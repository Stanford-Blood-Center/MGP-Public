#v 1.1.1

options(warn = 2) 

suppressPackageStartupMessages(library(lgr))

calcMatchGrade<-function(r_itl, d_itl, username){
  
  #set up log file
  log<-tempfile(tmpdir=tempdir(),pattern = paste("MG_", r_itl, '_', sep=""), fileext = ".log")
  lgr$add_appender(AppenderFile$new(log), name = "mg_log")
  
  lgr$info(paste('Executed by mTilda user: ', username, sep = ''))
  sunetID<-Sys.getenv('SHINYPROXY_USERNAME')
  lgr$info(paste('SUNetID: ', sunetID, sep = ''))
  
  tryCatch(
    {
      lgr$info('**********MATCH GRADE EVALUATION START**********')
      #version 3.55.0 as of 1/31/24
      alignments<<-readRDS('ref/alignments.rda')
      c2_antigen_ref<<-getC2RefAlleles()
      
      
      lgr$info(paste('Recipient ITL entered:', r_itl))
      
      con<-dbConn()
      lgr$info('Successfully connected to the database')
      
      #get Antigen table from DB
      antigen_ref<<-getAntigenTable(con)
      
      mg<-getMatchGrade(con, r_itl)
      
      if(nrow(mg)==0){
        mess<-paste('ITL', r_itl, 'was not found in Match Grade. Please make sure the ITL entered is correct.')
        lgr$info(mess)
        return()
      }
      
      lgr$info('Match Grade data extracted')
      
      r_mg<-mg %>%
        filter(donor_number == r_itl)
      recip_typing<-getTyping(r_mg, 'r')
      recip_hla<-recip_typing[[1]]
      
      #capture any recipient null alleles
      recip_null_allele<-recip_typing[[2]]
      
      lgr$info('Recipient typing extracted')
      
      #load NMDP file if there are any NMDP codes present in the typing data
      nmdpCheck<-mg %>%
        select(c(a_1, a_2, b_1, b_2, c_1, c_2, dr_1, dr_2, drp_1, drp_2, dqb_1, dqb_2, dqa_1, dqa_2,
                 dpa_1, dpa_2, dpb_1, dpb_2)) %>%
        #check if character directly after the colon is a letter
        mutate_all(~grepl(':[A-WYZ]', .))
      
      if(any(nmdpCheck)){
        nmdp_file<<-loadNMDP()
        lgr$info('NMDP typing detected. NMDP conversion file loaded.')
      }
      
      #get recipient's called antibodies
      sample_num<-getSampleNumber(con, r_itl)
      
      calculateDSA<-TRUE
      
      #if the sample number is NA, 'Recipient DSA Date' was not populated
      #field is not to be populated if only C1Q was ordered for a patient
      if(is.na(sample_num)){
        calculateDSA<-'Unknown'
      } else{
        #get IgG test numbers to get MFI values and AB screening results
        #for IgG specific tests
        testNums<-paste(getIgGTestNums(con, sample_num), collapse=',')
        
        ab_results<-getAbResults(con, testNums)
        ab_results$called_antibodies<-str_trim(ab_results$called_antibodies)
        
        #if called antibodies for both classes are 'Negative', DSA = N
        if(all(ab_results$called_antibodies == 'Negative')){
          calculateDSA<-FALSE
        } else{
          #recipient positive antigens 
          positive_antigens<-ab_results %>%
            filter(called_antibodies != 'Negative') %>%
            pull(called_antibodies) %>%
            strsplit(., ' ') %>%
            unlist()
          
          #split any haplotypes into separate alleles
          positive_antigens<-unlist(sapply(positive_antigens, function(x) if(grepl('/', x)) unlist(strsplit(x, '/')) else x), use.names = F)

          mfi_vals<-getMFIvals(con, r_itl, testNums)
          
          #if value is blank for average value, all beads in that antigen group have
          #reactivity 0; make blanks and NA 0 
          mfi_vals<-mfi_vals %>%
            mutate(average_value = case_when(average_value == '' ~ 0, 
                                             .default = as.numeric(average_value)))
          
          mfi_vals$average_value[which(is.na(mfi_vals$average_value))]<-0
        }
        
      }
      

      #will need to change this based on user input; will not be doing for llgr$oop
      #to go over all possibilities. 
      for(d in d_itl){
        
        lgr$info(paste('*****Calculating Match Grade for donor ITL ', d, '*****', sep=''))
        
        d_mg<-filterDonor(mg, d)
        donor_hla<-getTyping(d_mg, 'd')[[1]]
        
        ##### MATCH EVALUATIONS #####
        
        lgr$info('Donor typing extracted')
        
        lgr$info('Assessing A, B, and C...')
        ABC<-calcABCDRB('ABC', donor_hla, recip_hla)
        lgr$info('Finished assessing A, B, and C!')
        lgr$info('Assessing DRB1...')
        DRB1<-calcABCDRB('DRB1', donor_hla, recip_hla)
        lgr$info('Finished assessing DRB1!')
        lgr$info('Assessing DRB3/4/5...')
        DRB345<-calcABCDRB('DRB', donor_hla, recip_hla)
        lgr$info('Finished assessing DRB3/4/5!')
        lgr$info('Assessing DP...')
        DP<-calcDQDP('DP', donor_hla, recip_hla)
        lgr$info('Finished assessing DP!')
        lgr$info('Assessing DQ...')
        DQ<-calcDQDP('DQ', donor_hla, recip_hla)
        lgr$info('Finished assessing DQ!')
        
        #ABC-DR-DQ
        #match, total
        #ABC category should always have something
        lgr$info('Summing values for the ABC-DR-DQ category...')
        
        ABCDRDQ<-ABC[[1]]+DRB1[[1]]+DQ[[1]]
        
        print(ABCDRDQ[c(1,2)])
        
        #ABC-DRB1
        #match, total, mm gvh, mm hvg
        lgr$info('Summing values for the ABC-DRB1 category...')
        
        ABCDRB1<-ABC[[1]]+DRB1[[1]]
        
        print(ABCDRB1)
        
        #DRB3/4/5-DQ-DP
        #match, total, mm gvh, mm hvg
        if(sum(DRB345[[1]])==0 & sum(DQ[[1]])==0 & sum(DP[[1]])==0){
          DRB345DQDP<-c(NA, NA, NA, NA)
        } else{
          lgr$info('Summing values for the DRB3/4/5-DQ-DP category...')
          DRB345DQDP<-DRB345[[1]]+DQ[[1]]+DP[[1]]
        }
        
        print(DRB345DQDP)
        
        #7 LOCI TOTAL
        #match, total
        #return NA if any of the 7 loci have not been sequenced
        if(ABC[[1]][2]==0 | DRB1[[1]][2]==0 | DRB345[[1]][2]==0 | DP[[1]][2]==0 | DQ[[1]][2]==0){
          lgr$info('Skipping value summation for the 7 loci category, since not all loci have been sequenced')
          allLoci<-c(NA, NA)
        } else{
          lgr$info('Summing values for the 7 loci category...')
          allLoci<-c(ABCDRB1[1] + DRB345DQDP[1], ABCDRB1[2] + DRB345DQDP[2])
        }
        
        print(allLoci)
        
        #disable this for real values while testing
        lgr$info('Updating match values in Match Grade')
        updateMGtable(con, 'match', c(ABCDRDQ[c(1,2)], ABCDRB1, DRB345DQDP, allLoci, d))
        lgr$info('Match values successfully updated!')
        
        ##### TCE EVALUATIONS #####
        lgr$info('Assessing TCE permissibility...')
        tce<-getTCE(donor_hla, recip_hla)
        
        print(tce)
        if(!is.null(tce)){
          lgr$info('Finished assessing TCE permissibility!')
          updateMGtable(con, 'tce', c(tce,d))
          lgr$info('TCE value successfully updated!')
        }
        
        ##### DSA EVALUATIONS #####
        mm_alleles<-c(ABC[[2]], DRB1[[2]], DRB345[[2]], DQ[[2]], DP[[2]])
        
        #remove N suffix and combine with mm_alleles to evaluate DSA for non-null
        #variant
        if(!is.null(recip_null_allele)){
          str_sub(recip_null_allele, -1)<-''
          mm_alleles<-c(mm_alleles, recip_null_allele)
        }
        
        lgr$info('Evaluating DSA...')
        if(calculateDSA == TRUE){
          DSA<-calcDSA(con, mm_alleles, positive_antigens, mfi_vals)
        } else if(calculateDSA == FALSE){
          DSA<-'N'
        } else if (calculateDSA == 'Unknown'){
          DSA<-'U'
        }
        
        lgr$info(paste('DSA:', DSA))
        
        lgr$info('Updating DSA in Match Grade')
        updateMGtable(con, 'dsa', c(DSA, d))
        lgr$info('Finished updating DSA in Match Grade')
        
        lgr$info(paste('*****Finished calculating Match Grade for donor ITL ', d, '*****', sep=''))
      }
      return(c(TRUE, log))
      
    },
    error = function(e){
      lgr$fatal(e)
      return(c(FALSE, log))
    }
  )
}