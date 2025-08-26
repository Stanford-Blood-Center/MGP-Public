#Match Grade Populator © Stanford Blood Center, LLC.
#v 1.12.9

options(warn = 2) 

suppressPackageStartupMessages(library(lgr))

calcMatchGrade<-function(r_itl, d_itl, credentials, recip_hla, donor_hla, synnonsynList, d_mg, d_filter, r_filter){
  
  tryCatch(
    {
      lgr$info('**********MATCH GRADE EVALUATION START**********')

      lgr$info(paste('Recipient ITL entered:', r_itl))
      errorMessage<-NULL
      
      con<-dbConn()

      #get Antigen table from DB
      antigen_ref<<-getAntigenTable(con)
      
      #get recipient's called antibodies
      sample_num<-getSampleNumber(con, r_itl)
      
      lgr$info(sprintf('Sample Number: %s', sample_num))
      
      calculateDSA<-TRUE
      
      #if the sample number is NA, 'Recipient DSA Date' was not populated
      #field is not to be populated if only C1Q was ordered for a patient
      if(is.na(sample_num)){
        calculateDSA<-'Unknown'
      } else{
        #get IgG test numbers to get MFI values and AB screening results
        #for IgG specific tests
        testNums<-paste(getIgGTestNums(con, sample_num), collapse=',')
        
        lgr$info(sprintf('Test Numbers: %s', testNums))
        
        if(testNums == ""){
          errorMessage<-'The selected DSA date in the Match Grade software does not have any associated IgG tests. Please check the selected DSA date.'
          stop('Selected DSA date does not have associated IgG tests')
        }
        
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
          #positive_antigens<-unlist(sapply(positive_antigens, function(x) if(grepl('/', x)) unlist(strsplit(x, '/')) else x), use.names = F)
          
          mfi_vals<-getMFIvals(con, r_itl, testNums)

          #if value is blank for average value, all beads in that antigen group have
          #reactivity 0; make blanks and NA 0 
          mfi_vals<-mfi_vals %>%
            mutate(average_value = case_when(average_value == '' ~ 0, 
                                             .default = as.numeric(average_value)))
          
          mfi_vals$average_value[which(is.na(mfi_vals$average_value))]<-0
          
          #filter out Bw4, Bw6
          #if DP or DQ, only keep heterodimers 
          mfi_vals<-mfi_vals %>%
            filter(!antigen %in% c('Bw4', 'Bw6')) %>% 
            filter(!is.na(allele) & allele != '') %>%
            filter(!is.na(antigen)) %>%
            filter(case_when(grepl('DP|DQ', antigen)~grepl(',', probe_id),
                             TRUE ~ TRUE))
          
          mfi_vals$probe_id<-gsub(',', '/', mfi_vals$probe_id)
        }
      }

      lgr$info(paste('*****Calculating Match Grade for donor ITL ', d_itl, '*****', sep=''))
      
      ##### MATCH EVALUATIONS #####
      
      lgr$info('Donor typing extracted')
      lgr$info('Assessing A, B, and C...')
      ABC<-calcABCDRB('ABC', donor_hla, recip_hla, synnonsynList, d_filter, r_filter)
      lgr$info('Finished assessing A, B, and C!')
      lgr$info('Assessing DRB1...')
      DRB1<-calcABCDRB('DRB1', donor_hla, recip_hla, synnonsynList, d_filter, r_filter)
      lgr$info('Finished assessing DRB1!')
      lgr$info('Assessing DRB3/4/5...')
      DRB345<-calcABCDRB('DRB', donor_hla, recip_hla, synnonsynList, d_filter, r_filter)
      lgr$info('Finished assessing DRB3/4/5!')
      lgr$info('Assessing DP...')
      DP<-calcDQDP('DP', donor_hla, recip_hla, synnonsynList, d_filter, r_filter)
      lgr$info('Finished assessing DP!')
      lgr$info('Assessing DQ...')
      DQ<-calcDQDP('DQ', donor_hla, recip_hla, synnonsynList, d_filter, r_filter)
      lgr$info('Finished assessing DQ!')
      
      #ABC-DR-DQ
      #match, total
      #ABC category should always have something
      lgr$info('Summing values for the ABC-DR-DQ category...')
      
      ABCDRDQ<-ABC[[1]]+DRB1[[1]]+DQ[[1]]
      
      #print(ABCDRDQ[c(1,2)])
      
      #ABC-DRB1
      #match, total, mm gvh, mm hvg
      lgr$info('Summing values for the ABC-DRB1 category...')
      ABCDRB1<-ABC[[1]]+DRB1[[1]]
      
      #print(ABCDRB1)
      #DRB3/4/5-DQ-DP
      #match, total, mm gvh, mm hvg
      if(sum(DRB345[[1]])==0 & sum(DQ[[1]])==0 & sum(DP[[1]])==0){
        DRB345DQDP<-c(NA, NA, NA, NA)
      } else{
        lgr$info('Summing values for the DRB3/4/5-DQ-DP category...')
        DRB345DQDP<-DRB345[[1]]+DQ[[1]]+DP[[1]]
      }
      
      #print(DRB345DQDP)
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
      
      #print(allLoci)
      
      #disable this for real values while testing
      if(credentials == '30'){
        lgr$info('Updating match values in Match Grade')
        updateMGtable(con, 'match', c(ABCDRDQ[c(1,2)], ABCDRB1, DRB345DQDP, allLoci, d_itl))
        lgr$info('Match values successfully updated!')
      } else if(credentials %in% c('50', '60')){
        d_mg$ABCDRDQ_match<-ABCDRDQ[1]
        d_mg$ABCDRDQ_alleles<-ABCDRDQ[2]
        
        d_mg$ABCDRB1_match<-ABCDRB1[1]
        d_mg$ABCDRB1_alleles<-ABCDRB1[2]
        d_mg$ABCDRB1_mm_GVH<-ABCDRB1[3]
        d_mg$ABCDRB1_mm_HVG<-ABCDRB1[4]
        
        d_mg$DRB345DQDP_match<-DRB345DQDP[1]
        d_mg$DRB345DQDP_alleles<-DRB345DQDP[2]
        d_mg$DRB345DQDP_mm_GVH<-DRB345DQDP[3]
        d_mg$DRB345DQDP_mm_HVG<-DRB345DQDP[4]
        
        d_mg$seven_loci_match<-allLoci[1]
        d_mg$Seven_loci_alleles<-allLoci[2]
      }
      
      ##### TCE EVALUATIONS #####
      lgr$info('Assessing TCE permissibility...')
      tce<-getTCE(donor_hla, recip_hla)

      if(!is.null(tce)){
        if(credentials == '30'){
          lgr$info('Finished assessing TCE permissibility!')
          updateMGtable(con, 'tce', c(tce,d_itl))
          lgr$info('TCE value successfully updated!')
        } else if(credentials %in% c('50', '60')){
          d_mg$DRB345DQDP_mm_TCE<-tce
        }
      }
      
      gvh_mm_alleles<-c(ABC[[3]], DRB1[[3]], DRB345[[3]], DQ[[3]], DP[[3]])
      
      ##### DSA EVALUATIONS #####
      hvg_mm_alleles<-c(ABC[[2]], DRB1[[2]], DRB345[[2]], DQ[[2]], DP[[2]])
      hvg_mm_alleles_eval<-c(unlist(hvg_mm_alleles, use.names = FALSE))
      
      lgr$info('Evaluating DSA...')
      if(calculateDSA == TRUE){
        DSAresults<-calcDSA(con, hvg_mm_alleles_eval, positive_antigens, mfi_vals, unlist(donor_hla, use.names = F), unlist(recip_hla, use.names = F))
        DSA<-DSAresults[[1]]
        DSAmessage<-DSAresults[[2]]
      } else if(calculateDSA == FALSE){
        DSA<-'N'
        DSAmessage<-data.frame()
      } else if (calculateDSA == 'Unknown'){
        DSA<-'U'
        DSAmessage<-data.frame()
      }
      
      lgr$info(paste('DSA:', DSA))
      
      if(credentials == '30'){
        lgr$info('Updating DSA in Match Grade')
        updateMGtable(con, 'dsa', c(DSA, d_itl))
        lgr$info('Finished updating DSA in Match Grade')
      } else if(credentials %in% c('50', '60')){
        d_mg$DSA<-DSA
      }
      
      missing_alleles<-paste(c(ABC[[4]], DRB1[[4]], DRB345[[4]], DQ[[4]], DP[[4]]), collapse = ', ')
      final_missing_message<-NULL
      
      if(missing_alleles != ""){
        final_missing_message<-paste('The following alleles do not have a complete sequence in the IMGT reference alignment for the ARD: <b>', missing_alleles, '</b>. Mismatches may be incorrect. Please check counts.', sep = '')
      }
      A_mm<-B_mm<-C_mm<-DRB1_mm<-DRB345_mm<-DQA1_mm<-DQB1_mm<-DPA1_mm<-DPB1_mm<-NULL
      
      if('A' %in% names(ABC[[3]])){
        A_mm<-mismatchFormatter('A', ABC[[3]][[1]], ABC[[2]][[1]])
      }
      if('B' %in% names(ABC[[3]])){
        B_mm<-mismatchFormatter('B', ABC[[3]][[2]], ABC[[2]][[2]])
      }
      if('C' %in% names(ABC[[3]])){
        C_mm<-mismatchFormatter('C', ABC[[3]][[3]], ABC[[2]][[3]])
      }
      
      if('DRB1' %in% names(DRB1[[3]])){
        DRB1_mm<-mismatchFormatter('DRB1', DRB1[[3]][[1]], DRB1[[2]][[1]])
      }

      if('DRB' %in% names(DRB345[[3]])){
        DRB345_mm<-mismatchFormatter('DRB3/4/5', DRB345[[3]][[1]], DRB345[[2]][[1]])
      }

      if('DQA1' %in% names(DQ[[3]])){
        DQA1_mm<-mismatchFormatter('DQA1', DQ[[3]][[1]], DQ[[2]][[1]])
      }
      if('DQB1' %in% names(DQ[[3]])){
        DQB1_mm<-mismatchFormatter('DQB1', DQ[[3]][[2]], DQ[[2]][[2]])
        
      }
      if('DPA1' %in% names(DP[[3]])){
        DPA1_mm<-mismatchFormatter('DPA1', DP[[3]][[1]], DP[[2]][[1]])
      }
      if('DPB1' %in% names(DP[[3]])){
        DPB1_mm<-mismatchFormatter('DPB1', DP[[3]][[2]], DP[[2]][[2]])
      }

      lgr$info(paste('*****Finished calculating Match Grade for donor ITL ', d_itl, '*****', sep=''))
      
      return(list('TRUE', d_mg[c(36:49)], final_missing_message, errorMessage, A_mm, B_mm, C_mm, DRB1_mm, DRB345_mm, DQA1_mm, DQB1_mm, DPA1_mm, DPB1_mm, DSAmessage))
      
    },
    error = function(e){
      lgr$fatal(e)
      return(list('FALSE', NULL, NULL, errorMessage))
    }
  )
}
