library(lgr)
source('functions/functions.R')

#version 3.55.0 as of 1/31/24
alignments<-readRDS('ref/alignments.rda')
c2_antigen_ref<<-getC2RefAlleles()

con<-dbConn()
r<-160049
mg<-getMatchGrade(con, r)

antigen_ref<<-getAntigenTable(con)

#need to display donors for evaluation 
donors<-getDonors(mg)

r_mg<-mg %>%
  filter(donor_number == r)
recip_typing<-getTyping(r_mg, 'r')
recip_hla<-recip_typing[[1]]

#capture any recipient null alleles
recip_null_allele<-recip_typing[[2]]

#remove N from allele name if any null alleles
if(!is.null(recip_null_allele)){
  str_sub(recip_null_allele, -1)<-''
}

 #load NMDP file if there are any NMDP codes present in the typing data
nmdpCheck<-mg %>%
  select(c(a_1, a_2, b_1, b_2, c_1, c_2, dr_1, dr_2, drp_1, drp_2, dqb_1, dqb_2, dqa_1, dqa_2,
           dpa_1, dpa_2, dpb_1, dpb_2)) %>%
  #check if character directly after the colon is a letter
  mutate_all(~grepl(':[A-WYZ]', .))


if(any(nmdpCheck)){
  nmdp_file<-loadNMDP()
}

d<-160866
for(d in donors){
  print(d)
  
  d_mg<-filterDonor(mg, d)
  donor_hla<-getTyping(d_mg, 'd')[[1]]
  
  ##### MATCH EVALUATIONS #####

  ABC<-calcABCDRB('ABC', donor_hla, recip_hla)

  DRB1<-calcABCDRB('DRB1', donor_hla, recip_hla)

  DRB345<-calcABCDRB('DRB', donor_hla, recip_hla)

  DP<-calcDQDP('DP', donor_hla, recip_hla)

  DQ<-calcDQDP('DQ', donor_hla, recip_hla)

  #ABC-DR-DQ
  #match, total
  #ABC category should always have something

  ABCDRDQ<-ABC[[1]]+DRB1[[1]]+DQ[[1]]
  
  #print(ABCDRDQ[c(1,2)])
  
  #ABC-DRB1
  #match, total, mm gvh, mm hvg

  ABCDRB1<-ABC[[1]]+DRB1[[1]]
  
  #print(ABCDRB1)
  
  #DRB3/4/5-DQ-DP
  #match, total, mm gvh, mm hvg
  if(sum(DRB345[[1]])==0 & sum(DQ[[1]])==0 & sum(DP[[1]])==0){
    DRB345DQDP<-c(NA, NA, NA, NA)
  } else{
    DRB345DQDP<-DRB345[[1]]+DQ[[1]]+DP[[1]]
  }
  
  #print(DRB345DQDP)
  
  #7 LOCI TOTAL
  #match, total
  #return NA if any of the 7 loci have not been sequenced
  if(ABC[[1]][2]==0 | DRB1[[1]][2]==0 | DRB345[[1]][2]==0 | DP[[1]][2]==0 | DQ[[1]][2]==0){
    allLoci<-c(NA, NA)
  } else{
    allLoci<-c(ABCDRB1[1] + DRB345DQDP[1], ABCDRB1[2] + DRB345DQDP[2])
  }
  
 #print(allLoci)
  
  #disable this for real values while testing
  #updateMGtable(con, 'match', c(ABCDRDQ[c(1,2)], ABCDRB1, DRB345DQDP, allLoci, d))

  ##### TCE EVALUATIONS #####
  tce<-getTCE(donor_hla, recip_hla  )
  
  #print(tce)
  if(!is.null(tce)){
    #updateMGtable(con, 'tce', c(tce,d))
  }
  
  ##### DSA EVALUATIONS #####
  mm_alleles<-c(ABC[[2]], DRB1[[2]], DRB345[[2]], DQ[[2]], DP[[2]])
  
  #append non-null variant to mm_alleles if the recipient has a NULL allele
  if(!is.null(recip_null_allele)){
    mm_alleles<-c(mm_alleles, recip_null_allele)
  }
  
  DSA<-calcDSA(con, mm_alleles, r)
  print(DSA)
  
  #updateMGtable(con, 'dsa', c(DSA, d))
}





