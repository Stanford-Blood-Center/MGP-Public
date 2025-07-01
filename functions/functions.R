#external functions
#Match Grade Populator © Stanford Blood Center, LLC.
#v 1.12.6

suppressPackageStartupMessages(library(odbc))
suppressPackageStartupMessages(library(tidyverse))
library(stringr) 
library(httr)
suppressPackageStartupMessages(library(rvest))

dbConn <- function(){
  
  con <- dbConnect(odbc(),
                   Driver = Sys.getenv('DRIVER'),
                   Server = Sys.getenv('SERVER'),
                   Database = Sys.getenv('DB'),
                   UID = Sys.getenv('DB_USERNAME'),
                   PWD = Sys.getenv('DB_PW'))
  
  return(con)
}

#extract Match Grade data for recipient ITL entered
getMatchGrade<-function(con, r_itl){
  
  mg_res<-dbGetQuery(con, paste("SELECT * 
                                    FROM dbo.Match_grades 
                                    WHERE recipient_number = ", r_itl, sep = ""))
  
  return(mg_res)
}

#extract available donors for evaluation
getDonors<-function(mg_df){
  donors <- mg_df %>%
    filter(recipient_number != donor_number) %>%
    select(donor_number) %>%
    pull()
  
  return(donors)
}

#get HLA typing
#transform loci columns into rows, populate homozygous alleles
getTyping<-function(con, mgDF, type){
  
  alignments<<-readRDS('ref/alignments.rda')

  hla_cols <- sort(do.call(paste, c(expand.grid(c('a', 'b', 'c', 'dr', 'drp', 'dqa', 'dqb', 'dpa', 'dpb'), c(1,2)), sep = '_')))
  
  if(type == 'r'){
    prefix <- 'Recipient'
  } else{
    prefix <- 'Donor'
  }
  
  t_df <- mgDF %>%
    select(all_of(hla_cols)) %>%
    mutate(id = row_number()) %>%
    pivot_longer(-id) %>%
    mutate(name = toupper(sub("^(.*)_.*", "\\1", name))) %>%
    pivot_wider(names_from = 'name',
                values_from = 'value',
                values_fn = list) %>%
    unnest(cols = colnames(.)) %>%
    select(-id) %>%
    as.data.frame(stringsAsFactors = F)
  
  #paste standard locus names to allele name
  colnames(t_df)[names(t_df) == 'DPA'] <- 'DPA1'
  colnames(t_df)[names(t_df) == 'DPB'] <- 'DPB1'
  colnames(t_df)[names(t_df) == 'DQA'] <- 'DQA1'
  colnames(t_df)[names(t_df) == 'DQB'] <- 'DQB1'
  colnames(t_df)[names(t_df) == 'DRP'] <- 'DRB'
  colnames(t_df)[names(t_df) == 'DR'] <- 'DRB1'
  
  
  #find homozygous alleles (X or NA)
  ds_allele<-t_df %>%
    select(where(~any(. == "X" | is.na(.)))) %>%
    names() 
  
  #populate X or NA w extracted homozygous allele
  for(i in ds_allele){
    if(t_df[,i][1]=='X'|is.na(t_df[,i][1])){
      t_df[,i][1]<-t_df[,i][2]
    } else{
      t_df[,i][2]<-t_df[,i][1]
    }
  }
  
  #get any alleles with fully sequenced, null alleles
  null_alleles<-t_df %>%
    select(where(~any(str_sub(., -1) == 'N' & !is.na(.) & str_count(., '[[:alpha:]]') == 1))) %>%
    names()
  
  filter<-c()
  
  for(a in null_alleles){
    getNullAllele<-t_df %>% 
      select(all_of(a)) %>%
      pull()
    
    nullVars<-unique(keep(getNullAllele, ~ str_sub(.x, -1) == "N" & str_count(.x, "[[:alpha:]]") == 1))
    
    if(a == 'DRB'){
      delim = ''
    } else{
      delim = '*'
    }
    nullAppend<-paste(a, delim, nullVars, sep = '')
    lgr$info(paste(prefix, '- non-expression variant', nullAppend, 'detected', sep = ' '))
    filter<-append(filter, nullAppend)
  }
  
  r_null_allele<-c()
  
  if(type == 'r'){
    prefix <- 'Recipient'
    
    if(length(null_alleles)!=0){
      for(i in 1:length(null_alleles)){
        if(null_alleles[[i]] == 'DRB'){
          delim = ''
        } else{
          delim = '*'
        }
        r_null_allele<-append(r_null_allele, paste(null_alleles[[i]], t_df[,null_alleles[[i]]][str_sub(t_df[,null_alleles[[i]]], -1) == 'N'], sep = delim))
      }
    }
  } else{
    prefix <- 'Donor'
  }
  
  #if a locus has heterozygous alleles and does not have NMDP alleles, 
  #check if they have the same P group
  #skip NMDP, NAs, and novel alleles (denoted by @ suffix)
  residue<-NULL
  novelAllelesAll<-c()
  
  for(j in colnames(t_df)){
    if(any(grepl('@', t_df[,j]))){
      if(j == 'DRB'){
        delimit <- ''
      } else{
        delimit <- '*'
      }
      
      novelAlleles<-paste('HLA-', j, delimit, t_df[,j][grepl('@', t_df[,j])], sep = '')
      lgr$info(sprintf('%s: Novel allele %s detected', prefix, novelAlleles))
      #capture all novel alleles 
      novelAllelesAll<-append(novelAllelesAll, paste(prefix, novelAlleles, sep=': '))
    }
    
    if(!any(grepl(':[A-WYZ]', t_df[,j])) & !any(grepl('@', t_df[,j])) & !all(is.na(t_df[,j]))){
      
      if(t_df[1,j] != t_df[2,j]){
        
        delim<-'*'
        
        align_locus<-j
        
        residue<-90
        
        if(j %in% c('A', 'B', 'C')){
          residue<-182
        } 
        
        if(j == 'DRB'){
          drb_loci<-gsub('^(.*?)\\*.*$', '\\1', t_df[,j])
          #move onto next column if DRB loci are not the same; no need to compare
          if(drb_loci[1] != drb_loci[2]){
            next
          }
          delim<-''
          align_locus<-paste('DRB', unique(drb_loci), sep ='')
        } 
        
        s_pg_alleles<-alignments[[align_locus]] %>%
          filter(trimmed_allele %in% c(paste(j, t_df[,j], sep=delim))) %>%
          distinct(trimmed_allele, .keep_all = TRUE) %>%
          select(as.character(seq(1, residue))) %>%
          unite(., 'all', sep='') %>%
          select('all') %>%
          summarise(equal = all[1] == all[2]) %>%
          pull()
        
        #if two alleles in the same locus have the same P group, it is effectively
        #homozygous. append it to the filter variable to remove one during mismatch
        #evaluation
        #change align_locus back to DRB if locus is DRB
        if(j == 'DRB'){
          align_locus<-'DRB'
          delim<-''
        } else{
          delim<-'*'
        }
        
        if(s_pg_alleles){
          lgr$info(paste(prefix, ' - ', align_locus, delim, t_df[2,j], ' and ', align_locus, delim, t_df[1,j], ' have the same P-group', sep = ''))
          #t_df[2,j]<-paste(t_df[1,j], t_df[2,j], sep = ',')
          filter<-append(filter, paste(align_locus, delim, t_df[2,j], sep = ''))
        }
      }
    }
  }
  
  t_df<-t_df %>%
    rowwise() %>%
    mutate(across(everything(), ~ifelse(!is.na(.), ifelse(cur_column()=="DRB", paste(cur_column(), ., sep = ""), paste(cur_column(), ., sep = "*")), .)))
  
  return(list(t_df, r_null_allele, novelAllelesAll, filter))
}


#convert NMDP code to subtype 
convertNMDP<-function(n_allele){
  
  f1<-gsub('.*?(.*?)\\*', "\\1", sapply(strsplit(n_allele, ':'), '[[', 1))
  f2<-sapply(strsplit(n_allele, ':'), '[[', 2) 
  
  donor_st<-sapply(f1, function(x) NULL)
  
  subtypes<-GET(
    url=sprintf('https://hml.nmdp.org/mac/api/codes/%s', f2), 
    add_headers(Accept = 'text/plain')
  )
  
  if(subtypes$status_code != 200){
    stop(sprintf('Status code %s returned from the MAC API for %s. Terminating session...', subtypes$status_code, n_allele))
  }
  
  #get subtypes
  subtype<-strsplit(strsplit(trimws(gsub('\t|\\*', ' ', content(subtypes, encoding = 'UTF-8'))), ' ')[[1]][[2]], '/')[[1]]
  
  #if subtype does not already have first two fields, paste it
  if(any(!grepl(':', subtype))){
    donor_st<-paste(names(donor_st), subtype, sep=':')
  } else{
    donor_st<-subtype
  }
  
  return(unlist(donor_st, use.names=F))
  
}

#update dbo.Match_grades table with calculated values
updateMGtable<-function(connection, type, vals){
  
  #updating match columns
  if(type=='match'){
    
    res<-dbSendStatement(connection, 'UPDATE dbo.Match_grades
                        SET dbo.Match_grades.ABCDRDQ_match=?, 
                            dbo.Match_grades.ABCDRDQ_alleles=?,
                            dbo.Match_grades.ABCDRB1_match=?,
                            dbo.Match_grades.ABCDRB1_alleles=?,
                            dbo.Match_grades.ABCDRB1_mm_GVH=?,
                            dbo.Match_grades.ABCDRB1_mm_HVG=?,
                            dbo.Match_grades.DRB345DQDP_match=?,
                            dbo.Match_grades.DRB345DQDP_alleles=?,
                            dbo.Match_grades.DRB345DQDP_mm_GVH=?,
                            dbo.Match_grades.DRB345DQDP_mm_HVG=?,
                            dbo.Match_grades.Seven_loci_match=?,
                            dbo.Match_grades.Seven_loci_alleles=?
                        WHERE dbo.Match_grades.donor_number=?')
    
  }
  
  
  #updating tce column
  if(type=='tce'){
    
    res<-dbSendStatement(connection, 'UPDATE dbo.Match_grades
                        SET dbo.Match_grades.DRB345DQDP_mm_TCE=?
                        WHERE dbo.Match_grades.donor_number=?')
    
  }
  
  if(type == 'dsa'){
    
    res<-dbSendStatement(connection, 'UPDATE dbo.Match_grades
                        SET dbo.Match_grades.DSA=?
                        WHERE dbo.Match_grades.donor_number=?')
    
  }
  
  dbBind(res, c(vals))
  dbClearResult(res)
}

#get TCE permissibility for TCE
getTCE<-function(d_hla, r_hla){
  
  d_dpb<-d_hla %>%
    select(DPB1) %>%
    pull() %>%
    strsplit('\\*') %>%
    map(2) %>%
    unlist()
  
  
  if(all(is.na(d_dpb))){
    lgr$info('Donor has no DPB1 alleles to assess')
    return()
  }
  
  r_dpb<-r_hla %>%
    select(DPB1) %>%
    pull() %>%
    strsplit('\\*') %>%
    map(2) %>%
    unlist()
  
  #replace any null variants with the non-null DPB variant
  if(any(grepl('N', r_dpb))){
    r_dpb[which(grepl('N', r_dpb))]<-r_dpb[!grepl('N', r_dpb)]
  } else if(any(grepl('N', d_dpb))){
    d_dpb[which(grepl('N', d_dpb))]<-d_dpb[!grepl('N', d_dpb)]
  }
    
  legend<-list('Non-Permissive'='N', 'Permissive'='P','Unknown'='U')
  
  #get differences between recipient DP alleles and donor DP alleles
  rDPdiff<-setdiff(r_dpb, d_dpb)
  dDPdiff<-setdiff(d_dpb, r_dpb)
  
  dpDiff<-union(rDPdiff, dDPdiff)
  
  if(length(dpDiff)!=0){
    
    DPalignment<-alignments[['DPB1']]
    removeDiff <- c()
    
    #evaluate if mismatched donor alleles have the same p group as the recip alleles
    recipAlignments<-DPalignment %>%
      filter(trimmed_allele %in% c(paste('DPB1*', r_dpb, sep = ''))) %>%
      distinct(trimmed_allele, .keep_all = TRUE) %>%
      select(c(trimmed_allele,as.character(seq(1, 90)))) %>%
      unite(.,'all', -trimmed_allele, sep='')
    
    for(d in dDPdiff){
      donorEvalAlignments<-DPalignment %>%
        filter(trimmed_allele == paste('DPB1*', d, sep = '')) %>%
        slice(1) %>%
        select(as.character(seq(1, 90))) %>%
        paste(collapse = "")
      
      if(any(donorEvalAlignments == recipAlignments$all)){
        d_samePGroup<-paste(recipAlignments$trimmed_allele[donorEvalAlignments == recipAlignments$all], collapse = ', ')
        lgr$info(sprintf('Donor - %s removed from mismatched DPs as it has the same P-group as the following recipient alleles: %s', paste(paste('DPB1*', d, sep =''), collapse = ','), d_samePGroup))
        removeDiff<-append(removeDiff, d)
      }
    }
    
    #evaluate if mismatched recip alleles have the same p group as the donor alleles
    donorAlignments<-DPalignment %>%
      filter(trimmed_allele %in% c(paste('DPB1*', d_dpb, sep = ''))) %>%
      distinct(trimmed_allele, .keep_all = TRUE) %>%
      select(c(trimmed_allele,as.character(seq(1, 90)))) %>%
      unite(.,'all', -trimmed_allele, sep='')
    
    for(r in rDPdiff){
      recipEvalAlignments<-DPalignment %>%
        filter(trimmed_allele == paste('DPB1*', r, sep = '')) %>%
        slice(1) %>%
        select(as.character(seq(1, 90))) %>%
        paste(collapse = "")
      
      if(any(recipEvalAlignments == donorAlignments$all)){
        r_samePGroup<-paste(donorAlignments$trimmed_allele[recipEvalAlignments == donorAlignments$all], collapse = ', ')
        lgr$info(sprintf('Recipient - %s removed from mismatched DPs as it has the same P-group as the following donor alleles: %s', paste(paste('DPB1*', r, sep =''), collapse = ','), r_samePGroup))
        removeDiff<-append(removeDiff, r)
      }
    }
    dpDiff<-dpDiff[!dpDiff %in% removeDiff]
  }
  
  if(length(dpDiff)==0){
    tce<-'M'
  } else{
    
    api_body<-sprintf('pid=1&patdpb1=%s&patdpb2=%s&did=2&dondpb1=%s&dondpb2=%s',
                      r_dpb[1], r_dpb[2], d_dpb[1], d_dpb[2])
    
    res_text<-content(POST(
      url='https://www.ebi.ac.uk/cgi-bin/ipd/pl/hla/dpb_v2.cgi?',
      body=api_body
    )) %>%   
      html_nodes("td") %>%
      html_text() %>%
      .[17]
    
    tce<-strsplit(strsplit(res_text, ': ')[[1]][[2]], " ")[[1]][[1]]
    tce<-legend[[tce]]
  }
  
  return(tce)
}

#calculate matches for ABC, DRB3/4/5, and DRB1
#categories include the following:
#'ABC' - HLA-A, B, C
#'DRP' - HLA-DRB3/4/5
#'DR' - HLA-DRB1
calcABCDRB<-function(cat, d_hla, r_hla, synqList, filter_d, filter_r){
  
  residue<-90
  
  if(cat == 'ABC'){
    group<-strsplit(cat, '')[[1]]
    residue<-182
  } else{
    group<-cat
  }
  
  d_cat<-d_hla %>%
    select(all_of(group))
  
  #if category is DRP and all entries are NA, check if DRB1 has been sequenced
  #DRB1*01, DRB1*08, and DRB1*10 do not have associated DRB3/4/5 alleles
  ###DONOR###
  if(cat=='DRB'){
    if(all(is.na(d_cat[,cat]))){
      DRB1<-d_hla %>% 
        select('DRB1') %>%
        pull()
      
      #if DRB1 has been sequenced is not NMDP, then DRB3/4/5 is absent due to
      #not being expressed, and not due to lack of typing
      #change NA to ''
      if(all(!grepl(':[A-WYZ]', DRB1))){
        lgr$info('Donor: DRB3/4/5 empty and DRB1 is fully sequenced. HvG will be 0')
        d_cat[,cat]<-c('', '')
      } 
    }
  }
  
  #return 0 if selected category has not been sequenced in donor genotype
  #this could happen for DRB3/4/5 or HLA-C
  if(all(is.na(d_cat))){
    return(list(c(0,0,0,0), NULL, NULL, NULL))
  }
  
  #get loci that have not been sequenced
  na_cols<-names(which(colSums(is.na(d_cat)) > 0))
  
  if(length(na_cols)!=0){
    group<-setdiff(group, na_cols)
  }
  
  #only filter loci in recipient typing that have been sequenced in potential donor
  r_cat<-r_hla %>%
    select(all_of(group))
  
  #DRB3/4/5 check
  ###RECIPIENT###
  if(cat=='DRB'){
    if(all(is.na(r_cat[,cat]))){
      DRB1<-r_hla %>% 
        select('DRB1') %>%
        pull()
      
      #if DRB1 has been sequenced is not NMDP, then DRB3/4/5 is absent due to
      #not being expressed, and not due to lack of typing
      #change NA to ''
      if(all(!grepl(':[A-WYZ]', DRB1))){
        lgr$info('Recipient: DRB3/4/5 empty and DRB1 is fully sequenced. GvH will be 0')
        r_cat[,cat]<-c('', '')
      } 
    }
  }
  
  #return 0 if selected category for recipient only has NAs (i.e two non-expressed
  #alleles)
  if(all(is.na(r_cat))){
    return(list(c(0,0,0,0), NULL, NULL, NULL))
  }
  
  d_alleles<-na.omit(unlist(d_cat, use.names=F))
  r_alleles<-na.omit(unlist(r_cat, use.names=F))
  
  #novel allele handling - recipient
  if(any(grepl('@', r_alleles))){
    r_alleles<-processNovelAlleles(r_alleles, synqList)
  }
  
  #novel allele handling - donor
  if(any(grepl('@', d_alleles))){
    d_alleles<-processNovelAlleles(d_alleles, synqList)
  }
  
  #filter out one allele if 2 alleles have the same p group and filter out
  #null expression variant alleles for recipient and donor
  if(!is.null(filter_d)){
    d_filtered_alleles<-d_alleles[!d_alleles %in% filter_d]
  } else{
    d_filtered_alleles<-d_alleles
  }
  
  if(!is.null(filter_r)){
    r_filtered_alleles<-r_alleles[!r_alleles %in% filter_r]
  } else{
    r_filtered_alleles<-r_alleles
  }
  
  total<-matches<-length(d_alleles)
  
  gvh<-hvg<-sapply(group, function(x) 0)
  
  if(all(d_alleles %in% r_alleles) & all(r_alleles %in% d_alleles)){
    gvh<-hvg<-sapply(group, function(x) NULL)
    return(list(c(total, matches, 0, 0), hvg, gvh, NULL))
    
  } else{
    
    nmdp_flag<-reg_flag<-FALSE
    
    d_mm_alleles<-unique(d_filtered_alleles[which(!d_filtered_alleles %in% r_filtered_alleles)])
    r_mm_alleles<-unique(r_filtered_alleles[which(!r_filtered_alleles %in% d_filtered_alleles)])
    
    all_nmdp_alleles<-d_mm_alleles[isNMDP(d_mm_alleles)]
    nmdp_translated<-sapply(all_nmdp_alleles, function(x) NULL)
    
    if(length(nmdp_translated)!=0){
      lgr$info('NMDP typing detected. Converting NMDP alleles...')
      for(n in 1:length(nmdp_translated)){
        nmdpRes<-convertNMDP(names(nmdp_translated)[[n]])
        nmdp_translated[[n]]<-paste(gsub('^(.*?)\\*.*$', '\\1', names(nmdp_translated)[[n]]), nmdpRes, sep="*")
        
      }
    }
    
    missing_message<-c()
    gvh_alleles<-sapply(group, function(x) c())
    
    ##GvH calculation
    for(i in r_mm_alleles){
      
      nmdp_flag<-reg_flag<-FALSE
      
      if(i==''){
        gvh[['DRB']]<-0
        next
      } 
      
      align_locus<-mm_locus<-gsub('^(.*?)\\*.*$', '\\1', i)
      
      if(align_locus %in% c('DRB3', 'DRB4', 'DRB5')){
        mm_locus<-'DRB'
      }
      
      #filter to locus specific alleles for mismatched allele 
      #filtered alleles
      d_locus_alleles<-d_filtered_alleles[grepl(mm_locus,gsub('^(.*?)\\*.*$', '\\1', d_filtered_alleles))]
      #all locus alleles 
      all_d_locus_alleles<-d_alleles[grepl(mm_locus,gsub('^(.*?)\\*.*$', '\\1', d_alleles))]
      
      if(align_locus %in% c('DRB3', 'DRB4', 'DRB5')){
        #if DRB locus for recipient mismatch allele is not present in donor locus
        #alleles, add 1 to GvH and move onto next iteration
        if(all(align_locus != gsub('^(.*?)\\*.*$', '\\1', d_locus_alleles)) & !all(d_locus_alleles == '')){
          gvh_alleles[[mm_locus]]<-append(gvh_alleles[[mm_locus]], i)
          gvh[[mm_locus]]<-gvh[[mm_locus]]+1
          lgr$info(paste('GvH MM-', i, sep = ''))
          next
        }
        if(all(d_locus_alleles == '')){
          gvh_alleles[[mm_locus]]<-append(gvh_alleles[[mm_locus]], i)
          gvh[[mm_locus]]<-length(r_mm_alleles)
          lgr$info(paste('GvH MM-', i, sep = ''))
          next
        }
      } 
      
      #split up into nmdp and regular allele assessment, since some donors
      #can have nmdp and regular sequenced alleles
      #get locus specific translated NMDP alleles
      nmdp_alleles<-names(nmdp_translated)[grepl(mm_locus,gsub('^(.*?)\\*.*$', '\\1', names(nmdp_translated)))]
      reg_alleles<-setdiff(d_locus_alleles, nmdp_alleles)
      
      if(length(nmdp_alleles)!=0){
        
        nmdp_flag<-!i%in% unlist(nmdp_translated[c(nmdp_alleles)])
        
        #if mismatched recipient allele is in translated NMDP alleles, remove mismatched NMDP allele
        #from donor mismatched alleles
        if(!nmdp_flag){
          d_mm_alleles<-d_mm_alleles[!d_mm_alleles %in% names(nmdp_translated)[map(nmdp_translated, ~i%in% .)==TRUE]]
        }
        
        #assess for DRB1*04:07/DRB1*04:92 exception for NMDP codes
        #if category is DRB1, check if mm allele
        #is DRB1*04:07 or DRB1*04:92. if yes, check if the other allele
        #is present in the list of translated NMDP codes
        #if yes, change nmdp_flag = FALSE and remove the NMDP allele
        #from the donor mismatched alleles w that code
        if(cat=='DRB1'){
          if(i== 'DRB1*04:07'){
            if('DRB1*04:92' %in% unlist(nmdp_translated)){
              nmdp_flag<-FALSE
              d_mm_alleles<-d_mm_alleles[d_mm_alleles!=names(nmdp_translated)[map(nmdp_translated, ~'DRB1*04:92' %in% .)==TRUE]]
            }
            if(i== 'DRB1*04:92'){
              if('DRB1*04:07' %in% unlist(nmdp_translated)){
                nmdp_flag<-FALSE
                d_mm_alleles<-d_mm_alleles[d_mm_alleles!=names(nmdp_translated)[map(nmdp_translated, ~'DRB1*04:07' %in% .)==TRUE]]
              }
            }
          }
        }
      }
      
      if(length(reg_alleles)!=0){
        
        align<-alignments[[align_locus]]
        
        r_gvh_prot<-align %>%
          filter(trimmed_allele == i) %>%
          slice(1) %>%
          select(as.character(seq(1, residue))) %>%
          paste(collapse = "")
        
        if(grepl('*', r_gvh_prot, fixed = TRUE)){
          missing_message<-append(paste('RECIPIENT-', i, sep = ''), missing_message)
        }
        
        #exception for DRB1*04:92 and DRB1*04:07
        #check in all donor loci alleles, just in case the donor has 
        #DRB1 alleles with the same p group and one gets filtered out
        if(cat == 'DRB1'){
          if(i== 'DRB1*04:92'){
            if('DRB1*04:07' %in% all_d_locus_alleles){
              next 
            }
          } else if(i== 'DRB1*04:07'){
            if('DRB1*04:92' %in% all_d_locus_alleles){
              next
            }
          }
        }
        
        d_gvh_prot<-align %>%
          filter(trimmed_allele %in% c(d_locus_alleles)) %>%
          distinct(trimmed_allele, .keep_all = TRUE) %>%
          select(c(trimmed_allele, as.character(seq(1, residue)))) %>%
          unite(.,'all', -trimmed_allele, sep='')
        
        reg_flag<-!any(r_gvh_prot == d_gvh_prot$all)
        
        #if the recipient mature protein sequence is equivalent to one of the 
        #sequences for the mismatched donor allele, remove it from the list of
        #mismatched donor alleles to evaluate
        if(!reg_flag){
          filter<-d_gvh_prot$trimmed_allele[r_gvh_prot == d_gvh_prot$all]
          lgr$info(paste(i, ' and ', filter, ' have the same P-group', sep = ''))
          d_mm_alleles<-d_mm_alleles[!d_mm_alleles %in% filter]
        }
      }
      
      if(reg_flag | nmdp_flag){
        #if donor alleles consist of 1 nmdp and 1 regular allele, assess
        #if the recipient allele has FALSE for reg_flag or nmdp_flag
        #move onto next iteration without adding to the counter bc there is a match
        if(length(reg_alleles)==1 & length(nmdp_alleles)==1){
          if(!reg_flag | !nmdp_flag){
            next 
          }
        }
        gvh_alleles[[mm_locus]]<-append(gvh_alleles[[mm_locus]], i)
        gvh[[mm_locus]]<-gvh[[mm_locus]]+1
        lgr$info(paste('GvH MM-', i, sep = ''))
      }
    }
    
    hvg_alleles<-sapply(group, function(x) c())
    
    #HvG calculations
    for(j in d_mm_alleles){
      
      #add 0 to HvG for DRP if empty
      if(j==''){
        hvg[['DRB']]<-0
        next
      }
      
      align_locus<-mm_locus<-gsub('^(.*?)\\*.*$', '\\1', j)
      
      if(align_locus %in% c('DRB3', 'DRB4', 'DRB5')){
        mm_locus<-'DRB'
      }
      
      r_locus_alleles<-r_filtered_alleles[grepl(mm_locus,gsub('^(.*?)\\*.*$', '\\1', r_filtered_alleles))]
      all_r_locus_alleles<-r_alleles[grepl(mm_locus,gsub('^(.*?)\\*.*$', '\\1', r_alleles))]
      
      if(align_locus %in% c('DRB3', 'DRB4', 'DRB5')){
        if(all(r_locus_alleles == '')){
          hvg_alleles[[mm_locus]]<-append(hvg_alleles[[mm_locus]], j)          
          hvg[[mm_locus]]<-length(d_mm_alleles)
          lgr$info(paste('HvG MM-', j, sep = ''))
          next
        }
        if(all(align_locus != gsub('^(.*?)\\*.*$', '\\1', r_locus_alleles)) & !all(r_locus_alleles == '')){
          hvg_alleles[[mm_locus]]<-append(hvg_alleles[[mm_locus]], j)
          hvg[[mm_locus]]<-hvg[[mm_locus]]+1
          lgr$info(paste('HvG MM-', j, sep = ''))
          next
        }
      }
     
      #use all recipient alleles for evaluating NMDP codes
      if(isNMDP(j)){
        if(!any(nmdp_translated[[j]] %in% all_r_locus_alleles)){
          hvg_alleles[[mm_locus]]<-append(hvg_alleles[[mm_locus]], j)
          hvg[[mm_locus]]<-hvg[[mm_locus]]+1
          lgr$info(paste('HvG MM-', j, sep = ''))
        }
      } else{
        
        align<-alignments[[align_locus]]
        
        d_hvg_prot<-align %>%
          filter(trimmed_allele == j) %>%
          slice(1) %>%
          select(as.character(seq(1, residue))) %>%
          paste(collapse = "")
        
        if(grepl('*', d_hvg_prot, fixed = TRUE)){
          missing_message<-append(paste('DONOR-', j, sep = ''), missing_message)
        }
        
        if(cat == 'DRB1'){
          if(j == 'DRB1*04:92'){
            if('DRB1*04:07' %in% all_r_locus_alleles){
              next 
            }
          } else if(j == 'DRB1*04:07'){
            if('DRB1*04:92' %in% all_r_locus_alleles){
              next
            }
          }
        }
        
        r_hvg_prot<-align %>%
          filter(trimmed_allele %in% c(r_locus_alleles)) %>%
          distinct(trimmed_allele, .keep_all = TRUE) %>%
          select(as.character(seq(1, residue))) %>%
          unite(., 'all', sep='') %>%
          select('all')
        
        if(!any(d_hvg_prot == r_hvg_prot)){
          hvg_alleles[[mm_locus]]<-append(hvg_alleles[[mm_locus]], j)
          hvg[[mm_locus]]<-hvg[[mm_locus]]+1
          lgr$info(paste('HVG MM-', j, sep = ''))
        }
      }
    }
    
    gvh_total<-sum(gvh)
    hvg_total<-sum(hvg)
    
    mismatches<-0
    
    #get the greater of the two mismatches for each locus and count # of mismatches
    #to subtract from the total
    for(a in 1:length(gvh)){
      mismatches<-mismatches+max(gvh[[a]], hvg[[a]])
    }
    
    matches<-total-mismatches
    
    return(list(c(matches, total, gvh_total, hvg_total), hvg_alleles, gvh_alleles, missing_message))
    
  }
}

#calculate matches for DQ and DP
calcDQDP<-function(cat, d_hla, r_hla, synqList, filter_d, filter_r){
  
  group<-paste(cat, c('A1', 'B1'), sep = "")
  
  d_cat<-d_hla %>%
    select(all_of(group)) 
  
  #return 0 if selected category has not been sequenced in donor genotype
  if(all(is.na(d_cat))){
    return(list(c(0,0,0,0), NULL, NULL, NULL))
  }
  
  r_cat<-r_hla %>%
    select(all_of(group)) 
  
  #return 0 if selected category for recipient only has NAs (i.e two non-expressed
  #alleles)
  if(all(is.na(r_cat))){
    return(list(c(0,0,0,0), NULL, NULL, NULL))
  }
  
  d_alleles<-na.omit(unlist(d_cat, use.names=F))
  r_alleles<-na.omit(unlist(r_cat, use.names=F))
  
  #novel allele handling - recipient
  if(any(grepl('@', r_alleles))){
    r_alleles<-processNovelAlleles(r_alleles, synqList)
  }
  
  #novel allele handling - donor
  if(any(grepl('@', d_alleles))){
    d_alleles<-processNovelAlleles(d_alleles, synqList)
  }
  
  if(!is.null(filter_d)){
    d_filtered_alleles<-d_alleles[!d_alleles %in% filter_d]
  } else{
    d_filtered_alleles<-d_alleles
  }
  
  if(!is.null(filter_r)){
    r_filtered_alleles<-r_alleles[!r_alleles %in% filter_r]
  } else{
    r_filtered_alleles<-r_alleles
  }
  
  total<-matches<-length(d_cat[!is.na(d_cat)])
  
  if(total==4){
    total<-matches<-total/2
  }
  
  gvh<-hvg<-sapply(group, function(x) 0)
  
  if(all(d_alleles %in% r_alleles) & all(r_alleles %in% d_alleles)){
    gvh<-hvg<-sapply(group, function(x) NULL)
    return(list(c(total, matches, 0, 0), hvg, gvh, NULL))
  } else{
    
    nmdp_flag<-reg_flag<-FALSE
    
    d_mm_alleles<-unique(d_filtered_alleles[which(!d_filtered_alleles %in% r_filtered_alleles)])
    r_mm_alleles<-unique(r_filtered_alleles[which(!r_filtered_alleles %in% d_filtered_alleles)])
    
    all_nmdp_alleles<-d_mm_alleles[isNMDP(d_mm_alleles)]
    nmdp_translated<-sapply(all_nmdp_alleles, function(x) NULL)
    
    if(length(nmdp_translated)!=0){
      for(n in 1:length(nmdp_translated)){
        nmdpRes<-convertNMDP(names(nmdp_translated)[[n]])
        nmdp_translated[[n]]<-paste(gsub('^(.*?)\\*.*$', '\\1', names(nmdp_translated)[[n]]), nmdpRes, sep="*")
      }
    }
    
    missing_message<-c()
    gvh_alleles<-sapply(group, function(x) c())
    
    ##GvH calculation
    for(i in r_mm_alleles){
      
      mm_locus<-gsub('^(.*?)\\*.*$', '\\1', i)
      
      d_locus_alleles<-d_filtered_alleles[grepl(mm_locus,gsub('^(.*?)\\*.*$', '\\1', d_filtered_alleles))]
      
      nmdp_alleles<-names(nmdp_translated)[grepl(mm_locus,gsub('^(.*?)\\*.*$', '\\1', names(nmdp_translated)))]
      reg_alleles<-setdiff(d_locus_alleles, nmdp_alleles)
      
      if(length(nmdp_alleles)!=0){
        nmdp_flag<-!i %in% unlist(nmdp_translated[c(nmdp_alleles)])
      }
      
      if(length(reg_alleles)!=0){
        
        align<-alignments[[mm_locus]]
        
        r_gvh_prot<-align %>%
          filter(trimmed_allele == i) %>%
          slice(1) %>%
          select(as.character(seq(1, 90))) %>%
          paste(collapse = "")
        
        if(grepl('*', r_gvh_prot, fixed = TRUE)){
          missing_message<-append(paste('RECIPIENT-', i, sep =''), missing_message)
        }
        
        d_gvh_prot<-align %>%
          filter(trimmed_allele %in% c(d_locus_alleles)) %>%
          distinct(trimmed_allele, .keep_all = TRUE) %>%
          select(c(trimmed_allele,as.character(seq(1, 90)))) %>%
          unite(.,'all', -trimmed_allele, sep='')
        
        reg_flag<-!any(r_gvh_prot == d_gvh_prot$all)
        
        if(!reg_flag){
          filter<-d_gvh_prot$trimmed_allele[r_gvh_prot == d_gvh_prot$all]
          d_mm_alleles<-d_mm_alleles[!d_mm_alleles %in% filter]
        }
      }
      
      if(reg_flag | nmdp_flag){
        #if donor alleles consist of 1 nmdp and 1 regular allele, assess
        #if the recipient allele has FALSE for reg_flag or nmdp_flag -- if yes,
        #move onto next iteration without adding to the counter bc there is a match
        if(length(reg_alleles)==1 & length(nmdp_alleles)==1){
          if(!reg_flag | !nmdp_flag){
            next 
          }
        }
        gvh_alleles[[mm_locus]]<-append(gvh_alleles[[mm_locus]], i)
        gvh[[mm_locus]]<-gvh[[mm_locus]]+1
        lgr$info(paste('GvH MM-', i, sep = ''))
        
      }
    }
    
    hvg_alleles<-sapply(group, function(x) c())
    
    for(j in d_mm_alleles){
      
      mm_locus<-gsub('^(.*?)\\*.*$', '\\1', j)
      
      r_locus_alleles<-r_filtered_alleles[grepl(mm_locus,gsub('^(.*?)\\*.*$', '\\1', r_filtered_alleles))]
      all_r_locus_alleles<-r_alleles[grepl(mm_locus,gsub('^(.*?)\\*.*$', '\\1', r_alleles))]
      
      if(isNMDP(j)){
        if(!any(nmdp_translated[[j]] %in% all_r_locus_alleles)){
          hvg_alleles[[mm_locus]]<-append(hvg_alleles[[mm_locus]], j)
          hvg[[mm_locus]]<-hvg[[mm_locus]]+1
          lgr$info(paste('HvG MM-', j, sep = ''))
        }
        
      } else{
        
        align<-alignments[[mm_locus]]
        
        d_hvg_prot<-align %>%
          filter(trimmed_allele == j) %>%
          slice(1) %>%
          select(as.character(seq(1, 90))) %>%
          paste(collapse = "")
        
        if(grepl('*', d_hvg_prot, fixed = TRUE)){
          missing_message<-append(paste('DONOR-', j, sep =''), missing_message)
        }
        
        r_hvg_prot<-align %>%
          filter(trimmed_allele %in% c(all_r_locus_alleles)) %>%
          distinct(trimmed_allele, .keep_all = TRUE) %>%
          select(as.character(seq(1, 90))) %>%
          unite(., 'all', sep='') %>%
          select('all')
        
        if(!any(d_hvg_prot == r_hvg_prot)){
          hvg_alleles[[mm_locus]]<-append(hvg_alleles[[mm_locus]], j)
          hvg[[mm_locus]]<-hvg[[mm_locus]]+1
          lgr$info(paste('HvG MM-', j, sep = ''))
        } 
      }
    }
    
    gvh_total<-max(gvh)
    hvg_total<-max(hvg)
    
    matches<-total-max(gvh_total, hvg_total)
    
    return(list(c(matches, total, gvh_total, hvg_total), hvg_alleles, gvh_alleles, missing_message))
  }
}

#get sample number  
getSampleNumber<-function(con, r_itl){
  
  res<-dbGetQuery(con, sprintf('SELECT sample_number 
                                       FROM dbo.Match_grades 
                                       WHERE recipient_number = %s and donor_number = %s', r_itl, r_itl))$sample_number
  
  return(res)
  
}

#get test numbers for IgG where tests are billable
getIgGTestNums<-function(con, s_num){
  
  res<-dbGetQuery(con, sprintf("SELECT test_number
                           FROM dbo.Patient_tests
                           WHERE sample_number = %s AND test_type_code in ('LSAB1', 'LSAB2') AND billing_flag = 'B'", s_num))
  testNums<-res %>%
    pull(test_number)
  
  return(testNums)
  
}

#get antibody screening results
getAbResults <- function(con, testNumbers) {
  
  res<-dbGetQuery(con, sprintf('SELECT called_antibodies
                           FROM dbo.Screening_results
                           WHERE test_number in (%s)', testNumbers))
  
  #use str_trim to make sure there is no white space before or after if 'Negative'
  #is present in the called_antibodies column 
  res$called_antibodies<-str_trim(res$called_antibodies)
  
  return(res)
}

#get MFI values
getMFIvals<-function(con, p_itl, testNums){
  
  res<-dbGetQuery(con, sprintf('SELECT antigen, probe_id, allele, average_value
                                FROM dbo.Luminex_SA_bead_detail
                                WHERE patient_number = %s AND
                                test_number in (%s)', p_itl, testNums))
  
  #replace ag with specific antigen group for distinct alleles that are called
  #differently
  agSpecific<- list(
    'A*02:03' = 'A203',
    'A*02:10' = 'A210',
    'A*24:03' = 'A2403',
    'B*27:08' = 'B2708',
    'B*39:02'= 'B3902',
    'B*40:05'= 'B4005',
    'B*51:02' = 'B5102',
    'DRB1*14:03' = 'DR1403',
    'DRB1*14:04' = 'DR1404',
    'DPA1*01:03,DPB1*02:01' = 'DPB1*02:01',
    'DPA1*01:03,DPB1*02:02' = 'DPB1*02:02',
    'DPA1*01:03,DPB1*04:01' = 'DPB1*04:01',
    'DPA1*01:03,DPB1*04:02' = 'DPB1*04:02'
  )

  res<-res %>%
    mutate(antigen = recode(allele, !!!agSpecific, .default = antigen)) %>%
    filter(probe_id != '')

  return(res)
  
}

#get antigen table
getAntigenTable<-function(con){
  
  res<-dbGetQuery(con, "SELECT antigen_name, locus, sero_eq
                FROM   dbo.Antigens")
  
  res<-data.frame(sapply(res[1:ncol(res)], function(x) str_trim(x)))
  
  res$antigen_name<-sapply(res$antigen_name, function(x) BIGDAWG::GetField(x, 2))
  
  res<-res %>% 
    filter(!is.na(sero_eq) & !sero_eq %in% c('None', 'Unknown'))
  
  res<-res %>%
    distinct(antigen_name, .keep_all = TRUE)
  
  return(res)
}

#get DPDQ serological reference file
#made manually from 'serologically indistinguishable alleles' file
getC2RefAlleles<-function(){
  
  c2_sero_alleles<-readxl::read_excel('ref/DPDQ_serological_reference.xlsx')
  
  return(c2_sero_alleles)
}

#determine if allele is NMDP
isNMDP<-function(allele){
  return(grepl('[A-Z]', substr(gsub('.*?:', "", allele), 1,1)))
  
}

#calculate if DSA is Y or N
calcDSA<-function(db_con, mismatched_alleles, called_antibodies, mfi_vals, donorTyping, recipTyping){

  call_dsa<-'N'
  
  dsaDF<-data.frame(
    allele = character(), 
    surrogate = character(), 
    alleles = character(), 
    mfi = character()
  )
  
  #DSA = No if no mismatched alleles or called_antibodies = Negative
  if(length(mismatched_alleles)==0 | any(called_antibodies == 'Negative')){
    return(list(call_dsa, dsaDF))
  }
  
  #create vector for mismatched_alleles in case there are surrogates used for C2
  #surrogate alleles need to be appended to mismatched alleles while evaluating
  #if both subunits are present in heterodimers 
  mmAllelesAppend<-donorTyping
  
  for(t in mismatched_alleles){

    surrogate<-NA
    nmdp_allele<-NULL
    
    lgr$info(sprintf('Evaluating HvG mismatched allele %s...', t))
    
    heterodimerSurrogate = FALSE 
    
    #if any novel alleles, 
    if(grepl('@', t)){
      t<-gsub('@', '', t)
    }
    
    #skip if 'Q' suffix and not NMDP allele 
    if(str_sub(t, -1) =='Q' & !isNMDP(t)){
      lgr$info(sprintf('Skipping %s for DSA analysis', t))
      next
    }

    #replace suffix if low expression or protein level ambiguity
    if(grepl('L|P', str_sub(t, -1)) & !isNMDP(t)){
      t<-substr(t,1, nchar(t)-1)
    }
    
    locus<-gsub('^(.*?)\\*.*$', '\\1', t)
    
    #nmdp conversion; change iterator to subtype values
    if(isNMDP(t)){
      
      nmdp_allele<-t
      convertedAllele<-convertNMDP(t)[[1]]
      #use first subtype, which is the most common one
      t<-paste(locus, convertedAllele[[1]], sep='*')
    }
    
    #if heterodimer locus, find alleles that contain the mismatched allele
    if(locus %in% c('DPA1', 'DPB1', 'DQA1', 'DQB1')){
      allele_mfi<-mfi_vals %>%
        filter(grepl(t, allele, fixed = TRUE)) %>%
        filter(grepl(',', allele))
      
    } else{
      allele_mfi<-mfi_vals %>%
        filter(allele %in% t) %>%
        distinct(allele, .keep_all = T)
    }
    
    aspBeads<-keep(mfi_vals$allele, ~ !is.na(.x) & .x != "")
    
    #if allele is not tested by ab screening and is A, B, C, DRB1, DRB3/4/5, use 
    #antigen table to find serological equivalent
    if(nrow(allele_mfi)==0){
      
      lgr$info(sprintf('Using a serological surrogate for %s', paste(t, collapse = ', ')))
      
      if(locus %in% c('A', 'B', 'C', 'DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQB1')){
        
        search<-sprintf('%s_serotype', locus)
        if(locus %in% c('DRB1','DRB3','DRB4','DRB5')){
          search<-sprintf('%s_serotype', 'DRB1345')
        }
        sero_table<-readRDS(list.files('ref/', pattern = search, full.names = T))
        
        surrogateQuery<-sero_table %>%
          filter(Allele == t)
        
        surrogate<-surrogateQuery %>%
          select(Serotype) %>% 
          pull()
        
        #if allele is not found in KO table
        if(length(surrogate)==0){
          lgr$info(sprintf('No surrogate found for %s', t))
          next
        }
        
        #if asp, check if the bead is tested in ab testing
        if(grepl(':', surrogate)){
          #if asp surrogate is not in beads, check the sero_table again to see if 
          #the surrogate bead was used as a surrogate for a bead that is present in
          #ab testing
          #ie surrogate for A*02:14 is A*02:02, but we don't have A*02:02. however, A*02:02
          #is also a surrogate for A*02:05, which we do have
          
          if(!surrogate %in% aspBeads){
          
            aspSearch<-sero_table %>%
              filter(Serotype == surrogate) %>%
              filter(Allele %in% aspBeads)
            
            #if there is still no ASP surrogate, then use the ag from the ImputedAntigen col
            if(nrow(aspSearch)==0){
              surrogate<-surrogateQuery %>%
                select(ImputedAntigen) %>%
                pull()
              lgr$info(sprintf('Unable to find an allele specific surrogate for %s... using imputed antigen %s as surrogate', t, surrogate))
            } 
          
            if(nrow(aspSearch)!=0){ 
              surrogate<-aspSearch %>% pull(Allele)
            } 
          }
        }
        
        #DRB35 and serotype is not antigen
        if(locus %in% c('DRB3', 'DRB5') & grepl('DR-5', surrogate)){
            DRmapping<-list(
              #DRB3
              'DR-5201' = 'DRB3*01:01',
              'DR-5202' = c('DRB3*02:01', 'DRB3*01:02'), 
              'DR-5203' = 'DRB3*03:01',
              #DRB5
              'DR-5101' = c('DRB5*01:01', 'DRB5*01:02'), 
              'DR-5102' = 'DRB5*02:02',
              #use ag for DR-5103 since we don't have the bead (DRB5*01:03)
              'DR-5103' = 'DR51'
            )
            surrogate<-DRmapping[[surrogate]]
            #in case IgG test does not include new ExPlex beads (i.e. DRB3*01:02)
            surrogate<-surrogate[surrogate %in% aspBeads]
        }
        
        #use any(grepl()), just in case multiple surrogates are found for ASP beads
        allele_mfi<-mfi_vals %>%
          filter(if(any(grepl(':', surrogate))) allele %in% surrogate else antigen == surrogate)
        
        if(length(surrogate)==0){
          lgr$info(sprintf('No surrogate found for %s', t))
          surrogate<-NA
          next
        } else {
          lgr$info(sprintf('Surrogate used: %s', paste(surrogate, collapse = ', ')))
        }
        
        #alleles in this category should not be used as surrogates for their respective
        #antigen groups because they are serologically distinct
        #filter these alleles out
        if(length(surrogate) == 1){
          if(surrogate %in% c('B35', 'B62', 'Cw16')){
            filter_map <- list(
              "B35" = c('B*35:12', 'B*35:02'),
              "B62" = c('B*15:24'),
              "Cw16" = c('C*16:02')
            )
            if (surrogate %in% names(filter_map)) {
              allele_mfi <- allele_mfi %>% filter(!allele %in% filter_map[[surrogate]])
              }
            } 
        }
      } else if(locus %in% c('DQA1', 'DPA1', 'DPB1')){
        
        heterodimerSurrogate = TRUE
        
        surrogate<-c2_antigen_ref %>%
          filter(allele == t) %>%
          select(panelbead) %>%
          pull()
        
        #if bead is in DPDQ_serological_reference.xlsx, use serologically identical
        #bead
        if(length(surrogate)!=0){
          allele_mfi<-mfi_vals %>%
            filter(allele==surrogate)
        } else{
          #else, check Kazu's table for serotype surrogate, or direct allele surrogate
          sero_table<-readRDS(list.files('ref/', pattern = sprintf('%s_serotype', locus), full.names = TRUE))
          
          surrogate<-sero_table %>%
            filter(Allele == t) %>%
            select(Serotype) %>%
            pull()
         
          if(locus == 'DPB1'){
            if(length(surrogate)!=0){
              allele_mfi<-mfi_vals %>%
                filter(antigen == surrogate)
            } else{
              lgr$info(sprintf('No surrogate found for %s', t))
              surrogate<-NA
              next
            }
          } else{
            #DQA1, DPA1
            #can't use fixed = TRUE with later grepl with '|', so have to 
            #escape asterisks
            surrogate<-sero_table %>% 
              filter(Serotype == surrogate) %>%
              mutate(Allele = gsub('*', '\\*', Allele, fixed = TRUE)) %>%
              pull(Allele)
            
            if(length(surrogate)!=0){
              surrogate<-paste(surrogate, collapse = '|')
              
              #get serologically equivalent alpha alleles in beads
              allele_mfi<-mfi_vals %>%
                filter(grepl(surrogate, allele))
            } else{
              lgr$info(sprintf('No surrogate found for %s', t))
              surrogate<-NA
              next
            }
          }
        }
        lgr$info(sprintf('Surrogate used: %s', gsub('\\\\', '', surrogate)))
      }
    }  
    
    #evaluate if allele(alleles if there are multiple serological surrogates) 
    #have an average mfi value greater than 2000
    mfi_eval<-allele_mfi %>%
      mutate(bool=average_value>2000) 
    
    #if DQ or DP, check if both subunits are present in the donor's phenotype
    if(locus %in% c('DQA1', 'DQB1', 'DPA1', 'DPB1')){
      
      possibleDSA<-mfi_eval %>%
        filter(bool == TRUE) %>%
        pull(allele)
      
      if(length(possibleDSA) != 0){
        #if a heterodimer surrogate was used, add subunits to the donor's phenotype
        if(heterodimerSurrogate){
          #beta
          if(grepl('B', locus)){ 
            mmAllelesAppend<-append(mmAllelesAppend, unique(unlist(lapply(strsplit(allele_mfi$allele, ','), '[[', 2))))
          } else{
            #alpha
            mmAllelesAppend<-append(mmAllelesAppend, unique(unlist(lapply(strsplit(allele_mfi$allele, ','), '[[', 1))))
          }
        }
        for(i in 1:length(possibleDSA)){
          heterodimerSplit<-strsplit(possibleDSA[[i]], ',')[[1]]
          
          #check if both subunits are present in donor's phenotype; if not, remove 
          #from list
          if(!all(heterodimerSplit %in% mmAllelesAppend)){
            mfi_eval<-mfi_eval %>%
              filter(allele != possibleDSA[[i]])
          }
        } 
      }
    }
    
    mfi_bool<-mfi_eval%>%
      filter(!is.na(allele)) %>%
      filter(bool == TRUE)
    
    if(any(mfi_bool$bool)){ 
      if(!is.null(nmdp_allele)){
        t<-nmdp_allele  
      }
      
      #check if antigen with MFI > 2000 is in called_antibodies list
      #for surrogates, there can be multiple probe_ids
      if(any(unique(mfi_eval$antigen) %in% called_antibodies) | any(mfi_eval$probe_id %in% called_antibodies)){
        
        call_dsa<-'Y'
        
        #if there is one or more bead > 2000 MFI
        if(nrow(mfi_bool)>=1){
          #if surrogate was used, multiple beads can be present; use the min and max of all bead data
          if(nrow(mfi_eval) > 1){
            mfi_value<-sprintf('%s to %s',min(mfi_eval$average_value),max(mfi_eval$average_value))
          } else{
            mfi_value<-mfi_eval$average_value
          }
        }
        
        mfi_mess<-sprintf('%s (MFI = %s)', t, mfi_value)
        lgr$info(paste('DSA = Y due to ', mfi_mess, sep=''))
        
        dsaDF<-rbind(dsaDF, data.frame(allele = t, surrogate = surrogate, alleles = paste(mfi_bool$allele, collapse = ', '), mfi = mfi_value))
        
      } else{
        #may need to modify later due to mm allele not being found in positive abs
        #due to alpha, beta, or combo asp
        lgr$info(sprintf('False positive reported for %s ', t))
        
        #if call_dsa is not already Y, then change to N
        if(call_dsa != 'Y'){
          call_dsa<-'N'
        }
      }
    }
  }
  
  return(list(call_dsa, dsaDF))
}

#determine if mutation for novel position is in mature protein sequence
determinePosition<-function(position, locus){
  
  #class 1
  if(locus %in% c('A', 'B', 'C')){
    if(position %in% seq(1, 182)){
      return(TRUE)
    } else{
      return(FALSE)
    }
    #class 2
  } else{
    if(position %in% seq(1, 90)){
      return(TRUE)
    } else{
      return(FALSE)
    }
  }
}

#process novel alleles
processNovelAlleles<-function(novelAlleles, synqList){
  
  novel_handle<-novelAlleles[grepl('@', novelAlleles)]
  
  for(i in novel_handle){
    
    novelFull<-paste('HLA-', i, sep ='')
    
    #if the mutation is in the ARD and is SYNONYMOUS ('y'),
    #move onto the next novel allele and remove @ sign
    
    #if the mutation is in the ARD and is NON-SYNONYMOUS ('n'), keep @ sign
    
    #if the mutation is not in the ARD (i.e. it's not in the syn-nonsyn list), 
    #move onto the next allele and remove @ sign from novel allele
    if(novelFull %in% names(synqList)){
      if(synqList[[novelFull]] == 'y'){
        browser()
        novelAlleles[novelAlleles == i]<-gsub('@', '', i)
        next
      } else{
        next
      }
    } else{
      novelAlleles[novelAlleles == i]<-gsub('@', '', i)
      next
    }
  }
  return(novelAlleles)
}

#process phenotypes with the same P group or null
processGrouped<-function(alleles){
  
  same_group<-c()
  pAlleles<-alleles[grepl(',', alleles)]
  
  for(i in 1:length(pAlleles)){
    split<-strsplit(pAlleles[[i]], ',')[[1]]
    allele1<-split[[1]]
    allele2<-split[[2]]
    
    locus<-gsub('^(.*?)\\*.*$', '\\1', allele1)
    delim = '*'
    if(locus %in% c('DRB3', 'DRB4', 'DRB5')){
      locus<-'DRB'
      delim <- ''
    }
    
    fullAllele2<-paste(locus, allele2, sep = delim)
    alleles<-append(fullAllele2, alleles)
    alleles<-alleles[!grepl(',', alleles)]
    same_group<-append(fullAllele2, same_group)
  }
  
  return(list(alleles, same_group))
}

mismatchFormatter<-function(locus, gvh_alleles, hvg_alleles){
  
  gvh_alleles<-paste(gvh_alleles, collapse = ', ')
  
  if(gvh_alleles == ''){
    gvh_alleles<-'None'
  } 
  
  hvg_alleles<-paste(hvg_alleles, collapse = ', ')
  
  if(hvg_alleles == ''){
    hvg_alleles<-'None'
  } 
  
  return(paste(
    '<h3 style="font-weight: bold;"><span style="background-color: yellow;">', locus, '</span></h3>',
    '<b><u>GvH</b></u>: ', gvh_alleles, '<br>',
    '<b><u>HvG</b></u>: ', hvg_alleles,
    sep = ''
  ))
}

#generate module ID
generateID<-function(count){
  module_id<-paste0("module_", count)
}

#pop up modal for missing sequences
missingModal<-function(message){
  modalDialog(
    HTML(message),
    footer = modalButton('OK!')
  )
}

#get maintainer e-mail
getMaintainerEmail <- function(){
  
  maintainerEmail <- Sys.getenv('MAINTAINER_EMAIL')
  
  if(maintainerEmail == ""){
    maintainerEmail<-'maintaneremail@placeholder.com'
  }
  
  return(maintainerEmail)
}
