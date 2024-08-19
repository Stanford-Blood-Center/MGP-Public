#external functions
#v 1.3.1

suppressPackageStartupMessages(library(odbc))
suppressPackageStartupMessages(library(tidyverse))
library(stringr) 
library(httr)
suppressPackageStartupMessages(library(rvest))

dbConn <- function(){
  
  con <- dbConnect(odbc(),
                   Driver = "SQL Server",
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

#load NMDP file downloaded from https://hml.nmdp.org/MacUI/
#in the future, find a way to automate downloading new releases of the NMDP file 
#and updating the version in /ref
loadNMDP<-function(){
  
  nmdp_file<-readRDS(list.files('ref', pattern = 'nmdp', full.names = T))
  
  return(nmdp_file)
}

#extract available donors for evaluation
getDonors<-function(mg_df){
  donors <- mg_df %>%
    filter(recipient_number != donor_number) %>%
    select(donor_number) %>%
    pull()
  
  return(donors)
}

#filter to specific donor
filterDonor<-function(mg_df, d_itl){
  d_all<-mg_df %>%
    filter(donor_number == d_itl) 
  return(d_all)
}

#get HLA typing
#transform loci columns into rows, populate homozygous alleles
getTyping<-function(input_df, type){
  hla_cols <- sort(do.call(paste, c(expand.grid(c('a', 'b', 'c', 'dr', 'drp', 'dqa', 'dqb', 'dpa', 'dpb'), c(1,2)), sep = '_')))
  
  t_df <- input_df %>%
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
  
  r_null_allele<-NULL
  
  if(type == 'r'){
    if(length(null_alleles)!=0){
      if(null_alleles == 'DRB'){
        delim = ''
      } else{
        delim = '*'
      }
      r_null_allele<-paste(null_alleles, t_df[,null_alleles][str_sub(t_df[,null_alleles], -1) == 'N'], sep = delim)
    }
  }
  
  #null allele scenarios
  for(a in null_alleles){
    
    #if both alleles are null, replace both with NA
    if(str_sub(t_df[,a][1], -1) == 'N' & str_sub(t_df[,a][2], -1) == 'N'){
      t_df[,a][1]<-NA
      t_df[,a][2]<-NA
      lgr$info('Two non-expressed alleles detected')
      next
    }
    
    #if one non-expression allele is present, replace it with the other
    #expressed allele 
    if(str_sub(t_df[,a][1], -1) == 'N'){
      t_df[,a][1]<-t_df[,a][2]
    } else{
      t_df[,a][2]<-t_df[,a][1]
    }
    lgr$info('One non-expressed allele detected')
  }
  
  #if a locus has heterozygous alleles and does not have NMDP alleles, 
  #check if they have the same P group
  residue<-NULL
  for(j in colnames(t_df)){
    if(!any(grepl(':[A-WYZ]', t_df[,j])) & !all(is.na(t_df[,j]))){
      
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
        #homozygous. replace the second allele with the first one 
        if(s_pg_alleles){
          lgr$info(paste(align_locus, 'has 2 alleles with the same P-group. Replacing', t_df[2,j], 'with', t_df[1,j]))
          t_df[2,j]<-t_df[1,j]
        }
      }
    }
  }
  
  t_df<-t_df %>%
    rowwise() %>%
    mutate(across(everything(), ~ifelse(!is.na(.), ifelse(cur_column()=="DRB", paste(cur_column(), ., sep = ""), paste(cur_column(), ., sep = "*")), .)))
  
  return(list(t_df, r_null_allele))
}

#convert NMDP code to subtype 
convertNMDP<-function(n_alleles){
  
  f1<-gsub('.*?(.*?)\\*', "\\1", sapply(strsplit(n_alleles, ':'), '[[', 1))
  f2<-sapply(strsplit(n_alleles, ':'), '[[', 2) 
  
  donor_st<-sapply(f1, function(x) NULL)
  nonexistentNMDP<-c()
  
  for(i in 1:length(f2)){
    
    nmdpFilter<-nmdp_file %>%
      filter(code %in% c(f2[[i]])) %>%
      pull() %>%
      strsplit('/') %>%
      unlist()
    
    if(!is.null(nmdpFilter)){
      donor_st[[i]]<-nmdpFilter
      
      #paste first field to subtype if not present in NDMP codelist
      if(any(!grepl(':', donor_st[[i]]))){
        donor_st[[i]]<-paste(names(donor_st)[[i]], donor_st[[i]], sep =":")
      }
    } else{
      lgr$info(paste(f2[[i]], 'was not found in the NMDP reference file'))
      nonexistentNMDP<-append(f2[[i]], nonexistentNMDP)
    }
  }
  
  if(!is.null(nonexistentNMDP)){
    stop('NMDP alleles not found in the NMDP reference file')
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
  
  legend<-list('Non-Permissive'='N', 'Permissive'='P','Unknown'='U')
  
  if(length(union(setdiff(r_dpb, d_dpb), setdiff(d_dpb, r_dpb)))==0){
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
calcABCDRB<-function(cat, d_hla, r_hla){
  
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
    return(list(c(0,0,0,0), NULL))
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
    return(list(c(0,0,0,0), NULL))
  }
  
  d_alleles<-na.omit(unlist(d_cat, use.names=F))
  r_alleles<-na.omit(unlist(r_cat, use.names=F))
  
  total<-matches<-length(d_alleles)
  
  gvh<-hvg<-sapply(group, function(x) 0)
  
  if(all(d_alleles %in% r_alleles) & all(r_alleles %in% d_alleles)){
    return(list(c(total, matches, 0, 0), NULL))
    
  } else{
    
    d_mm_alleles<-unique(d_alleles[which(!d_alleles %in% r_alleles)])
    r_mm_alleles<-unique(r_alleles[which(!r_alleles %in% d_alleles)])
    nmdp_flag<-reg_flag<-FALSE
    
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
      
      d_locus_alleles<-d_cat %>%
        select(all_of(mm_locus)) %>%
        pull()
      
      if(align_locus %in% c('DRB3', 'DRB4', 'DRB5')){
        #if DRB locus for recipient mismatch allele is not present in donor locus
        #alleles, add 1 to GvH and move onto next iteration
        if(all(align_locus != gsub('^(.*?)\\*.*$', '\\1', d_locus_alleles)) & !all(d_locus_alleles == '')){
          gvh[[mm_locus]]<-gvh[[mm_locus]]+1
          next
        }
        if(all(d_locus_alleles == '')){
          gvh[[mm_locus]]<-length(r_mm_alleles)
          next
        }
      } 
      
      #split up into nmdp and regular allele assessment, since some donors
      #can have nmdp and regular sequenced alleles
      nmdp_alleles<-d_locus_alleles[grepl('[A-Z]', substr(gsub('.*?:', "", d_locus_alleles), 1,1))]
      reg_alleles<-setdiff(d_locus_alleles, nmdp_alleles)
      
      if(length(nmdp_alleles)!=0){
        
        d_nmdp<-convertNMDP(nmdp_alleles)
        nmdp_flag<-!i %in% paste(align_locus, d_nmdp, sep="*")
        
      }
      
      if(length(reg_alleles)!=0){
        
        align<-alignments[[align_locus]]
        r_gvh_prot<-align %>%
          filter(trimmed_allele == i) %>%
          slice(1) %>%
          select(as.character(seq(1, residue))) %>%
          paste(collapse = "")
        
        d_gvh_prot<-align %>%
          filter(trimmed_allele %in% c(d_locus_alleles)) %>%
          distinct(trimmed_allele, .keep_all = TRUE) %>%
          select(as.character(seq(1, residue))) %>%
          unite(., 'all', sep='') %>%
          select('all')
        
        reg_flag<-!any(r_gvh_prot == d_gvh_prot)
      }
      
      if(reg_flag | nmdp_flag){
        #if donor alleles consist of 1 nmdp and 1 regular allele, assess
        #if the recipient allele has FALSE for reg_flag or nmdp_flag move on to
        #move onto next iteration without adding to the counter bc there is a match
        if(length(reg_alleles)==1 & length(nmdp_alleles)==1){
          if(!reg_flag | !nmdp_flag){
            next 
          }
        }
        gvh[[mm_locus]]<-gvh[[mm_locus]]+1
      }
    }
    
    hvg_alleles<-NULL
    
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
      
      r_locus_alleles<-r_cat %>%
        select(all_of(mm_locus)) %>%
        pull()
      
      if(align_locus %in% c('DRB3', 'DRB4', 'DRB5')){
        if(all(r_locus_alleles == '')){
          hvg_alleles<-append(hvg_alleles, j)
          hvg[[mm_locus]]<-length(d_mm_alleles)
          next
        }
        if(all(align_locus != gsub('^(.*?)\\*.*$', '\\1', r_locus_alleles)) & !all(r_locus_alleles == '')){
          hvg_alleles<-append(hvg_alleles, j)
          hvg[[mm_locus]]<-hvg[[mm_locus]]+1
          next
        }
      }
      
      if(grepl('[A-Z]', substr(gsub('.*?:', "", j), 1,1))){
        
        d_nmdp<-convertNMDP(j)
        
        if(!any(paste(align_locus, d_nmdp, sep="*") %in% r_locus_alleles)){
          hvg_alleles<-append(hvg_alleles, j)
          hvg[[mm_locus]]<-hvg[[mm_locus]]+1
        }
      } else{
        
        align<-alignments[[align_locus]]
        
        d_hvg_prot<-align %>%
          filter(trimmed_allele == j) %>%
          slice(1) %>%
          select(as.character(seq(1, residue))) %>%
          paste(collapse = "")
        
        
        r_hvg_prot<-align %>%
          filter(trimmed_allele %in% c(r_locus_alleles)) %>%
          distinct(trimmed_allele, .keep_all = TRUE) %>%
          select(as.character(seq(1, residue))) %>%
          unite(., 'all', sep='') %>%
          select('all')
        
        if(!any(d_hvg_prot == r_hvg_prot)){
          hvg_alleles<-append(hvg_alleles, j)
          hvg[[mm_locus]]<-hvg[[mm_locus]]+1
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
    
    return(list(c(matches, total, gvh_total, hvg_total), hvg_alleles))
    
  }
}

#calculate matches for DQ and DP
calcDQDP<-function(cat, d_hla, r_hla){
  
  group<-paste(cat, c('A1', 'B1'), sep = "")
  
  d_cat<-d_hla %>%
    select(all_of(group)) 
  
  #return 0 if selected category has not been sequenced in donor genotype
  if(all(is.na(d_cat))){
    return(list(c(0,0,0,0), NULL))
  }
  
  r_cat<-r_hla %>%
    select(all_of(group)) 
  
  #return 0 if selected category for recipient only has NAs (i.e two non-expressed
  #alleles)
  if(all(is.na(r_cat))){
    return(list(c(0,0,0,0), NULL))
  }
  
  d_alleles<-na.omit(unlist(d_cat, use.names=F))
  r_alleles<-na.omit(unlist(r_cat, use.names=F))
  
  total<-matches<-length(d_cat[!is.na(d_cat)])
  
  if(total==4){
    total<-matches<-total/2
  }
  
  gvh<-hvg<-sapply(group, function(x) 0)
  
  if(all(d_alleles %in% r_alleles) & all(r_alleles %in% d_alleles)){
    return(list(c(total, matches, 0, 0), NULL))
  } else{
    
    d_mm_alleles<-unique(d_alleles[which(!d_alleles %in% r_alleles)])
    r_mm_alleles<-unique(r_alleles[which(!r_alleles %in% d_alleles)])
    nmdp_flag<-reg_flag<-FALSE
    
    ##GvH calculation
    for(i in r_mm_alleles){
      
      mm_locus<-gsub('^(.*?)\\*.*$', '\\1', i)
      
      d_locus_alleles<-d_cat %>%
        select(all_of(mm_locus)) %>%
        pull()
      
      nmdp_alleles<-d_locus_alleles[grepl('[A-Z]', substr(gsub('.*?:', "", d_locus_alleles), 1,1))]
      reg_alleles<-setdiff(d_locus_alleles, nmdp_alleles)
      
      if(length(nmdp_alleles)!=0){
        
        d_nmdp<-convertNMDP(nmdp_alleles)
        nmdp_flag<-!i %in% paste(mm_locus, d_nmdp, sep="*")
        
      }
      
      if(length(reg_alleles)!=0){
        
        align<-alignments[[mm_locus]]
        
        r_gvh_prot<-align %>%
          filter(trimmed_allele == i) %>%
          slice(1) %>%
          select(as.character(seq(1, 90))) %>%
          paste(collapse = "")
        
        d_gvh_prot<-align %>%
          filter(trimmed_allele %in% c(d_locus_alleles)) %>%
          distinct(trimmed_allele, .keep_all = TRUE) %>%
          select(as.character(seq(1, 90))) %>%
          unite(., 'all', sep='') %>%
          select('all')
        
        reg_flag<-!any(r_gvh_prot == d_gvh_prot)
        
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
        gvh[[mm_locus]]<-gvh[[mm_locus]]+1
      }
    }
    
    hvg_alleles<-NULL
    
    for(j in d_mm_alleles){
      
      mm_locus<-gsub('^(.*?)\\*.*$', '\\1', j)
      
      r_locus_alleles<-r_cat %>%
        select(all_of(mm_locus)) %>%
        pull()
      
      if(grepl('[A-Z]', substr(gsub('.*?:', "", j), 1,1))){
        
        d_nmdp<-convertNMDP(j)
        
        if(!any(paste(mm_locus, d_nmdp, sep="*") %in% r_locus_alleles)){
          hvg_alleles<-append(hvg_alleles, j)
          hvg[[mm_locus]]<-hvg[[mm_locus]]+1
        }
      } else{
        
        align<-alignments[[mm_locus]]
        
        d_hvg_prot<-align %>%
          filter(trimmed_allele == j) %>%
          slice(1) %>%
          select(as.character(seq(1, 90))) %>%
          paste(collapse = "")
        
        
        r_hvg_prot<-align %>%
          filter(trimmed_allele %in% c(r_locus_alleles)) %>%
          distinct(trimmed_allele, .keep_all = TRUE) %>%
          select(as.character(seq(1, 90))) %>%
          unite(., 'all', sep='') %>%
          select('all')
        
        if(!any(d_hvg_prot == r_hvg_prot)){
          
          hvg_alleles<-append(hvg_alleles, j)
          hvg[[mm_locus]]<-hvg[[mm_locus]]+1
        }
      }
    }
    
    gvh_total<-max(gvh)
    hvg_total<-max(hvg)
    
    matches<-total-max(gvh_total, hvg_total)
    
    return(list(c(matches, total, gvh_total, hvg_total), hvg_alleles))
  }
}

#get sample number  
getSampleNumber<-function(con, r_itl){
  
  res<-dbGetQuery(con, sprintf('SELECT sample_number 
                                       FROM dbo.Match_grades 
                                       WHERE recipient_number = %s and donor_number = %s', r_itl, r_itl))$sample_number
  
  return(res)
  
}

#get test numbers for IgG 
getIgGTestNums<-function(con, s_num){
  
  res<-dbGetQuery(con, sprintf("SELECT test_number
                           FROM dbo.Patient_tests
                           WHERE sample_number = %s and test_type_code in ('LSAB1', 'LSAB2')", s_num))
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


#calculate if DSA is Y or N
calcDSA<-function(db_con, mismatched_alleles, called_antibodies, mfi_vals){
  
  call_dsa<-'N'
  nmdp_allele<-NULL
  
  #DSA = No if no mismatched alleles or called_antibodies = Negative
  if(length(mismatched_alleles)==0 | any(called_antibodies == 'Negative')){
    return(call_dsa)
  }
  
  for(t in mismatched_alleles){
    
    #skip if 'Q' suffix
    if(str_sub(t, -1) =='Q' & str_count('Q', '[[:alpha:]]') == 1){
      lgr$info(sprintf('Skipping %s for DSA analysis', t))
      next
    }
    
    locus<-gsub('^(.*?)\\*.*$', '\\1', t)
    
    #nmdp conversion; change iterator to subtype values
    if(grepl('[A-Z]', substr(gsub('.*?:', "", t), 1,1))){
      
      nmdp_allele<-t
      
      #use first subtype, which is the most common one
      t<-paste(locus, convertNMDP(t)[[1]], sep='*')
    }
    
    allele_mfi<-mfi_vals %>%
      filter(allele %in% t) %>%
      distinct(allele, .keep_all = T)
    
    #if allele is not tested by ab screening and is A, B, C, DRB1, DRB3/4/5, use 
    #antigen table to find serological equivalent
    if(nrow(allele_mfi)==0 & locus %in% c('A', 'B', 'C', 'DRB1', 'DRB3', 'DRB4', 'DRB5', 'DQB1')){
      
      lgr$info(sprintf('Using a serological surrogate for %s', paste(t, collapse = ', ')))
      
      surrogate<-antigen_ref %>%
        filter(antigen_name == t) %>%
        select(sero_eq) %>% 
        distinct() %>%
        pull()
      
      if(length(surrogate)==0){
        lgr$info(sprintf('No surrogate found for %s', t))
        next
      }
      
      allele_mfi<-mfi_vals %>%
        filter(antigen == surrogate)
    } 
    
    if(nrow(allele_mfi)==0 & locus %in% c('DQA1', 'DPA1', 'DPB1')){
      
      lgr$info(sprintf('Using a serological surrogate for %s', paste(t, collapse = ', ')))
      
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
        sero_table<-readRDS(paste('ref/', locus, '_serotype.rda', sep=''))
        
        surrogate<-sero_table %>%
          filter(Allele == t) %>%
          select(Serotype) %>%
          pull()
        
        if(length(surrogate)!=0){
          #if ':' is present in surrogate, filter by allele column
          if(grepl(':', surrogate)){
            allele_mfi<-mfi_vals %>%
              filter(allele==surrogate)
          } else{
            #else filter by antigen group
            allele_mfi<-mfi_vals %>%
              filter(antigen==surrogate)
          }
        } else{
          lgr$info(sprintf('No surrogate found for %s', t))
          next
        }
      }
    }
    
    #evaluate if allele(alleles if there are multiple serological surrogates) 
    #have an average mfi value greater than 2000; if any do, DSA = Yes
    mfi_eval<-allele_mfi %>%
      mutate(bool=average_value>2000) 
    
    mfi_bool<-mfi_eval%>%
      pull(bool)
    
    if(any(mfi_bool)){
      if(!is.null(nmdp_allele)){
        t<-nmdp_allele  
      }
      #check if antigen with MFI > 2000 is in called_antibodies list
      #for surrogates, there can be multiple probe_ids
      if(unique(mfi_eval$antigen) %in% called_antibodies | any(mfi_eval$probe_id %in% called_antibodies)){
        call_dsa<-'Y'
        lgr$info(sprintf('Calling DSA Y due to %s', t))
        #if DSA = Y, break out of loop
        break
      } else{
        lgr$info(sprintf('False positive reported for %s ', t))
        call_dsa<-'N'
      }
    }
  }
  
  return(call_dsa)
}
