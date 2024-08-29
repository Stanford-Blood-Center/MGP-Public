#v 1.4.0

library(shiny)

executeMGPUI<- function(id, type) {
  
  ns <- NS(id)
  
  if(type == 'spinner'){
    uiOutput(ns('spinner'))
  } else{
    uiOutput(ns('log'))
  }
}


executeMGPServer<-function(id, creds, patient, donor, log, hlaRecipient, hlaRecipientNull, hlaDonor, syn_nonsyn, donorMatchGrade, clearButton) {
 
   moduleServer(
     id,
     function(input, output, session) {
      #browser()
       ns <- session$ns
       
       mgp<-reactiveValues(result=NULL, mgDBTable=NULL, mgpTable=NULL)
       final<-reactiveValues(log_path=NULL, result=NULL, df=NULL, mgpTable=NULL, mgDBTable=NULL)
        
                                                ##### OUTPUT DIFFERENCES #####
              #if technologist level, only output finish text to indicate success or failure
              #else, output comparison text, data from the Match_grades table in DB, and what MGP calculated
       
       if(creds == '30'){
        output$spinner <- renderUI({
            uiOutput(ns('finish_text'))
          })
       } else if(creds %in% c('50', '60')){
        output$spinner <-renderUI({
          tagList(uiOutput(ns('comparisonText')), uiOutput(ns('mgTable')), uiOutput(ns('mgPopulationTable')))
        })
       }
      
                                                     ##### CALCULATION #####
      
      if(creds == '30'){ 
      output$finish_text<-renderUI({
        
        run_res<-isolate({calcMatchGrade(patient, donor, log, creds, hlaRecipient, hlaRecipientNull, hlaDonor, syn_nonsyn, donorMatchGrade)})
        
        mgp$result<-run_res[[1]]
        final$log_path<-run_res[[2]]
        
        if(mgp$result){
          message<-paste('Match Grade evaluations for recipient ITL', patient, 'completed!')
          lgr$info(message)
          lgr$info('**********MATCH GRADE EVALUATION END**********')
          hide_spinner()
          message
        } else{
          message<-'Match Grade evaluations for recipient ITL %s. Please download the log and e-mail it to %s to troubleshoot.'
          lgr$info(sprintf(message, paste(patient, 'failed'), 'livtran@stanford.edu'))
          lgr$info('**********MATCH GRADE EVALUATION END**********')
          hide_spinner()
          HTML(sprintf(message, paste(patient, ' <b><span style=color:red;>failed</b></span>', sep=""), '<b>livtran@stanford.edu</b>'))
        } 
      })
    } else if(creds %in% c('50', '60')){
        
        #difference comparison output
        output$comparisonText<-renderUI({
          req(mgp$result)
          if(mgp$result == 'TRUE'){
            
            comparison<-all.equal(mgp$mgDBTable, mgp$mgpTable)
            
            baseText<-"<b>Mismatched columns: </b>"
            
            if(length(comparison)==1 & any(!grepl('Component', comparison))){
              comparison<-paste(baseText, ' None', sep = '')
            } else{
              #[\u201C\u201D] = curly quotes
              comparison<-paste(baseText, paste(gsub("[\u201C\u201D]", "", str_extract(comparison, "(?<=Component )[^:]+")), collapse = ', '), sep = '')
            }
            HTML(paste0('<div style="text-align: center;">', comparison, '</div>'))
          } else{
            NULL
          }
        })
        
        #MGP output
        output$mgTable<-renderUI({
          
          run_res<-isolate({calcMatchGrade(patient, donor, log, creds, hlaRecipient, hlaRecipientNull, hlaDonor, syn_nonsyn, donorMatchGrade)})
          
          mgp$result<-run_res[[1]]
          final$log_path<-run_res[[2]]
          df<-run_res[[3]]
          
          if(mgp$result){
            mgp$mgpTable<-df %>%
              replace(is.na(.), '') %>%
              relocate(ABCDRDQ_match, .before = ABCDRDQ_alleles) %>%
              relocate(ABCDRB1_match, .before = ABCDRB1_alleles) %>%
              relocate(DRB345DQDP_match, .before= DRB345DQDP_alleles) %>%
              relocate(Seven_loci_alleles, .before= seven_loci_match)
            
            hide_spinner()
            tagList(
              hr(),
              h3('MGP Results', style = "text-align: center; font-weight: bold;"),
              DT::renderDT(mgp$mgpTable, options = list(dom = 't', scrollX = TRUE, ordering = FALSE), rownames = FALSE),
              br())
          } else{
            message<-'Match Grade evaluations for recipient ITL %s. Please download the log and e-mail it to %s to troubleshoot.'
            lgr$info(sprintf(message, paste(patient, 'failed'), 'livtran@stanford.edu'))
            lgr$info('**********MATCH GRADE EVALUATION END**********')
            hide_spinner()
            HTML(sprintf(message, paste(patient, ' <b><span style=color:red;>failed</b></span>', sep=""), '<b>livtran@stanford.edu</b>'))
          }
        })
        
        #MG DB
        output$mgPopulationTable<-renderUI({
          req(mgp$result)
          if(mgp$result == 'TRUE'){
            con <- dbConn()
            
            mgPopulationTable <- dbGetQuery(con, sprintf("SELECT DSA, ABCDRDQ_alleles, ABCDRDQ_match,
                                                   ABCDRB1_alleles, ABCDRB1_match, ABCDRB1_mm_GVH, ABCDRB1_mm_HVG,
                                                   DRB345DQDP_alleles, DRB345DQDP_match, DRB345DQDP_mm_GVH, DRB345DQDP_mm_HVG,
                                                   DRB345DQDP_mm_TCE, Seven_loci_alleles, seven_loci_match
                                                   FROM dbo.Match_grades
                                                   WHERE recipient_number = %s and donor_number = %s", patient, donor))
            
            dbDisconnect(con)
            
            #remove any trailing spaces
            mgPopulationTable<-as.data.frame(lapply(mgPopulationTable, str_trim), stringsAsFactors = F)
            #make read in data numeric types numeric
            mgPopulationTable[,c(2:11, 13:length(mgPopulationTable))]<-lapply(mgPopulationTable[, c(2:11, 13:length(mgPopulationTable))], as.numeric)
            #replace any NAs with ''
            mgp$mgDBTable<-mgPopulationTable %>%
              replace(is.na(.), '') %>%
              relocate(ABCDRDQ_match, .before = ABCDRDQ_alleles) %>%
              relocate(ABCDRB1_match, .before = ABCDRB1_alleles) %>%
              relocate(DRB345DQDP_match, .before= DRB345DQDP_alleles) %>%
              relocate(Seven_loci_alleles, .before= seven_loci_match)
            
            tagList(
              hr(),
              h3('Match Grade Results', style = "text-align: center; font-weight: bold;"),
              DT::renderDT(mgp$mgDBTable, options = list(dom = 't', scrollX = TRUE, ordering = FALSE), rownames = FALSE))
          } else{
            NULL
          }
        })
      }
      
      #download log UI
      observe({
        req(mgp$result)
        output$log<-renderUI(
          tagList(
            hr(),
            div(style="text-align:center;",
                downloadButton(ns("downloadlog"), label = "Download log")
            )
          ))
      })

      #download log client
      output$downloadlog <- downloadHandler(
        filename <- function() {
          basename(final$log_path)},

        content <- function(file) {
          file.copy(final$log_path, file)
        }
      )
    }
  )
}

