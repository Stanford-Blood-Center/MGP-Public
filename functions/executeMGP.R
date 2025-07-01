#Match Grade Populator © Stanford Blood Center, LLC.
#v 1.12.0

library(shiny)

executeMGPUI<- function(id, type) {
  
  ns <- NS(id)
  
  if(type == 'spinner'){
    uiOutput(ns('spinner'))
  } else{
    uiOutput(ns('log'))
  }
}


executeMGPServer<-function(id, creds, patient, donor, hlaRecipient, hlaDonor, syn_nonsyn, donorMatchGrade, log, donor_filter, recip_filter) {
 
   moduleServer(
     id,
     function(input, output, session) {
       
       ns <- session$ns
       
       mgp<-reactiveValues(result=NULL, mgDBTable=NULL, mgpTable=NULL, missing_seq=NULL, A_mm = NULL, B_mm = NULL, C_mm = NULL, DRB1_mm = NULL, DRB345_mm= NULL, DQA1_mm = NULL,
                           DQB1_mm = NULL, DPA1_mm = NULL, DPB1_mm = NULL, fail_message = NULL, df = NULL)
       final<-reactiveValues(log_path = log, result=NULL, df=NULL, mgpTable=NULL, mgDBTable=NULL)
        
                                                ##### OUTPUT DIFFERENCES #####
              #if technologist level, only output finish text to indicate success or failure
              #else, output comparison text, data from the Match_grades table in DB, and what MGP calculated
       
       if(creds == '30'){
        output$spinner <- renderUI({
            tagList(uiOutput(ns('errorText')), uiOutput(ns('dsaMess')), uiOutput(ns('mmAlleles')))
          })
       } else if(creds %in% c('50', '60')){
        output$spinner <-renderUI({
          tagList(uiOutput(ns('errorText')), uiOutput(ns('comparisonText')), uiOutput(ns('dsaMess')), uiOutput(ns('mgTable')), uiOutput(ns('mgPopulationTable')), uiOutput(ns('mmAlleles')))
        })
       }
       
       observe({
         
         run_res<-isolate({calcMatchGrade(patient, donor, creds, hlaRecipient, hlaDonor, syn_nonsyn, donorMatchGrade, donor_filter, recip_filter)})
         
         mgp$result<-run_res[[1]]
         mgp$df<-run_res[[2]]
         mgp$missing_seq<-run_res[[3]]
         mgp$fail_message<-run_res[[4]]

         if(mgp$result){
           mgp$A_mm<-run_res[[5]]
           mgp$B_mm<-run_res[[6]]
           mgp$C_mm<-run_res[[7]]
           mgp$DRB1_mm<-run_res[[8]]
           mgp$DRB345_mm<-run_res[[9]]
           mgp$DQA1_mm<-run_res[[10]]
           mgp$DQB1_mm<-run_res[[11]]
           mgp$DPA1_mm<-run_res[[12]]
           mgp$DPB1_mm<-run_res[[13]]
           mgp$dsainfo<-run_res[[14]]

           message<-paste('Match Grade evaluations for recipient ITL', patient, 'completed!')
           if(!is.null(mgp$missing_seq)){
             showModal(missingModal(mgp$missing_seq))
           }
           lgr$info(message)
           lgr$info('**********MATCH GRADE EVALUATION END**********')
           hide_spinner()
         } else{
           output$errorText<-renderUI({
             
               message<-'Match Grade evaluations for recipient ITL %s %s'
               if(!is.null(mgp$fail_message)){
                 #no IgG tests selected error
                 message<-paste(message, mgp$fail_message, sep = '')
                 lgr$info(sprintf(message, patient, 'failed. '))
                 uiOutMess<-HTML(sprintf(message, patient, ' <b><span style=color:red;>failed</b></span>. '))
               } else{
                 message<-paste(message, ' Please download the log and e-mail it to %s to troubleshoot.', sep = '')
                 
                 email<-getMaintainerEmail()
                 
                 lgr$info(sprintf(message, patient, 'failed. ', email))
                 uiOutMess<-HTML(sprintf(message, patient, ' <b><span style=color:red;>failed</b></span>.', sprintf('<b>%s</b>', email)))
               }
               lgr$info('**********MATCH GRADE EVALUATION END**********')
               hide_spinner()
               uiOutMess
          })
         }
       })
      
                                                     ##### REVIEW MODE OUTPUTS #####
    if(creds %in% c('50', '60')){
      
      #difference comparison output
      output$comparisonText<-renderUI({
        req(mgp$result)
        if(mgp$result == 'TRUE'){
          
          comparison<-all.equal(mgp$mgDBTable, mgp$mgpTable)
          
          baseText<-"<b>Mismatched columns: </b>"
          
          if(length(comparison)==1 & any(!grepl('Component', comparison))){
            comparison<-'<h4 style="font-weight:bold; text-align:center;">None</h4>'
          } else{
            #[\u201C\u201D] = curly quotes
            comparison<-paste(gsub("[\u201C\u201D]", "", str_extract(comparison, "(?<=Component )[^:]+")), collapse = ', ')
          }
          tagList(
            h3('Mismatched columns', style = 'text-align: center; font-weight: bold;'),
            HTML(paste0('<div style="text-align: center;">', comparison, '</div>'))
          )
        } else{
          NULL
        }
      })
      
      #MGP output
      output$mgTable<-renderUI({
        
        if(mgp$result){
          
          mgp$mgpTable<-mgp$df %>%
            replace(is.na(.), '') %>%
            relocate(ABCDRDQ_match, .before = ABCDRDQ_alleles) %>%
            relocate(ABCDRB1_match, .before = ABCDRB1_alleles) %>%
            relocate(DRB345DQDP_match, .before= DRB345DQDP_alleles) %>%
            relocate(Seven_loci_alleles, .before= seven_loci_match)
          
          hide_spinner()
          
          if(!is.null(mgp$missing_seq)){
            showModal(missingModal(mgp$missing_seq))
          }
          
          tagList(
            hr(),
            h3('MGP Results', style = "text-align: center; font-weight: bold;"),
            DT::renderDT(mgp$mgpTable, options = list(dom = 't', scrollX = TRUE, ordering = FALSE), rownames = FALSE),
            br())
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
              DT::renderDT(mgp$mgDBTable, options = list(dom = 't', scrollX = TRUE, ordering = FALSE), rownames = FALSE), 
              hr())
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
        
        output$mmAlleles<-renderUI({
          req(mgp$result)

          if(mgp$result == 'TRUE'){
            
            mm_mess<-paste0(mgp$A_mm, mgp$B_mm, mgp$C_mm, mgp$DRB1_mm, mgp$DRB345_mm, mgp$DQA1_mm, mgp$DQB1_mm, mgp$DPA1_mm, mgp$DPB1_mm)
            
            #if the number of loci (indicated by style yellow due to highlighting) * 2 is equal to the 
            #number of Nones (None per GvH and HvG), return None for mismatches
            if(str_count(mm_mess, 'None') == str_count(mm_mess, 'yellow')*2){
              mm_mess<-'<h4 style="font-weight:bold; text-align:center;">None</h4>'
            }

            tagList(
              h3('Mismatches', style = 'font-weight:bold; text-align:center; margin-top: 1px;'),
              HTML(mm_mess))
          } else{
            NULL
          }
        })
        
        #DSA messaging 
        output$dsaMess<-renderUI({
          req(mgp$result)
          
          if(mgp$result == 'TRUE'){
            if(nrow(mgp$dsainfo)!=0){
              
              dsaTable<-DT::renderDT(mgp$dsainfo, options = list(dom = 't', scrollX = TRUE, ordering = FALSE,  columnDefs = list(
                list(targets = "_all", className = "dt-center")
              )), rownames = FALSE)
              
              #output at top; needs separation line after
              if(creds == '30'){
                tagOut<-tagList(
                  h3('DSA', style = 'font-weight:bold; text-align:center; margin-top: 1px;'),
                  dsaTable,
                  hr()
                )
              } else{
                tagOut<-tagList(
                  hr(),
                  h3('DSA', style = 'font-weight:bold; text-align:center; margin-top: 1px;'),
                  dsaTable,
                )
              }
              tagOut
            } else{
              NULL
            }
          }
        })
        
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

