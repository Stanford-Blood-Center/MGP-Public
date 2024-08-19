#server
#by: Livia Tran
#v1.3.1

suppressPackageStartupMessages(library(odbc))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(DT))

source('functions/functions.R')
source('functions/calcMatchGrade.R')
source('functions/instructions.R')

server <- function(input, output, session) {
  
  validate_modal <- modalDialog(
    textInput('username', 'mTilda username', ''),
    textOutput('username_not_found'),
    easyClose = F,
    footer = tagList(
      actionButton("validate", "Validate username")
    )
  )
  
  showModal(validate_modal)
  
  disable('query_donor')
  #require patient itl to enable donor query choice widget
  observe({
    if(input$p_itl != ""){
      enable('query_donor')
    } else{
      disable('query_donor')
    }
  })
  
  disable('run')
  
  username<-reactiveValues(log=NULL, creds=NULL)
  
  #validate mTilda username; capture mTilda username to use in log
  observeEvent(input$validate, {
    
    con <- dbConn()
    
    username_validation <- dbGetQuery(con, sprintf("select * from dbo.Staff where user_id = '%s'", input$username))
    username$creds <- username_validation$security_level
    
    if(nrow(username_validation)!=0){
      username$log<-input$username
      removeModal()
      
    } else{
      output$username_not_found<-renderText('The entered mTtilda username was not found. Please try again.')
    }
  })
  
  #clear button logic
  observeEvent(input$clear, {
    disable('run')
    enable('p_itl')
    updateTextInput(session, "p_itl", value="")
    output$d_check <- NULL
    output$p_check <- NULL
    output$q_check <- NULL
    output$donor_itl <- NULL
    output$p_itl_oops <- NULL
    output$log <- NULL 
    output$finish_text <- NULL
    output$mgTable<-NULL
    output$mgPopulationTable<-NULL
    output$comparisonText<-NULL
    run$result<-NULL
    run$log_path<-NULL
    run$df<-NULL
    run$mgpTable<-NULL
    run$mgDBTable<-NULL
    output$spinner <- NULL
  })
  
  donor<-reactiveValues(selection=NULL, donorDF=NULL)
  
  #donor selection
  observeEvent(input$query_donor,{
    
    disable('query_donor')
    disable('p_itl')
    
    if(grepl('[A-Za-z]', input$p_itl)){
      output$p_itl_oops<-renderText('The patient ITL entered has letters. ITLs should only be composed of numbers. Please re-enter a corrected ITL.')
      disable('p_itl')
      return()
    }
    
    con <- dbConn()
    
    mg_res <- dbGetQuery(con, paste("SELECT mg.*, p.last_name, p.first_name FROM dbo.Match_grades mg
                                     LEFT JOIN dbo.Patients p
                                     ON mg.donor_number = p.patient_number
                                     WHERE mg.recipient_number = ", input$p_itl))
    
    if(nrow(mg_res)==0){
      output$p_check <- renderText({'No patients were found with the input ITL. Please make sure the input patient ITL is correct and that it was imported into Match Grade.'})
      disable('p_itl') 
      dbDisconnect(con)
      return()
    } else{ 
      
      donor$donorDF <- mg_res %>%
        filter(recipient_number != donor_number) %>%
        mutate(name=paste(last_name, first_name, sep = ',')) %>%
        mutate(full=paste(name, donor_number, sep = ' - ')) %>%
        select(donor_number, name, full)
      
      donorEntry<- donor$donorDF %>%
        pull(full)
      
      output$donor_itl <- renderUI({
        tagList(
          selectInput('d_itl', 'Donor Selection', choices=c("Nothing selected", donorEntry), multiple = FALSE)
        )})
    }
    
    dbDisconnect(con)
    
  })
  
  #save donor selected ITLs to a reactive value
  observeEvent(input$d_itl, {
    
    #logic for enabling Run MGP button
    if(input$d_itl != 'Nothing selected'){
      enable('run')
    } else{
      disable('run')
    }
    
    donor$selection<-donor$donorDF %>%
      filter(full == input$d_itl) %>%
      select(donor_number) %>%
      pull()
    
    donor$selection<-as.character(donor$selection)
    
  })
  
  run<-reactiveValues(result=NULL, log_path=NULL, df=NULL, mgpTable=NULL, mgDBTable=NULL)
  
  #run MGP button logic
  observeEvent(input$run, {
    
    if(username$creds == '30'){
      output$spinner <- renderUI({
      withSpinner(uiOutput('finish_text'), color='#DC143C')
    })
    
    } else if(username$creds %in% c('50', '60')){
      output$spinner <-renderUI({
        withSpinner(tagList(uiOutput('comparisonText'), uiOutput('mgTable'), uiOutput('mgPopulationTable')), color='#DC143C')
      })
    }
    
    if(username$creds == '30'){
      
      output$finish_text<-renderUI({
        run_res<-calcMatchGrade(input$p_itl, donor$selection, username$log, username$creds)
        
        run$result<-run_res[[1]]
        run$log_path<-run_res[[2]]
        
        if(run$result){
          message<-paste('Match Grade evaluations for recipient ITL', input$p_itl, 'completed!')
          lgr$info(message)
          lgr$info('**********MATCH GRADE EVALUATION END**********')
          message
        } else{
          message<-'Match Grade evaluations for recipient ITL %s. Please download the log and e-mail it to %s to troubleshoot.'
          lgr$info(sprintf(message, paste(input$p_itl, 'failed'), 'livtran@stanford.edu'))
          lgr$info('**********MATCH GRADE EVALUATION END**********')
          HTML(sprintf(message, paste(input$p_itl, ' <b><span style=color:red;>failed</b></span>', sep=""), '<b>livtran@stanford.edu</b>'))
        }
      })
    } else if(username$creds %in% c('50', '60')){
      
        #MGP output
        output$mgTable<-renderUI({
          
          run_res<-calcMatchGrade(input$p_itl, donor$selection, username$log, username$creds)
          
          run$result<-run_res[[1]]
          run$log_path<-run_res[[2]]
          run$df<-run_res[[3]]
          
          if(run$result){
            run$mgpTable<-run$df %>% replace(is.na(.), '')
            
            tagList(
              hr(),
              h3('MGP Results', style = "text-align: center; font-weight: bold;"),
              DT::renderDT(run$mgpTable, options = list(dom = 't', scrollX = TRUE), rownames = FALSE), 
              br())
          } else{
            message<-'Match Grade evaluations for recipient ITL %s. Please download the log and e-mail it to %s to troubleshoot.'
            lgr$info(sprintf(message, paste(input$p_itl, 'failed'), 'livtran@stanford.edu'))
            lgr$info('**********MATCH GRADE EVALUATION END**********')
            HTML(sprintf(message, paste(input$p_itl, ' <b><span style=color:red;>failed</b></span>', sep=""), '<b>livtran@stanford.edu</b>'))
          }
        })
        
          #MG DB
          output$mgPopulationTable<-renderUI({
            req(run$result)
            if(run$result == 'TRUE'){
              con <- dbConn()
              
              mgPopulationTable <- dbGetQuery(con, sprintf("SELECT DSA, ABCDRDQ_alleles, ABCDRDQ_match, 
                                                   ABCDRB1_alleles, ABCDRB1_match, ABCDRB1_mm_GVH, ABCDRB1_mm_HVG, 
                                                   DRB345DQDP_alleles, DRB345DQDP_match, DRB345DQDP_mm_GVH, DRB345DQDP_mm_HVG,
                                                   DRB345DQDP_mm_TCE, Seven_loci_alleles, seven_loci_match
                                                   FROM dbo.Match_grades 
                                                   WHERE recipient_number = %s and donor_number = %s", input$p_itl, donor$selection))
              
              dbDisconnect(con)
              
              #remove any trailing spaces
              mgPopulationTable<-as.data.frame(lapply(mgPopulationTable, str_trim), stringsAsFactors = F)
              #make read in data numeric types numeric
              mgPopulationTable[,c(2:11, 13:length(mgPopulationTable))]<-lapply(mgPopulationTable[, c(2:11, 13:length(mgPopulationTable))], as.numeric)
              #replace any NAs with ''
              run$mgDBTable<-mgPopulationTable %>% replace(is.na(.), '')
              
              tagList(
                hr(),
                h3('Match Grade Results', style = "text-align: center; font-weight: bold;"),
                DT::renderDT(run$mgDBTable, options = list(dom = 't', scrollX = TRUE), width = '25%', rownames = FALSE))
            } else{
              NULL
            }
          })
          
          #difference comparison output
          output$comparisonText<-renderUI({
            req(run$result)
            if(run$result == 'TRUE'){
              
              comparison<-all.equal(run$mgDBTable, run$mgpTable)
              
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
        } 
    
  })
  
  
  #download log UI
  observe({
    req(run$result)
    output$log<-renderUI(
      tagList(
        hr(),
        div(style="text-align:center;", 
            downloadButton("downloadlog", label = "Download log")
        )
      ))
  })
  
  #download log client 
  output$downloadlog <- downloadHandler(
    filename <- function() {
      basename(run$log_path)},
    
    content <- function(file) {
      file.copy(run$log_path, file)
    }
  )
  
  #instructions
  inst<-getInstructions()
  
  output$instructions <- renderUI({
    inst
  })
  
  
}