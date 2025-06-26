#server
#by: Livia Tran
#v1.12.5

suppressPackageStartupMessages(library(odbc))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(shinybusy))

source('functions/functions.R')
source('functions/calcMatchGrade.R')
source('functions/instructions.R')
source('functions/refreshModal.R')
source('functions/executeMGP.R')

server <- function(input, output, session) {
  
  disable('validate')
  
  observe({
    req(!is.null(input$username))
    
    if(input$username != ''){
      enable('validate')
    } else{
      disable('validate')
    }
  })

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
    isolate(updateTextInput(session, "p_itl", value=""))
    output$d_check <- NULL
    output$p_check <- NULL
    output$q_check <- NULL
    output$donor_itl <- NULL
    output$p_itl_oops <- NULL
    output$typingException <- NULL
    main$error <- NULL
    novel$syn <- NULL
    novel$position<-NULL
    novel$click<-1
    novel$selection_syn<-NULL
    output$synonymousQuestion<-NULL
    output$mainPanel<-NULL
    output$sideLog<-NULL
    
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
  
  patient<-reactiveValues(itl=NULL)
  
  observeEvent(input$p_itl, {
    patient$itl<-isolate(input$p_itl)  
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
  
  novel<-reactiveValues(syn=NULL, position=NULL, click=1, selection_syn = list(), region_list = list(), totalAlleles = NULL, listNames = NULL, region = NULL)
  hla<-reactiveValues(recip=NULL, donor=NULL, recipNull=NULL, dMG=NULL, log=NULL, filter_donor = NULL, filter_recip = NULL)
  counter <- reactiveVal(1)
  main<-reactiveValues(error=NULL)
  
  #run MGP button logic
  observeEvent(input$run, {
    
    show_spinner()
    
    con<-dbConn()
    
    #set up log file
    hla$log<-tempfile(tmpdir=tempdir(),pattern = paste("MG_", patient$itl, '_', sep=""), fileext = ".log")
    lgr$add_appender(AppenderFile$new(hla$log), name = "mg_log")

    lgr$info(paste('Executed by mTilda user: ', username$log, sep = ''))
    
    spUsername<-Sys.getenv('SHINYPROXY_USERNAME')
    institution<-Sys.getenv('INSTITUTION_ID')
    
    #default to generic 'Institution ID' if not set as env variable
    if(institution == ""){
      institution<-"Institution ID"
    }
    
    #do not print institution id log line if there isn't a SHINYPROXY_USERNAME set 
    if(spUsername != ""){
      lgr$info(paste(institution, ': ', spUsername, sep = ''))
    }

    if(username$creds == '30'){
      lgr$info('***** CALCULATION MODE *****')
    } else if(username$creds %in% c('50', '60')){
      lgr$info('***** REVIEW MODE *****')
    }
    
    lgr$info(sprintf('Patient ITL: %s', patient$itl))
    
    tryCatch(
      {
    lgr$info('Extracting Match Grade data...')
    mg<-getMatchGrade(con, patient$itl)
    r_mg<-mg %>%
      filter(donor_number == patient$itl)
    
    lgr$info('Extracting recipient typing...')
    recip_Typing<-getTyping(con, r_mg, type = 'r')
    hla$recip<-recip_Typing[[1]]
    recip_novel_alleles<- recip_Typing[[3]]
    hla$filter_recip<-recip_Typing[[4]]
    
    hla$dMG<-mg %>%
      filter(donor_number == donor$selection)
    
    lgr$info('Extracting donor typing...')
    
    donor_Typing<-getTyping(con, hla$dMG, type = 'd')
  },
    error = function(e){
      hide_spinner()
      lgr$fatal(e)
      
      message<-'Match Grade evaluations for recipient ITL %s %s Please download the log and e-mail it to %s to troubleshoot.'
      
      email<-getMaintainerEmail()
      
      lgr$info(sprintf(message, patient$itl, 'failed. ', email))
      main$error<-HTML(sprintf(message, patient$itl, ' <b><span style=color:red;>failed</b></span>.', sprintf('<b>%s</b>', email)))

      }
  )
    
    if(!is.null(main$error)){
      
      output$typingException<-renderUI(main$error)
      
      output$sideLog<-renderUI(
        tagList(
          hr(),
          div(style="text-align:center;",
              downloadButton("downloadlog", label = "Download log")
          )
        ))
      
      #download log client
      output$downloadlog <- downloadHandler(
        filename <- function() {
          basename(hla$log)},
        
        content <- function(file) {
          file.copy(hla$log, file)
        }
      )
      
      return()
    }
    
    hla$donor<-donor_Typing[[1]]

    donor_novel_alleles<-donor_Typing[[3]]
    hla$filter_donor<-donor_Typing[[4]]
    
    novel$totalAlleles<-c(recip_novel_alleles, donor_novel_alleles)
    
    if(!is.null(novel$totalAlleles)){
      hide_spinner()
      novel$listNames<-gsub('Recipient: |Donor: ', '', novel$totalAlleles)
      refreshModal(novel$totalAlleles, novel$click)
      
    } else{
      
      module_id<-generateID(counter())
      counter(counter() + 1)
      
      output$mainPanel <- renderUI({
        executeMGPUI(module_id, 'spinner')
      })
      
      output$sideLog <- renderUI({
        executeMGPUI(module_id, 'log')
      })
      executeMGPServer(module_id, username$creds, patient$itl, donor$selection, hla$recip, hla$donor, novel$selection_syn, hla$dMG, hla$log, hla$filter_donor, hla$filter_recip)
    
      }
  })
  
  observeEvent(input$finish,{
    removeModal()

    show_spinner()
    module_id<-generateID(counter())
    counter(counter() + 1)
    
    output$mainPanel <- renderUI({
      executeMGPUI(module_id, 'spinner')
    })
    
    output$sideLog <- renderUI({
      executeMGPUI(module_id, 'log')
    })
    
    executeMGPServer(module_id, username$creds, input$p_itl, donor$selection, hla$recip, hla$donor, novel$selection_syn, hla$dMG, hla$log, hla$filter_donor, hla$filter_recip)

  })
  
  #next button
  observeEvent(input$nextButton, {
    novel$position <- NULL
    novel$click <- novel$click + 1
    removeModal()
    output$synonymousQuestion <- NULL
    refreshModal(novel$totalAlleles, novel$click)
  })

  observeEvent(input$syn_nonsyn,{
    novel$selection_syn[[novel$listNames[novel$click]]]<-input$syn_nonsyn
  })
  
  # #enable/disable button logic for 'Submit' button 
  observe({
    req(!is.null(novel$region))

    if(novel$region == ""){
      disable('submit')
    } else{
      enable('submit')
    }
  })
  
  #enable/disable button logic for 'Next' and 'Close' buttons
  observe({
    req(!is.null(novel$position))
    
    if(novel$position==FALSE){
      shinyjs::enable('nextButton')
      shinyjs::enable('finish')
    } else{
      if(!is.null(input$syn_nonsyn)){
        shinyjs::enable('nextButton')
        shinyjs::enable('finish')
      }
    }
  })
  
   observeEvent(input$region,{
     novel$region<-input$region
   })
  
  observeEvent(input$submit,{
    locus<-str_extract(novel$listNames[novel$click], '(?<=-)[^//*]+')
    novel$position<-determinePosition(novel$region, locus)

    shinyjs::hide('submit')
    
    output$synonymousQuestion<-renderUI({
      if(novel$position == TRUE){
        br()
        br()
        radioButtons('syn_nonsyn', 'Is the mutation in the novel allele synonymous?', c('Yes' = 'y', 'No' = 'n'), selected = character(0))
      }
    })
    
    #prevent flickering after clear button
    outputOptions(output, "synonymousQuestion", suspendWhenHidden = FALSE)
    
  })

  #instructions
  inst<-getInstructions()
  
  output$instructions <- renderUI({
    inst
  })
  
  
}