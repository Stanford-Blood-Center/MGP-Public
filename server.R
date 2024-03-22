#server
#by: Livia Tran
#v1.0

suppressPackageStartupMessages(library(odbc))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(dplyr))
library(DBI)

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
  
  username<-reactiveValues(log=NULL)
  
  #validate mTilda username; capture mTilda username to use in log
  observeEvent(input$validate, {
    
    con <- dbConn()
    
    username_validation <- dbGetQuery(con, sprintf("select * from dbo.Staff where user_id = '%s'", input$username))

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
    output$spinner <- NULL
    run$result<-NULL
    run$log_path<-NULL
  })
  
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
    
    mg_res <- dbGetQuery(con, paste("select * from dbo.Match_grades where recipient_number = ", input$p_itl))
    
    if(nrow(mg_res)==0){
      output$p_check <- renderText({'No patients were found with the input ITL. Please make sure the input patient ITL is correct and that it was imported into Match Grade.'})
      disable('p_itl')
      dbDisconnect(con)
      return()
    } else{ 
      
       donors <- mg_res %>%
          filter(recipient_number != donor_number) %>%
          select(donor_number) %>%
          pull()
       
         output$donor_itl <- renderUI({
           tagList(
            #div(style = "margin-bottom:15px"),
            pickerInput('d_itl', 'Donor Selection', 
                        choices = c(donors), multiple = TRUE, 
                        options = pickerOptions(actionsBox = TRUE, noneSelectedText = 'Select Donors'), 
                        selected = NULL
                        )
         )})
    }
    
    dbDisconnect(con)
    
    })
  
  
  donor<-reactiveValues(selection=NULL)
  
  #save donor selected ITLs to a reactive value
  observeEvent(input$d_itl, {
    
    donor$selection<-input$d_itl
    
  })
  
  #logic for enabling Run MGP button
  observe({
    if(!is.null(input$d_itl)){
      enable('run')
    } else{
      disable('run')
    }
    
  })
    
  run<-reactiveValues(result=NULL, log_path=NULL)
  
    #run MGP button logic
    observeEvent(input$run, {
      
      output$spinner <- renderUI({
        withSpinner(uiOutput('finish_text'), color='#DC143C')
      })
      
      output$finish_text<-renderUI({

        run_res<-calcMatchGrade(input$p_itl, donor$selection, username$log)

        run$result<-run_res[1]
        run$log_path<-run_res[2]

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