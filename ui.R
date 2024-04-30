#ui
#by: Livia Tran
#v1.1.0

suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinythemes))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(shinyWidgets))
suppressPackageStartupMessages(library(shinycssloaders))

ui <- navbarPage(title = div('Match Grade Populator (MGP) v 1.1.0', collapsible = FALSE, img(src = "sbc_logo.png", height = "45px",width = "100px",style = "position: fixed; right: 10px; top: 5px;")),theme = shinytheme('cosmo'),
                 
                 
                 #main tab 
                 tabPanel("Main",
                          
                          #enable shinyjs
                          useShinyjs(),
                          
                          #sidebar layout 
                          sidebarLayout(
                            sidebarPanel(
                              width =4,
                              textInput('p_itl', 'Patient ITL', ''),
                              column(
                                12,
                                actionButton('query_donor', 'Query Available Donors', style = "margin-bottom:15px; display:inline-block; text-align:center;"),
                                align='center'
                              ),
                              uiOutput('donor_itl'),
                              fixedRow(
                                column(
                                  12,
                                  actionButton('run', 'Run MGP'),
                                  actionButton('clear', 'Clear Inputs'),
                                  align='center')), 
                              uiOutput('log')
                            ),
                            
                            #main panel
                            mainPanel(
                              #style=paste('font-size: 20px; "flex-grow:1; resize:horizontal; overflow: hidden"'),
                              textOutput('p_itl'),
                              textOutput('p_check'),
                              textOutput('p_itl_oops'),
                              uiOutput('spinner')
                            )
                          )),
                 
                 
                 #about tab
                 tabPanel("Instructions", 
                          fluidPage(
                            tags$head(
                              tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
                            ),
                            fluidRow(
                              column(
                                width = 12,
                                div(
                                  class = "container",
                                  div(
                                    class = "instructions-box",
                                    uiOutput("instructions")
                                  )
                                )
                              )
                            )
                          )
                        )
                      )

