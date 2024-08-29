#pop up modal for novel alleles

refreshModal <- function(novelAlleles, clicks) {
  
  showModal(modalDialog(
    title = div(class='modal-title', novelAlleles[clicks]),
    textInput('region', 'What position is the mutation in?'),
    actionButton('submit', 'Submit'),
    isolate(uiOutput('synonymousQuestion')),
    footer = tagList(
      if(clicks == length(novelAlleles)){
        actionButton('finish', 'Close')
      } else{
        actionButton('nextButton', 'Next')
      }
    )
  ))
  disable('submit')
  disable('nextButton')
  disable('finish')
}