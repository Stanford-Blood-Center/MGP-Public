#Instructions section 


getInstructions<-function(){
  HTML(
      "
  <body>

    <h1><b>Instructions</b></h1>
    <br>

    <div class='instruction'>
        <b>Before beginning, please make sure donors for the recipient of interest have been moved from the 'Phenotypes'
            section of Match Grade to the middle window of the program.</b>
        <br>
        <u>Note:</u> 'Clear Inputs' can be clicked at any time to reset all fields back to their initial states.
    </div>

    <div class='instruction'>
        <b>1)</b> Enter the patient ITL to calculate Match Grade results for in the text box labeled 'Patient ITL'.
    </div>

    <div class='instruction'>
        <b>2)</b> Click on 'Query Available Donors'. This button will not be active until a patient ITL is provided. If
        the patient ITL entered is not found in the database or any letters are detected in the patient ITL entered,
        an error will appear in the main panel.
    </div>

    <div class='instruction'>
        <b>3)</b> A drop-down box populated with donors available for evaluation for the entered patient ITL. Select
        desired donors for Match Grade evaluation.
        <br>
        <b><u>Please make sure any patients with manually entered typings from other facilities are excluded.
                The application will calculate values for those donors, but they may be incorrect.</b></u>
    </div>

    <div class='instruction'>
        <b>4)</b> Click on the 'Run MGP' button to begin evaluation. Oscillating red lines in the main panel will appear
        while the program is running. If there are no errors and the run was successful, a message indicating completion
        will appear. Results will be available for review in the Match Grade application. If the job failed, an error message will appear.
    </div>

    <div class='instruction'>
        <b>5)</b> Regardless of the job status, a 'Download log' button will appear in the side panel. If the job failed,
        please download the log and e-mail it to <b>livtran@stanford.edu</b> for troubleshooting. If the job succeeded,
        downloading the log is optional. It is recommended to download and review the log to ensure everything ran smoothly.
        To view updated match values upon job success, click on 'Load Match Grade Donors'
        in the Match Grade application.
    </div>

    <h6><b>© 2024 Stanford Blood Center LLC. All rights reserved.</b></h6>

</body>"
    )
}
