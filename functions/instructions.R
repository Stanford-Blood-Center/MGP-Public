#Instructions section 


getInstructions<-function(){
  HTML(
      "
  <body>

    <h1><b>Instructions</b></h1>
    <h6><i><b>Version 1.8.0</b></i></h6>
    <br>

    <div class='instruction'>
        <b>Before beginning, please make sure donors for the recipient of interest have been moved from the 'Phenotypes'
            section of Match Grade to the middle window of the program. If the recipient of interest does not have IgG tests, 
            do <u>NOT</u> populate the 'Recipient DSA Date' field. Otherwise, please make sure the 'Recipient DSA Date' is populated 
            in the Match Grade application; this is required for DSA evaluation.</b>
        <br>
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
        the desired donor for evaluation.
        <br>
        <b><u>Please make sure any patients with manually entered typings from other facilities are excluded.
                The application will calculate values for those donors, but they may be incorrect.</b></u>
    </div>

    <div class='instruction'>
        <b>4)</b> Click on the 'Run MGP' button to begin evaluation. A red, busy spinner will appear in the main panel
        while the program is running. If there are no errors and the run was successful, a screen highlighting allele specific mismatches by locus will appear. 
        Results will be available for review in the Match Grade application. If the job failed, an error message will appear.
    </div>
    
    <div class='instruction'>
        <b>4a)</b> A pop-up window will appear at the end of the job if the recipient and/or donor has an allele that does not have a complete sequence for the ARD in the IMGT protein alignments. 
        Double check counts for the categories that the loci of the alleles appear in. Click 'OK' to exit.
    </div>
    
    <div class='instruction'>
        <b>4b)</b> If the recipient and/or donor has novel alleles, a window with some questions will appear for each novel allele. The window title will indicate
        which novel allele is being referenced. The first question is what position the mutation in the novel allele is at. Click 'Submit' after inputting
        the position. If the position is in the antigen recognition domain, a second question asking if the mutation is synonymous or non-synonymous will appear. A 'Next' button
        will be present if there are multiple novel alleles. A 'Finish' button will be present for the final novel allele. 
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
