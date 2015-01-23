$(".filter-selection *").prop('disabled', false);
		$(document).ready(function() {	
			/* --------------------------------------------------------------------------------------------------------- 
														GENE INPUT SECTION
			----------------------------------------------------------------------------------------------------------*/

			/** 
			** Show textarea to insert patient id when the checkbox is clicked
			**/
			$('#patientid').on('click',function () { 
				if($('#patientid').is(":checked")){
					$('#patientidTextarea').show();
					$("#plotBtn").hide()
				}
				else{
					$('#patientidTextarea').hide();
					$("#patientAlert").hide();

					$(".filter-selection *").prop('disabled', false);
					
					if($('.div-gene').hasClass("success"))
						$("#plotBtn").show();
					else
						$("#plotBtn").hide();
				}
				
			});


			/**
			**	This function is responsible for validating what is inserted in the textarea gene (id #geneInput).
			**	The server will verify that the gene or the list of genes exist in the bioMart table (Ensembl database v62).
			** 	If the gene(s) is/are valid, the server will return true, and the textarea will be coloured green, otherwise
			**	it will return false, and the text area is coloured red, and an alert box is displayed informing the user
			**	of what went wrong.
			**/
			Shiny.addCustomMessageHandler("verifyGeneInput", function(status){
				if (status == "BLANK") {
					// for if no genes were inputted
					$('.div-gene').removeClass('success error');	
					$("#geneAlert").hide();
					$('#plotBtn').hide();

					$(".filter-selection *").prop('disabled', true);
				}
				else {
					// for if genes were inputted
					if (status == "TRUE"){
						$('.div-gene').removeClass('error');
						$('.div-gene').addClass('success');
						$("#geneAlert").hide();
						$('#plotBtn').show();

					
						var patientidClasses = $(".patient-id").attr('class')
						if(patientidClasses.indexOf("error") == -1)
							$(".filter-selection *").prop('disabled', false);
					}
					else {
						$('.div-gene').removeClass('success');
						$('.div-gene').addClass('error');

						$(".filter-selection *").prop('disabled', true);

						// Create an error message
						var alertMessage = "The gene(s): ";
						if(status instanceof Array) {
							for(var i=0;i<status.length;i++){
								alertMessage += '"' + status[i] + '"';

								if(i+2 == status.length) {
									alertMessage += " and ";
								}
								else {
									if(i+1 < status.length)
										alertMessage += ", ";
								}
							}
						}
						else {
							alertMessage += '"' + status + '"';
						}
						alertMessage += " do(es) not exist at Ensembl database v62. "+
										"Please, insert one gene per line. Make sure that there isn't space after each gene<br><br>"+
										"Ensembl database v62: Ensembl Gene ID, Associated Gene Name, HGNC Symbol, "+ 
										"UCSC ID, RefSeq DNA ID, RefSeq Predicted DNA ID, RefSeq Protein ID, RefSeq "+
										"Predicted Protein ID and RefSeq Genomic IDs.";
						$('#alertGene').html(alertMessage);
						$("#geneAlert").show();
						$('#plotBtn').hide();
					}				
				}
			});


			/**
			**	This function is responsible for validating what is inserted in the textarea patientid (id #patientidTextarea).
			**	The server verifies if the patient id or the list of patient ids exist in the file lib_header.txt.
			**	If the list is valid, the server returns true and the textarea is coloured green, else the server 
			**	return false, so the textarea is coloured red and alert text is shown to the user, describing what is wrong.
			**/	
			Shiny.addCustomMessageHandler("verifyPatientId", function(status){
				
				if (status == "BLANK") {
					$('.patient-id').removeClass('success error');	
					$("#patientAlert").hide();
					$('#plotBtn').hide();

					$(".filter-selection *").prop('disabled', true);
				}
				else {
					if (status == "TRUE") {
						$('.patient-id').removeClass('error');
						$('.patient-id').addClass('success');
						$("#patientAlert").hide();
						$('#plotBtn').show();

						var geneClasses = $(".div-gene").attr('class')
						if(geneClasses.indexOf("error") == -1)
							$(".filter-selection *").prop('disabled', false);
					}
					else {
						$('.patient-id').removeClass('success');
						$('.patient-id').addClass('error');

						$(".filter-selection *").prop('disabled', true);

						// Create an error message
						var alertMessage = "The patient(s) id: ";
						if(status instanceof Array){
							for(var i=0;i<status.length;i++){
								alertMessage += '"' + status[i] + '"';

								if(i+2 == status.length)
									alertMessage += " and ";
								else
									if(i+1 < status.length)
										alertMessage += ", ";
							}
						}
						else {
							alertMessage += '"' + status + '"';
						}
						alertMessage += " do(es) not exist in our database."+
										"Please, insert one id per line. Make sure that there isn't space after each patient id";
						$('#alertPatient').html(alertMessage);
						$("#patientAlert").show();
						$('#plotBtn').hide();
					}
				}
			})


			/* --------------------------------------------------------------------------------------------------------- 
														SELECT FILTER SECTION
			----------------------------------------------------------------------------------------------------------*/
			/**
			**	Create the field that will be the 'class'.
			**/
			createDivSelectOptions("filter_0", filters);



			/**
			**	Show/hide the age when it is clicked.
			**/
			$('.div-filters').on('click','[name="age"]',function () { 
				if($('input[name=age]:checked').attr('value') == "American"){
					$("#americanAge").show();
					$("#europeanAge").hide();
					$('#europeanAge input:checked').each(function() {
					    $(this).prop('checked', false);
					});
				}
				else{
					$("#americanAge").hide();
					$("#europeanAge").show();
					$('#americanAge input:checked').each(function() {
					    $(this).prop('checked', false);
					});
					
				}
			});

			
			/**
			**	Add a new filter when the 'Add new filter' button is pressed
			**/
			$("#addFilterBtn").on("click", function(){
				addNewFilter();		
			});

			/**
			**	Reset all the select filters.
			**/
			$('#resetFilterBtn').on('click',function () {
				$('.div-filters').empty();
				numberOfFilters = 0;
				createDivSelectOptions("filter_0", filters);
				$("#addFilterBtn").show();	
				
				var geneClasses = $(".div-gene").attr('class')
				if(geneClasses.indexOf("success") != -1)
					$(".filter-selection *").prop('disabled', false);
				else {
					var patientidClasses = $(".patient-id").attr('class')
					if(patientidClasses.indexOf("success") != -1)
						$(".filter-selection *").prop('disabled', false);
				}
			});


			/**
			**	This functions is responsible for getting which options the user selected, appending it to the div and 
			**	showing the checkbox form that was already created and had its data updated.
			**/
			$(".div-filters").on('change','select',function () { 
				$(this).next("div").html($('#form_'+ $(this).val()).clone());
				$('#form_'+ $(this).val()).show();

				$("#"+$(this).parent().attr('id')).removeClass('error');
				$("#alert_"+$(this).parent().attr('id')).hide();
			});

			
			/**
			**	Every time options are selected, a json string is sent to the server containing all the selected filters
			**	and all the select options. The server is responsible for managing the json.
			**	In the server the json is split into 2 in the reactive function 'filter'. The data is split into filters which 
			**	contains all the selected filters and filter_x (x is a number) which  and contains the selected options.
			**	Example: if the filter Subgroups is at the position 1 one the vector filters them filter_1 contains all the
			** 	select options that correspond to the filter Subgroups
			**/
			$('.div-filters').on('click','input',function () {
				var json = []
				var agePattern = null
				for(var i=0; i<=numberOfFilters; i++){
					var selectedFilter = $("#div_filter_"+i).find(":selected").text();
					json.push(selectedFilter)

					 var options = "";
					 $('#form_'+ selectedFilter +' input:checked').each(function() {
					 	if($(this).attr('value') == "American" || $(this).attr('value') == "European")
					 		agePattern = $(this).attr('value')
					 	else{
					 		if( $(this).attr('value') == "SHH2")
					 			options = "SHH?" + " " + options;	
					 		else
						 		options = $(this).attr('value') + " " + options;
					 	}
					 });

					json.push(options)				
				}

			 	Shiny.onInputChange("json", json);

			 	if(agePattern != null)
			 		Shiny.onInputChange("agePattern", agePattern);

			 	var filter = $(this).closest('div').attr('id').split("_")[2]
			  	$("#div_filter_"+filter).removeClass('error');
				$("#alert_div_filter_"+filter).hide();

				setAuxiliaryVar($("#select_filter_"+filter).val());
			});


			/**
			**	Receive a json from the server. The json contains the data required to create the checkboxes
			**/
			Shiny.addCustomMessageHandler("getCheckboxData", function(data){
				createFormCheckbox(data);
			});


			/**
			**	Receive a json from the server. The json contains the data required to update the checkboxes in 
			**	accord to the select options
			**/
			Shiny.addCustomMessageHandler("udpateCheckboxOptions", function(data){
				updateFormCheckbox(data)
			});

			/* --------------------------------------------------------------------------------------------------------- 
														GENERAL FUNCTIONS
			----------------------------------------------------------------------------------------------------------*/
			/** 
			**	Handles confirming that a page has actually changed, such that the modal may be raised 
			**	if no changes were made (as none of the server.R code would be hit)
			**/
			$(".alterable").change(function() {
				pageAltered = true;
			});
			
			
			/** 
			**	Handles confirming whether or not a change has been made to the gene list 
			**	such that the heatmap page's error message will still display if only one gene is 
			**	selected multiple times - i.e. if the user attempts to plot the same gene but with
			**	different filters, they should still all have the error message displayed. 
			**/
			$("#geneInput").change(function() {
				geneAltered = true;
			});
			
			
			/**
			**	Receive a message from the server as soon as the data is loaded for the boxplot tab,
			**	in order to hide the loading modal, and make sure the page starts from the top
			**/
			Shiny.addCustomMessageHandler("hideBoxplotModal", function(status){
				if(status == "boxplot") {
					plotsLoaded = true;
				}
				else if(status == "tableData") {
					tablesLoaded = true;
				}
				
				if(plotsLoaded&& tablesLoaded) {
					$('html,body').scrollTop(0);
					$('#myModal').modal('hide');
					plotsLoaded = false;
					tablesLoaded = false;
				}
			});
			
			/**
			**	Receive a message from the server as soon as the data is loaded  for the heatmap tab, 
			**	in order to hide the loading modal and make sure the page starts from the top
			**/
			Shiny.addCustomMessageHandler("hideHeatmapModal", function(status){
				if(status == "heatmap") {
					heatmapLoaded = true;
				}
				else if(status == "tableData") {
					tablesLoaded = true;
				}
				
				if((heatmapLoaded&& tablesLoaded) || (status == "heatmapError")) {
					$('html,body').scrollTop(0);
					$('#myModal').modal('hide');
					heatmapLoaded = false;
					tablesLoaded = false;
				}
			});
			
			/** 
			*** This code is just used to display or hide the p value table checkbox according to 
			*** whether or not the p value can acutally be calculated for the dataset selected.
			**/
			Shiny.addCustomMessageHandler("displayPValueCheckbox", function(display) {
				if(display == "TRUE"){
					$('#pvalue-checkbox').show();
				}
				else {
					$('#pvalue-checkbox').hide();
				}
			});
			
			
			/**
			**	Defining the alert dialog box which will open when an error is caught by the server
			**/
			$( "#alertDialog" ).dialog({
				autoOpen: false,
				show: {
				effect: "blind",
				duration: 100
				},
				hide: {
				effect: "explode",
				duration: 100
				}
			});
			
			
			/** 
			*** Error handling code - currently the implementation is that an alert window will be displayed with the
			*** R error message, and the page will revert back to the original starting page (instead of being left in 
			*** an inbetween state/half completed state.
			**/
			Shiny.addCustomMessageHandler("caughtError", function(condition) {
				$('#myModal').modal('hide');
				
				$('#newPlot').hide();

				pageAltered = true;
				
				$("#myTab li").hide();
				$('#myTab a[href="#about"]').parent().show();
				$('#myTab a[href="#home"]').tab('show');
				$("#myTab .active").show();	
				
				$("#alertMessage").html(condition);
				$("#alertDialog").dialog("open");
				//window.alert(condition);
			});
			
			
			Shiny.addCustomMessageHandler("displayWarningMessage", function(message) {
				if(message.length > 0) 
				{
					$("#warningMessage").html(message);
					$("#warningArea").show();
				}
			});
			
			
			/**
			**	When Plot button is pressed, the tabs for the boxplot, heatmap and data are shown, and the tab Home is hidden.
			**	This function is also responsible for calculating the size of the div which the boxplot will be plotted on, 
			**	and generating the datetime value to be used in filenames.
			** 	The data will only be plotted if at least one option is selected or checked, otherwise, an error message
			**	will be displayed.
			**/
			$('#plotBtn').on('click',function () {				
				// Pass back the current date-time value to be used in the downloadable file names
				// (this should only be generated when the plot button has been clicked, such that ALL downloads from that 
				// session can be linked together by the date within the name)
				var currentDate = new Date();
				var formattedDate = currentDate.getFullYear() + '-' + (currentDate.getMonth()+1) + '-' + currentDate.getDate() + '_' + currentDate.toLocaleTimeString().replace(/:/gi, ".").replace(" ", "-") + '_';
				
				Shiny.onInputChange("currentDateTime", formattedDate);
				
				// make sure the warning message div is hidden to begin with
				// such that it really only displays if there is a proper warning
				hideWarningMessage();
				
				// Make sure the modal comes down once the plot button has been clicked
				$('#myModal').modal('show');
			
				if(checkSelectedFilters() == true){
					// Hide home tab, show boxplot, data and heatmap tabs and make boxplot tab the active one
					$("#myTab li").show();
					$("#myTab .active").hide();
					$('#myTab a[href="#about"]').parent().hide();
					$('#myTab a[href="#boxplot"]').tab('show');
					$('#newPlot').show();

					// Get the number of unique, non-empty, genes to calculate the size of the div
					var genes = $("#geneInput").val().toLowerCase().replace(/ /g, "\n").split("\n");
					genes = genes.filter(Boolean).unique();
					var numberOfGene = 0, i = genes.length;

					while (i--) {
						if (genes[i] != "")
							numberOfGene++;
					}
					
					// declare the variables for setting heights of the boxplot div on the page
					var heightBoxplotArea = 0;
					
					if($("#plotEachFilter").is(":checked")){
						// Send a message to the server warning that plot each filter is checked
						Shiny.onInputChange("plotEachFilter", true);

						// Get the number of selected options if plot for each filter is checked
					
						var numberOfSelectedOptions = 0;
					
						var filter = $("#select_filter_1").val();
						$('#form_'+filter+ ' input:checked').each(function() {
							numberOfSelectedOptions++;
						});

						// Subtract by one since the age pattern does not mattern at the sum
						if(filter == "Age")
							numberOfSelectedOptions--;

						// Calculate the size of the div which the boxplot will be plotted
						heightBoxplotArea = numberOfGene*numberOfSelectedOptions*500;
					}
					else{
						// Send a message to the server warning that plot each filter is not checked
						Shiny.onInputChange("plotEachFilter", false);

						// Since plot for each filter is not checked, the number of gene is the one that
						// matter at the sum.
						heightBoxplotArea = numberOfGene*500;
					}

					$('.bp-style').css("width", "90%");
					$('.bp-style').css("height", heightBoxplotArea+"px");
				}
				else{
					for(var i=0; i<checkSelectedFilters().length; i++){
						$("#"+checkSelectedFilters()[i]).addClass('error');
						$("#alert_"+checkSelectedFilters()[i]).show();
					}
					
					$('#myModal').modal('hide');
				}
				
				if(pageAltered == false) {
					$('#myModal').modal('hide');
				}
			})
			
			/**
			**	When the heatmap tab is selected, check to see if there are a valid number of genes
			**	being plotted (so that the correct sizes and divs can be selected for displaying.
			**	If there is only one gene, base the lowering of the modal on whether or not the gene had changed.
			**/
			$('#tab-heatmap').on('click', function(){
				// Get the number of unique, non-empty, genes to calculate the size of the div
				var genes = $("#geneInput").val().replace(/ /g, "\n").split("\n");
				genes = genes.filter(Boolean).unique();
				var numberOfGene = 0, i = genes.length;

				while (i--) {
					if (genes[i] != "")
						numberOfGene++;
				}
				
				if(numberOfGene > 1 ) {
					$(".toggle-on-error").show();
					$("#heatmapErrorMessage").hide();
					$("#heatmapWarningArea").hide();
					
					var heightHeatMapArea = numberOfGene*100 + 300;
					var legendHeight = 80 * Math.ceil((numberOfFilters+1)/3);
					
					// Set the width and the height of the div which the heatmap and legend will be plotted to
					$('.hm-style').css("width", "100%");
					$('.hm-style').css("height", heightHeatMapArea+"px");
					$('.hm-legend-style').css("width", "100%");
					$('.hm-legend-style').css("height", legendHeight+"px");
			
					if(pageAltered == true) {
						$('#myModal').modal('show');
					}
					else
					{
						$('#myModal').modal('hide');
					}
					
					// reset the check
					pageAltered = false;
				}
				else {
					$("#heatmapErrorMessage").show();
					$("#heatmapWarningArea").show();
					$(".toggle-on-error").hide();
					
					if(geneAltered == true) {
						$('#myModal').modal('show');
					}
					else
					{
						$('#myModal').modal('hide');
					}
					
					//reset the check
					geneAltered = false;
				}
			})


			/**
			** 	Reset everything and head back to the home page
			**/
			$('#newPlot').on('click',function () {
				$('#newPlot').hide();

				pageAltered = false;
				geneAltered = false;
				hideWarningMessage();
				
				$("#myTab li").hide();
				$('#myTab a[href="#about"]').parent().show();
				$('#myTab a[href="#home"]').tab('show');
				$("#myTab .active").show();					
			})
		})