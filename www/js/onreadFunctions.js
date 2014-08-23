$(".filter-selection *").prop('disabled', false);
		$(document).ready(function() {	
			/* --------------------------------------------------------------------------------------------------------- 
														GENE INPUT SECTION
			----------------------------------------------------------------------------------------------------------*/

			/** 
			** Show textarea to insert patient id
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
			**	This function is responsible to validate what is inserted in the textarea gene (id #geneInput).
			**	The server verify is the gene or the list of gene exist in the file bioMart_table (Ensembl database v62),
			** 	if the list is valid the server return true and the textarea is colorized in green, else the server 
			**	return false, the textare is colorized in red and a alert box is shown to the user, showing what is wrong.
			**/
			Shiny.addCustomMessageHandler("verifyGeneInput", function(status){
				if (status == "BLANK"){
					$('.div-gene').removeClass('success error');	
					$("#geneAlert").hide();
					$('#plotBtn').hide();

					$(".filter-selection *").prop('disabled', true);
				}
				else{
					if (status == "TRUE"){
						$('.div-gene').removeClass('error');
						$('.div-gene').addClass('success');
						$("#geneAlert").hide();
						$('#plotBtn').show();

					
					var patientidClasses = $(".patient-id").attr('class')
					if(patientidClasses.indexOf("error") == -1)
						$(".filter-selection *").prop('disabled', false);
				
						// $(".filter-selection *").prop('disabled', false);
					}
					else{
						$('.div-gene').removeClass('success');
						$('.div-gene').addClass('error');

						$(".filter-selection *").prop('disabled', true);

						// Create an error message
						var alertMessage = "The gene(s): ";
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
						else{
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
			**	This function is responsible to validate what is inserted in the textarea patientid (id #patientidTextarea).
			**	The server verify is the patient id or the list of patient id exist in the file lib_header,
			**	if the list is valid the server return true and the textarea is colorized in green, else the server 
			**	return false, the textare is colorized in red and a alert box is shown to the user, showing what is wrong.
			**/	
			Shiny.addCustomMessageHandler("verifyPatientId", function(status){
				
				if (status == "BLANK"){
					$('.patient-id').removeClass('success error');	
					$("#patientAlert").hide();
					$('#plotBtn').hide();

					$(".filter-selection *").prop('disabled', true);
				}
				else{
					if (status == "TRUE"){
						$('.patient-id').removeClass('error');
						$('.patient-id').addClass('success');
						$("#patientAlert").hide();
						$('#plotBtn').show();

						var geneClasses = $(".div-gene").attr('class')
						if(geneClasses.indexOf("error") == -1)
							$(".filter-selection *").prop('disabled', false);
				

						
					}
					else{
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
						else{
							alertMessage += '"' + status + '"';
						}
						alertMessage += " do(es) not exist in our database."+
										"Please, insert one id per line. Make sure that there isn't space after each patient id";
						$('#alertPatient').html(alert);
						$("#patientAlert").show();
						$('#plotBtn').hide();
					}
					
				}
			})


			/* --------------------------------------------------------------------------------------------------------- 
														SELECT FILTER SECTION
			----------------------------------------------------------------------------------------------------------*/
			/**
			**	Create the field that correst to the class.
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
			**	Add new filter when Add new filter button is pressed
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
				else{
					var patientidClasses = $(".patient-id").attr('class')
					if(patientidClasses.indexOf("success") != -1)
						$(".filter-selection *").prop('disabled', false);
				}
			});


			/**
			**	This functions is responsible to get which options was select by the user, append it to the div and 
			**	show the checkbox form that was already created and had its data updated.
			**/
			$(".div-filters").on('change','select',function () { 
				$(this).next("div").html($('#form_'+ $(this).val()).clone());
				$('#form_'+ $(this).val()).show();


				$("#"+$(this).parent().attr('id')).removeClass('error');
				$("#alert_"+$(this).parent().attr('id')).hide();
			});

			
			/**
			**	Every time an options is select a json is sent to the server. The json contains all the selected filters
			**	and all the select options. The server is responsible to manage the json.
			**	In the server the json is split into 2 in the reactive function filter. The data is split into filters which 
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
			**	Receive a message from the server as son as the data is loaded in order to hide the loading modal
			**/
			Shiny.addCustomMessageHandler("hideModal", function(status){
			if(status == "TRUE")
				$('#myModal').modal('hide');
			});


			/**
			**	When Plot button is pressed, the tabs boxplot and data are shown and the tab Home is hidden.
			**	This function is also responsible to calculate the size of the div which the boxplot will be plotted.
			** 	The data only will be ploted if any options is selected or checked
			**/
			$('#plotBtn').on('click',function () {

				if(checkSelectedFilters() == true){
					// Hidde home tab, show boxplot and data tab and make boxplot tab as the active one
					$("#myTab li").show();
					$("#myTab .active").hide();
					$('#myTab a[href="#about"]').parent().hide();
					$('#myTab a[href="#boxplot"]').tab('show');
					$('#newPlot').show();

					// Get the number of gene to calculate the size of the div
					var numberOfGene = $("#geneInput").val().split("\n").length;
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

					// Set the width and the height of the div which the boxplot will be plotted
					$('.bp-style').css("width", "90%");
					$('.bp-style').css("height", heightBoxplotArea+"px");
				}
				else{
					for(var i=0; i<checkSelectedFilters().length; i++){
						$("#"+checkSelectedFilters()[i]).addClass('error');
						$("#alert_"+checkSelectedFilters()[i]).show();
					}
				}
			})


			/**
			** 	Show home and about tab in order to plot again
			**/
			$('#newPlot').on('click',function () {
				$('#newPlot').hide()

				$("#myTab li").hide();
				$('#myTab a[href="#about"]').parent().show();
				$('#myTab a[href="#home"]').tab('show');
				$("#myTab .active").show();					
			})
		})