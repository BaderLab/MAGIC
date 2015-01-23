/**
**	This functions is responsible to create the div with its legend, html select tag and its options.
** 	The function also create an alert box that will be shown if the user does not select or check any option.
**/
function createDivSelectOptions(filterName, optionsData){


	//	Determine the legend of each div
	//	Class: First div - this is the x of the boxplot
	//	Filter: Others div - this is extra filters to filter the data ant plot the boxplot	
	if(filterName == "filter_0"){
		var titleName = "Class";
		var disabled = "disabled";
	}
	else{
		var titleName = "Filter " + filterName.split("_")[1];
		var disabled = "";

	}


	// Create the the div with its legend and a select campus
	$(".div-filters").append(
		'<div id="' + "div_" + filterName + '" class="select-filter span12 control-group">'+
			'<h1 class="div-title">'+
				'<p class="div-text-title">' + titleName + '</p>'+
			'</h1>'	+
			'<select class="form-control" id="'+ "select_" + filterName + '"'+disabled+'> </select>'
		);


	// Create the options for the select input.
	$(optionsData).each(function(){
			$("#select_"+filterName).append('<option value="' + this + '">' + this + '</option>');
	});

	// Start with nothing selected
	$("#select_"+filterName).prop('selectedIndex',-1);

	// Create a div which will be placed checbboxes related to the selected filter
	$("#div_" + filterName).append('<div id="' + "divcb_" + filterName + '"></div>');


	// Create an alert box for each div.
	$('#div_'+filterName).append('<div class="alert alert-danger fade in" id="alert_div_'+filterName+'" style="display: none;">'+
									'<p>Select a filter or an options before plot.</p>'+
					    		'</div>');

	// Create a checkbox which the user can select if he wants to plot one boxplot for each select option	
	if (filterName == "filter_1"){		
		$("#div_filter_1").append('<hr>');
		$("#div_filter_1").append('<input type="checkbox"  id="plotEachFilter" value="FALSE"> Plot for each options <br>');
	}
};


/**
**	This function is responsible for creating a form with a checkbox for each select options in accordance with the data 
**	from the lib_header file (data/lib_header.txt).
**/
function createFormCheckbox(checkboxData){

	for(var i=0; i<checkboxData.length; i++){
		
		var filterOptions = checkboxData[i].split("/");
		var checkboxValue = filterOptions[0].split(" ");

		if(filterOptions.length>1)
			var checkboxField = filterOptions[1].split("-");

		$form = $('<form id="' + "form_" + filters[i] + '" class="form-control" style="display: none;"></form>')
		$form.appendTo("body")

		if(filters[i] == "Age"){
			createAgeCheckbox(checkboxValue)
		}
		else{

			for(var j=0; j<checkboxValue.length; j++){
				if(checkboxValue[j] == "SHH?")
					checkboxValue[j] = "SHH2"
				$("#form_" + filters[i]).append('<input type="checkbox"  id="' + filters[i] + "_" +checkboxValue[j] + '" value="' + checkboxValue[j] + '"><span>'+checkboxField[j]+ '</span><br>');
			}
		}
	}
}


/**
**	This function is responsible to update the forms in real time in accord to the select options.
**/
function updateFormCheckbox(checkboxData){
	
	for(var i=0; i<checkboxData.length; i++){
		
		if (filters[i] != getAuxiliaryVar()){
			var filterOptions = checkboxData[i].split("/");
			var checkboxValue = filterOptions[0].split(" ");

			if(filterOptions.length>1){
				var checkboxField = filterOptions[1].split("-");
				
				for(var j=0; j<checkboxValue.length; j++){
					if(checkboxValue[j] == "SHH?")
						checkboxValue[j] = "SHH2"

					$('#'+filters[i]+"_"+checkboxValue[j] +"+span").html(checkboxField[j]);
				}
			}
			else{
				updateAgeCheckbox(checkboxValue)
			}
			
		}
	
	}
}


/**
**	This function is responsible for creating a form with its checkbox for each age pattern (american and european)
**	in accordance with the data from lib_header (data/lib_header.txt).
**/
function createAgeCheckbox(checkboxOptions){

	$("#form_Age").append('<input type="radio" id="americanRadio" name="age" value="American">American criteria<br>');
		$("#form_Age").append('<form id="americanAge" style="display: none;">');
			$("#americanAge").append('<input type="checkbox" id="american_Infant" value="Infant"><span>Infant - Age < 3 ('+ checkboxOptions[0] +')'+'</span><br>');
			$("#americanAge").append('<input type="checkbox" id="american_Child" value="Child"><span>Child - 3 <= Age < 16 ('+ checkboxOptions[1] +')'+'</span><br>');
			$("#americanAge").append('<input type="checkbox" id="american_Adult" value="Adult"><span>Adult - Age >= 16 ('+ checkboxOptions[2] +')'+'</span><br>');
			$("#americanAge").append('<input type="checkbox" id="american_NA" value="NA"><span>NA ('+ checkboxOptions[3] +')'+'</span><br>');

	$("#form_Age").append('<input type="radio" id="europeanRadio" name="age" value="European">European criteria<br>');
		$("#form_Age").append('<form id="europeanAge" style="display: none;">');
			$("#europeanAge").append('<input type="checkbox" id="european_Infant" value="Infant"><span>Infant - Age < 4 ('+ checkboxOptions[4] +')'+'</span><br>');
			$("#europeanAge").append('<input type="checkbox" id="european_Child" value="Child"><span>Child - 4 <= Age < 18 ('+ checkboxOptions[5] +')'+'</span><br>');
			$("#europeanAge").append('<input type="checkbox" id="european_Adult" value="Adult"><span>Adult - Age >= 18 ('+ checkboxOptions[6] +')'+'</span><br>');
			$("#europeanAge").append('<input type="checkbox" id="european_NA" value="NA"><span>NA ('+ checkboxOptions[7] +')'+'</span><br>');	
}

/**
**	This function is responsible for updating all of the options that correspond to the filter 'Age'.	
**/
function updateAgeCheckbox(checkboxOptions){
	$('#american_Infant +span').html("Infant - Age < 3 ("+checkboxOptions[0]+")");
	$('#american_Child+span').html("Child - 3 <= Age < 16 ("+checkboxOptions[1]+")");
	$('#american_Adult+span').html("Adult - Age >= 16 ("+checkboxOptions[2]+")");	
	$('#american_NA+span').html("NA ("+checkboxOptions[3]+")");

	$('#european_Infant+span').html("Infant - Age < 4 ("+checkboxOptions[4]+")");
	$('#european_Child+span').html("Child - 4 <= Age < 18 ("+checkboxOptions[5]+")");
	$('#european_Adult+span').html("Adult - Age >= 18 ("+checkboxOptions[6]+")");
	$('#european_NA+span').html("NA ("+checkboxOptions[7]+")");
}

/**
**	This function will add a new filter to the filters div
**/
function addNewFilter(){

	numberOfFilters++;

	if (numberOfFilters < filters.length){

		// Get all the select filters to create a list of options of the next select tag
		var selectedFilter = [];
		var noneSelection = 0;
		for(var i=0; i<numberOfFilters; i++){
			if($("#select_filter_"+i)[0].selectedIndex == -1){
				noneSelection++;
			}
			else{
					selectedFilter.push($("#select_filter_"+i).find(":selected").text());
				}
		}
		

		// Make a copy of the filters vector to creater a new one for the next select tag
		var auxiliaryFilterArray = filters.slice();

		// Remove the previous selected filters
		for(var i=0; i<selectedFilter.length; i++){	
			auxiliaryFilterArray = jQuery.grep(auxiliaryFilterArray, function(value) {
				return value != selectedFilter[i];
			});
		}

		// Delete the first element of the list if none was select in the previous filter
		if(noneSelection > 0){
			for (var i=0; i<noneSelection; i++){
				auxiliaryFilterArray.splice(0,1);
			}
		}
						
		// Create the select tag with its options
		createDivSelectOptions("filter_"+numberOfFilters, auxiliaryFilterArray);
	}

	if(numberOfFilters == filters.length - 1){
		$("#addFilterBtn").hide();
	}
}

/**
**	This function is responsible for checking if each input has a select filter and at least one 
**	selected option.
**	If none of the requirements are met, then the plot button is disabled.
**/
function checkSelectedFilters(){

		var uncheckedFilters = [];
		for(var i=0; i<=numberOfFilters; i++){
			
			var numberOfcheckedOptions = 0;
			var filter = $("#select_filter_"+i).val();
			$('#form_'+ filter +' input:checked').each(function() {
				numberOfcheckedOptions++;
			});

			
			comparisonValue = 0;
			if (filter == "Age"){
				comparisonValue = 1;
			}

			if(numberOfcheckedOptions == comparisonValue){
				uncheckedFilters.push("div_filter_"+i);
			}

		}

		if (uncheckedFilters.length > 0){
			return uncheckedFilters;
		}
		else{
			return true;
		}
}

function setAuxiliaryVar(value){
	auxiliaryVar = value;
}

function getAuxiliaryVar(){
	return auxiliaryVar;
}

/** 
**	Hides the warning message on the boxplot tab
**/
function hideWarningMessage() {
	$("#warningMessage").html("");
	$("#warningArea").hide();
}

// extending the basic javascript array object to allow for only unique elements to be obtained
// NOTE: this code relies on jquery's $.inArray() function
Array.prototype.unique = function () {
    var arr = this;
    return $.grep(arr, function (v, i) {
        return $.inArray(v, arr) === i;
    });
}