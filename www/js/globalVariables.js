// Increased each time add filter button is pressed. The variable is used to manage the number of filter
// that is allowed to be created
var numberOfFilters = 0;

// Filters that can be choosed by the user.
var filters = [	"Age", "Gender", "Subgroup", "Tissue_Type"];

// The moment in which the user click in an option to filter the data (Filter Selection) this variable is set in accord
// to the input value.
// This value is used latter at updateFormCheckbox() in order to do not update the checbkboxes that the user checked.
var auxiliaryVar = null;

// These variables will be used make sure the gene data table and at least one plot has been loaded properly
// prior to the loading modal being removed
var plotsLoaded = false;
var tablesLoaded = false;
var heatmapLoaded = false;

// When a user clicks the plotBtn after selecting no new values, the modal will lower, but the R code is never entered (because 
// it's based on reactive functions / are monitoring changes in the html values. Thus, this will monitor if no changes have
// been made
var pageAltered = false;
var geneAltered = false;