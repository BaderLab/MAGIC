// Increased each time add filter button is pressed. The variable is used to manage the number of filter
// that is allowed to be created
var numberOfFilters = 0;

// Filters that can be choosed by the user.
var filters = [	"Age", "Gender", "Subgroup", "Tissue_Type"];

// The moment in which the user click in an option to filter the data (Filter Selection) this variable is set in accord
// to the input value.
// This value is used latter at updateFormCheckbox() in order to do not update the checbkboxes that the user checked.
var auxiliaryVar = null;