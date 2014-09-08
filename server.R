#required libraries
library(shiny)


# Define server logic
shinyServer(function(input, output, session) {

	#----------------------------------------------------------------------------------------------------------#
	# 												READ DATA												   #
	#----------------------------------------------------------------------------------------------------------#


	#--------------------------------------- READ DATA FROM FOLDER DATA ---------------------------------------#
												

	# Data loaded to get the number of sample for each options of each filter once it is selected.
	lib_header <- reactive({
		lib_header = read.table("data/lib_header.txt", sep="\t", header=TRUE)

		return (lib_header)		
	})


	# Data that contains the genes and the rnaseq value
	rnaseq_rpkm <- reactive({
		rnaseq_rpkm = read.table("data/rnaseq_rpkm_matrix.txt", sep="\t", header=TRUE)

		return (rnaseq_rpkm)
	})	



	# Load Ensembl data v62. Used to validated the inputted gene and to get the Ensembl gene ID to plot the graphs.
	bioMart_table = read.table("data/bioMart_table.txt", sep="\t", header=TRUE)



	#---------------------------------------- GET DATA FROM CLIENT SIDE ----------------------------------------#

	# Get the inserted gene(s). If there is more than one, split them into an array.
	gene <- reactive({

		genes <- strsplit(input$geneInput, "\n")
		ensemblIDArray <- NULL
		for(i in 1:length(genes[[1]])){
			ensemblID <- which (bioMart_table == tolower(genes[[1]][i]), arr.ind=TRUE)
			ensemblIDArray <- c(ensemblIDArray, as.character(bioMart_table[ensemblID[1],1]))
		}

		return (ensemblIDArray)
	})


	# Get the status of the field Patiend ID. If it is selected then $patientid become true else false.
	# This reactive "function" also get the ids and split it into an array.
	patient <- reactive({
		patient <- list(patientid = input$patientid)
		patient$patientidtext <- strsplit(input$patientidTextarea, "\n")

		return (patient)
	})


	# Get the json sent from the client. The json contains all the selected filter and all the selected
	# options. The json is splitted in two arrays: $filters which contains all selected filters (subgroup as example)
	# and $filter_x in which x is a number. Each $filter_ contains the selected options.
	# If $filter[1] contains the filter subgroup, so $filter_1 contains the selected options.
	# This function also get the status of the checkbox Plot Each Filter that can be true (checked) 
	# or false (not checked).
	filter <- reactive({

		if(is.null(input$plotEachFilter))
			fltr <- list(plotEachFilter = FALSE)
		else
			fltr <- list(plotEachFilter = input$plotEachFilter)

		if(!is.null(input$json)){
			filters <- NULL
			for(i in seq(1, length(input$json), 2))
				filters <- c(filters, input$json[i])

			fltr <- c( fltr, list(filters = filters))

			j <- 1;
			for(i in seq(2, length(input$json), 2)){
				options <- strsplit(input$json[i], " ")[[1]]
				fltr[[paste("filter_",j,sep="")]] <- options
				j <- j + 1
			}

			for(i in 1:length(fltr$filters)){
				filter <- paste("filter_",i,sep="")
				fltr[[filter]] <- orderfilter(fltr[[filter]], fltr$filters[[i]])
			}
		}

		if(!is.null(input$agePattern))
			fltr <- c(fltr, list(agePattern = input$agePattern))


		return (fltr)
	})
		


	#----------------------------------------------------------------------------------------------------------#
	#											OBSERVER FUNCTIONS 											   # 
	#			Receive data from the client and automaticaly execute an action in accord to each data 		   #
	#----------------------------------------------------------------------------------------------------------#

	

	# Send to the client a json which contains all the data needed to created the checkbox options for each filter.
	# The sent json is an array which each position correspond to a filter. According to the file globalVariables.js 
	# the position 1 correspond to the Age and so on.
	# Each position of the array contains all the options and the number of sample in accord what was select by the 
	# user.
	# As the page is loaded the data is get from the hole sample (lib_header()) to create the checkboxes. When any 
	# option is selected the data is get from the reactive functions r_lib_header() and sent to the function that will 
	# update the data instead of create it as when the page is loaded.
	observe({

			filters <- c("Age", "Gender", "Subgroup", "Tissue_Type");

			ageId <- NULL
			if(length(filter()$filters) > 0){
				lh <- r_lib_header()
				ageId <- which(filter()$filters == "Age")
			}
			else
				lh <- lib_header()
			
			finalData <- NULL
			for (i in 1:length(filters)){
				filter <- filters[i]

				# Age has a different pattern then the other filters due to the data
				if(filter == "Age"){
					if(length(filter()[[paste("filter_",ageId,sep="")]]) > 0)
						checkboxOptions <- getNumberOfSampleAge2(lh)
					else
						checkboxOptions <- getNumberOfSampleAge(lh)

					checkboxData <- paste(checkboxOptions, collapse=" ")
					finalData <- c(finalData, checkboxData)
				}
				else{
					# Get the values
					checkboxValue <- levels(lib_header()[[filter]])
					checkboxValue <- checkNA(checkboxValue, filter, lh)
					checkboxValue <- orderfilter(checkboxValue,filter)

					# Get the options for the checkbox
					checkboxOptions <- getNumberOfSample(checkboxValue, filter, lh)

					# Format the json that will be sent to the client
					checkboxValue <- paste(checkboxValue, collapse=" ")
					checkboxOptions <- paste(checkboxOptions, collapse="-")
					checkboxData <- paste(checkboxValue,checkboxOptions,sep="/")
					
					# Receive all the value and options for each filter
					finalData <- c(finalData, checkboxData)
				}
			}
			
			# Send data to the client side
			if(length(filter()$filters) > 0)
				session$sendCustomMessage(type = "udpateCheckboxOptions", finalData)
			else
				session$sendCustomMessage(type = "getCheckboxData", finalData)

				gene_input$resume()
				patientid_input$resume()
	},priority = 2)

	# Verify if the gene inputted by the user exist at Ensembl database v62 (file at Data/bioMart_table.txt).
	# If exist return true, else return the gene(s) that do(es) not exist in our database;
	gene_input <- observe({

		if(input$geneInput == ""){
			session$sendCustomMessage(type = "verifyGeneInput", "BLANK")
		}
		else{
			genes <- strsplit(input$geneInput, "\n")
			matchGene <- TRUE
			wrongGenes <- NULL

			for(i in 1:length(genes[[1]])){
				gene <- tolower(genes[[1]][i])
				if (length(which (bioMart_table == gene)) == 0){
					matchGene <- FALSE
					wrongGenes <- c(wrongGenes, genes[[1]][i])
				}
			}

			if(matchGene == TRUE)
				session$sendCustomMessage(type = "verifyGeneInput", "TRUE")
			else
				session$sendCustomMessage(type = "verifyGeneInput", wrongGenes)
		}
	},suspend = TRUE)

	# Verify if the filed Patient ID is checked. If it is, then verify if the inputted patient(s) 
	# id(s) exist in our database.
	patientid_input <-observe({
		
		if (input$patientid == TRUE){
			if(input$patientidTextarea == ""){
				session$sendCustomMessage(type = "verifyPatientId", "BLANK")
			}
			else{
				patientids <- strsplit(input$patientidTextarea, "\n")

				matchPatientId <- TRUE
				wrongPatientIds <- NULL
				for(i in 1:length(patientids[[1]])){
					patientid <- toupper(patientids[[1]][i])
					if (length(which (lib_header()$Patient_ID == patientid)) == 0){
						matchPatientId <- FALSE
						wrongPatientIds <- c(wrongPatientIds, patientids[[1]][i])
					}
				}				

				if(matchPatientId == TRUE)
					session$sendCustomMessage(type = "verifyPatientId", "TRUE")
				else
					session$sendCustomMessage(type = "verifyPatientId", wrongPatientIds)
			}
		}
	},suspend = TRUE)


	#----------------------------------------------------------------------------------------------------------#
	#							MANIPULATE THE DATA TO SEND TO THE CLIENT SIDE								   #	
	#----------------------------------------------------------------------------------------------------------#

	# This function order the filter in accord to the stipulated order at the order list.
	# filter is a vector containing the disordered elements
	# order is a vector containing the desired order
	orderfilter <- function(filter, filterName){
		# These are the order of the filter in which order that will be ploted.
		if(filterName == "Subgroup")
			order <- c("WNT", "SHH", "Group3", "Group4", "control", "Normal")
		else
			if(filterName == "Tissue_Type")
				order <- c("control", "metastasis", "null", "primary", "recurrence")
			else
				order <- c("F", "M")

		
		length <- length(filter)
		orderedfilter <- NULL

		j <- 1
		while(j <= length(order)){
			if(length(which(filter==order[j])) > 0){
				filter <- manageExistingElement(order[j], filter)
				orderedfilter <- c(orderedfilter, order[j])
				order <- order[-j]
			}
			else
				j <- j + 1
		}

		orderedfilter <- c(orderedfilter, filter)

		return (orderedfilter)

	}

	# manageExistingElement search in the vector for a element. If the element
	# is found, it is removed from its original list.
	manageExistingElement <- function(element, filter){
		matchId <- which(filter==element)
		filter <- filter[-matchId]
	
		return (filter)
	}


	# toUppercase capitalize all the first letter of any string.
	toUppercase <- function(data){
		for(i in 1:length(data)){
	        s <- strsplit(data[i], " ")[[1]]
	        data[i] <- paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
		}

		return (data)
	}

	
	# Get the number of sample for the age since it has a different pattern compared to the others filters.
	getNumberOfSampleAge <- function(lib_header){

		# American
		american <- length(which(lib_header$Age < 3))
		american <- c(american, length(which(lib_header$Age >= 3 & lib_header$Age < 16)))
		american <- c(american, length(which(lib_header$Age >= 16)))
		american <- c(american, length(which(is.na(lib_header$Age))))

		# European
		european <- length(which(lib_header$Age < 4))
		european <- c(european, length(which(lib_header$Age >= 4 & lib_header$Age < 18)))
		european <- c(european, length(which(lib_header$Age >= 18)))
		european <- c(european, length(which(is.na(lib_header$Age))))
		
		return(c(american, european))
	}


	# Get the number of sample for the age when any field is selected by the user.
	getNumberOfSampleAge2 <- function(lib_header){

		ageData <- length(which(lib_header$Age == "Infant"))
		ageData <- c(ageData, length(which(lib_header$Age == "Child")))
		ageData <- c(ageData, length(which(lib_header$Age == "Adult")))
		ageData <- c(ageData, length(which(is.na(lib_header$Age))))


		if (filter()$agePattern == "American")
			ageData <- c(ageData, 0,0,0,ageData[4])
		else
			ageData <- c(0,0,0,ageData[4],ageData)

	}


	# Check and add NA to the list in case it is found.
	checkNA <- function(filterData, column, lib_header){
		# if (length(which(is.na(lib_header[[column]]))) > 0)
			filterData <- c(filterData, "NA")

		return (filterData)
	}

	# Function reponsible to get the number of sample for each element in each sample and to create
	# a string which contain the options and the name of sample. It is used to create the options
	# for each filter.
	# return example: Female (20)
	getNumberOfSample <- function(filterData, filter, lib_header){
		for(i in 1:length(filterData)){
			if (filterData[i] == "NA")
				noSample <- length(which(is.na(lib_header[[filter]])))
			else
				noSample <- length(which(lib_header[[filter]] == filterData[i]))


			if(filterData[i] == "F")
				filterData[i] <- "Female"
			else
				if(filterData[i] == "M")
					filterData[i] <- "Male"
				else
					filterData[i] <- toUppercase(filterData[i])
						
			filterData[i] <- paste(filterData[i], " ", "(", noSample, ")", sep="")
		}

		return(filterData)
	}

	
	#---------------------------------------------------------------------------------------------------------#
	#                                  MANAGE THE DATA TO PLOT THE BOXPLOT                                    #
	#---------------------------------------------------------------------------------------------------------#
	
	# This functions arrange the data in accord to the selected age type. It transforme each numeric field
	# to it specific name such as infant, child or adult.
	arrangeData <- function(lib_header, agePattern){
		
  		value <- NULL
		if(agePattern == "American"){
  			value$Infant <- 3
  			value$Adult <- 16
  		}
  		else{
  			value$Infant <- 4
  			value$Adult <- 18
  		}

  		matchedFilterIds <- NULL
  		
  		matchedFilterIds$Infant <- which(lib_header$Age < value$Infant)
		matchedFilterIds$Child <- which(lib_header$Age >= value$Infant & lib_header$Age < value$Adult)
		matchedFilterIds$Adult <- which(lib_header$Age >= value$Adult)
		
					
  		for(j in 1:3)	{
  			val <- names(matchedFilterIds)[j]
  			lib_header$Age[matchedFilterIds[[j]]] <- val
		}

		return (lib_header)
	}


	# Reduce the rnaseq_rpkm data. The reduced data frame contains only the genes inserted by the user.
	# It is made to reduce the time to process the data.
	r_rnaseq_rpkm <- reactive({
		matchedGeneIds <- NULL

  		# Get the number of the lines that correspond to each inserted gene
  		for (i in 1:length(gene())){
  			matchedGeneIds <- c(matchedGeneIds, which(rnaseq_rpkm() == toupper(gene()[i])))
  		}

  		# Create a data.frame that contains only the inserted gene
  		reduced_rnaseq_rpkm <- rnaseq_rpkm()[sort(matchedGeneIds),]

  		return (reduced_rnaseq_rpkm)
	})


	# Reduce the lib_header data. The reduced data frame contains only the selected filter and options.
	# It is made to reduce the time to process the data.
	r_lib_header <- reactive({

		# Create a data.frame that contains only the selected filters and the essential ones to create the table
  		# filters <- c("LIB", "Tissue_ID", "Patient_ID", filter()$filters)
  		# reduced_lib_header <- lib_header()[,filters]
  		reduced_lib_header <- lib_header()

  		# Reduce the reduced_lib_header in accord to each selected options for each filter
  		for(i in 1:length(filter()$filters)){
  			
  			# This value will be used to store the lines that match with the selected options
  			matchedFilterIds <- NULL

  			# Form filter name
  			filter <- paste("filter_",i,sep="")
	  			
	  			if (length(filter()[[filter]]) > 0 ){

	  			# If some selected filter is age, then the column Age from reduced_lib_header is transformed to 
	  			# the accord to the age pattern (american or european), i.e, trade each age to their correspondend
	  			# one (Infant, Child or Adult)
	  			if(filter()$filters[i] == "Age" ){
	  				reduced_lib_header <- arrangeData(reduced_lib_header, filter()$agePattern)
	  			}


	  			# Filter the reducied_lib_header in accord to the selected filters and options
	  			for(j in 1:length(filter()[[filter]])){
	  				if(filter()[[filter]][j] == "NA"){
	  					matchedFilterIds <- c(matchedFilterIds, which(is.na(reduced_lib_header[[filter()$filters[i]]])))
	  				}
	  				else{
	  					matchedFilterIds <- c(matchedFilterIds, which(reduced_lib_header[[filter()$filters[i]]] == filter()[[filter]][j]))
	  				}
	  			}

	  			reduced_lib_header <- reduced_lib_header[sort(matchedFilterIds),]
  			}
  		}

  		# If patient id checkbox is checked then the data will be filtred in accord to the patient id as well.
  		if(patient()$patientid == TRUE){

  			matchedFilterIds <- NULL
  			for(i in 1:length(patient()$patientidtext)){
  				matchedFilterIds <- c(matchedFilterIds, which(reduced_lib_header$Patient_ID == patient()$patientidtext[i]))
  			}
  			reduced_lib_header <- reduced_lib_header[sort(matchedFilterIds),]
  		}
  		
  		return (reduced_lib_header)
	})

	# Manage the data, i.e, create a data frame that contains the rnaseq value for the specific inserted
	# gene(s), patient(s) id(s) and selected filter and options.
  	boxplotData <- reactive({
  		print("2")
  		# Get the data only with the inserted genes
  		reduced_rnaseq_rpkm <- r_rnaseq_rpkm()

  		# Send message to the client warning the data is already read.
  		# This message hide the loading modal
  		session$sendCustomMessage(type = "hideModal", "TRUE")

  		# Get the data already filtered
  		reduced_lib_header <- r_lib_header()

  		# Match tables obtaining a vector containg all coluns at reduced_rnaseq_rpkm in which all userIds (LIB) match.
  		matchedSampleId <- match(reduced_lib_header$LIB, colnames(reduced_rnaseq_rpkm))

  		# Define the gene(s) as name of the rows
  		rownames(reduced_rnaseq_rpkm) <- reduced_rnaseq_rpkm[,1]

  		# Create a new data.frame containing only the data selected by the user (gene, patient id and selected filters)
  		reduced_rnaseq_rpkm <- reduced_rnaseq_rpkm[,matchedSampleId]

  		# Transpose the previous data.frame to bind with the data.frame that contains the selected filters (reduced_lib_header)
  		reduced_rnaseq_rpkmT <- t(reduced_rnaseq_rpkm)

  		# Bind data.frames
  		finalData <- cbind(reduced_rnaseq_rpkmT, reduced_lib_header)
  		
  		return(finalData)  		
  	})


  	# Split the final data and plot the boxplot
  	plotBoxplot <- function(){
		print("1")  		
  		# numberOfBoxplot represent the number of graph that will be ploted
  		if (filter()$plotEachFilter){
  			numberOfBoxplot <- length(gene()) * length(filter()$filter_2)
  		}
  		else{
  			numberOfBoxplot <- 	length(gene())
  		}

  		# Determine the number of the lines the area the boxplot will be must have
		par(mfrow = c(numberOfBoxplot,1))

		# Determine the colour of each boxplot if the class is subgroup as each ootions has a specific one
		color <- defineColors(filter()$filter_1)

		# Get the data to plot the boxplot
		plotData <- boxplotData()
  		
  		# Plot one boxplot for each inserted gene
  		if(!filter()$plotEachFilter){
  			for(i in 1:numberOfBoxplot){
  				dataP <- split(plotData[,i], plotData[[filter()$filters[1]]])
	  			
	  			# If the NA option was checked the it is appended to the data in order to be ploted	
  				if(length(which(filter()$filter_1 == "NA"))>0){
	  				matchedNAIds <- which(is.na(plotData[[filter()$filters[1]]]))
	  				values <- plotData[matchedNAIds,i]
	  				dataP <- c(dataP, list('NA' = values))
  				}
  				
	  			dataP <- dataP[filter()$filter_1]

	  			# If some checked options has 0 sample then the value of the options is set to 0 in order the appear at the plot.
		  		for(k in 1:length(dataP)){
		  			if (length(dataP[[k]]) == 0){
						dataP[k] <- 0
					}
		  		}

	  			title <- getHGNCsymbol(gene()[i])

  				boxplot(
  					dataP,
  					col = color,
  					frame = FALSE,
  					main = title,
  					notch = input$notch,
  					varwidth = input$varwidth,
  					las = 1
  				)
  			}
  		}
  		 # Plot one boxplot for each inesrted gene and each selected options that was chosen in filter 1
  		 else{
  		 	
  		 	partialData <- split(plotData, plotData[[filter()$filters[2]]])


  		 	for(i in 1:length(gene())){

  		 		for(j in 1:length(partialData)){
  		 			dataP <- split(partialData[[j]][,i], partialData[[j]][filter()$filters[1]])
  		 			
  		 			# If the NA option was checked the it is appended to the data in order to be ploted
  		 			if(length(which(filter()$filter_1 == "NA"))>0){
		  				matchedNAIds <- which(is.na(plotData[[filter()$filters[1]]]))
		  				values <- plotData[matchedNAIds,i]
		  				dataP <- c(dataP, list('NA' = values))
	  				}
	  				
		  			dataP <- dataP[filter()$filter_1]

		  			# If some checked options has 0 sample then the value of the options is set to 0 in order the appear at the plot.
		  			for(k in 1:length(dataP)){
		  				if (length(dataP[[k]]) == 0){
							dataP[k] <- 0
						}
		  			}
		  			
		  			title <- getHGNCsymbol(gene()[i])

	  				boxplot(
	  					dataP,
	  					col = color,
	  					frame = FALSE,
	  					main = title,
	  					notch = input$notch,
	  					varwidth = input$varwidth,
	  					las = 1
	  				)
  		 		}
  		 	}
  		 }
  	}

  	# Get the HGNC symbol from bioMart_table in order to get create the legend to the boxplot
  	getHGNCsymbol <- function(ensemblID){
  		
		HGNC_ID <- which (bioMart_table == ensemblID, arr.ind=TRUE)
		HGNC_symbol <- as.character(bioMart_table[HGNC_ID[1],3])
		
		return (HGNC_symbol)						
  	}

  	# Function responsible to determine the colors of the plots in the boxplot.
	# WNT - blue / SHH - red / Group3 - yellow / Group 4 - green / Control - gray
	defineColors<-function(filter){
		color <- NULL

		for (i in 1:length(filter)){
			if(filter[i] == "WNT"){
				color <- c(color, "blue")
			}
			else 
				if (filter[i] == "SHH"){
					color <- c(color, "red")
				}
				else
					if (filter[i] == "Group3"){
						color <- c(color, "yellow")
					}
					else
						if(filter[i] == "Group4"){
							color <- c(color, "green")
						}
						else
							color <- c(color, "gray")
		}

		return(color)
	}


	# Render boxplot at Boxplot tab
	output$plot <- renderPlot({
		plotBoxplot()
	})

	# Create a table at boxplot tab to show all the gene inserted by the user
	output$geneData <- renderTable({
		ensemblID_M <- as.matrix(gene())

		aux <- NULL
		for(i in 1:length(gene())){
			aux <- c(aux, getHGNCsymbol(gene()[i]))
		}


		HGNCsymbol_M <- as.matrix(aux)

		table <- cbind(HGNCsymbol_M, ensemblID_M)
		colnames(table) <- c("HGNC Symbol", "Ensembl Gene ID")

		table

	},include.rownames=FALSE)


	# Create a table at boxplot tab to show all the patient id if some was inserted
	output$patientidData <- renderTable({
		
		if(patient()$patientid == TRUE){

			
			if (length(patient()$patientidtext[[1]]) <= 7){
				aux <- length(patient()$patientidtext[[1]])

				table <- matrix(patient()$patientidtext[[1]], ncol=aux, byrow=TRUE)
			}
			else{
				aux <- patient()$patientidtext[[1]]

				for (i in (length(aux)%%7):6){
					aux <- c(aux, " ")
				}

				table <- matrix(aux, ncol=7, byrow=TRUE)
			}
		}
		else{
			data <- "NONE"
			table <- matrix(data)
			colnames(table) <- "Patient ID"
		}

		table
	},include.rownames=FALSE, include.colnames=FALSE)



	# Create a table at boxplot tab to show all the select filter and options
	output$filterData <- renderTable({
				
		table <- NULL

		if (filter()$filters[1] == "Age")
			max <- length(filter()$filter_1) - 1
		else
			max <- length(filter()$filter_1)

		if(length(filter()$filters) > 1){
			for (i in 2:length(filter()$filters)){

				if (filter()$filters[i] == "Age")
					aux <- length(filter()[[paste("filter_",i,sep="")]]) -1
				else
					aux <- length(filter()[[paste("filter_",i,sep="")]])

				if (max < aux){
					max <- aux
				}
			}
		}

		for (i in 1:length(filter()$filters)){

			if(filter()$filters[i] == "Age")
				aux <- filter()[[paste("filter_",i,sep="")]][-length(filter()[[paste("filter_",i,sep="")]])]
			else
				aux <- aux <- filter()[[paste("filter_",i,sep="")]]
				
			if(length(aux) < max)
				for(j in length(aux):(max-1)){
					aux <- c(aux, "")
				}

			table <- rbind(table, aux)
		}


		rownames(table) <- filter()$filters

		table
	},include.colnames=FALSE)



	#---------------------------------------------------------------------------------------------------------#
	#                                  MANAGE THE DATA TO CREATE A DATA TABLE                                 #
	#---------------------------------------------------------------------------------------------------------#

	# Manage the final data getting only the columns that matter to create a table at Data tab
	finalData <- reactive({
		order <- c("Tissue_ID", "Patient_ID", filter()$filters, toupper(gene()))
		data <- boxplotData()
		data <- data[order]

		return(data)
	})

	# Render the final data at Data tab
	output$dataTable <- renderDataTable(lib_header(), options = list(iDisplayLength = 10))

	#---------------------------------------------------------------------------------------------------------#
	#			                                  DOWNLOAD BUTTONS                                 			  #
	#---------------------------------------------------------------------------------------------------------#

	output$downloadBoxplotPDF <- downloadHandler(
		filename <- function() { paste('Boxplot.pdf') },
		content <- function(file) {
			pdf(file, width = input$myWidth/72, height = input$myHeight/72)
			## ---------------
			plotBoxplot()
			## ---------------
			dev.off()
		},
		contentType = 'application/pdf' # MIME type of the image
	)

	output$downloadDataPDF <- downloadHandler(
		filename <- function() { paste('data.pdf') },
		content <- function(file) {
			pdf(file, width = input$myWidth/72, height = input$myHeight/72)
			## ---------------
			finalData()
			## ---------------
			dev.off()
		},
		contentType = 'application/pdf' # MIME type of the image
	)


	output$downloadDataTXT <- downloadHandler(
    	filename = function() { "boxplotData.txt" },
   		content = function(file) {
		write.table(finalData(), file, row.names=FALSE)
    })


})

