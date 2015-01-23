#required libraries
library(shiny)
library(DESeq2)
library(plotrix)
library(gridExtra)
library(gplots)

#required source file
source("heatmap_2_revised.R")

# Define server logic
shinyServer(function(input, output, session) {
  
  #----------------------------------------------------------------------------------------------------------#
  #                                            PSEUDO CONSTANTS                                              #
  #         These values will be replaced with simply "recurrence" for the filtering - more can be added.    #
  #----------------------------------------------------------------------------------------------------------#
  
  replaceWithRecurrence <- c("recurrentmetastasis", "localrecurrence2")

	#----------------------------------------------------------------------------------------------------------#
	# 								                       				READ DATA									                        			   #
	#----------------------------------------------------------------------------------------------------------#


	#--------------------------------------- READ DATA FROM FOLDER DATA ---------------------------------------#

	# Data loaded to get the number of samples for each option of each filter once it is selected.
  # Once the txt file has been read in (which takes a bit of time) it is stored as a .RData file 
  # to be read in more quickly for every other time.
	lib_header <- reactive({
    if(file.exists("data/lib_header.RData")) 
    {
      load("data/lib_header.RData")
    }
    else 
    {
      lib_header = read.table("data/lib_header.txt", sep="\t", header=TRUE)
      
      # make sure to replace the values with recurrence before saving it to /data
      for(i in 1:length(replaceWithRecurrence))
      {
        lib_header$Tissue_Type[lib_header$Tissue_Type == replaceWithRecurrence[i]] <- "recurrence"
      }
      
      lib_header$Tissue_Type <- factor(lib_header$Tissue_Type)
      
      save(lib_header, file="data/lib_header.RData")
    }
    
		return (lib_header)
	})
  
  # Returns a data frame containing genes and their associated rnaseq raw count data (for use in the deseq normalization)
  # To make the read more efficient, if the RData object has been made already, then the program will simply read that.
  deseq_raw <- reactive({
    
    if(file.exists("data/rnaseq_raw_matrix.RData")) {
      load("data/rnaseq_raw_matrix.RData")
    }
    else {
      raw_count = read.table("data/rnaseq_raw_count_matrix.txt", sep="\t", header=TRUE, row.names=1)
      
      #making sure to round up all values to integers (for the counts to be whole)
      raw_count[,] <-ceiling(raw_count[,])
      
      rnaseq_raw_init=as.matrix(raw_count)
      
      # remove rows containing 0 for all counts
      count_sum=apply(rnaseq_raw_init,1,sum)
      table(count_sum>0)
      rnaseq_raw_matrix=rnaseq_raw_init[-which(count_sum==0),]
      
      mode(rnaseq_raw_matrix) <- "integer"
      
      save(rnaseq_raw_matrix, file="data/rnaseq_raw_matrix.RData")
    }

    return(rnaseq_raw_matrix)
  })
  
  
  # runs deseq normalization on the raw count data and returns the results as a data frame
  rnaseq_deseq_norm <- reactive({
    # gets a matrix of the raw count values
    M_rnaseq_raw = deseq_raw()
    
    # this is just setting up for normalizing (not actually running the full DESeq pipeline)
    conds=factor(rep("MB",ncol(M_rnaseq_raw)))
    table(conds)
    
    nb_A=round(ncol(M_rnaseq_raw)/2)
    nb_B=ncol(M_rnaseq_raw)-nb_A
    
    colData=data.frame("condition"=c(rep("MB1",nb_A),rep("MB2",nb_B)))
    rownames(colData)=colnames(M_rnaseq_raw)
    
    dds <- DESeqDataSetFromMatrix(countData = M_rnaseq_raw, colData = colData, design = ~ condition)
    colData(dds)$condition <- factor(colData(dds)$condition,levels=c("MB1","MB2"))
    
    dds <- estimateSizeFactors(dds)
    M_rnaseq_raw_norm=counts( dds, normalized=TRUE )
    
    rnaseq_deseq_normalized <- data.frame(M_rnaseq_raw_norm)
    return (rnaseq_deseq_normalized)
  })

	# Returns a dataframe containing genes and their associated rnaseq rpkm values (for the rpkm plots)
	rnaseq_rpkm <- reactive({
    
    if(file.exists("data/rnaseq_rpkm.RData"))
    {
      load("data/rnaseq_rpkm.RData")
    }
    else 
    {
      rnaseq_rpkm = read.table("data/rnaseq_rpkm_matrix.txt", sep="\t", header=TRUE, row.names=1)
      save(rnaseq_rpkm, file="data/rnaseq_rpkm.RData")
    }
    
    return(rnaseq_rpkm)
	})

	# Load Ensembl data v62. Used to validated the inputted gene and to get the Ensembl gene ID to plot the graphs.
  bioMart_table <- reactive({
    if(file.exists("data/bioMart_table.RData"))
    {
      load("data/bioMart_table.RData")
    }
    else 
    {
      bioMart_table = read.table("data/bioMart_table.txt", sep="\t", header=TRUE)
      bioMart_table[bioMart_table==''] <- NA
      save(bioMart_table, file="data/bioMart_table.RData")
    }
    
    return (bioMart_table)
  })
  
  lc_bioMart_table <- reactive({
    bioMart <- as.data.frame(lapply(bioMart_table(), tolower))    
    return (bioMart)
  })


  #---------------------------------------- GET DATA FROM CLIENT SIDE ----------------------------------------#

	# Get the inserted gene(s). If there is more than one, split them into an array.
	# returns a Character vector containing the bioMart HUGO Gene symbols.
	gene <- reactive({

    # split the input string by both newline and space, and then make sure no entries were empty
		genes <- strsplit(tolower(input$geneInput), "[ \n]")
    genes[[1]] <- unique(genes[[1]][which(genes[[1]]!="")])

    # make sure that the biomart table used in comparisons is all lower case    
		ensemblIDArray <- NULL
		if(length(genes[[1]]) > 0) {
  		for(i in 1:length(genes[[1]])){
  			ensemblID <- which (lc_bioMart_table() == genes[[1]][i], arr.ind=TRUE)
  			ensemblIDArray <- c(ensemblIDArray, as.character(bioMart_table()[ensemblID[1],'Ensembl.Gene.ID']))
  		}
		}

		return (ensemblIDArray)
	})


	# Get the status of the field 'Patient ID'. The value of $patientid is determined by the checkbox value.
	# This reactive "function" also get the ids and split them into an array.
	patient <- reactive({
		patient <- list(patientid = input$patientid)
    
    # also allow for the patient ids to be separated by both newlines and spaces
		patient$patientidtext <- strsplit(input$patientidTextarea, "[ \n]")
		patient$patientidtext[[1]] <- unique(patient$patientidtext[[1]][which(patient$patientidtext[[1]]!="")])
		
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
	#										                       	OBSERVER FUNCTIONS 	                    										   # 
	#		    	Receive data from the client and automaticaly execute an action in accord to each data 	    	   #
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
					checkboxValue <- checkNA(checkboxValue)
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

	# Verify if the gene inputted by the user exists within the Ensembl database (data/bioMart_table.txt).
	# If exist return true, else return the gene(s) that do(es) not exist in our database;
	gene_input <- observe({

		if(input$geneInput == ""){
			session$sendCustomMessage(type = "verifyGeneInput", "BLANK")
		}
		else{
			genes <- strsplit(input$geneInput, "[ \n]")
			genes[[1]] <- unique(genes[[1]][which(genes[[1]]!="")])
			matchGene <- TRUE
			wrongGenes <- NULL
			
			if(length(genes[[1]] > 0)) {
  			for(i in 1:length(genes[[1]])){
  				gene <- tolower(genes[[1]][i])
  				if (length(which (lc_bioMart_table() == gene)) == 0){
  					matchGene <- FALSE
  					wrongGenes <- c(wrongGenes, genes[[1]][i])
  				}
  			}
  			
  			if(matchGene == TRUE)
  			  session$sendCustomMessage(type = "verifyGeneInput", "TRUE")
  			else
  			  session$sendCustomMessage(type = "verifyGeneInput", wrongGenes)
      }
		}
	},suspend = TRUE)

	# Verify that the Patient ID checkbox is checked. If it is, then verify if the inputted patient(s) 
	# id(s) exist within the patient list.
	patientid_input <-observe({
		
		if (input$patientid == TRUE){
			if(input$patientidTextarea == ""){
				session$sendCustomMessage(type = "verifyPatientId", "BLANK")
			}
			else{
				patientids <- strsplit(input$patientidTextarea, "[ \n]")
				patientids[[1]] <- unique(patientids[[1]][which(patientids[[1]]!="")])
				
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
	#					                    		MANIPULATE THE DATA TO SEND TO THE CLIENT SIDE						        		   #	
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
				order <- c("control", "metastasis", "null", "primary", "recurrence", "na")
    else
      if(filterName == "Age")
        order <- c("Infant", "Child", "Adult")
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
    if(data == "na")
      return (toupper(data))
    
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
	checkNA <- function(filterData){
    if(!("na" %in% filterData) && !("NA" %in% filterData))
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

  # Reduce the rnaseq_deseq_norm data. The reduced data frame contains only the genes inserted by the user.
  # It is made to reduce the time to process data
  r_rnaseq_deseq <- reactive({
    matchedGeneIds <- NULL
    
    deseq_norm <- rnaseq_deseq_norm()
    
    for(i in 1:length(gene())) {
      matchedGeneIds <- c(matchedGeneIds, toupper(gene()[i]))
    }
    
    # create a data frame with only the inputted genes
    reduced_rnaseq_deseq <- deseq_norm[matchedGeneIds,]
    
    return(reduced_rnaseq_deseq)
  })


	# Reduce the rnaseq_rpkm data. The reduced data frame contains only the genes inserted by the user.
	# It is made to reduce the time to process the data.
	r_rnaseq_rpkm <- reactive({
		matchedGeneIds <- NULL

  		# Get the number of the lines that correspond to each inserted gene
      rpkm <- rnaseq_rpkm()
    
      for(i in 1:length(gene())) {
	    	matchedGeneIds <- c(matchedGeneIds, toupper(gene()[i]))
      }
    
  		# Create a data.frame that contains only the inserted gene
  		reduced_rnaseq_rpkm <- rpkm[matchedGeneIds,]
    
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
		  			# the accord to the age pattern (american or european), i.e, convert each age to their correspondend
		  			# category (Infant, Child or Adult)
		  			
		  			if(filter()$filters[i] == "Age" ){
		  			  reduced_lib_header <- arrangeData(reduced_lib_header, filter()$agePattern);
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
  			for(i in 1:length(patient()$patientidtext[[1]])){
  				matchedFilterIds <- c(matchedFilterIds, which(tolower(reduced_lib_header$Patient_ID) == tolower(patient()$patientidtext[[1]][i])))
  			}
  			reduced_lib_header <- reduced_lib_header[sort(matchedFilterIds),]
  		}
  		
  		return (reduced_lib_header)
	})

	# Manage the data, i.e, create a data frame that contains the rnaseq value for the specific inserted
	# gene(s), patient(s) id(s) and selected filter and options.
  	boxplotData <- reactive({      
      reduced_rnaseq <- NULL
      if(input$dataSource == "rpkm") 
      {
  		  # Get the data only with the inserted genes
  		  reduced_rnaseq <- r_rnaseq_rpkm()
      }
      else #we will use the DESeq normalized data
      {
        reduced_rnaseq <- r_rnaseq_deseq()
      }
      
      
  		# Get the already filtered data
  		reduced_lib_header <- r_lib_header()
      
  		# Match tables obtaining a vector containg all colunms at reduced_rnaseq in which all userIds (LIB) match.
  		matchedSampleId <- match(reduced_lib_header$LIB, colnames(reduced_rnaseq))

  		# Create a new data.frame containing only the data selected by the user (gene, patient id and selected filters)
  		reduced_rnaseq <- reduced_rnaseq[,matchedSampleId]

  		# Transpose the previous data.frame to bind with the data.frame that contains the selected filters (reduced_lib_header)
  		reduced_rnaseqT <- t(reduced_rnaseq)

  		# Bind data.frames
  		finalData <- cbind(reduced_rnaseqT, reduced_lib_header)
      
  		return(finalData)  		
  	})


    # get the table filled with p values for display within the boxplot
    getPValueTable <- function(dataP, filters, zeroLengthSets){
      
      pValues <- NULL
      r_names <- NULL
      dataList <- list()
      
      setsToConsider <- setdiff(filters, zeroLengthSets)
      
      # for each of the different filters present (i.e. 'SHH', 'WNT', etc.), take the first 
      # group of data which can be filtered out (for the specific filter value) and then 
      # test it against each of the other groups present
      for(i in 1:length(setsToConsider))
      {        
        dataX <- dataP[[setsToConsider[i]]]
        
        if(i < length(setsToConsider))
        {
          for(j in (i+1):length(setsToConsider))
          {
              dataY <- dataP[[setsToConsider[j]]]
            
              # add the wilcox test between each of the two groups to the pValues vector
              pValues <- c(pValues, format(wilcox.test(dataX, dataY, exact=FALSE)$p.value, digits=3))
              r_names <- c(r_names, paste(setsToConsider[i], setsToConsider[j], sep=" - "))
          }
        }
        
        # store the data in a list for the kruskal test later
        dataList[[setsToConsider[i]]] <- dataX
      }
      
      # call kruskal.test() on all of the different vectors of data from each of the groupings present in the filter
      if(length(dataList)>1) {
        pValues <- c(pValues, format(kruskal.test(dataList)$p.value, digits=3))
        r_names <- c(r_names, "All groups (Kruskal-Wallis Rank Sum test)")
      }      
      
      pval_matrix <- as.matrix(pValues)
      rowname_matrix <- as.matrix(r_names)
      
      return(cbind(rowname_matrix, pval_matrix))
    }


  	# Split the final data and plot the boxplot
  	plotBoxplot <- function(plotPDF){
       tryCatch({ 
    		# numberOfBoxplot represent the number of the graph that will be ploted
    		if (filter()$plotEachFilter){
    			numberOfBoxplot <- length(gene()) * length(filter()$filter_2)
    		}
    		else{
    			numberOfBoxplot <- 	length(gene())
    		}

    		# Determine the number of the lines the area the boxplot will be must have for the web page
    		if(!plotPDF)
    		{
            par(mfrow = c(numberOfBoxplot,1))
    		}
        
    		# Determine the colour of each boxplot if the class is subgroup as each ootions has a specific one
    		color <- defineColors(filter()$filter_1)
    
    		# Get the data to plot the boxplot
    		plotData <- boxplotData()
    		
    		zeroLengthSets <- NULL
        
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
              
              setLength <- length(dataP)
              k <- 1
              
              ## TODO: refactor this into a separate method (to avoid later repeat)
    	  			# If some one of the checked filter options has 0 cases to plot
              # then remove the empty set from the list, such that it does not display
              # in the plot. Then make sure to handle the removal properly, and 
      				# add it to the list of filters to display within the warning message
      				while(k <= setLength)
              {
    		  			if (length(dataP[[k]]) == 0){
    		  			  zeroLengthSets <- c(zeroLengthSets, names(dataP)[k])
                  dataP[[k]] <- NULL
                  setLength <- setLength-1
    					  }
                else
                {
                  # if there still is data, check to see if we need to 
                  # convert it to log base 2 scale
      		  			if(input$logscale) {
                    
                    # for the rpkm data, make sure to add 1 to each value such that the
                    # zero values are handled accordingly
                    if(input$dataSource == "rpkm")
                      dataP[[k]] <- dataP[[k]]+1 
                    
      		  			  dataP[[k]] <- log2(dataP[[k]])
      		  			}
                  
      		  			k <- k + 1
                }
    		  		}                
      
    	  			title <- getHGNCsymbol(gene()[i])
              
      				boxplot(
      					dataP,
      					col = color,
      					frame = FALSE,
      					main = title,
                outline = input$outliers,
      					notch = input$notch,
      					varwidth = input$varwidth,
      					las = 1
              )
              
              if(length(filter()$filter_1) > 1) {
                
                session$sendCustomMessage(type = "displayPValueCheckbox", "TRUE")
                
                if(input$pvalues) {
                  tab <- (getPValueTable(dataP, filter()$filter_1, zeroLengthSets))
                  
                  tablePosition <- "topright"
                  
                  if(plotPDF)
                  {
                    plot.new()
                    tablePosition <- "center"
                  }                  
                  
                  addtable2plot(tablePosition,table=tab, title=paste("Pairwise p-values (Wilcoxon Rank Sum Test)", title, sep=" - "),bty="o", 
                                hlines=TRUE, vlines=TRUE, cex=1.15, lwd=1.5, display.rownames=FALSE, display.colnames=FALSE,
                               xpad=0.1, ypad=0.5)
                }
              }
              else
              {
                session$sendCustomMessage(type = "displayPValueCheckbox", "FALSE")
              }
      			}
      		}
      		 # Plot one boxplot for each inesrted gene and each selected options that was chosen in filter 1
      		 else{
      		 	
      		 	partialData <- split(plotData, plotData[[filter()$filters[2]]])
             
            # add in the NA list
      		 	if(length(which(filter()$filter_2 == "NA"))>0){
      		 	  matchedNAIds <- which(is.na(plotData[[filter()$filters[2]]]))
      		 	  values <- plotData[matchedNAIds,]
      		 	  partialData <- c(partialData, list('NA' = values))
      		 	}
             
            # remove all of the ones that aren't actually relevant to what was filtered
            partialData <- partialData[filter()$filter_2]
    
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
    
    	  				
    	  				setLength <- length(dataP)
    	  				k <- 1
    	  				
    	  				# If some one of the checked filter options has 0 cases to plot
    	  				# then remove the empty set from the list, such that it does not display
    	  				# in the plot. Then make sure to handle the removal properly, and 
    	  				# add it to the list of filters to display within the warning message
    	  				while(k <= setLength)
    	  				{
    	  				  if (length(dataP[[k]]) == 0){
    	  				    zeroLengthSets <- c(zeroLengthSets, names(dataP)[k])
    	  				    dataP[[k]] <- NULL
    	  				    setLength <- setLength-1
    	  				  }
    	  				  else
    	  				  {
    	  				    # if there still is data, check to see if we need to 
    	  				    # convert it to log base 2 scale
    	  				    if(input$logscale) {
    	  				      
    	  				      # for the rpkm data, make sure to add 1 to each value such that the
    	  				      # zero values are handled accordingly
    	  				      if(input$dataSource == "rpkm")
    	  				        dataP[[k]] <- dataP[[k]]+1 
    	  				      
    	  				      dataP[[k]] <- log2(dataP[[k]])
    	  				    }
    	  				    
    	  				    k <- k + 1
    	  				  }
    	  				} 
    		  			
    		  			title <- paste(getHGNCsymbol(gene()[i]), filter()$filter_2[j], sep=" - ", collapse=NULL)
                
    	  				boxplot(
    	  					dataP,
    	  					col = color,
    	  					frame = FALSE,
    	  					main = title,
    	  					notch = input$notch,
    	  					varwidth = input$varwidth,
                  outline = input$outliers,
    	  					las = 1
    	  				)
                
    	  				if(length(filter()$filter_1) > 1) {
                  
    	  				  session$sendCustomMessage(type = "displayPValueCheckbox", "TRUE")
                  
    	  				  if(input$pvalues) {
                    
                    tablePos <- "topright"
                    if(plotPDF) {
                      plot.new()
                      tablePos <- "center"
                    }
    	  				  
      	  				  tab <- (getPValueTable(dataP, filter()$filter_1, zeroLengthSets))
      	  				  
      	  				  addtable2plot(tablePos,y=NULL,table=tab, title=paste("Pairwise p-values (Wilcoxon Rank Sum Test)", title, sep=" - "),bty="o", 
      	  				                hlines=TRUE, vlines=TRUE, cex=1.15, lwd=1.5, display.rownames=FALSE, display.colnames=FALSE,
      	  				                xpad=0.1, ypad=0.5)
                  }
    	  				}
                else
                {
                  session$sendCustomMessage(type = "displayPValueCheckbox", "FALSE")
                }
                
      		 		}
      		 	}
      		 }        
    		
    		# either display or hide the warning message accordingly
    		if(length(zeroLengthSets) > 0){
          message <- paste("For ",length(unique(zeroLengthSets)), " group(s), there is no data after filtering. <br>Therefore, the following have been excluded from some of the plots : <b>", paste(unique(zeroLengthSets), collapse=", "), "</b>", sep="")
    		  session$sendCustomMessage(type = "displayWarningMessage", message)
    		}
      },
      warning=function(cond) {
        session$sendCustomMessage(type = "displayWarningMessage", conditionMessage(cond))
      },
      error=function(cond) {
        # TODO: fix the code such that it no longer has the "figure margins too large" error to begin with
        if(!grepl("figure margins too large", conditionMessage(cond)))
        {
          session$sendCustomMessage(type = "caughtError", conditionMessage(cond))
        }
      })
  	}

    plotHeatMap <- function(){
      # Get the data to plot the boxplot
      plotData <- boxplotData()
            
      # format both the data and the row names into a matrix
      data <- NULL
      rnames <- NULL
      
      for(i in 1:length(gene())) {
        data <- rbind(data, plotData[,i])
        rnames <- c(rnames, getHGNCsymbol(gene()[i]))
      }
      
      rownames(data) <- rnames
      
      # define the colours used for the bar ontop of the heat map, and also get the legend sorted
      colsidecols <- NULL
      cnames <- NULL
      
      for(j in 1:length(filter()$filters)) {
        colors <- defineHeatMapColors(as.character(plotData[[filter()$filters[j]]]))
        colsidecols <- rbind(colsidecols, colors)
        cnames <- c(cnames, filter()$filters[j])
      }
      
      rownames(colsidecols) <- cnames  
     
      keyXLabel <- paste(toupper(input$dataSource), "Values", sep=" ")
      
      scale <- "none"
      
      if(input$scaleRow) scale <- "row"
      
      # plot both the heatmap and the legend
      heatmap.2_revised(data, col="bluered", ColSideColors=t(colsidecols), trace="none", 
                        labCol=rep("", ncol(data)), key.xlab=keyXLabel, key.title=NA, 
                        keysize=1, cexRow=1.5, margins=c(0,8), scale=scale, NumColSideColors=length(filter()$filters))
    }

    plotHeatMapLegend <- function() {

      titles <- NULL
      categories <- list()
      colors <- list()
      
      for(j in 1:length(filter()$filters)) {
        titles <- c(titles, filter()$filters[j])
        
        
        # get the values for the legend
        categories <- c(categories, list(filter()[[paste("filter_",j,sep="")]]))
        colors <- c(colors, list(defineHeatMapColors(filter()[[paste("filter_",j,sep="")]])))
      }
      
      nc <- length(titles)
      if(nc > 3)
        nc <- 3
      
      nr <- c(ceiling(length(titles)/nc))
      
      par(mar=c(0, 0, 0, 0), mfrow=c(nr, nc))  
      
      for(k in 1:length(categories))
      {
        plot.new()
        ncol <- 3
        if(length(categories[[k]]) < ncol)
          ncol <- length(categories[[k]])
        
        legend("center", legend=categories[[k]], col=colors[[k]], lty=1, lwd=10, title=titles[k], xpd=TRUE, ncol=ncol)
      }
    }

    # Will output the formatted error message for the heatmap page (when there are too few genes selected)
    getHeatmapErrorMessage <- function() {
      currentGene <- gene()[1] # there should be only one gene
      message <- paste("No heatmap could be generated as only one gene - ", getHGNCsymbol(currentGene), " - was inputted. Please input at least 2 genes for heatmaps.")
      return(message)
    }
    
  	# Get the HGNC symbol from bioMart_table in order to get create the legend to the boxplot
  	getHGNCsymbol <- function(ensemblID){
  		HGNC_ID <- which (lc_bioMart_table() == tolower(ensemblID), arr.ind=TRUE)
  		HGNC_symbol <- as.character(bioMart_table()[HGNC_ID[1],'HGNC.symbol'])
  		
  		return (HGNC_symbol)						
  	}

    

  	# Function responsible to determine the colors of the plots in the boxplot.
	# WNT - blue / SHH - red / Group3 - yellow / Group 4 - green / Control - gray
	defineColors <- function(filter){
	  color <- NULL

		for (i in 1:length(filter)){
      color <- c(color, switch(filter[i],
                      "WNT"="blue",
                      "SHH"="red",
                      "Group3"="yellow",
                      "Group4"="green",
                      "gray"))
		}
    
		return(color)
	}

  defineHeatMapColors <- function(filter){
    color <- NULL
    # TODO: still have to define the gradient green colour scheme for age ranges
    for (i in 1:length(filter)){
      color <- c(color, switch(filter[i],
                               "WNT"="blue",
                               "SHH"="red",
                               "Group3"="yellow",
                               "Group4"="green",
                               "NOS"="gray45",
                               "M"="skyblue3",
                               "F"="pink",
                               "primary"="cyan",
                               "metastasis"="orange",
                               "recurrence"="mediumorchid",
                               "NA"="white",
                               "na"="white",
                               "Infant"="palegreen",
                               "Child"="limegreen",
                               "Adult"="darkgreen",
                               "gray"))
    }
    
    return(color)
  }


	# Render boxplot at Boxplot tab
	output$plot <- renderPlot({
		plotBoxplot(F)
    
    # sending a custom message to the client side to allow for the hiding of the load modal
		session$sendCustomMessage(type = "hideBoxplotModal", "boxplot")
	})

  # render the heatmap on the Heatmap tab
  output$heatmapPlot <- renderPlot({
    plotHeatMap()
    session$sendCustomMessage(type="hideHeatmapModal", "heatmap")
  })

  output$heatmapLegend <- renderPlot({
    plotHeatMapLegend()
  })

  output$heatmapErrorMessage <- renderText({
    session$sendCustomMessage(type="hideHeatmapModal", "heatmapError")
    getHeatmapErrorMessage()
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
	  session$sendCustomMessage(type = "hideBoxplotModal", "tableData")
	  
	  return(table)
  },include.rownames=FALSE)


  output$hm_geneData <- renderTable({
    ensemblID_M <- as.matrix(gene())
    
    aux <- NULL
    for(i in 1:length(gene())){
      aux <- c(aux, getHGNCsymbol(gene()[i]))
    }
    
    
    HGNCsymbol_M <- as.matrix(aux)
    
    table <- cbind(HGNCsymbol_M, ensemblID_M)
    
    colnames(table) <- c("HGNC Symbol", "Ensembl Gene ID")
    
    session$sendCustomMessage(type = "hideHeatmapModal", "tableData")
    
    return(table)
  },include.rownames=FALSE)


	# Create a table at boxplot tab to show all the patient id if some was inserted
	output$patientidData <- renderTable({
	  session$sendCustomMessage(type = "hideBoxplotModal", "tableData")
	  patientTable(3)
	},include.rownames=FALSE, include.colnames=FALSE)

  output$hm_patientidData <- renderTable({
    session$sendCustomMessage(type = "hideHeatmapModal", "tableData")
    patientTable(3)
  }, include.rownames=FALSE, include.colnames=FALSE)

	# Create a table at boxplot tab to show all the select filter and options
	output$filterData <- renderTable({
	  session$sendCustomMessage(type = "hideBoxplotModal", "tableData")
		filterTable()
	},include.colnames=FALSE)

  output$hm_filterData <- renderTable({
    session$sendCustomMessage(type = "hideHeatmapModal", "tableData")
    filterTable()
  }, include.colnames=FALSE)

  output$dataSourceData <- renderTable({
    table <- NULL
    aux <- c(toupper(input$dataSource))
    table <- rbind(table, aux)
    session$sendCustomMessage(type = "hideBoxplotModal", "tableData")
    table
  }, include.colnames=FALSE, include.rownames=FALSE)

  output$hm_dataSourceData <- renderTable({
    table <- NULL
    aux <- c(toupper(input$dataSource))
    table <- rbind(table, aux)
    session$sendCustomMessage(type = "hideHeatmapModal", "tableData")
    table
  }, include.colnames=FALSE, include.rownames=FALSE)


	#---------------------------------------------------------------------------------------------------------#
	#                                  MANAGE THE DATA TO CREATE DATA TABLES                                  #
	#---------------------------------------------------------------------------------------------------------#
  
  geneTable <- function() {
    ncols <- 6
    aux <- NULL
    geneCount <- length(gene())
    for(i in 1:geneCount){
      aux <- c(aux, getHGNCsymbol(gene()[i]))
    }
    
    excess <- geneCount%%ncols
    if(excess > 0 && (geneCount > ncols)) {
      for(i in excess:(ncols-1)) {
        aux <- c(aux, "")
      }
    }
    
    cols <- ncols
    
    if(geneCount < ncols)
    {
      cols=geneCount
    }
    
    table <- matrix(aux, ncol=cols, byrow=TRUE)
    colnames(table) <- c("HGNC Symbol(s)", rep("", (cols-1)))
    
    return(table)
  }

  filterTable <- function(){
    table <- NULL

    max <- length(filter()$filter_1)
    
    if(length(filter()$filters) > 1){
      for (i in 2:length(filter()$filters)){

        aux <- length(filter()[[paste("filter_",i,sep="")]])
        
        if (max < aux){
          max <- aux
        }
      }
    }
    
    lh <- r_lib_header()
    rnames <- NULL
    
    for (i in 1:length(filter()$filters))
    {
      aux <- NULL
      
      if(filter()$filters[i] == "Age")
      {
        rnames <- c(rnames, paste(filter()$filters[i], " [", filter()$agePattern, "]", sep=""))
        numberedAgeGroups <- getNumberOfSampleAge2(lh)
        
        selectedFilters <- filter()[[paste("filter_", i, sep="")]]
        
        # go through each of the age groups selected and match them to their corresponding counts
        for(j in 1:length(selectedFilters))
        {
          k <- switch(selectedFilters[j],
                 "Infant"=1,
                 "Child"=2,
                 "Adult"=3,
                 "NA"=4)

          if(filter()$agePattern == "European")
            k <- k+4
          
          aux <- c(aux, paste(selectedFilters[j], paste("(", numberedAgeGroups[k], ")", sep=""), sep=" "))
        }
      }
      else
      {
        rnames <- c(rnames, filter()$filters[i])
        aux <- getNumberOfSample(filter()[[paste("filter_", i, sep="")]], filter()$filters[i], lh)
      }
      
      if(length(aux) < max)
        for(j in length(aux):(max-1)){
          aux <- c(aux, "")
        }
      
      table <- rbind(table, aux)
    }
    
    rownames(table) <- rnames
    colnames(table) <- c("Filter Values", rep("", ncol(table)-1))
    
    return(table)
  }

  patientTable <- function(ncols){
    
    if(patient()$patientid == TRUE){
      aux <- patient()$patientidtext[[1]]
  
      excess <- length(aux)%%ncols
      if(excess > 0) {
        for(i in excess:(ncols-1)) {
          aux <- c(aux, "")
        }
      }
      
      table <- matrix(aux, ncol=ncols, byrow=TRUE)
      colnames(table) <- c("Patient ID(s)", rep("", (ncols-1)))
    }
    else{
      data <- "NONE"
      table <- matrix(data)
      colnames(table) <- "Patient ID(s)"
    }
    
    return(table)
  }

  optionsTable <- reactive({
    table <- NULL
    
    cnames <- c("Data Source", "Log2 Scale", "Outliers", "Variable Width", "Notches", "p Values")
    values <- c(input$dataSource, input$logscale, input$outliers, input$varwidth, input$notch, input$pvalues)
    table <- matrix(values, ncol=length(values), byrow=TRUE)
    
    colnames(table) <- cnames
    return (table)
  })
  
  heatMapOptionsTable <- reactive({
    table <- NULL
    cnames <-c("Data Source", "Scaled By Row")
    values <- c(input$dataSource, input$scaleRow)
    table <- matrix(values, ncol=length(values), byrow=TRUE)
    
    colnames(table) <- cnames
    return (table)
  })

	# Manage the final data getting only the columns that matter to create a table at Data tab
	finalData <- reactive({
	  order <- c("Tissue_ID", "Patient_ID", filter()$filters, toupper(gene()))
		data <- boxplotData()
		data <- data[order]

		return(data)
	})


	# Render the final data at Data tab
	output$dataTable <- renderDataTable({
  	  validate(
  	    need(nrow(r_lib_header()) > 0, "The filters selected leaves no data to graph")
  	  )
      r_lib_header()
      }, options = list(iDisplayLength = 10)
  )

	#---------------------------------------------------------------------------------------------------------#
	#			                                        DOWNLOAD BUTTONS                                     			  #
	#---------------------------------------------------------------------------------------------------------#

	output$downloadBoxplotPDF <- downloadHandler(
		filename <- function() { paste(input$currentDateTime, 'boxplot.pdf', sep="") },
		content <- function(file) {
			pdf(file, paper='A4', onefile=T)
			## ---------------
      # plot the boxplots
			plotBoxplot(T)
      
      # plot the table of selection information
			grid.arrange(tableGrob(filterTable()),  tableGrob(geneTable()), tableGrob(optionsTable()), tableGrob(patientTable(4)), ncol=1,as.table =TRUE)
			## ---------------
			dev.off()
		},
		contentType = 'application/pdf' # MIME type of the image
	)


	output$downloadDataTXT <- downloadHandler(
    	filename = function() { paste(input$currentDateTime, "boxplotData.txt", sep="") },
   		content = function(file) {
	    	write.table(finalData(), file, row.names=FALSE)
    })


  output$downloadHeatMapPDF <- downloadHandler(
    filename <- function() { paste(input$currentDateTime, 'heatmap.pdf', sep="") },
    content <- function(file) {
      pdf(file, paper='A4', onefile=T)
      ## ---------------
      # plot the heatmap/legend on separate pages
      plotHeatMap()
      plotHeatMapLegend()
      # plot the table of selection information
      grid.arrange(tableGrob(filterTable()),  tableGrob(geneTable()), tableGrob(heatMapOptionsTable()), tableGrob(patientTable(4)), ncol=1,as.table =TRUE)
      ## ---------------
      dev.off()
    },
    contentType = 'application/pdf' # MIME type of the image
  )
})

