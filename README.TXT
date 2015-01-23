=============================================================================================================================
=															MAGIC															=
=																															=
=			Genomics browser for the Medulloblastoma Advanced Genomics International Consortium (MAGIC) project				=
=																															=
=============================================================================================================================

NOTE: the MAGIC portal is primarily coded for the display of Medulloblastoma data - for other datasets, please alter the server.R file accordingly. 

=====================================
=       REQUIRED R LIBRARIES        =
=====================================

There are a number of required R libraries for both the server.R and heatmap_2_revised.R code to function correctly, including:

	- shiny
		* To install, please follow instructions found here - http://rstudio.github.io/shiny-server/latest/#install-shiny
		
	- DESeq2
		* To install, please follow instructions found here - http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html
		
	- plotrix
		* Run install.packages('plotrix')
		
	- gridExtra
		* Run install.packages('gridExtra')
	
	- gplots (for heatmap.2_revised)
		* Run install.packages('gplots')
		
	- gtools (for heatmap.2_revised)
		* Run install.packages('gtools')
		
	- stats (for heatmap.2_revised)
		* Run install.packages('stats')
		
=====================
=       DATA        =
=====================

The data used in the application should be located in the Data folder. The files specified below are required for the application to run correctly.

Data Folder Contents:
	- bioMart_table.txt
		* Based on Ensembl database v75 for Home Sapiens data. The table contains all the Ensembl Gene Id, their Associated Gene Name, HGNC Symbol, and Entrez Gene ID.
		* The data should have the same headers that are supplied by BioMart, and must include at least headers of 'HGNC symbol', and 'Ensembl Gene ID'. 

	- lib_header.txt
		* A tab-delimited file containing samples and their associated information to filter on (i.e. the Subgroup type of the tissue sample). The default filters which the code will use are Subgroup, Age, Gender, and Tissue Type. The application can also filter on the ID of the patient for a more specific selection process. 

	- rnaseq_rpkm_matrix.txt
		* A matrix of the rnaseq value for each gene for each sample. The application requires that the column names are the tissue sample ids, and the row names are the Ensembl gene IDs. 
		i.e. the header line should be Gene_ID<TAB>Sample_ID1<TAB>Sample_ID2<TAB>Sample_ID3<TAB>...
			and all subsequent lines should be in the format:
				Ensembl_GeneID<TAB>RPKM_Value1<TAB>RPKM_Value2<TAB>...
		
	- rnaseq_raw_count_matrix.txt
		* A matrix containing the rnaseq raw count data for each gene for each sample -- will be used in DESeq normalization. The matrix format should be the same as was specified above for the RPKM matrix. 


=====================================
=  APPLICATION FOLDERS / KEY FILES  =
=====================================

- server.R : contains the bulk of the backend code for the application.

- heatmap_2_revised.R : contains a revised version of the heatmap.2 code which allows for > 1 side rows/columns on the heatmap
	|
	-- please see https://gist.github.com/nachocab/3853004 and http://tolstoy.newcastle.edu.au/R/e9/help/10/02/6212.html for details of the patch which was applied to get this file

- data: contains the above specified data filles 

- www: the index html page and all plugins, images and style which are used.
	|
	 -- css: all css files such as bootstrap, datatable and the style.
	|
	 -- js: all the javascript files such as boostrap, datatable, jquery, shiny and the functions.
	|
	 -- img: page logo and loading modal gif.



=================================
=  LOCATION OF THE APPLICATION	=
=================================

The application is located in the default found in which shiny-server look for it: /srv/shiny-server/

accessing the application:
	$ cd /srv/shiny-server/MAGIC


=====================
=      SERVER       =
=====================

The server used to run the application as a web page is a local install of a shiny-server. Shiny Server is a server program that makes Shiny applications available over the web. It’s free, open source, and available from Github.

Configuration:
	- The folder "/srv/shiny-server/" is set as default as the place shiny-server will look for the application.
	- The file "/etc/shiny-server/shiny-server.conf" provide the configuration for the server

	The complete manual to configure the server can be found at: http://rstudio.github.io/shiny-server/latest/
	Shiny-server page: https://github.com/rstudio/shiny-server

	The installer automatically deploys the necessary scripts to ensure that Shiny Server is started automatically on boot. The upstart is a system used to automatically start, stop and manage services. The installer  automatically deployed an Upstart script to /etc/init/shiny-server.conf. This script will initialize shiny-server as soon as the network is activated on the machine and stop when the machine is being shut down.
	The Upstart script also ensures that shiny-server is respawned if the process is terminated unexpectedly.

	To start the server:	
		$ sudo start shiny-server

	To stop the server:
		$ sudo stop shiny-server

	To restart the server:
		$ sudo restart shiny-server

	To reload the configuration while the server keeps running:
		$ sudo reload shiny-server

	To check the status or retrieve the Process ID associated with shiny-server:
		$ status shiny-server

=====================
=   QUICK SETUP     =
=====================

The steps for installing the MAGIC portal are as follows:

	1. Download the entire MAGIC folder locally
	2. Create a MAGIC/data folder and store the files specified in the 'DATA' section in the folder
	3. Install the required R libraries specified above onto the server you wish to set MAGIC up for
	4. Install a shiny server onto the same machine
	5. Change any of the shiny server defaults to suit your application (optional)
	6. Copy your entire MAGIC folder into the /srv/shiny-server directory
	7. Restart the shiny server 
	8. Access the MAGIC portal through a web browser