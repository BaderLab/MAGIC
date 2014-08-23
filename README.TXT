=============================================================================================================================
=															MAGIC															=
=																															=
=			Genomics browser for the Medulloblastoma Advanced Genomics International Consortium (MAGIC) project				=
=																															=
=============================================================================================================================

=====================
=		DATA		=
=====================

The data used in the application are located at the folder Data.

Content:
	- bioMart_table.txt
		* Based on Ensembl database v62 for Home Sapiens data. The table contains all the Ensembl Gene Id and their Associated Gene Name, HGNC Symbol, UCSC Id and Refseq Ids (DNA, Predicted DNA, Protein, Predicted Protein, and Genomic).
		* All the data should be in lower case since the data is treat in the server to allow the user insert capitalized or uncapitalized texts.
		* In case the data have to be changed, a way to uncapitalize all the text of the exported data from Ensembl database is to use the following linux shell command: dd if=input.txt of=output.txt conv=lcase

	- lib_header.txt
		* Contains samples (LIB) with their id (Patient ID), Subgroup, Tissue ID, Tissue Type, Age and Gender

	- rnaseq_rpkm_matrix.txt
		* Contains the rnaseq value for each gene for each sample.


=========================
=  APPLICATION FOLDERS	=
=========================

- data: contains the above data (bioMart_table.txt, lib_header.txt and rnaseq_rpkm_matrix.txt)

- www: the index page and all used plugins, images and style.
	|
	 -- css: all css files such as bootstrap, datatable and the style created by me.
	|
	 -- js: all the javascript files such as boostrap, datatable, jquery, shiny and the functions created by me.
	|
	 -- img: page folder and loading modal gif.



=================================
=  LOCATION OF THE APPLICATION	=
=================================

The application is located in the default found in which shiny-server look for it: /srv/shiny-server/
To access the application is necessary to ssh to the virtual machine and drive through the directories to the folder in which the application is located.

ssh the virtual machine:
	$ ssh <username>@172.17.16.43

accessing the application:
	$ cd /srv/shiny-server/MAGIC


=====================
=		SERVER		=
=====================

The server used to run the application as web page is shiny-server. Shiny Server is a server program that makes Shiny applications available over the web. It’s free, open source, and available from Github.

Configuration:
	- The folder "/srv/shiny-server/" is set as default as the place shiny-server will look for the application.
	- The file "/etc/shiny-server/shiny-server.conf" provide the configuration for the server

	The complete manual to configure the server can be found at: http://rstudio.github.io/shiny-server/latest/
	Shiny-server page: https://github.com/rstudio/shiny-server

	The installer automatically deploy the necessary scripts to ensure that Shiny Server is started automatically on boot. The upstart is a system used to automatically start, stop and manage services. The installer  automatically deployed an Upstart script to /etc/init/shiny-server.conf. This script will initialize shiny-server as soon as the network is activated on the machine and stop when the machine is being shut down.
	The Upstart script also ensure that shiny-server is respawned if the process is terminated unexpectedly

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

	Before running the these commands you should ssh to the virtual machine:
		$ ssh <username>@172.17.16.43



=====================
=  AUTHENTICATION 	=
=====================

The authentication system is made by the Apache SSL encryption and basic HTTP authentication. Once the authentication is valid, the user is redirected to the MAGIC application.

To set up a username/password:
	$ sudo htpasswd /etc/httpd/htpasswd <username>

Before running the command, is necessary to ssh to the virtual machine.

=====================
=  REMOVE AN USER	=
=====================

To remove an user is necessary to ssh to the virtual machine and edit htpasswd file in order and delete the line that correspond to the user.

ssh the server:
	$ ssh <username>@172.17.16.43

edit the file (in this examples I am using the vim text editor but any other text editor can be used, such as nano):
	$ sudo vim /etc/httpd/htpasswd


=================================
=   ACCESSING THE APPLICATION 	=
=================================

Once the shiny-server is running, anyone that is connect to the network can access the application by the url http://https://172.17.16.43/MAGIC/





