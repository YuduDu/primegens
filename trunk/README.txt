
EXECUTION COMMAND:
PRIMEGESv2.0.5 has simple command for any type of case. All the complexity 
is contained into configuration file.

Command line:
$>PRIMEGENSv2.0.5.exe -q <query filename> -d <db filename> -p <db path>

In case of multilple databases (not developed yet), following command should
be used.
$>PRIMEGENSv2.0.5.exe -q <query filename> -l <db list filename> -p <db path>


			What is PRIMEGENSv2 ?
			---------------------------------
	(PRIMEGENS with additional feature and graphical user interface)

PRIMEGENS (PRIMEr Design Using GEN Specific Fragments) is a computer program to select
gene-specific fragments and then design primer pairs for PCR amplifications on high throughput scale.


			Acknowledgment
			-----------------------

The development of PRIMEGENS was sponsored by the Office of Biological and Environmental 
Research, U.S. Department of Energy, under Contract DE-AC05-00OR22725, managed by 
UT-Battelle, LLC. Dong Xu, Gary Li, Liyou Wu, Jizhong Zhou, and Ying Xu at Oak Ridge National 
Laboratory are involved in this project. We would like to thank our colleagues Manesh Shah, Joel A. 
Klappenbach, and Julia Stair for their help and insightful comments.


			Installation Instructioins
			-------------------------------

	Linux Version
	------------------

1)  After successfully obtaining primegens.tar.gz, copy the file into an installation directory, e.g. /usr/src/

2)  Unpack the compressed tar file:
			tar -xvzf primegens.tar.gz

3)  It will generate a directory called "PRIMEGENS", and all the related files are saved in the directory.

4)  Define the environments ($FINDSEG_PATH ) and put them in .cshrc or .login file, e.g.,
	setenv FINDSEG_PATH /usr/local/primegens/
	#(the path of the primegens package)

5)  Define an alias for the excutable and put it in .cshrc or .login file:
	alias primegens '/usr/local/primegens/bin/primegens.sun'


	Windows 2000/XP
	------------------------

1)	Download PRIMEGENS.zip file from the URL and save it on your PC
	at any location of your interest.

2)	Set Environment Variable for PRIMEGENS location on your computer.

	1) Go to CONTROL PANEL of your computer by clicking on START menu.
	2) Double click on SYSTEM icon. a new window SYSTEM PROPERTIES is open.
	3) Click on ADVANCED tab located on top in system properties window.
	4) Click on ENVIRONMENT VARIABLES located on bottom in advanced tab.
	5) Click NEW to add new environment variable under SYSTEM VARIABLES title.
	6) Type "FINDSEG_PATH" in Variable Name field.
	     and location of main PRIMEGENS directory (including itself) in Variable Value field
		for example type "C:\software\PRIMEGENS" as variable value if
		the PRIMEGENS directory is located in "C:\software" path.

IMPORTANT !
	few things are very essential to be certain.
	1) spaces are not allowed for the FINSDEG_PATH value. therefore be careful about 
	     location of PRIMEGENS directory. [c:\PRIMEGENS is good option]

	2) the FINDSEG_PATH value should be such that README.txt should be located 
	    in that path. 
	    By default winzip wizard will create new directory with name PRIMEGENS. in that case 
	    you will see README file in <PRIMEGENS.zip directory>\PRIMEGENS\PRIMEGENS.  
	
	3) After setting environment variables, please open a new MS Command Prompt in order
	    to get the FINDSEG_PATH activated.

3)	For successful testing of PRIMEGENS
	1) Open a MS Command Prompt on your computer.
	2) Goto TEST directory inside the PRIMEGENS directory.
	3) type "..\bin\PRIMEGENS.exe database.txt. once successfully executed,
	4) Type "..\bin\PRIMEGENS.exe -lf subset.txt database.txt".
	if the display shows "Job Finished Successfully !" means the program is ready to use.


Running Instruction 
-------------------------

Before running PRIMEGENS, select a location where you, want to store all PRIMEGENS out files.
Copy your database in that location also copy subset file according to your option. o

INPUT: Database File

The user need to specify the primer design requirement by selecting the DNA sequence whose primer 
pair need to be designed. the content of database should look like this -

   >TC216017
   GGCACGAGGAGATGGCTGAAGAGACAGTGAAAAGAA ……..
   >TC216017
   GGCACGAGGAGATGGCTGAAGAGACAGTGAAAAGAA ……...
   ACGACCATCACCCCTGCGTCGTGTGCCAGGCCANNTN ………
   NGAGGACGGCGACCAGGTCCGCATGCTCCCCTGCCGC ………
   >TC216017
   GGCACGAGGAGATGGCTGAAGAGACAGTGAAAAGAA ……....
   CACGTCTTCCACCGCCGCTGCTTCGACGGCTGGCTCCA ……...  
   >TC216069
   CAATNNNTCCNCCACCACCACGCCGGCGCCGGCGGCC ………
   CCACTACAAGTTTAACTGCCCACTCTGCCGCTCACCGC ………
   TCTTCTCCGACGAGCGCGTGGCCGTCACAGAGCGCCGC ……..


INPUT: Subset File

The subset file is an optional feature in PRIMEGENS, which specifies all sequence names amongst 
those present in database, whose primer need to be designed. In case when user wants primer pairs 
for all the sequences in database, there is no need for selecting the subset option. Once subset file is 
selected user has option to maintain FASTA format or not. Accordingly the user needs to select 
option for PRIMEGENS to read the subset input file. Following is the example for how subset file 
looks like.

1. )  FASTA Format Subset File Input
         	>TC216069
         	>TC216017
         	>TC216017
       	---------------
       	---------------

2.)  non-FASTA Format Subset File Input
        	TC216069
        	TC216017
        	TC216017
       	-------------
       	-------------	


INPUT: Command Line Options ( For CONSOLE application )
	
PRIMEGENS provide various options for the user to select according to his choice. Each option has 
its own feature in terms of extra modified processing with additional results for analysis. Here is a 
brief summary of the various user optional flags.

	• If the subset file is in a list (non-FASTA format) e.g.
		$> PRIMEGENS.exe -l  <subset-filename>  <database-filename>

	• If the subset file is in FASTA format e.g.
		$> PRIMEGENS.exe -lf  <subset-filename>  <database-filename>

	• If sequences in database file has no functions defined, one can use a separate function 
	  file (e.g., "function") in a format:  e.g.
		ORF1 function a
		ORF2 function b
		..................
		$> PRIMEGENS.exe  -name  <function-filename>  <database-filename>

	• Keep the BLAST outputs in the temporary directory (default: no) e.g.
		$> PRIMEGENS.exe -bf  <database-filename > 

	• Define the maximum expectation value (default: 1e-15) to get a BLAST alignment e.g.
		$> PRIMEGENS.exe -e 1e-10  <database-filename >

	• Define the minimum length (default: 100 bases) to select a unique segment. e.g.
		$> PRIMEGENS.exe -mseg 80  <database-filename >

	• Define the maximum sequence similarity (default: 75%) between the select segment 
	   and any other sequence. e.g.
		$> PRIMEGENS.exe -msim  <database-filename >

	• Define the minimum length (default: 50 bases) to filter out a BLAST alignment as a 
	   potential primer. e.g.
		$> PRIMEGENS.exe -mfil  <database-filename >

	• Define the faction (default: 0.8) of the length as the lower limit in the 
	  "PRIMER_PRODUCT_SIZE_RANGE" parameter for Primer3. e.g.
		$> PRIMEGENS.exe -f  0.6  <database-filename >
	
	• Define  "PRIMER_PRODUCT_SIZE_RANGE" independent of sequence length. e.g.
		 $> PRIMEGENS.exe -fz  80  120  <database-filename >

	• Calculate the unique sequences only without designing primers using Primer3. e.g.
		$> PRIMEGENS.exe -np  <database-filename >

	• Define the maximum length (-maxsep, default: 10000 bases) or the minimum length 
	  (-minsep, default: 0) for the primer product. e.g.
		$> PRIMEGENS.exe -f  -maxsep  1200  -minsep  100  <database-filename >

NOTE:-
You can combine several options together. The order of options does not matter as long as the 
database file is at the end of the argument for PRIMEGENS in case of command line execution. 

For graphical version of PRIMEGENS, creation of file is same as above but instead of giving 
command line option, we will have flexibity to choose any of the given options from the option
window. Here also you can combine several options togather.

  
OUTPUT : Various result files

PRIMEGENS software model is developed in a way which supports permanent storage of primer 
design for any experiment. In order to generate organized results, the model creates various files and 
directories with specific contents in those files. Here are brief description of all types of generated files 
and information contents. 

NOTE that for simplification, it is ASSUMED that user has selected database.txt as database file and 
subset.txt as subset file input to PRIMEGENS.

	• Database.txt_nohit.txt
	  This file contains all those query sequences from subset.txt, which are unique in 
	  database.txt. The content of the file is the sequence name, sequence length and the 
	  exact DNA sequence.

	• Database.txt_nohit.txt_back.txt
	  This file is the copy of the database.txt_nohit.txt file contents of previous experiment 
	  results executed and stored by PRIMEGENSv2.

	• Database.txt_nohit.txt_primer.txt
	  This file contains all the unique query sequences, whose primer is successfully designed 
	  by the PRIMEGENSv2. The file also contains the left and right primer sequences along 
	  with their index number on the main query sequence. 

	• Database.txt_nohit.txt_primer.txt_back.txt
	  This file contains all the record same as of Database.txt_nohit.txt_primer.txt for the 
	  previous experimental primer design results from created by PRIMEGENSv2.

	• Database.txt_primer.xls
	  This is an excel sheet file, which contains all designed query –specific primers along 
	  with various other information about the primer pairs. 

	• Database.txt_primer.xls_back.txt
	  This file contains the same contents as of Database.txt_primer.xls for the previous 
	  primer design results created by PRIMEGENSv2.

	• Database.txt_primer_undo.xls
	  This is an excel sheet file, which contains all the query sequences, whose primer pairs 
	  could not be designed by PRIMEGENSv2.

	• Database.txt_primer_undo.xls_back.txt
 	  This file contains all the information same as Database.txt_primer_undo.txt but for 
	  the previous results created by PRIMEGENSv2.

 	• Database.txt_seg.txt
	  This file contains all the query sequences, which are not unique in the database. But 
	  still have some that sequence-specific fragment, which are unique in the whole database. 
	  The file contains the name of the query sequence along with the sequence-specific 
	  fragments and its location on the main query sequence.

	• Database.txt_seg.txt_back.txt
	  This file contains same information as of Database.txt_seg.txt for the previous results 
	   created by the PRIMEGENSv2 software.

 	• Database.txt_seg.txt_primer.txt
	  This file contains all those query sequence, which are present in Database_seg.txt and 
	   their primer design is successful. The file contains the query-specific fragment along 
	   with the designed primer pairs and its location on the original sequence.  

	• Database.txt_seg.txt_primer.txt_back.txt
	  This file contains information similar to the Database.txt_seg.txt_primer.txt for the 
	   previous results generated by PRIMEGENS.

	• Database.txt_sim.txt
	  This file contains information about the query sequences, which are not unique in the 
	  database. This file shows which sequence in the database is closest to the query 
	  sequence with how much similarity.

	• Database.txt_sim.txt_back.txt
	  This file is same as Database.txt_sim.txt for the previous results generated by the 
	  PRIMEGENS.

	• Primer_plate_left [index]
	  This file contains query sequence and its left primer which should be kept in the primer 
	   plate number represented by the index digit.
 
	• Primer_plate_right [index]
	  This file contains query sequence and its left primer which should be kept in the primer 
	  plate number represented by the index digit.

Beside the various files there are some temporary directories which are created by PRIMEGENS, 
unlike files these directories are cleaned before the software execution. The information contents in 
these directories are described below.

	• TMP DIRECTORY
	  The storage of this directory is optional and is set according to the user input. 
	  This directory contains the BLAST output for each query sequences, which are 
	   used by PRIMEGENSv2 to design primer pairs. The file name is selected on the 
	   basis of query sequence, whose BLAST output is store in this file.

	• LOGS DIRECTORY
	  This file contains various logs, and global alignments for all those query sequence, 
	  which are not unique and global alignment were performed by PRIMEGENSv2 to 
	  select the sequence-specific fragment for that query sequence.
	  This directory is cleaned once the user executes the software again. 
	 


	
 
