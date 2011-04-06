/*****************************************************************************
This file is part of the Primegens Project.  

PUT COPYRIGHT INFO HERE

PUT VERSION INFO HERE
*****************************************************************************/



/*
 * This module contains following functionlities
 * 1- formatting various databases for using blast
 * 2- find unique oligos from the pair pairs 
 * 3- run blast for unique oligos and save results
 *
 */


/*formatting database for blast uses*/

/* 
Function: int format_chromosome_database(int bloop, int cloop)
Inputs: The number of blast loops and chromosome loops
Outputs: creates a virtual database from a chromosome database
THIS FUNCTION MAY CAUSE THE PROGRAM TO TERMINATE IF THE DATABASE FILES CANNOT BE OPENED
*/
int format_chromosome_database(int bloop, int cloop);

/* 
Function: int format_single_database(void)
Inputs: None
Outputs: creates a virtual database from a single database file
*/
int format_single_database(void);

/* 
Function: int format_multiple_database(void)
Inputs: None
Outputs: creates a virtual database from multiple database files
*/
int format_multiple_database(void);



/*find unique oligos from primer pairs*/

/* 
Function: int find_uniq_oligos(void)
Inputs: None
Outputs: Finds unique oligios from primer pairs and stores them in uoligos
*/
int find_uniq_oligos(void);

/* 
Function: int chk_if_exist(char *oligo, UOligo *tuoligo, int size)
Inputs: The oligo to check for, the temporary list of unique oligos, and the size of the list
Outputs: Returns the position of the oligo in the list, or -1 if not found
*/
int chk_if_exist(char *oligo, UOligo *tuoligo, int size);



/*run blast for unique oligos only*/

/* 
Function: int run_uo_blast_chromosome(int bloop, int cloop)
Inputs: The number of blast loops and chromosome loops
Outputs: Runs blast on the unique oligos from chromosomes and stores the results in the file stored in fn_bout
*/
int run_uo_blast_chromosome(int bloop, int cloop);

/* 
Function: int run_uo_blast_single(void)
Inputs: None
Outputs: Runs blast on the unique oligos from a single database file and stores the results in the file stored in fn_bout
*/
int run_uo_blast_single(void);

/* 
Function: int run_uo_blast_multiple(int dbcount)
Inputs: the number of database files
Outputs: Runs blast on the unique oligos from several database files and stores the results in the file stored in fn_bout
*/
int run_uo_blast_multiple(int dbcount);



/*save blast output for unique oligos only*/

/* 
Function: int save_blast_output_chromosome(int bloop, int cloop)
Inputs: the number of blast loops and chromosome loops
Outputs: Stores the output from the BLASTs from the chromosomes
*/
int save_blast_output_chromosome(int bloop, int cloop);

/* 
Function: int save_blast_output_single(void)
Inputs: None
Outputs: Stores the output from the BLASTs from a single database file
*/
int save_blast_output_single(void);

/* 
Function: int save_blast_output_multiple(int dbcount)
Inputs: The number of database files
Outputs: Stores the output from the BLASTs from multiple database files
*/
int save_blast_output_multiple(int dbcount);
