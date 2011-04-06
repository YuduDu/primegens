/*****************************************************************************
This file is part of the Primegens Project.  

PUT COPYRIGHT INFO HERE

PUT VERSION INFO HERE
*****************************************************************************/


/*
 * this module contains all functions related to reading input
 * either from the command line or from files
 */



/*
Function: short init_param(void)
Inputs: None
Outputs: Initialize the parameters used throughout the rest of the program.  Returns pass if successful, fail otherwise.
*/
short init_param(void);

/*
Function: short get_param (int argc, char * argv[])
Inputs: the number of command line parameters, and the command line parameters themselves
Outputs: parses the command line parameters and updates the global variables.  Returns pass if successful, fail otherwise.
THIS FUNCTION WILL CAUSE THE PROGRAM TO TERMINATE IF THE STATUS FILE CANNOT BE OPENED
*/
short get_param (int argc, char * argv[]);

/*
Function: short set_default_param(void)
Inputs: None
Outputs: parses the command line parameters and updates the global variables.  Returns pass if successful, fail otherwise.
THIS FUNCTION WILL CAUSE THE PROGRAM TO TERMINATE IF THE OUTPUT FILE OR THE LOG FILE CANNOT BE OPENED
*/
short set_default_param(void);

/*
Function: short read_config_file(void)
Inputs: None
Outputs: Looks for a config.txt file in the output location, the current directory, then the primegens package, opening the first one found.  Then set any global variables defined in the config file.
THIS FUNCTION WILL CAUSE THE PROGRAM TO TERMINATE IF THE CONFIG FILE CANNOT BE OPENED
*/
short read_config_file(void);

/*
Function: short validate_config_file(void)
Inputs: None
Outputs: Ensures that the config file didn't give contradictory parameters.  Returns pass if there are no conflicts.
THIS FUNCTION WILL CAUSE THE PROGRAM TO TERMINATE IF ANY CONFLICTS ARE DETECTED
*/
short validate_config_file(void);

/*
Function: short read_query(char *fn)
Inputs: The name of the query file
Outputs: Reads in the query file.  Returns pass if successful.
THIS FUNCTION WILL CAUSE THE PROGRAM TO TERMINATE IF ANY PROBLEMS OCCUR WHILE READING THE QUERY
*/
short read_query(char *fn);

/*
Function: short read_fspd_sb(char *fn)
Inputs: The name of the query file
Outputs: Reads the query file and fragments the queries for fspd 
THIS FUNCTION WILL CAUSE THE PROGRAM TO TERMINATE IF ANY PROBLEMS OCCUR WHILE READING THE QUERY
*/
short read_fspd_sb(char *fn);

/*
Function: short tokenize_query_info(void)
Inputs: None
Outputs: Splits each query line into its various parameters and stores each query in the query array.  Returns pass if successful.
THIS FUNCTION WILL CAUSE THE PROGRAM TO TERMINATE IF ANY PROBLEMS OCCUR WHILE TOKENIZING THE QUERIES
*/
short tokenize_query_info(void);

/*
Function: short read_dbl(char *fn)
Inputs: The name of the file containing a list of the database files.
Outputs: Counts the number of database files, then stores all the filenames in the array fn_db.  Also calls read_dbl_test, which seems to unconditionally pass?  Then returns pass if successful.
THIS FUNCTION WILL CAUSE THE PROGRAM TO TERMINATE IF ANY PROBLEMS OCCUR WHILE READING THE LIST OF DATABASES
*/
short read_dbl(char *fn);

/*
Function: short read_dbl(char *fn)
Inputs: The name of the file containing a list of the database files.
Outputs: Prints all the database files, the returns pass.
*/
short read_dbl_test(char *fn);
