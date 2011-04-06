/*****************************************************************************
This file is part of the Primegens Project.  

PUT COPYRIGHT INFO HERE

PUT VERSION INFO HERE
*****************************************************************************/


/*
 * primer design module, used for primer design using various different 
 * programs.
 *
 */
 
/*
Function: short find_maxCutSite_region(Query qry, int *start, int *length)
Inputs: a query, and integers by reference to hold the start and length of the max cut site
Outputs: the values pointed to by start and length have been updated appropriately
THIS FUNCTION MAY CAUSE THE PROGRAM TO TERMINATE IF THE XLS FILE CANNOT BE OPENED
*/
short find_maxCutSite_region(Query qry, int *start, int *length);

/*
Function: short run_primer3(Query qry)
Inputs: a query on which to run primer3 
Outputs: The primer3 output is stored in the file specified by fnp3_out.  Returns PASS if successful.
THIS FUNCTION MAY CAUSE THE PROGRAM TO TERMINATE IF THE PRIMER3 INPUT FILE CANNOT BE OPENED
*/
short run_primer3(Query qry);

/*
Function: short save_primer3_output(void)
Inputs: None
Outputs: Reads in the primer pairs from the file named by fnp3_out
THIS FUNCTION MAY CAUSE THE PROGRAM TO TERMINATE IF THE PRIMER3 OUTPUT FILE CANNOT BE OPENED
*/
short save_primer3_output(void);
