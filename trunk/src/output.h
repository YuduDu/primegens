/*****************************************************************************
This file is part of the Primegens Project.  

PUT COPYRIGHT INFO HERE

PUT VERSION INFO HERE
*****************************************************************************/


/*
 * this module contains all functions relating to output
 * of the final results; many other functions have optional
 * debugging outputs to standard output, however this module
 * is for outputing the final results only
 */

/*
Function: int rank_ppair(void)
Inputs: None
Outputs: Sorts the primer pairs by the number of hybridizations, returns PASS if successful
*/
int rank_ppair(void);

/*
Function: int sort_on_hybrid(const void *ptr1, const void *ptr2)
Inputs: Two pointers to primer pair structures
Outputs: Returns 1 if ptr1 has more hybridizations, -1 if ptr2 has more hybridizations, and 0 if equal.  Used to sort primer pairs.
*/
int sort_on_hybrid(const void *ptr1, const void *ptr2);

/*
Function: int print_output(int query_ndx)
Inputs: The index of the query to be printed to the files
Outputs: The output associated with the query is printed to the appropriate output files.
THIS FUNCTION WILL CAUSE THE PROGRAM TO TERMINATE IF THE OUTPUT FILES CANNOT BE OPENED
*/
int print_output(int query_ndx);
