/*****************************************************************************
This file is part of the Primegens Project.  

PUT COPYRIGHT INFO HERE

PUT VERSION INFO HERE
*****************************************************************************/


/*
 * this module contains all the basic utilities, which are
 * frequently used in the software for various purposes.
 */


/*
Function: int p_min (int x, int y)
Inputs: two integers, x and y
Outputs: returns the smaller integer 
*/
int p_min (int x, int y);


/*
Function: int p_max (int x, int y)
Input: Two integers, x and y
Output: Returns the larger integer
*/
int p_max (int x, int y);


/*
Function: int p_malloc (size_t x)
Input: A size_t which represents an integral number of bytes to be allocated
Output: A pointer to the allocated block of memory 
*/
void *p_malloc(size_t x);


/*
Function: char *T_S (void);
Input: None
Output: Returns a string representing the time
*/
char *T_S (void);


/*
Function: void p_exit(char *msg)
Input: A message to be printed before exiting the program.
Output: Writes msg to the log file, prints msg to the user, removes the tmp_dir directory, then terminates program execution.  
Will quit the program if unable to open log file
*/
void p_exit(char *msg);


/*
Function: void p_continue(char *msg)
Input: A message to be printed to the user
Output: Writes msg to the log file, then prints msg to the user 
THIS FUNCTION WILL CAUSE THE PROGRAM TO TERMINATE IF THE LOG FILE CANNOT BE OPENED
*/
void p_continue(char *msg);


/*
Function: FILE *p_fopen_error_exit(char *fn, char *m)
Input: A filename fn to be opened and a mode m.
Output: Returns a pointer to the opened file.  
THIS FUNCTION WILL CAUSE THE PROGRAM TO TERMINATE IF THE GIVEN FILE CANNOT BE OPENED
*/
FILE *p_fopen_error_exit(char *fn, char *m);


/*
Function: void p_exec(char *cmd);
Input: A string containing a command for the system to execute
Output: The given command is passed to the system; if the command returns a failure, then the program will exit.
THIS FUNCTION WILL CAUSE THE PROGRAM TO TERMINATE IF THE GIVEN COMMAND RETURNS A FAILURE
*/
void p_exec(char *cmd);


/*
Function: int p_oligo_strcmp(char *str1, char *str2)
Input: Two strings of oligos
Output: Compares the two strings starting at the end until it reaches the beginning of the shorter string.  Returns 0 if the two strings are the same, otherwise return non-zero
*/
int p_oligo_strcmp(char *str1, char *str2);


/*
Function: void p_oligo_strcpy(char *str1, char *str2)
Input: A destination oligo string and a source oligo string for copying
Output: Copies the str2 into str1
*/
void p_oligo_strcpy(char *str1, char *str2);


/*
Function: char p_clean(void)
Input: None
Output: Frees allocated memory, removes all temporary created files and directory.  Returns nothing, should be void?
*/
char p_clean(void);


/*
Function: short get_nseq_from_db(char *qid, char *nseq, int *len);
Input: The query id, a string to read the query sequence into, and an integer by reference to hold the length of the sqeuence.
Output: Finds the given query in the database, then copies it into nseq and updates the reference of len.  Returns PASS if the function competes successfully
THIS FUNCTION WILL CAUSE THE PROGRAM TO TERMINATE IF THE DATABASE FILE CANNOT BE OPENED
*/
short get_nseq_from_db(char *qid, char *nseq, int *len);


/*
Function: short get_hybd_nseq_from_db(char *qid, int ls, int le, int rs, int re, char *hlseq, char *hrseq, int *hllen, int *hrlen);
Input: The query id, integers representing the left sequence's start and end, integers representing the right sequence's start and end, strings to hold the left and right hybridization, as well as integers by reference to hold the length of the hybridized sequences.
Output: Searches the database for the given left and right sequences, stores them in the strings hlseq and hrseq, and updates the values referenced by hrlen and hllen.  Returns PASS if the function competes successfully
THIS FUNCTION WILL CAUSE THE PROGRAM TO TERMINATE IF THE DATABASE FILE CANNOT BE OPENED
*/
short get_hybd_nseq_from_db(char *qid, int ls, int le, int rs, int re, char *hlseq, char *hrseq, int *hllen, int *hrlen);


/*
Function: get_oligo_tm_for_hybd(char *hlseq,char *hrseq,char *lseq,char *rseq,int *potl)
Input: Strings representing the hybrid left and right sequences, as well as the left and right primer sequences, and an integer by reference to store the tm
Output: potl is set to 0 if there is a mismatch, otherwise 1.  Returns PASS if the function competes successfully
*/
short get_oligo_tm_for_hybd(char *hlseq,char *hrseq,char *lseq,char *rseq,int *potl);


/*
Function: short get_nseq_from_chromosome(int q)
Input: The query number q
Output: Reads the query's sequence into the query structure from the database.  Returns PASS if the function completes successfully
*/
short get_nseq_from_chromosome(int q);


/*
Function: int GetFitnessScore(char *tseq,int tposn,char *qseq, int qposn)
Input: Two strings of nucleotides, and indices of nucleotides to be compared
Output: Return -1 if the two nucleotides are equal, otherwise return 1
*/
int GetFitnessScore(char *tseq,int tposn,char *qseq, int qposn);


/*
Function: short read_query_toget_format(char *fn)
Input: The filename of the queries
Output: The global variable QUERY_FORMAT is set to the appropriate format.  Returns PASS if the function competes successfully
THIS FUNCTION WILL CAUSE THE PROGRAM TO TERMINATE IF THE DATABASE FILE CANNOT BE OPENED
*/
short read_query_toget_format(char *fn);
