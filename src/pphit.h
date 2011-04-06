/*****************************************************************************
This file is part of the Primegens Project.  

PUT COPYRIGHT INFO HERE

PUT VERSION INFO HERE
*****************************************************************************/


/*
 * this module will generate the blast hit record and potential 
 * amplicons for each primer pair for any query sequence.
 *
 */
 
/*
Function: int ppair_chit_record_chromosome(void)
Inputs: None
Outputs: Stores the blast hits for primer pairs from chromosomes
THIS FUNCTION MY CAUSE THE PROGRAM TO TERMINATE IF THERE IS A PROBLEM WITH THE BLAST HITS
*/
int ppair_chit_record_chromosome(void);

/*
Function: int ppair_chit_record_single(void)
Inputs: None
Outputs: Stores the blast hits for primer pairs from a single database
THIS FUNCTION MY CAUSE THE PROGRAM TO TERMINATE IF THERE IS A PROBLEM WITH THE BLAST HITS
*/
int ppair_chit_record_single(void);

/* not implemented */
int ppair_chit_record_multiple(void);

/*
Function: int chk_if_seqid_exist(char *tid, SeqID *seqid, int seqid_count)
Inputs: The temporary sequence ID, the list of sequence IDs, and the number of sequence IDs
Outputs: Returns the index of tid in SeqID, or INVALID if not in the list */
int chk_if_seqid_exist(char *tid, SeqID *seqid, int seqid_count);

