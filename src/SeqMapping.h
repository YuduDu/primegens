/*****************************************************************************
This file is part of the Primegens Project.  

PUT COPYRIGHT INFO HERE

PUT VERSION INFO HERE
*****************************************************************************/


/*
 * this module will determine the energy components given the alignments,
 * extract an alignment, and construct a similarity matrix.
 */

/*macros definition*/
#define MAX(a,b)  ((a) > (b) ? (a):(b))
#define MIN(a,b)  ((a) < (b) ? (a):(b))
																						
/*
Function:int SeqMapping (char * tseq, char * qseq, int * seqmapping,  
       	int * inversemapping, int tlength, int qlength, int * identical)
Inputs: 	CHAR *TSEQ - any of other homologue sequence string for pairwise global alignment 
			CHAR *QSEQ - the query sequence string, which is the 'star' in global alignment 
			INT *SEQMAPPING - buffer[query length] array for query sequence of its sequence size 
			INT *INVERSEMAPPING - same type buffer(int array) for the other sequence of its size 
			INT TLENGTH - length of other homologue sequence string	
			INT QLENGTH - length of the query sequence string			  
			INT *IDENTICAL - the number of matches in the alignment by reference
Outputs: The tmp score, and IDENTICAL is updated by reference
*/
int SeqMapping (char * tseq, char * qseq, int * seqmapping,  
        int * inversemapping, int tlength, int qlength, int * identical);

/*
Function:void ExtractAlignment (char * tseq, char * qseq, int **
			scorematrix, int ** insertmatrix, int ** deletematrix, int ndx1,
			int ndx2, int temp_start, int seq_start, int * seqmapping, 
			int * inversemapping, int flag)
Inputs: 	CHAR * TSEQ - the sequence who alignment is to be extracted
			CHAR * QSEQ - the query sequence string
			INT ** SCOREMATRIX - the score matrix, filled in dynamically above
			INT ** INSERTMATRIX - the insert matrix, filled in dynamically above
			INT ** DELETEMATRIX - the delete matrix, filled in dynamically above
			INT NDX1 - the index in the temp sequence
			INT NDX2 - the index in the query sequence
			INT TEMP_START - the start position of the temp sequence
			INT SEQ_START - the start position of the sequence
			INT * SEQMAPPING - the sequence mapping extracted from the matrices (an array of indices)
			INT * INVERSEMAPPING - the inverse of seqmapping
			INT FLAG - A flag determining which matrix to extract from
Outputs: seqmapping and inversemapping are filled in using the value from the appropriate matrix
*/
void ExtractAlignment (char * tseq, char * qseq, int ** scorematrix, 
  int ** insertmatrix, int ** deletematrix, int ndx1, int ndx2, 
  int temp_start, int seq_start, int * seqmapping, int * inversemapping,
  int flag);
  
/*
Function:void ConstructSim (char * tseq, char * qseq, int * seqmapping,
			int * inversemapping, int length, short ** similarity)
Inputs: 	CHAR * TSEQ - the sequence whose similarity is being extracted
			CHAR * QSEQ - the query sequence
			INT * SEQMAPPING - the sequence mapping of tseq onto qseq (an array of indices)
			INT * INVERSEMAPPING - the inverse of seqmapping
			INT LENGTH - the length of the query sequence
			SHORT ** SIMILARITY - a 2D array of similarity values
Outputs: the similarity matrix is filled in such that similarity[x][y] is 1 if the alignment from tseq onto qseq starting going from indices x through y is not more similar than the maximum similarity and is longer than the minimum segment length, and is 0 otherwise
*/
void ConstructSim (char * tseq, char * qseq, int * seqmapping, int * inversemapping, int length, short ** similarity);
