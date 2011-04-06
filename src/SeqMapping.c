/*
 * this module will determine the energy components given the alignments.
 *
 */

#include "lib.h"
#include "defs.h"
#include "SeqMapping.h"


/* determine each energy components given the alignments */

int SeqMapping (char * tseq, char * qseq, int * seqmapping,  
        int * inversemapping, int tlength, int qlength, int * identical)
  { 
		register	int	ndx, ndx1, ndx2, length1, length2;
		int	minscore, min, sim;
		int	** scorematrix, ** insertmatrix, ** deletematrix;
		minscore = LargeScore;
	
/***********************************************************************/
	/* three type of matrix required for dynamic programming */
	/* scorematrix = Matrix[ query length ]x[ homologue length ] */
	/* insertmatrix = Matrix[ query length ]x[ homologue length ] */
	/* deletematrix = Matrix[ query length ]x[ homologue length ] */
/***********************************************************************/

		scorematrix  = (int **) malloc ((tlength+1) * sizeof (int*));
		insertmatrix = (int **) malloc ((tlength+1) * sizeof (int*));
		deletematrix = (int **) malloc ((tlength+1) * sizeof (int*));
		for (ndx1 = 0; ndx1 <= tlength; ndx1++) 
		{  
			scorematrix[ndx1]  = (int *) malloc ((qlength+1) * sizeof (int));
			insertmatrix[ndx1] = (int *) malloc ((qlength+1) * sizeof (int));
			deletematrix[ndx1] = (int *) malloc ((qlength+1) * sizeof (int));
		}

		/* dynamic programming initialization */      
		/* do a global alignment */
		scorematrix[0][0]  = 0.0;    
		insertmatrix[0][0] = 0.0; 
		deletematrix[0][0] = 0.0;   
		for (ndx1 = 1; ndx1 <= tlength; ndx1++)					
		//matrix[i][0] = initializing
		{
			scorematrix[ndx1][0]  = OpenGapPenalty + ndx1 * OneIndelPenalty;
			insertmatrix[ndx1][0] = OpenGapPenalty + ndx1 * OneIndelPenalty;
			deletematrix[ndx1][0] = LargeScore;
		}   
		for (ndx1 = 1; ndx1 <= qlength; ndx1++)					
		//matrix[0][i] = initializing		
		{ 
			scorematrix[0][ndx1]  = OpenGapPenalty + ndx1 * OneIndelPenalty; 
			deletematrix[0][ndx1] = OpenGapPenalty + ndx1 * OneIndelPenalty;
			insertmatrix[0][ndx1] = LargeScore;
		}   
		/* dynamic programming recurrence */
		for (ndx1 = 1; ndx1 <= tlength; ndx1++)   
		{ 
			for (ndx2 = 1; ndx2 <= qlength; ndx2++)   
			{ 
				deletematrix[ndx1][ndx2] = MIN(deletematrix[ndx1-1][ndx2] 
					+ OneIndelPenalty, MIN (scorematrix[ndx1-1][ndx2] 
					+ OpenGapPenalty + OneIndelPenalty, insertmatrix[ndx1-1]
					[ndx2] + OpenGapPenalty + OneIndelPenalty));

				insertmatrix[ndx1][ndx2] = MIN(insertmatrix[ndx1][ndx2-1] 
					+ OneIndelPenalty, MIN (scorematrix[ndx1][ndx2-1]	
					+ OpenGapPenalty + OneIndelPenalty, deletematrix[ndx1]
					[ndx2-1] + OpenGapPenalty + OneIndelPenalty));

				scorematrix[ndx1][ndx2] = MIN (scorematrix[ndx1-1][ndx2-1], 					MIN (deletematrix[ndx1-1][ndx2-1], insertmatrix[ndx1-1]
					[ndx2-1]))+ GetFitnessScore(tseq, ndx1-1, qseq, ndx2-1);
      }
    }


   minscore = MIN (deletematrix[tlength][qlength],
              MIN (insertmatrix[tlength][qlength], scorematrix[tlength]
              [qlength]));

    /* retrive alignments  */                                             
	if (minscore == deletematrix[tlength][qlength])
		ExtractAlignment (tseq, qseq, scorematrix, insertmatrix,  
		deletematrix, tlength, qlength, -1, 0, seqmapping, inversemapping, 1);
	else if (minscore == insertmatrix[tlength][qlength])
		ExtractAlignment (tseq, qseq, scorematrix, insertmatrix,  
		deletematrix, tlength, qlength, -1, 0, seqmapping, inversemapping, 2);
	else
		ExtractAlignment (tseq, qseq, scorematrix, insertmatrix,  
		deletematrix, tlength, qlength, -1, 0, seqmapping, inversemapping, 0);

   sim = 0;
	for (ndx = 0; ndx < qlength; ndx++)
	{
		if (seqmapping[ndx] != INVALID && qseq[ndx] == tseq[seqmapping[ndx]])
		sim++;
	}

	*identical = sim;
	
	/* free the matrices */
	for (ndx1 = 0; ndx1 <= tlength; ndx1++) 
	{  
		free (scorematrix[ndx1]);
		free (insertmatrix[ndx1]);
		free (deletematrix[ndx1]);
	}
	free (scorematrix);
	free (insertmatrix);
	free (deletematrix);

	return minscore;
}

/* Extract the alignment from the dynamically filler in matrices */
void ExtractAlignment (char * tseq, char * qseq, int ** scorematrix, 
  int ** insertmatrix, int ** deletematrix, int ndx1, int ndx2, 
  int temp_start, int seq_start, int * seqmapping, int * inversemapping,
  int flag) 
  {  
    if (ndx1 == 0 || ndx2 == 0) return;

    /* use score matrix */
    if (flag == 0)
    { 
      seqmapping[seq_start+ndx2-1] = temp_start+ndx1; 
      inversemapping[temp_start+ndx1] = seq_start+ndx2-1;

      if (scorematrix[ndx1][ndx2] - scorematrix[ndx1-1][ndx2-1] ==
             GetFitnessScore(tseq, temp_start+ndx1, qseq, seq_start+ndx2-1))
      { 
        ExtractAlignment (tseq, qseq, scorematrix, insertmatrix, 
          deletematrix, ndx1-1, ndx2-1, temp_start, seq_start, 
          seqmapping, inversemapping, 0);
      } 
      else if (scorematrix[ndx1][ndx2] - deletematrix[ndx1-1][ndx2-1] ==
             GetFitnessScore(tseq, temp_start+ndx1, qseq, seq_start+ndx2-1))
      { 
        ExtractAlignment (tseq, qseq, scorematrix, insertmatrix, 
           deletematrix, ndx1-1, ndx2-1, temp_start, seq_start, 
           seqmapping, inversemapping, 1);
      }
      else if (scorematrix[ndx1][ndx2] - insertmatrix[ndx1-1][ndx2-1] ==
             GetFitnessScore(tseq, temp_start+ndx1, qseq, seq_start+ndx2-1))
      { 
        ExtractAlignment (tseq, qseq, scorematrix, insertmatrix, 
           deletematrix, ndx1-1, ndx2-1, temp_start, seq_start, 
           seqmapping, inversemapping, 2);
      }
      else 
      { 
        fprintf(stderr, "Unable to retrive alignment 1: ndx1=%d, ndx2=%d,\
         \ntemp_start=%d, seq_start=%d. Quit!\n",
           ndx1, ndx2, temp_start, seq_start);
        exit (10);
      }
    } 
    /* use delete matrix */
    else if (flag == 1)
    { 
      if (deletematrix[ndx1][ndx2] - deletematrix[ndx1-1][ndx2] ==
            OneIndelPenalty)
      { 
        ExtractAlignment (tseq, qseq, scorematrix, insertmatrix, 
           deletematrix, ndx1-1, ndx2, temp_start, seq_start, 
           seqmapping, inversemapping, 1);
      }
      else if (deletematrix[ndx1][ndx2] - insertmatrix[ndx1-1][ndx2] ==
             OpenGapPenalty + OneIndelPenalty)
      { 
        ExtractAlignment (tseq, qseq, scorematrix, insertmatrix, 
           deletematrix, ndx1-1, ndx2, temp_start, seq_start, 
           seqmapping, inversemapping, 2);
      }
      else if (deletematrix[ndx1][ndx2] - scorematrix[ndx1-1][ndx2] ==
             OpenGapPenalty + OneIndelPenalty) 
      { 
        ExtractAlignment (tseq, qseq, scorematrix, insertmatrix, 
          deletematrix, ndx1-1, ndx2, temp_start, seq_start, 
          seqmapping, inversemapping, 0);
      } 
      else 
      { 
        fprintf(stderr, "Unable to retrive alignment 1: ndx1=%d, ndx2=%d,\
         \ntemp_start=%d, seq_start=%d. Quit!\n",
           ndx1, ndx2, temp_start, seq_start);
        exit (11);
      }
    } 
    /* use insert matrix */
    else if (flag == 2)
    { 
      if (insertmatrix[ndx1][ndx2] - insertmatrix[ndx1][ndx2-1] ==
               OneIndelPenalty) 
      { 
        ExtractAlignment (tseq, qseq, scorematrix, insertmatrix, 
           deletematrix, ndx1, ndx2-1, temp_start, seq_start, 
           seqmapping, inversemapping, 2);
      }
      else if (insertmatrix[ndx1][ndx2] - deletematrix[ndx1][ndx2-1] ==
            OpenGapPenalty + OneIndelPenalty) 
      { 
        ExtractAlignment (tseq, qseq, scorematrix, insertmatrix, 
           deletematrix, ndx1, ndx2-1, temp_start, seq_start, 
           seqmapping, inversemapping, 1);
      }
      else if (insertmatrix[ndx1][ndx2] - scorematrix[ndx1][ndx2-1] ==
             OpenGapPenalty + OneIndelPenalty) 
      { 
        ExtractAlignment (tseq, qseq, scorematrix, insertmatrix, 
          deletematrix, ndx1, ndx2-1, temp_start, seq_start, 
          seqmapping, inversemapping, 0);
      } 
      else 
      { 
        fprintf(stderr, "Unable to retrive alignment 1: ndx1=%d, ndx2=%d,\
         \ntemp_start=%d, seq_start=%d. Quit!\n",
           ndx1, ndx2, temp_start, seq_start);
        exit (12);
      }
    } 
    else
    { 
      fprintf(stderr, "Wrong flag (%d): ndx1=%d, ndx2=%d Quit!\n",
        flag, ndx1, ndx2);
      exit (20);
    } 
    return;
}

/* construct a 2D similarity matrix for tseq onto qseq */
void ConstructSim (char * tseq, char * qseq, int * seqmapping, int * inversemapping, int length, short ** similarity)
  {
     register int ndx, ndx1, ndx2, ndx0, s_type, t_type;
     int base_shift = 0, shift = 0, sim, size;
		
		/* let shift be the last valid index in seqmapping */
     for (ndx = length - 1; ndx > 0; ndx--)
     {
       if (seqmapping[ndx] != INVALID)
       { shift = ndx; break; }
     }

		
     for (ndx0 = 0; ndx0 < length - PSPD_MINSEGLENGTH; ndx0++)
     { 
       sim = 0;
       for (ndx = ndx0; ndx < length; ndx++)
       { 
         ndx1 = ndx + base_shift;

         if (seqmapping[ndx] != INVALID && 
             qseq[ndx] == tseq[seqmapping[ndx]])
           sim++;

         ndx2 = seqmapping[ndx]+1;
			
			
         while (ndx < shift && inversemapping[ndx2] == INVALID)
         {
           ndx1++;
           ndx2++;
           base_shift++;
         }
         size = ndx1 - ndx0 + 1;
         /* update the similarity only if the subsequence is big enough */
         /* and not too similar */
         if(ndx-ndx0+1 >= PSPD_MINSEGLENGTH && (float)sim/(float)size<PSPD_MAX_SIMILARITY)
           similarity[ndx0][ndx]++;
// printf ("ndx0=%d, ndx1=%d, size=%d, sim=%d\n", ndx0, ndx1, size, sim);
       }
     }
  }
