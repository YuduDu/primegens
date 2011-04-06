
/*this module will generate the blast hit record and potential 
 * amplicons for each primer pair for any query sequence.
 */


#include "lib.h"
#include "defs.h"
#include "pphit.h"
#include "util.h"

/* Store the blast hits for primer pairs from chromosomes */
int ppair_chit_record_chromosome(void) {
	int pp=0, uo=0, br=0;
	int ppc=0, pmc=0;
	size_t t_size;
	char chr[BISULCHR];
	int bh=0, h=0, bht=0, bhs=0, bhe=0, c=0;
	BRecord  *tbrc1, *tbrc2;

	printf("%s chit records for all primer pairs...\n",T_S());
	t_size = BLASTRECORDS*sizeof(BRecord);
	tbrc1 = (BRecord *)p_malloc(t_size);
	tbrc2 = (BRecord *)p_malloc(t_size);

	/* for each primer pair */
	for(pp=0; pp < ppair_count; pp++) {
		/*find left primer blast records*/
		if (dlevel > 3)
			printf("%s find and copy left primer blast record.\n",T_S());
		for(uo=0; uo < uoligo_count; uo++) {
			if (dlevel > 3){
				printf("uoligo[%d].nseq = %s\t", uo, uoligo[uo].nseq);
				printf("ppair[%d].lseq = %s\n", pp, ppair[pp].lseq);
			}
			if(p_oligo_strcmp(uoligo[uo].nseq, ppair[pp].lseq)==0) {
				for(br=0; br < BLASTRECORDS; br++) {
					strcpy(chr, uoligo[uo].brecord[br].chr);
					strcpy(tbrc1[br].chr, chr);
					if (dlevel > 3)
						printf("%s cid %2d [%s]\tcopied\n",T_S(),br,chr);
					tbrc1[br].bh = uoligo[uo].brecord[br].bh;
					/* copy the hits into tbrc1 */
					if(uoligo[uo].brecord[br].bh > 0) {
						bh = uoligo[uo].brecord[br].bh;
						tbrc1[br].bht=(int *)p_malloc(bh*sizeof(int));
						tbrc1[br].bhs=(int *)p_malloc(bh*sizeof(int));
						tbrc1[br].bhe=(int *)p_malloc(bh*sizeof(int));
						for(h=0; h < bh; h++) {
							bht = uoligo[uo].brecord[br].bht[h];
							bhs = uoligo[uo].brecord[br].bhs[h];
							bhe = uoligo[uo].brecord[br].bhe[h];
							tbrc1[br].bht[h]= bht;
							tbrc1[br].bhs[h]= bhs;
							tbrc1[br].bhe[h]= bhe;
						}
					}
				}
				break;
			}
		}
		if(uo==uoligo_count)
			p_exit("left primer not in uoligo list?\n");
		
		/*find right primer blast records*/
		if (dlevel > 3)
			printf("%s find and copy right primer blast record.\n",T_S());
		for(uo=0; uo < uoligo_count; uo++) {
			if(p_oligo_strcmp(uoligo[uo].nseq, ppair[pp].rseq)==0) {
				for(br=0; br < BLASTRECORDS; br++) {
					strcpy(chr, uoligo[uo].brecord[br].chr);
					strcpy(tbrc2[br].chr, chr);
					if(dlevel > 3)
						printf("%s cid %2d [%s]\tcopied\n",T_S(),br,chr);
					tbrc2[br].bh = uoligo[uo].brecord[br].bh;
					/* copy the hits into tbrc2 */
					if(uoligo[uo].brecord[br].bh > 0) {
						bh = uoligo[uo].brecord[br].bh;
						tbrc2[br].bht=(int *)p_malloc(bh*sizeof(int));
						tbrc2[br].bhs=(int *)p_malloc(bh*sizeof(int));
						tbrc2[br].bhe=(int *)p_malloc(bh*sizeof(int));
						for(h=0; h < bh; h++) {
							bht = uoligo[uo].brecord[br].bht[h];
							bhs = uoligo[uo].brecord[br].bhs[h];
							bhe = uoligo[uo].brecord[br].bhe[h];
							tbrc2[br].bht[h]= bht;
							tbrc2[br].bhs[h]= bhs;
							tbrc2[br].bhe[h]= bhe;
						}
					}
				}
				break;
			}
		}
		if(uo==uoligo_count)
			p_exit("right primer not in uoligo list?\n");
		
		/*
		printf("%s description of tbrc1 and tbrc2 records...\n",T_S());
		for(c=0; c < BLASTRECORDS; c++) {
		}*/

		/* make sure tbrc1 matchs tbrc2 */
		for(c=0; c < BLASTRECORDS; c++) {
			if(dlevel >3){
				printf("%s checking brecord # %d\n",T_S(),c);
				printf("tbrc1[%d].chr = %s <-", c, tbrc1[c].chr);
				printf("-> tbrc2[%d].chr = %s\n", c, tbrc2[c].chr);
			}
			if(strcmp(tbrc1[c].chr, tbrc2[c].chr)!=0) {
				p_exit("tbrc1.chr and tbrc2.chr are not same!\n");
			}
		}
			
		if(dlevel > 3)
			printf("%s allocate memory for pair chit record.\n",T_S());
		t_size = BLASTRECORDS*sizeof(CHit);
		ppair[pp].chit = (CHit *)p_malloc(t_size);
		
		/* count the number of PP hits and PM hits */
		for(c=0; c < BLASTRECORDS; c++) {
			strcpy(ppair[pp].chit[c].cid, tbrc1[c].chr); /*=tbrc2.chr*/
			ppair[pp].chit[c].ppc = 0;
			ppair[pp].chit[c].pmc = 0;
			for(h=0; h < tbrc1[c].bh; h++) {
				if(tbrc1[c].bht[h]==PP)
					ppair[pp].chit[c].ppc ++;
				else if(tbrc1[c].bht[h]==PM)
					ppair[pp].chit[c].pmc ++;
				else
					p_exit("all hits should be either PP or PM!\n");
			}
			for(h=0; h < tbrc2[c].bh; h++) {
				if(tbrc2[c].bht[h]==PP)
					ppair[pp].chit[c].ppc ++;
				else if(tbrc2[c].bht[h]==PM)
					ppair[pp].chit[c].pmc ++;
				else
					p_exit("all hits should be either PP or PM!\n");
			}

			/* allocate space for the pp hits */
			if(ppair[pp].chit[c].ppc > 0) {
				t_size = (ppair[pp].chit[c].ppc)*sizeof(int);
				ppair[pp].chit[c].pps=(int *)p_malloc(t_size);
				ppair[pp].chit[c].ppe=(int *)p_malloc(t_size);
			}
			
			/* allocate space for the pm hits */
			if(ppair[pp].chit[c].pmc > 0) {
				t_size = (ppair[pp].chit[c].pmc)*sizeof(int);
				ppair[pp].chit[c].pme=(int *)p_malloc(t_size);
				ppair[pp].chit[c].pms=(int *)p_malloc(t_size);
			}
			
			ppc = 0;
			pmc = 0;

			/* store the tbrc1 hits */
			for(h=0; h < tbrc1[c].bh; h++) {
				if(tbrc1[c].bht[h]==PP){
                                      ppair[pp].chit[c].pps[ppc]=tbrc1[c].bhs[h];
                                      ppair[pp].chit[c].ppe[ppc]=tbrc1[c].bhe[h];
                                      ppc++;}
				else if(tbrc1[c].bht[h]==PM){
                                      ppair[pp].chit[c].pme[pmc]=tbrc1[c].bhe[h];
                                      ppair[pp].chit[c].pms[pmc]=tbrc1[c].bhs[h];
                                      pmc++;}

				else
					p_exit("all hits should be either PP or PM!\n");
			}
			/* store the tbrc2 hits */
			for(h=0; h < tbrc2[c].bh; h++) {
				if(tbrc2[c].bht[h]==PP){
					ppair[pp].chit[c].pps[ppc]=tbrc2[c].bhs[h];
					ppair[pp].chit[c].ppe[ppc]=tbrc2[c].bhe[h];
					ppc++;
                                }
				else if(tbrc2[c].bht[h]==PM){
					ppair[pp].chit[c].pme[pmc]=tbrc2[c].bhe[h];
					ppair[pp].chit[c].pms[pmc]=tbrc2[c].bhs[h];
					pmc++;
                                }
				else
					p_exit("all hits should be either PP or PM!\n");
			}
			if(ppair[pp].chit[c].ppc != ppc)
				p_exit("some problem during pp grouping.\n");
			if(ppair[pp].chit[c].pmc != pmc)
				p_exit("some problem during PM grouping.\n");
		}
		if(dlevel > 3)
			printf("%s free tmp memory from inside brecord\n",T_S());
		/* free the temporary blast records */
		for(br=0; br < BLASTRECORDS; br++) {
			if(tbrc1[br].bh > 0) {
				free(tbrc1[br].bht);
				free(tbrc1[br].bhs);
				free(tbrc1[br].bhe);
			}
			if(tbrc2[br].bh > 0) {
				free(tbrc2[br].bht);
				free(tbrc2[br].bhs);
				free(tbrc2[br].bhe);
			}
		}
		printf("%s CHit record for ppair %02d...done!\n",T_S(),pp+1);
	}
	free(tbrc1);
	free(tbrc2);
	printf("%s chit records for all primer pairs...done!\n",T_S());
	return PASS;
}

/* store blast records for primer pairs from a single database */
int ppair_chit_record_single(void) {
	int pp=0, uo=0, br=0, flag=0;
	int ppc=0, pmc=0, chit_count=0;
	size_t t_size;
	char chr[SEQUENCEID];
	SeqID *tseqid;	/*buffer to get unique sequence ID for chit record*/
	int tseqid_count;
	int bh=0, h=0, bht=0, bhs=0, bhe=0, c=0, c1=0;
	int tbrc1_count, tbrc2_count;
	BRecord  *tbrc1, *tbrc2;

	printf("%s chit records for all primer pairs...\n",T_S());
	/* for each primer pair */
	for(pp=0; pp < ppair_count; pp++) {
		/*find left primer blast records*/
		for(uo=0; uo < uoligo_count; uo++) {
			if(p_oligo_strcmp(uoligo[uo].nseq, ppair[pp].lseq)==0) {
				/* copy the hits into tbrc1 */
				if(uoligo[uo].br_count==0) {
					tbrc1_count = 0;
				} else {
					tbrc1_count = uoligo[uo].br_count;
					t_size = uoligo[uo].br_count*sizeof(BRecord);
					tbrc1 = (BRecord *)p_malloc(t_size);
					for(br=0; br < uoligo[uo].br_count; br++) {
						strcpy(chr, uoligo[uo].brecord[br].chr);
						strcpy(tbrc1[br].chr, chr);
						tbrc1[br].bh = uoligo[uo].brecord[br].bh;
						if(uoligo[uo].brecord[br].bh > 0) {
							bh = uoligo[uo].brecord[br].bh;
							tbrc1[br].bht=(int *)p_malloc(bh*sizeof(int));
							tbrc1[br].bhs=(int *)p_malloc(bh*sizeof(int));
							tbrc1[br].bhe=(int *)p_malloc(bh*sizeof(int));
							for(h=0; h < bh; h++) {
								bht = uoligo[uo].brecord[br].bht[h];
								bhs = uoligo[uo].brecord[br].bhs[h];
								bhe = uoligo[uo].brecord[br].bhe[h];
								tbrc1[br].bht[h]= bht;
								tbrc1[br].bhs[h]= bhs;
								tbrc1[br].bhe[h]= bhe;
								if (dlevel > 3)
									printf("left-primer: %d %d %d\n",tbrc1[br].bht[h],tbrc1[br].bhs[h],tbrc1[br].bhe[h]);
							}
						}
					}
				}
				break;
			}
		}
		/*find right primer blast records*/
		for(uo=0; uo < uoligo_count; uo++) {
			if(p_oligo_strcmp(uoligo[uo].nseq, ppair[pp].rseq)==0) {
				/* copy the hits into tbrc2 */
				if(uoligo[uo].br_count==0) {
					tbrc2_count = 0;
				} else {
					tbrc2_count = uoligo[uo].br_count;
					t_size = uoligo[uo].br_count*sizeof(BRecord);
					tbrc2 = (BRecord *)p_malloc(t_size);
					for(br=0; br < uoligo[uo].br_count; br++) {
						strcpy(chr, uoligo[uo].brecord[br].chr);
						strcpy(tbrc2[br].chr, chr);
						tbrc2[br].bh = uoligo[uo].brecord[br].bh;
						if(uoligo[uo].brecord[br].bh > 0) {
							bh = uoligo[uo].brecord[br].bh;
							tbrc2[br].bht=(int *)p_malloc(bh*sizeof(int));
							tbrc2[br].bhs=(int *)p_malloc(bh*sizeof(int));
							tbrc2[br].bhe=(int *)p_malloc(bh*sizeof(int));
							for(h=0; h < bh; h++) {
								bht = uoligo[uo].brecord[br].bht[h];
								bhs = uoligo[uo].brecord[br].bhs[h];
								bhe = uoligo[uo].brecord[br].bhe[h];
								tbrc2[br].bht[h]= bht;
								tbrc2[br].bhs[h]= bhs;
								tbrc2[br].bhe[h]= bhe;
								if (dlevel > 3)
									printf("right-primer: %d %d %d\n",tbrc2[br].bht[h],tbrc2[br].bhs[h],tbrc2[br].bhe[h]);
							}
						}
					}
				}
				break;
			}
		}
		
		if((tbrc1_count + tbrc2_count)==0) {
			ppair[pp].chit_count = 0;
			continue;
		}
		/* create a lits of sequence IDs */
		t_size = (tbrc1_count + tbrc2_count)*sizeof(SeqID);
		tseqid = (SeqID *)p_malloc(t_size);
		tseqid_count = 0;
		/* if the sequence ID isn't in the list of SeqIDs, add it */
		for(br=0; br < tbrc1_count; br++) {
			flag=chk_if_seqid_exist(tbrc1[br].chr,tseqid,tseqid_count);
			if(flag == INVALID) {
				strcpy(tseqid[tseqid_count].id, tbrc1[br].chr);
				tseqid_count++;
			}
		}
		/* if the sequence ID isn't in the list of SeqIDs, add it */
		for(br=0; br < tbrc2_count; br++) {
			flag=chk_if_seqid_exist(tbrc2[br].chr,tseqid,tseqid_count);
			if(flag == INVALID) {
				strcpy(tseqid[tseqid_count].id, tbrc2[br].chr);
				tseqid_count++;
			}
		}
		ppair[pp].chit_count = tseqid_count;
		t_size = ppair[pp].chit_count*sizeof(CHit);
		ppair[pp].chit = (CHit *)p_malloc(t_size);
		for(c=0; c < ppair[pp].chit_count; c++) {
			strcpy(ppair[pp].chit[c].cid, tseqid[c].id);
			ppair[pp].chit[c].ppc = 0;
			ppair[pp].chit[c].pmc = 0;
		}
		free(tseqid);
		/* count the PP hits and PM hits */
		for(c=0; c < ppair[pp].chit_count; c++) {
			for(br=0; br < tbrc1_count; br++) {
				if(strcmp(ppair[pp].chit[c].cid,tbrc1[br].chr)==0) {
					for(h=0; h < tbrc1[br].bh; h++) {
						if(tbrc1[br].bht[h]==PP)
							ppair[pp].chit[c].ppc ++;
						else if(tbrc1[br].bht[h]==PM)
							ppair[pp].chit[c].pmc ++;
					}
				}
			}
			/* count the PP hits and PM hits */
			for(br=0; br < tbrc2_count; br++) {
				if(strcmp(ppair[pp].chit[c].cid,tbrc2[br].chr)==0) {
					for(h=0; h < tbrc2[br].bh; h++) {
						if(tbrc2[br].bht[h]==PP)
							ppair[pp].chit[c].ppc ++;
						else if(tbrc2[br].bht[h]==PM)
							ppair[pp].chit[c].pmc ++;
					}
				}
			}
			/* allocate space for the PP hits and PM hits */
			if(ppair[pp].chit[c].ppc > 0) {
				t_size = (ppair[pp].chit[c].ppc)*sizeof(int);
				ppair[pp].chit[c].pps=(int *)p_malloc(t_size);
				ppair[pp].chit[c].ppe=(int *)p_malloc(t_size);
			}
			if(ppair[pp].chit[c].pmc > 0) {
				t_size = (ppair[pp].chit[c].pmc)*sizeof(int);
				ppair[pp].chit[c].pme=(int *)p_malloc(t_size);
				ppair[pp].chit[c].pms=(int *)p_malloc(t_size);
			}	
			ppc = 0;
			pmc = 0;
			/* store the PP and PM hits */
			for(br=0; br < tbrc1_count; br++) {
				if(strcmp(ppair[pp].chit[c].cid,tbrc1[br].chr)==0) {
					for(h=0; h < tbrc1[br].bh; h++) {
						if(tbrc1[br].bht[h]==PP){
							ppair[pp].chit[c].pps[ppc]=tbrc1[br].bhs[h];
							ppair[pp].chit[c].ppe[ppc]=tbrc1[br].bhe[h];
                                                        ppc++;}
						else if(tbrc1[br].bht[h]==PM){
							ppair[pp].chit[c].pme[pmc]=tbrc1[br].bhe[h];
							ppair[pp].chit[c].pms[pmc]=tbrc1[br].bhs[h];
                                                        pmc++;}
					}
				}
			}
			/* store the PP and PM hits */
			for(br=0; br < tbrc2_count; br++) {
				if(strcmp(ppair[pp].chit[c].cid,tbrc2[br].chr)==0) {
					for(h=0; h < tbrc2[br].bh; h++) {
						if(tbrc2[br].bht[h]==PP){
							ppair[pp].chit[c].pps[ppc]=tbrc2[br].bhs[h];
							ppair[pp].chit[c].ppe[ppc]=tbrc2[br].bhe[h];
                                                        ppc++;
                                                        }
						else if(tbrc2[br].bht[h]==PM){
							ppair[pp].chit[c].pme[pmc]=tbrc2[br].bhe[h];
							ppair[pp].chit[c].pms[pmc]=tbrc2[br].bhs[h];
                                                        pmc++;
                                                        }
					}
				}
			}
			if(ppair[pp].chit[c].ppc != ppc)
				p_exit("some problem during pp grouping.\n");
			if(ppair[pp].chit[c].pmc != pmc)
				p_exit("some problem during PM grouping.\n");
		}
		/* free the temporary blast records */
		for(br=0; br < tbrc1_count; br++) {
			if(tbrc1[br].bh > 0) {
				free(tbrc1[br].bht);
				free(tbrc1[br].bhs);
				free(tbrc1[br].bhe);
			}
		}
		for(br=0; br < tbrc2_count; br++) {
			if(tbrc2[br].bh > 0) {
				free(tbrc2[br].bht);
				free(tbrc2[br].bhs);
				free(tbrc2[br].bhe);
			}
		}
		printf("%s CHit record for ppair %02d...done!\n",T_S(),pp+1);
		if(tbrc1_count > 0)
			free(tbrc1);
		if(tbrc2_count > 0)
			free(tbrc2);
	}
	printf("%s chit records for all primer pairs...done!\n",T_S());
	return PASS;
}

/* check if tid is in the SeqID list */
int chk_if_seqid_exist(char *tid, SeqID *seqid, int seqid_count) {
	int ndx=0;

	for(ndx=0; ndx < seqid_count; ndx++) {
		if(strcmp(seqid[ndx].id, tid)==0)
			return ndx;
	}
	return INVALID;
}

/* Not implemented */
int ppair_chit_record_multiple(void) {

	return PASS;
}


