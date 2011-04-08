/*
 * this module contains all the functions for determining hybridizations.
 */


#include "lib.h"
#include "defs.h"
#include "oligotm.h"
#include "hbdzn.h"

/* check the hybridizations for a query against a single database */
int chk_hybridization_single(int q_ndx) {
	int pp=0;							/*primer pair index*/
	int c=0;								/*chromosome index*/
	int p=0;								/*plus-plus hit index*/
	int m=0;								/*plus-minus hit index*/
	int p_p=0, p_m=0, best=0;		/*temporary buffer*/
	int left=0, right=0;
	int ampl=0;							/*amplicon index*/
	size_t t_size;
	char cid[SEQUENCEID];
        int ret;
        int r1=0, r2=0;


	printf("%s find all hybridization for ppair...\n",T_S());
	/* for each primer pair */
	for(pp=0; pp < ppair_count; pp++) {
		ppair[pp].hbrdn_count = 0;
		/* skip if there are no BLAST hits */
		if(ppair[pp].chit_count == 0)
			continue;
		/* for each BLAST hit, count the potential hybridizations */
		for(c=0; c < ppair[pp].chit_count; c++) {
			if(ppair[pp].chit[c].ppc == 0)
				continue;				/*amplicons not possible*/
			if(ppair[pp].chit[c].pmc == 0)
				continue;				/*amplicons not possible*/
			/* for each plus-plus start */
			for(p=0; p < ppair[pp].chit[c].ppc; p++) {
				left = ppair[pp].chit[c].pps[p];
				best = right = ppair[pp].chit[c].pme[0];	
				/* for each plus-minus end */
				for(m=1; m < ppair[pp].chit[c].pmc; m++) {
					right = ppair[pp].chit[c].pme[m];
					if((left < right) && (right < best)){
						best = right;		
						r1=ppair[pp].chit[c].pms[m];
						r2=ppair[pp].chit[c].pme[m];
					}
				}
				right = best;
				if(((right - left) <= AMPSB)&&((right-left) >0)){
					if(dlevel >3) {
						printf("this could be potential hybridization.\n");
						printf("left:[%d,%d]\t",ppair[pp].chit[c].pps[p],ppair[pp].chit[c].ppe[p]);
						printf("right:[%d,%d]\n", r1, r2);
					}
					get_hybd_nseq_from_db(ppair[pp].chit[c].cid,ppair[pp].chit[c].pps[p],ppair[pp].chit[c].ppe[p],r1,r2,&ppair[pp].chit[c].hlseq,&ppair[pp].chit[c].hrseq,&ppair[pp].chit[c].hllen,&ppair[pp].chit[c].hrlen);	
					if(dlevel >3) { 
			      	printf("left-primer sequence: %s\nright-primer sequence: %s\n",ppair[pp].lseq,ppair[pp].rseq);
				   	printf("left-hit sequence: %s\nright-hit sequence: %s\n",ppair[pp].chit[c].hlseq,ppair[pp].chit[c].hrseq);
				   }
					get_oligo_tm_for_hybd(ppair[pp].chit[c].hlseq,ppair[pp].chit[c].hrseq,ppair[pp].lseq,ppair[pp].rseq,&ppair[pp].chit[c].potential);

					if(ppair[pp].chit[c].potential == 1){ 
				       ppair[pp].hbrdn_count ++;
					}
				}
			}
		}
		printf("%s total hybridizations for ppair %02d",T_S(),pp+1);
		printf(" [%s]\t[%s]",ppair[pp].lseq, ppair[pp].rseq);
		printf("...%d\n", ppair[pp].hbrdn_count);
		
		/* store the potential hybridizations from above */
		t_size = (ppair[pp].hbrdn_count)*sizeof(Hbrdn);
		ppair[pp].hbrdn=(Hbrdn *)p_malloc(t_size);
		ampl = 0;
		for(c=0; c < ppair[pp].chit_count; c++) {
			if(ppair[pp].chit[c].ppc == 0)
				continue;		/*amplicons not possible*/
			if(ppair[pp].chit[c].pmc == 0)
				continue;		/*amplicons not possible*/
			for(p=0; p < ppair[pp].chit[c].ppc; p++) {
				left = ppair[pp].chit[c].pps[p];
				best = right = ppair[pp].chit[c].pme[0];
				for(m=1; m < ppair[pp].chit[c].pmc; m++) {
					right = ppair[pp].chit[c].pme[m];
					if((left < right) && (right < best)){
						best = right;		
						r1=ppair[pp].chit[c].pms[m];
						r2=ppair[pp].chit[c].pme[m];
					}
				}
				right = best;
				if(((right - left) <= AMPSB)&&((right-left) >0)){
               /* REMOVE? */
               if(dlevel >3) {
						get_hybd_nseq_from_db(ppair[pp].chit[c].cid,ppair[pp].chit[c].pps[p],ppair[pp].chit[c].ppe[p],r1,r2,&ppair[pp].chit[c].hlseq,&ppair[pp].chit[c].hrseq,&ppair[pp].chit[c].hllen,&ppair[pp].chit[c].hrlen);	 
						printf("left: %s\nright: %s\n",ppair[pp].chit[c].hlseq,ppair[pp].chit[c].hrseq);
						printf("left-primer: %s\nright-primer: %s\n",ppair[pp].lseq,ppair[pp].rseq);
						get_oligo_tm_for_hybd(ppair[pp].chit[c].hlseq,ppair[pp].chit[c].hrseq,ppair[pp].lseq,ppair[pp].rseq,&ppair[pp].chit[c].potential);
					}
					if(ppair[pp].chit[c].potential == 1) {
						strcpy(cid, ppair[pp].chit[c].cid);	/*tmp var*/
						strcpy(ppair[pp].hbrdn[ampl].id, cid);
						ppair[pp].hbrdn[ampl].hs = left;
						ppair[pp].hbrdn[ampl].he = right;
						ppair[pp].hbrdn[ampl].hsize = right - left;
						ampl ++;
					}
				}
			}
		}
		if(ampl != ppair[pp].hbrdn_count)
			p_exit("repeat counting of hybridization not same\n");
	}
	printf("%s find all hybridization for ppair...done!\n",T_S());
	return PASS;
}

/* check the hybridizations for a query against chromosomes */
int chk_hybridization_chromosome(int q) {			/*index for qinfo*/
	int pp=0;												/*primer pair index*/
	int c=0;													/*chromosome index*/
	int p=0;													/*plus-plus hit index*/
	int m=0;                                		/*plus-minus hit index*/		
        int r1=0, r2=0;					
	int p_p=0, p_m=0, best=0;							/*temporary buffer*/
	int left=0, right=0;
	int ampl=0;												/*amplicon index*/
	size_t t_size;
	char cid[BISULCHR];

	printf("%s find all hybridization for ppair...\n",T_S());
	/* for each primer pair */
	for(pp=0; pp < ppair_count; pp++) {
		if (dlevel > 3)
			printf("processing primer pair # %d\n", pp+1);
		ppair[pp].hbrdn_count = 0;
		/* count the potential hybridizations */
		for(c=0; c < BLASTRECORDS; c++) {
			if(ppair[pp].chit[c].ppc == 0)
				continue;				/*amplicons not possible*/
			if(ppair[pp].chit[c].pmc == 0)
				continue;				/*amplicons not possible*/
			/* for each plus-plus start */
			for(p=0; p < ppair[pp].chit[c].ppc; p++) {
				left = ppair[pp].chit[c].pps[p];
				best = right = ppair[pp].chit[c].pme[0];
				r1=ppair[pp].chit[c].pms[0];
				r2=ppair[pp].chit[c].pme[0];
				/* for each plus-minus end */
				for(m=1; m < ppair[pp].chit[c].pmc; m++) {
					right = ppair[pp].chit[c].pme[m];
					if((left < right) && (right < best)){
						best = right;
						r1=ppair[pp].chit[c].pms[m];
						r2=ppair[pp].chit[c].pme[m];
					}
				}
				right = best;
				if(((right - left) <= AMPSB)&&((right-left) >0)){
					if(dlevel > 3) {
						/* REMOVE? */
						printf("this could be potential hybridization.\n");
						printf("left:[%d,%d]\t",ppair[pp].chit[c].pps[p],ppair[pp].chit[c].ppe[p]);
						printf("right:[%d,%d]\n", r1, r2);
					}
					get_hybd_from_chromosome(ppair[pp].chit[c].cid,ppair[pp].chit[c].pps[p],ppair[pp].chit[c].ppe[p],r1,r2,ppair[pp].chit[c].cid,&ppair[pp].chit[c].hlseq,&ppair[pp].chit[c].hrseq);
					get_oligo_tm_for_hybd(ppair[pp].chit[c].hlseq,ppair[pp].chit[c].hrseq,ppair[pp].lseq,ppair[pp].rseq,&ppair[pp].chit[c].potential);
					if(ppair[pp].chit[c].potential == 1) {
						ppair[pp].hbrdn_count ++;
					}
				}
			}
		}
		printf("%s total hybridizations for ppair %02d",T_S(),pp+1);
		printf(" [%s]\t[%s]",ppair[pp].lseq, ppair[pp].rseq);
		printf("...%d\n", ppair[pp].hbrdn_count);
		
		/* store the potential hybridizations from above */
		t_size = (ppair[pp].hbrdn_count)*sizeof(Hbrdn);
		ppair[pp].hbrdn=(Hbrdn *)p_malloc(t_size);
		ampl = 0;
		for(c=0; c < BLASTRECORDS; c++) {
			if(ppair[pp].chit[c].ppc == 0)
				continue;		/*amplicons not possible*/
			if(ppair[pp].chit[c].pmc == 0)
				continue;		/*amplicons not possible*/
			for(p=0; p < ppair[pp].chit[c].ppc; p++) {
				left = ppair[pp].chit[c].pps[p];
				best = right = ppair[pp].chit[c].pme[0];
				r1=ppair[pp].chit[c].pms[0];
				r2=ppair[pp].chit[c].pme[0];
				for(m=1; m < ppair[pp].chit[c].pmc; m++) {
					right = ppair[pp].chit[c].pme[m];
					if((left < right) && (right < best)){
						best = right;
						r1=ppair[pp].chit[c].pms[m];
						r2=ppair[pp].chit[c].pme[m];
					}
				}
				right = best;
				if(((right - left) <= AMPSB)&&((right-left) >0)){
					get_hybd_from_chromosome(ppair[pp].chit[c].cid,ppair[pp].chit[c].pps[p],ppair[pp].chit[c].ppe[p],r1,r2,ppair[pp].chit[c].cid,&ppair[pp].chit[c].hlseq,&ppair[pp].chit[c].hrseq);
					get_oligo_tm_for_hybd(ppair[pp].chit[c].hlseq,ppair[pp].chit[c].hrseq,ppair[pp].lseq,ppair[pp].rseq,&ppair[pp].chit[c].potential);
					/* REMOVE? */
					if(dlevel>3)
						printf("potl=%d\n",ppair[pp].chit[c].potential);
					if(ppair[pp].chit[c].potential == 1) {
						strcpy(cid, ppair[pp].chit[c].cid);	/*tmp var*/
						strcpy(ppair[pp].hbrdn[ampl].id, cid);
						ppair[pp].hbrdn[ampl].hs = left;
						ppair[pp].hbrdn[ampl].he = right;
						ppair[pp].hbrdn[ampl].hsize = right - left;
						ampl++;
					}
				}
			}
		}
		if(ampl != ppair[pp].hbrdn_count)
			p_exit("repeat counting of hybridization not same\n");
	}
	printf("%s find all hybridization for ppair...done!\n",T_S());
	return PASS;
}
