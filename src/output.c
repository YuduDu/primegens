/*
 * This module contains all the funtions dealing with output, including
 * sorting the results and printing them to the user specified files.  
 */

#include "lib.h"
#include "defs.h"
#include "output.h"
#include "util.h"

/* sort the primer pairs based on hybridization */
int rank_ppair(void) {
	
	qsort(ppair, ppair_count, sizeof(PPair), sort_on_hybrid);
	return PASS;
}

/* a function to compare the # of hybridations of two primer pairs */
int sort_on_hybrid(const void *ptr1, const void *ptr2) {
 	 	 PPair *p1, *p2;

		 p1= (PPair *)ptr1;
		 p2= (PPair *)ptr2;
		 if(p1->hbrdn_count > p2->hbrdn_count)
			 return 1;
		 else if(p1->hbrdn_count < p2->hbrdn_count)
			 return -1;
		 else 
			 return 0;
}

/*print query q output*/
int print_output(int q) { 
	int pp=0, h=0, b=0;
	FILE *fout1;
	FILE *fout2;
	FILE *fout3;
	FILE *fout5;		/*4 is already used for failed query*/

	fout1 = p_fopen_error_exit(fn_out1, "a");
	fout2 = p_fopen_error_exit(fn_out2, "a");
	fout3 = p_fopen_error_exit(fn_out3, "a");
	
	/*if no query sequence in input, generate fasta file in output*/
	if(QUERY_SEQUENCE==0) {
		fout5 = p_fopen_error_exit(fn_out5, "a");
		fprintf(fout5, "%s\n", query[q].id);
		for(b=0; b < query[q].len; b++) {
			fprintf(fout5, "%c", query[q].nseq[b]);
			if((b + 1) % FASTALINE == 0)
				fprintf(fout5, "\n");
		}
		if(query[q].len % FASTALINE != 0)
			fprintf(fout5, "\n");
		fclose(fout5);
		printf("%s adding query into %s...done!\n",T_S(),fn_out5);
	}

	if(query[q].res==PASS) {	/*after primer pair sorting*/
		/*first output format-excel sheet*/
if (CUT_SITE_COUNT==0){
		if((QUERY_FORMAT==3)||(QUERY_FORMAT==4)) {
			fprintf(fout1, "\n>%s:", qinfo[q].chr);
			fprintf(fout1, "%d-", qinfo[q].start);
			fprintf(fout1, "%d\t", qinfo[q].end);
			if(QUERY_FORMAT==4)
				fprintf(fout1, "%s\n", qinfo[q].desc);
			else 
				fprintf(fout1, "\n");
			fprintf(fout1, "%d\t", query[q].len);
		 } else if((QUERY_FORMAT==1)||(QUERY_FORMAT==2)) {
			fprintf(fout1, "\n%s\t", query[q].id);
		//	if(QUERY_FORMAT==2)
		//		fprintf(fout1, "%s\t", qinfo[q].desc);
		 }
}
		fprintf(fout1, "%s\t", ppair[0].lseq);
		fprintf(fout1, "%d\t", ppair[0].lst);
		fprintf(fout1, "%d\t", ppair[0].llen);
		fprintf(fout1, "%f\t", ppair[0].ltm);
		fprintf(fout1, "%f\t", ppair[0].lgc);
		fprintf(fout1, "%s\t", ppair[0].rseq);
		fprintf(fout1, "%d\t", ppair[0].rst);
		fprintf(fout1, "%d\t", ppair[0].rlen);
		fprintf(fout1, "%f\t", ppair[0].rtm);
		fprintf(fout1, "%f\t", ppair[0].rgc);
		fprintf(fout1, "%d\t", ppair[0].psize);
		fprintf(fout1, "%d", ppair[0].hbrdn_count);
		
		/*second output format-primer list*/
		if((QUERY_FORMAT==3)||(QUERY_FORMAT==4)) {
			fprintf(fout2, "%s:\t", qinfo[q].chr);
			fprintf(fout2, "%d-\t", qinfo[q].start);
			fprintf(fout2, "%d:\t", qinfo[q].end);
			if(QUERY_FORMAT==4)
				fprintf(fout2, "%s:\n", qinfo[q].desc);
			else 
				fprintf(fout2, "\n");
			fprintf(fout2, "%d\t", query[q].len);
		 } else if((QUERY_FORMAT==1)||(QUERY_FORMAT==2)) {
			fprintf(fout2, "%s\t", query[q].id);
			if(QUERY_FORMAT==2)
				fprintf(fout2, "%s", qinfo[q].desc);
			fprintf(fout2, "\n");
		 }
		for(pp=0; pp < p_min(ppair_count,PRIMER_DISPLAY_COUNT); pp++) {
			fprintf(fout2, "\t%2d) %-30s", pp+1, ppair[pp].lseq);
			fprintf(fout2, "[%5d]  ", ppair[pp].lst);
			fprintf(fout2, "%-30s", ppair[pp].rseq);
			fprintf(fout2, "[%5d]  ", ppair[pp].rst);
			fprintf(fout2, "psize %4d  ", ppair[pp].psize);
			fprintf(fout2, "hbrdn %4d\n", ppair[pp].hbrdn_count);
		}
		
		/*third output format-query primer list*/
		if((QUERY_FORMAT==3)||(QUERY_FORMAT==4)) {
			fprintf(fout3, "%s:\t", qinfo[q].chr);
			fprintf(fout3, "%d-\t", qinfo[q].start);
			fprintf(fout3, "%d:\t", qinfo[q].end);
			if(QUERY_FORMAT==4)
				fprintf(fout3, "%s:\n", qinfo[q].desc);
			else 
				fprintf(fout3, "\n");
			fprintf(fout3, "%d\t", query[q].len);
		 } else if((QUERY_FORMAT==1)||(QUERY_FORMAT==2)) {
			fprintf(fout3, "%s\t", query[q].id);
			if(QUERY_FORMAT==2)
				fprintf(fout3, "%s", qinfo[q].desc);
			fprintf(fout3, "\n");
		 }
		fprintf(fout3, "%s\n", query[q].nseq);
		for(pp=0; pp < p_min(ppair_count,PRIMER_DISPLAY_COUNT); pp++) {
			fprintf(fout3, "\t%2d) %-30s", pp+1, ppair[pp].lseq);
			fprintf(fout3, "[%5d]  ", ppair[pp].lst);
			fprintf(fout3, "%-30s", ppair[pp].rseq);
			fprintf(fout3, "[%5d]  ", ppair[pp].rst);
			fprintf(fout3, "psize %4d  ", ppair[pp].psize);
			if(HBRDN_DISPLAY_COUNT==0)
				fprintf(fout3, "hbrdn %4d\n", ppair[pp].hbrdn_count);
			else {
				fprintf(fout3, "hbrdn %4d\t", ppair[pp].hbrdn_count);
				/* REMOVE? */
				//for(h=0; h < p_min(ppair[pp].hbrdn_count, HBRDN_DISPLAY_COUNT); h++)
				for(h=0; h < ppair[pp].hbrdn_count; h++)
					fprintf(fout3, "%s;", ppair[pp].hbrdn[h].id);
				//if(ppair[pp].hbrdn_count > HBRDN_DISPLAY_COUNT)
					//fprintf(fout3, "...\n");
			//	else
					fprintf(fout3, "\n");
			}
		}
	} 
	
	fclose(fout1);
	fclose(fout2);
	fclose(fout3);

	return PASS;
}
