
/*
 * primer design module, used for primer design using various different 
 * programs.
 *
 */
 

#include "lib.h"
#include "defs.h"


short run_primer3(Query qry, int start, int end);
short run_primer3_test(Query qry, int start, int end);

short find_maxCutSite_region(Query qry, int *start, int *length);
short count_cutSite_region(char *nseq, int start, int length);

short save_primer3_output(void);
short save_primer3_output_test(void);


short find_maxCutSite_region(Query qry, int *start, int *length) {
	
	int length, b=0, c=0, currCount=0, maxCount=0i, prevCount=0;
	int target_start=0, target_length=0;
	short cutSiteFlag=0;
	char **cutSite, line[SLINE];
	char *tok;

	if(CUT_SITE_REGION==0)
		p_exit("cut-site region is not given!\n");
	target_length = CUT_SITE_REGION;	
	




	/*get all cutsite info*/
	cutSite = p_malloc(CUT_SITE_COUNT*sizeof(char *));
	strcpy(line, CUT_SITE);
	for(c=0; c < CUT_SITE_COUNT; c++) {
		tok = strtok(line, ",:; \t\n");
		if(tok==NULL)
			p_exit("cut-site is null!\n");	
		strcat(tok, "\0");
		if(strlen(tok) >4)
			p_exit("cut-site is incorrect!\n");
		strcpy(cutSite[c], tok);	
	}
	/*cut site count initialization*/
	currCount=0;
	for(b=0; b < target_length; b++) {
		cutSiteFlag=0;	
		for(c=0; c < CUT_SITE_COUNT; c++) {
			if(strncmp(qry.nseq+b,cutSite[c],strlen(cutSite[c]))==0) {
				cutSiteFlag = 1;
				break;	
			}	
		}
		if(cutSiteFlag==1)
			currCount++;	
	}
	maxCount = currCount;
	prevCount = currCount;
	
	for(b=1; b <= length-target_legnth; b++) {
		/*left side region going out*/
		cutSiteFlag=0;	
		for(c=0; c < CUT_SITE_COUNT; c++) {
			if(strncmp(qry.nseq-b-1,cutSite[c],strlen(cutSite[c]))==0) {
				cutSiteFlag = 1;
				break;	
			}	
		}
		if(cutSiteFlag==1)
			currCount = prevCount -1;	
		/*right side region coming in*/	
		cutSiteFlag=0;	
		for(c=0; c < CUT_SITE_COUNT; c++) {
			if(strncmp(qry.nseq+b+target_legnth-2,cutSite[c],strlen(cutSite[c]))==0) {
				cutSiteFlag = 1;
				break;	
			}	
		}
		if(cutSiteFlag==1)
			currCount = prevCount -1;	
	}
	return PASS;
}


short count_cutSite_region(char *nseq, int start, int length) {

	return PASS;	
}

short run_primer3(Query qry, int start, int end) {
	FILE *fp;
	int length, b=0, c=0, currCount=0, maxCount=0i, prevCount=0;
	int target_start=0, target_length=0;
	short cutSiteFlag=0;
	char cmd[COMMAND], **cutSite, line[SLINE];
	char *tok;

	printf("%s run primer3 for %.30s...\n", T_S(), qry.id);
	length = end - start +1;
	fp = p_fopen_error_exit(fn_p3in, "w");
	fprintf(fp,"PRIMER_SEQUENCE_ID=%s\n",qry.id);
	fprintf(fp,"SEQUENCE=%s\n",qry.nseq);
	fprintf(fp,"INCLUDED_REGION=%d,%d\n", start, length);
	fprintf(fp,"PRIMER_PRODUCT_SIZE_RANGE=%s\n", PPSR);
	if(strlen(TARGET)>2)
		fprintf(fp, "TARGET=%s\n", TARGET);
	fprintf(fp,"PRIMER_EXPLAIN_FLAG=%d\n", PEF);
	fprintf(fp,"PRIMER_FILE_FLAG=%d\n", PFF);
	fprintf(fp,"PRIMER_NUM_RETURN=%d\n", PNR);
	fprintf(fp,"PRIMER_MAX_POLY_X=%d\n", PMPX);
	fprintf(fp,"PRIMER_INTERNAL_OLIGO_MAX_POLY_X=%d\n", PIOMPX);
	fprintf(fp,"PRIMER_MIN_TM=%f\n", PMINT);
	fprintf(fp,"PRIMER_OPT_TM=%f\n", POT);
	fprintf(fp,"PRIMER_MAX_TM=%f\n", PMAXT);
	fprintf(fp,"PRIMER_MIN_SIZE=%d\n", PMINS);
	fprintf(fp,"PRIMER_OPT_SIZE=%d\n", POS);
	fprintf(fp,"PRIMER_MAX_SIZE=%d\n", PMAXS);
	fprintf(fp,"PRIMER_MAX_GC=%d\n", PMAXGC);
	fprintf(fp,"PRIMER_MIN_GC=%d\n", PMINGC);
	fprintf(fp,"PRIMER_SELF_END=%d\n", PSE);
	fprintf(fp,"PRIMER_SELF_ANY=%d\n", PSA);
	fprintf(fp,"PRIMER_MAX_END_STABILITY=%f\n", PMES);
	fprintf(fp,"=\n");
	fclose(fp);
	
	sprintf(cmd,"%s < %s > %s",primer3, fn_p3in, fn_p3out);
	p_exec(cmd);
	printf("%s run primer3 for %.30s...done!\n", T_S(), qry.id);
	return PASS; 
}


short run_primer3_test(Query qry, int start, int end) {

	return PASS;
}


short save_primer3_output(void) {
	FILE *fp;
	int p=0, count=0, b=0;		/*nucleotide ndx for primer filter*/
	int cg_count=0;
	int lst=0, rst=0, llen=0, rlen=0, psize=0;
	double ltm=0.0, rtm=0.0, lgc=0.0, rgc=0.0;
	char line[SLINE], *tok;
	char plseq[PRIMER3TAG];
	char prseq[PRIMER3TAG];
	char pl[PRIMER3TAG];
	char pr[PRIMER3TAG];
	char pltm[PRIMER3TAG];
	char prtm[PRIMER3TAG];
	char plgc[PRIMER3TAG];
	char prgc[PRIMER3TAG];
	char pps[PRIMER3TAG];
	char perror[PRIMER3TAG];
	PPair *tppair;
	int tppair_count;
	
	printf("%s read primer3 output from %s file",T_S(), fn_p3out);
	printf(" and save all primer pairs info into memory\n");
	strcpy(plseq, "PRIMER_LEFT_SEQUENCE=");
	strcpy(prseq, "PRIMER_RIGHT_SEQUENCE=");
	strcpy(pl, "PRIMER_LEFT=");
	strcpy(pr, "PRIMER_RIGHT=");
	strcpy(pltm, "PRIMER_LEFT_TM=");
	strcpy(prtm, "PRIMER_RIGHT_TM=");
	strcpy(plgc, "PRIMER_LEFT_GC_PERCENT=");
	strcpy(prgc, "PRIMER_RIGHT_GC_PERCENT=");
	strcpy(pps, "PRIMER_PRODUCT_SIZE=");
	strcpy(perror, "PRIMER_ERROR=");
	
	fp = p_fopen_error_exit(fn_p3out, "r");
	count = 0;
	while(fgets(line, SLINE, fp)) { 
		if(strncmp(line, pps, strlen(pps)-1)==0) 	/*next pair info*/
			count++;
		if(strncmp(line, perror, strlen(perror))==0) {
			fprintf(stdout, "%s primer design error:", T_S());
			tok = strtok(line, "=");				/*remove tag*/
			tok = strtok(NULL, "\n");
			strcat(tok, "\0");
			fprintf(stdout, " %s\n", tok);
		}
	}
	tppair_count = --count;				/*primer pair count*/
	printf("%s total primer pairs designed = %d\n", T_S(), count);
	rewind(fp);
	if(tppair_count <= 0) {
		printf("%s primer design failed!...query skipped!\n",T_S());
		return FAIL;
	}
	tppair = (PPair *)malloc(tppair_count*sizeof(PPair));
	count = 0;
	while(fgets(line, SLINE, fp)) {
		if(strncmp(line, plseq, strlen(plseq))==0) {
			tok = strtok(line, "=");				/*remove tag*/
			strcat(tok, "\0");
			strcpy(tppair[count].lid, tok);
			tok = strtok(NULL, "\n");
			strcat(tok, "\0");
			strcpy(tppair[count].lseq, tok);
		}
		if(strncmp(line, prseq, strlen(prseq))==0) {
			tok = strtok(line, "=");				/*remove tag*/
			strcat(tok, "\0");
			strcpy(tppair[count].rid, tok);
			tok = strtok(NULL, "\n");
			strcat(tok, "\0");
			strcpy(tppair[count].rseq, tok);
		}
		if(strncmp(line, pl, strlen(pl))==0) {
			tok = strtok(line, "=");				/*remove tag*/
			tok = strtok(NULL, ",");
			strcat(tok, "\0");
			tppair[count].lst = atoi(tok);
			tok = strtok(NULL, "\n");
			strcat(tok, "\0");
			tppair[count].llen = atoi(tok);
		}
		if(strncmp(line, pr, strlen(pr))==0) {
			tok = strtok(line, "=");				/*remove tag*/
			tok = strtok(NULL, ",");
			strcat(tok, "\0");
			tppair[count].rst = atoi(tok);
			tok = strtok(NULL, "\n");
			strcat(tok, "\0");
			tppair[count].rlen = atoi(tok);
		}
		if(strncmp(line, pltm, strlen(pltm))==0) {
			tok = strtok(line, "=");				/*remove tag*/
			tok = strtok(NULL, "\n");
			strcat(tok, "\0");
			tppair[count].ltm = atof(tok);
		}
		if(strncmp(line, prtm, strlen(prtm))==0) {
			tok = strtok(line, "=");				/*remove tag*/
			tok = strtok(NULL, "\n");
			strcat(tok, "\0");
			tppair[count].rtm = atof(tok);
		}
		if(strncmp(line, plgc, strlen(plgc))==0) {
			tok = strtok(line, "=");				/*remove tag*/
			tok = strtok(NULL, "\n");
			strcat(tok, "\0");
			tppair[count].lgc = atof(tok);
		}
		if(strncmp(line, prgc, strlen(prgc))==0) {
			tok = strtok(line, "=");				/*remove tag*/
			tok = strtok(NULL, "\n");
			strcat(tok, "\0");
			tppair[count].rgc = atof(tok);
		}
		if(strncmp(line, pps, strlen(pps))==0) {	/*next pair info*/
			tok = strtok(line, "=");				/*remove tag*/
			tok = strtok(NULL, "\n");
			strcat(tok, "\0");
			tppair[count].psize = atoi(tok);
			count++;
			sprintf(plseq,"PRIMER_LEFT_%d_SEQUENCE=",count);
			sprintf(prseq,"PRIMER_RIGHT_%d_SEQUENCE=",count);
			sprintf(pl,"PRIMER_LEFT_%d=",count);
			sprintf(pr,"PRIMER_RIGHT_%d=",count);
			sprintf(pltm,"PRIMER_LEFT_%d_TM=",count);
			sprintf(prtm,"PRIMER_RIGHT_%d_TM=",count);
			sprintf(plgc,"PRIMER_LEFT_%d_GC_PERCENT=",count);
			sprintf(prgc,"PRIMER_RIGHT_%d_GC_PERCENT=",count);
			sprintf(pps,"PRIMER_PRODUCT_SIZE_%d=",count);
		}
	}
	fclose(fp);
	
	/*filtering primer pairs based on PRIMER_CG_MAX parameters*/
	if(PRIMER_CG_MAX > 0) {
		printf("%s filtering primer pairs on CG count...\n",T_S());
		/*count total selected primer pair*/
		count=0;
		for(p=0; p < tppair_count; p++) {
			cg_count=0;
			for(b=0; b < tppair[p].llen; b++) { /*forward primer:CG*/
				if((toupper(tppair[p].lseq[b])=='C')
								&& (toupper(tppair[p].lseq[b+1])=='G'))
						cg_count++;
			}
			if(cg_count > PRIMER_CG_MAX)
				continue;
			cg_count=0;
			for(b=0; b < tppair[p].rlen; b++) { /*reverse primer:CG*/
				if((toupper(tppair[p].rseq[b])=='C')
								&& (toupper(tppair[p].rseq[b+1])=='G'))
						cg_count++;
			}
			if(cg_count > PRIMER_CG_MAX)
				continue;
			count++;
		}
		printf("%s total filtered primer pair is %d\n",T_S(),count);
		ppair_count = count;
		if(ppair_count <= 0) {
			printf("%s primer design failed!...query skipped!\n",T_S());
			return FAIL;
		}
		ppair = (PPair *)p_malloc(ppair_count*sizeof(PPair));
		count=0;
		for(p=0; p < tppair_count; p++) {
			cg_count=0;
			for(b=0; b < tppair[p].llen; b++) { /*forward primer:CG*/
				if((toupper(tppair[p].lseq[b])=='C')
								&& (toupper(tppair[p].lseq[b+1])=='G'))
						cg_count++;
			}
			if(cg_count > PRIMER_CG_MAX)
				continue;
			cg_count=0;
			for(b=0; b < tppair[p].rlen; b++) { /*reverse primer:CG*/
				if((toupper(tppair[p].rseq[b])=='C')
								&& (toupper(tppair[p].rseq[b+1])=='G'))
						cg_count++;
			}
			if(cg_count > PRIMER_CG_MAX)
				continue;
			strcpy(ppair[count].lid, tppair[p].lid);
			strcpy(ppair[count].rid, tppair[p].rid);
			strcpy(ppair[count].lseq, tppair[p].lseq);
			strcpy(ppair[count].rseq, tppair[p].rseq);
			ppair[count].lst = tppair[p].lst;
			ppair[count].rst = tppair[p].rst;
			ppair[count].llen = tppair[p].llen;
			ppair[count].rlen = tppair[p].rlen;
			ppair[count].ltm = tppair[p].ltm;
			ppair[count].rtm = tppair[p].rtm;
			ppair[count].lgc = tppair[p].lgc;
			ppair[count].rgc = tppair[p].rgc;
			ppair[count].psize = tppair[p].psize;
			count++;
		}
		printf("%s filtering primer pairs on CG count...done!\n",T_S());
	} else {
		ppair_count = tppair_count;
		ppair = (PPair *)p_malloc(ppair_count*sizeof(PPair));
		for(p=0; p < tppair_count; p++) {
			strcpy(ppair[p].lid, tppair[p].lid);
			strcpy(ppair[p].rid, tppair[p].rid);
			strcpy(ppair[p].lseq, tppair[p].lseq);
			strcpy(ppair[p].rseq, tppair[p].rseq);
			ppair[p].lst = tppair[p].lst;
			ppair[p].rst = tppair[p].rst;
			ppair[p].llen = tppair[p].llen;
			ppair[p].rlen = tppair[p].rlen;
			ppair[p].ltm = tppair[p].ltm;
			ppair[p].rtm = tppair[p].rtm;
			ppair[p].lgc = tppair[p].lgc;
			ppair[p].rgc = tppair[p].rgc;
			ppair[p].psize = tppair[p].psize;
		}
	}
	free(tppair);
	for(p=0; p < ppair_count; p++) {
		printf("%s %03d) ", T_S(), p+1);
		printf("%26s  %-30s ", ppair[p].lid, ppair[p].lseq);
		printf("%26s  %-30s", ppair[p].rid, ppair[p].rseq);
		printf("  psize %d\n", ppair[p].psize);
	}
	printf("%s read & save primers info from ", T_S());
	printf("%s into memory...done!\n", fn_p3out);
	//save_primer3_output_test();
	return PASS;
}


short save_primer3_output_test(void) {
	int p=0;

	for(p=0; p < ppair_count; p++) {
		printf("%2d) lseq\t%s\n", p+1, ppair[p].lseq);
		printf("%2d) rseq\t%s\n", p+1, ppair[p].rseq);
		printf("%2d) lst,len\t%d,%d\n",p+1,ppair[p].lst,ppair[p].llen);
		printf("%2d) rst,len\t%d,%d\n",p+1,ppair[p].rst,ppair[p].rlen);
		printf("%2d) ltm\t%f\n", p+1, ppair[p].ltm);
		printf("%2d) rtm\t%f\n", p+1, ppair[p].rtm);
		printf("%2d) lgc\t%f\n", p+1, ppair[p].lgc);
		printf("%2d) rgc\t%f\n", p+1, ppair[p].rgc);
		printf("%2d) psize\t%d\n\n", p+1, ppair[p].psize);
	}
	
	return PASS;
}


