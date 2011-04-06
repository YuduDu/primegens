/*
 * primer design module, used for primer design using various different 
 * programs.
 *
 */
 

#include "lib.h"
#include "defs.h"
#include "util.h"
#include "pd.h"

/* find the max cutsite in the query, storing the start and length by reference */
short find_maxCutSite_region(Query qry, int *start, int *length) {
	int c=0, b=0, cutSiteCounter=0, *cutSitePos;
	short cutSiteFlag=0;
	int currCount=0, maxCount=0, trgt_start=0, trgt_length=0, trgt_end=0;
	
	FILE *fxl;        /*to write in xls file*/
	fxl = p_fopen_error_exit(fn_out1, "a");

	fprintf(fxl, "\n%s\t", qry.id);

	if(dlevel >=3)
		printf("%s find max cut-site region in query sequence...\n",T_S());	
	/*get all cut-site counts in query sequence*/
	cutSiteCounter=0;
	for(b=0; b < qry.len; b++) {
		cutSiteFlag=0;	
		for(c=0; c < cutSite_count; c++) {
			if(strncmp(qry.nseq+b,cutSite[c].nseq,strlen(cutSite[c].nseq))==0) {
				cutSiteFlag = 1;
				break;	
			}	
		}
		if(cutSiteFlag==1)
			cutSiteCounter++;			
	}
	printf("%s total cutsites in query is %d.\n",T_S(),cutSiteCounter);
	cutSitePos = (int *)p_malloc(cutSiteCounter*sizeof(int));
	
	/*get all cut-site position in query sequence*/
	if(dlevel >=3)
		printf("%s cut-site positions:", T_S());	
	currCount = 0;
	for(b=0; b < qry.len; b++) {
		for(c=0; c < cutSite_count; c++) {
			if(strncmp(qry.nseq+b,cutSite[c].nseq,strlen(cutSite[c].nseq))==0) {
				cutSitePos[currCount] = b+1; 	
				currCount++;
				if(dlevel >=3)
					printf(" %d(%s)", b+1, cutSite[c].nseq);	
				break;	
			}	
		}
	}

	if(dlevel >=3)
		printf("...\n");
	if(currCount != cutSiteCounter)
		p_exit("problem while locating cut-sites in query sequence!\n");

	currCount = maxCount = 0;
	trgt_length = CUT_SITE_REGION;
	
	/* find the max number of cut sites */
	for(b=0; b < cutSiteCounter; b++) {
		trgt_start = cutSitePos[b];
		trgt_end = trgt_start + trgt_length-1;
		currCount=1;
		for(c=b+1; c < cutSiteCounter; c++) {
			if(cutSitePos[c] <= trgt_end)	
				currCount++;	
		}
		if(currCount > maxCount)
			maxCount = currCount;	
	}
	if(dlevel >=3)
		printf("%s max-cutsite region has %d sites.\n",T_S(), maxCount);		
		
	/* find the max cut site */
	for(b=0; b < cutSiteCounter; b++) {
		trgt_start = cutSitePos[b];
		trgt_end = trgt_start + trgt_length-1;
		currCount=1;
		for(c=b+1; c < cutSiteCounter; c++) {
			if(cutSitePos[c] <= trgt_end)	
				currCount++;	
		}
		if(currCount==maxCount) {
			*start = trgt_start;
			for(c=b+1; c < cutSiteCounter; c++) {
				if(cutSitePos[c] > trgt_end)	
					break;
			}	
			*length = cutSitePos[c-1] - trgt_start+1;
			break;
		}	
	}
	if(dlevel >=3){
		printf("%s target start=%d length=%d\n",T_S(), *start, *length);	
		printf("%s cut-sites within target region are:\n",T_S());	}
	
	fprintf(fxl,"start=%d length=%d\t",*start,*length);	
	fprintf(fxl,"%d\t",maxCount);	

	int st = *start;
	int ln = *length;
	int en = st+ln;

	printf("%s start:%d\n",T_S(),st);	
	printf("%s length:%d\n",T_S(),ln);	
        
	for(b=st-1;b<en;b++)
	{
          for(c=0; c < cutSite_count; c++) {
		      if(strncmp(qry.nseq+b,cutSite[c].nseq,strlen(cutSite[c].nseq))==0) {
			     cutSitePos[currCount] = b;
				 currCount++;
		      if(dlevel >=3)
			        printf(" %d(%s)", b, cutSite[c].nseq);
			        fprintf(fxl," %d(%s)", b, cutSite[c].nseq);
			     break;
	          }
	      }
    }
	fprintf(fxl,"\t");
	fclose(fxl);
	return PASS;
    }  
 
/* Run primer3 for the given query */
short run_primer3(Query qry) {
	FILE *fp;
	

	int length, b=0, c=0, currCount=0, maxCount=0i, prevCount=0;
	int target_start=0, target_length=0;
	short cutSiteFlag=0;
	char cmd[COMMAND], **cutSite, line[SLINE];
	char *tok;


	/* set up the parameters for primer3 */
	printf("%s run primer3 for %.30s...\n", T_S(), qry.id);
	fp = p_fopen_error_exit(fn_p3in, "w");
	fprintf(fp,"PRIMER_SEQUENCE_ID=%s\n",qry.id);
	fprintf(fp,"SEQUENCE=%s\n",qry.nseq);
	
	if(strcmp(PIR,"unassigned")!=0){
	   fprintf(fp,"INCLUDED_REGION=%s\n", PIR);}
	if(strcmp(PER,"unassigned")!=0){
	   fprintf(fp,"EXCLUDED_REGION=%s\n", PER);}
	if(strcmp(PSQ,"unassigned")!=0){
	   fprintf(fp,"PRIMER_SEQUENCE_QUALITY=%s\n", PSQ);}
	       
	fprintf(fp,"PRIMER_PRODUCT_SIZE_RANGE=%s\n", PPSR);
	
        if(strcmp(TARGET,"unassigned")!=0){
            if(strlen(TARGET)>2)
	    	fprintf(fp, "TARGET=%s\n", TARGET);}
        
        if(CUT_SITE_COUNT > 0) {
		find_maxCutSite_region(qry, &target_start, &target_length);
		fprintf(fp, "TARGET=%d,%d\n", target_start, target_length);
	    }
        
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
	fprintf(fp,"PRIMER_MAX_GC=%f\n", PMAXGC);
	fprintf(fp,"PRIMER_OPT_GC_PERCENT=%f\n", POPTGC);
	fprintf(fp,"PRIMER_MIN_GC=%f\n", PMINGC);
	fprintf(fp,"PRIMER_SELF_END=%d\n", PSE);
	fprintf(fp,"PRIMER_SELF_ANY=%d\n", PSA);
	fprintf(fp,"PRIMER_MAX_END_STABILITY=%f\n", PMES);
	
        fprintf(fp,"PRIMER_START_CODON_POSITION=%d\n", PSCP);

    fprintf(fp,"PRIMER_PRODUCT_MAX_TM=%f\n", PPMAXT);
    fprintf(fp,"PRIMER_PRODUCT_MIN_TM=%f\n", PPMINT);
    fprintf(fp,"PRIMER_GC_CLAMP=%d\n", PGC);
    fprintf(fp,"PRIMER_DEFAULT_SIZE=%d\n", PDS);
    fprintf(fp,"PRIMER_MAX_DIFF_TM=%f\n", PMDT);
    fprintf(fp,"PRIMER_TM_SANTALUCIA=%d\n", PTS);
    fprintf(fp,"PRIMER_SALT_CONC=%f\n", PSC);
    fprintf(fp,"PRIMER_DIVALENT_CONC=%f\n", PDIVC);
    fprintf(fp,"PRIMER_DNTP_CONC=%f\n", PDNTPC);
    fprintf(fp,"PRIMER_SALT_CORRECTIONS=%f\n", PSALTC);
    fprintf(fp,"PRIMER_LOWERCASE_MASKING=%d\n", PLM);
    fprintf(fp,"PRIMER_DNA_CONC=%f\n", PDNAC);
    fprintf(fp,"PRIMER_NUM_NS_ACCEPTED=%d\n", PNNSA);
    fprintf(fp,"PRIMER_LIBERAL_BASE=%d\n", PLB);
    fprintf(fp,"PRIMER_FIRST_BASE_INDEX=%d\n", PFBI);
    fprintf(fp,"PRIMER_MIN_QUALITY=%d\n", PMQ);
    fprintf(fp,"PRIMER_MIN_END_QUALITY=%d\n", PMEQ);
    fprintf(fp,"PRIMER_QUALITY_RANGE_MIN=%d\n", PQRMIN);
    fprintf(fp,"PRIMER_QUALITY_RANGE_MAX=%d\n", PQRMAX);
    fprintf(fp,"PRIMER_INSIDE_PENALTY=%f\n", PIP);
    fprintf(fp,"PRIMER_OUTSIDE_PENALTY=%f\n", POP);
    fprintf(fp,"PRIMER_PRODUCT_OPT_TM=%f\n", PPOTM);
    fprintf(fp,"PRIMER_PRODUCT_OPT_SIZE=%d\n", PPOS);


    fprintf(fp,"PRIMER_WT_TM_GT=%f\n", PWTG);
    fprintf(fp,"PRIMER_WT_TM_LT=%f\n", PWTL);
    fprintf(fp,"PRIMER_WT_SIZE_LT=%f\n", PWSL);
    fprintf(fp,"PRIMER_WT_SIZE_GT=%f\n", PWSG);
    fprintf(fp,"PRIMER_WT_GC_PERCENT_LT=%f\n", PWGPL);
    fprintf(fp,"PRIMER_WT_GC_PERCENT_GT=%f\n", PWGPG);
    fprintf(fp,"PRIMER_WT_COMPL_ANY=%f\n", PWCA);
    fprintf(fp,"PRIMER_WT_COMPL_END=%f\n", PWCE);
    fprintf(fp,"PRIMER_WT_NUM_NS=%f\n", PWNNS);
    fprintf(fp,"PRIMER_WT_REP_SIM=%f\n", PWRS);
    fprintf(fp,"PRIMER_WT_SEQ_QUAL=%f\n", PWSQ);
    fprintf(fp,"PRIMER_WT_END_QUAL=%f\n", PWEQ);
    fprintf(fp,"PRIMER_WT_POS_PENALTY=%f\n", PWPP);
    fprintf(fp,"PRIMER_WT_END_STABILITY=%f\n", PWES);
    fprintf(fp,"PRIMER_PAIR_WT_PR_PENALTY=%f\n", PPWPP);
    fprintf(fp,"PRIMER_PAIR_WT_DIFF_TM=%f\n", PPWDTM);
    fprintf(fp,"PRIMER_PAIR_WT_COMPL_ANY=%f\n", PPWCA);
    fprintf(fp,"PRIMER_PAIR_WT_COMPL_END=%f\n", PPWCE);
    fprintf(fp,"PRIMER_PAIR_WT_PRODUCT_TM_LT=%f\n", PPWPTL);
    fprintf(fp,"PRIMER_PAIR_WT_PRODUCT_TM_GT=%f\n", PPWPTG);
    fprintf(fp,"PRIMER_PAIR_WT_PRODUCT_SIZE_LT=%f\n", PPWPSL);
    fprintf(fp,"PRIMER_PAIR_WT_PRODUCT_SIZE_GT=%f\n", PPWPSG);
    fprintf(fp,"PRIMER_PAIR_WT_REP_SIM=%f\n", PPWRS);

    fprintf(fp,"=\n");

	fclose(fp);
	
	/* run primer3 */
	sprintf(cmd,"%s < %s > %s",primer3, fn_p3in, fn_p3out);
	p_exec(cmd);
	printf("%s run primer3 for %.30s...done!\n", T_S(), qry.id);
	return PASS; 
}

/* Read the primer3 output file */
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
	/* read the primer pairs */
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
	
	/* read and store primer pair info */
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
			/* store only pairs with less then the max # if CGs */
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
		/* don't care about CG, just store all pairs */
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
	return PASS;
}
