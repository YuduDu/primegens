/*
 * This module contains following functionlities
 * 1- formatting various databases for using blast
 * 2- find unique oligos from the pair pairs 
 * 3- run blast for unique oligos and save results
 *
 */

#include "lib.h"
#include "defs.h"
#include "uhit.h"
#include "util.h"

/* Put all the chromosome info into a virtual database */
int format_chromosome_database(int bloop, int cloop) {
	char tdbname[FILENAME];
	char cmd[FCOMMAND];
	int d=0, c=0;

	printf("%s build and format virtual database...\n",T_S());
	if(cloop != (bloop*CHRDBINLOOP))
		p_exit("format database misordered!\n");
	sprintf(cmd, "cat");
	for(d=0, c=cloop; d < CHRDBINLOOP; d++, c++) {
		printf("%s adding database %s...",T_S(),fn_db[c].fn);
		#ifdef _WIN32
			sprintf(tdbname, "%s\\%s", dbpath, fn_db[c].fn);
		#else 
			sprintf(tdbname, "%s/%s", dbpath, fn_db[c].fn);
		#endif
		sprintf(cmd, "%s %s", cmd, tdbname);
		printf("done!\n");
	}
	sprintf(cmd, "%s | %s -p F -i stdin ", cmd, formatdb);
	sprintf(cmd, "%s -n %s%d -l %s -o", cmd, fn_vdb, bloop, fn_fmtdblog);
	p_exec(cmd);						/*formatdb command*/
	printf("%s build and format virtual database...done!\n",T_S());
	return PASS;
}

/* put a single database file into a virtual database */
int format_single_database(void) {
	char tdbname[FILENAME];
	char cmd[FCOMMAND];

	printf("%s build and format virtual single database...\n",T_S());
	#ifdef _WIN32
		sprintf(tdbname, "%s\\%s", dbpath, fn_db0.fn);
	#else 
		sprintf(tdbname, "%s/%s", dbpath, fn_db0.fn);
	#endif
	sprintf(cmd, "%s -p F -i %s -n %s", formatdb, tdbname, fn_vdb);
	sprintf(cmd, "%s -l %s -o", cmd, fn_fmtdblog);
	p_exec(cmd);						/*formatdb command*/
	printf("%s build and format virtual single database...done!\n",T_S());
	return PASS;
}

/* put multiple databases into a virtual database */
int format_multiple_database(void) {
	char tdbname[FILENAME];
	char cmd[FCOMMAND];
	int d=0;

	printf("%s build/format multiple virtual database...\n",T_S());
	for(d=0; d < db_count; d++) {
		printf("%s virtual database for %s...\n",T_S(),fn_db[d].fn);
		#ifdef _WIN32
			sprintf(tdbname, "%s\\%s", dbpath, fn_db[d].fn);
		#else 
			sprintf(tdbname, "%s/%s", dbpath, fn_db[d].fn);
		#endif
		sprintf(cmd, "%s -p F -i %s -n %s%d", formatdb, tdbname, fn_vdb, d);
		sprintf(cmd, "%s -l %s -o", cmd, fn_fmtdblog);
		p_exec(cmd);						/*formatdb command*/
		printf("done!\n");
	}
	printf("%s build/format multiple virtual database...done!\n",T_S());
	return PASS;
}

/* Finds all the unique oligos form the primer pairs */
int find_uniq_oligos(void) {	/*find uniq oligos, save into memory*/
	int uo=0, pp=0, flag=0, t_count=0;
	UOligo *t_uoligo;
	
	printf("%s find unique oligos from primer pairs...\n",T_S());
	t_count = 2 * ppair_count;			/*default size*/
	t_uoligo=(UOligo *)p_malloc(t_count*sizeof(UOligo));
	
	/* go through each primer pair */
	for(pp=0, uo=0; pp < ppair_count; pp++) {
		flag = chk_if_exist(ppair[pp].lseq, t_uoligo, uo);
		/*If the left sequence isn't found, add it to the list */
		if(flag < 0) {
			strcpy(t_uoligo[uo].id, ppair[pp].lid);
			p_oligo_strcpy(t_uoligo[uo].nseq, ppair[pp].lseq);
			t_uoligo[uo].len = (int)strlen(t_uoligo[uo].nseq); 
			uo++;
		}
		flag = chk_if_exist(ppair[pp].rseq, t_uoligo, uo);
		/*If the right sequence isn't found, add it to the list */
		if(flag < 0) {
			strcpy(t_uoligo[uo].id, ppair[pp].rid);
			p_oligo_strcpy(t_uoligo[uo].nseq, ppair[pp].rseq);
			t_uoligo[uo].len = (int)strlen(t_uoligo[uo].nseq); 
			uo++;
		}
	}
	uoligo_count = uo;
	if(uoligo_count==0)
		p_exit("there is no unique oligo?\timpossible!\n");
	printf("%s unique oligo count = %d\n",T_S(),uoligo_count);
	
	/* Copy the unique oligos into uoligo */
	uoligo=(UOligo *)p_malloc(uoligo_count*sizeof(UOligo));
	for(uo=0; uo < uoligo_count; uo++) {
		strcpy(uoligo[uo].id, t_uoligo[uo].id);
		strcpy(uoligo[uo].nseq, t_uoligo[uo].nseq);
		uoligo[uo].len = t_uoligo[uo].len;
	}
	free(t_uoligo);
	printf("%s find unique oligos from primer pairs...", T_S());
	printf("saved into memory!\n");
	find_uniq_oligos_test();
	return PASS;
}

/* REMOVE? */
/* Print the unique oligos */
int find_uniq_oligos_test() {
	int uo=0;
	
	printf("%s find unique oligos...\n",T_S());
	for(uo=0; uo < uoligo_count; uo++) {
		printf("%6d) %25s [%s]",uo+1,uoligo[uo].id,uoligo[uo].nseq);
		printf("\tlength %d\n", strlen(uoligo[uo].nseq));
	}
	printf("%s find unique oligos...done!\n",T_S());
	return PASS;
}

/* Look for oligo in tuoligo */
int chk_if_exist(char *oligo, UOligo *tuoligo, int size) {
	int pos=0, res=0;
	
	for(pos=0; pos < size; pos++) {
		if(p_oligo_strcmp(tuoligo[pos].nseq, oligo)==0)
			return pos;
	}
	return -1;
}

/* run blast on the unique oligos from chromosomes */
int run_uo_blast_chromosome(int bloop, int cloop) {
	char cmd[BCOMMAND];
	
	printf("%s run megablast for unique oligos...", T_S());
	sprintf(cmd,"%s -d %s%d -i %s ",megablast,fn_vdb,bloop,fn_bin);
	sprintf(cmd,"%s -W %d -F F -D 2 -o %s",cmd,WORDSIZE,fn_bout);
	p_exec(cmd);						/*blast command*/
	printf("done!\n");
	return PASS;
}


/* run blast on the unique oligos from a single database */
int run_uo_blast_single(void) {
	char cmd[BCOMMAND];
	char cmd1[BCOMMAND];
	
	printf("%s run blast for input unique oligo from ", T_S());
	printf("%s against virtual databse %s...\n", fn_bin, fn_vdb);
	sprintf(cmd1,"%s -d %s -i %s ",megablast, fn_vdb, fn_bin);
	sprintf(cmd,"%s -W %d -F F -D 2 -o %s",cmd1,WORDSIZE,fn_bout);
	p_exec(cmd);						/*blast command*/
	printf("%s megablast for input oligos from %s", T_S(), fn_bin);
	printf(" against virtual databse %s...done!\n",fn_vdb);
	return PASS;
}

/* run blast on oligos from multiple databases */
int run_uo_blast_multiple(int dbcount) {
	char cmd[BCOMMAND];
	
	printf("%s run blast for input unique oligo from ", T_S());
	printf("%s against virtual databse %s...\n", fn_bin, fn_vdb);
	sprintf(cmd,"%s -d %s%d -i %s ",megablast,fn_vdb,dbcount,fn_bin);
	sprintf(cmd,"%s -W %d -F F -D 2 -p 100 -o %s",cmd,WORDSIZE,fn_bout);
	p_exec(cmd);						/*blast command*/
	printf("%s megablast for input oligos from %s", T_S(), fn_bin);
	printf(" against virtual databse %s%d...done!\n",fn_vdb,dbcount);
	return PASS;
}

/* Save the output from blast on the unique oligos taken from chromosomes */
int save_blast_output_chromosome(int bloop, int cloop) {
	FILE *fp;
	char line[MLINE], tmp[SLINE], chr[BISULCHR], *tok;
	int count=0, start_flag=0, res=0;
	int	br_ndx=0, br=0;		/*br_ndx:blast records index*/
	int	uo=0, bo=0, db=0, h=0, len=0;	/*bo:blast output for a uligo*/
	int qs=0, qe=0, ss=0, se=0;		/*query & sbjct*/
	size_t size;

	printf("%s read unique oligos blast hits from", T_S());
	printf(" %s and save blast hits into memory...\n", fn_bout);

	for(uo=0; uo < uoligo_count; uo++) {
		for(db=0; db < CHRDBINLOOP; db++) {
			br_ndx = cloop + db;
			uoligo[uo].brecord[br_ndx].bh = 0;
			strcpy(chr, fn_db[br_ndx].fn);
			if(strcmp(uoligo[uo].brecord[br_ndx].chr, chr)!=0)
				p_exit("problem in saving blast output!\n");
		}
	}	/*important; not all brecords have hits*/
	
	fp=p_fopen_error_exit(fn_bout, "r");
	bo=-1;
	while(fgets(line, MLINE, fp)) {			/*count blast hit*/
		if(strncmp(line, "Query=", 6)==0) {
			bo++;		/*blast output for next oligo*/
			br_ndx = cloop-1;	/*initialize for next uoligo*/
			sscanf(line,"%*s%s",tmp);
			for(uo=0; uo < uoligo_count; uo++) 
				if(strcmp(tmp,uoligo[uo].id)==0)
					break;
			if(bo!=uo) 
				p_continue("uoligo-blast order mismatch\n");
			continue;
		}
		/* found a hit */
		if(strncmp(line, "Query:", 6)==0) {	
			uoligo[uo].brecord[br_ndx].bh++;
			continue;
		}
		/* make sure the blast record chromosome matches the BLAST output file */
		if(strncmp(line, ">", 1)==0) {
			tok = line;
			tok++;
			sscanf(tok, "%s", chr);
			res=FAIL;
			for(br_ndx=cloop; br_ndx < cloop+CHRDBINLOOP; br_ndx++) {
				if(strncmp(chr, uoligo[uo].brecord[br_ndx].chr, strlen(chr))==0){
					res=PASS;
					break;
				}
			}
			if(res==FAIL)
				p_exit("Brecord.chr is not same as in BLAST output file!\n"); 
			continue;
		}
	}
	/* print the number of hits */
	for(uo=0; uo < uoligo_count; uo++) {
		for(db=0; db < CHRDBINLOOP; db++) {
			printf("%s unique oligo %02d) ", T_S(), uo+1);
			br_ndx = cloop + db;
			if(uoligo[uo].brecord[br_ndx].bh <= 0) {
				printf(" %-25s ",uoligo[uo].id);
				printf("%-30s ",uoligo[uo].nseq);
				printf("blast record#%3d ... ", br_ndx+1);
				printf("       no hits found!\n");
			}
			else {
				printf(" %-25s ",uoligo[uo].id);
				printf("%-30s ",uoligo[uo].nseq);
				printf("blast record#%3d ... ", br_ndx+1);
				printf("%-6d hits found!\n",uoligo[uo].brecord[br_ndx].bh);
				size = uoligo[uo].brecord[br_ndx].bh*sizeof(int);
				uoligo[uo].brecord[br_ndx].bht=(int *)p_malloc(size);
				uoligo[uo].brecord[br_ndx].bhs=(int *)p_malloc(size);
				uoligo[uo].brecord[br_ndx].bhe=(int *)p_malloc(size);
			}
		}
	}
	
	rewind(fp);
	bo=-1;
	
	/* now store the hits themselves */
	while(fgets(line, MLINE, fp)) {
		if(strncmp(line, "Query=", 6)==0) {
			bo++;				/*blast oligo ndx*/
			br_ndx = cloop -1;		/*br_ndx start for next uoligo*/
			sscanf(line,"%*s%s",tmp);
			for(uo=0; uo < uoligo_count; uo++) 
				if(strcmp(tmp,uoligo[uo].id)==0)
					break;
			if(bo!=uo) 
				p_continue("blast-oligo order mismatch\n");
			continue;
		}
		if(strncmp(line, " Strand = ", 10)==0) {
			++h;				/*new hit*/
			sscanf(line,"%*s%*s%*s%*s%s",tmp);
			if(strcmp(tmp, TPP)==0) 
				uoligo[uo].brecord[br_ndx].bht[h] = PP;
			else if(strcmp(tmp, TPM)==0) 
				uoligo[uo].brecord[br_ndx].bht[h] = PM;
			fgets(line, MLINE, fp);		/*junk*/
			fgets(line, MLINE, fp);		/*junk*/
			fgets(line, MLINE, fp);	
			if(strncmp(line, "Query:", 6)!=0) 
				p_exit("blast format error!\n");
			sscanf(line,"%*s%d%*s%d", &qs, &qe);
			fgets(line, MLINE, fp);		/*junk*/
			fgets(line, MLINE, fp);	
			/* this line should be the hit info */
			if(strncmp(line, "Sbjct:", 6)!=0) 
				p_exit("blast format error!\n");
			sscanf(line,"%*s%d%*s%d", &ss, &se);
			if(uoligo[uo].brecord[br_ndx].bht[h]==PP) {
				len = uoligo[uo].len;
				uoligo[uo].brecord[br_ndx].bhs[h]=ss-(qs-1);	/*start*/
				uoligo[uo].brecord[br_ndx].bhe[h]=se+(len-qe);	/*end*/
			}
			if(uoligo[uo].brecord[br_ndx].bht[h]==PM) {
				len = uoligo[uo].len;
				uoligo[uo].brecord[br_ndx].bhs[h]=se-(len-qe);	/*start*/
				uoligo[uo].brecord[br_ndx].bhe[h]=ss+(qs-1);	/*end*/
			}
			continue;
		}
		/* make sure the blast record chromosome matches the BLAST output file */
		if(strncmp(line, ">", 1)==0) {	
			tok = line;
			tok++;
			sscanf(tok, "%s", chr);
			res=FAIL;
			for(br_ndx=cloop; br_ndx < cloop+CHRDBINLOOP; br_ndx++) {
				if(strncmp(chr, uoligo[uo].brecord[br_ndx].chr, strlen(chr))==0){
					res=PASS;
					break;
				}
			}
			if(res==FAIL)
				p_exit("Brecord.chr is not same as in BLAST output file!\n"); 
			h = -1;				/*bhast hit ndx*/
			continue;
		}
	}
	fclose(fp);
	printf("%s unique oligos blast hit from %s for",T_S(), fn_bout);
	printf(" bloop %02d, cloop %02d\tsaved into memory!\n",bloop,cloop);
	res = PASS;
	return res;
}

/* Save the output from blast on the unique oligos taken from a single database */
int save_blast_output_single(void) {
	FILE *fp;
	char line[MLINE], tmp[SLINE];
	int count=0, start_flag=0, res=0;
	int	br_ndx=0;			/*br_ndx:blast records index*/
	int	uo=0, bo=0, db=0, h=0, len=0;	/*bo:blast output for a uligo*/
	int qs=0, qe=0, ss=0, se=0;			/*query & sbjct*/
	size_t size;

	printf("%s read unique oligos blast hits from", T_S());
	printf(" %s and save blast hits into memory...\n", fn_bout);
	
	fp=p_fopen_error_exit(fn_bout, "r");
	bo=-1;
	while(fgets(line, MLINE, fp)) {			/*count blast hit*/
		if(strncmp(line, "Query=", 6)==0) {
			bo++;		/*blast output for next oligo*/
			sscanf(line,"%*s%s",tmp);
			for(uo=0; uo < uoligo_count; uo++) 
				if(strcmp(tmp,uoligo[uo].id)==0)
					break;
			if(bo!=uo) 
				p_continue("uoligo-blast order mismatch\n");
			uoligo[uo].br_count = 0;
			continue;
		}
		/* count the hits for the current uligo */
		if(strncmp(line,BHITTAG,BHITCOMP)==0) {
			br_ndx = 0;						/*blast hit counter for current uoligo*/
			fgets(line, MLINE, fp);		/*empty line*/
			fgets(line, MLINE, fp);		/*first hit always exist*/
			do {
				br_ndx++;
				fgets(line, MLINE, fp);
			} while(strlen(line) >= 2);		
			size = br_ndx*sizeof(BRecord);
			uoligo[uo].br_count = br_ndx;
			uoligo[uo].brecord=(BRecord *)p_malloc(size);
			for(br_ndx=0; br_ndx < uoligo[uo].br_count; br_ndx++)
				uoligo[uo].brecord[br_ndx].bh=0;
			br_ndx = -1;
			continue;
		}
		if(strncmp(line, ">", 1)==0) {	
			br_ndx++;
			sscanf(line, ">%s", uoligo[uo].brecord[br_ndx].chr);
			continue;
		}
		/* found a hit */
		if(strncmp(line, "Query:", 6)==0) {	
			uoligo[uo].brecord[br_ndx].bh++;
			continue;
		}
	}
	
	/* print the hits */
	for(uo=0; uo < uoligo_count; uo++) {
		printf("%s unique oligo %02d) ", T_S(), uo+1);
		printf(" %-25s ",uoligo[uo].id);
		printf("%-30s ",uoligo[uo].nseq);
		if(uoligo[uo].br_count==0)
			printf("       no hits found!\n");
		else {
			printf("%-6d hits found!\n",uoligo[uo].br_count);
			for(br_ndx=0; br_ndx < uoligo[uo].br_count; br_ndx++) {
				if (dlevel > 3) {
					printf("%s\t %4d) ", T_S(), br_ndx+1);
					printf("%-30s",uoligo[uo].brecord[br_ndx].chr);
					printf("...%2d hits\n",uoligo[uo].brecord[br_ndx].bh);
				}
				size = uoligo[uo].brecord[br_ndx].bh*sizeof(int);
				uoligo[uo].brecord[br_ndx].bht=(int *)p_malloc(size);
				uoligo[uo].brecord[br_ndx].bhs=(int *)p_malloc(size);
				uoligo[uo].brecord[br_ndx].bhe=(int *)p_malloc(size);
			}
		}
	}

	rewind(fp);
	bo=-1;
	/* store the hits themselves */
	while(fgets(line, MLINE, fp)) {
		if(strncmp(line, "Query=", 6)==0) {
			bo++;				/*blast oligo ndx*/
			br_ndx = -1;
			sscanf(line,"%*s%s",tmp);
			for(uo=0; uo < uoligo_count; uo++) 
				if(strcmp(tmp,uoligo[uo].id)==0)
					break;
			if(bo!=uo) 
				p_continue("blast-oligo order mismatch\n");
			continue;
		}
		if(strncmp(line, ">", 1)==0) {	
			br_ndx++;
			h = -1;				/*bhast hit ndx*/
			continue;
		}
		if(strncmp(line, " Strand = ", 10)==0) {
			++h;				/*new hit*/
			sscanf(line,"%*s%*s%*s%*s%s",tmp);
			if(strcmp(tmp, TPP)==0) 
				uoligo[uo].brecord[br_ndx].bht[h] = PP;
			else if(strcmp(tmp, TPM)==0) 
				uoligo[uo].brecord[br_ndx].bht[h] = PM;
			fgets(line, MLINE, fp);		/*junk*/
			fgets(line, MLINE, fp);		/*junk*/
			fgets(line, MLINE, fp);	
			if(strncmp(line, "Query:", 6)!=0) 
				p_exit("blast format error!\n");
			sscanf(line,"%*s%d%*s%d", &qs, &qe);
			fgets(line, MLINE, fp);		/*junk*/
			fgets(line, MLINE, fp);	
			if(strncmp(line, "Sbjct:", 6)!=0) 
				p_exit("blast format error!\n");
			sscanf(line,"%*s%d%*s%d", &ss, &se);
                        if(uoligo[uo].brecord[br_ndx].bht[h]==PP) {
				len = uoligo[uo].len;
				uoligo[uo].brecord[br_ndx].bhs[h]=ss-(qs-1);	/*start*/
				uoligo[uo].brecord[br_ndx].bhe[h]=se+(len-qe);	/*end*/
				if (dlevel > 3) {
					printf("%d ",uoligo[uo].brecord[br_ndx].bhs[h]);
					printf("%d\n",uoligo[uo].brecord[br_ndx].bhe[h]);
				}
			}
			if(uoligo[uo].brecord[br_ndx].bht[h]==PM) {
				len = uoligo[uo].len;
				uoligo[uo].brecord[br_ndx].bhs[h]=se-(len-qe);	/*start*/
				uoligo[uo].brecord[br_ndx].bhe[h]=ss+(qs-1);	/*end*/
				if (dlevel > 3) {
					printf("%d ",uoligo[uo].brecord[br_ndx].bhs[h]);
					printf("%d\n",uoligo[uo].brecord[br_ndx].bhe[h]);
				}
			}
			continue;
		}
	}
	fclose(fp);
	printf("%s unique oligos blast hit from %s",T_S(), fn_bout);
	printf("\tsaved into memory!\n");
	return PASS;
}
