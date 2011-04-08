/*
 * this module contains the Probe Specific Primer Design function.
 * 
 */

/* 
 * Copyright (c) 2002, Digital Bilogy Laboratory, University of 
 * Missouri-Columbia. All rights reserved. Please see full software
 * agreement in Copyright.txt or type primegens -h.
 */


/*
 * Sequence-specific primer design
 * PRIMEGENS execution:
 * $primegens.exe -q <filename> -d <database> -p <pathname>
 * $primegens.exe -q <filename> -l <database> -p <pathname>
 * -q: query sequence input filename
 * -l: multiple database filename list filename
 * -d: single database filename
 * -p: database location pathname 
 *
 */

/* header files*/

#include "lib.h"
#include "defs.h"
#include "util.h"
#include "pspd.h"
#include "hbdzn.h"
#include "output.h"
#include "pd.h"
#include "pphit.h"
#include "SeqMapping.h"
#include "uhit.h"


/* use the pspd algorithm on the queries */
short pspd_loop (void) { 
	
	/*Local Variable for the function */					
	register int count, sdx, sdx1, sdx2;
	register int ndx4;
	int q=0, b=0, db=0, uo=0, pp=0, h=0, c=0;
	int res=0, bloop=0, cloop=0, flag=0;
	int *length, ndx1, ndx2, ndx3, nsim;
	int best_score, best_sim, best_ndx, tmp_score, sim;
	int start_flag, read, start, end, start_tmp;
	int max, resid_shift=0, identical, size;
	int *seqmap, *inversemap, left, right;
	short block[PSPD_NumBlastMatch][2];
	short **similarity;
	size_t t_size;
	char string[4][QNSIZE], *tok;
	char qseq_fn[FILENAME], qseq_blast_fn[FILENAME];
	char tqfn[FILENAME];		/*temporary filename*/
	char cmd[COMMAND];
	char Expect[CONFIGFORMAT], line[BLINE];
	char **seq;					/*seq[i]=ith homologue nseq*/ 
	char tmp0[FILENAME];
	FILE *fp, *fblast, *fqal, *funiseg, *fsim, *fnohit, *fout4;
	SeqID *Hitlist;

	printf("%s selecting probe-specific primer design algorithm...\n",T_S());

	res = read_query(fn_query);			/*saved in *query*/
	if(res==FAIL)
		p_exit("problem inside read_query function\n");

	sprintf(cmd, "mkdir %s", tmp_dir);		/*working direcotry*/
	p_exec(cmd);
	sprintf(alignDir, "%s%s", out_path, ALIGNMENTDIR);
	sprintf(cmd, "mkdir %s", alignDir);	
	p_exec(cmd);
	strcpy(Expect, PSPD_EXPECT);
	sprintf(qseq_fn, "%s%s", tmp_dir, QSEQFILENAME);
	sprintf(qseq_blast_fn, "%s%s", tmp_dir, QSEQBLASTFILENAME);

	/*formatting database based on input type*/
	if(DATABASE_TYPE==3) {
		for(bloop=0, cloop=0; bloop < BLASTLOOP; bloop++) {
			res = format_chromosome_database(bloop, cloop);
			if(res==FAIL)
				p_exit("problem while formatdb all databases\n");
			cloop = cloop + CHRDBINLOOP;	
		}
	} else if(DATABASE_TYPE==2) {
			res = format_multiple_database(); 
			if(res==FAIL)
				p_exit("problem while formatdb all databases\n");
	} else if(DATABASE_TYPE==1) {
			res = format_single_database(); 
			if(res==FAIL)
				p_exit("problem while formatdb all databases\n");
	} else 
		p_exit("problem in selecting database type!\n");

	fnohit = p_fopen_error_exit(fn_nohit, "w");
	if(PSPD_RUN_SIM_FLAG==0)
		fsim = p_fopen_error_exit(fn_sim, "w");
	funiseg = p_fopen_error_exit(fn_uniseg, "w");
	fout4 = p_fopen_error_exit(fn_out4, "w");   /*record failed queries*/
   
	
	for(q=0; q < q_count; q++) {
		printf("\n\n");
		printf("%s =====================================",T_S());
		printf("==========================================\n");
		printf("%s %05d) processing %.60s\n",T_S(),q+1,query[q].id);
		printf("%s =====================================",T_S());
		printf("==========================================\n");
        /*skip if no sequence available*/
        if((query[q].len < 10)||(strlen(query[q].nseq) < 10)) {
            printf("%s query has no DNA sequence\tskipped.\n", T_S());  
            continue;
        }
		
		#ifdef _WIN32
			sprintf(cmd, "del %s", qseq_fn);
		#else
			sprintf(cmd, "rm %s", qseq_fn);
		#endif
		p_exec(cmd);
		fp = p_fopen_error_exit(qseq_fn, "w");
		fprintf(fp, "%s\n", query[q].id);
		fprintf(fp, "%s\n", query[q].nseq);
		fclose(fp);
		#ifdef _WIN32
			sprintf(cmd, "%s -p blastn -i %s -d %s -F F -G %d -E %d \
	   					-e %s -o %s", blastall, qseq_fn, fn_vdb, 
	      				(BLAST_OGP + BLAST_OIP), 
	      				BLAST_OIP, Expect, qseq_blast_fn);
		#else 
			sprintf(cmd, "%s -p blastn -i %s -d %s -F F -G %d -E %d \
	   					-e %s -o %s", blastall, qseq_fn, fn_vdb, 
	      				(BLAST_OGP + BLAST_OIP), 
	      				BLAST_OIP, Expect, qseq_blast_fn);
		#endif 
		p_exec(cmd);
		/*
		 * Read Blast result file and count the number of Hits for the 
		 * current sequence and store in flag                       
		 */
		fblast = p_fopen_error_exit(qseq_blast_fn, "r");
		flag = -1;
		while(fgets(line, BLINE, fblast)) {
			if((flag > 0) && (line[0]=='\n')) 
				break; 
			else if(strncmp(line,BHITTAG,BHITCOMP)==0)
				flag = 0;
			else if((flag >=0) && (line[0] !='\n')) {
				flag++;
			}
		}
		printf("%s total query blast hit = %d\n", T_S(), flag);

		
		if(flag >= 1) {				
			Hitlist=(SeqID *)p_malloc(flag*sizeof(SeqID));
			rewind(fblast); 
			flag = -1;
			while(fgets(line, BLINE, fblast)) {
				if((flag > 0) && (line[0]=='\n')) 
					break; 
				else if(strncmp(line,BHITTAG,BHITCOMP)==0)
					flag = 0;
				else if((flag >=0) && (line[0] !='\n')) {
					tok = strtok(line, " .\n");
					strcat(tok, "\0");
					sprintf(Hitlist[flag].id, ">%s", tok);
					flag++;
				}
			}
		}
		fclose(fblast);
		if(flag < MAXHOMOLOGS) {
			for(h=0; h < flag; h++) {
				if(dlevel >=1)
					printf("%s %03d) %s\n",T_S(),h+1,Hitlist[h].id);
			}
		} else 
			printf("%s query has %d homologs\thidden!\n",T_S(),flag); 
		
		if(flag <= 1) {				
			fprintf(fnohit, "%s\n", query[q].id);
			for(b=0; b < query[q].len; b++) {	
				fprintf(fnohit, "%c", query[q].nseq[b]);	
				if((b + 1) % FASTALINE == 0) 
					fprintf(fnohit, "\n"); 
			}
			if(query[q].len % FASTALINE != 0) 
				fprintf(fnohit, "\n");	
			printf("%s adding query into %s...done!\n",T_S(),fn_nohit); 
			
			/*primer design with hybridization check*/
			
		/*STEP 1) primer design using primer3 ]**************/
			printf("%s\t\t STEP 1: primer design using primer3",T_S());
			printf(" for query # %06d [%.30s]...\n", q+1, query[q].id);
			run_primer3(query[q]);
			res = save_primer3_output();        /*returns pass-fail*/
			if((res==FAIL)||(ppair_count<=0)) {
				query[q].res = FAIL;
				if((QUERY_FORMAT==1)||(QUERY_FORMAT==2)) {
					fprintf(fout4, "%s\n", query[q].id);
					if(QUERY_FORMAT==2)
						fprintf(fout4, "%s\n", qinfo[q].desc);
					else
						fprintf(fout4, "\n");
				} else {
					p_exit("query format chromosome is not supported!\n");	
				}
				if(QUERY_SEQUENCE==1) {
					for(b=0; b < query[q].len; b++) {
						fprintf(fout4, "%c", query[q].nseq[b]);
						if((b + 1) % FASTALINE == 0)
							fprintf(fout4, "\n");
					}
					if(query[q].len % FASTALINE != 0)
						fprintf(fout4, "\n");
				}
				continue;
			} else
				query[q].res = PASS;        /*primer design successful*/
			sprintf(cmd, "rm %s %s", fn_p3in, fn_p3out);
			p_exec(cmd);                    /*delete primer3 tmp files*/
			
		/*STEP 2) find unique oligos and save blast hits ]********/
			printf("%s\t\tSTEP 2: find uoilgos blast records",T_S());
			printf(" for query # %06d [%.30s]...\n", q+1, query[q].id);
			res = find_uniq_oligos();
			if(res==FAIL)
				p_exit("some problem inside find_uniq_oligos\n");
			printf("%s write uoligos into %s file...\n",T_S(),fn_bin);
			fp = p_fopen_error_exit(fn_bin, "w");
			for(uo=0; uo < uoligo_count; uo++) {
				fprintf(fp,">%s\n", uoligo[uo].id);
				fprintf(fp,"%s\n", uoligo[uo].nseq);
			}
			fclose(fp);
			printf("%s uoligos into %s file...saved!\n",T_S(),fn_bin);
			if(DATABASE_TYPE==1) {
				run_uo_blast_single();
				save_blast_output_single();
				sprintf(cmd, "rm %s", fn_bout);
				p_exec(cmd);
			} else
				p_exit("Other database type is not supported in PSPD!\n");
			sprintf(cmd, "rm %s", fn_bin);
			p_exec(cmd);
			
		/*STEP 3) build primer pairs chit records ]**************/
			printf("%s\t\t STEP 3: find primer pair blast records",T_S());
			printf(" for query # %06d [%.30s]...\n", q+1, query[q].id);
			if(DATABASE_TYPE==1) {
				res = ppair_chit_record_single();
				if(res==FAIL)
					p_exit("problem in saving blast records for pairs\n");
				printf("%s free uoligos blast records...\n", T_S());
				for(uo=0; uo < uoligo_count; uo++) {
					if(uoligo[uo].br_count > 0) {
						for(b=0; b < uoligo[uo].br_count; b++) {
							free(uoligo[uo].brecord[b].bht);
							free(uoligo[uo].brecord[b].bhs);
							free(uoligo[uo].brecord[b].bhe);
						}
						free(uoligo[uo].brecord);
					}
				}
				free(uoligo);
			} else {
				p_exit("only single database type is supported for PSPD\n");
			}
			printf("%s free uoligos blast records...done!\n", T_S());

		/*STEP 4) build potential hybridization records ]********/
			printf("%s\t\t STEP 4: find potential hybridization",T_S());
			printf(" for query # %06d [%.30s]...\n", q+1, query[q].id);
			if(DATABASE_TYPE==1) {
				res = chk_hybridization_single(q);      /*qinfo is required*/
				if(res==FAIL)
					p_exit("problem inside chk_hybirization function\n");
				printf("%s free ppair chit records...\n", T_S()); 
				for(pp=0; pp < ppair_count; pp++) { 
					for(c=0; c < ppair[pp].chit_count; c++) {
						if(ppair[pp].chit[c].ppc > 0)
							free(ppair[pp].chit[c].pps);
						if(ppair[pp].chit[c].pmc > 0)
							free(ppair[pp].chit[c].pme);
					}
					if(ppair[pp].chit_count > 0)
						free(ppair[pp].chit);
				}
			} else {
				p_exit("only single database type is supported for PSPD\n");
			}
			printf("%s free ppair chit records...done!\n", T_S()); 
		
		/*STEP 5) rank all primer pairs and save into file ]*****/
			printf("%s\t\t STEP 5: primer ranking and save",T_S());
			printf(" for query # %06d [%.30s]...\n", q+1, query[q].id);
			res = rank_ppair();             /*rank primer pairs*/
			if(res==FAIL)
				p_exit("problem inside primer pair ranking\n");
			res = print_output(q);          /*save output*/
			if(res==FAIL)
				p_exit("problem inside print output\n");
			for(pp=0; pp < ppair_count; pp++) {
				if(ppair[pp].hbrdn_count > 0) {
					free(ppair[pp].hbrdn);
				}
			}
			free(ppair);
			printf("\n\n");
			printf("%s %05d) query # %2d [%.30s]",T_S(),q+1,q+1,query[q].id);
			printf("...processed\t\tHave a great time!\n\n");
				
		} else {				
			sprintf(tqfn,"%s%s%d%s", alignDir, qinfo[q].id1, q+1, HOMOA);
			fqal = p_fopen_error_exit(tqfn, "w");	   		
			fprintf(fqal, "\n%s\t\tGLOBAL ALIGNMENTS\n",T_S());
			for(h=0; h < flag; h++) {
				uo = p_min(BLASTSEQID, strlen(query[q].id));
				if(strncmp(Hitlist[h].id, query[q].id, uo)==0){
					if(h != 0) {
						strcpy (line, Hitlist[h].id);
						strcpy (Hitlist[h].id, Hitlist[0].id);
						strcpy (Hitlist[0].id, line);          
					}
					break;
				}
			}
			length = (int *)p_malloc(flag*sizeof(int));  
			seq = (char **)p_malloc(flag*sizeof(char *)); 
			for(h=0; h < flag; h++) {	 
				seq[h] = (char *)p_malloc(QNSIZE*sizeof(char)); 
				get_nseq_from_db(Hitlist[h].id, seq[h], &length[h]);
				fprintf(fqal,"%s (length=%d)\n",Hitlist[h].id,length[h]);
				fprintf(fqal, "%s\n", seq[h]);
			}	
			if(dlevel >=2)
				printf("%s read blast to get aligned blocks...\n",T_S());
	   		//read blast file again to find aligned blocks
			fblast = p_fopen_error_exit(qseq_blast_fn, "r");
			read = -2; //initialization 
			ndx1 = 1;  	//initialization
	   		while(fgets(line, BLINE, fblast)) {
				if(strncmp(line, ">", 1)==0) {				
						start_flag = 0;
						if (read >= 0) {	
							/* aligned portion of query sequence */
							if(end - start >= PSPD_MINFILTERLENGTH) {
								block[read][0] = start; 
								block[read][1] = end;
							} else
								read--;
						}
						read++;	
						tok = strtok(line, "\n");
						strcat(tok, "\0");
						strcpy(tmp0, tok);	
						uo = p_min(BLASTSEQID, strlen(tmp0));
				}	
				if(read>=0 && strncmp(Hitlist[0].id, tmp0, uo)!=0) {	
	         		if((strncmp(Hitlist[ndx1].id, tmp0, uo)!=0) && 
	            			(strncmp(Hitlist[ndx1+1].id, tmp0, uo)==0))	
						ndx1++;
					if(strncmp (line,"Query:",6)==0) {	
						sscanf(line,"%*s%d%*s%d",&start_tmp,&end);
						if(start_flag == 0) {		
							start_flag = 1;
							start = start_tmp;
						}
					}
				}
			}
			fclose(fblast);
			if(dlevel >1)
				printf("%s all aligned blocks saved!\n",T_S());
			if((read >= 0) && (end-start >= PSPD_MINFILTERLENGTH)) {
					block[read][0]=start;
					block[read][1] = end;
			}
			seqmap=(int *)p_malloc((length[0]+1)*sizeof(int)); 
			similarity=(short **)p_malloc(length[0]*sizeof(short *));
			for(ndx2=0; ndx2 < length[0]; ndx2++) {
				t_size = length[0]*sizeof(short);
				similarity[ndx2]=(short *)p_malloc(t_size);
				for (ndx3=0; ndx3 < length[0]; ndx3++)		
					similarity[ndx2][ndx3] = 0;
			}
			best_score = LARGESCORE;
			for(ndx1=1; ndx1 < flag; ndx1++) {	//for each homolog
				if(dlevel >1)
					printf("%s query-homolog(%d) alignment\n",T_S(),ndx1);
				for(ndx2=0; ndx2 <= length[0]; ndx2++)		
					seqmap[ndx2] = EMPTY;
				t_size = (length[ndx1]+1)*sizeof (int);
				inversemap = (int *)p_malloc(t_size); 
				for(ndx2=0; ndx2 <= length[ndx1]; ndx2++)	
					inversemap[ndx2] = EMPTY;
				tmp_score = SeqMapping(seq[ndx1],seq[0],seqmap, 
			               inversemap,length[ndx1],length[0],&sim);
				if(best_score > tmp_score) {
						best_score = tmp_score; 
						best_sim = sim;
						best_ndx = ndx1;
				}
				for(sdx=0;sdx < QNSIZE; sdx++) {
					string[0][sdx] = ' ';
					string[1][sdx] = ' ';
					string[2][sdx] = ' ';
					string[3][sdx] = ' ';
				}
				resid_shift = 0;
				identical = 0;
				for(sdx=0;sdx < length[0];sdx++) { 
					sdx1 = sdx + resid_shift;
					if (seqmap[sdx] != EMPTY) {
						string[3][sdx1] = seq[0][sdx]; 
						string[2][sdx1]=seq[ndx1][seqmap[sdx]];
						if(string[3][sdx1]==string[2][sdx1]) {
							string[1][sdx1] = '|';
							identical++;
						} else 
							string[1][sdx1] = ' ';
						sdx2 = seqmap[sdx]+1;
						while (sdx2 < length[ndx1] && \
									inversemap[sdx2] == EMPTY) {
							sdx1++;
							string[3][sdx1] = '-';
							string[2][sdx1] = seq[ndx1][sdx2];
							string[1][sdx1] = ' ';
							sdx2++;
							resid_shift++;
						}
					} else {
						string[3][sdx1] = seq[0][sdx];
						string[2][sdx1] = '-';
						string[1][sdx1] = ' ';
					}
				}
				size = length[0] + resid_shift;
				string[1][size] = '\0';
				string[2][size] = '\0';
				string[3][size] = '\0';
				fprintf(fqal, "\nGlobal Alignment between %s and %s:\n", 
							Hitlist[0].id, Hitlist[ndx1].id);
				fprintf(fqal, "\nSequence identity: %.2f(%d / %d)\n\n", 
							(float) identical*100.0/(float)length[0], 
							identical, length[0]);
				if(dlevel >2) {
					printf("\nGlobal Alignment between %s and %s:\n", 
							Hitlist[0].id, Hitlist[ndx1].id);
					printf("\nSequence identity: %.2f(%d / %d)\n\n", 
							(float) identical*100.0/(float)length[0], 
							identical, length[0]);
				}
				for(sdx=0; sdx < size; sdx += LINEWIDTH) { 
					fprintf(fqal,"%.60s\n",string[3]+sdx);
					fprintf(fqal,"%.60s\n",string[1]+sdx);
					fprintf(fqal,"%.60s\n",string[2]+sdx);
					fprintf(fqal, "\n");
				} 
				if(dlevel >2) {
					for(sdx=0; sdx < size; sdx += LINEWIDTH) { 
						printf("%.60s\n",string[3]+sdx);
						printf("%.60s\n",string[1]+sdx);
						printf("%.60s\n",string[2]+sdx);
						printf("\n");
					}
				}	
				ConstructSim(seq[ndx1], seq[0], seqmap,
								inversemap, length[0], similarity);
				free (inversemap);
			}	//for each homolog
			
			max = 0; 
			left = EMPTY; 
			right = EMPTY;
			for(ndx1=0; ndx1 < length[0]; ndx1++) {
				for(ndx2=ndx1+PSPD_MINSEGLENGTH; ndx2 < length[0]; ndx2++) {
					count = 0;
					for(ndx3 = 0; ndx3 <=read; ndx3++) {
						if(block[ndx3][0]>ndx2 || block[ndx3][1]<ndx1) 
							++count;
					}
					if(read < 0 || count==read+1) {
						count = 0;
						if(similarity[ndx1][ndx2]==flag-1 && ndx2-ndx1 > max) {
							left = ndx1; 
							right = ndx2; 
							max = ndx2 - ndx1;
						}
					}
				}
			}
			for(ndx2=0; ndx2 < length[0]; ndx2++)
				free(similarity[ndx2]);
			free (similarity);
			
			if(left==EMPTY) {
				max = 0; 
				left = EMPTY; 
				right = EMPTY; 
				ndx4 = flag - 1;
				while(ndx4 > 1 && left == EMPTY) {
					for(ndx1=0; ndx1 < length[0]; ndx1++) {
						for(h=ndx1+PSPD_MINSEGLENGTH;h<length[0];h++){
							count = 0;
							for(b=0; b <= read; b++) {
								if(block[b][0]>h||block[b][1]<ndx1) 
									++count;
							}
							if(read < 0 || count==read+1) {
								count = 0;
								if(similarity[ndx1][h]==ndx4 
											&& h-ndx1 > max) {
									left = ndx1; 
									right = h; 
									max = h - ndx1;
								}
							}
						}
					}
					ndx4--;
				}
				nsim = flag - ndx4 - 2;
				if(left == EMPTY) { 
					left = 0; 
					right = length[0]; 
				}
				fprintf(funiseg,"%s ",Hitlist[0].id); 
				fprintf(funiseg,"(length = %d, ", max+1); 
				fprintf(funiseg,"left=%d, right=%d,",left+1,right+1);
				fprintf(funiseg," original = %d)\n", length[0]);
				for(ndx1 = left; ndx1 < right; ndx1++) {
					fprintf(funiseg, "%c", seq[0][ndx1]);
					if((ndx1-left+1) % FASTALINE==0) 
						fprintf(funiseg, "\n");
				}
				if((ndx1-right+1) % FASTALINE !=0) 
					fprintf(funiseg, "\n");
			} else {
				fprintf(funiseg,"%s (length = %d, ", Hitlist[0].id, max+1); 
				fprintf(funiseg, "left = %d, right = %d,",left+1, right+1);
				fprintf(funiseg, " original = %d)\n", length[0]);
				for(ndx1 = left; ndx1 < right; ndx1++) {	
					fprintf(funiseg, "%c", seq[0][ndx1]);
					if((ndx1-left+1) % FASTALINE ==0) 
							fprintf(funiseg, "\n");
				}
				if((ndx1-right+1) % FASTALINE !=0)
					fprintf(funiseg, "\n");
			}
			printf("%s adding query into %s...done!\n",T_S(),fn_uniseg);	
			
			/*primer design with hybridization check*/
			
		/*STEP 1) primer design using primer3 ]**************/
			printf("%s\t\t STEP 1: primer design using primer3",T_S());
			printf(" for query # %06d [%.30s]...\n", q+1, query[q].id);
			run_primer3(query[q]);
			res = save_primer3_output();        /*returns pass-fail*/
			if((res==FAIL)||(ppair_count<=0)) {
				query[q].res = FAIL;
				if((QUERY_FORMAT==1)||(QUERY_FORMAT==2)) {
					fprintf(fout4, "%s\n", query[q].id);
					if(QUERY_FORMAT==2)
						fprintf(fout4, "%s\n", qinfo[q].desc);
					else
						fprintf(fout4, "\n");
				} else {
					p_exit("query format chromosome is not supported!\n");	
				}
				if(QUERY_SEQUENCE==1) {
					for(b=0; b < query[q].len; b++) {
						fprintf(fout4, "%c", query[q].nseq[b]);
						if((b + 1) % FASTALINE == 0)
							fprintf(fout4, "\n");
					}
					if(query[q].len % FASTALINE != 0)
						fprintf(fout4, "\n");
				}
				continue;
			} else
				query[q].res = PASS;        /*primer design successful*/
			sprintf(cmd, "rm %s %s", fn_p3in, fn_p3out);
			p_exec(cmd);                    /*delete primer3 tmp files*/
			
		/*STEP 2) find unique oligos and save blast hits ]********/
			printf("%s\t\tSTEP 2: find uoilgos blast records",T_S());
			printf(" for query # %06d [%.30s]...\n", q+1, query[q].id);
			res = find_uniq_oligos();
			if(res==FAIL)
				p_exit("some problem inside find_uniq_oligos\n");
			printf("%s write uoligos into %s file...\n",T_S(),fn_bin);
			fp = p_fopen_error_exit(fn_bin, "w");
			for(uo=0; uo < uoligo_count; uo++) {
				fprintf(fp,">%s\n", uoligo[uo].id);
				fprintf(fp,"%s\n", uoligo[uo].nseq);
			}
			fclose(fp);
			printf("%s uoligos into %s file...saved!\n",T_S(),fn_bin);
			if(DATABASE_TYPE==1) {
				run_uo_blast_single();
				save_blast_output_single();
				sprintf(cmd, "rm %s", fn_bout);
				p_exec(cmd);
			} else
				p_exit("Other database type is not supported in PSPD!\n");
			sprintf(cmd, "rm %s", fn_bin);
			p_exec(cmd);
			
		/*STEP 3) build primer pairs chit records ]**************/
			printf("%s\t\t STEP 3: find primer pair blast records",T_S());
			printf(" for query # %06d [%.30s]...\n", q+1, query[q].id);
			if(DATABASE_TYPE==1) {
				res = ppair_chit_record_single();
				if(res==FAIL)
					p_exit("problem in saving blast records for pairs\n");
				printf("%s free uoligos blast records...\n", T_S());
				for(uo=0; uo < uoligo_count; uo++) {
					if(uoligo[uo].br_count > 0) {
						for(b=0; b < uoligo[uo].br_count; b++) {
							free(uoligo[uo].brecord[b].bht);
							free(uoligo[uo].brecord[b].bhs);
							free(uoligo[uo].brecord[b].bhe);
						}
						free(uoligo[uo].brecord);
					}
				}
				free(uoligo);
			} else {
				p_exit("only single database type is supported for PSPD\n");
			}
			printf("%s free uoligos blast records...done!\n", T_S());

		/*STEP 4) build potential hybridization records ]********/
			printf("%s\t\t STEP 4: find potential hybridization",T_S());
			printf(" for query # %06d [%.30s]...\n", q+1, query[q].id);
			if(DATABASE_TYPE==1) {
				res = chk_hybridization_single(q);      /*qinfo is required*/
				if(res==FAIL)
					p_exit("problem inside chk_hybirization function\n");
				printf("%s free ppair chit records...\n", T_S()); 
				for(pp=0; pp < ppair_count; pp++) { 
					for(c=0; c < ppair[pp].chit_count; c++) {
						if(ppair[pp].chit[c].ppc > 0)
							free(ppair[pp].chit[c].pps);
						if(ppair[pp].chit[c].pmc > 0)
							free(ppair[pp].chit[c].pme);
					}
					if(ppair[pp].chit_count > 0)
						free(ppair[pp].chit);
				}
			} else {
				p_exit("only single database type is supported for PSPD\n");
			}
			printf("%s free ppair chit records...done!\n", T_S()); 
		
		/*STEP 5) rank all primer pairs and save into file ]*****/
			printf("%s\t\t STEP 5: primer ranking and save",T_S());
			printf(" for query # %06d [%.30s]...\n", q+1, query[q].id);
			res = rank_ppair();             /*rank primer pairs*/
			if(res==FAIL)
				p_exit("problem inside primer pair ranking\n");
			res = print_output(q);          /*save output*/
			if(res==FAIL)
				p_exit("problem inside print output\n");
			for(pp=0; pp < ppair_count; pp++) {
				if(ppair[pp].hbrdn_count > 0) {
					free(ppair[pp].hbrdn);
				}
			}
			free(ppair);
			printf("\n\n");
			printf("%s %05d) query # %2d [%.30s]",T_S(),q+1,q+1,query[q].id);
			printf("...processed\t\tHave a great time!\n\n");
			
			for(ndx1=0; ndx1 < flag; ndx1++) {
				free (seq[ndx1]);
			}
			fclose(fqal);
			free(seq);
			free(length);
			free(seqmap); 
		}
		if(flag >= 1)				
			free(Hitlist);
	}
	#ifdef _WIN32
		sprintf(cmd, "del %s", qseq_fn);
		p_exec(cmd);
		sprintf(cmd, "del %s", qseq_blast_fn);
		p_exec(cmd);
		sprintf(cmd, "del %s %s", fn_p3in, fn_p3out);
		p_exec(cmd);                    /*delete primer3 tmp files*/
	#else
		sprintf(cmd, "rm %s", qseq_fn);
		p_exec(cmd);
		sprintf(cmd, "rm %s", qseq_blast_fn);
		p_exec(cmd);
		sprintf(cmd, "rm %s %s", fn_p3in, fn_p3out);
		p_exec(cmd);                    /*delete primer3 tmp files*/
	#endif
	
	fclose(fnohit);
	fclose(funiseg);
	if(PSPD_RUN_SIM_FLAG==0)
		fclose(fsim);
	return PASS;
}
