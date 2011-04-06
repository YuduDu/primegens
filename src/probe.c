
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
#include "input.h"
#include "probe.h"
#include "uhit.h"

/* Find probes for the given queries */	
int probe_loop (void) { 
	
	/*Local Variables for the function */					
	int q=0, b=0, h=0, n=0;
	int res=0, bloop=0, cloop=0, flag=0, flag1=0;
	int frag_start, frag_end;
	float gc=0.0;
	char qseq_fn[FILENAME], qseq_blast_fn[FILENAME];
	
	char cmd[COMMAND];
	char Expect[CONFIGFORMAT], line[BLINE];
	FILE *fp, *fblast, *fprobe;
			
	res = read_query(fn_query);			/*saved in *query*/
	if(res==FAIL)
		p_exit("problem inside read_query function\n");

	sprintf(cmd, "mkdir %s", tmp_dir);		/*working directory*/
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
	
	fprobe = p_fopen_error_exit(fn_probe, "w");
   
   /* for each query */
	for(q=0; q < q_count; q++) {
		printf("\n\n");
		if(dlevel >=1) {
			printf("%s =====================================",T_S());
			printf("==========================================\n");
		}
		printf("%s %05d) processing %.60s\n",T_S(),q+1,query[q].id);
		if(dlevel >=1) {
			printf("%s =====================================",T_S());
			printf("==========================================\n");
		}
		if((query[q].len < 10)||(strlen(query[q].nseq) < 10)) {
			printf("%s query has no DNA sequence\tskipped.\n", T_S());
			continue;
		}
        frag_start = 0;
        /* if the query is shorter than the probe length */
        if(PROBE_LENGTH >= query[q].len)
            frag_end = query[q].len -1;
        else
            frag_end = PROBE_LENGTH -1;
            
        /* while we haven't hit the end of the query sequence */
        while(frag_end < query[q].len) {
        
        
         /* count the GCs and skip if there are too many or too few */
         gc = 0.00;
			for(n=frag_start; n <= frag_end; n++)
             	if((toupper(query[q].nseq[n])=='G')
						||(toupper(query[q].nseq[n])=='C'))
					gc = gc + 1.0;
			gc = (gc*100)/(float)(frag_end-frag_start+1);
			if((gc < PROBE_MIN_GC)||(gc > PROBE_MAX_GC)) {
				printf("%s start: %d\tend: %d",T_S(),frag_start,frag_end);
				printf("  GC-content %f...skipped!\n", gc);
				frag_start += PROBE_PERIOD;
            	frag_end += PROBE_PERIOD;
				continue;
			}
			
			/* copy the fragments into a temp file for BLASTing */
			fp = p_fopen_error_exit(qseq_fn, "w");
			fprintf(fp, "%s\n", query[q].id);
            for(n=frag_start; n <= frag_end; n++)
                fprintf(fp, "%c", query[q].nseq[n]);
			fprintf(fp, "\n");
			fclose(fp);
			if(DATABASE_TYPE==1) {
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
				
				/* Run BLAST, store the output in qseq_blast_fn */
				p_exec(cmd);
				fblast = p_fopen_error_exit(qseq_blast_fn, "r");
				flag = -1;
				
				/* Read BLAST output and count the hits */
				while(fgets(line, BLINE, fblast)) {
					if((flag > 0) && (line[0]=='\n')) 
						break; 
					else if(strncmp(line,BHITTAG,BHITCOMP)==0)
						flag = 0;
					else if((flag >=0) && (line[0] !='\n')) {
						flag++;
					}
				}
            	fclose(fblast);
			} else if(DATABASE_TYPE==2) {
				p_exit("database type not supported with this algorithm!\n");	
			} else if(DATABASE_TYPE==3) {
				flag=0;
				
				/* loop for each chromosome */
				for(bloop=0, cloop=0; bloop < BLASTLOOP; bloop++) {
					#ifdef _WIN32
						sprintf(cmd, "%s -p blastn -i %s -d %s%d -F F -G %d -E %d \
	   							-e %s -o %s", blastall, qseq_fn, fn_vdb, bloop, 
	      						(BLAST_OGP + BLAST_OIP), 
	      						BLAST_OIP, Expect, qseq_blast_fn);
					#else 
						sprintf(cmd, "%s -p blastn -i %s -d %s%d -F F -G %d -E %d \
	   							-e %s -o %s", blastall, qseq_fn, fn_vdb, bloop,
	      						(BLAST_OGP + BLAST_OIP), 
	      						BLAST_OIP, Expect, qseq_blast_fn);
					#endif 
					/* Run BLAST, store the output in qseq_blast_fn */
					p_exec(cmd);
					fblast = p_fopen_error_exit(qseq_blast_fn, "r");
					flag1 = -1;
					
					/* Read BLAST output and count the hits */
					while(fgets(line, BLINE, fblast)) {
						if((flag1 > 0) && (line[0]=='\n')) 
							break; 
						else if(strncmp(line,BHITTAG,BHITCOMP)==0)
							flag1 = 0;
						else if((flag1 >=0) && (line[0] !='\n')) {
							flag1++;
						}
					}
            	
            	fclose(fblast);
					cloop = cloop + CHRDBINLOOP;
					
					/*update the total hits for all chromosomes */	
					flag += flag1;
					if(dlevel >2) {
						printf("%s start: %d\tend: %d",T_S(),frag_start,frag_end);
						printf("\tbloop: %d cloop: %d\t", bloop, cloop);
						if(flag1 < 0)
							printf("GC-content %f\tno hit found!\n",gc);
						else		
							printf("GC-content %f\tblast hit=%d\n",gc,flag1);
					}
				}
			} else
				p_exit("wrong database format provided!\n");
			
			if(dlevel >=1) {
				printf("%s start: %d\tend: %d",T_S(),frag_start,frag_end);
				printf("\tGC-content %f", gc);
				if(flag < 1)
					flag = 1;		/*some cases blast output is not formed*/
				printf("\tblast hit = %d", flag);
			}
			
			/* output the probes with the few enough hits */
			if(flag <= (PROBE_BLAST_HIT+1)) {
				printf("\tselected!\n");
				fprintf(fprobe,"%s%05d.%04d%04d#", PROBEPREFIX,
													q+1, frag_start, frag_end);
				fprintf(fprobe, ">%s\n", query[q].id);
            	for(n=frag_start; n <= frag_end; n++)
                	fprintf(fprobe, "%c", query[q].nseq[n]);
				fprintf(fprobe, "\n");
				if(PROBE_UNIQUE ==1) {
					break;
				} else {
					frag_start += PROBE_REGION;
					frag_end += PROBE_REGION;
					continue;
				}
			} else
				printf("\n");
			frag_start += PROBE_PERIOD;
            frag_end += PROBE_PERIOD;
        }
	}
	fclose(fprobe);

	/* delete the temporary files */
    #ifdef _WIN32
        sprintf(cmd, "del %s", qseq_fn);
        p_exec(cmd);
        sprintf(cmd, "del %s", qseq_blast_fn);
        p_exec(cmd);
    #else
        sprintf(cmd, "rm %s", qseq_fn);
        p_exec(cmd);
        sprintf(cmd, "rm %s", qseq_blast_fn);
        p_exec(cmd);
    #endif

	return PASS;
}

