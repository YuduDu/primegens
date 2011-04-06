
/*fragment-specific primer design*/

/* header files*/

#include "lib.h"
#include "defs.h"
#include "util.h"
#include "fspd.h"
#include "hbdzn.h"
#include "output.h"
#include "pd.h"
#include "pphit.h"
#include "uhit.h"


short fspd_loop(void) {
	int q=0, b=0, db=0, uo=0, pp=0;
	int res=0, bloop=0, cloop=0, c=0;
	size_t t_size;
	char cmd[COMMAND];
	FILE *fp, *fout4;

	res = read_fspd_sb(fn_query);			/*saved in *query*/
	if(res==FAIL)
		p_exit("problem inside read_fspd_sb function\n");

	sprintf(cmd, "mkdir %s", tmp_dir);
	p_exec(cmd);

    if(DATABASE_TYPE==1) {  /*single database type*/
        res = format_single_database(); 
        if(res==FAIL)
            p_exit("problem while formatdb all databases\n");

    } else if(DATABASE_TYPE==2) { /*multiple databases*/
        res = format_multiple_database();
        if(res==FAIL)
            p_exit("problem while formatdb all databases\n");

    } else if(DATABASE_TYPE==3) { /*chromosome sequence from genome*/
        for(bloop=0, cloop=0; bloop < BLASTLOOP; bloop++) {
            res = format_chromosome_database(bloop, cloop);
            if(res==FAIL)
                p_exit("problem while formatdb all databases\n");
            cloop = cloop + CHRDBINLOOP;
        }
    }


	fout4 = p_fopen_error_exit(fn_out4, "w");	/*record failed queries*/
	
	for(q=0; q < q_count; q++) {
		printf("\n\n");
		printf("%s===========================================",T_S());
		printf("===========================================\n");
		printf("%s %05d) primerdesign for %s\n",T_S(),q+1,query[q].id);
		printf("%s===========================================",T_S());
		printf("===========================================\n");
        /*skip if no sequence available*/
		if((query[q].len < 10)||(strlen(query[q].nseq) < 10)) {
			printf("%s query has no DNA sequence\tskipped.\n", T_S());  
			continue;
		}
		
	/*STEP 1) primer design using primer3 ]**************/
		printf("%s\t\t_____________________________________",T_S());
		printf("_______________________________________________\n");
		printf("%s\t\t STEP 1: primer design using primer3",T_S());
		printf(" for query # %04d [%s]...\n", q+1, query[q].id);
		printf("%s\t\t_____________________________________",T_S());
		printf("_______________________________________________\n");
		run_primer3(query[q]);
		res = save_primer3_output(); 		/*returns pass-fail*/
		if((res==FAIL)||(ppair_count<=0)) { 
			query[q].res = FAIL;
			if((QUERY_FORMAT==3)||(QUERY_FORMAT==4)) {
            	fprintf(fout4, "%s:\t", qinfo[q].chr);
            	fprintf(fout4, "%d-\t", qinfo[q].start);
				fprintf(fout4, "%d\t", qinfo[q].end);
				if(QUERY_FORMAT==4)
                	fprintf(fout4, "%s\n", qinfo[q].desc);
				else 
                	fprintf(fout4, "\n");
			} else if((QUERY_FORMAT==1)||(QUERY_FORMAT==2)) {
            	fprintf(fout4, "%s\t", qinfo[q].id);
				if(QUERY_FORMAT==2)
                	fprintf(fout4, "%s\n", qinfo[q].desc);
				else 
                	fprintf(fout4, "\n");
			}
			if(QUERY_SEQUENCE==1){
				for(b=0; b < query[q].len; b++) {
					fprintf(fout4, "%c", query[q].nseq[b]);
					if((b + 1) % FASTALINE == 0)
						fprintf(fout4, "\n");
				}
				if(query[q].len % FASTALINE != 0)
					fprintf(fout4, "\n");
			}
			continue;						/*skip this query*/
		} else 
			query[q].res = PASS;		/*primer design successful*/
		sprintf(cmd, "rm %s %s", fn_p3in, fn_p3out);
		p_exec(cmd);					/*delete primer3 tmp files*/
	

	/*STEP 2) find unique oligos and save blast hits ]********/
		printf("%s\t\t_____________________________________",T_S());
		printf("_______________________________________________\n");
		printf("%s\t\tSTEP 2: find uoilgos blast records",T_S());
		printf(" for query # %04d [%s]...\n", q+1, query[q].id);
		printf("%s\t\t_____________________________________",T_S());
		printf("_______________________________________________\n");
		
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

        } else if(DATABASE_TYPE==2) {
            for(db=0; db , db_count; db++) { 
                run_uo_blast_multiple(db);
            }
        
        } else if(DATABASE_TYPE==3) {
            printf("%s malloc for uoligo BRecords...\n",T_S());
            for(uo=0; uo < uoligo_count; uo++) { /*brecords for uoligo*/ 
                t_size = BLASTRECORDS*sizeof(BRecord);
                uoligo[uo].brecord=(BRecord *)p_malloc(t_size);
                for(b=0; b < BLASTRECORDS; b++)
                    strcpy(uoligo[uo].brecord[b].chr, fn_db[b].fn);
            }
            printf("%s malloc for uoligo BRecords...done!\n",T_S());
            for(bloop=0, cloop=0; bloop < BLASTLOOP; bloop++) {
                printf("%s blast loop [%02d] with start ",T_S(),bloop);
                printf("chromosome index [%02d] execution...\n",cloop);
                run_uo_blast_chromosome(bloop, cloop); 
                save_blast_output_chromosome(bloop, cloop); 
                sprintf(cmd, "rm %s", fn_bout);
                p_exec(cmd);
                printf("%s blast loop [%02d] with start ",T_S(),bloop);
                printf("chromosome index [%02d] execution...done!\n",cloop);
                cloop = cloop + CHRDBINLOOP;    
            }
        }
		
		sprintf(cmd, "rm %s", fn_bin);
		p_exec(cmd);		
		if(res==FAIL) 
			p_exit("problem in uoligo_brecords_test function.\n");	


	/*STEP 3) build primer pairs chit records ]**************/
		printf("%s\t\t_____________________________________",T_S());
		printf("_______________________________________________\n");
		printf("%s\t\t STEP 3: find primer pair blast records",T_S());
		printf(" for query # %04d [%s]...\n", q+1, query[q].id);
		printf("%s\t\t_____________________________________",T_S());
		printf("_______________________________________________\n");

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
        } else if(DATABASE_TYPE==2) {
            res = ppair_chit_record_multiple(); 
            if(res==FAIL)
                p_exit("problem in saving blast records for pairs\n");
        } else if(DATABASE_TYPE==3) {
            res = ppair_chit_record_chromosome();   
            if(res==FAIL)
                p_exit("problem in saving blast records for pairs\n");
            printf("%s free uoligos blast records...\n", T_S()); 
            for(uo=0; uo < uoligo_count; uo++) { 
                for(cloop=0; cloop < BLASTRECORDS; cloop++) {
                    if(uoligo[uo].brecord[cloop].bh > 0) {
                        free(uoligo[uo].brecord[cloop].bht);
                        free(uoligo[uo].brecord[cloop].bhs);
                        free(uoligo[uo].brecord[cloop].bhe);
                    }
                }
                free(uoligo[uo].brecord);
            }
            free(uoligo);
        }
        
        printf("%s free uoligos blast records...done!\n", T_S()); 
    

	/*STEP 4) build potential hybridization records ]********/
		printf("%s\t\t_____________________________________",T_S());
		printf("_______________________________________________\n");
		printf("%s\t\t STEP 4: find potential hybridization",T_S());
		printf(" for query # %04d [%s]...\n", q+1, query[q].id);
		printf("%s\t\t_____________________________________",T_S());
		printf("_______________________________________________\n");

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
        } else if(DATABASE_TYPE==2) {
        		/* REMOVE? */
            //res = chk_hybridization_multiple(q);      /*qinfo is required*/
            if(res==FAIL)
                p_exit("problem inside chk_hybirization function\n");
        } else if(DATABASE_TYPE==3) {
            res = chk_hybridization_chromosome(q);  /*qinfo is required*/
            if(res==FAIL)
                p_exit("problem inside chk_hybirization function\n");
            printf("%s free ppair chit records...\n", T_S()); 
            for(pp=0; pp < ppair_count; pp++) { 
                for(cloop=0; cloop < BLASTRECORDS; cloop++) {
                    if(ppair[pp].chit[cloop].ppc > 0)
                        free(ppair[pp].chit[cloop].pps);
                    if(ppair[pp].chit[cloop].pmc > 0)
                        free(ppair[pp].chit[cloop].pme);
                }
                free(ppair[pp].chit);
            }
        }
        printf("%s free ppair chit records...done!\n", T_S()); 

	
	/*STEP 5) rank all primer pairs and save into file ]*****/
		printf("%s\t\t_____________________________________",T_S());
		printf("_______________________________________________\n");
		printf("%s\t\t STEP 5: primer ranking and save",T_S());
		printf(" for query # %04d [%s]...\n", q+1, query[q].id);
		printf("%s\t\t_____________________________________",T_S());
		printf("_______________________________________________\n");
		
		res = rank_ppair();				/*rank primer pairs*/ 
		if(res==FAIL)
			p_exit("problem inside primer pair ranking\n");
		
		res = print_output(q);			/*save output*/
		if(res==FAIL)
			p_exit("problem inside print output\n");

		for(pp=0; pp < ppair_count; pp++)
			if(ppair[pp].hbrdn_count > 0)
				free(ppair[pp].hbrdn);
		free(ppair);
		printf("\n\n");
		printf("%s %05d) query # %2d [%.30s]",T_S(),q+1,q+1,query[q].id);
		printf("...processed\t\tHave a great time!\n\n");
	
	}		/*primer design for query(q) done*/
	fclose(fout4);
	
	return PASS;
}

