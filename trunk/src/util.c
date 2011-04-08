
/*
 * this module contains all the basic utilities, which are
 * frequently used in the software for various purposes.
 * More info can be found in the .h file
 */

#include "lib.h"
#include "defs.h"
#include "oligotm.h"
#include "util.h"



/*local utilies*/

/* Get the min of two ints */
int p_min (int x, int y) {
	return ((x < y)?x:y);
}

/* Get the max of two ints */
int p_max (int x, int y) {
	return ((x > y)?x:y);
}

/* Allocate x bytes of memory */
void *p_malloc(size_t x) {
	void *r = malloc(x);
	
	if (r==NULL) {
		printf("\n\t\tmemory allocation failed!\n");
		exit(1);
	}
	return r;
}

/* Get the time as a string */
char *T_S (void) {
	time_t lt;
	struct tm *t;

	lt = time(NULL);
	t = localtime (&lt);
	sprintf(t_Stamp, " %02d:%02d:%02d:", t->tm_hour, 
			t->tm_min, t->tm_sec);
	return t_Stamp;
}

/* Exit the program, logging and displaying msg */
void p_exit(char *msg) {
	char cmd[COMMAND];
	FILE *fp;
	
	/* Open the log file and log msg */ 
	fp = fopen(fn_log, "a");
	if(!fp)
		printf("\n\tfailed to log into log file (%s)!\n", fn_log);	
	fprintf(fp,"\n%s %s\n",T_S(),msg);
	fprintf(fp,"\n\tPRIMEGENS\tFATAL ERROR:\tAbort!\n\n\n");
	fclose(fp);
	
	/* print the message to the user */
	printf("\n%s %s\n",T_S(),msg);
	printf("\n\tPRIMEGENS\tFATAL ERROR:\tAbort!\n\n\n");
	
	/* Delete the temporary directories and files */
	#ifdef _WIN32
		sprintf("del %s\\*", tmp_dir);
	#else
		sprintf("rm %s/*", tmp_dir);
	#endif
		p_exec(cmd); 
		sprintf(cmd, "rmdir %s", tmp_dir);
}

/* Log and display msg */
void p_continue(char *msg) {
	FILE *fp;
	
	/* Open the log file and log msg */
	fp = p_fopen_error_exit(fn_log, "a");
	if(!fp)
		printf("\n\tfailed to log int log file!\n");	
	fprintf(fp,"%s %s\n",T_S(),msg);
	fclose(fp);
	/* print msg */
	printf("%s %s\n",T_S(),msg);
}

/* Open the file, exiting the program if opening fails */
FILE *p_fopen_error_exit(char *fn, char *m) {
	FILE *fp;
	/* Open the file in mode m */
	fp=fopen(fn, m);
	/* ABORT if file opening failed */
	if(!fp) {
		printf("%s file %s\t open error.\n",T_S(),fn);
		p_exit("\n\t\tprimegens abort!\n\n");
		return NULL;
	}
	return fp;
}

/* Execute the give system cmd */
void p_exec(char *cmd) {
	int state;

	if(dlevel > 3)
		printf("%s executing %s...", T_S(), cmd);
	state = system(cmd); 
	/* If command fails, exit program */
	if(state < 0) {
		printf("%s command output = %d\tfailed!\n",T_S(),state);
		exit(10);
	}	
	if(dlevel > 3)
		printf("successful!\n");
	return;
}

/* Compare two oligonucleotide strings */
int p_oligo_strcmp(char *str1, char *str2) {
	int i=0, limit=0;
	int len1=0, len2=0;
	int ndx1=0, ndx2=0;
	
	/* Set the limit to min of len1, len2, MINBIND */
	len1 = (int)strlen(str1);
	len2 = (int)strlen(str2);
	if(p_min(len1, len2) > MINBIND) 
		limit = MINBIND;
	else 
		limit = p_min(len1, len2);
		
	/* Compare each nucleotide, return nonzero if different */	
	for(i=1; i <= limit; i++) {   
		ndx1 = len1 - i;	  /* char index start from 0*/
		ndx2 = len2 - i;	  /* char index start from 0*/
		if(str1[ndx1]!=str2[ndx2])
			return (str1[ndx1]-str2[ndx2]);
	}
	/* return zero since strings match */
	return 0;
}

/* Copy str2 into str1 */
void p_oligo_strcpy(char *str1, char *str2) {
	int i=0, ndx=0, len=0;

	len = (int)strlen(str2);
	if(len < MINBIND) {
		strcpy(str1, str2);
		return;
	}
	/* copy only the last MINBIND nucleotides */
	str1[MINBIND] = 0;	/*important*/
	for(i=MINBIND-1; i >= 0; i--) {
		ndx = len - MINBIND +i;
		str1[i] = str2[ndx];
	}

	return;
}


/*final cleaning*/
char p_clean(void) {
	char cmd[COMMAND], t_vdb[FILENAME];	
   	int bloop;
	
	printf("%s final cleaning...\n", T_S());
	/* deallocate memory */
	free(fn_db);
	free(query);
	free(qinfo);
	printf("%s free global variable memory done!\n", T_S());

	/* Delete temporary database files */
	if(DATABASE_TYPE==1) {
		sprintf(t_vdb, "%s", fn_vdb);
		sprintf(cmd,"rm %s.nhr %s.nin %s.nsd",t_vdb,t_vdb,t_vdb);
		p_exec(cmd);
		sprintf(cmd,"rm %s.nsi %s.nsq", t_vdb, t_vdb);
		p_exec(cmd);
	}
	else if(DATABASE_TYPE==3) {
		for(bloop=0; bloop < BLASTLOOP; bloop++) {
			sprintf(t_vdb, "%s%d", fn_vdb, bloop);
			sprintf(cmd,"rm %s.nhr %s.nin %s.nsd",t_vdb,t_vdb,t_vdb);
			p_exec(cmd);
			sprintf(cmd,"rm %s.nsi %s.nsq", t_vdb, t_vdb);
			p_exec(cmd);
		}
	}
	/* remove temporary directories */
	sprintf(cmd, "rmdir %s", tmp_dir);
	p_exec(cmd); 
}

		
/*get sequence for given query from database*/
short get_nseq_from_db(char *qid, char *nseq, int *len) {
	FILE *fp;
	int flag=0;
	char line[ELLINE], tnseq[QNSIZE], tdbname[FILENAME];
	char tid[SEQUENCEID], *tok;

	if(dlevel >=3)
		printf("%s get sequence of %s...\n",T_S(), qid);
	*nseq = '\0';
	if(DATABASE_TYPE==1) {
		#ifdef _WIN32
			sprintf(tdbname, "%s\\%s", dbpath, fn_db0.fn);
		#else 
			sprintf(tdbname, "%s/%s", dbpath, fn_db0.fn);
		#endif
		/* Open the database */
		fp = p_fopen_error_exit(tdbname, "r");
		/* Read lines until we break or hit EOF */
		while(fgets(line, ELLINE, fp)) {
			/* if we hit a new sequence... */
			if(strncmp(line, ">", 1)==0) {
				/* If we already found our sequence, then it just ended, so break */
				if(flag==1) {
					strcat(nseq, "\0");
					break;
				}
				/* Otherwise, check if this is our query sequence */
				if(strncmp(qid, ">", 1)==0) {
					if(strncmp(line, qid, strlen(qid))==0) {
						flag = 1;
					}
				} else {
				/* Query might not start with '>' */
					sscanf(line, ">%s", tid);
					if(strncmp(tid, qid, strlen(tid))==0) {
						flag = 1;
					}
				}
				/* We're in our query sequence, so append this line */
			} else if(flag==1) {
				sscanf(line, "%s", tnseq);
				strcat(nseq, tnseq);
			}
		}
		fclose(fp);
		*len = (int)strlen(nseq); 	
	} else if(DATABASE_TYPE==2) {
		p_exit("feature not developed yet!\n");
	}
	if(dlevel >=3) {
		printf("%s qid=[%s] length=%d\n",T_S(), qid, *len);
		if(dlevel==4)
			printf("%s\n", nseq);
		printf("%s get sequence of %s...done!\n",T_S(), qid);
	}
	return PASS;
}

/* Get hybrid sequences from the database */
short get_hybd_nseq_from_db(char *qid, int ls, int le, int rs, int re, char *hlseq, char *hrseq, int *hllen, int *hrlen) {
        FILE *fp;
        int flag=0,b=0,count=0;
        char line[MLINE], tnseq[QNSIZE], tdbname[FILENAME];
        char tid[SEQUENCEID], *tok,base;
        char hl[30],hr[30];
        if(dlevel >=3)
           printf("%s get hit  sequence for primer pair of %s...\n",T_S(), qid);
        *hlseq = '\0';
        *hrseq = '\0';
        if(DATABASE_TYPE==1) {
                #ifdef _WIN32
                        sprintf(tdbname, "%s\\%s", dbpath, fn_db0.fn);
                #else
                        sprintf(tdbname, "%s/%s", dbpath, fn_db0.fn);
                #endif
                /* Open the file */
                fp = p_fopen_error_exit(tdbname, "r");
                /*Read lines until we break or hit the EOF */
                while(fgets(line, MLINE, fp)) {
                			/* If we hit a new sequence... */
                        if(strncmp(line, ">", 1)==0) {
                        			/* If we already found our sequence, then break */
                                if(flag==1) {
                                        strcat(hlseq, "\0");
                                        break;
                                }
                                /* Otherwise, check if this is our sequence */
                               if(strncmp(qid, ">", 1)==0) {
                                        if(strncmp(line, qid, strlen(qid))==0) {
                                                flag = 1;
                                        }
                               /* the query might not start with '>' */
                               } else {
                                        sscanf(line, ">%s", tid);
                                        if(strncmp(tid, qid, strlen(tid))==0) {
                                                flag = 1;
                                        }
                                }
                        /* We're in our query sequence, so append this line */
                        } else if(flag==1) {
                                count=0;
                                /* Read until we get to the starting line */
                                while((count + CHR_FASTA_LINE)< ls){
                                     fgets(line, MLINE, fp);
                                     count += CHR_FASTA_LINE;
                                }
                                /* Read until we get to the starting base */
                                for(b=count+1; b < ls; b++){
                                    fscanf(fp,"%c", &base);
                                }
                                /* Read in the sequence to hl */
                                for(b=ls; b <= le;) {
                                    fscanf(fp,"%c", &base);
                                    if(base != '\n'){
                                       hl[b-ls] = toupper(base);
                                       b++;
                                    }
                                }
                                strcpy(hlseq,hl);
                                strcat(hlseq, "\0");
                                rewind(fp);
                                count=0;
                                /* Read until we get to the starting line */
                                while((count + CHR_FASTA_LINE)< rs){
                                     fgets(line, MLINE, fp);
                                     count += CHR_FASTA_LINE;
                                }
                                /* Read until we get to the starting base */
                                for(b=count+1; b < rs; b++){
                                    fscanf(fp,"%c", &base);
                                }
                                /* Read in the sequence to hr */
                                for(b=rs; b <= re;) {
                                    fscanf(fp,"%c", &base);
                                    if(base != '\n'){
                                       hr[b-rs] = toupper(base);
                                       b++;
                                    }
                                }
                                strcpy(hrseq,hr);
                                strcat(hrseq, "\0");
                            		/* Read both sequences, so break */
                                break;
                        }
                }
                fclose(fp);
                /* update the lengths by reference */
                *hllen = (int)strlen(hlseq);
                *hrlen = (int)strlen(hrseq);
        } 
        return PASS;
}

/* Calculate potential for two oligonucleotides */
short get_oligo_tm_for_hybd(char *hlseq,char *hrseq,char *lseq,char *rseq,int *potl){
        char h_lseq[HLEN],h_rseq[HLEN],p_lseq[OLIGO],p_rseq[OLIGO];
        char h_rc_lseq[HLEN],h_rc_rseq[HLEN],p_rc_lseq[OLIGO],p_rc_rseq[OLIGO];
        char h_check_lseq[6],h_check_rseq[6];
        int h_lseq_len=0, h_rseq_len=0,p_lseq_len=0,p_rseq_len=0;
        int i, j,k;    
        
        h_lseq_len=strlen(hlseq); 
        h_rseq_len=strlen(hrseq);
        p_lseq_len=strlen(lseq);
        p_rseq_len=strlen(rseq);

        strcpy(h_lseq,hlseq);
        strcpy(h_rseq,hrseq);
        strcpy(p_lseq,lseq);
        strcpy(p_rseq,rseq);

        strcat(h_lseq,"\0");
        strcat(h_rseq,"\0");
        strcat(p_lseq,"\0");
        strcat(p_rseq,"\0");

        j = 0;
        k = h_lseq_len-6;
        /* Check the matching of the left primer and left hybrid sequences, in FASTA format */
        for(i=(p_lseq_len-6);i<p_lseq_len;i++){
         if(toupper(p_lseq[i])=='A'){
         if(toupper(h_lseq[k])=='A'){
            h_check_lseq[j]='A';
         }else if (toupper(h_lseq[k])=='T'){
            h_check_lseq[j]='H';
         }else if (toupper(h_lseq[k])=='G'){
            h_check_lseq[j]='I';
         }else if (toupper(h_lseq[k])=='C'){
            h_check_lseq[j]='J';
         }
    }
    if(toupper(p_lseq[i])=='C'){
       if (toupper(h_lseq[k])=='C'){
          h_check_lseq[j]='C';
       }else if (toupper(h_lseq[k])=='G'){
          h_check_lseq[j]='K';
       }else if (toupper(h_lseq[k])=='T'){
          h_check_lseq[j]='L';
       }else if (toupper(h_lseq[k])=='A'){
          h_check_lseq[j]='M';
       }
    }
    if(toupper(p_lseq[i])=='G'){
       if (toupper(h_lseq[k])=='G'){
          h_check_lseq[j]='G';
       }else if (toupper(h_lseq[k])=='C'){
          h_check_lseq[j]='V';
       }else if (toupper(h_lseq[k])=='T'){
          h_check_lseq[j]='O';
       }else if (toupper(h_lseq[k])=='A'){
          h_check_lseq[j]='P';
       }   
    }
    if(toupper(p_lseq[i])=='T'){
       if (toupper(h_lseq[k])=='T'){
          h_check_lseq[j]='T';
       }else if (toupper(h_lseq[k])=='A'){
          h_check_lseq[j]='Q';
       }else if (toupper(h_lseq[k])=='G'){
          h_check_lseq[j]='R';
       }else if (toupper(h_lseq[k])=='C'){
          h_check_lseq[j]='S';
       }
   }
    j++;
k++;
}
 h_check_lseq[j]= '\0';

j=0;
k=0;

/* Check the matching of the right primer and right hybrid sequences, in FASTA format */
for(i=(p_rseq_len-6);i<p_rseq_len;i++){
    if(toupper(p_rseq[i])=='A'){
       if(toupper(h_rseq[k])=='T'){
          h_check_rseq[j]='A';
       }else if (toupper(h_rseq[k])=='A'){
          h_check_rseq[j]='H';
       }else if (toupper(h_rseq[k])=='C'){
          h_check_rseq[j]='I';
       }else if (toupper(h_rseq[k])=='G'){
          h_check_rseq[j]='J';
       }
    }
    if(toupper(p_rseq[i])=='C'){
       if (toupper(h_rseq[k])=='G'){
          h_check_rseq[j]='C';
       }else if (toupper(h_rseq[k])=='C'){
          h_check_rseq[j]='K';
       }else if (toupper(h_rseq[k])=='A'){
          h_check_rseq[j]='L';
       }else if (toupper(h_rseq[k])=='T'){
          h_check_rseq[j]='M';
       }
    }
    if(toupper(p_rseq[i])=='G'){
       if (toupper(h_rseq[k])=='C'){
          h_check_rseq[j]='G';
       }else if (toupper(h_rseq[k])=='G'){
          h_check_rseq[j]='N';
       }else if (toupper(h_rseq[k])=='A'){
          h_check_rseq[j]='O';
       }else if (toupper(h_rseq[k])=='T'){
          h_check_rseq[j]='P';
       }
    }
    if(toupper(p_rseq[i])=='T'){
       if (toupper(h_rseq[k])=='A'){
          h_check_rseq[j]='T';
       }else if (toupper(h_rseq[k])=='T'){
          h_check_rseq[j]='Q';
       }else if (toupper(h_rseq[k])=='C'){
          h_check_rseq[j]='R';
       }else if (toupper(h_rseq[k])=='G'){
          h_check_rseq[j]='S';
       }
    }
    j++;
    k++;
}
 h_check_rseq[j]= '\0';

	
	/* If there is a mismatch, set potl to 0 and return */
    if(((h_check_lseq[5] == 'H') || (h_check_lseq[5] == 'I') || (h_check_lseq[5] == 'J') || \
        (h_check_lseq[5] == 'K') || (h_check_lseq[5] == 'L') || (h_check_lseq[5] == 'M') || \
        (h_check_lseq[5] == 'V') || (h_check_lseq[5] == 'O') || (h_check_lseq[5] == 'P') || \
        (h_check_lseq[5] == 'Q') || (h_check_lseq[5] == 'R') || (h_check_lseq[5] == 'S'))|| \
       ((h_check_lseq[4] == 'H') || (h_check_lseq[4] == 'I') || (h_check_lseq[4] == 'J') || \
        (h_check_lseq[4] == 'K') || (h_check_lseq[4] == 'L') || (h_check_lseq[4] == 'M') || \
        (h_check_lseq[4] == 'V') || (h_check_lseq[4] == 'O') || (h_check_lseq[4] == 'P') || \
        (h_check_lseq[4] == 'Q') || (h_check_lseq[4] == 'R') || (h_check_lseq[4] == 'S'))){
       *potl=0;
       return PASS;
    }
    if(((h_check_rseq[0] == 'H') || (h_check_rseq[0] == 'I') || (h_check_rseq[0] == 'J') || \
        (h_check_rseq[0] == 'K') || (h_check_rseq[0] == 'L') || (h_check_rseq[0] == 'M') || \
        (h_check_rseq[0] == 'V') || (h_check_rseq[0] == 'O') || (h_check_rseq[0] == 'P') || \
        (h_check_rseq[0] == 'Q') || (h_check_rseq[0] == 'R') || (h_check_rseq[0] == 'S')) || \
       ((h_check_rseq[1] == 'H') || (h_check_rseq[1] == 'I') || (h_check_rseq[1] == 'J') || \
        (h_check_rseq[1] == 'K') || (h_check_rseq[1] == 'L') || (h_check_rseq[1] == 'M') || \
        (h_check_rseq[1] == 'V') || (h_check_rseq[1] == 'O') || (h_check_rseq[1] == 'P') || \
        (h_check_rseq[1] == 'Q') || (h_check_rseq[1] == 'R') || (h_check_rseq[1] == 'S'))){
       *potl=0;
       return PASS;
    }

    if(((h_check_lseq[1] == 'H') || (h_check_lseq[1] == 'I') || (h_check_lseq[1] == 'J') || \
        (h_check_lseq[1] == 'K') || (h_check_lseq[1] == 'L') || (h_check_lseq[1] == 'M') || \
        (h_check_lseq[1] == 'V') || (h_check_lseq[1] == 'O') || (h_check_lseq[1] == 'P') || \
        (h_check_lseq[1] == 'Q') || (h_check_lseq[1] == 'R') || (h_check_lseq[1] == 'S')) && \
       ((h_check_lseq[2] == 'H') || (h_check_lseq[2] == 'I') || (h_check_lseq[2] == 'J') || \
        (h_check_lseq[2] == 'K') || (h_check_lseq[2] == 'L') || (h_check_lseq[2] == 'M') || \
        (h_check_lseq[2] == 'V') || (h_check_lseq[2] == 'O') || (h_check_lseq[2] == 'P') || \
        (h_check_lseq[2] == 'Q') || (h_check_lseq[2] == 'R') || (h_check_lseq[2] == 'S'))&& \
       ((h_check_lseq[3] == 'H') || (h_check_lseq[3] == 'I') || (h_check_lseq[3] == 'J') || \
        (h_check_lseq[3] == 'K') || (h_check_lseq[3] == 'L') || (h_check_lseq[3] == 'M') || \
        (h_check_lseq[3] == 'V') || (h_check_lseq[3] == 'O') || (h_check_lseq[3] == 'P') || \
        (h_check_lseq[3] == 'Q') || (h_check_lseq[3] == 'R') || (h_check_lseq[3] == 'S'))){
       *potl=0;
       return PASS;
    }
    if(((h_check_rseq[2] == 'H') || (h_check_rseq[2] == 'I') || (h_check_rseq[2] == 'J') || \
        (h_check_rseq[2] == 'K') || (h_check_rseq[2] == 'L') || (h_check_rseq[2] == 'M') || \
        (h_check_rseq[2] == 'V') || (h_check_rseq[2] == 'O') || (h_check_rseq[2] == 'P') || \
        (h_check_rseq[2] == 'Q') || (h_check_rseq[2] == 'R') || (h_check_rseq[2] == 'S')) && \
       ((h_check_rseq[3] == 'H') || (h_check_rseq[3] == 'I') || (h_check_rseq[3] == 'J') || \
        (h_check_rseq[3] == 'K') || (h_check_rseq[3] == 'L') || (h_check_rseq[3] == 'M') || \
        (h_check_rseq[3] == 'V') || (h_check_rseq[3] == 'O') || (h_check_rseq[3] == 'P') || \
        (h_check_rseq[3] == 'Q') || (h_check_rseq[3] == 'R') || (h_check_rseq[3] == 'S'))&& \
       ((h_check_rseq[4] == 'H') || (h_check_rseq[4] == 'I') || (h_check_rseq[4] == 'J') || \
        (h_check_rseq[4] == 'K') || (h_check_rseq[4] == 'L') || (h_check_rseq[4] == 'M') || \
        (h_check_rseq[4] == 'V') || (h_check_rseq[4] == 'O') || (h_check_rseq[4] == 'P') || \
        (h_check_rseq[4] == 'Q') || (h_check_rseq[4] == 'R') || (h_check_rseq[4] == 'S'))){
       *potl=0;
       return PASS;
    }

    if((((h_check_lseq[2] == 'H') || (h_check_lseq[2] == 'I') || (h_check_lseq[2] == 'J') || \
         (h_check_lseq[2] == 'K') || (h_check_lseq[2] == 'L') || (h_check_lseq[2] == 'M') || \
         (h_check_lseq[2] == 'N') || (h_check_lseq[2] == 'O') || (h_check_lseq[2] == 'P') || \
         (h_check_lseq[2] == 'Q') || (h_check_lseq[2] == 'R') || (h_check_lseq[2] == 'S')) && \
        ((h_check_lseq[3] == 'H') || (h_check_lseq[3] == 'I') || (h_check_lseq[3] == 'J') || \
         (h_check_lseq[3] == 'K') || (h_check_lseq[3] == 'L') || (h_check_lseq[3] == 'M') || \
         (h_check_lseq[3] == 'N') || (h_check_lseq[3] == 'O') || (h_check_lseq[3] == 'P') || \
         (h_check_lseq[3] == 'Q') || (h_check_lseq[3] == 'R') || (h_check_lseq[3] == 'S')))|| \
       (((h_check_lseq[2] == 'H') || (h_check_lseq[2] == 'I') || (h_check_lseq[2] == 'J') || \
         (h_check_lseq[2] == 'K') || (h_check_lseq[2] == 'L') || (h_check_lseq[2] == 'M') || \
         (h_check_lseq[2] == 'N') || (h_check_lseq[2] == 'O') || (h_check_lseq[2] == 'P') || \
         (h_check_lseq[2] == 'Q') || (h_check_lseq[2] == 'R') || (h_check_lseq[2] == 'S'))&& \
        ((h_check_lseq[1] == 'H') || (h_check_lseq[1] == 'I') || (h_check_lseq[1] == 'J') || \
         (h_check_lseq[1] == 'K') || (h_check_lseq[1] == 'L') || (h_check_lseq[1] == 'M') || \
         (h_check_lseq[1] == 'N') || (h_check_lseq[1] == 'O') || (h_check_lseq[1] == 'P') || \
         (h_check_lseq[1] == 'Q') || (h_check_lseq[1] == 'R') || (h_check_lseq[1] == 'S')))||\
       (((h_check_lseq[0] == 'H') || (h_check_lseq[0] == 'I') || (h_check_lseq[0] == 'J') || \
         (h_check_lseq[0] == 'K') || (h_check_lseq[0] == 'L') || (h_check_lseq[0] == 'M') || \
         (h_check_lseq[0] == 'N') || (h_check_lseq[0] == 'O') || (h_check_lseq[0] == 'P') || \
         (h_check_lseq[0] == 'Q') || (h_check_lseq[0] == 'R') || (h_check_lseq[0] == 'S'))&& \
        ((h_check_lseq[1] == 'H') || (h_check_lseq[1] == 'I') || (h_check_lseq[1] == 'J') || \
         (h_check_lseq[1] == 'K') || (h_check_lseq[1] == 'L') || (h_check_lseq[1] == 'M') || \
         (h_check_lseq[1] == 'N') || (h_check_lseq[1] == 'O') || (h_check_lseq[1] == 'P') || \
         (h_check_lseq[1] == 'Q') || (h_check_lseq[1] == 'R') || (h_check_lseq[1] == 'S')))){
       *potl=0;
       return PASS;
    }
    if((((h_check_rseq[2] == 'H') || (h_check_rseq[2] == 'I') || (h_check_rseq[2] == 'J') || \
         (h_check_rseq[2] == 'K') || (h_check_rseq[2] == 'L') || (h_check_rseq[2] == 'M') || \
         (h_check_rseq[2] == 'N') || (h_check_rseq[2] == 'O') || (h_check_rseq[2] == 'P') || \
         (h_check_rseq[2] == 'Q') || (h_check_rseq[2] == 'R') || (h_check_rseq[2] == 'S')) && \
        ((h_check_rseq[3] == 'H') || (h_check_rseq[3] == 'I') || (h_check_rseq[3] == 'J') || \
         (h_check_rseq[3] == 'K') || (h_check_rseq[3] == 'L') || (h_check_rseq[3] == 'M') || \
         (h_check_rseq[3] == 'N') || (h_check_rseq[3] == 'O') || (h_check_rseq[3] == 'P') || \
         (h_check_rseq[3] == 'Q') || (h_check_rseq[3] == 'R') || (h_check_rseq[3] == 'S')))|| \
       (((h_check_rseq[3] == 'H') || (h_check_rseq[3] == 'I') || (h_check_rseq[3] == 'J') || \
         (h_check_rseq[3] == 'K') || (h_check_rseq[3] == 'L') || (h_check_rseq[3] == 'M') || \
         (h_check_rseq[3] == 'N') || (h_check_rseq[3] == 'O') || (h_check_rseq[3] == 'P') || \
         (h_check_rseq[3] == 'Q') || (h_check_rseq[3] == 'R') || (h_check_rseq[3] == 'S'))&& \
        ((h_check_rseq[4] == 'H') || (h_check_rseq[4] == 'I') || (h_check_rseq[4] == 'J') || \
         (h_check_rseq[4] == 'K') || (h_check_rseq[4] == 'L') || (h_check_rseq[4] == 'M') || \
         (h_check_rseq[4] == 'N') || (h_check_rseq[4] == 'O') || (h_check_rseq[4] == 'P') || \
         (h_check_rseq[4] == 'Q') || (h_check_rseq[4] == 'R') || (h_check_rseq[4] == 'S')))|| \
       (((h_check_lseq[5] == 'H') || (h_check_lseq[5] == 'I') || (h_check_lseq[5] == 'J') || \
         (h_check_lseq[5] == 'K') || (h_check_lseq[5] == 'L') || (h_check_lseq[5] == 'M') || \
         (h_check_lseq[5] == 'N') || (h_check_lseq[5] == 'O') || (h_check_lseq[5] == 'P') || \
         (h_check_lseq[5] == 'Q') || (h_check_lseq[5] == 'R') || (h_check_lseq[5] == 'S'))&& \
        ((h_check_lseq[4] == 'H') || (h_check_lseq[4] == 'I') || (h_check_lseq[4] == 'J') || \
         (h_check_lseq[4] == 'K') || (h_check_lseq[4] == 'L') || (h_check_lseq[4] == 'M') || \
         (h_check_lseq[4] == 'N') || (h_check_lseq[4] == 'O') || (h_check_lseq[4] == 'P') || \
         (h_check_lseq[4] == 'Q') || (h_check_lseq[4] == 'R') || (h_check_lseq[4] == 'S')))){
       *potl=0;
       return PASS;
    }
	float l,r;
	int PTS=1;
	/* Calculate the binding energies */
   l = oligodg(h_check_lseq,PTS);
   r = oligodg(h_check_rseq,PTS);
    
   /* If the energies are both less than Primer Max End Stability, then there is potential */
   if((l<PMES) && (r<PMES))
   {
      *potl =1;
   }
	else{
		*potl=0;
	}
	return PASS;
}


/*get sequence for given query from genome*/
short get_nseq_from_chromosome(int q) {
	FILE *fp;
	int flag=0, len=0, b=0, count=0;
	int tstart=0, tend=0; 
	char line[MLINE], base, tdbname[FILENAME];
	char tid[SEQUENCEID], *tok;
	short d=0;

	if(dlevel >=3) {
		printf("%s retrieve sequence [%s:",T_S(), qinfo[q].chr);
		printf("%d-%d]...\n", qinfo[q].start, qinfo[q].end);
	}
	/* If the sequence goes backwards, reverse it */
	if(qinfo[q].end < qinfo[q].start) {
		tstart = qinfo[q].end;
		tend = qinfo[q].start;
	}
	else {  
		tstart = qinfo[q].start;
		tend = qinfo[q].end;
	}
	for(b=0; b < QNSIZE; b++)
		qinfo[q].nseq[b] = '\0'; 
	/* Search the database files for the chromosome */
	for(d=0; d < db_count; d++) {
		strcpy(tdbname, fn_db[d].fn);
      tok = strtok(tdbname,". |\0");   /*find first name of file*/
      /* If we find the chromosome */
		if(strcmp(qinfo[q].chr, tok)==0) {
			#ifdef _WIN32
				sprintf(tdbname, "%s\\%s", dbpath, fn_db[d].fn);
			#else 
				sprintf(tdbname, "%s/%s", dbpath, fn_db[d].fn);
			#endif
			if(dlevel==4)
				printf("%s retrieve nseq from %s...\n", T_S(), tdbname);
			fp = p_fopen_error_exit(tdbname, "r");
			fgets(line, MLINE, fp); /*first line is identifier*/
			count=0;
			/* Read through the lines until we get to the start line */
			while((count + CHR_FASTA_LINE)< tstart){
				fgets(line, MLINE, fp);
				count += CHR_FASTA_LINE;
			}
			/* Read through the bases until we get to the start base */
			for(b=count+1; b < tstart; b++){
				fscanf(fp,"%c", &base);
			}
			/* If forward, read in normal */
			if(qinfo[q].start < qinfo[q].end) {
				for(b=tstart; b <= tend;) {
					fscanf(fp,"%c", &base);
					if(base != '\n'){
						qinfo[q].nseq[b-tstart] = toupper(base);
						b++;
					}
				}
			/* If backward, read in backwards and flip bases */
			} else {
				for(b=tstart; b <= tend;) {
					fscanf(fp,"%c", &base);
					if(base != '\n'){
						if(toupper(base)=='A')
							qinfo[q].nseq[(tend-tstart)-(-b-tstart)] = 'T';
						if(toupper(base)=='T')
							qinfo[q].nseq[(tend-tstart)-(-b-tstart)] = 'A';
						if(toupper(base)=='G')
							qinfo[q].nseq[(tend-tstart)-(-b-tstart)] = 'C';
						if(toupper(base)=='C')
							qinfo[q].nseq[(tend-tstart)-(-b-tstart)] = 'G';
						b++;
					}
				}
			}
			strcat(qinfo[q].nseq, "\0");
			fclose(fp);
			break;
		}
	}
	/* update the query struct */
	strcpy(query[q].nseq, qinfo[q].nseq);
	query[q].len = strlen(query[q].nseq);
	if(dlevel >=3) {
		printf("%s qid=[%.30s] length=%d\n",T_S(), qinfo[q].id, query[q].len);
		if(dlevel==4)
			printf("%s\n", qinfo[q].nseq);
		printf("%s get sequence of %.30s...done!\n",T_S(), qinfo[q].id);
	}

	return PASS;
}

/*get sequence for given query from genome*/
short get_hybd_from_chromosome(char *qid, int ls, int le, int rs, int re, char *chr, char *hlseq, char *hrseq) {
	FILE *fp;
	int flag=0, len=0, b=0, count=0;
	int tstart=0, tend=0;
	char line[MLINE], base, tdbname[FILENAME];
	char tid[SEQUENCEID], *tok,*tok1;
	short d=0;
	char hl[30],hr[30];
	*hlseq = '\0';
	*hrseq = '\0';

	/* Get the start and end of the left sequence, reversing if backwards */
	if(le < ls) {
		tstart = ls-14;
		tend = ls;
	}
	else {
		tstart = ls;
		tend = ls+14;
	}
	/* find the database file we want */
	for(d=0; d < db_count; d++) {
		strcpy(tdbname, fn_db[d].fn);
		tok = strtok(tdbname,". _|\0");  /*find first name of file*/
		tok1 = strtok(chr,". _|\0");  /*find first name of file*/

		if(strcmp(tok1, tok)==0) {
		#ifdef _WIN32
			sprintf(tdbname, "%s\\%s", dbpath, fn_db[d].fn);
		#else
			sprintf(tdbname, "%s/%s", dbpath, fn_db[d].fn);
		#endif
			if(dlevel==4)
				printf("%s retrieve nseq from %s...\n", T_S(), tdbname);
			fp = p_fopen_error_exit(tdbname, "r");
			fgets(line, MLINE, fp); /*first line is identifier*/
			count=0;
			/* Read lines until we hit the starting line */
			while((count + CHR_FASTA_LINE)< tstart){
				fgets(line, MLINE, fp);
				count += CHR_FASTA_LINE;
			}
			/* Read Bases until we hit the starting base */
			for(b=count+1; b < ls; b++){
				fscanf(fp,"%c", &base);
			}
			/* If forward, copy normally */
			if(ls < le) {
				for(b=tstart; b <= tend;) {
					fscanf(fp,"%c", &base);
					if(base != '\n'){
						hl[b-tstart] = toupper(base);
						b++;
					}
				}
			/* Otherwise, read in backwards and swap bases */
			} else {
				for(b=tstart; b <= tend;) {
					fscanf(fp,"%c", &base);
					if(base != '\n'){
						if(toupper(base)=='A')
							hl[(tend-tstart)-(-b-tstart)] = 'T';
						if(toupper(base)=='T')
							hl[(tend-tstart)-(-b-tstart)] = 'A';
						if(toupper(base)=='G')
							hl[(tend-tstart)-(-b-tstart)] = 'C';
						if(toupper(base)=='C')
							hl[(tend-tstart)-(-b-tstart)] = 'G';
						b++;
					}
				}
			}
			strcpy(hlseq,hl);
			strcat(hlseq,"\0");
			rewind(fp);

			/* Get the start and end of the right sequence, reversing if backwards */
			if(re < rs) {
				tstart = rs-14;
				tend = rs;
			}
			else {
				tstart = rs;
				tend = rs+14;
			}

			count=0;
			/* Read lines until we hit the starting line */
			while((count + CHR_FASTA_LINE)< tstart){
				fgets(line, MLINE, fp);
				count += CHR_FASTA_LINE;
			}
			/* Read Bases until we hit the starting base */
			for(b=count+1; b < rs; b++){
				fscanf(fp,"%c", &base);
			}
			/* If forward, copy normally */
			if(ls < le) {
				for(b=tstart; b <= tend;) {
					fscanf(fp,"%c", &base);
					if(base != '\n'){
						hr[b-tstart] = toupper(base);
						b++;
					}
				}
			/* Otherwise, read in backwards and swap bases */
			} else {
				for(b=tstart; b <= tend;) {
					fscanf(fp,"%c", &base);
					if(base != '\n'){
						if(toupper(base)=='A')
							hr[(tend-tstart)-(-b-tstart)] = 'T';
						if(toupper(base)=='T')
							hr[(tend-tstart)-(-b-tstart)] = 'A';
						if(toupper(base)=='G')
							hr[(tend-tstart)-(-b-tstart)] = 'C';
						if(toupper(base)=='C')
							hr[(tend-tstart)-(-b-tstart)] = 'G';
						b++;
					}
				}
			}
			strcpy(hrseq,hr);
			strcat(hrseq,"\0");

			int hllen = (int)strlen(hlseq);
			int hrlen = (int)strlen(hrseq);

			fclose(fp);
			break;
		}
	}
	return PASS;
}



/*pspd function*/
int GetFitnessScore(char *tseq,int tposn,char *qseq, int qposn) { 
	if (tseq[tposn] == qseq[qposn]) 
		return -1; 
	else 
		return 1;
}

/*read query input file to get its format*/
short read_query_toget_format(char *fn) {
	FILE *fp;
	int count=0, qf1=0, qf2=0, qf3=0, qf4=0;
	int start=0, end=0, c=0;
	char line[ELLINE], line1[ELLINE], *tok;
	
	printf("%s read query input %s to get its format...\n", T_S(), fn);
	fp = p_fopen_error_exit(fn, "r");
	
	/* count the lines, ignoring lines with less than 2 characters */
	for(count=0; fgets(line, ELLINE, fp);)
		if(strlen(line)>=2)
			count++;	
	rewind(fp);
	/* read lines until EOF */
	for(qf1=0, qf2=0, qf3=0, qf4=0; fgets(line, ELLINE, fp);) {
		if(strlen(line) < 2)
			continue;			/*ignore this line*/
		if(line[0]=='>') {
			qf1++;
			strcpy(line1, line);
			tok = strtok(line1, " #\t\n");
			/* if there is a ' ','#',tab, or newline before the end of the line */
			if(strlen(line)-strlen(tok) >2)
				qf2++;	
		}
		if(strncmp(line, ">chr", 4)==0) {
			strcpy(line1, line);
			tok = strtok(line1, ": #\t\n");
			tok = strtok(NULL,  "- #\t\n");
			strcat(tok, "\0");
			start = atoi(tok);
			/* ensure that the second token in the line is a valid start */
			if((start < 0)||(start >999999999)) 
				p_exit("query format must be >chrX:start-end!\n");	
			tok = strtok(NULL,  "- #\t\n");
			strcat(tok, "\0");
			end = atoi(tok);
			/* ensure that the second token in the line is a valid end */
			if((end < 0)||(end >999999999)) 
				p_exit("query format must be >chrX:start-end!\n");	
			qf3++;
			strcpy(line1, line);
			/* if there is a ' ','#',tab, or newline before the end of the line */
			tok = strtok(line1, " #\t\n");
			if(strlen(line)-strlen(tok) >2)
				qf4++;
		}
	}
	/* if all the lines are comments */
	if((count==qf1)||(count==qf2)||(count==qf3)||(count==qf4))
		QUERY_SEQUENCE=0;		/*there is no query sequence*/
	else {	
		rewind(fp);
		/* read lines until we break */
		while(fgets(line, ELLINE, fp)) {
			if(strlen(line) < 2)
				continue;			/*ignore this line*/
			if(line[0]=='>') {		/*gene ID line*/
				fgets(line, ELLINE, fp);
				/* If two ID lines in a row with no sequence in between */
				if(line[0]=='>')
					p_exit("input query format missing DNA sequence!\n");
				/* Require all sequences to be in FASTA format */
				for(c=0; c < strlen(line)-1;  c++) {
					if( (toupper(line[c])!='A') &&
						(toupper(line[c])!='C') &&
						(toupper(line[c])!='G') &&
						(toupper(line[c])!='T') &&
						(toupper(line[c])!='U') &&
						(toupper(line[c])!='R') &&
						(toupper(line[c])!='Y') &&
						(toupper(line[c])!='K') &&
						(toupper(line[c])!='M') &&
						(toupper(line[c])!='S') &&
						(toupper(line[c])!='W') &&
						(toupper(line[c])!='B') &&
						(toupper(line[c])!='D') &&
						(toupper(line[c])!='H') &&
						(toupper(line[c])!='V') &&
						(toupper(line[c])!='N')	&&
						(line[c]!='-')&&(line[c]!='\r'))
						p_exit("query input fasta character error!\n");	
				}
			} else {
			/* Require all sequences to be in FASTA format */
				for(c=0; c < strlen(line)-1;  c++) {
					if( (toupper(line[c])!='A') &&
						(toupper(line[c])!='C') &&
						(toupper(line[c])!='G') &&
						(toupper(line[c])!='T') &&
						(toupper(line[c])!='U') &&
						(toupper(line[c])!='R') &&
						(toupper(line[c])!='Y') &&
						(toupper(line[c])!='K') &&
						(toupper(line[c])!='M') &&
						(toupper(line[c])!='S') &&
						(toupper(line[c])!='W') &&
						(toupper(line[c])!='B') &&
						(toupper(line[c])!='D') &&
						(toupper(line[c])!='H') &&
						(toupper(line[c])!='V') &&
						(toupper(line[c])!='N')	&&
						(line[c]!='-')&&(line[c]!='\r'))
						p_exit("query input fasta character error!\n");	
				}
			}
		}
		QUERY_SEQUENCE=1;	/*sequence present for all id*/
	}
	if(dlevel >3)
		printf("count=%d qf1=%d qf2=%d qf3=%d qf4=%d\n",count,qf1,qf2,qf3,qf4);
	QUERY_FORMAT=1;		/*default*/
	if(qf1==qf2)	
		QUERY_FORMAT=2;	
	if(qf1==qf3)
		QUERY_FORMAT=3;
	if((qf1==qf3)&&(qf3==qf4))
		QUERY_FORMAT=4;	
			
	printf("%s read query input %s to get its format...done!\n", T_S(), fn);
	if(dlevel >2) {
		printf("%s QUERY_FORMAT=%d\n", T_S(), QUERY_FORMAT);
		printf("%s QUERY_SEQUENCE=%d\n", T_S(), QUERY_SEQUENCE);
	}
						
	fclose(fp);
	return PASS;
}
