/* 
 * This module contains all functions dealing with input (both 
 * commandline as well as file input, such as database files
 *	See input.h for more information
 */

#include "lib.h"
#include "defs.h"
#include "input.h"
#include "util.h"


/* initialize the parameters, most of which can be found in lib.h */
short init_param(void) {
	
	/*default parameters value from config file */
	HBRDN_DISPLAY_COUNT = 0;

	/*default count of top quality primer pair display*/
	PRIMER_DISPLAY_COUNT = 10;
	
	/*default number of primers designed by Primer3*/
	PNR = 20;
	
	/*default blast expected value*/
	strcpy(PSPD_EXPECT, "1e-15");
	
	BLAST_OGP = 4;
	BLAST_OIP = 2;

	QUERY_FORMAT = INVALID;				/*initialization*/
	QUERY_SEQUENCE = INVALID;			/*initialization*/
	DATABASE_TYPE = INVALID;	  		/*initialization*/
	ALGORITHM_TYPE = INVALID;			/*initialization*/
	PRIMER_DISPLAY_COUNT = INVALID;		/*initialization*/
	fspd_length = INVALID;				/*initialization*/
	fspd_period = INVALID;				/*initialization*/
	fspd_overlap = INVALID;				/*initialization*/
	CHRCOUNT = INVALID;					/*initialization*/
	BLASTRECORDS = INVALID;				/*initialization*/
	BLASTLOOP = INVALID;				/*initialization*/
	CHRDBINLOOP = INVALID;				/*initialization*/
	AMPSB = INVALID;					/*initialization*/
	MINBIND = INVALID;					/*initialization*/
		
	PEF=0;              /*PRIMER_EXPLAIN_FLAG=1*/
	PFF=0;              /*PRIMER_FILE_FLAG=0*/
	PNR=25;             /*PRIMER_NUM_RETURN=20*/
	PMPX=0;             /*PRIMER_MAX_POLY_X=0*/
	PIOMPX=0;           /*PRIMER_INTERNAL_OLIGO_MAX_POLY_X=0*/
	PMINT=56.0;         /*PRIMER_MIN_TM=50.0*/
	POT=60.0;           /*PRIMER_OPT_TM=58.0*/
	PMAXT=66.0;         /*PRIMER_MAX_TM=60.0*/
	PMINS=19;           /*PRIMER_MIN_SIZE=20*/
	POS=20;             /*PRIMER_OPT_SIZE=24*/
	PMAXS=23;           /*PRIMER_MAX_SIZE=28*/
	PMAXGC=80;          /*PRIMER_MAX_GC=80*/
	POPTGC=50;          /*PRIMER_OPT_GC=50*/
	PMINGC=20;          /*PRIMER_MIN_GC=20*/
	PSE=3;              /*PRIMER_SELF_END=100*/
	PSA=8;              /*PRIMER_SELF_ANY=100*/
	PMES=100.0;         /*PRIMER_MAX_END_STABILITY=100.0*/

        strcpy(PPSR, "250-450");                  /*PRIMER_PRODUCT_SIZE_RANGE=250-450*/
        strcpy(TARGET, "unassigned");                /*target region to specify*/
        strcpy(PIR, "unassigned");                   /*INCLUDED_REGION*/
        strcpy(PER, "unassigned");                   /*EXCLUDED_REGION*/
        strcpy(PSQ, "unassigned");                   /*PRIMER_SEQUENCE_QUALITY*/
        PSCP=-1000000;                         /*PRIMER_START_CODON_POSITION=-1000000*/
        PPA=0;                                 /*PRIMER_PICK_ANYWAY=0*/

        PPMAXT=1000000.0;                      /*PRIMER_PRODUCT_MAX_TM=1000000.0*/
        PPMINT=-1000000.0;                     /*PRIMER_PRODUCT_MIN_TM=-1000000.0*/
        PGC=0;                                 /*PRIMER_GC_CLAMP=0*/
        PDS=20;                                /*PRIMER_DEFAULT_SIZE=20*/
        PMDT=100.0;                            /*PRIMER_MAX_DIFF_TM=100.0C*/
        PTS=1;                                 /*PRIMER_TM_SANTALUCIA=0*/
        PSC=50.0;                              /*PRIMER_SALT_CONC=50.0 mM*/
        PDIVC=0.0;                             /*PRIMER_DIVALENT_CONC=0.0 mM*/
        PDNTPC=0.0;                            /*PRIMER_DNTP_CONC=0.0 mM*/
        PSALTC=1;                              /*PRIMER_SALT_CORRECTIONS=0*/
        PLM=0;                                 /*PRIMER_LOWERCASE_MASKING=0*/
        PDNAC=50.0;                            /*PRIMER_DNA_CONC=50.0nM*/
        PNNSA=0;                               /*PRIMER_NUM_NS_ACCEPTED=0*/
        PLB=0;                                 /*PRIMER_LIBERAL_BASE=0*/
        PFBI=0;                                /*PRIMER_FIRST_BASE_INDEX=0*/
        PMQ=0;                                 /*PRIMER_MIN_QUALITY=0*/
        PMEQ=0;                                /*PRIMER_MIN_END_QUALITY=0*/
        PQRMIN=0;                              /*PRIMER_QUALITY_RANGE_MIN=0*/
        PQRMAX=100;                            /*PRIMER_QUALITY_RANGE_MAX=100*/
        PIP=-1.0;                              /*PRIMER_INSIDE_PENALTY=-1.0*/
        POP=0.0;                               /*PRIMER_OUTSIDE_PENALTY=0.0*/
        PPOTM=0.0;                             /*PRIMER_PRODUCT_OPT_TM=0.0*/
        PPOS=0;                                /*PRIMER_PRODUCT_OPT_SIZE=0*/

        PWTG=1.0;                              /*PRIMER_WT_TM_GT=1.0*/
        PWTL=1.0;                              /*PRIMER_WT_TM_LT=1.0*/
        PWSL=1.0;                              /*PRIMER_WT_SIZE_LT=1.0*/
        PWSG=1.0;                              /*PRIMER_WT_SIZE_GT=1.0*/
        PWGPL=1.0;                             /*PRIMER_WT_GC_PERCENT_LT=1.0*/
        PWGPG=1.0;                             /*PRIMER_WT_GC_PERCENT_GT=1.0*/
        PWCA=0.0;                              /*PRIMER_WT_COMPL_ANY=0.0*/
        PWCE=0.0;                              /*PRIMER_WT_COMPL_END=0.0*/
        PWNNS=0.0;                             /*PRIMER_WT_NUM_NS=0.0*/
        PWRS=0.0;                              /*PRIMER_WT_REP_SIM=0.0*/
        PWSQ=0.0;                              /*PRIMER_WT_SEQ_QUAL=0.0*/
        PWEQ=0.0;                              /*PRIMER_WT_END_QUAL=0.0*/
        PWPP=0.0;                              /*PRIMER_WT_POS_PENALTY=0.0*/
        PWES=0.0;                              /*PRIMER_WT_END_STABILITY=0.0*/
        PPWPP=1.0;                             /*PRIMER_PAIR_WT_PR_PENALTY=1.0*/
        PPWDTM=0.0;                            /*PRIMER_PAIR_WT_DIFF_TM=0.0*/
        PPWCA=0.0;                             /*PRIMER_PAIR_WT_COMPL_ANY=0.0*/
        PPWCE=0.0;                             /*PRIMER_PAIR_WT_COMPL_END=0.0*/
        PPWPTL=0.0;                            /*PRIMER_PAIR_WT_PRODUCT_TM_LT=0.0*/
        PPWPTG=0.0;                            /*PRIMER_PAIR_WT_PRODUCT_TM_GT=0.0*/
        PPWPSL=0.0;                            /*PRIMER_PAIR_WT_PRODUCT_SIZE_LT=0.0*/
        PPWPSG=0.0;                            /*PRIMER_PAIR_WT_PRODUCT_SIZE_GT=0.0*/
        PPWRS=0.0;                             /*PRIMER_PAIR_WT_REP_SIM=0.0*/

	/*probe design only default setting*/
	PROBE_DESIGN=0;
	PROBE_MAX_GC=80.0;
	PROBE_MIN_GC=40.0;
	PROBE_BLAST_HIT=0;
	PROBE_LENGTH=150;
	PROBE_PERIOD=20;
	PROBE_REGION=500;
	PROBE_GAP=60;
	PROBE_UNIQUE=1;

	/*other variables*/
	CUT_SITE_COUNT=0;
	CUT_SITE_REGION=0;

	return PASS;
}

/* Deal with command line parameters */
short get_param (int argc, char * argv[]) {
	short	i=0;	/*q; d/l; p; o options*/
	char *tmp;		/*used temporarity for storing string*/
	char s_fn_log[SLINE], s_fn_status[SLINE];
	FILE *fp;
	const char *hmsg = 
	"\n\n\t\t\tPRIMEGNENS help\n"
	"\t $>primegens.exe ...\n"
	"\t -q: query sequence list filename\n"
	"\t -d/l: database/list filename\n"
	"\t -p: database location path\n"
	"\t -o: output result location\n"
	"\t -h: help\n\n" ;

	/* too few arguments */
	if(argc < 9) {
		fprintf(stderr, hmsg);
		fprintf(stderr, "%s arguments count: %d\n", T_S(), argc);
		p_exit("less no. of arguments\n");
	}
	
	/* too many arguments */
	if(argc != 9) {
		printf("%s total arguments: %d\n", T_S(), argc);
		printf("%s $exe -q <query_input_filename>", T_S());
		printf(" -d <database_filename> -p <database_path>\n");
		printf("%s\t\tOR\n", T_S());
		printf("%s $exe -q <query_input_filename>", T_S());
		printf(" -d <database_list_filename> -p <database_path>\n");
		p_exit(" wrong command line input\n");
	}
	
	/* exactly the right number of arguments! */
	
	/* ensure there are only '-q', '-d', '-l', '-p', '-o' arguments */
	for(i=0; i < 4; i++)
		argf[i] = 0;
	for(i=1; i < argc; i++) { 
		if(strcmp(argv[i], "-q")==0) {
			argf[0] = 1;						/*mandatory option*/
			strcpy(fn_query, argv[++i]);
			printf("%s query input: %s\n", T_S(), fn_query);
		} else if(strcmp(argv[i], "-d")==0) {
			argf[1] = 1;						/*mandatory option*/
			#ifdef _WIN32
				tmp = strrchr(argv[++i], '\\');
			#else 
				tmp = strrchr(argv[++i], '/');
			#endif
			if(tmp==NULL)
				strcpy(fn_db0.fn, argv[i]);
			else {
				tmp ++;
				strcpy(fn_db0.fn, tmp);
			}
			printf("%s database input: %s\n", T_S(), fn_db0.fn);
		} else if(strcmp(argv[i], "-l")==0) {
			argf[1] = 2;						/*mandatory option*/
			strcpy(fn_dbl, argv[++i]);
			printf("%s database list input: %s\n",T_S(),fn_dbl);
		} else if(strcmp(argv[i], "-p")==0) {
			argf[2] = 1;						/*mandatory option*/
			strcpy(dbpath, argv[++i]);
			printf("%s database location: %s\n", T_S(), dbpath);
		} else if(strcmp(argv[i], "-o")==0) {
			argf[3] = 1;						/*mandatory option*/
			strcpy(out_path, argv[++i]);
			printf("%s output location: %s\n", T_S(), out_path);
		} else { 
			fprintf(stderr, hmsg);
			printf("%s argv[%d]=%s\n", T_S(), i, argv[i]);
			p_exit("wrong argument option\n");
		}
	}
	
	/* missing an argument */
	if(argf[0]*argf[1]*argf[2]*argf[3] <= 0) 
		p_exit("some problem in command. check!\n");
    
	strcpy(s_fn_log, "prime design log filename");
	strcpy(s_fn_status, "prime design status filename");
	#ifdef _WIN32
		sprintf(fn_log, "%s\\%s", out_path, LOGFILENAME);
		sprintf(fn_status, "%s\\%s", out_path, STATUSFILE);
	#else
		sprintf(fn_log, "%s/%s", out_path, LOGFILENAME);
		sprintf(fn_status, "%s/%s", out_path, STATUSFILE);
	#endif
	printf("%s %-40s: %s\n", T_S(), s_fn_log, fn_log);
	printf("%s %-40s: %s\n", T_S(), s_fn_status, fn_status);
	fp = p_fopen_error_exit(fn_status, "w");
	fprintf(fp, "PRIMEGENSv2 is still running...\n");
	fclose(fp);
	return PASS;
}

/* Set various parameters to their defaul values */
short set_default_param(void) {
    int tmp=0;
    char cmd[COMMAND];
	char s_blast[SLINE];
    char s_megablast[SLINE];
    char s_formatdb[SLINE];
    char s_primer3[SLINE];
    char s_env_path[SLINE];
	char s_fn_fmtdblog[SLINE];
	char s_fn_vdb[SLINE];
    char s_fn_p3in[SLINE];
    char s_fn_p3out[SLINE];
    char s_fn_uoligo[SLINE];
    char s_fn_bin[SLINE];
    char s_fn_bout[SLINE];
	char s_fn_query[SLINE];
	FILE *fout1, *fp;
	
    strcpy(s_blast, "blast command");
    strcpy(s_megablast, "megablast command");
    strcpy(s_formatdb, "formatdb command");
    strcpy(s_primer3, "primer3 command");
    strcpy(s_env_path, "environment path");
	strcpy(s_fn_fmtdblog, "formatdb log filename");
    strcpy(s_fn_vdb, "virtual database filename");
    strcpy(s_fn_p3in, "primer3 input filename");
    strcpy(s_fn_p3out, "primer3 output filename");
    strcpy(s_fn_uoligo, "unique oligo list filename");
    strcpy(s_fn_bin, "blast intput filename");
    strcpy(s_fn_bout, "blast output filename");
    
	/* Deal with Chromosome info */
	if(DATABASE_TYPE==3) {
		if(CHRCOUNT==INVALID)
			p_exit("database type genome with no chromosome?\n");
		if(BLASTRECORDS==INVALID)
			BLASTRECORDS = CHRCOUNT;
		if(BLASTRECORDS==24) {
			BLASTLOOP = 4;
			CHRDBINLOOP = 6;
		} else if(BLASTRECORDS==48) {
			BLASTLOOP = 4;
			CHRDBINLOOP = 12;
		} 
		if((BLASTRECORDS==CHRCOUNT)&&(BLASTLOOP==-1)||(CHRDBINLOOP==-1)){
			BLASTLOOP = CHRCOUNT;
			CHRDBINLOOP = 1;
		}
	}
	
	/*some PSPD options default setting */
	PSPD_design_primer_flag=1;  
	PSPD_BLAST_OGP=2;
	PSPD_BLAST_OIP=2;
	PSPD_MINFILTERLENGTH=50;
	PSPD_RUN_SIM_FLAG=1;
	PSPD_MAX_SIMILARITY=0.75;
	PSPD_MINSEGLENGTH=100;
	strcpy(PSPD_EXPECT, "1e-15");
	
	#ifdef _WIN32
		strcpy(FORMATDBEXEC, "formatdb.exe");
		strcpy(PRIMERDESIGNEXEC, "primer3.exe");
		strcpy(BLASTEXEC, "megablast.exe");
		strcpy(BLASTALLEXEC, "blastall.exe");
		strcpy(ALIGNMENTEXEC, "clustalw.exe");
	#else
		strcpy(FORMATDBEXEC, "formatdb");
		strcpy(PRIMERDESIGNEXEC, "primer3_core");
		strcpy(BLASTEXEC, "megablast");
		strcpy(BLASTALLEXEC, "blastall");
		strcpy(ALIGNMENTEXEC, "clustalw");
	#endif
	
	
	if(MINBIND > PMINS)		/*min primer match > min primer size*/
		p_exit("minbind is greater than primer min size\n");

	if((DATABASE_TYPE==3)&&(BLASTRECORDS != CHRCOUNT))
		p_exit("blast records != chromosome count\n");
	
	tmp = BLASTLOOP*CHRDBINLOOP;
	if((DATABASE_TYPE==3)&&(BLASTRECORDS != tmp))
		p_exit("blast records != bloop*chr_in_loop\n");

	printf("%s setting default parameters...\n", T_S());
		
	#ifdef _WIN32
		/*env_path*/
        sprintf(formatdb, "%s\\blast\\%s", env_path, FORMATDBEXEC);
        sprintf(blast, "%s\\blast\\%s", env_path, BLASTEXEC);
        sprintf(blastall, "%s\\blast\\%s", env_path, BLASTALLEXEC);
        sprintf(megablast, "%s\\blast\\%s", env_path, BLASTEXEC);
        sprintf(primer3, "%s\\primer3\\%s", env_path, PRIMERDESIGNEXEC);
		sprintf(tmp_dir, "%s\\tmpDir_%d\\", tmpDir_path, pid);
		/*work_dir*/
		sprintf(fn_vdb, "%s\\%s", tmp_dir, VDATABASE);
		sprintf(fn_p3in, "%s\\%s", tmp_dir, PRIMER3IN);
		sprintf(fn_p3out, "%s\\%s", tmp_dir, PRIMER3OUT);
		sprintf(fn_uoligo, "%s\\%s", tmp_dir, UOLIGOTMP);
        sprintf(fn_bin, "%s\\%s", tmp_dir, BLASTIN);
        sprintf(fn_bout, "%s\\%s", tmp_dir, BLASTOUT);
		/*out_path*/
		if(strrchr(fn_query, '\\')!=NULL)
			strcpy(s_fn_query, strrchr(fn_query, '\\'));
		else 
			strcpy(s_fn_query, fn_query);
		sprintf(fn_fmtdblog, "%s\\%s", out_path, FORMATDBLOG);
    	sprintf(fn_out1, "%s\\%s%s", out_path, s_fn_query, EXCELSHEET);
    	sprintf(fn_out2, "%s\\%s%s", out_path, s_fn_query, PRIMERLIST);
    	sprintf(fn_out3, "%s\\%s%s", out_path, s_fn_query, QUERYPRIMER);
    	sprintf(fn_out4, "%s\\%s%s", out_path, s_fn_query, QUERYFAILED);
    	sprintf(fn_out5, "%s\\%s%s", out_path, s_fn_query, QUERYFASTA);
		sprintf(fn_nohit, "%s\\%s%s", out_path, s_fn_query, QUERYNOHIT);
		sprintf(fn_uniseg, "%s\\%s%s", out_path, s_fn_query, QUERYUNISEG);
		sprintf(fn_sim, "%s\\%s%s", out_path, s_fn_query, QUERYSIMFILE);
		sprintf(fn_probe, "%s\\%s%s", out_path, s_fn_query, QUERYPROBEA);
		sprintf(fn_probe_hit, "%s\\%s%s", out_path, s_fn_query, QUERYPROBEB);
    #else
		/*env_path*/
        sprintf(formatdb, "%s/blast/%s", env_path, FORMATDBEXEC);
        sprintf(blast, "%s/blast/%s", env_path, BLASTEXEC);
        sprintf(blastall, "%s/blast/%s", env_path, BLASTALLEXEC);
        sprintf(megablast, "%s/blast/%s", env_path, BLASTEXEC);
        sprintf(primer3, "%s/primer3/%s", env_path, PRIMERDESIGNEXEC);
		sprintf(tmp_dir, "%s/tmpDir_%d/", tmpDir_path, pid);
        /*work_dir*/
        sprintf(fn_vdb, "%s/%s", tmp_dir, VDATABASE);
		sprintf(fn_p3in, "%s/%s", tmp_dir, PRIMER3IN);
		sprintf(fn_p3out, "%s/%s", tmp_dir, PRIMER3OUT);
		sprintf(fn_uoligo, "%s/%s", tmp_dir, UOLIGOTMP);
        sprintf(fn_bin, "%s/%s", tmp_dir, BLASTIN);
        sprintf(fn_bout, "%s/%s", tmp_dir, BLASTOUT);
		/*out_path*/
		if(strrchr(fn_query, '/')!=NULL)
			strcpy(s_fn_query, strrchr(fn_query, '/'));
		else 
			strcpy(s_fn_query, fn_query);
		sprintf(fn_fmtdblog, "%s/%s", out_path, FORMATDBLOG);
    	sprintf(fn_out1, "%s/%s%s", out_path, s_fn_query, EXCELSHEET);
    	sprintf(fn_out2, "%s/%s%s", out_path, s_fn_query, PRIMERLIST);
    	sprintf(fn_out3, "%s/%s%s", out_path, s_fn_query, QUERYPRIMER);
    	sprintf(fn_out4, "%s/%s%s", out_path, s_fn_query, QUERYFAILED);
    	sprintf(fn_out5, "%s/%s%s", out_path, s_fn_query, QUERYFASTA);
		sprintf(fn_nohit, "%s/%s%s", out_path, s_fn_query, QUERYNOHIT);
		sprintf(fn_uniseg, "%s/%s%s", out_path, s_fn_query, QUERYUNISEG);
		sprintf(fn_sim, "%s/%s%s", out_path, s_fn_query, QUERYSIMFILE);
		sprintf(fn_probe, "%s/%s%s", out_path, s_fn_query, QUERYPROBEA);
		sprintf(fn_probe_hit, "%s/%s%s", out_path, s_fn_query, QUERYPROBEB);
    #endif

    printf("%s commands & temporary file names\n\n", T_S());
    printf("%s %-40s: %s\n", T_S(), s_blast, blast); 
    printf("%s %-40s: %s\n", T_S(), s_megablast, megablast); 
    printf("%s %-40s: %s\n", T_S(), s_formatdb, formatdb);
    printf("%s %-40s: %s\n", T_S(), s_primer3, primer3);
    printf("%s %-40s: %s\n", T_S(), s_env_path, env_path);
    printf("%s %-40s: %s\n", T_S(), s_fn_fmtdblog, fn_fmtdblog);
    printf("%s %-40s: %s\n", T_S(), s_fn_vdb, fn_vdb);
    printf("%s %-40s: %s\n", T_S(), s_fn_p3in, fn_p3in);
    printf("%s %-40s: %s\n", T_S(), s_fn_p3out, fn_p3out);
    printf("%s %-40s: %s\n", T_S(), s_fn_uoligo, fn_uoligo);
    printf("%s %-40s: %s\n", T_S(), s_fn_bin, fn_bin);
    printf("%s %-40s: %s\n", T_S(), s_fn_bout, fn_bout);

	/*fragment-specific primer design parameters*/
	#ifdef _WIN32
		sprintf(cmd, "del %s", fn_out1);
    	p_exec (cmd);
    	sprintf(cmd, "del %s", fn_out2);
    	p_exec (cmd);
    	sprintf(cmd, "del %s", fn_out3);
    	p_exec (cmd);
    	sprintf(cmd, "del %s", fn_out4);
    	p_exec (cmd);
    	sprintf(cmd, "del %s", fn_out5);
    	p_exec (cmd);
    	sprintf(cmd, "del %s", fn_nohit);
    	p_exec (cmd);
    	sprintf(cmd, "del %s", fn_uniseg);
    	p_exec (cmd);

	#else
		sprintf(cmd, "rm %s", fn_out1);
    	p_exec (cmd);
    	sprintf(cmd, "rm %s", fn_out2);
    	p_exec (cmd);
    	sprintf(cmd, "rm %s", fn_out3);
    	p_exec (cmd);
    	sprintf(cmd, "rm %s", fn_out4);
    	p_exec (cmd);
    	sprintf(cmd, "rm %s", fn_out5);
    	p_exec (cmd);
    	sprintf(cmd, "rm %s", fn_nohit);
    	p_exec (cmd);
    	sprintf(cmd, "rm %s", fn_uniseg);
    	p_exec (cmd);
	#endif

	if(PROBE_DESIGN!=1) {
		fout1 = p_fopen_error_exit(fn_out1, "w");
    	fprintf(fout1, "QUERY_NAME\t");
        if (CUT_SITE_COUNT>0){
    	    fprintf(fout1, "MAX_CUTSITE_REGION\t");
    	    fprintf(fout1, "# of CUT_SITE\t");
    	    fprintf(fout1, "CUT_SITE\t");
        }
    	fprintf(fout1, "LEFT_PRIMER\t");
    	fprintf(fout1, "LEFT_PRIMER_START_POSITION\t");
    	fprintf(fout1, "LEFT_PRIMER_LENGTH\t");
    	fprintf(fout1, "LEFT_PRIMER_TM\t");
    	fprintf(fout1, "LEFT_PRIMER_GC_CONTENT\t");
    	fprintf(fout1, "RIGHT_PRIMER\t");    
    	fprintf(fout1, "RIGHT_PRIMER_START_POSITION\t");
    	fprintf(fout1, "RIGHT_PRIMER_LENGTH\t");
    	fprintf(fout1, "RIGHT_PRIMER_TM\t");
    	fprintf(fout1, "RIGHT_PRIMER_GC_CONTENT\t");
    	fprintf(fout1, "PRODUCT_SIZE\t");
    	fprintf(fout1, "HYBRIDIZATION\n");
		fclose(fout1);
	}
	
	
	/*dumping into log file*/
	fp = p_fopen_error_exit(fn_log, "w");
	fprintf(fp, "\n\t\t\t Primer Design Report\n");
	fprintf(fp, "\t\t\t\t____________________\n\n");
	fprintf(fp, "\t\t\t\tPrimer Design Start Time:%s\n", T_S());
	fprintf(fp, "___________________________________________________\n");
	fprintf(fp, "Input data\n");
	fprintf(fp, "%s\t query input: %s\n", T_S(), fn_query);
	fprintf(fp, "%s\t database input: %s\n", T_S(), fn_db0.fn);
	fprintf(fp, "%s\t database location: %s\n", T_S(), dbpath);
	fprintf(fp, "%s\t output location: %s\n", T_S(), out_path);
	fprintf(fp, "___________________________________________________\n");
	fprintf(fp, "Primer design configurations\n");
	fprintf(fp, "\t query format : %d\n", QUERY_FORMAT);
	fprintf(fp, "\t query sequence : %d \n", QUERY_SEQUENCE);
	fprintf(fp, "\t database type : %d \n", DATABASE_TYPE);
	fprintf(fp, "\t probe design: %d\n", PROBE_DESIGN);
	fprintf(fp, "\t probe length: %d\n", PROBE_LENGTH);
	fprintf(fp, "\t probe period: %d\n", PROBE_PERIOD);
	fprintf(fp, "\t probe region: %d\n", PROBE_REGION);
	fprintf(fp, "\t probe gap: %d\n", PROBE_GAP);
	fprintf(fp, "\t max 'cg' in primer : %d \n", PRIMER_CG_MAX);
	fprintf(fp, "\t algorithm type : %d \n", ALGORITHM_TYPE);
	fprintf(fp, "\t Probe design : %d\n", PROBE_DESIGN);
	fprintf(fp, "\t Probe min. GC-content: %f\n", PROBE_MIN_GC);
	fprintf(fp, "\t Probe max. GC-content: %f\n", PROBE_MAX_GC);
	fprintf(fp, "\t probe uniqueness: %d\n", PROBE_UNIQUE);
	fprintf(fp, "\t fspd length : %d \n", fspd_length);
	fprintf(fp, "\t fspd overlap : %d \n", fspd_overlap);
	fprintf(fp, "\t chromosom count : %d \n", CHRCOUNT);
	fprintf(fp, "\t blast records : %d \n", BLASTRECORDS);
	fprintf(fp, "\t blast loop count : %d \n", BLASTLOOP);
	fprintf(fp, "\t chrosome in one loop : %d \n", CHRDBINLOOP);
	fprintf(fp, "\t smallest possible amplicon : %d\n", AMPSB);
	fprintf(fp, "\t minimum oligo binding from right : %d\n", MINBIND);
	fprintf(fp, "\t formatdb program : %s\n", FORMATDBEXEC);
	fprintf(fp, "\t blast program : %s\n", BLASTEXEC);
	fprintf(fp, "\t word size for blast : %d\n", WORDSIZE);
	fprintf(fp, "\t primer design program : %s\n", PRIMERDESIGNEXEC);
	fprintf(fp, "\t alignment program : %s\n", ALIGNMENTEXEC);
	fprintf(fp, "\t create log file : %d\n", log_flag);
	fprintf(fp, "\t design primer decision : %d \n", PSPD_design_primer_flag);
	fprintf(fp, "\t PRIMER_PRODUCT_SIZE_RANGE : %s\n", PPSR);
	fprintf(fp, "\t PRIMER_EXPLAIN_FLAG : %d \n", PEF);
	fprintf(fp, "\t PRIMER_FILE_FLAG : %d \n", PFF);
	fprintf(fp, "\t PRIMER_NUM_RETURN : %d \n", PNR);
	fprintf(fp, "\t PRIMER_MAX_POLY_X : %d \n", PMPX);
	fprintf(fp, "\t PRIMER_INTERNAL_OLIGO_MAX_POLY_X : %d \n", PIOMPX);
	fprintf(fp, "\t PRIMER_MIN_TM : %f \n", PMINT);
	fprintf(fp, "\t PRIMER_OPT_TM : %f \n", POT);
	fprintf(fp, "\t PRIMER_MAX_TM : %f \n", PMAXT);
	fprintf(fp, "\t PRIMER_MIN_SIZE : %d \n", PMINS);
	fprintf(fp, "\t PRIMER_OPT_SIZE : %d \n", POS);
	fprintf(fp, "\t PRIMER_MAX_SIZE : %d \n", PMAXS);
	fprintf(fp, "\t PRIMER_MAX_GC : %d \n", PMAXGC);
	fprintf(fp, "\t PRIMER_OPT_GC : %d \n", POPTGC);
	fprintf(fp, "\t PRIMER_MIN_GC : %d \n", PMINGC);
	fprintf(fp, "\t PRIMER_SELF_END : %d \n", PSE);
	fprintf(fp, "\t PRIMER_SELF_ANY : %d \n", PSA);
	fprintf(fp, "\t PRIMER_MAX_END_STABILITY %f \n", PMES);
	fprintf(fp, "\t HBRDN_DISPLAY_COUNT: %d\n",HBRDN_DISPLAY_COUNT);
       
        fprintf(fp, "\t PRIMER_PRODUCT_MAX_TM : %f \n", PPMAXT);
        fprintf(fp, "\t PRIMER_PRODUCT_MIN_TM : %f \n", PPMINT);
        fprintf(fp, "\t PRIMER_GC_CLAMP : %d \n", PGC);
        fprintf(fp, "\t PRIMER_DEFAULT_SIZE : %d \n", PDS);
        fprintf(fp, "\t PRIMER_MAX_DIFF_TM : %f \n", PMDT);
        fprintf(fp, "\t PRIMER_TM_SANTALUCIA : %d \n", PTS);
        fprintf(fp, "\t PRIMER_SALT_CONC : %f \n", PSC);
        fprintf(fp, "\t PRIMER_DIVALENT_CONC : %f \n", PDIVC);
        fprintf(fp, "\t PRIMER_DNTP_CONC : %f \n", PDNTPC);
        fprintf(fp, "\t PRIMER_SALT_CORRECTIONS : %d \n", PSALTC);
        fprintf(fp, "\t PRIMER_LOWERCASE_MASKING : %d \n", PLM);
        fprintf(fp, "\t PRIMER_DNA_CONC : %f \n", PDNAC);
        fprintf(fp, "\t PRIMER_NUM_NS_ACCEPTED : %d \n", PNNSA);
        fprintf(fp, "\t PRIMER_LIBERAL_BASE : %d \n", PLB);
        fprintf(fp, "\t PRIMER_FIRST_BASE_INDEX : %d \n", PFBI);
        fprintf(fp, "\t PRIMER_MIN_QUALITY : %d \n", PMQ);
        fprintf(fp, "\t PRIMER_MIN_END_QUALITY : %d \n", PMEQ);
        fprintf(fp, "\t PRIMER_QUALITY_RANGE_MIN : %d \n", PQRMIN);
        fprintf(fp, "\t PRIMER_QUALITY_RANGE_MAX : %f \n", PQRMAX);
        fprintf(fp, "\t PRIMER_INSIDE_PENALTY : %f \n", PIP);
        fprintf(fp, "\t PRIMER_OUTSIDE_PENALTY : %f \n", POP);
        fprintf(fp, "\t PRIMER_PRODUCT_OPT_TM : %f \n", PPOTM);
        fprintf(fp, "\t PRIMER_PRODUCT_OPT_SIZE : %d \n", PPOS);

        fprintf(fp, "\t PRIMER_WT_TM_GT : %f \n", PWTG);
        fprintf(fp, "\t PRIMER_WT_TM_LT : %f \n", PWTL);
        fprintf(fp, "\t PRIMER_WT_SIZE_LT : %f \n", PWSL);
        fprintf(fp, "\t PRIMER_WT_SIZE_GT : %f \n", PWSG);
        fprintf(fp, "\t PRIMER_WT_GC_PERCENT_LT : %f \n", PWGPL);
        fprintf(fp, "\t PRIMER_WT_GC_PERCENT_GT : %f \n", PWGPG);
        fprintf(fp, "\t PRIMER_WT_COMPL_ANY : %f \n", PWCA);
        fprintf(fp, "\t PRIMER_WT_COMPL_END : %f \n", PWCE);
        fprintf(fp, "\t PRIMER_WT_NUM_NS : %f \n", PWNNS);
        fprintf(fp, "\t PRIMER_WT_REP_SIM : %f \n", PWRS);
        fprintf(fp, "\t PRIMER_WT_SEQ_QUAL : %f \n", PWSQ);
        fprintf(fp, "\t PRIMER_WT_END_QUAL : %f \n", PWEQ);
        fprintf(fp, "\t PRIMER_WT_POS_PENALTY : %f \n", PWPP);
        fprintf(fp, "\t PRIMER_WT_END_STABILITY : %f \n", PWES);
        fprintf(fp, "\t PRIMER_PAIR_WT_PR_PENALTY : %f \n", PPWPP);
        fprintf(fp, "\t PRIMER_PAIR_WT_DIFF_TM : %f \n", PPWDTM);
        fprintf(fp, "\t PRIMER_PAIR_WT_COMPL_ANY : %f \n", PPWCA);
        fprintf(fp, "\t PRIMER_PAIR_WT_COMPL_END : %f \n", PPWCE);
        fprintf(fp, "\t PRIMER_PAIR_WT_PRODUCT_TM_LT : %f \n", PPWPTL);
        fprintf(fp, "\t PRIMER_PAIR_WT_PRODUCT_TM_GT : %f \n", PPWPTG);
        fprintf(fp, "\t PRIMER_PAIR_WT_PRODUCT_SIZE_LT : %f \n", PPWPSL);
        fprintf(fp, "\t PRIMER_PAIR_WT_PRODUCT_SIZE_GT : %f \n", PPWPSG);
        fprintf(fp, "\t PRIMER_PAIR_WT_REP_SIM : %f \n", PPWRS);

	fprintf(fp, "___________________________________________________\n");
	fclose(fp);

    /*
    printf("%s mblast command %60s\n", T_S(), blast);
    printf("%s formatdb command %60s\n", T_S(), formatdb);
    printf("%s primer3 command %s\n", T_S(), primer3);
    printf("%s primegens env_path %s\n", T_S(), env_path);
    printf("%s virtual database filename: %s\n",T_S(),fn_vdb);
    printf("%s primer3 input filename: %s\n",T_S(),fn_p3in);
    printf("%s primer3 output filename: %s\n",T_S(),fn_p3out);
    printf("%s unique oligo list filename: %s\n",T_S(),fn_uoligo);
    printf("%s blast intput filename: %s\n", T_S(), fn_bin);
    printf("%s blast output filename: %s\n", T_S(), fn_bout);
    */

	/*Primer3 parameters*/
    printf("\n\n");

	return PASS;
}

/* deal with configuration file parameters */
short read_config_file(void) {
	FILE *fp;
	short res, flag=0;
	char line[MLINE], *tok;
	char cmd[COMMAND],fn_cfg[FILENAME];

	#ifdef _WIN32
		sprintf(fn_cfg, "%s\\%s", out_path, FNCFG);
	#else
		sprintf(fn_cfg, "%s/%s", out_path, FNCFG);
	#endif

   	/*first search*/
	if(!(fp = fopen(fn_cfg, "r"))) {
		if(!(fp = fopen(FNCFG, "r"))) {
			#ifdef _WIN32
				sprintf(cmd,"copy %s\\include\\%s %s",env_path,FNCFG,fn_cfg);
			#else 
				sprintf(cmd,"cp %s/include/%s %s",env_path,FNCFG,fn_cfg);
			#endif
			printf("%s copying config file from PRIMEGENS package.\n", T_S());
		} else {
			fclose(fp);
			printf("%s using config file from current directoty.\n", T_S());
			#ifdef _WIN32
				sprintf(cmd,"copy %s %s", FNCFG, fn_cfg);
			#else 
				sprintf(cmd,"cp %s %s", FNCFG, fn_cfg);
			#endif
		}
		system(cmd);		/*because log file is still not defined*/
	} else {
		fclose(fp);
		printf("%s configuration file found in output location!\n", T_S());
	}
	/* read lines and check for parameters */
	printf("%s read configuration (%s)...\n", T_S(), fn_cfg);
	fp = p_fopen_error_exit(fn_cfg, "r");
   	while(fgets(line, MLINE, fp)) {
		if(strlen(line) < 3)
			continue;			/*bad/empty line to read*/
		if(strncmp(line, "#", 1)==0)
			continue;			/*commented line*/
		tok = strtok(line, "=\n");
		strcat(tok, "\0");
		if(strcmp(tok, "QUERY_FORMAT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			QUERY_FORMAT = atoi(tok);
			printf("%s QUERY_FORMAT=%d\n", T_S(), QUERY_FORMAT);
		} else if(strcmp(tok, "QUERY_SEQUENCE")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			QUERY_SEQUENCE = atoi(tok);
			printf("%s QUERY_SEQUENCE=%d\n", T_S(), QUERY_SEQUENCE);
		} else if(strcmp(tok, "DATABASE_TYPE")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			DATABASE_TYPE = atoi(tok);
			printf("%s DATABASE_TYPE=%d\n", T_S(), DATABASE_TYPE);
		} else if(strcmp(tok, "PRIMER_CG_MAX")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			if(strncmp(tok, NA, 1)==0)
				PRIMER_CG_MAX = INVALID;
			else 	
				PRIMER_CG_MAX = atoi(tok);
			printf("%s PRIMER_CG_MAX=%d\n", T_S(), PRIMER_CG_MAX);
		} else if(strcmp(tok, "ALGORITHM_TYPE")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			ALGORITHM_TYPE = atoi(tok);
			printf("%s ALGORITHM_TYPE=%d\n", T_S(), ALGORITHM_TYPE);
		} else if(strcmp(tok, "PROBE_DESIGN")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PROBE_DESIGN = atoi(tok);
			printf("%s PROBE_DESIGN=%d\n", T_S(), PROBE_DESIGN);
		} else if(strcmp(tok, "BISULFITE_MODE")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			BISULFITE_MODE = atoi(tok);
			printf("%s BISULFITE_MODE=%d\n", T_S(), BISULFITE_MODE);
		} else if(strcmp(tok, "PROBE_MAX_GC")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PROBE_MAX_GC = atof(tok);
			printf("%s PROBE_MAX_GC=%f\n", T_S(), PROBE_MAX_GC);
		} else if(strcmp(tok, "PROBE_MIN_GC")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PROBE_MIN_GC = atof(tok);
			printf("%s PROBE_MIN_GC=%f\n", T_S(),  PROBE_MIN_GC);
		} else if(strcmp(tok, "PROBE_BLAST_HIT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PROBE_BLAST_HIT = atoi(tok);
			printf("%s PROBE_BLAST_HIT=%d\n", T_S(), PROBE_BLAST_HIT);
		} else if(strcmp(tok, "PROBE_LENGTH")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PROBE_LENGTH = atoi(tok);
			printf("%s PROBE_LENGTH=%d\n", T_S(), PROBE_LENGTH);
		} else if(strcmp(tok, "PROBE_PERIOD")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PROBE_PERIOD = atoi(tok);
			printf("%s PROBE_PERIOD=%d\n", T_S(), PROBE_PERIOD);
		} else if(strcmp(tok, "PROBE_REGION")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PROBE_REGION = atoi(tok);
			printf("%s PROBE_REGION=%d\n", T_S(), PROBE_REGION);
		} else if(strcmp(tok, "PROBE_UNIQUE")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PROBE_UNIQUE = atoi(tok);
			printf("%s PROBE_UNIQUE=%d\n", T_S(),  PROBE_UNIQUE);
		} else if(strcmp(tok, "PROBE_GAP")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PROBE_GAP = atoi(tok);
			printf("%s PROBE_GAP=%d\n", T_S(),  PROBE_GAP);
		} else if(strcmp(tok, "DISPLAY_LEVEL")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			dlevel = atoi(tok);
			printf("%s DISPLAY_LEVEL=%d\n", T_S(), dlevel);
		} else if(strcmp(tok, "HBRDN_DISPLAY_COUNT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			HBRDN_DISPLAY_COUNT = atoi(tok);
			printf("%s HBRDN_DISPLAY_COUNT=%d\n",T_S(),HBRDN_DISPLAY_COUNT);
		} else if(strcmp(tok, "PRIMER_DISPLAY_COUNT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PRIMER_DISPLAY_COUNT = atoi(tok);
			printf("%s PRIMER_DISPLAY_COUNT=%d\n",T_S(),PRIMER_DISPLAY_COUNT);
		} else if(strcmp(tok, "FSPD_LENGTH")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			fspd_length = atoi(tok);
			printf("%s FSPD_LENGTH=%d\n", T_S(), fspd_length);
		} else if(strcmp(tok, "FSPD_PERIOD")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			fspd_period = atoi(tok);
			printf("%s FSPD_PERIOD=%d\n", T_S(), fspd_period);
		} else if(strcmp(tok, "FSPD_OVERLAP")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			fspd_overlap = atoi(tok);
			printf("%s FSPD_OVERLAP=%d\n", T_S(), fspd_overlap);
		} else if(strcmp(tok, "CUT_SITE_COUNT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			CUT_SITE_COUNT = atoi(tok);
			printf("%s CUT_SITE_COUNT=%d\n", T_S(), CUT_SITE_COUNT);
		} else if(strcmp(tok, "CUT_SITE_REGION")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			CUT_SITE_REGION = atoi(tok);
			printf("%s CUT_SITE_REGION=%d\n", T_S(), CUT_SITE_REGION);
		} else if(strcmp(tok, "CUT_SITE")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			strcpy(CUT_SITE, tok);
			printf("%s CUT_SITE=%s\n", T_S(), CUT_SITE);
		} else if(strcmp(tok, "CHRCOUNT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			CHRCOUNT = atoi(tok);
			printf("%s CHRCOUNT=%d\n", T_S(), CHRCOUNT);
		} else if(strcmp(tok, "BLASTRECORDS")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			BLASTRECORDS = atoi(tok);
			printf("%s  BLASTRECORDS=%d\n", T_S(),  BLASTRECORDS);
		} else if(strcmp(tok, "BLASTLOOP")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			BLASTLOOP = atoi(tok);
			printf("%s BLASTLOOP=%d\n", T_S(), BLASTLOOP);
		} else if(strcmp(tok, "CHRDBINLOOP")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			CHRDBINLOOP = atoi(tok);
			printf("%s CHRDBINLOOP=%d\n", T_S(), CHRDBINLOOP);
		} else if(strcmp(tok, "AMPSB")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			AMPSB = atoi(tok);
			printf("%s AMPSB=%d\n", T_S(), AMPSB);
		} else if(strcmp(tok, "MINBIND")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			MINBIND = atoi(tok);
			printf("%s MINBIND=%d\n", T_S(), MINBIND);
		} else if(strcmp(tok, "FORMATDB_PROGRAM")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			strcpy(FORMATDBEXEC, tok);
			printf("%s FORMATDB_PROGRAM=%s\n", T_S(), FORMATDBEXEC);
		} else if(strcmp(tok, "BLAST_PROGRAM")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			strcpy(BLASTEXEC, tok);
			printf("%s BLAST_PROGRAM=%s\n", T_S(), BLASTEXEC);
		} else if(strcmp(tok, "WORDSIZE")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			WORDSIZE = atoi(tok);
			printf("%s WORDSIZE=%d\n", T_S(), WORDSIZE);
		} else if(strcmp(tok, "PRIMERDESIGN_PROGRAM")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			strcpy(PRIMERDESIGNEXEC, tok);
			printf("%s PRIMERDESIGN_PROGRAM=%s\n",T_S(),PRIMERDESIGNEXEC);
		} else if(strcmp(tok, "ALIGNMENT_PROGRAM")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			strcpy(ALIGNMENTEXEC, tok);
			printf("%s ALIGNMENT_PROGRAM=%s\n",T_S(),ALIGNMENTEXEC);
		} else if(strcmp(tok, "LOGS")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			log_flag = atoi(tok);
			printf("%s LOGS=%d\n", T_S(), log_flag);
		} else if(strcmp(tok, "PSPD_DESIGN_PRIMER")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PSPD_design_primer_flag = atoi(tok);
			printf("%s DESIGN_PRIMER=%d\n",T_S(),PSPD_design_primer_flag);
		} else if(strcmp(tok, "BLAST_OGP")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			BLAST_OGP = atoi(tok);
			printf("%s BLAST_OGP=%d\n", T_S(),  BLAST_OGP);
		} else if(strcmp(tok, "BLAST_OIP")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			BLAST_OIP = atoi(tok);
			printf("%s BLAST_OIP=%d\n", T_S(), BLAST_OIP);
		} else if(strcmp(tok, "PSPD_EXPECT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			strcpy(PSPD_EXPECT, tok);
			printf("%s PSPD_EXPECT=%s\n", T_S(), PSPD_EXPECT);
		} else if(strcmp(tok, "PSPD_MINFILTERLENGTH")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PSPD_MINFILTERLENGTH = atoi(tok);
			printf("%s PSPD_MINFILTERLENGTH=%d\n",T_S(),PSPD_MINFILTERLENGTH);
		} else if(strcmp(tok, "PSPD_RUN_SIM_FLAG")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PSPD_RUN_SIM_FLAG = atoi(tok);
			printf("%s PSPD_RUN_SIM_FLAG=%d\n", T_S(), PSPD_RUN_SIM_FLAG);
		} else if(strcmp(tok, "PSPD_MAX_SIMILARITY")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PSPD_MAX_SIMILARITY = (float)atof(tok);
			printf("%s PSPD_MAX_SIMILARITY=%f\n",T_S(),PSPD_MAX_SIMILARITY);
		} else if(strcmp(tok, "PSPD_MINSEGLENGTH")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PSPD_MINSEGLENGTH = atoi(tok);
			printf("%s PSPD_MINSEGLENGTH=%d\n", T_S(), PSPD_MINSEGLENGTH);
		} else if(strcmp(tok, "PRIMER_PRODUCT_SIZE_RANGE")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			strcpy(PPSR, tok);
			printf("%s PRIMER_PRODUCT_SIZE_RANGE=%s\n", T_S(), PPSR);
		} else if(strcmp(tok, "TARGET")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			strcpy(TARGET, tok);
			printf("%s TARGET=%s\n", T_S(), TARGET);
		} else if(strcmp(tok, "PRIMER_SEQUENCE_QUALITY")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			strcpy(PSQ, tok);
			printf("%s PRIMER_SEQUENCE_QUALITY=%s\n", T_S(), PSQ);
		} else if(strcmp(tok, "PRIMER_START_CODON_POSITION")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PSCP = atoi(tok);
			printf("%s PRIMER_START_CODON_POSITION=%d\n", T_S(), PSCP);
		} else if(strcmp(tok, "PRIMER_PICK_ANYWAY")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PPA = atoi(tok);
			printf("%s PRIMER_PICK_ANYWAY=%d\n", T_S(), PPA);
		} else if(strcmp(tok, "PRIMER_EXPLAIN_FLAG")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PEF = atoi(tok);
			printf("%s PRIMER_EXPLAIN_FLAG=%d\n", T_S(), PEF);
		} else if(strcmp(tok, "PRIMER_FILE_FLAG")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PFF = atoi(tok);
			printf("%s PRIMER_FILE_FLAG=%d\n", T_S(), PFF);
		} else if(strcmp(tok, "INCLUDED_REGION")==0) {
			tok = strtok(NULL, " \t\n");
 			strcat(tok, "\0");
			strcpy(PIR, tok);
			printf("%s INCLUDED_REGION=%s\n", T_S(), PIR);
		} else if(strcmp(tok, "EXCLUDED_REGION")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			strcpy(PER, tok);
			printf("%s EXCLUDED_REGION=%s\n", T_S(), PER);	
		} else if(strcmp(tok, "PRIMER_NUM_RETURN")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PNR = atoi(tok);
			printf("%s PRIMER_NUM_RETURN=%d\n", T_S(), PNR);
		} else if(strcmp(tok, "PRIMER_MAX_POLY_X")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PMPX = atoi(tok);
			printf("%s PRIMER_MAX_POLY_X=%d\n", T_S(), PMPX);
		} else if(strcmp(tok, "PRIMER_MIN_TM")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PMINT = (float)atof(tok);
			printf("%s PRIMER_MIN_TM=%f\n", T_S(), PMINT);
		} else if(strcmp(tok, "PRIMER_OPT_TM")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			POT = (float)atof(tok);
			printf("%s PRIMER_OPT_TM=%f\n", T_S(),  POT);
		} else if(strcmp(tok, "PRIMER_MAX_TM")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PMAXT = (float)atof(tok);
			printf("%s PRIMER_MAX_TM=%f\n", T_S(), PMAXT);
		} else if(strcmp(tok, "PRIMER_MIN_SIZE")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PMINS = atoi(tok);
			printf("%s PRIMER_MIN_SIZE=%d\n", T_S(), PMINS);
		} else if(strcmp(tok, "PRIMER_OPT_SIZE")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			POS = atoi(tok);
			printf("%s PRIMER_OPT_SIZE=%d\n", T_S(), POS);
		} else if(strcmp(tok, "PRIMER_MAX_SIZE")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PMAXS = atoi(tok);
			printf("%s PRIMER_MAX_SIZE=%d\n", T_S(), PMAXS);
		} else if(strcmp(tok, "PRIMER_MAX_GC")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PMAXGC = atoi(tok);
			printf("%s PRIMER_MAX_GC=%d\n", T_S(), PMAXGC);
		} else if(strcmp(tok, "PRIMER_MIN_GC")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PMINGC = atoi(tok);
			printf("%s PRIMER_MIN_GC=%d\n", T_S(), PMINGC);
		} else if(strcmp(tok, "PRIMER_SELF_END")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PSE = atoi(tok);
			printf("%s PRIMER_SELF_END=%d\n", T_S(), PSE);
		} else if(strcmp(tok, "PRIMER_SELF_ANY")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PSA = atoi(tok);
			printf("%s PRIMER_SELF_ANY=%d\n", T_S(), PSA);
		} else if(strcmp(tok, "PRIMER_MAX_END_STABILITY")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PMES = (float)atof(tok);
			printf("%s PRIMER_MAX_END_STABILITY=%f\n", T_S(), PMES);
      } else if(strcmp(tok, "PRIMER_PRODUCT_MAX_TM")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PPMAXT = (float)atof(tok);
			printf("%s PRIMER_PRODUCT_MAX_TM=%f\n", T_S(), PPMAXT);
		} else if(strcmp(tok, "PRIMER_PRODUCT_MIN_TM")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
 			PPMINT = (float)atof(tok);
			printf("%s PRIMER_PRODUCT_MIN_TM=%f\n", T_S(), PPMINT);
		} else if(strcmp(tok, "PRIMER_GC_CLAMP")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PGC = atoi(tok);
			printf("%s PRIMER_GC_CLAMP=%d\n", T_S(), PGC);
		} else if(strcmp(tok, "PRIMER_DEFAULT_SIZE")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PMES = atoi(tok);
			printf("%s PRIMER_DEFAULT_SIZE=%d\n", T_S(), PDS);
		} else if(strcmp(tok, "PRIMER_MAX_DIFF_TM")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PMDT = (float)atof(tok);
			printf("%s PRIMER_MAX_DIFF_TM=%f\n", T_S(), PMDT);
		} else if(strcmp(tok, "PRIMER_TM_SANTALUCIA")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PTS = atoi(tok);
			printf("%s PRIMER_TM_SANTALUCIA=%d\n", T_S(), PTS);
		} else if(strcmp(tok, "PRIMER_SALT_CONC")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PSC = (float)atof(tok);
			printf("%s PRIMER_SALT_CONC=%f\n", T_S(), PSC);
		} else if(strcmp(tok, "PRIMER_DIVALENT_CONC")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PDIVC = (float)atof(tok);
			printf("%s PRIMER_DIVALENT_CONC=%f\n", T_S(), PDIVC);
		} else if(strcmp(tok, "PRIMER_DNTP_CONC")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PDNTPC = (float)atof(tok);
			printf("%s PRIMER_DNTP_CONC=%f\n", T_S(), PDNTPC);
		} else if(strcmp(tok, "PRIMER_SALT_CORRECTIONS")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PSALTC = atoi(tok);
			printf("%s PRIMER_SALT_CORRECTIONS=%d\n", T_S(), PSALTC);
		} else if(strcmp(tok, "PRIMER_LOWERCASE_MASKING")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PLM = atoi(tok);
			printf("%s PRIMER_LOWERCASE_MASKING=%d\n", T_S(), PLM);
		} else if(strcmp(tok, "PRIMER_DNA_CONC")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PDNAC = (float)atof(tok);
			printf("%s PRIMER_DNA_CONC=%f\n", T_S(), PDNAC);
		} else if(strcmp(tok, "PRIMER_NUM_NS_ACCEPTED")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PNNSA = atoi(tok);
			printf("%s PRIMER_NUM_NS_ACCEPTED=%d\n", T_S(), PNNSA);
		} else if(strcmp(tok, "PRIMER_LIBERAL_BASE")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PLB = atoi(tok);
			printf("%s PRIMER_LIBERAL_BASE=%d\n", T_S(), PLB);
		} else if(strcmp(tok, "PRIMER_FIRST_BASE_INDEX")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PFBI = atoi(tok);
			printf("%s PRIMER_FIRST_BASE_INDEX=%d\n", T_S(), PFBI);
		} else if(strcmp(tok, "PRIMER_MIN_QUALITY")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PMQ = atoi(tok);
			printf("%s PRIMER_MIN_QUALITY=%d\n", T_S(), PMQ);
		} else if(strcmp(tok, "PRIMER_MIN_END_QUALITY")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PMEQ = atoi(tok);
			printf("%s PRIMER_MIN_END_QUALITY=%d\n", T_S(), PMEQ);
		} else if(strcmp(tok, "PRIMER_QUALITY_RANGE_MIN")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PQRMIN = atoi(tok);
			printf("%s PRIMER_QUALITY_RANGE_MIN=%d\n", T_S(), PQRMIN);
		} else if(strcmp(tok, "PRIMER_QUALITY_RANGE_MAX")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PQRMAX = atoi(tok);
			printf("%s PRIMER_QUALITY_RANGE_MAX=%d\n", T_S(), PQRMAX);
		} else if(strcmp(tok, "PRIMER_INSIDE_PENALTY")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PIP = (float)atof(tok);
			printf("%s PRIMER_INSIDE_PENALTY=%f\n", T_S(), PIP);
		} else if(strcmp(tok, "PRIMER_OUTSIDE_PENALTY")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			POP = (float)atof(tok);
			printf("%s PRIMER_OUTSIDE_PENALTY=%f\n", T_S(), POP);
		} else if(strcmp(tok, "PRIMER_PRODUCT_OPT_TM")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PPOTM = (float)atof(tok);
			printf("%s PRIMER_PRODUCT_OPT_TM=%f\n", T_S(), PPOTM);
		} else if(strcmp(tok, "PRIMER_PRODUCT_OPT_SIZE")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PPOS = atoi(tok);
			printf("%s PRIMER_PRODUCT_OPT_SIZE=%d\n", T_S(), PPOS);
		} else if(strcmp(tok, "PRIMER_WT_TM_GT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PWTG = (float)atof(tok);
			printf("%s PRIMER_WT_TM_GT=%f\n", T_S(), PWTG);
		} else if(strcmp(tok, "PRIMER_WT_TM_LT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PWTL = (float)atof(tok);
			printf("%s PRIMER_WT_TM_LT=%f\n", T_S(), PWTL);
		} else if(strcmp(tok, "PRIMER_WT_SIZE_LT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PWSL = (float)atof(tok);
			printf("%s PRIMER_WT_SIZE_LT=%f\n", T_S(), PWSL);
		} else if(strcmp(tok, "PRIMER_WT_SIZE_GT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PWSG = (float)atof(tok);
			printf("%s PRIMER_WT_TM_LT=%f\n", T_S(), PWSG);
		} else if(strcmp(tok, "PRIMER_WT_GC_PERCENT_LT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PWGPL = (float)atof(tok);
			printf("%s PRIMER_WT_GC_PERCENT_LT=%f\n", T_S(), PWGPL);
		} else if(strcmp(tok, "PRIMER_WT_GC_PERCENT_GT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PWGPG = (float)atof(tok);
			printf("%s PRIMER_WT_GC_PERCENT_GT=%f\n", T_S(), PWGPG);
		} else if(strcmp(tok, "PRIMER_WT_COMPL_ANY")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PWCA = (float)atof(tok);
			printf("%s PRIMER_WT_COMPL_ANY=%f\n", T_S(), PWCA);
		} else if(strcmp(tok, "PRIMER_WT_COMPL_END")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PWCE = (float)atof(tok);
			printf("%s PRIMER_WT_COMPL_END=%f\n", T_S(), PWCE);
		} else if(strcmp(tok, "PRIMER_WT_NUM_NS")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PWNNS = (float)atof(tok);
			printf("%s PRIMER_WT_NUM_NS=%f\n", T_S(), PWNNS);
		} else if(strcmp(tok, "PRIMER_WT_REP_SIM")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PWRS = (float)atof(tok);
			printf("%s PRIMER_WT_REP_SIM=%f\n", T_S(), PWRS);
		} else if(strcmp(tok, "PRIMER_WT_SEQ_QUAL")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PWSQ = (float)atof(tok);
			printf("%s PRIMER_WT_SEQ_QUAL==%f\n", T_S(), PWSQ);
		} else if(strcmp(tok, "PRIMER_WT_END_QUAL")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PWEQ = (float)atof(tok);
			printf("%s PRIMER_WT_END_QUAL==%f\n", T_S(), PWEQ);
		} else if(strcmp(tok, "PRIMER_WT_POS_PENALTY")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PWPP = (float)atof(tok);
			printf("%s PRIMER_WT_POS_PENALTY==%f\n", T_S(), PWPP);
		} else if(strcmp(tok, "PRIMER_WT_END_STABILITY")==0) {
 			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PWES = (float)atof(tok);
			printf("%s PRIMER_WT_END_STABILITY==%f\n", T_S(), PWES);
		} else if(strcmp(tok, "PRIMER_PAIR_WT_PR_PENALTY")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PPWPP = (float)atof(tok);
			printf("%s PRIMER_PAIR_WT_PR_PENALTY==%f\n", T_S(), PPWPP);
		} else if(strcmp(tok, "PRIMER_PAIR_WT_DIFF_TM")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PPWDTM = (float)atof(tok);
			printf("%s PRIMER_PAIR_WT_DIFF_TM==%f\n", T_S(), PPWDTM);
		} else if(strcmp(tok, "PRIMER_PAIR_WT_COMPL_ANY")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PPWCA = (float)atof(tok);
			printf("%s PRIMER_PAIR_WT_COMPL_ANY==%f\n", T_S(), PPWCA);
		} else if(strcmp(tok, "PRIMER_PAIR_WT_COMPL_END")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PPWCE = (float)atof(tok);
			printf("%s PRIMER_PAIR_WT_COMPL_END==%f\n", T_S(), PPWCE);
		} else if(strcmp(tok, "PRIMER_PAIR_WT_PRODUCT_TM_LT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PPWPTL = (float)atof(tok);
			printf("%s PRIMER_PAIR_WT_PRODUCT_TM_LT==%f\n", T_S(), PPWPTL);
		} else if(strcmp(tok, "PRIMER_PAIR_WT_PRODUCT_TM_GT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PPWPTG = (float)atof(tok);
			printf("%s PRIMER_PAIR_WT_PRODUCT_TM_GT==%f\n", T_S(), PPWPTG);
		} else if(strcmp(tok, "PRIMER_PAIR_WT_PRODUCT_SIZE_LT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PPWPSL = (float)atof(tok);
			printf("%s PRIMER_PAIR_WT_PRODUCT_SIZE_LT==%f\n", T_S(), PPWPSL);
		} else if(strcmp(tok, "PRIMER_PAIR_WT_PRODUCT_SIZE_GT")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PPWPSG = (float)atof(tok);
			printf("%s PRIMER_PAIR_WT_PRODUCT_SIZE_GT==%f\n", T_S(), PPWPSG);
		} else if(strcmp(tok, "PRIMER_PAIR_WT_REP_SIM")==0) {
			tok = strtok(NULL, " \t\n");
			strcat(tok, "\0");
			PPWRS = (float)atof(tok);
			printf("%s PRIMER_PAIR_WT_REP_SIM==%f\n", T_S(), PPWRS);
		} else {
			printf("%s error: read line: %s\n", T_S(), line);
			printf("%s token: %s not recognized.\tskipped!\n", T_S(), tok);
			//p_exit("configuration file wrong format\n");
		}
	}
	fclose(fp);
	return PASS;
}

/*check consistency across config file*/
short validate_config_file(void) { 	
	int c=0, b=0, currCount=0, cutSiteCounter=0, *cutSitePos;
	int target_length=0;
	short cutSiteFlag=0;
	char line[SLINE], *tok;

	if((ALGORITHM_TYPE==3)&&(PROBE_DESIGN==0)) {
		if(DATABASE_TYPE !=1)
			p_exit("PSPD algorithm must use single type database!\n");
	}
	if((argf[1]==1)&&(DATABASE_TYPE !=1))
		p_exit("database type should be single\n");

	/* [-l] option with chromosome database*/
	if((argf[1]==2)&&(DATABASE_TYPE==1)) 
		p_exit("[-l] option not matching with database type in config file\n");

	if((argf[1]==2)&&(DATABASE_TYPE==2)) 
		p_exit("current version doesn't support this option!\n");

	/*primer pair count to display must be less than primer pair designed*/
	if(PNR < PRIMER_DISPLAY_COUNT)
		p_exit("total primer pair designed is less than display count!\n");

	/*oligo taken from 3' end shouldn't be less than min. primer size*/
	if(MINBIND > PMINS)
		p_exit("oligo taken from 3' end is larger than min primer size!\n");
   
   /*find all cut-site information based on configuration file*/
	if(CUT_SITE_COUNT > 0) {	
		printf("%s retrieve cut-site information...\n", T_S());
		if(CUT_SITE_REGION==0)
			p_exit("cut-site region is not given!\n");
		target_length = CUT_SITE_REGION;    
		cutSite_count = CUT_SITE_COUNT;
		/*get all cutsite info*/
		cutSite = (CUTSITE *)p_malloc(cutSite_count*sizeof(CUTSITE));
		strcpy(line, CUT_SITE);
		tok = strtok(line, ",:; \t\n");
		if(tok==NULL)
			p_exit("cut-site is null!\n");  
		strcat(tok, "\0");
		if(strlen(tok) >4) 
			p_exit("cut-site is incorrect!\n");
		strcpy(cutSite[0].nseq, tok);
		for(c=1; c < CUT_SITE_COUNT; c++) {
			tok = strtok(NULL, ",:; \t\n");
			if(tok==NULL)
				p_exit("cut-site is null!\n");  
			strcat(tok, "\0");
			if(strlen(tok) >4) 
				p_exit("cut-site is incorrect!\n");
			strcpy(cutSite[c].nseq, tok);
		}   
		if(dlevel >3) {
			printf("%s cut-sites(%d):", T_S(), cutSite_count);  
			for(c=0; c < cutSite_count; c++)
				printf(" %s,", cutSite[c].nseq);
			printf("\n");
		}
	}
	return PASS;
}


/*read input query sequences*/
short read_query(char *fn) {
	FILE *fp;
	char line[LLINE], *tok, tnseq[QNSIZE];
	int count=0, flag=0, q=0;
	short res=0;	/*function return value (pass, fail)*/
	
	/*read query input file and understand its format*/
	read_query_toget_format(fn);
	
	printf("%s read queries from input file %s...\n",T_S(),fn);
	fp = p_fopen_error_exit(fn, "r");
	
	/*testing file format*/
	fgets(line, LLINE, fp);
	if(strncmp(line, ">", 1)!=0)
		p_exit("query file shouldn't have header line!\n"); 	
	rewind(fp);

	/* count the query sequence identifier lines */
	count=0;
	while(fgets(line, LLINE, fp)) { 
		if(strncmp(line, ">", 1)==0) 
			count++;
	}
    rewind(fp);
    q_count = count;
	if(q_count ==0)
		p_exit("No query sequence found as input!\n");
    printf("%s total queries in %s : %d\n", T_S(), fn, count);
    
	/* allocate an array of queries */
	query=(Query *)p_malloc(q_count*sizeof(Query));
	for(q=0; q < q_count; q++) 
		query[q].nseq[0] = 0;
		
	/* If the sequence is present in the query, read it in */
	if(QUERY_SEQUENCE==1)	{
		q=-1;
		while(fgets(line, LLINE, fp)) {
			if(strlen(line) <= 1)
				continue;
			if(strncmp(line, ">", 1)==0) {
				if(q != -1) {
					strcat(query[q].nseq, "\0");
					query[q].len = (int)strlen(query[q].nseq);
				}
				++q;
				tok = strtok(line, "\n");
				strcat(tok, "\0");
				strcpy(query[q].id, tok);
			} else {
				sscanf(line, "%s", tnseq);
				strcat(query[q].nseq, tnseq);
			}
		}
		strcat(query[q].nseq, "\0");
		query[q].len = (int)strlen(query[q].nseq);
		if((q+1)!=q_count)
			p_exit("reading query input counting problem\n");
		for(q=0; q < q_count; q++) {
			printf("%s query [%03d] %s\t",T_S(),q+1,query[q].id);
			printf("length %d\n", query[q].len);
			if(dlevel > 2)
				printf("%s\n", query[q].nseq);
		}
	/* Otherwise, look up the query sequence */
	} else {
		/*check points*/
		if(DATABASE_TYPE==1)
			if((QUERY_FORMAT!=1)&&(QUERY_FORMAT!=2))
				p_exit("query sequence retrieval problem!\n");
		else if(DATABASE_TYPE==2)
			p_exit("feature is under development!\n");
		else if(DATABASE_TYPE==3) {
			if((QUERY_FORMAT!=3)&&(QUERY_FORMAT!=4))
				p_exit("query sequence retrieval problem!\n");
			else if(strncmp(query[q].id, ">chr",4)!=0)
				p_exit("input query name format error!\n");
		} 
		q=0;		/*reading again from start*/
		while(fgets(line, LLINE, fp)) {
			if(strlen(line) < 2)
				continue;			/*junk line*/
			tok = strtok(line, "\n");
			strcat(tok, "\0");
			strcpy(query[q].id, tok);
			printf("%s query [%03d] %s\n",T_S(),q+1,query[q].id);
			q++;
		}
		if(q != q_count)
			p_exit("rewind reading problem!\n");
	}
	fclose(fp);
	printf("%s read query input file %s...done!\n",T_S(),fn);

	/*tokenize query information*/
	res = tokenize_query_info();            /*saved in *qinfo*/
	if(res==FAIL)
		p_exit("problem inside tokenize_query_info function\n");

	return PASS;
}

/*read fspd subset file*/
short read_fspd_sb(char *fn) {
	FILE *fp;
	char line[ELLINE], *tok;
	short count=0, flag=0, q=0, n=0, res=0;
	int frag_start, frag_end, tq_count;
	Query *tquery;
	Qinfo *tqinfo;

	/*read query input file and understand its format*/
	read_query_toget_format(fn);
	
	printf("%s read fspd input file %s...\n",T_S(),fn);
	fp = p_fopen_error_exit(fn, "r");
	
	/* count the query sequence identifier lines */
	count=0;
	while(fgets(line, ELLINE, fp)) 
		if(line[0]=='>') 
			count++;
	q_count = count;
   rewind(fp);
   printf("%s total query in %s : %d\n",T_S(),fn, q_count);
	query=(Query *)p_malloc(q_count*sizeof(Query));
	/* If the query sequences are included in the file, read them in */
	if(QUERY_SEQUENCE==1) {
		q=-1;
		while(!feof(fp)) {
			fgets(line, ELLINE, fp);
			if(strncmp(line, ">", 1)==0) {
				if(q != -1) {
					strcat(query[q].nseq, "\0");
					query[q].len = (int)strlen(query[q].nseq);
				}
				++q;
				tok = strtok(line, "\n");
				strcat(tok, "\0");
				strcpy(query[q].id, tok);
			} else {
				sscanf(line, "%s", tok);
				strcat(query[q].nseq, tok);
			}
		}
		strcat(query[q].nseq, "\0");
		query[q].len = (int)strlen(query[q].nseq);
		if((q+1)!=q_count)
			p_exit("reading query input counting problem\n");
		for(q=0; q < q_count; q++) {
			printf("%s query [%03d] %s\t",T_S(),q+1,query[q].id);
			printf("length %d\n%s\n",query[q].len, query[q].nseq);
		}
	/* Otherwise just read in the queries */
	} else {
		for(q=0; q < q_count; q++) {	
			fgets(line, ELLINE, fp);
			tok = strtok(line, "\n");
			strcat(tok, "\0");
			strcpy(query[q].id, tok);
			printf("%s query [%03d] %s\n",T_S(), q+1, query[q].id);
		}
	}
	fclose(fp);
	printf("%s read fspd input file %s...done!\n",T_S(),fn);
	
	/*tokenize query information*/
	res = tokenize_query_info();            /*saved in *qinfo*/
	if(res==FAIL)
		p_exit("problem inside tokenize_query_info function\n");

	/*copying query data structure into buffer*/
	tq_count = q_count;
	tquery=(Query *)p_malloc(tq_count*sizeof(Query));
	tqinfo=(Qinfo *)p_malloc(tq_count*sizeof(Qinfo));
	for(q=0; q < tq_count; q++) {
		strcpy(tquery[q].id, query[q].id);
		strcpy(tquery[q].nseq, query[q].nseq);
		tquery[q].len = query[q].len;
		strcpy(tqinfo[q].id, qinfo[q].id);
		strcpy(tqinfo[q].id1, qinfo[q].id1);
		strcpy(tqinfo[q].nseq, qinfo[q].nseq);
		strcpy(tqinfo[q].chr, qinfo[q].chr);
		tqinfo[q].start = qinfo[q].start;
		tqinfo[q].end = qinfo[q].end;
		strcpy(tqinfo[q].desc, qinfo[q].desc);
	}
	free(query);
	free(qinfo);
	
	/*fragmenting query sequences*/
	printf("%s fragmenting query sequences with ", T_S());
	printf("length= %d\toverlap=%d...\n", fspd_length, fspd_overlap);
	if((fspd_period==INVALID)&&(fspd_overlap==INVALID))
		p_exit("fspd period or overlap must be defined!\n");	
	if(fspd_length==INVALID)
		p_exit("fspd length must be defined!\n");
	if(fspd_period==INVALID)
		fspd_period = fspd_length - fspd_overlap;
	count=0;
	for(q=0; q < tq_count; q++) {
		frag_start = 0;
		if(fspd_length >= tquery[q].len) 
			frag_end = tquery[q].len -1; 
		else 
			frag_end = fspd_length -1; 
		while(frag_end < tquery[q].len) {
            frag_start += fspd_period;
            frag_end += fspd_period;
			count++;
		}
	}
	q_count = count;
	fprintf(stdout,"%s total query fragments %d\n",T_S(),q_count);
	query=(Query *)p_malloc(q_count*sizeof(Query));
	qinfo=(Qinfo *)p_malloc(q_count*sizeof(Qinfo));
	count=0;
	for(q=0; q < tq_count; q++) {
		frag_start = 0;
		if(fspd_length >= tquery[q].len) 
			frag_end = tquery[q].len -1; 
		else 
			frag_end = fspd_length -1; 
		while(frag_end < tquery[q].len) {
			sprintf(query[count].id,"%s (%d-%d)",tquery[q].id,frag_start,frag_end);
			strcpy(qinfo[count].id, query[count].id);
			strcpy(qinfo[count].id1, tquery[q].id);
			if((QUERY_FORMAT==3)||(QUERY_FORMAT==4)) {
				strcpy(qinfo[count].chr, tqinfo[q].chr);	
				qinfo[count].start = tqinfo[q].start;
				qinfo[count].end = tqinfo[q].end;
			}
			if((QUERY_FORMAT==2)||(QUERY_FORMAT==4)) {
				strcpy(qinfo[count].desc, tqinfo[q].desc);	
			}
           	for(n=frag_start; n <= frag_end; n++) {
				query[count].nseq[n-frag_start]=tquery[q].nseq[n];
				qinfo[count].nseq[n-frag_start]=tquery[q].nseq[n];
			}
			query[count].nseq[n-frag_start] = '\0';
		   	query[count].len = (int)strlen(query[count].nseq);
			frag_start += fspd_period;
            frag_end += fspd_period;
			count++;
		}
	}
	if(count != q_count)
		p_exit("counting problem in fragmenting\n");
	for(q=0; q < q_count; q++) {
		printf("%s query [%03d] %s\t",T_S(),q+1,query[q].id);
		printf("length %d\n%s\n",query[q].len, query[q].nseq);
	}
	fprintf(stdout, "%s fragmenting query sequences...done!\n", T_S());
	free(tquery);
	return PASS;
}


/* split the queries into the various query parameters */
short tokenize_query_info(void) {
    char *tok;
    char line[MLINE];
    int q=0, c1=0, c2=0;
	short res = 0;	/*function return value*/

    printf("%s tokenize query name information...\n",T_S());
    qinfo=(Qinfo *)p_malloc(q_count*sizeof(Qinfo));
	if((QUERY_FORMAT==3)||(QUERY_FORMAT==4)) {
		for(q=0; q < q_count; q++) {
    	    strcpy(line, query[q].id);	/*original qid line*/
    	    tok=strtok(line,":# \t");
    	    strcat(tok,"\0");
			if(strncmp(tok, ">", 1)==0)
				tok++;
    	    strcpy(qinfo[q].chr, tok);
    	    tok=strtok(NULL,"-# \t");
    	    strcat(tok,"\0");
    	    qinfo[q].start = atoi(tok);
    	    tok=strtok(NULL,"# \t\n");
    	    strcat(tok,"\0");
    	    qinfo[q].end = atoi(tok);
			if(QUERY_FORMAT==4){
    	   		tok=strtok(NULL,"\n");
    	    	strcat(tok,"\0");
    	    	strcpy(qinfo[q].desc, tok);
			}
    	    strcpy(qinfo[q].id, query[q].id);
			qinfo[q].id[strlen(qinfo[q].id)-1]='\0';
			if(QUERY_SEQUENCE==1)
    	    	strcpy(qinfo[q].nseq, query[q].nseq);
			if(dlevel >=3) {
				printf("%s%6d) qinfo.chr: %s\n",T_S(),q,qinfo[q].chr);
				printf("%s%6d) qinfo.start: %d\n",T_S(),q,qinfo[q].start);
				printf("%s%6d) qinfo.end: %d\n",T_S(),q,qinfo[q].end);
				printf("%s%6d) qinfo.desc: %s\n",T_S(),q,qinfo[q].desc);
			}
    	}
	} else if((QUERY_FORMAT==1)||(QUERY_FORMAT==2))	{
		for(q=0; q < q_count; q++) {
    	    strcpy(line, query[q].id);	/*original qid line*/
    	    tok=strtok(line,"# \t\n");
    	    strcat(tok,"\0");
    	    strcpy(qinfo[q].id, tok);
			if(qinfo[q].id[0]=='>')
				c1=1;
			else 
				c1=0;
			for(c2=0; c1 < (short)strlen(qinfo[q].id); c1++, c2++) {
				if(!isalnum(qinfo[q].id[c1]))
					qinfo[q].id1[c2] = '_';
				else
					qinfo[q].id1[c2] = qinfo[q].id[c1];
			}
			qinfo[q].id1[c2] = '\0';
			if(QUERY_FORMAT==2) {
    	    	tok = strtok(NULL, "\n");
				strcat(tok,"\0");
    	    	strcpy(qinfo[q].desc, tok);
				if(dlevel >3)
					printf("\tqinfo[%d].desc=%s\n", q, qinfo[q].desc);	
			}
		}
	}
	/*format cross-checking*/
	if(QUERY_SEQUENCE==0) {
	if(((QUERY_FORMAT==1)||(QUERY_FORMAT==2))&&(DATABASE_TYPE!=1))
		p_exit("to retrieve sequence, chose correct database, query format.\n");
	if(((QUERY_FORMAT==3)||(QUERY_FORMAT==4))&&(DATABASE_TYPE!=3))
		p_exit("to retrieve sequence, chose correct database, query format.\n");
	}
			
	/*retrieving query sequences in case needed*/
	if(QUERY_SEQUENCE==0) {
		for(q=0; q < q_count; q++) {
			if(DATABASE_TYPE==1) {
				res=get_nseq_from_db(query[q].id,query[q].nseq,&query[q].len);
				if(res==FAIL)
					p_exit("problem in function get_nseq_from_db!\n");
			} else if(DATABASE_TYPE==2){
				p_exit("database type not supported yet!\n");	
			} else if(DATABASE_TYPE==3) {
				res = get_nseq_from_chromosome(q);
				if(res==FAIL)
					p_exit("problem in function get_nseq_from_chromosome!\n");
			} else
				p_exit("database type error!\n");	
		}
	}
    printf("%s tokenize query name information...done!\n",T_S());
    return PASS;
}

/* read the file containing the list of databases */
short read_dbl(char *fdl) {
	int count=0, i=0, res=0;
	FILE *fp;
	char line[SLINE];

    printf("%s read database list input file %s...\n",T_S(),fdl);
    count=0;
    fp=p_fopen_error_exit(fdl, "r");
    
    /* count the files */
    while(fgets(line, SLINE, fp))
        count++;
    rewind(fp);
    db_count = count;
    printf("%s total database count = %d\n",T_S(),count);
	if(DATABASE_TYPE==3) {
		if(db_count != BLASTRECORDS) 
			p_exit("input db list is not same as blast records!\n");
	}
	
	/* store the filenames in fn_db */
	fn_db=(Dfile *)p_malloc(count*sizeof(Dfile));
    for(i=0; i < db_count; i++)
        fscanf(fp, "%s",fn_db[i].fn);
    fclose(fp);
   /* REMOVE? */
	res = read_dbl_test(fdl);
	if(res==FAIL)
		p_exit("problem in read_dbl_test function\n");
    printf("%s read dbl input file %s...done!\n",T_S(),fdl);
	return PASS;
}


short read_dbl_test(char *fn) {
	short i=0;

    for(i=0; i<db_count; i++) {
        printf("%s Database %03d )  \t%s.\n",T_S(),i+1,fn_db[i].fn);
    }
	return PASS;
}

/************************* END ***********************************/


