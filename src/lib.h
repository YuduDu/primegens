

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

#define NA "?"
#define CONFIGFORMAT 200	/*configuration file format string*/
#define QNAME 1000			/*query name string size*/
#define ONAME 1000			/*oligo name string size*/
#define QNSIZE 99999		/*query sequence string size*/
#define SEQUENCEID 200		/*sequence name in general like single db*/
#define QDESCRIPTION 1000	/*query function and descriptions*/
#define FILENAME 300		/*filename string size*/
#define PATHNAME 300		/*pathname string size*/
#define CHROMOSOME 10		/*string like chr5*/
#define BISULCHR 50			/*bisulphite chromosome id*/
#define SLINE 100			/*small length file lines*/
#define MLINE 9999			/*medium length file lines*/
#define LLINE 99999			/*large length file lines*/
#define ELLINE 100000		/*extra large length file lines*/
#define BLINE 1000			/*blast output file line*/
#define FASTALINE 60		/*fasta format line size for output*/
#define CHR_FASTA_LINE 50	/*fasta line size in chromosome database*/
#define OLIGO 40			/*primer oligo lengh*/
#define CUTSITESIZE 6		/*cut sitesjstring size*/ 
#define ALIGNMENT 40		/*example plus/plus in blast*/
#define COMMAND 500			/*command line string*/
#define BCOMMAND 10000		/*blast commands string size*/
#define FCOMMAND 10000		/*formatdb commands string size*/
#define LARGESCORE 9999999	/*default value of blast alignment score*/
#define LINEWIDTH 80		/*sequence alingment display line width*/
#define BLASTSEQID 60		/*used for matching id with blast hit seqid*/
#define EMPTY -1			/*original primegens code*/
#define PP 1				/*plus/plus type*/
#define PM 0				/*plus/minus type*/
#define PASS 1
#define FAIL 0
#define INVALID -1			/*primer_cg_max not available*/
#define HLEN 100

/*PSPD options*/
#define OneIndelPenalty 2
#define OpenGapPenalty  4 
#define PSPD_NumBlastMatch 1000
#define MAXHOMOLOGS 100				/*max homolog to display*/
#define LargeScore 9999999
#define TPP "Plus"			/* right of Plus/Plus strand type*/
#define TPM "Minus"			/* right of Plus/Minus strand type*/


#define TPSPD "product-specific primer design "			
#define TSSPD "sequence-specific primer design"		
#define TGSPD "gene-specific primer design"	
#define TFSPD "fragment-specific primer design"

#define BHITTAG "Sequences producing significant alignments:"
#define BHITCOMP 43

#define QSEQFILENAME "qseq.tmp"	/*dump query sequence into file*/ 
#define QSEQBLASTFILENAME "qseq_blast.tmp"	/*query blast output*/

/*primer design parameters*/
#define PRIMER3TAG 60			/*primer3 output tag string size*/
#define LOGFILENAME "primegens.log"
#define FNCFG "config.txt"
#define PATHVAR "PRIMEGENS_PATH"
#define TMPDIRPATHVAR "TMPDIR_PATH"
#define VDATABASE "virtual_database"
#define PRIMER3IN "primer3.input"
#define PRIMER3OUT "primer3.output"
#define UOLIGOTMP "unique_oligos.tmp"
#define BLASTIN "blast_in.fasta"
#define BLASTOUT "blast_out.fasta"
#define FORMATDBLOG "format_database.log"
#define STATUSFILE "Run_Status.txt"

/*output results naming convention*/
#define EXCELSHEET "_primer.xls"
#define PRIMERLIST "_primers_list.txt"
#define QUERYPRIMER "_query_primer_list.txt"
#define QUERYFAILED "_query_failed.txt"
#define QUERYFASTA "_query_fasta.txt"
#define QUERYNOHIT "_nohit.txt"
#define QUERYUNISEG "_uniseg.txt"
#define QUERYSIMFILE "_sim.txt"
#define QUERYPROBEA "_probe.txt"
#define QUERYPROBEB "_sequence_hit_info.xls"
#define PROBEPREFIX "probe_"
#define HOMOA "_homologs.txt"

/*last character in workdir is mandatory*/
#ifdef _WIN32
	#define ALIGNMENTDIR "\\alignments\\"		/*for pspd*/
	#define BLASTALLEXEC "blastall.exe"		/*for pspd*/
#else
	#define ALIGNMENTDIR "/alignments/"		/*for pspd*/
	#define BLASTALLEXEC "blastall"			/*for pspd*/
#endif


#define TSTAMP 40


/*configuration file variables*/
short BLAST_OGP;
short BLAST_OIP;
short WORDSIZE;				/*used for running blast*/
short MINBIND;				/*min primer oligo bind from right*/ 
short BLASTRECORDS;
short BLASTLOOP;
short CHRDBINLOOP;
char FORMATDBEXEC[CONFIGFORMAT];
char BLASTEXEC[CONFIGFORMAT];
char PRIMERDESIGNEXEC[CONFIGFORMAT];
char ALIGNMENTEXEC[CONFIGFORMAT];
int AMPSB;
short QUERY_FORMAT;
short DATABASE_TYPE;			/*single or multiple sequences*/
short QUERY_SEQUENCE;			/*if sequence is given or not*/
short ALGORITHM_TYPE;
short PROBE_DESIGN;
short BISULFITE_MODE;
float PROBE_MAX_GC;
float PROBE_MIN_GC;
short PROBE_BLAST_HIT;
short PROBE_LENGTH;
short PROBE_PERIOD;
short PROBE_REGION;
short PROBE_GAP;
short PROBE_UNIQUE;
short dlevel;
short HBRDN_DISPLAY_COUNT;
short PRIMER_DISPLAY_COUNT;
short PRIMER_CG_MAX;
short CHRCOUNT;

char PPSR[PRIMER3TAG];	                /*PRIMER_PRODUCT_SIZE_RANGE=250-450*/
char TARGET[PRIMER3TAG];       	        /*target region to specify*/
char PIR[PRIMER3TAG];                   /*INCLUDED_REGION*/
char PER[PRIMER3TAG];                   /*EXCLUDED_REGION*/
char PSQ[PRIMER3TAG];                   /*PRIMER_SEQUENCE_QUALITY*/
int  PSCP;                                 /*PRIMER_START_CODON_POSITION=-1000000*/
short PPA;                                  /*PRIMER_PICK_ANYWAY=0*/

char PML[PRIMER3TAG];                                   /*PRIMER_MISPRIMING_LIBRARY*/
short PLACC;                                /*PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=1*/
short PMM;                                  /*PRIMER_MAX_MISPRIMING=12*/
float PMTM;                                 /*PRIMER_MAX_TEMPLATE_MISPRIMING=-1.00*/
short PPMM;                                 /*PRIMER_PAIR_MAX_MISPRIMING=24*/
float PPMTM;                                /*PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=-1.00*/
float PPMAXT;                               /*PRIMER_PRODUCT_MAX_TM=1000000.0*/
float PPMINT;                               /*PRIMER_PRODUCT_MIN_TM=-1000000.0*/
short PPIO;                                 /*PRIMER_PICK_INTERNAL_OLIGO=0*/
short PGC;                                  /*PRIMER_GC_CLAMP=0*/
short PDS;                                  /*PRIMER_DEFAULT_SIZE=20*/
float PMDT;                                 /*PRIMER_MAX_DIFF_TM=100.0C*/
short PTS;                                  /*PRIMER_TM_SANTALUCIA=0*/
float PSC;                                  /*PRIMER_SALT_CONC=50.0 mM*/
float PDIVC;                                /*PRIMER_DIVALENT_CONC=0.0 mM*/
float PDNTPC;                               /*PRIMER_DNTP_CONC=0.0 mM*/
short PSALTC;                               /*PRIMER_SALT_CORRECTIONS=0*/
short PLM;                                  /*PRIMER_LOWERCASE_MASKING=0*/
float PDNAC;                                /*PRIMER_DNA_CONC=50.0nM*/
short PNNSA;                                /*PRIMER_NUM_NS_ACCEPTED=0*/
short PLB;                                  /*PRIMER_LIBERAL_BASE=0*/
short PFBI;                                 /*PRIMER_FIRST_BASE_INDEX=0*/
short PMQ;                                  /*PRIMER_MIN_QUALITY=0*/
short PMEQ;                                 /*PRIMER_MIN_END_QUALITY=0*/
short PQRMIN;                               /*PRIMER_QUALITY_RANGE_MIN=0*/
short PQRMAX;                               /*PRIMER_QUALITY_RANGE_MAX=100*/
float PIP;                                  /*PRIMER_INSIDE_PENALTY=-1.0*/
float POP;                                  /*PRIMER_OUTSIDE_PENALTY=0.0*/
float PPOTM;                                /*PRIMER_PRODUCT_OPT_TM=0.0*/
short PPOS;                                 /*PRIMER_PRODUCT_OPT_SIZE=0*/
char PT[PRIMER3TAG];                                    /*PRIMER_TASK=pick_pcr_primers*/

float PWTG;                                  /*PRIMER_WT_TM_GT=1.0*/
float PWTL;                                 /*PRIMER_WT_TM_LT=1.0*/
float PWSL;                                 /*PRIMER_WT_SIZE_LT=1.0*/
float PWSG;                                 /*PRIMER_WT_SIZE_GT=1.0*/
float PWGPL;                                 /*PRIMER_WT_GC_PERCENT_LT=1.0*/
float PWGPG;                                 /*PRIMER_WT_GC_PERCENT_GT=1.0*/
float PWCA;                                 /*PRIMER_WT_COMPL_ANY=0.0*/
float PWCE;                                 /*PRIMER_WT_COMPL_END=0.0*/
float PWNNS;                                 /*PRIMER_WT_NUM_NS=0.0*/
float PWRS;                                 /*PRIMER_WT_REP_SIM=0.0*/
float PWSQ;                                 /*PRIMER_WT_SEQ_QUAL=0.0*/
float PWEQ;                                 /*PRIMER_WT_END_QUAL=0.0*/
float PWPP;                                 /*PRIMER_WT_POS_PENALTY=0.0*/
float PWES;                                 /*PRIMER_WT_END_STABILITY=0.0*/
float PWTM;                                 /*PRIMER_WT_TEMPLATE_MISPRIMING=0.0*/
float PPWPP;                                 /*PRIMER_PAIR_WT_PR_PENALTY=1.0*/
float PPWIOP;                                 /*PRIMER_PAIR_WT_IO_PENALTY=0.0*/
float PPWDTM;                                 /*PRIMER_PAIR_WT_DIFF_TM=0.0*/
float PPWCA;                                 /*PRIMER_PAIR_WT_COMPL_ANY=0.0*/
float PPWCE;                                 /*PRIMER_PAIR_WT_COMPL_END=0.0*/
float PPWPTL;                                 /*PRIMER_PAIR_WT_PRODUCT_TM_LT=0.0*/
float PPWPTG;                                 /*PRIMER_PAIR_WT_PRODUCT_TM_GT=0.0*/
float PPWPSL;                                 /*PRIMER_PAIR_WT_PRODUCT_SIZE_LT=0.0*/
float PPWPSG;                                 /*PRIMER_PAIR_WT_PRODUCT_SIZE_GT=0.0*/
float PPWRS;                                 /*PRIMER_PAIR_WT_REP_SIM=0.0*/
float PPWTM;                                 /*PRIMER_PAIR_WT_TEMPLATE_MISPRIMING=0.0*/


short PEF;				/*PRIMER_EXPLAIN_FLAG=1*/
short PFF;				/*PRIMER_FILE_FLAG=0*/
short PNR;				/*PRIMER_NUM_RETURN=20*/
short PMPX;				/*PRIMER_MAX_POLY_X=0*/
short PIOMPX;				/*PRIMER_INTERNAL_OLIGO_MAX_POLY_X=0*/
float PMINT;			        /*PRIMER_MIN_TM=50.0*/
float POT;				/*PRIMER_OPT_TM=58.0*/
float PMAXT;			        /*PRIMER_MAX_TM=60.0*/
short PMINS;				/*PRIMER_MIN_SIZE=20*/
short POS; 				/*PRIMER_OPT_SIZE=24*/
short PMAXS;				/*PRIMER_MAX_SIZE=28*/
float PMAXGC;				/*PRIMER_MAX_GC=80*/
float POPTGC;				/*PRIMER_OPT_GC_PERCENT=50*/
float PMINGC;				/*PRIMER_MIN_GC=20*/
short PSE;				/*PRIMER_SELF_END=100*/
short PSA;				/*PRIMER_SELF_ANY=100*/
float PMES;				/*PRIMER_MAX_END_STABILITY=100.0*/

/*fragment specific primer design parameters*/
short fspd_length;		/*only for fspd option*/
short fspd_period;		/*only for fspd option*/
short fspd_overlap;		/*only for fspd option*/
short CUT_SITE_COUNT;	/*count of cut-sites*/
int CUT_SITE_REGION;	/*region to calculate cut-site density*/
char CUT_SITE[CONFIGFORMAT];	/*cut-site list*/

/*valid only for original pspd*/
short PSPD_design_primer_flag;	
short PSPD_BLAST_OGP;
short PSPD_BLAST_OIP;
short PSPD_MINFILTERLENGTH;
short PSPD_RUN_SIM_FLAG;
float PSPD_MAX_SIMILARITY;
short PSPD_MINSEGLENGTH;
char PSPD_EXPECT[CONFIGFORMAT];

int pid;				/*process id*/
short log_flag;			/*generate log*/

char fn_fmtdblog[FILENAME];
char fn_log[FILENAME];
char fn_status[FILENAME];
char fn_vdb[FILENAME];
char fn_p3in[FILENAME];
char fn_p3out[FILENAME];
char fn_uoligo[FILENAME];
char fn_bin[FILENAME];
char fn_bout[FILENAME];
char fn_query[FILENAME];
char fn_blast[FILENAME];
char fn_dbl[FILENAME];

char fn_out1[FILENAME];			/*output type-1*/
char fn_out2[FILENAME];			/*output type-2*/
char fn_out3[FILENAME];			/*output type-3*/
char fn_out4[FILENAME];			/*output type-4*/
char fn_out5[FILENAME];			/*output type-4*/
char fn_nohit[FILENAME];
char fn_uniseg[FILENAME];
char fn_sim[FILENAME];
char fn_probe[FILENAME];
char fn_probe_hit[FILENAME];

char t_Stamp[TSTAMP];			/*dummy time stamp*/
char env_path[PATHNAME];		/*environment variable*/
char tmpDir_path[PATHNAME];		/*environment variable*/
char tmp_dir[PATHNAME];			/*temporary directory*/
char alignDir[PATHNAME];		/*alignment directory*/
char dbpath[PATHNAME];			/*database path location*/
char out_path[PATHNAME];		/*output result location*/
short argf[4];

char primer3[COMMAND];
char blast[COMMAND];
char blastall[COMMAND];		/*original pspd*/
char megablast[COMMAND];
char formatdb[COMMAND];
