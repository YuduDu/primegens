/*data structures used in whole program*/

typedef struct {
	char fn[FILENAME];
} Dfile;

Dfile *fn_db;
int db_count;		/*used in case of multiple/chromosome*/
Dfile fn_db0;		/*used in case if single database type*/


typedef struct {
	char nseq[CUTSITESIZE];
} CUTSITE;
CUTSITE *cutSite;
short cutSite_count;

typedef struct {
	char id[SEQUENCEID];
} SeqID;		/*used in case of single database type*/

typedef struct {
	char cid[SEQUENCEID];
	int ppc;				/*plus-plus hit count*/
	int pmc;				/*plus-minus hit count*/
	int *pps;			/*plus-plus start position*/
	int *pme;			/*plus-minus end position*/
	int *ppe;			/*plus-plus end position*/
	int *pms;			/*plus-minus start position*/
	char hlseq[HLEN];	/*hybridizing left sequence*/
	char hrseq[HLEN];	/*hybridizing right sequence*/
	int hllen;
	int hrlen;
	float lseq_oligo_tm;
	float rseq_oligo_tm;
	int potential;
} CHit;

CHit *chit;			/* buffer for each primer pair*/
int chit_count;


typedef struct {
	char id[BISULCHR];
} Chr;


typedef struct {
	char chr[SEQUENCEID];	/*bisulphite chromosome id*/	
	int bh;						/*blast hits count*/
	int *bht;					/*blast hit type*/
	int *bhs;					/*blast hit start*/
	int *bhe;					/*blast hit end*/
} BRecord;


typedef struct {
	char id[ONAME];
	char nseq[OLIGO];
	int len;
	BRecord *brecord;
	int br_count;			/*used in case of single database*/
} UOligo;

UOligo *uoligo;
int uoligo_count;


typedef struct {			/*super set of amplicon*/
	char id[SEQUENCEID];
	int hs;					/*hybridization start*/
	int he;					/*hybridization end*/
	int hsize;				/*hybridization size*/
} Hbrdn;


typedef struct {			/*save primer design result for Q*/
	char lid[ONAME];		/*left primer id*/
	char rid[ONAME];		/*right primer id*/
	char lseq[OLIGO];		/*left primer sequence*/
	char rseq[OLIGO];		/*right primer sequence*/
	int lst;					/*left primer start position*/
	int rst;					/*right primer start position*/
	int llen;				/*left primer length*/
	int rlen;				/*right primer length*/
	double ltm;				/*left primer melting temperature*/
	double rtm;				/*right primer melting temperature*/
	double lgc;				/*left primer %GC content*/
	double rgc;				/*right primer %GC content*/
	int psize;				/*product size due to primer pair*/
	CHit *chit;				/*blast hit info for each primer pair*/
	int chit_count;		/*blast record in case of single database*/
	int hbrdn_count;		/*total no. hybridization on genome*/
	Hbrdn *hbrdn;			/*decsription of each hybridization*/
} PPair;

PPair *ppair;
int ppair_count;


/*
 * Minimum information about query sequence used during the 
 * process of primer design.
 */
typedef struct {
	char id[MLINE];			/*query sequence name/identifier*/
	char nseq[QNSIZE];		/*nucleotide sequence of query*/
	int len;						/*query sequence length*/		
	short res;					/*primer design fail or pass*/
	PPair bestpp;				/*best primer pair designed*/
}Query;

Query *query;
int q_count;


/*
 * information about each of the query gene, independent 
 * of the query input format
 */
typedef struct {
	char id[QNAME];				/*query sequence name/identifier*/
	char id1[QNAME];				/*query name with removed wild chars*/
	char nseq[QNSIZE];			/*nucleotide sequence of the query*/
	char chr[BISULCHR];			/*chromosome location if applicable*/
	int start;						/*start position on genome if applicable*/
	int end;							/*end positiion on genome if applicable*/
	char desc[QDESCRIPTION];	/*functional description if applicable*/
}Qinfo;

Qinfo *qinfo;


typedef struct {
        char query_id[QNAME];
        char subject_id[QDESCRIPTION];
        int perc_identity;
        int align_length;
        int mismatches;
        int gap_open;
        int q_start;
        int q_end;
        int s_start;
        int s_end;
        int e_value;
        int bit_score;
} HOMOLOG;
HOMOLOG *b_hit;


typedef struct {
        int start;
        int end;
        int hit_no;
} Hit;
Hit *hit, hit_struct,overlap_arr[10000][10000],sorted[100],temp;

