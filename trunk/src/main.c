
/* header files*/

#include <unistd.h>
#include "lib.h"
#include "defs.h"
#include "input.h"
#include "util.h"
#include "probe.h"
#include "sspd.h"
#include "fspd.h"
#include "pspd.h"

const char *hmsg =
    "\n\n\t\t\tPRIMEGNENS-v2.0 help\n"
    "\t $>primegens.exe [-q] <?> [-d/l] <?> [-p] <?> [-o] <?>\n"
    "\t -q: query sequence list filename\n"
    "\t -d/l: database/list filename\n"
    "\t -p: database location path\n"
    "\t -o: output result location\n"
    "\t -h: help\n\n" ;

const char *pmsg =
	"\n\t\tEnvironment variable setup\n"
	"\tThe value of PRIMEGENS environment variable must the location of\n"
	"\tPRIMEGENS root folder. This is the folder, which contains other\n"
	"\tsubfolder like bin, blast, primer3, include etc.\n" 
	"\tFor more detail, please follow README.txt\n" ;
	
	
/*main program*/
int main(int argc, char *argv[]) {
	char *ptr;
	int res=0;
	FILE *fp;


	ptr = (char *)getenv(PATHVAR);
	if(!ptr) { 
		/* get the parent directory */
		chdir("..");
		char cwd[1024];
		if (getcwd(cwd, sizeof(cwd)) != NULL)
			fprintf(stdout, "Environment Variable for path not defined, using current parent dir: %s\n", cwd);
		else
			p_exit("getcwd() error");
		strcpy(env_path, cwd);
		/* move back into the bin */
		chdir("bin");
	}
	else {
		strcpy(env_path, ptr);
	}
	pid = getpid();
	ptr = (char *)getenv(TMPDIRPATHVAR);
	if(!ptr) { 
		printf("%s tmpDir environment variable not defined. Copying\n",T_S());	
		strcpy(tmpDir_path, env_path);
		printf("%s tmpDir environment variable=%s\n",T_S(), tmpDir_path);
	}
	else {
		strcpy(tmpDir_path, ptr);	
		printf("%s tmpDir environment variable=%s\n",T_S(), tmpDir_path);
	}
	
	//show;
    
	if((argc==2)&&(strcmp(argv[1],"-h")==0)) {
        fprintf(stderr, hmsg);
        exit(0);
    }
    if(argc < 9) {
        fprintf(stderr, hmsg);
        fprintf(stderr, "%s arguments count: %d\n", T_S(), argc);
        p_exit("less no. of arguments\n");
    }
	
	/*initialize required parameters*/
	printf("%s initializing required parameters...\n", T_S());
	res = init_param();
	if(res==FAIL)
		p_exit("problem inside init_param function\n");

	/*getting input parameters*/
	printf("%s get parameters from command line...\n", T_S());
	res = get_param(argc, argv);
	if(res==FAIL)
		p_exit("problem inside get_param function\n");


	/*read configuration file*/
	printf("%s read configuration file...\n", T_S());
	res = read_config_file();
	if(res==FAIL)
		p_exit("problem reading configuration file\n");
	/*validation*/

	
	res = validate_config_file();   /*consistency check*/
	if(res==FAIL)
		p_exit("Inconsistent configuration\n");

	/*setting default parameters*/
	res = set_default_param();
	if(res==FAIL)
		p_exit("problem inside set_default_param function\n");
	
	if(DATABASE_TYPE ==3) {
		res = read_dbl(fn_dbl);         /*saved in *fn_db*/
		if(res==FAIL)
			p_exit("problem inside read_dbl function\n");
	}
	
	/*select primer-probe design algorithm*/
	if(PROBE_DESIGN==1) {
		printf("%s starting probe design module...\n", T_S());	
                res = probe_loop();
		if(res==FAIL)
			p_exit("problem in probe design loop function\n");
	} else if(ALGORITHM_TYPE==1) {
		fprintf(stdout, "%s %s\n", T_S(), TSSPD);
		
		res = sspd_loop();
		if(res==FAIL)
			p_exit("problem in sspd loop function\n");
	} else if(ALGORITHM_TYPE==3) { 
		fprintf(stdout, "%s\n", T_S(), TPSPD);
		res = pspd_loop();
		if(res==FAIL)
			p_exit("problem in pspd loop function\n");
	} else if(ALGORITHM_TYPE==2) { 
		fprintf(stdout, "%s\n", T_S(), TFSPD);
		res = fspd_loop();
		if(res==FAIL)
			p_exit("problem in fspd loop function\n");
	} else { 
		p_exit("no algorithm selected!\n");
	}
	/*final cleaning*/
	p_clean();

    printf("\n\n");
    printf("**************************** ");
    printf("primegens execution successfully terminated");
    printf(" **************************** ");
    printf("\n\n");

	fp = p_fopen_error_exit(fn_status, "w");
	fprintf(fp,"PRIMEGENSv2 completed successfully!\n");
	fclose(fp);

	return 0;	/*special return for successful completion*/
}


