/* ppgap.c   27.2.95. 
 *
 * Read a kb-input file and make a preprocessor file for GAP
 *
 * SYNOPSIS:
 * ppgap  groupname
 *
 */

#include <stdio.h>

#include "defs.h"
#include "fsa.h"
#include "rws.h"
#include "definitions.h"

#define MAXEQNS		1024
#define MAXREDUCELEN	4096

int  	 maxeqns = MAXEQNS,
         maxreducelen = MAXREDUCELEN;
boolean  maxreducelenset = FALSE; /* not used in this program */
rewriting_system  rws;

static char gpname[100],  outf[100];
FILE *rfile, *wfile;

/* Functions defined in this file */
void badusage();

/* Functions used in this file defined in other files: */
void read_kbinput_simple();
void printbuffer();
void add_to_buffer();
void rws_clear();

void main(argc, argv)
        int             argc;
        char           *argv[];
{ int  i, l, ct;
  boolean first;

  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  if (argc!=2) badusage();
  strcpy(gpname,argv[1]);
  strcpy(outf,gpname); strcat(outf,".gap");
  
/* First read in the defining relations for the group. */
  if ((rfile = fopen(gpname,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",gpname);
      exit(1);
  }
  read_kbinput_simple(rfile,TRUE);
  fclose(rfile);

  wfile = fopen(outf,"w");

  buffer[0]='\0';
  add_to_buffer(0,"_RWS.gpMon := FreeGroup(");
  first = TRUE;

  if (rws.num_gens==0) add_to_buffer(0,"0");
  else for (i=1;i<=rws.num_gens;i++){
     l = strlen(rws.gen_name[i]);
     if (l<=3 || strcmp(rws.gen_name[i]+l-3,"^-1")) {
       if (!first)
          add_to_buffer(0,",");
       first = FALSE;
       sprintf(buffer+strlen(buffer),"\"%s\"",rws.gen_name[i]);
     }
  }
  add_to_buffer(0,");");
  printbuffer(wfile);

  ct = 0;
  for (i=1;i<=rws.num_gens;i++) {
    l = strlen(rws.gen_name[i]);
    if (l<=3 || strcmp(rws.gen_name[i]+l-3,"^-1")) {
       ct++;
       fprintf(wfile,"%s := _RWS.gpMon.%d;\n",rws.gen_name[i],ct);
    }
  }
  fprintf(wfile,"_ := IdWord;\n");
  fclose(wfile);
  rws_clear(&rws);
  exit(0);
}
 
void
badusage()
{ fprintf(stderr, "Usage: ppgap groupname.\n");
	exit(1);
}
