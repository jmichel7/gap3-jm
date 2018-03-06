/* gpmult.c 28/11/94.
 * This simply calculates minimizes and outputs all multipliers for a
 * short-lex automatic group, using the general multiplier automaton.
 * It assumes that kbprog with -wd option, gpwa, gpgenmult
 * and preferably gpcheckmult  have already been run of G.
 *
 * SYNOPSIS:
 * gpmult [-ip d/s] [-op d/s] [-silent] [-v] [-l/-h] groupname
 *
 * Input is from groupname.gm
 * Output is to groupname.mi for each generator number i.
 *
 * OPTIONS:
 * -ip d/s	input in dense or sparse format - sparse is default
 * -op d/s	output in dense or sparse format - sparse is default
 * -v		verbose
 * -silent	no diagnostics
 * -l/-h	large/huge hash-tables (for big examples)
 */

#include <stdio.h>
#include "defs.h"
#include "fsa.h"
#include "definitions.h"

FILE *rfile,
     *wfile,
     *fopen();

void  badusage();

/* Functions defined in other files used in this file */
void  fsa_read();
void  fsa_makemult();
void  fsa_minimize();
void  fsa_print();
void  fsa_clear();
void  fsa_table_dptr_init();

main(argc, argv)
        int             argc;
        char           *argv[];
{ int arg, g, nlabels;
  fsa genmult, *mult;
  char groupname[100], inf[100], outf[100], fsaname[100];
  storage_type ip_store = SPARSE;
  boolean op_format_set = FALSE;
  storage_type op_format = SPARSE;

  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  groupname[0] = '\0';
  arg = 1;
  while (argc > arg) {
    if (strcmp(argv[arg],"-ip")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      if (strcmp(argv[arg],"d")==0)
        ip_store = DENSE;
      else if (strcmp(argv[arg],"s")==0)
        ip_store = SPARSE;
      else
        badusage();
    }
    else if (strcmp(argv[arg],"-op")==0) {
      op_format_set = TRUE;
      arg++;
      if (arg >= argc)
        badusage();
      if (strcmp(argv[arg],"d")==0)
        op_format = DENSE;
      else if (strcmp(argv[arg],"s")==0)
        op_format = SPARSE;
      else
        badusage();
    }
    else if (strcmp(argv[arg],"-silent")==0)
      print_level = 0;
    else if (strcmp(argv[arg],"-v")==0)
      print_level = 2;
    else if (strcmp(argv[arg],"-vv")==0)
      print_level = 3;
    else if (strcmp(argv[arg],"-l")==0)
      large = TRUE;
    else if (strcmp(argv[arg],"-h")==0)
      huge = TRUE;
    else {
       if (argv[arg][0] == '-')
         badusage();
       if (strcmp(groupname,"")!=0)
         badusage();
       else
         strcpy(groupname,argv[arg]);
    }
    arg++;
  }
  if (strlen(groupname)==0)
    badusage();
  
  strcpy(inf,groupname);
  strcat(inf,".gm");

  if ((rfile = fopen(inf,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf);
      exit(1);
  }
  fsa_read(rfile,&genmult,ip_store,0,0,TRUE,fsaname);

  tmalloc(mult,fsa,1);
  nlabels = genmult.states->labels->size;
/* nlabels may either be the numebr of generators ngens of the group, or
 * ngens+1 (in which case label number ngens+1 gives the equality
 * multiplier). The latter occurs if gpgenmult was run with th -id option.
 */

  for (g=1;g<=nlabels;g++) {
 /* Form multiplier number g */
    fsa_copy(mult,&genmult);
    fsa_makemult(mult,g);
    fsa_minimize(mult);

    base_prefix(fsaname);
    sprintf(fsaname+strlen(fsaname),".m%d",g);
    sprintf(outf,"%s.m%d",groupname,g);

    if (op_format_set)
      mult->table->printing_format = op_format;
    wfile = fopen(outf,"w");
    fsa_print(wfile,mult,fsaname);
    fclose(wfile);

    if (print_level>0)
      printf(
       "#Multiplier number %d with %d states computed.\n",
        g,mult->states->size);
    fsa_clear(mult);
  }

  fsa_clear(&genmult);
  tfree(mult);
  exit(0);
}
 
void
badusage()
{
    fprintf(stderr,
      "Usage: gpmult [-ip d/s] [-op d/s] [-silent] [-v] [-l/-h] groupname.\n");
    exit(1);
}
