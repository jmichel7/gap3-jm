/* gpmult2.c 25/11/94.
 * Calculates a multiplier automaton of a short-lex automatic group,
 * for words in the generators of length 2, using the length-2
 * general multiplier.
 * This program assumes that kbprog with -wd option, gpwa, gpgenmult
 * gpgenmult2 and preferably gpcheckmult  have already been run of G.
 *
 * SYNOPSIS:
 * gpmult2 [-ip d/s] [-op d/s] [-silent] [-v] [-l/-h] g1 g2 groupname
 *
 * Input is from groupname.gm2 (and possibly groupname.gm2_ut).
 * and individual 2-generator multipliers to groupname.mult_g1_g2,
 * where  g1  and  g2  are the generator numbers. 
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
#include "rws.h"
#include "definitions.h"

FILE *rfile,
     *wfile,
     *fopen();

void  badusage();

/* Functions defined in other files used in this file */
void  fsa_read();
void  fsa_makemult2();
void  fsa_minimize();
void  fsa_print();
void  fsa_clear();

main(argc, argv)
        int             argc;
        char           *argv[];
{ int arg, i, g1, g2;
  fsa genmult2;
  char inf[100], outf[100], fsaname[100];
  storage_type ip_store = SPARSE;
  boolean op_format_set = FALSE;
  storage_type op_format = SPARSE;

  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  inf[0] = '\0';
  g1 = g2 = 0;
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
       if (g1==0 && isdigit(argv[arg][0]))
         g1 = atoi(argv[arg]);
       else if (g2==0 && isdigit(argv[arg][0]))
         g2 = atoi(argv[arg]);
       else if (strcmp(inf,"")!=0)
         badusage();
       else
         strcpy(inf,argv[arg]);
    }
    arg++;
  }
  if (g1==0 || g2==0 || strlen(inf)==0)
    badusage();
  
  sprintf(outf,"%s.m%d_%d",inf,g1,g2);
  strcat(inf,".gm2");

  if ((rfile = fopen(inf,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf);
      exit(1);
  }
  fsa_read(rfile,&genmult2,ip_store,0,0,TRUE,fsaname);

  if (print_level>1)
    printf("  #Number of states of genmult2 = %d.\n",genmult2.states->size);

  base_prefix(fsaname);
  sprintf(fsaname+strlen(fsaname),".m%d_%d",g1,g2);
  fsa_makemult2(&genmult2,g1,g2);
  fsa_minimize(&genmult2);
  if (print_level>0)
    printf("#Length-2 multiplier with %d states computed.\n",
        genmult2.states->size);

  if (op_format_set)
    genmult2.table->printing_format = op_format;

  wfile = fopen(outf,"w");
  fsa_print(wfile,&genmult2,fsaname);
  fclose(wfile);

  fsa_clear(&genmult2);
  exit(0);
}
 
void
badusage()
{
    fprintf(stderr,
"Usage: gpmult2 [-ip d/s] [-op d/s] [-silent] [-v] [-l/-h] g1 g2 groupname.\n");
    exit(1);
}
