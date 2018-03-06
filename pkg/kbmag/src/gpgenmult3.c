/* gpgenmult3.c 28/11/94.
 * Calculates the general multiplier automaton of a short-lex automatic group,
 * for words in the generators of length 3.
 * This program assumes that kbprog with -wd option, gpwa, gpgenmult
 * and gpgenmult2 have already been run of G.
 *
 * SYNOPSIS:
 * gpgenmult3 [-ip d/s[dr]] [-op d/s] [-silent] [-v] [-l/-h] [-f] groupname
 *
 * Input is from groupname.gm
 * Output is to groupname.gm3 and the unformatted transition
 * table goes to groupname.gm3_ut
 *
 * OPTIONS:
 * -ip d/s[dr]	input in dense or sparse format - sparse is default
 * -op d/s	output in dense or sparse format - sparse is default
 * -v		verbose
 * -silent	no diagnostics
 * -l/-h	large/huge hash-tables (for big examples)
 * -f           don't read back the transition table and mimimize.
 */

#include <stdio.h>
#include "defs.h"
#include "fsa.h"
#include "definitions.h"

FILE *rfile1, *rfile2,
     *wfile,
     *fopen();

void  badusage();

/* Functions defined in other files used in this file */
void  fsa_read();
fsa  *fsa_genmult3();
void  fsa_labeled_minimize();
void  fsa_print();
void  fsa_clear();
void  fsa_table_dptr_init();

main(argc, argv)
        int             argc;
        char           *argv[];
{ int arg, i, g1, g2;
  fsa genmult, genmult2, *genmult3ptr;
  char inf1[100], inf2[100], outf[100], fsaname[100], tablefilename[100];
  storage_type ip_store = DENSE;
  int dr = 0;
  storage_type op_store = SPARSE;
  boolean readback = TRUE;

  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  inf1[0] = '\0';
  arg = 1;
  while (argc > arg) {
    if (strcmp(argv[arg],"-ip")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      if (strcmp(argv[arg],"d")==0)
        ip_store = DENSE;
      else if (argv[arg][0] == 's') {
        ip_store = SPARSE;
        if (strlen(argv[arg]) > 1)
          dr = atoi(argv[arg]+1);
      }
      else
        badusage();
    }
    else if (strcmp(argv[arg],"-op")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      if (strcmp(argv[arg],"d")==0)
        op_store = DENSE;
      else if (strcmp(argv[arg],"s")==0)
        op_store = SPARSE;
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
    else if (strcmp(argv[arg],"-f")==0)
      readback = FALSE;
    else {
       if (argv[arg][0] == '-')
         badusage();
       if (strcmp(inf1,"")!=0)
         badusage();
       else
         strcpy(inf1,argv[arg]);
    }
    arg++;
  }
  if (strlen(inf1)==0)
    badusage();
  
  strcpy(tablefilename,inf1);
  strcat(tablefilename,".gm3_ut");

  strcpy(inf2,inf1);
  strcat(inf2,".gm2");
  strcpy(outf,inf1);
  strcat(outf,".gm3");
  strcat(inf1,".gm");

  if ((rfile1 = fopen(inf1,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf1);
      exit(1);
  }
  if ((rfile2 = fopen(inf2,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf2);
      exit(1);
  }
  fsa_read(rfile1,&genmult,ip_store,dr,0,TRUE,fsaname);
  fsa_read(rfile2,&genmult2,ip_store,dr,0,TRUE,fsaname);
  fclose(rfile1); fclose(rfile2);

  genmult3ptr =
      fsa_genmult3(&genmult,&genmult2,op_store,TRUE,tablefilename,readback);

  if (print_level>1)
    printf("  #Number of states of genmult3 = %d.\n",genmult3ptr->states->size);

  if (readback) {
    fsa_labeled_minimize(genmult3ptr);
    if (print_level>1)
      printf("  #Number of states of genmult3 after minimization = %d.\n",
               genmult3ptr->states->size);
  }
  base_prefix(fsaname);
  strcat(fsaname,".gm3");
  wfile = fopen(outf,"w");
  fsa_print(wfile,genmult3ptr,fsaname);
  fclose(wfile);

  if (print_level>0)
    printf("#Generalised length-3 multiplier with %d states computed.\n",
           genmult3ptr->states->size);

  fsa_clear(genmult3ptr);
  tfree(genmult3ptr);
  exit(0);
}
 
void
badusage()
{
    fprintf(stderr,"Usage: \n");
    fprintf(stderr,
"gpgenmult3 [-ip d/s[dr]] [-op d/s] [-silent] [-v] [-l/-h] [-f] groupname.\n");
    exit(1);
}
