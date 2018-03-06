/* gpminkb.c 8/11/94.
 * Calculates the fsa which accepts minimal reducible words in a short-lex
 * automatic group, and then the 2-variable fsa which accepts the set of
 * minimal KB-equations. These have the minimal reducible words as LHS
 * and the corresponding irreducible words as RHS.
 * Finally it uses the minimal KB-machine to calculate the correct
 * diff1 machine.
 * The input is the word-acceptor and word-difference machine diff2
 * (diff1 would sometimes work, but diff2 is safer).
 *
 * SYNOPSIS:  gpminkb [-op1 d/s] [-op2 d/s] [-silent] [-v] [-l/-h] groupname
 *
 * Input is from groupname.wa and groupname.diff2
 * Output is to groupname.minred, groupname.minkb and groupname.diff1c.
 *
 * OPTIONS:
 * -op1 d/s	output of minred in dense or sparse format - dense is default
 * -op2 d/s	output of minkb and diff1c in dense or sparse format -
 *                                                          sparse is default
 * -v		verbose
 * -silent	no diagnostics
 * -l/-h	large/huge hash-tables (for big examples)
 *		(former is default). This should be the same as was used
 *		to construct the word-acceptor.
 * 
 */

#include <stdio.h>
#include "defs.h"
#include "fsa.h"
#include "definitions.h"

FILE *rfile,
     *wfile,
     *fopen();

fsa   wd_fsa; /* for word reduction */
void  badusage();

/* Functions defined in other files used in this file */
void  fsa_read();
fsa  *fsa_minred();
fsa  *fsa_minkb();
fsa  *fsa_diff1c();
void  fsa_print();
void  fsa_clear();

main(argc, argv)
        int             argc;
        char           *argv[];
{ int arg;
  fsa wa, *minred, *minkb, *diff1c;
  char gpname[100], inf[100], outf[100], fsaname[100], tempfilename[100];
  storage_type op_store1 = DENSE;
  storage_type op_store2 = SPARSE;

  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  gpname[0] = '\0';
  arg = 1;
  while (argc > arg) {
    if (strcmp(argv[arg],"-op1")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      if (strcmp(argv[arg],"d")==0)
        op_store1=DENSE;
      else if (strcmp(argv[arg],"s")==0)
        op_store1=SPARSE;
      else
        badusage();
    }
    else if (strcmp(argv[arg],"-op2")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      if (strcmp(argv[arg],"d")==0)
        op_store2=DENSE;
      else if (strcmp(argv[arg],"s")==0)
        op_store2=SPARSE;
      else
        badusage();
    }
    else if (strcmp(argv[arg],"-silent")==0 || strcmp(argv[arg],"-s")==0)
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
       if (strcmp(gpname,""))
         badusage();
       strcpy(gpname,argv[arg]);
    }
    arg++;
  }
  if (strlen(gpname)==0)
    badusage();

  strcpy(inf,gpname);
  strcat(inf,".wa");
  strcpy(outf,gpname);
  strcat(outf,".minred");

  strcpy(tempfilename,gpname);
  strcat(tempfilename,"temp_min_XXX");

  if ((rfile = fopen(inf,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf);
    exit(1);
  }
  fsa_read(rfile,&wa,DENSE,0,0,TRUE,fsaname);
  fclose(rfile);

  minred = fsa_minred(&wa,op_store1,FALSE,tempfilename);

  if (print_level>1)
    printf("  #Number of states of minred before minimisation = %d.\n",
        minred->states->size);
  fsa_minimize(minred);
  if (print_level>1)
    printf("  #Number of states of minred after minimisation = %d.\n",
        minred->states->size);

  base_prefix(fsaname);
  strcat(fsaname,".minred");
  wfile = fopen(outf,"w");
  fsa_print(wfile,minred,fsaname);
  fclose(wfile);
  if (print_level>0)
    printf("#Minimal reducible word acceptor with %d states computed.\n",
            minred->states->size);


/* Now do the second part - construct the minimal KB-rules machine. */
  strcpy(inf,gpname);
  strcat(inf,".diff2");
  strcpy(outf,gpname);
  strcat(outf,".minkb");

  if ((rfile = fopen(inf,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf);
    exit(1);
  }
  fsa_read(rfile,&wd_fsa,DENSE,0,0,TRUE,fsaname);
  fclose(rfile);

  minkb = fsa_minkb(minred,&wa,&wd_fsa,op_store2,FALSE,tempfilename);

  fsa_clear(&wa);
  fsa_clear(minred);
  tfree(minred);

  if (print_level>1)
    printf("  #Number of states of minkb before minimisation = %d.\n",
        minkb->states->size);
  fsa_minimize(minkb);
  if (print_level>1)
    printf("  #Number of states of minkb after minimisation = %d.\n",
        minkb->states->size);

  base_prefix(fsaname);
  strcat(fsaname,".minkb");
  wfile = fopen(outf,"w");
  fsa_print(wfile,minkb,fsaname);
  fclose(wfile);

  if (print_level>0)
    printf("#Minimal Knuth-Bendix equation fsa with %d states computed.\n",
            minkb->states->size);

  fsa_clear(minkb);

/* Now the final part - construct the correct diff1 machine. */
  strcpy(inf,gpname);
  strcat(inf,".minkb");
  strcpy(outf,gpname);
  strcat(outf,".diff1c");

/* We re-read the minkb machine to get it into dense format. */
  if ((rfile = fopen(inf,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf);
    exit(1);
  }
  fsa_read(rfile,minkb,DENSE,0,0,TRUE,fsaname);
  fclose(rfile);

  diff1c = fsa_diff1c(minkb,op_store2);

  if (print_level>1)
    printf("  #Number of word-differences = %d.\n",diff1c->states->size);

  base_prefix(fsaname);
  strcat(fsaname,".diff1c");
  wfile = fopen(outf,"w");
  fsa_print(wfile,diff1c,fsaname);
  fclose(wfile);

  if (print_level>0)
    printf("#Minimal diff1 fsa with %d states computed.\n",
            diff1c->states->size);

  fsa_clear(minkb);
  tfree(minkb);
  fsa_clear(diff1c);
  tfree(diff1c);
  exit(0);
}

void
badusage()
{
    fprintf(stderr, "Usage:\n");
    fprintf(stderr,
      "gpminkb [-op1 d/s] [-op2 d/s] [-silent] [-v] [-l/-h] groupname\n");
    exit(1);
}
