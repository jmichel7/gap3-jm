/* gpgenmult.c 3/11/94.
 * Calculates the general multiplier automaton of a short-lex automatic group.
 * This program assumes that kbprog with -wd option and gpwa have already
 * been run of G.
 *
 * SYNOPSIS:
 * gpgenmult [-op d/s] [-silent] [-v] [-l/-h] [-c] [-mwd maxwdiffs]
 *           [-m maxeqns] [-ns] [-f] groupname
 *
 * Input is from groupname.wa and gpname.diff2
 *  (and possibly from groupname.diff1 if -c option is called).
 * Output is to groupname.gm
 *  (or possibly to groupname.diff1 if -c option is called).
 *
 * OPTIONS:
 * -op d/s	output in dense or sparse format - sparse is default
 * -v		verbose
 * -silent	no diagnostics
 * -l/-h	large/huge hash-tables (for big examples)
 * -c		(correction)
 *              If equations are discovered which proves the word-acceptor
 *              to be incorrect, then the first word-difference machine
 *              (which should be in file groupname.diff1) is updated by making
 *              it accept these equations.
 *              There is no point in calling this option if gpwa was run
 *              with input from groupname.diff2, since such an equation
 *              cannot occur in that case.
 * -mwd maxwdiffs
 *              At most maxwdiffs word-differences possible.
 * -m  maxeqns
 *		(only relevant if -c is called)
 *              Abort the multiplier checking process after finding maxeqns
 *              offending words w (see above) -  default is MAXEQNS
 * -ns		Don't stop if nontrivial equation found in word-acceptor
 *		language.
 * -f           read the transition table repeatedly from file while mimimizing.
 *              this avoids storing the large table, but is a little slower.
 *
 * EXIT STATUS:
 * If new equations are discovered, then the exit status is 2
 * (whether or not the correction option was called).
 * Otherwise, with normal output to groupname.genmult, status is 0.
 */

#include <stdio.h>
#include "defs.h"
#include "fsa.h"
#include "rws.h"
#include "definitions.h"
#define MAXWDIFFS 2048
#define MAXEQNS 512

FILE *rfile,
     *wfile,
     *fopen();

fsa wd_fsa; /* This is for the case that we correct the diff1 machine */

void  badusage();
void  (*reduce_word)();

/* Functions defined in other files used in this file */
void  fsa_read();
fsa  *fsa_triples();
void  fsa_print();
void  fsa_labeled_minimize();
void  fsa_ip_labeled_minimize();
void  fsa_clear();
void  diff_reduce();
boolean add_wd_fsa();
void  calculate_inverses();
void fsa_table_dptr_init();

main(argc, argv)
        int             argc;
        char           *argv[];
{ int arg, i, *inv, ngens, maxwdiffs, maxeqns;
  fsa wa, diff, *genmultptr;
  char inf1[100], inf2[100], inf3[100], outf[100], fsaname[100],
       tempfilename[100];
  reduction_equation *eqnptr;
  boolean correction = FALSE;
  storage_type op_store = SPARSE;
  boolean idlabel = TRUE; /* always true for the moment */
  boolean eqnstop = TRUE;
  boolean readback = TRUE;

  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  maxwdiffs = MAXWDIFFS;
  maxeqns = MAXEQNS;
  inf1[0] = '\0';
  arg = 1;
  while (argc > arg) {
    if (strcmp(argv[arg],"-op")==0) {
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
    else if (strcmp(argv[arg],"-c")==0)
      correction = TRUE;
    else if (strcmp(argv[arg],"-ns")==0)
      eqnstop = FALSE;
    else if (strcmp(argv[arg],"-f")==0)
      readback = FALSE;
    else if (strcmp(argv[arg],"-mwd")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      maxwdiffs = atoi(argv[arg]);
    }
    else if (strcmp(argv[arg],"-m")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      maxeqns = atoi(argv[arg]);
    }
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
  
  strcpy(tempfilename,inf1);
  strcat(tempfilename,"temp_triples_XXX");

  strcpy(inf2,inf1);
  strcat(inf2,".diff2");

  if (correction) {
    strcpy(inf3,inf1);
    strcat(inf3,".diff1");
  }

  strcpy(outf,inf1);
  strcat(outf,".gm");

  strcat(inf1,".wa");

  if ((rfile = fopen(inf1,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf1);
      exit(1);
  }
  fsa_read(rfile,&wa,DENSE,0,0,TRUE,fsaname);
  fclose(rfile);
  if ((rfile = fopen(inf2,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf2);
    exit(1);
  }
  fsa_read(rfile,&diff,DENSE,0,0,TRUE,fsaname);
  fclose(rfile);

  if (correction)
    tmalloc(eqnptr,reduction_equation,maxeqns)
  else {
    eqnptr = 0;
    maxeqns = 0;
  }

  genmultptr = fsa_triples(&wa,&diff,op_store,TRUE,tempfilename,eqnptr,maxeqns,
						idlabel,eqnstop,readback);

  if (genmultptr==0) { /*This is the case where new equations were found. */
    if (correction)  {
      if (print_level>1)
        printf("  #Altering wd-machine to make it accept new equations.\n");
      if ((rfile = fopen(inf3,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",inf3);
          exit(1);
      }
      fsa_read(rfile,&wd_fsa,DENSE,0,maxwdiffs,TRUE,fsaname);
      fclose(rfile);
/* We will use this machine for word-reduction even while correcting it! */
      reduce_word = diff_reduce;
      fsa_table_dptr_init(&wd_fsa);
/* We need to know the inverses of generators - let's just work them out! */
      ngens = wd_fsa.alphabet->base->size;
      calculate_inverses(&inv,ngens);

      i=0;
      while (eqnptr[i].lhs && i<maxeqns) {
        add_wd_fsa(&wd_fsa,eqnptr[i],inv,FALSE);
        i++;
      }

      if (print_level>1)
        printf("  #Word-difference machine now has %d states.\n",
                 wd_fsa.states->size);

      wfile = fopen(inf3,"w");
      fsa_print(wfile,&wd_fsa,fsaname);
      fclose(wfile);

      tfree(inv);
      fsa_clear(&wd_fsa);
      i=0;
      while (eqnptr[i].lhs && i<maxeqns) {
        tfree(eqnptr[i].lhs); tfree(eqnptr[i].rhs);
        i++;
      }
      tfree(eqnptr);
    }
    exit(2);
  }

  if (print_level>1)
    printf("  #Number of states of triples before minimisation = %d.\n",
        genmultptr->states->size);
  if (readback)
    fsa_labeled_minimize(genmultptr);
  else
    fsa_ip_labeled_minimize(genmultptr);
  if (print_level>1)
    printf("  #Number of states of triples after minimisation = %d.\n",
        genmultptr->states->size);

  base_prefix(fsaname);
  strcat(fsaname,".gm");
  wfile = fopen(outf,"w");
    fsa_print(wfile,genmultptr,fsaname);
  fclose(wfile);

  if (print_level>0)
    printf("#Generalised multiplier with %d states computed.\n",
            genmultptr->states->size);

  fsa_clear(genmultptr);
  tfree(genmultptr);
  if (correction)
    tfree(eqnptr);
  exit(0);
}
 
void
badusage()
{
    fprintf(stderr,
   "Usage: gpgenmult [-op d/s] [-silent] [-v] [-l/-h] [-c] [-mwd maxwdiffs]\n");
    fprintf(stderr,"		[-m maxeqns] [-ns] [-f] groupname.\n");
    exit(1);
}
