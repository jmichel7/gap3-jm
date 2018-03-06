/* gpcheckmult.c 18/11/94.
 * This program assumes that kbprog, gpwa, gpgenmult have already been
 * run on the group G. 
 * The theory of this checking process for the general multiplier is as follows.
 *
 * For each generator g of G, let M_g be the associated multiplier automaton,
 * and construct the fsa E_g, which accepts word w iff there exists a word x
 * such that M_g accepts (w,x). E_g should equal the word-acceptor W for G.
 * Construct E_g and W' to check this. It should have empty language. If not
 * find words w that it accepts, reduce w*g, to x, say, and update the
 * word-difference machine for G used to construct the multipliers, to make
 * it accept the pair (w*g,x).
 *
 * In practice, the function fsa_checkmult does this check on the general
 * multiplier for all g at the same time.
 * Furthermore, it does not need to use the word-acceptor explicitly, since
 * a pair of words (w,x) can only be accepted by the mutliplier if w and x
 * are accepted by the word-acceptor.
 *
 * SYNOPSIS:
 * gpcheckmult [-ip d/s[dr]] [-v] [-silent] [-l/-h] [-m maxeqns] 
 *		[-mwd maxwdiffs] groupname
 *
 * Input is from groupname.wa, groupname.gpgenmult, and possibly groupname.diff2
 * (the last only if it requires updating).
 * If no updating is required there is no output.
 * Otherwise output is to groupname.diff2.
 *
 * EXIT STATUS:
 * This is 2 if groupname.diff2 requires updating, otherwise 0.
 *
 * OPTIONS:
 * -ip d/s[dr]  input in dense or sparse format - dense is default
 * -v           verbose
 * -silent      no diagnostics
 * -l/-h        large/huge hash-tables (for big examples)
 * -m  maxeqns  Abort after finding maxeqns offending words w (see above)
 *              default is MAXEQNS
 * -mwd maxwdiffs
 *		At most maxwdiffs word-differences possible.
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

fsa wd_fsa; /* This is for the case that we correct the diff2 machine */

void  badusage();
void  (*reduce_word)();

/* Functions defined in other files used in this file */
void  fsa_read();
int   fsa_checkmult();
void  fsa_clear();
void  diff_reduce();
boolean add_wd_fsa();
void  make_full_wd_fsa();
void  calculate_inverses();
void fsa_table_dptr_init();

main(argc, argv)
        int             argc;
        char           *argv[];
{ int arg, i, *inv, old_ndiff, maxeqns, numeqns, ngens, maxwdiffs;
  fsa  genmult;
  char inf1[100], inf2[100], outf[100], fsaname[100];
  reduction_equation *eqnptr;
  storage_type ip_store = DENSE;
  int dr = 0;

  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  maxeqns = MAXEQNS;
  maxwdiffs = MAXWDIFFS;
  inf1[0] = '\0';
  inf2[0] = '\0';
  outf[0] = '\0';
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
    else if (strcmp(argv[arg],"-m")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      maxeqns = atoi(argv[arg]);
    }
    else if (strcmp(argv[arg],"-mwd")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      maxwdiffs = atoi(argv[arg]);
    }
    else {
       if (argv[arg][0] == '-')
         badusage();
       if (strcmp(inf1,"")!=0)
         badusage();
       strcpy(inf1,argv[arg]);
    }
    arg++;
  }
  if (strlen(inf1)==0)
    badusage();
  strcpy(inf2,inf1);
  strcat(inf1,".gm");

  if ((rfile = fopen(inf1,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf1);
      exit(1);
  }
  fsa_read(rfile,&genmult,ip_store,dr,0,TRUE,fsaname);
  fclose(rfile);

  tmalloc(eqnptr,reduction_equation,maxeqns)
  if ((numeqns = fsa_checkmult(&genmult,eqnptr,maxeqns)) > 0) {
 /* A multiplier was not valid, so groupname.diff2 will need updating. */

    fsa_clear(&genmult);
    strcat(inf2,".diff2");
    strcpy(outf,inf2);
    if (print_level>1)
      printf("  #Altering wd-machine to make it accept new equations.\n");
    if ((rfile = fopen(inf2,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",inf2);
          exit(1);
    }
    fsa_read(rfile,&wd_fsa,DENSE,0,maxwdiffs,TRUE,fsaname);
    fclose(rfile);
/* We will use this machine for word-reduction even while correcting it!  */
    reduce_word = diff_reduce;
    fsa_table_dptr_init(&wd_fsa);
/* We need to know the inverses of generators - let's just work them out!  */
    ngens = wd_fsa.alphabet->base->size;
    calculate_inverses(&inv,ngens);
    old_ndiff = wd_fsa.states->size;

/* Now add the new equations
 * The right hand side of the equation to be added will be the reduction of
 * the lhs times the generator which is currently in the lhs.
 */
    for (i=0;i<numeqns;i++) {
      strcat(eqnptr[i].lhs,eqnptr[i].rhs);
      tfree(eqnptr[i].rhs);
      tmalloc(eqnptr[i].rhs,char,1+strlen(eqnptr[i].lhs));
      strcpy(eqnptr[i].rhs,eqnptr[i].lhs);
      diff_reduce(eqnptr[i].rhs);
      add_wd_fsa(&wd_fsa,eqnptr+i,inv,TRUE);
    }
    
    make_full_wd_fsa(&wd_fsa,inv,old_ndiff+1);

    if (print_level>1)
      printf("  #Word-difference machine now has %d states.\n",
               wd_fsa.states->size);

    wfile = fopen(inf2,"w");
    fsa_print(wfile,&wd_fsa,fsaname);
    fclose(wfile);

    tfree(inv);
    fsa_clear(&wd_fsa);
    for (i=0;i<numeqns;i++) {
      tfree(eqnptr[i].lhs);
      tfree(eqnptr[i].rhs);
    }
    tfree(eqnptr);
    exit(2);
  }

  if (print_level>0)
    printf("#Validity test on general multiplier succeeded.\n");

  fsa_clear(&genmult);
  tfree(eqnptr);
  exit(0);
}
 
void
badusage()
{
    fprintf(stderr,
     "Usage: gpcheckmult [-ip d/s[dr]] [-silent] [-v] [-l/-h] [-m maxeqns]\n");
    fprintf(stderr,"		[-mwd maxwdiffs]  groupname\n");
    exit(1);
}
