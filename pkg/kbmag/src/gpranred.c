/* gpranred.c 3/11/94.
 * Calculates the general multiplier automaton of a short-lex automatic group.
 * This program assumes that kbprog with -wd option and gpwa have already
 * been run of G.
 *
 * SYNOPSIS:
 * gpranred [-silent] [-v] [-l/-h] groupname
 *
 * Input is from groupname.wa and gpname.diff2
 *  (and possibly from groupname.diff1 if -c option is called).
 * Output is to groupname.ranred
 *  (or possibly to groupname.diff1 if -c option is called).
 *
 * OPTIONS:
 * -v		verbose
 * -silent	no diagnostics
 * -l/-h	large/huge hash-tables (for big examples)
 *
 * EXIT STATUS:
 * If a new equation is discovered, then the exit status is 2
 * Otherwise, with normal output to groupname.ranred, status is 0.
 */

#include <stdio.h>
#include "defs.h"
#include "fsa.h"
#include "rws.h"
#include "definitions.h"
#define MAXWDIFFS 8192
#define NUMTRIALS 256
#define MAXLEN    24

FILE *rfile,
     *wfile,
     *fopen();

fsa wd_fsa; /* This is for the case that we correct the diff1 machine */

void  badusage();
void  (*reduce_word)();

/* Functions defined in other files used in this file */
void  fsa_read();
void  fsa_print();
void  fsa_labeled_minimize();
void  fsa_clear();
void  diff_reduce();
boolean add_wd_fsa();
void  calculate_inverses();
void fsa_table_dptr_init();

main(argc, argv)
        int             argc;
        char           *argv[];
{ int arg, i, r, g, s, *inv, trial, len, ngens, **table, *validgen, nvg;
  int sc, **tablec;
  fsa wa, diff1, *gpranred, wac;
  char inf1[100], inf2[100], inf3[100], fsaname[100], infc[100],
       word1[16384], word2[16384];
  int max_len = MAXLEN, num_trials = NUMTRIALS;
  reduction_equation *eqnptr;
  boolean changed;

  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  inf1[0] = '\0';
  arg = 1;
  while (argc > arg) {
    if (strcmp(argv[arg],"-silent")==0)
      print_level = 0;
    else if (strcmp(argv[arg],"-v")==0)
      print_level = 2;
    else if (strcmp(argv[arg],"-vv")==0)
      print_level = 3;
    else if (strcmp(argv[arg],"-l")==0)
      large = TRUE;
    else if (strcmp(argv[arg],"-h")==0)
      huge = TRUE;
    else if (strcmp(argv[arg],"-nt")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      num_trials = atoi(argv[arg]);
    }
    else if (strcmp(argv[arg],"-ml")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      max_len = atoi(argv[arg]);
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
  
  srand(num_trials+max_len);

  strcpy(inf2,inf1);
  strcat(inf2,".diff2");

  strcpy(inf3,inf1);
  strcat(inf3,".wa");

strcpy(infc,inf1);
strcat(infc,".wak");

  strcat(inf1,".diff1");

  if ((rfile = fopen(inf3,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf3);
      exit(1);
  }
  fsa_read(rfile,&wa,DENSE,0,TRUE,fsaname);
  fclose(rfile);
if ((rfile = fopen(infc,"r")) == 0) {
  fprintf(stderr,"Cannot open file %s.\n",infc);
    exit(1);
}
fsa_read(rfile,&wac,DENSE,0,TRUE,fsaname);
  fclose(rfile);

  if ((rfile = fopen(inf1,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf1);
    exit(1);
  }
  fsa_read(rfile,&diff1,DENSE,MAXWDIFFS,TRUE,fsaname);
  fclose(rfile);

  if ((rfile = fopen(inf2,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf2);
    exit(1);
  }
  fsa_read(rfile,&wd_fsa,DENSE,0,TRUE,fsaname);
  fclose(rfile);
/* We shall be using the diff2 machine for reduction. */

  reduce_word = diff_reduce;
  tmalloc(eqnptr,reduction_equation,1)
  eqnptr->lhs = word1; eqnptr->rhs = word2;
  fsa_table_dptr_init(&diff1);
  fsa_table_dptr_init(&wd_fsa);
  ngens = wd_fsa.alphabet->base->size;
  table = wa.table->table_data_ptr;
tablec = wac.table->table_data_ptr;
  tmalloc(validgen,int,ngens);
  calculate_inverses(&inv,ngens);
  changed = FALSE;

  for (trial=1;trial<=num_trials;trial++) {
    if (print_level>=2)
      printf ("  #Trial number %d.\n",trial);
    word1[0] = '\0';
    s = wa.initial[1];
sc = wac.initial[1];
    for (len=1;len<=max_len;len++) {
      nvg = 0;
      for (i=1;i<=ngens;i++) if (table[i][s]) {
        validgen[nvg++] = i;
if (tablec[i][sc]==0)
printf("WOTS ALL THIS THEN?");
}
      if (nvg==0)
        break;
      r = rand(); g = r%nvg;
      word1[len-1] = validgen[g];
      word1[len] = '\0';
      s = table[validgen[g]][s];
sc = tablec[validgen[g]][sc];
    }
    strcpy(word2,word1);
    (*reduce_word)(word2);
    if (strcmp(word1,word2)) {
      if (print_level>1)
        printf("  #Altering wd-machine to make it accept new equation.\n");
      changed = TRUE;
      add_wd_fsa(&diff1,*eqnptr,inv,FALSE);
      if (print_level>=2)
        printf("  #Word-difference machine now has %d states.\n",
               wd_fsa.states->size);
    }
  }

  if (changed && print_level>0)
    printf("#Word-difference machine now has %d states.\n",
           wd_fsa.states->size);
  fsa_clear(&wa);
  fsa_clear(&diff1);
  fsa_clear(&wd_fsa);

  tfree(validgen);
  tfree(eqnptr);

  if (changed) exit(2);
  exit(0);
}
 
void
badusage()
{
    fprintf(stderr,
   "Usage: gpranred [-silent] [-v] [-l/-h]  groupname.\n");
    exit(1);
}
