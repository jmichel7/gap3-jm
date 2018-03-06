/* gpmakefsa.c 11/1/95.
 * This program performs the combined functions of gpwa, gpgenmult and
 * gpcheckmult, to make the finite state automata of a short-lex automatic
 * structure of a group G. It runs them repeatedly, until they pass the
 * correctness test.
 * This program assumes that kbprog with -wd option has already
 * been run of G.
 *
 * SYNOPSIS:
 * gpmakefsa [-opwa d/s] [-opgm d/s] [-ip d/s[dr]] [-silent] [-v]
 *           [-l/-h] [-diff1/-diff2] [-m maxeqns] [-mwd maxwdiffs] [-ns] [-f]
 *		groupname
 *
 * Input is from groupname.diff2
 *  (and possibly from groupname.diff1 if -c option is called).
 * Output is to groupname.wa and groupname.gm
 * Also the files groupname.diff2, and possibly groupname.diff1 if -c called,
 * may be updated by the addition of extra required word-differences.
 *
 * OPTIONS:
 * -opwa d/s	output to groupname.wa in dense or sparse format -
 *              dense is default
 * -opgm d/s	output to groupname.gm in dense or sparse format -
 *              sparse is default
 * -ip d/s[dr]  input of groupname.gm in dense or sparse format -
 *              dense is default
 * -v		verbose
 * -silent	no diagnostics
 * -l/-h	large/huge hash-tables (for big examples)
 * -diff1/diff2	(diff2 is default).
 *		This determines whether groupname.diff1 or groupname.diff2
 * 		is used to construct the word-acceptor.
 *		If diff1 is called, then the multiplier is constructed in
 *		"correction" mode. This means that
 *              if an equation is discovered which proves the word-acceptor
 *              to be incorrect, then the first word-difference machine
 *              (which should be in file groupname.diff1) is updated by making
 *              it accept this equation.
 * -m  maxeqns  Abort the multiplier checking process after finding maxeqns
 * 		offending words w (see above) -  default is MAXEQNS
 * -mwd maxwdiffs
 *              At most maxwdiffs word-differences possible.
 * -ns		Don't stop if nontrivial equation found in word-acceptor
 *		language.
 * -f           read the transition table repeatedly from file while mimimizing.
 *              this avoids storing the large table, but is a little slower.
 *
 */

#include <stdio.h>
#include "defs.h"
#include "fsa.h"
#include "rws.h"
#include "definitions.h"
#define MAXEQNS 512
#define MAXWDIFFS 2048

FILE *rfile,
     *wfile,
     *fopen();

fsa wd_fsa; /* This is for the case that we correct the diff1 machine */

void  badusage();
void  (*reduce_word)();

/* Functions defined in other files used in this file */
void  fsa_read();
fsa  *fsa_wa();
fsa  *fsa_triples();
int   fsa_checkmult();
void  fsa_print();
void  fsa_minimize();
void  fsa_ip_minimize();
void  fsa_labeled_minimize();
void  fsa_ip_labeled_minimize();
void  fsa_clear();
void  diff_reduce();
void  make_full_wd_fsa();
boolean add_wd_fsa();
void  calculate_inverses();
void fsa_table_dptr_init();

void main(argc, argv)
        int             argc;
        char           *argv[];
{ int arg, i, *inv, ngens, maxwdiffs, numeqns, old_ndiff;
  fsa *wa, *genmultptr;
  char inf1[100], inf2[100], outf1[100], outf2[100], fsaname[100],
       tempfilename[100];
  int maxeqns1, maxeqns2;
  reduction_equation *eqnptr;
  boolean diff1_ip = FALSE;
  boolean correction = FALSE;
  storage_type op_store_wa = DENSE;
  storage_type op_store_gm = SPARSE;
  storage_type ip_store = DENSE;
  int dr = 0;
  boolean eqnstop = TRUE;
  boolean idlabel = TRUE; /* always true for the moment */
  boolean readback = TRUE;
  boolean calc_wa, calc_gm;

  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  maxeqns1 = 0;
  maxeqns2 = MAXEQNS;
  maxwdiffs = MAXWDIFFS;
  inf1[0] = '\0';
  arg = 1;
  while (argc > arg) {
    if (strcmp(argv[arg],"-opwa")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      if (strcmp(argv[arg],"d")==0)
        op_store_wa = DENSE;
      else if (strcmp(argv[arg],"s")==0)
        op_store_wa = SPARSE;
      else
        badusage();
    }
    else if (strcmp(argv[arg],"-opgm")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      if (strcmp(argv[arg],"d")==0)
        op_store_gm = DENSE;
      else if (strcmp(argv[arg],"s")==0)
        op_store_gm = SPARSE;
      else
        badusage();
    }
    else if (strcmp(argv[arg],"-ip")==0) {
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
    else if (strcmp(argv[arg],"-diff1")==0)
      correction = diff1_ip = TRUE;
    else if (strcmp(argv[arg],"-diff2")==0)
      correction = diff1_ip = FALSE;
    else if (strcmp(argv[arg],"-m")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      maxeqns2 = atoi(argv[arg]);
    }
    else if (strcmp(argv[arg],"-mwd")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      maxwdiffs = atoi(argv[arg]);
    }
    else if (strcmp(argv[arg],"-ns")==0)
      eqnstop = FALSE;
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
  
  strcpy(tempfilename,inf1);
  strcat(tempfilename,"temp_fsa_XXX");

  strcpy(inf2,inf1);
  strcpy(outf1,inf1);
  strcpy(outf2,inf1);
  if (diff1_ip)
    strcat(inf1,".diff1");
  else
    strcat(inf1,".diff2");
  strcat(inf2,".diff2");
  strcat(outf1,".wa");
  strcat(outf2,".gm");

  calc_wa = TRUE;
  while (calc_wa) {
  /* First read in word-difference machine for calculating word-acceptor. */
    if ((rfile = fopen(inf1,"r")) == 0) {
      fprintf(stderr,"Cannot open file %s.\n",inf1);
      exit(1);
    }
    fsa_read(rfile,&wd_fsa,DENSE,0,0,TRUE,fsaname);
    fclose(rfile);
  
    wa = fsa_wa(&wd_fsa,op_store_wa,TRUE,tempfilename);
  
    if (print_level>1)
      printf("  #Number of states of word-acceptor before minimisation = %d.\n",
          wa->states->size);
    fsa_minimize(wa);
    if (print_level>1)
      printf("  #Number of states of word-acceptor after minimisation = %d.\n",
          wa->states->size);
  
    base_prefix(fsaname);
    strcat(fsaname,".wa");
    wfile = fopen(outf1,"w");
    fsa_print(wfile,wa,fsaname);
    fclose(wfile);

    if (print_level>0)
      printf("#Word-acceptor with %d states computed.\n",wa->states->size);

    fsa_clear(wa);

    calc_wa = FALSE;
    calc_gm = TRUE;
    while (calc_gm) {
     /* Read in word-acceptor and diff2 machine for calculating multiplier */
      if ((rfile = fopen(outf1,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",outf1);
          exit(1);
      }
      fsa_read(rfile,wa,DENSE,0,0,TRUE,fsaname);
      fclose(rfile);
      if ((rfile = fopen(inf2,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",inf2);
        exit(1);
      }
      fsa_read(rfile,&wd_fsa,DENSE,0,0,TRUE,fsaname);
      fclose(rfile);
    
      if (correction) {
	maxeqns1 = maxeqns2;
        tmalloc(eqnptr,reduction_equation,maxeqns1)
      }
      else
        eqnptr = 0;
    
      genmultptr = fsa_triples(wa,&wd_fsa,op_store_gm,TRUE,tempfilename,
				   eqnptr,maxeqns1,idlabel,eqnstop,readback);
    
      if (genmultptr==0) { /*This is the case where new equations was found. */
        if (correction)  {
          if (print_level>1)
            printf(
             "  #Altering first wd-machine to make it accept new equations.\n");
          if ((rfile = fopen(inf1,"r")) == 0) {
            fprintf(stderr,"Cannot open file %s.\n",inf1);
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
          while (eqnptr[i].lhs && i<maxeqns1) {
            add_wd_fsa(&wd_fsa,eqnptr[i],inv,FALSE);
            i++;
          }
    
          if (print_level>1)
            printf("  #First word-difference machine now has %d states.\n",
                     wd_fsa.states->size);
    
          wfile = fopen(inf1,"w");
          fsa_print(wfile,&wd_fsa,fsaname);
          fclose(wfile);
    
          tfree(inv);
          fsa_clear(&wd_fsa);
          i=0;
          while (eqnptr[i].lhs && i<maxeqns1) {
            tfree(eqnptr[i].lhs); tfree(eqnptr[i].rhs);
            i++;
          }
          tfree(eqnptr);
        }
        calc_wa = TRUE;
        tfree(wa);
        break;
      }
      if (correction)
        tfree(eqnptr);
    
      if (print_level>1)
        printf("  #Number of states of multiplier before minimisation = %d.\n",
            genmultptr->states->size);
      if (readback)
        fsa_labeled_minimize(genmultptr);
      else
        fsa_ip_labeled_minimize(genmultptr);
      if (print_level>1)
        printf("  #Number of states of multiplier after minimisation = %d.\n",
            genmultptr->states->size);
    
      base_prefix(fsaname);
      strcat(fsaname,".gm");
      wfile = fopen(outf2,"w");
      fsa_print(wfile,genmultptr,fsaname);
      fclose(wfile);

      if (print_level>0)
        printf("#General multiplier with %d states computed.\n",
                genmultptr->states->size);
    
      fsa_clear(genmultptr);
      calc_gm = FALSE;

      /* Now perform the check on the multiplier - we read it back in, since
       * the storage type has usually changed.
       */

      if ((rfile = fopen(outf2,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",inf1);
          exit(1);
      }
      fsa_read(rfile,genmultptr,ip_store,dr,0,TRUE,fsaname);
      fclose(rfile);
    
      tmalloc(eqnptr,reduction_equation,maxeqns2)
      if ((numeqns = fsa_checkmult(genmultptr,eqnptr,maxeqns2)) > 0) {
     /* A multiplier was not valid, so groupname.diff2 will need updating. */
    
        fsa_clear(genmultptr);
        if (print_level>1)
          printf(
            "  #Altering second wd-machine to make it accept new equations.\n");
        if ((rfile = fopen(inf2,"r")) == 0) {
            fprintf(stderr,"Cannot open file %s.\n",inf2);
              exit(1);
        }
        fsa_read(rfile,&wd_fsa,DENSE,0,maxwdiffs,TRUE,fsaname);
        fclose(rfile);
    /* We will use this machine for word-reduction even while correcting it!  */
        reduce_word = diff_reduce;
        fsa_table_dptr_init(&wd_fsa);
    /* We need to know the inverses of generators - let's just work them out! */
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
          printf("  #Second word-difference machine now has %d states.\n",
                   wd_fsa.states->size);
    
        wfile = fopen(inf2,"w");
        fsa_print(wfile,&wd_fsa,fsaname);
        fclose(wfile);
    
        tfree(inv);
        fsa_clear(&wd_fsa);
        for (i=0;i<numeqns;i++) {
          tfree(eqnptr[i].lhs); tfree(eqnptr[i].rhs);
        }
        calc_gm = TRUE;
        tfree(genmultptr);
      }
      else {
        if (print_level>0)
           printf("#Validity test on general multiplier succeeded.\n");
        fsa_clear(genmultptr);
      }
      tfree(eqnptr);
    }
  }

  tfree(wa);
  tfree(genmultptr);
  exit(0);
}
 
void
badusage()
{
    fprintf(stderr,
   "Usage: gpmakefsa [-opwa d/s] [-opgm d/s] [-ip d/s[dr]] [-silent] [-v]\n");
    fprintf(stderr,
      "\t\t[-l/-h] [-diff1/-diff2] [-m maxeqns] [-mwd maxwdiffs]\n");
    fprintf(stderr,
      "\t\t[-ns] [-f] groupname.\n");
    exit(1);
}
