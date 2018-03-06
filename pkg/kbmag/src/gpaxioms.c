/* gpaxioms.c   14.12.94. 
 *
 * Carry out the axiom-checking for a short-lex automatic group.
 *
 * SYNOPSIS:
 * gpaxioms [-ip d/s[dr]] [-op d/s] [-silent] [-v] [-l] [-f] [-s] groupname
 *
 * OPTIONS:
 * -ip d/s[dr]  input in dense or sparse format - dense is default
 * -op d/s      output in dense or sparse format - sparse is default
 * -v           verbose
 * -silent      no diagnostics
 * -l/-h        large/huge hash-tables (for big examples)
 * -f           read the transition table repeatedly from file while mimimizing,
 *              in calls of fsa_genmult2 and fsa_composite.
 *              This avoids storing large tables, but is a little slower.
 * -s		Throw away files immediately after use whenever possible, to
 *              save disk-space. This will slow things down considerably.
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
fsa		wd_fsa;
rewriting_system  rws;
void (*reduce_word)(); /* This is the function that will be passed to
                        * worddiff.c for reduction when calculating
                        * the word-difference table.
                        */

static char gpname[100], inf[100], outf[100], fsaname[100],
       tablefilename[100];
static fsa genmult, genmult2, *genmult2ptr, mult1, mult2, *compmult;
static int ngens, neqns, *inv, *newgenno, nlabgm2, numstoredmult;
static reduction_equation *eqn;
static char *lhs, *rhs, **labgm2, **storedmult;
static boolean allshort;
static boolean  redundantgen = FALSE;
static storage_type ip_store = DENSE;
static int dr = 0;
static storage_type op_store = SPARSE;
static boolean readback = TRUE;
static boolean keepfiles = TRUE;
FILE *rfile, *wfile;

/* Functions defined in this file */
void redundant_gen_check();
void balance_equations();
void check_short_relation();
char *file_suffix();
void long_word_multiplier();
void badusage();

/* Functions used in this file defined in other files: */
void read_kbinput_simple();
void printbuffer();
void add_to_buffer();
int add_word_to_buffer();
int int_len();
void fsa_table_dptr_init();
fsa  *fsa_genmult2();
void  fsa_makemult();
void  fsa_makemult2();
fsa  *fsa_composite();
void  diff_reduce();
boolean fsa_equal();
void fsa_clear();
void rws_clear();

void main(argc, argv)
        int             argc;
        char           *argv[];
{ int arg, i, j, l, ct;
  char *suffl, *suffr;
  boolean gotl, gotr;

  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  gpname[0] = '\0';
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
    else if (strcmp(argv[arg],"-s")==0)
      keepfiles = FALSE;
    else {
       if (argv[arg][0] == '-')
         badusage();
       if (strcmp(gpname,"")!=0)
         badusage();
       strcpy(gpname,argv[arg]);
    }
    arg++;
  }
  if (strlen(gpname)==0)
    badusage();
  
  strcpy(tablefilename,gpname);
  strcat(tablefilename,"temp_axXXX");

/* First read in the defining relations for the group. */
  strcpy(inf,gpname);
  if ((rfile = fopen(inf,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf);
      exit(1);
  }
  read_kbinput_simple(rfile,TRUE);
  fclose(rfile);
  ngens = rws.num_gens;
  neqns = rws.num_eqns;
  inv = rws.inv_of;
  eqn = rws.eqns;

/* We now see if any of the generators have become redundant, and record
 * any such equalities or trivialities in the array newgenno
 * We use the diff2 word-reduction machine for reducing words.
 */
  redundant_gen_check();

/* Now remove any redundant generators from the relations, and balance them. */
  balance_equations();
             
/* See if all relations are short (length(lhs)=2 ) */
  allshort = TRUE;
  for (i=1;i<=neqns;i++) if (strlen(eqn[i].lhs) > 2) {
    allshort = FALSE;
    break;
  }

/* Now calculate the general multiplier for words of length 2.
 * If allshort is true, we don't need the transitions - only the state labels.
 */
  strcpy(inf,gpname);
  strcat(inf,".gm");
  if ((rfile = fopen(inf,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf);
      exit(1);
  }
  fsa_read(rfile,&genmult,ip_store,dr,0,TRUE,fsaname);
  fclose(rfile);

  if (allshort) {
    genmult2ptr = fsa_genmult2(&genmult,op_store,TRUE,"",FALSE);
    if (print_level>1)
      printf(
           "  #Number of states of genmult2 = %d.\n",genmult2ptr->states->size);
  }
  else {
    genmult2ptr = fsa_genmult2(&genmult,op_store,TRUE,tablefilename,readback);
    if (print_level>1)
      printf(
           "  #Number of states of genmult2 = %d.\n",genmult2ptr->states->size);
  
    if (readback)
      fsa_labeled_minimize(genmult2ptr);
    else
      fsa_ip_labeled_minimize(genmult2ptr);
    if (print_level>1)
      printf("  #Number of states of genmult2 after minimization = %d.\n",
               genmult2ptr->states->size);
    strcpy(fsaname,rws.name);
    strcat(fsaname,".gm2");
    strcpy(outf,gpname);
    strcat(outf,".gm2");
    wfile = fopen(outf,"w");
    fsa_print(wfile,genmult2ptr,fsaname);
    if (print_level>0)
      printf("#General length-2 multiplier with %d states computed.\n",
            genmult2ptr->states->size);
    fclose(wfile);
  }

/* Now we check axioms for inverse-relators and short-relators, by simply
 * checking that the state-labels for the lhs and rhs correspond in
 * *genmult2ptr.
 */
  if (print_level>0) {
    printf("#Checking inverse and short relations.\n");
  }
  nlabgm2 = genmult2ptr->states->labels->size;
  labgm2 = genmult2ptr->states->labels->names;
  lhs = rws.testword1; rhs=rws.testword2;
  for (i=1;i<=ngens;i++) if (newgenno[i]==i) {
    lhs[0]=i; lhs[1]=inv[i]; lhs[2]='\0'; rhs[0]='\0';
    check_short_relation();
  }
  for (i=1;i<=neqns;i++) if (strlen(eqn[i].lhs)==2) {
    lhs = eqn[i].lhs; rhs = eqn[i].rhs;
    check_short_relation();
  }
  fsa_clear(genmult2ptr);

/* Finally we deal with the other relations, if any */
  if (!allshort) {
/* If keepfiles is false, then we will throw away every multiplier immediately
 * after it has been used. Otherwise, we will keep them, in case they are
 * needed again.
 * We store a list of words for which we have got the multipliers in
 * storedmult. We first form a rough upper bound on how long this list
 * could get - ngens + total relator length - 1.
 */
    strcpy(fsaname,rws.name);
    strcat(fsaname,".mult"); /* this is unimportant, since file is temporary */
    if (keepfiles) {
      ct = ngens;
      for (i=1;i<=neqns;i++) if (strlen(eqn[i].lhs)>2)
        ct += (strlen(eqn[i].lhs) + strlen(eqn[i].rhs));
      tmalloc(storedmult,char *,ct);
      numstoredmult = 0;
    }

    for (i=1;i<=neqns;i++) if (strlen(eqn[i].lhs)>2) {
      lhs = eqn[i].lhs; rhs = eqn[i].rhs;
      if (print_level > 0) {
        buffer[0] = '\0';
        add_to_buffer(0,"#Checking relation:  ");
        add_word_to_buffer(stdout,lhs,rws.gen_name);
        add_to_buffer(0," = ");
        if ((l=strlen(buffer))>44) {
           printbuffer(stdout,buffer);
           add_to_buffer(21,"");
        }
        add_word_to_buffer(stdout,rhs,rws.gen_name);
        printbuffer(stdout,buffer);
      }
      suffl = file_suffix(lhs); suffr = file_suffix(rhs);
      gotl = gotr = FALSE;
      if (keepfiles) {
/* Check to see if we have either of these multipliers already */
        for (j=1;j<=numstoredmult;j++) {
          if (strcmp(suffl,storedmult[j])==0)
            gotl = TRUE;
          if (strcmp(suffr,storedmult[j])==0)
            gotr = TRUE;
        }
      }
      if (!gotl)
        long_word_multiplier(lhs,suffl);
      if (!gotr)
        long_word_multiplier(rhs,suffr);
/* Read in the two multipliers and compare them */
      sprintf(inf,"%s.m%s",gpname,suffl);
      if ((rfile = fopen(inf,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",inf);
          exit(1);
      }
      fsa_read(rfile,&mult1,ip_store,0,0,TRUE,fsaname);
      fclose(rfile);
      sprintf(inf,"%s.m%s",gpname,suffr);
      if ((rfile = fopen(inf,"r")) == 0) {
        fprintf(stderr,"Cannot open file %s.\n",inf);
          exit(1);
      }
      fsa_read(rfile,&mult2,ip_store,0,0,TRUE,fsaname);
      fclose(rfile);
      if (!fsa_equal(&mult1,&mult2)) {
        buffer[0] = '\0';
        add_to_buffer(0,"#Relation check fails for relation:  ");
        add_word_to_buffer(stdout,lhs,rws.gen_name);
        add_to_buffer(0," = ");
        add_word_to_buffer(stdout,rhs,rws.gen_name);
        printbuffer(stdout,buffer);
        exit(2);
      }
      if (keepfiles) {
        if (gotl) tfree(suffl)
        else storedmult[++numstoredmult] = suffl;
        if (gotr) tfree(suffr)
        else storedmult[++numstoredmult] = suffr;
      }
      else {
        sprintf(inf,"%s.m%s",gpname,suffl); unlink(inf);
        sprintf(inf,"%s.m%s",gpname,suffr); unlink(inf);
        tfree(suffl); tfree(suffr);
      }
    }
    if (keepfiles) {
      for (i=1;i<=numstoredmult;i++) {
        sprintf(outf,"%s.m%s",gpname,storedmult[i]);
        unlink(outf);
        tfree(storedmult[i]);
      }
      tfree(storedmult);
    }
    fsa_clear(&mult1); fsa_clear(&mult2);
    strcpy(outf,gpname);
    strcat(outf,".gm2");
    unlink(outf);
  }
  tfree(genmult2ptr);
  tfree(newgenno);
  rws_clear(&rws);
  if (print_level>0)
    printf("#Axiom checking succeeded.\n");
  exit(0);
}
 
void
redundant_gen_check()
{ int i,j;
  if (print_level > 1)
    printf("  #Checking if any generators are redundant.\n");
  tmalloc(newgenno,int,ngens+1);
  for (i=1;i<=ngens;i++)
    newgenno[i] = i;
  strcpy(inf,gpname);
  strcat(inf,".diff2");
  if ((rfile = fopen(inf,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf);
      exit(1);
  }
  fsa_read(rfile,&wd_fsa,DENSE,0,0,TRUE,fsaname);
  fsa_table_dptr_init(&wd_fsa);
  fclose(rfile);
  reduce_word = diff_reduce;
  for (i=1;i<=ngens;i++) {
    rws.testword1[0] = i; rws.testword1[1] = '\0';
    reduce_word(rws.testword1);
    j = rws.testword1[0]; 
    if (j == '\0') {
      if (print_level>0)
        printf("#Generator number %d is trivial.\n",i);
      newgenno[i] = 0;
      redundantgen = TRUE;
    }
    else if (j < i) {
      if (print_level>1)
        printf("  #Generator number %d equals generator number %d.\n",i,j);
      newgenno[i] = j;
      redundantgen = TRUE;
    }
  }
  fsa_clear(&wd_fsa);
  return;
}

void
balance_equations()
{ int i, j, l1, l2;
  char *lhs, *rhs, *ptr1, *ptr2;
  if (print_level > 1)
    printf("  #Simplifying and balancing inverse-list and relations.\n");
  if (redundantgen) for (i=1;i<=ngens;i++) 
    inv[i] = newgenno[inv[i]];

  for (i=1;i<=neqns;i++) {
    lhs = eqn[i].lhs; rhs = eqn[i].rhs;
    l1 = strlen(lhs); l2 = strlen(rhs);
    if (redundantgen) {
      ptr1 = lhs;
      while (*ptr1 != '\0') {
        if (newgenno[*ptr1]==0) {
          ptr2 = ptr1;
          do { *ptr2 = *(ptr2+1); ptr2++;
          } while (*ptr2 != '\0');
          l1--;
        }
        else { *ptr1 = newgenno[*ptr1]; ptr1++;
        }
      }
      ptr1 = rhs;
      while (*ptr1 != '\0') {
        if (newgenno[*ptr1]==0) {
          ptr2 = ptr1;
          do { *ptr2 = *(ptr2+1); ptr2++;
          } while (*ptr2 != '\0');
          l2--;
        }
        else { *ptr1 = newgenno[*ptr1]; ptr1++;
        }
      }
    }
    if (l1+l2 <= 2) { /*We no longer need keep this equation */
      if (print_level > 1)
        printf("  #Equation number %d is being discarded.\n",i);
      lhs[0] = rhs[0] = '\0';
    }
    else { /* balance the lengths of the lhs and rhs. */
      if (l2 > l1) { /* swap lhs and rhs */
        eqn[i].lhs = rhs; eqn[i].rhs = lhs;
        lhs = eqn[i].lhs; rhs = eqn[i].rhs;
        l1 = strlen(lhs); l2 = strlen(rhs);
      }
      while (l1-l2 > 1) { /* move a generator from lhs to rhs */
        rhs[l2] = inv[lhs[l1-1]]; l2++; rhs[l2] = '\0';
        l1--; lhs[l1] = '\0';
      }
    }
  }
}

void
check_short_relation()
/* Check that the general multiplier automaton *genmult2ptr for words of
 * length two satisfies the short-equation lhs = rhs.
 * This is done merely by looking at the state-labels of *genmult2ptr.
 * The transitions are not needed.
 */
{ int i, j, l;
  char *lab;
  boolean foundlhs, foundrhs; 
  if (print_level > 1) {
    buffer[0] = '\0';
    add_to_buffer(0,"  #Checking short relation:  ");
    add_word_to_buffer(stdout,lhs,rws.gen_name);
    add_to_buffer(0," = ");
    add_word_to_buffer(stdout,rhs,rws.gen_name);
    printbuffer(stdout,buffer);
  }
/* pad rhs to length 2 if necessary */
  if (strlen(rhs)==0){
    rhs[0]=rhs[1]=ngens+1; rhs[2]='\0';
  }
  else if (strlen(rhs)==1){
    rhs[1]=ngens+1; rhs[2]='\0';
  }
/* Search through the state-labels of genmult2ptr.
 * For the relation to be satisfied by the fsa, it must be the case that
 * lhs occurs in a label <=> rhs occurs in that label.
 */
  for (i=1;i<=nlabgm2;i++) {
    foundlhs = foundrhs = FALSE;
    lab = labgm2[i];
    l = strlen(lab);
    for (j=0;j<l;j+=2) {
      if (lab[j]==lhs[0] && lab[j+1]==lhs[1])
        foundlhs = TRUE;
      if (lab[j]==rhs[0] && lab[j+1]==rhs[1])
        foundrhs = TRUE;
    }
    if (foundlhs != foundrhs) {
      buffer[0] = '\0';
      add_to_buffer(0,"#Relation check fails for short relation:  ");
      add_word_to_buffer(stdout,lhs,rws.gen_name);
      add_to_buffer(0," = ");
      add_word_to_buffer(stdout,rhs,rws.gen_name);
      printbuffer(stdout,buffer);
      sprintf(buffer,"Label number %d:  ",i);
      add_word_to_buffer(stdout,lab,rws.gen_name);
      printbuffer(stdout,buffer);
      exit(2);
    }
  }
}

char *
file_suffix(w)
	char *w;
/* For a word w in the generators, this function returns a corresponding
 * string with the terms of w replaced by integers separated by '_'.
 * This is used as a suffix in the filenames used for storing the
 * corresponding multiplier fsa's.
 */
{ char *s, *p;
  boolean first;
  int len;
/* First work out the length of the required string. */
  len = strlen(w);
  p = w-1;
  while (*(++p) != 0)
    len += int_len((int)(*p));
  tmalloc(s,char,len);
  s[0] = '\0';
  first = TRUE;
  p = w-1;
  while (*(++p) != 0) {
    if (first)
      first = FALSE;
    else
      sprintf(s+strlen(s),"_");
    sprintf(s+strlen(s),"%d",*p);
  }
  return s;
}

void
long_word_multiplier(w,s)
	char *w, *s;
/* Calculate the multiplier associated with the word w.
 * s is the suffix of the file in which it will be stored.
 * (s has been derived from w by a call of file_suffix).
 */
{ int i, l;
  char *wl, *wlt, *wr, *wrt, *suffl, *sufflt, *suffr, *suffrt;
  boolean gotl, gotr, gotlt, gotrt;
  if (print_level >= 3) {
    buffer[0] = '\0';
    add_to_buffer(0,"  #Calculating multiplier for word:  ");
    add_word_to_buffer(stdout,w,rws.gen_name);
    printbuffer(stdout,buffer);
  }
  l = strlen(w);

  if (l==1) { /* Length 1 - use fsa_makemult */
    strcpy(inf,gpname);
    strcat(inf,".gm");
    if ((rfile = fopen(inf,"r")) == 0) {
      fprintf(stderr,"Cannot open file %s.\n",inf);
        exit(1);
    }
    fsa_read(rfile,&genmult,op_store,0,0,TRUE,fsaname);
    fclose(rfile);
    fsa_makemult(&genmult,w[0]);
    fsa_minimize(&genmult);
    sprintf(outf,"%s.m%s",gpname,s);
    wfile = fopen(outf,"w");
    fsa_print(wfile,&genmult,fsaname);
    fclose(wfile);
    fsa_clear(&genmult);
  }
  else if (l==2) { /* Length 2 - use fsa_makemult */
    strcpy(inf,gpname);
    strcat(inf,".gm2");
    if ((rfile = fopen(inf,"r")) == 0) {
      fprintf(stderr,"Cannot open file %s.\n",inf);
        exit(1);
    }
    fsa_read(rfile,&genmult2,op_store,0,0,TRUE,fsaname);
    fclose(rfile);
    fsa_makemult2(&genmult2,w[0],w[1]);
    fsa_minimize(&genmult2);
    sprintf(outf,"%s.m%s",gpname,s);
    wfile = fopen(outf,"w");
    fsa_print(wfile,&genmult2,fsaname);
    fclose(wfile);
    fsa_clear(&genmult2);
  }
  else { /* general case - we have to split w up */
    if (l%2 == 0) {
      tmalloc(wl,char,l/2 + 1); tmalloc(wr,char,l/2 + 1);
      for (i=0;i<l/2;i++) wl[i]=w[i];
      wl[l/2]='\0';
      strcpy(wr,w+l/2);
      suffl = file_suffix(wl); suffr = file_suffix(wr);
    }
    else  {
      tmalloc(wl,char,l/2 + 2); tmalloc(wr,char,l/2 + 1);
      for (i=0;i<=l/2;i++) wl[i]=w[i];
      wl[l/2 + 1]='\0';
      strcpy(wr,w+l/2+1);
      suffl = file_suffix(wl); suffr = file_suffix(wr);
    }
/* See whether we have either of them already */
    gotl = gotr = FALSE;
    if (keepfiles) {
      for (i=1;i<=numstoredmult;i++) {
        if (strcmp(suffl,storedmult[i])==0)
          gotl = TRUE;
        if (strcmp(suffr,storedmult[i])==0)
          gotr = TRUE;
      }
    }

    if (keepfiles && l%2==1 && (!gotl || !gotr)) {
/* In this case, there are two possible ways to split w up -
 * we see if the other way has more multipliers already stored.
 */
      tmalloc(wlt,char,l/2 + 1); tmalloc(wrt,char,l/2 + 2);
      for (i=0;i<l/2;i++) wlt[i]=w[i];
      wlt[l/2]='\0';
      strcpy(wrt,w+l/2);
      sufflt = file_suffix(wlt); suffrt = file_suffix(wrt);
      gotlt = gotrt = FALSE;
      for (i=1;i<=numstoredmult;i++) {
        if (strcmp(sufflt,storedmult[i])==0)
          gotlt = TRUE;
        if (strcmp(suffrt,storedmult[i])==0)
          gotrt = TRUE;
      }
      if ((gotlt && gotrt) || ((gotlt || gotrt) && !gotl && !gotr)) {
        tfree(wl); tfree(wr); tfree(suffl); tfree(suffr);
        wl=wlt; wr=wrt; suffl=sufflt; suffr=suffrt; gotl=gotlt; gotr=gotrt;
      }
      else {
        tfree(wlt); tfree(wrt); tfree(sufflt); tfree(suffrt);
      }
    }
    if (!gotl)
      long_word_multiplier(wl,suffl);
    if (!gotr && strcmp(wl,wr)!=0)
      long_word_multiplier(wr,suffr);
/* Read back in the two multipliers and form their composite */
    sprintf(inf,"%s.m%s",gpname,suffl);
    if ((rfile = fopen(inf,"r")) == 0) {
      fprintf(stderr,"Cannot open file %s.\n",inf);
        exit(1);
    }
    fsa_read(rfile,&mult1,ip_store,dr,0,TRUE,fsaname);
    fclose(rfile);
    sprintf(inf,"%s.m%s",gpname,suffr);
    if ((rfile = fopen(inf,"r")) == 0) {
      fprintf(stderr,"Cannot open file %s.\n",inf);
        exit(1);
    }
    fsa_read(rfile,&mult2,ip_store,dr,0,TRUE,fsaname);
    fclose(rfile);

    compmult =
        fsa_composite(&mult1,&mult2,op_store,TRUE,tablefilename,readback);
    if (readback)
      fsa_minimize(compmult);
    else
      fsa_ip_minimize(compmult);
    sprintf(outf,"%s.m%s",gpname,s);
    wfile = fopen(outf,"w");
    fsa_print(wfile,compmult,fsaname);
    fclose(wfile);
    fsa_clear(compmult);
    tfree(compmult);

    if (keepfiles) {
      if (gotl) tfree(suffl)
      else storedmult[++numstoredmult] = suffl;
      if (gotr) tfree(suffr)
      else storedmult[++numstoredmult] = suffr;
    }
    else {
      sprintf(inf,"%s.m%s",gpname,suffl); unlink(inf);
      sprintf(inf,"%s.m%s",gpname,suffr); unlink(inf);
      tfree(suffl); tfree(suffr);
    }
    tfree(wl); tfree(wr);
  }


  return;
}

void
badusage()
{
   fprintf(stderr,
"Usage: gpaxioms [-ip d/s[dr]] [-op d/s] [-silent] [-v] [-l] [-f] groupname.\n"
  );
	exit(1);
}
