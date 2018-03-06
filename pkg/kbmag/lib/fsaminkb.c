/* file fsaminkb.c  16. 12. 94.
 * This file contains the routines necessary to calculate a 2-variable fsa
 * accepting precisely the minimal KB-rules of a short-lex automatic group.
 *
 * The first routine fsa_minred uses the word-acceptor to construct the
 * minimal reducible words - i.e. those for which all subwords are
 * accepted by the word-acceptor - clearly it is enough for the two words
 * obtained by omitting the first and last letters to be acceptable, so
 * it is essentially an application of fsa_and.
 *
 * The second routine fsa_minkb is a 2-variable machine accepting precisely the
 * minimal KB-reduction equations - i.e it accepts a pair (w1,w2) if w1 and w2
 * are equal in G, w1 is accepted by the fsa constructed by fsa_minred, and
 * w2 is accepted by the word-accpetor. Its construciton is as in fsa_triples.
 *
 * The third routine, fsa_diff1c uses tme minkb machine to calculate the correct
 * diff1 machine. This is mainly for experimentation, but it may be the most
 * efficient machine for word-reduction in some cases.
 */

#include <stdio.h>
#include "defs.h"
#include "fsa.h"
#include "hash.h"
#include "externals.h"

/* Functions defined in this file: */
fsa *fsa_minred();
fsa *fsa_minkb();
fsa *fsa_diff1c();

extern fsa	wd_fsa;
void  (*reduce_word)();
static char     testword[4096]; /* Used for reducing words */

/* Functions used in this file and defined elsewhere */
int sparse_target();
void fsa_init();
void fsa_table_dptr_init();
void fsa_set_is_accepting();
void srec_copy();
void fsa_clear();
void hash_init();
void short_hash_init();
int  hash_locate();
int  short_hash_locate();
void hash_clear();
void short_hash_clear();
int* hash_rec();
unsigned short* short_hash_rec();
void diff_reduce();
void calculate_inverses();
int diff_no();

fsa *
fsa_minred(waptr,op_table_type,destroy,tempfilename)
	fsa *waptr;
	storage_type op_table_type;
	boolean destroy;
	char *tempfilename;
{
  int **table, ne, nsi, nsi1, ns, dr, *fsarow,
      nt, cstate, csa, csb, im, i, g, len, ct, *ht_ptr, *tab_ptr;
  boolean dense_ip, dense_op, accept;
  fsa *minred;
  hash_table ht;
  FILE *tempfile, *fopen();

  if (print_level>=3)
    printf("    #Calling fsa_minred.\n");
  if (!waptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_minred only applies to DFA's.\n");
    exit(1);
  }

  tmalloc(minred,fsa,1);
  fsa_init(minred);
  srec_copy(minred->alphabet,waptr->alphabet);
  minred->flags[DFA] = TRUE;
  minred->flags[ACCESSIBLE] = TRUE;
  minred->flags[BFS] = TRUE;

  ne = waptr->alphabet->size;
  nsi = waptr->states->size;
  nsi1 = nsi+1;

  table = waptr->table->table_data_ptr;

  dense_ip = waptr->table->table_type==DENSE;
  dr = waptr->table->denserows;
  dense_op = op_table_type==DENSE;

  minred->states->type = SIMPLE;
  minred->num_initial = 1;
  tmalloc(minred->initial,int,2);
  minred->initial[1] = 1;
  minred->table->table_type = op_table_type;
  minred->table->denserows = 0;
  minred->table->printing_format = op_table_type;
  
  hash_init(&ht,TRUE,2,0,0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = waptr->initial[1];
  ht_ptr[1] = waptr->initial[1];
  im = hash_locate(&ht,2);
  if (im!=1) {
    fprintf(stderr,"Hash-initialisation problem in fsa_binop.\n");
    exit(1);
  }
  if ((tempfile=fopen(tempfilename,"w"))==0){
    fprintf(stderr,"Error: cannot open file %s\n",tempfilename);
    exit(1);
  }
  if (dense_op)
    tmalloc(fsarow,int,ne)
  else
    tmalloc(fsarow,int,2*ne+1)
 
  cstate = 0;
  if (dense_op)
    len = ne; /* The length of the fsarow output. */
  nt = 0; /* Number of transitions in minred */

/* A state of *minred will be essentially a pair of states of *waptr,
 * with the transitions coming from those of *waptr.
 * The differences are that the left hand side will accept words  w 
 * which are not accepted by *waptr but whose maximal prefix is,
 * whereas the right hand side will accept words  w which are not accepted
 * by *waptr but whose maximal suffix is.
 * Thus, on the lhs, transitions to 0 in *waptr will go to a new accept-
 * state nsi1 instead (with no transitions from nsi1) whereas on the rhs
 * the first transition will be back to the intiial state.
 * The initial state itself is non-accept.
 */
  while (++cstate <= ht.num_recs) {
    if (print_level>=3) {
      if ((cstate<=1000 && cstate%100==0)||(cstate<=10000 && cstate%1000==0)||
          (cstate<=100000 && cstate%5000==0) || cstate%50000==0)
       printf("    #cstate = %d;  number of states = %d.\n",cstate,ht.num_recs);
    }
    ht_ptr = hash_rec(&ht,cstate);
    csa = ht_ptr[0]; csb = ht_ptr[1];
    if (!dense_op)
      len = 0;
    for (g=1;g<=ne;g++) {
/* Calculate action of generator g on state cstate */
      ht_ptr = ht.current_ptr;
      if (csa==0 || csa==nsi1)
         ht_ptr[0] = 0;
      else {
        ht_ptr[0] = target(dense_ip,table,g,csa,dr);
        if (ht_ptr[0]==0)
           ht_ptr[0] = nsi1;
      }
      if (cstate==1)
        ht_ptr[1] = csb;
      else
        ht_ptr[1] = csb ? target(dense_ip,table,g,csb,dr) : 0;
      if (ht_ptr[0]==0 || ht_ptr[1]==0)
        im = 0;
      else
        im = hash_locate(&ht,2);
      if (dense_op)
         fsarow[g-1] = im;
      else if (im>0) {
         fsarow[++len] = g;
         fsarow[++len] = im;
      }
      if (im>0)
        nt++;
    }
    if (!dense_op)
      fsarow[0] = len++;
    fwrite((void *)fsarow,sizeof(int),(size_t)len,tempfile);
  }
  fclose(tempfile);

  ns = minred->states->size = ht.num_recs;
  minred->table->numTransitions = nt;
/* Locate the accept states: first count them and then record them. */
  ct = 0;
  for (i=1;i<=ns;i++) {
    ht_ptr = hash_rec(&ht,i);
    if (ht_ptr[0]==nsi1)
        ct++;
  }
  minred->num_accepting = ct;
  if (ct==1 || ct != ns) {
    tmalloc(minred->accepting,int,ct+1);
    ct = 0;
    for (i=1;i<=ns;i++) {
      ht_ptr = hash_rec(&ht,i);
      if (ht_ptr[0]==nsi1)
        minred->accepting[++ct] = i;
    }
  }
  hash_clear(&ht);
  tfree(fsarow);
  tfree(waptr->is_accepting);

  if (destroy) 
    fsa_clear(waptr);

/* Now read the transition table back in */
  tempfile = fopen(tempfilename,"r");
  compressed_transitions_read(minred,tempfile);
  fclose(tempfile);

  unlink(tempfilename);

  return minred;
}

fsa *
fsa_minkb(minredptr,waptr,diffptr,op_table_type,destroy,tempfilename)
	fsa *minredptr, *waptr, *diffptr;
	storage_type op_table_type;
	boolean destroy;
	char *tempfilename;
/* *waptr is assumed to be the word-acceptor of an automatic group.
 * and *minredptr the fsa which accepts minimal reducible words 
 * (computed using fsa_minred above).
 * In particular, all states of *waptr should be accepting.
 * *diffptr is assumed to be a word-difference machine of the same automatic
 *  group. Both are assumed to be stored in dense-format.
 * This routine constructs the fsa of which the states are triples (s1,s2,d),
 * with s1 and s2 states of *minredptr and *waptr and d a state of *diffptr.
 * (More precisely, if *waptr has n states, then s2 may also be equal
 * to n+1, meaning that the end of string symbol has been read on  rhs.
 * Since lhs is never shorter than rhs in an accpeted pair, the end of
 * string symbolon the lhs always leads to failure.)
 * The alphabet is 2-variable with base the alphabet of *waptr
 * (i.e. the same alphabet as *diffptr).
 * The alphabet member (g1,g2) maps (s1,s2,d) to (s1^g1,s2^g2,d^(g1,g2))
 * if all three components are nonzero, and to zero otherwise.
 * The transition-table of the resulting fsa is output in the usual way to
 * file tempfilename with table-type specified by op_table_type, before
 * minimisation.
 * A state is accept is s1, s2 and d all are (but s2 always is).
 * Short hash-tables will be used, so this routine won't work if *waptr
 * or *diffptr has more than 65536 states.
 *
 */
{
  int **minredtable, **watable, ***difftable, identity, ngens, ngens1, nswa1,
      ne, ns, *fsarow,
      nt, cstate, cswa1, cswa2, csdiff, im, i, e, len, ct, *tab_ptr,
      **table, bstate, bigger;
  unsigned short *ht_ptr;
  boolean  dense_op;
  fsa *minkbptr;
  srec *labels;
  short_hash_table ht;
  FILE *tempfile, *fopen();
  char g1, g2, bg1, bg2;

  if (print_level>=3)
    printf("    #Calling fsa_minkbptr.\n");

  if (!minredptr->flags[DFA] || !waptr->flags[DFA] || !diffptr->flags[DFA]){
    fprintf(stderr,"Error: fsa__minkbptr only applies to DFA's.\n");
    exit(1);
  }
  if (minredptr->alphabet->type!=IDENTIFIERS ||
                                       waptr->alphabet->type!=IDENTIFIERS) {
    fprintf(stderr,"Error in fsa_minkbptr: an fsa has wrong type.\n");
    exit(1);
  }
  if (waptr->num_accepting != waptr->states->size) {
    fprintf(stderr,
       "Error in fsa_minkbptr: second fsa should be a word-acceptor.\n");
    exit(1);
  }
  if (diffptr->alphabet->type!=PRODUCT || diffptr->alphabet->arity!=2) {
    fprintf(stderr,"Error in fsa_minkbptr: third fsa must be 2-variable.\n");
    exit(1);
  }
  if (diffptr->states->type!=WORDS) {
    fprintf(stderr,
       "Error in fsa_minkbptr: third fsa must be word-difference type.\n");
    exit(1);
  }
  if (!srec_equal(diffptr->alphabet->base,waptr->alphabet)) {
    fprintf(stderr,"Error in fsa_minkbptr: fsa's alphabet's don't match.\n");
    exit(1);
  }
  if (minredptr->states->size>MAXSHORT || waptr->states->size>MAXSHORT
                                       || diffptr->states->size>MAXSHORT) {
    fprintf(stderr,
       "Error in fsa_minkbptr: one of the fsa's has too many states.\n");
    exit(1);
  }

  fsa_table_dptr_init(diffptr);

  tmalloc(minkbptr,fsa,1);
  fsa_init(minkbptr);
  srec_copy(minkbptr->alphabet,diffptr->alphabet);
  minkbptr->flags[DFA] = TRUE;
  minkbptr->flags[ACCESSIBLE] = TRUE;
  minkbptr->flags[BFS] = TRUE;

  ngens = waptr->alphabet->size;
  ngens1 = ngens+1;
  ne = diffptr->alphabet->size;
  nswa1 = waptr->states->size + 1;

  if (ne != ngens1*ngens1-1) {
   fprintf(stderr,
       "Error: in a 2-variable fsa, alphabet size should = ngens^2 - 1.\n");
    exit(1);
  }

  minredtable = minredptr->table->table_data_ptr;
  watable = waptr->table->table_data_ptr;
  difftable = diffptr->table->table_data_dptr;

  dense_op = op_table_type==DENSE;

  minkbptr->states->type = SIMPLE;

  minkbptr->num_initial = 1;
  tmalloc(minkbptr->initial,int,2);
  minkbptr->initial[1] = 1;
  minkbptr->table->table_type = op_table_type;
  minkbptr->table->denserows = 0;
  minkbptr->table->printing_format = op_table_type;
  
  short_hash_init(&ht,TRUE,3,0,0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = minredptr->initial[1];
  ht_ptr[1] = waptr->initial[1];
  ht_ptr[2] = identity = diffptr->initial[1];
  im = short_hash_locate(&ht,3);
  if (im!=1) {
    fprintf(stderr,"Hash-initialisation problem in fsa_minkbptr.\n");
    exit(1);
  }
  if ((tempfile=fopen(tempfilename,"w"))==0){
    fprintf(stderr,"Error: cannot open file %s\n",tempfilename);
    exit(1);
  }
  if (dense_op)
    tmalloc(fsarow,int,ne)
  else
    tmalloc(fsarow,int,2*ne+1)
 
  cstate = 0;
  if (dense_op)
    len = ne; /* The length of the fsarow output. */
  nt = 0; /* Number of transitions in minkbptr */

  while (++cstate <= ht.num_recs) {
    if (print_level>=3) {
      if ((cstate<=1000 && cstate%100==0)||(cstate<=10000 && cstate%1000==0)||
          (cstate<=100000 && cstate%5000==0) || cstate%50000==0)
       printf("    #cstate = %d;  number of states = %d.\n",cstate,ht.num_recs);
    }
    ht_ptr = short_hash_rec(&ht,cstate);
    cswa1 = ht_ptr[0]; cswa2 = ht_ptr[1];
    csdiff = ht_ptr[2];
    if (dense_op)
      for (i=0;i<ne;i++) fsarow[i] = 0;
    else
      len = 0;
    e = 0; /* e is the number of the edge corresponding to the pair (g1,g2) */
    for (g1=1;g1<=ngens;g1++) for (g2=1;g2<=ngens1;g2++) {
      e++;
/* Calculate action of generator-pair (g1,g2) on state cstate - since the
 * padding symbol cannot occur on the lhs, g1 only goes up to ngens.
 */
      ht_ptr = ht.current_ptr;
      ht_ptr[2] = dense_dtarget(difftable,g1,g2,csdiff);
      if (ht_ptr[2]==0)
        im=0;
      else {
        ht_ptr[0] = dense_target(minredtable,g1,cswa1);
        if (ht_ptr[0]==0)
          im=0;
        else {
          ht_ptr[1] = g2==ngens1 ? nswa1 : cswa2==nswa1 ? 0 :
                                          dense_target(watable,g2,cswa2);
          if (ht_ptr[1]==0)
            im=0;
          else
            im = short_hash_locate(&ht,3);
        }
      }

      if (dense_op)
         fsarow[e-1] = im;
      else if (im>0) {
         fsarow[++len] = e;
         fsarow[++len] = im;
      }
      if (im>0)
        nt++;
    }  /* for (g1=1;g1<=ngens1; ... */
    if (!dense_op)
      fsarow[0] = len++;
    fwrite((void *)fsarow,sizeof(int),(size_t)len,tempfile);
  }  /*while (++cstate <= ht.num_recs) */
  fclose(tempfile);

  ns = minkbptr->states->size = ht.num_recs;
  minkbptr->table->numTransitions = nt;

  if (print_level>=3)
    printf("    #Calculated transitions - %d states, %d transitions.\n",ns,nt);

/* Now locate the accept states */

  fsa_set_is_accepting(minredptr);
  ct = 0;
  for (i=1;i<=ns;i++) {
    ht_ptr = short_hash_rec(&ht,i);
    if (minredptr->is_accepting[ht_ptr[0]] && ht_ptr[2]==identity)
      ct++;
  }
  minkbptr->num_accepting = ct;
  if (ct==1 || ct != ns) {
    tmalloc(minkbptr->accepting,int,ct+1);
    ct = 0;
    for (i=1;i<=ns;i++) {
      ht_ptr = short_hash_rec(&ht,i);
    if (minredptr->is_accepting[ht_ptr[0]] && ht_ptr[2]==identity)
        minkbptr->accepting[++ct] = i;
    }
  }

  short_hash_clear(&ht);
  tfree(minredptr->is_accepting);
  tfree(fsarow);
  if (destroy) {
    fsa_clear(minredptr); fsa_clear(waptr); fsa_clear(diffptr);
  }
/* Now read the transition table back in */
  tempfile = fopen(tempfilename,"r");
  compressed_transitions_read(minkbptr,tempfile);
  fclose(tempfile);

  unlink(tempfilename);

  return minkbptr;
}

fsa *
fsa_diff1c(minkbptr,op_table_type)
	fsa  *minkbptr;
        storage_type op_table_type;
/* *minkbptr should be the minimal KB-reduction machine output by fsa_minkb.
 * The correct word-difference machine for the minimal set of KB-reduction
 * equations is calculated from this machine. Note that wd_fsa must be defined
 * externally to enable word-reduction.
 */
{ int i, n, ngens, ngens1, ne, ns, g1, g2, minkb_state, diff1c_state,
  *state_diff, ***minkb_table, ***wd_table, **table, minkb_im, diff1c_im, *inv;
  fsa *diff1cptr;

  if (print_level>=3)
    printf("    #Calling fsa_diff1c.\n");
  if (!minkbptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_diff1c only applies to DFA's.\n");
    exit(1);
  }

  if (minkbptr->alphabet->type != PRODUCT || minkbptr->alphabet->arity != 2) {
    fprintf(stderr,"Error in fsa_diff1c: fsa must be 2-variable.\n");
    exit(1);
  }

  ne = minkbptr->alphabet->size;
  ngens = minkbptr->alphabet->base->size;
  ngens1 = ngens+1;
  ns = minkbptr->states->size;
  fsa_table_dptr_init(minkbptr);
  minkb_table = minkbptr->table->table_data_dptr;

  if (ne != ngens1*ngens1-1) {
   fprintf(stderr,
       "Error: in a 2-variable fsa, alphabet size should = ngens^2 - 1.\n");
    exit(1);
  }

  tmalloc(diff1cptr,fsa,1);
  fsa_init(diff1cptr);
  srec_copy(diff1cptr->alphabet,minkbptr->alphabet);

/* Maximum possible number of states of diff1c is of course ns */
  diff1cptr->states->type = WORDS;
  tmalloc(diff1cptr->states->names,char *,ns+1);
  diff1cptr->states->alphabet_size = minkbptr->alphabet->base->size;
  for (i=1;i<=minkbptr->alphabet->base->size;i++) {
    tmalloc(diff1cptr->states->alphabet[i],char,
                                 strlen(minkbptr->alphabet->base->names[i])+1);
    strcpy(diff1cptr->states->alphabet[i],minkbptr->alphabet->base->names[i]);
  }
  diff1cptr->states->size = 1;
/* Set up first state, which is the empty word. */
  tmalloc(diff1cptr->states->names[1],char,1);
  diff1cptr->states->names[1][0] = '\0';

  diff1cptr->num_initial = 1;
  tmalloc(diff1cptr->initial,int,2);
  diff1cptr->initial[1] = 1;
  diff1cptr->num_accepting = 1;
  tmalloc(diff1cptr->accepting,int,2);
  diff1cptr->accepting[1] = 1; /* Only the initial state is accepting */

  diff1cptr->flags[DFA] = TRUE;
  diff1cptr->flags[TRIM] = TRUE;

  fsa_table_init(diff1cptr->table,ns,diff1cptr->alphabet->size);
  diff1cptr->table->printing_format = op_table_type;
  fsa_table_dptr_init(diff1cptr);
  wd_table = diff1cptr->table->table_data_dptr;
  table = diff1cptr->table->table_data_ptr;
  for (i=1;i<=ne;i++)
    table[i][1] = 0;

  reduce_word = diff_reduce;
  calculate_inverses(&inv,ngens);

/* Now build the transition-table of diff1cptr, using that of minkbptr.
 * Each state of minkbptr maps onto one of *diff1cptr, with the mapping
 * store din the list state_diff.
 */
  tmalloc(state_diff,int,ns+1);
  for (i=1;i<=ns;i++)
    state_diff[i] = 0;
  state_diff[1] = 1;
  for (minkb_state=1;minkb_state<=ns;minkb_state++) {
    diff1c_state = state_diff[minkb_state];
    for (g1=1;g1<=ngens1;g1++) for (g2=1;g2<=ngens1;g2++) {
      if (g1==ngens1 && g2==ngens1)
        continue;
      minkb_im = dense_dtarget(minkb_table,g1,g2,minkb_state);
      if (minkb_im==0)
        continue;
      diff1c_im = state_diff[minkb_im];
      if (diff1c_im==0) {
/* We have to work out what word-difference the state maps onto */
        if (g1==ngens1)
          sprintf(testword,"%s%c",diff1cptr->states->names[diff1c_state],g2);
        else if (g2==ngens1)
          sprintf(testword,"%c%s",
                                inv[g1],diff1cptr->states->names[diff1c_state]);
        else
          sprintf(testword,"%c%s%c",
                             inv[g1],diff1cptr->states->names[diff1c_state],g2);
        (*reduce_word)(testword);
        n=diff_no(diff1cptr,testword);
        if (n==0) {
          n = ++diff1cptr->states->size;
          tmalloc(diff1cptr->states->names[n],char,strlen(testword)+1);
          strcpy(diff1cptr->states->names[n],testword);
          for (i=1;i<=ne;i++)
            table[i][n] = 0;
        }
        state_diff[minkb_im] = diff1c_im = n;
      }
      set_dense_dtarget(wd_table,g1,g2,diff1c_state,diff1c_im);
    }
  }

  tfree(inv);
  tfree(state_diff);
  return diff1cptr;
}
