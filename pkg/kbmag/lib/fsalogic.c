/* file fsalogic.c  2. 11. 94.
 * 18.11.94. - added fsa_and_not
 * This file contains routines for performing logical operations on
 * finite state automata.
 * The functions in this file currently only deal with deterministic fsa's
 *
 * The general methodology of some functions in this file, such as
 * fsa_and, fsa_or and fsa_exists, and of several functions in different files
 * is as follows.
 * These functions construct a particular finite state automaton, of which
 * the states are some kind of sets of states of the input automata,
 * and the transition table is generated in order.
 * The states are stored as records in a hash-table, so that they can be
 * located easily.
 * The initial state is constructed by hand and numbered 1, and then the
 * transitions from each state are generated in order. If a new state is
 * generated as target of a transition, then it is added to the end of the
 * list of states. Eventually the process will terminate when the table closes,
 * and all transitions have been generated from each state.
 * The transitions from a particular state will not be needed again until
 * the generation process is complete, and so they are output(in unformatted
 * form) to a file (with name tempfilename).
 * When the table is complete, the hash-table data can normally be discarded
 * (once it has been used to identify accept-states in the constructed fsa),
 * and then the transition table can be read back in.
 * There is also an option to clear the space taken by the input automata
 * before re-reading.
 * This results in saving of core memory, because the transition table of
 * the fsa being constructed and the data in the hash-table are never held
 * in main memory at the same time.
 */

#include <stdio.h>
#include "defs.h"
#include "fsa.h"
#include "hash.h"
#include "externals.h"

typedef enum {AND, OR, AND_NOT} binop;

/* Functions defined in this file: */
boolean srec_equal();
boolean table_equal();
boolean fsa_equal();
fsa *fsa_and();
fsa *fsa_or();
fsa *fsa_and_not();
fsa *fsa_binop();
fsa *fsa_or();
fsa *fsa_not();
fsa *fsa_exists();
fsa *fsa_exists_short();
fsa *fsa_exists_int();
fsa *fsa_greater_than();

/* Functions used in this file and defined elsewhere */
int sparse_target();
void fsa_init();
void fsa_table_init();
void fsa_table_dptr_init();
void fsa_set_is_initial();
void fsa_set_is_accepting();
void srec_copy();
void fsa_copy();
void fsa_clear();
void hash_init();
void short_hash_init();
int  hash_locate();
int  short_hash_locate();
void hash_clear();
void short_hash_clear();
int* hash_rec();
unsigned short* short_hash_rec();
int hash_rec_len();
int short_hash_rec_len();

boolean
srec_equal(srptr1,srptr2)
	srec *srptr1, *srptr2;
/* Test equality of set records *srptr1 and *srptr2 */
{ int i;
  if (srptr1->type != srptr2->type)
    return FALSE;
  if (srptr1->size != srptr2->size)
    return FALSE;
  if (srptr1->type==IDENTIFIERS || srptr1->type==STRINGS ||
                                   srptr1->type==WORDS) {
    for (i=1;i<=srptr1->size;i++)
      if (strcmp(srptr1->names[i],srptr2->names[i]))
        return FALSE;
  }
  if (srptr1->type==WORDS) {
    for (i=1;i<=srptr1->alphabet_size;i++)
      if (strcmp(srptr1->alphabet[i],srptr2->alphabet[i]))
        return FALSE;
  }
  if (srptr1->type==PRODUCT)
    if (!srec_equal(srptr1->base,srptr2->base))
      return FALSE;
  if (srptr1->type==PRODUCT) {
     if (srptr1->arity != srptr2->arity)
       return FALSE;
     if (srptr1->padding != srptr2->padding)
       return FALSE;
  }

  return TRUE;
    
}

boolean
table_equal(tableptr1,tableptr2,ne,ns)
	table_struc *tableptr1, *tableptr2;
	int ne, ns;
/* Test equality of the transition tables *tableptr1 and *tableptr2.
 * The storage-types don't need to be equal.
 */
{ int **table1, **table2, i, j, dr1, dr2;
  boolean dense1, dense2;

  dense1 = tableptr1->table_type==DENSE;
  dense2 = tableptr2->table_type==DENSE;
  table1 = tableptr1->table_data_ptr;
  table2 = tableptr2->table_data_ptr;
  dr1 = tableptr1->denserows;
  dr2 = tableptr2->denserows;

  for (j=1;j<=ne;j++) for(i=1;i<=ns;i++)
    if (target(dense1,table1,j,i,dr1) != target(dense2,table2,j,i,dr2) )
      return FALSE;

  return TRUE;
}

boolean
fsa_equal(fsaptr1,fsaptr2)
	fsa *fsaptr1, *fsaptr2;
/* Test equality of the fsa's fsaptr1 and fsaptr2 */
{ int ns, ne, ni, na, i;

  if (print_level>=3)
    printf("    #Calling fsa_equal.\n");
  if (!srec_equal(fsaptr1->alphabet,fsaptr2->alphabet) )
    return FALSE;

  if (!srec_equal(fsaptr1->states,fsaptr2->states) )
    return FALSE;

  ns = fsaptr1->states->size;
  ne = fsaptr1->alphabet->size;
 
  if (fsaptr1->num_initial != fsaptr2->num_initial)
    return FALSE;
  ni = fsaptr1->num_initial;
  if (ni>0 && ni<ns) {
    fsa_set_is_initial(fsaptr1); fsa_set_is_initial(fsaptr2);
    for (i=1;i<=ns;i++)
      if (fsaptr1->is_initial[i] != fsaptr2->is_initial[i])
        return FALSE;
    tfree(fsaptr1->is_initial);
    tfree(fsaptr2->is_initial);
  }
 
  if (fsaptr1->num_accepting != fsaptr2->num_accepting)
    return FALSE;
  na = fsaptr1->num_accepting;
  if (na>0 && na<ns) {
    fsa_set_is_accepting(fsaptr1); fsa_set_is_accepting(fsaptr2);
    for (i=1;i<=ns;i++)
      if (fsaptr1->is_accepting[i] != fsaptr2->is_accepting[i])
        return FALSE;
    tfree(fsaptr1->is_accepting);
    tfree(fsaptr2->is_accepting);
  }

  /* The flag strings are additional information, and needn't correspond */

  if (!table_equal(fsaptr1->table,fsaptr2->table,ne,ns) )
    return FALSE;

  return TRUE;
}

/* The code for the next 3 functions is so similar, that we combine them
 * into a single function fsa_binop.
 * (The negatives of these would be more difficult, because (0,0) would
 *  be an accept state - so use fsa_not for these.)
 */

fsa *
fsa_and(fsaptr1,fsaptr2,op_table_type,destroy,tempfilename)
	fsa *fsaptr1, *fsaptr2;
	storage_type op_table_type;
	boolean destroy;
	char *tempfilename;
{
  if (print_level>=3)
    printf("    #Calling fsa_and.\n");
  return 
    fsa_binop(fsaptr1,fsaptr2,op_table_type,destroy,tempfilename,AND);
}

fsa *
fsa_or(fsaptr1,fsaptr2,op_table_type,destroy,tempfilename)
	fsa *fsaptr1, *fsaptr2;
	storage_type op_table_type;
	boolean destroy;
	char *tempfilename;
{
  if (print_level>=3)
    printf("    #Calling fsa_or.\n");
  return 
    fsa_binop(fsaptr1,fsaptr2,op_table_type,destroy,tempfilename,OR);
}

fsa *
fsa_and_not(fsaptr1,fsaptr2,op_table_type,destroy,tempfilename)
	fsa *fsaptr1, *fsaptr2;
	storage_type op_table_type;
	boolean destroy;
	char *tempfilename;
{
  if (print_level>=3)
    printf("    #Calling fsa_and_not.\n");
  return 
    fsa_binop(fsaptr1,fsaptr2,op_table_type,destroy,tempfilename,AND_NOT);
}

fsa *
fsa_binop(fsaptr1,fsaptr2,op_table_type,destroy,tempfilename,op)
	fsa *fsaptr1, *fsaptr2;
	storage_type op_table_type;
	boolean destroy;
	char *tempfilename;
	binop op;
{
  int **table1, **table2, ne, ns, dr1, dr2, *fsarow,
      nt, cstate, csa, csb, im, i, g, len, ct, *ht_ptr, *tab_ptr, **table;
  boolean dense_ip1, dense_ip2, dense_op, accept;
  fsa *and_or_not;
  hash_table ht;
  FILE *tempfile, *fopen();

  if (print_level>=3)
    printf("    #Calling fsa_binop.\n");
  if (!fsaptr1->flags[DFA] || !fsaptr2->flags[DFA]){
    fprintf(stderr,"Error: fsa_binop only applies to DFA's.\n");
    exit(1);
  }

  if (!srec_equal(fsaptr1->alphabet,fsaptr2->alphabet)) {
    fprintf(stderr,"Error in fsa_binop: fsa's have different alphabets.\n");
    exit(1);
  }

  if (fsaptr1->flags[RWS])
    fsa_clear_rws(fsaptr1);

  if (fsaptr2->flags[RWS])
    fsa_clear_rws(fsaptr2);

  tmalloc(and_or_not,fsa,1);
  if (fsaptr2->num_initial==0 && (op==AND_NOT || op==OR)) {
    fsa_copy(and_or_not,fsaptr1);
    and_or_not->table->printing_format = op_table_type;
    if (destroy) {
      fsa_clear(fsaptr1); fsa_clear(fsaptr2);
    }
    return and_or_not;
  }
  if (fsaptr1->num_initial==0 && op==OR) {
    fsa_copy(and_or_not,fsaptr2);
    and_or_not->table->printing_format = op_table_type;
    if (destroy) {
      fsa_clear(fsaptr1); fsa_clear(fsaptr2);
    }
    return and_or_not;
  }

  fsa_init(and_or_not);
  srec_copy(and_or_not->alphabet,fsaptr1->alphabet);
  and_or_not->flags[DFA] = TRUE;
  and_or_not->flags[ACCESSIBLE] = TRUE;
  and_or_not->flags[BFS] = TRUE;

  and_or_not->states->type = SIMPLE;

  and_or_not->table->table_type = op_table_type;
  and_or_not->table->denserows = 0;
  and_or_not->table->printing_format = op_table_type;

  if (fsaptr1->num_initial==0 || fsaptr2->num_initial==0) {
    and_or_not->num_initial = 0;
    and_or_not->num_accepting = 0;
    and_or_not->states->size = 0;
    if (destroy) {
      fsa_clear(fsaptr1); fsa_clear(fsaptr2);
    }
    return and_or_not;
  }

  ne = fsaptr1->alphabet->size;

  fsa_set_is_accepting(fsaptr1);
  fsa_set_is_accepting(fsaptr2);
  table1 = fsaptr1->table->table_data_ptr;
  table2 = fsaptr2->table->table_data_ptr;

  dense_ip1 = fsaptr1->table->table_type==DENSE;
  dense_ip2 = fsaptr2->table->table_type==DENSE;
  dr1 = fsaptr1->table->denserows;
  dr2 = fsaptr2->table->denserows;
  dense_op = op_table_type==DENSE;

  and_or_not->num_initial = 1;
  tmalloc(and_or_not->initial,int,2);
  and_or_not->initial[1] = 1;
  
  hash_init(&ht,TRUE,2,0,0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = fsaptr1->initial[1];
  ht_ptr[1] = fsaptr2->initial[1];
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
  nt = 0; /* Number of transitions in and_or_not */

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
      ht_ptr[0] = csa ? target(dense_ip1,table1,g,csa,dr1) : 0;
      ht_ptr[1] = csb ? target(dense_ip2,table2,g,csb,dr2) : 0;
      if ( (op==AND && (ht_ptr[0]==0 || ht_ptr[1]==0)) || 
           (op==AND_NOT && ht_ptr[0]==0) )
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

  ns = and_or_not->states->size = ht.num_recs;
  and_or_not->table->numTransitions = nt;
/* Locate the accept states: first count them and then record them. */
  ct = 0;
  for (i=1;i<=ns;i++) {
    ht_ptr = hash_rec(&ht,i);
    csa = ht_ptr[0]; csb = ht_ptr[1];
    accept =
      op==AND     ? fsaptr1->is_accepting[csa] && fsaptr2->is_accepting[csb] :
      op==OR      ? fsaptr1->is_accepting[csa] || fsaptr2->is_accepting[csb] :
      op==AND_NOT ? fsaptr1->is_accepting[csa] && !fsaptr2->is_accepting[csb] :
          FALSE; /* default cannot occur */
    if (accept)
        ct++;
  }
  and_or_not->num_accepting = ct;
  if (ct==1 || ct != ns) {
    tmalloc(and_or_not->accepting,int,ct+1);
    ct = 0;
    for (i=1;i<=ns;i++) {
      ht_ptr = hash_rec(&ht,i);
      csa = ht_ptr[0]; csb = ht_ptr[1];
      accept =
      op==AND     ? fsaptr1->is_accepting[csa] && fsaptr2->is_accepting[csb] :
      op==OR      ? fsaptr1->is_accepting[csa] || fsaptr2->is_accepting[csb] :
      op==AND_NOT ? fsaptr1->is_accepting[csa] && !fsaptr2->is_accepting[csb] :
          FALSE; /* default cannot occur */
      if (accept)
        and_or_not->accepting[++ct] = i;
    }
  }
  hash_clear(&ht);
  tfree(fsarow);
  tfree(fsaptr1->is_accepting);
  tfree(fsaptr2->is_accepting);

  if (destroy) {
    fsa_clear(fsaptr1); fsa_clear(fsaptr2);
  }
/* Now read the transition table back in */
  tempfile = fopen(tempfilename,"r");
  compressed_transitions_read(and_or_not,tempfile);
  fclose(tempfile);

  unlink(tempfilename);

  return and_or_not;
}

fsa *
fsa_not(fsaptr,op_table_type)
	fsa *fsaptr;
	storage_type op_table_type;
/* This function ought to be easy - just interchange accept and non-accept
 * states. In fact it is complicated slightly by the fact that we are
 * working with partial fsa's, so the 0 state has to become a new state,
 * and a former state may become 0.
 * The storage-type of the computed fsa *not will always be dense.
 */
{ int i, j, k, ns, ne, **oldtable, **newtable, newzero, dr;
  fsa *not;
  boolean dense, zero;

  if (print_level>=3)
    printf("    #Calling fsa_not.\n");

  if (!fsaptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_not only applies to DFA's.\n");
    exit(1);
  }

  if (fsaptr->flags[RWS])
    fsa_clear_rws(fsaptr);

  tmalloc(not,fsa,1);
  fsa_init(not);
/* Most of the information is simply copied from *fsaptr to *not-
 * but we allow room for an extra accept state of *not, to replace 0
 */
  ne = fsaptr->alphabet->size;
  ns = fsaptr->states->size;

  srec_copy(not->alphabet,fsaptr->alphabet);
  not->states->type = SIMPLE;
  not->states->size = ns;
  not->num_initial = 1;
  tmalloc(not->initial,int,2);
  not->initial[1] = fsaptr->num_initial==0 ? 1 : fsaptr->initial[1];

  fsa_table_init(not->table,ns+1,ne);
  oldtable = fsaptr->table->table_data_ptr;
  newtable = not->table->table_data_ptr;
  dense = fsaptr->table->table_type==DENSE;
  dr = fsaptr->table->denserows;
  not->table->printing_format = op_table_type;
  not->table->denserows = 0;

  zero = ns==0;
  for (j=1;j<=ne;j++) for (i=1;i<=ns;i++) {
    k = target(dense,oldtable,j,i,dr);
    set_dense_target(newtable,j,i,k);
    if (!zero && k==0)
      zero = TRUE;
  }

  tmalloc(not->is_accepting,boolean,ns+2);
  if (fsaptr->num_accepting==ns) {
    for (i=1;i<=ns;i++)
      not->is_accepting[i] = FALSE;
    not->is_accepting[ns+1] = TRUE;
  }
  else {
    for (i=1;i<=ns+1;i++)
      not->is_accepting[i] = TRUE;
    for (i=1;i<=fsaptr->num_accepting;i++)
      not->is_accepting[fsaptr->accepting[i]] = FALSE;
  }

/* See if there is a new zero-state */
  newzero = 0;
  for (i=1;i<=ns;i++) if (!not->is_accepting[i]) {
    newzero = i;
    for (j=1;j<=ne;j++) if (dense_target(newtable,j,i)!=i) {
      newzero = 0;
      break;
    }
    if (newzero)
      break;
  }

  if (newzero==0 && zero)
    not->states->size = ++ns;

  not->num_accepting = 0;
  for (i=1;i<=ns;i++) if (not->is_accepting[i])
    not->num_accepting++;
  if (ns==1 || not->num_accepting!=ns) {
    k = 0;
    tmalloc(not->accepting,int,not->num_accepting+1);
    for (i=1;i<=ns;i++) if (not->is_accepting[i])
      not->accepting[++k] = i;
  }
  tfree(not->is_accepting);

  if (zero && newzero) {
  /* swap zero with the new zero-state */
    for (j=1;j<=ne;j++) for (i=1;i<=ns;i++) {
      if (i==newzero)
        continue;
      if (dense_target(newtable,j,i)==0)
        newtable[j][i] = newzero;
      else if (dense_target(newtable,j,i)==newzero)
        set_dense_target(newtable,j,i,0);
    }
  }
  else if (newzero)
    fsa_delete_state(not,newzero);
  else if (zero) {
    for (j=1;j<=ne;j++)
      set_dense_target(newtable,j,ns,ns);
    for (j=1;j<=ne;j++) for (i=1;i<ns;i++)
      if (dense_target(newtable,j,i)==0)
        set_dense_target(newtable,j,i,ns);
  }

  for (i=0;i<num_flag_strings;i++)
    not->flags[i] = fsaptr->flags[i];
  not->flags[BFS] = FALSE;
     /* BFS property would be ruined by 0 state of *fsaptr */

  return not;
}

fsa *
fsa_exists(fsaptr,op_table_type,destroy,tempfilename)
	fsa *fsaptr;
	storage_type op_table_type;
	boolean destroy;
	char *tempfilename;
/* *fsaptr must be a 2-variable fsa.
 * The returned fsa accepts a word w_1 iff (w_1,w_2) is accepted by *fsaptr,
 * for some word w_2 (with the usual padding conventions).
 * In fact, fsa_exists calls fsa_exists_int or (more usually) fsa_exists_short,
 * the latter using the short hash table functions.
 */
{
  if (print_level>=3)
    printf("    #Calling fsa_exists.\n");
  if (fsaptr->states->size < MAXUSHORT)
    return fsa_exists_short(fsaptr,op_table_type,destroy,tempfilename);
  else
    return fsa_exists_int(fsaptr,op_table_type,destroy,tempfilename);
}

fsa *
fsa_exists_short(fsaptr,op_table_type,destroy,tempfilename)
	fsa *fsaptr;
	storage_type op_table_type;
	boolean destroy;
	char *tempfilename;
{
  int **table, ne, ngens, nsi, ns, dr, *fsarow, e, es, ef,
      nt, cstate, cs, csi, im, i, g1, len, ct, *tab_ptr;
  unsigned short *ht_ptr, *ht_chptr, *ht_ptrb, *ht_ptre,
                 *cs_ptr, *cs_ptre, *ptr;
  boolean dense_ip, dense_op, got;
  short_hash_table ht;
  fsa *exists;
  FILE *tempfile, *fopen();

  if (print_level>=3)
    printf("    #Calling fsa_exists_short.\n");
  if (!fsaptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_exists only applies to DFA's.\n");
    exit(1);
  }

  if (fsaptr->alphabet->type != PRODUCT || fsaptr->alphabet->arity != 2) {
    fprintf(stderr,"Error in fsa_exists: fsa must be 2-variable.\n");
    exit(1);
  }

  tmalloc(exists,fsa,1);
  fsa_init(exists);
  srec_copy(exists->alphabet,fsaptr->alphabet->base);
  exists->flags[DFA] = TRUE;
  exists->flags[ACCESSIBLE] = TRUE;
  exists->flags[BFS] = TRUE;

  exists->states->type = SIMPLE;

  exists->table->table_type = op_table_type;
  exists->table->denserows = 0;
  exists->table->printing_format = op_table_type;

  if (fsaptr->num_initial==0) {
    exists->num_initial = 0;
    exists->num_accepting = 0;
    exists->states->size = 0;
    return exists;
  }
  ne = fsaptr->alphabet->size;
  ngens = exists->alphabet->size;
  nsi = fsaptr->states->size;

  if (ne != (ngens+1)*(ngens+1)-1) {
   fprintf(stderr,
     "Error: in a 2-variable fsa, alphabet size should = (ngens+1)^2 - 1.\n");
    exit(1);
  }

  fsa_set_is_accepting(fsaptr);

  dense_ip = fsaptr->table->table_type==DENSE;
  dr = fsaptr->table->denserows;
  dense_op = op_table_type==DENSE;
  table = fsaptr->table->table_data_ptr;

  exists->num_initial = 1;
  tmalloc(exists->initial,int,2);
  exists->initial[1] = 1;
  
  short_hash_init(&ht,FALSE,0,0,0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = fsaptr->initial[1];
  im = short_hash_locate(&ht,1);
/* Each state in 'exists' will be represented as a subset of the set of states
 * of *fsaptr. The initial state is one-element set containing the initial state
 * of *fsaptr. A subset is an accept-state if it contains an accept state of
 * *fsaptr, or if one can get to an accept-state by applying elements
 * ($,x) where $ is the padding symbol.
 * The subsets will be stored as variable-length records in the hash-table,
 * always in increasing order.
 */
  if (im!=1) {
    fprintf(stderr,"Hash-initialisation problem in fsa_exists.\n");
    exit(1);
  }
  if ((tempfile=fopen(tempfilename,"w"))==0){
    fprintf(stderr,"Error: cannot open file %s\n",tempfilename);
    exit(1);
  }
  if (dense_op)
    tmalloc(fsarow,int,ngens)
  else
    tmalloc(fsarow,int,2*ngens+1)
 
  cstate = 0;
  if (dense_op)
    len = ngens; /* The length of the fsarow output. */
  nt = 0; /* Number of transitions in exists */

  while (++cstate <= ht.num_recs) {
    if (print_level>=3) {
      if ((cstate<=1000 && cstate%100==0)||(cstate<=10000 && cstate%1000==0)||
          (cstate<=100000 && cstate%5000==0) || cstate%50000==0)
       printf("    #cstate = %d;  number of states = %d.\n",cstate,ht.num_recs);
    }
    cs_ptr = short_hash_rec(&ht,cstate);
    cs_ptre = short_hash_rec(&ht,cstate) + short_hash_rec_len(&ht,cstate) - 1;
    if (!dense_op)
      len = 0;

    for (g1=1;g1<=ngens;g1++) {
/* Calculate action of generator g1 on state cstate  - to get the image, we
 * have to apply (g1,g2) to each element in the subset corresponding to cstate,
 * and this for each generator g2 of the base-alphabet (including the padding
 * symbol).
 */
      ht_ptrb = ht.current_ptr;
      ht_ptre = ht_ptrb-1;
      ptr = cs_ptr-1;
      es = (g1-1)*(ngens+1) + 1;
      ef = g1*(ngens+1);
/* As g2 ranges from 1 to ngens+1 in the pair (g1,g2), for fixed g1, the
 * corresponding edge number in the fsa ranges from es to ef.
 */
      while (++ptr <= cs_ptre) {
        cs = *ptr;
        for (e=es;e<=ef;e++){
          csi = target(dense_ip,table,e,cs,dr);
          if (csi==0)
            continue;
          if (ht_ptrb>ht_ptre || csi> *ht_ptre) {
/* We have a new state for the image subset to be added to the end */
            *(++ht_ptre) = csi;
          }
          else {
            ht_chptr = ht_ptrb;
            while (*ht_chptr < csi)
              ht_chptr++;
            if (csi < *ht_chptr) {
/* we have a new state for the image subset to be added in the middle */
              ht_ptr = ++ht_ptre;
              while (ht_ptr > ht_chptr) {
                *ht_ptr = *(ht_ptr-1);
                ht_ptr--;
              }
              *ht_ptr = csi;
            }
          }
        }
      }
      im = short_hash_locate(&ht,ht_ptre-ht_ptrb+1);
      if (dense_op)
         fsarow[g1-1] = im;
      else if (im>0) {
         fsarow[++len] = g1;
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

  ns = exists->states->size = ht.num_recs;
  exists->table->numTransitions = nt;

/* Locate the accept states. This is slightly cumbersome, since a state
 * is an accept state if either the corresponding subset contains an
 * accept state of *fsaptr, or we can get from some such state to an
 * accept state by applying elements ($,x).
 * We will need to use the array exists->is_accepting.
 */
  tmalloc(exists->is_accepting,boolean,ns+1);
  for (i=1;i<=ns;i++)
    exists->is_accepting[i] = FALSE;
  ct = 0;
  for (cstate=ns;cstate>=1;cstate--) {
/* We work backwards through the states, since we wish to add on additional
 * elements at the end of the list in the hash-table - this destroys
 * later entries, but that doesn't matter, since we have already dealt
 * with them.
 */
    cs_ptr = short_hash_rec(&ht,cstate);
    cs_ptre = short_hash_rec(&ht,cstate) + short_hash_rec_len(&ht,cstate) - 1;
/* First see if the set itself contains an accept-state */
    ptr = cs_ptr-1;
    while (++ptr <= cs_ptre) if (fsaptr->is_accepting[*ptr]) {
      exists->is_accepting[cstate] = TRUE;
      ct++;
      break;
    }
    if (exists->is_accepting[cstate])
      continue;
/* Next apply generators ($,x) and see if we get anything new, and
 * if it is an accept state.
 * Anything new is added to the end of the list - we don't need to
 * bother about having them in increasing order this time.
 */
    es = ngens*(ngens+1)+1;
    ef = (ngens+1)*(ngens+1)-1;
    ptr = cs_ptr-1;
    while (++ptr <= cs_ptre) {
      cs = *ptr;
      for (e=es;e<=ef;e++){
        csi = target(dense_ip,table,e,cs,dr);
        if (csi==0)
          continue;
        if (fsaptr->is_accepting[csi]) {
          exists->is_accepting[cstate] = TRUE;
          ct++;
          break;
        }
/* see if csi is new */
        ht_chptr = cs_ptr-1;
        got = FALSE;
        while (++ht_chptr <= cs_ptre) if (csi == *ht_chptr) {
          got = TRUE;
          break;
        }
        if (!got) 
/* add csi to the end */
          *(++cs_ptre) = csi;
      }
      if (exists->is_accepting[cstate])
        break;
    }
  }

  exists->num_accepting = ct;
  if (ct==1 || ct != ns) {
    tmalloc(exists->accepting,int,ct+1);
    ct = 0;
    for (i=1;i<=ns;i++) if (exists->is_accepting[i])
      exists->accepting[++ct] = i;
  }
  tfree(fsaptr->is_accepting);
  tfree(exists->is_accepting);
  short_hash_clear(&ht);
  tfree(fsarow);

  if (destroy)
    fsa_clear(fsaptr);

/* Now read the transition table back in */
  tempfile = fopen(tempfilename,"r");
  compressed_transitions_read(exists,tempfile);
  fclose(tempfile);

  unlink(tempfilename);

  return exists;
}

fsa *
fsa_exists_int(fsaptr,op_table_type,destroy,tempfilename)
	fsa *fsaptr;
	table_struc op_table_type;
	boolean destroy;
	char *tempfilename;
{
  fprintf(stderr,"Sorry - fsa_exists is not yet implemented for machines.\n");
  fprintf(stderr,"with more than 65536 states.\n");
  exit(1);
}

fsa *
fsa_greater_than(alphptr)
	srec *alphptr;
/* This constructs the two-variable fsa with base-alphabet *alphptr
 * that accepts (w_1,w_2) iff w_1 > w_2 in the shortlex ordering.
 * The shorter of the two words (if any) is padded with the padding-symbol.
 * Any occurrence of the padding symbol in mid-word leads to failure.
 * The table-type of the returned fsa will always be dense.
 */
{ int ngens, ***dtable, i, j, p;
  fsa *greater_than;

  if (print_level>=3)
    printf("    #Calling fsa_greater_than.\n");
  tmalloc(greater_than,fsa,1);
  fsa_init(greater_than);

  ngens = alphptr->size;
  greater_than->alphabet->type = PRODUCT;
  greater_than->alphabet->size = (ngens+1)*(ngens+1) - 1;
  greater_than->alphabet->arity = 2;
  greater_than->alphabet->padding = '_';
  tmalloc(greater_than->alphabet->base,srec,1);
  srec_copy(greater_than->alphabet->base,alphptr);

  greater_than->states->type = SIMPLE;
  greater_than->states->size = 4;
/* state 1 = initial (fail state),
 * state 2 = words of equal length, lhs leading (accept state),
 * state 3 = words of equal length, rhs leading (fail state),
 * state 4 = eos on rhs (accept state),
 * state 0 = eos on lhs (dead state)
 */

  greater_than->num_initial = 1;
  tmalloc(greater_than->initial,int,2);
  greater_than->initial[1] = 1;

  greater_than->num_accepting = 2;
  tmalloc(greater_than->accepting,int,3);
  greater_than->accepting[1] = 2;
  greater_than->accepting[2] = 4;

  greater_than->flags[DFA] = TRUE;
  greater_than->flags[MINIMIZED] = TRUE;

  fsa_table_init(greater_than->table,4,greater_than->alphabet->size);
  greater_than->table->printing_format = DENSE;
  greater_than->table->denserows = 0;
  fsa_table_dptr_init(greater_than);
  dtable = greater_than->table->table_data_dptr;
  p = ngens+1;
  for (i=1;i<=ngens;i++) {
    set_dense_dtarget(dtable,i,p,1,4);
    set_dense_dtarget(dtable,i,p,2,4);
    set_dense_dtarget(dtable,i,p,3,4);
    set_dense_dtarget(dtable,i,p,4,4);
    set_dense_dtarget(dtable,p,i,1,0);
    set_dense_dtarget(dtable,p,i,2,0);
    set_dense_dtarget(dtable,p,i,3,0);
    set_dense_dtarget(dtable,p,i,4,0);
    set_dense_dtarget(dtable,i,i,1,1);
    set_dense_dtarget(dtable,i,i,2,2);
    set_dense_dtarget(dtable,i,i,3,3);
    set_dense_dtarget(dtable,i,i,4,0);
  }
  for (i=2;i<=ngens;i++) for(j=1;j<i;j++) {
    set_dense_dtarget(dtable,i,j,1,2);
    set_dense_dtarget(dtable,j,i,1,3);
    set_dense_dtarget(dtable,i,j,2,2);
    set_dense_dtarget(dtable,j,i,2,2);
    set_dense_dtarget(dtable,i,j,3,3);
    set_dense_dtarget(dtable,j,i,3,3);
    set_dense_dtarget(dtable,i,j,4,0);
    set_dense_dtarget(dtable,j,i,4,0);
  }

  return greater_than;
}
