/* file fsanosubword.c 7/11/94 */

/* This file essentially consists of the function fsa_no_sub_exists
 * (although this function calls fsa_nosub_exists_int or (more usually)
 * fsa_nosub_exists_short, the latter using the short hash table functions).
 * For general methodology, see comment at top of file fsalogic.c.
 *
 * This function is intended to be called only in the construction of the
 * word-acceptor of a shortlex automatic group.
 * *fsaptr must be a 2-variable finite state automaton.
 * but it will always be called in the situation where *fsaptr is
 * "WD" and "GT", where "WD" is a word-difference machine, and "GT" is
 * the "lhs >= rhs" fsa, with shortlex ordering.
 * The returned fsa accepts a word w_1 iff  w_1 contains no subword such that
 * (w_1,w_2) is accepted by *fsaptr,
 * for some word w_2 (with the usual padding conventions).
 * It is based on the function fsa_exists.
 *
 * WARNING: Since this function is called is a specific situation, it makes
 * some specific assumptions.
 * (i) *fsaptr only accepts pairs (w_1,w_2) for which w_2 is not longer
 * than w_1.
 * (ii) For each generator g, there is some transition
 * from (g,h) (for some h) with source the initial state of *fsaptr.
 *
 */

#include <stdio.h>
#include "defs.h"
#include "fsa.h"
#include "hash.h"
#include "externals.h"

#define HTMARGIN        4096
   /* if less space than this in hash-table, re-allocate */

/* Functions defined in this file: */
fsa * fsa_nosub_exists();
fsa * fsa_nosub_exists_short();
fsa * fsa_nosub_exists_int();

/* Functions used in this file and defined elsewhere */
int sparse_target();
void fsa_init();
void fsa_table_dptr_init();
void fsa_set_is_initial();
void fsa_set_is_accepting();
void srec_copy();
void fsa_clear();
void short_hash_init();
int short_hash_locate();
void short_hash_clear();
unsigned short* short_hash_rec();
int short_hash_rec_len();

fsa *
fsa_nosub_exists(fsaptr,op_table_type,destroy,tempfilename)
	fsa *fsaptr;
	storage_type op_table_type;
	boolean destroy;
	char *tempfilename;
{
  if (print_level>=3)
    printf("    #Calling fsa_nosub_exists.\n");
  if (fsaptr->states->size < MAXUSHORT)
    return fsa_nosub_exists_short(fsaptr,op_table_type,destroy,tempfilename);
  else
    return fsa_nosub_exists_int(fsaptr,op_table_type,destroy,tempfilename);
}

fsa *
fsa_nosub_exists_short(fsaptr,op_table_type,destroy,tempfilename)
	fsa *fsaptr;
	storage_type op_table_type;
	boolean destroy;
	char *tempfilename;
{
  int ***dtable, ne, ngens, nsi, ns, *fsarow,
      nt, cstate, cs, csi, im, i, g1, g2, len, ct, *tab_ptr;
  unsigned short *ht_ptr, *ht_chptr, *ht_ptrb, *ht_ptre,
                 *cs_ptr, *cs_ptre, *ptr;
  boolean dense_op, no_trans;
  short_hash_table ht;
  fsa *nosub;
  FILE *tempfile, *fopen();

  if (print_level>=3)
    printf("    #Calling fsa_nosub_exists_short.\n");
  if (!fsaptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_nosub_exists only applies to DFA's.\n");
    exit(1);
  }

  if (fsaptr->alphabet->type != PRODUCT || fsaptr->alphabet->arity != 2) {
    fprintf(stderr,"Error in fsa_nosub_exists: fsa must be 2-variable.\n");
    exit(1);
  }

  tmalloc(nosub,fsa,1);
  fsa_init(nosub);
  srec_copy(nosub->alphabet,fsaptr->alphabet->base);
  nosub->flags[DFA] = TRUE;
  nosub->flags[ACCESSIBLE] = TRUE;
  nosub->flags[BFS] = TRUE;

  ne = fsaptr->alphabet->size;
  ngens = nosub->alphabet->size;
  nsi = fsaptr->states->size;

  if (ne != (ngens+1)*(ngens+1)-1) {
   fprintf(stderr,
       "Error: in a 2-variable fsa, alphabet size should = ngens^2 - 1.\n");
    exit(1);
  }

  fsa_set_is_accepting(fsaptr);

  if (fsaptr->table->table_type!=DENSE) {
     fprintf(stderr,
       "Error: fsanosubword can only be called with a densely-stored fsa.\n");
     exit(1);
  }

  dense_op = op_table_type==DENSE;
  fsa_table_dptr_init(fsaptr);
  dtable = fsaptr->table->table_data_dptr;

  nosub->states->type = SIMPLE;
  nosub->num_initial = 1;
  tmalloc(nosub->initial,int,2);
  nosub->initial[1] = 1;
  nosub->table->table_type = op_table_type;
  nosub->table->denserows = 0;
  nosub->table->printing_format = op_table_type;
  
  short_hash_init(&ht,FALSE,0,0,0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = fsaptr->initial[1];
  im = short_hash_locate(&ht,1);
/* Each state in 'nosub' will be represented as a subset of the set of states
 * of *fsaptr. The initial state is one-element set containing the initial state
 * of *fsaptr. If a subset state contains an accept
 * state of *fsaptr, then it immediately becomes a dead state,
 * and can be discarded.
 * The subsets will be stored as variable-length records in the hash-table,
 * always in increasing order.
 */
  if (im!=1) {
    fprintf(stderr,"Hash-initialisation problem in fsa_nosub_exists.\n");
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
 * Since we are excluding words that contain subwords w_1 s.t. (w_1,w_2) is
 * accepted by *fsaptr, we also have to apply (g1,g2) to the initial state
 * of *fsaptr.
 */
      ht_ptrb = ht.current_ptr;
      ht_ptre = ht_ptrb-1;
      ptr = cs_ptr-1;
      no_trans = FALSE;
      while (ptr <= cs_ptre) {
/* We add the initial state of *fsaptr to the subset representing cstate */
        cs = ptr<cs_ptr ? fsaptr->initial[1] : *ptr;
        ptr++;
        for (g2=1;g2<=ngens+1;g2++){
          csi = dense_dtarget(dtable,g1,g2,cs);
          if (csi==0)
            continue;
          if (fsaptr->is_accepting[csi]) { /* g1 maps cstate to 0 */
            no_trans = TRUE;
            break;
          }
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
        if (no_trans)
          break;
      }
      if (no_trans) {
        if (dense_op)
          fsarow[g1-1] = 0;
        continue;
      }
      im = short_hash_locate(&ht,ht_ptre-ht_ptrb+1);
/* Because of assumption (ii) (see top of file), i can never be zero. */
      if (dense_op)
         fsarow[g1-1] = im;
      else {
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

  short_hash_clear(&ht);
  tfree(fsarow);

  ns = nosub->states->size = ht.num_recs;
  nosub->table->numTransitions = nt;

/* All states of nosub will be accept states. */
  nosub->num_accepting = ns;
  if (ns==1) {
    tmalloc(nosub->accepting,int,2);
    nosub->accepting[1] = 1;
  }
  tfree(fsaptr->is_accepting);
  if (destroy)
    fsa_clear(fsaptr);

/* Now read the transition table back in */
  tempfile = fopen(tempfilename,"r");
  compressed_transitions_read(nosub,tempfile);
  fclose(tempfile);

  unlink(tempfilename);

  return nosub;
}

fsa *
fsa_nosub_exists_int(fsaptr,op_table_type,destroy,tempfilename)
	fsa *fsaptr;
	table_struc op_table_type;
	boolean destroy;
	char *tempfilename;
{
  fprintf(stderr,"Sorry - fsa_nosub_exists is not yet implemented.\n");
  fprintf(stderr,"for machines with more than 65536 states.\n");
  exit(1);
}
