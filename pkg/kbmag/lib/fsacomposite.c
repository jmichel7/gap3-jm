/* file fsacomposite.c  24. 11. 94.
 * This file contains the routines necessary for axiom checking in an
 * automatic group.
 */

#include <stdio.h>
#include "defs.h"
#include "fsa.h"
#include "hash.h"
#include "externals.h"

/* Functions defined in this file: */
fsa  *fsa_genmult2();
fsa  *fsa_genmult2_short();
fsa  *fsa_genmult2_int();
fsa  *fsa_genmult3();
fsa  *fsa_genmult3_short();
fsa  *fsa_genmult3_int();
void  fsa_makemult();
void  fsa_makemult2();
fsa  *fsa_composite();
fsa  *fsa_composite_short();
fsa  *fsa_composite_int();

/* Functions used in this file and defined elsewhere */
boolean fsa_equal();
int sparse_target();
void fsa_init();
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
int hash_rec_len();
int short_hash_rec_len();

fsa *
fsa_genmult2(genmultptr,op_table_type,destroy,genmult2filename,readback)
	fsa *genmultptr;
	storage_type op_table_type;
	boolean destroy;
	char *genmult2filename;
	boolean readback;
/* *genmultptr should be the general multiplier fsa of an automatic group.
 * This function calculates the transition table of the general multiplier for
 * products of two generators. This is output (in unformatted form) to the
 * file with name genmult2filename, followed by a list of states, which
 * specifies which states are accept-states for which length-two words.
 * It can then be minimised for any such word as required.
 * If readback is false, then
 * the fsa returned does not have its transition table stored internally
 * (since these are in the file genmult2filename) - instead, its table
 * component table->filename contains the name of this file.
 * Otherwise, the transition table is read back in as usual.
 * If genmult2filename is the empty string, then the transition table is
 * not stored at all.
 * (This can be done when all relators of the group have length at most 4.)
 */
{
  if (print_level>=3)
    printf("    #Calling fsa_genmult2.\n");
  if (genmultptr->states->size < MAXUSHORT)
    return
      fsa_genmult2_short(genmultptr,op_table_type,
                           destroy,genmult2filename,readback);
  else
    return fsa_genmult2_int(genmultptr,op_table_type,
                           destroy,genmult2filename,readback);
}

fsa *
fsa_genmult2_short(genmultptr,op_table_type,destroy,genmult2filename,readback)
	fsa *genmultptr;
	storage_type op_table_type;
	boolean destroy;
	char *genmult2filename;
	boolean readback;
{
  int **table, ne, ngens, ngens1, ns, *fsarow, e, e1, e2, es1, ef1, dr,
      nt, cstate, im, csa, csb, csima, csimb,
      i, j, g1, g2, g3, len, reclen, nlab;
  unsigned short *ht_ptr, *ht_chptr, *ht_ptrb, *ht_ptre,
                 *cs_ptr, *cs_ptre, *ptr;
  boolean dense_ip, dense_op, got, leftpad, rightpad, keeptable;
  setToLabelsType *label, l1, l2, *newlabel;
  short_hash_table ht, labelht;
  fsa *genmult2ptr;
  srec *labelset;
  FILE *tablefile, *fopen();

  if (print_level>=3)
    printf("    #Calling fsa_genmult2_short.\n");
  if (!genmultptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_genmult2 only applies to DFA's.\n");
    exit(1);
  }

  if (genmultptr->alphabet->type!=PRODUCT || genmultptr->alphabet->arity!=2) {
    fprintf(stderr,"Error in fsa_genmult2: fsa must be 2-variable.\n");
    exit(1);
  }

  if (genmultptr->states->type!=LABELED) {
    fprintf(stderr,"Error in fsa_genmult2: states of fsa must be labeled.\n");
    exit(1);
  }

  keeptable = strlen(genmult2filename)>0;
  if (readback && !keeptable) {
   fprintf(stderr,"Error in fsa_genmult2: readback true but empty filename.\n");
    exit(1);
  }

  tmalloc(genmult2ptr,fsa,1);
  fsa_init(genmult2ptr);
  srec_copy(genmult2ptr->alphabet,genmultptr->alphabet);
  genmult2ptr->num_initial = 1;
  tmalloc(genmult2ptr->initial,int,2);
  genmult2ptr->initial[1] = 1;
  genmult2ptr->num_accepting = 0;
  /* In fact there will be lots of different categories of accept-states -
   * one for each generator pair.
   */
  genmult2ptr->flags[DFA] = TRUE;
  genmult2ptr->flags[ACCESSIBLE] = TRUE;
  genmult2ptr->flags[BFS] = TRUE;

  genmult2ptr->table->table_type = op_table_type;
  genmult2ptr->table->denserows = 0;
  genmult2ptr->table->printing_format = op_table_type;
  if (!readback) {
    tmalloc(genmult2ptr->table->filename,char,strlen(genmult2filename)+1);
    strcpy(genmult2ptr->table->filename,genmult2filename);
  }

  ne = genmultptr->alphabet->size;
  ngens = genmultptr->alphabet->base->size;
  ngens1 = ngens+1;

  if (ne != ngens1*ngens1-1) {
   fprintf(stderr,
     "Error: in a 2-variable fsa, alphabet size should = (ngens+1)^2 - 1.\n");
    exit(1);
  }

  dense_ip = genmultptr->table->table_type==DENSE;
  dr = genmultptr->table->denserows;
  dense_op = op_table_type==DENSE;
  table = genmultptr->table->table_data_ptr;

  short_hash_init(&ht,FALSE,0,0,0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = genmultptr->initial[1];
  ht_ptr[1] = genmultptr->initial[1];
  im = short_hash_locate(&ht,2);
/* Each state in the composite transition table will be represented as a subset
 * of the set of ordered pairs of states of *genmultptr.
 * The initial state contains just one such pair, of which both components are
 * the initial state of *genmultptr.
 * The subsets will be stored as variable-length records in the hash-table,
 * as a list of pairs in increasing order.
 * If a state is reached from a transition ($,x) or (x,$) (with $ the padding
 * symbol), then it needs to be marked as such, since we do not allow a $
 * to be followed by any other generator.
 * We do this by adding a 1 or a 2 to the end of the statelist - this is
 * recognisable by the fact that the length of the statelist then becomes odd.
 */
  if (im!=1) {
    fprintf(stderr,"Hash-initialisation problem in fsa_genmult2.\n");
    exit(1);
  }
  if (keeptable)
    if ((tablefile=fopen(genmult2filename,"w"))==0){
      fprintf(stderr,"Error: cannot open file %s\n",genmult2filename);
      exit(1);
    }
  if (dense_op)
    tmalloc(fsarow,int,ne)
  else
    tmalloc(fsarow,int,2*ne+1)
 
  label = genmultptr->states->setToLabels;
  cstate = 0;
  if (dense_op)
    len = ne; /* The length of the fsarow output. */
  nt = 0; /* Number of transitions in genmult2ptr */

  while (++cstate <= ht.num_recs) {
    if (print_level>=3) {
      if ((cstate<=1000 && cstate%100==0)||(cstate<=10000 && cstate%1000==0)||
          (cstate<=100000 && cstate%5000==0) || cstate%50000==0)
       printf("    #cstate = %d;  number of states = %d.\n",cstate,ht.num_recs);
    }
    reclen = short_hash_rec_len(&ht,cstate);
    cs_ptr = short_hash_rec(&ht,cstate);
    cs_ptre = short_hash_rec(&ht,cstate) + reclen - 1;
    if (reclen % 2 == 1) {
      if (*cs_ptre==1) {
        leftpad=TRUE;  rightpad=FALSE;
      } else {
        leftpad=FALSE;  rightpad=TRUE;
      }
      cs_ptre--;
    }
    else {
      leftpad=FALSE; rightpad=FALSE;
    }
      
    if (dense_op)
      for (i=0;i<ne;i++) fsarow[i] = 0;
    else
      len = 0;

    e = 0; /* e is the edge number of generator pair (g1,g3) */
    for (g1=1;g1<=ngens1;g1++)  for (g3=1;g3<=ngens1;g3++) {
/* Calculate action of pair (g1,g3) on state cstate  - to get the image, we
 * have to apply ( (g1,g2), (g2,g3) ) to each ordered pair in the subset
 * corresponding to cstate, * and this for each generator g2 of the
 * base-alphabet (including the padding symbol).
 */
       
      e++;
      if (g1==ngens1 && g3==ngens1)
        continue;

      if ((leftpad && g1<=ngens) || (rightpad && g3<=ngens))
        continue;
      ht_ptrb = ht.current_ptr;
      ht_ptre = ht_ptrb-1;

      es1 = (g1-1)*ngens1 + 1;
      ef1 = g1*ngens1;
/* As g2 ranges from 1 to ngens+1 in the pair (g1,g2), for fixed g1, the
 * corresponding edge number in the fsa ranges from es1 to ef1.
 *
 * As g2 ranges from 1 to ngens+1 in the pair (g2,g3), for fixed g3, the
 * corresponding edge number in the fsa ranges from g3 upwards, with
 * increments of size ngens+1.
 */

      ptr = cs_ptr;
      while (ptr <= cs_ptre) {
        csa = *(ptr++); csb = *(ptr++);
        for (g2=1,e1=es1,e2=g3; e1<=ef1; g2++,e1++,e2+=ngens1){
          csima = g1==ngens1 && g2==ngens1 ? (label[csa]>0 ? csa : 0) :
                   target(dense_ip,table,e1,csa,dr);
          if (csima==0)
            continue;

          csimb = g2==ngens1 && g3==ngens1 ? (label[csb]>0 ? csb : 0) :
                   target(dense_ip,table,e2,csb,dr);
          if (csimb==0)
            continue;
          if (ht_ptrb>ht_ptre || csima> *(ht_ptre-1) ||
                               (csima== *(ht_ptre-1) &&  csimb > *ht_ptre) ){
/* We have a new state for the image subset to be added to the end */
            *(++ht_ptre) = csima; *(++ht_ptre) = csimb;
          }
          else {
            ht_chptr = ht_ptrb;
            while (*ht_chptr < csima)
              ht_chptr+=2;
            if (*ht_chptr == csima)
              while (*(ht_chptr+1) < csimb)
                ht_chptr+=2;
            if (csima < *ht_chptr || csimb < *(ht_chptr+1)) {
/* we have a new state for the image subset to be added in the middle */
              ht_ptr = ht_ptre;
              ht_ptre += 2;
              while (ht_ptr >= ht_chptr) {
                *(ht_ptr+2) = *ht_ptr;
                ht_ptr--;
              }
              *ht_chptr = csima; *(ht_chptr+1) = csimb;
            }
          }
        }
      }
      if (ht_ptre > ht_ptrb) {
        if (g1==ngens1)
          *(++ht_ptre) = 1;
        else if (g3==ngens1)
          *(++ht_ptre) = 2;
      }
      im = short_hash_locate(&ht,ht_ptre-ht_ptrb+1);
      if (im > 0) {
        if (dense_op)
          fsarow[e-1] = im;
        else  {
         fsarow[++len] = e;
         fsarow[++len] = im;
        }
        nt++;
      }
    }
    if (!dense_op)
      fsarow[0] = len++;
    if (keeptable)
      fwrite((void *)fsarow,sizeof(int),(size_t)len,tablefile);
  }

  if (keeptable)
    fclose(tablefile);
  tfree(fsarow);

  ns =
  genmult2ptr->states->size = ht.num_recs;
  genmult2ptr->table->numTransitions = nt;
  genmult2ptr->num_accepting = 0;

  if (print_level>=3) {
    printf("    #Calculated transitions - %d states, %d transitions.\n",ns,nt);
    printf("    #Now calculating state labels.\n");
  }

/* Locate the accept states for  Mult_(a*b)  each generator pair (a,b).
 * This is slightly cumbersome, since a state S
 * is an accept state if either the  subset representing S contains a
 * pair (s1,s2), where s1 is accept state for Mult_a and s2 for Mult_b,
 * or we can get from to such a state by applying ( ($,g), (g,$) ) to one
 * of the pairs in S, for some generator g.
 * It is most convenient to work out this information taking the states
 * S in reverse order.
 * The information on the accept-states will be stored as labels, which
 * will be words in the generators. Each such word will have the form
 * a1*b1*a2*b2*...*ar*br, where the (ai,bi) are the generator pairs for
 * which that particular state is an accept state. The labels will be
 * collected initially in a new hash-table.
 */

  genmult2ptr->states->type = LABELED;
  tmalloc(genmult2ptr->states->labels,srec,1);
  labelset = genmult2ptr->states->labels;
  labelset->type = WORDS;
  labelset->alphabet_size = genmultptr->states->labels->size;
  for (i=1;i<=genmultptr->states->labels->size;i++) {
    tmalloc(labelset->alphabet[i],char,
                            strlen(genmultptr->states->labels->names[i])+1);
    strcpy(labelset->alphabet[i],genmultptr->states->labels->names[i]);
  }
  tmalloc(genmult2ptr->states->setToLabels,setToLabelsType,ns+1);
  newlabel = genmult2ptr->states->setToLabels;

  short_hash_init(&labelht,FALSE,0,8192,16384);

  es1 = ngens*ngens1 + 1;
  ef1 = ngens1*ngens1-1;
  for (cstate=ns;cstate>=1;cstate--) {
/* We work backwards through the states, since we wish to add on additional
 * elements at the end of the list in the hash-table - this destroys
 * later entries, but that doesn't matter, since we have already dealt
 * with them.
 */
    cs_ptr = short_hash_rec(&ht,cstate);
    reclen = short_hash_rec_len(&ht,cstate);
    if (reclen % 2 ==1)
      reclen--; /* The last entry marks the fact that this is a "padded state"*/
    cs_ptre = short_hash_rec(&ht,cstate) + reclen - 1;
/* Apply generators ( ($,g2), (g2,$) ) and see if we get anything new.
 * We won't bother about having them in increasing order this time.
 */
    ptr = cs_ptr;
    while (ptr <= cs_ptre) {
      csa = *(ptr++); csb = *(ptr++);
      for (e1=es1,e2=ngens1; e1<=ef1; e1++,e2+=ngens1){
        csima = target(dense_ip,table,e1,csa,dr);
        if (csima==0)
          continue;
        csimb = target(dense_ip,table,e2,csb,dr);
        if (csimb==0)
          continue;
  /* see if (csima,csimb) is new */
        ht_chptr = cs_ptr;
        got = FALSE;
        while (ht_chptr < cs_ptre) {
          if (csima == ht_chptr[0] && csimb == ht_chptr[1]) {
            got = TRUE;
            break;
          }
          ht_chptr+=2;
        }
        if (!got) {
  /* add (csima,csimb) to the end */
          *(++cs_ptre) = csima; *(++cs_ptre) = csimb;
        }
      }
    }
  /* Now we see which pairs in the subset are of form (s,t), where s is
   * an accept state for a generator a, and t for a generator b.
   * The list of all such pairs (a,b) will form the label of the state, which
   * will be the word a1*b1*a2*b2*...ar*br, with the (ai,bi) coming in
   * lex-order.
   */

    ht_ptrb = labelht.current_ptr;
    ht_ptre = ht_ptrb-1;
    ptr = cs_ptr;
    while (ptr <= cs_ptre) {
      csa = *(ptr++); csb = *(ptr++);
      if (((l1=label[csa]) != 0) && ((l2=label[csb]) != 0)) {
        if (ht_ptrb>ht_ptre || l1> *(ht_ptre-1) ||
                             (l1== *(ht_ptre-1) &&  l2 > *ht_ptre) ){
/* We have a new generator pair to be added to the end */
          *(++ht_ptre) = l1; *(++ht_ptre) = l2;
        }
        else {
          ht_chptr = ht_ptrb;
          while (*ht_chptr < l1)
            ht_chptr+=2;
          if (*ht_chptr == l1)
            while (*(ht_chptr+1) < l2)
              ht_chptr+=2;
          if (l1 < *ht_chptr || l2 < *(ht_chptr+1)) {
/* we have a new generator pair to be added in the middle */
            ht_ptr = ht_ptre;
            ht_ptre += 2;
            while (ht_ptr >= ht_chptr) {
              *(ht_ptr+2) = *ht_ptr;
              ht_ptr--;
            }
            *ht_chptr = l1; *(ht_chptr+1) = l2;
          }
        }
      }
    }
/* That completes the calculation of the label for cstate */
    newlabel[cstate] = short_hash_locate(&labelht,ht_ptre-ht_ptrb+1);
  } 
  short_hash_clear(&ht);

/* Finally copy the records from the lable hash-table into the set of labels */
  nlab =
  labelset->size = labelht.num_recs;
  if (print_level>=3)
    printf("    #There are %d distinct labels.\n",nlab);
  tmalloc(labelset->names,char *,nlab+1);
  for (i=1;i<=nlab;i++) {
    len = short_hash_rec_len(&labelht,i);
    tmalloc(labelset->names[i],char,len+1);
    ht_ptr = short_hash_rec(&labelht,i);
    for (j=0;j<len;j++)
      labelset->names[i][j] = ht_ptr[j];
    labelset->names[i][len] = '\0';
  }

  short_hash_clear(&labelht);
  if (destroy)
    fsa_clear(genmultptr);

  if (readback) {
    tablefile = fopen(genmult2filename,"r");
    compressed_transitions_read(genmult2ptr,tablefile);
    fclose(tablefile);
    unlink(genmult2filename);
  }

  return genmult2ptr;
}

fsa *
fsa_genmult2_int(genmultptr,op_table_type,destroy,genmult2filename,readback)
	fsa *genmultptr;
	storage_type op_table_type;
	boolean destroy;
	char *genmult2filename;
	boolean readback;
{
  fprintf(stderr,"Sorry - fsa_genmult2 is not yet implemented for machines.\n");
  fprintf(stderr,"with more than 65536 states.\n");
  exit(1);
}

fsa *
fsa_genmult3(genmultptr,genmult2ptr,op_table_type,destroy,
                                                 genmult3filename,readback)
	fsa *genmultptr, *genmult2ptr;
	storage_type op_table_type;
	boolean destroy;
	char *genmult3filename;
	boolean readback;
/* *genmultptr should be the general multiplier fsa of an automatic group.
 * and *genmult2ptr the length 2-multiplier.
 * This function calculates the transition table of the general multiplier for
 * products of three generators. This is output (in unformatted form) to the
 * file with name genmult3filename, followed by a list of states, which
 * specifies which states are accept-states for which length-two words.
 * It can then be minimised for any such word as required.
 * If readback is false, then
 * the fsa returned does not have its transition table stored internally
 * (since these are in the file genmult3filename) - instead, its table
 * component table->filename contains the name of this file.
 * Otherwise, the transition table is read back in as usual.
 */
{
  if (print_level>=3)
    printf("    #Calling fsa_genmult3.\n");
  if (genmultptr->states->size < MAXUSHORT)
    return fsa_genmult3_short(genmultptr,genmult2ptr,op_table_type,destroy,
                                             genmult3filename,readback);
  else
    return fsa_genmult3_int(genmultptr,genmult2ptr,op_table_type,destroy,
                                             genmult3filename,readback);
}

fsa *
fsa_genmult3_short(genmultptr,genmult2ptr,op_table_type,destroy,
                                             genmult3filename,readback)
	fsa *genmultptr, *genmult2ptr;
	storage_type op_table_type;
	boolean destroy;
	char *genmult3filename;
	boolean readback;
{
  int **table1, **table2, ne, ngens, ngens1, ns, *fsarow,
      e, e1, e2, es1, ef1, dr1, dr2, nt, cstate, im, csa, csb, csima, csimb,
      i, j, k, g1, g2, g3, len, reclen, nlab, lablen;
  unsigned short *ht_ptr, *ht_chptr, *ht_ptrb, *ht_ptre,
                 *cs_ptr, *cs_ptre, *ptr;
  boolean dense_ip1, dense_ip2, dense_op, got, leftpad, rightpad;
  setToLabelsType *label1, *label2, l1, l2, *newlabel;
  char **labelword, *lni;
  short_hash_table ht, labelht;
  fsa *genmult3ptr;
  srec *labelset;
  FILE *tablefile, *fopen();

  if (print_level>=3)
    printf("    #Calling fsa_genmult3_short.\n");
  if (!genmultptr->flags[DFA] || !genmult2ptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_genmult3 only applies to DFA's.\n");
    exit(1);
  }

  if (genmultptr->alphabet->type!=PRODUCT || genmultptr->alphabet->arity!=2) {
    fprintf(stderr,"Error in fsa_genmult3: fsa must be 2-variable.\n");
    exit(1);
  }
  if (!srec_equal(genmultptr->alphabet,genmult2ptr->alphabet) ) {
    fprintf(stderr,"Error in fsa_genmult3: fsa's must have same alphabets.\n");
    exit(1);
  }

  if (genmultptr->states->type!=LABELED || genmult2ptr->states->type!=LABELED) {
    fprintf(stderr,"Error in fsa_genmult3: states of fsa must be labeled.\n");
    exit(1);
  }
  if (genmultptr->states->labels->type != IDENTIFIERS) {
    fprintf(stderr,
        "Error in fsa_genmult3: labels of first fsa must be identifiers.\n");
    exit(1);
  }
  if (genmult2ptr->states->labels->type != WORDS) {
    fprintf(stderr,
        "Error in fsa_genmult3: labels of first fsa must be words.\n");
    exit(1);
  }

  tmalloc(genmult3ptr,fsa,1);
  fsa_init(genmult3ptr);
  srec_copy(genmult3ptr->alphabet,genmultptr->alphabet);
  genmult3ptr->num_initial = 1;
  tmalloc(genmult3ptr->initial,int,2);
  genmult3ptr->initial[1] = 1;
  genmult3ptr->num_accepting = 0;
  /* In fact there will be lots of different categories of accept-states -
   * one for each generator pair.
   */
  genmult3ptr->flags[DFA] = TRUE;
  genmult3ptr->flags[ACCESSIBLE] = TRUE;
  genmult3ptr->flags[BFS] = TRUE;

  genmult3ptr->table->table_type = op_table_type;
  genmult3ptr->table->denserows = 0;
  genmult3ptr->table->printing_format = op_table_type;
  if (!readback) {
    tmalloc(genmult3ptr->table->filename,char,strlen(genmult3filename)+1);
    strcpy(genmult3ptr->table->filename,genmult3filename);
  }

  ne = genmultptr->alphabet->size;
  ngens = genmultptr->alphabet->base->size;
  ngens1 = ngens+1;

  if (ne != ngens1*ngens1-1) {
   fprintf(stderr,
     "Error: in a 2-variable fsa, alphabet size should = (ngens+1)^2 - 1.\n");
    exit(1);
  }

  dense_ip1 = genmultptr->table->table_type==DENSE;
  dr1 = genmultptr->table->denserows;
  dense_ip2 = genmult2ptr->table->table_type==DENSE;
  dr2 = genmult2ptr->table->denserows;
  dense_op = op_table_type==DENSE;
  table1 = genmultptr->table->table_data_ptr;
  table2 = genmult2ptr->table->table_data_ptr;

  short_hash_init(&ht,FALSE,0,0,0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = genmultptr->initial[1];
  ht_ptr[1] = genmult2ptr->initial[1];
  im = short_hash_locate(&ht,2);
/* Each state in the composite transition table will be represented as a subset
 * of the set of ordered pairs of states the first from *genmultptr
 * and the second from *genmult2ptr.
 * The initial state contains just one such pair, of which the components are
 * the initial states of *genmultptr and *genmult2ptr.
 * The subsets will be stored as variable-length records in the hash-table,
 * as a list of pairs in increasing order.
 * If a state is reached from a transition ($,x) or (x,$) (with $ the padding
 * symbol), then it needs to be marked as such, since we do not allow a $
 * to be followed by any other generator.
 * We do this by adding a 1 or a 2 to the end of the statelist - this is
 * recognisable by the fact that the length of the statelist then becomes odd.
 */
  if (im!=1) {
    fprintf(stderr,"Hash-initialisation problem in fsa_genmult3.\n");
    exit(1);
  }
  if ((tablefile=fopen(genmult3filename,"w"))==0){
    fprintf(stderr,"Error: cannot open file %s\n",genmult3filename);
    exit(1);
  }
  if (dense_op)
    tmalloc(fsarow,int,ne)
  else
    tmalloc(fsarow,int,2*ne+1)
 
  label1 = genmultptr->states->setToLabels;
  label2 = genmult2ptr->states->setToLabels;
  cstate = 0;
  if (dense_op)
    len = ne; /* The length of the fsarow output. */
  nt = 0; /* Number of transitions in genmult3ptr */

  while (++cstate <= ht.num_recs) {
    if (print_level>=3) {
      if ((cstate<=1000 && cstate%100==0)||(cstate<=10000 && cstate%1000==0)||
          (cstate<=100000 && cstate%5000==0) || cstate%50000==0)
       printf("    #cstate = %d;  number of states = %d.\n",cstate,ht.num_recs);
    }
    reclen = short_hash_rec_len(&ht,cstate);
    cs_ptr = short_hash_rec(&ht,cstate);
    cs_ptre = short_hash_rec(&ht,cstate) + reclen - 1;
    if (reclen % 2 == 1) {
      if (*cs_ptre==1) {
        leftpad=TRUE;  rightpad=FALSE;
      } else {
        leftpad=FALSE;  rightpad=TRUE;
      }
      cs_ptre--;
    }
    else {
      leftpad=FALSE; rightpad=FALSE;
    }
      
    if (dense_op)
      for (i=0;i<ne;i++) fsarow[i] = 0;
    else
      len = 0;
    e = 0; /* e is the edge number of generator pair (g1,g3) */
    for (g1=1;g1<=ngens1;g1++)  for (g3=1;g3<=ngens1;g3++) {
/* Calculate action of pair (g1,g3) on state cstate  - to get the image, we
 * have to apply ( (g1,g2), (g2,g3) ) to each ordered pair in the subset
 * corresponding to cstate, * and this for each generator g2 of the
 * base-alphabet (including the padding symbol).
 */
      e++;
      if (g1==ngens1 && g3==ngens1)
        continue;
      if ((leftpad && g1<=ngens) || (rightpad && g3<=ngens))
        continue;
      ht_ptrb = ht.current_ptr;
      ht_ptre = ht_ptrb-1;

      es1 = (g1-1)*ngens1 + 1;
      ef1 = g1*ngens1;
/* As g2 ranges from 1 to ngens+1 in the pair (g1,g2), for fixed g1, the
 * corresponding edge number in the fsa ranges from es1 to ef1.
 *
 * As g2 ranges from 1 to ngens+1 in the pair (g2,g3), for fixed g3, the
 * corresponding edge number in the fsa ranges from g3 upwards, with
 * increments of size ngens+1.
 */

      ptr = cs_ptr;
      while (ptr <= cs_ptre) {
        csa = *(ptr++); csb = *(ptr++);
        for (g2=1,e1=es1,e2=g3; e1<=ef1; g2++,e1++,e2+=ngens1){
          csima = g1==ngens1 && g2==ngens1 ? (label1[csa]>0 ? csa : 0) :
                   target(dense_ip1,table1,e1,csa,dr1);
          if (csima==0)
            continue;

          csimb = g2==ngens1 && g3==ngens1 ? (label2[csb]>0 ? csb : 0) :
                   target(dense_ip2,table2,e2,csb,dr2);
          if (csimb==0)
            continue;

          if (ht_ptrb>ht_ptre || csima> *(ht_ptre-1) ||
                               (csima== *(ht_ptre-1) &&  csimb > *ht_ptre) ){
/* We have a new state for the image subset to be added to the end */
            *(++ht_ptre) = csima; *(++ht_ptre) = csimb;
          }
          else {
            ht_chptr = ht_ptrb;
            while (*ht_chptr < csima)
              ht_chptr+=2;
            if (*ht_chptr == csima)
              while (*(ht_chptr+1) < csimb)
                ht_chptr+=2;
            if (csima < *ht_chptr || csimb < *(ht_chptr+1)) {
/* we have a new state for the image subset to be added in the middle */
              ht_ptr = ht_ptre;
              ht_ptre += 2;
              while (ht_ptr >= ht_chptr) {
                *(ht_ptr+2) = *ht_ptr;
                ht_ptr--;
              }
              *ht_chptr = csima; *(ht_chptr+1) = csimb;
            }
          }
        }
      }
      if (ht_ptre > ht_ptrb) {
        if (g1==ngens1)
          *(++ht_ptre) = 1;
        else if (g3==ngens1)
          *(++ht_ptre) = 2;
      }
      im = short_hash_locate(&ht,ht_ptre-ht_ptrb+1);
      if (im > 0) {
        if (dense_op)
          fsarow[e-1] = im;
        else  {
         fsarow[++len] = e;
         fsarow[++len] = im;
        }
        nt++;
      }
    }
    if (!dense_op)
      fsarow[0] = len++;
    fwrite((void *)fsarow,sizeof(int),(size_t)len,tablefile);
  }

  fclose(tablefile);
  tfree(fsarow);

  ns =
  genmult3ptr->states->size = ht.num_recs;
  genmult3ptr->table->numTransitions = nt;
  genmult3ptr->num_accepting = 0;

  if (print_level>=3) {
    printf("    #Calculated transitions - %d states, %d transitions.\n",ns,nt);
    printf("    #Now calculating state labels.\n");
  }

/* Locate the accept states for  Mult_(a*b*c) for each generator triple (a,b,c).
 * This is slightly cumbersome, since a state S
 * is an accept state if either the subset representing S contains a
 * pair (s1,s2), where s1 is accept state for Mult_a and s2 for Mult_b*c,
 * or we can get from to such a state by applying ( ($,g), (g,$) ) to one
 * of the pairs in S, for some generator g.
 * It is most convenient to work out this information taking the states
 * S in reverse order.
 * The information on the accept-states will be stored as labels, which
 * will be words in the generators. Each such word will have the form
 * a1*b1*c1*a2*b2*c2*...*ar*br*cr, where the (ai,bi,ci) are the generator pairs
 * for which that particular state is an accept state. The labels will be
 * collected initially in a new hash-table.
 */

  genmult3ptr->states->type = LABELED;
  tmalloc(genmult3ptr->states->labels,srec,1);
  labelset = genmult3ptr->states->labels;
  labelset->type = WORDS;
  labelset->alphabet_size = genmultptr->states->labels->size;
  for (i=1;i<=genmultptr->states->labels->size;i++) {
    tmalloc(labelset->alphabet[i],char,
                            strlen(genmultptr->states->labels->names[i])+1);
    strcpy(labelset->alphabet[i],genmultptr->states->labels->names[i]);
  }
  tmalloc(genmult3ptr->states->setToLabels,setToLabelsType,ns+1);
  newlabel = genmult3ptr->states->setToLabels;

  short_hash_init(&labelht,FALSE,0,8192,16384);

  es1 = ngens*ngens1 + 1;
  ef1 = ngens1*ngens1-1;
  for (cstate=ns;cstate>=1;cstate--) {
/* We work backwards through the states, since we wish to add on additional
 * elements at the end of the list in the hash-table - this destroys
 * later entries, but that doesn't matter, since we have already dealt
 * with them.
 */
    cs_ptr = short_hash_rec(&ht,cstate);
    reclen = short_hash_rec_len(&ht,cstate);
    if (reclen % 2 ==1)
      reclen--; /* The last entry marks the fact that this is a "padded state"*/
    cs_ptre = short_hash_rec(&ht,cstate) + reclen - 1;
/* Apply generators ( ($,g2), (g2,$) ) and see if we get anything new.
 * We won't bother about having them in increasing order this time.
 */
    ptr = cs_ptr;
    while (ptr <= cs_ptre) {
      csa = *(ptr++); csb = *(ptr++);
      for (e1=es1,e2=ngens1; e1<=ef1; e1++,e2+=ngens1){
        csima = target(dense_ip1,table1,e1,csa,dr1);
        if (csima==0)
          continue;
        csimb = target(dense_ip2,table2,e2,csb,dr2);
        if (csimb==0)
          continue;

  /* see if (csima,csimb) is new */
        ht_chptr = cs_ptr;
        got = FALSE;
        while (ht_chptr < cs_ptre) {
          if (csima == ht_chptr[0] && csimb == ht_chptr[1]) {
            got = TRUE;
            break;
          }
          ht_chptr+=2;
        }
        if (!got) {
  /* add (csima,csimb) to the end */
          *(++cs_ptre) = csima; *(++cs_ptre) = csimb;
        }
      }
    }
  /* Now we see which pairs in the subset are of form (s,t), where s is
   * an accept state for a generator a, and t for a generator pair (b,c).
   * The list of all such pairs (a,b,c) will form the label of the state, which
   * will be the word a1*b1*c1*a2*b2*c2...ar*br*cr, with the (ai,bi,ci) coming
   * in the lex-order of (ai,labelnumber(bi,ci)) (which is the easiest).
   */

    ht_ptrb = labelht.current_ptr;
    ht_ptre = ht_ptrb-1;
    ptr = cs_ptr;
    while (ptr <= cs_ptre) {
      csa = *(ptr++); csb = *(ptr++);
      if (((l1=label1[csa]) != 0) && ((l2=label2[csb]) != 0)) {
        if (ht_ptrb>ht_ptre || l1> *(ht_ptre-1) ||
                             (l1== *(ht_ptre-1) &&  l2 > *ht_ptre) ){
/* We have a new label pair to be added to the end */
          *(++ht_ptre) = l1; *(++ht_ptre) = l2;
        }
        else {
          ht_chptr = ht_ptrb;
          while (*ht_chptr < l1)
            ht_chptr+=2;
          if (*ht_chptr == l1)
            while (*(ht_chptr+1) < l2)
              ht_chptr+=2;
          if (l1 < *ht_chptr || l2 < *(ht_chptr+1)) {
/* we have a new label pair to be added in the middle */
            ht_ptr = ht_ptre;
            ht_ptre += 2;
            while (ht_ptr >= ht_chptr) {
              *(ht_ptr+2) = *ht_ptr;
              ht_ptr--;
            }
            *ht_chptr = l1; *(ht_chptr+1) = l2;
          }
        }
      }
    }
/* That completes the calculation of the label for cstate */
    newlabel[cstate] = short_hash_locate(&labelht,ht_ptre-ht_ptrb+1);
  } 
  short_hash_clear(&ht);

/* Finally copy the records from the lable hash-table into the set of labels */
  nlab =
  labelset->size = labelht.num_recs;
  if (print_level>=3)
    printf("    #There are %d distinct labels.\n",nlab);
  tmalloc(labelset->names,char *,nlab+1);
  labelword = genmult2ptr->states->labels->names;
  for (i=1;i<=nlab;i++) {
    len = short_hash_rec_len(&labelht,i);
    ht_ptr = short_hash_rec(&labelht,i);
    lablen = 0;
    for (j=1;j<len;j+=2)
      lablen += (3*strlen(labelword[ht_ptr[j]])/2);
    tmalloc(labelset->names[i],char,lablen+1);
    lni = labelset->names[i];
    lni[0] = '\0';
    for (j=1;j<len;j+=2)  for (k=0;k<strlen(labelword[ht_ptr[j]]);k+=2)
      sprintf(lni+strlen(lni),"%c%c%c",
            ht_ptr[j-1], labelword[ht_ptr[j]][k], labelword[ht_ptr[j]][k+1]);
  }

  short_hash_clear(&labelht);
  if (destroy) {
    fsa_clear(genmultptr); fsa_clear(genmult2ptr);
  }

  if (readback) {
    tablefile = fopen(genmult3filename,"r");
    compressed_transitions_read(genmult3ptr,tablefile);
    fclose(tablefile);
    unlink(genmult3filename);
  }

  return genmult3ptr;
}

fsa *
fsa_genmult3_int(genmultptr,genmult2ptr,op_table_type,
                                          destroy,genmult3filename,readback)
	fsa *genmultptr, *genmult2ptr;
	storage_type op_table_type;
	boolean destroy;
	char *genmult3filename;
	boolean readback;
{
  fprintf(stderr,"Sorry - fsa_genmult3 is not yet implemented for machines.\n");
  fprintf(stderr,"with more than 65536 states.\n");
  exit(1);
}

void
fsa_makemult(genmultptr,g)
        fsa *genmultptr;
        int g;
/* This procedure takes the fsa *genmultptr produced by fsa_triples, 
 * and builds a particular  multiplier Mult_g1.
 * This merely involves setting the accept states of *genmultptr
 * in accordance with the labels of its states.
 * g can be ngens+1, for the identity generator, which inevitably
 * produces the diagonal of the word-acceptor.
 * This procedure alters its arguments and does not return anything.
 */
{
  int ne, ngens, ns, i, ct;
  setToLabelsType *label;

  if (print_level>=3)
    printf("    #Calling fsa_makemult with generator number %d.\n",g);
  if (!genmultptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_makemult only applies to DFA's.\n");
    exit(1);
  }

  if (genmultptr->alphabet->type!=PRODUCT || genmultptr->alphabet->arity!=2) {
    fprintf(stderr,"Error in fsa_makemult: fsa must be 2-variable.\n");
    exit(1);
  }

  if (genmultptr->states->type!=LABELED) {
    fprintf(stderr,"Error in fsa_genmult2: states of fsa must be labeled.\n");
    exit(1);
  }

  ns = genmultptr->states->size;
  ngens = genmultptr->alphabet->base->size;
  if (g<0 || g>ngens+1) {
    fprintf(stderr,"#Error in fsa_makemult: Generator is out of range.\n");
    exit(1);
  }
  
  tfree(genmultptr->accepting);
  tfree(genmultptr->is_accepting);
  label = genmultptr->states->setToLabels;
  ct = 0;
  for (i=1;i<=ns;i++) if (label[i]==g) 
    ct++;
  genmultptr->num_accepting = ct;
  if (ct<ns || ns==1) { 
    tfree(genmultptr->accepting);
    tmalloc(genmultptr->accepting,int,ct+1);
    ct = 0;
    for (i=1;i<=ns;i++) if (label[i]==g)
      genmultptr->accepting[++ct] = i;
  }
/* The state labelling is no longer relevant, so we clear it */
  genmultptr->states->type=SIMPLE;
  srec_clear(genmultptr->states->labels);
  tfree(genmultptr->states->labels);
  tfree(genmultptr->states->setToLabels);
}

void
fsa_makemult2(genmult2ptr,g1,g2)
        fsa *genmult2ptr;
        int g1, g2;
/* This procedure takes the fsa *genmult2ptr produced by fsa_genmult2, 
 * and builds a particular length-2 multiplier Mult_g1g2.
 * This merely involves locating the accept states.
 * This procedure alters its arguments and does not return anything.
 */
{
  int ne, ngens, ns, nlabs, *staterow, i, j, ct, len;
  char **labelname;
  boolean *accept;
  setToLabelsType *labelnumber;

  if (print_level>=3)
    printf("    #Calling fsa_makemult2 with generators %d and %d.\n",g1,g2);
  if (!genmult2ptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_makemult2 only applies to DFA's.\n");
    exit(1);
  }

  if (genmult2ptr->alphabet->type!=PRODUCT || genmult2ptr->alphabet->arity!=2) {
    fprintf(stderr,"Error in fsa_makemult2: fsa must be 2-variable.\n");
    exit(1);
  }

  if (genmult2ptr->states->type!=LABELED) {
    fprintf(stderr,"Error in fsa_genmult2: states of fsa must be labeled.\n");
    exit(1);
  }

  if (genmult2ptr->states->labels->type!=WORDS) {
    fprintf(stderr,"Error in fsa_genmult2: labels must be words.\n");
    exit(1);
  }

  ns = genmult2ptr->states->size;
  nlabs = genmult2ptr->states->labels->size;
  ngens = genmult2ptr->states->labels->alphabet_size;
  if (g1<0 || g1>ngens || g2<0 || g2>ngens) {
    fprintf(stderr,"#Error in fsa_makemult2: A generator is out of range.\n");
    exit(1);
  }
  
  tmalloc(accept,boolean,nlabs+1);
  labelname = genmult2ptr->states->labels->names;
/* Now we see which labels are accept-labels for the pair (g1,g2) - the
 * label is an accept-label if the word which is its name conatins g1*g2.
 */
  accept[0] = FALSE;
  for (i=1;i<=nlabs;i++) {
    accept[i] = FALSE;
    len = strlen(labelname[i]);
    for (j=0;j<len;j+=2) if (labelname[i][j]==g1 && labelname[i][j+1]==g2) {
      accept[i] = TRUE;
      break;
    }
  }

/* Now we can see which states are accept-states. A state is an accept-state
 * iff its label is an accept-label.
 * First we count the number.
 */
  ct = 0;
  labelnumber = genmult2ptr->states->setToLabels;
  for (i=1;i<=ns;i++) if (accept[labelnumber[i]])
    ct++;

  genmult2ptr->num_accepting = ct;
  if (ct<ns || ns==1) { 
    tfree(genmult2ptr->accepting);
    tmalloc(genmult2ptr->accepting,int,ct+1);
    ct = 0;
    for (i=1;i<=ns;i++) if (accept[labelnumber[i]])
      genmult2ptr->accepting[++ct] = i;
  }

  tfree(accept);

/* The state labelling is no longer relevant, so we clear it */
  genmult2ptr->states->type=SIMPLE;
  srec_clear(genmult2ptr->states->labels);
  tfree(genmult2ptr->states->labels);
  tfree(genmult2ptr->states->setToLabels);
}

fsa *
fsa_composite(mult1ptr,mult2ptr,op_table_type,destroy,compfilename,readback)
	fsa *mult1ptr, *mult2ptr;
	storage_type op_table_type;
	boolean destroy;
	char *compfilename;
        boolean readback;
/* *mult1ptr and *mult2ptr should be multiplier fsa's of an automatic group.
 * This function calculates the composite of these two multipliers.
 * So if *mult1ptr is the mutiplier for the word w1 and *mult2ptr for w2,
 * then the retruned fsa is the multiplier for w1*w2.
 */
{
  if (print_level>=3)
    printf("    #Calling fsa_composite.\n");
  if (mult1ptr->states->size < MAXUSHORT && mult2ptr->states->size < MAXUSHORT)
    return fsa_composite_short(mult1ptr,mult2ptr,op_table_type,
                                                destroy,compfilename,readback);
  else
    return fsa_composite_int(mult1ptr,mult2ptr,op_table_type,
                                                destroy,compfilename,readback);
}

fsa *
fsa_composite_short(mult1ptr,mult2ptr,op_table_type,
                                                  destroy,compfilename,readback)
	fsa *mult1ptr, *mult2ptr;
	storage_type op_table_type;
	boolean destroy;
	char *compfilename;
        boolean readback;
{
  int **table1, **table2, ne, ngens, ngens1, ns, *fsarow,
      e, e1, e2, es1, ef1, dr1, dr2, nt, cstate, im, csa, csb, csima, csimb,
      i, ct, g1, g2, g3, len, reclen;
  unsigned short *ht_ptr, *ht_chptr, *ht_ptrb, *ht_ptre,
                 *cs_ptr, *cs_ptre, *ptr;
  boolean dense_ip1, dense_ip2, dense_op, got, leftpad, rightpad;
  short_hash_table ht;
  fsa *compositeptr;
  FILE *tempfile, *fopen();

  if (print_level>=3)
    printf("    #Calling fsa_composite_short.\n");
  if (!mult1ptr->flags[DFA] || !mult2ptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_composite only applies to DFA's.\n");
    exit(1);
  }

  if (mult1ptr->alphabet->type!=PRODUCT || mult1ptr->alphabet->arity!=2) {
    fprintf(stderr,"Error in fsa_composite: fsa must be 2-variable.\n");
    exit(1);
  }
  if (mult2ptr->alphabet->type!=PRODUCT || mult2ptr->alphabet->arity!=2) {
    fprintf(stderr,"Error in fsa_composite: fsa must be 2-variable.\n");
    exit(1);
  }

  if (!srec_equal(mult1ptr->alphabet,mult2ptr->alphabet)) {
    fprintf(stderr,"Error in fsa_composite: fsa's must have same alphabet.\n");
    exit(1);
  }

  tmalloc(compositeptr,fsa,1);
  fsa_init(compositeptr);
  srec_copy(compositeptr->alphabet,mult1ptr->alphabet);
  compositeptr->states->type = SIMPLE;
  compositeptr->num_initial = 1;
  tmalloc(compositeptr->initial,int,2);
  compositeptr->initial[1] = 1;
  compositeptr->flags[DFA] = TRUE;
  compositeptr->flags[ACCESSIBLE] = TRUE;
  compositeptr->flags[BFS] = TRUE;

  compositeptr->table->table_type = op_table_type;
  compositeptr->table->denserows = 0;
  compositeptr->table->printing_format = op_table_type;
  if (!readback) {
    tmalloc(compositeptr->table->filename,char,strlen(compfilename)+1);
    strcpy(compositeptr->table->filename,compfilename);
  }

  ne = mult1ptr->alphabet->size;
  ngens = mult1ptr->alphabet->base->size;
  ngens1 = ngens+1;

  if (ne != ngens1*ngens1-1) {
   fprintf(stderr,
     "Error: in a 2-variable fsa, alphabet size should = (ngens+1)^2 - 1.\n");
    exit(1);
  }

  fsa_set_is_accepting(mult1ptr);
  fsa_set_is_accepting(mult2ptr);

  dense_ip1 = mult1ptr->table->table_type==DENSE;
  dr1 = mult1ptr->table->denserows;
  dense_ip2 = mult2ptr->table->table_type==DENSE;
  dr2 = mult2ptr->table->denserows;
  dense_op = op_table_type==DENSE;
  table1 = mult1ptr->table->table_data_ptr;
  table2 = mult2ptr->table->table_data_ptr;

  short_hash_init(&ht,FALSE,0,0,0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = mult1ptr->initial[1];
  ht_ptr[1] = mult2ptr->initial[1];
  im = short_hash_locate(&ht,2);
/* Each state in the composite transition table will be represented as a subset
 * of the set of ordered pairs in which the first component is a state in
 * *multptr1 and the second a state in *multptr2.
 * The initial state contains just one such pair, of which the components are
 * the initial states of *mult1ptr and *multptr2.
 * The subsets will be stored as variable-length records in the hash-table,
 * as a list of pairs in increasing order.
 * If a state is reached from a transition ($,x) or (x,$) (with $ the padding
 * symbol), then it needs to be marked as such, since we do not allow a $
 * to be followed by any other generator.
 * We do this by adding a 1 or a 2 to the end of the statelist - this is
 * recognisable by the fact that the length of the statelist then becomes odd.
 */
  if (im!=1) {
    fprintf(stderr,"Hash-initialisation problem in fsa_composite.\n");
    exit(1);
  }
  if ((tempfile=fopen(compfilename,"w"))==0){
    fprintf(stderr,"Error: cannot open file %s\n",compfilename);
    exit(1);
  }
  if (dense_op)
    tmalloc(fsarow,int,ne)
  else
    tmalloc(fsarow,int,2*ne+1)
 
  cstate = 0;
  if (dense_op)
    len = ne; /* The length of the fsarow output. */
  nt = 0; /* Number of transitions in compositeptr */

  while (++cstate <= ht.num_recs) {
    if (print_level>=3) {
      if ((cstate<=1000 && cstate%100==0)||(cstate<=10000 && cstate%1000==0)||
          (cstate<=100000 && cstate%5000==0) || cstate%50000==0)
       printf("    #cstate = %d;  number of states = %d.\n",cstate,ht.num_recs);
    }
    reclen = short_hash_rec_len(&ht,cstate);
    cs_ptr = short_hash_rec(&ht,cstate);
    cs_ptre = short_hash_rec(&ht,cstate) + reclen - 1;
    /* for (i=0;i<reclen;i++) printf(" %d",cs_ptr[i]); printf("\n"); */
    if (reclen % 2 == 1) {
      if (*cs_ptre==1) {
        leftpad=TRUE;  rightpad=FALSE;
      } else {
        leftpad=FALSE;  rightpad=TRUE;
      }
      cs_ptre--;
    }
    else {
      leftpad=FALSE; rightpad=FALSE;
    }
      
    if (dense_op)
      for (i=0;i<ne;i++) fsarow[i] = 0;
    else
      len = 0;

    e = 0; /* e is the edge number of generator pair (g1,g3) */
    for (g1=1;g1<=ngens1;g1++)  for (g3=1;g3<=ngens1;g3++) {
/* Calculate action of pair (g1,g3) on state cstate  - to get the image, we
 * have to apply ( (g1,g2), (g2,g3) ) to each ordered pair in the subset
 * corresponding to cstate, * and this for each generator g2 of the
 * base-alphabet (including the padding symbol).
 */
      e++;
      if (g1==ngens1 && g3==ngens1)
        continue;
      if ((leftpad && g1<=ngens) || (rightpad && g3<=ngens))
        continue;
      ht_ptrb = ht.current_ptr;
      ht_ptre = ht_ptrb-1;
      es1 = (g1-1)*ngens1 + 1;
      ef1 = g1*ngens1;
/* As g2 ranges from 1 to ngens+1 in the pair (g1,g2), for fixed g1, the
 * corresponding edge number in the fsa ranges from es1 to ef1.
 *
 * As g2 ranges from 1 to ngens+1 in the pair (g2,g3), for fixed g3, the
 * corresponding edge number in the fsa ranges from g3 upwards, with
 * increments of size ngens+1.
 */

      ptr = cs_ptr;
      while (ptr <= cs_ptre) {
        csa = *(ptr++); csb = *(ptr++);
        for (g2=1,e1=es1,e2=g3; e1<=ef1; g2++,e1++,e2+=ngens1){
          csima = g1==ngens1 && g2==ngens1 ?
                   (mult1ptr->is_accepting[csa]>0 ? csa : 0) :
                   target(dense_ip1,table1,e1,csa,dr1);
          if (csima==0)
            continue;

          csimb = g2==ngens1 && g3==ngens1 ?
                   (mult2ptr->is_accepting[csb]>0 ? csb : 0) :
                   target(dense_ip2,table2,e2,csb,dr2);
          if (csimb==0)
            continue;

          if (ht_ptrb>ht_ptre || csima> *(ht_ptre-1) ||
                               (csima== *(ht_ptre-1) &&  csimb > *ht_ptre) ){
/* We have a new state for the image subset to be added to the end */
            *(++ht_ptre) = csima; *(++ht_ptre) = csimb;
          }
          else {
            ht_chptr = ht_ptrb;
            while (*ht_chptr < csima)
              ht_chptr+=2;
            if (*ht_chptr == csima)
              while (*(ht_chptr+1) < csimb)
                ht_chptr+=2;
            if (csima < *ht_chptr || csimb < *(ht_chptr+1)) {
/* we have a new state for the image subset to be added in the middle */
              ht_ptr = ht_ptre;
              ht_ptre += 2;
              while (ht_ptr >= ht_chptr) {
                *(ht_ptr+2) = *ht_ptr;
                ht_ptr--;
              }
              *ht_chptr = csima; *(ht_chptr+1) = csimb;
            }
          }
        }
      }
      if (ht_ptre > ht_ptrb) {
        if (g1==ngens1)
          *(++ht_ptre) = 1;
        else if (g3==ngens1)
          *(++ht_ptre) = 2;
      }
      im = short_hash_locate(&ht,ht_ptre-ht_ptrb+1);
      if (im > 0) {
        if (dense_op)
          fsarow[e-1] = im;
        else  {
          fsarow[++len] = e;
          fsarow[++len] = im;
        }
        nt++;
      }
    }
    if (!dense_op)
      fsarow[0] = len++;
    fwrite((void *)fsarow,sizeof(int),(size_t)len,tempfile);
  }
  fclose(tempfile);

  ns = compositeptr->states->size = ht.num_recs;
  compositeptr->table->numTransitions = nt;

  if (print_level>=3)
    printf("    #Calculated transitions - %d states, %d transitions.\n",ns,nt);

/* Locate the accept states. This is slightly cumbersome, since a state
 * is an accept state if either the corresponding subset contains a
 * a pair (s1,s2), where s1 is and accept state of *mult1ptr and s2 an
 * accept state of *mult2ptr, or we can get from some such state to an
 * accept state by applying elements ( ($,x), (x,$ ).
 * We will need to use the array compositeptr->is_accepting.
 */

  tmalloc(compositeptr->is_accepting,boolean,ns+1);
  for (i=1;i<=ns;i++)
    compositeptr->is_accepting[i] = FALSE;
  ct = 0;
  es1 = ngens*ngens1 + 1;
  ef1 = ngens1*ngens1-1;
  for (cstate=ns;cstate>=1;cstate--) {
/* We work backwards through the states, since we wish to add on additional
 * elements at the end of the list in the hash-table - this destroys
 * later entries, but that doesn't matter, since we have already dealt
 * with them.
 */
    cs_ptr = short_hash_rec(&ht,cstate);
    reclen = short_hash_rec_len(&ht,cstate);
    if (reclen % 2 ==1)
      reclen--; /* The last entry marks the fact that this is a "padded state"*/
    cs_ptre = short_hash_rec(&ht,cstate) + reclen - 1;
/* First see if the set itself contains an accept-state */
    ptr = cs_ptr;
    while (ptr <= cs_ptre) {
      csa = *(ptr++); csb = *(ptr++);
      if (mult1ptr->is_accepting[csa] && mult2ptr->is_accepting[csb]) {
        compositeptr->is_accepting[cstate] = TRUE;
        ct++;
        break;
      }
    }
    if (compositeptr->is_accepting[cstate])
      continue;
/* Next apply generators ( ($,g2), (g2,$) ) and see if we get anything new.
 * We won't bother about having them in increasing order this time.
 */
    ptr = cs_ptr;
    while (ptr <= cs_ptre) {
      csa = *(ptr++); csb = *(ptr++);
      for (e1=es1,e2=ngens1; e1<=ef1; e1++,e2+=ngens1){
        csima = target(dense_ip1,table1,e1,csa,dr1);
        if (csima==0)
          continue;
        csimb = target(dense_ip2,table2,e2,csb,dr2);
        if (csimb==0)
          continue;

        if (csima==0 || csimb==0)
          continue;
  /* see if (csima,csimb) is accepting */
        if (mult1ptr->is_accepting[csima] && mult2ptr->is_accepting[csimb]) {
          compositeptr->is_accepting[cstate] = TRUE;
          ct++;
          break;
        }
  /* now see if it is new */
        ht_chptr = cs_ptr;
        got = FALSE;
        while (ht_chptr < cs_ptre) {
          if (csima == ht_chptr[0] && csimb == ht_chptr[1]) {
            got = TRUE;
            break;
          }
          ht_chptr+=2;
        }
        if (!got) {
  /* add (csima,csimb) to the end */
          *(++cs_ptre) = csima; *(++cs_ptre) = csimb;
        }
      }
      if (compositeptr->is_accepting[cstate])
        continue;
    }
  }

  compositeptr->num_accepting = ct;
  if (ct==1 || ct != ns) {
    tmalloc(compositeptr->accepting,int,ct+1);
    ct = 0;
    for (i=1;i<=ns;i++) if (compositeptr->is_accepting[i])
      compositeptr->accepting[++ct] = i;
  }
  tfree(compositeptr->is_accepting);
  tfree(mult1ptr->is_accepting);
  tfree(mult2ptr->is_accepting);
  short_hash_clear(&ht);
  tfree(fsarow);
  if (destroy) {
    fsa_clear(mult1ptr); fsa_clear(mult2ptr);
  }

  if (readback) {
    tempfile = fopen(compfilename,"r");
    compressed_transitions_read(compositeptr,tempfile);
    fclose(tempfile);
    unlink(compfilename);
  }

  return compositeptr;
}

fsa *
fsa_composite_int(mult1ptr,mult2ptr,op_table_type,destroy,compfilename,readback)
	fsa *mult1ptr, *mult2ptr;
	storage_type op_table_type;
	boolean destroy;
	char *compfilename;
        boolean readback;
{
 fprintf(stderr,"Sorry - fsa_composite is not yet implemented for machines.\n");
  fprintf(stderr,"with more than 65536 states.\n");
  exit(1);
}
