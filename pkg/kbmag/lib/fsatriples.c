/* file fsatriples.c */
/* 14. 11. 94.
 * This file contains the routine fsa_triples, which is used to construct the
 * general product fsa in an automatic group.
 * For general methodology, see comment at top of file fsalogic.c.
*/

#include <stdio.h>
#include "defs.h"
#include "fsa.h"
#include "hash.h"
#include "rws.h"
#include "externals.h"

/* Functions defined in this file: */
fsa * fsa_triples();


/* Functions used in this file and defined elsewhere */
boolean srec_equal();
boolean table_equal();
boolean fsa_equal();
int sparse_target();
void fsa_init();
void fsa_set_is_initial();
void fsa_set_is_accepting();
void srec_copy();
void fsa_clear();
void short_hash_init();
int short_hash_locate();
void short_hash_clear();
unsigned short* short_hash_rec();

fsa *
fsa_triples(waptr,diffptr,op_table_type,destroy,tempfilename,eqnptr,maxeqns,
						     idlabel,eqnstop,readback)
	fsa *waptr, *diffptr;
	storage_type op_table_type;
	boolean destroy;
	char *tempfilename;
	reduction_equation *eqnptr;
	int maxeqns;
        boolean idlabel;
        boolean eqnstop;
	boolean readback;
/* *waptr is assumed to be the word-acceptor of an automatic group.
 * (In particular, all states should be accepting.)
 * *diffptr is assumed to be a word-difference machine of the same automatic
 *  group. Both are assumed to be stored in dense-format.
 * This routine constructs the fsa of which the states are triples (s1,s2,d),
 * with s1 and s2 states of *waptr and d a state of *diffptr.
 * (More precisely, if *waptr has n states, then s1 and s2 may also be equal
 * to n+1, meaning that the end of string symbol has been read on lhs or rhs.
 * The alphabet is 2-variable with base the alphabet of *waptr
 * (i.e. the same alphabet as *diffptr).
 * The alphabet member (g1,g2) maps (s1,s2,d) to (s1^g1,s2^g2,d^(g1,g2))
 * if all three components are nonzero, and to zero otherwise.
 * The transition-table of the resulting fsa is output in the usual way to
 * file tempfilename with table-type specified by op_table_type, before
 * minimisation.
 * Short hash-tables will be used, so this routine won't work if *waptr
 * or *diffptr has more than 65536 states.
 *
 * If during the construction, a nontrivial equation between two words is
 * discovered as a result of encountering the identity word-difference,
 * then the word-acceptor *waptr must be accepting both of these words
 * which represent the same group-element, and must therefore be wrong.
 * The procedure therefore aborts without returning an fsa.
 * If the maxeqns is greater than zero, then a maximum of maxeqns such
 * equations are returned as eqnptr[i] - in order to do this, it is necessary to
 * store the defining transitions of the states as we proceed.
 * idlabel is set true if the identity is to be used as a label.
 */
{
  int **watable, ***difftable, identity, ngens, ngens1, nswa1, ne, ns, *fsarow,
      nt, cstate, cswa1, cswa2, csdiff, im, i, e, len, ct, *tab_ptr,
      **table, bstate, bigger, numeqns;
  unsigned short *ht_ptr;
  boolean  dense_op;
  fsa *triples;
  srec *labels;
  short_hash_table ht;
  FILE *tempfile, *fopen();
  char g1, g2, bg1, bg2;
  int maxv = 65536;
  struct vertexd {
     char g1;
     char g2;
     int  state;
  } *definition, *newdef; 
/* This is used to store the defining transition for the states of *triples.
 * If definition[i] = v, then state i is defined by the transition from
 * state v.state, with generator (v.g1,v.g2).
 * State 1 does not have a definition.
 */

  if (print_level>=3)
    printf("    #Calling fsa_triples.\n");

  if (!waptr->flags[DFA] || !diffptr->flags[DFA]){
    fprintf(stderr,"Error: fsa__triples only applies to DFA's.\n");
    exit(1);
  }
  if (waptr->alphabet->type!=IDENTIFIERS) {
    fprintf(stderr,"Error in fsa_triples: first fsa has wrong type.\n");
    exit(1);
  }
  if (waptr->num_accepting != waptr->states->size) {
    fprintf(stderr,
       "Error in fsa_triples: first fsa should be a word-acceptor.\n");
    exit(1);
  }
  if (diffptr->alphabet->type!=PRODUCT || diffptr->alphabet->arity!=2) {
    fprintf(stderr,"Error in fsa_triples: second fsa must be 2-variable.\n");
    exit(1);
  }
  if (diffptr->states->type!=WORDS) {
    fprintf(stderr,
       "Error in fsa_triples: second fsa must be word-difference type.\n");
    exit(1);
  }
  if (!srec_equal(diffptr->alphabet->base,waptr->alphabet)) {
    fprintf(stderr,"Error in fsa_triples: fsa's alphabet's don't match.\n");
    exit(1);
  }
  if (waptr->states->size>MAXSHORT || diffptr->states->size>MAXSHORT) {
    fprintf(stderr,
       "Error in fsa_triples: one of the fsa's has too many states.\n");
    exit(1);
  }

  fsa_table_dptr_init(diffptr);

  if (maxeqns>0) {
 /* We need to remember vertex definitions */ 
    tmalloc(definition, struct vertexd, maxv);
    ns = 1;
  }

  tmalloc(triples,fsa,1);
  fsa_init(triples);
  srec_copy(triples->alphabet,diffptr->alphabet);
  triples->flags[DFA] = TRUE;
  triples->flags[ACCESSIBLE] = TRUE;
  triples->flags[BFS] = TRUE;
  triples->num_accepting=0;
/* In fact there will be lots of different categories of accept-states -
 * one for each generator - they will be recorded as labels of the states.
 */

  ngens = waptr->alphabet->size;
  ngens1 = ngens+1;
  ne = diffptr->alphabet->size;
  nswa1 = waptr->states->size + 1;

  if (ne != ngens1*ngens1-1) {
   fprintf(stderr,
       "Error: in a 2-variable fsa, alphabet size should = ngens^2 - 1.\n");
    exit(1);
  }

  watable = waptr->table->table_data_ptr;
  difftable = diffptr->table->table_data_dptr;

  dense_op = op_table_type==DENSE;

  triples->states->type = LABELED;
  tmalloc(triples->states->labels,srec,1);
  if (!idlabel)
    srec_copy(triples->states->labels,waptr->alphabet);
  else {
/* Here we want one extra label, for the identity generator "_",
 * so we can't just use srec_copy straight from waptr->alphabet
 */
    labels = triples->states->labels;
    labels->size = ngens+1;
    labels->type = IDENTIFIERS;
    tmalloc(labels->names,char *,ngens+2);
    for (i=1;i<=ngens;i++) {
      tmalloc(labels->names[i],char,strlen(waptr->alphabet->names[i])+1);
      strcpy(labels->names[i],waptr->alphabet->names[i]);
    }
    tmalloc(labels->names[ngens+1],char,2);
    strcpy(labels->names[ngens+1],"_");
  }

  triples->num_initial = 1;
  tmalloc(triples->initial,int,2);
  triples->initial[1] = 1;
  triples->table->table_type = op_table_type;
  triples->table->denserows = 0;
  triples->table->printing_format = op_table_type;
  if (!readback) {
    tmalloc(triples->table->filename,char,strlen(tempfilename)+1);
    strcpy(triples->table->filename,tempfilename);
  }
  
  short_hash_init(&ht,TRUE,3,0,0);
  ht_ptr = ht.current_ptr;
  ht_ptr[0] = waptr->initial[1];
  ht_ptr[1] = waptr->initial[1];
  ht_ptr[2] =identity = diffptr->initial[1];
  im = short_hash_locate(&ht,3);
  if (im!=1) {
    fprintf(stderr,"Hash-initialisation problem in fsa_triples.\n");
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
  nt = 0; /* Number of transitions in triples */

  numeqns = 0; /* this becomes nonzero when we have started collecting
                * equations of equal words both accepted by word-acceptor.
                */
  while (++cstate <= ht.num_recs) {
    if (print_level>=3) {
      if ((cstate<=1000 && cstate%100==0)||(cstate<=10000 && cstate%1000==0)||
          (cstate<=100000 && cstate%5000==0) || cstate%50000==0)
       printf("    #cstate = %d;  number of states = %d.\n",cstate,ht.num_recs);
    }
    ht_ptr = short_hash_rec(&ht,cstate);
    cswa1 = ht_ptr[0]; cswa2 = ht_ptr[1];
    csdiff = ht_ptr[2];
    if (!dense_op)
      len = 0;
    e = 0; /* e is the num,ber of the edge corresponding to the pair (g1,g2) */
    for (g1=1;g1<=ngens1;g1++) for (g2=1;g2<=ngens1;g2++) {
      e++;
/* Calculate action of generator-pair (g1,g2) on state cstate */
      if (g1==ngens1 && g2==ngens1)
        continue;
      ht_ptr = ht.current_ptr;
      ht_ptr[2] = dense_dtarget(difftable,g1,g2,csdiff);
      if (ht_ptr[2]==0)
        im=0;
      else {
        ht_ptr[0] = g1>ngens ? nswa1 : cswa1==nswa1 ? 0 :
                                          dense_target(watable,g1,cswa1);
        if (ht_ptr[0]==0)
          im=0;
        else {
          ht_ptr[1] = g2>ngens ? nswa1 : cswa2==nswa1 ? 0 :
                                          dense_target(watable,g2,cswa2);
          if (ht_ptr[1]==0)
            im=0;
          else {
            if (eqnstop && ht_ptr[2]==identity && csdiff!=identity) {
 /* This means that we have found a new equation between two distinct words 
  * accepted by the word-acceptor *gpwa, and so  *gpwa must have been wrong.
  */
             if (print_level>0 && numeqns==0)
               printf(
            "#Equation found between two words accepted by word-acceptor.\n");
              if (maxeqns>0) { /* We reconstruct the equation explicitly */
                /* First we calculate the length of the equation */
                if (print_level>=3)
                  printf("    #Calculating equation number %d.\n",numeqns+1);
                len = 1;
                bg1 = g1; bg2 = g2; bstate = cstate;
                bigger = g2>ngens ? 1 : g1>ngens ? 2 : 0;
                /* bigger=1 or 2 means resp. lhs/rhs larger in shortlex order */
                while (bstate != 1) {
                  len++;
                  bg1 = definition[bstate].g1;
                  bg2 = definition[bstate].g2;
                  bstate = definition[bstate].state;
                }
                if (bigger==0)
                  bigger = bg1>bg2 ? 1 : 2;

                /* Now we allocate space for it and store it */
                tmalloc(eqnptr[numeqns].lhs,char,len+1);
                tmalloc(eqnptr[numeqns].rhs,char,len+1);
                eqnptr[numeqns].lhs[len] = eqnptr[numeqns].rhs[len] = '\0';
                bg1 = g1; bg2 = g2; bstate = cstate;
                while (1) {
                  len--;
                  if (bigger==1) {
                     eqnptr[numeqns].lhs[len] = bg1>ngens ? '\0' : bg1;
                     eqnptr[numeqns].rhs[len] = bg2>ngens ? '\0' : bg2;
                  }
                  else {
                     eqnptr[numeqns].rhs[len] = bg1>ngens ? '\0' : bg1;
                     eqnptr[numeqns].lhs[len] = bg2>ngens ? '\0' : bg2;
                  }
                  if (bstate==1)
                    break;
                  bg1 = definition[bstate].g1;
                  bg2 = definition[bstate].g2;
                  bstate = definition[bstate].state;
                }
              }

              if (numeqns==0) {
    /* We are no longer constructing the fsa, so we no longer need the file. */
                fclose(tempfile);
                unlink(tempfilename);
              }
              numeqns++;
              if (numeqns>=maxeqns) { /* exit */
                if (print_level>=2 && maxeqns>0)
                    printf("  #Found %d new equations - aborting.\n",maxeqns);
                short_hash_clear(&ht);
                tfree(fsarow);
                fsa_clear(triples);
                tfree(triples);
                if (maxeqns>0)
                  tfree(definition);
                if (destroy) {
                  fsa_clear(waptr); fsa_clear(diffptr);
                }
                return 0;
              }
              else eqnptr[numeqns].lhs = 0; /* to mark how many we have later */
            }
            im = short_hash_locate(&ht,3);
            if (maxeqns>0 && im>ns) {
              ns++;
              if (ns==maxv) { /* need room for more definitions */
                if (print_level>=3)
                 printf("    #Allocating more space for vertex definitions.\n");
                tmalloc(newdef,struct vertexd,2*maxv);
                for (i=1;i<maxv;i++) 
                  newdef[i] = definition[i];
                tfree(definition);
                definition = newdef;
                maxv *= 2;
              }
              definition[ns].g1 = g1;
              definition[ns].g2 = g2;
              definition[ns].state = cstate;
            }
          }
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
    if (numeqns==0)
      fwrite((void *)fsarow,sizeof(int),(size_t)len,tempfile);
  }  /*while (++cstate <= ht.num_recs) */

  if (numeqns>0) {
    short_hash_clear(&ht);
    tfree(fsarow);
    fsa_clear(triples);
    tfree(triples);
    tfree(definition);
    if (destroy) {
      fsa_clear(waptr); fsa_clear(diffptr);
    }
    if (print_level>=2)
      printf("  #Found %d new equations - aborting with algorithm complete.\n",
                  numeqns);
    return 0;
  }
    
  fclose(tempfile);

  ns = triples->states->size = ht.num_recs;
  triples->table->numTransitions = nt;

  if (print_level>=3) {
    printf("    #Calculated transitions - %d states, %d transitions.\n",ns,nt);
    printf("    #Now calculating state labels.\n");
  }
  tmalloc(triples->states->setToLabels,setToLabelsType,ns+1);
  triples->states->setToLabels[0]=0;
  for (i=1;i<=ns;i++) {
    ht_ptr = short_hash_rec(&ht,i);
    csdiff = ht_ptr[2];
    len = strlen(diffptr->states->names[csdiff]);
    if (len==0 && idlabel)
      triples->states->setToLabels[i]=ngens+1;
    else if (len==1)
      triples->states->setToLabels[i]=diffptr->states->names[csdiff][0];
    else
      triples->states->setToLabels[i] = 0;
  }

  short_hash_clear(&ht);
  tfree(fsarow);
  if (maxeqns>0)
    tfree(definition);
/* Now read the transition table back in */
  if (readback) {
    tempfile = fopen(tempfilename,"r");
    compressed_transitions_read(triples,tempfile);
    fclose(tempfile);
    unlink(tempfilename);
  }
  if (destroy) {
    fsa_clear(waptr); fsa_clear(diffptr);
  }

  return triples;
}
