/* file worddiff.c 6/10/94
   Routines for adding equations to become acceptable by a word-difference
   fsa.  These will be used both by kbprog, and later programs (mult and
   newdiff) that correct the word-difference machines, but the word-reduction
   algorithm will be different in both cases - in this file it is called
   reduce_word.
*/

#include "defs.h"
#include "fsa.h"
#include "rws.h"
#include "externals.h"


extern void 	(*reduce_word)();
static char	testword[4096]; /* Used for reducing words */

/* Functions defined in this file: */
void	initialise_wd_fsa();
void 	build_wd_fsa();
boolean	add_wd_fsa();
void 	make_full_wd_fsa();
void 	clear_wd_fsa();
int	diff_no();
void    calculate_inverses();

/* Functions used in this file defined in other files: */
void	fsa_init();
void	fsa_table_init();
void	fsa_table_dptr_init();

void
initialise_wd_fsa(wd_fsaptr,alphptr,maxwdiffs)
	fsa  *wd_fsaptr;
	srec *alphptr;
	int maxwdiffs;
/* Initialise a word-difference automaton, using  *alphptr as base-alphabet.
 * First state is empty word, which is initial and only accepting state.
 */
{ int i;

  fsa_init(wd_fsaptr);

  wd_fsaptr->states->type = WORDS; 
  tmalloc(wd_fsaptr->states->names,char *,maxwdiffs+1);
  wd_fsaptr->states->alphabet_size = alphptr->size;
  for (i=1;i<=alphptr->size;i++) {
    tmalloc(wd_fsaptr->states->alphabet[i],char,strlen(alphptr->names[i])+1);
    strcpy(wd_fsaptr->states->alphabet[i],alphptr->names[i]);
  }
  wd_fsaptr->states->size=1;
/* Set up first state, which is the empty word. */
  tmalloc(wd_fsaptr->states->names[1],char,1);
  wd_fsaptr->states->names[1][0] = '\0';

  wd_fsaptr->alphabet->type = PRODUCT;
  wd_fsaptr->alphabet->size = (alphptr->size+1)*(alphptr->size+1) - 1;
  wd_fsaptr->alphabet->arity = 2;
  wd_fsaptr->alphabet->padding = '_';
  tmalloc(wd_fsaptr->alphabet->base,srec,1);
  srec_copy(wd_fsaptr->alphabet->base,alphptr);

  wd_fsaptr->num_initial = 1;
  tmalloc(wd_fsaptr->initial,int,2);
  wd_fsaptr->initial[1] = 1;
  wd_fsaptr->num_accepting = 1;
  tmalloc(wd_fsaptr->accepting,int,2);
  wd_fsaptr->accepting[1] = 1; /* Only the initial state is accepting */

  wd_fsaptr->flags[DFA] = TRUE;
  wd_fsaptr->flags[TRIM] = TRUE;

  fsa_table_init(wd_fsaptr->table,maxwdiffs,wd_fsaptr->alphabet->size);
  wd_fsaptr->table->printing_format = SPARSE;
  fsa_table_dptr_init(wd_fsaptr);
}

void
build_wd_fsa(wd_fsaptr,rws,new_wd)
	fsa *wd_fsaptr;
	rewriting_system 	rws;
	boolean *new_wd;
/* Make the word-difference machine from the rewriting system rws */
{ int i, **table;
  boolean new;
  wd_fsaptr->states->size = 1;
  table = wd_fsaptr->table->table_data_ptr;
  for (i=1;i<=wd_fsaptr->alphabet->size;i++)
    set_dense_target(table,i,1,0);

  if (new_wd) /* new_wd nonzero means that we wish to note which equations
               * give rise to changes in the table.
               */
    for (i=1;i<=rws.num_eqns;i++)
      new_wd[i] = FALSE;

/* Now add each equation to the fsa. */
  for (i=1;i<=rws.num_eqns;i++){
     new = add_wd_fsa(wd_fsaptr,rws.eqns+i,rws.inv_of,FALSE);
     if (new_wd && new)
        new_wd[i] = TRUE;
  }
  if (print_level>=2)
    printf("	#There are %d word-differences.\n",wd_fsaptr->states->size);
}

boolean
add_wd_fsa(wd_fsaptr,eqn,inv,reverse)
        fsa *wd_fsaptr;
	reduction_equation	*eqn;
	int			*inv;
	boolean			reverse;
/* Alter the word-difference machine to make it accept the equation *eqn
 * If reverse is true, then for all transitions added, the inverse transition
 * is also added.
 */
{ char *wd1, *wd2, g1, g2;
  int i, j, l1, l2, maxl, state, image, size_pba, **table, ***wd_table;
  boolean ans = FALSE;
  size_pba = 1 + wd_fsaptr->alphabet->base->size;
  wd1 = eqn->lhs; wd2 = eqn->rhs;
  l1 = strlen(wd1); l2 = strlen(wd2);
  maxl = l1>l2 ? l1 : l2;

  table = wd_fsaptr->table->table_data_ptr;
  wd_table = wd_fsaptr->table->table_data_dptr;
  state = wd_fsaptr->initial[1];
  for (i=0;i<maxl;i++){
    g1 = i>=l1 ? size_pba : wd1[i];
    g2 = i>=l2 ? size_pba : wd2[i];
    image = dense_dtarget(wd_table,g1,g2,state);
    if (image==0){
      if (g1==size_pba)
          sprintf(testword,"%s%c",wd_fsaptr->states->names[state],g2);
      else if (g2==size_pba)
          sprintf(testword,"%c%s",inv[g1],wd_fsaptr->states->names[state]);
      else
          sprintf(testword,"%c%s%c",inv[g1],wd_fsaptr->states->names[state],g2);
      (*reduce_word)(testword);
      image = diff_no(wd_fsaptr,testword);
      if (image==0){ /* new state */
        image = (++wd_fsaptr->states->size);
        if (image > wd_fsaptr->table->maxstates){
          fprintf(stderr,
            "Maximum number of word-differences exceeded. Cannot continue.\n");
          exit(1);
        }
        tmalloc(wd_fsaptr->states->names[image],char,strlen(testword)+1);
        strcpy(wd_fsaptr->states->names[image],testword);
        for (j=1;j<=wd_fsaptr->alphabet->size;j++)
           set_dense_target(table,j,image,0);
      }
      set_dense_dtarget(wd_table,g1,g2,state,image);
      ans = TRUE;
    }
    if (reverse)
      set_dense_dtarget(wd_table,inv[g1],inv[g2],image,state);
    state=image;
  }
  return ans;
}

void
make_full_wd_fsa(wd_fsaptr,inv,start_no)
	fsa *wd_fsaptr;
	int *inv;
        int start_no;
/* Close the set of word-differences under inversion, and add all possible
 * transitions, starting at state number start_no.
 */
{ int i, j, ns, n, g1, g2, size_pba, **table, ***wd_table;
  char **wdn, *ptr, *ptre, *ptri;

  if (print_level>=3)
    printf("    #Calling make_full_wd_fsa.\n");
  size_pba = 1 + wd_fsaptr->alphabet->base->size;
  ns = wd_fsaptr->states->size;
  wdn = wd_fsaptr->states->names;
  table = wd_fsaptr->table->table_data_ptr;
  wd_table = wd_fsaptr->table->table_data_dptr;

  for (i=1;i<=ns;i++){
    ptr = wdn[i];
    ptre = ptr + strlen(ptr);
/* invert this word-difference and put into testword. */
    ptri = testword;
    while (--ptre >= ptr)
      *(ptri++) = inv[*ptre];
    *ptri = '\0';
    (*reduce_word)(testword);
    
    if (diff_no(wd_fsaptr,testword)==0){ /* new state */
      n = (++wd_fsaptr->states->size);
      if (n > wd_fsaptr->table->maxstates){
        fprintf(stderr,"Too many word-differences. Increase maxwdiffs.\n");
        exit(1);
      }
      tmalloc(wd_fsaptr->states->names[n],char,strlen(testword)+1);
      strcpy(wd_fsaptr->states->names[n],testword);
      for (j=1;j<=wd_fsaptr->alphabet->size;j++)
         set_dense_target(table,j,n,0);
    }
  }

/* Now fill out table */
  ns = wd_fsaptr->states->size;
  if (start_no<1)
    start_no = 1;
  for (i=start_no;i<=ns;i++){
    for (g1=1;g1<=size_pba;g1++) for (g2=1;g2<=size_pba;g2++){
      if (g1==size_pba && g2==size_pba)
        continue; /* Don't want padding-symbol on both sides. */
      if (dense_dtarget(wd_table,g1,g2,i) == 0) {
        if (g1==size_pba)
          sprintf(testword,"%s%c",wd_fsaptr->states->names[i],g2);
        else if (g2==size_pba)
          sprintf(testword,"%c%s",inv[g1],wd_fsaptr->states->names[i]);
        else
          sprintf(testword,"%c%s%c",inv[g1],wd_fsaptr->states->names[i],g2);
        (*reduce_word)(testword);
        if (n=diff_no(wd_fsaptr,testword))
          set_dense_dtarget(wd_table,g1,g2,i,n);
        if (n>0 && n<start_no)
          set_dense_dtarget(wd_table,inv[g1],inv[g2],n,i);
      }
    }
  }
}


void
clear_wd_fsa(wd_fsaptr)
	fsa *wd_fsaptr;
/* Clear the state-names for all states except the first. */
{ int i, ns;
  char **wdn;
  ns = wd_fsaptr->states->size;
  wdn = wd_fsaptr->states->names;

  for (i=2;i<=ns;i++)
    tfree(wdn[i]);
}


int
diff_no(wd_fsaptr,w)
	fsa *wd_fsaptr;
	char *w;
/* See if w is in the list of word-differences (state-names of wd_diff).
 * If so, return the number of the state.
 * If not, return 0.
 * This is done by a simple search. If it turns out to be too slow, we can use
 * hashing later.
 */
{ int i, ns;
  char **wdn;
  ns = wd_fsaptr->states->size;
  wdn = wd_fsaptr->states->names;

  for (i=1;i<=ns;i++) if (strcmp(w,wdn[i])==0)
    return i;

  return 0;
}

void
calculate_inverses(invptr,ngens)
	int **invptr;
	int ngens;
/* Use reduction in the word-difference machine to calculate
 * inverses of generators.
 */
{ int i, j;
  tmalloc(*invptr,int,ngens+2);
  for (i=1;i<=ngens;i++) for (j=1;j<=ngens;j++) {
    testword[0] = i;
    testword[1] = j;
    testword[2] = '\0';
    (*reduce_word)(testword);
    if (strlen(testword)==0) {
      (*invptr)[i] = j;
      break;
    }
  }
  (*invptr)[ngens+1] = ngens+1; /* to handle the padding symbol */
}
