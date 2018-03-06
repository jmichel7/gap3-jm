/* file fsa.c - 16. 9. 94.
 * 17.11.94. - implemented labeled srec type and function fsa_labeled_minimize.
 * This file contains routines for manipulating finite state automata
 * The functions in this file currently only deal with deterministic fsa's
*/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "fsa.h"
#include "hash.h"
#include "externals.h"

char *type_names[] = {"simple", "identifiers", "words", "strings",
                      "labeled", "product" };
char *flag_names[] = {"DFA","NFA","minimized","BFS","accessible","trim","RWS"};

/* Functions defined in this file: */
int sparse_target();
void fsa_init();
void fsa_table_init();
void fsa_set_is_initial();
void fsa_set_is_accepting();
void fsa_table_dptr_init();
void srec_copy();
void table_copy();
void fsa_copy();
void srec_clear();
void table_clear();
void fsa_clear();
void fsa_delete_state();
void fsa_permute_states();
void fsa_clear_rws();
void fsa_make_accessible();
void fsa_minimize();
void fsa_labeled_minimize();
void fsa_bfs();
int fsa_count();
boolean fsa_enumerate();

/* Functions used in this file and defined elsewhere */
void hash_init();
void hash_clear();
int *hash_rec();
int hash_rec_len();
int hash_locate();
void add_word_to_buffer();
void add_iword_to_buffer();
int int_len();
fsa * fsa_and();

int
sparse_target(g,p1,p2)
	int g;
	int *p1, *p2;
/* p1 points to the beginning of the row of targets from the required state
 * in a sparsely stored fsa, and p2 to the location beyond the end.
 * g is the generator number for which we seek the target.
 * If we don't find one, then there isn't one, and we return 0.
 */
{ while (p1 < p2){
    if (*p1 == g)
      return *(p1+1);
    p1+=2;
  }
  return 0;
}

void
fsa_init(fsaptr)
	fsa *fsaptr;
/* Allocate space for the sub-structures of fsaptr, and zero all pointers. */
{ int i;
  fsaptr->initial = 0;
  fsaptr->accepting = 0;
  fsaptr->is_accepting = 0;
  fsaptr->is_initial = 0;

  tmalloc(fsaptr->alphabet,srec,1);
  fsaptr->alphabet->names = 0;
  fsaptr->alphabet->base = 0;
  fsaptr->alphabet->type = SIMPLE;
  fsaptr->alphabet->size = 0;
  fsaptr->alphabet->labels = 0;
  fsaptr->alphabet->setToLabels = 0;

  tmalloc(fsaptr->states,srec,1);
  fsaptr->states->names = 0;
  fsaptr->states->base = 0;
  fsaptr->states->type = SIMPLE;
  fsaptr->states->size = 0;
  fsaptr->states->labels = 0;
  fsaptr->states->setToLabels = 0;

  for (i=0;i<num_flag_strings;i++)
    fsaptr->flags[i] = FALSE;
  tmalloc(fsaptr->table,table_struc,1);
  fsaptr->table->filename = 0;
  fsaptr->table->table_data_ptr = 0;
  fsaptr->table->table_data_dptr = 0;
}

void
fsa_table_init(tableptr,maxstates,ne)
	table_struc *tableptr;
	int maxstates;
	int ne;
/* Intialize the table_data and table_data_ptr fields of fsaptr->table,
 * for maxstates states, with table-type dense.
 * ne is the size of the alphabet.
 */
{ int i;
  tableptr->maxstates = maxstates;
  tableptr->table_type =  DENSE;
  tmalloc(tableptr->table_data_ptr,int *,ne+1);
  tmalloc(tableptr->table_data_ptr[0],int,maxstates*ne);
  for (i=1; i <= ne; i++)
      tableptr->table_data_ptr[i] =
          tableptr->table_data_ptr[0] + (i-1)*maxstates - 1;
}

void
fsa_set_is_initial(fsaptr)
/* Set the is_initial field in fsa *fsaptr - should be freed after use  */
	fsa *fsaptr;
{ int ns, i;
  ns = fsaptr->states->size;
  tmalloc(fsaptr->is_initial,boolean,ns+1);
  fsaptr->is_initial[0] = FALSE;
  if (fsaptr->num_initial == ns) {
    for (i=1;i<=ns;i++)
      fsaptr->is_initial[i] = TRUE;
  }
  else {
    for (i=1;i<=ns;i++)
      fsaptr->is_initial[i] = FALSE;
    for (i=1;i<=fsaptr->num_initial;i++)
      fsaptr->is_initial[fsaptr->initial[i]] = TRUE;
  }
}
 
void
fsa_set_is_accepting(fsaptr)
/* Set the is_accepting field in fsa *fsaptr - should be freed after use */
	fsa *fsaptr;
{ int ns, i;
  ns = fsaptr->states->size;
  tmalloc(fsaptr->is_accepting,boolean,ns+1);
  fsaptr->is_accepting[0] = FALSE;
  if (fsaptr->num_accepting == ns) {
    for (i=1;i<=ns;i++)
      fsaptr->is_accepting[i] = TRUE;
  }
  else {
    for (i=1;i<=ns;i++)
      fsaptr->is_accepting[i] = FALSE;
    for (i=1;i<=fsaptr->num_accepting;i++)
      fsaptr->is_accepting[fsaptr->accepting[i]] = TRUE;
  }
}
 
void
fsa_table_dptr_init(fsaptr)
/* Set the field fsaptr->table->table_data_dptr.
 * This is used for fast access to table entries in 2-variable machines.
 */
	fsa *fsaptr;
{ int i, ngens;
  if (fsaptr->alphabet->type!=PRODUCT || fsaptr->alphabet->arity!=2) {
    fprintf(stderr,"Error in fsa_table_dptr_init: fsa must be 2-variable.\n");
    exit(1);
  }
  if (fsaptr->table->table_type != DENSE)  {
    fprintf(stderr,
            "Error in fsa_table_dptr_init: storage-type must be dense.\n");
    exit(1);
  }
  ngens = fsaptr->alphabet->base->size;
  if (fsaptr->table->table_data_dptr == 0) {
    tmalloc(fsaptr->table->table_data_dptr,int **,ngens+2);
    for (i=1;i<=ngens+1;i++)
      fsaptr->table->table_data_dptr[i] =
                         fsaptr->table->table_data_ptr + (i-1)*(ngens+1);
  }
}

void
srec_clear(srptr)
/* Deallocate all space used by the set-record fsaptr */
	srec *srptr;
{ int i;
  if (srptr==0)
    return;
  if ((srptr->type==IDENTIFIERS || srptr->type==STRINGS ||srptr->type==WORDS)
               && srptr->names){
    for (i=1;i<=srptr->size;i++)
      tfree(srptr->names[i]);
    tfree(srptr->names);
  }
  if (srptr->type==PRODUCT && srptr->base){
    srec_clear(srptr->base);
    tfree(srptr->base);
  }
  if (srptr->type==WORDS)
    for (i=1;i<=srptr->alphabet_size;i++)
      tfree(srptr->alphabet[i]);
  if (srptr->type==LABELED) {
    srec_clear(srptr->labels);
    tfree(srptr->labels);
    tfree(srptr->setToLabels);
  }
}

void
srec_copy(srptr1,srptr2)
	srec *srptr1, *srptr2;
/* Copy the information in *srptr2 to *srptr1 (same way round as strcpy!)
 * It is assumed that srptr1 points to a zeroed set-record.
 */
{ int i;
  srptr1->type = srptr2->type;
  srptr1->size = srptr2->size;
  if (srptr1->type==IDENTIFIERS || srptr1->type==STRINGS ||
                                   srptr1->type==WORDS) {
    tmalloc(srptr1->names,char *,srptr1->size+1);
    for (i=1;i<=srptr1->size;i++) {
      tmalloc(srptr1->names[i],char,strlen(srptr2->names[i])+1);
      strcpy(srptr1->names[i],srptr2->names[i]);
    }
  }
  if (srptr1->type==WORDS) {
    srptr1->alphabet_size = srptr2->alphabet_size;
    for (i=1;i<=srptr1->alphabet_size;i++) {
      tmalloc(srptr1->alphabet[i],char,strlen(srptr2->alphabet[i])+1);
      strcpy(srptr1->alphabet[i],srptr2->alphabet[i]);
    }
  }
  if (srptr1->type==PRODUCT) {
    tmalloc(srptr1->base,srec,1);
    srec_copy(srptr1->base,srptr2->base);
    srptr1->arity = srptr2->arity;
    srptr1->padding = srptr2->padding;
  }
  if (srptr1->type==LABELED) {
    tmalloc(srptr1->labels,srec,1);
    srec_copy(srptr1->labels,srptr2->labels);
    tmalloc(srptr1->setToLabels,setToLabelsType,srptr1->size+1);
    for (i=1;i<=srptr1->size;i++)
      srptr1->setToLabels[i] = srptr2->setToLabels[i];
  }
}


void
table_copy(tableptr1,tableptr2,ne,ns)
	table_struc *tableptr1, *tableptr2;
	int ne, ns;
/* Copy the information in *tableptr2 to *tableptr1 (same way round as strcpy!)
 * It is assumed that tableptr1 points to a zeroed table_struc.
 * ne and ns are the sizes of the alphabet and state-set.
 */
{ int i, j, maxstates, dr, space, *row1, *row2, *ptr1, *ptr2, *ptr2e;
  if (tableptr2->filename) {
    tmalloc(tableptr1->filename,char,strlen(tableptr2->filename)+1);
    strcpy(tableptr1->filename,tableptr2->filename);
  }
  tableptr1->table_type = tableptr2->table_type;
  tableptr1->printing_format = tableptr2->printing_format;
  tableptr1->numTransitions = tableptr2->numTransitions;
  dr = 
  tableptr1->denserows = tableptr2->denserows;
  maxstates =
  tableptr1->maxstates = tableptr2->maxstates;
  
  if (tableptr1->table_type == DENSE) {
    tmalloc(tableptr1->table_data_ptr,int *,ne+1);
    tmalloc(tableptr1->table_data_ptr[0],int,ne*maxstates);
    for (i=1; i <= ne; i++) {
      row2 = tableptr2->table_data_ptr[i];
      row1 = tableptr1->table_data_ptr[i] =
          tableptr1->table_data_ptr[0] + (i-1)*maxstates - 1;
      for (j=1;j<=ns;j++)
        row1[j] = row2[j];
    }
  }
  else {
    tmalloc(tableptr1->table_data_ptr,int *,ns+2);
    if (dr>0) {
      tmalloc(tableptr1->table_data_ptr[0],int,ne*dr);
      for (i=1; i <= dr; i++) {
        row2 = tableptr2->table_data_ptr[i];
        row1 = tableptr1->table_data_ptr[i] =
          tableptr1->table_data_ptr[0] + (i-1)*ne - 1;
        for (j=1;j<=ne;j++)
          row1[j] = row2[j];
      }
    }
    if (ns>dr) {
      space = tableptr2->table_data_ptr[ns+1]-tableptr2->table_data_ptr[dr+1];
      tmalloc(tableptr1->table_data_ptr[dr+1],int,space+1);
      ptr1 = tableptr1->table_data_ptr[dr+1];
      ptr2 = tableptr2->table_data_ptr[dr+1];
      for (i=dr+1;i<=ns;i++) {
        tableptr1->table_data_ptr[i] = ptr1;
        ptr2e = tableptr2->table_data_ptr[i+1];
        while (ptr2 < ptr2e)
          *(ptr1++) = *(ptr2++);
      }
      tableptr1->table_data_ptr[ns+1] = ptr1;
    }
  }
}


void
fsa_copy(fsaptr1,fsaptr2)
	fsa *fsaptr1, *fsaptr2;
/* Copy the information in *fsaptr2 to *fsaptr1 (same way round as strcpy!)
 * It is assumed that fsaptr1 points to a zeroed fsa.
 */
{ int i, j, ne, ns, maxstates, *row1, *row2, *ptr1, *ptr2;
  fsa_init(fsaptr1);
  srec_copy(fsaptr1->alphabet,fsaptr2->alphabet);
  srec_copy(fsaptr1->states,fsaptr2->states);
  ne = fsaptr1->alphabet->size;
  ns = fsaptr1->states->size;

  fsaptr1->num_initial = fsaptr2->num_initial;
  if (ns==1 || fsaptr1->num_initial!=ns) {
    tmalloc(fsaptr1->initial,int,fsaptr1->num_initial+1);
    for (i=1;i<=fsaptr1->num_initial;i++)
      fsaptr1->initial[i] = fsaptr2->initial[i];
  }
  fsaptr1->num_accepting = fsaptr2->num_accepting;
  if (ns==1 || fsaptr1->num_accepting!=ns) {
    tmalloc(fsaptr1->accepting,int,fsaptr1->num_accepting+1);
    for (i=1;i<=fsaptr1->num_accepting;i++)
      fsaptr1->accepting[i] = fsaptr2->accepting[i];
  }

  for (i=0;i<num_flag_strings;i++)
    fsaptr1->flags[i] = fsaptr2->flags[i];

  table_copy(fsaptr1->table,fsaptr2->table,ne,ns);
}

void
table_clear(tableptr,ns)
/* Deallocate all space used by the table-struc fsaptr */
	table_struc *tableptr;
	int ns;
{
  if (tableptr==0)
    return;
  tfree(tableptr->filename);
  if (tableptr->table_data_ptr) {
    tfree(tableptr->table_data_ptr[0]);
    if (tableptr->table_type==SPARSE &&
                            tableptr->denserows>0 && tableptr->denserows<ns)
    tfree(tableptr->table_data_ptr[tableptr->denserows+1]);
    tfree(tableptr->table_data_ptr);
    tfree(tableptr->table_data_dptr);
  }
}

void
fsa_clear(fsaptr)
/* Deallocate all space used by the fsa *fsaptr */
	fsa	*fsaptr;
{ int ns;
  if (fsaptr==0)
    return;
  srec_clear(fsaptr->states);
  tfree(fsaptr->states);
  srec_clear(fsaptr->alphabet);
  tfree(fsaptr->alphabet);
  tfree(fsaptr->initial);
  tfree(fsaptr->accepting);
  tfree(fsaptr->is_accepting);
  tfree(fsaptr->is_initial);
  ns = fsaptr->states ? fsaptr->states->size : 0;
  table_clear(fsaptr->table,ns);
  tfree(fsaptr->table);
}

void
fsa_delete_state(fsaptr,stateno)
/* Delete state number stateno from the fsa *fsaptr - works for all fsa's. */
	fsa     *fsaptr;
	int	stateno;
{ int ns, ne, **table, *ptr, *ptr2, *ptr2e, i, j, k, l;
  srec *srptr;

  if (fsaptr->table->table_type==SPARSE && fsaptr->table->denserows>0) {
    fprintf(stderr,
"Sorry: fsa_delete_state not available for sparse storage with dense rows.\n");
    exit(1);
  }

  if (print_level>=3)
    printf("    #Calling fsa_delete_state on state number %d.\n",stateno);
  ne = fsaptr->alphabet->size;
  srptr = fsaptr->states;
  ns = srptr->size;
  if (stateno<=0 || stateno>ns){
    fprintf(stderr,"Error in fsa_delete_stateno - invalid state number.\n");
    exit(1);
  }
  if (srptr->type==IDENTIFIERS || srptr->type==STRINGS ||srptr->type==WORDS){
    tfree(srptr->names[stateno]);
    for (i=stateno;i<ns;i++)
      srptr->names[i] = srptr->names[i+1];
    srptr->names[ns]=0;
  }
  if (srptr->type==LABELED)
    for (i=stateno;i<ns;i++)
      srptr->setToLabels[i] = srptr->setToLabels[i+1];

  if (fsaptr->num_initial==ns)
    fsaptr->num_initial--;
  else
    for (i=1;i<=fsaptr->num_initial;i++){
      k = fsaptr->initial[i];
      if (k==stateno){
        for (j=i;j<fsaptr->num_initial;j++) {
          l = fsaptr->initial[j+1];
          fsaptr->initial[j] = l>stateno ? l-1 : l;
        }
        fsaptr->num_initial--;
        break;
      }
      else if (k>stateno)
        fsaptr->initial[i] = k-1;
    }

  if (fsaptr->num_accepting==ns) {
    fsaptr->num_accepting--;
    if (ns==2) {
      tmalloc(fsaptr->accepting,int,2);
      fsaptr->accepting[1] = 1;
    }
  }
  else
    for (i=1;i<=fsaptr->num_accepting;i++) {
      k = fsaptr->accepting[i];
      if (k==stateno){
        for (j=i;j<fsaptr->num_accepting;j++){
          l = fsaptr->accepting[j+1];
          fsaptr->accepting[j] = l>stateno ? l-1 : l;
        }
        fsaptr->num_accepting--;
        break;
      }
      else if (k>stateno)
        fsaptr->accepting[i] = k-1;
    }
  
  fsaptr->flags[NFA] = FALSE;
  fsaptr->flags[ACCESSIBLE] = FALSE;
  fsaptr->flags[TRIM] = FALSE;
  fsaptr->flags[BFS] = FALSE;
  fsaptr->flags[MINIMIZED] = FALSE;

  table = fsaptr->table->table_data_ptr;
  if (fsaptr->table->table_type==DENSE){
    for (j=1;j<=ne;j++) {
      for (i=1;i<stateno;i++)
        if (table[j][i] > stateno)
          table[j][i]--;
        else if (table[j][i] == stateno)
          table[j][i]=0;
      for (i=stateno;i<ns;i++) {
         k=table[j][i+1];
         table[j][i] = k<stateno ? k : k==stateno ? 0 : k-1;
      }
    }
  }
  else {
    ptr = fsaptr->table->table_data_ptr[0];
    ptr2 = ptr;
    for (i=1;i<stateno;i++) {
      table[i] = ptr;
      ptr2e = table[i+1];
      while (ptr2<ptr2e){
        if (*(ptr2+1) != stateno) {
          *(ptr++) = *(ptr2++);
          *(ptr++) = *ptr2>stateno ? *ptr2-1 : *ptr2;
          *ptr2++;
        }
        else
          ptr2+=2;
      }
    }
    ptr2 = table[stateno+1];
    for (i=stateno;i<ns;i++) {
      table[i] = ptr;
      ptr2e = table[i+2];
      while (ptr2<ptr2e){
        if (*(ptr2+1) != stateno) {
          *(ptr++) = *(ptr2++);
          *(ptr++) = *ptr2>stateno ? *ptr2-1 : *ptr2;
          *ptr2++;
        }
        else
          ptr2+=2;
      }
    }
    table[ns] = ptr;
  }
  fsaptr->states->size--;
}

void
fsa_permute_states(fsaptr,perm)
/* permute the states of fsa *fsaptr, using perm - works for all fsa's
 * perm should be a permutation of the integers 1 to ns - this is not checked
 * ALSO - perm[0] should be 0 !!
 */
	fsa     *fsaptr;
	int	*perm;
{ int ns, ne, *newtable, **newtableptr, **table, *perminv,
      *ptr, *ptr2, *ptr2e, i, j;
  srec *srptr;
  char ** newnames;
  setToLabelsType *newsetToLabels;

  if (fsaptr->table->table_type==SPARSE && fsaptr->table->denserows>0) {
    fprintf(stderr,
"Sorry: fsa_permute_state not available for sparse storage with dense rows.\n");
    exit(1);
  }
  if (print_level>=3)
    printf("    #Calling fsa_permute_states.\n");
  perm[0] = 0; /* let's make sure! */
  ne = fsaptr->alphabet->size;
  srptr = fsaptr->states;
  ns = srptr->size;
  if (srptr->type==IDENTIFIERS || srptr->type==STRINGS ||srptr->type==WORDS){
    tmalloc(newnames,char *,ns+1);
    for (i=1;i<=ns;i++)
      newnames[perm[i]] = srptr->names[i];
    tfree(srptr->names);
    srptr->names = newnames;
  }
  if (srptr->type==LABELED) {
    tmalloc(newsetToLabels,setToLabelsType,ns+1);
    for (i=1;i<=ns;i++)
      newsetToLabels[perm[i]] = srptr->setToLabels[i];
    tfree(srptr->setToLabels);
    srptr->setToLabels = newsetToLabels;
  }

  if (fsaptr->num_initial!=ns)
    for (i=1;i<=fsaptr->num_initial;i++)
      fsaptr->initial[i] = perm[fsaptr->initial[i]];

  if (fsaptr->num_accepting!=ns)
    for (i=1;i<=fsaptr->num_accepting;i++)
      fsaptr->accepting[i] = perm[fsaptr->accepting[i]];

  fsaptr->flags[BFS] = FALSE;

  table = fsaptr->table->table_data_ptr;
  if (fsaptr->table->table_type==DENSE){
    perm[0] = 0; /* just to make sure! */
    tmalloc(newtable,int,ns*ne);
    for (j=1;j<=ne;j++) {
      ptr = newtable + (j-1)*ns - 1;
      for (i=1;i<=ns;i++)
        ptr[perm[i]] = perm[table[j][i]];
      table[j] = ptr;
    }
    tfree(fsaptr->table->table_data_ptr[0]);
    fsaptr->table->table_data_ptr[0]  = newtable;
  }
  else {
/* We need the inverse of perm in this case */
    tmalloc(perminv,int,ns+1);
    for (i=1;i<=ns;i++)
      perminv[perm[i]] = i;
    tmalloc(newtable,int,table[ns+1]-table[1]);
    tmalloc(newtableptr,int *,ns+2);
    ptr = newtable;
    for (i=1;i<=ns;i++) {
      newtableptr[i] = ptr;
      ptr2 = table[perminv[i]];
      ptr2e = table[perminv[i]+1];
      while (ptr2<ptr2e){
        *(ptr++) = *(ptr2++);
        *(ptr++) = perm[*(ptr2++)];
      }
    }
    newtableptr[ns+1] = ptr;
    tfree(perminv);
    tfree(fsaptr->table->table_data_ptr[0]);
    tfree(fsaptr->table->table_data_ptr);
    fsaptr->table->table_data_ptr = newtableptr;
    fsaptr->table->table_data_ptr[0] = newtable;
  }
}

void
fsa_clear_rws(fsaptr)
/* If *fsaptr is an rws (rewriting-system) automaton, then it may have
 * negative targets, which point to reduction equations.
 * This function simply replaces all of these tagets by 0, to make it
 * into a conventional fsa.
 */
	fsa *fsaptr;
{ int ns, ne, **table, *ptr, *ptr2, *ptr2e, i, j;

  if (fsaptr->table->table_type==SPARSE && fsaptr->table->denserows>0) {
    fprintf(stderr,
    "Sorry: fsa_clear_rws not available for sparse storage with dense rows.\n");
    exit(1);
  }
  if (fsaptr->flags[RWS]==FALSE)
    return;

  ne = fsaptr->alphabet->size;
  ns = fsaptr->states->size;

  table = fsaptr->table->table_data_ptr;
  if (fsaptr->table->table_type == DENSE) {
    for (j=1;j<=ne;j++) for (i=1;i<=ns;i++)
      if (table[j][i]<0)
        table[j][i] = 0;
  }
  else {
    ptr = table[1];
    for (i=1;i<=ns;i++) {
      ptr2 = table[i]; ptr2e = table[i+1];
      table[i] = ptr;
      while (ptr2 < ptr2e) {
        if (*(ptr2+1) > 0) {
          *(ptr++) = *(ptr2++); *(ptr++) = *(ptr2++);
        }
        else
          ptr2 += 2;
      }
    }
    table[ns+1] = ptr;
  }

  fsaptr->flags[RWS]=FALSE;
}

void
fsa_make_accessible(fsaptr)
/* Make the fsa *fsaptr accessible - i.e. all states reachable from start-state
 */
        fsa     *fsaptr;
{ int ns, ct, ne, *gotlist, **table, i, j, k, l, *ptr, *ptre; 
  boolean *got;
  storage_type st = fsaptr->table->table_type;
  fsa *temp;

  if (st==SPARSE && fsaptr->table->denserows>0) {
    fprintf(stderr,
"Sorry: fsa_make_accessible unavailable for sparse storage with dense rows.\n");
    exit(1);
  }
  if (print_level>=3)
    printf("    #Calling fsa_make_accessible.\n");
  if (fsaptr->flags[MINIMIZED] || fsaptr->flags[ACCESSIBLE]
                               || fsaptr->flags[TRIM])
    return;

  if (fsaptr->flags[RWS])
    fsa_clear_rws(fsaptr);

  ne = fsaptr->alphabet->size;
  ns = fsaptr->states->size;

  if (fsaptr->num_initial==ns){
    fsaptr->flags[ACCESSIBLE] = TRUE;
    return;
  }

/* In the deterministic case, we do it by and-ing it with itself - should write
 * separate code for this really, but it hardly seems worth it.
 * This works faster than deleting redundant states.
 */
  if (fsaptr->flags[DFA]) {
    temp = fsa_and(fsaptr,fsaptr,st,TRUE,tmpnam(0));
    fsa_copy(fsaptr,temp);
    fsa_clear(temp);
    tfree(temp);
    return;
  }
    
/* Now find the list of states accessible from the start-state. */
  tmalloc(gotlist,int,ns+1);
  tmalloc(got,boolean,ns+1);
  table = fsaptr->table->table_data_ptr;
  ct = 0;
  for (i=1;i<=ns;i++)
    got[i] = FALSE;
  for (i=1;i<=fsaptr->num_initial;i++) {
    got[fsaptr->initial[i]]=TRUE;
    gotlist[++ct] = fsaptr->initial[i];
  }
  for (i=1;i<=ct;i++) {
    k = gotlist[i];
    if (fsaptr->table->table_type == DENSE)
      for (j=1;j<=ne;j++){
        l = table[j][k];
        if (l>0 && !got[l]) {
          got[l]=TRUE; gotlist[++ct]=l;
        }
      }
    else {
      ptr=table[k]; ptre=table[k+1];
      while (ptr<ptre) {
        l = *(++ptr);
        if (l>0 && !got[l]) {
          got[l]=TRUE; gotlist[++ct]=l;
        }
        ptr++;
      }
    }
  }
  tfree(gotlist);

/* Now delete inaccessible states */
  for (i=ns;i>=1;i--) if (!got[i])
    fsa_delete_state(fsaptr,i);

  tfree(got);
}

void
fsa_minimize(fsaptr)
	fsa	*fsaptr;
/* Minimize the fsa *fsaptr. */
{ int *block_numa, *block_numb, *block_swap, i, j, k, l, len,
       *ptr, *ptr2, *ptr2e, *ht_ptr,
       ne, ns_orig, **table, ns_final, ns_new, num_iterations;
  hash_table ht;
  boolean fixed;

  if (fsaptr->table->table_type==SPARSE && fsaptr->table->denserows>0) {
    fprintf(stderr,
"Sorry: fsa_minimize unavailable for sparse storage with dense rows.\n");
    exit(1);
  }

  if (print_level>=3)
    printf("    #Calling fsa_minimize.\n");
  if (!fsaptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_minimize only applies to DFA's.\n");
    exit(1);
  }
  if (fsaptr->flags[MINIMIZED])
    return;

  if (fsaptr->flags[RWS])
    fsa_clear_rws(fsaptr);

  fsa_make_accessible(fsaptr);

  ns_orig = fsaptr->states->size;
  if (ns_orig==0) {
    fsaptr->flags[TRIM] = TRUE;
    fsaptr->flags[MINIMIZED] = TRUE;
    return;
  }

  /* First throw away any existing structure on the state-set. */
  srec_clear(fsaptr->states);
  fsaptr->states->type = SIMPLE;
  ne = fsaptr->alphabet->size;
  table = fsaptr->table->table_data_ptr;

  tmalloc(block_numa,int,ns_orig+1); tmalloc(block_numb,int,ns_orig+1);
/* Start with block_numa equal to the accept/reject division
 * Remember that state/block number 0 is always failure with no hope.
 */
  if (fsaptr->num_accepting == ns_orig) {
    block_numa[0] = 0;
    for (i=1;i<=ns_orig;i++)
      block_numa[i] = 1;
  }
  else {
    for (i=0;i<=ns_orig;i++)
      block_numa[i] = 0;
    for (i=1;i<=fsaptr->num_accepting;i++)
      block_numa[fsaptr->accepting[i]] = 1;
  }

  fixed = fsaptr->table->table_type==DENSE;

  ns_new = 1;
  num_iterations = 0;
/* The main refinement loop follows. */
  do {
    num_iterations++;
    ns_final = ns_new;
/* Turn off excessive printing at this point */
    j=print_level; print_level=1;
    hash_init(&ht,fixed,ne+1,0,0);
    print_level=j;
    if (print_level>=3)
      printf("    #Iterating - number of state categories = %d.\n",ns_new);
    block_numb[0] = 0;
    for (i=1;i<=ns_orig;i++){
  /* Insert the encoded form of the transitions from state i into the hashtable
   * preceded by the current block number of i.
   */
      len = fixed ? ne+1 : table[i+1]-table[i]+1;
      ptr = ht.current_ptr;
      *ptr = block_numa[i];
      if (fixed) {
        for (j=1;j<len;j++)
          ptr[j] = block_numa[table[j][i]];
          l = len;
      }
      else {
        l = 0;
        for (j=1;j<len;j+=2) {
          k = block_numa[table[i][j]];
          if (k > 0){
            ptr[++l] = table[i][j-1];
            ptr[++l] = k;
          }
        }
        if (l>0 || *ptr>0)
          l++;
/* For technical reasons, we want the zero record to be empty */
      }
      block_numb[i] = hash_locate(&ht,l);
    }
    ns_new = ht.num_recs;
    block_swap=block_numa; block_numa=block_numb; block_numb=block_swap;
    if (ns_new > ns_final)
      hash_clear(&ht);
  } while (ns_new > ns_final);
  
/* At this stage, either ns_final = ns_new, or the fsa has empty accepted
 * language, ns_new=0 and ns_final=1.
 */

  fsaptr->flags[TRIM] = TRUE;
  fsaptr->flags[MINIMIZED] = TRUE;

  if (ns_new==0) {
/* This is the awkward case of no states - always causes problems! */
    fsaptr->states->size=0;
    fsaptr->num_initial=0;
    tfree(fsaptr->initial);
    fsaptr->num_accepting = 0;
    tfree(fsaptr->accepting);
    tfree(fsaptr->table->table_data_ptr[0]);
    tfree(fsaptr->table->table_data_ptr);
  }
  else if (ns_final<ns_orig) {
/* Re-define the fsa fields  */
    fsaptr->states->size = ns_final;

    fsaptr->initial[1] = block_numa[fsaptr->initial[1]];

    if (fsaptr->num_accepting == ns_orig) {
      fsaptr->num_accepting = ns_final;
      if (ns_final==1) {
        tmalloc(fsaptr->accepting,int,2);
        fsaptr->accepting[1] = 1;
      }
    }
    else {
      tmalloc(fsaptr->is_accepting,boolean,ns_final+1);
      for (i=1;i<=ns_final;i++)
        fsaptr->is_accepting[i] = FALSE;
      for (i=1;i<=fsaptr->num_accepting;i++)
        fsaptr->is_accepting[block_numa[fsaptr->accepting[i]]] = TRUE;
      fsaptr->num_accepting = 0;
      for (i=1;i<=ns_final;i++) if (fsaptr->is_accepting[i])
        fsaptr->num_accepting++;
      tfree(fsaptr->accepting);
      tmalloc(fsaptr->accepting,int,fsaptr->num_accepting+1);
      j = 0;
      for (i=1;i<=ns_final;i++) if (fsaptr->is_accepting[i])
        fsaptr->accepting[++j] = i;
      tfree(fsaptr->is_accepting);
    }

/* Finally copy the transition table data from the hash-table back to the fsa */
    tfree(fsaptr->table->table_data_ptr[0]);
    tfree(fsaptr->table->table_data_ptr);
    if (fixed){
      fsa_table_init(fsaptr->table,ns_final,ne);
      table = fsaptr->table->table_data_ptr;
      for (i=1;i<=ns_final;i++) {
        ht_ptr = hash_rec(&ht,i);       
        for (j=1;j<=ne;j++) 
          table[j][i] = ht_ptr[j];
      }
    }
    else{
      tmalloc(fsaptr->table->table_data_ptr,int *,ns_final+2);
      tmalloc(fsaptr->table->table_data_ptr[0],int,ht.tot_space-ns_final);
      table = fsaptr->table->table_data_ptr;
      table[1] = ptr = table[0];
      for (i=1;i<=ns_final;i++){
        ht_ptr = hash_rec(&ht,i);
        ptr2 = ht_ptr+1;
        ptr2e = ht_ptr + hash_rec_len(&ht,i) - 1;
        while (ptr2 <= ptr2e)
          *(ptr++) = *(ptr2++);
        table[i+1] = ptr;
      }
    }
  }
  hash_clear(&ht);
  tfree(block_numa);
  tfree(block_numb);
  if (print_level>=3)
    printf("    #Number of iterations = %d.\n",num_iterations);
}

void
fsa_labeled_minimize(fsaptr)
	fsa	*fsaptr;
/* This is the minimization function for fsa's which misght have more than
 * two categories of states. 
 * We use the labeled set-record type to identify the categories, so *fsaptr
 * must have state-set of labeled type.
 */
{ int *block_numa, *block_numb, *block_swap, i, j, k, l, len,
       *ptr, *ptr2, *ptr2e, *ht_ptr,
       ne, ns_orig, **table, ns_final, ns_new, num_iterations;
  hash_table ht;
  boolean fixed, *occurs;

  if (fsaptr->table->table_type==SPARSE && fsaptr->table->denserows>0) {
    fprintf(stderr,
"Sorry: fsa_minimize unavailable for sparse storage with dense rows.\n");
    exit(1);
  }
  if (print_level>=3)
    printf("    #Calling fsa_labeled_minimize.\n");
  if (!fsaptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_labeled_minimize only applies to DFA's.\n");
    exit(1);
  }

  if (fsaptr->states->type != LABELED) {
    fprintf(stderr,
        "Error: in fsa_labeled_minimize state=set must have type labeled.\n");
    exit(1);
  }

  if (fsaptr->flags[RWS])
    fsa_clear_rws(fsaptr);

  fsa_make_accessible(fsaptr);

  ne = fsaptr->alphabet->size;
  table = fsaptr->table->table_data_ptr;
  ns_orig = fsaptr->states->size;
  if (ns_orig==0)
    return;

/* Start with block_numa equal to the subdivision defined by the labeling.  */
  tmalloc(occurs,boolean,fsaptr->states->labels->size+1);
  for (i=1;i<=fsaptr->states->labels->size;i++)
    occurs[i] = FALSE;
  tmalloc(block_numa,int,ns_orig+1);
  tmalloc(block_numb,int,ns_orig+1);
  ns_new = 0;
  block_numa[0] = 0; /* The zero state is always regarded as having label 0 */
  for (i=1;i<=ns_orig;i++) {
    j = fsaptr->states->setToLabels[i];
    if (j>0 && !occurs[j]) {
      occurs[j] = TRUE;
      ns_new++;
    }
    block_numa[i] = j;
  }
  tfree(occurs);
  if (ns_new==0) ns_new=1; /* a patch for the case when there are no labels */

  fixed = fsaptr->table->table_type==DENSE;

  num_iterations = 0;
/* The main refinement loop follows. */
  do {
    num_iterations++;
    ns_final = ns_new;
/* Turn off excessive printing at this point */
    j=print_level; print_level=1;
    hash_init(&ht,fixed,ne+1,0,0);
    print_level=j;
    if (print_level>=3)
      printf("    #Iterating - number of state categories = %d.\n",ns_new);
    block_numb[0] = 0;
    for (i=1;i<=ns_orig;i++){
  /* Insert the encoded form of the transitions from state i into the hashtable
   * preceded by the current block number of i.
   */
      len = fixed ? ne+1 : table[i+1]-table[i]+1;
      ptr = ht.current_ptr;
      *ptr = block_numa[i];
      if (fixed) {
        for (j=1;j<len;j++)
          ptr[j] = block_numa[table[j][i]];
        l = len;
      }
      else {
        l = 0;
        for (j=1;j<len;j+=2) {
          k = block_numa[table[i][j]];
          if (k > 0){
            ptr[++l] = table[i][j-1];
            ptr[++l] = k;
          }
        }
        if (l>0 || *ptr>0)
          l++;
/* For technical reasons, we want the zero record to be empty */
      }
      block_numb[i] = hash_locate(&ht,l);
    }
    ns_new = ht.num_recs;
    block_swap=block_numa; block_numa=block_numb; block_numb=block_swap;
    if (ns_new > ns_final)
      hash_clear(&ht);
  } while (ns_new > ns_final);
  
/* At this stage, either ns_final = ns_new, or the fsa has empty accepted
 * language, ns_new=0 and ns_final=1.
 */
  
  if (ns_new==0) {
/* This is the awkward case of no states - always causes problems! */
    fsaptr->states->size=0;
    fsaptr->num_initial=0;
    tfree(fsaptr->initial);
    fsaptr->num_accepting = 0;
    tfree(fsaptr->accepting);
    tfree(fsaptr->table->table_data_ptr[0]);
    tfree(fsaptr->table->table_data_ptr);
  }
  else if (ns_final<ns_orig){
/* Re-define the fsa fields  */
    fsaptr->states->size = ns_final;

    fsaptr->initial[1] = block_numa[fsaptr->initial[1]];

    for (i=1;i<=ns_orig;i++)
      block_numb[block_numa[i]] = fsaptr->states->setToLabels[i];
    for (i=1;i<=ns_final;i++)
      fsaptr->states->setToLabels[i] = block_numb[i];

/* Finally copy the transition table data from the hash-table back to the fsa */
    tfree(fsaptr->table->table_data_ptr[0]);
    tfree(fsaptr->table->table_data_ptr);
    if (fixed){
      fsa_table_init(fsaptr->table,ns_final,ne);
      table = fsaptr->table->table_data_ptr;
      for (i=1;i<=ns_final;i++) {
        ht_ptr = hash_rec(&ht,i);
        for (j=1;j<=ne;j++)
          table[j][i] = ht_ptr[j];
      }
    }
    else{
      tfree(fsaptr->table->table_data_ptr);
      tmalloc(fsaptr->table->table_data_ptr,int *,ns_final+2);
      tmalloc(fsaptr->table->table_data_ptr[0],int,ht.tot_space-ns_final);
      table = fsaptr->table->table_data_ptr;
      table[1] = ptr = table[0];
      for (i=1;i<=ns_final;i++){
        ht_ptr = hash_rec(&ht,i);
        ptr2 = ht_ptr+1;
        ptr2e = ht_ptr + hash_rec_len(&ht,i) - 1;
        while (ptr2 <= ptr2e)
          *(ptr++) = *(ptr2++);
        table[i+1] = ptr;
      }
    }
  }
  hash_clear(&ht);
  tfree(block_numa);
  tfree(block_numb);
  if (print_level>=3)
    printf("    #Number of iterations = %d.\n",num_iterations);
}


void
fsa_bfs(fsaptr)
	fsa	*fsaptr;
/* Put the fsa *fsapt into bfs form. */
{  int ns, ne, *perm, *perminv, **table, ct, i, j, s, t;
   boolean *got, dense;

  if (fsaptr->table->table_type==SPARSE && fsaptr->table->denserows>0) {
    fprintf(stderr,
"Sorry: fsa_bfs unavailable for sparse storage with dense rows.\n");
    exit(1);
  }
  if (print_level>=3)
    printf("    #Calling fsa_bfs.\n");
  if (!fsaptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_bfs only applies to DFA's.\n");
    exit(1);
  }
  if (fsaptr->flags[BFS])
    return;

  fsa_make_accessible(fsaptr);

  if (fsaptr->flags[RWS])
    fsa_clear_rws(fsaptr);

  ne = fsaptr->alphabet->size;
  ns = fsaptr->states->size;
  if (ns==0)
    return;
 
  dense = fsaptr->table->table_type == DENSE;
  table = fsaptr->table->table_data_ptr;

  tmalloc(got,boolean,ns+1);
  for (i=1;i<=ns;i++)
    got[i] = FALSE;
  tmalloc(perm,int,ns+1);
  perm[0] = 0;
  perm[1] = fsaptr->initial[1];
  got[perm[1]] = TRUE;
  ct = 1;
  i = 1;
  while (i<=ct) {
    s = perm[i];
    for (j=1;j<=ne;j++) {
      t = target(dense,table,j,s,0);
      if (t>0 && !got[t]) {
        perm[++ct] = t;
        got[t] = TRUE;
      }
    }
    i++;
  }

  tfree(got);
/* The inverse of perm is what is required */
  tmalloc(perminv,int,ns+1);
  for (i=0;i<=ns;i++)
    perminv[perm[i]] = i;
  tfree(perm);
  fsa_permute_states(fsaptr,perminv);
  tfree(perminv);
  fsaptr->flags[BFS]=TRUE;
}

int
fsa_count(fsaptr)
	fsa *fsaptr;
/* Count the size of the accepted language of the fsa *fsaptr.
 * Return this size, or -1 if infinite.
 */
{ int i, j, ct, ne, ns, **table, dr,
      *in_degree, *ordered_state_list, *num_acc_words, cstate, im;
  boolean dense;

  if (!fsaptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_count only applies to DFA's.\n");
    exit(1);
  }
  if (!fsaptr->flags[TRIM]){
    printf("#WARNING: fsa_count applied to fsa not known to be trim.\n");
    printf("#It will be minimized before proceeding.\n");
    fsa_minimize(fsaptr);
  }

  ne = fsaptr->alphabet->size;
  ns = fsaptr->states->size;
  if (ns==0)
    return 0;

  dr = fsaptr->table->denserows;

/* We first count the number of edges going into each state */
  tmalloc(in_degree,int,ns+1);
  for (i=1;i<=ns;i++)
    in_degree[i] = 0;
  dense = fsaptr->table->table_type==DENSE;
  table = fsaptr->table->table_data_ptr;
  for (i=1;i<=ns;i++) for (j=1;j<=ne;j++) {
    im = target(dense,table,j,i,dr);
    if (im>0)
      in_degree[im]++;
  }

/* Now we try to order the states so that if state s <= state t, then there
 * is no path from state t to state s. If this is not possible, then the
 * accepted language must be infinite.
 */
  cstate = fsaptr->initial[1];
  if (in_degree[cstate] > 0) {
    tfree(in_degree);
    return -1;
  }
  tmalloc(ordered_state_list,int,ns+1);
  ordered_state_list[1] = cstate;

  ct = 1;
  i = 0;
  while (++i<=ct) {
    cstate = ordered_state_list[i];
    for (j=1;j<=ne;j++) {
      im = target(dense,table,j,cstate,dr);
      if (im>0) {
        in_degree[im]--;
        if (in_degree[im]==0) 
          ordered_state_list[++ct] = im;
      }
    }
  }

  if (ct!=ns) {
    tfree(in_degree) tfree(ordered_state_list)
    return -1;
  }

/* We have built the list, so the accepted language is finite. Now we work
 * backwards through the list, calculating the number of accepted words
 * starting at that state.
 * We might as well use the same space as was used by in_degree, which
 * is no longer needed.
 */
  fsa_set_is_accepting(fsaptr);
  num_acc_words = in_degree;
  for (i=ns;i>=1;i--) {
    cstate = ordered_state_list[i];
    num_acc_words[cstate] = 0;
    for (j=1;j<=ne;j++) {
      im = target(dense,table,j,cstate,dr);
      if (im>0)
        num_acc_words[cstate] += num_acc_words[im];
    }
    if (fsaptr->is_accepting[cstate])
      num_acc_words[cstate]++;
  }

  ct = num_acc_words[fsaptr->initial[1]];
  tfree(in_degree) tfree(ordered_state_list) tfree(fsaptr->is_accepting)
  return ct;
}

boolean
fsa_enumerate(wfile,fsaptr,min,max,putcomma)
	FILE *wfile;
	fsa *fsaptr;
	int min, max;
	boolean putcomma;
/* Enumerate the subset of the language of the finite state automaton *fsaptr,
 * consisting of those words having length l satisfying min <= l <= max.
 * Since there is no point in storing these words currently, they are
 * printed out to the file wfile, with all but the last one to be printed
 * followed by a comma and carriage return.
 * TRUE is returned if some words are found - otherwise FALSE.
 * If putcomma is true, then the first word to be printed is preceded by comma
 * and carriage-return (to handle the case when this function is called
 * repeatedly).
 */
{ int i, g1, g2, ne, ngens, ngens1, ns, **table, dr,
      *cword, firste, clength, clength1, clength2, cstate, im, *statelist;
  boolean dense, done, backtrack, foundword;
  char *cword1, *cword2, prod, numbers, **names;

  if (!fsaptr->flags[DFA]){
    fprintf(stderr,"Error: fsa_enumerate only applies to DFA's.\n");
    exit(1);
  }

  if (fsaptr->num_initial==0)
    return FALSE;

  ne = fsaptr->alphabet->size;
  ns = fsaptr->states->size;

  dr = fsaptr->table->denserows;

  dense = fsaptr->table->table_type==DENSE;
  table = fsaptr->table->table_data_ptr;

  tmalloc(cword,int,max+1); /* to hold the current word in the enumeration */

/* "names" will be used to store the symbols used for printing. If no natural
 * names are available, then we just use the numbers of the edges.
 */
  prod = fsaptr->alphabet->type==PRODUCT;
  if (prod) {
    if (fsaptr->alphabet->arity!=2){
      fprintf(stderr,
	"Sorry  can only cope with two-variable product alphabets.\n");
      return FALSE;
    }
    ngens = fsaptr->alphabet->base->size;
    ngens1 = ngens+1;
    numbers = fsaptr->alphabet->base->type==SIMPLE ||
              fsaptr->alphabet->base->type==LABELED ||
              fsaptr->alphabet->base->type==WORDS;
    if (!numbers)
      names = fsaptr->alphabet->base->names;
    tmalloc(cword1,char,max+1); tmalloc(cword2,char,max+1);
   /* These are used for the two components of the word cword = [cword1,cword2].
    */
  }
  else {
    numbers = fsaptr->alphabet->type==SIMPLE ||
	      fsaptr->alphabet->type==LABELED ||
	      fsaptr->alphabet->type==WORDS;
    if (!numbers)
        names = fsaptr->alphabet->names;
  }
  if (numbers) { /* define the integral names. */
    tmalloc(names,char *,ne+1);
    for (i=1;i<=ne;i++) {
      tmalloc(names[i],char,int_len(i)+1);
      sprintf(names[i],"%d",i);
    }
  }

  *buffer = '\0';
  tmalloc(statelist,int,max+1);
 /* this is used to store the state history on scanning the current word. */
  for (i=0;i<=max;i++)
    cword[i] = 0;
  clength = 0;
  statelist[0] = fsaptr->initial[1];
  fsa_set_is_accepting(fsaptr);
  done = FALSE;
  foundword=FALSE;
/* Backtrack search can now begin */
  while (!done) {
  /* first see if we want the current word - if so, print it */
    if (clength>=min && fsaptr->is_accepting[statelist[clength]]) {
      if (putcomma || foundword)
        fprintf(wfile,",\n");
      foundword = TRUE;
      if (prod) { /* split word into components and print them */
        clength1 = clength2 = 0;
        for (i=0;i<clength;i++) {
          g1 = (cword[i]-1)/ngens1 + 1; g2 = (cword[i]-1)%ngens1 + 1;
          if (g1<ngens1) cword1[clength1++] = g1;
          if (g2<ngens1) cword2[clength2++] = g2;
        }
        cword1[clength1] = '\0'; cword2[clength2] = '\0';
        add_to_buffer(0," [");
        add_word_to_buffer(wfile,cword1,names);
        add_to_buffer(0,",");
        add_word_to_buffer(wfile,cword2,names);
        add_to_buffer(0,"]");
      }
      else {
        add_to_buffer(0," ");
        add_iword_to_buffer(wfile,cword,names);
      }
      fprintf(wfile,"%s",buffer);
      *buffer = '\0';
    }
  /* Now proceed to next word in the search */
    firste = 1;
    backtrack = TRUE;
    while (backtrack && !done) {
      if (clength<max) {
        cstate = statelist[clength];
        i = firste-1;
        while (backtrack && ++i<=ne) 
          if ((im =target(dense,table,i,cstate,dr)) > 0) { /* found next node */
            cword[clength++] = i;
            statelist[clength] = im;
            backtrack = FALSE;
          }
      }
      if (backtrack) {
        if (clength==0)
          done = TRUE;
        else {
          firste = cword[--clength]+1;
          cword[clength] = '\0';
        }
      }
    }
  }
  
  if (numbers) {
    for (i=1;i<=ne;i++) tfree(names[i]);
    tfree(names);
  }
  tfree(fsaptr->is_accepting);
  tfree(cword);
  if (prod) {
    tfree(cword1); tfree(cword2);
  }
  tfree(statelist);

  return foundword;
}
