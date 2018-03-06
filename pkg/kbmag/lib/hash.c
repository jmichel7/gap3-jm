/* file hash.c
 * 23/12/94 - changes to avoid copying whole records - now records are
 * saved in a series of independent blocks, which are allocated when
 * needed. pointers are still copied when number needs increasing.
 * 24/10/94.
 * This files contains functions for dealing with hash tables.
 * If ht is a hash-table, then new records should be inserted into the table
 * at location beginning at ht.current_ptr, and located by the function
 * hash_locate (or short_hash_locate).
 * If the record is new it remains where it is.
 */
#define MAXBLOCKS	1024
#define MAXRECSINC	16384
#define SPACEINC	131072
#define MAXRECSINCL	262144
#define SPACEINCL	1048576
#define MAXRECSINCH	2097152
#define SPACEINCH	16777216
#define HASH_MOD	2039
#define HASH_VAL	65521
#define HASH_VALL	524047 /* last three constants are all prime */
#include "defs.h"
#include "hash.h"
#include "externals.h"

#define HTMARGIN	4096
/* if less space than this in hash-table, re-allocate */

/* The following functions are defined in this file and visible externally */
void 	hash_init();
void	short_hash_init();
void	hash_clear();
void	short_hash_clear();
int	hash_rec_len();
int	short_hash_rec_len();
int	*hash_rec();
unsigned short	*short_hash_rec();
int	hash_locate();
int	short_hash_locate();

/* The following functions should be used only within this file */
void hash_morerecs();
void short_hash_morerecs();
void hash_morespace();
void short_hash_morespace();

void
hash_init(htptr,fixed,len,num_recs_inc,space_inc)
	hash_table *htptr;
	boolean fixed;
	int len;
	int num_recs_inc, space_inc;
/* Note: len = record-length is ignored if fixed is FALSE
 * If num_recs_inc and/or space_inc is zero, they are given default values,
 * depending on whether large or huge is true.
 */
{ int i;

  if (print_level>=3) {
    printf("    #Calling hash_init.");
    if (fixed) printf(" Records of fixed length %d.\n",len);
    else printf(" Records of variable length.\n");
  }
  if (num_recs_inc==0) 
    htptr->num_recs_inc = huge ? MAXRECSINCH : large ? MAXRECSINCL : MAXRECSINC;
  else
    htptr->num_recs_inc = num_recs_inc;
  if (space_inc==0) 
    htptr->space_inc = huge ? SPACEINCH/2 : large ? SPACEINCL/2 : SPACEINC;
  else
    htptr->space_inc = space_inc;
  htptr->fixed_len = fixed;
  htptr->num_recs = 0;
  htptr->maxrecs = htptr->num_recs_inc; 
  tmalloc(htptr->table_block,int *,MAXBLOCKS);
  htptr->num_blocks = 1;
  tmalloc(htptr->table_block[0],int,htptr->space_inc);
  
  if (htptr->fixed_len) {
    htptr->rec_len = len;
    htptr->table_data_ptr = 0; /* we won't use it for fixed length records */
    for (i=0;i<len;i++)
      htptr->table_block[0][i] = 0; /* This is record number 0 */
    htptr->block_space = htptr->tot_space = len;
    htptr->current_ptr = htptr->table_block[0]+len;
    htptr->recs_per_block = htptr->space_inc/len;
  }
  else{
    tmalloc(htptr->table_data_ptr,int *,htptr->maxrecs);
    htptr->table_data_ptr[0] = htptr->table_data_ptr[1] = htptr->table_block[0];
    htptr->block_space = htptr->tot_space = 0;
    htptr->current_ptr = htptr->table_block[0];
    tmalloc(htptr->block_start_rec,int,MAXBLOCKS);
    tmalloc(htptr->block_last_len,int,MAXBLOCKS);
  }

  htptr->modulus = HASH_MOD;
  htptr->hash_values = huge || large ? HASH_VALL : HASH_VAL;

  tmalloc(htptr->first_rec,int,htptr->hash_values);
  htptr->first_rec[0] = 0;
  for (i=1;i<htptr->hash_values;i++)
    htptr->first_rec[i] = -1;
  tmalloc(htptr->next_rec,int,htptr->maxrecs);
  htptr->next_rec[0] = -1;
}

void
short_hash_init(htptr,fixed,len,num_recs_inc,space_inc)
	short_hash_table *htptr;
	boolean fixed;
	int len;
	int num_recs_inc, space_inc;
/* Note: len = record-length is ignored if fixed is FALSE
 * If num_recs_inc and/or space_inc is zero, they are given default values,
 * depending on whether large or huge is true.
 */
{ int i;

  if (print_level>=3) {
    printf("    #Calling short_hash_init.");
    if (fixed) printf(" Records of fixed length %d.\n",len);
    else printf(" Records of variable length.\n");
  }
  if (num_recs_inc==0) 
    htptr->num_recs_inc = huge ? MAXRECSINCH : large ? MAXRECSINCL : MAXRECSINC;
  else
    htptr->num_recs_inc = num_recs_inc;
  if (space_inc==0) 
    htptr->space_inc = huge ? SPACEINCH : large ? SPACEINCL : SPACEINC;
  else
    htptr->space_inc = space_inc;
  htptr->fixed_len = fixed;
  htptr->num_recs = 0;
  htptr->maxrecs = htptr->num_recs_inc; 
  tmalloc(htptr->table_block,unsigned short *,MAXBLOCKS);
  htptr->num_blocks = 1;
  tmalloc(htptr->table_block[0],unsigned short,htptr->space_inc);

  if (htptr->fixed_len) {
    htptr->rec_len = len;
    htptr->table_data_ptr = 0; /* we won't use it for fixed length records */
    htptr->block_space = htptr->tot_space = len;
    for (i=0;i<len;i++)
      htptr->table_block[0][i] = 0; /* This is record number 0 */
    htptr->current_ptr = htptr->table_block[0]+len;
    htptr->recs_per_block = htptr->space_inc/len;
  }
  else{
    tmalloc(htptr->table_data_ptr,unsigned short *,htptr->maxrecs);
    htptr->table_data_ptr[0] = htptr->table_data_ptr[1] = htptr->table_block[0];
    htptr->current_ptr = htptr->table_block[0];
    htptr->block_space = htptr->tot_space = 0;
    tmalloc(htptr->block_start_rec,int,MAXBLOCKS);
    tmalloc(htptr->block_last_len,int,MAXBLOCKS);
  }

  htptr->modulus = HASH_MOD;
  htptr->hash_values = huge || large ? HASH_VALL : HASH_VAL;

  tmalloc(htptr->first_rec,int,htptr->hash_values);
  htptr->first_rec[0] = 0;
  for (i=1;i<htptr->hash_values;i++)
    htptr->first_rec[i] = -1;
  tmalloc(htptr->next_rec,int,htptr->maxrecs);
  htptr->next_rec[0] = -1;
}

void
hash_clear(htptr)
	hash_table *htptr;
{ int i;
  for (i=0;i<htptr->num_blocks;i++)
    tfree(htptr->table_block[i]);
  tfree(htptr->table_block);
  tfree(htptr->table_data_ptr);
  if (!htptr->fixed_len) {
    tfree(htptr->block_start_rec);
    tfree(htptr->block_last_len);
  }
  tfree(htptr->first_rec);
  tfree(htptr->next_rec);
}

void
short_hash_clear(htptr)
	short_hash_table *htptr;
{ int i;
  for (i=0;i<htptr->num_blocks;i++)
    tfree(htptr->table_block[i]);
  tfree(htptr->table_block);
  tfree(htptr->table_data_ptr);
  if (!htptr->fixed_len) {
    tfree(htptr->block_start_rec);
    tfree(htptr->block_last_len);
  }
  tfree(htptr->first_rec);
  tfree(htptr->next_rec);
}

int
hash_rec_len(htptr,n)
	hash_table  *htptr;
	int n;
/* The length of hash-record number n. Only useful for variable-length tables */
{ int *ptr, *ptre, *bn, ct;
  if (htptr->fixed_len)
    return htptr->rec_len;
  else {
    ptre = htptr->table_data_ptr[n+1];
    if (ptre==0) {
      bn = htptr->block_start_rec;
      ct = 1;
      while (bn[ct] != n+1)
        ct++;
      return htptr->block_last_len[ct-1];
    }
    else  {
      ptr = htptr->table_data_ptr[n];
      if (ptr==0) {
        bn = htptr->block_start_rec;
        ct = 1;
        while (bn[ct] != n)
          ct++;
        return ptre - htptr->table_block[ct];
      }
      else
        return ptre - ptr;
    }
  }
}

int
short_hash_rec_len(htptr,n)
	short_hash_table  *htptr;
	int n;
/* The length of hash-record number n. Only useful for variable-length tables */
{ unsigned short *ptr, *ptre;
  int  *bn, ct;
  if (htptr->fixed_len)
    return htptr->rec_len;
  else {
    ptre = htptr->table_data_ptr[n+1];
    if (ptre==0) {
      bn = htptr->block_start_rec;
      ct = 1;
      while (bn[ct] != n+1)
        ct++;
      return htptr->block_last_len[ct-1];
    }
    else  {
      ptr = htptr->table_data_ptr[n];
      if (ptr==0) {
        bn = htptr->block_start_rec;
        ct = 1;
        while (bn[ct] != n)
          ct++;
        return ptre - htptr->table_block[ct];
      }
      else
        return ptre - ptr;
    }
  }
}

int *
hash_rec(htptr,n)
	hash_table  *htptr;
        int n;
/* Pointer to record number n */
{ int *ptr, *bn, ct;
  if (htptr->fixed_len)
    return htptr->table_block[n/htptr->recs_per_block] +
           (n % htptr->recs_per_block)*htptr->rec_len;
  else {
    ptr = htptr->table_data_ptr[n];
    if (ptr)
      return ptr;
    else {
      bn = htptr->block_start_rec;
      ct = 1;
      while (bn[ct] != n)
        ct++;
      return htptr->table_block[ct];
    }
  }
}

unsigned short *
short_hash_rec(htptr,n)
	short_hash_table  *htptr;
        int n;
/* Pointer to record number n */
{ unsigned short *ptr;
  int *bn, ct;
  if (htptr->fixed_len)
    return htptr->table_block[n/htptr->recs_per_block] +
           (n % htptr->recs_per_block)*htptr->rec_len;
  else {
    ptr = htptr->table_data_ptr[n];
    if (ptr)
      return ptr;
    else {
      bn = htptr->block_start_rec;
      ct = 1;
      while (bn[ct] != n)
        ct++;
      return htptr->table_block[ct];
    }
  }
}

int
hash_locate(htptr,reclen)
	hash_table *htptr;
	int reclen;
/* This is the basic search function, using the hash-table.
 * It is assumed that the entry for which we are searching is already
 * in the table, and pointed at by htptr->current_ptr
 * and has length reclen.
 * If the item is not found in the existing part of the table,
 * then all we have to do is increment num_recs to add the item to the table.
 * (and relevant first_rec/next_rec pointers).
 * We then allocate more space if necessary.
 * In any case, we return the number of the record.
 */
{ int nr, *rec, *cand, candlen, *nextptr, candno, hashval, coeff,
      k, m, hv, i;
  boolean found, fixed;

  nr = htptr->num_recs;
  rec = htptr->current_ptr;
  fixed = htptr->fixed_len;

/* Now calculate the hashed value of the record. */
  m = htptr->modulus; hv = htptr->hash_values;
  hashval = 0; coeff = 1;
  for (i=0;i<reclen;i++){
    k = rec[i] % m;
    hashval += (coeff*k); hashval %= hv;
    coeff = coeff<<1; coeff %= hv;
  }

  nextptr =  htptr->first_rec + hashval; 
  candno = *nextptr;
  while (candno >= 0) {
    candlen = hash_rec_len(htptr,candno);
    if (candlen == reclen) {
      cand =  hash_rec(htptr,candno);
      found = TRUE;
      for (i=0;i<reclen;i++)
        if (cand[i] != rec[i]) {
           found = FALSE;
           break;
        }
      if (found)
        return(candno);
    }
    nextptr =  htptr->next_rec + candno; 
    candno = *nextptr;
  }
/* The record is now known to be new */
  htptr->current_ptr += reclen;
  htptr->num_recs = ++nr;
  if (!fixed)
    htptr->table_data_ptr[nr+1] = htptr->current_ptr;
  *nextptr = nr;
  htptr->next_rec[nr] = -1;
  htptr->block_space += reclen; htptr->tot_space += reclen;

  if (nr+2 >= htptr->maxrecs)
    hash_morerecs(htptr);
  if ((fixed && htptr->space_inc - htptr->block_space <= reclen) ||
     (!fixed && htptr->space_inc - htptr->block_space <= HTMARGIN) )
    hash_morespace(htptr); 

  return(nr);
}

int
short_hash_locate(htptr,reclen)
	short_hash_table *htptr;
	int reclen;
/* This is the basic search function, using the short_hash-table.
 * It is assumed that the entry for which we are searching is already
 * in the table, and pointed at by htptr->table_data[htptr->num_recs+1]
 * and has length reclen.
 * If the item is not found in the existing part of the table,
 * then all we have to do is increment num_recs to add the item to the table.
 * (and relevant first_rec/next_rec pointers).
 * We then allocate more space if necessary.
 * In any case, we return the number of the record.
 */
{ int nr, candlen, *nextptr, candno, hashval, coeff,
      k, m, hv, i;
  unsigned short *rec, *cand; 
  boolean found, fixed;

  nr = htptr->num_recs;
  rec = htptr->current_ptr;
  fixed = htptr->fixed_len;

/* Now calculate the hashed value of the record. */
  m = htptr->modulus; hv = htptr->hash_values;
  hashval = 0; coeff = 1;
  for (i=0;i<reclen;i++){
    k = rec[i] % m;
    hashval += (coeff*k); hashval %= hv;
    coeff = coeff<<1; coeff %= hv;
  }

  nextptr =  htptr->first_rec + hashval; 
  candno = *nextptr;
  while (candno >= 0) {
    candlen = short_hash_rec_len(htptr,candno);
    if (candlen == reclen) {
      cand =  short_hash_rec(htptr,candno);
      found = TRUE;
      for (i=0;i<reclen;i++)
        if (cand[i] != rec[i]) {
           found = FALSE;
           break;
        }
      if (found)
        return(candno);
    }
    nextptr =  htptr->next_rec + candno; 
    candno = *nextptr;
  }
/* The record is now known to be new */
  htptr->current_ptr += reclen;
  htptr->num_recs = ++nr;
  if (!fixed)
    htptr->table_data_ptr[nr+1] = htptr->current_ptr;
  *nextptr = nr;
  htptr->next_rec[nr] = -1;
  htptr->block_space += reclen; htptr->tot_space += reclen;

  if (nr+2 >= htptr->maxrecs)
    short_hash_morerecs(htptr);
  if ((fixed && htptr->space_inc - htptr->block_space <= reclen) ||
     (!fixed && htptr->space_inc - htptr->block_space <= HTMARGIN) )
    short_hash_morespace(htptr); 

  return(nr);
}

void
hash_morerecs(htptr)
	hash_table  *htptr;
/* Allocate space for more records. */
{ int *new, **newp, *ptr, *ptre, *ptrc, **dptr, **dptre, **dptrc;

  if (print_level>=3)
    printf("    #Calling hash_morerecs. maxrecs increased from %d to %d.\n",
            htptr->maxrecs,htptr->maxrecs + htptr->num_recs_inc);
  htptr->maxrecs += htptr->num_recs_inc;

  if (!htptr->fixed_len){
    tmalloc(newp,int *,htptr->maxrecs);
    dptr = htptr->table_data_ptr;
    dptre = dptr+htptr->num_recs+1;
    dptrc = newp;
    while (dptr <= dptre)
      *(dptrc++) = *(dptr++);
    tfree(htptr->table_data_ptr);
    htptr->table_data_ptr = newp;
  }

  tmalloc(new,int,htptr->maxrecs);
  ptr = htptr->next_rec;
  ptre = ptr+htptr->num_recs;
  ptrc = new;
  while (ptr <= ptre)
    *(ptrc++) = *(ptr++);
  tfree(htptr->next_rec);
  htptr->next_rec = new;
}

void
short_hash_morerecs(htptr)
	short_hash_table  *htptr;
/* Allocate space for more records. */
{ int *new, *ptr, *ptre, *ptrc;
  unsigned short **newp, **dptr, **dptre, **dptrc;

  if (print_level>=3) 
    printf(
        "    #Calling short_hash_morerecs. maxrecs increased from %d to %d.\n",
            htptr->maxrecs,htptr->maxrecs + htptr->num_recs_inc);
  htptr->maxrecs += htptr->num_recs_inc;

  if (!htptr->fixed_len){
    tmalloc(newp,unsigned short *,htptr->maxrecs);
    dptr = htptr->table_data_ptr;
    dptre = dptr+htptr->num_recs+1;
    dptrc = newp;
    while (dptr <= dptre)
      *(dptrc++) = *(dptr++);
    tfree(htptr->table_data_ptr);
    htptr->table_data_ptr = newp;
  }

  tmalloc(new,int,htptr->maxrecs);
  ptr = htptr->next_rec;
  ptre = ptr+htptr->num_recs;
  ptrc = new;
  while (ptr <= ptre)
    *(ptrc++) = *(ptr++);
  tfree(htptr->next_rec);
  htptr->next_rec = new;
}

void
hash_morespace(htptr)
	hash_table *htptr;
/* Allocate more table space */
{ int nb, nr;

  htptr->num_blocks++;
  nb = htptr->num_blocks;
  if (nb > MAXBLOCKS) {
    fprintf(stderr,"#Too much hash-space - try running with -l or -h.\n");
    exit(1);
  }
  if (print_level>=3)
    printf("    #Calling hash_morespace. maxspace increased by %d to %d.\n",
            htptr->space_inc,htptr->space_inc * nb);

  nr =  htptr->num_recs;
  tmalloc(htptr->table_block[nb-1],int,htptr->space_inc);
  htptr->block_space = 0;
  htptr->current_ptr = htptr->table_block[nb-1];
  if (!htptr->fixed_len) {
    htptr->block_start_rec[nb-1] = nr+1; 
    htptr->block_last_len[nb-2] =
                       htptr->table_data_ptr[nr+1] -  htptr->table_data_ptr[nr];
    htptr->table_data_ptr[nr+1] = 0;
  }
}

void
short_hash_morespace(htptr)
	short_hash_table *htptr;
/* Allocate more table space */
{ int  nb, nr;
  htptr->num_blocks++;
  nb = htptr->num_blocks;
  if (nb > MAXBLOCKS) {
    fprintf(stderr,"#Too much hash-space - try running with -l or -h.\n");
    exit(1);
  }
  if (print_level>=3)
    printf(
        "    #Calling short_hash_morespace. maxspace increased by %d to %d.\n",
            htptr->space_inc,htptr->space_inc * nb);

  nr =  htptr->num_recs;
  tmalloc(htptr->table_block[nb-1],unsigned short,htptr->space_inc);
  htptr->block_space = 0;
  htptr->current_ptr = htptr->table_block[nb-1];
  if (!htptr->fixed_len) {
    htptr->block_start_rec[nb-1] = nr+1; 
    htptr->block_last_len[nb-2] =
                       htptr->table_data_ptr[nr+1] -  htptr->table_data_ptr[nr];
    htptr->table_data_ptr[nr+1] = 0;
  }
}
