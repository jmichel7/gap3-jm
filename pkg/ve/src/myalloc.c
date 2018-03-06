
/* Validated malloc and free */
/* These routines are only for debugging purposes, they aim to
   trap the most common forms of heap error */

#define MYALLOC_C
#include <stdio.h>
#include <string.h>
#include "myalloc.h"


#define SMAGIC 0x1234
#define EMAGIC 0x5678
#define EMPTY  (char)0xE5
#define FILLER (char) 0xF6
typedef long align;

typedef struct mablock * maptr;

struct mablock
{
  int smagic; /* magic number to mark start of block obtained from myalloc */
  size_t len; /* length requested from myalloc */
  size_t len1; /* length actually allocated */
  char *file;
  int line;
  maptr next; /* pointers to doubly linked list of allocated blocks */
  maptr last;
};

/* head for the doubly linked list */
static struct mablock mahead = { SMAGIC, 0, 0, __FILE__, __LINE__,
				   &mahead, &mahead};
#define ROUNDING 4 /* minimum unit for allocation (in bytes) */

/* Allocate n bytes. If b is true then return NULL if no more room,
   otherwise generate a fatal error */

pointer myreallocf(p,n,b,file,line)
     pointer p;
     size_t n;
     bool b;
     char *file;
     int line;
{
  char *s;
  pointer p1;
  maptr m;
  if (!p)
    {
      fprintf(stderr,"Reallocating NULL\n");
      abort();
    }
  m = (maptr)((char *)p - sizeof(struct mablock));
  if (n < m->len) /* reducing */
    {
      for (s = (char *)p+n; s < (char *)p + m->len; s++)
	*s = FILLER;
      m->len = n;
      m->file = file;
      m->line =line;
      return(p);
    }
  if (n < m->len1)
    {
      for (s = (char *)p+m->len; s < (char *)p+n; s++)
	*s = EMPTY;
      m->len = n;
      m->file = file;
      m->line = line;
      return(p);
    }
  p1 = myallocf(n,b,file,line);
  if (p1)
    {
      memcpy(p1,p,m->len);
      myfreef(p,file,line);
    }
  return(p1);
}


pointer  myallocf(n,b,file,line)
     size_t n;
     bool b;
     char *file;
     int line;
{
  char * p;
  maptr m;
  int i;
  unsigned int n1;
  n1 = n+ROUNDING-1;		/* round up to a legal amount */
  n1 = n1 - n1%ROUNDING;
  p = malloc(n1+sizeof(struct mablock)+sizeof(int)); /* allow for header and
						        trailer*/
  if (!p)			/* can't get space */
    {
      if (b)
	return((pointer)p);
      else
	{
	  fprintf(stderr,"Myalloc failure - heap full\n");
	  fprintf(stderr,"Requested %u at line %d of %s\n",n,line,file);
	  abort();
	}
    }
  m = (maptr) p;		/* set up the header */
  m->smagic = SMAGIC;
  m->len = n;
  m->len1 = n1;
  m->file = file;
  m->line = line;
  m->last = mahead.last;
  m->next = &mahead;
  m->last->next = m;
  mahead.last = m;
  p += sizeof(struct mablock);	/* now point at the block to be returned */
  for (i=0; i< n; i++)		/* 'clear' it all */
    p[i] = EMPTY;
  for (;i<n1;i++)		/* mark any extra space as such */
    p[i] = FILLER;
  *(int *)(p+n1) = EMAGIC;
  return((pointer )p);			/* return the block */
}

/* Free a block allocated by myalloc, checking for various errors */

void myfreef(p,file,line)
     pointer p;
     char *file;
     int line;
{
  maptr m;
  bool failed = false;
  int i,tl;
  if (!p)			/* The commonest one first */
    {
      fprintf(stderr,"Attempt to free NULL at line %d of %s\n",line,file);
      abort();
    }
  m = (maptr) (((char *)p)-sizeof(struct mablock)); /* look at the header */
  if (m->smagic != SMAGIC)
    {
      fprintf(stderr,"Magic number violation at start of block\n");
      failed = true;
    }
  if (m->next->smagic != SMAGIC || m->last->smagic != SMAGIC)
    {
      fprintf(stderr,"Improperly linked block\n");
      failed = true;
    }
  if (*(int *)(((char *)p)+m->len1) != EMAGIC) /* and the trailer */
    {
      fprintf(stderr, "Magic number violation at end of block\n");
      failed = true;
    }
  for (i=m->len; i < m->len1; i++) /* and the filler, if any */
    if (((char *)p)[i] != FILLER)
      {
	fprintf(stderr, "Filler space overwritten\n");
	failed = true;
      }
  if (failed)			/* if there were any errors then die */
    {
      fprintf(stderr,"Myfree failed at line %d of %s\n",line,file);
      abort();
    }
  m->last->next = m->next;	/* unlink */
  m->next->last = m->last;
  tl = m->len1+sizeof(struct mablock)+sizeof(int); /* total length */
  for (i=0; i< tl; i++)		/* clear the block */
    ((char *) m)[i] = EMPTY;
  free((pointer )m);		/* free it */
  return;
}

#include <ctype.h>

void DumpHeap(output)
     int output;
{
  maptr m;
  unsigned char *p;
  int i,topr,j;
  char hexd[49];
  char ascd[16];
  unsigned int c;
  for (m = mahead.next; m != &mahead; m = m->next)
    {
      if (m->smagic != SMAGIC)
	{
	  fprintf(stderr,"Bad block found in block from %s line %d\n",
		  m->file, m->line);
	  abort();
	}
      p = ((unsigned char *)m)+sizeof(mahead);
      if (*(int *)(p+m->len1) != EMAGIC)
	{
	  fprintf(stderr,
		  "Bad end magic number found\n in block from %s line %d\n",
		  m->file, m->line);
	  abort();
	}
      if (output)
	{
	  printf("\naddress %p request size %u allocated size %u\n",
		 (pointer)p,m->len,m->len1);
	  printf("  allocated from %s line %d\n",m->file,m->line);
	  if (output > 1)
	    {
	      for (i=0, topr = m->len; i < m->len; i += 16, topr -= 16)
		{
		  for (j=0; j < min(16,topr); j++)
		    {
		      c = p[i+j];
		      sprintf(hexd+3*j,"%02X ",c);
		      ascd[j] = isprint(c) ? c : '?';
		    }
		  for (j = 0; j < min(48,3*topr); j++)
		    putchar(hexd[j]);
		  for (; j <60; j++)
		    putchar(' ');
		  for (j = 0; j < min(16,topr); j++)
		    putchar(ascd[j]);
		  putchar('\n');
		}
	    }
	}
    }
}

