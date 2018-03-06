/* ========================== C MeatAxe =============================
   sequence.c - The sequence data type.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: sequence.c,v 1.2 1997/09/11 15:43:22 gap Exp $
 *
 * $Log: sequence.c,v $
 * Revision 1.2  1997/09/11 15:43:22  gap
 * New version 2.2.3. AH
 *
 * Revision 1.6  1994/11/28  16:39:58  mringe
 * Neue Namen: SFOpen() und SFSeek()
 *
 * Revision 1.5  1994/06/29  07:51:40  mringe
 * Neu: dup und print
 *
 * Revision 1.4  1994/05/18  10:06:18  mringe
 * Benutze NALLOC
 *
 * Revision 1.3  1994/02/21  13:08:12  mringe
 * neu: seq_remove().
 *
 * Revision 1.2  1994/02/15  13:39:15  mringe
 * Benutze SFSeek().
 *
 * Revision 1.1  1994/02/12  04:10:13  mringe
 * Initial revision
 *
 */

#include <string.h>
#include <stdlib.h>
#include "meataxe.h"



/* ------------------------------------------------------------------
   seq_desc() - Describe
   ------------------------------------------------------------------ */
    
char *seq_desc(s)
sequence_t *s;

{
    static char buf[50];
    sprintf(buf,"sequence of length %ld",(long)s->len);
    return buf;
}



/* ------------------------------------------------------------------
   seq_print() - Print a sequence
   ------------------------------------------------------------------ */

void seq_print(name,s)
char *name;
sequence_t *s;

{
    int i;

    if (name != NULL) printf("%s=",name);
    printf("[");
    for (i = 0; i < s->len; ++i)
	mtxprint(s->buf[i]);
	if (i < s->len-1) printf(",");
    printf("]");
    if (name != NULL) printf("\n");
    fflush(stdout);
}



/* ------------------------------------------------------------------
   seq_alloc() - Allocate a new sequence
   seq_free() - Free a sequence
   ------------------------------------------------------------------ */

sequence_t *seq_alloc(size)
size_t size;

{
    sequence_t *n;

    n = (sequence_t *) malloc(sizeof(sequence_t));
    if (n != NULL)
    {
        n->id = T_SEQUENCE;
	if (size < 1) size = 1; /* avoid malloc(0) */
    	if ((n->buf = NALLOC(void *,size)) != NULL)
	{
	    n->len = 0;
	    n->max = size;
	    return n;
	}
	free(n);
    }
    MTXFAIL(ERR_NOMEM,(sequence_t *)NULL);
}


void seq_free(s)
sequence_t *s;

{
    int i;
    for (i = 0; i < s->len; ++i)
	mtxfree(s->buf[i]);
    free(s);
}

/* ------------------------------------------------------------------
   seq_dup() - Duplicate a sequence
   ------------------------------------------------------------------ */

sequence_t *seq_dup(s)
sequence_t *s;

{
    sequence_t *x;
    int i;

    if ((x = seq_alloc(s->len)) == NULL) return NULL;
    for (i = 0; i < s->len; ++i)
    {
	void *y = mtxdup(s->buf[i]);
	if (y == NULL)
	{
	    seq_free(x);
	    return NULL;
	}
	seq_insert(x,i,y);
    }
    return x;
}




/* ------------------------------------------------------------------
   seq_insert() - Insert an object at position <pos>, or append if
       <pos> is negative.
   ------------------------------------------------------------------ */

int seq_insert(s,pos,x)
sequence_t *s;
int pos;
void *x;

{
    if (s->len == s->max)
    {
	size_t newsize = s->max + 8 - s->max % 8;
	void *newbuf;
	if ((newbuf = realloc(s->buf,newsize)) == NULL)
	    MTXFAIL(ERR_NOMEM,-1);
	s->buf = newbuf;
    }
    if (pos < 0)
	s->buf[s->len] = x;
    else
    {
	int k;
	void **b = s->buf;
	for (k = s->len; k > pos; --k)
	    b[k] = b[k-1];
	b[pos] = x;
    }
    ++s->len;
    return 0;
}



/* ------------------------------------------------------------------
   seq_read() - Read a sequence
   ------------------------------------------------------------------ */

sequence_t *seq_read(f)
FILE *f;

{
    long l;
    int i;
    sequence_t *s;

    if (zreadlong(f,&l,1) != 1)
	MTXFAIL(ERR_FILEREAD,NULL);
    if ((s = seq_alloc((size_t)l)) == NULL)
	MTXFAIL(ERR_NOMEM,NULL);
    for (i = (int) l; i > 0; --i)
    {
	void *x = mtxread(f);
	if (x == NULL) return NULL;
	seq_insert(s,-1,x);
    }
    return s;
}



/* ------------------------------------------------------------------
   seq_write() - Write a sequence
   ------------------------------------------------------------------ */

int seq_write(f, s)
FILE *f;
sequence_t *s;

{
    long l = s->len;
    int i;

    if (zwritelong(f,&l,1) != 1)
	MTXFAIL(ERR_FILEWRITE,-1);
    for (i = 0; i < s->len; ++i)
    {
	if (mtxwrite(f,s->buf[i]) != 0)
	    return -1;
    }
    return 0;
}



/* ------------------------------------------------------------------
   seq_remove() - Remove an object from a sequence.
   ------------------------------------------------------------------ */

int seq_remove(s,pos)
sequence_t *s;
int pos;

{
    int i;

    if (pos < 0 || pos >= s->len)
	MTXFAIL(ERR_RANGE,-1);
    mtxfree(s->buf[pos]);
    for (i = pos;  i < s->len -1; ++i)
	s->buf[i] = s->buf[i+1];
    s->buf[i] = NULL;
    --s->len;
    return 0;
}

