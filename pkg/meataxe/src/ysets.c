/* ========================== C MeatAxe =============================
   ysets.c - Integer sets (= sorted lists)

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: ysets.c,v 1.2 1997/09/11 15:43:33 gap Exp $
 *
 * $Log: ysets.c,v $
 * Revision 1.2  1997/09/11 15:43:33  gap
 * New version 2.2.3. AH
 *
 * Revision 2.4  1994/06/29  07:51:40  mringe
 * Neu: dup und print
 *
 * Revision 2.4  1994/06/29  07:51:40  mringe
 * Neu: dup und print
 *
 * Revision 2.3  1994/05/20  09:09:32  mringe
 * Benutze Smalloc() und Srealloc().
 *
 * Revision 2.2  1994/03/13  13:27:01  mringe
 * Maschinenunabhaengiges Format fuer Permutationen.
 *
 * Revision 2.1  1994/02/12  04:10:13  mringe
 * UMFANGREICHE AENDERUNGEN AN VIELEN DATENTYPEN.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.7  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.6  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.5  1993/10/05  19:02:08  mringe
 * yerror eliminiert.
 *
 * Revision 1.4  1993/10/05  18:58:31  mringe
 * Neue fehlermeldungen.
 *
 * Revision 1.3  1993/08/11  08:24:01  mringe
 * Changed zerrno to yerrno
 *
 * Revision 1.2  1993/08/11  08:19:47  mringe
 * Sets von ZZZ in YYY Library verlagert.
 *
 * Revision 1.1  1993/08/11  08:10:59  mringe
 * Initial revision
 *
 */


#include <stdlib.h>
#include "meataxe.h"





/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static size_t init_max = 10;
static size_t blksize = 10;



/* ------------------------------------------------------------------
   set_desc() - Describe
   ------------------------------------------------------------------ */
    
char *set_desc(s)
set_t *s;

{
    static char buf[50];
    sprintf(buf,"set of %u integers",(unsigned) s->len);
    return buf;
}


/* ------------------------------------------------------------------
   set_print() - Print
   ------------------------------------------------------------------ */
    
void set_print(name,s)
char *name;
set_t *s;

{
    int i;
    if (name != NULL) printf("%s=",name);
    printf("{");
    for (i = 0; i < s->len; ++i)
    {
	if (i > 0) printf(",");
	printf("%ld",s->buf[i]);
    }
    printf("}");
    if (name != NULL) printf("\n");
    fflush(stdout);
}


/* ------------------------------------------------------------------
   set_allocstrategy() - Set allocation strategy
   ------------------------------------------------------------------ */

int set_allocstrategy(first,blocksize)
size_t first;
size_t blocksize;

{
    if (first < 1 || blocksize < 1)
    {
    	MTXFAIL(ERR_RANGE,-1);
    }
    init_max = first;
    blksize = blocksize;
    return 0;
}



/* ------------------------------------------------------------------
   set_alloc() - Allocate an integer set
   set_free() - Free an integer set
   ------------------------------------------------------------------ */

set_t *set_alloc()

{
    set_t *x;

    x = ALLOC(set_t);
    if (x == NULL)
 	MTXFAIL(ERR_NOMEM,NULL);
    x->id = T_SET;
    x->len = 0;
    x->max = init_max;
    x->buf = NALLOC(long,init_max);
    if (x->buf == NULL)
    {
	free(x);
 	MTXFAIL(ERR_NOMEM,NULL);
    }

    return x;
}


void set_free(x)
set_t *x;

{
    free(x->buf);
    free(x);
}


/* ------------------------------------------------------------------
   set_dup() - Duplicate
   ------------------------------------------------------------------ */

set_t *set_dup(set_t *s)

{
    set_t *x;
    int i;

    x = ALLOC(set_t);
    if (x == NULL)
 	MTXFAIL(ERR_NOMEM,NULL);
    x->id = T_SET;
    x->max = x->len = s->len;
    x->buf = NALLOC(long,x->max);
    if (x->buf == NULL)
    {
	free(x);
 	MTXFAIL(ERR_NOMEM,NULL);
    }
    for (i = 0; i < s->len; ++i)
	x->buf[i] = s->buf[i];
    return x;
}



/* ------------------------------------------------------------------
   set_insert() - Insert an element
   ------------------------------------------------------------------ */

int set_insert(set,elem)
set_t *set;
long elem;

{
    int i,k;
    long *l,*m;

    for (l = set->buf, i = set->len; i > 0 && *l < elem; --i, ++l);
    if (i > 0 && *l == elem) return 0;	/* Ok, already in set */
    if (set->len >= set->max - 1)	/* Need extend the set? */
    {
	size_t newmax = set->max + blksize;
	long *newbuf = (long *) Srealloc(set->buf,newmax*sizeof(long));
	if (newbuf == NULL)
	    MTXFAIL(ERR_NOMEM,-1);
	set->max = newmax;
	set->buf = newbuf;
    }
    for (m = set->buf+set->len-1, k = i; k > 0; --k,--m)
	m[1] = m[0];
    *l = elem;
    ++set->len;
    return 0;
}


/* ------------------------------------------------------------------
   set_contains() - Check if a set contains a given element
   ------------------------------------------------------------------ */

int set_contains(set,elem)
set_t *set;
long elem;

{
    int i;
    long *l = set->buf;
 
    for (i = set->len; i > 0 && *l < elem; --i, ++l);
    return (i > 0 && *l == elem);
}


/* ------------------------------------------------------------------
   set_read() - Read an integer set from a file
   set_write() - Write an integer set to a file
   ------------------------------------------------------------------ */

set_t *set_read(f)
FILE *f;

{
    set_t *x;
    long l;

    if (zreadlong(f,&l,1) != 1)
	MTXFAIL(-1,NULL);
    x = ALLOC(set_t);
    if (x == NULL)
	MTXFAIL(ERR_NOMEM,NULL);
    x->id = T_SET;
    x->max = x->len = (size_t)l;
    x->buf = NALLOC(long,x->len);
    if (x->buf == NULL)
    {
	free(x);
	MTXFAIL(ERR_NOMEM,NULL);
    }
    if (zreadlong(f,x->buf,x->len) != x->len)
    {
	set_free(x);
	MTXFAIL(ERR_FILEREAD,NULL);
    }
    return x;
}


int set_write(f,set)
FILE *f;
set_t *set;

{
    long l = set->len;

    if (zwritelong(f,&l,1) != 1)
	MTXFAIL(ERR_FILEWRITE,-1);
    if (zwritelong(f,set->buf,set->len) != set->len)
	MTXFAIL(ERR_FILEWRITE,-1);
    return 0;
}


