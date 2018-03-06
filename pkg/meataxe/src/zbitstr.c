/* ========================== C MeatAxe =============================
   zbitstr.c - Bit string operations

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zbitstr.c,v 1.2 1997/09/11 15:43:36 gap Exp $
 *
 * $Log: zbitstr.c,v $
 * Revision 1.2  1997/09/11 15:43:36  gap
 * New version 2.2.3. AH
 *
 * Revision 2.4  1994/06/29  07:50:19  mringe
 * Neu: bs_dup()
 *
 * Revision 2.4  1994/06/29  07:50:19  mringe
 * Neu: bs_dup()
 *
 * Revision 2.3  1994/05/20  09:06:16  mringe
 * Benutze Smalloc().
 *
 * Revision 2.2  1994/02/12  04:10:13  mringe
 * UMFANGREICHE AENDERUNGEN AN VIELEN DATENTYPEN.
 *
 * Revision 2.1  1993/12/02  17:59:28  mringe
 * Ersetze bitstring_t durch bitstring_t *.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.15  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.14  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.13  1993/10/05  23:35:02  mringe
 * zerrno eliminiert.
 *
 * Revision 1.12  1993/09/30  12:58:01  mringe
 * bs_size als globale Variable.
 *
 * Revision 1.11  1993/08/12  16:14:23  mringe
 * Include <string.h>.
 *
 * Revision 1.10  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.9  1993/07/06  13:59:43  mringe
 * Neu: bs_minus()
 *
 * Revision 1.8  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.7  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.6  1993/01/14  18:58:58  mringe
 * *** empty log message ***
 *
 * Revision 1.5  1992/07/22  07:10:30  mringe
 * Changed 'global.h' to 'lattice.h'
 *
 * Revision 1.4  1992/07/13  11:57:02  mringe
 * Neu: bs_and()
 *
 * Revision 1.3  1992/07/04  12:01:49  mringe
 * Benutze long int Operationen in bs_or() und bs_issub()
 *
 * Revision 1.2  1992/07/04  10:18:42  mringe
 * bs_copy() und bs_issub(): Benutze long int Operationen.
 *
 * Revision 1.1  1992/05/24  08:03:16  mringe
 * Initial revision
 *
 */

#include <string.h>
#include <stdlib.h>
#include "meataxe.h"




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static unsigned char mask[8] = {1,2,4,8,16,32,64,128};
static unsigned char clmask[8] =
	{(unsigned char)~1,(unsigned char)~2,(unsigned char)~4,
	(unsigned char)~8,(unsigned char)~16,(unsigned char)~32,
        (unsigned char)~64,(unsigned char)~128};

static int matchcount[256] =
	{	0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
		1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
		1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
		2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
		1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
		2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
		2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
		3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
		1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
		2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
		2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
		3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
		2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
		3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
		3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
		4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8
	};
size_t bs_size = -1;
size_t bs_len = -1;

/* ------------------------------------------------------------------
   bs_print() - Print a bit string
   ------------------------------------------------------------------ */
    
void bs_print(b)
bitstring_t *b;

{
    int i;
    for (i = 0; i < bs_len; ++i)
	putc(bs_test(b,i) ? '+' : '.',stdout);
}

/* ------------------------------------------------------------------
   bs_desc() - Describe a bit string
   ------------------------------------------------------------------ */
    
char *bs_desc(b)
bitstring_t *b;

{
    static char buf[20];
    strcpy(buf,"bit string");
    return buf;
}


/* ------------------------------------------------------------------
   bs_setlen() - Set bit string length
   bs_alloc() - Allocate a new bit string
   bs_free() - Free a bit string
   bs_reset() - Clear all bits in a bit string
   ------------------------------------------------------------------ */

void bs_setlen(l)
int l;

{
    bs_len = l;
    bs_size = ((l-1) / 8) + 1;
}


bitstring_t *bs_alloc()

{
    bitstring_t *n;

    n = (bitstring_t *) Smalloc(sizeof(bitstring_t) + bs_size - 1);
    if (n == NULL) MTXFAIL(ERR_NOMEM,(bitstring_t *)NULL);
    n->id = T_BITS;
    memset(n->buf,0,(size_t)bs_size);
    return n;
}

bitstring_t *bs_dup(s)
bitstring_t *s;

{
    bitstring_t *n;

    n = (bitstring_t *) Smalloc(sizeof(bitstring_t) + bs_size - 1);
    if (n == NULL) MTXFAIL(ERR_NOMEM,(bitstring_t *)NULL);
    n->id = T_BITS;
    memcpy(n->buf,s->buf,(size_t)bs_size);
    return n;
}


void bs_free(b)
bitstring_t *b;

{
    free(b);
}


void bs_reset(x)
bitstring_t *x;

{
    memset(x->buf,0,bs_size);
}


/* ------------------------------------------------------------------
   bs_cmp() - Compare two bit strings
   bs_cpy() - Copy a bit string
   ------------------------------------------------------------------ */

int bs_cmp(a,b)		/* Compare two bit strings */
bitstring_t *a, *b;

{
    return memcmp(a->buf,b->buf,bs_size);
}


void bs_cpy(dest, src)
bitstring_t *dest, *src;

{	memcpy(dest->buf,src->buf,bs_size);
}


/* ------------------------------------------------------------------
   bs_set() - Set a bit
   bs_clear() - Clear a bit
   bs_test() - Check if a bit is set
   ------------------------------------------------------------------ */

void bs_set(b,i)
bitstring_t *b;
int i;

{	b->buf[i / 8] |= mask[i % 8];
}


void bs_clear(b, i)
bitstring_t *b;
int i;

{	b->buf[i / 8] &= clmask[i % 8];
}



int bs_test(b, i)
bitstring_t *b;
int i;

{	return ((b->buf[i / 8] & mask[i % 8]) != 0);
}


/* ------------------------------------------------------------------
   bs_read() - Read a bit string from a file
   bs_write() - Write a bit string to a file
   ------------------------------------------------------------------ */

bitstring_t *bs_read(f)
FILE *f;

{
    bitstring_t *b;

    b = bs_alloc();
    if (fread(b->buf,sizeof(unsigned char),bs_size,f) != bs_size)
    {
	bs_free(b);
	MTXFAIL(ERR_FILEREAD,(bitstring_t *)NULL);
    }
    return b;
}

int bs_write(f,b)
FILE *f;
bitstring_t *b;

{
    if (fwrite(b->buf,sizeof(unsigned char),bs_size,f) != bs_size)
	MTXFAIL(ERR_FILEWRITE,-1);
    return 0;
}


/* ------------------------------------------------------------------
   bs_and() - Logical AND: dest &= src
   ------------------------------------------------------------------ */

void bs_and(dest, src)
bitstring_t *dest, *src;

{
    register int i;
    /* gefaehrlich: buf liegt evtl nicht auf einer longword grenze*/
    register unsigned long *dl = (unsigned long *)dest->buf;
    register unsigned long *sl = (unsigned long *)src->buf;
    register unsigned char *d;
    register unsigned char *s;

    for (i = bs_size; i > sizeof(long); i -= sizeof(long))
	*dl++ &= *sl++;
    d = (unsigned char *) dl;
    s = (unsigned char *) sl;
    for (; i > 0; --i)
	*d++ &= *s++;
}


/* ------------------------------------------------------------------
   bs_or() - Logical OR: dest |= src
   ------------------------------------------------------------------ */

void bs_or(dest, src)
bitstring_t *dest, *src;

{	register int i;
	register unsigned long *dl = (unsigned long *)dest->buf;
	register unsigned long *sl = (unsigned long *)src->buf;
	register unsigned char *d;
	register unsigned char *s;

	for (i = bs_size; i > sizeof(long); i -= sizeof(long))
		*dl++ |= *sl++;
	d = (unsigned char *) dl;
	s = (unsigned char *) sl;
	for (; i > 0; --i)
		*d++ |= *s++;
}

/* ------------------------------------------------------------------
   bs_minus() - dest &= ~src
   ------------------------------------------------------------------ */

void bs_minus(dest, src)
bitstring_t *dest, *src;

{	register int i;
	register unsigned long *dl = (unsigned long *)dest->buf;
	register unsigned long *sl = (unsigned long *)src->buf;
	register unsigned char *d;
	register unsigned char *s;

	for (i = bs_size; i > sizeof(long); i -= sizeof(long))
		*dl++ &= ~*sl++;
	d = (unsigned char *) dl;
	s = (unsigned char *) sl;
	for (; i > 0; --i)
		*d++ &= ~*s++;
}



/* ------------------------------------------------------------------
   bs_match() - Returns number of bits set in both bit strings
   ------------------------------------------------------------------ */

int bs_match(x, y)
bitstring_t *x, *y;

{
    register int i, k = 0;
    register unsigned char *d = (unsigned char *)x->buf;
    register unsigned char *s = (unsigned char *)y->buf;
    register unsigned char t;

    for (i = bs_size; i > 0 && k < 2; i--)
    {
	if ((t = (*s++ & *d++)) != 0)
	    k += matchcount[t];
    }
    return k;
}


/* ------------------------------------------------------------------
   bs_issub() - Returns 1 if 'a' is a subset of 'b'
   ------------------------------------------------------------------ */

int bs_issub(a,b)
bitstring_t *a, *b;

{	register int i;
	register unsigned long *al = (unsigned long *)a->buf;
	register unsigned long *bl = (unsigned long *)b->buf;
	register unsigned char *ac;
	register unsigned char *bc;

	for (i = bs_size; i > sizeof(long); i -= sizeof(long))
	{	if ((*al ^ (*al & *bl)) != 0)
			return 0;
		al++;
		bl++;
	}
	ac = (unsigned char *) al;
	bc = (unsigned char *) bl;
	for (; i > 0; --i)
	{	if ((*ac ^ (*ac & *bc)) != 0)
			return 0;
		ac++;
		bc++;
	}
	return 1;
}


