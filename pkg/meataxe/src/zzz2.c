/* ========================== C MeatAxe =============================
   Basic row operations.

   (C) Copyright 1994 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zzz2.c,v 1.2 1997/09/11 15:44:43 gap Exp $
 *
 * $Log: zzz2.c,v $
 * Revision 1.2  1997/09/11 15:44:43  gap
 * New version 2.2.3. AH
 *
 * Revision 1.4  1994/11/25  14:08:17  mringe
 * ANSI-C, neue Namen fuer time-Funktionen.
 *
 * Revision 1.3  1994/07/29  07:41:52  mringe
 * Benutze zmulrow() zum Initialisieren in zalloc().
 *
 * Revision 1.2  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 1.1  1994/07/25  18:33:53  mringe
 * Initial revision
 *
 * Revision 1.1  1994/07/20  22:19:39  mringe
 * Initial revision
 *
 *
 */


#include <string.h>
#include <stdlib.h>
#include "meataxe.h"

#define LPR (zrowsize/sizeof(long))

/* ------------------------------------------------------------------
   zalloc()  -  Allocate memory and initialize
   ------------------------------------------------------------------ */

PTR zalloc(long nrows)

{
    PTR p;
    char *q;
    register long i;
    size_t req = (size_t)(zrowsize*nrows);

    if ((p = (PTR) Smalloc(req)) == NULL)
    {
	fprintf(stderr,"zalloc(): Request=%ld Bytes\n",(long)req);
	FATAL("OUT OF MEMORY");
    }

    /* Initialize */
    q = (char *) p;
    for (i = nrows; i > 0; --i)
    {
	zmulrow((PTR) q,F_ZERO);
	q += zrowsize;
    }
    return p;
}




/* ------------------------------------------------------------------
   zsize()  -  Size of a matrix or permutation(s)
   ------------------------------------------------------------------ */

size_t zsize(long nrows)

{	
    return ((size_t)(zrowsize * nrows));
}



/* ------------------------------------------------------------------
   zmoverow()  -  Move row/vector

   return value: none
   ------------------------------------------------------------------ */

void zmoverow(PTR dest, PTR src)

{
    register long *s = (long *) src, *d = (long *) dest, i;
    for (i = LPR; i > 0; --i)
   	 *d++ = *s++;
}


/* ------------------------------------------------------------------
   zswaprow()  -  Swap two rows

   return value: none
   ------------------------------------------------------------------ */

void zswaprow(PTR dest, PTR src)

{	register long *p1 = (long *) src, *p2 = (long *) dest, i, l;

	for (i = LPR; i > 0; --i)
	{	l = *p1;
		*p1++ = *p2;
		*p2++ = l;
	}
}



/* ------------------------------------------------------------------
   zcmprow()  -  Compare two rows
   ------------------------------------------------------------------ */

int zcmprow(PTR dest, PTR src)

{	
    register long k;
    register long *d = (long *) dest, *s = (long *) src;

    for (k = LPR; k > 0; --k)
	if (*d++ != *s++) return 1;
    return 0;
}




/* ------------------------------------------------------------------
   zadvance()  -  advance pointer

   return value: none
   ------------------------------------------------------------------ */

void zadvance(PTR *ptr, long nrows)

{
	*(long **)ptr += nrows * LPR;
}



/* ------------------------------------------------------------------
   zpermrow()  -  multiply a row by a permutation
   ------------------------------------------------------------------ */

PTR zpermrow(PTR row, long *perm, PTR result)

{
    register FEL f;
    register long i;
    register long *p = (long *)perm;

    for (i = 1; i <= znoc; ++i)
    {
	f = zextract(row,i);
	zinsert(result,*p++,f);
    }
    return result;
}



/* ------------------------------------------------------------------
   mtxinit()  -  Initialize
   ------------------------------------------------------------------ */

int mtxinit()

{
    zchar = 0;
    zgen = 0;
    zfl = 0;		/* This forces zsetlen() to read the tables! */
    SInit();
    return (ZZZVERSION);
}

