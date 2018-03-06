/* ========================== C MeatAxe =============================
   zpseed.c - Produce seed.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zpseed.c,v 1.2 1997/09/11 15:44:17 gap Exp $
 *
 * $Log: zpseed.c,v $
 * Revision 1.2  1997/09/11 15:44:17  gap
 * New version 2.2.3. AH
 *
 * Revision 2.1  1994/05/19  11:37:13  mringe
 * Benutze Smalloc(). zmkpivot() erzeugt keinen Fehler,
 * wenn die Matrix nicht in Staffelform ist.
 *
 * Revision 2.1  1994/05/19  11:37:13  mringe
 * Benutze Smalloc(). zmkpivot() erzeugt keinen Fehler,
 * wenn die Matrix nicht in Staffelform ist.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.8  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.7  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.6  1993/10/05  23:35:02  mringe
 * zerrno eliminiert.
 *
 * Revision 1.5  1993/08/24  12:39:48  mringe
 * Benutze LONG_MAX aus limits.h. Dok. entfernt.
 *
 * Revision 1.4  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.3  1993/06/30  11:32:49  mringe
 * *** empty log message ***
 *
 * Revision 1.2  1993/06/29  07:46:36  mringe
 * Dokumentation
 *
 * Revision 1.1  1993/02/10  19:40:54  mringe
 * Initial revision
 *
 * Revision 1.1  1993/02/10  19:40:54  mringe
 * Initial revision
 *
 * Revision 1.2  1993/01/08  13:35:33  mringe
 * MAXVECTOR (=100000) eingebaut.
 *
 * Revision 1.1  1992/09/23  07:47:05  mringe
 * Initial revision
 */

#include <limits.h>
#include <stdlib.h>
#include "meataxe.h"





/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

PTR zpseed_vec = NULL;		/* The seed vector */

static long dim;		/* Dimensino of basis */
static PTR basis;		/* The basis vectors */
static FEL *coeff;		/* Coefficients */
static long lastnum;		/* Number of last seed vector */
static long maxnum;		/* Maximal valid vector number + 1 */

/* ------------------------------------------------------------------
   zpseed_free() - Clean up
   ------------------------------------------------------------------ */

void zpseed_free()

{
    if (zpseed_vec != NULL) free(zpseed_vec);
    if (coeff != NULL) free(coeff);
}


/* ------------------------------------------------------------------
   zpseed_init() - Initialize the seed vector generator
   ------------------------------------------------------------------ */

int zpseed_init(sdim,sbasis)
long sdim;
PTR sbasis;

{
    int i;

    zpseed_free();
    zpseed_vec = zalloc((long)1);
    coeff = NALLOC(FEL,sdim);
    if (coeff == NULL)
    {
	MTXFAIL(ERR_NOMEM,-1);
    }
    dim = sdim;
    basis = sbasis;
    for (i = 0; i < (int) dim; ++i) coeff[i] = F_ZERO;
    zmulrow(zpseed_vec,F_ZERO);
    lastnum = 0;
    maxnum = 2;
    for (i = 1; i < dim; ++i)
    {
	/*if (maxnum * zfl > ZPSEED_MAX || maxnum * zfl < 0)
	    break;*/
	if (maxnum > LONG_MAX / zfl) break;
	maxnum *= zfl;
    }
    return 0;
}


/* ------------------------------------------------------------------
   zpseed_make() - Make a seed vector given its number
   ------------------------------------------------------------------ */

long zpseed_make(num)
long num;

{
    int i;
    PTR x = basis;

    if (num < 0 || num >= maxnum)
    {
	MTXFAIL(ERR_RANGE,-1);
    }
    lastnum = num;
    for (i = 0; i < (int) dim; ++i)
    {
	FEL co = zitof(num % zfl);
	num /= zfl;
	if (co != coeff[i])
	{
	    FEL delta = zsub(co,coeff[i]);
	    zaddmulrow(zpseed_vec,x,delta);
	    coeff[i] = co;
	}
	zadvance(&x,(long)1);
    }
    return lastnum;
}


/* ------------------------------------------------------------------
   zpseed_next() - Make the next seed vector
   ------------------------------------------------------------------ */

long zpseed_next()

{
    long w = lastnum+1, x, i;

    for (x = 1; (i = w/x) >= zfl; x *= zfl);
    if (i != 1) w = x * zfl;

    if (w >= maxnum) return -1;	/* Keine Fehlerbedingung! */
    return zpseed_make(w);
}


