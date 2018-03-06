/* ========================== C MeatAxe =============================
   zgauss.c - Gauss elimination and related functions.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zgauss.c,v 1.2 1997/09/11 15:43:53 gap Exp $
 *
 * $Log: zgauss.c,v $
 * Revision 1.2  1997/09/11 15:43:53  gap
 * New version 2.2.3. AH
 *
 * Revision 2.5  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.4  1994/05/20  12:41:25  mringe
 * Prototypen.
 *
 * Revision 2.3  1994/05/19  16:02:24  mringe
 * Schnellere Version von zmkpivot().
 *
 * Revision 2.2  1994/04/09  13:51:23  mringe
 * memset()-Bug ???
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.10  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.9  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.8  1993/10/05  23:35:02  mringe
 * zerrno eliminiert.
 *
 * Revision 1.7  1993/08/09  12:36:38  mringe
 * Dokumentation entfernt.
 *
 * Revision 1.6  1993/08/06  15:09:31  mringe
 * ERR_NOECH in ERR_NOTECH umbenannt.
 *
 * Revision 1.5  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.4  1993/06/14  09:29:11  mringe
 * Returnwerte korrigiert. Dikumentation.
 *
 * Revision 1.3  1993/05/12  08:38:39  mringe
 * Prototypes, die zweite.
 *
 * Revision 1.3  1993/05/12  08:38:39  mringe
 * Prototypes, die zweite.
 *
 * Revision 1.2  1993/05/12  08:35:09  mringe
 * Prototypes korrigiert.
 *
 * Revision 1.1  1993/02/10  19:40:54  mringe
 * Initial revision
 *
 * Revision 1.2  1993/01/28  07:35:51  mringe
 * *** empty log message ***
 *
 * Revision 1.1  1993/01/27  15:38:38  mringe
 * Initial revision
 *
 */

#include <stdlib.h>
#include "meataxe.h"



/* ------------------------------------------------------------------
   zmkpivot() - Find the pivot columns of a matrix. <matrix> is a
	<nor> by <noc> matrix in echelon form. The pivot columns
	are stored in <piv>, beginning with piv[1]. Thus, <piv> must
	be large enough to hold at least <nor>+1 entries.

   Return value:
      0  = ok
      -1 = Error (out of memory or matrix not in echelon form)
   ------------------------------------------------------------------ */


int zmkpivot(matrix, nor, piv)
PTR matrix;
long nor;
long *piv;

{
    PTR x;
    long i;
    FEL f;
    static unsigned char *ispiv = NULL;
    static long maxnoc;

    if (ispiv == NULL || znoc > maxnoc)
    {
	if (ispiv != NULL) free(ispiv);
	ispiv = NALLOC(unsigned char,znoc+1);
	maxnoc = znoc;
    }
    for (i = 1; i <= znoc; ++i) ispiv[(int)i] = 0;
    for (i=1, x=matrix; i<=nor; ++i, zadvance(&x,(long)1))
    {
	long pivi = zfindpiv(x,&f);
	if (pivi == 0 || ispiv[pivi])
	    MTXFAIL(ERR_NOTECH,-1);
	piv[i] = pivi;
	ispiv[pivi] = 1;
    }
    return 0;
}



/* ------------------------------------------------------------------
   zcleanrow() - Clean a row
   zcleanrow2() - Clean a row and store the operations.
   ------------------------------------------------------------------ */

#define CLEAN_PROC(extra_cmds)\
    for (i=1, x=matrix; i<=nor; ++i, zadvance(&x,(long)1))\
    {\
	FEL f = zextract(row,piv[i]);\
	if (f == F_ZERO) continue;\
	f = zdiv(f,zextract(x,piv[i]));\
	zaddmulrow(row,x,zneg(f));\
	extra_cmds;\
    }


int zcleanrow(row, matrix, nor, piv)
PTR row;
PTR matrix;
long nor;
long *piv;

{
    long i;
    PTR x;

    /*CLEAN_PROC({});*/
    for (i=1, x=matrix; i<=nor; ++i, zadvance(&x,(long)1))
    {
	FEL f = zextract(row,piv[i]);
	if (f == F_ZERO) continue;
	f = zdiv(f,zextract(x,piv[i]));
	zaddmulrow(row,x,zneg(f));
    }

    return 0;
}


int zcleanrow2(row, matrix, nor, piv, row2)
PTR row;
PTR matrix;
long nor;
long *piv;
PTR row2;

{
    long i;
    PTR x;

    CLEAN_PROC(zinsert(row2,i,f));
    return 0;
}

/* ------------------------------------------------------------------
   zmkechelon() - Convert a matrix to echelon form.

   Return value: Rank (also stored in piv[0])
   ------------------------------------------------------------------ */

long zmkechelon(matrix, nor, piv)
PTR matrix;
long nor;
long *piv;

{
    PTR x, newrow;
    long i, rank;

    rank = 0;
    newrow = matrix;
    for (i = 1, x = matrix; i <= nor; ++i, zadvance(&x,(long)1))
    {
	long newpiv;
	FEL f;

	if (rank < i-1) zmoverow(newrow,x);
	zcleanrow(newrow,matrix,rank,piv);
	newpiv = zfindpiv(newrow,&f);
	if (newpiv != 0)
	{
	    ++rank;
	    piv[rank] = newpiv;
	    zadvance(&newrow,(long)1);
	}
    }
    piv[0] = rank;
    return rank;
}


