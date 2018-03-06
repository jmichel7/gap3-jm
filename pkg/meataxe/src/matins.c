/* ========================== C MeatAxe =============================
   polins.c - Insert a matrix into a polynomial

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: matins.c,v 1.2 1997/09/11 15:42:59 gap Exp $
 *
 * $Log: matins.c,v $
 * Revision 1.2  1997/09/11 15:42:59  gap
 * New version 2.2.3. AH
 *
 * Revision 1.4  1994/09/10  10:47:44  mringe
 * Neu: matins_() - Ueberschreibt das Argument mit dem Ergebnis.
 *
 * Revision 1.3  1994/07/23  16:46:18  mringe
 * Profiling.
 *
 * Revision 1.2  1994/06/16  14:16:17  mringe
 * Profiling.
 *
 * Revision 1.1  1994/04/12  13:44:13  mringe
 * Initial revision
 *
 */


#include <string.h>
#include <stdlib.h>
#include "meataxe.h"

long tmatinsert = 0;

/* ------------------------------------------------------------------
   matinsert_() - Insert a matrix into a polynomial (destructive)
   Returns its argument or NULL on error.
   ------------------------------------------------------------------ */
    

matrix_t *matinsert_(mat,pol)
matrix_t *mat;
poly_t *pol;

{
    matrix_t *x = NULL;
    int i;
    long nor = mat->nor;
    long l;
    PTR v;
    FEL f;
    PROFILE_BEGIN(t);

    if (nor != mat->noc) MTXFAIL(ERR_NOTSQUARE,NULL);
    if (mat->fl != pol->fl) MTXFAIL(ERR_INCOMPAT,NULL);

    zsetfield(mat->fl);
    zsetlen(mat->noc);

    /* Special case: p(x) = 0
       ---------------------- */
    if (pol->deg == -1)
    {
	free(mat->d);
	mat->d = zalloc(mat->nor);
    	PROFILE_END(t,MatInsert);
	return mat;
    }

    /* Evaluate p(A)
       ------------- */
    if (pol->deg > 1) x = matdup(mat);
    if ((f = pol->buf[pol->deg]) != F_ONE)
    {
	for (l = nor, v = mat->d; l > 0; --l, zadvance(&v,(long)1))
	    zmulrow(v,f);
    }
    for (i = pol->deg - 1; i >= 0; --i)
    {
        if ((f = pol->buf[i]) != F_ZERO)
        {
	    for (l = 1, v = mat->d; l <= nor; ++l, zadvance(&v,(long)1))
		zinsert(v,l,zadd(zextract(v,l),f));
        }
	if (i > 0) matmul(mat,x);
    }
    if (pol->deg > 1) matfree(x);
    PROFILE_END(t,MatInsert);
    return mat;
}


/* ------------------------------------------------------------------
   matinsert() - Insert a matrix into a polynomial. Non-destructive
   version.
   ------------------------------------------------------------------ */
    
matrix_t *matinsert(mat,pol)
matrix_t *mat;
poly_t *pol;

{
    matrix_t *x;
    int i;
    long nor = mat->nor;
    long l;
    PTR v;
    FEL f;
    PROFILE_BEGIN(t);

    if (nor != mat->noc) MTXFAIL(ERR_NOTSQUARE,NULL);
    if (mat->fl != pol->fl) MTXFAIL(ERR_INCOMPAT,NULL);

    /* Special case: p(x) = 0
       ---------------------- */
    if (pol->deg == -1)
    {
	x = matalloc(mat->fl,mat->nor,mat->nor);
    	PROFILE_END(t,MatInsert);
	return x;
    }

    /* Evaluate p(A)
       ------------- */
    x = matdup(mat);
    if ((f = pol->buf[pol->deg]) != F_ONE)
    {
	for (l = nor, v = x->d; l > 0; --l, zadvance(&v,(long)1))
	    zmulrow(v,f);
    }
    for (i = pol->deg - 1; i >= 0; --i)
    {
        if ((f = pol->buf[i]) != F_ZERO)
        {
	    for (l = 1, v = x->d; l <= nor; ++l, zadvance(&v,(long)1))
		zinsert(v,l,zadd(zextract(v,l),f));
        }
	if (i > 0) matmul(x,mat);
    }
    PROFILE_END(t,MatInsert);
    return x;
}

