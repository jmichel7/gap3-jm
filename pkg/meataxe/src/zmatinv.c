/* ========================== C MeatAxe =============================
   zmatinv.c -  Matrix inversion

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */

/* $Id: zmatinv.c,v 1.2 1997/09/11 15:44:00 gap Exp $
 *
 * $Log: zmatinv.c,v $
 * Revision 1.2  1997/09/11 15:44:00  gap
 * New version 2.2.3. AH
 *
 * Revision 1.4  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 1.3  1993/12/08  11:48:50  mringe
 * Compiler warnings.
 *
 * Revision 1.2  1993/10/21  08:37:53  mringe
 * Compiler warnings.
 *
 * Revision 1.1  1993/10/20  18:17:07  mringe
 * Initial revision
 *
 */


#include "meataxe.h"



/* ------------------------------------------------------------------
   zmatinv() - Matrix inversion. On return, the original matrix is
	destroyed and <result> contains the inverse.

   Return value: 0 = OK, 1 = Error (matrix is singular)
   ------------------------------------------------------------------ */

int zmatinv(mat,result)
PTR mat, result;

{
    PTR xj1, xj2, xk1, xk2;
    FEL f1 = F_ZERO, f2;
    long j, k;

    /* Initialize result = identity matrix
       ----------------------------------- */
    for (j = 1, xj1 = result; j <= znoc; ++j, zadvance(&xj1,1))
    {
	zmulrow(xj1,F_ZERO);
	zinsert(xj1,j,F_ONE);
    }

    /* Matrix inversion
       ---------------- */
    xj1 = mat;
    xj2 = result;
    for (j = 1; j <= znoc; ++j)
    {

        for (xk1 = xj1, k = j;
	     k <= znoc && (f1 = zextract(xk1,j)) == F_ZERO;
	     ++k, zadvance(&xk1,(long)1));
	if (f1 == F_ZERO) MTXFAIL(ERR_DIV0,-1);
	if (k > j)	/* Swap rows */
	{
	    zswaprow(xk1,xj1);
	    xk2 = xj2;
	    zadvance(&xk2,k-j);
	    zswaprow(xk2,xj2);
	}
	f2 = zinv(f1);
	zmulrow(xj1,f2);
	zmulrow(xj2,f2);
	xk1 = mat;
	xk2 = result;
	for (k = 1; k <= znoc; ++k)
	{
	    if (k != j)
	    {
		f1 = zneg(zextract(xk1,j));
		zaddmulrow(xk1,xj1,f1);
		zaddmulrow(xk2,xj2,f1);
	    }
	    zadvance(&xk1,(long)1);
	    zadvance(&xk2,(long)1);
	}
	zadvance(&xj1,(long)1);
	zadvance(&xj2,(long)1);
    }
    return 0;
}

