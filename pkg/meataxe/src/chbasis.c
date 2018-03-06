/* ========================== C MeatAxe =============================
   Change basis.

   (C) Copyright 1994 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: chbasis.c,v 1.2 1997/09/11 15:42:41 gap Exp $
 *
 * $Log: chbasis.c,v $
 * Revision 1.2  1997/09/11 15:42:41  gap
 * New version 2.2.3. AH
 *
 * Revision 1.1  1994/05/18  10:12:01  mringe
 * Initial revision
 *
 * Revision 1.1  1994/05/18  10:12:01  mringe
 * Initial revision
 *
 */


#include <stdlib.h>
#include "meataxe.h"


/* ------------------------------------------------------------------
   chbasis() - Change basis: newgen[i] = basis * gen[i] * basis^(-1)
   newgen = gen is ok!
   ------------------------------------------------------------------ */

int chbasis(basis,ngen,gen,newgen)
matrix_t *basis;
int ngen;
matrix_t *gen[];
matrix_t *newgen[];

{
    matrix_t *bi, *tmp;
    int i;

    if ((bi = matinv(basis)) == NULL) return -1;
    for (i = 0; i < ngen; ++i)
    {
	tmp = matdup(basis);
	matmul(tmp,gen[i]);
	matmul(tmp,bi);
	if (newgen == gen) matfree(gen[i]);
	newgen[i] = tmp;
    }
    matfree(bi);
    return 0;
}


