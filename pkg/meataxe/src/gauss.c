/* ========================== C MeatAxe =============================
   gauss.c -  Null-space, echelon form etc.

   (C) Copyright 199$ Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */

/* $Id: gauss.c,v 1.2 1997/09/11 15:42:51 gap Exp $
 *
 * $Log: gauss.c,v $
 * Revision 1.2  1997/09/11 15:42:51  gap
 * New version 2.2.3. AH
 *
 * Revision 1.3  1995/02/08  10:01:14  mringe
 * PL_ entfernt.
 *
 * Revision 1.2  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 1.1  1994/07/22  18:54:34  mringe
 * Initial revision
 *
 */


#include <string.h>
#include <stdlib.h>

#include "meataxe.h"



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static long *piv = NULL;		/* Pivot table */
static long maxnoc = -1;		/* Table size */


/* ------------------------------------------------------------------
   Function prototypes
   ------------------------------------------------------------------ */

static void pivalloc(long noc);


/* ------------------------------------------------------------------
   pivalloc() - Allocate pivot table. All functions in this module
	use one global pivot table (piv). pivalloc() checks first
	if the currently allocated table is large enough.
   ------------------------------------------------------------------ */

static void pivalloc(noc)
long noc;

{
    if (noc <= maxnoc)
	return;
    maxnoc = noc;
    if (piv != NULL)
	free(piv);
    piv = NALLOC(long,noc+1);
    if (piv == NULL)
	FATAL("pivalloc(): out of memory");
}


/* ------------------------------------------------------------------
   echelon() - Echelon form of a matrix (non-destructive).
   echelon_() - Convert a matrix to echelon form.
   ------------------------------------------------------------------ */

matrix_t *echelon_(mat)
matrix_t *mat;

{
    long rank;

    pivalloc(mat->nor);
    zsetfield(mat->fl);
    zsetlen(mat->noc);
    rank = zmkechelon(mat->d,mat->nor,piv);
    mat->nor = rank;
    mat->d = (PTR) Srealloc(mat->d,zsize(rank));
    return mat;
}

matrix_t *echelon(mat)
matrix_t *mat;

{
    matrix_t *tmp;

    if ((tmp = matdup(mat)) == NULL) return NULL;
    return echelon_(tmp);
}


/* ------------------------------------------------------------------
   nullity() - Return the nullity of a matrix, matrix is not changed
   nullity_() - Return the nullity of a matrix, matrix is changed
   nullity__() - Return the nullity of a matrix, matrix is deleted
   ------------------------------------------------------------------ */

long nullity(mat)
matrix_t *mat;

{
    matrix_t *ech = echelon(mat);
    long nul;
    if (ech == NULL) return -1;
    nul = ech->noc - ech->nor;
    matfree(ech);
    return nul;
}

long nullity__(mat)
matrix_t *mat;

{
    long nul;
    echelon_(mat);
    nul = mat->noc - mat->nor;
    matfree(mat);
    return nul;
}

long nullity_(mat)
matrix_t *mat;

{
    echelon_(mat);
    return mat->noc - mat->nor;
}


/* ------------------------------------------------------------------
   nullspace() - Null-space of a matrix.
   nullspace_() - Matrix null-space, puts mat into echelon form.
   nullspace__() - Matrix null-space, deletes mat.
   ------------------------------------------------------------------ */

matrix_t *nullspace_(mat)
matrix_t *mat;

{
    long dim;
    matrix_t *nsp;
    PROFILE_BEGIN(t);

    pivalloc(mat->nor);
    nsp = matalloc(mat->fl,mat->nor,mat->nor);
    if (nsp == NULL) return NULL;
    zsetlen(mat->noc);
    dim = znullsp(mat->d,mat->nor,piv,nsp->d);

    zsetlen(nsp->noc);
    nsp->nor = dim;
    nsp->d = (PTR) Srealloc(nsp->d,zsize(dim));

    zsetlen(mat->noc);
    mat->nor = mat->nor - dim;
    mat->d = (PTR) Srealloc(mat->d,zsize(mat->nor));

    PROFILE_END(t,NullSpace);
    return nsp;
}

matrix_t *nullspace(mat)
matrix_t *mat;

{
    matrix_t *tmp, *nsp;
    if ((tmp = matdup(mat)) == NULL)
	return NULL;
    nsp = nullspace_(tmp);
    matfree(tmp);
    return nsp;
}

matrix_t *nullspace__(mat)
matrix_t *mat;

{
    matrix_t *nsp = nullspace_(mat);
    matfree(mat);
    return nsp;
}


