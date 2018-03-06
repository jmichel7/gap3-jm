/* ========================== C MeatAxe =============================
   Standard basis, matrix_t version.

   (C) Copyright 1994 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: sbasis.c,v 1.2 1997/09/11 15:43:21 gap Exp $
 *
 * $Log: sbasis.c,v $
 * Revision 1.2  1997/09/11 15:43:21  gap
 * New version 2.2.3. AH
 *
 * Revision 1.1  1994/05/18  10:08:41  mringe
 * Initial revision
 *
 * Revision 1.1  1994/05/18  10:08:41  mringe
 * Initial revision
 *
 */


#include <stdlib.h>
#include "meataxe.h"


/* ------------------------------------------------------------------
   sbasis() - Spin up `canonically'
   ------------------------------------------------------------------ */

matrix_t *sbasis(seed, ngen, gen)
matrix_t *seed;		/* A seed vector */
int ngen;		/* Number of generators */
matrix_t *gen[];	/* The generators */

{
    matrix_t *basis;
    PTR *genptr;
    long fl = gen[0]->fl;
    long dim = gen[0]->nor;
    int i, result;
    long *piv;
    PTR workspace;

    /* Check arguments
       --------------- */
    for (i = 0; i < ngen; ++i)
    {
	if (gen[i]->nor != dim || gen[i]->noc != dim)
	    FATAL("sbasis(): generators incompatible or not square");
    }
    if (seed->nor < 1 || seed->noc != dim)
	FATAL("sbasis(): bad seed");

    /* Allocate space
       -------------- */
    basis = matalloc(fl,dim,dim);
    workspace = zalloc(dim);
    piv = NALLOC(long,dim+1);
    genptr = NALLOC(PTR,ngen);
    for (i = 0; i < ngen; ++i) genptr[i] = gen[i]->d;


    /* Spin up
       ------- */
    result = zsbasis(seed->d,(long)1,ngen,genptr,workspace,piv,
	basis->d);
    if (result != 0) FATAL("sbasis(): not irreducible");
    free(workspace);
    free(piv);
    free(genptr);
    return basis;
}

