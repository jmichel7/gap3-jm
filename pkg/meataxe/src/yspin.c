/* ========================== C MeatAxe =============================
   yspin.c - Spin up

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: yspin.c,v 1.1.1.1 1996/12/11 12:40:18 werner Exp $
 *
 * $Log: yspin.c,v $
 * Revision 1.1.1.1  1996/12/11 12:40:18  werner
 * Preparing 3.4.4 for release
 *
 * Revision 2.2  1993/12/13  08:25:53  mringe
 * Reihenfolge der Fkt.-argumente vereinheitlicht.
 *
 * Revision 2.1  1993/10/26  08:19:43  mringe
 * Neu: matspin_f().
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.3  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.2  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.1  1993/02/15  13:50:28  mringe
 * Initial revision
 *
 */




#include <stdlib.h>
#include <string.h>
#include "meataxe.h"
#include "lattice.h"


static int dontclean = 0;	/* Do not clean up after spin-up */
static long *piv;		/* Pivot table */
static long subdim;		/* Subspace dimension */


/* ----------------------------------------------------------------
   matspin() - Spin up using n matrices
   ------------------------------------------------------------------ */

matrix_t *matspin(seed,ngen,gen)
matrix_t *seed;
int ngen;
matrix_t *gen[];

{
    matrix_t *sub, *myseed;
    long fl;
    int i;
    PTR mygen[MAXGEN];
    long dim;

    /* Initialize
       ---------- */
    myseed = echelon(seed);
    dim = seed->noc;
    fl = seed->fl;
    for (i = 0; i < ngen; ++i)
    {
	if (gen[i]->noc != dim || gen[i]->nor != dim ||
	    gen[i]->fl != fl)
	    FATAL("MATRICES INCOMPATIBLE");
	mygen[i] = gen[i]->d;
    }
    piv = (long *)malloc((size_t)(dim + 1) * sizeof(long));
    sub = matalloc(fl,dim,dim);	/* Worst case */

    /* Spin up
       ------- */
    zsetlen(fl,dim);
    memcpy(sub->d,myseed->d,zsize(myseed->nor));
    subdim = zspinup(sub->d,myseed->nor,piv,mygen,ngen,T_MATRIX);
    sub->nor = subdim;
    if (!dontclean)
    {
    	free(piv);
    	sub->d = (PTR) realloc(sub->d,zsize(subdim));
    }
    return sub;
}




/* ----------------------------------------------------------------
   matspin_f() - Same as matspin(), but resulting matrix is filled
   out with extra rows to make a nonsingular square matrix. Returns
   the dimension of the subspace in sdim.
   ------------------------------------------------------------------ */

matrix_t *matspin_f(seed,ngen,gen,sdim)
matrix_t *seed;
int ngen;
matrix_t *gen[];
long *sdim;

{
    matrix_t *result;
    char *flag;
    PTR x;
    long i;

    dontclean = 1;
    result = matspin(seed,ngen,gen);
    dontclean = 0;
    *sdim = subdim;
    flag = (char *) malloc(seed->noc+1);
    memset(flag,0,seed->noc+1);
    for (i = 1; i <= subdim; ++i) flag[piv[i]] = 1;
    x = result->d;
    zadvance(&x,subdim);
    for (i = 1; i <= subdim; ++i)
    {
	if (flag[piv[i]] == 0)
	{
	    zmulrow(x,F_ZERO);
	    zinsert(x,i,F_ONE);
	    zadvance(&x,1);
	}
    }
    result->nor = seed->noc;
    free(piv);
    free(flag);
    return result;
}

