/* ========================== C MeatAxe =============================
   minpol.c -  Minimal polynomial of a matrix

   (C) Copyright 1994 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: minpol.c,v 1.2 1997/09/11 15:43:02 gap Exp $
 *
 * $Log: minpol.c,v $
 * Revision 1.2  1997/09/11 15:43:02  gap
 * New version 2.2.3. AH
 *
 * Revision 1.4  1995/02/08  10:14:52  mringe
 * _PL entfernt.
 *
 * Revision 1.4  1995/02/08  10:14:52  mringe
 * _PL entfernt.
 *
 * Revision 1.3  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 1.2  1994/07/23  16:48:07  mringe
 * Funktionierende Version (?)
 *
 * Revision 1.1  1994/07/22  08:15:38  mringe
 * Initial revision
 *
 */

#include <string.h>
#include <stdlib.h>
#include "meataxe.h"




/* ------------------------------------------------------------------
   Function prototypes
   ------------------------------------------------------------------ */

static void cleanup(void);
static int init(matrix_t *matrix);
static void spinup_cyclic(void);
static poly_t *mkpoly(void);

/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static long fl, nor;		/* Field and size */
static long *APiv = NULL;	/* Pivot table (A) */
static long *CPiv = NULL;	/* Pivot table (C) */
static char *CIsPiv = NULL;	/* Pivot flags */
static PTR mat = NULL,		/* The matrix */
    A = NULL,			/* Work space (spin-up) */
    B = NULL,			/* Work space II (coefficients) */
    C = NULL,			/* Work space III (basis) */
    seed = NULL;
static long CDim;		/* Dimension reached so far */
static long ADim;
static poly_t *mpol = NULL;

/* ------------------------------------------------------------------
   cleanup() - Free memory
   ------------------------------------------------------------------ */

static void cleanup()

{
    if (mat != NULL) free(mat);
    if (A != NULL) free(A);
    if (B != NULL) free(B);
    if (C != NULL) free(C);
    if (seed != NULL) free(seed);
    if (APiv != NULL) free(APiv);
    if (CPiv != NULL) free(CPiv);
    if (CIsPiv != NULL) free(CIsPiv);
    if (mpol != NULL) polfree(mpol);
    mpol = NULL;
    mat = A = B = C = NULL;
    CPiv = APiv = NULL;
    CIsPiv = NULL;
}


/* ------------------------------------------------------------------
   init()
   ------------------------------------------------------------------ */

static int init(matrix)
matrix_t *matrix;

{	
    /* Set some global variables
       ------------------------- */
    if (matrix->nor != matrix->noc) MTXFAIL(ERR_NOTSQUARE,-1);
    fl = matrix->fl;
    nor = matrix->nor;
    CDim = 0;

    /* Allocate memory
       --------------- */
    cleanup();
    zsetfield(fl);
    zsetlen(nor);
    mat = zalloc(nor);
    A = zalloc(nor+1);	/* We need 1 extra row !! */
    B = zalloc(nor+1);
    C = zalloc(nor+1);
    seed = zalloc(1);
    mpol = polalloc(fl,0);
    APiv = NALLOC(long,nor+2);
    CPiv = NALLOC(long,nor+2);
    CIsPiv = NALLOC(char,nor+2);
    if (CPiv == NULL || CIsPiv == NULL)
	MTXFAIL(ERR_NOMEM,-1);

    /* Initialize memory
       ----------------- */
    memcpy(mat,matrix->d,zsize(nor));
    memset(CIsPiv,0,(size_t)(nor+2));

    return 0;
}


/* ------------------------------------------------------------------
   mkpoly() - Make polynomial for the latest cyclic subspace
   ------------------------------------------------------------------ */

static poly_t *mkpoly()

{
    int k;
    PTR x;
    poly_t *pol;
    
    pol = polalloc(fl,ADim);
    x = B;
    zadvance(&x,ADim);
    for (k = 1; k <= ADim; ++k)
	pol->buf[k-1] = zextract(x,k);
    pol->buf[ADim] = F_ONE;
    return pol;
}


/* ------------------------------------------------------------------
   spinup_cyclic() - Spin up one cyclic subaspace
   ------------------------------------------------------------------ */

static void spinup_cyclic()

{   PTR a, b, c;
    long k, pv;
    FEL f;

    a = A;
    b = B;
    c = C;
    zadvance(&c,CDim);
    zmoverow(a,seed);
    zmulrow(b,F_ZERO);
    zinsert(b,1,F_ONE);
    ADim = 0;
    while ((pv = zfindpiv(a,&f)) != 0)
    {	
	PTR x, y;

    	/* Add new vector to basis A
	   ------------------------- */
	zmoverow(c,a);
	APiv[++ADim] = pv;
	zadvance(&a,(long)1);
	zadvance(&b,(long)1);

    	/* Add new vector to basis C
	   ------------------------- */
	zcleanrow(c,C,CDim,CPiv);
	if ((pv = zfindpiv(c,&f)) != 0)
	{
	    CPiv[++CDim] = pv;
	    CIsPiv[pv] = 1;
	    zadvance(&c,(long)1);
	}

	/* Calculate the next vector
	   ------------------------- */
	zmaprow(seed,mat,nor,a);
	zmoverow(seed,a);
	zmulrow(b,F_ZERO);
	zinsert(b,ADim+1,F_ONE);

	/* Clean with existing basis vectors
	   --------------------------------- */
	x = A; y = B;
	for (k = 1; k <= ADim; ++k)
	{   f = zneg(zdiv(zextract(a,APiv[k]),
		zextract(x,APiv[k])));
	    zaddmulrow(a,x,f);
	    zaddmulrow(b,y,f);
	    zadvance(&x,(long)1);
	    zadvance(&y,(long)1);
	}
    }
}

/* ------------------------------------------------------------------
   minpolfactor() - Spin up one cyclic constituent and
   calculate the next factor of m(x), resp.
   ------------------------------------------------------------------ */

poly_t *minpolfactor(mat)
matrix_t *mat;

{
    int i;
    poly_t *ggt, *l, *h;
    PROFILE_BEGIN(t);

    /* If called with mat != NULL, initialize everything
       ------------------------------------------------- */
    if (mat != NULL)
    	if (init(mat) != 0)
    	    return NULL;

    /* If there is nothing left to do, return NULL
       ------------------------------------------- */
    /* we could call cleanup() here... */
    if (CDim >= nor) return NULL;

    /* Prepare the next seed vector
       ---------------------------- */
    zsetfield(fl);
    zsetlen(nor);
    for (i = 1; i <= nor && CIsPiv[i] != 0; ++i);
    if (i > nor) FATAL("internal error");
    zmulrow(seed,F_ZERO);
    zinsert(seed,i,F_ONE);

    /* Spin up and return the polynomial 
       --------------------------------- */
    spinup_cyclic();
    h = mkpoly();    
    ggt = polgcd(h,mpol);
    l = poldivmod(h,ggt);
    polfree(h);
    polfree(ggt);
    polmul(mpol,l);
    PROFILE_END(t,CharPol);
    return l;
}


