/* ========================== C MeatAxe =============================
   charpol.c -  Characteristic polynomial of a matrix

   (C) Copyright 1994 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: charpol.c,v 1.2 1997/09/11 15:42:40 gap Exp $
 *
 * $Log: charpol.c,v $
 * Revision 1.2  1997/09/11 15:42:40  gap
 * New version 2.2.3. AH
 *
 * Revision 1.9  1995/02/08  10:02:58  mringe
 * PL_ entfernt
 *
 * Revision 1.8  1994/09/14  09:59:38  mringe
 * Neue Variable: CharPolSeed.
 *
 * Revision 1.7  1994/07/29  07:43:36  mringe
 * Dimensionslimit in charpol() entfernt.
 *
 * Revision 1.6  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 1.5  1994/07/23  16:49:35  mringe
 * Profiling
 *
 * Revision 1.4  1994/07/07  05:52:21  mringe
 * Compiler warning
 *
 * Revision 1.3  1994/07/04  11:47:36  mringe
 * Dimensionlimit bei charpol()
 *
 * Revision 1.2  1994/05/04  11:42:29  mringe
 * hmmm.
 *
 * Revision 1.2  1994/05/04  11:42:29  mringe
 * hmmm.
 *
 * Revision 1.1  1994/04/07  21:11:02  mringe
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

long CharPolSeed = 1;		/* Seed */

static long fl, nor;		/* Field and size */
static long *piv = NULL;	/* Pivot table */
static char *ispiv = NULL;	/* Pivot flags */
static PTR mat = NULL,		/* The matrix */
    A = NULL,			/* Work space (for spin-up) */
    B = NULL;			/* Work space II (coefficients) */
static long dim;		/* Dimension reached so far */
static long n;			/* Dimension of cyclic subspace */



/* ------------------------------------------------------------------
   cleanup() - Free memory
   ------------------------------------------------------------------ */

static void cleanup()

{
    if (mat != NULL) free(mat);
    if (A != NULL) free(A);
    if (B != NULL) free(B);
    if (piv != NULL) free(piv);
    if (ispiv != NULL) free(ispiv);
    mat = A = B = NULL;
    piv = NULL;
    ispiv = NULL;
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
    dim = 0;
    if (CharPolSeed < 1 | CharPolSeed > nor)
	CharPolSeed = 1;

    /* Allocate memory
       --------------- */
    cleanup();
    zsetfield(fl);
    zsetlen(nor);
    mat = zalloc(nor);
    A = zalloc(nor+1);
    B = zalloc(nor+1);
    piv = (long *) malloc((size_t)(nor+2)*sizeof(long));
    ispiv = (char *) malloc((size_t)(nor+2));
    if (piv == NULL || ispiv == NULL) MTXFAIL(ERR_NOMEM,-1);

    /* Initialize memory
       ----------------- */
    memcpy(mat,matrix->d,zsize(nor));
    memset(ispiv,0,(size_t)(nor+2));

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
    
    pol = polalloc(fl,nor);
    x = B;
    zadvance(&x,n);
    pol->deg = n;
    for (k = 1; k <= n; ++k)
	pol->buf[k-1] = zextract(x,k);
    pol->buf[n] = F_ONE;
    return pol;
}


/* ------------------------------------------------------------------
   spinup_cyclic() - Spin up one cyclic subaspace
   ------------------------------------------------------------------ */

static void spinup_cyclic()

{   PTR a, b;
    long pv, k;
    FEL f;

/*fprintf(stderr,"spinup_cyclic(): dim=%ld\n",dim);*/

    a = A;
    zadvance(&a,dim);
    b = B;
    zmulrow(b,F_ZERO);
    n = 0;
    while ((pv = zfindpiv(a,&f)) != 0)
    {	
	PTR x, y;

    	/* Add new vector to basis
	   ----------------------- */
	piv[dim+1+n] = pv;
/*fprintf(stderr,"piv[%ld]=%ld\n",dim+1+n,pv);*/
	ispiv[pv] = 1;
	++n;
	zinsert(b,n,F_ONE);

	/* Calculate the next vector
	   ------------------------- */
	x = a;
	zadvance(&a,(long)1);
	zmaprow(x,mat,nor,a);
/*
{long ii; fprintf(stderr," x ="); for (ii=1; ii < nor; ++ii)
 fprintf(stderr,"%ld",zftoi(zextract(x,ii))); fprintf(stderr,"\n xA=");
 for (ii=1; ii < nor; ++ii) fprintf(stderr,"%ld",zftoi(zextract(a,ii)));
 fprintf(stderr,"\n");
}
*/
	y = b;
	zadvance(&b,(long)1);
	zmulrow(b,F_ZERO);
	for (k = 1; k < nor; ++k)
	    zinsert(b,k+1,zextract(y,k));
    
	/* Clean with existing basis vectors
	   --------------------------------- */
	x = A;
	y = B;
	for (k = 1; k <= dim+n; ++k)
	{   f = zdiv(zextract(a,piv[k]),zextract(x,piv[k]));
	    zaddmulrow(a,x,zneg(f));
	    if (k > dim)
	    {	zaddmulrow(b,y,zneg(f));
		zadvance(&y,(long)1);
	    }
	    zadvance(&x,(long)1);
	}
    }
    dim += n;
}


/* ------------------------------------------------------------------
   charpolfactor() - Returns next (reducible) factor
   ------------------------------------------------------------------ */

poly_t *charpolfactor(matrix_t *mat)

{
    PTR seed;
    int i;
    PROFILE_BEGIN(t);

    /* If called with mat != NULL, initialize everything
       ------------------------------------------------- */
    if (mat != NULL)
    	if (init(mat) != 0)
    	    return NULL;

    /* If there is nothing left to do, return NULL
       ------------------------------------------- */
    /* we could call cleanup() here... */
    if (dim >= nor) return NULL;

    /* Prepare the next seed vector
       ---------------------------- */
    zsetfield(fl);
    zsetlen(nor);
    seed = A;
    zadvance(&seed,dim);
    if (dim == 0)
	i = CharPolSeed;
    else
    	for (i = 1; i <= nor && ispiv[i] != 0; ++i);
    zmulrow(seed,F_ZERO);
    zinsert(seed,i,F_ONE);

    /* Spin up and return the polynomial 
       --------------------------------- */
    spinup_cyclic();
    PROFILE_END(t,CharPol);
    return mkpoly();    
}



/* ------------------------------------------------------------------
   charpol() - Returns factorized characteristic polynomial.
   ------------------------------------------------------------------ */

fpoly_t *charpol(matrix_t *mat)

{
    fpoly_t *cpol = fpolalloc();
    poly_t *p;

    for (p = charpolfactor(mat); p != NULL; p = charpolfactor(NULL))
    {
	fpoly_t *fac = factorization(p);
    	polfree(p);
	fpolmul(cpol,fac);
	fpolfree(fac);
    }
    return cpol;
}
