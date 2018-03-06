/* ========================== C MeatAxe =============================
   fpoly.c - Factorized polynomials

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: fpoly.c,v 1.2 1997/09/11 15:42:50 gap Exp $
 *
 * $Log: fpoly.c,v $
 * Revision 1.2  1997/09/11 15:42:50  gap
 * New version 2.2.3. AH
 *
 * Revision 1.6  1994/06/29  07:49:58  mringe
 * Neu: fpoldup().
 *
 * Revision 1.6  1994/06/29  07:49:58  mringe
 * Neu: fpoldup().
 *
 * Revision 1.5  1994/06/25  13:35:30  mringe
 * fpolread() und fpolwrite() (Dummies).
 *
 * Revision 1.4  1994/05/20  09:09:32  mringe
 * Benutze Smalloc() und Srealloc().
 *
 * Revision 1.3  1994/05/04  07:43:44  mringe
 * Neu: fpolprint().
 *
 * Revision 1.2  1994/04/12  14:01:17  mringe
 * Sortieren nach Vfh wieder rausgeworfen.
 *
 * Revision 1.1  1994/04/07  21:11:02  mringe
 * Initial revision
 *
 */


#include <string.h>
#include <stdlib.h>
#include "meataxe.h"


#define BLOCKSIZE 5	/* Granularity */

/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

/* ------------------------------------------------------------------
   fpoldesc() - Describe
   ------------------------------------------------------------------ */
    
char *fpoldesc(p)
fpoly_t *p;

{
    static char buf[50];
    sprintf(buf,"factorized polynomial (%ld factors)",
	(long) p->len);
    return buf;
}


/* ------------------------------------------------------------------
   fpolalloc() - Allocate a polynomial
   fpolfree() - Free a polynomial
   ------------------------------------------------------------------ */

fpoly_t *fpolalloc()

{
    fpoly_t *x;

    x = ALLOC(fpoly_t);
    if (x != NULL)
    {
    	x->id = T_FPOLY;
    	x->len = 0;
    	x->max = BLOCKSIZE;
    	x->p = NALLOC(poly_t *,BLOCKSIZE);
   	if (x->p != NULL)
	{
	    x->e = NALLOC(long,BLOCKSIZE);
	    if (x->e != NULL)
	    {
		return x;
	    }
	    free(x->p);
	}
	free(x);
    }
    MTXFAIL(ERR_NOMEM,NULL);
}


void fpolfree(x)
fpoly_t *x;

{
    int i;
    for (i = 0; i < x->len; ++i)
	polfree(x->p[i]);
    free(x->p);
    free(x->e);
    free(x);
}



/* ------------------------------------------------------------------
   fpoldup() - Duplicate
   ------------------------------------------------------------------ */

fpoly_t *fpoldup(s)
fpoly_t *s;

{
    fpoly_t *x;

    x = ALLOC(fpoly_t);
    if (x != NULL)
    {
    	x->id = T_FPOLY;
    	x->len = x->max = s->len;
    	x->p = NALLOC(poly_t *,s->len);
   	if (x->p != NULL)
	{
	    x->e = NALLOC(long,s->len);
	    if (x->e != NULL)
	    {
		int i;
		for (i = 0; i < s->len; ++i)
		{
		    x->e[i] = s->e[i];
		    if ((x->p[i] = poldup(s->p[i])) == NULL)
			break;
		}
		if (i >= s->len)
		    return x;		/* Done */

		/* Clean up after error */
		while (--i > 0) polfree(x->p[i]);
	    }
	    free(x->p);
	}
	free(x);
    }
    MTXFAIL(ERR_NOMEM,NULL);
}


/* ------------------------------------------------------------------
   fpolmulp() - Multiply with a power of an an irreducible polynomial
   ------------------------------------------------------------------ */

fpoly_t *fpolmulp(dest, src, pwr)
fpoly_t *dest;
poly_t  *src;
long pwr;

{
    int i;

    for (i = 0; i < dest->len && polcmp(dest->p[i],src) < 0; ++i);
    if (i >= dest->len || polcmp(dest->p[i],src))
    {
	int k;
	if (dest->len >= dest->max)	/* extend */
	{
	    size_t newsize = dest->max + BLOCKSIZE;
	    poly_t **x = (poly_t**) Srealloc(dest->p,newsize *
			sizeof(poly_t*));
	    long *e = (long *) Srealloc(dest->e,newsize * sizeof(long));
	    if (e == NULL || x == NULL) MTXFAIL(ERR_NOMEM,NULL);
	    dest->p = x;
	    dest->e = e;
	    dest->max = newsize;
	}
	/* Platz machen fuer neuen Faktor */
	for (k = dest->len; k > i; --k)
	{
	    dest->p[k] = dest->p[k-1];
	    dest->e[k] = dest->e[k-1];
	}
	++dest->len;
	/* Neuen Faktor einfuegen */
	dest->p[i] = poldup(src);
	dest->e[i] = pwr;
    }
    else
    {
	dest->e[i] += pwr;		/* Increase multiplicity */
    }
    return dest;
}




/* ------------------------------------------------------------------
   fpolmul() - Multiply
   ------------------------------------------------------------------ */

fpoly_t *fpolmul(dest, src)
fpoly_t *dest;
fpoly_t  *src;

{
    int i;

    for (i = 0; i < src->len; ++i)
	fpolmulp(dest,src->p[i],src->e[i]);
    return dest;
}


/* -----------------------------------------------------------------
   fpolprint() - Print to stdout
   ----------------------------------------------------------------- */

void fpolprint(name,p)
char *name;
fpoly_t *p;

{
    int i;

    if (name != NULL) printf("%s =",name);
    for (i = 0; i < p->len; ++i)
    {   
	if (i > 0) printf("      * ");
	printf("(");
	polprint(NULL,p->p[i]);
	printf(")^%ld\n",p->e[i]);
    }
}

/* -----------------------------------------------------------------
   fpolread() - Read
   ----------------------------------------------------------------- */

fpoly_t *fpolread(f)
FILE *f;

{
    return NULL;
}


/* -----------------------------------------------------------------
   fpolwrite() - Write
   ----------------------------------------------------------------- */

int fpolwrite(f,p)
FILE *f;
fpoly_t *p;

{
    return -1;
}

