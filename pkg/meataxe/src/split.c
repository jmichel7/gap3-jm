/* ========================== C MeatAxe =============================
   Split a representation.

   (C) Copyright 1994 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: split.c,v 1.2 1997/09/11 15:43:24 gap Exp $
 *
 * $Log: split.c,v $
 * Revision 1.2  1997/09/11 15:43:24  gap
 * New version 2.2.3. AH
 *
 * Revision 2.10  1994/07/29  07:42:36  mringe
 * Fehlermeldung verbessert.
 *
 * Revision 2.9  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.8  1994/07/07  07:13:41  mringe
 * split_pivot ist global.
 *
 * Revision 2.7  1994/05/18  10:09:12  mringe
 * Funktionen neu geschrieben und auf zwei Files verteilt
 * (spinup.c und split.c).
 *
 * Revision 2.6  1994/05/11  18:23:49  mringe
 * VOELLIG NEUE VERSION!!!
 *
 * Revision 2.5  1994/05/06  15:23:30  mringe
 * Neu: tr=2 (make closure).
 *
 * Revision 2.4  1994/04/09  13:51:23  mringe
 * memset()-Bug ???
 *
 * Revision 2.3  1993/12/13  08:25:53  mringe
 * Reihenfolge der Fkt.-argumente vereinheitlicht.
 *
 * Revision 2.2  1993/12/07  14:34:24  mringe
 * Teste bei realloc() auf size==0.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.12  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.11  1993/02/12  17:06:59  mringe
 * Woerter mit N Erzeugern.
 *
 * Revision 1.10  1993/02/12  08:31:41  mringe
 * Benutze Library-Funktionen fuer sun und quot.
 *
 * Revision 1.9  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.8  1993/01/27  13:05:02  mringe
 * Benutze F_ZERO und F_ONE statt 0 und 1.
 *
 * Revision 1.7  1993/01/15  07:31:42  mringe
 * try_split() auf N Erzeuger umgestellt.
 *
 * Revision 1.6  1993/01/15  06:48:03  mringe
 * sbasis() auf N Erzeuger umgestellt.
 *
 * Revision 1.5  1993/01/09  14:50:26  mringe
 * Berechne Projektion des seed space auf quot.
 *
 * Revision 1.4  1992/10/02  15:50:28  mringe
 * Benutze F_ZERO statt 0.
 *
 * Revision 1.3  1992/10/01  13:50:10  mringe
 * Header eingef"ugt.
 *
 * Revision 1.2  1992/07/22  07:10:30  mringe
 * Changed 'global.h' to 'lattice.h'
 *
 * Revision 1.1  1992/05/24  08:29:47  mringe
 * Initial revision
 *
 */



#include <stdlib.h>
#include "meataxe.h"


long *split_pivot = NULL;	/* Pivot table */
static long maxpiv;


/* ------------------------------------------------------------------
   split() - Calculate the action of the generators on subspace
  	and quotient. 
   ------------------------------------------------------------------ */

int split(subspace,ngen,gen,sub,quot)
matrix_t *subspace;
int ngen;
matrix_t *gen[];
matrix_t **sub;
matrix_t **quot;

{
    int g;
    long i;
    PTR xi, yi;
    long fl = subspace->fl;
    long dim = subspace->noc;
    long sdim = subspace->nor;
    long qdim = dim - sdim;
    PTR mat = subspace->d;
    PTR tmp;

    zsetfield(fl);
    zsetlen(dim);
    tmp = zalloc((long)1);

    /* Make pivot table. Echelonize the subspace if necessary.
       ------------------------------------------------------- */
    if (split_pivot == NULL || maxpiv < dim)
    {
	if (split_pivot != NULL) free(split_pivot);
	split_pivot = NALLOC(long,dim+1);
	maxpiv = dim;
    }
    if (zmkpivot(subspace->d,subspace->nor,split_pivot) != 0)
    {
	sdim = zmkechelon(subspace->d,sdim,split_pivot);
	subspace->nor = sdim;
	qdim = dim - sdim;
    }

    /* Subspace
       -------- */
    for (g = 0; g < ngen; ++g)
    {
	sub[g] = matalloc(fl,sdim,sdim);
	xi = mat;
	yi = sub[g]->d;
	for (i = 1; i <= sdim; ++i)
	{
	    FEL f;
	    zsetlen(sdim);
	    zmulrow(yi,F_ZERO);
	    zsetlen(dim);
	    zmaprow(xi,gen[g]->d,dim,tmp);
	    zadvance(&xi,(long)1);
	    zcleanrow2(tmp,mat,sdim,split_pivot,yi);
	    if (zfindpiv(tmp,&f) != 0)
		FATAL("split(): Subspace not invariant");
	    zsetlen(sdim);
	    zadvance(&yi,(long)1);
	}
    }

    /* Quotient
       -------- */
    if (quot != NULL)
    {
    	zsetlen(dim);
    	zquotinit(mat,sdim,split_pivot);
    	for (g = 0; g < ngen; ++g)
    	{
    	    quot[g] = matalloc(fl,qdim,qdim);
    	    zquotop(gen[g]->d,quot[g]->d);
    	}
    }

    zsetlen(dim);
    free(tmp);
    return 0;
}




/* ------------------------------------------------------------------
   quotproj() - Calculate the projection of vectors onto the quotient.
   ------------------------------------------------------------------ */

matrix_t *quotproj(subspace,vectors)
matrix_t *subspace;
matrix_t *vectors;

{
    long fl = subspace->fl;
    long dim = subspace->noc;
    long sdim = subspace->nor;
    long qdim = dim - sdim;
    PTR mat = subspace->d;
    matrix_t *result;


    if (subspace->fl != vectors->fl || subspace->noc != vectors->noc)
    {
	MTXFAIL(ERR_INCOMPAT,NULL);
    }

    /* Initialize the zquot routines
       ----------------------------- */
    if (split_pivot == NULL || maxpiv < dim)
    {
	if (split_pivot != NULL) free(split_pivot);
	split_pivot = NALLOC(long,dim+1);
	maxpiv = dim;
    }
    zsetfield(fl);
    zsetlen(dim);
    if (zmkpivot(subspace->d,subspace->nor,split_pivot) < 0)
    {
	sdim = zmkechelon(subspace->d,subspace->nor,split_pivot);
	subspace->nor = sdim;
	qdim = dim - sdim;
    }

    zsetlen(dim);
    zquotinit(mat,sdim,split_pivot);
    
    /* Calculate the projection and reduce to echelon form
       --------------------------------------------------- */
    result = matalloc(fl,vectors->nor,qdim);
    zsetlen(dim);
    zquot(vectors->d,vectors->nor,result->d);
    zsetlen(qdim);
    result->nor = zmkechelon(result->d,result->nor,split_pivot);
    zsetlen(dim);
    return result;
}


