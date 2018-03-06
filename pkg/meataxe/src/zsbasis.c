/* ========================== C MeatAxe =============================
   Standard basis.

   (C) Copyright 1994 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zsbasis.c,v 1.2 1997/09/11 15:44:24 gap Exp $
 *
 * --- COMPLETELY REWRITTEN  for 2.2.2 ---
 *
 */


#include <stdlib.h>
#include "meataxe.h"


int zsbasis(PTR seed, long nseed, int ngen, PTR *gen, PTR space,
	long *piv, PTR basis)

{
    long i, k, seedcount;
    PTR xi, yi, xk, yk;
    FEL f;
    int igen;

    /* Initialize
       ---------- */
    i = 1; xi = space; yi = basis;
    k = 1; xk = space; yk = basis;
    seedcount = 1;
    igen = 0;

    /* Main loop
       --------- */
    while (1)
    {
	if (i >= k)	/* Next seed vector */
	{
	    if (seedcount > nseed) break;
	    zmoverow(yk,seed);
	    zmoverow(xk,seed);
	    zadvance(&seed,(long)1);
	    ++seedcount;
	}
	else		/* Apply next generator */
	{
	    zmaprow(yi,gen[igen],znoc,yk);
	    zmoverow(xk,yk);
	    if (++igen >= ngen)	/* Alle Erzeuger benutzt */
	    {
		igen = 0;
		++i;
	    	zadvance(&xi,(long)1);
	    	zadvance(&yi,(long)1);
	    }
	}

	/* Clean and check if we got a new vector
	   -------------------------------------- */
	zcleanrow(xk,space,k-1,piv);
	if ((piv[k] = zfindpiv(xk,&f)) != 0)
	{
	    if (k == znoc) break;		/* Finished */
	    ++k;
	    zadvance(&xk,(long)1);
	    zadvance(&yk,(long)1);
	}
    }
    if (k < znoc) return 2;		/* Failed */
    if (seedcount > 2) return 1;	/* Used more than one vector */
    return 0;				/* Success */
}


