/* ========================== C MeatAxe =============================
   Spin-up and split.

   (C) Copyright 1994 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: spinup.c,v 1.2 1997/09/11 15:43:23 gap Exp $
 *
 * $Log: spinup.c,v $
 * Revision 1.2  1997/09/11 15:43:23  gap
 * New version 2.2.3. AH
 *
 * Revision 1.5  1994/08/22  12:05:56  mringe
 * Bug in do_spin() behoben.
 *
 * Revision 1.5  1994/08/22  12:05:56  mringe
 * Bug in do_spin() behoben.
 *
 * Revision 1.4  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 1.3  1994/06/14  12:05:48  mringe
 * Bug behoben (malloc).
 *
 * Revision 1.2  1994/06/13  13:32:06  mringe
 * MTXFAIL() falsch aufgerufen. Korrigiert.
 *
 * Revision 1.1  1994/05/18  10:09:12  mringe
 * Initial revision
 *
 * Revision 1.1  1994/05/18  10:09:12  mringe
 * Initial revision
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


static int Options;
static matrix_t *Subspace = NULL;	/* Invariant subspace */
static long *spin_pivot = NULL;		/* Pivot table */
static long maxpiv;

static PTR *Gen = NULL;
static int NGen;
static long Fl;
static int Dim;
static long seedno;






static int do_spin(sp,ns)
PTR sp;
long ns;

{
    PTR x = sp, y = Subspace->d;
    long l;
    int success;

    if (ns > Dim) FATAL("do_spin(): Overflow");
    for (l = ns; l > 0; --l)
    {
	zmoverow(y,x);
	zadvance(&x,(long)1);
	zadvance(&y,(long)1);
    }
    l = zspinup(Subspace->d,ns,spin_pivot,NGen,Gen,T_MATRIX);
    Subspace->nor = l;
    success = (l > 0 && l < Dim);
    if (success)
	Subspace->d = (PTR) realloc(Subspace->d,zsize(l));
    return success;
}



int spinup(seed,ngen,gen,options,subspace)
matrix_t *seed;
int ngen;
matrix_t *gen[];
int options;
matrix_t **subspace;

{
    /* Step 1: Initialize
       ------------------ */
    zsetfield(seed->fl);
    zsetlen(seed->noc);
    if ((options & 0xF) != SPL_CONTINUE)
    {
	int i, g;
	if (ngen < 0 || seed->nor < 1) MTXFAIL(ERR_BADARG,-1);
	for (i = 0; i < ngen; ++i)
	{
	    if (gen[i]->fl != seed->fl || gen[i]->noc != seed->noc)
		MTXFAIL(ERR_INCOMPAT,-1);
	    if (gen[i]->noc != gen[i]->nor)
		MTXFAIL(ERR_NOTSQUARE,-1);
	}
    	NGen = ngen;
	if (Gen != NULL) free(Gen);
    	Gen = NALLOC(PTR,NGen+1);	/* Avoid malloc(0) */
    	for (g = 0; g < NGen; ++g)
	    Gen[g] = gen[g]->d;
    	Options = options;
	Dim = seed->noc;
	Fl = seed->fl;
        *subspace = Subspace = matalloc(Fl,Dim+1,Dim);
        if (spin_pivot == NULL || maxpiv < Dim)
        {
	    if (spin_pivot != NULL) free(spin_pivot);
	    spin_pivot = NALLOC(long,Dim+1);
	    maxpiv = Dim;
        }
    }


    /* Step 2: Spin up
       --------------- */
    switch (Options & 0xF)
    {
	case SPL_SEED_MAKE:	/* Generate seed vectors */
	    zpseed_init(seed->nor,seed->d);
	    for (seedno = zpseed_next(); seedno >= 0;
		seedno = zpseed_next())
	    {
	    	if (do_spin(zpseed_vec,(long)1)) return 1;
	    }
	    return 0;
    
	case SPL_SEED_EACH:	/* Try each row of seed */
	    for (seedno = 1; seedno <= seed->nor; ++seedno)
	    {
	        PTR x;
	        x = seed->d;
	        zadvance(&x,seedno-1);
	        if (do_spin(x,(long)1)) return 1;
	    }
	    return 0;
    
    	case SPL_SEED_FIRST:	/* Try only the first vector */
	    return do_spin(seed->d,(long)1);
    
	case SPL_SEED_SPACE:	/* Make closure */
	    seed->nor = zmkechelon(seed->d,seed->nor,spin_pivot);
	    return do_spin(seed->d,seed->nor);
    }

    FATAL("Unknown option in spinup()");
    return -1;
}


