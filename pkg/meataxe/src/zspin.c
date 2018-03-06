/* ========================== C MeatAxe =============================
   zspin.c - Spin-up of vectors and related functions.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zspin.c,v 1.2 1997/09/11 15:44:30 gap Exp $
 *
 * $Log: zspin.c,v $
 * Revision 1.2  1997/09/11 15:44:30  gap
 * New version 2.2.3. AH
 *
 * Revision 2.6  1994/11/28  16:46:52  mringe
 * Type cast in zpermrow().
 *
 * Revision 2.5  1994/11/28  16:45:09  mringe
 * ANSI-C
 *
 * Revision 2.4  1994/07/23  16:51:53  mringe
 * Profiling.
 *
 * Revision 2.3  1994/05/17  06:42:58  mringe
 * Reihenfolge der Argumente.
 *
 * Revision 2.2  1994/04/09  13:51:23  mringe
 * memset()-Bug ???
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.4  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.3  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.2  1993/05/12  10:31:47  mringe
 * Prototypes korrigiert.
 *
 * Revision 1.1  1993/02/10  19:40:54  mringe
 * Initial revision
 *
 * Revision 1.1  1993/01/28  07:35:51  mringe
 * Initial revision
 *
 */


#include "meataxe.h"

/* ------------------------------------------------------------------
   zspinup() - Spin up one or more vectors using a set of matrices
   or permutations.
   
   Arguments: `space' is the seed space where seed vectors are
   taken from, `nseed' is the number of rows in `seed'. `piv' must
   point to a pivot table for the seed space. `ngen' is the number of
   generators, `gen' is the array of generators, and `gentype' is
   either T_MARTEIX of T_PERM.

   Return value: Returns the subspace dimension or -1 on error.
   ------------------------------------------------------------------ */

long zspinup(PTR space, long nseed, long *piv, int ngen, PTR gen[],
	int gentype)


{
    PTR xget, xput;
    long iget;
    long dim;
    int igen;
    PROFILE_BEGIN(t);

    if (zmkechelon(space,nseed,piv) == -1) return -1;
    dim = piv[0];
    iget = 1;
    xput = xget = space;
    zadvance(&xput,dim);    
    igen = 0;

    if (dim > 0 && ngen > 0)
    {
        while (iget <= dim && dim < znoc)
        {
	    FEL f;


	    if (gentype == T_MATRIX)
	        zmaprow(xget,(PTR) gen[igen],znoc,xput);
	    else
	        zpermrow(xget,(long *) gen[igen],xput);
	    if (++igen >= ngen)
	    {
	        igen = 0;
	        ++iget;
	        zadvance(&xget,(long)1);
	    }
	    zcleanrow(xput,space,dim,piv);
	    if ((piv[dim+1] = zfindpiv(xput,&f)) != 0)
	    {
	        ++dim;
	        zadvance(&xput,1);
	    }
        }
    }
    piv[0] = dim;
    PROFILE_END(t,SpinUp);
    return dim;
}




