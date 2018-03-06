/* ========================== C MeatAxe =============================
   znullsp.c - Calculate null-space

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: znullsp.c,v 1.2 1997/09/11 15:44:06 gap Exp $
 *
 * $Log: znullsp.c,v $
 * Revision 1.2  1997/09/11 15:44:06  gap
 * New version 2.2.3. AH
 *
 * Revision 1.5  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 1.4  1994/07/22  18:52:59  mringe
 * znullsp() liefert jetzt den fertigen Nullraum.
 *
 * Revision 1.3  1994/05/29  11:01:32  mringe
 * Bug behoben: fehlendes zsetlen() eingefuegt.
 *
 * Revision 1.2  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 1.1  1993/10/19  09:30:59  mringe
 * Initial revision
 *
 *
 */



#include <stdlib.h>
#include "meataxe.h"



long znullsp(matrix, nor, piv, nsp)
PTR matrix, nsp;
long nor, *piv;

{
    PTR x, y, a, b;
    int i;
    long noc = znoc;
    long dim;
    FEL f;
    PROFILE_BEGIN(t);

    /* Make identity matrix in nsp
       --------------------------- */
    zsetlen(nor);
    x = nsp;
    for (i = 1; i <= nor; ++i)
    {
	piv[i] = 0;
	zmulrow(x,F_ZERO);
	zinsert(x,i,F_ONE);
	zadvance(&x,(long)1);
    }

    /* Gaussian elimination
       -------------------- */
    x = matrix;
    y = nsp;
    for (i = 1; i <= nor; ++i)
    {
	PTR xx = matrix, yy = nsp;
	long k, p;

	for (k = 1; k < i; ++k)
	{
	    if ((p = piv[k]) != 0 && (f = zextract(x,p)) != F_ZERO)
	    {
		f = zneg(zdiv(f,zextract(xx,p)));
		zsetlen(noc);
		zaddmulrow(x,xx,f);
		zsetlen(nor);
		zaddmulrow(y,yy,f);
	    }
	    zsetlen(noc);
	    zadvance(&xx,(long)1);;
	    zsetlen(nor);
	    zadvance(&yy,(long)1);
	}
	zsetlen(noc);
	piv[i] = p = zfindpiv(x,&f);
	zsetlen(noc);
	zadvance(&x,(long)1);
	zsetlen(nor);
	zadvance(&y,(long)1);
    }

    /* Step 2: Reduce the null space to echelon form.
       ---------------------------------------------- */
    dim = 0;
    x = y = nsp;
    a = b = matrix;
    for (i = 1; i <= nor; ++i)
    {
	if (piv[i] == 0)
	{
	    zsetlen(nor);
	    if (y != x) zmoverow(y,x);
	    zcleanrow(y,nsp,dim,piv);
	    piv[++dim] = zfindpiv(y,&f);
	    zadvance(&y,1);
	}
	else
	{
	    zsetlen(noc);
	    if (b != a) zmoverow(b,a);
	    zadvance(&b,1);
	}
	zsetlen(nor);
	zadvance(&x,1);
	zsetlen(noc);
	zadvance(&a,1);
    }
    piv[0] = dim;
    PROFILE_END(t,NullSpace);
    return dim;
}

