/* ========================== C MeatAxe =============================
   zgap.c -  Convert between MeatAxe and GAP format.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zgap.c,v 1.2 1997/09/11 15:43:52 gap Exp $
 *
 * $Log: zgap.c,v $
 * Revision 1.2  1997/09/11 15:43:52  gap
 * New version 2.2.3. AH
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.6  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.5  1993/05/12  08:10:45  mringe
 * Prototypes korrigiert.
 *
 * Revision 1.4  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.3  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.2  1993/01/06  21:11:59  mringe
 * *** empty log message ***
 *
 * Revision 1.1  1993/01/06  20:31:50  mringe
 * Initial revision
 *
 */


#include "meataxe.h"



/* ------------------------------------------------------------------
   zftogap() - Convert field element into GAP readable text
	string.
   ------------------------------------------------------------------ */

#if defined (__STDC__)
char *zftogap(FEL f)
#else
char *zftogap(f)
FEL f;
#endif

{
    static char buffer[40];

    if (zchar == zfl)	/* Prime field */
    {
	FEL f2 = F_ZERO;
   	long k = 0;
    	while (f2 != f)
    	{   f2 = zadd(f2,zgen);
	    ++k;
	}
	sprintf(buffer,"%ld*Z(%ld)",k,zfl);
    }
    else		/* Other field */
    {
	if (f == F_ZERO)
	    sprintf(buffer,"0*Z(%ld)",zfl);
	else
	{
	    FEL f2 = zgen;
	    long k = 1;
	    while (f2 != f)
	    {   f2 = zmul(f2,zgen);
		++k;
	    }
	    sprintf(buffer,"Z(%ld)^%ld",zfl,k);
	}
    }
    return buffer;
}






