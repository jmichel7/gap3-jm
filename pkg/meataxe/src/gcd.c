/* ========================== C MeatAxe =============================
   gcd.c - Greatest common divisor

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: gcd.c,v 1.2 1997/09/11 15:42:53 gap Exp $
 *
 * $Log: gcd.c,v $
 * Revision 1.2  1997/09/11 15:42:53  gap
 * New version 2.2.3. AH
 *
 * Revision 1.2  1994/11/25  13:31:37  mringe
 * Neu: lcm(); ANSI-C
 *
 * Revision 1.1  1993/12/14  22:27:11  mringe
 * Initial revision
 *
 * Revision 1.1  1993/12/14  22:27:11  mringe
 * Initial revision
 *
 */

#include "meataxe.h"


/* ------------------------------------------------------------------
   gcd() - Greatest common divisor.
   ------------------------------------------------------------------ */

long gcd(long a,long b)

{
    if (a == 0) return b;
    if (b == 0) return a;
    while (1)
    {
	if ((a %= b) == 0) return b;
	if ((b %= a) == 0) return a;
    }
}


/* ------------------------------------------------------------------
   lcm() - Least common multiple
   ------------------------------------------------------------------ */

long lcm(long a,long b)

{
    register long ggt = gcd(a,b);
    if (ggt == 0) return 0;
    return (a / ggt) * b;
}
