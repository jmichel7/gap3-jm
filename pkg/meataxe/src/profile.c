/* ========================== C MeatAxe =============================
   profile.c - Global data for profiling

   (C) Copyright 1994 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: profile.c,v 1.2 1997/09/11 15:43:16 gap Exp $
 *
 * $Log: profile.c,v $
 * Revision 1.2  1997/09/11 15:43:16  gap
 * New version 2.2.3. AH
 *
 * Revision 1.1  1994/07/23  16:53:04  mringe
 * Initial revision
 *
 * Revision 1.1  1994/07/23  16:53:04  mringe
 * Initial revision
 *
 */


#include "meataxe.h"


/* ------------------------------------------------------------------
   Variables
   ------------------------------------------------------------------ */

long ProfFileIO = 0;		/* File I/O */
long ProfNullSpace = 0;		/* Matrix null space */
long ProfSpinUp = 0;		/* Spin-up */
long ProfMakeWord = 0;		/* Make word */
long ProfMatInsert = 0;		/* Matrix insertion */
long ProfPolFactor = 0;		/* Polynomial factorization */
long ProfCharPol = 0;		/* Char. and minimal polynomial */

