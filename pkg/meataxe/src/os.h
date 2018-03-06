/***********************************************************************
 * This file is part of the C Meat-Axe
 * Written by Michael Ringe <mringe@tiffy.math.rwth-aachen.de>
 * (C) Copyright 1992:	    Lehrstuhl D fuer Mathematik
 *                          RWTH Aachen
 *                          Aachen, Germany
 ***********************************************************************
 * $Source: /gap/CVS/GAP/3.4.4/pkg/meataxe/src/os.h,v $
 * $Revision: 1.1.1.1 $
 * $Date: 1996/12/11 12:40:20 $
 ***********************************************************************
 * Header file containing OS dependent stuff
 **********************************************************************/

/* ------------------------------------------------------------------
   os_xxx Functions (os.c)
   ------------------------------------------------------------------ */

#define FM_READ 1
#define FM_CREATE 2
#define FM_APPEND 3
#define FM_TEXT 0x10
#define FM_LIB 0x20

void os_init _PL((void));
long os_timeused _PL((void));
void os_timelimit _PL((long nsec));

char *os_mkfilename _PL((char *name));
FILE *os_fopen _PL((char *name, int mode));
int os_fseek _PL((FILE *f,long pos));

