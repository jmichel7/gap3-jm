/***********************************************************************
 * This file is part of the C Meat-Axe                                
 * Written by Michael Ringe <mringe@tiffy.math.rwth-aachen.de>
 * (C) Copyright 1992:	    Lehrstuhl D fuer Mathematik
 *                          RWTH Aachen
 *                          Aachen, Germany
 ***********************************************************************
 * $Source: /gap/CVS/GAP/3.4.4/pkg/meataxe/src/words.h,v $
 * $Revision: 1.1.1.1 $
 * $Date: 1996/12/11 12:40:19 $
 ***********************************************************************
 * Header file for words.c
 **********************************************************************/


/* ------------------------------------------------------------------
   Constants
   ------------------------------------------------------------------ */

#define MAXFP 6
#define MAXBASIS 50


/* ------------------------------------------------------------------
   Typedefs
   ------------------------------------------------------------------ */

typedef struct
{
    matrix_t *b[MAXBASIS];
    int ngen;
    matrix_t *w;
    matrix_t *z4, *z5, *z7, *z9, *z10, *z11, *z12;
}
basis_t;


/* ------------------------------------------------------------------
   Functions prototypes
   ------------------------------------------------------------------ */

long nextword _PL((long w));
void initbasis _PL((matrix_t *gen[], int ngen, basis_t *basis));
void freebasis _PL((basis_t *basis));
void mkword _PL((basis_t *b, long n));
void makefp _PL((basis_t *basis, long fp[]));
char *nameof _PL((basis_t *b,long n));

