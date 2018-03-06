/***********************************************************************
 * This file is part of the C Meat-Axe                                
 * Written by Michael Ringe <mringe@tiffy.math.rwth-aachen.de>
 * (C) Copyright 1992:	    Lehrstuhl D fuer Mathematik
 *                          RWTH Aachen
 *                          Aachen, Germany
 ***********************************************************************
 * $Source: /gap/CVS/GAP/3.4.4/pkg/meataxe/src/split.h,v $
 * $Revision: 1.1.1.1 $
 * $Date: 1996/12/11 12:40:19 $
 ***********************************************************************
 * Header file for "split.c".
 **********************************************************************/


/* ------------------------------------------------------------------
   Constants
   ------------------------------------------------------------------ */

#define MAXSEED 10	/* Max. dimension of the seed space */


/* ------------------------------------------------------------------
   Types
   ------------------------------------------------------------------ */

typedef struct
	{	int result;		/* 0 = not split, 1 = split */
		matrix_t *sub[MAXGEN];	/* Subspace */
		matrix_t *quot[MAXGEN];	/* Quotient */
		matrix_t *proj;
	}
	split_t;	/* Returned by split() */


/* ------------------------------------------------------------------
   Function prototypes
   ------------------------------------------------------------------ */


split_t *split _PL((matrix_t *seed, int ngen, matrix_t *gen[], int tr));
void splfree _PL((split_t *s));
matrix_t *sbasis _PL((matrix_t *seed, int ngen, matrix_t *gen[]));
matrix_t *mkseed _PL((matrix_t *subsp));


