/*
  Written by: Eddie Lo
  Date started: January 20, 1996.

  Part of the polycyclic quotient algorithm package.

  This module handles memory for the Hermite normal form package.
*/

#include <stdlib.h>
#include "pcqa.h"

exponent *expterm = NULL;
generator *nzterm = NULL;
generator *corr = NULL;


/* This procedure initializes for the hnf module. */
void Init_HNF(nrow, ncol)
counter nrow, ncol;
{
  expterm = (exponent *) calloc(ncol, sizeof(exponent));
  nzterm = (generator *) calloc(nrow, sizeof(generator));
  corr = (generator *) calloc(ncol+1, sizeof(generator));
}

/* This procedure resets for the hnf module. */
void Reset_HNF()
{
  free(expterm); free(nzterm); free(corr);
}

