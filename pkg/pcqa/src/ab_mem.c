/*
  Written by: Eddie Lo
  Date started: October 24, 1995.

  Part of the polycyclic quotient algorithm package.

  This program manages memory for the abelian quotient module.
*/

#include <stdlib.h>
#include "pcqa.h"

extern counter Size_R, Size_S;
exponent *test_elm, **row, *rowmem;

/* This procedure initializes for the abelian module. */
void Init_Abelian()
{
  counter k, l;

  row = (exponent **) calloc(Size_R, sizeof(exponent *));
  rowmem = (exponent *) calloc(Size_R*Size_S, sizeof(exponent));
  for (k = 0; k < Size_R; k++) {
    row[k] = rowmem+k*Size_S;
    for (l = 0; l < Size_S; l++) row[k][l] = 0;
  };
  test_elm = (exponent *) calloc(Size_S, sizeof(exponent));
}

/* This procedure frees up memory for the abelian module. */
void Reset_Abelian()
{
  free(rowmem); free(row); free(test_elm);
}

