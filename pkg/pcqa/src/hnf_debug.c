/*
  Written by: Eddie Lo
  Date started: October 24, 1995.

  Part of the polycyclic quotient algorithm package.

  This program contains procedures to help debug the abelian module.
*/

#include "pcqa.h"

extern exponent **row;
extern counter Size_R, Size_S;


/* This procedure prints a row. */
void Print_Row(r, ncol)
exponent *r;
counter ncol;
{
  generator j;

  for (j = 0; j < ncol; j++) printf("%d\t", r[j]);
  printf("\n");
}

/* This procedure prints the relation matrix. */
void Print_Matrix(Matrix, nrow, ncol)
exponent **Matrix;
counter nrow, ncol;
{
  generator i;

  for (i = 0; i < nrow; i++) Print_Row(Matrix[i], ncol);
}

/* This procedure multiplies two matrices to give a third matrix.
   For debugging. */
void Print_Mult_Matrix(M1, r1, c1, M2, r2, c2)
exponent **M1, **M2;
counter r1, c1, r2, c2;
{
  generator i1, i2, j;
  exponent num;

  if (c1 != r2) {
    printf("Incompatible size.\n");
    return;
  };
  for (i1 = 0; i1 < r1; i1++) {
    for (i2 = 0; i2 < c2; i2++) {
      num = 0;
      for (j = 0; j < c1; j++) num += M1[i1][j]*M2[j][i2];
      printf("%ld\t", num);
    };
    printf("\n");
  }
}

