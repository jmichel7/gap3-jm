/*
  Written by: Eddie Lo
  Date started: March 29, 1996.

  Part of the polycyclic quotient algorithm program.

  This module takes as inputs a matrix and computes a Smith normal
  form of it. This module requires the Hermite normal form hnf module.
  A list of module elements, a total of Test_Gen, is required besides
  the matrix.
*/

#include <stdlib.h>
#include "pcqa.h"

exponent **Smith_Conv = NULL;
exponent **Smith_Inv = NULL;
exponent *col = NULL;
Poly **Smith_Base = NULL;

extern generator NumMod;


/* This procedure initializes for the Smith normal form module. */
void Init_Smith(nrow, ncol)
counter nrow, ncol;
{
  counter k;

  Smith_Conv = (exponent **) calloc(ncol, sizeof(exponent *));
  Smith_Conv[0] = (exponent *) calloc(ncol*ncol, sizeof(exponent));
  for (k = 1; k < ncol; k++) Smith_Conv[k] = Smith_Conv[k-1]+ncol;
  for (k = ncol*ncol-1; k > 0; k--) Smith_Conv[0][k] = 0;
  for (k = 0; k < ncol; k++) Smith_Conv[k][k] = 1;
  Smith_Inv = (exponent **) calloc(ncol, sizeof(exponent *));
  Smith_Inv[0] = (exponent *) calloc(ncol*ncol, sizeof(exponent));
  for (k = 1; k < ncol; k++) Smith_Inv[k] = Smith_Inv[k-1]+ncol;
  for (k = ncol*ncol-1; k > 0; k--) Smith_Inv[0][k] = 0;
  for (k = 0; k < ncol; k++) Smith_Inv[k][k] = 1;
  if (nrow > ncol) col = (exponent *) calloc(nrow, sizeof(exponent));
  else col = (exponent *) calloc(ncol, sizeof(exponent));
}

/* This procedure resets memory for Smith normal form module. */
void Reset_Smith()
{
  free(Smith_Conv[0]);
  free(Smith_Conv);
  free(Smith_Inv[0]);
  free(Smith_Inv);
  free(col);
}

/* This procedure initializes for a new basis which corresponds
   to the Smith normal form. */
void Init_SBase(ncol)
counter ncol;
{
  generator i;

  Smith_Base = (Poly **) calloc(ncol, sizeof(Poly *));
  Smith_Base[0] = (Poly *) calloc(ncol*NumMod, sizeof(Poly));
  for (i = 1; i < ncol; i++) Smith_Base[i] = Smith_Base[i-1]+NumMod;
}

/* This procedure resets the basis initialized for the Smith normal
   form. */
void Reset_SBase()
{
  generator i;

  free(Smith_Base[0]);
  free(Smith_Base);
}

/* This procedure computes the Smith normal form given a matrix. */
void Compute_Smith(Matrix, nrow, ncol)
exponent **Matrix;
counter nrow, ncol;
{
  flag test;
  generator i, ci, cj, ck, r;
  exponent s, s1, s2, s3, s4;

  Abelian_Row_Reduce(Matrix, nrow, ncol);
  test = Test_Single(Matrix, nrow, ncol, &ci, &cj, &r);
  while (test == NO) {
    for (ck = cj+1; ck < ncol; ck++) {
      if (Matrix[r][ck] != 0 && abs(Matrix[r][ck]) < abs(Matrix[r][cj]))
        cj = ck;
    };
    s = Extended_GCD(Matrix[r][ci], abs(Matrix[r][cj]), &s1, &s2);
    s3 = Matrix[r][cj]/s; s4 = Matrix[r][ci]/s;
    if (Matrix[r][cj] < 0) s2 = -s2;
    for (i = 0; i < ncol; i++) {
      col[i] = s1*Smith_Conv[i][ci]+s2*Smith_Conv[i][cj];
      Smith_Conv[i][cj] = -s3*Smith_Conv[i][ci]+s4*Smith_Conv[i][cj];
    };
    for (i = 0; i < ncol; i++) Smith_Conv[i][ci] = col[i];
    for (i = 0; i < ncol; i++) {
      col[i] = s4*Smith_Inv[ci][i]+s3*Smith_Inv[cj][i];
      Smith_Inv[cj][i] = -s2*Smith_Inv[ci][i]+s1*Smith_Inv[cj][i];
    };
    for (i = 0; i < ncol; i++) Smith_Inv[ci][i] = col[i];
    for (i = 0; i < nrow; i++) {
      col[i] = s1*Matrix[i][ci]+s2*Matrix[i][cj];
      Matrix[i][cj] = -s3*Matrix[i][ci]+s4*Matrix[i][cj];
    };
    for (i = 0; i < nrow; i++) Matrix[i][ci] = col[i];
    Abelian_Row_Reduce(Matrix, nrow, ncol);
    test = Test_Single(Matrix, nrow, ncol, &ci, &cj, &r);
  };
  Align_Matrix(Matrix, nrow, ncol);
  test = Test_Divide(Matrix, nrow, ncol, &ci, &cj, &r);
  while (test == NO) {
    s = Extended_GCD(Matrix[r-1][ci], Matrix[r][cj], &s1, &s2);
    s3 = Matrix[r-1][ci]/s; s4 = Matrix[r][cj]/s;
    for (i = 0; i < ncol; i++) {
      col[i] = s3*Smith_Conv[i][ci]-s4*Smith_Conv[i][cj];
      Smith_Conv[i][cj] = s2*Smith_Conv[i][ci]+s1*Smith_Conv[i][cj];
    };
    for (i = 0; i < ncol; i++) Smith_Conv[i][ci] = col[i];
    for (i = 0; i < ncol; i++) {
      col[i] = s1*Smith_Inv[ci][i]+s4*Smith_Inv[cj][i];
      Smith_Inv[cj][i] = -s2*Smith_Inv[ci][i]+s3*Smith_Inv[cj][i];
    };
    for (i = 0; i < ncol; i++) Smith_Inv[ci][i] = col[i];
    Matrix[r][cj] = Matrix[r][cj]/s*Matrix[r-1][ci];
    Matrix[r-1][ci] = s;
    test = Test_Divide(Matrix, nrow, ncol, &ci, &cj, &r);
  }
}

/* This procedure tests whether a matrix has two nonzero entries in the
   same row. If yes, it outputs the two columns with the nonzero entries
   together with the row they are in. This procedure assumes that Matrix
   is in Hermite normal form. */
flag Test_Single(Matrix, nrow, ncol, cx, cy, r)
exponent **Matrix;
counter nrow, ncol;
generator *cx, *cy, *r;
{
  generator ri, ci, cj;

  for (ri = 0; ri < nrow; ri++) {
    for (ci = 0; ci < ncol && Matrix[ri][ci] == 0; ci++);
    if (ci == ncol) return YES;
    else {
      for (cj = ci+1; cj < ncol && Matrix[ri][cj] == 0; cj++);
      if (cj != ncol) {
        *cx = ci;
        *cy = cj;
        *r = ri;
        return NO;
      }
    }
  };
  return YES;
}

/* This procedure tests whether a matrix is in Smith normal form. If
   not, it outputs the two columns and the first row in which the
   hypothesis is not satisfied. This procedure assumes that no two
   nonzero entries are in the same row. This procedure also assumes
   that Matrix is in Hermite normal form. */
flag Test_Divide(Matrix, nrow, ncol, cx, cy, r)
exponent **Matrix;
counter nrow, ncol;
generator *cx, *cy, *r;
{
  exponent num;
  generator ci, cj, ri;

  num = 0;
  cj = 0;
  for (ri = 0; ri < nrow; ri++) {
    while (cj < ncol && Matrix[ri][cj] == 0) cj++;
    if (cj == ncol) return YES;
    else if (num == 0) num = Matrix[ri][cj];
    else if (Matrix[ri][cj] % num != 0) {
      *cx = ci;
      *cy = cj;
      *r = ri;
      return NO;
    }
    else {
      num = Matrix[ri][cj];
      ci = cj;
      cj++;
    }
  };
  return YES;
}

/* This procedure multiplies a row by a matrix. Assume that nrow is
   at least ncol. */
void Mult_Row_Matrix(row, Matrix, nrow, ncol)
exponent *row, **Matrix;
counter nrow, ncol;
{
  generator c, r;
  exponent sum;

  for (c = 0; c < ncol; c++) {
    col[c] = 0;
    for (r = 0; r < nrow; r++) col[c] += row[r]*Matrix[r][c];
  };
  for (c = 0; c < ncol; c++) row[c] = col[c];
}

/* This procedure aligns the rows and columns in the matrix, so it looks
   like a diagonal form. */
void Align_Matrix(Matrix, nrow, ncol)
exponent **Matrix;
counter nrow, ncol;
{
  generator i, j, ci;
  exponent num;

  for (i = 0; i < nrow; i++) {
    for (ci = i; ci < ncol && Matrix[i][ci] == 0; ci++);
    if (ci == ncol) return;
    if (ci != i) {
      Matrix[i][i] = Matrix[i][ci];
      Matrix[i][ci] = 0;
      for (j = 0; j < ncol; j++) {
        num = Smith_Conv[j][i];
        Smith_Conv[j][i] = Smith_Conv[j][ci];
        Smith_Conv[j][ci] = num;
        num = Smith_Inv[i][j];
        Smith_Inv[i][j] = Smith_Inv[ci][j];
        Smith_Inv[ci][j] = num;
      }
    }
  }
}

