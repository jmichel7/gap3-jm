/*
  Written by: Eddie Lo
  Date started: January 20, 1996.

  Part of the polycyclic quotient algorithm package.

  This program consists of procedures to compute Hermite normal form
  given a matrix.
*/

#include "pcqa.h"

extern exponent *expterm;
extern generator *nzterm, *corr;


/* This procedure mutually reduces two rows r1 and r2, where i is
   the index where both r1[i], r2[i] are nonzero, and r1[j], r2[j]
   are 0 for all j < i. */
void Abelian_Mutual_Reduce(ncol, r1, r2, i)
counter ncol;
exponent *r1, *r2;
generator i;
{
  generator j;
  exponent *s1, *s2, *temp;
  int quot;

  s1 = r2; s2 = r1;
  while (s2[i] != 0) {
    if (s1[i] >= 0 || s2[i] == 1 || s2[i] == -1) quot = s1[i]/s2[i];
    else quot = (s1[i]-1)/s2[i];
    for (j = i; j < ncol; j++) s1[j] -= s2[j]*quot;
    temp = s1; s1 = s2; s2 = temp;
  }
}

/* This procedure row reduces a matrix. */
void Abelian_Row_Reduce(Matrix, nrow, ncol)
exponent **Matrix;
counter nrow, ncol;
{
  generator ri, rj, rk, ci, cj, i, j;
  exponent *temp;
  int quot;

  ri = ci = 0;
  while (ri < nrow && ci < ncol) {

/* Assumption: Matrix[rl][cl] = 0 for all rl >= ri, cl < ci. */
    for (rj = ri; rj < nrow && Matrix[rj][ci] == 0; rj++);
    if (rj < nrow) {
      temp = Matrix[ri]; Matrix[ri] = Matrix[rj]; Matrix[rj] = temp;
      for (rk = ri+1; rk < nrow; rk++) {
        if (Matrix[rk][ci] != 0) {
          Abelian_Mutual_Reduce(ncol, Matrix[ri], Matrix[rk], ci);
          if (Matrix[rk][ci] != 0) {
            temp = Matrix[ri]; Matrix[ri] = Matrix[rk]; Matrix[rk] = temp;
          }
        }
      };
      if (Matrix[ri][ci] < 0)
        for (i = ci; i < ncol; i++) Matrix[ri][i] = -Matrix[ri][i];
      for (i = 0; i < ri; i++) {
        quot = Matrix[i][ci]/Matrix[ri][ci];
        if (Matrix[i][ci] < 0 && Matrix[i][ci]%Matrix[ri][ci]) quot--;
        for (j = ci; j < ncol; j++) Matrix[i][j] -= quot*Matrix[ri][j];
      };
      ri++; ci++;
    }
    else ci++;
  }
}

/* This procedure reduces a row by Matrix. */
void Abelian_Reduce(Matrix, nrow, ncol, r)
exponent **Matrix;
counter nrow, ncol;
exponent *r;
{
  generator i, j1, j2;
  int quot;

  for (i = 0; i < nrow && (j1 = nzterm[i]) < ncol; i++) {
    if (r[j1] != 0) {
      quot = r[j1]/Matrix[i][j1];
      if (r[j1] < 0 && Matrix[i][j1]*quot != r[j1]) quot--;
      for (j2 = j1; j2 < ncol; j2++) r[j2] -= quot*Matrix[i][j2];
    }
  }
}

/* This procedure finds the first non-zero term in the row of the relation
   matrix. */
generator Find_misc(Matrix, nrow, ncol)
exponent **Matrix;
counter nrow, ncol;
{
  generator i, j;

  for (i = 0; i < ncol; i++) expterm[i] = 0;
  j = 0;
  for (i = 0; i < nrow; i++) {
    while (j < ncol && Matrix[i][j] == 0) j++;
    if (j != ncol) expterm[j] = Matrix[i][j];
    nzterm[i] = j;
  };
  for (i = 0, j = 1; i < ncol; i++)
    if (expterm[i] != 1) {
      corr[j] = i;
      j++;
    };
  return j-1;
}

/* This procedure converts a ncol vector to a newcol vector, backward. */
void Fill_Vec(newcol, v, u)
counter newcol;
exponent *v;
exponent *u;
{
  generator i, j;

  for (i = 1; i <= newcol; i++) u[-i] = v[corr[i]];
}

