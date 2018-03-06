/*
  Written by: Eddie Lo
  Date started: January 20, 1996.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the Hermite normal form module.
*/


/* In file hnf.c. */

void Abelian_Mutual_Reduce(counter n, exponent *r1, exponent *r2, generator i);
void Abelian_Row_Reduce(exponent **M, counter nrow, counter ncol);
generator Find_misc(exponent **M, counter nrow, counter ncol);
void Abelian_Reduce(exponent **M, counter nrow, counter ncol, exponent *r);
void Fill_Vec(counter n, exponent *v, exponent *u);


/* In file hnf_mem.c. */

void Init_HNF(counter nrow, counter ncol);
void Reset_HNF();


/* In file hnf_debug.c. */

void Print_Row(exponent *r, counter ncol);
void Print_Matrix(exponent **Matrix, counter nrow, counter ncol);
void Print_Mult_Matrix(exponent **M1, counter r1, counter c1,
                       exponent **M2, counter r2, counter c2);

