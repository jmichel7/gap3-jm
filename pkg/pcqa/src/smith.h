/*
  Written by: Eddie Lo
  Date started: March 29, 1996.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the Smith normal form module.
*/


/* In file smith.c. */

void Init_Smith(counter nrow, counter ncol);
void Reset_Smith();
void Init_SBase(counter ncol);
void Reset_SBase();
void Compute_Smith(exponent **Matrix, counter nrow, counter ncol);
flag Test_Single(exponent **Matrix, counter nrow, counter ncol,
                 generator *cx, generator *cy, generator *r);
flag Test_Divide(exponent **Matrix, counter nrow, counter ncol,
                 generator *cx, generator *cy, generator *r);
void Mult_Row_Matrix(exponent *row, exponent **Matrix,
                     counter nrow, counter ncol);
void Align_Matrix(exponent **Matrix, counter nrow, counter ncol);

