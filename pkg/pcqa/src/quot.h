/*
  Written by: Eddie Lo
  Date started: February 6, 1996.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the quotient module.
*/


/* In file quotient.c. */

void Consistency_Poly();
void Consistency_Homom();
counter Correct_Forward_Map();
generator Correct_Backward_Map();
generator Onto_Map();
generator Q_Reduce(element row);
void Q_Mutual_Reduce(element r1, element r2, generator i);
void Full_Matrix();
void Sift_Add(element row);
void Consistency_Collect(element u);
exponent Extended_GCD(exponent m1, exponent m2, exponent *r1, exponent *r2);


/* In file quot_io.c. */

flag Quotient_Control();
void Quotient_Menu();
void Quotient_Run();
void Print_Arr();


/* In file quot_mem.c. */

void Init_Quotient();
void Reset_Quotient();

