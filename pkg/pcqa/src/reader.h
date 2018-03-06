/*
  Written by: Eddie Lo
  Date started: February 25, 1996.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the reader module.
*/


/* In file reader.c. */

void Init_Many_Reader();
void Init_Reader();
void Reset_Reader();
char Real_Char();
char Input_Char();
void ReaderError(char *str);
void Symbol();
coeff Coefficient();
exponent Positive();
Poly Block();
Poly Pterm();
Poly Fterm();
Poly Product();
Poly Polynomial();
void Melement(Poly *m);
void Read_Mod_Grammar(FILE *fi, Poly *m);

