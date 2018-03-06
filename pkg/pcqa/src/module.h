/*
  Written by: Eddie Lo
  Date started: February 25, 1996.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the module module.
*/


/* In file module.c. */

void Init_Many_Mod();
void Reset_Mod();
generator Non_Null_Mod(Poly *m);
void Mult_Mod(Poly *m_ans, Poly *m, element u, coeff num, generator start);
generator Merge_Mod(Poly *m_ans, Poly *m1, Poly *m2, generator start1,
  generator start2);
void Negate_Mod(Poly *m, generator start);
void Copy_Mod(Poly *m_ans, Poly *m, generator start);
void Print_Mod(FILE *print_file, Poly *m, generator start);
void Read_Mod(FILE *input_file, Poly *m, generator *start);
counter NumTerm_Mod(Poly *m, generator start);
void Memory_Mod(Poly *m, counter *tterm, counter *tcell, counter *tint,
  generator start);
void Free_Mod_Cell(Poly *m, generator start);
void Mult_Poly_Mod(Poly *m_ans, Poly f, Poly *m, generator start);
void Mult_Mod_Mod(Poly *m_ans, Poly *m1, Poly *m2, generator start);
coeff Find_MHTP(Poly *m, element elm, generator *start);
coeff BigTerm_Mod(Poly *m, generator start);
void Output_Mod(FILE *fo, Poly *m, generator start);
