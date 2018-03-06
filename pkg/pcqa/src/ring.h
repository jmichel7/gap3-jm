/*
  Written by: Eddie Lo
  Date started: August 23,1994.

  Part of the polycyclic quotient algorithm program.

  This is the header file of the ring module.
*/


/* In file ring.c. */

cellptr Copy_List(cellptr p);
cellptr Insert_List(cellptr head, cellptr obj);
coeff BigCoeff_Poly(Poly f);
counter NumTerm_Poly(Poly f);
counter Short_Poly(Poly f, counter n);
counter Sym_Poly(Poly *m, generator start);
flag Full_Reduce(Poly *m, Poly *m_red, generator start, flag pos);
flag Modulo_Reduce_Mod(Poly *m, Poly *m_red, generator start, flag pos);
flag Modulo_Reduce(Poly *f, coeff num);
Poly Constant_Poly(coeff num, flag use);
Poly Erase_Extra_Poly(Poly f);
Poly Merge_Poly(Poly f1, Poly f2);
Poly Monomial_Poly(element u, coeff num, flag use);
Poly Mult_Poly(Poly f, element u, coeff num);
Poly Read_Poly(FILE *input_file);
void Fix_div(coeff num);
void Fix_Zero(cellptr head);
void Free_Poly_Cell(Poly f);
void Interactive_Check(cellptr p);
void Memory_Poly(Poly f, counter *numterm, counter *numcell, counter *memory);
void Negate_Poly(Poly f);
void Print_Poly(FILE *print_file, Poly f);
void Show_Poly(FILE *print_file, Poly f, flag show);


/* In file ring_adv.c. */

Poly Mult_Poly_Poly(Poly f, Poly g);
Poly Power_Poly(Poly f, exponent n);
Poly Aug_Poly(element u);
Poly Split_Poly(Poly *f);
coeff First_Term_Poly(Poly *f, element u);
Poly Class_Poly(generator *perm, generator limit);


/* In file ring_basic.c. */

void Est_Basic();
void Init_Many_Basic();
void Reset_Basic();
void Allocate_Cell();
cellptr New_Cell();
void Free_Cell(cellptr p);
void Insert_INT(coeff p);
void Allocate_MP_INT();
coeff New_MP_INT();
void Free_MP_INT(coeff p);
void Print_MP_INT(coeff p);
void Print_Basic();


/* In file ring_mem.c. */

void Est_Structure();
void Init_Many_Structure();
void Reset_Structure();
void Init_Structure();
