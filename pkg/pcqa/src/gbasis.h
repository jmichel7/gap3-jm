/*
  Written by: Eddie Lo
  Date started: December 9, 1995.

  Part of the polycyclic quotient algorithm package.

  This is the header file for the polynomial package.
*/


/* In file gbasis.c. */

counter Criterion(counter n1, counter n2, element u1, element u2, 
  flag delete, generator start);
counter Save_Basis_Mod(Poly *m, generator start);
counter Save_Mod(Poly *m, generator start);
flag Alpha_Mod(HTSpt hpt1, HTSpt hpt2, Poly *m, generator start);
flag Compute_Basis();
flag Confirm_Basis();
flag Erase_Computation();
flag Minimum_Mod(Poly *m, generator start, element u, coeff num);
flag Mod_Reduce_Mod(Poly *m1, Poly *m2);
flag Next_Pair(counter *n1, counter *n2, HTSpt *hpt1, HTSpt *hpt2);
flag Process_Mod(Poly *m, generator *start, flag mode);
flag Process_Pairs();
flag Reduce_Mod(Poly *m, HTSpt hpt, generator *start, flag pos);
void Add_Basis();
void Add_Rule(Rule_Set *R, counter k);
void Basis_Reduce_Mod(Poly *m, generator *start);
void Create_Pairs(counter last);
void Form_Basis();
void Form_Univ(coeff num);
void Groebner_Basis(Poly **Initial_List, counter Initial_Num);
void HT_Compute(HTSpt hpt1, HTSpt hpt2);
void Initial_Pairs();
void Minimum_Reduce(Poly *m, generator *start);
void Mutual_Reduce();
void Nice_Full_Reduce_Mod(Poly *m, generator *start);
void Randomize_Set(counter last);
void Recover_Basis();
void Set_Full_Reduce_Mod(Poly *m, generator *start);
void Set_Reduce_Mod(Poly *m, generator *start);
void Short_Full_Reduce(Poly *m, generator start);
void Sort_Mod();


/* In file gbasis_basic.c. */

Poly Get_Poly(counter n);
ModSpt Get_ModStorage(counter n);
counter Get_NumTerm(counter n);
HTSpt Get_HT(counter n);
HTSpt Get_Basis_HT(counter n);
void Info_HT(counter n, element *lead, element *pre, coeff *num, ModSpt *mpt);
void Info_Basis_HT(counter n, element *lead, element *pre, coeff *num,
  ModSpt *mpt);
counter Save_HT(element lead, element pre, ModSpt mpt, coeff num);
HTSpt Store_HT_Data(element lead, element pre, ModSpt mpt, coeff num);
counter Save_Basis_HT_Mod(Poly *m, generator start);
void Assign_HT(counter n, HTSpt hpt);
void Assign_Basis_HT(counter n, HTSpt hpt);
void Delete_HT(counter n, flag erase);
void Delete_Basis_HT(counter n, flag erase);
void Clear_HT();
void Clear_Basis_HT();


/* In file gbasis_io.c. */

flag Check_Basis(flag ind);
void Save_Basis(FILE *fo);
void Save_Temp_Basis(FILE *fo);
void Read_Basis(FILE *fi);
void Print_Mod_HT(Poly *m, generator start);
void Print_Original(HTSpt hpt);
void Print_All_HT();
void Print_All_Basis_HT();
void Save_Build(Poly *m, generator start);


/* In file gbasis_debug.c. */

void Print_Imm_Stack();
void Print_Late_Stack();
counter List_Memory();
void Print_Memory();
void Print_Stack();
void Print_HT(HTSpt hpt);


/* In file gbasis_mem.c. */

void Est_Tool();
void Init_Tool();
void Reset_Tool();
void Est_Memory();
void Init_Memory();
void Reset_Memory();
ModSpt New_Mod();
void Allocate_Mod();
void Free_Mod(ModSpt mpt, flag erase);
HTSpt New_HT();
HTSpt New_Basis_HT();
void Allocate_HT();
void Free_HT(HTSpt hpt);
void Free_Basis_HT(HTSpt hpt);
void Reallocate_Basis(counter n);


/* In file stack.c. */
void Push_Exp_Stack(counter n1, counter n2);
void Push_Imm_Stack(counter n1, counter n2);
void Push_Late_Stack(counter n1, counter n2);
flag Search_Imm_Stack(counter *n1, counter *n2, HTSpt *h1, HTSpt *h2);
flag Search_Exp_Stack(counter *n1, counter *n2, HTSpt *h1, HTSpt *h2);
flag Search_Late_Stack(counter *n1, counter *n2, HTSpt *h1, HTSpt *h2);
void Pop_Exp_Stack(counter *n1, counter *n2, HTSpt *h1, HTSpt *h2);
void Pop_Imm_Stack(counter *n1, counter *n2, HTSpt *h1, HTSpt *h2);
void Pop_Late_Stack(counter *n1, counter *n2, HTSpt *h1, HTSpt *h2);

