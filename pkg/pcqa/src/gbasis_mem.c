/*
  Written by: Eddie Lo
  Date started: October 24, 1995.

  Part of the polycyclic quotient algorithm package.

  This program manages memeory for the reduce module.

Introduction:

  This module manages memory for several data types.
  ModStorage: storage for a module element
  HTStorage:  storage for head terms
  stack:      storage for pairs of aligned head terms

  A ModStorage stores a module element. These module elements are
  candidates for the Groebner basis. To save memory, symmetrization
  of a module element is performed without the symmetrized elements
  actually calculated. For instance, if m is a module element,
  m g_1, m g_2, ..., m g_l is a list of the symmetrized forms of m,
  then only one, say m g_1, will be stored as a ModStorage. The rest
  will be computed when needed.
  There are two pieces of information which are common to m g_1, m g_2,
  ..., m g_l. First, the starting module generator is the same. Also,
  the number of terms is the same. Since some of the m g_i may be
  deleted during Groebner basis calculation, we would also want to
  know how many of the m g_i currently exist. When this number gets
  to zero, we will have the liberty to delete this ModStorage to save
  memory.

  A HTStorage is an extension to ModStorage. Every symmetrized form
  of a module element stored in ModStorage is accounted for by a
  HTStorage. Using the notation as above, to represent an m g_i,
  we need to store its leading term, the term which correspond to
  the leading term in m g_1 and the constant in the leading term.
  For instance, if u is the leading term of m g_i and v is the
  term which corrsponds to u in m g_1, then m g_i = m g_1 v^-1 u
  modulo a sign. So m g_i can be recovered easily by multiplication.

  A stack is a data structure containing pairs of indices in the list
  of head terms, HT_List. Contents of the stacks are the pairs which
  have to be coincided and reduced. In this implementation there are
  three stacks. The first is Exp_Stack, it contains pairs with the
  same head term and in which one coefficient divides another.
  The second is Imm_Stack contains pairs with one head term dominating
  the other and its coefficient a multiple of another. Late_Stack
  contains the remaining pairs. Only aligned pairs are pushed onto
  stack. It is not evident at this stage of experimentation that
  using the Exp_Stack will speed up the computation. In fact, it looks
  like computation is slowed down by the use of Exp_Stack. Imm_Stack,
  however, seems to help computation a lot. The user is given the
  choice to use 1, 2 or all 3 stacks.

Data structure:

  A ModStorage has 5 fields: total number of terms, starting module
  generator, remaining symmetrized forms, and the actual module
  element itself. The final field is a union of the size of polynomial
  and a pointer. In this implementation, since the actual module
  element is actually in a ModStorage, the size of ModStorage is
  variable. So ModStorage is declared only as a pointer to counter.

  A HTStorage has 4 fields: leading term, the term which corresponds
  to the leading term, the constant, and a pointer pointing to a
  ModStorage. Again, the actual terms are stored in HTStorage, so
  its size is a variable. HTStorage is declared only as a pointer
  to exponent. The final field can be used as a pointer via the use of
  a union.

  A stack is just a collection of pairs.

Data description:

  Free_Mod_L: List of free ModStorage
  Free_HT_L: List of free HTStorage
  Free_Frame: Storage for two free frames of pairs

  Exp_Stack, Imm_Stack, Late_Stack: stacks of pairs
  HT_List: List of head terms

  Mod_Chain and HT_Chain chains up all the frames allocated for
  ModStorage and HTStorage. It facilitates resetting memory.
*/

#include <stdlib.h>
#include "pcqa.h"

counter used_mod = 0;
counter used_HT = 0;
counter Num_HT = 0;
counter Num_Basis = 0;
counter Cur_Stack = 0;
counter Cur_HT = 0;
counter Cur_Mod = 0;
counter Cur_Basis = 0;
counter Num_Exp_Pair = 0;
counter Num_Imm_Pair = 0;
counter Num_Late_Pair = 0;
counter Num_HT_Frame = 0;
ModSpt Mod_Chain = NULL;
ModSpt Free_Mod_L = NULL;
HTSpt HT_Chain = NULL;
HTSpt Free_HT_L = NULL;
HTSpt **HT_List = NULL;
HTSpt *Basis_HT = NULL;
stackframe Free_Frame[2];
stack Exp_Stack = NULL;
stack Imm_Stack = NULL;
stack Late_Stack = NULL;
counter Useful_Red = 0;
counter Zero_Red = 0;
counter Delete_Poly = 0;

MP_INT MP_poly[17];
coeff univ_exp = NULL;
coeff univ_tmp = NULL;
element reduce_u = NULL;
element reduce_s = NULL;
coeff reduce_div = NULL;
Poly *reduce_m = NULL;
coeff crit_res = NULL;
coeff crit_gcd = NULL;
coeff crit_div = NULL;
coeff minimum_num = NULL;
coeff minimum_gcd = NULL;
coeff minimum_c1 = NULL;
coeff minimum_c2 = NULL;
coeff minimum_c3 = NULL;
element minimum_s = NULL;
element minimum_u = NULL;
Poly *minimum_m1 = NULL;
Poly *minimum_m2 = NULL;
Poly *minimum_m3 = NULL;
element grob_u1 = NULL;
element grob_u2 = NULL;
element grob_s1 = NULL;
element grob_s2 = NULL;
element grob_lcm = NULL;
flag *grob_flag = NULL;
coeff grob_n = NULL;
coeff grob_res1 = NULL;
coeff grob_res2 = NULL;
coeff grob_gcd = NULL;
coeff grob_leadc1 = NULL;
coeff grob_leadc2 = NULL;
Poly *grob_m1 = NULL;
Poly *grob_m2 = NULL;
Poly *grob_m = NULL;
Poly *full_m = NULL;
element create_u1 = NULL;
element create_u2 = NULL;
element initial_u1 = NULL;
element initial_u2 = NULL;
element form_u = NULL;
Poly *form_m = NULL;
element save_u = NULL;
Poly *save_m = NULL;
Poly *read_m = NULL;
element original_u = NULL;
Poly *original_m = NULL;
MP_INT one, neg_one;

extern generator NumGen, NumGenP, NumMod;
extern counter pair_length1, pair_length2, pair_length3;
extern counter Cur_Cell, Cur_MP_INT;
extern flag ExistsMemory;


/* This procedure initializes the polynomial module, used only once. */
void Est_Tool()
{
  counter index;

  for(index = 0; index < 17; index++) mpz_init(MP_poly+index);

  reduce_div = MP_poly;
  crit_res = MP_poly+1;
  crit_gcd = MP_poly+2;
  crit_div = MP_poly+3;
  grob_n = MP_poly+4;
  grob_res1 = MP_poly+5;
  grob_res2 = MP_poly+6;
  grob_leadc1 = MP_poly+7;
  grob_leadc2 = MP_poly+8;
  grob_gcd = MP_poly+9;
  univ_exp = MP_poly+10;
  univ_tmp = MP_poly+11;
  minimum_num = MP_poly+12;
  minimum_gcd = MP_poly+13;
  minimum_c1 = MP_poly+14;
  minimum_c2 = MP_poly+15;
  minimum_c3 = MP_poly+16;

  mpz_init_set_si(&one, 1);
  mpz_init_set_si(&neg_one, -1);

  Useful_Red = Zero_Red = Delete_Poly = 0;
}

/* This procedure initializes the polynomial module, used many times. */
void Init_Tool()
{
  reduce_u = (element) calloc(NumGen*2, sizeof(exponent))+NumGen;
  reduce_s = reduce_u+NumGen;
  reduce_m = (Poly *) calloc(NumMod, sizeof(Poly));

  minimum_u = (element) calloc(NumGen*2, sizeof(exponent))+NumGen;
  minimum_s = reduce_u+NumGen;
  minimum_m1 = (Poly *) calloc(NumMod*3, sizeof(Poly));
  minimum_m2 = minimum_m1+NumMod;
  minimum_m3 = minimum_m2+NumMod;

  grob_lcm = (element) calloc(NumGen*3, sizeof(exponent))+NumGen;
  grob_s1 = grob_lcm+NumGen;
  grob_s2 = grob_s1+NumGen;
  grob_flag = (flag *) calloc(NumGenP, sizeof(flag));
  grob_m = (Poly *) calloc(NumMod*3, sizeof(Poly));
  grob_m1 = grob_m+NumMod;
  grob_m2 = grob_m1+NumMod;

  full_m = (Poly *) calloc(NumMod, sizeof(Poly));

  create_u1 = (element) calloc(NumGen*2, sizeof(exponent))+NumGen;
  create_u2 = create_u1+NumGen;

  initial_u1 = (element) calloc(NumGen*2, sizeof(exponent))+NumGen;
  initial_u2 = initial_u1+NumGen;

  form_u = (element) calloc(NumGen, sizeof(exponent))+NumGen;
  form_m = (Poly *) calloc(NumMod, sizeof(Poly));

  save_u = (element) calloc(NumGen, sizeof(exponent))+NumGen;
  save_m = (Poly *) calloc(NumMod, sizeof(Poly));

  read_m = (Poly *) calloc(NumMod, sizeof(Poly));

  original_u = (element) calloc(NumGen, sizeof(exponent))+NumGen;
  original_m = (Poly *) calloc(NumMod, sizeof(Poly));
}

/* This procedure frees up memory in poly module. */
void Reset_Tool()
{
  free(reduce_u-NumGen); free(reduce_m);
  free(minimum_u-NumGen); free(minimum_m1);
  free(grob_lcm-NumGen); free(grob_flag); free(grob_m);
  free(full_m);
  free(create_u1-NumGen);
  free(initial_u1-NumGen);
  free(form_u-NumGen); free(form_m);
  free(save_u-NumGen); free(save_m);
  free(read_m);
  free(original_u-NumGen); free(original_m);
  Useful_Red = Zero_Red = Delete_Poly = 0;
}

/* This procedure initializes for the memory module, used only once. */
void Est_Memory()
{
  used_HT = Num_HT = Num_Basis = used_mod = Num_HT_Frame = 0;
  Cur_Stack = Cur_HT = Cur_Mod = Cur_Basis = 0;
  Cur_Basis = Begin_Basis_Size;
  Num_Exp_Pair = Num_Imm_Pair = Num_Late_Pair = 0;
  Mod_Chain = NULL; Free_Mod_L = NULL;
  HT_Chain = NULL; Free_HT_L = NULL;
  HT_List = (HTSpt **) calloc(HT_List_Size, sizeof(HTSpt *));
  Basis_HT = (HTSpt *) calloc(Begin_Basis_Size, sizeof(HTSpt));
  Exp_Stack = (stack) calloc(Stack_List_Size, sizeof(stackframe));
  Imm_Stack = (stack) calloc(Stack_List_Size, sizeof(stackframe));
  Late_Stack = (stack) calloc(Stack_List_Size, sizeof(stackframe));
  Free_Frame[0] = Free_Frame[1] = NULL;
}

/* This procedure initializes for the memory module, used many times. */
void Init_Memory()
{
  ExistsMemory = YES;
}

/* This procedure frees up memory used in memory module. */
void Reset_Memory()
{
  counter k;
  HTSpt hpt;
  ModSpt mpt;

  for (k = 1; k <= used_HT; k++) {
    mpt = (Get_HT(k))->pm;
    Free_Mod_Cell(mpt->list, mpt->start);
  };
  for (k = 1; k <= Num_Basis; k++) {
    mpt = (Get_Basis_HT(k))->pm;
    Free_Mod_Cell(mpt->list, mpt->start);
  };
  while (Mod_Chain != NULL) {
    mpt = Mod_Chain;
    Mod_Chain = mpt->u.mspt;
    free(mpt->list);
    free(mpt);
  };
  while (HT_Chain != NULL) {
    hpt = HT_Chain;
    HT_Chain = hpt[0].u.hspt;
    free(hpt[0].lead);
    free(hpt);
  };
  for (k = 0; k < Num_HT_Frame; k++) free(HT_List[k]);
  used_HT = Num_HT = Num_Basis = used_mod = Num_HT_Frame = 0;
  Cur_HT = Cur_Mod = 0;
  Free_Mod_L = NULL;
  Free_HT_L = NULL;
  ExistsMemory = NO;
}

/* This procedure gets a new polynomial from the free list.
   Require global variable: Free_Mod_L, the list of free module elements. */
ModSpt New_Mod()
{
  ModSpt mpt;

  if (Free_Mod_L == NULL) Allocate_Mod();
  mpt = Free_Mod_L;
  Free_Mod_L = mpt->u.mspt;
  used_mod++;
  return mpt;
}

/* This procedure allocates more polynomial storages to the free list.
   Require global variable: Free_Mod_L, the list of free polynomials,
                            Mod_Chain, chain of module storage,
                            Cur_Mod, count of number of frames. */
void Allocate_Mod()
{
  ModSpt mpt;
  Poly *ppt;
  counter index;

  mpt = (ModSpt) calloc(Mod_Frame_Size, sizeof(ModStorage));
  mpt[0].u.mspt = Mod_Chain;
  mpt[0].list = (Poly *) calloc((Mod_Frame_Size-1)*NumMod, sizeof(Poly));
  Mod_Chain = &(mpt[0]);
  Cur_Mod++;
  mpt[1].u.mspt = Free_Mod_L;
  mpt[1].list = mpt[0].list;
  Free_Mod_L = &(mpt[Mod_Frame_Size-1]);
  for (index = 2; index < Mod_Frame_Size; index++) {
    mpt[index].u.mspt = &(mpt[index-1]);
    mpt[index].list = mpt[index-1].list+NumMod;
  }
}

/* This procedure frees a polynomial to the free list.
   Require global variable: Free_Mod_L, the list of free cells. */
void Free_Mod(mpt, erase)
ModSpt mpt;
flag erase;
{
  Poly *m;
  generator index;

  if (erase == YES) {
    m = mpt->list;
    for (index = mpt->start; index < NumMod; index++) Free_Poly_Cell(m[index]);
  };
  mpt->u.mspt = Free_Mod_L;
  Free_Mod_L = mpt;
  used_mod--;
}

/* This procedure returns a free head term to the calling procedure.
   Require global variable: Free_HT_L, list of free head terms.
                            used_HT, number of head terms in use. */
HTSpt New_HT()
{
  HTSpt hpt;

  if (Free_HT_L == NULL) Allocate_HT();
  hpt = Free_HT_L;
  Free_HT_L = hpt->u.hspt;
  used_HT++;
  return hpt;
}

/* This procedure returns a free head term for basis.
   Require global variable: Free_HT_L, list of free head terms. */
HTSpt New_Basis_HT()
{
  HTSpt hpt;

  if (Free_HT_L == NULL) Allocate_HT();
  hpt = Free_HT_L;
  Free_HT_L = hpt->u.hspt;
  return hpt;
}

/* This procedure allocates memory for head term elements.
   Require global variable: HT_Chain, the chain of head term elements,
                            Free_HT_L, the list of free head term elements,
                            Cur_HT, count of number of frames. */
void Allocate_HT()
{
  HTSpt hpt;
  counter index;

  hpt = (HTSpt) calloc(HT_Frame_Size, sizeof(HTStorage));
  hpt[0].u.hspt = HT_Chain;
  hpt[0].lead = (element) calloc((HT_Frame_Size-1)*2*NumGen, sizeof(exponent));
  HT_Chain = &(hpt[0]);
  Cur_HT++;
  hpt[1].u.hspt = Free_HT_L;
  hpt[1].lead = hpt[0].lead+NumGen;
  hpt[1].pre = hpt[1].lead+NumGen;
  Free_HT_L = &(hpt[HT_Frame_Size-1]);
  for (index = 2; index < HT_Frame_Size; index++) {
    hpt[index].u.hspt = &(hpt[index-1]);
    hpt[index].lead = hpt[index-1].pre+NumGen;
    hpt[index].pre = hpt[index].lead+NumGen;
  }
}

/* This procedure frees a head term element.
   Require global variable: Free_HT_L, list of free head terms. */
void Free_HT(hpt)
HTSpt hpt;
{
  hpt->u.hspt = Free_HT_L;
  Free_HT_L = hpt;
  used_HT--;
}

/* This procedure frees a head term element.
   Require global variable: Free_HT_L, list of free head terms. */
void Free_Basis_HT(hpt)
HTSpt hpt;
{
  hpt->u.hspt = Free_HT_L;
  Free_HT_L = hpt;
}

/* This procedure makes sure that there is enough memory for the next
   Groebner basis. */
void Reallocate_Basis(n)
counter n;
{
  if (Cur_Basis < n) {
    free(Basis_HT);
    Basis_HT = (HTSpt *) calloc(n, sizeof(HTStorage));
    Cur_Basis = n;
  }
}

