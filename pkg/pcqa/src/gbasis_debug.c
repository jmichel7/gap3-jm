/*
  Written by: Eddie Lo
  Date started: August 26,1994.

  Part of the polycyclic quotient algorithm program.

  This program debugs for the reduce module.
*/

#include "pcqa.h"
#include <malloc.h>

extern counter Cur_Cell, Cur_MP_INT;
extern counter used_mod, used_HT, Num_HT;
extern counter Cur_Stack, Cur_HT, Cur_Mod;
extern counter Num_Exp_Pair, Num_Imm_Pair, Num_Late_Pair, Num_HT_Frame;
extern stack Exp_Stack, Imm_Stack, Late_Stack;


/* Print out the entries in the Immediate Stack.
   Require global variables: Imm_Stack, the immediate stack.
                             Num_Imm_Pair, number of immediate pairs. */
void Print_Imm_Stack()
{
  counter n, m1, m2;

  for (n = 0; n < Num_Imm_Pair; n++) {
    m1 = n >> Log_SFS;
    m2 = n^(m1 << Log_SFS);
    printf("%d, %d\n", Imm_Stack[m1][m2].smaller, Imm_Stack[m1][m2].bigger);
  }
}

/* Print out the entries in the Immediate Stack.
   Require global variables: Late_Stack, the immediate stack.
                             Num_Late_Pair, number of immediate pairs. */
void Print_Late_Stack()
{
  counter n, m1, m2;

  for (n = 0; n < Num_Late_Pair; n++) {
    m1 = n >> Log_SFS;
    m2 = n^(m1 << Log_SFS);
    printf("%d, %d\n", Late_Stack[m1][m2].smaller, Late_Stack[m1][m2].bigger);
  }
}

/* This procedure computes the amount of used in storing the polynomials
   in the list of polynomial.
   Require global variable: Num_HT, Number of polynomials in the list.
                            Poly_List, the list of polynomials. */
counter List_Memory()
{
  counter numterm, numcell, numint;
  counter tterm, tcell, tint, tmod, tHT;
  counter n, index;
  Poly f;
  ModSpt mpt, prev;
  HTSpt hpt;

  tterm = tcell = tint = tmod = tHT = 0;
  prev = NULL;
  for (n = 1; n <= Num_HT; n++) {
    if ((hpt = Get_HT(n)) != NULL) {
      tHT++;
      if ((mpt = hpt->pm) != prev) {
        tmod++;
        Memory_Mod(mpt->list, &numterm, &numcell, &numint, mpt->start);
        tterm += numterm;
        tcell += numcell;
        tint += numint;
        prev = mpt;
      }
    }
  };
  printf("%d terms, %d cells, %d integers,  %d head terms, %d mods.\n", tterm,
    tcell, tint, tHT, tmod);
  return (tterm*sizeof(MP_INT)+tcell*sizeof(cell)+tint*sizeof(long int)+
    tHT*sizeof(HTStorage)+
    tmod*sizeof(ModStorage));
}

/* Output some memory information. */
void Print_Memory()
{
  printf(
    "Memory in list = %d, %d Cell, %d Stack, %d HT, %d Mod, %d INTs, %d HTF\n",
    List_Memory(), Cur_Cell, Cur_Stack, Cur_HT, Cur_Mod, Cur_MP_INT,
    Num_HT_Frame);
/*  printf("Arena = %d\n", mallinfo().arena); */
}

/* Output some stack information. */
void Print_Stack()
{
  printf("%d Exp_Pair, %d Imm_Pair, %d Late_Pair, %d Num_HT, %d used_HT\n",
    Num_Exp_Pair, Num_Imm_Pair, Num_Late_Pair, Num_HT, used_HT);
}

/* Print information about a head term. */
void Print_HT(hpt)
HTSpt hpt;
{
  ModSpt mpt;

  printf("Leading: "); Save_Elm(stdout, hpt->lead);
  printf("High   : "); Save_Elm(stdout, hpt->pre);
  printf("Coeff  : "); mpz_out_str(stdout, 10, hpt->u.cff);
  mpt = hpt->pm;
  printf("\nnumterm: %d\nremain: %d\nstart: %d\nbig: ",
    mpt->numterm, mpt->remain, mpt->start);
  mpz_out_str(stdout, 10, mpt->u.weight);
  printf("\n");
}

