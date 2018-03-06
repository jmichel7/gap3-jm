/*
  Written by: Eddie Lo
  Date started: October 24, 1995.

  Part of the polycyclic quotient algorithm package.

  This program manages memory for the ring module.
*/

#include <stdlib.h>
#include "pcqa.h"

extern generator NumGen, NumGenP;

MP_INT MP_structure[6];

cellptr *copy_from_h = NULL;
cellptr *copy_to_h = NULL;

cellptr *negate_h = NULL;

cellptr *merge1_h = NULL;
cellptr *merge2_h = NULL;
cellptr *merge_prev = NULL;

cellptr *print_h = NULL;
element print_u = NULL;

element read_u = NULL;
coeff read_num = NULL;

cellptr *mult_from_h = NULL;
cellptr *mult_to_h = NULL;
cellptr *mult_head = NULL;
element mult_u = NULL;
element mult_tmp = NULL;
element *mult_upt = NULL;
generator *mult_gen = NULL;

cellptr *free_h = NULL;

cellptr *sym_pos = NULL;
cellptr *sym_neg = NULL;
element sym_tar_u = NULL;
element sym_high_u = NULL;
generator *sym_gen = NULL;

cellptr *full_h1 = NULL;
cellptr *full_h2 = NULL;
element full_mul_u = NULL;
element full_ht_u = NULL;
element full_tar_u = NULL;
generator *full_gen = NULL;
coeff full_div = NULL;
coeff full_rem = NULL;
coeff full_n1 = NULL;
coeff full_n2 = NULL;

cellptr *modulo_h = NULL;
cellptr *modulo_prev = NULL;
coeff modulo_div = NULL;

cellptr *memory_h = NULL;

element mult_pp_u = NULL;

element aug_u = NULL;

element split_u = NULL;


/* -----  Initialization  ----- */

/* This procedure initializes for the ring module. */
void Est_Structure()
{
  counter index;

  for (index = 0; index < 6; index++) mpz_init(MP_structure+index);
  read_num = MP_structure;
  full_div = &MP_structure[1];
  full_rem = &MP_structure[2];
  full_n1 = &MP_structure[3];
  full_n2 = &MP_structure[4];
  modulo_div = &MP_structure[5];
}

/* This procedure initializes for the ring module, used many times. */
void Init_Structure()
{
  generator i;
  element elm;

  copy_from_h = (cellptr *) calloc(NumGenP, sizeof(cellptr));
  copy_to_h = (cellptr *) calloc(NumGenP, sizeof(cellptr));

  negate_h = (cellptr *) calloc(NumGenP, sizeof(cellptr));

  merge1_h = (cellptr *) calloc(NumGenP, sizeof(cellptr));
  merge2_h = (cellptr *) calloc(NumGenP, sizeof(cellptr));
  merge_prev = (cellptr *) calloc(NumGenP, sizeof(cellptr));
  merge1_h[0] = (cellptr) malloc(sizeof(cell));

  print_h = (cellptr *) calloc(NumGenP, sizeof(cellptr));
  print_u = (element) calloc(NumGen, sizeof(exponent))+NumGen;
  for (i = 1; i <= NumGen; i++) print_u[-i] = 0;

  read_u = (element) calloc(NumGen, sizeof(exponent))+NumGen;

  mult_from_h = (cellptr *) calloc(NumGenP, sizeof(cellptr));
  mult_to_h = (cellptr *) calloc(NumGenP, sizeof(cellptr));
  mult_head = (cellptr *) calloc(NumGenP, sizeof(cellptr));
  mult_u = (element) calloc(2*NumGen, sizeof(exponent))+NumGen;
  mult_tmp = mult_u+NumGen;
  elm = (element) calloc(NumGen*NumGen, sizeof(exponent));
  mult_upt = (element *) calloc(NumGenP, sizeof(element));
  for (i = 1; i <= NumGen; i++) mult_upt[i] = elm+i*NumGen;
  mult_gen = (generator *) calloc(NumGenP, sizeof(generator));
  mult_to_h[0] = (cellptr) malloc(sizeof(cell));
  mult_to_h[0]->ind = 0;

  free_h = (cellptr *) calloc(NumGenP, sizeof(cellptr));

  sym_pos = (cellptr *) calloc(NumGenP, sizeof(cellptr));
  sym_neg = (cellptr *) calloc(NumGenP, sizeof(cellptr));
  sym_tar_u = (element) calloc(NumGenP*2, sizeof(exponent))+NumGenP;
  sym_high_u = sym_tar_u+NumGenP;
  sym_gen = (generator *) calloc(NumGenP, sizeof(generator));
  for (i = 1; i <= NumGenP; i++) sym_tar_u[-i] = sym_high_u[-i] = 0;

  full_h1 = (cellptr *) calloc(NumGenP, sizeof(cellptr));
  full_h2 = (cellptr *) calloc(NumGenP, sizeof(cellptr));
  full_mul_u = (element) calloc(NumGenP*3, sizeof(exponent))+NumGenP;
  full_ht_u = full_mul_u+NumGenP;
  full_tar_u = full_ht_u+NumGenP;
  full_gen = (generator *) calloc(NumGenP, sizeof(generator));

  modulo_h = (cellptr *) calloc(NumGenP, sizeof(cellptr));
  modulo_prev = (cellptr *) calloc(NumGenP, sizeof(cellptr));
  modulo_h[0] = (cellptr) malloc(sizeof(cell));

  memory_h = (cellptr *) calloc(NumGenP, sizeof(cellptr));

  mult_pp_u = (element) calloc(NumGen, sizeof(exponent))+NumGen;

  aug_u = (element) calloc(NumGen, sizeof(exponent))+NumGen;
  for (i = NumGen; i >= 1; i--) aug_u[-i] = 0;

  split_u = (element) calloc(NumGen, sizeof(exponent))+NumGen;
}


/* -----  Resetting  ----- */

/* This procedure frees up memory for the ring module. */
void Reset_Structure()
{
  free(copy_from_h); free(copy_to_h);
  free(negate_h);
  free(merge1_h[0]); free(merge1_h); free(merge2_h); free(merge_prev); 
  free(print_h); free(print_u-NumGen);
  free(read_u-NumGen);
  free(mult_from_h); 
  free(mult_to_h[0]); free(mult_to_h); free(mult_head);
  free(mult_u-NumGen); free(mult_upt[1]-NumGen); free(mult_upt);
  free(mult_gen); 
  free(free_h);
  free(sym_pos); free(sym_neg); free(sym_tar_u-NumGenP); free(sym_gen);
  free(full_h1); free(full_h2); free(full_mul_u-NumGenP); free(full_gen);
  free(modulo_h[0]);free(modulo_h); free(modulo_prev); 
  free(memory_h);
  free(mult_pp_u-NumGen);
  free(aug_u-NumGen);
  free(split_u-NumGen);
}

