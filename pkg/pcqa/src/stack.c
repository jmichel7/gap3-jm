/*
  Written by: Eddie Lo
  Date started: October 24, 1995

  Part of the polycyclic quotient algorithm package.

  This program operates on the stack for the reduce module.
*/

#include <stdlib.h>
#include "pcqa.h"

extern counter Num_Exp_Pair, Num_Imm_Pair, Num_Late_Pair;
extern counter Cur_Stack;
extern stackframe Free_Frame[2];
extern stack Exp_Stack, Imm_Stack, Late_Stack;
extern counter pair_length1, pair_length2, pair_length3;


/* -----  Late stack  ----- */

/* This procedure pushes a pair onto Late_Stack.
   Require global variable: Late_Stack, immediate stack.
                            Num_Late_Pair, number of immediate pairs.
                            Free_Frame, storage for free frame. */
void Push_Late_Stack(n1, n2)
counter n1, n2;
{
  counter m1, m2;

  m1 = Num_Late_Pair >> Log_SFS;
  m2 = Num_Late_Pair^(m1 << Log_SFS);
  if (!m2) {
    if (Free_Frame[1] != NULL) {
      Late_Stack[m1] = Free_Frame[1];
      Free_Frame[1] = NULL;
    }
    else if (Free_Frame[0] != NULL) {
      Late_Stack[m1] = Free_Frame[0];
      Free_Frame[0] = NULL;
    }
    else {
      Late_Stack[m1] = (stackframe) calloc(Stack_Sublist_Size, sizeof(pair));
      Cur_Stack++;
    }
  };
  Late_Stack[m1][m2].smaller = n1;
  Late_Stack[m1][m2].bigger = n2;
  Num_Late_Pair++;
}

/* This procedure pops a pair out of Late_Stack. Pairs are chosen so that
   the sum of the number of terms of the polynomials is smaller than
   n if there exists such pair.
   Require global variable: Late_Stack, immediate stack.
                            Num_Late_Pair, number of immediate pairs.
                            Free_Frame, storage for free frames. */
flag Search_Late_Stack(n1, n2, hpt1, hpt2)
counter *n1, *n2;
HTSpt *hpt1, *hpt2;
{
  counter m1, m2;
  counter index, size, new_size;
  counter in1, in2;
  flag nopair, shorter;

  nopair = 1;
  while (Num_Late_Pair > 0 && nopair) {
    Num_Late_Pair--;
    m1 = Num_Late_Pair >> Log_SFS;
    m2 = Num_Late_Pair^(m1 << Log_SFS);
    in1 = Late_Stack[m1][m2].smaller;
    in2 = Late_Stack[m1][m2].bigger;
    if ((*hpt1 = Get_HT(in1)) != NULL && (*hpt2 = Get_HT(in2)) != NULL)
      nopair = 0;
    if (!m2) {
      if (Free_Frame[1] == NULL) Free_Frame[1] = Late_Stack[m1];
      else if (Free_Frame[0] == NULL) Free_Frame[0] = Late_Stack[m1];
      else {
        free(Late_Stack[m1]);
        Cur_Stack--;
      }
    }
  };
  if (nopair) {
    *n1 = 0;
    shorter = 0;
  }
  else {
    if ((size = ((*hpt1)->pm)->numterm+((*hpt2)->pm)->numterm) > pair_length3) {
      nopair = 1;
      for (index = 0; index < Num_Late_Pair && nopair; index++) {
        m1 = index >> Log_SFS;
        m2 = index^(m1 << Log_SFS);
        *n1 = Late_Stack[m1][m2].smaller;
        *n2 = Late_Stack[m1][m2].bigger;
        if ((*hpt1 = Get_HT(*n1)) != NULL && (*hpt2 = Get_HT(*n2)) != NULL &&
          (new_size = ((*hpt1)->pm)->numterm+((*hpt2)->pm)->numterm) < size) {
          Late_Stack[m1][m2].smaller = in1;
          Late_Stack[m1][m2].bigger = in2;
          in1 = *n1; in2 = *n2;
          if ((size = new_size) <= pair_length3) nopair = 0;
        }
      }
    };
    if (size > pair_length3) shorter = 0;
    else shorter = 1;
    *n1 = in1; *n2 = in2;
    *hpt1 = Get_HT(*n1); *hpt2 = Get_HT(*n2);
  };
  return shorter;
}

/* This procedure pops a pair out of Late_Stack. Only non-null pairs
   are returned. *n1 = 0 if no non-null pair exists.
   Require global variable: Late_Stack, immediate stack.
                            Num_Late_Pair, number of late pairs.
                            Free_Frame, storage for free frames. */
void Pop_Late_Stack(n1, n2, hpt1, hpt2)
counter *n1, *n2;
HTSpt *hpt1, *hpt2;
{
  counter m1, m2;
  flag nopair;

  nopair = 1;
  while (Num_Late_Pair > 0 && nopair) {
    Num_Late_Pair--;
    m1 = Num_Late_Pair >> Log_SFS;
    m2 = Num_Late_Pair^(m1 << Log_SFS);
    *n1 = Late_Stack[m1][m2].smaller;
    *hpt1 = Get_HT(*n1);
    if (*hpt1 != NULL) {
      *n2 = Late_Stack[m1][m2].bigger;
      *hpt2 = Get_HT(*n2);
      if (*hpt2 != NULL) nopair = 0;
    };
    if (!m2) {
      if (Free_Frame[1] == NULL) Free_Frame[1] = Late_Stack[m1];
      else if (Free_Frame[0] == NULL) Free_Frame[0] = Late_Stack[m1];
      else {
        free(Late_Stack[m1]);
        Cur_Stack--;
      }
    }
  };
  if (nopair) *n1 = 0;
}


/* -----  Immediate stack  ----- */

/* This procedure pushes a pair onto Imm_Stack.
   Require global variables: Imm_Stack, immediate stack.
                             Num_Imm_Pair, number of immediate pairs.
                             Free_Frame, storage for free frame. */
void Push_Imm_Stack(n1, n2)
counter n1, n2;
{
  counter m1, m2;

  m1 = Num_Imm_Pair >> Log_SFS;
  m2 = Num_Imm_Pair^(m1 << Log_SFS);
  if (!m2) {
    if (Free_Frame[1] != NULL) {
      Imm_Stack[m1] = Free_Frame[1];
      Free_Frame[1] = NULL;
    }
    else if (Free_Frame[0] != NULL) {
      Imm_Stack[m1] = Free_Frame[0];
      Free_Frame[0] = NULL;
    }
    else {
      Imm_Stack[m1] = (stackframe) calloc(Stack_Sublist_Size, sizeof(pair));
      Cur_Stack++;
    }
  };
  Imm_Stack[m1][m2].smaller = n1;
  Imm_Stack[m1][m2].bigger = n2;
  Num_Imm_Pair++;
}

/* This procedure pops a pair out of Imm_Stack. Pairs are chosen so that
   the sum of the number of terms of the polynomials is smaller than
   n if there exists such pair.
   Require global variable: Imm_Stack, immediate stack.
                            Num_Imm_Pair, number of immediate pairs.
                            Free_Frame, storage for free frames. */
flag Search_Imm_Stack(n1, n2, hpt1, hpt2)
counter *n1, *n2;
HTSpt *hpt1, *hpt2;
{
  counter m1, m2;
  counter index, size, new_size;
  counter in1, in2;
  flag nopair, shorter;

  nopair = 1;
  while (Num_Imm_Pair > 0 && nopair) {
    Num_Imm_Pair--;
    m1 = Num_Imm_Pair >> Log_SFS;
    m2 = Num_Imm_Pair^(m1 << Log_SFS);
    in1 = Imm_Stack[m1][m2].smaller;
    in2 = Imm_Stack[m1][m2].bigger;
    if ((*hpt1 = Get_HT(in1)) != NULL && (*hpt2 = Get_HT(in2)) != NULL)
      nopair = 0;
    if (!m2) {
      if (Free_Frame[1] == NULL) Free_Frame[1] = Imm_Stack[m1];
      else if (Free_Frame[0] == NULL) Free_Frame[0] = Imm_Stack[m1];
      else {
        free(Imm_Stack[m1]);
        Cur_Stack--;
      }
    }
  };
  if (nopair) {
    *n1 = 0;
    shorter = 0;
  }
  else {
    if ((size = ((*hpt1)->pm)->numterm+((*hpt2)->pm)->numterm) > pair_length2) {
      nopair = 1;
      for (index = 0; index < Num_Imm_Pair && nopair; index++) {
        m1 = index >> Log_SFS;
        m2 = index^(m1 << Log_SFS);
        *n1 = Imm_Stack[m1][m2].smaller;
        *n2 = Imm_Stack[m1][m2].bigger;
        if ((*hpt1 = Get_HT(*n1)) != NULL && (*hpt2 = Get_HT(*n2)) != NULL &&
          (new_size = ((*hpt1)->pm)->numterm+((*hpt2)->pm)->numterm) < size) {
          Imm_Stack[m1][m2].smaller = in1;
          Imm_Stack[m1][m2].bigger = in2;
          in1 = *n1; in2 = *n2;
          if ((size = new_size) <= pair_length2) nopair = 0;
        }
      }
    };
    if (size > pair_length2) shorter = 0;
    else shorter = 1;
    *n1 = in1; *n2 = in2;
    *hpt1 = Get_HT(*n1); *hpt2 = Get_HT(*n2);
  };
  return shorter;
}

/* This procedure pops a pair out of Imm_Stack. Only non-null pairs
   are returned. *n1 = 0 if no non-null pair exists.
   Require global variable: Imm_Stack, immediate stack.
                            Num_Imm_Pair, number of immediate pairs.
                            Free_Frame, storage for free frames. */
void Pop_Imm_Stack(n1, n2, hpt1, hpt2)
counter *n1, *n2;
HTSpt *hpt1, *hpt2;
{
  counter m1, m2;
  flag nopair;

  nopair = 1;
  while (Num_Imm_Pair > 0 && nopair) {
    Num_Imm_Pair--;
    m1 = Num_Imm_Pair >> Log_SFS;
    m2 = Num_Imm_Pair^(m1 << Log_SFS);
    *n1 = Imm_Stack[m1][m2].smaller;
    *hpt1 = Get_HT(*n1);
    if (*hpt1 != NULL) {
      *n2 = Imm_Stack[m1][m2].bigger;
      *hpt2 = Get_HT(*n2);
      if (*hpt2 != NULL) nopair = 0;
    };
    if (!m2) {
      if (Free_Frame[1] == NULL) Free_Frame[1] = Imm_Stack[m1];
      else if (Free_Frame[0] == NULL) Free_Frame[0] = Imm_Stack[m1];
      else {
        free(Imm_Stack[m1]);
        Cur_Stack--;
      }
    }
  };
  if (nopair) *n1 = 0;
}


/* -----  Express stack ----- */

/* This procedure pushes a pair onto Imm_Stack.
   Require global variables: Exp_Stack, immediate stack.
                             Num_Exp_Pair, number of immediate pairs.
                             Free_Frame, storage for free frame. */
void Push_Exp_Stack(n1, n2)
counter n1, n2;
{
  counter m1, m2;

  m1 = Num_Exp_Pair >> Log_SFS;
  m2 = Num_Exp_Pair^(m1 << Log_SFS);
  if (!m2) {
    if (Free_Frame[1] != NULL) {
      Exp_Stack[m1] = Free_Frame[1];
      Free_Frame[1] = NULL;
    }
    else if (Free_Frame[0] != NULL) {
      Exp_Stack[m1] = Free_Frame[0];
      Free_Frame[0] = NULL;
    }
    else {
      Exp_Stack[m1] = (stackframe) calloc(Stack_Sublist_Size, sizeof(pair));
      Cur_Stack++;
    }
  };
  Exp_Stack[m1][m2].smaller = n1;
  Exp_Stack[m1][m2].bigger = n2;
  Num_Exp_Pair++;
}

/* This procedure pops a pair out of Exp_Stack. Pairs are chosen so that
   the sum of the number of terms of the polynomials is smaller than
   n if there exists such pair.
   Require global variable: Exp_Stack, immediate stack.
                            Num_Exp_Pair, number of immediate pairs.
                            Free_Frame, storage for free frames. */
flag Search_Exp_Stack(n1, n2, hpt1, hpt2)
counter *n1, *n2;
HTSpt *hpt1, *hpt2;
{
  counter m1, m2;
  counter index, size, new_size;
  counter in1, in2;
  flag nopair, shorter;

  nopair = 1;
  while (Num_Exp_Pair > 0 && nopair) {
    Num_Exp_Pair--;
    m1 = Num_Exp_Pair >> Log_SFS;
    m2 = Num_Exp_Pair^(m1 << Log_SFS);
    in1 = Exp_Stack[m1][m2].smaller;
    in2 = Exp_Stack[m1][m2].bigger;
    if ((*hpt1 = Get_HT(in1)) != NULL && (*hpt2 = Get_HT(in2)) != NULL)
      nopair = 0;
    if (!m2) {
      if (Free_Frame[1] == NULL) Free_Frame[1] = Exp_Stack[m1];
      else if (Free_Frame[0] == NULL) Free_Frame[0] = Exp_Stack[m1];
      else {
        free(Exp_Stack[m1]);
        Cur_Stack--;
      }
    }
  };
  if (nopair) {
    *n1 = 0;
    shorter = 0;
  }
  else {
    if ((size = ((*hpt1)->pm)->numterm+((*hpt2)->pm)->numterm) > pair_length1) {
      nopair = 1;
      for (index = 0; index < Num_Exp_Pair && nopair; index++) {
        m1 = index >> Log_SFS;
        m2 = index^(m1 << Log_SFS);
        *n1 = Exp_Stack[m1][m2].smaller;
        *n2 = Exp_Stack[m1][m2].bigger;
        if ((*hpt1 = Get_HT(*n1)) != NULL && (*hpt2 = Get_HT(*n2)) != NULL &&
          (new_size = ((*hpt1)->pm)->numterm+((*hpt2)->pm)->numterm) < size) {
          Exp_Stack[m1][m2].smaller = in1;
          Exp_Stack[m1][m2].bigger = in2;
          in1 = *n1; in2 = *n2;
          if ((size = new_size) <= pair_length1) nopair = 0;
        }
      }
    };
    if (size > pair_length1) shorter = 0;
    else shorter = 1;
    *n1 = in1; *n2 = in2;
    *hpt1 = Get_HT(*n1); *hpt2 = Get_HT(*n2);
  };
  return shorter;
}

/* This procedure pops a pair out of Imm_Stack. Only non-null pairs
   are returned. *n1 = 0 if no non-null pair exists.
   Require global variable: Exp_Stack, express stack.
                            Num_Exp_Pair, number of immediate pairs.
                            Free_Frame, storage for free frames. */
void Pop_Exp_Stack(n1, n2, hpt1, hpt2)
counter *n1, *n2;
HTSpt *hpt1, *hpt2;
{
  counter m1, m2;
  flag nopair;

  nopair = 1;
  while (Num_Exp_Pair > 0 && nopair) {
    Num_Exp_Pair--;
    m1 = Num_Exp_Pair >> Log_SFS;
    m2 = Num_Exp_Pair^(m1 << Log_SFS);
    *n1 = Exp_Stack[m1][m2].smaller;
    *hpt1 = Get_HT(*n1);
    if (*hpt1 != NULL) {
      *n2 = Exp_Stack[m1][m2].bigger;
      *hpt2 = Get_HT(*n2);
      if (*hpt2 != NULL) nopair = 0;
    };
    if (!m2) {
      if (Free_Frame[1] == NULL) Free_Frame[1] = Exp_Stack[m1];
      else if (Free_Frame[0] == NULL) Free_Frame[0] = Exp_Stack[m1];
      else {
        free(Exp_Stack[m1]);
        Cur_Stack--;
      }
    }
  };
  if (nopair) *n1 = 0;
}

