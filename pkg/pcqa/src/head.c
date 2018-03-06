/*
  Written by: Eddie Lo
  Date started: December 5, 1995.

  Part of the polycyclic quotient algorithm package.

  This is the "head" routine, as opposed to the "tail" routine
  in the ANU nilpotent quotient package. Module elements here are
  collected to the left.

  Data description can be found in extend.c.
*/

#include "pcqa.h"

extern Poly **ppmod, **pnmod, **npmod, **nnmod, **powermod, **npowermod;
extern Poly **defmod;
extern Poly *head_m1, *head_m2;
extern def *module_definition;
extern generator *ppstart, *pnstart, *npstart, *nnstart;
extern generator *powerstart, *npowerstart, *defstart;
extern generator NumGen, NumMod;
extern exponent *Power_Elm, *Power_Inv_Elm, *finite_index;
extern element Identity, head_u;
extern autlist *posaut, *negaut;
extern flag *gg_defined, *n_defined;
extern counter *table2;
extern MP_INT one, neg_one;
extern counter Size_S, choose2;
extern flag ExistsModDef;
extern flag DEBUG;


/* This procedure computes the number of head module elements required. */
counter Num_HeadMod()
{
  return (4*choose2+2*NumGen+Size_S);
}

/* This subprocedure of Make_Head computes the head of a pn pair.
   Require global variable: head_m1, head_m2, a module element. */
flag pn_Head(i, j, entry)
generator i, j;
counter entry;
{
  generator start1, start2, index;

  if ((start1 = ppstart[entry]) < NumMod)
    Mult_Mod(head_m1, ppmod[entry], ppaut(i, j), &neg_one, start1);
  else for (index = 0; index < NumMod; index++) head_m1[index] = NULL;
  Store_Elm_Inv(ppaut(i, j));
  if (!Collect(head_m2, &start2)) {
    printf("Abort collecting a%d A%d.\n", i, j);
    return 0;
  };
  gg_defined[entry] |= 2;
  start1 = Merge_Mod(head_m1, head_m1, head_m2, start1, start2);
  Basis_Reduce_Mod(head_m1, &start1);
  if (start1 < NumMod) {
    pnmod[entry] = New_HeadMod();
    for (index = 0; index < NumMod; index++)
      pnmod[entry][index] = head_m1[index];
  }
  else pnmod[entry] = NULL;
  pnstart[entry] = start1;
  return 1;
}

/* This subprocedure of Make_Head computes the head of a np pair.
   Require global variable: head_m1, a module element.
                            head_u, an element. */
flag np_Head(i, j, entry)
generator i, j;
counter entry;
{
  generator start1, start2, index;

  Store_Gen(j, -1);
  Store_Gen(i, 1);
  Store_Elm(npaut(i, j));
  if (!Collect(head_m1, &start1)) {
    printf("Abort collecting A%d a%d.\n", i, j);
    return 0;
  };
  gg_defined[entry] |= 4;
  if (start1 < NumMod) {
    Identity[-i] = 1; Identity[-j] = -1;
    Group_Pow_Check(head_u, Identity);
    Identity[-i] = Identity[-j] = 0;
    Mult_Mod(head_m2, head_m1, head_u, &neg_one, start1);
    start2 = start1;
    Basis_Reduce_Mod(head_m2, &start2);
    if (start2 < NumMod) {
      npmod[entry] = New_HeadMod();
      for (index = 0; index < NumMod; index++)
        npmod[entry][index] = head_m2[index];
    }
    else npmod[entry] = NULL;
    npstart[entry] = start2;
  }
  else {
    npmod[entry] = NULL;
    npstart[entry] = NumMod;
  };
  Free_Mod_Cell(head_m1, start1);
  return 1;
}

/* This subprocedure of Make_Head computes the head of a nn pair.
   Require global variable: head_m1, head_m2, module elements. */
flag nn_Head(i, j, entry)
generator i, j;
counter entry;
{
  Poly *fp;
  generator start1, start2, index;

  if (!(gg_defined[entry] & 4)) {
    printf("Abort collecting A%d A%d.\n", i, j);
    return 0;
  };
  if ((start1 = npstart[entry]) < NumMod)
    Mult_Mod(head_m1, npmod[entry], npaut(i, j), &neg_one, start1);
  else for (index = 0; index < NumMod; index++) head_m1[index] = NULL;
  Store_Elm_Inv(npaut(i, j));
  if (!Collect(head_m2, &start2)) {
    printf("Abort collecting A%d A%d.\n", i, j);
    return 0;
  };
  gg_defined[entry] |= 8;
  start1 = Merge_Mod(head_m1, head_m1, head_m2, start1, start2);
  Basis_Reduce_Mod(head_m1, &start1);
  if (start1 < NumMod) {
    nnmod[entry] = New_HeadMod();
    for (index = 0; index < NumMod; index++)
      nnmod[entry][index] = head_m1[index];
  }
  else nnmod[entry] = NULL;
  nnstart[entry] = start1;
  return 1;
}

/* This subprocedure of Make_Head computes the head of a npower pair.
   Require global variable: head_m1, head_m2, module elements. */
flag npower_Head(i)
generator i;
{
  Poly *fp;
  generator start1, start2, index;

  if ((start1 = powerstart[i]) < NumMod)
    Mult_Mod(head_m1, powermod[i], powelm(i), &neg_one, start1);
  else for (index = 0; index < NumMod; index++) head_m1[index] = NULL;
  Store_Elm_Inv(powelm(i));
  if (!Collect(head_m2, &start2)) {
    printf("Abort collecting A%d.\n", i);
    return 0;
  };
  n_defined[i] = 1;
  start1 = Merge_Mod(head_m1, head_m1, head_m2, start1, start2);
  Basis_Reduce_Mod(head_m1, &start1);
  if (start1 < NumMod) {
    npowermod[i] = New_HeadMod();
    for (index = 0; index < NumMod; index++)
      npowermod[i][index] = head_m1[index];
  }
  else npowermod[i] = NULL;
  npowerstart[i] = start1;
  return 1;
}

/* Given the collection rules for positive generators and power relations,
   find the other rules. */
void Make_Head()
{
  generator i, j;
  counter entry;

  for (i = NumGen; i >= 1; i--) {
    for (j = NumGen; j > i; j--) {
      entry = table2[i-1]+j-i-1;
      if (!finite_index[j] && !(gg_defined[entry] & 2)) pn_Head(i, j, entry);
    };
    for (j = NumGen; j > i; j--) {
      entry = table2[i-1]+j-i-1;
      if (!finite_index[i] && !(gg_defined[entry] & 4)) np_Head(i, j, entry);
      if (!finite_index[i] && !finite_index[j] && !(gg_defined[entry] & 8))
        nn_Head(i, j, entry);
    };
    if (finite_index[i] != 0 && n_defined[i] == 0) npower_Head(i);
  }
}

/* Store the definition of module generators. This procedure mirrors the
   procedure Define_Head. */
void Define_ModDef()
{
  generator i, j1, j2;
  counter entry;

  entry = 0;
  for (i = 1; i <= Size_S; i++) {
    module_definition[entry].ind = DEF;
    module_definition[entry].lower = i;
    entry++;
  };
  for (i = 1; i <= NumGen; i++) {
    if (finite_index[i]) {
      module_definition[entry].ind = EXP;
      module_definition[entry].lower = i;
      entry++;
    }
  };
  for (i = 1; i < NumGen; i++) {
    for (j1 = i+1; j1 <= NumGen; j1++) {
      module_definition[entry].ind = CONJ;
      module_definition[entry].lower = i;
      module_definition[entry].upper = j1;
      entry++;
    }
  };
  ExistsModDef = USED;
}

/* Form an initial definition of module generators. */
void Define_Head()
{
  generator i, j1, j2;
  counter k, entry, left;
  element u;

  entry = 0;
  for (i = 0; i < Size_S; i++) {
    defstart[i] = entry;
    defmod[i] = New_HeadMod();
    for (j1 = 0; j1 < NumMod; j1++) defmod[i][j1] = NULL;
    defmod[i][entry] = Constant_Poly(&one, 0);
    entry++;
  };
  for (i = 1; i <= NumGen; i++) {
    if (finite_index[i] != 0) {
      powerstart[i] = entry;
      powermod[i] = New_HeadMod();
      for (j1 = 0; j1 < NumMod; j1++) powermod[i][j1] = NULL;
      powermod[i][entry] = Constant_Poly(&one, 0);
      entry++;
    }
    else powerstart[i] = NumMod;
  };
  k = 0;
  for (i = 1; i < NumGen; i++) {
    for (j1 = i+1; j1 <= NumGen; j1++) {
      ppstart[k] = entry;
      ppmod[k] = New_HeadMod();
      for (j2 = 0; j2 < NumMod; j2++) ppmod[k][j2] = NULL;
      ppmod[k][entry] = Constant_Poly(&one, 0);
      k++; entry++;
    }
  }
}

