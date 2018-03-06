/*
  Written by: Eddie Lo
  Date started: August 26,1994.

  Part of the polycyclic quotient algorithm program.

  This program takes as inputs the "definition" of the module generators
  using positive generators and extends to rules involving both positive
  and negative generators.

Introduction:

  This module contains procedure like Extend_Rule in Sims' book.

Data description:

  ppmod, pnmod, npmod, nnmod, powermod, npowermod are the module
  elements which are appended to the polycyclic presentation of the
  previous quotient. ppmod corresponds to those will left sides
  a_i a_j, npmod to those with left sides a_i^(-1) a_j, etc.
  powermod and npowermod correspond to power relations and inverse
  power relations.

  ppprel, powerprel, ppowerrel, powerrel, drel, relrel are the
  module elements which are rules in determining the Groebner basis.
  ppprel is the relation from different ways of collecting
  a_i a_j a_k, powerprel a_i^m a_j, etc.
  drel relations are formed from definition of certain group
  generators. relrel are formed from original relations of the group.

  *start are the starting module elements of the *mod.

  table2 and table3 are indexing functions.

  module_definition stores the definitions of the module generators.

  gg_defined are flags which tell whether ppmod, pnmod, npmod and
  nnmod are defined.
  n_defined are flags which tell whether npowermod is defined.

  *_found are flags which tell whether the module relations are found.
*/

#include "pcqa.h"

extern Poly **ppmod, **pnmod, **npmod, **nnmod, **powermod, **npowermod;
extern Poly **defmod;
extern Poly **ppprel, **powerprel, **ppowerrel, **powerrel, **drel, **relrel;
extern generator *ppstart, *pnstart, *npstart, *nnstart;
extern generator *powerstart, *npowerstart;
extern counter *table2, *table3;
extern Poly *extend_m1, *extend_m2, *make_m1, *make_m2;
extern Rule_Set RUse;
extern def *module_definition;
extern flag *gg_defined, *n_defined, *d_found, *rel_found;
extern flag *ppp_found, *powerp_found, *ppower_found, *power_found;
extern autlist *negaut;

extern generator NumGen, NumGenP, NumMod;
extern exponent *Power_Elm, *Power_Inv_Elm, *finite_index;
extern flag DEBUG;
extern flag ExistsRule;
extern counter Size_R;
extern relterm **RelGen, **Poly_Pre;


/* This subprocedure makes rules of the form a1 a2 a3.
   Require global variables: make_m, a module element. */
flag ppp_Rule(i1, i2, i3, k)
generator i1, i2, i3;
counter k;
{
  generator start1, start2, j;
  counter entry;

  entry = table2[i1-1]+i2-i1-1;
  start1 = ppstart[entry];
  if (start1 < NumMod) Copy_Mod(make_m1, ppmod[entry], 0);
  else for (j = 0; j < NumMod; j++) make_m1[j] = NULL;
  Store_Elm(ppaut(i1, i2));
  Store_Gen(i1, 1);
  Store_Gen(i3, 1);
  if (Collect(make_m2, &start2) == NO) return NO;
  start1 = Merge_Mod(make_m1, make_m1, make_m2, start1, start2);
  Store_Gen(i1, 1);
  Store_Gen(i2, 1);
  Store_Gen(i3, 1);
  if (Collect(make_m2, &start2) == NO) return NO;
  Negate_Mod(make_m2, start2);
  start1 = Merge_Mod(make_m1, make_m1, make_m2, start1, start2);
  Basis_Reduce_Mod(make_m1, &start1);
  if (start1 < NumMod) {
    ppprel[k] = Store_Rule(&RUse, make_m1);
    return YES;
  };
  return USED;
}

/* This subprocedure makes rules of the form a1^n a2.
   Require global variables: make_m1, make_m2, module elements. */
flag powerp_Rule(i1, i2, k)
generator i1, i2;
counter k;
{
  generator start1, start2, j;

  start1 = powerstart[i1];
  if (start1 < NumMod) Copy_Mod(make_m1, powermod[i1], 0);
  else for (j = 0; j < NumMod; j++) make_m1[j] = NULL;
  Store_Elm(powelm(i1));
  Store_Gen(i2, 1);
  if (Collect(make_m2, &start2) == NO) return NO;
  start1 = Merge_Mod(make_m1, make_m1, make_m2, start1, start2);
  Store_Gen(i1, finite_index[i1]);
  Store_Gen(i2, 1);
  if (Collect(make_m2, &start2) == NO) return NO;
  Negate_Mod(make_m2, start2);
  start1 = Merge_Mod(make_m1, make_m1, make_m2, start1, start2);
  Basis_Reduce_Mod(make_m1, &start1);
  if (start1 < NumMod) {
    powerprel[k] = Store_Rule(&RUse, make_m1);
    return YES;
  };
  return USED;
}

/* This subprocedure makes rules of the form a1 a2^n.
   Require global variables: make_m1, make_m2, module elements. */
flag ppower_Rule(i1, i2, k)
generator i1, i2;
counter k;
{
  counter entry;
  generator start1, start2, j;

  entry = table2[i1-1]+i2-i1-1;
  start1 = ppstart[entry];
  if (start1 < NumMod) Copy_Mod(make_m1, ppmod[entry], 0);
  else for (j = 0; j < NumMod; j++) make_m1[j] = NULL;
  Store_Elm(ppaut(i1, i2));
  Store_Gen(i1, 1);
  if(finite_index[i2] > 1) Store_Gen(i2, finite_index[i2]-1);
  if (Collect(make_m2, &start2) == NO) return NO;
  start1 = Merge_Mod(make_m1, make_m1, make_m2, start1, start2);
  Store_Gen(i1, 1);
  Store_Gen(i2, finite_index[i2]);
  if (Collect(make_m2, &start2) == NO) return NO;
  Negate_Mod(make_m2, start2);
  start1 = Merge_Mod(make_m1, make_m1, make_m2, start1, start2);
  Basis_Reduce_Mod(make_m1, &start1);
  if (start1 < NumMod) {
    ppowerrel[k] = Store_Rule(&RUse, make_m1);
    return YES;
  };
  return USED;
}

/* This subprocedure makes rules of the form a1^(n+1).
   Require global variables: make_m1, make_m2, module elements. */
flag powerpower_Rule(i1)
generator i1;
{
  generator start1, start2, j;

  start1 = powerstart[i1];
  if (start1 < NumMod) Copy_Mod(make_m1, powermod[i1], 0);
  else for (j = 0; j < NumMod; j++) make_m1[j] = NULL;
  if (finite_index[i1] == 1) {
    Store_Elm(powelm(i1));
    Store_Gen(i1, 1);
    if (Collect(make_m2, &start2) == NO) return NO;
    start1 = Merge_Mod(make_m1, make_m1, make_m2, start1, start2);
  };
  Store_Gen(i1, 1);
  Store_Gen(i1, finite_index[i1]);
  if (Collect(make_m2, &start2) == NO) return NO;
  Negate_Mod(make_m2, start2);
  start1 = Merge_Mod(make_m1, make_m1, make_m2, start1, start2);
  Basis_Reduce_Mod(make_m1, &start1);
  if (start1 < NumMod) {
    powerrel[i1] = Store_Rule(&RUse, make_m1);
    return YES;
  };
  return USED;
}

/* This procedure generates rules from forcing coincidence. */
void Make_Rule()
{
  generator i1, i2, i3;
  counter index1, index2;

  index1 = index2 = 0;
  for (i1 = 1; i1 <= NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      for (i3 = i2+1; i3 <= NumGen; i3++) {
        if (ppp_found[index1] == NO)
          ppp_found[index1] = ppp_Rule(i1, i2, i3, index1);
        index1++;
      };
      if (finite_index[i1] != 0 && powerp_found[index2] == NO)
        powerp_found[index2] = powerp_Rule(i1, i2, index2);
      if (finite_index[i2] != 0 && ppower_found[index2] == NO)
        ppower_found[index2] = ppower_Rule(i1, i2, index2);
      index2++;
    };
    if (finite_index[i1] > 1 && power_found[i1] == NO)
      power_found[i1] = powerpower_Rule(i1);
  }
}

/* This subprocedure generates rule by inspecting definition.
   Require global variable: extend_m1, a module element. */
flag d_Rule(i)
generator i;
{
  generator start;

  Store_Rel(Poly_Pre[i-1]);
  if (Collect(extend_m1, &start) == NO) return NO;
  Basis_Reduce_Mod(extend_m1, &start);
  if (start < NumMod) {
    drel[i] = Store_Rule(&RUse, extend_m1);
    return YES;
  };
  return USED;
}

/* This procedure generates rule by inspecting definition. */
void Definition_Rule()
{
  generator i;

  for (i = 1; i <= NumGen; i++) if (d_found[i] == NO) d_found[i] = d_Rule(i);
}

/* This procedure generates rule by inspecting relation.
   Require global variable: extend_m1, a module element. */
flag r_Rule(i)
generator i;
{
  generator start;

  Store_Rel(RelGen[i]);
  if (Collect(extend_m1, &start) == NO) return NO;
  Basis_Reduce_Mod(extend_m1, &start);
  if (start < NumMod) {
    relrel[i] = Store_Rule(&RUse, extend_m1);
    return YES;
  };
  return USED;
}

/* This procedure generates rule by inspecting relation. */
void Relation_Rule()
{
  generator i;

  for (i = 0; i < Size_R; i++) if (!rel_found[i]) rel_found[i] = r_Rule(i);
}

