/*
  Written by : Eddie Lo
  Date started: February 21, 1996.

  Part of the polycyclic quotient algorithm package.

  This program processes special rules in the rule module.
*/

#include "pcqa.h"
#define avail(k) (!(RUse.Active[(k)] & (chosen | utilized | deleted)))

extern Rule_Set RUse;
extern generator NumGen, NumMod;


/* This procedure finds the "next" "permutation" given a permutation perm
   of length limit. For example, the next of 4 4 1 4 2 when limit = 5
   and NumGen = 4 is 1 1 2 4 2. */
flag Next_Permutation(perm, limit)
generator *perm, limit;
{
  generator i, j;

  for (i = 0; i < limit && perm[i] == NumGen; i++) ;
  if (i == limit) return NO;
  perm[i]++;
  for (j = 0; j < i; j++) perm[j] = 1;
  return YES;
}

/* This procedure inputs all the relations given by the permutation perm.
   Look also at the function Class_Poly in the ring module. */
void Get_Class_Mod(perm, limit)
generator *perm, limit;
{
  Poly *m, f;
  generator i, j;

  f = Class_Poly(perm, limit);
  m = New_Rule(&RUse);
  for (j = 1; j < NumMod; j++) m[j] = NULL;
  m[0] = f;
  Fill_Rule(&RUse, RUse.current);
  for (i = 1; i < NumMod; i++) {
    m = New_Rule(&RUse);
    for (j = 0; j < NumMod; j++) m[j] = NULL;
    m[i] = Copy_List(f);
    Fill_Rule(&RUse, RUse.current);
  }
}

/* This procedure uses the rules to reduce themselves. It returns the
   number of rules modified. */
counter Reduce_Among_Rules()
{
  generator i1, i2, start;
  Poly *m;
  flag newhead, modified;
  counter k;

  k = 0;
  i1 = 1;
  while (i1 <= RUse.current) {
    if (!(RUse.Active[i1-1] & (utilized | deleted))) {
      m = Get_Rule(&RUse, i1);
      modified = NO;
      i2 = 1;
      while (i2 <= RUse.current) {
        if (i1 != i2 && avail(i2-1)) {
          newhead = Mod_Reduce_Mod(m, Get_Rule(&RUse, i2));
          if (newhead == 3) {
            for (start = 0; start < NumMod && m[start] == NULL; start++);
            if (start == NumMod) {
              RUse.Active[i1-1] = deleted;
              i2 = RUse.current+1;
            }
            else i2 = 1;
          }
          else i2++;
          if (newhead != NO) modified = YES;
        }
        else i2++;
      };
      if (modified == YES) {
        k++;
        Fill_Rule(&RUse, i1);
        if (RUse.Length[i1-1] == 0) RUse.Active[i1-1] = deleted;
      }
    };
    i1++;
  };
  return k;
}

