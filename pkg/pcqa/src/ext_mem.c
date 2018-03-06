/*
  Written by: Eddie Lo
  Date started: November 10, 1995.

  Part of the polycyclic quotient algorithm package.

  This program manages memory for the extend module.
*/

#include <stdlib.h>
#include "pcqa.h"

Poly **ppprel = NULL;
Poly **powerprel = NULL;
Poly **ppowerrel = NULL;
Poly **powerrel = NULL;
Poly **drel = NULL;
Poly **relrel = NULL;
counter *table3 = NULL;
Poly *extend_m1 = NULL;
Poly *extend_m2 = NULL;
Poly *make_m1 = NULL;
Poly *make_m2 = NULL;
element extend_u = NULL;
counter used_rule = 0;
counter Cur_Rule = 0;
flag *ppp_found = NULL;
flag *powerp_found = NULL;
flag *ppower_found = NULL;
flag *power_found = NULL;
flag *d_found = NULL;
flag *rel_found = NULL;

extern generator NumGen, NumGenP, NumMod;
extern counter Size_R, choose2, choose3;
extern flag ExistsExtend, ExistsFound;
extern counter *table2;


/* This procedure initializes for the set of found flas. */
void Init_Found()
{
  generator i;

  ppp_found = (flag *) calloc(choose3, sizeof(flag));
  for (i = 0; i < choose3; i++) ppp_found[i] = NO;
  powerp_found = (flag *) calloc(choose2, sizeof(flag));
  for (i = 0; i < choose2; i++) powerp_found[i] = NO;
  ppower_found = (flag *) calloc(choose2, sizeof(flag));
  for (i = 0; i < choose2; i++) ppower_found[i] = NO;
  power_found = (flag *) calloc(NumGenP, sizeof(flag));
  for (i = 1; i <= NumGen; i++) power_found[i] = NO;
  d_found = (flag *) calloc(NumGenP, sizeof(flag));
  for (i = 1; i <= NumGen; i++) d_found[i] = NO;
  rel_found = (flag *) calloc(Size_R, sizeof(flag));
  for (i = 0; i < Size_R; i++) rel_found[i] = NO;
  ExistsFound = YES;
}

/* This procedure initializes for the extend module. */
void Init_Extend()
{
  generator i;

  if (NumGen > 2) {
    table3 = (counter *) calloc(NumGen-2, sizeof(counter));
    table3[0] = 0;
    for (i = 1; i < NumGen-2; i++)
      table3[i] = table3[i-1]+(NumGen-i)*(NumGen-i-1)/2;
  };
  extend_m1 = (Poly *) calloc(2*NumMod, sizeof(Poly));
  extend_m2 = extend_m1+NumMod;
  extend_u = (element) calloc(NumGen, sizeof(exponent))+NumGen;
  make_m1 = (Poly *) calloc(2*NumMod, sizeof(Poly));
  make_m2 = make_m1+NumMod;
  ppprel = (Poly **) calloc(choose3, sizeof(Poly *));
  powerprel = (Poly **) calloc(choose2, sizeof(Poly *));
  ppowerrel = (Poly **) calloc(choose2, sizeof(Poly *));
  powerrel = (Poly **) calloc(NumGenP, sizeof(Poly *));
  drel = (Poly **) calloc(NumGenP, sizeof(Poly *));
  relrel = (Poly **) calloc(Size_R, sizeof(Poly *));
  ExistsExtend = YES;
}

/* This procedure frees up the memory used for the found flags. */
void Reset_Found()
{
  free(ppp_found); free(powerp_found); free(ppower_found); free(power_found);
  free(d_found); free(rel_found);
  ExistsFound = NO;
}

/* This procedure frees up the memory used for the extend module. */
void Reset_Extend()
{
  counter k;

  if (NumGen > 2) free(table3);
  Cur_Rule = used_rule = 0;
  free(extend_m1); free(extend_u-NumGen);
  free(make_m1);
  free(ppprel); free(powerprel); free(ppowerrel); free(powerrel);
  free(drel); free(relrel);
  ExistsExtend = NO;
}

