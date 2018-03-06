/*
  Written by: Eddie Lo
  Date started: November 10, 1995.

  Part of the polycyclic quotient algorithm package.

  This program manages memory for the head module.
*/

#include <stdlib.h>
#include "pcqa.h"

Poly **ppmod = NULL;
Poly **pnmod = NULL;
Poly **npmod = NULL;
Poly **nnmod = NULL;
Poly **powermod = NULL;
Poly **npowermod = NULL;
Poly **defmod = NULL;
generator *ppstart = NULL;
generator *pnstart = NULL;
generator *npstart = NULL;
generator *nnstart = NULL;
generator *powerstart = NULL;
generator *npowerstart = NULL;
generator *defstart = NULL;
def *module_definition = NULL;
flag *gg_defined = NULL;
flag *n_defined = NULL;
Poly *head_m1 = NULL;
Poly *head_m2 = NULL;
element head_u = NULL;
Poly *HeadMod_List = NULL;
counter used_headmod = 0;

extern generator NumGen, NumGenP, NumMod;
extern counter Size_S;
extern counter choose2;
extern flag ExistsHead, ExistsModDef;


/* This procedure initializes for definition of module generators. */
void Init_ModDef()
{
  module_definition = (def *) calloc(NumMod, sizeof(def));
  ExistsModDef = YES;
}

/* This procedure initializes for the head module. */
void Init_Head()
{
  generator i;
  counter k;

  HeadMod_List = (Poly *) calloc(Num_HeadMod()*NumMod, sizeof(Poly));
  used_headmod = 0;

  ppmod = (Poly **) calloc(choose2, sizeof(Poly *));
  pnmod = (Poly **) calloc(choose2, sizeof(Poly *));
  npmod = (Poly **) calloc(choose2, sizeof(Poly *));
  nnmod = (Poly **) calloc(choose2, sizeof(Poly *));
  powermod = (Poly **) calloc(NumGenP, sizeof(Poly *));
  npowermod = (Poly **) calloc(NumGenP, sizeof(Poly *));
  defmod = (Poly **) calloc(Size_S, sizeof(Poly *));

  ppstart = (generator *) calloc(choose2, sizeof(generator));
  pnstart = (generator *) calloc(choose2, sizeof(generator));
  npstart = (generator *) calloc(choose2, sizeof(generator));
  nnstart = (generator *) calloc(choose2, sizeof(generator));
  powerstart = (generator *) calloc(NumGenP, sizeof(generator));
  npowerstart = (generator *) calloc(NumGenP, sizeof(generator));
  defstart = (generator *) calloc(Size_S, sizeof(generator));

  gg_defined = (flag *) calloc(choose2, sizeof(flag));
  for (k = 0; k < choose2; k++) gg_defined[k] = 1;
  n_defined = (flag *) calloc(NumGenP, sizeof(flag));
  for (i = 1; i <= NumGen; i++) n_defined[i] = 0;

  head_m1 = (Poly *) calloc(2*NumMod, sizeof(Poly));
  head_m2 = head_m1+NumMod;
  head_u = (element) calloc(NumGen, sizeof(exponent))+NumGen;

  ExistsHead = YES;
}

/* This procedure frees up memory used for the head module. */
void Reset_Head()
{
  counter k;

  for (k = 0; k < used_headmod; k++) Free_Mod_Cell(HeadMod_List+k*NumMod, 0);
  free(HeadMod_List);
  used_headmod = 0;
  free(ppmod); free(pnmod); free(npmod); free(nnmod);
  free(powermod); free(npowermod); free(defmod);
  free(ppstart); free(pnstart); free(npstart); free(nnstart);
  free(powerstart); free(npowerstart); free(defstart);
  free(gg_defined); free(n_defined);
  free(head_m1); free(head_u-NumGen);
  ExistsHead = NO;
}

/* This procedure deletes storage for definition of module generators. */
void Reset_ModDef()
{
  free(module_definition);
  ExistsModDef = NO;
}

/* This procedure returns an Extend element to the calling procedure. */
Poly *New_HeadMod()
{
  Poly *m;

  m = HeadMod_List+used_headmod*NumMod;
  used_headmod++;
  return m;
}

