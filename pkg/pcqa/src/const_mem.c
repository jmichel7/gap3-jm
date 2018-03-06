/*
  Written by: Eddie Lo
  Date started: October 24, 1995

  Part of the polycyclic quotient algorithm package.

  This program manages memory for the construct module.
*/

#include <stdlib.h>
#include "pcqa.h"

exponent *posmax = NULL;
exponent *negmax = NULL;
flag *check = NULL;
counter *range1 = NULL;
counter *range2 = NULL;

element construct_u = NULL;
Poly *construct_m = NULL;
element modelm_u = NULL;
Poly *modelm_m = NULL;
Poly *map_m = NULL;
coeff map_num = NULL;
cellptr *fill_h = NULL;
element fill_u = NULL;
Poly *fill_m = NULL;

exponent **RelPow = NULL;
exponent *RPmem = NULL;
exponent *Test_Row = NULL;
exponent *modelm_row = NULL;
counter Test_Gen = 0;
counter relpow = 0;
rel **ModStr = NULL;

generator *new_start = NULL;
generator *new_perm = NULL;
generator *new_invperm = NULL;
exponent *new_element = NULL;
exponent *new_power = NULL;
Poly **new_mod = NULL;
generator New_Gen = 0;
generator new_NumGen = 0;
generator new_NumGenP = 0;
/* exponent *new_power = NULL; */
autlist *new_posaut = NULL;
autlist *new_negaut = NULL;
exponent *new_Cpp = NULL;
exponent *new_Power_Elm = NULL;
exponent *new_Power_Inv_Elm = NULL;
exponent *new_finite_index = NULL;
counter new_Sum_H = 0;
element *new_Pre_Poly = NULL;
relterm **new_Poly_Pre = NULL;
rel **new_RelDef = NULL;

extern generator NumMod, NumGen, NumGenP;
extern autlist *posaut, *negaut;
extern exponent *Cpp;
extern exponent *Power_Elm, *Power_Inv_Elm, *finite_index;
extern counter Size_S, Sum_H, Num_Basis;
extern element *Pre_Poly;
extern relterm **Poly_Pre;

extern flag ExistsHead, ExistsCollect, ExistsNumMod, ExistsModDef;
extern flag ExistsRule, ExistsExtend, ExistsFound, ExistsConstruct;
extern flag ExistsParsed, ExistsPoly, ExistsHomom, ExistsNewPoly, ExistsMemory;
extern flag ExistsReader, ExistsQuotient;


/* This procedure establishes memory for the construct module. */
void Est_Construct()
{
  map_num = (coeff) malloc(sizeof(MP_INT));
  mpz_init(map_num);
}

/* This procedure initializes for the construct module, used many times. */
void Init_Construct()
{
  generator i;

  posmax = (exponent *) calloc(NumGen*NumMod, sizeof(exponent));
  negmax = (exponent *) calloc(NumGen*NumMod, sizeof(exponent));
  range1 = (counter *) calloc(NumMod+1, sizeof(counter));
  range2 = (counter *) calloc(NumMod+1, sizeof(counter));
  ModStr = (rel **) calloc(NumMod, sizeof(rel *));

  construct_u = (element) calloc(NumGen, sizeof(exponent))+NumGen;
  construct_m = (Poly *) calloc(NumMod, sizeof(Poly));
  map_m = (Poly *) calloc(NumMod, sizeof(Poly));
  fill_h = (cellptr *) calloc(NumGenP, sizeof(cellptr));
  fill_u = (element) calloc(NumGen, sizeof(exponent))+NumGen;
  for (i = 1; i <= NumGen; i++) fill_u[-i] = 0;
  fill_m = (Poly *) calloc(NumMod, sizeof(Poly));
  ExistsConstruct = YES;
}

/* This procedure frees up memory for the construct module. */
void Reset_Construct()
{
  free(posmax); free(negmax);
  free(range1); free(range2); free(ModStr);
  free(construct_u-NumGen); free(construct_m); free(modelm_m);
  free(fill_h); free(fill_u-NumGen); free(fill_m);
  ExistsConstruct = NO;
}

/* This procedure initializes for the ModElm module. */
void Init_ModElm()
{
  modelm_u = (element) calloc(NumGen, sizeof(exponent))+NumGen;
  modelm_m = (Poly *) calloc(NumMod, sizeof(Poly));
  modelm_row = (exponent *) calloc(Test_Gen, sizeof(exponent));
}

/* This procedure resets for the ModElm module. */
void Reset_ModElm()
{
  free(modelm_u-NumGen);
  free(modelm_m);
  free(modelm_row);
}

/* This procedure initializes the check flags. */
void Init_check()
{
  counter k;

  check = (flag *) calloc(Num_Basis+1, sizeof(flag));
  for (k = 1; k <= Num_Basis; k++) check[k] = 0;
}

/* This procedure resets the check flags. */
void Reset_check()
{
  free(check);
}

/* This procedure initializes a set of generators. */
void Init_Generator()
{
  generator i, j;

  new_start = (generator *) calloc(Test_Gen, sizeof(generator));
  new_element = (exponent *) calloc(Test_Gen*NumGen, sizeof(exponent));
  new_perm = (generator *) calloc(Test_Gen, sizeof(generator));
  new_invperm = (generator *) calloc(Test_Gen, sizeof(generator));
  new_mod = (Poly **) calloc(Test_Gen, sizeof(Poly *));
  new_mod[0] = (Poly *) calloc(Test_Gen*NumMod, sizeof(Poly));
  for (i = 0; i < Test_Gen; i++) {
    new_mod[i] = new_mod[0]+i*NumMod;
    for (j = 0; j < NumMod; j++) new_mod[i][j] = NULL;
  };
  new_power = (exponent *) calloc(Test_Gen, sizeof(exponent));
}

/* Reset the definition of the generators. */
void Reset_Generator()
{
  generator i;

  free(new_start);
  free(new_element);
  free(new_perm); free(new_invperm);
  for (i = 0; i < New_Gen; i++) Free_Mod_Cell(new_mod[i], 0);
  free(new_mod[0]);
  free(new_mod);
  free(new_power);
}

/* This procedure initializes a set of rows. */
void Init_RelPow()
{
  counter k;

  if (relpow > 0) {
    RelPow = (exponent **) calloc(relpow, sizeof(exponent *));
    RPmem = (exponent *) calloc(relpow*Test_Gen, sizeof(exponent));
    for (k = 0; k < relpow; k++) RelPow[k] = RPmem+k*Test_Gen;
  };
  Test_Row = (exponent *) calloc(Test_Gen, sizeof(exponent));
}

/* This procedure resets a set of rows. */
void Reset_RelPow()
{
  if (relpow > 0) {
    free(RPmem);
    free(RelPow);
  };
  free(Test_Row);
}

/* This procedure initializes memory for a new presentation. */
exponent *Init_New_Presentation()
{
  generator i;
  counter k, disp, sum;

  new_NumGen = NumGen+New_Gen;
  new_NumGenP = new_NumGen+1;
  new_finite_index = (exponent *) calloc(new_NumGenP, sizeof(exponent));
  new_posaut = (autlist *) calloc(new_NumGen-1, sizeof(autlist));
  new_negaut = (autlist *) calloc(new_NumGen-1, sizeof(autlist));
  for (i = 0; i < new_NumGen-1; i++) {
    new_posaut[i] = (autstorage *) calloc(MaxAutStore, sizeof(autstorage));
    new_negaut[i] = (autstorage *) calloc(MaxAutStore, sizeof(autstorage));
    for (k = 1; k < MaxAutStore; k++) {
      new_posaut[i][k] = NULL;
      new_negaut[i][k] = NULL;
    }
  };
  sum = 2*new_NumGen*new_NumGen*new_NumGen;
  new_Cpp = (exponent *) calloc(sum, sizeof(exponent));
  for (k = 0; k < sum; k++) new_Cpp[k] = 0;
  disp = 0;
  for (i = 1; i < NumGen; i++) {
    new_posaut[i-1][0] = (autstorage) new_Cpp+disp;
    disp += 2*(NumGen-i)*NumGen;
  };
  for (i = 1; i < NumGen; i++) {
    *(new_negaut[i-1]) = (autstorage) new_Cpp+disp;
    disp += 2*(NumGen-i)*NumGen;
  };
  new_Power_Elm = (exponent *) new_Cpp+2*(new_NumGen-1)*new_NumGen*new_NumGen;
  new_Power_Inv_Elm = (exponent *) new_Power_Elm+new_NumGen*new_NumGen;
  return new_Cpp;
}

/* Initialize storage for maps between the polycyclic quotient and the
   presentation. */
void Init_Map()
{
  counter k;

  new_Pre_Poly = (element *) calloc(Size_S, sizeof(element));
  new_Pre_Poly[0] = (element) calloc(Size_S*new_NumGen, sizeof(exponent))+
                      new_NumGen;
  for (k = 1; k < Size_S; k++)
    new_Pre_Poly[k] = new_Pre_Poly[k-1]+new_NumGen;
  new_Poly_Pre = (relterm **) calloc(new_NumGen, sizeof(relterm *));
  new_Poly_Pre[0] = (relterm *) calloc(new_Sum_H, sizeof(relterm));
}

/* Initialize storage for relators defining the polycyclic generators. */
void Init_RelDef()
{
  new_RelDef = (rel **) calloc(New_Gen, sizeof(rel *));
}

/* Free storage used for defining the polycyclic generators. */
void Reset_RelDef()
{
  generator i;

  for (i = 0; i < New_Gen; i++) Free_Str(new_RelDef[i]);
  free(new_RelDef);
}

/* Delete a new presentation. */
void Delete_New_Poly_Presentation()
{
  generator i;
  counter k;

  free(new_finite_index);
  free(new_Cpp);
  for (i = 0; i < new_NumGen-1; i++) {
    free(new_posaut[i]); free(new_negaut[i]);
  };
  free(new_posaut); free(new_negaut);
  free(new_Pre_Poly[0]-new_NumGen);
  free(new_Pre_Poly);
  free(new_Poly_Pre[0]);
  free(new_Poly_Pre);
  ExistsNewPoly = NO;
}

/* Delete all the old stuff once a new polycyclic presentation is
   initialized. */
void Delete_Old()
{
  if (ExistsHead == YES) Reset_Head();
  if (ExistsCollect == YES) Reset_Collect();
  if (ExistsNumMod != NO) {
    NumMod = 0;
    ExistsNumMod = NO;
  };
  if (ExistsModDef != NO) Reset_ModDef();
  if (ExistsRule == YES) Reset_Rule();
  if (ExistsExtend == YES) Reset_Extend();
  if (ExistsFound != NO) Reset_Found();
  if (ExistsConstruct == YES) Reset_Construct();
  if (ExistsHomom == YES) Reset_Homom();
  if (ExistsReader == YES) Reset_Reader();
  if (ExistsMemory == YES) {
    Reset_Tool();
    Reset_Memory();
  };
  if (ExistsQuotient == YES) Reset_Quotient();
  if (ExistsPoly == YES) Reset_Group();
}

/* Reinitialize variables for the new group. */
void Replace()
{
  generator i;
  counter k;

  NumGen = new_NumGen;
  NumGenP = new_NumGenP;
  finite_index = new_finite_index;
  posaut = new_posaut;
  negaut = new_negaut;
  Cpp = new_Cpp;
  Power_Elm = new_Power_Elm;
  Power_Inv_Elm = new_Power_Inv_Elm;
  Init_Group();
  Pre_Poly = new_Pre_Poly;
  Poly_Pre = new_Poly_Pre;
  Sum_H = new_Sum_H;
  ExistsHomom = YES;
  ExistsNewPoly = NO;
}
