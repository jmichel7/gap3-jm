/*
  Written by: Eddie Lo
  Date started: November 21,1994.

  Part of the polycyclic quotient algorithm program.

  This program finds the abelian quotient.

Introduction:

  Given some relators in RelGen, this module computes the abelian
  quotient and its polycyclic presentation. The variables
  posaut, negaut, Power_Elm, Power_Inv_Elm, finite_index and
  definition will be initialized accordingly.

Data description:

  row stores the pointers to the rows.
*/

#include <stdlib.h>
#include "pcqa.h"

extern exponent *test_elm;
extern exponent **row, *expterm;
extern generator *nzterm, *corr;
extern generator NumGen, NumGenP;
extern exponent *Cpp;
extern autlist *posaut, *negaut;
extern exponent *Power_Elm, *Power_Inv_Elm, *finite_index;
extern relterm **RelGen;
extern counter Size_R, Size_S, Sum_H;
extern element *Pre_Poly;
extern relterm **Poly_Pre;
extern flag ExistsPoly;


/* This procedure stores the relations. */
void Store_Row()
{
  generator i;
  relterm *p;

  for (i = 0; i < Size_R; i++) {
    p = RelGen[i];
    while (p->g != 0) {
      row[i][p->g-1] += p->pow; p++;
    }
  }
}

/* This procedure generates a presentation from the relation matrix. */
void Abelian_Presentation()
{
  generator i, j1, j2, j3;
  element u;
  counter disp, k, sum;
  relterm *p;

  NumGen = Find_misc(row, Size_R, Size_S);
  NumGenP = NumGen+1;
  finite_index = (exponent *) calloc(NumGenP, sizeof(exponent));
  for (i = 1; i <= NumGen; i++) finite_index[i] = expterm[corr[i]];
  if (NumGen > 1) {
    posaut = (autlist *) calloc(NumGen-1, sizeof(autlist));
    negaut = (autlist *) calloc(NumGen-1, sizeof(autlist));
  };
  for (i = 0; i < NumGen-1; i++) {
    posaut[i] = (autstorage *) calloc(MaxAutStore, sizeof(autstorage));
    negaut[i] = (autstorage *) calloc(MaxAutStore, sizeof(autstorage));
    for (k = 1; k < MaxAutStore; k++) {
      posaut[i][k] = NULL;
      negaut[i][k] = NULL;
    }
  };
  sum = 2*NumGen*NumGen*NumGen;
  Cpp = (exponent *) calloc(sum, sizeof(exponent));
  for (k = 0; k < sum; k++) Cpp[k] = 0;
  Power_Elm = (exponent *) Cpp+2*(NumGen-1)*NumGen*NumGen;
  Power_Inv_Elm = (exponent *) Power_Elm+NumGen*NumGen;
  for (j1 = 0, i = 1; i <= NumGen; i++)
    if (finite_index[i] > 1) {
      for (; nzterm[j1] != corr[i]; j1++) ;
      for (j2 = 0; j2 < Size_S; j2++) test_elm[j2] = -row[j1][j2];
      test_elm[corr[i]] = 0;
      Abelian_Reduce(row, Size_R, Size_S, test_elm);
      u = (element) Power_Elm+i*NumGen;
      Fill_Vec(NumGen, test_elm, (exponent *) u);
      u = (element) Power_Inv_Elm+i*NumGen;
      Fill_Vec(NumGen, row[j1], (exponent *) u);
      u[-i] = finite_index[i]-1;
      j1++;
    };
  disp = 0;
  for (i = 1; i < NumGen; i++) {
    *(posaut[i-1]) = (autstorage) Cpp+disp;
    if (!finite_index[i]) {
      u = (element) *(posaut[i-1])+(NumGen-i)*NumGen;
      for (j1 = i+1; j1 <= NumGen; j1++) {
        if (finite_index[j1] != 1) {
          for (j2 = NumGen; j2 >= 1; j2--) u[-j2] = 0;
          u[-j1] = 1;
        }
        else {
          for (j2 = 1; j2 <= j1; j2++) u[-j2] = 0;
          for (; j2 <= NumGen; j2++) u[-j2] = Power_Elm[j1*NumGen-j2];
        };
        u -= NumGen;
      }
    };
    disp += 2*(NumGen-i)*NumGen;
  };
  for (i = 1; i < NumGen; i++) {
    if (!finite_index[i]) {
      u = (element) *(posaut[i-1])+2*(NumGen-i)*NumGen;
      for (j1 = i+1; j1 <= NumGen; j1++) {
        for (j2 = NumGen; j2 >= 1; j2--) u[-j2] = 0;
        if (!finite_index[j1]) u[-j1] = -1;
        u -= NumGen;
      }
    }
  };
  for (i = 1; i < NumGen; i++) {
    *(negaut[i-1]) = (autstorage) Cpp+disp;
    u = (element) *(negaut[i-1])+(NumGen-i)*NumGen;
    for (j1 = i+1; j1 <= NumGen; j1++) {
      if (finite_index[j1] != 1) {
        for (j2 = NumGen; j2 >= 1; j2--) u[-j2] = 0;
        u[-j1] = 1;
      }
      else {
        for (j2 = 1; j2 <= j1; j2++) u[-j2] = 0;
        for (; j2 <= NumGen; j2++) u[-j2] = Power_Elm[j1*NumGen-j2];
      };
      u -= NumGen;
    };
    disp += 2*(NumGen-i)*NumGen;
  };
  for (i = 1; i < NumGen; i++) {
    u = (element) *(negaut[i-1])+2*(NumGen-i)*NumGen;
    for (j1 = i+1; j1 <= NumGen; j1++) {
      for (j2 = NumGen; j2 >= 1; j2--) u[-j2] = 0;
      if (!finite_index[j1]) u[-j1] = -1;
      u -= NumGen;
    }
  };
  Init_Group();

  if (Check_Homom(1) == YES) {
    Sum_H = 2*NumGen;
    Init_Homom();
    for (i = 1; i < NumGen; i++) Poly_Pre[i] = Poly_Pre[i-1]+2;
    for (i = 1, p = Poly_Pre[0]; i <= NumGen; i++) {
      p->g = corr[i]+1; p->pow = 1; p++;
      p->g = 0; p++;
    };
    for (i = 1; i < Size_S; i++) Pre_Poly[i] = Pre_Poly[i-1]+NumGen;
    for (i = 0, j1 = 1, j3 = 0; i < Size_S; i++) {
      u = Pre_Poly[i];
      if (expterm[i] == 1) {
        for (j2 = 0; j2 < Size_S; j2++) test_elm[j2] = -row[j3][j2];
        test_elm[i] = 0;
        Abelian_Reduce(row, Size_R, Size_S, test_elm);
        for (j2 = 1; j2 <= NumGen; j2++) u[-j2] = test_elm[corr[j2]];
        j3++;
      }
      else {
        for (j2 = 1; j2 < j1; j2++) u[-j2] = 0;
        u[-j1] = 1;
        for (j2 = j1+1; j2 <= NumGen; j2++) u[-j2] = 0;
        j1++;
        if (expterm[i] > 0) j3++;
      }
    }
  }
}

/* This procedure computes the abelian quotient. */
flag Get_Abelian()
{
  if (Check_Poly_Presentation(1) == NO) return NO;
  if (Check_Parsed_Presentation(1) == NO) return NO;
  if (Input_Raw_Presentation() == NO) return NO;
  Init_Abelian();
  Init_HNF(Size_R, Size_S);
  Store_Row();
  Abelian_Row_Reduce(row, Size_R, Size_S);
  Abelian_Presentation();
  Reset_HNF();
  Reset_Abelian();
  return YES;
}

