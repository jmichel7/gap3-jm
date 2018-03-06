/*
  Written by: Eddie Lo
  Date started: December 3, 1995.

  Part of the polycyclic quotient algorithm package.

  This consists of tools to work with a polycyclic presentation.
  For detail of a polycyclic presentation, please refer to group.c.
*/

#include <stdlib.h>
#include "pcqa.h"

extern generator NumGen, NumGenP;
extern autlist *posaut, *negaut;
extern exponent *Cpp;
extern exponent *Power_Elm, *Power_Inv_Elm, *finite_index;
extern flag ExistsPoly;
extern element id;


/* Print out a polycyclic presentation. */
void Print_Poly_Presentation(fo)
FILE *fo;
{
  generator i, j;
  element u;

  fprintf(fo, "Number of generators = %hu\n", NumGen);
  for (i = 1; i < NumGen; i++) {
    u = (element) *(negaut[i-1])+(NumGen-i)*NumGen;
    for (j = i+1; j <= NumGen; j++) {
      fprintf(fo, "a%hu a%hu = ", i, j);
      Print_Elm(fo, u);
      fprintf(fo, "a%hu\n", i);
      u -= NumGen;
    }
  };
  for (i = 1; i < NumGen; i++) {
    u = (element) *(negaut[i-1])+2*(NumGen-i)*NumGen;
    for (j = i+1; j <= NumGen; j++) {
      if (finite_index[j] == 0) {
        fprintf(fo, "a%hu A%hu = ", i, j);
        Print_Elm(fo, u);
        fprintf(fo, "a%hu\n", i);
      }
      u -= NumGen;
    }
  };
  for (i = 1; i < NumGen; i++) {
    if (finite_index[i] == 0) {
      u = (element) *(posaut[i-1])+(NumGen-i)*NumGen;
      for (j = i+1; j <= NumGen; j++) {
        fprintf(fo, "A%hu a%hu = ", i, j);
        Print_Elm(fo, u);
        fprintf(fo, "A%hu\n", i);
        u -= NumGen;
      }
    }
  };
  for (i = 1; i < NumGen; i++) {
    if (finite_index[i] == 0) {
      u = (element) *(posaut[i-1])+2*(NumGen-i)*NumGen;
      for (j = i+1; j <= NumGen; j++) {
        if (finite_index[j] == 0) {
          fprintf(fo, "A%hu A%hu = ", i, j);
          Print_Elm(fo, u);
          fprintf(fo, "A%hu\n", i);
        };
        u -= NumGen;
      }
    }
  };
  for (j = 1; j <= NumGen; j++) {
    if (finite_index[j] != 0) {
      u = (element) Power_Elm+j*NumGen;
      fprintf(fo, "a%hu^%ld = ", j, finite_index[j]);
      Print_Elm(fo, u);
      fprintf(fo, "\n");
      u = (element) Power_Inv_Elm+j*NumGen;
      fprintf(fo, "A%hu = ", j);
      Print_Elm(fo, u);
      fprintf(fo, "\n");
    }
  }
}

/* Save a polycyclic presentation to a file fo. */
void Save_Poly_Presentation(fo)
FILE *fo;
{
  generator i, j1, j2;
  element u;

  fprintf(fo, "%hu\n", NumGen);
  for (i = 1; i <= NumGen; i++) fprintf(fo, "%u ", finite_index[i]);
  fprintf(fo, "\n\n");
  for (i = 1; i < NumGen; i++) {
    u = (element) *(posaut[i-1])+(NumGen-i)*NumGen;
    for (j1 = i+1; j1 <= NumGen; j1++) {
      for (j2 = NumGen; j2 >= 1; j2--) fprintf(fo, "%ld ", u[-j2]);
      fprintf(fo, "\n");
      u -= NumGen;
    }
  };
  fprintf(fo, "\n");
  for (i = 1; i < NumGen; i++) {
    u = (element) *(posaut[i-1])+2*(NumGen-i)*NumGen;
    for (j1 = i+1; j1 <= NumGen; j1++) {
      for (j2 = NumGen; j2 >= 1; j2--) fprintf(fo, "%ld ", u[-j2]);
      fprintf(fo, "\n");
      u -= NumGen;
    }
  };
  fprintf(fo, "\n");
  for (i = 1; i < NumGen; i++) {
    u = (element) *(negaut[i-1])+(NumGen-i)*NumGen;
    for (j1 = i+1; j1 <= NumGen; j1++) {
      for (j2 = NumGen; j2 >= 1; j2--) fprintf(fo, "%ld ", u[-j2]);
      fprintf(fo, "\n");
      u -= NumGen;
    }
  };
  fprintf(fo, "\n");
  for (i = 1; i < NumGen; i++) {
    u = (element) *(negaut[i-1])+2*(NumGen-i)*NumGen;
    for (j1 = i+1; j1 <= NumGen; j1++) {
      for (j2 = NumGen; j2 >= 1; j2--) fprintf(fo, "%ld ", u[-j2]);
      fprintf(fo, "\n");
      u -= NumGen;
    }
  };
  fprintf(fo, "\n");
  for (j1 = 1; j1 <= NumGen; j1++)
    if (finite_index[j1]) {
      u = (element) Power_Elm+j1*NumGen;
      for (j2 = NumGen; j2 >= 1; j2--) fprintf(fo, "%ld ", u[-j2]);
      fprintf(fo, "\n");
      u = (element) Power_Inv_Elm+j1*NumGen;
      for (j2 = NumGen; j2 >= 1; j2--) fprintf(fo, "%ld ", u[-j2]);
      fprintf(fo, "\n\n");
    }
}

/* Save a polycyclic presentation to a file fo. */
void Save_GAP_Poly_Presentation(fo)
FILE *fo;
{
  generator i, j1, j2;
  element u;

  fprintf(fo, "PCQAcpp := rec (\nGenerators := %hu,\n", NumGen);
  fprintf(fo, "ExponentList := [");
  for (i = 1; i <= NumGen; i++) {
    fprintf(fo, "%u", finite_index[i]);
    if (i != NumGen) fprintf(fo, ",");
  };
  fprintf(fo, "],\nppRelations := [\n");
  for (i = 1; i < NumGen; i++) {
    fprintf(fo, "[\n");
    u = (element) *(posaut[i-1])+(NumGen-i)*NumGen;
    for (j1 = i+1; j1 <= NumGen; j1++) {
      if (j1 > i+1) fprintf(fo, ",");
      if (finite_index[i] == 0 || !Equal_Term(u, id)) Save_GAP_Elm(fo, u);
      u -= NumGen;
    };
    if (i < NumGen-1) fprintf(fo, "],\n");
    else fprintf(fo, "]\n");
  };
  fprintf(fo, "],\npnRelations := [\n");
  for (i = 1; i < NumGen; i++) {
    fprintf(fo, "[\n");
    u = (element) *(posaut[i-1])+2*(NumGen-i)*NumGen;
    for (j1 = i+1; j1 <= NumGen; j1++) {
      if (j1 > i+1) fprintf(fo, ",");
      if ((finite_index[i] == 0 && finite_index[j1] == 0) || !Equal_Term(u, id))
        Save_GAP_Elm(fo, u);
      u -= NumGen;
    };
    if (i < NumGen-1) fprintf(fo, "],\n");
    else fprintf(fo, "]\n");
  };
  fprintf(fo, "],\nnpRelations := [\n");
  for (i = 1; i < NumGen; i++) {
    fprintf(fo, "[\n");
    u = (element) *(negaut[i-1])+(NumGen-i)*NumGen;
    for (j1 = i+1; j1 <= NumGen; j1++) {
      if (j1 > i+1) fprintf(fo, ",");
      Save_GAP_Elm(fo, u);
      u -= NumGen;
    };
    if (i < NumGen-1) fprintf(fo, "],\n");
    else fprintf(fo, "]\n");
  };
  fprintf(fo, "],\nnnRelations := [\n");
  for (i = 1; i < NumGen; i++) {
    fprintf(fo, "[\n");
    u = (element) *(negaut[i-1])+2*(NumGen-i)*NumGen;
    for (j1 = i+1; j1 <= NumGen; j1++) {
      if (j1 > i+1) fprintf(fo, ",");
      if (finite_index[j1] == 0 || !Equal_Term(u, id)) Save_GAP_Elm(fo, u);
      u -= NumGen;
    };
    if (i < NumGen-1) fprintf(fo, "],\n");
    else fprintf(fo, "]\n");
  };
  fprintf(fo, "],\nPowerRelations := [\n");
  for (j1 = 1; j1 <= NumGen; j1++) {
    if (finite_index[j1]) {
      u = (element) Power_Elm+j1*NumGen;
      Save_GAP_Elm(fo, u);
      fprintf(fo, ",");
      u = (element) Power_Inv_Elm+j1*NumGen;
      Save_GAP_Elm(fo, u);
      fprintf(fo, ",");
    }
    else fprintf(fo, ",,");
  };
  fprintf(fo, "]\n);\n");
}

/* Read a poly presentation from a file fi. The file input has to be of the
   form:
   First line: NumGen = n.
   Second line: The orders k1,...,kn of the quotients. 0 is it is infinite.
   The conjugacy relations, in the order:
     a2^a1, a3^a1, ..., an^a(n-1),
     (a2^-1)^a1, (a3^-1)^a1, ..., (an^-1)^a(n-1)
     a2^(a1^-1), a3^(a1^-1), ..., an^(a(n-1)^-1)
     (a2^-1)^(a1^-1), (a3^-1)^(a1^-1), ..., (an^-1)^(a(n-1)^-1).
   The power relations: a1^k1, a2^k2, ..., an^kn, if exist. */
void Read_Poly_Presentation(fi)
FILE *fi;
{
  generator i, j1, j2;
  counter k, disp;
  element u;

  fscanf(fi, "%hu", &NumGen);
  NumGenP = NumGen+1;
  finite_index = (exponent *) calloc(NumGenP, sizeof(exponent));
  for (i = 1; i <= NumGen; i++) fscanf(fi, "%ld", finite_index+i);
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
  Cpp = (exponent *) calloc(2*NumGen*NumGen*NumGen, sizeof(exponent));
  disp = 0;
  for (i = 1; i < NumGen; i++) {
    *(posaut[i-1]) = (autstorage) Cpp+disp;
    u = (element) *(posaut[i-1])+(NumGen-i)*NumGen;
    for (j1 = i+1; j1 <= NumGen; j1++) {
      for (j2 = NumGen; j2 >= 1; j2--) fscanf(fi, "%ld", u-j2);
      u -= NumGen;
    };
    disp += 2*(NumGen-i)*NumGen;
  };
  for (i = 1; i < NumGen; i++) {
    u = (element) *(posaut[i-1])+2*(NumGen-i)*NumGen;
    for (j1 = i+1; j1 <= NumGen; j1++) {
      for (j2 = NumGen; j2 >= 1; j2--) fscanf(fi, "%ld", u-j2);
      u -= NumGen;
    }
  };
  for (i = 1; i < NumGen; i++) {
    *(negaut[i-1]) = (autstorage) Cpp+disp;
    u = (element) *(negaut[i-1])+(NumGen-i)*NumGen;
    for (j1 = i+1; j1 <= NumGen; j1++) {
      for (j2 = NumGen; j2 >= 1; j2--) fscanf(fi, "%ld", u-j2);
      u -= NumGen;
    };
    disp += 2*(NumGen-i)*NumGen;
  };
  for (i = 1; i < NumGen; i++) {
    u = (element) *(negaut[i-1])+2*(NumGen-i)*NumGen;
    for (j1 = i+1; j1 <= NumGen; j1++) {
      for (j2 = NumGen; j2 >= 1; j2--) fscanf(fi, "%ld", u-j2);
      u -= NumGen;
    }
  };
  Power_Elm = (exponent *) Cpp+2*(NumGen-1)*NumGen*NumGen;
  Power_Inv_Elm = (exponent *) Power_Elm+NumGen*NumGen;
  for (j1 = 1; j1 <= NumGen; j1++)
    if (finite_index[j1]) {
      u = (element) Power_Elm+j1*NumGen;
      for (j2 = NumGen; j2 >= 1; j2--) fscanf(fi, "%ld", u-j2);
      u = (element) Power_Inv_Elm+j1*NumGen;
      for (j2 = NumGen; j2 >= 1; j2--) fscanf(fi, "%ld", u-j2);
    }
}

/* Check whether a polycyclic presentation exists already and
   act appropriately. If ind is 1, it prompts the user whether
   to delete the polycyclic presentation if it exists. Otherwise,
   it prompts the user whether to read a polycyclic presentation
   if it doesn't already exists. If the answer is no in either
   case, it returns 0. Otherwise it returns 1. */
flag Check_Poly_Presentation(ind)
flag ind;
{
  if (ind == 1) {
    if (ExistsPoly == YES) {
      printf("Polycyclic presentation exists. Erase? ");
      if (answer()) {
        Delete_Old();
        return YES;
      }
    }
    else return YES;
  }
  else {
    if (ExistsPoly == NO) {
      printf("No polycyclic presentation exists. Read? ");
      if (answer()) {
        Input_Poly_Presentation();
        return YES;
      }
    }
    else return YES;
  };
  return NO;
}

/* Delete the storage for a consistent polycyclic presentation. */
void Delete_Poly_Presentation()
{
  generator i;
  counter k;

  free(finite_index);
  free(Cpp);
  for (i = 0; i < NumGen-1; i++) {
    for (k = 1; k < MaxAutStore && posaut[i][k] != NULL; k++)
      free(posaut[i][k]);
    for (k = 1; k < MaxAutStore && negaut[i][k] != NULL; k++)
      free(negaut[i][k]);
    free(posaut[i]); free(negaut[i]);
  };
  if (NumGen > 1) {
    free(posaut);
    free(negaut);
  }
}

