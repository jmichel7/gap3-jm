/*
  Written by: Eddie Lo
  Date started: December 3, 1995.

  Part of the polycyclic quotient package.

  This is the module for homomorphism. Here, homomorphism
  means an epimorphism from the group presentation to a polycyclic
  presentation. Suppose we have a finite presentation and a consistent
  polycyclic presentation. It is necessary to know how the generators
  in the finite presentation can be expressed as a collected word in
  the polycyclic presentation. Conversely, it would be helpful to know
  how a generator in the polycyclic presentation can be expressed in
  terms of the generators in the finite presentation. This module
  provides tools for expressing these relations.

Data definition:
  Pre_Poly: This specifies the epimorphism from the group presentation
            to the polycyclic presentation. Suppose there are Size_S
            generators, then Pre_Poly consists of Size_S slots of
            data type element.
  Poly_Pre: This specifies how the polycyclic generators can be expressed
            in terms of the group generators. The data structure for
            this is expressed as blocks of relterms. Please refer to
            the corresponding modules of relterms for details. For
            example, relator.c.
*/

#include <stdlib.h>
#include "pcqa.h"

counter Sum_H = 0;
element *Pre_Poly = NULL;
relterm **Poly_Pre = NULL;

extern counter Size_S;
extern generator NumGen;
extern flag ExistsHomom;


/* This procedure prints out a quotient. */
void Print_Homom(fo)
FILE *fo;
{
  generator i;
  counter k;
  relterm *p;

  for (k = 0; k < Size_S; k++) {
    fprintf(fo, "g%lu = ", k+1);
    Print_Elm(fo, Pre_Poly[k]);
    fprintf(fo, "\n");
  };
  fprintf(fo, "\n");
  for (i = 1; i <= NumGen; i++) {
    p = Poly_Pre[i-1];
    fprintf(fo, "a%hu = ", i);
    while (p->g != 0) {
      Print_relterm(fo, *p);
      p++;
    };
    if (p == Poly_Pre[i-1]) fprintf(fo, "1\n");
    else fprintf(fo, "\n");
  }
}

/* This procedure saves a quotient into the file fo. */
void Save_Homom(fo)
FILE *fo;
{
  counter k;
  relterm *p;

  for (k = 0; k < Size_S; k++) Save_Elm(fo, Pre_Poly[k]);
  fprintf(fo, "%lu\n", Sum_H);
  for (p = Poly_Pre[0], k = 0; k < Sum_H; p++, k++) {
    if (p->g == 0) fprintf(fo, "0\n");
    else fprintf(fo, "%hu %ld ", p->g, p->pow);
  };
  fprintf(fo, "\n");
}

/* This procedure saves a quotient into the file fo. */
void Save_GAP_Homom(fo)
FILE *fo;
{
  counter k;
  relterm *p;

  fprintf(fo, "PCQAhom := rec(\nEpimorphism := [\n");
  for (k = 0; k < Size_S; k++) {
    if (k != 0) fprintf(fo, ",");
    Save_GAP_Elm(fo, Pre_Poly[k]);
  };
  fprintf(fo, "],\nInverseMap := [\n");
  if (Sum_H > 0) fprintf(fo, "[");
  for (p = Poly_Pre[0], k = 0; k < Sum_H; p++, k++) {
    if (p->g == 0) {
      if (k < Sum_H-1) fprintf(fo, "],\n[");
      else fprintf(fo, "]\n");
    }
    else fprintf(fo, "[%hu,%ld],", p->g, p->pow);
  };
  fprintf(fo, "]\n);\n");
}

/* Read a quotient homomorphism into storage. */
void Read_Homom(fi)
FILE *fi;
{
  counter k;
  generator i;
  relterm *p;

  Pre_Poly = (element *) calloc(Size_S, sizeof(element));
  Pre_Poly[0] = (element) calloc(Size_S*NumGen, sizeof(exponent))+NumGen;
  Read_Elm(fi, Pre_Poly[0]);
  for (k = 1; k < Size_S; k++) {
    Pre_Poly[k] = Pre_Poly[k-1]+NumGen;
    Read_Elm(fi, Pre_Poly[k]);
  };
  fscanf(fi, "%lu", &Sum_H);
  Poly_Pre = (relterm **) calloc(NumGen, sizeof(relterm *));
  Poly_Pre[0] = (relterm *) calloc(Sum_H, sizeof(relterm));
  for (p = Poly_Pre[0], k = 0, i = 1; k < Sum_H; p++, k++) {
    fscanf(fi, "%hu", &(p->g));
    if (p->g == 0 && i < NumGen) {
      Poly_Pre[i] = p+1; i++;
    }
    else fscanf(fi, "%ld", &(p->pow));
  };
  ExistsHomom = YES;
}

/* Check whether a definition for homomorphism exists and act
   accordingly. If ind is 1 and a homomorphism exists, then the
   user is prompted whether to delete homomorphism. If ind is 0
   and a homomorphism does not exist, then the user is prompted
   whether to read a homomorphism. If the answer is no, 0 is
   returned. 1 is returned otherwise. */
flag Check_Homom(ind)
flag ind;
{
  if (ind == 1) {
    if (ExistsHomom == YES) {
      printf("Homomorphism exists. Erase? ");
      if (answer()) {
        Reset_Homom();
        return YES;
      }
    }
    else return YES;
  }
  else {
    if (ExistsHomom == NO) {
      printf("No homomorphism exists. Read? ");
      if (answer()) {
        Init_Homom();
        Input_Homom();
        return YES;
      }
    }
    else return YES;
  };
  return NO;
}

/* Delete storage for a quotient. */
void Reset_Homom()
{
  free(Pre_Poly[0]-NumGen);
  free(Pre_Poly);
  Sum_H = 0;
  free(Poly_Pre[0]);
  free(Poly_Pre);
  ExistsHomom = NO;
}

/* Initialize for the quotient. */
void Init_Homom()
{
  Pre_Poly = (element *) calloc(Size_S, sizeof(element));
  Pre_Poly[0] = (element) calloc(Size_S*NumGen, sizeof(exponent))+NumGen;
  Poly_Pre = (relterm **) calloc(NumGen, sizeof(relterm *));
  Poly_Pre[0] = (relterm *) calloc(Sum_H, sizeof(relterm));
  ExistsHomom = YES;
}

/* This procedure prints out the definition of the polycyclic generator with
   index i. */
void Print_Definition(fo, i)
FILE *fo;
generator i;
{
  relterm *p;

  fprintf(fo, "a%d = ", i);
  p = Poly_Pre[i];
  if (p->g == 0) {
    fprintf(fo, "1\n");
    return ;
  };
  for (; p->g != 0; p++) Print_relterm(fo, *p);
  fprintf(fo, "\n");
}

/* Print out the definitions of all the generators. */
void Print_All_Definition(fo)
FILE *fo;
{
  generator i;

  for (i = 1; i <= NumGen; i++) Print_Definition(fo, i);
}

