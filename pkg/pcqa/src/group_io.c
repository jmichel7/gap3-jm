/*
  Written by: Eddie Lo
  Date started: October 24, 1995

  Part of the polycyclic quotient algorithm package.

  This program contains procedures for interface in the group module.
*/

#include "pcqa.h"

extern generator NumGen;


/* Read an element from an input file. */
void Read_Elm(FILE *fi, element u)
{
  generator i;

  for (i = NumGen; i >= 1; i--) fscanf(fi, "%ld", u-i);
}

/* Save an element to a file fo. */
void Save_Elm(fo, u)
FILE *fo;
element u;
{
  generator i;

  for (i = NumGen; i >= 1; i--) fprintf(fo, "%ld ", u[-i]);
  fprintf(fo, "\n");
}

/* Save an element to a file fo in GAP list form. */
void Save_GAP_Elm(fo, u)
FILE *fo;
element u;
{
  generator i;

  fprintf(fo, "[");
  for (i = NumGen; i >= 1; i--) {
    fprintf(fo, "%ld", u[-i]);
    if (i != 1) fprintf(fo, ",");
  };
  fprintf(fo, "]\n");
}

/* Print out an element in a readable form. */
void Print_Elm(fo, print_elm_u)
FILE *fo;
element print_elm_u;
{
  generator i;
  exponent k;
  flag zero;

  zero = YES;
  for (i = NumGen; i >= 1; i--) {
    if (print_elm_u[-i] > 0) {
      fprintf(fo, "a%hu", i);
      if (print_elm_u[-i] != 1) fprintf(fo, "^%d ", print_elm_u[-i]);
      else fprintf(fo, " ");
      zero = NO;
    }
    else if (print_elm_u[-i] < 0) {
      fprintf(fo, "A%hu", i);
      if (print_elm_u[-i] != -1) fprintf(fo, "^%d ", -print_elm_u[-i]);
      else fprintf(fo, " ");
      zero = NO;
    }
  };
  if (zero == YES) fprintf(fo, "1");
}

/* Print out an automorphism. */
void Save_Aut(fo, A, i)
FILE *fo;
autstorage A;
generator i;
{
  generator j;
  element u;

  u = (element) A;
  for (j = NumGen; j > i; j--) {
    u += NumGen;
    Save_Elm(fo, u);
  };
  fprintf(fo, "\n");
  for (j = NumGen; j > i; j--) {
    u += NumGen;
    Save_Elm(fo, u);
  };
  fprintf(fo, "\n");
}

/* Prints out an automorphism. */
void Print_Aut(A, i)
autstorage A;
generator i;
{
  Save_Aut(stdout, A, i);
}

