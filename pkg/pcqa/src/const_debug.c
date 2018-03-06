/*
  Written by: Eddie Lo
  Date started: February 5, 1996.

  Part of the polycyclic quotient algorithm package.

  This program contains procedures for debugging in the construct module.
*/

#include "pcqa.h"

extern counter *range1, *range2;
extern generator NumMod;


/* This procedure prints out the arrays range1 and range2. */
void Print_ranges()
{
  generator i;

  printf("range1:");
  for (i = 0; i <= NumMod; i++) printf(" %lu", range1[i]);
  printf("\nrange2:");
  for (i = 0; i <= NumMod; i++) printf(" %lu", range2[i]);
  printf("\n");
}

