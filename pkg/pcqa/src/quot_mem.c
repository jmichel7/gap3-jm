/*
  Written by: Eddie Lo
  Date started: February 6, 1996.

  Part of the polycyclic quotient algorithm package.

  This program initalizes memory the quotient module.
*/

#include <stdlib.h>
#include "pcqa.h"

element quotient_u1 = NULL;
element quotient_u2 = NULL;
element quotient_u3 = NULL;
element sift_u = NULL;
element *Arr = NULL;
flag *occ = NULL;

extern generator NumGen;
extern flag ExistsQuotient;


/* Initialize for the quotient module. */
void Init_Quotient()
{
  generator i;

  quotient_u1 = (element) calloc(4*NumGen, sizeof(exponent))+NumGen;
  quotient_u2 = quotient_u1+NumGen;
  quotient_u3 = quotient_u2+NumGen;
  sift_u = quotient_u3+NumGen;
  Arr = (element *) calloc(NumGen, sizeof(element));
  Arr[0] = (element) calloc(NumGen*NumGen, sizeof(exponent))+NumGen;
  for (i = 1; i < NumGen; i++) Arr[i] = Arr[i-1]+NumGen;
  occ = (flag *) calloc(NumGen, sizeof(flag));
  ExistsQuotient = YES;
}

/* Reset memory for the quotient module. */
void Reset_Quotient()
{
  free(quotient_u1-NumGen);
  free(Arr[0]-NumGen);
  free(Arr);
  free(occ);
  ExistsQuotient = NO;
}

