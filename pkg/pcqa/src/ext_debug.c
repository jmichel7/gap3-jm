/*
  Written by: Eddie Lo
  Date started: December 14, 1995.

  Part of the polycyclic quotient algorithm package.

  This module debugs for the extend module.
*/

#include "pcqa.h"

extern Poly **ppmod, **pnmod, **npmod, **nnmod, **powermod, **npowermod;
extern generator *ppstart, *pnstart, *npstart, *nnstart;
extern generator *powerstart, *npowerstart;
extern counter *table2;
extern flag *gg_defined, *n_defined, *d_found, *rel_found;
extern flag *ppp_found, *powerp_found, *ppower_found, *power_found;
extern flag *d_found, *rel_found;
extern generator NumGen, NumMod;
extern exponent *finite_index;
extern counter Size_R;


/* Print the contents in the stored module element individually. */
void Print_Stored(i, j, ind)
generator i, j;
flag ind;
{
  counter k;

  if (ind & 15) {
    k = table2[i-1]+j-i-1;
    if (ind & 1 && ppstart[k] < NumMod) {
      printf("pp %d %d\n", i, j);
      Print_Mod(stdout, ppmod[k], ppstart[k]);
    };
    if (ind & 2 && pnstart[k] < NumMod) {
      printf("pn %d %d\n", i, j);
      Print_Mod(stdout, pnmod[k], pnstart[k]);
    };
    if (ind & 4 && npstart[k] < NumMod) {
      printf("np %d %d\n", i, j);
      Print_Mod(stdout, npmod[k], npstart[k]);
    };
    if (ind & 8 && nnstart[k] < NumMod) {
      printf("nn %d %d\n", i, j);
      Print_Mod(stdout, nnmod[k], nnstart[k]);
    };
  };
  if (ind & 16 && powerstart[i] < NumMod) {
    printf("power %d\n", i);
    Print_Mod(stdout, powermod[i], powerstart[i]);
  };
  if (ind & 32 && npowerstart[i] < NumMod) {
    printf("npower %d\n", i);
    Print_Mod(stdout, npowermod[i], npowerstart[i]);
  }
}

/* Print the contents in the stored module element. */
void Print_All_Stored(ind)
flag ind;
{
  generator i, j;
  counter k;

  k = 0;
  for (i = 1; i <= NumGen; i++) {
    for (j = i+1; j <= NumGen; j++) {
      if (ind & 1 && ppstart[k] < NumMod) {
        printf("pp %d %d\n", i, j);
        Print_Mod(stdout, ppmod[k], ppstart[k]);
      };
      if (ind & 2 && pnstart[k] < NumMod) {
        printf("pn %d %d\n", i, j);
        Print_Mod(stdout, pnmod[k], pnstart[k]);
      };
      if (ind & 4 && npstart[k] < NumMod) {
        printf("np %d %d\n", i, j);
        Print_Mod(stdout, npmod[k], npstart[k]);
      };
      if (ind & 8 && nnstart[k] < NumMod) {
        printf("nn %d %d\n", i, j);
        Print_Mod(stdout, nnmod[k], nnstart[k]);
      };
      k++;
    };
    if (ind & 16 && powerstart[i] < NumMod) {
      printf("power %d\n", i);
      Print_Mod(stdout, powermod[i], powerstart[i]);
    };
    if (ind & 32 && npowerstart[i] < NumMod) {
      printf("npower %d\n", i);
      Print_Mod(stdout, npowermod[i], npowerstart[i]);
    }
  }
}

/* This procedure prints out which module element hasn't been defined. */
counter Count_Defined(print)
flag print;
{
  generator j1, j2;
  counter entry, k, l;

  entry = k = l = 0;
  for (j1 = 1; j1 < NumGen; j1++) {
    for (j2 = j1+1; j2 <= NumGen; j2++) {
      if (!finite_index[j2]) {
        if (!(gg_defined[entry] & 2)) {
          k++;
          if (print) printf("a%d A%d not defined.\n", j1, j2);
        }
        else l++;
      };
      if (!finite_index[j1]) {
        if (!(gg_defined[entry] & 4)) {
          k++;
          if (print) printf("A%d a%d not defined.\n", j1, j2);
        }
        else l++;
      };
      if (!finite_index[j1] && !finite_index[j2]) {
        if (!(gg_defined[entry] & 8)) {
          k++;
          if (print) printf("A%d A%d not defined.\n", j1, j2);
        }
        else l++;
      };
      entry++;
    }
  };
  for (j1 = 1; j1 <= NumGen; j1++) {
    if (finite_index[j1]) {
      if (!n_defined[j1]) {
        k++;
        if (print) printf("A%d not defined.\n", j1);
      }
      else l++;
    }
  };
  if (print) printf("%d defined, %d undefined.\n", l, k);
  return k;
}

/* This procedure prints the status of collecting rules.
   It returns the number of base rules not found. */
counter Count_Found(print)
flag print;
{
  generator i1, i2, i3;
  counter index, k1, k2;

  k1 = k2 = 0;
  for (i1 = 1; i1 < NumGen-1; i1++) {
    for (i2 = i1+1; i2 < NumGen; i2++) {
      for (i3 = i2+1; i3 <= NumGen; i3++) {
        if (ppp_found[index] == 0) {
          if (print) printf("a%d a%d a%d relation not found.\n", i1, i2, i3);
          k1++;
        }
        else k2++;
        index++;
      }
    }
  };
  index = 0;
  for (i1 = 1; i1 < NumGen; i1++) {
    for (i2 = i1+1; i2 <= NumGen; i2++) {
      if (finite_index[i1]) {
        if (powerp_found[index] == 0) {
          if (print) printf("a%d^%d a%d relation not found.\n",
            i1, finite_index[i1], i2);
          k1++;
        }
        else k2++;
      };
      if (finite_index[i2]) {
        if (ppower_found[index] == 0) {
          if (print) printf("a%d a%d^%d relation not found.\n",
            i1, i2, finite_index[i2]);
          k1++;
        }
        else k2++;
      };
      index++;
    }
  };
  for (i1 = 1; i1 <= NumGen; i1++) {
    if (finite_index[i1] > 1) {
      if (power_found[i1] == 0) {
        if (print) printf("a%d^%d relation not found.\n",
          i1, finite_index[i1]+1);
        k1++;
      }
      else k2++;
    }
  };
  for (i1 = 1; i1 <= NumGen; i1++) {
    if (d_found[i1] == 0) {
      if (print) printf("a%d definition relation not found.\n", i1);
      k1++;
    }
    else k2++;
  };
  for (i1 = 0; i1 < Size_R; i1++) {
    if (rel_found[i1] == 0) {
      if (print) printf("Relation %d not found.\n", i1+1);
      k1++;
    }
    else k2++;
  };
  if (print) printf("%d relation found, %d relation not found.\n", k2, k1);
  return k1;
}

