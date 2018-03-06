/*
  Written by: Eddie Lo
  Date started: May 6, 1996.

  Part of the polycyclic quotient algorithm package.

  This module reads in the positive conjugacy relations and some power
  relations of a polycyclic presentation and completes it to make it
  a consistent polcycylic presentation.
*/

#include <stdlib.h>
#include "pcqa.h"

extern generator NumGen, NumGenP;
extern autlist *posaut, *negaut;
extern exponent *Cpp;
extern exponent *Power_Elm, *Power_Inv_Elm, *finite_index;


/* This procedure reads in an incomplete polycyclic presentation. */
void Read_Incomplete(fi)
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
    disp += 2*(NumGen-i)*NumGen;
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

/* This procedure completes an incomplete polycyclic presentation. */
void Complete_Presentation()
{
  generator i, j, ind;
  element u1, u2;

  for (i = NumGen-1; i >= 1; i--) {
    u1 = (element) *(negaut[i-1])+(NumGen-i)*NumGen;
    u2 = (element) *(negaut[i-1])+2*(NumGen-i)*NumGen;
    for (j = i+1; j <= NumGen; j++) {
      Group_Elm_Inv(u2, u1);
      u1 -= NumGen;
      u2 -= NumGen;
    };
    u1 = (element) *(negaut[i-1])+(NumGen-i)*NumGen;
    u2 = (element) *(posaut[i-1])+2*(NumGen-i)*NumGen;
    for (j = i+1; j <= NumGen; j++) {
      for (ind = 1; ind <= NumGen; ind++) u2[-ind] = u1[-ind];
      u1 -= NumGen;
      u2 -= NumGen;
    };
    Group_Aut_Inv(*(posaut[i-1]), *(posaut[i-1])+(NumGen-i)*NumGen, i);
    u1 = (element) *(posaut[i-1])+(NumGen-i)*NumGen;
    u2 = (element) *(posaut[i-1])+2*(NumGen-i)*NumGen;
    for (j = i+1; j <= NumGen; j++) {
      Group_Elm_Inv(u2, u1);
      u1 -= NumGen;
      u2 -= NumGen;
    }
  }
}

