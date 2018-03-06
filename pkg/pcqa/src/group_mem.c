/*
  Written by: Eddie Lo
  Date started: October 24, 1995

  Part of the polycyclic quotient algorithm package.

  This program manages memory for the group module.
*/

#include <stdlib.h>
#include "pcqa.h"

element Identity = NULL;
element id = NULL;
autlist *posaut = NULL;
autlist *negaut = NULL;
exponent *Cpp = NULL;
exponent *Power_Elm = NULL;
exponent *Power_Inv_Elm = NULL;
exponent *finite_index = NULL;
autstorage limited_A = NULL;
exponent *Stack_Elm = NULL;
frlistptr Frame_Chain = NULL;
frlistptr Frame_Elm = NULL;
int stack_ptr = 0;
int frame_size = 0;
int max_stack = 0;
int num_stack = 0;
generator MinComm = 0;
generator NumGen = 0;
generator NumGenP = 1;
counter *table2 = NULL;
counter choose2 = 0;
counter choose3 = 0;

extern flag ExistsPoly;


/* -----  Initialization ----- */

/* Initialize a group with its presentation, used many times. */
void Init_Group()
{
  generator i, j;
  flag commute;
  exponent *ep;

  frame_size = NumGen << Log_EFS;
  stack_ptr = 0; Stack_Elm = NULL; Frame_Elm = NULL; Frame_Chain = NULL;
  max_stack = num_stack = 0;
  Identity = Get_Elm();
  id = Get_Elm();
  for (i = 1; i <= NumGen; i++) {
    Identity[-i] = 0;
    id[-i] = 0;
  };
  choose2 = NumGen*(NumGen-1)/2;
  choose3 = choose2*(NumGen-2)/2;
  MinComm = NumGen;
  commute = 1;
  while (commute && MinComm > 1) {
    ep = (exponent *) *(negaut[MinComm-2]);
    for (i = NumGen; commute && i >= MinComm; i--) {
      ep += NumGen;
      Identity[-i] = 1;
      for (j = MinComm; commute && j <= NumGen; j++)
        if (Identity[-j] != ep[-j]) commute = 0;
      Identity[-i] = 0;
    };
    if (commute) MinComm--;
  };
  if (NumGen > 1) {
    limited_A = (autstorage) calloc(2*(NumGen-1)*NumGen, sizeof(exponent));
    table2 = (counter *) calloc(NumGen-1, sizeof(counter));
    table2[0] = 0;
    for (i = 1; i < NumGen-1; i++) table2[i] = table2[i-1]+NumGen-i;
  };
  Init_Structure();
  ExistsPoly = YES;
}


/* -----  Resetting  ----- */

/* Frees up memory used in the group module. */
void Reset_Group()
{
  generator i;
  counter k;

  Delete_Poly_Presentation();
  if (Frame_Elm != NULL) {
    free(Frame_Elm->list);
    free(Frame_Elm);
  };
  while (Frame_Chain != NULL) {
    Frame_Elm = Frame_Chain;
    Frame_Chain = Frame_Chain->next;
    free(Frame_Elm->list);
    free(Frame_Elm);
  };
  Frame_Elm = NULL;
  Stack_Elm = NULL;
  Reset_Structure();
  if (NumGen > 1) {
    free(table2);
    free(limited_A);
  };
  ExistsPoly = NO;
}


/* -----  Elements  ----- */

/* Get an element from the free element stack. */
element Get_Elm()
{
  if (stack_ptr == 0) {
    if (Frame_Elm != NULL) {
      Frame_Elm->next = Frame_Chain;
      Frame_Chain = Frame_Elm;
      Stack_Elm = Frame_Elm->list;
      Frame_Elm = NULL;
      stack_ptr = Elm_Frame_Size;
    }
    else Allocate_Elm();
    num_stack++;
  };
  stack_ptr--;
  if (num_stack*Elm_Frame_Size-stack_ptr > max_stack)
    max_stack = num_stack*Elm_Frame_Size-stack_ptr;
  return (Stack_Elm+(stack_ptr+1)*NumGen);
}

/* Allocate more elements to be used. This procedure assumes that
   both Frame_Elm and Stack_Elm are NULL. */
void Allocate_Elm()
{
  frlistptr lp;

  lp = Frame_Chain;
  Frame_Chain = (frlistptr) malloc(sizeof(frlist));
  Stack_Elm = (exponent *) calloc(frame_size, sizeof(exponent));
  Frame_Chain->list = Stack_Elm;
  Frame_Chain->next = lp;
  stack_ptr = Elm_Frame_Size;
}

/* Free n elements to the free element stack,
   assuming that n <= Elm_Frame_Size. */
void Free_Elm(n)
counter n;
{
  exponent *current;
  frlistptr lp;

  stack_ptr += n;
  while (stack_ptr >= Elm_Frame_Size) {
    if (Frame_Elm != NULL) {
      free(Frame_Elm->list);
      free(Frame_Elm);
    };
    Frame_Elm = Frame_Chain;
    Frame_Chain = Frame_Chain->next;
    if (Frame_Chain != NULL) Stack_Elm = Frame_Chain->list;
    stack_ptr -= Elm_Frame_Size;
    num_stack--;
  }
}


/* -----  Automorphisms  --- */

/* Get a storage of automorphism from the free element stack,
   assuming that n <= Elm_Frame_Size. */
autstorage Get_Aut(n, disp)
counter n;
int *disp;
{
  if (stack_ptr < n) {
    *disp = stack_ptr+n;
    if (Frame_Elm != NULL) {
      Frame_Elm->next = Frame_Chain;
      Frame_Chain = Frame_Elm;
      Stack_Elm = Frame_Elm->list;
      Frame_Elm = NULL;
      stack_ptr = Elm_Frame_Size;
    }
    else Allocate_Elm();
    num_stack++;
  }
  else *disp = n;
  stack_ptr -= n;
  if (num_stack*Elm_Frame_Size-stack_ptr > max_stack)
    max_stack = num_stack*Elm_Frame_Size-stack_ptr;
  return (autstorage) (Stack_Elm+stack_ptr*NumGen);
}


