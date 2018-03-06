/*
  Written by: Eddie Lo
  Date started: August 26,1994.

  Part of the polycyclic quotient algorithm program.

  This program manages memory for the polycyclic quotient algorithm.

Introduction:

  The basic building blocks in every polynomial are cells and MP_INTs.
  This procedure allocates and assigns cells and MP_INTs when needed.

Data structure:

  The basic unit for each polynomial is a cell. It consists of four
  fields. If gen = i, pow = p, then the element it represents is
  a_i^p. The field hpt is a pointer to the next cell and should
  be looked as a "+" sign. The field vpt/cpt is a union. It is
  another pointer to the next cell if the field ind = 0. It is
  a pointer to a MP_INT if the field ind = 1. If ind = 0, vpt should
  be looked as a "*".

Data description:

  Free_Cell_L is a list of free cells.
  Free_MP_INT_L is a list of free MP_INT.
*/

#include "pcqa.h"
#include <malloc.h>

#define Cell_Frame_Size 16384
#define MP_INT_Frame_Size 16384

#define vpt_Num(p)\
  ((p)->ind & 1)

extern generator NumGen, NumGenP;

counter used_cell = 0;
counter used_int = 0;
counter Cur_Cell = 0;
counter Cur_MP_INT = 0;
cellptr Free_Cell_L = NULL;
IntL *Free_MP_INT_L = NULL;
IntL *Int_Frame = NULL;
counter Int_count = 0;


/* This procedure initializes for the memory module, used only once. */
void Est_Basic()
{
  Free_Cell_L = NULL;
  Free_MP_INT_L = NULL;
  Int_Frame = NULL;
  Int_count = MP_INT_Frame_Size-1;
  used_cell = used_int = 0;
  Cur_Cell = Cur_MP_INT = 0;
}

/* This procedure allocates more cells to the list of free cells.
   Require global variable: Free_Cell_L, the list of free cells. */
void Allocate_Cell()
{
  cellptr p;
  counter index;

  p = (cellptr) calloc(Cell_Frame_Size, sizeof(cell));
  Cur_Cell++;
  for (index = 0; index < Cell_Frame_Size; index++, p++) {
    p->hpt = Free_Cell_L;
    Free_Cell_L = p;
  }
}

/* This procedure provides a new cell to the calling procedure.
   Require global variable: Free_Cell_L, the list of free cells. */
cellptr New_Cell()
{
  cellptr p;

  if (Free_Cell_L == NULL) Allocate_Cell();
  p = Free_Cell_L;
  Free_Cell_L = p->hpt;
  used_cell++;
  return p;
}

/* This procedure returns the cell pointed to by p to the free cell list.
   Require global variable: Free_Cell_L, the list of free cells. */
void Free_Cell(p)
cellptr p;
{
  if (vpt_Num(p)) Free_MP_INT(p->u.cpt);
  p->hpt = Free_Cell_L;
  Free_Cell_L = p;
  used_cell--;
}

/* This procedure allocates more MP_INTs to the list of free integers.
   Require global variable: Free_MP_INT_L, the list of free cells. */
void Allocate_MP_INT()
{
  coeff p;
  counter index;

  p = (coeff) calloc(MP_INT_Frame_Size, sizeof(MP_INT));
  Cur_MP_INT++;
  for (index = 0; index < MP_INT_Frame_Size; index++, p++) {
    mpz_init(p);
    Insert_INT(p);
  }
}

/* This procedure inserts a MP_INT into the list of free MP_INT.
   Require global variable: Int_count, the count of the free MP_INT,
                            Free_MP_INT_L, the list of free integers,
                            Int_Frame, a free frame of integer. */
void Insert_INT(p)
coeff p;
{
  Int_count++;
  if (Int_count == MP_INT_Frame_Size) {
    Int_count = 0;
    if (Int_Frame == NULL) {
      Int_Frame = (IntL *) malloc(sizeof(IntL));
      Int_Frame->frame = (coeff *) calloc(MP_INT_Frame_Size, sizeof(coeff));
    };
    Int_Frame->next = Free_MP_INT_L;
    Free_MP_INT_L = Int_Frame;
    Int_Frame = NULL;
  };
  (Free_MP_INT_L->frame)[Int_count] = p;
}

/* This procedure takes a MP_INT from the list of free MP_INT.
   Require global variable: Int_count, the count of the free MP_INT,
                            Free_MP_INT_L, the list of free integers,
                            Int_Frame, a free frame of integer. */
coeff New_MP_INT()
{
  used_int++;
  if (Free_MP_INT_L == NULL) Allocate_MP_INT();
  if (Int_count == 0) {
    if (Int_Frame != NULL) {
      free(Int_Frame->frame);
      free(Int_Frame);
    };
    Int_Frame = Free_MP_INT_L;
    Free_MP_INT_L = Free_MP_INT_L->next;
    Int_count = MP_INT_Frame_Size-1;
    return (Int_Frame->frame)[0];
  };
  Int_count--;
  return (Free_MP_INT_L->frame)[Int_count+1];
}

/* This procedure returns the integer pointed to by p to the free integer list.
   It is assumed that the alloc field of p has been cleared.
   Require global variable: Free_MP_INT_L, the list of free integers. */
void Free_MP_INT(p)
coeff p;
{
  used_int--;
  Insert_INT(p);
}

/* This procedure prints out the content of a coefficient pointer. */
void Print_MP_INT(p)
coeff p;
{
  counter k;

  printf("allocate = %d, size = %d\n", p->_mp_alloc, p->_mp_size);
  for (k = 0; k < p->_mp_alloc; k++) printf("%d ", p->_mp_d[k]);
  printf("\n");
}

/* This procedure prints out the memory usage for the cells and the
   MP_INTs. */
void Print_Basic()
{
  printf("used_cell = %d, used_int = %d\n", used_cell, used_int);
}

