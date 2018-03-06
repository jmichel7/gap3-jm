/* This program is Copyright (C) Stephen Linton 1991 and 1992  
  See the file copying for details */

/* This is the internal header of the lattice code. */

/* $Header: /gap/CVS/GAP/3.4.4/pkg/ve/src/ilatt.h,v 1.1.1.1 1996/12/11 12:39:38 werner Exp $
   $Log: ilatt.h,v $
   Revision 1.1.1.1  1996/12/11 12:39:38  werner
   Preparing 3.4.4 for release
   */

#include "me.h"		/* incidentally picks up the external */
				/* header for the lattice stuff */

/* Internal functions - these are global names, but not published as it were */

extern lattmen FreeLMs;


#define newlm(l)        {if (((l)=FreeLMs)) FreeLMs = FreeLMs->next; else (l) = lmalloc();}
#define killlm(l)       {(l)->next = FreeLMs; FreeLMs = (l);}

extern lattmen lmalloc PT((void));

typedef enum {LInUse, LOutOfUse} lstatus;

typedef struct s_lrow *lrow;

typedef struct s_lrow {
  vector v;			/* the lattice vector */
  basiselt piv;			/* the pivot column */
  lattpos next;			/* the ordering wrt which the v are echelonized */
  lattpos prev;			/* is given by this DLL */
  lstatus status;		/*  */
} _lrow;

extern _lrow LLHead;

extern lrow Lattice;
extern lattpos lattsize;


extern retcode ChangeLVec();
extern retcode InsLatt();
extern void StackGcdOne();
extern void FixGcd();
extern void RemLatt();
extern void lmcheck();

#define LATT_BLOCK 1000
#define LATT_MAX 100000






