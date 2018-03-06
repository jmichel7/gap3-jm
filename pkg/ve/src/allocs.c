/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */
#include <string.h>
#include "me.h"

pointer *AllBlocks = NULL;
int NBlocks = 0;

vector FreeVecs=NULL;
#define VecBlock 300

#ifdef DEBUG
unsigned int vecsinuse=0, freevecs=0, totalvecs=0;
vector AllVecs = NULL;
#endif


vector valloc()
{
  vector v,v1;
  int i;
  v = (vector)myalloc(sizeof(_vector)*VecBlock,false);
  if (!AllBlocks)
    AllBlocks = (pointer *)myalloc(sizeof(pointer)*BLOCKBLOCK,false);
  else if (!(NBlocks %BLOCKBLOCK))
    AllBlocks =
      (pointer *)myrealloc(AllBlocks,(NBlocks+BLOCKBLOCK)*sizeof(pointer), false);
  AllBlocks[NBlocks++] = (pointer)v;
#ifdef DEBUG
  totalvecs += VecBlock;
  freevecs += VecBlock-1;
  vecsinuse++;
  for (i=0; i< VecBlock-1; i++)
     v[i].dbnext = v+i+1;
  v[VecBlock-1].dbnext = AllVecs;
  AllVecs = v;
  v->dbflags = 0;
  if (inpusha) v->dbflags |= DBF_PUSHA;
  if (inpushg) v->dbflags |= DBF_PUSHG;
  if (doingstacks) v->dbflags |= DBF_DOSTACKS;
  if (insubgens) v->dbflags |= DBF_SUBGENS;
#endif
  FreeVecs = v1 = v+1;
  for (i=0; i<VecBlock-2; i++)
    v1 = NextFreeV(v1) = v1+1;
  NextFreeV(v1) = NULL;
  return(v);
}

#ifndef LEN
menlist FreeMenl=NULL;
#define MenBlock 100

#ifdef DEBUG
unsigned int mensinuse=0, freemens=0, totalmens=0;
#endif


menlist menalloc()
{
  menlist m,m1;
  int i;
  m = (menlist)myalloc(sizeof(_menl)*MenBlock,false);
  if (!AllBlocks)
    AllBlocks = (pointer *)myalloc(sizeof(pointer)*BLOCKBLOCK,false);
  else if (!(NBlocks %BLOCKBLOCK))
    AllBlocks = myrealloc(AllBlocks,(NBlocks+BLOCKBLOCK)*sizeof(pointer), false);
  AllBlocks[NBlocks++] = (pointer)m;
#ifdef DEBUG
  totalmens += MenBlock;
  freemens += MenBlock-1;
  mensinuse++;
#endif
  FreeMenl = m1 = m+1;
  for (i=0; i< MenBlock-2; i++)
    m1 = m1->next = m1+1;
  m1->next = NULL;
  return(m);
}

#endif

void CleanAllocs()
{
  int i;
  for (i=0; i < NBlocks; i++)
    myfree(AllBlocks[i]);
  myfree(AllBlocks);
  AllBlocks = NULL;
  NBlocks = 0;
}


#if defined(RATIONAL) || defined(INTEGRAL)

#define GMP_LIMB mp_limb_t


/* This all works on the assumption that an mp_limb is at least
   as big as a pointer */

#define WBLOCK 300

typedef union _w *wp;

union _w {
  wp p;
  GMP_LIMB w;
};

wp FreeWps = {NULL};

typedef union _dw *dwp;

union _dw {
  dwp p;
  GMP_LIMB w[2];
};

dwp FreeDwps = {NULL};

pointer alloc_for_gmp(n)
     size_t n;
{
  pointer p;
  wp w,w1;
  dwp dw,dw1;
  int i;
  if (n == sizeof(GMP_LIMB))
    {
      if (FreeWps)
	{
	  p = (pointer)FreeWps;
	  FreeWps = FreeWps->p;
	  return p;
	}
      else
	{
	  w1 = w = (wp) myalloc(sizeof(union _w)*WBLOCK,false);
	  if (!AllBlocks)
	    AllBlocks = (pointer *)myalloc(sizeof(pointer)*BLOCKBLOCK,false);
	  else if (!(NBlocks %BLOCKBLOCK))
	    AllBlocks = (pointer *)myrealloc(AllBlocks,(NBlocks+BLOCKBLOCK)*sizeof(pointer), false);
	  AllBlocks[NBlocks++] = (pointer)w++;
	  for (i=0; i < WBLOCK-2; i++)
	    {
	      w->p = w+1;
	      w++;
	    }
	  w->p = (wp)NULL;
	  FreeWps = w1+1;
	  return (pointer)w1;
	}
    }
  else if (n == 2*sizeof(GMP_LIMB))
    {
      if (FreeDwps)
	{
	  p = (pointer)FreeDwps;
	  FreeDwps = FreeDwps->p;
	  return p;
	}
      else
	{
	  dw1 = dw = (dwp) myalloc(sizeof(union _dw)*WBLOCK,false);
	  if (!AllBlocks)
	    AllBlocks = (pointer *)myalloc(sizeof(pointer)*BLOCKBLOCK,false);
	  else if (!(NBlocks %BLOCKBLOCK))
	    AllBlocks = (pointer *)myrealloc(AllBlocks,(NBlocks+BLOCKBLOCK)*sizeof(pointer), false);
	  AllBlocks[NBlocks++] = (pointer)dw++;
	  for (i=0; i < WBLOCK-2; i++)
	    {
	      dw->p = dw+1;
	      dw++;
	    }
	  dw->p = (dwp)NULL;
	  FreeDwps = dw1+1;
	  return (pointer)dw1;
	}
    }
  else
    return myalloc(n,false);
}

void free_for_gmp(p, n)
     pointer p;
     size_t n;
{
  if (n==sizeof(GMP_LIMB))
    {
      ((wp)p)->p = FreeWps;
      FreeWps = (wp)p;
    }
  else if (n == 2*sizeof(GMP_LIMB))
    {
      ((dwp)p)->p = FreeDwps;
      FreeDwps = (dwp)p;
    }
  else
    myfree(p);
}

pointer realloc_for_gmp(p, o, n)
     pointer p;
     size_t o,n;
{
  pointer p1;
  if (o == n)
    return p;
  /* first get the new space */
  p1 = alloc_for_gmp(n);
  /* now copy the data */
  if (o <= n)
    memcpy(p1,p,o);
  else
    memcpy(p1,p,n);
  /* now deallocate the old space */
  free_for_gmp(p,o);
  return p1;
}



#endif
