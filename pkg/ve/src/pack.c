/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */
#include <stdlib.h>
#include <string.h>
#include "me.h"

#define HAVE_QSORT

static int rct=0;		/* number of relators. Initalized first
				   time pack is called*/
static rell *rl = NULL;		/* space for flat list of relators  - likewise*/
static rell *lrl = NULL;		/* and another copy for sorting onm a different field */
				   
void CleanPack()
{
  if (rl)
    myfree((pointer)rl);
  if (lrl)
    myfree((pointer)lrl);
  rl = NULL;
  lrl = NULL;
  rct = 0;
}



static int cmprel(r1,r2)
	const pointer r1,r2;
{
    basiselt n1= (*(rell *)r1)->NexttoDo,n2=(*(rell*)r2)->NexttoDo;
    if (n1 < n2) return -1;
    else if (n1 == n2) return 0;
    else return 1;
 }

static int lcmprel(r1,r2)
	const pointer r1,r2;
{
    basiselt n1= (*(rell*)r1)->SavedNext,n2=(*(rell*)r2)->SavedNext;
    if (n1 < n2) return -1;
    else if (n1 == n2) return 0;
    else return 1;
 }

				   
void pack()				/* pack the table. */
{
#ifdef LEN
  vector v;
  basiselt *newname;
#else
  menlist m;
#ifdef PACKED_VEC
  pl p;
#endif
#endif
  basiselt b,b1,b2=0,last=0;
  rell r;
  int i,pt,lpt;
  weight w;
  gpgen g;
#ifdef LEN
  newname = (basiselt *)myalloc(sizeof(basiselt)*nextbe, false);
  for (b=0; b < nextbe; b++)
    {
      newname[b] = (basiselt)-1;
      if (table[b].status == InUse)
	for (g=0; g < nextgen; g++)
	  if ((v=tabent(b,g)))
	    {
	      v = vroot(v);
	      if (v != tabent(b,g))
		{
		  vfree(tabent(b,g));
		  tabent(b,g) = v;
		}
	    }
    }
#endif
  if (nalive == nextbe)		/* no progress made, nothing to pack */
    {
#ifdef LEN
      myfree((pointer)newname);
#endif
    return;
    }
  LOG(logpack,
    fprintf(logfile,"%sPacking %u to %u\n",logprefix,nextbe,nalive););
#ifdef INTEGRAL
  PackLatt();
#endif
  if (!rct)
    for (r = relators; r; r= r->next)
      rct++;			/* count the relators */
  if (!rl)
    {
      rl = (rell *)myalloc(rct*sizeof(rell),false);
      i = 0;
      for (r = relators; r; r=r->next) /* get the list of relators */
	rl[i++] = r;
    }
  qsort((pointer)rl,(size_t)rct,sizeof(rell),cmprel);
  pt = 0;
  if (LookAhead)		/* In lookahead we need a second list also */
    {
      if (!lrl)
	{
	  lrl = (rell *)myalloc(rct*sizeof(rell),false);
	  memcpy((pointer)lrl,(pointer)rl,rct*sizeof(rell));
	}
      qsort((pointer)lrl,(size_t)rct,sizeof(rell),lcmprel);
      lpt = 0;
    }
  w = 1;
  /* Now the main loop. */			   
  for (b1=0; b1 < nextbe; b1++) 
    {
      if (b1 >= startvecs && table[b1].status != InUse) /* got a deleted row */
	{
	  LOG(logpack>1,
	    fprintf(logfile,"%sClearing %d\n",logprefix,b1););
	  if (rowdel(b1))
	    vfree(replace(b1)); /* discard this pointer */
	                        /* look for undeleted row to copy */
	  for (b2 = max(b1+1,b2);
	       b2 < nextbe && (table[b2].status != InUse);
	       b2++)
	    ;
	  if (b2 == nextbe)	/* no more rows */
	    {
	      table[b1].status = OutOfUse; /* can't re-use this one */
	      break;		/* leave the b1 loop */
	    }

				/* now relocating b2 as b1 */
				/* first correct packed vectors */
	  last = b1;
	  LOG(logpack>1,
	    fprintf(logfile,"%sRelocating %d to %d\n",logprefix,b2,b1););
	  SC(scpk(b1,b2));
#ifdef LEN
	  newname[b2] = b1;
#endif
#ifdef INTEGRAL
	  RenameLatt(b1,b2);
#endif
#ifndef LEN
#ifdef PACKED_VEC
	  for (p = plisthead->next; p != plisthead; p = p->next)
	    {
	      vector vp = elookup(p->en);
	      fldelt *r = BodyP(vp);
	      if (Len(vp) > b2 && !IsZero(r[b2]))
		{
		  r[b1] = r[b2];
		  r[b2] = FZero;
		  if (Len(vp) == b2+1)
		    Len(vp) = b1+1;
		}
	    }
#endif
	  for (m = table[b2].mentions; m; m=m->next)
	    {			/* now the ones in sparse vectors */
	      vector v = elookup(m->en);
	      sent s;
	      for (s = BodyS(v); s->loc != b2; s++)
		;
	      s->loc = b1;
	      if (Len(v) == b2+1)
		Len(v) = b1+1;
	    }
	  /* Now we have to sort out other
	     peoples mention lists, and the pointers
	     up from any packed entries */
	  table[b1].mentions = table[b2].mentions;
	  for (g=0; g < nextgen; g++)
	    {
	      vector v;
	      if ((v = tabent(b2,g)))
		{
		  PS(v,
		    v->list.p->en.be = b1;,
		    {
		      sent s; sent send = BodyEndS(v);
		      for (s = BodyS(v); s<send; s++)
			{
			  for (m = table[s->loc].mentions;
			       m->en.be != b2 || m->en.g != g;
			       m = m->next)
			    ASSERT(m);
			  m->en.be = b1;
			}
		    });

		}
	    }
#endif
	  table[b1].status = InUse;
	  table[b1].v.entries = table[b2].v.entries;
	  table[b1].fmgen = table[b2].fmgen;
	  table[b1].wlen = table[b2].wlen;
	  table[b1].word = table[b2].word;
	  table[b2].status = OutOfUse;
#ifdef INTEGRAL
	  table[b1].latt = table[b2].latt;
#endif
	}
      else			/* b1 undeleted  or b1 is a starting vector*/
	{
	  b2 = b1;
	  last = b1;
#ifdef LEN
	  newname[b1] = b1;
#endif
	}
      while (pt < rct && b2 >= rl[pt]->NexttoDo)
	rl[pt++]->NexttoDo = b1;
      if (LookAhead)
	while (lpt < rct && b2 >= lrl[lpt]->SavedNext)
	  lrl[lpt++]->SavedNext = b1;
      while (FirstofWt[w] != (basiselt)-1 && b2 >= FirstofWt[w])
	FirstofWt[w++] = b1;
    }
  for (b1 = last+1; b1 < nextbe; b1++)
    if (rowdel(b1))
      {
	vfree(replace(b1));
	table[b1].status = OutOfUse;
      }
  nextbe = last+1;
  while (pt < rct)
    rl[pt++]->NexttoDo = nextbe;
  if (LookAhead)
    while (lpt < rct)
      lrl[lpt++]->SavedNext = nextbe;
  while (FirstofWt[w] != (basiselt)-1)
    FirstofWt[w++] = nextbe;
#if defined(PACKED_VEC) || defined(LEN)
  for (b1 = 0; b1 < nextbe; b1++) /* finally, tidy all entries */
    if (table[b1].status == InUse)
      {
        gpgen g;
        for (g=0; g < nextgen; g++)
	  {
	    vector v1;
#ifdef PACKED_VEC
	    vector v2;
#endif
#ifdef LEN
	    basiselt vlast;
#endif
	    if ((v1 = tabent(b1,g)))
	      {
#ifdef LEN
		PS(v1,
		   basiselt l=Len(v1);
		   basiselt i;
		   fldelt *r =BodyP(v1);
		   fldelt *r1 = r;
		   vlast = 0;
		   for (i=0; i < l; i++,r++)
		      if (!IsZero(*r))
		         {
			   ASSERT(newname[i] != (basiselt)-1);
			   if (newname[i] != i)
			     {
			       r1[newname[i]] = *r;
			       Fcopy(FZero,*r);
			     }
			   vlast = newname[i]+1;
			 },
		   sent s;
		   sent send = BodyEndS(v1);
		   vlast = 0;
		   for (s = BodyS(v1); s < send; s++)
		     {
		       ASSERT(newname[s->loc] != (basiselt)-1);
		       s->loc = newname[s->loc];
		       vlast = s->loc+1;
		     });
		Len(v1) = vlast;
#endif
#ifdef PACKED_VEC
		v1 = vcopy(v1);
		v2 = vtidy(v1);
		if (v1 != v2)
		  {
		    clrent(b1,g);
		    install(b1,g,v2);
		  }
		else 
		  {
		    vfree(v2);
		  }
#endif
	      }
	  }
      }
#ifdef LEN
  myfree((pointer)newname);
#endif
#endif
}




