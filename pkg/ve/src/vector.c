/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */
/* Subroutines for processing vectors and such-like */
#include "me.h"
#include <memory.h>

#ifdef PACKED_VEC
#define PACKCON       ((sizeof(_sent))/sizeof(fldelt)) 
			/* determines when to pack a vector*/

#define PACKERR 10		/* allowed variation */

vector vtidy(v)			/* replace packed vectors w. few non-zero */
				/* entries with sparse ones  etc. */
     vector v;
{
  sent s,send;
  basiselt i,last;
  fldelt *r;
  vector vt;
  if (Type(v) == sparse) /* no point in using PS format here */
    {
      if (1 < Wt(v) && Len(v) + PACKERR < PACKCON * Wt(v)) 
				/* wt 1 is always sparse */
	{			/* packing a sparse vector */
	  newpvec(vt,Len(v));
	  r = BodyP(vt);
	  s = BodyS(v);
	  send = BodyEndS(v);
	  for (i=0; i < Len(v); i++)
	    if (s < send && i == s->loc)
	      {
		r[i] = s->fac;
		s++;
	      }
	    else
	      r[i] = FZero;
	  Wt(vt) = Wt(v);
	  vfree(v);
	}
      else
	vt = v;
    }
  else
    {
      r = BodyP(v);
      last = Len(v);
      while (last && IsZero(r[last-1]))
	last--;
      if (Wt(v) == 1 || Wt(v)*PACKCON+PACKERR < last) 
				/* unpacking a packed vector */
	{
	  newsvec(vt,Wt(v));
	  s = BodyS(vt);
	  for (i=0; i<last; i++)
	    if (!IsZero(r[i]))
	      {
		s->loc = i;
		s->fac = r[i];
		s++;
	      }
	  ASSERT((s-BodyS(vt)) == Wt(v));
	  Len(vt) = last;
	  Wt(vt)  = Wt(v);
	  vfree(v);
	}
      else
	{
	  Len(v) = last;	/* fix up 'short' packed vectors
				 this is just so that they print nicely */
	  vt = v;
	}
    }
  return(vt);
}

#endif

#ifdef GFP /* this is speed critical, so copy P into a register for speed */
#undef myP
#endif

vector adds(v1,v2,f)		/* Compute v1+f*v2 */
     vector v1;
     vector v2;
     fldelt f;
{
#ifdef GFP
  register fldelt myP = P;
#endif
  vector v;
  weight w = 0;
  basiselt last = (basiselt)-1;
  fldelt f1,f2;		/* convenient to do v1.f1 + v2.f2 */
  f1 = FOne; f2 = f;
  if (IsZero(f))			/* might get this I suppose */
    return(vcopy(v1));
#ifdef PACKED_VEC
  if (Type(v1) == packed && Type(v2) == packed)
    {
      basiselt i,l1,l2;
      fldelt *p1,*p2,*p,t,t1;
      bool o1,o2;
      l1 = Len(v1);
      l2 = Len(v2);
      if (l1 < l2)		/* get the longest one first */
	{
	  basiselt l;
	  fldelt ft;
	  l = l1; l1 = l2; l2 = l;
	  v = v1; v1 = v2; v2 = v;
	  ft = f1; f1 = f2; f2 = ft;
	}
      o1 = IsOne(f1);
      o2 = IsOne(f2);
      ASSERT(o1 || o2);
      newpvec(v,l1);		/* allocate space for sum */
      p1 = BodyP(v1);
      p2 = BodyP(v2);
      p = BodyP(v);
      if (myP == 2) /* this is arithmetic-specific */
	{
	  for (i = 0; i< l2; i++)	/* add entries */
	    if ((*p++ = *p1++ ^ *p2++))
	      {
		w++;
		last = i;
	      }
	  for (; i< l1; i++)	/* copy entries */
	    if ((*p++ = *p1++))
	      {
		w++;
		last = i;
	      }
	  Wt(v) = w;
	}
      else
	{
	  if (o1 && o2)
	    {
	      for (i = 0; i< l2; i++)	/* add entries */
		{
		  if (IsZero(*p1))
		    {
		      Fcopy(*p2++,*p);
		      p1++;
		    }
		  else if IsZero(*p2)
		    {
		      Fcopy(*p1++,*p);
		      p2++;
		    }
		  else
		    {
		      Fadd(*p1++,*p2++,*p);
		    }
		  if (!IsZero(*p++))
		    {
		      w++;
		      last = i;
		    }
		}
	      for (; i< l1; i++)	/* copy entries */
		{
		  Fcopy(*p1++,*p);
		  if (!IsZero(*p++))
		    {
		      w++;
		      last = i;
		    }
		}
	    }
	  else
	    if (o1)
	      {
		for (i = 0; i< l2; i++)	/* add entries */
		  {
#ifdef GFP			/* avoid reducing mod P more than we need to */
		    *p = (fldelt)(((xfld)*p1++ +(xfld)*p2++*(xfld)f2)%(xfld)myP);
#else
		    Fmul(f2,*p2++,t1);
		    Fadd(*p1++,t1,*p);
#endif
		    if (!IsZero(*p++))
		      {
			w++;
			last = i;
		      }
		  }
		for (; i< l1; i++)	/* copy entries */
		  {
		    Fcopy(*p1++,*p);
		    if (!IsZero(*p++))
		      {
			w++;
			last = i;
		      }
		  }
	      }
	    else
	      {
		ASSERT(o2);
		for (i = 0; i< l2; i++)	/* add entries */
		  {
#ifdef GFP			/* avoid reducing mod P more than we need to */
		    *p = (fldelt)(((xfld)*p1++*(xfld)f1 +(xfld)*p2++)%(xfld)myP);
#else
		    Fmul(f1,*p1++,t)
		      Fadd(t,*p2++,*p);
#endif
		    if (!IsZero(*p++))
		      {
			w++;
			last = i;
		      }
		  }
		for (; i< l1; i++)	/* copy entries */
		  {
		    Fmul(*p1++,f1,*p);
		    if (!IsZero(*p++))
		      {
			w++;
			last = i;
		      }
		  }
	      }
	}
      Wt(v) = w;
      Len(v) = last+1;
      return(v);		/* return */
    }
  else if (Type(v1) == sparse && Type(v2) == sparse)
#endif
    {
      sent s1,s2,send1,send2;
      basiselt b1,b2;
      sent s;
      bool o;
      ASSERT(IsOne(f1));
      s1 = BodyS(v1);
      send1 = BodyEndS(v1);
      s2 = BodyS(v2);
      send2 = BodyEndS(v2);
      b1 = s1 < send1 ? s1->loc : (basiselt) -1;
      b2 = s2 < send2 ? s2->loc : (basiselt) -1;
      o = IsOne(f2);
      newsvec(v,Wt(v1)+Wt(v2));
      s = BodyS(v);
      while (s1 < send1 || s2 < send2 )	/* essentially a merge operation */
	{
	  if (b1 < b2)
	    {
	      s->loc = b1;
	      Fcopy(s1->fac,s->fac);
	      s1++;
	      b1 = s1 < send1 ? s1->loc : (basiselt) -1;
	      last = s->loc;
	      s++;
	      w++;
	    }
	  else if (b2 < b1)
	    {
	      s->loc = b2;
	      if (o)
		Fcopy(s2->fac,s->fac);
	      else
		Fmul(s2->fac,f2,s->fac);
	      s2++;
	      b2 = s2 < send2 ? s2->loc : (basiselt) -1;
	      ASSERT(!IsZero(s->fac))
	      last = s->loc;
	      s++;
	      w++;
	    }
	  else  /* b1 == b2 */
	    {
	      s->loc = b1;
#ifdef GFP
	      if (P==2)
		s->fac = s1->fac ^ s2->fac;
	      else
#endif
		if (o)
		  Fadd(s1->fac,s2->fac,s->fac);
		else
		  {
		    Fmul(s2->fac,f2,s->fac);
		    Fadd(s1->fac,s->fac,s->fac);
		  }
	      s1++;
	      b1 = s1 < send1 ? s1->loc : (basiselt) -1;
	      s2++;
	      b2 = s2 < send2 ? s2->loc : (basiselt) -1;
	      if (!IsZero(s->fac))	/* don't  write 0 entries */
		{
		  last = s->loc;
		  s++;
		  w++;
		}
	    }
	}
      Wt(v) = w;
      Len(v)=last+1;
      for (b1=w; b1 < Wt(v1)+Wt(v2); b1++)
	{
	  Fclear(s->fac);
	  s++;
	}
      return(v);
    }
#ifdef PACKED_VEC
  else				/* one packed one sparse */
    {
      basiselt l;
      if (Type(v1) == packed)	/* get v1 sparse, v2 packed */
	{
	  fldelt ft;
	  v = v1; v1 = v2; v2 = v;
	  ft = f1; f1 = f2; f2 = ft;
	}
      l = max(Len(v1),Len(v2));
      if ( (Wt(v1)+Wt(v2)) == 1 ||
	  (Wt(v1) + Wt(v2)) * PACKCON >= l) /* sparse + packed = packed */
	{
	  bool o1,o2;
	  basiselt l2 = Len(v2);
	  sent s,send;
	  fldelt *r,*ri,t;
	  basiselt i;
	  o1 = IsOne(f1);
	  o2 = IsOne(f2);
	  ASSERT(o1 || o2);
	  Finit(t);
	  newpvec(v,l);
	  r = BodyP(v);
	  ri = BodyP(v2);
	  if (o2)
	    for (i=0; i< l2; i++)
	      Fcopy(*ri++,*r++);
	  else
	    for (i=0; i< l2; i++)
	      Fmul(*ri++,f2,*r++);
	  for (;i<l;i++)
	    Fcopy(FZero,*r++);
	  r = BodyP(v);
	  w = Wt(v2);
	  send = BodyEndS(v1);
	  if (P==2)
	    for (s = BodyS(v1); s < send; s++)
	      {
		if (*(ri= r+s->loc)) w--;
		if ((*ri ^= s->fac)) w++;
	      }
	  else
	    if (o1)
	      {
		for (s = BodyS(v1); s < send; s++)
		  {
		    if (!IsZero(*(ri= r+s->loc))) w--;
		    Fadd(*ri,s->fac,*ri);
		    if (!IsZero(*ri)) w++;
		  }
	      }
	    else
	      {
		for (s = BodyS(v1); s < send; s++)
		  {
		    if (!IsZero(*(ri= r+s->loc))) w--;
#ifdef GFP
		    *ri = (fldelt)(((xfld)*ri + (xfld)(s->fac)* (xfld)f1)%(xfld)myP);
#else
		    Fmul(f1,s->fac,t);
		    Fadd(*ri,t,*ri);
#endif
		    if (!IsZero(*ri)) w++;
		  }
	      }

	  Wt(v) = w;
	  while (Len(v) && IsZero(r[Len(v)-1]))
	    Len(v)--;
	  Fclear(t);
	  return(v);
	}
      else			/* packed + sparse making sparse */
	{
	  basiselt i,b,l2=Len(v2);
	  fldelt *r;
	  sent s,s1,send1;
	  fldelt fac,t;
	  bool o1,o2;
	  Finit(fac);
	  Finit(t);
	  o1 = IsOne(f1);
	  o2 = IsOne(f2);
	  ASSERT(o1 || o2);
	  newsvec(v,Wt(v1)+Wt(v2));
	  r = BodyP(v2);
	  s1 = BodyS(v1);
	  send1 = BodyEndS(v1);
	  s = BodyS(v);
	  b = s1 < send1 ? s1->loc : (basiselt)-1;
	  if (o1)
	    {
	      for (i=0; i < l2; i++)
		{
		  Fcopy(*r++,fac);
		  if (IsZero(fac))
		    {
		      if (i == b)
			{
			  Fcopy(s1->fac,fac);
			  s1++;
			  b = s1 < send1 ? s1->loc : (basiselt) -1;
			}
		    }
		  else
		    {
		      if(!o2)
			Fmul(f2,fac,fac);
		      if (i == b)
			{
			  if (P == 2)
			    fac ^= s1->fac;
			  else
			    Fadd(fac,s1->fac,fac);
			  s1++;
			  b = s1 < send1 ? s1->loc : (basiselt) -1;
			}
		    }
		  if (!IsZero(fac))
		    {
		      s->loc = i;
		      s->fac = fac;
		      s++;
		      w++;
		      last = i;
		    }
		}
	    }
	  else
	    {
	      ASSERT(o2);
	      for (i=0; i < l2; i++)
		{
		  Fcopy(*r++,fac);
		  if (IsZero(fac))
		    {
		      if (i == b)
			{
			  Fmul(s1->fac,f1,fac);
			  s1++;
			  b = s1 < send1 ? s1->loc : (basiselt) -1;
			}
		    }
		  else
		    {
		      if (i == b)
			{
			  if (P == 2)
			    fac ^= s1->fac;
			  else
			    {
#ifdef GFP
			      fac = (fldelt)(((xfld)(s1->fac) * (xfld)f1 + (xfld)fac)%(xfld)myP);
#else
			      Fmul(s1->fac,f1,t);
			      Fadd(fac,t,fac);
#endif
			    }
			  s1++;
			  b = s1 < send1 ? s1->loc : (basiselt) -1;
			}
		    }
		  if (!IsZero(fac))
		    {
		      s->loc = i;
		      s->fac = fac;
		      s++;
		      w++;
		      last = i;
		    }
		}
	      
	    }
	  while(s1 < send1)
	    {
	      last  = s->loc = s1->loc;
	      if (o1)
		Fcopy(s1->fac,s->fac);
	      else
		Fmul(s1->fac,f1,s->fac);
	      w++;
	      s++;
	      s1++;
	    }
	  Wt(v) = w;
	  Len(v) = last+1;
	  Fclear(fac);
	  Fclear(t);
	  return(v);
	}
    }
#endif
}

#ifdef GFP
#define myP P
#endif

#if defined(RATIONAL) || defined(INTEGRAL)
#define CLEANENTS      sent s; sent send; send = BodyEndS(v); \
                       for (s = BodyS(v); s < send; s++)    Fclear(s->fac)
#else
#define CLEANENTS
#endif


void vfree(v)			/* free a vector */
     vector v;
{
  PS(v,
     myfree((pointer)BodyP(v));
     IFDEBUG(BodyP(v) = NULL);
     ,
     CLEANENTS;
     myfree((pointer)BodyS(v));
     IFDEBUG(BodyS(v) = NULL);
     );
  killvec(v);
}

#undef CLEANENTS

vector smul(v,f)		/* multiply a vector (in place)
				 by a scalar, returns 1st arg */
     vector v;
     fldelt f;
{
  if (!IsOne(f))
    {
      PS(v,
	  fldelt *r;
	  basiselt i; basiselt l;
	  r = BodyP(v);
	  l = Len(v);
	  for (i=0; i<l; i++)
	    {
	      Fmul(*r,f,*r);
	      r++;
	    }
	  if (IsZero(f)) Wt(v) = 0;,
	  sent s; sent send;
	  if (!IsZero(f))
	    {
	      send = BodyEndS(v);
	      for (s = BodyS(v); s < send; s++)
		Fmul(s->fac,f,s->fac);
	    }
	  else
	    {
	      Len(v) = 0;
	      Wt(v) = 0;
	    }
	 );
    }
  return(v);
}
	   

vector btov(b)			/* weight one vector corr. to basis elt b */
     basiselt b;
{
  vector v;
  sent s;
  newsvec(v,1);
  s = BodyS(v);
  Wt(v) = 1;
  Len(v) = b+1;
  s->loc = b;
  Fcopy(FOne,s->fac);
  return(v);
}

vector vcopy(v)			/* copy a vector */
     vector v;
{
  vector vc;
  PS(v,
      basiselt i; basiselt l;
      fldelt *r; fldelt *rc;
      l = Len(v);
      newpvec(vc,l);
      r = BodyP(v);
      rc= BodyP(vc);
      for (i=0; i< l; i++)
      {
        Fcopy(*r,*rc);
	r++;
	rc++;
      },
      sent s; sent sc; sent send;
      newsvec(vc,Wt(v));
      send = BodyEndS(v);
      sc = BodyS(vc);
      for (s = BodyS(v); s<send; s++)
	{
	  sc->loc = s->loc;
	  Fcopy(s->fac,sc->fac);
	  sc++;
	});
  Len(vc) = Len(v);
  Wt(vc) = Wt(v);
  return(vc);
}

#ifdef GFP

static fldelt *fldinvs = NULL;	/* field inverses look up table
				 initialized in fldinv */


fldelt fldinv(f)		/* invert a field elt. */
     fldelt f;
{
  ASSERT(!IsZero(f));
  if (!fldinvs) /* implementation specific stuff */
    {
      fldelt f1,f2;
      fldinvs = (fldelt *) myalloc((P-1)*sizeof(fldelt),false);
      for (f1 =0; f1 < P-1; f1++)	
	fldinvs[f1] = 0;
      for (f1 =1; f1 < P; f1++)	
	if (!fldinvs[f1-1])
	  for (f2 = f1; f2 < P; f2++)
	    {
	      if ((fldelt)((f2*(xfld)f1)%P) == 1)
		{
		  fldinvs[f1-1] = f2;
		  if (f1 != f2)
		    fldinvs[f2-1] = f1;
		  break;
		}
	    }
    }
  return(fldinvs[f-1]);
}



void ClearFldInv()
{
  myfree((pointer)fldinvs);
  fldinvs = NULL;
}

#endif

#ifdef RATIONAL

fldelt fldinv(f)
     fldelt f;
{
  fldelt g SAFE;
  Finit(g);
  mpq_inv(&g,&f);
  return(g);
}

void ClearFldInv()
{
  return;
}

#endif

sent FindEnt(v,b)
     vector v;
     basiselt b;
{
  long top,bottom,middle;
  basiselt b1;
  sent s;
#ifdef PACKED_VEC
  ASSERT(v->type == sparse);
#endif
  if (b > Len(v))
    return NULL;
  top = Wt(v);
  s = BodyS(v);
  bottom = 0;
  while (top >= bottom)
    {
      middle = (top + bottom)/2;
      b1 = s[middle].loc;
      if (b1 < b)
	bottom = middle+1;
      else if (b1 > b)
	top = middle -1;
      else
	return s+middle;
    }
  return NULL;
}


void DelSent(v,s)
     vector v;
     sent s;
{
#ifdef PACKED_VEC
  ASSERT(Type(v) == sparse);
#endif
  ASSERT(BodyS(v) <= s && s < BodyEndS(v));
  Fclear(s->fac);
  if (s-BodyS(v) < --Wt(v))
    memmove((pointer)(s),(pointer)(s+1),(Wt(v) - (s-BodyS(v)))*sizeof(_sent));
  else
    if (s > BodyS(v))
      Len(v) = 1+(s-1)->loc;
    else
      Len(v) = 0;
}

void DelLast(v)
     vector v;
{
#ifdef PACKED_VEC
  ASSERT(Type(v) == sparse);
#endif
  ASSERT(Wt(v));
  Fclear((BodyEndS(v)-1)->fac);
  if (--Wt(v))
    Len(v) = (BodyEndS(v)-1)->loc+1;
  else
    Len(v) = 0;
}



#ifndef PACKED_VEC

vector NewsVecFunc(w)
     basiselt w;
{
  vector v;
  sent s;
  basiselt i;
  newvec(v);
  Len(v) = 0;
  Wt(v) = 0;
  s = BodyS(v) = (sent) myalloc(w*sizeof(_sent),false);
  for (i=0; i < w; i++)
    {
      Finit(s->fac);
      s++;
    }
  return v;
}

#endif
