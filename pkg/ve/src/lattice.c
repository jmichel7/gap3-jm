/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */

/* Here are, at least some of the main parts of the lattice code */

/* $Header: /gap/CVS/GAP/3.4.4/pkg/ve/src/lattice.c,v 1.1.1.1 1996/12/11 12:39:38 werner Exp $
   $Log: lattice.c,v $
   Revision 1.1.1.1  1996/12/11 12:39:38  werner
   Preparing 3.4.4 for release
   */
#include "ilatt.h"
#define HAVE_QSORT
lrow Lattice=NULL;
lattpos lattsize = 0;
lattpos nextlat = 0;
lattpos *LFoW;
_lrow LLHead = {NULL,0,(lattpos)-1,(lattpos)-1,LOutOfUse};

#ifdef DEBUG
extern void CheckLatt();
#endif
     
#define LMBLOCK 100

lattmen FreeLMs = NULL;

lattmen lmalloc()
{
  lattmen m,m1;
  m = (lattmen) myalloc(sizeof(_lattmen)*LMBLOCK,false);
  if (!AllBlocks)
    AllBlocks = (pointer *)myalloc(sizeof(pointer)*BLOCKBLOCK,false);
  else if (!(NBlocks %BLOCKBLOCK))
    AllBlocks = myrealloc(AllBlocks,(NBlocks+BLOCKBLOCK)*sizeof(pointer), false);
  AllBlocks[NBlocks++] = (pointer)m;
  for (m1 = m+1; m1 < m+LMBLOCK-1; m1++)
    m1->next = m1+1;
  m1->next = NULL;
  FreeLMs = m+1;
  return(m);
}

retcode AddLatt(v)
     vector v;
{
  lattpos l;
  fldelt g,s,t;
  sent se,s1;
  vector v1,v2,v3;
#ifdef DEBUG
  lmcheck();
  CheckLatt();
#endif
  LOG(loglatt > 1,log("Adding to lattice: "); vprint(logfile,v););
  Finit(g);
  Finit(s);
  Finit(t);
  ASSERT(v == vroot(v));
  v1 = v;
  for (l = LLHead.next; l != (lattpos)-1; l = Lattice[l].next)
    {
      se = FindEnt(Lattice[l].v,Lattice[l].piv);
      ASSERT(se);
      s1 = BodyEndS(v1)-1; /* Look at the LAST non-zero in v1 */
      ASSERT(!IsZero(s1->fac));
      if (IsOne(s1->fac) || IsMOne(s1->fac))
	{
	  stackc(v1);
	  LOG(loglatt > 1, log("Reduced to Coincidence\n"));
	  Fclear(g);
	  Fclear(s);
	  Fclear(t);
	  return OK;
	}
      if (s1->loc == se->loc)
	{
	  LOG(loglatt > 2, log ("Stripping "); vprint(logfile,v1); 
	      fprintf(logfile,"%s with row %u, pivot %hu vector ",logprefix,l,Lattice[l].piv);
	      vprint(logfile,Lattice[l].v););
	  mpz_gcdext(&g,&t,&s,&(s1->fac),&(se->fac)); /* NB order is important here
							 to get the = and =- cases correct */
	  if (IsZero(t))
	    {
	      LOG(loglatt > 2, log("Clears exactly\n"));
	      mpz_div(&t,&(s1->fac),&(se->fac)); /* this must work exactly here */
	      mpz_neg(&t,&t);
	      v2 = adds(v1,Lattice[l].v,t);
	      vfree(v1);
	      v1 = v2;
	    }
	  else if (IsZero(s))
	    {
	      LOG(loglatt > 2, log("Replaces\n"););
	      v2 = Lattice[l].v;
	      ChangeLVec(Lattice+l,v1); /* s1->fac divides the pivot so replace the lattice vector
					   by v1. Updating mention lists and gcds */
	      mpz_div(&t,&(se->fac),&(s1->fac)); /* this must work exactly */
	      mpz_neg (&t,&t);
	      v1 = adds(v2,v1,t);
	      vfree(v2);
	    }
	  else /* Now the general case where neither divides the other */
	    {
	      LOG(loglatt > 2, log("General case\n"););
	      v2 = vcopy(v1);
	      smul(v2,t);
	      v3 = adds(v2,Lattice[l].v,s);
	      vfree(v2);
	      v2 = Lattice[l].v;
	      ChangeLVec(Lattice+l,v3);
	      vfree(v2);
	      mpz_div(&t,&(s1->fac),&g);
	      mpz_neg(&t,&t);
	      v2 = adds(v1,v3,t);
	      vfree(v1);
	      v1 = v2;
	    }
	  /* Here we have flattened the relevant entry of v1 */
	  if (zerop(v1)) 
	    break; /* stripped to zero, in the space */
	}
      else if (s1->loc > se->loc)
	break;
    }
  /* Here v1 is what is left of v, which might be nothing
   We must insert before l (at the end if l is -1 */
  if (!zerop(v1))
    InsLatt(v1,l);   /* If there is something put it into the lattice before the netry in position l
		      find it a pivot column and generally settle it in */
  else
    {
      LOG(loglatt > 1, log("Stripped exactly\n"););
      vfree(v1);
    }
  Fclear(g);
  Fclear(s);
  Fclear(t);
#ifdef DEBUG
  lmcheck();
  CheckLatt();
#endif
  return(OK);
}


retcode ChangeLVec(l,v) /* remove the vector l->v from the lattice, but do not free it
			insert v in the same e-place, but at the end of the w-order
			v will have the same pivot column. Since we have v to hand, we 
			can fix up all the gcds at once */
     lrow l;
     vector v;
{
  lrow l1;
  sent s,send;
  lattpos p,lp = l-Lattice;
  lattmen m,*m1;
#ifdef DEBUG
  lmcheck();
#endif
  LOG(loglatt > 1, log("Replacing vector ");vprint(logfile,l->v);
      fprintf(logfile,"%s in position %d by ",logprefix,l-Lattice);vprint(logfile,v););
  /* Lets begin by putting v into the lattice */
  if (nextlat == lattsize)
    {
      p = lattsize;
      lattsize = min(LATT_MAX,lattsize + max(LATT_BLOCK,lattsize/4));
      if (lattsize <= p)
	return(OutOfSpace);
      Lattice = (lrow)myrealloc((pointer)Lattice,lattsize*sizeof(_lrow),false);
      l = Lattice+lp;
    }
  l1 = Lattice+nextlat;
  l1->v = v;
  l1->piv = l->piv;
  l1->next = l->next;
  l1->prev = l->prev;
  if (l1->next != (lattpos)-1)
    Lattice[l1->next].prev = nextlat;
  else
    LLHead.prev = nextlat;
  if (l1->prev != (lattpos)-1)
    Lattice[l1->prev].next = nextlat;
  else
    LLHead.next = nextlat;
  nextlat++;
  l->status = LOutOfUse;
  l1->status = LInUse;
  /* Now we haveto remove the mentions of l->v and add those for v */
  /* we can sort out gcd's at the same time */
  send = BodyEndS(l->v);
  for (s = BodyS(l->v); s < send; s++)
    {
      m1 = &(table[s->loc].latt.mentions);
      for  (m = *m1; m->pos != lp ; m = *m1)
	{
	  ASSERT(m);
	  m1 = &(m->next);
	}
      *m1 = m->next;
      killlm(m);
      if (closed)
	mpz_set_ui(&(table[s->loc].latt.g),0);
    }
  send = BodyEndS(v);
  for (s = BodyS(v); s < send; s++)
    {
      newlm(m);
      m->next  = table[s->loc].latt.mentions;
      m->pos = nextlat-1;
      table[s->loc].latt.mentions = m;
      if (closed)
	{
	  if (!IsZero(table[s->loc].latt.g))
	    {
	      mpz_gcd(&(table[s->loc].latt.g),&(s->fac),&(table[s->loc].latt.g));
	      if (IsOne(table[s->loc].latt.g))
		StackGcdOne(s->loc);
	    }
	  else
	    FixGcd(s->loc);
	}
    }
#ifdef DEBUG
  lmcheck();
#endif
  return(OK);
}

  
void StackGcdOne(b)
     basiselt b;

     /* process the fact that a column of the lattice has GCD one */
{
  fldelt g,s,t;
  lattmen m;
  vector v,v1;
  sent se;
  Finit(g);
  Finit(s);
  Finit(t);
  m = table[b].latt.mentions;
  ASSERT(m);
  v1 = Lattice[m->pos].v;
  se = FindEnt(v1,b);
  ASSERT(se);
  mpz_abs(&g,&(se->fac));
  v = vcopy(v1);
  if (IsNegative(se->fac))
    smul(v,MinusOne);
  while ((m = m->next))
    {
      v1 = Lattice[m->pos].v;
      se = FindEnt(v1,b);
      mpz_gcdext(&g,&s,&t,&g,&(se->fac));
      if (IsZero(s))
	{
	  ASSERT(!IsZero(t));
	  vfree(v);
	  v = vcopy(v1);
	  if (!IsOne(t))
	    smul(v,t);
	}
      else
	{
	  if (!IsOne(s))
	    {
	      smul(v,s);
	      if (!IsZero(t))
		{
		  v1 = adds(v,v1,t);
		  vfree(v);
		  v = v1;
		}
	    }
	  else
	    {
	      ASSERT(IsZero(t));
	    }
	}
    }
  stackc(v);
  Fclear(g);
  Fclear(s);
  Fclear(t);
  return;
}
	      
void FixGcd(b)
     basiselt b;
     /* Called when the gcd of column b is unknown and needs to be known */
{
  fldelt *gp;
  lattmen m;
  vector v;
  sent s;
  gp = &(table[b].latt.g);
  m = table[b].latt.mentions;
  ASSERT(m);
  v = Lattice[m->pos].v;
  s = FindEnt(v,b);
  mpz_abs(gp,&(s->fac));
  while ((m = m->next))
    {
      v =  Lattice[m->pos].v;
      s = FindEnt(v,b);
      mpz_gcd(gp,gp,&(s->fac));
    }
  if (IsOne(*gp))
    StackGcdOne(b);
}


retcode InsLatt(v,before)	      
     vector v;
     lattpos before;
{
  fldelt a,b;
  sent s,send;
  lattpos p;
  lattmen m;
  lrow l;
  basiselt be;
  basiselt best;
  fldelt fac;
#ifdef DEBUG
  lmcheck();
#endif
  LOG(loglatt > 1, fprintf(logfile,"%sAdding new lattice vector in position %u ",logprefix,nextlat);
      vprint(logfile,v););
  ASSERT(!zerop(v));
  if (nextlat == lattsize)
    {
      p = lattsize;
      lattsize = min(LATT_MAX,lattsize + max(LATT_BLOCK,lattsize/4));
      if (lattsize <= p)
	return(OutOfSpace);
      if (Lattice)
	Lattice = (lrow)myrealloc((pointer)Lattice,lattsize*sizeof(_lrow),false);
      else
	Lattice = (lrow)myalloc(lattsize*sizeof(_lrow),false);
    }
  l = Lattice+nextlat;
  l->v = v;
  l->piv = (BodyEndS(v)-1)->loc;
  send = BodyEndS(v);
  l->status = LInUse;
  l->next = before;
  if (before == (lattpos)-1) 
    {
      l->prev = LLHead.prev;
      LLHead.prev = nextlat;
    }
  else
    {
      l->prev = Lattice[before].prev;
      Lattice[before].prev = nextlat;
    }
  if (l->prev != (lattpos)-1)
    Lattice[l->prev].next = nextlat;
  else
    LLHead.next = nextlat;
  nextlat++;
  for (s = BodyS(v); s < send; s++)
    {
      be = s->loc;
      newlm(m);
      m->pos = nextlat-1;
      m->next = table[be].latt.mentions;
      table[be].latt.mentions = m;
      if (closed)
	{
	  fac = s->fac;
	  if (!IsZero(table[be].latt.g))
	    {
	      mpz_gcd(&(table[be].latt.g),&(table[be].latt.g),&fac);
	      if (IsOne(table[be].latt.g))
		StackGcdOne(be);
	    }
	  else
	    FixGcd(be);
	}
      
    }
#ifdef DEBUG
  lmcheck();
#endif
  return(OK);
}
  
    
void RenameLatt(b1,b2)
     basiselt b1,b2;
     /* called from pack when basiselt b2 is being renamed as b1
	called BEFORE changing the main table 
	pack is responsible for moving the lattice column info, 
	so we just have to fix the actual lattice vectors */
{
  lattmen m;
  vector v;
  sent s;
#ifdef DEBUG
  lmcheck();
#endif
  ASSERT(b1<b2);
  for (m=table[b2].latt.mentions; m; m = m->next)
    {
      v = Lattice[m->pos].v;
      s = FindEnt(v,b2);
      ASSERT(s);
      s->loc = b1;
      if (Lattice[m->pos].piv == b2)
	Lattice[m->pos].piv = b1;
    }
}

typedef struct s_gl {
  gpgen g;
  lattpos l;
} _gl;

typedef _gl *gl;

static int glc(p1,p2)
     pointer p1;
     pointer p2;
{
  lattpos l1 = ((gl)p1)->l,l2=((gl)p2)->l;
  if (l1 < l2)
    return -1;
  else if (l1 == l2)
    return 0;
  else return 1;
}

#ifdef DEBUG
void CheckLatt()
{
  bool *InList;
  int i;
  lattpos p;
  InList = (bool *)myalloc(sizeof(bool)*nextlat,false);
  for (i=0; i < nextlat; i++)
    InList[i] = false;
  for (p = LLHead.next; p != (lattpos)-1;)
    {
      InList[p] = true;
      p = Lattice[p].next;
      if (p != (lattpos)-1 && InList[p])
	DIE();
    }
    for (i=0; i < nextlat; i++)
      if (InList[i] != (Lattice[i].status == LInUse))
	DIE();
  myfree((pointer)InList);
}
	
#endif
void PackLatt()
     /* reclaim space in the lattice table */
{
  lattpos l1,l2=0,last=0;
  int w=0;
  vector v;
  sent s,send;
  lattmen m;
  gl gens,lgens,ogens,olgens;
  gpgen g;
#ifdef DEBUG
  CheckLatt();
#endif
  if (!nextlat)
    return;
  ogens = gens = (gl)myalloc(sizeof(_gl)*(nextgen+1),false);
  for (g = 0; g < nextgen; g++)
    {
      gens[g].g = g;
      gens[g].l = NextToDo[g];
    }
  
  if (LookAhead)
    {
      olgens = lgens = (gl)myalloc(sizeof(_gl)*(nextgen+1),false);
      for (g = 0; g < nextgen; g++)
	{
	  lgens[g].g = g;
	  lgens[g].l = SavedNext[g];
	}
    }
#ifdef HAVE_QSORT
  qsort((pointer)gens,nextgen,sizeof(_gl),glc);
  if (LookAhead)
    qsort((pointer)lgens,nextgen,sizeof(_gl),glc);
#else
  DIE();
#endif
  gens[nextgen].l = (lattpos)-1; /* End stops */
  if (LookAhead)
    lgens[nextgen].l = (lattpos)-1;
  for (l1 = 0; l1 < nextlat; l1++)
    {
      if (Lattice[l1].status == LOutOfUse)
	{
	  for (l2 = max(l2,l1+1); Lattice[l2].status != LInUse && l2 < nextlat; l2++)
	    ;
	  if (l2 < nextlat)
	    {
	      /* Here we have to relabel l2 as l1 */
	      /* First lets do the mention list entries */
	      while(gens->l <= l2)
		{
		  NextToDo[gens->g] = l1;
		  gens++;
		}
	      if (LookAhead)
		while(lgens->l <= l2)
		  {
		    SavedNext[lgens->g] = l1;
		    lgens++;
		  }
	      v = Lattice[l2].v;
	      send = BodyEndS(v);
	      for (s = BodyS(v); s < send; s++)
		{
		  for (m = table[s->loc].latt.mentions; m->pos != l2; m = m->next)
		    ASSERT(m);
		  m->pos = l1;
		}
	      Lattice[l1] = Lattice[l2];
	      if (Lattice[l1].next == (lattpos)-1)
		LLHead.prev = l1;
	      else
		Lattice[Lattice[l1].next].prev = l1;
	      if (Lattice[l1].prev == (lattpos)-1)
		LLHead.next = l1;
	      else
		Lattice[Lattice[l1].prev].next = l1;
	      Lattice[l2].status = LOutOfUse;
	      last = l1;
	      /* now lets update the table that tells us about weights of lattice vectors */
	      while ( w < maxwt && LFoW[w] <= l2 && w < maxwt )
		{
		  LFoW[w] = l1;
		  w++;
		}
#ifdef DEBUG
	      CheckLatt();
#endif
	    }
	  else
	    break;
	}
      else
	last = l1;
    }
  while (gens->l != (lattpos)-1)
    NextToDo[gens++->g] = last+1;
  myfree(ogens);
  if (LookAhead)
    {
      while(lgens->l != (lattpos)-1)
	SavedNext[lgens++->g] = last+1;
      myfree(olgens);
    }
  DumpHeap(0);
  /* Here when we ran out of vectors to rename */
  LOG(loglatt,fprintf(logfile,"%sPacked lattice from %u to %u\n",logprefix,nextlat,last+1));
  nextlat = last+1;
#ifdef DEBUG
  CheckLatt();
#endif
}
		

void CoinLatt(b,v)
     basiselt b;
     vector v;
     /* apply the coincidence b=v to the lattice */
     /* must be called before the b row is cleared in the main table */
{
  lattmen m,m1;
  lattpos p,p1;
  vector v0,v1;
  sent s;
  fldelt fac;
  Finit(fac);

  LOG(loglatt > 2, fprintf(logfile,"%sApplying to lattice : e%hu = ",logprefix,b); vprint(logfile,v));
  for (m= table[b].latt.mentions; m ; m = m1)
    {
      m1 = m->next; /* have to do this as m will probably get deleted as we work */
      p = m->pos;
      LOG(loglatt > 2, fprintf(logfile,"%sAltering lattice entry %u\n",logprefix,p););
      ASSERT(Lattice[p].status == LInUse);
      v0 = vcopy(Lattice[p].v);
      s = FindEnt(v0,b);
      ASSERT(s);
      Fcopy(s->fac,fac);
      DelSent(v0,s);
      v1 = adds(v0,v,fac);
      vfree(v0);
      LOG(loglatt>2, log(""); vprint(logfile,Lattice[p].v); log(" becomes "); vprint(logfile,v1););
      /* Now we want to replace v0 by v1, but this might break the echelonization */
      /* Lets do this cleanly */
      RemLatt(p);
      stackc(v1);
    }
  ASSERT(!table[b].latt.mentions);
#ifdef DEBUG
  lmcheck();
#endif
  Fclear(fac);
}

void RemLatt(p)
     lattpos p;
{
  vector v;
  lrow l;
  sent s,send;
  lattmen m,*m1;
  l = Lattice+p;
  v = l->v;
  if (l->next != (lattpos)-1)
    Lattice[l->next].prev = l->prev;
  else
    LLHead.prev = l->prev;
  if (l->prev != (lattpos)-1)
    Lattice[l->prev].next = l->next;
  else
    LLHead.next = l->next;
  l->status = LOutOfUse;
  send = BodyEndS(v);
  for (s = BodyS(v); s < send; s++)
    {
      m1 = &(table[s->loc].latt.mentions);
      while ((m = *m1)->pos != p)
	{
	  m1 = &(m->next);
	  ASSERT(*m1);
	}
      *m1 = m->next;
      killlm(m);
      if (closed)
	mpz_set_ui(&(table[s->loc].latt.g),0);
    }
  vfree(v);
}

	  
void PrintLatt(f,pref)
     FILE *f;
     char *pref;
{
  lrow l;
  lattpos p;
  for (p = LLHead.next; p != (lattpos)-1; p = Lattice[p].next)
    {
      l = Lattice+p;
      fprintf(f,"%sPosition %d, pivot %hu vector ",pref,p,l->piv);
      vprint(f,l->v);
    }
}

void OutLatt(f,m)
     FILE *f;
     OutMode m;
{
  lattpos p;
  outstyle *s;
  s = outstyles+m;
  if (s->stmat)
    (s->stmat)(f,(gpgen)-2);
  for (p = LLHead.next; p != (lattpos)-1; p = Lattice[p].next)
    {
      if (p != LLHead.next)
	fputs(s->vsep,f);
      OutRow(f,Lattice[p].v,s);
    }
  fputs(s->emat,f);
  fputs(s->eims,f);
}

lattpos *NextToDo;
lattpos *SavedNext;
weight *GenWt,*GenLWt;

retcode PushLatt(w)
     weight w;
     /* This generates all the images needed at weight w */
{
  weight gw;
  gpgen g;
  retcode rc;
  lattpos p;
  vector v1,v2;
  LOG(loglatt,fprintf(logfile,"%sPushing lattice at weight %u\n",logprefix,w));
  for (g =0; g < nextgen; g++)
    if (inverse[g] == NOINVERSE || inverse[g] >= g)
      {
	gw = LookAhead ? GenWt[g] : GenLWt[g];
	while (NextToDo[g] < LFoW[w-gw+1])
	  {
	    p = NextToDo[g]++;
	    if (Lattice[p].status == LInUse)
	      {
		LOG(loglatt > 1, fprintf(logfile,"%sPushing lattice vector %u under %s\n",logprefix,p,gennames[g]));
		v1 = action(Lattice[p].v,g,!LookAhead,&rc);
		switch (rc)
		  {
		  case OK:
		    newsvec(v2,0);
		    stackeq(v1,v2);
		    vfree(v1);
		    vfree(v2);
		    dostacks();
		    break;
		  case NeedToDefine:
		    ASSERT(LookAhead);
		    break;
		  case OutOfSpace:
		    return(rc);
		  default:
		    DIE();
		  }
	      }
	  }
      }
  for (g = 0; g < nextgen; g++)
    if (inverse[g] == NOINVERSE || inverse[g] >= g)
      if (NextToDo[g] < nextlat)
	return(OK);
  return(NoMoreWork);
}

	    

		    
void lmcheck()
{
  basiselt b;
  lattmen m;
  lattpos p;
  sent s,send;
  vector v;
  for (b = 0; b < nextbe; b++)
    if (table[b].status == InUse)
	for (m = table[b].latt.mentions; m; m = m->next)
	  {
	    ASSERT(FindEnt(Lattice[m->pos].v,b));
	  }
  for (p = 0; p < nextlat; p++)
    if (Lattice[p].status == LInUse)
      {
	v = Lattice[p].v;
	send = BodyEndS(v);
	for (s = BodyS(v); s < send; s++)
	  {
	    ASSERT(!rowdel(s->loc));
	    for (m = table[s->loc].latt.mentions; m && m->pos != p; m = m->next)
	      ;
	    ASSERT(m);
	  }
      }
}


void LattClosing()
{
  basiselt b;
  LOG(loglatt,log("Table closing, lattice processing under more liberal rules\n"));
  for (b=0; b < nextbe; b++)
    if (table[b].status == InUse && table[b].latt.mentions)
      FixGcd(b);
  dostacks();
}


void CleanLatt()
{
  lattpos p;
  for (p=0; p < nextlat; p++)
    if (Lattice[p].status == LInUse)
      vfree(Lattice[p].v);
  if (Lattice)
    myfree((pointer)Lattice);
  lattsize = nextlat = 0;
  Lattice = NULL;
}


#ifdef DEBUG

void StampLVecs()
{
  lattpos p;
  for (p=0; p < nextlat; p++)
    if (Lattice[p].status == LInUse)
      Lattice[p].v->dbflags |= DBF_INLATT;
}


#endif
