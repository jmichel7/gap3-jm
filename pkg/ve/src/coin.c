/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */
/* Code to deal with coincidences */
#include "me.h"

coinl coincidences=NULL;	/* stacks */
deductl deductions=NULL;

void stackc(v)			/* Stack a coincidence */
     vector v;
{
  coinl c;
  LOG(logcoin,
      fprintf(logfile,"%sStacked coincidence 0 = ",logprefix);
      vprint(logfile,v););
  SC(sccheckc(v));
  c = (coinl)myalloc(sizeof(_coinl),false);
  c->body.v = v;
  c->next = coincidences;
  coincidences = c;
}

void stackd(b,g,v)			/* Stack a deduction */
     basiselt b;
     gpgen g;
     vector v;
{
  deductl d;
  LOG(logcoin,
      fprintf(logfile,"%sStacked deduction e%d.%s = ",logprefix,b,gennames[g]);
      vprint(logfile,v););
  d = (deductl)myalloc(sizeof(_deductl),false);
  d->body.b = b;
  d->body.g = g;
  d->body.v = v;
  d->next = deductions;
  deductions = d;
}

retcode stackxeq(v1,x,v2)	/* handle the info. v1.x=v2 */
     vector v1,v2;
     gpgen x;
{
  vector v1x;
#ifndef INTEGRAL
  vector v;
  basiselt b;
  fldelt f SAFE,fi SAFE;
  bool def =false;
  vector v0;
#endif
  retcode rc;
  SC(sccheckxeq(v1,x,v2));
  ASSERT(v1 == vroot(v1));
  ASSERT(v2 == vroot(v2));
#ifdef INTEGRAL /* then too hard to be clever  for the time being*/
  v1x = action(v1,x,DefinesCritical,&rc);
  if (!v1x)
    return rc;
  else
    {
       stackeq(v1x,v2);  
       vfree(v1x);
     }
#else 
  newsvec(v1x,0);
  PS(v1,
     {
       fldelt *r=BodyP(v1);
       basiselt i;
       basiselt l = Len(v1);
       for (i=0; i<l; i++,r++)
	 if (!IsZero(*r))
	   {
	     v0 = baction(i,x,def ? DefinesCritical : NoDefines,&rc);
	     if (v0)
	       {
		 v = adds(v1x,v0,*r);
/* 		 if (logcoin > 4)
		   {
		     log("+++");
		     vprint(logfile,v1x);
		     vprint(logfile,v0);
		     outfelt(logfile,*r);
		     fputc('\n',logfile);
		     vprint(logfile,v);
		     log("###");
		   } */
		 vfree(v1x);
		 v1x = v;
	       }
	     else if (!def && rc == NeedToDefine)
	       /* here on the first undefined only */
	       {
		 b = i;
		 f = *r;
		 def = true;
	       }
	     else
	       {
		 vfree(v1x);
		 return(rc);
	       }
	   }
     },
     {
       sent s;
       sent send;
       send = BodyEndS(v1);
       for (s=BodyS(v1); s < send ; s++)
	 {
	   v0 = baction(s->loc,x,def ? DefinesCritical : NoDefines,&rc);
	   if (v0)
	     {
	       v = adds(v1x,v0,s->fac);
/*	       if (logcoin > 4)
		 {
		   log("***");
		   vprint(logfile,v1x);
		   vprint(logfile,v0);
		   outfelt(logfile,s->fac);
		   fputc('\n',logfile);
		   vprint(logfile,v);
		   log("---");
		 } */
	       vfree(v1x);
	       v1x = v;
	     }
	   else if (!def && rc == NeedToDefine)
	     /* here on the first undefined only */
	     {
	       b = s->loc;
	       f = s->fac;
	       def = true;
	     }
	   else 
	     {
	       vfree(v1x);
	       return(rc);
	     }
	 }
     });
  if (!def)
    {
      stackeq(v1x,v2);
      vfree(v1x);
    }
  else
    {
      v = adds(v2,v1x,MinusOne);
      vfree(v1x);
      fi = fldinv(f);
      smul(v,fi);
      Fclear(fi);
      stackd(b,x,v);
    }
#endif
  return(OK);
}
	  
      
void stackeq(v1,v2)			/* find the coincidence in v1 = v2
				 and stack it */
     vector v1,v2;
{
  vector v;
  LOG(logcoin,fprintf(logfile,"%sStacking equation ",logprefix);
      vprint(logfile,v1);
      fprintf(logfile,"%s                = ",logprefix);
      vprint(logfile,v2));
  SC(sccheckeq(v1,v2));
				/* First see if we are in one of
				 two easy cases */
  ASSERT(v1 == vroot(v1));
  ASSERT(v2 == vroot(v2));

  if (zerop(v2)) v = vcopy(v1);	/* two special cases */
  else if (zerop(v1)) v = vcopy(v2);
  else v = adds(v1,v2,MinusOne);
  if (zerop(v))
    vfree(v);
  else
    stackc(v);
  return;
}

retcode processc(v)			/* process the coincidence e_b = v */
				/* may free v */
     vector v;
{
  row r;
  vector v0 = vroot(v);		/* go to undeleted image */
  retcode rc,xrc;
#ifndef LEN
  menlist m,m1;
#ifdef PACKED_VEC
  pl p,p2;
#endif
#endif
  gpgen g;
  vector v1;
  basiselt b;
  fldelt f SAFE;
#ifdef INTEGRAL
  sent s;
  bool unit;
#else
  fldelt fac SAFE,ifac SAFE;
  Finit(fac);
#endif
  Finit(f);
  SC(sccheckc(v0));
  LOG(logcoin,
      fprintf(logfile,"%sProcessing 0 = ",logprefix);
      vprint(logfile,v);
      if (v0 != v)
      {
	log("Root converts to: ");
	vprint(logfile,v0);
      });
  if (v0 != v) vfree(v);
  if (zerop(v0))
    {
      vfree(v0);
      return(OK);
    }
#ifdef INTEGRAL
  if (!closed) /* can only get a coincidence if the unit is in last place */
    {
      s = BodyEndS(v0)-1;
      unit = false;
      if  (((s->fac).d[0]) == 1)
	if ((s->fac).size == 1)
	  {
	    smul(v0,MinusOne);
	    unit = true;
	  }
	else if ((s->fac).size == -1)
	  unit = true;
      if (unit)
	{
	  ASSERT(IsMOne(s->fac));
	  b = s->loc;
	  DelSent(v0,s);
	}
    }
  else /* can get a coincidence from any unit */
    {
      for (s = BodyEndS(v0)-1; s >= BodyS(v0); s--)
	{
	  unit = false;
	  if  (*((s->fac).d) == 1)
	    if ((s->fac).size == 1)
	      {
		smul(v0,MinusOne);
		unit = true;
	      }
	    else if ((s->fac).size == -1)
	      unit = true;
	  if (unit)
	    {
	      ASSERT(IsMOne(s->fac));
	      b = s->loc;
	      DelSent(v0,s);
	      break;
	    }
	}
    }
  /* Here, hopefully, either unit=true and b and v0 give a coincidence
     to apply, or unit=false and v0 needs to be added to the lattice */
	      
  if (!unit)
    {
      AddLatt(v0);
      return (OK);
    }
#else
  PS(v0,
     for (b = Len(v0)-1; b != (basiselt)-1; b--)
     if (!IsZero(BodyP(v0)[b]))
     {
       Fcopy(BodyP(v0)[b],fac);
       Fcopy(FZero,BodyP(v0)[b]);
       if (--Wt(v0))
	 {
	   Len(v0) = b;
	   while (IsZero(BodyP(v0)[--Len(v0)]))
	     ;
	   Len(v0)++;
	 }
       else
	 Len(v0)= 0;
       break;
     },
     b = Len(v0)-1;
     Fcopy((BodyEndS(v0)-1)->fac,fac);
     DelLast(v0);
     );
				 /* b is the one we're going to eliminate */
  if (b == (basiselt)-1)			/* trivial coincidence */
    {
      vfree(v0);
      return (OK);
    }

				/* now fac is the coefficient of the b entry
				 which has been removed, so that 
				 v = (v1-v2)-e_b.fac */
  Negate(fac,fac);
  ifac = fldinv(fac);
  Fclear(fac);
  v0 = smul(vtidy(v0),ifac); /* compute v/-fac */
  Fclear(ifac);
#endif /* OK here, regardless of field we ned to delete row b and replace it by v0 */
  ASSERT(!rowdel(b));
#ifdef INTEGRAL
  CoinLatt(b,v0); 
#endif
  r = table+b;
  rc = OK;
  for (g=0; g < nextgen; g++)
    if ((v1 = r->v.entries[g]))
      {
#ifdef LEN
	v1 = vroot(v1);
	if (v1 != r->v.entries[g])
	  {
	    vfree(r->v.entries[g]);
	    r->v.entries[g] = v1;
	  }
#endif
	xrc = stackxeq(v0,g,v1);
	if (xrc == CriticalOutOfSpace)
	  {
	    vfree(v0);
	    return xrc;
	  }
	ASSERT(xrc == OK);
	clrent(b,g);
	r = table+b;
      }
  if (rc != OK)		/* don't delete if we couldn't empty */
  {
    vfree(v0);
    return(rc);
  }
  r->status = Deleted;	/* mark it deleted */
  blanks -= nextgen;
  nalive--;
  myfree((pointer)r->v.entries);
  if (TrackPreIm && r->wlen)
      myfree((pointer)r->word);
  r->v.replacement = v0;
#ifdef INTEGRAL
  Fclear(r->latt.g);
#endif
#ifndef LEN
  /* now the search - start w the mentions */
  for (m = r->mentions; m; m=m1) 
    {
      sent s;
      basiselt bb;
      gpgen gg;
      v1 = elookup(m->en);
#ifdef PACKED_VEC
      ASSERT(Type(v1) == sparse);
#endif
      s = FindEnt(v1,b);
      ASSERT(s);
      Fcopy(s->fac,f);
      DelSent(v1,s);
      v1 = adds(v1,v0,f); /* now add in the replacement */
      m1 = m->next;
      bb = m->en.be;
      gg = m->en.g;
      clrent(bb, gg);	/* and install it */
      install(bb,gg,vtidy(v1));
      killment(m);
    }
  table[b].mentions = NULL;
#ifdef PACKED_VEC
  /* now the packed list */
  for (p=plisthead->next; p != plisthead; p= p2)
    {
      fldelt f SAFE;
      p2 = p->next;
      v1 = elookup(p->en);
      ASSERT(Type(v1) == packed);
      if (Len(v1) > b && !IsZero(BodyP(v1)[b])) /* a mention! */
	{
	  vector v3;
	  Fcopy(BodyP(v1)[b],f);
	  BodyP(v1)[b] = 0;	/* so remove it */
	  Wt(v1)--;
	  v3 = vtidy(adds(v1,v0,f));
	  if (Type(v3) == packed) /* and install the result */
	    {
	      v3->list.p = p;
	      elookup(p->en) = v3;
	      vfree(v1);
	    }
	  else
	    {
	      p->next->prev = p->prev;
	      p->prev->next = p->next;
	      elookup(p->en) = NULL;
	      blanks++;
	      install(p->en.be,p->en.g,v3);
	      killpl(p);
	      vfree(v1);
	    }
	}
    }
#endif
#endif
  Fclear(f);
  return(rc);
}


void clrent(b,g)			/* clear an entry */
     basiselt b;
     gpgen g;
{
  vector v = tabent(b,g);
#ifndef LEN
  PS(v,
    {
      pl p = v->list.p;
      p->next->prev = p->prev;
      p->prev->next = p->next;
      killpl(p);
    },
    {
      sent s; sent send;
      menlist m; menlist *m1;
      send = BodyEndS(v);
      for (s = BodyS(v); s < send; s++)
	{
	  m1 = &(table[s->loc].mentions);
	  while ((m = *m1))
	    {
	      if (m->en.be == b && m->en.g == g)
		{
		  *m1 = m->next;
		  killment(m);
		  break;
		}
	      m1 = &(m->next);
	    }
	}
    });
#endif
  vfree(v);
  tabent(b,g) = NULL;	/* kill the entry */
  blanks++;
}

vector vroot(v)			/* find the undeleted name of v */
      vector v;
{
  vector v3,v1 = v;
  bool moved = false;
  basiselt b;
  fldelt fac SAFE;
  Finit(fac);
  while (true)
    PS(v1,
       {
	 register basiselt l=Len(v1);
	 register fldelt *r=BodyP(v1);
	 register row rw=table;
	 for (; l; r++,rw++,l--)
	   if (!IsZero(*r) && rw->status == Deleted)	/* something to do */
	     {
	       b = rw-table;
	       LOG(logcoin > 1,
		   fprintf(logfile,"%sIn root at %d changing ",logprefix,b);
		   vprint(logfile,v1); );
	       if (!moved) 
		 {
		   v1 = vcopy(v);
		   moved = true;
		   r = BodyP(v1)+b;
		 }
	       Fcopy(*r,fac);
	       Fcopy(FZero,*r);
	       Wt(v1)--;
	       v3 = vroot(replace(b));
	       if (v3 != replace(b))
		 {
		   vfree(replace(b));
		   replace(b) = v3;
		 }
	       v3 = adds(v1,v3,fac);
	       vfree(v1);
	       v1 = v3;
	       LOG(logcoin > 1,
		   log("to ");
		   vprint(logfile,v1););
	       break;
	     }
	 if (!l) break;
       },
       {
	 register sent s;
	 register sent send;
	 register row mytable=table;
	 send = BodyEndS(v1);
	 for (s = BodyS(v1); s < send; s++)
	  {
	    if (mytable[s->loc].status == Deleted) /* something to do */
	      {
		LOG(logcoin > 1,
		    fprintf(logfile,"%sIn root at %d changing ",logprefix,s->loc);
		    vprint(logfile,v1););
		if (!moved)
		  {
		    v1 = vcopy(v);
		    s = BodyS(v1)+(s-BodyS(v));
		    moved = true;
		    send = BodyEndS(v1);
		  }
		Fcopy(s->fac,fac);
		b = s->loc;
		DelSent(v1,s);
		v3 = vroot(replace(b));
                if (v3 != replace(b))
                {
                   vfree(replace(b));
                   replace(b) = v3;
                } 
		v3 = adds(v1,v3,fac);
		vfree(v1);
		v1 = v3;
		LOG(logcoin > 1, log("to "); vprint(logfile,v1););
		break;
	      }
	  }
	 if (s == send) break;
       });
  Fclear(fac);
  return(v1);
}

retcode processd(b,g,v)			/* process the deduction b.g = v */
     basiselt b;
     gpgen g;
     vector v;
{
  vector v1,v0 = vroot(v);
  retcode rc;
  SC(sccheckd(b,g,v));
  LOG(logcoin,
      fprintf(logfile,"%sProcessing e%d.%s = ",logprefix,b,gennames[g]);
      vprint(logfile,v);
      if (v0 != v)
	{
	  log("Root converts to: ");
	  vprint(logfile,v0);
	});
  if (v0 != v) vfree(v);
  if (rowdel(b))		/* row deleted already */
    {
      v1 = vroot(replace(b));
      if (v1 != replace(b))
	{
	  vfree(replace(b));
	  replace(b) = v1;
	}
      rc = stackxeq(v1,g,v0);
      vfree(v0);
      return(rc);
    }
  if ((v1 = tabent(b,g)))	/* entry already full */
    {
#ifdef LEN
      v1 = vroot(tabent(b,g));
      if (v1 != tabent(b,g))
	{
	  vfree(tabent(b,g));
	  tabent(b,g) = v1;
	}
#endif
      stackeq(v0,v1);
      vfree(v0);
     }
  else
    install(b,g,vtidy(v0));
  return(OK);
}

void install(b,g,v)			/* install a known trouble-free
				  entry */
     basiselt b;
     gpgen g;
     vector v;
{
  blanks--;
  tabent(b,g) = v;
#ifndef LEN
  PS(v,{
    pl p;
    newpl(p);
    p->next = plisthead->next;
    p->prev = plisthead;
    p->next->prev = p;
    p->prev->next = p;
    p->en.be = b;
    p->en.g = g;
    v->list.p = p;
  }, {
    sent s; sent send;
    menlist m;
    send = BodyEndS(v);
    for (s=BodyS(v); s < send; s++)
      {
	newment(m);
	m->en.be = b;
	m->en.g = g;
	m->next = table[s->loc].mentions;
	table[s->loc].mentions = m;
      }
  });
#endif
}

static void DelCoin(c)
     coinl c;
{
  coinl cn;
  if (!c)
    return;
  vfree(c->body.v);
  cn = c->next;
  myfree((pointer)c);
  DelCoin(cn);
}

static void DelDeduct(d)
     deductl d;
{
  deductl dn;
  if (!d)
    return;
  vfree(d->body.v);
  dn = d->next;
  myfree((pointer)d);
  DelDeduct(dn);
}
  


  
retcode dostacks()			/* clear both stacks */
{
  retcode xrc,rc=OK;
#ifdef DEBUG
	doingstacks = true;
#endif
  while (1)
    {
      if (coincidences)		/* coincidences first */
	{
	  coinl c = coincidences;
	  vector v = c->body.v;
	  coincidences = c->next;
	  myfree((pointer)c);
	  xrc = processc(v);
	  if (xrc != OK)
	    {
	      rc = xrc;
	      if (xrc == CriticalOutOfSpace)
		break;
	    }
	}
      else if (deductions)
	{
	  deductl d = deductions;
	  basiselt b = d->body.b;
	  gpgen g = d->body.g;
	  vector v = d->body.v;
	  deductions = d->next;
	  myfree((pointer)d);
	  xrc = processd(b,g,v);
	  if (xrc != OK)
	    {
	      rc = xrc;
	      if (xrc == CriticalOutOfSpace)
		break;
	    }
	}
      else break;
    }
#ifdef DEBUG
	doingstacks = false;
#endif
  if (rc == CriticalOutOfSpace)
    {
      DelCoin(coincidences);
      DelDeduct(deductions);
    }
  return(rc);
}
  
      
	  
	
  


