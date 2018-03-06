/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */
#include <memory.h>
#include "me.h"

retcode pushg(gw,be)			/* Push a group word */
     gpwd gw;
     basiselt be;
{
  vector vf=btov(be),vb=btov(be),v1; /* like SCE, front and back ptrs */
  gpgen g;
  int f,b,l=gw->len;
  retcode rc;
#ifdef DEBUG
  inpushg = true;
#endif
  LOG(logpush,
      char s[120];
      gpwdtos(s,gw);
      fprintf(logfile,"%sPushing %s from %d\n",logprefix,s,(int) be););
  if (!gw->len)
    return OK;
  for (f=0; f < l; f++)		/* go forwards */
    {
      v1 = action(vf,gw->word[f],NoDefines,&rc);
      if (!v1) 
	{
	  if (rc != NeedToDefine)
	    {
	      vfree(vf);
	      vfree(vb);
#ifdef DEBUG
	      inpushg = false;
#endif	     
	      return(rc);
	    }
	  else
	    break;
	}
      vfree(vf);
      vf = v1;
    }
  if (f == l)			/* we got round */
    stackeq(vf,vb);
  else
    {
      for (b = l-1; b > f; b--)	/* go backwards */
	{
	  v1 = action(vb,inverse[gw->word[b]],LookAhead ? NoDefines : DefinesOK,&rc);
	  /* define allowed now  unless we're in Lkahd*/
	  if (!v1)		/* Failed in Lkahd or table full */
	    {
	      vfree(vf);
	      vfree(vb);
#ifdef DEBUG
	      inpushg = false;
#endif
	      return(rc);
	    }
	  vfree(vb);
	  vb = v1;
	}
      g = gw->word[f];
      v1 = action(vb,inverse[g],NoDefines,&rc);	/* see if we can close the 
						   last gap without defining */
      if (v1)
	{
	  stackeq(vf,v1);
	  vfree(v1);
	}
      else if (rc == NeedToDefine)
	  rc = stackxeq(vf,g,vb);
    }
  vfree(vf);
  vfree(vb);
  if (rc == OK) 
    rc = dostacks();		/* clear the stacks */
#ifdef DEBUG
	inpushg = false;
#endif
  return(rc);
}
  

retcode pushrl(rl,w)		/* push a complete list of relators 
				   at a given weight. Returns true
				   if there are no pushes left to do at
				   this or any higher weight */
     rell rl;
     weight w;
{
  weight wt;
  rell r;
  retcode rc;
  LOG(logwt,
    fprintf(logfile,
	    "%sStarting weight %d in %s mode, %u alive out of %u\n",
	    logprefix,
	    w,
	    LookAhead ? "lookahead" : "define",
	    nalive,
	    nextbe););
  for (r = rl; r; r=r->next)
    {
      wt = LookAhead ? r->lwt : r->wt;
      if (wt < w)
	{
	  while (r->NexttoDo < FirstofWt[w-wt+1])
	    if (!rowdel(r->NexttoDo++))
	      {
		if (r->type == algebra)
		  rc = push(r->entry.a.l,r->entry.a.r,r->NexttoDo-1, false);
		else
		  rc = pushg(r->entry.g,r->NexttoDo-1);
		switch (rc)
		  {
		  case RowDeleted:
		  case OK:
		  case NeedToDefine: /* must be in lookahead */
		    break;
		  case OutOfSpace:
		    r->NexttoDo--;
		    if (nalive >= tablesize || nalive >= UnreservedLimit)
		      return(OutOfSpace);
		    pack();
		    break;
		  case CriticalOutOfSpace:
		    return CriticalOutOfSpace;
		  default:
		    fprintf(stderr,"Bad rc from push %d\n",rc);
		    abort();
		  }
	      }
	}

    }
  for (r = relators; r; r = r->next)
    if (r->NexttoDo < nextbe)
      return(OK);
  return(NoMoreWork);
}


retcode push(w1,w2,b, inuniv)			/* push a g.a.w. */
     gaw w1,w2;
     basiselt b;
     bool inuniv;
{
  vector v1,v2,v3;
  retcode rc;
#ifdef DEBUG
  inpusha = true;
#endif
  if (rowdel(b) && !inuniv)
  {
#ifdef DEBUG
    inpusha = false;
#endif
    return(RowDeleted);
  }
  LOG(logpush,
    {
      log("Pushing ");
      printgaw(logfile,w1);
      if (w2)
      {
      	fputc('=',logfile);
      	printgaw(logfile,w2);
      }
      fprintf(logfile," from %d\n",b);
    });
  v1 = btov(b);			/* easy this way, as we just have to define
				 everywhere, so there's no complication */
  if (inuniv && rowdel(b))
    {
      v2 = vroot(v1);
      vfree(v1);
      v1 = v2;
    }
  v2 = image(v1,w1,&rc);
  if (rc != OK)
  {
    vfree(v1);
#ifdef DEBUG
	inpusha = false;
#endif
    return(rc);
  }
  if (w2)
  {
    v3 = image(v1,w2,&rc);
    vfree(v1);
    if (rc != OK)
    {
       vfree(v2);
#ifdef DEBUG
	inpusha = false;
#endif
      return(rc);
   }
  }
  else
    v3 = v1;
  LOG(logpush > 1,
    {
      log("Images ");
      vprint(logfile,v2);
      log("       ");
      vprint(logfile,v3);
    });
  stackeq(v2,v3);
  vfree(v3);
  vfree(v2);
  rc = dostacks();
#ifdef DEBUG
        inpusha = false;
#endif
  return(rc);
}


vector image(v,w,rc)		/* image of v under gp alg. word w */
     vector v;
     gaw w;
     retcode *rc;
{
  int i;
  gpwd g;
  vector v1,v2,v3,vr;
  LOG(logact,
      log("");
      printgaw(logfile,w);
      fputs(" acting on ",logfile);
      vprint(logfile,v););
  switch(w->type)
    {
    case scalar:
      if (IsZero(w->body.scalar))
	{
	  newsvec(vr,0);
	}
      else
	{
	  v1 = vcopy(v);
	  vr = smul(v1,w->body.scalar);
	}
      break;

    case grp:
      g = w->body.group;
      if (g->len)
	{
	  v1 = v;
	  for (i = 0; i < g->len; i++)
	    {
	      v2 = action(v1,g->word[i],LookAhead ? NoDefines : DefinesOK,rc);
	      if (i) vfree(v1);
	      if (!v2)
		return(NULL);
	      v1 = v2;
	    }
	  vr = v1;
	}
      else
	vr = vcopy(v);
      break;
    case sum:
      v1 = image(v,w->body.sum.a,rc);
      if (!v1)
	return(NULL);
      v2 = image(v,w->body.sum.b,rc);
      if (!v2)
      {
      	vfree(v1);
	return(NULL);
      }
      v3 = add(v1,v2);
      vfree(v1);
      vfree(v2);
      vr = v3;
      break;
    case product:
      v1 = image(v,w->body.prod.c,rc);
      if (!v1) return(NULL);
      v2 = image(v1,w->body.prod.d,rc);
      vfree(v1);
      if (!v2)
	return(NULL);
      vr = v2;
      break;
    }
  LOG(logact,
      log("Producing ");
      vprint(logfile,vr););
  *rc = OK;
  return(vr);
}

vector action(v,g,def,rc)	/* Action of group generator g on v
				 defining new basis elts is allowed if
				 def is true */
     vector v;
     gpgen g;
     DefineStatus def;
     retcode *rc;
{
  vector vg,vg1,v1;
  int status = 0;
#ifdef DEBUG
  inact = true;
#endif
  LOG(logact,
      fprintf(logfile,"%s%s acting on ",logprefix,gennames[g]);
      vprint(logfile,v);
      );
  PS(v,
      basiselt i;
      basiselt l;
      fldelt *r;
      *rc = OK;
      l = Len(v);
      r = BodyP(v);
      for (i = 0; i<l; r++,i++)
	if (!IsZero(*r))
	  {
	    v1 = baction(i,g,def,rc);
	    if (!v1)		/* couldn't do it -- def false */
	      {
		if (status == 2) vfree(vg);
		IFDEBUG(inact = false);
		return(NULL);
	      }
	    switch(status)
	      {
	      case 0:
		if (IsOne(*r))
		  {
		    vg = v1;
		    status++;
		  }
		else
		  {
		    vg = vcopy(v1);
		    smul(vg,*r);
		    status += 2;
		  }
		break;

	      case 1:
		vg = adds(vg,v1,*r);
		status++;
		break;
		
	      case 2:
		vg1 = adds(vg,v1,*r);
		vfree(vg);
		vg = vg1;
		break;
		
	      default:
		DIE();
	      }
	  },
      sent s; sent send;
      send = BodyEndS(v);
      for (s = BodyS(v); s < send; s++)
	{
	  v1 = baction(s->loc,g,def,rc);
	  if (!v1)		/* couldn't do it again */
	    {
	      if (status == 2)
		vfree(vg);
	      IFDEBUG(inact = false);
	      return(NULL);
	    }
	  switch(status)
	    {
	    case 0:
	      if (IsOne(s->fac))
		{
		  vg = v1;
		  status++;
		}
	      else
		{
		  vg = vcopy(v1);
		  smul(vg,s->fac);
		  status = 2;
		}
	      break;

	    case 1:
	      vg = adds(vg,v1,s->fac);
	      status++;
	      break;

	    case 2:
	      vg1 = adds(vg,v1,s->fac);
	      vfree(vg);
	      vg = vg1;
	      break;
	      
	    default:
	      DIE();
	    }
	});
  switch(status)
    {
    case 0:
      newsvec(vg,0);
      break;
    case 1:
      vg = vcopy(vg);
      break;
    case 2:
      break;
    default:
      DIE();
    }
  LOG(logact,
      log("Producing ");
      vprint(logfile,vg););
  IFDEBUG( inact = false);
  return(vg);
}
  

vector baction(b,g,def,rc)	/* action of generator on basis elt
				 again def allows or forbids new basis elts*/
     basiselt b;
     gpgen g;
     DefineStatus def;
     retcode *rc;
{
  vector v;
  basiselt b1;
  LOG(logact,
      fprintf(logfile,
	      "%sBasic action of %s on %u\n",
	      logprefix,
	      gennames[g],
	      (unsigned int) b););

  if ((v=tabent(b,g)))
    {
#ifdef LEN
      v =vroot(v);
      if (v != tabent(b,g))
	{
	  vfree(tabent(b,g));
	  tabent(b,g) = v;
	}
#endif
      LOG(logact,
	  log("Returning ");
	  vprint(logfile,v););
      *rc = OK;
      return(v);
    }
  else if (def != NoDefines)
    {
      *rc = initrow(nextbe,def);
      if (*rc != OK)
	return(NULL);
      SC(scdef(b,g,nextbe));
      LOG(loginit > 1,
	  fprintf(logfile,
		  "%s  as image of %u under %s\n",
		  logprefix,
		  (unsigned int)b,
		  gennames[g]););
      v = btov(nextbe);
      install(b,g,v);	/* set up new row */
      if (TrackPreIm)
      {
      	b1 = table[b].wlen;
      	table[nextbe].fmgen = table[b].fmgen;
      	table[nextbe].wlen = b1+1;
      	table[nextbe].word = (gpgen *)myalloc(sizeof(gpgen)*(b1+1),false);
      	memcpy((pointer)table[nextbe].word,
               (pointer)table[b].word,
               sizeof(gpgen)*b1);
        table[nextbe].word[b1] = g;       
	LOG(loginit > 2,int ii;
	    log("      ");
	    for (ii=0; ii < table[nextbe].wlen; ii++) {
	    fputs(gennames[table[nextbe].word[ii]],logfile);
	    fputc(' ',logfile);}
	    fputc('\n',logfile);
	    );
      }
      if (inverse[g] != NOINVERSE)
	 install(nextbe,inverse[g],btov(b));
      nextbe++;
      LOG(logact,
	  log("Returning ");
	  vprint(logfile,v););
      return(v);
    }
  else
    {
      LOG(logact,log("Returning null\n"););
      *rc = NeedToDefine;
      return(NULL);
    }
}
	
retcode initrow(b,stat)			/* initialize a new row */
     basiselt b;
     DefineStatus stat;
{
  gpgen g;
  basiselt ots,b1;
  row ntab;
  LOG(loginit,fprintf(logfile,"%sInitializing row %d\n",logprefix,b););
  ASSERT(stat != NoDefines);
  if (stat == DefinesOK && b > UnreservedLimit)
    return OutOfSpace;      
  if (b >= tablesize)
    {
      if (tablesize >=  TableLimit )
	return CriticalOutOfSpace; /* stat must be CriticalDefines here */
      ots = tablesize;
      tablesize = min(tablesize+TABLE_BLOCK,TableLimit);
      ntab = (row)myrealloc(table,sizeof(_row)*tablesize,true);
      if (!ntab)
	return(OutOfSpace);
      table = ntab;
      for (b1 = ots; b1 < tablesize; b1++)
	{
	  table[b1].v.entries = NULL;
	  table[b1].status = OutOfUse;
	}
    }
  table[b].status = InUse;
#ifndef LEN
  table[b].mentions = NULL;
#endif
  table[b].v.entries = (vector *) myalloc(nextgen * sizeof(vector),false);
#ifdef INTEGRAL
  table[b].latt.mentions = NULL;
  Finit(table[b].latt.g);
#endif
  for (g=0; g < nextgen; g++)
      tabent(b,g) = NULL;
  blanks += nextgen;
  nalive++;
  return OK;
}



