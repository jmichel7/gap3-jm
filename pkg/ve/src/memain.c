/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */
#include "me.h"
#include <ctype.h>


vector MultiImage(wl)
	wordl wl;
{	
  vector v1,v2,v3;
  basiselt b;
  retcode rc;
  swl s;
  LOG(logact > 1,
      log("Computing multi-image of ");
	  printwordl(logfile,wl);
	  fputc('\n',logfile););
  newsvec(v3,0);
  if (wl->type == wl_packed)
    {
      for (b=0; b < startvecs; b++)
	{
	  if (wl->body.p.wds[b]->type != scalar || !IsZero(wl->body.p.wds[b]->body.scalar))
	    {
	      v2 = btov(b);
	      v1 = vroot(v2);
	      if (v2 != v1) vfree(v2);
	      v2 = image(v1,wl->body.p.wds[b],&rc);
	      switch(rc)
		{
		case OK:
		  break;
		case OutOfSpace:
		case CriticalOutOfSpace:
		  fprintf(stderr,"Out of space processing submodule\n");
		  vfree(v1);
		  vfree(v3);
		  return NULL;
		default:
		  DIE();
		}
	      vfree(v1);
	      v1 = add(v3,v2);
	      LOG(logact > 2,
		  fprintf(logfile,"%sPartial image %u: ",logprefix,b);
		  vprint(logfile,v2);
		  log("Running total: ");
		  vprint(logfile,v1););
	      vfree(v3);
	      vfree(v2);
	      v3 = v1;
	    }
	}
    }
  else
    {
      for (s = wl->body.s; s; s = s->next)
	{
	  v2 = btov(s->loc);
	  v1 = vroot(v2);
	  if (v2 != v1) vfree(v2);
	  v2 = image(v1,s->w,&rc);
	  switch(rc)
	    {
	    case OK:
	      break;
	    case OutOfSpace:
	    case CriticalOutOfSpace:
	      fprintf(stderr,"Out of space processing submodule\n");
	      vfree(v1);
	      vfree(v3);
	      return NULL;
	    default:
	      DIE();
	    }
	  vfree(v1);
	  v1 = add(v3,v2);
	  LOG(logact > 2,
	      fprintf(logfile,"%sPartial image at %u: ",logprefix,s->loc);
	      vprint(logfile,v2);
	      log("Running total: ");
	      vprint(logfile,v1););
	  vfree(v3);
	  vfree(v2);
	  v3 = v1;
	}
    }
  LOG(logact > 1,
    log("Returning from multi-image: ");
    vprint(logfile,v3););
  return (v3);
}

retcode DoSubGens()
{
  rell r;
  smgl s;
  vector v1,v2;
  retcode rc;
  basiselt b;
#ifdef DEBUG
  insubgens = true;
#endif
  if (subgens->type == rels)
    for (r = subgens->list.r; r; r=r->next)
      {
	if (r->type == algebra)
	  rc = push(r->entry.a.l,r->entry.a.r,0, false);
	else
	  rc = pushg(r->entry.g,0);
	switch(rc)
	  {
	  case OK:
	    break;
	  case OutOfSpace:
	  case CriticalOutOfSpace:
	    fprintf(stderr,"Out of space processing submodule\n");
	    return(rc);
	  default:
	    DIE();
	  }
	if (TimedOut)
	  return(OutOfTime);
      }
  else
    {
      for (s = subgens->list.s; s; s = s->next)
	{
	  if (s->type == normal)
	    {
	      v1 = MultiImage(s->body.norm.lhs);
	      if (s->body.norm.rhs)
	   	v2 = MultiImage(s->body.norm.rhs);
	      else
		{
		  newsvec(v2,0);
		}
	      stackeq(v1,v2);
	      vfree(v1);
	      vfree(v2);	   
	      rc = dostacks();
	      switch(rc)
		{
		case OK:
		  break;
		case OutOfSpace:
		case CriticalOutOfSpace:
		  fprintf(stderr,"Out of space processing submodule\n");
		  return(rc);
		default:
		  DIE();
		}
	    }
	  else
	    {
	     for (b = s->body.univ.start; b <= s->body.univ.end; b++)
	       {
		 rc = push(s->body.univ.w,GawZero,b, true);
		 switch(rc)
		   {
		   case OK:
		     break;
		   case OutOfSpace:
		   case CriticalOutOfSpace:
		     fprintf(stderr,"Out of space processing submodule\n");
		     return(rc);
		   default:
		     DIE();
		   }
	       }
	   }
	  if ((nextbe > 2*nalive && nextbe > 3*startvecs) ||
	      TABLE_CEILING - nextbe < nextbe - nalive)
	    pack();
	  if (TimedOut)
	    return(OutOfTime);
	}
    }
#ifdef DEBUG
	insubgens = false;
#endif
  return OK;
}

retcode  DoLookahead(wmin,wmax)
     weight wmin,wmax;
{
  weight w;
  rell r;
  retcode rc;
#ifdef INTEGRAL
  gpgen g;
#endif
  LOG(logwt,
    log("Looking ahead ...\n"););
  for (r = relators; r; r= r->next)
    r->SavedNext = r->NexttoDo;
#ifdef INTEGRAL
  for (g = 0; g < nextgen; g++)
    SavedNext[g] = NextToDo[g];
#endif
  LookAhead = true;
  for (w = wmin; w < wmax; w++)
    {
      rc = pushrl(relators,w);
      if (rc == NoMoreWork)
	break;
      if (rc == CriticalOutOfSpace)
	return rc;
      if (TimedOut)
	return(OutOfTime);
      ASSERT (rc == OK);
      LOG(true, fflush(logfile););
      FirstofWt[w+1] = nextbe;
#ifdef INTEGRAL
      LFoW[w+1] = nextlat;
#endif
      if (!blanks)
	{
#ifdef INTEGRAL
	  if (!closed)
	    {
	      closed = true;
	      LattClosing();
	    }
#endif
	  if (EC)
	    {
	      if (nalive > ecmin && nalive < ecmax)
		{
		  pack();
		  break;
		}
	    }
	}
    }
  for (w = wmin+1; w <= maxwt; w++)
    {
      FirstofWt[w] = (basiselt)-1;
#ifdef INTEGRAL
      LFoW[w] = (lattpos)-1;
#endif
    }
  for (r = relators; r; r = r->next)
    r->NexttoDo = r->SavedNext;
#ifdef INTEGRAL
  for (g=0; g < nextgen; g++)
    NextToDo[g] = SavedNext[g];
#endif
  threshold = 2*nalive;
  LookAhead = false;
  LOG(logwt, log("  ...done\n"););
  return OK;
}


static retcode runme()
{
  bool norels;
  weight w;
  retcode rc;
  basiselt b;
  bool inec = false;
  
/* Can't tidy this up because it's two layered */
#ifdef LOGGING
  if (logstages>1)
    {
      log("Submodule gens\n");  
      printsmg(logfile,logprefix);
      log("Relator\n");
      printrl(logfile,logprefix,relators);
    }
  else if (logstages)
    log("Read Input\n");
#endif
#ifdef DEBUG
  LOG(logmem, ReportAllocs(); vcheck(););
#endif
  table = (row) myalloc(tablesize*sizeof(_row),false); /* initialize table */
  for (b=0; b<tablesize; b++)
    {
      table[b].v.entries = NULL;
      table[b].status = OutOfUse;
    }
  for (nextbe = 0; nextbe < startvecs; nextbe++)
  {
    if (TrackPreIm)
    {
    	table[nextbe].fmgen = nextbe;
    	table[nextbe].wlen = 0;
    	table[nextbe].word = NULL;
    }
    initrow(nextbe,DefinesOK);
  }
  FirstofWt = (basiselt *) myalloc((maxwt+1)*sizeof(basiselt),false);
  FirstofWt[1] = 0;
  FirstofWt[2] = startvecs;
  FirstofWt[3] = nextbe;
  for (w=4; w< maxwt; w++)
    FirstofWt[w] = (basiselt)-1;
#ifdef INTEGRAL
  LFoW = (lattpos *) myalloc((maxwt+1)*sizeof(lattpos),false);
  LFoW[1] = 0;
  LFoW[2] = 0;
  LFoW[3] = nextlat;
  for (w=4; w< maxwt; w++)
    LFoW[w] = (lattpos)-1;
#endif
  rc = DoSubGens();
  if (rc != OK)
    return rc;
  LOG(logstages,log("Done submodule generators\n"););
  LOG(logstages > 3, printtable(logfile,logprefix););
#ifdef DEBUG
  LOG(logmem, ReportAllocs(); vcheck(););
#endif
  if (!threshold) threshold = 2*nextbe;
  for (w=2; w<maxwt; w++)
    {
    TryAgain:
      rc = pushrl(relators,w);
      if (rc == OutOfSpace)
	{
	  rc = DoLookahead(w,maxwt);
	  if (rc != OK)
	    {
	      ASSERT(rc == OutOfTime || rc == CriticalOutOfSpace);
	      return(rc);
	    }
	  if (nalive >= tablesize || nalive >= UnreservedLimit)
	    {
	      rc = OutOfSpace;
	      break;
	    }
	  goto TryAgain;
	}
      norels = (rc == NoMoreWork);
#ifndef INTEGRAL
      if (norels)
	break;
      if (rc == OutOfSpace || rc == CriticalOutOfSpace)
	return(rc);
#else
    LTryAgain:
      rc = PushLatt(w);
      if (rc == OutOfSpace)
	{
	  DoLookahead(w,maxwt);
	  if (nalive >= tablesize)
	    break;
	  goto LTryAgain;
	}
      if (rc == CriticalOutOfSpace)
	return rc;
      if (norels && rc == NoMoreWork)
	break;
      LFoW[w+1] = nextlat;
#endif
      LOG(true,fflush(logfile));
      FirstofWt[w+1] = nextbe;
#ifdef INTEGRAL
      LFoW[w+1] = nextlat;
#endif
      if (!blanks)
	{
#ifdef INTEGRAL
	  if (!closed)
	    {
	      closed = true;
	      LOG(logstages > 2, log("Table closed, more flexible rules apply\n"));
	      LattClosing();
	    }
#endif
	  if (EC)
	    {
	      if (nalive > ecmin && nalive < ecmax)
		{
		  pack();
		  LOG(logstages, log ("Early Closing\n"););
		  inec = true;
		  break;
		}
	    }
	}
      
      if (LA && nalive >= threshold)
	{
	  rc = DoLookahead(w+1,min(maxwt,w+LookWeight+1));
	  if (rc != OK)
	    {
	      ASSERT(rc == OutOfTime || rc == CriticalOutOfSpace);
	      return(rc);
	    }
	  pack();
	}
#ifdef DEBUG
      LOG(logmem > 1, ReportAllocs(); vcheck(););
#endif	 
    }
#ifdef DEBUG
  LOG(logmem > 1, ReportAllocs(); vcheck(););
#endif	 
  if (!inec && rc != NoMoreWork)
    return(rc);
  LOG(logstages, fprintf(logfile,"%sClosed, %d rows defined\n",logprefix,nextbe ););
  pack();
#ifdef DEBUG
  LOG(logmem > 1, ReportAllocs(); vcheck(););
#endif
  LOG(logstages, fprintf(logfile,"%s%d live dimensions\n",logprefix,nalive););
  LOG(logstages > 2, printtable(logfile,logprefix););
  return OK;
}

#define LOG_INTERVAL 100

static retcode runquot()
{
  basiselt b,nsmg;
  vector v;
  unsigned int p,x;
  gpgen g;
  fscanf(ifp,"%u%u",&p,&x); /* get the characteristic and dimension */
#ifdef GFP
  if (!ForceCharacteristic)
    {
      if (!p || !IsPrime(p))
	{
	  fprintf(stderr,"Characteristic %u not permitted\n",p);
	  DIE();
	}
      P = p;
    }
#else
  if (!ForceCharacteristic && p != 0)
    {
      fprintf(stderr,"Characteristic must be zero with this program, not %u\n",p);
      DIE();
    }
#endif
  nextbe = x;
  startvecs = nextbe;
  readgens();
  inverse = (gpgen *)myalloc(sizeof(gpgen)*nextgen,false);
  for (g = 0; g < nextgen; g++)
  	inverse[g] = NOINVERSE;
  table = (row)myalloc(sizeof(_row)*nextbe,false);
  nalive = blanks = 0;
  for (b=0; b < nextbe; b++)
      initrow(b,DefinesOK);
  fscanf(ifp,"%u",&x);
  nsmg = x;
  for (b = 0; b < nsmg; b++)
   {
      v = readvec(ifp);
      if (!v)
	DIE();
      if (!zerop(v))
	 stackc(v);
   }
   dostacks();
  LOG(logstages,
      fprintf(logfile,
	      "%sRead submodule generators %u live dimensions\n",logprefix,nalive);
      fflush(logfile););
  for (b=0; b < nextbe; b++)
  {
     for (g=0; g < nextgen; g++)
     {
         v = readvec(ifp);
	 if (!v)
	   DIE();
         stackd(b,g,v);
     }
     dostacks();
     if (TimedOut)
       return(OutOfTime);
     LOG(logstages > 1 && !(b % LOG_INTERVAL),
	 fprintf(logfile,
		 "%sRead %u sets of images\n\t%5u alive%6u blanks\n",
		 logprefix,b,nalive,blanks); fflush(logfile););
  }
  LOG(logstages,fprintf(logfile,"%sAll read in %u alive\n",logprefix,nalive);fflush(logfile););
  return(OK);
}


FILE *yyin;
int main(argc,argv)		       
     int argc;
     char ** argv;
{
  retcode rc;
/*   clock(); */
/*  mallopt(M_MXFAST,256);
  mallopt(M_GRAIN,16); */
  yyin = stdin;
#if defined(RATIONAL) || defined(INTEGRAL)
  mp_set_memory_functions(alloc_for_gmp,realloc_for_gmp,free_for_gmp);
#endif
#ifdef RATIONAL
  Finit(FOne);
  mpq_set_ui(&FOne,1,1);
  Finit(FZero);
  mpq_set_ui(&FZero,0,1);
  Finit(MinusOne);
  Negate(FOne,MinusOne);
#endif
#ifdef INTEGRAL
  mpz_init_set_ui(&(FOne),1);
  mpz_init_set_ui(&(FZero),0);
  mpz_init_set_si(&(MinusOne),-1);
#endif
  SetDefaults();
  comline(argc,argv);
  if (!IsQuot)
    {
      getinput();
      SC(scsetup());
      rc = runme();
    }
  else
    rc = runquot();
  switch(rc)
    {
    case OK:
      if (ofn)
	mxsout();
	  break;
    case OutOfSpace:
      LOG(logstages,log("Out of space\n"));
      break;
    case CriticalOutOfSpace:
      LOG(logstages,log("Out of space for critical define\n"));
      break;
    case OutOfTime:
      LOG(logstages,log("Out of time\n"));
      break;
    case OutOfWeights:
      LOG(logstages,log("Out of weights\n"));
      break;
    default:
      DIE();
    }
  CleanUp();
  return(rc != OK);
}


