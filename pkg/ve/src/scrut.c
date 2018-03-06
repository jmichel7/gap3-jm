/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */
/* this is debugging code */



#include "me.h"
#include <string.h>
#ifdef SCRUT			/* compile nothing unless this is defined */
vector **sctable;		/* target table */
vector *trans;			/* translate table for the basis vectors
				 of the enumerator */
bool scrut = false;		/* scrut on/off */
char *sfile = "meout";		/* file name to read matrices */

extern vector scact();		/* two subroutines which are not exported */
extern vector scv();

static basiselt scdimen=0;
static basiselt sctabsize;


#ifdef GFP
static fldelt scrdfelt(f,b)
     FILE *f;
     bool *b;
{
  fldelt x;
  *b = 1 == fread((pointer)&x,sizeof(fldelt),1,f);
  return x;
}
#endif

#ifdef RATIONAL
static fldelt scrdfelt(f,b)
     FILE *f;
     bool *b;
{
  fldelt x;
  Finit(x);
  mpz_inp_raw(&(x.num),f);
  mpz_inp_raw(&(x.den),f);
  *b = true; /* can't tell whether we succeeded */
  return x;
}

#endif

#ifdef INTEGRAL
static fldelt scrdfelt(f,b)
     FILE *f;
     bool *b;
{
  fldelt x;
  Finit(x);
  mpz_inp_raw(&(x),f);
  *b = true; /* can't tell whether we succeeded */
  return x;
}

#endif

static vector scrdv(f)
     FILE *f;
{
  fldelt x;
  bool readok;
  sent s;
  basiselt b,last=(basiselt)-1;
  vector v;
  newsvec(v,scdimen);
  s = BodyS(v);
  for (b=0; b < scdimen; b++)
    {
      x = scrdfelt(f,&readok);
      if (!readok)
	{
	  fprintf(stderr,"Can't read vector from scrut file\n");
	  vfree(v);
	  return(NULL);
	}
      if (!IsZero(x))
	{
	  s->loc = b;
	  Fclear(s->fac);
	  s->fac = x;
	  last = b;
	  Wt(v)++;
	  s++;
	}
      else
	Fclear(x);
    }
  Len(v) = last+1;
  return v;
}


#ifndef INTEGRAL
static void scinvert(g,h)
     gpgen g,h;
{
  vector *extra,v,v1;
  basiselt b,b1;
  sent s;
  fldelt x,x1,xi;
  /* fill in sctable for g from that for h (g = h^-1) */
  extra = (vector *)myalloc(sizeof(vector)*scdimen,false);
  for (b=0; b < scdimen; b++)
    {
      extra[b] = vcopy(sctable[b][h]);
      sctable[b][g] = btov(b);
    }
  for (b=0; b < scdimen; b++)
    {
      v = extra[b];
      PS(v,
	 x = (Len(v) > b ? BodyP(v)[b] : FZero);,
	 s = FindEnt(v,b);
	 x = (s) ? s->fac : FZero;);
      if (IsZero(x)) /* got a zero on the diagonal so go hunting
		 for a non-zero in this column */
	{
	  for (b1 = b+1; b1 < scdimen; b1++)
	    {
	      v1 = extra[b1];
	      PS(v1,
		 x = (Len(v1) > b ? BodyP(v1)[b] : FZero);,
		 s = FindEnt(v1,b);
		 x = (s) ? s->fac : FZero;);
	      if (!IsZero(x))
		break;
	    }
	  ASSERT(b1 < scdimen);
	  extra[b] = v1;
	  extra[b1] = v;
	  v = sctable[b1][g];
	  sctable[b1][g] = sctable[b][g];
	  sctable[b][g] = v;
	  v = v1;
	}
      /* Now x <> 0 and x is the b component of v,
	 which is extra[b] */
      xi = fldinv(x);
      smul(v,xi);
      smul(sctable[b][g],xi);
      Fclear(xi);
      Finit(x1);
      for (b1 = 0; b1 < scdimen; b1++)
	if (b1 != b)
	  {
	    v1 = extra[b1];
	    PS(v1,
	       if (Len(v1) > b)
	          Fcopy(BodyP(v1)[b],x1);
	       else
	          Fcopy(FZero,x1);,
	       if ((s = FindEnt(v1,b)))
	         Fcopy(s->fac,x1);
	       else
	          Fcopy(FZero,x1););
	    Negate(x1,x1);
	    v1 = adds(v1,v,x1);
	    vfree(extra[b1]);
	    extra[b1] = v1;
	    v1 = adds(sctable[b1][g],sctable[b][g],x1);
	    vfree(sctable[b1][g]);
	    sctable[b1][g] = v1;
	  }
    }
  Fclear(x1);
  for (b = 0; b < scdimen; b++)
    vfree(extra[b]);
  myfree((pointer)extra);
#ifdef CHECKINV
  for (b = 0; b < scdimen; b++)
    {
      vector v2;
      v1 = sctable[b][h];
      v2 = scact(v1,g);
      if (b == 331 && h == 1)
	puts("This time");
      if (Wt(v2) != 1)
	{
	  fprintf(stderr,"Bad inverse for %s at %u\n",gennames[h],b);
	  DIE();
	}
      v1 = sctable[b][g];
      v2 = scact(v1,h);
      if (Wt(v2) != 1)
	{
	  fprintf(stderr,"Bad inverse for %s at %u\n",gennames[g],b);
	  DIE();
	}
    }
#endif
}

#endif

#define readx fread(&x,sizeof(unsigned int),1,f)
#define readxerr(s,v) if (!readx) {fputs(s,stderr);\
				     fputs("Turning scrut off\n",stderr);\
				       scrut = false; return;} v = x

void scsetup()			/* initialize */
{
  FILE *f;
  int l = strlen(sfile);
  unsigned r;
  char *fn = myalloc(l+5,false);
  basiselt b;
  gpgen g;
  int i;
  unsigned int x,ntg;
  bool preims,ims;
  strcpy(fn,sfile);
  strcat(fn,".pb");
  if (!(f = fopen(fn,"rb")))
    {
      fprintf(stderr,"Unable to open %s to read scrut. matrix\n",fn);
      scrut = false;
      myfree(fn);
      return;
    }
  myfree(fn);
  readxerr("Can't read from scrut file",scdimen);
  readxerr("Can't read from scrut file",ntg);
  readxerr("Can't read from scrut file",ims);
  readxerr("Can't read from scrut file",preims);
  if ((!ims && startvecs > 1) || ims != startvecs)
    {
      fprintf(stderr,
	      "Scrut input file does not contain images or contains the wrong number\n");
      fputs("Turning scrut off\n",stderr);
      scrut = false;
      fclose(f);
      return;
    }
  sctabsize = tablesize;
  trans = (vector *)myalloc(sctabsize*sizeof(vector),false);
  for (i=0; i < sctabsize; i++)
    trans[i] = NULL;
  sctable = (vector **)myalloc(scdimen*sizeof(vector *),false);
  for (i = 0; i < scdimen; i++)
    sctable[i] = (vector *)myalloc(nextgen*sizeof(vector),false);
  if (ims)
    {
      for (b=0; b < startvecs; b++)
	if (!(trans[b] = scrdv(f)))
	  {
	    fputs("Can't read images from scrut file\nTurning scrut off\n",stderr);
	    scrut = false;
	    while (b > 0)
	      vfree(trans[--b]);
	    for (i = 0; i < scdimen; i++)
	      myfree((pointer)sctable[i]);
	    myfree((pointer)trans);
	    myfree((pointer)sctable);
	    return;
	  }
    }
  else
    {
      if (scdimen)
	trans[0] = btov(0);
      else
	{
	  newsvec(trans[0],0);
	}
    }
  if (preims)			/* skip over the pre-images in the file */
    for (i = 0; i < scdimen; i++)
      {
	if (!readx) goto SkipErr;
	do {
	  if (!readx) goto SkipErr;
	} while (x);
      }  
  for (g = 0; g < nextgen; g++)
#ifndef INTEGRAL
    if (inverse[g] == NOINVERSE || g <= inverse[g])
      {
#endif
	for (r =0; r < scdimen; r++)
	  if (!( sctable[r][g] = scrdv(f))
#ifdef INTEGRAL
	      || (inverse[g] != NOINVERSE && inverse[g] < g)
#endif
	    )
	{
	  gpgen g1;
	  fputs("Can't read matrices from scrut file\nTurning scrut off\n",stderr);
#ifdef INTEGRAL
	  fputs("Note that invertible non-involutary generators break scrut over Z\n",stderr);
#endif
	  scrut = false;
	  for (b = 0; b < startvecs; b++)
	    vfree(trans[b]);
	  for (g1 =0; g1 < g; g1++)
	    for (b=0; b < scdimen; b++)
	      vfree(sctable[b][g1]);
	  for (b = 0; b < r; b++)
	    vfree(sctable[b][g]);
	  for (b = 0; b < scdimen; b++)
	    myfree((pointer)sctable[b]);
	  myfree((pointer)trans);
	  myfree((pointer)sctable);
	  return;
	}
#ifndef INTEGRAL
     }  
   else
     scinvert(g,inverse[g]);
#endif
  fclose(f);
  return;
  
   SkipErr:
  fprintf(stderr,"Can't read from scrut file while skipping preims\n");
   scrut = false;
    for (b = 0; b < startvecs; b++)
    vfree(trans[b]);
	for (i = 0; i < scdimen; i++)
	  {
	    for (g=0; g < nextgen; g++)
	      vfree(sctable[i][g]);
	    myfree((pointer)sctable[i]);
	  }

  myfree((pointer)trans);
  myfree((pointer)sctable);
  return;

}

void scpk(b1,b2)			/* Note that b2 is being moved to b1 */
     basiselt b1,b2;
{
  if (trans[b1])
    vfree(trans[b1]);
  trans[b1] = trans[b2];
  trans[b2] = NULL;
}


void scdef(b1,g,b2)			/* process the defn b2 = b1.g */
     basiselt b1,b2;
     gpgen g;
{
  basiselt b;
  if (b2 >= sctabsize)
    {
      sctabsize += TABLE_BLOCK;
      trans = (vector *)myrealloc(trans,sctabsize*sizeof(vector),false);
      for (b = sctabsize-TABLE_BLOCK; b < sctabsize; b++)
	trans[b] = NULL;
      }
  else
    if (trans[b2])
      vfree(trans[b2]);
  trans[b2] = scact(trans[b1],g);
}

vector scact(s,g)		/* (real) action of generator on a vector */
     vector s;
     gpgen g;
{
  vector s1,s2;
#ifdef PANIC
  if (Wt(sctable[337][8]) > 1)
    {
      DIE();
    }
#endif
  newsvec(s1,0);
  PS(s,
     fldelt *r = s->body.p;
     basiselt i;
     basiselt l=s->len;
     for (i=0; i<l; i++)
     if (r[i])
     {
       s2 = adds(s1,sctable[i][g],r[i]);
       vfree(s1);
       s1 = s2;
     },
     sent se; sent send;
     send = BodyEndS(s);
     for (se = BodyS(s); se < send; se++)
     {
       s2 = adds(s1,sctable[se->loc][g],se->fac);
       vfree(s1);
       s1 = s2;
     });
  return(s1);
}

vector scv(v)			/* convert a vector to its 
				 translation */
     vector v;
{
  vector s,s1;
  newsvec(s,0);
  PS(v,
      basiselt i; basiselt l = v->len;
      fldelt *r = BodyP(v);
      for (i=0; i<l; i++,r++)
	if (*r)
	  {
	    s1 = adds(s,trans[i],*r);
	    vfree(s);
	    s = s1;
	  },
      sent se; sent send = BodyEndS(v);
      for (se=BodyS(v); se < send; se++)
	{
	  s1 = adds(s,trans[se->loc],se->fac);
	  vfree(s);
	  s = s1;
	});
  return(s);
}

void sccheckc(v)		/* check that b & v are equal */
     vector v;
{
  vector s2=scv(v);
  if (!zerop(s2)) 
    {
      fflush(logfile);
      fprintf(stderr,"Scrut. error : ");
      vprint(stderr,v);
      vprint(stderr,s2);
      DIE();
    }
  vfree(s2);
  return;
}

void sccheckd(b,g,v)		/* check that eb.g = v */
     basiselt b;
     gpgen g;
     vector v;
{
  vector s1,s2,s3,s;
  s1 = trans[b];
  s2 = scact(s1,g);
  s3 = scv(v);
  s = adds(s2,s3,MinusOne);
  if (!zerop(s))
    {
      fflush(logfile);
      fprintf(stderr,"Scrut. error %d.%s = ",b,gennames[g]);
      vprint(stderr,v);
      vprint(stderr,s1);
      vprint(stderr,s2);
      vprint(stderr,s3);
      DIE();
    }
  vfree(s2);
  vfree(s3);
  vfree(s);
  return;
}

void sccheckeq(v1,v2)		/* check that v1 = v2 */
     vector v1,v2;
{
  vector s1 = scv(v1);
  vector s2 = scv(v2);
  vector s = adds(s1,s2,MinusOne);
  if (!zerop(s))
    {
      fflush(logfile);
      fprintf(stderr,"Scrut. error ");
      vprint(stderr,v1);
      fprintf(stderr," = ");
      vprint(stderr,v2);
      vprint(stderr,s1);
      vprint(stderr,s2);
      DIE();
    }
  vfree(s1);
  vfree(s2);
  vfree(s);
}

void sccheckxeq(v1,g,v2)
     vector v1,v2;
     gpgen g;
{
  vector s1,s2,s3,s;
  s1 = scv(v1);
  s2 = scv(v2);
  s3 = scact(s1,g);
  s = adds(s2,s3,MinusOne);
  if (!zerop(s))
    {
      fflush(logfile);
      fprintf(stderr,"Scrut Error %s.",gennames[g]);
      vprint(stderr,v1);
      fprintf(stderr," = ");
      vprint(stderr,v2);
      vprint(stderr,s1);
      vprint(stderr,s2);
      vprint(stderr,s3);
      DIE();
      }
  vfree(s1);
  vfree(s2);
  vfree(s3);
  vfree(s);
  return;
}

void scClean()
{
  basiselt b;
  gpgen g;
  for (b = 0; b < scdimen; b++)
    {
      for (g = 0; g < nextgen; g++)
	vfree(sctable[b][g]);
      myfree((pointer)sctable[b]);
    }
  myfree((pointer)sctable);
  for (b=0; b < sctabsize; b++)
    if (trans[b])
      vfree(trans[b]);
  myfree((pointer)trans);
}
  

#endif
