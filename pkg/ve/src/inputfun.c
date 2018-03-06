/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */
/* C routines to support the input parser */
#include "input.h"
#include <string.h>

gpgen nextgen=0;		/* next free gen. number */
gpgen truegens;
gpgen *inverse;		        /* inverses */
char **gennames;
bool NotGroup = false;		/* true when there is a non-invertible generator */
int WarnLevel;

bool ErrorSpotted = false; /* set by recoverable parse errors */
static  struct _gaw _GawZero;
const gaw GawZero = &_GawZero;

#define GENBLOCK 50


gpgen *newgl()
{
   gpgen *p;
   int i;
   p = (gpgen *)myalloc(sizeof(gpgen)*(nextgen+1),false);
   for (i=0; i <= nextgen; i++)
     p[i] = 0;
   return (p);
}

static void makenotinv(g)
     gpgen g;
{
  int i;
  if (inverse[g] !=  g)
    {
      fprintf(stderr,"Generator %s marking non-invertible and not an involution\n", gennames[g]);
      fprintf(stderr,"treated as non-invertible\n");
    }
  else
    {
      inverse[g] = nextgen;
      inverse[nextgen] = g;
      i = strlen(gennames[g]);
      if (!(nextgen%GENBLOCK))
	gennames = (char **)myrealloc(gennames,
				      sizeof(char *)*(nextgen+GENBLOCK),
				      false);
      gennames[nextgen] = (char *)myalloc(i+5,false);
      strcpy(gennames[nextgen],gennames[g]);
      strcpy(gennames[nextgen]+i,"^-1");
    }
}

void setupinverses( noinv, notinv)
     gpgen *noinv;
     gpgen *notinv;
{
  int i;
  gpgen g;
  gpgen *s;
  truegens = nextgen;
  if (noinv)
    {
      for (i=0; i < nextgen; i++)
	inverse[i] = i;
      for (s = noinv; *s; s++)
	{
	  NotGroup = true;
	  inverse[*s-1] = NOINVERSE;
	}
      myfree((pointer)noinv);
    }
  else /* NULL here implies that there was a '*' in the presentation */
    {
      for (i =0; i < nextgen; i++)
	inverse[i] = NOINVERSE;
      NotGroup = true;
    }
  if (notinv)
    {
      for (s = notinv; *s; s++,nextgen++)
	makenotinv(*s-1);
      myfree((pointer)notinv);
    }
  else
    for (g=0; g < truegens; g++,nextgen++)
      makenotinv(g);
}	

void addugen(l)
     char l;
{
  if (!(nextgen%GENBLOCK))
    if (nextgen)
      gennames = (char **)myrealloc(gennames,
				    sizeof(char *)*(nextgen+GENBLOCK),
				    false);
    else
      gennames = (char **)myalloc(sizeof(char *)*GENBLOCK,false);
  gennames[nextgen] = (char *)myalloc(2,false);
  gennames[nextgen][0] = l;
  gennames[nextgen][1] = '\0';
  nextgen++;
}

void addlgen(s)
     char *s;
{
  if (!(nextgen%GENBLOCK))
    if (nextgen)
      gennames = (char **)myrealloc(gennames,
				    sizeof(char *)*(nextgen+GENBLOCK),
				    false);
    else
      gennames = (char **)myalloc(sizeof(char *)*GENBLOCK,
				  false);
  gennames[nextgen] = s;
  nextgen++;
}

				
gpwd ginverse(w)		/* invert a gp word */
     gpwd w;
{
  gpwd w1;
  int i;
  newgpwd(w1,w->len);
  for (i=0; i<w->len; i++)
    if ((w1->word[w->len-1-i] = inverse[w->word[i]])==NOINVERSE)
    {
    	fprintf(stderr,"Can't invert generator %c\n",w->word[i]);
    	exit(2);
    }
  return(w1);
}

gpwd gpower(w,n)		/* power one */
     gpwd w;
     int n;
{
  gpwd w1;
  int i,j,k;
  newgpwd(w1,n*w->len);
  k = 0;
  for (i = 0; i < n; i++)
    for (j=0; j<w->len; j++)
      w1->word[k++] = w->word[j];
  wfree(w);
  return(w1);
}

gpwd gconcat(w1,w2)		/* Concatenate two */
     gpwd w1,w2;
{
  gpwd w;
  int i,j;
  newgpwd(w, w1->len+w2->len);
				/* first one */
  for (i=0; i<w1->len; i++)
    w->word[i] = w1->word[i];
  j = 0;
				/* do any cancellation */
  while (i > 0 && w2->word[j] == inverse[w->word[i-1]] && j < w2->len)
    {
      j++;
      i--;
      w->len -= 2;
    }
				/* rest of the second one */
  for (; j<w2->len; j++)
    w->word[i++] = w2->word[j];
  return(w);
}

gpwd gconj(w1,w2)
     gpwd w1,w2;
{
  gpwd w,w3,w4;
  w3 = gconcat(w1,w2);
  w4 = ginverse(w2);
  w = gconcat(w4,w3);
  wfree(w1);
  wfree(w2);
  wfree(w3);
  wfree(w4);
  return w;
}

gpwd gcomm(w1,w2)
     gpwd w1,w2;
{
  gpwd w,w3,w4,w5,w6;
  w3 = gconcat(w1,w2);
  w4 = ginverse(w1);
  w5 = ginverse(w2);
  w6 = gconcat(w4,w5);
  w = gconcat(w6,w3);
  wfree(w1);
  wfree(w2);
  wfree(w3);
  wfree(w4);
  wfree(w5);
  wfree(w6);
  return w;
}


#define AlgRelWt 3		/* Standard weight for algebra relators */
#define AlgRelLwt 3		/* lookahead weight for same */
#define InvRelWt 3		/* Standard weight for Gg relators */
#define InvRelLwt 6		/* lookahead weight for same */
#define GpRelDiv 2		/* Gp rels get weight len/this */
#define AbWt 2
#define AbLWt 2


rell neweq(lft,rt)
     gaw lft,rt;
{
  rell r;
  r = (rell)myalloc(sizeof(struct _rell),false);
  r->type = algebra;
  r->entry.a.l = lft;
  r->entry.a.r = rt;
  r->lwt = AlgRelLwt;
  r->wt = AlgRelWt;
  r->next = NULL;
  r->NexttoDo = 0;
  return(r);
}

rell gtorel(g)
     gpwd g;
{
  rell r;
  int i;
  r = (rell)myalloc(sizeof(struct _rell),false);
  r->type = group;
  r->entry.g = g;
  for (i=0; i < g->len; i++)
     if (inverse[g->word[i]] == NOINVERSE)
     {
  	fprintf(stderr,"Non invertible generator in group relator\n");
	DIE();
     }
  r->lwt = r->wt = (g->len+GpRelDiv-1)/GpRelDiv;
  r->next = NULL;
  r->NexttoDo = 0;
  return(r);
}

rell relconcat(r1,r2)		/* in situ concatenate two relator lists */
     rell r1,r2;
{
  rell r;
  if (!r1)
    return(r2);
  for (r = r1; r->next; r = r->next)
    ;
  r->next = r2;
  return(r1);
}

void cleanrels(rp)
     rell *rp;
{
  rell r;
  while ((r = *rp))
    {
      if (r->type == group && r->entry.g->len == 0)
	{
	  *rp = r->next;
	  myfree(r->entry.g);
	  myfree(r);
	}
      else
	rp = &(r->next);
    }
}


void getinput()			/* drive the parser */
{
  gpgen g,g1;
  gaw w,w1;
  gpwd gw;
  rell r;
  GawZero->type = scalar;
  Finit(GawZero->body.scalar);
  Fcopy(FZero,GawZero->body.scalar);
  yyin = ifp;
  yyparse();			/* let Bison do its work */
  if (ErrorSpotted)
    {
      fprintf(stderr,"Error in presentation file\n");
      DIE();
    }
  cleanrels(&relators);
  if (subgens->type == rels)
    cleanrels(&(subgens->list.r));
#ifdef INTEGRAL
  GenWt = (weight *)myalloc(sizeof(weight)*nextgen,false);
  GenLWt = (weight *)myalloc(sizeof(weight)*nextgen,false);
  NextToDo = (lattpos *)myalloc(sizeof(lattpos)*nextgen,false);
  SavedNext = (lattpos *)myalloc(sizeof(lattpos)*nextgen,false);
#endif
  for (g = 0; g < nextgen; g++)	/* put in Gg relations */
    {
      if (inverse[g] != NOINVERSE)
	{
	  newgpwd(gw,2);
	  gw->word[0] = g;
	  gw->word[1] = inverse[g];
	  r = gtorel(gw);
	  r->next = relators;
	  r->wt = InvRelWt;
	  r->lwt = InvRelLwt;
	  relators = r;
	}
      else if (!Abelian && !SuppressRelators)
	{
	  newgpwd(gw,1);
	  gw->word[0] = g;
	  w = newgaw;
	  w->type = grp;
	  w->body.group = gw;
	  r = neweq(w,w);
	  r->wt = InvRelWt;
	  r->lwt = InvRelLwt;
	  r->next = relators;
	  relators = r;
	}
      if (Abelian)
	{
	  for (g1 = g+1; g1 < nextgen; g1++)
	    {
	      newgpwd(gw,2);
	      gw->word[0] = g;
	      gw->word[1] = g1;
	      w = newgaw;
	      w->type = grp;
	      w->body.group = gw;
	      newgpwd(gw,2);
	      gw->word[0] = g1;
	      gw->word[1] = g;
	      w1 = newgaw;
	      w1->type = grp;
	      w1->body.group = gw;
	      r = neweq(w,w1);
	      r->wt = AbWt;
	      r->lwt = AbLWt;
	      r->next = relators;
	      relators = r;
	    }
	}
#ifdef INTEGRAL
      GenWt[g] = STD_GWT;
      GenLWt[g] = STD_GLWT;
      NextToDo[g] = 0;
      SavedNext[g] = 0;
#endif
    }
}

int yyerror(s)			/* minimal error handler */
     char *s;
{
  ErrorSpotted = true;
  return true;
}

swl newswl(loc,w,next)
     basiselt loc;
     gaw w;
     swl next;
{
  swl s;
  s = (swl) myalloc(sizeof(_swl),false);
  s->loc = loc;
  s->w = w;
  s->next = next;
  return s;
}

wordl newwordl(w)
	gaw w;
{	
    wordl wl;
    wl = (wordl)myalloc(sizeof(struct _wordl),false);
    if (w)
      {
	wl->type = wl_packed;
	wl->body.p.wds = (gaw *) myalloc(sizeof(gaw)*startvecs,false);
	wl->body.p.wds[0] = w;
	wl->body.p.len = 1;
      }
    else
      {
	wl->type = wl_sparse;
	wl->body.s = NULL;
      }
    return(wl);
}	

smgl newsmgl(lft,rt)
     wordl lft,rt;
{
  smgl s;
  s = (smgl) myalloc(sizeof(struct _smgl),false);
  s->type = normal;
  s->body.norm.lhs = lft;
  s->body.norm.rhs = rt;
  s->next = NULL;
  return(s);
}

smgl newuniv(start,end,w)
     basiselt start,end;
     gaw w;
{
  smgl s;
  s = (smgl) myalloc(sizeof(struct _smgl),false);
  s->type = universal;
  s->body.univ.start = start;
  s->body.univ.end = end;
  s->body.univ.w = w;
  s->next = NULL;
  return s;
}

#define newsmg(s) s = (smg) myalloc(sizeof(struct _smg),false);

smg relltosmg(r)
     rell r;
{
  smg s;
  newsmg(s);
  s->type = rels;
  s->list.r = r;
  return(s);
}

smg smgltosmg(sl)
     smgl sl;
{
  smg s;
  newsmg(s);
  s->type = smgens;
  s->list.s = sl;
  return(s);
}


void setfield(p)
     unsigned int p;
{
#ifdef GFP
  unsigned int top,bottom,mid;
  fldelt pm;
  if (!ForceCharacteristic)
    {
      if (p)
	{
	  if (!IsPrime(p))
	    Warn("Characteristic not prime\n",true);
	  P = p;
	}
      else
	{
	  Warn("Use qme or zme for characteristic 0\n",true);
	}
    }
#else
  if (!ForceCharacteristic && p)
    {
      fprintf(stderr,"Use me for positive characteristic %d\n", (int)p);
      ErrorSpotted = true;
    }
#endif
}


#ifdef GFP

int inttorat(x)
     int x;
{
  int r=x;
  while (r < 0)
    r += P;
  return r;
}

int makerat(x,y)
     int x,y;
{
  fldelt z,xx,yy;
  xx = rattofelt(inttorat(x));
  yy = rattofelt(inttorat(y));
  z = fldinv(yy);
  Fmul(z,xx,z);
  return (int)z;
}

#endif

#if defined(RATIONAL) || defined(INTEGRAL)

MP_INT startint(x)
     int x;
{
  MP_INT t;
  mpz_init_set_si(&t,x);
  return t;
}

MP_INT contint(t,x)
     MP_INT t;
     int x;
{
  mpz_mul_ui(&t,&t,(unsigned)10);
  mpz_add_ui(&t,&t,(unsigned)x);
  return t;
}

MP_INT negint(t)
     MP_INT t;
{
  mpz_neg(&t,&t);
  return t;
}

#endif

#ifdef RATIONAL
MP_RAT inttorat(t)
     MP_INT t;
{
  MP_RAT q;
  mpq_init(&q);
  mpq_set_ui(&q,0,1);
  mpq_set_num(&q,&t);
  mpz_clear(&t);
  return q;
}

MP_RAT makerat(x,y)
     MP_INT x;
     MP_INT y;
{
  MP_RAT q;
  mpq_init(&q);
  mpq_set_num(&q,&x);
  mpq_set_den(&q,&y);
  mpz_clear(&x);
  mpz_clear(&y);
  return q;
}

#endif


#ifdef INTEGRAL

MP_INT makerat(x,y)
     MP_INT x;
     MP_INT y;
{
  MP_INT q;
  MP_INT r;
  mpz_init(&q);
  mpz_init(&r);
  mpz_divmod(&q,&r,&x,&y);
  if (mpz_cmp_ui(&r,0))
    Warn("Inexact division",true);
  mpz_clear(&r);
  mpz_clear(&x);
  mpz_clear(&y);
  return q;
}

#endif
