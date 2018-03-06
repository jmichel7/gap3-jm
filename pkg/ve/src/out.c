/* This program is Copyright (C) Stephen Linton 1991 and 1992 
   See the file copying for details */

#include "me.h"
#include <ctype.h>
#include <string.h>
/* Various output routines */

char *ofn;		/* prefix for matrix output */
bool cayfmt;		/* output selectors */
#ifdef GFP
bool mtxfmt;
#endif
bool gapfmt;
bool gapfmt2;
bool gapfmt3;
bool mcrfmt;
bool mcbfmt;
bool axiomfmt;
bool PrintImages;
char *gap2var;


#ifdef GFP
void outfelt(f,x)
     FILE *f;
     fldelt x;
{
  fprintf(f,"%u",(unsigned int)x);
}

void outbin(f,x)
     FILE *f;
     fldelt x;
{
    fwrite(&x,sizeof(fldelt),1,f);
}

char *felttos(s,x)
     char *s;
     fldelt x;
{
  sprintf(s,"%u",(unsigned int)x);
  while (*s)
    s++;
  return(s);
}

#endif

#ifdef RATIONAL

void outfelt(f,q)
     FILE *f;
     fldelt q;
{
  MP_INT t;
  mpz_init(&t);
  mpq_get_num(&t,&q);
  mpz_out_str(f,10,&t);
  if (mpz_cmp_ui(&t,0))
    {
      mpq_get_den(&t,&q);
      if (mpz_cmp_ui(&t,1))
	{
	  fputc('/',f);
	  mpz_out_str(f,10,&t);
	}
    }
  mpz_clear(&t);
}

char * felttos(s,q)
     char *s;
     fldelt q;
{
  MP_INT t;
  mpz_init(&t);
  mpq_get_num(&t,&q);
  mpz_get_str(s,10,&t);
  while (*s)
    s++;
  if (mpz_cmp_ui(&t,0))
    {
      mpq_get_den(&t,&q);
      if (mpz_cmp_ui(&t,1))
	{
	  *s++ = '/';
	  mpz_get_str(s,10,&t);
	  while (*s)
	    s++;
	}
    }
  mpz_clear(&t);
  return(s);
}

void outbin(f,x)
     FILE *f;
     fldelt x;
{
     mpz_out_raw(f,mpq_numref(&x));
     mpz_out_raw(f,mpq_denref(&x));
}

#endif

#ifdef INTEGRAL

void outfelt(f,z)
     FILE *f;
     fldelt z;
{
  mpz_out_str(f,10,&z);
}

void outbin(f,z)
     FILE *f;
     fldelt z;
{
  mpz_out_raw(f,&z);
}

char *felttos(s,z)
     char *s;
     fldelt z;
{
  mpz_get_str(s,10,&z);
  while (*s)
    s++;
  return s;
}

#endif

static void CayHead(f)
     FILE *f;
{
  char *s = ofn;
  while (*++s)		/* find end of file name */
    ;
  while (*s != '/' && s >= ofn) /* search for start of last component */
    s--;
  s++;
  fprintf(f,"LIBRARY %s;\n",s);
#ifdef GFP
  fprintf(f,"fld : field(%u);\n",P);
#endif
#ifdef RATIONAL
  fprintf(f,"fld : rationals;\n");
#endif
#ifdef INTEGRAL
  fprintf(f,"fld : integers;\n");
#endif
  fprintf(f,"vs : vector space (%u,fld);\n",nalive);
}

static void CayStMat(f,g)
     FILE *f;
     gpgen g;
{
#ifdef INTEGRAL
  if (g == (gpgen) -2)
    fprintf(f,"Lattice := MAT(\n\t\t");
  else
#endif
    if (g != (gpgen) -1)
      fprintf(f,"\t%s = MAT(\n\t\t",gennames[g]);
    else
      fprintf(f,"Images := MAT(\n\t\t");
}

static void CayStGens(f)
     FILE *f;
{
  fputs("grp : matrix group (vs);\ngrp.generators :\n",f);
}

static void CayEGens(f)
     FILE *f;
{
  fputs(";\n",f);
}

static void CayTail(f)
     FILE *f;
{
  fprintf(f,"VEruntime = %lu;\n",RunTime());
  fputs("FINISH;",f);
}

static void GAPHead(f)
     FILE *f;
{
#ifdef GFP
  fprintf(f,"field := GF(%u);\n",(unsigned int)P);
#endif
#ifdef RATIONAL
  fputs("field := Rationals;\n",f);
#endif
#ifdef INTEGRAL
  fputs("field := Integers;\n",f);
#endif
  fprintf(f,"VErunTime := %lu;\n",RunTime());
}

static void GAPStMat(f,g)
     FILE *f;
     gpgen g;
{
#ifdef INTEGRAL
  if (g == (gpgen)-2)
    fputs("lattice_mat",f);
  else
#endif
    if (g == (gpgen)-1)
      fputs("images_mat",f);
    else
      fprintf(f,"gen_%s",gennames[g]);
  fputs(" := field.one*[",f);
}

static void GAPEGens(f)
     FILE *f;
{
  bool started = false;
  gpgen g;
  fputs("gens := [",f);
  for (g=0; g < nextgen; g++)
    if (inverse[g] == NOINVERSE || g <= inverse[g])
      {
	if (started)
	  fputc(',',f);
	fprintf(f,"gen_%s",gennames[g]);
	started = true;
      }
  fputs("];\n",f);
}

static void MCRHead(f)
     FILE *f;
{
  fprintf(f,
	  "%u %u %u %u",
	  nalive,
	  truegens,
	  PrintImages ? (unsigned int)startvecs : 0,
	  (unsigned int)TrackPreIm);
#ifdef INTEGRAL
  fprintf(f," %u",nextlat);
#endif
  fputc('\n',f);
}

static void MCRTail(f)
     FILE *f;
{
  fprintf(f,"%lu\n", RunTime());
}

static void MCBHead(f)
     FILE *f;
{
  unsigned int t;
  t = (unsigned int)nalive;
  fwrite(&t,sizeof(unsigned int),1,f);
  t = (unsigned int) truegens;
  fwrite(&t,sizeof(unsigned int),1,f);
  t = PrintImages ? (unsigned int)startvecs : 0;
  fwrite(&t,sizeof(unsigned int),1,f);
  t = (unsigned int)TrackPreIm;
  fwrite(&t,sizeof(unsigned int),1,f);
#ifdef INTEGRAL
  t = (unsigned int)nextlat;
  fwrite(&t,sizeof(unsigned int),1,f);
#endif
}

static void MCBTail(f)
     FILE *f;
{
  unsigned long t;
  t = RunTime();
  fwrite(&t,sizeof(unsigned long),1,f);
}


static void AXHead(f)
     FILE *f;
{
#ifdef GFP
  fprintf(f,"fld := FF(%u,1)\n",(unsigned int) P);
#endif
#ifdef RATIONAL
  fputs("fld := FRAC INT\n",f);
#endif
#ifdef INTEGRAL
  fputs("fld := INT\n",f);
#endif
  fprintf(f,"VErunTime := %lu\n",RunTime());  
}

static void AXStMat(f,g)
     FILE *f;
     gpgen g;
{
#ifdef INTEGRAL
  if (g == (gpgen)-2)
    fputs("lattice_mat",f);
  else
#endif
    if (g == (gpgen)-1)
      fputs("images_mat",f);
    else
      fprintf(f,"gen_%s",gennames[g]);
  fputs(" := matrix([",f);
}

static void AXEGens(f)
     FILE *f;
{
  bool started = false;
  gpgen g;
  fputs("gens := [",f);
  for (g=0; g < nextgen; g++)
    if (inverse[g] == NOINVERSE || g <= inverse[g])
      {
	if (started)
	  fputc(',',f);
	fprintf(f,"gen_%s",gennames[g]);
	started = true;
      }
  fputs("]\n",f);
}

void GAPStPims(f)
     FILE *f;
{
  gpgen g;
  for (g = 0; g < nextgen; g++)
    if (inverse[g] == NOINVERSE || inverse[g] >= g)
      fprintf(f,"g%s := AbstractGenerator(\"%s\");\n",
	      gennames[g],gennames[g]);
  fputs("preImages := [\n",f);
}

void GAPStPim(f,fmg)
     FILE *f;
     basiselt fmg;
{
  fprintf(f,"rec( modGen := %u, word := ",fmg);
}

void GAPWrGen(f,g)
     FILE *f;
     gpgen g;
{
  fprintf(f,"g%s",gennames[g]);
}

void GAPEPim(f)
     FILE *f;
{
  fputs(")",f);
}

void MCRStPim(f,fmg)
     FILE *f;
     basiselt fmg;
{
  fprintf(f,"%u ",(unsigned int)fmg);
}

void MCRWrGen(f,g)
     FILE *f;
     gpgen g;
{
  if (inverse[g] == NOINVERSE || inverse[g] >= g)
    fprintf(f,"%d",(int)g+1);
  else
    fprintf(f,"%d",-((int)inverse[g]+1));
}

void MCREPim(f)
     FILE *f;
{
  fputs(" 0",f);
}

void MCBStPim(f,fmg)
     FILE *f;
     basiselt fmg;
{
  unsigned int t = (unsigned int)fmg;
  fwrite(&t,sizeof(unsigned int),1,f);
}

void MCBWrGen(f,g)
     FILE *f;
     gpgen g;
{
  int t;
  if (inverse[g] == NOINVERSE || inverse[g] >= g)
    t = (int)g+1;
  else
    t = -((int)inverse[g]+1);
  fwrite(&t,sizeof(unsigned int),1,f);
}

void MCBEPim(f)
     FILE *f;
{
  unsigned int t = 0;
  fwrite(&t,sizeof(unsigned int),1,f);
}

void AXStPims(f)
     FILE *f;
{
  gpgen g;
  bool started = false;
  fputs("rec := Record(modGen : PositiveInteger, word : Polynomial Integer)\n",f);
  fputs("gen_names : LIST Symbol := [",f);
  for (g=0; g < nextgen; g++)
    if (inverse[g] == NOINVERSE || inverse[g] >= g)
      {
	if (started)
	  fputc(',',f);
	fputs(gennames[g],f);
	started = true;
      }
  fputs("]\npreImagesL : LIST LIST Any := [_\n",f);
}

void AXStPim(f,fmg)
     FILE *f;
     basiselt fmg;
{
  fprintf(f,"[ %u, ",(unsigned int)fmg);
}

void AXWrGen(f,g)
     FILE *f;
     gpgen g;
{
  if (inverse[g] == NOINVERSE || inverse[g] >= g)
    fputs(gennames[g],f);
  else
    fprintf(f,"%s ** -1",gennames[inverse[g]]);
}

void AXEPim(f)
     FILE *f;
{
  fputs("]",f);  
}

static void GAP2Head(f)
     FILE *f;
{
  fputs(gap2var,f);
  fputs(" := rec(\n",f);
#ifdef GFP
  fprintf(f,"field := GF(%u),\n",(unsigned int)P);
#endif
#ifdef RATIONAL
  fputs("field := Rationals,\n",f);
#endif
#ifdef INTEGRAL
  fputs("field := Integers,\n",f);
#endif
  fprintf(f,"dim := %u,\n",nalive);
  fprintf(f,"statistics := rec( time := %lu/1000 ),\n",RunTime());
}


static void GAP2EGens(f)
     FILE *f;
{
  bool started = false;
  gpgen g;
  fputs("gens := [",f);
  for (g=0; g < nextgen; g++)
    if (inverse[g] == NOINVERSE || g <= inverse[g])
      {
	if (started)
	  fputc(',',f);
	fprintf(f,"~.gen_%s",gennames[g]);
	started = true;
      }
  fputs("],\n",f);
}

void GAP2StPims(f)
     FILE *f;
{
  gpgen g;
  for (g = 0; g < nextgen; g++)
    if (inverse[g] == NOINVERSE || inverse[g] >= g)
      fprintf(f,"g%s := AbstractGenerator(\"%s\"),\n",
	      gennames[g],gennames[g]);
  fputs("preImages := [\n",f);
}

static void GAP2Tail(f)
     FILE *f;
{
  fputs("isNMEOutput := true);\n",f);
}

void GAP2WrGen(f,g)
     FILE *f;
     gpgen g;
{
  fprintf(f,"~.g%s",gennames[g]);
}


static void GAP2StMat(f,g)
     FILE *f;
     gpgen g;
{
#ifdef INTEGRAL
  if (g == (gpgen)-2)
    fputs("lattice_mat",f);
  else
#endif
  if (g == (gpgen)-1)
    fputs("images_mat",f);
  else
    fprintf(f,"gen_%s",gennames[g]);
  fputs(" := ~.field.one*[",f);
}


outstyle outstyles[] = {
  {"","w",outfelt, FuncPt CayHead,	/* CAYLEY */
     ",","","",":\n\t\t",FuncPt CayStMat,")",
     ",\n",FuncPt CayStGens, FuncPt CayEGens,FuncPt CayTail,";\n",
   false},
  {".g","w",outfelt,FuncPt GAPHead,	/* GAP */
     ",","[","]",",\n",FuncPt GAPStMat,"];\n",
     "",NULL,FuncPt GAPEGens, NULL,"",
     true, FuncPt GAPStPims, FuncPt GAPStPim, FuncPt GAPWrGen,"*","IdWord",FuncPt GAPEPim,
     ",\n","];\n"},
  {".pa","w",outfelt,FuncPt MCRHead,	/* Machine readable - ASCII */
     " ","","","\n",NULL,"\n",
     "",NULL,NULL, FuncPt MCRTail,"",
   true, NULL, FuncPt MCRStPim, FuncPt MCRWrGen, " ", "", FuncPt MCREPim, "\n", "\n" },
  {".pb","wb",outbin,FuncPt MCBHead,	/* Machine readable - binary */
     "","","","",NULL,"",
     "",NULL,NULL, FuncPt MCBTail,"",
   true, NULL, FuncPt MCBStPim, FuncPt MCBWrGen, "", "", FuncPt MCBEPim, "",""},
  {".input","w",outfelt,FuncPt AXHead,	/* AXIOM */
     ",","[","]",",_\n",FuncPt AXStMat,"]) :: MATRIX fld\n",
     "",NULL,FuncPt AXEGens, NULL,"",
   true,FuncPt AXStPims, FuncPt AXStPim, FuncPt AXWrGen, "*", "1", FuncPt AXEPim,",_\n",
     "]\npreImages := preImagesL :: LIST rec\n" },
  {".g","w",outfelt,FuncPt GAP2Head,	/* GAP */
     ",","[","]",",\n",FuncPt GAP2StMat,"],\n",
     "",NULL,FuncPt GAP2EGens, FuncPt GAP2Tail,"",
     true, FuncPt GAP2StPims, FuncPt GAPStPim, FuncPt GAP2WrGen,"*","IdWord",FuncPt GAPEPim,
     ",\n","],\n"},

};


void OutRow(f,v,s)
     FILE *f;
     vector v;
     outstyle *s;
{
  basiselt i,bb;
  sent se, send;
#ifdef PACKED_VEC
  fldelt *r;
  basiselt l;
#endif
  bool started=false;
  fldelt z;
  z = FZero;
  fputs(s->stvec,f);
  PS(v,
     l = Len(v);
     r = BodyP(v);
     for (i=0; i<nextbe; i++)
     {
       if (!rowdel(i))
	 {
	   if (started)
	     fputs(s->sep,f);
	   ASSERT(s->wrfelt);
	   (*s->wrfelt)(f,i < l ? *r : z);
	   started = true;
	 }
       r++;
     }
     ,
     {
       send = BodyEndS(v);
       se = BodyS(v);
       if (se < send)
	 bb = se->loc;
       for (i = 0; i< nextbe; i++)
	 if (se < send && i==bb)
	   {
	     ASSERT(!rowdel(i));
	     if (started)
	       fputs(s->sep,f);
	     ASSERT(s->wrfelt);
	     (s->wrfelt)(f,se->fac);
	     
	     started = true;
	     if (++se < send)
	       bb = se->loc;
	   }
	 else {
	   if (!rowdel(i))
	     {
	       if (started)
		 fputs(s->sep,f);
	       ASSERT(s->wrfelt);
	       (s->wrfelt)(f,z);
	       started = true;
	     }
	 }});
  fputs(s->evec,f);
}


static void DoOutput(m)
     OutMode m;
{
  char *fn;
  outstyle *s = outstyles+m;
  FILE *f;
  basiselt b;
  vector v1,v;
  gpgen g;
  int i;
  bool started,gstarted;
  fn = myalloc(sizeof(char)+(strlen(ofn)+strlen(s->suffix)+1),false);
  strcpy(fn,ofn);
  strcat(fn,s->suffix);
  f = fopen(fn,s->ft);
  if (!f)
    {
      fprintf(stderr,"Can't open output file %s in mode %s\n",fn,s->ft);
      myfree(fn);
      return;
    }
  myfree(fn);
  if (s->wrhead)
    (s->wrhead)(f); /* First do the header according to the style */
  if (PrintImages)
    {
      if (s->stmat)
	(s->stmat)(f,(gpgen)-1); /* Now set up to write images */
      for (b=0; b < startvecs; b++)
      	{
	  if (b) fputs(s->vsep,f);
	  v = btov(b);
	  v1 = vroot(v);
	  if (v1 != v) vfree(v);
	  OutRow(f,v1,s);
	  vfree(v1);
	}
      fputs(s->emat,f);
      fputs(s->eims,f);
    }
  if (TrackPreIm && s->PreImsPossible)
    {
      if (s->stpims)
	(s->stpims)(f);
      started = false;
      for (b=0; b < nextbe; b++)
	if (!rowdel(b))
	  {
	    if (started)
	      fputs(s->pimsep,f);
	    started = true;
	    if (s->stpim)
	      (s->stpim)(f,table[b].fmgen+1);
	    for (i=0; i < table[b].wlen; i++)
	      {
		if (i)
		  fputs(s->gensep,f);
		(s->wrgen)(f,table[b].word[i]);
	      }
	    if (!i)
	      fputs(s->nullword,f);
	    if (s->epim)
	      (s->epim)(f);
	  }
      fputs(s->epims,f);
    }
  if (s->stgens)
    (s->stgens)(f);
  gstarted = false;
  for (g=0; g < nextgen; g++)
    if (inverse[g] == NOINVERSE || inverse[g] >= g)
      {
	if (gstarted)
	  fputs(s->matsep,f);
	gstarted = true;
	if (s->stmat)
	  (s->stmat)(f,g);
	started = false;
	for (b = 0; b < nextbe; b++)
	  if (!rowdel(b))
	    {
	      if (started)
		fputs(s->vsep,f);
	      started = true;
	      v = tabent(b,g);
	      ASSERT(v);
	      OutRow(f,v,s);
	    }
	fputs(s->emat,f);
      }
  if (s->egens)
    (s->egens)(f);
#ifdef INTEGRAL
  OutLatt(f,m);
#endif
  if (s->wrtail)
    (s->wrtail)(f);
  fclose(f);
}

#ifdef GFP
static void MTXoutrow(f,v,mode,len,fmt)
     FILE *f;
     vector v;
     int mode,len;
     char *fmt;
{
  basiselt l,i,bb;
  sent se, send;
  fldelt *r;
  int skipped=0;
  bool started=false;
  PS(v,
     l = Len(v);
     r = BodyP(v);
     for (i=0; i<nextbe; i++)
     {
       if (!rowdel(i))
	 {
	   fprintf(f,fmt,(i<l) ? *r : 0);
	   if ((i-skipped)%len == len-1) fputc('\n',f);
	   started = true;
	 }
       else
	 skipped++;
       r++;
     }
     ,
     {
       send = BodyEndS(v);
       se = BodyS(v);
       if (se < send)
	 bb = se->loc;
       for (i = 0; i< nextbe; i++)
	 if (se < send && i==bb)
	   {
	     ASSERT(!rowdel(i));
	     fprintf(f,fmt,se->fac);
	     if ((i-skipped)%len == len-1)
	       fputc('\n',f);
	     started = true;
	     if (++se < send)
	       bb = se->loc;
	   }
	 else {
	   if (!rowdel(i))
	     {
	       fprintf(f,fmt,0);
	       if ((i-skipped)%len == len-1)
		 fputc('\n',f);
	       started = true;
	     }
	   else
	     skipped++;
	 }});
  if (nalive%len) fputc('\n',f);
}


static void DoMTXOutput()
{
  char *fn,*fmt;
  int mode,len;
  basiselt b;
  gpgen g;
  FILE *f;
  vector v,v1;
  if (P < 10)
    {
      mode = 1;
      len = 80;
      fmt = "%1u";
    }
  else if (P<100)
    {
      mode = 3;
      len = 20;
      fmt = "%3u";
    }
  else
    {
      mode = 4;
      len = 5;
      fmt = "%8u";
    }
  for (g=0; g< nextgen; g++)
    if (inverse[g] == NOINVERSE || inverse[g] >= g)
      {
	fn = myalloc(strlen(ofn)+3+strlen(gennames[g]),false);
	strcpy(fn,ofn);
	strcat(fn,".");
	strcat(fn,gennames[g]);
	f = fopen(fn,"w");
	if (!f)
	  {
	    fprintf(stderr,"Cannot open output file %s\n",fn);
	    exit(12);
	  }
	fprintf(f,
		"%2u%6u%6u%6u\n",
		(unsigned int)mode,
		(unsigned int)P,
		(unsigned int)nalive,
		(unsigned int)nalive);
	for (b=0; b < nextbe; b++)
	  if (!rowdel(b))
	    MTXoutrow(f,tabent(b,g),mode,len,fmt);
	fclose(f);
	myfree(fn);
      }
  if (PrintImages)
    {
      fn = myalloc(strlen(ofn)+4,false);
      strcpy(fn,ofn);
      strcat(fn,".IM");
      f = fopen(fn,"w");
      if (!f)
	{
	  fprintf(stderr,"Cannot open output file %s\n",fn);
	  exit(12);
	}
      fprintf(f,
	      "%2u%6u%6u%6u\n",
	      (unsigned int)mode,
	      (unsigned int)P,
	      (unsigned int)startvecs,
	      (unsigned int)nalive);
      for (b=0; b < startvecs; b++)
      	{
	  v = btov(b);
	  v1 = vroot(v);
	  if (v1 != v) vfree(v);
	  MTXoutrow(f,v1,mode,len,fmt);
	  vfree(v1);
	}
      fclose(f);
      myfree(fn);
    }	
}
#endif

void DoLimittedOutput()
{
  char *fn;
  FILE *f;
  fn = (char *)myalloc(sizeof(char)*(strlen(ofn)+3),false);
  strcpy(fn,ofn);
  strcat(fn,".g");
  f = fopen(fn,"w");
  if (!f)
    {
      fprintf(stderr,"Can't open file %s for writing\n",fn);
      myfree(fn);
      return;
    }
  myfree(fn);
  fprintf(f,"%s := rec( dim := %u,\n",gap2var,nalive);
  fprintf(f,"statistics := rec( time := %lu/1000 )\n",RunTime());
  fputs(");\n",f);
  fclose(f);
  return;
}

void mxsout()
{
  if (cayfmt)
    DoOutput(Cayley);
  if (gapfmt)
    DoOutput(GAP);
  if (gapfmt2)
    DoOutput(GAP2);
  if (gapfmt3)
    DoLimittedOutput();
  if (mcrfmt)
    DoOutput(MCR);
  if (mcbfmt)
    DoOutput(MCB);
  if (axiomfmt)
    DoOutput(AXIOM);
#ifdef GFP
  if (mtxfmt) /* have to do this one specially as it all happens in different files */
    DoMTXOutput();
#endif
}


#ifdef LOGGING			/* don't need these unless loggging  */

void printgaw(f,w)			/* print group algebra word w on file f*/
     FILE *f;
     gaw w;
{
  char s[1000];
  gawtos(w,s,false);
  fputs(s,f);
}

void gpwdtos(s,g)			/* put the string for g into buffer at s */
     char *s;
     gpwd g;
{
  char *s1;
  int  i,l = g->len;
  *s = '\0';
  for (i=0; i<l; i++)
    {
      s1 = gennames[g->word[i]];
      if (i && islower(*s1))
	*s++ = '*';
      strcpy(s,s1);
      while (*++s)
	;
    }
  return;
}

void gawtos(w,s,f)			/* put the string for w into buffer s
				 if f then an addition needs brackets */
     gaw w;
     char *s;
     bool f;
{
  switch(w->type)
    {
    case scalar:
      s = felttos(s,w->body.scalar);
      return;
    case grp:
      gpwdtos(s,w->body.group);
      return;
    case sum:
      if (f) *s++ = '(';
      gawtos(w->body.sum.a,s,false);
      while (*s++);
      *(s-1) = '+';
      gawtos(w->body.sum.b,s,false);
      if (f)
	{
	  while (*s++);
	  *(s-1) = ')';
	  *s = '\0';
	}
      return;
    case product:
      gawtos(w->body.prod.c,s,true);
      while (*s++);
      *(s-1) = '*';
      gawtos(w->body.prod.d,s,true);
      return;
    }
}


void printtable(f,pref)			/* Print the entire coset table */
     FILE *f;
     char *pref;
{
  basiselt b,b1;
  gpgen g;
  vector v;
  for (b=0; b < nextbe; b++)
    {
      fprintf(f,"%s%4d: ",pref,(int)b);
      if (rowdel(b))
	{
	  fprintf(f,"Deleted - replaced by ");
	  vprint(f,replace(b));
	}
      else
      {
	for (g=0; g < nextgen; g++)
	  {
	    if (g) fprintf(f,"%s      ",pref);
	    fputs(gennames[g],f);
	    fputc(' ',f);
	    if ((v=tabent(b,g)))
	      vprint(f,v);
	    else
	      fprintf(f,"Undefined\n");
	  }
	if (TrackPreIm)
	{  
	  fprintf(f,"%sPre-image %u",pref,table[b].fmgen);
	  for (b1=0; b1 < table[b].wlen; b1++)
	    fputs(gennames[table[b].word[b1]],f);
	  fputc('\n',f);
	}
      }
    }
#ifdef INTEGRAL
  PrintLatt(f,pref);
#endif
}

void printrel(f,r)
	FILE *f;
	rell r;
{
   char s[1000];
   fprintf(f,
	   "Weight %2d %s type \t",
	   r->wt,
	   (r->type == algebra)?"algebra":"group");
   if (r->type == algebra)
    {
       printgaw(f,r->entry.a.l);
       if (r->entry.a.r)
	 {
	    fputc('=',f);
	    printgaw(f,r->entry.a.r);
	 }
    }
   else
    {
        gpwdtos(s,r->entry.g);
	fputs(s,f);
     }
}


void printrl(f,pref,rl)
     FILE *f;
     char *pref;
     rell rl;
{
  rell r;
  for (r = rl; r; r=r->next)
    {
      fprintf(f,"%s",pref);
      printrel(f,r);
      fputc('\n',f);
    }
}

#endif

#if (defined(LOGGING) || defined(SCRUT))

void vprint(f,v)			/* print a vector */
     FILE *f;
     vector v;
{
  PS(v,
     basiselt i;
     fprintf(f,"Packed vector, length %u:",(unsigned int) Len(v));
     for (i=0; i< Len(v); i++)
     {
       putc(' ',f);
       outfelt(f,BodyP(v)[i]);
     },
     sent s; sent send;
     send = BodyEndS(v);
     fprintf(f,"Sparse vector, weight %u:",(unsigned int)Wt(v));
     for (s=BodyS(v); s < send; s++)
     {
       fprintf(f," e%u.",s->loc);
       outfelt(f,s->fac);
     });
  fputc('\n',f);

}
#endif

#ifdef LOGGING
void printstacks(f)			/* print stacks */
     FILE *f;
{
  coinl c;
  deductl d;
  fprintf(f,"Coincidence Stack: \n");
  for (c=coincidences; c; c= c->next)
    {
      fprintf(f,"  0 = ");
      vprint(f,c->body.v);
    }
  fprintf(f,"Deduction Stack: \n");
  for (d=deductions; d; d= d->next)
    {
      fprintf(f,"  %d.%s = ",d->body.b,gennames[d->body.g]);
      vprint(f,d->body.v);
    }
}

void printwordl(f,wl)
     FILE *f;
     wordl wl;
{
  basiselt b;
  swl s;
  bool started = false;
  if (wl->type == wl_packed)
    {
      fputc('(',f);
	for (b = 0; b < startvecs; b++)
	  {
	    if (b) fputc(',',f);
	    printgaw(f,wl->body.p.wds[b]);
	  }
      fputc(')',f);
    }
  else
    {
      fputc('[',f);
      for (s = wl->body.s; s; s = s->next)
	{
	  if (!started)
	    started = true;
	  else
	    fputc(',',f);
	  fprintf(f,"%u,", s->loc+1);
	  printgaw(f,s->w);
	}
      fputc(']',f);
    }
}

void printsmg(f,pref)
     FILE *f;
     char *pref;
{
  smgl s;

  if (subgens->type == rels)
    printrl(f,pref,subgens->list.r);
  else
    for (s=subgens->list.s; s; s=s->next)
      {
	fprintf(f,"%s",pref);
	if (s->type == normal)
	  {
	    printwordl(f,s->body.norm.lhs);
	    if (s->body.norm.rhs)
	      {
		fputc('=',f);
		printwordl(f,s->body.norm.rhs);
	      }    
	  }
	else
	  {
	    fprintf(f,"[%u-%u,",s->body.univ.start+1,s->body.univ.end+1);
	    printgaw(f,s->body.univ.w);
	    fputc(']',f);
	  }
	    
	fputc('\n',f);
      }
  return;
} 
#endif










