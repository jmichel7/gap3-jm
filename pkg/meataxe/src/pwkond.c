/* ========================== C MeatAxe =============================
   This program finds the peak words and does the generalized
   kondensation for each composition factor.

   (C) Copyright 1994 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: pwkond.c,v 1.2 1997/09/11 15:43:18 gap Exp $
 *
 * $Log: pwkond.c,v $
 * Revision 1.2  1997/09/11 15:43:18  gap
 * New version 2.2.3. AH
 *
 * Revision 1.9  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 1.8  1994/09/10  10:36:49  mringe
 * Option -p: Benutze Polynome, nicht nur Linearfaktoren.
 *
 * Revision 1.7  1994/09/08  07:52:57  mringe
 * Benutzt matinsert_() statt matinsert()
 *
 * Revision 1.6  1994/08/29  17:03:00  mringe
 * Peakwoerter: Benutze Paar aus WortNr,Polynom.
 *
 * Revision 1.5  1994/08/25  12:53:58  mringe
 * Suche nach beliebigen Eigenwerten.
 *
 * Revision 1.4  1994/07/29  07:43:08  mringe
 * Neuer Wortgenerator.
 *
 * Revision 1.3  1994/07/22  18:54:57  mringe
 * Neue Funktionen nullity(), echelon(), nullspace()
 *
 * Revision 1.2  1994/07/07  09:40:33  mringe
 * *** empty log message ***
 *
 * Revision 1.1  1994/07/07  09:31:03  mringe
 * Initial revision, aus mkpeak.c und gkond.c zusammengesetzt.
 *
 */


#include <string.h>
#include <stdlib.h>

#include "meataxe.h"
#include "lattice.h"
#include "files.h"



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

long dim;
long quotdim;
static matrix_t *gen[MAXGEN];	/* Generators */
static wgdata_t *basis;
int opt_G = 0;
int opt_n = 0;			/* No kondensation, PW only */
int opt_p = 0;			/* Use full polynomials in PW search */
matrix_t *peakword;

static char *helptext[] = {
"SYNTAX",
"    pwkond [<Options>] <Name>",
"",
"OPTIONS",
"    -Q            Quiet, no messages.",
"    -V            Verbose, more messages.",
"    -G            GAP output (implies -Q).",
"    -T <Seconds>  Set CPU time limit",
"    -n            Find peak words only, do not kondense",
"    -p            Use full polynomials in peak word search",
"    -e <List>     Exclude words from search. <List> is a",
"                  comma-separated list of numbers or ranges (X-Y).",
"",
"FILES",
"    Input files are produced by CHOP (see `chop -help').",
"    Output is written to the following files:",
"       <Name><Dim><a>.<N>k        Kondensed generators",
"       <Name><Dim><a>.np          Kondensed peak word",
"       <Name><Dim><a>.im          Image used for kondenstion",
"       <Name><Dim><a>.k           Unkondense matrix (used by mkinc)",
NULL};

static proginfo_t pinfo =
   { "pwkond", "Peakword Kondensation",
     "$Revision: 1.2 $", helptext };



/* ------------------------------------------------------------------
   init()
   ------------------------------------------------------------------ */

static void init2()

{
    char fn[MAXBASENAME+10];
    int i;
    long fl;

    for (i = 0; i < ngen; ++i)
    {
	sprintf(fn,"%s.%d",cfbasename,i+1);
	gen[i] = matload(fn);
    }
    fl = gen[0]->fl;
    basis = WGInit(ngen,gen);
    dim = gen[0]->nor;
}


/* ------------------------------------------------------------------
   power() - Find the stable power of the peak word
   ------------------------------------------------------------------ */

static matrix_t *power(i)
int i;

{
    long nul1, nul2, k0;
    matrix_t *tmp;

    k0 = 1;
    tmp = matdup(peakword);
    do
    {   matrix_t *tmp2 = matmul(matdup(tmp),tmp);
	nul1 = nullity__(tmp);
	tmp = matdup(tmp2);
	nul2 = nullity__(tmp2);
	k0 *= 2;
    }
    while (nul2 != nul1);

    MESSAGE(0,("pwr=%ld, nul=%ld, mult=%ld, ",k0/1,nul2,
	cfinfo[i].mult));
    if (nul2 != cfinfo[i].mult * cfinfo[i].spl)
	FATAL("Something is wrong...");
    quotdim = nul2;
    return tmp;
}



/* ------------------------------------------------------------------
   gkond()
   ------------------------------------------------------------------ */

static void gkond(i,k,w,name)
int i;
matrix_t *k, *w;
char *name;

{
    char fn[MAXBASENAME+10];
    matrix_t *x1, *x2;

    x1 = matdup(k);
    matmul(x1,w);

    x2 = matalloc(x1->fl,x1->nor,quotdim);
    zquot(x1->d,x1->nor,x2->d);

    sprintf(fn,"%s%s.%s",cfbasename,cfname(i),name);
    matsave(x2,fn);
    matfree(x1);
    matfree(x2);
}


/* ------------------------------------------------------------------
   kond() - Generalized kondensation for one irreducible
   ------------------------------------------------------------------ */

static void kond(cf)
int cf;

{
    char fn[MAXBASENAME+10];
    matrix_t *kern, *bild, *m, *k;
    int j;
		
    /* Make the peak word, find its stable power,
       calculate both kernel and image.
       ------------------------------------------- */
    MESSAGE(0,("Kondensing %s%s: word=%ld, ",cfbasename,cfname(cf),
    	cfinfo[cf].peakword));
    peakword = MakeWord(basis,cfinfo[cf].peakword);
    matinsert_(peakword,cfinfo[cf].peakpol);
    bild = power(cf);
    kern = nullspace_(bild);
    echelon_(bild);
    zquotinit(bild->d,bild->nor,NULL);

    /* Write out the image
       ------------------- */
    sprintf(fn,"%s%s.im",cfbasename,cfname(cf));
    matsave(bild,fn);

    /* Write out the `unkondense matrix'
       --------------------------------- */
    m = matalloc(kern->fl,kern->nor,quotdim);
    zquot(kern->d,kern->nor,m->d);
    k = matinv(m);
    matmul(k,kern);
    sprintf(fn,"%s%s.k",cfbasename,cfname(cf));
    matsave(k,fn);

    /* Kondense all generators
       ----------------------- */
    MESSAGE(0,("("));
    for (j = 0; j < ngen; ++j)
    {
	sprintf(fn,"%dk",j+1);
	gkond(cf,k,gen[j],fn);
	MESSAGE(0,(" %d",j+1));
    }
    MESSAGE(0,(" )\n"));

    /* Kondense the peak word
       ---------------------- */
    gkond(cf,k,peakword,"np");

    matfree(m); matfree(k);
    matfree(kern);
    matfree(bild);
    matfree(peakword);
}




#define MAXLOCK 100
wgdata_t *cfbasis[MAXCF];
long exclude[MAXLOCK][2];
int nexclude=0;



static int isexcluded(w)
long w;

{	int i;

	for (i = 0; i < nexclude; ++i)
		if (w >= exclude[i][0] && w <= exclude[i][1])
			return 1;
	return 0;
}

static void parselist(c)
char *c;

{	long a, b;

	while (*c != 0)
	{	a = b = 0;
		while (*c >= '0' && *c <= '9')
			a = a * 10 + (*c++ - '0');
		if (*c == '-')
		{	++c;
			while (*c >= '0' && *c <= '9')
				b = b * 10 + (*c++ - '0');
		}
		else
			b = a;
		if (a == 0 || b == 0 || a > b) 
			FATAL("BAD ARGUMENTS");
		exclude[nexclude][0] = a;
		exclude[nexclude][1] = b;
		++nexclude;
		if (*c == ',') ++c;
	}
}


static void init(argc, argv)
int argc;
char *argv[];

{
    char fn[50];
    int i;
    matrix_t *cfgen[MAXGEN];

    readcfinfo();

    /* Read the generators for each composition factor
       ----------------------------------------------- */
    for (i = 0; i < ncf; ++i)
    {
	int k;
	for (k = 0; k < ngen; ++k)
	{
	    sprintf(fn,"%s%s.%d",cfbasename,cfname(i),k+1);
	    MESSAGE(1,("Reading %s\n",fn));
	    cfgen[k] = matload(fn);
	}
	cfbasis[i] = WGInit(ngen,cfgen);
	cfinfo[i].peakword = -1;
    }
}


/* ------------------------------------------------------------------
   finished() - Find out if a peak word has been found for each
	module.
   ------------------------------------------------------------------ */

static int finished()

{
    int i;

    for (i = 0; i < ncf; ++i)
	if (cfinfo[i].peakword == -1) return 0;
    return 1;
}


/* ------------------------------------------------------------------
   try() - Try another word
   ------------------------------------------------------------------ */

static void addid(m,f)
matrix_t *m;
FEL f;
{
    PTR x;
    long i;
    if (f == F_ZERO) return;
    zsetlen(m->noc);
    for (i = 1, x = m->d; i <= m->nor; ++i, zadvance(&x,1))
	zinsert(x,i,zadd(zextract(x,i),f));
}


static int try2(w,f)
long w;
FEL f;

{
    int i, ppos = -1;
    long nul;

    if (isexcluded(w)) return -1;
    MESSAGE(3,("Word %ld+%ldI:",w,zftoi(f)));
    for (i = 0; i < ncf; ++i)  /* For each composition factor... */
    {
	matrix_t *word;
	
	word = MakeWord(cfbasis[i],w);
	addid(word,f);
	nul = nullity__(matdup(word));
        MESSAGE(3,(" %ld",nul));
	if (nul != 0 && nul != cfinfo[i].spl)
	{
	    matfree(word);
	    MESSAGE(3,("failed\n"));
	    return -1;
	}
	if (nul == cfinfo[i].spl)
	{
	    if (ppos >= 0 || cfinfo[i].peakword > 0)
	    {
	    	MESSAGE(3,("failed\n"));
		return -1;  /* Nullity should be 0 */
	    }
	    nul = nullity__(matmul(matdup(word),word));
	    if (nul != cfinfo[i].spl)
	    {
		matfree(word);
	    	MESSAGE(3,("failed (nullity not stable)\n"));
		return -1;  /* Nullity is not stable */
	    }
	    ppos = i;
	}
        matfree(word);
    }
    if (ppos > -1) /* we have found a new peak word */
    {
	poly_t *pp;
	cfinfo[ppos].peakword = w;
	pp = polalloc(zfl,1);
	pp->buf[0] = f;
	cfinfo[ppos].peakpol = pp;
    }
    MESSAGE(3,("\n"));
    return ppos;
}

static int try(w)
long w;

{
    long f;
    for (f = 0; f < zfl; ++f)
    {
        int result = try2(w,zitof(f));
        if (result >= 0) return result;
    }
    return -1;
}


/* ------------------------------------------------------------------
   MinPol() - Returns the minimal polynomial in factorized form.
   ------------------------------------------------------------------ */

static fpoly_t *MinPol(matrix_t *m)

{
    poly_t *f = minpolfactor(m);
    fpoly_t *mp = fpolalloc();
    while (f != NULL)
    {
    	fpoly_t *ff = factorization(f);
    	fpolmul(mp,ff);
    	polfree(f);
    	fpolfree(ff);
	f = minpolfactor(NULL);
    }
    return mp;
}

/* ------------------------------------------------------------------
   tryp2()
   ------------------------------------------------------------------ */

static int tryp2(long w, int cf, poly_t *pol)

{
    int i;

for (i = 0; i < ncf; ++i)
{
    matrix_t *word, *wordp;
    long nul;

    if (i == cf)
    {
    	continue;
    }
    word = MakeWord(cfbasis[i],w);
    wordp = matinsert(word,pol);
    matfree(word);
    nul = nullity__(wordp);
    if (nul != 0) return -1;
}
return 0;
}

/* ------------------------------------------------------------------
   try_p() - Try another word using polynomials
   ------------------------------------------------------------------ */

static int try_p(w)
long w;

{
    int i;

    if (isexcluded(w)) return -1;

for (i = 0; i < ncf; ++i)
{
    matrix_t *word;
    fpoly_t *mp;
    int k;

    if (cfinfo[i].peakword > 0) continue;
    word = MakeWord(cfbasis[i],w);
    mp = MinPol(word);
    if (MSG3)
    {
	int j;
	printf("%s, minpol =\n",cfname(i));
	for (j = 0; j < mp->len; ++j)
	{
	    polprint(NULL,mp->p[j]);
	    printf(" ^ %ld\n",mp->e[j]);
	}
    }
    for (k = 0; k < mp->len; ++k)
    {
	if (mp->p[k]->deg * mp->e[k] == cfinfo[i].spl)
	{
	   matrix_t *wp, *wp2;
	   long nul;

	   if (MSG3)
	   {
	       printf("%s, ",cfname(i));
	       polprint("factor",mp->p[k]);
	   }
	   if (tryp2(w,i,mp->p[k]) == -1) continue;
	   wp = matinsert(word,mp->p[k]);
	   wp2 = matmul(matdup(wp),wp);
	   matfree(wp);
	   nul = nullity__(wp2);
	   if (nul != cfinfo[i].spl) continue;
	   break;
	}
    }
    matfree(word);
    if (k < mp->len)
    {
	cfinfo[i].peakword = w;
	cfinfo[i].peakpol = mp->p[k];
	return i;
    }
    else
    	fpolfree(mp);
}

return -1;

}

/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */


int main(argc, argv)
int argc;
char *argv[];

{
    int i;
    long w;

    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("Gpne:")) != OPT_END)
    {
	switch (i)
	{
	    case 'G': opt_G = 1; msg_level = -100; break;
	    case 'n': opt_n = 1; break;
	    case 'p': opt_p = 1; break;
	    case 'e':
		parselist(opt_text);
		break;
	}
    }
    if (opt_ind != argc-1) errexit(ERR_BADUSAGE,pinfo.name);
    setbasename(argv[argc-1]);
    MESSAGE(0,("\n*** PEAK WORD KONDENSATION ***\n\n"));
    init(argc,argv);
    init2();

    for (w = 1; !finished(); ++w)
    {
	if ((MSG1 && (w % 50 == 0)) || MSG2)
	{
	    printf("Word %ld\n",w);
	    fflush(stdout);
	}
	i = opt_p ? try_p(w) : try(w);
	if (i >= 0)
	{
	    if (MSG0)
	    {
	    	printf("Peak word for %s%s is %ld (",cfbasename,
		    cfname(i), cfinfo[i].peakword);
	    	polprint(NULL,cfinfo[i].peakpol);
	    	printf("), nullity=%ld\n",cfinfo[i].spl);
	    	fflush(stdout);
	    }
	    if (!opt_n) kond(i);
	}
    }
    writecfinfo();
    if (opt_G)
    {
	printf("MeatAxe.PeakWords := [");
	for (i = 0; i < ncf; ++i)
	{
	    if (i > 0) printf(",");
	    printf("%ld",cfinfo[i].peakword);
	}
	printf("];\n");
    }
    if (msg_level >= 0)
    {
	printf("\n");
	prtimes();
    }
    return 0;
}
