/* ========================== C MeatAxe =============================
   mkpeak.c - Calculates peak words.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: mkpeak.c,v 1.1.1.1 1996/12/11 12:40:21 werner Exp $
 *
 * $Log: mkpeak.c,v $
 * Revision 1.1.1.1  1996/12/11 12:40:21  werner
 * Preparing 3.4.4 for release
 *
 * Revision 2.6  1994/03/26  06:34:04  mringe
 * basename umbenannt wg. Namenskonflikt.
 *
 * Revision 2.5  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.4  1994/01/24  08:39:41  mringe
 * Helptext: -e.
 *
 * Revision 2.3  1993/12/08  11:33:02  mringe
 * Neue CPU time - Funktionen.
 *
 * Revision 2.2  1993/10/22  16:08:19  mringe
 * Neues Numerierungsschema fuer irreduzible.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.18  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.17  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.16  1993/10/02  16:23:02  mringe
 * matread() und matwrite() in matload() bzw. matsave() umbenannt.
 *
 * Revision 1.15  1993/09/20  19:05:30  mringe
 * Doppelte Aufrufe entfernt: setbasename(), mtxinit()
 *
 * Revision 1.14  1993/09/20  17:40:14  mringe
 * help, Optionen -T, -G, -Q, -V.
 *
 * Revision 1.13  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.12  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.11  1993/02/16  18:32:46  mringe
 * string.h und stdio.h werden jetzt in meataxe.h included.
 *
 * Revision 1.10  1993/02/12  17:08:07  mringe
 * Woerter mit N Erzeugern.
 *
 * Revision 1.9  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.8  1992/10/01  13:50:10  mringe
 * Header eingef"ugt.
 *
 * Revision 1.7  1992/07/22  07:10:30  mringe
 * Changed 'global.h' to 'lattice.h'
 *
 * Revision 1.6  1992/07/15  09:27:50  mringe
 * Some minor changes.
 *
 * Revision 1.5  1992/07/10  16:03:13  mringe
 * Finde den Koerper beim Einlesen der Irred. heraus
 *
 * Revision 1.4  1992/07/02  15:00:13  mringe
 * Fange bei Wort 1 (statt 0) an.
 *
 * Revision 1.3  1992/07/01  14:37:55  mringe
 * Argumente aus nicht-ANSI Prototypes entfernt.
 *
 * Revision 1.2  1992/05/25  17:19:14  mringe
 * *** empty log message ***
 *
 */




#include "meataxe.h"
#include "lattice.h"
#include "words.h"
#include "files.h"


#define MAXLOCK 100


/* ------------------------------------------------------------------
   Function prototypes
   ------------------------------------------------------------------ */

static int isexcluded _PL((long w));
static void list _PL((char *c));
static void init _PL((int argc, char *argv[]));
int finished _PL((void));
int try _PL((long w));


/* ------------------------------------------------------------------
   Global variables
   ------------------------------------------------------------------ */

static char *helptext[] = {
"SYNTAX",
"    mkpeak [<Options>] <Name>",
"",
"OPTIONS",
"    -Q            Quiet, no messages.",
"    -V            Verbose, more messages.",
"    -G            GAP output (implies -Q).",
"    -T <Seconds>  Set CPU time limit",
"    -e <List>     Exclude words from search. <List> is a",
"                  comma-separated list of numbers or ranges (X-Y).",
"",
"FILES",
"    Input files are produced by CHOP (see `chop -help').",
"    Output is written to <Name>.cfinfo.",
NULL};

static proginfo_t pinfo =
   { "mkpeak", "Find Peak Words", "$Revision: 1.1.1.1 $", helptext };


int opt_G = 0;
basis_t basis[MAXCF];
matrix_t *tmp[MAXCF];
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

static void list(c)
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
    matrix_t *gen[MAXGEN];

    MESSAGE(0,("\n*** FIND PEAK WORDS ***\n\n"));
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
	    gen[k] = matload(fn);
	}
	initbasis(gen,ngen,basis+i);
	tmp[i] = matalloc(gen[0]->fl,cfinfo[i].dim,cfinfo[i].dim);
	for (k = 0; k < ngen; ++k)
	    matfree(gen[k]);
	cfinfo[i].peakword = -1;
    }
}


/* ------------------------------------------------------------------
   finished() - Find out if a peak word has been found for each
	module.
   ------------------------------------------------------------------ */

int finished()

{
    int i;

    for (i = 0; i < ncf; ++i)
	if (cfinfo[i].peakword == -1) return 0;
    return 1;
}


/* ------------------------------------------------------------------
   try() - Try another word
   ------------------------------------------------------------------ */

int try(w)
long w;

{
    int i, ppos = -1;
    long nul;

    if (isexcluded(w)) return -1;
    MESSAGE(1,("Word %ld:",w));
    for (i = 0; i < ncf; ++i)  /* For each composition factor... */
    {
	mkword(basis+i,w);
	matmove(tmp[i],basis[i].w);
	nul = nullity(tmp[i]);
        MESSAGE(1,(" %ld",nul));
	if (nul != 0 && nul != cfinfo[i].spl)
	    return -1;
	if (nul == cfinfo[i].spl)
	{
	    if (ppos >= 0 || cfinfo[i].peakword > 0)
		return -1;  /* Nullity should be 0 */
	    matmove(tmp[i],basis[i].w);
	    matmul(tmp[i],basis[i].w);
	    nul = nullity(tmp[i]);
	    if (nul != cfinfo[i].spl)
		return -1;  /* Nullity is not stable */
	    ppos = i;
	}
    }
    if (ppos > -1) /* we have found a new peak word */
	cfinfo[ppos].peakword = w;
    return ppos;
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
    while ((i = zgetopt("Ge:")) != OPT_END)
    {
	switch (i)
	{
	    case 'G': opt_G = 1; msg_level = -100; break;
	    case 'e':
		list(opt_text);
		break;
	}
    }

    if (opt_ind != argc-1) errexit(ERR_BADUSAGE,"mkpeak");
    mtxinit();
    setbasename(argv[argc-1]);

    init(argc,argv);
    for (w = 1; !finished(); w = nextword(w))
    {
	i = try(w);
	MESSAGE(1,("\n"));
	if (i >= 0)
	{
	    MESSAGE(0,("Peak word for %s%s is %ld, nullity = %ld\n",
		cfbasename,cfname(i),cfinfo[i].peakword,
		cfinfo[i].spl));
	    fflush(stdout);
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


