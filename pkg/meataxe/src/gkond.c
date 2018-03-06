/* ========================== C MeatAxe =============================
   This program performs a generalized kondensation for each
   composition factor using the peak words found by mkpeak.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: gkond.c,v 1.1.1.1 1996/12/11 12:40:22 werner Exp $
 *
 * $Log: gkond.c,v $
 * Revision 1.1.1.1  1996/12/11 12:40:22  werner
 * Preparing 3.4.4 for release
 *
 * Revision 2.8  1994/03/26  06:34:04  mringe
 * basename umbenannt wg. Namenskonflikt.
 *
 * Revision 2.7  1994/02/19  15:18:49  mringe
 * unistd.h raus.
 *
 * Revision 2.6  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.5  1993/12/08  11:33:02  mringe
 * Neue CPU time - Funktionen.
 *
 * Revision 2.4  1993/12/02  17:56:26  mringe
 * Include stdlib.h fuer free().
 *
 * Revision 2.3  1993/10/22  16:08:19  mringe
 * Neues Numerierungsschema fuer irreduzible.
 *
 * Revision 2.2  1993/10/21  08:37:53  mringe
 * Compiler warnings.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.16  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.15  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.14  1993/10/02  16:23:02  mringe
 * matread() und matwrite() in matload() bzw. matsave() umbenannt.
 *
 * Revision 1.13  1993/09/20  20:27:36  mringe
 * *** empty log message ***
 *
 * Revision 1.12  1993/09/20  19:47:59  mringe
 * *** empty log message ***
 *
 * Revision 1.11  1993/09/20  19:30:22  mringe
 * Optionen -GQV, help.
 *
 * Revision 1.10  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.9  1993/02/16  18:32:46  mringe
 * string.h und stdio.h werden jetzt in meataxe.h included.
 *
 * Revision 1.8  1993/02/12  17:08:07  mringe
 * Woerter mit N Erzeugern.
 *
 * Revision 1.7  1993/02/10  20:46:56  mringe
 * Bug in gkond() behoben.
 *
 * Revision 1.6  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.5  1993/01/21  13:29:45  mringe
 * Filenamen der Erzeuger korrigiert
 *
 * Revision 1.4  1993/01/15  15:12:36  mringe
 * Auf n Erzeuger erweitert.
 *
 * Revision 1.3  1992/07/22  07:10:30  mringe
 * Changed 'global.h' to 'lattice.h'
 *
 * Revision 1.2  1992/07/15  09:25:55  mringe
 * Some minor changes.
 *
 * Revision 1.1  1992/05/26  07:29:08  mringe
 * Initial revision
 *
 */


#include <string.h>
#include <stdlib.h>

#include "meataxe.h"
#include "lattice.h"
#include "files.h"
#include "words.h"


/* ------------------------------------------------------------------
   Function prototypes
   ------------------------------------------------------------------ */

static void init _PL((void));
static void power _PL((int i));
static void gkond _PL((int i, matrix_t *k, matrix_t *w, char *name));
static void kond _PL((int cf));


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

long dim;
long quotdim;
static matrix_t *gen[MAXGEN];		/* Generators */
static matrix_t *tmp1, *tmp2;
static basis_t basis;
int opt_G = 0;

static char *helptext[] = {
"SYNTAX",
"    gkond [<Options>] <Name>",
"",
"OPTIONS",
"    -Q            Quiet, no messages.",
"    -V            Verbose, more messages.",
"    -G            GAP output (implies -Q).",
"    -T <Seconds>  Set CPU time limit",
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
   { "gkond", "Generalized Kondensation",
     "$Revision: 1.1.1.1 $", helptext };



/* ------------------------------------------------------------------
   init()
   ------------------------------------------------------------------ */

static void init()

{
    char fn[MAXBASENAME+10];
    int i;
    long fl;

    mtxinit();
    readcfinfo();

    for (i = 0; i < ngen; ++i)
    {
	sprintf(fn,"%s.%d",cfbasename,i+1);
	gen[i] = matload(fn);
    }
    fl = gen[0]->fl;
    initbasis(gen,ngen,&basis);
    dim = gen[0]->nor;
    tmp1 = matalloc(fl,dim,dim);
    tmp2 = matalloc(fl,dim,dim);
}


/* ------------------------------------------------------------------
   power()
   ------------------------------------------------------------------ */

static void power(i)
int i;

{
    long nul1, nul2, k0;

    mkword(&basis,cfinfo[i].peakword);
    matmove(tmp1,basis.w);
    matmove(tmp2,tmp1);
    k0 = 1;
    nul2 = nullity(tmp2);
    nul1 = -1;
    while (nul2 != nul1)	/* while nullity not stable */
    {
	k0 <<= 1;
	nul1 = nul2;
	matmove(tmp2,tmp1);
	matmul(tmp2,tmp1);
	matmove(tmp1,tmp2);
	matmul(tmp1,basis.w);
	matmove(tmp2,tmp1);
	nul2 = nullity(tmp2);
    }
    MESSAGE(0,("word=%ld, pwr=%ld, nul=%ld, mult=%ld\n",
	  cfinfo[i].peakword,k0,nul2,cfinfo[i].mult));
    fflush(stdout);

    /* Consistency check
       ----------------- */
    if (nul2 != cfinfo[i].mult * cfinfo[i].spl)
	FATAL("Something is wrong...");
    quotdim = nul2;
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
    free(x2->d);
    x2->d = zquot(x1->d,x1->nor);

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

    MESSAGE(0,("%s%s: ",cfbasename,cfname(cf)));
    fflush(stdout);
		
    /* Make the peak word, find its stable power,
       calculate both kernel and image.
       ------------------------------------------- */
    power(cf);
    kern = nullsp(tmp1);
    bild = echelon(tmp1);
    zquotinit(bild->d,bild->nor,NULL);

    /* Write out the image
       ------------------- */
    sprintf(fn,"%s%s.im",cfbasename,cfname(cf));
    matsave(bild,fn);

    /* Write out the `unkondense matrix'
       --------------------------------- */
    m = matalloc(kern->fl,kern->nor,quotdim);
    free(m->d);
    m->d = zquot(kern->d,kern->nor);
    k = matinv(m);
    matmul(k,kern);
    sprintf(fn,"%s%s.k",cfbasename,cfname(cf));
    matsave(k,fn);

    /* Kondense all generators
       ----------------------- */
    for (j = 0; j < ngen; ++j)
    {
	sprintf(fn,"%dk",j+1);
	gkond(cf,k,gen[j],fn);
    }

    /* Kondense the peak word
       ---------------------- */
    gkond(cf,k,basis.w,"np");

    matfree(m); matfree(k);
    matfree(kern); matfree(bild);
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{
    int i;


    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("G")) != OPT_END)
    {
	switch (i)
	{
	    case 'G': opt_G = 1; msg_level = -100; break;
	}
    }
    if (opt_ind != argc-1) errexit(ERR_BADUSAGE,"gkond");

    MESSAGE(0,("\n*** GENERALIZED KONDENSATION ***\n\n"));
    setbasename(argv[argc-1]);
    init();
    for (i = 0; i < ncf; ++i)
	kond(i);
    if (msg_level >= 0)
    {
        printf("\n");
        prtimes();
    }
    return 0;
}





