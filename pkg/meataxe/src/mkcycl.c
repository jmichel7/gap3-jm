/* ========================== C MeatAxe =============================
   mkcycl.c - This program calculates a representative for each
   cyclic submodule of the kondensed modules.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: mkcycl.c,v 1.2 1997/09/11 15:43:03 gap Exp $
 *
 * $Log: mkcycl.c,v $
 * Revision 1.2  1997/09/11 15:43:03  gap
 * New version 2.2.3. AH
 *
 * Revision 2.11  1995/02/08  10:01:14  mringe
 * PL_ entfernt.
 *
 * Revision 2.10  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.9  1994/07/07  06:57:14  mringe
 * Include files.
 *
 * Revision 2.8  1994/05/18  10:13:51  mringe
 * matspin() und ycompare() durch neue Lib.funktionen ersetzt.
 *
 * Revision 2.7  1994/05/11  08:43:53  mringe
 * spinup() in Spinup() umbenannt.
 *
 * Revision 2.6  1994/03/26  06:34:04  mringe
 * basename umbenannt wg. Namenskonflikt.
 *
 * Revision 2.5  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.4  1994/01/28  08:17:50  mringe
 * Bug behoben, Aufruf von matspin() jetzt richtig.
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
 * Revision 1.13  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.12  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.11  1993/10/02  16:23:02  mringe
 * matread() und matwrite() in matload() bzw. matsave() umbenannt.
 *
 * Revision 1.10  1993/09/20  20:36:01  mringe
 * Neu: -TGQV, help verbessert.
 *
 * Revision 1.9  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.8  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.7  1993/07/05  16:51:57  mringe
 * N Erzeuger, help().
 *
 * Revision 1.6  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.5  1993/02/15  13:49:24  mringe
 * Funktionen aus cyclic.c -> yyy-Lib verschoben.
 *
 * Revision 1.4  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.3  1992/10/01  13:50:10  mringe
 * Header eingef"ugt.
 *
 * Revision 1.2  1992/07/22  07:10:30  mringe
 * Changed 'global.h' to 'lattice.h'
 *
 * Revision 1.1  1992/05/24  08:22:56  mringe
 * Initial revision
 *
 */


#include <string.h>
#include "meataxe.h"
#include "lattice.h"
#include "files.h"



/* ------------------------------------------------------------------
   Function prototypes
   ------------------------------------------------------------------ */

static void init(int argc, char *argv[]);
void cyclic(int cf);
void Spinup(matrix_t *seed);
void writeresult(int cf);


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static char *helptext[] = {
"SYNTAX",
"    mkcycl [<Options>] <Name>",
"",
"OPTIONS",
"    -Q            Quiet, no messages.",
"    -V            Verbose, more messages.",
"    -G            GAP output (implies -Q).",
"    -T <Seconds>  Set CPU time limit",
"",
"FILES",
"    <Name><Dim><a>.<N>k    i    Generators in kondensed modules",
"    <Name><Dim><a>.np      i    Kondensed peak word",
"    <Name><Dim><a>.v       o    Cyclic submodules",
NULL,
};

static proginfo_t pinfo =
   { "mkcycl", "Find Cyclic Submodules",
     "$Revision: 1.2 $", helptext };


int opt_G = 0;

matrix_t *gen[MAXGEN+1];	/* Generators and peak word*/
long kdim;			/* Dimension of kondensed module */
int ncycl;			/* Number of cyclic submodules found */
long count;			/* Number of vectors tried */
matrix_t *cmod[MAXCYCL];	/* Cyclic submodules */


static void init(argc, argv)
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
    if (opt_ind != argc-1) errexit(ERR_BADUSAGE,"mkcycl");
    setbasename(argv[argc-1]);
    readcfinfo();
}




int main(argc, argv)
int argc;
char *argv[];

{
    int i;

    init(argc, argv);
    MESSAGE(0,("\n*** FIND CYCLIC SUBMODULES ***\n\n"));
    for (i = 0; i < ncf; ++i)
	cyclic(i);
    if (msg_level >= 0)
    {
        printf("\n");
        prtimes();
    }
    return 0;
}



void cyclic(cf)
int cf;

{
    matrix_t *seed;
    char fn[200], fn1[200];
    long k, maxk;
    FEL f;
    int i;

    MESSAGE(0,("%s%s: ",cfbasename,cfname(cf)));
    fflush(stdout);

    /* Read generators
       --------------- */
    for (i = 0; i < ngen; ++i)
    {
	sprintf(fn,"%s%s.%dk",cfbasename,cfname(cf),i+1);
	gen[i] = matload(fn);
	if (gen[i]->nor != gen[i]->noc)
	    errexit(ERR_NOTSQUARE,fn);
	if (i == 0)
	{
	    kdim = gen[0]->nor;
	    strcpy(fn1,fn);
	}
	else
	{
	    char buf[200];
	    if (gen[i]->nor != kdim)
	    {   sprintf("%s and %s",fn,fn1);
		errexit(ERR_INCOMPAT,buf);
	    }
	}
    }

    /* Read the kondensed peak word
       ---------------------------- */
    sprintf(fn,"%s%s.np",cfbasename,cfname(cf));
    gen[ngen] = matload(fn);

    /* Initialize the seed vector generator
       ------------------------------------ */
    seed = matalloc(zfl,(long)1,kdim);
    zmulrow(seed->d,(FEL)0);
    ncycl = 0;
    count = 0;
    maxk = 1;

    /* Spin up each seed vector
       ------------------------ */
    while (1)
    {
	for (k=1; 1; ++k)
	{
	    /* This won't work with bigzzz! */

	    f = zextract(seed->d,k) + 1;
	    if (k == maxk && f == 2)
	    {
		if (k == kdim)
		{	writeresult(cf);
			matfree(seed);
			return;
		}
		else
		{	f = 0;
			++maxk;
		}
	    }
	    if (f == (FEL) zfl)
		f = F_ZERO;
	    zinsert(seed->d,k,f);
	    if (f != 0)
	    {
		Spinup(seed);
		break;
	    }
	}
    }
}



void Spinup(seed)
matrix_t *seed;

{
    matrix_t *sub;
    int i;

    ++count;
    spinup(seed,ngen+1,gen,SPL_SEED_FIRST,&sub);
    for (i = 0; i < ncycl; ++i)
    {
	if (sub->nor == cmod[i]->nor &&
	    spccomp(cmod[i],sub,(long)1) == 0)
	{
	    matfree(sub);	/* Module already in list */
	    return;
	}
    }
    if (ncycl >= MAXCYCL)
	FATAL("TOO MANY CYCLIC SUBMODULES");
    cmod[ncycl] = sub;
    ++ncycl;
}



void writeresult(cf)
int cf;

{
    int i;
    char fn[100];
    FILE *f;

    MESSAGE(0,("%d cyclic submodule%s (%ld vectors tried)\n",
	ncycl,ncycl == 1 ? " " : "s",count));
    sprintf(fn,"%s%s.v",cfbasename,cfname(cf));
    zsetlen(kdim);
    if ((f = zwritehdr(fn,zfl,(long)ncycl,kdim)) == NULL)
	errexit(-1,fn);
    for (i = 0; i < ncycl; ++i)
	zwritevec(f,cmod[i]->d,1);
    fclose(f);

    /* Clean up
       -------- */
    for (i = 0; i < ncycl; ++i)
	matfree(cmod[i]);
    for (i = 0; i <= ngen; ++i)
	matfree(gen[i]);
}


