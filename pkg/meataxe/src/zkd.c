/* ========================== C MeatAxe =============================
   zkd.c - Kondense a permutation

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zkd.c,v 1.2 1997/09/11 15:43:59 gap Exp $
 *
 * $Log: zkd.c,v $
 * Revision 1.2  1997/09/11 15:43:59  gap
 * New version 2.2.3. AH
 *
 * Revision 2.7  1995/02/09  14:03:44  mringe
 * Neues Format fuer Bahnlaengendatei
 *
 * Revision 2.6  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.5  1994/03/12  13:03:53  mringe
 * Verschiedene Bugs beseitigt.
 *
 * Revision 2.4  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.3  1993/10/26  10:47:35  mringe
 * Compiler Warnings.
 *
 * Revision 2.2  1993/10/21  21:57:35  mringe
 * Permutationen.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.6  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.5  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.4  1993/08/09  07:15:48  mringe
 * help() eingebaut.
 *
 * Revision 1.3  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.2  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.1  1992/05/26  17:39:21  mringe
 * Initial revision
 *
 */


#include <stdlib.h>
#include "meataxe.h"




/* ------------------------------------------------------------------
   Global variables
   ------------------------------------------------------------------ */

static char *helptext[] = {
"SYNTAX",
"    zkd <Field> <OrbSz> <Perm> <Kond>",
"",
"ARGUMENTS",
"    <Field>   The field to use for kondensation.",
"",
"FILES",
"    <OrbSz>   i  Orbit sizes (produced by ZMO)",
"    <Perm>    i  Permutation to be kondensed.",
"    <Kond>    o  The result, a matrix",
NULL};
static proginfo_t pinfo =
   { "zkd", "Kondense A Permutation", "$Revision: 1.2 $", helptext };

char *orbname, *permname, *kondname;
FILE *orbfile, *permfile, *kondfile;
long fl;		/* Field order */
long prime;		/* Characteristic */
long ppow;		/* l.c.m. of the p-parts of orbit sizes */
long degree;		/* Degree of the permutation */
long norb;		/* Number of orbits */

long *orbsize;		/* Array of orbit sizes */
PTR hsz;
PTR perm;		/* The permutation to be kondensed */
PTR m1;			/* One row of the output matrix */



/* ------------------------------------------------------------------
   init() - Initialize everything
   ------------------------------------------------------------------ */

static void init(int argc, char *argv[])

{
    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while (zgetopt("") != OPT_END);
    if (opt_ind != argc-4) errexit(ERR_NARGS,"zkd");
    opt_text_ptr = argv[opt_ind];
    fl = getint();
    if (fl == GETINT_ERR || *opt_text_ptr != 0)
	errexit(ERR_USAGE,"<Field>");
    orbname=argv[opt_ind+1];
    permname=argv[opt_ind+2];
    kondname=argv[opt_ind+3];

}


/* ------------------------------------------------------------------
   readdata() - Open and read input files
   ------------------------------------------------------------------ */

static void readdata()

{
    long fl1, n, nperm;
    int i;

    /* Open input files and check
       -------------------------- */
    if ((orbfile = zreadhdr(orbname,&fl1,&n,&norb)) == NULL)
	errexit(-1,orbname);
    if (n != 1 || fl1 != -T_IMAT)
	errexit(ERR_FILEFMT,orbname);
    if ((permfile = zreadhdr(permname,&fl1,&degree,&nperm)) == NULL)
	errexit(-1,permname);
    if (fl1 != -1) errexit(ERR_NOTPERM,permname);

    /* Allocate memory, read orbit sizes and permutation
       ------------------------------------------------- */
    orbsize = (long *) malloc(sizeof(long) * (norb+1));
    perm = (PTR) malloc(sizeof(long) * degree);
    if (orbsize == NULL || perm == NULL) errexit(ERR_NOMEM,"zkd");
    if (zreadlong(orbfile,orbsize,norb) != norb)
	errexit(ERR_FILEREAD,orbname);
    if (zreadlong(permfile,(long *)perm,degree) != degree)
	errexit(ERR_FILEREAD,permname);
    fclose(orbfile);
    fclose(permfile);

    /* Allocate memory for output (1 row) and hsz
       ------------------------------------------ */
    zsetfield(fl);
    zsetlen(norb);
    m1 = zalloc((long)1);
    hsz = zalloc((long)1);

    /* Find the largest power of the characteristic
       which divides any of the orbit sizes
       -------------------------------------------- */
    prime = zchar;
    if (MSG1)
	printf("Kondensation over GF(%ld), Characteristic is %ld\n",
		fl,zchar);
    ppow = prime;
    for (i = 0; i < (int) norb; ++i)
    {
	n = orbsize[i];
	while (n % ppow == 0)
	    ppow *= prime;
    }
    ppow /= prime;
    if (MSG0) printf("p-part taken has order %ld\n",ppow);
    fflush(stdout);
}


/* ------------------------------------------------------------------
   init2() - Initialize hsz, calculate starting points
   ------------------------------------------------------------------ */

static void init2()

{
    int i;
    long l, m;
    FEL f;

    zmulrow(hsz,F_ZERO);
    for (i = 1; i <= (int) norb; ++i)
    {	l = orbsize[i-1];
	m = l / ppow;
	if (m * ppow == l)
	    f = zinv(zitof(m % prime));
	else
	    f = F_ZERO;
	zinsert(hsz,i,f);
    }

    /* Calculate one point of each orbit
       --------------------------------- */
    m = 1;
    for (i = 0; i < (int) norb; ++i)
    {	l = orbsize[i];
	orbsize[i] = m;
	m += l;
    }
    orbsize[norb] = degree + 1;
}



int main(int argc, char *argv[])

{
    int orbit;	/* Current orbit */
    long pt;
    long ee;
    long orb;
    FEL f;

    init(argc,argv);
    readdata();
    init2();

    if ((kondfile = zwritehdr(kondname,fl,norb,norb)) == NULL)
	errexit(-1,kondname);;
    for (orbit = 1; orbit <= (int) norb; ++orbit)
    {	zmulrow(m1,F_ZERO);	/* Clear row */
	for (pt = orbsize[orbit-1]; pt < orbsize[orbit]; ++pt)
	{
	    ee = ((long *)perm)[pt-1];
	    for (orb = 0; orbsize[orb] <= ee; ++orb);
	    f = zextract(m1,orb);
	    f = zadd(f,zextract(hsz,orb));
	    zinsert(m1,orb,f);
	}
	zwritevec(kondfile,m1,1);
    }
    fclose(kondfile);
    return (EXIT_OK);
}



