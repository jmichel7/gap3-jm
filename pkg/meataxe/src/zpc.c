/* ========================== C MeatAxe =============================
   zpc.c - Permutation chop.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zpc.c,v 1.2 1997/09/11 15:44:12 gap Exp $
 *
 * $Log: zpc.c,v $
 * Revision 1.2  1997/09/11 15:44:12  gap
 * New version 2.2.3. AH
 *
 * Revision 2.6  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.5  1994/07/07  06:57:14  mringe
 * Include files.
 *
 * Revision 2.4  1994/03/13  13:27:01  mringe
 * Maschinenunabhaengiges Format fuer Permutationen.
 *
 * Revision 2.3  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
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
 * Revision 1.9  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.8  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.7  1993/08/26  17:31:35  mringe
 * Helptext.
 *
 * Revision 1.6  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.5  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.4  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.3  1992/07/10  15:17:10  mringe
 * Argumente aus nicht-ANSI-prototypes entfernt.
 *
 * Revision 1.2  1992/06/01  08:08:43  mringe
 * CL warnings beseitigt.
 *
 * Revision 1.1  1992/05/26  07:56:18  mringe
 * Initial revision
 */

#include <string.h>
#include <stdlib.h>
#include "meataxe.h"


#define MAXGEN 20


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

long *seed;			/* Read from G3 */
long npoints;			/* Degree of input permutations */
long blksize;			/* Block size */
long nblocks;			/* Number of blocks */
long *perm[MAXGEN+1];		/* Generators */
int ngen;			/* Number of generators */
int nperm;			/* Number of permutations to chop */
long *orb;
long *num;
long start;
long orblen;

char *genname[MAXGEN+1];
char *subname[MAXGEN+1];
char *quotname[MAXGEN+1];
char *seedname;

int opt_b = 0;			/* Option -b (permute blocks) */
int opt_g = 0;			/* Option -g (additional generators) */


static char *helptext[] = {
"SYNTAX",
"    zpc [-b] <Perm1> <Perm2> <Seed> <S1> <S2> <Q1> <Q2>",
"    zpc [-b] -g <#Perm>[.<#Gen>] <Perm> <Seed> <Sub> <Quot>",
"",
"OPTIONS",
"    -b   Find action on block system",
"    -g   Multiple permutation mode. <#Perm> is the number of",
"         permutations to chop, <#Gen> is the number of generators.",
"",
"FILES",
"    <Perm1,2> i  Generators: permutations of equal degree",
"    <Seed>    i  Seed: One point or (with -b) one block",
"    <S1,2>    o  Action on the orbit or block system",
"    <Q1,2>    o  Action on the rermaining points",
"",
"If the -g option is used, actual file names are derived from <Perm>,",
"<Sub> and <Quot>, respectively, by appending `1', `2', ...",
"",
NULL};

static proginfo_t pinfo =
   { "zpc", "Permutation Chop", "$Revision: 1.2 $", helptext };




/* ------------------------------------------------------------------
   err() - Print error message and exit
   ------------------------------------------------------------------ */

static void err(int c)

{   fprintf(stderr,"ZPC ERROR - ");
    switch (c)
    {	
	case 'b':
	    fprintf(stderr,"BLOCK SIZE DOES NOT DIVIDE DEGREE\n");
	    break;
	case 'x':
	    fprintf(stderr,"NOT A BLOCK SYSTEM\n");
	    break;
    }
    exit(EXIT_ERR);
}



/* ------------------------------------------------------------------
   mkname()
   ------------------------------------------------------------------ */

static char *mkname(char *basename,int i)

{   
    char *c;

    c = malloc(strlen(basename)+10);
    sprintf(c,"%s%d",basename,i);
    return c;
}


/* ------------------------------------------------------------------
   setnperm()
   ------------------------------------------------------------------ */

static void setnperm(char *c)

{
    nperm = -1;
    opt_g = 1;
    if (sscanf(c,"%d.%d",&nperm,&ngen) == 1)
	ngen = nperm;
    if (nperm > MAXGEN || ngen > nperm || ngen < 2)
	errexit(ERR_OPTION,"-g");
    if (ngen != nperm && MSG1)
        printf("%d GENERATORS, %d PERMUTATIONS",ngen,nperm);
}


/* ------------------------------------------------------------------
   parseargs() - Process command line options and arguments
   ------------------------------------------------------------------ */

static void parseargs(int argc, char **argv)

{
    int i;

    ngen = 2;
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("Gb")) != OPT_END)
    {
	switch (i)
	{
	    case 'b': opt_b = 1; break;
	    case 'g': setnperm(opt_text_ptr); break;
	}
    }

    /* Get file names
       -------------- */
    if (opt_g)
    {
	if (argc - opt_ind != 4) errexit(ERR_NARGS,"zpc");
	seedname = argv[opt_ind+1];
	for (i = 1; i <= nperm; ++i)
	{	genname[i] = mkname(argv[opt_ind+0],i);
		subname[i] = mkname(argv[opt_ind+2],i);
		quotname[i] = mkname(argv[opt_ind+3],i);
	}
    }
    else
    {	nperm = 2;
	if (argc - opt_ind != 7) errexit(ERR_NARGS,"zpc");
	quotname[2] = argv[opt_ind+6];
	quotname[1] = argv[opt_ind+5];
	subname[2] = argv[opt_ind+4];
	subname[1] = argv[opt_ind+3];
	seedname = argv[opt_ind+2];
	genname[2] = argv[opt_ind+1];
	genname[1] = argv[opt_ind+0];
    }
}



/* ------------------------------------------------------------------
   init() - Initialize everything, read input files
   ------------------------------------------------------------------ */

static void init(int argc, char **argv)

{
    long fl, nop, nor2;
    long *m;
    int i;
    FILE *f;

    mtxinit();
    parseargs(argc, argv);



    /* Read the Permutation
       -------------------- */
    for (i = 1; i <= nperm; ++i)
    {
	if ((f = zreadhdr(genname[i],&fl,&nor2,&nop)) == NULL)
	    errexit(-1,genname[i]);
	if (fl != -1 || nop != 1)
    	    errexit(ERR_NOTPERM,genname[i]);
	if (i == 1)
	    npoints = nor2;
	else if (nor2 != npoints)
	{
    	    char buf[200];
	    sprintf(buf,"%s and %s",genname[1],genname[i]);
	    errexit(ERR_INCOMPAT,buf);
	}
    	m = (long *) malloc(sizeof(long)*npoints);
	if (m == NULL) errexit(ERR_NOMEM,"zpc: init()");
	if (zreadlong(f,m,npoints) != npoints)
	    errexit(-1,genname[i]);
        perm[i] = m - 1;
	fclose(f);
    }

    /* Seed
       ---- */
    if ((f = zreadhdr(seedname,&fl,&blksize,&nop)) == NULL)
	errexit(-1,seedname);
    if (!opt_b)
    {	blksize = 1;
	nblocks = npoints;
    }
    else
    {	nblocks = npoints / blksize;
	if (npoints % blksize != 0) err('b');
    }
    if (fl != -1) errexit(ERR_FILEFMT,seedname);

    seed = (long *) malloc(sizeof(long)*blksize);
    if (seed == NULL) errexit(ERR_NOMEM,"zpc");
    if (zreadlong(f,seed,blksize) != blksize)
	errexit(-1,seedname);
    fclose(f);

    /* Allocate tables
       --------------- */
    orb = (long *) malloc(((size_t)npoints+1) * sizeof(long));
    num = (long *) malloc(((size_t)npoints+1) * sizeof(long));
    for (i = 1; i <= (int) npoints; ++i)
    {
	orb[i] = 0;
	num[i] = 0;
    }
}


/* ------------------------------------------------------------------
   writeresult() - Calculate action on orbits and write output files.
   ------------------------------------------------------------------ */

static void writeresult()

{   
    long cosize;
    long *s, *p;
    int i;
    long k, im;
    FILE *f;

    cosize = npoints - orblen;
    if (cosize == 0)
	printf("Transitive on %ld %s\n",orblen/blksize,
		opt_b ? "blocks" : "points");
    else
	printf("Intransitive - 'sub' %ld  'quot' %ld\n",
		orblen/blksize,cosize/blksize);

    /* Calculate action on first orbit
       ------------------------------- */
    if (orblen % blksize != 0) err('x');
    s = (long*) malloc(((size_t)(orblen/blksize)+1) * sizeof(long));
    for (i = 1; i <= nperm; ++i)
    {	p = perm[i];
	for (k = 1; k <= orblen/blksize; ++k)
		s[k] = -1;
	for (k = 1; k <= npoints; ++k)
	{	long x, y;
			if (num[k] > orblen) continue;
			x = (num[k] - 1) / blksize + 1;
			y = (num[p[k]] - 1) / blksize + 1;
			if (s[x] == -1)
				s[x] = y;
			else
				if (s[x] != y) err('x');
	}
	if ((f=zwritehdr(subname[i],(long)-1,orblen/blksize,
			 (long)1)) == NULL)
	    errexit(-1,subname[i]);
	if (zwritelong(f,s+1,orblen/blksize) != orblen/blksize)
	    errexit(-1,subname[i]);
	fclose(f);
    }
    free(s);
    if (opt_b || cosize <= 0) return;

    /* Calculate action on other orbits (not with -b option)
       ----------------------------------------------------- */
    s = (long*) malloc(((size_t)cosize+1) * sizeof(long));
    for (i = 1; i <= nperm; ++i)
    {	p = perm[i];
	for (k = 1; k <= npoints; ++k)
	{
	    if (num[k] <= orblen) continue;
	    im = p[k];
	    s[num[k]-orblen] = num[im]-orblen;
	}
	if ((f = zwritehdr(quotname[i],(long)-1,cosize,(long)1))==NULL)
	    errexit(-1,quotname[i]);
	if (zwritelong(f,s+1,cosize) != cosize)
	    errexit(ERR_FILEWRITE,subname[i]);
	fclose(f);
    }
    free(s);
}


/* ------------------------------------------------------------------
   chop() - Make orbit
   ------------------------------------------------------------------ */

static void chop()

{	long level, i, k, newpt;
	int gen;

	/* Set up the orbit to contain only the seed point/block
	   ----------------------------------------------------- */
	for (i = 1; i <= blksize; ++i)
	{	orb[i] = seed[i-1];
		num[orb[i]] = i;
	}
	orblen = blksize;
	level = 1;

	/* Make orbit
	   ---------- */
	while (level <= orblen)
	{	/* Apply all generators
		   -------------------- */
		for (gen = 1; gen <= ngen; ++gen)
		{	if (num[perm[gen][orb[level]]] == 0)
				/* New point ? */
			{	for (i = 1; i <= blksize; ++i)
				{	newpt=perm[gen][orb[level+i-1]];
					num[newpt] = orblen+i;
					orb[orblen+i] = newpt;
				}
				orblen += blksize;
			}
		}
		level += blksize;
	}

	if (opt_b)	/* No `quotient' when permuting blocks */
		return;

	/* Renumber the remaining points
	   ----------------------------- */
	k = orblen;
	for (i = 1; i <= npoints; ++i)
	{	if (num[i] == 0)
			num[i] = ++k;
	}
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char *argv[])

{	init(argc,argv);
	chop();
	writeresult();
	return EXIT_OK;
}
