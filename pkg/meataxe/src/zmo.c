/* ========================== C MeatAxe =============================
   zmo.c - Make orbits under permutations

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zmo.c,v 1.2 1997/09/11 15:44:01 gap Exp $
 *
 * $Log: zmo.c,v $
 * Revision 1.2  1997/09/11 15:44:01  gap
 * New version 2.2.3. AH
 *
 * Revision 2.17  1995/02/09  14:03:44  mringe
 * Neues Format fuer Bahnlaengendatei
 *
 * Revision 2.16  1994/12/21  13:40:48  mringe
 * Tippfehler
 *
 * Revision 2.15  1994/12/21  11:06:13  mringe
 * Benutze ErrorExit().
 *
 * Revision 2.14  1994/12/14  07:10:29  mringe
 * ANSI C
 *
 * Revision 2.13  1994/08/23  14:18:49  mringe
 * Syntax jetzt wie bei zsp.
 *
 * Revision 2.12  1994/08/22  08:41:23  mringe
 * Bug behoben.
 *
 * Revision 2.11  1994/08/22  08:28:35  mringe
 * Neu: Option -s (Set seed).
 *
 * Revision 2.10  1994/03/12  13:03:53  mringe
 * Verschiedene Bugs beseitigt.
 *
 * Revision 2.9  1994/03/09  14:03:42  mringe
 * Warnung, wenn bei -g > 1 permutation in einem File steht.
 *
 * Revision 2.8  1994/03/09  13:53:30  mringe
 * Option -g
 *
 * Revision 2.7  1994/03/09  12:14:33  mringe
 * Fehler im Sizes-File behoben.
 *
 * Revision 2.6  1994/02/15  10:28:33  mringe
 * MSDOS_BCC entfernt.
 *
 * Revision 2.5  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.4  1993/10/26  10:47:35  mringe
 * Compiler Warnings.
 *
 * Revision 2.3  1993/10/21  21:57:35  mringe
 * Permutationen.
 *
 * Revision 2.2  1993/10/21  08:57:15  mringe
 * Compiler-Warnings.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.8  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.7  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.6  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.5  1993/08/06  12:53:45  mringe
 * help() eingebaut.
 *
 * Revision 1.4  1993/07/23  13:46:27  mringe
 * OS-Symbole neu (SYS_xxx)
 *
 * Revision 1.3  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.2  1992/06/27  17:20:11  mringe
 * Sagt jetzt, wenn es mehr als 200 verschiedene
 * Bahnl"angen gibt.
 *
 * Revision 1.1  1992/05/26  17:39:25  mringe
 * Initial revision
 *
 */


#include <stdlib.h>
#include <string.h>
#include "meataxe.h"




#define MAXPERM 50		/* Max. number of permutations */
#define MAXLEV 100000
#define MAXNORB 20000


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

long *perm[MAXPERM];			/* Permutations */
long sstt[201], ors[MAXNORB+1];
long fl, nor;
char *orbname, *sizname;
int nperm = 0;

static char *helptext[] = {
"SYNTAX",
"    zmo [<Options>] <Perm1> <Perm2> <Orbits> <OrbSz>",
"    zmo [<Options>] [-g <#Perms>] <Perms> <Orbits> <OrbSz>",
"",
"OPTIONS",
"    -Q         Quiet, no messages",
"    -V         Verbose output",
"    -g         Set number of permutations.",
"    -s <Seed>  Set seed (default:1)",
"",
"FILES",
"    <Perm1,2>  i  Permutations",
"    <Orbits>   o  Orbits under the permutations (Permutation format)",
"    <OrbSz>    o  Orbit sizes (Permutation format)",
"",
NULL};

static proginfo_t pinfo =
   { "zmo", "Make Orbits", "$Revision: 1.2 $", helptext };




/* ------------------------------------------------------------------
   err() - Print error message and exit
   ------------------------------------------------------------------ */

static void err(int c)

{
    fprintf(stderr,"ZMO ERROR - ");
    switch (c)
    {
	case 'o':
	    fprintf(stderr,"MORE THAN %ld ORBITS\n",(long)MAXNORB);
	    break;
	}
	exit(EXIT_ERR);
}


/* ------------------------------------------------------------------
   ReadPerm() - Read a permutation
   ------------------------------------------------------------------ */

static void ReadPerm(int i, char *fn)

{
    FILE *f;
    static char n1[200];
    static int first = 1;
    static long nor1;
    long nop;

    if ((f = zreadhdr(fn,&fl,&nor,&nop)) == NULL)
	ErrorExit("%s: %E",fn,mtxerrno);
    if (fl != -1)
	ErrorExit("%s: %E",fn,ERR_NOTPERM);
    MESSAGE(1,("Reading %s: Degree = %ld\n",fn,nor));
    if (first)
    {
	first = 0;
	strcpy(n1,fn);
	nor1 = nor;
    }
    else if (nor != nor1)
	ErrorExit("%s and %s: %E",fn,n1,ERR_INCOMPAT);
    if (nop > 1)
	fprintf(stderr,"%s: Only 1 of %ld permutations read\n",fn,nop);
    
    perm[i] = NALLOC(long,nor);
    if (perm[i] == NULL)
    	ErrorExit("%s: %E",pinfo.name,ERR_NOMEM);
    if (zreadlong(f,perm[i],nor) != nor)
	ErrorExit("%s: %E",fn,ERR_FILEREAD);
    fclose(f);
    perm[i] = perm[i] - 1;
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char *argv[])

{
    long i, seed = 1, nout, siz, dis, norb;
    long *stk, *mout;
    long f1, lev;
    FILE *f;

    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("g:s:")) != OPT_END)
    {
	switch (i)
	{
	    case 'g':
		nperm = (int) getint();
		if (nperm < 1 || nperm > MAXPERM || *opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-g");
		break;
	    case 's':
		seed = getint();
		if (seed == GETINT_ERR || *opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-s");
		break;
	}
    }

    if (nperm == 0)	/* Two generators */
    {
    	if (opt_ind != argc-4) errexit(ERR_NARGS,pinfo.name);
	nperm = 2;
	for (i = 0; i < 2; ++i)
	    ReadPerm(i,argv[opt_ind+i]);
    }
    else		/* N generators */
    {
	char n1[200];
    	if (opt_ind != argc-3) errexit(ERR_NARGS,pinfo.name);
	for (i = 0; i < nperm; ++i)
	{
	    sprintf(n1,"%s.%ld",argv[opt_ind],i+1);
	    ReadPerm(i,n1);
	}
    }
    orbname = argv[argc-2];
    sizname = argv[argc-1];

    if (seed < 0 || seed > nor)
    {
	char buf[100];
	sprintf(buf,"Seed is out of range (1-%ld)",nor);
	errexit(-1,buf);
    }

    stk = NALLOC(long,MAXLEV);
    mout = NALLOC(long,250);
    if (stk == NULL || mout == NULL) errexit(ERR_NOMEM,pinfo.name);

    if ((f = zwritehdr(orbname,fl,nor,(long)1)) == NULL)
	errexit(-1,orbname);
    norb = 0;
    nout = 0;
    while (seed <= nor)
    {
	if (norb >= MAXNORB) err('o');
	f1 = perm[0][seed];
	if (f1 < 0) errexit(-1,"internal error");
	perm[0][seed] = -f1;	/* Mark as done */
	siz = 0;
	stk[lev = 1] = seed;	/* Put first point on stack */

	while (lev > 0)		/* Stack not empty */
	{
	    f1 = stk[lev--];	/* Take from stack */
	    if (nout == 250)	/* If buffer full, write out */
	    {	
		if (zwritelong(f,mout,250) != 250)
		    errexit(ERR_FILEWRITE,orbname);
		nout = 0;
	    }
	    MESSAGE(2,("mout[%ld]=%ld\n",(long)nout,(long)f1));
	    mout[nout++] = f1;	/* Put into buffer */
	    ++siz;

	    /* Apply all permutations and collect new points
	       --------------------------------------------- */
	    for (i = 0; i < nperm; ++i)
	    {	
	    	register long f3, f5;
		f3 = perm[i][f1];
		if (f3 < 0) f3 = -f3;
		f5 = perm[0][f3];
		if (f5 <= 0) continue;	/* Already done, ignore */
		perm[0][f3] = -f5;	/* Mark as done */
		if (lev >= MAXLEV - 1)
		    errexit(ERR_NOMEM,"zmo (stack check)");
		stk[++lev] = f3;	/* Put on stack */
	    }
	}
	ors[norb++] = siz;

	/* find next unused point */
	for (seed = 1; seed <= nor && perm[0][seed] < 0; ++seed);
    }

    /* Flush buffer and close file
       --------------------------- */
    if (zwritelong(f,mout,nout) != nout)
   	 errexit(ERR_FILEWRITE,orbname);
    fclose(f);
 

    /* Write the orbit sizes file
       -------------------------- */
    MESSAGE(1,("Writing orbit sizes to %s\n",sizname));
    if ((f = zwritehdr(sizname,-T_IMAT,1,norb)) == NULL)
	errexit(-1,sizname);
    if (zwritelong(f,ors,norb) != norb)
   	 errexit(ERR_FILEWRITE,sizname);
    fclose(f);

    /* Print orbit sizes
       ----------------- */
    if (msg_level >= 0)
    {
	dis = 0;
	for (i = 0; dis < 200 && i < norb; ++i)
	{
	    int k;
	    siz = ors[i];

	    for (k = 0; k < dis && ors[k] != siz; ++k);
	    if (k < dis) ++sstt[k];
	    else if (dis++ < 200)
	    {
		sstt[k] = 1;
		ors[k] = siz;
	    }
	}
	if (dis > 200)
	    printf("More than 200 distinct orbit sizes\n");
	else if (dis > 15)
	    printf("%ld distinct orbit sizes\n", (long) dis);
	else
	    for (i = 0; i < dis; ++i)
		printf("%6ld ORBIT%c OF SIZE %6ld\n",
		    sstt[i],sstt[i] > 1 ? 'S':' ',ors[i]);
    }
    return (EXIT_OK);
}


