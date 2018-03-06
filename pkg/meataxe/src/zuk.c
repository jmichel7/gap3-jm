/* ========================== C MeatAxe =============================
   zuk - Unkondense

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zuk.c,v 1.2 1997/09/11 15:44:39 gap Exp $
 *
 * $Log: zuk.c,v $
 * Revision 1.2  1997/09/11 15:44:39  gap
 * New version 2.2.3. AH
 *
 * Revision 2.10  1995/06/29  16:28:53  mringe
 * Kosmetik
 *
 * Revision 2.9  1995/06/16  14:16:52  mringe
 * Lies Orbitlaengen im neuen FOrmat.
 *
 * Revision 2.8  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.7  1994/09/08  07:51:12  mringe
 * Bug behoben. Offenbar hat bis jetzt niemand dieses Programm benutzt...
 *
 * Revision 2.6  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.5  1994/03/13  13:27:01  mringe
 * Maschinenunabhaengiges Format fuer Permutationen.
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
 * Revision 1.7  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.6  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.5  1993/08/10  14:51:42  mringe
 * Include string.h
 *
 * Revision 1.4  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.3  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.2  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.1  1993/01/08  14:40:20  mringe
 * Initial revision
 *
 * Revision 1.1  1993/01/08  14:29:30  mringe
 * Initial revision
 *
 */


#include <stdlib.h>
#include <string.h>
#include "meataxe.h"





/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static char *helptext[] ={
"SYNTAX",
"    zuk <Vectors> <Orbits> <Result>",
"",
"FILES",
"    <Vectors>    i  A matrix (#columns = #orbits)",
"    <Orbits>     i  An orbit sizes file",
"    <Result>     o  Unkondensed vectors, a matrix",
NULL};

static proginfo_t pinfo =
   { "zuk", "Unkondense Vectors", "$Revision: 1.2 $", helptext };

static char *vecname, *orbname, *resname;


/* ------------------------------------------------------------------
   err() - Exit with error message
   ------------------------------------------------------------------ */

static void err(ch)
int ch;

{
    fprintf(stderr,"ZUK ERROR - ");
    switch (ch)
    {
	case '=':
	    fprintf(stderr,"MATRIX AND ORBIT SIZES ARE INCOMPATIBLE\n");
	    break;
    }
    exit(EXIT_ERR);
}




/* ------------------------------------------------------------------
   unkondense()
   ------------------------------------------------------------------ */

static void unkondense()

{
    long fl, nor, noc, fl2, norb, noc2, noc3;
    PTR inp, perm, out;
    long *os, *x;
    long i;
    FILE *vecfile, *orbfile, *resfile;

    /* Open the vector file
       -------------------- */
    if ((vecfile = zreadhdr(vecname,&fl,&nor,&noc)) == NULL)
	errexit(-1,vecname);
    if (fl < 2) errexit(ERR_NOTMATRIX,vecname);
    zsetfield(fl);
    zsetlen(noc);
    inp = zalloc((long)1);

    /* Read the orbit sizes
       -------------------- */
    if ((orbfile = zreadhdr(orbname,&fl2,&noc2,&norb)) == NULL)
	errexit(-1,orbname);
    if (fl2 != -T_IMAT || noc2 != 1)
	errexit(ERR_FILEFMT,orbname);
    if (norb != noc) err('=');
    perm = (PTR) malloc(sizeof(long)*norb);
    if (perm == NULL) errexit(ERR_NOMEM,pinfo.name);
    if (zreadlong(orbfile,(long *)perm,norb) != norb)
	errexit(ERR_FILEREAD,orbname);
    fclose(orbfile);
    os = (long *)perm - 1;

    /* Calculate length of output vectors
       ---------------------------------- */
    x = os + 1;
    noc3 = 0;
    for (i = 1; i <= norb; ++i)
    {
	noc3 += *x;
	++x;
    }
    MESSAGE(0,("OUTPUT IS %ld x %ld\n",nor,noc3));
    zsetlen(noc3);
    out = zalloc((long)1);
    if ((resfile = zwritehdr(resname,fl,nor,noc3)) == NULL)
	errexit(-1,resname);

    /* Unkondense row by row
       --------------------- */
    for (i = 1; i <= nor; ++i)
    {
	long k, outpos = 1;

	zsetlen(noc);
	zreadvec(vecfile,inp,1);
	zsetlen(noc3);
	zmulrow(out,F_ZERO);
	x = os + 1;
	for (k = 1; k <= norb; ++k)
	{
	    FEL f;
	    long l;

	    zsetlen(noc);
	    f = zextract(inp,k);
	    zsetlen(noc3);
	    for (l = 1; l <= *x; ++l)
		zinsert(out,outpos++,f);  
	    ++x;
	}
	zsetlen(noc3);
	zwritevec(resfile,out,1);
    }
    fclose(resfile);
    fclose(vecfile);
}




/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];


{
    mtxinit();

    /* Parse command line
       ------------------ */
    initargs(argc, argv, &pinfo);
    while (zgetopt("") != OPT_END);
    if (opt_ind != argc-3) errexit(ERR_NARGS,pinfo.name);
    vecname = argv[opt_ind];
    orbname = argv[opt_ind+1];
    resname = argv[opt_ind+2];
    unkondense();
    return (EXIT_OK);
}


