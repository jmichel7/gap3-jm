/* ========================== C MeatAxe =============================
   zfr.c - Apply the Frobenius automorphism to a matrix.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zfr.c,v 1.2 1997/09/11 15:43:51 gap Exp $
 *
 * $Log: zfr.c,v $
 * Revision 1.2  1997/09/11 15:43:51  gap
 * New version 2.2.3. AH
 *
 * Revision 2.3  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.2  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
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
 * Revision 1.5  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.4  1993/08/05  15:48:29  mringe
 * help() eingebaut.
 *
 * Revision 1.3  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.2  1992/07/13  12:59:10  mringe
 * Falsche Zeilenzahl beim Lesen korrigiert
 *
 * Revision 1.1  1992/05/26  17:39:08  mringe
 * Initial revision
 *
 */


#include "meataxe.h"


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static char *helptext[] = {
"SYNTAX",
"    zfr [-QV] <Matrix> <Result>",
"",
"OPTIONS",
"    -Q   Quiet, no messages",
"    -V   Verbose, more messages",
"",
"FILES",
"    <Matrix>    i  A matrix",
"    <Result>    o  The result, a matrix of the same size",
"",
NULL};
static proginfo_t pinfo =
   { "zfr", "Frobenius Map", "$Revision: 1.2 $", helptext };

char *iname, *oname;
FILE *ifile, *ofile;


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{
    long fl;			/* Field */
    PTR m1;				/* One row of the matrix */
    long nor, noc;			/* Number of rows and columns */
    long i, k, n, p;
    FEL f1, f2;

    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while (zgetopt("") != OPT_END);
    if (argc - opt_ind != 2) errexit(ERR_NARGS,"zfr");
    iname = argv[opt_ind];
    oname = argv[opt_ind+1];

    /* Open the input file
       ------------------- */
    if ((ifile = zreadhdr(iname,&fl,&nor,&noc)) == NULL)
	errexit(-1,iname);
    if (fl < 2) errexit(ERR_NOTMATRIX,iname);
    zsetfield(fl); zsetlen(noc);
    p = zchar;
    MESSAGE(0,("CHARACTERISTIC IS %ld\n",p));

    /* Open output file, allocate memory
       --------------------------------- */
    if ((ofile = zwritehdr(oname,fl,nor,noc)) == NULL)
	errexit(-1,oname);
    m1 = zalloc((long)1);

    /* Apply the frobenius map to each entry
       ------------------------------------- */
    for (i = 1; i <= nor; ++i)
    {
	zreadvec(ifile,m1,(long) 1);
	for (k = 1; k <= noc; ++k)
	{
	    f1 = zextract(m1,k);
	    f2 = f1;
	    for (n = p - 1; n != 0; --n)
		f2 = zmul(f1,f2);
	    zinsert(m1,k,f2);
	}
	zwritevec(ofile,m1,(long) 1);
    }

    fclose(ifile);
    fclose(ofile);
    return (EXIT_OK);
}


