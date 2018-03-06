/* ========================== C MeatAxe =============================
   ziv.c - Invert a matrix or permutation.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: ziv.c,v 1.2 1997/09/11 15:43:57 gap Exp $
 *
 * $Log: ziv.c,v $
 * Revision 1.2  1997/09/11 15:43:57  gap
 * New version 2.2.3. AH
 *
 * Revision 2.7  1995/02/09  14:04:19  mringe
 * ANSI C
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
 * Revision 1.4  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.3  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.2  1992/07/15  09:27:50  mringe
 * Some minor changes.
 *
 * Revision 1.1  1992/05/26  07:56:07  mringe
 * Initial revision
 *
 */


#include <stdlib.h>
#include "meataxe.h"



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static char *iname, *oname;
static FILE *ifile, *ofile;
static long fl;
static long nor, noc;

static char *helptext[] = {
"SYNTAX",
"    ziv <File> <Result>",
"",
"FILES",
"    <File>    i  A square matrix or permutation",
"    <Result>  o  The result",
NULL};

static proginfo_t pinfo =
   { "ziv", "Invert Matrix Or Permutation",
     "$Revision: 1.2 $", helptext };



/* ------------------------------------------------------------------
   invperm() - Invert a permutation.
   ------------------------------------------------------------------ */

void invperm()

{
    long *p1, *p2;
    long i;

    if (noc > 1)
    {   fprintf(stderr,
	     "ZIV WARNING - ONLY 1 OF %ld ELEMENTS INVERTED\n",noc);
    }
    p1 = (long *) malloc(sizeof(long) * nor);
    p2 = (long *) malloc(sizeof(long) * nor);
    if (p1 == NULL || p2 == NULL) errexit(ERR_NOMEM,"ziv");
    if (zreadlong(ifile,p1,nor) != nor)
	errexit(-1,iname);

    /* Calculate inverse in p2
       ----------------------- */
    --p1;
    --p2;
    for (i = 1; i <= nor; ++i)
	p2[p1[i]] = i;
    if ((ofile = zwritehdr(oname,fl,nor,(long)1)) == NULL)
	errexit(-1,oname);
    if (zwritelong(ofile,p2+1,nor) != nor)
	errexit(-1,oname);
    fclose(ofile);
}


/* ------------------------------------------------------------------
   invmat() - Invert a matrix.
   ------------------------------------------------------------------ */

void invmat()

{
    PTR mat, result;

    /* Read matrix
       ----------- */
    zsetfield(fl); zsetlen(noc);
    if (nor != noc) errexit(ERR_NOTSQUARE,iname);
    mat = zalloc(nor);
    result = zalloc(nor);
    if (zreadvec(ifile,mat,nor) != nor) errexit(-1,iname);

    /* Make inverse
       ------------ */
    if (zmatinv(mat,result)) errexit(-1,iname);

    /* Write result
       ------------ */
    if ((ofile = zwritehdr(oname,fl,nor,noc)) == NULL)
	errexit(-1,oname);
    zwritevec(ofile,result,nor);
    fclose(ofile);
}



int main(argc,argv)
int argc;
char *argv[];

{

    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while (zgetopt("") != OPT_END);
    if (opt_ind != argc-2) errexit(ERR_NARGS,"ziv");
    iname = argv[opt_ind];
    oname = argv[opt_ind+1];
    if ((ifile = zreadhdr(iname,&fl,&nor,&noc)) == NULL)
	errexit(-1,iname);
    if (fl < 0)
       	invperm();
    else
       	invmat();
    fclose(ifile);
    return (EXIT_OK);
}



