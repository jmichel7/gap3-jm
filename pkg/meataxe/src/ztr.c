/* ========================== C MeatAxe =============================
   ztr.c - Transpose a matrix.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: ztr.c,v 1.2 1997/09/11 15:44:38 gap Exp $
 *
 * $Log: ztr.c,v $
 * Revision 1.2  1997/09/11 15:44:38  gap
 * New version 2.2.3. AH
 *
 * Revision 2.4  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.3  1994/06/25  13:34:18  mringe
 * Benutze zextractcol().
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
 * Revision 1.8  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.7  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.6  1993/10/05  23:44:33  mringe
 * err() eliminiert.
 *
 * Revision 1.5  1993/08/06  15:11:49  mringe
 * Bug behoben.
 *
 * Revision 1.4  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.3  1993/08/05  15:52:49  mringe
 * File header
 *
 * Revision 1.2  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.1  1992/05/26  18:17:58  mringe
 * Initial revision
 *
 */


#include "meataxe.h"

/* ------------------------------------------------------------------
   Global Data
   ------------------------------------------------------------------ */

static char *iname, *oname;

static char *helptext[] = {
"SYNTAX",
"    ztr [-QV] <Mat> <Result>",
"",
"OPTIONS",
"    -Q    Quiet, no messages",
"    -V    Verbatim, more messages",
"",
"FILES",
"    <Mat>      i  Any matrix",
"    <Result>   o  The transposed matrix",
NULL};

static proginfo_t pinfo =
   { "ztr", "Transpose", "$Revision: 1.2 $", helptext };



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{
    long fl;
    PTR m1, m2;
    long nor, noc, j;
    FILE *f;


    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while (zgetopt("") != OPT_END);
    if (opt_ind != argc-2) errexit(ERR_NARGS,"ztr");
    iname = argv[opt_ind];
    oname = argv[opt_ind+1];

    /* Read matrix
       ----------- */
    if ((f = zreadhdr(iname,&fl,&nor,&noc)) == NULL)
	errexit(-1,iname);
    if (fl < 2) errexit(ERR_NOTMATRIX,argv[opt_ind]);
    zsetfield(fl);
    zsetlen(noc);
    m1 = zalloc(nor);
    zreadvec(f,m1,nor);
    fclose(f);

    /* Transpose
       --------- */
    zsetlen(nor);
    m2 = zalloc((long)1);
    if ((f = zwritehdr(oname,fl,noc,nor)) == NULL)
	errexit(-1,oname);
    for (j = 1; j <= noc; ++j)
    {
	zsetlen(nor);
	zmulrow(m2,F_ZERO);
	zsetlen(noc);
	zextractcol(m1,nor,j,m2);
	zsetlen(nor);
	zwritevec(f,m2,(long)1);
    }
    fclose(f);
    return (EXIT_OK);
}


