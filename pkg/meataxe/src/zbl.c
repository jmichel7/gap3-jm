/* ========================== C MeatAxe =============================
   zbl.c - Make matrix lower triangular (keeping bottom left).

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zbl.c,v 1.2 1997/09/11 15:43:37 gap Exp $
 *
 * $Log: zbl.c,v $
 * Revision 1.2  1997/09/11 15:43:37  gap
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
 * Revision 1.5  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.4  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.3  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.2  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.1  1992/05/26  18:24:49  mringe
 * Initial revision
 *
 */

#include "meataxe.h"



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static char *helptext[] = {
"SYNTAX",
"    zbl <Matrix> <Result>",
"",
"FILES",
"    <Matrix>   i  Any matrix",
"    <Result>   o  The result, a matrix of the same size",
NULL};

static proginfo_t pinfo =
   { "zbl", "Bottom Left of a Matrix", "$Revision: 1.2 $", helptext };



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{
    long fl;
    PTR m1;
    long nor, noc, i, j;
    char *iname, *oname;
    FILE *inp, *out;

    mtxinit();

    /* Parse command line
       ------------------ */
    initargs(argc, argv, &pinfo);
    while (zgetopt("") != OPT_END);
    if ((i = opt_ind) != argc-2) errexit(ERR_NARGS,"zbl");
    iname = argv[opt_ind];
    oname = argv[opt_ind+1];
    if ((inp = zreadhdr(iname,&fl,&nor,&noc)) == NULL)
    {
	perror(iname);
	errexit(-1,iname);
    }
    if (fl < 1) errexit(ERR_NOTMATRIX,iname);
    zsetfield(fl);
    zsetlen(noc);
    m1 = zalloc((long)1);
    if ((out = zwritehdr(oname,fl,nor,noc)) == NULL)
    {
	perror(oname);
	errexit(-1,oname);
    }
    for (i = 1; i <= nor; ++i)
    {
	zreadvec(inp,m1,1);
	for (j = i + 1; j <= noc; ++j)
	    zinsert(m1,j,F_ZERO);
	zwritevec(out,m1,1);
    }
    fclose(inp);
    fclose(out);
    return (EXIT_OK);
}
