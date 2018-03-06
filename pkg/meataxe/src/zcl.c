/* ========================== C MeatAxe =============================
   zcl.c - Clean matrix.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zcl.c,v 1.2 1997/09/11 15:43:40 gap Exp $
 *
 * $Log: zcl.c,v $
 * Revision 1.2  1997/09/11 15:43:40  gap
 * New version 2.2.3. AH
 *
 * Revision 2.4  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.3  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
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
 * Revision 1.10  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.9  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.8  1993/08/10  14:51:42  mringe
 * Include string.h
 *
 * Revision 1.7  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.6  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.5  1993/07/19  13:21:54  mringe
 * help(), message library.
 *
 * Revision 1.4  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.3  1992/07/15  09:25:55  mringe
 * Some minor changes.
 *
 * Revision 1.2  1992/05/26  18:30:18  mringe
 * Schreibfehler im USAGE korrigiert.
 *
 * Revision 1.1  1992/05/26  18:26:16  mringe
 * Initial revision
 *
 */

#include <string.h>
#include <stdlib.h>
#include "meataxe.h"




/* ------------------------------------------------------------------
   Variables
   ------------------------------------------------------------------ */

static char *helptext[] = {
"SYNTAX",
"    zcl <Subsp> <Mat> <Cleaned mat> <Ops>",
"",
"FILES",
"    <Subsp>    i  The subspace to clean with",
"    <Mat>      i  The matrix to be cleaned.",
"    <Cleaned>  o  The cleaned matrix",
"    <Ops>      o  Contains the operations that were performed, i.e.",
"                  <Mat> = <Ops> * <Subsp> + <Cleaned>",
NULL};

static proginfo_t pinfo =
   { "zcl", "Clean Matrix", "$Revision: 1.2 $", helptext };


static char *matname, *subname, *clname, *opname;



int main(argc, argv)
int argc;
char *argv[];

{
    PTR mat, row, op;
    long fl, fl2;
    long nor1, noc, nor2, noc2;
    long j;
    long *piv;
    int i;
    char buf[200];
    FILE *f, *fcl, *fop;


    /* Parse command line
       ------------------ */
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("")) != OPT_END)
    {
    }

    if (opt_ind != argc-4) errexit(ERR_NARGS,"zcl");
    subname = argv[opt_ind];
    matname = argv[opt_ind+1];
    clname = argv[opt_ind+2];
    opname = argv[opt_ind+3];

    /* Open files and initialize
       ------------------------- */
    if ((f = zreadhdr(subname,&fl,&nor1,&noc)) == NULL)
	errexit(-1,subname);
    if (fl <= 1) errexit(ERR_NOTMATRIX,subname);
    piv = (long *) malloc((size_t)(nor1+1) * sizeof(long));
    if (piv == NULL) errexit(ERR_NOMEM,"zcl");
    zsetfield(fl);
    zsetlen(noc);
    mat = zalloc(nor1);
    row = zalloc((long)1);
    zreadvec(f,mat,nor1);
    fclose(f);
    if (zmkpivot(mat,nor1,piv) == -1) /* Make pivot table */
	errexit(-1,subname);
    if ((f = zreadhdr(matname,&fl2,&nor2,&noc2)) == NULL)
	errexit(-1,matname);
    if (noc2 != noc || fl2 != fl)
    {
	sprintf(buf,"%s and %s",subname,matname);
        errexit(ERR_INCOMPAT,buf);
    }
    zsetlen(nor1);
    op = zalloc((long)1);
    zmulrow(op,F_ZERO);
    if ((fcl = zwritehdr(clname,fl,nor2,noc)) == NULL)
	errexit(-1,clname);
    if ((fop = zwritehdr(opname,fl,nor2,nor1)) == NULL)
	errexit(-1,opname);

    /* Process matrix row by row
       ------------------------- */
    for (j = 1; j <= nor2; ++j)
    {
	zsetlen(noc);
	zreadvec(f,row,1);		/* Read next row */
	zcleanrow2(row,mat,nor1,piv,op);/* Clean with matrix */
	zwritevec(fcl,row,1);		/* Write cleaned row */
	zsetlen(nor1);
	zwritevec(fop,op,1);		/* Write operation */
	zmulrow(op,F_ZERO);
    }
    fclose(f);
    fclose(fop);
    fclose(fcl);
    return (EXIT_OK);
}

