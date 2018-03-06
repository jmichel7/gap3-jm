/* ========================== C MeatAxe =============================
   zqt.c - Projection on quotient.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zqt.c,v 1.2 1997/09/11 15:44:20 gap Exp $
 *
 * $Log: zqt.c,v $
 * Revision 1.2  1997/09/11 15:44:20  gap
 * New version 2.2.3. AH
 *
 * Revision 2.5  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.4  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.3  1994/05/10  16:28:30  mringe
 * Neue Versionen von zquot() und zquotop().
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
 * Revision 1.11  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.10  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.9  1993/08/10  14:51:42  mringe
 * Include string.h
 *
 * Revision 1.8  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.7  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.6  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.5  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.4  1993/01/28  21:52:18  mringe
 * Use yyy library.
 *
 * Revision 1.3  1993/01/28  15:30:25  mringe
 * Benutze yyy library.
 *
 * Revision 1.2  1992/06/01  08:10:11  mringe
 * CL warnings beseitigt.
 *
 * Revision 1.1  1992/05/26  07:56:29  mringe
 * Initial revision
 *
 */

#include <string.h>
#include <stdlib.h>
#include "meataxe.h"




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static int opt_i = 0;
static PTR subsp, mat, tmp, quot;
static long fl, nor, noc, sdim, qdim;
static char *mname, *sname, *oname;
static FILE *mfile, *sfile, *ofile;


static char *helptext[] = {
"SYNTAX",
"    zqt [-i] [<Subsp> <Matrix> <Quotient>]",
"",
"OPTIONS",
"    -i   Take only insignificant rows of <Matrix>. <Quotient> will be",
"         the action of <Matrix> on the quotient by <Subspace>.",
"",
"FILES",
"    <Subsp>    i  The invariant subspace, in semi-echelon form",
"    <Matrix>   i  The matrix, must have the same number of columns",
"    <Quotient> o  Insignificant columns of <Matrix>, cleaned with <Subspace>",
"",
NULL};

static proginfo_t pinfo =
   { "zqt", "Clean And Quotient", "$Revision: 1.2 $", helptext };


/* ------------------------------------------------------------------
   err() - Print error message and exit
   ------------------------------------------------------------------ */

static void err(c)
int c;

{
    fprintf(stderr,"ZQT ERROR - ");
    switch (c)
    {
	case 'm':
	    fprintf(stderr,"MATRICES NOT COMPATIBLE\n");
	    break;
	case 'e':
	    fprintf(stderr,"SUBSPACE NOT IN ECHELON FORM\n");
	    break;
	default:
	    fprintf(stderr,"UNKNOWN ERROR `%c'\n",c);
	    break;
    }
    exit(1);
}


/* ------------------------------------------------------------------
   init() - Process command line options and arguments
   ------------------------------------------------------------------ */

static void init(argc, argv)
int argc;
char *argv[];

{   int i;

    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("i")) != OPT_END)
    {
	switch (i)
	{
	    case 'i': opt_i = 1; break;
	}
    }
    if (opt_ind != argc-3) errexit(ERR_NARGS,"zqt");
    sname = argv[opt_ind];
    mname = argv[opt_ind+1];
    oname = argv[opt_ind+2];
}


/* ------------------------------------------------------------------
   readfiles() - Read the input files, allocate tables
   ------------------------------------------------------------------ */

void readfiles()

{
    long fls, nors, nocs, flm, norm, nocm;

    if ((sfile = zreadhdr(sname,&fls,&nors,&nocs)) == NULL)
	errexit(-1,sname);
    if ((mfile = zreadhdr(mname,&flm,&norm,&nocm)) == NULL)
	errexit(-1,mname);
    if (nocs != nocm || fls != flm) err('m');
    if (opt_i && norm != nocm) errexit(ERR_NOTSQUARE,mname);
    zsetfield(fls);
    zsetlen(nocs);
    subsp = zalloc(nors);
    zreadvec(sfile,subsp,nors);
    mat = zalloc(norm);
    tmp = zalloc((long)1);
    zreadvec(mfile,mat,norm);
    fclose(sfile);
    fclose(mfile);
    fl = flm;
    noc = nocm;
    nor = norm;
    sdim = nors;
    qdim = noc - nors;
}


/* ------------------------------------------------------------------
   doit()
   ------------------------------------------------------------------ */
void doit()

{
    long *piv = (long*)malloc((size_t)(sdim+1)*sizeof(long));
    zsetfield(fl); zsetlen(noc);
    if (zmkpivot(subsp,sdim,piv) != 0) err('e');
    zquotinit(subsp,sdim,piv);
    if (opt_i)
    {
        zsetlen(qdim);
	quot = zalloc(qdim);
	zquotop(mat,quot);
    }
    else
    {
	zsetlen(sdim);
	quot = zalloc(nor);
	zquot(mat,nor,quot);
    }
}


/* ------------------------------------------------------------------
   writefile()
   ------------------------------------------------------------------ */

void writefile()

{
    long resnor = opt_i ? qdim : nor;
	
    zsetlen(qdim);
    if ((ofile = zwritehdr(oname,fl,resnor,qdim)) == NULL)
    	errexit(-1,oname);
    zwritevec(ofile,quot,resnor);
    fclose(ofile);
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{
    init(argc, argv);
    readfiles();
    doit();
    writefile();
    return 0;
}

