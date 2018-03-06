/* ========================== C MeatAxe =============================
   zad.c - Add two matrices.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zad.c,v 1.2 1997/09/11 15:43:34 gap Exp $
 *
 * $Log: zad.c,v $
 * Revision 1.2  1997/09/11 15:43:34  gap
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
 * Revision 1.12  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.11  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.10  1993/10/02  15:33:21  mringe
 * Neue File I/O-Library.
 *
 * Revision 1.9  1993/08/10  14:51:42  mringe
 * Include string.h
 *
 * Revision 1.8  1993/08/09  12:39:42  mringe
 * Dokumentation entfernt
 *
 * Revision 1.7  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.6  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.5  1993/06/29  06:44:57  mringe
 * Dokumentation eingefuegt.
 *
 * Revision 1.4  1993/05/13  07:25:28  mringe
 * Hilfe verbessert.
 *
 * Revision 1.3  1993/05/13  07:23:43  mringe
 * help() eingebaut, Prototypes korrigiert.
 *
 * Revision 1.2  1992/05/30  09:21:17  mringe
 * Matrizen zeilenweise einlesen und rausschreiben.
 *
 * Revision 1.1  1992/05/26  07:55:51  mringe
 * Initial revision
 *
 */


#include <string.h>
#include "meataxe.h"



/* ------------------------------------------------------------------
   Variables
   ------------------------------------------------------------------ */

static char *helptext[] = {
"SYNTAX",
"    zad <Mat1> <Mat2> <Result>",
"",
"FILES",
"    <Mat1>    i  First matrix",
"    <Mat2>    i  Second matrix. Must be of the same type as <Mat1>.",
"    <Result>  o  The result",
NULL};

static proginfo_t pinfo =
   { "zad", "Add Matrices", "$Revision: 1.2 $", helptext };


/* ------------------------------------------------------------------
   Function prototypes
   ------------------------------------------------------------------ */


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{
    char buf[200];
    long nor, noc, nor2, noc2, fl, fl2;
    FILE *file1, *file2, *file3;
    PTR m1, m2;
    register long i;

    mtxinit();

    /* Parse command line
       ------------------ */
    initargs(argc, argv, &pinfo);
    while (zgetopt("") != OPT_END);
    if ((i = opt_ind) != argc-3) errexit(ERR_BADUSAGE,"zad");

    /* Open input files
       ---------------- */
    if ((file1 = zreadhdr(argv[opt_ind],&fl,&nor,&noc)) == NULL)
    {
	perror(argv[opt_ind]);
	errexit(mtxerrno,argv[opt_ind]);
    }
    if (fl < 2) errexit(ERR_NOTMATRIX,argv[opt_ind]);
    if ((file2 = zreadhdr(argv[opt_ind+1],&fl2,&nor2,&noc2)) == NULL)
    {
	perror(argv[opt_ind+1]);
	errexit(ERR_FILEREAD,argv[opt_ind+1]);
    }
    if (fl2 < 2) errexit(ERR_NOTMATRIX,argv[opt_ind+1]);
    if (fl != fl2 || nor != nor2 || noc != noc2)
    {
	sprintf(buf,"%s and %s",argv[opt_ind],argv[opt_ind+1]);
    	errexit(ERR_INCOMPAT,buf);
    }
    zsetfield(fl);
    zsetlen(noc);
    file3 = zwritehdr(argv[opt_ind+2],fl,nor,noc);
    if (file3 == NULL)
    {
	perror(argv[opt_ind+2]);
	errexit(mtxerrno,argv[opt_ind+2]);
    }

    /* Add matrices row by row
       ----------------------- */
    m1 = zalloc((long)1);
    m2 = zalloc((long)1);
    for (i = nor; i > 0; --i)
    {
	if (zreadvec(file1,m1,1) != 1) errexit(-1,argv[1]);
	if (zreadvec(file2,m2,1) != 1) errexit(-1,argv[2]);
	zaddrow(m1,m2);
	if (zwritevec(file3,m1,1) != 1) errexit(-1,argv[3]);
    }
    fclose(file1);
    fclose(file2);
    fclose(file3);
    return EXIT_OK;
}

