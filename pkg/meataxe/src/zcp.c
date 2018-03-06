/* ========================== C MeatAxe =============================
   zcp.c -  Characteristic polynomial of a matrix

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zcp.c,v 1.2 1997/09/11 15:43:41 gap Exp $
 *
 * $Log: zcp.c,v $
 * Revision 1.2  1997/09/11 15:43:41  gap
 * New version 2.2.3. AH
 *
 * Revision 2.11  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.10  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.9  1994/07/23  17:00:57  mringe
 * Minimalpolynom.
 *
 * Revision 2.8  1994/07/21  11:17:11  mringe
 * *** empty log message ***
 *
 * Revision 2.7  1994/07/20  12:17:04  mringe
 * Minimalpolynome
 *
 * Revision 2.6  1994/05/04  11:43:21  mringe
 * charpol() ind charpolfactor() umbenannt.
 *
 * Revision 2.5  1994/04/07  21:11:57  mringe
 * Hauptteil nach charpol.c verlagtert.
 *
 * Revision 2.4  1994/03/14  12:59:28  mringe
 * Option -f.
 *
 * Revision 2.3  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.2  1993/10/26  10:47:35  mringe
 * Compiler Warnings.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.15  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.14  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.13  1993/08/10  14:51:42  mringe
 * Include string.h
 *
 * Revision 1.12  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.11  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.10  1993/07/23  11:40:20  mringe
 * GAP-Output: Benutze MeatAxe.CharPol.
 *
 * Revision 1.9  1993/07/17  09:14:30  mringe
 * Bugs korrigiert.
 *
 * Revision 1.8  1993/07/17  08:59:05  mringe
 * Benutze message library. Optionen: -G, -Q, -V.
 *
 * Revision 1.7  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.6  1993/05/21  07:55:49  mringe
 * Newlines auf dem Hilfstext entfernt.
 *
 * Revision 1.5  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.4  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.3  1993/01/07  15:02:51  mringe
 * GAP-Output
 *
 */

#include <string.h>
#include <stdlib.h>
#include "meataxe.h"




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

char *fname;			/* Matrix file name */
matrix_t *mat;			/* The matrix */
int opt_G = 0;			/* GAP output */
int opt_f = 0;			/* Factorization */
int opt_m = 0;			/* Minimal polynomial */
fpoly_t *cpol;			/* Char. polynomial (with -f) */

static char *helptext[] = {
    "SYNTAX",
    "    zcp [-GQVfm] <File>",
    "",
    "OPTIONS",
    "    -G   GAP output",
    "    -Q   Quiet, no messages",
    "    -V   Verbose, more messages",
    "    -m   Calculate the minimal polynomial",
    "    -f   Factorize the polynomial",
    "",
    "FILES",
    "    <File>    i  A square matrix",
    NULL
    };

static proginfo_t pinfo =
   { "zcp", "Characteristic and Minimal Polynomial",
     "$Revision: 1.2 $", helptext };




/* ------------------------------------------------------------------
   init()
   ------------------------------------------------------------------ */

static void init(argc, argv)
int argc;
char *argv[];

{	
    int i;
    FILE *f;

    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("Gmf")) != OPT_END)
    {
	switch (i)
	{
	    case 'G': opt_G = 1; msg_level = -100; break;
	    case 'f': opt_f = 1; break;
	    case 'm': opt_m = 1; break;
	}
    }
    if (opt_ind != argc - 1) errexit(ERR_NARGS,pinfo.name);
    fname = argv[opt_ind];

    /* Read the matrix
       --------------- */
    f = SFOpen(fname,FM_READ);
    if (f == NULL) errexit(-1,fname);
    mat = matread(f);
    fclose(f);
    if (mat == NULL) errexit(-1,fname);
    if (mat->nor != mat->noc) errexit(ERR_NOTSQUARE,fname);
}


/* ------------------------------------------------------------------
   write_init() - Called once before factors are written
   write_factor() - Write one factor
   write_end() - Called once after all factors have been written
   ------------------------------------------------------------------ */

static int first;

static void write_init()

{
    first = 1;
    if (opt_m)
    {
    	if (opt_G)
	    printf("MeatAxe.MinPol:=[\n");
    	else
	    printf("MINIMAL POLYNOMIAL:\n");
    }
    else
    {
    	if (opt_G)
	    printf("MeatAxe.CharPol:=[\n");
    	else
	    printf("CHARACTERISTIC POLYNOMIAL:\n");
    }
}


static void write_end()

{
    int i;

    if (opt_G)
	printf("];\n");
    else
    {
	if (opt_f)
	    for (i = 0; i < cpol->len; ++i)
    	    {
		printf("( ");
		polprint(NULL,cpol->p[i]);
		printf(" )^%ld\n",cpol->e[i]);
    	    }
    }
}



static void write_one(pol)
poly_t *pol;

{ 
    if (opt_G)
    {
	int i;
	if (!first)
	    printf(",\n");
	printf("[");
	for (i = 0; i < pol->deg; ++i)
	    printf("%s,",zftogap(pol->buf[i]));
	printf("%s]",zftogap(F_ONE));
    }
    else if (!opt_f)
    {
	polprint(NULL,pol);
	printf("\n");
    }
    else
    {
	fpoly_t *f = factorization(pol);
	fpolmul(cpol,f);
	fpolfree(f);
    }
    first = 0;
}




/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{	
    poly_t *p;

    init(argc,argv);
    if (opt_f) cpol = fpolalloc();
    write_init();
    if (opt_m)
    {
    	for (p = minpolfactor(mat); p != NULL; p = minpolfactor(NULL))
	    write_one(p);
    }
    else
    {
    	for (p = charpolfactor(mat); p != NULL; p = charpolfactor(NULL))
	    write_one(p);
    }
    write_end();
    return EXIT_OK;
}


