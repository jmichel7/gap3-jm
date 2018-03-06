/* ========================== C MeatAxe =============================
   zev - Calculate eigenvalues and multiplitcities.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zev.c,v 1.2 1997/09/11 15:43:48 gap Exp $
 *
 * $Log: zev.c,v $
 * Revision 1.2  1997/09/11 15:43:48  gap
 * New version 2.2.3. AH
 *
 * Revision 2.6  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.5  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.4  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.3  1994/07/23  16:52:19  mringe
 * gauss() in Gauss() umenannt.
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
 * Revision 1.10  1993/08/10  14:29:19  mringe
 * Include string.h
 *
 * Revision 1.9  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.8  1993/07/29  14:47:00  mringe
 * GAP output (Option -G).
 *
 * Revision 1.7  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.6  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.5  1992/10/02  09:01:32  mringe
 * Limits f"ur L"ange der Namen auf 1023 erh"oht.
 *
 * Revision 1.4  1992/06/28  10:47:35  mringe
 * Ausgabe bei -e auf GAP-Format gebracht.
 *
 * Revision 1.3  1992/06/23  07:01:42  mringe
 * help() eingebaut. Verschiedene kleinere Korrekturen.
 * Neu: Option -e (GAP output)
 *
 * Revision 1.2  1992/06/16  14:53:23  mringe
 * F"ur nicht-Primk"orper erweitert. Benutze zitof().
 *
 * Revision 1.1  1992/05/26  07:56:01  mringe
 * Initial revision
 *
 */


#include <string.h>
#include <stdlib.h>
#include "meataxe.h"



#define MAXDEG 200		/* Max. degree */
#define LINEWIDTH 2048		/* Line width for input file */


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

long fl;
long nor, noc;
FILE *src = NULL;		/* Input file */
FILE *matfile;
char *matname;
char grpname[LINEWIDTH] = "";	/* Input selector (empty = ALL) */
FEL poly[MAXDEG+1];		/* Polynomial */
long deg;			/* Degree */
char name[LINEWIDTH];		/* Value (atlas notation) */
char thisgrp[LINEWIDTH];
int opt_G = 0;			/* -G option (GAP output) */
long *piv;
PTR mat, m2, m3;	/* Matrices */
PTR v2;			/* One row */


static char *helptext[] = {
"SYNTAX",
"    zev [-G] [<Matrix> [<Poly> [<Group>]]]",
"",
"OPTIONS",
"    -G      GAP output",
"",
"ARGUMENTS",
"    <Matrix>   A square matrix (default: G1)",
"    <Poly>     Data file with polynomials (default: standard input)",
"    <Group>    Selects a group of polynomials (default: all)",
"",
NULL};


static proginfo_t pinfo =
   { "zev", "Eigenvalues", "$Revision: 1.2 $", helptext };



/* ------------------------------------------------------------------
   init() - Initialze everything
   ------------------------------------------------------------------ */

void init(int argc, char **argv)

{
    int i;
    src = stdin;		/* Input file */

    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("G")) != OPT_END)
    {
	switch (i)
	{
	    case 'G': opt_G = 1; msg_level = -100; break;
	}
    }

    i = opt_ind;
    switch (argc-i)
    {
	case 3:
	    strcpy(grpname,argv[i+2]);
	case 2:
	    if (strcmp(argv[i+1],"-"))
	    {
		src = SFOpen(argv[i+1],FM_READ);
		if (src == NULL)
		{
		    perror(argv[i+1]);
		    errexit(ERR_FILEOPEN,argv[i+1]);
		}
	    }
	case 1:
	    matname = argv[i];
	case 0:
	    break;
	default:
	    errexit(ERR_NARGS,"zev");
    }
    if ((matfile = zreadhdr(matname,&fl,&nor,&noc)) == NULL)
	errexit(-1,matname);
    if (nor != noc) errexit(ERR_NOTSQUARE,argv[i]);
    zsetfield(fl); zsetlen(noc);
    mat = zalloc(nor);
    m2 = zalloc(nor);
    m3 = zalloc(nor);
    v2 = zalloc((long)1);
    zreadvec(matfile,mat,nor);
    fclose(matfile);
    piv = (long*) malloc((size_t)(nor+1) * sizeof(long));
    if (piv == NULL) errexit(ERR_NOMEM,"zev");
}





/* ------------------------------------------------------------------
   Gauss() - Calculate nullity and print the result
   ------------------------------------------------------------------ */

void Gauss(void)

{
    long dim;
    static int first = 1;

    dim = zmkechelon(m3,nor,piv);

    if (opt_G)	/* GAP output */
    {	long mult = (nor - dim) / deg;

    	if (mult > 0)
    	{
	    if (!first)
		printf(" + ");
	    else
	    {
		printf("MeatAxe.BrauerChar := ");
		first = 0;
	    }
	    printf("%ld*(%s)",mult,name);
	}
	if (mult*deg != (nor-dim))
	    fprintf(stderr,"Non-integer multiplicity for %s\n",name);
    }
    else	/* Table output */
    {
	printf("%10s %20s %4ld %4ld",thisgrp,name,deg,(nor-dim)/deg);
	if (((nor-dim)/deg)*deg != (nor-dim))
	    printf(" (non-integer)\n");
	else
	    printf("\n");
	fflush(stdout);
    }
}


/* ------------------------------------------------------------------
   readln() - Read one input line, strip comments
   ------------------------------------------------------------------ */

int readln(char *buf)

{
    char *c;
    int flag = 0;
    while (!flag)
    {	if (feof(src) || ferror(src)) return 0;
	*buf = 0;
	fgets(buf,LINEWIDTH,src);
	if (*buf == '#') continue;	/* comment */
	for (c = buf; *c != 0 & *c != '\n'; ++c)
		if (*c != ' ' && *c != '\t')
			flag = 1;		/* */
	*c = 0;					/* remove EOL */
    }
    return 1;
}


/* ------------------------------------------------------------------
   getnextpol() - Read the next polynomial
   ------------------------------------------------------------------ */

int getnextpol(void)

{
    char line[LINEWIDTH];
    FEL coeff[MAXDEG+1];
    char *c;
    int i;
    FEL f;

    while (1)
    {	if (!readln(line)) return 0;
	if (*line != ' ')	/* new group */
	{	strcpy(thisgrp,line);
		continue;
	}
	if (grpname[0] == 0 || !strcmp(thisgrp,grpname))
		break;
    }
    c = strtok(line," \t");
    strcpy(name,c);
    for (deg = 0; (c = strtok(NULL," \t")) != NULL; ++deg)
    {	if (!strcmp(c,"-1"))
		f = zneg(F_ONE);
	else
		f = atof(c);
	coeff[deg] = f;
    }
    --deg;
    for (i = 0; i <= (int) deg; ++i)
	poly[deg-i] = coeff[i];
    return 1;
}


/* ------------------------------------------------------------------
   eval() - Calculate f(A)
   ------------------------------------------------------------------ */

void eval(void)

{
    FEL f;
    long i, k;
    PTR x2, x3;

    memcpy(m2,mat,zsize(nor));	/* m2 = m1 */
    f = poly[0];
    x3 = m3;
    for (i = 1; i <= nor; ++i)	/* m3 = poly[0] * I */
    {	zmulrow(x3,F_ZERO);
    	zinsert(x3,i,f);
    	zadvance(&x3,(long)1);
    }
    for (k = 1; k <= deg; ++k)
    {	f = poly[k];
    	x2 = m2;		/* m3 = m3 + poly[k] * m2 */
    	x3 = m3;
    	for (i = 1; i <= nor; ++i)
    	{	zaddmulrow(x3,x2,f);
    		zadvance(&x2,(long)1);
    		zadvance(&x3,(long)1);
    	}
    	if (k == deg) break;
    	x2 = m2;		/* m2 = m2 * mat */
    	for (i = 1; i <= nor; ++i)
    	{	zmaprow(x2,mat,nor,v2);
    		zmoverow(x2,v2);
    		zadvance(&x2,(long)1);
    	}
    }
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char *argv[])

{
    init(argc,argv);
    while (getnextpol())
    {
	eval();
	Gauss();
    }
    if (opt_G) printf(";\n");
    return 0;
}


