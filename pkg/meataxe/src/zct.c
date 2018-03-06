/* ========================== C MeatAxe =============================
   zct.c - Cut file.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zct.c,v 1.2 1997/09/11 15:43:43 gap Exp $
 *
 * $Log: zct.c,v $
 * Revision 1.2  1997/09/11 15:43:43  gap
 * New version 2.2.3. AH
 *
 * Revision 2.10  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.9  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.8  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.7  1994/05/20  12:41:25  mringe
 * Prototypen.
 *
 * Revision 2.6  1994/03/13  13:27:01  mringe
 * Maschinenunabhaengiges Format fuer Permutationen.
 *
 * Revision 2.5  1994/02/15  14:09:19  mringe
 * zsetlen() erweitert fuer Permutationen.
 *
 * Revision 2.4  1994/02/15  13:12:35  mringe
 * Benutze zseek() statt SFSeek().
 *
 * Revision 2.3  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
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
 * Revision 1.9  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.8  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.7  1993/08/10  14:44:19  mringe
 * Include string.h
 *
 * Revision 1.6  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.5  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.4  1993/06/30  16:12:13  mringe
 * Bug bei der Erkennung der Zeilen-/Spaltenliste behoben.
 *
 * Revision 1.3  1993/06/30  11:32:49  mringe
 * help() verbessert, benutze ':' statt ';'
 *
 * Revision 1.2  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.1  1992/05/26  07:55:57  mringe
 * Initial revision
 *
 */


#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "meataxe.h"


#define MAXPIECES 10


/* ------------------------------------------------------------------
   Variables
   ------------------------------------------------------------------ */


static char *helptext[] = {
"SYNTAX",
"    zct <Rows>[:<Columns>] <Input> <Output>",
"",
"    <Rows> and <Columns> are lists of integers or ranges",
"    separated by commas. Example: 1-5,6,8-9:4,10-11",
"",
"FILES",
"    <Input>    The input file (a matrix or permutations)",
"    <Output>   The output file",
NULL};

static proginfo_t pinfo =
   { "zcd", "Cut Matrices Or Permutations",
     "$Revision: 1.2 $", helptext };


	
int nrows = 0;
long rowlist[MAXPIECES][2];
int ncols = 0;
long collist[MAXPIECES][2];
char *savedpos;
char *ifilename;
char *ofilename;
long fl, nor, noc;
long onor, onoc;

FILE *inp, *out;



/* ------------------------------------------------------------------
   err() - Print error message and exit
   ------------------------------------------------------------------ */

static void err(c)
int c;

{   fprintf(stderr,"ZCT ERROR - ");
    switch (c)
    {	
	case 'a':
	    fprintf(stderr,"BAD RANGE\n");
	    break;
	case 'x':
	    fprintf(stderr,"TOO MANY PIECES\n");
	    break;
	case 'P':
	    fprintf(stderr,"CANNOT CUT COLUMNS FROM A PERMUTATION\n");
	    break;
	case 'o':
	    fprintf(stderr,"CANNOT CUT THE REQUESTED PIECE\n");
	    break;
    }
    exit(EXIT_ERR);
}


/* ------------------------------------------------------------------
   list() - Scan a list of rows/columns
   ------------------------------------------------------------------ */

static int list(x, c)
long x[MAXPIECES][2];
char *c;

{   
    int n = 0;

    while (isdigit(*c))
    {
	if (n >= MAXPIECES) err('x');
	x[n][0] = atol(c);
	while (isdigit(*c)) ++c;
	if (*c == '-')
	{	++c;
		if (!isdigit(*c))
    	    	    errexit(ERR_BADUSAGE,"zct");
		x[n][1] = atol(c);
		while (isdigit(*c)) ++c;
	}
	else
		x[n][1]= x[n][0];
	if (*c == ',') ++c;
	if (x[n][1] < x[n][0]) err('a');
	++n;
    }
    savedpos = c;
    return n;
}




/* ------------------------------------------------------------------
   parseargs()
   ------------------------------------------------------------------ */

static void parseargs(argc,argv)
int argc;
char *argv[];

{
    int i;
    char *c;

    /* Parse command line
       ------------------ */
    initargs(argc, argv, &pinfo);
    while (zgetopt("") != OPT_END);

    /* Process the <Range> argument
       ---------------------------- */
    i = opt_ind;
    if (i >= argc) errexit(ERR_BADUSAGE,"zct");
    c = argv[i];
    nrows = list(rowlist,c);
    c = savedpos;
    if (*c != 0)
    {
        if (*c != ':' && *c != ';')
    	    errexit(ERR_BADUSAGE,"zct");
        ++c;
        ncols = list(collist,c);
    }
    ++i;

    /* Process file names
       ------------------ */
    if (i != argc - 2) errexit(ERR_NARGS,"zct");
    ifilename = argv[i];
    ofilename = argv[i+1];
}



/* ------------------------------------------------------------------
   init()
   ------------------------------------------------------------------ */

static void init()

{
    int i;
	
    mtxinit();
    inp = zreadhdr(ifilename,&fl,&nor,&noc);
    if (inp == NULL) errexit(-1,ifilename);
    if (fl == -1)	/* Is it a permutation? */
    {	if (ncols > 0) err('P');
    }
    else
    {	if (ncols == 0)
	{	ncols = 1;
		collist[0][0] = 1;
		collist[0][1] = noc;
	}
    }

    if (nrows == 0)
	onor = nor;
    else
    {	onor = 0;
	for (i = 0; i < nrows; ++i)
		onor += rowlist[i][1]-rowlist[i][0]+1;
    }
    if (ncols == 0)
	onoc = noc;
    else
    {	onoc = 0;
	for (i = 0; i < ncols; ++i)
		onoc += collist[i][1]-collist[i][0]+1;
    }
}


/* ------------------------------------------------------------------
   cutperm()
   ------------------------------------------------------------------ */

static void cutperm()

{
    int i;
    long *x, *y;

    if (rowlist[nrows-1][1] > noc) err('o');

    /* Read permutations
       ----------------- */
    y = x = (long *) malloc(nor * onor * sizeof(long));
    if (x == NULL) errexit(ERR_NOMEM,"zct");
    for (i = 0; i < nrows; ++i)
    {
	size_t n = nor * (rowlist[i][1]-rowlist[i][0]+1);
	SFSeek(inp,12+(rowlist[i][0]-1)*nor*4);
	if (zreadlong(inp,y,n) != n)
	    errexit(-1,ifilename);
	y += n;
    }

    /* Write output
       ------------ */
    if ((out = zwritehdr(ofilename,(long)-1,nor,onor)) == NULL)
	errexit(-1,ofilename);
    if (zwritelong(out,x,onor*nor) != onor*nor)
	errexit(-1,ofilename);
}


/* ------------------------------------------------------------------
   cutmatrix()
   ------------------------------------------------------------------ */

static void cutmatrix()

{   int i, ii;
    long k, kk, pos;
    PTR x,y,row;

    if (rowlist[nrows-1][1] > nor) err('o');
    if (collist[ncols-1][1] > noc) err('o');
    zsetfield(fl);
    zsetlen(noc);
    row = zalloc((long) 1);
    zsetlen(onoc);
    x = zalloc(onor);
    y = x;
    for (i = 0; i < nrows; ++i)
    {	zsetlen(noc);
	zseek(inp,rowlist[i][0]);
	for (k = 0; k <= rowlist[i][1]-rowlist[i][0]; ++k)
	{   zsetlen(noc);
	    zreadvec(inp,row,1);
	    pos = 1;
	    zsetlen(onoc);
	    zmulrow(y,F_ZERO);
	    for (ii = 0; ii < ncols; ++ii)
		for (kk = collist[ii][0]; kk <= collist[ii][1]; ++kk)
		{   FEL f = zextract(row,kk);
		    zinsert(y,pos,f);
		    ++pos;
		}
	    zadvance(&y,(long)1);
	}
    }

    /* Write output
       ------------ */
    zsetlen(onoc);
    if ((out = zwritehdr(ofilename,fl,onor,onoc)) == NULL)
	errexit(-1,ofilename);
    zwritevec(out,x,onor);
}



int main(argc, argv)
int argc;
char *argv[];

{


    parseargs(argc,argv);

    init();
    if (fl == -1)
	cutperm();
    else
	cutmatrix();
    fclose(inp);
    fclose(out);
    return EXIT_OK;
}



