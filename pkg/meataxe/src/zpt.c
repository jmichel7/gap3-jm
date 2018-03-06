/* ========================== C MeatAxe =============================
   zpt.c - Paste matrices or permutations.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zpt.c,v 1.2 1997/09/11 15:44:19 gap Exp $
 *
 * $Log: zpt.c,v $
 * Revision 1.2  1997/09/11 15:44:19  gap
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
 * Revision 2.3  1993/10/28  19:15:38  mringe
 * Erlaube DImension 0 bei Teilstuecken.
 *
 * Revision 2.2  1993/10/21  21:57:35  mringe
 * Permutationen.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
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
 * Revision 1.10  1993/08/10  14:52:15  mringe
 * Include string.h
 *
 * Revision 1.9  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.8  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.7  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.6  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.5  1993/01/06  20:59:57  mringe
 * getopt in() in zgetopt() umbenannt.
 *
 * Revision 1.4  1992/11/04  07:24:45  mringe
 * Behandlung von Permutationen berichtigt (Thanks to Klaus)
 *
 * Revision 1.3  1992/10/02  11:49:02  mringe
 * Erlaube "-" als Filenamen (= Nullmatrix).
 *
 * Revision 1.2  1992/07/29  10:41:38  mringe
 * Fixed some bugs; use zgetopt().
 *
 * Revision 1.1  1992/05/26  17:39:40  mringe
 * Initial revision
 *
 */


#include <string.h>
#include <stdlib.h>
#include "meataxe.h"




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

int nrows = 0;
int ncols = 0;
int opt_r = 0, opt_R = 0;
int opt_c = 0, opt_C = 0;
char *ofilename;
char **ifilename;
FILE *ifile, *ofile;
long fl, nor, noc, maxnor;
long *width, *height;

static char *helptext[] = {
"SYNTAX",
"    zpt [-{c|C} <NCols>] [-{r|R} <NRows>] <Out> [<Inp> ...]",
"",
"OPTIONS",
"    -c|-C   Set number of columns.",
"    -r|-R   Set number of rows.",
"",
"    -c/-r expect input files specified separately, -C/-R expect",
"    a base name and append one or two numbers to build the input",
"    file name. `-' as input file name represents the zero matrix.",
"",
"FILES",
"    <Out>   o   The result goes here.",
"    <Inp>   i   The pieces.",
NULL};

static proginfo_t pinfo =
   { "zpt", "Paste Matrices Or Permutations",
     "$Revision: 1.2 $", helptext };



/* ------------------------------------------------------------------
   err() - Print error message and exit
   ------------------------------------------------------------------ */

static void err(int c)

{   fprintf(stderr,"ZPT ERROR - ");
    switch (c)
    {   
	case 'p':
	    fprintf(stderr,"CANNOT PASTE PERMUTATIONS THIS WAY\n");
	    break;
	case '0':
	    fprintf(stderr,"UNDETERMINED SIZE(S)\n");
	    break;
    }
    exit(EXIT_ERR);
}


/* ------------------------------------------------------------------
   myitoa()
   ------------------------------------------------------------------ */

static void myitoa(i, c)
int i;
char *c;

{	char buf[10];
	int l = 0;

	if (i == 0)
	{	l = 1;
		buf[0] = '0';
	}
	else
		while (i != 0)
		{	buf[l++] = (char)((i % 10) + '0');
			i /= 10;
		}
	while (l > 0)
		*c++ = buf[--l];
	*c = 0;
}



/* ------------------------------------------------------------------
   init()
   ------------------------------------------------------------------ */

static void init(argc,argv)
int argc;
char *argv[];

{   int i, names_needed;


    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("c:C:r:R:")) != OPT_END)
    {
	switch (i)
	{
	    case 'r':
	    case 'R':
		if (opt_r || opt_R) errexit(ERR_OPTION,"-r");
		if (i == 'r') opt_r = 1; else opt_R = 1;
		nrows = getint();
		if (nrows == GETINT_ERR || *opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-r");
		break;
	    case 'c':
	    case 'C':
		if (opt_c || opt_C) errexit(ERR_OPTION,"-c");
		if (i == 'c') opt_c = 1; else opt_C = 1;
		ncols = getint();
		if (ncols == GETINT_ERR || *opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-c");
		break;
	}
    }
    i = opt_ind;
    if (i >= argc-1) errexit(ERR_NARGS,"zpt");
    ofilename = argv[i++];
    ifilename = argv + i;

    if (!opt_r && !opt_R && !opt_c && !opt_C)
    {	opt_c = opt_r = 1;
	ncols = 1;
	nrows = argc - i;
	names_needed = nrows;
    }
    else if (!opt_r && !opt_R)
    {	nrows = 1;
	opt_r = 1;
	names_needed = opt_C ? 1 : ncols;
    }
    else if (!opt_C && !opt_c)
    {	ncols = 1;
	opt_c = 1;
	names_needed = opt_R ? 1 : nrows;
    }
    else
	names_needed = (opt_C ? 1 : ncols) * (opt_R ? 1 : nrows);

    if (argc - i != names_needed) errexit(ERR_NARGS,"zpt");
    if (nrows < 1 || ncols < 1) errexit(ERR_BADUSAGE,"zpt");
    width = (long *) malloc((size_t)(ncols+1)*sizeof(long));
    height = (long *) malloc((size_t)(nrows+1)*sizeof(long));
}



/* ------------------------------------------------------------------
   mkname()
   ------------------------------------------------------------------ */

static char *mkname(r, c)
int r, c;

{
    int i;
    char *x;
    static char buf[200];

    if (opt_R && opt_C)
	i = 0;
    else if (opt_R)
	i = c - 1;
    else if (opt_C)
	i = r - 1;
    else
	i = (r-1) * ncols + (c-1);
    strcpy(buf,ifilename[i]);
    for (x = buf; *x != 0; ++x);
    if (opt_R)
    {	myitoa(r,x);
	while (*x != 0) ++x;
    }
    if (opt_C)
    	myitoa(c,x);
    return buf;
}


/* ------------------------------------------------------------------
   checksizes()
   ------------------------------------------------------------------ */

static void checksizes()

{
    int i, k;
    long fl2, nor2, noc2;
    char *c;

    fl = 0;
    for (i = 1; i <= nrows; ++i) height[i] = -1;
    for (k = 1; k <= ncols; ++k) width[k] = -1;

    for (i = 1; i <= nrows; ++i)
    {
	for (k = 1; k <= ncols; ++k)
	{
	    if (!strcmp(c=mkname(i,k),"-"))
		fl2 = nor2 = noc2 = 0;
	    else
	    {
	        char buf[200];
		char *inam;

		inam = mkname(i,k);
		if ((ifile = zreadhdr(inam,&fl2,&nor2,&noc2)) == NULL)
		    errexit(-1,inam);
		fclose(ifile);
	    
		if (fl == 0) fl = fl2;
    		if (fl == -1 && ncols != 1) err('p');
		else if (fl != fl2)
		{
		    sprintf(buf,"%s and %s",mkname(1,1),inam);
		    errexit(ERR_INCOMPAT,buf);
		}

		if (fl == -1)
		{
		    int tmp;
		    tmp = nor2;
		    nor2 = noc2;
		    noc2 = tmp;
		}
	    
		if (height[i] == -1) height[i] = nor2;
		else if (height[i] != nor2)
		{
		    sprintf(buf,"%s and %s",mkname(i,k),mkname(i,1));
		    errexit(ERR_INCOMPAT,buf);
		}

		if (width[k] == -1) width[k] = noc2;
		else if (width[k] != noc2)
		{
		    sprintf(buf,"%s and %s",mkname(i,k),mkname(1,k));
		    errexit(ERR_INCOMPAT,buf);
		}

	    }
	}
    }

    /* Calculate nor, noc and maxnor
       ----------------------------- */
    noc = nor = maxnor = 0;
    for (i = 1; i <= nrows; ++i)
    {
	if (height[i] == -1) err('0');
	if (height[i] > maxnor) maxnor = height[i];
	nor += height[i];
    }
    for (k = 1; k <= ncols; ++k)
    {
	if (width[k] == -1) err('0');
	noc += width[k];
    }
}


/* ------------------------------------------------------------------
   pasteperm() - Concatenate permutations
   ------------------------------------------------------------------ */

static void pasteperm()

{
    long *p;
    int i;
    long fl2, npts, nperms;

    p = (long *) malloc(sizeof(long)*noc);
    if (p == NULL) errexit(ERR_NOMEM,"zpt");
    if ((ofile = zwritehdr(ofilename,fl,noc,nor)) == NULL)
	errexit(-1,ofilename);
    for (i = 1; i <= nrows; ++i)
    {
	char *inam = mkname(i,1);
	if ((ifile = zreadhdr(inam,&fl2,&npts,&nperms)) == NULL)
	    errexit(-1,inam);
	while (nperms-- > 0)
	{
	    if (zreadlong(ifile,p,noc) != noc)
		errexit(-1,inam);
	    if (zwritelong(ofile,p,noc) != noc)
		errexit(-1,ofilename);
	}
	fclose(ifile);
    }
}



/* ------------------------------------------------------------------
   pastemat() - Paste matrices
   ------------------------------------------------------------------ */

static void pastemat()

{   
    PTR m, piece, x, y;
    int i, k;
    long l, j;
    long fl2, nor2, noc2, pos;

    zsetfield(fl);
    zsetlen(noc);
    m = zalloc(maxnor);
    if ((ofile = zwritehdr(ofilename,fl,nor,noc)) == NULL)
	errexit(-1,ofilename);
    for (i = 1; i <= nrows; ++i)
    {	zsetlen(noc);
	x = m;
	for (l = maxnor; l > 0; --l)
	{	zmulrow(x,F_ZERO);
		zadvance(&x,(long)1);
	}
	pos = 0;
	for (k = 1; k <= ncols; ++k)
	{
	    char *c = mkname(i,k);

	    if (strcmp(c,"-"))
	    {
		if ((ifile = zreadhdr(c,&fl2,&nor2,&noc2)) == NULL)
		    errexit(-1,c);
		zsetlen(noc2);
		piece = zalloc(nor2);
		zreadvec(ifile,piece,nor2);
		fclose(ifile);
		x = m;
		y = piece;
		for (l = 1; l <= nor2; ++l)
		{	for (j = 1; j <= noc2; ++j)
				zinsert(x,j+pos,zextract(y,j));
			zsetlen(noc2);
			zadvance(&y,(long)1);
			zsetlen(noc);
			zadvance(&x,(long)1);
		}
		free(piece);
	    }
	    pos += width[k];
	}
	zsetlen(noc);
	zwritevec(ofile,m,height[i]);
    }
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{
    init(argc,argv);
    checksizes();
    if (fl == -1)
	pasteperm();
    else
	pastemat();
    fclose(ofile);
    return EXIT_OK;
}

