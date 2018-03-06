/* ========================== C MeatAxe =============================
   zmu.c - Multiply matrices or permutations.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zmu.c,v 1.2 1997/09/11 15:44:03 gap Exp $
 *
 * $Log: zmu.c,v $
 * Revision 1.2  1997/09/11 15:44:03  gap
 * New version 2.2.3. AH
 *
 * Revision 2.9  1994/11/30  10:28:58  mringe
 * ANSI-C, Neue Random-Funktionen.
 *
 * Revision 2.8  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.7  1994/03/13  13:28:13  mringe
 * Benutze zreadlong()/zwritelong()
 *
 * Revision 2.6  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.5  1994/01/18  08:22:12  mringe
 * Bug behoben.
 *
 * Revision 2.4  1994/01/18  08:09:10  mringe
 * Erste Matrix zeilenweise einlesen.
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
 * Revision 1.19  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.18  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.17  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.16  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.15  1993/03/30  14:09:31  mringe
 * Benutze help-Funktion in der utils-Library.
 *
 * Revision 1.14  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.13  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.12  1993/01/06  20:59:57  mringe
 * getopt in() in zgetopt() umbenannt.
 *
 * Revision 1.11  1992/12/17  20:07:39  mringe
 * Einige Fehlermeldungen verbessert.
 *
 * Revision 1.10  1992/11/04  09:11:47  mringe
 * permrow() nach zzz.c verlagert.
 *
 * Revision 1.9  1992/11/02  09:15:44  mringe
 * Matrix*Perm: "Uberflu"ssiges zadvance() entfernt.
 *
 * Revision 1.8  1992/10/13  18:24:57  mringe
 * Pr"ufe, ob Filenamen verschieden sind.
 *
 * Revision 1.7  1992/10/13  18:17:11  mringe
 * Matrix*Perm: Lese die Matrix zeilenweise ein.
 *
 * Revision 1.6  1992/10/13  18:11:17  mringe
 * Perm*Matrix: Lese die Matrix zeilenweise.
 *
 * Revision 1.5  1992/07/29  10:47:51  mringe
 * Fixed a bug (default for #Rows, #Cols is 2)
 *
 * Revision 1.4  1992/07/29  10:04:07  mringe
 * Hilfstext verbessert.
 *
 * Revision 1.3  1992/07/23  18:33:07  mringe
 * Bug in multpp() behoben. Neu: multmp(), help().
 *
 * Revision 1.2  1992/07/16  14:58:00  mringe
 * Multiply any matrix with 1 by 1 matrix.
 *
 * Revision 1.1  1992/05/26  07:56:12  mringe
 * Initial revision
 *
 */

#include <stdlib.h>
#include "meataxe.h"




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static char *helptext[] = {
"SYNTAX",
"    zmu [-r <Row>[.<#Rows>]] [-c <Col>[.<#Cols>]] <A> <B> <Result>",
"",
"OPTIONS",
"    -r and -c can be used for piecewise matrix multiplication.",
"    E.g., `-r1 -c2' selects the upper right part. By default,",
"    <#Rows> = <#Cols> = 2, i.e., the result is divided into",
"    four pieces.",
"",
"FILES",
"    <A> and <B> are the objects to be multiplied. Their product",
"    (A*B) is written to <Result>. Compatible data types are:",
"",
"        M(a,b) * M(b,c)                   = M(a,c)",
"        M(1,1) * M(a,b) = M(a,b) * M(1,1) = M(a,b)",
"        P(a) * P(b)                       = P(max {a,b})",
"        M(a,b) * P(b)                     = M(a,b)",
"        P(a) * M(a,b)                     = M(a,b)",
"",
"    M(a,b) means `a by b matrix' and P(a) `Permutation of degree a'",
NULL};

static proginfo_t pinfo =
   { "zmu", "Multiply", "$Revision: 1.2 $", helptext };


char *aname, *bname, *cname;
FILE *afile, *bfile, *cfile;
static long fl1, fl2;
static long nor1, noc1, nor2, noc2;
int nrows = 1, thisrow = 1;	/* Arguments to -r option */
int ncols = 1, thiscol = 1;	/* Arguments to -c option */




/* ------------------------------------------------------------------
   err() - Print error message and exit.
   ------------------------------------------------------------------ */

static void err(int ch)

{   fprintf(stderr,"ZMU ERROR - ");
    switch (ch)
    {   case 'c':
	    fprintf(stderr,"INCOMPATIBLE OBJECTS\n");
	    break;
	case 'a':
	    fprintf(stderr,"BAD ARGUMENT TO -r OR -c\n");
	    break;
	case 'T':
	    fprintf(stderr,"MATRIX TOO SMALL\n");
	    break;
 	case '=':
	    fprintf(stderr,"IDENTICAL FILE NAMES NOT ALLOWED\n");
	    break;
 	case 'P':
	    fprintf(stderr,"MONOMIALS NOT SUPPORTED\n");
	    break;
 	case 'p':
	    fprintf(stderr,"PERMUTATIONS INCOMPATIBLE\n");
	    break;
 	case 'r':
	    fprintf(stderr,"ILLEGAL RANGE\n");
	    break;
    }
	exit(EXIT_ERR);
}


/* ------------------------------------------------------------------
   multpm() - Multiply permutation by matrix
   ------------------------------------------------------------------ */

static void multpm(void)

{
    PTR m2;
    long *p, i;

    if (nor1 != nor2) err('c');

    /* Allocate memory and read the permutation
       ---------------------------------------- */
    p = (long *) malloc(sizeof(long)*nor1);
    if (zreadlong(afile,p,nor1) != nor1)
	errexit(-1,aname);
    --p;
    zsetfield(fl2);
    zsetlen(noc2);
    m2 = zalloc((long) 1);

    /* Write out the rows in permuted order
       ------------------------------------ */
    if ((cfile = zwritehdr(cname,fl2,nor2,noc2)) == NULL)
	errexit(-1,cname);
    for (i = 1; i <= nor1; ++i)
    {	
	zseek(bfile,p[i]);
	zreadvec(bfile,m2,(long) 1);
	zwritevec(cfile,m2,(long)1);
    }
}


/* ------------------------------------------------------------------
   multmp() - Multiply matrix by permutation
   ------------------------------------------------------------------ */

static void multmp(void)

{
    PTR m1,m3;
    long *m2;

    if (noc1 != nor2) err('c');

    zsetfield(fl1);
    zsetlen(noc1);
    m1 = zalloc((long) 1);
    m3 = zalloc((long) 1);
    m2 = NALLOC(long,nor2);
    if (zreadlong(bfile,m2,nor2) != nor2)
	errexit(-1,bname);

    if ((cfile = zwritehdr(cname,fl1,nor1,noc1)) == NULL)
	errexit(-1,cname);
    zsetlen(noc1);
    while (nor1--)
    {
        zreadvec(afile,m1,(long) 1);
	zpermrow(m1,m2,m3);
	zwritevec(cfile,m3,(long) 1);
    }
}



/* ------------------------------------------------------------------
   multsm() - Multiply scalar with matrix
   ------------------------------------------------------------------ */

static void multsm(FILE *s,FILE *m,long nor3,long noc3)

{
    PTR ms, mm;
    FEL f;
    long i;

    zsetfield(fl1);
    zsetlen(1);
    ms = zalloc(1);
    zreadvec(s,ms,1);
    f = zextract(ms,1);
    zsetlen(noc3);
    mm = zalloc((long)1);
	
    if ((cfile = zwritehdr(cname,fl1,nor3,noc3)) == NULL)
	errexit(-1,cname);
    for (i = 1; i <= nor3; ++i)
    {	
	zreadvec(m,mm,(long)1);
	zmulrow(mm,f);
	zwritevec(cfile,mm,(long)1);
    }
}



/* ------------------------------------------------------------------
   multmm() - Multiply two matrices
   ------------------------------------------------------------------ */

static void multmm(void)

{
    PTR m1,m2,tmp,out,x;
    long row1, row2;	/* Range of rows */
    long col1, col2;	/* Range of columns */
    long i, k;

    if (fl1 != fl2) err('c');
    if (noc1 == 1 && nor1 == 1)
    {	multsm(afile,bfile,nor2,noc2);
	return;
    }
    else if (noc2 == 1 && nor2 == 1)
    {	multsm(bfile,afile,nor1,noc1);
	return;
    }
    if (noc1 != nor2) err('c');
    if ((long) nrows > nor1 || (long) ncols > noc2) err('T');
    row1 = (nor1 / nrows) * (thisrow - 1) + 1;
    row2 = (nor1 / nrows) * thisrow;
    col1 = (noc2 / ncols) * (thiscol - 1) + 1;
    col2 = (noc2 / ncols) * thiscol;

    /* First matrix
       ------------ */
    zsetfield(fl1);
    zsetlen(noc1);
    m1 = zalloc((long)1);
    zseek(afile,row1);

    /* Read second matrix
       ------------------ */
    zsetlen(noc2);
    tmp = zalloc((long)1);
    zsetlen(col2-col1+1);
    m2 = zalloc(nor2);
    out = zalloc((long)1);
    if (col2 - col1 + 1 < noc2)
    {
	x = m2;
	for (i = 1; i <= nor2; ++i)
	{
	    zsetlen(noc2);
	    zreadvec(bfile,tmp,(long)1);
	    for (k = col1; k <= col2; ++k)
		zinsert(x,k - col1 + 1,zextract(tmp,k));
	    zsetlen(col2-col1+1);
	    zadvance(&x,(long)1);
	}
    }
    else
    {
	zsetlen(noc2);
	zreadvec(bfile,m2,nor2);
    }

    /* Multiply and write result
       ------------------------- */
    if ((cfile = zwritehdr(cname,fl1,row2-row1+1,col2-col1+1)) == NULL)
	errexit(-1,cname);
    for (i = 1; i <= row2 - row1 + 1; ++i)
    {
        zsetlen(noc1);
        zreadvec(afile,m1,(long)1);	/* Read one row */
	zsetlen(col2-col1+1);
	zmaprow(m1,m2,nor2,out);
	zwritevec(cfile,out,(long)1);
    }
}


/* ------------------------------------------------------------------
   multpp() - Multiply two permutations
   ------------------------------------------------------------------ */

static void multpp(void)

{
    long i;
    long *p1, *p2;
    long maxdeg;

    if (fl1 != fl2) err('p');
    if (fl1 != -1) err('P');
    maxdeg = (nor1 > nor2) ? nor1 : nor2;

    /* Allocate memory and read the permutations
       ----------------------------------------- */
    p1 = (long *) malloc(sizeof(long)*maxdeg);
    p2 = (long *) malloc(sizeof(long)*maxdeg);
    if (p1 == NULL || p2 == NULL) errexit(ERR_NOMEM,"zmu");
    if (zreadlong(afile,p1,nor1) != nor1)
	errexit(-1,aname);
    if (zreadlong(bfile,p2,nor2) != nor2)
	errexit(-1,bname);
    --p1;
    --p2;

    /* If the permutations have different
       degrees, extend the smaller one.
       ---------------------------------- */
    for (i = nor1 + 1; i <= maxdeg; ++i)
	p1[i] = i;
    for (i = nor2 + 1; i <= maxdeg; ++i)
	p2[i] = i;

    /* Calculate the product
       --------------------- */
    for (i = 1; i <= maxdeg; ++i)
	p1[i] = p2[p1[i]];

	/* Write out the result
	   -------------------- */
    if ((cfile = zwritehdr(cname,fl1,maxdeg,(long)1)) == NULL)
	errexit(-1,cname);
    if (zwritelong(cfile,p1+1,maxdeg) != maxdeg)
	errexit(-1,cname);
}



int main(int argc, char *argv[])

{   int i;


    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("r:c:")) != OPT_END)
    {
	switch (i)
	{
	    case 'c':
		thiscol = (int) getint();
		if (*opt_text_ptr == '.')
		{	++opt_text_ptr;
			ncols = (int) getint();
		}
		else
			ncols = 2;
		if (*opt_text_ptr != 0) errexit(ERR_OPTION,"-c");
		if (thiscol < 1 || ncols < thiscol) err('a');
		break;
	    case 'r':
		thisrow = (int) getint();
		if (*opt_text_ptr == '.')
		{	++opt_text_ptr;
			nrows = (int) getint();
		}
		else
			nrows = 2;
		if (*opt_text_ptr != 0) errexit(ERR_OPTION,"-r");
		if (thisrow < 1 || nrows < thisrow) err('a');
		break;
	}
    }

    /* Set file names
       -------------- */
    switch (argc-opt_ind)
    {
	case 3:
	    aname = argv[opt_ind];
	    bname = argv[opt_ind+1];
	    cname = argv[opt_ind+2];
	    break;
	default:
	    errexit(ERR_BADUSAGE,"zmu");
    }
   if (!strcmp(aname,cname) || !strcmp(bname,cname))
	err('=');	/* Don't overwrite input files */

    /* Read file headers and call the 
       appropriate function for multiplying.
       ------------------------------------- */
    afile = zreadhdr(aname,&fl1,&nor1,&noc1);
    if (afile == NULL) errexit(-1,aname);
    bfile = zreadhdr(bname,&fl2,&nor2,&noc2);
    if (bfile == NULL) errexit(-1,bname);
    if (fl1 < 0 && fl2  < 0) multpp(); else
    if (fl1 > 1 && fl2 > 1) multmm(); else
    if (fl1 > 1 && fl2 == -1) multmp(); else
    if (fl1 == -1 && fl2  > 1) multpm(); else
	err('c');
    fclose(afile);
    fclose(bfile);
    fclose(cfile);
    return EXIT_OK;
}



