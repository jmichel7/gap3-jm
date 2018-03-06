/* ========================== C MeatAxe =============================
   zsi.c -  Sum and intersection of two subspaces.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zsi.c,v 1.2 1997/09/11 15:44:25 gap Exp $
 *
 * $Log: zsi.c,v $
 * Revision 1.2  1997/09/11 15:44:25  gap
 * New version 2.2.3. AH
 *
 * Revision 2.9  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.8  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.7  1994/07/07  06:58:50  mringe
 * Include files.
 *
 * Revision 2.6  1994/02/25  11:21:27  mringe
 * Bug in -s entfernt.
 *
 * Revision 2.5  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.4  1994/02/04  16:19:45  mringe
 * Option -s berechnet nur die Summe.
 *
 * Revision 2.3  1993/12/08  11:48:50  mringe
 * Compiler warnings.
 *
 * Revision 2.2  1993/10/28  19:10:15  mringe
 * BUG behoben (2. Argument wird jetzt eingelesen).
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
 * Revision 1.7  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.6  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.5  1993/07/19  15:57:53  mringe
 * Optionen -Q, -V, -G.
 *
 * Revision 1.4  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.3  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.2  1992/07/15  09:25:55  mringe
 * Some minor changes.
 *
 * Revision 1.1  1992/05/26  18:14:12  mringe
 * Initial revision
 *
 */


#include <stdlib.h>
#include <string.h>
#include "meataxe.h"



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

long fl, nor1, nor2, noc;
long *piv;
FILE *afile, *bfile;
long sumdim, secdim;
PTR sum, sec, intersection;
char *aname, *bname, *sumname, *intname;

static char *helptext[] = {
"SYNTAX",
"    zsi [-GQV] <Space1> <Space2> <Sum> <Int>",
"    zsi -s [-GQV] <Space1> <Space2> <Sum>",
/*"    zsi -i [-GQV] <Space1> <Space2> <Int>",*/
"",
"OPTIONS",
"    -s     Calculate sum only",
/*"    -i     Calculate intersection only",*/
"    -G     GAP output",
"    -Q     Quiet, no messages",
"    -V     Verbose, more messages",
"",
"FILES",
"    <Spc1>    i  The first space, a matrix.",
"    <Spc2>    i  The second space, a matrix.",
"    <Sum>     o  The sum of the two spaces.",
"    <Int>     o  The intersection of the two spaces.",
NULL};

static proginfo_t pinfo =
   { "zsi", "Sum And Intersection", "$Revision: 1.2 $", helptext };


int opt_G = 0;
int opt_s = 0;

/* ------------------------------------------------------------------
   openinp() - Open input files
   ------------------------------------------------------------------ */

static void openinp()

{
    long fl1, fl2, noc1, noc2;

    if ((afile = zreadhdr(aname,&fl1,&nor1,&noc1)) == NULL)
	errexit(-1,aname);
    MESSAGE(1,("%s: %ld vectors, row size = %ld\n",aname,nor1,noc1));
    if (fl1 < 2) errexit(ERR_NOTMATRIX,aname);
    if ((bfile = zreadhdr(bname,&fl2,&nor2,&noc2)) == NULL)
	errexit(-1,bname);
    MESSAGE(1,("%s: %ld vectors, row size = %ld\n",bname,nor2,noc2));
    if (fl2 < 2) errexit(ERR_NOTMATRIX,bname);
    if (fl1 != fl2 || noc1 != noc2)
    {
        char buf[200];
	sprintf(buf,"%s and %s",aname,bname);
    	errexit(ERR_INCOMPAT,buf);
    }

    fl = fl1;
    noc = noc1;
    zsetfield(fl); zsetlen(noc);
}


/* ------------------------------------------------------------------
   sumonly() - Calculate the sum
   ------------------------------------------------------------------ */

static void sumonly()

{
    FILE *cfile;
    PTR a, b;
    long maxrows;

    /* Read first space
       ---------------- */
    openinp();
    maxrows = (noc < nor1 + nor2) ? noc : nor1 + nor2;
    a = zalloc(maxrows);
    zreadvec(afile,a,nor1);
    fclose(afile);

    /* Make echelon form of first space
       -------------------------------- */
    piv = (long *) malloc(sizeof(long) * (size_t)(maxrows+1));
    if (piv == NULL) errexit(ERR_NOMEM,"zsi");
    nor1 = zmkechelon(a,nor1,piv);

    /* Clean vectors from 2nd space and extend first space
       --------------------------------------------------- */
    b = a;
    zadvance(&b,nor1);
    while (nor1 < maxrows && nor2 > 0)	/* Any vectors left? */
    {
	long nread;
	PTR x;

	/* Read as many vectors as possible
	   -------------------------------- */
	nread = maxrows - nor1;
	if (nread > nor2) nread = nor2;
        zreadvec(bfile,b,nread);
	nor2 -= nread;

	/* Clean the new vectors
	   --------------------- */
	for (x = b; nread > 0; --nread, zadvance(&x,(long)1))
	{
	    long p;
	    FEL f;
	    zcleanrow(x,a,nor1,piv);
	    if ((p = zfindpiv(x,&f)) > 0)  /* A new vector! */
	    {
		++nor1;
		piv[nor1] = p;
		if (x != b) zmoverow(b,x);
		zadvance(&b,1);
	    }
	}
    }
    fclose(bfile);

    MESSAGE(0,("SUM %ld\n",nor1));
    if ((cfile = zwritehdr(sumname,fl,nor1,noc)) == NULL)
	errexit(-1,sumname);
    zwritevec(cfile,a,nor1);
    fclose(cfile);
}


/* ------------------------------------------------------------------
   readfiles() - Read the input files, allocate tables
   ------------------------------------------------------------------ */

static void readfiles()

{
    PTR x;

    openinp();
    sum = zalloc(nor1+nor2);
    sec = zalloc(nor1+nor2);
    zreadvec(afile,sum,nor1);
    fclose(afile);
    x = sum;
    zadvance(&x,nor1);
    zreadvec(bfile,x,nor2);
    fclose(bfile);
    memcpy(sec,sum,zsize(nor1));
    piv = (long *) malloc(sizeof(long) * (size_t)(nor1+nor2+1));
    if (piv == NULL) errexit(ERR_NOMEM,"zsi");
}


/* ------------------------------------------------------------------
   suminter() - Calculate sum and intersection (Zassenhaus algorithm)
   ------------------------------------------------------------------ */

static void suminter()

{
    long i, k, p = 0;
    PTR xi, xk, yi, yk;
    FEL f;

    /* Step 1
       ------ */
    xi = sum;
    yi = sec;
    for (i = 1; i <= nor1+nor2; )
    {
    	xk = xi;
    	yk = yi;
    	for (k = i; k <= nor1+nor2; ++k)
    	{	p = zfindpiv(xk,&f);
		if (p != 0) break;
		zadvance(&xk,(long)1);
		zadvance(&yk,(long)1);
	}
	if (p == 0) break;
	if (k != i)
	{	zswaprow(xk,xi);
		zswaprow(yk,yi);
	}
	xk = sum;
	yk = sec;
	for (k = 1; k < i; ++k)
	{	f = zneg(zextract(xi,piv[k]));
		zaddmulrow(xi,xk,f);
		zaddmulrow(yi,yk,f);
		zadvance(&xk,(long)1);
		zadvance(&yk,(long)1);
	}
	p = zfindpiv(xi,&f);
	if (p != 0)
	{	piv[i++] = p;
		zmulrow(xi,zinv(f));
		zmulrow(yi,zinv(f));
		zadvance(&xi,(long)1);
		zadvance(&yi,(long)1);
	}
    }
    sumdim = i-1;
    intersection = yi;
    MESSAGE(1,("Step 1 completed, sumdim = %ld\n",sumdim));

    /* Step 2
       ------ */
    while (1)
    {	if (i > nor1+nor2) break;
		yk = yi;
		for (k = i; k <= nor1+nor2; ++k)
		{	p = zfindpiv(yk,&f);
			if (p != 0) break;
			zadvance(&yk,(long)1);
		}
		if (p == 0) break;
		if (k != i)
			zswaprow(yk,yi);
		yk = intersection;
		for (k = sumdim+1; k < i; ++k)
		{	f = zneg(zextract(yi,piv[k]));
			zaddmulrow(yi,yk,f);
			zadvance(&yk,(long)1);
		}
		p = zfindpiv(yi,&f);
		if (p != 0)
		{	piv[i++] = p;
			zmulrow(yi,zinv(f));
			zadvance(&yi,(long)1);
		}
    }
    secdim = i - sumdim - 1;

    /* Print result
       ------------ */
    if (opt_G)
    {
        printf("MeatAxe.SumDimension := %ld;\n",sumdim);
        printf("MeatAxe.IntDimension := %ld;\n",secdim);
    }
    else
        MESSAGE(0,("SUM %ld INTERSECTION %ld\n",sumdim,secdim));
}



/* ------------------------------------------------------------------
   writefiles() - Write out the result
   ------------------------------------------------------------------ */

static void writefiles()

{	
    FILE *f;
    
    if ((f = zwritehdr(sumname,fl,sumdim,noc)) == NULL)
	errexit(-1,sumname);
    zwritevec(f,sum,sumdim);
    fclose(f);
    if ((f = zwritehdr(intname,fl,secdim,noc)) == NULL)
	errexit(-1,intname);
    zwritevec(f,intersection,secdim);
    fclose(f);
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{
    int i;

    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("sG")) != OPT_END)
    {
	switch (i)
	{
	    case 'G': opt_G = 1; msg_level = -100; break;
	    case 's': opt_s = 1; break;
	}
    }
    if (opt_ind != argc - (opt_s ? 3 : 4)) errexit(ERR_NARGS,"zsi");
    aname = argv[opt_ind];
    bname = argv[opt_ind+1];
    sumname = argv[opt_ind+2];
    if (!opt_s) intname = argv[opt_ind+3];
	
    if (opt_s)
	sumonly();
    else
    {
	readfiles();
    	suminter();
    	writefiles();
    }
    return (EXIT_OK);
}

