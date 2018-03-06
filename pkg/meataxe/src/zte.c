/* ========================== C MeatAxe =============================
   zte.c - Tensor two matrices or permutations.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zte.c,v 1.2 1997/09/11 15:44:34 gap Exp $
 *
 * $Log: zte.c,v $
 * Revision 1.2  1997/09/11 15:44:34  gap
 * New version 2.2.3. AH
 *
 * Revision 2.9  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.8  1994/09/19  15:46:02  mringe
 * Bug bei Permutationen behoben.
 *
 * Revision 2.7  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.6  1994/03/13  13:27:01  mringe
 * Maschinenunabhaengiges Format fuer Permutationen.
 *
 * Revision 2.5  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.4  1993/12/08  11:33:02  mringe
 * Neue CPU time - Funktionen.
 *
 * Revision 2.3  1993/10/26  10:47:35  mringe
 * Compiler Warnings.
 *
 * Revision 2.2  1993/10/21  21:58:39  mringe
 * Permutationen.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.8  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.7  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.6  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.5  1993/08/05  15:52:49  mringe
 * File header
 *
 * Revision 1.4  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.3  1993/07/17  12:52:05  mringe
 * Message lib, help, time check.
 *
 * Revision 1.2  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.1  1992/05/26  18:17:04  mringe
 * Initial revision
 *
 */


#include <stdlib.h>
#include "meataxe.h"


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

long fl, nor1, noc1, nor2, noc2;
char *aname, *bname, *cname;	/* File names */
FILE *afile, *bfile;

static char *helptext[] = {
"SYNTAX",
"    zte [-T <MaxTime>] <A> <B> <Result>",
"",
"OPTIONS",
"    -T <NSeconds>    Set CPU time limit in seconds.",
"",
"FILES",
"    <A>       i  A matrix or permutation",
"    <B>       i  A matrix or permutation",
"    <Result>  o  The tensor product (A*B)",
NULL};

static proginfo_t pinfo =
   { "zte", "Tensor Product", "$Revision: 1.2 $", helptext };





/* ------------------------------------------------------------------
   tensormatrices() - Calculate the tensor product of two matrices
   ------------------------------------------------------------------ */

static void tensormatrices()

{
    PTR m1, m2, m3, x1, x2;
    long nor3, noc3;
    long j1,j2,j3,i1,i2;
    FEL f1,f2,f3;
    FILE *f;

    zsetfield(fl);
    zsetlen(noc1);
    m1 = zalloc(nor1);
    if (zreadvec(afile,m1,nor1) != nor1) errexit(-1,aname);
    zsetlen(noc2);
    m2 = zalloc(nor2);
    if (zreadvec(bfile,m2,nor2) != nor2) errexit(-1,bname);
    nor3 = nor1 * nor2;
    noc3 = noc1 * noc2;
    MESSAGE(1,("Computing matrix tensor product:"));
    MESSAGE(1,("(%ld,%ld)*(%ld,%lx)=(%ld,%ld)\n", nor1,noc1,nor2,noc2,
    nor3,noc3));
    zsetlen(noc3);
    m3 = zalloc((long)1);
    if ((f = zwritehdr(cname,fl,nor3,noc3)) == NULL)
	errexit(-1,cname);
    x1 = m1;
    for (i1 = 1; i1 <= nor1; ++i1)
    {
		x2 = m2;
		for (i2 = 1; i2 <= nor2; ++i2)
		{	j3 = 1;
			zmulrow(m3,F_ZERO);
			for (j1 = 1; j1 <= noc1; ++j1)
			{	f1 = zextract(x1,j1);
				for (j2 = 1; j2 <= noc2; ++j2)
				{	f2 = zextract(x2,j2);
					f3 = zmul(f1,f2);
					zinsert(m3,j3,f3);
					++j3;
				}
			}
			zsetlen(noc3);
			zwritevec(f,m3,(long)1);
			zsetlen(noc2);
			zadvance(&x2,(long)1);
		}
		zsetlen(noc1);
		zadvance(&x1,(long)1);
    }
    fclose(f);
}


/* ------------------------------------------------------------------
    tensorperms() - Calculate the tensor product of two permutations.
   ------------------------------------------------------------------ */

static void tensorperms(void)

{
    long nor3;
    long i, k, ii, kk, l;
    long *m1, *m2;
    FILE *f;
    long *p1, *p2;

    /* Read permutations
       ----------------- */
    m1 = (long *) malloc(sizeof(long) * nor1);
    m2 = (long *) malloc(sizeof(long) * nor2);
    if (m1 == NULL || m2 == NULL) errexit(ERR_NOMEM,"zte");
    p1 = m1-1;
    p2 = m2-1;
    if (zreadlong(afile,m1,nor1) != nor1)
	errexit(-1,aname);
    if (zreadlong(bfile,m2,nor2) != nor2)
	errexit(-1,bname);

    /* Open output file
       ---------------- */
    nor3 = nor1 * nor2;
    if ((f = zwritehdr(cname,(long)-1,nor3,(long)1)) == NULL)
	errexit(-1,cname);
    MESSAGE(1,("Computing permutation tensor product:"));
    MESSAGE(1,("%ld*%ld=%ld\n", nor1,nor2,nor3));

    /* Calculate the tensor product
       ---------------------------- */
    for (i = 1; i <= nor1; ++i)
    {
	ii = p1[i] - 1;
	for (k = 1; k <= nor2; ++k)
	{
	    kk = p2[k];
	    l = ii * nor2 + kk;
	    if (zwritelong(f,&l,1) != 1)
		errexit(-1,cname);
	}
    }
    fclose(f);
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char *argv[])

{
    long fl2;
    char buf[200];

    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while (zgetopt("") != OPT_END);
    if (opt_ind != argc-3) errexit(ERR_BADUSAGE,"zte");
    aname = argv[opt_ind];
    bname = argv[opt_ind+1];
    cname = argv[opt_ind+2];

    if ((afile = zreadhdr(aname,&fl,&nor1,&noc1)) == NULL)
	errexit(-1,aname);
    if ((bfile = zreadhdr(bname,&fl2,&nor2,&noc2)) == NULL)
	errexit(-1,bname);
    if (fl2 != fl)
    {
	sprintf(buf,"%s and %s",aname,bname);
	errexit(ERR_INCOMPAT,buf);
    }
    if (fl == -1)
	tensorperms();
    else
	tensormatrices();
    fclose(afile);
    fclose(bfile);
    return EXIT_OK;
}

