/* ========================== C MeatAxe =============================
   zzztest.c - Test the arithmetic module.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zzztest.c,v 1.2 1997/09/11 15:44:45 gap Exp $
 *
 * $Log: zzztest.c,v $
 * Revision 1.2  1997/09/11 15:44:45  gap
 * New version 2.2.3. AH
 *
 * Revision 2.13  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.12  1994/11/30  10:28:58  mringe
 * ANSI-C, Neue Random-Funktionen.
 *
 * Revision 2.11  1994/11/25  14:08:57  mringe
 * ANSI-C, neue Namen fuer time-Funktionen
 *
 * Revision 2.10  1994/11/25  13:33:55  mringe
 * Neue Random...-Funktionen
 *
 * Revision 2.9  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.8  1994/07/21  11:35:43  mringe
 * *** empty log message ***
 *
 * Revision 2.7  1994/07/20  17:45:11  mringe
 * timings.
 *
 * Revision 2.6  1994/07/07  07:04:54  mringe
 * Zeitangabe korrigiert.
 *
 * Revision 2.5  1994/06/19  16:15:34  mringe
 * Timings.
 *
 * Revision 2.4  1994/06/16  18:40:45  mringe
 * *** empty log message ***
 *
 * Revision 2.3  1994/06/16  17:06:46  mringe
 * Fehlermeldung erweitert.
 *
 * Revision 2.2  1994/06/16  14:17:07  mringe
 * Extract/insert verbessert.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.12  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.11  1993/02/16  18:32:46  mringe
 * string.h und stdio.h werden jetzt in meataxe.h included.
 *
 * Revision 1.10  1992/07/28  13:28:15  mringe
 * *** empty log message ***
 *
 * Revision 1.9  1992/07/28  12:58:41  mringe
 * *** empty log message ***
 *
 * Revision 1.8  1992/07/28  12:51:40  mringe
 * Fixed more bugs.
 *
 * Revision 1.7  1992/07/28  12:50:08  mringe
 * Fixed some bugs with BIG version.
 *
 * Revision 1.6  1992/07/28  09:10:55  mringe
 * Bug in prtab() behoben.
 *
 * Revision 1.5  1992/07/23  06:45:17  mringe
 * zgetgen() durch zgen ersetzt.
 *
 * Revision 1.4  1992/07/15  09:31:47  mringe
 * *** empty log message ***
 *
 * Revision 1.3  1992/07/10  15:23:28  mringe
 * *** empty log message ***
 *
 * Revision 1.2  1992/07/08  08:49:19  mringe
 * Alle Funktionen auch f"ur 'big' version.
 *
 * Revision 1.2  1992/07/08  08:49:19  mringe
 * Alle Funktionen auch f"ur 'big' version.
 *
 * Revision 1.1  1992/05/26  18:28:43  mringe
 * Initial revision
 *
 */


#include <stdlib.h>
#include "meataxe.h"


#if defined (BIG)
#define ISFEL(f) ((f) == 0xffff ||\
	(unsigned short)(f) < (unsigned short)fld-1)
#else
#define ISFEL(f) ((unsigned int)(f) < (unsigned int)fld)
#endif

#define TNOC 20	  /* Number of colums for testing row operations */



/* ------------------------------------------------------------------
   Global variables
   ------------------------------------------------------------------ */

int zzzvers;
long fld;
FEL *ftab;
int opt_t = 0;
int opt_b = 0;

char amsg[] = "zzztest aborted";

static char *helptext[] = {
"SYNTAX",
"    zzztest [-t] [-b] <Field> ...",
"",
"OPTIONS",
"    -t     Print tables",
"    -b     Benchmark",
NULL};

static proginfo_t pinfo =
   { "zzztest", "Test finite field functions", "$Revision: 1.2 $",
      helptext };











#define ERR(msg) err(__LINE__,msg)

void err(int line,char *msg)

{	printf("Line %d: %s\n",line,msg);
	fflush(stdout);
	exit(1);
}




void testitof()

{
    long l, m;

    if (!opt_b && !opt_t)
	printf("Testing int <--> FEL conversion\n");
    if (zitof((long)0) != F_ZERO) ERR("Error 1");
    if (zitof((long)1) != F_ONE) ERR("Error 2");
    if (zftoi(F_ZERO) != 0) ERR("Error 3");
    if (zftoi(F_ONE) != 1) ERR("Error 4");
    for (l = 0; l < fld; ++l)
    {	ftab[l] = zitof(l);
	if (!ISFEL(ftab[l])) ERR("Error 5");
	for (m = 0; m < l; ++m)
	    if (ftab[m] == ftab[l]) ERR("Error 6");
	if (zftoi(ftab[l]) != l) ERR("Error 7");
    }
}



void testfield()

{   FEL a,b,c;
    int ai, bi, ci;
	
    printf("Testing field arithmetic\n");

    /* Test one and zero elements
       -------------------------- */
    for (ai = 0; ai < (int)fld; ++ai)
    {	a = ftab[ai];
	if (zadd(a,F_ZERO) != a)
	{   printf("%d+0=%d\n",(int)a,(int)zadd(a,F_ZERO));
	    ERR(amsg);
	}
	if (zmul(a,F_ONE) != a)
	{	printf("%d*1=%d\n",(int)a,(int)zmul(a,F_ONE));
		ERR(amsg);
	}
    }

    /* Test negative and inverse
       ------------------------- */
    for (ai = 0; ai < (int)fld; ++ai)
    {	a = ftab[ai];
	b = zneg(a);
	if (!ISFEL(b) || zadd(a,b) != F_ZERO) 
	    ERR("Illegal negative");
	if (a != F_ZERO)
	{   b = zinv(a);
	    if (!ISFEL(b) || zmul(a,b) != F_ONE) 
		ERR("Illegal inverse");
	}
    }

    /* Test arithmetic
       --------------- */
    for (ai = 0; ai < (int)fld; ++ai)
    {	a = ftab[ai];
	for (bi = ai; bi < (int)fld; ++bi)
	{   b = ftab[bi];
	    c = zadd(a,b);
	    if (!ISFEL(c)) ERR("zadd() error");
	    if (c != zadd(b,a)) ERR("'+' not commutative");
	    c = zmul(a,b);
	    if (!ISFEL(c)) ERR("zmul() error");
	    if (c != zmul(b,a)) ERR("'*' not commutative");
	    
	    for (ci = 0; ci < (int)fld; ++ci)
	    {	c = ftab[ci];
		if (zadd(a,zadd(b,c)) != zadd(zadd(a,b),c))
		    ERR("'+' not associative");
		if (zmul(a,zmul(b,c)) != zmul(zmul(a,b),c))
		    ERR("'*' not associative");
		if (zmul(a,zadd(b,c)) != zadd(zmul(a,b),zmul(a,c)))
		    ERR("a*(b+c) != a*b+a*c");
	    }
	}
    }
}


void testgen()
{
    FEL a, b;
    int i;

    MESSAGE(0,("Testing zgen\n"));
    a = zgen;
    MESSAGE(1,("Generator = %d\n",(int) zgen));
    b = a;
    for (i = 1; i < (int)fld-1; ++i)
    {	
	if (b == F_ONE)
	{
	    ERR("Generator test failed");
	    break;
	}
	b = zmul(a,b);
    }
    if (b != F_ONE)
    	ERR("Generator test failed");
}


void testextractinsert()

{
    long noc;

    printf("Testing zinsert()/zextract()\n");
    for (noc = 1; noc < 35; ++noc)
    {
	PTR x;
	unsigned i;

	zsetfield(fld);
	zsetlen(noc);
	x = zalloc((long)1);
  	for (i = 0; i < 50; ++i)
	{
	    long col;
	    RandInit(i);
	    for (col = 1; col <= noc; ++col)
		zinsert(x,col,ftab[RandInt(fld)]);
	    RandInit(i);
	    for (col = 1; col <= noc; ++col)
		if (zextract(x,col) != ftab[RandInt(fld)])
		    ERR("insert/extract test failed");
	}
	free(x);
    }
}


void testfindmark()

{	PTR x;
	FEL a;
	long i;

	zsetfield(fld);
	zsetlen(TNOC);
	printf("Testing zfindpiv()\n");
	x = zalloc((long)1);

	for (i = 1; i <= TNOC; ++i)
		zinsert(x,i,ftab[i%(fld-1)+1]);
	for (i = 1; i <= TNOC; ++i)
	{	if (zfindpiv(x,&a) != i) ERR("Error 1");
		if (a != ftab[i%(fld-1)+1]) ERR("Error 2");
		zinsert(x,i,F_ZERO);
	}
	free(x);
}


void testaddrow()

{	PTR x, y;
	int a, b;
	FEL c;
	long i;

	zsetfield(fld);
	zsetlen(TNOC);
	printf("Testing zaddrow()\n");
	x = zalloc((long)1);
	y = zalloc((long)1);
	for (a = 0; a < (int)fld; ++a)
	    for (b = 0; b < (int)fld; ++b)
	    {	for (i = 1; i <= TNOC; ++i)
		{	zinsert(x,i,ftab[a]);
			zinsert(y,i,ftab[b]);
		}
		zaddrow(x,y);
		c = zadd(ftab[a],ftab[b]);
                for (i = 1; i <= TNOC; ++i)
		{	if (zextract(x,i) != c)
				ERR("Test failed");
		}
	    }
	free(x);
	free(y);
}



void testrowops()

{	PTR m1,m2,m3;
	long i, k;
	int d;

	zsetfield(fld);
	zsetlen(TNOC);
	printf("Testing row operations\n");
	m1 = zalloc((long)1);
	m2 = zalloc((long)1);
	m3 = zalloc((long)1);

	/* Test zmoverow() and zcmprow()
	   ----------------------------- */
	for (i = 1; i <= TNOC; ++i)
		zinsert(m1,i,ftab[i%fld]);
	for (i = 1; i <= TNOC; ++i)
	{	zmoverow(m2,m1);
		if (zcmprow(m2,m1)) ERR("Error 1");
		if (memcmp(m2,m1,zsize((long)1))) ERR("Error 2");
		zinsert(m1,i,ftab[(i+1)%fld]);
		if (!zcmprow(m2,m1)) ERR("Error 3");
	}

	/* Test zaddrow() and zaddmulrow()
	   ------------------------------- */
	for (d = 0; d < 6; ++d)
	{	for (i = 1; i <= TNOC; ++i)
		{	zinsert(m1,i,ftab[i%fld]);
			zinsert(m2,i,ftab[(i+d)%fld]);
		}
		zmoverow(m3,m1);
		zaddrow(m3,m2);
		for (i = 1; i <= TNOC; ++i)
		{	if (zextract(m3,i) !=
			    zadd(ftab[i%fld],ftab[(i+d)%fld]))
				ERR("Error 4");
		}
		for (k = 0; k < fld; ++k)
		{	zmoverow(m3,m1);
			zaddmulrow(m3,m2,ftab[k]);
			for (i = 1; i <= TNOC; ++i)
			{	if (zextract(m3,i) != zadd(ftab[i%fld],
				    zmul(ftab[(i+d)%fld],ftab[k])))
					ERR("Error 5");
			}
		}
	}
}


void prtables()
{	int a, b;

	printf(" + ");
	for (a = 0; a < (FEL)fld; ++a)
		printf("%3d", a);
	printf("\n");
	for (a = 0; a < (int)fld; ++a)
	{       printf("%3d",a);
		for (b = 0; b < (int)fld; ++b)
		    printf("%3ld",zftoi(zadd(ftab[a],ftab[b])));
		printf("\n");
	}

	printf("\n * ");
	for (a = 0; a < (FEL)fld; ++a)
		printf("%3d", a);
	printf("\n");
	for (a = 0; a < (int)fld; ++a)
	{       printf("%3d",a);
		for (b = 0; b < (int)fld; ++b)
			printf("%3ld",zftoi(zmul(ftab[a],ftab[b])));
		printf("\n");
	}

}



static void timing()

{
    static int firsttime = 1;
    long noc = 2000;
    long count = 0;
    PTR x, y, z;
    long i, t;

    if (firsttime)
    {
	firsttime = 0;
	printf("ZZZ version: %s\n",zzzversion);
	printf("ZZZ compile info: %s\n",zzz_cc);
	printf("Field  Integer Extr/Ins    RowOp\n");
    }
    printf("%5ld",fld);

    zsetlen(noc);
    x = zalloc(1);
    y = zalloc(1);
    z = zalloc(1);
    RandInit(1);
    for (i = 1; i <= noc; ++i)
    {
	zinsert(x,i,ftab[RandInt(fld)]);
	zinsert(y,i,ftab[RandInt(fld)]);
    }

    for (count = 0, t = STimeUsed(); STimeUsed() - t < 50; ++count)
    {
	register long i, k;
	for (i = k = 0; i < 10000; ++i) k += i / 17;
    }
    t = STimeUsed() - t;
    printf("%9ld",count); 
    fflush(stdout);


    /* extract/insert */
    for (count = 0, t = STimeUsed(); STimeUsed() - t < 100; ++count)
    {
        for (i = 1; i <= noc; ++i)
	    zinsert(z,i,zextract(y,i));
    }
    printf("%9ld",count); 
    fflush(stdout);


    for (count = 0, t = STimeUsed(); STimeUsed() - t < 100; ++count)
    {
	zaddrow(z,x);
	zaddmulrow(z,y,ftab[RandInt(fld)]);
    }
    t = STimeUsed() - t;
    printf("%9ld",count/10); 
    printf("\n");
}



static void zzztest()

{
    MESSAGE(0,("Testing with GF(%ld)\n",fld));
    MESSAGE(0,("ZZZ version: %d\n",zzzvers));
    testfield();
    testgen();
    testextractinsert();
    testfindmark();
    testaddrow();
    testrowops();
} 


int main(argc, argv)

int argc;
char *argv[];

{   int k;

    zzzvers = mtxinit();
    initargs(argc, argv, &pinfo);
    while ((k = zgetopt("bt")) != OPT_END)
    {   switch (k)
	{ 
	    case 'b': opt_b = 1; break;
	    case 't': opt_t = 1; break;
	}
    }
    
    for (; opt_ind < argc; ++opt_ind)
    {
	fld = atol(argv[opt_ind]);
	zsetfield(fld);
	zsetlen(TNOC);
	ftab = NALLOC(FEL,fld);
    	testitof();

    	if (opt_t)
	    prtables();
	else if (opt_b)
	    timing();
	else
	    zzztest();
	free(ftab);
    }
    return 0;
}




