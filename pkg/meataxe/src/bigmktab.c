/* ========================== C MeatAxe =============================
   This program generates the MEAT-AXE table files for large fields
   (q < 65535).

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */

/* $Id: bigmktab.c,v 1.2 1997/09/11 15:42:36 gap Exp $
 *
 * $Log: bigmktab.c,v $
 * Revision 1.2  1997/09/11 15:42:36  gap
 * New version 2.2.3. AH
 *
 * Revision 2.3  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.2  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.1  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.6  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.5  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.4  1992/07/28  12:52:28  mringe
 * Fixed some bugs.
 *
 * Revision 1.3  1992/06/30  07:56:14  mringe
 * Conway-Polynome von GAP eingetragen
 *
 * Revision 1.2  1992/06/27  17:21:11  mringe
 * Einige kleinere Tippfehler behoben.
 *
 * Revision 1.1  1992/05/25  18:07:50  mringe
 * Initial revision
 *
 */



#if !defined(BIG)
#define BIG
#endif

#include "meataxe.h"

#define MAXGRAD 25		/* Maximum degree of polynomials */

typedef unsigned short POLY[MAXGRAD+1];



/* -----------------------------------------------------------------
   Global data
   ----------------------------------------------------------------- */

char rcsid[] = "$Revision: 1.2 $";
char rcsdate[] = "$Date: 1997/09/11 15:42:36 $";

unsigned short *inc;
unsigned short Minusone;

unsigned short 	P;	/* Characteristic of the field */
unsigned short	Q;	/* Order of the field */
unsigned short	N;	/* Q = P^N */
FILE	*fd;		/* Output file */
unsigned short Gen;	/* Generator */

POLY irred;		/*  Polynomial which defines the field */
unsigned short *indx;	/*  Index i of a field element g,g=X^i */
unsigned short ppwr[MAXGRAD+1];
unsigned short ppindex[MAXGRAD+1];



/* -----------------------------------------------------------------
   printpol() - Print a polynomial
   ----------------------------------------------------------------- */

void printpol(a)

POLY a;

{	int i,flag = 0;

	for (i = MAXGRAD; i >= 0; i--)
		if (a[i]!=0)
		{	if (flag) printf("+");
			if (a[i] != 1) printf("%d",(int)a[i]);
			printf("x^%d",i);
			flag=1;
        	}
        printf("\n");
}


/* -----------------------------------------------------------------
   number() - Convert polynomial to number
   ----------------------------------------------------------------- */

unsigned short number(a)
POLY a;			/* polynomial a0..an -> number sum(ai*p^i) */

{	unsigned short k;
	int i;

	k = 0;
	for (i = MAXGRAD; i >= 0; i--)
		k = k * (unsigned short) P + a[i];
	return (k);
}


/* -----------------------------------------------------------------
   polmultx() - Multiply a polynomial by X.
   ----------------------------------------------------------------- */

void polmultx(a)
POLY a;

{	int i;

	for (i = MAXGRAD; i > 0; --i)
		a[i] = a[i-1];
	a[0] = 0;
}


/* -----------------------------------------------------------------
   polmod() - Reduce the polynomial a modulo b. b ist assumed to
	be normalized.
   ----------------------------------------------------------------- */

void polmod(a,b)
POLY a,b;

{	int i, l, dl, f;

	/* l= index of leading coeff. in b (must be 1) */
	for (l = MAXGRAD; b[l]==0 && l>0; l--);
	for (dl = MAXGRAD; dl>=l; dl--)
	{	f = (int) a[dl];
		if (f == 0) continue;
		f = (int)P - f;
		for (i = 0; i <= l; ++i)
			a[i+dl-l] =
				(unsigned short) ((f*b[i] + a[i+dl-l]) % (int)P);
        }
}


/* -----------------------------------------------------------------
   testprim() - Test for primitivity.
   ----------------------------------------------------------------- */

void testprim()

{	unsigned short i, *a;

	a = malloc(sizeof(short)*(size_t)Q);
	memset(a,0,sizeof(short)*(size_t)Q);
	for (i = 1; i < Q; i++)
		++a[indx[i]];
	for (i = 0; i < Q-1; i++)
        	if(a[i] != 1)
		{	fprintf(stderr,
			 "bigmktab: Polynome is not primitive!\n");
			exit(1);
		}
	free(a);
}


/* -----------------------------------------------------------------
   initarith() - Initialize index table.
   ----------------------------------------------------------------- */

void initarith()

{	int i, elem;
	POLY a;

	memset(a,0,sizeof(POLY));
	a[0] = 1;			/* a=X^0 (=generator) */
	indx[0] = (unsigned short) 0xFFFF;
	for (i = 0; i < (int)Q-1; i++)
	{	elem = number(a);
		indx[elem] = (unsigned short) i;
		polmultx(a);
		polmod(a,irred);
        }
	testprim();
	Gen = P;
}

/* -----------------------------------------------------------------
   initarithP() - Initialize index table for prime field
   ----------------------------------------------------------------- */

void initarithP()

{	unsigned short i;
	unsigned short a, gen, x;


	/* Find generator
	   -------------- */
	for (gen = 1; gen < P; ++gen)
	{	x = gen;
		for (i = 1; x != 1; ++i)
			x = (unsigned short) (((long)x * gen) % P);
		if (i == P - 1) break;
	}
	Gen = gen;

	/* Make index table
	   ---------------- */
	indx[0] = (unsigned short) 0xFFFF;
	a = 1;
	for (i = 0; i < P-1; i++)
	{	indx[a] = i;
		a = (unsigned short) (((long) a * gen) % P);
        }
	testprim();
	irred[0] = 0;
	irred[1] = 1;
}


/* ------------------------------------------------------------------
   comptables() - Compute tables
   ------------------------------------------------------------------ */

void comptables()

{	unsigned short i, j;

	for (i = 1; i <= Q-1; ++i)
	{	j = (int)((i%P)==P-1 ? i+1-P : i+1);
		if (j == 0)
			Minusone = indx[i];
		inc[indx[i]] = indx[j];
	}
	for (i = 0; i < N; ++i)
		ppindex[i] = indx[ppwr[i]];
}


/* ------------------------------------------------------------------
   writeheader() - Write file header, select polynomial, and
	initialize tables.
   ------------------------------------------------------------------ */

void writeheader()

{	char fname[50];
	char rev[10], date[20];
	unsigned short header[5];

	sprintf(fname,"p%5.5ld.zzz",(long) Q);
	fd = SFOpen(fname,FM_CREATE);
	if (fd == NULL)
	{	perror(fname);
		fprintf(stderr,"bigmktab: Cannot create table file!\n");
		exit(EXIT_ERR);
	}
	header[0] = ZZZVERSION;
	header[1] = P;
	header[2] = Q;
	header[3] = N;
	header[4] = Gen;

	sscanf(rcsid,"$Revision:%s",rev);
	sscanf(rcsdate,"$Date:%s",date);
	printf("BIGMKTAB Revision %s (%s)\n",rev,date);

	printf("ZZZ version : %hd\n",header[0]);
	printf("Field order : %hd=%hd^%hd\n",header[2],header[1],N);
	if (P != Q)
	{	printf("Polynomial  : ");
		printpol(irred);
		printf("Generator   : x\n");
	}
	else
		printf("Generator   : %hd\n",Gen);
	if (fwrite(header,sizeof(short),5,fd) != 5)
	{   perror(fname);
	    fprintf(stderr,"bigmktab: Error writing file header!\n");
	    exit(EXIT_ERR);
	}
}


/* -----------------------------------------------------------------
   getpol() - Get the polynomial which defines the field
   ----------------------------------------------------------------- */

typedef struct
	{	unsigned short p, n;
		POLY pol;
	}
	poltab_t;


/* Conway polynomials. Created from GAP's Pols[] array */

poltab_t poltab[] =
{   
    {  2, 2,{1,1,1}},				/* GF(4) */
    {  2, 3,{1,1,0,1}},				/* GF(8) */
    {  2, 4,{1,1,0,0,1}},			/* GF(16) */
    {  2, 5,{1,0,1,0,0,1}},			/* GF(32) */
    {  2, 6,{1,1,0,1,1,0,1}},			/* GF(64) */
    {  2, 7,{1,1,0,0,0,0,0,1}},			/* GF(128) */
    {  2, 8,{1,0,1,1,1,0,0,0,1}},		/* GF(256) */
    {  2, 9,{1,0,0,0,1,0,0,0,0,1}},		/* GF(512) */
    {  2,10,{1,1,1,1,0,1,1,0,0,0,1}},		/* GF(1024) */
    {  2,11,{1,0,1,0,0,0,0,0,0,0,0,1}},		/* GF(2048) */
    {  2,12,{1,1,0,1,0,1,1,1,0,0,0,0,1}},	/* GF(4096) */
    {  2,13,{1,1,0,1,1,0,0,0,0,0,0,0,0,1}},	/* GF(8192) */
    {  2,14,{1,0,0,1,0,1,0,1,0,0,0,0,0,0,1}},	/* GF(16384) */
    {  2,15,{1,0,1,0,1,1,0,0,0,0,0,0,0,0,0,1}},	/* GF(32768) */
    {  2,16,{1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,1}},/* GF(65536) */
    {  3, 2,{2,2,1}},				/* GF(9) */
    {  3, 3,{1,2,0,1}},				/* GF(27) */
    {  3, 4,{2,0,0,2,1}},			/* GF(81) */
    {  3, 5,{1,2,0,0,0,1}},			/* GF(243) */
    {  3, 6,{2,2,1,0,2,0,1}},		 	/* GF(729) */
    {  3, 7,{1,0,2,0,0,0,0,1}},			/* GF(2187) */
    {  3, 8,{2,2,2,0,1,2,0,0,1}},		/* GF(6561) */
    {  3, 9,{1,1,2,2,0,0,0,0,0,1}},		/* GF(19683) */
    {  3,10,{2,1,0,0,2,2,2,0,0,0,1}},		/* GF(59049) */
    {  5, 2,{2,4,1}},				/* GF(25) */
    {  5, 3,{3,3,0,1}},				/* GF(125) */
    {  5, 4,{2,4,4,0,1}},			/* GF(625) */
    {  5, 5,{3,4,0,0,0,1}},		   	/* GF(3125) */
    {  5, 6,{2,0,1,4,1,0,1}},		 	/* GF(15625) */
    {  7, 2,{3,6,1}},				/* GF(49) */
    {  7, 3,{4,0,6,1}},		       		/* GF(343) */
    {  7, 4,{3,4,5,0,1}},		     	/* GF(2401) */
    {  7, 5,{4,1,0,0,0,1}},		   	/* GF(16807) */
    { 11, 2,{2,7,1}},				/* GF(121) */
    { 11, 3,{9,2,0,1}},				/* GF(1331) */
    { 11, 4,{2,10,8,0,1}},			/* GF(14641) */
    { 13, 2,{2,12,1}},		        	/* GF(169) */
    { 13, 3,{11,2,0,1}},	      		/* GF(2197) */
    { 13, 4,{2,12,3,0,1}},		    	/* GF(28561) */
    { 17, 2,{3,16,1}},		        	/* GF(289) */
    { 17, 3,{14,1,0,1}},	      		/* GF(4913) */
    { 19, 2,{2,18,1}},		        	/* GF(361) */
    { 19, 3,{17,4,0,1}},	      		/* GF(6859) */
    { 23, 2,{5,21,1}},		        	/* GF(529) */
    { 23, 3,{18,2,0,1}},	      		/* GF(12167) */
    { 29, 2,{2,24,1}},		        	/* GF(841) */
    { 29, 3,{27,2,0,1}},	      		/* GF(24389) */
    { 31, 2,{3,29,1}},		        	/* GF(961) */
    { 31, 3,{28,1,0,1}},	      		/* GF(29791) */
    { 37, 2,{2,33,1}},		       	 	/* GF(1369) */
    { 37, 3,{35,6,0,1}},	      		/* GF(50653) */
    { 41, 2,{6,38,1}},		       		/* GF(1681) */
    { 43, 2,{3,42,1}},				/* GF(1849) */
    { 47, 2,{5,45,1}},				/* GF(2209) */
    { 53, 2,{2,49,1}},				/* GF(2809) */
    { 59, 2,{2,58,1}},				/* GF(3481) */
    { 61, 2,{2,60,1}},				/* GF(3721) */
    { 67, 2,{2,63,1}},				/* GF(4489) */
    { 71, 2,{7,69,1}},				/* GF(5041) */
    { 73, 2,{5,70,1}},				/* GF(5329) */
    { 79, 2,{3,78,1}},				/* GF(6241) */
    { 83, 2,{2,82,1}},				/* GF(6889) */
    { 89, 2,{3,82,1}},				/* GF(7921) */
    { 97, 2,{5,96,1}},				/* GF(9409) */
    {101, 2,{2,97,1}},				/* GF(10201) */
    {103, 2,{5,102,1}},		       		/* GF(10609) */
    {107, 2,{2,103,1}},		       		/* GF(11449) */
    {109, 2,{6,108,1}},				/* GF(11881) */
    {113, 2,{3,101,1}},				/* GF(12769) */
    {127, 2,{3,126,1}},				/* GF(16129) */
    {131, 2,{2,127,1}},				/* GF(17161) */
    {137, 2,{3,131,1}},				/* GF(18769) */
    {139, 2,{2,138,1}},				/* GF(19321) */
    {149, 2,{2,145,1}},				/* GF(22201) */
    {151, 2,{6,149,1}},				/* GF(22801) */
    {157, 2,{5,152,1}},				/* GF(24649) */
    {163, 2,{2,159,1}},				/* GF(26569) */
    {167, 2,{5,166,1}},				/* GF(27889) */
    {173, 2,{2,169,1}},				/* GF(29929) */
    {179, 2,{2,172,1}},				/* GF(32041) */
    {181, 2,{2,177,1}},				/* GF(32761) */
    {191, 2,{19,190,1}},			/* GF(36481) */
    {193, 2,{5,192,1}},				/* GF(37249) */
    {197, 2,{2,192,1}},				/* GF(38809) */
    {199, 2,{3,193,1}},				/* GF(39601) */
    {211, 2,{2,207,1}},				/* GF(44521) */
    {223, 2,{3,221,1}},				/* GF(49729) */
    {227, 2,{2,220,1}},				/* GF(51529) */
    {229, 2,{6,228,1}},				/* GF(52441) */
    {233, 2,{3,232,1}},				/* GF(54289) */
    {239, 2,{7,237,1}},				/* GF(57121) */
    {241, 2,{7,238,1}},				/* GF(58081) */
    {251, 2,{6,242,1}},				/* GF(63001) */
    { 0,0,{0}}					/* End of table */
};

void getpol()

{	poltab_t *x;
	for (x = poltab; x->p != 0; ++x)
		if (x->p == P && x->n == N)
			break;
	if (x->p != 0)
		memcpy(irred,x->pol,sizeof(POLY));
	else
	{	fprintf(stderr,"bigmktab: Polynomial not found\n");
		exit(EXIT_ERR);
	}
}



/* -----------------------------------------------------------------
   init() - Initialize
   ----------------------------------------------------------------- */

void init(l)

long l;

{	unsigned short q, d;

	if (l < 2 || l > 65535)
	{	fprintf(stderr,
		 "bigmktab: Field order out of range (2-65535)\n");
		exit(EXIT_ERR);
	}
	Q = (unsigned short) l;
	q = Q;
	for (d = 2; q % d != 0; ++d);
	P = d;
	for (N = 0; (q % d) == 0; ++N)
        	q /= d;
	if (q != 1)
	{	fprintf(stderr,"bigmktab: Think first\n");
        	exit(EXIT_ERR);
        }
	ppwr[0] = 1;
	for (d = 1; d < N; ++d)
		ppwr[d] = ppwr[d-1] * P;

	inc = malloc((size_t)Q * sizeof(short));
	indx = malloc((size_t)Q * sizeof(short));
	if (N != 1)
	{	getpol();
		initarith();
	}
        else
	{	initarithP();
	}
}



/* -----------------------------------------------------------------
   writetables()
   ----------------------------------------------------------------- */

void writetables()

{	printf("Writing tables to `");
	printf(TABFILENAME,(long)Q);
	printf("'...");
	fflush(stdout);
	if (fwrite(irred,sizeof(short),(size_t)N+1,fd) != (size_t)N+1 ||
	  fwrite(ppwr,sizeof(short),(size_t) N,fd) != (size_t)N ||
	  fwrite(ppindex,sizeof(short),(size_t) N,fd) != (size_t)N ||
	  fwrite(&Minusone,sizeof(short),1,fd) != 1 ||
	  fwrite(inc,sizeof(short),(size_t)Q-1,fd) != (size_t) Q-1)
	{	fprintf(stderr,"bigmktab: Error writing tables\n");
		perror("fwrite");
		exit(EXIT_ERR);
	}
	fclose(fd);
	printf("Ok\n");
}


/* -----------------------------------------------------------------
   main()
   ----------------------------------------------------------------- */

int main(argc, argv)
int argc;
char *argv[];

{	long int l;

	if (argc != 2 || sscanf(argv[1],"%ld",&l) != 1)
	{	fprintf(stderr,"Usage: bigmktab <field order>\n");
		return(EXIT_ERR);
	}
	init(l);
	writeheader();
	comptables();
	writetables();
	return(0);
}


