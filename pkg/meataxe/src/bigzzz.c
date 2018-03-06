/* ========================== C MeatAxe =============================
   Finite field arithmetic and common functions. `Big' version for
   large fields (q<=65536).

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: bigzzz.c,v 1.2 1997/09/11 15:42:38 gap Exp $
 *
 * $Log: bigzzz.c,v $
 * Revision 1.2  1997/09/11 15:42:38  gap
 * New version 2.2.3. AH
 *
 * Revision 2.5  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.4  1994/02/15  10:28:33  mringe
 * MSDOS_BCC entfernt.
 *
 * Revision 2.3  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
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
 * Revision 1.16  1993/10/05  23:35:02  mringe
 * zerrno eliminiert.
 *
 * Revision 1.15  1993/08/06  14:04:46  mringe
 * SYS_MSDOS durch OS_MSDOS ersetzt
 *
 * Revision 1.14  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.13  1993/07/29  14:47:00  mringe
 * OS_xxxx Symbole.
 *
 * Revision 1.12  1993/07/23  13:46:27  mringe
 * OS-Symbole neu (SYS_xxx)
 *
 * Revision 1.11  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.10  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.9  1992/11/04  09:10:16  mringe
 * Neu: zpermrow()
 *
 */


#if !defined(BIG)
#define BIG 1
#endif

#include "meataxe.h"

#define MAXPWR 20


/* ------------------------------------------------------------------
   Data
   ------------------------------------------------------------------ */

FEL minusone;			/* -1 */
FEL *inc = NULL;		/* a+1 = inc[a] */

static unsigned short P = 0;		/* Characteristic */
static unsigned short Q = 0;		/* Field order */
static unsigned short Q1 = 0;		/* Q-1 */
static unsigned short N;
static unsigned short Gen;
static unsigned short ppwr[MAXPWR];
static unsigned short ppindex[MAXPWR];
static unsigned short poly[MAXPWR];

static long NOC = 0;            /* No. of columns per row */
static long IPR = 0;            /* No. of long ints per row */
static long rowsize = 0;        /* Row size (bytes)=IPR*sizeof(long) */

static char iname[MAXFILES+1][50], oname[MAXFILES+1][50];
static FILE *ifile[MAXFILES+1], *ofile[MAXFILES+1];


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

long znoc, zfl = -1, zchar;
FEL zgen = 1;


/* ------------------------------------------------------------------
   Argument checking macros
   ------------------------------------------------------------------ */

#if defined(DEBUG)

#define CHECKRANGE(x,lo,hi) if ((x)<(lo)||(x)>(hi)) {\
	fprintf(stderr,"%ld <= %ld <= %ld ?\n",(long)(lo),\
		(long)(x),(long)(hi));\
	fatal("RANGE CHECK ERROR",__LINE__);}
#define CHECKFILE(x)  CHECKRANGE(file,0,MAXFILES)
#define CHECKCOL(x)  CHECKRANGE(x,1,NOC)
#define CHECKFEL(x) { \
	if ((x) != 0xFFFF && ((x) > Q-2)) \
		fatal("range check error",__LINE__); \
	}

#else

#define CHECKRANGE(x,lo,hi)
#define CHECKCOL(x)
#define CHECKFILE(x)
#define CHECKFEL(x)

#endif



/* ------------------------------------------------------------------
   fatal() - fatal error
   ------------------------------------------------------------------ */

#define FATAL(msg) fatal(msg,__LINE__)

#if defined(__STDC__)
static void fatal(char *msg, int line);
#endif

static void fatal(msg, line)
char *msg;
int line;

{	fprintf(stderr,"ZZZ ERROR (%s, line %d) - %s\n",__FILE__,
		line,msg);
	exit(EXIT_ERR);
}



/* ------------------------------------------------------------------
   mtxinit()  -  Initialize

   return value: version number
   ------------------------------------------------------------------ */

int mtxinit()

{	int i;
	char c[10];

	zchar = 0;
	zfl = 0;
	for (i = 0; i <= MAXFILES; ++i)
	{	ifile[i] = ofile[i] = NULL;
		sprintf(c,"G%d",i);
		zname(i,0,c);
		sprintf(c,"P%d",i);
		zname(i,1,c);
	}
	return (ZZZVERSION);
}

/* ------------------------------------------------------------------
   zsetlen()  -  Set field parameter and row length
   ------------------------------------------------------------------ */

void zsetlen(field, ncols)
long field;
long ncols;

{	static char filename[50];
	FILE *fd;
	unsigned short info[5];

	if (field != zfl)
	{	sprintf(filename,"p%5.5ld.zzz",field);
		if ((fd = SFOpen(filename,FM_READ|FM_LIB)) == NULL ||
		    fread((char *)info,sizeof(short),5,fd) != 5)
		{	perror(filename);
			FATAL("Error opening table file");
		}
		P = info[1];
		Q = info[2];
		N = info[3];
		Q1 = Q - 1;
		Gen = info[4];
		zfl = (long) Q;
		zchar = (long) P;
		if (Q == 2) zgen = 0; else zgen = 1;
		if ((info[0]) != ZZZVERSION || N >= MAXPWR || Q < 2 ||
		    P < 2 || P > Q || Q % P != 0)
			FATAL("ERROR IN TABLE FILE HEADER");
		if (inc != NULL) free(inc);
		inc = malloc((size_t) (Q-1) * sizeof(FEL));
		if (
		    fread(poly,sizeof(short),(size_t) N+1,fd) != (size_t)N+1 ||
		    fread(ppwr,sizeof(short),(size_t) N,fd) != (size_t)N ||
		    fread(ppindex,sizeof(short),(size_t) N,fd) != (size_t)N ||
		    fread(&minusone,sizeof(short),1,fd) != 1 ||
		    fread(inc,sizeof(short),(size_t) Q-1,fd) != (int) Q-1
		   )
		{	perror(filename);
			exit(EXIT_ERR);
		}
		if (poly[N] != 1) FATAL("ERROR IN TABLE FILE");
		fclose(fd);
	}
	NOC = ncols;       			/* columns per row */
	IPR = (NOC*sizeof(FEL)-1)/(sizeof(long))+1;/* ints per row */
	rowsize = IPR * sizeof(long);	/* total bytes per row */
}


/* ------------------------------------------------------------------
   zadd(), zsub(), zmul(), zdiv() - Field arithmetic
   ------------------------------------------------------------------ */

#if defined (__STDC__)
FEL zadd(FEL a, FEL b)
#else
FEL zadd(a, b)
FEL a, b;
#endif

{	register FEL x;

	CHECKFEL(a);
	CHECKFEL(b);

	if (b == F_ZERO) return a;
	if (a == F_ZERO) return b;
	if (a >= b)
	{	x = inc[a-b];
		if (x == F_ZERO) return F_ZERO;
		if ((x += b) >= Q1) x -= Q1;
		return x;
	}
	x = inc[b-a];
	if (x == F_ZERO) return F_ZERO;
	if ((x += a) >= Q1) x -= Q1;
	return x;
}


#if defined (__STDC__)
FEL zsub(FEL a, FEL b)
#else
FEL zsub(a, b)
FEL a, b;
#endif

{	register FEL x, bb;

	CHECKFEL(a);
	CHECKFEL(b);

	if (b == F_ZERO) return a;
	if (b == a) return F_ZERO;
	if ((bb = b + minusone) >= Q1) bb -= Q1;	/* bb = -b */
	if (a == F_ZERO) return bb;
	if (a >= bb)
	{	x = inc[a-bb];
		if ((x += bb) >= Q1) x -= Q1;
		return x;
	}
	x = inc[bb-a];
	if ((x += a) >= Q1) x -= Q1;
	return x;
}


#if defined (__STDC__)
FEL zmul(FEL a, FEL b)
#else
FEL zmul(a, b)
FEL a, b;
#endif

{	register FEL c;

	CHECKFEL(a);
	CHECKFEL(b);

	if (b == F_ZERO || a == F_ZERO) return F_ZERO;
	if ((c = a + b) >= Q1) c -= Q1;
	return (c);
}


#if defined (__STDC__)
FEL zdiv(FEL a, FEL b)
#else
FEL zdiv(a, b)
FEL a, b;
#endif

{

	CHECKFEL(a);
	CHECKFEL(b);
	if (b == F_ZERO) FATAL("Division by zero");
	if (a == F_ZERO) return F_ZERO;
	if (a >= b)
		return (a-b);
	else
		return ((a+Q1)-b);
}


/* ------------------------------------------------------------------
   zembed() - Embed a subfield element.

   return value: Embedded field element.
   ------------------------------------------------------------------ */

#if defined (__STDC__)
FEL zrestrict(FEL a, long subfield)
#else
FEL zrestrict(a, subfield)
FEL a;
long subfield;
#endif

{	int i;
	long l;

	/* If a is zero, return zero.
	   -------------------------- */
	if (a == F_ZERO) return F_ZERO;

	/* If a is non-zero, find out the degree of the subfield
	   relative to the current field, and divide by this
	   number.
	   ------------------------------------------------------ */
	for (l = subfield, i = 1; l < Q; l *= subfield, ++i);
	if (l != Q)
		FATAL("ILLEGAL SUBFIELD");
	return (a / i);
}


/* ------------------------------------------------------------------
   zembed() - Embed a subfield element.

   return value: Embedded field element.
   ------------------------------------------------------------------ */

#if defined (__STDC__)
FEL zembed(FEL a, long subfield)
#else
FEL zembed(a, subfield)
FEL a;
long subfield;
#endif

{	int i;
	long l;

	/* If a is zero, return zero.
	   -------------------------- */
	if (a == F_ZERO) return F_ZERO;

	/* If a is non-zero, find out the degree of the subfield
	   relative to the current field, and multiply by this
	   number.
	   ------------------------------------------------------ */
	for (l = subfield, i = 1; l < Q; l *= subfield, ++i);
	if (l != Q)
		FATAL("ILLEGAL SUBFIELD");
	return (a * i);
}


/* ------------------------------------------------------------------
   zsize()  -  Size of a matrix or permutation(s)

   return value: size in bytes
   ------------------------------------------------------------------ */

size_t zsize(nrows)
long nrows;

{	return ((size_t)(rowsize * nrows));
}



/* ------------------------------------------------------------------
   zalloc()  -  Allocate memory and initialize
   ------------------------------------------------------------------ */

PTR zalloc(nrows)
long nrows;			/* number of rows */

{	PTR p;

	if ((p = (PTR) malloc((size_t)(rowsize * nrows))) == NULL)
	{	FATAL("NOT ENOUGH MEMORY");
	}
	memset(p,0xFF,(size_t)(rowsize*nrows));
	return p;
}

/* ------------------------------------------------------------------
   zadvance()  -  advance pointer

   return value: none
   ------------------------------------------------------------------ */

void zadvance(ptr, nrows)
PTR *ptr;
long nrows;

{	*(long **)ptr += nrows * IPR;
}



/* ------------------------------------------------------------------
   zinsert()  -  Insert mark into row
   zextract()  -  Extract mark from row
   ------------------------------------------------------------------ */

#if defined (__STDC__)
void zinsert(PTR row, long col, FEL mark)
#else
void zinsert(row, col, mark)
PTR row;
long col;
FEL mark;
#endif

{	/*CHECKCOL(col);*/
	CHECKFEL(mark);
	row[col-1] = mark;
}

FEL zextract(row, col)
PTR row;
long col;

{	/*CHECKCOL(col);*/

	return (row[col-1]);
}


/* ------------------------------------------------------------------
   zfindpiv()  -  Find first non-zero mark in row

   return value: column (0 if row is zero)
   ------------------------------------------------------------------ */

long zfindpiv(row, mark)
PTR row;
FEL *mark;

{	register long i;
	register PTR p = row;

	for (i = 1; i <= NOC; ++i, ++p)
	{	if (*p != F_ZERO)
		{	*mark = *p;
			return (i);
		}
	}
	return 0;
}


/* ------------------------------------------------------------------
   zmoverow()  -  Move row

   return value: none
   ------------------------------------------------------------------ */

void zmoverow(dest, src)
PTR dest, src;

{	register long *d = (long *) dest, *s = (long *) src;
	register long i = IPR;

	while (i--)
		*d++ = *s++;
}

/* ------------------------------------------------------------------
   zswaprow()  -  Swap two rows

   return value: none
   ------------------------------------------------------------------ */

void zswaprow(dest, src)
PTR dest, src;

{	register long *d = (long *) dest, *s = (long *) src;
	register long i = IPR;
	register long x;

	while (i--)
	{	x = *d;
		*d++ = *s;
		*s++ = x;
	}
}


/* ------------------------------------------------------------------
   zcmprow()  -  compare rows

   return value:
	0	if rows are equal
	1 or -1	if rows are different
   ------------------------------------------------------------------ */

int zcmprow(dest, src)
PTR dest, src;

{	return memcmp(dest,src,(size_t)NOC * sizeof(FEL));
}


/* -----------------------------------------------------------------------
   zaddrow()  -  Add two rows (row1 += row2)

   return value: none
   ----------------------------------------------------------------------- */

void zaddrow(row1, row2)
PTR row1, row2;

{	register long i;
	register FEL *p1 = row1;
	register FEL *p2 = row2;

	for (i = NOC; i != 0; --i)
	{	*p1++ = zadd(*p1,*p2++);
	}
}


/* ------------------------------------------------------------------
   zmulrow()  -  multiply row by field element

   return value: none
   ------------------------------------------------------------------ */

#if defined (__STDC__)
void zmulrow(PTR row, FEL mark)
#else
void zmulrow(row, mark)
PTR row;
FEL mark;
#endif

{	register FEL *m;
	register long i;

	CHECKFEL(mark);
	if (mark == F_ZERO)
	{	m = row;
		for (i = NOC; i != 0; --i)
		{	*m++ = F_ZERO;
		}
	}
	else
	{	m = row;
		for (i = NOC; i != 0; --i)
		{	if (*m != F_ZERO)
			{	if ((*m += mark) >= Q1)
					*m -= Q1;
			}
			++m;
		}
	}
}


/* ------------------------------------------------------------------
   zaddmulrow() - row1 += f * row2

   return value: none
   ------------------------------------------------------------------ */

#if defined (__STDC__)
void zaddmulrow(PTR row1, PTR row2, FEL f)
#else
void zaddmulrow(row1, row2, f)
PTR row1, row2;
FEL f;
#endif

{	register long i;
	register FEL *p1, *p2;

	CHECKFEL(f);
	if (f == F_ZERO) return;
	if (f == F_ONE)
	{	zaddrow(row1,row2);
		return;
	}
	p1 = row1;
	p2 = row2;
	for (i = NOC; i != 0; --i)
	{	*p1 = zadd(*p1,zmul(*p2,f));
		++p1;
		++p2;
	}
}



/* ------------------------------------------------------------------
   zmaprow()  -  Multiply a row by a matrix

   noc(result) = noc(matrix) = current NOC value
   noc(row) = nor(matrix) = nor parameter

   return value: none
   ------------------------------------------------------------------ */

void zmaprow(row, matrix, nor, result)
PTR row, matrix, result;
long nor;			/* number of rows in matrix */

{	register long i;
	register FEL f;
	PTR m = matrix;

	zmulrow(result,F_ZERO);
	for (i = 1; i <= nor; ++i)
	{	f = row[i-1];
		if (f != F_ZERO)
		{	register FEL *v = m;
			register FEL *r = result;
			register long k = NOC;
			if (f == F_ONE)
			{	for (; k != 0; --k)
				{	*r = zadd(*r,*v++);
					++r;
				}
			}
			else
			{	for (; k != 0; --k)
				{	*r = zadd(*r,zmul(*v++,f));
					++r;
				}
			}
		}
		m += NOC;		/* next row */
	}
}


/* ------------------------------------------------------------------
   zpermrow()  -  multiply a row by a permutation

	If 'perm' maps i to k, then the ith mark of 'row' is stored
	in the kth position of 'result'.
	'result' must not eqal 'row'!
   ------------------------------------------------------------------ */

void zpermrow(row, perm, result)
PTR row, result;
PTR perm;

{
    register FEL f;
    register long i;
    register long *p = (long *)perm;

    for (i = 1; i <= znoc; ++i)
    {
	f = zextract(row,i);
	zinsert(result,*p++,f);
    }
}



/* ------------------------------------------------------------------
   zitof(), zftoi() - Convert between long int and FEL
   ------------------------------------------------------------------ */

static void polmul(a, b)
unsigned short *a, *b;

{	unsigned short x[2*MAXPWR], f;
	int i, k;

	memset(x,0,sizeof(x));
	for (i = 0; i < (int) N; ++i)
		for (k = 0; k < (int)N; ++k)
			x[i+k] = (unsigned short)
				((x[i+k] + (long)a[i]*b[k]) % P);
	for (i = 2*(int)N-2; i >= (int)N; --i)
	{	if ((f = x[i]) == 0) continue;
		f = (unsigned short) P - f;
		for (k = (int)N; k >= 0; --k)
		{	x[i-(int)N+k] = (unsigned short)
			    ((x[i-(int)N+k] + (long) f * poly[k]) % P);
		}
	}
	memcpy(a,x,(size_t)N * sizeof(short));
}






FEL zitof(l)
long l;

{	register FEL f = F_ZERO;
	register int i;

	l = l % Q;
	if (l < 0) l += Q;
	for (i = (int)N-1; i >= 0; --i)
	{	while (l >= ppwr[i])
		{	f = zadd(f,ppindex[i]);
			l -= ppwr[i];
		}
	}
	return f;
}


#if defined (__STDC__)
long zftoi(FEL f)
#else
long zftoi(f)
FEL f;
#endif

{	unsigned short i;
	long l, m;
	unsigned short a[MAXPWR], b[MAXPWR];

	if (f == F_ZERO) return 0;
	if (N == 1)
	{	m = Gen;
		l = 1;
		for (i = 1; f != 0; i <<= 1)
		{	if ((f & i) != 0)
			{	l = (l * m) % P;
				f &= ~i;
			}
			m = (m * m) % P;
		}
	}
	else
	{	memset(a,0,sizeof(a));
		a[0] = 1;			/* a(x) = 1 */
		memset(b,0,sizeof(b));
		b[1] = 1;			/* b(x) = x */
		for (i = 1; f != 0; i <<= 1)
		{	if ((f & i) != 0)
			{	polmul(a,b);
				f &= ~i;
			}
			polmul(b,b);
		}
		l = 0;
		m = 1;
		for (i = 0; i < (short) N; ++i)
		{	l += a[i] * m;
			m *= P;
		}
	}
	return l;
}


/* ==================================================================
   Common functions
   ================================================================== */

