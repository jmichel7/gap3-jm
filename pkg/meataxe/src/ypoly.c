/* ========================== C MeatAxe =============================
   ypoly.c - Polynomials over finite fields

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: ypoly.c,v 1.2 1997/09/11 15:43:32 gap Exp $
 *
 * $Log: ypoly.c,v $
 * Revision 1.2  1997/09/11 15:43:32  gap
 * New version 2.2.3. AH
 *
 * Revision 2.11  1995/02/08  10:14:52  mringe
 * _PL entfernt.
 *
 * Revision 2.11  1995/02/08  10:14:52  mringe
 * _PL entfernt.
 *
 * Revision 2.10  1994/08/22  19:30:33  mringe
 * polread()/polwrite(): Benutzt standard fileformat.
 *
 * Revision 2.9  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.8  1994/05/04  11:41:37  mringe
 * hmmm....
 *
 * Revision 2.7  1994/04/12  13:35:32  mringe
 * BUG behoben: polgcd(x,0) liefert jetzt x statt 1 zurueck.
 *
 * Revision 2.6  1994/03/14  12:59:28  mringe
 * Geaendert: vec<-->pol, polprint().
 *
 * Revision 2.5  1994/03/13  13:28:13  mringe
 * Benutze zreadlong()/zwritelong()
 *
 * Revision 2.4  1994/02/12  04:10:13  mringe
 * UMFANGREICHE AENDERUNGEN AN VIELEN DATENTYPEN.
 *
 * Revision 2.3  1994/02/11  10:23:33  mringe
 * gcd(p,0) = 0.
 *
 * Revision 2.2  1994/02/11  09:49:32  mringe
 * Diverse Korrekturen.
 *
 * Revision 2.1  1994/02/04  16:40:25  mringe
 * polderive(), polpack(),...
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
 * Revision 1.6  1993/10/05  18:58:31  mringe
 * Neue fehlermeldungen.
 *
 * Revision 1.5  1993/10/05  14:20:58  mringe
 * Bug in poldup() behoben.
 *
 * Revision 1.4  1993/10/02  16:06:10  mringe
 * Neue Polynomfunktionen.
 *
 * Revision 1.3  1993/09/29  10:48:40  mringe
 * Neu: Polynome.
 *
 * Revision 1.2  1993/09/21  07:17:44  mringe
 * *** empty log message ***
 *
 * Revision 1.1  1993/08/31  12:58:35  mringe
 * Initial revision
 *
 */


#include <string.h>
#include <stdlib.h>
#include "meataxe.h"


#define INITIAL_SIZE 10

/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */


/* ------------------------------------------------------------------
   poldesc() - Describe a matrix
   ------------------------------------------------------------------ */
    
char *poldesc(p)
poly_t *p;

{
    static char buf[50];
    sprintf(buf,"polynomial of degree %ld over GF(%ld)",
	p->deg, p->fl);
    return buf;
}


/* ------------------------------------------------------------------
   polalloc() - Allocate a polynomial, initialize to X^(degree)
   polfree() - Free a polynomial
   poldup() - Duplicate a polynomial
   ------------------------------------------------------------------ */

poly_t *polalloc(fl,degree)
long fl;
long degree;

{
    poly_t *x;
    size_t s;
    int i;

    if (degree < 0) degree = -1;
    if ((s = degree + 1) < 1) s = 1;
    zsetfield(fl);
    x = (poly_t *) malloc(sizeof(poly_t));
    if (x == NULL) MTXFAIL(ERR_NOMEM,NULL)
    x->id = T_POLY;
    x->fl = fl;
    x->deg = degree;
    x->size = s;
    x->buf = (FEL *) malloc(s * sizeof(FEL));
    if (x->buf == NULL)
    {
	free(x);
        MTXFAIL(ERR_NOMEM,NULL)
    }
    for (i = 0; i < s-1; ++i) x->buf[i] = F_ZERO;
    x->buf[s-1] = F_ONE;
    return x;
}


void polfree(x)
poly_t *x;

{
    free(x->buf);
    free(x);
}


poly_t *poldup(x)
poly_t *x;

{
    poly_t *y;


    if ((y = malloc(sizeof(poly_t))) == NULL)
	MTXFAIL(ERR_NOMEM,NULL);
    y->id = T_POLY;
    y->buf = (FEL *) malloc((x->deg+1) * sizeof(FEL));
    if (y->buf == NULL)
    {
	free(y);
        MTXFAIL(ERR_NOMEM,NULL)
    }
    memcpy(y->buf,x->buf,x->deg+1);
    y->fl = x->fl;
    y->deg = x->deg;
    y->size = x->deg + 1;
    return y;
}


/* ------------------------------------------------------------------
   resize() - Extend buffer.
   ------------------------------------------------------------------ */

static int resize(p, newdeg)
poly_t *p;
size_t newdeg;

{
    int i;
    FEL *x;

    if (p->deg < newdeg)
    {
	if (p->size < newdeg+1)	/* Allocate new buffer */
    	{
	    x = (FEL *) malloc((newdeg+1) * sizeof(FEL));
	    if (x == NULL) return 1;
	    memcpy(x,p->buf,(p->deg+1)*sizeof(FEL));
	    free(p->buf);
	    p->buf = x;
	    p->size = newdeg + 1;
    	}
    	for (i = p->deg + 1; i <= newdeg; ++i)
	    p->buf[i] = F_ZERO;
        p->deg = newdeg;
    }
    return 0;
}


/* ------------------------------------------------------------------
   normalize()
   ------------------------------------------------------------------ */

static void normalize(p)
poly_t *p;

{
    int i = p->deg;

    while (i >= 0 && p->buf[i] == F_ZERO) --i;
    p->deg = i;
}



/* ------------------------------------------------------------------
   poladd() - Add polynomials: dest += src
   ------------------------------------------------------------------ */

poly_t *poladd(dest, src)
poly_t *dest, *src;

{
    FEL *s, *d;
    int i;

    if (dest->fl != src->fl) MTXFAIL(ERR_INCOMPAT,NULL);
    if ((i = src->deg) == -1) return dest;	/* src = 0 */
    zsetfield(src->fl);
    if (resize(dest,i)) MTXFAIL(ERR_NOMEM,NULL)
    s = src->buf;
    d = dest->buf;
    for (; i >= 0; --i)
    {
       *d = zadd(*d,*s++);
       ++d;
    }
    normalize(dest);
    return dest;
}


/* ------------------------------------------------------------------
   polmul() - Multiply polynomials: dest *= src
   ------------------------------------------------------------------ */

poly_t *polmul(dest, src)
poly_t *dest, *src;

{
    FEL *x, *y, *d = dest->buf, *s = src->buf;
    int di, si;
    size_t xdeg = src->deg + dest->deg;


    if (dest->fl != src->fl) MTXFAIL(ERR_INCOMPAT,NULL);
    zsetfield(src->fl);
    if (dest->deg == -1)	/* dest = 0 */
	return dest;
    if (src->deg == -1)		/* src = 0 */
    {
	dest->deg = -1;
	return dest;
    }
    x = (FEL *) malloc((xdeg+1) * sizeof(FEL));
    if (x == NULL) MTXFAIL(ERR_NOMEM,NULL);
    for (di = xdeg, y = x; di >= 0; --di) *y++ = F_ZERO;

    for (di = 0; di <= dest->deg; ++di)
    	for (si = 0; si <= src->deg; ++si)
	    x[si+di] = zadd(x[si+di],zmul(s[si],d[di]));
    free(dest->buf);
    dest->buf = x;
    dest->deg = xdeg;
    dest->size= xdeg+1;
    return dest;
}


/* ------------------------------------------------------------------
   poldivmod() - Polynomial division. Returns the quotient,
	remainder is stored in a;
   ------------------------------------------------------------------ */

poly_t *poldivmod(a,b)
poly_t *a, *b;

{
    poly_t *q;

    if (a->fl != b->fl) MTXFAIL(ERR_INCOMPAT,NULL);
    zsetfield(a->fl);
    if (b->deg <= -1) MTXFAIL(ERR_DIV0,NULL);
    if (a->deg < b->deg)
    {
	q = polalloc(a->fl,-1);	/* Trivial case: Quotient = 0 */
    }
    else
    {
	FEL lead = b->buf[b->deg];
	int i, k;

	if (lead == F_ZERO) MTXFAIL(ERR_DIV0,NULL);
	q = polalloc(zfl,a->deg - b->deg);
	if (q == NULL) return NULL;
	for (i = a->deg; i >= b->deg; --i)
	{
	    FEL qq = zneg(zdiv(a->buf[i],lead));
	    for (k = 0; k <= b->deg; ++k)
		a->buf[i-k] = zadd(a->buf[i-k],
		    zmul(qq,b->buf[b->deg - k]));
	    q->buf[i-b->deg] = zneg(qq);
	}
        normalize(a);
    }
    return q;
}


/* ------------------------------------------------------------------
   polmod() - Reduce polynomial `a' modulo `b'
   ------------------------------------------------------------------ */

poly_t *polmod(a,b)
poly_t *a, *b;

{
    if (a->fl != b->fl) MTXFAIL(ERR_INCOMPAT,NULL);
    zsetfield(a->fl);
    if (b->deg <= -1) MTXFAIL(ERR_DIV0,NULL);
    if (a->deg >= b->deg)
    {
	FEL lead = b->buf[b->deg];
	int i, k;

	if (lead == F_ZERO) MTXFAIL(ERR_DIV0,NULL);
	for (i = a->deg; i >= b->deg; --i)
	{
	    FEL qq = zneg(zdiv(a->buf[i],lead));
	    for (k = 0; k <= b->deg; ++k)
		a->buf[i-k] = zadd(a->buf[i-k],
		    zmul(qq,b->buf[b->deg - k]));
	}
        normalize(a);
    }
    return a;
}

/* -----------------------------------------------------------------
   polshiftmod() - Multiply p(x) by x^n and reduce mod q(x)
   ----------------------------------------------------------------- */

poly_t *polshiftmod(p,n,q)
poly_t *p;
long n;
poly_t *q;

{
    FEL qlead = q->buf[q->deg];
    FEL *qbuf;
    FEL *pbuf;
    int deg;
    if (p->deg >= q->deg) polmod(p,q);
    resize(p,q->deg-1);
    qbuf = q->buf;
    pbuf = p->buf;
    deg = q->deg;
    while (n > 0)
    {
	int i,k;
	FEL f;

	for (i = deg-1; i >= 0 && pbuf[i] != F_ZERO; --i);
	if (i < 0) break;
	if (deg - i > n) break;
	f = zdiv(pbuf[i],qlead);
	for (k = deg - 1; k >= deg-i; --k)
	    pbuf[k] = zsub(pbuf[k-deg+i],zmul(qbuf[k],f));
	for (; k >= 0; --k)
	    pbuf[k] = zneg(zmul(qbuf[k],f));
	n -= deg - i;
    }
    normalize(p);
    return p;
}


/* -----------------------------------------------------------------
   polprint() - Print a polynomial to stdout
   ----------------------------------------------------------------- */

void polprint(name,p)

char *name;
poly_t *p;

{
    int i,flag = 0;

    if (name != NULL) printf("%s=",name);
    zsetfield(p->fl);
    if (p->deg == -1) 
    {
	printf("0x^0");
    }
    for (i = p->deg; i >= 0; i--)
    {
	if (p->buf[i] != F_ZERO)
	{
	    if (flag) printf("+");
	    if (p->buf[i] != F_ONE || i == 0)
		printf("%ld",zftoi(p->buf[i]));
	    if (i > 1) printf("x^%d",i);
	    else if (i == 1) printf("x");
	    flag=1;
       	}
    }
    if (name != NULL) printf("\n");
}


/* -----------------------------------------------------------------
   polread() - Read a polynomial from a file.
   ----------------------------------------------------------------- */

static long tmpfl = 0;
static long tmpdeg = 0;
static PTR tmpvec = NULL;

static void mktmp(long fl, long deg)

{
    zsetfield(fl);
    if (deg > 0) zsetlen(deg+1);
    if (tmpfl != fl || tmpdeg < deg)
    {
	if (tmpvec != NULL) free(tmpvec);
	tmpvec = zalloc((long)1);
	tmpdeg = deg;
	tmpfl = fl;
    }
}

poly_t *polread(f)
FILE *f;

{
    poly_t *p;
    long hdr[3];
    
    if (zreadlong(f,hdr,3) != 3)
	MTXFAIL(ERR_FILEREAD,NULL);
    if (hdr[0] != -2)
	MTXFAIL(ERR_BADTYPE,NULL);
    mktmp(hdr[1],hdr[2]);
    if ((p = polalloc(hdr[1],hdr[2])) == NULL)
	return NULL;
    if (p->deg > 0)
    {
	int i;
	if (zreadvec(f,tmpvec,1) != 1)
	{
	    polfree(p);
	    MTXFAIL(ERR_FILEREAD,NULL);
	}
	for (i = 0; i <= p->deg; ++i)
	    p->buf[i] = zextract(tmpvec,i+1);
    }
    return p;
}


/* -----------------------------------------------------------------
   polwrite() - Write a polynomial to a file.
   ----------------------------------------------------------------- */

int polwrite(f,p)

FILE *f;
poly_t *p;

{
    long hdr[3];

    mktmp(p->fl,p->deg);
    hdr[0] = -2;
    hdr[1] = p->fl;
    hdr[2] = p->deg;
    if (zwritelong(f,hdr,3) != 3) 
	MTXFAIL(ERR_FILEWRITE,-1);
    if (p->deg >= 0)
    {
	int i;
    	for (i = 0; i <= p->deg; ++i)
	    zinsert(tmpvec,i+1,p->buf[i]);
        if (zwritevec(f,tmpvec,1) != 1)
	    MTXFAIL(ERR_FILEWRITE,-1);
    }
    return 0;
}



/* -----------------------------------------------------------------
   polcmp() - Compare polynomials.
   ----------------------------------------------------------------- */

int polcmp(a,b)
poly_t *a, *b;

{
    int i;

    if (a->fl > b->fl) return 1;
    if (a->fl < b->fl) return -1;

    if (a->deg > b->deg) return 1;
    if (a->deg < b->deg) return -1;

    for (i = a->deg; i >= 0; --i)
    {
	if (a->buf[i] > b->buf[i]) return 1;
	if (a->buf[i] < b->buf[i]) return -1;
    }

    return 0;
}


/* -----------------------------------------------------------------
   polgcd() - Calculate the g.c.d. of two Polynomials.
   ----------------------------------------------------------------- */

poly_t *polgcd(a,b)
poly_t *a, *b;

{
    poly_t *p,*q,*tmp;
    FEL f;

    if (a->fl != b->fl) MTXFAIL(ERR_INCOMPAT,NULL);

    /* Handle special cases
       -------------------- */
    if (a->deg == -1)
    {
	if (b->deg == -1) MTXFAIL(ERR_DIV0,NULL);
	return poldup(b);
    }
    else if (b->deg == -1)
	return poldup(a);

    /* a,b != 0: Find gcd
       ------------------ */
    zsetfield(a->fl);
    if (a->deg < b->deg)
    {
	p=poldup(b);
	q=poldup(a);
    }
    else
    {
	p=poldup(a);
	q=poldup(b);
    }
   
    while (q->deg >= 0)
    {
	polmod(p,q);
	tmp=p;
	p=q;
	q=tmp;
    }
    polfree(q);
    if ((f = p->buf[p->deg]) != F_ONE)
    {
	FEL *buf = p->buf;
	int i = p->deg;
	while (i-- > 0)
	{
	    if (*buf != F_ZERO) *buf = zdiv(*buf,f);
	    ++buf;
	}
	*buf = F_ONE;
    }
    return p;
}


/* -----------------------------------------------------------------
   pol2vec() - Convert polynomial to packed vector. Sets field and
	row size.
   vec2pol() - Convert packed vector to polynomial. Use zfl, znoc.

   These functions are for char. pols. Normalization is done
   automatically.
   ----------------------------------------------------------------- */

int pol2vec(p,vec)
poly_t *p;
PTR vec;

{
    long i;
    FEL f;

    if (p->deg < 1) MTXFAIL(ERR_BADARG,-1);
    zsetfield(p->fl);
    zsetlen(p->deg);
    zmulrow(vec,F_ZERO);
    for (i = 0; i < znoc; ++i)
	zinsert(vec,i+1,p->buf[i]);
    if ((f = p->buf[p->deg]) != F_ONE)
	zmulrow(vec,zinv(f));
    return 0;
}


int vec2pol(vec,pol)
PTR vec;
poly_t *pol;

{
    long i;

    if (resize(pol,znoc))
	MTXFAIL(ERR_NOMEM,-1);
    pol->fl = zfl;
    pol->deg = znoc;
    for (i = 1; i <= znoc; ++i)
	pol->buf[i-1] = zneg(zextract(vec,i));
    pol->buf[znoc] = F_ONE;
    return 0;
}

/* -----------------------------------------------------------------
   polpack() - Pack polynomial into row. Sets field and row length.
   ----------------------------------------------------------------- */

int polpack(p,vec)
poly_t *p;
PTR vec;

{
    long i;

    zsetfield(p->fl);
    zsetlen(p->deg+1);
    zmulrow(vec,F_ZERO);
    for (i = 0; i <= p->deg; ++i)
	zinsert(vec,i+1,p->buf[i]);
    return 0;
}


/* -----------------------------------------------------------------
   polderive() - Derive a polynomial. Returns pol, makes no copy of
       the original polynomial.
   ----------------------------------------------------------------- */

poly_t *polderive(pol)
poly_t *pol;

{
    int i, maxdeg = -1;
    FEL *buf = pol->buf;
    FEL f = F_ZERO;

    zsetfield(pol->fl);
    for (i = 0; i < pol->deg; ++i)
    {
	f = zadd(f,F_ONE);
	buf[i] = zmul(buf[i+1],f);
	if (buf[i] != F_ZERO) maxdeg = i;
    }
    pol->deg = maxdeg;
    return pol;
}
      

