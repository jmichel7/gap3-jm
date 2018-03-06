/* ========================== C MeatAxe =============================
   ymatrix.c -  Matrix operations

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */

/* $Id: ymatrix.c,v 1.2 1997/09/11 15:43:29 gap Exp $
 *
 * $Log: ymatrix.c,v $
 * Revision 1.2  1997/09/11 15:43:29  gap
 * New version 2.2.3. AH
 *
 * Revision 2.16  1995/07/05  10:33:53  mringe
 * BUG in matadd() behoben.
 *
 * Revision 2.15  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.14  1994/11/25  13:32:13  mringe
 * ANSI-C
 *
 * Revision 2.13  1994/09/29  14:31:05  mringe
 * Bug in matorder() behoben -- Thanks to Steve Linton.
 *
 * Revision 2.12  1994/08/22  08:24:01  mringe
 * matpower(x,0) gibt jetzt die Einheitsmatrix zurueck
 *
 * Revision 2.11  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.10  1994/07/23  16:48:52  mringe
 * Nullraum- und Echelonfunktionen nach gauss.c ausgelagert.
 *
 * Revision 2.9  1994/06/25  13:34:18  mringe
 * Benutze zextractcol().
 *
 * Revision 2.9  1994/06/25  13:34:18  mringe
 * Benutze zextractcol().
 *
 * Revision 2.8  1994/06/16  14:19:10  mringe
 * Profiling.
 *
 * Revision 2.7  1994/05/20  07:18:26  mringe
 * Benutze Smalloc() und Srealloc().
 *
 * Revision 2.6  1994/05/18  05:16:48  mringe
 * matquot() entfernt.
 *
 * Revision 2.5  1994/05/11  08:43:03  mringe
 * Neue Version von zquot().
 *
 * Revision 2.4  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.3  1994/02/12  04:10:13  mringe
 * UMFANGREICHE AENDERUNGEN AN VIELEN DATENTYPEN.
 *
 * Revision 2.2  1993/12/07  14:34:24  mringe
 * Teste bei realloc() auf size==0.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.17  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.16  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.15  1993/10/05  18:57:28  mringe
 * Neue Fehlermeldungen.
 *
 * Revision 1.14  1993/10/05  11:54:02  mringe
 * Benutze zreadvec()/zwritevec()
 *
 * Revision 1.13  1993/10/02  16:23:02  mringe
 * matread() und matwrite() in matload() bzw. matsave() umbenannt.
 *
 * Revision 1.12  1993/10/02  16:09:23  mringe
 * prmat() in matprint umbenannt.
 *
 * Revision 1.11  1993/08/31  12:57:03  mringe
 * ERR_DIV0 ersetzt ERR_SINGULAR
 *
 * Revision 1.10  1993/08/12  16:14:23  mringe
 * Include <string.h>.
 *
 * Revision 1.9  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 */


#include <string.h>
#include <stdlib.h>

#include "meataxe.h"





#define FOREACHROW(m,i,x)\
    for ((i)=1,(x)=(m)->d;(i)<=(m)->nor;++(i),zadvance(&(x),(long)1))



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static long *piv = NULL;		/* Pivot table */
static long maxnoc = -1;		/* Table size */


/* ------------------------------------------------------------------
   matdesc() - Describe a matrix
   ------------------------------------------------------------------ */
    
char *matdesc(matrix_t *m)

{
    static char buf[50];
    sprintf(buf,"%ld x %ld matrix over GF(%ld)",
	m->nor, m->noc, m->fl);
    return buf;
}

/* ------------------------------------------------------------------
   matprint() - Print a matrix to stdout
   ------------------------------------------------------------------ */

void matprint(char *name, matrix_t *m)

{
    PTR x;
    long i, k;

    zsetfield(m->fl);
    zsetlen(m->noc);
    x = m->d;
    if (name != NULL) printf("%s=\n",name);
    for (i = 1; i <= m->nor; ++i)
    {
	for (k = 1; k <= m->noc; ++k)
	    printf("%d",zextract(x,k));
	printf("\n");
	zadvance(&x,(long)1);
    }
}



/* ------------------------------------------------------------------
   pivalloc() - Allocate pivot table. All functions in this module
	use one global pivot table (piv). pivalloc() checks first
	if the currently allocated table is large enough.
   ------------------------------------------------------------------ */

static void pivalloc(long noc)

{
    if (noc <= maxnoc)
	return;
    maxnoc = noc;
    if (piv != NULL)
	free(piv);
    piv = NALLOC(long,noc+1);
    if (piv == NULL)
	FATAL("pivalloc(): out of memory");
}


/* ------------------------------------------------------------------
   matalloc() - Allocate a matrix
   matfree() - Free a matrix
   ------------------------------------------------------------------ */

matrix_t *matalloc(long fl, long nor, long noc)

{
    matrix_t *m;

    if (fl == 0) fl = zfl;	/* Default value */
    zsetfield(fl);
    zsetlen(noc);
    m = ALLOC(matrix_t);
    if (m == NULL)
	MTXFAIL(ERR_NOMEM,NULL);
    m->id = T_MATRIX;
    m->d = zalloc(nor);
    if (m->d == NULL)
    {
	free(m);
	MTXFAIL(ERR_NOMEM,NULL);
    }
    m->fl = fl;
    m->nor = nor;
    m->noc = noc;
    return m;
}


void matfree(matrix_t *m)

{
    if (m->d != NULL)
	free(m->d);
    free(m);
}




/* ------------------------------------------------------------------
   matid() - Identity matrix
   ------------------------------------------------------------------ */

matrix_t *matid(long fl, long nor)

{
    matrix_t *m;
    PTR x;
    long i;

    m = matalloc(fl,nor,nor);
    if (m == NULL) return NULL;
    FOREACHROW(m,i,x)
    {
	zmulrow(x,F_ZERO);
	zinsert(x,i,F_ONE);
    }
    return m;
}




/* ------------------------------------------------------------------
   matdup() - Duplicate a matrix
   matmove() - Move a matrix
   ------------------------------------------------------------------ */

matrix_t *matdup(matrix_t *src)

{
    matrix_t *m;

    m = matalloc(src->fl,src->nor,src->noc);
    if (m == NULL) return NULL;
    memcpy(m->d,src->d,zsize(src->nor));
    return m;
}


int matmove(matrix_t *dest, matrix_t *src)

{
    if (dest->fl != src->fl || dest->nor != src->nor ||
	dest->noc != src->noc)
    {
	MTXFAIL(ERR_INCOMPAT,-1);
    }
    zsetfield(src->fl);
    zsetlen(src->noc);
    memcpy(dest->d,src->d,zsize(src->nor));
    return 0;
}



/* ------------------------------------------------------------------
   matextract() - Extract <n> rows out of a matrix, starting with
   <first>.
   ------------------------------------------------------------------ */

matrix_t *matextract(matrix_t *src,long first,long n)

{
    matrix_t *m;
    PTR p;

    zsetfield(src->fl);
    zsetlen(src->noc);
    if (first < 1 || first+n-1 > src->nor || n < 0)
    {
	MTXFAIL(ERR_RANGE,NULL);
    }
    m = matalloc(src->fl,n,src->noc);
    if (m == NULL) return NULL;
    p = src->d;
    zadvance(&p,first-1);
    memcpy(m->d,p,zsize(n));
    return m;
}


/* ------------------------------------------------------------------
   matread(), matwrite() - Read and write matrices
   ------------------------------------------------------------------ */

matrix_t *matread(FILE *f)

{
    matrix_t *m;
    long hdr[3];

    if (zreadlong(f,hdr,3) != 3) MTXFAIL(ERR_FILEREAD,NULL);
    if (hdr[0] < 2) MTXFAIL(ERR_BADTYPE,NULL);
    m = matalloc(hdr[0],hdr[1],hdr[2]);
    if (m == NULL) return NULL;
    if (zreadvec(f,m->d,(size_t)m->nor) != m->nor)
    {
	matfree(m);
	return NULL;
    }
    return m;
}

int matwrite(FILE *f, matrix_t *mat)

{
    long hdr[3];

    hdr[0] = mat->fl;
    hdr[1] = mat->nor;
    hdr[2] = mat->noc;
    if (zwritelong(f,hdr,3) != 3) MTXFAIL(ERR_FILEWRITE,-1);
    zsetfield(mat->fl);
    zsetlen(mat->noc);
    if (zwritevec(f,mat->d,(size_t)mat->nor) != mat->nor)
	return -1;
    return 0;
}



/* ------------------------------------------------------------------
   matload(), matsave() - Read and write matrices
   ------------------------------------------------------------------ */

matrix_t *matload(char *fn)

{
    FILE *f;
    matrix_t *m;

    if ((f = SFOpen(fn,FM_READ)) == NULL)
    {
	perror(fn);
	errexit(ERR_FILEREAD,fn);
    }
    m = matread(f);
    fclose(f);
    if (m == NULL)
    {
	mtxerror(fn);
	errexit(ERR_FILEREAD,fn);
    }
    return m;
}


int matsave(matrix_t *mat, char *fn)

{
    FILE *f;
    int i;

    if ((f = SFOpen(fn,FM_CREATE)) == NULL)
    {
	perror(fn);
	errexit(ERR_FILEWRITE,fn);
    }
    i = matwrite(f,mat);
    fclose(f);
    if (i != 0)
    {
	mtxerror(fn);
	errexit(ERR_FILEWRITE,fn);
    }
    return i;
}



/* ------------------------------------------------------------------
   matadd(), matmul() - Basic matrix arithmetic
   ------------------------------------------------------------------ */

matrix_t *matadd(matrix_t *dest, matrix_t *src)

{
    PTR s, d;
    long i;

    if (src->fl != dest->fl || src->noc != dest->noc ||
	src->nor != dest->nor)
	MTXFAIL(ERR_INCOMPAT,NULL);
    zsetfield(src->fl);
    zsetlen(src->noc);
    d = dest->d;
    s = src->d;
    for (i = src->nor; i != 0; --i)
    {
	zaddrow(d,s);
	zadvance(&d,(long)1);
	zadvance(&s,(long)1);
    }
    return dest;
}



matrix_t *matmul(matrix_t *dest, matrix_t *src)

{
    PTR x, tmp, result;
    long i;

    if (src->fl != dest->fl || src->nor != dest->noc)
    {
	MTXFAIL(ERR_INCOMPAT,NULL);
    }
    zsetfield(src->fl);
    zsetlen(src->noc);
    result = tmp = zalloc(dest->nor);
    if (result == NULL)
    {
	MTXFAIL(ERR_NOMEM,NULL);
    }
    x = dest->d;
    for (i = dest->nor; i != 0; --i)
    {
	zsetlen(src->noc);
	zmaprow(x,src->d,src->nor,tmp);
	zadvance(&tmp,(long)1);
	zsetlen(dest->noc);
	zadvance(&x,(long)1);
    }
    free(dest->d);
    dest->d = result;
    dest->noc = src->noc;
    return dest;
}


/* ------------------------------------------------------------------
   mattr() - Transpose a matrix
   ------------------------------------------------------------------ */

matrix_t *mattr(matrix_t *src)

{
    PTR s, d;
    long i;
    matrix_t *dest;

    dest = matalloc(src->fl,src->noc,src->nor);
    if (dest == NULL) return NULL;
    d = dest->d;
    s = src->d;
    for (i = 1; i <= src->noc; ++i)
    {
	zsetlen(dest->noc);
	zmulrow(d,F_ZERO);
	zsetlen(src->noc);
	zextractcol(s,src->nor,i,d);
	zsetlen(dest->noc);
	zadvance(&d,(long)1);
    }
    return dest;
}


/* ------------------------------------------------------------------
   matinv() - Matrix inversion
   ------------------------------------------------------------------ */

matrix_t *matinv(matrix_t *src)

{
    PTR tmp = NULL;	/* Workspace */
    matrix_t *dest;

    if (src->nor != src->noc)
    {
	MTXFAIL(ERR_NOTSQUARE,NULL);
    }
    dest = matid(src->fl,src->nor);
    if (dest == NULL) return NULL;

    /* Copy matrix into workspace
       -------------------------- */
    tmp = zalloc(src->nor);
    memcpy(tmp,src->d,zsize(src->nor));

    /* Inversion
       --------- */
    if (zmatinv(tmp,dest->d) != 0) 
    {
	matfree(dest);
	return NULL;
    }
    return dest;
}


/* ------------------------------------------------------------------
   chkechelon() - Check if a matrix is in full echelon form
   ------------------------------------------------------------------ */

int chkechelon(matrix_t *mat)

{
    PTR x, y;
    long i, k, piv;
    FEL f;

    zsetfield(mat->fl);
    zsetlen(mat->noc);
    x = mat->d;
    for (i = 1; i <= mat->nor; ++i)
    {	piv = zfindpiv(x,&f);
    	if (f != F_ONE) return 1;
    	y = mat->d;
    	for (k = 1; k <= mat->nor; ++k, zadvance(&y,(long)1))
    	{
	    if (k == i) continue;
    	    if (zextract(y,piv) != F_ZERO) return 1;
    	}
    	zadvance(&x,(long)1);
    }
    return 0;
}


/* ------------------------------------------------------------------
   matpower() - Calculate the <n>-th power of <mat> using the
	binary method.
   ------------------------------------------------------------------ */

static void matpwr_(long n, PTR inp, PTR out, PTR tmp2)

{
    PTR x, y;
    long i;
    int first = 1;

    while (n > 0)
    {	if (n % 2 == 1)
	{
	    if (first)
	    {
		memcpy(out,inp,zsize(znoc));
		first = 0;
	    }
	    else
	    {
		x = out;
		for (i = 1; i <= znoc; ++i)
		{
		    zmaprow(x,inp,znoc,tmp2);
		    zmoverow(x,tmp2);
		    zadvance(&x,(long)1);
		}
	    }
	}
	x = inp;
	y = tmp2;
	for (i = 1; i <= znoc; ++i)
	{
	    zmaprow(x,inp,znoc,y);
	    zadvance(&x,(long)1);
	    zadvance(&y,(long)1);
	}
	memcpy(inp,tmp2,zsize(znoc));
	n /= 2;
    }
}


matrix_t *matpower(matrix_t *mat, long n)

{
    matrix_t *result;
    PTR tmp, tmp2;
    if (mat->nor != mat->noc)
    {
	MTXFAIL(ERR_NOTSQUARE,NULL);
    }
    if (n == 0) return matid(mat->fl,mat->nor);
    zsetfield(mat->fl);
    zsetlen(mat->noc);
    tmp = zalloc(znoc);
    memcpy(tmp,mat->d,zsize(znoc));
    tmp2 = zalloc(znoc);
    result = matalloc(mat->fl,mat->nor,mat->noc);
    if (result == NULL)
    {
	free(tmp);
	free(tmp2);
	return NULL;
    }
    matpwr_(n,tmp,result->d,tmp2);
    free(tmp);
    free(tmp2);
    return result;
}


/* ------------------------------------------------------------------
   matorder() - Order of a matrix. Returns -1 on failure.
   ------------------------------------------------------------------ */

long matorder(matrix_t *mat)

{
    PTR m1, v1, v2, v3;
    PTR basis, bend, bptr;
    char *done;
    long dim;
    long j1, i;
    FEL f1;
    long tord;
    long ord;
    int flag;

    if (mat->nor != mat->noc)
    {
	MTXFAIL(ERR_NOTSQUARE,-1);
    }
    zsetfield(mat->fl);
    zsetlen(mat->noc);
    m1 = zalloc(mat->nor);
    memcpy(m1,mat->d,zsize(mat->nor));
    bend = basis = zalloc(mat->nor+1);
    pivalloc(mat->nor+1);
    done = NALLOC(char,mat->nor+1);
    memset(done,0,(size_t)((mat->nor+1)));
    v1 = zalloc((long)1);
    v2 = zalloc((long)1);
    v3 = zalloc((long)1);
    tord = ord = 1;
    dim = 0;
    j1 = 1;
    while (dim < mat->nor && tord <= 1000 && ord <= 1000000)
    {
	/* Get next start vector
	   --------------------- */
	for (j1 = 1; j1 <= mat->nor && done[j1]; ++j1);
	if (j1 > mat->nor) break;	/* Done! */
	zmulrow(v1,F_ZERO);
	zinsert(v1,j1,F_ONE);

	/* Calculate order on cyclic subspace
	   ---------------------------------- */
	tord = 0;
	flag = 1;
	zmoverow(v3,v1);
	do
	{   zmoverow(v2,v3);
	    if (flag)
	    {   zmoverow(bend,v3);
	        bptr = basis;
	        for (i = 1; i <= dim; ++i)
	        {   f1 = zextract(bend,piv[i]);
		    if (f1 != 0)
			zaddmulrow(bend,bptr,zneg(
			zdiv(f1,zextract(bptr,piv[i]))));
		    zadvance(&bptr,(long)1);
	        }
		if ((piv[dim+1] = zfindpiv(bend,&f1)) != 0)
		{   done[piv[++dim]] = 1;
		    zadvance(&bend,(long)1);
		}
		else
		    flag = 0;
	    }
	    zmaprow(v2,m1,mat->nor,v3);
	    ++tord;
	}
	while (tord <= 1000 && zcmprow(v3,v1) != 0);

	/* Order = lcm(orders on cyclic subspaces)
	   --------------------------------------- */
	ord = lcm(ord,tord);
    }
    free(done);
    free(v1);
    free(v2);
    free(v3);
    free(m1);
    free(basis);
    if (tord > 1000 || ord > 1000000)
	return -1;
    return ord;
}


