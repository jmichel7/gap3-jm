/* ========================== C MeatAxe =============================
   yperm.c - Operations with permutations

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: yperm.c,v 1.2 1997/09/11 15:43:30 gap Exp $
 *
 * $Log: yperm.c,v $
 * Revision 1.2  1997/09/11 15:43:30  gap
 * New version 2.2.3. AH
 *
 * Revision 2.9  1995/06/22  07:53:46  mringe
 * permorder() optimiert - jetzt deutlich schneller
 *
 * Revision 2.8  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.7  1994/07/05  16:33:54  mringe
 * *** empty log message ***
 *
 * Revision 2.6  1994/06/29  07:50:56  mringe
 * Neu: permprint()
 *
 * Revision 2.5  1994/06/27  10:25:37  mringe
 * Bug in permdup(): zsize() nicht mehr fuer Perm. geeignet!
 *
 * Revision 2.4  1994/03/12  13:38:48  mringe
 * Benutze zreadlong() und zwritelong().
 *
 * Revision 2.3  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.3  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.2  1994/02/12  04:10:13  mringe
 * UMFANGREICHE AENDERUNGEN AN VIELEN DATENTYPEN.
 *
 * Revision 2.1  1993/10/21  21:57:35  mringe
 * Permutationen.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.10  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.9  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.8  1993/10/05  19:04:39  mringe
 * *** empty log message ***
 *
 * Revision 1.7  1993/10/05  19:02:08  mringe
 * yerror eliminiert.
 *
 * Revision 1.6  1993/10/05  18:57:28  mringe
 * Neue Fehlermeldungen.
 *
 * Revision 1.5  1993/10/05  17:33:26  mringe
 * permsave() und permload().
 *
 * Revision 1.4  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.3  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.2  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.1  1993/02/10  19:40:54  mringe
 * Initial revision
 *
 * Revision 1.5  1992/07/22  07:10:30  mringe
 * Changed 'global.h' to 'lattice.h'
 *
 * Revision 1.4  1992/07/15  09:25:55  mringe
 * Some minor changes.
 *
 * Revision 1.3  1992/07/10  15:21:55  mringe
 * Typdeklarationen aus Arugumentlisten entfernt.
 *
 * Revision 1.2  1992/05/29  07:29:53  mringe
 * Neu: permpower()
 *
 * Revision 1.1  1992/05/26  18:34:08  mringe
 * Initial revision
 *
 */

#include <stdlib.h>
#include <string.h>
#include "meataxe.h"


/* ------------------------------------------------------------------
   permdesc() - Describe
   ------------------------------------------------------------------ */
    
char *permdesc(p)
perm_t *p;

{
    static char buf[50];
    sprintf(buf,"permutation of degree %ld",p->deg);
    return buf;
}

/* ------------------------------------------------------------------
   permprint() - Print
   ------------------------------------------------------------------ */

#define SIZE(i) ((i)<9?2 : (i)<99?3 : (i)<999?4 : (i)<9999?5 : \
	(i)<99999?6 : (i)>999999?7 : (i)<9999999?8 : 9)
    
void permprint(name,perm)
char *name;
perm_t *perm;

{
    long *p = (long *)(perm->d) - 1;
    int cycle = 1;
    int empty = 1;
    int count = 0, k;
    int i;

    if (name != NULL) printf("%s=",name);
    while (1)
    {
	/* Find next cycle
	   --------------- */
	while (cycle <= perm->deg && p[cycle] < 0) ++cycle;
	if (cycle > perm->deg) break;	/* Done */
	i = cycle;

	     /* Check if it is a fixed point (we don't
	        print fixed points, GAP doesn't like them)
	        ----------------------------------------- */
	    if (p[i] == i)
	    {
		p[i] = -p[i];	/* Mark it as done */
		continue;
	    }

	    /* Print cycle
	       ----------- */
	    if ((count += SIZE(i)) > 77)
	    {
		printf("\n    (%ld",(long)i);
	        count = 5 + SIZE(i);
	    }
	    else
		printf("(%ld",(long)i);
	    empty = 0;

	    while (1)
	    {
		k = i;
		i = p[i];
		p[k] = -p[k];
		if (p[i] < 0) break;

	        if ((count += SIZE(i)) > 77)
	        {
		    printf(",\n     %ld",(long)i);
	            count = 4 + SIZE(i);
	        }
	        else
		    printf(",%ld",(long)i);
	    }
	    printf(")");
	    ++count;
    }
    if (empty) printf("()");
    if (name != NULL) printf("\n");
    for (i = 1; i <= perm->deg; ++i) p[i] = -p[i];
}



/* ------------------------------------------------------------------
   permalloc() - Allocate a permutation
   permfree() - Free a permutation
   ------------------------------------------------------------------ */

perm_t *permalloc(deg)
long deg;

{
    perm_t *p;

    p = malloc(sizeof(perm_t));
    if (p == NULL)
    {
	MTXFAIL(ERR_NOMEM,NULL);
    }
    p->id = T_PERM;
    p->deg = deg;
    p->d = (PTR) malloc(deg * sizeof(long));
    if (p == NULL)
    {
	free(p);
	MTXFAIL(ERR_NOMEM,NULL);
    }
    return p;
}


void permfree(p)
perm_t *p;

{	free(p->d);
	free(p);
}


/* ------------------------------------------------------------------
   permdup() - Duplicate a permutation
   permmove() - Move a permutation
   ------------------------------------------------------------------ */

perm_t *permdup(src)
perm_t *src;

{
    perm_t *p;

    p = permalloc(src->deg);
    if (p == NULL) return NULL;
    memcpy(p->d,src->d,(size_t) src->deg * sizeof(long));
    return p;
}

perm_t *permmove(dest, src)
perm_t *dest, *src;

{
    if (dest->deg != src->deg)
	MTXFAIL(ERR_INCOMPAT,NULL);
    memcpy(dest->d,src->d,(size_t) src->deg * sizeof(long));
    return dest;
}


/* ------------------------------------------------------------------
   permread() - Read a permutation
   permwrite() - Write a permutation
   ------------------------------------------------------------------ */
   
perm_t *permread(f)
FILE *f;

{
    perm_t *p;
    long hdr[3];

    if (zreadlong(f,hdr,3) != 3) MTXFAIL(ERR_FILEREAD,NULL);
    if (hdr[0] != -1) MTXFAIL(ERR_NOTPERM,NULL);
    p = permalloc(hdr[1]);
    if (p == NULL) return NULL;
    if (zreadlong(f,(long *)p->d,hdr[1]) != hdr[1])
    {
	permfree(p);
	MTXFAIL(ERR_FILEREAD,NULL);
    }
    return p;
}


int permwrite(f, perm)
FILE *f;
perm_t *perm;

{
    long hdr[3];

    hdr[0] = -1;
    hdr[1] = perm->deg;
    hdr[2] = 1;
    if (zwritelong(f,hdr,3) != 3) 
	MTXFAIL(ERR_FILEWRITE,-1);
    if (zwritelong(f,(long *)perm->d,hdr[1]) != hdr[1])
    	MTXFAIL(ERR_FILEWRITE,-1);
    return 0;
}



/* ------------------------------------------------------------------
   permload() - Read a permutation
   permsave() - Write a permutation
   ------------------------------------------------------------------ */

perm_t *permload(fn)
char *fn;

{
    FILE *f;
    perm_t *p;

    if ((f = SFOpen(fn,FM_READ)) == NULL)
    {
	perror(fn);
	errexit(ERR_FILEREAD,fn);
    }
    p = permread(f);
    if (p == NULL) errexit(ERR_FILEREAD,fn);
    fclose(f);
    return p;
}


int permsave(perm, fn)
perm_t *perm;
char *fn;

{
    FILE *f;
    int result;

    if ((f = SFOpen(fn,FM_CREATE)) == NULL)
    {
	perror(fn);
	errexit(ERR_FILEREAD,fn);
    }
    result = permwrite(f,perm);
    fclose(f);
    if (result != 0) errexit(ERR_FILEREAD,fn);
    return result;
}


/* ------------------------------------------------------------------
   permmul() - Multiply permutations
   ------------------------------------------------------------------ */

perm_t *permmul(dest, src)
perm_t *dest;
perm_t *src;

{
    long i;
    long *d = (long *)dest->d;
    long *s = ((long *)src->d)-1;

    if (dest->deg != src->deg)
    {
	MTXFAIL(ERR_INCOMPAT,NULL);
    }
    for (i = dest->deg; i > 0; --i)
    {
	*d = s[*d];
	++d;
    }
    return dest;
}


/* ------------------------------------------------------------------
   permorder() - Order of a permutation
   ------------------------------------------------------------------ */

long permorder(perm)
perm_t *perm;

{
    long ord = 1;
    long deg = perm->deg;
    long *p = ((long *)perm->d) -1;
    long *end = p + (deg+1);
    long start = 1;
    register long *seed = p, *x, tord;

    while (1)
    {
	/* Find next starting point
	   ------------------------ */
	while (*++seed < 0 && seed != end);
	if (seed == end) break;

	/* Find orbit size
	   --------------- */
	tord = 0;
	x = seed;
	while (1)
	{
	    register long tmp = *x;
	    if (tmp < 0) break;
	    *x = -tmp;
	    x = p + tmp;
	    ++tord;
	}

	/* Calculate l.c.m of orbit sizes
	   ------------------------------ */
	{
	    register long i = ord, j = tord, k;
	    while (1)
	    {	if ((i %= j) == 0)
		{	k = j; break;
		}
		if ((j %= i) == 0)
		{	k = i; break;
		}
	    }
	    ord = ord / k * tord;
	}
    }

    /* Restore the permutation
       ----------------------- */
    for (x = p + 1; x != end; ++x)
	*x = -*x;

    return ord;
}


/* ------------------------------------------------------------------
   permpower() - Power of a permutation
   ------------------------------------------------------------------ */

perm_t *permpower(p,n)
perm_t *p;
long n;

{	perm_t *q = permalloc(p->deg);
	long *xp = (long *)(p->d) - 1;
	long *xq = (long *)(q->d) - 1;
	long i, k, l;

	for (i = 1; i <= p->deg; ++i)
	{	for (k = i, l = n; l > 0; --l)
			k = xp[k];
		xq[i] = k;
	}
	return q;
}




