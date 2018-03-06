/* ========================== C MeatAxe =============================
   words.c - The word generator.

   (C) Copyright 1994 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: words.c,v 1.2 1997/09/11 15:43:26 gap Exp $
 *
 * $Log: words.c,v $
 * Revision 1.2  1997/09/11 15:43:26  gap
 * New version 2.2.3. AH
 *
 * Revision 2.17  1995/02/08  10:14:52  mringe
 * _PL entfernt.
 *
 * Revision 2.17  1995/02/08  10:14:52  mringe
 * _PL entfernt.
 *
 * Revision 2.16  1994/08/29  17:04:26  mringe
 * MakeWordF() wieder rausgeworfen.
 *
 * Revision 2.15  1994/08/25  12:53:36  mringe
 * Neu: MakeWordF()
 *
 * Revision 2.14  1994/08/25  07:36:53  mringe
 * Speicherleck in mkbasis() beseitigt.
 *
 * Revision 2.13  1994/07/29  07:43:08  mringe
 * Neuer Wortgenerator.
 *
 * Revision 2.12  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.11  1994/07/23  16:55:11  mringe
 * Profiling, makeword() gibt jetzt das Wort zurueck.
 *
 * Revision 2.10  1994/06/27  14:11:37  mringe
 * Bug in nameof() behoben.
 *
 * Revision 2.9  1994/06/25  13:55:41  mringe
 * *** empty log message ***
 *
 * Revision 2.8  1994/06/25  13:45:09  mringe
 * Memory leak in nameof().
 *
 * Revision 2.7  1994/06/19  16:14:24  mringe
 * words.h eliminiert.
 *
 * Revision 2.6  1994/06/16  14:16:28  mringe
 * Profiling
 *
 * Revision 2.5  1994/05/11  18:23:00  mringe
 * initbasis(): Reihenfolge der Parameter.
 *
 * Revision 2.4  1994/05/10  12:37:26  mringe
 * Unbenutzte variablen entfernt.
 *
 * Revision 2.3  1994/05/05  14:16:37  mringe
 * Multipl. mit Erzeuger wieder rausgenommen.
 *
 * Revision 2.2  1994/05/04  09:20:12  mringe
 * Bilde Linearkomb. mit Koeffizienten aus dem ganzen Koerper,
 * nicht nur aus dem Primkoerper.
 *
 * Revision 2.1  1994/04/11  14:35:06  mringe
 * Verbesserte Erzeugerbildung.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.15  1993/08/12  16:14:23  mringe
 * Include <string.h>.
 *
 * Revision 1.14  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.13  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.12  1993/02/13  11:46:41  mringe
 * Bug in mkbasisname() behoben.
 *
 * Revision 1.11  1993/02/13  11:34:45  mringe
 * nameof() von ZSM hierher verlagert.
 *
 * Revision 1.10  1993/02/12  17:06:59  mringe
 * Woerter mit N Erzeugern.
 *
 * Revision 1.9  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.8  1993/01/18  16:54:27  mringe
 * Fehler beseitigt.
 *
 * Revision 1.7  1993/01/18  16:51:36  mringe
 * z3 entfernt. Berechne Fingerprint via mkword().
 *
 * Revision 1.6  1993/01/08  00:24:45  mringe
 * makeword() verbessert. Macht jetzt keine Annahme
 * ueber vorangegangene Aufrufe mehr.
 *
 * Revision 1.5  1992/09/03  12:53:44  mringe
 * Bug in mkword() behoben.
 *
 * Revision 1.4  1992/08/31  14:40:32  mringe
 * Worte -1..-6 (=Standard-Fingerprint)
 *
 * Revision 1.3  1992/07/22  07:10:30  mringe
 * Changed 'global.h' to 'lattice.h'
 *
 * Revision 1.2  1992/06/01  08:34:25  mringe
 * CL warnings entfernt.
 *
 * Revision 1.1  1992/05/24  09:43:36  mringe
 * Initial revision
 *
 */

#include <string.h>
#include <stdlib.h>

#include "meataxe.h"



static int mkbasis(wgdata_t *b, int max);
static char *basisname(int i, int ngen);



/* ------------------------------------------------------------------
   WGInit() - Initialize the word generator. Returns a pointer to
   a basis_t struncture or NULL on error.
   ------------------------------------------------------------------ */

wgdata_t *WGInit(ngen,gen)
int ngen;
matrix_t *gen[];

{
    int i, k;
    long fl, dim;
    wgdata_t *b;

    /* Check the generators
       -------------------- */
    if (ngen < 1) MTXFAIL(ERR_RANGE,NULL);
    fl = gen[0]->fl;
    dim = gen[0]->nor;
    for (i = 0; i < ngen; ++i)
    {
    	if (gen[i]->fl != fl || gen[i]->nor != dim)
    	    MTXFAIL(ERR_INCOMPAT,NULL);
    	if (gen[i]->nor != gen[i]->noc)
    	    MTXFAIL(ERR_NOTSQUARE,NULL);
    }

    /* Create a new basis_t structure
       ------------------------------  */
    b = ALLOC(wgdata_t);
    if (b == NULL) MTXFAIL(ERR_NOMEM,NULL);
    b->b = NALLOC(matrix_t*,ngen + 10);
    if (b->b == NULL) MTXFAIL(ERR_NOMEM,NULL);
    b->max = ngen + 10;
    b->len = ngen;
    b->ngen = ngen;
    for (k = 0; k < ngen; ++k)
    	b->b[k] = gen[k];
    return (void *) b;
}


/* ------------------------------------------------------------------
   WGFree() - Free a basis_t structure
   ------------------------------------------------------------------ */

int WGFree(b)
wgdata_t *b;

{
    int k;

    /* Check the handle
       ---------------- */
    if (b == NULL || b->ngen < 0 || b->ngen > 1000 || b->len < b->ngen)
    	MTXFAIL(ERR_BADARG,-1);

    /* Free all matrices except the generators
       --------------------------------------- */
    for (k = b->ngen; k < b->len; ++k)
    	matfree(b->b[k]);
    free(b->b);
    free(b);
    return 0;
}


/* ------------------------------------------------------------------
   mkbasis() - Calculate basis elements up to number <max>
   ------------------------------------------------------------------ */

#define MAXCP 30
static int cptab[MAXCP] =
      {0 ,1 ,2,2,3,3,2,1,1, 1,5,1,3,7,3,3,3,3,3,3,3,2,4,3,3,1,4,4,3,1};
static int mul[MAXCP] =
      {-1,-1,1,7,6,3,1,5,1, 0,2,5,1,4,7,3,5,1,4,7,0,6,3,9,1,4,3,2,7,1};
 
static int cp(i)
int i;
{
    return i < MAXCP ? i - cptab[i] : (i * 87) % i;
}

static int mu(i,ngen)
int i;
int ngen;
{
    return (i < MAXCP ? mul[i] : (i * 23)) % ngen;
}


static int mkbasis(b,max)
wgdata_t *b;
int max;

{
    int i;

    if (max >= b->max)
    {
    	matrix_t **x = NREALLOC(b->b,matrix_t *,max+5);
    	if (x == NULL) return -1;
	b->b = x;
	b->max = max + 5;
    }
    for (i = b->len; i <= max; ++i)
	b->b[i] = matmul(matdup(b->b[cp(i)]),b->b[mu(i,b->ngen)]);
    if (max >= b->len)
    	b->len = max + 1;
    return 0;
}


/* ------------------------------------------------------------------
   mapnumber() - Map word number (to get some special words first)
   ------------------------------------------------------------------ */

static long fpnum[MAXFP+1] = {-1,7,39,217,219,223,249};

static long mapnumber(n)
long n;

{
    int i;
    if (n <= MAXFP && n >= 1)
	n = fpnum[n];
    else
	for (i = 1; i <= MAXFP; ++i)
	    if (n == fpnum[i]) n = i;
    return n;
}


/* ------------------------------------------------------------------
   MakeWord() - Calculate a word
   ------------------------------------------------------------------ */

matrix_t *MakeWord(b, n)
wgdata_t *b;
long n;

{
    long l;
    int i;
    matrix_t *w;
    PROFILE_BEGIN(t);

    /* Initialize
       ---------- */
    if (n <= 0) MTXFAIL(ERR_RANGE,NULL);
    n = mapnumber(n);
    w = NULL;

    /* Calculate the word
       ------------------ */
    for (i = 0, l = n; l != 0; ++i, l >>= 1)
    {
	if ((l % 2) == 0) continue;
	if (i >= b->len) mkbasis(b,i);
	if (w != NULL)
	    matadd(w,b->b[i]);
	else
	    if ((w = matdup(b->b[i])) == NULL) return NULL;
    }
    PROFILE_END(t,MakeWord);
    return w;
}


/* ------------------------------------------------------------------
   makefp() - Calculate standard fingerprint
   ------------------------------------------------------------------ */
#if (MAXFP != 6)
	************* FEHLER *****************
#endif

void makefp(b,fp)
wgdata_t *b;
long fp[];

{	
    int i;

    for (i = 1; i <= MAXFP; ++i)
    	fp[i-1] = nullity__(MakeWord(b,i));
}



/* ------------------------------------------------------------------
   basisname() - Symbolic representation of basis elements.
   ------------------------------------------------------------------ */

static char *basisname(i,ngen)
int i;
int ngen;

{
    static char buf[50];
    char b[2];

    if (i < ngen)
    {
    	buf[0] = (char) (i + 'a');
    	buf[1] = 0;
    	return buf;
    }

    basisname(cp(i),ngen);
    b[0] = (char) (mu(i,ngen) + 'a');
    b[1] = 0;
    strcat(buf,b);
    return buf;
}


/* ------------------------------------------------------------------
   SymbolicName() - Returns the name of a word (given its number).
   ------------------------------------------------------------------ */

char *SymbolicName(b,n)
wgdata_t *b;
long n;

{
    static char name[200];
    char *c = name;
    int i;
    long digit;

    /* Check the arguments
       ------------------- */
    if (n <= 0)
    	MTXFAIL(ERR_RANGE,NULL);

    n = mapnumber(n);
    for (i = 0; n != 0; ++i)
    {	digit = n % 2;	/* Get next digit */
	n /= 2;
	if (digit == 0) continue;
	if (c != name) *c++ = '+';
	strcpy(c,basisname(i,b->ngen));
	while (*c != 0) ++c;
    }
    return name;
}


