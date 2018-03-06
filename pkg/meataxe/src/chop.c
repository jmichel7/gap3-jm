/* ========================== C MeatAxe =============================
   chop.c - Chop module (Calculate composition series)

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: chop.c,v 1.2 1997/09/11 15:42:43 gap Exp $
 *
 * $Log: chop.c,v $
 * Revision 1.2  1997/09/11 15:42:43  gap
 * New version 2.2.3. AH
 *
 * Revision 2.43  1995/05/29  14:10:17  mringe
 * Fehlermeldung modifiziert.
 *
 * Revision 2.42  1995/02/08  10:01:14  mringe
 * PL_ entfernt.
 *
 * Revision 2.41  1994/09/22  08:19:36  mringe
 * Berechne char.Pol in zwei Schritten.
 *
 * Revision 2.40  1994/09/14  09:59:38  mringe
 * Benutze variablen Seed beim spinup.
 *
 * Revision 2.39  1994/08/22  08:27:03  mringe
 * Berechne den Seed fuer dual split per spin-up.
 *
 * Revision 2.39  1994/08/22  08:27:03  mringe
 * Berechne den Seed fuer dual split per spin-up.
 *
 * Revision 2.38  1994/08/15  09:05:46  mringe
 * *** empty log message ***
 *
 * Revision 2.37  1994/07/29  12:48:34  mringe
 * Neuer Algorithmus, vermeidet Einsetzen der Martrix in das Polynom.
 *
 * Revision 2.35  1994/07/29  07:43:08  mringe
 * Neuer Wortgenerator.
 *
 * Revision 2.34  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.33  1994/07/23  16:53:40  mringe
 * Profiling, neue Nullraum- und Echelonfunktionen.
 *
 * Revision 2.32  1994/07/07  05:36:12  mringe
 * ULTRIX_42 in ULTRIX umbenannt
 *
 * Revision 2.31  1994/07/05  15:27:10  mringe
 * Option -d und -n
 *
 * Revision 2.30  1994/06/27  14:11:05  mringe
 * Erlaube -g 1.
 *
 * Revision 2.29  1994/06/22  14:38:40  mringe
 * Option -p: keine Polynome.
 *
 * Revision 2.28  1994/06/19  16:15:02  mringe
 * words.h eliminiert.
 *
 * Revision 2.27  1994/06/16  19:19:01  mringe
 * Statistic verbessert.
 *
 * Revision 2.26  1994/06/16  14:20:17  mringe
 * Profiling, diverse kleinere Aenderungen.
 *
 * Revision 2.25  1994/05/19  11:35:54  mringe
 * Bad words via set_t implementiert.
 *
 * Revision 2.24  1994/05/18  11:52:50  mringe
 * Nach vielen Aenderungen... lauffaehige Version.
 *
 * Revision 2.23  1994/05/11  18:22:18  mringe
 * Einige Aenderungen in chop(): Reihenfolge der Faktoren, etc.
 *
 * Revision 2.22  1994/05/06  15:21:36  mringe
 * Benutze den gesamten Nullraum als seed.
 * Falls ein irred. Faktor von c(x) vorliegt, braucht man
 * beim Norton-Test nur einen Vektor (?)
 *
 * Revision 2.21  1994/05/05  14:17:05  mringe
 * Bug in chop() behoben. Test auf c(x) irred. jetzt o.k.
 *
 * Revision 2.20  1994/05/04  15:25:44  mringe
 * Bug in equiv() behoben.
 *
 * Revision 2.19  1994/05/04  09:29:26  mringe
 * Teste, ob die Erzeuger ueber dem gleichen Koerper sind.
 *
 * Revision 2.18  1994/05/04  08:46:31  mringe
 * Polynome.
 *
 * Revision 2.17  1994/04/09  13:51:23  mringe
 * memset()-Bug ???
 *
 * Revision 2.16  1994/03/26  06:34:04  mringe
 * basename umbenannt wg. Namenskonflikt.
 *
 * Revision 2.15  1994/03/02  11:06:14  mringe
 * Vermeide mehrfache Berechnung von Nullitaeten.
 *
 * Revision 2.14  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.13  1994/01/28  08:24:38  mringe
 * Schreibe die Erzeuger sofort raus, nicht erst am Ende.
 *
 * Revision 2.12  1993/12/14  22:38:12  mringe
 * Neu: Funktion gcd().
 *
 * Revision 2.11  1993/12/13  08:27:53  mringe
 * split(): Reihenfolge der Arguimente.
 *
 * Revision 2.10  1993/12/13  08:25:53  mringe
 * Reihenfolge der Fkt.-argumente vereinheitlicht.
 *
 * Revision 2.9  1993/12/08  11:48:50  mringe
 * Compiler warnings.
 *
 * Revision 2.8  1993/12/08  11:33:02  mringe
 * Neue CPU time - Funktionen.
 *
 * Revision 2.7  1993/12/07  16:44:24  mringe
 * Schreibe erst .cfinfo, dann die irreduziblen.
 * (Sonst funktioniert cfname() nicht richtig).
 *
 * Revision 2.6  1993/10/27  09:34:47  mringe
 * IdWords.
 *
 * Revision 2.5  1993/10/22  16:50:49  mringe
 * *** empty log message ***
 *
 * Revision 2.4  1993/10/22  16:08:19  mringe
 * Neues Numerierungsschema fuer irreduzible.
 *
 * Revision 2.3  1993/10/21  08:37:53  mringe
 * Compiler warnings.
 *
 * Revision 2.2  1993/10/20  18:26:42  mringe
 * Bug bei matload() behoben.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.44  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.43  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.42  1993/10/02  16:23:02  mringe
 * matread() und matwrite() in matload() bzw. matsave() umbenannt.
 *
 * Revision 1.41  1993/09/20  19:47:59  mringe
 * *** empty log message ***
 *
 * Revision 1.40  1993/09/20  17:40:14  mringe
 * msg_prefix gesetzt.
 *
 * Revision 1.39  1993/08/27  15:27:26  mringe
 * Option -T
 *
 * Revision 1.38  1993/08/25  15:57:32  mringe
 * Verwaltung der Worte verbessert.
 *
 * Revision 1.37  1993/08/10  15:52:06  mringe
 * Test nullword(n) statt nullword(n->parent).
 *
 * Revision 1.36  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.35  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.34  1993/07/28  06:06:12  mringe
 * Probiere erst alle Woerter, die schon einmal zum Nachweis der
 * Irreduzibilitaet benutzt wurden (nicht nur bei passender Dimension).
 *
 * Revision 1.33  1993/07/23  13:46:27  mringe
 * OS-Symbole neu (SYS_xxx)
 *
 * Revision 1.32  1993/07/19  15:00:29  mringe
 * Option -G (GAP output).
 *
 * Revision 1.31  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.30  1993/02/16  18:32:46  mringe
 * string.h und stdio.h werden jetzt in meataxe.h included.
 *
 * Revision 1.29  1993/02/15  13:49:24  mringe
 * Funktionen aus cyclic.c -> yyy-Lib verschoben.
 *
 * Revision 1.28  1993/02/12  17:06:59  mringe
 * Woerter mit N Erzeugern.
 *
 * Revision 1.27  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.26  1993/02/10  19:23:03  mringe
 * !comp
 *
 * Revision 1.25  1993/01/28  07:35:51  mringe
 * Benutze "utils" Bibliothek
 *
 * Revision 1.24  93/01/19  07:21:10  07:21:10  mringe (  Michael Ringe)
 * Bug in checkspl() entfernt.
 * 
 * Revision 1.23  1993/01/15  14:50:01  mringe
 * ngen in nach files.h ausgelagert.
 *
 * Revision 1.22  1993/01/15  07:30:49  mringe
 * Erweiterung auf N Erzeuger.
 *
 * Revision 1.21  1993/01/09  14:50:26  mringe
 * `Save vectors' implementiert.
 *
 * Revision 1.20  1993/01/08  19:50:08  mringe
 * Benutze auch den Fingerprint bei der Suche nach einem
 * Wort mit minimaler Nullitaet.
 *
 * Revision 1.19  1993/01/08  10:53:51  mringe
 * Suchalgorithmus verbessert: (1) beruechsichtige die schon
 * gefundenen irreduziblen und (2) benutze die Worte, mit denen
 * vorher schon einal gesplittet wurde.
 *
 * Revision 1.18  1993/01/07  20:12:44  mringe
 * Teste zuerst, ob es einer der schon bekannten Irred. ist.
 *
 * Revision 1.17  1993/01/06  20:59:57  mringe
 * getopt in zgetopt() umbenannt.
 *
 * Revision 1.16  1992/10/13  18:00:59  mringe
 * No Changes
 *
 * Revision 1.15  1992/10/12  10:40:30  hiss
 * Neu: Option -s
 *
 * Revision 1.14  1992/10/02  16:40:00  mringe
 * Merke nullit"at=0 auch bei W"ortern -1..-6.
 *
 * Revision 1.13  1992/10/01  13:50:10  mringe
 * Header eingef"ugt.
 *
 * Revision 1.12  1992/09/01  13:01:11  mringe
 * Bug in chop() beseitigt.
 *
 * Revision 1.11  1992/08/31  14:40:32  mringe
 * Benutze Worte -1..-6 zum Splitten.
 *
 * Revision 1.10  1992/07/29  08:27:07  mringe
 * Debug messages in extendbasis() entfernt.
 *
 * Revision 1.9  1992/07/28  08:33:55  mringe
 * Cleaned up checkspl() and findw1().
 *
 * Revision 1.8  1992/07/28  08:20:29  mringe
 * Option -m erweitert (Startwert.Maximum)
 *
 * Revision 1.7  1992/07/22  07:10:30  mringe
 * Changed 'global.h' to 'lattice.h'
 *
 * Revision 1.6  1992/07/15  09:25:55  mringe
 * Some minor changes.
 *
 * Revision 1.5  1992/07/08  15:28:31  mringe
 * Erkennung des Zerf.k"orpers verbessert.
 *
 * Revision 1.4  1992/07/02  15:42:07  mringe
 * Sichere Erkennung des Zerf"allungsk"orpers.
 *
 * Revision 1.3  1992/07/01  14:43:09  mringe
 * Einige type casts eingef"ugt.
 *
 * Revision 1.2  1992/06/01  08:33:30  mringe
 * CL warnings entfernt.
 *
 * Revision 1.1  1992/05/23  16:46:00  mringe
 * Initial revision
 */


#include <string.h>
#include <stdlib.h>
#include "meataxe.h"
#include "lattice.h"
#include "files.h"


/* ------------------------------------------------------------------
   Some constants
   ------------------------------------------------------------------ */

#define MAXTRIES 2500	/* Number of tries before findw1() fails */


/* ------------------------------------------------------------------
   typedefs
   ------------------------------------------------------------------ */

typedef struct nodestruct
{
    struct nodestruct *sub, *quot, *parent;
    long dim, num;		/* Dimension and number */
    matrix_t *gen[MAXGEN];	/* Generators */
    long mult;			/* Multiplicity */
    long spl;			/* Degree of splitting field */
    long idword;		/* Word used for std basis */
    poly_t *idpol;		/* Polynomial "  "    "    */
    poly_t *f1, *f2;		/* characteristic Polynomial c=f1*f2*/
    long fprint[MAXFP];		/* Fingerprint */
    set_t *badwords;
    matrix_t *nsp, *nsptr;	/* Null space */
    matrix_t *subsp;		/* Subspace */
    matrix_t **gentr;
    long ggt;			/* g.c.d. of nullities */
    void *wg;			/* Used by word generator */
    long wnum;
    matrix_t *word;
}
node_t;


/* ------------------------------------------------------------------
   Function prototypes
   ------------------------------------------------------------------ */

static void init(void);
static int isbadword(long w,node_t *n);
static void writetree(node_t *n);
static void writeinfo(node_t *root);
static void splitnode(node_t *n, int tr);
static node_t *chop(node_t *n);
static matrix_t *extendbasis(matrix_t *basis, matrix_t *space);
static int checkspl(matrix_t *gen[], matrix_t *nsp);
static int findidword(node_t *n);
static int equiv(node_t *n1, node_t *n2);
static void newirred(node_t *n);


/* ------------------------------------------------------------------
   Global variables
   ------------------------------------------------------------------ */

node_t *root;
long opt_deglimit = -1;
long opt_nullimit;
node_t *irred[MAXCF];			/* List of irreducibles */
long firstword = -1;
int opt_G = 0;				/* GAP output */
set_t *goodwords;			/* List of `good' words */

static long stat_svsplit = 0;
static long stat_cpirred = 0;
static long stat_nssplit = 0;
static long stat_dlsplit = 0;
static long stat_irred = 0;



static char *helptext[] = {
"SYNTAX",
"    chop [<Options>] <Name>",
"",
"OPTIONS",
"    -G               GAP output (implies -Q).",
"    -Q               Quiet, no messages.",
"    -V               Verbose, more messages.",
"    -g <NGen>        Set number of generators (default is 2).",
"    -n <MaxNul>      Set limit on nullity",
"    -d <MaxDeg>      Set limit on degrees of polynomials",
"    -T <MaxTime>     Set CPU time limit",
"",
"FILES",
"    The input must be <NGen> square matrices. They are read from",
"    <Name>.1, <Name>.2, ... <Name>.<NGen>.",
"",
"    Output is written to various files:",
"      <Name>.cfinfo        Information about the composition factors.",
"      <Name><Dim><a>.<N>   Generators on the composition factors.",
"                           <Dim> is the dimension, <a> is a letter",
"                           discriminating different factors, and <N>",
"                           runs from 1 to <NGen>.",
NULL};

static proginfo_t pinfo =
   { "chop", "Composition Series", "$Revision: 1.2 $", helptext };


/* ------------------------------------------------------------------
   newnode() - Create a new node
   ------------------------------------------------------------------ */

static node_t *newnode(gen,parent)
matrix_t *gen[];
node_t *parent;

{
    node_t *n = ALLOC(node_t);
    n->parent = parent;
    n->sub = n->quot = NULL;
    n->dim = gen[0]->nor;
    n->num = -1;
    n->mult = -1;
    memcpy(n->gen,gen,(size_t)ngen*sizeof(matrix_t *));
    n->spl = -1;
    n->idword = -1;
    n->idpol = NULL;
    n->badwords = set_alloc();
    n->nsp = n->nsptr = NULL;
    n->subsp = NULL;
    n->ggt = 0;
    n->f1 = n->f2 = NULL;

    n->wg = WGInit(ngen,gen);
    if (n->wg == NULL) errexit(-1,"WGInit");
    n->wnum = -1;
    n->word = NULL;
    n->gentr = NULL;

    return n;
}

#define MATFREE(x) ( ((x) != NULL) ? (matfree(x), (x)=NULL) : NULL )


/* ------------------------------------------------------------------
   makeword() - Make a word
   ------------------------------------------------------------------ */

static void makeword(n,w)
node_t *n;
long w;

{
    if (n->word != NULL) matfree(n->word);
    n->word = MakeWord(n->wg,w);
    n->wnum = w;
}

/* ------------------------------------------------------------------
   insertword() - Insert the word into a polynomial and calculate the
   null-space (nsp) and the null-space of p(A)^T (nsptr).
   ------------------------------------------------------------------ */

static void insertword(n,p)
node_t *n;
poly_t *p;

{
    matrix_t *m, *mt;

    m = matinsert(n->word,p);
    mt = mattr(m);
    if (n->nsp != NULL) matfree(n->nsp);
    n->nsp = nullspace__(m);
    if (n->nsptr != NULL) matfree(n->nsptr);
    n->nsptr = nullspace__(mt);
}


/* ------------------------------------------------------------------
   init () - Read generators
   ------------------------------------------------------------------ */

static void init()

{
    char fn[100];
    char fn0[200];
    int i;
    matrix_t *x;
    matrix_t *gen[MAXGEN];

    for (i = 0; i < ngen; ++i)
    {
	sprintf(fn,"%s.%d",cfbasename,i+1);
	if (i == 0) strcpy(fn0,fn);
	x = gen[i] = matload(fn);
	if (x->nor != x->noc)
	    errexit(ERR_NOTSQUARE,fn);
	if (i > 0 && (x->nor != gen[0]->nor ||
		x->fl != gen[0]->fl))
	    errexit(ERR_INCOMPAT,strcat(strcat(fn0," and "),fn));
    }
    root = newnode(gen,NULL);
}


/* ------------------------------------------------------------------
   isbadword() - Look up a word in the 'badwords' list
   ------------------------------------------------------------------ */

static int isbadword(w,n)
long w;
node_t *n;

{
    if (n == NULL) return 0;
    if (n->badwords != NULL && set_contains(n->badwords,w))
	return 1;
    return isbadword(w,n->parent);
}



/* ------------------------------------------------------------------
   writetree() - Write the composition series to stdout
   ------------------------------------------------------------------ */

static void writetree(n)
node_t *n;

{
    static int count = 0;	/* Characters written in one line */
    static int first = 1;
    int i;

    if (n->sub == NULL)		/* Irreducible */
    {
	if (count >= 20)
	{   printf("\n");
	    count = 0;
	}
	for (i = 0; cfinfo[i].dim != n->dim || cfinfo[i].num != n->num;
	    ++i);
	if (opt_G)
	{
	    if (!first) printf(","); else { printf(" "); first = 0; }
	    printf("%d",i+1);
	}
	else
	    printf("%s ",cfname(i));
	++count;
    }
    else
    {
	writetree(n->sub);
	writetree(n->quot);
    }
}


/* ------------------------------------------------------------------
   writeinfo() - Write some information to stdout and create the
	.cfinfo file.
   ------------------------------------------------------------------ */

static void writeinfo(root)
node_t *root;

{
    int i, k;

    MESSAGE(0,(
    "\n\nChopping completed: %d different composition factors\n",ncf));

    /* Write the cfinfo file
       --------------------- */
    MESSAGE(0,("Writing %s.cfinfo\n",cfbasename));
    for (i = 0; i < ncf; ++i)
    {
    	cfinfo[i].dim = irred[i]->dim;
	cfinfo[i].num = irred[i]->num;
	cfinfo[i].mult = irred[i]->mult;
	cfinfo[i].spl = irred[i]->spl;
	cfinfo[i].idword = irred[i]->idword;
	cfinfo[i].idpol = irred[i]->idpol;
    }
    writecfinfo();


    /* Write composition factors
       ------------------------- */
    if (opt_G)
    {
	printf("MeatAxe.CompositionFactors := [\n");
    	for (i = 0; i < ncf; ++i)
    	{
	    printf("  [ \"%s\", %ld, %ld ]",cfname(i),
	      irred[i]->mult,irred[i]->spl);
	    if (i < ncf-1) printf(",");
	    printf("\n");
	}
	printf("\n];\n");
    }
    else
    {
    	printf("\nName   Mult  SF  Fingerprint\n");
    	for (i = 0; i < ncf; ++i)
    	{
	    printf("%-6s %4ld  %2ld  ",cfname(i),
	      irred[i]->mult,irred[i]->spl);
	    for (k = 0; k < MAXFP; ++k)
		printf("%ld%s",irred[i]->fprint[k],k==MAXFP-1?"\n":",");
	}
    }

    /* Write the composition series
       ---------------------------- */
    if (opt_G)
	printf("MeatAxe.CompositionSeries := [\n");
    else
	printf("\nAscending composition series:\n");
    writetree(root);
    if (opt_G)
	printf("];\n");
    else
	printf("\n");

    /* Write statistics
       ---------------- */
    if (MSG3)
    {
	printf("\nStatistics\n----------\n");

	printf("Saved vectors split: %4ld\n",stat_svsplit);
	printf("c(x) irreducible:    %4ld\n",stat_cpirred);
	printf("Normal split:        %4ld\n",stat_nssplit);
	printf("Dual split:          %4ld\n",stat_dlsplit);
	printf("Irreducible:         %4ld\n",stat_irred);

	printf("\nTimings: makeword=%8ld\n",ProfMakeWord);
	printf("         chpol   =%8ld\n",ProfCharPol);
	printf("         spin    =%8ld\n",ProfSpinUp);
	printf("         nullsp  =%8ld\n",ProfNullSpace);
	printf("         insert  =%8ld\n",ProfMatInsert);
    }
}


static void cleannode(node_t *n)

{   int i;
    if (n->gentr != NULL)
    {
        for (i = 0; i < ngen; ++i)
	    matfree(n->gentr[i]);
	free(n->gentr);
	n->gentr = NULL;
    }
    if (n->subsp != NULL) { matfree(n->subsp); n->subsp = NULL; }
    if (n->wg != NULL) { WGFree(n->wg); n->wg = NULL; }
    if (n->nsp != NULL) { matfree(n->nsp); n->nsp = NULL; }
    if (n->nsptr != NULL) { matfree(n->nsptr); n->nsptr = NULL; }
    if (n->word != NULL) { matfree(n->word); n->word = NULL; }
    if (n->f1 != NULL) { polfree(n->f1); n->f1 = NULL; }
    if (n->f2 != NULL) { polfree(n->f2); n->f2 = NULL; }
}

static void matstat(matrix_t *gen[])
{
    int i;
    long count = 0, count0 = 0;

    for (i = 0; i < ngen; ++i)
    {
	long dim = gen[i]->nor;
	long k;
	PTR x = gen[i]->d;
	zsetlen(dim);
	for (k = 1; k <= dim; ++k)
	{
	    long j;
	    for (j = 1; j <= dim; ++j)
	    	if (zextract(x,j) == F_ZERO) ++count0;
	    count += dim;
	}
	zadvance(&x,1);
    }
    MESSAGE(1,("count = %ld, count0 = %ld\n",count,count0));

}

/* ------------------------------------------------------------------
   splitnode() - Split a node
   ------------------------------------------------------------------ */

static void splitnode(n,tr)
node_t *n;		/* Node to split */
int tr;			/* Indicates that it was a `dual split' */

{	
    matrix_t **sub, **quot;
    int i;

    /* Split the module
       ---------------- */
    MESSAGE(0,("Split: Subspace=%ld, Quotient=%ld\n",
	n->subsp->nor,n->dim - n->subsp->nor));
    sub = NALLOC(matrix_t *,ngen);
    quot = NALLOC(matrix_t *,ngen);
    split(n->subsp,ngen,tr ? n->gentr : n->gen,sub,quot);
matstat(sub);
matstat(quot);

    /* If it was a dual split, subspace and quotient have been
       calculated in the dual module. To get back to the original
       module, transpose again and exchange sub and quot.
       ---------------------------------------------------------- */
    if (tr)
    {	int i;
	matrix_t *x, *y;
	for (i = 0; i < ngen; ++i)
	{
	    x = mattr(sub[i]);
	    matfree(sub[i]);
	    y = mattr(quot[i]);
	    matfree(quot[i]);
	    sub[i] = y;
	    quot[i] = x;
	}
    }
    
    /* Make new nodes for subspace and quotient
       ---------------------------------------- */
    n->sub = newnode(sub,n);
    n->quot = newnode(quot,n);
    free(sub);
    free(quot);

    /* Project saved vectors on the quotient
       ------------------------------------- */
    if (!tr && n->nsp != NULL)
	n->quot->nsp = quotproj(n->subsp,n->nsp);

    /* Clean up
       -------- */
    cleannode(n);
    for (i = 0; i < ngen; ++i)
    {
	matfree(n->gen[i]);
	n->gen[i] = NULL;
    }

    /* Chop the subspace and quotient
       ------------------------------ */
    chop(n->sub);
    chop(n->quot);
}



/* ------------------------------------------------------------------
   chop_saved_vectors() - Chop a module using saved vectors. Returns
   1 on success, 0 otherwise.
   ------------------------------------------------------------------ */

static int chop_saved_vectors(node_t *n)

{
    int result;

    if (n->nsp == NULL || n->nsp->nor == 0) return 0;
    MESSAGE(1,("Trying saved vectors..."));

    if (n->subsp != NULL) matfree(n->subsp);
    result = spinup(n->nsp,ngen,n->gen,SPL_SEED_EACH,&(n->subsp));
    if (result == SPL_SUCCESS)
    {
	MESSAGE(1,("\n"));
	++stat_svsplit;
	splitnode(n,0);
    }
    else
	MESSAGE(1,("failed\n"));
    return (result == SPL_SUCCESS);
}


/* ------------------------------------------------------------------
   polymap() - Calculate the image of a single vector -- or a set
   of vectors -- under p(A).
   ------------------------------------------------------------------ */

matrix_t *polymap(matrix_t *v, matrix_t *m, poly_t *p)

{
    matrix_t *result, *tmp;
    int i;

    result = matalloc(v->fl,v->nor,v->noc);
    tmp = matdup(v);
    for (i = 0; i <= p->deg; ++i)
    {
	FEL f = p->buf[i];
	PTR x = result->d, y = tmp->d;
	long i;
	for (i = v->nor; i > 0; --i)
	{   zaddmulrow(x,y,f);
	    zadvance(&x,1);
	    zadvance(&y,1);
	}
	matmul(tmp,m);
    }
    matfree(tmp);
    return result;
}


/* ------------------------------------------------------------------
   make_kern() - Calculate a vector in the null-space of p(A). p(x)
   is assumed to be a factor of f1(x), i.e., p(x) must occur in the
   first cyclic subspaces generated in make_charpol().
   ------------------------------------------------------------------ */

void make_kern(node_t *n, poly_t *p)

{
    poly_t *cof, *f;
    matrix_t *seed;

    f = poldup(n->f1);
    cof = poldivmod(f,p);
    if (f->deg != -1) printf("ERROR 1");
    polfree(f);
    if (n->nsp != NULL) matfree(n->nsp);
    seed = matalloc(zfl,1,n->dim);
    zinsert(seed->d,CharPolSeed,F_ONE);
    n->nsp = polymap(seed,n->word,cof);
    matfree(seed);
    polfree(cof);
}


/* ------------------------------------------------------------------
   make_trkern() - Calculate a vector in the null-space of p(A^T). We
   consider only the first cyclic subspace, so make_trkern() may not
   find a vector.

   Returns 0 if a vector has been found, -1 else.
   ------------------------------------------------------------------ */

int make_trkern(node_t *n, poly_t *p)

{
    matrix_t *mt;
    poly_t *pt, *cof;
    int result = -1;

    mt = mattr(n->word);
    pt = charpolfactor(mt);
    cof = poldivmod(pt,p);
    if (pt->deg == -1)		/* Glueck gehabt! */
    {
	matrix_t *seed = matalloc(zfl,1,n->dim);
    	zinsert(seed->d,CharPolSeed,F_ONE);
    	if (n->nsptr != NULL) matfree(n->nsptr);
    	n->nsptr = polymap(seed,mt,cof);
	matfree(seed);
	result = 0;
    }
    polfree(cof);
    polfree(pt);
    matfree(mt);
    return result;
}


/* ------------------------------------------------------------------
   make_charpol() - Calculate the characteristic polynomial.
   ------------------------------------------------------------------ */

fpoly_t *make_charpol(n)
node_t *n;

{
    fpoly_t *cpol;
    int i, k;

    if (n->f2 != NULL) polfree(n->f2);
    n->f2 = NULL;
    if (n->f1 != NULL) polfree(n->f1);
    n->f1 = charpolfactor(n->word);
    cpol = factorization(n->f1);

    /* Rearrange the factors
       --------------------- */
    for (i = 0; i < cpol->len; ++i)
    {
    	for (k = i+1; k < cpol->len; ++k)
    	{
	    long val1 = cpol->e[i];
	    long val2 =  cpol->e[k];
	    if (val1 > val2)
	    {
	    	poly_t *p = cpol->p[i];
	    	long e = cpol->e[i];
	    	cpol->p[i] = cpol->p[k];
	    	cpol->e[i] = cpol->e[k];
	    	cpol->p[k] = p;
	    	cpol->e[k] = e;
	    }
    	}
    }
    if (MSG3) fpolprint("  f1(x)",cpol);
    return cpol;
}

void make_charpol2(n)	/* Make f2 */
node_t *n;

{
    poly_t *f;

    if (n->f2 != NULL) return;
    n->f2 = polalloc(n->f1->fl,0);
    while ((f = charpolfactor(NULL)) != NULL)
    {
	polmul(n->f2,f);
	polfree(f);
    }
    if (MSG3)
    { fpoly_t *x = factorization(n->f2);
      fpolprint("  f2(x)",x);
      fpolfree(x);
    }
}


static int tryit(n,pol,vfh)
node_t *n;
poly_t *pol;
long vfh;

{
    long nul;
    int result;
    int chkirred = 1;

    make_kern(n,pol);
    nul = n->nsp->nor;

    /* Try to split
       ------------ */
    if (n->subsp != NULL) matfree(n->subsp);
    result = spinup(n->nsp,ngen,n->gen,SPL_SEED_FIRST,&(n->subsp));
    if (result == SPL_SUCCESS) 		/* Module was split */
    {
	++stat_nssplit;
	set_insert(goodwords,n->wnum);
	splitnode(n,0);
	return 1;
    }
    MESSAGE(3,("Split failed\n"));

    if (n->f2 == NULL) make_charpol2(n);	/* Make f2 */
{
 poly_t *ggt = polgcd(pol,n->f2);
 long deg = ggt->deg;
 polfree(ggt);
 if (vfh != 1 || deg != 0)
 {
    chkirred = 0;
 }
}


    /* Not split, try dual split
       ------------------------- */
    MESSAGE(3,("Try dual split...\n"));
    if (make_trkern(n,pol) != 0)
    {
	MESSAGE(3,(" No seed vector found, dual split skipped\n"));
	return 0;
    }
    if (n->gentr == NULL)	/* Transpose generators */
    {   int i;
	n->gentr = NALLOC(matrix_t *,ngen);
	for (i = 0; i < ngen; ++i)
	    n->gentr[i] = mattr(n->gen[i]);
    }
    result = spinup(n->nsptr,ngen,n->gentr,SPL_SEED_FIRST,&(n->subsp));
    if (result == SPL_SUCCESS) 	/* Dual split */
    {
	++stat_dlsplit;
	splitnode(n,1);
	return 1;
    }

    /* Dual split failed*/
    MESSAGE(3,("Dual split failed\n"));
    if (!chkirred)
	return 0;
    newirred(n);
    set_insert(goodwords,n->wnum);
    ++stat_irred;
    return 1;
}


/* ------------------------------------------------------------------
   chop_word() - Chop a module using a given word. Returns 1 on
   success, 0 otherwise.
   ------------------------------------------------------------------ */

static int chop_word(node_t *n, long wn)

{
    long dlimit = opt_deglimit;		/* Limit on degree */
    int pi, done;
    fpoly_t *cpol;

    CharPolSeed = (CharPolSeed + 1) % n->dim + 1;
    makeword(n,wn);
    cpol = make_charpol(n);

    /* If c(x) is irreducible, then the module is irreducible
       ------------------------------------------------------ */
    if (cpol->p[0]->deg == n->dim)
    {
	MESSAGE(2,("c(x) is irreducible\n"));
	newirred(n);
	++stat_cpirred;
	set_insert(goodwords,wn);
	fpolfree(cpol);
	return 1;
    }

    /* Try all factors of c(x)
       ----------------------- */
    for (done = pi = 0; !done && pi < cpol->len; ++pi)
    {    
	if (MSG3) {
	    printf("Next factor: ");
	    polprint(NULL,cpol->p[pi]);
	    printf("^%ld\n",cpol->e[pi]);
	}
	if (dlimit > 0 && cpol->p[pi]->deg > dlimit)
	{
	   MESSAGE(3,("deg > %ld -- discarded\n",dlimit));
	   continue;
	}
	done = tryit(n,cpol->p[pi],cpol->e[pi]);
	MESSAGE(3,("done=%d\n",done));
    }
    fpolfree(cpol);
    return done;
}


/* ------------------------------------------------------------------
   chop() - Chop a module
   ------------------------------------------------------------------ */

static node_t *chop(n)
node_t *n;			/* Parent node */

{
    long w;

    MESSAGE(0,("Chop: Dim=%ld\n",n->dim));

    /* Check if the module is one-dimensional
       -------------------------------------- */
    if (n->dim == 1)
    {
	MESSAGE(1,("Dimension is one -- irreducible!\n"));
	newirred(n);
	++stat_irred;
	return n;
    }

    /* Try saved vectors
       ----------------- */
    if (chop_saved_vectors(n)) return n;

    /* Try words
       --------- */
    for (w = 0; ; )
    {
	long i;
	if ((int)-w < goodwords->len)
	    i = goodwords->buf[(int) -(w--)];
	else if (w <= 0)
	    w = i = 1;
	else
	    i = ++w;
	if (isbadword(i,n)) continue;
	MESSAGE(2,("Next word is %ld (=%s)\n",i,SymbolicName(n->wg,i)));
	if (chop_word(n,i))
	    break;
	/*set_insert(n->badwords,i);*/
    }
    return n;
}


/* ------------------------------------------------------------------
   extendbasis() - Find a vector in 'space' which is not in the
	(linear) span of 'basis'. 'basis' must be linearly
	independent.
   ------------------------------------------------------------------ */

static matrix_t *extendbasis(basis, space)
matrix_t *basis, *space;

{
    long i, j, piv;
    long dimb = basis->nor;
    long dims = space->nor;
    PTR tmp, x, y;
    FEL f;


    /* Concatenate basis and space
       --------------------------- */
    zsetlen(basis->noc);
    tmp = zalloc(dimb+dims);
    memcpy(tmp,basis->d,zsize(dimb));
    x = tmp;
    zadvance(&x,dimb);
    memcpy(x,space->d,zsize(dims));

    /* Clean with basis
       ---------------- */
    for (i = 1, x = tmp; i <= dimb; zadvance(&x,(long)1), ++i)
    {
	piv = zfindpiv(x,&f);
	if (piv == 0) FATAL("extendbasis(): zero vector in basis");
	y = x;
	for (j = i+1; j <= dimb+dims; ++j)
	{
	    zadvance(&y,(long)1);
	    zaddmulrow(y,x,zsub(F_ZERO,zdiv(zextract(y,piv),f)));
	}
    }


    /* Find the first non-zero row
       --------------------------- */
    x = tmp;
    zadvance(&x,dimb);
    for (j = 1; zfindpiv(x,&f) == 0; ++j, zadvance(&x,(long)1));
    free(tmp);
    if (j > dims)
    {
	FATAL("extendbasis() failed");
    }
    return matextract(space,j,(long)1);
}



/* ------------------------------------------------------------------
   checkspl() - Checks if a given representation's splitting field
   has degree [E:F] = dim(V), where V is a given subspace (usually
   the kernel of an algebra element).
   
   Returns 1 if [E:F]=dim(V) or zero otherwise.
   ------------------------------------------------------------------ */

#define MAXENDO 10	/* Max. dimension of endomorphism ring */


static int checkspl(gen,nsp)
matrix_t *gen[];
matrix_t *nsp;

{
    matrix_t *sb1, *sb2;	/* Standard bases */
    matrix_t *g1[MAXGEN];	/* Generators in standard basis sb1 */
    matrix_t *g2[MAXGEN];	/* Generators in standard basis sb2 */
    matrix_t *endo[MAXENDO];	/* Endomorphisms */
    int nendo = 0;		/* # of endomorphisms obtained so far */
    long dim = gen[0]->nor;
    int i, result;

    /* Take the first vector from nsp and change to standard basis.
       ------------------------------------------------------------ */
    sb1 = sbasis(nsp,ngen,gen);		/* Spin up to std basis */
    chbasis(sb1,ngen,gen,g1);

    sb2 = NULL;	/* Mark as unused */
    while (1)
    {  
	matrix_t *v2, *subsp;

	/* Spin up v1 under all endomorphisms found so far. If this
	   yields the whole null-space, we know that the endomorphism
	   ring has at least dimension dim(nsp).
	   ---------------------------------------------------------- */
	spinup(nsp,nendo,endo,SPL_SEED_FIRST,&subsp);
	if (subsp->nor == nsp->nor)
	{
	    matfree(subsp);
	    result = 1;	/* Successfull! */
	    break;
	}

	/* Take a vector which is not in span(v1)
	   and make again the standard basis.
	   --------------------------------------- */
	if (sb2 != NULL)	/* Clean up first */
	{
	    int j;
	    matfree(sb2);
	    for (j = 0; j < ngen; ++j) matfree(g2[j]);
	}
	v2 = extendbasis(subsp,nsp);	/* Get vector */
	matfree(subsp);
	sb2 = sbasis(v2,ngen,gen);	/* Spin up */
	matfree(v2);
	chbasis(sb2,ngen,gen,g2);

	/* Compare the two representations. If they are different,
	   we know that the splitting field degree must be smaller
	   than dim(nsp).
	   -------------------------------------------------------- */
	result = 1;
	for (i = 0; result && i < ngen; ++i)
	    if (memcmp(g1[i]->d,g2[i]->d,zsize(dim)))
	        result = 0;
	if (result == 0) break;	/* Not successfull */

	/* They are identical, i.e., we have found an endomorphism.
	   Put it into the list and try the next vector.
	   -------------------------------------------------------- */
	if (nendo >= MAXENDO) FATAL("Too many endomorphisms");
	endo[nendo] = matinv(sb2);
	matmul(endo[nendo],sb1);
	++nendo;
    }

    /* Clean up
       -------- */
    matfree(sb1); 
    for (i = 0; i < ngen; ++i) matfree(g1[i]);
    if (sb2 != NULL)
    {
	matfree(sb2);
        for (i = 0; i < ngen; ++i) matfree(g2[i]);
    }
    while (nendo > 0)
	matfree(endo[--nendo]);
    MESSAGE(3,("checkspl(): result = %d\n",result));
    return result;
}



/* ------------------------------------------------------------------
   findidword() - Find an identifying word for a module. This is a
   word in the generators with minimal nullity, i.e., the nullity
   equals the splitting field degree [E:F]. The word is stored in
   n->idword, the polynomial in n->idpol, and its null-space in
   n->nsp.
   ------------------------------------------------------------------ */

static int findidword(n)
node_t *n;

{
    long i;
    int k;
    long count = 0;
    static fpoly_t *cpol = NULL;

    /* Main loop: Try all words
       ------------------------ */
    MESSAGE(1,("Searching idword, dim=%ld\n",n->dim));
    for (i = 1; count <= MAXTRIES; ++i)
    {
	if (isbadword(i,n)) continue;
	n->wnum = i;

	/* Make the word and its characteristic polynomial
	   ----------------------------------------------- */
	MESSAGE(2,("  Word %ld  gcd=%ld\n",i,n->ggt));
	makeword(n,i);
	if (cpol != NULL) fpolfree(cpol);
	cpol = charpol(n->word);
	if (MSG3) fpolprint("  c(x)",cpol);
	for (k = 0; k < cpol->len; ++k)
	    n->ggt = gcd(n->ggt,cpol->e[k] * cpol->p[k]->deg);

	/* Try all factors with degree<=gcd
	   ------------------------------- */
	for (k = 0; k < cpol->len; ++k)
	{
	    if (cpol->p[k]->deg > n->ggt) continue;
	    ++count;
	    if (MSG3) polprint("  factor",cpol->p[k]);
	    insertword(n,cpol->p[k]);
	    n->ggt = gcd(n->ggt,n->nsp->nor);
	    if (checkspl(n->gen,n->nsp))
	    {
		MESSAGE(2,("  idword = %ld\n",i));
		n->idword = i;
		n->idpol = poldup(cpol->p[k]);
		return 0;
	    }
	}
    }

    /* Failed...
       --------- */
    FATAL("findidword() failed");
    return 0;
}



/* ------------------------------------------------------------------
   equiv() - Check if two modules are isomorphic. n2 is assumed to
   be in standard basis.
   ------------------------------------------------------------------ */

static int equiv(n1, n2)
node_t *n1, *n2;

{
    matrix_t *m, *b, *bi, *seed;
    int result = 0, j;
    size_t datasize;
    matrix_t *word;

    word = MakeWord(n1->wg,n2->idword);
    m = matinsert(word,n2->idpol);
    matfree(word);
    seed = nullspace__(m);
    if (seed->nor != n2->spl)
    {	matfree(seed);
	return 0;
    }

    datasize = zsize(n1->dim);
    b = sbasis(seed,ngen,n1->gen);
    matfree(seed);
    bi = matinv(b);
    for (j = 0; result == 0 && j < ngen; ++j)
    {
	matrix_t *g = matdup(b);
	matmul(g,n1->gen[j]);
	matmul(g,bi);
	if (memcmp(g->d,n2->gen[j]->d,datasize))
	    result = 1;
	matfree(g);
    }
    matfree(b);
    matfree(bi);
    return (result == 0);
}




/* ------------------------------------------------------------------
   newirred() - Check if a given irreducible module is already
	contained in the list of composition factors. If yes, return
	the index. If not, insert the new irreducible and return its
	index.
   ------------------------------------------------------------------ */


static void newirred(n)
node_t *n;

{
    int i, k;
    matrix_t *b;
    
    /* Check if the module is already in the list
       ------------------------------------------ */
    makefp(n->wg,n->fprint);
    for (i = 0; i < ncf && n->dim >= irred[i]->dim; ++i)
    {
	/* Compare dimensions and fingerprints
	   ----------------------------------- */
	if (n->dim != irred[i]->dim ||
	    memcmp(n->fprint,irred[i]->fprint,sizeof(n->fprint)))
	    continue;

	/* Check for equivalence
	   --------------------- */
	if (equiv(n,irred[i]))
	{
	    ++irred[i]->mult;
	    n->num = irred[i]->num;
    	    cfinfo[i].dim = irred[i]->dim;  /* cfname() needs this */
    	    cfinfo[i].num = irred[i]->num;
	    MESSAGE(0,("Irreducible (%s)\n",cfname(i)));
	    cleannode(n);
	    return;
	}
    }

    /* It's a new irreducible!
       ----------------------- */
    if (ncf >= MAXCF)
	FATAL("TOO MANY COMPOSITION FACTORS");
    for (k = ncf-1; k >= i; --k)
	irred[k+1] = irred[k];
    irred[i] = n;
    ++ncf;

    /* Set the fields
       -------------- */
    if (i == 0 || irred[i]->dim != irred[i-1]->dim)
	irred[i]->num = 0;
    else
	irred[i]->num = irred[i-1]->num + 1;
    irred[i]->mult = 1;
    cfinfo[i].dim = irred[i]->dim;  /* cfname() needs this */
    cfinfo[i].num = irred[i]->num;
    MESSAGE(0,("New irreducible (%s)\n",cfname(i)));

    /* Make idword and change to std basis 
       ----------------------------------- */
    if (n->word != NULL) matfree(n->word);
    n->word = NULL;
    if (n->nsp != NULL) matfree(n->nsp);
    n->nsp = NULL;
    if (n->idword == -1)
	findidword(n);
    else
    {
	matrix_t *m;
    	n->word = MakeWord(n->wg,n->idword);
    	m = matinsert(n->word,irred[i]->idpol);
    	n->nsp = nullspace__(m);
    }
    n->spl = n->nsp->nor;
    b = sbasis(n->nsp,ngen,n->gen);
    chbasis(b,ngen,n->gen,n->gen);
    matfree(b);

    /* Write out the generators
       ------------------------ */
    if (irred[i]->spl > 1)
	MESSAGE(1,("Splitting field has degree %ld\n",irred[i]->spl));
    for (k = 0; k < ngen; ++k)
    {
	char fn[200];
	sprintf(fn,"%s%s.%d",cfbasename,cfname(i),k+1);
   	matsave(irred[i]->gen[k],fn);
    }
    cleannode(n);
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

#define MALLOC_CKCHAIN	103	/* turn on chain checking	*/
#define MALLOC_FILLAREA	104	/* turn off area filling	*/
#define MALLOC_CKDATA	106	/* turn off/on data checking	*/
#define MALLOC_ERRFILE	112

union dbmalloptarg
{
	int	  i;
	char	* str;
} mo;

int  dbmallopt(int, union dbmalloptarg *);



int main(argc, argv)
int argc;
char *argv[];

{   
    int i;

#if  0 && defined(OS_ULTRIX)
mo.str = "dbmalloc.err";
dbmallopt(MALLOC_ERRFILE,&mo);
mo.i = 1;
dbmallopt(MALLOC_CKCHAIN,&mo);
dbmallopt(MALLOC_FILLAREA,&mo);
dbmallopt(MALLOC_CKDATA,&mo);
malloc_chain_check(1);
#endif

    mtxinit();
    mtxerraction = 2;	/* Print error messages */
    goodwords = set_alloc();
    opt_nullimit = zfl > 10 ? 2 : 3;

    /* Parse command line
       ------------------ */
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("Gs:g:n:d:")) != OPT_END)
    {
	switch (i)
	{
	    case 'G': opt_G = 1; msg_level = -100; break;
	    case 's':
		firstword = getint();
		if (firstword < 1 || *opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-s");
		break;
	    case 'g':
		ngen = getint();
		if (ngen < 1 || ngen > MAXGEN || *opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-g");
		break;
	    case 'd':
		opt_deglimit = getint();
		if (opt_deglimit < 0 || *opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-d");
		break;
	    case 'n':
		opt_nullimit = getint();
		if (opt_nullimit < 1 || *opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-n");
		break;
	}
    }

    if (opt_ind != argc-1) errexit(ERR_NARGS,"chop");
    setbasename(argv[argc-1]);
    MESSAGE(0,("*** CHOP MODULE ***\n\n"));
    init();
    chop(root);
    writeinfo(root);
    if (msg_level >= 0)
    {
        printf("\n");
        prtimes();
    }
    return 0;
}



