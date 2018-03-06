/* ========================== C MeatAxe =============================
   mksub.c - Calculate submodule lattice.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: mksub.c,v 1.2 1997/09/11 15:43:11 gap Exp $
 *
 * $Log: mksub.c,v $
 * Revision 1.2  1997/09/11 15:43:11  gap
 * New version 2.2.3. AH
 *
 * Revision 2.15  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.14  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.13  1994/06/25  12:13:39  mringe
 * Sortiere die Bloecke, sonst stimmt die Zuordnung
 * der Mountains zu den irred. nicht.
 *
 * Revision 2.12  1994/03/26  06:34:04  mringe
 * basename umbenannt wg. Namenskonflikt.
 *
 * Revision 2.11  1994/03/13  13:27:01  mringe
 * Maschinenunabhaengiges Format fuer Permutationen.
 *
 * Revision 2.10  1994/02/21  12:47:24  mringe
 * Debug-printfs raus.
 *
 * Revision 2.9  1994/02/19  15:18:49  mringe
 * Debug-printfs entfernt.
 *
 * Revision 2.8  1994/02/15  10:28:33  mringe
 * MSDOS_BCC entfernt.
 *
 * Revision 2.7  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.6  1994/02/12  04:11:38  mringe
 * Neuer bitstring_t.Typ.
 *
 * Revision 2.5  1993/12/08  11:33:02  mringe
 * Neue CPU time - Funktionen.
 *
 * Revision 2.4  1993/12/02  18:33:51  mringe
 * Ersetze bitstring_t durch bitstring_t *.
 *
 * Revision 2.3  1993/10/27  11:12:13  mringe
 * printf() formats.
 *
 * Revision 2.2  1993/10/22  16:08:19  mringe
 * Neues Numerierungsschema fuer irreduzible.
 *
 * Revision 2.2  1993/10/22  16:08:19  mringe
 * Neues Numerierungsschema fuer irreduzible.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.40  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.39  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.38  1993/10/02  16:05:36  mringe
 * Schreibe Bitstrinfs nach xxx.sub.
 *
 * Revision 1.37  1993/08/27  15:27:26  mringe
 * Option -T
 *
 * Revision 1.36  1993/08/10  14:29:19  mringe
 * Include string.h
 *
 * Revision 1.35  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.34  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.33  1993/07/28  13:34:49  mringe
 * Gap-Output: Fange Indizs mit 1 an.
 *
 * Revision 1.32  1993/07/28  08:24:05  mringe
 * Neues Format f"ur .lat-File.
 *
 * Revision 1.31  1993/07/26  08:43:37  mringe
 * Bug in extend behoben (MAXDOTL statt MAXCYCL)
 *
 * Revision 1.30  1993/07/23  13:46:27  mringe
 * OS-Symbole neu (SYS_xxx)
 *
 * Revision 1.29  1993/07/19  13:43:05  mringe
 * extend() verbessert - schneller!!!
 *
 * Revision 1.28  1993/07/19  09:27:38  mringe
 * Pruefe auf Fehler nach bs_read().
 *
 * Revision 1.27  1993/07/17  19:13:05  mringe
 * Aenderungen fuer Borland C.
 *
 * Revision 1.26  1993/07/16  13:27:00  mringe
 * Option -O zur Kontrolle des Output-Formats.
 *
 * Revision 1.1  1993/07/13  20:30:59  mringe
 * Initial revision
 *
 * Revision 1.25  1993/07/02  08:00:35  mringe
 * Berechne Dimensionen und maximale bei vorzeitigem
 * Abbruch (gibt sonst core dump, wenn MAXNSUB  erreicht wird).
 *
 * Revision 1.24  1993/07/02  00:58:01  mringe
 * Debug-Output rausgenommen.
 *
 * Revision 1.23  1993/07/01  23:56:51  mringe
 * Ueberlappende D.L. nicht zusammenwerfen --- leider
 * falsch gedacht!
 *
 * Revision 1.22  1993/07/01  22:39:25  mringe
 * Lese Dimensionen der Mountains von XXX.mnt
 *
 * Revision 1.21  1993/06/30  11:26:00  mringe
 * Code vereinfacht.
 *
 * Revision 1.20  1993/03/03  14:26:05  mringe
 * Prototypes ergaenzt.
 *
 * Revision 1.20  1993/03/03  14:26:05  mringe
 * Prototypes ergaenzt.
 *
 * Revision 1.19  1993/03/03  14:05:01  mringe
 * Option -b (blocks) implementiert.
 *
 * Revision 1.18  1993/03/03  10:08:35  mringe
 * help() eingebaut.
 *
 * Revision 1.17  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.16  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.15  1993/01/09  12:58:35  mringe
 * *** empty log message ***
 *
 * Revision 1.14  1993/01/06  21:12:56  mringe
 * ???
 *
 * Revision 1.13  1992/10/01  13:50:10  mringe
 * Header eingef"ugt.
 *
 * Revision 1.12  1992/09/07  10:33:26  mringe
 * Radikalreihe korrigiert (o geh"ort dazu)
 *
 * Revision 1.11  1992/09/04  10:07:38  mringe
 * Schreibe Anzahl der Teilmoduln ins .gra-File.
 *
 * Revision 1.10  1992/09/04  09:41:39  mringe
 * Schreibe .gra-File f"ur mkgraph & Diverse "Anderungen.
 *
 * Revision 1.9  1992/08/26  09:51:04  mringe
 * Schreibe Isomorphietypen der einfachen Faktoren raus.
 *
 * Revision 1.8  1992/08/15  12:35:41  mringe
 * Schreibe Liste der Mountains raus.
 *
 * Revision 1.7  1992/07/23  19:15:00  mringe
 * Schreibe Inzidenzen name X.lat.
 *
 * Revision 1.6  1992/07/22  08:22:56  mringe
 * Removed 'zzz.h'
 *
 * Revision 1.5  1992/07/22  07:10:30  mringe
 * Changed 'global.h' to 'lattice.h'
 *
 * Revision 1.4  1992/07/13  11:59:16  mringe
 * Sockel- und Radikalreihe
 *
 * Revision 1.3  1992/07/13  06:00:30  mringe
 * Sockelreihe
 *
 * Revision 1.2  1992/07/04  12:49:55  mringe
 * Benutze, da"s dlspan[i] abgeschl. unter Inzidenzen ist
 *
 * Revision 1.1  1992/05/25  17:43:40  mringe
 * Initial revision
 */

#include <string.h>
#include <stdlib.h>

#include "meataxe.h"
#include "lattice.h"
#include "files.h"



#define BITS (bnmount < 100)


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

#define O_MOUNTAINS		0x01
#define O_SUBMODULES		0x02
#define O_DOTTEDLINES		0x04
#define O_EXTFILES		0x08
#define O_RADICAL		0x10
#define O_SOCLE			0x20
#define O_INCIDENCES		0x40
#define O_ALL (O_MOUNTAINS|O_SUBMODULES|O_DOTTEDLINES|O_EXTFILES|\
	       O_RADICAL|O_SOCLE|O_INCIDENCES)


int opt_b = 0;			/* -b option (blocks) */
int opt_o = O_ALL;
int done[MAXCF];		/* */
int blnum;			/* Number of current block */
int blsize;			/* Block size */
int block[MAXCF];		/* Block members */

int firstm[MAXCF+1];		/* First mountain */
int firstdl[MAXCF+1];		/* First dotted line */

/* Data read from input files
   -------------------------- */
int xnmount = 0;		/* Number of mountains */
int xndotl = 0;			/* Number of dotted lines */
bitstring_t *xsubof[MAXCYCL];	/* Incidence matrix */
bitstring_t *xdotl[MAXDOTL];	/* Dotted lines */
long xmdim[MAXCYCL];		/* Mountain dimensions */

/* Data for current block
   ---------------------- */
int bnmount, bndotl;		/* As above, but for current block */
bitstring_t *bsubof[MAXCYCL];
bitstring_t *bsupof[MAXCYCL];	/* Transposed incidence matrix */
bitstring_t *bdotl[MAXDOTL];
bitstring_t *bdlspan[MAXDOTL];
long bmdim[MAXCYCL];

/* Data used during computation
   ---------------------------- */
bitstring_t *sub[MAXNSUB];		/* Submodules */
int nsub = 0;			/* Number of submodules */
int lastgen = 0;		/* First submodule of last generation */
int generation;			/* Generation count */
int nadd;			/* Number of calls to addtolist() */
int oldnsub;			/* */
bitstring_t *y;			/* Temporary bit string */
long *subdim;			/* Submodule dimensions */
char *israd;			/* Radical series */
char *issoc;			/* Socle series */
char *ismount;			/* Mountains */
int **max;			/* List of maximal submodules */


static char *helptext[] = {
"SYNTAX",
"    mksub [<Options>] <Name>",
"",
"OPTIONS",
"    -b       Find blocks",
"    -o <Fmt>",
"    -n <Fmt>",
"             Output format. <Fmt> is any combination of m (mountains),",
"             d (dotted lines), i (incidence matrix), e (.lat and .gra files),",
"             s (submodule list), r (radical series), and o (socle series).",
"             -o includes, -n excludes a component.",
"    -T <MaxTime>]    Set CPU time limit",
"",
"FILES",
"        <Name>.cfinfo  i/o    Composition factor data ",
"        <Name>.inc     i      Incidence matrix generated by MKINC",
"        <Name>.dot     i      Dotted lines generated by MKDOTL",
"        <Name>.mnt     i      Mountain dimensions",
"        <Name>.out     o      Submodule lattice",
"        <Name>.lat     o      Incidence matrix of the submodule",
"                              lattice (in GAP format)",
"        <Name>.gra     o      Description of the lattice in a format",
"                              suitable as input for other programs",
"",
"    If -b is used, output files are produces for each block, and a",
"    block number is appended to the file names (e.g., `psl27.out.1').",
"",
NULL};

static proginfo_t pinfo =
   { "mksub", "Find Submodules",
     "$Revision: 1.2 $", helptext };



/* -----------------------------------------------------------------
   init() - Read .cfinfo and .inc file
   ----------------------------------------------------------------- */

static void init()

{	int i;
	long l[2];
	char fn[40];
	FILE *f;

	readcfinfo();

	/* Read incidence matrix
	   --------------------- */
	f = SFOpen(strcat(strcpy(fn,cfbasename),".inc"),FM_READ);
	if (f == NULL) FATAL("CANNOT OPEN .inc FILE");
	printf("Reading %s: ",fn);
	fflush(stdout);
	zreadlong(f,l,1);
	xnmount = (int) l[0];
	printf("%d mountain%s\n",xnmount,xnmount == 1 ? "" : "s");
	fflush(stdout);
	if (xnmount > MAXCYCL) FATAL("TOO MANY MOUNTAINS");
	bs_setlen(xnmount);
	for (i = 0; i < xnmount; ++i)
	{
		if ((xsubof[i] = bs_read(f)) == NULL)
		    FATAL("ERROR READING INCIDENCE MATRIX");
	}
	fclose(f);

	/* Read dotted lines
	   ----------------- */
	f = SFOpen(strcat(strcpy(fn,cfbasename),".dot"),FM_READ);
	if (f == NULL) FATAL("CANNOT OPEN .dot FILE");
	printf("Reading %s: ",fn);
	fflush(stdout);
	zreadlong(f,l,1);
	xndotl = (int) l[0];
	printf("%d dotted line%s\n",xndotl,xndotl == 1 ? "" : "s");
	fflush(stdout);
	if (xndotl > MAXDOTL) FATAL("TOO MANY DOTTED LINES");
	for (i = 0; i < xndotl; ++i)
	{
	    if ((xdotl[i] = bs_read(f)) == NULL)
		FATAL("ERROR READING DOTTED LINES");
	}
	fclose(f);

	y = bs_alloc();

    /* Read dimensions
       --------------- */
    f = SFOpen(strcat(strcpy(fn,cfbasename),".mnt"),FM_READ|FM_TEXT);
    if (f == NULL) FATAL("ERROR OPENING .mnt FILE");
    printf("Reading %s\n",fn);
    fflush(stdout);
    for (i = 0; i < xnmount; ++i)
    {
	long mno, mdim;
	if (fscanf(f,"%ld%ld",&mno,&mdim) != 2 || mno != i || mdim < 1)
	    FATAL("ERROR IN .md FILE");
	xmdim[i] = mdim;
	while (getc(f) != '\n')	/* Skip class */
	    if (ferror(f) || feof(f))
	    	FATAL("ERROR IN .md FILE");
    }
    fclose(f);
}


/* -----------------------------------------------------------------
   init2()
   ----------------------------------------------------------------- */

static void init2()

{
    int i;

    /* Set firstm and firstdl
       ---------------------- */
    firstm[0] = 0;
    firstdl[0] = 0;
    for (i = 0; i < ncf; ++i)
    {
	firstm[i+1] = firstm[i] + cfinfo[i].nmount;
	firstdl[i+1] = firstdl[i] + cfinfo[i].ndotl;
    }

    /* Initialize done[]
       ----------------- */
    for (i = 0; i < ncf; ++i) done[i] = 0;
}


/* -----------------------------------------------------------------
   isotype() - Finde den Isomorphietyp eines gegebenen Mountains.
   ----------------------------------------------------------------- */

static int isotype(int mnt)

{
    int m;

    for (m = 0; (mnt -= cfinfo[block[m]].nmount) >= 0; ++m);
    return (block[m]);
}


/* -----------------------------------------------------------------
   findrsm() - Find the radical and socle series, and the
	mountains.
   ----------------------------------------------------------------- */

void findrsm()

{
    char *flag = (char *) malloc((size_t) nsub);
    bitstring_t *bs = bs_alloc();
    int i, k;


    /* Berechne maximale Teilmoduln und Dimensionen
       -------------------------------------------- */
    ismount = (char *) malloc((size_t) nsub);
    max = (int **) malloc((size_t) nsub * sizeof(int *));
    subdim = (long *) malloc((size_t) nsub * sizeof(long));

    for (i = 0; i < nsub; ++i)
    {
	int maxcount = 0, *lp;
        memset(flag,0,(size_t) nsub);

	/* Bestimme alle maximalen Teilmoduln
	   ---------------------------------- */
	for (k = i-1; k >= 0; --k)
	{
	    if (flag[k] != 0) continue;
	    if (bs_issub(sub[k],sub[i]))
	    {
		int l;
		flag[k] = 1;
		++maxcount;
		for (l = k-1; l >= 0; --l)
		    if (bs_issub(sub[l],sub[k])) flag[l] = 2;
	    }
	}

	/* Baue eine Liste mit den Nummern der maximalen Teilmoduln
	   und den Isomorphietypen der einfachen Faktoren auf.
	   -------------------------------------------------------- */
	lp = max[i] =
	 (int *) malloc((size_t)(2*maxcount+1)*sizeof(int));
	for (k = 0; k < i; ++k)
	    if (flag[k] == 1)
	    {
		int l;
		*lp++ = k;
		for (l = 0; !bs_test(sub[i],l)||bs_test(sub[k],l); ++l);
		*lp++ = isotype(l);
	    }
	*lp = -1;
	ismount[i] = (maxcount == 1);

	/* Berechne die Dimension. Da wir aufsteigend vorgehen,
	   ist die Dimension der Maximalen hier schon bekannt.
	   ---------------------------------------------------- */
	if (maxcount == 0)
	    subdim[i] = 0;
	else
	{
	    subdim[i] = subdim[max[i][0]] + cfinfo[max[i][1]].dim;
	}
    }

    /* Berechne die Radikalreihe
       -------------------------- */
    israd = (char *) malloc((size_t) nsub);
    memset(israd,0,(size_t)nsub);
    for (i = nsub-1; i > 0; )
    {	int *lp;

	bs_cpy(bs,sub[i]);
	for (lp = max[i]; *lp >= 0; lp += 2)
	    bs_and(bs,sub[*lp]);
	for (i = nsub-1; !bs_issub(sub[i],bs); --i);
	israd[i] = 1;
    }

    /* Berechne die Sockelreihe
       ------------------------ */
    issoc = (char *) malloc((size_t) nsub);
    memset(issoc,0,(size_t)nsub);
    for (i = 0; i < nsub-1; )
    {
	/* Find the simple submodules
	   -------------------------- */
	memset(flag,0,(size_t) nsub);
	for (k = i+1; k < nsub; ++k)
	{
	    if (flag[k] != 0) continue;
	    if (bs_issub(sub[i],sub[k]))
	    {
		int l;
		flag[k] = 1;
		for (l = k + 1; l < nsub; ++l)
		    if (bs_issub(sub[k],sub[l])) flag[l] = 2;
	    }
        }

        /* Calculate the sum of all simple submodules (=Socle)
	   --------------------------------------------------- */
	bs_cpy(bs,sub[i]);
	for (k = i; k < nsub; ++k)
	    if (flag[k] == 1) bs_or(bs,sub[k]);

	for (i = 0; !bs_issub(bs,sub[i]); ++i);
	issoc[i] = 1;
    }
}


/* -----------------------------------------------------------------
   extend() - Add one mountain to a given bitstring_t
   ----------------------------------------------------------------- */

static char dlflag[MAXDOTL] = {0};

static void extend(x, i, nextend)
bitstring_t *x;
int i;
int nextend;

{
    int k;
    int changed;

    bs_or(x,bsupof[i]);		/* Add mountain and its subspaces */
    if (nextend) bs_clear(x,i);	/* Add the radical only */

    /* Make closure
       ------------ */
    memset(dlflag,0,sizeof(dlflag));
    for (changed = 1; changed; )
    {
	changed = 0;
        for (k = 0; k < bndotl; ++k)
        {
	    if (!dlflag[k] && bs_match(x,bdotl[k]) >= 2)
	    {
	        bs_or(x,bdlspan[k]);
	        dlflag[k] = 1;
	        changed = 1;
	    }
	}
    }
}


/* -----------------------------------------------------------------
   writeresult() - Write output files
   ----------------------------------------------------------------- */

FILE *openout(name)
char *name;

{
    FILE *f;
    char fn[200];

    sprintf(fn,opt_b ? "%s%s.%d" : "%s%s",cfbasename,name,blnum);
    printf("Writing %s\n",fn);
    fflush(stdout);
    f = SFOpen(fn,FM_CREATE|FM_TEXT);
    if (f == NULL)
    {
	perror(fn);
	FATAL("CANNOT OPEN OUTPUT FILE");
    }
    return f;
}

#define NDIG(x) (x>99999?6:x>9999?5:x>999?4:x>99?3:x>9?2:1)

static int printbs(f, b)
FILE *f;
bitstring_t *b;

{
    int k, k1, k2, len, flag;

    if (BITS)
    {
	for (k = 0; k < bnmount;  ++k)
	    putc(bs_test(b,k) ? '+' : '.',f);
	len = bnmount;
    }
    else
    {
	len = k = 0;
	flag = 0;
	while (1)
	{
	    while (k < bnmount && !bs_test(b,k)) ++k;
	    if (k >= bnmount) break;
	    k1 = k;
	    while (k < bnmount && bs_test(b,k)) ++k;
	    k2 = k-1;
	    if (flag)
	    {
		putc(',',f);
		++len;
	    }
	    else
		flag = 1;
	    if (k2 > k1)
	    {
		fprintf(f,"%d-%d",k1,k2);
		len += NDIG(k1)+NDIG(k2)+1;
	    }
	    else
	    {
		fprintf(f,"%d",k1);
		len += NDIG(k1);
	    }
	}
    }
    return len;
}



void writeresult()

{
    FILE *f, *g;
    int i, k;
    char tmp[100];
    bitstring_t *b = bs_alloc();

    f = openout(".out");

    /* Write irreducibles
       ------------------ */
    fprintf(f,"Irreducibles:\n");
    fprintf(f,"    Type   Mult   SF   Mountains           Dotted lines\n");
    for (i = 0; i < blsize; ++i)
    {
	sprintf(tmp,"%s",cfname(block[i]));
	fprintf(f,"    %-7s%-7ld%-5ld",tmp,cfinfo[block[i]].mult,
	    cfinfo[block[i]].spl);
	sprintf(tmp,"%ld (%ld-%ld)",cfinfo[block[i]].nmount,
	    (long) firstm[block[i]],
	    cfinfo[block[i]].nmount+firstm[block[i]]-1);
	fprintf(f,"%-20s",tmp);
	if (cfinfo[block[i]].ndotl > 0)
	    sprintf(tmp,"%ld (%ld-%ld)",cfinfo[block[i]].ndotl,
	        (long) firstdl[block[i]],
	        cfinfo[block[i]].ndotl+firstdl[block[i]]-1);
	else
	    sprintf(tmp,"0");
	fprintf(f,"%-20s\n",tmp);
    }
    fprintf(f,"\n");

    /* Write mountains
       --------------- */
    if (opt_o & O_MOUNTAINS)
    {
	fprintf(f,"Mountains:\n");
	fprintf(f,"    No     Dim    Maximal Submountains\n");
	for (i = 0; i < bnmount; ++i)
	{
	    fprintf(f,"    %-7d%-7ld",i,bmdim[i]);
	    bs_cpy(b,bsupof[i]);
	    bs_clear(b,i);
	    for (k = 0; k < bnmount; ++k)
	    {
		if (!bs_test(b,k)) continue;
		bs_minus(b,bsupof[k]);
		bs_set(b,k);
	    }
	    for (k = 0; k < bnmount; ++k)
		if (bs_test(b,k)) fprintf(f,"%d ",k);
	    fprintf(f,"\n");
	}
	fprintf(f,"\n");
    }

    /* Write incidence matrix
       ---------------------- */
    if (opt_o & O_INCIDENCES)
    {
	printf("  Incidence matrix (%d by %d)\n",bnmount,bnmount);
	fflush(stdout);
	fprintf(f,"Incidence matrix:\n");
	for (i = 0; i < bnmount; ++i)
	{
	    fprintf(f,"    %3d: ",i);
	    printbs(f,bsupof[i]);
	    fprintf(f,"\n");
	}
	fprintf(f,"\n");
    }

    /* Write dotted lines
       ------------------ */
    if (opt_o & O_DOTTEDLINES)
    {
	printf("  Dotted lines (%d)\n",bndotl);
	fflush(stdout);
	fprintf(f,"Dotted lines:\n");
	for (i = 0; i < bndotl; ++i)
	{
		fprintf(f,"    ");
		printbs(f,bdotl[i]);
		fprintf(f,"\n");
	}
	fprintf(f,"\n");
    }

    /* Write submodules
       ---------------- */
    if (opt_o & O_SUBMODULES)
    {
	printf("  Submodules (%d)\n",nsub);
	fflush(stdout);
    	g = openout(".sub");
	fprintf(f,"Submodules:\n");
	fprintf(f,"    No    Dim  Flags  Ident                           Max\n");
	for (i = 0; i < nsub; ++i)
	{
	    int *lp;

	    fprintf(f,"    %-6d%-5ld",i,subdim[i]);
	    putc(ismount[i] ? 'M' : ' ',f);
	    putc(israd[i] ? 'R' : ' ',f);
	    putc(issoc[i] ? 'S' : ' ',f);
	    fprintf(f,"    ");
	    k = printbs(f,sub[i]);
	    for (; k < 30; ++k) putc(' ',f);
	    fprintf(f,"  ");
	    for (lp = max[i]; *lp >= 0; lp += 2)
	    {
		if (lp != max[i]) putc(',',f);
		fprintf(f,"%d",*lp);
	    }
	    fprintf(f,"\n");
	    bs_write(g,sub[i]);
	}
	fprintf(f,"\n");
	fclose(g);
    }

    /* Radikal- und Sockelreihe
       ------------------------ */
    if (opt_o & O_RADICAL)
    {
	static int mult[MAXCF];
	long rdim;
	int layer;
	bitstring_t
	    *rad = bs_alloc(),
	    *newrad = bs_alloc(),
	    *x = bs_alloc(),
	    *zero = bs_alloc();

	MESSAGE(0,("  Radical series\n"));
	fflush(stdout);
	fprintf(f,"Radical series:\n");
	for (i = 0; i < bnmount; ++i) bs_set(rad,i);
	bs_reset(zero);
	for (i = 0, rdim = 0; i < ncf; ++i)
	    rdim += cfinfo[i].dim * cfinfo[i].mult;
	for (layer = 1; bs_cmp(rad,zero); ++layer)
	{
	    bs_reset(x);
	    bs_reset(newrad);
	    MESSAGE(1,("Starting layer %d\n",layer));

	    /* Extend the zero module = x by all those mountains
	       which are contained in the radical and nextend y
	       ------------------------------------------------- */
	    for (i = 0; i < bnmount && bs_cmp(rad,x); ++i)
	    {
		if ( bs_test(rad,i) && !(bs_test(x,i)) )
		{
		    MESSAGE(2,("extend(%i)\n",i));
		    extend(x,i,0);
		    MESSAGE(2,("nextend(%i)\n",i));
		    extend(newrad,i,1);
		}
	    }

	    /* Find the irreducible factors in this layer
	       ------------------------------------------ */
	    memset(mult,0,sizeof(mult));
	    bs_cpy(x,newrad);
	    for (i = 0; i < bnmount && bs_cmp(rad,x); ++i)
	    {
	      if ( bs_test(rad,i) && !(bs_test(x,i)) )
	      {
		  extend(x,i,0);
		  k = isotype(i);
		  ++mult[k];
		  rdim -= cfinfo[k].dim;
	      }
	    }
	    fprintf(f,"    Layer %d: Dim=%-4ld  ",layer,rdim);
	    for (i = 0; i < ncf; ++i)
		for (; mult[i] > 0; --mult[i])
		    fprintf(f,"%s ",cfname(i));
	    fprintf(f,"\n");
	    bs_cpy(rad,newrad);
	}
	fprintf(f,"\n");
    }
    fclose(f);


    if ((opt_o & O_EXTFILES) && (opt_o & O_SUBMODULES))
    {

	/* Write the .lat file
	   ------------------- */
	f = openout(".lat");
#if 0
	fprintf(f,"MeatAxe.IncidenceMatrix := [\n");
	for (i = 0; i < nsub; ++i)
	{
	    fprintf(f,"[");
	    for (k = 0; k < nsub; ++k)
	    {
		fprintf(f,"%d",bs_issub(sub[i],sub[k]) ? 1 : 0);
		if (k < nsub-1)
		    fprintf(f,",");
	    }
	    fprintf(f,"]");
	    if (i < nsub-1)
		fprintf(f,",");
	    fprintf(f,"\n");
	}
	fprintf(f,"];\n");
#endif

	fprintf(f,"MeatAxe.Lattice := [\n");
	for (i = 0; i < nsub; ++i)
	{
	    int *lp = max[i];

	    fprintf(f,"[%ld,[",subdim[i]);
	    for (k = 0, lp = max[i]; *lp >= 0; lp += 2, ++k)
	    {
	    	fprintf(f,"[%d,%d]",lp[0]+1,lp[1]+1);
		if (lp[2] >= 0)
		{
		    fprintf(f,",");
		    if (k % 10 == 9) fprintf(f,"\n");
		}
	    }
	    if (i < nsub-1)
		fprintf(f,"]],\n");
	    else
		fprintf(f,"]]\n");
	}
	fprintf(f,"];\n");


	fclose(f);

	/* Write the .gra file
	   ------------------- */
	f = openout(".gra");
	fprintf(f,"%d\n",nsub);
	for (i = 0; i < nsub; ++i)
	{
	int *lp;

	putc(ismount[i] ? 'm' : '.',f);
	putc(israd[i] ? 'r' : '.',f);
	putc(issoc[i] ? 's' : '.',f);
	for (k = 0, lp = max[i]; *lp >= 0; lp += 2, ++k);
	fprintf(f," %2d",k);
	for (lp = max[i]; *lp >= 0; lp += 2)
	    fprintf(f," %d %d",lp[0],lp[1]);
	fprintf(f,"\n");
	}
	fclose(f);
    }
}



/* -----------------------------------------------------------------
   sameblock() - Test if irred. #i and #k belong to the same block
	(Actually, we check only if there is any mountain of
	irred. i contained in any mountain of irred. k and vice
	versa)
   ----------------------------------------------------------------- */

static int sameblock(i,k)
int i, k;

{
    int ii, kk;

    for (ii = firstm[i]; ii < firstm[i+1]; ++ii)
	for (kk = firstm[k]; kk < firstm[k+1]; ++kk)
	{
	    if (bs_test(xsubof[ii],kk)) return 1;
	    if (bs_test(xsubof[kk],ii)) return 1;
	}
    return 0;
}




/* -----------------------------------------------------------------
   nextblock() - Build next block; returns 0 if no more
   	constituents remain
   ----------------------------------------------------------------- */

static int nextblock()

{
    int i, k;

    ++blnum;
    for (i = 0; i < ncf && done[i]; ++i);
    if (i >= ncf) return 0;

    /* If -b was not used, build one single
       block containing all irreducibles.
       ------------------------------------ */
    if (!opt_b)
    {
	blsize = ncf;
	for (i = 0; i < ncf; ++i)
	{
	    block[i] = i;
	    done[i] = 1;
	}
	return 1;
    }

    /* Otherwise, make the next block
       ------------------------------ */
    MESSAGE(2,("Making next block (%d)\n",blnum));
    done[i] = 1;
    blsize = 1;
    block[0] = i;
    i = 0;
    while (i < blsize)
    {
	for (k = 0; k < ncf; ++k)
    	    if (!done[k] && sameblock(block[i],k))
    	    {
	        done[k] = 1;
	        block[blsize++] = k;
	    }
	++i;
    }

    /* Sort block
       ---------- */
    MESSAGE(2,("Sorting\n"));
    for (i = 0; i < blsize; ++i)
	for (k = i+1; k < blsize; ++k)
	{
	    if (block[i] > block[k])
	    {
		int tmp = block[i];
		block[i] = block[k];
		block[k] = tmp;
	    }
	}
    if (MSG0)
    {
        printf("\nBlock %d: ",blnum);
	for (i = 0; i < blsize; ++i)
	    printf(" %s%s",cfbasename,cfname(block[i]));
	printf("\n");
        fflush(stdout);
    }
    return 1;
}



/* -----------------------------------------------------------------
   sort()
   ----------------------------------------------------------------- */

void sort()

{
    int i, k;
    bitstring_t *x;

    printf("Sorting\n");
    fflush(stdout);
    for (i = 0; i < nsub; ++i)
    {
	for (k = i+1; k < nsub; ++k)
	{
	    if (bs_issub(sub[k],sub[i]))
	    {
		x = sub[i];
		sub[i] = sub[k];
		sub[k] = x;
	    }
	}
    }
}

/* -----------------------------------------------------------------
   addtolist() - Find out if a given submodule is new. If yes, add
	the new sub module to the list
   ----------------------------------------------------------------- */

static void addtolist(x)
bitstring_t *x;

{
    int i;

    ++nadd;
    for (i = nsub-1; i >= 0; --i)
    {
	if (!bs_cmp(sub[i],x))
	    return;
    }
    if (nsub >= MAXNSUB)
    {
	fprintf(stderr,"Too many submodules (> %ld)\n",(long)MAXNSUB);
	sort();
	findrsm();
	writeresult();	/* Write out what we have found so far */
	FATAL("TOO MANY SUBMODULES - PROGRAM ABORTED");
    }
    sub[nsub] = bs_alloc();
    bs_cpy(sub[nsub++],x);
}





/* -----------------------------------------------------------------
   initblock()
   ----------------------------------------------------------------- */

static void initblock()

{
    int i, k, ii, kk, row, col;

    /* Find out the number of mountains in this block
       ---------------------------------------------- */
    bnmount = 0;
    for (i = 0; i < blsize; ++i)
	bnmount += cfinfo[block[i]].nmount;
    bs_setlen(bnmount);

    /* Build the incidence matrix
       -------------------------- */
    printf("Building incidence matrix\n");
    fflush(stdout);
    for (i = 0; i < bnmount; ++i)
    {
	bsubof[i] = bs_alloc();
	bsupof[i] = bs_alloc();
    }
    row = 0;
    for (i = 0; i < blsize; ++i)
    {
        for (ii = firstm[block[i]]; ii < firstm[block[i]+1]; ++ii)
	{
	    col = 0;
	    for (k = 0; k < blsize; ++k)
	    {
		for (kk=firstm[block[k]]; kk<firstm[block[k]+1]; ++kk)
		{
		    if (bs_test(xsubof[ii],kk))
		    {
			bs_set(bsubof[row],col);
		        bs_set(bsupof[col],row);
		    }
		    ++col;
		}
	    }
	    bmdim[row] = xmdim[ii];
	    ++row;
	}
    }

    /* Build the dotted lines for one block
       ------------------------------------ */
    printf("Building dotted lines\n");
    fflush(stdout);
    bndotl = 0;
    for (i = 0; i < blsize; ++i)
    {
        for (ii = firstdl[block[i]]; ii < firstdl[block[i]+1]; ++ii)
	{
	    bdotl[bndotl] = bs_alloc();
	    bdlspan[bndotl] = bs_alloc();
	    col = 0;
	    for (k = 0; k < blsize; ++k)
	    {
		for (kk=firstm[block[k]]; kk<firstm[block[k]+1]; ++kk)
		{
		    if (bs_test(xdotl[ii],kk))
		    {
			bs_or(bdlspan[bndotl],bsupof[col]);
			bs_set(bdotl[bndotl],col);
		    }
		    ++col;
		}
	    }
	    ++bndotl;
	}
    }

    /* Initialize global variables
       --------------------------- */
    generation = 0;
    nsub = 0;		/* Number of submodules */
    lastgen = 0;	/* First submodule of previous generation */
    nadd = 0;
    oldnsub = 0;
    addtolist(bs_alloc());	/* Null module */
}



/* -----------------------------------------------------------------
   cleanupblock() - Clean up after each block
   ----------------------------------------------------------------- */

static void cleanupblock()

{
    int i;

    for (i = 0; i < bnmount; ++i)
    {
	free(bsubof[i]);
	free(bsupof[i]);
    }
    for (i = 0; i < bndotl; ++i)
    {
	free(bdotl[i]);
    }
    if (opt_o & O_SUBMODULES)
    {
        for (i = 0; i < nsub; ++i)
        {
	    free(sub[i]);
	    free(max[i]);
        }
        free(ismount);
        free(israd);
        free(issoc);
        free(subdim);
        free(max);
    }
}



/* -----------------------------------------------------------------
   nextgen() - Make next generation (modules generated by n+1
	mountains)
   ----------------------------------------------------------------- */

void nextgen()

{
    int i, k, oldnsub = nsub;
    bitstring_t *x;

    x = bs_alloc();
    for (i = lastgen; i < oldnsub; ++i)
    {
	for (k = 0; k < bnmount; ++k)
	{
	   bs_cpy(x,sub[i]);
	   if (!bs_test(x,k))
	   {	extend(x,k,0);
		addtolist(x);
	   }
	}
    }
    lastgen = oldnsub;
    ++generation;
}


/* -----------------------------------------------------------------
   main()
   ----------------------------------------------------------------- */

int main(argc, argv)
int argc;
char *argv[];

{
    int opt_o_found = 0;
    int i;


    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("bo:n:")) != OPT_END)
    {
	switch (i)
	{
	    case 'b':
		opt_b = 1;
		break;
	    case 'o':
		if (!opt_o_found)
		{
		    opt_o = 0;
		    opt_o_found = 1;
		}
		for (; *opt_text_ptr != 0; ++opt_text_ptr)
		switch (*opt_text_ptr)
		    {
		    case 'm': opt_o |= O_MOUNTAINS; break;
		    case 's': opt_o |= O_SUBMODULES; break;
		    case 'd': opt_o |= O_DOTTEDLINES; break;
		    case 'e': opt_o |= O_EXTFILES; break;
		    case 'r': opt_o |= O_RADICAL; break;
		    case 'o': opt_o |= O_SOCLE; break;
		    case 'i': opt_o |= O_INCIDENCES; break;
		    default: errexit(ERR_OPTION,"-o");
		    }
		break;
	    case 'n':
		opt_o_found = 1;
		for (; *opt_text_ptr != 0; ++opt_text_ptr)
		switch (*opt_text_ptr)
		    {
		    case 'm': opt_o &= ~O_MOUNTAINS; break;
		    case 's': opt_o &= ~O_SUBMODULES; break;
		    case 'd': opt_o &= ~O_DOTTEDLINES; break;
		    case 'e': opt_o &= ~O_EXTFILES; break;
		    case 'r': opt_o &= ~O_RADICAL; break;
		    case 'o': opt_o &= ~O_SOCLE; break;
		    case 'i': opt_o &= ~O_INCIDENCES; break;
		    default: errexit(ERR_OPTION,"-u");
		    }
		break;
	}
    }

    if (opt_ind != argc-1)
	errexit(ERR_NARGS,"mksub");
    printf("\n*** CALCULATE ALL SUBMODULES ***\n\n");
    setbasename(argv[opt_ind]);
    init();
    init2();

    while (nextblock())
    {
	initblock();

	if (opt_o & O_SUBMODULES)
	{
	    do
	    {
		printf("Generation %d: ",generation);
		printf("%d tr%s, %d new submodule%s\n",
		    nadd,nadd == 1 ? "y":"ies",
		    nsub-oldnsub,nsub-oldnsub == 1 ? "" : "s");
		fflush(stdout);
		nadd = 0;
		oldnsub = nsub;
		nextgen();
	    }
	    while (oldnsub != nsub);
	    sort();
	    findrsm();
	}
	else
	    printf("Submodules not calculated\n");

	writeresult();
	printf("\n");
	cleanupblock();
    }
    prtimes();
    return 0;
}
	
	
