/* ========================== C MeatAxe =============================
   zpr.c - Print a matrix or permutaion.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zpr.c,v 1.2 1997/09/11 15:44:15 gap Exp $
 *
 * $Log: zpr.c,v $
 * Revision 1.2  1997/09/11 15:44:15  gap
 * New version 2.2.3. AH
 *
 * Revision 2.24  1995/03/01  10:14:40  mringe
 * Bug in primat().
 *
 * Revision 2.23  1995/02/09  13:59:17  mringe
 * Ganzzahlige Matrizen.
 *
 * Revision 2.22  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.21  1994/09/18  14:00:32  mringe
 * Polynome
 *
 * Revision 2.20  1994/09/17  16:33:37  mringe
 * Header fuer Polynome.
 *
 * Revision 2.19  1994/09/14  07:23:54  mringe
 * Benutze zchar um auf Primkoerper zu testen.
 *
 * Revision 2.18  1994/08/28  01:29:24  mringe
 * *** empty log message ***
 *
 * Revision 2.17  1994/08/27  12:48:15  mringe
 * Neues Dateiformat fuer Permutationen.
 *
 * Revision 2.16  1994/08/26  09:08:57  mringe
 * Neuer Fileheader fuer Matrizen.
 *
 * Revision 2.15  1994/08/23  14:31:48  mringe
 * Zwei kleine Bugs behoben (Gap format).
 *
 * Revision 2.14  1994/08/22  19:28:44  mringe
 * Print polynomials
 *
 * Revision 2.13  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.12  1994/07/22  09:24:54  mringe
 * Bug beseitigt, -r geht wieder.
 *
 * Revision 2.11  1994/05/18  09:25:05  mringe
 * stdio.h entfernt.
 *
 * Revision 2.10  1994/03/30  08:26:06  mringe
 * Gap-Permutationen korrigiert fuer p=1.
 *
 * Revision 2.9  1994/03/13  12:37:35  mringe
 * Maschinenunabhaengiges Dateiformat (Perm)
 *
 * Revision 2.8  1994/02/15  14:09:19  mringe
 * zsetlen() erweitert fuer Permutationen.
 *
 * Revision 2.7  1994/02/15  13:28:01  mringe
 * Benutze zseek() statt fseek().
 *
 * Revision 2.6  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.5  1993/12/08  12:00:41  mringe
 * malloc() prototype.
 *
 * Revision 2.4  1993/12/08  11:48:50  mringe
 * Compiler warnings.
 *
 * Revision 2.3  1993/10/26  10:47:35  mringe
 * Compiler Warnings.
 *
 * Revision 2.2  1993/10/21  21:57:35  mringe
 * Permutationen.
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
 * Revision 1.16  1993/10/06  06:09:00  mringe
 * Option -s.
 *
 * Revision 1.15  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.14  1993/08/10  13:46:18  mringe
 * Ueberfluessige klammern beim Perm.-Output weglassen.
 *
 * Revision 1.13  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.12  1993/07/19  17:12:36  mringe
 * Option -G: Benutze MeatAxe-Record.
 *
 * Revision 1.11  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.11  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.10  1993/06/28  18:13:00  mringe
 * GAP Perm.utationen: Benutze Filenamen fuer die Elemente.
 *
 * Revision 1.9  1993/03/04  16:24:27  mringe
 * GAP-Output verbessert (maximal 80 Zeichen/Zeile).
 *
 * Revision 1.8  1993/03/04  15:53:38  mringe
 * GAP-Output von Permutationen beschleunigt.
 * Benutze utils-Library fuer help().
 *
 * Revision 1.7  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.6  1992/07/22  14:06:38  mringe
 * Use zgen
 *
 * Revision 1.5  1992/07/15  09:25:55  mringe
 * Some minor changes.
 *
 * Revision 1.4  1992/05/26  18:48:37  mringe
 * help() eingebaut.
 *
 * Revision 1.3  1992/05/26  18:07:05  mringe
 * Ausgabe von Permutationen im GAP-Format: Zeilenl"ange korrigiert.

 * Revision 1.2  1992/05/26  17:59:05  mringe
 * Bei Permutationen im GAP-Format keine Fixpunkte ausgeben,
 * da GAP diese am Anfang einer Permutation nicht einlesen
 * kann.
 *
 * Revision 1.1  1992/05/26  17:41:49  mringe
 * Initial revision
 *
 */

#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include "meataxe.h"

/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static long HdrPos;
static long hdr[3];
static FILE *dest = NULL;	/* Output file */
static char *inpname = NULL;
static FILE *inpfile = NULL;
static long first = 1;		/* Output range (option -r) */
static long last = 10000000;

static char *helptext[] = {
"SYNTAX",
"    zpr [-Gs] [-r <Range>] [<Binfile> [<Textfile>]]",
"",
"OPTIONS",
"    -G   GAP output",
"    -r   Print a range of rows (or permutations). <Range> may be a",
"         single number, or of the form <First>-<Last>.",
"    -s   Print summary only\n",
"",
"FILES",
"    <Binfile>   i  A matrix or permutation in binary format",
"    <Textfile>  i  The output in text format (default: stdout)",
"",
NULL};

static proginfo_t pinfo =
   { "zpr", "Print Permutations Or Matrices", "$Revision: 1.2 $",
      helptext };



/* ------------------------------------------------------------------
   err() - Print error message and exit.
   ------------------------------------------------------------------ */

static void err(c)
int c;

{   fprintf(stderr,"ZPR ERROR - ");
    switch (c)
    {
	case 'f':
	    fprintf(stderr,"FILE I/O ERROR\n");
	    break;
	case 't':
	    fprintf(stderr,"CANNOT PRINT THIS FILE TYPE\n");
	    break;
    }
    exit(EXIT_ERR);
}


/* ------------------------------------------------------------------
   PrString(), PrLong() - Prettty printer
   ------------------------------------------------------------------ */

static void PrString(char *c)
{
    static int pos = 0;
    int l = strlen(c);

    if (l + pos >= 78)
    {
	fprintf(dest,"\n");
	pos = 0;
    }
    fprintf(dest,"%s",c);
    for (; *c != 0; ++c)
	if (*c == '\n') pos = 0; else ++pos;
}


static void PrLong(long l)

{
    char tmp[20];
    sprintf(tmp,"%ld",l);
    PrString(tmp);
}


/* ------------------------------------------------------------------
   prmatrix() - Print a matrix in standard format.
   ------------------------------------------------------------------ */

static void prmatrix()

{
    PTR m1;
    FEL f1;
    long loop1, j1;
    int md, mx, iv;

    last = hdr[1];
    first = 1;
    zsetfield(hdr[0]);
    zsetlen(hdr[2]);
    m1 = zalloc((long)1);

    if (hdr[0] < 10) {md = 1; mx = 80;}
    else if (hdr[0] < 100) {md = 3; mx = 25;}
    else if (hdr[0] < 1000) {md = 4; mx = 20;}
    else if (hdr[0] < 10000) {md = 5; mx = 15;}
    else {md = 6; mx = 12;}

    fprintf(dest,"matrix field=%ld rows=%ld cols=%ld\n",hdr[0],hdr[1],hdr[2]);
    for (loop1 = 1; loop1 <= last; ++loop1)
    {	if (zreadvec(inpfile,m1,1) != 1) errexit(-1,inpname);
	if (loop1 < first) continue;
	iv = 1;
	for (j1 = 1; j1 <= hdr[2]; ++j1)
	{
	    f1 = zextract(m1,j1);
	    switch (md)
	    {	case 1: fprintf(dest,"%1ld",zftoi(f1)); break;
		case 2: fprintf(dest,"%2ld",zftoi(f1)); break;
		case 3: fprintf(dest,"%3ld",zftoi(f1)); break;
		case 4: fprintf(dest,"%4ld",zftoi(f1)); break;
		case 5: fprintf(dest,"%5ld",zftoi(f1)); break;
		case 6: fprintf(dest,"%6ld",zftoi(f1)); break;
	    }
	    if (iv++ >= mx)
	    {	fprintf(dest,"\n");
		iv = 1;
	    }
	}
	if (iv > 1)
		fprintf(dest,"\n");
    }
}


/* ------------------------------------------------------------------
   prgapmat() - Print a matrix in GAP format.
   ------------------------------------------------------------------ */

static void prgapmat()

{   PTR m1;
    FEL f1;
    FEL gen;
    long loop1, j1;
    int cnt, isprimefield;


    zsetfield(hdr[0]);
    zsetlen(hdr[2]);
    gen = zgen;		/* Generator */
    isprimefield = (zchar == zfl);

    m1 = zalloc((long)1);
    PrString("MeatAxe.Matrix := [\n");
    for (loop1 = 1; loop1 <= hdr[1]; ++loop1)
    {	
	if (zreadvec(inpfile,m1,1) != 1) errexit(-1,inpname);
	cnt = 0;
	fprintf(dest,"[");
	for (j1 = 1; j1 <= hdr[2]; ++j1)
	{   if (cnt > 75)
	    {	fprintf(dest,"\n ");
		cnt = 0;
	    }
	    f1 = zextract(m1,j1);
	    if (isprimefield)
	    {   FEL f2=F_ZERO;
	   	long k=0;
	    	while (f2 != f1)
	    	{   f2 = zadd(f2,gen);
		    ++k;
		}
		fprintf(dest,"%ld",k);
		cnt += k>9999?5:k>999?4:k>99?3:k>9?2:1;
	    }
	    else
	    {   if (f1 == F_ZERO)
	    	{   fprintf(dest,"0*Z(%ld)",hdr[0]);
		    cnt += 5;
		    cnt += hdr[0]>9999?5:hdr[0]>999?4:hdr[0]>99?3:hdr[0]>9?2:1;
		}
		else
		{   FEL f2 = gen;
		    long k = 1;
		    while (f2 != f1)
		    {   f2 = zmul(f2,gen);
			++k;
		    }
		    fprintf(dest,"Z(%ld)^%ld",hdr[0],k);
		    cnt += 4;
		    cnt += hdr[0]>9999?5:hdr[0]>999?4:hdr[0]>99?3:hdr[0]>9?2:1;
		    cnt += k>9999?5:k>999?4:k>99?3:k>9?2:1;
		}
	    }
	    if (j1 < hdr[2])
	    {	fprintf(dest,",");
		++cnt;
	    }
	}
	fprintf(dest,"]");
	if (loop1 < hdr[1])
	    fprintf(dest,",");
	fprintf(dest,"\n");
    }
    fprintf(dest,"]");
    if (isprimefield)
	fprintf(dest,"*Z(%ld)",hdr[0]);
    fprintf(dest,";\n");
}


/* ------------------------------------------------------------------
   prgapperm() - Print a permutation in GAP format.
   ------------------------------------------------------------------ */

#define SIZE(i) ((i)<9?2 : (i)<99?3 : (i)<999?4 : (i)<9999?5 : \
	(i)<99999?6 : (i)>999999?7 : (i)<9999999?8 : 9)


static void prgapperm()

{
    long i, k, *p, pos;
    long cycle;
    long *perm;
    int empty;

    perm = (long *) malloc(sizeof(long) * hdr[1]);
    if (perm == NULL) errexit(ERR_NOMEM,"zpr: prgapperm()");

    PrString("MeatAxe.Perms := [\n");
    for (pos = 1; pos <= hdr[2]; ++pos)
    {
	if (zreadlong(inpfile,perm,hdr[1]) != hdr[1])
	    errexit(ERR_FILEREAD,inpname);
	PrString("    ");
	p = (long *) perm;
	cycle = 0;
	empty = 1;
	while (1)
	{
	    /* Find next cycle
	       --------------- */
	    while (cycle < hdr[1] && p[cycle] == 0) ++cycle;
	    if (cycle >= hdr[1]) break;	/* Done */
	    i = cycle;

	     /* Check if it is a fixed point (we don't
	        print fixed points, GAP doesn't like them)
	        ----------------------------------------- */
	    if (p[i] == i+1)
	    {
		p[i] = 0;	/* Mark it as done */
		continue;
	    }

	    /* Print cycle
	       ----------- */
	    empty = 0;
	    PrString("(");
	    PrLong(i+1);

	    while (1)
	    {
		k = i;
		i = p[i] - 1;
		p[k] = 0;
		if (p[i] == 0) break;
	    	PrString(",");
	    	PrLong(i+1);
	    }
	    PrString(")");
	}
	if (empty) PrString("()");
	if (pos < hdr[2]) PrString(",");
	PrString("\n");
    }
    PrString("];\n");
}


/* ------------------------------------------------------------------
   prgap() - Print a matrix or permutation in GAP format.
   ------------------------------------------------------------------ */

static void prgap()

{
    if (hdr[0] == -1)
	prgapperm();
    else if (hdr[0] >= 2)
	prgapmat();
    else
	err('t');
}



/* ------------------------------------------------------------------
   prpol() - Print a polynomial
   ------------------------------------------------------------------ */

static void prpol()

{
    int i;
    poly_t *p;
    
    SFSeek(inpfile,HdrPos);
    if ((p = polread(inpfile)) == NULL)
	errexit(-1,inpname);
    

    fprintf(dest,"polynomial field=%ld degree=%ld\n",hdr[1],hdr[2]);
    for (i = 0; i <= p->deg; ++i)
    {
	PrLong(zftoi(p->buf[i]));
	PrString(" ");
    }
    PrString("\n");
}


/* ------------------------------------------------------------------
   prperm() - Print a permutation in standard format.
   ------------------------------------------------------------------ */

static void prperm()

{
    long f1, i, count;
    long *perm;

    perm = NALLOC(long,hdr[1]);
    if (perm == NULL) errexit(ERR_NOMEM,"zpr: prperm()");

    for (count = 1; count <= hdr[2]; ++count)
    {
        fprintf(dest,"permutation degree=%ld\n",hdr[1]);
	if (zreadlong(inpfile,perm,hdr[1]) != hdr[1])
	    errexit(ERR_FILEREAD,inpname);
	for (i = 0; i < hdr[1]; ++i)
	{
	    f1 = perm[i];
	    PrLong(f1);
	    PrString(" ");
	}
	PrString("\n");
    }
}

/* ------------------------------------------------------------------
   primat() - Print an integer matrix
   ------------------------------------------------------------------ */

static void primat()

{
    long k, i;
    long *row;

    row = NALLOC(long,hdr[2]);
    if (row == NULL) errexit(ERR_NOMEM,"zpr: primat()");

    fprintf(dest,"integer matrix rows=%ld cols=%ld\n",hdr[1],
	hdr[2]);
    for (i = 0; i < hdr[1]; ++i)
    {
	if (zreadlong(inpfile,row,hdr[2]) != hdr[2])
    	    errexit(ERR_FILEREAD,inpname);
    	for (k = 0; k < hdr[2]; ++k)
    	{
    	    PrLong(row[k]);
    	    PrString(" ");
    	}
    	PrString("\n");
    }
}



/* ------------------------------------------------------------------
   prmtx() - Print a matrix or permutation in standard format.
   ------------------------------------------------------------------ */

static void prmtx()

{	if (hdr[0] > 1)
		prmatrix();
	else if (hdr[0] == -1)
		prperm();
	else if (hdr[0] == -2)
		prpol();
	else if (hdr[0] == -T_IMAT)
		primat();
	else
		err('t');
}


/* ------------------------------------------------------------------
   setrange() - Set output range
   ------------------------------------------------------------------ */

static void setrange(c)
char *c;

{	
    if ((first = getint()) == GETINT_ERR) errexit(ERR_OPTION,"-r");
    if (*opt_text_ptr == 0)
    {	last = first;
    }
    else if (*opt_text_ptr == '-')
    {   ++opt_text_ptr;
        if ((last = getint()) == GETINT_ERR) 
    	    errexit(ERR_OPTION,"-r");
    }
    if (*opt_text_ptr != 0 || first < 1 || last < first)
       	errexit(ERR_OPTION,"-r");
}

/* ------------------------------------------------------------------
   printsummary()
   ------------------------------------------------------------------ */

static void printsummary()

{
    printf("%s: ",inpname);
    if (hdr[0] == -1)
    {
	printf("%ld Permutation%s of degree %ld\n",hdr[2],
		hdr[2] == 1 ? "" : "s", hdr[1]);
	SFSeek(inpfile,ftell(inpfile) + hdr[1]*hdr[2]*4);
    }
    else if (hdr[0] >= 2)
    {
	PTR x;
	long l;
	printf("%ld x %ld matrix over GF(%ld)\n",hdr[1],hdr[2],hdr[0]);
	zsetfield(hdr[0]);
	zsetlen(hdr[2]);
	x = zalloc(1);
	for (l = hdr[1]; l > 0; --l)
	    zreadvec(inpfile,x,1);
	free(x);
    }
    else if (hdr[0] == -2)
    {
	printf("Polynomial of degree %ld over GF(%ld)\n",hdr[2],hdr[1]);
	zsetfield(hdr[1]);
	zsetlen(hdr[2]+1);
	SFSeek(inpfile,ftell(inpfile) + zsize(1));
    }
    else if (hdr[0] == -T_IMAT)
    {
	printf("%ld x %ld integer matrix\n",hdr[1],hdr[2]);
	SFSeek(inpfile,ftell(inpfile)+hdr[1]*hdr[2]*4);
    }
    else
	printf("Unknown file format.\n");
}

/* ------------------------------------------------------------------
   ReadHeader() - Read file header. Returns 1 on success, 0 on error
   or EOF.
   ------------------------------------------------------------------ */

static int ReadHeader(void)

{
    HdrPos = ftell(inpfile);
    if (feof(inpfile)) return 0;
    if (zreadlong(inpfile,hdr,3) != 3) return 0;
    /* Check the header */
    if (hdr[0] > 65536 || hdr[0] < -20 || hdr[1] < 0 || hdr[2] < 0)
    {
	fprintf(stderr,"Bad header in data file (%ld %ld %ld)\n",
	hdr[0],hdr[1],hdr[2]);
	errexit(ERR_FILEFMT,inpname);
    }
    return 1;
}

/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{
    enum {MTX,GAP,SUMMARY} mode = MTX;
    int i;

    dest = stdout;	/* Output file */
    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("Gsr:")) != OPT_END)
    {
	switch (i)
	{
	    case 'G': mode = GAP; msg_level = -100; break;
	    case 's': mode = SUMMARY; break;
	    case 'r': setrange(opt_text); break;
	}
    }
    switch (argc - opt_ind)
    {  
    	case 2:
	    if ((dest = SFOpen(argv[opt_ind+1],FM_CREATE|FM_TEXT))
							== NULL)
	    {	perror(argv[opt_ind+1]);
		err('f');
	    }
	    /* no break here! */
        case 1:
	    inpname = argv[opt_ind];
	    if ((inpfile = SFOpen(inpname,FM_READ)) == NULL)
		errexit(-1,inpname);
	    break;
    	default:
	    errexit(ERR_BADUSAGE,pinfo.name);
    }

    while (ReadHeader())
    {
    	switch (mode)
    	{
	    case MTX: prmtx(); break;
	    case GAP: prgap(); break;
	    case SUMMARY: printsummary(); break;
    	}
    }
    fclose(inpfile);
    fclose(dest);
    return (EXIT_OK);
}


