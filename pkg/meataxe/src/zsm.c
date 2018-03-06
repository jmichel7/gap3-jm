/* ========================== C MeatAxe =============================
   zsm.c - Small matrices package.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zsm.c,v 1.2 1997/09/11 15:44:27 gap Exp $
 *
 * $Log: zsm.c,v $
 * Revision 1.2  1997/09/11 15:44:27  gap
 * New version 2.2.3. AH
 *
 * Revision 2.18  1994/11/30  10:28:58  mringe
 * ANSI-C, Neue Random-Funktionen.
 *
 * Revision 2.17  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.16  1994/11/25  13:33:08  mringe
 * Neue Random...-Funktionen
 *
 * Revision 2.15  1994/08/22  19:28:44  mringe
 * Bug behoben, mwXid jetzt korrekt.
 *
 * Revision 2.14  1994/08/22  08:41:49  mringe
 * `Make word' erweitert: mw5id erzeugt 5*Einheitsmatrix.
 *
 * Revision 2.13  1994/07/29  07:43:08  mringe
 * Neuer Wortgenerator.
 *
 * Revision 2.12  1994/07/23  17:03:09  mringe
 * nullspace()__
 *
 * Revision 2.12  1994/07/23  17:03:09  mringe
 * nullspace()__
 *
 * Revision 2.11  1994/06/27  14:11:05  mringe
 * Erlaube -g 1.
 *
 * Revision 2.10  1994/06/19  16:14:24  mringe
 * words.h eliminiert.
 *
 * Revision 2.9  1994/05/18  10:11:38  mringe
 * initbasis(): Reihenfolge der Argumente geaendert.
 *
 * Revision 2.8  1994/05/10  09:42:14  mringe
 * Pruefe vor dem Einlesen der Erzeuger, ob Argumente vorhanden sind.
 *
 * Revision 2.7  1994/05/04  11:42:50  mringe
 * Neuer Wortalgorithmus und neue Numerierung.
 *
 * Revision 2.6  1994/03/12  11:51:13  mringe
 * mw: Teste, ob Filename fuer das Wort angegeben wurde.
 *
 * Revision 2.5  1994/02/15  10:28:33  mringe
 * MSDOS_BCC entfernt.
 *
 * Revision 2.4  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.3  1993/12/08  11:48:50  mringe
 * Compiler warnings.
 *
 * Revision 2.2  1993/12/08  10:30:36  mringe
 * Option -T (Time limit).
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.24  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.23  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.22  1993/10/05  17:33:26  mringe
 * permsave() und permload().
 *
 * Revision 1.21  1993/10/02  16:23:02  mringe
 * matread() und matwrite() in matload() bzw. matsave() umbenannt.
 *
 * Revision 1.20  1993/09/20  20:36:33  mringe
 * -Q, -V funktionieren jetzt, -G fuer Random Orders.
 *
 * Revision 1.19  1993/08/10  14:51:42  mringe
 * Include string.h
 *
 * Revision 1.18  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.17  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.16  1993/07/28  13:34:49  mringe
 * *** empty log message ***
 *
 * Revision 1.15  1993/07/23  13:46:27  mringe
 * OS-Symbole neu (SYS_xxx)
 *
 * Revision 1.14  1993/07/19  16:55:20  mringe
 * Optionen -G, -Q, -V.
 *
 * Revision 1.13  1993/07/17  19:13:05  mringe
 * Aenderungen fuer Borland C.
 *
 * Revision 1.12  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.11  1993/05/21  08:59:18  mringe
 * Lese Erzeuger von <Name>.1, <Name>.2, ..., wenn -g benutzt wird.
 *
 * Revision 1.10  1993/05/12  11:02:53  mringe
 * Option -g im Helptext aufgefuehrt.
 *
 * Revision 1.9  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.8  1993/02/15  13:20:50  mringe
 * Auf N Erzeuger umgestellt. Random words benutzt jetzt random.c
 *
 * Revision 1.7  1993/02/13  11:34:45  mringe
 * Auf N Erzeuger umgestellt.
 * nameof() nach words.c ausgelagert.
 *
 * Revision 1.6  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.5  1992/07/22  08:13:02  mringe
 * Function prototypes
 *
 * Revision 1.4  1992/05/31  12:34:09  mringe
 * Diverse kleine Bugs entfernt.
 * Syntax der Dokumentation angepasst.
 *
 * Revision 1.3  1992/05/30  13:39:01  mringe
 * Alten Code entfernt.
 *
 * Revision 1.2  1992/05/29  07:42:33  mringe
 * Power f"ur Permutationen
 *
 * Revision 1.1  1992/05/26  07:38:23  mringe
 * Initial revision
 */

#include <ctype.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>

#include "meataxe.h"
#include "lattice.h"	/* For MAXGEN */



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

char **args;
int arg_ind, nargs;
int ngen=2;
int opt_g = 0;
int opt_G = 0;
matrix_t *gen[MAXGEN];
perm_t *genp[MAXGEN];
int gentype;
static char *helptext[] = {
"SYNTAX",
"    zsm [-GQV] [-T <#Secs>] [-g <#Gen>] <Command> <Files>",
"",
"    zsm pwr<Num> <Mat|Perm> <Result>",
"      Raise a matrix or permutation to the <Num>-th power",
"",
"    zsm ro[<Num>] <Mat|Perm1> <Mat|Perm2>",
"      Calculate <Num> random orders. Default is <Num>=120.",
"",
"    zsm fp <Mat1> <Mat2>",
"       Calculate the standard fingerprint (6 Nullities).",
"",
"    zsm fp[<Start>-]<End> <Mat1> <Mat2>",
"      Calculate nullities from word number <Start> (default = 1) up",
"      to word number <End>.",
"",
"    zsm mw<Num> <Mat1> <Mat2> <Word> [<Nullspace>]",
"      Write word number <Num> to <Word> and its nullspace, if",
"      desired, to <Nullspace>. `mw<Num>id' makes <Num> * Id.",
"",
"OPTIONS",
"    -G    GAP output.",
"    -Q    Quiet, no messages.",
"    -V    Verbose, more messages.",
"    -T    Set CPU time limit.",
NULL
};

static proginfo_t pinfo =
   { "zsm", "Small Matrices Package", "$Revision: 1.2 $", helptext };




/* ------------------------------------------------------------------
   error() - Print error message and exit
   ------------------------------------------------------------------ */

static void error(char *fmt, ...)

{
    va_list vl;

    va_start(vl,fmt);
    fprintf(stderr,"zsm: ");
    vfprintf(stderr,fmt,vl);
    va_end(vl);
    fprintf(stderr,"\n");
    exit(EXIT_ERR);
}

/* ------------------------------------------------------------------
   filetype() - Get field parameter from file header
   ------------------------------------------------------------------ */

static long filetype(char *fn)

{
    long fl;
    FILE *f;

    f = SFOpen(fn,FM_READ);
    if (f == NULL) errexit(-1,fn);
    if (zreadlong(f,&fl,1) != 1) errexit(-1,fn);
    fclose(f);
    return fl;
}


/* ------------------------------------------------------------------
   readgens() - Read matrices or permutations
   ------------------------------------------------------------------ */

static void readgens(void)

{
    int i;
    char name[100];
    char name0[100];
    char buf[200];
    long ft;

    for (i = 0; i < ngen; ++i)
    {
	if (opt_g)
	{
	    if (arg_ind >= nargs)
		error("missing generator name");
	    sprintf(name,"%s.%d",args[arg_ind],i+1);
	}
	else
	{
	    if (arg_ind+i >= nargs)
		error("missing generator name");
	    strcpy(name,args[arg_ind+i]);
	}
	ft = filetype(name);
	if (i == 0)
	{   gentype = (ft >= 2) ? T_MATRIX : T_PERM;
	    strcpy(name0,name);
	}
	switch (gentype)
	{
	    case T_MATRIX:
		if (ft < 2) errexit(ERR_NOTMATRIX,name);
		gen[i] = matload(name);
	        if (gen[i]->nor != gen[i]->noc)
	    	    errexit(ERR_NOTSQUARE,name);
	       	if (i > 0)
	    	{
		    if (gen[i]->fl != gen[0]->fl ||
		        gen[i]->nor != gen[0]->nor)
			{
			    sprintf(buf,"%s and %s",name0,name);
		            errexit(ERR_INCOMPAT,buf);
			}
	    	}
		break;
	    case T_PERM:
	    	if (filetype(name) != -1) errexit(ERR_NOTPERM,name);
	    	genp[i] = permload(name);
	    	if (i > 0)
	    	{
		    if (genp[i]->deg != genp[0]->deg)
			{
			    sprintf(buf,"%s and %s",name0,name);
		            errexit(ERR_INCOMPAT,buf);
			}
	    	}
		break;
	}
    }
    if (opt_g)
	++arg_ind;
    else
	arg_ind += 2;
}


/* ------------------------------------------------------------------
   Makeword() - Calculate a word and its null-space
   ------------------------------------------------------------------ */

static void Makeword(long l)

{
    matrix_t *nsp, *word;
    wgdata_t *b;

    readgens();
    if (arg_ind >= nargs) error("Missing output file name");
    if (l >= 0)
    {
    	b = WGInit(ngen,gen);
    	word = MakeWord(b,l);
    	MESSAGE(0,("WORD %ld (%s)",l,SymbolicName(b,l)));
    }
    else
    {
	long dim = gen[0]->nor, i;
	FEL f = zitof(-l);
	PTR x;
	word = matalloc(gen[0]->fl,dim,dim);
	for (x = word->d, i = 1; i <= dim; ++i, zadvance(&x,1))
	    zinsert(x,i,f);
    	MESSAGE(0,("WORD %ldid (%ld*Id)",-l,-l));
    }
    matsave(word,args[arg_ind++]);
    if (arg_ind < nargs)
    {
	nsp = nullspace__(word);
	matsave(nsp,args[arg_ind++]);
	MESSAGE(0,(" HAS NULLITY %ld",nsp->nor));
    }
    MESSAGE(0,("\n"));
}


/* ------------------------------------------------------------------
   fingerprint() - Calculate nullities
   ------------------------------------------------------------------ */

static void fingerprint(long l1,long l2)

{
    wgdata_t *b;
    long w;
    int count = 0;

    readgens();
    b = WGInit(ngen,gen);

    if (opt_G)
	printf("MeatAxe.Fingerprint := [");
    else
        printf("NUMBER NULLITY WORD\n");
    for (w = l1; w <= l2; ++w)
    {	matrix_t *word = MakeWord(b,w);
        if (opt_G)
	{
	    if (w > l1) printf(",");
	    if (++count % 5 == 0) printf("\n  ");
	    printf("[%ld,%ld]",w,nullity__(word));
	}
    	else
	    printf("%6ld %7ld %s\n",w,nullity__(word),
		SymbolicName(b,w));
    }
    if (opt_G)
	printf("];\n");
}


/* ------------------------------------------------------------------
   std_fingerprint() - Standard fingerprint (6 nullities)
   ------------------------------------------------------------------ */

static void std_fingerprint(void)

{
    wgdata_t *b;
    int i;
    long fp[MAXFP];

    readgens();
    b = WGInit(ngen,gen);
    makefp(b,fp);
    if (opt_G)
	printf("MeatAxe.Fingerprint := [ ");
    else
        printf("FINGERPRINT: ");
    for (i = 0; i < MAXFP; ++i)
    {
	if (opt_G)
	{
	    if (i > 0) printf(",");
	    printf("%ld",fp[i]);
	}
	else
	    printf("%3ld ",fp[i]);
    }

    if (opt_G)
	printf("];\n");
    else
        printf("\n");
}



/* ------------------------------------------------------------------
   random_orders() - Calculate random orders
   ------------------------------------------------------------------ */

static void random_orders(long n)

{
    matrix_t *m = NULL;
    perm_t *p = NULL;
    long count;

    readgens();
    if (gentype == T_MATRIX)
	m = matdup(gen[0]);
    else
	p = permdup(genp[0]);

    if (opt_G)
	printf("MeatAxe.RandomOrders := [");

    for (count = 0; count < n; ++count)
    {
	long o;
	if (gentype == T_MATRIX)
	    o = matorder(m);
	else
	    o = permorder(p);
	if (count % 12 == 0) printf("\n    ");
	printf("%4ld",o);
    	if (opt_G && count < n-1) printf(",");

	if (gentype == T_MATRIX)
	    matmul(m,gen[RandInt(ngen)]);
	else
	    permmul(p,genp[RandInt(ngen)]);
    }
    if (opt_G)
	printf("];\n");
    else
	printf("\n");
}



/* ------------------------------------------------------------------
   power() - Power of a matrix
   ------------------------------------------------------------------ */

static void power(long n)

{
    matrix_t *x, *y;
    perm_t *p, *q;

    if (filetype(args[0]) == -1)
    {	p = permload(args[0]);
	q = permpower(p,n);
	permsave(q,args[1]);
    }
    else
    {	x = matload(args[0]);
   	if (x->nor != x->noc)
	    errexit(ERR_NOTSQUARE,args[0]);
	y = matpower(x,n);
	matsave(y,args[1]);
    }
    arg_ind += 2;
}


/* ------------------------------------------------------------------
   mygetint() - Convert string to integer
   match() - Compare string with pattern
   ------------------------------------------------------------------ */

char *cp;

static int mygetint(long *l)

{	if (!isdigit(*cp))
		return 1;
	for (*l = 0; isdigit(*cp); ++cp)
		*l = 10 * *l + (*cp - '0');
	return 0;
}


static int match(char *pattern)

{	char *tmp = cp;

	while (*pattern != 0)
	{	switch (*pattern)
		{	case '$':
				if (*tmp != 0) return 1;
				break;
			case '#':
				if (!isdigit(*tmp)) return 1;
				while (isdigit(*tmp)) ++tmp;
				break;
			default:
				if (*tmp != *pattern) return 1;
				++tmp;
				break;
		}
		++pattern;
	}
	cp = tmp;
	return 0;
}




/* ------------------------------------------------------------------
   init() - Process command line options and arguments
   ------------------------------------------------------------------ */

static void init(int argc, char *argv[])

{   int i;

    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("Gg:")) != OPT_END)
    {
	switch (i)
	{
	    case 'G': opt_G = 1; msg_level = -100; break;
	    case 'g':
		opt_g = 1;
		ngen = getint();
		if (ngen < 1 || ngen > MAXGEN || *opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-g");
		break;
	}
    }

    /* Skip command
       ------------ */
    if (opt_ind >= argc) errexit(ERR_NARGS,"zsm");
    cp = argv[opt_ind++];


    /* Set file names
       -------------- */
    args = argv + opt_ind;
    nargs = argc - opt_ind;
    arg_ind = 0;
}




/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char *argv[])

{

    init(argc,argv);

    if (!match("pwr"))	/* pwr - Power */
    {	long l;
	if (mygetint(&l) || match("$") || nargs != 2)
	    error("Usage: `zsm pwr<Num> <Input> <Result>'");
	power(l);
    }
    else if (!match("ro"))	/* ro - Random orders */
    {	long l = 120;
	if ((match("$") && (mygetint(&l) || match("$"))))
	    error("Usage: `zsm ro<Num> <Gen1> <Gen2>'");
	random_orders(l);
    }
    else if (!match("fp"))	/* fp - Fingerprint */
    {
	if (!match("$"))
	    std_fingerprint();
	else
	{   long l1, l2;
	    if (mygetint(&l1))
	    	error("Usage: fp[<Num>[-<Num>]]");
	    if (!match("$"))
		fingerprint(1,l1);
	    else
	    {	if (match("-") || mygetint(&l2) || match("$"))
	    	    error("Usage: fp[<Num>[-<Num>]]");
		fingerprint(l1,l2);
	    }
	}
    }
    else if (!match("mw"))	/* mw - Make word */
    {	long l;

	if (mygetint(&l)) error("Usage: mw<Number>");
	if (!match("id$")) l = -l;
	else if (match("$")) error("Usage: mw<Number>");
	Makeword(l);
    }
    else errexit(ERR_BADUSAGE,"zsm");
    if (arg_ind != nargs)
	error("Extra arguments ignored: %s ...",args[arg_ind]);
    return EXIT_OK;
}


