/* ========================== C MeatAxe =============================
   zsp.c - Spin up and split.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zsp.c,v 1.2 1997/09/11 15:44:29 gap Exp $
 *
 * $Log: zsp.c,v $
 * Revision 1.2  1997/09/11 15:44:29  gap
 * New version 2.2.3. AH
 *
 * Revision 2.16  1995/06/16  14:34:02  mringe
 * Bug bei Permutationen behoben.
 *
 * Revision 2.15  1995/03/15  14:55:06  mringe
 * Option -t eingebaut.
 *
 * Revision 2.14  1994/11/30  10:28:58  mringe
 * ANSI-C, Neue Random-Funktionen.
 *
 * Revision 2.13  1994/11/15  10:23:52  mringe
 * Messages erweitert.
 * Bug behoben: Erkennt jetzt sofort, wenn dim >= nor wird.
 *
 * Revision 2.12  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.11  1994/07/22  09:24:54  mringe
 * Erlaube -g 1.
 *
 * Revision 2.10  1994/06/14  08:22:13  mringe
 * Bug behoben: Split mit quotient=0.
 *
 * Revision 2.9  1994/05/18  12:21:31  mringe
 * Fehler beim Testen der Seedpace-Dimension behoben.
 *
 * Revision 2.8  1994/05/16  08:28:00  mringe
 * Benutze zcleanrow(), zquotXX(), zmkechelon().
 * Optino -ss setzt jetzt NICHT mehr l.u. seed-Vektoren voraus.
 *
 * Revision 2.7  1994/05/11  08:44:27  mringe
 * split in dosplit umbenannt [Konflikt mit split()]
 *
 * Revision 2.6  1994/03/13  13:27:01  mringe
 * Maschinenunabhaengiges Format fuer Permutationen.
 *
 * Revision 2.5  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.4  1993/12/08  11:33:02  mringe
 * Neue CPU time - Funktionen.
 *
 * Revision 2.3  1993/10/21  21:57:35  mringe
 * Permutationen.
 *
 * Revision 2.2  1993/10/21  08:37:53  mringe
 * Compiler warnings.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.39  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.38  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.37  1993/08/24  12:50:16  mringe
 * Kleine Korrektur am helptext.
 *
 * Revision 1.36  1993/08/24  12:48:34  mringe
 * Option -T
 *
 * Revision 1.35  1993/08/24  08:40:43  mringe
 * Bilde Filenamen mit '.', wie CHOP. MAXDIM entfernt.
 *
 * Revision 1.34  1993/08/10  14:51:42  mringe
 * Include string.h
 *
 * Revision 1.33  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.32  1993/08/05  16:14:53  mringe
 * Prototype von err() korrigiert.
 *
 * Revision 1.31  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.30  1993/07/28  13:34:49  mringe
 * *** empty log message ***
 *
 * Revision 1.29  1993/07/23  13:46:27  mringe
 * OS-Symbole neu (SYS_xxx)
 *
 * Revision 1.28  1993/07/20  06:18:06  mringe
 * GAP-Output: benutze MeatAxe-Record.
 *
 * Revision 1.27  1993/07/17  19:13:05  mringe
 * Aenderungen fuer Borland C.
 *
 * Revision 1.26  1993/07/17  08:38:14  mringe
 * msg_t durch message_t ersetzt
 *
 * Revision 1.25  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.24  1993/07/08  10:09:42  mringe
 * Belege nur so viel Speicher, wie mit -d angegeben.
 *
 * Revision 1.23  1993/07/02  10:49:48  mringe
 * Deklaration von err() fuer Nicht-ANSI-Compiler
 *
 * Revision 1.22  1993/04/07  08:21:56  mringe
 * Verbesserte Fehlermeldung bei inkompatiblen Files
 *
 * Revision 1.21  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.20  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.19  1993/01/21  10:59:00  mringe
 * fflush(stdout) bei verbose output.
 *
 * Revision 1.18  1993/01/11  16:08:37  hiss
 * Pruefe, ob seed vector und Permutationen kompatibel sind.
 *
 * Revision 1.17  1993/01/08  12:57:59  hiss
 * Berechne Operation auf dem Teilraum bei Permutationen.
 *
 * Revision 1.16  1993/01/08  12:31:46  hiss
 * Option -v (verbose).
 *
 * Revision 1.15  1993/01/06  20:59:57  mringe
 * getopt in zgetopt() umbenannt.
 *
 * Revision 1.14  1992/11/04  09:45:31  mringe
 * Behandlung von Permutationen verbessert.
 *
 * Revision 1.13  1992/11/04  09:11:47  mringe
 * Erlaube Permutationen als Erzeuger.
 *
 * Revision 1.12  1992/10/02  10:05:53  mringe
 * Benutze pseed.c zum Generieren der Vektoren.
 *
 * Revision 1.11  1992/09/27  08:58:11  mringe
 * help() verbessert.
 *
 * Revision 1.10  1992/09/27  08:50:57  mringe
 * help() ge"andert.
 *
 * Revision 1.9  1992/08/14  10:52:03  mringe
 * Bug bei -ss korrigiert.
 *
 * Revision 1.8  1992/07/15  13:32:01  mringe
 * Prototypes mit _PLi, benutze args.c
 *
 * Revision 1.7  1992/07/15  13:02:09  mringe
 * Optionen -sg<Num> und -sr<Num>
 *
 * Revision 1.6  1992/07/15  09:25:55  mringe
 * Some minor changes.
 *
 * Revision 1.5  1992/06/29  16:55:06  mringe
 * Hilfstext verk"urzt.
 * Fehler bei der Auswertung der -s Option behoben.
 *
 * Revision 1.4  1992/06/28  10:40:34  mringe
 * Syntax vereinfacht.
 *
 * Revision 1.3  1992/06/28  08:59:12  mringe
 * Neue Version.
 *
 * interactive() eingebaut
 *
 * Revision 1.1  1992/05/26  18:33:33  mringe
 * Initial revision
 *
 */

#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include "meataxe.h"






/* ------------------------------------------------------------------
   Constants
   ------------------------------------------------------------------ */

#define EXIT_SPLIT	1  /* Return value if space was split */
#define EXIT_NOSPLIT	0  /* Return value if space was not split */

#define MAXGEN 20	/* Maximal number of generators */



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

enum {DONTKNOW,MATRIX,PERMUTATION}
	gentype = DONTKNOW;	/* Type of generators */
enum {GENERATE,READ,READALL}
	smode = GENERATE;	/* Seed mode */
int opt_G = 0;			/* GAP output */
int dosplit = 1;		/* Split after spin-up */
int output = 1;			/* Enable file ouput */
long seedvecno = 0;		/* Selected vector (0 = all) */
long iseed = 0;			/* Current vector */
long maxdim = -1;		/* Limit on susbspace dimension */
long maxtry = -1;		/* How often do we try to apply gens */
int nmat = 2;			/* Number of matrices to split */
int ngen = 2;			/* Number of generators */
char *genname[MAXGEN] = {"G1","G2"};
char *subname[MAXGEN] = {"P4","P5"};
char *quotname[MAXGEN] = {"P6","P7"};
char *seedname = "G3";
char *vname = "V";
FILE *seedfile;

PTR	gen[MAXGEN];		/* Matrices to split */
PTR	m3, 			/* Subspace */
	tmp1, tmp2;		/* Two extra rows */

long	nor,			/* Degree of representation */
	dim,			/* Dimension of invariant subspace */
	ssdim,			/* Dimension of seed vector space */
	fl;			/* Field order */

long *piv;		/* Pivot table for the subspace */
char *ins;		/* Used to find the insignificant columns */


/* ------------------------------------------------------------------
   help()
   ------------------------------------------------------------------ */

static char *helptext[] = {
"SYNTAX",
"    zsp [<Options>] [<G1> <G2> <Seed> [{<S1> <S2> <Q1> <Q2> | <Spc>}]]",
"    zsp [<Options>] -g <N>[.<M>] [<Gen> [<Seed> [{<Sub> <Quot> | <Spc>}]]]",
"",
"OPTIONS",
"    -s<Mode>      Select seed for spin-up:",
"                    -sg[<Num>]   Generate seed vectors (Default)",
"                    -sr[<Num>]   Read seed vectors one by one",
"                    -ss          Read a whole space and make its closure",
"                  <Num>, if specified, selects a single vector",
"    -n            Spin up and write out invariant subspace, but do not split",
"    -o            No output to files, messages only",
"    -V            Verbose",
"    -Q            Quiet, no messages",
"    -d <Dim>      Set limit on dimension of subspace",
"    -t <Num>      Assume subspace is invariant after <Num> tries",
"    -T <Maxtime>  Set limit on CPU time (in seconds)",
"    -g <N>[.<M>]  Read <M> matrices and treat the first <N> as generators",
"                  (default: <N>=<M>). Read matrices from <Gen>.i (i=1..<M>),",
"                  and write output to <Sub>.i and <Quot>.i.",
"",
"FILES",
"    <G1,2>   i  Generators: square matrices of equal size",
"    <Seed>   i  Seed vectors: a matrix with NOC(Seed)=NOC(G1)",
"    <S1,2>   o  Action on subspace",
"    <Q1,2>   o  Action on quotient",
"    <Spc>    o  Invariant subspace (with -n)",
NULL};

static proginfo_t pinfo =
   { "zsp", "Spin Up And Split", "$Revision: 1.2 $", helptext };



/* ------------------------------------------------------------------
   err() - Terminate program with error message
   ------------------------------------------------------------------ */

void err(int code)

{
    fprintf(stderr,"ZSP ERROR - ");
    switch (code)
    {
	case 'b':
	    fprintf(stderr,"INTERNAL ERROR\n");
	    break;
	case '0':
	    fprintf(stderr,"CANNOT SPIN UP ZERO VECTOR\n");
	    break;
	case 'v':
	    fprintf(stderr,"CANNOT MAKE THIS SEED VECTOR\n");
	    break;
	case 'S':
	    fprintf(stderr,
		  "SEED SPACE IS GREATER THAN DIMENSION LIMIT\n");
	    break;
	}
	exit(EXIT_ERR);
}



/* ------------------------------------------------------------------
   parseargs()
   ------------------------------------------------------------------ */

static void setname(char **name,int i,char *template)

{
    char *c;

    c = malloc(strlen(template) + 10);
    sprintf(c,"%s.%d",template,i+1);
    name[i] = c;
}



void parseargs(int argc, char *argv[])

{   int k, j;
    int opt_g = 0;


    /* Parse command line
       ------------------ */
    initargs(argc, argv, &pinfo);
    while ((k = zgetopt("ionvGt:d:g:s:")) != OPT_END)
    {   switch (k)
	{ 
	    case 'o':
		output = 0;
		break;
	    case 'G':
		msg_level = -100;
		opt_G = 1;
		break;
	    case 'n':
		dosplit = 0;
		break;
	    case 't':
		maxtry = getint();
		if (maxtry <= 0 || *opt_text_ptr != 0)
		    ErrorExit("Bad argument after -t");
		break;
	    case 'd':
		maxdim = getint();
		if (maxdim <= 0 || *opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-d");
		break;
	    case 'g':
		opt_g = 1;
		nmat = (int) getint();
		if (*opt_text_ptr == '.')
		{	++opt_text_ptr;
			ngen = (int) getint();
		}
		else
		    ngen = nmat;
		if (ngen > nmat || ngen < 1 || nmat > MAXGEN ||
				*opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-g");
		break;
	    case 's':	/* Seed mode */
		switch (opt_text[0])
		{   case 'g': smode = GENERATE; ++opt_text_ptr; break;
		    case 'r': smode = READ; ++opt_text_ptr; break;
		    case 's': smode = READALL; ++opt_text_ptr; break;
		    default: break; /* Allow -s<Num> for -sg<Num> */
		}
		if (smode == GENERATE || smode == READ)
		{   if (*opt_text_ptr != 0)
		    {	seedvecno = getint();
			if (seedvecno <= 0) errexit(ERR_OPTION,"-s");
		    }
		}
		if (*opt_text_ptr != 0) errexit(ERR_OPTION,"-s");
		break;
	}
    }

    /* Set file names
       -------------- */
    k = opt_ind;	/* Index of first file arg */
    if (!opt_g)		/* Old Syntax */
    {	switch (argc - k)
	{   case 0:
		break;
	    case 7:
		quotname[0] = argv[k+5];
		quotname[1] = argv[k+6];
	    case 5:
		subname[0] = argv[k+3];
		subname[1] = argv[k+4];
		if (!dosplit) errexit(ERR_BADUSAGE,"zsp");
	    case 4:
		if (!dosplit)
			vname = argv[k+3];
	    case 3:
		genname[0] = argv[k+0];
		genname[1] = argv[k+1];
		seedname = argv[k+2];
		break;
	    default:
		errexit(ERR_BADUSAGE,"zsp");
	}
    }
    else		/* The -g option was used */
    {
    	for (j = 0; j < nmat; ++j)
	{
	    setname(genname,j,"G");
	    setname(subname,j,"S");
	    setname(quotname,j,"Q");
	}
	seedname = "G0";
	switch (argc - k)
	{   case 4:
		if (!dosplit) errexit(ERR_BADUSAGE,"zsp");
		for (j = 0; j < nmat; ++j)
		    setname(quotname,j,argv[k+3]);
	    case 3:
		if (!dosplit)
		    vname = argv[k+2];
		else
		    for (j = 0; j < nmat; ++j)
		        setname(subname,j,argv[k+2]);
	    case 2:
		seedname = argv[k+1];
	    case 1:
		for (j = 0; j < nmat; ++j)
		    setname(genname,j,argv[k]);
	    case 0:
		break;
	    default:
		errexit(ERR_BADUSAGE,"zsp");
	}
    }
}




/* ------------------------------------------------------------------
   init() - Read input files and initialize
   ------------------------------------------------------------------ */

static char wm1[] =
   "zsp: Warning: Cannot calculate quotient for Permutations\n";
static char wm2[] =
   "zsp: Warning: Only one permutation read from %s\n";

void init(void)

{
    long fl2, noc, nor2, noc2;
    int j;
    FILE *f;

    /* Open the <Seed> file
       -------------------- */
    if (MSG2) printf("Open seed file (%s)\n",seedname);
    if ((seedfile = zreadhdr(seedname,&fl,&ssdim,&noc)) == NULL)
	errexit(-1,seedname);
    if (MSG1) printf("Found %ld seed vectors\n",ssdim);
    if (fl < 2) errexit(ERR_NOTMATRIX,seedname);
    nor = noc;


    /* Read the generators
       ------------------- */
    for (j = 0; j < nmat; ++j)	/* Read the generators */
    {
	if ((f = zreadhdr(genname[j],&fl2,&nor2,&noc2)) == NULL)
	    errexit(-1,genname[j]);
	if (fl2 == -1)	/* Permutations */
	{
	    if (gentype == MATRIX)
		errexit(ERR_NOTPERM,genname[j]);
	    if (gentype == DONTKNOW && dosplit) fprintf(stderr,wm1);
	    gentype = PERMUTATION;
	    if (noc2 > 1)
		fprintf(stderr,wm2,genname[j]);
	    if (nor2 != noc)
	    {
		char buf[200];
		sprintf(buf,"%s and %s",genname[j],seedname);
		errexit(ERR_INCOMPAT,buf);
	    }
	    gen[j] = (PTR) malloc(sizeof(long)*noc);
	    if (gen[j] == NULL) errexit(ERR_NOMEM,"zsp");
	    if (zreadlong(f,(long *)(gen[j]),noc) != noc)
		errexit(-1,genname[j]);
	}
	else
	{
	    if (gentype == PERMUTATION)
		errexit(ERR_NOTMATRIX,genname[j]);
	    gentype = MATRIX;
	    if (fl2 != fl || noc2 != noc)
	    {
		char buf[200];
		sprintf(buf,"%s and %s",genname[j],seedname);
		errexit(ERR_INCOMPAT,buf);
	    }
	    if (noc2 != nor2) errexit(ERR_NOTSQUARE,genname[j]);
	    zsetfield(fl);
	    zsetlen(nor);
	    gen[j] = zalloc(nor);
	    zreadvec(f,gen[j],nor);
	}
	fclose(f);
    }

    /* Set default value of maxdim
       --------------------------- */
    if (maxdim == -1) maxdim = nor - 1; /* We test for dim > maxdim! */
    if (smode == READALL && maxdim < ssdim) err('S');

    /* Allocate memory
       --------------- */
    zsetfield(fl);
    zsetlen(nor);
    m3 = zalloc(maxdim + 2);
    tmp1 = zalloc((long)1);
    tmp2 = zalloc((long)1);
    piv = (long *)malloc((size_t)(nor+1)*sizeof(long));
    ins = (char *)malloc((size_t)(nor+1));
    if (piv == NULL || ins == NULL)
	errexit(ERR_NOMEM,"zsp");
}




/* ------------------------------------------------------------------
   do_split() - Calculate the action of g1 and g2 on subspace and
	quotient.
   ------------------------------------------------------------------ */

void do_split(void)

{
    int g;
    long j, qdim;
    PTR xa, quotop;
    FILE *ofile;
    FEL dummy;


    /* Calculate the action on the subspace
       ------------------------------------ */
    MESSAGE(1,("Splitting..."));
    for (g = 0; g < nmat; ++g)
    {	
	if ((ofile = zwritehdr(subname[g],fl,dim,dim)) == NULL)
	    errexit(-1,subname[g]);
        MESSAGE(2,("%s...",subname[g]));
	xa = m3;
	for (j = 1; j <= dim; ++j)
	{
	    zsetlen(dim);
	    zmulrow(tmp2,F_ZERO);

	    /* Map the j-th basis vector
               ------------------------- */
	    zsetlen(nor);
	    switch (gentype)
	    {
	    	case MATRIX:
		    zmaprow(xa,gen[g],nor,tmp1);
		    break;
	    	case PERMUTATION:
		    zpermrow(xa,(long *)gen[g],tmp1);
		    break;
		case DONTKNOW:
		    err('b');
	    }
	    zadvance(&xa,(long)1);

	    /* Clean it with the basis and store the
               coefficients in tmp2.
               ------------------------- */
	    zcleanrow2(tmp1,m3,dim,piv,tmp2);

	    /* Check if the row is clean
	       ------------------------- */
	    if (zfindpiv(tmp1,&dummy) != 0)
	        fprintf(stderr,"WARNING: SUBSPACE IS NOT INVARIANT\n");

	    /* Write one row of output
	       ----------------------- */
	    zsetlen(dim);
		zwritevec(ofile,tmp2,(long)1);
	}
	fclose(ofile);
    }


    /* Calculate the action on the quotient
       ------------------------------------ */
    if (gentype == PERMUTATION) return;
    qdim = nor - dim;		/* Dimension of the quotient */
    zsetlen(qdim);
    quotop = zalloc(qdim);
    zsetlen(nor);
    zquotinit(m3,dim,piv);
    for (g = 0; g < nmat; ++g)
    {
        MESSAGE(2,("%s...",quotname[g]));
	if ((ofile = zwritehdr(quotname[g],fl,qdim,qdim)) == NULL)
	    errexit(-1,quotname[g]);
    	zsetlen(nor);
	zquotop(gen[g],quotop);
	zsetlen(qdim);
	zwritevec(ofile,quotop,qdim);
	fclose(ofile);
    }
    MESSAGE(1,("Done\n"));
    zsetlen(nor);
    free(quotop);
}


/* ------------------------------------------------------------------
   wrmod() - Write out the submodule
   ------------------------------------------------------------------ */

void wrmod(long cnt)

{
    char buf[200];
    FILE *ofile;

    if (smode == READALL)
	strcpy(buf,vname);
    else
    {
	sprintf(buf,"%s.%4.4ld",vname,cnt);
    }
    if ((ofile = zwritehdr(buf,fl,dim,nor)) == NULL)
	errexit(-1,buf);
    zwritevec(ofile,m3,dim);
    fclose(ofile);
}



/* ------------------------------------------------------------------
   spin_up() - Spin up the seed vector (s).

   returns 1 if a proper subspace was found, 0 otherwise
   ------------------------------------------------------------------ */

int spin_up(void)

{
    int g;
    long count, nmu, pv;
    PTR xput, xget;
    FEL mk;
    int message = 1;

    xput = m3;
    dim = 0;

    /* If smode == READALL, echelonize the seed space.
       ----------------------------------------------- */
    if (smode == READALL)
    {
	dim = zmkechelon(m3,ssdim,piv);
	if (dim <= 0) err('0');
	MESSAGE(0,("SEED HAS DIMENSION %ld\n",dim));
	--dim;
	xput = m3;
	zadvance(&xput,dim);
    }

    /* Spin up
       ------- */
    g = -1;
    xget = m3;
    count = nmu = 0;
    while (1)
    {
	pv = zfindpiv(xput,&mk);
	if (pv != 0)
	{
	    count = 0;
	    message = 1;
	    piv[(int)++dim] = pv;
	    zadvance(&xput,(long)1);
	    if (dim > maxdim)
		return 0; /* Spans whole Space */
	}
    	if (message)
	{
	    message = 0;
	    if (MSG3)
	    	printf("dim = %ld, stack = %ld\n",dim,dim-nmu);
	    else if (MSG2 && dim % 10 == 0)
	    	printf("dim = %ld, stack = %ld\n",dim,dim-nmu);
	    else if (MSG1 && dim % 100 == 0)
		printf("dim=%ld\n",dim);
	}

	/* Apply the next generator
	   ------------------------ */
	if (++g >= ngen)
	{
	    g = 0;
	    zadvance(&xget,(long)1);
	    ++count;
	    if (++nmu >= dim ||
		(maxtry > 0 && count >= maxtry))
	    {
		if (dim == 0) err('0');
		if (dim < nor)
		    return 1;	/* Split! */
		else
		    return 0;    /* Spans whole space */
	    }
	    message = 1;
	}
	switch (gentype)
	{
	    case MATRIX:
		zmaprow(xget,gen[g],nor,xput);
		break;
	    case PERMUTATION:
		zpermrow(xget,(long *)gen[g],xput);
		break;
	    case DONTKNOW:
		err('b');
		break;
	}

	/* Clean the image with the basis.
	   ------------------------------- */
	zcleanrow(xput,m3,dim,piv);
    }
}




/* ------------------------------------------------------------------
   try_split()
   ------------------------------------------------------------------ */

void try_split(void)

{
    int hassplit;

    hassplit = spin_up();
    if (opt_G)
	printf("MeatAxe.SplitDimensions := [%ld,%ld];\n",dim,nor-dim);

    if (!dosplit)
    {
	if (smode != READALL)
	    MESSAGE(0,("VECTOR %5ld: ",iseed));
	if (!hassplit)
	{
	    if (dim < nor)
		MESSAGE(0,("DIMENSION LIMIT EXCEEDED\n"));
	    else
		MESSAGE(0,("SPANS WHOLE SPACE\n"));
	}
	else
	    MESSAGE(0,("SUB=%6ld QUOT=%6ld\n",dim,nor-dim));
	if (output)
	    wrmod(iseed);
    }
    else
    {	if (hassplit)
	{
	    if (smode != READALL)
	        MESSAGE(0,("VECTOR %ld SPLITS: ",iseed));
	    MESSAGE(0,("SUBSPACE %ld  QUOTIENT %ld\n",dim,nor-dim));
	    if (output)
		do_split();
	    exit(EXIT_OK);
	}
    }
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char *argv[])

{
    PTR base;

    mtxinit();
    parseargs(argc,argv);
    init();

    switch (smode)
    {
	case GENERATE:
	    base = zalloc(ssdim);
	    zreadvec(seedfile,base,ssdim);
	    zpseed_init(ssdim,base);
	    if (seedvecno == 0)
	    {
		iseed = zpseed_make((long)1);
		while (iseed >= 0)
		{
		    zmoverow(m3,zpseed_vec);
		    try_split();
		    iseed = zpseed_next();
		}
	    }
	    else
	    {
		iseed = zpseed_make(seedvecno);
		if (iseed < 0) err('v');
		zmoverow(m3,zpseed_vec);
		try_split();
	    }
	    break;
	case READ:
	    if (seedvecno == 0)
	    {
		for (iseed = 1; iseed <= ssdim; ++iseed)
		{
		    zreadvec(seedfile,m3,(long)1);
		    try_split();
		}
	    }
	    else
	    {	if (seedvecno > ssdim) err('v');
		zseek(seedfile,seedvecno);
		zreadvec(seedfile,m3,(long)1);
		iseed = seedvecno;
		try_split();
	    }
	    break;

	case READALL:
	    zreadvec(seedfile,m3,ssdim);
	    try_split();
	    break;
    }

    /* The representation was not split. Print an
       appropriate message depending on the options
       used and exit with the right exit code.
       --------------------------------------------- */
    if (!dosplit)
	return EXIT_OK;
    if (seedvecno == 0 && smode != READALL)
	MESSAGE(0,("EVERY SEED VECTOR "));
    MESSAGE(0,("SPANS WHOLE SPACE\n"));
    return EXIT_NOSPLIT;
}

