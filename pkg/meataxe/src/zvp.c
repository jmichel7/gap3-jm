/* ========================== C MeatAxe =============================
   zvp.c - Vector permute (make permutations from matrices)

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zvp.c,v 1.2 1997/09/11 15:44:40 gap Exp $
 *
 * $Log: zvp.c,v $
 * Revision 1.2  1997/09/11 15:44:40  gap
 * New version 2.2.3. AH
 *
 * Revision 2.10  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.9  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.8  1994/07/23  19:46:32  mringe
 * Benutze MESSAGE
 *
 * Revision 2.7  1994/07/07  06:57:14  mringe
 * Include files.
 *
 * Revision 2.6  1994/05/10  16:27:57  mringe
 * sbasis in seedbasis umbenannt.
 *
 * Revision 2.5  1994/03/26  06:35:08  mringe
 * zseek(): "Anderung f"ur 64-Bit-Architekturen.
 *
 * Revision 2.4  1994/03/12  13:34:43  mringe
 * Benutze zwritelong().
 * Teste seed == 0.
 *
 * Revision 2.3  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
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
 * Revision 1.16  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.15  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.14  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.13  1993/06/09  18:55:22  mringe
 * Behandle -l <Num> richtig als maximale Bahnl"ange
 *
 * Revision 1.12  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.11  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.10  1993/01/06  20:59:57  mringe
 * getopt in zgetopt() umbenannt.
 *
 * Revision 1.9  1992/12/08  08:36:35  mringe
 * Fehlermeldung "out of seed vectors" geaendert.
 *
 * Revision 1.8  1992/11/27  10:18:48  mringe
 * Kleinere Bugs behoben.
 *
 * Revision 1.7  1992/11/27  09:54:25  mringe
 * Option -g fertig
 *
 * Revision 1.6  1992/11/20  09:07:00  mringe
 * Option -g angefangen; noch nicht voll implementiert.
 *
 */


#include <stdlib.h>
#include <string.h>
#include "meataxe.h"



#define MAXVEC 100000	/* Default value for maxvec */
#define MAXMATRIX 20


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

long maxvec = MAXVEC;	/* Max. orbit size */
long tabsize;		/* Size of hash table */
long nvec;		/* Total number of vectors obtained so far */
long nfinished;		/* Number of finished vectors */
long *vecpos;		/* Position of vectors in table */
long *vecno;		/* Numbers of vectors in table */
char *isfree;
PTR vtable;		/* Holds the vectors */
long *perm[MAXMATRIX];	/* Permutations */

long fl, dim;		/* Field and dimension */
PTR gen[MAXMATRIX];	/* Generators (Matrices) */
long nseed;		/* Number of vectors in <Seed> */
PTR seedbasis;		/* Vectors read from <Seed> */
PTR tmp;

int nmat = 2;		/* Number of matrices */
int ngen = 2;		/* Number of generators */
long startvec = 1;	/* First seed vector */
long iseed = -1;	/* Current seed vector */
enum {GENERATE,READ} seedmode = GENERATE;
			/* Seed mode */
int proj = 0;		/* Operate on 1-spaces rather than vectors */
int vecout = 0;		/* Output vectors, too */
int noout = 0;		/* No output, print orbit sizes only */
int opt_g = 0;		/* -g was used */
long hasm, hasp, neh;	/* Parameters for the hash function */

char *matname[MAXMATRIX];	/* Input file names (with -g) */
char *permname[MAXMATRIX];	/* Output file names (with -g) */
char *seedname, *orbname;

static char *helptext[] = {
"SYNTAX",
"    zvp [<Options>] <Mat1> <Mat2> <Seed> <Perm1> <Perm2> [<Orbit>]",
"    zvp -g <#Mat> [<Options>] <Mat> <Seed> <Perm> [<Orbit>]",
"",
"OPTIONS",
"    -g <#Mat>   Set number of matrices",
"    -n          No Output",
"    -p          Permute 1-spaces rather than vectors",
"    -v          Write Vectors to <Orbit>",
"    -l <Limit>  Set maximal orbit size",
"    -sg[<N>]    Generate seed, start with vector <N> (Default)",
"    -sr[<N>]    Read seed vectors, beginning with vector <N>.",
"",
"FILES",
"    <Mat1,2>    i  Square matrices of equal size",
"    <Seed>      i  A Matrix, NOC(Seed)=NOC(Mat1)",
"    <Perm1,2>   o  Permutations",
"    <Orbit>     o  The orbit, a matrix",
NULL};

static proginfo_t pinfo =
   { "zvp", "Vector Permute", "$Revision: 1.2 $", helptext };




/* ------------------------------------------------------------------
   err()
   ------------------------------------------------------------------ */

static void err(ch)
int ch;

{
    fprintf(stderr,"ZVP ERROR - ");
    switch (ch)
    {  
	case 'c':
	    fprintf(stderr,"MATRICES INCOMPATIBLE\n");
	    break;
    	case 'o':
	    fprintf(stderr,"TABLE OVERFLOW\n");
	    break;
	case 'x':
	    fprintf(stderr,"NO ORBIT < %ld WAS FOUND\n",maxvec);
	    break;
    }
    exit(EXIT_ERR);
}


/* ------------------------------------------------------------------
   init() - Process arguments, allocate memory
   ------------------------------------------------------------------ */

static void setname(name,i,template)
char **name;
int i;
char *template;

{
    char *c;

    c = malloc(strlen(template) + 10);
    sprintf(c,"%s.%d",template,i+1);
    name[i] = c;
}


static void init(argc, argv)
int argc;
char *argv[];

{   int i;

    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("g:pvnl:s:")) != OPT_END)
    {
	switch (i)
	{
	    case 'n': noout = 1; break;
	    case 'v': vecout = 1; break;
	    case 'p': proj = 1; break;
	    case 'g':
		opt_g = 1;
		ngen = nmat = getint();
		if (nmat == GETINT_ERR || nmat < 1 || nmat > MAXMATRIX)
		    errexit(ERR_OPTION,"-g");
		if (*opt_text_ptr == '.')
		{
		    ngen = getint();
		    if (ngen == GETINT_ERR || ngen < 1 || ngen > nmat)
		        errexit(ERR_OPTION,"-g");
		}
		if (*opt_text_ptr != 0) 
		    errexit(ERR_OPTION,"-g");
		break;
	    case 'l':
		maxvec = getint();
		if (maxvec == GETINT_ERR || maxvec < 1)
		    errexit(ERR_OPTION,"-l");
		break;
	    case 's':
		switch (*opt_text_ptr++)
		{
		    case 'g': seedmode = GENERATE; break;
		    case 'r': seedmode = READ; break;
		    default: errexit(ERR_OPTION,"-s");
		}
		if (*opt_text_ptr != 0)
		    startvec = getint();
		if (startvec < 1 || startvec == GETINT_ERR)
		    errexit(ERR_OPTION,"-s");
		break;
	}
    }

    /* File names
       ---------- */
    if (!opt_g)
    {
	switch (argc - opt_ind)
    	{
	    case 6:
	    	orbname = argv[opt_ind+5];
	    case 5:
	    	permname[0] = argv[opt_ind+3];
	    	permname[1] = argv[opt_ind+4];
    	    case 3:
	    	matname[0] = argv[opt_ind+0];
	    	matname[1] = argv[opt_ind+1];
	    	seedname = argv[opt_ind+2];
	    	break;
	    default:
    		errexit(ERR_BADUSAGE,"zvp");
    	}
    }
    else /* -g was used */
    {
	switch (argc - opt_ind)
    	{
	    case 4:
	    	orbname = argv[opt_ind+3];		/* Orbit */
	    case 3:
		for (i = 0; i < nmat; ++i)
		{
		    setname(permname,i,argv[opt_ind+2]);
		    setname(matname,i,argv[opt_ind]);
		}
	    	seedname = argv[opt_ind+1];
	    	break;
	    default:
    		errexit(ERR_BADUSAGE,"zvp");
    	}
    }


    /* Allocate tables
       --------------- */
    tabsize = maxvec + maxvec / 10 + 1;
    vecpos = (long *) malloc((size_t)((tabsize+1)*sizeof(long)));
    vecno = (long *) malloc((size_t)((tabsize+1)*sizeof(long)));
    isfree = (char *) malloc((size_t)((tabsize+1)*sizeof(char)));
    for (i = 0; i < nmat; ++i)
    {
	perm[i] = (long *) malloc((size_t)((maxvec+1)*sizeof(long)));
	if (perm[i] == NULL) errexit(ERR_NOMEM,"zvp");
    }
}



/* ------------------------------------------------------------------
   readfiles() - Read input files, allocate more memory.
   ------------------------------------------------------------------ */

static void readfiles()

{   long fl2, nor2, noc2;
    int i;
    FILE *f;

    /* Read the matrices
       ----------------- */
    for (i = 0; i < nmat; ++i)
    {
    	if ((f = zreadhdr(matname[i],&fl2,&nor2,&noc2)) == NULL)
	    errexit(-1,matname[i]);
    	if (fl2 < 2) errexit(ERR_NOTMATRIX,matname[i]);
    	if (nor2 != noc2) errexit(ERR_NOTSQUARE,matname[i]);
	if (i == 0)
	{
	    dim = nor2;
	    fl = fl2;
    	    zsetfield(fl);
    	    zsetlen(dim);
	}
	else
	{
	    if (dim != nor2 || fl != fl2) err('c');
	}

	/* Allocate memory and read the matrix
	   ----------------------------------- */
	gen[i] = zalloc(dim);
    	zreadvec(f,gen[i],dim);
    	fclose(f);
    }

    /* Allocate some more memory
       -------------------------- */
    tmp = zalloc((long)1);
    vtable = zalloc(tabsize+1);

    /* Read the seed vectors
       --------------------- */
    if ((f = zreadhdr(seedname,&fl2,&nseed,&noc2)) == NULL)
	errexit(-1,seedname);
    if (fl2 != fl || noc2 != dim) err('c');
    seedbasis = zalloc(nseed);
    zreadvec(f,seedbasis,nseed);
    fclose(f);

    if (seedmode == GENERATE)
        zpseed_init(nseed,seedbasis);
}


/* ------------------------------------------------------------------
   inithash() - Initialize the hash parameters
   ------------------------------------------------------------------ */

static void inithash()

{	long tod, i;
	FEL f1, f2;

	neh = (dim > (long) 25) ? (long) 25 : dim;
	tod = 100;
	hasp = tabsize;
	if (hasp < 110) tod = 7;
	if (hasp > 11)
	{	i = 1;
		while (++i <= tod)
		{	if (hasp % i == 0)
			{	--hasp;
				i = 1;
			}
		}
	}
	if (hasp < 100)
		hasm = 3;
	else
	{	f1 = F_ONE;
		f2 = F_ZERO;
		for (i = 1; i <= 89; ++i)
			f2 = zadd(f1,f2);
		if (f2 == F_ZERO)
			hasm = 83;
		else
			hasm = 89;
	}
}


/* ------------------------------------------------------------------
   normalize() - Normalize a vector
   ------------------------------------------------------------------ */

static void normalize(row)
PTR row;

{	FEL f;

	zfindpiv(row,&f);
	zmulrow(row,zinv(f));
}


/* ------------------------------------------------------------------
   hash() - The hash function
   ------------------------------------------------------------------ */

static long hash(row)
PTR row;

{	register long pos = 0, i;

	for (i = 1; i <= neh; ++i)
		pos = (pos * hasm + zftoi(zextract(row,i))) % hasp;
	return pos;
}



/* ------------------------------------------------------------------
   mkorbit() - Make the orbit
   ------------------------------------------------------------------ */

static void mkorbit()

{
    long pos, im, pos1;
    PTR x;
    int igen;	/* Apply which generator */

    igen = 0;
    while (nfinished < nvec && nvec <= maxvec)
    {
	x = vtable;
	zadvance(&x,vecpos[nfinished]);
	zmaprow(x,gen[igen],dim,tmp);
	if (proj) normalize(tmp);
	pos1 = pos = hash(tmp);
	x = vtable;
	zadvance(&x,pos);
	while (!isfree[pos] && zcmprow(tmp,x))
	{
	    if (++pos == tabsize)
	    {	pos = 0;
			x = vtable;
	    }
	    else
	    	zadvance(&x,(long)1);
	    /* The following should never happen as the
	       hash table is always larger than maxvec */
	    if (pos == pos1)	/* Table full */
		err('o');
	}

	if (isfree[pos])	/* New vector */
	{
	    zmoverow(x,tmp);
	    im = nvec;
	    isfree[pos] = 0;
	    vecpos[nvec] = pos;
	    vecno[pos] = nvec;
	    nvec++;
	}
	else
	    im = vecno[pos];
	
	perm[igen][nfinished] = im;

	/* Next generator
	   -------------- */
	if (++igen == ngen)
	{
	    igen = 0;
	    ++nfinished;
	}
    }
}


/* ------------------------------------------------------------------
   makeseed() - Make the next seed vector
   ------------------------------------------------------------------ */

PTR makeseed()

{
    PTR x;

    if (seedmode == GENERATE)
    {
	if (iseed == -1)
	    iseed = zpseed_make(startvec);
	else
	    iseed = zpseed_next();
	if (iseed == -1) err('x');
	return zpseed_vec;
    }
    else
    {
	if (iseed == -1)
	    iseed = startvec;
	else
	    ++iseed;
	if (iseed > nseed) err('x');
	x = seedbasis;
	zadvance(&x,iseed-1);
	return x;
    }
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{   
    int i;
    long pos;
    PTR x, seed;
    FILE *f;
    FEL a;

    init(argc, argv);
    readfiles();
    inithash();

    /* Main loop
       --------- */
    while (1)
    {
	    /* Read (or generate) next seed vector
	       ----------------------------------- */
	    do
	        seed = makeseed();
	    while (zfindpiv(seed,&a) == 0);
	    if (proj) normalize(seed);

	    /* Inititalize tables
	       ------------------ */
	    for (i = 0; i < (int) tabsize; ++i)
    		isfree[i] = 1;
	    nvec = 1;
	    nfinished = 0;
	    pos = hash(seed);
	    x = vtable;
	    zadvance(&x,pos);
	    zmoverow(x,seed);
	    isfree[pos] = 0;
	    vecpos[0] = pos;
	    vecno[pos] = 0;

	    /* Make the orbit and check for success
	       ------------------------------------ */
	    mkorbit();
	    if (nfinished < nvec)
	    {	MESSAGE(0,("ORBIT %ld HAS MORE THAN %ld VECTORS\n",
			iseed,maxvec));
		continue;
	    }
	    else
	    {	MESSAGE(0,("ORBIT %ld HAS %ld VECTORS\n",
			iseed,nvec));
		break;
	    }
	}


	/* Write vectors
	   ------------- */
	if (nvec > maxvec) return EXIT_ERR;
	if (noout) return (EXIT_OK);
	if (vecout)
	{	if ((f = zwritehdr(orbname,fl,nvec,dim)) == NULL)
		    errexit(-1,orbname);
		for (i = 0; i < (int) nvec; ++i)
		{	x = vtable;
			zadvance(&x,vecpos[i]);
			zwritevec(f,x,(long)1);
		}
		fclose(f);
	}

    /* Write permutations
       ------------------ */
    for (i = 0; i < nmat; ++i)
    {   long k, *pk;
	pk = perm[i];
	for (k = nvec; k != 0; --k) ++*(pk++);
	if ((f = zwritehdr(permname[i],(long)-1,nvec,(long)1)) == NULL)
	    errexit(-1,permname[i]);;
	zwritelong(f,perm[i],nvec);
	fclose(f);
    }
    return EXIT_OK;
}



