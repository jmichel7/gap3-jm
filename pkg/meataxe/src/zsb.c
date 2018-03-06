/* ========================== C MeatAxe =============================
   zsb.c - Standard basis

   (C) Copyright 1993-95 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zsb.c,v 1.2 1997/09/11 15:44:23 gap Exp $
 *
 * $Log: zsb.c,v $
 * Revision 1.2  1997/09/11 15:44:23  gap
 * New version 2.2.3. AH
 *
 * Revision 2.7  1995/06/20  12:09:56  mringe
 * Output operations.
 *
 * Revision 2.6  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.5  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.4  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.3  1994/01/28  08:25:41  mringe
 * Optionen -G, -Q arbeiten jetzt.
 *
 * Revision 2.2  1993/12/08  11:33:02  mringe
 * Neue CPU time - Funktionen.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.11  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.10  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.9  1993/09/22  08:29:23  mringe
 * Optionen -GQVT. Options -g verhaelt sich jetzt wie bei den
 * anderen programmen.
 *
 * Revision 1.8  1993/08/10  14:51:42  mringe
 * Include string.h
 *
 * Revision 1.7  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.6  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.5  1993/05/13  09:31:55  mringe
 * Fehler bei der Basisnumerierung korrigiert.
 *
 * Revision 1.4  1993/05/13  09:06:49  mringe
 * --
 *
 * Revision 1.3  1993/05/13  09:05:48  mringe
 * Fehlermeldungen korrigiert.
 *
 * Revision 1.2  1993/05/13  08:40:29  mringe
 * Fast komplett neu geschrieben: N Erzeuger, help(),
 * verbesserte Fehlermeldungen, Z-Library...
 *
 * Revision 1.1  1992/05/26  18:14:06  mringe
 * Initial revision
 *
 */

#include <string.h>
#include <stdlib.h>

#include "meataxe.h"
#include "lattice.h"	/* for MAXGEN */



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

char *genname[MAXGEN];	/* File name of generators */
char *seedname;		/* Seed vector file name */
char *outname;		/* Output file name */
char *OpName = NULL;	/* Name of <op> file */
FILE *seedfile, *outfile;
long fl, nor;		/* Field and size of the generators */
long *piv;		/* Pivot table */
long *OpTable;		/* Operations, written to <Op> */
#define OPVEC(i) OpTable[2*(i)]
#define OPGEN(i) OpTable[2*(i)+1]
int ngen = 2;		/* Number of generators */
PTR gen[MAXGEN];	/* Generators */
PTR seed, A, B;		/* Workspace */
long seeddim;		/* Number of seed vectors in input file */
long nseed;		/* Number of seed vectors read */
long *cycldim = NULL;	/* Dimension of cyclic subspaces */
int opt_G = 0;		/* GAP output */

static char *helptext[] = {
"SYNTAX",
"    zsb [<Options>] <Gen1> <Gen2> <Seed> <Basis> [<Op>]",
"    zsb [<Options>] [-g <#Gen>] <Gen> <Seed> <Basis> [<Op>]",
"",
"OPTIONS",
"    -Q            Quiet, no messages.",
"    -V            Verbose, more messages.",
"    -G            Produce GAP output.",
"    -T <Seconds>  Set CPU time limit.",
"    -Q            Quiet, no messages.",
"    -g <#Gen>     Set number of generators.",
"",
"FILES",
"    <Gen1/2>   i  The generators (square matrices of equal size)",
"    <Seed>     i  Seed vectors (a matrix over the same field and",
"                  with the same row size as the generators)",
"    <Basis>    o  Standard basis",
"    <Op>       o  Description of <Basis> in terms of seed vectors",
"                  and generators (integer matrix)",
"",
"    If -g is used, Generators are read from <Gen>.1, <Gen>.2, ...",
NULL};

static proginfo_t pinfo =
   { "zsb", "Standard Basis", "$Revision: 1.2 $", helptext };



/* ------------------------------------------------------------------
   init()
   ------------------------------------------------------------------ */

static void init()

{
    int i;
    long fl2, nor2, noc2;
    FILE *f;

    mtxinit();

    /* Read generators
       --------------- */
    for (i = 0; i < ngen; ++i)
    {
	if ((f = zreadhdr(genname[i],&fl2,&nor2,&noc2)) == NULL)
	    errexit(-1,genname[i]);
	if (fl2 < 2) errexit(ERR_NOTMATRIX,genname[i]);
	if (nor2 != noc2) errexit(ERR_NOTSQUARE,genname[i]);
	if (i == 0)
	{
	    fl = fl2;
	    nor = nor2;
	}
	else
	    if (fl != fl2 || nor != nor2)
	    {
		char buf[200];
		sprintf(buf,"%s and %s",genname[0],genname[i]);
		errexit(ERR_INCOMPAT,buf);
	    }
	zsetfield(fl);
	zsetlen(nor);
	gen[i] = zalloc(nor);
	zreadvec(f,gen[i],nor);
	fclose(f);
    }

    /* Allocate workspace
       ------------------ */
    A = zalloc(nor+1);
    B = zalloc(nor+1);
    seed = zalloc((long)1);
    piv = NALLOC(long,nor+2);
    OpTable = NALLOC(long,2*(nor + 2));

    /* Open seed vector file
       --------------------- */
    if ((seedfile = zreadhdr(seedname,&fl2,&seeddim,&noc2)) == NULL)
	errexit(-1,seedname);
    cycldim = (long *) malloc((size_t)seeddim * sizeof(long));
    if (cycldim == NULL) errexit(ERR_NOMEM,"zsb");
    if (fl2 < 2) errexit(ERR_NOTMATRIX,seedname);
    if (noc2 != nor || fl2 != fl)
    {
	char buf[200];
	sprintf(buf,"%s and %s",genname[0],seedname);
	errexit(ERR_INCOMPAT,buf);
    }
    nseed = 0;
}


/* ------------------------------------------------------------------
   init_args()
   ------------------------------------------------------------------ */

static void init_args(int argc, char **argv)

{
    int i;
    int opt_g = 0;

    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("Gg:")) != OPT_END)
    {
	switch (i)
	{
	    case 'G': opt_G = 1; msg_level = -1000; break;
	    case 'g':
		ngen = getint();
		opt_g = 1;
		if (ngen < 2 || ngen > MAXGEN || *opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-g");
		break;
	}
    }

    if (opt_g)
    {
	char buf[200];
	switch (argc - opt_ind)
	{
	    case 4:
		OpName = argv[opt_ind + 3];
	    case 3:
    		seedname = argv[opt_ind + 1];
    		outname = argv[opt_ind + 2];
		break;
	    default:
	    	errexit(ERR_NARGS,"zsb");
	}
	for (i = 0; i < ngen; ++i)
	{
	    sprintf(buf,"%s.%d",argv[opt_ind],i+1);
	    genname[i] = malloc(strlen(buf)+1);
	    strcpy(genname[i],buf);
	}
    }
    else	/* Two generators */
    {
	switch (argc - opt_ind)
	{
	    case 5:
		OpName = argv[opt_ind+4];
	    case 4:
		genname[0] = argv[opt_ind+0]; 
		genname[1] = argv[opt_ind+1]; 
    		seedname = argv[opt_ind+2];
    		outname = argv[opt_ind+3];
		break;
	    default:
	    	errexit(ERR_NARGS,"zsb");
	}
    }
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char *argv[])

{
    long i, dim, src, olddim;
    PTR a, b, aget, bget;
    int g;
    FEL f;

    init_args(argc,argv);
    init();

    a = A; b = B;
    dim = 0;
    while (dim < nor && nseed < seeddim)
    {

	/* Read next seed vector. Check if it is
	   contained in the subspace found so far.
  	   --------------------------------------- */
	zreadvec(seedfile,a,(long)1);
	zmoverow(b,a);
	cycldim[nseed] = 0;
	++nseed;
	OPVEC(dim+1) = nseed;
	OPGEN(dim+1) = 0;	/* It's a seed vector */
	zcleanrow(a,A,dim,piv);
	i = zfindpiv(a,&f);
	if (i == 0) continue;	/* Discard this vector */

	/* Initialize everything for spin-up.
   	   ---------------------------------- */
	olddim = dim;
	aget = a;
	bget = b;
	src = dim;
	zadvance(&a,(long)1);
	zadvance(&b,(long)1);
	piv[++dim] = i;
	g = 0;

	/* Spin up using the standard basis procedure
   	   ------------------------------------------ */
	while (dim < nor && aget != a)
	{
	    zmaprow(aget,gen[g],nor,a);
	    zcleanrow(a,A,dim,piv);
	    i = zfindpiv(a,&f);
	    if (i != 0)		/* We have a new vector */
	    {
		piv[++dim] = i;
		zmaprow(bget,gen[g],nor,b);
		zadvance(&a,(long)1);
		zadvance(&b,(long)1);
		OPVEC(dim) = src;
		OPGEN(dim) = g + 1;
	    }
	    if (++g >= ngen)
	    {
		g = 0;
		zadvance(&aget,(long)1);
		zadvance(&bget,(long)1);
		++src;
	    }
	}

	/* Print a message for this seed vector
  	   ------------------------------------ */
	MESSAGE(0,("VECTOR %ld SPANS %ld\n",nseed,dim-olddim));
	cycldim[nseed-1] = dim-olddim;
    }

    /* Done! Check if we have found a basis and print a
       warning message if not.
       ------------------------------------------------- */
    fclose(seedfile);
    if (dim < nor && msg_level >= 0)
	fprintf(stderr,"ZSB WARNING - SPAN IS ONLY %ld OF SPACE %ld\n",
	    dim,nor);
    if (nseed < seeddim && msg_level >= 0)
	fprintf(stderr,
	    "ZSB WARNING - ONLY %ld OF %ld SEED VECTORS WERE USED\n",
	    nseed,seeddim);

    /* Write out the result
       -------------------- */
    if (opt_G)
    {
	int k;

	printf("MeatAxe.Dimension := %ld;\n",dim);
	printf("MeatAxe.CyclicDimensions := [");
	for (k = 0; k < nseed; ++k)
	{
	    if (k > 0)
	    {
		printf(",");
		if (k % 10 == 0) printf("\n    ");
	    }
	    printf("%ld",cycldim[k]);
	}
	printf("];\n");
    }

    /* Write Standard basis
       -------------------- */
    if ((outfile = zwritehdr(outname,fl,dim,nor)) == NULL)
	errexit(-1,outname);
    zwritevec(outfile,B,dim);
    fclose(outfile);

    /* Write <Op> file
       --------------- */
    if (OpName != NULL)
    {
    	if ((outfile = zwritehdr(OpName,-T_IMAT,dim,2)) == NULL)
	    errexit(-1,OpName);
	zwritelong(outfile,OpTable+2,2*dim);
    	fclose(outfile);
    }
    return (EXIT_OK);
}



