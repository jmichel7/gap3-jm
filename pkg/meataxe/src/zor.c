/* ========================== C MeatAxe =============================
   zor.c - Order of a matrix or permutations.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zor.c,v 1.2 1997/09/11 15:44:08 gap Exp $
 *
 * $Log: zor.c,v $
 * Revision 1.2  1997/09/11 15:44:08  gap
 * New version 2.2.3. AH
 *
 * Revision 2.7  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.6  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.5  1994/03/13  13:27:01  mringe
 * Maschinenunabhaengiges Format fuer Permutationen.
 *
 * Revision 2.4  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
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
 * Revision 1.14  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.13  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.12  1993/08/26  17:31:35  mringe
 * Helptext.
 *
 * Revision 1.11  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.10  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.9  1993/07/19  16:32:16  mringe
 * help(), Optionen -G, -Q, -V
 *
 * Revision 1.8  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.7  1993/02/16  18:32:46  mringe
 * string.h und stdio.h werden jetzt in meataxe.h included.
 *
 * Revision 1.6  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.5  1993/01/06  20:59:57  mringe
 * getopt in() in zgetopt() umbenannt.
 *
 * Revision 1.4  1992/07/29  10:02:38  mringe
 * Neu: Option -q (Quick) und -m <Max>
 *
 * Revision 1.3  1992/07/23  06:38:03  mringe
 * Bug in matorder() beseitigt.
 *
 * Revision 1.2  1992/07/22  06:28:04  mringe
 * Betrachte zyklische Teilraeume.
 *
 * Revision 1.1  1992/05/26  17:39:36  mringe
 * Initial revision
 */


#define MAXORDER 100000		/* Maximal order */
#define MAXORDER_C 1000		/* Maximal order on cyclic subspaces */


#include <stdlib.h>
#include "meataxe.h"




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

char *iname;
FILE *ifile;
long fl, nor, noc;
long maxord = -1;	/* Limit set with -m option */
int opt_q = 0;		/* Quick mode */
int opt_G = 0;		/* GAP output */

static char *helptext[] = {
"SYNTAX",
"    zor [-m <MaxOrder>] [-GQVq] <File>",
"",
"OPTIONS",
"    -m   Set highest possible order.",
"    -q   Quick mode: Find a lower bound for the order.",
"    -G   GAP output.",
"    -Q   Quiet, no messages.",
"    -V   Verbose, more messages.",
"",
"FILES",
"    <File>   i   A matrix or permutation(s).",
NULL};
static proginfo_t pinfo =
   { "zor", "Order of Matrices Or Permutations",
     "$Revision: 1.2 $", helptext };



/* ------------------------------------------------------------------
   ordmat() - Order of a matrix
   ------------------------------------------------------------------ */

static void ordmat()

{   PTR m1, v;
    PTR base, bend;
    long dim;
    long *piv;
    char *ispiv;
    long ord, i;

    zsetfield(fl); zsetlen(noc);
    if (nor != noc) errexit(ERR_NOTSQUARE,iname);
    m1 = zalloc(nor);
    base = zalloc(nor+1);
    piv = (long *) malloc((size_t)(nor+1)*sizeof(long));
    ispiv = (char *) malloc((size_t)(nor+1));
    for (i = 1; i <= nor; ++i) ispiv[i] = 0;
    v = zalloc((long)1);
    zreadvec(ifile,m1,nor);
    ord = 1;
    bend = base;

    for (dim = 0; dim < nor; )
    {	PTR start = bend;
	long tord = 0;
	int closed = 0;

	/* Find the next seed vector
	   ------------------------- */
	for (i = 1; i <= nor && ispiv[i]; ++i);
	zmulrow(bend,F_ZERO);
	zinsert(bend,i,F_ONE);

	/* Calculate order on the cyclic subspace
	   -------------------------------------- */
	do
	{   PTR b;
	    FEL f;
	    long pv;

	    /* Save the vector and extend the basis,
	       if the vector is linearly independent.
	       -------------------------------------- */
	    zmoverow(v,bend);
	    if (!closed)
	    {	b = base;
	    	for (i = 1; i <= dim; ++i)
	    	{   f = zextract(bend,piv[i]);
		    if (f != F_ZERO)
		    {	zaddmulrow(bend,b,zneg(
			zdiv(f,zextract(b,piv[i]))));
		    }
		    zadvance(&b,(long)1);
		}
		pv = zfindpiv(bend,&f);
		if (pv != 0)
		{   piv[++dim] = pv;
		    ispiv[piv[dim]] = 1;
		    zadvance(&bend,(long)1);
		}
		else
		    closed = 1;
	    }

	    /* Apply the matrix.
	       ----------------- */
	    if (++tord > MAXORDER_C)
	    {  
		fprintf(stderr,"zor: Partial order is over %ld\n",
		    (long) MAXORDER_C);
		exit(1);
	    }
	    zmaprow(v,m1,nor,bend);
	}
	while (zcmprow(bend,start));

	/* Calculate l.c.m. of all tord's
	   ------------------------------ */
	for (i = ord; ord % tord != 0; ord += i);
	if (ord > MAXORDER)
	{
	    fprintf(stderr,"zor: Order is over %ld\n",(long)MAXORDER);
	    exit(1);
	}
	if (opt_q && dim > nor/10) break;
	if (maxord > 1)
	{   if (ord > maxord)
	    {
	    	fprintf(stderr,"zor: Order is over %ld\n",maxord);
	    	exit(1);
	    }
	    if (ord == maxord) break;
	}
    }

    if (opt_q && ord != maxord)
	MESSAGE(0,("ORDER IS A MULTIPLE OF %ld\n",(long int) ord));
    else
    {
	if (opt_G)
	    printf("MeatAxe.Order := %ld;\n",ord);
	else
	    printf("ORDER IS %ld\n",ord);
    }
}


/* ------------------------------------------------------------------
   ordperm() - Order of permutations
   ------------------------------------------------------------------ */

static void ordperm()

{
    long order, iper;
    perm_t *perm;

    if (fl != -1) errexit(ERR_NOTPERM,iname);	/* No monomials */
    if ((perm = permalloc(nor)) == NULL)
	errexit(-1,"zor: ordperm()");
    if (opt_G) printf("MeatAxe.Orders := [");
    for (iper = 1; iper <= noc; ++iper)
    {
	if (zreadlong(ifile,(long*)(perm->d),nor) != nor)
	    errexit(-1,iname);
	order = permorder(perm);
	if (order < 0)
	{
	    mtxerror("zor");
	    continue;
	}

	if (opt_G)
	{
	    if (iper > 1) printf(",");
	    printf("%ld",order);
	}
	else
	    printf("ELEMENT %ld HAS ORDER %ld\n",iper, order);
    }
    if (opt_G) printf("];\n");
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{
    int i;

    /* Parse command line
       ------------------ */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("Gqm:")) != OPT_END)
    {
	switch (i)
	{
	    case 'G': opt_G = 1; msg_level = -100; break;
	    case 'q': opt_q = 1; break;
	    case 'm':
		maxord = getint();
		if (maxord < 1 || maxord == GETINT_ERR)
		    errexit(ERR_OPTION,"-m");
		break;
	}
    }
    if (opt_ind != argc-1) 
	errexit(ERR_BADUSAGE,"zor");
    iname = argv[opt_ind];

    /* Open input file, call the appropriate function
       ---------------------------------------------- */
    if ((ifile = zreadhdr(iname,&fl,&nor,&noc)) == NULL)
	errexit(-1,iname);
    if (fl < 0)
	ordperm();
    else
	ordmat();

    fclose(ifile);
    return (EXIT_OK);
}


