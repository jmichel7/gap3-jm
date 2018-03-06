/* ========================== C MeatAxe =============================
   ztc.c - Trace of a matrix or permutation.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: ztc.c,v 1.2 1997/09/11 15:44:33 gap Exp $
 *
 * $Log: ztc.c,v $
 * Revision 1.2  1997/09/11 15:44:33  gap
 * New version 2.2.3. AH
 *
 * Revision 2.8  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.7  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.6  1994/03/13  13:27:01  mringe
 * Maschinenunabhaengiges Format fuer Permutationen.
 *
 * Revision 2.5  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
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
 * Revision 1.16  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.15  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.14  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.13  1993/07/23  13:46:27  mringe
 * OS-Symbole neu (SYS_xxx)
 *
 * Revision 1.12  1993/07/16  14:09:18  mringe
 * Option -G
 *
 * Revision 1.11  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.10  1993/07/09  14:41:34  mringe
 * help() verbessert, -Q und -V eingefuehrt
 *
 * Revision 1.9  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.8  1993/01/07  14:52:49  mringe
 * Semikolon bei GAP-Output
 *
 * Revision 1.7  1993/01/07  13:26:29  mringe
 * *** empty log message ***
 *
 * Revision 1.6  1993/01/06  21:12:20  mringe
 * GAP-Output f"ur Matrizen.
 *
 * Revision 1.5  1993/01/06  20:59:57  mringe
 * getopt in() in zgetopt() umbenannt.
 *
 * Revision 1.4  1993/01/06  20:39:17  mringe
 * help(), GAP output f"ur Permutationen.
 *
 */


#include <stdlib.h>
#include "meataxe.h"



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

long fl, nor, noc;
int opt_G = 0;		/* GAP output */
int msg = 0;
FILE *inp;
char *inpname;

static char *helptext[] = {
"SYNTAX",
"    ztc [-GQV] <File>",
"",
"OPTIONS",
"    -G   GAP output",
"    -V   Verbose output (repeat to increase number of messages)",
"    -Q   Quiet, no messages",
"",
"FILES",
"    <File>    i  A square matrix or permutation(s)",
"",
NULL};

static proginfo_t pinfo =
   { "ztc", "Trace", "$Revision: 1.2 $", helptext };




/* ------------------------------------------------------------------
   trmat() - Trace of a matrix
   ------------------------------------------------------------------ */

static void trmat()

{
    FEL tr;
    long i, max;
    PTR m1;

    zsetfield(fl); zsetlen(noc);
    m1 = zalloc((long)1);
    tr = F_ZERO;
    if ((max = nor) > noc) max = noc;
    for (i = 1; i <= max; ++i)
    {
	zreadvec(inp,m1,1);
	tr = zadd(tr,zextract(m1,i));
    }
    if (!opt_G)		/* Standard output */
	printf("TRACE IS %ld\n",zftoi(tr));
    else		/* GAP output */
        printf("MeatAxe.Trace := %s;\n",zftogap(tr));
}


/* ------------------------------------------------------------------
   trperm() - Trace of a permutation (= number of fixed points)
   ------------------------------------------------------------------ */

static void trperm()

{
    long *m1;
    long i,k,tr=0;

    m1 = (long *) malloc(sizeof(long)*nor);
    if (m1 == NULL) errexit(ERR_NOMEM,"ztc");
    if (opt_G)
	printf("MeatAxe.Trace := [");
    for (i = 1; i <= noc; ++i)
    {
	if (zreadlong(inp,m1,nor) != nor)
	    errexit(-1,inpname);
	tr = 0;
	--m1;
	for (k = 1; k <= nor; ++k)
	    if (m1[k] == k) ++tr;

	if (!opt_G)	/* Standard output */
	    printf("ELEMENT %3ld HAS TRACE %6ld\n",i,tr);
	else		/* GAP output */
	{
	    if (i > 1) printf(",");
	    if (i % 10 == 0) printf("\n   ");
	    printf("%ld",tr);
	}
    }
    if (opt_G)
	printf("];\n");
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
    while ((i = zgetopt("G")) != OPT_END)
    {
	switch (i)
	{
	    case 'G': opt_G = 1; msg_level = -100; break;
	}
    }
    if (opt_ind != argc-1) errexit(ERR_NARGS,"ztc");
    inpname = argv[opt_ind];
    if ((inp = zreadhdr(inpname,&fl,&nor,&noc)) == NULL)
	errexit(-1,inpname);
    if (msg >= 1)
    {
	printf("Input is ");
	if (fl == -1)
	    printf("%ld permutation(s) of degree %ld\n",noc,nor);
	else
	    printf("a %ldx%ld matrix over GF(%ld)\n",nor,noc,fl);
    }
    if (fl < -1) errexit(ERR_BADTYPE,argv[opt_ind]);
    if (fl < 0)
	trperm();
    else
	trmat();
    fclose(inp);
    return (EXIT_OK);
}




