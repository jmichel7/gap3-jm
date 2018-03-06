/* ========================== C MeatAxe =============================
   znu.c - This program calculates the null space of a matrix.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: znu.c,v 1.2 1997/09/11 15:44:05 gap Exp $
 *
 * $Log: znu.c,v $
 * Revision 1.2  1997/09/11 15:44:05  gap
 * New version 2.2.3. AH
 *
 * Revision 2.4  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.3  1994/07/22  18:52:59  mringe
 * znullsp() liefert jetzt den fertigen Nullraum.
 *
 * Revision 2.2  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.8  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.7  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.6  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.5  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.4  1993/07/19  15:26:52  mringe
 * Optionen -G, -Q, -V
 *
 * Revision 1.3  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.2  1992/08/13  13:13:43  mringe
 * Output in semi-echelon form.
 *
 * Revision 1.1  1992/05/26  17:39:30  mringe
 * Initial revision
 *
 */



#include <stdlib.h>
#include "meataxe.h"





/* ------------------------------------------------------------------
   Variables
   ------------------------------------------------------------------ */

static char *helptext[] = {
"SYNTAX",
"    znu [-GQV] <Matrix> [<NullSpace>]",
"",
"OPTIONS",
"    -G     GAP output",
"    -Q     Quiet, no messages",
"    -V     Verbose, more messages",
"",
"FILES",
"    <Matrix>       i  Any matrix.",
"    <NullSpace>    o  The null space.",
NULL};

static proginfo_t pinfo =
   { "znu", "Matrix Null-Space", "$Revision: 1.2 $", helptext };


int opt_G = 0;
static char *matname;
static char *nspname = NULL;
static FILE *matfile, *nspfile;



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{   long nor, noc, dim;
    long *piv;
    long fl;
    int i;
    PTR m1, m2;

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
    switch (argc - opt_ind)
    {
        case 2:
	    nspname = argv[opt_ind+1];
        case 1:
	    matname = argv[opt_ind];
	    break;
	default:
	    errexit(ERR_NARGS,"znu");
    }

    /* Read input file, check sizes
       ---------------------------- */
    if ((matfile = zreadhdr(matname,&fl,&nor,&noc)) == NULL)
	errexit(-1,matname);
    if (fl < 0) errexit(ERR_NOTMATRIX,matname);
    piv = (long *) malloc((size_t)(sizeof(long)*(nor+1)));
    zsetfield(fl);
    zsetlen(nor);
    m2 = zalloc(nor);
    zsetlen(noc);
    m1 = zalloc(nor);
    zreadvec(matfile,m1,nor);
    fclose(matfile);

    /* Calculate null-space
       -------------------- */
    dim = znullsp(m1,nor,piv,m2);
    if (opt_G)
    {
	printf("MeatAxe.Nullity := %ld;\n",dim);
	fflush(stdout);
    }
    else
        if (MSG0) printf("NULLITY %ld\n",dim);

    /* Write out the null space
       ------------------------ */
    if (nspname != NULL)
    {
        if (MSG1) printf("Writing null-space to %s\n",nspname);
    	if ((nspfile = zwritehdr(nspname,fl,dim,nor)) == NULL)
	    errexit(-1,nspname);
    	zsetlen(nor);
	zwritevec(nspfile,m2,dim);
        fclose(nspfile);
    }
    return 0;
}



