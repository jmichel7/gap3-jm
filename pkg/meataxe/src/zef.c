/* ========================== C MeatAxe =============================
   zef.c - Reduce to (normalized) echelon form

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zef.c,v 1.2 1997/09/11 15:43:46 gap Exp $
 *
 * $Log: zef.c,v $
 * Revision 1.2  1997/09/11 15:43:46  gap
 * New version 2.2.3. AH
 *
 * Revision 2.3  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
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
 * Revision 1.12  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.11  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.10  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.9  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.8  1993/07/28  13:34:49  mringe
 * *** empty log message ***
 *
 * Revision 1.7  1993/07/20  08:44:45  mringe
 * Semikolon im GAP-Output.
 *
 * Revision 1.6  1993/07/19  14:09:50  mringe
 * Option -G
 *
 * Revision 1.5  1993/07/17  08:38:14  mringe
 * msg_t durch message_t ersetzt
 *
 * Revision 1.4  1993/07/16  14:05:19  mringe
 * Helptext, message library.
 *
 * Revision 1.3  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.2  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.1  1992/05/26  18:27:41  mringe
 * Initial revision
 *
 */




#include <stdlib.h>
#include "meataxe.h"




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static long *piv;
static int opt_G = 0;	/* GAP output */
static char *iname, *oname;
FILE *ifile, *ofile;

static char *helptext[] = {
"SYNTAX",
"    zef [-GQV] <Inp> <Out>",
"",
"OPTIONS",
"    -G   GAP Output format",
"    -Q   Quiet, no messages",
"    -V   Verbose, more messages",
"",
"FILES",
"    <Inp>    i  A matrix",
"    <Out>    o  Matrix reduced to echelon form",
NULL};

static proginfo_t pinfo =
   { "zef", "Echelon form", "$Revision: 1.2 $", helptext };



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{
    int i;
    long k, rank, fl, nor, noc;
    PTR m1, x;


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

    if (opt_ind != argc-2)
	errexit(ERR_NARGS,"zef");
    iname = argv[opt_ind];
    oname = argv[opt_ind+1];
	
    if ((ifile = zreadhdr(iname,&fl,&nor,&noc)) == NULL)
	errexit(-1,iname);
    if (fl < 1) errexit(ERR_NOTMATRIX,iname);
    if (MSG1) printf("%s: %ldx%ld-Matrix over GF(%ld)\n",
	iname,nor,noc,fl);
    piv = (long *) malloc((size_t)(sizeof(long) * (nor+1)));
    if (piv == NULL) errexit(ERR_NOMEM,"zef");
    zsetfield(fl); zsetlen(noc);
    m1 = zalloc(nor);
    zreadvec(ifile,m1,nor);
    fclose(ifile);

    /* Reduce to echelon form and normalize
       ------------------------------------ */
    rank = zmkechelon(m1,nor,piv);
    for (x=m1,k=1; k <= rank; ++k, zadvance(&x,(long)1))
	zmulrow(x,zinv(zextract(x,piv[k])));

    /* Write result
       ------------ */
    if ((ofile = zwritehdr(oname,fl,rank,noc)) == NULL)
	errexit(-1,oname);
    zwritevec(ofile,m1,rank);
    fclose(ofile);
    if (opt_G)
	printf("MeatAxe.Rank := %ld;\n",rank);
    else
        MESSAGE(0,("RANK %ld\n",rank));
    return (EXIT_OK);
}

