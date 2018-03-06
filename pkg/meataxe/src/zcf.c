/* ========================== C MeatAxe =============================
   Change field.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zcf.c,v 1.2 1997/09/11 15:43:38 gap Exp $
 *
 * $Log: zcf.c,v $
 * Revision 1.2  1997/09/11 15:43:38  gap
 * New version 2.2.3. AH
 *
 * Revision 2.6  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.5  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.4  1994/03/13  13:27:01  mringe
 * Maschinenunabhaengiges Format fuer Permutationen.
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
 * Revision 1.11  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.10  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.9  1993/08/31  08:16:41  mringe
 * Optionen -Q, -V, neues help().
 *
 * Revision 1.8  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.7  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.6  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.5  1992/07/22  14:06:38  mringe
 * Restrict to subfield.
 *
 * Revision 1.4  1992/07/15  16:33:36  mringe
 * Convert permutation to matrix
 *
 * Revision 1.3  1992/07/15  09:25:55  mringe
 * Some minor changes.
 *
 * Revision 1.2  1992/07/03  09:18:06  mringe
 * Prototypes korrigiert.
 *
 * Revision 1.1  1992/06/30  12:49:24  mringe
 * Initial revision
 *
 */


#include <stdlib.h>
#include "meataxe.h"





/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static char *iname, *oname;
FILE *inp, *out;
long fl;	/* Field parameter of input file */
long fl2;	/* Field to convert to */
long nor, noc;	/* Parameters of input file */

static char *helptext[] = {
"SYNTAX",
"    zcf [-QV] <Field> <Input> <Output>",
"",
"OPTIONS",
"    -Q     Quiet, no messages",
"    -V     Verbose, more messages",
"",
"ARGUMENTS",
"    <Field>    The field order of the output matrix",
"    <Input>    The input file, a Permutation or matrix",
"    <Output>   The output matrix",
NULL};


static proginfo_t pinfo =
   { "zcf", "Change Format", "$Revision: 1.2 $", helptext };




/* ------------------------------------------------------------------
   err() - Print error message and exit
   ------------------------------------------------------------------ */

static void err(c)
int c;

{   fprintf(stderr,"ZCF ERROR - ");
    switch (c)
    {
	case 'm':
	    errexit(ERR_NOMEM,"zcf");
	case '=':
	    fprintf(stderr,"MATRIX IS ALREADY OVER GF(%ld)\n",fl);
	    break;
	case 'e':
	    fprintf(stderr,"CANNOT CHANGE FROM GF(%ld) TO GF(%ld)\n",
		fl,fl2);
	    break;
	case 'f':
	    fprintf(stderr,"Invalid field order\n");
	    break;
    }
    exit(EXIT_ERR);
}





/* ------------------------------------------------------------------
   checkfl() - Check if GF(fl) < GF(fl2)
   ------------------------------------------------------------------ */

static int checkfl()

{	long f;
	
	if (fl2 < 2) err('f');
	if (fl == -1) return 3;
	if (fl == fl2) err('=');
	else if (fl < fl2)
	{	for (f = fl2; f % fl == 0; f /= fl);
		if (f == 1) return 1;
		else err('e');
	}
	else if (fl > fl2)
	{	for (f = fl; f % fl2 == 0; f /= fl2);
		if (f == 1) return 2;
		else err('e');
	}
	return 0;
}


/* ------------------------------------------------------------------
   permmat() - Convert permutation -> matrix
   ------------------------------------------------------------------ */

void permmat()

{	PTR m2;
	long *p, i;

    /* Read the permutation
       -------------------- */
    p = (long *) malloc(sizeof(long) * (size_t) nor);
    if (p == NULL) errexit(ERR_NOMEM,"zcf");
    if (zreadlong(inp,p,nor) != nor)
	errexit(-1,iname);
	fclose(inp);
	--p;

	/* Allocate one row
	   ---------------- */
	zsetfield(fl2);
	zsetlen(nor);
	m2 = zalloc((long)1);

	/* Create the output matrix
	   ------------------------ */
	if ((out = zwritehdr(oname,fl2,nor,nor)) == NULL)
	    errexit(-1,oname);
	for (i = 1; i <= nor; ++i)
	{	zmulrow(m2,F_ZERO);
		zinsert(m2,p[i],F_ONE);
		zwritevec(out,m2,1);
	}
	fclose(out);
	MESSAGE(0,("Converted to GF(%ld)\n",fl2));
}




/* ------------------------------------------------------------------
   restrictfield() - Restrict to subfield
   ------------------------------------------------------------------ */

void restrictfield()

{	PTR m1, x;
	long j, k;
	FEL *fp, *tmp;

	checkfl();
	zsetfield(fl);
	zsetlen(noc);
	m1 = zalloc(nor);
	tmp = (FEL *) malloc((size_t)(sizeof(FEL)*nor*noc));
	if (tmp == NULL) err('m');
	zreadvec(inp,m1,(size_t) nor);
	fclose(inp);

	fp = tmp;
	x = m1;
	for (j = 1; j <= nor; ++j)
	{	for (k = 1; k <= nor; ++k)
		{	*fp++ = zrestrict(zextract(x,k),fl2);
		}
		zadvance(&x,(long)1);
	}

	free(m1);

    zsetfield(fl2);
    zsetlen(noc);
    if ((out = zwritehdr(oname,fl2,nor,noc)) == NULL)
	errexit(-1,oname);
    m1 = zalloc((long)1);

    fp = tmp;
    for (j = 1; j <= nor; ++j)
    {
	zmulrow(m1,F_ZERO);
	for (k = 1; k <= nor; ++k)
	    zinsert(m1,k,*fp++);
	zwritevec(out,m1,(long)1);
    }
    fclose(out);
    MESSAGE(0,("Restricted to GF(%ld)\n",fl2));
}



/* ------------------------------------------------------------------
   matmat() - Convert matrix -> matrix
   ------------------------------------------------------------------ */

void matmat()

{
    PTR m1, x;
    long j, k;
    FEL *fp, *tmp;

    checkfl();
    zsetfield(fl); zsetlen(noc);
    m1 = zalloc(nor);
    tmp = (FEL *) malloc((size_t)(sizeof(FEL)*nor*noc));
    if (tmp == NULL) err('m');
    zreadvec(inp,m1,nor);
    fclose(inp);

    fp = tmp;
    x = m1;
    for (j = 1; j <= nor; ++j)
    {	for (k = 1; k <= nor; ++k)
		{	*fp++ = zextract(x,k);
		}
		zadvance(&x,(long)1);
    }

    free(m1);
    zsetfield(fl2);
    zsetlen(noc);
    if ((out = zwritehdr(oname,fl2,nor,noc)) == NULL)
	errexit(-1,oname);
    m1 = zalloc((long)1);

    fp = tmp;
    for (j = 1; j <= nor; ++j)
    {	zmulrow(m1,F_ZERO);
	for (k = 1; k <= nor; ++k)
		zinsert(m1,k,zembed(*fp++,fl));
	zwritevec(out,m1,1);
    }
    fclose(out);
    MESSAGE(0,("Embedded into GF(%ld)\n",fl2));
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
    while ((i = zgetopt("")) != OPT_END)
    {
    }
    if (opt_ind != argc - 3) errexit(ERR_NARGS,"zcf");
    iname = argv[opt_ind+1];
    oname = argv[opt_ind+2];
    fl2 = atol(argv[opt_ind+0]);
    if ((inp = zreadhdr(iname,&fl,&nor,&noc)) == NULL)
	errexit(-1,iname);

    /* Convert
       ------- */
    switch (checkfl())
    {	case 1:	/* Embed into larger field */
		matmat();
		break;
	case 2:
		restrictfield();
		break;
	case 3:
		permmat();
		break;
	default:
		err('e');
    }
    return EXIT_OK;
}

