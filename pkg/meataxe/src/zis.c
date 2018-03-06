/* ========================== C MeatAxe =============================
   zis.c - Make invariant subspace.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zis.c,v 1.2 1997/09/11 15:43:55 gap Exp $
 *
 * $Log: zis.c,v $
 * Revision 1.2  1997/09/11 15:43:55  gap
 * New version 2.2.3. AH
 *
 * Revision 2.3  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.2  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.4  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.3  1993/08/05  15:52:49  mringe
 * File header
 *
 * Revision 1.2  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.1  1992/05/26  17:39:14  mringe
 * Initial revision
 *
 */

#include "meataxe.h"


/* ------------------------------------------------------------------
   err() - Print error message and exit.
   ------------------------------------------------------------------ */

static void err(int c)

{	fprintf(stderr,"ZIS ERROR - ");
	switch (c)
	{	case 'c':
			fprintf(stderr,"MATRICES INCOMPATIBLE\n");
		case '0':
			fprintf(stderr,"CANNOT SPIN UP ZERO VECTOR\n");
		case 's':
			fprintf(stderr,"Usage: `%s'\n",USAGE);
			break;
        }
	exit(EXIT_ERR);
}


int main(argc, argv)
int argc;
char *argv[];


{	long fl, fl2;
	long nor, noc, nor2, noc2, *piv;
	PTR m1, m2, m3, m4, xp, x3, x8, xe, *sad;
	long dim, nmu, pv, j, x7;
	FEL mk, f1, f2;
	int gen;

	mtxinit();
	switch (argc)
	{	case 1:
			break;
		case 5:
			zname(1,0,argv[1]);
			zname(2,0,argv[2]);
			zname(3,0,argv[3]);
			zname(4,1,argv[4]);
			break;
		default:
			err('s');
	}
	zReadhdr(1,&fl,&nor,&noc);
	if (nor != noc) err('c');
	zsetfield(fl); zsetlen(noc);
	m1 = zalloc(nor);
	m2 = zalloc(nor);
	m3 = zalloc(nor);
	m4 = zalloc((long)1);
	piv = (long *) malloc((size_t)nor * sizeof(long));
	sad = (PTR *) malloc((size_t)nor * sizeof(PTR));
	zRead(1,m1,nor);
	ziclose(1);
	zReadhdr(2,&fl2,&nor2,&noc2);
	if (fl2 != fl || nor2 != nor || noc2 != noc)
		err('c');
	zRead(2,m2,nor);
	ziclose(2);
	zReadhdr(3,&fl2,&nor2,&noc2);
	if (noc2 != noc || fl2 != fl)
		err('c');
	zRead(3,m3,(long)1);
	ziclose(3);

	x3 = m3;
	dim = 0;
	nmu = 0;
	gen = 2;

	while (1)
	{

	pv = zfindpiv(x3,&mk);
	if (pv != 0)
	{	if (mk != F_ONE)
		{	f1 = zinv(mk);
			zmulrow(x3,f1);
		}
		++dim;
		piv[dim] = pv;
		sad[dim] = x3;
		zadvance(&x3,(long)1);
		if (dim >= noc)
		{	fprintf(stderr,"SPANS WHOLE SPACE\n");
			return (EXIT_OK);
		}
	}
	++gen;
	if (gen > 2)
	{	gen = 1;
		if (++nmu > dim)
		{	if (dim == 0) err('0');
			fprintf(stderr,
				"SUBSPACE HAS DIMENSION %ld\n",dim);
			zWritehdr(4,fl,dim,noc);
			zWrite(4,m3,dim);
			zoclose(4);
			return (EXIT_OK);
		}
	}
	xp = m1;
	if (gen == 2) xp = m2;
	xe = sad[nmu];

	zmaprow(xe,xp,noc,x3);

	for (j = 1; j <= dim; ++j)
	{	x7 = piv[(int)j];
		f1 = zextract(x3,x7);
		if (f1 != F_ZERO)
		{	f2 = zneg(f1);
			x8 = sad[j];
			if (f2 == F_ONE)
				zaddrow(x3,x8);
			else
			{	zmoverow(m4,x8);
				zmulrow(m4,f2);
				zaddrow(x3,m4);
			}
		}
	}

	}
}

