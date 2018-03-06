/* ========================== C MeatAxe =============================
   zsy.c - Symmetrized tensor product.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zsy.c,v 1.2 1997/09/11 15:44:31 gap Exp $
 *
 * $Log: zsy.c,v $
 * Revision 1.2  1997/09/11 15:44:31  gap
 * New version 2.2.3. AH
 *
 * Revision 2.7  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.6  1994/08/01  13:31:57  mringe
 * Pruefe, ob die Matrix hinreichend gross ist.
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
 * Revision 1.13  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.12  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.11  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.10  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.9  1993/07/23  13:46:27  mringe
 * OS-Symbole neu (SYS_xxx)
 *
 * Revision 1.8  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.7  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.6  1992/07/07  14:44:06  mringe
 * Bug in ze4() beseitigt.
 *
 * Revision 1.5  1992/07/06  17:40:25  mringe
 * E3 f"ur Permutationen.
 *
 * Revision 1.4  1992/07/06  17:08:15  mringe
 * S2 f"ur Permutationen.
 *
 * Revision 1.3  1992/07/06  09:16:09  mringe
 * E2 f"ur Permutationen.
 *
 * Revision 1.2  1992/06/01  08:34:25  mringe
 * CL warnings entfernt.
 *
 * Revision 1.1  1992/05/26  18:14:23  mringe
 * Initial revision
 *
 */

#include <stdlib.h>
#include "meataxe.h"




/* ------------------------------------------------------------------
   Global Data
   ------------------------------------------------------------------ */

static int opt_G = 0;		/* GAP output */
static char *iname, *oname;	/* File names */
static FILE *ofile;
enum {M_E2,M_E3,M_E4,M_S2,M_M3} mode;
static long fl;		/* Field */
static long nor, noc;	/* Input sizes */
static long nor2, noc2;	/* Output sizes */
static PTR m1;		/* Input */
static PTR m2;		/* One row of the output matrix */
static PTR *row;	/* Pointers to the rows of the input matrix */


static char *helptext[] = {
"SYNTAX",
"    zsy [-GQV] <Mode> <Inp> <Out>",
"",
"    Valid modes are",
"        e2       Antisymmetric (exterior) square",
"        e3       Antisymmetric (exterior) cube",
"        e4       Antisymmetric fourth power",
"        s2       Symmetric square",
"        m3       Tensor cube, mixed symmetry",
"",
"OPTIONS",
"    -G     GAP output",
"    -Q     Quiet, no messages",
"    -V     Verbose, more messages",
"",
"FILES",
"    <Inp>   i  The input matrix (must be square)",
"    <Out>   o  The output matrix",
"",
NULL};

static proginfo_t pinfo =
   { "zsy", "Symmetrized Tensor Product",
     "$Revision: 1.2 $", helptext };



/* ------------------------------------------------------------------
   parseargs() - Process command line arguments
   ------------------------------------------------------------------ */

static void parseargs(int argc, char **argv)

{
    int i;
    char *c;

    /* Parse command line
       ------------------ */
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("G")) != OPT_END)
    {
	switch (i)
	{
	    case 'G': opt_G = 1; msg_level = -100; break;
	}
    }
    if (opt_ind != argc - 3) errexit(ERR_NARGS,pinfo.name);
    iname = argv[opt_ind+1];
    oname = argv[opt_ind+2];

    /* Process the <Mode> argument
       --------------------------- */
    c = argv[opt_ind];
    if (!strcmp(c,"e2")) mode = M_E2;
    else if (!strcmp(c,"e3")) mode = M_E3;
    else if (!strcmp(c,"e4")) mode = M_E4;
    else if (!strcmp(c,"s2")) mode = M_S2;
    else if (!strcmp(c,"m3")) mode = M_M3;
    else errexit(ERR_BADUSAGE,pinfo.name);
}



/* ------------------------------------------------------------------
   init() - Initialize everything, read input files, allocate memory
   ------------------------------------------------------------------ */

static char wm1[] =
    "zsy: Warning: %s contains more than one permutation";

static void init()

{
    long i;
    FILE *f;

    /* Read input
       ---------- */
    if ((f = zreadhdr(iname,&fl,&nor,&noc)) == NULL)
	errexit(-1,iname);
    if (mode != M_E2 && mode != M_S2 && mode != M_E3 && fl < 2)
	errexit(ERR_NOTMATRIX,iname);
    if (mode == M_M3 && nor != noc)
	errexit(ERR_NOTSQUARE,iname);

    if (fl >= 2)	/* Matrix */
    {
	zsetfield(fl); zsetlen(noc);
	m1 = zalloc(nor);
	zreadvec(f,m1,nor);

	/* Set up pointers to the rows of the input matrix
   	   ----------------------------------------------- */
	row = NALLOC(PTR,nor+1);
	if (row == NULL)
	    errexit(ERR_NOMEM,pinfo.name);
	for (i = 1; i <= nor; ++i)
	{	row[i] = m1;
		zadvance(&m1,(long)1);
	}
    }
    else		/* Permutation */
    {	
	if (noc != 1) fprintf(stderr,wm1,iname);
	noc = nor;
	m1 = (PTR) Smalloc(sizeof(long)*(size_t)nor);
	if (m1 == NULL) errexit(ERR_NOMEM,pinfo.name);
	if (zreadlong(f,(long*) m1,nor) != nor)
	    errexit(ERR_FILEREAD,iname);
    }
    fclose(f);

    /* Calculate nor2 and noc2
       ----------------------- */
    switch (mode)
    {	case M_S2:
	    nor2 = (nor * (nor+1)) / 2;
	    noc2 = (noc * (noc+1)) / 2;
	    break;
	case M_E2:
	    nor2 = (nor * (nor-1)) / 2;
	    noc2 = (noc * (noc-1)) / 2;
	    break;
	case M_E3:
	    nor2 = nor*(nor-1)/2*(nor-2)/3;
	    noc2 = noc*(noc-1)/2*(noc-2)/3;
	    break;
	case M_E4:
	    nor2 = nor*(nor-1)/2*(nor-2)/3*(nor-3)/4;
	    noc2 = noc*(noc-1)/2*(noc-2)/3*(noc-3)/4;
	    break;
	case M_M3:
	    noc2 = nor2 = nor*(nor-1)*(nor+1)/3;
	    break;
	default:
	    errexit(ERR_BADUSAGE,pinfo.name);
    }
    if (nor2 < 0 || noc2 <= 0)
    {
	fprintf(stderr,"%s: Matrix too small\n",pinfo.name);
	exit(EXIT_ERR);
    }


    /* Allocate the output buffer
       -------------------------- */
    if (fl >= 2)
    {
	if (MSG0) printf("Output is %ld x %ld\n",nor2,noc2);
	fflush(stdout);
	zsetlen(noc2);
	m2 = zalloc((long)1);
    }
    else
    {	
	if (MSG0) printf("Output has degree %ld\n",nor2);
	fflush(stdout);
	noc2 = 1;
    	m2 = (PTR) Smalloc(sizeof(long)*(size_t)nor2);
	if (m2 == NULL) errexit(ERR_NOMEM,pinfo.name);
    }
}



/* ------------------------------------------------------------------
   zs2() - Symmetric square
   ------------------------------------------------------------------ */

static void zs2()

{   long i1, i2, j1, j2, j3;
    FEL f11, f12, f21, f22;
    FEL w1,w2,f1,f2;

    MESSAGE(1,("Mode S2, part 1\n"));
    for (i1 = 1; i1 < nor; ++i1)
    {	for (i2 = i1 + 1; i2 <= nor; ++i2)
	{   zmulrow(m2,F_ZERO);
	    j3 = 0;
	    for (j1 = 1; j1 < noc; ++ j1)
	    {	f11 = zextract(row[i1],j1);
		f21 = zextract(row[i2],j1);
		for (j2 = j1+1; j2 <= noc; ++j2)
		{   ++j3;
		    f12 = zextract(row[i1],j2);
		    f22 = zextract(row[i2],j2);
		    w1 = zmul(f11,f22);
		    w2 = zmul(f12,f21);
		    zinsert(m2,j3,zadd(w1,w2));
		}
	    }
	    for (j2 = 1; j2 <= noc; ++j2)
	    {	f1 = zextract(row[i1],j2);
		f2 = zextract(row[i2],j2);
		++j3;
		zinsert(m2,j3,zmul(f1,f2));
	    }
	    zwritevec(ofile,m2,1);
	}
    }

    MESSAGE(1,("Mode S2, part 2\n"));
    for (i1 = 1; i1 <= nor; ++i1)
    {	j3 = 0;
	zmulrow(m2,F_ZERO);
	for (j1 = 1; j1 < noc; ++j1)
	{   f1 = zextract(row[i1],j1);
	    for (j2 = j1+1; j2 <= noc; ++j2)
	    {	++j3;
		f2 = zextract(row[i1],j2);
		w2 = zmul(f1,f2);
		zinsert(m2,j3,zadd(w2,w2));
	    }
	}
	for (j2 = 1; j2 <= noc; ++j2)
	{   f1 = zextract(row[i1],j2);
	    ++j3;
	    zinsert(m2,j3,zmul(f1,f1));
	}
	zwritevec(ofile,m2,1);
    }
}


/* ------------------------------------------------------------------
   zs2p() - Antisymmetric square (permutations)
   ------------------------------------------------------------------ */

static long maps2(i,k)
long i, k;

{	if (i <= k)
		return ((k*(k-1))/2 + i);
	else
		return ((i*(i-1))/2 + k);
}



static void zs2p()

{	long *p1, *p2;
	register long i, k;

	p1 = (long *) m1 - 1;
	p2 = (long *) m2 - 1;

	for (i = 1; i <= nor; ++i)
	{	for (k = 1; k <= i; ++k)
		{	p2[maps2(i,k)] = maps2(p1[i],p1[k]);
		}
	}
	zwritelong(ofile,(long*) m2,nor2);
}



/* ------------------------------------------------------------------
   ze2p() - Antisymmetric square (permutations)
   ------------------------------------------------------------------ */

static long mape2(i,k)
long i, k;

{	if (i < k)
		return (((k-1)*(k-2))/2 + i);
	else
		return (((i-1)*(i-2))/2 + k);
}



static void ze2p()

{	long *p1, *p2;
	register long i, k;

	p1 = (long *) m1 - 1;
	p2 = (long *) m2 - 1;

	for (i = 1; i <= nor; ++i)
	{	for (k = 1; k < i; ++k)
		{	p2[mape2(i,k)] = mape2(p1[i],p1[k]);
		}
	}
	zwritelong(ofile,(long*) m2,nor2);
}




/* ------------------------------------------------------------------
   ze2() - Antisymmetric square
   ------------------------------------------------------------------ */

static void ze2()

{
    long i1, i2, j1, j2, j3;
    FEL f12, f22, w1, w2, w3;

    for (i1 = 1; i1 < nor; ++i1)
    {	
        MESSAGE(1,("i1 = %ld\n",i1)); 
	for (i2 = i1+1; i2 <= nor; ++i2)
	{
            MESSAGE(2,("i2 = %ld\n",i2)); 
	    zmulrow(m2,F_ZERO);
	    j3 = 1;
	    for (j1 = 1; j1 < noc; ++j1)
	    {	register FEL f11, f21;
		f11 = zextract(row[i1],j1);
		f21 = zextract(row[i2],j1);
		for (j2 = j1+1; j2 <= noc; ++j2)
		{
		    f12 = zextract(row[i1],j2);
		    f22 = zextract(row[i2],j2);
		    w1 = zmul(f11,f22);
		    w2 = zmul(f12,f21);
		    w3 = zsub(w1,w2);
		    zinsert(m2,j3,w3);
		    ++j3;
		}
	    }
	    zwritevec(ofile,m2,1);
	}
    }
}


/* ------------------------------------------------------------------
   ze3() - Antisymmetric cube
   ------------------------------------------------------------------ */

static void ze3()

{   FEL f11,f12,f13,f21,f22,f23,f31,f32,f33;
    FEL e,g12,g13,g23;
    long i1, i2, i3, j1, j2, j3, jins;

    for (i1 = 1; i1 < nor-1; ++i1)
    {	
   	MESSAGE(1,("i1 = %ld\n",i1)); 
	for (i2 = i1+1; i2 < nor; ++i2)
	{   
   	    MESSAGE(2,("i2 = %ld\n",i2)); 
	    for (i3 = i2+1; i3 <= nor; ++i3)
	    {	
   	       	MESSAGE(3,("i3 = %ld\n",i3)); 
		zmulrow(m2,F_ZERO);
		jins = 1;
		for (j1 = 1; j1 < noc-1; ++j1)
		{   f11 = zextract(row[i1],j1);
		    f21 = zextract(row[i2],j1);
		    f31 = zextract(row[i3],j1);
		    for (j2 = j1+1; j2 < noc; ++j2)
		    {	f12 = zextract(row[i1],j2);
			f22 = zextract(row[i2],j2);
			f32 = zextract(row[i3],j2);
			g12 = zsub(zmul(f11,f22),zmul(f21,f12));
			g13 = zsub(zmul(f31,f12),zmul(f11,f32));
			g23 = zsub(zmul(f21,f32),zmul(f31,f22));
			for (j3 = j2+1; j3 <= noc; ++j3)
			{   f13 =zextract(row[i1],j3);
			    f23 =zextract(row[i2],j3);
			    f33 =zextract(row[i3],j3);
			    e = zadd(zadd(zmul(g12,f33),zmul(g13,f23)),
			    	zmul(g23,f13));
			    zinsert(m2,jins,e);
			    ++jins;
			}
		    }
		}
		zwritevec(ofile,m2,1);
	    }
	}
    }
}


/* ------------------------------------------------------------------
   ze3p() - Antisymmetric cube (permutations)
   ------------------------------------------------------------------ */

#define SWAP(x,y) {tmp=x; x=y; y=tmp;}

static long mape3(i,k,l)
long i, k, l;

{	register long tmp;
	if (i < k) SWAP(i,k);
	if (i < l) SWAP(i,l);
	if (k < l) SWAP(k,l);
	return ((i-1)*(i-2)/2*(i-3)/3 + (k-1)*(k-2)/2 + l);
}



static void ze3p()

{	long *p1, *p2;
	register long i, k, l;

	p1 = (long *) m1 - 1;
	p2 = (long *) m2 - 1;

	for (i = 3; i <= nor; ++i)
	{	for (k = 2; k < i; ++k)
			for (l = 1; l < k; ++l)
				p2[mape3(i,k,l)] =
					mape3(p1[i],p1[k],p1[l]);
	}
	zwritelong(ofile,(long*) m2,nor2);
}



/* ------------------------------------------------------------------
   ze4() - Antisymmetric fourth power
   ------------------------------------------------------------------ */

static void ze4()

{   FEL f11,f12,f13,f14,f21,f22,f23,f24,f31,f32,f33,f34,f41,f42,f43,f44;
    FEL e,g12,g13,g14,g23,g24,g34,g123,g124,g134,g234;
    long i1, i2, i3, i4, j1, j2, j3, j4, jins;

    for (i1 = 1; i1 < nor-1; ++i1)
    {	
        MESSAGE(1,("i1 = %ld\n",i1)); 
	for (i2 = i1+1; i2 < nor; ++i2)
	{   
            MESSAGE(2,("i2 = %ld\n",i2)); 
	    for (i3 = i2+1; i3 <= nor; ++i3)
	    {   
            	MESSAGE(3,("i3 = %ld\n",i3)); 
		for (i4 = i3+1; i4 <= nor; ++i4)
	    	{   zmulrow(m2,F_ZERO);
		    jins = 1;

		    for (j1 = 1; j1 < noc-1; ++j1)
		    {   f11 = zextract(row[i1],j1);
		    	f21 = zextract(row[i2],j1);
		    	f31 = zextract(row[i3],j1);
		    	f41 = zextract(row[i4],j1);

		    	for (j2 = j1+1; j2 < noc; ++j2)
		    	{   f12 = zextract(row[i1],j2);
			    f22 = zextract(row[i2],j2);
			    f32 = zextract(row[i3],j2);
			    f42 = zextract(row[i4],j2);

			    g12 = zsub(zmul(f11,f22),zmul(f21,f12));
			    g13 = zsub(zmul(f11,f32),zmul(f31,f12));
			    g14 = zsub(zmul(f11,f42),zmul(f41,f12));
			    g23 = zsub(zmul(f21,f32),zmul(f31,f22));
			    g24 = zsub(zmul(f21,f42),zmul(f41,f22));
			    g34 = zsub(zmul(f31,f42),zmul(f41,f32));

			    for (j3 = j2+1; j3 <= noc; ++j3)
			    {   f13 =zextract(row[i1],j3);
			    	f23 =zextract(row[i2],j3);
			    	f33 =zextract(row[i3],j3);
			    	f43 =zextract(row[i4],j3);

				g123 = zmul(f13,g23);
				g123 = zsub(g123,zmul(f23,g13));
				g123 = zadd(g123,zmul(f33,g12));
				g124 = zmul(f13,g24);
				g124 = zsub(g124,zmul(f23,g14));
				g124 = zadd(g124,zmul(f43,g12));
				g134 = zmul(f13,g34);
				g134 = zsub(g134,zmul(f33,g14));
				g134 = zadd(g134,zmul(f43,g13));
				g234 = zmul(f23,g34);
				g234 = zsub(g234,zmul(f33,g24));
				g234 = zadd(g234,zmul(f43,g23));

			    	for (j4 = j3+1; j4 <= noc; ++j4)
				{   f14 =zextract(row[i1],j4);
				    f24 =zextract(row[i2],j4);
				    f34 =zextract(row[i3],j4);
				    f44 =zextract(row[i4],j4);

				    e = zmul(f24,g134);
				    e = zsub(e,zmul(f14,g234));
				    e = zadd(e,zmul(f44,g123));
				    e = zsub(e,zmul(f34,g124));
			    	    zinsert(m2,jins,e);
			    	    ++jins;
				}
			    }
			}
		    }
		    zwritevec(ofile,m2,1);
		}
	    }
	}
    }
}


/* ------------------------------------------------------------------
   zm3() - (2,1) part of tensor cube
   ------------------------------------------------------------------ */

#if defined(__STDC__)
static void addtorow(long i1, long i2, long i3, FEL f)
#else
static void addtorow(i1,i2,i3,f) long i1, i2, i3;
FEL f;
#endif

{	long pos;

/*printf("addtorow(%ld,%ld,%ld)\n",i1,i2,i3);fflush(stdout);*/
	if (i1 >= i2 || i1 > i3)
	{	fprintf(stderr,"internal error (%s, line %d)\n",
			__FILE__,__LINE__);
	}

	pos = (i2 - i1 - 1) * (nor - i1 + 1) + i3 - i1 + 1;
	for (; i1 > 1; --i1)
		pos += (nor-i1+1)*(nor-i1+2);
	zinsert(m2,pos,zadd(zextract(m2,pos),f));
}



static void mapvector(i1,i2,i3)
long i1, i2, i3;

{   long k1, k2, k3;
    FEL f1, f2, f3;
    long l1, l2;

    for (k1 = 1; k1 <= nor; ++k1)
    {   f1 = zextract(row[i1],k1);
	for (k2 = 1; k2 <= nor; ++k2)
	{   if (k2 == k1) continue;
	    f2 = zmul(f1,zextract(row[i2],k2));
	    if (k2 < k1)
	    {	l1 = k2;
		l2 = k1;
		f2 = zadd(F_ZERO,f2);
	    }
	    else
	    {	l1 = k1;
		l2 = k2;
	    }
	    for (k3 = 1; k3 <= nor; ++k3)
	    {	f3 = zmul(f2,zextract(row[i3],k3));
		if (k3 < l1)
		{	f3 = zneg(f3);
			addtorow(k3,l1,l2,f3);
			addtorow(k3,l2,l1,f3);
		}
		else
		{	addtorow(l1,l2,k3,f3);
		}
	    }
	}
    }
}


static void zm3()

{
    long i,j,k;

    for (i = 1; i < nor; ++i)
    {	
	for (j = i + 1; j <= nor; ++j)
	{
	    for (k = i; k <= nor; ++k)
	    {	
		zmulrow(m2,F_ZERO);
		mapvector(i,j,k);
	    	zwritevec(ofile,m2,1);
	    }
	}
    }
}



/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char *argv[])

{
    mtxinit();
    parseargs(argc,argv);
    init();
    if (fl >= 2)
    {
	if ((ofile = zwritehdr(oname,fl,nor2,noc2)) == NULL)
	    errexit(-1,oname);;
    }
    else
    {
	if ((ofile = zwritehdr(oname,fl,nor2,(long)1)) == NULL)
	    errexit(-1,oname);
    }
    switch (mode)
    {	case M_S2:
		if (fl >= 2)
			zs2();
		else
			zs2p();
		break;
	case M_E2:
		if (fl >= 2)
			ze2();
		else
			ze2p();
		break;
	case M_E3:
		if (fl >= 2)
			ze3();
		else
			ze3p();
		break;
	case M_E4:
		ze4();
		break;
	case M_M3:
		zm3();
		break;
	default:
		errexit(ERR_BADUSAGE,pinfo.name);
		break;
	}
	fclose(ofile);
	return EXIT_OK;
}


