/* ========================== C MeatAxe =============================
   mkinc.c - Calculate the incidences between mountains.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: mkinc.c,v 1.2 1997/09/11 15:43:09 gap Exp $
 *
 * $Log: mkinc.c,v $
 * Revision 1.2  1997/09/11 15:43:09  gap
 * New version 2.2.3. AH
 *
 * Revision 2.22  1995/06/27  08:13:05  mringe
 * GAP output (written by Max Bisani).
 *
 * Revision 2.21  1995/02/08  10:01:14  mringe
 * PL_ entfernt.
 *
 * Revision 2.20  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.19  1994/08/22  12:08:03  mringe
 * *** empty log message ***
 *
 * Revision 2.18  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.17  1994/07/07  06:57:36  mringe
 * mount in mountlist umbenannt.
 *
 * Revision 2.16  1994/07/07  05:35:16  mringe
 * *** empty log message ***
 *
 * Revision 2.15  1994/06/16  19:20:02  mringe
 * *** empty log message ***
 *
 * Revision 2.14  1994/05/18  09:48:40  mringe
 * matspin(), matquot(), issubspace() durch die
 * entsprechenden neuen Libraryfunktionen ersetzt.
 *
 * Revision 2.13  1994/03/26  06:34:04  mringe
 * basename umbenannt wg. Namenskonflikt.
 *
 * Revision 2.12  1994/03/15  10:37:59  mringe
 * Allokiere proj[] dynamisch.
 *
 * Revision 2.11  1994/03/13  13:28:13  mringe
 * Benutze zreadlong()/zwritelong()
 *
 * Revision 2.10  1994/02/21  12:41:41  mringe
 * mtxini() am ANFANG!
 *
 * Revision 2.9  1994/02/15  10:28:33  mringe
 * MSDOS_BCC entfernt.
 *
 * Revision 2.8  1994/02/15  09:42:27  mringe
 * Kleiner Bug (%d statt %ld).
 *
 * Revision 2.7  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.6  1994/02/12  04:10:13  mringe
 * UMFANGREICHE AENDERUNGEN AN VIELEN DATENTYPEN.
 *
 * Revision 2.5  1993/12/13  08:25:53  mringe
 * Reihenfolge der Fkt.-argumente vereinheitlicht.
 *
 * Revision 2.4  1993/12/08  11:48:50  mringe
 * Compiler warnings.
 *
 * Revision 2.3  1993/12/08  11:33:02  mringe
 * Neue CPU time - Funktionen.
 *
 * Revision 2.2  1993/10/22  17:03:01  mringe
 * Neues Numerierungsschema.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.18  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.17  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.16  1993/10/02  16:23:02  mringe
 * matread() und matwrite() in matload() bzw. matsave() umbenannt.
 *
 * Revision 1.15  1993/09/21  07:15:15  mringe
 * Optionen -GQVT, help verbessert.
 *
 * Revision 1.14  1993/08/10  14:29:19  mringe
 * Include string.h
 *
 * Revision 1.13  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.12  1993/07/17  19:13:05  mringe
 * Aenderungen fuer Borland C.
 *
 * Revision 1.11  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.10  1993/07/01  22:39:25  mringe
 * WICHTIGE AENDERUNG: Rechne nicht mehr Aeq.klassen von zyklischen
 * Teilmoduln im kondensierten aus, die alle zum gleichen Mountain
 * aufspinnen. Statt dessen fasse alle Vektoren zu einer Klasse
 * zusammen, die in einem kondensierten Mountain liegen. Das sind
 * keine Aequ.klassen mehr, da sie sich ueberlappen koennen.
 *
 * Revision 1.9  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.9  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.8  1993/02/15  13:49:24  mringe
 * Funktionen aus cyclic.c -> yyy-Lib verschoben.
 *
 * Revision 1.7  1993/02/10  19:40:54  mringe
 * Libraries angelegt (YYY und ZZZ).
 *
 * Revision 1.6  1992/10/01  13:50:10  mringe
 * Header eingef"ugt.
 *
 * Revision 1.5  1992/07/22  07:10:30  mringe
 * Changed 'global.h' to 'lattice.h'
 *
 * Revision 1.4  1992/07/15  09:25:55  mringe
 * Some minor changes.
 *
 * Revision 1.3  1992/07/14  10:03:17  mringe
 * step2(): Namen der Irred. nur einmal ausgeben
 *
 * Revision 1.2  1992/06/16  15:15:19  mringe
 * Bug in step2() behoben.
 *
 * Revision 1.1  1992/05/25  17:31:56  mringe
 * Initial revision
 *
 */

#include <string.h>
#include <stdlib.h>



#include "meataxe.h"
#include "lattice.h"
#include "files.h"



/* -----------------------------------------------------------------
   Function prototypes
   ----------------------------------------------------------------- */

static void init(void);
static void writemountains(void);
static int newmountain(matrix_t *vec, int cf);
static void makeclass(int mnt, int cf, matrix_t *vectors);
static void step1(void);
static void writeresult(void);
static void step2(void);


/* -----------------------------------------------------------------
   Global data
   ----------------------------------------------------------------- */

matrix_t *gen[MAXGEN];			/* Generators */
matrix_t *bild[MAXCF];			/* Image of peak word (gkond) */
long *piv;				/* Pivot table */
int nmount = 0;				/* Number of mountains */
matrix_t * mountlist[MAXCYCL];		/* List of mountains */
long mountdim[MAXCYCL];			/* Dim. of mountains */
matrix_t **proj[MAXCYCL];		/* Projections of mountains */
int moffset[MAXCF];			/* Number of first mountain */
long *class[MAXCYCL];			/* Classes of vectors */
bitstring_t * subof[MAXCYCL];		/* Incidence matrix */

static char *helptext[] = {
"SYNTAX",
"    mkinc [<Options>] <Name>",
"",
"OPTIONS",
"    -Q            Quiet, no messages.",
"    -V            Verbose, more messages.",
"    -G            GAP output (implies -Q).",
"    -T <Seconds>  Set CPU time limit",
"",
"FILES",
"    <Name><Dim><a>.<N>   i    Generators",
"    <Name><Dim><a>.<N>k  i    Kondensed generators",
"    <Name><Dim><a>.v     i    Cyclic submodules",
"    <Name><Dim><a>.im    i    Image used for kondensation",
"    <Name><Dim><a>.k     i    Unkondense matrix",
"    <Name>.v             o    Mountains",
"    <Name>.mnt           o    Dimensions of mountains and classes",
"                              of cyclic subspaces",
"    <Name>.inc           o    Incidence matrix",
NULL};

static proginfo_t pinfo =
   { "mkinc", "Mountains And Incidence Matrix",
     "$Revision: 1.2 $", helptext };


int opt_G = 0;


/* -----------------------------------------------------------------
   init() - Read generators and images of peak words
   ----------------------------------------------------------------- */

static void init()

{
    char fn[200];
    int i;

    readcfinfo();
    for (i = 0; i < ngen; ++i)
    {
	sprintf(fn,"%s.%d",cfbasename,i+1);
	gen[i] = matload(fn);
	if (gen[i] == NULL)
	{
	    perror(fn);
	    FATAL("Error reading file");
	}
    }
    for (i = 0; i < ncf; ++i)
    {
	sprintf(fn,"%s%s.im",cfbasename,cfname(i));
	bild[i] = matload(fn);
    }
    piv = (long *) malloc((size_t)(gen[0]->nor+1)*sizeof(long));
}


/* -----------------------------------------------------------------
   writemountains() - Write mountains, dimensions and classes.
   ----------------------------------------------------------------- */

static void writemountains()

{
    FILE *f;
    char fn[200];
    int i;
    long *p;

    /* Write dimensions and classes
       ---------------------------- */
    f= SFOpen(strcat(strcpy(fn,cfbasename),".mnt"),FM_CREATE|FM_TEXT);
    if (f == NULL) errexit(-1,fn);
    MESSAGE(0,("Writing dimensions and classes to %s\n",fn));
    for (i = 0; i < nmount; ++i)
    {
	fprintf(f,"%4d %4ld  ",i,mountdim[i]);
	for (p = class[i]; *p > 0; ++p)
	    fprintf(f,"%ld ",*p);
	fprintf(f,"-1\n");
    }
    fclose(f);

    /* Write mountains
       --------------- */
    strcat(strcpy(fn,cfbasename),".v");
    MESSAGE(0,("Writing mountains to %s\n",fn));
    zsetlen(gen[0]->noc);
    if ((f = zwritehdr(fn,zfl,nmount,gen[0]->noc)) == NULL)
	errexit(-1,fn);
    for (i = 0; i < nmount; ++i)
    {
	zwritevec(f,mountlist[i]->d,(long)1);
	matfree(mountlist[i]);	/* We don't need them for step 2*/
    }
    fclose(f);
}


/* -----------------------------------------------------------------
   newmountain() - Take a vector, spin it up and check if it is a
	new mountain. If yes, add it to the mountain list and
	calculate the projection on each kondensed module.
	Returns 1 if a new mountain has been found, 0 if not.
   ----------------------------------------------------------------- */

static int newmountain(vec, cf)
matrix_t *vec;
int cf;

{
    matrix_t *span, *backproj;
    int i;

    /* Spin up the vector and project back onto
       the kondensed module where it came from.
       ----------------------------------------- */
    spinup(vec,ngen,gen,SPL_SEED_FIRST,&span);
    MESSAGE(2,("Next vector spins up to %ld\n",span->nor));
    backproj = quotproj(bild[cf],span);

    /* Check if it is a new mountain
       ----------------------------- */
    for (i = moffset[cf]; i < nmount; ++i)
    {
	if (backproj->nor == proj[i][cf]->nor &&
		spccontains(backproj,proj[i][cf]))
	    break;
    }

    /* If it is new, add it to the list and calculate
       the other projections. Otherwise just forget it.
       ------------------------------------------------ */
    if (i == nmount)
    {
	int k;

	if (nmount >= MAXCYCL)
	    FATAL("TOO MANY MOUNTAINS, INCREASE MAXCYCL");
	proj[nmount] = (matrix_t **) malloc((size_t)ncf *
		sizeof(matrix_t *));

    	MESSAGE(2,("New Mountain %d\n",(int)i));
	for (k = 0; k < ncf; ++k)
	{
    	    MESSAGE(3,("Projecting on %d\n",k));
	    if (k == cf)
	    	proj[nmount][cf] = backproj;
	    else
	    {
		proj[nmount][k] = quotproj(bild[k],span);
	    }
	}
	mountlist[nmount] = vec;
	mountdim[nmount] = span->nor;
	++nmount;
	matfree(span);


	return 1;
    }
    else
    {
	matfree(backproj);
	matfree(span);
	return 0;
    }
}


/* -----------------------------------------------------------------
   makeclass() - Find all cyclic subspaces of the projection of a
   mountain onto `its' irreducible.
   ----------------------------------------------------------------- */

static void makeclass(mnt,cf,vectors)
int mnt;
int cf;
matrix_t *vectors;

{
    char *tmp = (char *) malloc((size_t)vectors->nor+2);
    matrix_t *vec;
    long k, *p;
    size_t nvec;

    nvec = 0;
    MESSAGE(2,("Making equivalence class\n"));
    for (k = 1; k <= vectors->nor; ++k)
    {
	vec = matextract(vectors,k,(long)1);
	tmp[k] = 0;
        MESSAGE(3,(" %ld",k));
	if (spccontains(proj[mnt][cf],vec))
	{
	    tmp[k] = 1;
	    ++nvec;
	}
	matfree(vec);
    }
    MESSAGE(3,("\n"));

    p = class[mnt] = NALLOC(long,nvec+2);
    *p++ = nvec;
    for (k = 1; k <= vectors->nor; ++k)
	if (tmp[k]) *p++ = k;
    *p = -1;
    free(tmp);
}


/* -----------------------------------------------------------------
   step1() - Make all mountains and calculate the projections of
	mountains on kondensed modules.
   ----------------------------------------------------------------- */

static void step1()

{
    matrix_t *vectors, *vec, *uk;
    char fn[200];
    int cf;
    long i;

    MESSAGE(0,("Step 1 (Mountains)\n"));
    nmount = 0;
    for (cf = 0; cf < ncf; ++cf)	/* For each irreducible */
    {
	MESSAGE(0,("  %s%s: ",cfbasename,cfname(cf)));

	/* Read the vectors and the unkondense matrix
	   ------------------------------------------ */
	sprintf(fn,"%s%s.v",cfbasename,cfname(cf));
	vectors = matload(fn);
	sprintf(fn,"%s%s.k",cfbasename,cfname(cf));
	uk = matload(fn);

	/* Try each vector
	   --------------- */
	moffset[cf] = nmount;
	for (i = 1; i <= vectors->nor; ++i)
	{
	    if (i % 50 == 0)
		MESSAGE(1,("[%ld] ",i));
	    vec = matextract(vectors,i,(long)1);
	    matmul(vec,uk);	/* Unkondense */
	    if (newmountain(vec,cf))
		makeclass(nmount-1,cf,vectors);
	    else
	        matfree(vec);
	}
	cfinfo[cf].nmount = nmount - moffset[cf];

	matfree(vectors);
	matfree(uk);
	MESSAGE(0,("%ld mountain%s\n",cfinfo[cf].nmount,
		cfinfo[cf].nmount != 1 ? "s" : ""));

    }
    MESSAGE(0,("Total: %d mountain%s\n",nmount,nmount != 1 ? "s" : ""));
}




/* -----------------------------------------------------------------
   writeresult() - Write the incidence matrix
   ----------------------------------------------------------------- */

static void writeresult()

{
    FILE *f;
    char fn[200];
    int i;
    long l;

    /* Write the incidence matrix
       -------------------------- */
    f = SFOpen(strcat(strcpy(fn,cfbasename),".inc"),FM_CREATE);
    if (f == NULL) errexit(-1,fn);
    MESSAGE(0,("Writing incidence matrix (%s)\n",fn));
    l = (long) nmount;
    zwritelong(f,&l,1);
    for (i = 0; i < nmount; ++i)
	bs_write(f,subof[i]);
    fclose(f);

    /* Write the .cfinfo file
       ---------------------- */
    writecfinfo();
}



/* -----------------------------------------------------------------
   step2() - Calculate all incidences
   ----------------------------------------------------------------- */

static void step2()

{
    int i, k;
    int cfi, cfk;	/* Comp. factor corresponding to mountain i,j */

    MESSAGE(0,("Step 2 (Incidences)\n"));

    /* Allocate memory for the incidence matrix
       ---------------------------------------- */
    bs_setlen(nmount);
    for (i = 0; i < nmount; ++i)
	subof[i] = bs_alloc();


    /* Calculate the incidences
       ------------------------ */
    for (cfi=0, i = 0; i < nmount; ++i)
    {
	if (i == moffset[cfi])
	    MESSAGE(0,("  %s%s\n",cfbasename,cfname(cfi)));
	for (cfk=0, k = 0; k < nmount; ++k)
	{
	    if (spccontains(proj[i][cfk],proj[k][cfk]))
		bs_set(subof[k],i);
	    if (spccontains(proj[k][cfi],proj[i][cfi]))
		bs_set(subof[i],k);
	    if (cfk < ncf && k == moffset[cfk+1]-1)
		++cfk;
    	}
    	if (cfi < ncf && i == moffset[cfi+1]-1)
			++cfi;
    }
}


/* -----------------------------------------------------------------
   writemountainsGAP() - Write mountains, dimensions and classes.
   ----------------------------------------------------------------- */

static void writemountainsGAP()

{
    int i;
    long *p;

    /* Write dimensions and classes
       ---------------------------- */
    printf( "MeatAxe.Dimensions := [" ) ;
    for (i = 0; i < nmount-1 ; ++i)
	printf( "%ld,", mountdim[i] );
    printf( "%ld] ;\n", mountdim[nmount-1] ) ;

    printf( "MeatAxe.Classes := [\n" ) ;
    for (i = 0; i < nmount; ++i)
    {
        printf( "[" ) ;
	for (p = class[i] ; *p > 0 ; ++p)
	{
	    printf( "%ld", *p ) ;
	    printf( p[1] > 0 ? "," : "]" ) ;
	}
	printf( i < nmount-1 ? ",\n" : "] ;\n");
    }
}


/* -----------------------------------------------------------------
   writeresultGAP() - Write the incidence matrix
   ----------------------------------------------------------------- */

static void writeresultGAP()

{	int i,  j;

	/* Write the incidence matrix
	   -------------------------- */
	printf( "MeatAxe.NMount := %ld;\n", (long) nmount );
	printf( "MeatAxe.Incidences := [\n" );
	for (i = 0; i < nmount; ++i)
        {
	    printf( "BlistList([" );
	    for (j = 0 ; j < nmount ; j++)
      	        printf( j < (nmount - 1) ? "%s," : "%s], [1])" , 
		       bs_test(subof[i],j) ? "1" : "0" ) ;

	    printf( i < nmount-1 ? ",\n" : "] ;\n" ) ;
	}
}




/* -----------------------------------------------------------------
   main()
   ----------------------------------------------------------------- */

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
    if (opt_ind != argc-1) errexit(ERR_BADUSAGE,"mkinc");
    MESSAGE(0,("\n*** INCIDENCE MATRIX ***\n\n"));
    setbasename(argv[argc-1]);
    init();
    step1();
    writemountains();
    step2();
    writeresult();
    if (opt_G) 
    {
	writeresultGAP();
	writemountainsGAP();
    }
    if (msg_level >= 0)
    {
    	printf("\n");
    	prtimes();
    }
    return 0;
}


