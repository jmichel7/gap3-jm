/* ========================== C MeatAxe =============================
   ztm.c - Multiply tensor product of two matrices to a vector

   (C) Copyright 1994 Christoph Jansen, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <cjansen@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: ztm.c,v 1.2 1997/09/11 15:44:36 gap Exp $
 *
 * $Log: ztm.c,v $
 * Revision 1.2  1997/09/11 15:44:36  gap
 * New version 2.2.3. AH
 *
 * Revision 1.3  1995/05/12  10:28:08  mringe
 * Initialisiere mat1 und mat2 (Compiler warning)
 *
 * Revision 1.2  1994/12/17  20:33:08  mringe
 * Neu: Message() und ErrorExit().
 *
 * Revision 1.1  1994/12/14  14:58:24  mringe
 * Initial revision
 *
 */

#include <stdlib.h>
#include <string.h>
#include "meataxe.h"


/* ------------------------------------------------------------------
   Variables
   ------------------------------------------------------------------ */

static char *helptext[] = {
"SYNTAX",
"    ztm <A> <B> <Vectors> <Result>",
"",
"FILES",
"    <A> and <B> are square matrices of dimension m and n, resp.",
"    <Vectors> is a matrix with mn columns and any number of rows.",
"    The result has the same size as <Vectors>.\n",
NULL};

static proginfo_t pinfo =
   { "ztm", "Tensor Multiply", "$Revision: 1.2 $", helptext };

static char *AName, *BName;	/* Matrices */
static char *VName, *RName;	/* Vectors, Result */



/* ------------------------------------------------------------------
   VectorToMatrix() - Convert vector to matrix
   ------------------------------------------------------------------ */

matrix_t *VectorToMatrix(PTR vec, long fl, long nrows, long ncol)

{
    matrix_t *mat;
    PTR d;
    long i, k;

    mat = matalloc(fl, nrows, ncol);
    d = mat->d;
    zsetlen(ncol);
    for (i = 0; i < nrows; i++) 
    {
        for (k = 1; k <= ncol; k++) 
        {
            FEL f = zextract(vec,i*ncol+k);
            zinsert(d,k,f);
        }
        zadvance(&d,(long)1);
    }
    return mat;
}




/* ------------------------------------------------------------------
   MatrixToVector() - Convert matrix to vector
   ------------------------------------------------------------------ */

void MatrixToVector(matrix_t *mat, PTR vec)

{
    PTR y;
    long i, k, nrows, ncol, l;

    ncol = mat->noc;
    nrows = mat->nor;

    /* Clear the vector 
       ---------------- */
    zsetlen(ncol * nrows);
    zmulrow(vec,F_ZERO);

    /* Convert
       ------- */
    zsetlen(ncol);
    y = mat->d;
    l = 1;
    for (i=0; i < nrows; i++) 
    {
        for (k=1; k <= ncol; k++) 
        {
            FEL f = zextract(y,k);
            zinsert(vec,l++,f);
        }
        zadvance(&y,(long)1);
    }
}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(int argc, char *argv[])

{  
    matrix_t *mat1=NULL, *mat2=NULL, *mat1tr;
    FILE *f, *vecfile, *resultfile;
    PTR tmp;
    long i, ssdim, fl, nocol;

    /* Initialize everything, get options
       ---------------------------------- */
    mtxinit();
    initargs(argc, argv, &pinfo);
    while (zgetopt("") != OPT_END);
    if ((i = opt_ind) != argc-4)
	ErrorExit("%s: %E",pinfo.name,ERR_BADUSAGE);
    AName = argv[opt_ind];
    BName = argv[opt_ind+1];
    VName = argv[opt_ind+2];
    RName = argv[opt_ind+3];

    /* Read matrices
       ------------- */
    if ((f = SFOpen(AName,FM_READ)) == NULL ||
        (mat1 = matread(f)) == NULL)
	errexit(-1,AName);
    fclose(f);
    if ((f = SFOpen(BName,FM_READ)) == NULL ||
        (mat2 = matread(f)) == NULL)
	errexit(-1,BName);
    fclose(f);

    /* Test if matrices are square and over the same field
       --------------------------------------------------- */
    if (mat1->nor != mat1->noc)
	ErrorExit("%s: %E",AName,ERR_NOTSQUARE);
    if (mat2->nor != mat2->noc)
	ErrorExit("%s: %E",BName,ERR_NOTSQUARE);
    if (mat1->fl != mat2->fl)
	ErrorExit("%s and %s: %E",AName,BName,ERR_INCOMPAT);
    MESSAGE(1,("%s: %ldx%ld matrix over GF(%ld)\n",
	AName,mat1->nor,mat1->noc,mat1->fl));
    MESSAGE(1,("%s: %ldx%ld matrix over GF(%ld)\n",
	BName,mat2->nor,mat2->noc,mat2->fl));

    /* Open the <Vectors> file and check the header
       -------------------------------------------- */
    if ((vecfile = zreadhdr(VName,&fl,&ssdim,&nocol)) == NULL)
	errexit(-1,VName);
    if (fl < 2) errexit(ERR_NOTMATRIX,VName);
    MESSAGE(1,("%s: %ldx%ld matrix over GF(%ld)\n",VName,ssdim,
	nocol,fl));
    if (fl != mat1->fl || nocol != mat1->nor * mat2->nor)
	ErrorExit("%s x %s and %s: %E",AName,BName,VName,ERR_INCOMPAT);
    

    /* Write beginning of output file
       ------------------------------ */
    if ((resultfile = zwritehdr(RName,mat1->fl,ssdim,nocol)) == NULL)
        errexit(-1,RName);

    /* Transpose first matrix 
       ---------------------- */
    mat1tr = mattr(mat1);
    matfree(mat1);

    /* Allocate buffer for one vector
       ------------------------------ */
    zsetlen(nocol);
    tmp = zalloc((long)1);

    /* Process <Vectors> one by one
       ---------------------------- */
    for ( i=1; i <= ssdim; ++i)
    {
	matrix_t *mat3, *newmat;

        /* Read on evector and convert to matrix
           ------------------------------------- */
        zsetlen(nocol);
        zreadvec(vecfile,tmp,1);

        /* Turn vector into matrix
           ----------------------- */
        mat3 = VectorToMatrix(tmp,mat1tr->fl,mat1tr->nor,mat2->nor);
       
        /* Multiply
           -------- */
	newmat = matdup(mat1tr);
	matmul(newmat,mat3);
	matmul(newmat,mat2);

        /* Turn matrix into vector and write out
           ------------------------------------- */
        MatrixToVector(newmat,tmp);
        zsetlen(nocol);
        zwritevec(resultfile,tmp,1);

        /* Free memory
           ----------- */
        matfree(mat3);
        matfree(newmat);
    }
   
    /* Close files
       ---------- */
    fclose(vecfile);
    fclose(resultfile);
    return 0;
}


