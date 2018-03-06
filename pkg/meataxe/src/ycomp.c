/* ========================== C MeatAxe =============================
   ycomp.c - Compare vector spaces

   (C) Copyright 1994 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */

/* $Id: ycomp.c,v 1.2 1997/09/11 15:43:27 gap Exp $
 *
 * $Log: ycomp.c,v $
 * Revision 1.2  1997/09/11 15:43:27  gap
 * New version 2.2.3. AH
 *
 * Revision 2.5  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.4  1994/06/16  19:20:14  mringe
 * Bug in spccontains() beseitigt.
 *
 * Revision 2.3  1994/05/18  11:53:31  mringe
 * Neue Funktionen spXXXX().
 *
 * Revision 2.2  1993/12/08  11:48:50  mringe
 * Compiler warnings.
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
 * Revision 1.6  1993/10/05  19:17:50  mringe
 * *** empty log message ***
 *
 * Revision 1.5  1993/10/05  19:02:08  mringe
 * yerror eliminiert.
 *
 * Revision 1.4  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.3  1993/08/06  12:54:10  mringe
 * ERR_NOECH in ERR_NOIECH umbenannt.
 *
 * Revision 1.2  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.1  1993/02/15  13:50:28  mringe
 * Initial revision
 *
 */

#include <stdlib.h>
#include "meataxe.h"


static long *piv = NULL;
static long maxpiv;


/* ------------------------------------------------------------------
   spccomp() - Compare two vector spaces. If n == 0, m1 and m2 are
   full vector spaces. The spaces must be in echelon form.
   If n > 0, the first n rows of each matrix are considered to
   generate the whole space.

   Return value:
          -1  m1 < m2
	   0  m1 = m2
	   1  m1 > m2
	   9  else
	  -9  error
   ------------------------------------------------------------------ */

int spccomp(m1,m2,n)
matrix_t *m1, *m2;
long n;

{
    PTR tmp, y;
    long i;
    FEL f;
    long dim1 = m1->nor, dim2 = m2->nor;  /* Dimensions */
    long ngen1 = dim1, ngen2 = dim2;	  /* Number of generators */

    /* Check the arguments
       ------------------- */
    if (n < 0) MTXFAIL(ERR_BADARG,-9);
    if (m1->noc != m2->noc || m1->fl != m2->fl)
	MTXFAIL(ERR_INCOMPAT,-9);
    zsetfield(m1->fl);
    zsetlen(m1->noc);

    /* Set ngen1 and ngen2
       ------------------- */
    if (n > 0)
    {
	if (n < dim1) ngen1 = n;
	if (n < dim2) ngen2 = n;
    }

    /* Check if m1 <= m2. We don't try this, of
       course, if dim(m1) > dim(m2).
       ---------------------------------------- */
    if (dim1 <= dim2)
    {
	if (piv == NULL || maxpiv < dim2)
	{
	    if (piv != NULL) free(piv);
	    piv = NALLOC(long,dim2+1);
	    maxpiv = dim2;
	}
	if (zmkpivot(m2->d,dim2,piv) != 0) MTXFAIL(ERR_NOTECH,-9);
    	tmp = zalloc((long)1);
	for (i = ngen1, y = m1->d; i > 0; --i, zadvance(&y,1))
	{
	    zmoverow(tmp,y);
	    zcleanrow(tmp,m2->d,dim2,piv);
	    if (zfindpiv(tmp,&f) > 0) break;
	}
	free(tmp);
	if (i == 0)
	{
	    if (dim1 == dim2) return 0;	/* m1 = m2 */
	    return -1;			/* m1 <= m2 */
	}
    }

    /* Now, check if m1 > m2.
       ---------------------- */
    if (dim1 > dim2)
    {
	if (piv == NULL || maxpiv < dim1)
	{
	    if (piv != NULL) free(piv);
	    piv = NALLOC(long,dim1+1);
	    maxpiv = dim1;
	}
	if (zmkpivot(m1->d,dim1,piv) != 0) MTXFAIL(ERR_NOTECH,-9);
    	tmp = zalloc((long)1);
	for (i = ngen2, y = m2->d; i > 0; --i, zadvance(&y,1))
	{
	    zmoverow(tmp,y);
	    zcleanrow(tmp,m1->d,dim1,piv);
	    if (zfindpiv(tmp,&f) > 0) break;
	}
	free(tmp);
	if (i == 0)
	    return 1;			/* m1 > m2 */
    }

    /* Neither m1<=m2 nor m1>=m2.
       -------------------------- */
    return 9;	/* */
}


/* ------------------------------------------------------------------
   spcequal() - Check if m1==m2. m1 must be in echelon form, m2
	may be any matrix.
   Return value:  0  m1!=m2
		  1  m1=m2
		  -1 error
   ------------------------------------------------------------------ */

int spcequal(m1,m2)
matrix_t *m1, *m2;

{
    PTR tmp, y;
    long dim, i;
    FEL f;

    /* Check the arguments
       ------------------- */
    if (m1->noc != m2->noc || m1->fl != m2->fl)
	MTXFAIL(ERR_INCOMPAT,-1);
    zsetfield(m1->fl);
    zsetlen(m1->noc);

    /* Compare the dimensions
       ---------------------- */
    if (m1->nor != m2->nor)
	return 0;

    /* Make the pivot table for m1
       --------------------------- */
    dim = m1->nor;
    if (piv == NULL || maxpiv < dim)
    {
	if (piv != NULL) free(piv);
	piv = NALLOC(long,dim+1);
	maxpiv = dim;
    }
    if (zmkpivot(m1->d,dim,piv) != 0)
	MTXFAIL(ERR_NOTECH,-1);

    /* Clean the rows of m2 with m1
       ---------------------------- */
    tmp = zalloc((long)1);
    for (i = dim, y = m2->d; i > 0; --i, zadvance(&y,1))
    {
	zmoverow(tmp,y);
	zcleanrow(tmp,m1->d,dim,piv);
	if (zfindpiv(tmp,&f) > 0) break;
    }
    free(tmp);
    if (i == 0)
	return 1;			/* m1 = m2 */
    else
	return 0;			/* m1 != m2 */
}


/* ------------------------------------------------------------------
   spccontains() - Check if m1>=m2. m1 must be in echelon form, m2
	may be any matrix.
   Return value:  0  no
		  1  yes
		  -1 error
   ------------------------------------------------------------------ */

int spccontains(m1,m2)
matrix_t *m1, *m2;

{
    PTR tmp, y;
    long dim, i;
    FEL f;

    /* Check the arguments
       ------------------- */
    if (m1->noc != m2->noc || m1->fl != m2->fl)
	MTXFAIL(ERR_INCOMPAT,-1);
    zsetfield(m1->fl);
    zsetlen(m1->noc);

    /* Compare the dimensions
       ---------------------- */
    if (m1->nor < m2->nor)
	return 0;

    /* Make the pivot table for m1
       --------------------------- */
    dim = m1->nor;
    if (piv == NULL || maxpiv < dim)
    {
	if (piv != NULL) free(piv);
	piv = NALLOC(long,dim+1);
	maxpiv = dim;
    }
    if (zmkpivot(m1->d,dim,piv) != 0)
	MTXFAIL(ERR_NOTECH,-1);

    /* Clean the rows of m2 with m1
       ---------------------------- */
    tmp = zalloc((long)1);
    for (i = m2->nor, y = m2->d; i > 0; --i, zadvance(&y,1))
    {
	zmoverow(tmp,y);
	zcleanrow(tmp,m1->d,dim,piv);
	if (zfindpiv(tmp,&f) > 0) break;
    }
    free(tmp);
    if (i == 0)
	return 1;			/* m1 >= m2 */
    else
	return 0;			/* 'no' */
}

