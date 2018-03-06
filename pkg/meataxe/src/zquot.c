/* ========================== C MeatAxe =============================
   zqt.c - Projection on quotient.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zquot.c,v 1.2 1997/09/11 15:44:21 gap Exp $
 *
 * $Log: zquot.c,v $
 * Revision 1.2  1997/09/11 15:44:21  gap
 * New version 2.2.3. AH
 *
 * Revision 2.3  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.2  1994/05/20  09:06:16  mringe
 * Benutze Smalloc().
 *
 * Revision 2.1  1994/05/10  16:29:15  mringe
 * NEU: Der Speicher muss jetzt vom aufrufenden Programm
 * allokiert werden.
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
 * Revision 1.6  1993/10/05  23:35:02  mringe
 * zerrno eliminiert.
 *
 * Revision 1.5  1993/08/12  16:14:23  mringe
 * Include <string.h>.
 *
 * Revision 1.4  1993/08/06  15:09:31  mringe
 * ERR_NOECH in ERR_NOTECH umbenannt.
 *
 * Revision 1.3  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.2  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.1  1993/02/10  19:40:54  mringe
 * Initial revision
 *
 * Revision 1.1  1993/01/28  13:42:30  mringe
 * Initial revision
 *
 */

#include <string.h>
#include <stdlib.h>
#include "meataxe.h"




/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static long *mypiv;
static int ownpiv = 0;
static long *locn = NULL;
static long vdim, qdim, sdim;
static PTR mysubsp;


/* ------------------------------------------------------------------
   zquotinit() - Initialize. Finds the insignificant rows.

   Arguments: <subspace> must be a matrix in semi-echelon form with
	<znoc> columns and <dim> rows. <piv> must point to a pivot
	table for this matrix. If <piv>==NULL, a pivot table is
	built by a call to zmkpivot().

   Return value: Returns 0 on success, -1 on error.
   ------------------------------------------------------------------ */

int zquotinit(subspace,dim,piv)
PTR subspace;
long dim;
long *piv;

{
    long i;
    char *ispiv;

    if (locn != NULL)
    {
	free(locn);
	locn = NULL;
    }
    if (ownpiv)
    {
	free(mypiv);
	ownpiv = 0;
    }

    if (piv != NULL)
    {
    	mypiv = piv;
	ownpiv = 0;
    }
    else
    {
	mypiv = NALLOC(long,dim+1);
	if (mypiv == NULL)
	{
	    MTXFAIL(ERR_NOMEM,-1);
	}
	ownpiv = 1;
	if (zmkpivot(subspace,dim,mypiv) != 0)
	    return -1;
    }
    vdim = znoc;
    sdim = dim;
    mysubsp = subspace;

    ispiv = NALLOC(char,znoc+1);
    memset(ispiv,0,(size_t)(znoc+1));
    locn = NALLOC(long,vdim-sdim+1);

    for (i = 1; i <= sdim; ++i)
    {
	if (ispiv[mypiv[i]])
	{
	    free(ispiv);
	    MTXFAIL(ERR_NOTECH,-1);
	}
	ispiv[mypiv[i]] = 1;
    }
    qdim = 0;
    for (i = 1; i <= znoc; ++i)
    {
	if (!ispiv[i])
	    locn[++qdim] = i;
    }
    free(ispiv);
    return 0;
}


/* ------------------------------------------------------------------
   zquot() - Projection of rows on the quotient.

   Parameters: <space> must be a matrix with <vdim> columns and
	<dim> rows. <quot> must be a matrix with <dim> rows and
	<qdim> columns.

   Return value: 0 = ok, -1=error
   ------------------------------------------------------------------ */

int zquot(space,dim,quot)
PTR space;
long dim;
PTR quot;

{
    long i, k;
    PTR qx, sx, tmp;

    qx = quot;
    zsetlen(vdim);
    sx = space;
    tmp = zalloc((long)1);

    for (i = 1; i <= dim; ++i)
    {
	zsetlen(vdim);
	zmoverow(tmp,sx);
	zadvance(&sx,(long)1);
	zcleanrow(tmp,mysubsp,sdim,mypiv);
	zsetlen(qdim);
	zmulrow(qx,F_ZERO);
	for (k = 1; k <= qdim; ++k)
	    zinsert(qx,k,zextract(tmp,locn[k]));
	zadvance(&qx,(long)1);
    }
    free(tmp);
    zsetlen(vdim);
    return 0;
}


/* ------------------------------------------------------------------
   zquotop() - Calculate the action of a square matrix on the
	quotient.

   Parameters: <matrix> must be a square matrix of dimension <vdim>.
	<quot> must be a square matrix of dimension <qdim>.

   Return value: 0=ok, -1=error.
   ------------------------------------------------------------------ */

int zquotop(matrix,quot)
PTR matrix;
PTR quot;

{
    long i, k;
    PTR qx, sx, tmp;

    qx = quot;
    zsetlen(vdim);
    tmp = zalloc((long)1);

    for (i = 1; i <= qdim; ++i)
    {
	zsetlen(vdim);
        sx = matrix;
	zadvance(&sx,locn[i]-1);
	zmoverow(tmp,sx);
	zcleanrow(tmp,mysubsp,sdim,mypiv);
	zsetlen(qdim);
	zmulrow(qx,F_ZERO);
	for (k = 1; k <= qdim; ++k)
	    zinsert(qx,k,zextract(tmp,locn[k]));
	zadvance(&qx,(long)1);
    }
    free(tmp);
    zsetlen(vdim);
    return 0;
}
