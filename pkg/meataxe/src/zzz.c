/* ========================== C MeatAxe =============================
   Finite field arithmetic and common functions. `Small' version for
   field orders q <= 256. Originally based on the `hprout.c' written
   by Klaus Lux.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zzz.c,v 1.2 1997/09/11 15:44:42 gap Exp $
 *
 * $Log: zzz.c,v $
 * Revision 1.2  1997/09/11 15:44:42  gap
 * New version 2.2.3. AH
 *
 * Revision 2.20  1995/06/22  13:19:45  mringe
 * Teste, ob MeatAxeBinDir != 0
 *
 * Revision 2.19  1995/05/12  10:07:20  mringe
 * Benutze MeatAxeBinDir beim Aufruf von maketab.
 *
 * Revision 2.18  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.17  1994/11/25  14:08:17  mringe
 * ANSI-C, neue Namen fuer time-Funktionen.
 *
 * Revision 2.16  1994/09/20  20:43:08  mringe
 * *** empty log message ***
 *
 * Revision 2.15  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.15  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.14  1994/07/23  16:50:01  mringe
 * Prototypen.
 *
 * Revision 2.13  1994/07/22  18:55:37  mringe
 * BUG in maprow behoben.
 *
 * Revision 2.12  1994/07/21  11:17:35  mringe
 * zzz_cc.
 *
 * Revision 2.11  1994/06/25  13:36:01  mringe
 * Neu: zextractcol().
 *
 * Revision 2.10  1994/06/14  09:58:21  mringe
 * Teste ob subfield==zfl.
 *
 * Revision 2.9  1994/05/19  11:34:50  mringe
 * Benutze Smalloc().
 *
 * Revision 2.8  1994/05/18  12:22:21  mringe
 * Vermeide malloc(0).
 *
 * Revision 2.7  1994/04/09  13:27:43  mringe
 * memset() durch Schleifen ersetzt.
 *
 * Revision 2.6  1994/03/30  08:52:13  mringe
 * zextract() ist jetzt ein makro.
 *
 * Revision 2.5  1994/02/15  14:09:19  mringe
 * zsetlen() erweitert fuer Permutationen.
 *
 * Revision 2.4  1994/02/15  10:28:33  mringe
 * MSDOS_BCC entfernt.
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
 */


#include <string.h>
#include <stdlib.h>
#include "meataxe.h"


/* ------------------------------------------------------------------
   Gobal data
   ------------------------------------------------------------------ */

char zzz_cc[] = ZZZ_CC;
char zzzversion[] = "conway-qadic-1.0";
typedef unsigned char BYTE;
BYTE
	tmult[256][256],
	tadd[256][256],
	taddinv[256],
	tmultinv[256];

static BYTE tffirst[256][2],
	textract[8][256],
	tnull[8][256],
	tinsert[8][256];
static long embedord[4];	/* Subfield orders */
static BYTE embed[4][16];	/* Embeddings of subfields */
static BYTE restrict[4][256];	/* Restrictions to subfields */

long zchar = 0;		/* Characteristic */
long zfl = -1;		/* Field order */
FEL zgen = -1;		/* Generator */
long znoc = 0;		/* No. of columns for row operations */
size_t zrowsize = 0;	/* Row size in memory */
size_t zrowsize_io = 0;	/* Row size for file i/o */

static size_t MPB = 0;	/* No. of marks per byte */
static size_t LPR = 0;	/* Long ints per row */


/* ------------------------------------------------------------------
   Argument checking (for debbugging)
   ------------------------------------------------------------------ */

#if defined(DEBUG)
#define CHECKRANGE(x,lo,hi) if ((x)<(lo)||(x)>(hi)) {\
	fprintf(stderr,"%ld <= %ld <= %ld ?\n",(long)(lo),\
		(long)(x),(long)(hi));\
	FATAL("ZZZ RANGE CHECK ERROR");}
#else
#define CHECKRANGE(x,lo,hi)
#endif

#define CHECKFEL(x) CHECKRANGE(x,0,zfl-1)
#define CHECKCOL(x)  CHECKRANGE(x,1,znoc)



/* ------------------------------------------------------------------
   zinsert()  -  insert mark into row

   return value: none
   ------------------------------------------------------------------ */

void zinsert(PTR row, long col, FEL mark)

{
    register BYTE *loc = (BYTE *)row + ((col - 1) / MPB);
    register long idx = (col - 1) % MPB;

    CHECKFEL(mark);
    *loc = tnull[idx][*loc] + tinsert[idx][mark];
}



/* ------------------------------------------------------------------
   zextract()  -  extract mark from row

   return value: mark extracted from row[col]
   ------------------------------------------------------------------ */

FEL zextract(PTR row, long col)

{
    register BYTE *loc = (BYTE *)row + ((col - 1) / MPB);
    return (textract[(col - 1) % MPB][*loc]);
}


/* ------------------------------------------------------------------
   zfindpiv()  -  Find first non-zero mark in row

   return value: column (0 if row is zero)
   ------------------------------------------------------------------ */

long zfindpiv(PTR row, FEL *mark)

{
    register long *l = (long *) row;
    register long idx;
    register BYTE *m;

    for (idx = 0; *l == 0 && idx < LPR; ++idx, ++l);
    if (*l == 0) return 0;
    idx = idx * sizeof(long) * MPB + 1;
    m = (BYTE *)l;
    while (*m == 0)
    {	++m;
	idx += MPB;
    }
    idx += tffirst[*m][1];
    if (idx <= znoc)
    {	*mark = tffirst[*m][0];
	return idx;
    }
    return 0;
}


/* -----------------------------------------------------------------
   zsetfiel()  -  Select the field. Read table file if necessary.
   ------------------------------------------------------------------ */

int zsetfield(long field)

{   static char filename[50];
    FILE *fd = NULL;
    long hdr[5];
    int count;

    if (field == zfl || field < 2) return 0;
    
    /* Open table file
       --------------- */
    for (count = 2; count >= 0; --count)
    {   sprintf(filename,"p%3.3ld.zzz",(long)field);
	if ((fd = SFOpen(filename,FM_READ|FM_LIB)) != NULL)
	    break;	/* O.k. */
	switch (count)
	{
	    case 0:
	    	perror(filename);
	        FATAL("ERROR OPENING TABLE FILE");
		break;
	    case 1:
		sprintf(filename,"maketab -Q %ld",field);
		system(filename);
		break;
	    case 2:
		if (*MeatAxeBinDir != 0)
		    sprintf(filename,"%s/maketab -Q %ld",
			MeatAxeBinDir,field);
		else
		    sprintf(filename,"maketab -Q %ld",field);
		system(filename);
		break;
	}
    }

    /* Read header, perform some checks
       -------------------------------- */
    if (zreadlong(fd,hdr,5) != 5)
    {
	perror(filename);
	FATAL("ERROR READING FILE HEADER");
    };
    if (hdr[2] != field || hdr[1] < 0 || hdr[1] > field ||
        hdr[0] <= 1 || hdr[2] % hdr[0] != 0 || hdr[3] < 1 || hdr[3] > 8)
    {	fprintf(stderr,"%s: Invalid header\n", filename);
	FATAL("TABLE FILE CORRUPTED");
    }
    zchar = hdr[0];
    zgen = (FEL) hdr[1];
    MPB = hdr[3];
    if (hdr[4] != (long) ZZZVERSION)
	fprintf(stderr,"BAD VERSION NUMBER IN TABLE FILE");

    /* Read tables
       ----------- */
    if (fread(tmult,4,sizeof(tmult)/4,fd) != sizeof(tmult)/4 ||
        fread(tadd,4,sizeof(tadd)/4,fd) != sizeof(tadd)/4 ||
        fread(tffirst,1,sizeof(tffirst),fd) != sizeof(tffirst) ||
        fread(textract,1,sizeof(textract),fd) != sizeof(textract) ||
        fread(taddinv,1,sizeof(taddinv),fd) != sizeof(taddinv) ||
        fread(tmultinv,1,sizeof(tmultinv),fd) != sizeof(tmultinv) ||
        fread(tnull,1,sizeof(tnull),fd) != sizeof(tnull) ||
        fread(tinsert,1,sizeof(tinsert),fd) != sizeof(tinsert) ||
        zreadlong(fd,embedord,4) != 4 ||
        fread(embed,16,4,fd) != 4 ||
        fread(restrict,256,4,fd) != 4
       )
    {
	perror(filename);
        FATAL("ERROR READING FILE");
    }
    fclose(fd);
    zfl = field;
    zsetlen(zfl);
    return 0;
}


/* -----------------------------------------------------------------
   zsetlen()  -  Set row length.
   ------------------------------------------------------------------ */

int zsetlen(long ncols)

{ 
    znoc = ncols;
    zrowsize_io = (znoc-1) / MPB + 1;		/* Used bytes per row */
    LPR = (zrowsize_io-1) / (sizeof(long)) + 1;	/* Ints per row */
    zrowsize = LPR * sizeof(long);	/* Total bytes per row */
    return 0;
}



/* ------------------------------------------------------------------
   zembed() - Embed a subfield element.

   return value: Embedded field element.
   ------------------------------------------------------------------ */

FEL zembed(FEL a, long subfield)

{   int i;

    if (subfield == zfl) return a;
    for (i = 0; embedord[i] != subfield && i < 4; ++i);
    if (i >= 4)
	FATAL("NO SUCH SUBFIELD");
    return embed[i][a];
}




/* ------------------------------------------------------------------
   zrestrict() - Restrict to a subfield.

   return value: Embedded field element.
   ------------------------------------------------------------------ */

FEL zrestrict(FEL a, long subfield)

{
    int i;

    if (subfield == zfl) return a;
    for (i = 0; embedord[i] != subfield && i < 4; ++i);
    if (i >= 4)
	FATAL("NO SUCH SUBFIELD");
    return restrict[i][a];
}



/* ------------------------------------------------------------------
   zaddrow()  -  Add two rows (row1 += row2)
   ------------------------------------------------------------------ */

PTR zaddrow(PTR dest, PTR src)

{
    register long i;

    if (zchar == 2)	/* characteristic 2 is simple... */
    {	register long *l1 = (long *) dest;
	register long *l2 = (long *) src;
	for (i = LPR; i != 0; --i)
	{
	    register long x = *l2++;
	    if (x != 0) *l1 ^= x;
	    l1++;
	}
    }
    else		/* any other characteristic */
    {	register BYTE *p1 = dest;
	register BYTE *p2 = src;
	for (i = zrowsize_io; i != 0; --i)
	{
	    register int x = *p2++;
	    if (x != 0) *p1 = tadd[*p1][x];
	    p1++;
	}
    }
    return dest;
}



/* ------------------------------------------------------------------
   zmulrow() - Multiply row by field element
   ------------------------------------------------------------------ */

PTR zmulrow(PTR row, FEL mark)

{
    register BYTE *m;
    register BYTE *multab;
    register long i;

    CHECKFEL(mark);
    if (mark == F_ZERO)
    {
        register long *l = (long *)row;
        for (i = LPR; i > 0; --i) *l++ = 0;
    }
    else if (mark != F_ONE)
    {
	multab = tmult[mark];
	m = row;
	for (i = zrowsize_io; i != 0; --i)
	{
	    register int x = *m;
	    if (x != 0) *m = multab[x];
	    ++m;
	}
    }
    return row;
}



/* ------------------------------------------------------------------
   zaddmulrow() - row1 += f * row2

   return value: none
   ------------------------------------------------------------------ */

PTR zaddmulrow(PTR dest, PTR src, FEL f)

{
    register long i;
    register BYTE *p1, *p2, *multab;

    CHECKFEL(f);
    if (f == F_ZERO) return dest;
    if (f == F_ONE) return zaddrow(dest,src);
    multab = tmult[f];
    p1 = dest;
    p2 = src;
    for (i = zrowsize_io; i != 0; --i)
    {
	*p1 = tadd[*p1][multab[*p2++]];
	++p1;
    }
    return dest;
}


/* -----------------------------------------------------------------
   zmaprow()  -  multiply a row by a matrix
   ------------------------------------------------------------------ */

PTR zmaprow(PTR row, PTR matrix, long nor, PTR result)

{
    register int i;
    register FEL f;
    BYTE *m = (BYTE *) matrix;
    register long *l = (long *)result;

    for (i = LPR; i > 0; --i) *l++ = 0;

    if (zfl == 2)	/* GF(2) is a special case */
    {
	register long *x1 = (long *) matrix;
	register BYTE *r = (BYTE *) row;
	register BYTE mask;

	for (i = nor; i > 0; ++r)
    	{
	    if (*r == 0)	/* Skip zero bytes */
	    {
		x1 += 8 * LPR;
		i -= 8;
		continue;
	    }

	    for (mask = 0x80; mask != 0 && i > 0; --i, mask >>= 1)
	    {
	        if ((*r & mask) != 0)
	        {
		    register long *x2 = (long *)result;
		    register long k;
		    for (k = LPR; k; --k)
			*x2++ ^= *x1++;
	    	}
	    	else
		    x1 += LPR;	/* Skip that row */
	    }
	}
    }
    else		/* Any other field */
    {
	register BYTE *brow = (BYTE *) row;
	register int pos = 0;

	for (i = nor; i > 0; --i)
	{
	    f = textract[pos][*brow];
	    if (++pos == MPB)
	    {
		pos = 0;
		++brow;
	    }
	    if (f != F_ZERO)
	    {
		register BYTE *v = m;
		register BYTE *r = result;
		register long k = zrowsize_io;
		if (f == F_ONE)
		{
		    for (; k != 0; --k)
		    {
			*r = tadd[*r][*v++];
			++r;
		    }
		}
		else
		{
		    register BYTE *multab = tmult[f];
		    for (; k != 0; --k)
		    {
			*r = tadd[*r][multab[*v++]];
			++r;
		    }
		}
	    }
	    m += zrowsize;		/* next row */
	}
    }
    return result;
}


/* ------------------------------------------------------------------
   zitof(), zftoi() - Convert between long and FEL
   ------------------------------------------------------------------ */

FEL zitof(long l)

{
    register long f;
    f = l % zfl;
    if (f < 0) f += zfl;
    return (FEL) f;
}

long zftoi(FEL f)

{	return (long) f;
}



/* ------------------------------------------------------------------
   zextractcol() - Extract a column. The size of result must be at
   least <nor>. <mat>=<result> is NOT allowed.
   ------------------------------------------------------------------ */

PTR zextractcol(PTR mat,long nor,long col,PTR result)

{
    register BYTE *x = (BYTE *)mat + ((col - 1) / MPB);
    register BYTE *extab = textract[(col - 1) % MPB];
    register BYTE a = 0;
    register int ind = 0;
    register BYTE *y = result;
    register long count;
    
    for (count = nor; count > 0; --count)
    {
	a += tinsert[ind][extab[*x]];
	if (++ind == MPB)
	{
	    *y++ = a;
	    a = 0;
	    ind = 0;
	}
	x += zrowsize;
    }
    if (ind != 0) *y = a;
    return result;
}

