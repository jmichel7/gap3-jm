/* ========================== C MeatAxe =============================
   zfile.c - Stream based i/o for simple data types.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: zfile.c,v 1.2 1997/09/11 15:43:50 gap Exp $
 *
 * $Log: zfile.c,v $
 * Revision 1.2  1997/09/11 15:43:50  gap
 * New version 2.2.3. AH
 *
 * Revision 2.14  1995/02/08  10:14:52  mringe
 * _PL entfernt.
 *
 * Revision 2.14  1995/02/08  10:14:52  mringe
 * _PL entfernt.
 *
 * Revision 2.13  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.12  1994/07/28  06:04:43  mringe
 * zsetfield() und zsetlen() als getrennte Funktionen.
 *
 * Revision 2.11  1994/07/23  17:01:36  mringe
 * Profiling
 *
 * Revision 2.10  1994/05/20  12:41:25  mringe
 * Prototypen.
 *
 * Revision 2.9  1994/03/26  06:36:37  mringe
 * zseek() fuer 64-Bit-Rechner ver"andert.
 *
 * Revision 2.8  1994/03/13  13:48:19  mringe
 * Teste, ob Byteordnung stimmt.
 *
 * Revision 2.7  1994/02/15  11:23:50  mringe
 * SFSeek() durch SFSeek() ersetzt.
 *
 * Revision 2.6  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.5  1994/02/12  04:10:13  mringe
 * UMFANGREICHE AENDERUNGEN AN VIELEN DATENTYPEN.
 *
 * Revision 2.4  1994/01/18  14:32:05  mringe
 * *** empty log message ***
 *
 * Revision 2.3  1993/12/11  06:12:45  mringe
 * Tippfehler.
 *
 * Revision 2.2  1993/12/11  04:13:25  mringe
 * File I/O-Funktionen: kleine Korrekturen.
 *
 * Revision 2.1  1993/10/20  18:17:07  mringe
 * MeatAxe-2.0, Phase II.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.16  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.15  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.14  1993/10/05  23:35:02  mringe
 * zerrno eliminiert.
 *
 * Revision 1.13  1993/10/05  17:24:55  mringe
 * Kleine Schoenheitskorrekturen.
 *
 * Revision 1.12  1993/10/05  14:17:48  mringe
 * zreadvec()/zwritevec(): Bugs behoben.
 *
 * Revision 1.11  1993/10/05  11:54:02  mringe
 * Neu: zreadvec(), zwritevec()
 *
 * Revision 1.10  1993/08/14  11:50:26  mringe
 * Benutze HAS_UNSTD_H
 *
 * Revision 1.9  1993/08/10  14:51:42  mringe
 * Include string.h
 *
 * Revision 1.8  1993/08/06  14:40:42  mringe
 * Dokumentation entfernt, ist jetzt in zlib.tex.
 *
 * Revision 1.7  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.6  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.5  1993/07/28  13:34:49  mringe
 * *** empty log message ***
 *
 * Revision 1.4  1993/07/23  13:46:27  mringe
 * OS-Symbole neu (SYS_xxx)
 *
 * Revision 1.3  1993/07/13  20:30:59  mringe
 * Neue File i/o library.
 *
 * Revision 1.2  1993/05/25  13:08:07  mringe
 * Deklarationen von BPR und LPR in zzz.h
 *
 * Revision 1.1  1993/05/25  12:49:01  mringe
 * Initial revision
 *
 */

#include <string.h>
#include <stdlib.h>
#include "meataxe.h"


#define ZZZFATAL(msg) fatal(msg,__LINE__)


/* ------------------------------------------------------------------
   Prototypes
   ------------------------------------------------------------------ */

static void fatal(char *msg, int line);
static void swapbytes(long *l);


/* -----------------------------------------------------------------
   zreadlong() - Read n long integers from a file.
   zwritelong() - Write n long ints to a file.
   ------------------------------------------------------------------ */

static long test = 1;

size_t zreadlong(f,buf,n)
FILE *f;
long *buf;
size_t n;

{
    unsigned char a[4];
    unsigned nread;
    PROFILE_BEGIN(t);

    if (sizeof(long) == 4 && *(char *)&test == 1)
	nread = fread((char *)buf,sizeof(long),n,f);
    else
    {
    	for (nread = 0; nread < n; ++nread)
    	{
	    if (fread(a,1,4,f) != 4) break;
	    *buf = ((unsigned long)a[0]|((unsigned long)a[1] << 8)|
		((unsigned long)a[2] << 16)|((long)(char)a[3] << 24));
	    ++buf;
    	}
    }
    PROFILE_END(t,FileIO);
    return nread;
}

size_t zwritelong(f,buf,n)
FILE *f;
long *buf;
size_t n;

{
    unsigned char a[4];
    unsigned nwritten;
    PROFILE_BEGIN(t);

    if (sizeof(long) == 4 && *(char *)&test == 1)
	nwritten = fwrite((char *)buf,sizeof(long),n,f);
    else
    {
    	for (nwritten = 0; nwritten < n; ++nwritten)
    	{
    	    a[0] = (unsigned char) *buf;
    	    a[1] = (unsigned char) (*buf >> 8);
    	    a[2] = (unsigned char) (*buf >> 16);
    	    a[3] = (unsigned char) (*buf >> 24);
    	    if (fwrite(a,1,4,f) != 4) break;
	    ++buf;
    	}
    }
    PROFILE_END(t,FileIO);
    return nwritten;
}


/* ------------------------------------------------------------------
   zreadvec() - Read packed vectors from a file
   zwritevec() - Write packed vectors to a file
   ------------------------------------------------------------------ */

size_t zreadvec(f,buf,n)
FILE *f;
PTR  buf;
size_t n;

{
    size_t i;
    register char *b = (char *) buf;
    PROFILE_BEGIN(t);

    for (i = 0; i < n; ++i)
    {
        if (fread(b,zrowsize_io,1,f) != 1) break;
	b += zrowsize;
    }
    if (ferror(f)) mtxerrno = ERR_FILEREAD;
    PROFILE_END(t,FileIO);
    return i;
}


size_t zwritevec(f,buf,n)
FILE *f;
PTR  buf;
size_t n;

{
    size_t i;
    register char *b = (char *) buf;
    PROFILE_BEGIN(t);

    for (i = 0; i < n; ++i)
    {
        if (fwrite(b,zrowsize_io,1,f) != 1) break;
	b += zrowsize;
    }
    if (ferror(f)) mtxerrno = ERR_FILEWRITE;
    PROFILE_END(t,FileIO);
    return i;
}


/* ------------------------------------------------------------------
   zseek()  -  Set read/write position
   ------------------------------------------------------------------ */

int zseek(f, pos)
FILE *f;
long pos;

{
    long addr;

    if (zfl != -1)
	addr = zrowsize_io * (pos - 1) + 12;
    else
	addr = znoc * 4 * (pos - 1) + 12;
    if (SFSeek(f,addr) == -1)
	MTXFAIL(ERR_FILEREAD,-1);
    return 0;
}


/* ------------------------------------------------------------------
   zreadhdr() - Open file and read file header
   zwritehdr() - Open file and write file header

   Return value: File handle or NULL on error
   ------------------------------------------------------------------ */

FILE *zreadhdr(name, field, nrows, ncols)
char *name;
long *field, *nrows, *ncols;

{
    FILE *fd;
    long header[3];

    fd = SFOpen(name,FM_READ);
    if (fd == NULL)
	MTXFAIL(ERR_FILEOPEN,NULL);
    if (zreadlong(fd,header,3) != 3)
    {
	fclose(fd);
	MTXFAIL(ERR_FILEREAD,NULL);
    }

if (header[0] > 2)
{
long flen;
SFSeek(fd,(long)-1);	/* Seek to eof */
flen = ftell(fd)-12;
SFSeek(fd,(long)12);

if (
	flen != header[1] * header[2] &&
	flen != header[1] * ((header[2]-1)/2+1) &&
	flen != header[1] * ((header[2]-1)/3+1) &&
	flen != header[1] * ((header[2]-1)/4+1) &&
	flen != header[1] * ((header[2]-1)/5+1) &&
	flen != header[1] * ((header[2]-1)/8+1) 
    )
    {fprintf(stderr,"%s: Bad file size, run conv\n",name);
    ZZZFATAL("ERROR READING FILE HEADER");
    }
}

	/* Try to correct wrong byte ordering.
	   ----------------------------------- */
	if (header[0] > 256 || header[1] < 0 || header[2] < 0)
	{	swapbytes(header);
		swapbytes(header+1);
		swapbytes(header+2);
		if (header[0] > 256 || header[1] < 0 || header[2] < 0)
		{	fprintf(stderr,"%s: Invalid header\n",name);
			ZZZFATAL("ERROR IN FILE HEADER");
		}
	}

	*field = header[0];
	*nrows = header[1];
	*ncols = header[2];

    return fd;
}


FILE *zwritehdr(name, field, nrows, ncols)
char *name;
long field;
long nrows, ncols;

{
    FILE *fd;
    long header[3];

    header[0] = (long) field;
    header[1] = (long) nrows;
    header[2]= (long) ncols;
    fd = SFOpen(name,FM_CREATE);
    if (fd == NULL) MTXFAIL(ERR_FILEOPEN,NULL);
    if (zwritelong(fd,header,3) != 3)
    {
	fclose(fd);
	MTXFAIL(ERR_FILEWRITE,NULL);
    }
    return fd;
}


/* ------------------------------------------------------------------
   swapbytes() - Swap bytes in a 4-byte block: ABCD -> DCBA
   ------------------------------------------------------------------ */

static void swapbytes(l)
long *l;

{	char inp[4], out[4];

	memcpy(inp,l,4);
	out[0] = inp[3];
	out[1] = inp[2];
	out[2] = inp[1];
	out[3] = inp[0];
	memcpy(l,out,4);
}

/* ------------------------------------------------------------------
   fatal() - Fatal error: Print message and exit
   ------------------------------------------------------------------ */

static void fatal(msg, line)
char *msg;
int line;

{	fprintf(stderr,"ZZZ ERROR (%s, line %d) - %s\n",
		__FILE__,line,msg);
	exit(EXIT_ERR);
}


