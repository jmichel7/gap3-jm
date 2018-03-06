/* ========================== C MeatAxe =============================
   conv.c - Convert machine dependent binary files.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: conv.c,v 1.2 1997/09/11 15:42:44 gap Exp $
 *
 * $Log: conv.c,v $
 * Revision 1.2  1997/09/11 15:42:44  gap
 * New version 2.2.3. AH
 *
 * Revision 2.9  1995/02/11  21:58:17  mringe
 * ANSI C
 *
 * Revision 2.8  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.7  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.6  1994/07/07  06:57:14  mringe
 * Include files.
 *
 * Revision 2.5  1994/03/30  08:52:13  mringe
 * Option -o.
 *
 * Revision 2.4  1994/03/26  06:33:05  mringe
 * Wieder auf die alte Version zur"uckgegriffen,
 * 2.3 lief nicht auf 64-Bit-Maschinen.
 *
 * Revision 2.2  1994/02/15  13:39:15  mringe
 * Benutze SFSeek().
 *
 * Revision 2.1  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.17  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.16  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.15  1993/08/14  11:50:26  mringe
 * Test des Fileheaders verbessert.
 *
 * Revision 1.14  1993/08/10  14:29:19  mringe
 * Include string.h
 *
 * Revision 1.13  1993/08/10  06:51:38  mringe
 * Hole SEEK_xxx aus unistd.h, falls nicht definiert.
 *
 * Revision 1.12  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.11  1993/07/28  13:37:40  mringe
 * *** empty log message ***
 *
 * Revision 1.10  1993/07/23  13:46:27  mringe
 * OS-Symbole neu (SYS_xxx)
 *
 * Revision 1.9  1993/06/09  18:36:53  mringe
 * singed/unsigned Bug behoben.
 *
 * Revision 1.8  1993/06/09  13:37:27  mringe
 * Bug beim Erkennen der richtigen Bytefolge behoben.
 *
 * Revision 1.7  1993/05/18  12:17:50  mringe
 * Neues Fileformat. Benutze zreadlong(). Konvertiere nur fuer die
 *  Maschine, auf dem das Programm laeuft.
 *
 * Revision 1.6  1993/02/16  18:32:46  mringe
 * string.h und stdio.h werden jetzt in meataxe.h included.
 *
 * Revision 1.5  1993/02/16  17:33:20  mringe
 * Symbole SYS_BSD und SYS_SYSV eingefuehrt.
 *
 * Revision 1.4  1992/11/19  09:47:44  mringe
 * Diverse kleine "Anderungen; Neu: Option -c
 *
 * Revision 1.3  1992/11/19  09:34:26  mringe
 * Kleine "Anderungen, help().
 *
 */


#include <stdlib.h>
#include <string.h>
#include "meataxe.h"



#define BUFSIZE 4096


/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

static char filename[100];
static char tmpname[100] = "conv.tmp";
static FILE *f, *tmp;
static long flen;
static long hdr[3];
static int opt_o = 0;

static char *helptext[] = {
"SYNTAX",
"    conv [<Options>] <File> ...",
"",
"OPTIONS",
"    -o    Convert back to old format.",
"FILES",
"    <File> is any MeatAxe binary data file (matrix or permutation).",
"    The converted file replaces the old one.",
NULL};


static proginfo_t pinfo =
   { "conv", "Convert Binary Files", "$Revision: 1.2 $", helptext };





/* ------------------------------------------------------------------
   mywritelong() - Write long integers
  
   ------------------------------------------------------------------ */

static size_t mywritelong(FILE *f,long *buf,size_t n)

{
    if (opt_o)
    	return fwrite(buf,sizeof(long),n,f);
    else
    	return zwritelong(f,buf,n);
}




/* ------------------------------------------------------------------
   swapl() - Swap bytes in a long integer
   ------------------------------------------------------------------ */

static long swapl(long l)

{
    unsigned char a,b,c,d;

    a = (unsigned char) (l & 0xFF);
    b = (unsigned char) ((l>>8) & 0xFF);
    c = (unsigned char) ((l>>16) & 0xFF);
    d = (unsigned char) ((l>>24) & 0xFF);
    return ((unsigned long) d |
            ((unsigned long) c << 8) |
            ((unsigned long) b << 16) |
            ((unsigned long) a << 24));
}


/* ------------------------------------------------------------------
   examine() - Decide if swapping is necessary.
   ------------------------------------------------------------------ */

static int examine(char *buf,long wordsize,long n)

{
    unsigned char *c = (unsigned char *) buf;
    int i, result1, result2;
    long l, count;

    for (count = n; count > 0; --count, c += wordsize)
    {
	for (l = 0, i = 0; i < wordsize; ++i)
	    l = (l << 8) | c[i];
	result1 = (l > 0 && l <= n);
	for (l = 0, i = wordsize-1; i >= 0; --i)
	    l = (l << 8) | c[i];
	result2 = (l > 0 && l <= n);

	if (result1 && !result2) return 1;	/* Swap */
	if (result2 && !result1) return 0;	/* Don't swap */
    }
    return -1;	/* Failed */
}



/* ------------------------------------------------------------------
   intcopy() - Copy integers with different word sizes and
	byte orders.
   ------------------------------------------------------------------ */

static void intcopy(long *dest, char *src, long wordsize, int swap, 
	long n)

{
    unsigned char *c;
    int i;
    long count;
    long *d, l;

    d = dest;
    c = (unsigned char *) src;
    for (count = 0; count < n; ++count)
    {
	l = 0;
	if (swap)
            for (i = 0; i < wordsize; ++i)
		l = (l << 8) | c[i];
	else
            for (i = wordsize; i >= 0; --i)
		l = (l << 8) | c[i];
        *d++ = l;
	c += wordsize;
    }
}


/* ------------------------------------------------------------------
   chkhdr() - Check file header. Returns 1 if hdr is a valid matrix
	header, 2 for a valid permutation header.
   ------------------------------------------------------------------ */

static int chkhdr(long *hdr)

{	if (hdr[0] > 1 && hdr[0] <= 256 && hdr[1] >= 0 && hdr[2] >= 0)
		return 1;	/* Matrix */
	if (hdr[0] == -1 && hdr[1] > 0 && hdr[2] > 0 &&
	    hdr[1] < 0x20000000 && hdr[2] < 0x20000000)
		return 2;		/* Permutation */
	return 0;
}


/* ------------------------------------------------------------------
   convperm() - Convert a permutation
   ------------------------------------------------------------------ */

static int convperm(void)


{
    char *rbuf;
    long *wbuf;
    long i;
    long wordsize;
    int swap = -1;
    static long test = 1;

    if (hdr[1]*hdr[2] == 0)
	wordsize = sizeof(long);	/* Egal, was hier steht! */
    else
	wordsize = (flen - 12) / (hdr[1]*hdr[2]);
    if (wordsize*hdr[1]*hdr[2] != flen - 12)
    {
	printf("Bad file size");
	return 1;
    }
    rbuf = (char *) malloc((size_t)hdr[1]*wordsize);
    wbuf = (long *) malloc((size_t)hdr[1]*sizeof(long));

    for (i = 0; i < hdr[2]; ++i)
    {
	if (fread(rbuf,(size_t)wordsize,(size_t)hdr[1],f) != hdr[1])
	{
	    printf("Error reading file");
	    return 1;
	}

	if (swap == -1)
	{
	    swap = examine(rbuf,wordsize,hdr[1]);
	    if (swap == -1)
	    {
	    	printf("Cannot decide if swapping is necessary");
	    	return 1;
	    }
	}
	intcopy(wbuf,rbuf,wordsize,swap,hdr[1]);
	if (mywritelong(tmp,wbuf,(size_t)hdr[1]) != hdr[1])
	{
	    printf("Error writing temporary file");
	    return 1;
	}
    }
    printf("%ld permutation(s) on %ld points",hdr[2],hdr[1]);
    if (wordsize != 4)
	printf(" (word size changed from %d to 4)",(int)wordsize);

    if (swap ^ (*((char*)&test) != 1))
	printf(" (swapped)");
    return 0;
}


/* ------------------------------------------------------------------
   convmatrix() - Convert a matrix
   ------------------------------------------------------------------ */

static int convmatrix()

{
    long bpr, newbpr, i;
    char *buf;
	
    if (hdr[1] == 0) bpr = 0;
    else bpr = (flen - 12) / hdr[1];
    if ((flen - 12) % hdr[1])
    {
	printf("Bad file size");
	return 1;
    }
    for (newbpr = 0, i = hdr[0]; i <= 256; ++newbpr, i*=hdr[0]);
    newbpr = (hdr[2] - 1) / newbpr + 1;
    if (newbpr > bpr || bpr - newbpr > 8)
    {
	printf("Bad file size");
	return 1;
    }
    buf = (char *) malloc((size_t)bpr);
    for (i = 0; i < hdr[1]; ++i)
    {
	if (fread(buf,(size_t)bpr,1,f) != 1)
	{
	    printf("Error reading file");
	    return 1;
	}
	if (fwrite(buf,(size_t)newbpr,1,tmp) != 1)
	{
	    printf("Error writing temporary file");
	    return 1;
	}
    }
    printf("%ldx%ld matrix over GF(%ld)",hdr[1],hdr[2],hdr[0]);
    return 0;
}


/* ------------------------------------------------------------------
   convmert() - Convert a matrix or permutation
   ------------------------------------------------------------------ */

static int convert(void)

{
    /* Lese File-Header, pruefe auf richtige Byte-Ordnung
       und schreibe den Header in die temoraere Datei.
       -------------------------------------------------- */
    SFSeek(f,(long)0);
    if (zreadlong(f,hdr,3) != 3)
    {
	printf("Error reading file header");
	return 1;
    }
    if (chkhdr(hdr) == 0)
    {
	hdr[0] = swapl(hdr[0]);
	hdr[1] = swapl(hdr[1]);
	hdr[2] = swapl(hdr[2]);
	if (chkhdr(hdr) == 0)
	{
	    printf("Not a MeatAxe file");
	    return 1;
	}
    }
    if (mywritelong(tmp,hdr,3) != 3)
    {
	printf("Error writing file header");
	return 1;
    }

    /* Konvertiere den Rest des Files
       ------------------------------ */
    if (hdr[0] >= 2) /* Matrix */
	return convmatrix();
    else	/* Permutation */
	return convperm();
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
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("o")) != OPT_END)
    {
    	switch (i)
    	{
    	    case 'o': opt_o = 1; break;
    	}
    }

    for (i = opt_ind; i < argc; ++i)
    {
	/* Open file, get file size
	   ------------------------ */
	strcpy(filename,argv[i]);
	f = SFOpen(filename,FM_READ);
	if (f == NULL)
	{	perror(filename);
		continue;
	}
	printf("%s: ",filename);
	SFSeek(f,(long)-1);	/* End of file */
	flen = ftell(f);

	/* Open temporary file
	   ------------------- */
	tmp = SFOpen(tmpname,FM_CREATE);
	if (tmp == NULL)
	{	perror(tmpname);
		fclose(f);
		continue;
	}

	/* Convert the file
	   ---------------- */
	if (convert() == 0)	/* Success */
	{
	    fclose(f);
	    fclose(tmp);
	    remove(filename);
	    rename(tmpname,filename);
	}
	else			/* Discard tmp file */
	{
	    fclose(f);
	    fclose(tmp);
	    remove(tmpname);
	    printf(" - Not converted.");
	}
	printf("\n");
    }
    return EXIT_OK;
}

