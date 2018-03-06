/* ========================== C MeatAxe =============================
   files.c - Functions for reading and writing the .cfinfo file.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: files.c,v 1.2 1997/09/11 15:42:47 gap Exp $
 *
 * $Log: files.c,v $
 * Revision 1.2  1997/09/11 15:42:47  gap
 * New version 2.2.3. AH
 *
 * Revision 2.11  1995/02/09  14:04:19  mringe
 * ANSI C
 *
 * Revision 2.10  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 2.9  1994/08/29  17:03:00  mringe
 * Peakwoerter: Benutze Paar aus WortNr,Polynom (wie bei den IdWords in chop).
 *
 * Revision 2.8  1994/05/04  11:40:31  mringe
 * Polynome.
 *
 * Revision 2.7  1994/03/26  06:34:04  mringe
 * basename umbenannt wg. Namenskonflikt.
 *
 * Revision 2.6  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.5  1994/01/24  08:44:36  mringe
 * Benutze FM_READ|FM_TEXT und FM_CREATE|FM_TEXT.
 *
 * Revision 2.4  1993/10/27  09:34:00  mringe
 * Benutze 0 als Wert f"ur 'nicht bekannt'.
 *
 * Revision 2.3  1993/10/27  09:23:36  mringe
 * Neues cfinfo-Format.
 *
 * Revision 2.2  1993/10/22  17:00:14  mringe
 * *** empty log message ***
 *
 * Revision 2.1  1993/10/22  16:08:19  mringe
 * Neues Numerierungsschema fuer irreduzible.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.10  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.9  1993/08/10  14:51:42  mringe
 * Include string.h
 *
 * Revision 1.8  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.7  1993/07/23  13:46:27  mringe
 * OS-Symbole neu (SYS_xxx)
 *
 * Revision 1.6  1993/02/16  18:32:46  mringe
 * string.h und stdio.h werden jetzt in meataxe.h included.
 *
 * Revision 1.5  1993/02/16  17:33:20  mringe
 * Symbole SYS_BSD und SYS_SYSV eingefuehrt.
 *
 * Revision 1.4  1993/01/15  14:49:07  mringe
 * Speichere Anzahl der Erzeuger im cfinfo file
 *
 * Revision 1.3  1992/07/22  07:10:30  mringe
 * Changed 'global.h' to 'lattice.h'
 *
 * Revision 1.2  1992/07/15  09:25:55  mringe
 * Some minor changes.
 *
 * Revision 1.1  1992/05/26  07:08:08  mringe
 * Initial revision
 *
 */


#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "meataxe.h"
#include "lattice.h"
#include "files.h"


char cfbasename[MAXBASENAME];
cfinfo_t cfinfo[MAXCF];
int ncf=0;		/* Number of composition factors */
int ngen=2;		/* Number of generators */

static char fnbuf[MAXBASENAME+10];	/* File name buffer */





static FILE *myfopen(name, mode)
char *name;
int mode;

{
    FILE *f;

    f = SFOpen(name,mode);
    if (f == NULL)
    {
	perror(name);
	FATAL("CANNOT OPEN FILE");
    }
    return f;
}



void setbasename(bn)
char *bn;

{
    int i;

    if (strlen(bn) >= MAXBASENAME)
	FATAL("NAME TOO LONG");
    strcpy(cfbasename,bn);
    for (i = 0; i < MAXCF; ++i)
    {
	cfinfo[i].dim = 0;
	cfinfo[i].num = -1;
	cfinfo[i].mult = -1;
	cfinfo[i].idword = 0;
	cfinfo[i].idpol = NULL;
	cfinfo[i].peakword = 0;
	cfinfo[i].peakpol = NULL;
	cfinfo[i].spl = 0;
	cfinfo[i].ndotl = -1;
	cfinfo[i].nmount = -1;
    }
}


static char *fn_info()

{
    sprintf(fnbuf,"%s.cfinfo",cfbasename);
    return fnbuf;
}


FILE *f;
char buf[2000];
char *bp;

void nextline()

{
    if (feof(f))
	errexit(ERR_FILEFMT,fn_info());
    fgets(buf,sizeof(buf),f);
    buf[sizeof(buf)-1] = 0;
    if (ferror(f))
	errexit(ERR_FILEREAD,fn_info());
    bp = buf;
}


static int match(c)
char *c;

{
    char *b = bp;
    for (b = bp; *b != 0 && *c != 0; ++c)
    {
	if (*c == ' ')
	    while (*b != 0 && isspace(*b)) ++b;
	else
	    if (*b == *c) ++b;
	    else return 0;
    }
    if (*c == 0)
    {	
	bp = b;
	return 1;	/* Success */
    }
    return 0;		/* Failed */
}


static long readint()

{
    long l;
    int neg = 0;

    while (isspace(*bp)) ++bp;
    if (*bp == '-') { neg = 1; ++bp; }
    if (!isdigit(*bp)) errexit(ERR_FILEFMT,fn_info());
    l = atol(bp);
    while (isdigit(*bp)) ++ bp;
    return neg ? -l : l;
}

static long scanint()

{
    long l;

    l = readint();
    if (*bp == ',')
	++bp;
    else if (*bp == ']')
	if (bp[1] != ';') errexit(ERR_FILEFMT,fn_info());
    return l;
}


/* ------------------------------------------------------------------
   writeword(), readword() - Write and read words
   ------------------------------------------------------------------ */

static void writeword(FILE *f, long w, poly_t *p)

{
    int i;
    if (p == NULL)
    {
        fprintf(f,"[%ld,%ld,-1]",w,zfl);
	return;
    }
    fprintf(f,"[%ld,%ld,%ld",w,(long)p->fl,(long)p->deg);
    for (i = 0; i <= p->deg; ++i)
	fprintf(f,",%ld",zftoi(p->buf[i]));
    fprintf(f,"]");
}

static int readword(long *w, poly_t **p)

{
    long fl, deg;  
    int i;

    if (!match(" [")) errexit(ERR_FILEFMT,fn_info());
    *w = readint();
    if (!match(",")) errexit(ERR_FILEFMT,fn_info());
    fl = readint();
    if (!match(",")) errexit(ERR_FILEFMT,fn_info());
    deg = readint();
    if (deg == -1)
	*p = NULL;
    else
    {
        *p = polalloc(fl,deg);
    	for (i = 0; i <= deg; ++i)
    	{
            if (!match(",")) errexit(ERR_FILEFMT,fn_info());
	    (*p)->buf[i] = zitof(readint());
    	}
    }
    if (!match("]")) errexit(ERR_FILEFMT,fn_info());
    return 1;
}


/* ------------------------------------------------------------------
   readcfinfo() - Read the .cfinfo file
   ------------------------------------------------------------------ */

void readcfinfo()

{
    int i;

    f = myfopen(fn_info(),FM_READ|FM_TEXT);
    ncf = 0;
    ngen = 2;
    nextline();
    if (!match("CFInfo := rec();\n"))
	errexit(ERR_FILEFMT,fn_info());
    while (!feof(f))
    {
	nextline();
	if (match("CFInfo.NCF := "))
	    sscanf(bp,"%d",&ncf);
	else if (match("CFInfo.Field := "))
	{   long fl = scanint();
	    zsetfield(fl);
	}
	else if (match("CFInfo.NGen := "))
	    sscanf(bp,"%d",&ngen);
	else if (match("CFInfo.Dimension := ["))
	    for (i = 0; i < ncf; ++i) cfinfo[i].dim = scanint();
	else if (match("CFInfo.Number := ["))
	    for (i = 0; i < ncf; ++i) cfinfo[i].num = scanint();
	else if (match("CFInfo.Multiplicity := ["))
	    for (i = 0; i < ncf; ++i) cfinfo[i].mult = scanint();
	else if (match("CFInfo.SplittingField := ["))
	    for (i = 0; i < ncf; ++i) cfinfo[i].spl = scanint();
	else if (match("CFInfo.IdWord := ["))
	    for (i = 0; i < ncf; ++i)
	    {
		readword(&(cfinfo[i].idword),&(cfinfo[i].idpol));
		if (!match(i < ncf - 1 ? "," : "];"))
			errexit(ERR_FILEFMT,fn_info());
	    }
	else if (match("CFInfo.PeakWord := ["))
	    for (i = 0; i < ncf; ++i)
	    {
		readword(&(cfinfo[i].peakword),&(cfinfo[i].peakpol));
		if (!match(i < ncf - 1 ? "," : "];"))
			errexit(ERR_FILEFMT,fn_info());
	    }
	else if (match("CFInfo.NMountains := ["))
	    for (i = 0; i < ncf; ++i) cfinfo[i].nmount = scanint();
	else if (match("CFInfo.NDottedLines := ["))
	    for (i = 0; i < ncf; ++i) cfinfo[i].ndotl = scanint();
	else
	    errexit(ERR_FILEFMT,fn_info());
    }
    fclose(f);
    if (MSG1)
	printf("Read %s: %d composition factors\n",fn_info(),ncf);
}


/* ------------------------------------------------------------------
   writecfinfo() - Write the .cfinfo file
   ------------------------------------------------------------------ */

void writecfinfo()

{
    FILE *f;
    int i;

    f = myfopen(fn_info(),FM_CREATE|FM_TEXT);

    fprintf(f,"CFInfo := rec();\n");
    fprintf(f,"CFInfo.NGen := %d;\n",ngen);
    fprintf(f,"CFInfo.Field := %ld;\n",zfl);
    fprintf(f,"CFInfo.NCF := %d;\n",ncf);
    fprintf(f,"CFInfo.Dimension := [");
    for (i = 0; i < ncf; ++i)
        fprintf(f,"%ld%s",cfinfo[i].dim,i == ncf-1 ? "];\n":",");
    fprintf(f,"CFInfo.Number := [");
    for (i = 0; i < ncf; ++i)
        fprintf(f,"%ld%s",cfinfo[i].num,i == ncf-1 ? "];\n":",");
    fprintf(f,"CFInfo.Multiplicity := [");
    for (i = 0; i < ncf; ++i)
        fprintf(f,"%ld%s",cfinfo[i].mult,i == ncf-1 ? "];\n":",");
    fprintf(f,"CFInfo.SplittingField := [");
    for (i = 0; i < ncf; ++i)
        fprintf(f,"%ld%s",cfinfo[i].spl,i == ncf-1 ? "];\n":",");
    fprintf(f,"CFInfo.NMountains := [");
    for (i = 0; i < ncf; ++i)
        fprintf(f,"%ld%s",cfinfo[i].nmount,i == ncf-1 ? "];\n":",");
    fprintf(f,"CFInfo.NDottedLines := [");
    for (i = 0; i < ncf; ++i)
        fprintf(f,"%ld%s",cfinfo[i].ndotl,i == ncf-1 ? "];\n":",");

    fprintf(f,"CFInfo.PeakWord := [");
    for (i = 0; i < ncf; ++i)
    {
        writeword(f,cfinfo[i].peakword,cfinfo[i].peakpol);
	if (i < ncf-1) fprintf(f,",");
    }
    fprintf(f,"];\n");

    fprintf(f,"CFInfo.IdWord := [");
    for (i = 0; i < ncf; ++i)
    {
        writeword(f,cfinfo[i].idword,cfinfo[i].idpol);
	if (i < ncf-1) fprintf(f,",");
    }
    fprintf(f,"];\n");

    fclose(f);
}



/* ------------------------------------------------------------------
   cfname() - Returns the name of a composition factor.
   ------------------------------------------------------------------ */

char *cfname(cf)
int cf;

{
    char buf[9];
    int num = (char) (cfinfo[cf].num);

    if (num < 26)
	sprintf(buf,"%c",(char)num+'a');
    else if (num < 26*26)
	sprintf(buf,"%c%c",(char)(num%26)+'a',(char)(num/26-1)+'a');
    else
	sprintf(buf,"cf%d",num);
    sprintf(fnbuf,"%ld%s",cfinfo[cf].dim,buf);
    return fnbuf;
}

