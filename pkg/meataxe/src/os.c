/* ========================== C MeatAxe =============================
   os.c - OS dependent stuff

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */

/* $Id: os.c,v 1.2 1997/09/11 15:43:14 gap Exp $
 *
 * $Log: os.c,v $
 * Revision 1.2  1997/09/11 15:43:14  gap
 * New version 2.2.3. AH
 *
 * Revision 1.16  1995/02/08  10:14:52  mringe
 * _PL entfernt.
 *
 * Revision 1.16  1995/02/08  10:14:52  mringe
 * _PL entfernt.
 *
 * Revision 1.15  1994/11/28  16:38:00  mringe
 * Neue Namen: SFOpen() und SFSeek().
 *
 * Revision 1.14  1994/11/25  13:57:13  mringe
 * Neue Namen fuer os_clock(), timeused(), ...; ANSI-C
 *
 * Revision 1.13  1994/07/04  11:47:00  mringe
 * Bug in timeused() behoben.
 *
 * Revision 1.12  1994/06/24  12:20:05  mringe
 * Prototyp fuer setitimer().
 *
 * Revision 1.11  1994/06/22  15:44:49  mringe
 * Prototypen fuer SunOS.
 *
 * Revision 1.10  1994/06/22  14:38:40  mringe
 * Korrekturen fuer HP/UX.
 *
 * Revision 1.9  1994/06/16  14:17:36  mringe
 * timeused benutzt jetzt 1/10 Sekunden. 'os_' bei einigen
 * Funktionen entfernt.
 *
 * Revision 1.8  1994/06/14  12:06:34  mringe
 * Neue OS_xxxx Konstanten. Prototypen in meataxe.h
 *
 * Revision 1.7  1994/05/20  12:41:25  mringe
 * Prototypen.
 *
 * Revision 1.6  1994/05/19  11:34:23  mringe
 * Neu: Smalloc() und Srealloc().
 *
 * Revision 1.5  1994/05/10  09:42:51  mringe
 * SFOpen() setzt jetzt mtxerrno.
 *
 * Revision 1.4  1994/04/07  16:17:00  mringe
 * Deklaration ov getrusage fuer Ultrix eingefuegt.
 *
 * Revision 1.3  1994/02/15  11:08:20  mringe
 * *** empty log message ***
 *
 * Revision 1.1  1994/02/13  18:26:56  mringe
 * Initial revision
 *
 */
 


/* ------------------------------------------------------------------
   Include files
   ------------------------------------------------------------------ */

#if defined(OS_HPUX)

#  define _INCLUDE_POSIX_SOURCE
#  include <sys/types.h>
#  define _INCLUDE_HPUX_SOURCE
#  include <sys/time.h>
#  include <signal.h>
#  undef _INCLUDE_HPUX_SOURCE
#  include <unistd.h>
#  include <sys/times.h>
#  undef _INCLUDE_POSIX_SOURCE

#else

#  include <sys/time.h>
#  include <sys/times.h>
#  include <signal.h>
#  include <unistd.h>

#endif


#include <sys/resource.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdarg.h>
#include "meataxe.h"

#if defined(OS_IBMVM)
#   include <ctype.h>
#else
#   include <unistd.h>
#endif


/* ------------------------------------------------------------------
   Add some missing declarations
   ------------------------------------------------------------------ */

#if !defined(SIGFUNTYPE)
#define SIGFUNTYPE void
#endif

#if defined(PROTO_GETRUSAGE)
extern int getrusage(int who, struct rusage *r);
#endif

#if defined(PROTO_CLOCK)
extern long clock(void);
#endif

#if defined(PROTO_SETITIMER)
extern int setitimer(int which, struct itimerval *val, struct
	itimerval *oval);
#endif



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

#if defined(OS_MSDOS)
time_t zinittime = 0;
#endif


/* ------------------------------------------------------------------
   fopen() modes
   ------------------------------------------------------------------ */

static char *fmodes[6] = { 
#if !defined(FMODES)
 "at","rt","wt","ab","rb","wb"
#else
 FMODES
#endif
 };



/* ------------------------------------------------------------------
   SInit() - Initialize OS dependent stuff. Called from mtxinit().
   ------------------------------------------------------------------ */

void SInit()

{
#if defined(NO_GETRUSAGE) && defined(NO_TIMES)
    zinittime = time(NULL);
#endif
}


/* ------------------------------------------------------------------
   SClock() - Returns CPU time used (in system dependent units).
   ------------------------------------------------------------------ */

long SClock()
{
    return (long) clock();
}


/* ------------------------------------------------------------------
   STimeUsed() - Returns CPU time used in 1/10 seconds.
   ------------------------------------------------------------------ */

/* Case 1: Systems without CPU time information
   -------------------------------------------- */

#if defined(NO_GETRUSAGE) && defined(NO_TIMES)

long STimeUsed(void)
{
    return (time(NULL) - zinittime) * 10;
}

/* Case 2: Systems without getrusage(). Use sysconf()
   and times().
   ---------------------------------------------------- */

#elif defined(NO_GETRUSAGE)


long STimeUsed(void)

{
    struct tms t;
    static long clk_tck = 0;

    if (clk_tck == 0) clk_tck = sysconf(_SC_CLK_TCK);
    times(&t);
    return ((long)((t.tms_utime + t.tms_stime) * 10 / clk_tck ));
}


/* Case 3: Systems with getrusage()
   -------------------------------- */

#else

long STimeUsed(void)

{
    static struct rusage ru;
    getrusage(RUSAGE_SELF,&ru);
    return (ru.ru_utime.tv_sec * 10 + ru.ru_utime.tv_usec / 100000);
}

#endif



/* ------------------------------------------------------------------
   STimeLimit() - Set CPU time limit in seconds.
   ------------------------------------------------------------------ */

#if !defined(NO_ITIMER)


SIGFUNTYPE vtalarm(int i)

{
    errexit(ERR_GAMEOVER,"I'm sorry, but");
#if (SIGFUNTYPE != void)
    return (SIGFUNTYPE)0;
#endif
}

void STimeLimit(long nsecs)

{
    struct itimerval tv;

    tv.it_interval.tv_sec = 0;
    tv.it_interval.tv_usec = 0;
    tv.it_value.tv_sec = nsecs;
    tv.it_value.tv_usec = 0;
    setitimer(ITIMER_VIRTUAL,&tv,NULL);
    signal(SIGVTALRM,vtalarm);
}

#else /* No interval timer, sorry... */

void STimeLimit(long nsecs)

{
}

#endif


/* ------------------------------------------------------------------
   os_mkfilename() - Convert any string into a valid file name.
   ------------------------------------------------------------------ */

char *os_mkfilename(name)
char *name;

{

#if defined(OS_MSDOS)

    static char buf[30];
    char *d;
    int i, dp, nd = 0;

    for (i = 0; name[i] != 0; ++i)
	if (name[i] == '.')
	    if (++nd == 1) dp = i;
    if (nd == 1 && dp <= 8 && i <= 12)
	return name;
    d = buf+12;
    *d = 0;
    for (--i; d != buf && i >= 0; --i)
    {
	if (name[i] != '.')
	{
	    if (d == buf+8) *--d = '.';
	    *--d = name[i];
	}
    }
    while (d != buf) *--d = 'x';
    return buf;

#elif defined(OS_IBMVM)

    static char buf[30];
    char *d = buf;
    int i;

    for (i = 0; name[i] != 0 && name[i] != '.'; ++i);
    if (i > 8) i -= 8; else i = 0;
    while (name[i] != 0 && name[i] != '.'; ++i)
    {
	if (isalnum(name[i]))
	    *d++ = name[i];
	else
	    *d++ = '-';
    }
    *d++ = ' ';
    if (name[i] == '.')
    {
	int i0 = i+1;
        for (i = i0; name[i] != 0; ++i);
        if (i > 8) i -= 8; else i = i0;
        while (name[i] != 0; ++i)
	{
	    if (isalnum(name[i]))
	        *d++ = name[i];
	    else
		*d++ = '-';
	}
    }
    else 
    {
	strcpy(d,"MTX");
	d += 3;
    }
    strcpy(d," A");
    return buf;

#else	/* Default: no translation */

    return name;

#endif

}


/* ------------------------------------------------------------------
   SFOpen() - OS independent fopen()
   ------------------------------------------------------------------ */

FILE *SFOpen(name, mode)
char *name;
int mode;

{
    char *canonical_name = os_mkfilename(name);
    char buf[200];
    int m;
    FILE *f;
    char *c;

    m = mode & 0x0F;			/* Append, read or create? */
    if ((mode & FM_TEXT) == 0) m += 3;	/* Binary mode */
    if (m < 0 || m >= 6)
	MTXFAIL(ERR_BADARG,NULL);		/* Invalid mode */
    f = fopen(canonical_name,fmodes[m]);
    if (f != NULL) return f;

    /* Search library directory
       ------------------------ */
    if ((mode & FM_LIB) == 0) MTXFAIL(ERR_FILEOPEN,NULL);
    if ((c = getenv("MTXLIB")) != NULL)
    {	strcpy(buf,c);
	strcat(buf,canonical_name);
	if ((f = fopen(buf,fmodes[m])) != NULL)
	    return f;
    }
    strcpy(buf,MTXLIB);
    strcat(buf,canonical_name);
    if ((f = fopen(buf,fmodes[m])) != NULL)
	return f;
    MTXFAIL(ERR_FILEOPEN,NULL);
}


/* ------------------------------------------------------------------
   SFSeek() - OS independent fseek()
   pos >= 0: absolute seek
   pos < 0: seek to end of file
   ------------------------------------------------------------------ */

#if !defined(SEEK_SET)
#define SEEK_SET 0
#endif
#if !defined(SEEK_END)
#define SEEK_END 2
#endif

int SFSeek(file,pos)
FILE *file;
long pos;

{
    if (pos < 0)
	return fseek(file,(long) 0,SEEK_END);
    else
	return fseek(file,pos,SEEK_SET);
}



/* ------------------------------------------------------------------
   Smalloc() - malloc() wrapper.
   Srealloc() - realloc() wrapper.
   
   These functions avoid allocation of zero length segments, because
   this is not portable.
   ------------------------------------------------------------------ */
   

void *Smalloc(nbytes)
size_t nbytes;

{
    if (nbytes == 0) nbytes = 1;
    return malloc(nbytes);
}


void *Srealloc(buf,nbytes)
void *buf;
size_t nbytes;

{
    if (nbytes == 0) nbytes = 1;
    return realloc(buf,nbytes);
}


