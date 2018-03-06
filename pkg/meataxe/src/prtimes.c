/* ========================== C MeatAxe =============================
   prtimes.c - CPU time functions.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: prtimes.c,v 1.2 1997/09/11 15:43:17 gap Exp $
 *
 * $Log: prtimes.c,v $
 * Revision 1.2  1997/09/11 15:43:17  gap
 * New version 2.2.3. AH
 *
 * Revision 2.6  1994/11/25  13:57:13  mringe
 * Neue Namen fuer os_clock(), timeused(), ...; ANSI-C
 *
 * Revision 2.5  1994/06/16  14:19:35  mringe
 * Gibt die Zeit bis auf 1/10 Sekunde an.
 *
 * Revision 2.4  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.3  1993/12/08  12:00:15  mringe
 * Argumente fuer signal handler korrigiert.
 *
 * Revision 2.2  1993/12/08  11:49:35  mringe
 * Time limit voellig neu geschrieben.
 *
 * Revision 2.1  1993/11/25  13:52:02  mringe
 * Verwende _NO_CYSCONF und HAS_UNSTD_H.
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
 * Revision 1.14  1993/10/05  14:17:48  mringe
 * GETRU
 *
 * Revision 1.13  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.12  1993/08/05  15:48:54  mringe
 * Neues message.c
 *
 * Revision 1.11  1993/07/29  14:47:00  mringe
 * Benutze getrusage() ODER times(), abh. vom OS.
 *
 * Revision 1.10  1993/07/28  13:34:49  mringe
 * Benutze sysconf.
 *
 * Revision 1.9  1993/07/23  13:46:27  mringe
 * OS-Symbole neu (SYS_xxx)
 *
 * Revision 1.8  1993/07/23  12:54:52  mringe
 * Bug behoben: Default timelimit ist jetzt unendlich.
 *
 * Revision 1.7  1993/07/18  15:28:15  mringe
 * timeused() fuer SYS_MSDOS.
 *
 * Revision 1.6  1993/07/17  12:52:05  mringe
 * Neu: time_setlimit() und time_checklimit().
 *
 * Revision 1.5  1993/02/16  18:32:46  mringe
 * string.h und stdio.h werden jetzt in meataxe.h included.
 *
 * Revision 1.4  1993/02/16  17:33:20  mringe
 * Symbole SYS_BSD und SYS_SYSV eingefuehrt.
 *
 * Revision 1.3  1993/01/28  07:35:51  mringe
 * In Bibliothek integriert.
 *
 * Revision 1.2  1993/01/22  12:16:28  hiss
 * NOPRTIMES f"ur Maschinen ohne sysconf ...
 *
 * Revision 1.1  1992/05/25  17:29:24  mringe
 * Initial revision
 */
 

#include "meataxe.h"


/* ------------------------------------------------------------------
   prtimes() - Print the CPU time to stdout.
   ------------------------------------------------------------------ */


void prtimes()

{
    long t = STimeUsed();
    printf("Time used: %ld.%ld seconds\n",t/10,t%10);
    fflush(stdout);
}

