/* ========================== C MeatAxe =============================
   args.c - Command line options and arguments.

   (C) Copyright 1993 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */


/* $Id: args.c,v 1.2 1997/09/11 15:42:31 gap Exp $
 *
 * $Log: args.c,v $
 * Revision 1.2  1997/09/11 15:42:31  gap
 * New version 2.2.3. AH
 *
 * Revision 2.5  1995/05/12  10:07:20  mringe
 * Neue Variable MeatAxeBinDir.
 *
 * Revision 2.4  1994/11/25  13:57:13  mringe
 * Neue Namen fuer os_clock(), timeused(), ...; ANSI-C
 *
 * Revision 2.3  1994/06/19  15:59:17  mringe
 * Fix sscanf bug (?) in HP/UX.
 *
 * Revision 2.2  1994/06/16  14:18:34  mringe
 * os_timelimit() in timelimit() umbenannt.
 *
 * Revision 2.1  1994/02/13  18:26:56  mringe
 * Neu: os.c, os.h.
 *
 * Revision 2.0  1993/10/14  18:54:18  mringe
 * MeatAxe-2.0, Phase I
 *
 * Revision 1.12  1993/10/11  19:05:28  mringe
 * Neue Library-Struktur.
 *
 * Revision 1.11  1993/10/06  04:41:05  mringe
 * utils Library eliminiert.
 *
 * Revision 1.10  1993/08/06  14:01:59  mringe
 * Neuer File-header.
 *
 * Revision 1.9  1993/06/14  19:23:32  mringe
 * Tippfehler korrigiert.
 *
 * Revision 1.8  1993/06/14  19:20:52  mringe
 * Dokumentation
 *
 * Revision 1.7  1993/02/17  11:16:12  mringe
 * Include-Files...
 *
 * Revision 1.6  1993/01/28  07:35:51  mringe
 * In Bibliothek integriert.
 *
 * Revision 1.5  1993/01/06  20:59:57  mringe
 * getopt in zgetopt() umbenannt.
 *
 * Revision 1.4  1992/07/28  08:16:47  mringe
 * getint() verbessert.
 *
 * Revision 1.3  1992/07/28  08:09:05  mringe
 * getint() verbessert
 *
 * Revision 1.2  1992/07/28  08:06:47  mringe
 * getint() verbessert
 *
 * Revision 1.1  1992/07/15  09:33:16  mringe
 * Initial revision
 */

#include "meataxe.h"
#include <ctype.h>
#include <string.h>
#include <stdlib.h>



/* ------------------------------------------------------------------
   Global data
   ------------------------------------------------------------------ */

char opt_char;
char opt_text[50];
char *opt_text_ptr;
int opt_ind = 0;
char MeatAxeBinDir[250] = "";

static int origargc;
static char **origargv;
static proginfo_t *pinfo;


/* ------------------------------------------------------------------
   initargs()
   ------------------------------------------------------------------ */

void initargs(int argc, char **argv, proginfo_t *pi)

{
    char hrev[20];
    char *c;

    origargc = argc;
    origargv = argv;
    pinfo = pi;
#if defined(OS_HPUX)	/* sscanf broken in 9.0 ????? */
    { char tmp[20]; strncpy(tmp,pi->rcsrev,sizeof(tmp));
      sscanf(tmp,"$Revision: %s",hrev);
    }
#else
    sscanf(pi->rcsrev,"$Revision: %s",hrev);
#endif

    /* Setze MeatAxeBinDir
       ------------------- */
    if ((c = getenv("MTXBIN")) != NULL)
	strcpy(MeatAxeBinDir,c);
    else
    {
	strcpy(MeatAxeBinDir,argv[0]);
    	for (c = MeatAxeBinDir; *c != 0; ++c);
    	while (c != MeatAxeBinDir && *c != '/' && *c != '\\') --c;
    	*c = 0;
    }

    /* Werte die Option -help aus
       -------------------------- */
    if ((argc == 2 && !strcmp(argv[1],"help"))
	|| (argc >= 2 && !strcmp(argv[1],"-help")))
    {
	char **ht = pi->helptext;

	printf("NAME\n    %s - %s\n    Revision %s\n\n",
	  pi->name,pi->shortdesc,hrev);
	if (ht != NULL)
	    for (; *ht != NULL; ++ht)
	   	 printf("%s\n",*ht);
	else
	    printf("Sorry, no help available\n");
	exit(EXIT_OK);
    }
}


/* ------------------------------------------------------------------
   zgetopt()
   ------------------------------------------------------------------ */

int zgetopt(pattern)
char *pattern;

{
    static char pat[50];
    static char *c = NULL;
    char *d;

    strcpy(pat,pattern);
    strcat(pat,"QVT:");

    while (1)
    {
    	opt_char = 0;
    	opt_text[0] = 0;
    	opt_text_ptr = opt_text;
    	if (c == NULL)
    	{
	    if (++opt_ind >= origargc || *(c = origargv[opt_ind]) !=
		'-' || *c == 0)
	    	return OPT_END;		/* No more options */
	    c = origargv[opt_ind] + 1;
    	}
    	opt_char = *c;
    	for (d = pat; *d != 0 && *d != opt_char; ++d);
    	if (*d != opt_char)
	    errexit(ERR_BADUSAGE,pinfo->name);
    	if (d[1] == ':')		/* Option needs an argument */
    	{
	    if (*++c == 0)
	    {
		if (++opt_ind >= origargc)
		    errexit(ERR_BADUSAGE,pinfo->name);
		c = origargv[opt_ind];
	    }
	    strcpy(opt_text,c);
	    c = NULL;
    	}
    	else
    	{
	    if (*++c == 0) c = NULL;
    	}
    
        /* Check standard options
           ---------------------- */
        if (opt_char == 'Q')
	   msg_level = -1000;
        else if (opt_char == 'V')
	   ++msg_level;
        else if (opt_char == 'T')
        {
	    long tl = getint();
	    if (tl < 0) errexit(ERR_OPTION,"-T");
	    STimeLimit(tl);
        }
        else
            return opt_char;
    }
}


/* ------------------------------------------------------------------
   getint()
   ------------------------------------------------------------------ */

long getint()

{   long k;

    if (!isdigit(*opt_text_ptr))
	return GETINT_ERR;
    for (k = 0; isdigit(*opt_text_ptr); ++opt_text_ptr)
    {   switch (*opt_text_ptr)
	{   case '0': k = k * 10 + 0; break;
	    case '1': k = k * 10 + 1; break;
	    case '2': k = k * 10 + 2; break;
	    case '3': k = k * 10 + 3; break;
	    case '4': k = k * 10 + 4; break;
	    case '5': k = k * 10 + 5; break;
	    case '6': k = k * 10 + 6; break;
	    case '7': k = k * 10 + 7; break;
	    case '8': k = k * 10 + 8; break;
	    case '9': k = k * 10 + 9; break;
	}
    }
    return k;
}



