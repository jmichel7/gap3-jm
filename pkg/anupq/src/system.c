#include "pq_defs.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/times.h>

/* system and operating system dependent pieces of code */

/* return CPU time in CLOCK TICKS -- the program should report 
   correct CPU times for each of SunOS and Solaris if compiled 
   and run under that operating system; under Solaris, 
   CLK_TCK is defined in <limits.h>; if compiled under SunOS
   and run under Solaris, then multiply reported times by
   3/5 to get correct user time */

int runTime ()
{
   struct tms buffer;

   times(&buffer);
   return buffer.tms_utime + buffer.tms_cutime;
}

/* print startup message */

int print_message (int work_space)
{
   time_t now;
   char *id;
   char string[100];

#if defined (HOSTNAME) 
   char s[100];
#else 
   char *s;
#endif

#ifdef VMS
   int time_val;
   struct tm *x;
   time (&time_val);
   x = localtime (&time_val);
   printf ("ANU p-Quotient Program Version 1.4 running with workspace %d\n", 
	   work_space);
   printf ("Date is %d:%d:%d %d-%d-19%d\n", x->tm_hour, x->tm_min, x->tm_sec,  
	   x->tm_mday, x->tm_mon + 1, x->tm_year);
#else 

#if defined (HOSTNAME) 
   gethostname (s, 100);
#else 
   s = (char *) getenv ("HOST");
   if (s == NULL) s = "unknown";
#endif 

#if defined (LIE) 
   id = "Lie Ring Program Version 1.0";
#endif
#if defined (GROUP) 
   id = "p-Quotient Program Version 1.4";
#endif 

   printf ("ANU %s running with workspace %d on %s\n", 
	   id, work_space, s);
   now = time (NULL);
#ifndef HAS_NO_STRFTIME
   strftime (string, 100, "%a %b %d %H:%M:%S %Z %Y", localtime (&now));
   printf ("%s\n", string);
#else
   printf ("\n");
#endif
#endif
}
