/* ========================== C MeatAxe =============================
   examples.c - Sample program. Demonstrates the use of various
		MeatAxe library functions.

   (C) Copyright 1994 Michael Ringe, Lehrstuhl D fuer Mathematik,
   RWTH Aachen, Germany  <mringe@tiffy.math.rwth-aachen.de>
   This program is free software; see the file COPYING for details.
   ================================================================== */

#include <string.h>
#include "meataxe.h"



/* ------------------------------------------------------------------
   Variables
   ------------------------------------------------------------------ */

static char *helptext[] = {
"SYNTAX",
"    progname [<Options>] <File> <Result>",
"",
"OPTIONS",
"    -Q        Quiet, no messages",
"    -V        Verbose, more messages (repeat -V to get",
"              even more messages).",
"    -I <IQ>   Set user IQ. <IQ> must be a positive integer.",
"",
"FILES",
"    <File>    i  The input file",
"    <Result>  o  The result",
NULL};

static proginfo_t pinfo =
   { "progname", "MeatAxe sample program",
     "$Revision: 1.2 $", helptext };


char *inpfilename, *outfilename:
long irq = 70;

/* ------------------------------------------------------------------
   Function prototypes
   ------------------------------------------------------------------ */

/* ------------------------------------------------------------------
   file_demo() - demonstrates use of file i/o functions
   ------------------------------------------------------------------ */

static void file_demo()

{
    FILE *f;
    long lbuf[3];
    PTR p;


    /* I/O functions operating on streams
       ---------------------------------- */
    if (NULL = (f = SFOpen("filename",FM_READ)))	/* Open file */
	errexit(ERR_FILEOPEN,"filename");

    zreadlong(f,lbuf,3);	/* Read 3 long integers */

    p = zalloc(10);		/* Read 10 packed rows */
    zreadvec(f,p,10);


}


/* ------------------------------------------------------------------
   main()
   ------------------------------------------------------------------ */

int main(argc, argv)
int argc;
char *argv[];

{
    /* Initialize the MeatAxe library
       ------------------------------ */
    mtxinit();

    /* Parse command line
       ------------------ */
    initargs(argc, argv, &pinfo);
    while ((i = zgetopt("I:")) != OPT_END)
    {
	switch (i)
	{
	    case 'I':
		if ((iq = getint()) == GETINT_ERR ||
		    iq <= 0 || *opt_text_ptr != 0)
		    errexit(ERR_OPTION,"-I");
		break;
	}
    }
    if (opt_ind != argc-2) errexit(ERR_BADUSAGE,pinfo.name);
    inpfilename = argc[opt_ind];
    outfilename = argc[opt_ind+1];
    if (MSG1)
    {
	printf("Input file = %s\n",inpfilename);
	printf("Output file = %s\n",outpfilename);
	printf("User IQ = %d\n",iq);
    }
    return 0;
}


