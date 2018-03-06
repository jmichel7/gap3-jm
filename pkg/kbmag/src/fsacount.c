/* fsacount.c 27/10/94.
 * Counts the accepted language of a fsa
 * This outputs to stdout only.
 *
 * SYNOPSIS: fsacount [-ip d/s] [-silent] [-v] [filename]\n");
 *
 * Input is from filename or stdin, which should contain a fsa.
 *
 * OPTIONS:
 * -ip d/s      input in dense or sparse format - dense is default
 * -v           verbose
 * -silent      no diagnostics
 */

#include <stdio.h>
#include "defs.h"
#include "fsa.h"
#include "definitions.h"

FILE *rfile,
     *wfile,
     *fopen();

void  badusage();

/* Functions defined in other files used in this file */
void  fsa_read();
int   fsa_count();
void  fsa_clear();

main(argc, argv)
        int             argc;
        char           *argv[];
{ int arg, ct;
  fsa testfsa;
  char inf[100],fsaname[100];
  storage_type ip_store = DENSE;

  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  inf[0] = '\0';
  arg = 1;
  while (argc > arg) {
    if (strcmp(argv[arg],"-ip")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      if (strcmp(argv[arg],"d")==0)
        ip_store = DENSE;
      else if (argv[arg][0] == 's')
        ip_store = SPARSE;
      else
        badusage();
    }
    else if (strcmp(argv[arg],"-silent")==0)
      print_level = 0;
    else if (strcmp(argv[arg],"-v")==0)
      print_level = 2;
    else if (strcmp(argv[arg],"-vv")==0)
      print_level = 3;
    else {
       if (argv[arg][0] == '-')
         badusage();
       if (strcmp(inf,""))
         badusage();
       strcpy(inf,argv[arg]);
    }
    arg++;
  }
  if (strlen(inf)==0)
    rfile = stdin;
  else if ((rfile = fopen(inf,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf);
    exit(1);
  }
  fsa_read(rfile,&testfsa,ip_store,0,0,TRUE,fsaname);
  if (strlen(inf))
    fclose(rfile);

  ct = fsa_count(&testfsa);
  if (ct == -1)
    printf("#The language accepted is infinite.\n");
  else
    printf("#The language accepted has size %d.\n",ct);

  fsa_clear(&testfsa);
  exit(0);
}
 
void
badusage()
{
    fprintf(stderr,
      "Usage: fsacount [-ip d/s] [-silent] [-v] [filename]\n");
    exit(1);
}
