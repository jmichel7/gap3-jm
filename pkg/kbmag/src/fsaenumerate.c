/* fsaenumerate.c 16/1/95
 * Enumerates the accepted language of a fsa
 *
 * SYNOPSIS: fsaenumerate [-ip d/s] [-bfs/-dfs] min max filename\n");
 *
 * Input is from filename, which should contain a fsa.
 * Output is to filename.enumerate, and contains a list of the accepted words
 * of the finite state automaton with lengths l satisfying min <= l <= max.
 * Note that the integers min and max must be specified.
 *
 * OPTIONS:
 * -ip d/s      input in dense or sparse format - dense is default
 * -bfs/-dfs	specifies listing according to breadth-first-search of
 *		depth-first-search, resepctively. The latter is default,
 *		and is somewhat quicker, since the bradth-first-search involves
 *		repeated calls of the procedure, for each individual length.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "fsa.h"
#include "definitions.h"

FILE *rfile,
     *wfile,
     *fopen();

void  badusage();

/* Functions defined in other files used in this file */
void  fsa_read();
boolean   fsa_enumerate();
void  fsa_clear();
boolean is_int();

main(argc, argv)
        int             argc;
        char           *argv[];
{ int arg, min, max, i;
  fsa testfsa;
  char inf[100], outf[100], fsaname[100];
  storage_type ip_store = DENSE;
  boolean bfs, minset, maxset, putcomma;
  

  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  minset = maxset = FALSE;
  bfs = FALSE;
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
    else if (strcmp(argv[arg],"-dfs")==0)
      bfs = FALSE;
    else if (strcmp(argv[arg],"-bfs")==0)
      bfs = TRUE;
    else {
       if (argv[arg][0] == '-')
         badusage();
       if (strcmp(inf,""))
         badusage();
       if (!minset) {
         if (!is_int(argv[arg])) badusage();
         min = atoi(argv[arg]);
         minset = TRUE;
       }
       else if (!maxset) {
         if (!is_int(argv[arg])) badusage();
         max = atoi(argv[arg]);
         maxset = TRUE;
       }
       else strcpy(inf,argv[arg]);
    }
    arg++;
  }
  if (strlen(inf)==0)
    badusage();
  strcpy(outf,inf);
  strcat(outf,".enumerate");

  if ((rfile = fopen(inf,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf);
    exit(1);
  }
  fsa_read(rfile,&testfsa,ip_store,0,0,TRUE,fsaname);
  fclose(rfile);
  strcat(fsaname,".words");

  wfile = fopen(outf,"w");
  fprintf(wfile,"%s = [\n",fsaname);

  putcomma = FALSE;
  if (bfs)
    for (i=min;i<=max;i++)
      putcomma = fsa_enumerate(wfile,&testfsa,i,i,putcomma) || putcomma;
  else
    fsa_enumerate(wfile,&testfsa,min,max,putcomma);

  fprintf(wfile,"\n];\n");
  fclose(wfile);

  fsa_clear(&testfsa);
  exit(0);
}
 
void
badusage()
{
    fprintf(stderr,
      "Usage: fsaenumerate [-ip d/s] [-bfs/-dfs] min max filename\n");
    exit(1);
}
