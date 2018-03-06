/* fsaand.c 3/11/94.
 * calculates "and" of two finite state automata
 *
 * SYNOPSIS:  fsaand [-ip d/s[dr]] [-op d/s] [-silent] [-v] [-l/-h]
                filename1 filename2 outfilename
 *
 * Input is from filename1 and filename2, which should contain fsa's.
 * Output is to outfilename
 *
 * OPTIONS:
 * -ip d/s[dr]  input in dense or sparse format - dense is default
 * -op d/s      output in dense or sparse format - default is as in current
 *               value of table->printing_format, in the fsa.
 * -v           verbose
 * -silent      no diagnostics
 * -l/-h        large/huge hash-tables (for big examples)
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
fsa  *fsa_and();
void  fsa_print();
void  fsa_clear();

main(argc, argv)
        int             argc;
        char           *argv[];
{ int arg;
  fsa fsa1, fsa2, *fsaand;
  char inf1[100], inf2[100], outf[100], fsaname1[100], fsaname2[100],
       tempfilename[100];
  storage_type ip_store = DENSE;
  int dr = 0;
  storage_type op_store = DENSE;


  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  inf1[0] = '\0';
  inf2[0] = '\0';
  outf[0] = '\0';
  arg = 1;
  while (argc > arg) {
    if (strcmp(argv[arg],"-ip")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      if (strcmp(argv[arg],"d")==0)
        ip_store = DENSE;
      else if (argv[arg][0] == 's') {
        ip_store = SPARSE;
        if (strlen(argv[arg]) > 1)
          dr = atoi(argv[arg]+1);
      }
      else
        badusage();
    }
    else if (strcmp(argv[arg],"-op")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      if (strcmp(argv[arg],"d")==0)
        op_store = DENSE;
      else if (strcmp(argv[arg],"s")==0)
        op_store = SPARSE;
      else
        badusage();
    }
    else if (strcmp(argv[arg],"-silent")==0)
      print_level = 0;
    else if (strcmp(argv[arg],"-v")==0)
      print_level = 2;
    else if (strcmp(argv[arg],"-vv")==0)
      print_level = 3;
    else if (strcmp(argv[arg],"-l")==0)
      large = TRUE;
    else if (strcmp(argv[arg],"-h")==0)
      huge = TRUE;
    else {
       if (argv[arg][0] == '-')
         badusage();
       if (strcmp(outf,""))
         badusage();
       if (strcmp(inf1,"")==0)
         strcpy(inf1,argv[arg]);
       else if (strcmp(inf2,"")==0)
         strcpy(inf2,argv[arg]);
       else
         strcpy(outf,argv[arg]);
    }
    arg++;
  }
  if (strlen(inf1)==0 || strlen(inf1)==0 || strlen(outf)==0)
    badusage();

  if ((rfile = fopen(inf1,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf1);
    exit(1);
  }
  fsa_read(rfile,&fsa1,ip_store,dr,0,TRUE,fsaname1);
  fclose(rfile);

  if ((rfile = fopen(inf2,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf2);
    exit(1);
  }
  fsa_read(rfile,&fsa2,ip_store,dr,0,TRUE,fsaname2);
  fclose(rfile);
  
  strcpy(tempfilename,inf1);
  strcat(tempfilename,"temp_and_XXX");
  fsaand = fsa_and(&fsa1,&fsa2,op_store,TRUE,tempfilename);

  if (print_level>1)
    printf("  #Number of states of fsaand before minimisation = %d.\n",
        fsaand->states->size);
  fsa_minimize(fsaand);
  if (print_level>1)
    printf("  #Number of states of fsaand after minimisation = %d.\n",
        fsaand->states->size);

  strcat(fsaname1,"_and_");
  strcat(fsaname1,fsaname2);
  wfile = fopen(outf,"w");
  fsa_print(wfile,fsaand,fsaname1);
  fclose(wfile);

  if (print_level>0)
    printf("#\"And\" fsa with %d states computed.\n",fsaand->states->size);

  fsa_clear(fsaand);
  tfree(fsaand);

  exit(0);
}
 
void
badusage()
{
    fprintf(stderr,
      "Usage: fsaand [-ip d/s[dr]] [-op d/s] [-silent] [-v] [-l/-h]\n");
    fprintf(stderr,
      "       filename1 filename2 outfilename\n");
    exit(1);
}
