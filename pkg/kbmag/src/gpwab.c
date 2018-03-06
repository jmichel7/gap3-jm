/* gpwab.c 7/11/94.
 * Calculates word acceptor of an automatic group.
 * This is now redundant, since gpwa seems to be faster.
 * The input is a word-difference fsa for the group
 * (from a run of kbprog with the -wd option).
 *
 * SYNOPSIS:  gpwa [-op d/s] [-silent] [-v] [-l/-h] [-diff1/-diff2] groupname
 *
 * Input is from groupname.diff1 or groupname.diff2
 * Output is to groupname.wab
 *
 * OPTIONS:
 * -op d/s      output in dense or sparse format - dense is default
 * -v           verbose
 * -silent      no diagnostics
 * -l/-h        large/huge hash-tables (for big examples)
 * -diff1/diff2 take input from groupname.diff1 or groupname.diff2
 *              (former is default). Latter needs more space and is
 *              slower, but can sometimes save time later.
 *
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
fsa  *fsa_greater_than();
fsa  *fsa_and();
fsa  *fsa_nosub_exists();
void  fsa_print();
void  fsa_clear();

main(argc, argv)
        int             argc;
        char           *argv[];
{ int arg;
  fsa *fsawd, *fsagt, *diffgt, *gpwab;
  char inf[100], outf[100], fsaname[100], tempfilename[100];
  storage_type op_store = DENSE;
  boolean diff1_ip = TRUE;

  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  inf[0] = '\0';
  arg = 1;
  while (argc > arg) {
    if (strcmp(argv[arg],"-op")==0) {
      arg++;
      if (arg >= argc)
        badusage();
      if (strcmp(argv[arg],"d")==0)
        op_store=DENSE;
      else if (strcmp(argv[arg],"s")==0)
        op_store=SPARSE;
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
    else if (strcmp(argv[arg],"-diff1")==0)
      diff1_ip = TRUE;
    else if (strcmp(argv[arg],"-diff2")==0)
      diff1_ip = FALSE;
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
    badusage();

  strcpy(outf,inf);
  strcat(outf,".wab");

  strcpy(tempfilename,inf);
  strcat(tempfilename,"temp_wab_XXX");

  if (diff1_ip)
    strcat(inf,".diff1");
  else
    strcat(inf,".diff2");

  tmalloc(fsawd,fsa,1);

/* We always use dense format for the word-difference machine - this
 * is much faster, and the machine is usually not too big.
 */
  if ((rfile = fopen(inf,"r")) == 0) {
    fprintf(stderr,"Cannot open file %s.\n",inf);
    exit(1);
  }
  fsa_read(rfile,fsawd,DENSE,0,0,TRUE,fsaname);
  fclose(rfile);

  if (fsawd->alphabet->type != PRODUCT || fsawd->alphabet->arity != 2) {
    fprintf(stderr,"Error: fsa must be 2-variable.\n");
    exit(1);
  }
  
  fsagt = fsa_greater_than(fsawd->alphabet->base);

  diffgt = fsa_and(fsawd,fsagt,DENSE,TRUE,tempfilename);
  tfree(fsawd) tfree(fsagt)
  if (print_level>1)
    printf("  #Number of states of diffgt = %d.\n",diffgt->states->size);

  fsa_minimize(diffgt);
  if (print_level>1)
    printf("  #Number of states of diffgt after minimization = %d.\n",
      diffgt->states->size);

  gpwab = fsa_nosub_exists(diffgt,op_store,TRUE,tempfilename);
 
  tfree(diffgt);
  if (print_level>1)
    printf("  #Number of states of gpwab before minimisation = %d.\n",
        gpwab->states->size);
  fsa_minimize(gpwab);
  if (print_level>1)
     printf("  #Number of states of gpwab after minimisation = %d.\n",
        gpwab->states->size);

  base_prefix(fsaname);
  strcat(fsaname,".wab");
  wfile = fopen(outf,"w");
    fsa_print(wfile,gpwab,fsaname);
  fclose(wfile);

  if (print_level>0)
    printf("#Word-acceptor with %d states computed.\n",gpwab->states->size);

  fsa_clear(gpwab);
  tfree(gpwab);
  exit(0);
}

void
badusage()
{
    fprintf(stderr,
   "Usage: gpwab [-op d/s] [-silent] [-v] [-l/-h] [-diff1/-diff2] groupname\n");
    exit(1);
}
