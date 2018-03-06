/* wordreduce.c 1/11/94.
 * 18/1/95 - added options to reduce words using rewriting system.
 * Reduce words using a rewriting system or a word-difference machine.
 *
 * SYNOPSIS:
 *      wordreduce [-kbprog/-diff1/-diff2/-diff1c] [-mrl maxreducelen]
 *							groupname [filename]
 *
 * Input is from groupname.kbprog, groupname.diff1, groupname.diff2 or
 * groupname.diff1c, and stdin or filename.
 * Output is to stdout or to filename.reduced.
 *
 * OPTIONS:
 * -kbprog/diff1/diff2/diff1c
 *		This determines which finite state automaton is used for
 *		the word-reduction.
 *		If the kbprog flag is called, then a rewriting system is
 *		input from groupname.kbprog.
 *		Otherwise, a woord-difference machine is input from the
 *		appropriate file.
 *		The default is groupname.kbprog if that file is present,
 *		and otherwise groupname.diff2.
 *		
 * If the optional "filename" argument is present, then the file filename should
 * contain an assignment of a list of words to be reduced
 * (e.g. words := [a^3,(a*b)^7];), and output will be a lsit of reduced words
 * to filename.reduced.
 * Otherwise words are input interactively to stdin and output to stdout.
 */

#include <stdio.h>
#include "defs.h"
#include "fsa.h"
#include "rws.h"
#include "definitions.h"

#define MAXEQNS         32767
#define MAXREDUCELEN    32767

int      maxeqns = MAXEQNS,
         maxreducelen = MAXREDUCELEN,
	 maxslowhistoryspace = 0; /*not used in this program */
boolean  maxreducelenset = FALSE;

char	*gen_word;

rewriting_system rws;
fsa	reduction_fsa;
fsa	wd_fsa;

void  badusage();

/* Functions defined in other files used in this file */
void  read_kbinput_simple();
void  fsa_read();
void  fsa_table_dptr_init();
void  read_word();
void  diff_reduce();
void  rws_reduce();
void  add_word_to_buffer();
void  printbuffer();
void  process_names();
void  fsa_clear();
   
FILE *rfile,
     *wfile,
     *fopen();

main(argc, argv)
        int             argc;
        char           *argv[];
{ int delim, arg, num_gens;
  char gpname[100], inf1[100], fsaname[100], inf2[100], outf[100],
	**names;
  boolean rws_ip, diff1_ip, diff2_ip, diff1c_ip, file_ip, open, first;
  void (*reduce_word)();

  setbuf(stdout,(char*)0);
  setbuf(stderr,(char*)0);

  gpname[0] = '\0';
  inf2[0] = '\0';
  arg = 1;
  rws_ip = diff1_ip = diff2_ip = diff1c_ip = file_ip = FALSE;
  while (argc > arg) {
    if (strcmp(argv[arg],"-kbprog")==0)
      rws_ip = TRUE;
    else if (strcmp(argv[arg],"-diff1")==0)
      diff1_ip = TRUE;
    else if (strcmp(argv[arg],"-diff2")==0)
      diff2_ip = TRUE;
    else if (strcmp(argv[arg],"-diff1c")==0)
      diff1c_ip = TRUE;
    else if (strcmp(argv[arg], "-mrl") == 0) {
      maxreducelenset = TRUE;
      arg++;
      if (arg >= argc)
        badusage();
      maxreducelen = atoi(argv[arg]);
    }
    else {
       if (argv[arg][0] == '-')
         badusage();
       if (strcmp(inf2,""))
         badusage();
       if (strcmp(gpname,"")==0)
         strcpy(gpname,argv[arg]);
       else {
         file_ip = TRUE;
         strcpy(inf2,argv[arg]);
	 strcpy(outf,inf2); strcat(outf,".reduced");
       }
    }
    arg++;
  }

  if (strlen(gpname)==0)
    badusage();

  strcpy(inf1,gpname);
  open = FALSE;
  if (rws_ip)
    strcat(inf1,".kbprog");
  else if (diff1_ip)
    strcat(inf1,".diff1");
  else if (diff2_ip)
    strcat(inf1,".diff2");
  else if (diff1c_ip)
    strcat(inf1,".diff1c");
  else {
    strcat(inf1,".kbprog");
    if (rfile = fopen(inf1,"r")) {
      rws_ip = TRUE;
      open = TRUE;
    }
    else {
      diff2_ip = TRUE;
      strcpy(inf1,gpname); strcat(inf1,".diff2");
    }
  }

  if (!open) 
    if ((rfile = fopen(inf1,"r")) == 0) {
      fprintf(stderr,"Cannot open file %s.\n",inf1);
      exit(1);
    }
 
  if (rws_ip) {
    read_kbinput_simple(rfile,FALSE);
    fclose(rfile);
    strcpy(inf1,gpname);
    strcat(inf1,".reduce");
    if ((rfile = fopen(inf1,"r")) == 0) {
      fprintf(stderr,"Cannot open file %s.\n",inf1);
      exit(1);
    }
    fsa_read(rfile,&reduction_fsa,DENSE,0,0,TRUE,fsaname);
    fclose(rfile);
    rws.reduction_machine = &reduction_fsa;
    num_gens = rws.num_gens;
    process_names(rws.gen_name,num_gens);
    tmalloc(rws.history,int,maxreducelen+1);
  }
  else {
    fsa_read(rfile,&wd_fsa,DENSE,0,0,TRUE,fsaname);
    fclose(rfile);
    num_gens = wd_fsa.alphabet->base->size;
    process_names(wd_fsa.alphabet->base->names,num_gens);
/* Set the pointers in the word-difference machine */
    fsa_table_dptr_init(&wd_fsa);
  }

/* Set the generator names and the reduction algorithm in the two cases. */
  names = rws_ip ? rws.gen_name : wd_fsa.alphabet->base->names;
  reduce_word = rws_ip ? rws_reduce : diff_reduce;

  tmalloc(gen_word,char,maxreducelen);
  if (file_ip) {
/* open file containing list of words, and deal with assignment */
    if ((rfile = fopen(inf2,"r")) == 0) {
      fprintf(stderr,"Cannot open file %s.\n",inf2);
      exit(1);
    }
    wfile = fopen(outf,"w");
    read_ident(rfile,buffer,&delim,FALSE);
    fprintf(wfile,"%s.reduced := [\n",buffer);
    if (delim != ':'){
      fprintf(stderr,
        "#Input error: file must contain a record assignment\n");
      exit(1);
    }
    check_next_char(rfile,'=');
    read_delim(rfile,&delim);
    if (delim != '['){
      fprintf(stderr,
        "#Input error: file must contain a record assignment\n");
      exit(1);
    }
    first = TRUE;
  }
  else {
    printf("#Input words for reduction.\n#Separate words with ','.\n");
    printf("#Terminate input with ';'.\n");
    rfile=stdin; wfile=stdout;
  }

  delim = 0;
  while (delim != ';' && delim != ']') {
    read_word(rfile,gen_word,gen_word+maxreducelen,&delim,
		names,num_gens,TRUE);
    reduce_word(gen_word);
    if (file_ip) {
      if (first) first = FALSE;
      else fprintf(wfile,",\n");
    }
    else
      fprintf(wfile,"reduces to:\n");
    strcpy(buffer,"  ");
    add_word_to_buffer(wfile,gen_word,names);
    fprintf(wfile,"%s",buffer);
    if (!file_ip) fprintf(wfile,"\n");
  }

  if (file_ip) {
    check_next_char(rfile,';');
    fprintf(wfile,"\n];\n");
    fclose(wfile);
  }

  if (rws_ip) {
    rws_clear(&rws);
    fsa_clear(&reduction_fsa);
  }
  else {
    fsa_clear(&wd_fsa);
  }
  tfree(gen_word);
  exit(0);
}
 
void
badusage()
{
    fprintf(stderr,
   "Usage: wordreduce [-kbprog/-diff1/-diff2/-diff1c] [-mrl maxreducelen]\n");
    fprintf(stderr,"\tgroupname [filename]\n");
    exit(1);
}
