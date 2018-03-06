/* file rwsio.c -  11. 10. 94.
 * 14.3.95. introduced read_extra_kbinput.
 * 12. 12. 94. implemented format type c) (see below).
 *
 * This file contains input and output routines to read and write GAP syntax
 * files for the Knuth-Bendix program.
 * (from 18/1/95) other such routines are in the file rwsio2.c.
 * The basic functions, which are shared by the finite state automata programs,
 * are in the file miscio.c.
 * Since words will be stored as ASCII strings, a maximum of 127
 * ( = MAXCHAR) monoid generators is allowed.
 * To help fast reading of large sets of words, two specials formats
   * for monoid generators are recognised:
 * a) single letter characters (usually 'A' will be inverse of 'a', etc.)
 * b) names of form <prefix><posint>, for a fixed prefix, where posint should
 *    be at most 255.
 * In all cases, rws.num_gens is the total number of monoid generators,
 * and the generator names are stored in the array rws.gen_name.
 * In case a), the variable algno is set equal to 0, and the array gen_no
 * is used to translate from rws.gen_name back to the gneerator number.
 * In case b), algno is set equal to the length of <prefix> (which must be
 * strictly positive), and gen_no is defined on the <posint> suffixes to
 * get the generator number back from the rws.gen_name.
 * If neither case a) nor case b) applies then algno is set equal to -1,
 * and names are located in the list by a straightforward linear search - of
 * course this will be considerably slower for long lists of words.
 */

#include "defs.h"
#include "fsa.h"
#include "rws.h"
#include "externals.h"

extern int 	tidyint,
		maxeqns,
		maxstates,
		current_maxstates,
		init_fsaspace,	
		maxwdiffs,
		confnum,
		maxslowhistoryspace,
		maxreducelen,
		num_states,
		**table;
extern boolean  tidyon,
                tidyintset,
                maxeqnsset,
                maxstatesset,
                maxreducelenset,
                orderingset,
                maxwdiffsset,
		confnumset,
		silentset,
		verboseset,
		resume_with_orig,
		outputprefixes;
extern fsa      reduction_fsa;
extern fsa      wd_fsa;
extern rewriting_system rws;

int 		*eqn_no;
		/* eqn_no[i]=j means the i-th equation read in becomes
		 * equation number j. (usually i=j) - needed only if the
		 * "done" values are read in.
		 */

/* Functions defined in this file: */
void initialise_reduction_fsa();
void initialise_eqns();
void read_eqns();
void read_done();
void read_kbinput();
void read_extra_kbinput();
void print_kboutput();
void print_wdoutput();

/* Functions defined in other files called from this file: */
boolean isdelim();
boolean isvalid();
void read_delim();
void skip_gap_expression();
void read_ident();
void read_int();
void read_string();
void read_gens();
void read_inverses();
void process_names();
void build_quicktable();
int modify_table();
int insert();
void printbuffer();
void add_to_buffer();
int add_word_to_buffer();
int int_len();
void fsa_init();
void fsa_table_init();
void srec_copy();
void fsa_print();

void
initialise_reduction_fsa()
/* We initialise the reduction automaton as an fsa type.
 * Some of its components like fsa.states->size are stored elsewhere
 * (as num_states) and will only be updated at the end, before printing.
 */
{ int i;
  fsa_init(&reduction_fsa);

  reduction_fsa.states->type = SIMPLE;
  reduction_fsa.states->size = 0;

  reduction_fsa.alphabet->type = IDENTIFIERS;
  reduction_fsa.alphabet->size = rws.num_gens;
  tmalloc(reduction_fsa.alphabet->names,char *,rws.num_gens+1);
  for (i=1;i<=rws.num_gens;i++) {
    tmalloc(reduction_fsa.alphabet->names[i],char,strlen(rws.gen_name[i])+1);
    strcpy(reduction_fsa.alphabet->names[i],rws.gen_name[i]);
  }
  
  reduction_fsa.num_initial = 1;
  tmalloc(reduction_fsa.initial,int,2);
  reduction_fsa.initial[1] = 1;
/* All states will be accepting, so reduction_fsa.num_accepting will
 * be equal to states.size.
 * This is manipulated as num_states in the program, so we do not set it until
 * the end.
 */
  reduction_fsa.flags[DFA] = TRUE;
  reduction_fsa.flags[TRIM] = TRUE;
  reduction_fsa.flags[RWS] = TRUE;
  current_maxstates = maxstates==0 ? rws.num_gens==0 ? 1 :
			init_fsaspace/rws.num_gens : maxstates;
  fsa_table_init(reduction_fsa.table,current_maxstates,rws.num_gens);
  reduction_fsa.table->printing_format = DENSE;
  table = reduction_fsa.table->table_data_ptr; /* for fast access to table */
  tmalloc(rws.preflen,int,current_maxstates+1);
  tmalloc(rws.prefno,int,current_maxstates+1);
  build_quicktable();
}

void
initialise_eqns()
/* This function is called after the generators and inverses have been read
 * in. We initialise the list of equations by making one equation
 * with lhs a*inv(a) and rhs identity for each generator a.
 * (These equations will probably be in the list that we read in as well,
 *  but that does not matter.)
 * We also initialise the finite state automaton that is used for reducing
 * words.
 */
{  int i;
/* maxeqns, maxstates and maxreducelen should have been set by now */
  tmalloc(rws.eqns,reduction_equation,maxeqns+1);
  tmalloc(eqn_no,int,maxeqns+1);
  tmalloc(rws.testword1,char,maxreducelen+1);
  tmalloc(rws.testword2,char,maxreducelen+1);
  rws.num_eqns = 0;
  for (i = 1; i <= rws.num_gens; i++) if (rws.inv_of[i]) {
     rws.num_eqns++;
     tmalloc(rws.eqns[rws.num_eqns].lhs, char, 3);
     tmalloc(rws.eqns[rws.num_eqns].rhs, char, 1);
     rws.eqns[rws.num_eqns].lhs[0] = i;
     rws.eqns[rws.num_eqns].lhs[1] = rws.inv_of[i];
     rws.eqns[rws.num_eqns].lhs[2] = '\0';
     rws.eqns[rws.num_eqns].rhs[0] = '\0';
     rws.eqns[rws.num_eqns].done = TRUE;
  }
  rws.num_inveqns = rws.num_eqns;
  rws.reduction_machine = &reduction_fsa;
  initialise_reduction_fsa();
}

void
read_eqns(rfile,check)
	FILE *rfile;
	boolean check;
/* Read the initial reduction equations and install them. */
{ int delim, i, ct, iv;
  char *test1 = rws.testword1, *test2 = rws.testword2;
  check_next_char(rfile,'[');
  read_delim(rfile,&delim);
  if (delim==']')
    return;
  if (delim!='['){
    fprintf(stderr,"#Input error: '[' expected.\n");
    exit(1);
  }
  ct = 0;
  while (1){
    ct++;
    read_word(rfile,test1,test1+maxreducelen,
                                        &delim,rws.gen_name,rws.num_gens,check);
    if (delim!=','){
      fprintf(stderr,"#Input error: ',' expected.\n");
      exit(1);
    }
    read_word(rfile,test2,test2+maxreducelen,
					&delim,rws.gen_name,rws.num_gens,check);
    if (delim!=']'){
      fprintf(stderr,"#Input error: ']' expected.\n");
      exit(1);
    }
    if (strlen(test1)!=2 || strlen(test2)!=0 || test1[1]!=rws.inv_of[test1[0]]){
/* this is  NOT an inverse equation, which we already know about. */
      if ((iv=insert( &(rws.eqns[rws.num_eqns+1].lhs),
                  &(rws.eqns[rws.num_eqns+1].rhs) )) > 0 ){
          i = modify_table(rws.num_eqns+1);
          if (i == -1){
            fprintf(stderr,"#maxstates is too small. Cannot get started.\n");
            exit(1);
          }
          if (i == 1){
            rws.num_eqns++;
            if (rws.num_eqns >= maxeqns){
             printf("#Too many equations - increase maxeqns to get started.\n");
              exit(1);
            }
            rws.eqns[rws.num_eqns].done = FALSE;
	    eqn_no[ct] = rws.num_eqns;
          }
      }
      else if (iv== -1) {
        fprintf(stderr,
             "#Error: input equation has lhs or rhs longer than limit set.\n");
        exit(1);
      }
    }

    read_delim(rfile,&delim);
    if (delim==']')
      break;
    if (delim!=','){
      fprintf(stderr,"#Input error: ',' expected.\n");
      exit(1);
    }
    check_next_char(rfile,'[');
  }
}

void
read_done(rfile)
	FILE *rfile;
/* Read the list of equation numbers that have already been processed. */
{ int delim, i, j, n;
  check_next_char(rfile,'[');
  read_delim(rfile,&delim);
  if (delim==']')
    return;
  if (delim!='['){
    fprintf(stderr,"#Input error: '[' expected.\n");
    exit(1);
  }
  while (1){
    read_int(rfile,&i,&delim);
    if (delim=='.'){
      check_next_char(rfile,'.');
      read_int(rfile,&j,&delim);
      for (n=i;n<=j;n++) if (n>0 && n<=maxeqns && eqn_no[n]>0)
        rws.eqns[eqn_no[n]].done = TRUE;
    }
    else if (i>0 && i<=maxeqns && eqn_no[i]>0)
      rws.eqns[eqn_no[i]].done = TRUE;

    read_delim(rfile,&delim);
    if (delim==']')
      break;
    if (delim!=','){
      fprintf(stderr,"#Input error: ',' expected.\n");
      exit(1);
    }
    check_next_char(rfile,'[');
  }

}

void
read_kbinput(rfile,check)
	FILE *rfile;
	boolean check;
/* This function reads the full input for the Knuth-Bendix program 
 * from the file rfile, which should already be open..
 * The rewriting system is read into the externally defined rws.
 * If check is true, then the words in the equations are checked for
 * validity - this could make input slower if there are many equations
 */
{ int delim, n, m;
  boolean isRWS=FALSE, seengens=FALSE, seeneqns=FALSE;

  read_ident(rfile,rws.name,&delim,FALSE);
  if (delim != ':'){
    fprintf(stderr,
        "#Input error: file must contain a record assignment\n");
    exit(1);
  }
  check_next_char(rfile,'=');
  read_ident(rfile,buffer,&delim,FALSE);
  if (delim != '(' || strcmp(buffer,"rec")!=0){
    fprintf(stderr,
        "#Input error: file must contain a record assignment\n");
    exit(1);
  }

/* main loop reading the fields of the record follows. */
  do{
    read_ident(rfile,buffer,&delim,FALSE);
    if (delim != ':'){
      fprintf(stderr,
          "#Input error: bad record field assignment\n");
      exit(1);
    }
    check_next_char(rfile,'=');
    if (strcmp(buffer,"isRWS")==0){
      isRWS = TRUE;
      read_ident(rfile,buffer,&delim,FALSE);
      if (strcmp(buffer,"true")!=0){
        fprintf(stderr,
            "#Input error: isRWS field must equal \"true\"\n");
        exit(1);
      }
    }
    else if (strcmp(buffer,"isConfluent")==0){
      read_ident(rfile,buffer,&delim,FALSE);
      if (strcmp(buffer,"true")==0 && !resume_with_orig){
        fprintf(stderr,
            "#System is already confluent!\n");
        exit(0);
      }
    }
    else if (strcmp(buffer,"tidyint")==0){
      read_int(rfile,&n,&delim);
/* Parameters on the command-line override those in the input file */
      if (!tidyintset && n>0)
         tidyint = n;
    }
    else if (strcmp(buffer,"maxeqns")==0){
      read_int(rfile,&n,&delim);
      if (seeneqns)
        fprintf(stderr,
            "Warning: maxeqns field should precede equations field.\n");
/* We will exclude ridiculously small values of limit parameters. */
      else if (!maxeqnsset && !resume_with_orig && n>16)
         maxeqns = n;
    }
    else if (strcmp(buffer,"maxstates")==0){
      read_int(rfile,&n,&delim);
      if (seeneqns)
        fprintf(stderr,
            "Warning: maxstates field should precede equations field.\n");
      else if (!maxstatesset && n>128)
         maxstates = n;
    }
    else if (strcmp(buffer,"confnum")==0){
      read_int(rfile,&n,&delim);
      if (!confnumset && n>0)
         confnum = n;
    }
    else if (strcmp(buffer,"sorteqns")==0){
      read_ident(rfile,buffer,&delim,FALSE);
      if (strcmp(buffer,"true")==0)
         rws.sorteqns = TRUE;
    }
    else if (strcmp(buffer,"maxoplen")==0){
      read_int(rfile,&n,&delim);
      if (rws.maxoplen==0 && n>=0){
         rws.sorteqns = TRUE;
         rws.maxoplen = n;
      }
    }
    else if (strcmp(buffer,"maxstoredlen")==0){
      check_next_char(rfile,'[');
      read_int(rfile,&n,&delim);
      if (delim != ','){
        fprintf(stderr,
          "#Input error:  format of maxstoredlen field wrong\n");
        exit(1);
      }
      read_int(rfile,&m,&delim);
      if (delim != ']'){
        fprintf(stderr,
          "#Input error:  format of maxstoredlen field wrong\n");
        exit(1);
      }
      if (n>0 && m>0 && rws.maxlenleft==0 && rws.maxlenright==0){
	rws.maxlenleft = n;
	rws.maxlenright = m;
      }
      read_delim(rfile,&delim);
    }
    else if (strcmp(buffer,"maxoverlaplen")==0){
      read_int(rfile,&n,&delim);
      if (n>0)
         rws.maxoverlaplen = n;
    }
    else if (strcmp(buffer,"maxwdiffs")==0){
      read_int(rfile,&n,&delim);
      if (!maxwdiffsset && n>16)
         maxwdiffs = n;
    }
    else if (strcmp(buffer,"maxreducelen")==0){
      read_int(rfile,&n,&delim);
      if (n>4096 && !maxreducelenset)
        maxreducelen = n;
    }
    else if (strcmp(buffer,"silent")==0){
      read_ident(rfile,buffer,&delim,FALSE);
      if (!silentset && !verboseset && strcmp(buffer,"true")==0)
         print_level=0;
    }
    else if (strcmp(buffer,"verbose")==0){
      read_ident(rfile,buffer,&delim,FALSE);
      if (!silentset && !verboseset && strcmp(buffer,"true")==0)
         print_level=2;
    }
    else if (strcmp(buffer,"veryVerbose")==0){
      read_ident(rfile,buffer,&delim,FALSE);
      if (!silentset && !verboseset && strcmp(buffer,"true")==0)
         print_level=3;
    }
    else if (strcmp(buffer,"RabinKarp")==0){
      check_next_char(rfile,'[');
      read_int(rfile,&n,&delim);
      if (delim != ','){
        fprintf(stderr,
          "#Input error:  format of maxstoredlen field wrong\n");
        exit(1);
      }
      read_int(rfile,&m,&delim);
      if (delim != ']'){
        fprintf(stderr,
          "#Input error:  format of maxstoredlen field wrong\n");
        exit(1);
      }
      if (n>0 && m>0 && rws.rkminlen==0 && rws.rkmineqns==0){
	rws.rkminlen = n;
	rws.rkmineqns = m;
      }
      read_delim(rfile,&delim);
    }
    else if (strcmp(buffer,"ordering")==0){
      read_string(rfile,buffer,&delim);
      if (!orderingset){
          if (strcmp(buffer,"shortlex")==0)
            rws.ordering = SHORTLEX;
          else if  (strcmp(buffer,"recursive")==0)
            rws.ordering = RECURSIVE;
          else if (strcmp(buffer,"rt_recursive")==0)
            rws.ordering = RT_RECURSIVE;
          else if (strcmp(buffer,"wtlex")==0)
            rws.ordering = WTLEX;
          else if (strcmp(buffer,"wreathprod")==0)
            rws.ordering = WREATHPROD;
          else{
            fprintf(stderr,
                "#Input error: invalid string for ordering field\n");
            exit(1);
          }
      }
    }
    else if (strcmp(buffer,"generatorOrder")==0){
      read_gens(rfile);
      process_names(rws.gen_name,rws.num_gens);
      seengens = TRUE;
      read_delim(rfile,&delim);
    }
    else if (strcmp(buffer,"weight")==0){
      if (rws.ordering != WTLEX)
        skip_gap_expression(rfile,&delim);
      else {
        if (!seengens){
          fprintf(stderr,
            "#Input error: generator field must precede weight field\n");
          exit(1);
        }
        tmalloc(rws.weight,int,rws.num_gens+1);
        check_next_char(rfile,'[');
        for (n=1;n<=rws.num_gens;n++) {
            read_int(rfile,rws.weight+n,&delim);
            if ((n<rws.num_gens && delim!=',') ||
                 (n==rws.num_gens && delim!=']')){
              fprintf(stderr,"#Input error: ',' or ']' expected.\n");
              exit(1);
            }
        }
        read_delim(rfile,&delim);
      }
    }
    else if (strcmp(buffer,"level")==0){
      if (rws.ordering != WREATHPROD)
        skip_gap_expression(rfile,&delim);
      else {
        if (!seengens){
          fprintf(stderr,
            "#Input error: generator field must precede level field\n");
          exit(1);
        }
        tmalloc(rws.level,int,rws.num_gens+1);
        check_next_char(rfile,'[');
        for (n=1;n<=rws.num_gens;n++) {
            read_int(rfile,rws.level+n,&delim);
            if ((n<rws.num_gens && delim!=',') ||
                 (n==rws.num_gens && delim!=']')){
              fprintf(stderr,"#Input error: ',' or ']' expected.\n");
              exit(1);
            }
        }
        read_delim(rfile,&delim);
      }
    }
    else if (strcmp(buffer,"inverses")==0){
      if (!seengens){
        fprintf(stderr,
          "#Input error: generator field must precede inverses field\n");
        exit(1);
      }
      read_inverses(rfile);
      initialise_eqns();
      read_delim(rfile,&delim);
    }
    else if (strcmp(buffer,"equations")==0){
      if (rws.ordering==WTLEX && rws.weight==0) {
         fprintf(stderr,
                 "Input error: weight field missing (for this ordering).\n");
         exit(1);
      }
      if (rws.ordering==WREATHPROD && rws.level==0) {
         fprintf(stderr,
                 "Input error: level field missing (for this ordering).\n");
         exit(1);
      }
      if (rws.num_gens!=0 && rws.inv_of==0){
        fprintf(stderr,
            "#Input error: record must have 'inverses' field\n");
        exit(1);
      }
      read_eqns(rfile,check);
      read_delim(rfile,&delim);
      seeneqns = TRUE;
    }
    else if (strcmp(buffer,"done")==0){
      if (!seeneqns){
        fprintf(stderr,
            "#Input error: 'equations' field must precede 'done' field\n");
        exit(1);
      }
      read_done(rfile);
      read_delim(rfile,&delim);
    }
    else {
      fprintf(stderr,
          "#Input error: unknown record field: %s\n",buffer);
      exit(1);
    }
    if (delim != ')' && delim != ','){
      fprintf(stderr,
          "#Input error:  field %s assignment must end ',' or ')', not %c\n",
               buffer,delim);
      exit(1);
    }
  } while (delim != ')');

  check_next_char(rfile,';');
  if (!isRWS){
    fprintf(stderr,
        "#Input error: record must have 'isRWS' field\n");
    exit(1);
  }
  if (rws.num_gens!=0 && rws.inv_of==0){
    fprintf(stderr,
        "#Input error: record must have 'inverses' field\n");
    exit(1);
  }
  tfree(eqn_no);
}

void
read_extra_kbinput(rfile,check)
        FILE *rfile;
        boolean check;
/* This function reads the additional equations (from the original
 * input file), when these are to be re-adjoined to the output equations
 * (which is what happens under the -ro option of kbprog).
 * The file rfile should already be open.
 * If check is true, then the words in the equations are checked for
 * validity
 */
{ int delim;

  read_ident(rfile,rws.name,&delim,FALSE);
  if (delim != ':'){
    fprintf(stderr,
        "#Input error: file must contain a record assignment\n");
    exit(1);
  }
  check_next_char(rfile,'=');
  read_ident(rfile,buffer,&delim,FALSE);
  if (delim != '(' || strcmp(buffer,"rec")!=0){
    fprintf(stderr,
        "#Input error: file must contain a record assignment\n");
    exit(1);
  }

/* main loop reading the fields of the record follows. */
  do {
    read_ident(rfile,buffer,&delim,FALSE);
    if (delim != ':'){
      fprintf(stderr,
          "#Input error: bad record field assignment\n");
      exit(1);
    }
    check_next_char(rfile,'=');
    if (strcmp(buffer,"equations")==0){
      tmalloc(eqn_no,int,maxeqns+1);
      read_eqns(rfile,check);
      read_delim(rfile,&delim);
    }
    else
      skip_gap_expression(rfile,&delim);
  } while (delim != ')');
  check_next_char(rfile,';');
  tfree(eqn_no);
}


void
print_kboutput(wfile)
	FILE *wfile;
/* This function prints the output from the KB program to the file named
 * wfile, which should already be open for writing.
 * Note that the rewriting system rws and its reduction-fsa are
 * defined externally.
 * This starts with the rewriting system, in the same format as the
 * input, and is followed by a new field, the reduction automaton.
 * The fsa function fsa_print is used to print this.
 */
{ int i, j, i1, l, n;
  boolean in,first;
  char **pref, *w;
  time_t t;

  t = time( (time_t *)0 );
  sprintf(buffer,"%s := rec(",rws.name);
  printbuffer(wfile);
  add_to_buffer(16,"isRWS");
  sprintf(buffer+strlen(buffer)," := true,");
  printbuffer(wfile);

  add_to_buffer(16,"isConfluent");
  if (rws.confluent)
    sprintf(buffer+strlen(buffer)," := true,");
  else
    sprintf(buffer+strlen(buffer)," := false,");
  printbuffer(wfile);

  if (rws.confluent) {
   /* Since the number of equations won't be increased again, we specify this
    * number here as the maxeqns component. This is useful for programs that
    * re-read the system later.
    */
    add_to_buffer(16,"maxeqns");
    sprintf(buffer+strlen(buffer)," := %d,",rws.num_eqns);
    printbuffer(wfile);
  }

  add_to_buffer(16,"generatorOrder");
  sprintf(buffer+strlen(buffer)," := [");
  for (i=1;i<=rws.num_gens;i++){
    if (i==1 || strlen(buffer)+strlen(rws.gen_name[i]) <= 76){
       if (i>1)
          add_to_buffer(0,",");
       sprintf(buffer+strlen(buffer),"%s",rws.gen_name[i]);
    }
    else {
      add_to_buffer(0,",");
      printbuffer(wfile);
      add_to_buffer(21,"");
      sprintf(buffer+strlen(buffer),"%s",rws.gen_name[i]);
    }
  }
  add_to_buffer(0,"],");
  printbuffer(wfile);

  add_to_buffer(16,"ordering");
  if (rws.ordering==SHORTLEX)
    sprintf(buffer+strlen(buffer)," := \"shortlex\",");
  else if (rws.ordering==RECURSIVE)
    sprintf(buffer+strlen(buffer)," := \"recursive\",");
  else if (rws.ordering==RT_RECURSIVE)
    sprintf(buffer+strlen(buffer)," := \"rt_recursive\",");
  else if (rws.ordering==WTLEX)
    sprintf(buffer+strlen(buffer)," := \"wtlex\",");
  else if (rws.ordering==WREATHPROD)
    sprintf(buffer+strlen(buffer)," := \"wreathprod\",");
  printbuffer(wfile);

  if (rws.ordering==WTLEX) {
    add_to_buffer(16,"weight");
    sprintf(buffer+strlen(buffer)," := [");
    for (i=1;i<=rws.num_gens;i++){
      if (i>1) add_to_buffer(0,",");
      sprintf(buffer+strlen(buffer),"%d",rws.weight[i]);
    }
    add_to_buffer(0,"],");
    printbuffer(wfile);
  }

  if (rws.ordering==WREATHPROD) {
    add_to_buffer(16,"level");
    sprintf(buffer+strlen(buffer)," := [");
    for (i=1;i<=rws.num_gens;i++){
      if (i>1) add_to_buffer(0,",");
      sprintf(buffer+strlen(buffer),"%d",rws.level[i]);
    }
    add_to_buffer(0,"],");
    printbuffer(wfile);
  }

  add_to_buffer(16,"inverses");
  sprintf(buffer+strlen(buffer)," := [");
  for (i=1;i<=rws.num_gens;i++){
    if (i>1)
      add_to_buffer(0,",");
    if (rws.inv_of[i] != 0){
      if (strlen(buffer)+strlen(rws.gen_name[rws.inv_of[i]]) > 76){
        printbuffer(wfile);
        add_to_buffer(21,"");
      }
      sprintf(buffer+strlen(buffer),"%s",rws.gen_name[rws.inv_of[i]]);
    }
  }
  add_to_buffer(0,"],");
  printbuffer(wfile);

  add_to_buffer(16,"equations");
  sprintf(buffer+strlen(buffer)," := [");
  for (i=1;i<=rws.num_eqns;i++){
    printbuffer(wfile);
    add_to_buffer(10,"[");
    n=add_word_to_buffer(wfile,rws.eqns[i].lhs,rws.gen_name);
    sprintf(buffer+strlen(buffer),",");
    if (n>0 || strlen(buffer)>40){
      printbuffer(wfile);
      add_to_buffer(12,"");
    }
    add_word_to_buffer(wfile,rws.eqns[i].rhs,rws.gen_name);
    if (i==rws.num_eqns)
      sprintf(buffer+strlen(buffer),"]");
    else
      sprintf(buffer+strlen(buffer),"],");
  }
  printbuffer(wfile);
  if (rws.confluent)
    add_to_buffer(8,"]");
  else
    add_to_buffer(9,"],");
  printbuffer(wfile);

  if (!rws.confluent){
/* print the list of equations that have been processed, in case the run is
 * continued later.
 */
    add_to_buffer(16,"done");
    sprintf(buffer+strlen(buffer)," := [");
    first = TRUE;
    in = FALSE;
    for (i=1;i<=rws.num_eqns+1;i++){
      if (i<=rws.num_eqns && !in && rws.eqns[i].done){
        in = TRUE; i1 = i;
      }
      else if (in && (i>rws.num_eqns || !rws.eqns[i].done)){
        in = FALSE;
        if (!first)
          sprintf(buffer+strlen(buffer),",");
        printbuffer(wfile);
        first = FALSE;
        add_to_buffer(10,"[");
        if (i==i1+1)
          sprintf(buffer+strlen(buffer),"%d]",i1);
        else
          sprintf(buffer+strlen(buffer),"%d..%d]",i1,i-1);
      }
    }
    printbuffer(wfile);
    add_to_buffer(8,"]");
    printbuffer(wfile);
  }

  sprintf(buffer,");");
  printbuffer(wfile);

  strcat(rws.name,".reductionFSA");
/* We complete the definition of the fsa. before printing it */
  reduction_fsa.states->size = num_states;
  reduction_fsa.num_accepting = num_states;
  if (num_states==1) {
    tmalloc(reduction_fsa.accepting,int,2);
    reduction_fsa.accepting[1]=1;
  }
  if (outputprefixes){
/* In this case we want to output the prefixes of the left-hand sides of
 * equations associated with the states. First we must calculate them.
 */
    reduction_fsa.states->type = WORDS;
    reduction_fsa.states->alphabet_size = reduction_fsa.alphabet->size;
    for (i=1;i<=reduction_fsa.alphabet->size;i++) {
     tmalloc(reduction_fsa.states->alphabet[i],char,
                                   strlen(reduction_fsa.alphabet->names[i])+1);
     strcpy(reduction_fsa.states->alphabet[i],reduction_fsa.alphabet->names[i]);
    }
    tmalloc(reduction_fsa.states->names,char *,num_states+1);
    pref = reduction_fsa.states->names;
    for (i=1;i<=num_states;i++){
      l = rws.preflen[i];
      w = rws.eqns[rws.prefno[i]].lhs;
/* Copy prefix from w to pref[i] */
      tmalloc(pref[i],char,l+1);
      for (j=0;j<l;j++)
         pref[i][j] = w[j];
      pref[i][l] = '\0';
    }
  }
}

void
print_wdoutput(wfile,suffix)
	FILE *wfile;
	char *suffix;
/* This function prints the word-difference from the KB program to the file
 * named wfile, which should already be open for writing.
 * The fsa function fsa_print is used to print this.
 */
{
  time_t t;
  char diffname[128];
	
  strcpy(diffname,rws.name);
  strcat(diffname,suffix);
  fsa_print(wfile,&wd_fsa,diffname);
}
