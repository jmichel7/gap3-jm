/* rwsio2.c 18/1/95
 * This is the second file (first if rwsio.c) containing i/o functions for
 * re-writing systems. Some of the routines here are simplified versions
 * of those in rwsio.c, used when we are merely reading in a re-writing system
 * and not applying Knuth-Bendix to it.
 * These are used, for example, by the progrmas gpaxioms and reducewords.
 */

#include <stdio.h>

#include "defs.h"
#include "fsa.h"
#include "rws.h"
#include "externals.h"


extern int  	maxeqns,
         	maxreducelen;
extern boolean  maxreducelenset;
extern rewriting_system rws;

/* Functions defined in this file */
void read_kbinput_simple();
void read_gens();
void read_inverses();
void read_eqns_simple();
void rws_clear();

/* Functions used in this file defined in other files: */
void read_ident();
void read_delim();
void skip_gap_expression();
void read_int();
void read_string();
void process_names();
void read_words();

void
read_kbinput_simple(rfile,check)
        FILE *rfile;
        boolean check;
/* This is a simplified version of read_kbinput() in rwsio.c.
 * It is used when a rewriting system will be read in but not altered.
 * If check is true, then the words in the equations are checked for
 * validity - this could make input slower if there are many equations
 */
{ int delim, n, m, i, j;
  boolean isRWS=FALSE, seengens=FALSE, seeneqns=FALSE, maxeqnsset=FALSE;

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
    else if (strcmp(buffer,"isConfluent")==0) {
      read_ident(rfile,buffer,&delim,FALSE);
      rws.confluent = strcmp(buffer,"true")==0;
    }
    else if (strcmp(buffer,"maxeqns")==0){
      read_int(rfile,&n,&delim);
      if (seeneqns)
        fprintf(stderr,
            "Warning: maxeqns field should precede equations field.\n");
/* We will exclude ridiculously small values of limit parameters. */
      else if (n>16) {
         maxeqns = n;
	 maxeqnsset = TRUE;
      }
    }
    else if (strcmp(buffer,"maxreducelen")==0) {
      read_int(rfile,&n,&delim);
      if (n>4096 && !maxreducelenset)
        maxreducelen = n;
    }
    else if (strcmp(buffer,"maxstates")==0 || strcmp(buffer,"tidyint")==0 ||
             strcmp(buffer,"confnum")==0 || strcmp(buffer,"sorteqns")==0 ||
             strcmp(buffer,"maxoplen")==0 || strcmp(buffer,"maxstoredlen")==0 ||
	     strcmp(buffer,"maxwdiffs")==0 || strcmp(buffer,"silent")==0 ||
	     strcmp(buffer,"verbose")==0)
	skip_gap_expression(rfile,&delim);
    else if (strcmp(buffer,"generatorOrder")==0){
      read_gens(rfile);
      process_names(rws.gen_name,rws.num_gens);
      seengens = TRUE;
      read_delim(rfile,&delim);
    }
    else if (strcmp(buffer,"inverses")==0){
      if (!seengens){
        fprintf(stderr,
          "#Input error: generator field must precede inverses field\n");
        exit(1);
      }
      read_inverses(rfile);
      read_delim(rfile,&delim);
    }
    else if (strcmp(buffer,"ordering")==0)
      read_string(rfile,buffer,&delim);
    else if (strcmp(buffer,"weight")==0){
      if (!seengens){
        fprintf(stderr,
          "#Input error: generator field must precede weight field\n");
        exit(1);
      }
      check_next_char(rfile,'[');
      for (n=1;n<=rws.num_gens;n++) {
        read_int(rfile,&m,&delim);
        if (!(n<rws.num_gens && delim==',') && !(n==rws.num_gens&& delim==']')){
          fprintf(stderr,
            "#Input error:  format of weight field wrong\n");
          exit(1);
        }
      }
      read_delim(rfile,&delim);
    }
    else if (strcmp(buffer,"level")==0){
      if (!seengens){
        fprintf(stderr,
          "#Input error: generator field must precede level field\n");
        exit(1);
      }
      check_next_char(rfile,'[');
      for (n=1;n<=rws.num_gens;n++) {
        read_int(rfile,&m,&delim);
        if (!(n<rws.num_gens && delim==',') && !(n==rws.num_gens&& delim==']')){
          fprintf(stderr,
            "#Input error:  format of level field wrong\n");
          exit(1);
        }
      }
      read_delim(rfile,&delim);
    }
    else if (strcmp(buffer,"equations")==0){
      if (rws.num_gens!=0 && rws.inv_of==0){
        fprintf(stderr,
            "#Input error: record must have 'inverses' field\n");
        exit(1);
      }
      if (maxeqnsset) 
       /* increase it by rws.num_gens, to avoid a possible difficulty with
        * the inverse equations.
        */
         maxeqns += rws.num_gens;
      tmalloc(rws.eqns,reduction_equation,maxeqns+1);
      tmalloc(rws.testword1,char,maxreducelen+1);
      tmalloc(rws.testword2,char,maxreducelen+1);
      rws.num_eqns=0;
      read_eqns_simple(rfile,check);
      read_delim(rfile,&delim);
      seeneqns = TRUE;
    }
    else if (strcmp(buffer,"done")==0){
      if (!seeneqns){
        fprintf(stderr,
            "#Input error: 'equations' field must precede 'done' field\n");
        exit(1);
      }
      check_next_char(rfile,'[');
      read_delim(rfile,&delim);
      if (delim!='[' && delim!=']'){
        fprintf(stderr,"#Input error: '[' or ']' expected.\n");
        exit(1);
      }
      if (delim=='[') while (1){
        read_int(rfile,&i,&delim);
        if (delim=='.'){
          check_next_char(rfile,'.');
          read_int(rfile,&j,&delim);
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
      read_delim(rfile,&delim);
    }
    else {
      fprintf(stderr,
          "#WARNING: unknown record field: %s\n",buffer);
      skip_gap_expression(rfile,&delim);
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
}

void
read_gens(rfile)
        FILE * rfile;
/* Read the list of generator names into the array rws.gen_name, and
 * set rws.num_gens equal to the number of generators.
 */
{ int delim;
  rws.num_gens=0;
  check_next_char(rfile,'[');
  tmalloc(rws.gen_name,char *,MAXCHAR+1);
  do{
    read_ident(rfile,buffer,&delim,TRUE);
    if (delim!=',' && delim!=']'){
      fprintf(stderr,"#Input error: ',' or ']' expected.\n");
      exit(1);
    }
    else if (strlen(buffer)==0 && rws.num_gens==0 && delim==']')
      return; /* no generators */
    else{
      rws.num_gens++;
      if (rws.num_gens>MAXCHAR){
        fprintf(stderr,
            "Sorry -  no more than 127 monoid generators allowed.\n");
        exit(1);
      }
      tmalloc(rws.gen_name[rws.num_gens],char,strlen(buffer)+1);
      strcpy(rws.gen_name[rws.num_gens],buffer);
    }
  } while (delim != ']');
/* Add the padding-symbol to the end of the list - not always needed,
 * but no harm.
 */
  tmalloc(rws.gen_name[rws.num_gens+1],char,2);
  strcpy(rws.gen_name[rws.num_gens+1],"_");
}

void
read_inverses(rfile)
        FILE * rfile;
/* The list of inverses is read.
 * An inverse inv(g) of a generator g satisfies g.inv(g) = inv(g).g = 1 -
 * i.e. this lists only 2-sided inverses
 * so the inv function must be involutory.
 * If the monoid is a group, then of course all generators have inverses.
 * In general, the list may have gaps.
 */
{ int g, i, ino, delim;
  check_next_char(rfile,'[');
  tmalloc(rws.inv_of,int,MAXCHAR+2);
  for (i=1;i<=rws.num_gens;i++)
    rws.inv_of[i] = 0;
  g = 0;
  do{
    read_ident(rfile,buffer,&delim,TRUE);
    if (delim!=',' && delim!=']'){
      fprintf(stderr,"#Input error: ',' or ']' expected.\n");
      exit(1);
    }
    else if (rws.num_gens==0 && delim==']') /* no generators */
      return;
    else{
      g++;
      if (g>rws.num_gens){
          fprintf(stderr,"#Input error: inverse list too long.\n");
        exit(1);
      }
      if (strlen(buffer)==0)
        ino=0;
      else{
        if (algno == -1) {
          ino = 0;
          for (i=1;i<=rws.num_gens;i++) if (strcmp(buffer,rws.gen_name[i])==0) {
            ino = i;
            break;
          }
          if (ino==0) {
            fprintf(stderr,
                      "#Input error: unknown generator in inverses list.\n");
            exit(1);
          }
        }
        else {
          ino = (algno==0) ? gen_no[*buffer] : gen_no[atoi(buffer+algno)];
          if (ino==0 || strcmp(buffer,rws.gen_name[ino]) != 0){
            fprintf(stderr,
                "#Input error: unknown generator in inverses list.\n");
            exit(1);
          }
        }
        if (rws.inv_of[ino]!=0 && rws.inv_of[ino]!=g){
          fprintf(stderr,"#Input error: inverse list must be involutory.\n");
          exit(1);
        }
      }
      rws.inv_of[g]=ino;
    }
  } while (delim != ']');
  rws.inv_of[rws.num_gens+1] = rws.num_gens+1;
    /* In case there is a padding symbol */
}

void
read_eqns_simple(rfile,check)
        FILE *rfile;
        boolean check;
/* Read the initial reduction equations  -
 * simplified version of "read_eqns" in rwsio.c
 * Here we simply read them in and store them without any processing.
 */
{ int delim, i, ct, l1, l2, max;
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
    ct++;
    if (ct >= maxeqns){
      printf("#Too many equations - increase maxeqns to get started.\n");
      exit(1);
    }
    l1 = strlen(test1); l2 = strlen(test2);
    max = l1>l2 ? l1 : l2;
/* We assign the same space for both sides, since we may be balancing them
 * later.
 */
    tmalloc(rws.eqns[ct].lhs,char,max+1);
    strcpy(rws.eqns[ct].lhs,test1);
    tmalloc(rws.eqns[ct].rhs,char,max+1);
    strcpy(rws.eqns[ct].rhs,test2);

    read_delim(rfile,&delim);
    if (delim==']')
      break;
    if (delim!=','){
      fprintf(stderr,"#Input error: ',' expected.\n");
      exit(1);
    }
    check_next_char(rfile,'[');
  }

  rws.num_eqns += ct;
}

void
rws_clear(rwsptr)
        rewriting_system *rwsptr;
/* Free the assigned space associated with a rewriting-system,
 * excluding its reduction_automaton, which should be freed separately
 * with fsa_clear, if it is present.
 */
{ int i;
        for (i=1;i<=rwsptr->num_eqns;i++) {
                tfree((rwsptr->eqns)[i].lhs);
                tfree((rwsptr->eqns)[i].rhs);
        }
        tfree(rwsptr->eqns);
        for (i=1;i<=rwsptr->num_gens+1;i++)
                tfree(rwsptr->gen_name[i]);
        tfree(rwsptr->gen_name);
        tfree(rwsptr->weight);
        tfree(rwsptr->level);
        tfree(rwsptr->inv_of);
        tfree(rwsptr->history);
        tfree(rwsptr->slowhistory);
        tfree(rwsptr->slowhistorysp);
        tfree(rwsptr->preflen);
        tfree(rwsptr->prefno);
        tfree(rwsptr->testword1);
        tfree(rwsptr->testword2);
}
