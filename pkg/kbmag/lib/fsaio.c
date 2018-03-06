/* file fsaio.c - 11. 10. 94.
 * 5.12.94. - changes for sparse storage type with some dense rows.
 * 17.11.94. Added parts to do with labeled srec type.
 * This file contains io-routines for manipulating finite state automata
 * They can deal with both deterministic and non-deterministic automata,
 * but the functions in fsa.c currently only deal with dfa's.
 */

#define MAXWORDLEN	65535

#include <stdlib.h>
#include "defs.h"
#include "fsa.h"
#include "externals.h"

char 		*wbuffer; /* used only for calls of read_word - this can
                           * be allocated temporarily.
                           */

/* Functions defined in this file: */
void srec_print();
void table_print();
void fsa_print();
void srec_read();
void table_read();
void fsa_read();
void compressed_transitions_read();


/* Functions used in this file defined outside this file: */
void printbuffer();
void add_to_buffer();
int add_word_to_buffer();
void read_ident();
void read_string();
void read_int();
void read_word();
void read_delim();
int int_len();
void fsa_init();

void
srec_print(wfile,srptr,name,offset,endstring)
/* Print the set record *srptr. Follows corresponding GAP routine.
 * Currently, rather arbitrarily, identifiers and strings names are
 * printed in dense format, and words in sparse format.
 */
	FILE *wfile;
	srec *srptr;
	char *name;
	int  offset;
	char *endstring;
	
{ int ct;
  srec sr;
  boolean first;

  sr = *srptr;
  *buffer = '\0';
  if (strlen(name)==0)
     add_to_buffer(0,"rec(");
  else{
     add_to_buffer(12+offset,name);
     add_to_buffer(0," := rec(");
  }

  printbuffer(wfile);
  add_to_buffer(16+offset,"type");
  sprintf(buffer+strlen(buffer)," := \"%s\",",type_names[sr.type]);
  printbuffer(wfile);

  add_to_buffer(16+offset,"size");
  sprintf(buffer+strlen(buffer)," := %d",sr.size);
  if (sr.type!=SIMPLE)
    add_to_buffer(0,",");
  printbuffer(wfile);

  if (sr.type == PRODUCT){
    add_to_buffer(16+offset,"arity");
    sprintf(buffer+strlen(buffer)," := %d,",sr.arity);
    printbuffer(wfile);

    add_to_buffer(16+offset,"padding");
    sprintf(buffer+strlen(buffer)," := %c,",sr.padding);
    printbuffer(wfile);

    srec_print(wfile,sr.base,"base",offset+4," ");
  }

  else if (sr.type==IDENTIFIERS){
    add_to_buffer(16+offset,"format");
    sprintf(buffer+strlen(buffer)," := \"dense\",");
    printbuffer(wfile);

    add_to_buffer(16+offset,"names");
    sprintf(buffer+strlen(buffer)," := [");
    ct = 1;
    while (ct <= sr.size){
      if (ct==1 || strlen(buffer)+strlen(sr.names[ct]) <= 76){
         if (ct>1)
            add_to_buffer(0,",");
          sprintf(buffer+strlen(buffer),"%s",sr.names[ct]);
      }
      else {
        add_to_buffer(0,",");
        printbuffer(wfile);
        add_to_buffer(21+offset,"");
        sprintf(buffer+strlen(buffer),"%s",sr.names[ct]);
      }
      ct++;
    }
    add_to_buffer(0,"]");
    printbuffer(wfile);
  }

  else if (sr.type==STRINGS){
    add_to_buffer(16+offset,"format");
    sprintf(buffer+strlen(buffer)," := \"dense\",");
    printbuffer(wfile);

    add_to_buffer(16+offset,"names");
    sprintf(buffer+strlen(buffer)," := [");
    ct = 1;
    while (ct <= sr.size){
      if (ct==1 || strlen(buffer)+strlen(sr.names[ct])+2 <= 76){
         if (ct>1)
            add_to_buffer(0,",");
          sprintf(buffer+strlen(buffer),"\"%s\"",sr.names[ct]);
      }
      else {
        add_to_buffer(0,",");
        printbuffer(wfile);
        add_to_buffer(21+offset,"");
        sprintf(buffer+strlen(buffer),"\"%s\"",sr.names[ct]);
      }
      ct++;
    }
    add_to_buffer(0,"]");
    printbuffer(wfile);
  }

  else if (sr.type==WORDS){
    add_to_buffer(16+offset,"alphabet");
    sprintf(buffer+strlen(buffer)," := [");
    ct = 1;
    while (ct <= sr.alphabet_size){
      if (ct==1 || strlen(buffer)+strlen(sr.alphabet[ct]) <= 76){
         if (ct>1)
            add_to_buffer(0,",");
          sprintf(buffer+strlen(buffer),"%s",sr.alphabet[ct]);
      }
      else {
        add_to_buffer(0,",");
        printbuffer(wfile);
        add_to_buffer(21+offset,"");
        sprintf(buffer+strlen(buffer),"%s",sr.alphabet[ct]);
      }
      ct++;
    }
    add_to_buffer(0,"],");
    printbuffer(wfile);

    add_to_buffer(16+offset,"format");
    sprintf(buffer+strlen(buffer)," := \"sparse\",");
    printbuffer(wfile);

    add_to_buffer(16+offset,"names");
    sprintf(buffer+strlen(buffer)," := [");
    printbuffer(wfile);
    for (ct=1;ct <= sr.size; ct++){ 
/* If the set is large, we will format less nicely, to save space. */
      if (sr.size > 128)
        add_to_buffer(0,"\t");
      else
        add_to_buffer(12+offset,"");
      sprintf(buffer+strlen(buffer),"[%d,",ct);
      add_word_to_buffer(wfile,sr.names[ct],sr.alphabet);
      if (ct==sr.size)
        add_to_buffer(0,"]");
      else
        add_to_buffer(0,"],");
      printbuffer(wfile);
    }
    add_to_buffer(11+offset,"]");
    printbuffer(wfile);
  }

  else if (sr.type==LABELED){
    srec_print(wfile,sr.labels,"labels",offset+4,",");

    add_to_buffer(16+offset,"format");
    sprintf(buffer+strlen(buffer)," := \"sparse\",");
    printbuffer(wfile);

    add_to_buffer(16+offset,"setToLabels");
    sprintf(buffer+strlen(buffer)," := [");
    printbuffer(wfile);
    first = TRUE;
    for (ct=1;ct<=sr.size;ct++) if (sr.setToLabels[ct] != 0) {
      if (first)
        first = FALSE;
      else {
        strcat(buffer,",");
        printbuffer(wfile);
      }
      if (sr.size <= 128)
        add_to_buffer(21+offset,"");
      else
        add_to_buffer(0,"\t");
      sprintf(buffer+strlen(buffer),"[%d,%d]",ct,sr.setToLabels[ct]);
    }
    printbuffer(wfile);
    add_to_buffer(21+offset,"]");
    printbuffer(wfile);
  }

  if (strlen(name)==0)
     add_to_buffer(0,")");
  else
     add_to_buffer(12+offset,")");
  add_to_buffer(0,endstring);
  printbuffer(wfile);
}

void
table_print(wfile,tableptr,name,offset,endstring,ns,ne)
/* Print the table record *tableptr. */
	FILE *wfile;
	table_struc  *tableptr;
	char *name;
	int  offset;
	char *endstring;
	int ns, ne;
	
{ int **table, ct, g, i, j, k, nl, dr, *ptr, *ptre;
  boolean dense, densepf, first, firstg;
  
  table = tableptr->table_data_ptr;
  dense = tableptr->table_type == DENSE;
  densepf = tableptr->printing_format == DENSE;
  dr = tableptr->denserows;

  *buffer = '\0';
  if (strlen(name)==0)
     add_to_buffer(0,"rec(");
  else{
     add_to_buffer(12+offset,name);
     add_to_buffer(0," := rec(");
  }

  printbuffer(wfile);

/* Calculate number of nonzero transitions - this ought to be stored in
 * tableptr->numTransitions, but we won't rely on that (unless the table
 * is stored externally).
 */

  add_to_buffer(offset+16,"format");
  if (densepf)
    sprintf(buffer+strlen(buffer)," := \"dense deterministic\",");
  else
    sprintf(buffer+strlen(buffer)," := \"sparse\",");
  printbuffer(wfile);

  if (tableptr->filename == 0) {
    if (dense) {
      ct = 0;
      for (g=1;g<=ne;g++) for (i=1;i<=ns;i++)
        if (target(dense,table,g,i,dr) != 0)
          ct++;
    }
    else {
      ct = 0;
      for (g=1;g<=ne;g++) for (i=1;i<=dr;i++)
        if (target(dense,table,g,i,dr) != 0)
          ct++;
      if (dr < ns)
        ct += (table[ns+1]-table[dr+1])/2;
    }
    tableptr->numTransitions = ct;
  }

  add_to_buffer(offset+16,"numTransitions");
  sprintf(buffer+strlen(buffer)," := %d,",tableptr->numTransitions);
  printbuffer(wfile);

  if (tableptr->filename) {
    add_to_buffer(offset+16,"filename");
    sprintf(buffer+strlen(buffer)," := \"%s\"",tableptr->filename);
    printbuffer(wfile);
    if (strlen(name)==0)
       add_to_buffer(0,")");
    else
       add_to_buffer(12+offset,")");
    add_to_buffer(0,endstring);
    printbuffer(wfile);
    return;
  }
  
  add_to_buffer(offset+16,"transitions");
  sprintf(buffer+strlen(buffer)," := [");
  if (ns == 0 || ns > 128)
      printbuffer(wfile);
  first = TRUE;
  for (i=1;i<=ns;i++){
    if (!first && ns<=128)
      add_to_buffer(offset+21,"");
    sprintf(buffer+strlen(buffer),"[");
    first = FALSE;
    if (!dense && !densepf && i>dr) {
/* In this case we can print directly */
      ptr = table[i]; ptre = table[i+1];
      firstg = TRUE;
      while (ptr < ptre){
        if (!firstg)
          sprintf(buffer+strlen(buffer),",");
        firstg = FALSE;
        g = *(ptr++); k = *(ptr++);
        nl = int_len(k) + 3 + int_len((int)g);
        if (strlen(buffer)+nl >= 76){
          printbuffer(wfile);
/* If the state-set is large, we format less nicely to save space. */
          if (ns<=128)
            add_to_buffer(offset+22,"");
          else
            add_to_buffer(0," ");
        }
        sprintf(buffer+strlen(buffer),"[%d,%d]",g,k);
      }
    }
    else {
      firstg = TRUE;
      for (g=1;g<=ne;g++){
        k = target(dense,table,g,i,dr);
        if (densepf || k!=0){
          if (!firstg)
            sprintf(buffer+strlen(buffer),",");
          firstg = FALSE;
          nl = densepf ? 
                int_len(k) :
                int_len(k) + 3 + int_len((int)g);
          if (strlen(buffer)+nl >= 76){
            printbuffer(wfile);
            if (ns<=128)
              add_to_buffer(offset+22,"");
            else
              add_to_buffer(0," ");
          }
          if (densepf)
            sprintf(buffer+strlen(buffer),"%d",k);
          else
            sprintf(buffer+strlen(buffer),"[%d,%d]",g,k);
        }
      }
    }
    if (i < ns)
      sprintf(buffer+strlen(buffer),"],");
    else
      sprintf(buffer+strlen(buffer),"]");
    printbuffer(wfile);
  }
  add_to_buffer(offset+20,"");
  sprintf(buffer+strlen(buffer),"]");
  printbuffer(wfile);

  if (strlen(name)==0)
     add_to_buffer(0,")");
  else
     add_to_buffer(12+offset,")");
  add_to_buffer(0,endstring);
  printbuffer(wfile);
}

void
fsa_print(wfile,fsaptr,name)
	FILE *wfile;
	fsa *fsaptr;
	char *name;
{ int i, ns, ne;
  boolean first;
  
  if (print_level>=3)
    printf("    #Calling fsa_print.\n");
  ns = fsaptr->states->size;
  ne = fsaptr->alphabet->size;

  *buffer = '\0';
  if (strlen(name)==0)
     add_to_buffer(0,"rec(");
  else{
     sprintf(buffer,"%s := rec(",name);
  }
  printbuffer(wfile);

  add_to_buffer(16,"isFSA");
  sprintf(buffer+strlen(buffer)," := true,");
  printbuffer(wfile);

  srec_print(wfile,fsaptr->alphabet,"alphabet",4,",");

  srec_print(wfile,fsaptr->states,"states",4,",");

  add_to_buffer(16,"flags");
  sprintf(buffer+strlen(buffer)," := [");
  first = TRUE;
  for (i=0;i<num_flag_strings;i++) if (fsaptr->flags[i]){
    if (!first)
      sprintf(buffer+strlen(buffer),",");
    first = FALSE;
    sprintf(buffer+strlen(buffer),"\"%s\"",flag_names[i]);
  }
  sprintf(buffer+strlen(buffer),"],");
  printbuffer(wfile);

  add_to_buffer(16,"initial");
  sprintf(buffer+strlen(buffer)," := [");
  if (ns>1 && ns==fsaptr->num_initial)
    sprintf(buffer+strlen(buffer),"1..%d",fsaptr->states->size);
  else{
    for (i=1;i<=fsaptr->num_initial;i++){
      if (i>1)
        sprintf(buffer+strlen(buffer),",");
      if (strlen(buffer)+int_len(fsaptr->initial[i]) > 76){
        printbuffer(wfile);
        add_to_buffer(21,"");
      }
      sprintf(buffer+strlen(buffer),"%d",fsaptr->initial[i]);
    }
  }
  sprintf(buffer+strlen(buffer),"],");
  printbuffer(wfile);

  add_to_buffer(16,"accepting");
  sprintf(buffer+strlen(buffer)," := [");
  if (ns>1 && ns==fsaptr->num_accepting)
    sprintf(buffer+strlen(buffer),"1..%d",fsaptr->states->size);
  else{
    for (i=1;i<=fsaptr->num_accepting;i++){
      if (i>1)
        sprintf(buffer+strlen(buffer),",");
      if (strlen(buffer)+int_len(fsaptr->accepting[i]) > 76){
        printbuffer(wfile);
        add_to_buffer(21,"");
      }
      sprintf(buffer+strlen(buffer),"%d",fsaptr->accepting[i]);
    }
  }
  sprintf(buffer+strlen(buffer),"],");
  printbuffer(wfile);

/* If the fsa is not known to be deterministic, we will print it sparsely */
  if (!fsaptr->flags[DFA])
    fsaptr->table->printing_format = SPARSE;
  table_print(wfile,fsaptr->table,"table",4,"",ns,ne);

  sprintf(buffer,");");
  printbuffer(wfile);
}

void
srec_read(rfile,srptr,maxsize)
/* Read the set record *srptr from rfile, assigning space as required.
 * If maxsize is larger than srptr->size, and space is allocated for
 * names or labels, then space is allocated for maxsize of these.
 * This allows for possible later expansion.
 */
	FILE *rfile;
	srec *srptr;
	int maxsize;
	
{ int delim, i, n, r;
  boolean typeset, sizeset, baseset, arityset, paddingset, namesset, formatset,
          labelsset, alphabetset, setToLabelsset;
  storage_type st;
  read_ident(rfile,buffer,&delim,FALSE);
  if (delim != '(' || strcmp(buffer,"rec")!=0){
    fprintf(stderr,
        "#Input error: \"rec(\" expected\n");
    exit(1);
  }

/* main loop reading the fields of the record follows. */
  typeset = FALSE;
  sizeset = FALSE;
  baseset = FALSE;
  arityset = FALSE;
  paddingset = FALSE;
  namesset = FALSE;
  formatset = FALSE;
  labelsset = FALSE;
  alphabetset = FALSE;
  setToLabelsset = FALSE;
  do{
    read_ident(rfile,buffer,&delim,FALSE);
    if (delim != ':'){
      fprintf(stderr,
          "#Input error: bad record field assignment\n");
      exit(1);
    }
    check_next_char(rfile,'=');
    if (strcmp(buffer,"type")==0){
      typeset = TRUE;
      read_string(rfile,buffer,&delim);
      if (strcmp(buffer,"simple")==0)
        srptr -> type = SIMPLE;
      else if (strcmp(buffer,"identifiers")==0)
        srptr -> type = IDENTIFIERS;
      else if (strcmp(buffer,"words")==0)
        srptr -> type = WORDS;
      else if (strcmp(buffer,"strings")==0)
        srptr -> type = STRINGS;
      else if (strcmp(buffer,"product")==0)
        srptr -> type = PRODUCT;
      else if (strcmp(buffer,"labeled")==0 || strcmp(buffer,"labelled")==0)
        srptr -> type = LABELED;
      else{
       fprintf(stderr,"Error: unknown set-record type %s.\n",buffer);
       exit(1);
      }
    }
    else if (strcmp(buffer,"size")==0){
      sizeset = TRUE;
      read_int(rfile,&(srptr->size),&delim);
    }
    else if (strcmp(buffer,"arity")==0){
      arityset = TRUE;
      if (!typeset || srptr->type != PRODUCT){
        fprintf(stderr,"Error: arity field only used for product type.\n");
        exit(1);
      }
      read_int(rfile,&(srptr->arity),&delim);
    }
    else if (strcmp(buffer,"padding")==0){
      paddingset = TRUE;
      if (!typeset || srptr->type != PRODUCT){
        fprintf(stderr,"Error: padding field only used for product type.\n");
        exit(1);
      }
      read_ident(rfile,buffer,&delim,FALSE);
      srptr->padding = buffer[0];
    }
    else if (strcmp(buffer,"base")==0){
      baseset = TRUE;
      if (!typeset || srptr->type != PRODUCT){
        fprintf(stderr,
                "Error: base field only used for product type.\n");
        exit(1);
      }
      tmalloc(srptr->base,srec,1);
      srec_read(rfile,srptr->base,0);
      read_delim(rfile,&delim);
    }
    else if (strcmp(buffer,"format")==0){
      formatset = TRUE;
      if (srptr->type!=IDENTIFIERS && srptr->type!=WORDS &&
          srptr->type!=STRINGS && srptr->type!=LABELED){
       fprintf(stderr,"Error: set-record type doesn't require format field.\n");
       exit(1);
      }
      read_string(rfile,buffer,&delim);
      if (strcmp(buffer,"dense")==0)
        st = DENSE;
      else if (strcmp(buffer,"sparse")==0)
        st = SPARSE;
      else{
       fprintf(stderr,"Error: unknown storage type %s.\n",buffer);
       exit(1);
      }
    }
    else if (strcmp(buffer,"alphabet")==0){
      alphabetset = TRUE;
      if (srptr->type!=WORDS){
       fprintf(stderr,
                   "Error: set-record type doesn't require alphabet field.\n");
       exit(1);
      }
      check_next_char(rfile,'[');
      i=0;
      do {
         read_ident(rfile,buffer,&delim,TRUE);
         if (strlen(buffer)>0) {
           i++;
           if (i>MAXCHAR) {
              fprintf(stderr,"Error: alphabet can have size at most 127.\n");
              exit(1);
           }
           tmalloc(srptr->alphabet[i],char,strlen(buffer)+1);
           strcpy(srptr->alphabet[i],buffer);
         }
      } while (delim!=']');
      srptr->alphabet_size = i;
      check_next_char(rfile,',');
    }
    else if (strcmp(buffer,"names")==0){
      if (!typeset || !sizeset || !formatset){
        fprintf(stderr,
                "Error: type, size and format fields must precede names.\n");
        exit(1);
      }
      if (srptr->type!=IDENTIFIERS && srptr->type!=WORDS &&
          srptr->type!=STRINGS){
        fprintf(stderr,"Error: set-record type doesn't require names field.\n");
        exit(1);
      }
      namesset = TRUE;
      check_next_char(rfile,'[');
      if (maxsize < srptr->size)
        maxsize = srptr->size;
      tmalloc(srptr->names,char *,maxsize + 1);
      if (srptr->type==WORDS){
        if (!alphabetset){
          fprintf(stderr,"Error: alphabet field must precede names.\n");
          exit(1);
        }
/* We now have to work out the algorithm for converting words to generator
 * numbers, based on the alphabet names - similar to that in rwsio.c
 */
        process_names(srptr->alphabet,srptr->alphabet_size);
/* Finally allocate some space temporarily to read the words into. */
        tmalloc(wbuffer,char,MAXWORDLEN+1);
      }
      if (st==DENSE){
        if (srptr->size==0)
           read_delim(rfile,&delim);
        for (i=1;i<=srptr->size;i++){
          if (srptr->type==IDENTIFIERS){
            read_ident(rfile,buffer,&delim,TRUE);
            tmalloc(srptr->names[i],char,strlen(buffer)+1);
            strcpy(srptr->names[i],buffer);
          }
          else if (srptr->type==STRINGS){
            read_string(rfile,buffer,&delim);
            tmalloc(srptr->names[i],char,strlen(buffer)+1);
            strcpy(srptr->names[i],buffer);
          }
          else if (srptr->type==WORDS){
            read_word(rfile,wbuffer,wbuffer+MAXWORDLEN,&delim,
                      srptr->alphabet,srptr->alphabet_size,TRUE);
            tmalloc(srptr->names[i],char,strlen(wbuffer)+1);
            strcpy(srptr->names[i],wbuffer);
          }
        }
      }
      else {
        read_delim(rfile,&delim);
        while (delim != ']') {
          read_int(rfile,&i,&delim);
          if (i<=0 || i>srptr->size){
            fprintf(stderr,"Error: name-number out of range.\n");
            exit(1);
          }
          if (srptr->type==IDENTIFIERS){
            read_ident(rfile,buffer,&delim,TRUE);
            tmalloc(srptr->names[i],char,strlen(buffer)+1);
            strcpy(srptr->names[i],buffer);
          }
          else if (srptr->type==STRINGS){
            read_string(rfile,buffer,&delim);
            tmalloc(srptr->names[i],char,strlen(buffer)+1);
            strcpy(srptr->names[i],buffer);
          }
          else if (srptr->type==WORDS){
            read_word(rfile,wbuffer,wbuffer+MAXWORDLEN,&delim,
                      srptr->alphabet,srptr->alphabet_size,TRUE);
            tmalloc(srptr->names[i],char,strlen(wbuffer)+1);
            strcpy(srptr->names[i],wbuffer);
          }
          read_delim(rfile,&delim);
          if (delim == ',')
            read_delim(rfile,&delim);
        }
      }
      read_delim(rfile,&delim);
/* don't forget to de-allocate space for reading words if necessary */
      if (srptr->type==WORDS)
        tfree(wbuffer);
    }
    else if (strcmp(buffer,"labels")==0){
      labelsset = TRUE;
      if (!typeset || (srptr->type != LABELED)){
        fprintf(stderr,
                "Error: labels field only used for labeled type.\n");
        exit(1);
      }
      tmalloc(srptr->labels,srec,1);
      srec_read(rfile,srptr->labels,0);
      read_delim(rfile,&delim);
      if (srptr->labels->size > MAXUSHORT){
        fprintf(stderr,"Error: label-set can have size at most 65535.\n");
        exit(1);
      }
    }
    else if (strcmp(buffer,"setToLabels")==0){
      setToLabelsset = TRUE;
      if (!typeset || (srptr->type != LABELED)){
        fprintf(stderr,
                "Error: setToLabels field only used for labeled type.\n");
        exit(1);
      }
      if (!labelsset){
        fprintf(stderr,"Error: labels field must precede setToLabels.\n");
        exit(1);
      }
      if (!formatset){
        fprintf(stderr,"Error: format field must precede setToLabels.\n");
        exit(1);
      }
      check_next_char(rfile,'[');
      if (maxsize < srptr->size)
        maxsize = srptr->size;
      tmalloc(srptr->setToLabels,setToLabelsType,maxsize+1);
      for (i=0;i<=srptr->size;i++)
        srptr->setToLabels[i] = 0;
      if (st==DENSE){
        for (i=1;i<=srptr->size;i++){
          read_int(rfile,&n,&delim);
          if (n<0 || n>srptr->labels->size) {
            fprintf(stderr,"Error: label out of range.\n");
            exit(1);
          }
          srptr->setToLabels[i] = n;
          if (delim==']')
            break;
        }
      }
      else {
        read_delim(rfile,&delim);
        while (delim != ']') {
          read_int(rfile,&i,&delim);
          if (i<=0 || i>srptr->size){
            fprintf(stderr,"Error: label-number out of range.\n");
            exit(1);
          }
          read_int(rfile,&n,&delim);
          if (n<0 || n>srptr->labels->size) {
            fprintf(stderr,"Error: label out of range.\n");
            exit(1);
          }
          srptr->setToLabels[i] = n;
          read_delim(rfile,&delim);
          if (delim == ',')
            read_delim(rfile,&delim);
        }
      }
      read_delim(rfile,&delim);
    }
    else { 
      fprintf(stderr,"#Warning: Invalid field name %s.\n",buffer);
      exit(1);
    }
  } while (delim != ')');

  if (!typeset || !sizeset){
    fprintf(stderr,"Error: type and size fields must be set.\n");
    exit(1);
  }
  if (srptr->type==IDENTIFIERS || srptr->type==WORDS ||
     srptr->type==STRINGS){
    if (!namesset){
       fprintf(stderr,"Error: set-record type requires names field.\n");
       exit(1);
    }
  }
  if (srptr->type==PRODUCT) {
    if (!arityset || !paddingset || !baseset){
       fprintf(stderr,
         "Error: set-record type require arity, padding and base fields.\n");
       exit(1);
    }
    /* The size of the set-record should be (n+1)^r or (n+1)^r-1, where
     * n is the size of the base set, and r the arity.
     */
       n = srptr->base->size + 1;
       for (i=1;i<srptr->arity;i++)
         n *= (srptr->base->size + 1);
       if (srptr->size != n && srptr->size != n-1){
         fprintf(stderr,
           "Error: set-record size wrong for product type.\n");
         exit(1);
       }
  }
  if (srptr->type==LABELED) {
    if (!setToLabelsset){
       fprintf(stderr,
         "Error: set-record type require setToLabels field.\n");
       exit(1);
    }
  }
}

void
table_read(rfile,tableptr,table_storage_type,dr,ns,maxstates,ne)
/* Read the table_struc *tableptr from rfile, assigning space as required.
 * dr is the number of rows stored densely if storage_type=SPARSE.
 * ns and ne are the sizes of the state-set and alphabet-set.
 * maxstates is the maximum number of states that we allocate space for
 * in dense-storage mode.
 */
	FILE *rfile;
	table_struc *tableptr;
	storage_type table_storage_type;
	int dr;
	int ns, maxstates, ne;
	
{ int delim, i, j, k, ntleft, dt = 0;
  boolean  filenameset, formatset, numTransitionsset, transitionsset;
  int **table, *sparseptr;

  formatset = FALSE;
  numTransitionsset = FALSE;
  transitionsset = FALSE;
  filenameset = FALSE;

  read_ident(rfile,buffer,&delim,FALSE);
  if (delim != '(' || strcmp(buffer,"rec")!=0){
    fprintf(stderr,
        "#Input error: \"rec(\" expected\n");
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
    if (strcmp(buffer,"filename")==0){
/* In this case, the transition table (and possibly other info) is stored
 * externally in the file "filename"
 */
      if (transitionsset) {
        fprintf(stderr,
         "#Input error - filename and transitions fields are both present.\n");
        exit(1);
      } 
      filenameset = TRUE;
      read_string(rfile,buffer,&delim);
      tmalloc(tableptr->filename,char,strlen(buffer)+1);
      strcpy(tableptr->filename,buffer);
    }
    else if (strcmp(buffer,"format")==0){
      formatset = TRUE;
      read_string(rfile,buffer,&delim);
      if (strcmp(buffer,"dense deterministic")==0)
        tableptr->printing_format = DENSE;
      else if (strcmp(buffer,"sparse")==0)
        tableptr->printing_format = SPARSE;
      else {
        fprintf(stderr,"Error: invalid format value %s\n",buffer);
        exit(1);
      }
    }
    else if (strcmp(buffer,"numTransitions")==0){
      numTransitionsset = TRUE;
      read_int(rfile,&(tableptr->numTransitions),&delim);
    }
    else if (strcmp(buffer,"defaultTarget")==0){
      read_int(rfile,&dt,&delim);
      if (dt != 0)
/* Our code doesn't currently cope with this, so we will enforce storage-type
 * dense in this case.
 */
        table_storage_type = DENSE;
    }
    else if (strcmp(buffer,"transitions")==0){
      if (filenameset) {
        fprintf(stderr,
         "#Input error - filename and transitions fields are both present.\n");
        exit(1);
      } 
      transitionsset = TRUE;
      if (!numTransitionsset){
         /* In this case, we have currently to use dense storage */
         table_storage_type = DENSE;
      }
      tableptr->table_type = table_storage_type;
      if (table_storage_type==SPARSE) {
        if (dr>ns) dr=ns;
        tableptr->denserows = dr;
      }
      tableptr->maxstates = maxstates;
/*  assign space according to storage type */
      if (table_storage_type == DENSE){
        fsa_table_init(tableptr,maxstates,ne);
        table = tableptr->table_data_ptr;
      }
      else {
          tmalloc(tableptr->table_data_ptr,int *,ns+2);
          table = tableptr->table_data_ptr;
          if (dr > 0) {
/* Here we have to read in the data in two separate chunks. */
            tmalloc(table[0],int,dr*ne);
            for (i=1;i<=dr;i++)
              table[i] = table[0] + (i-1)*ne - 1;
            ntleft = tableptr->numTransitions;
          }
          else {
            tmalloc(table[0],int,2*tableptr->numTransitions+1);
            sparseptr = table[0];
          }
      }
/* Now read the transitions  */
      check_next_char(rfile,'[');
      for (i=1;i<=ns;i++){
        check_next_char(rfile,'[');
        if (table_storage_type == SPARSE) {
          if (dr>0 && i==dr+1) { /* allocate remaining storage */
            tmalloc(table[i],int,2*ntleft+1);
            sparseptr = table[i];
          }
          else if (i > dr)
            table[i] = sparseptr;
        }
        if (tableptr->printing_format == DENSE){
          if (ne==0)
            read_delim(rfile,&delim);
          else for (j=1;j<=ne;j++){
            read_int(rfile,&k,&delim);
            if (table_storage_type == DENSE)
              table[j][i] = k;
            else if (i<=dr) {
              table[i][j] = k;
              if (k!=0) ntleft--;
            }
            else if (k!=0) {
              *(sparseptr++)=j; *(sparseptr++)=k;
            }
          }
          if (delim != ']'){
            fprintf(stderr,"Error: format error in table.\n");
            exit(1);
          }
        }
        else {
          if (table_storage_type == DENSE)
            for (j=1;j<=ne;j++)
              table[j][i] = dt; /* in case there is a default target */
          else if (i<=dr)
            for (j=1;j<=ne;j++)
              table[i][j] = 0;
          read_delim(rfile,&delim);
          while (delim != ']') {
            read_int(rfile,&j,&delim);
            read_int(rfile,&k,&delim);
            if (table_storage_type == DENSE)
              table[j][i] = k;
            else if (i<=dr) {
              table[i][j] = k;
              if (k!=0) ntleft--;
            }
            else if (k!=0) {
              *(sparseptr++)=j; *(sparseptr++)=k;
            }
            read_delim(rfile,&delim);
            if (delim == ',')
              read_delim(rfile,&delim);
          }
        } 
        read_delim(rfile,&delim);
      }
      if (ns==0)
        read_delim(rfile,&delim);
      if (table_storage_type == SPARSE && ns>dr)
/* Set extra pointer to show end of data */
        table[ns+1] = sparseptr;
      read_delim(rfile,&delim);
    }
    else 
      fprintf(stderr,"#Warning: Invalid field name %s.\n",buffer);
  } while (delim != ')');

  if (!transitionsset && !filenameset){
    fprintf(stderr,"transitions and filename fields are not set.\n");
    exit(1);
  }
  if (filenameset) 
/* In this case, the storage_type has to be as specified by the format field */
    tableptr->table_type = tableptr->printing_format;
}

void
fsa_read(rfile,fsaptr,table_storage_type,dr,maxstates,assignment,name)
/* Read the fsa *fsaptr from rfile, assigning space as required.
 * If maxstates is greater than the actual number of states, then
 * space will be assigned for maxstates states, to allow more states
 * to be added later (this only makes sense with dense storage_type).
 * If assignment is true, an assignment to an fsa is read.
 * The name is returned in *name, which is assumed to have enough space.
 */
	FILE *rfile;
	fsa *fsaptr;
	storage_type table_storage_type;
	int dr;
	int maxstates;
	boolean assignment;
	char *name;
	
{ int delim, i, j, k, ct, ns, ne, *ia, *buf1, *buf2, *swapbuf, bufsize;
  boolean isFSAset, statesset, alphabetset, initialset, acceptingset,
          flagsset, tableset, in, compressed;

  if (print_level>=3)
    printf("    #Calling fsa_read.\n");
  isFSAset = FALSE;
  isFSAset = FALSE;
  statesset = FALSE;
  alphabetset = FALSE;
  initialset = FALSE;
  acceptingset = FALSE;
  flagsset = FALSE;
  tableset = FALSE;

  fsa_init(fsaptr);

  if (assignment){
    read_ident(rfile,buffer,&delim,FALSE);
    if (delim != ':'){
       fprintf(stderr,
        "#Input error: file must contain a record assignment\n");
       exit(1);
     }
     strcpy(name,buffer);
     check_next_char(rfile,'=');
  }

  read_ident(rfile,buffer,&delim,FALSE);
  if (delim != '(' || strcmp(buffer,"rec")!=0){
    fprintf(stderr,
        "#Input error: \"rec(\" expected\n");
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
    if (strcmp(buffer,"isFSA")==0){
      isFSAset = TRUE;
      read_ident(rfile,buffer,&delim,FALSE);
      if (strcmp(buffer,"true")!=0){
          fprintf(stderr,
              "#Input error: isFSA field must equal \"true\"\n");
          exit(1);
      }
    }
    else if (strcmp(buffer,"states")==0){
      statesset = TRUE;
      srec_read(rfile,fsaptr->states,maxstates);
      ns = fsaptr->states->size;
      read_delim(rfile,&delim);
    }
    else if (strcmp(buffer,"alphabet")==0){
      alphabetset = TRUE;
      srec_read(rfile,fsaptr->alphabet,0);
      ne = fsaptr->alphabet->size;
      read_delim(rfile,&delim);
    }
    else if (strcmp(buffer,"flags")==0){
      flagsset = TRUE;
      check_next_char(rfile,'[');
      do {
        read_string(rfile,buffer,&delim);
        for (i=0;i<num_flag_strings;i++)
          if (strcmp(buffer,flag_names[i])==0)
            fsaptr->flags[(flag_strings)i]=TRUE;
      } while (delim != ']');
      read_delim(rfile,&delim);
/* If the fsa is not known to be deterministic, then we should store is in
 * sparse format
 */
      if (!fsaptr->flags[DFA]) {
        table_storage_type = SPARSE;
        dr = 0;
      }
    }
    else if (strcmp(buffer,"initial")==0 || strcmp(buffer,"accepting")==0){
/* We'll handle these together since they are so similar. */
      compressed = FALSE;
      bufsize = 1024;
      in = strcmp(buffer,"initial")==0;
      if (!statesset){
        fprintf(stderr,
          "Error: states field must precede initial and accepting fields.\n");
        exit(1);
      }
      ct = 0;
      tmalloc(buf1,int,bufsize);
/* We don't know how long the list will be in advance */
      check_next_char(rfile,'[');
      do {
        read_int(rfile,&i,&delim);
        if (i<0 || i>ns){
         fprintf(stderr,"Invalid state number in initial or accepting list.\n");
          exit(1);
        }
        if (i>0){
          ct++;
          if (delim == '.'){
             compressed = TRUE;
             if (ct!=1){
               fprintf(stderr,
                        "Error in format of initial or accepting field.\n");
               exit(1);
             }
             check_next_char(rfile,'.');
             read_int(rfile,&j,&delim);
             if (j>=i && j<=ns){
               if (delim!=']'){
                 fprintf(stderr,
                         "Error in format of initial or accepting field.\n");
                 exit(1);
               }
               ct = j-i+1;
               if (ct==1 || ct!=ns){
                   /* if ct = ns, we don't need to store the list. */
                 tmalloc(ia,int,ct+1);
                 if (in)
                   fsaptr->initial = ia;
                 else
                    fsaptr->accepting = ia;
                 for (k=1;k<=ct;k++)
                   ia[k] = k+i-1;
               }
             }
             else {
               fprintf(stderr,
                       "Error in format of initial or accepting field.\n");
               exit(1);
             }
          }
          else {
             compressed = FALSE;
             if (ct>=bufsize){
                bufsize*=2;
                tmalloc(buf2,int,bufsize);
                for (i=1;i<ct;i++)
                  buf2[i] = buf1[i];
                tfree(buf1);
                swapbuf=buf1; buf1=buf2; buf2=swapbuf;
             }
             buf1[ct] = i;
          }
        } /* if (i>0) */
      } while (delim != ']');
      if (!compressed){
        tmalloc(ia,int,ct+1);
        if (in)
          fsaptr->initial = ia;
        else
          fsaptr->accepting = ia;
        for (i=1;i<=ct;i++)
          ia[i] = buf1[i];
      }
      tfree(buf1);
      if (in){
        initialset = TRUE;
        fsaptr->num_initial = ct;
      }
      else{
        acceptingset = TRUE;
        fsaptr->num_accepting = ct;
      }
      read_delim(rfile,&delim);
    }
    else if (strcmp(buffer,"table")==0){
      if (!statesset){
        fprintf(stderr,
             "Error: states field must precede table fields.\n");
        exit(1);
      }
      if (!alphabetset){
        fprintf(stderr,
             "Error: alphabet field must precede table fields.\n");
        exit(1);
      }
      tableset = TRUE;
      if (fsaptr->alphabet->type == PRODUCT)
        k = fsaptr->alphabet->base->size; 
      else
        k = -1;
      if (maxstates < ns)
        maxstates = ns;
      table_read(rfile,fsaptr->table,table_storage_type,dr,ns,maxstates,ne,k);
      read_delim(rfile,&delim);
    }
    else 
      fprintf(stderr,"#Warning: Invalid field name %s.\n",buffer);
  } while (delim != ')');

  if (!isFSAset || !statesset || !alphabetset || !initialset ||
      !acceptingset || !flagsset || !tableset){
    fprintf(stderr,"One of the compulsory fsa-fields is not set.\n");
    exit(1);
  }
  if (assignment)
    read_delim(rfile,&delim);
}

void
compressed_transitions_read(fsaptr,rfile)
	fsa *fsaptr;
        FILE *rfile;
/* Read the transition-table of the fsa *fsaptr which is stored in the
 * file with *file in unformatted form.
 * *file should be opened and closed externally to the function.
 * This is used by the functions in fsalogic.c
 */
{ int ns, ne, nt, **table, *tab_ptr, len, i, j; 

  ns = fsaptr->states->size;
  ne = fsaptr->alphabet->size;
  nt = fsaptr->table->numTransitions;

  if (fsaptr->table->table_type==DENSE) {
/* it is a pity that the table was output in transposed form! */
    fsa_table_init(fsaptr->table,ns,ne);
    table = fsaptr->table->table_data_ptr;
    for (i=1;i<=ns;i++) for(j=1;j<=ne;j++)
      fread((void *)(table[j]+i),sizeof(int),(size_t)1,rfile);
  }
  else {
    tmalloc(fsaptr->table->table_data_ptr,int *,ns+2);
    tmalloc(fsaptr->table->table_data_ptr[0],int,2*nt+1);
    tab_ptr = fsaptr->table->table_data_ptr[0];
    table = fsaptr->table->table_data_ptr;
    for (i=1;i<=ns;i++) {
      table[i] = tab_ptr;
      fread((void *)&len,sizeof(int),(size_t)1,rfile);
      if (len>0)
        fread((void *)tab_ptr,sizeof(int),(size_t)len,rfile);
      tab_ptr += len;
    }
    table[ns+1] = tab_ptr;
  }
}
