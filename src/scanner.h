/****************************************************************************
**
*A  scanner.h                   GAP source                   Martin Schoenert
**
*A  @(#)$Id: scanner.h,v 1.1.1.1 1996/12/11 12:43:59 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file declares the functions of the scanner, which is responsible for
**  all input and output processing.
**
**  The scanner  exports two very  important abstractions.  The  first is the
**  concept that an input file is  a stream of symbols,  such nasty things as
**  <space>,  <tab>,  <newline> characters or  comments (they are worst  :-),
**  characters making  up identifiers  or  digits that  make  up integers are
**  hidden from the rest of GAP.
**
**  The second is  the concept of  a current input  and output file.   In the
**  main   module   they are opened  and   closed  with the  'OpenInput'  and
**  'CloseInput' respectively  'OpenOutput' and 'CloseOutput' calls.  All the
**  other modules just read from the  current input  and write to the current
**  output file.
**
**  The scanner relies on the functions  provided  by  the  operating  system
**  dependent module 'system.c' for the low level input/output.
**
*H  $Log: scanner.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:59  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.4  1993/05/05  11:10:12  fceller
*H  added 'LogInputTo'
*H
*H  Revision 3.3  1992/12/16  19:46:21  martin
*H  added character constants
*H
*H  Revision 3.2  1992/12/08  11:50:26  martin
*H  added '<list>{<positions>}'
*H
*H  Revision 3.1  1991/04/30  16:12:45  martin
*H  initial revision under RCS
*H
*/


/****************************************************************************
**
*V  Symbol  . . . . . . . . . . . . . . . . .  current symbol read from input
**
**  The  variable 'Symbol' contains the current  symbol read from  the input.
**  It is represented as an unsigned long integer.
**
**  The possible values for 'Symbol' are defined in the  definition  file  of
**  this package as follows:
*/
#define S_ILLEGAL       (0L)

#define S_IDENT         ((1L<< 3))
#define S_INT           ((1L<< 4))
#define S_CHAR          ((1L<< 5)+0)
#define S_STRING        ((1L<< 5)+1)

#define S_DOT           ((1L<< 6))
#define S_LBRACK        ((1L<< 7)+0)
#define S_RBRACK        ((1L<< 8)+0)
#define S_LBRACE        ((1L<< 7)+1)
#define S_RBRACE        ((1L<< 8)+1)
#define S_LPAREN        ((1L<< 9))
#define S_RPAREN        ((1L<<10))
#define S_COMMA         ((1L<<11)+0)
#define S_DOTDOT        ((1L<<11)+1)

#define S_IF            ((1L<<12)+0)
#define S_THEN          ((1L<<13))
#define S_ELIF          ((1L<<14)+0)
#define S_ELSE          ((1L<<14)+1)
#define S_FI            ((1L<<15))

#define S_FOR           ((1L<<12)+1)
#define S_DO            ((1L<<16))
#define S_OD            ((1L<<17))

#define S_REPEAT        ((1L<<12)+2)
#define S_UNTIL         ((1L<<18))
#define S_WHILE         ((1L<<12)+3)

#define S_ASSIGN        ((1L<<19))
#define S_SEMICOLON     ((1L<<20))

#define S_FUNCTION      ((1L<<21))
#define S_LOCAL         ((1L<<22))
#define S_END           ((1L<<23))
#define S_RETURN        ((1L<<12)+4)
#define S_MAPTO         ((1L<<24))

#define S_NOT           ((1L<<25)+0)
#define S_AND           ((1L<<25)+1)
#define S_OR            ((1L<<25)+2)

#define S_EQ            ((1L<<26)+0)
#define S_LT            ((1L<<26)+1)
#define S_GT            ((1L<<26)+2)
#define S_NE            ((1L<<26)+3)
#define S_LE            ((1L<<26)+4)
#define S_GE            ((1L<<26)+5)
#define S_IN            ((1L<<26)+6)

#define S_PLUS          ((1L<<27)+0)
#define S_MINUS         ((1L<<27)+1)

#define S_MULT          ((1L<<28)+0)
#define S_DIV           ((1L<<28)+1)
#define S_MOD           ((1L<<28)+2)
#define S_POW           ((1L<<28)+3)

#define S_QUIT          ((1L<<29))
#define S_EOF           ((1L<<30))

extern  unsigned long   Symbol;


/****************************************************************************
**
*T  TypSymbolSet  . . . . . . . . . . . . . . . . . . type of sets of symbols
**
**  'TypSymbolSet' is the type of sets of symbols.  Sets  of symbols are used
**  in the error recovery of the  parser  to specify that 'Match' should skip
**  all symbols until finding one in a specified set.
**
**  If there were less than 32 different symbols  things would be  very easy.
**  We could  simply assign   the  symbolic constants   that are the possible
**  values for 'Symbol' values 1, 2, 4, 8, 16, ...  and so on.  Then making a
**  set  would  simply mean  or-ing the  values, as in  'S_INT|S_STRING', and
**  checking whether a symbol is in a set would be '(<symbol> & <set>) != 0'.
**
**  There  are however more  than 32 different  symbols, so  we must  be more
**  clever.  We  group some  symbols that  are syntactically  equivalent like
**  '*', '/' in a class. We use the least significant 3 bits to differentiate
**  between members in one class.  And now  every symbol class, many of which
**  contain   just  one  symbol,  has exactely  one   of  the  remaining most
**  significant 29  bits  set.   Thus   sets  of symbols  are  represented as
**  unsigned long integers, which is typedef-ed to 'TypSymbolSet'.
**
**  The classes are as follows, all other symbols are in a class themself:
**      if, for, repeat, while, return
**      elif, else
**      not, and, or
**      =, <>, <, >=, <=, >, in
**      +, -
**      *, /, mod, ^
**
**  'TypSymbolSet'  is defined in the   definition  file of  this  package as
**  follows:
*/
typedef unsigned long   TypSymbolSet;


/****************************************************************************
**
*F  IS_IN( <symbol>, <set> )  . . . . . . . . is a symbol in a set of symbols
**
**  'IS_IN' returns 1 if the symbol <symbol> is in the symbol set <set> and 0
**  otherwise.  Due to the grouping into classes some symbol sets may contain
**  more than mentioned, for  example 'IS_IN(S_POW,S_MULT|S_DIV|S_MOD)' is 1.
**
**  'IS_IN' is defined in the definition file of this package as follows:
*/
#define IS_IN(SYMBOL,SET)       ((SYMBOL) & ((SET) & ~7))


/****************************************************************************
**
*V  EXPRBEGIN . . . . . . . . . . . . set of symbols that start an expression
*V  STATBEGIN . . . . . . . . . . . . . set of symbols that start a statement
**
**  'EXPRBEGIN' is  the set   of symbols   that might  start   an expression.
**  'STATBEGIN' is the set of symbols that might  start a stament, this  is a
**  superset of 'EXPRBEGIN', since expressions are themselfs statments.
**
**  'EXPRBEGIN' and 'STATBEGIN'  are defined in  the definition  file of this
**  package as follows:
*/
#define EXPRBEGIN  (S_IDENT|S_INT|S_STRING|S_LPAREN|S_FUNCTION)
#define STATBEGIN  (EXPRBEGIN|S_IF|S_FOR|S_WHILE|S_REPEAT|S_RETURN)


/****************************************************************************
**
*V  Value . . . . . . . . . . . .  value of the identifier, integer or string
**
**  If 'Symbol' is 'S_IDENT', 'S_INT' or 'S_TRING' the variable 'Value' holds
**  the name of the identifier, the digits of the integer or the value of the
**  string constant.
**
**  Note that the size of  'Value' limits the  maximal number of  significant
**  characters  of an identifier,   the maximal size  of  an  integer and the
**  maximal length of a  string.   'GetIdent', 'GetInt' and 'GetStr' truncate
**  identifier, integers or strings after that many characters.
*/
#define MAXVALUE "2048"
char            Value [2048];

/****************************************************************************
**
*V  NrError . . . . . . . . . . . . . . . .  number of errors in current expr
*V  NrErrLine . . . . . . . . . . . . . . .  number of errors on current line
**
**  'NrError' is an integer whose value is the number of errors already found
**  in the current expression.  It is set to 0 at the beginning of 'Read' and
**  incremented with each 'SyntaxError' call, including those  from  'Match'.
**
**  If 'NrError' is greater than zero the parser functions  will  not  create
**  new bags.  This prevents the parser from creating new bags after an error
**  occured.
**
**  'NrErrLine' is an integer whose value is the number of  errors  found  on
**  the current line.  It is set to 0 in 'GetLine' and incremented with  each
**  'SyntaxError' call, including those from 'Match'.
**
**  If 'NrErrLine' is greater  than  zero  'SyntaxError' will  not  print  an
**  error message.  This prevents the printing of multiple error messages for
**  one line, since they  probabely  just reflect  the  fact that the  parser
**  has not resynchronized yet.
*/
extern  long            NrError;
extern  long            NrErrLine;


/****************************************************************************
**
*V  Prompt  . . . . . . . . . . . . . . . . . . . . . .  prompt to be printed
**
**  'Prompt' holds the string that is to be printed if a  new  line  is  read
**  from the interactive files '*stdin*' or '*errin*'.
**
**  It is set to 'gap> ' or 'brk> ' in the  read-eval-print loops and changed
**  to the partial prompt '> ' in 'Read' after the first symbol is read.
*/
extern  char            * Prompt;


/****************************************************************************
**
*F  SyntaxError( <msg> )  . . . . . . . . . . . . . . .  raise a syntax error
**
**  'SyntaxError' prints the current line, followed by the error message:
**
**      ^ syntax error, <msg> in <current file name>
**
**  with the '^' pointing to the current symbol on the current line.  If  the
**  <current file name> is '*stdin*' it is not printed.
**
**  'SyntaxError' is called from the parser to print error messages for those
**  errors that are not cought by 'Match',  for example if the left hand side
**  of an assignment is not a variable, a list element or a record component,
**  or if two formal arguments of a function have the same identifier.  It is
**  also called for warnings, for example if a statement has no effect.
**
**  'SyntaxError' first increments 'NrError' by   1.  If 'NrError' is greater
**  than zero the parser functions  will not create  new bags.  This prevents
**  the parser from creating new bags after an error occured.
**
**  'SyntaxError'  also  increments  'NrErrLine'  by  1.  If  'NrErrLine'  is
**  greater than zero  'SyntaxError' will not print an  error  message.  This
**  prevents the printing of multiple error messages for one line, since they
**  probabely  just reflect the  fact  that the parser has not resynchronized
**  yet.  'NrErrLine' is reset to 0 if a new line is read in 'GetLine'.
*/
void            SyntaxError ( char * msg );


/****************************************************************************
**
*F  Match( <symbol>, <msg>, <skipto> )  . match current symbol and fetch next
**
**  'Match' is the main  interface between the  scanner and the  parser.   It
**  performs the  4 most common actions in  the scanner  with  just one call.
**  First it checks that  the current symbol stored  in the variable 'Symbol'
**  is the expected symbol  as passed in the  argument <symbol>.  If  it  is,
**  'Match' reads the next symbol from input  and returns.  Otherwise 'Match'
**  first prints the current input line followed by the syntax error message:
**  '^ syntax error, <msg> expected' with '^' pointing to the current symbol.
**  It then  skips symbols up to one  in the resynchronisation  set <skipto>.
**  Actually 'Match' calls 'SyntaxError' so its comments apply here too.
**
**  One kind of typical 'Match' call has the form
**
**      'Match( Symbol, "", 0L );'.
**
**  This is used if the parser knows that the current  symbol is correct, for
**  example in 'RdReturn'  the   first symbol must be 'S_RETURN',   otherwise
**  'RdReturn' would not have been  called.  Called this  way 'Match' will of
**  course never raise an syntax error,  therefore <msg>  and <skipto> are of
**  no concern, they are passed nevertheless  to please  lint.  The effect of
**  this call is merely to read the next symbol from input.
**
**  Another typical 'Match' call is in 'RdIf' after we read the if symbol and
**  the condition following, and now expect to see the 'then' symbol:
**
**      Match( S_THEN, "then", STATBEGIN|S_ELIF|S_ELSE|S_FI|follow );
**
**  If the current symbol  is 'S_THEN' it is  matched  and the next symbol is
**  read.  Otherwise 'Match'  prints the  current line followed by the  error
**  message: '^ syntax error, then expected'.  Then 'Match' skips all symbols
**  until finding either  a symbol  that can begin  a statment,  an 'elif' or
**  'else' or 'fi' symbol, or a symbol that is  contained in the set <follow>
**  which is passed to  'RdIf' and contains  all symbols allowing  one of the
**  calling functions to resynchronize, for example 'S_OD' if 'RdIf' has been
**  called from 'RdFor'.  <follow>  always contain 'S_EOF', which 'Read' uses
**  to resynchronise.
**
**  If 'Match' needs to  read a  new line from  '*stdin*' or '*errin*' to get
**  the next symbol it prints the string pointed to by 'Prompt'.
*/
void Match ( unsigned long symbol, char * msg, TypSymbolSet skipto );

/****************************************************************************
**
*F  Pr( <format>, <arg1>, <arg2> )  . . . . . . . . .  print formatted output
**
**  'Pr' is the output function. The first argument is a 'printf' like format
**  string containing   up   to 2  '%'  format   fields,   specifing  how the
**  corresponding arguments are to be  printed.  The two arguments are passed
**  as  'long'  integers.   This  is possible  since every  C object  ('int',
**  'char', pointers) except 'float' or 'double', which are not used  in GAP,
**  can be converted to a 'long' without loss of information.
**
**  The function 'Pr' currently support the following '%' format  fields:
**  '%c'    the corresponding argument represents a character,  usually it is
**          its ASCII or EBCDIC code, and this character is printed.
**  '%s'    the corresponding argument is the address of  a  null  terminated
**          character string which is printed.
**  '%d'    the corresponding argument is a signed integer, which is printed.
**          Between the '%' and the 'd' an integer might be used  to  specify
**          the width of a field in which the integer is right justified.  If
**          the first character is '0' 'Pr' pads with '0' instead of <space>.
**  '%>'    increment the indentation level.
**  '%<'    decrement the indentation level.
**  '%%'    can be used to print a single '%' character. No argument is used.
**
**  You must always  cast the arguments to  '(long)' to avoid  problems  with
**  those compilers with a default integer size of 16 instead of 32 bit.  You
**  must pass 0L if you don't make use of an argument to please lint.
*/
void            Pr ( char * format, long arg1, long arg2 );


/****************************************************************************
**
*F  OpenInput( <filename> ) . . . . . . . . . .  open a file as current input
**
**  'OpenInput' opens  the file with  the name <filename>  as  current input.
**  All  subsequent input will  be taken from that  file, until it is  closed
**  again  with 'CloseInput'  or  another file  is opened  with  'OpenInput'.
**  'OpenInput'  will not  close the  current  file, i.e., if  <filename>  is
**  closed again, input will again be taken from the current input file.
**
**  'OpenInput'  returns 1 if  it   could  successfully open  <filename>  for
**  reading and 0  to indicate  failure.   'OpenInput' will fail if  the file
**  does not exist or if you do not have permissions to read it.  'OpenInput'
**  may  also fail if  you have too  many files open at once.   It  is system
**  dependent how many are  too many, but  16  files should  work everywhere.
**
**  Directely after the 'OpenInput' call the variable  'Symbol' has the value
**  'S_ILLEGAL' to indicate that no symbol has yet been  read from this file.
**  The first symbol is read by 'Read' in the first call to 'Match' call.
**
**  You can open  '*stdin*' to  read  from the standard  input file, which is
**  usually the terminal, or '*errin*' to  read from the standard error file,
**  which  is  the  terminal  even if '*stdin*'  is  redirected from  a file.
**  'OpenInput' passes those  file names  to  'SyFopen' like any other  name,
**  they are  just  a  convention between the  main  and the system  package.
**  'SyFopen' and thus 'OpenInput' will  fail to open  '*errin*' if the  file
**  'stderr'  (Unix file  descriptor  2)  is  not a  terminal,  because  of a
**  redirection say, to avoid that break loops take their input from a file.
**
**  It is not neccessary to open the initial input  file, 'InitScanner' opens
**  '*stdin*' for  that purpose.  This  file on   the other   hand can not be
**  closed by 'CloseInput'.
*/
long            OpenInput ( char * filename );


/****************************************************************************
**
*F  CloseInput()  . . . . . . . . . . . . . . . . .  close current input file
**
**  'CloseInput'  will close the  current input file.   Subsequent input will
**  again be taken from the previous input file.   'CloseInput' will return 1
**  to indicate success.
**
**  'CloseInput' will not close the initial input file '*stdin*', and returns
**  0  if such  an  attempt is made.   This is  used in  'Error'  which calls
**  'CloseInput' until it returns 0, therebye closing all open input files.
**
**  Calling 'CloseInput' if the  corresponding  'OpenInput' call failed  will
**  close the current output file, which will lead to very strange behaviour.
*/
long            CloseInput ( void );


/****************************************************************************
**
*F  OpenOutput( <filename> )  . . . . . . . . . open a file as current output
**
**  'OpenOutput' opens the file  with the name  <filename> as current output.
**  All subsequent output will go  to that file, until either   it is  closed
**  again  with 'CloseOutput' or  another  file is  opened with 'OpenOutput'.
**  The file is truncated to size 0 if it existed, otherwise it  is  created.
**  'OpenOutput' does not  close  the  current file, i.e., if  <filename>  is
**  closed again, output will go again to the current output file.
**
**  'OpenOutput'  returns  1 if it  could  successfully  open  <filename> for
**  writing and 0 to indicate failure.  'OpenOutput' will fail if  you do not
**  have  permissions to create the  file or write   to it.  'OpenOutput' may
**  also   fail if you   have  too many files   open  at once.   It is system
**  dependent how many are too many, but 16 files should work everywhere.
**
**  You can open '*stdout*'  to write  to the standard output  file, which is
**  usually the terminal, or '*errout*' to write  to the standard error file,
**  which is the terminal  even   if '*stdout*'  is  redirected to   a  file.
**  'OpenOutput' passes  those  file names to 'SyFopen'  like any other name,
**  they are just a convention between the main and the system package.
**
**  It is not neccessary to open the initial output file, 'InitScanner' opens
**  '*stdout*' for that purpose.  This  file  on the other hand   can not  be
**  closed by 'CloseOutput'.
*/
long            OpenOutput ( char * filename );


/****************************************************************************
**
*F  CloseOutput() . . . . . . . . . . . . . . . . . close current output file
**
**  'CloseOutput' will  first flush all   pending output and  then  close the
**  current  output  file.   Subsequent output will  again go to the previous
**  output file.  'CloseOutput' returns 1 to indicate success.
**
**  'CloseOutput' will  not  close the  initial output file   '*stdout*', and
**  returns 0 if such attempt is made.  This  is  used in 'Error' which calls
**  'CloseOutput' until it returns 0, thereby closing all open output files.
**
**  Calling 'CloseOutput' if the corresponding 'OpenOutput' call failed  will
**  close the current output file, which will lead to very strange behaviour.
**  On the other  hand if you  forget  to call  'CloseOutput' at the end of a
**  'PrintTo' call or an error will not yield much better results.
*/
long            CloseOutput ( void );


/****************************************************************************
**
*F  OpenAppend( <filename> )  . . open a file as current output for appending
**
**  'OpenAppend' opens the file  with the name  <filename> as current output.
**  All subsequent output will go  to that file, until either   it is  closed
**  again  with 'CloseAppend' or  another  file is  opened with 'OpenOutput'.
**  Unlike 'OpenOutput' 'OpenAppend' does not truncate the file to size 0  if
**  it exists.  Appart from that 'OpenAppend' is equal to 'OpenOutput' so its
**  description applies to 'OpenAppend' too.
*/
long            OpenAppend ( char * filename );


/****************************************************************************
**
*F  CloseAppend() . . . . . . . . . . . . . . . . . close current output file
**
**  'CloseAppend' will  first flush all   pending output and  then  close the
**  current  output  file.   Subsequent output will  again go to the previous
**  output file.  'CloseAppend' returns 1 to indicate success.  'CloseAppend'
**  is exactely equal to 'CloseOutput' so its description applies.
*/
long            CloseAppend ( void );


/****************************************************************************
**
*F  OpenLog( <filename> ) . . . . . . . . . . . . . log interaction to a file
**
**  'OpenLog'  instructs  the scanner to   echo  all  input   from  the files
**  '*stdin*' and  '*errin*'  and  all  output to  the  files '*stdout*'  and
**  '*errout*' to the file with  name <filename>.  The  file is truncated  to
**  size 0 if it existed, otherwise it is created.
**
**  'OpenLog' returns 1 if it could  successfully open <filename> for writing
**  and 0  to indicate failure.   'OpenLog' will  fail if  you do  not   have
**  permissions  to create the file or   write to  it.  'OpenOutput' may also
**  fail if you have too many files open at once.  It is system dependent how
**  many   are too   many, but  16   files should  work everywhere.   Finally
**  'OpenLog' will fail if there is already a current logfile.
*/
long            OpenLog ( char * filename );


/****************************************************************************
**
*F  CloseLog()  . . . . . . . . . . . . . . . . . . close the current logfile
**
**  'CloseLog' closes the current logfile again, so that input from '*stdin*'
**  and '*errin*' and output to '*stdout*' and '*errout*' will no  longer  be
**  echoed to a file.  'CloseLog' will return 1 to indicate success.
**
**  'CloseLog' will fail if there is no logfile active and will return  0  in
**  this case.
*/
long            CloseLog ( void );


/****************************************************************************
**
*F  OpenInputLog( <filename> )	. . . . . . . . . . . . . log input to a file
**
**  'OpenInputLog'  instructs the  scanner  to echo  all input from the files
**  '*stdin*' and  '*errin*' to the file  with  name <filename>.  The file is
**  truncated to size 0 if it existed, otherwise it is created.
**
**  'OpenInputLog' returns 1  if it  could successfully open  <filename>  for
**  writing  and  0 to indicate failure.  'OpenInputLog' will fail  if you do
**  not have  permissions to create the file  or write to it.  'OpenInputLog'
**  may also fail  if you  have  too many  files open  at once.  It is system
**  dependent  how many are too many,  but 16 files  should work  everywhere.
**  Finally 'OpenInputLog' will fail if there is already a current logfile.
*/
long            OpenInputLog ( char* );



/****************************************************************************
**
*F  CloseInputLog() . . . . . . . . . . . . . . . . close the current logfile
**
**  'CloseInputLog'  closes  the current  logfile again,  so  that input from
**  '*stdin*'  and   '*errin*'  will  no  longer   be  echoed   to  a   file.
**  'CloseInputLog' will return 1 to indicate success.
**
**  'CloseInputLog' will fail if there is no logfile active and will return 0
**  in this case.
*/
long            CloseInputLog ( void );


/****************************************************************************
**
*F  OpenTest( <filename> )  . . . . . . . .  open an input file for test mode
**
**  'OpenTest'  opens the file with the  name <filename> as current input for
**  test mode.  All subsequent input will  be taken  from that file, until it
**  is closed   again with  'CloseTest'   or another  file is   opened   with
**  'OpenInput'.   'OpenTest' will  not  close the   current file,  i.e.,  if
**  <filename> is  closed again, input will be  taken again from  the current
**  input file.
**
**  Test mode works as follows.  If the scanner  is about to  print a line to
**  the  current output file  (or to be more precise  to the output file that
**  was current when  'OpenTest' was called) this line  is  compared with the
**  next line from the test input  file, i.e., the  one opened by 'OpenTest'.
**  If this line starts with '#>' and the rest of it  matches the output line
**  the output line is  not printed and the input  comment line is discarded.
**  Otherwise the  scanner prints the output line   and does not  discard the
**  input line.
**
**  On the other hand if an input line is encountered on  the test input that
**  starts with '#>' the scanner assumes that this is an expected output line
**  that did not appear and echoes this line to the current output file.
**
**  The upshot is that  you can write  test files that consist of alternating
**  input and,  as  '#>' test  comment  lines the  expected  output.   If GAP
**  behaves normal and produces the expected  output then nothing is printed.
**  But if something  goes wrong you see  what actually was printed  and what
**  was expected instead.
**
**  As a convention GAP test files should end with a  print  statement  like:
**
**    Print("prime   3.002   06-Jul-90 ",417000000/Runtime()," GAPstones\n");
**
**  without a matching '#>' comment line.  This tells the user that the  test
**  file completed and also how much time it took.  The  constant  should  be
**  such that a VAX 11/780 gets roughly 1000 GAPstones.
**
**  'OpenTest' returns 1 if it could successfully open <filename> for reading
**  and  0 to indicate failure.  'OpenTest'  will fail if   the file does not
**  exist or if you have no permissions to read it.  'OpenTest' may also fail
**  if you have too many files open at once.  It is system dependent how many
**  are too may, but 16 files shoule work everywhere.
**
**  Directely after the 'OpenTest'  call the variable  'Symbol' has the value
**  'S_ILLEGAL' to indicate that no symbol has yet been  read from this file.
**  The first symbol is read by 'Read' in the first call to 'Match' call.
*/
long            OpenTest ( char * filename );


/****************************************************************************
**
*F  CloseTest() . . . . . . . . . . . . . . . . . . close the test input file
**
**  'CloseTest'  closes the  current test  input  file and ends  test   mode.
**  Subsequent  input   will again be taken   from  the previous  input file.
**  Output will no longer be compared with  comment lines from the test input
**  file.  'CloseTest' will return 1 to indicate success.
**
**  'CloseTest' will not close a non test input file and returns 0 if such an
**  attempt is made.
*/
long            CloseTest ( void );


/****************************************************************************
**
*F  InitScanner() . . . . . . . . . . . . . .  initialize the scanner package
**
**  'InitScanner' initializes  the  scanner  package.  This  justs  sets  the
**  current input file to '*stdin*' and current output  file  to  '*stdout*'.
*/
void            InitScanner ( void );



