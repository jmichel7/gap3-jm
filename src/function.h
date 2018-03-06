/****************************************************************************
**
*A  function.h                  GAP source                   Martin Schoenert
**
*A  @(#)$Id: function.h,v 1.1.1.1 1996/12/11 12:43:53 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This package delcares the functions  that  mainly  deal  with  functions.
**
*H  $Log: function.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:53  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.1  1991/04/30  16:12:19  martin
*H  initial revision under RCS
*H
*/


/****************************************************************************
**
*V  HdExec  . . . . . . . . . . . . . . . handle of the topmost execution bag
**
**  'HdExec' is the handle of the topmost execution bag on the execution  bag
**  linked list.  For every active function there is  one  execution  bag  on
**  this list.  'HdExec' is the execution bag of the current function.
**
**  The execution bag list is the equivalent  of  the  stack  in  programming
**  languages like C or Pascal.  'HdExec' is therefor the equivalent  of  the
**  Stackpointer in those languages.
**
**  To be precise a function execute bag has the following form:
**
**      handle 0:               handle of enclosing environment.
**      handle 1:               handle 'HdTrue' if this frame is current.
**      handle 2:               handle of function definition.
**      handle 3:               handle of function call bag, debug only.
**      handle 4:               handle of calling enviroment, debug only.
**      handle 5..n+4:          handles of old arguments.
**      handle n+5..n+m+4:      handles of old local variables.
**      handle n+m+5..2*n+m+4:  handles of new arguments.
*/
extern  TypHandle       HdExec;


/****************************************************************************
**
*F  ChangeEnv( <hdEnv> )  . . . .  change the environment for a function call
**
**  'ChangeEnv' changes the environment for a function call.  A *environment*
**  is  the set of  bindings of  identifiers to variables.   GAP has  lexical
**  binding, i.e., the environment in effect  when a function is created, the
**  so  called *definition  environment*,  determines the  variable bindings.
**  Thus when a  function is called its   definition environment is made  the
**  current environment.   When the function  terminates the old environment,
**  the so called *execution environment* is made current again.
**
**  An environment is stored  as a linked list  of exec-bags.  Every exec-bag
**  contains, among other things, the changes that a certain function made to
**  the environment.  I.e., when a function is  called it introduces a set of
**  new  arguments and local variables.   If the function was already active,
**  i.e., if the call was recursive, the exec-bag remembers the old values of
**  the  variables.  If  the  function was  not  already active, the exec-bag
**  remembers the fact that the variables had no values prior to the call.
**
**  The following picture should make the  operation  of  'ChangeEnv'  clear:
**
**      <hdEnv>  -> <exec 1> -> <exec 2> -\
**                                         \
**      'HdExec' -> <exec 3> -> <exec 4> ---+-> <exec 5> -> <exec 6> ... -> 0
**
**  'HdExec' is the handle of the current environment.  <hdEnv> is the handle
**  of an environment of a function that is  just beeing called.  'ChangeEnv'
**  must now change  the environment from 'HdExec' to  <hdEnv>.  To do  so it
**  must undo the changes stored in <exec 3> and  <exec 4> and then must redo
**  the changes stored in <exec 2> and <exec 1>, in that order.
**
**  Note that functions which are defined globally  can not access non-local,
**  non-global  variables.  Therefor it makes no  difference in  which such a
**  function  is executed.  In this  case  'EvFunccall'  does not change  the
**  environment at all.  Thus instead of:
**
**      <hdEnv> -> <exec 1> -------------\
**                                        \
**      'HdExec' -> <exec 3> -> <exec 4> --+-> 0
**
**  'EvFuncall' acts as if the situation was:
**
**      <hdEnv> -> <exec 1> -\
**                            \
**      'HdExec' --------------+-> <exec 3> -> <exec 4> -> 0
*/
void            ChangeEnv ( TypHandle hdEnv );



/****************************************************************************
**
*F  EvFunccall( <hdCall> )  . . . . . . . . . . . . evaluates a function call
**
**  'EvFunccall' evaluates the function call with  the  handle  <hdCall>  and
**  returns the value returned by the function or 'HdVoid'  if  the  function
**  did not return any value at all.
**
**  The function call bag <hdCall> has the following form:
**
**      handle 0:               handle of the function definition bag.
**      handle 1.. :            handles of arguments (not yet evaluated).
**
**  'EvFunccall' first creates a new  execute  bag.  The  new  arguments  are
**  evaluated and the values are put in the execute bag.  Then the old values
**  of the arguments and local variables are saved in the  execute bag.  Then
**  'ChangeEnv' is called to copy the new values from the  execute  bag  into
**  the variables.  Now the binding is complete and the statement sequence is
**  executed.  After that 'ChangeEnv' is called  again  to  restore  the  old
**  values from the execute bag.
*/
TypHandle       EvFunccall ( TypHandle hdCall );


/****************************************************************************
**
*F  EvFunction( <hdFun> ) . . . . . . . . . . . . . . . . evaluate a function
**
**  'EvFunction' returns the value of the function <hdFun>.  Since  functions
**  are constants and thus selfevaluating it just returns <hdFun>.
*/
TypHandle       EvFunction ( TypHandle hdDef );


/****************************************************************************
**
*F  EvMakefunc( <hdFun> ) . . . . . . . . . . . . . . . . . . make a function
**
**  'EvMakefunc' makes a function, i.e., turns a  variable  function  into  a
**  constant one.  GAP has lexical binding.  This means that the binding from
**  identifiers to variables is determined by the environment that was active
**  when a function was created and not by the one active  when  the function
**  is  executed.  'ChangeEnv'  performs  the  task  of  switching  from  the
**  active execution environment to the definition environment of a function.
**  But in order to do this it needs to know the  definition  environment  of
**  a function.  'EvMakefunc' copies the function  definition  bag  and  adds
**  the handle of the current  environment  to  that  bag.  This  process  is
**  usually called closing the function and the result is called  a  closure.
**
**  To be precise, the make-function bag created by the parser has the form:
**
**      handle 0:               handle of the statement sequence.
**      handle 1..n:            handles of the arguments.
**      handle n+1..n+m:        handles of the local variables.
**      handle n+m+1:           0.
**      data 1:                 (short) number of arguments (n).
**      data 2:                 (short) number of local variables (m).
**
**  And 'EvMakefunc' makes a copy of the form:
**
**      handle 0:               handle of the statement sequence.
**      handle 1..n:            handles of the arguments.
**      handle n+1..n+m:        handles of the local variables.
**      handle n+m+1:           handle of the definition environment.
**      data 1:                 (short) number of arguments (n).
**      data 2:                 (short) number of local variables (m).
*/
TypHandle       EvMakefunc ( TypHandle hdFun );


/****************************************************************************
**
*F  EvReturn( <hdRet> ) . . . . . . . . . . . . . evaluate a return-statement
**
**  'EvReturn' executes the return-statement withthe handle <hdRet>.
**
**  This   is sort of tricky.  'EvReturn'   actually  does nothing but return
**  <hdRet>.   Any  of the  functions  that execute  statements,  for example
**  'EvIf'  look for return bags as  the result of  executing a substatement.
**  If  they encounter  such  a bag they  pass  it to their calling function.
**  Thus sooner  or later it will get  to 'EvFunccall' as result of executing
**  the   function body.   'EvFunccall'    then either  returns   'HdVoid' or
**  evaluates the expression of the 'return' and returns its value.
*/
TypHandle       EvReturn ( TypHandle hdRet );


/****************************************************************************
**
*F  FunTrace( <hdCall> )  . . . . . . . . . . . . . . internal function Trace
**
**  'FunTrace' implements the internal function 'Trace'.
**
**  'Trace( <function>... )'
**
**  'Trace' switches on  tracing  for  the  functions  passed  as  arguments.
**  Whenever such a function is called GAP prints a message of the form:
**
**      <function1>( <arg1>, <arg2>, ... )
**       <function2>()
**        ...
**       <function2> returns
**      <function1> returns <value>
**
**  Where <function1>, <function2>, <arg1>, <arg2> and <value>  are  replaced
**  by the respective values.
**
**  'Untrace' switches this off again.
*/
TypHandle       FunTrace ( TypHandle hdCall );


/****************************************************************************
**
*F  FunUntrace( <hdCall> )  . . . . . . . . . . . . internal function Untrace
**
**  'FunUntrace' implements the internal function 'Untrace'.
**
**  'Untrace( <function>... )'
**
**  'Untrace' switches of the tracing for the functions passed as  arguments.
*/
TypHandle       FunUntrace ( TypHandle hdCall );


/****************************************************************************
**
*F  FunProfile( <hdCall> )  . . . . . . . . . . . . internal function Profile
**
**  'FunProfile' implements the internal function 'Profile'.
**
**  'Profile( true )'
**  'Profile( false )'
**  'Profile()'
**
**  'Profile' controls the function profiling.
**
**  In the first form,  with  the  argument  'true', 'Profile'  switches  the
**  profiling on.  From that moment on for every function GAP  remembers  the
**  number of times  this  function  was  called,  the  time  spent  in  this
**  function without its children, i.e., the functions it  called  and  their
**  children, and the time spent in this function together with them.  If the
**  profiling was already on, 'Profile' clears the profiling information.
**
**  In the second form, with the  argument  'false', 'Profile'  switches  the
**  profiling off again.  Note that programs run faster without profiling.
**
**  In the third form, without  arguments,  'Profile'  prints  the  profiling
**  information.
*/
TypHandle       FunProfile ( TypHandle hdCall );


/****************************************************************************
**
*F  PrFuncint( <hdFun> )  . . . . . . . . . . . .  print an internal function
**
**  'PrFuncint' prints the internal function with the handle  <hdFun> in  the
**  short form: 'function (...) internal; end'.
*/
void            PrFuncint ( TypHandle hdFun );


/****************************************************************************
**
*F  PrFunction( <hdFun> ) . . . . . . . . . . . . . . . . .  print a function
**
**  'PrFunction' prints the function with the handle <hdFun>  either  in  the
**  short format:
**
**      function ( <args> ) ... end
**
**  if 'prFull' is 0, or in the long format:
**
**      function ( <args> )
**          local <locals>;
**          <statements>
**      end
**
**  otherwise.
*/
void            PrFunction ( TypHandle hdFun );


/****************************************************************************
**
*F  PrintFunction( <hdFun> )  . . . . . . . print a function in the full form
**
**  'PrintFunction' prints the function with the handle <hdFun> in  the  full
**  form, i.e., with the statement sequence.  It is called from 'Print'.
*/
void            PrintFunction ( TypHandle hdFun );


/****************************************************************************
**
*F  PrFunccall( <hdCall> )  . . . . . . . . . . . . . . print a function call
**
**  'PrFunccall' prints the function call with the  handle  <hdCall>  in  the
**  usual form:  '<function>( <args> )'.
**
**  Linebreaks are preffered after the opening  parenthesis  and  the  commas
**  between the arguments.
*/
void            PrFunccall ( TypHandle hdCall );


/****************************************************************************
**
**  PrReturn( <hdRet> ) . . . . . . . . . . . . . .  print a return statement
**
**  'PrReturn' prints the return statement with the  handle  <hdRet>  in  the
**  usual form 'return;' or 'return <expr>;'.
*/
void            PrReturn ( TypHandle hdRet );


/****************************************************************************
**
*F  InitFunc()  . . . . . . . . . . .  initialize function evaluation package
**
**  'InitFunc' initializes the function evaluation package.
*/
void            InitFunc ( void );



