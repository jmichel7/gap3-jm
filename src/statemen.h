/****************************************************************************
**
*A  statemen.h                  GAP source                   Martin Schoenert
**
*A  @(#)$Id: statemen.h,v 1.1.1.1 1996/12/11 12:43:59 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This module defines the functions for executing the  various  statements.
**  Assignments are dealed with in 'eval.c' and functions  are  in  'func.c'.
**
*H  $Log: statemen.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:59  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.1  1991/04/30  16:12:48  martin
*H  initial revision under RCS
*H
*/


/****************************************************************************
**
*V  StrStat . . . . . . .  beginning text of the currently evluated statement
*V  HdStat  . . . . . . . . . . . handle of the currently evaluated statement
**
**  'StrStat' is  the  beginnnig text of   the currently evaluated statement.
**  'HdStat' is the handle of the rest of  the currently evaluated statement.
**  If an error occurs 'Error' prints 'StrStat' and 'HdStat' to give the user
**  an idea where the error occured.
*/
extern  char *          StrStat;
extern  TypHandle       HdStat;


/****************************************************************************
**
*F  EvStatseq( <hdSSeq> ) . . . . . . . . . . .  execute a statement sequence
**
**  'EvStatseq' executes the statement sequence <hdSSeq>.
**
**  This is   done by executing  the   <statement> one  after another.   If a
**  'return <expr>;' is executed inside the  statement sequence the execution
**  of the statement sequence  terminates and <expr>  is returned.  Otherwise
**  'HdVoid' is returned after execution of the last statement.
*/
TypHandle       EvStatseq ( TypHandle hdSSeq );


/****************************************************************************
**
*F  EvIf( <hdIf> )  . . . . . . . . . . . . . . . . . execute an if statement
**
**  'EvIf' executes the 'if' statement <hdIf>.
**
**  This  is done  by   evaluating the <conditions>  until  one  evaluates to
**  'true'.   Then the    corresponding <statements>  are  executed.   If  no
**  <condition> evaluates to 'true' the 'else'  <statements> are executed, if
**  present.  If  a 'return  <expr>;' statement is   executed inside the 'if'
**  statement the execution of the 'if' statement is terminated and <expr> is
**  returned, otherwise 'HdVoid' is returned.
*/
TypHandle       EvIf ( TypHandle hdIf );


/****************************************************************************
**
*F  EvFor( <hdFor> )  . . . . . . . . . . . . . . . . execute a for statement
**
**  'EvFor' executes the 'for' loop <hdFor>.
**
**  This is  done by evaluating   <list> and executing the <statements>  with
**  <variable> bound to the  first element  of the  list, then executing  the
**  <statements> with <variable> bound to the next element of the list and so
**  on.  Unbound entries in the list are skipped.  If  new elements are added
**  to the end of  <list> during a loop  iteration they will be iterated over
**  too.  If   a  'return <expr>;'  is executed  inside   the 'for'  loop the
**  execution of the 'for' loop terminates and <expr> is returned.  Otherwise
**  'HdVoid' is returned after execution of the last statement.
*/
TypHandle       EvFor ( TypHandle hdFor );


/****************************************************************************
**
*F  EvWhile( <hdWhile> )  . . . . . . . . . . . . . execute a while statement
**
**  'EvWhile' executes the 'while' loop <hdWhile>.
**
**  This is   done by   executing  the  <statements> while  the   <condition>
**  evaluates  to  'true'.  If   a 'return  <expr>;' is executed   inside the
**  'while' loop the  execution of the  'while' loop terminates and <expr> is
**  returned.  Otherwise  'HdVoid' is returned  after  execution of the  last
**  statement.
*/
TypHandle       EvWhile ( TypHandle hdWhile );


/****************************************************************************
**
*F  EvRepeat( <hdRep> ) . . . . . . . . . . . . . . . . execute a repeat loop
**
**  'EvRepeat' executes the 'repeat until' loop <hdRep>.
**
**  This   is done by executing    the  <statements>  until the   <condition>
**  evaluates   to  'true'.  If a   'return <expr>;'  is executed  inside the
**  'repeat' loop the execution of the 'repeat' loop terminates and <expr> is
**  returned.    Otherwise 'HdVoid' is returned after   execution of the last
**  statement.
*/
TypHandle       EvRepeat ( TypHandle hdRep );


/****************************************************************************
**
*F  PrStatseq( <hdSSeq> ) . . . . . . . . . . . .  print a statement sequence
**
**  'PrStatseq' prints the statement sequence <hdSSeq>.
**
**  A linebreak is forced after each <statement> except the last one.
*/
void            PrStatseq ( TypHandle hdSSeq );


/****************************************************************************
**
*F  PrIf( <hdIf> )  . . . . . . . . . . . . . . . . . . print an if statement
**
**  'PrIf' prints the 'if' statement <hdIf>.
**
**  A Linebreak is forced after the 'then'  and  <statements>.  If  necessary
**  it is preferred immediately before the 'then'.
*/
void            PrIf ( TypHandle hdIf );


/****************************************************************************
**
*F  PrFor( <hdFor> )  . . . . . . . . . . . . . . . . . . .  print a for loop
**
**  'PrFor' prints the 'for' loop <hdFor>.
**
**  A linebreak is forced after the 'do' and the <statements>.  If  necesarry
**  it is preferred immediately before the 'in'.
*/
void            PrFor ( TypHandle hdFor );


/****************************************************************************
**
*F  PrWhile( <hdWhile> )  . . . . . . . . . . . . . . . .  print a while loop
**
**  'PrWhile' prints the 'while' loop <hdWhile>.
**
**  A linebreak is forced after the 'do' and the <statements>.  If  necessary
**  it is preferred immediately before the 'do'.
*/
void            PrWhile ( TypHandle hdWhile );


/****************************************************************************
**
*F  PrRepeat( <hdRep> ) . . . . . . . . . . . . . . . . . print a repeat loop
**
**  'PrRepeat' prints the 'repeat until' loop <hdRep>.
**
**  A linebreak is forced after the 'repeat' and the <statements>.
*/
void            PrRepeat ( TypHandle hdRep );


/****************************************************************************
**
*F  InitStat()  . . . . . . . . . . . . . . . initialize the statement module
**
**  Is called from 'InitEval' to initialize the statement evaluation  module.
*/
void            InitStat ( void );




