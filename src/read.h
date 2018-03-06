/****************************************************************************
**
*A  read.h                      GAP source                   Martin Schoenert
**
*A  @(#)$Id: read.h,v 1.1.1.1 1996/12/11 12:43:58 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This module declares the functions to read  expressions  and  statements.
**
**      <Ident>         :=  a|b|..|z|A|B|..|Z { a|b|..|z|A|B|..|Z|0|..|9|_ }
**
**      <Var>           :=  <Ident>
**                      |   <Var> '.' <Ident>
**                      |   <Var> '[' <Expr> ']'
**                      |   <Var> '{' <Expr> '}'
**                      |   <Var> '(' [ <Expr> { ',' <Expr> } ] ')'
**
**      <List>          :=  '[' [ <Expr> ] {',' [ <Expr> ] } ']'
**                      |   '[' <Expr> '..' <Expr> ']'
**
**      <Record>        :=  'rec( [ <Ident>:=<Expr> {, <Ident>:=<Expr> } ] )'
**
**      <Permutation>   :=  ( <Expr> {, <Expr>} ) { ( <Expr> {, <Expr>} ) }
**
**      <Function>      :=  'function (' [ <Ident> {',' <Ident>} ] ')'
**                              [ 'local'  <Ident> {',' <Ident>} ';' ]
**                              <Statments>
**                          'end'
**
**      <String>        :=  " { <any character> } "
**
**      <Int>           :=  0|1|..|9 { 0|1|..|9 }
**
**      <Atom>          :=  <Int>
**                      |   <Var>
**                      |   '(' <Expr> ')'
**                      |   <Permutation>
**                      |   <String>
**                      |   <Function>
**                      |   <List>
**                      |   <Record>
**
**      <Factor>        :=  {'+'|'-'} <Atom> [ '^' {'+'|'-'} <Atom> ]
**
**      <Term>          :=  <Factor> { '*'|'/'|'mod' <Factor> }
**
**      <Arith>         :=  <Term> { '+'|'-' <Term> }
**
**      <Rel>           :=  { 'not' } <Arith> { '=|<>|<|>|<=|>=|in' <Arith> }
**
**      <And>           :=  <Rel> { 'and' <Rel> }
**
**      <Log>           :=  <And> { 'or' <Rel> }
**
**      <Expr>          :=  <Log>
**                      |   <Var> [ '->' <Log> ]
**
**      <Statment>      :=  <Expr>
**                      |   <Var> ':=' <Expr>
**                      |   'if'   <Expr>  'then' <Statments>
**                        { 'elif' <Expr>  'then' <Statments> }
**                        [ 'else'                <Statments> ] 'fi'
**                      |   'for' <Var>  'in' <Expr>  'do' <Statments>  'od'
**                      |   'while' <Expr>  'do' <Statments>  'od'
**                      |   'repeat' <Statments>  'until' <Expr>
**                      |   'return' [ <Expr> ]
**                      |   'quit'
**
**      <Statments>     :=  { <Statment> ; }
**                      |   ;
**
*H  $Log: read.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:58  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.2  1992/12/08  11:50:26  martin
*H  added '<list>{<positions>}'
*H
*H  Revision 3.1  1991/04/30  16:12:41  martin
*H  initial revision under RCS
*H
*H  Revision 3.0  1990/12/27  12:00:00  martin
*H  changed the precedence of the logical operators
*H
*/


/****************************************************************************
**
*F  ReadIt()  . . . . . . . . . . . . . . . . . read a statement interactivly
**
**  'ReadIt' reads a single statement, returning the handle to the  new  bag.
**  This is the only reading function that doesn't expect the first symbol of
**  its input already read and wont read the first symbol of the  next input.
**  This is the main interface function for the  various  ReadEvalPrintLoops.
**
**  It has this funny name, because 'Read' would give name clash with  'read'
**  from the C library on the stupid VAX, which turns all names to uppercase.
*/
TypHandle       ReadIt ( void );




