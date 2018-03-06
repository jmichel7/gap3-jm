/****************************************************************************
**
*A  rational.h                  GAP source                   Martin Schoenert
**
*A  @(#)$Id: rational.h,v 1.1.1.1 1996/12/11 12:43:58 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file declares  the  functions  for  the  artithmetic  of  rationals.
**
**  A  rational is represented  as a pair of  two integers.  The first is the
**  numerator and the second  is  the  denominator.  This representation   is
**  always reduced, i.e., numerator and denominator  are relative prime.  The
**  denominator is always  positive and  greater than 1.    If it were 1  the
**  rational would be an integer and would be represented as  integer.  Since
**  the denominator is always positive the numerator carries the sign of  the
**  rational.
**
**  It  is very easy to  see  that for every  rational   there is one  unique
**  reduced representation.  Because of   this comparisons of  rationals  are
**  quite easy,  we just compare  numerator  and denominator.  Also numerator
**  and denominator are as small as possible,  reducing the effort to compute
**  with them.   Of course  computing  the reduced  representation comes at a
**  cost.   After every arithmetic operation we  have to compute the greatest
**  common divisor of numerator and denominator, and divide them by the gcd.
**
**  In  the following   descriptions when  we talk  about  rationals  we will
**  usually include integers.  For example the function to compute the sum of
**  two rationals also handles the case that one operand is an integer.  This
**  is done by writing an integer $i$ as rational  $i / 1$ and then computing
**  with this rational.
**
*H  $Log: rational.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:58  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.1  1991/04/30  16:12:39  martin
*H  initial revision under RCS
*H
*/


/****************************************************************************
**
*F  EvRat( <hdRat> )  . . . . . . . . . . . . . . . . . . evaluate a rational
**
**  'EvRat' returns the value of the rational <hdRat>.  Because rationals are
**  constants and thus selfevaluating this just returns <hdRat>.
*/
TypHandle       EvRat ( TypHandle hdRat );


/****************************************************************************
**
*F  SumRat( <hdL>, <hdR> )  . . . . . . . . . . . . . .  sum of two rationals
**
**  'SumRat'  returns the   sum of two  rationals  <hdL>  and <hdR>.   Either
**  operand may also be an integer.  The sum is reduced.
**
**  Is called from the 'Sum' binop, so both operands are already evaluated.
*/
TypHandle       SumRat ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  DiffRat( <hdL>, <hdR> ) . . . . . . . . . . . difference of two rationals
**
**  'DiffRat' returns the  difference  of  two  rationals  <hdL>  and  <hdR>.
**  Either operand may also be an integer.  The difference is reduced.
**
**  Is called from the 'Diff' binop, so both operands are already evaluated.
*/
TypHandle       DiffRat ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  ProdRat( <hdL>, <hdR> ) . . . . . . . . . . . .  product of two rationals
**
**  'ProdRat' returns the  product of two rationals <hdL> and  <hdR>.  Either
**  operand may also be an integer.  The product is reduced.
**
**  Is called from the 'Prod' binop, so both operands are already evaluated.
*/
TypHandle       ProdRat ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  QuoRat( <hdL>, <hdR> )  . . . . . . . . . . . . quotient of two rationals
**
**  'QuoRat'  returns the quotient of two rationals <hdL> and  <hdR>.  Either
**  operand may also be an integer.  The quotient is reduced.
**
**  Is called from the 'Quo' binop, so both operands are already evaluated.
*/
TypHandle       QuoRat ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  ModRat( <hdL>, <hdL> )  . . . . . . . . remainder of rational mod integer
**
**  'ModRat' returns the remainder  of the rational  <hdL> modulo the integer
**  <hdR>.  The remainder is always an integer.
**
**  '<r>  / <s> mod  <n>' yields  the remainder of   the rational '<r> / <s>'
**  modulo the integer '<n>'.
**
**  The  modular  remainder of  $r  / s$  mod $n$  is defined  as  a $l$ from
**  $0..n-1$ such that $r = l s$ mod $n$.  As a special  case $1 / s$ mod $n$
**  is the modular inverse of $s$ modulo $n$.
**
**  Note  that the remainder will  not exist if $s$  is not relative prime to
**  $n$.  However note that $4 / 6$  mod $32$ does  exist (and is $22$), even
**  though $6$ is not invertable modulo $32$, because the $2$ cancels.
**
**  Another possible  definition of $r/s$ mod $n$  would be  a rational $t/s$
**  such that $0 \<= t/s \< n$ and $r/s - t/s$ is a multiple of $n$.  This is
**  rarely needed while computing modular inverses is very useful.
**
**  Is called from the 'Mod' binop, so both operands are already evaluated.
*/
TypHandle       ModRat ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  PowRat( <hdL>, <hdR> )  . . . . . .  raise a rational to an integer power
**
**  'PowRat' raises the rational <hdL> to the  power  given  by  the  integer
**  <hdR>.  The power is reduced.
**
**  Is called from the 'Pow' binop, so both operands are already evaluated.
*/
TypHandle       PowRat ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  EqRat( <hdL>, <hdR> ) . . . . . . . . . . . . . . test if <ratL> = <ratR>
**
**  'EqRat' returns 'true' if the two rationals <ratL> and <ratR>  are  equal
**  and 'false' otherwise.
**
**  Is called from 'EvEq' binop, so both operands are already evaluated.
*/
TypHandle       EqRat ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  LtRat( <hdL>, <hdR> ) . . . . . . . . . . . . . . test if <ratL> < <ratR>
**
**  'LtRat' returns 'true'  if  the  rational  <ratL>  is  smaller  than  the
**  rational <ratR> and 'false' otherwise.  Either operand may be an integer.
**
**  Is called from 'EvLt' binop, so both operands are already evaluated.
*/
TypHandle       LtRat ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  PrRat( <hdRat> )  . . . . . . . . . . . . . . . . . . .  print a rational
**
**  'PrRat' prints a rational in the standard form:
**
**      <numerator> / <denominator>
*/
void            PrRat ( TypHandle hdRat );


/****************************************************************************
**
*F  FunIsRat( <hdCall> )  . . . . . . . . .  internal function IsRat( <obj> )
**
**  'IsRat'  returns 'true' if the object <obj> is  a  rational  and  'false'
**  otherwise.  May cause an error if <obj> is an unbound variable.
*/
TypHandle       FunIsRat ( TypHandle hdCall );


/****************************************************************************
**
*F  FunNumerator( <hdCall> )  . . . . . . . . . . internal function Numerator
**
**  'Numerator' returns the numerator of the rational argument.
*/
TypHandle       FunNumerator ( TypHandle hdCall );


/****************************************************************************
**
*F  FunDenominator( <hdCall> )  . . . . . . . . internal function Denominator
**
**  'Denominator' returns the denominator of the rational argument.
*/
TypHandle       FunDenominator ( TypHandle hdCall );


/****************************************************************************
**
*F  InitRat() . . . . . . . . . . . . . . . . initialize the rational package
**
**  'InitRat' initializes the rational package.
*/
void            InitRat ( void );



