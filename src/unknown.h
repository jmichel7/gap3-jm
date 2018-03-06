/****************************************************************************
**
*A  unknown.h                   GAP source                   Martin Schoenert
**
*A  @(#)$Id: unknown.h,v 1.1.1.1 1996/12/11 12:44:01 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This  file  defines  the  arithmetic  for unknown  values,  unknowns  for
**  short.  Unknowns  are written as 'Unknown(<n>)'  where  <n> is an integer
**  that distuingishes different unknowns.  Every unknown stands for a fixed,
**  well defined, but unknown  scalar value,  i.e., an  unknown  integer,  an
**  unknown rational, or an unknown cyclotomic.
**
**  Being unknown is a contagious property.  That is  to say  that the result
**  of  a scalar operation involving an  unknown is   also  unknown, with the
**  exception of multiplication by 0,  which  is  0.  Every scalar  operation
**  involving an  unknown operand is  a  new  unknown, with  the exception of
**  addition of 0 or multiplication by 1, which is the old unknown.
**
**  Note that infinity is not regarded as a well defined scalar value.   Thus
**  an unknown never stands for infinity.  Therefor division by 0 still gives
**  an  error, not an unknown.  Also  division by an  unknown gives an error,
**  because the unknown could stand for 0.
**
*H  $Log: unknown.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:44:01  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.2  1992/03/19  18:56:32  martin
*H  changed unknowns, they can no longer stand for finite field elements
*H
*H  Revision 3.1  1991/04/30  16:12:55  martin
*H  initial revision under RCS
*H
*H  Revision 3.0  1990/10/09  12:00:00  martin
*H  added unknown package
*H
*/


/****************************************************************************
**
*F  EvUnknown( <hdUnd> )  . . . . . . . . . . . . . . . . evaluate an unknown
**
**  'EvUnknown' returns the value of the unknown <hdUnd>.  Since unknowns are
**  constants and thus selfevaluating this simply returns <hdUnd>.
*/
extern  TypHandle       EvUnknown ( TypHandle hdUnk );


/****************************************************************************
**
*F  SumUnknown( <hdL>, <hdR> )  . . . . . . . . . . . . . sum of two unknowns
**
**  'SumUnknown' returns  the  sum  of  the  two  unknowns <hdL>  and  <hdR>.
**  Either operand may also be a known scalar value.
**
**  Is called from the 'Sum' binop, so both operands are already evaluated.
*/
extern  TypHandle       SumUnknown ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  DiffUnknown( <hdL>, <hdR> ) . . . . . . . . .  difference of two unknowns
**
**  'DiffUnknown' returns the difference of the two unknowns <hdL> and <hdR>.
**  Either operand may also be a known scalar value.
**
**  Is called from the 'Diff' binop, so both operands are already evaluated.
*/
extern  TypHandle       DiffUnknown ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  ProdUnknown( <hdL>, <hdR> ) . . . . . . . . . . . product of two unknowns
**
**  'ProdUnknown' returns the product of the two  unknowns  <hd>  and  <hdR>.
**  Either operand may also be a known scalar value.
**
**  Is called from the 'Prod' binop, so both operands are already evaluated.
*/
extern  TypHandle       ProdUnknown ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  QuoUnknown( <hdL>, <hdR> )  . . . . . . . . . .  quotient of two unknowns
**
**  'QuoUnknown' returns the quotient of the unknown  <hdL>  and  the  scalar
**  <hdR>.  <hdR> must not be zero, and must not be an unknown,  because  the
**  unknown could stand for zero.
**
**  Is called from the 'Quo' binop, so both operands are already evaluated.
*/
extern  TypHandle       QuoUnknown ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  PowUnknown( <hdL>, <hdR> )  . . . . . . . . . . . . . power of an unknown
**
**  'PowUnknown' returns the unknown <hdL> raised to the integer power <hdR>.
**  If <hdR> is 0, the result is the integer 1.  If <hdR> must  not  be  less
**  than 0, because <hdL> could stand for 0.
**
**  Is called from the 'Pow' binop, so both operands are already evaluted.
*/
extern  TypHandle       PowUnknown ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  EqUnknown( <hdL>, <hdR> ) . . . . . .  . . test if two unknowns are equal
**
**  'EqUnknown' returns 'true' if the two unknowns <hdL> and <hdR>  are equal
**  and 'false' otherwise.
**
**  Note that 'EqUnknown' assumes that two unknowns with  different  <n>  are
**  different.  I dont like this at all.
**
**  Is called from 'EvEq' binop, so both operands are already evaluated.
*/
extern  TypHandle       EqUnknown ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  LtUnknown( <hdL>, <hdR> ) . . .  test if one unknown is less than another
**
**  'LtUnknown' returns 'true' if the unknown <hdL> is less than the  unknown
**  <hdR>  are equal and 'false' otherwise.
**
**  Note that 'LtUnknown' assumes that two unknowns with  different  <n>  are
**  different.  I dont like this at all.
**
**  Is called from 'EvLt' binop, so both operands are already evaluated.
*/
extern  TypHandle       LtUnknown ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  PrUnknown( <hdUnk> )  . . . . . . . . . . . . . . . . .  print an unknown
**
**  'PrUnknown' prints the unknown <hdUnk> in the form 'Unknown(<n>)'.
*/
extern  void            PrUnknown ( TypHandle hdUnk );


/****************************************************************************
**
*F  FunUnknown( <hdCall> )  . . . . . . . . . . . . . . . . create an unknown
**
**  'FunUnknown' implements the internal function 'Unknown'.
**
**  'Unknown()'\\
**  'Unknown(<n>)'
**
**  In the first form 'Unknown' returns a new unknown 'Unknown(<n>)'  with  a
**  <n> that was not previously used.
**
**  In the second form 'Unknown' returns the unknown 'Unknown(<n>)'.
*/
extern  TypHandle       FunUnknown ( TypHandle hdCall );


/****************************************************************************
**
*F  FunIsUnknown( <hdCall> )  . . . . . . . .  test if an object is a unknown
**
**  'FunIsUnknown' implements the internal function 'IsUnknown'.
**
**  'IsUnknown( <obj> )'
**
**  'IsUnknown' returns 'true' if the object <obj> is an unknown and  'false'
**  otherwise.  Will cause an error if <obj> is an unbound variable.
*/
extern  TypHandle       FunIsUnknown ( TypHandle hdCall );


/****************************************************************************
**
*F  InitUnknown() . . . . . . . . . . . . . .  initialize the unknown package
**
**  'InitUnknown' initializes the unknown package.
*/
extern void             InitUnknown ( void );



