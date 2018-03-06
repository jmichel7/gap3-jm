/****************************************************************************
**
*A  cyclotom.h                  GAP source                   Martin Schoenert
**
*A  @(#)$Id: cyclotom.h,v 1.1.1.1 1996/12/11 12:43:52 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file implements the arithmetic for elements from  cyclotomic  fields
**  $Q(e^{{2 \pi i}/n}) = Q(e_n)$,  which  we  call  cyclotomics  for  short.
**
**  For a full description of cyclotimics  see chapter  "Cyclotomics"  in the
**  \GAP\ Manual.  Read   also  section "More  about  Cyclotomics" about  the
**  choosen base.
**
*H  $Log: cyclotom.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:52  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.6  1991/04/30  16:12:14  martin
*H  initial revision under RCS
*H
*H  Revision 3.5  1990/12/07  12:00:00  martin
*H  removed 'FunTraceCyc'
*H
*H  Revision 3.4  1990/09/30  12:00:00  martin
*H  added 'TraceCyc'
*H
*H  Revision 3.3  1990/09/29  12:00:00  martin
*H  added 'CoeffsCyc'
*H
*H  Revision 3.2  1990/09/28  13:00:00  martin
*H  added 'NofCyc' and 'IsCycInt'
*H
*H  Revision 3.1  1990/09/28  12:00:00  martin
*H  added 'GaloisCyc'
*H
*H  Revision 3.0  1990/09/27  12:00:00  martin
*H  added new cyclotomics package
*H
*/


/****************************************************************************
**
*F  EvCyc( <hdCyc> )  . . . . . . . . . . . . . . . . . evaluate a cyclotomic
**
**  'EvCyc'   returns   the  value   of the    cyclotomic  <hdCyc>.   Because
**  cyclomtomics are constants  and  thus  selfevaluating  this just  returns
**  <hdCyc>.
*/
extern  TypHandle       EvCyc ( TypHandle hdCyc );


/****************************************************************************
**
*F  SumCyc( <hdL>, <hdR> )  . . . . . . . . . . . . .  sum of two cyclotomics
**
**  'SumCyc' returns  the  sum  of  the  two  cyclotomics  <hdL>  and  <hdR>.
**  Either operand may also be an integer or a rational.
**
**  Is called from the 'Sum' binop, so both operands are already evaluated.
*/
extern  TypHandle       SumCyc ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  DiffCyc( <hdL>, <hdR> ) . . . . . . . . . . difference of two cyclotomics
**
**  'DiffCyc' returns the difference of the two cyclotomic <hdL>  and  <hdR>.
**  Either operand may also be an integer or a rational.
**
**  Is called from the 'Diff' binop, so both operands are already evaluated.
*/
TypHandle       DiffCyc ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  ProdCycI( <hdL>, <hdR> )  . . . .  product of a cyclotomic and an integer
**
**  'ProdCycI' returns the product of a cyclotomic and a integer or rational.
**  Which operand is the cyclotomic and wich the integer does not matter.
**
**  Called from the 'Prod' binop, so both operands are already evaluated.
**
**  This is a special case, because if the integer is not 0, the product will
**  automatically be base reduced.  So we dont need to  call  'ConvertToBase'
**  or 'Reduce' and directly write into a result bag.
*/
extern  TypHandle       ProdCycI ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  ProdCyc( <hdL>, <hdR> ) . . . . . . . . . . .  product of two cyclotomics
**
**  'ProdCyc' returns the product of the two  cyclotomics  <hdL>  and  <hdR>.
**  Either operand may also be an integer or a rational.
**
**  Is called from the 'Prod' binop, so both operands are already evaluated.
*/
extern  TypHandle       ProdCyc ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  QuoCyc( <hdL>, <hdR> )  . . . . . . . . . . . quotient of two cyclotomics
**
**  'QuoCyc' returns the quotient of the  cyclotomic  <hdL>  divided  by  the
**  cyclotomic <hdR>.  Either operand may also be an integer or a rational.
**
**  Is called from the 'Quo' binop, so both operands are already evaluated.
**
**  'QuoCyc' computes the inverse of <hdR> by computing the  product  $p$  of
**  nontrivial galois conjugates of <hdR>.  Then  $hdR * (p / (hdR * p)) = 1$
**  so $p / (hdR * p)$ is the  inverse  of  $hdR$.  Because  the  denominator
**  $hdR*p$ is the norm of $hdR$ over the rationals it is rational so we  can
**  compute the quotient $p / (hdL * p)$ with 'ProdCycI'.
*/
extern  TypHandle       QuoCyc ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  PowCyc( <hdL>, <hdR> )  . . . . . . . . . . . . . . power of a cyclotomic
**
**  'PowCyc' returns the <hdR>th, which must be  an  integer,  power  of  the
**  cyclotomic <hdL>.  The left operand may also be an integer or a rational.
**
**  Is called from the 'Pow' binop, so both operands are already evaluated.
*/
extern  TypHandle       PowCyc ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  EqCyc( <hdL>, <hdR> ) . . . . . . . . . test if two cyclotomics are equal
**
**  'EqCyc' returns 'true' if the two cyclotomics <hdL>  and <hdR>  are equal
**  and 'false' otherwise.
**
**  'EqCyc'  is  pretty  simple because   every    cyclotomic  has a   unique
**  representation, so we just have to compare the terms.
**
**  Is called from 'EvEq' binop, so both operands are already evaluated.
*/
extern  TypHandle       EqCyc ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  LtCyc( <hdL>, <hdR> ) . . . . test if one cyclotomic is less than another
**
**  'LtCyc'  returns  'true'  if  the  cyclotomic  <hdL>  is  less  than  the
**  cyclotomic <hdR> and 'false' otherwise.
**
**  Cyclotomics are first sorted according to the order of the primitive root
**  they are written in.  That means that the rationals  are  smallest,  then
**  come cyclotomics from $Q(e_3)$ followed by cyclotomics from $Q(e_4)$ etc.
**  Cyclotomics from the same field are sorted lexicographicaly with  respect
**  to their representation in the base of this field.  That means  that  the
**  cyclotomic with smaller coefficient for the first base root  is  smaller,
**  for cyclotomics with the same first coefficient the second decides  which
**  is smaller, etc.
**
**  'LtCyc'  is  pretty  simple because   every    cyclotomic  has a   unique
**  representation, so we just have to compare the terms.
**
**  Is called from 'EvLt' binop, so both operands are already evaluated.
*/
extern  TypHandle       LtCyc ( TypHandle hdL, TypHandle hdR );


/****************************************************************************
**
*F  PrCyc( <hdCyc> )  . . . . . . . . . . . . . . . . . .  print a cyclotomic
**
**  'PrCyc' prints the cyclotomic <hdCyc> in the standard form.
*/
extern void            PrCyc ( TypHandle hdCyc );


/****************************************************************************
**
*F  FunE( <hdCall> )  . . . . . . . . . . . . . . create a new primitive root
**
**  'FunE' implements the internal function 'E'.
**
**  'E( <n> )'
**
**  'E' return a the primitive root of order <n>, which must  be  a  positive
**  integer, represented as cyclotomic.
*/
extern  TypHandle       FunE ( TypHandle hdCall );


/****************************************************************************
**
*F  FunIsCyc( <hdCall> )  . . . . . . . . .  test if an object is a cyclomtic
**
**  'FunIsCyc' implements the internal function 'IsCyc'.
**
**  'IsCyc( <obj> )'
**
**  'IsCyc' returns 'true' if the object <obj> is a  cyclotomic  and  'false'
**  otherwise.  Will cause an error if <obj> is an unbound variable.
*/
extern  TypHandle       FunIsCyc ( TypHandle hdCall );


/****************************************************************************
**
*F  FunIsCycInt( <hdCall> ) . . . .  test if an object is a cyclomtic integer
**
**  'FunIsCycInt' implements the internal function 'IsCycInt'.
**
**  'IsCycInt( <obj> )'
**
**  'IsCycInt' returns 'true' if the object <obj> is a cyclotomic integer and
**  'false' otherwise.  Will cause an error if <obj> is an unbound variable.
**
**  'IsCycInt' relies on the fact that the base is an integral base.
*/
extern  TypHandle       FunIsCycInt ( TypHandle hdCall );


/****************************************************************************
**
*F  FunNofCyc( <hdCall> ) . . . . . . . . . . . . . . . . . N of a cyclotomic
**
**  'FunNofCyc' implements the internal function 'NofCyc'.
**
**  'NofCyc( <cyc> )'
**
**  'NofCyc' returns the N of the cyclotomic <cyc>, i.e., the  order  of  the
**  roots of which <cyc> is written as a linear combination.
*/
extern  TypHandle       FunNofCyc ( TypHandle hdCall );


/****************************************************************************
**
*F  FunCoeffsCyc( <hdCall> )  . . . . . . . . .  coefficients of a cyclotomic
**
**  'FunCoeffCyc' implements the internal function 'CoeffsCyc'.
**
**  'CoeffsCyc( <cyc> )'
**
**  'CoeffsCyc' returns a list of the coefficients of the  cyclotomic  <cyc>.
**  The list has lenght <n> if <n> is the order of the primitive  root  $e_n$
**  of which <cyc> is written as a linear combination.  The <i>th element  of
**  the list is the coefficient of $e_l^{i-1}$.
*/
extern  TypHandle       FunCoeffsCyc ( TypHandle hdCall );


/****************************************************************************
**
*F  FunGaloisCyc(<hdCall>)  image of a cyclotomic under a galois automorphism
**
**  'FunGaloisCyc' implements the internal function 'GaloisCyc'.
**
**  'GaloisCyc( <cyc>, <ord> )'
**
**  'GaloisCyc' computes the image of the cyclotomic <cyc>  under  the galois
**  automorphism given by <ord>, which must be an integer.
**
**  The galois automorphism is the mapping that  takes  $e_n$  to  $e_n^ord$.
**  <ord> may be any integer, of course if it is not relative prime to  $ord$
**  the mapping will not be an automorphism, though still an endomorphism.
*/
extern  TypHandle       FunGaloisCyc ( TypHandle hdCall );


/****************************************************************************
**
*F  InitCyc() . . . . . . . . . . . . . . . initialize the cyclotomic package
**
**  'InitCyc' initializes the cyclotomic package.
*/
extern  void            InitCyc ( void );



