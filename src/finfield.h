/****************************************************************************
**
*A  finfield.h                  GAP source                      Werner Nickel
**                                                         & Martin Schoenert
**
*H  @(#)$Id: finfield.h,v 1.2 1999/10/11 15:53:49 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file declares the functions to compute with  finite fields elements.
**
**  Finite fields  are an  important   domain in computational   group theory
**  because the classical matrix groups are defined over those finite fields.
**  In GAP we support finite fields with  up to 65536 elements, larger fields
**  can be realized as polynomial domains over smaller fields.
**
*H  $Log: finfield.h,v $
*H  Revision 1.2  1999/10/11 15:53:49  werner
*H  This is a fix of the bug which causes wrong results if the primitive
*H  element of a finite field has a high exponent and is raised to a high
*H  exponent again.  First noticed in GAP 4 but already there in GAP 3.
*H
*H  Revision 1.1.1.1  1996/12/11 12:43:53  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.4  1994/01/27  15:15:55  fceller
*H  removed 'FunIntFFE',  added 'FunIntFFE2'.
*H
*H  Revision 3.3  1993/02/04  10:51:10  martin
*H  changed the interface slightly
*H
*H  Revision 3.2  1992/04/03  16:07:03  martin
*H  fixed 'IsVector' and 'IsMatrix'
*H
*H  Revision 3.1  1991/04/30  16:12:18  martin
*H  initial revision under RCS
*H
*H  Revision 3.0  1991/01/06  12:00:00  martin
*H  improved the finite field package
*H
*/


/****************************************************************************
**
*T  TypFFE  . . . . . . . . . . . type of the values of finite field elements
**
**  'TypFFE' is the type used to store the values of finite  field  elements,
**  i.e., the logarithm of the element with respect to  the  root  plus  one.
**
**  Finite fields are restricted to contain at most 65536 elements, so we can
**  put the logarithm of any finite field element into an 'unsigned short'.
**
**  It is possible to change this to 'unsigned long'  to  allow  fields  with
**  more than than 65536 elements.  All the below macros and have been  coded
**  in such a way that they work without problems.  The exception is 'POW_FF'
**  which will only work if the product of integers of type 'TypFFE' does not
**  cause an overflow.  And of course the successor table stored for a finite
**  field will become quite large for fields with more than 65536 elements.
*/
typedef unsigned short          TypFFE;


/****************************************************************************
**
*F  SUCC_FF(<hdFF>) . . . . . . . . . . . . . successor table of finite field
**
**  'SUCC_FF' returns a pointer to the successor table  of  the  finite field
**  <hdFF>.
**
**  Note that 'SUCC_FF' is a macro, so do not call  it  with  arguments  that
**  sideeffects.
*/
#define SUCC_FF(FF)             ((TypFFE*)PTR(FF))


/****************************************************************************
**
*F  SIZE_FF(<hdFF>) . . . . . . . . . . . . . . . . .  size of a finite field
**
**  'SIZE_FF' returns the size of the finite field <hdFF>.
**
**  Note that 'SIZE_FF' is a macro, so do not call  it  with  arguments  that
**  have sideeffects.
*/
#define SIZE_FF(FF)             (*SUCC_FF(FF)+1)


/****************************************************************************
**
*F  FLD_FFE(<hdFFE>)  . . . . . . . . . . . . field of a finite field element
**
**  'FLD_FFE' returns the handle of the finite field over  which  the  finite
**  field element <hdFFE> is represented.
**
**  Note that 'FLD_FFE' is a macro, so do not call  it  with  arguments  that
**  have sideeffects.
*/
#define FLD_FFE(FFE)            (PTR(FFE)[0])


/****************************************************************************
**
*F  SET_FLD_FFE(<hdFFE>,<hdFF>) . . . set the field of a finite field element
**
**  'SET_FLD_FFE' sets the field of  the  finite  field  element  <hdFFE>  to
**  <hdFF>.
**
**  Note that 'SET_FLD_FFE' is a macro, so do not call it with arguments that
**  have sideeffects.
*/
#define SET_FLD_FFE(FFE,FF)     (FLD_FFE(FFE)=(FF))


/****************************************************************************
**
*F  VAL_FFE(<hdFFE>)  . . . . . . . . . . . . value of a finite field element
**
**  'VAL_FFE' returns the value of the finite field element <hdFFE>.  That is
**  if <hdFFE> is $0_F$,'VAL_FFE' returns 0;  if <hdFFE> is $1_F$,  'VAL_FFE'
**  returns 1; if <hdFFE> is the primitive generator $z$,  'VAL_FFE'  returns
**  2; and otherwise if <hdFFE> is $z^i$, 'VAL_FFE' returns '<i>+1'.
**
**  Note that 'VAL_FFE' is a macro, so do not call  it  with  arguments  that
**  have sideeffects.
*/
#define VAL_FFE(FFE)            (*(TypFFE*)(PTR(FFE)+1))


/****************************************************************************
**
*F  SET_VAL_FFE(<hdFFE>,<val>)  . . . set the value of a finite field element
**
**  'SET_VAL_FFE' sets the value of the finite field element <hdFFE>  to  the
**  value <val>.  Note that no checking is performed, whether <val>  lies  in
**  the allowed range.
**
**  Note that 'SET_VAL_FFE' is a macro, so do not call it with arguments that
**  have sideeffects.
*/
#define SET_VAL_FFE(FFE,VAL)    (VAL_FFE(FFE)=(VAL))


/****************************************************************************
**
*F  SUM_FF(<a>,<b>,<f>) . . . . . . . . . . . . .  sum of finite field values
**
**  'SUM_FF' returns the sum of the two  finite  field  values  <a>  and  <b>
**  from the finite field pointed to by the pointer <f>.
**
**  Note that  'SUM_FF'  may only be used if  the operands are represented in
**  the same finite field.  If you want to add two elements where one lies in
**  a subfield of the other use 'SumFFE'.
**
**  Use 'SUM_FF' only with arguments that are  variables  or  array elements,
**  because it is a macro and arguments with sideeffects will behave strange,
**  and because it is a complex macro so most C compilers will  be  upset  by
**  complex arguments.  In particular, do not use 'SUM_FF(a,NEG_FF(b,f),f)'.
**
**  If either operand is 0, the sum is just the other operand.
**  If $a <= b$ we have
**  $a + b ~ z^{a-1}+z^{b-1} = z^{a-1} * (z^{(b-1)-(a-1)}+1) ~ a * f[b-a+1]$,
**  otherwise we have
**  $a + b ~ z^{b-1}+z^{a-1} = z^{b-1} * (z^{(a-1)-(b-1)}+1) ~ b * f[a-b+1]$.
*/
#define SUM_FF(a,b,f)   ( (a)==0 || (b)==0 ? (a)+(b) :\
                          ( (a)<=(b) ? PROD_FF(a,(f)[(b)-(a)+1],f) :\
                                       PROD_FF(b,(f)[(a)-(b)+1],f) ) )


/****************************************************************************
**
*F  NEG_FF(<a>,<f>) . . . . . . . . . . . . .  negative of finite field value
**
**  'NEG_FF' returns the negative of the finite  field  value  <a>  from  the
**  finite field pointed to by the pointer <f>.
**
**  Use 'NEG_FF' only with arguments that are  variables  or  array elements,
**  because it is a macro and arguments with sideeffects will behave strange,
**  and because it is a complex macro so most C compilers will  be  upset  by
**  complex arguments.  In particular, do not use 'NEG_FF(PROD_FF(a,b,f),f)'.
**
**  If the characteristic is 2, every element is its  own  additive  inverse.
**  Otherwise note that $z^{o-1} = 1 = -1^2$ so $z^{(o-1)/2} = 1^{1/2} = -1$.
**  If $a <= (o-1)/2$ we have
**  $-a ~ -1 * z^{a-1} = z^{(o-1)/2} * z^{a-1} = z^{a+(o-1)/2-1} ~ a+(o-1)/2$
**  otherwise we have
**  $-a ~ -1 * z^{a-1} = z^{a+(o-1)/2-1} = z^{a+(o-1)/2-1-(o-1)} ~ a-(o-1)/2$
*/
#define NEG_FF(a,f)     ( (a)==0 ? 0 :\
                          ( *(f)%2==1 ? a :\
                            ( (a)<=*(f)/2 ? (a)+*(f)/2 : (a)-*(f)/2 ) ) )


/****************************************************************************
**
*F  PROD_FF(<a>,<b>,<f>)  . . . . . . . . . . . product of finite field value
**
**  'PROD_FF' returns the product of the two finite field values <a> and  <b>
**  from the finite field pointed to by the pointer <f>.
**
**  Note that  'PROD_FF'  may only be used if the operands are represented in
**  the same finite field.  If you want to multiply two  elements  where  one
**  lies in a subfield of the other use 'ProdFFE'.
**
**  Use 'PROD_FF' only with arguments that are variables  or  array elements,
**  because it is a macro and arguments with sideeffects will behave strange,
**  and because it is a complex macro so most C compilers will  be  upset  by
**  complex arguments.  In particular, do not use 'NEG_FF(PROD_FF(a,b,f),f)'.
**
**  If one of the values is 0 the product is 0.
**  If $a+b <= o$ we have $a * b ~ z^{a-1} * z^{b-1} = z^{(a+b-1)-1} ~ a+b-1$
**  otherwise   we   have $a * b ~ z^{(a+b-2)-(o-1)} = z^{(a+b-o)-1} ~ a+b-o$
*/
#define PROD_FF(a,b,f)  ( (a)==0 || (b)==0 ? 0 :\
                          ( (a)-1<=*(f)-(b) ? (a)-1+(b) : (a)-1-(*(f)-(b)) ))


/****************************************************************************
**
*F  QUO_FF(<a>,<b>,<f>) . . . . . . . . . . . quotient of finite field values
**
**  'QUO_FF' returns the quotient of the two finite field values <a> and  <b>
**  from the finite field pointed to by the pointer <f>.
**
**  Note that  'QUO_FF'  may only be used if the operands are represented  in
**  the same finite field.  If you want to  divide  two  elements  where  one
**  lies in a subfield of the other use 'QuoFFE'.
**
**  Use 'QUO_FF' only with arguments that are  variables  or  array elements,
**  because it is a macro and arguments with sideeffects will behave strange,
**  and because it is a complex macro so most C compilers will  be  upset  by
**  complex arguments.  In particular, do not use 'NEG_FF(PROD_FF(a,b,f),f)'.
**
**  A division by 0 is an error,  and dividing 0 by a nonzero value gives  0.
**  If $0 <= a-b$ we have  $a / b ~ z^{a-1} / z^{b-1} = z^{a-b+1-1} ~ a-b+1$,
**  otherwise   we   have  $a / b ~ z^{a-b+1-1}  =  z^{a-b+(o-1)}   ~ a-b+o$.
*/
#define QUO_FF(a,b,f)   ( (a)==0 ? 0 :\
                          ( (b)<=(a) ? (a)-(b)+1 : (*f)-(b)+1+(a) ) )


/****************************************************************************
**
*F  POW_FF(<a>,<n>,<f>) . . . . . . . . . . . . power of a finite field value
**
**  'POW_FF' returns the <n>th power of the finite field value <a>  from  the
**  the finite field pointed to by the pointer <f>.
**
**  Note that 'POW_FF' may only be used if the right operand is an integer in
**  the range $0..order(f)-1$.
**
**  Finally 'POW_FF' may only be used if the product of two integers  of  the
**  size of 'TypFFE' does not cause an overflow, i.e.  only  if  'TypFFE'  is
**  'unsigned short'.
**
**  Note that 'POW_FF' is a macro, so do not  call  it  with  arguments  that
**  have sideeffects.  For optimal performance put the operands in  registers
**  before calling 'POW_FF'.
**
**  If the finite field element is 0 the power is also 0, otherwise  we  have
**  $a^n ~ (z^{a-1})^n = z^{(a-1)*n} = z^{(a-1)*n % (o-1)} ~ (a-1)*n % (o-1)$
*/
#define POW_FF(a,n,f) ( (n)==0 ? 1 : ( (a)==0 ? 0 : \
                        ((((unsigned long)(a))-1) * ((unsigned long)(n))) \
                                                  % (unsigned long)*(f) + 1))


/****************************************************************************
**
*F  RootFiniteField( <q> )  . .  construct the finite field with <q> elements
**
**  'RootFiniteField' returns the handle of the  primitive root of the finite
**  field with <q> elements.  If <q> is not  a power of  a prime or is larger
**  than $2^{16}$  'RootFiniteField' returns the handle 0.    If the field is
**  already constructed,   i.e., is  in  the list   'HdFields' it  is  simply
**  returned.  Otherwise 'RootFiniteField' constructs  this finite  field and
**  remembers it in 'HdFields'.
*/
extern  TypHandle       RootFiniteField (
            unsigned long       q );


/****************************************************************************
**
*F  EvFFE( <hdFFE> )  . . . . . . . . . . . . evaluate a finite field element
**
**  'EvFFE' returns  the value of  the finite field  element  <hdFFE>.  Since
**  finite field elements  are constants  and  thus selfevaluating this  just
**  returns <hdFFE>.
*/
extern  TypHandle       EvFFE (
            TypHandle           hdFFE );


/****************************************************************************
**
*F  SumFFE( <hdL>, <hdR> )  . . . . . . . . . .  sum of finite field elements
**
**  'SumFFE' returns the sum of the  two  finite  field  elements  <hdL>  and
**  <hdR>.  The sum is represented over the field  over  which  the  operands
**  are represented, even if it lies in a much smaller field.
**
**  If one of the operands  is an integer  it is converted  into the field of
**  the  other  operand  before  the  addition.  If  one of the  elements  is
**  represented over a subfield of the field over which  the other element is
**  represented it is lifted into the larger field before the addition.
**
**  Is called from the 'Sum' binop, so both operands are already  evaluated.
*/
extern  TypHandle       SumFFE (
            TypHandle           hdL,
            TypHandle           hdR );


/****************************************************************************
**
*F  DiffFFE( <hdL>, <hdR> ) . . . . . . . difference of finite field elements
**
**  'DiffFFE' returns the difference of the two finite field  elements  <hdL>
**  and <hdR>.  The difference is represented over the field over  which  the
**  operands are represented, even if it lies in a much smaller field.
**
**  If one of the operands  is an integer  it is converted  into the field of
**  the other  operand before the subtraction.  If one  of  the  elements  is
**  represented over a subfield of the field over which  the other element is
**  represented it is lifted into the larger field before the subtraction.
**
**  Is called from the 'Diff' binop, so both operands are already  evaluated.
*/
extern  TypHandle       DiffFFE (
            TypHandle           hdL,
            TypHandle           hdR );


/****************************************************************************
**
*F  ProdFFE( <hdL>, <hdR> ) . . . . . . . .  product of finite field elements
**
**  'ProdFFE' returns the product of the two finite field elements <hdL>  and
**  <hdR>.  The product is represented over the field over which the operands
**  are represented, even if it lies in a much smaller field.
**
**  If one of the operands  is an integer  it is converted  into the field of
**  the other  operand before the multiplication.  If  one of the elements is
**  represented over a subfield of the field over which  the other element is
**  represented it is lifted into the larger field before the multiplication.
**
**  Is called from the 'Prod' binop, so both operands are already  evaluated.
*/
extern  TypHandle       ProdFFE (
            TypHandle           hdL,
            TypHandle           hdR );


/****************************************************************************
**
*F  QuoFFE( <hdL>, <hdR> ) . . . . . . . .  quotient of finite field elements
**
**  'QuoFFE' returns the quotient of the two finite field elements <hdL>  and
**  <hdR>. The quotient is represented over the field over which the operands
**  are represented, even if it lies in a much smaller field.
**
**  If one of the operands  is an integer  it is converted  into the field of
**  the other  operand before  the  division.  If  one  of  the  elements  is
**  represented over a subfield of the field over which  the other element is
**  represented it is lifted into the larger field before the division.
**
**  Is called from the 'Quo' binop, so both operands are already  evaluated.
*/
extern  TypHandle       QuoFFE (
            TypHandle           hdL,
            TypHandle           hdR );


/****************************************************************************
**
*F  PowFFE( <hdL>, <hdR> )  . . . . . . . . . power of a finite field element
**
**  'PowFFE' returns the power of the finite  field  element  <hdL>  and  the
**  integer <hdR>.  The power is represented over the field  over  which  the
**  left operand is represented, even if it lies in a much smaller field.
**
**  Is called from the 'Pow' binop, so both operands are already  evaluated.
*/
extern  TypHandle       PowFFE (
            TypHandle           hdL,
            TypHandle           hdR );


/****************************************************************************
**
*F  EqFFE( <hdL>, <hdR> ) . . . . . . test if finite field elements are equal
**
**  'EqFFE' returns 'HdTrue' if the two finite field elements <hdL> and <hdR>
**  are equal and 'HdFalse' othwise.
**
**  Is called from the 'Eq' binop, so both operands are already  evaluated.
*/
extern  TypHandle       EqFFE (
            TypHandle           hdL,
            TypHandle           hdR );


/****************************************************************************
**
*F  LtFFE( <hdL>, <hdR> ) . . . . .  test if finite field elements is smaller
**
**  'LtFFE' returns 'HdTrue' if the finite field element  <hdL>  is  strictly
**  less than the finite field element <hdR> and 'HdFalse' otherwise.
**
**  Is called from the 'Lt' binop, so both operands are already  evaluated.
*/
extern  TypHandle       LtFFE (
            TypHandle           hdL,
            TypHandle           hdR );


/****************************************************************************
**
*F  PrFFE( <hdFFE> )  . . . . . . . . . . . . .  print a finite field element
**
**  'PrFFE' prints the finite field element <hdFFE>.
*/
extern  void            PrFFE (
            TypHandle           hdFFE );


/****************************************************************************
**
*F  PrFF( <hdField>, <value> )  . . . . . . . . .  print a finite field value
**
**  'PrFF' prints the value <value> from the finite field <hdField>.
**
**  This procedure is called by the 'PrVector' printing procedure, which  can
**  not call 'PrFFE' because it would have to create  finite  field  elements
**  to do so and calling 'NewBag' from a printing procedure is forbidden.
*/
extern  void            PrFF (
            TypHandle           hdField,
            unsigned int        value );


/****************************************************************************
**
*F  FunIsFFE( <hdCall> )  . . . . . . . . . .  test for finite field elements
**
**  'FunIsFFE' implements the internal function 'IsFFE( <obj> )'.
**
**  'IsFFE' returns  'true' if its argument  <obj> is a finite  field element
**  and 'false' otherwise.   'IsFFE' will cause  an  error if  called with an
**  unbound variable.
*/
extern  TypHandle       FunIsFFE (
            TypHandle           hdCall );


/****************************************************************************
**
*F  CharFFE(<hdFFE>)  . . . . . . .  characteristic of a finite field element
**
**  'CharFFE' returns the characteristic of the field  in  which  the  finite
**  field element <hdFFE> lies.
*/
extern  long            CharFFE (
            TypHandle           hdFFE );


/****************************************************************************
**
*F  DegreeFFE(<hdFFE>)  . . . . . . . . . .  degree of a finite field element
**
**  'DegreeFFE' returns the degree of the smallest finite field in which  the
**  finite field element <hdFFE> lies.
*/
extern  long            DegreeFFE (
            TypHandle           hdFFE );


/****************************************************************************
**
*F  FunLogFFE( <hdCall> ) . . . . . . .  logarithm of a finite field constant
**
**  'FunLogFFE' implements the internal function 'LogFFE( <x> )'.
**
**  If called with one argument 'LogFFE' returns the logarithm of the  finite
**  field element <x> with respect to the generator of the finite field  over
**  which <x> is represented.  If called with two arguments 'LogFFE'  returns
**  the logarithm of the finite field element <x> with respect to the  second
**  argument <r> which must lie in the same field like <x>.
**
**  If <x> is 0 'LogFFE' causes an error.
*/
extern  TypHandle       FunLogFFE (
            TypHandle           hdCall );


/****************************************************************************
**
*F  FunIntFFE( <hdCall> ) . . .  convert a finite field element to an integer
**
**  'FunIntFFE' implements the internal function 'IntFFE( <z> )'.
**
**  'IntFFE'  returns  the integer  that  corresponds  to  the  finite  field
**  element <z>, which must of course be  an element  of a prime field, i.e.,
**  the smallest integer <i> such that '<i> * <z>^0 = <z>'.
*/
extern TypHandle ConvTabIntFFE ( long );
extern TypHandle FunIntFFE ( TypHandle hdCall );


/****************************************************************************
**
*F  FunZ( <hdCall> )  . . . . . . . .  return the generator of a finite field
**
**  'FunZ' implements the internal function 'Z( <q> )'.
**
**  'Z' returns the generators of the finite field  with  <q>  elements.  <q>
**  must be a positive prime power.
**
**  'Z' remembers all finite fields that exist in  the list  with the  handle
**  'HdFields' and will not create an already existing field.
*/
extern  TypHandle       FunZ (
            TypHandle           hdCall );


/****************************************************************************
**
*F  InitFF()  . . . . . . . . . . . . . . . . initialize finite field package
**
**  'InitFF' initializes the finite field package.
*/
extern  void            InitFF ( void );
