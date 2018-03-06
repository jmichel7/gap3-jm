/****************************************************************************
**
*A  eval.c                      GAP source                   Martin Schoenert
**
*H  @(#)$Id: eval.c,v 1.1.1.1 1996/12/11 12:43:52 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file contains the main evaluation functions.
**
*H  $Log: eval.c,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:52  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.30  1994/06/10  01:05:05  mschoene
*H  added coding theory package
*H
*H  Revision 3.29  1993/12/17  08:00:36  mschoene
*H  changed 'PrVar' to insert appropriate escape sequences
*H
*H  Revision 3.28  1993/10/26  12:22:51  martin
*H  fixed 'PrBinop' to print '(-1)^n' correctly
*H
*H  Revision 3.27  1993/02/12  17:50:28  martin
*H  added large permutations
*H
*H  Revision 3.26  1993/02/04  10:51:10  martin
*H  changed to the new list interface
*H
*H  Revision 3.25  1992/11/16  19:03:28  martin
*H  added packages 'costab' and 'tietze'
*H
*H  Revision 3.24  1992/08/14  15:32:29  fceller
*H  added 'EnterKernel' and 'ExitKernel' in 'Copy'
*H
*H  Revision 3.23  1992/05/25  08:45:24  fceller
*H  added polynomial package
*H
*H  Revision 3.22  1992/04/27  13:24:48  martin
*H  fixed 'IsBound' for blists
*H
*H  Revision 3.21  1992/01/23  16:21:58  martin
*H  changed the printing, so '~' is only used to break recursion
*H
*H  Revision 3.20  1992/01/09  16:14:10  martin
*H  fixed a minor problem in 'Print' (this is getting out of hand)
*H
*H  Revision 3.19  1992/01/08  11:13:19  martin
*H  changed 'IsBound' and 'Unbind' to accept '<rec>.(<int-expr>)'
*H
*H  Revision 3.18  1992/01/07  11:34:23  martin
*H  improved 'Print' to print '~'
*H
*H  Revision 3.17  1992/01/02  11:09:39  martin
*H  added '<rec>.(<name>)' construct
*H
*H  Revision 3.16  1991/10/01  11:36:05  martin
*H  changed 'Copy' to work with objects that are not trees
*H
*H  Revision 3.15  1991/07/31  13:09:02  fceller
*H  "pcpresen.h" is now included.
*H
*H  Revision 3.14  1991/04/30  16:12:15  martin
*H  initial revision under RCS
*H
*H  Revision 3.13  1991/03/05  12:00:00  martin
*H  added 'Unbind'
*H
*H  Revision 3.12  1991/01/30  12:00:00  martin
*H  improved the permutation package considerably
*H
*H  Revision 3.11  1991/01/21  12:00:00  fceller
*H  added 'LeftQuotient'
*H
*H  Revision 3.10  1991/01/18  12:00:00  martin
*H  fixed 'ShallowCopy' for vectors again
*H
*H  Revision 3.9  1991/01/16  12:00:00  martin
*H  improved the undefined binary operation errors
*H
*H  Revision 3.8  1990/12/12  12:00:00  martin
*H  fixed 'Copy' so that it copies the last data too
*H
*H  Revision 3.7  1990/12/07  13:00:00  martin
*H  changed shifts to please TurboC
*H
*H  Revision 3.6  1990/12/07  12:00:00  martin
*H  added prototypes for function tables
*H
*H  Revision 3.5  1990/12/06  12:00:00  martin
*H  added yet another list package
*H
*H  Revision 3.4  1990/11/20  12:00:00  martin
*H  added new list package
*H
*H  Revision 3.3  1990/11/01  12:00:00  martin
*H  improved comparisons with special case
*H
*H  Revision 3.2  1990/10/02  12:00:00  martin
*H  added 'quit'
*H
*H  Revision 3.1  1990/09/03  12:00:00  martin
*H  fixed 'Copy' and 'ShallowCopy' for vectors
*H
*H  Revision 3.0  1990/08/28  12:00:00  martin
*H  fixed 'IsBound' for vectors and ranges
*H
*/

#include        "system.h"              /* system dependent functions      */
#include        "gasman.h"              /* dynamic storage manager         */
#include        "scanner.h"             /* reading of tokens and printing  */
#include        "idents.h"              /* symbol table managment          */
#include        "integer.h"             /* 'InitInt', large integers       */

#include        "rational.h"            /* 'InitRat'                       */
#include        "cyclotom.h"            /* 'InitCyc'                       */
#include        "unknown.h"             /* 'InitUnknown'                   */
#include        "finfield.h"            /* 'InitFF'                        */
#include        "polynom.h"             /* 'InitPolynom'                   */

#include        "permutat.h"            /* 'InitPermutat'                  */
#include        "word.h"                /* 'InitWord'                      */
#include        "costab.h"              /* 'InitCostab'                    */
#include        "tietze.h"              /* 'InitTietze'                    */
#include        "aggroup.h"             /* 'InitAg'                        */
#include        "pcpresen.h"            /* 'InitPcPres'                    */

#include        "list.h"                /* 'InitList', generic list funcs  */
#include        "plist.h"               /* 'InitPlist', 'LEN_PLIST', ..    */
#include        "set.h"                 /* 'InitSet'                       */
#include        "vector.h"              /* 'InitVector'                    */
#include        "vecffe.h"              /* 'InitVecFFE'                    */
#include        "blister.h"             /* 'InitBlist'                     */
#include        "range.h"               /* 'InitRange'                     */
#include        "string.h"              /* 'InitString', 'IsString'        */

#include        "record.h"              /* 'InitRec'                       */
#include        "statemen.h"            /* 'InitStat'                      */
#include        "function.h"            /* 'InitFunc'                      */
#include        "coding.h"              /* 'InitCoding'                    */

#include        "eval.h"                /* definition part of this package */


/****************************************************************************
**
*V  HdVoid  . . . . . . . . . . . . . . . . . . . . .  handle of the void bag
**
**  'HdVoid' is the handle of the void back, which is returned by procedures,
**  i.e., functions that when viewed at the  GAP  level do not return a value
**  at all.  This plays a role similar to  '*the-non-printing-object*'  which
**  exists in some lisp systems.
*/
TypHandle       HdVoid;


/****************************************************************************
**
*V  HdReturn  . . . . . . . . . . . . . . . . . . .  handle of the return bag
**
**  'HdReturn' is the handle of the bag where 'EvReturn' puts the value of a
**  'return' statement.  This bag is then passed through all  the  statement
**  execution functions all the  way  back  to  'EvFunccall'.  For  'return'
**  statements without an expression 'EvReturn' puts 'HdVoid' into this bag.
*/
TypHandle       HdReturn;


/****************************************************************************
**
*V  HdTrue  . . . . . . . . . . . . . . . . . . . . .  handle of the true bag
*V  HdFalse   . . . . . . . . . . . . . . . . . . . . handle of the false bag
**
**  'HdTrue' is the handle of the unique bag that represents the value 'true'
**  and 'HdFalse' is likewise the unique handle of the  bag  that  represents
**  the value 'HdFalse'.
*/
TypHandle       HdTrue,  HdFalse;


/****************************************************************************
**
*F  EVAL( <hd> )  . . . . . . . . . . . . . . . . . . . .  evaluate an object
**
**  'EVAL' evaluates the bag <hd>  by  calling  the  corresponding  function.
**
**  It is defined in the definition file of this package as followings:
**
#define EVAL(hd)        ((long)(hd)&T_INT ? (hd) : (* EvTab[TYPE(hd)])((hd)))
*/


/****************************************************************************
**
*V  EvTab[<type>] . . . . . . . . evaluation function for bags of type <type>
**
**  Is the main dispatching table that contains for every type a  pointer  to
**  the function that should be executed if a bag  of  that  type  is  found.
*/
TypHandle       (* EvTab[ T_ILLEGAL ]) ( TypHandle hd );


/****************************************************************************
**
*F  CantEval( <hd> )  . . . . . . . . . . . . illegal bag evaluation function
**
**  Is called if a illegal bag should be evaluated, it  generates  an  error.
**  If this is actually ever executed in GAP it  indicates  serious  trouble,
**  for  example  that  the  type  field  of  a  bag  has  been  overwritten.
*/
TypHandle       CantEval ( hd )
    TypHandle           hd;
{
    return Error("Panic: can't eval bag of type %d",(long)TYPE(hd),0L);
}


/****************************************************************************
**
*F  Sum( <hdSum> )  . . . . . . . . . . . . . . . . . . . . .  evaluate a sum
*F  SUM(<hdL>,<hdR>)  . . . . . . . . . . . . . . . . . . . .  evaluate a sum
*V  TabSum[<typeL>][<typeR>]  . . . . . . . . . . table of addition functions
*F  CantSum(<hdL>,<hdR>)  . . . . . . . . . . . . . . . . . . . undefined sum
**
**  'Sum' returns the sum of the two objects '<hdSum>[0]'  and  '<hdSum>[1]'.
**  'Sum' is called from 'EVAL' to eval bags of type 'T_SUM'.
**
**  'Sum' evaluates the operands and then calls the 'SUM' macro.
**
**  'SUM' finds the types of the two operands and uses  them  to  index  into
**  the table 'TabSum' of addition functions.
**
**  At places where performance really matters one should  copy  the  special
**  code from 'Sum' which checks for the addition of two  immediate  integers
**  and computes their sum without calling 'SUM'.
**
**  'SUM' is defined in the header file of this package as follows:
**
#define SUM(hdL,hdR)    ((*TabSum[TYPE(hdL)][TYPE(hdR)])((hdL),(hdR)))
*/
TypHandle       (*TabSum[T_VAR][T_VAR]) ( TypHandle, TypHandle );

TypHandle Sum(TypHandle hd)
{
    TypHandle           hdL,  hdR,  hdS;

    hdL = EVAL( PTR(hd)[0] );  hdR = EVAL( PTR(hd)[1] );

    /* add two small integers with a small sum                             */
    /* check that no overflow occured      */
    if ( ISINT(hdL) && ISINT(hdR)) {
        hdS = CheckedSum(hdL,hdR);
        if(hdS) return hdS;
    }

    return SUM( hdL, hdR );
}

TypHandle       CantSum ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    return Error("operations: sum of %s and %s is not defined",
                 (long)NameType[TYPE(hdL)], (long)NameType[TYPE(hdR)] );
}


/****************************************************************************
**
*F  Diff( <hdDiff> )  . . . . . . . . . . . . . . . . . evaluate a difference
*F  DIFF(<hdL>,<hdR>)   . . . . . . . . . . . . . . . . evaluate a difference
*V  TabDiff[<typeL>][<typeR>] . . . . . . . .  table of subtraction functions
*F  CantDiff(<hdL>,<hdR>) . . . . . . . . . . . . . . .  undefined difference
**
**  'Diff' returns  the  difference of  the  two  objects  '<hdDiff>[0]'  and
**  '<hdDiff>[1]'.  'Diff'  is  called from  'EVAL'  to  eval  bags  of  type
**  'T_DIFF'.
**
**  'Diff' evaluates the operands and then calls the 'DIFF' macro.
**
**  'DIFF' finds the types of the two operands and uses them  to  index  into
**  the table 'TabDiff' of subtraction functions.
**
**  At places where performance really matters one should  copy  the  special
**  code from 'Diff'  which  checks for  the  subtraction  of  two  immediate
**  integers and computes their difference without calling 'DIFF'.
**
**  'DIFF' is defined in the header file of this package as follows:
**
#define DIFF(hdL,hdR)   ((*TabDiff[TYPE(hdL)][TYPE(hdR)])((hdL),(hdR)))
*/
TypHandle       (*TabDiff[T_VAR][T_VAR]) ( TypHandle, TypHandle );

TypHandle Diff(TypHandle hd)
{
    TypHandle  hdL, hdR, hdS;

    hdL = EVAL( PTR(hd)[0] );  hdR = EVAL( PTR(hd)[1] );

    /* subtract two small integers with a small difference                 */
    /* sub and compare top two bits to check that no overflow occured      */
    if ( ISINT(hdL) && ISINT(hdR) ) {
        hdS = CheckedDiff(hdL,hdR);
        if(hdS)return hdS;
    }

    return DIFF(hdL,hdR);
}

TypHandle       CantDiff ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    return Error("operations difference of %s and %s is not defined",
                 (long)NameType[TYPE(hdL)], (long)NameType[TYPE(hdR)] );
}

/****************************************************************************
**
*F  Prod( <hdProd> )  . . . . . . . . . . . . . . . . . .  evaluate a product
*F  PROD(<hdL>,<hdR>)   . . . . . . . . . . . . . . . . .  evaluate a product
*V  TabProd[<typeL>][<typeR>] . . . . . . . table of multiplication functions
*F  CantProd(<hdL>,<hdR>) . . . . . . . . . . . . . . . . . undefined product
**
**  'Prod'  returns   the  product  of  the  two  objects  '<hdProd>[0]'  and
**  '<hdProd>[1]'.  'Prod'  is  called from  'EVAL'  to  eval  bags  of  type
**  'T_PROD'.
**
**  'Prod' evaluates the operands and then calls the 'PROD' macro.
**
**  'PROD' finds the types of the two operands and uses them  to  index  into
**  the table 'TabProd' of multiplication functions.
**
**  At places where performance really matters one should  copy  the  special
**  code from 'Prod'  which  checks for  the  subtraction  of  two  immediate
**  integers and computes their product without calling 'PROD'.
**
**  'PROD' is defined in the header file of this package as follows:
**
#define PROD(hdL,hdR)   ((*TabProd[TYPE(hdL)][TYPE(hdR)])((hdL),(hdR)))
*/
TypHandle       (*TabProd[T_VAR][T_VAR]) ( TypHandle, TypHandle );

TypHandle Prod(TypHandle hd)
{
    TypHandle           hdL,  hdR, hdRR;

    hdL = EVAL( PTR(hd)[0] );  hdR = EVAL( PTR(hd)[1] );

    /* multiply two small integers with a small product                    */
    if ( ISINT(hdL) && ISINT(hdR)) {
        hdRR=CheckedProd(hdL,hdR);
        if(hdRR)return hdRR;
    }

    return PROD( hdL, hdR );
}

TypHandle       CantProd ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    return Error("operations: product of %s and %s is not defined",
                 (long)NameType[TYPE(hdL)], (long)NameType[TYPE(hdR)] );
}

/****************************************************************************
**
*F  Quo( <hdQuo> )  . . . . . . . . . . . . . . . . . . . evaluate a quotient
*F  QUO(<hdL>,<hdR>)  . . . . . . . . . . . . . . . . . . evaluate a quotient
*V  TabQuo[<typeL>][<typeR>]  . . . . . . . . . . table of division functions
*F  CantQuo(<hdL>,<hdR>)  . . . . . . . . . . . . . . . .  undefined quotient
**
**  'Quo'  returns   the   quotient  of  the  two  objects  '<hdQuo>[0]'  and
**  '<hdQuo>[1]'.  'Quo' is called from 'EVAL' to eval bags of type 'T_QUO'.
**
**  'Quo' evaluates the operands and then calls the 'QUO' macro.
**
**  'QUO' finds the types of the two operands and uses  them  to  index  into
**  the table 'TabQuo' of division functions.
**
**  'QUO' is defined in the header file of this package as follows:
**
#define QUO(hdL,hdR)    ((*TabQuo[TYPE(hdL)][TYPE(hdR)])((hdL),(hdR)))
*/
TypHandle       (*TabQuo[T_VAR][T_VAR]) ( TypHandle, TypHandle );

TypHandle Quo(TypHandle hd)
{
    TypHandle           hdL,  hdR;

    hdL = EVAL( PTR(hd)[0] );  hdR = EVAL( PTR(hd)[1] );

    return QUO( hdL, hdR );
}

TypHandle       CantQuo ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    return Error("operations: quotient of %s and %s is not defined",
                 (long)NameType[TYPE(hdL)], (long)NameType[TYPE(hdR)] );
}


/****************************************************************************
**
*F  Mod( <hdMod> )  . . . . . . . . . . . . . . . . . .  evaluate a remainder
*F  MOD(<hdL>,<hdR>)  . . . . . . . . . . . . . . . . .  evaluate a remainder
*V  TabMod[<typeL>][<typeR>]  . . . . . . . . . . table of division functions
*F  CantMod(<hdL>,<hdR>)  . . . . . . . . . . . . . . . . undefined remainder
**
**  'Mod' returns   the  remainder   of  the  two  objects  '<hdMod>[0]'  and
**  '<hdMod>[1]'.  'Mod' is called from 'EVAL' to eval bags of type 'T_MOD'.
**
**  'Mod' evaluates the operands and then calls the 'MOD' macro.
**
**  'MOD' finds the types of the two operands and uses  them  to  index  into
**  the table 'TabMod' of remainder functions.
**
**  'MOD' is defined in the header file of this package as follows:
**
#define MOD(hdL,hdR)    ((*TabMod[TYPE(hdL)][TYPE(hdR)])((hdL),(hdR)))
*/
TypHandle       (*TabMod[T_VAR][T_VAR]) ( TypHandle, TypHandle );

TypHandle       Mod ( hd )
    TypHandle           hd;
{
    TypHandle           hdL,  hdR;

    hdL = EVAL( PTR(hd)[0] );  hdR = EVAL( PTR(hd)[1] );

    return MOD( hdL, hdR );
}

TypHandle       CantMod ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    return Error("operations: remainder of %s and %s is not defined",
                 (long)NameType[TYPE(hdL)], (long)NameType[TYPE(hdR)] );
}


/****************************************************************************
**
*F  Pow( <hdPow> )  . . . . . . . . . . . . . . . . . . . .  evaluate a power
*F  POW(<hdL>,<hdR>)  . . . . . . . . . . . . . . . . . . .  evaluate a power
*V  TabPow[<typeL>][<typeR>]  . . . . . . . table of exponentiation functions
*F  CantPow(<hdL>,<hdR>)  . . . . . . . . . . . . . . . . . . undefined power
**
**  'Pow' returns the power of the two objects '<hdPow>[0]' and '<hdPow>[1]'.
**  'Pow' is called from 'EVAL' to eval bags of type 'T_POW'.
**
**  'Pow' evaluates the operands and then calls the 'POW' macro.
**
**  'POW' finds the types of the two operands and uses  them  to  index  into
**  the table 'TabPow' of powering functions.
**
**  'POW' is defined in the header file of this package as follows:
**
#define POW(hdL,hdR)    ((*TabPow[TYPE(hdL)][TYPE(hdR)])((hdL),(hdR)))
*/
TypHandle       (*TabPow[T_VAR][T_VAR]) ( TypHandle, TypHandle );

TypHandle       Pow ( hd )
    TypHandle           hd;
{
    TypHandle           hdL,  hdR;

    hdL = EVAL( PTR(hd)[0] );  hdR = EVAL( PTR(hd)[1] );

    return POW( hdL, hdR );
}

TypHandle       CantPow ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    return Error("operations: power of %s and %s is not defined",
                 (long)NameType[TYPE(hdL)], (long)NameType[TYPE(hdR)] );
}


/****************************************************************************
**
*F  FunComm( <hdCall> ) . . . . . . . . . . . . . . . . evaluate a commutator
**
**  'FunComm' implements the internal function 'Comm'.
**
**  'Comm( <expr1>, <expr2> )'
**
**  'Comm' returns the commutator of  the  two  group  elements  <expr1>  and
**  <expr2>, i.e., '<expr1>^-1 * <expr2>^-1 * <expr1> * <expr2>'.
**
**  This is a hack to replace the commutator operator until I have fixed  the
**  parser to read something like '(a & b)'
*/
TypHandle       (*TabComm[T_VAR][T_VAR]) ( TypHandle, TypHandle );

TypHandle       IntComm ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdL, hdR;

    /* check the arguments                                                 */
    if ( SIZE(hdCall) != 3 * SIZE_HD )
        return Error("usage: Comm( <expr>, <expr> )",0L,0L);

    /* evaluate the arguments and jump through the function table          */
    hdL = EVAL( PTR(hdCall)[1] );  hdR = EVAL( PTR(hdCall)[2] );
    return (* TabComm[ TYPE(hdL) ][ TYPE(hdR) ]) ( hdL, hdR );
}

TypHandle       CantComm ( hdL, hdR )
    TypHandle           hdL;
    TypHandle           hdR;
{
    return Error("operations: commutator of %s and %s is not defined",
                 (long)NameType[TYPE(hdL)], (long)NameType[TYPE(hdR)] );
}


/****************************************************************************
**
*F  FunLeftQuotient( <hdCall> ) . . . . . . . . . .  evaluate a left quotient
**
**  'FunLeftQuotient' implements the internal function 'LeftQuotient'.
**
**  'LeftQuotient( <expr1>, <expr2> )'
**
**  'LeftQuotient'  returns  the  left  quotient  of  the  two group elements
**  <expr1> and <expr2>, i.e., '<expr1>^-1 * <expr2>'.
*/
TypHandle       FunLeftQuotient ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdL, hdR;

    /** check the arguments ************************************************/
    if ( SIZE( hdCall ) != 3 * SIZE_HD )
        return Error( "usage: LeftQuotient( <expr>, <expr> )", 0L, 0L );

    /** evaluate the arguments and jump through the function table *********/
    hdL = EVAL( PTR( hdCall )[ 1 ] );
    hdR = EVAL( PTR( hdCall )[ 2 ] );
    return ( * TabMod[ TYPE(hdL) ][ TYPE(hdR) ] ) ( hdL, hdR );
}


/****************************************************************************
**
*F  Eq( <hdEq> )  . . . . . . . . . . . . . . . . .  test if <objL> =  <objR>
*F  EQ(<hdL>,<hdR>) . . . . . . . . . . . . . . . .  test if <objL> =  <objR>
*V  TabEq[<typeL>][<typeR>] . . . . . . . . . . table of comparison functions
**
**  'Eq' returns 'HdTrue' if the object '<hdEq>[0]' is equal  to  the  object
**  '<hdEq>[1]' and 'HdFalse'  otherwise.  'Eq'  is  called  from  'EVAL'  to
**  evaluate bags of type 'T_EQ'.
**
**  'Eq' evaluates the operands and then calls the 'EQ' macro.
**
**  'EQ' finds the types of the two operands and  uses  them  to  index  into
**  the table 'TabEq' of comparison functions.
**
**  At places where performance really matters one should  copy  the  special
**  code from 'Eq'  which  checks for the comparison  of  immediate  integers
**  and computes the result without calling 'EQ'.
**
**  'EQ' is defined in the header file of this package as follows:
**
#define EQ(hdL,hdR)     ((*TabEq[TYPE(hdL)][TYPE(hdR)])((hdL),(hdR)))
*/
TypHandle       (*TabEq[T_VAR][T_VAR]) ( TypHandle, TypHandle );

TypHandle       Eq ( hd )
    TypHandle           hd;
{
    TypHandle           hdL,  hdR;

    hdL = EVAL( PTR(hd)[0] );  hdR = EVAL( PTR(hd)[1] );

    /* if the handles are equal the objects certainly will be equal too    */
    if ( hdL == hdR )
        return HdTrue;

    /* Special code to compare two immediate integers.                     */
    if ( ((intptr_t)hdL & (intptr_t)hdR & T_INT) ) {
        if ( HD_TO_INT(hdL) == HD_TO_INT(hdR) )  return HdTrue;
        else                                     return HdFalse;
    }

    return EQ( hdL, hdR );
}


/****************************************************************************
**
*F  Lt( <hdLt> )  . . . . . . . . . . . . . . . . .  test if <objL> <  <objR>
*F  LT(<hdL>,<hdR>) . . . . . . . . . . . . . . . .  test if <objL> <  <objR>
*V  TabLt[<typeL>][<typeR>] , . . . . . . . . . table of comparison functions
**
**  'Lt'  returns 'HdTrue' if  the object '<hdLt>[0]' is less than the object
**  '<hdLt>[1]' and  'HdFalse'  otherwise.  'Lt'  is  called from  'EVAL'  to
**  evaluate bags of type 'T_LT'.
**
**  'Lt' evaluates the operands and then calls the 'LT' macro.
**
**  'LT' finds the types of the two operands and  uses  them  to  index  into
**  the table 'TabLt' of comparison functions.
**
**  At places where performance really matters one should  copy  the  special
**  code  from 'Lt' which  checks for the comparison  of  immediate  integers
**  and computes the result without calling 'LT'.
**
**  'LT' is defined in the header file of this package as follows:
**
#define LT(hdL,hdR)     ((*TabLt[TYPE(hdL)][TYPE(hdR)])((hdL),(hdR)))
*/
TypHandle       (*TabLt[T_VAR][T_VAR]) ( TypHandle, TypHandle );

TypHandle       Lt ( hd )
    TypHandle           hd;
{
    TypHandle           hdL,  hdR;

    hdL = EVAL( PTR(hd)[0] );  hdR = EVAL( PTR(hd)[1] );

    /* if the handles are equal the objects certainly will be equal too    */
    if ( hdL == hdR )
        return HdFalse;

    /* Special code to compare two immediate integers.                     */
    if ( ((intptr_t)hdL & (intptr_t)hdR & T_INT) ) {
        if ( HD_TO_INT(hdL) < HD_TO_INT(hdR) )  return HdTrue;
        else                                    return HdFalse;
    }

    return LT( hdL, hdR );
}


/****************************************************************************
**
*F  Ne( <hdNe> )  . . . . . . . . . . . . . . . . .  test if <objL> <> <objR>
**
**  'Ne'  return 'HdTrue' if  the object <objL> is not equal  to  the  object
**  <objR>.  'Ne' is called from 'EVAL' to evaluate bags of type 'T_NE'.
**
**  'Ne' is simply implemented as 'not <objL> = <objR>'.
*/
TypHandle       Ne ( hd )
    TypHandle           hd;
{
    TypHandle           hdL,  hdR;

    hdL = EVAL( PTR(hd)[0] );  hdR = EVAL( PTR(hd)[1] );

    /* if the handles are equal the objects certainly will be equal too    */
    if ( hdL == hdR )
        return HdFalse;

    /* Special code to compare two immediate integers.                     */
    if ( ((intptr_t)hdL & (intptr_t)hdR & T_INT) ) {
        if ( HD_TO_INT(hdL) != HD_TO_INT(hdR) )  return HdTrue;
        else                                     return HdFalse;
    }

    /* compute 'not <objL> = <objR>' and return it                         */
    if ( EQ(hdL,hdR) == HdTrue )  hdL = HdFalse;
    else                          hdL = HdTrue;
    return hdL;
}


/****************************************************************************
**
*F  Le( <hdLe> )  . . . . . . . . . . . . . . . . .  test if <objL> <= <objR>
**
**  'Le' returns 'HdTrue' if the object <objL>  is  less  than  or  equal  to
**  the object <objR> and 'HdFalse' otherwise.  'Le' is  called  from  'EVAL'
**  to evaluate bags of type 'T_LE'.
**
**  'Le' is simply implemented as 'not <objR> < <objL>'.
*/
TypHandle       Le ( hd )
    TypHandle           hd;
{
    TypHandle           hdL,  hdR;

    hdL = EVAL( PTR(hd)[0] );  hdR = EVAL( PTR(hd)[1] );

    /* if the handles are equal the objects certainly will be equal too    */
    if ( hdL == hdR )
        return HdTrue;

    /* Special code to compare two immediate integers.                     */
    if ( ((intptr_t)hdL & (intptr_t)hdR & T_INT) ) {
        if ( HD_TO_INT(hdL) <= HD_TO_INT(hdR) )  return HdTrue;
        else                                     return HdFalse;
    }

    /* compute 'not <objR> < <objL>' and return it                         */
    if ( LT( hdR, hdL ) == HdTrue )  hdL = HdFalse;
    else                             hdL = HdTrue;
    return hdL;
}


/****************************************************************************
**
*F  Gt( <hdLe> )  . . . . . . . . . . . . . . . . .  test if <objL> >  <objR>
**
**  'Gt' returns 'HdTrue' if the object <objL>  is greater  than  the  object
**  <objR> and 'HdFalse' otherwise.  'Gt' is called from 'EVAL'  to  evaluate
**  bags of type 'T_GT'.
**
**  'Gt' is simply implemented as '<objR> < <objL>'.
*/
TypHandle       Gt ( hd )
    TypHandle           hd;
{
    TypHandle    hdL,  hdR;

    hdL = EVAL( PTR(hd)[0] );  hdR = EVAL( PTR(hd)[1] );

    /* if the handles are equal the objects certainly will be equal too    */
    if ( hdL == hdR )
        return HdFalse;

    /* Special code to compare two immediate integers.                     */
    if ( ((intptr_t)hdL & (intptr_t)hdR & T_INT) ) {
        if ( HD_TO_INT(hdL) >  HD_TO_INT(hdR) )  return HdTrue;
        else                                     return HdFalse;
    }

    return LT( hdR, hdL );
}


/****************************************************************************
**
*F  Ge( <hdLe> )  . . . . . . . . . . . . . . . . .  test if <objL> >= <objR>
**
**  'Ge' returns 'HdTrue' if the object  <objL>  is  greater  or  equal  than
**  the object <objR> and 'HdFalse' otherwise.  'Le' is  called  from  'EVAL'
**  to evaluate bags of type 'T_GE'.
**
**  'Ge' is simply implemented as 'not <objL> < <objR>'.
*/
TypHandle       Ge ( hd )
    TypHandle           hd;
{
    TypHandle           hdL,  hdR;

    hdL = EVAL( PTR(hd)[0] );  hdR = EVAL( PTR(hd)[1] );

    /* if the handles are equal the objects certainly will be equal too    */
    if ( hdL == hdR )
        return HdTrue;

    /* Special code to compare two immediate integers.                     */
    if ( ((intptr_t)hdL & (intptr_t)hdR & T_INT) ) {
        if ( HD_TO_INT(hdL) >= HD_TO_INT(hdR) )  return HdTrue;
        else                                     return HdFalse;
    }

    /* compute 'not <objL> < <objR>' and return it                         */
    if ( LT( hdL, hdR ) == HdTrue )  hdL = HdFalse;
    else                             hdL = HdTrue;
    return hdL;
}


/****************************************************************************
**
*F  IsTrue( <hdL>, <hdR> )  . . . . . . . .  default function for comparisons
**
**  'IsTrue' always returns  'HdTrue'  no  matter  what  the  arguments  are.
**  Is is used for those comparisons where already the types of the  operands
**  determines the outcome.  E.g., it is  used  above  the  diagonal  of  the
**  'TabLt' table.
*/
/*ARGSUSED*/
TypHandle       IsTrue ( hdL, hdR )
    TypHandle           hdL, hdR;
{
    return HdTrue;
}


/****************************************************************************
**
*F  IsFalse( <hdL>, <hdR> ) . . . . . . . .  default function for comparisons
**
**  'IsFalse' always returns 'HdFalse' no  matter  what  the  arguments  are.
**  Is is used for those comparisons where already the types of the  operands
**  determines the outcome.  E.g., it is  used  below  the  diagonal  of  the
**  'TabLt' table.
*/
/*ARGSUSED*/
TypHandle       IsFalse ( hdL, hdR )
    TypHandle           hdL, hdR;
{
    return HdFalse;
}


/****************************************************************************
**
*F  EvVar( <hdVar> )  . . . . . . . . . . . . . . . . . . evaluate a variable
**
**  'EvVar' returns the value  of  the  variable  with  the  handle  <hdVar>.
**
**  The value is the only subobject <hdVar>.  If this has the handle '0' then
**  no value has been assigned to the variable and  an  error  is  generated.
*/
TypHandle       EvVar ( hdVar )
    TypHandle           hdVar;
{
    if ( PTR(hdVar)[0] == 0 )
        return Error("Variable: '%s' must have a value",
                     (long)(PTR(hdVar)+1), 0L );
    return PTR(hdVar)[0];
}


/****************************************************************************
**
*F  EvVarAuto( <hdVar> )  . . . . . . . . . . . . . eval an autoread variable
*/
TypHandle       EvVarAuto ( hdVar )
    TypHandle           hdVar;
{
    /* evaluate the value cell, unless it is already a constant            */
    if ( T_VAR <= TYPE( PTR(hdVar)[0] ) ) {
        EVAL( PTR(hdVar)[0] ); /* return value ignored */
        if ( T_VAR <= TYPE( PTR(hdVar)[0] ) ) {
            return Error("AUTO: '%s' must be defined by the evaluation",
                         (long)(PTR(hdVar)+1), 0L );
        }
    }

    /* convert the autoread variable to a normal one                       */
    Retype( hdVar, T_VAR );

    /* return the value                                                    */
    return PTR(hdVar)[0];
}


/****************************************************************************
**
*F  EvVarAss( <hdAss> ) . . . . . . . . . . . . . . . . execute an assignment
**
**  'EvVarAss' assigns the value of '<hdAss>[1]' to the variable '<hdAss>[0]'
**  and returns the value so that it can be printed in the ReadEvalPrintLoop.
**
**  'EvVarAss' is called from 'EVAL' for bags of type 'T_VARASS'.
*/
TypHandle       EvVarAss ( hdAss )
    TypHandle           hdAss;
{
    TypHandle           hdVal;

    hdVal = EVAL( PTR(hdAss)[1] );
    if ( hdVal == HdVoid )
        return Error("Assignment: function must return a value",0L,0L);
    PTR( PTR(hdAss)[0] )[0] = hdVal;
    return hdVal;
}


/****************************************************************************
**
*F  EvBool( <hdBool> )  . . . . . . . . . . . . . .  evaluate a boolean value
**
**  'EvBool' returns the value of the boolean value <hdBool>.  Since  boolean
**  values are constants and thus selfevaluating it just returns <hdBool>.
*/
TypHandle       EvBool ( hdBool )
    TypHandle           hdBool;
{
    return hdBool;
}


/****************************************************************************
**
*F  EvNot( <hdBool> ) . . . . . . . . . . . . . . . .  negate a boolean value
**
**  'EvNot' returns the boolean negation of the boolean value <hdBool>, i.e.,
**  it returns 'HdTrue' if <hdBool> is 'HdFalse' and vica versa.
*/
TypHandle       EvNot ( hdBool )
    TypHandle           hdBool;
{
    /* evaluate the operand                                                */
    hdBool = EVAL( PTR(hdBool)[0] );

    /* check that it is 'true' or 'false' and return the negation          */
    if ( hdBool == HdTrue )
        return HdFalse;
    else if ( hdBool == HdFalse )
        return HdTrue;
    else
        return Error("not: <expr> must evaluate to 'true' or 'false'",0L,0L);
}


/****************************************************************************
**
*F  EvAnd( <hdAnd> )  . . . . . . . . . . .  evaluate a boolean and operation
**
**  'EvAnd' returns the logical and  of  the  two  operand  '<hdAnd>[0]'  and
**  '<hdAnd>[1]' which must be boolean values.
**
**  If '<hdAnd>[0]' is already  'false'  'EvAnd'  returns  'HdFalse'  without
**  evaluating '<hdAnd>[1]'.  This allows constructs like
**
**      if index <= max  and list[index] = 0  then ... fi;
*/
TypHandle       EvAnd ( hd )
    TypHandle           hd;
{
    TypHandle           hd1;

    /* evaluate and check the left operand                                 */
    hd1 = EVAL( PTR(hd)[0] );
    if ( hd1 == HdFalse )
        return HdFalse;
    else if ( hd1 != HdTrue )
        return Error("and: <expr> must evaluate to 'true' or 'false'",0L,0L);

    /* evaluate and check the right operand                                */
    hd1 = EVAL( PTR(hd)[1] );
    if ( hd1 == HdFalse )
        return HdFalse;
    else if ( hd1 != HdTrue )
        return Error("and: <expr> must evaluate to 'true' or 'false'",0L,0L);

    return HdTrue;
}


/****************************************************************************
**
*F  EvOr( <hdOr> )  . . . . . . . . . . . . . evaluate a boolean or operation
**
**  'EvOr' returns the  logical  or  of  the  two  operands  '<hdOr>[0]'  and
**  '<hdOr>[1]' which must be boolean values.
**
**  If '<hdOr>[0]' is already 'true' 'EvOr' returns 'true' without evaluating
**  '<hdOr>[1]'.  This allows constructs like
**
**      if index > max  or list[index] = 0  then ... fi;
*/
TypHandle       EvOr ( hd )
    TypHandle           hd;
{
    TypHandle           hd1;

    /* evaluate and check the left operand                                 */
    hd1 = EVAL( PTR(hd)[0] );
    if ( hd1 == HdTrue )
        return HdTrue;
    else if ( hd1 != HdFalse )
        return Error("or: <expr> must evaluate to 'true' or 'false'",0L,0L);

    /* evaluate and check the right operand                                */
    hd1 = EVAL( PTR(hd)[1] );
    if ( hd1 == HdTrue )
        return HdTrue;
    else if ( hd1 != HdFalse )
        return Error("or: <expr> must evaluate to 'true' or 'false'",0L,0L);

    return HdFalse;
}


/****************************************************************************
**
*F  EqBool( <hdL>, <hdR> )  . . . . . . . . . . .  test if <boolL> =  <boolR>
**
**  'EqBool' returns 'HdTrue' if the  two  boolean  values  <hdL>  and  <hdR>
**  are equal, and 'HdFalse' otherwise.
*/
TypHandle       EqBool ( hdL, hdR )
    TypHandle           hdL,  hdR;
{
    if ( hdL == hdR )  return HdTrue;
    else               return HdFalse;
}


/****************************************************************************
**
*F  LtBool( <hdL>, <hdR> )  . . . . . . . . . . .  test if <boolL> <  <boolR>
**
**  'LtBool' return 'HdTrue' if  the  boolean value <hdL> is  less  than  the
**  boolean value <hdR> and 'HdFalse' otherwise.
*/
TypHandle       LtBool ( hdL, hdR )
    TypHandle           hdL, hdR;
{
    if ( hdL == HdTrue && hdR == HdFalse )  return HdTrue;
    else                                    return HdFalse;
}


/****************************************************************************
**
*F  PrBool( <hdBool> )  . . . . . . . . . . . . . . . . print a boolean value
**
**  'PrBool' prints the boolean value <hdBool>.
*/
void            PrBool ( hd )
    TypHandle           hd;
{
    if ( hd == HdTrue )  Pr("true",0L,0L);
    else                 Pr("false",0L,0L);
}


/****************************************************************************
**
*F  FunIsBool( <hdCall> ) . . . . . . . . . internal function IsBool( <obj> )
**
**  'IsBool' returns 'true' if the object <obj>  is  a  boolean  and  'false'
**  otherwise.  May cause an error if <obj> is an unbound variable.
*/
TypHandle       FunIsBool ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdObj;

    /* evaluate and check the argument                                     */
    if ( SIZE(hdCall) != 2 * SIZE_HD )
        return Error("usage: IsBool( <obj> )",0L,0L);
    hdObj = EVAL( PTR(hdCall)[1] );
    if ( hdObj == HdVoid )
        return Error("IsBool: function must return a value",0L,0L);

    /* return 'true' if <obj> is a boolean and 'false' otherwise           */
    if ( TYPE(hdObj) == T_BOOL )
        return HdTrue;
    else
        return HdFalse;
}


/****************************************************************************
**
*F  FunShallowCopy( <hdCall> )  . . . . . .  make a shallow copy of an object
**
**  'FunShallowCopy' implements the internal functin 'ShallowCopy( <obj> )'.
**
**  'ShallowCopy' makes a copy of the object  <obj>.  If <obj> is not a  list
**  or a record, 'ShallowCopy' simply returns <obj>, since those objects  can
**  never be modified there is no way to distinguish the original object from
**  any copy, so we might as well not copy it.  If < obj>  is  a  list  or  a
**  record 'ShallowCopy' makes a copy of this object,  but does not copy  the
**  subobjects.
*/
TypHandle       FunShallowCopy ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdOld;          /* handle of the old object        */
    TypHandle           * ptOld;        /* pointer to the old object       */
    TypHandle           hdNew;          /* handle of the new object        */
    TypHandle           * ptNew;        /* pointer to the new object       */
    unsigned long       i;              /* loop variable                   */

    /* check the argument                                                  */
    if ( SIZE(hdCall) != 2*SIZE_HD )
        return Error("usage: ShallowCopy( <obj> )",0L,0L);

    /* evaluate the argument                                               */
    hdOld = EVAL( PTR(hdCall)[1] );

    /* for mutable objects copy the bag                                    */
    if ( T_LIST <= TYPE(hdOld) && TYPE(hdOld) < T_VAR ) {
        hdNew = NewBag( TYPE(hdOld), SIZE(hdOld) );
        ptOld = PTR(hdOld);
        ptNew = PTR(hdNew);
        for ( i = (SIZE(hdOld)+SIZE_HD-1)/SIZE_HD; 0 < i; i-- )
            *ptNew++ = *ptOld++;
    }

    /* otherwise return the original object                                */
    else {
        hdNew = hdOld;
    }

    return hdNew;
}


/****************************************************************************
**
*F  Copy( <hdObj> ) . . . . . . . . . . . . . . . .  make a copy of an object
**
**  'Copy' makes a copy of the  object <hdObj>.  If <obj>  is not a list or a
**  record, 'Copy' simply  returns  <obj>, since those  objects can  never be
**  modified there  is no way  to  distinguish the  original object  from any
**  copy, so we might as  well not copy  it.  If <obj>  is a list or a record
**  'Copy' makes a copy of this object,  and calls itself recursively to copy
**  the subobjects.
*/
TypHandle       CopyShadow ( hdOld )
    TypHandle           hdOld;
{
    TypHandle           hdNew;          /* shadow of <hdOld>               */
    TypHandle           hdTmp;          /* shadow of element of <hdOld>    */
    unsigned long       n;              /* number of handles of <hdOld>    */
    unsigned long       i;              /* loop variable                   */

    /* make a shadow of the old bag                                        */
    hdNew = NewBag( TYPE(hdOld), SIZE(hdOld) );
    hdOld->name[2] = 0;

    /* and make recursively shadows of the subobjects                      */
    n = NrHandles( TYPE(hdOld), SIZE(hdOld) );
    for ( i = n; 0 < i; i-- ) {
        if ( PTR(hdOld)[i-1] != 0
          && T_LIST <= TYPE(PTR(hdOld)[i-1])
          && TYPE(PTR(hdOld)[i-1]) < T_VAR
          && PTR(hdOld)[i-1]->name[2] != 0 ) {
            hdTmp = CopyShadow( PTR(hdOld)[i-1] );
            PTR(hdNew)[i-1] = hdTmp;
        }
    }

    /* return the shadow                                                   */
    return hdNew;
}

void            CopyForward ( hdOld, hdNew )
    TypHandle           hdOld;          /* old bag                         */
    TypHandle           hdNew;          /* shadow of <hdOld>               */
{
    unsigned long       n;              /* number of handles of <hdOld>    */
    unsigned long       i;              /* loop variable                   */

    /* set the forward pointer for <hdOld>                                 */
    PTR(hdOld)[-1] = hdNew;

    /* and do that recursively for all subobjects of <hdOld>               */
    n = NrHandles( TYPE(hdOld), SIZE(hdOld) );
    for ( i = n; 0 < i; i-- ) {
        if ( PTR(hdOld)[i-1] != 0 && PTR(hdNew)[i-1] != 0 )
            CopyForward( PTR(hdOld)[i-1], PTR(hdNew)[i-1] );
    }

}

void            CopyCopy ( hdOld, hdNew )
    TypHandle           hdOld;          /* old bag                         */
    TypHandle           hdNew;          /* shadow of <hdOld>               */
{
    unsigned long       n;              /* number of handles of <hdOld>    */
    unsigned long       i;              /* loop variable                   */

    /* copy the data area                                                  */
    n = NrHandles( TYPE(hdOld), SIZE(hdOld) );
    for ( i = (SIZE(hdOld)+SIZE_HD-1)/SIZE_HD; n < i; i-- ) {
        PTR(hdNew)[i-1] = PTR(hdOld)[i-1];
    }

    /* copy the handles area                                               */
    for ( i = n; 0 < i; i-- ) {
        if ( PTR(hdOld)[i-1] != 0 && PTR(hdNew)[i-1] != 0 )
            CopyCopy( PTR(hdOld)[i-1], PTR(hdNew)[i-1] );
        else if ( PTR(hdOld)[i-1] != 0 && TYPE(PTR(hdOld)[i-1]) != T_INT )
            PTR(hdNew)[i-1] = PTR( PTR(hdOld)[i-1] )[-1];
        else
            PTR(hdNew)[i-1] = PTR(hdOld)[i-1];
    }

}

void            CopyCleanup ( hdOld )
    TypHandle           hdOld;
{
    unsigned long       n;              /* number of handles of <hdOld>    */
    unsigned long       i;              /* loop variable                   */

    /* clean up this bag                                                   */
    PTR(hdOld)[-1] = hdOld;
    hdOld->name[2] = 1;

    /* and recursively clean up the rest                                   */
    n = NrHandles( TYPE(hdOld), SIZE(hdOld) );
    for ( i = n; 0 < i; i-- ) {
        if ( PTR(hdOld)[i-1] != 0
          && T_LIST <= TYPE(PTR(hdOld)[i-1])
          && TYPE(PTR(hdOld)[i-1]) < T_VAR
          && PTR(hdOld)[i-1]->name[2] == 0 )
            CopyCleanup( PTR(hdOld)[i-1] );
    }

}

TypHandle       Copy ( hdOld )
    TypHandle           hdOld;
{
    TypHandle           hdNew;          /* copy of <hdOld>                 */

    /* copy mutable objects                                                */
    EnterKernel();
    if ( T_LIST <= TYPE(hdOld) && TYPE(hdOld) < T_VAR ) {
        hdNew = CopyShadow( hdOld );
        CopyForward( hdOld, hdNew );
        CopyCopy( hdOld, hdNew );
        CopyCleanup( hdOld );
    }

    /* for other objects simply return the object                          */
    else {
        hdNew = hdOld;
    }

    /* return the copy                                                     */
    ExitKernel(hdNew);
    return hdNew;
}


/****************************************************************************
**
*F  FunCopy( <hdCall> ) . . . . . . . . . . . . . .  make a copy of an object
**
**  'FunCopy' implements the internal function 'Copy( <obj> )'.
**
**  'Copy' makes a copy of the  object <hdObj>.  If <obj>  is not a list or a
**  record, 'Copy' simply  returns  <obj>, since those  objects can  never be
**  modified there  is no way  to  distinguish the  original object  from any
**  copy, so we might as  well not copy  it.  If <obj>  is a list or a record
**  'Copy' makes a copy of this object,  and calls itself recursively to copy
**  the subobjects.
*/
TypHandle       FunCopy ( hdCall )
    TypHandle           hdCall;
{
    /* check the argument                                                  */
    if ( SIZE(hdCall) != 2 * SIZE_HD )
        return Error("usage: Copy( <obj> )",0L,0L);

    /* return a copy of the object                                         */
    return Copy( EVAL( PTR(hdCall)[1] ) );
}


/****************************************************************************
**
*F  FunIsBound( <hdCall> )  . . . .  test if a variable has an assigned value
**
**  'FunIsBound' implements the internal function 'IsBound( <expr> )'.
**
*/
TypHandle       FunIsBound ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hd, hdList, hdInd, hdRec, hdNam, Result;
    unsigned long       i;              /* loop variable                   */
    char                value [16];     /* <i> as a string                 */
    char                * p;            /* beginning of <i> in <value>     */

    /* check the argument                                                  */
    if ( SIZE(hdCall) != 2*SIZE_HD )
        return Error("usage: IsBound( <obj> )",0L,0L);
    hd = PTR(hdCall)[1];
    if ( TYPE(hd) != T_VAR     && TYPE(hd) != T_VARAUTO
      && TYPE(hd) != T_LISTELM && TYPE(hd) != T_RECELM )
        return Error("IsBound: <obj> must be a variable",0L,0L);

    /* easiest case first                                                  */
    if ( TYPE(hd) == T_VAR ) {
        if ( PTR(hd)[0] != 0 )
            Result = HdTrue;
        else
            Result = HdFalse;
    }

    /* an variable that autoreads a file is considered bound               */
    else if ( TYPE(hd) == T_VARAUTO ) {
        Result = HdTrue;
    }

    /* is a list element bound                                             */
    else if ( TYPE(hd) == T_LISTELM ) {
        hdList = EVAL( PTR(hd)[0] );
        if ( ! IS_LIST(hdList) )
            return Error("IsBound: <list> must be a list",0L,0L);
        hdInd = EVAL( PTR(hd)[1] );
        if ( TYPE(hdInd) != T_INT || HD_TO_INT(hdInd) <= 0 )
            return Error("IsBound: <index> must be positive int",0L,0L);
        if ( HD_TO_INT(hdInd) <= LEN_LIST(hdList)
          && ELMF_LIST(hdList,HD_TO_INT(hdInd)) != 0 )
            Result = HdTrue;
        else
            Result = HdFalse;
    }

    /* is a record element bound                                           */
    else {
        hdRec = EVAL( PTR(hd)[0] );
        hdNam = PTR(hd)[1];
        if ( TYPE(hdNam) != T_RECNAM ) {
            hdNam = EVAL(hdNam);
            if ( IsString( hdNam ) ) {
                hdNam = FindRecname( (char*)PTR(hdNam) );
            }
            else if ( TYPE(hdNam) == T_INT && 0 <= HD_TO_INT(hdNam) ) {
                i = HD_TO_INT(hdNam);
                p = value + sizeof(value);  *--p = '\0';
                do { *--p = '0' + i % 10; } while ( (i /= 10) != 0 );
                hdNam = FindRecname( p );
            }
            else {
                return Error("<rec>.(<name>) <name> must be a string",0L,0L);
            }
        }
        if ( TYPE(hdRec) != T_REC )
            return Error("IsBound: <record> must be a record",0L,0L);
        for ( i = 0; i < SIZE(hdRec)/(2*SIZE_HD); ++i )
            if ( PTR(hdRec)[2*i] == hdNam )
                break;
        if ( i < SIZE(hdRec)/(2*SIZE_HD) )
            Result = HdTrue;
        else
            Result = HdFalse;
    }

    return Result;
}


/****************************************************************************
**
*F  FunUnbind( <hdCall> ) . . . . . . . . . . . . . . . . unassign a variable
**
**  'FunUnbind' implements the internal function 'Unbind( <expr> )'.
*/
TypHandle       FunUnbind ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hd, hdList, hdInd, hdRec, hdNam;
    unsigned long       i;              /* loop variable                   */
    char                value [16];     /* <i> as a string                 */
    char                * p;            /* beginning of <i> in <value>     */

    /* check the argument                                                  */
    if ( SIZE(hdCall) != 2*SIZE_HD )
        return Error("usage: Unbind( <obj> )",0L,0L);
    hd = PTR(hdCall)[1];
    if ( TYPE(hd) != T_VAR     && TYPE(hd) != T_VARAUTO
      && TYPE(hd) != T_LISTELM && TYPE(hd) != T_RECELM )
        return Error("Unbind: <obj> must be a variable",0L,0L);

    /* easiest case first                                                  */
    if ( TYPE(hd) == T_VAR ) {
        PTR(hd)[0] = 0;
    }

    /* an variable that autoreads a file is considered bound               */
    else if ( TYPE(hd) == T_VARAUTO ) {
        Retype( hd, T_VAR );
        PTR(hd)[0] = 0;
    }

    /* is a list element bound                                             */
    else if ( TYPE(hd) == T_LISTELM ) {
        hdList = EVAL( PTR(hd)[0] );
        if ( ! IS_LIST( hdList ) )
            return Error("Unbind: <list> must be a list",0L,0L);
        PLAIN_LIST( hdList );
        Retype( hdList, T_LIST );
        hdInd = EVAL( PTR(hd)[1] );
        if ( TYPE(hdInd) != T_INT || HD_TO_INT(hdInd) <= 0 )
            return Error("Unbind: <index> must be positive int",0L,0L);
        i = HD_TO_INT(hdInd);
        if ( i < LEN_PLIST(hdList) ) {
            SET_ELM_PLIST( hdList, i, 0 );
        }
        else if ( i == LEN_PLIST( hdList ) ) {
            SET_ELM_PLIST( hdList, i, 0 );
            while ( 0 < i && ELM_PLIST( hdList, i ) == 0 )
                i--;
            SET_LEN_PLIST( hdList, i );
        }
    }

    /* is a record element bound                                           */
    else {
        hdRec = EVAL( PTR(hd)[0] );
        hdNam = PTR(hd)[1];
        if ( TYPE(hdNam) != T_RECNAM ) {
            hdNam = EVAL(hdNam);
            if ( IsString( hdNam ) ) {
                hdNam = FindRecname( (char*)PTR(hdNam) );
            }
            else if ( TYPE(hdNam) == T_INT && 0 <= HD_TO_INT(hdNam) ) {
                i = HD_TO_INT(hdNam);
                p = value + sizeof(value);  *--p = '\0';
                do { *--p = '0' + i % 10; } while ( (i /= 10) != 0 );
                hdNam = FindRecname( p );
            }
            else {
                return Error("<rec>.(<name>) <name> must be a string",0L,0L);
            }
        }
        if ( TYPE(hdRec) != T_REC )
            return Error("Unbind: <record> must be a record",0L,0L);
        for ( i = 0; i < SIZE(hdRec)/(2*SIZE_HD); i++ ) {
            if ( PTR(hdRec)[2*i] == hdNam )
                break;
        }
        if ( i < SIZE(hdRec)/(2*SIZE_HD) ) {
            while ( i < SIZE(hdRec)/(2*SIZE_HD)-1 ) {
                PTR(hdRec)[2*i] = PTR(hdRec)[2*i+2];
                PTR(hdRec)[2*i+1] = PTR(hdRec)[2*i+3];
                i++;
            }
            Resize( hdRec, SIZE(hdRec)-2*SIZE_HD );
        }
    }

    return HdVoid;
}


/****************************************************************************
**
*V  PrTab[<type>] . . . . . . .  printing function for objects of type <type>
**
**  is the main dispatching table that contains for every type a  pointer  to
**  the function that should be executed if a bag  of  that  type  is  found.
*/
void            (* PrTab[ T_ILLEGAL ] ) ( TypHandle hd );


/****************************************************************************
**
*F  Print( <hd> ) . . . . . . . . . . . . . . . . . . . . . . print an object
**
**  'Print'  prints  the  object  with  handle  <hd>.  It dispatches   to the
**  appropriate function stored in 'PrTab[TYPE(<hd>)]'.
*/
TypHandle       HdTildePr;

void            Print ( hd )
    TypHandle           hd;
{
    unsigned long       len;            /* hdObj[1..<len>] are a path from */
    TypHandle           hdObj[256];     /* '~' to <hd>, where hdObj[<i>+1] */
    unsigned long       index[256];     /* is PTR(hdObj[<i>])[index[<i>]]  */
    TypHandle           cur;            /* current object along that path  */
    unsigned long       i;              /* loop variable                   */

    /* check for interrupts                                                */
    if ( SyIsIntr() ) {
        Pr( "%c", (long)'\03', 0L );
        /*N 19-Jun-90 martin do something about the current indent         */
        Error("user interrupt while printing",0L,0L); 
	     /* ignore return value of 'Error'  */
    }

    /* print new objects                                                   */
    if ( TYPE(hd) == T_INT || hd->name[2] != 0 ) {

        /* assign the current object to '~' if this is it                  */
        if ( PTR(HdTildePr)[0] == 0 )
            PTR(HdTildePr)[0] = hd;

        /* mark objects for '~...' detection                               */
        if ( (T_LIST <= TYPE(hd) && TYPE(hd) < T_VAR)
          || TYPE(hd) == T_PERM16 || TYPE(hd) == T_PERM32 )
            hd->name[2] = 0;

        /* dispatch to the appropriate method                              */
        (* PrTab[ TYPE(hd) ] ) (hd);

        /* unmark object again                                             */
        if ( (T_LIST <= TYPE(hd) && TYPE(hd) < T_VAR)
          || TYPE(hd) == T_PERM16 || TYPE(hd) == T_PERM32 )
            hd->name[2] = 1;

        /* unassign '~' again                                              */
        if ( hd == PTR(HdTildePr)[0] )
            PTR(HdTildePr)[0] = 0;

    }

    /* handle common subobject                                             */
    else {

        /* find the subobject in the object again by a backtrack search    */
        len = 0;
        hdObj[0] = HdTildePr;
        index[0] = 0;
        cur = PTR( hdObj[len] )[ index[len] ];
        while ( hd != cur ) {
            for ( i = 0; i <= len && hdObj[i] != cur; i++ ) ;
            if ( cur != 0
              && (TYPE(cur)==T_LIST || TYPE(cur)==T_SET || TYPE(cur)==T_REC)
              && SIZE(cur) != 0
              && len < i ) {
                len++;
                hdObj[len] = cur;
                index[len] = 0;
                cur = PTR( hdObj[len] )[ index[len] ];
            }
            else if ( index[len] < SIZE(hdObj[len])/SIZE_HD-1 ) {
                index[len]++;
                cur = PTR( hdObj[len] )[ index[len] ];
            }
            else {
                if ( len != 0 )  len--;
                cur = 0;
            }
        }

        /* print the path just found                                       */
        for ( i = 0; i <= len; i++ ) {
            if ( TYPE(hdObj[i]) == T_VAR )
                Pr("~",0L,0L);
            else if ( TYPE(hdObj[i])==T_LIST || TYPE(hdObj[i])==T_SET )
                Pr("[%d]",index[i],0L);
            else
                Pr(".%s",(long)PTR(PTR(hdObj[i])[index[i]-1]),0L);
        }

    }

}


/****************************************************************************
**
*F  CantPrint( <hd> ) . . . . . . . . . . . . . illegal bag printing function
**
**  Is called if a illegal bag should be  printed,  it  generates  an  error.
**  If this is actually ever executed in GAP it  indicates  serious  trouble,
**  for  example  that  the  type  field  of  a  bag  has  been  overwritten.
*/
void            CantPrint ( hd )
    TypHandle           hd;
{
    Error("Panic: can't print bag of type %d",(long)TYPE(hd),0L);
}


/****************************************************************************
**
*F  PrVar( <hdVar> )  . . . . . . . . . . . . . . . . . . .  print a variable
**
**  'PrVar' prints  the variable <hdVar>, or precisly  the identifier of that
**  variable.
*/
void            PrVar ( hdVar )
    TypHandle           hdVar;
{
    char *              name;

    /* check for a keyword                                                 */
    name = (char*)(PTR(hdVar)+1);
    if ( !SyStrcmp(name,"and")      || !SyStrcmp(name,"do")
      || !SyStrcmp(name,"elif")     || !SyStrcmp(name,"else")
      || !SyStrcmp(name,"end")      || !SyStrcmp(name,"fi")
      || !SyStrcmp(name,"for")      || !SyStrcmp(name,"function")
      || !SyStrcmp(name,"if")       || !SyStrcmp(name,"in")
      || !SyStrcmp(name,"local")    || !SyStrcmp(name,"mod")
      || !SyStrcmp(name,"not")      || !SyStrcmp(name,"od")
      || !SyStrcmp(name,"or")       || !SyStrcmp(name,"repeat")
      || !SyStrcmp(name,"return")   || !SyStrcmp(name,"then")
      || !SyStrcmp(name,"until")    || !SyStrcmp(name,"while")
      || !SyStrcmp(name,"quit") ) {
        Pr("\\",0L,0L);
    }

    /* print the name                                                      */
    for ( name = (char*)(PTR(hdVar)+1); *name != '\0'; name++ ) {
        if ( IsAlpha(*name) || IsDigit(*name) || *name == '_' )
            Pr("%c",(long)(*name),0L);
        else
            Pr("\\%c",(long)(*name),0L);
    }

}


/****************************************************************************
**
*F  PrVarAss( <hdAss> ) . . . . . . . . . . . . . . . . . print an assignment
**
**  'PrVarAss' prints an assignment to a variable: '<Var> := <Expr>;'
**
**  Linebreaks are preffered before the ':='.
*/
void            PrVarAss ( hdAss )
    TypHandle           hdAss;
{
    Pr("%2>",0L,0L);
    Print(PTR(hdAss)[0]);
    Pr("%< %>:= ",0L,0L);
    Print(PTR(hdAss)[1]);
    Pr("%2<",0L,0L);
}


/****************************************************************************
**
*V  prPrec  . . . . . . . . . . . . . . . . . . . . current preceedence level
**
**  This variable contains the current preceedence  level,  i.e.  an  integer
**  that indicates the binding  power  of  the  currently  printed  operator.
**  If one of the operands is an operation that has lower binding power it is
**  printed in parenthesis.  If the right operand has the same binding  power
**  it is put in parenthesis,  since all the operations are left associative.
**  Preceedence: 12: ^; 10: mod,/,*; 8: -,+; 6: in,=; 4: not; 2: and,or.
**  This sometimes puts in superflous  parenthesis:  2 * f( (3 + 4) ),  since
**  it doesn't know that a  function  call  adds  automatically  parenthesis.
*/
long            prPrec;


/****************************************************************************
**
*F  PrNot( <hdNot> )  . . . . . . . . . . . . .  print a boolean not operator
**
**  'PrNot' print a not operation in the following form: 'not <expr>'.
*/
void            PrNot ( hdNot )
    TypHandle           hdNot;
{
    long                oldPrec;

    oldPrec = prPrec;  prPrec = 4;
    Pr("not%> ",0L,0L);  Print( PTR(hdNot)[0] );  Pr("%<",0L,0L);
    prPrec = oldPrec;
}


/****************************************************************************
**
*F  PrBinop( <hdOp> ) . . . . . . . . . . . . . . .  prints a binary operator
**
**  This prints any of the binary operator using  prPrec  for parenthesising.
*/
void            PrBinop ( hdOp )
    TypHandle           hdOp;
{
    long                oldPrec;
    char                * op;

    oldPrec = prPrec;

    switch ( TYPE(hdOp) ) {
    case T_AND:    op = "and";  prPrec = 2;   break;
    case T_OR:     op = "or";   prPrec = 2;   break;
    case T_EQ:     op = "=";    prPrec = 6;   break;
    case T_LT:     op = "<";    prPrec = 6;   break;
    case T_GT:     op = ">";    prPrec = 6;   break;
    case T_NE:     op = "<>";   prPrec = 6;   break;
    case T_LE:     op = "<=";   prPrec = 6;   break;
    case T_GE:     op = ">=";   prPrec = 6;   break;
    case T_IN:     op = "in";   prPrec = 6;   break;
    case T_SUM:    op = "+";    prPrec = 8;   break;
    case T_DIFF:   op = "-";    prPrec = 8;   break;
    case T_PROD:   op = "*";    prPrec = 10;  break;
    case T_QUO:    op = "/";    prPrec = 10;  break;
    case T_MOD:    op = "mod";  prPrec = 10;  break;
    case T_POW:    op = "^";    prPrec = 12;  break;
    default:       op = "<bogus-operator>";   break;
    }

    if ( oldPrec > prPrec )  Pr("%>(%>",0L,0L);
    else                     Pr("%2>",0L,0L);
    if ( TYPE(hdOp) == T_POW
      && ((TYPE(PTR(hdOp)[0]) == T_INT && HD_TO_INT(PTR(hdOp)[0]) < 0)
        || TYPE(PTR(hdOp)[0]) == T_INTNEG) )
        Pr("(",0L,0L);
    Print( PTR(hdOp)[0] );
    if ( TYPE(hdOp) == T_POW
      && ((TYPE(PTR(hdOp)[0]) == T_INT && HD_TO_INT(PTR(hdOp)[0]) < 0)
        || TYPE(PTR(hdOp)[0]) == T_INTNEG) )
        Pr(")",0L,0L);
    Pr("%2< %2>%s%> %<",(long)op,0L);
    ++prPrec;
    Print( PTR(hdOp)[1] );
    --prPrec;
    if ( oldPrec > prPrec )  Pr("%2<)",0L,0L);
    else                     Pr("%2<",0L,0L);
    prPrec = oldPrec;
}


/****************************************************************************
**
*F  PrComm( <hdComm> )  . . . . . . . . . . . . . . . . .  print a commutator
**
**  This prints a commutator.
*/
void            PrComm ( hd )
    TypHandle           hd;
{
    Pr("%>Comm(%> ",0L,0L);
    Print(PTR(hd)[0]);
    Pr("%<,%>",0L,0L);
    Print(PTR(hd)[1]);
    Pr("%2<)",0L,0L);
}


/****************************************************************************
**
*F  InstEvFunc( <type>, <func> ) . . . . . . .  install a evaluation function
**
**  Installs the function  <func> as evaluation function for bags of  <type>.
*/
void            InstEvFunc ( type, func )
    unsigned int        type;
    TypHandle           (* func) ();
{
    EvTab[ type ] = func;
}


/****************************************************************************
**
*F  InstBinOp( <tab>, <typeL>, <typeR>, <func> )  .  install binary operation
**
**  Installs the function  <func>  as  evaluation  function  for  the  binary
**  operation with the table <tab> for operands of type  <typeL> and <typeR>.
*/
void            InstBinOp ( table, leftType, rightType, func )
    TypHandle           (* table [T_VAR][T_VAR]) ();
    unsigned int        leftType,  rightType;
    TypHandle           (* func) ();
{
    table[ leftType ][ rightType ] = func;
}


/****************************************************************************
**
*F  InstPrFunc( <type>, <func> )  . . . . . . . . install a printing function
**
**  Installs the function <func> as printing function  for  bags  of  <type>.
*/
void            InstPrFunc ( type, func )
    unsigned int        type;
    void                (* func) ();
{
    PrTab[ type ] = func;
}


/****************************************************************************
**
*F  InstVar( <name>, <hdVal> )  . . . . . . . . . . . . . installs a variable
**
**  Installs the value <hdVal> ar value of the new variable with name <name>.
*/
void            InstVar ( name, hdVal )
    char                * name;
    TypHandle           hdVal;
{
    TypHandle           hdVar;

    hdVar = FindIdent( name );
    if ( PTR(hdVar)[0] != 0 )
        Error("Panic: symbol clash %s during initialization",(long)name,0L);
    PTR(hdVar)[0] = hdVal;
}


/****************************************************************************
**
*F  InstIntFunc( <name>, <func> ) . . . . . . .  install an internal function
**
**  Installs the function <func> as internal function with the  name  <name>.
*/
void            InstIntFunc ( name, func )
    char                name [];
    TypHandle           (* func) ();
{
    TypHandle           hdDef,  hdVar;

    /* nice casts, aren't they?                                            */
    hdDef = NewBag( T_FUNCINT, sizeof(TypHandle(**)()) );
    * (TypHandle(**)())PTR(hdDef) = func;

    hdVar = FindIdent( name );
    if ( PTR(hdVar)[0] != 0 )
        Error("Panic: symbol clash %s during initialization",(long)name,0L);
    PTR(hdVar)[0] = hdDef;
}


/****************************************************************************
**
*F  InitEval  . . . . . . . . . . . . . initialize the evaluator main package
**
**  This is called relative lately during the initialization from  InitGap().
*/
void            InitEval ()
{
    unsigned int        type,  typeL,  typeR;

    /* clear the tables for the evaluation dispatching                     */
    for ( type = 0; type < T_ILLEGAL; ++type ) {
        EvTab[type] = CantEval;
        PrTab[type] = CantPrint;
    }
    for ( typeL = 0; typeL < T_VAR; ++typeL ) {
        for ( typeR = 0; typeR < T_VAR; ++typeR ) {
            TabSum[typeL][typeR]  = CantSum;
            TabDiff[typeL][typeR] = CantDiff;
            TabProd[typeL][typeR] = CantProd;
            TabQuo[typeL][typeR]  = CantQuo;
            TabMod[typeL][typeR]  = CantMod;
            TabPow[typeL][typeR]  = CantPow;
            TabComm[typeL][typeR] = CantComm;
        }
    }
    for ( typeL = 0; typeL < T_VAR; ++typeL ) {
        for ( typeR = 0; typeR <= typeL; ++typeR ) {
            TabEq[typeL][typeR] = IsFalse;
            TabLt[typeL][typeR] = IsFalse;
        }
        for ( typeR = typeL+1; typeR < T_VAR; ++typeR ) {
            TabEq[typeL][typeR] = IsFalse;
            TabLt[typeL][typeR] = IsTrue;
        }
    }

    /* install the evaluators main evaluation functions                    */
    InstEvFunc( T_SUM,      Sum      );
    InstEvFunc( T_DIFF,     Diff     );
    InstEvFunc( T_PROD,     Prod     );
    InstEvFunc( T_QUO,      Quo      );
    InstEvFunc( T_MOD,      Mod      );
    InstEvFunc( T_POW,      Pow      );
    /*N hack to replace commutator operator until I fix the parser         */
    /*N InstEvFunc( T_COMM,     Comm     );                                */
    InstIntFunc( "Comm",  IntComm  );
    InstIntFunc( "LeftQuotient",  FunLeftQuotient );
    InstEvFunc( T_EQ,       Eq       );
    InstEvFunc( T_LT,       Lt       );
    InstEvFunc( T_LE,       Le       );
    InstEvFunc( T_NE,       Ne       );
    InstEvFunc( T_GT,       Gt       );
    InstEvFunc( T_GE,       Ge       );

    /* install the main printing functions.                                */
    InstPrFunc( T_SUM,      PrBinop    );
    InstPrFunc( T_DIFF,     PrBinop    );
    InstPrFunc( T_PROD,     PrBinop    );
    InstPrFunc( T_QUO,      PrBinop    );
    InstPrFunc( T_MOD,      PrBinop    );
    InstPrFunc( T_POW,      PrBinop    );
    InstPrFunc( T_COMM,     PrComm     );
    InstPrFunc( T_EQ,       PrBinop    );
    InstPrFunc( T_LT,       PrBinop    );
    InstPrFunc( T_GT,       PrBinop    );
    InstPrFunc( T_NE,       PrBinop    );
    InstPrFunc( T_LE,       PrBinop    );
    InstPrFunc( T_GE,       PrBinop    );
    InstPrFunc( T_IN,       PrBinop    );

    /* variables and assignments                                           */
    InstEvFunc( T_VAR,      EvVar      );
    InstEvFunc( T_VARAUTO,  EvVarAuto  );
    InstEvFunc( T_VARASS,   EvVarAss   );
    InstPrFunc( T_VAR,      PrVar      );
    InstPrFunc( T_VARAUTO,  PrVar      );
    InstPrFunc( T_VARASS,   PrVarAss   );

    /* void bag                                                            */
    HdVoid  = NewBag( T_VOID, 0L );

    /* boolean operations                                                  */
    HdTrue  = NewBag(T_BOOL,0L);  InstVar( "true",  HdTrue  );
    HdFalse = NewBag(T_BOOL,0L);  InstVar( "false", HdFalse );
    InstEvFunc( T_BOOL,     EvBool     );
    InstEvFunc( T_NOT,      EvNot      );
    InstEvFunc( T_AND,      EvAnd      );
    InstEvFunc( T_OR,       EvOr       );
    InstPrFunc( T_BOOL,     PrBool     );
    InstPrFunc( T_NOT,      PrNot      );
    InstPrFunc( T_AND,      PrBinop    );
    InstPrFunc( T_OR,       PrBinop    );
    TabEq[ T_BOOL ][ T_BOOL ] = EqBool;
    TabLt[ T_BOOL ][ T_BOOL ] = LtBool;
    InstIntFunc( "IsBool",      FunIsBool      );

    /* install main evaluator internal functions.                          */
    InstIntFunc( "ShallowCopy", FunShallowCopy );
    InstIntFunc( "Copy",        FunCopy        );
    InstIntFunc( "IsBound",     FunIsBound     );
    InstIntFunc( "Unbind",      FunUnbind      );

    /* install the printing tilde                                          */
    HdTildePr = FindIdent( "~~" );

    /* initialize the evaluators subpackages                               */
    InitInt();                          /* init integer package            */
    InitRat();                          /* init rational package           */
    InitCyc();                          /* init cyclotomic integer package */
    InitUnknown();                      /* init unknown package            */
    InitFF();                           /* init finite field package       */
    InitPolynom();                      /* init polynomial package         */
    InitPermutat();                     /* init permutation package        */
    InitWord();                         /* init word package               */
    InitCosTab();                       /* init coset table package        */
    InitTietze();                       /* init tietze package             */
    InitAg();                           /* init soluable group package     */
    InitPcPres();                       /* init polycyclic pres            */
    InitList();                         /* init list package               */
    InitPlist();                        /* init plain list package         */
    InitSet();                          /* init set package                */
    InitVector();                       /* init vector package             */
    InitVecFFE();                       /* init finite fld vector package  */
    InitBlist();                        /* init boolean list package       */
    InitRange();                        /* init range package              */
    InitString();                       /* init string package             */
    InitRec();                          /* init record package             */
    InitStat();                         /* init statment package           */
    InitFunc();                         /* init function package           */
    InitCoding();                       /* init coding package             */

    /* initialization of further evaluation packages goes here !!!         */
}
