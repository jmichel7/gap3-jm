/****************************************************************************
**
*A  record.c                    GAP source                   Martin Schoenert
**
*H  @(#)$Id: record.c,v 1.1.1.1 1996/12/11 12:43:58 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file contains the functions for computing with records.
**
**  A record with <n> elements is stored as a bag with 2  * <n> entries.  The
**  odd entries are the record names of the elements and the even entries are
**  the corresponding values.
**
*H  $Log: record.c,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:58  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.29  1994/07/06  10:21:44  mschoene
*H  changed 'EvRecAss' to check for void value
*H
*H  Revision 3.28  1993/12/17  09:04:45  mschoene
*H  changed 'PrRec' und 'PrRecElm' to insert appropriate escape backquotes
*H
*H  Revision 3.27  1993/10/26  11:24:27  martin
*H  removed unused variable 'i'
*H
*H  Revision 3.26  1993/10/26  10:49:04  martin
*H  fixed 'MakeRec' from using 'i' for two different purposes
*H
*H  Revision 3.25  1993/02/09  10:40:29  martin
*H  fixed 'LtRec', it must not use 'LT' to compare record names
*H
*H  Revision 3.24  1993/02/04  10:51:10  martin
*H  changed 'T_STRING' to 'T_RECNAME' and to use 'plist' interface
*H
*H  Revision 3.23  1992/11/16  17:24:00  martin
*H  changed 'MakeRec' to check for return value 'HdVoid'
*H
*H  Revision 3.22  1992/04/28  14:08:40  martin
*H  changed a few things to silence GCC
*H
*H  Revision 3.21  1992/02/03  17:05:25  martin
*H  changed record operations function call bags to have a proper name
*H
*H  Revision 3.20  1992/01/28  16:14:54  martin
*H  improved the printing of records
*H
*H  Revision 3.19  1992/01/08  11:56:51  martin
*H  fixed the printing of record elements with evaluation
*H
*H  Revision 3.18  1992/01/08  11:13:19  martin
*H  changed 'EvRecElm', etc. to accept '<rec>.(<int-expr>)'
*H
*H  Revision 3.17  1992/01/06  14:15:15  martin
*H  changed the implementation of '~' slightly
*H
*H  Revision 3.16  1992/01/02  14:44:34  martin
*H  added magic variable '~'
*H
*H  Revision 3.15  1992/01/02  13:12:28  martin
*H  changed 'Backtrace' to handle binary operators
*H
*H  Revision 3.14  1992/01/02  11:52:59  martin
*H  fixed comparison of records
*H
*H  Revision 3.13  1992/01/02  11:09:39  martin
*H  added '<rec>.(<name>)' construct
*H
*H  Revision 3.12  1991/10/14  11:22:32  martin
*H  fixed printing of assignments
*H
*H  Revision 3.11  1991/08/12  09:15:15  fceller
*H  renamed 'GetRecField' to 'RecField'.
*H
*H  Revision 3.10  1991/07/30  14:20:53  martin
*H  fixed 'EqRec' to check for identical objects before calling 'EQ'
*H
*H  Revision 3.9  1991/06/21  17:32:16  fceller
*H  '~.operations.print' renamed to '~.operations.Print',
*H  '~.operations.comm' renmaed to '~.operations.Comm'.
*H
*H  Revision 3.8  1991/06/03  07:01:29  martin
*H  removed unused variables in 'InRec'
*H
*H  Revision 3.7  1991/05/29  15:15:51  martin
*H  added record operation 'in'
*H
*H  Revision 3.6  1991/05/29  14:51:39  martin
*H  changed the record operations function names to '+', '-', etc.
*H
*H  Revision 3.5  1991/05/29  14:46:16  martin
*H  changed the record operations to look at the right operand first
*H
*H  Revision 3.4  1991/04/30  16:12:42  martin
*H  initial revision under RCS
*H
*H  Revision 3.3  1990/12/07  12:00:00  goetz
*H  added 'FunIsRecField' and 'FunRecFields'
*H
*H  Revision 3.2  1990/12/04  12:00:00  martin
*H  improved record ops to make their args collectable
*H
*H  Revision 3.1  1990/09/29  12:00:00  martin
*H  renamed 'GetRec' to 'GetRecField'
*H
*H  Revision 3.0  1990/09/28  12:00:00  martin
*H  added 'GetRec' and 'SetRec'
*H
**
*N  05-Jun-90 martin 'PrRec' should be capable of ignoring elements
*N  05-Jun-90 martin 'PrRec' should support '~.<path>'
*/

#include        "system.h"              /* system dependent functions      */
#include        "gasman.h"              /* dynamic storage management      */
#include        "scanner.h"             /* reading of tokens and printing  */
#include        "eval.h"                /* evaluator main dispatcher       */
#include        "integer.h"             /* arbitrary size integers         */

#include        "idents.h"              /* 'FindRecname'                   */
#include        "plist.h"               /* 'LEN_PLIST', 'SET_LEN_PLIST',.. */
#include        "string.h"              /* 'IsString'                      */

#include        "record.h"              /* declaration part of the package */


/****************************************************************************
**
*F  EvRec( <hdRec> )  . . . . . . . . . . . . . . . . . . . evaluate a record
**
**  'EvRec' evaluates  the record <hdRec>.   Since records are constants  and
**  thus selfevaluating this simply returns <hdRec>.
*/
TypHandle       EvRec ( hdRec )
    TypHandle           hdRec;
{
    return hdRec;
}


/****************************************************************************
**
*V  HdTilde . . . . . . . . . . . . . . . . . . . . . . . . . .  variable '~'
**
**  'HdTilde' is the handle of the variable bag of  the  variable  '~'.  This
**  variable can be used inside record and list literals to refer to the list
**  or record that is currently being created.  So for example '[ ~ ]' is the
**  list that contains itself as only element.
*/
TypHandle       HdTilde;


/****************************************************************************
**
*F  EvMakeRec(<hdLiteral>)  . . .  evaluate record literal to record constant
**
**  'EvMakeRec' evaluates the record literal, i.e., not yet evaluated, record
**  <hdLiteral> to a constant record.
**
**  'EvMakeRec' just calls 'MakeRec' telling it that the result goes into the
**  variable '~'.  Thus expressions in the variable  record can refer to this
**  variable and its subobjects to create objects that are not trees.
*/
TypHandle       EvMakeRec ( hdLiteral )
    TypHandle           hdLiteral;
{
    TypHandle           hdRec;          /* handle of the result            */

    /* top level literal, make the list into '~'                           */
    if ( PTR(HdTilde)[0] == 0 ) {
        hdRec = MakeRec( HdTilde, 0, hdLiteral );
        PTR(HdTilde)[0] = 0;
    }

    /* not top level, do not write the result somewhere                    */
    else {
        hdRec = MakeRec( 0, 0, hdLiteral );
    }

    /* return the result                                                   */
    return hdRec;
}


/****************************************************************************
**
*F  MakeRec(<hdDst>,<ind>,<hdLiteral>)  . evaluate record literal to constant
**
**  'MakeRec' evaluates the record literal <hdLiteral>  to a constant one and
**  puts  the result  into the bag  <hdDst>  at position <ind>.    <hdDst> is
**  either the variable '~', a list, or a record.
**
**  Because of literals like 'rec( a := rec( b := 1, c  := ~.a.b ) )' we must
**  enter the handle  of the result  into the superobject  before we begin to
**  evaluate the record literal.
**
**  A record  literal  is very much  like a  record, except that  at the even
**  places   the do not    yet contain  the  values  of  the  components, but
**  expressions, which yield the elements.   Evaluating a record literal thus
**  means looping over  the components, copying  the names at the odd entries
**  and evaluating the values at the even entries.
*/
TypHandle       MakeRec ( hdDst, ind, hdLiteral )
    TypHandle           hdDst;
    unsigned long       ind;
    TypHandle           hdLiteral;
{
    TypHandle           hdRec;          /* handle of the result            */
    TypHandle           hdNam;          /* handle of component name bag    */
    TypHandle           hdVal;          /* handle of component value       */
    unsigned long       i;              /* loop variable                   */
    unsigned long       k;              /* value from <rec>.(<int>)        */
    char                value [16];     /* <k> as a string                 */
    char                * p;            /* beginning of <k> in <value>     */

    /* make the result bag and enter its handle in the superobject         */
    hdRec = NewBag( T_REC, SIZE(hdLiteral) );
    if ( hdDst != 0 )  PTR(hdDst)[ind] = hdRec;

    /* loop over the components                                            */
    for ( i = 0; i < SIZE(hdLiteral)/SIZE_HD/2; i++ ) {

        /* evaluate the name of the component if it is not constant        */
        hdNam = PTR(hdLiteral)[2*i];
        if ( TYPE(hdNam) != T_RECNAM ) {
            hdNam = EVAL(hdNam);
            if ( IsString( hdNam ) ) {
                hdNam = FindRecname( (char*)PTR(hdNam) );
            }
            else if ( TYPE(hdNam) == T_INT && 0 <= HD_TO_INT(hdNam) ) {
                k = HD_TO_INT(hdNam);
                p = value + sizeof(value);  *--p = '\0';
                do { *--p = '0' + k % 10; } while ( (k /= 10) != 0 );
                hdNam = FindRecname( p );
            }
            else {
                return Error("<rec>.(<name>) <name> must be a string",0L,0L);
            }
        }
        PTR(hdRec)[2*i] = hdNam;

        /* evaluate and enter the value of this component                  */
        if ( TYPE( PTR(hdLiteral)[2*i+1] ) == T_MAKELIST ) {
            MakeList( hdRec, 2*i+1, PTR(hdLiteral)[2*i+1] );
        }
        else if ( TYPE( PTR(hdLiteral)[2*i+1] ) == T_MAKEREC ) {
            MakeRec( hdRec, 2*i+1, PTR(hdLiteral)[2*i+1] );
        }
        else {
            hdVal             = EVAL( PTR(hdLiteral)[2*i+1] );
            while ( hdVal == HdVoid )
                Error("Record: function must return a value",0L,0L);
            PTR(hdRec)[2*i+1] = hdVal;
        }
    }

    /* return the record                                                   */
    return hdRec;
}


/****************************************************************************
**
*F  EvRecElm( <hdElm> ) . . . . . . . . . . . . . . . .  get a record element
**
**  'EvRecElm' returns the value of 'PTR(<hdElm>)[0] . PTR(<hdElm>)[1]'.
**
**  '<record> . <name>'
**
**  This evaluates to the value of the record element with the name <name> in
**  the record <record>.  It is an  error  if  the  record  <record>  has  no
**  element with the name <name>.
**
**  'EvRecElm' simply iterates over  the record looking  for the record name.
**  If it is found the corresponding  value  is returned.  Note  that  it  is
**  not  neccessary  to compare  the strings of  the record names since every
**  record name is  store in a  unique record name  bag, i.e.,  no two record
**  name bags have the same string.  Thus we can simply compare handles.
*/
TypHandle     EvRecElm ( hdElm )
    TypHandle           hdElm;
{
    TypHandle           hdRec,  hdNam,  * ptRec,  *ptEnd;
    unsigned long       k;              /* value from <rec>.(<int>)        */
    char                value [16];     /* <k> as a string                 */
    char                * p;            /* beginning of <k> in <value>     */

    /* first get the record                                                */
    hdRec = EVAL( PTR(hdElm)[0] );
    if ( TYPE(hdRec) != T_REC )
        return Error("Record: left operand must be a record",0L,0L);

    /* then get the right operand, this is by construction a record name   */
    hdNam = PTR(hdElm)[1];
    if ( TYPE(hdNam) != T_RECNAM ) {
        hdNam = EVAL(hdNam);
        if ( IsString( hdNam ) ) {
            hdNam = FindRecname( (char*)PTR(hdNam) );
        }
        else if ( TYPE(hdNam) == T_INT && 0 <= HD_TO_INT(hdNam) ) {
            k = HD_TO_INT(hdNam);
            p = value + sizeof(value);  *--p = '\0';
            do { *--p = '0' + k % 10; } while ( (k /= 10) != 0 );
            hdNam = FindRecname( p );
        }
        else {
            return Error("<rec>.(<name>) <name> must be a string",0L,0L);
        }
    }

    /* find the record name                                                */
    ptRec = PTR(hdRec);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdRec));
    while ( ptRec < ptEnd && ptRec[0] != hdNam )  ptRec += 2;

    /* it is an error if the record name is not found                      */
    if ( ptRec == ptEnd ) {
        return Error("Record: element '%s' must have an assigned value",
                     (long)PTR(hdNam), 0L );
    }

    /* return the value                                                    */
    return ptRec[1];
}


/****************************************************************************
**
*F  EvRecAss( <hdAss> ) . . . . . . . . .  assign a value to a record element
**
**  'EvRecAss' assigns the value  'EVAL( <hdAss>[1] )'  to the  element  with
**  the name '<hdAss>[0][1]' in the record '<hdAss>[0][0]'.
**
**  '<record>.<name> := <expr>;'
**
**  This assigns the value  of the expression  <expr> to the element with the
**  name <name> in  the record <record>.  Further  references to this element
**  will  return this new  value, until another assignment is  performed.  If
**  the record has  no element  with the  name   <name> it is   automatically
**  extended.
*/
TypHandle       EvRecAss ( hdAss )
    TypHandle           hdAss;
{
    TypHandle           hdRec,  hdNam,  hdVal;
    TypHandle           * ptRec,  * ptEnd;
    unsigned long       k;              /* value from <rec>.(<int>)        */
    char                value [16];     /* <k> as a string                 */
    char                * p;            /* beginning of <k> in <value>     */

    /* get the record                                                      */
    hdRec = EVAL( PTR(PTR(hdAss)[0])[0] );
    if ( TYPE(hdRec) != T_REC ) {
        return Error(
          "Record Assignment: left operand must be a record",
                     0L, 0L );
    }

    /* then get the right operand, this is by construction a record name   */
    hdNam = PTR(PTR(hdAss)[0])[1];
    if ( TYPE(hdNam) != T_RECNAM ) {
        hdNam = EVAL(hdNam);
        if ( IsString( hdNam ) ) {
            hdNam = FindRecname( (char*)PTR(hdNam) );
        }
        else if ( TYPE(hdNam) == T_INT && 0 <= HD_TO_INT(hdNam) ) {
            k = HD_TO_INT(hdNam);
            p = value + sizeof(value);  *--p = '\0';
            do { *--p = '0' + k % 10; } while ( (k /= 10) != 0 );
            hdNam = FindRecname( p );
        }
        else {
            return Error("<rec>.(<name>) <name> must be a string",0L,0L);
        }
    }

    /* evaluate the expression                                             */
    hdVal = EVAL( PTR(hdAss)[1] );
    if ( hdVal == HdVoid ) {
        return Error(
          "Record Assignment: function must return a value",
                     0L, 0L );
    }

    /* find the record name                                                */
    ptRec = PTR(hdRec);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdRec));
    while ( ptRec < ptEnd && ptRec[0] != hdNam )  ptRec += 2;

    /* extend the record if the record name is not found                   */
    if ( ptRec == ptEnd ) {
        Resize( hdRec, SIZE(hdRec) + 2*SIZE_HD );
        ptRec = PTR(hdRec) + SIZE(hdRec)/SIZE_HD - 2;
        ptRec[0] = hdNam;
    }

    /* assign the value                                                    */
    ptRec[1] = hdVal;
    return hdVal;
}


/****************************************************************************
**
*V  HdRnOp  . . . . . . . . handle of the 'operations' record name bag, local
*V  HdCall1 . . . . . . . . . handle of a function call bag with 1 arg, local
*V  HdCall2 . . . . . . . . . handle of a function call bag with 2 arg, local
*/
TypHandle       HdRnOp;
TypHandle       HdCall1,  HdCall2;


/****************************************************************************
**
*F  SumRec( <hdL>, <hdR> )  . . . . . . . . . . . . . . . . sum of two record
*V  HdRnSum . . . . . . . . . . .  handle of the 'sum' record name bag, local
*V  HdCallSum . . . . . . . . . . . . . handle of the 'sum' function call bag
**
**  'SumRec' returns the sum of the two operands <hdL> and <hdR>, of which at
**  least one must be a record.
**
**  '<left> + <right>'
**
**  The sum of two records or an object and a record is defined as follows:
**
**  If the right operand <right> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element '+', which is a function,
**  or if the left operand <left> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element '+' which is a function,
**  then this function is called with <left> and <right> as arguments
**    and '<left> + <right>' is the value returned by this function.
**
**  In all other cases an error is raised.
*/
TypHandle       HdRnSum;

TypHandle       HdStrSum;

TypHandle       HdCallSum;

TypHandle       SumRec ( hdL, hdR )
    TypHandle           hdL, hdR;
{
    TypHandle           hdOp;
    TypHandle           * ptRec,  * ptEnd;

    /* if the right operand is a record look for the 'operations' element  */
    if ( TYPE(hdR) != T_REC )  goto l1;
    ptRec = PTR(hdR);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdR));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the '+' element            */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l1;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnSum )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l1;
    PTR(HdCallSum)[0] = ptRec[1];
    PTR(HdCallSum)[1] = hdL;
    PTR(HdCallSum)[2] = hdR;
    hdOp = EVAL( HdCallSum );
    PTR(HdCallSum)[0] = HdStrSum;
    PTR(HdCallSum)[1] = 0;
    PTR(HdCallSum)[2] = 0;
    return hdOp;

l1:
    /* if the left operand is a record look for the 'operations' element   */
    if ( TYPE(hdL) != T_REC )  goto l2;
    ptRec = PTR(hdL);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdL));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the '+' element            */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l2;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnSum )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l2;
    PTR(HdCallSum)[0] = ptRec[1];
    PTR(HdCallSum)[1] = hdL;
    PTR(HdCallSum)[2] = hdR;
    hdOp = EVAL( HdCallSum );
    PTR(HdCallSum)[0] = HdStrSum;
    PTR(HdCallSum)[1] = 0;
    PTR(HdCallSum)[2] = 0;
    return hdOp;

l2:
    return Error("Record: one operand must have '~.operations.+'",0L,0L);
}


/****************************************************************************
**
*F  DiffRec( <hdL>, <hdR> ) . . . . . . . . . . . .  difference of two record
*V  HdRnDiff  . . . . . . . . . . handle of the 'diff' record name bag, local
*V  HdCallDiff  . . . . . . . . . . .  handle of the 'diff' function call bag
**
**  'DiffRec' returns the difference the of two operands  <hdL> and <hdR>, of
**  which at least one must be a record.
**
**  '<left> - <right>'
**
**  The difference of two records or an object and a  record  is  defined  as
**  follows:
**
**  If the right operand <right> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element '-', which is a function,
**  or if the left operand <left> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element '-' which is a function,
**  then this function is called with <left> and <right> as arguments
**    and '<left> - <right>' is the value returned by this function.
**
**  In all other cases an error is raised.
*/
TypHandle       HdRnDiff;

TypHandle       HdStrDiff;

TypHandle       HdCallDiff;

TypHandle       DiffRec ( hdL, hdR )
    TypHandle           hdL, hdR;
{
    TypHandle           hdOp;
    TypHandle           * ptRec,  * ptEnd;

    /* if the right operand is a record look for the 'operations' element  */
    if ( TYPE(hdR) != T_REC )  goto l1;
    ptRec = PTR(hdR);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdR));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the '-' element           */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l1;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnDiff )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l1;
    PTR(HdCallDiff)[0] = ptRec[1];
    PTR(HdCallDiff)[1] = hdL;
    PTR(HdCallDiff)[2] = hdR;
    hdOp = EVAL( HdCallDiff );
    PTR(HdCallDiff)[0] = HdStrDiff;
    PTR(HdCallDiff)[1] = 0;
    PTR(HdCallDiff)[2] = 0;
    return hdOp;

l1:
    /* if the left operand is a record look for the 'operations' element   */
    if ( TYPE(hdL) != T_REC )  goto l2;
    ptRec = PTR(hdL);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdL));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the '-' element            */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l2;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnDiff )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l2;
    PTR(HdCallDiff)[0] = ptRec[1];
    PTR(HdCallDiff)[1] = hdL;
    PTR(HdCallDiff)[2] = hdR;
    hdOp = EVAL( HdCallDiff );
    PTR(HdCallDiff)[0] = HdStrDiff;
    PTR(HdCallDiff)[1] = 0;
    PTR(HdCallDiff)[2] = 0;
    return hdOp;

l2:
    return Error("Record: one operand must have '~.operations.-'",0L,0L);
}


/****************************************************************************
**
*F  ProdRec( <hdL>, <hdR> ) . . . . . . . . . . . . . . product of two record
*V  HdRnProd  . . . . . . . . . . handle of the 'prod' record name bag, local
*V  HdCallProd  . . . . . . . . . . .  handle of the 'prod' function call bag
**
**  'ProdRec' returns the product of the two operands  <hdL>  and  <hdR>,  of
**  which at least one must be a record.
**
**  '<left> * <right>'
**
**  The product of two records or an  object  and  a  record  is  defined  as
**  follows:
**
**  If the right operand <right> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element '*', which is a function,
**  or if the left operand <left> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element '*' which is a function,
**  then this function is called with <left> and <right> as arguments
**    and '<left> * <right>' is the value returned by this function.
**
**  In all other cases an error is raised.
*/
TypHandle       HdRnProd;

TypHandle       HdStrProd;

TypHandle       HdCallProd;

TypHandle       ProdRec ( hdL, hdR )
    TypHandle           hdL, hdR;
{
    TypHandle           hdOp;
    TypHandle           * ptRec,  * ptEnd;

    /* if the right operand is a record look for the 'operations' element  */
    if ( TYPE(hdR) != T_REC )  goto l1;
    ptRec = PTR(hdR);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdR));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the '*' element            */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l1;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnProd )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l1;
    PTR(HdCallProd)[0] = ptRec[1];
    PTR(HdCallProd)[1] = hdL;
    PTR(HdCallProd)[2] = hdR;
    hdOp = EVAL( HdCallProd );
    PTR(HdCallProd)[0] = HdStrProd;
    PTR(HdCallProd)[1] = 0;
    PTR(HdCallProd)[2] = 0;
    return hdOp;

l1:
    /* if the left operand is a record look for the 'operations' element   */
    if ( TYPE(hdL) != T_REC )  goto l2;
    ptRec = PTR(hdL);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdL));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the '*' element            */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l2;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnProd )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l2;
    PTR(HdCallProd)[0] = ptRec[1];
    PTR(HdCallProd)[1] = hdL;
    PTR(HdCallProd)[2] = hdR;
    hdOp = EVAL( HdCallProd );
    PTR(HdCallProd)[0] = HdStrProd;
    PTR(HdCallProd)[1] = 0;
    PTR(HdCallProd)[2] = 0;
    return hdOp;

l2:
    return Error("Record: one operand must have '~.operations.*'",0L,0L);
}


/****************************************************************************
**
*F  QuoRec( <hdL>, <hdR> )  . . . . . . . . . . . . .  quotient of two record
*V  HdRnQuo . . . . . . . . . . .  handle of the 'quo' record name bag, local
*V  HdCallQuo . . . . . . . . . . . . . handle of the 'quo' function call bag
**
**  'QuoRec' returns the quotient of the two operands  <hdL>  and  <hdR>,  of
**  which at least one must be a record.
**
**  '<left> / <right>'
**
**  The quotient of two records or an object  and  a  record  is  defined  as
**  follows:
**
**  If the right operand <right> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element '/', which is a function,
**  or if the left operand <left> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element '/' which is a function,
**  then this function is called with <left> and <right> as arguments
**    and '<left> / <right>' is the value returned by this function.
**
**  In all other cases an error is raised.
*/
TypHandle       HdRnQuo;

TypHandle       HdStrQuo;

TypHandle       HdCallQuo;

TypHandle       QuoRec ( hdL, hdR )
    TypHandle           hdL, hdR;
{
    TypHandle           hdOp;
    TypHandle           * ptRec,  * ptEnd;

    /* if the right operand is a record look for the 'operations' element  */
    if ( TYPE(hdR) != T_REC )  goto l1;
    ptRec = PTR(hdR);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdR));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the '/' element            */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l1;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnQuo )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l1;
    PTR(HdCallQuo)[0] = ptRec[1];
    PTR(HdCallQuo)[1] = hdL;
    PTR(HdCallQuo)[2] = hdR;
    hdOp = EVAL( HdCallQuo );
    PTR(HdCallQuo)[0] = HdStrQuo;
    PTR(HdCallQuo)[1] = 0;
    PTR(HdCallQuo)[2] = 0;
    return hdOp;

l1:
    /* if the left operand is a record look for the 'operations' element   */
    if ( TYPE(hdL) != T_REC )  goto l2;
    ptRec = PTR(hdL);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdL));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the '/' element            */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l2;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnQuo )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l2;
    PTR(HdCallQuo)[0] = ptRec[1];
    PTR(HdCallQuo)[1] = hdL;
    PTR(HdCallQuo)[2] = hdR;
    hdOp = EVAL( HdCallQuo );
    PTR(HdCallQuo)[0] = HdStrQuo;
    PTR(HdCallQuo)[1] = 0;
    PTR(HdCallQuo)[2] = 0;
    return hdOp;

l2:
    return Error("Record: one operand must have '~.operations./'",0L,0L);
}


/****************************************************************************
**
*F  ModRec( <hdL>, <hdR> )  . . . . . . . . . . . . . remainder of two record
*V  HdRnMod . . . . . . . . . . .  handle of the 'mod' record name bag, local
*V  HdCallMod . . . . . . . . . . . . . handle of the 'mod' function call bag
**
**  'ModRec' returns the remainder of the two operands <hdL>  and  <hdR>,  of
**  which at least one must be a record.
**
**  '<left> mod <right>'
**
**  The remainder of two records or an object and  a  record  is  defined  as
**  follows:
**
**  If the right operand <right> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element 'mod', which is a function,
**  or if the left operand <left> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element 'mod' which is a function,
**  then this function is called with <left> and <right> as arguments
**    and '<left> mod <right>' is the value returned by this function.
**
**  In all other cases an error is raised.
*/
TypHandle       HdRnMod;

TypHandle       HdStrMod;

TypHandle       HdCallMod;

TypHandle       ModRec (TypHandle hdL, TypHandle hdR )
{
    TypHandle           hdOp;
    TypHandle           * ptRec,  * ptEnd;

    /* if the right operand is a record look for the 'operations' element  */
    if ( TYPE(hdR) != T_REC )  goto l1;
    ptRec = PTR(hdR);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdR));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the 'mod' element          */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l1;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnMod )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l1;
    PTR(HdCallMod)[0] = ptRec[1];
    PTR(HdCallMod)[1] = hdL;
    PTR(HdCallMod)[2] = hdR;
    hdOp = EVAL( HdCallMod );
    PTR(HdCallMod)[0] = HdStrMod;
    PTR(HdCallMod)[1] = 0;
    PTR(HdCallMod)[2] = 0;
    return hdOp;

l1:
    /* if the left operand is a record look for the 'operations' element   */
    if ( TYPE(hdL) != T_REC )  goto l2;
    ptRec = PTR(hdL);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdL));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the 'mod' element          */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l2;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnMod )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l2;
    PTR(HdCallMod)[0] = ptRec[1];
    PTR(HdCallMod)[1] = hdL;
    PTR(HdCallMod)[2] = hdR;
    hdOp = EVAL( HdCallMod );
    PTR(HdCallMod)[0] = HdStrMod;
    PTR(HdCallMod)[1] = 0;
    PTR(HdCallMod)[2] = 0;
    return hdOp;

l2:
    return Error("Record: one operand must have '~.operations.mod'",0L,0L);
}


/****************************************************************************
**
*F  PowRec( <hdL>, <hdR> )  . . . . . . . . . . . . . . . power of two record
*V  HdRnPow . . . . . . . . . . .  handle of the 'pow' record name bag, local
*V  HdCallPow . . . . . . . . . . . . . handle of the 'pow' function call bag
**
**  'PowRec' returns the power of the two operands <hdL> and <hdR>, of  which
**  at least one must be a record.
**
**  '<left> ^ <right>'
**
**  The power of two records or an object and a record is defined as follows:
**
**  If the right operand <right> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element '^', which is a function,
**  or if the left operand <left> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element '^' which is a function,
**  then this function is called with <left> and <right> as arguments
**    and '<left> ^ <right>' is the value returned by this function.
**
**  In all other cases an error is raised.
*/
TypHandle       HdRnPow;

TypHandle       HdStrPow;

TypHandle       HdCallPow;

TypHandle       PowRec ( hdL, hdR )
    TypHandle           hdL, hdR;
{
    TypHandle           hdOp;
    TypHandle           * ptRec,  * ptEnd;

    /* if the right operand is a record look for the 'operations' element  */
    if ( TYPE(hdR) != T_REC )  goto l1;
    ptRec = PTR(hdR);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdR));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the '^' element            */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l1;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnPow )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l1;
    PTR(HdCallPow)[0] = ptRec[1];
    PTR(HdCallPow)[1] = hdL;
    PTR(HdCallPow)[2] = hdR;
    hdOp = EVAL( HdCallPow );
    PTR(HdCallPow)[0] = HdStrPow;
    PTR(HdCallPow)[1] = 0;
    PTR(HdCallPow)[2] = 0;
    return hdOp;

l1:
    /* if the left operand is a record look for the 'operations' element   */
    if ( TYPE(hdL) != T_REC )  goto l2;
    ptRec = PTR(hdL);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdL));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the '^' element            */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l2;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnPow )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l2;
    PTR(HdCallPow)[0] = ptRec[1];
    PTR(HdCallPow)[1] = hdL;
    PTR(HdCallPow)[2] = hdR;
    hdOp = EVAL( HdCallPow );
    PTR(HdCallPow)[0] = HdStrPow;
    PTR(HdCallPow)[1] = 0;
    PTR(HdCallPow)[2] = 0;
    return hdOp;

l2:
    return Error("Record: one operand must have '~.operations.^'",0L,0L);
}


/****************************************************************************
**
*F  CommRec( <hdL>, <hdR> ) . . . . . . . . . . . .  commutator of two record
*V  HdRnComm  . . . . . . . . . . handle of the 'comm' record name bag, local
*V  HdCallComm  . . . . . . . . . . .  handle of the 'comm' function call bag
**
**  'CommRec' returns the commutator of the two operands <hdL> and  <hdR>, of
**  which at least one must be a record.
**
**  'Comm( <left>, <right> )'
**
**  The  commutator  of  two  records or an object and a record is defined as
**  follows:
**
**  If the right operand <right> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element 'comm', which is a function,
**  or if the left operand <left> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element 'comm' which is a function,
**  then this function is called with <left> and <right> as arguments
**    and '<left> + <right>' is the value returned by this function.
**
**  In all other cases an error is raised.
*/
TypHandle       HdRnComm;

TypHandle       HdStrComm;

TypHandle       HdCallComm;

TypHandle       CommRec ( hdL, hdR )
    TypHandle           hdL, hdR;
{
    TypHandle           hdOp;
    TypHandle           * ptRec,  * ptEnd;

    /* if the right operand is a record look for the 'operations' element  */
    if ( TYPE(hdR) != T_REC )  goto l1;
    ptRec = PTR(hdR);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdR));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the 'comm' element         */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l1;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnComm )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l1;
    PTR(HdCallComm)[0] = ptRec[1];
    PTR(HdCallComm)[1] = hdL;
    PTR(HdCallComm)[2] = hdR;
    hdOp = EVAL( HdCallComm );
    PTR(HdCallComm)[0] = HdStrComm;
    PTR(HdCallComm)[1] = 0;
    PTR(HdCallComm)[2] = 0;
    return hdOp;

l1:
    /* if the left operand is a record look for the 'operations' element   */
    if ( TYPE(hdL) != T_REC )  goto l2;
    ptRec = PTR(hdL);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdL));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the 'comm' element         */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l2;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnComm )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l2;
    PTR(HdCallComm)[0] = ptRec[1];
    PTR(HdCallComm)[1] = hdL;
    PTR(HdCallComm)[2] = hdR;
    hdOp = EVAL( HdCallComm );
    PTR(HdCallComm)[0] = HdStrComm;
    PTR(HdCallComm)[1] = 0;
    PTR(HdCallComm)[2] = 0;
    return hdOp;

l2:
    return Error("Record: one operand must have '~.operations.comm'",0L,0L);
}


/****************************************************************************
**
*F  EqRec( <hdL>, <hdR> ) . . . . . . . . . . .  test if two record are equal
*V  HdRnEq  . . . . . . . . . . . . handle of the 'eq' record name bag, local
*V  HdCallEq  . . . . . . . . . . . . .  handle of the 'eq' function call bag
**
**  'EqRec' returns 'HdTrue' two operands <hdL> and <hdR>, of which at  least
**  one must be a record, are equal and 'HdFalse' otherwise.
**
**  '<left> = <right>'
**
**  The equal of two records or an object and a record is defined as follows:
**
**  If the right operand <right> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element '=', which is a function,
**  or if the left operand <left> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element '=' which is a function,
**  then this function is called with <left> and <right> as arguments
**    and '<left> = <right>' is the value returned by this function.
**
**  In all other cases the operands are considered equal if both are  records
**  and they have the same names and all corresponding elements equal.
*/
TypHandle       HdRnEq;

TypHandle       HdStrEq;

TypHandle       HdCallEq;

TypHandle       EqRec ( hdL, hdR )
    TypHandle           hdL, hdR;
{
    TypHandle           hdOp;
    TypHandle           * ptRec,  * ptEnd;
    unsigned long       i, k;

    /* if the right operand is a record look for the 'operations' element  */
    if ( TYPE(hdR) != T_REC )  goto l1;
    ptRec = PTR(hdR);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdR));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the '=' element            */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l1;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnEq )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l1;
    PTR(HdCallEq)[0] = ptRec[1];
    PTR(HdCallEq)[1] = hdL;
    PTR(HdCallEq)[2] = hdR;
    hdOp = EVAL( HdCallEq );
    PTR(HdCallEq)[0] = HdStrEq;
    PTR(HdCallEq)[1] = 0;
    PTR(HdCallEq)[2] = 0;
    return hdOp;

l1:
    /* if the left operand is a record look for the 'operations' element   */
    if ( TYPE(hdL) != T_REC )  goto l2;
    ptRec = PTR(hdL);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdL));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the '=' element            */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l2;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnEq )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l2;
    PTR(HdCallEq)[0] = ptRec[1];
    PTR(HdCallEq)[1] = hdL;
    PTR(HdCallEq)[2] = hdR;
    hdOp = EVAL( HdCallEq );
    PTR(HdCallEq)[0] = HdStrEq;
    PTR(HdCallEq)[1] = 0;
    PTR(HdCallEq)[2] = 0;
    return hdOp;

l2:
    /* if no function is applied both operands must be records             */
    if ( TYPE(hdL) != T_REC || TYPE(hdR) != T_REC || SIZE(hdL) != SIZE(hdR) )
        return HdFalse;

    /* loop over all names of the left record                              */
    for ( i = 0; i < SIZE(hdL)/(2*SIZE_HD); ++i ) {

        /* look for that name in the right record                          */
        for ( k = 0; k < SIZE(hdR)/(2*SIZE_HD); ++k ) {

            /* if found compare the elements                               */
            if ( PTR(hdL)[2*i] == PTR(hdR)[2*k] ) {
                if ( PTR(hdL)[2*i+1] != PTR(hdR)[2*k+1]
                  && EQ( PTR(hdL)[2*i+1], PTR(hdR)[2*k+1] ) != HdTrue )
                    return HdFalse;
                break;
            }
        }

        /* if not found the record are not equal                           */
        if ( k == SIZE(hdR)/(2*SIZE_HD) )
            return HdFalse;
    }

    /* everything matched, the records are equal                           */
    return HdTrue;
}


/****************************************************************************
**
*F  LtRec( <hdL>, <hdR> ) . . . . . . test if one record is less than another
*V  HdRnLt  . . . . . . . . . . . . handle of the 'lt' record name bag, local
*V  HdCallLt  . . . . . . . . . . . . .  handle of the 'lt' function call bag
**
**  'LtRec' returns 'HdTrue' the operand <hdL> is less than the operand <hdR>
**  and 'HdFalse' otherwise .  At least one of the operands must be a record.
**
**  '<left> < <right>'
**
**  The lt of two records or an object and a record is defined as follows:
**
**  If the right operand <right> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element '<', which is a function,
**  or if the left operand <left> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element '<' which is a function,
**  then this function is called with <left> and <right> as arguments
**    and '<left> < <right>' is the value returned by this function.
**
**  In all other cases the operands are considered ltual if both are  records
**  and they have the same names and all corresponding elements ltual...
*/
TypHandle       HdRnLt;

TypHandle       HdStrLt;

TypHandle       HdCallLt;

TypHandle       LtRec ( hdL, hdR )
    TypHandle           hdL, hdR;
{
    TypHandle           hdOp;           /* handle of operations record     */
    TypHandle           * ptRec;        /* pointer into the record         */
    TypHandle           * ptEnd;        /* pointer to end of the record    */
    TypHandle           hdNam;          /* handle of a component name      */
    TypHandle           hdVal;          /* handle of a component value     */
    unsigned long       h;              /* gap width in shellsort          */
    unsigned long       i, k;           /* loop variables                  */

    /* if the right operand is a record look for the 'operations' element  */
    if ( TYPE(hdR) != T_REC )  goto l1;
    ptRec = PTR(hdR);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdR));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the '<' element            */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l1;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnLt )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l1;
    PTR(HdCallLt)[0] = ptRec[1];
    PTR(HdCallLt)[1] = hdL;
    PTR(HdCallLt)[2] = hdR;
    hdOp = EVAL( HdCallLt );
    PTR(HdCallLt)[0] = HdStrLt;
    PTR(HdCallLt)[1] = 0;
    PTR(HdCallLt)[2] = 0;
    return hdOp;

l1:
    /* if the left operand is a record look for the 'operations' element   */
    if ( TYPE(hdL) != T_REC )  goto l2;
    ptRec = PTR(hdL);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdL));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the '<' element            */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l2;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnLt )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l2;
    PTR(HdCallLt)[0] = ptRec[1];
    PTR(HdCallLt)[1] = hdL;
    PTR(HdCallLt)[2] = hdR;
    hdOp = EVAL( HdCallLt );
    PTR(HdCallLt)[0] = HdStrLt;
    PTR(HdCallLt)[1] = 0;
    PTR(HdCallLt)[2] = 0;
    return hdOp;

l2:
    /* if no function is applied both operands must be records             */
    if ( TYPE(hdL) < TYPE(hdR) )
        return HdTrue;
    else if ( TYPE(hdR) < TYPE(hdL) )
        return HdFalse;

    /* sort the left record with a shellsort                               */
    h = 1;  while ( 9*h + 4 < SIZE(hdL)/(2*SIZE_HD) )  h = 3*h + 1;
    while ( 0 < h ) {
        for ( i = h+1; i <= SIZE(hdL)/(2*SIZE_HD); i++ ) {
            hdNam = PTR(hdL)[2*i-2];
            hdVal = PTR(hdL)[2*i-1];
            k = i;
            while ( h < k
                 && SyStrcmp( (char*)PTR(hdNam),
                              (char*)PTR( PTR(hdL)[2*(k-h)-2] ) ) < 0 ) {
                PTR(hdL)[2*k-2] = PTR(hdL)[2*(k-h)-2];
                PTR(hdL)[2*k-1] = PTR(hdL)[2*(k-h)-1];
                k -= h;
            }
            PTR(hdL)[2*k-2] = hdNam;
            PTR(hdL)[2*k-1] = hdVal;
        }
        h = h / 3;
    }

    /* sort the right record with a shellsort                              */
    h = 1;  while ( 9*h + 4 < SIZE(hdR)/(2*SIZE_HD) )  h = 3*h + 1;
    while ( 0 < h ) {
        for ( i = h+1; i <= SIZE(hdR)/(2*SIZE_HD); i++ ) {
            hdNam = PTR(hdR)[2*i-2];
            hdVal = PTR(hdR)[2*i-1];
            k = i;
            while ( h < k
                 && SyStrcmp( (char*)PTR(hdNam),
                              (char*)PTR( PTR(hdR)[2*(k-h)-2] ) ) < 0 ) {
                PTR(hdR)[2*k-2] = PTR(hdR)[2*(k-h)-2];
                PTR(hdR)[2*k-1] = PTR(hdR)[2*(k-h)-1];
                k -= h;
            }
            PTR(hdR)[2*k-2] = hdNam;
            PTR(hdR)[2*k-1] = hdVal;
        }
        h = h / 3;
    }

    /* now see what differs                                                */
    for ( i = 1; i <= SIZE(hdR)/(2*SIZE_HD); i++ ) {
        if ( i > SIZE(hdL)/(2*SIZE_HD) ) {
            return HdTrue;
        }
        else if ( PTR(hdL)[2*i-2] != PTR(hdR)[2*i-2] ) {
            if ( SyStrcmp( (char*)PTR( PTR(hdR)[2*i-2] ),
                           (char*)PTR( PTR(hdL)[2*i-2] ) ) < 0 ) {
                return HdTrue;
            }
            else {
                return HdFalse;
            }
        }
        else if ( EQ( PTR(hdL)[2*i-1], PTR(hdR)[2*i-1] ) == HdFalse ) {
            return LT( PTR(hdL)[2*i-1], PTR(hdR)[2*i-1] );
        }
    }

    /* the records are equal, or the right is a proper prefix of the left  */
    return HdFalse;
}


/****************************************************************************
**
*F  InRec( <hdL>, <hdR> ) . . . . . . . . . .  test if a record is in another
*V  HdRnIn  . . . . . . . . . . . . handle of the 'in' record name bag, local
*V  HdCallIn  . . . . . . . . . . . . .  handle of the 'in' function call bag
**
**  'InRec' returns 'HdTrue' the operand <hdL> is in the  operand  <hdR>  and
**  'HdFalse' otherwise .  At least the right operand must be a record.
**
**  '<left> in <right>'
**
**  The 'in' of two records or an object and a record is defined as follows:
**
**  If the right operand <right> is a record
**    and is has a element with the name 'operations', which is a record,
**    and this record has a element 'in', which is a function,
**  then this function is called with <left> and <right> as arguments
**    and '<left> in <right>' is the value returned by this function.
**
**  In all other cases an error is raised.
*/
TypHandle       HdRnIn;

TypHandle       HdStrIn;

TypHandle       HdCallIn;

TypHandle       InRec ( hdL, hdR )
    TypHandle           hdL, hdR;
{
    TypHandle           hdOp;
    TypHandle           * ptRec,  * ptEnd;

    /* if the right operand is a record look for the 'operations' element  */
    if ( TYPE(hdR) != T_REC )  goto l1;
    ptRec = PTR(hdR);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdR));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the 'in' element           */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l1;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnIn )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l1;
    PTR(HdCallIn)[0] = ptRec[1];
    PTR(HdCallIn)[1] = hdL;
    PTR(HdCallIn)[2] = hdR;
    hdOp = EVAL( HdCallIn );
    PTR(HdCallIn)[0] = HdStrIn;
    PTR(HdCallIn)[1] = 0;
    PTR(HdCallIn)[2] = 0;
    return hdOp;

l1:
    /* if no function is applied both operands must be records             */
    return Error("Record: right operand must have '~.operations.in'",0L,0L);
}


/****************************************************************************
**
*F  PrRec( <hdRec> )  . . . . . . . . . . . . . . . . . . . .  print a record
*V  HdRnPrint . . . . . . . . .  handle of the 'print' record name bag, local
*V  HdCallPrint . . . . . . . . . . . handle of the 'print' function call bag
**
**  'PrRec' prints the record with the handle <hdRec>.
**
**  If <hdRec> has an element 'operations' which is a record
**    and this record as an element 'print' which is a function
**  then this function is called with <hdRec> as argument and should print it
**
**  In all other cases the record is printed in the following form:
**
**  'rec( <name> := <expr>,... )'
**
**  'PrRec' is also called to print variable records, i.e., records that have
**  not yet been evaluated.  They are always printed in the second form.
*/
TypHandle       HdRnPrint;

TypHandle       HdStrPrint;

TypHandle       HdCallPrint;

void            PrRec ( hdRec )
    TypHandle           hdRec;
{
    TypHandle           hdOp;
    TypHandle           * ptRec,  * ptEnd;
    unsigned long       i;
    char                * name;

    /* if the left operand is a record look for the 'operations' element   */
    if ( TYPE(hdRec) != T_REC )  goto l1;
    ptRec = PTR(hdRec);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdRec));
    while ( ptRec < ptEnd && ptRec[0] != HdRnOp )  ptRec += 2;

    /* if it was found and is a record look for the 'print' element        */
    if ( ptRec == ptEnd || TYPE(ptRec[1]) != T_REC )  goto l1;
    hdOp = ptRec[1];
    ptRec = PTR(hdOp);
    ptEnd = (TypHandle*)((char*)ptRec + SIZE(hdOp));
    while ( ptRec < ptEnd && ptRec[0] != HdRnPrint )  ptRec += 2;

    /* if it was found and is a function then apply it                     */
    if ( ptRec == ptEnd )  goto l1;
    PTR(HdCallPrint)[0] = ptRec[1];
    PTR(HdCallPrint)[1] = hdRec;
    EVAL( HdCallPrint ); /* value ignored */
    PTR(HdCallPrint)[0] = HdStrPrint;
    PTR(HdCallPrint)[1] = 0;
    return;

l1:
    /* otherwise print the record in the usual form                        */
    /*N 05-Jun-90 martin 'PrRec' should be capable of ignoring elements    */
    /*N 05-Jun-90 martin 'PrRec' should support '~.<path>'                 */
    Pr("%2>rec(\n%2>",0L,0L);
    for ( i = 0; i < SIZE(hdRec)/(2*SIZE_HD); ++i ) {

        /* print an ordinary record name                                   */
        if ( TYPE( PTR(hdRec)[2*i] ) == T_RECNAM ) {

            /* check for a keyword                                         */
            name = (char*)PTR(PTR(hdRec)[2*i]);
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

            /* print the name                                              */
            for ( name = (char*)PTR(PTR(hdRec)[2*i]); *name!='\0'; name++ ) {
                if ( IsAlpha(*name) || IsDigit(*name) || *name == '_' )
                    Pr("%c",(long)(*name),0L);
                else
                    Pr("\\%c",(long)(*name),0L);
            }

        }

        /* print an evaluating record name                                 */
        else {
            Pr(" (",0L,0L);
            Print( PTR(hdRec)[2*i] );
            Pr(")",0L,0L);
        }

        /* print the component                                             */
        Pr("%< := %>",0L,0L);
        Print( PTR(hdRec)[2*i+1] );
        if ( i < SIZE(hdRec)/(2*SIZE_HD)-1 )
            Pr("%2<,\n%2>",0L,0L);

    }
    Pr(" %4<)",0L,0L);
}


/****************************************************************************
**
*F  PrRecElm( <hdElm> ) . . . . . . . . . . . . . . .  print a record element
**
**  'PrRecElm' prints the record element in the following form:
**
**  '<record> . <name>'
*/
void            PrRecElm ( hdElm )
    TypHandle           hdElm;
{
    char *              name;

    /* print the record                                                    */
    Pr( "%>", 0L, 0L );
    Print( PTR(hdElm)[0] );

    /* print an ordinary record name                                       */
    if ( TYPE( PTR(hdElm)[1] ) == T_RECNAM ) {

        /* print the dot                                                   */
        Pr("%<.%>",0L,0L);

        /* check for a keyword                                             */
        name = (char*)PTR(PTR(hdElm)[1]);
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

        /* print the name                                                  */
        for ( name = (char*)PTR(PTR(hdElm)[1]); *name != '\0'; name++ ) {
            if ( IsAlpha(*name) || IsDigit(*name) || *name == '_' )
                Pr("%c",(long)(*name),0L);
            else
                Pr("\\%c",(long)(*name),0L);
        }

        /* thats it                                                        */
        Pr("%<",0L,0L);
 
    }

    /* print an evaluating record name                                     */
    else {
        Pr( "%<.%>(", 0L, 0L );
        Print( PTR(hdElm)[1] );
        Pr( ")%<", 0L, 0L );
    }

}


/****************************************************************************
**
*F  PrRecAss( <hdAss> ) . . . . . . . . . . . . . . print a record assignment
**
**  'PrRecAss' prints the record assignment in the form:
**
**  '<record>.<name> := <expr>;'
*/
void            PrRecAss ( hdAss )
    TypHandle           hdAss;
{
    Pr( "%2>", 0L, 0L );
    Print( PTR(hdAss)[0] );
    Pr( "%< %>:= ", 0L, 0L );
    Print( PTR(hdAss)[1] );
    Pr( "%2<", 0L, 0L );
}


/****************************************************************************
**
*F  FunIsRec( <hdCall> )  . . . . . . . . . . . . . internal function 'IsRec'
**
**  'IsRec'  returns 'true' if the object  <obj>  is  a  record  and  'false'
**  otherwise.  May cause an error if <obj> is an unbound variable.
*/
TypHandle       FunIsRec ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdObj;

    /* evaluate and check the argument                                     */
    if ( SIZE(hdCall) != 2 * SIZE_HD )
        return Error("usage: IsRec( <obj> )",0L,0L);
    hdObj = EVAL( PTR(hdCall)[1] );
    if ( hdObj == HdVoid )
        return Error("IsRec: function must return a value",0L,0L);

    /* return 'true' if <obj> is a rational and 'false' otherwise          */
    if ( TYPE(hdObj) == T_REC )
        return HdTrue;
    else
        return HdFalse;
}


/****************************************************************************
**
*F  FunRecFields( <hdCall> )  . . . . . . . . . . . . . .  list record fields
**
**  'FunRecFields' implements the internal function 'RecFields'.
**
**  'RecFields( <rec> )'
**
**  'RecFields' returns a list of strings representing all record  fields  of
**  the record <rec>.
**
**  You must use 'RecFields' if you want to make a selective copy of a record
**  for example to delete record fields arising from 'SetRecField'.
**
**  |    gap> r := rec();;
**      gap> i := 2^3-1;;
**      gap> s := ConcatenationString( "r", String(i) );;
**      gap> SetRecField( r, s, 0 );;
**      gap> r;
**      rec( r7 := 0 )
**      gap> RecFields( r );
**      [ "r7" ] |
*/
TypHandle       FunRecFields ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdRec,  hdNam;
    TypHandle           hdStr;
    long                i;              /* loop variable                   */

    /* get and check the arguments                                         */
    if ( SIZE(hdCall) != 2*SIZE_HD )
        return Error("usage: RecFields( <rec> )",0L,0L);
    hdRec = EVAL( PTR(hdCall)[1] );
    if ( TYPE(hdRec) != T_REC )
        return Error("RecFields: <rec> must be a record",0L,0L);
    hdNam = NewBag( T_LIST, SIZE_PLEN_PLIST( SIZE(hdRec) / SIZE_HD / 2 ) );
    SET_LEN_PLIST( hdNam, SIZE(hdRec) / SIZE_HD / 2 );

    /* list the record names                                               */
    for ( i = 1; i <= LEN_PLIST( hdNam ); i++ ) {
        hdStr = NewBag( T_STRING, SIZE( PTR(hdRec)[2*i-2] ) );
        SyStrncat( (char*)PTR(hdStr),
                   (char*)PTR( PTR(hdRec)[2*i-2] ),
                   SIZE( PTR(hdRec)[2*i-2] ) );
        SET_ELM_PLIST( hdNam, i, hdStr );
    }

    /* return the list                                                     */
    return hdNam;
}


/****************************************************************************
**
*F  InitRec() . . . . . . . . . . . . . . . .  intitialize the record package
**
**  'InitRec' initializes the record package.
*/
void            InitRec ()
{
    long                i;

    /* install the evaluation and printing functions                       */
    InstEvFunc( T_REC,     EvRec     );
    InstEvFunc( T_MAKEREC, EvMakeRec );
    InstEvFunc( T_RECELM,  EvRecElm  );
    InstEvFunc( T_RECASS,  EvRecAss  );
    InstPrFunc( T_REC,     PrRec     );
    InstPrFunc( T_MAKEREC, PrRec     );
    InstPrFunc( T_RECELM,  PrRecElm  );
    InstPrFunc( T_RECASS,  PrRecAss  );

    /* install the magic variable '~'                                      */
    HdTilde = FindIdent( "~" );

    /* find the record name of the operations record                       */
    HdRnOp     = FindRecname( "operations" );

    /* find the record names of operators and create function call bags    */
    HdRnSum     = FindRecname( "+"     );
    HdCallSum   = NewBag( T_FUNCCALL, 3 * SIZE_HD );
    HdStrSum    = FindIdent( "<rec1> + <rec2>" );
    PTR(HdCallSum)[0] = HdStrSum;
    HdRnDiff    = FindRecname( "-"     );
    HdCallDiff  = NewBag( T_FUNCCALL, 3 * SIZE_HD );
    HdStrDiff   = FindIdent( "<rec1> - <rec2>" );
    PTR(HdCallDiff)[0] = HdStrDiff;
    HdRnProd    = FindRecname( "*"     );
    HdCallProd  = NewBag( T_FUNCCALL, 3 * SIZE_HD );
    HdStrProd   = FindIdent( "<rec1> * <rec2>" );
    PTR(HdCallProd)[0] = HdStrProd;
    HdRnQuo     = FindRecname( "/"     );
    HdCallQuo   = NewBag( T_FUNCCALL, 3 * SIZE_HD );
    HdStrQuo    = FindIdent( "<rec1> / <rec2>" );
    PTR(HdCallQuo)[0] = HdStrQuo;
    HdRnMod     = FindRecname( "mod"   );
    HdCallMod   = NewBag( T_FUNCCALL, 3 * SIZE_HD );
    HdStrMod    = FindIdent( "<rec1> mod <rec2>" );
    PTR(HdCallMod)[0] = HdStrMod;
    HdRnPow     = FindRecname( "^"     );
    HdCallPow   = NewBag( T_FUNCCALL, 3 * SIZE_HD );
    HdStrPow    = FindIdent( "<rec1> ^ <rec2>" );
    PTR(HdCallPow)[0] = HdStrPow;
    HdRnComm    = FindRecname( "Comm"  );
    HdCallComm  = NewBag( T_FUNCCALL, 3 * SIZE_HD );
    HdStrComm   = FindIdent( "Comm( <rec1>, <rec2> )" );
    PTR(HdCallComm)[0] = HdStrComm;
    HdRnEq      = FindRecname( "="     );
    HdCallEq    = NewBag( T_FUNCCALL, 3 * SIZE_HD );
    HdStrEq     = FindIdent( "<rec1> = <rec2>" );
    PTR(HdCallEq)[0] = HdStrEq;
    HdRnLt      = FindRecname( "<"     );
    HdCallLt    = NewBag( T_FUNCCALL, 3 * SIZE_HD );
    HdStrLt     = FindIdent( "<rec1> < <rec2>" );
    PTR(HdCallLt)[0] = HdStrLt;

    /* note that 'in' is special because it is not implemented as binop    */
    HdRnIn      = FindRecname( "in"    );
    HdCallIn    = NewBag( T_FUNCCALL, 3 * SIZE_HD );
    HdStrIn     = FindIdent( "<obj> in <rec>" );
    PTR(HdCallIn)[0] = HdStrIn;

    /* note that 'print is special because it is a function not a binop    */
    HdRnPrint   = FindRecname( "Print" );
    HdCallPrint = NewBag( T_FUNCCALL, 2 * SIZE_HD );
    HdStrPrint  = FindIdent( "Print( <rec> )" );
    PTR(HdCallPrint)[0] = HdStrPrint;

    for ( i = T_VOID; i < T_VAR; ++i ) {
        TabSum[  i ][ T_REC ] = SumRec;
        TabSum[  T_REC ][ i ] = SumRec;
        TabDiff[ i ][ T_REC ] = DiffRec;
        TabDiff[ T_REC ][ i ] = DiffRec;
        TabProd[ i ][ T_REC ] = ProdRec;
        TabProd[ T_REC ][ i ] = ProdRec;
        TabQuo[  i ][ T_REC ] = QuoRec;
        TabQuo[  T_REC ][ i ] = QuoRec;
        TabMod[  i ][ T_REC ] = ModRec;
        TabMod[  T_REC ][ i ] = ModRec;
        TabPow[  i ][ T_REC ] = PowRec;
        TabPow[  T_REC ][ i ] = PowRec;
        TabEq[   i ][ T_REC ] = EqRec;
        TabEq[   T_REC ][ i ] = EqRec;
        TabLt[   i ][ T_REC ] = LtRec;
        TabLt[   T_REC ][ i ] = LtRec;
        TabComm[ i ][ T_REC ] = CommRec;
        TabComm[ T_REC ][ i ] = CommRec;
    }

    /* install the internal functions                                      */
    InstIntFunc( "IsRec",       FunIsRec       );
    InstIntFunc( "RecFields",   FunRecFields   );
}
