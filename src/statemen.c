/****************************************************************************
**
*A  statemen.c                  GAP source                   Martin Schoenert
**
*H  @(#)$Id: statemen.c,v 1.1.1.1 1996/12/11 12:43:59 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This module contains the functions for executing the various  statements.
**  Assignments are dealed with in 'eval.c' and functions  are  in  'func.c'.
**
*H  $Log: statemen.c,v $
*H  Revision 1.1.1.1  1996/12/11 12:43:59  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.8  1993/03/01  20:38:06  martin
*H  fixed the printing of 'repeat'-loops
*H
*H  Revision 3.7  1993/02/04  10:51:10  martin
*H  changed to new use new list interface
*H
*H  Revision 3.6  1991/04/30  16:12:47  martin
*H  initial revision under RCS
*H
*H  Revision 3.5  1991/01/17  12:00:00  martin
*H  improved 'for' loop for range constants
*H
*H  Revision 3.4  1990/12/20  12:00:00  martin
*H  added the boolean list package
*H
*H  Revision 3.3  1990/12/19  12:00:00  martin
*H  improved the list like objects package interface
*H
*H  Revision 3.2  1990/12/06  12:00:00  martin
*H  added yet another list package
*H
*H  Revision 3.1  1990/11/20  12:00:00  martin
*H  added new list package
*H
*H  Revision 3.0  1990/10/02  12:00:00  martin
*H  added 'quit'
*H
*/

#include        "system.h"              /* system dependent functions      */
#include        "gasman.h"              /* dynamic storage manager         */
#include        "scanner.h"             /* reading of tokens and printing  */
#include        "eval.h"                /* evaluator main dispatcher       */
#include        "integer.h"             /* arbitrary size integers         */

#include        "list.h"                /* generic list package            */

#include        "statemen.h"            /* definition part of this module  */


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
char *          StrStat;
TypHandle       HdStat;


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
**
**  A statement sequence with <n> statements is represented by a bag with <n>
**  handles, the first is the handle of the first <statement>, the second  is
**  the handle of the second <statement>, and so on.
*/
TypHandle       EvStatseq ( hdSSeq )
    TypHandle           hdSSeq;
{
    TypHandle           hdRes;
    unsigned long       k;

    /* execute the <statement> one after the other                         */
    for ( k = 0; k < SIZE(hdSSeq)/SIZE_HD; ++k ) {
        StrStat = "";  HdStat = PTR(hdSSeq)[k];
        hdRes = EVAL( HdStat );
        if ( hdRes == HdReturn )  return hdRes;
    }

    return HdVoid;
}


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
**
**  An 'if' statement is represented by  a bag where the  first handle is the
**  handle of the  <condition> following the  'if', the second  handle is the
**  handle of the corresponding <statements>,  the third handle is the handle
**  of the <condition> following the  first 'elif', the  fourth handle is the
**  handle  of  the  corresponding <statements>,  and   so  on.  If  the 'if'
**  statement  has no 'else'  part  the bag has   an even number of  handles,
**  otherwise the number of handles is odd and the last  handle is the handle
**  of the <statements> following the 'else'.
*/
TypHandle       EvIf ( hdIf )
    TypHandle           hdIf;
{
    TypHandle           hdRes,  hdSSeq;
    unsigned long       i,  k;

    /* handle the 'if'/'elif' branches in order                            */
    for ( i = 0; i < SIZE(hdIf) / (2*SIZE_HD); ++i ) {

        /* evaluate the <condition>                                        */
        if ( i == 0 ) StrStat = "if ";  else StrStat = "elif ";
        HdStat = PTR(hdIf)[2*i];
        hdRes = EVAL( HdStat );
        while ( hdRes != HdTrue && hdRes != HdFalse )
          hdRes=Error("if: <expr> must evaluate to 'true' or 'false'",0L,0L);

        /* if 'true', execute the <statements> and terminate               */
        if ( hdRes == HdTrue ) {
            hdSSeq = PTR(hdIf)[2*i+1];
            if ( TYPE(hdSSeq) == T_STATSEQ ) {
                for ( k = 0; k < SIZE(hdSSeq)/SIZE_HD; ++k ) {
                    StrStat = "";  HdStat = PTR(hdSSeq)[k];
                    hdRes = EVAL( HdStat );
                    if ( hdRes == HdReturn )  return hdRes;
                }
            }
            else {
                StrStat = "";  HdStat = hdSSeq;
                hdRes = EVAL( HdStat );
                if ( hdRes == HdReturn )  return hdRes;
            }
            return HdVoid;
        }
    }

    /* if present execute the 'else' <statements> and return               */
    if ( SIZE(hdIf) % (2*SIZE_HD) != 0 ) {
        hdSSeq = PTR(hdIf)[SIZE(hdIf)/SIZE_HD-1];
        if ( TYPE(hdSSeq) == T_STATSEQ ) {
            for ( k = 0; k < SIZE(hdSSeq)/SIZE_HD; ++k ) {
                StrStat = "";  HdStat = PTR(hdSSeq)[k];
                hdRes = EVAL( HdStat );
                if ( hdRes == HdReturn )  return hdRes;
            }
        }
        else {
            StrStat = "";  HdStat = hdSSeq;
            hdRes = EVAL( HdStat );
            if ( hdRes == HdReturn )  return hdRes;
        }
        return HdVoid;
    }

    return HdVoid;
}


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
**
**  The 'for' loop  is  represented by a  bag  with three handles, the  first
**  handle is the handle  of the <variable>,  the second handle is the handle
**  of the <list> and the third is the handle of the <statements>.
*/
TypHandle       EvFor ( hdFor )
    TypHandle           hdFor;
{
    TypHandle           hdList,  hdRes,  hdVar,  hdSSeq,  hdElm;
    long                j;
    unsigned long       i,  k;

    /* first evaluate the <list> we are to loop over                       */
    hdVar  = PTR(hdFor)[0];
    hdSSeq = PTR(hdFor)[2];
    StrStat = "for <var>  in ";  HdStat = PTR(hdFor)[1];

    /* special case for a range that appear as constant in the text        */
    /*N 1992/12/16 martin handle general range literals                    */
    if ( TYPE(HdStat) == T_MAKERANGE && SIZE(HdStat) == 2*SIZE_HD ) {

        /* get the low and the high value of the range                     */
        hdList = HdStat;
        hdElm = EVAL( PTR(hdList)[0] );
        while ( TYPE(hdElm) != T_INT )
            hdElm = Error("Range: <low> must be an integer",0L,0L);
        hdList = EVAL( PTR(hdList)[1] );
        while ( TYPE(hdList) != T_INT )
            hdList = Error("Range: <high> must be an integer",0L,0L);

        /* loop over the range                                             */
        for ( j = HD_TO_INT(hdElm); j <= HD_TO_INT(hdList); j++ ) {

            /* assign the element of the range to the variable             */
            PTR(hdVar)[0] = INT_TO_HD( j );

            /* now execute the <statements>                                */
            EnterKernel();
            if ( TYPE(hdSSeq) == T_STATSEQ ) {
                for ( k = 0; k < SIZE(hdSSeq)/SIZE_HD; ++k ) {
                    StrStat = "";  HdStat = PTR(hdSSeq)[k];
                    hdRes = EVAL( HdStat );
                    if ( hdRes == HdReturn ) {
                        ExitKernel( hdRes );
                        return hdRes;
                    }
                }
            }
            else {
                StrStat = "";  HdStat = hdSSeq;
                hdRes = EVAL( HdStat );
                if ( hdRes == HdReturn ) {
                    ExitKernel( hdRes );
                    return hdRes;
                }
            }

            /* give the user the chance to interrupt this loop             */
            StrStat = "for ";  HdStat = hdVar;
            if ( SyIsIntr() )  Error("user interrupt",0L,0L);
            ExitKernel( (TypHandle)0 );

        }

    }

    /* general case                                                        */
    else {

        /* evaluate the list                                               */
        hdList = EVAL( HdStat );
        while ( ! IS_LIST( hdList ) )
            hdList = Error("for: <list> must evaluate to a list",0L,0L);

        /* protect <list> from being removed by the garbage collection     */
        EnterKernel();  ExitKernel( hdList );

        /* loop over all elements in the list                              */
        /* note that the type of the list may dynamically change in loop   */
        for ( i = 1; 1; ++i ) {

            /* get the <i>th element, break if we have reached the end     */
            if ( LEN_LIST(hdList) < i )  break;
            hdElm = ELMF_LIST( hdList, i );
            if ( hdElm == 0 )  continue;
            PTR(hdVar)[0] = hdElm;

            /* now execute the <statements>                                */
            EnterKernel();
            if ( TYPE(hdSSeq) == T_STATSEQ ) {
                for ( k = 0; k < SIZE(hdSSeq)/SIZE_HD; ++k ) {
                    StrStat = "";  HdStat = PTR(hdSSeq)[k];
                    hdRes = EVAL( HdStat );
                    if ( hdRes == HdReturn ) {
                        ExitKernel( hdRes );
                        return hdRes;
                    }
                }
            }
            else {
                StrStat = "";  HdStat = hdSSeq;
                hdRes = EVAL( HdStat );
                if ( hdRes == HdReturn ) {
                    ExitKernel( hdRes );
                    return hdRes;
                }
            }

            /* give the user the chance to interrupt this loop             */
            StrStat = "for ";  HdStat = hdVar;
            if ( SyIsIntr() )  Error("user interrupt",0L,0L);
            ExitKernel( (TypHandle)0 );

        }

    }

    /* and thats it                                                        */
    return HdVoid;
}


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
**
**  The 'while' loop is represented by a bag  with  two  handles,  the  first
**  handle is the handle of the <condition>, the second handle is the  handle
**  of the <statements>.
*/
TypHandle       EvWhile ( hdWhile )
    TypHandle           hdWhile;
{
    TypHandle           hdRes,  hdCond,  hdSSeq;
    unsigned long       k;

    /* get the handles                                                     */
    hdCond = PTR(hdWhile)[0];
    hdSSeq = PTR(hdWhile)[1];

    /* evaluate the <condition> for the first iteration                    */
    StrStat = "while ";  HdStat = hdCond;
    hdRes = EVAL( hdCond );
    while ( hdRes != HdTrue && hdRes != HdFalse )
      hdRes = Error("while: <expr> must evalate to 'true' or 'false'",0L,0L);
    if ( SyIsIntr() )  Error("user interrupt",0L,0L);

    /* while <condition>  do                                               */
    while ( hdRes == HdTrue ) {

        /* execute the <statements>                                        */
        EnterKernel();
        if ( TYPE(hdSSeq) == T_STATSEQ ) {
            for ( k = 0; k < SIZE(hdSSeq)/SIZE_HD; ++k ) {
                StrStat = "";  HdStat = PTR(hdSSeq)[k];
                hdRes = EVAL( HdStat );
                if ( hdRes == HdReturn ) {
                    ExitKernel( hdRes );
                    return hdRes;
                }
            }
        }
        else {
            StrStat = "";  HdStat = hdSSeq;
            hdRes = EVAL( HdStat );
            if ( hdRes == HdReturn ) {
                ExitKernel( hdRes );
                return hdRes;
            }
        }

        /* evaluate the <condition> for the next iteration                 */
        StrStat = "while ";  HdStat = hdCond;
        hdRes = EVAL( hdCond );
        while ( hdRes != HdTrue && hdRes != HdFalse )
            hdRes=Error("while: <expr> must evaluate to 'true' or 'false'",
                        0L,0L);
        if ( SyIsIntr() )  Error("user interrupt",0L,0L);
        ExitKernel( (TypHandle)0 );

    }

    return HdVoid;
}


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
**
**  The 'repeat' loop is represented by a bag  with two  handles,  the  first
**  handle is the handle of the <condition>, the second handle is the  handle
**  of the <statements>.
*/
TypHandle       EvRepeat ( hdRep )
    TypHandle           hdRep;
{
    TypHandle           hdRes, hdCond, hdSSeq;
    unsigned long       k;

    /* get the handles                                                     */
    hdCond = PTR(hdRep)[0];
    hdSSeq = PTR(hdRep)[1];

    /* repeat the <statements> until the <condition> is 'true'             */
    do {

        /* execute the <statements>                                        */
        EnterKernel();
        if ( TYPE(hdSSeq) == T_STATSEQ ) {
            for ( k = 0; k < SIZE(hdSSeq)/SIZE_HD; ++k ) {
                StrStat = "";  HdStat = PTR(hdSSeq)[k];
                hdRes = EVAL( HdStat );
                if ( hdRes == HdReturn ) {
                    ExitKernel( hdRes );
                    return hdRes;
                }
            }
        }
        else {
            StrStat = "";  HdStat = hdSSeq;
            hdRes = EVAL( HdStat );
            if ( hdRes == HdReturn ) {
                ExitKernel( hdRes );
                return hdRes;
            }
        }

        /* evaluate the <condition>                                        */
        StrStat = "until ";  HdStat = hdCond;
        hdRes = EVAL( hdCond );
        while ( hdRes != HdTrue && hdRes != HdFalse )
            hdRes=Error("repeat: <expr> must evaluate to 'true' or 'false'",
                        0L,0L);
        if ( SyIsIntr() )  Error("user interrupt",0L,0L);
        ExitKernel( (TypHandle)0 );

    } while ( hdRes != HdTrue );

    return HdVoid;
}


/****************************************************************************
**
*F  PrStatseq( <hdSSeq> ) . . . . . . . . . . . .  print a statement sequence
**
**  'PrStatseq' prints the statement sequence <hdSSeq>.
**
**  A linebreak is forced after each <statement> except the last one.
*/
void            PrStatseq ( hdSSeq )
    TypHandle           hdSSeq;
{
    unsigned long       k;

    /* print the <statements> one after another, separated by linebreaks   */
    for ( k = 0; k < SIZE(hdSSeq)/SIZE_HD; ++k ) {
        Print( PTR(hdSSeq)[k] );
        if ( k < SIZE(hdSSeq)/SIZE_HD-1 )
            Pr(";\n",0L,0L);
    }
}


/****************************************************************************
**
*F  PrIf( <hdIf> )  . . . . . . . . . . . . . . . . . . print an if statement
**
**  'PrIf' prints the 'if' statement <hdIf>.
**
**  A Linebreak is forced after the 'then'  and  <statements>.  If  necessary
**  it is preferred immediately before the 'then'.
*/
void            PrIf ( hdIf )
    TypHandle           hdIf;
{
    unsigned long       i;

    /* print the 'if' and 'elif' parts                                     */
    for ( i = 0; i < SIZE(hdIf)/SIZE_HD/2; ++i ) {
        if ( i == 0 ) Pr("if%4> ",0L,0L);  else Pr("elif%4> ",0L,0L);
        Print( PTR(hdIf)[2*i] );
        Pr("%2<  then%2>\n",0L,0L);
        Print( PTR(hdIf)[2*i+1] );
        Pr(";%4<\n",0L,0L);
    }

    /* print the 'else' part if it exists                                  */
    if ( SIZE(hdIf)/SIZE_HD % 2 != 0 ) {
        Pr("else%4>\n",0L,0L);
        Print( PTR(hdIf)[ SIZE(hdIf)/SIZE_HD -1 ] );
        Pr(";%4<\n",0L,0L);
    }

    /* print the 'fi'                                                      */
    Pr("fi",0L,0L);
}


/****************************************************************************
**
*F  PrFor( <hdFor> )  . . . . . . . . . . . . . . . . . . .  print a for loop
**
**  'PrFor' prints the 'for' loop <hdFor>.
**
**  A linebreak is forced after the 'do' and the <statements>.  If  necesarry
**  it is preferred immediately before the 'in'.
*/
void            PrFor ( hdFor )
    TypHandle           hdFor;
{
    Pr("for%4> ",0L,0L);       Print( PTR(hdFor)[0] );
    Pr("%2<  in%2> ",0L,0L);   Print( PTR(hdFor)[1] );
    Pr("%2<  do%2>\n",0L,0L);  Print( PTR(hdFor)[2] );
    Pr(";%4<\nod",0L,0L);
}


/****************************************************************************
**
*F  PrWhile( <hdWhile> )  . . . . . . . . . . . . . . . .  print a while loop
**
**  'PrWhile' prints the 'while' loop <hdWhile>.
**
**  A linebreak is forced after the 'do' and the <statements>.  If  necessary
**  it is preferred immediately before the 'do'.
*/
void            PrWhile ( hdWhile )
    TypHandle           hdWhile;
{
    Pr("while%4> ",0L,0L);     Print( PTR(hdWhile)[0] );
    Pr("%2<  do%2>\n",0L,0L);  Print( PTR(hdWhile)[1] );
    Pr(";%4<\nod",0L,0L);
}


/****************************************************************************
**
*F  PrRepeat( <hdRep> ) . . . . . . . . . . . . . . . . . print a repeat loop
**
**  'PrRepeat' prints the 'repeat until' loop <hdRep>.
**
**  A linebreak is forced after the 'repeat' and the <statements>.
*/
void            PrRepeat ( hdRep )
    TypHandle           hdRep;
{
    Pr("repeat%4>\n",0L,0L);
    Print( PTR(hdRep)[1] );
    Pr(";%4<\nuntil%2> ",0L,0L);
    Print( PTR(hdRep)[0] );
    Pr("%2<",0L,0L);
}


/****************************************************************************
**
*F  InitStat()  . . . . . . . . . . . . . . . initialize the statement module
**
**  Is called from 'InitEval' to initialize the statement evaluation  module.
*/
void            InitStat ()
{
    InstEvFunc( T_STATSEQ,  EvStatseq );
    InstEvFunc( T_IF,       EvIf      );
    InstEvFunc( T_FOR,      EvFor     );
    InstEvFunc( T_WHILE,    EvWhile   );
    InstEvFunc( T_REPEAT,   EvRepeat  );

    InstPrFunc( T_STATSEQ,  PrStatseq );
    InstPrFunc( T_IF,       PrIf      );
    InstPrFunc( T_FOR,      PrFor     );
    InstPrFunc( T_WHILE,    PrWhile   );
    InstPrFunc( T_REPEAT,   PrRepeat  );
}
