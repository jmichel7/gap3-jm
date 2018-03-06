/****************************************************************************
**
*A  function.c                  GAP source                   Martin Schoenert
**
*H  @(#)$Id: function.c,v 1.1.1.1 1996/12/11 12:44:02 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This package contains the functions  that  mainly  deal  with  functions.
**
*H  $Log: function.c,v $
*H  Revision 1.1.1.1  1996/12/11 12:44:02  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 3.12  1994/01/31  11:54:28  fceller
*H  added 'Apply'Func'
*H
*H  Revision 3.11  1993/12/20  13:23:47  mschoene
*H  fixed 'EvFunccall' to check for void arguments in arg functions
*H
*H  Revision 3.10  1993/03/19  17:30:30  martin
*H  changed to of stack frames to 'T_EXEC'
*H
*H  Revision 3.9  1993/02/04  10:51:10  martin
*H  changed to use 'plist' interface
*H
*H  Revision 3.8  1992/04/03  16:08:23  martin
*H  renamed 'Trace' to 'TraceFunc'
*H
*H  Revision 3.7  1992/01/06  14:25:43  martin
*H  changed 'EvFuncCall', '~' must not be changed to early
*H
*H  Revision 3.6  1992/01/06  14:03:27  martin
*H  changed the implementation of '~' slightly
*H
*H  Revision 3.5  1991/06/17  08:05:09  martin
*H  fixed 'PrReturn' to print return statements without expressions
*H
*H  Revision 3.4  1991/04/30  16:12:19  martin
*H  initial revision under RCS
*H
*H  Revision 3.3  1990/12/07  12:00:00  martin
*H  added prototypes for function tables
*H
*H  Revision 3.2  1990/10/02  12:00:00  martin
*H  added 'quit'
*H
*H  Revision 3.1  1990/10/02  12:00:00  martin
*H  fixed 'Trace' for functions without local vars
*H
*H  Revision 3.0  1990/08/31  12:00:00  martin
*H  changed profiling of local functions
*/

#include        "system.h"              /* system dependent functions      */
#include        "gasman.h"              /* dynamic storage manager         */
#include        "scanner.h"             /* reading of tokens and printing  */
#include        "eval.h"                /* evaluator main dispatcher       */
#include        "integer.h"             /* arbitrary size integers         */

#include        "list.h"                /* generic list package            */
#include        "plist.h"               /* plain list package              */
#include        "record.h"              /* 'HdTilde'                       */
#include        "statemen.h"            /* 'HdStat', 'StrStat'             */

#include        "function.h"            /* declaration file of the package */


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
TypHandle       HdExec;


/****************************************************************************
**
*V  IsProfiling . . . . . . . . . . . . . .  is 1 if profiling is switched on
*V  HdTimes . . .  handle of the list that contains the profiling information
*V  Timesum . .  total time spent in all functions that have completed so far
**
**  'IsProfiling' is 1 if profiling is enabled and 0 otherwise.
**
**  'HdTimes'  is  the  handle  of  the  list  that  contains  the  profiling
**  information.  This list contains for every function the following entries
**
**      <function>      handle of the function bag (or body of the function)
**      <name>          handle of the first entry in the function call
**      <count>         number of times this function was called
**      <time>          time spent in this function without its children
**      <total>         time spent in this function with its childer
**
**  'Timesum' is the total time spent in all functions that have completed so
**  far.  When a function  is  called  the  current  value  of  'Timesum'  is
**  remembered.  When the function completes 'Timesum' - <old>  is  the  time
**  spent in all childeren of this function.  If this is subtracted from  the
**  total time spent in this function we have the time spent in this function
**  without its children.
*/
long            IsProfiling;
TypHandle       HdTimes;
unsigned long   Timesum;


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
void            ChangeEnv ( hdEnv )
    TypHandle           hdEnv;
{
    TypHandle  hdDo, hdComm, hdTmp, hdUndo;
    TypHandle  * ptUndo,  * ptDef,  * ptDo;
    short      nr,  i;

    /* first walk down the new chain until we find a active exec bag       */
    /* we reverse the links, so we can walk back later                     */
    hdDo   = 0;
    hdComm = hdEnv;
    while ( hdComm != 0 && PTR(hdComm)[1] != HdTrue ) {
        hdTmp          = PTR(hdComm)[0];
        PTR(hdComm)[0] = hdDo;
        hdDo           = hdComm;
        hdComm         = hdTmp;
    }

    /* then we undo all changes from the topmost down to the common exec   */
    hdUndo = HdExec;
    while ( hdUndo != hdComm ) {
        ptUndo = PTR(hdUndo);
        ptDef  = PTR(ptUndo[2]);
        nr = (SIZE(ptUndo[2])-2*sizeof(short)-2*SIZE_HD)/SIZE_HD;
        for ( i = 1; i <= nr; ++i ) {
            hdTmp = ptUndo[i+4];
            ptUndo[i+4] = PTR(ptDef[i])[0];
            PTR(ptDef[i])[0] = hdTmp;
        }
        ptUndo[1] = HdFalse;
        hdUndo    = ptUndo[0];
    }

    /* then we redo all changes from the common up to the new topmost exec */
    while ( hdDo != 0 ) {
        ptDo   = PTR(hdDo);
        ptDef  = PTR(ptDo[2]);
        nr = (SIZE(ptDo[2])-2*sizeof(short)-2*SIZE_HD)/SIZE_HD;
        for ( i = 1; i <= nr; ++i ) {
            hdTmp = ptDo[i+4];
            ptDo[i+4] = PTR(ptDef[i])[0];
            PTR(ptDef[i])[0] = hdTmp;
        }
        ptDo[1] = HdTrue;
        hdTmp   = ptDo[0];
        ptDo[0] = hdComm;
        hdComm  = hdDo;
        hdDo    = hdTmp;
    }

    /* reflect the new environment in HdExec                               */
    HdExec = hdComm;
}


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
**  'EvFunccall' first creates a new execute bag.  Then it evaluates the  new
**  arguments, and puts the values in the execute bag.  Then it saves the old
**  values of the arguments and local variables in the execute bag.  Then  it
**  calles 'ChangeEnv' to copy the new values from the execute bag  into  the
**  variables.  Now the binding is complete, and  'EvFunccall'  executes  the
**  statement sequence.  After that 'EvFunccall' calls 'ChangeEnv'  again  to
**  restore the old values from the execute bag.
*/
TypHandle       EvFunccall ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdDef,  hdExec,  hdRes = 0,  hdOld;
    TypHandle           hdStat,  hd;
    TypHandle           hdTilde;
    short               nrArg,  nrLoc,  i,  trace;
    unsigned long       time = 0,  sime = 0;

    /* remember the old value of HdStat to recover later                   */
    hdStat = HdStat;

    /* Get the handle of the function definition                           */
    hdDef = EVAL( PTR(hdCall)[0] );
    if ( TYPE(hdDef) != T_FUNCTION && TYPE(hdDef) != T_FUNCINT )
        return Error("Function: <function> must be a function",0L,0L);

    /* treat the special case of internal functions                        */
    if ( TYPE(hdDef) == T_FUNCINT ) {
        if ( IsProfiling ) {
            IsProfiling = 1;
            time = SyTime()-Timesum;
            for ( i = 0; i < SIZE(HdTimes)/SIZE_HD; i += 5 ) {
                if ( PTR(HdTimes)[i] == hdDef ) {
                    sime = SyTime() - HD_TO_INT( PTR(HdTimes)[i+4] );
                    break;
                }
            }
            if ( i == SIZE(HdTimes)/SIZE_HD ) {
                sime = SyTime();
            }
        }
        hdRes = (** (TypHandle(**)())PTR(hdDef)) ( hdCall );
        if ( IsProfiling == 1 ) {
            time = SyTime()-Timesum-time;  Timesum += time;
            for ( i = 0; i < SIZE(HdTimes)/SIZE_HD; i += 5 ) {
             if ( PTR(HdTimes)[i] == hdDef ) {
              PTR(HdTimes)[i+2]=INT_TO_HD(HD_TO_INT(PTR(HdTimes)[i+2])+1);
              PTR(HdTimes)[i+3]=INT_TO_HD(HD_TO_INT(PTR(HdTimes)[i+3])+time);
              PTR(HdTimes)[i+4]=INT_TO_HD(SyTime()-sime);
              break;
             }
            }
            if ( i == SIZE(HdTimes)/SIZE_HD ) {
                Resize( HdTimes, SIZE(HdTimes) + 5*SIZE_HD );
                PTR(HdTimes)[i]   = hdDef;
                PTR(HdTimes)[i+1] = PTR(hdCall)[0];
                PTR(HdTimes)[i+2] = INT_TO_HD(1);
                PTR(HdTimes)[i+3] = INT_TO_HD(time);
                PTR(HdTimes)[i+4] = INT_TO_HD(SyTime()-sime);
            }
            while ( 0 < i
                && (long)PTR(HdTimes)[i-2] < (long)PTR(HdTimes)[i+3] ) {
                hd = PTR(HdTimes)[i-5];
                PTR(HdTimes)[i-5] = PTR(HdTimes)[i];
                PTR(HdTimes)[i]   = hd;
                hd = PTR(HdTimes)[i-4];
                PTR(HdTimes)[i-4] = PTR(HdTimes)[i+1];
                PTR(HdTimes)[i+1] = hd;
                hd = PTR(HdTimes)[i-3];
                PTR(HdTimes)[i-3] = PTR(HdTimes)[i+2];
                PTR(HdTimes)[i+2] = hd;
                hd = PTR(HdTimes)[i-2];
                PTR(HdTimes)[i-2] = PTR(HdTimes)[i+3];
                PTR(HdTimes)[i+3] = hd;
                hd = PTR(HdTimes)[i-1];
                PTR(HdTimes)[i-1] = PTR(HdTimes)[i+4];
                PTR(HdTimes)[i+4] = hd;
                i -= 5;
            }
        }
        return( hdRes );
    }

    /* tell Gasman that we are entering the kernel                         */
    EnterKernel();

    /* compute the number of arguments and locals                          */
    trace = 0;
    nrArg = ((short*)((char*)PTR(hdDef) + SIZE(hdDef)))[ -2 ];
    nrLoc = ((short*)((char*)PTR(hdDef) + SIZE(hdDef)))[ -1 ];
    if ( nrArg == -1 ) {
        hdRes = NewBag( T_LIST, SIZE_PLEN_PLIST( SIZE(hdCall)/SIZE_HD-1 ) );
        SET_LEN_PLIST( hdRes, SIZE(hdCall)/SIZE_HD-1 );
        for ( i = 1; i < SIZE(hdCall)/SIZE_HD; i++ ) {
            hd = EVAL( PTR(hdCall)[i] );
            if ( TYPE(hd) == T_VOID )
                hd = Error("illegal void argument",0L,0L);
            SET_ELM_PLIST( hdRes, i, hd );
        }
        nrArg = 1;
        trace = 2;
    }
    else if ( nrArg != SIZE(hdCall) / SIZE_HD - 1 )
        return Error("Function: number of args must be %d",(long)nrArg,0L);

    /* check if the function is to be traced                               */
    if ( nrLoc < 0 ) {
        trace |= 1;  nrLoc = -nrLoc-1;
        Pr("\n%2>",0L,0L);  Print( PTR(hdCall)[0] );  Pr("%<( ",0L,0L);
    }

    /* Now create the new execute bag                                      */
    hdExec = NewBag( T_EXEC, SIZE_HD*(2*nrArg+nrLoc+5) );

    /* enter all relevant information into the execbag                     */
    if ( PTR(hdDef)[nrArg+nrLoc+1] == 0 )  PTR(hdExec)[0] = HdExec;
    else PTR(hdExec)[0] = PTR(hdDef)[nrArg+nrLoc+1];
    PTR(hdExec)[1] = HdFalse;           /* this frame is not yet current   */
    PTR(hdExec)[2] = hdDef;             /* function definition             */
    PTR(hdExec)[3] = hdCall;            /* function call, for debug only   */
    PTR(hdExec)[4] = HdExec;            /* calling environment, dbg only   */

    /* enter the new evaluated arguments in the execbag                    */
    for ( i = 1; i <= nrArg; ++i ) {
        if ( ! (trace & 2) )
            hdRes = EVAL( PTR(hdCall)[i] );
        if ( TYPE(hdRes) == T_VOID )
            hdRes = Error("illegal void argument",0L,0L);
        PTR(hdExec)[i+4] = hdRes;
        PTR(hdExec)[nrArg+nrLoc+i+4] = hdRes;
        if ( trace & 1 ) {
            Pr("%>",0L,0L);  Print( hdRes );
            if ( i < nrArg )  Pr("%<, ",0L,0L);
            else              Pr("%< )",0L,0L);
        }
    }

    /* And now change the environment                                      */
    hdOld = HdExec;
    ChangeEnv( hdExec );

    /* If there are timed functions compute the timing                     */
    if ( IsProfiling ) {
        IsProfiling = 1;
        time = SyTime()-Timesum;
        for ( i = 0; i < SIZE(HdTimes)/SIZE_HD; i += 5 ) {
            if ( PTR(HdTimes)[i] == PTR(hdDef)[0] ) {
                sime = SyTime() - HD_TO_INT( PTR(HdTimes)[i+4] );
                break;
            }
        }
        if ( i == SIZE(HdTimes)/SIZE_HD ) {
            sime = SyTime();
        }
    }
    StrStat = "";  HdStat = PTR(hdDef)[0];

    /* remember the old value of '~' to recover later                      */
    hdTilde = PTR(HdTilde)[0];
    PTR(HdTilde)[0] = 0;

    /* well here's what all is about                                       */
    hdRes =  EVAL( PTR(hdDef)[0] );
    if ( hdRes == HdReturn )
        hdRes = PTR(hdRes)[0];
    else
        hdRes = HdVoid;

    /* If there are timed functions compute the timing                     */
    if ( IsProfiling == 1 ) {
        time = SyTime()-Timesum-time;  Timesum += time;
        for ( i = 0; i < SIZE(HdTimes)/SIZE_HD; i += 5 ) {
            if ( PTR(HdTimes)[i] == PTR(hdDef)[0] ) {
              PTR(HdTimes)[i+2]=INT_TO_HD(HD_TO_INT(PTR(HdTimes)[i+2])+1);
              PTR(HdTimes)[i+3]=INT_TO_HD(HD_TO_INT(PTR(HdTimes)[i+3])+time);
              PTR(HdTimes)[i+4]=INT_TO_HD(SyTime()-sime);
              break;
            }
        }
        if ( i == SIZE(HdTimes)/SIZE_HD ) {
            Resize( HdTimes, SIZE(HdTimes) + 5*SIZE_HD );
            PTR(HdTimes)[i]   = PTR(hdDef)[0];
            PTR(HdTimes)[i+1] = PTR(hdCall)[0];
            PTR(HdTimes)[i+2] = INT_TO_HD(1);
            PTR(HdTimes)[i+3] = INT_TO_HD(time);
            PTR(HdTimes)[i+4] = INT_TO_HD(SyTime()-sime);
        }
        while ( 0 < i
             && (long)PTR(HdTimes)[i-2] < (long)PTR(HdTimes)[i+3] ) {
            hd = PTR(HdTimes)[i-5];
            PTR(HdTimes)[i-5] = PTR(HdTimes)[i];
            PTR(HdTimes)[i]   = hd;
            hd = PTR(HdTimes)[i-4];
            PTR(HdTimes)[i-4] = PTR(HdTimes)[i+1];
            PTR(HdTimes)[i+1] = hd;
            hd = PTR(HdTimes)[i-3];
            PTR(HdTimes)[i-3] = PTR(HdTimes)[i+2];
            PTR(HdTimes)[i+2] = hd;
            hd = PTR(HdTimes)[i-2];
            PTR(HdTimes)[i-2] = PTR(HdTimes)[i+3];
            PTR(HdTimes)[i+3] = hd;
            hd = PTR(HdTimes)[i-1];
            PTR(HdTimes)[i-1] = PTR(HdTimes)[i+4];
            PTR(HdTimes)[i+4] = hd;
            i -= 5;
        }
    }

    /* restore old environment                                             */
    ChangeEnv( hdOld );

    /* If the function is traced, print the return value                   */
    if ( trace & 1 ) {
        Pr("\n%>",0L,0L); Print( PTR(hdCall)[0] );  Pr("%< returns",0L,0L);
        if ( hdRes != HdVoid )  Print( hdRes );
        Pr("%<",0L,0L);
    }

    /* recover the value of HdStat to enable debugging                     */
    HdStat = hdStat;

    /* recover the value of '~'                                            */
    PTR(HdTilde)[0] = hdTilde;

    /* tell Gasman that we exit the kernel again                           */
    ExitKernel(hdRes);
    return hdRes;
}


/****************************************************************************
**
*F  EvFunction( <hdFun> ) . . . . . . . . . . . . . . . . evaluate a function
**
**  'EvFunction' returns the value of the function <hdFun>.  Since  functions
**  are constants and thus selfevaluating it just returns <hdFun>.
*/
TypHandle       EvFunction ( hdDef )
    TypHandle           hdDef;
{
    return hdDef;
}


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
TypHandle       EvMakefunc ( hdFun )
    TypHandle           hdFun;
{
    TypHandle           Result;
    short               nrArg,  nrLoc, i;

    Result = NewBag( T_FUNCTION, SIZE(hdFun) );

    /* copy the info about the number of arguments and locals              */
    nrArg = ((short*)((char*)PTR(hdFun) + SIZE(hdFun)))[ -2 ];
    nrLoc = ((short*)((char*)PTR(hdFun) + SIZE(hdFun)))[ -1 ];
    ((short*)((char*)PTR(Result) + SIZE(Result)))[ -2 ] = nrArg;
    ((short*)((char*)PTR(Result) + SIZE(Result)))[ -1 ] = nrLoc;

    /* now copy the formal arguments and locals                            */
    if ( nrArg == -1 )  nrArg = 1;
    for ( i = 0; i <= nrArg+nrLoc; ++i )
        PTR(Result)[i] = PTR(hdFun)[i];

    /* add the environment, i.e., close the function                       */
    PTR(Result)[nrArg+nrLoc+1] = HdExec;

    /* return the new function                                             */
    return Result;
}


/****************************************************************************
**
*F  EvReturn( <hdRet> ) . . . . . . . . . . . . . evaluate a return-statement
**
**  'EvReturn' executes the return-statement with the handle <hdRet>.
**
**  'EvReturn' evaluates the expression in the return bag and puts the  value
**  in the 'HdReturn' bag.  This bag is then  passed  back  through  all  the
**  statement execution functions, until  it  finally  reaches  'EvFunccall'.
**  'EvFunccall' then returns the value in the 'HdResult' bag.
**
**  Note that a quit statement is implemented as a return bag with the  value
**  'HdReturn' in it.  When 'EvReturn' sees this it does not try to  evaluate
**  it but just puts it into the 'HdReturn' bag.  The rules for  'EvFunccall'
**  now say that the function call will return 'HdReturn', thus it will make
**  its way back to the mail loop.
*/
TypHandle       EvReturn ( hdRet )
    TypHandle           hdRet;
{
    TypHandle           hd;

    if ( PTR(hdRet)[0] == HdReturn )
        hd = HdReturn;
    else if ( PTR(hdRet)[0] == HdVoid )
        hd = HdVoid;
    else
        hd = EVAL( PTR(hdRet)[0] );
    PTR(HdReturn)[0] = hd;
    return HdReturn;
}


/****************************************************************************
**
*F  FunIsFunc( <hdCall> ) . . . . . . . . . . . .  internal function 'IsFunc'
**
**  'IsFunc'  returns 'true' if the object <obj> is  a  function and  'false'
**  otherwise.  May cause an error if <obj> is an unbound variable.
*/
TypHandle       FunIsFunc ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdObj;

    /* evaluate and check the argument                                     */
    if ( SIZE(hdCall) != 2 * SIZE_HD )
        return Error("usage: IsFunc( <obj> )",0L,0L);
    hdObj = EVAL( PTR(hdCall)[1] );
    if ( hdObj == HdVoid )
        return Error("IsFunc: function must return a value",0L,0L);

    /* return 'true' if <obj> is a rational and 'false' otherwise          */
    if ( TYPE(hdObj) == T_FUNCTION || TYPE(hdObj) == T_FUNCINT )
        return HdTrue;
    else
        return HdFalse;
}


/****************************************************************************
**
*F  FunTrace( <hdCall> )  . . . . . . . . . . . . . internal function 'Trace'
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
TypHandle       FunTrace ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdDef;
    short               nrLoc,  i;

    for ( i = 1; i < SIZE(hdCall)/SIZE_HD; ++i ) {
        hdDef = EVAL( PTR(hdCall)[i] );
        if ( TYPE(hdDef) == T_FUNCINT )
            return Error("sorry I can not trace internal function",0L,0L);
        if ( TYPE(hdDef) != T_FUNCTION )
            return Error("usage: Trace( <function>... )",0L,0L);
        nrLoc = ((short*)((char*)PTR(hdDef)+SIZE(hdDef)))[-1];
        if ( 0 <= nrLoc )  nrLoc = -nrLoc-1;
        ((short*)((char*)PTR(hdDef)+SIZE(hdDef)))[-1] = nrLoc;
    }
    return HdVoid;
}


/****************************************************************************
**
*F  FunUntrace( <hdCall> )  . . . . . . . . . . . internal function 'Untrace'
**
**  'FunUntrace' implements the internal function 'Untrace'.
**
**  'Untrace( <function>... )'
**
**  'Untrace' switches of the tracing for the functions passed as  arguments.
*/
TypHandle       FunUntrace ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdDef;
    short               nrLoc, i;

    for ( i = 1; i < SIZE(hdCall)/SIZE_HD; ++i ) {
        hdDef = EVAL( PTR(hdCall)[i] );
        if ( TYPE(hdDef) != T_FUNCTION )
            return Error("usage: Untrace( <function>... )",0L,0L);
        nrLoc = ((short*)((char*)PTR(hdDef)+SIZE(hdDef)))[-1];
        if ( nrLoc < 0 )  nrLoc = -nrLoc-1;
        ((short*)((char*)PTR(hdDef)+SIZE(hdDef)))[-1] = nrLoc;
    }
    return HdVoid;
}


/****************************************************************************
**
*F  FunProfile( <hdCall> )  . . . . . . . . . . . internal function 'Profile'
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
TypHandle       FunProfile ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdArg;
    short               i;
    long                total;

    /* check argument count                                                */
    if ( 2 * SIZE_HD < SIZE(hdCall) ) {
        return Error("usage: Profile( true|false ) or Profile()",0L,0L);
    }

    /* switch profiling on or off                                          */
    else if ( SIZE(hdCall) == 2 * SIZE_HD ) {
        hdArg = EVAL( PTR(hdCall)[1] );
        if ( hdArg == HdTrue ) {
            IsProfiling = 2;
            Resize( HdTimes, 0 * SIZE_HD );
        }
        else if ( hdArg == HdFalse ) {
            IsProfiling = 0;
        }
        else {
            return Error("usage: Profile( true|false ) or Profile()",0L,0L);
        }
    }

    /* print profiling information, this should be formatted much nicer    */
    else {
        total = 0;
        for ( i = 0; i < SIZE(HdTimes)/SIZE_HD; i += 5 )
            total = total + HD_TO_INT( PTR(HdTimes)[i+3] );
        if ( total == 0 )  total = 1;
        Pr(" count    time percent time/call child function\n",0L,0L);
        for ( i = 0; i < SIZE(HdTimes)/SIZE_HD; i += 5 ) {
            Pr("%6d  ", HD_TO_INT( PTR(HdTimes)[i+2] ), 0L );
            Pr("%6d  ", HD_TO_INT( PTR(HdTimes)[i+3] ), 0L );
            Pr("%6d  ", 100 * HD_TO_INT(PTR(HdTimes)[i+3]) / total, 0L );
            Pr("%6d  ", HD_TO_INT( PTR(HdTimes)[i+3] ) /
                        HD_TO_INT( PTR(HdTimes)[i+2] ), 0L );
            Pr("%6d  ", HD_TO_INT( PTR(HdTimes)[i+4] ), 0L );
            Print( PTR(HdTimes)[i+1] );
            Pr("\n",0L,0L);
        }
        Pr("        %6d     100                  TOTAL\n",total-1,0L);
    }

    return HdVoid;
}


/****************************************************************************
**
*F  FunApplyFunc( <hdCall> )  . . . . . . . . . . . . .  internal 'ApplyFunc'
*/
TypHandle FunApplyFunc ( hdCall )
    TypHandle		hdCall;
{
    TypHandle           hdNew;		/* the new function call bag       */
    TypHandle           hdFunc;         /* the function                    */
    TypHandle           hdList;         /* and the list                    */
    long                i;              /* loop                            */

    /* check arguments                                                     */
    if ( SIZE(hdCall) != 3*SIZE_HD )
	return Error( "usage: ApplyFunc( <func>, <list> )", 0L, 0L );
    hdFunc = EVAL(PTR(hdCall)[1]);
    hdList = EVAL(PTR(hdCall)[2]);
    if ( ! IS_DENSE_LIST(hdList) )
	return Error( "<list> must be a dense list", 0L, 0L );

    /* create a new function call bag                                      */
    hdNew = NewBag( T_FUNCCALL, SIZE_HD*(1+LEN_LIST(hdList)) );
    PTR(hdNew)[0] = hdFunc;

    /* copy arguments into it                                              */
    for ( i = LEN_LIST(hdList);  0 < i;  i-- )
	PTR(hdNew)[i] = ELMF_LIST( hdList, i );

    /* evaluate this call                                                  */
    return EVAL(hdNew);
}


/****************************************************************************
**
*F  PrFuncint( <hdFun> )  . . . . . . . . . . . .  print an internal function
**
**  'PrFuncint' prints the internal function with the handle  <hdFun> in  the
**  short form: 'function (...) internal; end'.
*/
/*ARGSUSED*/
void            PrFuncint ( hdFun )
    TypHandle           hdFun;
{
    Pr("%2>function%< %>(...)%< %>internal;%< %>end%2<",0L,0L);
}


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
long            prFull;

void            PrFunction ( hdFun )
    TypHandle           hdFun;
{
    short               nrArg,  nrLoc,  i;

    Pr("%5>function%< ( %>",0L,0L);
    nrArg = ((short*)((char*)PTR(hdFun) + SIZE(hdFun)))[-2];
    if ( nrArg == -1 )  nrArg = 1;
    for ( i = 1; i <= nrArg; ++i ) {
        Print( PTR(hdFun)[i] );
        if ( i != nrArg )  Pr("%<, %>",0L,0L);
    }
    Pr(" %<)",0L,0L);

    if ( prFull == 0 ) {
        Pr(" ...%4< ",0L,0L);
    }
    else {
        Pr("\n",0L,0L);
        nrLoc = ((short*)((char*)PTR(hdFun) + SIZE(hdFun)))[-1];
        if ( nrLoc < 0 )  nrLoc = -nrLoc-1;
        if ( nrLoc >= 1 ) {
            Pr("%>local  ",0L,0L);
            for ( i = 1; i <= nrLoc; ++i ) {
                Print( PTR(hdFun)[i+nrArg] );
                if ( i != nrLoc )  Pr("%<, %>",0L,0L);
            }
            Pr("%<;\n",0L,0L);
        }
        Print( PTR(hdFun)[0] );
        Pr(";%4<\n",0L,0L);
    }

    Pr("end",0L,0L);
}


/****************************************************************************
**
*F  PrintFunction( <hdFun> )  . . . . . . . print a function in the full form
**
**  'PrintFunction' prints the function with the handle <hdFun> in  the  full
**  form, i.e., with the statement sequence.  It is called from 'Print'.
*/
void            PrintFunction ( hdFun )
    TypHandle           hdFun;
{
    prFull = 1L;
    PrFunction( hdFun );
    prFull = 0L;
}


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
void            PrFunccall ( hdCall )
    TypHandle           hdCall;
{
    long                i;

    Pr("%2>",0L,0L);  Print( PTR(hdCall)[0] ); Pr("%<( %>",0L,0L);
    for ( i = 1; i < SIZE(hdCall)/SIZE_HD; ++i ) {
	if ( PTR(hdCall)[i] )
        Print( PTR(hdCall)[i] );
	else
	    Pr( "%>...%<", 0L, 0L );
        if ( i != SIZE(hdCall)/SIZE_HD-1 )
            Pr("%<, %>",0L,0L);
    }
    Pr(" %2<)",0L,0L);
}


/****************************************************************************
**
**  PrReturn( <hdRet> ) . . . . . . . . . . . . . .  print a return statement
**
**  'PrReturn' prints the return statement with the  handle  <hdRet>  in  the
**  usual form 'return;' or 'return <expr>;'.
*/
void            PrReturn ( hdRet )
    TypHandle           hdRet;
{
    if ( PTR(hdRet)[0] == HdReturn ) {
        Pr("quit",0L,0L);
    }
    else if ( PTR(hdRet)[0] == HdVoid ) {
        Pr("return",0L,0L);
    }
    else {
        Pr("%2>return%< %>",0L,0L);
        Print( PTR(hdRet)[0] );
        Pr("%2<",0L,0L);
    }
}


/****************************************************************************
**
*F  InitFunc()  . . . . . . . . . . .  initialize function evaluation package
**
**  'InitFunc' initializes the function evaluation package.
*/
void            InitFunc ()
{
    InstEvFunc( T_FUNCCALL, EvFunccall  );
    InstEvFunc( T_FUNCTION, EvFunction  );
    InstEvFunc( T_FUNCINT,  EvFunction  );
    InstEvFunc( T_MAKEFUNC, EvMakefunc  );
    InstEvFunc( T_RETURN,   EvReturn    );

    InstPrFunc( T_FUNCCALL, PrFunccall  );
    InstPrFunc( T_FUNCTION, PrFunction  );
    InstPrFunc( T_FUNCINT,  PrFuncint   );
    InstPrFunc( T_MAKEFUNC, PrFunction  );
    InstPrFunc( T_RETURN,   PrReturn    );

    InstIntFunc( "IsFunc",      FunIsFunc    );
    InstIntFunc( "TraceFunc",   FunTrace     );
    InstIntFunc( "UntraceFunc", FunUntrace   );
    InstIntFunc( "ApplyFunc",   FunApplyFunc );

    HdTimes = NewBag( T_LIST, 0 );
    InstIntFunc( "Profile", FunProfile );
    HdReturn = NewBag( T_RETURN, SIZE_HD );
}
