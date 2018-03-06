/****************************************************************************
**
*A  gap.c                       GAP source                   Martin Schoenert
**
*H  @(#)$Id: gap.c,v 1.4 1997/01/25 18:28:50 werner Exp $
**
*Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file contains the various read-eval-print loops and  related  stuff.
**
*H  $Log: gap.c,v $
*H  Revision 1.4  1997/01/25 18:28:50  werner
*H  Fixed a couple of revision comments.
*H
*H  Revision 1.3  1997/01/10 14:13:02  sal
*H  Revised the tilde bug fix. ~ is now unset before entry to the break
*H  loop but kept locally and restored if the break loop is left by
*H  return. I believe that the local copy is safe from garbage collection
*H  due to an enclosing EnterKernel() ExitKernel() pair, but I am not sure
*H  that I have understood this correctly.
*H
*H  	SL 10/1/97
*H
*H  Revision 1.2  1997/01/10  12:08:13  sal
*H  Reset global variable ~ just before returning from an error. This
*H  fixes the bug Rob Wainwright found.  It would be a little safer to
*H  reset it BEFORE entering a break loop rather than after leaving it
*H  (since it contains an incomplete data structure) but it might be
*H  useful in debugging.
*H 
*H  	SL
*H 
*H  Revision 1.1.1.1  1996/12/11  12:43:53  werner
*H  Preparing 3.4.4 for release
*H 
*H  Revision 3.28  1994/06/10  16:22:28  mschoene
*H  changed 'VERSRC' for 3.4
*H
*H  Revision 3.27  1993/11/08  23:12:40  martin
*H  corrected the date in 'VERSRC'
*H
*H  Revision 3.26  1993/10/15  09:21:46  martin
*H  changed 'SyLinelength' to 'SyNrRows', added 'VERSRC' and 'VERSYS'
*H
*H  Revision 3.25  1993/05/25  08:20:40  fceller
*H  added 'IsIdentical'
*H
*H  Revision 3.24  1993/05/05  11:10:12  fceller
*H  added 'LogInputTo'
*H
*H  Revision 3.23  1993/03/11  13:00:03  fceller
*H  added package mode
*H
*H  Revision 3.22  1993/02/15  14:39:42  fceller
*H  fixed a typo in 'AppendTo'
*H
*H  Revision 3.21  1993/02/12  11:20:18  martin
*H  strings should print as strings in the main loop
*H
*H  Revision 3.20  1993/02/09  10:22:29  martin
*H  changed 'Print' to print empty string literals empty
*H
*H  Revision 3.19  1993/02/04  14:07:46  martin
*H  changed 'Print' to stop printing the empty list as a string
*H
*H  Revision 3.18  1993/02/04  13:06:02  martin
*H  fixed wrong test in 'CoefficientsInt'
*H
*H  Revision 3.17  1993/02/04  10:51:10  martin
*H  changed to use new list interface
*H
*H  Revision 3.16  1992/11/06  08:52:19  fceller
*H  Added 'FunTmpName' and 'SyTmpname'
*H
*H  Revision 3.15  1992/04/02  17:26:05  martin
*H  added the banner
*H
*H  Revision 3.14  1992/04/02  15:21:08  martin
*H  replaced 'Linelength' with 'SizeScreen'
*H
*H  Revision 3.13  1992/03/04  18:26:23  goetz
*H  restricted 'Linelength' to 256.
*H
*H  Revision 3.12  1992/03/01  20:48:33  jmnich
*H  added statistic functions
*H
*H  Revision 3.11  1992/01/02  13:12:28  martin
*H  changed 'Backtrace' to handle binary operators
*H
*H  Revision 3.10  1991/09/24  14:42:55  fceller
*H  'Coefficients' is now 'CoefficientsInt'.
*H
*H  Revision 3.9  1991/08/19  12:49:34  fceller
*H  Added Prototype for 'FunIgnore'.
*H
*H  Revision 3.8  1991/07/29  13:00:45  fceller
*H  'FunIgnore' added.
*H
*H  Revision 3.7  1991/04/30  16:12:20  martin
*H  initial revision under RCS
*H
*H  Revision 3.6  1991/01/23  12:00:00  martin
*H  improved 'Pr' to accept field width for strings
*H
*H  Revision 3.5  1991/01/15  12:00:00  martin
*H  added statistics to Gasman
*H
*H  Revision 3.4  1991/01/15  12:00:00  martin
*H  added 'Coefficients' temporary
*H
*H  Revision 3.3  1991/01/14  12:00:00  martin
*H  added expert functions 'TYPE' and 'SIZE'
*H
*H  Revision 3.2  1990/12/07  12:00:00  goetz
*H  added 'FunLinelength'
*H
*H  Revision 3.1  1990/10/02  12:00:00  martin
*H  added 'quit'
*H
*H  Revision 3.0  1990/09/30  12:00:00  martin
*H  changed '-2^2' to '-4'
*H
*/

#include        <setjmp.h>              /* definition of setjmp buffer     */

#include        "system.h"              /* system dependent functions      */
#include        "gasman.h"              /* dynamic storage manager         */
#include        "scanner.h"             /* reading of single tokens        */
#include        "eval.h"                /* evaluator main dispatcher       */
#include        "integer.h"             /* arbitrary size integers         */

#include        "idents.h"              /* 'InitIdents', 'FindIdent'       */
#include        "read.h"                /* 'ReadIt'                        */

#include        "list.h"                /* generic list package            */
#include        "plist.h"               /* 'LEN_PLIST', 'SET_LEN_PLIST',.. */
#include        "string.h"              /* 'IsString', 'PrintString'       */

#include        "statemen.h"            /* 'HdStat', 'StrStat'             */
#include        "function.h"            /* 'HdExec', 'ChangeEnv', 'PrintF' */
#include        "record.h"              /* 'HdCall*', 'HdTilde'            */


/****************************************************************************
**
*V  ErrRet  . . . . . . . . . . . . . . . environment for error return, local
**
**  'ErrRet' is  the  saved  environment  from  the  beginning  of  the  main
**  read-eval-print loop.  If you 'quit;' a break  loop  control  returns  to
**  that place.
*/
jmp_buf         ErrRet;


/****************************************************************************
**
*V  HdLast  . . . . . . . . . . . . . .  handle of the variable 'last', local
*V  HdLast2 . . . . . . . . . . . . . . handle of the variable 'last2', local
*V  HdLast3 . . . . . . . . . . . . . . handle of the variable 'last3', local
**
**  'HdLast' is the handle of  the  variable  'last'.  This  global  variable
**  holds the result of the last evaluation in the main read-eval-print loop.
**  'HdLast2' likewise holds the next to last result, and 'HdLast3' holds the
**  result before that.
*/
TypHandle       HdLast, HdLast2, HdLast3;


/****************************************************************************
**
*V  HdTime  . . . . . . . . . . . . . .  handle of the variable 'time', local
**
**  'HdTime' is the handle of the variable 'time'.
**
**  'time' holds the time in milliseconds that  the  execution  of  the  last
**  statement took.  This variable is set at the end of  the  read-eval-print
**  cycle.
*/
TypHandle       HdTime;


/****************************************************************************
**
*F  main( <argc>, <argv> )  . . . . . . .  main program, read-eval-print loop
**
**  'main' is the entry point of GAP.  The operating sytem transfers  control
**  to this point when GAP is started.  'main' calls 'InitGap' to  initialize
**  everything.  Then 'main' starts the read-eval-print loop, i.e.,  it reads
**  expression, evaluates it and prints the value.  This continues until  the
**  end of the input file.
*/
int  main (int argc, char *argv[])
{
    TypHandle           hd;
    unsigned long       start;
    extern void         InitGap ();
    extern char         * In;

    /* initialize everything                                               */
    InitGap( argc, argv );

    /* set the environment to return in case of an error                   */
    setjmp( ErrRet );

    /* repeat the read-eval-print cycle until end of input                 */
    while ( Symbol != S_EOF ) {

        /* read an expression                                              */
        Prompt = "gap> ";
        EnterKernel();
        NrError = 0;
        hd = ReadIt();

        /* if there we no syntax error evaluate the expression             */
        if ( hd != 0 ) {
            SyIsIntr();
            start = SyTime();
            hd = EVAL( hd );
            if ( hd == HdReturn && PTR(hd)[0] != HdReturn )
                Error("'return' must not be used in main loop",0L,0L);
            else if ( hd == HdReturn ) {
                hd = HdVoid;
                Symbol = S_EOF;
            }
            PTR(HdTime)[0]  = INT_TO_HD( SyTime() - start );

            /* assign the value to 'last' and then print it                */
            if ( TYPE(hd) != T_VOID ) {
                PTR(HdLast3)[0] = PTR(HdLast2)[0];
                PTR(HdLast2)[0] = PTR(HdLast)[0];
                PTR(HdLast)[0]  = hd;
                if ( *In != ';' ) {
                    IsString( hd );
                    Print( hd );
                    Pr("\n",0L,0L);
                }
            }

        }

        ExitKernel( (TypHandle)0 );
    }

    /* exit to the operating system, the return is there to please lint    */
    SyExit( 0 );
    return 0;
}


/****************************************************************************
**
*F  FunBacktrace( <hdCall> )  . . . . . . . . . internal function 'Backtrace'
**
**  'FunBacktrace' implements the internal function 'Backtrace'.
**
**  'Backtrace()' \\
**  'Backtrace( <level> )'
**
**  'Backtrace' can be used inside a break loop to print  a  history  of  the
**  computation.  'Backtrace' prints a list of  all  active  functions,  most
**  recent first, up to maximal <level>  nestings.  If  <level>  is  positive
**  the names of the formal arguments of the  functions  calls  are  printed,
**  otherwise the  values  of  the  actual  arguments  are  printed  instead.
**  <level> default to 5, i.e., calling 'Backtrace'  with  no  argument  will
**  print the 5 most recent functions with the names of the formal arguments.
**
**  When a break loop (see "Break Loops") is entered  'Backtrace'  is  called
**  automatically.
*/
TypHandle FunBacktrace(TypHandle hdCall)
{
    short               level,  nrArg,  nrLoc,  i;
    TypHandle           hdExec, hdDef;

    /* get the value of <level>                                            */
    if ( hdCall == 0 || SIZE(hdCall) == SIZE_HD ) {
        level = 5;
    }
    else if ( SIZE(hdCall) == 2 * SIZE_HD ) {
        hdDef = EVAL( PTR(hdCall)[1] );
        if ( TYPE(hdDef) != T_INT )
            return Error("usage: Backtrace( <level> )",0L,0L);
        else
            level = HD_TO_INT( hdDef );
    }
    else {
        return Error("usage: Backtrace( <level> )",0L,0L);
    }

    /* for <level> frames                                                  */
    for ( hdExec=HdExec; hdExec!=0 && level!=0; hdExec=PTR(hdExec)[4] ) {

        /* if <level> is positive print only the names of the formal args  */
        if ( 0 < level ) {
            if ( PTR(hdExec)[3] == HdCallSum )
                Pr("<rec1> + <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallDiff )
                Pr("<rec1> - <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallProd )
                Pr("<rec1> * <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallQuo )
                Pr("<rec1> / <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallMod )
                Pr("<rec1> mod <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallPow )
                Pr("<rec1> ^ <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallComm )
                Pr("Comm(<rec1>,<rec2>) called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallEq )
                Pr("<rec1> = <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallLt )
                Pr("<rec1> < <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallIn )
                Pr("<elm> in <rec> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallPrint )
                Pr("Print( <rec> ) called from\n",0L,0L);
            else {
                Print( PTR(hdExec)[3] );
                Pr(" called from\n",0L,0L);
            }
            --level;
        }

        /* if <level> is negative print the values of the arguments        */
        else {
            if ( PTR(hdExec)[3] == HdCallSum )
                Pr("<rec1> + <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallDiff )
                Pr("<rec1> - <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallProd )
                Pr("<rec1> * <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallQuo )
                Pr("<rec1> / <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallMod )
                Pr("<rec1> mod <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallPow )
                Pr("<rec1> ^ <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallComm )
                Pr("Comm(<rec1>,<rec2>) called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallEq )
                Pr("<rec1> = <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallLt )
                Pr("<rec1> < <rec2> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallIn )
                Pr("<elm> in <rec> called from\n",0L,0L);
            else if ( PTR(hdExec)[3] == HdCallPrint )
                Pr("Print( <rec> ) called from\n",0L,0L);
            else {
                Print( PTR( PTR(hdExec)[3] )[0] );
                Pr("%>( %>",0L,0L);
                hdDef = EVAL( PTR( PTR(hdExec)[3] )[0] );
                nrArg = ((short*)((char*)PTR(hdDef)+SIZE(hdDef)))[-2];
                if ( nrArg == -1 )  nrArg = 1;
                nrLoc = ((short*)((char*)PTR(hdDef)+SIZE(hdDef)))[-1];
                if ( nrLoc < 0 )  nrLoc = -nrLoc-1;
                for ( i = 1; i <= nrArg; ++i ) {
                    Print( PTR(hdExec)[i+nrArg+nrLoc+4] );
                    if ( i < nrArg )  Pr("%<, %>",0L,0L);
                }
                Pr(" %2<) called from\n",0L,0L);
            }
            ++level;
        }

    }

    /* print the bottom of the function stack                              */
    if ( hdExec == 0 ) {
        Pr("main loop\n",0L,0L);
    }
    else {
        Pr("...\n",0L,0L);
    }

    return HdVoid;
}


/****************************************************************************
**
*F  Error( <msg>, <arg1>, <arg2> )  . . . . . . . . . . . . . . error handler
**
**  'Error' is the GAP kernel error handler.
*/
TypHandle       Error ( msg, arg1, arg2 )
    char                msg [];
    long                arg1, arg2;
{
    TypHandle           hd;
    TypHandle           hdTilde; 
    long                i;
    extern TypHandle    FunPrint();
    extern char         * In;
    extern TypHandle    HdStack;
    extern unsigned long        TopStack;

    /* open the standard error output file [return code of OpenOutPut ignored]*/
    OpenOutput( "*errout*" );

    /* print the error message, special if called from 'FunError'          */
    if ( SyStrcmp( msg, "FunError" ) != 0 ) {
        Pr("Error, ",0L,0L);  Pr( msg, arg1, arg2 );
    }
    else {
        Pr("Error, ",0L,0L);  FunPrint( (TypHandle)arg1 );
    }

    /* print the stack traceback                                           */
    if ( HdExec != 0 ) {
        if ( HdStat != 0 && SyStrcmp( msg, "FunError" ) != 0 ) {
            Pr(" at\n%s", (long)StrStat, 0L );
            Print( HdStat );
            Pr(" ...",0L,0L);
        }
        Pr(" in\n",0L,0L);
        FunBacktrace( (TypHandle)0 );
    }
    else {
        Pr("\n",0L,0L);
    }

    /* we must disable tilde during break loop processing, since it points to
       a broken object. We may need to restore it if the user 'return's from the
       break loop, so we keep a local copy */
    hdTilde = PTR(HdTilde)[0];
    PTR(HdTilde)[0] = 0;
    
    /* if requested enter a break loop                                     */
    if ( HdExec != 0 && OpenInput( "*errin*" ) ) {

        /* first enter all funcdef bags from the exec list onto the        */
        /* stack, so that we can access args and locals in the loop        */
        hd = HdExec;  TopStack = 0;
        while ( hd != 0 && TopStack+1 < SIZE(HdStack)/SIZE_HD ) {
            ++TopStack;
            hd = PTR(hd)[4];
        }
        hd = HdExec;  i = 0;
        while ( hd != 0 && TopStack-i+1 > 0 ) {
            ++i;
            PTR(HdStack)[TopStack-i+1] = PTR(hd)[2];
            hd = PTR(hd)[4];
        }

        /* now enter a read-eval-print loop, just as in main               */
        while ( Symbol != S_EOF ) {

            /* read an expression                                          */
            Prompt = "brk> ";
            EnterKernel();
            NrError = 0;
            hd = ReadIt();

            /* if there we no syntax error evaluate the expression         */
            if ( hd != 0 ) {
                SyIsIntr();
                hd = EVAL( hd );
                if ( hd == HdReturn && PTR(hd)[0] != HdReturn ) {
                    while ( TopStack >= 1 )
                        PTR(HdStack)[TopStack--] = 0;
		    PTR(HdTilde)[0] = hdTilde; /* restore ~ */
                    ExitKernel( hd );
                    CloseInput();CloseOutput(); /* error codes ignored */
                    return PTR(hd)[0];
                }
                else if ( hd == HdReturn ) {
                    hd = HdVoid;
                    Symbol = S_EOF;
                }

                /* assign the value to 'last' and then print it            */
                if ( TYPE(hd) != T_VOID ) {
                    PTR(HdLast)[0]  = hd;
                    if ( *In != ';' ) {
                        Print( hd );
                        Pr("\n",0L,0L);
                    }
                }

            }

            ExitKernel( (TypHandle)0 );
        }

        /* remove function definitions from the stack and close "*errin*"  */
        while ( TopStack >= 1 )
            PTR(HdStack)[TopStack--] = 0;
        CloseInput(); /* error code ignored */
    }

    /* call ExitKernel(2) to clear new handles bag                         */
    ExitKernel( (TypHandle)2 );
    while ( HdExec != 0 )  ChangeEnv( PTR(HdExec)[4] );

    
    /* close "*errout*" and return to the main read-eval-print loop        */
    while ( CloseOutput() ) ;
    while ( CloseInput() ) ;
    longjmp( ErrRet, 1 );
    return 0;                           /* just to please lint ...         */
}


/****************************************************************************
**
*F  FunIgnore( <hdCall> ) . . . . . . . . . . . .  internal function 'Ignore'
**
**  'FunIgnore' implements the internal function 'Ignore'.
**
**  'Ignore( <arg1>, <arg2>, ... )'
**
**  'Ignore' ignores all its arguments,  it does not even evaluate  them.  So
**  for tracing a GAP function,  use a function 'InfoSomething'  which either
**  has value 'Print' and prints its arguments or has value 'Ignore' and does
**  nothing at all.
*/
TypHandle       FunIgnore( hdCall )
    TypHandle       hdCall;
{
    return HdVoid;
}


/****************************************************************************
**
*F  FunError( <hdCall> )  . . . . . . . . . . . . . internal function 'Error'
**
**  'FunError' implements the internal function 'Error'.
**
**  'Error( <arg1>, <arg2>,... )'
**
**  raises an error.
**  ...A lot of bla about errors and break loops...
**
**  'FunError' simply calls the GAP  kernel  function  'Error',  which  knows
**  that it has been called from 'FunError' because the  format  argument  is
**  'FunError'.  'FunError' passes <hdCall> as the first extra argument.
*/
TypHandle       FunError ( hdCall )
    TypHandle           hdCall;
{
    return Error("FunError", (long)hdCall, 0L );
}


/****************************************************************************
**
*F  FunWindowCmd( <hdCall> )  . . . . . . . . . . .  execute a window command
*/
TypHandle	FunWindowCmd ( hdCall )
    TypHandle	    hdCall;
{
    TypHandle       hdStr;
    TypHandle       hdTmp;
    TypHandle       hdCmd;
    TypHandle       hdLst;
    long            len;
    long            n,  m;
    long            i;
    char          * ptr;
    char          * qtr;

    /* check arguments                                                     */
    if ( SIZE(hdCall) != 2*SIZE_HD )
	return Error( "usage: WindowCmd( <cmds> )", 0L, 0L );
    hdCmd = EVAL(PTR(hdCall)[1]);
    if ( !IsList(hdCmd) )
	return Error( "usage: WindowCmd( <cmds> )", 0L, 0L );
    hdTmp = ELM_LIST(hdCmd,1);
    if ( TYPE(hdTmp) != T_STRING )
	return Error( "<cmd> must be a string", 0L, 0L );
    if ( SIZE(hdTmp) != 4 )
	return Error( "<cmd> is not a valid command", 0L, 0L );

    /* compute size needed to store argument string                        */
    len   = 13;
    hdLst = NewBag( T_LIST, (LEN_LIST(hdCmd)+1)*SIZE_HD );
    for ( i = LEN_LIST(hdCmd);  1 < i;  i-- )
    {
	hdTmp = ELM_LIST(hdCmd,i);
	if ( TYPE(hdTmp) != T_INT && ! IsString(hdTmp) )
	    return Error("%d.th argument must be a string or integer",i,0L);
	PTR(hdLst)[i] = hdTmp;
	if ( TYPE(hdTmp) == T_INT )
	    len += 12;
	else
	    len += 5 + 2*SIZE(hdTmp);
    }

    /* convert <hdCall> into an argument string                            */
    hdStr  = NewBag( T_STRING, len + 13 );
    ptr    = (char*) PTR(hdStr);
    *ptr   = '\0';

    /* first the command name                                              */
    SyStrncat( ptr, (char*)PTR(ELM_LIST(hdCmd,1)), 3 );
    ptr += 3;

    /* and at last the arguments                                           */
    for ( i = 2;  i < SIZE(hdLst)/SIZE_HD;  i++ )
    {
	hdTmp = PTR(hdLst)[i];
	if ( TYPE(hdTmp) == T_INT )
	{
	    *ptr++ = 'I';
	    m = HD_TO_INT(hdTmp);
	    for ( m = (m<0)?-m:m;  0 < m;  m /= 10 )
		*ptr++ = (m%10) + '0';
	    if ( HD_TO_INT(hdTmp) < 0 )
		*ptr++ = '-';
	    else
		*ptr++ = '+';
	}
	else
	{
	    *ptr++ = 'S';
	    m = SIZE(hdTmp)-1;
	    for ( n = 7;  0 <= n;  n--, m /= 10 )
		*ptr++ = (m%10) + '0';
	    qtr = (char*) PTR(hdTmp);
	    for ( m = SIZE(hdTmp)-1;  0 < m;  m-- )
		*ptr++ = *qtr++;
	}
    }
    *ptr = 0;

    /* compute correct length of argument string                           */
    qtr = (char*) PTR(hdStr);
    len = (long)(ptr - qtr);

    /* now call the window front end with the argument string              */
    ptr = SyWinCmd( qtr, len );
    len = SyStrlen(ptr);

    /* now convert result back into a list                                 */
    hdLst = NewBag( T_LIST, SIZE_PLEN_PLIST(11) );
    SET_LEN_PLIST( hdLst, 0 );
    i = 1;
    while ( 0 < len )
    {
	if ( *ptr == 'I' )
	{
	    ptr++;
	    for ( n=0,m=1; '0' <= *ptr && *ptr <= '9'; ptr++,m *= 10,len-- )
		n += (*ptr-'0') * m;
	    if ( *ptr++ == '-' )
		n *= -1;
	    len -= 2;
	    AssPlist( hdLst, i, INT_TO_HD(n) );
	}
	else if ( *ptr == 'S' )
	{
	    ptr++;
	    for ( n = 0, m = 7;  0 <= m;  m-- )
		n = n*10 + (ptr[m]-'0');
	    hdTmp = NewBag( T_STRING, n+1 );
	    *(char*)PTR(hdTmp) = '\0';
	    ptr += 8;
	    SyStrncat( (char*)PTR(hdTmp), ptr, n );
	    ptr += n;
	    len -= n+9;
	    AssPlist( hdLst, i, hdTmp );
	}
	else
	    return Error( "unknown return value '%s'", (long)ptr, 0 );
	i++;
    }

    /* if the first entry is one signal an error */
    if ( ELM_LIST(hdLst,1) == INT_TO_HD(1) )
    {
	hdStr = NewBag( T_STRING, 30 );
	SyStrncat( (char*) PTR(hdStr), "window system: ", 15 );
	SET_ELM_PLIST( hdLst, 1, hdStr );
	Resize( hdLst, i*SIZE_HD );
	return Error( "FunError", (long)hdLst, 0L );
    }
    else
    {
	for ( m = 1;  m <= i-2;  m++ )
	    SET_ELM_PLIST( hdLst,m, ELM_LIST(hdLst,m+1) );
	SET_LEN_PLIST( hdLst, i-2 );
	return hdLst;
    }
}


/****************************************************************************
**
*F  FunREAD( <hdCall> ) . . . . . . . . . . . . . .  internal function 'READ'
**
**  'FunREAD' implements the internal function 'READ'.
**
**  'READ( <filename> )'
**
**  'READ' instructs GAP to read from the file with the  name  <filename>. If
**  it is not found or could not be opened for reading  'false'  is returned.
**  If the file is found GAP reads all expressions and statements  from  this
**  file and evaluates respectively executes them and finally returns 'true'.
**  Then GAP continues evaluation or execution of what it was  doing  before.
**  'READ' can be nested, i.e., it is legal to execute a 'READ' function call
**  in a file that is read with 'READ'.
**
**  If a syntax error is found 'READ' continues reading the  next  expression
**  or statement, just  as  GAP  would  in  the  main  read-eval-print  loop.
**  If an evaluation error occurs, 'READ' enters a break loop.  If you 'quit'
**  this break loop, control returns to the  main  read-eval-print  loop  and
**  reading of <filename> terminates.
**
**  Note that this function is a helper function for  'Read',  which  behaves
**  similar, but causes an error if a file is not found.  'READ'  could  also
**  be used for a 'ReadLib' which searches for a file in various directories.
*/
TypHandle       FunREAD ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hd,  hdName,  hdOld;

    /* check the number and type of arguments                              */
    if ( SIZE(hdCall) != 2*SIZE_HD )
        return Error("usage: READ( <filename> )",0L,0L);
    hdName = EVAL( PTR(hdCall)[1] );
    if ( ! IsString(hdName) )
        return Error("usage: READ( <filename> )",0L,0L);

    /* try to open the given file, if the file is not found return 'false' */
    if ( ! OpenInput( (char*)PTR(hdName) ) )
        return HdFalse;

    /* now comes a read-eval-noprint loop, similar to the one in 'main'    */
    hdOld = HdExec;  HdExec = 0;
    while ( Symbol != S_EOF ) {
        EnterKernel();
        hd = ReadIt();
        if ( hd != 0 )  hd = EVAL( hd );
        if ( hd == HdReturn && PTR(hd)[0] != HdReturn )
            return Error("READ: 'return' must not be used here",0L,0L);
        else if ( hd == HdReturn )
            return Error("READ: 'quit' must not be used here",0L,0L);
        ExitKernel( (TypHandle)0 );
    }
    HdExec = hdOld;

    /* close the input file again, and return 'true'                       */
    if ( ! CloseInput() )
        Error("READ: can not close input, this should not happen",0L,0L);
    return HdTrue;
}


/****************************************************************************
**
*F  FunAUTO( <hdCall> ) . . . . . . . . . . . . . .  internal function 'AUTO'
**
**  'FunAUTO' implements the internal function 'AUTO'.
**
**  'AUTO( <expression>, <var1>, <var2>,... )'
**
**  'AUTO' associates the expression <expression> with the variables <var1>,
**  <var2> etc.  Whenever one those variables is evaluated, i.e.,  when  its
**  value is required, <expression> is automatically  evaluated.  This  must
**  assign a new value to the variable, otherwise an error  is  raised.  The
**  new value is then returned.
**
**  Here is an example of the most important special usage of 'AUTO':
**
**  |    AUTO( ReadLib("integer"), Int, Abs, Sign, Maximum, Minimum ); |
**
**  When one of the variables, 'Int', 'Abs', etc., is  evaluated  the  libary
**  file 'integer.g' is automatically read.  This then defines the functions.
**  This makes it possible to load the library function only on demand.
**
**  'AUTO' is a procedure, i.e., does not return a value.
*/
TypHandle       FunAUTO ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdExpr,  hdVar;
    long                i;

    /* check the number of arguments                                       */
    if ( SIZE(hdCall) < 3 * SIZE_HD )
        return Error("usage: AUTO( <expr>, <var>, <var>... )",0L,0L);

    /* get the expression                                                  */
    hdExpr = PTR(hdCall)[1];

    /* for all remaining arguments                                         */
    for ( i = 2; i < SIZE(hdCall)/SIZE_HD; ++i ) {
        hdVar = PTR(hdCall)[i];

        /* check that they are variables                                   */
        if ( TYPE(hdVar) != T_VAR && TYPE(hdVar) != T_VARAUTO )
            return Error("usage: AUTO( <expr>, <var>, <var>... )",0L,0L);

        /* turn them into automatic variables and bind them to <expr>      */
        Retype( hdVar, T_VARAUTO );
        PTR(hdVar)[0] = hdExpr;

    }

    return HdVoid;
}


/****************************************************************************
**
**  FunPrint( <hdCall> )  . . . . . . . . . . . . . internal function 'Print'
**
**  'FunPrint' implements the internal function 'Print'.
**
**  'Print( <obj1>, <obj2>... )'
**
**  'Print' prints the objects <obj1>, <obj2>,  etc.  one  after  the  other.
**  Strings are printed without the double quotes and special characters  are
**  not escaped, e.g., '\n' is printed as <newline>.  This  makes  a  limited
**  amount of formatting possible.  Functions are printed in the  full  form,
**  i.e., with the function body, not in the abbreviated form.
**
**  'Print' is a procedure, i.e., does not return a value.
**
**  Note that an empty string literal '""' prints empty (remember strings are
**  printed without the double quotes), while an empty list  '[]'  prints  as
**  '[ ]'.
**
**      gap> s := "";;  l := [];;  s = l;
**      gap> Print( s, "\n", l, "\n" );
**      
**      [ ]
**
**  To achieve this 'Print' must be able to distinguish between empty  string
**  literals and other empty lists.  For that it relies on  'IsString'  *not*
**  to convert empty lists to type 'T_STRING'.  This is ugly.
*/
TypHandle       FunPrint ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hd;
    long                i;

    /* print all the arguments, take care of strings and functions         */
    for ( i = 1; i < SIZE(hdCall)/SIZE_HD; ++i ) {
        hd = EVAL( PTR(hdCall)[i] );
        if ( IsString( hd ) && TYPE(hd) == T_STRING )  PrintString( hd );
        else if ( TYPE( hd ) == T_MAKEFUNC )           PrintFunction( hd );
        else if ( TYPE( hd ) == T_FUNCTION )           PrintFunction( hd );
        else if ( TYPE( hd ) != T_VOID )               Print( hd );
        else  hd = Error("function must return a value",0L,0L);
    }

    return HdVoid;
}


/****************************************************************************
**
*F  FunPrntTo( <hdCall> ) . . . . . . . . . . . . internal function 'PrintTo'
**
**  'FunPrntTo' implements the internal function 'PrintTo'.  The stupid  name
**  is neccessary to avoid a name conflict with 'FunPrint'.
**
**  'PrintTo( <filename>, <obj1>, <obj2>... )'
**
**  'PrintTo' prints the objects <obj1>, <obj2>, etc. to the  file  with  the
**  name <filename>.
**
**  'PrintTo' works as follows.  It opens the file with the name  <filename>.
**  If the file does not exist it is  created,  otherwise  it  is  truncated.
**  If you do not want to truncate the file use 'AppendTo'  (see "AppendTo").
**  After opening the file 'PrintTo' evaluates  its  arguments  in  turn  and
**  prints the values to  <filename>.  Finally  it  closes  the  file  again.
**  During evaluation of the arguments <filename> is the current output file.
**  This means that output printed with 'Print' during  the  evaluation,  for
**  example to inform the user about the progress, also goes  to  <filename>.
**  To make this feature more useful 'PrintTo' will silently ignore if one of
**  the arguments is a procedure call, i.e., does not return a value.
**
**  'PrintTo' is a procedure, i.e., does not return a value.
**
**  See the note about empty string literals and empty lists in 'Print'.
*/
TypHandle       FunPrntTo ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hd;
    long                i;

    /* check the number and type of the arguments, nothing special         */
    if ( SIZE(hdCall) == SIZE_HD )
        return Error("usage: PrintTo( <file>, <obj>, <obj>... )",0L,0L);
    hd = EVAL( PTR(hdCall)[1] );
    if ( ! IsString(hd) )
        return Error("usage: PrintTo( <file>, <obj>, <obj>... )",0L,0L);

    /* try to open the given output file, raise an error if you can not    */
    if ( OpenOutput( (char*)PTR(hd) ) == 0 )
        return Error("PrintTo: can not open the file for writing",0L,0L);

    /* print all the arguments, take care of strings and functions         */
    for ( i = 2; i < SIZE(hdCall)/SIZE_HD; ++i ) {
        hd = EVAL( PTR(hdCall)[i] );
        if ( IsString( hd ) && TYPE(hd) == T_STRING )  PrintString( hd );
        else if ( TYPE( hd ) == T_MAKEFUNC )           PrintFunction( hd );
        else if ( TYPE( hd ) == T_FUNCTION )           PrintFunction( hd );
        else if ( TYPE( hd ) != T_VOID )               Print( hd );
        else                                           Pr("",0L,0L);
    }

    /* close the output file again, and return nothing                     */
    if ( ! CloseOutput() )
        Error("PrintTo: can not close output, this should not happen",0L,0L);
    return HdVoid;
}


/****************************************************************************
**
*F  FunAppendTo( <hdCall> ) . . . . . . . . . .  internal function 'AppendTo'
**
**  'FunAppendTo' implements the internal function 'AppendTo'.
**
**  'AppendTo( <filename>, <obj1>, <obj2>... )'
**
**  'AppendTo' appends the obejcts <obj1>, <obj2>, etc. to the file with  the
**  name <filename>.  'AppendTo' works like 'PrintTo' (see "PrintTo")  except
**  that it does not truncate the file if it exists.
**
**  'AppendTo' is a procedure, i.e., does not return a value.
**
**  See the note about empty string literals and empty lists in 'Print'.
*/
TypHandle       FunAppendTo ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hd;
    long                i;

    /* check the number and type of the arguments, nothing special         */
    if ( SIZE(hdCall) == SIZE_HD )
        return Error("usage: AppendTo( <file>, <obj>, <obj>... )",0L,0L);
    hd = EVAL( PTR(hdCall)[1] );
    if ( ! IsString(hd) )
        return Error("usage: AppendTo( <file>, <obj>, <obj>... )",0L,0L);

    /* try to open the given output file, raise an error if you can not    */
    if ( OpenAppend( (char*)PTR(hd) ) == 0 )
        return Error("AppendTo: can not open the file for appending",0L,0L);

    /* print all the arguments, take care of strings and functions         */
    for ( i = 2; i < SIZE(hdCall)/SIZE_HD; ++i ) {
        hd = EVAL( PTR(hdCall)[i] );
        if ( IsString( hd ) && TYPE(hd) == T_STRING )  PrintString( hd );
        else if ( TYPE( hd ) == T_MAKEFUNC )           PrintFunction( hd );
        else if ( TYPE( hd ) == T_FUNCTION )           PrintFunction( hd );
        else if ( TYPE( hd ) != T_VOID )               Print( hd );
        else                                           Pr("",0L,0L);
    }

    /* close the output file again, and return nothing                     */
    if ( ! CloseOutput() )
       Error("AppendTo: can not close output, this should not happen",0L,0L);
    return HdVoid;
}


/****************************************************************************
**
*F  FunLogTo( <hdCall> )  . . . . . . . . . . . . . internal function 'LogTo'
**
**  'FunLogTo' implements the internal function 'LogTo'.
**
**  'LogTo( <filename> )' \\
**  'LogTo()'
**
**  'LogTo' instructs GAP to echo all input from the  standard  input  files,
**  '*stdin*' and '*errin*' and all output  to  the  standard  output  files,
**  '*stdout*'  and  '*errout*',  to  the  file  with  the  name  <filename>.
**  The file is created if it does not  exist,  otherwise  it  is  truncated.
**
**  'LogTo' called with no argument closes the current logfile again, so that
**  input   from  '*stdin*'  and  '*errin*'  and  output  to  '*stdout*'  and
**  '*errout*' will no longer be echoed to a file.
*/
TypHandle       FunLogTo ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdName;

    /* 'LogTo()'                                                           */
    if ( SIZE(hdCall) == SIZE_HD ) {
        if ( ! CloseLog() )
            return Error("LogTo: can not close the logfile",0L,0L);
    }

    /* 'LogTo( <filename> )'                                               */
    else if ( SIZE(hdCall) == 2 * SIZE_HD ) {
        hdName = EVAL( PTR(hdCall)[1] );
        if ( ! IsString(hdName) )
            return Error("usage: LogTo() or LogTo( <string> )",0L,0L);
        if ( ! OpenLog( (char*)PTR(hdName) ) )
            return Error("LogTo: can not log to %s",(long)PTR(hdName),0L);
    }

    return HdVoid;
}


/****************************************************************************
**
*F  FunLogInputTo( <hdCall> ) . . . . . . . .  internal function 'LogInputTo'
**
**  'FunLogInputTo' implements the internal function 'LogInputTo'.
**
**  'LogInputTo( <filename> )' \\
**  'LogInputTo()'
**
**  LogInputTo'  instructs  GAP  to echo  all  input from the  standard input
**  files, '*stdin*' and  '*errin*',  to the file  with the  name <filename>.
**  The file is created if it does not exist, otherwise it is truncated.
**
**  'LogInputTo' called with no argument closes the current logfile again, so
**  that input  from '*stdin*' and '*errin*' will  no longer  be echoed  to a
**  file.
*/
TypHandle       FunLogInputTo ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdName;

    /* 'LogInputTo()'                                                      */
    if ( SIZE(hdCall) == SIZE_HD ) {
        if ( ! CloseInputLog() )
            return Error("LogInputTo: can not close the logfile",0L,0L);
    }

    /* 'LogInputTo( <filename> )'                                          */
    else if ( SIZE(hdCall) == 2 * SIZE_HD ) {
        hdName = EVAL( PTR(hdCall)[1] );
        if ( ! IsString(hdName) )
           return Error("usage: LogInputTo() or LogTo( <string> )",0L,0L);
        if ( ! OpenInputLog( (char*)PTR(hdName) ) )
           return Error("LogInputTo: cannot log to %s",(long)PTR(hdName),0L);
    }

    return HdVoid;
}


/****************************************************************************
**
*F  FunReadTest( <hdCall> ) . . . . . . . . . .  internal function 'ReadTest'
**
**  'FunReadTest' implements the internal function 'ReadTest'.
**
**  'ReadTest( <filename> )'
**
**  'ReadTest' instructs GAP to  read test input from the  file with the name
**  <filename>.   If it is  not found or could not  be  opened for reading an
**  error is raised.  If  the  file is found  GAP  reads all expressions  and
**  statements  from this  file  and  evaluates  respectively executes  them.
**  After  that GAP continues  evaluation or  execution of  what it was doing
**  before.  'ReadTest' can be not nested, i.e., it is not legal to execute a
**  'ReadTest' function call in a file that is read with 'ReadTest'.
**
**  Test mode works as follows.  If GAP is about  to  print  a  line  to  the
**  current  output  file  (or  to  be  more precise  to the output file that
**  was current when  'ReadTest' was called) this line  is  compared with the
**  next line from the test input  file, i.e., the  one opened by 'ReadTest'.
**  If this line starts with '#>' and the rest of it  matches the output line
**  the output line is  not printed and the input  comment line is discarded.
**  Otherwise GAP prints the output line and does not discard the input line.
**
**  On the other hand if an input line is encountered on  the test input that
**  starts with '#>' the GAP assumes that this is  an  expected  output  line
**  that did not appear and echoes this line to the current output file.
**
**  The upshot is that  you can write  test files that consist of alternating
**  input and,  as  '#>' test  comment  lines the  expected  output.   If GAP
**  behaves normal and produces the expected  output then nothing is printed.
**  But if something  goes wrong you see  what actually was printed  and what
**  was expected instead.
**
**  As a convention GAP test files should end with a  print  statement  like:
**
**    Print("prime   3.002   06-Jul-90 ",Quo(417000000,time)," GAPstones\n");
**
**  without a matching '#>' comment line.  This tells the user that the  test
**  file completed and also how much time it took.  The  constant  should  be
**  such that a VAX 11/780 gets roughly 1000 GAPstones.
**
**  If a syntax error is found 'ReadTest' continues reading a next expression
**  or statement, just  as  GAP  would  in  the  main  read-eval-print  loop.
**  If an evaluation error occurs, 'ReadTest' enters a break  loop,  but  the
**  input for this break loop is taken from the test input file.
*/
TypHandle       FunReadTest ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hd,  hdName,  hdOld;
    unsigned long       start;
    extern char         * In;

    /* check the number and type of arguments                              */
    if ( SIZE(hdCall) != 2*SIZE_HD )
        return Error("usage: ReadTest( <filename> )",0L,0L);
    hdName = EVAL( PTR(hdCall)[1] );
    if ( ! IsString(hdName) )
        return Error("usage: ReadTest( <filename> )",0L,0L);

    /* try to open the given file, if the file is not found return 'false' */
    if ( ! OpenTest( (char*)PTR(hdName) ) )
        return Error("ReadTest: file '%s' must exist and be readable\n",
                     (long)PTR(hdName), 0L );
    start = SyTime();

    /* now comes a read-eval-print loop, similar to the one in 'main'      */
    hdOld = HdExec;  HdExec = 0;
    while ( Symbol != S_EOF ) {
        EnterKernel();
        hd = ReadIt();
        if ( hd != 0 ) {
            hd = EVAL( hd );
            if ( hd == HdReturn && PTR(hd)[0] != HdReturn )
                return Error("ReadTest: 'return' must not be used",0L,0L);
            else if ( hd == HdReturn )
                return Error("ReadTest: 'quit' must not be used",0L,0L);
            PTR(HdTime)[0]  = INT_TO_HD( SyTime() - start );
            if ( TYPE(hd) != T_VOID ) {
                if ( *In != ';' ) {
                    IsString( hd );
                    Print( hd );
                    Pr("\n",0L,0L);
                }
            }
        }
        ExitKernel( (TypHandle)0 );
    }
    HdExec = hdOld;

    /* close the input file again, and return 'true'                       */
    if ( ! CloseTest() )
        Error("ReadTest: can not close input, this should not happen",0L,0L);
    return HdVoid;
}


/****************************************************************************
**
*F  FunHelp( <hdCall> ) . . . . . . . . . . . . . .  internal function 'Help'
**
**  'FunHelp' implements the internal function 'Help'.
**
**  'Help( <topic> )'
**
**  'Help' prints a section from the on-line documentation about <topic>.
*/
TypHandle       FunHelp ( hdCall )
    TypHandle           hdCall;
{
    return Error("Help: not yet implemented",0L,0L);
}


/****************************************************************************
**
*F  FunExec( <hdCall> ) . . . . . . . . . . . . . .  internal function 'Exec'
**
**  'FunExec' implements the internal function 'Exec'.
**
**  'Exec( <command> )'
**
**  'Exec' passes the string <command> to  the  command  interpreter  of  the
**  operating system.  The precise mechanismen of this is  system  dependent.
**  Also operating system dependent are the possible commands.
**
**  'Exec' is a procedure, i.e., does not return a value.
*/
TypHandle       FunExec ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdCmd;
    if ( SIZE(hdCall) != 2*SIZE_HD )
        return Error("usage: Exec( <command> )",0L,0L);
    hdCmd = EVAL( PTR(hdCall)[1] );
    if ( ! IsString(hdCmd) )
        return Error("usage: Exec( <command> )",0L,0L);
    SyExec( (char*)PTR(hdCmd) );
    return HdVoid;
}


/****************************************************************************
**
*F  FunRuntime( <hdCall> )  . . . . . . . . . . . internal function 'Runtime'
**
**  'FunRuntime' implements the internal function 'Runtime'.
**
**  'Runtime()'
**
**  'Runtime' returns the time spent since the start of GAP in  milliseconds.
**  How much time execution of statements take is of course system dependent.
**  The accuracy of this number is also system dependent.
*/
TypHandle       FunRuntime ( hdCall )
    TypHandle           hdCall;
{
    if ( SIZE(hdCall) != SIZE_HD )
        return Error("usage: Runtime()",0L,0L);
    return INT_TO_HD( SyTime() );
}


/****************************************************************************
**
*F  FunSizeScreen( <hdCall> ) . . . . . . . .  internal function 'SizeScreen'
**
**  'FunSizeScreen' implements the internal function 'SizeScreen' to  get  or
**  set the actual screen size.
**
**  'SizeScreen()'
**
**  In this form 'ScreeSize' returns the size of the screen as  a  list  with
**  two entries.  The first is the length of each line,  the  second  is  the
**  number of lines.
**
**  'SizeScreen( [ <x>, <y> ] )'
**
**  In this form 'SizeScreen' sets the size of the screen.  <x> is the length
**  of each line, <y> is the number of lines.  Either value may  be  missing,
**  to leave this value unaffected.  Note that those parameters can  also  be
**  set with the command line options '-x <x>' and '-y <y>'.
*/
TypHandle       FunSizeScreen ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdSize;         /* argument and result list        */
    long                len;            /* length of lines on the screen   */
    long                nr;             /* number of lines on the screen   */

    /* check the arguments                                                 */
    if ( SIZE(hdCall) != SIZE_HD && SIZE(hdCall) != 2*SIZE_HD )
        return Error("usage: SizeScreen() or SizeScreen([<x>,<y>])",0L,0L);

    /* no argument is equivalent to the empty list                         */
    if ( SIZE(hdCall) == SIZE_HD ) {
        hdSize = NewBag( T_LIST, SIZE_PLEN_PLIST(0) );
        SET_LEN_PLIST( hdSize, 0 );
    }

    /* otherwise check the argument                                        */
    else {
        hdSize = EVAL( PTR(hdCall)[1] );
        if ( ! IS_LIST(hdSize) || 2 < LEN_LIST(hdSize) )
          return Error("usage: SizeScreen() or SizeScreen([<x>,<y>])",0L,0L);
    }

    /* extract the length                                                  */
    if ( LEN_LIST(hdSize) < 1 || ELMF_LIST(hdSize,1) == 0 ) {
        len = SyNrCols;
    }
    else {
        if ( TYPE( ELMF_LIST(hdSize,1) ) != T_INT )
            return Error("SizeScreen: <x> must be an integer",0L,0L);
        len = HD_TO_INT( ELMF_LIST(hdSize,1) );
        if ( len < 20  )  len = 20;
        if ( 256 < len )  len = 256;
    }

    /* extract the number                                                  */
    if ( LEN_LIST( hdSize ) < 2 || ELMF_LIST(hdSize,2) == 0 ) {
        nr = SyNrRows;
    }
    else {
        if ( TYPE( ELMF_LIST(hdSize,2) ) != T_INT )
            return Error("SizeScreen: <y> must be an integer",0L,0L);
        nr = HD_TO_INT( ELMF_LIST(hdSize,2) );
        if ( nr < 10 )  nr = 10;
    }

    /* set length and number                                               */
    SyNrCols = len;
    SyNrRows = nr;

    /* make and return the size of the screen                              */
    hdSize = NewBag( T_LIST, SIZE_PLEN_PLIST(2) );
    SET_LEN_PLIST( hdSize, 2 );
    SET_ELM_PLIST( hdSize, 1, INT_TO_HD(len) );
    SET_ELM_PLIST( hdSize, 2, INT_TO_HD(nr) );
    return hdSize;
}


/****************************************************************************
**
*F  FunTmpName( <hdCall> )  . . . . . . . . . . . internal function 'TmpName'
**
**  'TmpName()' returns a file names that can safely be used for a temporary
**  file.  It returns 'false' in case of failure.
*/
TypHandle	FunTmpName ( hdCall )
    TypHandle       hdCall;
{
    TypHandle       hdStr;
    char          * str;

    if ( SIZE(hdCall) != SIZE_HD )
	return Error( "usage: TmpName()", 0L, 0L );
    str = SyTmpname();
    if ( str == (char*)0 )
	return HdFalse;
    hdStr = NewBag( T_STRING, SyStrlen(str)+1 );
    *((char*)PTR(hdStr)) = 0;
    SyStrncat( (char*)PTR(hdStr), str, SyStrlen(str) );
    return hdStr;
}


/****************************************************************************
**
*F  FunIsIdentical( <hdCall> )  . . . . . . . internal function 'IsIdentical'
**
**  'FunIsIdentical' implements 'IsIdentical'
*/
TypHandle       FunIsIdentical ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdL;
    TypHandle           hdR;

    if ( SIZE(hdCall) != 3*SIZE_HD )
        return Error( "usage: IsIdentical( <l>, <r> )", 0L, 0L );
    hdL = EVAL( PTR(hdCall)[1] );
    hdR = EVAL( PTR(hdCall)[2] );
    if ( TYPE(hdL) < T_LIST && TYPE(hdR) < T_LIST )
	return EQ( hdL, hdR );
    else if ( TYPE(hdL) < T_LIST || TYPE(hdR) < T_LIST )
	return HdFalse;
    else
	return ( hdL == hdR ) ? HdTrue : HdFalse;
}


/****************************************************************************
**
*F  FunHANDLE( <hdCall> ) . . . . . . . . . . . . .  expert function 'HANDLE'
**
**  'FunHANDLE' implements the internal function 'HANDLE'.
**
**  'HANDLE( <obj> )'
**
**  'HANDLE' returns the handle  of  the  object  <obj>  as  an  integer.  It
**  exists only for debugging purposes and should only be  used  by  experts.
*/
TypHandle       FunHANDLE ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdHD;
    TypHandle           hdObj;

    if ( SIZE(hdCall) != 2*SIZE_HD )
        return Error("usage: HANDLE( <obj> )",0L,0L);
    hdObj = EVAL( PTR(hdCall)[1] );
    hdHD  = INT_TO_HD( (intptr_t)hdObj );
    if ( HD_TO_INT(hdHD) != (intptr_t)hdObj )
        return Error("HANDLE: %d does not fit into 28 bits",(long)hdObj,0L);

    return hdHD;
}


/****************************************************************************
**
*F  FunOBJ( <hdCall> )  . . . . . . . . . . . . . . . . expert function 'OBJ'
**
**  'FunOBJ' implements the internal function 'OBJ'.
**
**  'OBJ( <int> )'
**
**  'OBJ' returns the object with the handle given by the integer  <int>.  It
**  is the inverse function to 'HD'.  Note that passing an integer  to  'OBJ'
**  which is not a valid handle is likely to crash GAP.  Thus  this  function
**  is only there for debugging purposes and should only be used by experts.
*/
TypHandle       FunOBJ ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdObj;
    TypHandle           hdHD;

    if ( SIZE(hdCall) != 2*SIZE_HD )
        return Error("usage: OBJ( <hd> )",0L,0L);
    hdHD = EVAL( PTR(hdCall)[1] );
    if ( TYPE( hdHD ) != T_INT )
        return Error("OBJ: <hd> must be a small integer",0L,0L);
    hdObj = (TypHandle)HD_TO_INT( hdHD );

    return hdObj;
}


/****************************************************************************
**
*F  FunTYPE( <hdCall> ) . . . . . . . . . . . . . . .  expert function 'TYPE'
**
**  'FunTYPE' implements the internal function 'TYPE'.
**
**  'TYPE( <obj> )'
**
**  'TYPE' returns the type of the object <obj> as a string.
*/
TypHandle       FunTYPE ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdType;
    TypHandle           hdObj;

    if ( SIZE(hdCall) != 2*SIZE_HD )
        return Error("usage: TYPE( <obj> )",0L,0L);
    hdObj  = EVAL( PTR(hdCall)[1] );
    if ( hdObj == 0 ) {
        hdType = NewBag( T_STRING, 5 );
        SyStrncat( (char*)PTR(hdType), "null", 4 );
    }
    else {
        hdType = NewBag( T_STRING, SyStrlen(NameType[TYPE(hdObj)])+1 );
        SyStrncat( (char*)PTR(hdType), NameType[TYPE(hdObj)],
                   SyStrlen(NameType[TYPE(hdObj)])+1 );
    }

   return hdType;
}


/****************************************************************************
**
*F  FunSIZE( <hdCall> ) . . . . . . . . . . . . . . .  expert function 'SIZE'
**
**  'FunSIZE' implements the internal function 'SIZE'.
**
**  'SIZE( <obj> )'
**
**  'SIZE' returns the size of the object <obj> including all its subobjects.
**
**  First the  all   bags of  the object  are marked by  'MarkObj' by  adding
**  'T_ILLEGAL'  to their type.  Then 'SizeObj'   only counts marked bags and
**  unmarks them before recursing to subobjects.   This way every bag is only
**  counted once, even if it  appear several times in  the object.  This also
**  helps  to   avoid  infinite recursion if    an  object contains itself as
**  subobject.
*/
void            MarkObj ( hdObj )
    TypHandle           hdObj;
{
    unsigned long       i;

    /* void and small integers do not have a handle structure              */
    if ( hdObj == 0 || TYPE(hdObj) == T_INT )
        return;

    /* do not mark a bag twice                                             */
    if ( T_ILLEGAL <= TYPE(hdObj) )
        return;

    /* mark this bag                                                       */
    hdObj->type += T_ILLEGAL;

    /* mark the subobjects                                                 */
    for ( i = NrHandles( (hdObj->type)-T_ILLEGAL, SIZE(hdObj) ); 0 < i; i-- )
        MarkObj( PTR(hdObj)[i-1] );
}

unsigned long   SizeObj ( hdObj )
    TypHandle           hdObj;
{
    unsigned long       size;
    unsigned long       i;

    /* void and small integers do not use any memory at all                */
    if ( hdObj == 0 || TYPE(hdObj) == T_INT )
        return 0L;

    /* do not count unmarked bags                                          */
    if ( TYPE(hdObj) < T_ILLEGAL )
        return 0L;

    /* unmark this bag                                                     */
    hdObj->type -= T_ILLEGAL;

    /* start with the size of this bag                                     */
    size = SIZE( hdObj );

    /* add the sizes of the subobjects                                     */
    for ( i = NrHandles( TYPE(hdObj), SIZE(hdObj) ); 0 < i; i-- )
        size += SizeObj( PTR(hdObj)[i-1] );

    /* return the size                                                     */
    return size;
}

TypHandle       FunSIZE ( hdCall )
    TypHandle           hdCall;
{
    unsigned long       size;
    TypHandle           hdObj;

    if ( SIZE(hdCall) != 2*SIZE_HD )
        return Error("usage: SIZE( <obj> )",0L,0L);
    hdObj  = EVAL( PTR(hdCall)[1] );
    MarkObj( hdObj );
    size = SizeObj( hdObj );

    return INT_TO_HD( size );
}


/****************************************************************************
**
*F  FunGASMAN( <hdCall> ) . . . . . . . . . . . . .  expert function 'GASMAN'
**
**  'FunGASMAN' implements the internal function 'GASMAN'
**
**  'GASMAN( "display" | "clear" | "collect" | "message" )'
*/
TypHandle       FunGASMAN ( hdCall )
    TypHandle           hdCall;
{
    TypHandle           hdCmd;          /* handle of an argument           */
    unsigned long       i,  k;          /* loop variables                  */

    /* check the argument                                                  */
    if ( SIZE(hdCall) == SIZE_HD )
        return Error(
            "usage: GASMAN( \"display\"|\"clear\"|\"collect\"|\"message\" )",
                     0L,0L);

    /* loop over the arguments                                             */
    for ( i = 1; i < SIZE(hdCall)/SIZE_HD; i++ ) {

        /* evaluate and check the command                                  */
        hdCmd = EVAL( PTR(hdCall)[i] );
        if ( ! IsString(hdCmd) )
           return Error(
            "usage: GASMAN( \"display\"|\"clear\"|\"collect\"|\"message\" )",
                        0L,0L);

        /* if request display the statistics                               */
        if ( SyStrcmp( (char*)PTR(hdCmd), "display" ) == 0 ) {
            Pr("\t\t    type     alive     size     total     size\n",0L,0L);
            for ( k = T_VOID; k < T_ILLEGAL-1; k++ ) {
                Pr("%24s  ", (long)NameType[k], 0L );
                Pr("%8d %8d  ",(long)GasmanStatAlive[k],
                               (long)GasmanStatASize[k]);
                Pr("%8d %8d\n",(long)GasmanStatTotal[k],
                               (long)GasmanStatTSize[k]);
            }
        }

        /* if request display the statistics                               */
        else if ( SyStrcmp( (char*)PTR(hdCmd), "clear" ) == 0 ) {
            for ( k = T_VOID; k < T_ILLEGAL; k++ ) {
                GasmanStatTotal[k] = GasmanStatAlive[k];
                GasmanStatTSize[k] = GasmanStatASize[k];
            }
        }

        /* or collect the garbage                                          */
        else if ( SyStrcmp( (char*)PTR(hdCmd), "collect" ) == 0 ) {
            CollectGarb();
        }

        /* or finally toggle Gasman messages                               */
        else if ( SyStrcmp( (char*)PTR(hdCmd), "message" ) == 0 ) {
            SyGasman = ! SyGasman;
        }

        /* otherwise complain                                              */
        else {
           return Error(
            "usage: GASMAN( \"display\"|\"clear\"|\"collect\"|\"message\" )",
                        0L,0L);
        }
    }

    /* return nothing, this function is a procedure                        */
    return HdVoid;
}


/****************************************************************************
**
*F  FunCoefficients( <hdCall> ) . . . . . .  internal function 'Coefficients'
**
**  'FunCoefficients' implements the internal function 'Coefficients'.
**
**  'Coefficients( <list>, <number> )'
**
*N  15-Jan-91 martin this function should not be here
*N  15-Jan-91 martin this function should not be called 'Coefficients'
*/
TypHandle       FunCoefficients ( hdCall )
    TypHandle           hdCall;
{
    long                pos, num, val;
    TypHandle           hdRes, hdList, hdInt;


    if ( SIZE( hdCall ) != 3 * SIZE_HD )
        return Error("usage: Coefficients( <list>, <int> )",0L,0L);

    hdList = EVAL( PTR(hdCall)[1] );
    hdInt  = EVAL( PTR(hdCall)[2] );
    if ( ! IS_LIST(hdList) || TYPE(hdInt) != T_INT)
        return Error("usage: Coefficients( <list>, <int> )",0L,0L);

    pos   = LEN_LIST( hdList );
    hdRes = NewBag( T_LIST, SIZE_PLEN_PLIST( pos ) );
    SET_LEN_PLIST( hdRes, pos );

    num = HD_TO_INT( hdInt );
    if ( num < 0 )
        return Error("Coefficients: <int> must be non negative",0L,0L);

    while ( 0 < num && 0 < pos ) {
        hdInt = ELMF_LIST( hdList, pos );
        if ( hdInt == 0 || TYPE( hdInt ) != T_INT )
          return Error("Coefficients: <list>[%d] must be a positive integer",
                       (long)pos,0L);
        val = HD_TO_INT(hdInt);
        if ( val <= 0 )
          return Error("Coefficients: <list>[%d] must be a positive integer",
                        (long)pos,0L);
        SET_ELM_PLIST( hdRes, pos, INT_TO_HD( num % val ) );
        pos--;
        num /= val;
    }

    while ( 0 < pos ) {
        SET_ELM_PLIST( hdRes, pos, INT_TO_HD( 0 ) );
        pos--;
    }

    return hdRes;
}


/****************************************************************************
**
*F  FunNUMBERHANDLES( <hdCall> )  . . . .  internal function 'NUMBER_HANDLES'
**
**  'FunNUMBERHANDLES' implements the internal function 'NUMBER_HANDLES'.
**
**  'NUMBER_HANDLES( <type> )'
*/
TypHandle       FunNUMBERHANDLES ( hdCall )
    TypHandle           hdCall;
{
    long                typ;
    TypHandle           hdTyp;


    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error("usage: NUMBER_HANDLES( <type> )",0L,0L);

    hdTyp = EVAL( PTR(hdCall)[1] );
    if (TYPE(hdTyp) != T_INT)
        return Error("usage: NUMBER_HANDLES( <type> )",0L,0L);

    typ = HD_TO_INT( hdTyp );
    if (typ < 0 || typ >= T_ILLEGAL)
        return Error("NUMBER_HANDLES: <type> must lie in [%d,%d]",
                     0L,(long)(T_ILLEGAL-1));

    return INT_TO_HD( GasmanStatTotal[typ] );
}


/****************************************************************************
**
*F  FunSIZEHANDLES( <hdCall> )  . . . . . .  internal function 'SIZE_HANDLES'
**
**  'FunSIZEHANDLES' implements the internal function 'SIZE_HANDLES'.
**
**  'SIZE_HANDLES( <type> )'
*/
TypHandle       FunSIZEHANDLES ( hdCall )
    TypHandle           hdCall;
{
    long                typ;
    TypHandle           hdTyp;


    if ( SIZE( hdCall ) != 2 * SIZE_HD )
        return Error("usage: SIZE_HANDLES( <type> )",0L,0L);

    hdTyp = EVAL( PTR(hdCall)[1] );
    if (TYPE(hdTyp) != T_INT)
        return Error("usage: SIZE_HANDLES( <type> )",0L,0L);

    typ = HD_TO_INT( hdTyp );
    if (typ < 0 || typ >= T_ILLEGAL)
        return Error("SIZE_HANDLES: <type> must lie in [%d,%d]",
                     0L,(long)(T_ILLEGAL-1));

    return INT_TO_HD( GasmanStatTSize[typ] );
}


/****************************************************************************
**
*F  InitGap( <argc>, <argv> ) . . . . . . . . . . . . . . . . initializes GAP
**
**  'InitGap' initializes GAP.
*/
void            InitGap (int argc, char *argv[])
{
    TypHandle           hd;
    long                i;
    char *              version;

    /* initialize all subpackages of GAP                                   */
    InitSystem( argc, argv );
    InitScanner();
    InitGasman();
    InitIdents();
    InitEval();

    /* create the variables last, last2, last3                             */
    HdLast  = FindIdent( "last"  );
    HdLast2 = FindIdent( "last2" );
    HdLast3 = FindIdent( "last3" );
    HdTime  = FindIdent( "time"  );

    hd = FindIdent( "VERSRC" );
    version = " 6 Feb 2018"; 
    PTR(hd)[0] = NewBag( T_STRING, SyStrlen(version)+1 );
    SyStrncat( (char*)PTR(PTR(hd)[0]), version, SyStrlen(version)+1 );
    hd = FindIdent( "VERSYS" ); version = VERSYS;
    PTR(hd)[0] = NewBag( T_STRING, SyStrlen(version)+1 );
    SyStrncat( (char*)PTR(PTR(hd)[0]), version, SyStrlen(version)+1 );

    hd = FindIdent( "LIBNAME" );
    PTR(hd)[0] = NewBag( T_STRING, (unsigned long)(SyStrlen(SyLibname)+1) );
    SyStrncat( (char*)PTR(PTR(hd)[0]), SyLibname, SyStrlen(SyLibname) );

    hd = FindIdent( "QUIET" );
    if ( SyQuiet )  PTR(hd)[0] = HdTrue;
    else            PTR(hd)[0] = HdFalse;

    hd = FindIdent( "BANNER" );
    if ( SyBanner )  PTR(hd)[0] = HdTrue;
    else             PTR(hd)[0] = HdFalse;

    /* install all internal function from this package                     */
    InstIntFunc( "Ignore",     FunIgnore     );
    InstIntFunc( "Error",      FunError      );
    InstIntFunc( "Backtrace",  FunBacktrace  );
    InstIntFunc( "WindowCmd",  FunWindowCmd  );

    InstIntFunc( "READ",       FunREAD       );
    InstIntFunc( "AUTO",       FunAUTO       );
    InstIntFunc( "Print",      FunPrint      );
    InstIntFunc( "PrintTo",    FunPrntTo     );
    InstIntFunc( "AppendTo",   FunAppendTo   );
    InstIntFunc( "LogTo",      FunLogTo      );
    InstIntFunc( "LogInputTo", FunLogInputTo );
    InstIntFunc( "ReadTest",   FunReadTest   );

    InstIntFunc( "Help",        FunHelp        );
    InstIntFunc( "Exec",        FunExec        );
    InstIntFunc( "Runtime",     FunRuntime     );
    InstIntFunc( "SizeScreen",  FunSizeScreen  );
    InstIntFunc( "TmpName",     FunTmpName     );
    InstIntFunc( "IsIdentical", FunIsIdentical );
    InstIntFunc( "HANDLE",      FunHANDLE      );
    InstIntFunc( "OBJ",         FunOBJ         );
    InstIntFunc( "TYPE",        FunTYPE        );
    InstIntFunc( "SIZE",        FunSIZE        );
    InstIntFunc( "GASMAN",      FunGASMAN      );

    InstIntFunc( "NUMBER_HANDLES",   FunNUMBERHANDLES );
    InstIntFunc( "SIZE_HANDLES",     FunSIZEHANDLES   );

    /*N  15-Jan-91 martin this function should not be here                 */
    InstIntFunc( "CoefficientsInt", FunCoefficients );

    /* read all init files, stop doing so after quiting from Error         */
    if ( ! setjmp( ErrRet ) ) {
        for ( i=0; i<sizeof(SyInitfiles)/sizeof(SyInitfiles[0]); ++i ) {
            if ( SyInitfiles[i][0] != '\0' ) {
                if ( OpenInput( SyInitfiles[i] ) ) {
                    while ( Symbol != S_EOF ) {
                        EnterKernel();
                        hd = ReadIt();
                        if ( hd != 0 )  hd = EVAL( hd );
                        if ( hd == HdReturn && PTR(hd)[0] != HdReturn )
                            Error("Read: 'return' must not be used",0L,0L);
                        else if ( hd == HdReturn )
                             Error("Read: 'quit' must not be used",0L,0L);
                        ExitKernel( (TypHandle)0 );
                   }
                   CloseInput(); /* error code ignored */
                }
                else {
                    Error("can't read from \"%s\"",(long)SyInitfiles[i],0L);
                }
            }
        }
    }

}



