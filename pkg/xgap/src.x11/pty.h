/****************************************************************************
**
*A  pty.h                       XGAP source                      Frank Celler
**
*H  @(#)$Id: pty.h,v 1.1.1.1 1996/12/11 12:40:09 werner Exp $
**
*Y  Copyright 1995-1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
*H  $Log: pty.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:40:09  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 1.3  1995/07/24  09:28:30  fceller
*H  reworked most parts to use nice typedefs
*H
*H  Revision 1.2  1994/06/06  08:57:24  fceller
*H  added database
*H
*H  Revision 1.1  1994/06/03  10:50:09  fceller
*H  fixed exec problem (again)
*H
*H  Revision 1.0  1993/04/05  11:42:18  fceller
*H  Initial revision
*/
#ifndef _pty_h
#define _pty_h


/****************************************************************************
**
*V  QuitGapCtrlD  . . . . . . . . . . . . . . . . . . . . . . . quit on CTR-D
*/
extern Boolean QuitGapCtrlD;


/****************************************************************************
**
*V  ScreenSizeBuffer  . . . . . . . . . . . . . .  screen size change command
*/
extern char ScreenSizeBuffer[];


/****************************************************************************
**
*V  ExecRunning . . . . . . . . . . . . . . . . . .  external program running
*/
extern Boolean ExecRunning;


/****************************************************************************
**
*P  Prototypes  . . . . . . . . . . . . . . . . . . . . . function prototypes
*/
extern Int  CheckCaretPos P(( Int, Int ));
extern int  StartGapProcess P(( String, String argv[] ));
extern void GapOutput P(( XtPointer, Int*,  XtInputId ));
extern void InterruptGap P(( void ));
extern void KeyboardInput P(( String, Int ));
extern void KillGap P(( void ));
extern void StoreInput P(( String, Int ));
extern void ProcessStoredInput P(( Int ));


/****************************************************************************
**
*D  ReadGap( <buf>, <len> ) . . . . . . . . . . . . . . . read bytes from gap
*D  WriteGap( <buf>, <len> )  . . . . . . . . . . . . . .  write bytes to gap
*/
#ifdef DEBUG_ON
    extern Int              READ_GAP  P(( String, Int, String, Int ));
    extern void             WRITE_GAP P(( String, Int, String, Int ));
#   define ReadGap(a,b)	    READ_GAP ( __FILE__, __LINE__, a, b )
#   define WriteGap(a,b)    WRITE_GAP( __FILE__, __LINE__, a, b )
#else
    extern Int              ReadGap  P(( String, Int ));
    extern void             WriteGap P(( String, Int ));
#endif

#endif
