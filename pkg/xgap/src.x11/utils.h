/****************************************************************************
**
*A  utils.h                     XGAP Source                      Frank Celler
**
*H  @(#)$Id: utils.h,v 1.1.1.1 1996/12/11 12:40:10 werner Exp $
**
*Y  Copyright 1995-1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This  file contains the  utility  functions  and  macros  used in   XGAP,
**  basically  the list functions  ('ELM', 'LEN',  'AddList', and 'List') and
**  the  debug  macro ('DEBUG').  This  file also  includes all the necessary
**  system and X11 include files and defines the following data types:
**
**      Boolean         "True" or "False" (defined by X11)
**      Char            a "Char" will be able to hold one character
**      Int             a 32-bit signed integer
**      Long            an integer able to hold a pointer
**      Pointer         a generic pointer
**      Short           a 16-bit signed integer
**      String          an array of chars
**      UChar           unsigned version of "Char"
**      UInt            a 32-bit unsigned integer
**      ULong           unsigned version of "Long"
**      UShort          a 16-bit unsigned integer
**
**  List( <len> )
**  -------------
**  'List' creates a new list able to hold <len> pointers of type 'Pointer'.
**
**  AddList( <lst>, <elm> )
**  -----------------------
**  'AddList' appends  the  new  element  <elm> of type  'Pointer'  to <lst>,
**  enlarging the list if necessary.
**
**  ELM( <lst>, <i> )
**  -----------------
**  'ELM' returns the <i>.th element of <lst>.
**
**  LEN( <lst> )
**  ------------
**  'LEN' returns the length of <lst>.
**
**  DEBUG( <type>, ( <debug-text>, ... ) )
**  --------------------------------------
**  'DEBUG' uses 'printf' to print the  <debug-text> in case that  'Debug' &&
**  <type> is true.  The text  is preceded by the line number  and the source
**  file name.  The following types are available:  D_LIST, D_XCMD, D_COMM.
**
*H  $Log: utils.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:40:10  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 1.7  1995/08/08  16:32:27  fceller
*H  IRIX doesn't like declaration of ioctl
*H
*H  Revision 1.6  1995/08/08  15:51:43  fceller
*H  fixed a type (HAVE_GETPTY instead of HAVE__GETPTY)
*H
*H  Revision 1.5  1995/07/28  10:37:24  fceller
*H  don't include "sys/ioctl.h" if including "termio.h"
*H
*H  Revision 1.4  1995/07/24  09:28:30  fceller
*H  reworked most parts to use nice typedefs
*H
*H  Revision 1.3  1994/06/06  08:57:24  fceller
*H  added database
*H
*H  Revision 1.2  1993/12/23  08:47:41  fceller
*H  removed malloc debug functions
*H
*H  Revision 1.1  1993/10/18  11:04:47  fceller
*H  added fast updated,  fixed timing problem
*H
*H  Revision 1.0  1993/04/05  11:42:18  fceller
*H  Initial revision
*/
#ifndef _utils_h
#define _utils_h


/****************************************************************************
**
*F  P(( <arg> ))  . . . . . . . . . . . . . . . . . . . . . . . .  prototypes
*/
#ifdef  __STDC__
# define P(ARGS) ARGS
#else
# define P(ARGS) ()
#endif


/****************************************************************************
**
*F  Autoconfig  . . . . . . . . . . . . check defines from "configure" script
*/
#if STDC_HEADERS
# define SYS_HAS_STDLIB
# define SYS_HAS_STDARG
#endif

#if HAVE_LIBC_H
# define SYS_HAS_LIBC
#else
# if HAVE_UNISTD_H
#  define SYS_HAS_UNISTD
# endif
#endif

#if HAVE_TERMIO_H
# define SYS_HAS_TERMIO
# define SYS_HAS_TCSETAW
#endif

#if HAVE_SIGNAL_H
# define SYS_HAS_SIGNAL
#endif

#if HAVE_SYS_SIGNAL_H
# define SYS_HAS_SYS_SIGNAL
#endif

#if TIME_WITH_SYS_TIME
# include <sys/time.h>
# include <time.h>
#else
# if defined(HAVE_SYS_TIME_H) || defined(SYS_HAS_SYS_TIME)
#  include <sys/time.h>
# else
#  include <time.h>
# endif
#endif

#if HAVE_GETPSEUDOTTY
# define SYS_HAS_GETPSEUDOTTY
#endif

#if HAVE__GETPTY
# define SYS_HAS_GETPTY
#endif

#if !defined(HAVE_SYS_WAIT_H) && defined(AUTOCONF)
# define SYS_HAS_UNION_WAIT
#endif


/****************************************************************************
**
*F  Include . . . . . . . . . . . . . . . . . . . . . .  system include files
*/
#include    <stdio.h>                   /* standard C i/o library          */
#include    <fcntl.h>
#include    <pwd.h>
#include    <sys/errno.h>
#include    <sys/stat.h>
#include    <sys/types.h>
#include    <sys/resource.h>
#include    <sys/wait.h>
#include    <sys/param.h>


#ifdef SYS_HAS_STDLIB
# include   <stdlib.h>                  /* standard C library              */
#endif

#ifdef SYS_HAS_UNISTD
# include   <unistd.h>                  /* another standard C library      */
#endif

#ifdef SYS_HAS_LIBC
# include   <libc.h>                    /* standard NeXT C library         */
#endif

#ifdef SYS_HAS_STDARG
# include   <stdarg.h>                  /* variable argument list          */
#endif

#ifdef SYS_HAS_TERMIO
# include   <termio.h>
#else
# include   <sgtty.h>
# include   <sys/ioctl.h>
#endif

#ifdef SYS_HAS_SIGNAL
# include   <signal.h>
#endif

#ifdef SYS_HAS_SYS_SIGNAL
# include   <sys/signal.h>
#endif


/****************************************************************************
**
*F  Include . . . . . . . . . . . . . . . . . . . . . . . . X11 include files
*/
#include    <X11/X.h>                   /* X11 basic definition            */
#include    <X11/Xos.h>
#include    <X11/Xatom.h>
#include    <X11/Xlib.h>
#include    <X11/StringDefs.h>
#include    <X11/keysym.h>

#include    <X11/Intrinsic.h>           /* X Intrinsic                     */
#include    <X11/IntrinsicP.h>
#include    <X11/CoreP.h>
#include    <X11/Composite.h>
#include    <X11/Shell.h>

#include    <X11/cursorfont.h>          /* cursor font                     */

#include    <X11/Xaw/AsciiText.h>       /* Athena widgets                  */
#include    <X11/Xaw/Box.h>
#include    <X11/Xaw/Cardinals.h>
#include    <X11/Xaw/Command.h>
#include    <X11/Xaw/Dialog.h>
#include    <X11/Xaw/Form.h>
#include    <X11/Xaw/Label.h>
#include    <X11/Xaw/List.h>
#include    <X11/Xaw/MenuButton.h>
#include    <X11/Xaw/Paned.h>
#include    <X11/Xaw/Scrollbar.h>
#include    <X11/Xaw/SimpleMenu.h>
#include    <X11/Xaw/SmeBSB.h>
#include    <X11/Xaw/SmeLine.h>
#include    <X11/Xaw/Text.h>
#include    <X11/Xaw/TextP.h>
#include    <X11/Xaw/TextSink.h>
#include    <X11/Xaw/TextSrc.h>
#include    <X11/Xaw/TextSrcP.h>
#include    <X11/Xaw/Viewport.h>
#include    <X11/Xaw/ViewportP.h>
#include    <X11/Xaw/XawInit.h>


/****************************************************************************
**
*F  Prototypes  . . . . . . . . . . . . . . . . . . . . . . system prototypes
*/
#if !defined(SYS_HAS_UNISTD) && !defined(SYS_HAS_LIBC)
extern int write();
#endif

#if !defined(SYS_HAS_TIME_PROTO) && !defined(AUTOCONF)
extern int time();
#endif

#if defined(SYS_HAS_PID_T) || defined(AUTOCONF)
extern pid_t wait3();
#else
extern int wait3();
#endif

extern int select();

/* IRIX System V.4 running IRIX Release 5.3 already defines ioctl and  */
/* therefore doesn't like the declaration of ioctl                     */

/* extern int ioctl(); */


/****************************************************************************
**

*T  Boolean . . . . . . . . . . . . . . . . . . . . . . . . . a boolean value
*/
#if 0
typedef enum { False, True } Boolean;
#endif


/****************************************************************************
**
*T  Char  . . . . . . . . . . . . . . . . . . . . . . . . . . . . a character
*/
typedef char Char;


/****************************************************************************
**
*T  Int . . . . . . . . . . . . . . . . . . . . . . . . . . . a signed 32-bit
*/
typedef int Int;


/****************************************************************************
**
*T  Long  . . . . . . . . . . . . . . a signed integer able to hold a pointer
*/
typedef long Long;


/****************************************************************************
**
*T  Pointer . . . . . . . . . . . . . . . . . . . . . . . . a generic pointer
*/
typedef void * Pointer;


/****************************************************************************
**
*T  Short . . . . . . . . . . . . . . . . . . . . . . . . . . a signed 16-bit
*/
typedef short Short;


/****************************************************************************
**
*T  String  . . . . . . . . . . . . . . . . . . . . . . . . an array of chars
*/
#if 0
typedef char * String;
#endif


/****************************************************************************
**
*T  UChar . . . . . . . . . . . . . . . . . . . . . . . an unsigned character
*/
typedef unsigned char UChar;


/****************************************************************************
**
*T  UInt  . . . . . . . . . . . . . . . . . . . . . . . .  an unsigned 32-bit
*/
typedef unsigned int UInt;


/****************************************************************************
**
*T  ULong . . . . . . . . . . . .  an unsigned integer able to hold a pointer
*/
typedef unsigned long ULong;


/****************************************************************************
**
*T  UShort  . . . . . . . . . . . . . . . . . . . . . . .  an unsigned 16-bit
*/
typedef unsigned short UShort;


/****************************************************************************
**

*F  DEBUG(( <str> ))  . . . . . . . . . . . . . . . print <str> as debug info
*/
extern Int Debug;

#define D_LIST		1
#define D_XCMD          2
#define D_COMM          4

#define DEBUG(a,b) {                                       \
            if ( Debug & a ) {                             \
                printf( "%04d:%s: ", __LINE__, __FILE__ ); \
                printf b;                                  \
            }                                              \
        } while(0)


/****************************************************************************
**
*F  MAX( <a>, <b> ) . . . . . . . . . . . . . . . . .  maximum of <a> and <b>
*/
#undef  MAX
#define MAX(a,b)        (((a) < (b)) ? (b) : (a))


/****************************************************************************
**
*F  MIN( <a>, <b> ) . . . . . . . . . . . . . . . . .  minimum of <a> and <b>
*/
#undef  MIN
#define MIN(a,b)        (((a) < (b)) ? (a) : (b))


/****************************************************************************
**

*T  TypeList  . . . . . . . . . . . . . . . . . . . . . . . .  list structure
*/
typedef struct _list
{
    UInt        size;
    UInt        len;
    Pointer   * ptr;
}
* TypeList;


/****************************************************************************
**
*F  ELM( <lst>, <i> ) . . . . . . . . . . . . . . . . <i>th element of a list
*/
#define ELM(lst,i)      (lst->ptr[i])


/****************************************************************************
**
*F  LEN( <lst> )  . . . . . . . . . . . . . . . . . . . . .  length of a list
*/
#define LEN(lst)        (lst->len)


/****************************************************************************
**
*F  AddList( <lst>, <elm> ) . . . . . . . .  add list element <elm> to <list>
*/
#ifdef DEBUG_ON
    extern void 	ADD_LIST P(( String, Int, TypeList, Pointer ));
#   define AddList(a,b)	ADD_LIST( __FILE__, __LINE__, a, b )
#else
    extern void 	AddList P(( TypeList, Pointer ));
#endif


/****************************************************************************
**
*F  List( <len> )   . . . . . . . . . . . . . . . . . . .   create a new list
*/
#ifdef DEBUG_ON
    extern TypeList 	LIST P(( String, Int, UInt ));
#   define List(a) 	LIST( __FILE__, __LINE__, a )
#else
    extern TypeList 	List P(( UInt ));
#endif

#endif
