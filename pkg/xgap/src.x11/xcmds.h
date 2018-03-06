/****************************************************************************
**
*A  xcmds.h                     XGAP Source                      Frank Celler
**
*H  @(#)$Id: xcmds.h,v 1.1.1.1 1996/12/11 12:40:11 werner Exp $
**
*Y  Copyright 1995-1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
*H  $Log: xcmds.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:40:11  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 1.5  1995/07/24  09:28:30  fceller
*H  reworked most parts to use nice typedefs
*H
*H  Revision 1.4  1993/12/23  08:49:29  fceller
*H  changed some comments
*H
*H  Revision 1.3  1993/10/18  11:04:47  fceller
*H  added fast updated,  fixed timing problem
*H
*H  Revision 1.2  1993/04/13  07:30:44  fceller
*H  added text selectors
*H
*H  Revision 1.1  1993/04/13  07:16:39  fceller
*H  added missing __STDC__ define
*H
*H  Revision 1.0  1993/04/05  11:42:18  fceller
*H  Initial revision
*/
#ifndef _xcmds_h
#define _xcmds_h

#include "utils.h"


/****************************************************************************
**
*T  TypeGapWindow . . . . . . . . . . . . . . . . description of a gap window
*/
typedef struct _gap_window
{
    Widget              top;
    Widget              box;
    Widget              viewport;
    Widget              draw;
    Widget              text;
    Boolean             used;
    TypeList            menus;
    Int                 height;
    Int                 width;
    Int                 line_width;
    Int                 color;
    Boolean             fast_update;
}
TypeGapWindow;


/****************************************************************************
**
*T  TypeTextSelector  . . . . . . . . . . . .  description of a text selector
*/
typedef struct _text_selector
{
    Widget              top;
    Widget              list;
    String            * text;
    TypeList            buttons;
}
TypeTextSelector;


/****************************************************************************
**
*T  TypeArg . . . . . . . . . . . . . . . . . . . . . . . . . . . . arguments
*/
#define MAX_ARG		10

typedef struct _arg
{
    TypeGapWindow     * win;
    TypeTextSelector  * sel;
    XFontStruct       * font;
    Int                 iargs[MAX_ARG];
    String              sargs[MAX_ARG];
    String              opts;
}
TypeArg;


/****************************************************************************
**
*T  TypeWindowCommand . . . . . . . . . . . . . . . . description of commands
*/
typedef struct _window_command
{
    String	name;
    String      args;
    Boolean     (*func) P(( TypeArg* ));
}
TypeWindowCommand;


/****************************************************************************
**
*T  TypeMenu  . . . . . . . . . . . . . . . . . . . . . . description of menu
*/
typedef struct _menu
{
    Widget      button;
    Widget      shell;
    TypeList	entries;
    String      name;
    String      string;
}
TypeMenu;


/****************************************************************************
**
*T  TypeMenuData  . . . . . . . . . . . . . . . . . . . . . . . .  menu entry
*/
typedef struct _menu_data
{
    Widget      shell;
    Int         window;
    Int         popup;
    Int         pane;
}
TypeMenuData;


/****************************************************************************
**
*T  TypePaneData  . . . . . . . . . . . . . . . . . . . . .  popup menu entry
*/
typedef struct _pane_data
{
    Widget      shell;
    Int         popup;
    Int         pane;
}
TypePaneData;


/****************************************************************************
**
*P  Prototypes	. . . . . . . . . . . . . . . . . . . . . function prototypes
*/
extern void     InitXCMDS P(( void ));
extern void     ExitXMCDS P(( void ));
extern void     UpdateXCMDS P(( Int ));
extern Boolean  GapWindowCmd P(( String, Int ));

#endif
