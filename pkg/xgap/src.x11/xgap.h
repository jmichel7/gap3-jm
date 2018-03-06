/****************************************************************************
**
*A  xgap.h                      XGAP Source                      Frank Celler
**
*H  @(#)$Id: xgap.h,v 1.1.1.1 1996/12/11 12:40:11 werner Exp $
**
*Y  Copyright 1995-1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
*H  $Log: xgap.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:40:11  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 1.3  1995/08/08  16:32:52  fceller
*H  added missing extern before CursorTL
*H
*H  Revision 1.2  1995/07/24  09:28:30  fceller
*H  reworked most parts to use nice typedefs
*H
*H  Revision 1.1  1993/10/18  11:04:47  fceller
*H  added fast updated,  fixed timing problem
*H
*H  Revision 1.0  1993/04/05  11:42:18  fceller
*H  Initial revision
*H
*/
#ifndef _xgap_h
#define _xgap_h

/****************************************************************************
**
*T  TypeMenuItem  . . . . . . . . . . . . . . . . . . . . .  menu description
*/
#define S_ALWAYS          0
#define S_INPUT_ONLY	  1
#define S_ERROR_ONLY      2
#define S_NORMAL_ONLY     3
#define S_RUNNING_ONLY    4
#define S_HELP_ONLY       5

typedef struct _menu_item
{
  char 	  * label;
  void      (*click)();
  int       sensitive;
  Widget    entry;
}
TypeMenuItem;


/****************************************************************************
**

*V  AppContext	. . . . . . . . . . . . . . . . . . . . .  aplication context
*/
extern XtAppContext AppContext;


/****************************************************************************
**
*V  GapDisplay  . . . . . . . . . . . . . . . . . . . . . . . current display
*/
extern Display * GapDisplay;


/****************************************************************************
**
*V  GapScreen . . . . . . . . . . . . . . . . . . . . . . . .  current screen
*/
extern long GapScreen;


/****************************************************************************
**
*V  GapTalk . . . . . . . . . . . . . . . . . . . . . . . . . gap text window
*/
extern Widget GapTalk;


/****************************************************************************
**
*V  GapState  . . . . . . . . . . . . . . . . . . . . . . . . . status of gap
*/
#define GAP_NOGAP       0
#define	GAP_RUNNING	1
#define GAP_INPUT       2
#define GAP_ERROR       3
#define GAP_HELP        4

extern Int GapState;


/****************************************************************************
**
*V  MyRootWindow  . . . . . . . . . . . . . . . . . . . . current root window
*/
extern Drawable MyRootWindow;


/****************************************************************************
**
*V  SpyMode . . . . . . . . . . . . . . . . . . . . copy GAP output to stderr
*/
extern Boolean SpyMode;


/****************************************************************************
**
*V  WmDeleteWindowAtom  . . . . . . .  window manager "delete window" request
*/
extern Atom WmDeleteWindowAtom;


/****************************************************************************
**
*V  XGap  . . . . . . . . . . . . . . . . . . . . . . . . . .  toplevel shell
*/
extern Widget XGap;


/****************************************************************************
**

*V  CheckMarkSymbol . . . . . . . . . . . . . symbol for checked menu entries
*/
extern Pixmap CheckMarkSymbol;


/****************************************************************************
**
*V  CursorTL  . . . . . . . . . . . . . . . . . . . . . . . .  top left arrow
*/
extern Cursor CursorTL;


/****************************************************************************
**
*V  EmptyMarkSymbol . . . . . . . . . . . . symbol for unchecked menu entries
*/
extern Pixmap EmptyMarkSymbol;


/****************************************************************************
**
*V  MenuSymbol	. . . . . . . . . . . . . . . . .  symbol for drop down menus
*/
extern Pixmap MenuSymbol;


/****************************************************************************
**
*V  ExMarkSymbol  . . . . . . . . . . . . . . . . . . . . .  exclamation mark
*/
extern Pixmap ExMarkSymbol;


/****************************************************************************
**

*P  Prototypes  . . . . . . . . . . . prototypes of public gap text functions
*/
extern void SimulateInput P(( String ));
extern void UpdateMenus P(( Int ));
extern void UpdateMemoryInfo P(( Int, Int ));


#endif
