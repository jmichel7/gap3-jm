/****************************************************************************
**
*A  popdial.h			XGAP source	                 Frank Celler
**
*H  @(#)$Id: popdial.h,v 1.1.1.1 1996/12/11 12:40:08 werner Exp $
**
*Y  Copyright 1995-1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
**  This file contains functions for popping up dialogs.
**
*H  $Log: popdial.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:40:08  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 1.4  1995/07/24  09:28:30  fceller
*H  reworked most parts to use nice typedefs
*H
*H  Revision 1.3  1995/01/24  11:15:18  fceller
*H  a dot cursor is show if the cursor is inside a button
*H
*H  Revision 1.2  1994/06/06  08:57:24  fceller
*H  added database
*H
*H  Revision 1.1  1993/10/18  11:04:47  fceller
*H  added fast updated,  fixed timing problem
*H
*H  Revision 1.0  1993/04/05  11:42:18  fceller
*H  Initial revision
*/
#ifndef _popdial_h
#define _popdial_h

#include    "utils.h"			/* utility functions */


/****************************************************************************
**
*D  PD_YES  . . . . . . . . . . . . . . . . . . . . . . . . exit button "yes"
*D  PD_NO   . . . . . . . . . . . . . . . . . . . . . . . .  exit button "no"
*D  PD_OK   . . . . . . . . . . . . . . . . . . . . . . . .  exit button "OK"
*D  PD_CANCEL   . . . . . . . . . . . . . . . . . . . .  exit button "cancel"
*D  PD_ABORT  . . . . . . . . . . . . . . . . . . . . . . exit button "abort"
*D  PD_RETRY  . . . . . . . . . . . . . . . . . . . . . . exit button "retry"
*D  PD_APPEND   . . . . . . . . . . . . . . . . . . . .  exit button "append"
*D  PD_OVERWRITE  . . . . . . . . . . . . . . . . . . exit button "overwrite"
*/
#define	PD_YES	        0x0001
#define PD_NO           0x0002
#define PD_OK           0x0004
#define PD_CANCEL       0x0008
#define PD_ABORT        0x0010
#define PD_RETRY        0x0020
#define PD_APPEND  	0x0040
#define PD_OVERWRITE	0x0080


/****************************************************************************
**
*T  TypePopupDialog . . . . . . . . . . . . . . . pointer to dialog structure
*/
typedef struct _popup_dialog
{
    Widget	    topLevel;
    Widget          popupShell;
    Widget          dialog;
    XtAppContext    context;
    Int             result;
    Int             button;
    Int             defaultButton;
    Widget          buttons[PD_OVERWRITE+1];
}
* TypePopupDialog;


/****************************************************************************
**
*P  Prototypes  . . . . . . . . . . . prototypes of public gap text functions
*/
extern TypePopupDialog CreatePopupDialog P(( XtAppContext, Widget, String,
					     Int, Int, Pixmap ));
extern Int  PopupDialog P(( TypePopupDialog, String, String, String* ));
extern void PopupDialogBrokenWM P(( void ));

#endif
