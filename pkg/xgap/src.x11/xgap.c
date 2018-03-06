/****************************************************************************
**
*A  xgap.c                      XGAP Source                      Frank Celler
**
*H  @(#)$Id: xgap.c,v 1.1.1.1 1996/12/11 12:40:11 werner Exp $
**
*Y  Copyright 1995-1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
*H  $Log: xgap.c,v $
*H  Revision 1.1.1.1  1996/12/11 12:40:11  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 1.10  1995/08/16  11:16:18  fceller
*H  verbose signal handler if DEBUG_ON is set
*H
*H  Revision 1.9  1995/07/24  09:28:30  fceller
*H  reworked most parts to use nice typedefs
*H
*H  Revision 1.8  1994/06/11  12:41:43  fceller
*H  added 'SYS_HAS_TIME_PROTO'
*H
*H  Revision 1.7  1994/06/06  08:56:07  fceller
*H  updated database
*H
*H  Revision 1.6  1994/06/03  10:50:09  fceller
*H  fixed exec problem (again)
*H
*H  Revision 1.5  1993/12/23  08:45:16  fceller
*H  added "gap drop prompt"
*H  renamed 'SignalHandler' to 'MySignalHandler' to avoid linux conflict
*H
*H  Revision 1.4  1993/10/18  11:04:47  fceller
*H  added fast updated,  fixed timing problem
*H
*H  Revision 1.3  1993/10/06  16:17:09  fceller
*H  new gap package mode
*H
*H  Revision 1.2  1993/08/12  13:48:06  fceller
*H  added help menu
*H
*H  Revision 1.1  1993/07/29  07:02:30  fceller
*H  added 'Exec' hack
*H
*H  Revision 1.0  1993/04/05  11:42:18  fceller
*H  Initial revision
*/
#include    "utils.h"			/* utility functions */

#include    "popdial.h"			/* popup dialogs */
#include    "gapgraph.h"		/* gap graphic sheet */
#include    "gaptext.h"			/* gap text sheet */
#include    "popdial.h"                 /* popup dialogs */
#include    "xcmds.h"
#include    "pty.h"
#include    "selfile.h"

#include    "xgap.h"


/* * * * * * * * * * * * * * * global variables  * * * * * * * * * * * * * */


/****************************************************************************
**

*V  AppContext	. . . . . . . . . . . . . . . . . . . . .  aplication context
*/
XtAppContext AppContext;


/****************************************************************************
**
*V  GapDisplay  . . . . . . . . . . . . . . . . . . . . . . . current display
*/
Display * GapDisplay;


/****************************************************************************
**
*V  GapScreen . . . . . . . . . . . . . . . . . . . . . . . .  current screen
*/
long GapScreen;


/****************************************************************************
**
*V  GapState  . . . . . . . . . . . . . . . . . . . . . . . . . status of gap
*/
#define GAP_NOGAP       0
#define	GAP_RUNNING	1
#define GAP_INPUT       2
#define GAP_ERROR       3
#define GAP_HELP        4

Int GapState = GAP_NOGAP;


/****************************************************************************
**
*V  GapTalk . . . . . . . . . . . . . . . . . . . . . . . . . gap text window
*/
Widget GapTalk;


/****************************************************************************
**
*V  MyRootWindow  . . . . . . . . . . . . . . . . . . . . current root window
*/
Drawable MyRootWindow;


/****************************************************************************
**
*V  SpyMode . . . . . . . . . . . . . . . . . . . . copy GAP output to stderr
*/
Boolean SpyMode = False;


/****************************************************************************
**
*V  WmDeleteWindowAtom  . . . . . . .  window manager "delete window" request
*/
Atom WmDeleteWindowAtom;


/****************************************************************************
**
*V  XGap  . . . . . . . . . . . . . . . . . . . . . . . . . .  toplevel shell
*/
Widget XGap;


/* * * * * * * * * * * * * * * various symbols * * * * * * * * * * * * * * */



/****************************************************************************
**

*V  CheckMarkSymbol . . . . . . . . . . . . . symbol for checked menu entries
*/
Pixmap CheckMarkSymbol;


/****************************************************************************
**
*V  CursorTL  . . . . . . . . . . . . . . . . . . . . . . . .  top left arrow
*/
Cursor CursorTL;


/****************************************************************************
**
*V  EmptyMarkSymbol . . . . . . . . . . . . symbol for unchecked menu entries
*/
Pixmap EmptyMarkSymbol;


/****************************************************************************
**
*V  ExMarkSymbol  . . . . . . . . . . . . . . . . . . . . .  exclamation mark
*/
Pixmap ExMarkSymbol;


/****************************************************************************
**
*V  MenuSymbol	. . . . . . . . . . . . . . . . .  symbol for drop down menus
*/
Pixmap MenuSymbol;


/* * * * * * * * * * * * * * * local variables * * * * * * * * * * * * * * */


/****************************************************************************
**

*V  CommandOptions  . . . . . . . . . . . . . . . . . .  command line options
*/
static XrmOptionDescRec CommandOptions[] =
{
    { "-colorModel",    "*colorModel",      XrmoptionSepArg,  0 },
    { "-colors",        "*colors",          XrmoptionSepArg,  0 },
    { "-huge",          "*hugeFont",        XrmoptionSepArg,  0 },
    { "-hugeFont",      "*hugeFont",        XrmoptionSepArg,  0 },
    { "-large",         "*largeFont",       XrmoptionSepArg,  0 },
    { "-largeFont",     "*largeFont",       XrmoptionSepArg,  0 },
    { "-normal", 	"*normalFont",      XrmoptionSepArg,  0 },
    { "-normalFont", 	"*normalFont",      XrmoptionSepArg,  0 },
    { "-positionTitle", "*titlePosition",   XrmoptionSepArg,  0 },
    { "-small",         "*smallFont",       XrmoptionSepArg,  0 },
    { "-smallFont",     "*smallFont",       XrmoptionSepArg,  0 },
    { "-tiny",		"*tinyFont",	    XrmoptionSepArg,  0 },
    { "-tinyFont",      "*tinyFont",        XrmoptionSepArg,  0 },
    { "-titlePosition", "*titlePosition",   XrmoptionSepArg,  0 },
    { "-tp",            "*titlePosition",   XrmoptionSepArg,  0 },
};


/****************************************************************************
**
*V  FallbackResources . . . . . . . . . . . . . . . . . . . default resources
*/
static char *FallbackResources[] =
{
    "*menu.line.height:                       10",
    "*xgapMenu*shapeStyle:                    ShapeOval",
    "*xgapDialog*shapeStyle:                  ShapeOval",

    /* gap talk window */
    "*xgapTalk.height:                        600",
    "*xgapTalk.width:                         600",
    "*xgapMenu.showGrip:                      False",
    "*xgapTalk.showGrip:                      False",
    "*xgapTalk.quitGapCtrD:                   True",
    "*xgapTalk.pasteGapPrompt:                True",

    /* gap menu */
    "*xgapMenu.gapButton.label:               Gap",
    "*xgapMenu.gapButton*pastePrompt.label:   Paste 'gap>'",
    "*xgapMenu.gapButton*quitGapCTRD.label:   Quit on CTR-D",
    "*xgapMenu.gapButton*editFile.label:      Edit File ...",
    "*xgapMenu.gapButton*readFile.label:      Read File ...",
    "*xgapMenu.gapButton*reloadLib.label:     Reload Library",
    "*xgapMenu.gapButton*changeLib.label:     Change Library ...",
#ifdef DEBUG_ON
    "*xgapMenu.gapButton*resyncGap.label:     Resync with GAP",
#endif
    "*xgapMenu.gapButton*quit.label:          Quit Gap",
    "*xgapMenu.gapButton*kill.label:          Kill Gap",

    /* run menu */
    "*xgapMenu.runButton.label:               Run",
    "*xgapMenu.runButton*quitBreak.label:     Leave Breakloop",
    "*xgapMenu.runButton*contBreak.label:     Continue Execution",
    "*xgapMenu.runButton*interrupt.label:     Interrupt",
    "*xgapMenu.runButton*garbColl.label:      Collect Garbage",
    "*xgapMenu.runButton*garbMesg.label:      Toggle GC Messages",
    "*xgapMenu.runButton*infoRead.label:      Toggle Library Read Mesg",

    /* help menu */
    "*xgapMenu.helpButton.label:              Help",
    "*xgapMenu.helpButton*copyHelp.label:     Copyright",
    "*xgapMenu.helpButton*helpHelp.label:     Helpsystem",
    "*xgapMenu.helpButton*abotHelp.label:     About GAP",
    "*xgapMenu.helpButton*chpsHelp.label:     Chapters",
    "*xgapMenu.helpButton*secsHelp.label:     Sections",
    "*xgapMenu.helpButton*nchpHelp.label:     Next Chapter",
    "*xgapMenu.helpButton*pchpHelp.label:     Previous Chapter",
    "*xgapMenu.helpButton*nextHelp.label:     Next Help Section",
    "*xgapMenu.helpButton*prevHelp.label:     Previous Help Section",
    "*xgapMenu.helpButton*contHelp.label:     Continue Help",
    "*xgapMenu.helpButton*quitHelp.label:     Quit Help",

    /* gap graphic window */
    "*xgapWindowViewport.width:               400",
    "*xgapWindowViewport.height:              600",

    /* query a input file name */
#ifdef NO_FILE_SELECTOR
    "*queryFileName.xgapDialog.icon:          Term",
#else
    "*selFileCancel*ShapeStyle:               Oval",
    "*selFileCancel*label:                    Cancel",
    "*selFileOK*ShapeStyle:                   Oval",
    "*selFileOK*label:                        OK",
    "*selFileHome*ShapeStyle:                 Oval",
    "*selFileHome*label:                      Home",
#endif
    0
};



/* * * * * * * * * * * * * gap talk window menus * * * * * * * * * * * * * */


/****************************************************************************
**

*V  GapMenu . . . . . . . . . . . . . . . . . . . . . . . . xgap's "Gap" menu
**
*/
static void MenuQuitGap ()   { KeyboardInput( "@C@A@Kquit;\nquit;\n", 18 ); }
static void MenuKillGap ()   { KillGap();                                   }
static void MenuReloadLib () { SimulateInput( "ReadLib( \"init\" );\n" );   }

#ifdef DEBUG_ON
static void MenuResyncGap ()
{
    ExecRunning = 0;
    GapState = GAP_INPUT;
    UpdateMenus(GapState);
    UpdateXCMDS(True);
    ProcessStoredInput(0);
}
#endif

static void MenuPastePrompt ( item )
    TypeMenuItem *      item;
{
    static Boolean	paste = False;

    paste = !paste;
    GTDropGapPrompt( GapTalk, !paste );
    if ( paste )
	XtVaSetValues( item->entry, XtNrightBitmap, CheckMarkSymbol, 0 );
    else
	XtVaSetValues( item->entry, XtNrightBitmap, EmptyMarkSymbol, 0 );
}

static void MenuQuitGapCTRD ( item )
    TypeMenuItem *      item;
{
    QuitGapCtrlD = !QuitGapCtrlD;
    if ( QuitGapCtrlD )
	XtVaSetValues( item->entry, XtNrightBitmap, CheckMarkSymbol, 0 );
    else
	XtVaSetValues( item->entry, XtNrightBitmap, EmptyMarkSymbol, 0 );
}

#ifndef NO_FILE_SELECTOR
void MenuReadFile ( item )
    TypeMenuItem *	item;
{
    Int			res;
    String              str;
    String              input;

    res = XsraSelFile( XGap, "Select a File", 0, 0, &str );
    if ( res && str[0] )
    {
	input = XtMalloc( strlen(str)+30 );
	strcpy( input, "Read( \"" );
	strcat( input, str );
	strcat( input, "\" );\n" );
	SimulateInput(input);
	XtFree(input);
	XtFree(str);
    }
}
#endif

static TypeMenuItem GapMenu[] =
{
    { "pastePrompt",    MenuPastePrompt,        S_ALWAYS,	0 },
    { "quitGapCTRD",    MenuQuitGapCTRD,        S_ALWAYS,       0 },
    { "-----------",    0,                      0,              0 },
    { "reloadLib",  	MenuReloadLib, 		S_NORMAL_ONLY,	0 },
#ifndef NO_FILE_SELECTOR
    { "readFile",       MenuReadFile,		S_NORMAL_ONLY,  0 },
#endif
    { "-----------",    0,                      0,              0 },
#ifdef DEBUG_ON
    { "resyncGap",      MenuResyncGap,          S_ALWAYS,       0 },
#endif
    { "quit",           MenuQuitGap,            S_ALWAYS,       0 },
    { "kill",           MenuKillGap,            S_ALWAYS,       0 },
    { 0,                0,             	        0,              0 }
};


/****************************************************************************
**
*V  HelpMenu  . . . . . . . . . . . . . . . . . . . . . .  xgap's "Help" menu
**
*/
static void MenuChapters ()     { SimulateInput( "?Chapters\n" ); }
static void MenuSections ()     { SimulateInput( "?Sections\n" ); }
static void MenuCopyright ()    { SimulateInput( "?Copyright\n" );}
static void MenuAboutGAP ()     { SimulateInput( "?About GAP\n" );}
static void MenuHelp ()         { SimulateInput( "?Help\n" );     }
static void MenuNextHelp ()     { SimulateInput( "?>\n" );        }
static void MenuNextChapter ()  { SimulateInput( "?>>\n" );       }
static void MenuPrevChapter ()  { SimulateInput( "?<<\n" );       }
static void MenuPrevHelp ()     { SimulateInput( "?<\n" );        }
static void MenuContHelp ()     { KeyboardInput( " ", 1 );        }
static void MenuQuitHelp ()     { KeyboardInput( "q", 1 );        }


static TypeMenuItem HelpMenu[] =
{
    { "copyHelp",       MenuCopyright,          S_INPUT_ONLY,   0 },
    { "helpHelp",       MenuHelp,               S_INPUT_ONLY,   0 },
    { "abotHelp",       MenuAboutGAP,           S_INPUT_ONLY,   0 },
    { "---------",      0,                      0,              0 },
    { "chpsHelp",       MenuChapters,           S_INPUT_ONLY,   0 },
    { "secsHelp",       MenuSections,           S_INPUT_ONLY,   0 },
    { "---------",      0,                      0,              0 },
    { "nchpHelp",       MenuNextChapter,        S_INPUT_ONLY,   0 },
    { "pchpHelp",       MenuPrevChapter,        S_INPUT_ONLY,   0 },
    { "nextHelp",       MenuNextHelp,           S_INPUT_ONLY,   0 },
    { "prevHelp",       MenuPrevHelp,           S_INPUT_ONLY,   0 },
    { "contHelp",       MenuContHelp,           S_HELP_ONLY,    0 },
    { "quitHelp",       MenuQuitHelp,           S_HELP_ONLY,    0 },
    { 0,              	0,               	0,              0 }
};


/****************************************************************************
**
*V  RunMenu . . . . . . . . . . . . . . . . . . . . . . . . xgap's "Run" menu
**
*/
static void MenuInterrupt () { InterruptGap();                            }
static void MenuQuitBreak () { SimulateInput( "quit;\n" );                }
static void MenuContBreak () { SimulateInput( "return;\n" );              }
static void MenuGarbColl ()  { SimulateInput( "GASMAN(\"collect\");\n" ); }
static void MenuGarbMesg ()  { SimulateInput( "GASMAN(\"message\");\n" ); }
static void MenuInfoRead ()  { SimulateInput(
"if InfoRead1=Print then InfoRead1:=Ignore; else InfoRead1:=Print; fi;\n"); }

static TypeMenuItem RunMenu[] =
{
    { "interrupt",      MenuInterrupt,          S_RUNNING_ONLY, 0 },
    { "---------",      0,                      0,              0 },
    { "quitBreak",    	MenuQuitBreak,    	S_ERROR_ONLY,	0 },
    { "contBreak",      MenuContBreak,          S_ERROR_ONLY,   0 },
    { "---------",      0,                      0,              0 },
    { "garbColl",       MenuGarbColl,           S_INPUT_ONLY,   0 }, 
    { "garbMesg",       MenuGarbMesg,           S_INPUT_ONLY,   0 },
    { "infoRead",       MenuInfoRead,           S_INPUT_ONLY,   0 },
    { 0,                0,             	        0,              0 }
};


/****************************************************************************
**

*F  CreateMenu( <button>, <items> ) . . . . . . . . . . . . create a pop menu
**
**  RESOURCES
**    *menu.line.height
**        height of menu line separator, default 10
*/
static TypeList ListInputOnly   = 0;
static TypeList ListErrorOnly   = 0;
static TypeList ListNormalOnly  = 0;
static TypeList ListRunningOnly = 0;
static TypeList ListHelpOnly    = 0;

static void CreateMenu P(( Widget, TypeMenuItem* ));
static void MenuSelected P(( Widget, TypeMenuItem*, caddr_t ));

static void CreateMenu ( button, items )
    Widget	    button;
    TypeMenuItem  * items;
{
    Widget          menu;

    /* if this is the first call,  create lists */
    if ( ListInputOnly == 0 )
    {
	ListErrorOnly   = List(0);
	ListHelpOnly    = List(0);
	ListInputOnly   = List(0);
	ListNormalOnly  = List(0);
	ListRunningOnly = List(0);
    }

    /* create new simple menu */
    menu = XtCreatePopupShell( "menu", simpleMenuWidgetClass, button, 0, 0 );

    /* and add menu buttons */
    for ( ;  items->label != 0;  items++ )
    {
	if ( *(items->label) == '-' )
	    (void) XtVaCreateManagedWidget( "line",
		       smeLineObjectClass, menu, NULL );
	else
	{
	    items->entry = XtVaCreateManagedWidget(
			       items->label, smeBSBObjectClass, menu,
			       XtNrightMargin, 14,
                               XtNrightBitmap, EmptyMarkSymbol,
			       0 );
	    XtAddCallback( items->entry, XtNcallback,
			   (XtCallbackProc)MenuSelected, items );
	    switch ( items->sensitive )
	    {
		case S_INPUT_ONLY:
		    AddList( ListInputOnly, items->entry );
		    XtVaSetValues( items->entry, XtNsensitive, False, 0 );
		    break;
		case S_ERROR_ONLY:
		    AddList( ListErrorOnly, items->entry );
		    XtVaSetValues( items->entry, XtNsensitive, False, 0 );
		    break;
		case S_NORMAL_ONLY:
		    AddList( ListNormalOnly, items->entry );
		    XtVaSetValues( items->entry, XtNsensitive, False, 0 );
		    break;
		case S_RUNNING_ONLY:
		    AddList( ListRunningOnly, items->entry );
		    XtVaSetValues( items->entry, XtNsensitive, False, 0 );
		    break;
		case S_HELP_ONLY:
		    AddList( ListHelpOnly, items->entry );
		    XtVaSetValues( items->entry, XtNsensitive, False, 0 );
		    break;
		case S_ALWAYS:
		    break;
	    }
	}
    }
}

static void MenuSelected ( w, item, dummy )
    Widget 	    w;
    TypeMenuItem *  item;
    caddr_t         dummy;
{
    if ( item->click != 0 )
	(*(item->click))(item);
    else
    {
	fputs( "Warning: menu item ", stderr   );
	fputs( XtName(w), stderr               );
	fputs( " has been selected.\n", stderr );
    }
}


/****************************************************************************
**
*F  UpdateMenus( <state> )  . . . . . .  update menus in case of state change
*/
void UpdateMenus ( state )
    Int         state;
{
    TypeList	l;
    Int         i;

    /* menu entry active only in break loop */
    l = ListErrorOnly;
    for ( i = 0;  i < l->len;  i++ )
    {
	if ( state == GAP_ERROR )
	    XtVaSetValues( (Widget)l->ptr[i], XtNsensitive, True, 0 );
	else
	    XtVaSetValues( (Widget)l->ptr[i], XtNsensitive, False, 0 );
    }

    /* menu entry active only during input */
    l = ListInputOnly;
    for ( i = 0;  i < l->len;  i++ )
    {
	if ( state == GAP_ERROR || state == GAP_INPUT )
	    XtVaSetValues( (Widget)l->ptr[i], XtNsensitive, True, 0 );
	else
	    XtVaSetValues( (Widget)l->ptr[i], XtNsensitive, False, 0 );
    }

    /* menu entry active only during normal input */
    l = ListNormalOnly;
    for ( i = 0;  i < l->len;  i++ )
    {
	if ( state == GAP_INPUT )
	    XtVaSetValues( (Widget)l->ptr[i], XtNsensitive, True, 0 );
	else
	    XtVaSetValues( (Widget)l->ptr[i], XtNsensitive, False, 0 );
    }

    /* menu entry active only while gap is running */
    l = ListRunningOnly;
    for ( i = 0;  i < l->len;  i++ )
    {
	if ( state == GAP_RUNNING )
	    XtVaSetValues( (Widget)l->ptr[i], XtNsensitive, True, 0 );
	else
	    XtVaSetValues( (Widget)l->ptr[i], XtNsensitive, False, 0 );
    }

    /* menu entry active only while gap is helping */
    l = ListHelpOnly;
    for ( i = 0;  i < l->len;  i++ )
    {
	if ( state == GAP_HELP )
	    XtVaSetValues( (Widget)l->ptr[i], XtNsensitive, True, 0 );
	else
	    XtVaSetValues( (Widget)l->ptr[i], XtNsensitive, False, 0 );
    }
}


/* * * * * * * * * * * * * * * gap talk window * * * * * * * * * * * * * * */


/****************************************************************************
**

*F  UpdateMemoryInfo( <type>, <val> )	. . . . . . . . . update memory usage
*/
static Widget LabelBytesFree;
static Widget LabelBytesTotal;
static Widget LabelHandleUsed;
static Widget LabelHandleDead;

void UpdateMemoryInfo ( type, val )
    Int		type;
    Int         val;
{
    char        tmp[30];

    
    switch ( type )
    {
	case 1:
	    sprintf( tmp, "Handle used: %-5d ", val );
  	    XtVaSetValues( LabelHandleUsed, XtNlabel, tmp, NULL );
	    break;
	case 2:
            sprintf( tmp, "Handle dead: %-5d ", val );
            XtVaSetValues( LabelHandleDead, XtNlabel, tmp, NULL );
	    break;
	case 3:
            sprintf( tmp, "KBytes free: %-4d ", val/1024 );
            XtVaSetValues( LabelBytesFree, XtNlabel, tmp, NULL );
	    break;
	case 4:
            sprintf( tmp, "KBytes total: %-4d ", val/1024 );
            XtVaSetValues( LabelBytesTotal, XtNlabel, tmp, NULL );
	    break;
    }
}


/****************************************************************************
**
*F  GapTalkResized( <talk>, <cd>, <evt>, <ctd> )  . . . . . . resize callback
*/
static void GapTalkResized P(( Widget, XtPointer, XEvent*, Boolean* ));

static void GapTalkResized ( talk, cd, evt, ctd )
    Widget	            talk;
    XtPointer               cd;
    XEvent                * evt;
    Boolean               * ctd;
{
    Int                     i;
    String                  ptr;
    Widget                  snk;
    XFontStruct           * font;
    char                    buf[128];
    static UInt             h,  h1 = 0;
    static UInt             w,  w1 = 0;

    /* is this a resize event */
    if ( evt->type == ConfigureNotify )
    {
	
	/* compute a sensible size */
	XtVaGetValues( talk, XtNtextSink, &snk,  0 );
	XtVaGetValues( snk,  XtNfont,     &font, 0 );
	w = evt->xconfigure.width / font->max_bounds.width - 3;
	h = evt->xconfigure.height / ( font->max_bounds.ascent
	    + font->max_bounds.descent ) - 2;
	if ( w < 2 )  w = 2;
	if ( h < 2 )  h = 2;
	if ( w == w1 && h == h1 )
	    return;
	w1 = w;
	h1 = h;

	/* construct gap command */
	strcpy( buf, "SizeScreen([ " );
	ptr = buf + strlen(buf);
	for ( i = 3;  0 <= i;  i--, w = w / 10 )
	    ptr[i] = w%10 + '0';
	ptr += 4;
	*ptr++ = ',';
	*ptr++ = ' ';
	for ( i = 3;  0 <= i;  i--, h = h / 10 )
	    ptr[i] = h%10 + '0';
        ptr += 4;
        strcpy( ptr, " ]);;\n" );

	/* if gap is waiting for input, do it */
	if ( GapState == GAP_INPUT || GapState == GAP_ERROR )
	    SimulateInput( buf );
	else
	    strcpy( ScreenSizeBuffer, buf );
    }
}


/****************************************************************************
**
*F  CreateGapWindow() . . . . . . . . . . . . . . create communication window
**
**  RESOURCES
**    *xgapMenu*shapeStyle
**        style of the menu buttons, default "ShapeOval"
**    *xgap.height
**    *xgap.width
**        start size of the communication text window
*/
static void CreateGapWindow P(( void ));

static void CreateGapWindow ()
{
    Widget	paned;
    Widget      box;
    Widget      button;
    Pixmap      symbol;
    Display   * display;
    Boolean     flag;
    Int         i;

    /* create a "paned" for the menu and text window */
    paned = XtVaCreateManagedWidget( "paned", panedWidgetClass,
				      XGap, NULL );

    /* create a menu box for the menu buttons */
    box = XtVaCreateManagedWidget( "xgapMenu", boxWidgetClass,
	      paned,
	      XtNx,	    	        0,
	      XtNy,                     0,
 	      XtNresizeToPreferred,	True,
	      NULL );

    /* create a menu button drop down symbol */
    display = XtDisplay(box);
    symbol = XCreateBitmapFromData( display,
	     DefaultRootWindow(display),
             "\376\3\2\2\2\6\162\6\2\6\162\6\2\6\162\6\2\6\2\6\376\7\370\7",
	     12, 12 );

    /* create file menu button and file menu */
    button = XtVaCreateManagedWidget( "gapButton", menuButtonWidgetClass,
	         box,
	         XtNleftBitmap,     symbol,
		 XtNx,              0,
		 NULL );
    CreateMenu( button, GapMenu );

    /* create run menu button and run menu */
    button = XtVaCreateManagedWidget( "runButton", menuButtonWidgetClass,
	         box,
		 XtNleftBitmap,     symbol,
		 XtNx,              10,
		 0 );
    CreateMenu( button, RunMenu );

    /* create help menu button and help menu */
    button = XtVaCreateManagedWidget( "helpButton", menuButtonWidgetClass,
	         box,
		 XtNleftBitmap,     symbol,
		 XtNx,              10,
		 0 );
    CreateMenu( button, HelpMenu );

    /* create the communication window */
    GapTalk = XtVaCreateManagedWidget( "xgapTalk", gapTextWidgetClass,
                  paned,
		  XtNinputCallback,    KeyboardInput,
		  XtNcheckCaretPos,    CheckCaretPos,
                  XtNscrollHorizontal, XawtextScrollWhenNeeded,
                  XtNscrollVertical,   XawtextScrollAlways,
                  XtNeditType,         XawtextEdit,
                  XtNbottomMargin,     15,
                  XtNx,                0,
                  XtNy,                10,
                  XtNdisplayCaret,     True,
                  NULL );
    XtAddEventHandler(GapTalk,StructureNotifyMask,False,GapTalkResized,0);
    GTDropGapPrompt( GapTalk, True );

    /* to quit or not do quit on CTR-D */
    XtVaGetValues( GapTalk, XtNquitGapCtrD, &flag,  0 );
    if ( flag )
    {
	for ( i = 0;  GapMenu[i].label;  i++ )
	    if ( !strcmp( GapMenu[i].label, "quitGapCTRD" ) )
		break;
	if ( GapMenu[i].label && GapMenu[i].click )
	    GapMenu[i].click(&(GapMenu[i]));
    }

    /* paste GAP prompt into talk window? */
    XtVaGetValues( GapTalk, XtNpasteGapPrompt, &flag,  0 );
    if ( flag )
    {
	for ( i = 0;  GapMenu[i].label;  i++ )
	    if ( !strcmp( GapMenu[i].label, "pastePrompt" ) )
		break;
	if ( GapMenu[i].label && GapMenu[i].click )
	    GapMenu[i].click(&(GapMenu[i]));
    }

    /* create a box and labels for garbage info */
    box = XtVaCreateManagedWidget( "xgapInfo", boxWidgetClass,
	      paned,
	      XtNx, 	                0,
	      XtNy,                     20,
	      XtNskipAdjust,            True,
	      XtNresizeToPreferred, 	True,
	      NULL );
    LabelBytesFree = XtVaCreateManagedWidget( "bytesFree",
			 labelWidgetClass, box,
			 XtNborderWidth, 0,
			 NULL );
    LabelBytesTotal = XtVaCreateManagedWidget( "bytesTotal",
		          labelWidgetClass, box,
 			  XtNborderWidth, 0,
			  NULL );
    LabelHandleUsed = XtVaCreateManagedWidget( "handleUsed",
		          labelWidgetClass, box,
			  XtNborderWidth, 0,
			  NULL );
    LabelHandleDead = XtVaCreateManagedWidget( "handleDead",
			  labelWidgetClass, box,
			  XtNborderWidth, 0,
			  NULL );
    UpdateMemoryInfo( 1, 0 );
    UpdateMemoryInfo( 2, 0 );
    UpdateMemoryInfo( 3, 0 );
    UpdateMemoryInfo( 4, 0 );
}


/* * * * * * * * * * * * * * * error handler * * * * * * * * * * * * * * * */


/****************************************************************************
**

*F  MyErrorHandler(<dis>) . . . . . . . . . . . . kill gap in case of X error
*/
static int (*OldErrorHandler)();

static int MyErrorHandler ( dis, evt )
    Display       * dis;
    XErrorEvent	    evt;
{
#   ifdef DEBUG_ON
        fputs( "killing gap because of X error\n", stderr );
#   endif
    KillGap();
    return OldErrorHandler( dis, evt );
}


/****************************************************************************
**
*F  MyIOErrorHandler(<dis>) . . . . . . . . . . . kill gap in case of X error
*/
static int (*OldIOErrorHandler)();

static int MyIOErrorHandler ( dis )
    Display   * dis;
{
#   ifdef DEBUG_ON
        fputs( "killing gap because of X IO error\n", stderr );
#   endif
    KillGap();
    return OldIOErrorHandler(dis);
}


/****************************************************************************
**
*F  MySignalHandler() . . . . . . . . . . . . . .  kill gap in case of signal
*/
#ifdef DEBUG_ON

static void (*OldSignalHandlerHUP)();
static void (*OldSignalHandlerINT)();
static void (*OldSignalHandlerQUIT)();
static void (*OldSignalHandlerILL)();
static void (*OldSignalHandlerIOT)();
static void (*OldSignalHandlerBUS)();
static void (*OldSignalHandlerSEGV)();

static void MySignalHandlerHUP ()
{
    fputs( "killing gap because of signal HUP\n", stderr );
    KillGap();
    OldSignalHandlerHUP();
    exit(1);
}
static void MySignalHandlerINT ()
{
    fputs( "killing gap because of signal INT\n", stderr );
    KillGap();
    OldSignalHandlerINT();
    exit(1);
}
static void MySignalHandlerQUIT ()
{
    fputs( "killing gap because of signal QUIT\n", stderr );
    KillGap();
    OldSignalHandlerQUIT();
    exit(1);
}
static void MySignalHandlerILL ()
{
    fputs( "killing gap because of signal ILL\n", stderr );
    KillGap();
    OldSignalHandlerILL();
    exit(1);
}
static void MySignalHandlerIOT ()
{
    fputs( "killing gap because of signal IOT\n", stderr );
    KillGap();
    OldSignalHandlerIOT();
    exit(1);
}
static void MySignalHandlerBUS ()
{
    fputs( "killing gap because of signal BUS\n", stderr );
    KillGap();
    OldSignalHandlerBUS();
    exit(1);
}

static void MySignalHandlerSEGV ()
{
    fputs( "killing gap because of signal SEGV\n", stderr );
    KillGap();
    OldSignalHandlerSEGV();
    exit(1);
}

#else

static void MySignalHandler ()
{
    KillGap();
    exit(1);
}

#endif


/* * * * * * * * * * * * * * * * main program  * * * * * * * * * * * * * * */


/****************************************************************************
**

*F  ParseArgs( <argc>, <argv> ) . . . create argument arry for gap subprocess
*/
static char * nargv[1024];

static void ParseArgs ( argc, argv )
    Int         argc;
    char     ** argv;
{
    Int		nargc;
    Int         i,  j;
    String      p;

    /* at first assume that "gap" is started with "gap", append "-p" */
    nargc = 0;
    nargv[nargc++] = "gap";
    nargv[nargc++] = "-p";

    /* parse arguments, collect arguments for gap in <nargv> */
    for ( argv++, argc--;  0 < argc;  argv++, argc-- )
    {
	if ( *argv[0] == '-' )
	{

	    /* don't group any options */
	    if ( strlen(*argv) != 2 )
	    {
		fputs("Gap: sorry options must not be grouped '", stderr);
		fputs(*argv, stderr);
		fputs("'.\n", stderr);
		exit(1);
	    }
	    switch( argv[0][1] )
	    {

		/* toggle debug */
		case 'D':
#ifdef DEBUG_ON
		    if ( argc-- < 0 )
		    {
			fputs( "XGAP: option '-D' must have an argument.\n",
			       stderr );
		    }
		    Debug = atoi(*++argv);
#else
		    fputs( "XGAP: compile XGAP using 'COPTS=-DDEBUG_ON'.\n",
			   stderr );
#endif
		    break;

                /* broken window manager */
		case 'W':
		    PopupDialogBrokenWM();
		    break;

		/* copy GAP output to stderr */
	        case 'E':
		    SpyMode = !SpyMode;
		    break;

		/* get name of gap subprocess */
		case 'G':
		    if ( argc-- < 2 )
		    {
			fputs( "XGAP: option '-G' must have an argument.\n",
			       stderr );
			exit(1);
		    }
		    p = *++argv;
		    nargv[0] = p;
		    j = 0;
		    while ( *++p )
		    {
			if ( *p == ' ' )
			{
			    *p = '\0';
			    j  = j + 1;
			}
		    }
		    if ( 0 < j )
		    {
			for ( i = nargc-1;  0 < i;  i-- )
			    nargv[i+j] = nargv[i];
			nargc = nargc + j;
			while ( 0 < j )
			{
			    while ( *--p ) ;
			    nargv[j--] = p+1;
			}
		    }
		    break;

		/* let GAP handle this argument */
		default:
		    nargv[nargc++] = *argv;
		    break;
	    }
	}
	else
	    break;
    }

    /* copy any remaining arguments */
    for ( ; 0 < argc;  argv++, argc-- )
	nargv[nargc++] = *argv;
    nargv[nargc] = 0;
    return;
}


/****************************************************************************
**
*F  main( <argc>, <argv> )  . . . . . . . . . . . . . . . .   main event loop
*/
#include "checksym.bm"
#include "emptymk.bm"
#include "exmark.bm"
#include "menusym.bm"

int main ( argc,  argv )
    int         argc;
    char     ** argv;
{
    String      color;
    String      colors;
    int         fromGap;
    Int         mod = -1;
    Int         len;


    /* create a new top level shell and an applictation context */
    XGap = XtVaAppInitialize( &AppContext, "XGap",
			      CommandOptions, XtNumber(CommandOptions),
			      &argc, argv, FallbackResources, 0 );
    GapDisplay   = XtDisplay(XGap);
    GapScreen    = DefaultScreen(GapDisplay);
    MyRootWindow = RootWindow( GapDisplay, GapScreen );

    /* create top left arrow cusor */
    CursorTL = XCreateFontCursor( XtDisplay(XGap), XC_top_left_arrow );

    /* create menu symbol */
    MenuSymbol = XCreateBitmapFromData( GapDisplay, MyRootWindow,
				        menusym_bits, menusym_width,
				        menusym_height );

    /* create check mark and empty mark */
    CheckMarkSymbol = XCreateBitmapFromData( GapDisplay, MyRootWindow,
					     checksym_bits, checksym_width,
					     checksym_height );
    EmptyMarkSymbol = XCreateBitmapFromData( GapDisplay, MyRootWindow,
					     emptymk_bits, emptymk_width,
					     emptymk_height );

    /* exclamation mark */
    ExMarkSymbol = XCreateBitmapFromData( GapDisplay, MyRootWindow,
					  exmark_bits, exmark_width,
					  exmark_height );

    /* WM_DELETE_WINDOW atom */
    WmDeleteWindowAtom = XInternAtom(GapDisplay, "WM_DELETE_WINDOW", False);

    /* install our error handler, we have to kill gap in this case */
    OldIOErrorHandler = XSetIOErrorHandler( MyIOErrorHandler );
    OldErrorHandler   = XSetErrorHandler  ( MyErrorHandler   );

    /* install your signal handler, we have to kill gap in this case */
#   ifdef DEBUG_ON
        OldSignalHandlerHUP  = signal( SIGHUP,  MySignalHandlerHUP );
        OldSignalHandlerINT  = signal( SIGINT,  MySignalHandlerINT );
        OldSignalHandlerQUIT = signal( SIGQUIT, MySignalHandlerQUIT );
        OldSignalHandlerILL  = signal( SIGILL,  MySignalHandlerILL );
        OldSignalHandlerIOT  = signal( SIGIOT,  MySignalHandlerIOT );
        OldSignalHandlerBUS  = signal( SIGBUS,  MySignalHandlerBUS );
        OldSignalHandlerSEGV = signal( SIGSEGV, MySignalHandlerSEGV );
#   else
        signal( SIGHUP,  MySignalHandler );
        signal( SIGINT,  MySignalHandler );
        signal( SIGQUIT, MySignalHandler );
        signal( SIGILL,  MySignalHandler );
        signal( SIGIOT,  MySignalHandler );
        signal( SIGBUS,  MySignalHandler );
        signal( SIGSEGV, MySignalHandler );
#   endif

    /* create the gap talk window */
    CreateGapWindow();
    XtRealizeWidget(XGap);

    /* initialize window commands */
    InitXCMDS();

    /* get color model */
    XtVaGetValues( GapTalk, XtNcolorModel, &color, 0 );
    len = strlen(color);
    if ( !strncmp( color, "black&white", len ) )
	mod = CM_BW;
    else if ( !strncmp( color, "monochrome", len ) )
	mod = CM_BW;
    else if ( !strncmp( color, "grey", len ) )
	mod = CM_GRAY;
    else if ( !strncmp( color, "gray", len ) )
	mod = CM_GRAY;
    else if ( !strncmp( color, "color5", len ) )
	mod = CM_COLOR5;
    else if ( !strncmp( color, "color3", len ) )
	mod = CM_COLOR3;
    else if ( !strncmp( color, "default", len ) )
	mod = -1;
    else
    {
	fputs( "Gap: unkown color model '", stderr );
	fputs( color, stderr );
	fputs( "'\n", stderr );
	exit(1);
    }
    if ( mod != -1 ) {
	XtVaGetValues( GapTalk, XtNcolors, &colors, 0 );
	GCSetColorModel( GapDisplay, mod, colors );
    }

    /* parse remaining arguments */
    ParseArgs( argc, argv );
    fromGap = StartGapProcess( nargv[0], nargv );

    /* add callback for output from gap*/
    XtAppAddInput( AppContext,  fromGap,  (XtPointer) XtInputReadMask,
                   (XtInputCallbackProc) GapOutput,  (XtPointer) 0 );

    /* force a garbage collection in the beginning */
    StoreInput( "GASMAN(\"collect\");\n", 19 );

    /* talk window is drawn only partial during start up otherwise (why?) */
    XFlush( GapDisplay );
    sleep(1);
    XFlush( GapDisplay );
    sleep(1);
    XFlush( GapDisplay );
    sleep(1);

    /* enter main read-eval loop */
    XtAppMainLoop(AppContext);
    return 0;
}

