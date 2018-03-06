/****************************************************************************
**
*A  gapgraph.h                  XGAP source                      Frank Celler
**
*H  @(#)$Id: gapgraph.h,v 1.1.1.1 1996/12/11 12:40:08 werner Exp $
**
*Y  Copyright 1995-1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
**
*H  $Log: gapgraph.h,v $
*H  Revision 1.1.1.1  1996/12/11 12:40:08  werner
*H  Preparing 3.4.4 for release
*H
*H  Revision 1.6  1995/07/24  09:28:30  fceller
*H  reworked most parts to use nice typedefs
*H
*H  Revision 1.5  1995/01/24  16:29:14  fceller
*H  'GCSetColorModel' has a new parameter which allows to set the colors
*H
*H  Revision 1.4  1995/01/24  12:58:52  fceller
*H  added color support
*H
*H  Revision 1.3  1994/06/06  08:57:24  fceller
*H  added database
*H
*H  Revision 1.2  1993/10/20  12:52:02  fceller
*H  fixed prototype for non stdc
*H
*H  Revision 1.1  1993/10/18  11:04:47  fceller
*H  added fast updated,  fixed timing problem
*H
*H  Revision 1.0  1993/04/05  11:42:18  fceller
*H  Initial revision
*/
#ifndef _gapgraph_h
#define _gapgraph_h


/****************************************************************************
**
*D  T_XXX . . . . . . . . . . . . . . . . . . . . . . . .  the graphic object
*/
#define T_LINE          1
#define T_CIRCLE        2
#define T_DISC          3
#define T_TEXT          4
#define T_RECT          5
#define T_BOX           6


/****************************************************************************
**
*D  C_XXX . . . . . . . . . . . . . . . . . . . . . . . . . .  the color XXXX
*/
#define C_BLACK		0
#define C_WHITE		1
#define C_LIGHT_GRAY	2
#define C_DARK_GRAY     3
#define C_LAST_GRAY     3
#define C_RED		4
#define C_BLUE		5
#define C_GREEN		6
#define C_LAST          6


/****************************************************************************
**
*D  CM_BW . . . . . . . . . . . . . . . . . . . . .  color model: black/white
*D  CM_GRAY . . . . . . . . . . . . . . . . . . . . . .  color model: 2 grays
*D  CM_COLOR3 . . . . . . . . . . . . . . . . . . . . . color model: 3 colors
*D  CM_COLOR5 . . . . . . . . . . . . . . . .  color model: 3 colors, 2 grays
*/
#define CM_BW		1
#define CM_GRAY		2
#define CM_COLOR3       3
#define CM_COLOR5       4


/****************************************************************************
**

*T  TypeGapGraphicObject  . . . . . . . . . . . . .  graphic object of widget
*/
typedef struct _gap_graphic_obj
{
    Short   type;
    Short   color;
    Int     x, y, w, h;
    union
    {
        struct { Int x1, x2, y1, y2, w;                }   line;
        struct { Int x, y, r, w;                       }   circle;
        struct { Int x, y, r;                          }   disc;
        struct { Int x, y, len; String str; Font font; }   text;
	struct { Int x1, x2, y1, y2, w;                }   rect;
	struct { Int x1, x2, y1, y2;                   }   box;
    } desc;
}
TypeGapGraphicObject;


/****************************************************************************
**
*T  GapGrahpicClassRec	. . . . . . . . . . . . . .  gap graphic class record
*/
typedef struct {int empty;} GapGraphicClassPart;

typedef struct _GapGraphicClassRec
{
    CoreClassPart       core_class;
    GapGraphicClassPart gap_graphic_class;
}
GapGraphicClassRec;

extern GapGraphicClassRec gapGraphicClassRec;


/****************************************************************************
**
*T  GapGarphicRec . . . . . . . . . . . . . . . . . gap graphic widget record
*/
typedef struct {

    /* dimension of window */
    UInt                width;
    UInt                height;

    /* reference number */
    Int                 number;

    /* list of graphic objects */
    TypeList            objs;

    /* display information */
    Display           * display;
    Pixel               black;
    Pixel               white;
    GC                  gc;

    /* bounding box to update */
    Boolean             update;
    Boolean             fast_update;
    Int                 lx,  hx;
    Int                 ly,  hy;
}
GapGraphicPart;

typedef struct _GapGraphicRec
{
    CorePart            core;
    GapGraphicPart      gap_graphic;
}
GapGraphicRec;


/****************************************************************************
**
*T  GapGraphicWidgetClass . . . . . . . . . . . . . . . . . .  class datatype
*/
typedef struct _GapGraphicClassRec    * GapGraphicWidgetClass;


/****************************************************************************
**
*T  GapGraphicWidget  . . . . . . . . . . . . . . . . . . . instance datatype
*/
typedef struct _GapGraphicRec * GapGraphicWidget;


/****************************************************************************
**
*V  gapGraphicWidgetClass . . . . . . . . . . . . . . . . .  class definition
*/
extern WidgetClass gapGraphicWidgetClass;


/****************************************************************************
**

*P  Prototypes  . . . . . . . . . . . prototypes of public gap text functions
*/
extern Int     GCColorModel P(( Display* ));
extern void    GCSetColorModel P(( Display*, Int, String ));
extern Int     GGAddObject P(( Widget, TypeGapGraphicObject* ));
extern void    GGFreeAllObjects P(( Widget ));
extern void    GGFreeGapGraphicObjects P(( Widget ));
extern void    GGFreeObject P(( TypeGapGraphicObject* ));
extern Boolean GGRemoveObject P(( Widget, Int ));
extern void    GGResize P(( Widget, Int, Int ));
extern void    GGStartRemove P(( Widget ));
extern void    GGStopRemove P(( Widget ));

#if defined(__GNUC__) && !defined(__GNUC_MINOR__)
extern void    GGDrawObject ();
extern void    GGFastUpdate ();
#else
extern void    GGDrawObject P(( Widget, TypeGapGraphicObject*, Boolean ));
extern void    GGFastUpdate P(( Widget, Boolean ));
#endif

#endif

