/*****************************************************************************
*                                                                            *
* This is the header file for versions 1.6 of nautaux.c.                     *
*                                                                            *
*   Copyright (1984-1989) Brendan McKay.  All rights reserved.               *
*   Subject to the waivers and disclaimers in nauty.h.                       *
*                                                                            *
*   CHANGE HISTORY                                                           *
*       26-Apr-89 : initial creation for version 1.5.                        *
*       14-Oct-90 : renamed as version 1.6 (no changes to this file)         *
*                                                                            *
*****************************************************************************/

#include "nauty.h"           /* which includes stdio.h */

EXTPROC(boolean equitable,(graph*,nvector*,nvector*,int,int,int))
EXTPROC(long ptncode,(graph*,nvector*,nvector*,int,int,int))
EXTPROC(int component,(graph*,int,set*,int,int))
