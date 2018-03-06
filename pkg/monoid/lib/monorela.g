#############################################################################
##
#A  GAP                                                 Goetz.Pfeiffer@UCG.IE
##
#A  $Id: monorela.g,v 2.0 1997/05/05 16:35:40 goetz Exp $
##
#Y  Copyright (C) 1997, Mathematics Dept, University College Galway, Ireland.
##
##  This file defines the functions for monoids of binary relations.
##

#############################################################################
##
#F  InfoMono? . . . . . . . . . . . . . . . . . . . . . . . . info functions.
##
if not IsBound(InfoMono1) then InfoMono1:= Ignore; fi;
if not IsBound(InfoMono2) then InfoMono2:= Ignore; fi;

#############################################################################
##
#V  RelMonoidOps . . . . . . . . . . . . . . . . . . . . . operations record.
##
##  This is the operations record for monoids of binary relations.  This is a
##  special kind of monoid, thus 'RelMonoidOps' inherits from 'MonoidOps'.
##
RelMonoidOps:= OperationsRecord("RelMonoidOps", MonoidOps);

#############################################################################
##
#F  RelMonoidOps.Degree( <M> ) . . . . . . . . . . . . . . . . . . .  degree.
##
##  The *degree* of a monoid  of binary relations  is the number of points it
##  acts upon.
##
RelMonoidOps.Degree:= M-> Degree(M.identity);

#############################################################################
##
#E  Emacs  . . . . . . . . . . . . . . . . . . . . . . local emacs variables.
##
##  Local Variables:
##  mode:               outline
##  outline-regexp:     "#F\\|#V\\|#A\\|#E"
##  fill-column:        77
##  fill-prefix:        "##  "
##  eval:               (hide-body)
##  End:
##
