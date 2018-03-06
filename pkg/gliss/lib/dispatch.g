###############################################################################
##
##  dispatch.g  GLISSANDO ver 1.0  Christof Noebauer   1995, 1996
##                                                       
#############################################################################
##  Dispatchers
#############################################################################
##
#F  DisplayCayleyTable( <D> ). . .nicely print a Cayley table of a domain <D> 
##
##  Implemented for <D> a semigroup or a group or a nearring.
##
DisplayCayleyTable:=Dispatcher("DisplayTable");

#############################################################################
##
#F  IsCommutative( <D> ) . . . . . . . . . test commutativity of a domain <D>
##
##  Dispatcher function to test if a domain <D> is commutative where <D> is
##  a nearring or a semigroup.
##
IsCommutative:=AttributeDispatcher("IsCommutative");
