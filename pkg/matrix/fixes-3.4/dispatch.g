#############################################################################
##
#W  fixes_dispatch.g	   	Matrix Packages                  Frank Celler
#
#H  @(#)$Id: dispatch.g,v 1.1 1997/03/10 13:51:40 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file  fixes functions from  the  file "dispatch.g".
##
RevisionMatrix.fixes_dispatch_g :=
    "@(#)$Id: dispatch.g,v 1.1 1997/03/10 13:51:40 gap Exp $";

# read in the files to patch
#############################################################################
##

#F  BlowupVec( <V>, <v> )
##
BlowupVec := function( V, v )
    BlowupVec := CRecSL.BlowupVec;
    return BlowupVec( V, v );
end;


#############################################################################
##
#F  Rewrite( <struct>, <obj> )
##
Rewrite := function( struct, obj )
    return struct.operations.Rewrite( struct, obj );
end;


#############################################################################
##
#F  IsPrimitive( <obj> )  . . . . . . . . . . . . check if <obj> is primitive
##
MatGroupOps.IsPrimitive := function ( arg )
     return ApplyFunc(MatrixMTX.IsPrimitive,arg);
end;
