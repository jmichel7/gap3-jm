#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @ (#)$Id: cdualmod.g,v 1.1 1997/03/10 13:52:23 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: cdualmod.g,v $
#H  Revision 1.1  1997/03/10 13:52:23  gap
#H  VERSION 1.0
#H
#H  Revision 1.2  1997/01/05 10:49:18  fceller
#H  added Eamonn's new version to the reprository
#H
#H  Revision 1.1  1996/12/25 09:03:31  fceller
#H  changed long filenames to MS-DOS conform filenames,
#H  the init files are *NOT* yet updated
#H
#H  Revision 1.1  1996/11/28 13:14:40  fceller
#H  added "smash" and "reducible" to the repository
#H
##
#dualmod.g
#
###############################################################################
##
#F  DualGModule ( module ) . . . . . dual of a G-module
##
## DualGModule computes the dual of the module m1. 
## This is done simply by transposing and inverting the matrices
## 
DualGModule := function ( module )
    local i, mats, dmats, dmodule, F; 

   F := FieldFlag (module);
   mats := GeneratorsFlag (module);
   dmats := [];
   for i in [1..Length (mats)] do
     dmats[i] := TransposedMat (mats[i]^-1);
   od;
   dmodule := GModule (dmats, F);

   return dmodule;
end;
