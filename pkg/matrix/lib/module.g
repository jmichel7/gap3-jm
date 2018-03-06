#############################################################################
##
#W  module.g	 	   	Matrix Packages                    Derek Holt
#W                                                   &  Charles Leedham-Green
#W                                                           & Eamonn O'Brien
#W                                                               & Sarah Rees 
##
#H  @(#)$Id: module.g,v 1.1 1997/03/10 13:52:05 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  The function GModule  creates a  module as a    record, out of  a set  of
##  matrices,  or a group (generated  by a set  of matrices).  The compulsory
##  flags .field, .dimension, .generators  and .isGmodule are set at creation
##  time. Other optional flags appear below.  The field  of the module may be
##  specified as an argument.   If not the field  is  set to be the  smallest
##  field over which the matrices are defined.
##
RevisionMatrix.recclass_g :=
    "@(#)$Id: module.g,v 1.1 1997/03/10 13:52:05 gap Exp $";


GModuleOps:=OperationsRecord("GModuleOps");

GModuleOps.Identity:=G->GeneratorsFlag(G)[1]^0;

GModuleOps.IsIrreducible:=M->MatrixMTX.isIrreducible(M);

GModuleOps.IsPrimitive:=function(arg)
  return ApplyFunc(MatrixMTX.IsPrimitive,arg);end;

GModuleOps.CompositionFactors:=M->MatrixMTX.CompositionFactors(M);
            
GModuleOps.IsAbsolutelyIrreducible:=M->MatrixMTX.IsAbsolutelyIrreducible(M);

#############################################################################
##

#F  GModule( <matrices> [,<F>] )
##
##  Make  a  module from  a list  of matrices  or  from a  matrix  group. The
##  underlying field may be specified as an an optional argument.
##
GModule := function( arg )
    local   matrices,  deg,  i,  j,  char,  F,  dimension;

    if Number(arg) = 1  then
        if IsGroup(arg[1])  then 
            matrices := Generators (arg[1]);
        else 
            matrices := arg[1];
        fi;
        deg := 1;
        for i  in matrices  do
            for j  in i  do
                deg := LcmInt( deg, DegreeFFE(j) );
            od;
        od;
        char := CharFFE(matrices[1][1][1]);
        F := GF(char^deg);
    elif Number (arg) = 2  then
        if IsGroup (arg[1]) then 
            matrices := Generators(arg[1]);
        else 
            matrices := arg[1];
        fi;
        F := arg[2];
    else
        Error("usage: GModule( <matrices> [,<F>] )");
    fi;
    if 0 = Length(matrices)  then
        if IsRec(arg[1]) and IsBound(arg[1].identity)  then
            dimension := Length(arg[1].identity);
        elif IsRec(arg[1]) and IsBound(arg[1].dimension)  then
            dimension := arg[1].dimension;
        else
            Error( "dimension is not known" );
        fi;
    else
        dimension := Length(matrices[1]);
    fi;
    for i  in matrices  do
        if not IsMatrix(i)  then
            Error( "<matrices> must be a list of matrices" );
        fi;
        if Length(i) <> dimension  then
            Error( "<matrices> must be a list of ", dimension,
                   "dimensional matrices" );
        fi;
    od;

    return rec( field := F, 
                dimension := dimension,
                generators := matrices, 
                isGModule := true,
		operations:=GModuleOps);
end;


#############################################################################
##
#F  IsGModule( <module> )
##
IsGModule := function (module)
   if IsRec (module) = false or IsBound (module.isGModule) = false then
      return false;
   fi;
   return module.isGModule;
end;


#############################################################################
##
#F  DualFrobeniusGModule( <module> )
##
DualFrobeniusGModule := function( module )
    local   F,  k,  dim,  mats,  dmats,  qq,  i,  j,  l;

    F := FieldFlag(module);
    k := LogInt( Size(F), Characteristic(F) );
    if k mod 2 = 1  then
        Error( "field <F> is not a square" );
    fi;
    dim   := DimensionFlag(module);
    mats  := GeneratorsFlag(module);
    dmats := Copy(mats);
    qq    := Characteristic(F) ^ ( k / 2 );
    for i  in [ 1 .. Length(mats) ]  do
        for j  in [ 1 .. dim ]  do
            for l  in [ 1 .. dim ]  do
                dmats[i][j][l] := mats[i][l][j]^qq;
            od;
        od;
    od;
    return GModule( List( dmats, x -> x^-1 ), F );
end;



#############################################################################
##

#E  module.g  . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
##
