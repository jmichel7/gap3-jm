#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @ (#)$Id: ctensorp.g,v 1.1 1997/03/10 13:52:27 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: ctensorp.g,v $
#H  Revision 1.1  1997/03/10 13:52:27  gap
#H  VERSION 1.0
#H
#H  Revision 1.2  1997/01/05 10:49:23  fceller
#H  added Eamonn's new version to the reprository
#H
#H  Revision 1.1  1996/12/25 09:03:39  fceller
#H  changed long filenames to MS-DOS conform filenames,
#H  the init files are *NOT* yet updated
#H
#H  Revision 1.1  1996/11/28 13:14:42  fceller
#H  added "smash" and "reducible" to the repository
#H
##
#c_tensorprod.g
#
###############################################################################
##
#F  TensorProductGModule ( m1, m2 ) . . . . . tensor product of two G-modules
##
## TensorProductGModule calculates the tensor product of modules m1 and m2. 
## They are assumed to be modules for the same group so, in particular, they
## should have the same number of generators.
## 
TensorProductGModule := function ( m1, m2)

   local mat1, mat2, F1, F2,  gens, i, l;

   mat1 := GeneratorsFlag (m1); mat2 := GeneratorsFlag (m2);
   F1 := FieldFlag (m1); F2 := FieldFlag (m2);
   if (F1 <> F2) then
      Error ("GModules are defined over different fields.\n");
   fi;
   l := Length (mat1);
   if (l <> Length (mat2)) then
      Error ("GModules have different numbers of generators.");
   fi;

   gens := [];
   for i in [1..l] do
      gens[i] := KroneckerProduct (mat1[i], mat2[i]);
   od;

   return GModule (gens, F1);
end;

###############################################################################
##
#F  WedgeGModule ( module ) . . . . . wedge product of a G-module
##
## WedgeGModule calculates the wedge product of a G-module.
## That is the action on antisymmetrix tensors.
## 
WedgeGModule := function ( module)
   local mats, mat, newmat, row, F, gens, dim, nmats, i, j, k, m, n, x;

   mats := GeneratorsFlag (module);
   F := FieldFlag (module);
   nmats := Length (mats);
   dim := Length (mats[1]);

   gens := [];
   for i in [1..nmats] do
      mat := mats[i];
      newmat := [];
      for j in [1..dim] do
         for k in [1..j - 1] do
            row := [];
            for m in [1..dim] do
               for n in [1..m - 1] do
                  x := mat[j][m] * mat[k][n] - mat[j][n] * mat[k][m];
                  Add (row, x);
               od;
            od;
            Add (newmat, row);
         od;
      od;
      Add (gens, newmat);
   od;

   return GModule (gens, F);
end;
