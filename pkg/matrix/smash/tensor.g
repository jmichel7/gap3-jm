#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @ (#)$Id: tensor.g,v 1.1 1997/03/10 13:52:46 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: tensor.g,v $
#H  Revision 1.1  1997/03/10 13:52:46  gap
#H  VERSION 1.0
#H
#H  Revision 1.3  1997/01/10 11:48:47  fceller
#H  put Eamonn's (final) version into the repository,
#H  made a small fix to 'PermGroupRepresentation',
#H  made a small fix to alternating recognition
#H
#H  Revision 1.2  1997/01/05 10:49:40  fceller
#H  added Eamonn's new version to the reprository
#H
#H  Revision 1.1  1996/12/22 07:45:11  fceller
#H  added new files to the repository
#H
#H  Revision 1.2  1996/12/12 10:51:33  fceller
#H  new version by Eamonn
#H
##
############################################################################
##
#F  InfoTensorProduct  (...)  . . . . . . . . . . . for debugging assistance
##
##
if not IsBound (InfoTensorProduct)  then InfoTensorProduct := Ignore;  fi;
#############################################################################
##
#F  TensorProductDecomposition ( module, basis, dim1, dim2)  . . 
##  test for tensor product of modules of dimensions dim1, dim2.
##
TensorProductDecomposition := function  (module, basis, dim1, dim2) 

   local F, invbasis, g, x, factors, matrices, matrices1, matrices2;

   if IsGModule (module) = false then
      Error ("usage: TensorProductDecomposition (<module>, <basis>, <dim1>, <dim2>)");
   fi;

   if dim1 * dim2 <> DimensionFlag (module) then return false; fi;

   F := FieldFlag (module);
   invbasis := basis^-1;
   matrices1 := [];
   matrices2 := [];
   matrices := GeneratorsFlag (module);
   for g in matrices do
      x := basis * g * invbasis;
      factors := KroneckerFactors (x, dim1, dim2, F);
      if factors = false then 
         return false;
      else
         Add (matrices1, factors[1]);
         Add (matrices2, factors[2]);
      fi;
   od;

   SetTensorProductFlag (module, true);
   SetTensorBasisFlag (module, basis);
   SetTensorFactorsFlag (module, 
           [GModule (matrices1, F), GModule (matrices2, F)]);
   InfoTensorProduct ("#I Module is a tensor product of modules of dimension ", 
                      dim1, " and ", dim2, ".\n");
   return true;

end;

#############################################################################
##
#F  KroneckerFactors ( x, dim1, dim2 [,F]) . . test to see if x is a Kronecker
## product of 2 matrices of dimensions dim1, dim2 respectively, over F.  
## More precisely, we try to find A, a dim1 x dim1 matrix, and B, 
## a dim2 x dim2 matrix, so that x decomposes into dim1 x dim1 blocks, 
## with the k, l-th block equal to  A[k][l] * B, 
## i.e. x is the Kronecker product of A and B.
## If we can find such matrices we return the pair [A, B], 
## otherwise we return false.
## 
KroneckerFactors := function  ( arg ) 
   local zero, x, dim1, dim2, F, r, s, r0, s0, i, j, k, l, y, A, B;

   x := arg[1];

   dim1 := arg[2];
   dim2 := arg[3];
   if dim1 * dim2 <> Length (x) then return false; fi;

   if Number (arg) = 4 then 
      F := arg[4];
   else 
      F := Field (Flat (x));
   fi;

   zero := Zero (F);
   A := [];
   B := [];

   # first find a position where there's a non-zero entry
   i := 1; j := 1;
   while  x[i][j] = zero do
      if j < dim2 then j := j + 1; else j := 1; i := i + 1; fi;
   od;

   # so x[i][j] <> 0
   y := x[i][j];
   r := (i - 1) mod dim2 + 1; r0 := i - r;
   s := (j - 1) mod dim2 + 1; s0 := j - s;
   for i in [1..dim2] do
      B[i] := [];
      for j in [1..dim2] do
         B[i][j] := x[r0 + i][s0 + j];
      od;
   od;
   for k in [1..dim1] do
      A[k] := [];
      for l in [1..dim1] do
         A[k][l] := x[(k - 1) * dim2 + r][ (l - 1) * dim2 + s] / y;  
      od;
   od;

   if x <> KroneckerProduct (A, B) then 
      return false; 
   else 
      return [A, B];
   fi;

end;

#############################################################################
##
#F  UndoTensorProductFlags ( module)  . . undo flags set by TensorProduct 
##
## 
UndoTensorProductFlags := function  (module) 

   UndoTensorProductFlag (module);
   UndoTensorBasisFlag (module);
   UndoTensorFactorsFlag (module);

end;

#############################################################################
##
#F  SymTensorProductDecomposition ( module, Smodule)  . . 
## test to see if we have a symmetric tensor product
##
## module is a module for a finite matrix group G over a finite field.
## and Smodule is the module corresponding to the action of a subgroup 
## <S> of G on the same vector space. 
## G and <S> are assumed to act absolutely irreducibly.
## The function returns true if Smodule can be decomposed as a 
## tensor product of spaces all of the same dimension, and if further 
## these tensor factors are permuted by the action of G. 
## In that case components of the record module record the tensor 
## decomposition and the action of G permuting the factors. 
## If no such decomposition is found, the function returns false.
## Since the function uses random elements to try and find the tensor
## decomposition, a negative answer is *not* conclusive.
## 
SymTensorProductDecomposition := function  (module, Smodule) 

   local d, F, matrices, ngens, permutes, numTries, divisors, 
         poss, pair, dd, n, tenpow, basis, invbasis, permaction, 
         pi_ij, pi_g, factors, g, h, i, j, k;


   if IsGModule (module) = false then
      Error ("usage: SymTensorProductDecomposition (<module>, <smodule>)");
   fi;

   numTries := 20;  

   d := DimensionFlag (module);
   F := FieldFlag (module);
   matrices := GeneratorsFlag (module);
   ngens := Length (matrices);
   permaction := [];

   divisors := DivisorsInt (d);
   poss := [];
   for dd in divisors do
      if dd <> 1 and dd <> d then
         n := IntPower (d, dd);
         if n <> false then Add (poss, [dd, n]); fi;
      fi;
   od;
   if poss = [] then 
      InfoTensorProduct ("#I Dimension is not a proper power.\n");
      return false; 
   fi;

   for pair in poss do
      InfoTensorProduct ("#I Trying pair ", pair, " in SymTensorProductDecomposition.\n");
      dd := pair[1]; n := pair[2];
      tenpow := MultipleTensorProductDecomposition (Smodule, dd, n, numTries);
      if tenpow <> false then
         InfoTensorProduct ("#I Found a tensor power decomposition.\n");
         InfoTensorProduct ("#I GModule is ", n, "-th power of a ", 
                                dd, "-dimensional module.\n");
         basis := tenpow[1];
         invbasis := basis^-1;
         k := 1;
         permutes := true;
         while permutes and k <= ngens do
            g := basis * matrices[k] * invbasis;
            pi_g := ();
            i := 1;
            while permutes = true and i < n do
               j := i;
               factors := KroneckerFactors (g, dd^(n - i), dd, F);
               if factors = false then
                  repeat
                     j := j + 1; 
                     if j <= n then
                        pi_ij := SwapFactors (1, j + 1 - i, dd, n + 1 - i, F);  
                        factors := KroneckerFactors (g * pi_ij, dd^(n - i), dd, F);
                        if factors <> false then pi_g :=  (i, j) * pi_g; fi;
                     fi;
                  until j > n or factors <> false;
               fi;
               if factors = false then 
                  permutes := false;
               else 
                  g := factors[1]; i := i + 1; 
               fi;
            od;
            if permutes = true then 
               InfoTensorProduct (k, "-th generator acts as permutation ", 
                       pi_g, " on factors.\n");
               Add (permaction, pi_g); 
               k := k + 1;
            else 
               InfoTensorProduct (k, "-th generator does not permute factors.\n");
            fi; 
         od;
         if permutes = true then 
            SetSymTensorProductFlag (module, true);
            SetSymTensorBasisFlag (module, basis);
            SetSymTensorFactorsFlag (module, tenpow[2]);
            SetSymTensorPermFlag (module, permaction);
            return true;
         fi; 
      fi;
   od;

   return false;

end;

#############################################################################
##
#F  MultipleTensorProductDecomposition (module, d, n, numTries)  . . 
##
## The function uses random methods to try and decompose the module
## `module' as a tensor product of n spaces of dimension d.
## The method is iterative; at each stage of a successful
## decomposition a space W of dimension a power of d is written 
## as tensor product of two such spaces W1, W2 of smaller dimension.  
## This is done using a randomly generated element of an appropriate
## order as input for a call of the function SmashGModule.
## At most numTries random elements are tried at each stage of the 
## of the tensor decomposition.
## 
MultipleTensorProductDecomposition := function (module, d, n, numTries) 

   local h, pair, po, poprimes, r, hh, SS, nn, 
         F, q, N, GLprimes, i, basis, P, factors, 
         module1, module2, dim1, dim2, tenpow1, tenpow2, try, numTries;

   InfoTensorProduct ("#I Trying to decompose as a ", n, "-th tensor power of a ", 
           d, "-dimensional module.\n");

   if d^n <> DimensionFlag (module) then return false; fi;

   F := FieldFlag (module);
   q := Size (F);
   if not IsIrreducible (module) or not IsAbsolutelyIrreducible (module) then
      return false;
   fi;

   GLprimes := [Characteristic (F)];
   for i in [1..d] do 
      if q <> 2 or i <> 1 then  # exclude q^i - 1 = 1
         GLprimes := Concatenation (GLprimes, Set (FactorsInt (q^i - 1))); 
         GLprimes  := Set (GLprimes);
      fi;
   od;

   try := 1;
   while try <= numTries do
      InfoTensorProduct ("#I In MultipleTensorProductDecomposition loop, try ", 
                         try, " of ", numTries, ".\n");
      h := PseudoRandom (module);
      pair := ProjectiveOrderMat (h); po := pair[1];
      if po <> 1 then
         poprimes := Set (FactorsInt (po));
      else 
         poprimes := [];
      fi;

      for r in poprimes do
         if not r in GLprimes then
            InfoTensorProduct ("#I Projective order ", po, " of element incompatible with ", n, "-fold tensor power of V (", d, ", ", q, ").\n");
            InfoTensorProduct ("#I GLprimes = ", GLprimes, "\n");
            InfoTensorProduct ("#I poprimes = ", poprimes, "\n");
            return false;
         fi;
      od;

      for r in poprimes do
         hh := h^(po / r);
         SS := [hh];
         UndoTensorProductFlags (module);
         if SmashGModule (module, SS, "PartialSmash") = true and 
            TensorProductFlag (module) = true then 
            basis := TensorBasisFlag (module);
            module1 := TensorFactorsFlag (module)[1];
            module2 := TensorFactorsFlag (module)[2];
            dim1 := DimensionFlag (module1);
            dim2 := DimensionFlag (module2);
            nn := IntPower (dim1, d);
            if nn <> false then
               if nn <> 1 then
                  tenpow1 := MultipleTensorProductDecomposition 
                             (module1, d, nn, numTries);
                  # We may want to change this: return false if recursive 
                  # call to  MultipleTensorProductDecomposition returns false
                  if tenpow1 = false then return false; fi;
               else 
                  tenpow1 := [IdentityMat (dim1, F), [module1]]; 
               fi;
               nn := n - nn;
               if nn <> 1 and tenpow1 <> false then
                  tenpow2 := MultipleTensorProductDecomposition 
                             (module2, d, nn, numTries);
                  # We may want to change this: return false if recursive call 
                  # to MultipleTensorProductDecomposition returns false
                  if tenpow2 = false then return false; fi;
               else 
                  tenpow2 := [ IdentityMat (dim2, F), [ module2]]; 
               fi;

               if tenpow1 <> false and tenpow2 <> false then 
                  P := KroneckerProduct (tenpow1[1], tenpow2[1]);
                  factors := Concatenation (tenpow1[2], tenpow2[2]);
                  return [P * basis, factors]; 
               fi;
            fi;
         fi;
      od;
      try := try + 1;
   od;

   InfoTensorProduct ("#I Failed to decompose as a ", n, "-th tensor power of a ", 
           d, "-dimensional module.\n");
   return false;

end;

#############################################################################
##
#F  SwapFactors (i, j, dd, n, F)   . . .
## set pi_ij to be a d = dd^n by d matrix over F that swaps 
## the i-th and j-th factors in the tensor product of n
## dd-dimensional spaces. We assume that i < j.
##
## 
SwapFactors := function  ( i, j, dd, n, F) 

   local zero, one, d, k, l, A, B, C, a, b, c, Mb, Mc, r, s, Mr, Ms, pi_ij;

   one := One (F); zero := Zero (F);
   d := dd^n;
   pi_ij := [];
   for k in [1..d] do pi_ij[k] := []; od;
   for k in [1..d] do
      for l in [1..d] do
         pi_ij[k][l] := zero;
      od;
   od;

   Mb := dd^i;
   Mc := dd^j;
   Mr := dd^(i - 1);
   Ms := dd^(j - 1);

   A := dd^(i - 1) - 1;
   B := dd^(j - i - 1) - 1;
   C := dd^(n - j) - 1;

   for a in [0..A] do
      for b in [0..B] do
         for c in [0..C] do
            for r in [0..dd - 1] do
               for s in [0..dd - 1] do
                  k := a + Mr * r + Mb * b + Ms * s + Mc * c + 1;
                  l := a + Mr * s + Mb * b + Ms * r + Mc * c + 1;
                  pi_ij[k][l] := one; pi_ij[l][k] := one;
               od;
            od;
         od;
      od;
   od;    

   return pi_ij;

end;
