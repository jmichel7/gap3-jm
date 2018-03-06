#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @(#)$Id: semilinr.g,v 1.1 1997/03/10 13:52:39 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: semilinr.g,v $
#H  Revision 1.1  1997/03/10 13:52:39  gap
#H  VERSION 1.0
#H
#H  Revision 1.2  1997/01/05 10:49:35  fceller
#H  added Eamonn's new version to the reprository
#H
#H  Revision 1.1  1996/12/25 09:07:37  fceller
#H  changed long filenames to MS-DOS conform filenames,
#H  the init files are *NOT* yet updated
#H
#H  Revision 1.2  1996/12/22 07:48:15  fceller
#H  new smash version
#H
#H  Revision 1.1  1996/11/28 13:14:56  fceller
#H  added "smash" and "reducible" to the repository
#H
##
############################################################################
##
#F  InfoSemiLinear (...)  . . . . . . . . . . . for debugging assistance
##
##
if not IsBound(InfoSemiLinear)  then InfoSemiLinear := Ignore;  fi;

#############################################################################
##
#F  SemiLinearDecomposition( module, S, C, e )  . . 
## applied when S acts irreducibly.
## module is the module of the group G, and <S> is a subgroup of G.
## C centralises the action of <S> on the underlying vector space and
## acts as multiplication by a scalar x
## in GF(q^e) for some embedding of a d/e-dimensional vector
## space over GF(q^e) in the d-dimensional space, where x is a field generator
## of GF(q^e). Thus, provided C centralises the action of <S>^G, the normal
## closure of <S>, <S>^G embeds in GL(d/e,q^e), and  G 
## embeds in GammaL(d/e, q^e). We test for that by trying to 
## construct a map from G to Aut(GF(q^e).    
## We check to see if G can be embedded in GammaL(d/e,q ^e), using the function
## PowerMaps. If PowerMaps returns false, it must be because C doesn't 
## centralise the action of <S>^G, so there is some conjugate of an 
## element of S which is not centralised by C. 
## We find such a conjugate and return false.
## If PowerMaps returns true we return true.
## If true is returned then the SemiLinearFlag is set to true, the 
## DegreeFieldExtFlag is set to e, the LinearPartFlag is set to S, and the
## FrobeniusAutomorphismsFlag is set to the sequence of integers i such that
## multiplication of C by the generator g corresponds to the action of the 
## field automorphism x -> x^q^i(g) on the corresponding element of GF(q^e). 
## 
SemiLinearDecomposition := function (module, S, C, e) 

   local G, matrices, g, powermaps;

   InfoSemiLinear ("#I Looking for powermaps.\n");
   powermaps := PowerMaps (module, C, e);

   if powermaps <> false then
      SetSemiLinearFlag (module, true);
      SetDegreeFieldExtFlag (module, e);
      SetLinearPartFlag (module, S);
      SetCentMatFlag (module, C);
      SetFrobeniusAutomorphismsFlag (module, powermaps);
      return true;
   else
      # Enlarge S by a random conjugate (of an element of S) that doesn't
      # commute with C
      matrices := GeneratorsFlag (module);
      G := Group (matrices, matrices[1]^0);
      repeat 
         g := RandomConjugate (G, S);
      until g * C <> C * g;
      Add (S, g);
   fi;

end;

#############################################################################
##
#F  PowerMaps (module, C, e)  . . part of test to see if 
##  GL (d/e, q^e) < G <= GammaL (d/e, q^e)
## module is d-dimensional, acted on by G, and C is a d x d matrix, which
## acts as multiplication by a scalar x  (a field generator of GF (q^e))
## for some embedding of a d/e-dimensional vector
## space over GF (q^e) in the d-dimensional space.
## G acts as a semilinear group of automorphisms on the d/e-dimensional
## space if and only if, for each generator g of G, there is an integer 
## i = i (g) such that Cg = gC^{q^i}, i.e. g corresponds to the
## field automorphism x -> x^(q^i). 
## Then we have a map from G to the  (cyclic) group Aut (GF (q^e), 
## and C centralises
## the the action of the kernel of this map, which thus lies in GL (d, q^e)
## We test this by first, if possible, finding such i=i (g) 
## such that wCg = wgC^(q^i) for a single vector w of the d-dimensional space
## (in fact the first vector of the standard basis) and then checking that 
## vCg = vgC^ (q^i) for all other vector v in the basis.
## This function returns a list, powermaps, consisting of the integers
## found, or false if no such integers can be found.
## 
PowerMaps := function  (module, C, e) 
   local zero, one, powermaps, g, matrices, found, 
         v, F, q, dim, L, M, N, R, s, i;

   matrices := GeneratorsFlag (module);
   dim := DimensionFlag (module);
   F := FieldFlag (module);
   q := Size (F);
   zero := Zero (F);
   one := One (F);

   powermaps := [];
   for g in matrices do
      found := false;

      v := zero * C[1]; 
      v[1] := one; 

      L := v * C * g;
      M := v * g;
      N := C;
      s := 0;
      repeat 
         R := M * N;
         if L = R then found := true; else N := N^q; s := s + 1; fi;
      until found = true or s = e; 
      if s = e then
         InfoSemiLinear ("#I No powermap found.\n"); return false;
      fi;
      for i in [2..dim] do
         v := zero * C[1]; 
         v[i] := one; 
         M := v * g;
         L := v * C * g;
         R := M * N; 
         if L <> R then
            InfoSemiLinear ("#I No consistent powermap found.\n"); 
            return false;
         fi;
      od;
      Add (powermaps, s);
   od;
   return powermaps;
end;

#############################################################################
##
#F  IsSemiLinear (G)  . . . . . . .  decide if G or GModule is semi-linear 
## 
##  
IsSemiLinear := function (G)

   local module, Result, gens, S, AllScalar, i; 

   if IsMatGroup (G) = false and IsGModule (G) = false then
      return Error ("Argument must be a matrix group or GModule\n");
   fi;

   if IsMatGroup (G) then
      module := GModule (G);
   else
       module := G;
   fi;

   InfoSemiLinear ("#I Input G-module has dimension ", 
        DimensionFlag (module), " over ", FieldFlag (module), "\n"); 

   if DimensionFlag (module) = 1 then return [false, module]; fi;

   #is the module irreducible?
   if IsIrreducible (module) = false then 
      Print ("#I GModule is not irreducible\n");
      return [false, module];
   fi;

   #is the module absolutely irreducible?
   if IsAbsolutelyIrreducible (module) = false then 
      Print ("#I GModule is not absolutely irreducible\n");
      return [true, module];
   fi;

   SemiLinearTest (module);

   #we may discover that module acts imprimitively or semilinearly;
   #otherwise we can conclude that it is not semilinear

   #have we found that it is semilinear?
   if SemiLinearFlag (module) = true then return [true, module]; fi;

   #if we have not found that it's imprimitive, 
   # then we now know that it's not semilinear

   if ImprimitiveFlag (module) <> true then      
      SetSemiLinearFlag (module, false);
      return [false, module]; 
   else 
      #otherwise we have found out that it is imprimitive 
      #and failed to decide semilinearity 
      return ["unknown", module];
   fi;

end; #IsSemiLinear
