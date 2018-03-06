#############################################################################
##
#A  Matrix package                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#############################################################################
# 
# given irreducible matrix group G and collection of subspaces of
# associated vector space V, write V as direct sum of subspaces of V

DirectSumSpaces := function (G, Spaces)

   local S, d, k, n, SumSpaces, bound, g, I, degS, sum, D, F, V, i;

   if Length (Spaces) = 0 then return []; fi;

   S := Spaces[1];
   d := Dimension (S);

   # force all of Spaces to have parent V
   F := BaseRing (G);
   V := VectorSpace (Identity (G), F);
   for i in [1..Length (Spaces)] do
      if Parent (Spaces[i]) <> V then 
         Spaces[i] := Subspace (V, Base (Spaces[i]));
      fi;
   od;

   degS := Length (GeneratorsFlag (S)[1]);
   if RemInt (degS, d) <> 0 then
      k := QuoInt (degS, d) + 1;
   else
      k := QuoInt (degS, d);
   fi;

   while Length (Spaces) < k do

      n := Length (Spaces);
      SumSpaces := Sum (Spaces);

      bound := n * d;

      repeat
         g := PseudoRandom (G);
         I := S^g;
         I := Subspace (V, Base (I));
         sum := SumSpaces + I;
         D := Dimension (sum);
      until D > bound;

      if (D <> (n + 1) * d) then
         Spaces := [Intersection (SumSpaces, I)];
         Spaces := DirectSumSpaces (G, Spaces);
         return Spaces;
      fi;

      Add (Spaces, I);

   od;

   return Spaces;

end;
