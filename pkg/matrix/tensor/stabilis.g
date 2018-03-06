#############################################################################
##
#A  Matrix package                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
##
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
## 
############################################################################
##
#F  InfoTensor1  (...)  . . . . . . . . . . . . . .  for debugging assistance
##
##
if not IsBound (InfoTensor1)  then InfoTensor1 := Ignore;  fi;
#############################################################################

# return GAP integer seed  
GetSeed := function ()
   return [R_N, Copy (R_X)];
end;

# set GAP integer seed  
SetSeed := function (IntSeed)
   R_N := IntSeed[1];
   R_X := Copy (IntSeed[2]);
end;

# generate those elements whose indices are in IndexList

GenerateElements := function (G, IndexList)
   local x, i, Good;

   Good := [];
   for i in [1..Maximum (IndexList)] do
      x := PseudoRandom (G);
      if i in IndexList then Add (Good, x); fi;
   od;

   return Good;
end;

# generate elements of set stabiliser of collection of spaces
# do this by generating NmrImages images of spaces under
# random elements of G

StabiliserOfSet := function (G, Spaces, NmrImages)
   local W, MatrixSeed, IntSeed, M, Collisions, i, g, I, Images,
         index, Stab, Good, L, IndexList, EltNmr, x, xinv, Im, j, y, tmp;

   if Length (Spaces) = 0 then return G; fi;
   W := Spaces[1];

   if Dimension (W) = 0 then return G; fi;

   MatrixSeed := Copy (RandomSeedFlag (G));
   IntSeed := GetSeed ();

   Images := [];

   # we're getting too many coincidences
   M := NextPrimeInt (QuoInt (8 * NmrImages, 3));
   Collisions := [];

   for i in [1..NmrImages] do
      g := PseudoRandom (G);
      I := Set (List (Spaces, W -> W^g));
      index := HashSet (I, M);
      if IsBound (Collisions[index]) = false then
         Collisions[index] := [];
      fi;
      Add (Collisions[index], i);
   od;

   Stab := Set ([]);

   # which random elements do we consider ?

   IndexList := Concatenation (
                Filtered (Collisions, i -> IsBound (i) and Length (i) > 1));
   Sort (IndexList);

   # get these elements
   SetSeed (IntSeed);
   SetRandomSeedFlag (G, MatrixSeed);
   Good := GenerateElements (G, IndexList);
   InfoTensor1 ("#I We consider ", Length (Good), " elements \n");

   # now test them
   for i in [1..Length (Collisions)] do
      if IsBound (Collisions[i]) then
         L := Collisions[i];
         if Length (L) <> 1 then

            EltNmr := Position (IndexList, L[1]);
            x := Good[EltNmr];
            xinv := x^-1;
            Im := Set (List (Spaces, W -> W^x));
            for j in [2..Length (L)] do
               EltNmr := Position (IndexList, L[j]);
               y := Good[EltNmr];
               Images := Set (List (Spaces, W -> W^y));
               if Im = Images then
                  tmp := y * xinv;
                  if not tmp in Stab then
                     Add (Stab, tmp);
                  fi;
               fi;
            od;
         fi;
      fi;
   od;

   return Group (Stab, IdentityFlag (G)); 

end;

# O is either a single space or set of spaces; compute some
# elements of its stabiliser in G; the Boolean returned
# indicates whether the whole of G stabilises O

SpaceStabiliser := function (G, O, Nmr, NmrTries)
   local id, Gens, Spaces, N, fixed, j, W, MinNmrElements, i, K, H;

   if IsVectorSpace (O) then 
      Spaces := Set ([O]);
   else
      Spaces := O;
   fi;

   Gens := GeneratorsFlag (G);
   N := Length (Gens);
   fixed := Filtered (Gens, j -> Set (List (Spaces, W -> W^j)) = Spaces);
   if Length (fixed) = N then return [G, true]; fi;

   id := Identity (G);
   H := Group (fixed, id);

   # we want to find at least this many elements in one try
   MinNmrElements := 2;

   i := 1;
   repeat
      K := StabiliserOfSet (G, Spaces, NmrTries * i);
      InfoTensor1 ("#I Ngens (K) ", Length (Generators (K)), "\n");
      H := Group (Concatenation (GeneratorsFlag (H), GeneratorsFlag (K)), id);
      i := i + 1;
   until Length (GeneratorsFlag (K)) > MinNmrElements or i > Nmr;

   InfoTensor1 ("#I We found ", Length (GeneratorsFlag (H)), 
           " elements of the set stabiliser of the spaces\n");

   return [H, false];

end;
