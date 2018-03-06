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
# generate some elements having powers which have order RequiredOrder;
# in practice, we use this where RequiredOrder is the characteristic;
# we sort in decreasing order and seek to remove powers of elements
# from the list 

pLocalElements := function (G, RequiredOrder, NmrTries)

   local e, M, N, o, i, rem, x, p, Orders, Elts;

   p := Characteristic (BaseRing (G));

   Orders := [];
   Elts := [];

   M := RootInt (NmrTries);
   for i in [1..M] do
      N := 1;
      repeat
         x := PseudoRandom (G);
         if RequiredOrder = p then 
            o := OrderMat (x);
         else 
            o := ProjectiveOrder (x);
         fi;
         rem := (o mod RequiredOrder = 0);
         N := N + 1;
      until rem or (N = M);
      if rem and not (o in Orders) then
         Add (Orders, o);
         Add (Elts, x);
      fi;
   od;

   SortParallel (Orders, Elts);
   Orders := Reversed (Orders);
   Elts := Reversed (Elts);
   EliminateRepetitions (Elts, Orders);

   return [Elts, Orders];

end;

ExtractBlock := function (A, row, col, rowdim, coldim)

   local B, i, j;

   B := [];
   for i in [row..row + rowdim - 1] do 
      Add (B, []);
      for j in [col..col + coldim - 1] do 
         Add (B[Length (B)], A[i][j]);
      od;
   od;

   return B;

end;

# find action of elements of Gens on W and on V / W, 
# given change of basis matrix, C 

SubQuotAction := function (Gens, W, C)

   local d, A, S, Q, dimS, dimW, dimQ, x, i;

   d := Length (Gens[1]);

   A := [];
   for x in Gens do 
      Add (A, C * x * C^-1);
   od;

   S := [];
   Q := [];
   dimW := DimensionFlag (W);
   dimQ := d - dimW;

   for i in [1..Length (A)] do
      S[i] := ExtractBlock (A[i], 1, 1, dimW, dimW);
      Q[i] := ExtractBlock (A[i], dimW + 1, dimW + 1, dimQ, dimQ);
   od;

   return [S, Q];

end ;

# given G is a matrix group of characteristic p
# X is a collection of elements of order dividing p;
#
# Is <X>^G a p-group or trivial?
#
# We need to find a chain of G-submodules
#
# 0 = V_0 <= V_1 <= ... <= V_k = V
#
# of the natural module such that x centralises
# V_i / V_(i - 1) for all x in X and all i 

IspNormal := function (G, X)

   local R, M, Sub, Quo, C, s, q, resultS, resultQ;

   if IsGroup (G) = true then 
      M := GModule (G);
   elif IsGModule (G) then 
      M := G;
   else 
      return Error ("Input must be G-module or group\n");
   fi;

   if Set (X) = [Identity (G)] then return true; fi;

   if IsIrreducible (M) = true then return false; fi;

   R := InducedAction (M, SubbasisFlag (M));
   Sub := R[1]; Quo := R[2]; C := R[4];

   R := SubQuotAction (X, Sub, C);

   s := R[1]; q := R[2];

   resultS := IspNormal (Sub, s);
   resultQ := IspNormal (Quo, q);

   return (resultS and resultQ);

end;

# compute chain of subspaces  
#      0 = V_0 <= V_1 <= ... <= V_k = V
# of the natural module such that x centralises 
# V_i / V_(i - 1) for all x in X and all i 

ComputeSpaces := function (G, X)

   local d, F, W, V, D, Spaces, i;

   d := DimensionFlag (G);
   F := BaseRing (G);
   W := NullSpace (G.1, F); #just to get vector space with basis zero

   Spaces := [];
   i := 0;
   repeat 
      W := CentralisedSpaceSetMod (G, X, W);
      i := i + 1;
      Spaces[i] := W;
   until Dimension (W) = d;

   # Spaces := [Spaces[j] : j in [1..i]];
   D := List (Spaces, x -> Dimension (x));

   return [Spaces, D];

end;

# X is a collection of non-scalar elements of G and generates a 
# p-subgroup of G, W is an <X>-invariant space; return a p-local 
# subgroup H which normalises a p-subgroup containing X 

pLocal := function (G, X, W, Nmr, NmrTries)

   local F, p, U, R, x, H, gens;

   F := BaseRing (G);
   p := Characteristic (F);

   # do all elements have order dividing p? 
   if Length (Filtered (List (X, x -> OrderMat (x)), 
              x -> x <> p and x <> 1)) <> 0 then 
      InfoTensor1 ("#W Supplied set contains element of order different from ", p, "\n");
      return false;
   fi;

   W := NullSpace (G.1, F); #just to get vector space with basis zero
   U := CentralisedSpaceSetMod (G, X, W);

   R := SpaceStabiliser (G, U, Nmr, NmrTries);
   H := R[1];
   gens := GeneratorsFlag (H);
   Append (gens, X);
   H := Group (gens, gens[1]^0);

   x := IspNormal (H, X);
   if IspNormal (H, X) then return H; fi;

   return pLocal (H, X, U, Nmr, NmrTries);

end;

# can we find element of order p in H which commutes with every 
# element of X? we also require that its chain of subspaces  
# differs from those stored in Spaces 

CommutingElement := function (G, H, X, g, Spaces, p, NmrTries)

   local S, D, R, N, h, s, i;

   #how hard do we try to find new element? 
   N := 10;

   s := Group (g, g^0);
   i := 0;
   repeat
      i := i + 1;
      h := ElementOfOrder (H, p, N);
      if IsBool (h) = false then 
         if not (h in s) and ElementCommutes (h, X) then 
            R := ComputeSpaces (H, [h]);
            S := R[1]; D := R[2];
            if not (S in Spaces) then 
               InfoTensor1 ("#I New element found on attempt ", i, "\n");
               InfoTensor1 ("#I Dimensions of associated spaces is ", D, "\n");
               return [true, h, S, D]; 
            fi;
         fi;
      fi;
   until i = N;

   InfoTensor1 ("#I No new element found after ", N, " attempts\n");
   return [false, Identity (H), false, false];

end; 

# construct a list of p-local subgroups where each p-group 
# contains g and may have rank at most Rank 

pLocalSubgroups := function (G, parent, g, Nmr, NmrTries)

   local Subs, Settled, MaxRank, p, X, W, H, Bin, Spaces, 
         Dims, R, found, h, S, D, i, gens;

   X := [g];
   W := [];

   # compute p-local subgroup containing g 
   H := pLocal (G, X, W, Nmr, NmrTries);
   gens := GeneratorsFlag (H);
   Add (gens, parent);
   H := Group (gens, gens[1]^0);
   # InfoTensor1 ("#I Order is ", Size (H), "\n");

   Subs := [H]; 
   Settled := HasShortLattice (H);

   # if H is cyclic, we will not find a new element 
   if Length (GeneratorsFlag (H)) = 1 or Settled then 
      return [Subs, Settled]; 
   fi;

   MaxRank := 3;
   Bin := [X];

   Spaces := []; Dims := [];
   # compute information about g 
   R := ComputeSpaces (G, Bin[1]);
   Spaces[1] := R[1]; Dims[1] := R[2];
   InfoTensor1 ("#I Dimensions of associated spaces is ", Dims, "\n");

   p := Characteristic (BaseRing (G));

   # now look for element h in p-local subgroup H which commutes with 
   # all elements in the bin of H; if the dimensions of its associated 
   # spaces are the same as those of elements in that bin, we expect that 
   # the elements are conjugate; we add h to the bin for H and compute 
   # p-local containing this larger bin; otherwise, we create new bin 
   # with single element h  

   i := 0;
   repeat
      i := i + 1;
      H := Subs[i];
      X := Bin[i];

      if Length (GeneratorsFlag (H)) <> 1 then 

         R := CommutingElement (G, H, X, g, Spaces, p, NmrTries);
         found := R[1]; h := R[2]; S := R[3]; D := R[4];

         # if we find a new element, its construction ensures that the 
         # new element h and the elements of X are not conjugate in H 

         # did we find new element? 
         if found then 

            # if the new element has the same dimensions as the elements 
            # of the bin, then they are probably conjugate in G (we know
            # by construction that these are not conjugate in H); if so,
            # compute local subgroup containing existing bin + new element; 
            # otherwise compute local subgroup for new element 

            if D = Dims[i] then
               X := Concatenation (Bin[i], [h]);
            else 
               X := [h];
            fi;

            InfoTensor1 ("#I i = ", i, " Length (X) is now ", Length (X), "\n");

            Add (Bin, X);
            Add (Spaces, S);
            Add (Dims, D);

            H := pLocal (G, X, W, Nmr, NmrTries);

            # can this subgroup settle the computation? 
            if HasShortLattice (H) then 
               Subs := [H]; Settled := true; 
               return [Subs, Settled]; 
            fi;

            Add (Subs, H);
         fi;

      fi; #H is not cyclic

   until (i = Minimum (MaxRank, Length (Subs))); 

   Settled := false;

   return [Subs, Settled]; 

end;
