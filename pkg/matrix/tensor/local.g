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
if not IsBound (InfoTensor)  then InfoTensor := Ignore;  fi;
#############################################################################
#we need the Lattice command from meataxe package
RequirePackage ("meataxe");

# we believe that the best local subgroup is that 
# whose G-Module has shortest composition length 

BestCompLength := function (Subs)

   local Len, CF, i, Len, M, m, index;

   M := List (Subs, x -> GModule (x));
   CF := List (M, x -> CompositionFactors (x));
   Len := CompositionSeriesLength (CF);

   InfoTensor1 ("#I Composition lengths are  ", Len, "\n");
   m := Minimum (Len);
   index := Position (Len, m);
   # InfoTensor1 ("#I Shortest composition length is ", m, "\n");

   return [m, index]; 

end;

# g is a power of parent; find a p-local subgroup of G which
# contains the normaliser of g  

LocalSubgroup := function (G, parent, g, Nmr, NmrTries)

   local p, H;

   p := Characteristic (BaseRing (G));

   if ProjectiveOrder (g) = p then 
      H := pLocalSubgroups (G, parent, g, Nmr, NmrTries);
   else
      H := pDashSubgroup (G, parent, g, Nmr, NmrTries);
   fi;

   return H;

end;

SubmoduleLattice := function (M, Max)

   local F, A, L, mgens;

   F := FieldFlag (M);
   mgens := List (GeneratorsFlag (M), x -> MeatAxeMat (x));

   A := NaturalModule (UnitalAlgebra (F, mgens));
   L := Lattice (A);
   InfoTensor1 ("#I The number of submodules is ", Length (L.dimensions), "\n"); 
   return L;

end;

SubmoduleLatticeAbort := function (M, Max)

   InfoTensor1 ("#I Submodule lattice abort -- This is a dummy function\n");
   return [false, false];

end;

# does the submodule of H have a short lattice? 

HasShortLattice := function (H)

   local MaxLen, M, MaxNmrSubmodules, R, AllFound;

   # maximum composition length 
   MaxLen := 4;

   M := GModule (H);
   #cl := CompositionSeriesLength (CompositionFactors (M));
   #if cl > MaxLen then 
   #   InfoTensor1 ("#I Lattice has composition length ", cl, "\n");
   #   return false; 
   #fi;

   # maximal number of submodules to compute 
   MaxNmrSubmodules := 500;

   R := SubmoduleLatticeAbort (M, MaxNmrSubmodules);
   AllFound := R[1];

   return AllFound;

end; 

# compute lattice of H-submodules; if ComputePartial is true, 
# then return partial lattice 

ComputeLattice := function (H, MaxNmrSubmodules, ComputePartial)

   local MaxLen, M, MaxNmrSubmodules, R, AllFound, L;

   M := GModule (H);

   if ComputePartial then 
      L  := SubmoduleLattice (M, MaxNmrSubmodules);
      AllFound := true;
   else 
      R := SubmoduleLatticeAbort (M, MaxNmrSubmodules);
      AllFound := R[1]; L := R[2];
      if not AllFound then L := []; fi;
   fi; 

   if AllFound then 
      InfoTensor1 ("#I Length of submodule lattice is ", Length (L.dimensions), "\n");
   else 
      InfoTensor1 ("#I Lattice construction incomplete; it has length at least ", 
              MaxNmrSubmodules, "\n");
   fi;

   return [L, AllFound];

end;

# construct local subgroups for G which may act reducibly on 
# some of the factors whose dimensions are recorded in Dims 

LocalSubgroups := function (G, Nmr, NmrTries, Dims)

   local ShortCompLength, NmrpLocal, TotalLocal, NmrSearches, Subs, p, x,
         R, Settled, NewSubs, i, j, r, parent, g, o, Orders, D, pDashSubs, Elts;

   # parameters 
   ShortCompLength := 4;
   NmrpLocal := 5;
   TotalLocal := 10;
   NmrSearches := 25;

   p := Characteristic (BaseRing (G));

   InfoTensor1 ("\n#I First look for ", p, "-local subgroups", "\n");

   R := pLocalElements (G, p, NmrSearches);
   Elts := R[1]; Orders := R[2];

   # distinct orders of parents 
   InfoTensor1 ("#I Orders of parents of ",  p, "-local elements are ", Orders, "\n");

   Subs := [];

   i := 0;
   while (i < Length (Orders) and Length (Subs) < NmrpLocal) do 
      i := i + 1;
      o := Orders[i];
      parent := Elts[i];
      g := parent^(o / p);
      R := pLocalSubgroups (G, parent, g, Nmr, NmrTries);
      NewSubs := R[1]; Settled := R[2];

      # have we found a p-local subgroup with a small submodule lattice? 
      if Settled then 
         r := rec (H := NewSubs[1], pLocal := true, 
                   pLocalElement := g, ReducibleDims := Dims);
         Subs := [r];
         return Subs; 
      fi;

      for j in [1..Length (NewSubs)] do 
         r := rec (H := NewSubs[j], pLocal := true, 
                   pLocalElement := g, ReducibleDims := Dims);
         Add (Subs, r);
      od;

   od;

   R := pDashLocalElements (G, NmrSearches);
   Elts := R[1]; Orders := R[2];

   InfoTensor1 ("#I Orders of parents of ", p, "-dash elements are ", Orders, "\n");

   # distinct primes obtainable from Orders 
   D := Set (Flat (List (Orders, x -> DistinctPrimes (x))));
   InfoTensor1 ("#I Distinct primes are ", D, "\n");

   InfoTensor1 ("#W Currently p-dash code is incomplete\n");
   return Subs;

   # now compute local subgroups for different primes 
   i := 0;
   while (i < Length (Orders) and Length (Subs) < TotalLocal) do 
      i := i + 1; 
      o := Orders[i];
      parent := Elts[i];
      Settled := false;
      pDashSubs := pDashLocals (G, parent, o, Nmr, NmrTries, Dims, Settled);
      Append (Subs, pDashSubs);
      if (Settled) then return Subs; fi;
   od;

   return Subs;

end;

# investigate each submodule in the lattice and decide if it 
# is a flat in a u-projective geometry for some u in DimU 

InvestigateLattice := function (G, L, DimU)

   local F, d, V, j, B, R, S, TestDim, Spaces, u, Status, CB, x;

   F := BaseRing (G);
   d := DimensionFlag (G);
   V := VectorSpace (IdentityMat (d, F), F);

   #is 1 the trivial and Length (L) the whole?
   for j in [2..Length (L.dimensions) - 1] do 

      B := GeneratorsSubmodule (L, j);
      S := Subspace (V, B.operations.GapObject (B));

      TestDim := Set (Filtered (DimU, x -> Dimension (S) mod x = 0));
      if Length (TestDim) <> 0 then 

         Spaces := [S];
         InfoTensor1 ("#I Processing space ", j, " of dimension ", Dimension (S), "\n");
         InfoTensor1 ("#I First call DirectSumSpaces\n");
         Spaces := DirectSumSpaces (G, Spaces);
         for u in TestDim do 
            InfoTensor1 ("#I Investigating geometry for u = ", u, "\n");
            R := FindPoint (G, Spaces, u);
            Status := R[1]; CB := R[2];
            if Status = true then 
               if IsList (CB) then 
                  InfoTensor ("#I Found geometry for u = ", CB[2], "\n"); 
                  return [Status, CB]; 
               else 
                  InfoTensor ("#I Probably found tensor decomposition over larger field", "\n");
               fi;
            fi;
         od;
      fi;
   od;

   Status := false; CB := "undefined";

   return [Status, CB]; 

end;


CompositionSeriesLength := function (CF)

   local i, j, L, Len;

   Len := [];
   for i in [1..Length (CF)] do
      L := [];
      for j in [1..Length (CF[i])] do
         Add (L, CF[i][j][2]);
      od;
      Len[i] := Sum (L);
   od;

   return Len;

end;

# run local subgroup test; Dims is possible dimensions of factorisations 

LocalTest := function (G, Dims, Nmr, NmrTries) 

   local Result, M, Len, Subs, Trial, d, MaxNmrSubmodules, NmrTrials,
         CF, RedH, Status, DimU, Search, i, H, m, p, R, AllFound, L, x;

   Result := false;
   Subs := LocalSubgroups (G, Nmr, NmrTries, Dims);

   #o := List (Subs, x ->  x.H);
   #InfoTensor1 ("#I Orders for these subgroups is ", o, "\n");

   # we use composition lengths to sort local subgroups 

   M := List (Subs, x -> GModule (x.H));
   CF := List (M, x -> CompositionFactors (x));
   Len := CompositionSeriesLength (CF);
   InfoTensor1 ("#I Composition lengths are ", Len, "\n");

   d := DimensionFlag (G);
   MaxNmrSubmodules := [1000, 10000];
   NmrTrials := 2; 

   SortParallel (Len, Subs);

   Trial := 0;
   repeat 
      Trial := Trial + 1;
      for i in [1..Length (Subs)] do
         InfoTensor1 ("#I Possible dimensions of U are now ", Dims, "\n");
         H := Subs[i].H;
         m := Len[i];
         InfoTensor1 ("#I Shortest composition length for local subgroup H is ", m, "\n");
         p := ProjectiveOrder (Subs[i].pLocalElement);
         InfoTensor1 ("#I The subgroup is ", p, "-local", "\n");

         if Subs[i].pLocal then 
            # here, H acts reducibly in all dimensions but there may be a 
            # more general version later 
            RedH := Subs[i].ReducibleDims;
            InfoTensor1 ("#I H is guaranteed to act reducibly in dimension ", RedH, "\n");
         fi;

         #R := ComputeLattice (H, MaxNmrSubmodules[Trial], Trial = NmrTrials);
         R := ComputeLattice (H, MaxNmrSubmodules[Trial], true);
         L := R[1]; AllFound := R[2];

         # we may change this later 
         if (AllFound = true or Trial = NmrTrials) then 

            # consider projective geometries for relevant dimensions 

            Search := Dims;

            InfoTensor1 ("#I Search for geometry in following dimensions ", Search, "\n");
            R := InvestigateLattice (G, L, Search);
            Status := R[1]; Result := R[2];
            if Status = true then return Result; fi;

            # there may be a more general version where pLocal is not true 
            if AllFound and Subs[i].pLocal then 
               DimU := Set ([RedH, x -> d / x]);
               SubtractSet (Dims, DimU);
            fi;

            if Length (Dims) = 0 then Result := false; return Result; fi;

         fi; #if (AllFound = true or Trial = NmrTrials) 

      od;
   until Trial = NmrTrials; 

   if Length (Dims) = 0 then 
      Result := false; 
   else 
      Result := "unknown"; 
   fi;

   return Result;

end;
