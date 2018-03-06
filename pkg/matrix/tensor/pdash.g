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

JordanBlocks := Print;

# generate some elements having orders which
# are not powers of the characteristic 

pDashLocalElements := function (G, NmrTries)

   local Orders, Elts, p, x, i, o;

   Orders := [];
   Elts := [];

   p := Characteristic (BaseRing (G));

   for i in [1..NmrTries] do 
      x := PseudoRandom (G);
      o := ProjectiveOrder (x);
      if not IsAPower (o, p) and not (o in Orders) then 
         Add (Elts, x);
         Add (Orders, o);
      fi;
   od;

   SortParallel (Orders, Elts);
   Orders := Reversed (Orders);
   Elts := Reversed (Elts);
   EliminateRepetitions (Elts, Orders);

   return [Elts, Orders];

end;

# given g of projective order p <> characteristic of F, 
# find a subgroup H of G which normalises a p-subgroup of G 

pDash := function (G, g, Nmr, NmrTries) 

   local F, p, o, d, V, Spaces, Degrees, R, X, x, LargeField,
         T, H, fixed, i, NewSum;

   F := BaseRing (G);
   p := Characteristic (F);

   o := ProjectiveOrder (g);

   if not IsPrime (o) or Gcd (p, o) <> 1 then 
      InfoTensor1 ("#I Element has order ", o, "\n");
      return [false, false];
   fi;

   if IsScalar (g) then 
      InfoTensor1 ("#I Element is scalar\n");
      return [G, false];
   fi;

   d := Dimension (G);

   #V := VectorSpace (F, d);

   R := JordanBlocks (g);
   Spaces := R[1]; Degrees := R[2];
   InfoTensor1 ("#I The number of Jordan blocks is ", Length (Spaces), "\n");

   X := Set (Degrees);
   InfoTensor1 ("#I Distinct degrees in Jordan form are ", X, "\n");

   # are all the eigen values present in the field? 
   LargeField := X <> [1];

   H := G;

   if Length (X) > 1 then 
      for x in X do
         NewSum := Sum (Filtered (Spaces, i -> Degrees[i] = x));
         R := SpaceStabiliser (H, NewSum, Nmr, NmrTries);
         H := R[1]; fixed := R[2];

         #if we get a subgroup of H as new stabiliser, 
         # then we must initialise the seed for H  ??? probably not 
         #necessary to worry in GAP context 

      od;
   fi; 

   for x in X do
      T := Set (Filtered (Spaces, i -> Degrees[i] = x));
      R := SpaceStabiliser (H, T, Nmr, NmrTries);
      H := R[1]; fixed := R[2];
      # if we got a proper subgroup of H as the stabiliser of T, 
      #  then we must initialise the seed for H 
   od;

   H := Group (UniteSet (GeneratorsFlag (H), g), g^0);

   return [H, LargeField];

end;

# write H defined over GF (q) over larger field GF(q^n) 

EmbedLargerField := function (H, q, n) 

   local Large;
   Large := Copy (H);
   SetFieldFlag (Large, GF(q^n));

   return Large;

   #Large := GL (DimensionFlag (H), q^n);
   #Large := Large.operations.Subgroup (Large, GeneratorsFlag (H));
   #return Group (GeneratorsFlag (H), H.identity);

end;

# the larger field of H is an extension of a base field and 
# we write H over the base field 

RestrictSmallerField := function (H, q)

   local Small;
   Small := Copy (H);
   SetFieldFlag (Small, GF(q));

   return Small;

   #Small := GL (DimensionFlag (H), q);
   #Small := Small.operations.Subgroup (Small, GeneratorsFlag (H));
   #return Group (GeneratorsFlag (H), H.identity);

end;

# compute subgroup H which is p-local; in particular <g>^H is an 
# elementary abelian p-group; g has order p^k, projective order p 

pDashSubgroup := function (G, parent, gg, Nmr, NmrTries) 

   local R, g, H, LargeField, q, n; 

   g := Copy (gg);

   R := pDash (G, g, Nmr, NmrTries); 
   H := R[1]; LargeField := R[2];

   H := Group (UniteSet (GeneratorsFlag (H), parent), parent^0);

   if LargeField then 

      q := Size (BaseRing (G));

      # degree of extension 
      n := OrderMod (q, OrderMat (g));

      H := EmbedLargerField (H, q, n);

      R := pDash (H, g, Nmr, NmrTries); 
      H := R[1]; g := R[2];

      if LargeField then return Error ("** ERROR in pDashSubgroup **\n"); fi; 

      H := RestrictSmallerField (H, q);

   fi;

   return H;

end;

# parent is element of order o; generate g, an element of prime order 
# p from parent; construct p-local subgroup for this prime and 
# store it in Subs 

pDashLocals := function (G, parent, o, Nmr, NmrTries, List, Settled)

   local Subs, p, D, d, Settled, H, r, g, x;

   D := DistinctPrimes (o); 
   d := DimensionFlag (G);

   Settled := false;
   Subs := [];

   # don't consider p-local subgroups here 
   p := Characteristic (BaseRing (G));
   D := Filtered (D, x -> x <> p);

   for p in D do

      g := parent^(o / p);
      if not IsScalar (g) then 
         InfoTensor1 ("#I Consider element of projective order ", p, "\n");

         # compute local subgroup containing g 
         InfoTensor1 ("#I Look for local subgroup for prime ", p, "\n");
         H := pDashSubgroup (G, parent, g, Nmr, NmrTries);
         r := rec (H := H, pLocal := false, pLocalElement := g);
         Add (Subs, r);
      fi;

   od;

   return Subs;

end;
