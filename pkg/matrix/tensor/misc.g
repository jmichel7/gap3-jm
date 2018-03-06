#############################################################################
##
#A  Matrix package                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
##
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
## 
#############################################################################

ProjectiveOrder := function (g)

   return ProjectiveOrderMat (g)[1];

end;

LogQ := function (F, x)

   if x = Zero (F) then return Size (F) - 1; fi;
   return LogFFE (x);

end;

#remove element with index position index from L

Remove := function (L, index)
   if index > Length (L) then return; fi;
   RemoveSet (L, L[index]);

end;

HashSpace := function (I, M)

   local F, E, L, p, B, d, k, e, x, lenx, i, j, h, y;

   F := BaseRing (I);
   E := Elements (F);
   L := List (E, f -> LogQ (F, f));

   p := Characteristic (F);

   B := Base (I);
   d := Length (B);
   k := Int (3 * d / 4) + 1;

   x := [];
   for i in [1..Length (B)] do
      e := B[i];
      Append (x, List ([k + 1..d], j -> L[ Position (E, e[j]) ]));
   od;

   h := 0; lenx := Length (x);
   for i in [lenx, lenx - 1..1] do
      h := (h + x[i]) * p;
   od;

   h := h mod M + 1;

   return h;

end;

# compute a hash value for a collection of objects S,
# and return result modulo prime M 

HashSet := function (S, M)

   local h, x;

   h := Sum (List (S, x -> HashSpace (x, M)));
   return h mod M + 1;

end;

# does h commute with all elements of X? 

ElementCommutes := function (h, X)

   local id;

   id := h^0;
   return ForAll (X, x -> Comm (x, h) = id);

end;

# first non-zero entry of B 

FirstNonZeroEntry := function (B)

   local E, zero, index;

   # this may be expensive if B is large 
   E := Flat (B);
   zero := 0 * B[1][1]^0;

   if Set (E) = [zero] then return [false, false]; fi; 
   index := 0;
   repeat index := index + 1; until E[index] <> zero;

   return [E[index], index];

end;

# is f a n-th power of some polynomial? 

IsPowerOfPolynomial := function (f, n)

   local factors, mults, i;

   factors := Collected (Factors (f));

   mults := List (factors, i -> i[2]);

   return [ForAll (mults, i -> i mod n = 0), factors];

end;

DeletePair := function (L, d, r)
   return Difference (L, [r, d/r]);
end;

Swap := function (x, i, j)

   local temp;

   temp := x[j];
   x[j] := x[i];
   x[i] := temp;
   return x;

end;

# find space of V centralised by g mod W 

CentralisedMod := function (G, g, W) 

   local bs, e, i, V, F, d, A, Q, x, K;

   F := BaseRing (G);
   d := DimensionFlag (G);

   V := VectorSpace (IdentityMat (d, F), F);
   A := g - g^0;

   bs := BaseSteinitz(V, W);
   bs := Concatenation (bs.subspace, bs.factorspace);
   d := Dimension (V);
   e := Dimension (W);
   bs := bs^-1;
   bs := bs{[1..d]}{[e+1..d]};

   x := A * bs;
   K := NullSpace (x, F);

   return K;

end;

NullSpace := function (g, F)

   local N, zero, d, x;

   N := NullspaceMat (g);

   if N = [] then 
      zero := 0 * g[1][1]^0; 
      d := Length (g);
      N := [List ([1..d], x -> x * zero)];
   fi;

   return VectorSpace (N, F);

end;

# find space centralised by g 

CentralisedSpace := function (g, F)

   local N;
   N := NullSpace (g - g^0, F);
   return N;

end;

# find space centralised by collection of elements, S 

CentralisedSpaceSet := function (F, S)

   local g;
   return Intersection (List (S, g -> CentralisedSpace (F, g)));

end;

# find space centralised by collection of elements, S, mod space W 

CentralisedSpaceSetMod := function (G, S, W)

   local g;
   return Intersection (List (S, g -> CentralisedMod (G, g, W)));

end;

ComplementSpace := function (G, g)

   local N, d, F, x, f, facs, factors, i, pos, rem, h, A, a;

   F := BaseRing (G);
   x := Indeterminate (F);
   d := DimensionFlag (G);
   f := MinimalPolynomial (g);

   facs := Collected (Factors (f));
   factors := List (facs, i -> i[1]);

   pos := Position (factors, x - 1);
   if pos = false then 
      rem := f;
   else 
      rem := f / (facs[pos][1]^facs[pos][2]);
   fi;

   h := Value (rem, g);

   N := NullSpace (h, F);

   return N;

end;

ComplementSpaceSet := function (G, S)

   local g;
   return Sum (S, g -> ComplementSpace (G, g));

end;

# if power of an element in Elts is also in Elts, remove it 

EliminateRepetitions := function (Elts, Orders)

   local i, g, o, D, x, h, index;

   if Length (Elts) <= 1 then return; fi;

   i := 0;
   repeat
      i := i + 1;
      g := Elts[i];
      o := Orders [i];
      D := Set (Factors (o));
      for x in D do
         h := g^x;
         index := Position (Elts, h);
         if index <> 0 then
            Remove (Elts, index);
            Remove (Orders, index);
         fi;
      od;
   until i = Length (Elts);

end;
