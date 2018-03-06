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

# return partitions of set S into subsets

SetPartitions := function (S)
   local x, X, T, Q, P;

   if Length (S) = 1 then
      return Set ([Set ([S])]);
   else
      x := Random (S);
      RemoveSet (S, x);
      P := SetPartitions (S);
      Q := Set ([]);
      for X in P do
         for T in X do
            Q := Union (Q, Set ([Union (Difference (X, Set ([T])), 
                         Set ([Union (T, Set ([x]))]))]));
         od;
         Q := Union (Q, Set ([Union (X, Set ([Set ([x])]))]));
      od;
      return Q;
   fi;

end;

# return the sum of the lcms of the elements of X

Score := function (X)
   local T;

   return Sum (List (X, T -> Lcm (List (T))));

end;

# returns the least d such that GL (d, q) has an element
# of order n and the corresponding factorisation of n;
# Gcd (n, q) = 1

LeastLinearSemiSimple := function (n, q)
   local p, f, orders, x, y, S, scores, least, nparts, i;

   p := Collected (Factors (q))[1][1];

   f := Collected (Factors (n));

   # for each prime-power factor of n, store the least d
   # such that GL (d, q) has an element of this order
   orders := Set (List ([1..Length (f)], i -> OrderMod (q, f[i][1]^f[i][2])));

   # also need to remove divisors
   for x in orders do
      for y in orders do
         if x < y and y mod x = 0 then RemoveSet (orders, x); fi;
      od;
   od;

   if Length (orders) = 0 then return [1, 1]; fi;

   S := SetPartitions (orders);
   S := List (S);

   scores := List (S, x -> Score (x));

   # minimise partition score subject to maximising number of parts

   least := scores[1]; nparts := Length (S[1]);
   for i in [2..Length (scores)] do
      if (scores[i] <= least) and (Length (S[i]) > nparts) then
         least := scores[i];
         nparts := Length (S[i]);
      fi;
   od;

   return [least, nparts];

end;

LeastProjectiveSemiSimple := function (n, q)
   local R, m, u;

   if Gcd (n, q) <> 1 then return false; fi;

   if n = 1 then return 1; fi;

   R := LeastLinearSemiSimple (n, q);
   u := R[1]; m := R[2];

   if m > 1 then return u; fi;

   if QuoInt (q^u - 1, q - 1) mod n = 0 then return u; fi;

   return u + 1;

end;

# find smallest d such that PGL (d, q) can contain an element
# of projective order n

LeastProjective := function (n, q)
   local R, p, f, primes, index, alpha, factor, m;

   # write n = p^alpha * m
   p := Collected (Factors (q))[1][1];

   # find p'-part of projective order n
   f := Collected (Factors (n));
   primes := List ([1..Length (f)], i -> f[i][1]);

   index := Position (primes, p);
   if IsBool (index) then
      alpha := 0;
      factor := 1;
   else
      alpha := f[index][2];
      factor := 1 + p^(alpha - 1);
      f := List (Concatenation ([1..index - 1], [index + 1..Length (f)]), 
                 x -> f[x]);
   fi;

   # p'-part of projective order
   m := Product (List ([1..Length (f)], i -> f[i][1]^f[i][2]));

   if alpha = 0 then
      return LeastProjectiveSemiSimple (m, q);
   elif m = 1 then
      return factor;
   else
      R := LeastLinearSemiSimple (m, q);
      return factor + R[1];
   fi;

end;

# take prime powers of n; compute scores for elements having
# these prime power order; choose prime with largest score
# as best prime

FindBestPrime := function (n, q)
   local D, Score, max, index, p, s;

   D := PrimePowersInt (n);

   Score := List ([1, 3 .. Length (D) - 1], 
                  i -> LeastProjective (D[i] * D[i + 1], q));

   max := Maximum (Score);
   index := Position (Score, max);

   p := D[2 * (index - 1) + 1];
   s := D[2 * index];

   return [p, s];

end;

# is there a co-prime factorisation of n as k * l such
# that Score (k * m) <= DimU and Score (l * m) <= DimW ?

ExistsFactorisation := function (n, m, d, q, DimU, DimW)
   local P, x, u, y;

   # find co-prime factorisations of n
   P := CoPrimeFactorisations (n);

   # is there a valid co-prime factorisation of n?
   # -- that is, one whose components fit into each side

   for x in P do
      InfoTensor1 ("#I Processing order factorisation ", x, "\n");
      u := List ([1..2], i -> LeastProjective (x[i] * m, q));

      y := List ([1..2], i -> x[i] * m);
      InfoTensor1 ("#I Score for ", y, " = ", u, "\n");

      # is DimU >= u[1] and DimW >= u[2]?
      if (u[1] <= DimU) and (u[2] <= DimW) then
         return true;
      fi;
   od;

   return false;

end;

# can an element g of order n rule out possible tensor
# factorisation DimU x DimW of a subgroup of GL (d, q) ?
# TestedPrimes records the prime order of elements obtained
# as powers of g which are not projectivities

PossibleFactorisation := function (G, g, n, d, q, DimU, DimW, D, TestedPrimes)

   local R, m, Result, CB, h, flag, s, p;

   m := 1; 
   flag := ExistsFactorisation (n, m, d, q, DimU, DimW);

   while flag do

      p := FindBestPrime (n, q);
      s := p[2];
      p := p[1];
      h := g^QuoInt (m * n, p);
      InfoTensor ("#I Projective order of possible scalar element is ", 
              ProjectiveOrder (h), "\n");

      if not (h in TestedPrimes) then
         R := IsProjectivity (G, h, D);
         Result := R[1]; CB := R[2];
         if Result = true then return CB; fi;
         if Result = "unknown" then return Result; fi;
         AddSet (TestedPrimes, h);
      fi;

      # we can now conclude that if there is such a tensor decomposition, 
      # then an element of order m acts as a non-scalar on both factors
      n := QuoInt (n, p^s); m := m * p^s;
      InfoTensor1 ("#I n is now ", n, " m is now ", m, "\n");
      flag := ExistsFactorisation (n, m, d, q, DimU, DimW);
   od;  

   return flag;

end;

# generate random elements and try to decide whether an element
# of projective order n rules out any possible tensor factorisation
# of a subgroup GL (d, q)

OrderTest := function (G, N, L, Record)
   local F, q, Result, NmrElts, Tested, TestedPrimes, 
         D, Facs, R, g, n, i, u, w, MinScore, d, Fit;

   F := BaseRing (G);
   q := Size (F);
   d := DimensionFlag (G);

   Result := false;
   NmrElts := 0;
   Tested := [];

   TestedPrimes := Set ([]);
   Facs := L;

   # generate N random elements and compute their scores
   repeat
      g := PseudoRandom (G);
      NmrElts := NmrElts + 1;
      n := ProjectiveOrder (g);

      if not (n in Tested) then
         InfoTensor1 ("\n#I Processing Element #", NmrElts, 
                 " of projective order ", n, "\n");
         # what is the smallest dimension which can contain
         # element of projective order n?
         MinScore := LeastProjective (n, q);
         Add (Record, [g, n]);
         Add (Tested, n);

         # now consider each possible factorisation of d as u x w
         # and decide whether such a factorisation of d is compatible
         # with an element of projective order n

         L := Copy (Facs);
         i := 1; Fit := false;
         while i <= Length (L) and Fit = false do 
            u := L[i][1]; w := L[i][2];
            InfoTensor1 ("\n#I Consider dimension factorisation u = ", u, 
                    ", w = ", w, "\n");
            Fit := (MinScore <= u) and (MinScore <= w);
            if Fit then 
               InfoTensor1 ("#I Element of projective order ", n, 
                       " fits into both factors\n");
            else 
               # the element doesn't fit into both factors;
               # however, there may exist a coprime factorisation;
               # we may also be able to conclude that the element can't
               # act in desired manner by calls to IsProjectivity

               D := Set (Concatenation (Facs));
               Result := PossibleFactorisation (G, g, n, d, q, u, 
                                 w, D, TestedPrimes);

               if Result = false then 
                  InfoTensor ("#I No valid score exists for dimension factorisation ", L[i], "\n");
                  RemoveSet (Facs, L[i]);
                  InfoTensor ("#I Valid dimensions are now ", Facs, "\n");

               elif Result = "unknown" then 
                  InfoTensor1 ("#I Could not settle this element\n");

               elif not IsBool (Result) then
                  InfoTensor ("#I  Found tensor decomposition\n");
                  return Result;

               fi;

            fi; #not Fit

            i := i + 1;
         od;
      fi;
   until (Length (Facs) = 0) or (NmrElts >= N) or (Result = true);

   return [Result, Facs];

end;
