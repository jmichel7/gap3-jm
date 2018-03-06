#############################################################################
##
#A  Matrix package                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
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

IsSystemSolvable := function (A, v)

   return [false, false, false];

end;

NonNegativeSolution := function (A, t)

   InfoTensor1 ("#I NonNegativeSolution -- This is a dummy function\n");
   return false;

end; 

# find tensor product of polynomials f and g 
# ensure that it's embedded in F 

PolynomialTensorProduct := function(F, f, g)

   local PR;

   PR := PolynomialRing (F);

   return EmbeddedPolynomial (PR, CharacteristicPolynomial (
                  KroneckerProduct (CompanionMatrix (PR, f), CompanionMatrix (PR, g))));

end;

# take some power of g, an element of (projective) order n 
# to obtain an element of (projective) order at most Limit 

PowerOfSmallOrder := function (g, n, Limit)

   local f, powers, power, newg;

   f := Collected (Factors (n));
   powers := List (f, x -> x[1]^x[2]);
   powers := Filtered (powers, x -> x <= Limit);

   if Length (powers) > 0 then 
      power := Random (powers);
      newg := g^(n / power);
   else
      newg := g^(n / f[1][1]);
   fi;

   return newg;

end;

# use the inherent symmetry of a left-hand factor to write 
# down permutations to reduce the number of possible solutions 

ApplySymmetry := function (F, n, PolyBasis) 

   local fixed, j, P, lp, q, factor, omega, lambda, f, image, o, p, Perms, i;

   q := Size (F);
   factor := Gcd (q - 1, n);
   InfoTensor1 ("#I factor is ", factor, "\n");

   if factor = 1 then return []; fi;

   omega := Root (F);
   lambda := omega^((q - 1) / factor);

   f := Indeterminate (F) - lambda;
   image := [];
   for i in [1..Length (PolyBasis)] do
      image[i] := Position (PolyBasis, 
                          PolynomialTensorProduct (F, PolyBasis[i], f));
   od;

   InfoTensor1 ("#I image = ", image, "\n");

   lp := Length (PolyBasis) + 1;

   p := PermList (image);
   o := OrderPerm (p);
   Perms := [];
   for i in [o - 1, o - 2 .. 1] do 
      Perms[i] := ListPerm (p^i);
      #permutation must have length PolyBasis + 1
      fixed := [Length (Perms[i]) + 1..lp];
      Append (Perms[i], fixed);
   od; 

   #InfoTensor1 ("#I Perms = ", Perms, "\n");
   return Perms;

end;

# P is the list of permutations, u is one possible left-hand side; 
# if some image of u under an element of P occurs later in an 
# ordering, we don't need to process u 

ProcessVector := function (P, u, lenu)

   local i, v, Im, j;

   i := 1; 
   repeat 
      v := P[i];
      Im := [];
      for j in [1..lenu] do 
         Add (Im, u[v[j]]);
      od;
      if Im > u then return false; fi;
      i := i + 1;
   until i > Length (P);

   return true;

end;

# given sequence of polynomials, some product of which is f,
# find which exponents occur in f 

ExponentsOfFactors := function (R, f)

   local fac, exponents, i, factor, j;

   fac := Collected (Factors (f));

   exponents := List ([1..Length (R)], x -> 0);

   for i in [1..Length (fac)] do
      factor := fac[i];
      j := Position (R, factor[1]);
      exponents[j] := factor[2];
   od;

   return exponents;

end;

# setup the basis matrices and write them over the integers 

SetupMatrices := function (Table)

   local n, m, M, i, x, y;

   n := Length (Table[1]);
   m := Length (Table[1][1]);

   M := [];
   for i in [1..Length (Table)] do
      x := Table[i];
      #I don't know that we need to care 
      #      y := &cat[x[j] : j in [1..Length (x)]];
      M[i] := x;
   od;

   return M;

end;

# Perms is list of possible permutations which can be 
# used to reduce number of possible left-hand sides; 
# M is the set of matrices; t is the right-hand side;
# Degrees is list of degrees of the factors;
# DimU is the degree of the u factor; 
# build up left-hand side and solve system 

FindFactorisation := function (Perms, M, t, Degrees, DimU)

   local Outstanding, WorkHard, Resolved, m, M, h, s, K,
         NonNegative, zm, zs,  t, R, flag, K, WorkHard, 
         tot, n, lenm, m, x, index, A;

   tot := 0;

   n := Length (M);
   lenm := n + 1;
   m := List ([1..n + 1], x -> 0);
   x := 0;

   WorkHard := false;
   # can we settle the question for this element? 
   Outstanding := false; 

   repeat 
      index := 1;
      m[index] := m[index] + 1;
      x := x + Degrees[index];

      while (index <= n and x > DimU) do
         x := x - m[index] * Degrees[index];
         m[index] := 0;
         index := index + 1;
         m[index] := m[index] + 1;
         if index <= n then 
            x := x + Degrees[index];
         fi;
      od;

      if x = DimU then 
         if Length (Perms) = 0 or ProcessVector (Perms, m, lenm) then 
            t := [];
            for h in [1..n] do 
               t[h] := m[h] * M[h];
            od;
            A := Sum (t);
            tot := tot + 1; 
            R := IsSystemSolvable (A, t); 
            flag := R[1]; s := R[2]; K := R[3];

            if flag then 
               InfoTensor1 ("#I A solution over Z was found after testing", 
                            tot, "vectors of correct weight\n");
               InfoTensor1 ("#I s = ", s, "\n");
               InfoTensor1 ("#I Kernel has dimension ", Dimension (K), "\n");
               if Dimension (K) > 0 then 
                  InfoTensor1 ("#I Kernel has dimension ", Dimension (K), "\n");
               fi;

               Resolved := true;
               if ForAny (s, x -> x < 0) then
                  if Dimension (K) > 0 then 
                     if WorkHard then 
                        # we should test if some translate of this solution 
                        # is non-negative; this may be very expensive 
                        InfoTensor1 ("#I Now try for a solution over N\n");
                        R := NonNegativeSolution (A, t);
                        NonNegative := R[1]; s := R[2];
                     else 
                        Resolved := false;
                        # record one possible solution over Z 
                        NonNegative := false;
                        zm := m; zs := s;
                     fi; # if WorkHard 
                  else 
                     InfoTensor1 ("#I Solution is unique\n");
                     NonNegative := false;
                  fi; # Dimension (K) > 0 
               else
                  InfoTensor1 ("#I Our existing solution is over N\n");
                  NonNegative := true;
               fi; # if exists 

               if not Resolved then Outstanding := true; fi;

               if NonNegative then 
                  InfoTensor1 ("#I A solution over N found after testing", tot,  
                          "vectors of correct weight\n");
                  InfoTensor1 ("#I m = ", m, "s = ", s, "\n");
                  return [true, m, s, true];
               fi;

            fi; # if flag 

         fi; # if ProcessVector 

      fi; # if x = DimU 

   until (index > n);

   if Outstanding then 
      InfoTensor1 ("#I ** Existence of non-negative solution for some u unresolved\n"); 
      return [true, zm, zs, false]; 
   fi;

   InfoTensor1 ("#I Number of vectors of correct weight tested is ", tot, "\n");

   return [false, false, false, true];

end;

# compute factors of x^n - theta 

ListFactors := function (F, n, theta)

   local P, R, PolyBasis, Degrees, x;

   P := Indeterminate (F)^n - theta;
   R := Collected (Factors (P));
   R := Reversed (R);

   PolyBasis := List (R, x -> x[1]);
   Degrees := List (R, x -> Degree (x[1]));

   return [PolyBasis, Degrees];

end;

# compute tensor product of each element of PolyBasis1 with 2
# and record what combination of elements of PolyBasis
# is equal to this product    

ComputeTensorTable := function (F, PolyBasis, PolyBasis1, PolyBasis2)

   local T, i, j, tp;

   T := [];
   for i in [1..Length (PolyBasis1)] do
      T[i] := [];
      for j in [1..Length (PolyBasis2)] do
         tp := PolynomialTensorProduct (F, PolyBasis1[i], PolyBasis2[j]);
         T[i][j] := ExponentsOfFactors (PolyBasis, tp);
      od;
   od;

   return T;

end;

# f is the characteristic polynomial of an element of order n;
# does it have a tensor factorisation with a factor of degree dimU? 

DecideFactorisation := function (F, f, n, phi, n0, PolyBasis, dimU, t) 

   # run over theta where theta is an element of F^* / (F^*)^n  
   # and F^* is the multiplicative group of the field 

   local PolyBasis1, PolyBasis2, Degrees1, Degrees2,  Resolved, M, 
         omega, x, E, Outstanding, theta, Perms, R, found, u, v, Resolved;

   omega := PrimitiveElement (F);
   E := List ([0..n0 - 1], x -> omega^x);

   Outstanding := false;
   for theta in E do

      R := ListFactors (F, n, theta);
      PolyBasis1 := R[1]; Degrees1 := R[2];
      R := ListFactors (F, n, phi * theta^-1);
      PolyBasis2 := R[1]; Degrees2 := R[2];
      M := ComputeTensorTable (F, PolyBasis, PolyBasis1, PolyBasis2);
      #M := SetupMatrices (Table);

      Perms := ApplySymmetry (F, n, PolyBasis1);
      InfoTensor1 ("#I theta = ", theta,  " Degrees1 is ", Degrees1, "\n");
      R := FindFactorisation (Perms, M, t, Degrees1, dimU);
      found := R[1]; u := R[2]; v := R[3]; Resolved := R[4];
      InfoTensor1 ("#I u = ", u, " v = ", v, "\n");
      if found and Resolved then return [found, u, v, Resolved]; fi;
      if not Resolved then Outstanding := true; fi; 

   od;

   return [false, false, false, not Outstanding];

end;

# try to tensor factorise the characteristic polynomials
# of N random elements of G; Outstanding records if the 
# possible factorisation over the natural nuumbers
# is not conclusively decided 

FactorisePolynomials := function (G, N, L) 

   local MaxOrder, MaxNmrFactors, F, q, Tested, NmrElts, pair, Outstanding, 
         found, u, v, Resolved, g, R, n, phi,  f, n0, PolyBasis, Degrees, t;

   MaxOrder := 40;
   MaxNmrFactors := 24;

   F := BaseRing (G);
   q := Size (F);

   # examine tensor factorisation of characteristic polynomial 
   Tested := [];
   NmrElts := 0;

   Outstanding := false;

   repeat

      NmrElts := NmrElts + 1;
      g := PseudoRandom (G);
      R := ProjectiveOrderMat (g);
      n := R[1]; phi := R[2];

      # if the order of g is too large, replace g by some power 
      if n > MaxOrder then 
         g := PowerOfSmallOrder (g, n, MaxOrder);
         R := ProjectiveOrderMat (g);
         n := R[1]; phi := R[2];
      fi;

      f := CharacteristicPolynomial (g);

      if not (f in Tested) then
         AddSet (Tested, f);
         n0 := Gcd (n, q - 1);

         # get factors of x^n - phi and express f as a product 
         # of its irreducible factors 
         R := ListFactors (F, n, phi);
         PolyBasis := R[1]; Degrees := R[2];

         if Length (Degrees) > MaxNmrFactors then 
            InfoTensor1 ("#I ** Too many factors -- choose another element **\n"); 
         else 
            InfoTensor1 ("#I Projective order of element is ", n, "\n");

            t := ExponentsOfFactors (PolyBasis, f);
            #         t := RSpace (Integers (), Length (t)) ! t;

            for pair in L do
               InfoTensor1 ("#I Processing dimension factorisation ", pair, "\n");
               # do we find tensor product factorisation? 
               R := DecideFactorisation (F, f, n, phi, n0, 
                            PolyBasis, pair[1], t);
               found := R[1]; u := R[2]; v := R[3]; Resolved := R[4];

               InfoTensor1 ("#I Resolved is ", Resolved, "\n");

               #temporary -- we have no IsConsistent 
               #            if (not found) and Resolved then 
               #               L := Filtered (L, x -> x <> pair);
               #            fi;

               if not Resolved then Outstanding := true; fi;
            od;
         fi;
      fi;

   until Length (L) = 0 or NmrElts >= N;

   return [Outstanding, L];

end;
