#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @(#)$Id: stabilis.g,v 1.1 1997/03/10 13:52:42 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: stabilis.g,v $
#H  Revision 1.1  1997/03/10 13:52:42  gap
#H  VERSION 1.0
#H
#H  Revision 1.2  1997/01/05 10:49:38  fceller
#H  added Eamonn's new version to the reprository
#H
#H  Revision 1.1  1996/12/25 09:07:39  fceller
#H  changed long filenames to MS-DOS conform filenames,
#H  the init files are *NOT* yet updated
#H
#H  Revision 1.3  1996/12/22 07:48:18  fceller
#H  new smash version
#H
#H  Revision 1.2  1996/12/10 12:05:46  fceller
#H  added new versions to the repository
#H
#H  Revision 1.1  1996/11/28 13:14:59  fceller
#H  added "smash" and "reducible" to the repository
#H
##
#############################################################################
#
#code for block stabiliser test -- part of primitivity testing
#
#note that all references to order are to projective order
#############################################################################
##
#F  InfoPrim (...)  . . . . . . . . . . . information function for package
#F  InfoPrim2 (...)  . . . . . . . . . . . additional information function
##
##
if not IsBound(InfoPrim)  then InfoPrim := Ignore;  fi;
if not IsBound(InfoPrim2)  then InfoPrim2 := Ignore;  fi;

#############################################################################
##
#F FindExpressions (s, d, m)  . . . . . . find all linear combinations of the
## integer  s  in  terms of the  d[i],  where each coefficient of the d[i] 
## is bounded above by m[i] 
##
## it is assumed that m[i] is the bound on the number of occurrences of d[i]
## 
FindExpressions := function (s, d, m)

   local R, len, solutions, i, C, j, good;

   R := RestrictedPartitions (s, d);

   #check whether the number of occurrences of each entry in a
   #partition is at most the corresponding entry in the sequence m

   solutions := []; 

   for i in [1..Length (R)] do
      C := Collected (R[i]);

      len := Length (C);
      InfoPrim2 ("#I A possible combination is ", C, "\n");
      j := 1;

      repeat
         good := (C[j][2] <= m[ Position (d, C[j][1]) ]);
         j := j + 1;
      until good = false or j > len;

      if good then 
         Add (solutions, C);
      fi;
   od;

   return solutions;

end; #FindExpressions

#############################################################################
##
#F SortCFs (CF) .  . . . . sort the composition factors of each dimension by
## decreasing multiplicity
## 
SortCFs := function (CF) 

   local Result, R, d, m, dd, mult, x, i, j;

   d := [];
   m := [];
  
   R := InvariantsOfCF (CF);
   for i in [1..Length (R)] do
     d[i] := R[i][1];
     m[i] := R[i][2];
   od;

   Result := [];

   dd := Set (d); 
   for i in [1..Length (dd)] do
      x := [];
      mult := [];
      for j in [1..Length (CF)] do
         if DimensionFlag (CF[j][1]) = dd[i] then
            Add (x, CF[j]);
            Add (mult, m[j]);
         fi;
      od;
      SortParallel (mult, x);
      x := Reversed (x);
      Append (Result, x);
   od;

   return Result;

end;

#############################################################################
##
#F InvariantsOfCF (CF) . . . return list of invariants of composition factors 
## 
## for each isomorphism type, a list containing 
## [dimension, number of factors of this dimension]
## 
InvariantsOfCF := function (CF)

   local invariants, i, dim, mult;

   invariants := [];
   for i in [1..Length (CF)] do
      dim := DimensionFlag (CF[i][1]);
      mult := CF[i][2];
      invariants[i] := [dim, mult];
   od;

   return invariants;

end; #InvariantsOfCF

#############################################################################
##
#F FindSolutions (CF, s)  . . . . . .   
## 
## CF is the list of composition factors 
## set up a list of invariants -- for each composition factor,
## its dimension and the number of occurrences of this factor 
##
## at the same time, set up sequence d which stores the distinct 
## dimensions of factors; and a sequence m which counts number of 
## factors of each distinct dimension 
## 
## now use FindExpressions to find all solutions a[i] to the equation
##
##       sum of a[i] * d[i] = s
## 
## where a[i] <= m[i]
## 
## return invariants, set of solutions, d, m
## 
FindSolutions := function (CF, s)

   local R, m, d, i, invariants, dim, mult, pos;

   m := []; d := []; invariants := [];
   for i in [1..Length (CF)] do
      d[i] := 0;
      m[i] := 0;
      dim := Dimension (CF[i][1]);
      mult := CF[i][2];
      invariants[i] := [dim, mult];
      pos := Position (d, dim); 
      if pos <> false then 
         m[pos] := m[pos] + mult;
      else 
         d[i] := dim;
         m[i] := mult;
      fi;
   od;

   d := Filtered (d, x -> x <> 0);
   m := Filtered (m, x -> x <> 0);

   InfoPrim2 ("#I d = ", d, "\n");
   InfoPrim2 ("#I m = ", m, "\n");
   InfoPrim ("#I Invariants = ", invariants, "\n");

   R := FindExpressions (s, d, m);
   InfoPrim ("#I Result of find_solutions is ", R, "\n");

   return [invariants, R, d, m];

end; #FindSolutions

########################################################################
##
#F TestSolution (invariants, solution, d, m, combined, check) . . . . .
##
## run two tests on each component of the solution 
## if no component of the solution passes either test, then we must 
## retain the submodule
##
## if any component of the solution passes either test, we note relevant
## composition factors and later obtain a minimal submodule containing
## this composition factor to hand to MinBlocks; since the only
## requirement for MinBlocks is that we have a space contained in it,
## it is sufficient to find any space contained in the block;
## we do not need to find all
## 
TestSolution := function (invariants, solution, d, m, combined, check)

   local NmrComponents, index, retain, j, D, C, k, len, 
         case2, reqnmr, nmr, dim;

   #for each solution, take each component 
   #
   #     [D, C] 
   #
   #and see if it satisfies either of the following:
   #
   #  (i) are there C factors of dimension D?
   # (ii) do all factors of dimensions D occur with multiplicity 1?
   #
   #if the component of the solution satisfies either of these tests, 
   #then note relevant composition factors
   #
   #if no component of the solution satisfies either of these
   #tests, we cannot settle definitely whether the associated 
   #minimal submodule may be contained in a block so we must 
   #save this submodule for later processing  
   
   len := Length (invariants);

   retain := true;
   NmrComponents := Length (solution); 

   j := 0;
   while j < NmrComponents and retain do 

      j := j + 1;

      #consider component [D, C]
      D := solution[j][1];
      C := solution[j][2];

      #first check if the coefficent C equals the sum of the
      #multiplicities for all factors of this dimension D

      if Position (combined, solution[j]) <> false then  
         retain := false; 
         InfoPrim ("#I Component [", D, ", ", C, "] passed Test 1\n");
         #note index of any one of the composition factors of dimension D
         k := 0;
         repeat
            k := k + 1; 
         until (invariants[k][1] = D) or k = len;
         AddSet (check, k);
          
      else 

         #check if the multiplicity of each factor of dimension D is 1
         case2 := true;
         k := 1;
         repeat
            dim := invariants[k][1];
            if dim = D then 
               case2 := (invariants[k][2] = 1);
            fi;
            k := k + 1;
         until case2 = false or dim > D or k > len;

         #if so, we need to take at most 
         #
         #  (number of factors of dimension D) - C + 1 
         #
         #CFs before finding some minimal submodule containing this block 

         if case2 = true then 
            InfoPrim ("#I Component [", D, ", ", C, "] passed Test 2\n");
            retain := false;
            index := Position (d, D);
            reqnmr := m[index] - C + 1; 
            k := 1; nmr := 0;
            repeat 
               if invariants[k][1] = D then 
                  AddSet (check, k);
                  nmr := nmr + 1;
               fi;
               k := k + 1;
            until nmr = reqnmr;
         fi; 

      fi; #else Position = false

      if retain = true then 
         InfoPrim ("#I Component [", D, ", ", C, "] fails both tests\n");
      fi;

   od; #while j 

   return retain;

end; #TestSolution

#############################################################################
##
#F ProcessLattice (M, module, CF, MaxNmr) . . .
##
## search for at most MaxNmr minimal submodules of module isomorphic to CF 
## and hand each to MinBlocks 
## 
ProcessLattice := function (M, module, CF, MaxNmr)

    local NmrSubmods, Submods, NmrBlocks, BT, i;

    Submods := MinimalSubGModules (CF, module, MaxNmr);
    NmrSubmods := Length (Submods);

    #did we find all minimal submodules isomorphic to CF?
    #if not, we do not process these
    if (NmrSubmods >= MaxNmr) then return false; fi;

    #hand each minimal submodule to MinBlocks 
    InfoPrim2 ("#I Number of submodules constructed ", NmrSubmods, "\n");
    for i in [1..NmrSubmods] do
       TriangulizeMat (Submods[i]);
       BT := MinBlocks (M, Submods[i]);
       NmrBlocks := NumberBlocksFlag (BT);
       if NmrBlocks > 1 then
          return BT;
       #else
       #   InfoPrim ("#I Number of blocks is 1\n");
       fi;
    od;

    return true;

end; #ProcessLattice 

########################################################################
##
#F FindMinimalSubGModules (M, module, CF, s, iteration)  . . . . . .   
## 
## M is irreducible module for matrix group 
## module is the reducible module for subgroup
## CF is result of CompositionFactors for the reducible module
## s is the block size 
## 
## if non-trivial block system found, it is returned,
## else the procedure returns value of RetainSubgroup
## 
FindMinimalSubGModules := function (M, module, CF, s, iteration)

    local R, invariants, solutions, d, m, combined, retain, i, j, S, check,
          MaxNmr, index, NmrBlocks, RetainSubgroup, BT, AllAtMostOne, HM,
          HomBasis, BasisLength, AllSubmodulesExhausted, SolnSizes, q,
          Remain, Resolved;

    RetainSubgroup := false; #do not keep this subgroup 

    R := FindSolutions (CF, s);
   
    invariants := R[1];
    solutions := R[2];

    if Length (solutions) = 0 then return RetainSubgroup; fi;

    #set up combined 
    d := R[3];
    m := R[4];
    combined := [];
    for i in [1..Length (d)] do
       combined[i] := [];
       combined[i][1] := d[i];
       combined[i][2] := m[i];
    od;

    #when all solutions are processed, take the set of invariants noted,
    #and for each such invariant, find a minimal submodule of 
    #module containing this composition factor and hand it to MinBlocks 

    #indices for invariants to note
    check := [];

    #record resolved solutions 
    Resolved := [];

    for i in [1..Length (solutions)] do 

       retain := TestSolution (invariants, solutions[i], d, m, combined, check);

       #does any component in the solution fail both case 1 and case 2?
       #if so, we try to apply test 3 to complete solution

       if retain = true then

          #if block size s is (d_1, 1) then apply case 3 -- 
          #compute a basis for the homomorphisms from each composition
          #factor of dimension s to the module;
          #if the length of the basis is at most one for each such 
          #composition factor, we may apply MinBlocks to each basis;
          #if no block system found, we rule out module 

          if (s = solutions[i][1][1]) then
             HM := [];
             AllAtMostOne := true;

             InfoPrim ("#I Applying test 3\n");

             for j in [1..Length (CF)] do

                if (invariants[j][1]) = s then 

                   HomBasis := HomGModule(CF[j][1], module);
                   BasisLength := Length (HomBasis);

                   #store the length of the basis as part of the 
                   #invariant sequence for this submodule 
                   invariants[j][3] := BasisLength;

                   InfoPrim ("#I Length of HomBasis for CF[", j, "] is ",
                              BasisLength, "\n");
                   if BasisLength > 1 then AllAtMostOne := false; fi;
                   if BasisLength = 1 then Add (HM, HomBasis); fi;
                fi;
             od;

             #if all bases for the homomorphisms from each CF[i] to module
             #have dimension at most one, we run MinBlocks on each basis 

             if AllAtMostOne then 
                InfoPrim ("#I HomBasis Length is <= 1 for all relevant CFs\n");
                InfoPrim ("#I Hence run MinBlocks on each HomBasis\n");

                AddSet (Resolved, i);

                #apply MinBlocks to each HM[j] 
                for j in [1..Length (HM)] do
                   BT := MinBlocks (M, HM[j][1]);
                   NmrBlocks := NumberBlocksFlag (BT);
                   if NmrBlocks > 1 then 
                      return BT;
                   else 
                      InfoPrim ("#I Number of blocks from HM[", j,"] is 1\n");
                   fi; 
                od; 
             fi;
          else 
             InfoPrim ("#I No component of solution passed test 1 or 2\n");
             InfoPrim ("#I Complete solution did not pass test 3\n");
          fi; #test 3
       else 
          AddSet (Resolved, i);
       fi; #if retain 
          
    od; #for i 

    #run through all composition factors whose indices occur in check;
    #compute a minimal submodule of module containing this composition
    #factor and hand this submodule to MinBlocks 

    check := Set (check);

    for i in [1..Length (check)] do
       index := check[i];
       InfoPrim ("#I Now computing submodule for CF[", index, "] \n");
       Distinguish (CF, index);
       S := MinimalSubGModule (module, CF, index);
       BT := MinBlocks (M, S);
       NmrBlocks := NumberBlocksFlag (BT);
       if NmrBlocks > 1 then 
          return BT;
       else 
          InfoPrim ("#I Number of blocks is 1\n");
       fi; 
    od; 

    #now, if necessary, apply Desperate strategy;
    #for each composition factor whose dimension occurs in a solution, 
    #we compute some minimal submodules of module which are isomorphic 
    #to this factor; if for each factor we compute all such minimal 
    #submodules and hand each to MinBlocks without finding a block 
    #system we do not need to retain this module
    
    Remain := Difference ([1..Length (solutions)], Resolved);

    if Length (Remain) <> 0 then 

       InfoPrim ("#I Now trying ProcessLattice\n");

       #what dimensions occur in solutions?
       SolnSizes := [];
       for i in [1..Length (solutions)] do
          for j in [1..Length (solutions[i])] do
             AddSet (SolnSizes, solutions[i][j][1]);
          od;
       od;
       InfoPrim2 ("#I Dimensions in solutions are ", SolnSizes, "\n");

       #set the maximum number of minimal submodules to compute for 
       #each factor to be q + 1 or 50 whichever is greater
       q := Size (FieldFlag (M));
       MaxNmr := 50 * iteration;
       MaxNmr := Maximum (q + 2, MaxNmr + 1);

       #process all composition factors of dimension an element of SolnSizes; 
       #abort the process if for any composition factor we have failed to 
       #compute all the minimal submodules isomorphic to this factor

       for j in [1..Length (CF)] do
          if (DimensionFlag (CF[j][1]) in SolnSizes) then 
             InfoPrim2 ("#I Calling ProcessLattice with CF[", j, "]\n");
             R := ProcessLattice (M, module, CF[j][1], MaxNmr);

             #ProcessLattice might find a block system or exhaust all 
             #minimal submodules isomorphic to this composition factor

             if IsBlockSystem (R) then return R; fi;
             if R = false then return invariants; fi;
          fi;
       od;

    fi; #Length (Remain) ne 0
             
    return false;

end; #FindMinimalSubGModules 

#############################################################################
##
#F ExtractLargestPrimeOrderElement (Elts, Orders) . . . . 
## 
## check which element of prime order can be obtained from powers of the 
## supplied elements; return element of largest prime order found 
## 
ExtractLargestPrimeOrderElement := function (Elts, Orders)

   local i, p, prime, max, pos, z, NmrTries;
   
   prime := [];

   NmrTries := Length (Elts);
   for i in [1..NmrTries] do 
      #store maximum prime which occurs in factorisation of this order 
      prime[i] := Maximum (Set (FactorsInt (Orders[i])));
   od;

   InfoPrim2 ("#I Order sequence is ", Orders, "\n");
   InfoPrim2 ("#I Prime sequence is ", prime, "\n");

   #maximum prime encountered 
   p := Maximum (prime);
   pos := Position (prime, p);
    
   z := Elts[pos]^(Orders[pos] / p);
   return [z, p];

end; #ExtractLargestPrimeOrderElement

#############################################################################
##
#F SetupFactor (F, n, a)  . . . . . . . . . . .  set up factor x^n - a over F
##
## 
SetupFactor := function (F, n, a)

   local zero, one, factor, i;

   zero := Zero (F);
   one := One (F);
   factor := [];
   factor[1] := -one * a;
   for i in [2..n] do
      factor[i] := zero;
   od;
   factor[n + 1] := one;

   return Polynomial (F, factor);

end; #SetupFactor

#############################################################################
##
#F SetupElements (y, z, w, index, p) .  . . . . . .
##
## given z an element of order p and w 
## set up w[index] * z^k * w[j]^-1 for j = 1 .. i - 1 and k = 1..p 
## and add them to y
##
SetupElements := function (y, z, w, index, p)

   local power, inv, j, k;

   if index < 2 then return y; fi;

   #set up powers of z 
   power := [z];
   for k in [2..p] do
      power[k] := z^k;
   od;

   for j in [1..index - 1] do
      inv := w[j]^-1;
      for k in [1..p] do
         Add (y, w[index] * power[k] * inv);
      od;
   od;

   return y;

end; #SetupElements

#############################################################################
##
#F FixedPointFreeElement (M, z, p, r) . . . 
##                                         
## can we deduce that z of order p acts fixed-point freely on r blocks?
## we can if p divides r, p does not divide exp GL ( d / r, q)
## and char pol of z is (x^p - 1)^(d / p) 
## 
FixedPointFreeElement := function (M, z, p, r)

   local exp, d, F, q, pol, f, factor;

   if Mod (r, p) <> 0 then return false; fi;

   d := DimensionFlag (M);
   F := FieldFlag (M);
   q := Size (F);

   #find exponent of GL (d / r, q)
   exp := FactorsToInt (ExponentGL (d / r, q));

   InfoPrim2 ("#I d = ", d, " p = ", p, " and expGL is ", exp, "\n");

   if Mod (exp, p) = 0 then return false; fi;

   #set up pol = (x^p - 1)^(d / p) 
   factor := SetupFactor (F, p, 1);
   pol := factor^Int (d / p);

   #has element z a characteristic polynomial equal to pol?
   f := CharacteristicPolynomial (z);

   return (f = pol);

end; #FixedPointFreeElement

#############################################################################
##
#F ExtractLargestPrimePowerOrderElement (Elts, Orders, r) . . . . 
## 
## take the supplied elements and check which element of prime power
## order can be obtained from powers of each; among these, find and return 
## the element of largest prime-power order co-prime to r -- 
## if you do not want to impose this restriction, supply 1 as the value of r 
##
ExtractLargestPrimePowerOrderElement := function (Elts, Orders, r)

   local new, i, pow, powers, p, prime, max, pos, z, j, NmrTries;
   
   prime := [];

   NmrTries := Length (Elts);
   for i in [1..NmrTries] do 

      powers := PrimePowersInt (Orders[i]);

      max := 1;
      for j in [1..Length (powers) / 2] do
         pow := powers[2 * j - 1]^powers[2 * j];
         if pow > max then 
            max := pow;
         fi;
      od;
         
      #store maximum prime power 
      prime[i] := max;

   od;

   InfoPrim2 ("#I order sequence is ", Orders, "\n");
   InfoPrim2 ("#I prime sequence is ", prime, "\n");

   #non-trivial prime powers encountered which are coprime to r
   new := Filtered (prime, x -> x <> 1 and Gcd (x, r) = 1);

   InfoPrim2 ("#I Filtered prime sequence is ", new, "\n");

   if Length (new) = 0 then return false; fi;

   #largest non-trivial prime power encountered which is coprime to r
   p := Maximum (new);

   #find its position in original list 
   pos := Position (prime, p);
    
   z := Elts[pos]^(Orders[pos] / prime[pos]);
   return [z, p, pos];

end; #ExtractLargestPrimePowerOrderElement

#############################################################################
##
#F ProcessSubGModule (M, H, s, iteration) . . . 
## 
## Compute composition factors for submodule H; 
## if reducible, then find appropriate minimal submodules  
## 
## the function returns a block system, or a list of invariants if we 
## must keep submodule, or false if no solution nor block system found 
## 
ProcessSubGModule := function (M, H, s, iteration) 

   local R, CF;

   #compute composition factors for this module 
   CF := CompositionFactors(H);

   InfoPrim2 ("#I Length of CF is ", Length (CF), "\n");

   #if the module is reducible, can we find block system?
   if Length (CF) > 1 or CF[1][2] > 1 then 
      InfoPrim ("#I ** GModule is reducible\n");
      CF := SortCFs (CF);
      R := FindMinimalSubGModules (M, H, CF, s, iteration);
   else 
      R := false;
   fi; 

   return R;

end; #ProcessSubGModule

#############################################################################
##
#F ExamineCompositionFactors (M, S, y, s) . . . 
## 
## apply MeatAxe to various modules to get their composition factors
## 
## M irreducible module
## S collection of submodules 
## z element of order p a prime 
## s block size 
## 
## if block system found 
##    return [true, block system]
## else 
##    return [false, list T of subgroups not eliminated, list of invariants]
##
ExamineCompositionFactors := function (M, S, SInvariants, 
                                      ConjugatesRequired, z, p, s, iteration)

   local n, i, j, k, H, CF, T, RetainSubgroup, NmrTries, gens, NmrOfGens,
         w, y, R, BT, t, ReplaceElement, TInvariants, RetainGModules,
         RetainInvariants, Current;

   NmrTries := Int ((DimensionFlag (M) / s) / p) + 1;
   InfoPrim ("#I We must select ", NmrTries, " elements w[n]\n");

   T := [];
   TInvariants := [];

   w := [];

   #if we decide that a particular element w[n] is not suitable, then 
   #ReplaceElement is set true and we choose new w[n]

   ReplaceElement := false;

   n := 0;
   while n < NmrTries do
 
      n := n + 1;
      InfoPrim2 ("#I Now working with element w[", n, "]\n");

      w[n] := PseudoRandom (M);
      ReplaceElement := false;

      if ConjugatesRequired then 
         y := [w[n] * z * w[n]^-1];
      else 
         y := [];
      fi;

      t := Runtime ();
      SetupElements (y, z, w, Length (w), p);
      InfoPrim ("#I Computed ", Length (y), " elements y[j]\n");
      InfoPrim2 ("#I Time to find new y[j] is ", Runtime () - t, " ms\n");

      #modules and their invariants retained for current w[n]
      RetainGModules := [];
      RetainInvariants := [];

      i := 0;
      while i < Length (S) and ReplaceElement = false do 
         i := i + 1;
         InfoPrim ("#I Invariants of S[", i, "] is ", SInvariants[i], "\n");

         #initialise Current to contain invariants of S[i] and add 
         #to it the list of invariants retained for <S[i], y[j]>;
         #if two submodules containing S[i] have the same set of 
         #invariants, we choose a new w[n]

         Current := [ SInvariants[i] ]; 

         gens := Copy (GeneratorsFlag (S[i]));
         NmrOfGens := Length (gens);

         j := 0;
         while j < Length (y) and ReplaceElement = false do 
            j := j + 1;

            #set up H to be the module generated by S[i] and y[j] 
            gens[NmrOfGens + 1] := y[j];
            H := GModule (gens);

            InfoPrim ("#I Processing submodule ", n, " ", i, " ", j, "\n");

            R := ProcessSubGModule (M, H, s, iteration);

            #did we find block system? if so, return 
            if IsBlockSystem (R) then 
               return [true, R];
            fi;

            #do we retain the subgroup?
            if IsList (R) then 
               Sort (R);
               ReplaceElement := (R in Current);
               if ReplaceElement = false then 
                  InfoPrim2 ("#I Invariants of new are not in Current\n");
                  InfoPrim ("#I Adding subgroup to T\n");
                  Add (RetainGModules, Copy (H));
                  Add (RetainInvariants, R); 
               else 
                  InfoPrim2 ("#I Invariants of new are in Current -- choose again\n");
                  n := n - 1;
                  RetainGModules := [];
                  RetainInvariants := [];
               fi;
            fi;

         od; #for j 
      od; #for i 

      #update contents of T by adding to it modules retained  
      Append (T, RetainGModules);
      Append (TInvariants, RetainInvariants);

   od; #for n 

   InfoPrim ("#I Length of list T is ", Length (T), "\n");

   return [false, T, TInvariants];

end; #ExamineCompositionFactors
      
#############################################################################
##
#F SampleOfElements (a, b, m) . . . .
##
## select m random elements from the range [a..b]
## 
SampleOfElements := function (a, b, m)

    local Sample, i;

    Sample := [];
    if m = 0 then return Sample; fi;

    repeat 
       AddSet (Sample, RandomList ([a..b]));
    until Length (Sample) = m;
    
    return Sample;

end; #SampleOfElements

#############################################################################
##
#F SampleOfCommutators (Elts, m) . . . .
##
## calculate commutators of a random sample of m elements from Elts
## 
SampleOfCommutators := function (Elts, m)

   local Sample;

   Sample := SampleOfElements (1, Length (Elts), m);
   #InfoPrim2 ("#I Sample for commutators is ", Sample, "\n");

   if Sample = [] then return []; fi;

   Elts := Filtered (Elts, x -> Position (Elts, x) in Sample);
   if Elts = [] then return []; fi;

   return Commutators (Elts);

end; #SampleOfCommutators 

#############################################################################
##
#F SomePowerFixesBlock (M, w, o, z) . . . .
##
## w is an element of order o = p^n, a prime-power, and w^o = z; if the 
## characteristic polynomial of w is not a power of (x^o - z), then at 
## least one orbits of w does not have length o and so w^(p^(n - 1)) 
## will fix a block
##
SomePowerFixesBlock := function (M, w, o, z)

   local factor, F, f, d;

   if IsPrimePowerInt (o) = false then return false; fi;

   F := FieldFlag (M);
   factor := SetupFactor (F, o, z);

   f := CharacteristicPolynomial (w);

   d := DimensionFlag (M);
   return (f <> factor^Int( d / o));
  
end; #SomePowerFixesBlock 

#############################################################################
##
#F ChooseFirstElement (M, Elts, Orders, r) . . .
##
## find an element which fixes at least one of the r blocks 
## 
## in BlockStabiliserTest we build up the stabiliser of the block
## and want to find a generating set for this stabiliser as quickly
## as possible; hence, we also want an element which lies in as few 
## proper subgroups of the block stabiliser as we can find; we try 
## to find an element of large prime-power order which fixes a block
## 
## a non-trivial element with this property is very useful to
## BlockStabiliserTest strategy -- thus we work hard to try to find one 
## 
## if we can find an element *whose order is prime-power* and *coprime to r*
## it will fix at least one block;
##
## we first try to find such an element among powers of our currently 
## selected random elements; 
##
## if we don't find one, we first select new random elements; 
## if we still don't find one, we compute commutators of some 
## existing elements and examine these; 
## if we don't find one using any of these strategies, we then examine
## elements of prime-power order obtained as powers of the random
## elements to see whether we can deduce from characteristic polynomial 
## structure that a power of this element of prime-power order 
## fixes a block
##
## we iterate these searches for a preset number of times determined by the
## value of NmrIterations 
##
## if after all this, we cannot find one, we return the identity 
##
ChooseFirstElement := function (M, Elts, Orders, r)

   local i, j, w, o, R, powers, p, n, pow, x, g, order, factors, pos, t, z,
         NewElts, EltComms, NewOrders, CommOrders, 
         NmrElements, NmrSearches, NmrIterations, offset, 
         SampleSize, MaxSampleSize, Sample, SampleList, SampleOrders; 

   #can we find one among powers of our existing random elements?

   R := ExtractLargestPrimePowerOrderElement (Elts, Orders, r);
   if R <> false then return R; fi;

   InfoPrim ("#I No suitable w among supplied collection of elements\n");

   #these parameters determine length of search for element 
   MaxSampleSize := 10;
   NmrIterations := 3;

   NmrSearches := 0;

   repeat 
      NmrSearches := NmrSearches + 1;

      #current number of random elements 
      NmrElements := Length (Elts);

      #extend the collection of random elements 
      SampleSize := Minimum (Length (Elts), MaxSampleSize);
      NewElts := ChooseRandomElements (M, SampleSize);
      NewOrders := [];
      for i in [1..Length (NewElts)] do
         NewOrders[i] := ProjectiveOrderMat (NewElts[i])[1];
      od;
      Append (Elts, NewElts);
      Append (Orders, NewOrders);

      R := ExtractLargestPrimePowerOrderElement (NewElts, NewOrders, r);
      if R <> false then R[3] := R[3] + NmrElements; return R; fi;

      InfoPrim ("#I After selecting new random elements, no suitable w found\n");

      NmrElements := Length (Elts); 

      #if we have not found suitable first element, add in 
      #commutators of some existing random elements and try again

      SampleSize := Minimum (NmrElements, RootInt (2 * MaxSampleSize, 2) + 1);
      #t := Runtime ();
      EltComms := SampleOfCommutators (Elts, SampleSize);
      #InfoPrim2 ("#I Get ", Length (EltComms), " commutators in ", 
      #           Runtime () - t, " ms\n");

      CommOrders := [];
      for i in [1..Length (EltComms)] do
         CommOrders[i] := OrderMat (EltComms[i]);
      od;
      Append (Elts, EltComms);
      Append (Orders, CommOrders);

      R := ExtractLargestPrimePowerOrderElement (EltComms, CommOrders, r);
      if R <> false then R[3] := R[3] + NmrElements; return R; fi;

      InfoPrim ("#I After adding commutators, no suitable w found\n");

      #if we have still not found one, look at characteristic 
      #polynomial structure of elements of prime-power order 
      #and try to deduce that power of one does act with some fixed point

      #if first iteration, look at current collection of random elements 
      #else look at the new elements selected  

      if NmrSearches = 1 then
         offset := 0;
         SampleList := Elts;
         SampleOrders := Orders;
      else 
         Append (NewElts, EltComms);
         Append (NewOrders, CommOrders);
         offset := Length (Elts) - Length (NewElts);
         SampleList := NewElts;
         SampleOrders := NewOrders;
      fi;
         
      NmrElements := Length (SampleList);

      SampleSize := Minimum (NmrElements, MaxSampleSize);
      Sample := SampleOfElements (1, NmrElements, SampleSize);
      #InfoPrim2 ("#I Sample size for polynomial test is ", Length (Sample), "\n");
 
      #initialise o 
      o := 1;

      for i in [1..SampleSize] do

         InfoPrim ("#I Examining powers of element ", i, "\n");
         g := SampleList[Sample[i]];
         order := SampleOrders[Sample[i]];

         #look at prime-powers 
         powers := PrimePowersInt (order);
         for j in [1..Length (powers) / 2] do
            p := powers[2 * j - 1];
            if p > o then 
               n := powers[2 * j];
               pow := p^n;
               InfoPrim ("#I Examing characteristic polynomials of element of order ", pow, "\n");
               x := g^(order / pow);
               z := x^pow; #note scalar 
               if SomePowerFixesBlock (M, x, pow, z[1][1]) = true then 
                  w := x^(p^(n - 1));
                  o := p;
                  pos := Sample[i];
               fi;
            fi;
         od;

      od; 
         
      if o <> 1 then
         return [w, o, offset + pos];
      fi;

      InfoPrim ("#I Characteristic polynomials of elements of ppo give no suitable w\n");

   until NmrSearches = NmrIterations;

   #if all else fails, return the identity element
   InfoPrim ("#I Failed to find non-trivial element which fixes at least one of ", 
           r, " blocks \n");
   R := [Elts[1]^0, 1, 0];

   return R;

end; #ChooseFirstElement 

#############################################################################
##
#F BlockStabiliserTest (M, Elts, Orders, r) . . . 
##
## seek to eliminate block number r 
## return true if successful; otherwise return false
## if block system found, it is stored as component of M 
##
BlockStabiliserTest := function (M, Elts, Orders, r)

   local w, o, n, new, s, T, R, result, y, z, p, pos, fpf, CF,
         Invariants, iteration;

   #find element w which fixes at least one block 
   R := ChooseFirstElement (M, Elts, Orders, r);

   #set up w and note its order 
   w := R[1];
   o := R[2]; 
   InfoPrim ("#I First element w has order ", o, "\n");

   #replace the element used to provide w by a random conjugate to try 
   #to ensure that the next element z selected is different from w 
   if R[3] <> 0 then 
      pos := R[3];
      n := 0;
      repeat 
         n := n + 1;
         new := Elts[pos]^PseudoRandom (M); 
      until new <> Elts[pos] or n > 10;
      Elts[pos] := new;
   fi;

   #search for elements of prime order and return the largest, z
   R := ExtractLargestPrimeOrderElement (Elts, Orders);
   
   #element z of order p
   z := R[1];
   p := R[2];

   InfoPrim2 ("#I The element z chosen has order ", p, "\n");

   #is it possible that w = z?
   if w = z then 
      Error ("w = z in BlockStabiliserTest\n");
   fi;

   #can we deduce that z acts fixed-point freely on the r blocks?
   fpf := FixedPointFreeElement (M, z, p, r);

   #if not, does it fix every block?

   if fpf = false then 

      #call SmashGModule 
      InfoPrim2 ("#I Calling SmashGModule with element z ...\n");
      R := SmashElement (M, z, false);
      if R <> false then return R; fi;

      #not all of the r blocks are fixed; thus, if p = r, z is fpf
      fpf := (p = r);

   fi;

   #note the block size 
   s := DimensionFlag (M) / r;

   T := [ GModule ([w]) ];

   #set up invariants of first entry in T 
   InfoPrim ("#I Compute composition factors for start module ...\n");
   CF := CompositionFactors(T[1]);
   Invariants := [];
   Invariants[1] := InvariantsOfCF (CF);
   Sort (Invariants[1]);
   InfoPrim ("#I Invariants of start submodule is ", Invariants[1], "\n");

   #in ExamineCompositionFactors, find collection of elements y[i], 
   #at least one of which fixes the same block as does w; 
   #we must add to y conjugates of z by the random elements selected if we 
   #have not deduced that z acts fixed-point freely on the blocks 

   InfoPrim ("#I Conjugates required? ", not fpf,"\n");

   #now apply meataxe to various modules, iterating this until
   #the set T is empty or until we find a block system

   iteration := 0;
   repeat 
      iteration := iteration + 1;
      InfoPrim ("#I Calling MeatAxe with T of length ", Length (T), "\n");
      R := ExamineCompositionFactors (M, T, Invariants, not fpf, 
                                      z, p, s, iteration);
  
      #if first is true, second is block system
      if R[1] = true then 
         SetBlockSystemFlag (M, R[2]);
         return false;
      fi;

      #otherwise second is list of outstanding subgroups and third
      #is collection of invariants of these subgroups 
      T := R[2];
      Invariants := R[3];

   until Length (T) = 0;

   return true;

end; #BlockStabiliserTest
