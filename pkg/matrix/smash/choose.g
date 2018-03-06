#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @(#)$Id: choose.g,v 1.1 1997/03/10 13:52:24 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: choose.g,v $
#H  Revision 1.1  1997/03/10 13:52:24  gap
#H  VERSION 1.0
#H
#H  Revision 1.3  1997/01/05 10:49:19  fceller
#H  added Eamonn's new version to the reprository
#H
#H  Revision 1.2  1996/12/10 12:05:42  fceller
#H  added new versions to the repository
#H
#H  Revision 1.1  1996/11/28 13:14:44  fceller
#H  added "smash" and "reducible" to the repository
#H
##
#############################################################################
#
# function to select various kinds of random elements 
#
# G can be a group or a G-module 
# 
#############################################################################
##
#F ElementOfOrder (G, RequiredOrder, NmrTries) . . . . . 
##
## try to find element of RequiredOrder in NmrTries attempts 
## 
ElementOfOrder := function (G, RequiredOrder, NmrTries)

   local g, o, i, rem;

   i := 0;
   repeat
      g := PseudoRandom (G); 
      o := OrderMat (g); 
      i := i + 1;
      rem := Mod (o, RequiredOrder) = 0; 
   until rem = true or i > NmrTries;
   
   if rem then 
      return g^(o / RequiredOrder);
   else 
      return false;
   fi;

end; #ElementOfOrder

#############################################################################
##
#F LargestPrimeOrderElement (G, NmrTries) . . . . 
## 
## generate NmrTries random elements and check which element of prime order 
## can be obtained from powers of each; return largest and its order 
##
LargestPrimeOrderElement := function (G, NmrTries)

   local i, g, o, p, prime, max, pos, z;
   
   g := [];
   o := [];
   prime := [];

   for i in [1..NmrTries] do 
      g[i] := PseudoRandom (G);
      o[i] := OrderMat (g[i]);
      #store maximum prime which occurs in factorisation of this order 
      prime[i] := Maximum (Set (FactorsInt (o[i])));
   od;

   #maximum prime encountered 
   p := Maximum (prime);
   pos := Position (prime, p);
    
   z := g[pos]^(o[pos] / p);
   return [z, p];

end; #LargestPrimeOrderElement
 
#############################################################################
##
#F LargestPrimePowerOrderElement (G, NmrTries) . . . . 
## 
## generate NmrTries random elements and check which element of prime power
## order can be obtained from powers of each; return largest and its order 
##
LargestPrimePowerOrderElement := function (G, NmrTries)

   local i, g, pow, powers, o, prime, max, pos, z, j;
   
   g := [];
   o := [];
   prime := [];

   for i in [1..NmrTries] do 
      g[i] := PseudoRandom (G);
      o[i] := OrderMat (g[i]);

      powers := PrimePowersInt (o[i]);

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

   #find position of element of largest prime-power order in list 
   pos := Position (prime, Maximum (prime));
    
   z := g[pos]^(o[pos] / prime[pos]);
   return [z, prime[pos]];

end; #LargestPrimePowerOrderElement
 
#############################################################################
##
#F ChooseRandomElements (G, NmrElts) . . . choose NmrElts random elements
##
ChooseRandomElements := function (G, NmrElts)

   local w, i;

   w := [];
   for i in [1..NmrElts] do 
      w[i] := PseudoRandom (G);
   od;

   return w;

end; #ChooseRandomElements

#############################################################################
##
#F ElementWithCharPol (G, o, pol, Limit) . . 
##
## find element of order o with characteristic polynomial pol; 
## Limit is number of tries to find such element
##
ElementWithCharPol := function (G, o, pol, Limit)

   local NmrTries, z, cf, ValidElement;

   NmrTries := 0;
   repeat
      z := ElementOfOrder (G, o, Limit);
      if z = false then
         Print (Limit, " tries failed to find element of order ", o, "\n");
         return false;
      fi;
      cf := CharacteristicPolynomial (z);
      ValidElement := (pol = cf);
      NmrTries := NmrTries + 1;
   until ValidElement or NmrTries > Limit;

   if ValidElement = false then
      Print (Limit, " attempts failed to find an element of order ", o,
                    "with given characteristic polynomial\n");
      return false;
   fi;

   return z;

end; #ElementWithCharPol

##############################################################################
##
#F  RandomConjugate (matrices, S) . . .
## 
##  return the conjugate of a random element of the set S by 
## a random element of the group G  
## 
RandomConjugate := function (G, S)

   local g, h;

   if Length (S) = 0 then Error ("S is empty"); fi; 

   g := PseudoRandom (G);
   h := S[Random ([1..Length (S)])];
   return g^-1 * h * g;

end;

##############################################################################
##
#F AddRandomTranslatingConjugate (G, S, W, F) . . .
##
## Find a random conjugate  of a random element of smatrices by a random
## element of G that doesn't fix the subspace W, and add it to smatrices
## 
AddRandomTranslatingConjugate := function (G, S, W, F)

   local dimW, g;

   if Length (S) = 0 then Error ("S is empty"); fi;

   dimW := Length (W);
   repeat 
      g := RandomConjugate (G, S);
   until Length (Base (RowSpace (Concatenation (W, W * g), F))) > dimW;

   Add (S, g);

end;
   
#############################################################################
##
#F Commutators (matrices) . . . . . form non-trivial commutators of matrices 
##
Commutators := function (matrices)

   local matrices, S, i, j, k, g, id; 

   S := [];
   k := Length (matrices);
   if k < 2 then return S; fi;
   id := matrices[1]^0;
   for i in [1..k] do
     for j in [i + 1..k] do
       g := Comm (matrices[i], matrices[j]);
       if g <> id then AddSet (S, g); fi;
     od;
   od;
 
   return S;
end;
