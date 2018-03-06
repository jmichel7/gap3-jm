#############################################################################
##
#A  Matrix package                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
##
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
## 
#############################################################################

# return set of proper divisors of n 
ProperDivisors := function (n)

   return Difference (DivisorsInt (n), [1, n]);

end;

# is n a prime-power? 

IsPrimePower := function (n)

   return Length (PrimePowersInt (n)) = 2;

end; 

# is n a power of m? 

IsAPower := function (n, m)

   return m^(LogInt (n, m)) = n;

end;

# return distinct primes which occur in factorisation of n 

DistinctPrimes := function (n)

   return Set (Factors (n));

end;

# compute the inverse mod d of each element of S 

InverseSet := function (d, S)

   local x;

   return Set (List (S, x -> d / x));

end;

# given n and prime p, write n = p^alpha * m, return alpha and m 

PairOfFactors := function (n, p)

   local m, primes, index, alpha, i, f;

   # find p'-part of n 

   f := Collected (Factors (n));
   primes := List (f, i -> i[1]);

   index := Position (primes, p);
   if index > 0 then
      alpha := f[index][2];
      Remove (f, index);
   else
      alpha := 0;
   fi;

   # p'-part of n 
   m := FactorsToInt (f);

   return [alpha, m];

end;

# return factorisations of n into ordered pairs 

FactorList := function (n)

   local Pairs, D, x, y;

   Pairs := [];
   D := DivisorsInt (n);
   for x in D do for y in D do
      if x * y = n and x <= y then 
         Add (Pairs, [x, y]);
      fi;
   od; od;
   Sort (Pairs);
   return Pairs;

end;

# find co-prime factorisations of n 

CoPrimeFactorisations := function (n)

   local L, x;

   # find co-prime factorisations of n 
   L := FactorList (n);
   return Filtered (L, x -> Gcd (x[1], x[2]) = 1);

end;


