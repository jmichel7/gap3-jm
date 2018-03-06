#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @(#)$Id: misc.g,v 1.1 1997/03/10 13:52:35 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: misc.g,v $
#H  Revision 1.1  1997/03/10 13:52:35  gap
#H  VERSION 1.0
#H
#H  Revision 1.4  1997/01/05 10:49:31  fceller
#H  added Eamonn's new version to the reprository
#H
#H  Revision 1.3  1996/12/22 07:48:12  fceller
#H  new smash version
#H
#H  Revision 1.2  1996/12/12 10:51:25  fceller
#H  new version by Eamonn
#H
#H  Revision 1.1  1996/11/28 13:14:51  fceller
#H  added "smash" and "reducible" to the repository
#H
##
###############################################################################
##
## Various functions 
## 
###############################################################################
##
#F  FactorsToInt (A) . . .    given [a, b, c, d, ...] return a^b * c ^d * ...
##
FactorsToInt := function (A)

   local result, i;

   result := 1;
   for i in [1..Int (Length (A) / 2)] do
      result := result * A[2 * i - 1]^A[2 * i];
   od;

   return result;

end; #FactorsToInt

#############################################################################
##
#F  IsScalar (A)  . . . . . . . . . . .                   is matrix A scalar?
## 
IsScalar := function(A)local i;

   if IsMat (A) = false then 
      Error ("Input should be a matrix\n");
   fi;

   if IsDiagonalMat (A) = false then return false; fi;

   for i in [2..Length (A)] do
      if A[i][i] <> A[i - 1][i - 1] then return false; fi;
   od;

   return true;

end; #IsScalar

#############################################################################
##
#F  InverseMod (u, p)  . .
##  p is prime, 0 < u < p
##  The function returns the positive integer x less than p such that
##  ux = 1 mod p.
##
InverseMod := function (u,p)
   local x;
   x := GcdRepresentation(u,p)[1];
   if x > 0 then return x; else return p + x; fi;
end;


#############################################################################
##
#F  IntPower ( A, a)  . . . if A = a^k, for k >= 1 return k, otherwise false
## Also return false if either A or a is not a positive integer.
##
IntPower := function  ( A, a)

  local N, k;

  if A < 0  or a < 0 or A < a then return false; fi;
  N := A; k := 0;
  repeat N := QuoInt (N, a); k := k + 1; until RemInt (N, a) <> 0;
  if N = 1 then return k; else return false; fi;

end;

#############################################################################
##
#F  Identity (G)  . . . . . . . . . . . . return identity element of G 
#  JM suppress by changing function GModule() and using generic Identity
##
