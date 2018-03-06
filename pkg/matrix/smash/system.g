#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
## 
###############################################################################
##
#F  LargestMovedPoint ( P )  . . . . . . . . .  degree of permutation group P
##
LargestMovedPoint := function (P) 
   
   if IsList (P) and IsPerm (P[1]) then
      P := Group (P, P[1]^0);
   fi;

   if IsPermGroup (P) then 
      return DegreeOperation (P, [1..PermGroupOps.LargestMovedPoint (P)]);
   else
      return 0;
   fi;

end; #LargestMovedPoint

###############################################################################
##
#F  Root ( F )  . . . . . . . . . . . . . . . . . . . . . . .  root of a field 
##

Root := function (F)

   return F.root;

end; 

PrimitiveElement := Root;

PolCoefficients := function (f)

   return f.coefficients;

end;

BaseRingZero := function (f)

   return f.baseRing.zero;

end;

OrderKnownDividend := function (R, genpol, minpol, pp)

   return R.operations.OrderKnownDividend (R, genpol, minpol, pp);

end; 
