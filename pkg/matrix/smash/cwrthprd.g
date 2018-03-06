#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @ (#)$Id: cwrthprd.g,v 1.1 1997/03/10 13:52:28 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: cwrthprd.g,v $
#H  Revision 1.1  1997/03/10 13:52:28  gap
#H  VERSION 1.0
#H
#H  Revision 1.2  1997/01/05 10:49:24  fceller
#H  added Eamonn's new version to the reprository
#H
#H  Revision 1.1  1996/12/25 09:03:40  fceller
#H  changed long filenames to MS-DOS conform filenames,
#H  the init files are *NOT* yet updated
#H
#H  Revision 1.2  1996/12/22 07:48:05  fceller
#H  new smash version
#H
#H  Revision 1.1  1996/11/28 13:14:43  fceller
#H  added "smash" and "reducible" to the repository
#H
#H  Revision 3.1 1996/08/20				Derek Holt
#H  ImprimitiveWreathProduct (G, P) and PowerWreathProduct (G, P) altered to 
##  work correctly where the permutation group P is intransitive.
#H  This error was found and partially corrected by J.D. Gilbey
##
#############################################################################
##
#F  ImprimitiveWreathProduct (G, P) . . 
##                         wreath product of a group and a permutation group
## 
##  G can be either  (i) a matrix group, GModule of list of matrices, 
##  or               (ii) a permutation group or list of permutations.
##  P  must be a permutation group or list of permutations.
##  ImprimitiveWreathProduct constructs the wreath product of 
##  G by P in its natural action. This is as a matrix or permutation 
##  group of degree d * t, where d is the degree of G and t the degree of P.
##  NOTE: There is a GAP library function WreathProduct, which uses the
##  regular representation of P by default.
## 
ImprimitiveWreathProduct := function (G, P)
   local i, j, k, error, matgroup, F, Ggens, nGgens, Pgens, nPgens, 
         one, d, t, gen, Wgens, Wgen, orb, orbrep;

   # First sort out the what the first argument is.
   if IsRec (G) and IsGModule (G) then
      matgroup := true;
      Ggens := GeneratorsFlag (G);
      nGgens := Length (Ggens);
      F := FieldFlag (G);
      d := DimensionFlag (G);
   else
      if IsGroup (G) then
         Ggens := Generators (G);
      else
         Ggens := G;
      fi;
      error := false;
      if IsList (Ggens) = false or Length (Ggens) = 0 then
         error := true;
      else
         nGgens := Length (Ggens);
         if IsMat (Ggens[1]) then
            matgroup := true;
            F := Field (Flat (Ggens));
            d := Length (Ggens[1]);
         elif IsPerm (Ggens[1]) then
            matgroup := false;
            d := LargestMovedPoint (Ggens);
         else
            error := true;
         fi;     
      fi;
      if error then
         return Error (
      "First argument must be a permutation or matrix group or list of permutations or matrices");
      fi;
   fi;

   # Now check and work out the degree of the second argument.
   if IsGroup (P) then
      Pgens := Generators (P);
   else
      Pgens := P;
   fi;
   error := false;
   if IsList (Pgens) = false or Length (Pgens) = 0 then
      error := true;
   else
      nPgens := Length (Pgens);
      if IsPerm (Pgens[1]) then
         t := LargestMovedPoint (Pgens);
      else     
         error := true;
      fi;
   fi;
   if error then
      return Error (
      "Second argument must be a permutation group or list of permutations");
   fi;

   one := One (F);

   # Now we can start the calculation. First do mat group case then perm group
   Wgens := [];
   if matgroup then
      # The generators consist of those of G on one representative of each
      # orbit of P on blocks, and filled out with 1's on the diagonal, 
      # and those of P as permutation matrices on blocks
      for orb in Orbits (P, [1..t]) do
         orbrep := orb[1] - 1;
         for gen in Ggens do
            Wgen := IdentityMat (d * t, F);
            Wgen{[d * orbrep + 1..d * (orbrep + 1)]}{[d * orbrep + 1..d * (orbrep + 1)]}  := 
              Copy (gen);
            Add (Wgens, Wgen);
         od;
      od;
      for gen in Pgens do
         Wgen := NullMat (d * t, d * t, F);
         for i in [1..t] do for j in [1..d] do
            Wgen[ (i - 1) * d + j][(i^gen - 1) * d + j] := one;
         od; od;
         Add (Wgens, Wgen);
      od;
   else
      # The generators consist of those of G on one representative of each
      # orbit of P on blocks, and those of P on blocks
      for orb in Orbits (P, [1..t]) do
         orbrep := orb[1] - 1;
         for gen in Ggens do
            Wgen := [1..d * t];
            for j in [1..d] do
               Wgen[orbrep * d + j] := orbrep * d + j^gen;
            od;
            Add (Wgens, PermList (Wgen));
         od;
      od;
      for gen in Pgens do
         Wgen := [];
         for i in [1..t] do for j in [1..d] do
            Wgen[ (i - 1) * d + j] := (i^gen - 1) * d + j;
         od; od;
         Add (Wgens, PermList (Wgen));
      od;
   fi;

   return Group (Wgens, Wgens[1]^0);

end;

#############################################################################
##
#F  PowerWreathProduct (G, P)  . 
##  wreath power of a group and a permutation group
## 
##  G can be either  (i) a matrix group, GModule or list of matrices, 
##  or               (ii) a permutation group or list of permutations.
##  P  must be a permutation group or list of permutations.
##  PowerWreathProduct constructs the wreath product of G by P in its power
##  action. This is as a matrix or permutation group of degree d^t, where
##  d  is the degree of  G  and  t  the degree of  P.
## 
PowerWreathProduct := function (G, P)

   local i, j, k, l, error, matgroup, F, Ggens, nGgens, Pgens, nPgens, 
        one, d, t, gen, Wgens, Wgen, orb, orbrep, dpowerrep, dpowerrep1;

   # First sort out the what the first argument is.
   if IsRec (G) and IsGModule (G) then
      matgroup := true;
      Ggens := GeneratorsFlag (G);
      nGgens := Length (Ggens);
      F := FieldFlag (G);
      d := DimensionFlag (G);
   else
      if IsGroup (G) then
         Ggens := Generators (G);
      else
         Ggens := G;
      fi;
      error := false;
      if IsList (Ggens) = false or Length (Ggens) = 0 then
         error := true;
      else
         nGgens := Length (Ggens);
         if IsMat (Ggens[1]) then
            matgroup := true;
            F := Field (Flat (Ggens));
            d := Length (Ggens[1]);
         elif IsPerm (Ggens[1]) then
            matgroup := false;
            d := LargestMovedPoint (Ggens);
         else
            error := true;
         fi;     
      fi;
      if error then
         return Error (
      "First argument must be a permutation or matrix group or list of permutations or matrices");
      fi;
   fi;

   # Now check and work out the degree of the second argument.
   if IsGroup (P) then
      Pgens := Generators (P);
   else
      Pgens := P;
   fi;
   error := false;
   if IsList (Pgens) = false or Length (Pgens) = 0 then
      error := true;
   else
      nPgens := Length (Pgens);
      if IsPerm (Pgens[1]) then
         t := LargestMovedPoint (Pgens);
      else     
         error := true;
      fi;
   fi;
   if error then
      return Error (
      "Second argument must be a permutation group or list of permutations");
   fi;

   one := One (F);

   # Now we can start the calculation. First do matrix group case, 
   # then permutation group
   Wgens := [];
   if matgroup then
      # The generators consist of those of G  repeated in blocks of d on the
      # d^t basis vectors and those of P in the product action, 
      # taking care to do this once for each orbit of P on blocks
      for orb in Orbits (P, [1..t]) do
         orbrep := orb[1];
         dpowerrep1 := d^(orbrep - 1);
         dpowerrep  := d^orbrep;
         for gen in Ggens do
            Wgen := NullMat (d^t, d^t, F);
            for i in [1..d^(t - orbrep)] do for j in [1..d] do
               for k in [1..dpowerrep1] do for l in [1..d] do
                  Wgen[ (i - 1) * dpowerrep + (j - 1) * dpowerrep1 + k]
                 [(i - 1) * dpowerrep + (l - 1) * dpowerrep1 + k] := gen[j][l];
               od; od;
            od; od;
            Add (Wgens, Wgen);
         od;
      od;
      for gen in Pgens do
         Wgen := NullMat (d^t, d^t, F);
         for i in [1..d^t] do
            Wgen[i][PowerPerm (gen, d, t, i)] := one;
         od;
         Add (Wgens, Wgen);
      od;
   else
      # The generators consist of those of G  repeated in blocks of d on the
      # d^t points and those of P in the product action, taking care to do
      # this once for each orbit of P on blocks
      for orb in Orbits (P, [1..t]) do
         orbrep := orb[1];
         dpowerrep1 := d^(orbrep - 1);
         dpowerrep  := d^orbrep;
         for gen in Ggens do
            Wgen := [];
            for i in [1..d^(t - orbrep)] do for j in [1..d] do
               for k in [1..dpowerrep1] do
                  Wgen[ (i - 1) * dpowerrep + (j - 1) * dpowerrep1 + k]  := 
                    (i - 1) * dpowerrep + (j^gen - 1) * dpowerrep1 + k;
               od;
            od; od;
            Add (Wgens, PermList (Wgen));
         od;
      od;
      for gen in Pgens do
         Wgen := [];
         # We construct the generator as a list. The images of the points are
         # calculated using PowerPerm.
         for i in [1..d^t] do
            Wgen[i] := PowerPerm (gen, d, t, i);
         od;
         Add (Wgens, PermList (Wgen));
      od;
   fi;

   return Group (Wgens, Wgens[1]^0);

end;

#############################################################################
##
#F  PowerPerm (p, d, t, pt) . induced power action of a permutation on a point
## 
##  p is a permutation on t point. Let pp be the induced action of p on
##  the set of t - tuples of d points (by permuting coordinates).
##  PowerPerm returns the action of pp on the point pt.
## 
PowerPerm := function (p, d, t, pt)
   local  i, tup, ttup, x;
   if pt <= 0 or pt > d^t then
      return Error ("pt must be in range [1..d^t].");
   fi;
   # Find the expression of pt as a t - tuple of d points.
   tup := [];
   x := pt - 1;
   for i in [1..t] do
      tup[i] := RemInt (x, d) + 1;
      x := QuoInt (x, d);
   od;
   # Now apply the permutation to tup to get ttup
   ttup := [];
   for i in [1..t] do
      ttup[i^p] := tup[i];
   od;
   # Finally calculate point to return.
   x := 1;
   for i in [1..t] do
      x := x + (ttup[i] - 1) * d^ (i - 1);
   od;
   return x;

end;
