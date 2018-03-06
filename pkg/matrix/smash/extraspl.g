#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @ (#)$Id: extraspl.g,v 1.1 1997/03/10 13:52:29 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: extraspl.g,v $
#H  Revision 1.1  1997/03/10 13:52:29  gap
#H  VERSION 1.0
#H
#H  Revision 1.2  1997/01/05 10:49:25  fceller
#H  added Eamonn's new version to the reprository
#H
#H  Revision 1.1  1996/12/25 09:03:43  fceller
#H  changed long filenames to MS-DOS conform filenames,
#H  the init files are *NOT* yet updated
#H
#H  Revision 1.3  1996/12/22 07:48:08  fceller
#H  new smash version
#H
#H  Revision 1.2  1996/12/12 10:51:23  fceller
#H  new version by Eamonn
#H
#H  Revision 1.1  1996/11/28 13:14:47  fceller
#H  added "smash" and "reducible" to the repository
#H
##
##
############################################################################
##
#F  InfoExtraSpecial  (...)  . . . . . . . . . . . for debugging assistance
##
##
if not IsBound (InfoExtraSpecial)  then InfoExtraSpecial := Ignore;  fi;
#############################################################################
##
#F  ExtraSpecialDecomposition (module, S)  . .  
##  module is a module for a finite matrix group G over a field of 
##  characteristic p.
##  S is a set of invertible matrices, assumed to act absolutely irreducibly 
##  on the underlying vector space of module.
##  Extraspecial returns true if (mod scalars)
##  S is an extraspecial r-group, for some prime r, 
##  or a 2-group of symplectic type (that is, the central product of an 
##  extra-special 2 group with a cyclic group of order 4), normalised by G.
##  Otherwise it returns false.
##
##  It attempts to prove that S is extra-special or of symplectic type by
##  construction. That is, it tries to find elements 
##  x[1], y[1], x[2], y[2], ..., x[k], y[k], z, which generate
##  S (mod scalars) so that the x[i]'s and y[i]'s 
##  are non-central in <S>, z is central and scalar of order r, 
##  x[i] and x[j] commute for all i, j, 
##  y[i] and y[j] commute for all i, j, 
##  x[i] and y[j] commute for distinct i, j, 
##  but the commutator of x[i] and y[i] is equal to z for all i.
##  Such generators  are found by adjusting the set S. 
##  x[1] and y[1] can be chosen as any two non-commuting non-scalar elements of
##  S  (if there are none then S cannot be extraspecial/of symplectic type), 
##  and z is set equal to their commutator.
##  Each other generator x[i] or y[i] is equal to a different non-scalar 
##  element of S multiplied by a word in those x[j]'s and y[j]'s 
##  already selected, or a power of such. 
##  More specifically, the function Stripped, applied to a generator in S, 
##  multiplies it by a word in x[1], ..x[i-1], y[1], ..y[i-1] until the 
##  product commutes with all of those generators. We find x[i] and y[i], 
##  if possible, by finding two such elements, both non-scalar, 
##  and non-commuting. x[i] is the first, and y[i] is a non-scalar power 
##  of the second such that [x[i], y[i]] = z.
##  If the construction fails at any stage, the function
##  returns false. Otherwise, the construction stops with S exhausted. 
##  Finally conjugates of each x[i] and y[i] by generators of G are computed, 
##  and checked for membership of <S> using the function Stripped.
##  If some conjugate is not in S false is returned, otherwise the function
##  returns true.
##  If the function returns true, the flag .extraSpecialGroup is set to S, 
##  and .extraSpecialPrime to r.
## 
ExtraSpecialDecomposition := function (module, S) 

   local d, F, matrices, r, m, primes, p, x, y, z, doneGen, 
         ngensG, ngensS, i, j, id, g, h, comm, u;

   if IsGModule (module) = false then
      Error ("usage: ExtraSpecialDecomposition (<module>, <matrices>)");
   fi;

   d := DimensionFlag (module);

   primes := PrimePowersInt (d);
   if Length (primes) > 2 then 
      InfoExtraSpecial ("#I Dimension ", d, " is not a prime power.\n");
      return false;
   fi;

   if Length (S) = [] then return false; fi;

   #d = r^m 
   r := primes[1];
   m := primes[2];

   F := FieldFlag (module);
   matrices := GeneratorsFlag (module);

   ngensG := Length (matrices);
   ngensS := Length (S);
   x := []; y := []; doneGen := [];

   for j in [1..ngensS] do 
      h := S[j];
      if IsScalar (h) then 
         doneGen[j] := true; 
      elif IsScalar (h^r) then 
         if x = [] then 
            Add (x, h); doneGen[j] := true; 
         else 
            doneGen[j] := false; 
         fi;
      else
         InfoExtraSpecial ("#I ", r, "-th power of matrix is not scalar.\n");
         InfoExtraSpecial ("#I ExtraSpecialDecomposition returns false.\n");
         return false;
      fi;
   od;

   if x = [] then
      Error ("All matrices in S are scalar.");
   fi;

   j := 1;
   id := x[1]^0;
   while y = [] and j <= ngensS do
      if doneGen[j] = false then
         z := Comm (x[1], S[j]);
         if z <> id then
            if z^r <> id then 
               InfoExtraSpecial (r, "-th power of commutator is not identity.\n");
               return false;
            else 
               Add (y, S[j]);
               doneGen[j] := true;
            fi;
         fi;
      fi;
      j := j + 1;
   od;

   if y = [] then
      Error ("Z(S) contains non-scalar, S can't act absolutely irreducibly");
   fi;

   for i in [2..m] do
      j := 1;
      while Length (x) < i and j <= ngensS do
         if doneGen[j] = false then
            h := Stripped (S[j], x, y, z, r);
            if h = false then 
               InfoExtraSpecial ("#I Strip failed, group isn't extraspecial.\n");
               return false;
            elif IsScalar (h) = false then
               if IsScalar (h^r) then
                  Add (x, h); doneGen[j] := true;
               else 
                  InfoExtraSpecial (r, "-th power of matrix is not scalar.\n");
                  return false;
               fi;
            fi;
         fi;
         j := j + 1;
      od;

      if Length (x) < i then
         Error ("Couldn't construct x[i].");
      fi;

      j := 1;
      while Length (y) < i and j <= ngensS do
         if doneGen[j] = false then
            h := Stripped (S[j], x, y, z, r);
            if h = false then 
               InfoExtraSpecial ("#I Strip failed, group isn't extraspecial.\n");
               return false;
            else
               comm := Comm (x[i], h);
               if comm <> id then 
                  u := 0;
                  while u < r and comm <> z^u do u := u + 1; od;
                  if u = r then
                     InfoExtraSpecial ("#I Commutator isn't a power of z.\n");
                     return false;
                  fi;
                  if u <> 0 then 
                     h := h^InverseMod (u, r); 
                     if IsScalar (h^r) then
                        Add (y, h); doneGen[j] := true;
                     else 
                        InfoExtraSpecial (r, "-th power of matrix is not scalar.\n");
                        return false;
                     fi;
                  fi;
               fi;
            fi;
         fi;
         j := j + 1;
      od;

      if Length (y) < i then
         Error ("Couldn't construct y[i].");
      fi;

   od;

   for j in [1..ngensS] do
      if doneGen[j] = false then
         h := Stripped (S[j], x, y, z, r);
         if h = false then 
            InfoExtraSpecial ("#I Group isn't extraspecial.\n");
            return false;
         elif IsScalar (h) = false then
            InfoExtraSpecial ("#I Group isn't extraspecial.\n");
            return false;
         fi;
      fi;
   od;

   for j in [1..m] do
      for g in matrices do
         h := Stripped (x[j]^g, x, y, z, r);
         if h = false then 
            InfoExtraSpecial ("#I Strip failed, group isn't extraspecial.\n");
            return false;
         elif IsScalar (h) = false then
            InfoExtraSpecial ("#I Group isn't normalised by G.\n");
            return false;
         fi;
         h := Stripped (y[j]^g, x, y, z, r);
         if h = false then 
            InfoExtraSpecial ("#I Strip failed, group isn't normalised by G.\n");
            return false;
         elif IsScalar (h) = false then
            InfoExtraSpecial ("#I Group isn't normalised by G.\n");
            return false;
         fi;
      od;
   od;

   SetExtraSpecialFlag (module, true);
   SetExtraSpecialGroupFlag (module, Group (S, S[1]^0));
   SetExtraSpecialPrimeFlag (module, r);
   return true;

end;

#############################################################################
##
#F  Stripped (h, x, y, z, r)  . .  
## 
## r is a prime, h and z are invertible non-scalar matrices and x, y lists 
## of invertible matrices all of the same dimension over some finite field.
## y has length t, and x length t or t + 1.
## The function returns a matrix k equal to a product of h with a word in
## powers of x[1], y[1], ..., x[t], y[t]  (in that order), and with the 
## property that k commutes with each of x[1], y[1], ..., x[t], y[t].
## 
Stripped := function (h, x, y, z, r) 

   local k, t, u, v, c1, c2, i;

   t := Length (y);
   k := Copy (h);

   for i in [1..t] do
      c1 := Comm (x[i], k);
      c2 := Comm (y[i], k);
      u := 0;
      while u < r and c1 <> z^u do u := u + 1; od;
      if u = r then
         InfoExtraSpecial ("#I Commutator isn't a power of z.\n");
         return false;
      fi;
      v := 0;
      while v < r and c2 <> z^v do v := v + 1; od;
      if v = r then
         InfoExtraSpecial ("#I Commutator isn't a power of z.\n");
         return false;
      fi;

      k := k * x[i]^v * y[i]^(r - u);
   od;

   return k;

end;
