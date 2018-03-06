#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @ (#)$Id: composit.g,v 1.1 1997/03/10 13:52:26 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: composit.g,v $
#H  Revision 1.1  1997/03/10 13:52:26  gap
#H  VERSION 1.0
#H
#H  Revision 1.2  1997/01/05 10:49:22  fceller
#H  added Eamonn's new version to the reprository
#H
#H  Revision 1.1  1996/12/25 09:03:34  fceller
#H  changed long filenames to MS-DOS conform filenames,
#H  the init files are *NOT* yet updated
#H
#H  Revision 1.2  1996/12/22 07:48:07  fceller
#H  new smash version
#H
#H  Revision 1.1  1996/12/10 12:05:43  fceller
#H  added new versions to the repository
#H
#H  Revision 1.1  1996/11/28 13:14:45  fceller
#H  added "smash" and "reducible" to the repository
#H
##
#composition.g
#
############################################################################
##
#F  InfoComposition  (...)  . . . . . . . . . . . for debugging assistance
##
##
if not IsBound (InfoComposition)  then InfoComposition := Ignore;  fi;
###############################################################################
##
#F  CompositionFactors ( module ) . . . . find composition factors of a module
##
## CompositionFactors calls IsIrreducible repeatedly to find the composition 
## factors of the GModule `module'. It also calls IsomorphismGModule to 
## determine which are isomorphic.
## It returns a list [f1, f2, ..fr], where each fi is a list [m, n], 
## where m is an irreducible composition factor of module, and n is the
## number of times it occurs in module.
## 
Smash.CompositionFactors := function ( module )
   local dim, factors, factorsout, queue, cmod, new, d, i, j, l, q;

   if IsGModule (module) = false then
      return Error ("Argument is not a module.");
   fi;

   dim := DimensionFlag (module);
   factors := [];
   for i in [1..dim] do
      factors[i] := [];
   od;
   #factors[i] will contain a list [f1, f2, ..., fr] of the composition factors
   #of module of dimension i. Each fi will have the form [m, n], where m is
   #the module, and n its multiplicity.

   queue := [module];
   #queue is the list of modules awaiting processing.

   while Length (queue) > 0 do
      l := Length (queue);
      cmod := queue[l];
      Unbind (queue[l]);
      InfoComposition ("#I Length of queue = ", l, ", dimension = ", 
                 DimensionFlag (cmod), ".\n");

      if IsIrreducible (cmod) then
         InfoComposition ("#I Irreducible: ");
         #module is irreducible. See if it is already on the list.
         d := DimensionFlag (cmod);
         new := true;
         l := Length (factors[d]);
         i := 1;
         while new and i <= l do
            if IsomorphismGModule (factors[d][i][1], cmod) <> false then
               new := false;
               factors[d][i][2] := factors[d][i][2] + 1;
            fi;
            i := i + 1;
         od;
         if new then
            InfoComposition ("#I  new.\n");
            factors[d][l + 1] := [cmod, 1];
         else 
            InfoComposition ("#I  old.\n");
         fi;
      else
         InfoComposition ("#I Reducible.\n");
         #module is reducible. Add sub- and quotient-modules to queue.
         l := Length (queue);
         q := InducedAction (cmod, SubbasisFlag (cmod));
         queue[l + 1] := q[1]; queue[l + 2] := q[2];
      fi;
   od;

   #Now repack the sequence for output.
   l := 0;
   factorsout := [];
   for i in [1..dim] do
      for j in [1..Length (factors[i])] do
         l := l + 1;
         factorsout[l] := factors[i][j];
      od;
   od;

   return factorsout;

end;

###############################################################################
##
#F  Distinguish ( cf, i ) . .  distinguish a composition factor of a module
##
## cf is assumed to be the output of a call to CompositionFactors, 
## and i is the number of one of the cf.
## Distinguish tries to find a group-algebra element for factor[i]
## which gives nullity zero when applied to all other cf.
## Once this is done, it is easy to find submodules containing this
## composition factor.
## 
Distinguish := function ( cf, i )
   local el, genpair, ngens, orig_ngens, mat, matsi, mats, M, 
         dimi, dim, F, fac, sfac, p, q, oldp, found, extdeg, j, k, 
         lcf, lf, x, y, wno, deg, trying, N, fact, R;

   lcf := Length (cf);
   ngens := Length (GeneratorsFlag (cf[1][1]));
   orig_ngens := ngens;
   F := FieldFlag (cf[1][1]);
   R := PolynomialRing (F);
   matsi := GeneratorsFlag (cf[i][1]);
   dimi := DimensionFlag (cf[i][1]);

   # First check that the existing nullspace has dimension 1 
   # over centralising field. 
   GoodElementGModule (cf[i][1]);

   # First see if the existing element is OK
   # Apply the alg. el. of factor i to every other factor and see if the
   # matrix is nonsingular.
   found := true;
   el := AlgElFlag (cf[i][1]);
   fact := AlgElCharPolFacFlag (cf[i][1]);
   for j in [1..lcf] do
      if j <> i and found then
         mats := GeneratorsFlag (cf[j][1]);
         dim := DimensionFlag (cf[j][1]);
         for genpair in el[1] do
            ngens := ngens + 1;
            mats[ngens] := mats[genpair[1]] * mats[genpair[2]];
         od;
         M := NullMat (dim, dim, F);
         for k in [1..ngens] do
            M := M + el[2][k] * mats[k];
         od;
         #Now throw away extra generators of module
         for k in [orig_ngens + 1..ngens] do
            Unbind (mats[k]);
         od;
         ngens := orig_ngens;
         mat := Value (fact, M);
         if RankMat (mat) < dim then
            found := false;
            InfoComposition ("#I Current element failed on factor ", j, "\n");
         fi;
      fi;
   od;

   if found then
      InfoComposition ("#I Current element worked.\n");
      return;
   fi;

   # That didn't work, so we have to try new random elements.
   wno := 0;
   el := []; el[1] := [];
   extdeg := DegreeFieldExtFlag (cf[i][1]);

   while found = false do
      InfoComposition  ("#I Trying new one. ");
      wno := wno + 1;
      # Add a new generator if there are less than 8 or if wno mod 10 = 0.
      if  ngens < 8 or wno mod 10 = 0 then
         x := Random ([1..ngens]);
         y := x;
         while y = x and ngens > 1 do y := Random ([1..ngens]); od;
         Add (el[1], [x, y]);
         ngens := ngens + 1;
         matsi[ngens] := matsi[x] * matsi[y];
      fi;
      # Now take the new random element
      el[2] := [];
      for j in [1..ngens] do el[2][j] := Random (F); od;
      # First evaluate on cf[i][1].
      M := NullMat (dimi, dimi, F);
      for k in [1..ngens] do
         M := M + el[2][k] * matsi[k];
      od;
      p := CharacteristicPolynomial (M);
      p := EmbeddedPolynomial (R, p);
      # That is necessary in case p is defined over a smaller field that F.
      oldp := Copy (p);
      # extract irreducible factors
      deg := 0;
      fac := [];
      trying := true;
      while deg <= extdeg and trying do
         repeat
            deg := deg + 1;
            if deg > extdeg then
               fac := [p];
            else
               fac := FactorsPolDeg (R, p, deg);
               sfac := Set (fac);
            fi;
         until fac <> [];
         lf := Length (fac);
         if trying and deg <= extdeg then
            j := 1;
            while j <= lf and trying do
               mat := Value (fac[j], M);
               N := NullspaceMat (mat);
               if Length (N) = extdeg then
                  trying := false;
                  SetAlgElFlag (cf[i][1], el);
                  SetAlgElMatFlag (cf[i][1], M);
                  SetAlgElCharPolFlag (cf[i][1], oldp);
                  SetAlgElCharPolFacFlag (cf[i][1], fac[j]);
                  SetAlgElNullspaceVecFlag (cf[i][1], N[1]);
               fi;
               j := j + 1;
            od;
         fi;

         if trying then
            for q in fac do
               p := Quotient (R, p, q);
            od;
         fi;
      od;

      # Now see if it works against the other factors of cf
      if trying = false then
         InfoComposition ("#I Found one. ");
         found := true;
         fact := AlgElCharPolFacFlag (cf[i][1]);
         # Apply the alg. el. of factor i to every other factor and 
         # see if the matrix is nonsingular.
         for j in [1..lcf] do
            if j <> i and found then
               mats := GeneratorsFlag (cf[j][1]);
               dim := DimensionFlag (cf[j][1]);
               ngens := orig_ngens;
               for genpair in el[1] do
                  ngens := ngens + 1;
                  mats[ngens] := mats[genpair[1]] * mats[genpair[2]];
               od;
               M := NullMat (dim, dim, F);
               for k in [1..ngens] do
                  M := M + el[2][k] * mats[k];
               od;
               # Now throw away extra generators of module
               for k in [orig_ngens + 1..ngens] do
                  Unbind (mats[k]);
               od;
               mat := Value (fact, M);
               if RankMat (mat) < dim then
                  found := false;
                  InfoComposition ("#I Failed on factor ", j, "\n");
               fi;
            fi;
         od;
      fi;
      if found then
         InfoComposition ("#I It worked!\n");
      fi;
   od;
   # Finally throw away extra generators of s[i][1]
   for k in [orig_ngens + 1..ngens] do
      Unbind (matsi[k]);
   od;

end;

###############################################################################
##
#F  MinimalSubGModule ( module, cf, i ) . .  find minimal submodule containing 
##                                                 a given composition factor.
##
## cf is assumed to be the output of a call to CompositionFactors (module), 
## and i is the number of one of the cf.
## It is assumed that Distinguish (cf, i) has already been called.
## A basis of a minimal submodule of module containing the composition factor
## cf[i][1] is calculated and returned - i.e. if cf[i][2] = 1.
##
MinimalSubGModule := function ( module, cf, i )
   local el, genpair, ngens, orig_ngens, mat, mats, M, dim, F, 
         j, k, N, fact;

   if IsGModule (module) = false then
      return Error ("First argument is not a module.");
   fi;

   ngens := Length (GeneratorsFlag (module));
   orig_ngens := ngens;
   F := FieldFlag (module);

   #Apply the alg. el. of factor i to module
   el := AlgElFlag (cf[i][1]);
   mats := GeneratorsFlag (module);
   dim := DimensionFlag (module);
   for genpair in el[1] do
      ngens := ngens + 1;
      mats[ngens] := mats[genpair[1]] * mats[genpair[2]];
   od;
   M := NullMat (dim, dim, F);
   for k in [1..ngens] do
      M := M + el[2][k] * mats[k];
   od;
   #Now throw away extra generators of module
   for k in [orig_ngens + 1..ngens] do
      Unbind (mats[k]);
   od;
   ngens := orig_ngens;
   fact := AlgElCharPolFacFlag (cf[i][1]);
   mat := Value (fact, M);
   N := NullspaceMat (mat);
   return (SpinBasis (N[1], mats, ngens));

end;
