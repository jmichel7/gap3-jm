#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @ (#)$Id: homomorp.g,v 1.1 1997/03/10 13:52:31 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: homomorp.g,v $
#H  Revision 1.1  1997/03/10 13:52:31  gap
#H  VERSION 1.0
#H
#H  Revision 1.2  1997/01/05 10:49:27  fceller
#H  added Eamonn's new version to the reprository
#H
#H  Revision 1.1  1996/12/25 09:03:45  fceller
#H  changed long filenames to MS-DOS conform filenames,
#H  the init files are *NOT* yet updated
#H
#H  Revision 1.2  1996/12/22 07:48:09  fceller
#H  new smash version
#H
#H  Revision 1.1  1996/11/28 13:14:48  fceller
#H  added "smash" and "reducible" to the repository
#H
##
#homomorphism.g
#
############################################################################
##
#F  InfoHomGModule (...)  . . . . . . . . . . . for debugging assistance
##
##
if not IsBound (InfoHomGModule)  then InfoHomGModule:= Ignore;  fi;
#############################################################################
##
#F  IsomomorphismGModule (module1, module2) . . . . 
##  decide whether two irreducible modules are isomorphic.
## 
## If the 2 modules are not isomorphic, this function returns false;
## if they are isomorphic it returns the matrix B, whose rows form the 
## basis of module2  which is the image of the standard basis for module1.
## Thus if X and Y are corresponding matrices in the generating sets
## for module1 and module2 respectively, Y = BXB^-1
## It is assumed that the same group acts on both modules.
## Otherwise who knows what will happen?
## 
IsomorphismGModule := function (module1, module2)
   local matrices, matrices1, matrices2, F, R, dim, swapmodule, genpair,
         swapped, orig_ngens, i, j, el, p, fac, ngens, M, mat, v1, v2, v, 
         N, basis, basis1, basis2;

   if IsGModule (module1) = false then 
      return Error ("Argument is not a module.");
   elif IsGModule (module2) = false then 
      return Error ("Argument is not a module.");
   elif FieldFlag (module1) <> FieldFlag (module2) then 
      return Error ("GModules are defined over different fields.");
   fi;

   swapped := false;
   if ReducibleFlag (module1) <> false then
      if ReducibleFlag (module2) <> false then
         return Error ("Neither module is known to be irreducible.");
      else
         # The second module is known to be irreducible, so swap arguments.
         swapmodule := module2; module2 := module1; module1 := swapmodule;
         swapped := true;
         InfoHomGModule ("#I Second module is irreducible. Swap them round.\n");
      fi;
   fi;

   #At this stage, module1 is known to be irreducible
   dim := DimensionFlag (module1);
   if dim <> DimensionFlag (module2) then
      InfoHomGModule ("#I GModules have different dimensions.\n");
      return false;
   fi;
   F := FieldFlag (module1);
   R := PolynomialRing (F);

   #First we must check that our nullspace is 1-dimensional over the
   #centralizing field.

   InfoHomGModule ("#I Checking nullspace 1-dimensional over centralising field.\n");
   GoodElementGModule (module1);
   matrices1 := GeneratorsFlag (module1); matrices2 := GeneratorsFlag (module2);
   ngens := Length (matrices1);
   orig_ngens := ngens;
   if ngens <> Length (matrices2) then
      return Error ("GModules have different numbers of defining matrices.");
   fi;

   # Now we calculate the element in the group algebra of module2 that 
   # corresponds to that in module1. This is done using the AlgEl flag 
   # for module1. We first extend the generating set in the same way as 
   # we did for module1, and then calculate the group alg. element as 
   # a linear sum in the generators.

   InfoHomGModule ("#I Extending generating set for second module.\n");
   el := AlgElFlag (module1);
   for genpair in el[1] do
      ngens := ngens + 1;
      matrices2[ngens] := matrices2[genpair[1]] * matrices2[genpair[2]];
   od;
   M := NullMat (dim, dim, F);
   for i in [1..ngens] do
      M := M + el[2][i] * matrices2[i];
   od;
   # Having done that, we no longer want the extra generators of module2, 
   # so we throw them away again.
   for i in [orig_ngens + 1..ngens] do
      Unbind (matrices2[i]);
   od;

   InfoHomGModule ("#I Calculating characteristic polynomial for second module.\n");
   p := CharacteristicPolynomial (M);
   p := EmbeddedPolynomial (R, p);
   if p <> AlgElCharPolFlag (module1) then
      InfoHomGModule ("#I Characteristic polynomial different.\n");
      return false;
   fi;
   fac := AlgElCharPolFacFlag (module1);
   mat := Value (fac, M);
   InfoHomGModule ("#I Calculating nullspace for second module.\n");
   N := NullspaceMat (mat);
   if Length (N) <> AlgElNullspaceDimensionFlag (module1) then
      InfoHomGModule ("#I Null space dimensions different.\n");
      return false;
   fi;

   # That concludes the easy tests for nonisomorphism. Now we must proceed
   # to spin up. We first form the direct sum of the generating matrices.
   InfoHomGModule ("#I Spinning up in direct sum.\n");
   matrices := MatrixSum (matrices1, matrices2);
   v1 := AlgElNullspaceVecFlag (module1);
   v2 := N[1];
   v := Concatenation (v1, v2);
   basis := SpinBasis (v, matrices);
   if Length (basis) = dim then
      basis1 := []; basis2 := [];
      for i in [1..dim] do
         basis1[i] := []; basis2[i] := [];
         for j in [1..dim] do
            basis1[i][j] := basis[i][j];
            basis2[i][j] := basis[i][j + dim];
         od;
      od;
      if swapped then
         return basis2^-1 * basis1;
      else
         return basis1^-1 * basis2;
      fi;
   else
      return false;
   fi;
end;

#############################################################################
##
#F  MatrixSum (matrices1, matrices2)  . . direct sum of two lists of matrices
##
MatrixSum := function (matrices1, matrices2) 
   local dim, dim1, dim2, matrices, zero, nmats, i, j, k;
   dim1 := Length (matrices1[1]); dim2 := Length (matrices2[1]);
   dim := dim1 + dim2;
   zero:= 0 * matrices1[1][1][1];
   matrices := [];
   nmats := Length (matrices1);
   for i in [1..nmats] do
      matrices[i] := NullMat (dim, dim, zero);
      for j in [1..dim1] do for k in [1..dim1] do
         matrices[i][j][k] := matrices1[i][j][k];
      od; od;
      for j in [1..dim2] do for k in [1..dim2] do
         matrices[i][j + dim1][k + dim1] := matrices2[i][j][k];
      od; od;
   od;

   return  matrices;
end;

#############################################################################
##
#F  HomGModule ( m1, m2) . . . . homomorphisms from an irreducible
##                         . . . GModule to an arbitrary GModule
##
## It is assumed that m1 is a module that has been proved irreducible
##  (using IsIrreducible), and m2 is an arbitrary module for the same group.
## A basis of the space of G-homomorphisms from m1 to m2 is returned.
## Each homomorphism is given as a list of base images.
## 
HomGModule := function (m1, m2)

   local F, ngens, orig_ngens, mats1, mats2, dim1, dim2, m1bas, imbases, 
         el, genpair, fac, mat, N, imlen, subdim, leadpos, vec, imvecs, 
         numrels, rels, leadposrels, newrels, bno, genno, colno, rowno, 
         zero, looking, ans, i, j, k;

   if IsGModule (m1) = false then 
      return Error ("First argument is not a module.");
   elif ReducibleFlag (m1) <> false then
      return Error ("First module is not known to be irreducible.");
   fi;

   if IsGModule (m2) = false then 
      return Error ("Second argument is not a module.");
   fi;

   mats1 := GeneratorsFlag (m1); mats2 := GeneratorsFlag (m2);
   ngens := Length (mats1);
   if ngens <> Length (mats2) then
      return Error ("GModules have different numbers of generators.");
   fi;

   F := FieldFlag (m1);
   if F <> FieldFlag (m2) then
      return Error ("GModules are defined over different fields.");
   fi;
   zero := Zero (F);

   dim1 := DimensionFlag (m1); dim2 := DimensionFlag (m2);

   m1bas := [];
   m1bas[1] :=  AlgElNullspaceVecFlag  (m1);

   # In any homomorphism from m1 to m2, the vector in the nullspace of the
   # algebraic element that was used to prove irreducibility  (which is now
   # m1bas[1]) must map onto a vector in the nullspace of the same algebraic
   # element evaluated in m2. We therefore calculate this nullspaces, and
   # store a basis in imbases.

   InfoHomGModule ("#I Extending generating set for second module.\n");
   orig_ngens := ngens;
   el := AlgElFlag (m1);
   for genpair in el[1] do
      ngens := ngens + 1;
      mats2[ngens] := mats2[genpair[1]] * mats2[genpair[2]];
   od;
   mat := NullMat (dim2, dim2, F);
   for i in [1..ngens] do
      mat := mat + el[2][i] * mats2[i];
   od;
   # Having done that, we no longer want the extra generators of m2, 
   # so we throw them away again.
   for i in [orig_ngens + 1..ngens] do
      Unbind (mats2[i]);
   od;
   ngens := orig_ngens;

   fac := AlgElCharPolFacFlag (m1);
   mat := Value (fac, mat);
   InfoHomGModule ("#I Calculating nullspace for second module.\n");
   N := NullspaceMat (mat);
   imlen := Length (N);
   InfoHomGModule ("#I Dimension = ", imlen, ".\n");
   if imlen = 0 then
      return [];
   fi;

   imbases := [];
   for i in [1..imlen] do
      imbases[i] := [N[i]];
   od;

   # Now the main algorithm starts. We are going to spin the vectors in m1bas
   # under the action of the module generators, norming as we go. Every
   # operation that we perform on m1bas will also be performed on each of the 
   # vectors in  imbas[1], ..., imbas[imlen].
   # When we find a vector that norms to zero in m1bas, then the image of this
   # under a homomorphism must be zero. This leads to a linear relation
   # amongst some vectors in imbas. We store up such relations, echelonizing as
   # we go. At the end, if we have numrels subch independent relations, then
   # there will be imlen - numrels independent homomorphisms from m1 to m2, 
   # which we can then calculate.

   subdim := 1; # the dimension of module spanned by m1bas
   numrels := 0;
   rels := [];

   #leadpos[j] will be the position of the first nonzero entry in m1bas[j]
   leadpos := [];
   vec := m1bas[1];
   j := 1;
   while j <= dim1 and vec[j] = zero do j := j + 1; od;
   leadpos[1] := j;
   k := vec[j]^-1;
   m1bas[1] := k * vec;
   for i in [1..imlen] do
      imbases[i][1] := k * imbases[i][1];
   od;

   leadposrels := [];
   #This will play the same role as leadpos but for the relation matrix.
   InfoHomGModule ("#I Starting spinning.\n");
   bno := 1;
   while bno <= subdim do
      for genno in [1..ngens] do
         # apply generator no. genno to submodule generator bno
         vec := m1bas[bno] * mats1[genno];
         # and do the same to the images
         imvecs := [];
         for i in [1..imlen] do
            imvecs[i] := imbases[i][bno] * mats2[genno];
         od;
         # try to express w in terms of existing submodule generators
         # make same changes to images
         j := 1;
         for  j in [1..subdim] do
            k := vec[leadpos[j]];
            if k <> zero then
               vec := vec - k * m1bas[j];
               for i in [1..imlen] do
                  imvecs[i] := imvecs[i] - k * imbases[i][j];
               od;
            fi;
         od;

         j := 1;
         while j <= dim1 and vec[j] = zero do j := j + 1; od;
         if j <= dim1 then
            #we have found a new generator of the submodule
            subdim := subdim + 1;
            leadpos[subdim] := j;
            k := vec[j]^-1;
            m1bas[subdim] := k * vec;
            for i in [1..imlen] do
               imbases[i][subdim] := k * imvecs[i];
            od;
         else
            # vec has reduced to zero. We get relations among the imvecs.
            # (these are given by the transpose of imvec)
            # reduce these against any existing relations.
            newrels := TransposedMat (imvecs);
            for i in [1..Length (newrels)] do
               vec := newrels[i];
               for j in [1..numrels] do
                  k := vec[leadposrels[j]];
                  if k <> zero then
                     vec := vec - k * rels[j];
                  fi;
               od;
               j := 1;
               while j <= imlen and vec[j] = zero do j := j + 1; od;
               if j <= imlen then
                  # we have a new relation
                  numrels := numrels + 1;
                  # if we have imlen relations, there can be no homomorphisms
                  # so we might as well give up immediately
                  if numrels = imlen then
                     return [];
                  fi;
                  k := vec[j]^-1;
                  rels[numrels] := k * vec; 
                  leadposrels[numrels] := j;
               fi;
            od;
         fi;
      od;
      bno := bno + 1;
   od;

   # That concludes the spinning. Now we do row operations on the im1bas to
   # make it the identity, and do the same operations to the imvecs.
   # Then the homomorphisms we output will be the basis images.
   InfoHomGModule ("#I Done. Reducing spun up basis.\n");

   for colno in [1..dim1] do
      rowno := colno;
      looking := true;
      while rowno <= dim1 and looking do
         if m1bas[rowno][colno] <> zero then
            looking := false;
            if rowno <> colno then
               #swap rows rowno and colno
               vec := m1bas[rowno]; m1bas[rowno] := m1bas[colno]; 
               m1bas[colno] := vec;
               #and of course the same in the images
               for i in [1..imlen] do
                  vec := imbases[i][rowno];
                  imbases[i][rowno] := imbases[i][colno]; 
                  imbases[i][colno] := vec;
               od;
            fi;
            # and then clear remainder of column
            for j in [1..dim1] do
               if j <> colno and m1bas[j][colno] <> zero then
                  k := m1bas[j][colno];
                  m1bas[j] := m1bas[j] - k * m1bas[colno];
                  for i in [1..imlen] do
                     imbases[i][j] := imbases[i][j] - k * imbases[i][colno];
                  od;
               fi;
            od;
         fi;
         rowno := rowno + 1;
      od;
   od;

   #Now we are ready to compute and output the linearly independent 
   #homomorphisms.  The coefficients for the solution are given by 
   #the basis elements of the nullspace of the transpose of rels.

   InfoHomGModule ("#I Done. Calculating homomorphisms.\n");
   if rels = [] then
      rels := NullMat (imlen, 1, F);
   else
      rels := TransposedMat (rels);
   fi;
   N := NullspaceMat (rels);
   ans := [];
   for k in [1..Length (N)] do
      vec := N[k];
      mat := NullMat (dim1, dim2, F);
      for i in [1..imlen] do
         mat := mat + vec[i] * imbases[i];
      od;
      ans[k] := mat;
   od;

   return ans;
end;

#############################################################################
##
#F  SortHomGModule ( m1, m2, homs)  . . 
##  sort output of HomGModule according to their images
##
## It is assumed that m1 is a module that has been proved irreducible
## (using IsIrreducible), and m2 is an arbitrary module for the same group, 
## and that homs is the output of a call HomGModule (m1, m2).
## Let e be the degree of the centralising field of m1.
## If e = 1 then SortHomGModule does nothing. If e > 1, then it replaces 
## the basis contained in homs by a new basis arranged in the form
## b11, b12, ..., b1e, b21, b22, ...b2e, ..., br1, br2, ...bre,  where each
## block of  e  adjacent basis vectors are all equivalent under the
## centralising field of m1, and so they all have the same image in  m2.
## A complete list of the distinct images can then be obtained with a call
## to DistinctIms (m1, m2, homs).
## 
SortHomGModule := function (m1, m2, homs)
   local e, F, ngens, mats1, mats2, dim1, dim2, centmat, fullimbas, oldhoms, 
         hom, homno, dimhoms, newdim, subdim, leadpos, vec, nexthom, 
         i, j, k, zero;

   if IsAbsolutelyIrreducible (m1) then return; fi;

   e := DegreeFieldExtFlag (m1);
   F := FieldFlag (m1);
   zero := Zero (F);

   mats1 := GeneratorsFlag (m1);  mats2 := GeneratorsFlag (m2);
   dim1 := DimensionFlag (m1);  dim2 := DimensionFlag (m2);
   ngens := Length (mats1);
   centmat := CentMatFlag (m1);

   fullimbas := [];
   subdim := 0;
   leadpos := [];

   # fullimbas will contain an echelonised basis for the submodule of m2
   # generated by all images of the basis vectors of hom that we have found
   # so far; subdim is its length.

   # We go through the existing basis of homs. 
   # For each hom in the basis, we first check whether the first vector in 
   # the image  of hom is in the space spanned by fullimbas. 
   # If so, we reject hom. If not, then hom is adjoined to the new
   # basis of homs, as are the other e-1 linearly independent homomorphisms
   # that are equivalent to hom by a multiplication by centmat. The
   # resulting block of e homomorphisms all have the same image in m2.

   # first make a copy of homs.

   oldhoms := Copy (homs);
   dimhoms := Length (homs);

   homno := 0; newdim := 0;

   while homno < dimhoms and newdim < dimhoms do
      homno := homno + 1;
      nexthom := oldhoms[homno];
      vec := nexthom[1];

      #Now check whether vec is in existing submodule spanned by fullimbas   
      j := 1;
      for j in [1..subdim] do
         k := vec[leadpos[j]];
         if k <> zero then
            vec := vec - k * fullimbas[j];
         fi;
      od;

      j := 1;
      while j <= dim2 and vec[j] = zero do j := j + 1; od;

      if j <= dim2 then
         #vec is not in the image, so we adjoin this homomorphism to the list;
         #first adjoin vec and all other basis vectors in the image to fullimbas
         subdim := subdim + 1;
         leadpos[subdim] := j;
         k := vec[j]^-1;
         fullimbas[subdim] := k * vec;
         for i in [2..dim1] do
            vec := nexthom[i];
            j := 1;
            for  j in [1..subdim] do
               k := vec[leadpos[j]];
               if k <> zero then
                  vec := vec - k * fullimbas[j];
               fi;
            od;

            j := 1;
            while j <= dim2 and vec[j] = zero do j := j + 1; od;
            subdim := subdim + 1;
            leadpos[subdim] := j;
            k := vec[j]^-1;
            fullimbas[subdim] := k * vec;
         od;

         newdim := newdim + 1;
         homs[newdim] := nexthom;

         #Now add on the other e - 1 homomorphisms equivalent to 
         #newhom by centmat.
         for k in [1..e - 1] do
            nexthom := centmat * nexthom;
            newdim := newdim + 1;
            homs[newdim] := nexthom;
         od;
      fi;
   od;

end;

#############################################################################
##
#F MinimalSubGModules (m1, m2, [max]) . . 
## minimal submodules of m2 isomorphic to m1
##
## It is assumed that m1 is a module that has been proved irreducible
##  (using IsIrreducible), and m2 is an arbitrary module for the same group.
## MinimalSubGModules computes and outputs a list of normed bases for all of 
## the distinct minimal submodules of m2 that are isomorphic to m1.
## max is an optional maximal number - if the total number of submodules
## exceeds max, then the procedure aborts.
## First HomGModule is called and then SortHomGModule to get a basis for
## the homomorphisms from m1 to m2 in the correct order.
## It is then easy to write down the list of distinct images.
## 
MinimalSubGModules := function (arg)

   local m1, m2, max, e, homs, coeff,  dimhom, edimhom, F, elF, q, 
         submodules, sub, adno, more, count, sr, er, i, j, k ;

   if Number (arg) < 2 or Number (arg) > 3 then
      Error ("Number of arguments to MinimalSubGModules must be 2 or 3.");
   fi;

   m1 := arg[1]; m2 := arg[2];
   if Number (arg) = 2 then max := 0; else max := arg[3]; fi;

   InfoHomGModule ("#I Calculating homomorphisms from m1 to m2.\n");
   homs := HomGModule (m1, m2);
   InfoHomGModule ("#I Sorting them.\n");
   SortHomGModule (m1, m2, homs);

   F := FieldFlag (m1);
   e := DegreeFieldExtFlag (m1);
   dimhom := Length (homs);
   edimhom := dimhom / e;
   submodules := [];
   count := 0;
   coeff := [];
   elF := Elements (F);
   q := Length (elF);
   for i in [1..dimhom] do coeff[i] := 1; od;

   #coeff[i] will be an integer in the range [1..q] corresponding to the
   #field element elF[coeff[i]].
   #Each submodule will be calculated as the image of the homomorphism
   #elF[coeff[1]] * homs[1] +...+  elF[coeff[dimhom]] * homs[dimhom]
   #for appropriate field elements elF[coeff[i]]. 
   #We get each distinct submodule
   #exactly once by making the first nonzero elF[coeff[i]] to be 1, 
   #and all other elF[coeff[i]]'s in that block equal to zero.

   InfoHomGModule ("#I Done. Calculating submodules.\n");

   for i in Reversed ([1..edimhom]) do
      j := e * (i - 1) + 1;
      coeff[j] := 2;  #giving field element 1.
      for k in [j + 1..dimhom] do coeff[k] := 1; od; # field element 0.
      sr := j + e; er := dimhom;
      #coeff[i] for i in [sr..er] ranges over all field elements.

      more := true;
      adno := er;
      while more do
         count := count + 1;
         if max > 0 and count > max then
            InfoHomGModule ("#I Number of submodules exceeds ", max, ". Aborting.\n");
            return submodules;
         fi;

         # Calculate the next submodule
         sub := Copy (homs[j]);
         for k in [sr..er] do
            sub := sub + elF[coeff[k]] * homs[k];
         od;
         TriangulizeMat (sub);
         Add (submodules, sub);

         #Move on to next set of coefficients if any
         while adno >= sr and coeff[adno]=q do
            coeff[adno] := 1;
            adno := adno - 1;
         od;
         if adno < sr then
            more := false;
         else
            coeff[adno] := coeff[adno] + 1;
            adno := er;
         fi;
      od;

   od;

   return submodules;

end;
