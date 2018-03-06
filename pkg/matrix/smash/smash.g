#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @ (#)$Id: smash.g,v 1.1 1997/03/10 13:52:40 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: smash.g,v $
#H  Revision 1.1  1997/03/10 13:52:40  gap
#H  VERSION 1.0
#H
#H  Revision 1.4  1997/01/05 10:49:36  fceller
#H  added Eamonn's new version to the reprository
#H
#H  Revision 1.3  1996/12/22 07:48:16  fceller
#H  new smash version
#H
#H  Revision 1.2  1996/12/12 10:51:31  fceller
#H  new version by Eamonn
#H
#H  Revision 1.1  1996/11/28 13:14:58  fceller
#H  added "smash" and "reducible" to the repository
#H
##
############################################################################
##
#F  InfoSmash  (...)  . . . . . . . . . . . for debugging assistance
##
##
if not IsBound (InfoSmash)  then InfoSmash := Ignore;  fi;

#############################################################################
##
#F  SmashGModule ( module, S[, flag])  . . 
##
## S is a set of matrices, SmashGModule attempts to find some way of decomposing
## the module with respect to S. It return true if some breakdown is found, 
## false otherwise.
##
SmashGModule := function ( arg )

   local module, S, PartialSmash, matrices, F, d, ranSub, Smodule, irreds, 
         submodules, summands, t, blocks, W, Wmodule, dimW, AllScalar, one,
         A, B, i, m, mm, g, c, C, e, s, basis;

   module := arg[1];
   if IsGModule (module) = false then
      return Error ("First argument must be a GModule");
   fi;

   S := arg[2];
   if IsList (S) = false or ForAll (S, x -> IsMat (x)) = false then 
      return Error ("Second argument must be a list of matrices");
   fi;

   # When PartialSmash is true, we terminate as soon as we discover
   # that the normal subgroup must act absolutely irreducibly 
   # -- hence, we do not seek to decide whether or not the group is 
   # is a normaliser of a p-group  or there is symmetric tensor 
   # product decomposition 

   if Number (arg) = 3 and arg[3] = "PartialSmash" then 
      PartialSmash := true; 
   else 
      PartialSmash := false; 
   fi;

   InfoSmash ("#I Starting call to SmashGModule\n");

   if S = [] then
      Print ("#I Set S is empty.\n");
      return false;
   elif ForAll (S, x -> IsScalar (x))  then
      InfoSmash ("#I S contains only scalar matrices.\n");
      return false; 
   fi;

   if IsIrreducible (module) = false then 
      InfoSmash ("#I GModule is not irreducible.\n");
      return false;
   fi;

   if IsAbsolutelyIrreducible (module) = false then 
      InfoSmash ("#I GModule is not absolutely irreducible.\n");
      return false;
   fi;

   F := FieldFlag (module);
   one := One (F);
   d := DimensionFlag (module);

   while 0 <> 1 do
      InfoSmash ("#I At top of main SmashGModule loop, S has ", Length (S), 
              " elements.\n");
      Smodule := GModule (S, F);
      matrices := GeneratorsFlag (module);
      ranSub := RandomIrreducibleSubGModule (Smodule);

      if ranSub <> false then  
         # <S> acts reducibly.
         W := ranSub[1];
         Wmodule := ranSub[2];

         irreds := SemiSimpleDecomposition (matrices, S, W, Wmodule, F);

         #S may well have changed during this function call, whatever the 
         #outcome.  If irreds is false then we have to go back up to the 
         #top of the loop and reselect a random irreducible submodule with 
         #our new enlarged S.  If irreds isn't false, then it's a pair 
         #of lists. The first list contains a set of bases for irreducible 
         #submodules, the second list the irreducible submodules themselves.

         if irreds <> false then 

            #the module decomposes as a sum of irreducibles W1, .., Wt 
            #under S, where t > 1.
            summands := irreds[1];
            submodules := irreds[2];
            t := Length (summands);

            InfoSmash ("#I GModule is decomposed as a sum of ", t, 
                    " irreducibles of dimension ", Length (W), "\n");
            dimW := Length (W);
            blocks := MinBlocks (module, W);
            if NumberBlocksFlag (blocks) > 1 then
               SetImprimitiveFlag (module, true);
               SetBlockSystemFlag (module, blocks);
               InfoSmash ("#I Found a block system.\n");
               InfoSmash ("#I SmashGModule returns true.\n");
               return true;
            fi;

            IsAbsolutelyIrreducible (Wmodule);
            basis := W;
            i := 2;
            while i <= t do
               B := IsomorphismGModule(Wmodule, submodules[i]);
               # Testing for isomorphism.
               if B <> false then 
                  basis := Concatenation (basis, B * summands[i]);
                  i := i + 1;
               else 
                  i := t + 1;
               fi;
            od;

            if B = false then 
               #the modules are not isomorphic.
               #That must mean we don't have the right decomposition, 
               #or we would have found blocks of imprimitivity. 
               #So we enlarge S and loop around again. 

               InfoSmash ("#I SubGModules, non-isomorphic.\n");
               InfoSmash ("#I Adding random translating conjugate to S.\n");
               AddRandomTranslatingConjugate (module, S, W, F);

            elif AbsolutelyReducibleFlag (Wmodule) = false then 

               #If the modules are isomorphic and absolutely irreducible,
               #we look for a tensor decomposition. If we fail to find 
               #one it must be because we don't have the right decomposition. 
               #So we enlarge S and loop around again.

               if TensorProductDecomposition (module, basis, t, dimW) then
                  InfoSmash ("#I Found a tensor product decomposition.\n");
                  InfoSmash ("#I SmashGModule returns true.\n");
                  return true;
               else
                  InfoSmash ("#I Failed to find tensor decomposition.\n");
                  InfoSmash ("#I Adding random translating conjugate to S.\n");
                  AddRandomTranslatingConjugate (module, S, W, F);
               fi;
            else 
               #we have a centralising field of degree e
               #We hope to find the group inside GammaL (d/e, q^e). 
               #If we fail to find the embedding, it must be because we 
               #have the wrong centraliser. So we enlarge S by something 
               #which doesn't commute with the matrix C which we thought 
               #was the centraliser, and loop around again.

               c := CentMatFlag (Wmodule);  
               #c is now the centraliser of Wmodule
               #and is a dimW x dimW matrix, wrt the standard basis.
               #We find the centraliser of the full S-module using the 
               #isomorphisms between the submodules. 
               #First set all the entries to be zero

               C := basis^-1 * 
                    KroneckerProduct (IdentityMat (t, one), c) * basis;

               # The Kronecker product gives us our centralising matrix with 
               # respect to the basis which is the concatenation of 
               # summands[1], B[1]summands[2], .. , B[t-1]summands[t]
               # We check that it commutes with everything in S.

               for s in S do
                  if s * C <> C * s then 
                     Error ("ERROR: centraliser isn't right.\n");
                  else 
                     InfoSmash ("#I Centraliser commutes ok.\n");
                  fi;
               od;

               e := DegreeFieldExtFlag (Wmodule);
               if SemiLinearDecomposition  (module, S, C, e) = true then
                  InfoSmash ("#I Group embeds in GammaL (", d/e, ", ", 
                          F, "^", e, ").\n");
                  InfoSmash ("#I SmashGModule returns true.\n");
                  return true;
               fi;
            fi;
         fi;
      else 
         # <S> acts irreducibly on the module.
         # So either <S> acts absolutely irreducibly or 
         # <S> <= GL (d/e, q^e) for some E, and 
         # G <= GammaL (d/e, q^e), for some e > 1. 
         repeat
            if IsAbsolutelyIrreducible (Smodule) = true then
               InfoSmash ("#I S acts absolutely irreducibly on the module.\n");
               if PartialSmash = true then return false; fi;

               if ExtraSpecialDecomposition (module, S) = true then
                  InfoSmash ("#I G normalises a p-group.\n");
                  InfoSmash ("#I SmashGModule returns true.\n");
                  return true;
               elif SymTensorProductDecomposition (module, Smodule) = true then
                  InfoSmash ("#I G preserves a symmetric tensor product.\n");
                  InfoSmash ("#I SmashGModule returns true.\n");
                  return true;
               else
                  InfoSmash ("#I SmashGModule returns false.\n");
                  return false;
               fi;
            else
               C := CentMatFlag (Smodule);
               e := DegreeFieldExtFlag (Smodule);
               if SemiLinearDecomposition (module, S, C, e) = true then
                  InfoSmash ("#I Group embeds in GammaL (", d/e, ", ", 
                          F, "^", e, ").\n");
                  InfoSmash ("#I SmashGModule returns true.\n");
                  return true;
               else
                  # When SemiLinearDecomposition returns false, 
                  # it must be the case that C doesn't centralise <S>^G. 
                  # So S is enlarged by the addition of an element which 
                  # doesn't centralise C. Then IsAbsolutelyIrreducible will 
                  # be reapplied. 
                  # We don't want to recompute the whole module and retest 
                  # irreducibility, so we change some flags on the existing
                  # module, those corresponding to the set of matrices that 
                  # act on it, and those set during IsAbsolutelyIrreducible. 
                  InfoSmash ("#I Group doesn't embed in GammaL (", d/e, 
                          ", ", F, "^", e, ").\n");
                  UndoAbsolutelyIrreducibleFlags (Smodule);
                  mm := Length (GeneratorsFlag (Smodule));
                  m := Length (S) - mm;
                  for i in [1..m] do
                     EnlargeIrreducibleGModule (Smodule, S[mm + i]);
                  od;
               fi;
            fi;
         until 0 = 1;
      fi;

   od;  

end;

#############################################################################
##
#F SemiSimpleDecomposition (matrices, S, W, Wmodule, F) . . . 
##
## S generates a subgroup of G = <matrices>, which acts on a vector space 
## over F.  W  spans an irreducible S-module, Wmodule. We hope to
## decompose the underlying vector space as a direct sum of 
## irreducible <S'>-modules
## which are translates of W under G, where S <= S' <= S^G. If we can do this, 
## the function returns a pair of lists. The first list contains
## a set of bases for irreducible submodules, the second list the irreducible 
## submodules themselves. At the same time, S is replaced by S'
## Otherwise the function returns false, and S is replaced by S' 
## which no longer fixes the subspace W.
## S may change during this function, whatever the outcome.
## If it does, and the function return true, Wmodule will have the extra 
## matrices added to it.

SemiSimpleDecomposition := function ( matrices, S, W, Wmodule, F)
   local s, translates, M;

   translates := TranslatesDirectSum (matrices, W, F);

   if translates = false then 
      InfoSmash ("#I Adding random conjugate to S that translates W.\n");
      M := GModule (matrices, F);  
      AddRandomTranslatingConjugate (M, S, W, F);
      return false;
   else 
      s := TranslatesSGModules (matrices, S, translates, F);
      if s <> true then 
         InfoSmash ("#I Translates of W aren't modules.\n");
         Add (S, s); 
         return false;
      else 
         return TranslatesIrreducible (matrices, S, Wmodule, translates, F);
      fi;
   fi;

end;

#############################################################################
##
#F TranslatesDirectSum (matrices, W, F) . . .
##
## W spans a subspace of a vector space over F on which G = <matrices> acts.
## We attempt to find the space as a direct sum of translates of <W>. 
## If we succeed we return a pair of lists, [summands, info].
## Summands is a list of bases of translates of W which span the space as
## a direct sum. The second list identifies the translation; 
## if info[k]=[i, j] then summands[k] spans a translate of an earlier 
## summand summand[i] by the element matrices[j].
## If we fail we return false. 
## 
TranslatesDirectSum := function ( matrices, W, F)
   local socle, newsocle, dimsoc, summands, info, i, j, Wi, dimW, g, Wig;

   socle := W;
   dimW := Length (W);
   dimsoc := dimW;
   summands := [W];
   info := [[]]; # the first entry is meaningless

   i := 1;
   while i <= Length (summands) do
      Wi := summands[i];
      j := 1;
      while j <= Length (matrices) do
         g := matrices[j];

         # If Wig turns out to span an submodule, it's essential 
         # (for the SubGModulefunction) that it's a normed basis. 
         # The function Base ensures that it is. 
         # We might as well call it here, since it'll make all the 
         # RowSpace function calls below run faster, even if
         # it turns out not to span a submodule.  

         Wig := Base (RowSpace (Wi * g, F)); 
         newsocle := Base (RowSpace (Concatenation (socle, Wig), F));

         # the dimension should be dimsoc or dimsoc + dimW. 
         # If it isn't, we return false.
         if Length (newsocle) = dimsoc + dimW then
            socle := newsocle;
            dimsoc := Length (newsocle);
            Add (summands, Wig);
            Add (info, [i, j]);
            if dimsoc = Length (matrices[1]) then return [summands, info]; fi;
         elif Length (newsocle) <> dimsoc then
            return false;
         fi;

         j := j + 1;
      od;
      i := i + 1;
   od;

   return false;

end;

############################################################################
##
#F TranslatesSGModules (matrices, S, translates, F) . . .
##
## G = <matrices> and S, a subgroup of G, act on the same space.
## translates is a pair [summands, info] as returned by TranslatesDirectSum, 
## so summands is a list of normed bases, which together span the space.
## The first basis spans a module W for S, the remaining are
## translates of that and the list of pairs info identifies this translation.
## We check to see if the translates are also submodules of Smodule. If not
## then some conjugate of S, which we can identify, does not fix W, 
## and we return that conjugate. Otherwise we return true.
## 
TranslatesSGModules := function (matrices, S, translates, F)

   local summands, info, k, Wk, dimW, s, g; 

   summands := translates[1];
   info := translates[2]; 
   k := 2;
   dimW := Length (summands[1]);
   while k <= Length (summands) do
      Wk := summands[k];
      for s in S do
         if Length (Base (RowSpace (Concatenation (Wk, Wk * s), F))) > dimW then
            g := matrices[info[k][2]];
            k := info[k][1];
            # SR - it's not clear to me that we really need to get the g that 
            # translates W to Wk. We might do as well to use simply the g that 
            # translated a previous Wi to give Wk. Then gsg^-1 doesn't fix Wi, 
            # and is not yet in S. 

            while k <> 1 do 
               g := matrices[info[k][2]] * g; 
               k := info[k][1];
            od;
            return g * s * g^-1;
         fi;
      od; 
      k := k + 1;
   od;
   return true;

end;

###############################################################################
##
#F TranslatesIrreducible (matrices, S, Wmodule, translates, F) . . .
##
## translates is a pair [summands, info] as returned by TranslatesDirectSum, 
## such that each of the subspaces in the list summands has already
## been shown to be a module for S.
## We check to see if the translates are also irreducible as S modules. If
## so we return true. If not we  replace S by S' < S^G, until either one of the
## translates is no longer a module for S  (in which case we return false)
## or all are irreducible  (in which case we return  a list of bases for the 
## submodules, together with a list of the submodules.).  
## S may change in this function, whatever it returns, and so may Wmodule.
## 
TranslatesIrreducible := function (matrices, S, Wmodule, translates, F)
   local summands, info, k, W, dimW, Wk, Wkmodule, g, s, s1, s2, i, U, 
         subSmodules, subS;

   summands := translates[1];
   info := translates[2];
   W := summands[1];
   dimW := Length (W);
   subSmodules := [Wmodule];
   k := 2;
   while k <= Length (summands) do
      Wk := summands[k];
      Wkmodule := SubGModule(GModule (S, F), Wk);
      while IsIrreducible (Wkmodule)=false do
         g := matrices[info[k][2]];

         # g translated a previous irreducible Wi to give Wk. 
         # Wk has a submodule which is a translate Ug of a subspace U of Wi.
         # Some s in S doesn't fix U, so g^-1 * s * g doesn't fix Ug. 
         # Add this conjugate to S.
         # A basis for Ug in terms of the basis for Wk is stored as a flag 
         # for Wkmodule. The matrix SubbasisFlag (Wgmodule) * Wk expresses a 
         # basis for Ug in terms of the standard basis for the full space.

         U := SubbasisFlag (Wkmodule) * Wk * g^-1;
         i := 1;
         while i <= Length (S) do
            s := S[i];
            if Length (Base (RowSpace (Concatenation (U, U * s), F))) > Length (U) then
               # s doesn't fix U
               s1 := g^-1 * s * g;
               i := Length (S) + 1;
            else 
               i := i + 1;
            fi;
         od;
         Add (S, s1);

         # Now if this s1 preserves each of the translates we can 
         # carry on, otherwise we have to restart the decomposition.
         if Length (Base (RowSpace (Concatenation (W, W * s1), F))) > dimW then
            return false;
         fi;
         s2 := TranslatesSGModules (matrices, [s1], translates, F); 
         if s2 <> true then
            Add (S, s2);
            return false;
         else
            subS := GeneratorsFlag (Wkmodule);
            s := SubGModuleAction (Wk, s1);
            Add (subS, s);
            Wkmodule := GModule (subS, F);
            # We remake Wkmodule completely because its reducibility may 
            # have changed.  
            # For the remaining modules we just add in the action of the 
            # additional matrix - irreducibility has already been proved.
            for i in [1..k - 1] do
               s := SubGModuleAction (summands[i], s1);
               EnlargeIrreducibleGModule (subSmodules[i], s);
            od;
            # I think that this updates Wmodule too, because Wmodule is simply
            # subSmodules[1]. Put a check in.
            if Length (GeneratorsFlag (subSmodules[1])) <> 
               Length (GeneratorsFlag (Wmodule)) then
               Error ("ERROR: in TranslatesIrreducible.\n");
            fi;
         fi;
      od;
      Add (subSmodules, Wkmodule);
      k := k + 1;
   od;

   return [summands, subSmodules];
end;

#############################################################################
