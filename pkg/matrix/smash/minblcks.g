#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @(#)$Id: minblcks.g,v 1.1 1997/03/10 13:52:34 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: minblcks.g,v $
#H  Revision 1.1  1997/03/10 13:52:34  gap
#H  VERSION 1.0
#H
#H  Revision 1.2  1997/01/05 10:49:30  fceller
#H  added Eamonn's new version to the reprository
#H
#H  Revision 1.1  1996/12/25 09:04:21  fceller
#H  changed long filenames to MS-DOS conform filenames,
#H  the init files are *NOT* yet updated
#H
#H  Revision 1.4  1996/12/23 21:50:55  fceller
#H  added Eamonn's new versions
#H
#H  Revision 1.2  1996/12/10 12:05:45  fceller
#H  added new versions to the repository
#H
##

if not IsBound(InfoMinBlocks)  then InfoMinBlocks := Ignore;  fi;

#############################################################################
## compute block system of matrix group which contains supplied vector 
## 
#############################################################################
##
#F  FilterVector (v, W, A, one) . . . . . filter vector v through echelonised 
##  basis W and update coefficient matrix A
## 
##
FilterVector := function (v, W, A, one)

   local depth, dim, a, c;

   depth := DepthVector (v);
   a := 0 * v;
   dim := Length (v);

   while IsBound (W[depth]) do
      c := W[depth][depth]^-1 * v[depth];
      v := v - c * W[depth];
      a := a - c * A[depth];
      depth := DepthVector (v);
   od;

   if depth <= dim then  
      W[depth] := v;
      a[depth] := one;
      A[depth] := a;
      return [true, v];
   else
      return [false, a];
   fi;

end; #FilterVector

#############################################################################
##
#F  Rep (i, BT) . . .                   find representative of vector i in BT
## 
##
Rep := function (i, BT)

   local j;

   if i <> 0 then 
      j := BT[i][Length (BT[i])];
      if j < 0 then return Rep (-j, BT); fi;
   fi;

   return i;

end; #Rep

#############################################################################
##
#F  SetRep (i, j, BT) . . .      set representative of vector i to be j in BT
## 
##
SetRep := function (i, j, BT)

   BT[i][Length (BT[i])] := -j;

end; #SetRep 

#############################################################################
##
#F  SetImage (i, j, k, BT) . . .  set image of vector i under generator j to 
##                                  be in BT[k]
## 
##
SetImage := function (i, j, k, BT)

   BT[i][j] := k;

end; #SetImage 

#############################################################################
##
#F  InverseColumn (i, NmrColumns) . . .  given column for generator i, 
##                                       what is column for generator i^-1? 
## 
##
InverseColumn := function (i, NmrColumns)

   local NmrGens;

   NmrGens := Int (NmrColumns / 2);
   if i <= NmrGens then 
      return i + NmrGens;
   else
      return i - NmrGens;
   fi;

end; #InverseColumn

#############################################################################
##
#F  EquateBlocks (c, d, BT) . . . equate blocks c and d  of block table BT
##                                and clean up all consequences 
## 
##
EquateBlocks := function (c, d, BT)

   local NmrColumns, i, a, b, x;

   if c = d then return; fi;

   SetRep (d, c, BT);

   NmrColumns := Length (BT[1]);

   for i in [1..NmrColumns - 1] do
      a := Rep (BT[d][i], BT);
      if a <> 0 then 
         c := Rep (c, BT);
         x := Rep (BT[a][InverseColumn (i, NmrColumns)], BT); 
         if x <> 0 then 
            EquateBlocks (c, x, BT);
            a := Rep (a, BT);
            c := Rep (c, BT); 
         fi;
         SetImage (a, InverseColumn (i, NmrColumns), c, BT);
         b := Rep (BT[c][i], BT);
         if b = 0 then 
            SetImage (c, i, a, BT);
         elif a <> b then 
            EquateBlocks (a, b, BT);
         fi;
      fi;
   od;

end; #EquateBlocks

#############################################################################
##
#F  CountBlocks (BT) . . .                 count number of blocks in table BT 
## 
##
CountBlocks := function (BT)

   local i, NmrColumns, NmrBlocks;

   NmrColumns := Length (BT [1]);

   NmrBlocks := 0;
   for i in [1..Length (BT)] do
      if (BT[i][NmrColumns] > 0) then
         NmrBlocks := NmrBlocks + 1;
      fi;
   od;

   return NmrBlocks;

end; #CountBlocks

#############################################################################
##
#F  AmalgamateBlocks (S, Bk, BT) . . .    amalgamate those blocks of BT whose 
##                              indices are listed in the set S and update Bk
## 
##
AmalgamateBlocks := function (S, Bk, BT)

   local i;

   for i in [2..Length (S)] do
      EquateBlocks (Rep (S[1], BT), Rep (S[i], BT), BT);
   od;

   for i in [1..Length (Bk)] do
      Bk[i] := Rep (Bk[i], BT);
   od;

end; #AmalgamateBlocks


#############################################################################
##
#F  InitialiseBlock (b, Width, BT) . . .   initialise block b of table BT to
##                                         consist of 0 with last entry b 
## 
##
InitialiseBlock := function (b, Width, BT)

   local i;

   BT[b] := [];
   for i in [1..Width - 1] do
      BT[b][i] := 0;
   od;

   BT[b][Width] := b;

end; #InitialiseBlock


#############################################################################
##
#F  BlockImage (M, i, j, SpaceBasis, W, A, Bk, BT, Map, NextCoset) . . .     
##  find image of block i under generator j
## 
##
BlockImage := function (M, i, j, SpaceBasis, W, A, Bk, BT, Map, NextCoset)

   local matrices, b, c, v, w, x, S, d, Len, F, zero, one,
         flag, Result, NmrColumns;

   F := FieldFlag (M);
   NmrColumns := Length (BT[1]);

   #Bk[i] is the block containing the ith vector in SpaceBasis
   c := Bk[i];

   #let b the image of block c under generator j 
   b := BT[c][j];

   v := SpaceBasis[i];
   matrices := GeneratorsFlag (M);
   w := v * matrices[j];

   one := One (F);
   zero := Zero (F);

   Result := FilterVector (w, W, A, one);

   flag := Result[1];
   x := Result[2];

   if flag <> true then 
      #w lies in the space spanned by SpaceBasis 

      #set S to be the mimimum set of blocks whose direct sum contains w
      S := [];
      for d in [1..Length (x)] do
         if x[d] <> zero then 
            AddSet (S, Bk[Map[d]]);
         fi;
      od;

      #if w lies in just one block and the image of c under j is undefined
      # then define the image of c under j to be this block
      if (Length (S) = 1 and b = 0) then
         SetImage (c, j, S[1], BT);
         SetImage (S[1], InverseColumn (j, NmrColumns), c, BT);
      else 
         if b <> 0 then
            AddSet (S, b);
         fi;
         #now amalgamate all of the blocks in S into one
         if Length (S) > 1 then 
            AmalgamateBlocks (S, Bk, BT);
         fi;
      fi;
   else
      #w does not lie in the space spanned by SpaceBasis; 
      #add in new basis element 
      #the image of c is undefined
      if b = 0 then 
         b := NextCoset;
         NextCoset := NextCoset + 1;
         InitialiseBlock (b, Length (BT[1]), BT);
         SetImage (c, j, b, BT);
         SetImage (b, InverseColumn (j, NmrColumns), c, BT);
      fi;
      #increment length of basis 
      Len := Length (SpaceBasis) + 1;
      Map[DepthVector (x)] := Len;
      SpaceBasis[Len] := w;
      Bk[Len] := b;
   fi;

   return NextCoset;

end; #BlockImage    


#############################################################################
##
#F  StandardVector (Pos, Len, F) . . .    set up standard basis vector e_Pos
##                                        of length Len over F
## 
##
StandardVector := function (Pos, Len, F) 

   local e, i, zero, one;

   zero := Zero (F);
   one := One (F);

   e := [];
   for i in [1..Len] do
      e[i] := zero;
   od;
   e[Pos] := one;

   return e;

end; #StandardVector


#############################################################################
##
#F  SetupBasis (SpaceBasis) . . . .  set up the elements of SpaceBasis in 
##  array W where W[i] is a vector of depth i  -- there may be some gaps in W 
## 
##
SetupBasis := function (SpaceBasis)

   local i, depth, W;

   W := [];
   for i in [1..Length (SpaceBasis)] do
      depth := DepthVector (SpaceBasis[i]);
      if IsBound (W[depth]) then
         return false;
      else
         W[depth] := SpaceBasis[i];
      fi;
   od;

   return W;

end; #SetupBasis


#############################################################################
##
#F  SetupBlockPermGroup (BT, B, Bk) . . set up permutation group representing 
##  action on blocks; also find the block containing the input basis and the 
##  number of blocks 
##
##
SetupBlockPermGroup := function (BT, B, Bk)

   local LastColumn, BlockNumber, i, j, number, renumber, livenumber,
     BlockSystem, NmrBlocks, Block, Perms, PermGroup, row, Maps, nonid;

   LastColumn := Length (BT[1]);
   renumber := []; livenumber := [];

   BlockNumber := 0;
   for i in [1..Length (BT)] do
      number := BT[i][LastColumn];
      if number > 0 then 
         BlockNumber := BlockNumber + 1;
         renumber[i] := BlockNumber;
         livenumber[BlockNumber] := i;
      fi;
   od;

   NmrBlocks := Length (livenumber);

   #block containing the given collection of vectors 
   BlockNumber := Bk[1]; 
   Block := [B[1]];

   for i in [2..Length (B)] do
      if Bk[i] = BlockNumber then 
         Add (Block, B[i]);
      fi;
   od;

   #now write down the permutations 
   Maps := []; nonid := 0;
   Perms := [];
   for i in [1..(LastColumn - 1) / 2] do  
      Perms[i] := [];
      for j in [1..NmrBlocks] do
         row := livenumber[j];
         Perms[i][j] := renumber[BT[row][i]];
      od;
      Perms[i] := PermList (Perms[i]); 
      if Perms[i] = () then 
         Maps[i] := 0;
      else 
         nonid := nonid + 1;
         Maps[i] := nonid;
      fi;
   od;

   #set up the permutation group 
   PermGroup := Group (Perms, Perms[1]^0);

   #now set up block system as record 
   BlockSystem := rec (numberBlocks := NmrBlocks, block := Block,
                       maps := Maps,
                       permGroup := PermGroup, isBlockSystem := true); 

   return BlockSystem;

end; #SetupBlockPermGroup 

#############################################################################
##
#F  MinBlocks (M, B) . . .     find smallest block containing the echelonised 
##                             basis B under the action of module M 
## 
## return number of blocks; the block containing the supplied basis B
## and the permutation group which describes the action on the blocks   
## 
## It is assumed that M is irreducible and that B is a basis 
##
MinBlocks := function (M, B)

   local Map, d, i, j, F, A, BT, Bk, W, NextCoset, SpaceBasis, NmrGens;

   #Bk[i] is the block containing the ith vector in SpaceBasis

   #W is the echelonised basis with vectors arranged according to depth

   #A is coefficient matrix 

   #BT[i][j] is the image of block i under the action of generator j 
   #where the images under inverse of generator i is stored in 
   #column i + NmrGens

   #Map[i] = j says that element j of SpaceBasis is stored as vector i of W 
   #(that is, the depth of j is i)

   SpaceBasis := Copy (B);

   W := SetupBasis (SpaceBasis);
   if (W = false) then
      Print ("#I Supplied collection of vectors is not an echelonised basis\n");
      return false;
   fi;

   Map := [];
   for i in [1..Length (B)] do 
      Map[DepthVector (B[i])] := i; 
   od;

   d := Length (B[1]);
   F := FieldFlag (M);

   A := [];
   Bk := [];

   for i in [1..Length (W)] do
      if IsBound (W[i]) then 
         A[i] := StandardVector (i, d, F);
      fi;
   od;

   for i in [1..Length (SpaceBasis)] do 
      Bk[i] := 1;
   od;

   NmrGens := Length (GeneratorsFlag (M));

   #length of row in block table is the number of generators * 2 + 1 
   BT := [];
   InitialiseBlock (1, 2 * NmrGens + 1, BT);

   NextCoset := 2;
   for i in [1..d] do
      for j in [1..NmrGens] do
         NextCoset := BlockImage (M, i, j, SpaceBasis, W, A, Bk, BT, Map, NextCoset);
      od;
   od;

   InfoMinBlocks ("#I The number of blocks is ",  CountBlocks (BT), "\n");

   return SetupBlockPermGroup (BT, SpaceBasis, Bk);

end; #MinBlocks
