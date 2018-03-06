#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @(#)$Id: isprimit.g,v 1.1 1997/03/10 13:52:32 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: isprimit.g,v $
#H  Revision 1.1  1997/03/10 13:52:32  gap
#H  VERSION 1.0
#H
#H  Revision 1.3  1997/01/10 11:48:45  fceller
#H  put Eamonn's (final) version into the repository,
#H  made a small fix to 'PermGroupRepresentation',
#H  made a small fix to alternating recognition
#H
#H  Revision 1.2  1997/01/05 10:49:28  fceller
#H  added Eamonn's new version to the repository
#H
#H  Revision 1.1  1996/12/25 09:03:51  fceller
#H  changed long filenames to MS-DOS conform filenames,
#H  the init files are *NOT* yet updated
#H
#H  Revision 1.1  1996/12/22 07:45:09  fceller
#H  added new files to the repository
#H
#H  Revision 1.2  1996/12/12 10:51:27  fceller
#H  new version by Eamonn
#H
#H  Revision 1.1  1996/11/28 13:14:53  fceller
#H  added "smash" and "reducible" to the repository
#H
##
#############################################################################
##
#F  InfoPrim (...)  . . . . . . . . . . . information function for package
#F  InfoPrim2 (...)  . . . . . . . . . . . additional information function 
##
##
if not IsBound(InfoPrim)  then InfoPrim := Ignore;  fi;
if not IsBound(InfoPrim2)  then InfoPrim2 := Ignore;  fi;

#############################################################################
##
#F  BlockSizes (M) . . . . . . . . . return block sizes field of M
##
BlockSizes := function (M)

   if IsBound (M.blockSizes) = false then return "unknown"; fi;
   return M.blockSizes;

end; #BlockSizes

#############################################################################
##
#F  SetBlockSizes (M) . set block sizes of M to SizesOfBlocks
##
SetBlockSizes := function (M, SizesOfBlocks)

   M.blockSizes := SizesOfBlocks;

end; #SetBlockSizes

#############################################################################
##
#F  BlockNumbers (M)  . .  . . . . return block numbers of M
##
BlockNumbers := function (M)

   if IsBound (M.blockNumbers) = false then return "unknown"; fi;
   return M.blockNumbers;

end; #BlockNumbers

#############################################################################
##
#F  SetBlockNumbers (M) . set block numbers of M to NmrOfBlocks
##
SetBlockNumbers := function (M, NmrOfBlocks)

   M.blockNumbers := NmrOfBlocks;

end; #SetBlockNumbers

############################################################################
##
#F  IsBlockSystem (BlockSystem) . . . . . . . . is object a block system? 
##
IsBlockSystem := function (BlockSystem)

   return IsRec (BlockSystem) and IsBound (BlockSystem.isBlockSystem)
          and BlockSystem.isBlockSystem; 

end; #IsBlockSystem

#############################################################################
##
#F  DeletePrimitivityComponents (M)  . . . . . . . . delete components from M 
##
DeletePrimitivityComponents := function (M)

   Unbind (M.blockNumbers);
   Unbind (M.blockSizes);

end; #DeletePrimitivityComponents 

#############################################################################
##
#F  ConstructBlock (M, TensorFactor) . . . M has a tensor decomposition 
##  with first component TensorFactor which is imprimitive; 
##  use the stored block for its block system to write down 
##  a block under the action of M 
##
ConstructBlock := function (M, TensorFactor)

   local r, s, e, f, B, i, j, k, SmallBlockSystem, SmallBlock, LargeBlock, el;

   #dimension of second factor
   e := DimensionFlag (TensorFactorsFlag (M)[2]);

   #dimension of first factor
   f := DimensionFlag (TensorFactor);

   #natural basis for the decomposition
   B := TensorBasisFlag (M);

   #extract block system for action of M on second tensor factor 
   SmallBlockSystem := BlockSystemFlag (TensorFactor);

   #the stored representative block for the small system
   SmallBlock := BlockFlag (SmallBlockSystem);

   #number of blocks in small system
   r := NumberBlocksFlag (SmallBlockSystem);

   #size of each block in the small system
   s := f / r;

   #now use this small block to construct block for larger system
   LargeBlock := [];
   for i in [1..s] do
      for k in [1..e] do
         el := SmallBlock[i][1] * B[(1 - 1)* e + k];
         for j in [2..f] do
            el := el + SmallBlock[i][j] * B[(j - 1)* e + k];
         od;
         Add (LargeBlock, el);
      od;
   od;

   return LargeBlock;

end; #ConstructBlock 

#############################################################################
##
#F  ResolveTensor (M)  . . . . . 
##  M has a tensor decomposition found by SmashGModule 
##  check whether first component of tensor product is imprimitive;
##  if so, use the block found for the first component to construct a
##  block for M, hand this block to MinBlocks (this is strictly unnecessary 
##  but the call sets up proper fields in M), and return TRUE; 
##  if primitive return FALSE; if unknown return "unknown"
##
ResolveTensor := function (M)

   local R, Result, TensorFactor, LargeBlock;

   TensorFactor := TensorFactorsFlag (M)[1];
   InfoPrim ("#I ** About to call StartPrimitivityTest again **\n");
   R := StartPrimitivityTest (TensorFactor, false);
   Result := R[1];

   if (Result = false) then

      InfoPrim2 ("#I ResolveTensor found component is imprimitive\n");

      #use block under first factor to construct a block under action of M  
      LargeBlock := ConstructBlock (M, TensorFactor);

      #now call MinBlocks and set appropriate flags for M 
      TriangulizeMat (LargeBlock);
      SetBlockSystemFlag (M, MinBlocks (M, LargeBlock));

      #temporary code -- remove later 
      if (LargeBlock <> BlockFlag (BlockSystemFlag (M))) then 
         Error ("Blocks are not equal in ResolveTensor\n");
      fi;

      SetBlockSizes (M, [NumberBlocksFlag (BlockSystemFlag (M))]);
      SetBlockNumbers (M, InverseSet (DimensionFlag (M), BlockSizes (M)));

      return true;
   elif (Result = true) then 
      return false;
   else 
      return Result;
   fi;

end; #ResolveTensor

#############################################################################
##
#F  InverseSet (d, Set) . . compute the inverse mod d of each element of Set 
## 
InverseSet := function (d, L)

    local InverseSet;
    InverseSet := List (L, x -> d / x);
    Sort (InverseSet);
    return InverseSet;

end; #InverseSet 

#############################################################################
##
#F  GcdSeq (A, B) . . . .          compute gcd of two numbers given by their
##                                 prime factorisations A and B 
## 
GcdSeq := function (A, B)

   local found, prime, C, D, i, j, lenA, lenB;

   D := Copy (A);

   lenA := Int (Length (A) / 2);
   lenB := Int (Length (B) / 2);

   for i in [1..lenA] do

      found := false;
      prime := D[2 * i - 1];
      j := 1;
      while (j <= lenB) and (not found) do
         if B[2 * j - 1] = prime then
            D[2 * i] := Minimum (D[2 * i], B[2 * j]); 
            found := true;
         fi;
         j := j + 1;
      od;
      if not found then 
         D[2 * i] := 0;
      fi;
   od;

   C := [1, 1];
   j := 1;
   for i in [1..lenA] do
      if D[2 * i] <> 0 then
         C[j] := D[2 * i - 1];
         C[j + 1] := D[2 * i];
         j := j + 2;
      fi; 
   od; 

   return C;

end; #GcdSeq

#############################################################################
##
#F  LcmSeq (A, B) . . . .          compute lcm of two numbers given by their
##                                 prime factorisations A and B 
## 
LcmSeq := function (A, B) 

   local found, D, i, j, prime, lenA, lenB;

   lenA := Int (Length (A) / 2);
   lenB := Int (Length (B) / 2);

   D := Copy (A);
   for i in [1..lenA] do
      found := false;
      prime := D[2 * i - 1];
      j := 1;
      while (j <= lenB) and not found do
         if B[2 * j - 1] = prime then
            D[2 * i] := Maximum (D[2 * i], B[2 * j]); 
            found := true;
         fi; 
         j := j + 1;
      od; 
   od;

   #add in any primes which occur in B but not in A

   for i in [1..lenB] do
      found := false;
      prime := B[2 * i - 1];
      j := 1;
      while (j <= lenA) and not found do
         found := (D[2 * j - 1] = prime);
         j := j + 1;
      od; 

      if (not found) then 
         D[2 * lenA + 1] := prime;
         D[2 * lenA + 2] := B[2 * i];
         lenA := lenA + 1;
      fi; 
   od; 

   return D;

end; #LcmSeq

#############################################################################
##
#F  ExponentGL (d, q) . . . . . . . . . . . .   compute exponent of GL (d, q) 
## 
ExponentGL := function (d, q) 

   local x, n, i, p, t, pow, ExpSeq;

   n := PrimePowersInt (q - 1);
   for i in [2..d] do
      n := LcmSeq (n, PrimePowersInt (q^i - 1));
   od;

   p :=  FactorsInt (q)[1];
   t := 1;
   pow := p;
   while (pow < d) do
      pow := pow * p;
      t := t + 1;
   od;

   ExpSeq := Concatenation ([p, t], n);

   return ExpSeq;

end; #ExponentGL

#############################################################################
##
#F  IsValidSymOrder (Order, r) . . . .         is Order a valid order for an
##                                                 element in Symmetric (r)?
## 
##  an element of order p1^n1 * p2^n2 * ... * pk^nk needs at least 
##  p1^n1 + p2^n2 + ... + pk^nk points 
## 
IsValidSymOrder := function (order, r) 

   local S, total, i;

   S := PrimePowersInt (order);

   total := 0;
   for i in [1..Int (Length (S) / 2)] do
      total := total + S[2 * i - 1]^S[2 * i];
   od;

   return total <= r;

end; #IsValidSymOrder

#############################################################################
##
#F  GcdOrderGL (d, q, Order) . . .      compute the gcd of the order of the 
##                                      element and the exponent of GL(d, q)
##
GcdOrderGL := function (d, q, Order)

   local i, Char, Factors, p, Result, power;

   Factors := Set (FactorsInt (Order));
   Char := FactorsInt (q)[1];

   Result := 1;
   for p in Factors do 
      if p <> Char then 
         i := 1;
         while (OrderMod (q, p^i) <= d) and (Order mod p^i = 0) do
            i := i + 1;
         od; 
         Result := Result * p^(i - 1);
      else
         power := 1;
         while (power < p * d) and (Order mod power = 0) do
            power := power * p;
         od;
         Result := Result * power / p;
      fi;
   od;

   return Result;

end; #GcdOrderGL    

#############################################################################
##
#F  IsOrderValid (d, q, r, Order) . . . 
##  given element in GL(d, q) of order 'Order'; 
##  if Order divides exponent of GL(s,q) wr Sym(r), return true, else false 
##
IsOrderValid := function (d, q, r, Order)

   local matord, permord;

   #compute the gcd of the order of the element and the exponent of GL (d/r, q)
   matord := GcdOrderGL (d / r, q, Order);
   InfoPrim2 ("#I Value of matorder is ", matord, "\n");

   permord := Order / matord;
   InfoPrim2 ("#I Value of permorder is ", permord, "\n");

   return IsValidSymOrder (permord, r);

end; #IsOrderValid

#############################################################################
##
#F OrderOfElement (M, Order) . . . . does an element of this Order
##                                eliminate any of the possible block sizes? 
##
OrderOfElement := function (M, Order)

   local d, q, r, j, NmrOfBlocks;

   if Order = 1 then return; fi;

   NmrOfBlocks := BlockNumbers (M);

   d := DimensionFlag (M);
   q := Length (Elements (FieldFlag (M)));

   for j in [1..Length (NmrOfBlocks)] do
      r := NmrOfBlocks[j];
      InfoPrim2 ("#I r is ", r, "\n");
      if IsOrderValid (d, q, r, Order) = false then 
         InfoPrim ("#I r = ", r, " is invalid\n");
         NmrOfBlocks[j] := 0;
      fi;
   od; 

   NmrOfBlocks := Difference (NmrOfBlocks, [0]);

   SetBlockNumbers (M, NmrOfBlocks);

end; #OrderOfElement

############################################################################
##
#F PolynomialQuotient (F, f, g) . . .
## return quotient of polynomials f by g as polynomial over F 
## if f is not divisible by g then return false
## 
PolynomialQuotient := function (F, f, g)

   local P;

   P := PolynomialRing (F);
   f := EmbeddedPolynomial (P, f);
   g := EmbeddedPolynomial (P, g);

   return Quotient (f, g);

end; #PolynomialQuotient


#############################################################################
##
#F PolynomialCoefficients (f) . . . . . . return coefficients of polynomial f 
## where coeff[i] is the power of x^i in f
##
## leading zeroes are suppressed in f.coefficients;
## f.valuation is the largest power of x which divides f;
## hence we prepend f.valuation zeros to the list of coefficients 
## 
PolynomialCoefficients := function (f)

   local coeffs, i, zero, NmrZeros;

   coeffs := [];
   zero := BaseRingZero (f);
   NmrZeros := Valuation (f);
   coeffs := List ([1..NmrZeros], i -> 0);

   Append (coeffs, PolCoefficients (f));

   return coeffs;

end; #PolynomialCoefficients

#############################################################################
##
#F PolynomialRemainder (F, f, g) . . . return remainder of polynomials f by g 
##                                     as a polynomial over F 
## 
PolynomialRemainder := function (F, f, g)

   local r;

   r := RemainderCoeffs (PolynomialCoefficients (f), 
                PolynomialCoefficients (g));
   return Polynomial (F, r);

end; #PolynomialRemainder

#############################################################################
##
#F  FindLargestPower (TestElement, f, p, n, a, F) . . .     
## 
## if p <> char of F 
## 
## find largest power of x^(p^n) - a in f 
## 
## else 
## 
## find Rank ((x - a)^(p^n - 1)) evaluated for TestElement
## 
FindLargestPower := function (TestElement, p, n, a, f, F)

   local ZeroPol, i, factor, t, quotient, remainder;

   ZeroPol := Zero (PolynomialRing (F));

   if p <> Characteristic (F) then

      #set up x^(p^n) - a over F
      factor := SetupFactor (F, p^n, a);

      #find the largest power, t, of x^(p^n) - a which divides f 
      t := -1;
      repeat 
         quotient := PolynomialQuotient (F, f, factor);
         remainder := PolynomialRemainder (F, f, factor);
         if remainder = ZeroPol then 
            f := quotient;
         fi;
         t := t + 1;
      until remainder <> ZeroPol;

      return [t, f];

   else
      #set up f = (x - a)^(p^n - 1) 
      f := SetupFactor (F, 1, 1)^(p^n - 1); 

      #find d - dimension of NS (f(TestElement)) = rank of f(TestElement) 
      t := RankMat (Value (f, TestElement));

      return [t];

   fi; #p = Char

end; #FindLargestPower

#############################################################################
##
#F  FreeRank (TestElement, p, f, F) . . .     
##  find largest t such that V contains a free F[C_p] module of rank t
## 
FreeRank := function (TestElement, p, f, F)

   local Result, Residue, rank, i, one, factor, t, h, Excess;

   one := One (F);

   if p <> Characteristic (F) then

      #find the largest power of x^p - 1 which divides f 
      #where f = (x^p - 1)^t * Residue

      Result := FindLargestPower (TestElement, p, 1, 1, f, F);
      t := Result[1];
      Residue := Result[2];
      InfoPrim2 ("#I Residue = ", Residue, "\n");

      #now check whether Residue is simply (x - 1)^Degree(Residue)
      factor := Polynomial (F, [-one, one])^Degree (Residue);
      InfoPrim2 ("#I Factor = ", factor, "\n");
      Excess := (factor <> Residue);

      return [t, Excess];

   else

      #set up x - 1
      factor := Polynomial (F, [-one, one]);

      Result := FindLargestPower (TestElement, p, 1, 1, f, F);
      t := Result[1];

      #now find height, h, of residue
      #do this by finding smallest h such that Rank ((x - 1)^h) = t * (p - h)

      f := factor;
      h := 0;
      repeat 
         rank := RankMat (Value (f, TestElement));
         f := f * factor;
         h := h + 1;
      until rank = t * (p - h);

      return [t, h];

   fi; #p = Char

end; #FreeRank

#############################################################################
##
#F  ExamineSmashResult (M, report) . . . . . .  examine results of Smash 
## 
ExamineSmashResult := function (M, report)

   #the following are possible outcomes from the call to SmashGModule: 
   #
   #(a) all commutators of generators are scalar (so the group 
   #    must be central by abelian) and SmashGModule exits with an error 
   #(b) group is imprimitive
   #(c) group is a tensor product
   #(d) group is semilinear 
   #(e) the commutator subgroup acts absolutely irreducibly 
   #
   #(a) may result in an Error; 
   #if (b) we return true; 
   #if (c) we return false; 
   #if (d) we return "unknown"
   #(e) has no impact and we return false;

   if ImprimitiveFlag (M) = true then
      if report then 
         Print ("#I  Matrix group is imprimitive\n");
      fi;
      return true;
   elif TensorProductFlag (M) = true then 
      if report then 
         Print ("#I  GModule is a tensor product\n");
      fi;
      return false;
   elif SemiLinearFlag (M) = true then 
      if report then 
         Print ("#I  GModule is semilinear\n");
      fi;
      return "unknown";
   fi;

   return false;

end; #ExamineSmashResult

#############################################################################
##
#F  AddSmashQueue (TestElement, t, Queue) . . add TestElement to smash queue 
## 
## keep track of element and the rank of the free module
## 
AddSmashQueue := function (TestElement, t, Queue)

   if IsScalar (TestElement) = false then 
      Add (Queue, [TestElement, t]);
   fi;

end; #AddSmashQueue

#############################################################################
##
#F  IndexMinimumRankElement (Queue) . . . 
##  index in Queue of element of smallest rank 
## 
IndexMinimumRankElement := function (Queue)

   local t, x;

   if Length (Queue) = 0 then return 0; fi;

   t := List (Queue, x -> x[2]);

   return Position (t, Minimum (t));

end; #IndexMinimumRankElement

#############################################################################
##
#F  SmashElement (M, TestElement, TensorTest) . . call Smash with TestElement 
## 
##  if TensorTest is true, we are carrying out test for tensor product;
##
##  if SmashGModule finds system of imprimitivity, return true;
##  if SmashGModule finds a tensor product, then it may discover that the 
##  first component is semilinear -- in this case return "unknown";
##  else return false
## 
SmashElement := function (M, TestElement, TensorTest) 

   local SizesOfBlocks, Result, pos;

   SmashGModule (M, [TestElement], "PartialSmash");

   #SmashGModule may find system of imprimitivity or a tensor product 
   #or it may find that the module is semilinear and hence
   #we set the imprimitive flag to "unknown"

   if ImprimitiveFlag (M) = true then 
      SetBlockSizes (M, [NumberBlocksFlag (BlockSystemFlag (M))]);
      return true;
   elif TensorProductFlag (M) = true then 
      if TensorTest = true then return; fi;

      Result := ResolveTensor (M);

      #if our recursive call to PrimitiveTest found that the second
      #component of the tensor decomposition is imprimitive, then 
      #our group is imprimitive; alternately, we may have obtained 
      #semilinear as the result; in either case, we set flag;

      if Result <> false then
         SetImprimitiveFlag (M, Result); 
      fi;

      return Result;
   fi;

   return false;

end; #SmashElement

#############################################################################
##
#F  CallSmashGModule (M, TestElement, t) . call SmashGModule with TestElement 
##
##  if SmashGModule finds system of imprimitivity, return true;
##  if it finds component is semilinear, return "unknown";
##  else return false
## 
CallSmashGModule := function (M, TestElement, t) 

   local SizesOfBlocks, Result, pos;

   Result := SmashElement (M, TestElement, false);

   #SmashGModule may find system of imprimitivity or 
   #may not be able to settle case 

   if (Result <> false) then 
      return Result;
   fi;

   #if SmashGModule did not find system, we can rule out all block sizes > t 
   SizesOfBlocks := BlockSizes (M);
   pos := PositionProperty (SizesOfBlocks, n -> n > t); 
   if pos <> false then 
      SizesOfBlocks := Sublist (SizesOfBlocks, [1..pos - 1]);
      SetBlockSizes (M, SizesOfBlocks);
   fi;

   return false;

end; #CallSmashGModule

#############################################################################
##
#F  CharPolPrimeOrder (M, TestElement, p, Queue) . . .
## 
##  test characteristic polynomial structure of elements of prime order
## 
CharPolPrimeOrder := function (M, TestElement, p, Queue)

   local Remove, first, f, Result, q, t, h, i, s, OrdMod, F, 
         SizesOfBlocks, Excess, Char, x;

   F := FieldFlag (M);
   f := CharacteristicPolynomial (TestElement);

   SizesOfBlocks := BlockSizes (M);

   #   x := Runtime ();
   Result := FreeRank (TestElement, p, f, F);
   #   Print ("#I Time to run FreeRank is ", Runtime () - x, "\n");

   t := Result[1];

   Char := Characteristic (F);

   if p <> Char then

      Excess := Result[2];
      q := Length (Elements (F));
      #smallest integer m such that q^m - 1 is divisible by p
      OrdMod := OrderMod (q, p);

      InfoPrim ("#I In CharPolPrimeOrder, p = ", p, " Excess? ", Excess, 
              " OrdMod = ", OrdMod, " t = ", t, "\n");

      h := 0;
   else 
      h := Result[2];
      Excess := false;
      OrdMod := 0;
      InfoPrim ("#I In CharPolPrimeOrder, p = ", p, " t = ", t, " h = ", h, "\n");
   fi;

   #run over and seek to eliminate block sizes; 
   #store eliminated block sizes in Remove  

   first := true; Remove := [];

   for s in SizesOfBlocks do 
      if (p <> Char and s < OrdMod and Excess) or (p = Char and s < h) then 
         Add (Remove, s);

      elif p = Char and s < p and Mod (t, s) <> 0 then  	 
         Add (Remove, s);

      elif s > t and first then 
         AddSmashQueue (TestElement, t, Queue);
         first := false;
      fi;
   od;      

   SizesOfBlocks := Difference (SizesOfBlocks, Remove);

   SetBlockSizes (M, SizesOfBlocks);

end; #CharPolPrimeOrder

#############################################################################
##
#F  CharPolPrimePowerOrder (M, TestElement, p, n, Queue) . . .
## 
##  order of TestElement is p^n 
##  test characteristic polynomial structure of elements of prime power order
## 
##  compute projective order, p^m, of TestElement
##  TestElement^(p^m) is scalar in a, say
##  now find largest power, t, of x^(p^m) - a 
##  which occurs in char poly of TestElement 
##
CharPolPrimePowerOrder := function (M, TestElement, p, n, Queue)

   local F, f, Result, t, x, o, m, scalar, a, SizesOfBlocks;

   F := FieldFlag (M);
   SizesOfBlocks := BlockSizes (M);

   f := CharacteristicPolynomial (TestElement);

   #note both projective order of TestElement and the scalar a
   Result := ProjectiveOrderMat (TestElement);
   o := PrimePowersInt (Result[1]);
   m := o[2];
   a := Result[2];

   #   x := Runtime ();
   Result := FindLargestPower (TestElement, p, m, a, f, F);
   #   Print ("#I Time to find t is ", Runtime () - x, "\n");

   t := Result[1];

   InfoPrim ("#I In CharPolPrimePowerOrder, p^n = ", p^n, 
             " projective order = ", p^m, " scalar = ", a, " t = ", t, "\n");

   #can we eliminate block sizes? 

   if t < Maximum (SizesOfBlocks) then 
      #keep track of element of prime order for possible later processing
      AddSmashQueue (TestElement^(p^(m - 1)), t, Queue);
   fi;

end; #CharPolPrimePowerOrder

#############################################################################
##
#F  CharPolStructure (M, g, Order, Queue) . . .     
## 
##  examine characteristic polynomial of elements of prime-power order which 
##  can be obtained as powers of g, an element of order Order
## 
CharPolStructure := function (M, g, Order, Queue)

   local F, factors, p, i, TestElement, powers, p, n, SizesOfBlocks;

   F := FieldFlag (M);

   SizesOfBlocks := BlockSizes (M);

   #consider elements of prime order 

   factors := Set (FactorsInt (Order));
   for p in factors do

      InfoPrim2 ("#I Call CharPolPrimeOrder with element of order ", p, "\n");

      TestElement := g^(Order / p);

      if IsScalar (TestElement) = false then 

         InfoPrim2 ("#I Before call, Blocksizes are ", SizesOfBlocks, "\n");
         CharPolPrimeOrder (M, TestElement, p, Queue);
         SizesOfBlocks := BlockSizes (M);

         InfoPrim ("#I Blocksizes after CharPolPrimeOrder: ", SizesOfBlocks, "\n");
         if SettleComputation (M, Queue) then return; fi;

      else 
         InfoPrim ("#I Element of order ", p, " is scalar\n");
      fi;

   od; 

   #now consider elements of prime-power order, p^n, where n > 1

   powers := PrimePowersInt (Order);
   for i in [1..Length (powers) / 2] do

      p := powers[2 * i - 1];
      n := powers[2 * i];

      if n > 1 then 

         InfoPrim2 ("#I Call CharPolPrimePowerOrder with element of order ", p^n, "\n");

         TestElement := g^(Order / p^n);

         if IsScalar (TestElement) = false then 
            CharPolPrimePowerOrder (M, TestElement, p, n, Queue);
            if SettleComputation (M, Queue) then return; fi;

         else 
            InfoPrim ("#I Element of order ", p^n, " is scalar\n");
         fi;

      fi; #if n > 1

   od; 

end; #CharPolStructure

#############################################################################
##
#F  CompositeOrders (M, Elts, Orders) . . .     
## 
## given list of elements and their projective orders
## select an element x of composite order o
## if IsValidSymOrder (o, r) = false then
## the element x cannot act faithfully on all of the blocks
## let p run over the primes dividing o
## then one of x^(o / p) must fix all blocks
## hand each of these elements to SmashGModule 
## if we do not find a decomposition we have ruled out 
## all number of blocks <= r 
## 
## hence, we find the element which fails IsValidSymOrder
## for the largest possible value of r and hand all of its
## prime constituents to SmashGModule 
##
CompositeOrders := function (M, Elts, Orders)

   local Result, o, rem, l, i, j, best, pos, g, Order, p, x, r, NmrOfBlocks;

   #first find set of composite orders 
   o := Set (Filtered (Orders, x -> Length (Set (FactorsInt (x))) > 1));
   InfoPrim2 ("#I Projective orders for CompositeOrders  is ", Orders, "\n");
   InfoPrim ("#I Input orders for CompositeOrders test is ", o, "\n");

   NmrOfBlocks := BlockNumbers (M);

   #now find which values of r fail IsValidSymOrder 
   rem := [];
   for i in [1..Length (o)] do
      l := [];
      for j in [1..Length (NmrOfBlocks)] do
         r := NmrOfBlocks[j];
         if IsValidSymOrder (o[i], r) = false then
            Add (l, r);
         fi;
      od;
      InfoPrim ("#I o is ", o[i], " invalid r is ", l, "\n");
      if l <> [] then 
         rem[i] := Maximum (l);
      fi;
   od;

   InfoPrim2 ("#I Rem is ", rem, "\n");

   if Length (rem) = 0 then return false; fi;

   #maximum r encountered which fails required test
   best := Maximum (rem);
   pos := Position (rem, best);
   pos := Position (Orders, o[pos]);

   #note the element g which fails for maximum r and its order 
   g := Elts[pos];
   Order := Orders[pos];

   InfoPrim ("#I Use element of order ", Order, " to rule out r <= ", best, "\n");

   #apply SmashGModule to all elements of prime order obtained from g
   for p in Set (FactorsInt (Order)) do
      InfoPrim ("#I Now calling SmashGModule with element of order ", p, "\n");
      x := g^(Order / p);
      Result := SmashElement (M, x, false);
      if Result <> false then
         return Result;
      fi;
   od;

   #we can now rule out all number of blocks <= best 
   SetBlockNumbers (M, Filtered (NmrOfBlocks, x -> x > best));

   return false;

end; #CompositeOrders

#############################################################################
##
#F  SemiLinearTest (M)  . . . . . . . . .  test if the module is semi-linear 
## 
SemiLinearTest := function (M)

   local gens, S, x, i; 

   gens := GeneratorsFlag (M);

   S := Commutators (gens);

   #does S consist entirely of scalars?
   #if so, add a non-scalar generator to S; if the group 
   #does not have a non-scalar generator, it is reducible 
   #and should be eliminated earlier in BasicReductionTests 

   if ForAll (S, x -> IsScalar (x)) then 
      i := PositionProperty (gens, x -> IsScalar (x) = false);
      if i <> false then 
         Add (S, gens[i]);
      fi;
   fi;

   SmashGModule (M, S, "PartialSmash");

end; #SemiLinearTest

#############################################################################
##
#F  BasicReductionTests (M) . . .     
##  carry out tests for basic reductions of the module M 
##
BasicReductionTests := function (M)

   local Result;

   #is the module irreducible?
   if IsIrreducible (M) = false then 
      Print ("#I  GModule is not irreducible\n");
      return true;
   fi;

   #is the module absolutely irreducible?
   if IsAbsolutelyIrreducible (M) = false then 
      Print ("#I  GModule is not absolutely irreducible\n");
      return true;
   fi;

   #test for semilinearity
   SemiLinearTest (M);

   Result := ExamineSmashResult (M, true);

   if TensorProductFlag (M) = true then 
      Result := ResolveTensor (M);

      #if our recursive call to PrimitiveTest found that a component of 
      #the tensor decomposition is imprimitive, we can deduce that our 
      #group is imprimitive; alternately, it may have found that the 
      #associated module is semilinear

      if Result <> false then
         SetImprimitiveFlag (M, Result); 
      fi;

   fi;

   return Result;

end; #BasicReductionTests

#############################################################################
##
#F  FinishComputation (M, Queue, ProcessLeast) . . . . . . . 
##
##  if there is an element of rank < the smallest remaining block size, 
##  then finish the computation by calling SmashGModule;
##  even if this is not so, if ProcessLeast is true, then call SmashGModule 
##  with minimum rank element from SmashGModule Queue;
##  set appropriate flags; return true if block system found or
##  primitivity proved, else false
## 
FinishComputation := function (M, Queue, ProcessLeast)

   local Index, Result, TestElement, t, SizesOfBlocks;

   SizesOfBlocks := BlockSizes (M);

   #do we already know the answer?
   Result := ImprimitiveFlag (M);
   if Result = true or Result = false then
      return true;
   fi;

   if Length (SizesOfBlocks) <> 0 then 

      #note index position of element of least rank in SmashGModule queue
      Index := IndexMinimumRankElement (Queue);

      #if this is smaller than remaining valid block sizes or
      #ProcessLeast is true, call SmashGModule

      if Index <> 0 then 
         t := Queue[Index][2];
         if ProcessLeast or (Minimum (SizesOfBlocks) > t) then 

            InfoPrim ("#I Call SmashGModule with element of rank t = ", t, "\n");

            TestElement := Queue[Index][1];
            Result := CallSmashGModule (M, TestElement, t); 
            if (Result <> true) and (Result <> false) then
               return true;
            fi;

            SetBlockNumbers (M, InverseSet (DimensionFlag (M), BlockSizes (M)));
            if Result = true then 
               return true;
            fi;
         fi;
      fi;

   fi;

   if Length (BlockSizes (M)) = 0 then 
      SetImprimitiveFlag (M, false);
      return true; 
   fi;

   return false;

end; #FinishComputation 

#############################################################################
##
#F  SettleComputation (M, Queue) . . . . . . . 
##
##  have we eliminated all valid block sizes?
##  can we settle computation via a single call to SmashGMod?
##  if either is the case, return true, else return false
##
SettleComputation := function (M, Queue)

   local Index, SizesOfBlocks;

   SizesOfBlocks := BlockSizes (M);

   if Length (SizesOfBlocks) = 0 then return true; fi;

   #note index position of element of least rank in SmashGModule queue
   Index := IndexMinimumRankElement (Queue);

   #is this smaller than existing valid block sizes?
   if Index <> 0 then 
      return (Minimum (SizesOfBlocks) > Queue[Index][2]);
   fi;

   return false;

end; #SettleComputation 

#############################################################################
##
#F  PrimitiveTest (M, NmrElements) . . .  test module M for primitivity 
## 
##  NmrElements is the number of random elements to choose
##  function sets ImprimitiveFlag (M) to be true or false
##  if true, BlockSystemFlag (M) contains a block system
##
PrimitiveTest := function (M, NmrElements)

   local Orders, o, ProjectiveOrders, R, r, t, MaxDegree, E, 
         i, d, g, Result, Elts, NmrElts, Queue;

   d := DimensionFlag (M);

   Elts := []; Orders := []; ProjectiveOrders := [];
   Queue := []; NmrElts := 0;

   repeat 
      g := PseudoRandom (M);
      NmrElts := NmrElts + 1;

      o := OrderMat (g);
      
      if not (o in Orders) and (o <> 1) then

         InfoPrim ("#I Order of element is ", o, "\n");

         #keep this element for later tests 
         Add (Elts, g);
         Add (Orders, o);
         Add (ProjectiveOrders, ProjectiveOrderMat (g)[1]);

         #seek to eliminate possible numbers of blocks 
         OrderOfElement (M, o);
         SetBlockSizes (M, InverseSet (d, BlockNumbers (M)));

         InfoPrim ("#I After order test, block sizes are ", BlockSizes (M), "\n");

         #are we finished?
         if FinishComputation (M, Queue, false) then return; fi;

         #seek to eliminate possible sizes of blocks 
         CharPolStructure (M, g, o, Queue);
         SetBlockNumbers (M, InverseSet (d, BlockSizes (M)));

         #are we finished?
         if FinishComputation (M, Queue, false) then return; fi;
      fi;

   until NmrElts = NmrElements;

   #eliminate all block sizes > minimum value of rank 
   if FinishComputation (M, Queue, true) then return; fi;

   #can we use composite elements to rule out certain block numbers?
   Result := CompositeOrders (M, Elts, ProjectiveOrders);

   #a call to SmashGModule may occur during the call to CompositeOrders  
   #resulting in a "unknown"
   if Result <> true and Result <> false then
      return Result;
   fi;

   SetBlockSizes (M, InverseSet (d, BlockNumbers (M)));

   #are we finished?
   if FinishComputation (M, Queue, false) then return; fi;

   InfoPrim ("#I Possible block numbers after Order, CharPolStructure, CompositeOrder tests are ", 
           BlockNumbers (M), "\n");

   #for each remaining number r of blocks 
   #apply BlockStabiliserTest to either find blocks or 
   #show that there are no block systems of r blocks 

   for r in BlockNumbers (M) do

      InfoPrim ("#I Calling BlockStabiliserTest with r = ", r, "\n");

      Result := BlockStabiliserTest (M, Elts, ProjectiveOrders, r);

      if Result = true then 
         InfoPrim ("#I Successfully eliminated r = ", r, "\n");
         SetBlockNumbers (M, Difference (BlockNumbers (M), [r]));
      elif Result = false then  
         r := NumberBlocksFlag (BlockSystemFlag (M));
         SetBlockNumbers (M, [r]);
         SetBlockSizes (M, [d / r]);
         SetImprimitiveFlag (M, true);
         return;
      else 
         #SmashGModule call found module is semilinear 
         return;
      fi;
   od;

   SetBlockSizes (M, InverseSet (d, BlockNumbers (M)));

   FinishComputation (M, Queue, false);

end; #PrimitiveTest

#############################################################################
##
#F  ReportResult (M)  . . . . . . have we discovered that M is (im)primitive? 
##
ReportResult := function (M, PrintFlag)

   local Result;

   #is the group imprimitive?
   Result := ImprimitiveFlag (M);

   #do we know the answer?
   if Result <> true and Result <> false then return Result; fi;

   #is the group primitive?
   Result := not Result;

   if PrintFlag then 
      if Result = false then 
         Print ("#I  Number of blocks is ", 
                 NumberBlocksFlag (BlockSystemFlag (M)),"\n");
      fi;
      Print ("#I  Group primitive? ", Result, "\n");
   fi;

   #delete those components now unnecessary 
   DeletePrimitivityComponents (M);

   return Result;

end; #ReportResult

#############################################################################
##
#F  StartPrimitivityTest (M, PrintFlag)  . . . . . . begin test on module M 
##    
StartPrimitivityTest := function (M, PrintFlag)

   local Result, NmrElements, d;

   #can we reduce the problem?
   if BasicReductionTests (M) <> false then 
      Result := ReportResult (M, PrintFlag);
      return [Result, M]; 
   fi;

   d := DimensionFlag (M);

   #initialise these components 
   if BlockNumbers (M) = "unknown" then 
      SetBlockNumbers (M, Difference (DivisorsInt (d), [1]));
      SetBlockSizes (M, Difference (DivisorsInt (d), [d]));
   fi;

   #now test for (im)primitivity
   NmrElements := 20;
   PrimitiveTest (M, NmrElements);

   Result := ReportResult (M, PrintFlag);

   return [Result, M]; 

end; #StartPrimitivityTest 

#############################################################################
##
#F  IsPrimitive (G [, <L>])   . .  is matrix group or G-module primitive? 
##    
##  L is a list of integer factorisations of dimension of G 
##  the function returns list containing a boolean and a module M 
##  if boolean is false, then BlockSystemFlag (M) returns the block system
##
Smash.IsPrimitive := function (arg)

   local G, M, L, Result, d;

   if Number (arg) = 1 then
      G := arg[1];
   elif Number (arg) >= 2  then
      G := arg[1];
      L := arg[2];
   fi;

   if IsMatGroup (G) = false and IsGModule (G) = false then
      return Error ("First argument must be a matrix group or GModule\n");
   fi;

   d := DimensionFlag (G);
   if IsGModule (G) = false then
      M := GModule (G);
   else
      M := G;
   fi;

   # possible dimensions of tensor factors
   if IsBound (L) = true then
      if IsList (L) = false or ForAll (L, x -> IsList (x) and
                 Length (x) = 2 and x[1] * x[2] = d) = false then return
          Error ("Second argument must be a list of factorisations of ", d, "\n");
      else 
         SetBlockNumbers (M, Difference (List (L, x -> x[1]), [1]));
         SetBlockSizes (M, Difference (List (L, x -> x[2]), [d]));
      fi;
   fi;

   InfoPrim ("#I Input GModule has dimension ", DimensionFlag (M), " over ", 
          FieldFlag (M), "\n");

   Result := StartPrimitivityTest (M, true);

   return Result;

end; #IsPrimitive
