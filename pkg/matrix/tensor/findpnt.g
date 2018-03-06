#############################################################################
##
#A  Matrix package                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
############################################################################
##
#F  InfoTensor1  (...)  . . . . . . . . . . . . . .  for debugging assistance
##
##
if not IsBound (InfoTensor1)  then InfoTensor1 := Ignore;  fi;
#############################################################################
# do the projectivities generate a field? this is not a conclusive test;
# we check that they commute and that we can find a generating element

ProjectivitiesGenerateField := function (S, order)
    local y, grp, g;
    
    for y in S do
        if ForAny (S, x -> x * y <> y * x) then
            InfoTensor1 ("#I Non-commuting\n");
            return [false, false];
        fi;
    od;
    
    grp := Group (S, S[1]^0);
    g := ElementOfOrder (grp, order, 100);
    
    if not IsMat (g) then
        InfoTensor1 ("#I Didn't find element of appropriate order\n");
        return [false, false];
    fi;
    
    return [true, g];
    
end;

# X is a matrix, return the non-zero blocks of size k in X

BlocksOfMatrix := function (X, k)
    local F, d, Nmr, Blocks, i, j, A;
    
    d := Length (X);

    Nmr := QuoInt (d, k);
    Blocks := [];
    
    for i in [1..Nmr] do
        for j in [1..Nmr] do
            A := X{[(i - 1) * k + 1..i * k]}{[(j - 1) * k + 1..j * k]};
            if not IsScalar (A) then 
                if not A in Blocks then
                    Add (Blocks, A); 
                fi;
            fi;
        od;
    od;
    
    return Blocks;
    
end;

# rewrite basis for P wrt N

ConstructNewFlat := function (N, P)
    local Coeffs, B, v, F, x, y;
    
    Coeffs := Base (N);
    B := Base (P);
    v := List ([1..Length (Coeffs)], 
                x -> Sum (List ([1..Length (B)], y -> B[y] * Coeffs[x][y])));
    F := BaseRing (N);        #this is P in Magma code 
    return VectorSpace (v, F);
    
end;

# can we use the potential projectivity C to give us
# a singular element?

FoundSingularElement := function (C, DimU, P, F)
    local R, f, flag, factors, h, N, m;
    
    # compute its characteristic polynomial
    f := CharacteristicPolynomial (C);
    
    # is f the DimU power of some polynomial ?
    R := IsPowerOfPolynomial (f, DimU);
    flag := R[1];
    factors := R[2];
    if not flag then return [true, true]; fi;
    
    # is f a power of an irreducible polynomial ?
    if Length (factors) > 1 then
        InfoTensor1 ("#I CP ... compute eigenspace \n");
        # compute generalised eigenspace 
        h := factors[1][1];
        N := NullSpace (Value (h, C), F);
        return [true, ConstructNewFlat (N, P)];
    else
        m := MinimalPolynomial (C);
        factors := Collected (Factors (m));
        if factors[1][2] > 1 then
            InfoTensor1 ("#I MP...compute eigenspace\n");
            # compute generatlised eigenspace of C
            h := factors[1][1];
            N := NullSpace (Value (h, C), F);
            return [true, ConstructNewFlat (N, P)];
        else
            return [false, false];
        fi;
    fi;
    
end;

# Y collection of matrices written wrt geometric basis;
# P is a potential flat; we are searching for a point
# of dimension DimU

SetupBlocks := function (Y, P, DimU, F)
    local R, zero, Proj, Blocks, DimP, i, B, N, Ainv, C, a, b, j;

    Proj := Set ([]);
    Blocks := [];
    DimP := Dimension (P);
    zero := Zero (F);
    
    # look through DimP x DimP in matrix for each element of Y
    # for one with determinant 0
    for i in [1..Length (Y)] do
        if not IsScalar (Y[i]) then
            B := BlocksOfMatrix (Y[i], DimP);
            for b in B do
                if Determinant (b) = zero then
                    InfoTensor1 ("#I Block A has nullspace\n");
                    N := NullSpace (b, F);
                    return [true, ConstructNewFlat (N, P)];
                fi;
            od;
            Add (Blocks, B);
        fi;
    od;
        
    # construct potential projectivities
    for i in [1..Length (Blocks)] do
        if Length (Blocks[i]) <> 0 then
            Ainv := Blocks[i][1]^-1;
            for j in [2..Length (Blocks[i])] do
                C := Blocks[i][j] * Ainv;
                if not IsScalar (C) then 
                    R := FoundSingularElement (C, DimU, P, F);
                    a := R[1];
                    b := R[2];
                    if a then return [true, b]; fi;
                    if not C in Proj then
                        Add (Proj, C);
                    fi;
                fi;
            od;
        else
            InfoTensor1 ("#I Number of Blocks from A is 0");
        fi;
    od;

    return [false, Proj];
    
end;
    
# search for singular element in algebra generated by collection
# of projectivities CC

SearchForSingularElement := function (C, P, DimU, NmrTries, F)
    local R, C, module, i, x, y, z, a;

    module := GModule (C); #just to allow us to get random elements 
    i := 0;
    repeat
        x := PseudoRandom (module);
        y := PseudoRandom (module);
        z := x + y;
        if not IsScalar (z) then
            R := FoundSingularElement (z, DimU, P, F);
            a := R[1];
            if a then return R; fi;
            if not z in C then Add (C, z); fi;
        fi;
        i := i + 1;
    until i = NmrTries;
    
    return [false, C];
    
end;

# try to find singular element in algebra of projectivities

InvestigateMatrices := function (Y, P, DimU, NmrTries, F)
    local R, flag, C;
    
    R := SetupBlocks (Y, P, DimU, F);
    flag := R[1]; 
    if flag = true then return R; fi;

    C := R[2];
    if Length (C) = 0 then return [false, false]; fi;
    
    R := SearchForSingularElement (C, P, DimU, NmrTries, F);
    return R;
    
end;

# set up matrix whose rows are the vectors of the 
# bases for each subspace in SumSpaces

SetupMatrix := function (SumSpaces)
   local x, A;
   A := [];
   for x in SumSpaces do
      Append (A, Base (x));
   od;
   return A;
end;

# find minimal set of components of SumSpaces whose direct sum contains P;
# U is the inclusion of P into V wrt the given basis for P and
# for the basis of V obtained by concatenating the given
# bases for the direct summands
IdentifySubset := function (A, P, SumSpaces)
    local zero, Ainv, U, DegU, k, Dim, Indices, i, j, ind, v;
    
    Ainv := A^-1;
    zero := 0 * A[1][1];
    
    U := List (Base (P), v -> v * Ainv);
    DegU := Length (U[1]);
    k := Length (SumSpaces);
    
    Dim := Dimension (SumSpaces[1]);
    Indices := Set ([]);
    for i in [1..Length (U)] do
        for j in [1..DegU] do
            if U[i][j] <> zero then
                ind := QuoInt (j - 1, Dim) + 1;
                if not ind in Indices then
                    Add (Indices, ind);
                fi;
                if Length (Indices) = k then return [Indices, U]; fi;
            fi;
        od;
    od;
    
    return [Indices, U];
    
end;

ComponentsOfSum := function (P, SumSpaces, index, M)
    local B, i, j;
    
    B := Base (SumSpaces[index]);
    return List ([1..Length (M[1])], i -> Sum (List ([1..Length (B)], 
                  j -> M[i][j] * B[j])));
    
end;

# extract the (index)th matrix of dimension DimP x DimP from U

BasisMatrix := function (P, U, index)
    local DimP, m, i, j;
    
    DimP := Dimension (P);
    
    m := List ([1..Length (U)], 
        i -> List ([(index - 1) * DimP + 1 .. index * DimP], j -> U[i][j]));
    return m;
    
end;
    
# construct the change of basis matrix and return it
# together with the generators of G wrt the new basis

ConstructMatrices := function (G, P, Equated, A)
    local k, DimP, i, j, y, x, C, z;
    
    k := Length (Equated);
    DimP := Dimension (P);
    
    x := [];
    for i in [1..Length (Equated)] do
        y := [];
        for j in [1..Length (Equated[i])] do
            y[j] := Sum (List ([1..Length (Equated[i][j])], 
                         k -> Equated[i][j][k] * A[(i - 1) * DimP + k]));
        od;
        x := Concatenation (x, y);
    od;
    
    C := x^-1;
    
    # write down the generators of G wrt to new basis
    return [List (GeneratorsFlag (G), z -> z^C), C];

end;

# are the matrices in X composed of k x k blocks which differ
# only by scalars? if so, return the decomposition

AreProportional := function (X, k)
    local d, zero, Nmr, TensorFactors, x, first, Component, 
          row, i, j, C, B, First, Second, entryB, index, entryC, alpha;

    if IsMat (X) then X := [X]; fi;
    
    d := Length (X[1]);
    zero := X[1][1][1] * 0;
    Nmr := QuoInt (d, k);
    
    TensorFactors := [];
    
    for x in X do
        
        first := true;
        Component := [];
        
        for i in [1..Nmr] do
            row := [];
            for j in [1..Nmr] do
                C := x{[(i - 1) * k + 1..i * k]}{[(j - 1) * k + 1..j * k]};
                if C <> 0 * C and first then
                    B := Copy (C);
                    First := Copy (C);
                    entryB := FirstNonZeroEntry (B);
                    index := entryB[2];
                    entryB := entryB[1];
                    first := false;
                fi;
                
                # is block C a multiple of block B?
                if C <> 0 * C then
                    entryC := Concatenation (C)[index];
                    alpha := entryC / entryB;
                    if alpha * B <> C then return [false, false]; fi;
                else
                    alpha := zero;
                fi;
                Add (row, alpha);
            od;
            Add (Component, row);
       od;
       
       Second := Component;
       Add (TensorFactors, [First, Second]);
    od;
    
    return [true, TensorFactors];

end;

# construct the isomorphisms from one of the equated spaces
# to each of the rest

FindIsom := function (G, P, I, U, SumSpaces, E, Equated)
    local F, StartDim, Common, fixed, Mfixed, Component, 
          NewP, MfixedInv, i, P, M, Isom;
    
    F := BaseRing (G);
    
    StartDim := Dimension (P);
    
    Common := Intersection (I, E);
    
    if Length (Common) <> 0 and Length (I) <> Length (Common) then
        
        #fixed := Representative (Common);
        fixed := Common[1];
        Mfixed := BasisMatrix (P, U, fixed);
        
        # is Mfixed a basis for SumSpaces[fixed]?
        Component := ComponentsOfSum (P, SumSpaces, fixed, Mfixed);
        NewP := VectorSpace (Component, F);

        if Dimension (NewP) < StartDim then
            P := NewP;
            return P;
        fi;

        MfixedInv := Mfixed^-1;

        for i in Difference (I, E) do

            M := BasisMatrix (P, U, i);
            Component := ComponentsOfSum (P, SumSpaces, i, M);
            
            NewP := VectorSpace (Component, F);
            
            if Dimension (NewP) < StartDim then
                P := NewP;
                return P;
            fi;
            
            # compute the isomorphism matrix from space Mfixed to M
            Isom := MfixedInv * M;
            Equated[i] := Equated[fixed] * Isom;
            if not i in E then
                Add (E, i);
            fi;
        od;
    fi;

    return P;
end;

# S = SumSpaces[1]
# S is a potential flat in a projective geometry of dimension DimU
# find a point in the geometry or decide that S is not a flat

FindPoint := function (G, SumSpaces, DimU)

    local k, R, S, StartDim, A, F, Idnn, Equated, g, P, I, U, DimP, E,
          NmrElts, ExtField, flag, Status, CB, Y, NmrTries, Result, NewP;
    
    Status := false;
    CB := "undefined";

    NmrElts := 0; NmrTries := 100;

    if Dimension (SumSpaces[1]) mod DimU <> 0 then
        return [Status, CB];
    fi;
    
    # need to equate k spaces
    k := Length (SumSpaces);

    S := SumSpaces[1];
    StartDim := Dimension (S);
    
    A := SetupMatrix (SumSpaces);
    
    F := BaseRing (G);
    Idnn := IdentityMat (StartDim, F);
    Equated := List ([1..k], i -> Idnn);
    
    E := Set ([1]);
    
    InfoTensor1 ("#I Dimension of spaces after DirectSum is ", Dimension (S), "\n");
    
    repeat
        
        g := PseudoRandom (G);
        
        P := S^g;
        R := IdentifySubset (A, P, SumSpaces);
        I := R[1]; U := R[2];
        P := FindIsom (G, P, I, U, SumSpaces, E, Equated);
        
        DimP := Dimension (P);
        
        if DimP mod DimU <> 0 then return [Status, CB]; fi;
        
        if DimP < StartDim then
            S := P;
            SumSpaces := [S];
            SumSpaces := DirectSumSpaces (G, SumSpaces);
            R := FindPoint (G, SumSpaces, DimU);
            return R;
        else 
           NmrElts := NmrElts + 1;
           if NmrElts > NmrTries then
              Print ("#I  This group probably acts imprimitively\n");
              Print ("#I  Cannot decide if it also preserves a tensor product\n");
              return ["unknown", CB];
           fi;
        fi;
   
    until Length (E) = k;
    
    InfoTensor1 ("#I After setting up points in general position Dim = ", 
           Dimension (P), "\n");
    
    R := ConstructMatrices (G, S, Equated, A);
    Y := R[1]; CB := R[2];
    
    # are the generators of G proportional ?
    R := AreProportional (Y, DimP);
    Status := R[1];
    InfoTensor1 ("#I Are matrices proportional for dim ", DimP, "? ", Status, "\n");
    
    if Status = false then CB := "undefined"; fi;
    if Status then CB := [CB, DimP]; fi;
    if Status or (Dimension (P) = DimU) then
        return [Status, CB];
    fi;
    
    R := InvestigateMatrices (Y, P, DimU, NmrTries, F);
    Result := R[1];
    NewP := R[2];
    InfoTensor1 ("#I Result of InvestigateMatrices is ", Result, "\n");
    
    # did we construct a possible new flat ?
    if IsVectorSpace (NewP) then
        InfoTensor1 ("#I Found a posible new flat of dimension ", 
                Dimension (NewP), "\n");
        if Dimension (NewP) mod DimU <> 0 then
            Status := false;
            return [Status, CB];
        fi;
        S := NewP;
        SumSpaces := [S];
        SumSpaces := DirectSumSpaces (G, SumSpaces);
        R := FindPoint (G, SumSpaces, DimU);
        return R;
    elif Result = false then
        ExtField := Size (F)^QuoInt (DimP, DimU) - 1;
        R := ProjectivitiesGenerateField (NewP, ExtField);
        flag := R[1]; 
        if flag = true then 
            Print ("#I ** Probably found tensor decomposition as ", 
                   DimP, " x ", DimensionFlag (G) / DimP, 
                   " over GF(", ExtField + 1, ") **\n");
            Status := "unknown";
            CB := NewP;
            return [Status, CB];
        else
            # I don't really know what to do at this point --
            # I expect we will never hit this
            Error ("** CRISIS -- projectivities do not generate field **\n");
        fi;
    fi;

    return [Status, CB];
    
end;

# S is a potential flat in a projective geometry of dimension DimU
# find a point in the geometry or decide that S is not a flat

GeneralFindPoint := function (G, S, DimU)
    local R, SumSpaces;

    SumSpaces := [S];
    
    SumSpaces := DirectSumSpaces (G, SumSpaces);
    
    R := FindPoint (G, SumSpaces, DimU);
    
    return R;
end;

# Decide if S is a point; Status is set to true if
# we verify S is a point

IsPoint := function (G, S, DimU, Status, CB)
    local R, SumSpaces, k, StartDim, A, F, Idnn, Equated, E,
          NmrElts, NmrTries, g, P, I, U, DimP, Y;
    
    Status := false;

    SumSpaces := [S];
    SumSpaces := DirectSumSpaces (G, SumSpaces);
    
    if Dimension (SumSpaces[1]) mod DimU <> 0 then
        return [Status, CB];
    fi;
    
    # need to equate k spaces
    k := Length (SumSpaces);
    
    S := SumSpaces[1];
    StartDim := Dimension (S);
    
    A := SetupMatrix (S);
    
    F := BaseRing (G);
    
    Idnn := IdentityMat (StartDim, F);
    Equated := List ([1..k], i -> Idnn);
    
    E := Set ([1]);
    
    NmrElts := 0; NmrTries := 100;

    repeat
        
        g := PseudoRandom (G);
        
        P := S^g;
        I := IdentifySubset (A, P, SumSpaces);
        U := I[2];
        I := I[1];
        
        P := FindIsom (G, P, I, U, SumSpaces, E, Equated);
        DimP := Dimension (P);
        if DimP < DimU then return [Status, CB]; fi;

        NmrElts := NmrElts + 1;
        if NmrElts > NmrTries then
           Print ("#I  This group probably acts imprimitively\n");
           Error ("#I  Cannot decide if it also preserves a tensor product\n");
        fi;
        
    until Length (E) = k;
    
    Y := ConstructMatrices (G, S, Equated, A);
    
    # are the generators of G of the correct shape ?
    R := AreProportional (Y, DimU);
    Status := R[1];
    Print ("#I Are matrices proportional? ", Status, "\n");
    
end;
