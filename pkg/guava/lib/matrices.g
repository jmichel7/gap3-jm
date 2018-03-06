#############################################################################
##
#A  matrices.g              GUAVA library                       Reinald Baart
#A                                                        &Jasper Cramwinckel
#A                                                           &Erik Roijackers
##
##  This file contains functions for generating matrices
##
#H  $Log: matrices.g,v $
#H  Revision 1.2  1997/01/20 15:06:18  werner
#H  Upgrade from Guava 1.2 to Guava 1.3 for GAP release 3.4.4.
#H
#H  Revision 1.8  1994/11/11  15:57:39  jcramwin
#H  changed the way IsStandardForm takes a few columns of a matrix
#H
#H  Revision 1.7  1994/10/19  10:00:55  rbaart
#H  Changes in comments
#H
#H  Revision 1.6  1994/10/13  15:19:18  rbaart
#H  Changed codeword functions
#H
#H  Revision 1.5  1994/09/28  11:50:26  jcramwin
#H  just some comments
#H
#H  Revision 1.4  1994/09/28  11:26:24  jcramwin
#H  something small in IsLatinSquare
#H
#H  Revision 1.3  1994/09/28  10:10:10  jcramwin
#H  changed some comments
#H
#H  Revision 1.2  1994/09/27  19:22:21  rbaart
#H  Layout adjusted to my humble demands
#H
#H  Revision 1.1  1994/09/27  19:17:13  rbaart
#H  Initial revision
#H
##

#############################################################################
##
#F  KrawtchoukMat( <n> [, <q>] )  . . . . . . .  matrix of Krawtchouk numbers
##
KrawtchoukMat := function(arg)
    local res,i,k,n,q;
    if Length(arg) = 1 then
        n := arg[1];
        q := 2;
    elif Length(arg) = 2 then
        n := arg[1];
        q := arg[2];
    else
        Error("usage: KrawtchoukMat( <n> [, <q>] )");
    fi;
    if not IsPrimePowerInt(q) then
        Error("q must be prime power");
    fi;
    res := NullMat(n+1,n+1);
    for k in [0..n] do
        res[1][k+1]:=1;
    od;
    for k in [0..n] do
        res[k+1][1] := Binomial(n,k)*(q-1)^k;
    od;
    for i in [2..n+1] do
        for k in [2..n+1] do
            res[k][i] := res[k][i-1] - (q-1)*res[k-1][i] - res[k-1][i-1];
        od;
    od;
    return res;
end;

#############################################################################
##
#F  GrayMat( <n> [, <F>] )  . . . . . . . . .  matrix of Gray-ordered vectors
##
##  GrayMat(n [, F]) returns a matrix in which rows a(i) have the property
##  d( a(i), a(i+1) ) = 1 and a(1) = 0.
##
GrayMat := function (arg) 
    local n, M, result, row, series, column, elem, q, elementnr, line,
          goingup; 
    if Length(arg) = 1 then
        n:=arg[1];
        elem:=Elements(GF(2));
        q := 2;
    elif Length(arg) = 2 then
        n:=arg[1];
        elem:=Elements(arg[2]);
        q := Length(elem);
    else
        Error("usage: GrayMat( <n> [, <F>] )");
    fi;
    M:=q^n;
    result:=NullMat(M,n);
    for column in [1..n] do
        goingup := true;
        row:=0;
        for series in [1..q^(column-1)] do
            for elementnr in [1..q] do
                for line in [1..q^(n-column)] do
                    row:=row+1;
                    if goingup then
                        result[row][column]:= elem[elementnr];
                    else
                        result[row][column]:= elem[q+1-elementnr];
                    fi;
                od;
            od;
            goingup:= not goingup;
        od;
    od;
    return result;
end;

#############################################################################
##
#F  SylvesterMat( <n> ) . . . . . . . . . . . . Sylvester matrix of order <n>
##
SylvesterMat := function(n)
    local result, syl;
    if n = 1 then
        return [[1]];
    elif Mod(n,2)=0 then
        syl:=SylvesterMat(n/2);
        result:=List(syl, x->Concatenation(x,x));
        Append(result,List(syl,x->Concatenation(x,-x)));
        return result;
    else
        Error("n must be a power of 2");
    fi;
end;

#############################################################################
##
#F  HadamardMat( <n> )  . . . . . . . . . . . .  Hadamard matrix of order <n>
##
HadamardMat := function(n)
    local result, had, i, j;
    if n = 1 then
        return [[1]];
    elif (n=2) or (n=4) or (Mod(n,8)=0) then
        had:=HadamardMat(n/2);
        result:=List(had, x->Concatenation(x,x));
        Append(result,List(had,x->Concatenation(x,-x)));
        return result;
    elif IsPrimeInt(n-1) and Mod(n,4)=0 then
        result := NullMat(n,n)+1;
        for i in [2..n] do
            result[i][i]:=-1;
            for j in [i+1..n] do
                result[i][j]:=Legendre(j-i,n-1);
                result[j][i]:=-result[i][j];
            od;
        od;
        return result;
    elif Mod(n,4)=0 then
        Error("The Hadamard matrix of order ",n," is not yet implemented");
    else
        Error("The Hadamard matrix of order ",n," does not exist");
    fi;
end;

#############################################################################
##
#F  IsLatinSquare( <M> )  . . . . . . .  determines if matrix is latin square
##
##  IsLatinSquare determines if M is a latin square, that is a q*q array whose
##  entries are from a set of q distinct symbols such that each row and each
##  column of the array contains each symbol exactly once
##
IsLatinSquare := function(M)
    local i, j, s, n, isLS, MT;
    if IsMat(M) then
        n:=Length(M);
        s:=Set(M[1]);
        isLS:= (Length(s) = n);
    else
        isLS := false;
    fi;
    i:=2;
    if isLS then
        MT:=TransposedMat(M);
    fi;
    while isLS and i<=n do
        isLS:= (Set(M[i]) = s);
        i:=i+1;
    od;
    i := 1;
    while isLS and i<=n do
        isLS:= (Set(MT[i]) = s);
        i:=i+1;
    od;
    return isLS;
end;

#############################################################################
##
#F  AreMOLS( <matlist> )  . . . . . . . . .  determines if arguments are MOLS
##
##  AreMOLS(M1, M2, ...) determines if the arguments are mutually orthogonal
##  latin squares.
##
AreMOLS := function(arg)
    local i, j, s, M, n, q2, first, second, max, fast;
    if Length(arg) = 1 then
        M:=arg[1];
    else
        M:=List([1..Length(arg)],i->arg[i]);
    fi;
    n:=Length(M);
    if ( n >= Length(M[1]) ) or not ForAll(M, i-> IsLatinSquare(i)) then
        return false; #this is right
    fi;
    q2 := Length(M[1])^2;
    max := Maximum(M[1][1]) + 1;
    M := List(M, i -> Flat(i));
    fast := (DefaultField(Flat(M)) = Rationals);
    first := 1;
    repeat
        second := first+1;
        if fast then
            repeat
                s := Set( M[first] * max + M[second] );
                second := second + 1;
            until (Length(s) < q2) or (second > n);
        else
            repeat
                s:=Set([]);
                for i in [1 .. q2] do
                    AddSet(s, [  M[first][i], M[second][i] ]);
                od;
                second := second + 1;
            until (Length(s) < q2) or (second > n);
        fi;
        first:=first + 1;
    until (Length(s) < q2) or (first >= n);
    return Length(s) = q2;
end;

#############################################################################
##
#F  MOLS( <q> [, <n>] ) . . . . . . . . . .  list of <n> MOLS of size <q>*<q>
##
##  MOLS( q [, n]) returns a list of n Mutually Orthogonal Latin
##  Squares of size q * q. If n is omitted, MOLS will return a list
##  of two MOLS. If it is not possible to return n MOLS of size q,
##  MOLS will return a boolean false.
##
MOLS := function(arg)
    local facs, res, Merged, Squares, nr, S, ToInt, q, n;

    ToInt := function(M)
        local res, els, q, i, j;
        q:=Length(M);
        els:=Elements(GF(q));
        res :=NullMat(q,q)+1;
        for i in [1..q] do
            for j in [1..q] do
                while els[res[i][j]] <> M[i][j] do
                    res[i][j]:=res[i][j]+1;
                od;
            od;
        od;
        return res-1;
    end;

    Squares := function(q, n)
        local els, res, i, j, k;
        els:=Elements(GF(q));
        res:=List([1..n], x-> NullMat(q,q,GF(q)));
        for i in [1..q] do
            for j in [1..q] do
                for k in [1..n] do
                    res[k][i][j] := els[i] + els[k+1] * els[j];
                od;
            od;
        od;
        return List([1..n],x -> ToInt(res[x]));
    end;

    Merged := function(A, B)
        local i, j, q1, q2, res;
        q1:=Length(A);
        q2:=Length(B);
        res:=KroneckerProduct(A,NullMat(q2,q2)+1);
        for i in [1 .. q1*q2] do
            for j in [1 .. q1*q2] do
                res[i][j]:= res[i][j] + q1 *
                            B[((i-1) mod q2)+1][((j-1) mod q2)+1];
            od;
        od;
        return res;
    end;

    if Length(arg) < 1 or Length(arg) > 2 then
        Error("usage: MOLS( <size> [, <# of MOLS>]");
    elif (arg[1] < 3) or (arg[1] = 6) or (arg[1] mod 4) = 2 then
        return false; #this must be so
    elif Length(arg) = 2 and (arg[2] <> 2) then
        q:=arg[1];
        n:=arg[2];
        if (not IsPrimePowerInt(q)) or (n >= q) then
            return false; #this is right
        elif IsPrimeInt(q) then
            return List([1..n],i -> List([0..q-1], y -> List([0..q-1], x -> 
                           (x+i*y) mod q)));
        else
            return Squares(q,n);
        fi;
    else
        q:=arg[1];
        res:=[[[0]],[[0]]];
        facs:=Collected(Factors(q));
        for nr in facs do
            if nr[2] = 1 then 
                S:= List([1..2], i -> List([0..nr[1]-1],y ->
                            List([0..nr[1]-1], x -> (x+i*y) mod nr[1])));
            else
                S:=Squares(nr[1]^nr[2],2);
            fi;
            res:=[Merged(res[1],S[1]),Merged(res[2],S[2])];
        od;
        return res;
    fi;
end;

#############################################################################
##
#F  VerticalConversionFieldMat( <M> ) . . . . . . .  converts matrix to GF(q)
##
##  VerticalConversionFieldMat (M) converts a matrix over GF(q^m) to a matrix
##  over GF(q) with vertical orientation of the tuples
##
VerticalConversionFieldMat := function(arg)
    local res, M, F, q, Fq, m, n, r, ConvTable, x, temp, i, j, k, zero;
    if Length(arg) = 1 then
    	M := arg[1]; 
    	F := DefaultField(Flat(M));
    elif Length(arg) = 2 then
        M := arg[1];
        F := arg[2];
    else Error("usage: VerticalConversionFieldMat( <M> [, <F> ]");
    fi;
    q := F.char;
    Fq := GF(q);
    zero := Fq.zero;
    m := F.degree;
    n := Length(M[1]);
    r := Length(M);

    ConvTable := [];
    x := Indeterminate(Fq);
    temp := Polynomial(Fq, MinPol(Z(q^m)));
    for i in [1.. q^m - 1] do
        ConvTable[i] := VectorCodeword(x^(i-1) mod temp, m+1);
    od;

    res := NullMat(r * m, n, Fq);
    for i in [1..r] do
        for j in [1..n] do
            if M[i][j] <> zero then
                temp := ConvTable[LogFFE(M[i][j]) + 1];
                for k in [1..m] do
                    res[(i-1)*m + k][j] := temp[k];
                od;
            fi;
        od;
    od;
    return res;
end;

#############################################################################
##
#F  HorizontalConversionFieldMat( <M>, <F> )  . . .  converts matrix to GF(q)
##
##  HorizontalConversionFieldMat (M, F) converts a matrix over GF(q^m) to a
##  matrix over GF(q) with horizontal orientation of the tuples
##
HorizontalConversionFieldMat := function (arg)
    local M, F, res, vec, k, n, coord, i, p, q, m, zero, g, Nul, ConvTable,
          x, F;
    if Length(arg) = 1 then
    	M := arg[1]; 
    	F := DefaultField(Flat(M));
    elif Length(arg) = 2 then
        M := arg[1];
        F := arg[2];
    else Error("usage: HorizontalConversionFieldMat( <M> [, <F> ]");
    fi;
    q := F.char;
    m := F.degree;
    zero := F.zero;
    g := Polynomial(GF(q), MinPol(F.root));
    Nul := List([1..m], i -> zero);
    ConvTable := [];
    x := Indeterminate(GF(q));
    for i in [1..Size(F) - 1] do
        ConvTable[i] := VectorCodeword(x^(i-1) mod g, m);
    od;
    res := [];
    n := Length(M[1]);
    k := Length(M);
    for vec in [0..k-1] do
        for i in [1..m] do res[m*vec+i] := []; od;
        for coord in [1..n] do
            if M[vec+1][coord] <> zero then
                p := LogFFE(M[vec+1][coord]);
                for i in [1..m] do
                    if (p+i) mod q^m = 0 then p := p+1; fi;
                    Append(res[m*vec+i], ConvTable[(p + i) mod (q^m)]);
                od;
            else
                for i in [1..m] do
                    Append(res[m*vec+i], Nul);
                od;
            fi;
        od;
    od;
    return res;
end;

#############################################################################
##
#F  IsInStandardForm( <M> [, <boolean>] ) . . . . is matrix in standard form?
##
##  IsInStandardForm(M [, identityleft]) determines if M is in standard form;
##  if identityleft = false, the identitymatrix must be at the right side
##  of M; otherwise at the left side.
##
IsInStandardForm := function(arg)
    local M, l;
    if Length(arg) > 2 or Length(arg) < 1 then
        Error("usage: IsInStandardForm( <M> [, <identityleft>] )");
    fi;
    M := arg[1];
    l := Length(M);
    if Length(arg) = 2 and arg[2] = false then
        return IdentityMat(l, DefaultField(Flat(M))) =
               M{[1..l]}{[Length(M[1])-l+1..Length(M[1])]};
    else    
        return IdentityMat(l, DefaultField(Flat(M))) = M{[1..l]}{[1..l]};
    fi;
end;

#############################################################################
##
#F  PutStandardForm( <M> [, <boolean>] [, <F>] )  .  put <M> in standard form
##
##  PutStandardForm(Mat [, idleft] [, F]) puts matrix Mat in standard form,
##  the size of Mat being (n x m). If idleft is true or is omitted, the
##  the identity matrix is put left, else right. The permutation is returned.
##
PutStandardForm := function(arg)
    local Mat, n, m, idleft, F, row, i, j, h, hp, s, zero, P;
    if Length(arg) > 3 then
        Error("usage: PutStandardForm(<mat> [, <idleft>] [, <F>])");
    fi;
    Mat := arg[1];
    n := Length(Mat);   # not the word length!
    m := Length(Mat[1]);
    if Length(arg) = 1 then
        idleft := true;
        F := DefaultField(Flat(Mat));
    elif Length(arg) = 2 then
        if IsField(arg[2]) then
            F := arg[2];
            idleft := true;
        else
            idleft := arg[2];
            F := DefaultField(Flat(Mat));
        fi;
    else
        idleft := arg[2];
        F := arg[3];
    fi;
    if idleft then 
        s := 0;
    else	
        s := m-n;
    fi;
    zero := F.zero;
    P := ();
    for j in [1..n] do
        if Mat[j][j+s] =zero then
            i := j+1;
            while (i <= n) and (Mat[i][j+s] = zero) do
                i := i + 1;
            od;
            if i <= n then
                row := Mat[j];
                Mat[j] := Mat[i];
                Mat[i] := row;
            else
                h := j+s;
                while Mat[j][h] = zero do
                    h := h + 1;
                    if h > m then h := 1; fi;
                od;
                for i in [1..n] do
                    Mat[i] := Permuted(Mat[i],(j+s,h));  
                od;   
                P := P*(j+s,h);
            fi;
        fi;
        Mat[j] := Mat[j]/Mat[j][j+s];
        for i in [1..n] do
            if i <> j then
                if Mat[i][j+s] <> zero then
                    Mat[i] := Mat[i]-Mat[i][j+s]*Mat[j];
                fi;
            fi;
        od;
    od;
    return P;
end;

