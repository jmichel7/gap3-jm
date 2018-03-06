###############################################################################
##
##  grim.g  (Version 1.0)   GRIM Library        Robert Beals
##
##
##  Copyright 1997 Robert Beals, Department of Mathematics,
##                 University of Arizona, Tucson, AZ 85721, USA
##
 
RandomIntegerComb := function(v)
# v is a list of vectors.
# RandomIntegerComb returns a linear combination of
# the vectors on the list with random integer coefficients.

    local n, i, comb;
    n := Length(v);
    comb := Random(Integers)*v[1];

    for i in [2..n] do
        comb := comb + Random(Integers)*v[i];
    od;

    return comb;
end;


HermiteNormalForm_grim := function(A,M)
# A is an m x n matrix of full column rank and M is a positive
# integer such that any row vector whose entries are multiples
# of M lies in the row lattice of A.  For example, M could be
# the absolute value of an n x n nonsingular subdeterminant of A.
#
# This algorithm is from  A. Schrijver: Theory of Linear and Integer 
# Programming, Wiley, New York, 1986.
# jm 4/2016:  renamed from HermiteNormalForm to solve conflict with cryst

    local C, D, m, n, i, j, k, tmp, newrows, x, y, q, d, MyEuclideanRemainder;

    MyEuclideanRemainder := function(x,y)
       local r;    
       r := EuclideanRemainder(x,y);
       if r < 0 then r := r + AbsInt(y); fi;
       return r;
    end;


    m := Length(A);
    n := Length(A[1]);

    C := [];
    D := List(A, row -> List(row, entry -> MyEuclideanRemainder(entry, M)));
    A := [];

    for k in [0..n-1] do
        D[m+1] := List([k+1..n],i->0);
        D[m+1][1] := M;

        i := 1;
        while D[i][1] = 0 do
            i := i + 1;
        od;

        if i > 1 then
            tmp := D[i];
            D[i] := D[1];
            D[1] := tmp;
        fi;

        i := i + 1;
        while i < m + 2 do
            if D[i][1] = 0 then
                i := i + 1;
            else
                x := D[1][1];
                y := D[i][1];
                C[1] := GcdRepresentation(x,y);
                d := C[1]*[x,y];
                C[2] := [-y/d,x/d];
                newrows := List(C*[D[1],D[i]],
                                row -> List(row,
                                            entry->MyEuclideanRemainder(
                                                         entry,M)));
                if newrows[1][1] = 0 then newrows[1][1] := M; fi;

                D[1] := newrows[1];
                D[i] := newrows[2];
                i := i + 1;
            fi;
        od;

        A[k+1] := List([1..k],i->0);
        Append(A[k+1],D[1]);
        D := D{[2..m+1]}{[2..n-k]};
    od;

    for j in [2..n] do
        for i in [1..j-1] do
            q := EuclideanQuotient(A[i][j],A[j][j]);

            if A[i][j] < 0 then q := q - 1; fi;

            if q <> 0 then A[i] := A[i] - q*A[j]; fi;
        od;
    od;

    return A;
end;
                                              



CycPolIndices:=[1,2];      # CycPolIndices is a list of all i,
CycPolsDegreeBound := 1;   # in increasing order of phi(i),
                           # such that phi(i) is at most
                           # CycPolsDegreeBound.



CycPolsWithSmallDegrees:= function(n)
# returns a list of i, in increasing order of phi(i),
# including at least all i such that phi(i) is at most n.

    local phi_m,m,phi_i,i,j,indices,philess;

    if n <= CycPolsDegreeBound then
	return CycPolIndices;
    fi;

    philess:= function(i,j) return Phi(i) < Phi(j); end;
    m := 1;
    phi_m := 1;
    i := 1;

    while phi_m < n do
        m := m * Primes[i];
        phi_m := phi_m * (Primes[i] -1);
        i := i + 1;
    od;

    j := 1;
    indices := [];
    for i in [CycPolsDegreeBound+2..m] do
        phi_i := Phi(i);
        if phi_i <= n and phi_i > CycPolsDegreeBound then
            indices[j] := i;
            j := j + 1;
        fi;
    od;

    Sort(indices,philess);
    CycPolsDegreeBound := n;
    Append(CycPolIndices, indices);
    return CycPolIndices;
end;

PositiveDefinite:= function(A)
# A is a symmetric n by n rational matrix.
# PositiveDefinite returns true if A is positive definite,
# false otherwise.

    local n, pvt, i, j;

    A := Copy(A);
    n := Length(A);
    if n <> Length(A[1])  then
        Error("PositiveDefinite: <mat> must be a square matrix.");
    fi;

    if Field(Flat(A)) <> Rationals then
        Error("PositiveDefinite: <mat> must have rational entries.");
    fi;

    for i in [1..n-1] do
        for j in [i+1,n] do
            if A[i][j] <> A[j][i] then
                Error("PositiveDefinite: <mat> must be symmetric.");
            fi;
        od;
    od;

    pvt := 1;
    i := 1;


    while i <= n and pvt > 0 do
        pvt := A[i][i];
        if pvt > 0 and i < n then
            for j in [i+1..n] do
                if A[j][i] <> 0 then
                    A[j] := A[j] - (A[j][i]/pvt)*A[i];
                fi;
            od;
        fi;
        i := i + 1;
    od;
    return pvt>0;
end;

MinPolyMonteCarlo := function( A )
# This returns the minimal polynomial of the rational matrix A.
# This is a Monte Carlo algorithm; it is guaranteed to return
# a divisor of the minimal polynomial of A.  This is faster, in
# some cases, than the deterministic version, from which this is 
# adapted.

    local   F,  L,  d,  p,  M,  N,  one,  R,  h,  v,  w,  j, range;

    # try to find the field of <A>
    F := Field( Flat(A) );

    if  F <> Rationals then
    Error("MinPolyMonteCarlo: <mat> must have rational entries.");
    fi;

    # set up the other variables
    d := Length( A );
    range := [-d .. d];
    N := F.zero * A[1];
    M := ShallowCopy( N );
    Add( M, N[1] );


            # clear right and left sides <R>
            L := [];
            R := [];

            # start random vector
            v := ShallowCopy( N );
            v := List(v,function(x) return RandomList(range); end);


            # <j>-1 gives the power of <A> we are looking at
            j := 1;

            # spin vector around and construct polynomial
            repeat

                # compute the head of <v>
                h := 1;
                while h <= d and v[h] = F.zero  do
                    h := h + 1;
                od;

              # start with appropriate polynomial x^(<j>-1)
                p := ShallowCopy( M );
                p[j] := F.one;

                # divide by known left sides
                w := v;
                while h <= d and IsBound( L[h] )  do
                  p := p - w[h] * R[h];
                    w := w - w[h] * L[h];
                    while h <= d and w[h] = F.zero  do
                        h := h + 1;
                    od;
                od;

                # if <v> is not the zero vector try next power
                if h <= d  then
                    R[h] := p * w[h]^-1;
                    L[h] := w * w[h]^-1;                   
                    j := j + 1;
                    v := v * A;

                fi;
            until h > d;


    return Polynomial(F,p);

end;

OrderRatMat := function(A)
# Returns the order of the rational matrix A.
# A has finite order iff the minimal polynomial
# of A is a product of distinct cyclotomic polynomials.
# This function uses MinPolyMonteCarlo, so it is possible 
# (though unlikely) that a proper divisor of the order of A 
# will be returned, or that a finite number will be returned 
# even though A has infinite order.


    local p, q, n, i, j, div, indices, order, x, one, zero;

    n := Length(A);
    if n <> Length(A[1])  then
        Error("OrderRatMat: <mat> must be a square matrix.");
    fi;
#    if RankMat( A ) <> n  then
#        Error("OrderRatMat: <mat> must be invertible.");
#    fi;
#
## Since this is only called on group elements, we don't need
## to check invertibility.

    if Field([A[1][1]]) <> Rationals then
        Error("OrderRatMat: <mat> must have rational entries.");
    fi;

    if AbsInt(TraceMat(A)) > n then
        return "infinity";
    fi;

    x := Indeterminate(Rationals);
    one := x^0;
    zero := 0 * x^0;
    p := MinPolyMonteCarlo(A);



    if Ring(p.coefficients).name = "Integers"
       and (Degree(p) = 0 or AbsInt(p.coefficients[Degree(p)]) <= Degree(p))
       and AbsInt(p.coefficients[1]) = 1
       and Gcd(p,Derivative(p)) = one then
        indices := CycPolsWithSmallDegrees(n);
        i := 1;
        order := 1;

        while i <= Length(indices) and Phi(indices[i])<=Degree(p) do

            j := indices[i];
            q := CyclotomicPolynomial(Rationals,j);
            div := QuotientRemainder(p,q);

            if div[2] = zero then

                p:= div[1];
                order := Lcm(order,j);
            fi;

            i := i + 1;
        od;

        if p = one then
            return order;
        else
            return "infinity";
        fi;
    else
        return "infinity";
    fi;
end;



SimplestRational := function(min, max)
# Uses continued fractions to find the simplest rational
# number s in the range min <= s <= max.

    local u, v, w, x, sgn, alpha, beta, tmp, a_1,b_1,a_2,b_2;

    a_1 := Numerator(min);
    b_1 := Denominator(min);
    a_2 := Numerator(max);
    b_2 := Denominator(max);


    sgn := SignInt(a_1)*SignInt(b_1);
    if sgn < SignInt(a_2)*SignInt(b_2) then return 0; fi;
    if sgn > SignInt(a_2)*SignInt(b_2) then return "undefined"; fi;

    a_1 := AbsInt(a_1);
    a_2 := AbsInt(a_2);
    b_1 := AbsInt(b_1);
    b_2 := AbsInt(b_2);

    if sgn = -1 then
        tmp := a_1;
        a_1 := a_2;
        a_2 := tmp;

        tmp := b_1;
        b_1 := b_2;
        b_2 := tmp;
    fi;

    if a_1*b_2 > a_2*b_1 then return "undefined"; fi;
    if a_1*b_2 = a_2*b_1 then return sgn*a_1/b_1; fi;


    u := 1;
    v := 0;
    w := 0;
    x := 1;

    repeat

        alpha := QuotientRemainder(a_1,b_1);
        beta := QuotientRemainder(a_2,b_2);
        if alpha[1] + 1 = beta[1] then
            return sgn*(u*beta[1] + w)/(v*beta[1]+x);
        fi;
        if alpha[1] < beta[1] then
            tmp := EuclideanQuotient(alpha[1]+1+beta[1],2);
            return sgn*(u*tmp + w)/(v*tmp + x);
        fi;
        if alpha[1] > beta[1] then
            Error("Fatal error in SimplestRational.");
        fi;

        tmp := w;
        w := u;
        u := alpha[1] * u + tmp;

        tmp := x;
        x := v;
        v := alpha[1] * v + tmp;

        a_1 := b_2;
        a_2 := b_1;
        b_1 := beta[2];
        b_2 := alpha[2];

    until b_1 = 0 or b_2 = 0;
    return sgn*u/v;
end;


        



InvariantLattice := function(G)
# G is a rational matrix group.
#
# InvariantLattice finds a lattice L (represented by a basis) which is
# G-invariant.  That is, L*G*L^(-1) is an integer matrix group.
# 
# The lattice L is returned, and stored in G.invariantLattice.
#
# The group L*G*L^(-1) is stored in G.integerMatrixGroup.
#
# An L is calculated unless G has elements of noninteger trace,
# in which case no such L exists and false is returned.
#
# This algorithm is Las Vegas: randomness is used,
# but it only affects the running time, not the correctness
# of the output.

    local gens, A, L, n, M, v, w, X, i, m, det, prev, rand_gens, Y, Z, l;


    if IsBound(G.invariantLattice) then

        return G.invariantLattice;
    fi;

    gens := ShallowCopy(G.generators);
    rand_gens := ShallowCopy(G.generators);
    l := Length(rand_gens);
    n := Length(gens[1]);
    L := IdentityMat(n);

    while Maximum(List(Flat(gens), Denominator)) > 1 do

        A := IdentityMat(n);
        M := 1;
        det := 1;



        repeat

            X := RandomList(rand_gens);
            Y := RandomList(rand_gens);
            Z := X * Y;
            l := l + 1;
            rand_gens[l] := Z;
            if Denominator(TraceMat(Z)) > 1 then
                G.invariantLattice := false;

                return false;
            fi;


            prev := det;

            for i in [n+1, n+n] do
                w := RandomIntegerComb(A);
                X := RandomIntegerComb(gens);
                v := w * X;

                m := Lcm(List(v,Denominator));
                A := m*A;
                M := m*M;
                v := m*v;

                Add(A,v);
            od;
    
            A := HermiteNormalForm_grim(A,M);
            det := Product([1..n],i->A[i][i]);
            M := Minimum(M, det);


        until det = prev;

        L := A*L;
        gens := List(gens, M->A*M*A^(-1));
    od;
    G.invariantLattice := L;
    G.integerMatrixGroup := Group(gens, IdentityMat(n));

    return L;
end;







IntegerMatrixGroupIsFinite := function (G)
#  G is an integer matrix group.
#
#  This function tests whether or not G is finite.
#
#  If G is finite, then a positive definite G-invariant
#  quadratic form is stored in G.quadraticForm, and
#  G.isFinite is set to true, and true is returned.
#
#  If G is not finite, then the algorithm finds an element
#  of G of infinite order.  This element is not stored or
#  returned.  G.isFinite is set to false, and false is returned.
#
#  This algorithm is Las Vegas: randomness is used, but it
#  only affects the running time, not the correctness of the
#  output.
#
#  The quadratic form B is found by approximating the average
#  of TransposedMat(A)*A for A in G, and by rounding the approximation.

local gens, s, rand_gens, n, B, N, v, maxexps, l, p, sum, x, count, X, Y, A,
      B_round, i, j, w, m, automorph, e, entries, max, min, round, ok,
      k, Last10, diff, maxnorm, err, runtime, checktime;


    if not (IsBound(G.isMatGroup) and G.isMatGroup and
        Ring(Flat(G.generators)).name = "Integers") then
        Error("IntegerMatrixGroupIsFinite: <G> must be an integer matrix group.");
    fi;

    checktime := 0;
    err := 1;
    gens := G.generators;
    s := Length(gens);
    if s=0 then return true;fi;
    rand_gens := ShallowCopy(gens);
    n := Length(gens[1]);
    B := IdentityMat(n);
    B_round := IdentityMat(n);
    Last10 := [B,B,B,B,B,B,B,B,B,B];
    N := 1;


    v := [];
    for i in [1..n] do v[i] := RandomList([-n,n]); od;


## The following could be used in the rounding process, since
## the exact average of TransposedMat(A)*A will have denominators
## which divide the product of Primes[i]^maxeps[i].
##
## In practice, this seems to slow things down.
##
#    maxexps := [];  # maxexps[i] is the highest exponenent e
#    l := 1;         # such that a group of order Primes[i]^e
#    p := 2;         # has a representation in dimension n.
#    repeat
#
#        sum := 0;
#        x := p-1;
#        while n >= x do
#            sum := sum + EuclideanQuotient(n,x);
#            x := x * p;
#        od;
#        maxexps[l] := sum;
#        l := l + 1;
#        p := Primes[l];        
#    until p > n + 1;
#



    count := 1;
    l := Length(rand_gens);

## Compute the maximum frobenius norm of the generators:
    
#    maxnorm := Maximum(List(gens, A -> Sum(List(Flat(A),x->x*x))));
#
##    If G is finite, then for any A in G we must have:
##    Sum(List(Flat(A),x->x*x)) <= maxnorm^(n^2)*n^(4*n^2 + 7)
#
#    maxnorm := maxnorm^(n^2)*n^(4*n^2+7);
#
##    This is really too big to work with, and in practice
##    an element of infinite order is found long before any
##    matrix is found which violates the norm bound.    

    if Maximum(List(gens, A -> OrderRatMat(A))) = "infinity" then
        G.isFinite := false;
        return false;
    fi;


    repeat

##############################################
# GENERATE A RANDOM MATRIX A

        A := RandomList(rand_gens);
        if RandomList([true,false]) then
            l := l + 1;
            repeat
                A := A * RandomList(rand_gens);
            until RandomList([true,false,true]);
            rand_gens[l] := A;
        fi;
##############################################

        if AbsInt(TraceMat(A)) > n then
            G.isFinite := false;
            return false;
        fi;


        ##  See comment above.  maxnorm is too big to be useful.
        ##  
        # if Sum(List(Flat(A),x->x*x)) > maxnorm then
        #     G.isFinite := false;
        #     return false;
        # fi;

        diff := TransposedMat(A) * B * A - B;
        B := 2*B + diff;
        Last10 := 2 * Sublist(Last10, [2..10]);
        Add(Last10, B);
        
        N := N+N;
        #
        # (1/N)*B is converges to the average of
        # TransposedMat(A)*A for A in G.  It seems
        # to work better to keep B an integer matrix.


        if  Runtime() > checktime then
            runtime := Runtime();



            if OrderRatMat(A) = "infinity" then
                G.isFinite := false;
                return false;



            elif 2*err*Maximum(List(Flat(diff), AbsInt)) < N then

    
                if Maximum(List(Flat(diff),AbsInt)) > 0 then
                    err := err + err;
                fi;
    
                i := 1;
                j := 1;
    
                repeat
    
                    entries := List(Last10, X -> X[i][j]);
                    max := Maximum(entries);
                    min := Minimum(entries);
                    round := SimplestRational((B[i][j] + min - max)/N,
                                              (B[i][j] - min + max)/N);
                    ok := round <> "undefined";

##  This attempts to check the rounding by seeing if the 
##  denominators are plausible.  It seems to take too long
##  to be useful.
##
#                if ok then
#                    m := Denominator(round);
#        
#                    k := 1;
#                    while Primes[k] - 1 <= n and m > 1 do
#                
#                        p := Primes[k];
#                        e := 0;
#            
#            
#                        while EuclideanRemainder(m,p) = 0 do
#        
#                            e := e + 1;
#                            m := m / p;
#                        od;
#               
#                        if e > maxexps[k] then
#                            m := 0;
#                        fi;
#                        k := k + 1;
#                   od;
#                
#                    ok := m = 1;
#                fi;

                    if ok then
                        B_round[i][j] := round;
                        B_round[j][i] := round;
                        if j < i then
                            j := j + 1;
                        else 
                            i := i + 1;
                            j := 1;
                        fi;
                    fi;
                until i > n or not ok;
    
                if ok then
    
    
                    i := 1;
                    repeat
                        A := gens[i];
                        automorph := (TransposedMat(A) * (B_round * (A * v))) =
                                     B_round * v;
                        i := i + 1;
                    until i > s or not automorph;
                
                    i := 1;
                    while automorph and i <= s do
                        A := gens[i];
                        automorph := TransposedMat(A) * B_round * A = B_round;
    
    
                        i := i + 1;
                    od;
                    
                    if PositiveDefinite(B_round) then
                        if automorph then
                            G.quadraticForm := B_round;
                            G.isFinite := true;
                            return true;

##
##  I haven't experimented enough with this to see
##  if it is useful.  It maintains that B is positive definite... 
##                        
#                    else   # This might save some space
#                        B := B_round;
#                        N := Lcm(List(Flat(B),Denominator));
#                        B := B * N;
#                        Last10 := [B,B,B,B,B,B,B,B,B,B];
                        fi;
                    fi;
                fi;
            fi;
            checktime := 11 * Runtime() - 10 * runtime;


        fi;

	count := count + 1;

    until false;
end;


#############################################################################
##
#F  MatGroupOps.IsFinite(<G>) . . . . . . .  test if a matrix group is finite
##
MatGroupOps.IsFinite := function ( G )
    if IsFinite( G.field )  then
        return true;
    elif Field(Flat(G.generators)) = Rationals then
        if G.generators=[] or 
	   Maximum(List(Flat(G.generators), Denominator)) = 1 then
            return IntegerMatrixGroupIsFinite(G);
        elif InvariantLattice(G) <> false and
             IntegerMatrixGroupIsFinite(G.integerMatrixGroup) then

            G.isFinite := true;
            G.quadraticForm := TransposedMat(G.invariantLattice)*
                               G.integerMatrixGroup.quadraticForm*
                               G.invariantLattice;
            return true;
        else
            G.isFinite := false;
            return false;
        fi;

    else
        return GroupOps.IsFinite( G );
    fi;
end;

SymmetricSquareMat := function (M)
#  M is a square matrix.
#
#  A matrix representing the action X -> TransposedMat(M)*X*M
#  on symmetric matrices X is returned.
#
#  This function, together with the following two, satisfies
#
#            VecToSymMat(SymMatToVec(X)*SymmetricSquareMat(M),n)
#              = TransposedMat(M)*X*M
#
#  for X a symmetric n by n matrix and M an n by n matrix.

local n, m, i, j, k, l, x, y, T;

    n := Length(M);
    m := n*(n+1)/2;

    T := NullMat(m,m);

    for i in [1..n] do
        x := i*(i-1)/2;
        for j in [1..i] do
            x := x + 1;
            for k in [1..n] do
                y := k*(k-1)/2;
                for l in [1..k] do
                    y := y + 1;

                    T[x][y] := M[i][k]*M[j][l] + M[j][k]*M[i][l];
                    if k = l then T[x][y] := T[x][y]/2;fi;

                od;
            od;
        od;
    od;    

    return T;
end;

SymMatToVec := function (B)
# B is a symmetric matrix.
# B is returned in the form of a vector.

    local i, j, x, v, n;

    n := Length(B);
    v := [];

    for i in [1..n] do
        x := i*(i-1)/2;
        for j in [1..i] do
            x := x + 1;
            v[x] := B[i][j];
            if i = j then
                v[x] := v[x] / 2;
            fi;
        od;
    od;

    return v;
end;

VecToSymMat := function (v, n)
# v is a vector of length n*(n-1)/2 representing
# a symmetric matrix, which is returned in matrix form.
local i, j, x, B;

    B := NullMat(n,n);

    for i in [1..n] do
        x := i*(i-1)/2;
        for j in [1..i] do
            x := x + 1;
            B[i][j] := v[x];
            if i = j then
                B[i][j] := v[x] * 2;
            else
                B[j][i] := v[x];
            fi;
        od;
    od;

    return B;
end;


IsFiniteDeterministic := function (G)
# G is an integer matrix group.
# 
# G is tested for finiteness.  If G is
# finite, then G.quadraticForm is set to 
# some G-invariant positive definite quadratic form,
# G.isFinite is set to true, and true is returned.
# Otherwise, G.isFinite is set to false and
# false is returned.  This algorithm is deterministic.

local n, m, i, j, r, M, C, Cinverse, I, L, Y, X, Z, B;


    if not (IsBound(G.isMatGroup) and G.isMatGroup and
        Ring(Flat(G.generators)).name = "Integers") then
        Error("IntegerMatrixGroupIsFinite: <G> must be an integer matrix group.");
    fi;

    n := Length(G.generators[1]);
    m := n*(n+1)/2;
    I := IdentityMat(m);

    M := Product(List(G.generators, A->(SymmetricSquareMat(A)+I)/2));

    C := NullspaceMat(M-I);
    r := Length(C);
    i := 1;
    j := 1;

    if r = 0 then
        G.isFinite := false;
        return false;
    fi;

    while j <= m do
        if not IsBound(C[i]) or C[i][j] = 0 then
            Add(C,I[j]);
        else
            i := i + 1;
        fi;
        j := j + 1;
    od;

    Cinverse := C^-1;

    L := C*M*Cinverse;
    Y := L{[r+1..m]}{[r+1..m]};
    X := L{[r+1..m]}{[1..r]};
    Z := ((IdentityMat(m-r) - Y)^-1)*X;

    L{[r+1..m]}{[r+1..m]} := NullMat(m-r,m-r);
    L{[r+1..m]}{[1..r]} := Z;

    B := VecToSymMat(((SymMatToVec(IdentityMat(n))*Cinverse)*L)*C,n);

    if PositiveDefinite(B) then
        G.isFinite := true;
        G.quadraticForm := B;
        return true;
    else
        G.isFinite := false;
        return false;
    fi;

end;
