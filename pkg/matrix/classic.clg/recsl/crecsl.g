#############################################################################
##
#A  crecsl.g                    GAP library                      Frank Celler
##
#H  @(#)$Id: crecsl.g,v 1.1 1997/03/10 13:49:10 gap Exp $
##
#Y  Copyright (C) 1995,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This  file  contains  functions  which  will  help  to  recognize  groups
##  containing SL(n,q) in a constructive way.
##
Revision_crecsl_g :=
    "@(#)$Id: crecsl.g,v 1.1 1997/03/10 13:49:10 gap Exp $";


#############################################################################
##
#F  InfoRecSL?(...) . . . . . . . . . . . . . . . . . . . . . . .  infomation
##
if not IsBound(InfoRecSL1)   then InfoRecSL1  := Ignore;  fi;
if not IsBound(InfoRecSL2)   then InfoRecSL2  := Ignore;  fi;
if not IsBound(InfoRecSL3)   then InfoRecSL3  := Ignore;  fi;


#############################################################################
##

#V  CRecSL  . . . . . . . . .  functions to (constructivly) recognise SL(n,q)
##
CRecSL := rec();

CRecSL.STAT_TRANSVECTION := 1;
CRecSL.STAT_BASIS        := 2;
CRecSL.STAT_PRESERVING   := 3;
CRecSL.STAT_SPINNING     := 4;
CRecSL.STAT_SPIN_ELEMENT := 5;
CRecSL.STAT_SPANNING     := 6;
CRecSL.STAT_DIAGONAL     := 7;
CRecSL.STAT_DIAGONAL_TWO := 8;

CRecSL.STAT_NAMES := [ "finding first transvection",
  "finding transvection basis", "preserving via conjugating",
  "preserving via spinning", "elements fixing 1st vector",
  "spanning transvections", "obsolete", "obsolete" ];

CRecSL.TIME_NAMES := [ "runtime", "finding first transvection",
  "finding transvection basis", "spanning transvections" ];

CRecSL.TIME_ALL := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[1] := pos.timing[1] + R;
end;

CRecSL.TIME_FIRST_TRANS := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[2] := pos.timing[2] + R;
end;

CRecSL.TIME_TRANS_BASIS := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[3] := pos.timing[3] + R;
end;

CRecSL.TIME_SPANNING := function( pos, flag )
    local   R;
    if flag  then R := -Runtime();  else R := Runtime();  fi;
    pos.timing[4] := pos.timing[4] + R;
end;


#############################################################################
##

#F  CRecSL.BlowupVec( <V>, <v> )  . . .  write <v> as vector over prime field
##
CRecSL.BlowupVec := function( V, v )
    local   x,  c;

    # for degree one return a shallow copy
    if V.degree = 1  then
        return ShallowCopy(v);

    # store a table for small fields
    elif Size(V) <= 257  then
        if not IsBound(V.blowupTable)  then
            V.blowupTable := [Coefficients(V,V.zero)];
            for x  in [ 0 .. Size(V)-2 ] do
                V.blowupTable[x+2] := Coefficients(V,V.root^x);
            od;
        fi;
        c := [];
        for x  in v  do
            if x = V.zero  then
                Append( c, V.blowupTable[1] );
            else
                Append( c, V.blowupTable[LogFFE(x,V.root)+2] );
            fi;
        od;

    # do it using 'Coeffients'
    else
        c := [];
        for x  in v  do
            Append( c, Coefficients(V,x) );
        od;
    fi;
    return c;

end;

#############################################################################
##
#F  CRecSL.Diagonalize( <M> ) . . . . . . . . . . . . . .  diagonalize matrix
##
CRecSL.Diagonalize := function( M )
    local   n,  z,  C,  R,  e,  i,  j,  tmp,  d,  c;

    # get dimension of matrix and zero element
    n := Length(M);
    z := M[1][1] * 0;

    # store row and column operations
    C := [];
    R := [];
    e := z^0;

    # bring it into upper triangular form
    for i  in [ 1 .. n ]  do

        # clear rubbish to the left
        InfoRecSL3( "#I  clearing left part of row ", i, "\n" );
        for j  in [ 1 .. i-1 ]  do
            tmp := M[i][j];
            if tmp <> z  then
                M[i] := M[i] - tmp / M[j][j] * M[j];
                Add( R, [ i, j, tmp/M[j][j] ] );
            fi;
        od;

        # if we have a zero on the diagonal,  we are in trouble
        if M[i][i] = z  then
            d := DepthVector(M[i]);
            InfoRecSL3( "#I  adding column ", d, " to ", i, "\n" );
            c := 1/M[i][d];
            for j  in [ 1 .. n ]  do
                M[j][i] := M[j][i] + c * M[j][d];
            od;
            Add( C, [ d, i, -c ] );

        # try to get a one into the diagonal
        elif M[i][i] <> e  then
            d := i+1;
            while d <= n and M[i][d] = z  do
                d := d + 1;
            od;
            if d <= n  then
                InfoRecSL3( "#I  adding column ", d, " to ", i, "\n" );
                c := (1-M[i][i])/M[i][d];
                for j  in [ 1 .. n ]  do
                    M[j][i] := M[j][i] + c * M[j][d];
                od;
                Add( C, [ d, i, -c ] );
            fi;
        fi;
    od;

    # bring matrix into diagonal form
    for i  in [ 1 .. n-1 ]  do
        InfoRecSL3( "#I  clearing right part of row ", i, "\n" );
        for j  in [ i+1 .. n ]  do
            tmp := M[i][j];
            if tmp <> z  then
                M[i] := M[i] - tmp / M[j][j] * M[j];
                Add( R, [ i, j, tmp/M[j][j] ] );
            fi;
        od;
    od;

    # make sure that we have ones on the diagonal
    InfoRecSL3( "#I  normalizing diagonal\n" );
    i := 1;
    while i < n  do
        if M[i][i] <> e  then
            j := i+1;
            while M[j][j] = e  do
                j := j + 1;
            od;
            InfoRecSL3( "#I  normalizing [ ", i, ", ", j, " ]\n" );
            Add( R, [ i, j, -1/M[j][j] ] );
            Add( C, [ j, i, M[i][i]-e ] );
            Add( R, [ j, i, (e-M[i][i])*M[j][j] ] );
            Add( R, [ i, j, 1/(M[i][i]*M[j][j]) ] );
            M[j][j] := M[j][j] * M[i][i];
            M[i][i] := e;
            i := j;
        else
            i := i + 1;
        fi;
    od;

    # combine the elementary row and column operations
    R := Concatenation( R, Reversed(C) );

    # now combine common rows
    C := [];
    i := 1;
    n := Length(R);
    while i <= n  do
        j := i;
        z := 0 * M[1];
        while j <= n and R[i][1] = R[j][1]  do
            z[R[j][2]] := R[j][3];
            j := j + 1;
        od;
        Add( C, [ R[i][1], z ] );
        i := j;
    od;

    # and return
    return C;

end;


#############################################################################
##
#F  CRecSL.LcmRepresentation( <lst> ) . . . . . . . . . representation of lcm
##
##  if X_i  are elements in a  abelian group of  order o_i, return  a list of
##  integers n_i, such that Prod(X_i^n_i) is of order Lcm(o_i).
##
CRecSL.LcmRepresentation := function( arg )
    local   lcm,  pow,  n,  p,  c;

    # unravel args
    if Length(arg) = 1 and IsList(arg)  then
        arg := arg[1];
    fi;

    # compute the lcm and factorise it
    lcm := Set(Collected(Factors(Lcm(arg))));

    # loop over <arg>
    pow := [];
    for n  in arg  do
        p := 1;
        for c  in Collected(Factors(n))  do
            if c  in lcm  then
                p := p * c[1]^c[2];
                RemoveSet( lcm, c );
            fi;
        od;
        Add( pow, n/p );
    od;

    # and return
    return pow;

end;


#############################################################################
##

#F  CRecSL.ApplyInverseTT( <v>, <T> ) . . . . . . . . . .  compute <v>*<T>^-1
##
CRecSL.ApplyInverseTT := function( v, T )
    return v - (T[1]*v)/(T[1]*T[2]) * T[3];
end;


#############################################################################
##
#F  CRecSL.ApplyTT( <v>, <T> )  . . . . . . . . . . . . . . . compute <v>*<T>
##
CRecSL.ApplyTT := function( v, T )
    return v + (T[1]*v)/(T[1]*T[2]) * T[3];
end;


#############################################################################
##
#F  CRecSL.TT2TV( <tr> )  . . . . . . . . . . . . . . . .  inverse of 'TV2TT'
##
CRecSL.TT2TV := function( tr )
    local   e;

    # construct the identity
    e := List( tr[1], x -> tr[1] )^0;

    # apply the transvection to each vector
    return List( e, x -> CRecSL.ApplyTT( x, tr ) );

end;


#############################################################################
##
#F  CRecSL.TV2TT( <t> ) . . . . . . . . . .  convert transvection into triple
##
CRecSL.TV2TT := function( t )
    local   w, b, p;

    # compute the left nullspace of <t> - 1
    w := LeftNullspaceMat( t - t^0 );

    # compute the missing basis vector
    b := w[1] * 0;
    p := Difference( [1..Length(w[1])], List(w,DepthVector) )[1];
    b[p] := b[p]^0;

    # compute the image of <b> under <t> and return the triple
    return [ RightNullspaceMat(w)[1], b, b*t-b ];

end;


#############################################################################
##

#F  CRecSL.EnlargeBasis( <b>, <v> ) . . . . . enlarge basis <b> by vector <v>
##
CRecSL.EnlargeBasis := function( b, v )
    local   n,  d,  i,  j;

    n := Length(v);
    d := DepthVector(v);
    i := 1;
    while d <= n  do
        while i <= Length(b) and DepthVector(b[i]) < d  do
            i := i + 1;
        od;

        # if we reach the end of our basis append <v> and return
        if Length(b) < i  then
            Add( b, v*v[d]^-1 );
            return true;

        # if we found a vector with the same depth,  subtract it
        elif DepthVector(b[i]) = d  then
            v := v - v[d] * b[i];
            d := DepthVector(v);

        # otherwise add <v> at position <i>
        else
            for j  in [ Length(b)+1, Length(b) .. i+1 ]  do
                b[j] := b[j-1];
            od;
            b[i] := v*v[d]^-1;
            return true;
        fi;
    od;
    return false;
end;


#############################################################################
##
#F  CRecSL.IsIndependent( <b>, <v> )  . . . . . . check if <v> depends on <b>
##
CRecSL.IsIndependent := function( b, v )
    local   n,  d,  i;

    n := Length(v);
    d := DepthVector(v);
    i := 1;
    while d <= n  do
        while i <= Length(b) and DepthVector(b[i]) < d  do
            i := i + 1;
        od;

        # if we reach the end of our basis <v> must be independent
        if Length(b) < i  then
            return true;

        # if we found a vector with the same depth,  subtract it
        elif DepthVector(b[i]) = d  then
            v := v - v[d] * b[i];
            d := DepthVector(v);

        # otherwise add <v> must be independent
        else
            return true;
        fi;
    od;
    return false;
end;


#############################################################################
##

#F  CRecSL.Random( <pos> )  . . . . .  create a random element of <pos.group>
##
##  if  <pos.smallRandom>  is  true  construct  only small words  in the free
##  group to  find the corresponding "random"  word.  Otherwise  use our semi
##  random algorithm
##
CRecSL.Random := function( pos )
    local  i,  j;
    
    # catch trivial case
    if Length(pos.generators) = 0  then
        return pos.identity;
    fi;

    # if <pos.smallRandom> is true  construct words in the free group
    if pos.smallRandom  then
        if not IsBound(pos.randomSeed)  then

            # use the old random to construct a start
            pos.smallRandom := false;
            CRecSL.Random(pos);

            # and shake it well
            for i  in [ 1 .. 50 ]  do
                CRecSL.Random(pos);
            od;
            pos.smallRandom := true;

            # construct structures for cyclic reduced words
            pos.randomGens  := Concatenation( Reversed(List(
                                  pos.randomSeed, x->x^-1)), [0],
                                  pos.randomSeed );
            pos.randomAgens := Concatenation( Reversed(List(
                                  pos.abstractSeed, x->x^-1)), [0],
                                  pos.abstractSeed );
            pos.randomSeed  := InitCRWord(Length(pos.randomSeed));
            Unbind(pos.abstractSeed);
        fi;
        i := NextCRWord(pos.randomSeed) + (1+pos.randomSeed.rank);
        return [ Product(pos.randomGens{i}), Product(pos.randomAgens{i}) ];

    # use semi random algorithm
    else

        # if '<pos>.randomSeed' is unbound, construct a start
        if not IsBound(pos.randomSeed)  then
            pos.randomSeed   := ShallowCopy(pos.generators);
            pos.abstractSeed := ShallowCopy(pos.abstractGenerators);
            for i  in [ Length(pos.randomSeed)+1 .. 10 ]  do
                j := RandomList( [ 1 .. Length(pos.generators) ] );
                pos.randomSeed[i]   := pos.generators[j];
                pos.abstractSeed[i] := pos.abstractGenerators[j];
            od;
        fi;
    
        # get two random position and a random L/R
        i := Random( [ 1 .. Length(pos.randomSeed) ] );
        repeat
            j := Random( [ 1 .. Length(pos.randomSeed) ] );
        until i <> j;
        if Random( [ true, false ] )  then
            pos.randomSeed[j]   := pos.randomSeed[j]*pos.randomSeed[i];
            pos.abstractSeed[j] := pos.abstractSeed[j]*pos.abstractSeed[i];
        else
            pos.randomSeed[j]   := pos.randomSeed[i]*pos.randomSeed[j];
            pos.abstractSeed[j] := pos.abstractSeed[i]*pos.abstractSeed[j];
        fi;
        return [ pos.randomSeed[j], pos.abstractSeed[j] ];
    fi;
    
end;


#############################################################################
##
#F  CRecSL.Order( <pos>, <pol> )  . . . . . . . . . . . . . .  order of <pol>
##
CRecSL.Order := function( pos, pol )
    local   po;

    po := RecSL.OrderScalar( pos, pol );
    return OrderFFE(po[2]) * po[1];
end;


#############################################################################
##

#F  CRecSL.DiagonalMat( <pos> ) . . . . . . . . . . . .  find diagonal matrix
##
CRecSL.DiagonalMat := function( pos )
    local   f,  i,  a,  m,  o,  c,  b,  d,  w,  v;

    # use <f> to check that minimal polynomial contains one linear factor
    f := Indeterminate(pos.field)^(pos.q-1) - 1;
    for i  in [ 1 .. 10*pos.d ]  do

        # construct a random element
        repeat
            a := pos.operations.Random(pos);
        until a <> pos.group.identity;
        m := MinimalPolynomial( FiniteFieldMatrices, a[1] );
        m.baseRing := pos.field;

        # test if it is almost irreducible
        if pos.d <> Degree(m)  then
            InfoRecSL3( "#I  degree of minimal polynomial is too small\n" );
        elif 0 < Degree(Gcd(m,Derivative(m)))  then
            InfoRecSL3( "#I  minimal polynomial is not square free\n" );
        elif 2 < pos.d and 1 <> Degree(Gcd(m,f))  then
            InfoRecSL3( "#I  minimal polynomial has no linear factor\n" );
        elif 2 = pos.d and Degree(Gcd(m,f)) < 1  then
            InfoRecSL3( "#I  minimal polynomial has no linear factor\n" );
        elif 2 <> Length(List(Factors(m),Degree))  then
            InfoRecSL3( "#I  element is not almost irreducible\n" );

        # it is almost irreducible,  make it similar to a diagonal matrix
        else
            o := Order( FiniteFieldMatrices, a[1] );
            o := o / Gcd( o, pos.q-1 );
            a := [ a[1]^o, a[2]^o ];

            # recompute the minimal polynomial
            m := MinimalPolynomial( FiniteFieldMatrices, a[1] );
            m.baseRing := pos.field;

            # test if the (two) roots are distinct
            if 2 = Length(Set(Factors(m)))  then
                InfoRecSL2( "#I  found diagonal after ", i, " tries\n" );

                # chose root with higher power
                c := CharacteristicPolynomial( FiniteFieldMatrices, a[1] );
                m := Factors(m);
                if c mod m[1]^2 = 0*c  then
                    b := -m[1].coefficients[1];
                    d := -m[2].coefficients[1];
                else
                    b := -m[2].coefficients[1];
                    d := -m[1].coefficients[1];
                fi;

                # is there a condition?
                if true  then

                    # compute dual of eigen space
                    w := LeftNullspaceMat(a[1]-b*pos.group.identity);
                    v := RightNullspaceMat(w)[1];
                    v := v[DepthVector(v)]^-1 * v;
                    f := LeftNullspaceMat(a[1]-d*pos.group.identity)[1];

                    # and return
                    pos.statistic[CRecSL.STAT_DIAGONAL] := i;
                    a[2].name := "d1";
                    return rec( diagonal       := a,
                                dualEigenSpace := v,
                                eigenSpace     := w,
                                eigenValue     := b,
                                eigenVector    := f );
                fi;

            # sorry,  we cannot use this
            else
                InfoRecSL3( "#I  power is a scalar\n" );
            fi;
        fi;
    od;
    InfoRecSL2( "#I  failed to find diagonal after ", i, " tries\n" );
    pos.statistic[CRecSL.STAT_DIAGONAL] := -i;
    return false;

end;


#############################################################################
##
#F  CRecSL.IndependentDiagonalMat( <pos>, <diag> )  . . . second diagonal mat
##
CRecSL.IndependentDiagonalMat := function( pos, diag )
    local   i,  r,  w,  b,  es,  vs,  b1,  b2,  basis,  fac,  id,  pol;

    # don't use generators,  because we need need various ones to try
    for i  in [ 1 .. 50 ]  do

        # conjugate with a random element in order to get different ones
        r := pos.operations.Random(pos);

        # try to find a diagonal matrix with a different eigen space
        w := diag.dualEigenSpace * TransposedMat(r[1]^-1);
        w := w[DepthVector(w)]^-1 * w;

        # if the standard duals are different,  we found our second diag
        if w <> diag.dualEigenSpace  then
            b := [ diag.diagonal[1]^r[1], diag.diagonal[2]^r[2] ];

            # compute the common eigen space
            es := RightNullspaceMat( [ w, diag.dualEigenSpace ] );

            # get two missing basis vector as eigen vectors for es of dim 1
            b1 := diag.eigenVector;
            b2 := diag.eigenVector * r[1];

            # put together a common bases
            basis := Concatenation( [b1,b2], es );

            # compute representation of on the factor space < <b1>, <b2> >
            fac := [ [], [] ];
            fac[1][1] := SolutionMat( basis, b1*diag.diagonal[1] ){[1,2]};
            fac[1][2] := SolutionMat( basis, b2*diag.diagonal[1] ){[1,2]};
            fac[2][1] := SolutionMat( basis, b1*b[1] ){[1,2]};
            fac[2][2] := SolutionMat( basis, b2*b[1] ){[1,2]};

            # do a few sanity checks
            id  := fac[1]^0;
            pol := CharacteristicPolynomial(fac[1]*fac[2]).coefficients;
            if    pol[2] = pos.field.zero
               or DegreeFFE(pol[1]/pol[2]^2) < pos.k
            then
                InfoRecSL3("#I  group is definable over smaller field\n");
            elif    Comm(fac[1]^2,fac[2]^2) = id
                and Comm(fac[1]^2,(fac[1]*fac[2])^2) = id
                and Comm(fac[2]^2,(fac[1]*fac[2])^2) = id
            then
                InfoRecSL3( "#I  squares commute\n" );
            else
                InfoRecSL2("#I  found next diagonal after ",i," tries\n");
                pos.statistic[CRecSL.STAT_DIAGONAL_TWO] := i;
                b[2].name := "d2";
                return rec( diagonal1  := diag.diagonal,
                            diagonal2  := b,
                            quotient   := fac,
                            eigenSpace := es,
                            basis      := basis,
                            fixedSpace := diag.dualEigenSpace,
                            eigenValue := diag.eigenValue,
                            eigenBasis := Concatenation([b1],diag.eigenSpace)
                          );
            fi;
        fi;
    od;

    # ops,  it might not contains SL
    pos.statistic[CRecSL.STAT_DIAGONAL_TWO] := -i;
    InfoRecSL2( "#I  failed to find next diagonal after ",i," tries\n" );
    return false;

end;


#############################################################################
##
#F  CRecSL.TransvectionSmall( <pos> ) . . . . find a transvection (small fld)
##
CRecSL.TransvectionSmall := function( pos )

    local   g,          # our group
            x,          # indeterminate over the field of <g>
            q,          # size of the field
            p,          # char of the field
            i,          # number of tries so far
            a,          # current random element
            m,          # minimal polynomial of <a>
            c,          # characteristic polynomial of <a>
            e,          # root of <m>
            l,          # <x> - <e>
            zero,       # zero polynomial
            tmp;

    # if we already know a transvection,  return
    if IsBound(pos.transvection)  then return pos.transvection;  fi;

    # start the timer
    CRecSL.TIME_FIRST_TRANS(pos,true);

    # get the group
    g := pos.group;

    # start generating random elements
    x := Indeterminate(pos.field);
    q := Size(pos.field);
    p := pos.field.char;
    zero := pos.field.zero * x;
    for i  in [ 1 .. Maximum(q,21)*3 ]  do
        InfoRecSL3( "#I  searching for a transvection, ", i, ".th try\n" );

        # compute a (pseudo) random group element and its minimal polynomial
        if IsBound(pos.useTransvection)  then
            a := pos.useTransvection;
            Unbind(pos.useTransvection);
        else
            repeat
                a := pos.operations.Random(pos);
            until a[1] <> g.identity;
        fi;
        m := MinimalPolynomial( FiniteFieldMatrices, a[1] );
        m.baseRing := pos.field;

        # compute the gcd with its derivative
        tmp := Gcd( m, Derivative(m) );

        # it should be of degree one or two
        if     Degree(tmp) < 3
           and Length(Set(Factors(tmp))) = 1
           and Degree(Factors(tmp)[1]) = 1
        then
            e := - Factors(tmp)[1].coefficients[1];
            l := x - e;
            InfoRecSL3( "#I  found root ", e, "\c" );

            # check that the multiplicity is two
            tmp := QuotientRemainder( m, l^2 );
            if Value(tmp[1],e) <> pos.field.zero  then
                InfoRecSL3( " with multiplicity 2\n" );
                c := CharacteristicPolynomial( FiniteFieldMatrices, a[1] );
                c.baseRing := pos.field;
                if m = c or ( c mod l^4 ) <> zero  then
                    InfoRecSL2("#I  found transvection after ",i," tries\n");

                    # store the diagonal element
                    tmp := a[1]^p;
                    pos.almostElement := [tmp,a[2]^p,TransposedMat(tmp^-1)];

                    # compute the order (to avoid big powers)
                    tmp := pos.operations.Order( pos, m ) / p;
                    a := [ a[1]^tmp, a[2]^tmp ];
                    pos.transvection := [ CRecSL.TV2TT(a[1]), a[2] ];
                    pos.statistic[CRecSL.STAT_TRANSVECTION] := i;
                    CRecSL.TIME_FIRST_TRANS(pos,false);
                    return pos.transvection;
                fi;

            # if the Characteristic is two,  square the element
            elif p = 2  then
                InfoRecSL3( " with multiplicity > 2,  using square\n" );
                a := [ a[1]^2, a[2]^2 ];
                m := MinimalPolynomial( FiniteFieldMatrices, a[1] );
                c := CharacteristicPolynomial( FiniteFieldMatrices, a[1] );
                m.baseRing := pos.field;
                c.baseRing := pos.field;
                if m = c or ( c mod l^4 ) <> zero  then

                    # store the diagonal element
                    pos.almostElement := [ a[1]^p, a[2]^p ];
                    pos.almostElement[3] :=
                        TransposedMat( pos.almostElement[1]^-1 );

                    # compute the order (to avoid big powers)
                    InfoRecSL2("#I  found transvection after ",i," tries\n");
                    tmp := pos.operations.Order( pos, m ) / p;
                    a := [ a[1]^tmp, a[2]^tmp ];
                    pos.transvection := [ CRecSL.TV2TT(a[1]), a[2] ];
                    pos.statistic[CRecSL.STAT_TRANSVECTION] := i;
                    CRecSL.TIME_FIRST_TRANS(pos,false);
                    return pos.transvection;
                fi;
            else
                InfoRecSL3( " with multiplicity > 2\n" );
            fi;
        fi;
    od;
    pos.statistic[CRecSL.STAT_TRANSVECTION] := -i;
    InfoRecSL2("#I  failed to find a transvection after ", i, " tries\n");
    CRecSL.TIME_FIRST_TRANS(pos,false);
    return false;
end;


#############################################################################
##
#F  CRecSL.TransvectionBasis( <pos> ) . . . .  find independent transvections
##
##  this function is used for small fields
##
CRecSL.TransvectionBasis := function( pos )

    local   tt,         # one transvection
            diff,       # transvection part of <tt>
            basis,      # basis over the prime field found so far
            ws,         # normalized dual space of fixed spaces of <tt>
            found,      # transvection fixing <ws> found so far
            conj,       # next conjugate of <tt>
            ws2,        # normalized dual space of fixed spaces of <new>
            stat,       # number of tries so far
            prev,       # number of preserving transvections found so far
            tmp,
            c;

    # if we know a basis,  return
    if IsBound(pos.transvectionBasis) then return pos.transvectionBasis;  fi;

    # get the first transvection
    tt := pos.operations.TransvectionSmall(pos);
    if tt = false  then return false;  fi;
    CRecSL.TIME_TRANS_BASIS(pos,true);

    # blow up the difference
    diff := tt[1][3];
    if 1 < pos.field.degree  then
        diff := CRecSL.BlowupVec( pos.field, diff );
    fi;
    basis := [];
    CRecSL.EnlargeBasis( basis, diff );

    # check if we want more
    if pos.field.degree*(pos.d-1) = 1  then
        pos.transvectionBasis := [tt];
        CRecSL.TIME_TRANS_BASIS(pos,false);
        return pos.transvectionBasis;
    fi;

    # normalize fixed space
    ws := tt[1][1];
    ws := ws[DepthVector(ws)]^-1 * ws;

    # collect centralizing transvections in <found>
    found := [tt];

    # start conjugating around
    stat  := 0;
    prev  := 0;
    for stat  in [ 1 .. 15*Size(pos.field)*pos.field.degree*(pos.d-1) ]  do

        # conjugate transvection with a random element
        tmp  := pos.operations.Random(pos);
        conj := [ [ tt[1][1]*TransposedMat(tmp[1]^-1), 
                    tt[1][2]*tmp[1],
                    tt[1][3]*tmp[1] ],
                  tt[2]^tmp[2] ];

        # test if the transvection perserves the fixed space of <tt>
        ws2 := conj[1][1][DepthVector(conj[1][1])]^-1*conj[1][1];
        if ws2=ws or tt[1][1]*conj[1][3]=pos.field.zero  then
            InfoRecSL3( "#I  found preserving transvection\n" );
            prev := prev + 1;

            # if the fixed spaces are not equal compute the conjugate
            if ws2 <> ws  then
                c := CRecSL.ApplyInverseTT( tt[1][2], conj[1] );
                c := CRecSL.ApplyTT( c, tt[1] );
                c := CRecSL.ApplyTT( c, conj[1] ) - tt[1][2];
                tmp := [ [ tt[1][1], tt[1][2], c ], tt[2] ^ conj[2] ];

            # otherwise use transvection <conj>
            else
                c := CRecSL.ApplyTT( tt[1][2], conj[1] ) - tt[1][2];
                tmp := [ [ tt[1][1], tt[1][2], c ], conj[2] ];
            fi;

            # blow up difference
            diff := tmp[1][3];
            if 1 < pos.field.degree  then
                diff := CRecSL.BlowupVec( pos.field, diff );
            fi;

            # add it if we don't know it already
            if CRecSL.EnlargeBasis( basis, diff )  then
                Add( found, tmp );
                if Length(found) = pos.field.degree*(pos.d-1)  then
                    InfoRecSL2( "#I  found transvection basis after ",
                                stat, " (", prev, ") tries\n" );
                    pos.transvectionBasis := found;
                    pos.statistic[CRecSL.STAT_BASIS] := stat;
                    pos.statistic[CRecSL.STAT_PRESERVING] := prev;
                    CRecSL.TIME_TRANS_BASIS(pos,false);
                    return pos.transvectionBasis;
                fi;
                InfoRecSL3( "#I  transvection is indepent, found ",
                            Length(found), "\n" );
            fi;
        else
            InfoRecSL3( "#I  conjugate does not preserve space, ",
                        stat, ".th try\n" );
        fi;
    od;

    # return failure
    pos.statistic[CRecSL.STAT_BASIS] := -stat;
    pos.statistic[CRecSL.STAT_PRESERVING] := prev;
    InfoRecSL2( "#I  failed to find transvection basis after ", stat-1,
                " (", prev, ") tries\n" );
    CRecSL.TIME_TRANS_BASIS(pos,false);
    return false;

end;


#############################################################################
##
#F  CRecSL.TransvectionBasisD4( <pos> ) . find independent transvections, d>3
##
CRecSL.TransvectionBasisD4 := function( pos )
    local   tt,  diff,  basis,  ws,  found,  stab,  spin,  spvc,  tmp,  
            ud,  stat,  conj,  ws2,  c,  next,  ns,  vcs;

    # dimension must be 4
    if pos.d < 4  then
        Error( "dimension must be at least 4" );
    fi;

    # if we know a basis,  return
    if IsBound(pos.transvectionBasis) then return pos.transvectionBasis;  fi;

    # get the first transvection
    tt := pos.operations.TransvectionSmall(pos);
    if tt = false  then return false;  fi;
    CRecSL.TIME_TRANS_BASIS(pos,true);

    # blow up the difference
    diff := tt[1][3];
    if 1 < pos.field.degree  then
        diff := CRecSL.BlowupVec( pos.field, diff );
    fi;
    basis := [];
    CRecSL.EnlargeBasis( basis, diff );

    # normalize fixed space
    ws := tt[1][1];
    ws := ws[DepthVector(ws)]^-1 * ws;

    # collect centralizing transvections in <found>
    found := [tt];

    # collect conjugating elements where the conjugate fixes first vector
    stab := [];

    # remember next vector to spin around in <spin>
    spin := 1;
    spvc := pos.almostElement;

    # find dual of the image of the transvection
    tmp := RightNullspaceMat([tt[1][1]]);
    tmp := Filtered( tmp, x -> DepthVector(x) <> DepthVector(tt[1][3]) );
    Add( tmp, tt[1][2] );
    ud := RightNullspaceMat(tmp);

    # start conjugating around
    for stat  in [ 1 .. (pos.q+pos.d*pos.k)*3 ]  do

        # if we can spin around do it
        if spin < Length(found)  then
            InfoRecSL3( "#I  using spinned instead of conjugate\n" );
            spin := spin + 1;
            tmp  := found[spin];
            conj := [ [ tmp[1][1] * spvc[3], tmp[1][2] * spvc[1],
                        tmp[1][3] * spvc[1] ], tmp[2] ^ spvc[2] ];
            ws2 := conj[1][1][DepthVector(conj[1][1])]^-1*conj[1][1];
            if ws2 <> ws  then
                Error( "this should not happen" );
            fi;
            c := CRecSL.ApplyTT( tt[1][2], conj[1] ) - tt[1][2];
            conj := [ [ tt[1][1], tt[1][2], c ], conj[2] ];
            pos.statistic[CRecSL.STAT_SPINNING] := 
              pos.statistic[CRecSL.STAT_SPINNING] + 1;
            

        # maybe we have a conjugate stabilising the first vector
        elif 2 < Length(found) and 0 < Length(stab)  then
            pos.statistic[CRecSL.STAT_SPIN_ELEMENT] := 
              pos.statistic[CRecSL.STAT_SPIN_ELEMENT] + 1;

            # try next element in <stab>
            next := stab[Length(stab)];
            Unbind(stab[Length(stab)]);

            # map image vectors of the transvections
            tmp  := List( found{[2..Length(found)]}, x -> CRecSL.BlowupVec(
                         pos.field, [tt[1][1]*(x[1][3]*next[1])] ) );

            # get a basis for the solution
            ns := LeftNullspaceMat(tmp);

            # and check that they don't fix the image of the 1st transvection
            if ns <> false  then
                vcs  := List( found{[2..Length(found)]}, x -> x[1][3] );
                tmp := false;
                for c  in ns  do
                    if tmp = false and ud*(c*vcs) <> pos.field.zero  then
                        tmp := c;
                    fi;
                od;
                if tmp <> false  then
                    InfoRecSL3( "#I  found another spinning element\n" );
                    tmp  := List( tmp, Int );
                    spvc := [ [ tt[1][1]*TransposedMat(next[1]^-1),
                                tt[1][2]*next[1], Sum( [2..Length(found)],
                              x -> tmp[x-1]*found[x][1][3] )*next[1] ],
                              Product( [2..Length(found)], x ->
                              found[x][2]^tmp[x-1] ) ^ next[2] ];
                    spvc[1] := CRecSL.TT2TV(spvc[1]);
                    spvc[3] := TransposedMat(spvc[1]^-1);
                    spin := 1;
                fi;
            fi;
            conj := false;

        # conjugate transvection with a random element
        else
            tmp  := pos.operations.Random(pos);
            conj := [ [ tt[1][1]*TransposedMat(tmp[1]^-1), 
                        tt[1][2]*tmp[1],
                        tt[1][3]*tmp[1] ],
                      tt[2]^tmp[2] ];

            # maybe we found a transvection stabilising the first basis vec
            if tt[1][2] * conj[1][1] = pos.field.zero  then
                InfoRecSL3("#I  found conjugate fixing first vector\n");
                Add( stab, tmp );
            fi;

            # test if the transvection perserves the fixed space of <tt>
            ws2 := conj[1][1][DepthVector(conj[1][1])]^-1*conj[1][1];
            if ws2=ws or tt[1][1]*conj[1][3]=pos.field.zero  then
                InfoRecSL3( "#I  found preserving transvection\n" );
                pos.statistic[CRecSL.STAT_PRESERVING] := 
                  pos.statistic[CRecSL.STAT_PRESERVING] + 1;

                # if the fixed spaces are not equal compute the conjugate
                if ws2 <> ws  then
                    c := CRecSL.ApplyInverseTT( tt[1][2], conj[1] );
                    c := CRecSL.ApplyTT( c, tt[1] );
                    c := CRecSL.ApplyTT( c, conj[1] ) - tt[1][2];
                    conj := [ [ tt[1][1], tt[1][2], c ], tt[2] ^ conj[2] ];

                # otherwise use transvection <conj>
                else
                    c := CRecSL.ApplyTT( tt[1][2], conj[1] ) - tt[1][2];
                    conj := [ [ tt[1][1], tt[1][2], c ], conj[2] ];
                fi;

            # sorry, it doesn't preserve the fix space
            else
                InfoRecSL3( "#I  conjugate does not preserve space, ",
                            stat, ".th try\n" );
                conj := false;
            fi;
        fi;

        # if we found a conjugate use it
        if conj <> false  then

            # blow up difference
            diff := conj[1][3];
            if 1 < pos.field.degree  then
                diff := CRecSL.BlowupVec( pos.field, diff );
            fi;

            # add it if we don't know it already
            if CRecSL.EnlargeBasis( basis, diff )  then
                Add( found, conj );
                if Length(found) = pos.field.degree*(pos.d-1)  then
                    tmp := pos.statistic[CRecSL.STAT_PRESERVING];
                    InfoRecSL2( "#I  found transvection basis after ",
                                stat, " (", tmp, ") tries\n" );
                    pos.transvectionBasis := found;
                    pos.statistic[CRecSL.STAT_BASIS] := stat;
                    CRecSL.TIME_TRANS_BASIS(pos,false);
                    return pos.transvectionBasis;
                fi;
                InfoRecSL3( "#I  transvection is indepent, found ",
                            Length(found), "\n" );
            fi;
        fi;
    od;

    # return failure
    pos.statistic[CRecSL.STAT_BASIS] := -stat;
    InfoRecSL2( "#I  failed to find transvection basis after ", stat-1,
                " (", pos.statistic[CRecSL.STAT_PRESERVING], ") tries\n" );
    CRecSL.TIME_TRANS_BASIS(pos,false);
    return false;

end;


#############################################################################
##
#F  CRecSL.TransvectionBasisLarge( <pos> )  .  transvection basis (large fld)
##
CRecSL.TransvectionBasisLarge := function( pos )
    local   a,  tt,  ind,  i,  j,  b,  h,  qos,  x,  f,  L,  z,  m,  
            r,  e1,  tr,  st,  diff;

    # if we already know a transvection,  return
    if IsBound(pos.transvectionBasis) then return pos.transvectionBasis; fi;

    # try to find a diagonal matrix, if this fails there is no hope
    a := pos.operations.DiagonalMat(pos);
    if a = false  then return false;  fi;

    # store independent transvection triples in <tt>
    tt := [];

    # store basis computed so far in <ind>
    ind := [];

    # this diagonal determines the fixed space,  but we need a second
    for i  in [ 1 .. 3*pos.d ]  do

        # find a second one which generators GL(2,q)
        j := 0;
        repeat
            j := j + 1;

            # if we cannot find a second one  give up
            b := pos.operations.IndependentDiagonalMat( pos, a );
            if b = false  then return false;  fi;

            # lets hope that they generate SL
            h := Group( b.quotient[1], b.quotient[2] );
            h.field := pos.field;

            # check if <h> contains SL,  use small random algorithm
            qos := CRecSL.Setup( h, b.quotient );
            qos.smallRandom := true;
            qos := CRecognizeSL(qos);

        until j = 5 or qos.containsSL;

        # if we failed to find a second diagonal,  give up
        if not qos.containsSL  then
            return false;
        fi;

        # compute the rep of the new basis vector using the eigen basis
        x := SolutionMat( b.eigenBasis, b.basis[2] )[1];

        # compute a basis for the field over the prime field
        f := Base( pos.field / pos.primeField );

        # and generate the transvections in the small group
        L := [];
        for z  in f  do
            m := [ [ 1-z*x, z ], [ -x^2*z, z*x+1 ] ];
            r := qos.operations.Rewrite( qos, m );
            j := Value( r, [ b.eigenValue, b.eigenValue ] );
            if j = pos.field.one  then
                Add( L, [ m, r ] );
            else
                j := OrderFFE(j);
                Add( L, [ m^j, r^j ] );
            fi;
        od;

        # construct the transvection triples
        e1 := b.eigenBasis[1] * b.basis^-1;
        e1 := e1{[1,2]};
        j  := b.eigenBasis[1] * 0;
        tr := [];
        for m  in L  do
            r := Copy(j);
            r{[1,2]} := e1 * ( m[1] - m[1]^0 );
            r := r * b.basis;
            Add( tr, [ b.fixedSpace, b.eigenBasis[1], r ] );
        od;

        # warning: a simple 'Value' will destory subtrees
        st := L[1][2].operations.Values( List( L, x -> x[2] ),
                  [ b.diagonal1[2], b.diagonal2[2] ] );

        # collect independent ones
        for j  in [ 1 .. Length(st) ]  do
            diff := tr[j][3];
            if 1 < pos.field.degree  then
                diff := CRecSL.BlowupVec( pos.field, diff );
            fi;
            if CRecSL.EnlargeBasis( ind, diff )  then
                Add( tt, [ tr[j], st[j] ] );
                if Length(tt) = pos.field.degree*(pos.d-1)  then
                    pos.transvectionBasis := tt;
                    return tt;
                fi;
                InfoRecSL3( "#I  transvection is indepent, found ",
                            Length(tt), "\n" );
            fi;
        od;
    od;
    InfoRecSL2( "#I  failed to find a basis after ", i, " tries\n" );
    return false;

end;


#############################################################################
##
#F  CRecSL.SetSpanningBasis( <pos>, <tt> )  . . . . . . . find a common basis
##
CRecSL.SetSpanningBasis := function( pos, tt )
    pos.spanningBasisInverse := TransposedMat(tt);
    pos.spanningBasis := pos.spanningBasisInverse^-1;
end;


#############################################################################
##
#F  CRecSL.SpanningTransvections( <pos> ) . . . find independent fixed spaces
##
CRecSL.SpanningTransvections := function( pos )
    local   tr,  tt,  cc,  b,  t,  stat,  i,  new;

    # do we already know a set of conjugating elements
    if IsBound(pos.spanningConjugatingElements)  then
        return pos.spanningConjugatingElements;
    fi;
    CRecSL.TIME_SPANNING(pos,true);

    # compute a transvection
    tr := pos.operations.TransvectionSmall(pos);
    if tr = false  then return false;  fi;

    # use only the fixed space of the transvection
    tt := [ tr[1][1] ];
    cc := [ [ pos.identity, pos.abstractIdentity ] ];
    b  := [];
    CRecSL.EnlargeBasis( b, tt[1] );

    # and start conjugating it around
    t := 1;
    stat := 0;
    while t <= Length(tt)  do
        for i  in [ 1 .. Length(pos.generators) ]  do
            new  := tt[t] * pos.dualGenerators[i];
            stat := stat + 1;
            if CRecSL.EnlargeBasis(b,new)  then
                Add( tt, new );
                Add( cc, [ cc[t][1]*pos.generators[i],
                           cc[t][2]*pos.abstractGenerators[i] ] );
                InfoRecSL3( "#I  found ", Length(tt),
                            " spanning transvection\n" );
                if Length(tt) = pos.d  then
                    InfoRecSL2( "#I  found spanning transvections\n" );
                    pos.operations.SetSpanningBasis( pos, tt );
                    pos.spanningConjugatingElements := cc;
                    pos.statistic[CRecSL.STAT_SPANNING] := stat;
                    CRecSL.TIME_SPANNING(pos,false);
                    return pos.spanningConjugatingElements;
                fi;
            else
                InfoRecSL3( "#I  transvection is dependent\n" );
            fi;
        od;
        t := t + 1;
    od;
    pos.statistic[CRecSL.STAT_SPANNING] := -stat;
    InfoRecSL2( "#I  found only ", Length(tt[1]), " transvection\n" );

    # bind it and return
    CRecSL.TIME_SPANNING(pos,false);
    return false;

end;


#############################################################################
##

#F  CRecSL.in . . . . . . . . . . . . . . . . . . . . . . . . . 'in' operator
##
CRecSL.\in := function( obj, pos )
    local   ord;

    # check the dimension and the field (cheap test, sorry)
    if not IsMat(obj)  then
        Error( "<obj> must be a matrix" );
    elif Length(obj) <> pos.d or Length(obj[1]) <> pos.d  then
        Error( "dimension must be ", pos.d );
    elif CharFFE(obj[1][1]) <> pos.p  then
        Error( "field must be ", pos.field );
    fi;

    # <pos> must contain SL
    if not pos.containsSL  then
        Error( "this will only work if the group contains SL" );
    fi;

    # if <pos> contains GL return, otherwise test the derminant
    if pos.isGL  then
        return true;
    else
        ord := Order( pos.field, DeterminantMat(obj) );
        return pos.ext mod ord = 0;
    fi;

end;


#############################################################################
##
#F  CRecSL.AddGenerators( <pos>, <mats> ) . . . . . . . . . .  add new matrix
##
CRecSL.AddGenerators := function( pos, mats )
    local   mat,  changed;

    # if <mat> is already in (this will only work for 'containsSL'), return
    if IsMatrix(mats)  then
        mats := [ mats ];
    fi;

    # append <mat> to the generators and recompute the root
    changed := false;
    for mat  in mats  do
        if not mat in pos  then
            Add( pos.generators, mat );
            Add( pos.abstractGenerators, pos.moreAbstractGens[1] );
            pos.moreAbstractGens := 
              pos.moreAbstractGens{[2..Length(pos.moreAbstractGens)]};
            Add( pos.detOrders, Order( pos.field, DeterminantMat(mat) ) );
            pos.operations.RootDeterminant(pos);
            pos.size := pos.ext * pos.sizeSL;
            changed := true;
        fi;
    od;

    # return 'true' showing that we have succeed
    return changed;

end;


#############################################################################
##
#F  CRecSL.RootDeterminant( <pos> ) . . . . . . . .  find expression for G/SL
##
CRecSL.RootDeterminant := function( pos )
    local   lcm,  t,  i;

    # check if it is SL
    pos.isSL := ForAll( pos.detOrders, x -> x = 1 );
    pos.ext  := Lcm(pos.detOrders);
    pos.isGL := pos.ext = Size(pos.field)-1;

    # return if it is SL
    if pos.isSL  then
        return;
    fi;

    # compute an element with a generators of G/SL as determinant
    lcm := CRecSL.LcmRepresentation(pos.detOrders);

    # ok,  construct a word
    t := [ pos.group.identity, pos.abstractIdentity ];
    for i  in [ 1 .. Length(pos.generators) ]  do
        if pos.detOrders[i] <> lcm[i]  then
            t[1] := t[1] * pos.generators[i]^lcm[i];
            t[2] := t[2] * pos.abstractGenerators[i]^lcm[i];
        fi;
    od;
    t[2].fixed          := true;
    pos.root            := t;
    pos.rootDeterminant := DeterminantMat(pos.root[1]);
end;


#############################################################################
##
#F  CRecSL.Rewrite( <pos>, <g> )  . . . . . . . .  write <g> in transvections
##
CRecSL.Rewrite := function( pos, g )
    local   det,  w,  i,  l,  v,  r,  e;

    # compute the determinant of <pos>
    det := DeterminantMat(g);

    # its order must divide <pos.ext>
    if pos.ext mod Order(pos.field,det) <> 0  then
        return false;
    fi;

    # create product of <g> with <pos.root> of determinant one
    if det = pos.field.one  then
        w := pos.abstractIdentity;
    else
        i := LogFFE( det, pos.rootDeterminant );
        g := pos.root[1]^(-i) * g;
        w := pos.root[2]^i;
    fi;

    # do basis change
    g := pos.spanningBasis * g * pos.spanningBasisInverse;

    # decompose into transvections
    l := CRecSL.Diagonalize(g);

    # blow up difference and find solution
    for i  in [ 1 .. Length(l) ]  do
        v := l[i][2]
             * pos.spanningBasis
             * pos.spanningConjugatingElements[l[i][1]][1]^-1;
        if 1 < pos.field.degree  then
            v := CRecSL.BlowupVec( pos.field, v );
        fi;    
        l[i][2] := IntVecFFE( SolutionMat( pos.rewritingRules, v ) );
    od;

    # now create an expression tree
    for r  in l  do
        v := pos.abstractIdentity;
        for e  in [ 1 .. Length(r[2]) ]  do
            if r[2][e] <> 0  then
                v := v * pos.transvectionBasis[e][2]^r[2][e];
            fi;
        od;
        w := w * v^pos.spanningConjugatingElements[r[1]][2];
    od;

    # and return
    return w;

end;


#############################################################################
##
#F  CRecSL.Print( <pos> ) . . . . . . . . . . . . . . . . . . .  pretty print
##
CRecSL.Print := function( pos )
    local   max,  log,  i;

    if 0 < pos.printLevel  then
        if pos.isSL  then
            Print( "#I  <G> is SL( ", pos.d, ", ", pos.q, " )\n" );
        elif pos.isGL  then
            Print( "#I  <G> is GL( ", pos.d, ", ", pos.q, " )\n" );
        elif pos.containsSL  then
            Print( "#I  <G> is SL( ",pos.d,", ",pos.q," ) . ",pos.ext,"\n" );
        fi;
        if not pos.containsSL  then
            if IsBound(pos.transvection)  then
                Print( "#I  transvection known\n" );
            fi;
            if IsBound(pos.transvectionBasis)  then
                Print( "#I  transvection basis known\n" );
            fi;
            if IsBound(pos.spanningConjugatingElements)  then
                Print( "#I  spanning transvections known\n" );
            fi;
        fi;
    fi;
    if 1 < pos.printLevel  then
        Print( "#I  statistics:\n" );
        max := Maximum( List( CRecSL.STAT_NAMES, Length ) ) + 1;
        log := Maximum(List(pos.statistic,x->LogInt(Maximum(1,x),10)))+1;
        for i  in [ 1 .. Length(CRecSL.STAT_NAMES) ]  do
            if pos.statistic[i] <> 0   then
                Print( "#I    ", String( Concatenation( CRecSL.STAT_NAMES[i],
                       ":" ), -max ), " ", String( pos.statistic[i], log ),
                       "\n" );
            fi;
        od;
        Print( "#I  timings:\n" );
        max := Maximum( List( CRecSL.TIME_NAMES, Length ) ) + 1;
        log := Maximum(List(pos.timing,x->LogInt(Maximum(1,x),10)))+1;
        for i  in [ 1 .. Length(CRecSL.TIME_NAMES) ]  do
            if pos.timing[i] <> 0   then
                Print( "#I    ", String( Concatenation( CRecSL.TIME_NAMES[i],
                       ":" ), -max ), " ", String( pos.timing[i], log ),
                       "\n" );
            fi;
        od;
    fi;
    Print( "<< constructive SL recognition record >>" );

end;

#############################################################################
##
#F  CRecSL.SetPrintLevel( <pos>, <lev> )  . . . . . . . . . . set printl evel
##
CRecSL.SetPrintLevel := function( obj, lev )
    if 2 < lev or lev < 0  then
        Error( "invalid print level" );
    fi;
    obj.printLevel := lev;
end;


#############################################################################
##
#F  CRecSL.Setup( <G>, <gens> ) . . . . . . . . . set up a possibility record
##
CRecSL.Setup := function( G, gens )
    local   pos,  d1,  d2,  tmp,  size,  qi,  i;

    # set up pos recgonition record
    pos             := rec();
    pos.d           := Length(G.identity);
    pos.q           := Size(G.field);
    pos.p           := G.field.char;
    pos.k           := LogInt( pos.q, pos.p );
    pos.group       := G;
    pos.identity    := G.identity;
    pos.containsSL  := false;
    pos.isSL        := false;
    pos.isGL        := false;
    pos.field       := G.field;
    pos.primeField  := GF(G.field.char);
    pos.statistic   := List( CRecSL.STAT_NAMES, x -> 0 );
    pos.timing      := List( CRecSL.TIME_NAMES, x -> 0 );
    pos.printLevel  := 1;
    pos.isCRecSL    := true;
    pos.smallRandom := false;
    pos.operations  := CRecSL;

    # some "ranges" we will need a few times to remove an entry in a vector
    pos.ranges := List([1..pos.d],x->Concatenation([1..x-1],[x+1..pos.d]));

    # compute the transposed inverse of the generators
    pos.generators     := gens;
    pos.dualGenerators := List( gens, x -> TransposedMat(x^-1) );

    # compute the abstract generators (a few extract ones for the det)
    d1  := Length(pos.generators);
    d2  := Length(DivisorsInt(pos.q-1));
    tmp := RexpTree(d1+d2+100);
    pos.abstractGenerators := tmp.generators{[1..d1]};
    pos.abstractIdentity   := tmp.identity;
    pos.moreAbstractGens   := tmp.generators{[d1+1..d1+d2+100]};

    # add the size of SL
    size := 1;
    qi   := pos.q;
    for i  in [ 2 .. pos.d ]  do
        qi   := qi * pos.q;
        size := size * (qi-1);
    od;
    pos.sizeSL := pos.q^(pos.d*(pos.d-1)/2) * size;

    # and return
    return pos;

end;


#############################################################################
##

#F  CRecognizeSL( <G> ) . . . . . . . . . . . . . . . . . . recognize SL(n,q)
##
CRecognizeSL := function( arg )
    local   G,  pos,  t,  i,  j,  tmp;

    # if <G> is already a crecsl record, do not set it up again
    if not IsBound(arg[1].isCRecSL) or not arg[1].isCRecSL  then

        # get the group
        G := arg[1];

        # if the group is trivial raise an error
        if Length(G.generators) = 0  then
            Error("<G> must be non-trivial");
        fi;

        # set up pos recgonition record
        if Length(arg) = 1  then
            pos := CRecSL.Setup( G, Set(G.generators) );
        else
            pos := CRecSL.Setup( G, arg[2] );
        fi;
    else
        pos := arg[1];
        G   := pos.group;
        if 2 = Length(arg)  then
            Error( "use 'AddGenerator' to add a new generator" );
        fi;
    fi;
    CRecSL.TIME_ALL(pos,true);

    # our strategy for finding a basis for the fixed space depends on <d>
    # old: if pos.d < 3 or pos.q < 80 or pos.q < pos.d+2  then
    if pos.d < 4  then
        if pos.operations.TransvectionBasis(pos) = false  then
            CRecSL.TIME_ALL(pos,false);
            return pos;
        fi;
    else
        if pos.operations.TransvectionBasisD4(pos) = false  then
            CRecSL.TIME_ALL(pos,false);
            return pos;
        fi;
    fi;

    # try to find a set of spanning transvections
    if pos.operations.SpanningTransvections(pos) = false  then
        CRecSL.TIME_ALL(pos,false);
        return pos;
    fi;

    # convert transvections to use spanning basis and give them a name
    t := pos.transvectionBasis;
    for j in [ 1 .. Length(t) ]  do
        tmp := CRecSL.ApplyTT( pos.spanningBasis[1], t[j][1] );
        t[j][1][2]    := pos.spanningBasis[1];
        t[j][1][3]    := tmp - pos.spanningBasis[1];
        t[j][2].name  := Concatenation( "t", String(j) );
        t[j][2].fixed := true;
    od;
    tmp := List( t, x -> CRecSL.BlowupVec(pos.field,x[1][3]) );
    pos.rewritingRules := tmp;

    # now <G> must contains SL
    pos.containsSL := true;

    # check if <G> is SL
    tmp           := List( pos.generators, x -> DeterminantMat(x) );
    pos.detOrders := List( tmp, x -> Order( pos.field, x ) );

    # if <G> is not SL,  compute | <G> / SL |
    pos.operations.RootDeterminant(pos);
    pos.size := pos.ext * pos.sizeSL;

    # and return
    CRecSL.TIME_ALL(pos,false);
    return pos;

end;

CRecogniseSL := CRecognizeSL;
