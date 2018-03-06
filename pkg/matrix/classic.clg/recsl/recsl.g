#############################################################################
##
#A  recsl.g                     GAP library                      Frank Celler
##
#H  @(#)$Id: recsl.g,v 1.1 1997/03/10 13:49:13 gap Exp $
##
#Y  Copyright (C) 1995,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This  file  contains  functions  which  will  help  to  recognize  groups
##  containing SL(n,q).
##
Revision_recsl_g :=
    "@(#)$Id: recsl.g,v 1.1 1997/03/10 13:49:13 gap Exp $";


#############################################################################
##

#V  InfoRecSL?(...) . . . . . . . . . . . . . . . . . . . . . . .  infomation
##
##  InfoRecSL3:  information about number of possible chevalley groups
##  InfoRecSL4:  runtime information
##  InfoRecSL5:  order algorithm
##
if not IsBound(InfoRecSL1)   then InfoRecSL1  := Ignore;  fi;
if not IsBound(InfoRecSL2)   then InfoRecSL2  := Ignore;  fi;
if not IsBound(InfoRecSL3)   then InfoRecSL3  := Ignore;  fi;
if not IsBound(InfoRecSL4)   then InfoRecSL4  := Ignore;  fi;
if not IsBound(InfoRecSL5)   then InfoRecSL5  := Ignore;  fi;


#############################################################################
##

#V  RecSL . . . . . . . . . . . . . . . . . .  functions to recognise SL(n,q)
##
RecSL := rec();
RecSL.STAT_ALT       :=  1;
RecSL.STAT_CHEV      :=  2;
RecSL.STAT_CLASSICAL :=  3;
RecSL.STAT_LARGER    :=  4;
RecSL.STAT_PGROUP    :=  5;
RecSL.STAT_POWER     :=  6;
RecSL.STAT_PRIMITIVE :=  7;
RecSL.STAT_PRODUCT   :=  8;
RecSL.STAT_REDUCIBLE :=  9;
RecSL.STAT_SMALLER   := 10;
RecSL.STAT_SPORADIC  := 11;


#############################################################################
##

#V  RecSL.GL.<func>( <d>, <p>, <k> )  . . . . . . . . .  general linear group
##
##  <RecSL.GL>  contains  a  collection of  functions to  extract information
##  about GL(<d>,<p>^<k>).
##
RecSL.GL := rec();

RecSL.GL.uexponent := function( d, p, k )
    local   q,  exp,  qi,  i;

    q := p^k;
    if d = 1  then
        return q-1;
    else
        exp := 1;
        qi  := 1;
        for i  in [ 1 .. d ]  do
            qi  := qi * q;
            exp := LcmInt( exp, qi-1 );
        od;
        return p^(1+LogInt(d-1,p)) * exp;
    fi;

end;


#############################################################################
##
#V  RecSL.PGL.<func>( <d>, <p>, <k> ) . . . . projective general linear group
##
##  <RecSL.PGL>  contains  a collection of  functions to  extract information
##  about PGL(<d>,<p>^<k>).
##
RecSL.PGL := rec();

RecSL.PGL.uexponent := function( d, p, k )
    local   q,  exp,  qi,  i;

    q := p^k;
    if d = 1  then
        return 1;
    else
        exp := 1;
        qi  := 1;
        for i  in [ 1 .. d ]  do
            qi  := qi * q;
            exp := LcmInt( exp, qi-1 );
        od;
       return p^(1+LogInt(d-1,p)) * exp;
    fi;

end;


#############################################################################
##
#V  RecSL.GU.<func>( <d>, <p>, <k> )  . . . . . . . . . general unitary group
##
##  <RecSL.GL>  contains  a  collection of  functions to  extract information
##  about GU(<d>,<p>^<k>) < GL(<d>,<p>^2<k>).
##
RecSL.GU := rec();

RecSL.GU.order := function( d, p, k )
    local   i,  o,  q;

    q := p^k;
    o := 1;
    for i  in [ 1 .. d ]  do
        o := o * ( q^i - (-1)^i );
    od;
    return o*q^(d*(d-1)/2);

end;

RecSL.GU.uexponent := function( d, p, k )
    local   i,  e,  q,  qi;

    q := p^k;
    if d = 1  then
        return q+1;
    fi;
    e  := p^(1+LogInt(d-1,p));
    qi := 1;
    for i  in [ 1 .. d ]  do
        qi := qi * q;
        e  := LcmInt( e, qi - (-1)^i );
    od;
    return GcdInt( e, RecSL.GU.order( d, p, k ) );
end;


#############################################################################
##
#V  RecSL.O.<func>( <s>, <d>, <p>, <k> )  . . . . . . . . .  orthogonal group
##
##  <RecSL.O>  contains  a  collection of  functions to  extract information
##  about O<s>(<d>,<p>^<k>).
##
RecSL.O := rec();

RecSL.O.order := function( s, d, p, k )
    local   q,  o,  q2,  qi,  i;

    if s = 0 and d mod 2 = 0  then
        Error( "sign zero but dimension even" );
    fi;

    q  := p^k;
    o  := 1;
    q2 := q^2;
    qi := 1;
    if s = 0  then
        for i  in [ 1 .. (d-1)/2 ]  do
            qi := qi * q2;
            o  := o * (qi-1);
        od;
        o := 2*q^((d-1)^2/4)*o;
    else
        for i  in [ 1 .. d/2-1 ]  do
            qi := qi * q2;
            o  := o * (qi-1);
        od;
        o := 2*q^(d*(d-2)/4)*(q^(d/2)-s)*o;
    fi;
    return o;

end;

RecSL.O.uexponent := function( s, d, p, k )
    local   q,  e,  qi,  i;

    if d = 1  then return 2;  fi;
    q  := p^k;
    e  := p^(1+LogInt(d-1,p));
    qi := 1;
    if s = 0  then
        for i  in [ 1 .. (d-1)/2 ]  do
            qi := qi * q;
            e := LcmInt( e, LcmInt( qi-1, qi+1 ) );
        od;
    else
        for i  in [ 1 .. d/2 ]  do
            qi := qi * q;
            e := LcmInt( e, LcmInt( qi-1, qi+1 ) );
        od;
    fi;
    return GcdInt( e, RecSL.O.order( s, d, p, k ) );
end;


#############################################################################
##
#V  RecSL.SP.<func>( <d>, <p>, <k> )  . . . . . . . . . . .  symplectic group
##
##  <RecSL.SP>  contains  a  collection of  functions to  extract information
##  about SP(<d>,<p>^<k>).
##
RecSL.SP := rec();

RecSL.SP.order := function( d, p, k )
    local   i,  o,  q,  qi,  q2;

    d  := d/2;
    q  := p ^ k;
    o  := 1;
    q2 := q^2;
    qi := 1;
    for i  in [ 1 .. d ]  do
        qi := qi * q2;
        o  := o * ( qi - 1 );
    od;
    o  := o * q^(d^2);
    return o;

end;

RecSL.SP.uexponent := function( d, p, k )
    local   e,  i,  q,  qi;

    e  := p ^ ( 1 + LogInt( d-1, p ) );
    q  := p ^ k;
    qi := 1;
    for i  in [ 1 .. d/2 ]  do
        qi := qi * q;
        e  := Lcm( e, qi - 1 );
        e  := Lcm( e, qi + 1 );
    od;
    return Gcd( e, RecSL.SP.order( d, p, k ) );

end;


#############################################################################
##

#F  RecSL.UpperBoundOrder( <pos>, <f> ) . . . . . . . . . . . . . upper bound
##
##  Computes the  irreducible factors f_i  of a polynomial  <f>  over a field
##  with  p^n  elements.  It returns a list  l of triples (f_i,a_i,pp_i) such
##  that the p-part  of x  mod  f_i is p^a_i and  the p'-part divides d_i for
##  which the semi prime powers pp_i are given.
##
RecSL.UpperBoundOrder := function( pos, f )
    local   F,  fs,  a,  pp,  f,  d,  phi,  B,  i,  L;

    # factorize <f> into irreducible factors
    InfoRecSL5( "#I  UpperBoundOrder:\n" );
    fs := Collected(Factors(f));

    # get the field over which the polynomials are written
    F := pos.field;

    # store a list of semi prime powers in <pos>
    if not IsBound(pos.phiList)  then
        pos.phiList := [ PrimePowersInt(F.char-1) ];
    fi;
    if not IsBound(pos.semiPrimes)  then
        pos.semiPrimes := [];
    fi;

    # <phi>(m) gives ( minpol of 1^(1/m) )( F.char )
    L   := pos.phiList;
    phi := function( m )
        local   x,  d,  pp,  i,  tmp;
        if not IsBound( L[m] )  then
            x := F.char^m-1;
            for d  in Difference( DivisorsInt( m ), [m] )  do
                pp := phi(d);
                for i  in [ 1 .. Length(pp)/2 ]  do
                    x := x / pp[2*i-1]^pp[2*i];
                od;
            od;
            tmp := SemiPrimePowersInt( x, pos.p, pos.d*pos.k );
            L[m] := Concatenation(tmp);
            AddSet( pos.semiPrimes, tmp[2]{([1..Length(tmp[2])/2]*2)-1} );
        fi;
        return L[m];
    end;

    # compute a_i and pp_i
    B := [];
    for i  in [ 1 .. Length(fs) ]  do

        # p-part is p^Roof(Log_p(e_i)) where e_i is the multiplicity of f_i
        a := 0;
        if fs[i][2] > 1  then
            a := 1+LogInt(fs[i][2]-1,F.char);
        fi;

        # p'-part: (p^n)^d_i-1/(p^n-1) where d_i is the degree of f_i
        d := Degree(fs[i][1]);
        InfoRecSL5( "#I    irreducible factor of degree ", d, "\n" );
        pp := [];
        for f  in DivisorsInt( d*F.degree )  do
            if F.degree mod f <> 0  then
                pp := ProductPP( phi(f), pp );
            fi;
        od;

        # add <a> and <pp> to <B>
        Add( B, [fs[i][1],a,pp] );
    od;

    # OK, that's it
    InfoRecSL5( "#I  UpperBoundOrder returns ", B, "\n" );
    return B;

end;


#############################################################################
##
#F  RecSL.OrderScalar( <pos>, <f> ) semi projective order of a polynomial <f>
##
##  Return an integer n and  a finite  field element  const  such that x^n  =
##  const mod <f>.  <n> is as small as possible as  long as the factorisation
##  works, otherwise it the  smallest product of semi  primes.  Note that q-1
##  must not be a semi prime.
##
RecSL.OrderScalar := function( pos, f )
    local   U,  R,  x,  O,  n,  g,  q,  o;

    # if degree is zero, return
    if 0 = Degree(f)  then
        return [ 1, f.coefficients[1] ];
    fi;

    # use 'UpperBoundOrder' to split <f> into irreducibles
    U := RecSL.UpperBoundOrder( pos, f );

    # run through the irrducibles and compute their order
    R := PolynomialRing(pos.field);
    x := Polynomial( R.baseRing, [ R.baseRing.zero, R.baseRing.one ] );
    O := [];
    n := 1;
    for g  in U  do
        o := R.operations.OrderKnownDividend( R, x mod g[1], g[1], g[3] );
        q := R.baseRing.char^g[2];
        n := Lcm( n, o[1]*q );
        Add( O, [ o[1]*q, o[2]^q ] );
    od;

    # try to get the same constant in each block
    U := [];
    q := Size( R.baseRing ) - 1;
    for g  in O  do
        AddSet( U, g[2]^((n/g[1]) mod q) );
    od;

    # return the order <n> times the order of <U>
    o := OrderKnownDividendList( U, PrimePowersInt(q) );
    return [ n*o[1], o[2] ];

end;


#############################################################################
##
#F  RecSL.Random( <G> ) . . . . . . . . . . .  return a random element of <G>
##
RecSL.Random := function( G )
    local  i,  j;
    
    # catch trivial case
    if Length(G.generators) = 0  then
        return G.identity;
    fi;

    # if '<G>.randomSeed' is unbound, construct a start
    if not IsBound(G.randomSeed)  then
        G.randomSeed := Set(G.generators);
        for i  in [ Length(G.randomSeed)+1 .. 10 ]  do
            G.randomSeed[i] := Random(G.generators);
        od;
        for i  in [ 1 .. 10*Length(G.generators) ]  do
            RecSL.Random(G);
        od;
    fi;
    
    # get two random position and a random L/R
    i := Random( [ 1 .. Length(G.randomSeed) ] );
    repeat
        j := Random( [ 1 .. Length(G.randomSeed) ] );
    until i <> j;
    if Random( [ true, false ] )  then
        G.randomSeed[j] := G.randomSeed[j]*G.randomSeed[i];
    else
        G.randomSeed[j] := G.randomSeed[i]*G.randomSeed[j];
    fi;
    return G.randomSeed[j];
    
end;


#############################################################################
##

#F  RecSL.SetAlternating( <pos> ) . . . . . . . . . . . set alternating group
##
RecSL.SetAlternating := function( pos )
    local   timer,  n,  o,  x;

    # start timer
    timer := Runtime();

    # store possible alternating groups in <pos.alternating>
    pos.alternating := [];

    # get alternating groups of degree <= d+2 (4<n)
    n := 5;
    o := Factorial(5)/2;
    if IsBound(pos.orderGL) and pos.orderGL mod o <> 0  then 
        pos.isAlternating := false;
        InfoRecSL3( "#I  0 alternating groups, A5 !< G\n" );
        return;
    fi;
    if IsBound(pos.orderGL)  then x := pos.orderGL/12;  else x := 0;  fi;
    while ( n <= pos.d+2 ) and ( x mod n = 0 ) and ( n < 9 )  do
        x := x / n;
        if n = 5 and ( 2 < pos.d or ( pos.q <> 4 and pos.q <> 5 ) )  then
            if 2 = pos.d and 9 = pos.q  then
                InfoRecSL3( "#I  A_6 = PSL(2,9) > A_5 = PSL(2,4), ",
                            "check occurs in Chevalley\n" );
            else
                Add( pos.alternating, n );
            fi;
        elif n = 6 and ( 2 < pos.d or pos.q <> 9 )  then
            Add( pos.alternating, n );
        elif n = 7 and 3 <= pos.d  then
            if pos.d = 3 and pos.q = 4  then
                InfoRecSL3( "#I  avoiding |PSL(3,4)| !> |A_7|\n" );
            else
                Add( pos.alternating, n );
            fi;
        elif n = 8 and 4 <= pos.d and ( 4 < pos.d or pos.q <> 2 )  then
            if pos.d = 3 and pos.q = 4  then
                InfoRecSL3( "#I  avoiding |PSL(3,4)| = |A_8|\n" );
            else
                Add( pos.alternating, n );
            fi;
        fi;
        n := n + 1;
        o := o * n;
    od;
    while ( n <= pos.d+2 ) and ( x mod n = 0 )  do
        Add( pos.alternating, n );
        x := x / n;
        n := n + 1;
        o := o * n;
    od;

    # fix the exceptions
    if pos.d = 2  then
        if pos.p = 2  then
            if pos.q <> 4  then
                Add( pos.alternating, 5 );
            else
                InfoRecSL3( "#I  avoiding SL(2,4) = A_5\n" );
            fi;
        fi;
    fi;
    if pos.d = 3  then
        if pos.p = 3  then
            Add( pos.alternating, 6 );
        fi;
    fi;
    if pos.d = 4  then
        if pos.p = 2  then
            Add( pos.alternating, 7 );
            if pos.q <> 2  then
                Add( pos.alternating, 8 );
            else
                InfoRecSL3( "#I  avoiding SL(4,2) = A_8\n" );
            fi;
        fi;
    fi;

    pos.alternating := Set(pos.alternating);
    InfoRecSL3( "#I  ", Length(pos.alternating), " alternating groups\n" );
    pos.isAlternating := 0 < Length(pos.alternating);

    # show runtime
    InfoRecSL4( "#I  alternating: ", Runtime()-timer, " msec\n" );

end;


#############################################################################
##
#F  RecSL.CheckAlternating( <pos>, <pord>, <psod> ) . . .  element order test
##
RecSL.CheckAlternating := function( pos, pord, psod )
    local   z,  new,  n,  z2;

    # are alternating groups still possible?
    if not pos.isAlternating then return;  fi;
    pos.statistic[RecSL.STAT_ALT] := pos.tries;

    # avoid big primes
    new := [];
    for n  in pos.alternating  do
        z := 2*Factorial(n);
        if z mod pord = 0  then
            if 1 = psod or 1 < GcdInt(psod,n)  then
                Add( new, n );
            fi;
        fi;
    od;
    pos.alternating := new;
    if 0 = Length(pos.alternating)  then
        pos.isAlternating := false;
        InfoRecSL2( "#I  <G> is not an alternating group,  ",
                    "group order criteria failed\n" );
        return;
    fi;

    # compute minimal cycle length
    z := Sum( List( Collected(Factors(pord)), x -> x[1]^x[2] ) );

    # check element orders
    new := [];
    for n  in pos.alternating  do
        if n = 6  then
            z2 := Sum( Collected(Factors(pord/Gcd(pord,2))), x->x[1]^x[2] );
            if z2 <= n  then
                Add( new, n );
            fi;
        elif z <= n  then
            Add( new, n );
        fi;
    od;
    pos.alternating := new;
            
    if 0 = Length(pos.alternating)  then
        pos.isAlternating := false;
        InfoRecSL2( "#I  <G> is not an alternating group,  ",
                    "element order criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSL.SetChevalley( <pos> ) . . . . . . . . . . possible Chevalley groups
##
RecSL.SetChevalley := function( pos )
    local   chev,  par,  t,  timer,  name;

    # start timer
    timer := Runtime();

    # we use determinants and orders for GL(2,9)
    if pos.d = 2 and pos.q = 9  then
        pos.isChevalley := true;
        pos.expsChev := [ [ 0, "A", 1, 4, 1 ], [ 0, "A", 1, 5, 1 ] ];
        pos.squares  := Set( List( Elements(pos.field), x -> x^2 ) );
        InfoRecSL4( "#I  chevalley: ", Runtime()-timer, " msec\n" );
        return;
    fi;

    # run through the simple Chevalley groups
    pos.isChevalley := true;
    pos.expsChev    := [];
    for chev in [ ChevA, ChevB, ChevC, ChevD, ChevG,
                  ChevF, ChevE, Chev2A, Chev2B, Chev2G,
                  Chev2F, Chev2D, Chev3D, Chev2E ]
    do

        # extract name and start timer
        name  := chev.name{[5..Length(chev.name)]};
        timer := Runtime();

        # check the order in small cases
        if IsBound(pos.orderGL)  then
            par := chev.possibilities( pos.d, pos.q3d, pos.orderGL, pos.p );
        else
            par := chev.possibilities( pos.d, pos.q3d, 0, pos.p );
        fi;
        InfoRecSL4( "#I  chevalley: possibilities ", chev.name, ": ",
                    Runtime()-timer, "\n" );
        timer := Runtime();
        if 0 < Length(par[1])+Length(par[2])  then
            InfoRecSL3("#I  ", chev.name, ": ", Length(par[1]), " over ", 
                       pos.p, " and ", Length(par[2]), " over coprime ",
                       "characteristic\n");
        fi;
        for t  in Concatenation(par[1], par[2])  do
            if chev = ChevA  then
               if t[1]+1 = pos.d and t[2] = pos.p and t[3] = pos.k  then
                   InfoRecSL3( "#I  ChevA: avoiding SL(d,q)\n" );
               elif pos.d=2 and pos.q=4 and t[1]=1 and t[2]=5 and t[3]=1 then
                   InfoRecSL3( "#I  ChevA: avoiding L2(5)=L2(4)\n" );
               elif pos.d=2 and pos.q=5 and t[1]=1 and t[2]=2 and t[3]=2 then
                   InfoRecSL3( "#I  ChevA: avoiding L2(5)=L2(4)\n" );
               elif pos.d=3 and pos.q=2 and t[1]=1 and t[2]=7 and t[3]=1 then
                   InfoRecSL3( "#I  ChevA: avoiding L3(2)=L2(7)\n" );
               elif pos.d=3 and pos.q=4 and t[1]=3 and t[2]=2 and t[3]=1 then
                   InfoRecSL3( "#I  ChevA: avoiding L3(4)<>L4(2)\n" );
               else
                   Add( pos.expsChev,
                        [ chev.uexponent(t[1],t[2],t[3]) 
                          * chev.automorphism(t[1],t[2],t[3]),
                          name, t[1], t[2], t[3], chev ] );
               fi;
            else
                Add( pos.expsChev,
                     [ chev.uexponent(t[1],t[2],t[3]) 
                       * chev.automorphism(t[1],t[2],t[3]),
                       name, t[1], t[2], t[3], chev ] );
            fi;
        od;
        InfoRecSL4( "#I  chevalley: exponents ", chev.name, ": ",
                    Runtime()-timer, "\n" );
    od;

    # show runtime
    InfoRecSL4( "#I  chevalley: ", Runtime()-timer, " msec\n" );

end;


#############################################################################
##
#F  RecSL.CheckChevalley( <pos>, <pord>, <psod>, <elm> ) . . .  element order
##
RecSL.CheckChevalley := function( pos, pord, psod, elm )
    local   d;

    # are Chevalley groups still possible?
    if not pos.isChevalley  then return;  fi;
    pos.statistic[RecSL.STAT_CHEV] := pos.tries;

    # special case: PSL(2,9) = A_6 > A_5 = PSL(2,4)
    if pos.d = 2 and pos.q = 9  then
        if pord = 8 or pord = 4  then
            d := DeterminantMat(elm);
            if d in pos.squares  then
                pos.isChevalley := false;
                InfoRecSL2( "#I  <G> is not a Chevalley group,  ",
                            "non-det-square element of porder 4/8\n" );
            fi;
        fi;

    # check element orders
    else
        pos.expsChev := Filtered(pos.expsChev, x -> x[1] mod pord = 0);
        if 1 < psod  then
            pos.expsChev := Filtered(pos.expsChev, x -> 1<GcdInt(x,psod));
        fi;
        if 0 = Length(pos.expsChev)  then
            pos.isChevalley := false;
            InfoRecSL2( "#I  <G> is not a Chevalley group,  ",
                        "element order criteria failed\n" );
        fi;
    fi;

end;


#############################################################################
##
#F  RecSL.SetClassicalGroups( <pos> ) . . . . . . . possible classical groups
##
RecSL.SetClassicalGroups := function( pos )
    local   gens,  id,  i,  j;

    # orthogonal group: if <q> is even, SP > GO
    pos.isOrthogonal := (2 < pos.p) or (1 = pos.d mod 2) or (2 = pos.d);
    if pos.q = 2 and pos.d = 2  then
        pos.expsOrthogonal := [ 2, 3 ];
    fi;
    if not pos.isOrthogonal  then
        InfoRecSL2( "#I  symplectic > orthogonal, ignoring orthogonal\n" );
    fi;

    # symplectic group, dimension must be even
    pos.isSymplectic := 0 = pos.d mod 2;
    if not pos.isSymplectic  then
        InfoRecSL2( "#I  <G> is not symplectic, dimension is odd\n" );
    fi;

    # unitary group, <q> must be a square
    pos.isUnitary := 0 = pos.k mod 2;
    if not pos.isUnitary  then
        InfoRecSL2( "#I  <G> is not unitary, field has no square root\n" );
    fi;

    # in dimension 2, SL_2 = SP_2 = SU_2
    if pos.d = 2  then
        InfoRecSL2( "#I  dimension is 2, SL = SP = SU\n" );
        pos.isSymplectic := false;
        pos.isUnitary    := false;
    fi;

    # or all cases together
    pos.isClassical := pos.isSymplectic or pos.isOrthogonal or pos.isUnitary;

end;


#############################################################################
##
#F  RecSL.CheckClassicalGroups( <pos>, <cpol>, <pord>, <elm> )  check classic
##
RecSL.CheckClassicalGroups := function( pos, cpol, pord, elm )
    local   c,  I,  g,  d,  e,  t,  i0,  a,  l;

    # classical groups still possible?
    if not pos.isClassical  then return;  fi;
    pos.statistic[RecSL.STAT_CLASSICAL] := pos.tries;

    # is <G> orthogonal in small dimenstion
    if pos.isOrthogonal and 2 = pos.d  then
        if 2 = pos.q  then
            pos.expsOrthogonal := Filtered( pos.expsOrthogonal,
                                            x -> x mod pord = 0 );
            pos.isOrthogonal   := 0 < Length(pos.expsOrthogonal);
        else
            pos.isOrthogonal := pos.isOrthogonal 
                                and ( ( 0 = (2*(pos.q-1)) mod pord )
                                   or ( 0 = (2*(pos.q+1)) mod pord ) );
        fi;
        if not pos.isOrthogonal  then
            InfoRecSL2("#I  <G> is not orthogonal,  element order failed\n");
        else
            c := elm^2;
            I := pos.identity;
            if c <> I  then
                for g  in pos.generators  do
                    if pos.isOrthogonal  then
                        d := c^g;
                        e := d / c;
                        if e <> e[1][1]*I  then
                            e := d * c;
                            if e <> e[1][1]*I  then
                                pos.isOrthogonal := false;
                            fi;
                        fi;
                    fi;
                od;
                if not pos.isOrthogonal  then
                    InfoRecSL2( "#I  <G> is not orthogonal,  ",
                                "dihedral test failed\n" );
                fi;
            fi;
        fi;
    fi;

    # get coefficient list of the characteristic polynomial
    c := cpol.coefficients;
    d := pos.d;

    # is <G> orthogonal or symplectic
    if pos.isOrthogonal or pos.isSymplectic  then
        I := Filtered( [0..d], x -> c[x+1] <> cpol.baseRing.zero );
        if ForAny( I, x -> c[d-x+1] = cpol.baseRing.zero )  then
            pos.isSymplectic := false;
            pos.isOrthogonal := false;
        else
            t  := GcdRepresentation(I);
            i0 := I*t;
            a  := c[1];
            l  := List( [1..Length(I)], x ->(a*c[d-I[x]+1]/c[I[x]+1]) );
            g  := Product( [1..Length(I)],x -> l[x]^t[x] );
            if ForAny( [1..Length(I)], x -> l[x] <> g^(I[x]/i0) )  then
                pos.isSymplectic := false;
                pos.isOrthogonal := false;
            fi;
        fi;
        if not ( pos.isOrthogonal or pos.isSymplectic )  then
            InfoRecSL2( "#I  <G> is not orthogonal/symplectic,  ",
                        "characteristic polynomial failed\n" );
        fi;
    fi;

    # is <G> unitary
    if pos.isUnitary  then
        I := Filtered( [0..d], x -> c[x+1] <> cpol.baseRing.zero );
        if ForAny( I, x -> c[d-x+1] = cpol.baseRing.zero )  then
            pos.isUnitary := false;
        else
            t  := GcdRepresentation(I);
            i0 := I*t;
            a  := c[1];
            l  := List([1..Length(I)], x ->(a*c[d-I[x]+1]^pos.qq/c[I[x]+1]));
            g  := Product( [1..Length(I)],x -> l[x]^t[x] );
            if ForAny( [1..Length(I)], x -> l[x] <> g^(I[x]/i0) )  then
                pos.isUnitary := false;
            fi;
        fi;
        if not pos.isUnitary  then
            InfoRecSL2( "#I  <G> is not unitary,  characteristic ",
                        "polynomial failed\n" );
        fi;
    fi;

    # or all cases together
    pos.isClassical := pos.isSymplectic or pos.isOrthogonal or pos.isUnitary;

end;


#############################################################################
##
#F  RecSL.SetLargerField( <pos> ) . . . . . . . . . . . possible larger field
##
##  d = e*f, G < TL(e,q^f) = GL(e,q^f).f
##
RecSL.SetLargerField := function( pos )
    local   d;

    # no reduction is possible if dimension is one
    pos.isLarger := 1 < pos.d;
    if not pos.isLarger  then
        InfoRecSL2("#I  <G> is not definable over a larger field,  ",
                   "dimension is one\n" );

    # representation that preserve a larger field
    elif 2 < pos.d  then
        d := Filtered( DivisorsInt(pos.d), x -> 1 < x );
        d := List( d, x -> RecSL.GL.uexponent(pos.d/x,pos.p,x*pos.k)*x );
        pos.expsLarger := d;

    # catch GL(2,2)
    elif 2 = pos.q  then
        pos.isLarger := false;
    fi;

end;


#############################################################################
##
#F  RecSL.CheckLargerField( <pos>, <ord>, <elm> ) . . . . check element order
##
RecSL.CheckLargerField := function( pos, ord, elm )
    local   c,  I,  g,  d,  e;

    # are larger fields still possible?
    if not pos.isLarger  then return;  fi;
    pos.statistic[RecSL.STAT_LARGER] := pos.tries;

    # check element orders (if d > 2)
    if 2 < pos.d  then
        pos.expsLarger := Filtered(pos.expsLarger, x -> x mod ord = 0);
        if 0 = Length(pos.expsLarger)  then
            pos.isLarger := false;
            InfoRecSL2( "#I  <G> is not definable over a larger field,  ",
                        "element order criteria failed\n" );
        fi;

    # dihedral test
    else
        c := elm^2;
        I := pos.identity;
        if c <> I  then
            for g  in pos.generators  do
                if pos.isLarger  then
                    d := c^g;
                    e := d / c;
                    if e <> e[1][1]*I  then
                        e := d / c^pos.q;
                        if e <> e[1][1]*I  then
                            pos.isLarger := false;
                        fi;
                    fi;
                fi;
            od;
            if not pos.isLarger  then
                InfoRecSL2( "#I  <G> is not defined over TL(q^2),  ",
                            "dihedral test failed\n" );
            fi;
        fi;
    fi;

end;


#############################################################################
##
#F  RecSL.SetMysteriousPGroup( <pos> )  . . . . . .  is mysterious p possible
##
##  R < G < N,  where N is the normalizer of an extraspecial r-group R.
##
RecSL.SetMysteriousPGroup := function( pos )
    local   r,  m,  e;

    # mysterious p group: dimension must be a prime power
    pos.isMysteriousP := false;
    if pos.k mod 2 = 1 and IsPrimePowerInt(pos.d)  then
        r := SmallestRootInt( pos.d );
        m := LogInt( pos.d, r );
        if r <> pos.p  then
            e := 1;
            while pos.p^e mod r <> 1  do e := e + 1;  od;
            if e = pos.k  then
                if e*r mod 2 = 1  then
                    pos.isMysteriousP  := true;
                    pos.expMysteriousP := r^2*RecSL.SP.uexponent(2*m,r,1);
                elif e = 1 and r = 2 and 2 = pos.d and 3 <> pos.q  then
                    pos.isMysteriousP  := true;
                    pos.expMysteriousP := 4*6;
                fi;
            fi;
            if r = 2  and 4 <= pos.d and pos.k = 1  then
                e := 1;
                while pos.p^e mod 4 <> 1  do e := e + 1;  od;
                if e = 1  then
                    pos.isMysteriousP  := true;
                    pos.expMysteriousP := 4*RecSL.SP.uexponent(2*m,2,1);
                fi;
            fi;
        fi;
    fi;
    if not pos.isMysteriousP  then
        InfoRecSL2( "#I  <G> is no mysterious p-group,  dimension ",
                    "is not a prime power\n" );
    fi;
    
end;


#############################################################################
##
#F  RecSL.CheckMysteriousPGroup( <pos>, <pord> )  . . . . check element order
##
##  <pord> is the correct part of the projective order.
##
RecSL.CheckMysteriousPGroup := function( pos, pord )

    # mysterious p-group still possible?
    if not pos.isMysteriousP  then return;  fi;
    pos.statistic[RecSL.STAT_PGROUP] := pos.tries;

    # check element order
    pos.isMysteriousP := pos.expMysteriousP mod pord = 0;
    if not pos.isMysteriousP  then
        InfoRecSL2( "#I  <G> is no mysterious p-group,  element order ",
                    "criteria failed\n" );
    fi;
    
end;


#############################################################################
##
#F  RecSL.SetImprimitive( <pos> ) . . . . . . . . possible imprimitive groups
##
##  d = e*f,  G < GL(e,q) wr Sym(f)
##
RecSL.SetImprimitive := function( pos )
    local   d;

    # store the various pairs <e>, <f>
    d := Filtered( DivisorsInt(pos.d), x -> x <> pos.d );
    pos.dimsImprimitive := List( d, x -> [ x, pos.d/x ] );
    pos.isImprimitive := true;
    
end;


#############################################################################
##
#F  RecSL.CheckImprimitive( <pos>, <ord> )  . . . . . . . check element order
##
##  Assume that all integer <= <pos.d> are no semi primes.
##
RecSL.CheckImprimitive := function( pos, ord )
    local   new,  p,  m;

    # imprimitivity impossible?
    if not pos.isImprimitive  then return;  fi;
    pos.statistic[RecSL.STAT_PRIMITIVE] := pos.tries;

    # check exponents
    new := [];
    for p  in pos.dimsImprimitive  do
        m := ord / Gcd( ord, pos.expsGL[p[1]] );
        if Factorial(p[2]) mod m = 0 
           and Sum(Collected(Factors(m)),x->x[1]^x[2]) <= p[2]
        then
            Add( new, p );
        fi;
    od;
    pos.dimsImprimitive := new;
    
    # if <new> is trivial no reduction is possible
    if 0 = Length(new)  then
        pos.isImprimitive := false;
        InfoRecSL2( "#I  <G> is not imprimitive,  element order ",
                    "criteria failed\n" );
    fi;
    
end;


#############################################################################
##
#F  RecSL.SetReducible( <pos> ) . . . . . . . . . . possible reducible groups
##
##  <pos>.dimsReducible is eitehr empty (all dimensions are possible) or is a
##  list of all possible dimensions.  At first it is empty.
##
RecSL.SetReducible := function( pos )
    pos.dimsReducible := [];
    pos.isReducible := true;
end;


#############################################################################
##
#F  RecSL.CheckReducible( <pos>, <cpol> ) . . . . . . . . . . . .  check cpol
##
##  Compute the  degrees  of  the irreducible factors of  <cpol>  and  update
##  <pos>.dimsReducible.
##
RecSL.CheckReducible := function( pos, cpol )
    local   deg,  dims,  g;

    # reducible groups still possible?
    if not pos.isReducible  then return;  fi;
    pos.statistic[RecSL.STAT_REDUCIBLE] := pos.tries;

    # compute the degrees of the irreducible factors
    deg := List(Factors(cpol), EuclideanDegree);
    
    # compute all possible dims (one could restrict this to 2s <=d)
    dims := [0];
    for g  in deg  do
        UniteSet( dims, dims+g );
    od;
    
    # and intersect it with <pos>.dimsReducible
    if 0 = Length(pos.dimsReducible)  then
        pos.dimsReducible := dims;
    else
        IntersectSet( pos.dimsReducible, dims );
    fi;
    
    # G acts irreducible if only 0 and d are possible
    if 2 = Length(pos.dimsReducible)  then
        pos.isReducible := false;
        InfoRecSL2("#I  <G> acts irreducible,  block criteria failed\n");
    fi;
    
end;


#############################################################################
##
#F  RecSL.SetSmallerField( <pos> )  . . . . . . . . .  possible smaller field
##
##  If <q> = <p>^<k> is a prime, then the matrices are already written over a
##  prime  field,  no  reduction  is possible in this case.   Otherwise upper
##  bounds for the exponents  of PGL(<n>,<p>^<i>),  <i> T <k>,  are computed.
##  The record component <pos>.smallerField containes the smaller field still
##  possible.
##
RecSL.SetSmallerField := function( pos )
    local   i,  j,  q,  exp;

    # if <field> is the prime field,  no reduction is possible
    pos.isSmaller := 1 < pos.k;
    if not pos.isSmaller  then
        InfoRecSL2( "#I  <G> is not definable over a smaller field,  ",
                    "field is the prime field\n" );

    # loop over the maximal divisors of <k>
    else
        pos.expsSmaller := [];
        for i in DivisorsInt(pos.k) do
            if i < pos.k and IsPrime(pos.k/i)  then
                AddSet( pos.expsSmaller,RecSL.GL.uexponent(pos.d,pos.p,i));
            fi;
        od;
        pos.smallerField := GF(pos.p);
    fi;         

end;        


#############################################################################
##
#F  RecSL.CheckSmallerField( <pos>, <cpol>, <pord> )  check order and charpol
##
##  <pord> is the projective order  of the  group  element  A, <cpol>  is the
##  characteristic polynomial  of A.  The function  first  checks  if  <pord>
##  divides  at least one  of the exponents stored in  <pos>.expsSmaller.  It
##  then computes the smallest field which contains <cpol> * zeta.
##
RecSL.CheckSmallerField := function( pos, cpol, pord )
    local   c,  d,  t,  p,  i0,  I;

    # are smaller fields still possible?
    if not pos.isSmaller  then return;  fi;
    pos.statistic[RecSL.STAT_SMALLER] := pos.tries;

    # first check the projective order
    pos.expsSmaller := Filtered( pos.expsSmaller, x -> x mod pord = 0 );
    if 0 = Length(pos.expsSmaller)  then
        pos.isSmaller := false;
        InfoRecSL2( "#I  <G> is not definable over a smaller field,  ",
                    "element order criteria failed\n" );

    # otherwise compute the smallest field containing <cpol> * zeta
    else
        c  := cpol.coefficients;
        d  := pos.d;
        I  := Filtered( [0..d], x -> c[d-x+1] <> cpol.baseRing.zero );
        t  := GcdRepresentation(I);
        i0 := I*t;
        p  := Product( [1..Length(I)], x -> c[d-I[x]+1]^(-t[x]) );
        c  := List( [1..Length(I)], x -> c[d-I[x]+1]*p^(I[x]/i0) );
        
        # compute the closure of <pos>.smallerField and <c>
        Add( c, pos.smallerField.root );
        pos.smallerField := Field(c);
        
        # if the degree is <k>,  stop
        if pos.smallerField.degree = pos.k  then
            pos.isSmaller := false;
            InfoRecSL2( "#I  <G> is not definable over a smaller field, ",
                        "char polynomial failed\n" );
        fi;
    fi;

end;


#############################################################################
##
#V  SporadicGroupsInfo  . . . . . . . . . . information about sporadic groups
##
##  names[i]
##      Atlas name (usable by 'CharTable')
##  outer[i]
##      order of outer automorphism group
##  schur[i]
##      order of the Schur-multiplier
##  order[i]
##      order of the sporadic group
##  orders[i]
##      possible element orders of the automorphisms group
##  minDeg0[i]
##      minimal degree of a faithful (projective) representation in char 0 if
##      known, one otherwise.
##  minDegP[i][p]
##      minimal degree of a faithful (proj) representation in char  p | order
##      if known, one otherwise
##
##  The following  program will create a file containing 'SporadicGroupsInfo'
##  by using the information stored in the ordinary and modular atlas in gap.
##  The file name is "sporadic.inf".
##
##  SporadicGroupInfo := [
##    [ "M11", 1, 1 ],    [ "M12", 2, 2 ],   [ "M22", 2, 12 ],
##    [ "M23", 1, 1 ],    [ "M24", 1, 1 ],   [ "J2", 2, 2 ],
##    [ "Suz", 2, 6 ],    [ "HS", 2, 2 ],    [ "McL", 2, 3 ],
##    [ "Co3", 1, 1 ],    [ "Co2", 1, 1 ],   [ "Co1", 1, 2 ],
##    [ "He", 2, 1 ],     [ "Fi22", 2, 6 ],  [ "Fi23", 1, 1 ],
##    [ "F3+", 2, 3 ],    [ "HN", 2, 1 ],    [ "Th", 1, 1 ],
##    [ "B", 1, 2 ],      [ "M", 1, 1 ],     [ "J1", 1, 1 ],
##    [ "ON", 2, 3 ],     [ "J3", 2, 3 ],    [ "Ly", 1, 1 ],
##    [ "Ru", 1, 2 ],     [ "J4", 1, 1 ]
##  ];
##  
##  SporadicGroups         := rec();
##  SporadicGroups.names   := [];
##  SporadicGroups.outer   := [];
##  SporadicGroups.schur   := [];
##  SporadicGroups.orders  := [];
##  SporadicGroups.order   := [];
##  SporadicGroups.minDeg0 := [];
##  SporadicGroups.minDegP := [];
##  
##  for i  in [ 1 .. 26 ]  do
##  
##      # read in standard table, get possible element orders
##      name := SporadicGroupInfo[i][1];
##      Print( "#I  Group: ", name, "\n" );
##      tbl   := CharTable( name );
##      ords  := Set( tbl.orders );
##      schur := SporadicGroupInfo[i][3];
##      if SporadicGroupInfo[i][2] = 2  then
##         UniteSet( ords, 2 * ords );
##      fi;
##      tmp := [];
##      for j  in DivisorsInt(schur)  do
##          UniteSet( tmp, j*ords );
##      od;
##      ords := tmp;
##  
##      # Add the name, orders, outer to our collection
##      Add( SporadicGroups.names,  name );
##      Add( SporadicGroups.outer,  SporadicGroupInfo[i][2] );
##      Add( SporadicGroups.orders, ords );
##      Add( SporadicGroups.order,  tbl.order );
##      Add( SporadicGroups.schur,  schur );
##  
##      # Minimial matrix representation in different characteristic
##      min := Set( List( tbl.irreducibles, x->x[1] ) )[2];
##      for d  in Difference( DivisorsInt(schur), [1] )  do
##              if 1 < min  then
##                  tbl := CharTable(ConcatenationString(String(d),".",name));
##                  if tbl = false  then
##                      Print( "#I      ", d, "*", name, " failed\n" );
##                      min := 1;
##                  else
##                      Print( "#I      ", d, "*", name, " ok\n" );
##                      mn2 := Set( List( tbl.irreducibles, x->x[1] ) )[2];
##                      min := Minimum( mn2, min );
##                  fi;
##              fi;
##      od;
##      Add( SporadicGroups.minDeg0, min );
##  
##      # Minimal matrix representation in same characteristic
##      minDegP := [];
##      for p  in Set( Factors( tbl.order ) )  do
##          Unbind(min);
##          for d  in DivisorsInt(schur)  do
##              if not IsBound(min) or 1 < min  then
##                      namep := ConcatenationString( name, "mod", String(p) );
##                      if d = 1  then
##                          tbl := CharTable(namep);
##                          tmp := namep;
##                      else
##                          tmp := ConcatenationString(String(d),".",namep);
##                          tbl := CharTable(tmp);
##                  fi;
##                      if tbl = false  then
##                          Print( "#I      ", tmp, " failed\n" );
##                          min := 1;
##                      else
##                          Print( "#I      ", tmp, " ok\n" );
##                          mn2 := Set( List( tbl.irreducibles, x->x[1] ) )[2];
##                          if d = 1  then
##                              min := mn2;
##                          else
##                              min := Minimum( mn2, min );
##                          fi;
##                      fi;
##                  fi;
##              od;
##              minDegP[p] := min;
##      od;
##      Add( SporadicGroups.minDegP, minDegP );
##  od;
##
##  # Fix unknown entries (see P.Kleidman & M.Liebeck)
##  Fix := [ [ "M11", 5 ],   [ "M12", 6 ],    [ "M22", 6 ],
##       [ "M23", 11 ],  [ "M24", 11 ],   [ "J1", 7 ],
##       [ "J2", 6 ],    [ "J3", 9 ],     [ "J4", 110 ],
##       [ "HS", 20 ],   [ "McL", 21 ],   [ "He", 18 ],
##       [ "Ru", 28 ],   [ "Suz", 12 ],   [ "ON", 31 ],
##       [ "Co1", 24 ],  [ "Co2", 22 ],   [ "Co3", 23 ],
##       [ "Fi22", 27 ], [ "Fi23", 234 ], [ "F3+", 702 ],
##       [ "HN", 56 ],   [ "Ly", 110 ],   [ "Th", 48 ],
##       [ "B", 234 ],   [ "M", 729 ] ];
##  
##  for i  in Fix  do
##      p := Position( SporadicGroups.names, i[1] );
##      for j  in [ 1 .. Length(SporadicGroups.minDegP[p]) ]  do
##          if IsBound(SporadicGroups.minDegP[p][j])  then
##              d := SporadicGroups.minDegP[p][j];
##              if d < i[2]  then
##                  Print("#I  fixing ",i[1],", p=",j,":",d,"->",i[2],"\n" );
##                  SporadicGroups.minDegP[p][j] := i[2];
##              fi;
##          fi;
##      od;
##  od;
##  
##  # print <SporadicGroups> to file "sporadic.inf"
##  SizeScreen( [ 70, 24 ] );
##  PrintTo("sporadic.inf","SporadicGroupsInfo := ",SporadicGroups,";\n");
##  SizeScreen( [ 80, 24 ] );
##  
SporadicGroupsInfo := rec(
  names := [ "M11", "M12", "M22", "M23", "M24", "J2", "Suz", "HS", 
      "McL", "Co3", "Co2", "Co1", "He", "Fi22", "Fi23", "F3+", 
      "HN", "Th", "B", "M", "J1", "ON", "J3", "Ly", "Ru", "J4" ],
  outer := [ 1, 2, 2, 1, 1, 2, 2, 2, 2, 1, 1, 1, 2, 2, 1, 2, 2, 1, 
      1, 1, 1, 2, 2, 1, 1, 1 ],
  schur := [ 1, 2, 12, 1, 1, 2, 6, 2, 3, 1, 1, 2, 1, 6, 1, 3, 1, 1, 
      2, 1, 1, 3, 3, 1, 2, 1 ],
  orders := 
   [ [ 1, 2, 3, 4, 5, 6, 8, 11 ], [ 1, 2, 3, 4, 5, 6, 8, 10, 11, 
          12, 16, 20, 22, 24, 32, 40, 44 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 18, 20, 
          21, 22, 24, 28, 30, 32, 33, 36, 40, 42, 44, 48, 56, 60, 
          64, 66, 72, 84, 88, 96, 120, 132, 144, 168, 192, 264 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 11, 14, 15, 23 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 14, 15, 21, 23 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 15, 16, 20, 24, 28, 30, 
          32, 40, 48, 60 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 
          20, 21, 22, 24, 26, 27, 28, 30, 32, 33, 36, 39, 40, 42, 
          44, 45, 48, 52, 54, 56, 60, 63, 66, 72, 78, 80, 84, 90, 
          96, 108, 120, 126, 132, 144, 156, 168, 180, 216, 240, 
          252, 288 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 14, 15, 16, 20, 22, 24, 
          28, 30, 32, 40, 44, 48, 60, 80 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 18, 20, 
          21, 22, 24, 27, 28, 30, 33, 36, 42, 45, 48, 54, 60, 66, 
          72, 84, 90, 180 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 18, 20, 21, 
          22, 23, 24, 30 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 18, 20, 
          23, 24, 28, 30 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 
          20, 21, 22, 23, 24, 26, 28, 30, 32, 33, 35, 36, 39, 40, 
          42, 44, 46, 48, 52, 56, 60, 66, 70, 72, 78, 80, 84, 120 ],
      [ 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14, 15, 16, 17, 20, 21, 24, 
          28, 30, 34, 42, 56 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 18, 
          20, 21, 22, 24, 26, 27, 28, 30, 32, 33, 36, 39, 40, 42, 
          44, 45, 48, 52, 54, 56, 60, 63, 64, 66, 72, 78, 80, 84, 
          88, 90, 96, 108, 120, 126, 132, 144, 156, 168, 180, 192, 
          216, 240, 252, 264, 288, 360 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
          18, 20, 21, 22, 23, 24, 26, 27, 28, 30, 35, 36, 39, 42, 
          60 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
          18, 20, 21, 22, 23, 24, 26, 27, 28, 29, 30, 32, 33, 34, 
          35, 36, 39, 40, 42, 44, 45, 46, 48, 51, 52, 54, 56, 58, 
          60, 63, 66, 69, 70, 72, 78, 81, 84, 87, 90, 96, 99, 102, 
          105, 108, 117, 120, 126, 132, 135, 138, 144, 156, 162, 
          168, 174, 180, 198, 210, 216, 234, 252, 270, 360 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 18, 19, 
          20, 21, 22, 24, 25, 28, 30, 35, 38, 40, 42, 44, 50, 60, 
          70, 80 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 18, 19, 20, 
          21, 24, 27, 28, 30, 31, 36, 39 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
          18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 30, 31, 32, 
          33, 34, 35, 36, 38, 39, 40, 42, 44, 46, 47, 48, 50, 52, 
          54, 55, 56, 60, 62, 64, 66, 68, 70, 72, 76, 78, 80, 84, 
          88, 92, 94, 96, 104, 110, 112, 120, 132, 140 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 
          18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 
          32, 33, 34, 35, 36, 38, 39, 40, 41, 42, 44, 45, 46, 47, 
          48, 50, 51, 52, 54, 55, 56, 57, 59, 60, 62, 66, 68, 69, 
          70, 71, 78, 84, 87, 88, 92, 93, 94, 95, 104, 105, 110, 
          119 ], [ 1, 2, 3, 5, 6, 7, 10, 11, 15, 19 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 18, 19, 
          20, 21, 22, 24, 28, 30, 31, 32, 33, 36, 38, 40, 42, 45, 
          48, 56, 57, 60, 62, 66, 72, 84, 90, 93, 96, 114, 120, 
          168, 186 ], 
      [ 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15, 16, 17, 18, 19, 20, 24, 
          27, 30, 34, 36, 38, 45, 48, 51, 54, 57, 60, 72, 90, 102, 
          114 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 18, 20, 21, 
          22, 24, 25, 28, 30, 31, 33, 37, 40, 42, 67 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13, 14, 15, 16, 20, 24, 26, 
          28, 29, 30, 32, 40, 48, 52, 58 ], 
      [ 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 14, 15, 16, 20, 21, 22, 
          23, 24, 28, 29, 30, 31, 33, 35, 37, 40, 42, 43, 44, 66 ] ]
   ,
  order := [ 7920, 95040, 443520, 10200960, 244823040, 604800, 
      448345497600, 44352000, 898128000, 495766656000, 
      42305421312000, 4157776806543360000, 4030387200, 
      64561751654400, 4089470473293004800, 
      1255205709190661721292800, 273030912000000, 90745943887872000,
      4154781481226426191177580544000000, 
      808017424794512875886459904961710757005754368000000000, 
      175560, 460815505920, 50232960, 51765179004000000, 
      145926144000, 86775571046077562880 ],
  minDeg0 := [ 10, 10, 10, 22, 23, 6, 12, 22, 22, 23, 23, 24, 51, 
      78, 782, 783, 133, 248, 4371, 196883, 56, 342, 18, 2480, 28, 
      1333 ],
  minDegP := [ [ , 10, 5,, 10,,,,,, 9 ], [ , 6, 6,, 10,,,,,, 10 ], 
      [ , 6, 6,, 10,, 10,,,, 10 ], 
      [ , 11, 22,, 22,, 22,,,, 22,,,,,,,,,,,, 21 ], 
      [ , 11, 22,, 23,, 23,,,, 23,,,,,,,,,,,, 23 ], 
      [ , 6, 6,, 6,, 6 ], [ , 12, 12,, 12,, 12,,,, 12,, 12 ], 
      [ , 20, 22,, 21,, 22,,,, 22 ], [ , 22, 21,, 21,, 22,,,, 22 ], 
      [ , 23, 23,, 23,, 23,,,, 23,,,,,,,,,,,, 23 ], 
      [ , 22, 23,, 23,, 23,,,, 23,,,,,,,,,,,, 23 ], 
      [ , 24, 24,, 24,, 24,,,, 24,, 24,,,,,,,,,, 24 ], 
      [ , 51, 51,, 51,, 50,,,,,,,,,, 51 ], 
      [ , 27, 27,, 27,, 78,,,, 78,, 78 ], 
      [ , 234, 234,, 782,, 782,,,, 782,, 782,,,, 234,,,,,, 782 ], 
      [ , 702, 702,, 702,, 702,,,, 783,, 702,,,, 702,,,,,, 783,,,,,
          , 702 ], [ , 56, 56,, 56,, 133,,,, 133,,,,,,,, 133 ], 
      [ , 48, 48,, 48,, 48,,,,,, 248,,,,,, 48,,,,,,,,,,,, 248 ], 
      [ , 234, 234,, 234,, 234,,,, 4371,, 234,,,, 4371,, 234,,,, 
          4371,,,,,,,, 234,,,,,,,,,,,,,,,, 234 ], 
      [ , 729, 729,, 729,, 729,,,, 729,, 729,,,, 196883,, 196883,,,
          , 196883,,,,,, 729,, 196883,,,,,,,,,, 729,,,,,, 729,,,,,,,
          ,,,,, 729,,,,,,,,,,,, 729 ], 
      [ , 20, 56,, 56,, 31,,,, 7,,,,,,,, 22 ], 
      [ , 153, 31,, 342,, 31,,,, 31,,,,,,,, 342,,,,,,,,,,,, 31 ], 
      [ , 9, 9,, 18,,,,,,,,,,,, 18,, 18 ], 
      [ , 110, 110,, 110,, 2480,,,, 2480,,,,,,,,,,,,,,,,,,,, 2480,,,
          ,,, 110,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, 110 ], 
      [ , 28, 28,, 28,, 28,,,,,, 28,,,,,,,,,,,,,,,, 28 ], 
      [ , 110, 110,, 1333,, 1333,,,, 110,,,,,,,,,,,, 110,,,,,, 110,
          , 110,,,,,, 1333,,,,,, 110 ] ] );


#############################################################################
##
#F  RecSL.SetSporadicGroups( <pos> )  . . . . possible groups in PGL(<d>,<q>)
##
##  'SetSporadicGroups' sets the component '<pos>.sporadicGroups' to  a  list
##  of integers (corresponding to groups in the list  'SporadicGroupsInfo.*')
##  describing the possible sporadic groups in PGL(<n>,<q>).
##
RecSL.SetSporadicGroups := function( pos )
    local   l,  p;

    # at first check the order of the sporadic groups
    if IsBound(pos.orderGL)  then
        l := Filtered( [1..26], x ->
                       pos.orderGL mod SporadicGroupsInfo.order[x] = 0
                       and SporadicGroupsInfo.order[x] < pos.q3d );
    else
        l := Filtered( [1..26], x -> SporadicGroupsInfo.order[x] < pos.q3d );
    fi;
    if 0 = Length(l)  then
        InfoRecSL2( "#I  <G> is not a sporadic group,  ",
                    "group order criteria failed\n" );
    fi;

    # check possible minimal degrees
    if 0 < Length(l)  then
        p := SmallestRootInt( pos.q );
        pos.sporadicGroups := Filtered( l, 
                function(x)
                   if SporadicGroupsInfo.order[x] mod p = 0  then
                       return SporadicGroupsInfo.minDegP[x][p] <= pos.d+1;
                   else
                       return SporadicGroupsInfo.minDeg0[x] <= pos.d+1;
                   fi;
                end );
        if 0 = Length(pos.sporadicGroups)  then
            InfoRecSL2( "#I  <G> is not a sporadic group,  ",
                        "minimal degree criteria failed\n" );
        fi;
    else
        pos.sporadicGroups := l;
    fi;

    # are same sporadic groups still possible?
    pos.isSporadic := 0 < Length(pos.sporadicGroups);

end;


#############################################################################
##
#F  RecSL.CheckSporadicGroups( <pos>, <pord> )  . . . . .  element order test
##
RecSL.CheckSporadicGroups := function( pos, pord )
    
    # sporadic groups still possible?
    if not pos.isSporadic  then return;  fi;
    pos.statistic[RecSL.STAT_SPORADIC] := pos.tries;

    # check the possible element orders of the automorphism groups
    pos.sporadicGroups := Filtered( pos.sporadicGroups, x -> pord in
                                    SporadicGroupsInfo.orders[x] );

    # if the list is trivial reset '<pos>.isSporadic'
    if 0 = Length(pos.sporadicGroups)  then
        pos.isSporadic := false;
        InfoRecSL2( "#I  <G> is not a sporadic group,  ",
                    "element order criteria failed\n" );
    fi;

    # if the exponent is too big stop
    if pos.exponent >= pos.q3d  then
        pos.isSporadic := false;
        InfoRecSL2( "#I  <G> is not a sporadic group, ",
                    "group exponent criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSL.SetTensorPowers( <pos> )  . . . . . . . . .  possible tensor powers
##
##  d =  e^f, G < GL(e,q)  wr Sym(f).  Assume  that there are  no semi primes
##  smaller then f+1.
##
RecSL.SetTensorPowers := function( pos )
    local   e,  f,  k,  d,  exp;

    # no reduction is possible if f is 1
    e := SmallestRootInt( pos.d );
    if e = pos.d  then
        pos.isTensorPower := false;
        InfoRecSL2( "#I  <G> is no tensor power,  dimension is ",
                    "not a power\n" );
        
    # store the various pairs <e>, <f>
    else
        pos.isTensorPower := true;
        f := LogInt( pos.d, e );
        d := Filtered( DivisorsInt(f), x -> x <> f );
        pos.dimsTensorPowers := List( d, x -> [ e^x, f/x ] );
    fi;

end;  


#############################################################################
##
#F  RecSL.CheckTensorPowers( <pos>, <ord> ) . . . . . . . check element order
##
##  <ord> is the order of the  group element A. If A lies in GL(e,q) X Sym(f)
##  then gcd(<ord>,exp(GL(e,q))) must be a valid element order in Sym(f).
##        
##  Assume that all integer <= <pos.d> are no semi primes.
##
RecSL.CheckTensorPowers := function( pos, ord )
    local   new,  p,  m;
    
    # are tensor powers still possible?
    if not pos.isTensorPower  then return;  fi;
    pos.statistic[RecSL.STAT_POWER] := pos.tries;

    # check exponents
    new := [];
    for p  in pos.dimsTensorPowers  do
        m := ord / Gcd( ord, pos.expsGL[p[1]] );
        if 2*Factorial(p[2]) mod m = 0
           and Sum(Collected(Factors(m)),x->x[1]^x[2]) <= p[2]
        then
            Add( new, p );
        fi;
    od;
    pos.dimsTensorPowers := new;
    
    # if <new> is trivial no reduction is possible
    if 0 = Length(new)  then
        pos.isTensorPower := false;
        InfoRecSL2( "#I  <G> is no tensor power,  element order ",
                    "criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSL.SetTensorProducts( <pos> )  . . . . . . .  possible tensor products
##
##  G < GL(d1,q) x GL(d2,q) with 1 < d1 < d2 < d and d1*d2 = d.
##
RecSL.SetTensorProducts := function( pos )
    local   i,  exp;

    # if <pos>.d is a prime, no tensor products are possible
    if IsPrime(pos.d)  then
        pos.isTensorProduct := false;
        InfoRecSL2( "#I  <G> is no tensor product,  dimension is ",
                    "a prime\n" );
        return;
    else
        pos.isTensorProduct := true;
    fi;
        
    # compute the exponents of the tensor products
    pos.expsTensorProducts := [];
    for i  in DivisorsInt(pos.d)  do
        if 1 < i and i^2 < pos.d  then
            exp := LcmInt( pos.expsGL[i], pos.expsGL[pos.d/i] );
            Add( pos.expsTensorProducts, [ exp, i ] );
        fi;
    od;

end;        


#############################################################################
##
#F  RecSL.CheckTensorProducts( <pos>, <ord> ) . . . . . . check element order
##
RecSL.CheckTensorProducts := function( pos, ord )

    # are tensor products still possible?
    if not pos.isTensorProduct  then return;  fi;
    pos.statistic[RecSL.STAT_PRODUCT] := pos.tries;

    # check if <ord> divides any exponent of a tensor product
    pos.expsTensorProducts := Filtered( pos.expsTensorProducts,
                                        x -> x[1] mod ord = 0 );
    
    # if <pos>.expsTensorProducts is empty no tps are possible
    if 0 = Length(pos.expsTensorProducts)  then
        pos.isTensorProduct := false;
        InfoRecSL2( "#I  <G> is no tensor product,  element order ",
                    "criteria failed\n" );
    fi;

end;


#############################################################################
##

#F  RecSL.InvariantForm( <pos> )  . . . . . . . . .  find the invariant forms
##
RecSL.InvariantForm := function( pos )
    local   clas,  orthogonal,  form;

    # try to find a form
    if (pos.isOrthogonal or pos.isSymplectic) and pos.isUnitary  then
        if IsBound(pos.unitaryForm)  then
            clas := ClassicalForms( pos.group, "dual" );
        elif IsBound(pos.symmetricForm) or IsBound(pos.symplecticForm)  then
            clas := ClassicalForms( pos.group, "frobenius" );
        else
            clas := ClassicalForms( pos.group, "dual", "frobenius" );
        fi;
    elif (pos.isOrthogonal or pos.isSymplectic) and not pos.isUnitary  then
        clas := ClassicalForms( pos.group, "dual" );
    elif not (pos.isOrthogonal or pos.isSymplectic) and pos.isUnitary  then
        clas := ClassicalForms( pos.group, "frobenius" );
    else
        return;
    fi;

    # check which forms where found
    orthogonal := ["orthogonalminus","orthogonalplus","orthogonalcircle"];
    for form  in clas  do
        if form[1] = "symplectic"  then
            pos.symplecticForm    := form[2];
            pos.symplecticScalars := form[3];
            if pos.p = 2  then
                pos.symmetricForm    := form[2];
                pos.symmetricScalars := form[3];
            fi;
        elif form[1] = "unitary"  then
            pos.unitaryForm    := form[2];
            pos.unitaryScalars := form[3];
        elif form[1] in orthogonal  then
            pos.symmetricForm    := form[2];
            pos.symmetricScalars := form[3];
            pos.quadraticForm := form[4];
            if 4 < Length(form)  then
                pos.squareDiscriminat := form[5];
            fi;
            if pos.p = 2  then
                pos.symplecticForm    := form[2];
                pos.symplecticScalars := form[3];
            fi;
            if form[1] = "orthogonalminus"  then
                pos.signum := -1;
            elif form[1] = "orthogonalplus"  then
                pos.signum := +1;
            else
                pos.signum := 0;
            fi;
        fi;
    od;
                
end;


#############################################################################
##

#F  RecSL.Print . . . . . . . . . . . . . . . . . . . . . . nice pretty print
##
RecSL.Print := function( obj )
    local   name;

    if obj.printLevel = 0  then
        Print( "<< SL recognition record >>" );
        return;
    fi;

    Print( "#I  field: ", obj.q, ", dimension: ", obj.d,
           ", number of generators: ", Length(obj.generators), "\n" );
    if IsBound(obj.isChevalley) and obj.isChevalley  then
        Print( "#I  <G> could be almost simple" );
        if 1 < obj.printLevel  then
            Print( ": " );
            for name in obj.expsChev  do
                Print(name[2],"_",name[3],"(",name[4],"^",name[5],") ");
            od;
        fi;
        Print( "\n");
    fi;
    
    # sporadic groups
    if IsBound(obj.isSporadic) and obj.isSporadic  then
        Print( "#I  <G> could be an almost sporadic group" );
        if 1 < obj.printLevel  then
            Print( ": " );
            for name in SporadicGroupsInfo.names{obj.sporadicGroups}  do
                Print( name, " " );
            od;
        fi;
        Print( "\n");
    fi;

    if IsBound(obj.isAlternating) and obj.isAlternating  then
        Print( "#I  <G> could be an almost alternating group" );
        if 1 < obj.printLevel  then
            Print( ": ");
            for name in obj.alternating  do
                Print( "A", name, " " );
            od;
        fi;
        Print( "\n" );
    fi;
    if IsBound(obj.isClassical) and obj.isClassical  then
        Print( "#I  <G> could be a classical group" );
        if 1 < obj.printLevel  then
            Print( ": " );
            if obj.isUnitary  then
                Print( "SU " );
            fi;
            if obj.isSymplectic  then
                Print( "SP " );
            fi;
            if ( obj.isSymplectic and obj.p = 2 ) or obj.isOrthogonal   then
                Print( "Omega " );
            fi;
        fi;
        Print( "\n" );
    fi;
    if IsBound(obj.isLarger) and obj.isLarger  then
        Print( "#I  <G> could be definable over a larger field\n" );
    fi;
    if IsBound(obj.isSmaller) and obj.isSmaller  then
        Print( "#I  <G> could be definable over ", obj.smallerField, "\n" );
    fi;
    if IsBound(obj.isMysteriousP) and obj.isMysteriousP  then
       Print( "#I  <G> could be acting on a mysterious p-group\n" );
    fi;
    if IsBound(obj.isImprimitive) and obj.isImprimitive  then
        Print( "#I  <G> could be imprimitive\n" );
    fi;
    if IsBound(obj.isReducible) and obj.isReducible  then
        Print( "#I  <G> could be reducible" );
        if 1 < obj.printLevel  then
            Print( ": " );
            if 0 = Length(obj.dimsReducible)  then
                Print( "1 .. ", obj.d-1 );
            else
                for name  in obj.dimsReducible  do
                    if name < obj.d  then
                        Print( name, " " );
                    fi;
                od;
            fi;
        fi;
        Print( "\n" );
    fi;
    if IsBound(obj.isTensorPower) and obj.isTensorPower  then
        Print( "#I  <G> could be a tensor power\n" );
    fi;
    if IsBound(obj.isTensorProduct) and obj.isTensorProduct  then
        Print( "#I  <G> could be a tensor product\n" );
    fi;
    if IsBound(obj.isSL) and obj.isSL  then
        Print( "#I  <G> is SL( ", obj.d, ", ", obj.q, " )\n" );
    elif IsBound(obj.isGL) and obj.isGL  then
        Print( "#I  <G> is GL( ", obj.d, ", ", obj.q, " )\n" );
    elif IsBound(obj.containsSL) and obj.containsSL  then
        Print( "#I  <G> is SL( ",obj.d,", ",obj.q," ) . ",obj.ext,"\n" );
    else
        Print( "#I  <G> could still contain SL( ",obj.d,", ",obj.q," )\n" );
    fi;
    Print( "<< SL recognition record >>" );
end;


#############################################################################
##
#F  RecSL.SetPrintLevel( <obj>, <lev> ) . . . . . . . . . . . set print level
##
RecSL.SetPrintLevel := function( obj, lev )
    if 2 < lev or lev < 0  then
        Error( "invalid print level" );
    fi;
    obj.printLevel := lev;
end;


#############################################################################
##
#F  RecSL.AddGenerator( <pos>, <mat> )  . . . . . . . . . add a new generator
##
RecSL.AddGenerator := function( pos, mat )

    # add <mat> to the generator
    Add( pos.generators, mat );

    # update the random elements
    if IsBound(pos.randomSeed)  then
        pos.randomSeed := pos.randomSeed * mat;
        Add( pos.randomSeed, mat );
    fi;

    # if <pos> contains SL update information
    if pos.containsSL  then
        if not pos.isGL  then
            Unbind(pos.isGL);
            Unbind(pos.ext);
        fi;
        RecognizeSL( pos, 0 );
    fi;
    
end;


#############################################################################
##
#F  RecSL.Setup( <G>, <gens> )  . . . . . . . . . . . . . . set up pos record
##
RecSL.Setup := function( G, gens )

    local   F,                  # field
            R,                  # full polynomial ring over <F>
            q,                  # size of <F>
            p,                  # char of <F>
            qi,                 # q^i
            pos,                # possible answers
            timer,              # runtime storage
            i, j;               # loop and tmp

    # get the field, its size and characteristic
    F := G.field;
    R := PolynomialRing(F);
    q := Size(F);
    p := F.char;

    # enter important information about the group
    pos            := rec();
    pos.p          := p;
    pos.k          := LogInt(q, p);
    pos.q          := q;
    pos.d          := Length(G.identity);
    pos.qq         := p^QuoInt(pos.k, 2);
    pos.q3d        := q^(3*pos.d);
    pos.field      := F;
    pos.group      := G;
    pos.identity   := G.identity;
    pos.generators := ShallowCopy(gens);
    pos.semiPrimes := [];
    pos.statistic  := [ 1 .. 11 ] * 0;
    pos.printLevel := 1;
    pos.isRecSL    := true;
    pos.containsSL := false;
    pos.operations := RecSL;
    pos.setupTime  := -Runtime();
    InfoRecSL1( "#I  field: ", pos.q, ", dimension: ", pos.d, 
                ", number of generators: ", Length(pos.generators), "\n" );

    # compute the exponent of GL(i,q) for 1 <= i < d, i|d
    timer  := Runtime();
    pos.expsGL := [ q-1 ];
    j  := q-1;
    qi := q;
    i  := 2;
    if 1 < pos.d  then
        while pos.d mod i <> 0  do i := i + 1;  od;
        for i  in [ 2 .. pos.d/i ]  do
            qi := qi * q;
            j  := LcmInt( j, qi-1 );
            if pos.d mod i = 0  then
                pos.expsGL[i] := j * p^(1+LogInt(i-1, p));
            fi;
        od;
    fi;
    InfoRecSL4( "#I  exponents: ", Runtime()-timer, " msec\n" );

    # compute the order of GL( d, q ) for small dimensions
    if pos.d < 10  then
        timer := Runtime();
        pos.orderGL := 1;
        qi := 1;
        for i  in [ 1 .. pos.d ]  do
            qi := qi * q;
            pos.orderGL := pos.orderGL * (qi-1);
        od;
        pos.orderGL := pos.orderGL * q^(pos.d*(pos.d-1)/2);
        InfoRecSL4( "#I  order gl: ", Runtime()-timer, " msec\n" );
    fi;

    # set the possible groups contained in <G>
    timer := Runtime();
    pos.operations.SetAlternating     ( pos );
    pos.operations.SetChevalley       ( pos );
    pos.operations.SetClassicalGroups ( pos );
    pos.operations.SetLargerField     ( pos );
    pos.operations.SetMysteriousPGroup( pos );
    pos.operations.SetImprimitive     ( pos );
    pos.operations.SetReducible       ( pos );
    pos.operations.SetSmallerField    ( pos );
    pos.operations.SetSporadicGroups  ( pos );
    pos.operations.SetTensorPowers    ( pos );
    pos.operations.SetTensorProducts  ( pos );
    InfoRecSL4( "#I  subgroup setup: ", Runtime()-timer, " msec\n" );
    pos.setupTime := Runtime() + pos.setupTime;

    # exponent is at least 1
    pos.exponent := 1;
    pos.orders   := [];

    # and return
    return pos;

end;


#############################################################################
##
#F  RecognizeSL( <G>, <gens>, <tries> ) . . . . . . . . . .  regonize SL(n,q)
##
RecognizeSL := function( arg )

    local   G,                  # the group or a supergroup
            gens,               # generators
            tries,              # number of tries
            p,                  # o semi prime
            A,                  # random matrix of <G>
            o,                  # order of <A> = <oo> * <so>
            oe,                 # order of a finite field element
            oo,                 # semi prime free part of <o>
            so,                 # semi prime part of <o>
            po,                 # projective order of <A>
            poo,                # projective semi prime free part of <po>
            pso,                # projective semi prime part of <po>
            cpol,               # characteristic polynomial of <A>
            facs,               # factors of <cpol>
            mpol,               # minimal polynomial of <A>
            fac,                # factor from <mpol>
            pos,                # possible answers
            quo,                # quotient of <cpol>
            tmp;

    # check the arguments
    if 2 = Length(arg)  then
        G     := arg[1];
        gens  := Set(G.generators);
        tries := arg[2];
    elif 3 = Length(arg)  then
        G     := arg[1];
        gens  := arg[2];
        tries := arg[3];
    else
        Error( "usage: RecognizeSL( <G>, <gens>, <tries> )" );
    fi;
    if not IsInt(tries)  then
        Error( "<tries> must be an integer" );
    fi;
    if IsBound(G.isRecSL) and 3 = Length(arg)  then
        Error( "use 'AddGenerator' to add a new generator" );
    fi;

    # if the group is trivial return
    if IsBound(G.isRecSL) and G.isRecSL  then
        if 0 = Length(G.generators)  then
            return G;
        fi;
    elif 0 = Length(gens)  then
        pos := RecSL.Setup( G, gens );
        if pos.d = 1  then
            pos.isAlternating   := false;
            pos.isChevalley     := false;
            pos.isImprimitive   := false;
            pos.isSporadic      := false;
            pos.isTensorProduct := false;
            pos.isTensorPower   := false;
            pos.isLarger        := false;
            pos.isSmaller       := false;
            pos.isMysteriousP   := false;
            pos.isReducible     := false;
            pos.isClassical     := false;
            pos.containsSL      := true;
            pos.ext             := 1;
            pos.isSL            := true;
            pos.isGL            := pos.q = 2;
        else
            pos.containsSL := false;
            pos.isSL       := false;
            pos.isGL       := false;
        fi;
        return pos;
    fi;

    # if <d> = 1  we are done
    if IsBound(G.isRecSL) and G.isRecSL  then
        if G.d = 1  then
            return G;
        fi;

    # set up dummy record
    elif Length(G.identity) = 1  then
        pos            := rec();
        pos.d          := Length(G.identity);
        pos.q          := Size(G.field);
        pos.field      := G.field;
        pos.generators := ShallowCopy(gens);
        pos.statistic  := [ 1 .. 10 ] * 0;
        pos.printLevel := 1;
        pos.containsSL := true;
        pos.isRecSL    := true;
        pos.operations := RecSL;

        # check which subgroup if F* it is
        o := List( pos.generators, x -> DeterminantMat(x) );
        pos.isSL := ForAll( o, x -> x = pos.field.one );
        if not pos.isSL  then
            pos.ext  := Lcm(List( o, x -> Order( pos.field, x ) ));
            pos.isGL := pos.ext = pos.q-1;
        else
            pos.isGL := false;
        fi;

        # and return
        return pos;
    fi;

    # set up pos recgonition record
    if IsBound(G.isRecSL) and G.isRecSL  then
        pos := G;
        G   := G.group;
    else
        pos := RecSL.Setup( G, gens );
    fi;

    # initialize orders
    o := 0;  po := 0;  poo := 0;  pso := 0;  cpol := 0;

    # start trying random non-one group elements
    pos.runTime := -Runtime();
    pos.tries   := 0;
    while not pos.containsSL and pos.tries < tries  do

        # find a non-trivial random element of <G>
        repeat
            A := pos.operations.Random(pos);
        until A <> pos.identity;
        pos.tries := pos.tries + 1;
        InfoRecSL2("#I  trying ", pos.tries, ".th element of <G>\n");

        # compute the minimal polynomial of <A> and the projective order
        if    pos.isChevalley
           or pos.isAlternating
           or pos.isImprimitive
           or pos.isLarger
           or pos.isMysteriousP
           or pos.isSmaller
           or pos.isSporadic
           or pos.isTensorPower
           or pos.isTensorProduct
           or ( pos.isOrthogonal and 2 = pos.d )
        then
            mpol := MinimalPolynomial( FiniteFieldMatrices, A );
            mpol.baseRing := pos.field;
            po := pos.operations.OrderScalar( pos, mpol );
            oe := OrderFFE(po[2]);
            po := po[1];
            InfoRecSL2( "#I  projective order = ", 
                    SemiPrimePowersInt(po,pos.p,pos.k*pos.d), "\n" );

            # remove semi primes from order
            poo := po;
            pso := 1;
            for p in pos.semiPrimes  do
                while poo mod p = 0  do poo := poo / p;  pso := pso * p;  od;
            od;

            # <o> is the order, <so> the semi order part, <o> = <so> * <oo>
            o  := po  * oe;
            oo := poo * oe;
            so := pso;

            # the exponent is at least the lcm of <exp> and <oo>
            pos.exponent := LcmInt( pos.exponent, oo );
            Add( pos.orders, [ poo, pso, oe ] );
        else
            mpol := false;
        fi;

        # compute the characteristic polynomial of <A>
        if pos.isReducible or pos.isClassical or pos.isSmaller  then

            # if the minimal polynomial has degree <pos.d> it is the char
            if mpol <> false and Degree(mpol) = pos.d  then
                cpol := mpol;

            # otherwise we have to compute it
            else
                cpol := CharacteristicPolynomial( FiniteFieldMatrices, A );
                cpol.baseRing := pos.field;
                if mpol <> false and pos.isReducible  then
                    quo  := cpol;
                    facs := [];
                    for fac in Set(Factors(mpol))  do
                        tmp := Quotient( quo, fac );
                        while tmp <> false  do
                            Add( facs, fac );
                            quo := tmp;
                            tmp := Quotient( quo, fac );
                        od;
                    od;
                    Sort(facs);
                    cpol.factors := facs;
                fi;
            fi;
        else
            cpol := false;
        fi;

        # check all possibilities
        pos.operations.CheckAlternating     ( pos,       poo, pso    );
        pos.operations.CheckChevalley       ( pos,       poo, pso, A );
        pos.operations.CheckClassicalGroups ( pos, cpol, poo,      A );
        pos.operations.CheckLargerField     ( pos,       po,       A );
        pos.operations.CheckMysteriousPGroup( pos,       poo         );
        pos.operations.CheckImprimitive     ( pos,       po          );
        pos.operations.CheckReducible       ( pos, cpol              );
        pos.operations.CheckSmallerField    ( pos, cpol, po          );
        pos.operations.CheckSporadicGroups  ( pos,       po          );
        pos.operations.CheckTensorPowers    ( pos,       po          );
        pos.operations.CheckTensorProducts  ( pos,       po          );

        # set 'containsSL'
        pos.containsSL :=     not pos.isReducible
                          and not pos.isImprimitive
                          and not pos.isMysteriousP
                          and not pos.isTensorProduct
                          and not pos.isTensorPower
                          and not pos.isLarger
                          and not pos.isSmaller
                          and not pos.isAlternating
                          and not pos.isChevalley
                          and not pos.isSporadic
                          and not pos.isClassical;
    od;

    # set 'containsSL'
    pos.containsSL :=     not pos.isReducible
                      and not pos.isImprimitive
                      and not pos.isMysteriousP
                      and not pos.isTensorProduct
                      and not pos.isTensorPower
                      and not pos.isLarger
                      and not pos.isSmaller
                      and not pos.isAlternating
                      and not pos.isChevalley
                      and not pos.isSporadic
                      and not pos.isClassical;
    pos.runTime := Runtime() + pos.runTime;

    # if <G> contains the special linear group compute the determinants
    pos.postTime := -Runtime();
    if pos.containsSL  then
        if pos.q = 2  then
            pos.isSL := true;
            pos.isGL := true;
            pos.ext  := 1;
        else
            o        := List( pos.generators, x -> DeterminantMat(x) );
            pos.isSL := ForAll( o, x -> x = pos.field.one );

            # if <G> is not SL,  compute | <G> / SL |
            if not pos.isSL  then
                pos.ext  := Lcm(List( o, x -> Order( pos.field, x ) ));
                pos.isGL := pos.ext = pos.q-1;
            else
                pos.isGL := false;
                pos.ext  := 1;
            fi;
        fi;
        pos.orderSL := 1;
        o := 1;
        for tmp  in [ 1 .. pos.d ]  do
            o := o * pos.q;
            pos.orderSL := pos.orderSL * (o-1);
        od;
        pos.orderSL := pos.orderSL * pos.q^(pos.d*(pos.d-1)/2) / (pos.q-1);
        pos.size := pos.ext * pos.orderSL;
    else
        pos.isSL := false;
        pos.isGL := false;
    fi;
    pos.postTime := Runtime() + pos.postTime;

    # return what we have found so far
    return pos;

end;

RecogniseSL := RecognizeSL;
