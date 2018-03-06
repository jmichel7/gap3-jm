#############################################################################
##
#A  chevgrp.g                   GAP library                      Frank Celler
##
#H  @(#)$Id: chevgrp.g,v 1.1 1997/03/10 13:49:08 gap Exp $
##
#Y  Copyright (C) 1995,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains information about the Chevalley-groups.  They are used
##  in "recsl.g" to exclude subgroups of GL(q,n).
##
Revision_chevgrp_g :=
    "@(#)$Id: chevgrp.g,v 1.1 1997/03/10 13:49:08 gap Exp $";


#############################################################################
##
#F  ChevOps . . . . . . . . . . . . . . . . . . . . . . . . operations record
##
ChevOps := rec();

ChevOps.Print := function(obj) Print(obj.name);  end;
ChevOps.\= := function(a,b) return a.name = b.name;  end;


#############################################################################
##
##  'ChevX.lowerBound( <n>, <p>, <k> )':
##      gives  an  lower  bound  for  the dimension of the  representation of
##      'ChevX' over  a prime other than <p>. See Vicente Landazuri and  Gary
##      M. Seitz, "On the  Minimal Degrees of Projective presentations of the
##      Finite  Chevalley  Groups"  and   P.Kleiderman  and  M.Liebeck,  "The
##      Subgroup Structure of the finite Classical groups".
##
##  'ChevX.order( <n>, <p>, <k> )':
##      gives the order  of the adjoint  Chevalley  group. Returns 'false' if
##      the group is not simple for the given parameter.
##
##  'ChevX.automorphism( <n>, <p>, <k> )':
##      gives the order of the outer automorphism group
##
##  'ChevX.multiplier( <n>, <p>, <k> )':
##      gives the order of the Schur multiplier
##
##  'ChevX.uexponent( <n>, <p>, <k> )':
##      gives an multiple of the exponent. This might be the order!
##
##  'ChevX.possibilities( <d>, <m>, <o>, <p> )'
##      returns  a  pair  [l1,l2] of lists such  that l2 contains  triples of
##      parameters for simple Chevalley-X-groups of dimension at most <d> and
##      order  at  most  <m> over a  field with characteristic different than
##      <p>, l1 contains  triples for simple  Chevalley-X-groups of dimension
##      at most <d> for fields with characteristic <p>. The order must divide
##      <o> in both cases.
##


#############################################################################
##

#V  ChevA . . . . . . . . . . . . . . . . . . . . . .  PSL(n+1,q) information
##
##  Exceptions:
##    A_1(2) = S_3
##    A_1(3) = A_4
##
ChevA := rec( name := "ChevA", operations := ChevOps );

ChevA.order := function( n, p, k )
    local   i,  o,  q,  qi;

    if n < 1 or ( n = 1 and k = 1 and ( p = 2 or p = 3 ) )  then
        return false;
    fi;
    q  := p ^ k;
    o  := q ^ (n*(n+1)/2);
    qi := q;
    for i  in [ 2 .. n+1 ]  do
        qi := qi * q;
        o  := o * ( qi - 1 );
    od;
    return o / GcdInt( n+1, q-1 );

end;

ChevA.possibilities := function( n, m, ord, p )
    local   l1,  l2,  t,  o,  d,  exc;

    # at first no possibilities at all
    l1 := [];
    l2 := [];

    # at first characteristic <p>
    t := [1,p,1];
    if p = 2 or p = 3  then
        t[3] := 2;
    fi;
    o := ChevA.order( t[1], t[2], t[3] );
    while o <= m  do
        while o <= m and ord mod o = 0  do
            Add( l1, Copy(t) );
            t[3] := t[3] + 1;
            o := ChevA.order( t[1], t[2], t[3] );
        od;
        t[1] := t[1] + 1;
        t[3] := 1;
        o := ChevA.order( t[1], t[2], t[3] );
    od;

    # exceptions in lower bound formula or order
    exc := [ [1,2,2], [1,3,2], [2,2,1], [2,2,2], [2,2,1], [2,3,1] ];

    # now different characteristic
    t := [1,2,1];
    if p = 2  then
        t[2] := 3;
    fi;
    if t[2] = 2 or t[2] = 3  then
        t[1] := 2;
    fi;
    o := ChevA.order( t[1], t[2], t[3] );
    d := ChevA.lowerBound( t[1], t[2], t[3] );
    while ( t in exc ) or ( o <= m and d <= n )  do
        while ( t in exc ) or ( o <= m and d <= n )  do
            while (t in exc) or (o <= m and d <= n and ord mod o = 0)  do
                if o<=m and d<=n and (not t in exc or ord mod o=0)  then
                    Add( l2, Copy(t) );
                fi;
                t[1] := t[1] + 1;
                o := ChevA.order( t[1], t[2], t[3] );
                d := ChevA.lowerBound( t[1], t[2], t[3] );
            od;
            t[1] := 1;
            t[3] := t[3] + 1;
            o := ChevA.order( t[1], t[2], t[3] );
            d := ChevA.lowerBound( t[1], t[2], t[3] );
        od;
        t[1] := 1;
        t[2] := NextPrimeInt(t[2]);
        t[3] := 1;
        if t[2] = p  then
            t[2] := NextPrimeInt(t[2]);
        fi;
        if t[2] = 3  then
            t[1] := 2;
        fi;
        o := ChevA.order( t[1], t[2], t[3] );
        d := ChevA.lowerBound( t[1], t[2], t[3] );
    od;

    # return the possibilities
    return [ l1, l2 ];

end;

ChevA.lowerBound := function( n, p, k )
    local   q;

    q := p^k;
    if n = 1 and q = 4  then
        return 2;
    elif n = 1 and q = 9  then
        return 3;
    elif n = 2 and q = 2  then
        return 2;
    elif n = 2 and q = 4  then
        return 4;
    elif n = 1  then
        return ( q - 1 ) / GcdInt( 2, q - 1 );
    else
        return q^n - 1;
    fi;

end;

ChevA.automorphism := function( n, p, k )

    if n = 1  then
        return GcdInt( 2, p^k - 1 ) * k;
    else
        return GcdInt( n+1, p^k-1 ) * k * 2;
    fi;

end;

ChevA.multiplier := function( n, p, k )
    local   q;

    q := p^k;
    if n = 1 and q = 4 then
        return GcdInt( 2, q-1 ) * 2;
    elif n = 1 and q = 9  then
        return GcdInt( 2, q-1 ) * 3;
    elif n = 2 and q = 2  then
        return GcdInt( n+1, q-1 ) * 2;
    elif n = 2 and q = 4  then
        return GcdInt( n+1, q-1 ) * 4^2;
    elif n = 3 and q = 2  then
        return GcdInt( n+1, q-1 ) * 2;
    else
        return GcdInt( n+1, q-1 );
    fi;

end;    

ChevA.uexponent := function( n, p, k )
    local   e,  i,  q,  qi;

    q  := p ^ k;
    e  := p ^ ( 1 + LogInt( n, p ) );
    qi := 1;
    for i  in [ 1 .. n+1 ]  do
        qi := qi * q;
        e  := LcmInt( e, qi - 1 );
    od;
    return e / GcdInt( n+1, q+1 );

end;
        

#############################################################################
##
#V  ChevB . . . . . . . . . . . . . . . . . . . . . . . O_2n+1(q) information
##
ChevB := rec( name := "ChevB", operations := ChevOps );

ChevB.order := function( n, p, k )
    local   i,  o,  q;

    if n < 3 and 2 < p  then
        return false;
    fi;
    q := p ^ k;
    o := q ^ ( n^2 );
    for i  in [ 1 .. n ]  do
        o := o * ( q^(2*i) - 1 );
    od;
    return o / 2;

end;

ChevB.possibilities := function( n, m, ord, p )
    local   l1,  l2,  t,  o,  d,  exc;

    # at first no possibilities at all
    l1 := [];
    l2 := [];

    # at first characteristic <p>
    if 2 < p  then
        t := [3,p,1];
        o := ChevB.order( t[1], t[2], t[3] );
        while o <= m  do
            while o <= m and ord mod o = 0  do
                Add( l1, Copy(t) );
                t[3] := t[3] + 1;
                o := ChevB.order( t[1], t[2], t[3] );
            od;
            t[1] := t[1] + 1;
            t[3] := 1;
            o := ChevB.order( t[1], t[2], t[3] );
        od;
    fi;

    # exceptions in lower bound formula or order
    exc := [ [3,3,1], [3,5,1] ];

    # now different characteristic
    t := [3,3,1];
    if p = 3  then
        t[2] := 5;
    fi;
    o := ChevB.order( t[1], t[2], t[3] );
    d := ChevB.lowerBound( t[1], t[2], t[3] );
    while ( t in exc ) or ( o <= m and d <= n )  do
        while ( t in exc ) or ( o <= m and d <= n )  do
            while (t in exc) or (o <= m and d <= n and ord mod o = 0)  do
                if o<=m and d<=n and (not t in exc or ord mod o=0)  then
                    Add( l2, Copy(t) );
                fi;
                t[1] := t[1] + 1;
                o := ChevB.order( t[1], t[2], t[3] );
                d := ChevB.lowerBound( t[1], t[2], t[3] );
            od;
            t[1] := 3;
            t[3] := t[3] + 1;
            o := ChevB.order( t[1], t[2], t[3] );
            d := ChevB.lowerBound( t[1], t[2], t[3] );
        od;
        t[1] := 3;
        t[2] := NextPrimeInt(t[2]);
        t[3] := 1;
        if t[2] = p  then
            t[2] := NextPrimeInt(t[2]);
        fi;
        o := ChevB.order( t[1], t[2], t[3] );
        d := ChevB.lowerBound( t[1], t[2], t[3] );
    od;

    # return the possibilities
    return [ l1, l2 ];

end;

ChevB.lowerBound := function( n, p, k )
    local   q;

    q := p ^ k;
    if n = 3 and q = 3  then
        return 27;
    elif q < 7  then
        return q^(n-1) * ( q^(n-1) - 1 );
    else
        return q^( 2*(n-1) ) - 1;
    fi;

end;

ChevB.automorphism := function( n, p, k )
    return 2 * k;
end;

ChevB.multiplier := function( n, p, k )
    local   q;

    q := p ^ k;
    if n = 3 and q = 3  then
        return 6;
    else
        return 2;
    fi;

end;    

ChevB.uexponent := function( n, p, k )
    local   e,  i,  q;

    #N this is too big, there should be a lower bound for the <p>-part
    e := p ^ ( 1 + LogInt( 2*n, p ) );
    q := p ^ k;
    for i  in [ 1 .. n ]  do
        e := LcmInt( e, q^i - 1 );
        e := LcmInt( e, q^i + 1 );
    od;
    return GcdInt( e, ChevB.order( n, p, k ) );

end;
        

#############################################################################
##
#V  ChevC . . . . . . . . . . . . . . . . . . . . . . .  Sp_2n(q) information
##
##  Exceptions:
##    Sp_4(2) = S_6
##
ChevC := rec( name := "ChevC", operations := ChevOps );

ChevC.order := function( n, p, k )
    local   i,  o,  q;

    if n < 2 or ( n = 2 and p = 2 and k = 1 )  then
        return false;
    fi;
    q := p ^ k;
    o := q ^ ( n^2 );
    for i  in [ 1 .. n ]  do
        o := o * ( q^(2*i) - 1 );
    od;
    return o / GcdInt( 2, q-1 );

end;

ChevC.possibilities := function( n, m, ord, p )
    local   l1,  l2,  t,  o,  d,  exc;

    # at first no possibilities at all
    l1 := [];
    l2 := [];

    # at first characteristic <p>
    t := [2,p,1];
    if p = 2  then
        t[3] := 2;
    fi;
    o := ChevC.order( t[1], t[2], t[3] );
    while o <= m  do
        while o <= m and ord mod o = 0  do
            Add( l1, Copy(t) );
            t[3] := t[3] + 1;
            o := ChevC.order( t[1], t[2], t[3] );
        od;
        t[1] := t[1] + 1;
        t[3] := 1;
        o := ChevC.order( t[1], t[2], t[3] );
    od;

    # exceptions in lower bound formula or order
    exc := [ [2,2,1], [3,2,1] ];

    # now different characteristic
    t := [2,2,1];
    if p = 2  then
        t[2] := 3;
    fi;
    if t[2] = 2  then
        t[1] := 2;
    fi;
    o := ChevC.order( t[1], t[2], t[3] );
    d := ChevC.lowerBound( t[1], t[2], t[3] );
    while ( t in exc ) or ( o <= m and d <= n )  do
        while ( t in exc ) or ( o <= m and d <= n )  do
            while (t in exc) or (o <= m and d <= n and ord mod o = 0)  do
                if o<=m and d<=n and (not t in exc or ord mod o=0)  then
                    Add( l2, Copy(t) );
                fi;
                t[1] := t[1] + 1;
                o := ChevC.order( t[1], t[2], t[3] );
                d := ChevC.lowerBound( t[1], t[2], t[3] );
            od;
            t[1] := 2;
            t[3] := t[3] + 1;
            o := ChevC.order( t[1], t[2], t[3] );
            d := ChevC.lowerBound( t[1], t[2], t[3] );
        od;
        t[1] := 2;
        t[2] := NextPrimeInt(t[2]);
        t[3] := 1;
        if t[2] = p  then
            t[2] := NextPrimeInt(t[2]);
        fi;
        o := ChevC.order( t[1], t[2], t[3] );
        d := ChevC.lowerBound( t[1], t[2], t[3] );
    od;

    # return the possibilities
    return [ l1, l2 ];

end;

ChevC.lowerBound := function( n, p, k )
    local   q;

    q := p ^ k;
    if n = 3 and q = 2  then
        return 7;
    elif p = 2  then
        return q^(n-1) * ( q^(n-1) - 1 ) * ( q - 1 ) / 2;
    else
        return ( q^n - 1 ) / 2;
    fi;

end;

ChevC.automorphism := function( n, p, k )
    return GcdInt( 2, p^k-1 ) * k;
end;

ChevC.multiplier := function( n, p, k )

    if n = 3 and p = 2 and k = 1  then
        return 2;
    else
        return GcdInt( 2, p^k-1 );
    fi;

end;    

ChevC.uexponent := function( n, p, k )
    local   e,  i,  q,  qi,  t;

    #N this is too big, there should be a lower bound for the <p>-part
    e  := p ^ ( 1 + LogInt( 2*n-1, p ) );
    q  := p ^ k;
    qi := 1;
    for i  in [ 1 .. n ]  do
        qi := qi * q;
        t  := LcmInt( qi-1, qi+1 );
        e  := LcmInt( e, t );
    od;
    return GcdInt( e, ChevC.order( n, p, k ) );

end;
        

#############################################################################
##
#V  ChevD . . . . . . . . . . . . . . . . . . . . . . .  O_2n+(q) information
##
ChevD := rec( name := "ChevD", operations := ChevOps );

ChevD.order := function( n, p, k )
    local   i,  o,  q;

    if n < 4  then
        return false;
    fi;
    q := p ^ k;
    o := q ^ (n^2-n) * ( q^n - 1 );
    for i  in [ 1 .. n-1 ]  do
        o := o * ( q^(2*i) - 1 );
    od;
    return o / GcdInt( 4, q^n-1 );

end;

ChevD.possibilities := function( n, m, ord, p )
    local   l1,  l2,  t,  o,  d,  exc;

    # at first no possibilities at all
    l1 := [];
    l2 := [];

    # at first characteristic <p>
    t := [4,p,1];
    o := ChevD.order( t[1], t[2], t[3] );
    while o <= m  do
        while o <= m  do
            if ord mod o = 0  then
                Add( l1, Copy(t) );
            fi;
            t[3] := t[3] + 1;
            o := ChevD.order( t[1], t[2], t[3] );
        od;
        t[1] := t[1] + 1;
        t[3] := 1;
        o := ChevD.order( t[1], t[2], t[3] );
    od;

    # exceptions in lower bound formula or order
    exc := [ [4,2,1], [4,3,1], [4,5,1] ];

    # now different characteristic
    t := [4,2,1];
    if p = 2  then
        t[2] := 3;
    fi;
    o := ChevD.order( t[1], t[2], t[3] );
    d := ChevD.lowerBound( t[1], t[2], t[3] );
    while ( o <= m and d <= n ) or ( t in exc )  do
        while ( o <= m and d <= n ) or ( t in exc )  do
            while ( o <= m and d <= n ) or ( t in exc )  do
                if o <= m and d <= n and ord mod o = 0  then
                    Add( l2, Copy(t) );
                fi;
                t[1] := t[1] + 1;
                o := ChevD.order( t[1], t[2], t[3] );
                d := ChevD.lowerBound( t[1], t[2], t[3] );
            od;
            t[1] := 4;
            t[3] := t[3] + 1;
            o := ChevD.order( t[1], t[2], t[3] );
            d := ChevD.lowerBound( t[1], t[2], t[3] );
        od;
        t[1] := 4;
        t[2] := NextPrimeInt(t[2]);
        t[3] := 1;
        if t[2] = p  then
            t[2] := NextPrimeInt(t[2]);
        fi;
        o := ChevD.order( t[1], t[2], t[3] );
        d := ChevD.lowerBound( t[1], t[2], t[3] );
    od;

    # return the possibilities
    return [ l1, l2 ];

end;

ChevD.lowerBound := function( n, p, k )
    local   q;

    q := p ^ k;
    if n = 4 and q = 2  then
        return 8;
    elif q < 7  then
        return q^(n-2) * ( q^(n-1) - 1 );
    else
        return ( q^(n-1) - 1 ) * ( q^(n-2) + 1 );
    fi;

end;

ChevD.automorphism := function( n, p, k )

    if n = 4  then
        return GcdInt( 2, p^k-1 )^2 * k * 6;
    elif n mod 2 = 0  then
        return GcdInt( 2, p^k-1 )^2 * k * 2;
    else
        return GcdInt( 4, p^(k*n)-1 ) * k * 2;
    fi;

end;

ChevD.multiplier := function( n, p, k )

    if n = 4 and p = 2 and k = 1  then
        return 4;
    elif n mod 2 = 0  then
        return GcdInt( 2, p^k-1 )^2;
    else
        return GcdInt( 4, p^(k*n)-1 );
    fi;

end;    

ChevD.uexponent := function( n, p, k )
    local   e,  i,  q;

    #N this is too big, there should be a lower bound for the <p>-part
    e := p ^ ( 1 + LogInt( 2*n-1, p ) );
    q := p ^ k;
    for i  in [ 1 .. n ]  do
        e := LcmInt( e, q^i - 1 );
        e := LcmInt( e, q^i + 1 );
    od;
    return GcdInt( e, ChevD.order( n, p, k ) );

end;
        

#############################################################################
##
#V  ChevG . . . . . . . . . . . . . . . . . . . . . . . .  G_2(q) information
##
ChevG := rec( name := "ChevG", operations := ChevOps );

ChevG.order := function( n, p, k )
    local   q;

    q := p ^ k;
    if n <> 2  then
        return false;
    fi;
    return q^6 * ( q^6-1 ) * ( q^2-1 );

end;

ChevG.possibilities := function( n, m, ord, p )
    local   l1,  l2,  t,  o,  d,  exc;

    # at first no possibilities at all
    l1 := [];
    l2 := [];

    # at first characteristic <p>
    t := [2,p,1];
    o := ChevG.order( t[1], t[2], t[3] );
    while o <= m and ord mod o = 0  do
        Add( l1, Copy(t) );
        t[3] := t[3] + 1;
        o := ChevG.order( t[1], t[2], t[3] );
    od;

    # exceptions in lower bound formula or order
    exc := [ [2,3,1], [2,2,1], [2,2,2] ];

    # now different characteristic
    t := [2,2,1];
    if p = 2  then
        t[2] := 3;
    fi;
    o := ChevG.order( t[1], t[2], t[3] );
    d := ChevG.lowerBound( t[1], t[2], t[3] );
    while ( o <= m and d <= n ) or ( t in exc )  do
        while ( o <= m and d <= n ) or ( t in exc )  do
            if o <= m and d <= n and ord mod o = 0  then
                Add( l2, Copy(t) );
            fi;
            t[3] := t[3] + 1;
            o := ChevG.order( t[1], t[2], t[3] );
            d := ChevG.lowerBound( t[1], t[2], t[3] );
        od;
        t[2] := NextPrimeInt(t[2]);
        t[3] := 1;
        if t[2] = p  then
            t[2] := NextPrimeInt(t[2]);
        fi;
        o := ChevG.order( t[1], t[2], t[3] );
        d := ChevG.lowerBound( t[1], t[2], t[3] );
    od;

    # return the possibilities
    return [ l1, l2 ];

end;

ChevG.lowerBound := function( n, p, k )
    local   q;

    q := p ^ k;
    if n = 2 and q = 3  then
        return 14;
    elif n =2 and q = 4  then
        return 12;
    else
        return q * ( q^2 - 1 );
    fi;

end;

ChevG.automorphism := function( n, p, k )

    if p = 3  then
        return k * 2;
    else
        return k;
    fi;

end;

ChevG.multiplier := function( n, p, k )

    if p = 3 and k = 1  then
        return 3;
    elif p = 2 and k = 2  then
        return 2;
    else
        return 1;
    fi;

end;    

ChevG.uexponent := ChevG.order; #N  What is the exponent of G_n(n)?
        

#############################################################################
##
#V  ChevF . . . . . . . . . . . . . . . . . . . . . . . .  F_4(q) information
##
ChevF := rec( name := "ChevF", operations := ChevOps );

ChevF.order := function( n, p, k )
    local   q;

    q := p ^ k;
    if n <> 4  then
        return false;
    fi;
    return q^24 * (q^12-1) * (q^8-1) * (q^6-1) * (q^2-1);

end;

ChevF.possibilities := function( n, m, ord, p )
    local   l1,  l2,  t,  o,  d,  exc;

    # at first no possibilities at all
    l1 := [];
    l2 := [];

    # at first characteristic <p>
    t := [4,p,1];
    o := ChevF.order( t[1], t[2], t[3] );
    while o <= m and ord mod o = 0  do
        Add( l1, Copy(t) );
        t[3] := t[3] + 1;
        o := ChevF.order( t[1], t[2], t[3] );
    od;

    # exceptions in lower bound formula or order
    exc := [ [4,3,1], [4,2,1] ];

    # now different characteristic
    t := [4,2,1];
    if p = 2  then
        t[2] := 3;
    fi;
    o := ChevF.order( t[1], t[2], t[3] );
    d := ChevF.lowerBound( t[1], t[2], t[3] );
    while ( o <= m and d <= n ) or ( t in exc )  do
        while ( o <= m and d <= n ) or ( t in exc )  do
            if o <= m and d <= n and ord mod o = 0  then
                Add( l2, Copy(t) );
            fi;
            t[3] := t[3] + 1;
            o := ChevF.order( t[1], t[2], t[3] );
            d := ChevF.lowerBound( t[1], t[2], t[3] );
        od;
        t[2] := NextPrimeInt(t[2]);
        t[3] := 1;
        if t[2] = p  then
            t[2] := NextPrimeInt(t[2]);
        fi;
        o := ChevF.order( t[1], t[2], t[3] );
        d := ChevF.lowerBound( t[1], t[2], t[3] );
    od;

    # return the possibilities
    return [ l1, l2 ];

end;

ChevF.lowerBound := function( n, p, k )
    local   q;

    q := p ^ k;
    if q = 2  then
        return 44;
    elif q mod 2 = 0  then
        return q^7 * ( q^3 - 1 ) * ( q - 1 ) / 2;

    else
        return q^6 * ( q^2 - 1 );
    fi;

end;

ChevF.automorphism := function( n, p, k )

    if p = 2  then
        return k * 2;
    else
        return k;
    fi;

end;

ChevF.multiplier := function( n, p, k )

    if p = 2 and k = 1 then
        return 2;
    else
        return 1;
    fi;

end;    

ChevF.uexponent := ChevF.order; #N What is the exponent of F_4(q)?


#############################################################################
##
#V  ChevE . . . . . . . . . . . . . . . . . . . . . . . .  E_n(q) information
##
ChevE := rec( name := "ChevE", operations := ChevOps );

ChevE.order := function( n, p, k )
    local   q;

    if 8 < n or n < 6  then
        return false;
    fi;
    q := p ^ k;
    if n = 6  then
        return q^36*(q^12-1)*(q^9-1)*(q^8-1)*(q^6-1)*(q^5-1)*(q^2-1)
               / GcdInt( 3, q-1 );
    elif n = 7  then
        return q^63*(q^18-1)*(q^14-1)*(q^12-1)*(q^10-1)*(q^8-1)
               *(q^6-1)*(q^2-1) / GcdInt( 2, q-1 );
    elif n = 8  then
        return q^120*(q^30-1)*(q^24-1)*(q^20-1)*(q^18-1)*(q^14-1)
               *(q^12-1)*(q^8-1)*(q^2-1);
    fi;

end;

ChevE.possibilities := function( n, m, ord, p )
    local   l1,  l2,  t,  o,  d,  exc,  i;

    # at first no possibilities at all
    l1 := [];
    l2 := [];

    # at first characteristic <p>
    t := [6,p,1];
    for i  in [ 6, 7, 8 ]  do
        t[1] := i;
        t[3] := 1;
        o := ChevE.order( t[1], t[2], t[3] );
        while o <= m  do
            if ord mod o = 0  then
                Add( l1, Copy(t) );
            fi;
            t[3] := t[3] + 1;
            o := ChevE.order( t[1], t[2], t[3] );
        od;
    od;

    # exceptions in lower bound formula or order
    exc := [];

    # now different characteristic
    t := [6,2,1];
    if p = 2  then
        t[2] := 3;
    fi;
    for i  in [ 6, 7, 8 ]  do
        t[1] := i;
        o := ChevE.order( t[1], t[2], t[3] );
        d := ChevE.lowerBound( t[1], t[2], t[3] );
        while ( o <= m and d <= n ) or ( t in exc )  do
            while ( o <= m and d <= n ) or ( t in exc )  do
                if o <= m and d <= n and ord mod o = 0  then
                    Add( l2, Copy(t) );
                fi;
                t[3] := t[3] + 1;
                o := ChevE.order( t[1], t[2], t[3] );
                d := ChevE.lowerBound( t[1], t[2], t[3] );
            od;
            t[2] := NextPrimeInt(t[2]);
            t[3] := 1;
            if t[2] = p  then
                t[2] := NextPrimeInt(t[2]);
            fi;
            o := ChevE.order( t[1], t[2], t[3] );
            d := ChevE.lowerBound( t[1], t[2], t[3] );
        od;
    od;

    # return the possibilities
    return [ l1, l2 ];

end;

ChevE.lowerBound := function( n, p, k )
    local   q;

    q := p ^ k;
    if n = 6  then
        return q^9 * ( q^2-1 );
    elif n = 7  then
        return q^15 * ( q^2-1 );
    elif n = 8  then
        return q^27 * ( q^2-1 );
    fi;

end;

ChevE.automorphism := function( n, p, k )

    if n = 6  then
        return GcdInt( 3, p^k-1 ) * k * 2;
    elif n = 7  then
        return GcdInt( 2, p^k-1 ) * k;
    elif n = 8  then
        return k;
    fi;

end;

ChevE.multiplier := function( n, p, k )

    if n = 6  then
        return GcdInt( 3, p^k-1 );
    elif n = 7  then
        return GcdInt( 2, p^k-1 );
    elif n = 8  then
        return 1;
    fi;

end;    

ChevE.uexponent := ChevE.order; #N What is the exponent of E_n(q)?


#############################################################################
##
#V  Chev2A  . . . . . . . . . . . . . . . . . . . . .  PSU(n+1,q) information
##
Chev2A := rec( name := "Chev2A", operations := ChevOps );

Chev2A.order := function( n, p, k )
    local   i,  o,  q;

    if n < 2 or ( n = 2 and p = 2 and k = 1 )  then
        return false;
    fi;
    q := p ^ k;
    o := q ^ ( n * ( n+1 ) / 2 );
    for i  in [ 1 .. n ]  do
        o := o * ( q^(i+1) - (-1)^(i+1) );
    od;
    return o;

end;

Chev2A.possibilities := function( n, m, ord, p )
    local   l1,  l2,  t,  o,  d,  exc;

    # at first no possibilities at all
    l1 := [];
    l2 := [];

    # at first characteristic <p>
    t := [2,p,1];
    if p = 2  then
        t[3] := 2;
    fi;
    o := Chev2A.order( t[1], t[2], t[3] );
    while o <= m  do
        while o <= m  do
            if ord mod o = 0  then
                Add( l1, Copy(t) );
            fi;
            t[3] := t[3] + 1;
            o := Chev2A.order( t[1], t[2], t[3] );
        od;
        t[1] := t[1] + 1;
        t[3] := 2;
        o := Chev2A.order( t[1], t[2], t[3] );
    od;

    # exceptions in lower bound formula or order
    exc := [ [2,2,1], [3,2,1], [3,3,1] ];

    # now different characteristic
    t := [2,2,1];
    if p = 2  then
        t[2] := 3;
    fi;
    if t[2] = 2  then
        t[1] := 2;
    fi;
    o := Chev2A.order( t[1], t[2], t[3] );
    d := Chev2A.lowerBound( t[1], t[2], t[3] );
    while ( o <= m and d <= n ) or ( t in exc )  do
        while ( o <= m and d <= n ) or ( t in exc )  do
            while ( o <= m and d <= n ) or ( t in exc )  do
                if o <= m and d <= n and ord mod o = 0  then
                    Add( l2, Copy(t) );
                fi;
                t[1] := t[1] + 1;
                o := Chev2A.order( t[1], t[2], t[3] );
                d := Chev2A.lowerBound( t[1], t[2], t[3] );
            od;
            t[1] := 2;
            t[3] := t[3] + 1;
            o := Chev2A.order( t[1], t[2], t[3] );
            d := Chev2A.lowerBound( t[1], t[2], t[3] );
        od;
        t[1] := 2;
        t[2] := NextPrimeInt(t[2]);
        t[3] := 1;
        if t[2] = p  then
            t[2] := NextPrimeInt(t[2]);
        fi;
        o := Chev2A.order( t[1], t[2], t[3] );
        d := Chev2A.lowerBound( t[1], t[2], t[3] );
    od;

    # return the possibilities
    return [ l1, l2 ];

end;

Chev2A.lowerBound := function( n, p, k )
    local   q;

    q := p ^ k;
    if n = 3 and q = 2  then
        return 4;
    elif n = 3 and q = 3  then
        return 6;
    elif n mod 2 = 0  then
        return ( q^(n+1) - 1 ) / ( q+1 );
    else
        return q * ( q^n - 1 ) / ( q+1 );
    fi;

end;

Chev2A.automorphism := function( n, p, k )
    return GcdInt( n+1, p^k+1 ) * k;
end;

Chev2A.multiplier := function( n, p, k )

    if n = 3 and p = 2 and k = 1  then
        return GcdInt( n+1, p^k+1 ) * 2;
    elif n = 3 and p = 3 and k = 1  then
        return GcdInt( n+1, p^k+1 ) * 9;
    elif n = 3 and p = 2 and k = 2  then
        return GcdInt( n+1, p^k+1 ) * 4;
    else
        return GcdInt( n+1, p^k+1 );
    fi;

end;    

Chev2A.uexponent := function( n, p, k )
    local   i,  e,  q,  qi;

    #N this is too big, there should be a lower bound for the <p>-part
    e  := p ^ ( 1 + LogInt( n, p ) );
    q  := p ^ k;
    qi := 1;
    for i  in [ 1 .. n+1 ]  do
        qi := qi * q;
        e  := LcmInt( e, qi - (-1)^i );
    od;
    return GcdInt( e, Chev2A.order( n, p, k ) );

end;

#############################################################################
##
#V  Chev2B  . . . . . . . . . . . . . . . . . . . . . . . . Sz(q) information
##
Chev2B := rec( name := "Chev2B", operations := ChevOps );

Chev2B.order := function( n, p, k )
    local   q;

    if n <> 2 or p <> 2 or k mod 2 = 0  then
        return false;
    fi;
    q := p ^ k;
    return q^2 * ( q^2+1 ) * ( q-1 );

end;

Chev2B.possibilities := function( n, m, ord, p )
    local   l1,  l2,  t,  o,  d,  exc;

    # at first no possibilities at all
    l1 := [];
    l2 := [];

    # at first characteristic <p>
    if p = 2  then
        t := [2,p,3];
        o := Chev2B.order( t[1], t[2], t[3] );
        while o <= m  do
            if ord mod o = 0  then
                Add( l1, Copy(t) );
            fi;
            t[3] := t[3] + 2;
            o := Chev2B.order( t[1], t[2], t[3] );
        od;
    fi;

    # exceptions in lower bound formula or order
    exc := [ [2,2,3] ];

    # now different characteristic
    if 2 < p  then
        t := [2,2,3];
        o := Chev2B.order( t[1], t[2], t[3] );
        d := Chev2B.lowerBound( t[1], t[2], t[3] );
        while ( o <= m and d <= n ) or ( t in exc )  do
            if o <= m and d <= n and ord mod o = 0  then
                Add( l2, Copy(t) );
            fi;
            t[3] := t[3] + 2;
            o := Chev2B.order( t[1], t[2], t[3] );
            d := Chev2B.lowerBound( t[1], t[2], t[3] );
        od;
    fi;

    # return the possibilities
    return [ l1, l2 ];

end;

Chev2B.lowerBound := function( n, p, k )
    local   q;

    q := p ^ k;
    if q = 8  then
        return 8;
    else
        return p ^ ((k-1)/2) * ( q - 1 );
    fi;

end;

Chev2B.automorphism := function( n, p, k )
    return k;
end;

Chev2B.multiplier := function( n, p, k )

    if p = 2 and k = 2  then
        return 4;
    else
        return 1;
    fi;

end;    

Chev2B.uexponent := Chev2B.order; #N What is the exponent of Sz(q)?


#############################################################################
##
#V  Chev2D  . . . . . . . . . . . . . . . . . . . . . .  O_2n-(q) information
##
Chev2D := rec( name := "Chev2D", operations := ChevOps );

Chev2D.order := function( n, p, k )
    local   o,  i,  q;

    if n < 4  then
        return false;
    fi;
    q := p ^ k;
    o := q ^ ( n*(n-1) ) * ( q^n+1 );
    for i  in [ 1 .. n-1 ]  do
        o := o * ( q^(2*i)-1 );
    od;
    return o / GcdInt( 4, q^n+1 );

end;

Chev2D.possibilities := function( n, m, ord, p )
    local   l1,  l2,  t,  o,  d,  exc;

    # at first no possibilities at all
    l1 := [];
    l2 := [];

    # at first characteristic <p>
    t := [4,p,1];
    o := Chev2D.order( t[1], t[2], t[3] );
    while o <= m  do
        while o <= m  do
            if ord mod o = 0  then
                Add( l1, Copy(t) );
            fi;
            t[3] := t[3] + 1;
            o := Chev2D.order( t[1], t[2], t[3] );
        od;
        t[1] := t[1] + 1;
        t[3] := 1;
        o := Chev2D.order( t[1], t[2], t[3] );
    od;

    # exceptions in lower bound formula or order
    exc := [];

    # now different characteristic
    t := [4,2,1];
    if p = 2  then
        t[2] := 3;
    fi;
    o := Chev2D.order( t[1], t[2], t[3] );
    d := Chev2D.lowerBound( t[1], t[2], t[3] );
    while ( o <= m and d <= n ) or ( t in exc )  do
        while ( o <= m and d <= n ) or ( t in exc )  do
            while ( o <= m and d <= n ) or ( t in exc )  do
                if o <= m and d <= n and ord mod o = 0  then
                    Add( l2, Copy(t) );
                fi;
                t[1] := t[1] + 1;
                o := Chev2D.order( t[1], t[2], t[3] );
                d := Chev2D.lowerBound( t[1], t[2], t[3] );
            od;
            t[1] := 4;
            t[3] := t[3] + 1;
            o := Chev2D.order( t[1], t[2], t[3] );
            d := Chev2D.lowerBound( t[1], t[2], t[3] );
        od;
        t[1] := 4;
        t[2] := NextPrimeInt(t[2]);
        t[3] := 1;
        if t[2] = p  then
            t[2] := NextPrimeInt(t[2]);
        fi;
        o := Chev2D.order( t[1], t[2], t[3] );
        d := Chev2D.lowerBound( t[1], t[2], t[3] );
    od;

    # return the possibilities
    return [ l1, l2 ];

end;

Chev2D.lowerBound := function( n, p, k )
    local   q;

    q := p ^ k;
    return ( q^(n-1) + 1 ) * ( q^(n-2) - 1 );

end;

Chev2D.automorphism := function( n, p, k )
    return GcdInt( 4, p^(k*n)+1 ) * 2 * k;
end;

Chev2D.multiplier := function( n, p, k )
    return GcdInt( 4, p^(k*n)+1 );
end;    

Chev2D.uexponent := function( n, p, k )
    local   e,  i,  q;

    #N this is too big, there should be a lower bound for the <p>-part
    e := p ^ ( 1 + LogInt( 2*n-1, p ) );
    q := p ^ k;
    for i  in [ 1 .. n ]  do
        e := LcmInt( e, q^i - 1 );
        e := LcmInt( e, q^i + 1 );
    od;
    return GcdInt( e, Chev2D.order( n, p, k ) );

end;


#############################################################################
##
#V  Chev3D  . . . . . . . . . . . . . . . . . . . . . . . 3D_4(q) information
##
Chev3D := rec( name := "Chev3D", operations := ChevOps );

Chev3D.order := function( n, p, k )
    local   q;

    if n <> 4  then
        return false;
    fi;
    q := p ^ k;
    return q^12 * ( q^8+q^4+1 ) * ( q^6-1 ) * ( q^2-1 );

end;

Chev3D.possibilities := function( n, m, ord, p )
    local   l1,  l2,  t,  o,  d,  exc;

    # at first no possibilities at all
    l1 := [];
    l2 := [];

    # at first characteristic <p>
    t := [4,p,1];
    o := Chev3D.order( t[1], t[2], t[3] );
    while o <= m  do
        if ord mod o = 0  then
            Add( l1, Copy(t) );
        fi;
        t[3] := t[3] + 1;
        o := Chev3D.order( t[1], t[2], t[3] );
    od;

    # exceptions in lower bound formula or order
    exc := [];

    # now different characteristic
    t := [4,2,1];
    if p = 2  then
        t[2] := 3;
    fi;
    o := Chev3D.order( t[1], t[2], t[3] );
    d := Chev3D.lowerBound( t[1], t[2], t[3] );
    while ( o <= m and d <= n ) or ( t in exc )  do
        while ( o <= m and d <= n ) or ( t in exc )  do
            if o <= m and d <= n and ord mod o = 0  then
                Add( l2, Copy(t) );
            fi;
            t[3] := t[3] + 1;
            o := Chev3D.order( t[1], t[2], t[3] );
            d := Chev3D.lowerBound( t[1], t[2], t[3] );
        od;
        t[2] := NextPrimeInt(t[2]);
        t[3] := 1;
        if t[2] = p  then
            t[2] := NextPrimeInt(t[2]);
        fi;
        o := Chev3D.order( t[1], t[2], t[3] );
        d := Chev3D.lowerBound( t[1], t[2], t[3] );
    od;

    # return the possibilities
    return [ l1, l2 ];

end;

Chev3D.lowerBound := function( n, p, k )
    local   q;

    q := p ^ k;
    return q^3 * ( q^2-1 );

end;

Chev3D.automorphism := function( n, p, k )
    return 3 * k;
end;

Chev3D.multiplier := function( n, p, k )
    return 1;
end;    

Chev3D.uexponent := Chev3D.order; #N What is the exponent of 3D(q)?


#############################################################################
##
#V  Chev2G  . . . . . . . . . . . . . . . . . . . . . . . 2G_2(q) information
##
Chev2G := rec( name := "Chev2G", operations := ChevOps );

Chev2G.order := function( n, p, k )
    local   q;

    if n <> 2 or p <> 3 or k mod 2 = 0  then
        return false;
    fi;
    q := p ^ k;
    return q^3 * ( q^3+1 ) * ( q-1 );

end;

Chev2G.possibilities := function( n, m, ord, p )
    local   l1,  l2,  t,  o,  d,  exc;

    # at first no possibilities at all
    l1 := [];
    l2 := [];

    # at first characteristic <p>
    if p = 3  then
        t := [2,p,3];
        o := Chev2G.order( t[1], t[2], t[3] );
        while o <= m  do
            if ord mod o = 0  then
                Add( l1, Copy(t) );
            fi;
            t[3] := t[3] + 2;
            o := Chev2G.order( t[1], t[2], t[3] );
        od;
    fi;

    # exceptions in lower bound formula or order
    exc := [];

    # now different characteristic
    if 3 <> p  then
        t := [2,3,3];
        o := Chev2G.order( t[1], t[2], t[3] );
        d := Chev2G.lowerBound( t[1], t[2], t[3] );
        while ( o <= m and d <= n ) or ( t in exc )  do
            if o <= m and d <= n and ord mod o = 0  then
                Add( l2, Copy(t) );
            fi;
            t[3] := t[3] + 2;
            o := Chev2G.order( t[1], t[2], t[3] );
            d := Chev2G.lowerBound( t[1], t[2], t[3] );
        od;
    fi;

    # return the possibilities
    return [ l1, l2 ];

end;

Chev2G.lowerBound := function( n, p, k )
    local   q;

    q := p ^ k;
    return q * ( q - 1 );

end;

Chev2G.automorphism := function( n, p, k )
    return k;
end;

Chev2G.multiplier := function( n, p, k )
    return 1;

end;    

Chev2G.uexponent := Chev2G.order; #N What is the exponent of 2G_2(q)?


#############################################################################
##
#V  Chev2F  . . . . . . . . . . . . . . . . . . . . . . . 2F_4(q) information
##
Chev2F := rec( name := "Chev2F", operations := ChevOps );

Chev2F.order := function( n, p, k )
    local   q,  o;

    if n <> 4 or p <> 2 or k mod 2 = 0  then
        return false;
    fi;
    q := p ^ k;
    if q = 2  then
        return q^12 * ( q^6+1 ) * ( q^4-1 ) * ( q^3+1 ) * ( q-1 ) / 2;
    else
        return q^12 * ( q^6+1 ) * ( q^4-1 ) * ( q^3+1 ) * ( q-1 );
    fi;

end;

Chev2F.possibilities := function( n, m, ord, p )
    local   l1,  l2,  t,  o,  d,  exc;

    # at first no possibilities at all
    l1 := [];
    l2 := [];

    # at first characteristic <p>
    if p = 2  then
        t := [4,p,1];
        o := Chev2F.order( t[1], t[2], t[3] );
        while o <= m  do
            if ord mod o = 0  then
                Add( l1, Copy(t) );
            fi;
            t[3] := t[3] + 2;
            o := Chev2F.order( t[1], t[2], t[3] );
        od;
    fi;

    # exceptions in lower bound formula or order
    exc := [ [4,2,1] ];

    # now different characteristic
    if 2 < p  then
        t := [4,2,1];
        o := Chev2F.order( t[1], t[2], t[3] );
        d := Chev2F.lowerBound( t[1], t[2], t[3] );
        while ( o <= m and d <= n ) or ( t in exc )  do
            if o <= m and d <= n and ord mod o = 0  then
                Add( l2, Copy(t) );
            fi;
            t[3] := t[3] + 2;
            o := Chev2F.order( t[1], t[2], t[3] );
            d := Chev2F.lowerBound( t[1], t[2], t[3] );
        od;
    fi;

    # return the possibilities
    return [ l1, l2 ];

end;

Chev2F.lowerBound := function( n, p, k )
    local   q;

    q := p ^ k;
    return p ^ ( (k-1)/2 ) * q^4 * ( q-1 );

end;

Chev2F.automorphism := function( n, p, k )

    if k = 1  then
        return 2 * k;
    else
        return k;
    fi;

end;

Chev2F.multiplier := function( n, p, k )
    return 1;
end;    

Chev2F.uexponent := Chev2F.order; #N What is the exponent?


#############################################################################
##
#V  Chev2E  . . . . . . . . . . . . . . . . . . . . . . . 2E_6(q) information
##
Chev2E := rec( name := "Chev2E", operations := ChevOps );

Chev2E.order := function( n, p, k )
    local   q;

    if n <> 6  then
        return false;
    fi;
    q := p ^ k;
    return q^36*(q^12-1)*(q^9+1)*(q^8-1)*(q^6-1)*(q^5+1)*(q^2-1)/GcdInt(3,q+1);

end;

Chev2E.possibilities := function( n, m, ord, p )
    local   l1,  l2,  t,  o,  d,  exc;

    # at first no possibilities at all
    l1 := [];
    l2 := [];

    # at first characteristic <p>
    t := [6,p,1];
    o := Chev2E.order( t[1], t[2], t[3] );
    while o <= m  do
        if ord mod o = 0  then
            Add( l1, Copy(t) );
        fi;
        t[3] := t[3] + 1;
        o := Chev2E.order( t[1], t[2], t[3] );
    od;

    # exceptions in lower bound formula or order
    exc := [];

    # now different characteristic
    t := [6,2,1];
    if p = 2  then
        t[2] := 3;
    fi;
    o := Chev2E.order( t[1], t[2], t[3] );
    d := Chev2E.lowerBound( t[1], t[2], t[3] );
    while ( o <= m and d <= n ) or ( t in exc )  do
        while ( o <= m and d <= n ) or ( t in exc )  do
            if o <= m and d <= n and ord mod o = 0  then
                Add( l2, Copy(t) );
            fi;
            t[3] := t[3] + 1;
            o := Chev2E.order( t[1], t[2], t[3] );
            d := Chev2E.lowerBound( t[1], t[2], t[3] );
        od;
        t[2] := NextPrimeInt(t[2]);
        t[3] := 1;
        if t[2] = p  then
            t[2] := NextPrimeInt(t[2]);
        fi;
        o := Chev2E.order( t[1], t[2], t[3] );
        d := Chev2E.lowerBound( t[1], t[2], t[3] );
    od;

    # return the possibilities
    return [ l1, l2 ];

end;

Chev2E.lowerBound := function( n, p, k )
    local   q;

    q := p ^ k;
    return q^9 * ( q^2-1 );

end;

Chev2E.automorphism := function( n, p, k )
    local   q;

    q := p ^ k;
    return GcdInt( 3, q+1 ) * 2 * k;

end;

Chev2E.multiplier := function( n, p, k )

    if p = 2 and k = 1  then
        return GcdInt( 3, p^k+1 ) * 4;
    else
        return GcdInt( 3, p^k+1 );
    fi;

end;    

Chev2E.uexponent := Chev2E.order; #N What is the exponent?
