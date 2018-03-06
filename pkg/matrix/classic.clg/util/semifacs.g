#############################################################################
##
#A  semifacs.g                  GAP library                      Frank Celler
##
#H  @(#)$Id: semifacs.g,v 1.1 1997/03/10 13:49:44 gap Exp $
##
#Y  Copyright (C) 1995,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains functions which compute semi factors.
##
Revision_semifacs_g :=
    "@(#)$Id: semifacs.g,v 1.1 1997/03/10 13:49:44 gap Exp $";


#############################################################################
##
#F  SemiFactorsInt( <n>, <q>, <d> ) . . . . . . . . . . . semi factors of <n>
##
SemiFactorsInt := function ( n, q, d )
    local   sign,  factors,  composite,  p,  i,  tmp,  new,  m,  qi;

    # make <n> positive and handle trivial cases
    sign := 1;
    if n < 0  then sign := -sign;  n := -n;  fi;
    if n < 4  then return [ [ sign * n ], [] ];  fi;

    # <factors> holds the prime factors of <n>,  <composite> the semi ones
    factors   := [];
    composite := [];

    # do trial divisions by the primes less than 1000
    for p  in Primes  do
        while n mod p = 0  do Add( factors, p );  n := n / p;  od;
        if n < (p+1)^2 and 1 < n  then Add(factors,n);  n := 1;  fi;
        if n = 1  then
            factors[1] := sign*factors[1];
            return [ factors, composite ];
        fi;
    od;

    # do trial divisions by known factors
    for p  in Primes2  do
        while n mod p = 0  do Add( factors, p );  n := n / p;  od;
        if n = 1  then
            factors[1] := sign*factors[1];
            return [ factors, composite ];
        fi;
    od;

    # handle perfect powers
    p := SmallestRootInt( n );
    if p < n  then
        while 1 < n  do
            Append( factors, FactorsInt(p) );
            n := n / p;
        od;
        Sort( factors );
        factors[1] := sign * factors[1];
        return [ factors, composite ];
    fi;

    # factor gcds with <n> and <q>^<i>-1 using 'FactorsRho'
    qi := 1;
    for i  in [ 1 .. d ]  do
        if 1 < n  then
            qi := qi * q;
            tmp := GcdInt( n, qi-1 );
            if 1 < tmp  then
                tmp := FactorsRho( tmp, 1, 16, 1024 );
                Append( factors, tmp[1] );
                Append( composite, tmp[2] );
                n := n / Product(tmp[1]);
                n := n / Product(tmp[2]);
            fi;
        fi;
    od;

    # factor remaining <n>
    if 1 < n  then
        tmp := FactorsRho( n, 1, 16, 1024 );
        Append( factors, tmp[1] );
        Append( composite, tmp[2] );
    fi;

    # make sure there is no weird common non-prime factor
    Sort(composite);
    if 1 < Length(composite)  then
        repeat
            new := [];
            for n  in composite  do
                for m  in Difference( composite, [n] )  do
                    if 1 < n  then
                        tmp := GcdInt( n, m );
                        if 1 < tmp  then
                            Add( new, tmp );
                            n := n / tmp;
                        fi;
                    fi;
                od;
                if 1 < n  then
                    Add( new, n );
                fi;
            od;
            Sort(new);
            tmp := composite;
            composite := new;
        until composite = tmp;
    fi;

    # sort factors
    Sort(factors);

    # add sign
    if 0 < Length(factors)  then
        factors[1] := sign * factors[1];
    else
        composite[1] := sign * composite[1];
    fi;

    # and return
    return [ factors, composite ];

end;


#############################################################################
##
#F  SemiPrimePowersInt( <n>, <q>, <d> ) . . . . . .  semi prime powers of <n>
##
SemiPrimePowersInt := function( n, q, d )
    local   res,  lst,  pows,  p;

    if n = 1  then
        return [];
    elif n = 0  then
        Error( "<n> must be non zero" );
    elif n < 0  then
        n := -1 * n;
    fi;
    res := [];
    for lst  in SemiFactorsInt( n, q, d )  do
        pows := [];
        for p  in Set(lst)  do
            Add( pows, p );
            Add( pows, Number( lst, x -> x = p ) );
        od;
        Add( res, pows );
    od;
    return res;

end;
