#############################################################################
##
#W  ppd.g                       GAP library                      Frank Celler
##
#H  @(#)$Id: ppd.g,v 1.1 1997/03/10 13:49:41 gap Exp $
##
#Y  Copyright (C) 1995,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains functions to deal with PPD elements.
##
RevisionMatrix.classic_util_ppd_g :=
    "@(#)$Id: ppd.g,v 1.1 1997/03/10 13:49:41 gap Exp $";


#############################################################################
##

#F  PPDPartPDM1( <d>, <p> ) . . . . . . . . compute the ppd part in <p>^<d>-1
##
PPDPartPDM1B := function( d, p )
    local   n,  q,  i,  m,  x,  y;

    # compute the (repeated) gcd with p^d-1
    n := p^d - 1;
    x := 1;
    q := 1;
    for i  in [ 1 .. d-1 ]  do
        q := q * p;
        if d mod i = 0  then
            repeat
                m := GcdInt( n, q-1 );
                n := n / m;
                x := x * m;
            until m = 1;
        fi;
    od;

    # compute the possible gcd with <d>+1
    y := 1;
    if IsPrimeInt(d+1) and (n mod (d+1)) = 0 and (n mod (d+1)^2) <> 0  then
        y := d+1;
        n := n / (d+1);
    fi;

    # and return
    return rec( ppd := y,  lppd := n,  quo := x );

end;


#############################################################################
##
#F  PPDIrreducibleFactor( <R>, <f>, <d>, <q> )  . . . .  large factors of <f>
##
PPDIrreducibleFactor := function ( R, f, d, q )
    local   d,  px,  pow,  i,  cyc,  gcd,  a;

    # handle trivial case
    if Degree(f) <= 2  then
        return false;
    fi;

    # compute the deriviative
    a := Derivative( R, f );

    # if the derivative is nonzero then $f / Gcd(f,a)$ is squarefree
    if a <> R.zero  then

        # compute the gcd of <f> and the derivative <a>
        f := Quotient( R, f, Gcd( R, f, a ) );

        # $deg(f) <= d/2$ implies that there is no large factor
        if Degree(f) <= d/2  then
            return false;
        fi;

        # remove small irreducible factors
        px  := Indeterminate(R.baseRing);
        pow := PowerMod( R, px, q, f );
        for i  in [ 1 .. QuoInt(d,2) ]  do

            # next cyclotomic polynomial x^(q^i)-x
            cyc := pow - px;

            # compute the gcd of <f> and <cyc>
            gcd := Gcd( R, f, cyc );
            if 0 < Degree(gcd)  then
                f := Quotient( R, f, gcd );
                if Degree(f) <= d/2  then
                    return false;
                fi;
            fi;

            # replace <pow> by x^(q^(i+1))
            pow := PowerMod( R, pow, q, f );
        od;
        return StandardAssociate( R, f );

    # otherwise <f> is the <p>-th power of another polynomial <r>
    else
        return false;
    fi;

end;


#############################################################################
##
#F  IsPpdElement( <F>, <m>, <d>, <p>, <a> )
##
IsPpdElement := function( F, m, d, p, a )
    local   c,  R,  pm,  g;

    # compute the characteristic polynomial
    if IsMatrix(m)  then
        c := CharacteristicPolynomial( FiniteFieldMatrices, m );
        c.baseRing := F;
    else
        c := m;
    fi;

    # try to find a large factor
    R := PolynomialRing(F);
    c := PPDIrreducibleFactor( R, c, d, p^a );

    # return if we failed to find one
    if c = false  then
        return false;
    fi;

    # find the ppd and lppd parts
    pm := PPDPartPDM1B( Degree(c)*a, p );

    # get rid of the non-ppd part
    g := PowerMod( Indeterminate(F), pm.quo, c );

    # if it is one there is no ppd involved
    if g = R.one  then
        return false;
    fi;

    # check if there is a non-large ppd involved
    if 1 < pm.ppd  then
        g := PowerMod( g, pm.ppd, c );
        if g = R.one  then
            return [ Degree(c), false ];
        else
            return [ Degree(c), true ];
        fi;
    elif 1 < pm.lppd  then
        return [ Degree(c), true ];
    else
        Error( "should not happen" );
    fi;

end;


#############################################################################
##

#E  classic_util_ppd.g  . . . . . . . . . . . . . . . . . . . . . . ends here
##
