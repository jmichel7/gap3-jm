#############################################################################
##
#A  Matrix package                                      Derek Holt
#A                                                      Charles Leedham-Green
#A                                                      Eamonn O'Brien
#A                                                      Sarah Rees 
##
#A  @(#)$Id: polyfdeg.g,v 1.1 1997/03/10 13:52:36 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log: polyfdeg.g,v $
#H  Revision 1.1  1997/03/10 13:52:36  gap
#H  VERSION 1.0
#H
#H  Revision 1.2  1997/01/05 10:49:32  fceller
#H  added Eamonn's new version to the reprository
#H
#H  Revision 1.1  1996/12/25 09:04:25  fceller
#H  changed long filenames to MS-DOS conform filenames,
#H  the init files are *NOT* yet updated
#H
#H  Revision 1.1  1996/11/28 13:14:52  fceller
#H  added "smash" and "reducible" to the repository
#H
##
#############################################################################
##
#F  FactorsSquarefreePolDeg( <R>, <f>, <deg> )  ... fixed degree factors
##                                                 of a squarefree polynom
##
##   'FactorsSquarefreePolDeg' returns a list of the normed irreducible factors
##   of degree deg  of  the  polynom <f>, which must be squarefree,
##   in the polynomial ring <R> over finite field of order  q.
##   f must have leading coefficient 1. It is based on
##   FiniteFieldPolynomialRingOps.FactorsSquarefree in polyfin.g
##   and it calls  FiniteFieldPolynomialRingOps.FactorsCommonDegree
##   IT CAN ONLY BE USED IF f HAS NO IRREDUCIBLE FACTORS OF DEGREE < DEG.
## 
FactorsSquarefreePolDeg := function ( R, f, deg )
    local       facs,           ## list of factors (result)
                cyc,            ## cyclotomic polynomial $x^{q^deg} - x$
                gcd,            ## gcd of $f$ and $cyc$
		 px,
                  d;

    # if <f> has a trivial constant term signal an error
    if f.valuation <> 0  then
        Error( "<f> must have a non-trivial constant term" );
    fi;

    # <facs> will contain factorisation
    facs := [];

    ## handle trivial case
    if Degree( f ) = 1  then
        return [ f ];
    fi;

    ## in the following $cyc = x^{q^deg} - x$
    px :=  Polynomial( R.baseRing, [ R.baseRing.zero, R.baseRing.one ] );
    cyc := Copy(px);

    for d in [1..deg] do
       cyc := PowerMod( R, cyc, R.baseRing.size, f );
    od;
    cyc := cyc - px;

    ## compute the gcd of $f$ and $x^{q^d} - x$
    gcd := Gcd(R, f, cyc );

    ## split the gcd with 'R.operations.FactorsCommonDegree'
    if 0 < Degree( gcd )  then
            Append( facs, R.operations.FactorsCommonDegree( R, gcd, deg ) );
    fi;

    ## return the factorization
    return facs;
end;

#############################################################################
##
#F  FactorsPolDeg(<R>, <f>, <deg> )  . .  fixed degree factors of a polynom
##
##  'FactorsPolDeg' is based on FiniteFieldPolynomialRingOps.Factors.
##  It returns a list of the
##  normed irreducible factors of degree deg  of the polynomial f.
##  It first reduces to the square-free case and then calls
##  FactorsSquarefreePolDeg.
##  IT CAN ONLY BE USED IF f HAS NO IRREDUCIBLE FACTORS OF DEGREE < DEG.
##  <f> must be a polynomial
##   in the polynomial ring <R> over finite field of order  q.
##
FactorsPolDeg := function ( R, f, deg )
    local  facs, d, g, h, r, i, v, k, l;

    ## handle trivial cases
    if Degree(f) < 2  then
        return [ f ];
    elif Length(f.coefficients) = 1  then
        l := List( [ 1 .. f.valuation ], x -> Indeterminate(f.baseRing) );
        l[1] := l[1] * f.coefficients[1];
        return l;
    fi;

    # make the polynomial normed
    g := R.operations.StandardAssociate( R, f );
    v := g.valuation;
    k := Polynomial( R.baseRing, g.coefficients );

    ## compute the deriviative
    d := R.operations.Derivative( R, k );

    ## if the derivative is nonzero then $k / GcdPol(k,d)$ is squarefree
    if d <>  R.zero  then

        ## compute the gcd of <k> and the derivative <d>
        g := Gcd( R, k, d );

        ## find factors of the squarefree quotient and the remainder
        facs := FactorsSquarefreePolDeg(R, Quotient(R,k,g), deg );
        for h in ShallowCopy(facs) do
            while 0 = Length( R.operations.EuclideanRemainder( R, g, h ).coefficients )  do
                 Add(facs,h);
                 g := Quotient(R,g,h);
             od;
        od;
        if 0 < Degree( g )  then
            Append( facs, FactorsPolDeg( R, g, deg ) );
        fi;

    ## otherwise <k> is the <p>-th power of another polynom <r>
    else

        ## compute the <p>-th root of <f>
        r := R.operations.RootsRepresentative( R, k, R.baseRing.char );

        ## factor this polynom
        h := FactorsPolDeg( R,r,deg );

        ## each factor appears <p> times in <f>
        facs := [];
        for i  in [1..R.baseRing.char]  do
            Append(facs,h);
        od;

    fi;

    # Sort the factorization
    Append( facs, List( [1..v], x -> Indeterminate(R.baseRing) ) );
    Sort(facs);

    return facs;
end;
