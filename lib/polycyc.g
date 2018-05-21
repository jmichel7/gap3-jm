#############################################################################
##
#A  polycyc.g     modified from GAP library                      Frank Celler
#A                                                         & Alexander Hulpke
##  
##
#Y  Copyright 1990-1995,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains functions for polynomials over the Cyclotomics.
##  Now faster internal functions are used for the arithmetic.
##  Thanks to the authors for sending the file. 

#############################################################################
##
#V  CyclotomicsPolynomialOps  . . . . . . . . polynomial over the cyclotomics
##
CyclotomicsPolynomialOps := OperationsRecord( "CyclotomicsPolynomialOps",
                                            PolynomialOps );

#############################################################################
##
#F  CyclotomicsPolynomialOps.\+ . . . . . . . . . . .  sum of two polynomials
##
CyclotomicsPolynomialOps.\+ := function( l, r ) local   sum,  val,  vdf;

    # handle the case that one argument is a list
    if IsList(l)  then return List( l, x -> x+r );
    elif IsList(r)  then return List( r, x -> l+x );
    fi;

    # handle the case <scalar> + <polynomial>
    if not (IsRec(l) and IsBound(l.isPolynomial) and l.isPolynomial)  then

        # <r> must have the cyclotomics as base ring
        if Cyclotomics <> r.baseRing  then
            Error( "<r> must have the cyclotomics as base ring" );
        fi;

        # <l> must lie in the base ring of <r>
        if not IsCyc(l)  then
            Error( "<l> must lie in the base ring of <r>" );
        fi;

        # if <l> is trivial return <r>
        if l = 0  then
            return r;
        fi;
 
        # otherwise convert <l> into a polynomial
        l := CyclotomicsOps.FastPolynomial( Cyclotomics, [l], 0 );
    fi;

    # handle the case <polynomial> + <scalar>
    if not (IsRec(r) and IsBound(r.isPolynomial) and r.isPolynomial)  then

        # <l> must have the cyclotomics as base ring
        if Cyclotomics <> l.baseRing  then
            Error( "<l> must have the cyclotomics as base ring" );
        fi;

        # <r> must lie in the base ring of <l>
        if not IsCyc(r)  then
            Error( "<r> must lie in the base ring of <l>" );
        fi;

        # if <r> is trivial return <l>
        if r = 0  then
            return l;
        fi;
 
        # otherwise convert <r> into a polynomial
        r := CyclotomicsOps.FastPolynomial( Cyclotomics, [r], 0 );
    fi;

    # depth greater than one are handle by our superclass
    if not IsBound(l.depth)  then  l.operations.Depth(l);  fi;
    if not IsBound(r.depth)  then  r.operations.Depth(r);  fi;
    if 1 <> l.depth or 1 <> r.depth  then
        return PolynomialOps.\+( l, r );

    # give up if we have rings other then the cyclotomics
    elif Cyclotomics <> l.baseRing  then
        Error( "<l> must have the cyclotomics as base ring" );
    elif Cyclotomics <> r.baseRing  then
        Error( "<r> must have the cyclotomics as base ring" );

    # if <l> is the null polynomial return <r>
    elif Length(l.coefficients) = 0  then
        return r;

    # if <r> is the null polynomial return <l>
    elif Length(r.coefficients) = 0  then
        return l;

    # sum of two polynomials
    else

        # get the valuation minimum;
        vdf := r.valuation - l.valuation;

        # if <r>.valuation is the minimum shift <l>
        if r.valuation < l.valuation  then
            val := r.valuation;
            sum := ShiftedCoeffs( l.coefficients, -vdf );
            AddCoeffs( sum, r.coefficients );

        # if <l>.valuation is the minimum shift <r>
        elif l.valuation < r.valuation  then
            val := l.valuation;
            sum := ShiftedCoeffs( r.coefficients, vdf );

            # the cyclotomics are commutative
            AddCoeffs( sum, l.coefficients );

        # otherwise they are equal
        else
            sum := SumCoeffs( l.coefficients, r.coefficients );
            val := l.valuation;
        fi;

        # return the sum
        sum := CyclotomicsOps.FastPolynomial( Cyclotomics, sum, val );
        sum.depth := 1;
        return sum;
    fi;

end;


#############################################################################
##
#F  CyclotomicsPolynomialOps.\- . . . . . . . . . . . . diff of two polynomials
##
CyclotomicsPolynomialOps.\- := function( l, r ) local   dif,  val,  vdf;

    # handle the case that one argument is a list
    if IsList(l)  then return List( l, x -> x-r );
    elif IsList(r)  then return List( r, x -> l-x );
    fi;

    # handle the case <scalar> - <polynomial>
    if not (IsRec(l) and IsBound(l.isPolynomial) and l.isPolynomial)  then

        # <r> must have the cyclotomics as base ring
        if Cyclotomics <> r.baseRing  then
            Error( "<r> must have the cyclotomics as base ring" );
        fi;

        # <l> must lie in the base ring of <r>
        if not IsCyc(l)  then
            Error( "<l> must lie in the base ring of <r>" );
        fi;

        # if <l> is trivial return -<r>
        if l = 0  then
            return CyclotomicsOps.FastPolynomial(
                       Cyclotomics, (-1) * r.coefficients, r.valuation );
        fi;
 
        # otherwise convert <l> into a polynomial
        l := CyclotomicsOps.FastPolynomial( Cyclotomics, [l], 0 );
    fi;

    # handle the case <polynomial> - <scalar>
    if not (IsRec(r) and IsBound(r.isPolynomial) and r.isPolynomial)  then

        # <l> must have the cyclotomics as base ring
        if Cyclotomics <> l.baseRing  then
            Error( "<l> must have the cyclotomics as base ring" );
        fi;

        # <r> must lie in the base ring of <l>
        if not IsCyc(r)  then
            Error( "<r> must lie in the base ring of <l>" );
        fi;

        # if <r> is trivial return <l>
        if r = 0  then
            return l;
        fi;
 
        # otherwise convert <r> into a polynomial
        r := CyclotomicsOps.FastPolynomial( Cyclotomics, [r], 0 );
    fi;

    # depth greater than one are handle by our superclass
    if not IsBound(l.depth)  then  l.operations.Depth(l);  fi;
    if not IsBound(r.depth)  then  r.operations.Depth(r);  fi;
    if 1 <> l.depth or 1 <> r.depth  then
        return PolynomialOps.\-( l, r );

    # give up if we have rings other then the cyclotomics
    elif Cyclotomics <> l.baseRing  then
        Error( "<l> must have the cyclotomics as base ring" );
    elif Cyclotomics <> r.baseRing  then
        Error( "<r> must have the cyclotomics as base ring" );

    # if <l> is the null polynomial return -<r>
    elif Length(l.coefficients) = 0  then
        return -r;

    # if <r> is the null polynomial return <l>
    elif Length(r.coefficients) = 0  then
        return l;

    # difference of two polynomials
    else

        # get the valuation minimum;
        vdf := r.valuation - l.valuation;

        # if <r>.valuation is the minimum shift <l>
        if r.valuation < l.valuation  then
            val := r.valuation;
            dif := ShiftedCoeffs( l.coefficients, -vdf );
            AddCoeffs( dif, r.coefficients, -1 );

        # if <l>.valuation is the minimum shift <r>
        elif l.valuation < r.valuation  then
            val := l.valuation;
            dif := (-1)*ShiftedCoeffs( r.coefficients, vdf );

            # the cyclotomics are commutative
            AddCoeffs( dif, l.coefficients );

        # otherwise they are equal
        else
            val := l.valuation;
            dif := Copy(l.coefficients);
            AddCoeffs( dif, r.coefficients, -1 );
        fi;

        # return the difference
        dif := CyclotomicsOps.FastPolynomial( Cyclotomics, dif, val );
        dif.depth := 1;
        return dif;
    fi;

end;


#############################################################################
##
#F  CyclotomicsPolynomialOps.\*  . . . . . . . . .  product of two polynomials
##
CyclotomicsPolynomialOps.\* := function( l, r ) local   R,  prd,  val;

    # handle the case that one argument is a list
    if IsList(l)  then return List( l, x -> x*r );
    elif IsList(r)  then return List( r, x -> l*x );

    # handle the case <scalar> * <polynomial>
    elif not (IsRec(l) and IsBound(l.isPolynomial) and l.isPolynomial)  then

        # <r> must have the cyclotomics as base ring
        if Cyclotomics <> r.baseRing  then
            Error( "<r> must have the cyclotomics as base ring" );
        fi;

        # <l> must lie in the base ring of <r>
        if not IsCyc(l)  then
            Error( "<l> must lie in the base ring of <r>" );
        fi;

        # compute the product
        if l = 0 or r.coefficients = []  then
            prd := [];
            val := 0;
        else
            prd := l * r.coefficients;
            val := r.valuation;
        fi;
        prd            := CyclotomicsOps.FastPolynomial( Cyclotomics, prd, val );
        prd.depth      := 1;
        prd.groundRing := Cyclotomics;

        # and return
        return prd;

    # handle the case <polynomial> * <scalar>
    elif not (IsRec(r) and IsBound(r.isPolynomial) and r.isPolynomial)  then

        # <l> must have the cyclotomics as base ring
        if Cyclotomics <> l.baseRing  then
            Error( "<l> must have the cyclotomics as base ring" );
        fi;

        # <r> must lie in the base ring of <r>
        if not IsCyc(r)  then
            Error( "<r> must lie in the base ring of <l>" );
        fi;

        # compute the product
        if r = Zero( l.baseRing ) or l.coefficients = []  then
            prd := [];
            val := 0;
        else
            prd := l.coefficients * r;
            val := l.valuation;
        fi;
        prd            := CyclotomicsOps.FastPolynomial( Cyclotomics, prd, val );
        prd.depth      := 1;
        prd.groundRing := Cyclotomics;

        # and return
        return prd;
    fi;

    # our superclass will handle different depth
    if not IsBound(l.depth)  then  l.operations.Depth(l);  fi;
    if not IsBound(r.depth)  then  r.operations.Depth(r);  fi;
    if 1 <> l.depth or 1 <> r.depth  then
        return PolynomialOps.\*( l, r );

    # give up if we have rings other then the cyclotomics
    elif Cyclotomics <> l.baseRing  then
        Error( "<l> must have the cyclotomics as base ring" );
    elif Cyclotomics <> r.baseRing  then
        Error( "<r> must have the cyclotomics as base ring" );

    # if <l> is the null polynomial return <l>
    elif Length(l.coefficients) = 0  then
        return l;

    # if <r> is the null polynomial return <r>
    elif Length(r.coefficients) = 0  then
        return r;

    # multiply two polynomials
    else

        # get a common ring
        R := l.baseRing;

        # use 'ProductCoeffs' in order to fold product
        prd := ProductCoeffs( l.coefficients, r.coefficients );
        val := l.valuation + r.valuation;

        # compute the product
        prd            := R.operations.FastPolynomial( R, prd, val );
        prd.depth      := 1;
        prd.groundRing := Cyclotomics;

        # and return
        return prd;
    fi;

end;


#############################################################################
##
#F  CyclotomicsPolynomialOps.\mod  . . . . . . .  remainder of two polynomials
##
CyclotomicsPolynomialOps.\mod := function( l, r ) local  R,  rem,  val,  vdf;

    # <l> must be a polynomial
    if not (IsRec(l) and IsBound(l.isPolynomial) and l.isPolynomial)  then
        Error( "<l> must be a polynomial" );
    fi;

    # if <r> is a integer reduce the coefficients of <l>
    if IsInt(r)  then
        rem := Copy(l.coefficients);
        ReduceCoeffsMod( rem, r );
        return CyclotomicsOps.FastPolynomial( l.baseRing, rem, l.valuation );
    fi;


    # otherwise <r> must be a non-zero polynomial
    if not (IsRec(r) and IsBound(r.isPolynomial) and r.isPolynomial)  then
        Error( "<r> must be a polynomial" );
    fi;
    if Length(r.coefficients) = 0  then
        Error( "<r> must be non zero" );
    fi;

    # our superclass will handle different depth
    if not IsBound(l.depth)  then  l.operations.Depth(l);  fi;
    if not IsBound(r.depth)  then  r.operations.Depth(r);  fi;
    if 1 <> l.depth or 1 <> r.depth  then
        return PolynomialOps.\mod( l, r );

    # give up if we have different rings
    elif l.baseRing <> r.baseRing  then
        Error( "polynomials must have the same ring" ); 

    # reduce the polynomial <l> by <r>
    else

        # if one is a Laurent polynomial use 'EuclideanRemainder'
        if l.valuation < 0 or r.valuation < 0  then
            return EuclideanRemainder( DefaultRing(l,r), l, r );
        fi;

        # get a common ring and the value difference
        R   := l.baseRing;
        vdf := r.valuation - l.valuation;

        # if <r>.valuation is the minimum shift <l>
        if r.valuation < l.valuation  then
            val := r.valuation;
            rem := ShiftedCoeffs( l.coefficients, -vdf );
            ReduceCoeffs( rem, r.coefficients );

        # if <l>.valuation is the minimum shift <r>
        elif l.valuation < r.valuation  then
            r   := ShiftedCoeffs( r.coefficients, vdf );
            rem := RemainderCoeffs( l.coefficients, r );
            val := l.valuation;

        # otherwise they are equal
        else
            rem := RemainderCoeffs( l.coefficients, r.coefficients );
            val := l.valuation;
        fi;

        # return the remainder
        rem := R.operations.FastPolynomial( R, rem, val );
        rem.depth := 1;
        return rem;
    fi;
end;


#############################################################################
##
#F  CyclotomicsPolynomialOps.\^  . . . . . . . . . . .  power of a polynomials
##
CyclotomicsPolynomialOps.\^ := function( l, r ) local   R,  pow, val;

    # <l> must be a polynomial over the cyclotomics and <r> an integer
    if not (IsRec(l) and IsBound(l.isPolynomial) and l.isPolynomial)  then
        Error( "<l> must be a polynomial" );
    elif l.baseRing <> Cyclotomics  then
        Error( "<l> must be a polynomial over the cyclotomics" );
    fi;
    if not IsInt(r)  then
        Error( "<r> must be an integer" );
    fi;

    # invert <l> if necessary
    if r < 0  then
        R := LaurentPolynomialRing( l.baseRing );
        l := R.operations.Quotient( R, One( R ), l );
        r := -r;
    fi;

    # if <r> is zero, return x^0
    if r = 0  then
        return CyclotomicsOps.FastPolynomial( Cyclotomics, [ 1 ], 0 );

    # if <r> is one return <l>
    elif r = 1  then
        return l;

    # if <l> is trivial return
    elif Length(l.coefficients) = 0  then
        return l;

    # if <l> is of degree less than 2, return
    elif Length(l.coefficients) = 1  then
        return CyclotomicsOps.FastPolynomial(
                   Cyclotomics,
                   [l.coefficients[1]^r],
                   l.valuation*r );
    fi;

    # use repeated squaring
    val := l.valuation * r;
    pow := [ One( l.baseRing ) ];
    l   := l.coefficients;
    while 0 < r  do
        if r mod 2 = 1  then
            pow := ProductCoeffs( pow, l );
            r   := r - 1;
        fi; 
        if 1 < r  then
            l := ProductCoeffs( l, l );
            r := r / 2;
        fi;
    od;

    # return the power
    return CyclotomicsOps.FastPolynomial( Cyclotomics, pow, val );

end;


#############################################################################
##
#F  CyclotomicsPolynomialOps.String( <f> )  . . . . construct a pretty string
##
CyclotomicsPolynomialOps.String := function( f )
    local   x,  i,  d,  v,  s,  l, hi;

    # find a name for the indeterminate
    x := Indeterminate(f.baseRing);
    if IsBound(x.name)  then x := x.name;  else x := "x";  fi;

    # run through the coefficients of <f>
    v := f.valuation-1;
    l := Length(f.coefficients);
    for i  in Reversed([ 1 .. l ])  do
      d := f.coefficients[i];
      if 0 <> d  then
	if i = l and d = 1 and i+v <> 0  then s := "";
	elif i = l and d = 1  then s := "1";
	elif i = l and d = -1 and i+v <> 0  then s := "-";
	elif i = l and d = -1  then s := "-1";
	elif i = l  then s := String(d);
	elif d = 1 and i+v <> 0  then s := ConcatenationString( s, "+" );
	elif d = 1  then s := ConcatenationString( s, "+1" );
	elif d = -1 and i+v <> 0  then s := ConcatenationString( s, "-" );
	elif d = -1  then s := ConcatenationString( s, "-1" );
	elif d < 0  then s := ConcatenationString( s, String(d) );
	elif 0 < d  then s := ConcatenationString( s, "+", String(d) );
	else Error( "internal error in 'CyclotomicsPolynomialOps.String'" );
	fi;
# some changes such that we get "2*q+1" instead of "2q+1":
	if i+v < 0 or 1 < i+v then
	    hi := ConcatenationString(  x, "^", String(i+v) );
	elif i+v = 1  then hi := x;
	else hi := "";
	fi;
	if s<>"" and hi<>"" and 
	  s[Length(s)] in "0123456789" then
	  hi := ConcatenationString( "*", hi );
	fi;
	s:= ConcatenationString( s, hi ); 
      fi;
    od;

    # catch a special case
    if l = 0  then s := "0";  fi;
    return s;
end;

CyclotomicsPolynomialOps.New := function( R, coeffs, val )local res;
    if 0 = Length(coeffs)  then val := 0;  fi;
    res:=rec( coefficients := coeffs,
	      baseRing     := R,
	      isPolynomial := true,
	      valuation    := val,
	      operations   := CyclotomicsPolynomialOps );
    if val < 0  then res.domain:=LaurentPolynomials;
    else             res.domain:=CyclotomicsPolynomials;
    fi;
    return res;
end;

#############################################################################
##

#V  CyclotomicsPolynomials  . . . . .  domain of polynomials over the cyclotomics
##
CyclotomicsPolynomials            := Copy( Polynomials );
CyclotomicsPolynomials.name       := "CyclotomicsPolynomials";
CyclotomicsPolynomialsOps         := OperationsRecord
    ( "CyclotomicsPolynomialsOps", FieldPolynomialRingOps );
CyclotomicsPolynomials.operations := CyclotomicsPolynomialsOps;

# show that this a polynomial ring
CyclotomicsPolynomials.isPolynomialRing := true;

# cyclotomics polynomials form a ring
CyclotomicsPolynomials.isDomain := true;
CyclotomicsPolynomials.isRing   := true;

# set known properties
CyclotomicsPolynomials.isFinite := false;
CyclotomicsPolynomials.size     := "infinity";

# add properties of polynom ring over a field
CyclotomicsPolynomials.isCommutativeRing         := true;
CyclotomicsPolynomials.isIntegralRing            := true;
CyclotomicsPolynomials.isUniqueFactorizationRing := true;
CyclotomicsPolynomials.isEuclideanRing           := true;

# set one, zero and base ring
CyclotomicsPolynomials.one  := Polynomial( Cyclotomics, [1] );
CyclotomicsPolynomials.zero := Polynomial( Cyclotomics, [] );
CyclotomicsPolynomials.baseRing := Cyclotomics;


#############################################################################
##
#F  CyclotomicsPolynomialsOps.\= . . . . . . . . . . . . . . . . equaltity test
##
CyclotomicsPolynomialsOps.\= := function( R, S )
    
    # both rings must be full polynomial rings
    if not IsPolynomialRing(R) or not IsPolynomialRing(S)  then
        return RingOps.\=( S, R );
        
    # compare the base rings in this case
    else
        return R.baseRing = S.baseRing;
    fi;
    
end;

#############################################################################
##
#F  CyclotomicsPolynomialsOps.\in  . . . . . . . . . . . . . .  membership test
##
CyclotomicsPolynomialsOps.\in := function( p, CyclotomicsPolynomials )
    return     IsRec( p )
           and IsBound( p.isPolynomial )
           and p.isPolynomial
           and IsField( p.baseRing )
           and 0 <= p.valuation
           and p.baseRing = Cyclotomics;
end;

#############################################################################
##
#F  CyclotomicsPolynomialsOps.DefaultRing( <L> )  . . . . . . . .  default ring
##
CyclotomicsPolynomialsOps.DefaultRing := PolynomialsOps.DefaultRing;

#############################################################################
##
#F  CyclotomicsPolynomialsOps.EuclideanRemainder( <R>, <f>, <g> ) . . . . . rem
##
CyclotomicsPolynomialsOps.EuclideanRemainder := function( R, f, g )
    return f mod g;
end;
