#############################################################################
##
#A  decoders.g              GUAVA library                       Reinald Baart
#A                                                        &Jasper Cramwinckel
#A                                                           &Erik Roijackers
##
##  This file contains functions for decoding codes
##
#H  $Log: decoders.g,v $
#H  Revision 1.2  1997/01/20 15:06:10  werner
#H  Upgrade from Guava 1.2 to Guava 1.3 for GAP release 3.4.4.
#H
#H  Revision 1.2  1994/10/13  15:17:58  rbaart
#H  Changed codeword functions
#H
#H  Revision 1.1  1994/09/28  09:50:09  jcramwin
#H  Initial revision
#H
##

#############################################################################
##
#F  BCHDecoder( <C>, <r> )  . . . . . . . . . . . . . . . . decodes BCH codes
##
BCHDecoder := function (C, r)
    local F, q, n, m, ExtF, x, a, t, ri_1, ri, rnew, si_1, si, snew,
          ti_1, ti, qi, sigma, i, cc, cl, mp, ErrorLocator, zero,
          Syndromes, null, pol, ExtSize, ErrorEvaluator, Fp;
    F := Field(C);
    q := Size(F);
    n := WordLength(C);
    m := OrderMod(q,n);
    t := QuoInt(C.designedDistance - 1, 2);
    ExtF := GF(q^m);
    x := X(ExtF);
    a := PrimitiveUnityRoot(q,n);
    zero := ExtF.zero;
    r := PolyCodeword(r, F, n);
    if Value(GeneratorPol(C), a) <> zero then
        return CycCodeOps.Decode(C, r);
    fi;
    # Calculate syndrome: this simple line is faster than using minimal pols.
    Syndromes :=  List([1..2*QuoInt(C.designedDistance - 1,2)],
                       i->Value(r, a^i));
    if Maximum(Syndromes) = F.zero then # no errors
        return Codeword(r / GeneratorPol(C), C);
    fi;
    # Use Euclidean algorithm:
    ri_1 := x^(2*t); ri := Polynomial(ExtF, Syndromes); rnew := Copy(ri);
    si_1 := x^0; si := 0*x; snew := 0*x;
    ti_1 := 0*x; ti := x^0; sigma := x^0;
    while Degree(rnew) >= t do
        rnew := Mod(ri_1,ri);
        qi := (ri_1 - rnew) / ri;
        snew := si_1 - qi*si;
        sigma := ti_1 - qi*ti;
        ri_1 := ri; ri := rnew;
        si_1 := si; si := snew;
        ti_1 := ti; ti := sigma;
    od;
    # Chien search for the zeros of error locator polynomial:
    ErrorLocator := []; 
    null := a^0;
    ExtSize := q^m-1;
    for i in [0..ExtSize-1] do
        if Value(sigma, null) = zero then
            AddSet(ErrorLocator, (ExtSize-i) mod n);
        fi;
        null := null * a;
    od;
    # And decode:
    if Length(ErrorLocator) = 0 then
        Error("not decodable");
    fi;
    x := X(F);
    if q = 2 then # error locator is not necessary
        pol := Sum(List([1..Length(ErrorLocator)], i->x^ErrorLocator[i]));
        return Codeword((r - pol) / GeneratorPol(C), C);
    else
        pol := Derivative(sigma);
        Fp := F.one*(x^n-1);
        ErrorEvaluator := List(ErrorLocator,i->
                              Value(rnew,a^-i)/Value(pol, a^-i));
        pol := Sum(List([1..Length(ErrorLocator)], i->
                       -ErrorEvaluator[i]*x^ErrorLocator[i]));
        return Codeword((r - pol) / GeneratorPol(C), C);
    fi;
end;

#############################################################################
##
#F  HammingDecoder( <C>, <r> )  . . . . . . . . . . . . decodes Hamming codes
##
##  Generator matrix must have all unit columns
###############################################################
HammingDecoder := function(C, r)
    local H, S,p, F, fac, e,z,x,ind, i,Sf;
    S := VectorCodeword(Syndrome(C,r));
    r := VectorCodeword(r);
    F := Field(C);
    p := PositionProperty(S, s->s<>F.zero);
    if p <> false then
        z := Z(CharFFE(S[p]))^0;
        if z = S[p] then
            fac := F.one;
        else
            fac := S[p]/z;
        fi;
        Sf := S/fac;
        H := CheckMat(C);
        ind := [1..WordLength(C)];
        for i in [1..Redundancy(C)] do
            ind := Filtered(ind, j-> H[i][j] = Sf[i]);
        od;
        e := ind[1];
        r[e] := r[e]-fac;     # correct error
    fi;
    x := SolutionMat(GeneratorMat(C), r);
    return Codeword(x);
end;	

