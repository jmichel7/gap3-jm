#############################################################################
##
#A  recsp.g                     GAP library                      Frank Celler
##
#H  @(#)$Id: recsp.g,v 1.1 1997/03/10 13:49:29 gap Exp $
##
#Y  Copyright (C) 1996,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This  file  contains functions  which will  help to recognize irreduzible
##  groups containing SP(n,q).
##
Revision_recsp_g :=
    "@(#)$Id: recsp.g,v 1.1 1997/03/10 13:49:29 gap Exp $";


#############################################################################
##

#V  InfoRecSP?(...) . . . . . . . . . . . . . . . . . . . . . . .  infomation
##
##  InfoRecSP4: runtime information
##  InfoRecSP5: MeatAxe
##
if not IsBound(InfoRecSP1)   then InfoRecSP1  := Ignore;  fi;
if not IsBound(InfoRecSP2)   then InfoRecSP2  := Ignore;  fi;
if not IsBound(InfoRecSP3)   then InfoRecSP3  := Ignore;  fi;
if not IsBound(InfoRecSP4)   then InfoRecSP4  := Ignore;  fi;
if not IsBound(InfoRecSP5)   then InfoRecSP5  := Ignore;  fi;


#############################################################################
##

#V  RecSP . . . . . . . . . . . . . . . . . .  functions to recognise SP(n,q)
##
RecSP := rec();
RecSP.STAT_ALT       :=  1;
RecSP.STAT_CHEV      :=  2;
RecSP.STAT_LARGER    :=  3;
RecSP.STAT_PGROUP    :=  4;
RecSP.STAT_POWER     :=  5;
RecSP.STAT_PRIMITIVE :=  6;
RecSP.STAT_PRODUCT   :=  7;
RecSP.STAT_SMALLER   :=  8;
RecSP.STAT_SPORADIC  :=  9;
RecSP.STAT_FORM      := 10;


#############################################################################
##

#F  RecSP.Random( <G> ) . . . . . . . . . . .  return a random element of <G>
##
RecSP.Random := RecSL.Random;


#############################################################################
##
#F  RecSP.SymplecticForm( <slpos>, <pos> )  . . try to find a symplectic form
##
RecSP.SymplecticForm := function( slpos, pos )
    if not (IsBound(slpos.symmetricForm) or IsBound(slpos.symplecticForm))
    then
        RecSL.InvariantForm(slpos);
    fi;
    if IsBound(slpos.quadraticForm)  then
        pos.quadraticForm := slpos.quadraticForm;
    fi;
    if IsBound(slpos.symplecticForm)  then
        pos.symplecticForm  := slpos.symplecticForm;
        pos.gScalars        := slpos.symplecticScalars;
        pos.isMeataxeLarger := false;
        return pos.symplecticForm;
    else
        return false;
    fi;
end;


#############################################################################
##

#F  RecSP.SetAlternating( <pos>, <slpos> )  . . . . . . set alternating group
##
RecSP.SetAlternating := function( pos, slpos )
    local   grp;

    # use alternating groups still possible
    if not slpos.isAlternating  then
        pos.alternating   := [];
        pos.isAlternating := false;
        return;
    fi;

    # catch the special case SP(4,2) = S_6
    if pos.d = 4 and pos.q = 2  then
        if IsGModule(pos.group)  then
            grp := ApplyFunc( Group, GeneratorsFlag(pos.group) );
        else
            grp := pos.group;
        fi;
        if Size(grp) = 720  then
            pos.alternating   := [];
            pos.isAlternating := false;
            return;
        fi;
    fi;

    # check group order
    if IsBound(pos.orderSP)  then
        pos.alternating := Filtered( slpos.alternating, x -> 2*pos.orderSPS
                                     mod Factorial(x) = 0 );
    else
        pos.alternating := slpos.alternating;
    fi;
    pos.isAlternating := 0 < Length(pos.alternating);

end;


#############################################################################
##
#F  RecSP.CheckAlternating( <pos>, <pord>, <psod> ) . . .  element order test
##
RecSP.CheckAlternating := function( pos, pord, psod )
    local   z,  new,  n,  z2;

    # are alternating groups still possible?
    if not pos.isAlternating then return;  fi;
    pos.statistic[RecSP.STAT_ALT] := pos.tries;

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
        InfoRecSP2( "#I  <G> is not an alternating group,  ",
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
        InfoRecSP2( "#I  <G> is not an alternating group,  ",
                    "element order criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSP.SetChevalley( <pos>, <slpos> )  . . . . . possible Chevalley groups
##
RecSP.SetChevalley := function( pos, slpos )
    local   new,  x;

    # use Chevalley groups still possible
    if not slpos.isChevalley  then
        pos.expsChev    := [];
        pos.isChevalley := false;
        return;
    fi;

    # check orders
    if IsBound(pos.orderSP)  then
        new := Filtered( slpos.expsChev, x -> pos.orderSPS mod
                                  x[6].order(x[3],x[4],x[5]) = 0 );
    else
        new := slpos.expsChev;
    fi;

    # avoid certain small cases
    pos.expsChev := [];
    for x  in new  do

        # get rid of ChevC
        if x[6] = ChevC  then
            if 2*x[3] <> pos.d or x[4] <> pos.p or x[5] <> pos.k  then
                Add( pos.expsChev, x );
            fi;

        # A5 < S6=SP(4,2) is check in alternating test
        elif x[6] = ChevA and pos.d = 4 and pos.q = 2  then
            InfoRecSP3( "#I  ignoring Chevalley groups in S6\n" );

        # get rid of SL(d/2,q^2) (MEATAXE test!)
        elif x[6] = ChevA  then
            if 2*(x[3]+1)<>pos.d or x[4]<>pos.p or x[5]<>2*pos.k  then
                Add( pos.expsChev, x );
            fi;

        # S4(3) = U4(2)
        elif x[6] = Chev2A and pos.q = 3 and pos.d = 4  then
            if x[3] <> 3 or x[4] <> 2 or x[5] <> 1  then
                Add( pos.expsChev, x );
            fi;

        # ok, add me
        else
            Add( pos.expsChev, x );
        fi;

    od;

    # is it possible
    pos.isChevalley := 0 < Length(pos.expsChev);

end;


#############################################################################
##
#F  RecSP.CheckChevalley( <pos>, <pord>, <psod> ) . . . . . . . element order
##
RecSP.CheckChevalley := function( pos, pord, psod )

    # are Chevalley groups still possible?
    if not pos.isChevalley  then return;  fi;
    pos.statistic[RecSP.STAT_CHEV] := pos.tries;

    pos.expsChev := Filtered(pos.expsChev, x -> x[1] mod pord = 0);
    if 1 < psod  then
        pos.expsChev := Filtered(pos.expsChev, x -> 1 < GcdInt(x,psod));
    fi;
    if 0 = Length(pos.expsChev)  then
        pos.isChevalley := false;
        InfoRecSP2( "#I  <G> is not a Chevalley group,  ",
                    "element order criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSP.SetLargerField( <pos>, <slpos> )  . . . . . . possible larger field
##
##  d = e*f, G < TSp(e,q^f) = Sp(e,q^f).f
##  or
##  GU(d/2,q^2).2 if q is odd
##
RecSP.SetLargerField := function( pos, slpos )
    local   d;

    # check if we already know it
    if not slpos.isLarger   then
        pos.isLarger := false;
        return;
    fi;

    # no reduction is possible if dimension is one
    pos.isLarger := 1 < pos.d;
    if not pos.isLarger  then
        InfoRecSP2("#I  <G> is not definable over a larger field,  ",
                   "dimension is one\n" );

    # representation that preserve a larger field
    else
        d := Filtered( DivisorsInt(pos.d), x -> 1<x and (pos.d/x) mod 2=0 );
        d := List( d, x -> RecSL.SP.uexponent(pos.d/x,pos.p,x*pos.k)*x );
        pos.expsLargerSp := d;
        if pos.q mod 2 = 1  then
            pos.expLargerGU := RecSL.GU.uexponent(pos.d/2,pos.p,pos.k)*2;
        else
            pos.expLargerGU := false;
        fi;
    fi;

end;


#############################################################################
##
#F  RecSP.CheckLargerField( <pos>, <pord> ) . . . . . . . check element order
##
RecSP.CheckLargerField := function( pos, pord )
    local   c,  I,  g,  d,  e;

    # are larger fields still possible?
    if not pos.isLarger  then return;  fi;
    pos.statistic[RecSP.STAT_LARGER] := pos.tries;

    # check element orders
    pos.expsLargerSp := Filtered(pos.expsLargerSp, x -> x mod pord = 0);
    if pos.expLargerGU <> false  then
        if pos.expLargerGU mod pord <> 0  then
            pos.expLargerGU := false;
        fi;
    fi;
    if 0 = Length(pos.expsLargerSp) and pos.expLargerGU = false  then
        pos.isLarger := false;
        InfoRecSP2( "#I  <G> is not definable over a larger field,  ",
                    "element order criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSP.SetMysteriousPGroup( <pos>, <slpos> ) . .  is mysterious p possible
##
##  R < G < N,  where N is the normalizer of an extraspecial r-group R.
##
RecSP.SetMysteriousPGroup := function( pos, slpos )
    local   r,  m,  e;

    # mysterious p group: dimension must be a prime power
    pos.isMysteriousP := false;
    if pos.k = 1 and IsPrimePowerInt(pos.d)  then
        r := SmallestRootInt( pos.d );
        m := LogInt( pos.d, r );
        if r = 2 and r <> pos.p  then
            e := 1;
            while pos.p^e mod r <> 1  do e := e + 1;  od;
            if e = 1  then
                pos.isMysteriousP  := true;
                pos.expMysteriousP := 4*RecSL.O.uexponent(-1,2*m,r,1);
            fi;
        fi;
    fi;
    if not pos.isMysteriousP  then
        InfoRecSP2( "#I  <G> is no mysterious p-group,  dimension ",
                    "is not a prime power\n" );
    fi;
    
end;


#############################################################################
##
#F  RecSP.CheckMysteriousPGroup( <pos>, <pord> )  . . . . check element order
##
##  <pord> is the correct part of the projective order.
##
RecSP.CheckMysteriousPGroup := function( pos, pord )

    # mysterious p-group still possible?
    if not pos.isMysteriousP  then return;  fi;
    pos.statistic[RecSP.STAT_PGROUP] := pos.tries;

    # check element order
    pos.isMysteriousP := pos.expMysteriousP mod pord = 0;
    if not pos.isMysteriousP  then
        InfoRecSP2( "#I  <G> is no mysterious p-group,  element order ",
                    "criteria failed\n" );
    fi;
    
end;


#############################################################################
##
#F  RecSP.SetImprimitive( <pos>, <slpos> )  . . . possible imprimitive groups
##
##   d = e*f,  G < Sp(e,q) wr Sym(f) or
##  or
##   GL(d/2,q).2
##
RecSP.SetImprimitive := function( pos, slpos )
    local   d;

    # check if we already know it
    if not slpos.isImprimitive   then
        pos.isImprimitive := false;
        return;
    fi;

    # store the various pairs <e>, <f> for symplectic decomposition
    d := DivisorsInt(pos.d);
    if pos.q = 2  then
        d := Filtered( d, x -> x <> pos.d and x mod 2 = 0 and x <> 2 );
    else
        d := Filtered( d, x -> x <> pos.d and x mod 2 = 0 );
    fi;
    pos.dimsImprimitive := List( d, x -> [ x, pos.d/x ] );

    # store exponent of totally singular decomposition
    if pos.p <> 2  then
        pos.expImprimitive := RecSL.PGL.uexponent(pos.d/2,pos.p,pos.k)
                              *(pos.q-1);
    else
        pos.expImprimitive := false;
    fi;
    pos.isImprimitive := true;

end;


#############################################################################
##
#F  RecSP.CheckImprimitive( <pos>, <pord> ) . . . . . . . check element order
##
##  Assume that all integer <= <pos.d> are no semi primes.
##
RecSP.CheckImprimitive := function( pos, pord )
    local   new,  p,  m;

    # imprimitivity impossible?
    if not pos.isImprimitive  then return;  fi;
    pos.statistic[RecSP.STAT_PRIMITIVE] := pos.tries;

    # check totally singular decomposition
    if pos.expImprimitive <> false  then
        if pos.expImprimitive mod pord <> 0  then
            pos.expImprimitive := false;
        fi;
    fi;
    
    # check exponents of symplectic decomposition
    new := [];
    for p  in pos.dimsImprimitive  do
        m := pord / Gcd( pord, pos.expsSP[p[1]] );
        if Factorial(p[2]) mod m = 0 
           and Sum(Collected(Factors(m)),x->x[1]^x[2]) <= p[2]
        then
            Add( new, p );
        fi;
    od;
    pos.dimsImprimitive := new;

    # if <new> is trivial no reduction is possible
    if 0 = Length(new) and pos.expImprimitive = false  then
        pos.isImprimitive := false;
        InfoRecSP2( "#I  <G> is not imprimitive,  element order ",
                    "criteria failed\n" );
    fi;
    
end;


#############################################################################
##
#F  RecSP.SetSmallerField( <pos>, <slpos> ) . . . . .  possible smaller field
##
##  If <q> = <p>^<k> is a prime, then the matrices are already written over a
##  prime  field,  no  reduction  is possible in this case.   Otherwise upper
##  bounds for the exponents  of PSP(<n>,<p>^<i>),  <i> T <k>,  are computed.
##  The record component <pos>.smallerField containes the smaller field still
##  possible.
##
RecSP.SetSmallerField := function( pos, slpos )
    local   i,  j,  q,  exp;

    # check if we already know it
    if not slpos.isSmaller   then
        pos.isSmaller:= false;
        return;
    fi;

    # if <field> is the prime field,  no reduction is possible
    pos.isSmaller := 1 < pos.k;
    if not pos.isSmaller  then
        InfoRecSP2( "#I  <G> is not definable over a smaller field,  ",
                    "field is the prime field\n" );

    # loop over the maximal divisors of <k>
    else
        pos.expsSmaller := [];
        for i in DivisorsInt(pos.k) do
            if i < pos.k and IsPrime(pos.k/i)  then
                AddSet( pos.expsSmaller,RecSL.SP.uexponent(pos.d,pos.p,i));
            fi;
        od;
        pos.smallerField := GF(pos.p);
    fi;         

end;        


#############################################################################
##
#F  RecSP.CheckSmallerField( <pos>, <cpol>, <pord> )  check order and charpol
##
##  <pord> is the projective order  of the  group  element  A, <cpol>  is the
##  characteristic polynomial  of A.  The function  first  checks  if  <pord>
##  divides  at least one  of the exponents stored in  <pos>.expsSmaller.  It
##  then computes the smallest field which contains <cpol> * zeta.
##
RecSP.CheckSmallerField := function( pos, cpol, pord )
    local   c,  d,  t,  p,  i0,  I;

    # are smaller fields still possible?
    if not pos.isSmaller  then return;  fi;
    pos.statistic[RecSP.STAT_SMALLER] := pos.tries;

    # first check the projective order
    pos.expsSmaller := Filtered( pos.expsSmaller, x -> x mod pord = 0 );
    if 0 = Length(pos.expsSmaller)  then
        pos.isSmaller := false;
        InfoRecSP2( "#I  <G> is not definable over a smaller field,  ",
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
            InfoRecSP2( "#I  <G> is not definable over a smaller field, ",
                        "char polynomial failed\n" );
        fi;
    fi;

end;


#############################################################################
##
#F  RecSP.SetSporadicGroups( <pos>, <slpos> ) . . .  possible sporadic groups
##
RecSP.SetSporadicGroups := function( pos, slpos )

    # use sporadics still possible
    if not slpos.isSporadic  then
        pos.sporadicGroups := [];
        pos.isSporadic     := false;
        return;
    fi;
    pos.sporadicGroups := slpos.sporadicGroups;
    pos.isSporadic     := true;

end;


#############################################################################
##
#F  RecSP.CheckSporadicGroups( <pos>, <pord> )  . . . . .  element order test
##
RecSP.CheckSporadicGroups := function( pos, pord )
    
    # sporadic groups still possible?
    if not pos.isSporadic  then return;  fi;
    pos.statistic[RecSP.STAT_SPORADIC] := pos.tries;

    # check the possible element orders of the automorphism groups
    pos.sporadicGroups := Filtered( pos.sporadicGroups, x -> pord in
                                    SporadicGroupsInfo.orders[x] );

    # if the list is trivial reset '<pos>.isSporadic'
    if 0 = Length(pos.sporadicGroups)  then
        pos.isSporadic := false;
        InfoRecSP2( "#I  <G> is not a sporadic group,  ",
                    "element order criteria failed\n" );
    fi;

    # if the exponent is too big stop
    if pos.exponent >= pos.q3d  then
        pos.isSporadic := false;
        InfoRecSP2( "#I  <G> is not a sporadic group, ",
                    "group exponent criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSP.SetTensorPowers( <pos>, <slpos> ) . . . . .  possible tensor powers
##
##  d  = m^t, G <  Sp(m,q) wr Sym(t).   qt odd, m  even, t >  2  and (q,m) !=
##  (3,2). Assume that there are no semi primes smaller then t+1.
##
RecSP.SetTensorPowers := function( pos, slpos )
    local   m,  t,  i;

    # check if we already know it
    if not slpos.isTensorPower   then
        pos.isTensorPower:= false;
        return;
    fi;

    # no reduction is possible if <t> is 1
    m := SmallestRootInt( pos.d );
    if m = pos.d  then
        pos.isTensorPower := false;
        InfoRecSP2( "#I  <G> is no tensor power,  dimension is ",
                    "not a power\n" );
        
    # if <m> is odd no reduction is possible
    elif m mod 2 = 1  then
        pos.isTensorPower := false;
        InfoRecSP2( "#I  <G> is no tensor power, dimension is an ",
                    "odd power\n" );

    # if <q> is even no reduction is possible
    elif pos.q mod 2 = 0  then
        pos.isTensorPower := false;
        InfoRecSP2( "#I  <G> is no tensor power, characteristic is 2\n" );

    # store the various pairs <t>, <m>
    else
        pos.isTensorPower    := true;
        pos.dimsTensorPowers := [];
        t := LogInt( pos.d, m );
        for i  in DivisorsInt(t)  do
            if 2<t/i and pos.q*t/i mod 2 = 1 and (pos.q<>3 or m^i<>2)  then
                Add( pos.dimsTensorPowers, [ m^i, t/i ] );
            fi;
        od;
    fi;

end;  


#############################################################################
##
#F  RecSP.CheckTensorPowers( <pos>, <ord> ) . . . . . . . check element order
##
##  <ord> is the order of the  group element A. If A lies in GL(e,q) X Sym(f)
##  then gcd(<ord>,exp(GL(e,q))) must be a valid element order in Sym(f).
##        
##  Assume that all integer <= <pos.d> are no semi primes.
##
RecSP.CheckTensorPowers := function( pos, ord )
    local   new,  p,  m;
    
    # are tensor powers still possible?
    if not pos.isTensorPower  then return;  fi;
    pos.statistic[RecSP.STAT_POWER] := pos.tries;

    # check exponents
    new := [];
    for p  in pos.dimsTensorPowers  do
        m := ord / Gcd( ord, pos.expsSP[p[1]] );
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
        InfoRecSP2( "#I  <G> is no tensor power,  element order ",
                    "criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSP.SetTensorProducts( <pos>, <slpos> ) . . .  possible tensor products
##
##  G < Sp(d1,q) x O(d2,q) with d=d1*d2, d1 even, d2 > 2, q odd
##
RecSP.SetTensorProducts := function( pos, slpos )
    local   i,  s,  exp;

    # check if we already know it
    if not slpos.isTensorProduct   then
        pos.isTensorProduct := false;
        return;
    fi;

    # if <pos.d> is a prime, no tensor products are possible
    if IsPrime(pos.d)  then
        pos.isTensorProduct := false;
        InfoRecSP2( "#I  <G> is no tensor product,  dimension is ",
                    "a prime\n" );
        return;

    # if <pos.q> is even, no tensor product are possible
    elif pos.q mod 2 = 0  then
        pos.isTensorProduct := false;
        InfoRecSP2( "#I  <G> is no tensor product, characteristic is ",
                    "two\n" );
        return;
    else
        pos.isTensorProduct := true;
    fi;
        
    # compute the exponents of the tensor products
    pos.expsTensorProducts := [];
    for i  in DivisorsInt(pos.d)  do
        if i mod 2 = 0 and 2 < pos.d/i  then
            if pos.d/i mod 2 = 1  then
                exp := LcmInt( pos.expsSP[i],
                               RecSL.O.uexponent(0,pos.d/i,pos.p,pos.k) );
                Add( pos.expsTensorProducts, [ exp, i ] );
            else
                for s  in [ -1, 1 ]  do
                    exp := LcmInt( pos.expsSP[i],
                               RecSL.O.uexponent(s,pos.d/i,pos.p,pos.k) );
                    Add( pos.expsTensorProducts, [ exp, i ] );
                od;
            fi;
        fi;
    od;

end;        


#############################################################################
##
#F  RecSP.CheckTensorProducts( <pos>, <pord> )  . . . . . check element order
##
RecSP.CheckTensorProducts := function( pos, pord )

    # are tensor products still possible?
    if not pos.isTensorProduct  then return;  fi;
    pos.statistic[RecSP.STAT_PRODUCT] := pos.tries;

    # check if <ord> divides any exponent of a tensor product
    pos.expsTensorProducts := Filtered( pos.expsTensorProducts,
                                        x -> x[1] mod pord = 0 );
    
    # if <pos>.expsTensorProducts is empty no tps are possible
    if 0 = Length(pos.expsTensorProducts)  then
        pos.isTensorProduct := false;
        InfoRecSP2( "#I  <G> is no tensor product,  element order ",
                    "criteria failed\n" );
    fi;

end;


#############################################################################
##

#F  RecSP.Print . . . . . . . . . . . . . . . . . . . . . . nice pretty print
##
RecSP.Print := function( obj )
    local   name;

    if obj.printLevel = 0  then
        Print( "<< SP recognition record >>" );
        return;
    fi;
    Print( "#I  field: ", obj.q, ", dimension: ", obj.d,
           ", number of generators: ", Length(obj.generators), "\n" );

    if IsBound(obj.symplecticForm)  then
        Print( "#I  symplectic form is known\n" );
    else
        Print( "#I  no symplectic form is known\n" );
        Print( "<< SP recognition record >>" );
        return;
    fi;
    if obj.p = 2 and IsBound(obj.quadraticForm)  then
        Print( "#I  quadratic form is known\n" );
        Print( "<< SP recognition record >>" );
        return;
    fi;
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
    if IsBound(obj.isLarger) and obj.isLarger and obj.isMeataxeLarger  then
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
    if IsBound(obj.isTensorPower) and obj.isTensorPower  then
        Print( "#I  <G> could be a tensor power\n" );
    fi;
    if IsBound(obj.isTensorProduct) and obj.isTensorProduct  then
        Print( "#I  <G> could be a tensor product\n" );
    fi;
    if IsBound(obj.containsSP) and obj.containsSP  then
        Print( "#I  <G> contains SP( ",obj.d,", ",obj.q," )\n" );
    else
        Print( "#I  <G> could still be an symplectic group\n" );
    fi;
    Print( "<< SP recognition record >>" );
end;


#############################################################################
##
#F  RecSP.SetPrintLevel( <obj>, <lev> ) . . . . . . . . . . . set print level
##
RecSP.SetPrintLevel := function( obj, lev )
    if 2 < lev or lev < 0  then
        Error( "invalid print level" );
    fi;
    obj.printLevel := lev;
end;


#############################################################################
##
#F  RecSP.Setup( <slpos> )  . . . . . . . . . . . . . . . . set up pos record
##
RecSP.Setup := function( slpos )
    local   pos,  timer,  j,  qi,  i;

    # set up a SP possibility record
    pos            := rec();
    pos.q          := slpos.q;
    pos.p          := slpos.p;
    pos.k          := slpos.k;
    pos.d          := slpos.d;
    pos.q3d        := slpos.q3d;
    pos.group      := slpos.group;
    pos.identity   := slpos.identity;
    pos.field      := slpos.field;
    pos.exponent   := slpos.exponent;
    pos.generators := slpos.generators;
    pos.orders     := Copy(slpos.orders);
    pos.statistic  := [1..10]*0;
    pos.printLevel := 1;
    pos.containsSP := false;
    pos.isRecSP    := true;
    pos.operations := RecSP;
    pos.setupTime  := -Runtime();

    # store <slpos>
    pos.recognizeSL := slpos;
        
    InfoRecSP1( "#I  field: ", pos.q, ", dimension: ", pos.d, 
                ", number of generators: ", Length(pos.generators), "\n" );

    # compute the exponent of SP(i,q) for 1 <= i < d, i|d, i even
    timer  := Runtime();
    pos.expsSP := [];
    j  := 1;
    qi := 1;
    for i  in [ 2, 4 .. pos.d ]  do
        qi := qi * pos.q;
        j  := LcmInt( j, LcmInt( qi-1, qi+1 ) );
        if pos.d mod i = 0  then
            pos.expsSP[i] := j * pos.p^(1+LogInt(i-1, pos.p));
        fi;
    od;
    InfoRecSP4( "#I  exponents: ", Runtime()-timer, " msec\n" );

    # compute the order of SP( d, q ) for small dimensions
    if pos.d < 10  then
        timer := Runtime();
        pos.orderSP  := RecSL.SP.order( pos.d, pos.p, pos.k );
        pos.orderSPS := pos.orderSP * (pos.q-1) / Gcd( 2, pos.q-1 );
        InfoRecSP4( "#I  order sp: ", Runtime()-timer, " msec\n" );
    fi;

    # the current meataxe will rule out larger fields
    pos.isMeataxeLarger := true;

    # set the possible groups contained in <G>
    timer := Runtime();
    pos.operations.SetAlternating     ( pos, slpos );
    pos.operations.SetChevalley       ( pos, slpos );
    pos.operations.SetLargerField     ( pos, slpos );
    pos.operations.SetMysteriousPGroup( pos, slpos );
    pos.operations.SetImprimitive     ( pos, slpos );
    pos.operations.SetSmallerField    ( pos, slpos );
    pos.operations.SetSporadicGroups  ( pos, slpos );
    pos.operations.SetTensorPowers    ( pos, slpos );
    pos.operations.SetTensorProducts  ( pos, slpos );
    InfoRecSP4( "#I  subgroup setup: ", Runtime()-timer, " msec\n" );
    pos.setupTime := Runtime() + pos.setupTime;

    # and return
    return pos;

end;


#############################################################################
##
#F  RecognizeSP( <slpos>, <tries> ) . . . . . . . . . . . .  regonize SP(n,q)
##
RecognizeSP := function( slpos, tries )
    local   pos,  form,  ord,  poo,  o,  A,  mpol,  po,  oe,  pso,  p,  
            oo,  so,  cpol;

    # is this a reentry
    if not IsBound(slpos.isRecSP) or not slpos.isRecSP  then

        # dimension should be at least 4 (otherwise SL=SP)
        if slpos.d = 2  then
            Error( "dimension must be at least 4,  otherwise SP=SL" );
        fi;

        # <G> must have a chance to be the symplectic group
        if not slpos.isSymplectic  then
            Error( "SL symplectic test for the group already failed" );
        fi;

        # dimension must be even
        if slpos.d mod 2 = 1  then
            Error( "dimension must be even" );
        fi;

        # set up a SP possibility record
        pos := RecSP.Setup(slpos);

        # use old orders to rule out possibilities
        pos.tries := 0;
        InfoRecSP2( "#I  checking old element orders\n" );
        for ord  in slpos.orders  do

            # construct order,  projective order and semi order
            poo := ord[1];  pso := ord[2];  oe  := ord[3];  po := poo * pso;
            o   := po*oe;   oo  := poo*oe;  so  := pso;

            # check possibilities (don't check smaller field)
            pos.operations.CheckAlternating     ( pos,    poo, pso );
            pos.operations.CheckChevalley       ( pos,    poo, pso );
            pos.operations.CheckLargerField     ( pos, po          );
            pos.operations.CheckMysteriousPGroup( pos,    poo      );
            pos.operations.CheckImprimitive     ( pos, po          );
            pos.operations.CheckTensorPowers    ( pos, po          );
            pos.operations.CheckTensorProducts  ( pos, po          );
        od;

    # use old recognition record
    else
        pos   := slpos;
        slpos := pos.recognizeSL;
    fi;

    # try to find an invariant symplectic form
    if not IsBound(pos.symplecticForm)  then
        form := pos.operations.SymplecticForm( slpos, pos );

        # return if we failed to find a symplectic form
        if form = false  then
            InfoRecSP1( "#I  failed to find a symplectic form\n" );
            return pos;
        else
            InfoRecSP1( "#I  found invariant symplectic form\n" );
        fi;
    fi;

    # if <q> is even  try to find a invariant quadratic form
    if IsBound(pos.quadraticForm)  then
        InfoRecSP1( "#I  quadratic form is already known\n" );
        return pos;
    fi;

    # set 'containsSP'
    pos.containsSP :=     not pos.isAlternating
                      and not pos.isChevalley
                      and not pos.isImprimitive
                      and not ( pos.isLarger and pos.isMeataxeLarger )
                      and not pos.isMysteriousP
                      and not pos.isSmaller
                      and not pos.isSporadic
                      and not pos.isTensorProduct
                      and not pos.isTensorPower;

    # start trying random non-trivial group elements
    cpol := false;
    po   := false;
    pos.runTime := -Runtime();
    pos.tries   := 0;
    while not pos.containsSP and pos.tries < tries  do

        # find a non-trivial random element of <G>
        repeat
            A := pos.operations.Random(pos.group);
        until A <> pos.identity;
        pos.tries := pos.tries + 1;
        InfoRecSP2("#I  trying ", pos.tries, ".th element of <G>\n");

        # compute the minimal polynomial of <A> and the projective order
        if    pos.isAlternating
           or pos.isChevalley
           or pos.isImprimitive
           or pos.isLarger
           or pos.isMysteriousP
           or pos.isSmaller
           or pos.isSporadic
           or pos.isTensorPower
           or pos.isTensorProduct
        then
            mpol := MinimalPolynomial( FiniteFieldMatrices, A );
            mpol.baseRing := pos.field;
            po := RecSL.OrderScalar( pos, mpol );
            oe := OrderFFE(po[2]);
            po := po[1];
            InfoRecSP2( "#I  projective order = ", 
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

        fi;

        # compute the characteristic polynomial of <A>
        if pos.isSmaller  then
            cpol := CharacteristicPolynomial( FiniteFieldMatrices, A );
            cpol.baseRing := pos.field;
        fi;

        # check all possibilities
        pos.operations.CheckAlternating     ( pos,       poo, pso );
        pos.operations.CheckChevalley       ( pos,       poo, pso );
        pos.operations.CheckLargerField     ( pos,       po       );
        pos.operations.CheckMysteriousPGroup( pos,       poo      );
        pos.operations.CheckImprimitive     ( pos,       po       );
        pos.operations.CheckSmallerField    ( pos, cpol, po       );
        pos.operations.CheckSporadicGroups  ( pos,       po       );
        pos.operations.CheckTensorPowers    ( pos,       po       );
        pos.operations.CheckTensorProducts  ( pos,       po       );

        # set 'containsSP'
        pos.containsSP :=     not pos.isAlternating
                          and not pos.isChevalley
                          and not pos.isImprimitive
                          and not ( pos.isLarger and pos.isMeataxeLarger )
                          and not pos.isMysteriousP
                          and not pos.isSmaller
                          and not pos.isSporadic
                          and not pos.isTensorProduct
                          and not pos.isTensorPower;
    od;

    return pos;

end;

RecogniseSP := RecognizeSP;
