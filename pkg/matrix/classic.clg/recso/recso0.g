#############################################################################
##
#A  recso0.g                    GAP library                      Frank Celler
##
#H  @(#)$Id: recso0.g,v 1.1 1997/03/10 13:49:20 gap Exp $
##
#Y  Copyright (C) 1995,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This  file  contains functions  which will  help to recognize irreduzible
##  groups containing O0(d,q).
##
Revision_recso0_g :=
    "@(#)$Id: recso0.g,v 1.1 1997/03/10 13:49:20 gap Exp $";


#############################################################################
##
#V  RecSO0 . . . . . . . . . . . . . . . . . .  functions to recognise O0(d,q)
##
RecSO0 := Copy(RecSO);


#############################################################################
##

#F  RecSO0.SetAlternating( <pos>, <slpos> ) . . . . . . set alternating group
##
RecSO0.SetAlternating := function( pos, slpos )
    local   timer,  n,  o,  x;

    # use alternating groups still possible
    if not slpos.isAlternating  then
        pos.alternating   := [];
        pos.isAlternating := false;
        return;
    fi;

    # check group order
    if IsBound(pos.orderO)  then
        pos.alternating := Filtered( slpos.alternating, x -> pos.orderOS
                                     mod Factorial(x) = 0 );
    else
        pos.alternating := slpos.alternating;
    fi;
    pos.isAlternating := 0 < Length(pos.alternating);

end;


#############################################################################
##
#F  RecSO0.CheckAlternating( <pos>, <pord>, <psod> )  . .  element order test
##
RecSO0.CheckAlternating := function( pos, pord, psod )
    local   z,  new,  n,  z2;

    # are alternating groups still possible?
    if not pos.isAlternating then return;  fi;
    pos.statistic[RecSO0.STAT_ALT] := pos.tries;

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
        InfoRecSO2( "#I  <G> is not an alternating group,  ",
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
        InfoRecSO2( "#I  <G> is not an alternating group,  ",
                    "element order criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSO0.SetChevalley( <pos>, <slpos> ) . . . . . possible Chevalley groups
##
RecSO0.SetChevalley := function( pos, slpos )
    local   new,  x;

    # use Chevalley groups still possible
    if not slpos.isChevalley  then
        pos.expsChev    := [];
        pos.isChevalley := false;
        return;
    fi;

    # check orders
    if IsBound(pos.orderO)  then
        new := Filtered( slpos.expsChev, x -> pos.orderOS mod
                                  x[6].order(x[3],x[4],x[5]) = 0 );
    else
        new := slpos.expsChev;
    fi;

    # avoid certain small cases
    pos.expsChev := [];
    for x  in new  do

        # get rid of O0
        if x[6] = ChevB  then
            if 2*x[3]+1 <> pos.d or x[4] <> pos.p or x[5] <> pos.k  then
                Add( pos.expsChev, x );
            fi;

        # avoid |SP| = |O0|
        elif x[6] = ChevC  then
            if 2*x[3]+1 <> pos.d or x[4] <> pos.p or x[5] <> pos.k  then
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
#F  RecSO0.CheckChevalley( <pos>, <pord>, <psod> )  . . . . . . element order
##
RecSO0.CheckChevalley := function( pos, pord, psod )

    # are Chevalley groups still possible?
    if not pos.isChevalley  then return;  fi;
    pos.statistic[RecSO0.STAT_CHEV] := pos.tries;

    pos.expsChev := Filtered(pos.expsChev, x -> x[1] mod pord = 0);
    if 1 < psod  then
        pos.expsChev := Filtered(pos.expsChev, x -> 1 < GcdInt(x,psod));
    fi;
    if 0 = Length(pos.expsChev)  then
        pos.isChevalley := false;
        InfoRecSO2( "#I  <G> is not a Chevalley group,  ",
                    "element order criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSO0.SetLargerField( <pos>, <slpos> ) . . . . . . possible larger field
##
##  d = e*f, G < TO(e,q^f)
##
RecSO0.SetLargerField := function( pos, slpos )
    local   new,  d;

    # check if we already know it
    if not slpos.isLarger   then
        pos.isLarger := false;
        return;
    fi;

    # no reduction is possible if dimension is one
    pos.isLarger := 1 < pos.d;
    if not pos.isLarger  then
        InfoRecSO2("#I  <G> is not definable over a larger field,  ",
                   "dimension is one\n" );

    # representation that preserve a larger field
    else
        new := [];
        for d  in Filtered( DivisorsInt(pos.d), x -> 1 < x )  do
            if d mod 2 = 0  then
                Add( new, d*RecSL.O.uexponent(-1,pos.d/d,pos.p,d*pos.k));
                Add( new, d*RecSL.O.uexponent(+1,pos.d/d,pos.p,d*pos.k));
            else
                Add( new, d*RecSL.O.uexponent(0,pos.d/d,pos.p,d*pos.k));
            fi;
        od;
        pos.expsLarger := new;
    fi;

end;


#############################################################################
##
#F  RecSO0.CheckLargerField( <pos>, <pord> )  . . . . . . check element order
##
RecSO0.CheckLargerField := function( pos, pord )

    # are larger fields still possible?
    if not pos.isLarger  then return;  fi;
    pos.statistic[RecSO0.STAT_LARGER] := pos.tries;

    # check element orders
    pos.expsLarger := Filtered( pos.expsLarger, x -> x mod pord = 0 );
    if 0 = Length(pos.expsLarger)  then
        pos.isLarger := false;
        InfoRecSO2( "#I  <G> is not definable over a larger field,  ",
                    "element order criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSO0.SetMysteriousPGroup( <pos>, <slpos> )  .  is mysterious p possible
##
##  there is no mysterious p group in O0
##
RecSO0.SetMysteriousPGroup := function( pos, slpos )
    pos.isMysteriousP := false;
end;


#############################################################################
##
#F  RecSO0.SetImprimitive( <pos>, <slpos> ) . . . possible imprimitive groups
##
##   d = e*f,  1 < e, G < O(e,q) wr Sym(f) or
##  or
##   q=p, G < O(1,q) wr Sym(d)
##
RecSO0.SetImprimitive := function( pos, slpos )
    local   d;

    # check if we already know it
    if not slpos.isImprimitive   then
        pos.isImprimitive := false;
        return;
    fi;

    # store the various pairs <e>, <f> for orthogonal decomposition
    pos.dimsImprimitive := [];
    for d  in Filtered( DivisorsInt(pos.d), x -> 1 < x and x < pos.d )  do
        if d mod 2 = 1  then
            Add( pos.dimsImprimitive, [ d, 0, pos.d/d ] );
        else
            Add( pos.dimsImprimitive, [ d, +1, pos.d/d ] );
            Add( pos.dimsImprimitive, [ d, -1, pos.d/d ] );
        fi;
    od;
    if pos.q = pos.p  then
        Add( pos.dimsImprimitive, [ 1, 0, pos.d ] );
    fi;
    pos.isImprimitive := true;

end;


#############################################################################
##
#F  RecSO0.CheckImprimitive( <pos>, <pord> )  . . . . . . check element order
##
##  Assume that all integer <= <pos.d> are no semi primes.
##
RecSO0.CheckImprimitive := function( pos, pord )
    local   new,  p,  m;

    # imprimitivity impossible?
    if not pos.isImprimitive  then return;  fi;
    pos.statistic[RecSO0.STAT_PRIMITIVE] := pos.tries;

    # check exponents of orthogonal decomposition
    new := [];
    for p  in pos.dimsImprimitive  do
        m := pord / Gcd( pord, pos.expsO[p[2]+2][p[1]] );
        if Factorial(p[3]) mod m = 0
           and Sum(Collected(Factors(m)),x->x[1]^x[2]) <= p[3]
        then
            Add( new, p );
        fi;
    od;
    pos.dimsImprimitive := new;

    # if <new> is trivial no reduction is possible
    if 0 = Length(new)  then
        pos.isImprimitive := false;
        InfoRecSO2( "#I  <G> is not imprimitive,  element order ",
                    "criteria failed\n" );
    fi;
    
end;


#############################################################################
##
#F  RecSO0.SetSmallerField( <pos>, <slpos> )  . . . .  possible smaller field
##
##  If <q> = <p>^<k> is a prime, then the matrices are already written over a
##  prime  field,  no  reduction  is possible in this case.
##
RecSO0.SetSmallerField := function( pos, slpos )
    local   i,  j,  q,  exp;

    # check if we already know it
    if not slpos.isSmaller   then
        pos.isSmaller:= false;
        return;
    fi;

    # if <field> is the prime field,  no reduction is possible
    pos.isSmaller := 1 < pos.k;
    if not pos.isSmaller  then
        InfoRecSO2( "#I  <G> is not definable over a smaller field,  ",
                    "field is the prime field\n" );

    # loop over the maximal divisors of <k>
    else
        pos.expsSmaller := [];
        for i in DivisorsInt(pos.k) do
            if i < pos.k and IsPrime(pos.k/i)  then
                AddSet( pos.expsSmaller, GcdInt(2,pos.k/i)*
                        RecSL.O.uexponent(0,pos.d,pos.p,i) );
            fi;
        od;
        pos.smallerField := GF(pos.p);
    fi;         

end;        


#############################################################################
##
#F  RecSO0.CheckSmallerField( <pos>, <cpol>, <pord> ) check order and charpol
##
##  <pord> is the projective order  of the  group  element  A, <cpol>  is the
##  characteristic polynomial  of A.  The function  first  checks  if  <pord>
##  divides  at least one  of the exponents stored in  <pos.expsSmaller>.  It
##  then computes the smallest field which contains <cpol> * zeta.
##
RecSO0.CheckSmallerField := function( pos, cpol, pord )
    local   c,  d,  t,  p,  i0,  I;

    # are smaller fields still possible?
    if not pos.isSmaller  then return;  fi;
    pos.statistic[RecSO0.STAT_SMALLER] := pos.tries;

    # first check the projective order
    pos.expsSmaller := Filtered( pos.expsSmaller, x -> x mod pord = 0 );
    if 0 = Length(pos.expsSmaller)  then
        pos.isSmaller := false;
        InfoRecSO2( "#I  <G> is not definable over a smaller field,  ",
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
            InfoRecSO2( "#I  <G> is not definable over a smaller field, ",
                        "char polynomial failed\n" );
        fi;
    fi;

end;


#############################################################################
##
#F  RecSO0.SetSporadicGroups( <pos>, <slpos> )  . .  possible sporadic groups
##
RecSO0.SetSporadicGroups := function( pos, slpos )

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
#F  RecSO0.CheckSporadicGroups( <pos>, <pord> ) . . . . .  element order test
##
RecSO0.CheckSporadicGroups := function( pos, pord )
    
    # sporadic groups still possible?
    if not pos.isSporadic  then return;  fi;
    pos.statistic[RecSO0.STAT_SPORADIC] := pos.tries;

    # check the possible element orders of the automorphism groups
    pos.sporadicGroups := Filtered( pos.sporadicGroups, x -> pord in
                                    SporadicGroupsInfo.orders[x] );

    # if the list is trivial reset '<pos>.isSporadic'
    if 0 = Length(pos.sporadicGroups)  then
        pos.isSporadic := false;
        InfoRecSO2( "#I  <G> is not a sporadic group,  ",
                    "element order criteria failed\n" );
    fi;

    # if the exponent is too big stop
    if pos.exponent >= pos.q3d  then
        pos.isSporadic := false;
        InfoRecSO2( "#I  <G> is not a sporadic group, ",
                    "group exponent criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSO0.SetTensorPowers( <pos>, <slpos> )  . . . .  possible tensor powers
##
##  d = m^t, m odd, G < O(m,q) wr Sym(t) Assume that there are no semi primes
##  smaller then t+1.
##
RecSO0.SetTensorPowers := function( pos, slpos )
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
        InfoRecSO2( "#I  <G> is no tensor power,  dimension is ",
                    "not a power\n" );
        
    # if <m> is even no reduction is possible
    elif m mod 2 = 2  then
        pos.isTensorPower := false;
        InfoRecSO2( "#I  <G> is no tensor power, dimension is an ",
                    "even power\n" );

    # store the various pairs <t>, <m>
    else
        pos.isTensorPower    := true;
        pos.dimsTensorPowers := [];
        t := LogInt( pos.d, m );
        for i  in DivisorsInt(t)  do
            if 1 < t/i and (pos.q<>3 or m^i<>3)  then
                Add( pos.dimsTensorPowers, [ m^i, t/i ] );
            fi;
        od;
    fi;

end;  


#############################################################################
##
#F  RecSO0.CheckTensorPowers( <pos>, <ord> )  . . . . . . check element order
##
##  <ord> is the order of the  group element A. If A lies in O(e,q) w Sym(f)
##  then gcd(<ord>,exp(GL(e,q))) must be a valid element order in Sym(f).
##        
##  Assume that all integer <= <pos.d> are no semi primes.
##
RecSO0.CheckTensorPowers := function( pos, ord )
    local   new,  p,  m;
    
    # are tensor powers still possible?
    if not pos.isTensorPower  then return;  fi;
    pos.statistic[RecSO0.STAT_POWER] := pos.tries;

    # check exponents
    new := [];
    for p  in pos.dimsTensorPowers  do
        m := ord / Gcd( ord, pos.expsO[0+2][p[1]] );
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
        InfoRecSO2( "#I  <G> is no tensor power,  element order ",
                    "criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSO0.SetTensorProducts( <pos>, <slpos> )  . .  possible tensor products
##
##  G < O(d1,q) x O(d2,q) with d=d1*d2
##
RecSO0.SetTensorProducts := function( pos, slpos )
    local   i;

    # check if we already know it
    if not slpos.isTensorProduct   then
        pos.isTensorProduct := false;
        return;
    fi;

    # if <pos.d> is a prime, no tensor products are possible
    if IsPrime(pos.d)  then
        pos.isTensorProduct := false;
        InfoRecSO2( "#I  <G> is no tensor product,  dimension is ",
                    "a prime\n" );
        return;
    else
        pos.isTensorProduct := true;
    fi;
        
    # compute the exponents of the tensor products
    pos.expsTensorProducts := [];
    for i  in DivisorsInt(pos.d)  do
        if 2 < i and i*i < pos.d  then
            Add( pos.expsTensorProducts, [ LcmInt(pos.expsO[0+2][pos.d/i],
                    pos.expsO[0+2][i] ), i ] );
        fi;
    od;

end;        


#############################################################################
##
#F  RecSO0.CheckTensorProducts( <pos>, <pord> ) . . . . . check element order
##
RecSO0.CheckTensorProducts := function( pos, pord )

    # are tensor products still possible?
    if not pos.isTensorProduct  then return;  fi;
    pos.statistic[RecSO0.STAT_PRODUCT] := pos.tries;

    # check if <ord> divides any exponent of a tensor product
    pos.expsTensorProducts := Filtered( pos.expsTensorProducts,
                                        x -> x[1] mod pord = 0 );
    
    # if <pos>.expsTensorProducts is empty no tps are possible
    if 0 = Length(pos.expsTensorProducts)  then
        pos.isTensorProduct := false;
        InfoRecSO2( "#I  <G> is no tensor product,  element order ",
                    "criteria failed\n" );
    fi;

end;


#############################################################################
##

#F  RecSO0.Setup( <pos>, <slpos> )  . . . . . . . . . . . . . . finish set up
##
RecSO0.Setup := function( pos, slpos )
    local   timer;

    # set signum
    pos.signum := 0;

    # compute the order of O( d, q ) for small dimensions
    if pos.d < 10  then
        timer := Runtime();
        pos.orderO  := RecSL.O.order( 0, pos.d, pos.p, pos.k );
        pos.orderOS := pos.orderO * (pos.q-1) / 2;
        InfoRecSO4( "#I  order so: ", Runtime()-timer, " msec\n" );
    fi;

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
    InfoRecSO4( "#I  subgroup setup: ", Runtime()-timer, " msec\n" );
    pos.setupTime := Runtime() - timer + pos.setupTime;

end;


#############################################################################
##
#F  RecognizeSO0( <slpos>, <pos>, <tries> ) . . . . . . . . regonize SO0(d,q)
##
RecognizeSO0 := function( slpos, pos, tries )
    local   ord,  poo,  o,  cpol,  po,  A,  mpol,  oe,  pso,  p,  oo,  
            so;

    # if <q> is even <d> must be even (otherwise the group acts reducible)
    if pos.p = 2  then
        Error( "characteristic must be odd,  if dimension if odd" );
    fi;

    # finish set up in case this is not a reentry
    if not IsBound(slpos.isRecSO) or not slpos.isRecSO  then
        pos.operations := RecSO0;
        pos.operations.Setup( pos, slpos );
    fi;

    # use old orders to rule out possibilities
    pos.tries := 0;
    InfoRecSO2( "#I  checking old element orders\n" );
    for ord  in slpos.orders  do

        # construct order,  projective order and semi order
        poo := ord[1];  pso := ord[2];  oe  := ord[3];  po := poo * pso;
        o   := po*oe;   oo  := poo*oe;  so  := pso;

        # check possibilities (don't check smaller field)
        pos.operations.CheckAlternating     ( pos,    poo, pso );
        pos.operations.CheckChevalley       ( pos,    poo, pso );
        pos.operations.CheckLargerField     ( pos, po          );
        pos.operations.CheckImprimitive     ( pos, po          );
        pos.operations.CheckTensorPowers    ( pos, po          );
        pos.operations.CheckTensorProducts  ( pos, po          );
    od;

    # set 'containsSO'
    pos.containsSO :=     not pos.isAlternating
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
    while not pos.containsSO and pos.tries < tries  do

        # find a non-trivial random element of <G>
        repeat
            A := pos.operations.Random(pos.group);
        until A <> pos.identity;
        pos.tries := pos.tries + 1;
        InfoRecSO2("#I  trying ", pos.tries, ".th element of <G>\n");

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
            InfoRecSO2( "#I  projective order = ", 
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
        pos.operations.CheckImprimitive     ( pos,       po       );
        pos.operations.CheckSmallerField    ( pos, cpol, po       );
        pos.operations.CheckSporadicGroups  ( pos,       po       );
        pos.operations.CheckTensorPowers    ( pos,       po       );
        pos.operations.CheckTensorProducts  ( pos,       po       );

        # set 'containsSO'
        pos.containsSO :=     not pos.isAlternating
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
