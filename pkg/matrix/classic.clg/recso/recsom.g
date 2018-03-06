#############################################################################
##
#A  recsom.g                    GAP library                      Frank Celler
##
#H  @(#)$Id: recsom.g,v 1.1 1997/03/10 13:49:22 gap Exp $
##
#Y  Copyright (C) 1995,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This  file  contains functions  which will  help to recognize irreduzible
##  groups containing O-(d,q).
##
Revision_recsom_g :=
    "@(#)$Id: recsom.g,v 1.1 1997/03/10 13:49:22 gap Exp $";


#############################################################################
##

#V  RecSOm . . . . . . . . . . . . . . . . . .  functions to recognise O-(d,q)
##
RecSOm := Copy(RecSO);


#############################################################################
##

#F  RecSOm.SetAlternating( <pos>, <slpos> ) . . . . . . set alternating group
##
RecSOm.SetAlternating := RecSO0.SetAlternating;


#############################################################################
##
#F  RecSOm.CheckAlternating( <pos>, <pord>, <psod> )  . .  element order test
##
RecSOm.CheckAlternating := RecSO0.CheckAlternating;


#############################################################################
##
#F  RecSOm.SetChevalley( <pos>, <slpos> ) . . . . . possible Chevalley groups
##
RecSOm.SetChevalley := function( pos, slpos )
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

        # get rid of O-
        if x[6] = Chev2D  then
            if 2*x[3] <> pos.d or x[4] <> pos.p or x[5] <> pos.k  then
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
#F  RecSOm.CheckChevalley( <pos>, <pord>, <psod> )  . . . . . . element order
##
RecSOm.CheckChevalley := RecSO0.CheckChevalley;


#############################################################################
##
#F  RecSOm.SetLargerField( <pos>, <slpos> ) . . . . . . possible larger field
##
RecSOm.SetLargerField := function( pos, slpos )
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
        if pos.d/2 mod 2 = 1  then
            Add( new, 8*pos.d/2*RecSL.GU.uexponent(pos.d/2,pos.p,pos.k) );
        fi;
        for d  in Filtered( DivisorsInt(pos.d), x -> 1 < x )  do
            if 2 < pos.d/d and pos.d/d mod 2 = 0 and IsPrimeInt(d)  then
                Add( new, d*RecSL.O.uexponent(-1,pos.d/d,pos.p,d*pos.k) );
            fi;
        od;
        if 9 < pos.d and pos.d*pos.q/2 mod 2 = 1  then
            Add( new, 2*RecSL.O.uexponent(0,pos.d/2,pos.p,2*pos.k) );
        fi;
        pos.expsLarger := new;
    fi;

end;


#############################################################################
##
#F  RecSOm.CheckLargerField( <pos>, <pord> )  . . . . . . check element order
##
RecSOm.CheckLargerField := RecSO0.CheckLargerField;


#############################################################################
##
#F  RecSOm.SetMysteriousPGroup( <pos>, <slpos> )  .  is mysterious p possible
##
##  there is no mysterious p group in O-
##
RecSOm.SetMysteriousPGroup := function( pos, slpos )
    pos.isMysteriousP := false;
end;


#############################################################################
##
#F  RecSOm.SetImprimitive( <pos>, <slpos> ) . . . possible imprimitive groups
##
RecSOm.SetImprimitive := function( pos, slpos )
    local   d;

    # check if we already know it
    if not slpos.isImprimitive   then
        pos.isImprimitive := false;
        return;
    fi;

    # store the various pairs <e>, <f> for orthogonal decomposition
    pos.dimsImprimitive := [];
    for d  in Filtered( DivisorsInt(pos.d), x -> 1 < x and 2*x <= pos.d )  do
        if d*pos.q mod 2 = 1  then
            Add( pos.dimsImprimitive, [ d, 0, pos.d/d ] );
        fi;
        if d mod 2 = 0 and pos.d/d mod 2 = 1  then
            Add( pos.dimsImprimitive, [ d, -1, pos.d/d ] );
        fi;
    od;
    if pos.q = pos.p and 2 < pos.p  then
        Add( pos.dimsImprimitive, [ 1, 0, pos.d ] );
    fi;
    if pos.d/2 mod 2 = 1 and pos.q mod 4 = 1  then
        Add( pos.dimsImprimitive, [ pos.d/2, 0, 2 ] );
    fi;
    pos.isImprimitive := true;

end;


#############################################################################
##
#F  RecSOm.CheckImprimitive( <pos>, <pord> )  . . . . . . check element order
##
##  Assume that all integer <= <pos.d> are no semi primes.
##
RecSOm.CheckImprimitive := function( pos, pord )
    local   new,  p,  m;

    # imprimitivity impossible?
    if not pos.isImprimitive  then return;  fi;
    pos.statistic[RecSOm.STAT_PRIMITIVE] := pos.tries;

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
#F  RecSOm.SetSmallerField( <pos>, <slpos> )  . . . .  possible smaller field
##
##  If <q> = <p>^<k> is a prime, then the matrices are already written over a
##  prime  field,  no  reduction  is possible in this case.
##
RecSOm.SetSmallerField := function( pos, slpos )
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
            if 2*i < pos.k and IsPrime(pos.k/i)  then
               AddSet(pos.expsSmaller,2*RecSL.O.uexponent(-1,pos.d,pos.p,i));
            fi;
        od;
        pos.smallerField := GF(pos.p);
    fi;         

end;        


#############################################################################
##
#F  RecSOm.CheckSmallerField( <pos>, <cpol>, <pord> ) check order and charpol
##
RecSOm.CheckSmallerField := RecSO0.CheckSmallerField;


#############################################################################
##
#F  RecSOm.SetSporadicGroups( <pos>, <slpos> )  . .  possible sporadic groups
##
RecSOm.SetSporadicGroups := RecSO0.SetSporadicGroups;


#############################################################################
##
#F  RecSOm.CheckSporadicGroups( <pos>, <ord> )  . . . . .  element order test
##
RecSOm.CheckSporadicGroups := RecSO0.CheckSporadicGroups;


#############################################################################
##
#F  RecSOm.SetTensorPowers( <pos>, <slpos> )  . . . .  possible tensor powers
##
##  there are no tensor powers in O-
##
RecSOm.SetTensorPowers := function( pos, slpos )
    pos.isTensorPower := false;
end;  


#############################################################################
##
#F  RecSOm.SetTensorProducts( <pos>, <slpos> )  . .  possible tensor products
##
##  G < O0(d1,q) x O-(d2,q) with d=d1*d2, 3 < d2, q*d1 odd
##
RecSOm.SetTensorProducts := function( pos, slpos )
    local   i;

    # check if we already that it cannot be a tensor product
    if not slpos.isTensorProduct   then
        pos.isTensorProduct := false;

    # if the dimension is a prime no tensor products are possible
    elif IsPrime(pos.d)  then
        pos.isTensorProduct := false;
        InfoRecSO2( "#I  <G> is no tensor product,  dimension is ",
                    "a prime\n" );
        
    # there are no tensor products in characteristic two
    elif pos.p = 2  then
        pos.isTensorProduct := false;
        InfoRecSO2( "#I  <G> is no tensor product,  characteristic is 2\n" );

    # compute the exponents of the tensor products
    else
        pos.expsTensorProducts := [];
        for i  in DivisorsInt(pos.d)  do
            if 1 < i and i mod 2 = 1 and 3 < pos.d/i  then
                Add( pos.expsTensorProducts, [ LcmInt(
                    pos.expsO[0+2][i], pos.expsO[-1+2][pos.d/i] ), i ] );
            fi;
        od;
        pos.isTensorProduct := 0 < Length(pos.expsTensorProducts);
    fi;
    

end;        


#############################################################################
##
#F  RecSOm.CheckTensorProducts( <pos>, <pord> ) . . . . . check element order
##
RecSOm.CheckTensorProducts := RecSO0.CheckTensorProducts;


#############################################################################
##

#F  RecSOm.Setup( <pos>, <slpos> )  . . . . . . . . . . . . . . finish set up
##
RecSOm.Setup := function( pos, slpos )
    local   timer;

    # compute the order of O( d, q ) for small dimensions
    if pos.d < 10  then
        timer := Runtime();
        pos.orderO  := RecSL.O.order( -1, pos.d, pos.p, pos.k );
        pos.orderOS := pos.orderO * (pos.q-1) / 2;
        InfoRecSO4( "#I  order o: ", Runtime()-timer, " msec\n" );
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
    pos.operations.SetTensorProducts  ( pos, slpos );
    pos.operations.SetTensorPowers    ( pos, slpos );
    InfoRecSO4( "#I  subgroup setup: ", Runtime()-timer, " msec\n" );
    pos.setupTime := Runtime() - timer + pos.setupTime;

end;


#############################################################################
##
#F  RecognizeSOm( <slpos>, <pos>, <tries> ) . . . . . . . . regonize SOm(d,q)
##
RecognizeSOm := function( slpos, pos, tries )
    local   ord,  poo,  o,  cpol,  po,  A,  mpol,  oe,  pso,  p,  oo,  
            so;

    # finish set up in case this is not a reentry
    if not IsBound(slpos.isRecSO) or not slpos.isRecSO  then
        pos.operations := RecSOm;
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
        pos.operations.CheckLargerField     ( pos,       po     );
        pos.operations.CheckImprimitive     ( pos,       po     );
        pos.operations.CheckSmallerField    ( pos, cpol, po     );
        pos.operations.CheckSporadicGroups  ( pos,       po     );
        pos.operations.CheckTensorProducts  ( pos,       po     );

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
