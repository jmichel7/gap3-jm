#############################################################################
##
#A  recsop.g                    GAP library                      Frank Celler
##
#H  @(#)$Id: recsop.g,v 1.1 1997/03/10 13:49:24 gap Exp $
##
#Y  Copyright (C) 1995,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This  file  contains functions  which will  help to recognize irreduzible
##  groups containing O+(d,q).
##
Revision_recsop_g :=
    "@(#)$Id: recsop.g,v 1.1 1997/03/10 13:49:24 gap Exp $";


#############################################################################
##

#V  RecSOp . . . . . . . . . . . . . . . . . .  functions to recognise O+(d,q)
##
RecSOp := Copy(RecSO);


#############################################################################
##

#F  RecSOp.SetAlternating( <pos>, <slpos> ) . . . . . . set alternating group
##
RecSOp.SetAlternating := RecSO0.SetAlternating;


#############################################################################
##
#F  RecSOp.CheckAlternating( <pos>, <pord>, <psod> )  . .  element order test
##
RecSOp.CheckAlternating := RecSO0.CheckAlternating;


#############################################################################
##
#F  RecSOp.SetChevalley( <pos>, <slpos> ) . . . . . possible Chevalley groups
##
RecSOp.SetChevalley := function( pos, slpos )
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

        # get rid of O+
        if x[6] = ChevD  then
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
#F  RecSOp.CheckChevalley( <pos>, <pord>, <psod> )  . . . . . . element order
##
RecSOp.CheckChevalley := RecSO0.CheckChevalley;


#############################################################################
##
#F  RecSOp.SetLargerField( <pos>, <slpos> ) . . . . . . possible larger field
##
RecSOp.SetLargerField := function( pos, slpos )
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
        if pos.d mod 4 = 0  then
            Add( new, pos.d/2*RecSL.GU.uexponent(pos.d/2,pos.p,pos.k) );
        fi;
        for d  in Filtered( DivisorsInt(pos.d), x -> IsPrimeInt(x) )  do
            if 3 < pos.d/d and pos.d/d mod 2 = 0  then
                Add( new, 2*d*RecSL.O.uexponent(+1,pos.d/d,pos.p,d*pos.k) );
            fi;
        od;
        if pos.d*pos.q/2 mod 2 = 1  then
            Add( new, 2*RecSL.O.uexponent(0,pos.d/2,pos.p,2*pos.k) );
        fi;
        pos.expsLarger := new;
    fi;

end;


#############################################################################
##
#F  RecSOp.CheckLargerField( <pos>, <pord> )  . . . . . . check element order
##
RecSOp.CheckLargerField := RecSO0.CheckLargerField;


#############################################################################
##
#F  RecSOp.SetMysteriousPGroup( <pos>, <slpos> )  .  is mysterious p possible
##
##  R < G < N,  where N is the normalizer of an extraspecial r-group R.
##
RecSOp.SetMysteriousPGroup := function( pos, slpos )
    local   r,  m,  e;

    # mysterious p group: dimension must be a prime power
    pos.isMysteriousP := false;
    if pos.p = pos.q and 2 < pos.p and IsPrimePowerInt(pos.d)  then
        r := SmallestRootInt(pos.d);
        m := LogInt( pos.d, r );
        if r = 2  then
            e := 1;
            while pos.p^e mod r <> 1  do e := e + 1;  od;
            if pos.k = e  then
                pos.isMysteriousP := true;
                pos.expMysteriousP:=4*RecSL.O.uexponent(+1,2*m,2,1);
            fi;
        fi;
    fi;
    if not pos.isMysteriousP  then
        InfoRecSO2( "#I  <G> is no mysterious p-group,  dimension ",
                    "is not a prime power\n" );
    fi;
    
end;


#############################################################################
##
#F  RecSOp.CheckMysteriousPGroup( <pos>, <pord> ) . . . . check element order
##
##  <pord> is the correct part of the projective order.
##
RecSOp.CheckMysteriousPGroup := function( pos, pord )

    # mysterious p-group still possible?
    if not pos.isMysteriousP  then return;  fi;
    pos.statistic[RecSO.STAT_PGROUP] := pos.tries;

    # check element order
    pos.isMysteriousP := pos.expMysteriousP mod pord = 0;
    if not pos.isMysteriousP  then
        InfoRecSO2( "#I  <G> is no mysterious p-group,  element order ",
                    "criteria failed\n" );
    fi;
    
end;


#############################################################################
##
#F  RecSOp.SetImprimitive( <pos>, <slpos> ) . . . possible imprimitive groups
##
RecSOp.SetImprimitive := function( pos, slpos )
    local   d;

    # check if we already know it
    if not slpos.isImprimitive   then
        pos.isImprimitive := false;
        return;
    fi;

    # store the various pairs <e>, <f> for orthogonal decomposition
    pos.dimsImprimitive := [];
    for d  in Filtered( DivisorsInt(pos.d), x -> 1 < x )  do
        if pos.d/d mod 2 = 0  then
            Add( pos.dimsImprimitive, [ pos.d/d, +1, d ] );
            if d mod 2 = 0  then
                Add( pos.dimsImprimitive, [ pos.d/d, -1, d ] );
            fi;
        else
            if 2 < pos.p and pos.squareDiscriminat  then
                Add( pos.dimsImprimitive, [ pos.d/d, 0, d ] );
            fi;
        fi;
    od;
    if pos.q = pos.p and 2 < pos.p and pos.squareDiscriminat then
        Add( pos.dimsImprimitive, [ 1, 0, pos.d ] );
    fi;
    if pos.q*pos.d/2 mod 2 = 1 and pos.q mod 4 = 3  then
        Add( pos.dimsImprimitive, [ pos.d/2, 0, 2 ] );
    fi;

    # store exponent of totally singular decomposition
    pos.expImprimitive := RecSL.GL.uexponent(pos.d/2,pos.p,pos.k)*2*(pos.q-1);

    pos.isImprimitive := true;

end;


#############################################################################
##
#F  RecSOp.CheckImprimitive( <pos>, <pord> )  . . . . . . check element order
##
##  Assume that all integer <= <pos.d> are no semi primes.
##
RecSOp.CheckImprimitive := function( pos, pord )
    local   new,  p,  m;

    # imprimitivity impossible?
    if not pos.isImprimitive  then return;  fi;
    pos.statistic[RecSOp.STAT_PRIMITIVE] := pos.tries;

    # check totally singular decomposition
    if pos.expImprimitive <> false  then
        if pos.expImprimitive mod pord <> 0  then
            pos.expImprimitive := false;
        fi;
    fi;

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
    if 0 = Length(new) and pos.expImprimitive = false  then
        pos.isImprimitive := false;
        InfoRecSO2( "#I  <G> is not imprimitive,  element order ",
                    "criteria failed\n" );
    fi;
    
end;


#############################################################################
##
#F  RecSOp.SetSmallerField( <pos>, <slpos> )  . . . .  possible smaller field
##
##  If <q> = <p>^<k> is a prime, then the matrices are already written over a
##  prime  field,  no  reduction  is possible in this case.
##
RecSOp.SetSmallerField := function( pos, slpos )
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
            if IsPrime(pos.k/i)  then
               AddSet(pos.expsSmaller,2*RecSL.O.uexponent(+1,pos.d,pos.p,i));
            fi;
        od;
        if pos.k mod 2 = 0  then
            AddSet( pos.expsSmaller,
                    2*RecSL.O.uexponent(-1,pos.d,pos.p,pos.k/2) );
        fi;
        pos.smallerField := GF(pos.p);
    fi;         

end;        


#############################################################################
##
#F  RecSOp.CheckSmallerField( <pos>, <cpol>, <pord> ) check order and charpol
##
RecSOp.CheckSmallerField := RecSO0.CheckSmallerField;


#############################################################################
##
#F  RecSOp.SetSporadicGroups( <pos>, <slpos> )  . .  possible sporadic groups
##
RecSOp.SetSporadicGroups := RecSO0.SetSporadicGroups;


#############################################################################
##
#F  RecSOp.CheckSporadicGroups( <pos>, <pord> ) . . . . .  element order test
##
RecSOp.CheckSporadicGroups := RecSO0.CheckSporadicGroups;


#############################################################################
##
#F  RecSOp.SetTensorPowers( <pos>, <slpos> )  . . . .  possible tensor powers
##
RecSOp.SetTensorPowers := function( pos, slpos )
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
        
    # store the various pairs <t>, <m>
    else
        pos.isTensorPower    := true;
        pos.dimsTensorPowers := [];
        t := LogInt( pos.d, m );
        for i  in DivisorsInt(t)  do
            if i < t and 2 < pos.p and 4 <= m^i  then
                Add( pos.dimsTensorPowers, [ m^i, -1, t/i ] );
                Add( pos.dimsTensorPowers, [ m^i, +1, t/i ] );
            fi;
            if i < t and 2 = pos.p and i mod 2 = 0 and m mod 2 = 0  then
                Add( pos.dimsTensorPowers, [ m^i, t/i ] );
            fi;
        od;
    fi;

end;  


#############################################################################
##
#F  RecSOp.CheckTensorPowers( <pos>, <ord> )  . . . . . . check element order
##
RecSOp.CheckTensorPowers := function( pos, pord )
    local   new,  p,  m;
    
    # are tensor powers still possible?
    if not pos.isTensorPower  then return;  fi;
    pos.statistic[RecSO.STAT_POWER] := pos.tries;

    # check exponents
    new := [];
    for p  in pos.dimsTensorPowers  do
        if Length(p) = 2  then
            m := pord / Gcd( pord, RecSL.SP.uexponent(p[1],pos.p,pos.k) );
            if 2*Factorial(p[2]) mod m = 0
               and Sum(Collected(Factors(m)),x->x[1]^x[2]) <= p[2]
            then
                Add( new, p );
            fi;
        else
            m := pord / Gcd( pord, 96*pos.expsO[p[2]+2][p[1]] );
            if 2*Factorial(p[2]) mod m = 0
               and Sum(Collected(Factors(m)),x->x[1]^x[2]) <= p[3]
            then
                Add( new, p );
            fi;
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
#F  RecSOp.SetTensorProducts( <pos>, <slpos> )  . .  possible tensor products
##
RecSOp.SetTensorProducts := function( pos, slpos )
    local   i;

    # check if we already that it cannot be a tensor product
    if not slpos.isTensorProduct   then
        pos.isTensorProduct := false;

    # if the dimension is a prime no tensor products are possible
    elif IsPrime(pos.d)  then
        pos.isTensorProduct := false;
        InfoRecSO2( "#I  <G> is no tensor product,  dimension is ",
                    "a prime\n" );
        
    # compute the exponents of the tensor products
    else
        pos.expsTensorProducts := [];
        for i  in DivisorsInt(pos.d)  do
            if 1<i and i<pos.d and i mod 2 = 0 and pos.d/i mod 2 = 0  then
                Add( pos.expsTensorProducts, [ LcmInt(
                  RecSL.SP.uexponent(i,pos.p,pos.k),
                  RecSL.SP.uexponent(pos.d/i,pos.p,pos.k) ), i ] );
            fi;
            if 2 < pos.p and i mod 2 = 0 and pos.d/i mod 2 = 0 then
                if 3 < i and 3*i < pos.d  then
                    Add( pos.expsTensorProducts, [ 8*LcmInt(
                      pos.expsO[1+2][i], pos.expsO[1+2][pos.d/i] ), i ] );
                    Add( pos.expsTensorProducts, [ 8*LcmInt(
                      pos.expsO[1+2][i], pos.expsO[-1+2][pos.d/i] ), i ] );
                    Add( pos.expsTensorProducts, [ 4*LcmInt(
                      pos.expsO[-1+2][i], pos.expsO[-1+2][pos.d/i] ), i ] );
                fi;
            elif i<pos.d and 2<pos.p and i mod 2=0 and pos.d/i mod 2=1  then
                Add( pos.expsTensorProducts, [ LcmInt(
                  pos.expsO[+1+2][i], pos.expsO[0+2][pos.d/i] ), i ] );
            fi;
        od;
        pos.isTensorProduct := 0 < Length(pos.expsTensorProducts);
    fi;
    

end;        


#############################################################################
##
#F  RecSOp.CheckTensorProducts( <pos>, <pord> ) . . . . . check element order
##
RecSOp.CheckTensorProducts := RecSO0.CheckTensorProducts;


#############################################################################
##

#F  RecSOp.Setup( <pos>, <slpos> )  . . . . . . . . . . . . . . finish set up
##
RecSOp.Setup := function( pos, slpos )
    local   timer;

    # compute the order of O( d, q ) for small dimensions
    if pos.d < 10  then
        timer := Runtime();
        pos.orderO  := RecSL.O.order( +1, pos.d, pos.p, pos.k );
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
#F  RecognizeSOp( <slpos>, <pos>, <tries> ) . . . . . . . . regonize SO+(d,q)
##
RecognizeSOp := function( slpos, pos, tries )
    local   ord,  poo,  o,  cpol,  po,  A,  mpol,  oe,  pso,  p,  oo,  
            so;

    # finish set up in case this is not a reentry
    if not IsBound(slpos.isRecSO) or not slpos.isRecSO  then
        pos.operations := RecSOp;
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
        pos.operations.CheckMysteriousPGroup( pos,    poo      );
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
        pos.operations.CheckMysteriousPGroup( pos,       poo      );
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
