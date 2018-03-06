#############################################################################
##
#A  recsu.g                     GAP library                      Frank Celler
##
#H  @(#)$Id: recsu.g,v 1.1 1997/03/10 13:49:35 gap Exp $
##
#Y  Copyright (C) 1995,   Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This  file  contains functions  which will  help to recognize irreduzible
##  groups containing GU(n,q).
##
Revision_recsu_g :=
    "@(#)$Id: recsu.g,v 1.1 1997/03/10 13:49:35 gap Exp $";


#############################################################################
##
#F  InfoRecSU?(...) . . . . . . . . . . . . . . . . . . . . . . .  infomation
##
##  InfoRecSU4: runtime information
##  InfoRecSU5: MeatAxe
##
if not IsBound(InfoRecSU1)   then InfoRecSU1  := Ignore;  fi;
if not IsBound(InfoRecSU2)   then InfoRecSU2  := Ignore;  fi;
if not IsBound(InfoRecSU3)   then InfoRecSU3  := Ignore;  fi;
if not IsBound(InfoRecSU4)   then InfoRecSU4  := Ignore;  fi;
if not IsBound(InfoRecSU5)   then InfoRecSU5  := Ignore;  fi;


#############################################################################
##

#V  RecSU . . . . . . . . . . . . . . . . . .  functions to recognise SU(n,q)
##
RecSU := rec();
RecSU.STAT_ALT       :=  1;
RecSU.STAT_CHEV      :=  2;
RecSU.STAT_LARGER    :=  3;
RecSU.STAT_PGROUP    :=  4;
RecSU.STAT_POWER     :=  5;
RecSU.STAT_PRIMITIVE :=  6;
RecSU.STAT_PRODUCT   :=  7;
RecSU.STAT_SMALLER   :=  8;
RecSU.STAT_SPORADIC  :=  9;
RecSU.STAT_FORM      := 10;


#############################################################################
##
#F  RecSU.Random( <G> ) . . . . . . . . . . .  return a random element of <G>
##
RecSU.Random := RecSL.Random;


#############################################################################
##
#F  RecSU.UnitaryForm( <slpos>, <pos> ) . . . . .  try to find a unitary form
##
RecSU.UnitaryForm := function( slpos, pos )
    if not IsBound(slpos.unitaryForm)  then
        RecSL.InvariantForm(slpos);
    fi;
    if IsBound(slpos.unitaryForm)  then
        pos.unitaryForm     := slpos.unitaryForm;
        pos.gScalars        := slpos.unitaryScalars;
        pos.isMeataxeLarger := false;
        return pos.unitaryForm;
    else
        return false;
    fi;
end;


#############################################################################
##

#F  RecSU.SetAlternating( <pos>, <slpos> )  . . . . . . set alternating group
##
RecSU.SetAlternating := function( pos, slpos )
    local   timer,  n,  o,  x;

    # use alternating groups still possible
    if not slpos.isAlternating  then
        pos.alternating   := [];
        pos.isAlternating := false;
        return;
    fi;

    # check group order
    if IsBound(pos.orderGU)  then
        pos.alternating := Filtered( slpos.alternating, x -> 2*pos.orderGUS
                                     mod Factorial(x) = 0 );
    else
        pos.alternating := slpos.alternating;
    fi;
    pos.isAlternating := 0 < Length(pos.alternating);

end;


#############################################################################
##
#F  RecSU.CheckAlternating( <pos>, <pord>, <psod> ) . . .  element order test
##
RecSU.CheckAlternating := function( pos, pord, psod )
    local   z,  new,  n,  z2;

    # are alternating groups still possible?
    if not pos.isAlternating then return;  fi;
    pos.statistic[RecSU.STAT_ALT] := pos.tries;

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
        InfoRecSU2( "#I  <G> is not an alternating group,  ",
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
        InfoRecSU2( "#I  <G> is not an alternating group,  ",
                    "element order criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSU.SetChevalley( <pos>, <slpos> )  . . . . . possible Chevalley groups
##
RecSU.SetChevalley := function( pos, slpos )
    local   new,  x;

    # use Chevalley groups still possible
    if not slpos.isChevalley  then
        pos.expsChev    := [];
        pos.isChevalley := false;
        return;
    fi;

    # check orders
    if IsBound(pos.orderGU)  then
        new := Filtered( slpos.expsChev, x -> pos.orderGUS mod
                                  x[6].order(x[3],x[4],x[5]) = 0 );
    else
        new := slpos.expsChev;
    fi;

    # avoid certain small cases
    pos.expsChev := [];
    for x  in new  do

        # avoid 2A_<d>-1(<qq>)
        if x[6] = Chev2A  then
            if x[3]+1 <> pos.d or x[4] <> pos.p or 2*x[5] <> pos.k  then
                Add( pos.expsChev, x );
            fi;

        # avoid U4(2) = PSp4(3) = C2(3)
        elif x[6] = ChevC and pos.d = 4 and pos.q = 4 then
            if x[3] <> 2 or x[4] <> 3 or x[5] <> 1  then
                Add( pos.expsChev, x );
            fi;

        # set tighter bound for PSL(3,4) < SU(4,3)
        elif x[6] = ChevA and pos.d = 4 and pos.q = 9  then
            if x[3] = 2 and x[4] = 2 and x[5] = 2  then
                x := ShallowCopy(x);
                x[1] := 1680;
            fi;
            Add( pos.expsChev, x );


        # set element orders for PSL(2,7) < SU(3,3)
        elif x[6] = ChevA and pos.d = 3 and pos.q = 9  then
            if x[3] = 1 and x[4] = 7 and x[5] = 1  then
                x := ShallowCopy(x);
                x[1] := [ 1, 2, 3, 4, 7 ];
            fi;
            Add( pos.expsChev, x );

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
#F  RecSU.CheckChevalley( <pos>, <pord>, <psod> ) . . . . . . . element order
##
RecSU.CheckChevalley := function( pos, pord, psod )
    local   new,  x;

    # are Chevalley groups still possible?
    if not pos.isChevalley  then return;  fi;
    pos.statistic[RecSU.STAT_CHEV] := pos.tries;

    # filtered chevalley groups
    new := [];
    for x  in pos.expsChev  do
        if IsInt(x[1])  then
            if x[1] mod pord = 0 and ( 1 = psod or 1 < GcdInt(x,psod) )  then
                Add( new, x );
            fi;
        elif pord in x[1]  then
            Add( new, x );
        fi;
    od;
    pos.expsChev := new;

    if 0 = Length(pos.expsChev)  then
        pos.isChevalley := false;
        InfoRecSU2( "#I  <G> is not a Chevalley group,  ",
                    "element order criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSU.SetLargerField( <pos>, <slpos> )  . . . . . . possible larger field
##
##  GU(d/r,qq^r).r for r > 2
##
RecSU.SetLargerField := function( pos, slpos )
    local   d;

    # check if we already know it
    if not slpos.isLarger   then
        pos.isLarger := false;
        return;
    fi;

    # no reduction is possible if dimension is one
    pos.isLarger := 1 < pos.d;
    if not pos.isLarger  then
        InfoRecSU2("#I  <G> is not definable over a larger field,  ",
                   "dimension is one\n" );

    # representation that preserve a larger field
    else
        d:=Filtered( DivisorsInt(pos.d), x -> 1 < x and x mod 2 = 1 );
        d:=List(d,x->(pos.qq-1)*RecSL.GU.uexponent(pos.d/x,pos.p,x*pos.k)*x);
        pos.expsLarger := d;
    fi;

end;


#############################################################################
##
#F  RecSU.CheckLargerField( <pos>, <pord> ) . . . . . . . check element order
##
RecSU.CheckLargerField := function( pos, pord )
    local   c,  I,  g,  d,  e;

    # are larger fields still possible?
    if not pos.isLarger  then return;  fi;
    pos.statistic[RecSU.STAT_LARGER] := pos.tries;

    # check element orders
    pos.expsLarger := Filtered(pos.expsLarger, x -> x mod pord = 0);
    if 0 = Length(pos.expsLarger)  then
        pos.isLarger := false;
        InfoRecSU2( "#I  <G> is not definable over a larger field,  ",
                    "element order criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSU.SetMysteriousPGroup( <pos>, <slpos> ) . .  is mysterious p possible
##
##  R < G < N,  where N is the normalizer of an extraspecial r-group R.
##
RecSU.SetMysteriousPGroup := function( pos, slpos )
    local   r,  m,  e;

    # mysterious p group: dimension must be a prime power
    pos.isMysteriousP := false;
    if pos.k mod 2 = 0 and IsPrimePowerInt(pos.d)  then
        r := SmallestRootInt(pos.d);
        m := LogInt( pos.d, r );
        if pos.p <> r  then
            if r <> 2  then
                e := 1;
                while pos.p^e mod r <> 1  do e := e + 1;  od;
                if pos.k = e and r mod 2 = 1  then
                    if pos.d = 3 and pos.q = 4  then
                        pos.isMysteriousP := false;
                    else
                        pos.isMysteriousP := true;
                        pos.expMysteriousP:=r^2*RecSL.SP.uexponent(2*m,r,1);
                    fi;
                fi;
            else
                e := 1;
                while pos.p^e mod 4 <> 1  do e := e + 1;  od;
                if pos.k = e and e = 2  then
                    pos.isMysteriousP  := true;
                    pos.expMysteriousP := 4*RecSL.SP.uexponent(2*m,2,1);
                fi;
            fi;
        fi;
    fi;
    if not pos.isMysteriousP  then
        InfoRecSU2( "#I  <G> is no mysterious p-group,  dimension ",
                    "is not a prime power\n" );
    fi;
    
end;


#############################################################################
##
#F  RecSU.CheckMysteriousPGroup( <pos>, <pord> )  . . . . check element order
##
##  <pord> is the correct part of the projective order.
##
RecSU.CheckMysteriousPGroup := function( pos, pord )

    # mysterious p-group still possible?
    if not pos.isMysteriousP  then return;  fi;
    pos.statistic[RecSU.STAT_PGROUP] := pos.tries;

    # check element order
    pos.isMysteriousP := pos.expMysteriousP mod pord = 0;
    if not pos.isMysteriousP  then
        InfoRecSU2( "#I  <G> is no mysterious p-group,  element order ",
                    "criteria failed\n" );
    fi;
    
end;


#############################################################################
##
#F  RecSU.SetImprimitive( <pos>, <slpos> )  . . . possible imprimitive groups
##
##   d = e*f,  G < GU(e,q) wr Sym(f) or
##  or
##   GL(d/2,q).2
##
RecSU.SetImprimitive := function( pos, slpos )
    local   d;

    # check if we already know it
    if not slpos.isImprimitive   then
        pos.isImprimitive := false;
        return;
    fi;

    # store the various pairs <e>, <f> for unitary decomposition
    d := DivisorsInt(pos.d);
    d := Filtered( d, x -> x <> pos.d );
    pos.dimsImprimitive := List( d, x -> [ x, pos.d/x ] );

    # store exponent of totally singular decomposition
    if pos.d mod 2 = 0  then
        pos.expImprimitive := RecSL.GL.uexponent(pos.d/2,pos.p,pos.k)*2;
    else
        pos.expImprimitive := false;
    fi;
    pos.isImprimitive  := true;

end;


#############################################################################
##
#F  RecSU.CheckImprimitive( <pos>, <pord> ) . . . . . . . check element order
##
##  Assume that all integer <= <pos.d> are no semi primes.
##
RecSU.CheckImprimitive := function( pos, pord )
    local   new,  p,  m;

    # imprimitivity impossible?
    if not pos.isImprimitive  then return;  fi;
    pos.statistic[RecSU.STAT_PRIMITIVE] := pos.tries;

    # check totally singular decomposition
    if pos.expImprimitive <> false  then
        if pos.expImprimitive mod pord <> 0  then
            pos.expImprimitive := false;
        fi;
    fi;
    
    # check exponents of unitary decomposition
    new := [];
    for p  in pos.dimsImprimitive  do
        m := pord / Gcd( pord, pos.expsGU[p[1]] );
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
        InfoRecSU2( "#I  <G> is not imprimitive,  element order ",
                    "criteria failed\n" );
    fi;
    
end;


#############################################################################
##
#F  RecSU.SetSmallerField( <pos>, <slpos> ) . . . . .  possible smaller field
##
##  If <q> = <p>^<k> is a prime, then the matrices are already written over a
##  prime  field,  no  reduction  is possible in this case.   Otherwise upper
##  bounds for the exponents  of PSU(<n>,<p>^<i>),  <i> T <k>,  are computed.
##  The record component <pos>.smallerField containes the smaller field still
##  possible.
##
RecSU.SetSmallerField := function( pos, slpos )
    local   i,  j,  q,  exp;

    # check if we already know it
    if not slpos.isSmaller   then
        pos.isSmaller:= false;
        return;
    fi;

    # if <field> is the prime field,  no reduction is possible
    pos.isSmaller := 1 < pos.k;
    if not pos.isSmaller  then
        InfoRecSU2( "#I  <G> is not definable over a smaller field,  ",
                    "field is the prime field\n" );

    # loop over the maximal divisors of <k>
    else
        pos.expsSmaller := [];
        for i in DivisorsInt(pos.k) do
            if 2 < pos.k/i and IsPrime(pos.k/i)  then
                AddSet( pos.expsSmaller,RecSL.GU.uexponent(pos.d,pos.p,i));
            fi;
        od;
        if pos.q*pos.d mod 2 = 1  then
            AddSet(pos.expsSmaller,RecSL.O.uexponent(0,pos.d,pos.q,pos.k));
        fi;
        if pos.q mod 2 = 1 and pos.d mod 2 = 0  then
            AddSet(pos.expsSmaller,RecSL.O.uexponent(-1,pos.d,pos.q,pos.k));
            AddSet(pos.expsSmaller,RecSL.O.uexponent(+1,pos.d,pos.q,pos.k));
        fi;
        if pos.d mod 2 = 0  then
            AddSet( pos.expsSmaller,RecSL.SP.uexponent(pos.d,pos.q,pos.k) );
        fi;
        pos.smallerField := GF(pos.p);
    fi;         

end;        


#############################################################################
##
#F  RecSU.CheckSmallerField( <pos>, <cpol>, <pord> )  check order and charpol
##
##  <pord> is the projective order  of the  group  element  A, <cpol>  is the
##  characteristic polynomial  of A.  The function  first  checks  if  <pord>
##  divides  at least one  of the exponents stored in  <pos>.expsSmaller.  It
##  then computes the smallest field which contains <cpol> * zeta.
##
RecSU.CheckSmallerField := function( pos, cpol, pord )
    local   c,  d,  t,  p,  i0,  I;

    # are smaller fields still possible?
    if not pos.isSmaller  then return;  fi;
    pos.statistic[RecSU.STAT_SMALLER] := pos.tries;

    # first check the projective order
    pos.expsSmaller := Filtered( pos.expsSmaller, x -> x mod pord = 0 );
    if 0 = Length(pos.expsSmaller)  then
        pos.isSmaller := false;
        InfoRecSU2( "#I  <G> is not definable over a smaller field,  ",
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
            InfoRecSU2( "#I  <G> is not definable over a smaller field, ",
                        "char polynomial failed\n" );
        fi;
    fi;

end;


#############################################################################
##
#F  RecSU.SetSporadicGroups( <pos>, <slpos> ) . . .  possible Sporadic groups
##
RecSU.SetSporadicGroups := function( pos, slpos )

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
#F  RecSU.CheckSporadicGroups( <pos>, <pord> )  . . . . .  element order test
##
RecSU.CheckSporadicGroups := function( pos, pord )
    
    # Sporadic groups still possible?
    if not pos.isSporadic  then return;  fi;
    pos.statistic[RecSU.STAT_SPORADIC] := pos.tries;

    # check the possible element orders of the automorphism groups
    pos.sporadicGroups := Filtered( pos.sporadicGroups, x -> pord in
                                    SporadicGroupsInfo.orders[x] );

    # if the list is trivial reset '<pos>.isSporadic'
    if 0 = Length(pos.sporadicGroups)  then
        pos.isSporadic := false;
        InfoRecSU2( "#I  <G> is not a Sporadic group,  ",
                    "element order criteria failed\n" );
    fi;

    # if the exponent is too big stop
    if pos.exponent >= pos.q3d  then
        pos.isSporadic := false;
        InfoRecSU2( "#I  <G> is not a Sporadic group, ",
                    "group exponent criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSU.SetTensorPowers( <pos>, <slpos> ) . . . . .  possible tensor powers
##
##  d = m^t, G < GU(m,q) wr Sym(t).   (q,m) != (2,3)  and m >= 2. Assume that
##  there are no semi primes smaller then t+1.
##
RecSU.SetTensorPowers := function( pos, slpos )
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
        InfoRecSU2( "#I  <G> is no tensor power,  dimension is ",
                    "not a power\n" );
        
    # store the various pairs <t>, <m>
    else
        pos.isTensorPower    := true;
        pos.dimsTensorPowers := [];
        t := LogInt( pos.d, m );
        for i  in DivisorsInt(t)  do
            if i < t and 2 < m^i and not ( pos.q = 2 and m^i = 3 )  then
                Add( pos.dimsTensorPowers, [ m^i, t/i ] );
            fi;
        od;
    fi;

end;  


#############################################################################
##
#F  RecSU.CheckTensorPowers( <pos>, <ord> ) . . . . . . . check element order
##
##  <ord> is the order of the  group element A. If A lies in GU(e,q) X Sym(f)
##  then gcd(<ord>,exp(GL(e,q))) must be a valid element order in Sym(f).
##        
##  Assume that all integer <= <pos.d> are no semi primes.
##
RecSU.CheckTensorPowers := function( pos, ord )
    local   new,  p,  m;
    
    # are tensor powers still possible?
    if not pos.isTensorPower  then return;  fi;
    pos.statistic[RecSU.STAT_POWER] := pos.tries;

    # check exponents
    new := [];
    for p  in pos.dimsTensorPowers  do
        m := ord / Gcd( ord, pos.expsGU[p[1]] );
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
        InfoRecSU2( "#I  <G> is no tensor power,  element order ",
                    "criteria failed\n" );
    fi;

end;


#############################################################################
##
#F  RecSU.SetTensorProducts( <pos>, <slpos> ) . . .  possible tensor products
##
##  G < GU(d1,q) x GU(d2,q) with d=d1*d2
##
RecSU.SetTensorProducts := function( pos, slpos )
    local   i,  s,  exp;

    # check if we already know it
    if not slpos.isTensorProduct   then
        pos.isTensorProduct := false;
        return;
    fi;

    # if <pos.d> is a prime, no tensor products are possible
    if IsPrime(pos.d)  then
        pos.isTensorProduct := false;
        InfoRecSU2( "#I  <G> is no tensor product,  dimension is ",
                    "a prime\n" );
        return;
    else
        pos.isTensorProduct := true;
    fi;
        
    # compute the exponents of the tensor products
    pos.expsTensorProducts := [];
    for i  in DivisorsInt(pos.d)  do
        if 1 < i and i*i < pos.d  then
            exp := LcmInt( pos.expsGU[i], pos.expsGU[pos.d/i] );
            Add( pos.expsTensorProducts, [ exp, i ] );
        fi;
    od;

end;        


#############################################################################
##
#F  RecSU.CheckTensorProducts( <pos>, <pord> )  . . . . . check element order
##
RecSU.CheckTensorProducts := function( pos, pord )

    # are tensor products still possible?
    if not pos.isTensorProduct  then return;  fi;
    pos.statistic[RecSU.STAT_PRODUCT] := pos.tries;

    # check if <ord> divides any exponent of a tensor product
    pos.expsTensorProducts := Filtered( pos.expsTensorProducts,
                                        x -> x[1] mod pord = 0 );
    
    # if <pos>.expsTensorProducts is empty no tps are possible
    if 0 = Length(pos.expsTensorProducts)  then
        pos.isTensorProduct := false;
        InfoRecSU2( "#I  <G> is no tensor product,  element order ",
                    "criteria failed\n" );
    fi;

end;


#############################################################################
##

#F  RecSU.Print . . . . . . . . . . . . . . . . . . . . . . nice pretty print
##
RecSU.Print := function( obj )
    local   name;

    if obj.printLevel = 0  then
        Print( "<< SU recognition record >>" );
        return;
    fi;
    Print( "#I  field: ", obj.q, ", dimension: ", obj.d,
           ", number of generators: ", Length(obj.generators), "\n" );

    if IsBound(obj.unitaryForm)  then
        Print( "#I  unitary form is known\n" );
    else
        Print( "#I  no unitary form is known\n" );
        Print( "<< SU recognition record >>" );
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
    
    # Sporadic groups
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
    if IsBound(obj.containsSU) and obj.containsSU  then
        Print( "#I  <G> contains SU( ",obj.d,", ",obj.qq," )\n" );
    else
        Print( "#I  <G> could still be an unitary group\n" );
    fi;
    Print( "<< SU recognition record >>" );
end;


#############################################################################
##
#F  RecSU.SetPrintLevel( <obj>, <lev> ) . . . . . . . . . . . set print level
##
RecSU.SetPrintLevel := function( obj, lev )
    if 2 < lev or lev < 0  then
        Error( "invalid print level" );
    fi;
    obj.printLevel := lev;
end;


#############################################################################
##
#F  RecSU.Setup( <slpos> )  . . . . . . . . . . . . . . . . set up pos record
##
RecSU.Setup := function( slpos )
    local   pos,  i,  j,  k,  timer,  qi,  q2;

    # set up a SU possibility record
    pos            := rec();
    pos.q          := slpos.q;
    pos.p          := slpos.p;
    pos.k          := slpos.k;
    pos.d          := slpos.d;
    pos.q3d        := slpos.q3d;
    pos.qq         := slpos.qq;
    pos.group      := slpos.group;
    pos.identity   := slpos.identity;
    pos.field      := slpos.field;
    pos.exponent   := slpos.exponent;
    pos.generators := slpos.generators;
    pos.orders     := Copy(slpos.orders);
    pos.statistic  := [1..1]*0;
    pos.printLevel := 1;
    pos.containsSU := false;
    pos.isRecSU    := true;
    pos.operations := RecSU;
    pos.setupTime  := -Runtime();
    if IsBound(slpos.unitaryForm)  then
        pos.unitaryForm := slpos.unitaryForm;
    fi;
    InfoRecSU1( "#I  field: ", pos.q, ", dimension: ", pos.d, 
                ", number of generators: ", Length(pos.generators), "\n" );

    # store <slpos>
    pos.recognizeSL := slpos;
        
    # compute the exponent of GU(i,q) for 1 <= i < d, i|d
    timer  := Runtime();
    pos.expsGU := [ pos.qq+1 ];
    j  := pos.qq+1;
    qi := pos.qq;
    i  := 2;
    while pos.d mod i <> 0  do i := i + 1;  od;
    for i  in [ 2 .. pos.d/i ]  do
        qi := qi * pos.qq;
        j  := LcmInt( j, qi - (-1)^i );
        if pos.d mod i = 0  then
            pos.expsGU[i] := j * pos.p^(1+LogInt(i-1, pos.p));
        fi;
    od;
    InfoRecSU4( "#I  exponents: ", Runtime()-timer, " msec\n" );

    # compute the order of GU( d, q ) for small dimensions
    if pos.d < 10  then
        timer := Runtime();
        pos.orderGU  := RecSL.GU.order( pos.d, pos.p, pos.k/2 );
        pos.orderGUS := pos.orderGU * (pos.qq-1);
        InfoRecSU4( "#I  order SU: ", Runtime()-timer, " msec\n" );
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
    InfoRecSU4( "#I  subgroup setup: ", Runtime()-timer, " msec\n" );
    pos.setupTime := Runtime() + pos.setupTime;

    # and return
    return pos;

end;


#############################################################################
##
#F  RecognizeSU( <slpos>, <tries> ) . . . . . . . . . . . .  regonize SU(n,q)
##
RecognizeSU := function( slpos, tries )
    local   pos,  form,  ord,  poo,  o,  A,  mpol,  po,  oe,  pso,  p,  
            oo,  so,  cpol;

    # is this a reentry
    if not IsBound(slpos.isRecSU) or not slpos.isRecSU  then

        # dimension should be at least 3 (otherwise SL=SU)
        if slpos.d = 2  then
            Error( "dimension must be at least 3,  otherwise SU=SL" );
        fi;

        # <G> must have a chance to be the unitary group
        if not slpos.isUnitary  then
            Error( "SL unitary test for the group already failed" );
        fi;

        # field must be a square
        if slpos.k mod 2 = 1  then
            Error( "field must be a square" );
        fi;

        # set up a SU possibility record
        pos := RecSU.Setup(slpos);

        # use old orders to rule out possibilities
        pos.tries := 0;
        InfoRecSU2( "#I  checking old element orders\n" );
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

    # try to find an invariant unitary form
    if not IsBound(pos.unitaryForm)  then
        form := pos.operations.UnitaryForm( slpos, pos );

        # return if we failed to find a unitary form
        if form = false  then
            InfoRecSU1( "#I  failed to find a unitary form\n" );
            return pos;
        else
            slpos.unitaryForm := form;
            InfoRecSU1( "#I  found invariant unitary form\n" );
        fi;
    fi;

    # set 'containsSU'
    pos.containsSU :=     not pos.isAlternating
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
    while not pos.containsSU and pos.tries < tries  do

        # find a non-trivial random element of <G>
        repeat
            A := pos.operations.Random(pos.group);
        until A <> pos.identity;
        pos.tries := pos.tries + 1;
        InfoRecSU2("#I  trying ", pos.tries, ".th element of <G>\n");

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
            InfoRecSU2( "#I  projective order = ", 
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

        # set 'containsSU'
        pos.containsSU :=     not pos.isAlternating
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

RecogniseSU := RecognizeSU;
