#############################################################################
##
#A  codeword.g             GUAVA library                       Reinald Baart
#A                                                        &Jasper Cramwinckel
#A                                                           &Erik Roijackers
##
##  This file contains functions for working with codewords
##  Codeword is a record with the following fields:
##  .vector
##  .isCodeword
##  .isDomain
##  .q
##  .weight (not used yet)
##  .wordLength
##  .treatAsPoly
##
#H  $Log: codeword.g,v $
#H  Revision 1.1  1997/01/20 15:15:46  werner
#H  Upgrade from Guava 1.2 to Guava 1.3 for GAP release 3.4.4.
#H
#H  Revision 1.9  1994/10/28  10:42:10  rbaart
#H  WeightCodeword: added zero to local vars
#H
#H  Revision 1.8  1994/10/28  10:24:42  jcramwin
#H  Changed a DistanceCodeword in DistanceVecFFE
#H
#H  Revision 1.7  1994/10/28  09:21:45  jcramwin
#H  WeightCodeword and DistanceCodeword now use DistanceVecFFE
#H
#H  Revision 1.6  1994/10/26  11:02:59  rbaart
#H  DistanceCodeword: fixed bug
#H
#H  Revision 1.5  1994/10/17  10:37:46  rbaart
#H  Changed name of file
#H
#H  Revision 1.4  1994/10/13  15:21:51  jcramwin
#H  changed some names, to make ik more Gappy (if that is at all
#H  possible with the Codeword record
#H
#H  Revision 1.3  1994/10/04  13:33:33  rbaart
#H  Support: fixed bug with .weight
#H
#H  Revision 1.2  1994/09/28  12:55:54  jcramwin
#H  changed intFFE in IntVecFFE
#H
#H  Revision 1.1  1994/09/28  09:51:23  jcramwin
#H  Initial revision
#H
##

#############################################################################
##
#F  IsCodeword( <v> ) . . . . . . . . .  determines if <v> is a codeword type
##
IsCodeword := function(v) 
    return IsRec(v) and IsDomain(v) and IsBound(v.isCodeword) and v.isCodeword;
end;

#############################################################################
##
#F  Codeword( <list> [, <F>] or . . . . . . . . . . . .  creates new codeword
#F  Codeword( <P> [, <n>] [, <F>] ) . . . . . . . . . . . . . . . . . . . . .
##
Codeword := function(arg)
    local v, i, p, l, zero, defFieldSize, elms, done;
    l := Length(arg);
    if l>3 or l<1 then
        Error("Give me a break guys, I am only doing my job here");
    fi;
    p := Copy(arg[1]);
    if l = 2 and IsCode(arg[2]) then
        arg :=[,WordLength(arg[2]),Field(arg[2])];
        l := 3;
    fi;
    if IsList(p) and not ( IsString(p) or IsRat(p[1]) or IsFFE(p[1]) ) then
            if l = 1 then
                return List(p, i->Codeword(i));
            elif l = 2 then
                return List(p,i->Codeword(i,arg[2]));
            else
                return List(p,i->Codeword(i,arg[2], arg[3]));
            fi;
    elif IsCodeword(p) then
        if l > 1 and IsInt(arg[2]) then
            if arg[2] > p.wordLength then
                Append(p.vector,List([p.wordLength..arg[2]-1], i->Z(p.q)*0));
                p.wordLength := arg[2];
            elif arg[2] < p.wordLength then
                p.vector := p.vector{[1..arg[2]]};
                p.wordLength := arg[2];
            fi;
            i := 3;
        else
            i := 2;
        fi;
        if l = i and IsField(arg[l]) and Size(arg[l]) <> p.q  then
            #if Size(arg[l]) mod p.q = 0 then
            #    p.q :=  Size(arg[l]);
            #else
                Error("illegal field conversion");
            #fi;
        fi;
        return p;
    else
        v:= rec(isDomain:=true,isCodeword:=true,operations:=CodewordOps);
        if IsList(p) then
            if IsString(p) then
                v.vector := [];
                for i in p do
                    if i = '1' then
                        Add(v.vector, 1 );
                    elif i = '2' then
                        Add(v.vector, 2 );
                    elif i = '3' then
                        Add(v.vector, 3 );
                    elif i = '4' then
                        Add(v.vector, 4 );
                    elif i = '5' then
                        Add(v.vector, 5 );
                    elif i = '6' then
                        Add(v.vector, 6 );
                    elif i = '7' then
                        Add(v.vector, 7 );
                    elif i = '8' then
                        Add(v.vector, 8 );
                    elif i = '9' then
                        Add(v.vector, 9 );
                    else
                        Add(v.vector, 0 );
                    fi;
                od;
                defFieldSize := false;
            else
                v.vector := p;
                if IsRat(p[1]) then
                    defFieldSize := false;
                elif IsFFE(p[1]) then
                    defFieldSize := CharFFE(p)^DegreeFFE(p);
                fi;
            fi;

            done := false;
            if IsField(arg[l]) then
                v.q := Size(arg[l]);
                if IsInt(defFieldSize) then
                    if defFieldSize^LogInt(v.q, defFieldSize) = v.q then
                        done := true;
                    else
                        Error("illegal field conversion");
                    fi;
                fi;
            else
                if defFieldSize = false then
                    v.q := Maximum(v.vector) + 1;
                    while not IsPrimePowerInt( v.q ) do
                        v.q := v.q + 1;
                    od;
                else
                    v.q := defFieldSize;
                    done := true;
                fi;
            fi;
            if not done then
                if IsPrime(v.q) then
                    v.vector := v.vector*Z(v.q)^0;
                else
                    elms := Elements(GF(v.q));
                    for i in [1..Length(v.vector)] do
                        v.vector[i] := elms [ Int(v.vector[i]) mod v.q + 1];
                    od;
                fi;
            fi;

            if l > 1 and IsInt(arg[2]) then
                i := Length(v.vector);
                if arg[2] > i then
                    Append(v.vector,List([i+1..arg[2]],i->0*Z(v.q)));
                elif arg[2] < i then
                    v.vector := v.vector{[1..arg[2]]};
                fi;
            fi;
            v.treatAsPoly := false;
            v.wordLength := Length(v.vector);
        elif IsPolynomial(p) then
            if l > 1 and IsInt(arg[2]) then
                v.wordLength := arg[2];
                i := 3;
            else
                v.wordLength := Degree(p)+1;
                i := 2;
            fi;
            if l = i and arg[l] <> p.baseRing then
                Error("illegal field conversion");
            fi;
            v.q := Size(p.baseRing);
            v.treatAsPoly := true;
            v.vector := List([1..p.valuation], x->0*Z(v.q));
            Append(v.vector,p.coefficients);
            Append(v.vector, List([1..v.wordLength-Length(v.vector)],
                    x->0*Z(v.q)));
        else
            Error("I don't know how to convert this into a codeword type");
        fi;
        return v;
    fi;
end;

#############################################################################
##
#F  VectorCodeword( <arg> ) . . . . . . . . . . . .  converts input to vector
##
## Input may be codeword, polynomial, vector or a list of those
##
VectorCodeword := function(arg)
    local v;
    if Length(arg) = 1 then
    	if not IsCodeword(arg[1]) then
        	v := Codeword(arg[1]);
        else
        	v := arg[1];
        fi;
    elif Length(arg) = 2 then
        v := Codeword(arg[1], arg[2]);
    else
        v := Codeword(arg[1], arg[2], arg[3]);
    fi;
    if IsList(v) then
        return List(v, i -> VectorCodeword(i));
    else
        return v.vector;
    fi;
end;

#############################################################################
##
#F  PolyCodeword( <arg> ) . . . . . . . . . . converts input to polynomial(s)
##
## Input may be codeword, polynomial, vector or a list of those
##
PolyCodeword := function(arg)
    local v;
    if Length(arg) = 1 then
    	if not IsCodeword(arg[1]) then
        	v := Codeword(arg[1]);
        else
        	v := arg[1];
        fi;
    elif Length(arg) = 2 then
        v := Codeword(arg[1], arg[2]);
    else
        v := Codeword(arg[1], arg[2], arg[3]);
    fi;
    if IsList(v) then
        return List(v, i -> PolyCodeword(i));
    else
        return Polynomial(GF(v.q),v.vector);
    fi;
end;

#############################################################################
##
#F  DistanceCodeword( <a>, <b> )  . the distance between codeword <a> and <b>
##
DistanceCodeword := function( a, b )
    return DistanceVecFFE(a.vector, b.vector);
end;

#############################################################################
##
#F  WeightCodeword( <v> ) . . . . . . . calculates the weight of codeword <v>
##
WeightCodeword := function(v)
    local i, zero;
    if not IsBound(v.weight) then
        zero := 0*Z(v.q);
        v.weight := DistanceVecFFE( List([1..v.wordLength], i->zero),
                            v.vector);
    fi;
    return v.weight;
end;

#############################################################################
##
#F  Support( <v> )  . . . . . . . set of coordinates in which <v> is not zero
##
Support := function (v)
    local i, S, zero;
    S:=[];
    zero := 0 * Z(v.q);
    v.weight := 0;
    for i in [1..v.wordLength] do
        if v.vector[i] <> zero then
            Add(S, i);
            v.weight := v.weight + 1;
        fi;
    od;
    return S;
end;

#############################################################################
##
#F  TreatAsVector( <v> )  . . . . . . . . . . . .  treat codeword as a vector
##
##  The codeword <v> will be treated as a vector
##
TreatAsVector := function(v)
    local i;
    if IsList(v) then
        for i in v do
            i.treatAsPoly := false;
        od;
    else
        v.treatAsPoly := false;
    fi;
end;

#############################################################################
##
#F  TreatAsPoly( <v> )  . . . . . . . . . . . .  treat codeword as polynomial
##
##  The codeword <v> will be treated as a polynomial
##
TreatAsPoly := function(v)
    local i;
    if IsList(v) then
        for i in v do
            if IsCodeword(i) then
                i.treatAsPoly := true;
            fi;
        od;
    else
        if IsCodeword(v) then
            v.treatAsPoly := true;
        fi;
    fi;
end;
        
#############################################################################
##
#F  CodewordOps . . . . . . . . . . . . . . . . . .  operations for codewords
##
CodewordOps := rec( name := "CodewordOps" );

#############################################################################
##
#F  CodewordOps.Print( <v> )  . . . . . . . . . . . . . . . prints a codeword
##
CodewordOps.Print := function(v)
    local l, isclear, power, i;
    if v.treatAsPoly then
        isclear := true;
        for power in Reversed([0..v.wordLength-1]) do
            if v.vector[power + 1] <> 0*Z(v.q) then
                if not isclear then
                    Print(" + ");
                fi;
                isclear := false;
                if power = 0 or v.vector[power + 1] <> Z(v.q)^0 then
                    if IsPrime(v.q) then
                        Print(String(Int( v.vector[power + 1] ) ));
                    else
                        i := LogFFE( v.vector[power + 1], Z(v.q));
                        if i = 0 then
                            Print("1");
                        elif i = 1 then
                            Print("a");
                        else
                            Print("(a^",String(i),")");
                        fi;  
                    fi;
                fi;
                if power > 0 then
                    Print("x");
                    if power > 1 then
                        Print("^", String( power ));
                    fi;
                fi;
            fi;
        od;
        if isclear then
            Print("0");
        fi;
    else
        if not IsPrime(v.q) then
            Print("[ ");
            for i in v.vector do
                if i = 0*Z(v.q) then
                    Print("0 ");
                else
                    l := LogFFE(i, Z(v.q));
                    if l = 0 then
                        Print("1 ");
                    elif l = 1 then
                        Print("a ");
                    else
                        Print("a^",String(l)," ");
                    fi;
                fi;
            od;
            Print("]");
        else
            Print("[ ");
            for i in IntVecFFE(v.vector) do
                Print(i," ");
            od;
            Print("]");
        fi;
    fi;
end;

#############################################################################
##
#F  CodewordOps.\+( <l>, <r> )  . . . . . . . . . . . . . .  sum of codewords
##
CodewordOps.\+ := function(l,r)
    local left;
    if not IsCodeword(r)  then
        return r + l;
    elif IsCode(l) then
    	return CosetCode(l, r);
    elif IsFFE(l) or IsRat(l) then
        return Codeword( l + r.vector);
    else
        return Codeword(VectorCodeword(l) + r.vector);
    fi;
end;

#############################################################################
##
#F  CodewordOps.\*( <l>, <r> )  . . . . . . . . . . . .  product of codewords
##
CodewordOps.\* := function(l,r)
    local left;
    if not IsCodeword(r)  then
        return r * l;
    elif IsFFE(l) or IsRat(l) or IsMat(l) then
        return Codeword( l * r.vector);
    elif IsCodeword(l) then
    	return l.vector * r.vector;
    else
    	Error("<l> has incompatible type");
    fi;
end;

#############################################################################
##
#F  CodewordOps.\-( <l>, <r> )  . . . . . . . . . . . difference of codewords
##
CodewordOps.\- := function(l,r)
    if not IsCodeword(r)  then
        return -(r - l);
    elif IsFFE(l) or IsRat(l) then
        return Codeword( l - r.vector);
    else
        return Codeword(VectorCodeword(l) - r.vector);
    fi;
end;

#############################################################################
##
#F  CodewordOps.\=( <l>, <r> )  . . . . . . . . . . . . equality of codewords
##
CodewordOps.\= := function(l,r)
    if not IsCodeword(r) then
        return r = l;
    else
        return VectorCodeword(l) = r.vector;
    fi;
end;

#############################################################################
##
#F  CodewordOps.Field( <v> )  . . . . . . . . . . . . . . . field of codeword
##
CodewordOps.Field := function (v)
    return GF(v.q);
end;

#############################################################################
##
#F  NullWord( <C> ) or NullWord( <n>, <F> ) . . . . . . . . . . all zero word
##
NullWord := function(arg)
    local n, F;
    if Length(arg) = 1 and IsCode(arg[1]) then
        n := WordLength(arg[1]);
        F := Field(arg[1]);
    elif Length(arg) = 1 then
        n := arg[1];
        F := GF(2);
    elif Length(arg) = 2 then
        n := arg[1];
        F := arg[2];
    else
        Error("usage: NullWord( n [, F] )");
    fi;
    return Codeword(List([1..n], i-> 0), F);
end;
