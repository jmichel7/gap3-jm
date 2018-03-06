#############################################################################
##
#A  util.g                  GUAVA library                       Reinald Baart
#A                                                        &Jasper Cramwinckel
#A                                                           &Erik Roijackers
##
##  This file contains miscellaneous functions
##
#H  $Log: util.g,v $
#H  Revision 1.1  1997/01/20 15:16:11  werner
#H  Upgrade from Guava 1.2 to Guava 1.3 for GAP release 3.4.4.
#H
#H  Revision 1.7  1994/11/03  11:17:28  rbaart
#H  NullVector: fixed double quotes
#H
#H  Revision 1.6  1994/11/03  11:08:39  jcramwin
#H  changed some comments and NullVector
#H
#H  Revision 1.5  1994/10/19  10:11:20  rbaart
#H  Changed some comments
#H
#H  Revision 1.4  1994/10/13  15:20:48  rbaart
#H  Changed codeword functions
#H
#H  Revision 1.3  1994/10/06  15:22:00  rbaart
#H  PrimitiveUnityRoot: changed GAP 3.3 to GAP 3
#H
#H  Revision 1.2  1994/09/28  12:53:44  jcramwin
#H  IntPPFFE removed
#H
#H  Revision 1.1  1994/09/28  09:44:04  jcramwin
#H  Initial revision
#H
##

#############################################################################
##
#F  SphereContent( <n>, <e> [, <F>] ) . . . . . . . . . . .  contents of ball
##
##  SphereContent(n, e [, F]) calculates the contents of a ball of radius e in 
##  the space (GF(q))^n
##
SphereContent := function (arg)
    local res, num, den, i,n,e,q_1;
    if Length(arg) = 2 then
        n:=arg[1];
        e:=arg[2];
        q_1:=1;
    elif Length(arg) = 3 then
        n:=arg[1];
        e:=arg[2];
        if IsInt(arg[3]) then
            q_1:=arg[3] - 1;
        else
            q_1:=Size(arg[3]) - 1;
        fi;
    else
        Error("usage: SphereContent( <n>, <radius> [, <F>] )");
    fi;
    res := 0;
    num := 1;
    den := 1;
    for i in [0..e] do
        res := res + (num * den);
        num := num * q_1;
        den := (den * (n-i)) / (i+1);
    od;
    return res;
end;

#############################################################################
##
#F  Krawtchouk( <k>, <i>, <n> [, <F>] ) . . . . . .  Krwatchouk number K_k(i)
##
##  Krawtchouk(k, i, n [, F]) calculates the Krawtchouk number K_k(i) 
##  over field of size q (or 2), wordlength n.
##  Pre: 0 <= k <= n
##
Krawtchouk := function(arg)
    local k,i,n,q;
    if Length(arg) = 3 then
        k := arg[1];
        i := arg[2];
        n := arg[3];
        q := 1;
    elif Length(arg) = 4 then
        k := arg[1];
        i := arg[2];
        n := arg[3];
        if IsInt(arg[4]) then
            q := arg[4] - 1;
        else
            q := Size(arg[4]) - 1;
        fi;
    else
        Error("usage: Krawtchouk( <k>, <i>, <n> [, <F>] )");
    fi;
    if k > n or k < 0 then
        Error("0 <= k <= n");
    elif not IsPrimePowerInt(q+1) then
        Error("q must be a prime power");
    fi;
    return Sum([0..k],j->Binomial(i,j)*Binomial(n-i,k-j)*(-1)^j*q^(k-j));
end;

#############################################################################
##
#F  PermutedCols( <M>, <P> )  . . . . . . . . . .  permutes columns of matrix
##
PermutedCols := function(M, P)
    if P = () then
        return M;
    else
        return List(M, i -> Permuted(i,P));
    fi;
end;

#############################################################################
##
#F  ReciprocalPolynomial( <p> [, <n>] ) . . . . . .  reciprocal of polynomial
##
ReciprocalPolynomial := function (arg)
    local n, p, cl;
    if Length(arg) = 1 then
        p := arg[1];
        cl := VectorCodeword(p);
    elif Length(arg) = 2 then
        p := arg[1];
        n := arg[2];
        cl := VectorCodeword(p, n+1);
    else
        Error("usage: ReciprocalPolynomial( <polynomial> [, <length>] )");
    fi;
    return Polynomial(p.baseRing, Reversed(cl));
end;

#############################################################################
##
#F  CyclotomicCosets( [<q>, ] <n> ) . . . .  cyclotomic cosets of <q> mod <n>
##
CyclotomicCosets := function (arg)
    local addel, set, res, nrelements, elements, start, q, n;
    if Length(arg) = 1 then
        q := 2;
        n := arg[1];
    elif Length(arg) = 2 then
        q := arg[1];
        n := arg[2];
    else
        Error("usage: CyclotomicCosets( <q>, <n> )");
    fi;
    if Gcd(q,n) <> 1 then
        Error("q and n must be relative primes");
    fi;
    res := [[0]];
    nrelements := 1;
    elements := Set([1..n-1]);
    repeat
        start := elements[1];
        addel := start;
        set := [];
        repeat
            Add(set, addel);
            RemoveSet(elements, addel);
            addel := addel * q mod n;
            nrelements := nrelements + 1;
        until addel = start;
        Add(res, set);
    until nrelements >= n;
    return res;
end;

#############################################################################
##
#F  PrimitiveUnityRoot( [<q>, ] <n> ) . .  primitive n'th power root of unity
##
PrimitiveUnityRoot := function(arg)
    local n,q,qm;
    if Length(arg) = 1 then
        n := arg[1];
        q := 2;
    elif Length(arg) = 2 then
        q := arg[1];
        n := arg[2];
        if not IsInt(q) then
            q := Size(q);
        fi;
    else
        Error("usage: PrimitiveUnityRoot( [ <q> ,] <n> )");
    fi;
    qm := q ^ OrderMod(q,n);
    if not qm in [2..65536] then
    Error("GAP 3 cannot compute in a finite field of size larger than 2^16");
    fi;
    return Z(qm)^((qm - 1) / n);
end;

#############################################################################
##
#F  RemoveFiles( <arglist> )  . . . . . . . .  removes all files in <arglist>
##
##  used for functions which use external programs (like Leons stuff)
##
RemoveFiles := function(arg)
    local f;
    for f in arg do
        Exec(Concatenation("rm -f ",f));
    od;
end;

#############################################################################
##
#F  NullVector( <n> [, <F> ] )  . .  vector consisting of <n> coordinates <o>
##
NullVector := function(arg)
    local i, res, n, zero;
    if Length(arg) = 1 then
        n := arg[1];
        zero := 0;
    elif Length(arg) = 2 then
        n := arg[1];
        zero := arg[2].zero;
    else
        Error("usage: NullVector( <n> [, <F>] )");
    fi;
    res := [];
    for i in [1..n] do
        res[i] := zero;
    od;
    return res;
end;

