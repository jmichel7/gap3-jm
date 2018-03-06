###########################################################################
##
#A  autrec.g                 autag package                 Michael J Smith
##
##  November 1996
##
##  This file forms part of a package for computing automorphism groups of
##  finite soluble groups which are given in terms of special soluble group
##  presentations.
##
###########################################################################


## Changes:

# 21:10 Mon 29 Jul 1996: powers of auts are done using binary expansion. 

# 13:55 Sun 28 Jul 1996: made a domain for auts so that the aut group
# record is a proper gap group. Also removed dependance on GrpHomByIms
# except for inversion --- they are no longer carried around in the aut
# records.


IsAut := function (a)
    return IsRec(a) and IsBound(a.isAut) and a.isAut;
end;
AutOps.IsAut := IsAut;


# domain of all automorphism records
#
AutGroupElements.isDomain := true;
AutGroupElements.name := "AutGroupElements";

AutGroupElements.isFinite := false;
AutGroupElements.size := "infinity";

AutGroupElements.operations := Copy( GroupElementsOps );
AutGroupElementsOps := AutGroupElements.operations;


# function that makes an automorphism group
#
AutGroupElementsOps.Group := function ( AutGroupElements, gens, id )
    local A;
    if not ForAll( gens, gen -> gen.group = gens[1].group ) then
        Error("the automorphism group generators must act on the same group");
    fi;

    # let the default function do the main work
    A := GroupElementsOps.Group( AutGroupElements, gens, id );
    
    A.isAutGroup := true;
    A.operations := AutGroupOps;

    # now store some information about what the group acts on etc
    A.group := id.group;
    
    return A;
end;



AutGroupOps.\in := function (a, A)
    #
    # NB to avoid having to compute a full list of elements in order to
    # define a subgroup, we check whether we know <A> is all of the
    # aut group of the group, and if so just return yes. Otherwise, we
    # have to check the elements of <A>.

    if a.group <> A.group then
        return false;
    fi;
    if IsBound(A.isFullAutGroup) and A.isFullAutGroup then
        return true;
    fi;
    return a in Elements(A);
end;



AutOps.MakeAut := function (arg)
    #
    # Make an automorphism record
    #
    #     Aut( <G>, <ims> [, <w>]) -->  aut specified by images of gen set
    #     Aut( <G>, <g> [, <w>])   -->  inner aut given by <g> in <G>
    #
    # The optional argument <w> is aut weight, which indicates where on
    # the LG-series-factors the aut first starts acting.

    local G, images, r;
    G := arg[1];
    images := arg[2];

    # make the group element
    r := rec();
    r.isGroupElement := true;
    r.isAut := true;
    r.domain := AutGroupElements;
    r.operations := AutOps;

    # now the aut specific information
    if IsList(images) then
        r.images := images;
    else
        r.images := List(G.generators, g -> g^images);
        if r.images = G.generators then
            r.inner := G.identity;
        else
            r.inner := images;
        fi;
    fi;
    r.group := G;
    
    # the weight (indicating which layers of LG-series it acts on)
    if Length(arg) > 2 then
        r.weight := arg[3];
    fi;
    return r;
end;


Aut := function (G, images)
    #
    # Make and aut record of <G> defined by <images>, and check whether it
    # actually defines an automorphism.
    local a;
    a := AutOps.MakeAut(G, images);
    if not IsAutomorphism(AutOps.hom(a, false)) then
        Error("images do not define an automorphism");
    fi;
    return a;
end;


AutOps.MakeInnerAut := function (G, g)
    if not g in G then
        Error("<g> not an element of <G> in InnerAut\n");
    fi;
    return AutOps.MakeAut(G, g);
end;
InnerAut := AutOps.MakeInnerAut;


AutOps.Print := function (r)
    if IsBound(r.inner) and r.inner <> false then
        Print("InnerAut(",r.group,", ", r.inner,")");
    else
        Print("Aut(",r.group,", ",r.images,")");
    fi;
end;


AutOps.\= := function (r, s)
    return r.images = s.images;
end;


AutOps.\< := function (r, s)
    return r.images < s.images;
end;


AutOps.\* := function (a, b)
    local prd;
    if IsAut(a) and IsAut(b) then
        if AutOps.IsKnownInner(a) and AutOps.IsKnownInner(b) then
            prd := AutOps.MakeAut(a.group, a.inner*b.inner);
        else
            prd := AutOps.MakeAut(a.group, AutOps.Compose(a, b));
            if IsBound(a.inverse) and IsBound(b.inverse) then
                prd.inverse := AutOps.MakeAut(a.group,
                                     AutOps.Compose(b.inverse,a.inverse));
                prd.inverse.inverse := prd;
            fi;
        fi;
    elif IsList(a) then
        prd := List( a, x -> x * b);
    elif IsList(b) then
        prd := List( b, x -> a * x);
    else
        Error("product of <a> and <b> is not defined");
    fi;

    return prd;
end;


AutOps.\^ := function (a,b)
    local r;
    if AutOps.IsAut(a) and IsInt(b) then
        if AutOps.IsKnownInner(a) then
            r := AutOps.MakeAut(a.group, a.inner^b);
        else
            if b = 0 then
                r := AutOps.trivialAut(a.group);
            elif b > 0 then
                r := AutOps.Power(a, b);
            elif b < 0 then
                r := AutOps.Power(AutOps.inverse(a),(-b));
            fi;
        fi;
    elif AutOps.IsAut(a) and AutOps.IsAut(b) then
        r :=  b^-1 * a * b;
    elif IsAgWord(a) and AutOps.IsAut(b) and (a in b.group) then
        r := AutGroupOps.Map(b,a);
    elif IsList(a) then
        r := List(a, x -> x^b);
    elif AutOps.IsAut(b) and IsSubgroup(b.group, a) then
        r := Subgroup(a.parent, a.generators^b);
    else
        Error(" power of <lft> by <rgt> must be defined");
    fi;

    return r;
end;


AutOps.Power := function (a, n)
    # 
    # the following powering algorithm is directly from Knuth (4.6.3)

    local y, z, odd;

    y := AutOps.trivialAut(a.group);
    z := a;
    repeat
        odd := n mod 2 = 1;
        n := Int(n/2);
        if odd then
            y := y * z;
            if n <> 0 then
                z := z * z;
            fi;
        else
            z := z * z;
        fi;
    until n = 0;
    return  y;
end;


AutOps.hom := function (arg)
    # 
    #     AutOps.hom(a)    
    #     AutOps.hom(a, true)       <- assume images define valid aut
    #     AutOps.hom(a, false)      <- don't assume it
    # 
    # Convert aut record into a GroupHomomorphismByImages record.  If
    # second argument <flag> is true (the default), set the fields of the
    # result to tell GAP that it is actually a valid automorphism (and save
    # some runtime later when it needs to know the answer to that question).
    
    local a, G, flag, hom;

    a := arg[1];
    G := a.group;
    if Length(arg) > 1 then
        flag := arg[2];
    else
        flag := true;
    fi;
    hom := GroupHomomorphismByImages(G, G, G.generators, a.images);
    if flag and not AutGroupOps.Debugging then
        # we set all of these fields because we wish to avoid corresponding
        # checks during the computation.
        hom.isMapping       := true;
        hom.isInjective     := true;
        hom.isSurjective    := true;
        hom.isBijection     := true;
        hom.isHomomorphism  := true;
        hom.isMonomorphism  := true;
        hom.isEpimorphism   := true;
        hom.isIsomorphism   := true;
        hom.isEndomorphism  := true;
        hom.isAutomorphism  := true;
    fi;
    return hom;
end;    


AutOps.inverse := function (a)
    local hom, t;
    if not IsBound(a.inverse) then
        if AutOps.IsKnownInner(a) then
            a.inverse := AutOps.MakeAut(a.group, a.inner^-1);
            a.inverse.inverse := a;
        else
            hom := AutOps.hom(a)^-1;
            t := AutOps.MakeAut(a.group, hom.genimages);
            a.inverse := t;
            t.inverse := a;
        fi;
    fi;
    return a.inverse;
end;


AutOps.trivialAut := function (G)
    return AutOps.MakeAut(G, G.identity);
end;


AutOps.Compose := function (aut, aut2)
    #
    # Have an automorphism of an AgGroup, given in terms of a list of
    # images of the generators, and wish to compose it with some
    # homomorphism whose domain is the AgGroup.
    #
    # NB this is substantially faster than composing GrHomByIm's in GAP3r4p3

    local images, n, ans, i, exp, j;

    if AutOps.IsAut(aut2) then
        images := aut2.images;
    else
        images := aut2;
    fi;
    n := Length(aut.images);
    ans := [1..n];
    for i in ans do
        exp := ExponentsAgWord(aut.images[i]);
        ans[i] := images[1]^exp[1];
        for j in [2..n] do
            if exp[j] <> 0 then
                ans[i] := ans[i]*images[j]^exp[j];
            fi;
        od;
    od;
    return ans;
end;


AutOps.IsNontrivial := function (a)
    return not a.images = a.group.generators;
end;

AutOps.IsTrivial := function (a)
    return a.images = a.group.generators;
end;


AutOps.IsInner := function (a)
    # 
    # A proactive function: if not already determined, check whether or
    # not <a> is inner. Returns true if <a> is inner, false otherwise.
    # Also sets the <a.inner> to reflect this new knowledge.

    if IsBound(a.inner) then
        if a.inner <> false then
            return true;
        else
            return false;
        fi;
    else
        a.inner := AutGroupOps.CheckInner(a.group, a.images);
        return a.inner <> false;
    fi;
end;
IsInnerAut := AutOps.IsInner;


AutOps.IsKnownInner := function (a)
    #
    # Return true if <a> is known to be an inner aut. Returns false if we
    # know <a> is not inner or if we haven't checked yet and don't know.

    return IsBound(a.inner) and a.inner <> false;
end;


AutOps.IsOuter := function (a)
    return not AutOps.IsInner(a);
end;


AutOps.IsKnownOuter := function (a)
    return IsBound(a.inner) and a.inner = false;
end;


AutOps.MappedAut := function (G, a, H) 
    #
    # If G and H are two copies of a group, so that an isomorphism maps the
    # generating set of G directly to that of H, then this function takes
    # an aut <a> of <G>, and convert it to an aut of H. It can also be
    # called with <G> and <H> generating sets instead of groups.

    local b;
    if AutOps.IsKnownInner(a) then
        b := AutOps.MakeAut(H, AutGroupOps.Map(H, a.inner));
    else
        b := AutOps.MakeAut(H, List(a.images, g -> 
                     AutGroupOps.Map(H, g)));
    fi;
    if IsBound(a.weight) then 
        b.weight := a.weight;
    fi;
    return b;
end;



AutGroupOps.ComputeCentraliserChain := function (G)
    #
    # compute a chain of centralisers in <G> for determining whether or not
    # an aut is inner.
    
    local rt, gens, n, perm, index, C, i;

    InfoAutgroup("#I Computing centraliser chain");
    rt := Runtime();

    gens := G.generators;
    n := Length(gens);

    # Unresolved: whether to compute centraliser chain from group generator
    # 1 up to n, or in opposite order (or indeed any other order).
    # 
    index := [n, n-1 .. 1];
    #index := [1..n];

    C := [G];
    for i in [1 .. n] do
        C[i+1] := Centralizer(C[i], gens[index[i]]);
    od;

    InfoAutgroup(" ( ", Time(Runtime()-rt), ")\n");

    G.centraliserChain := C;
    G.centraliserIndex := index;
end;


AutGroupOps.CheckInner := function (G, aut)
    #
    # <aut> is the image of the generating sequence of <G> under an
    # automorphism.  check whether or not the automorphism is inner.
    #
    # Returns false if not inner, and an element inducing the inner
    # automorphism if it is inner.
    
    local gens, n, index, w, outer, i, b;
    
    gens := G.generators;
    n := Length(gens);

    if not IsBound(G.centraliserChain) then
        AutGroupOps.ComputeCentraliserChain(G);
    fi;
    
    index := G.centraliserIndex;

    w := gens[1]^0;
    outer := false;
    i := 1;

    while i <= Length(gens) and not outer do
        
        b := RepresentativeOperation(G.centraliserChain[i], 
                     gens[index[i]], aut[index[i]]^(w^-1));
        
        outer := b = false;
        if not outer then
            w := b * w;
        fi;
        i := i + 1;
    od;
    
    if outer then
        return false;
    else
        return w;
    fi;
end;


AutGroupOps.Map := function (arg)
    # 
    #     Map(gens, g)
    #     Map(G, g)
    #
    #
    # where <g> is an ag-word. Returns the image of <g> under the
    # action taking each generator of its group to the corresponding
    # element of <gens> (or to <G.generators> when called with a group).

    local gens, exp, ans, j, domain, word, codomain, dgens, cgens, image, 
          i, g, k;

    if Length(arg) = 2 then
        gens := arg[1];
        exp := ExponentsAgWord(arg[2]);
    fi;

    if IsGroup(gens) then
        gens := gens.generators;
    elif AutOps.IsAut(gens) then
        gens := gens.images;
    fi;
    ans := gens[1]^exp[1];
    for j in [2..Minimum(Length(exp),Length(gens))] do
        if exp[j] <> 0 then
            ans := ans * gens[j]^exp[j];
        fi;
    od;
    return ans;

end;
