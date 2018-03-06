###########################################################################
##
#A  autgrp.g                 autag package                 Michael J Smith
##
##  November 1996
##
##  This file forms part of a package for computing automorphism groups of
##  finite soluble groups which are given in terms of special soluble group
##  presentations.
##
###########################################################################


# This is the main file for the autag package. The algorithm used in these
# files is an improvement of the algorithm described in my PhD thesis. My
# thesis is available by anonymous ftp from:
#
#    ftp://wwwmaths.anu.edu.au/pub/smith/thesis.dvi.gz
#
# However, the complete algorithm will be described in a forthcoming paper
# about this work.


## changes:

# 22:41 Mon 29 Jul 1996: extracted the code to compute a description of a
# single factor in the normal series of an aut group, into function
# AutGroupOps.AutGroupFactor. This can be used to prune unnecessary
# generators in the main loop...


ReadLocalDir(LOCALNAMEaut, "autrec");
ReadLocalDir(LOCALNAMEaut, "extend");
ReadLocalDir(LOCALNAMEaut, "nonsplit");
ReadLocalDir(LOCALNAMEaut, "orbit");
ReadLocalDir(LOCALNAMEaut, "split");
ReadLocalDir(LOCALNAMEaut, "pairs");
ReadLocalDir(LOCALNAMEaut, "eqns");


AutGroupOps.Debugging := false;


AutGroupSagGroup := function (arg)
    #
    # Usage: AutomorphismGroup(G)
    #        AutomorphismGroup(G, l)
    #
    # Compute the automorphism group of a finite soluble group that is
    # given by a special AG-presentation. If optional argument <l> is
    # given, compute the automorphism group of the quotient of <G> by the
    # <l+1>-st term of its LG-series (ie a quotient of <G> which includes
    # the <l>-th term).

    local rt, G, lim, order, nlayers, n, E, A, p, gl, innerauts, i, j, inn, 
          E2, ans;

    rt := Runtime();

    G := arg[1];
    if not IsSagGroup(G) then
        Error("G is not a special ag group\n");
    fi;

    lim := false;
    if Length(arg) > 1 then lim := arg[2]; fi;

    order := 1;
    nlayers := Maximum(G.layers);
    n := Number(G.layers, x -> x = 1);

    if nlayers > 1 then
        E := AutGroupOps.MakeExtension(G, 2);
    else
        E := rec(group := G, extension := G);
    fi;

    # A will hold aut group generators (minus inner auts)
    A := [];
    
    # Compute matrices acting on the top layer of the LG-series.
    p := G.weights[1][3];
    gl := GL(n,p);
    A := List(gl.generators, x -> List(x, y -> Product([1..n], 
                 i -> E.group.generators[i]^Int(y[i]))));
    order := gl.size;

    
    # Convert the matrices acting on the top layer to aut records
    # 
    A := List(A, a -> AutOps.MakeAut(E.group, a, 1)); # NB: weight 1
    
    E.autOrder := order;
    E.autChain := [ order, 1 ];

    # now keep a list of inner automorphisms corresponding to various
    # weights for use in the aut chain
    innerauts := [];

    # Start the lifting to the 2-nd layer of the LG-series.
    i := 2;

    while i <= nlayers and (lim = false or i <= lim) do
        
        if G.weights[G.first[i]][2] = 1 then

            InfoAutgroup("\n#I SPLIT EXTENSION");
            InfoAutgroup(",  module dimension ", E.module.dim, 
                    ",  prime ", E.module.field.size,
                    ",  LG-layer ",i,"\n");

            A := AutGroupOps.LiftSplit(E, A);

        else

            InfoAutgroup("\n#I NONSPLIT EXTENSION");
            InfoAutgroup(",  module dimension ", E.module.dim, 
                    ",  prime ", E.module.field.size,
                    ",  LG-layer ",i,"\n");
            A := AutGroupOps.LiftNonsplit(E, A);

        fi;
        order := E.autOrder;


        # Try to trim some redundant generators from the list. We do this
        # by looking at those weights with a matrix group representation,
        # and removing those that act trivially or in the same way as
        # another in this representation.

        A := AutGroupOps.TrimGenerators(E, A);


        # Remove inner automorphisms from list.
        # Keep them so that inner auts with weights can be used later.
        A := AutGroupOps.Set(A);
        for j in [1..Length(A)] do
            if AutOps.IsInner(A[j]) and AutOps.IsNontrivial(A[j]) then
                inn := AutOps.MakeAut(G, AutGroupOps.Map(G, A[j].inner),
                               A[j].weight);
                AutGroupOps.Add(innerauts, inn);
            fi;
        od;
        A := Filtered(A, x -> not AutOps.IsInner(x));

        InfoAutgroup("#I A",i," = ",A,"\n  |A| = ",order,"\n");
        InfoAutgroup("#I Length(A) = ", Length(A), "\n" );
        InfoAutgroup("#I ",Time(Runtime()-rt)," so far\n");

        if AutGroupOps.Debugging then
            #
            # check automorphisms - debugging
            for j in [1..Length(A)] do
                if not IsAutomorphism(AutOps.hom(A[j],false)) then
                    Error("\n\n NB: problem with automorphism ",j,"\n\n");
                fi;
            od;
        fi;

        # Increment pointer to layer and compute new extension record
        # if necessary.
        i := i + 1;
        if i <= nlayers and (lim = false or i <= lim) then
            E2 := AutGroupOps.NextExtension(E);
            for j in [1..Length(A)] do
                A[j] := AutOps.MappedAut(E.extension, A[j], E2.group);
            od;
            E2.autChain := E.autChain;
            E := E2;
            E.autOrder := order;
        fi;
    od;

    # add the inner automorphisms before returning

    if G = E.extension then
        Append(A, innerauts);
    else
        Append(A, List(innerauts, a -> AutOps.MakeAut(E.extension,
                AutGroupOps.Map(E.extension, a.inner), a.weight)));
    fi;

    # sort the generators before returning -- by increasing weight
    Sort(A, function(a,b) return a.weight < b.weight; end);

    InfoAutgroup("#I AutGrp Total time = ",Time(Runtime()-rt),"\n");

    ans := Group(A, AutOps.trivialAut(E.extension));
    ans.size := order;
    ans.time := Runtime()-rt;
    ans.ordersChain := E.autChain;

    if lim = false or lim >= nlayers then
        ans.isFullAutGroup := true;
        G.autGroup := ans;
    fi;

    return ans;
end;


AutGroupOps.PrintAutGroup := function (A)
    Print("AutGroup(", A.group, ", ", A.generators, ")");
end;    


AutGroupOps.Set := function (set)
    #
    # Sort into descending weight, then trim redundant generators. Note
    # that the sorting is crucial, since we need to keep the higher weight
    # version if there is a duplicate.

    local newset, x;
    newset := [];
    Sort(set, function(a,b) return a.weight > b.weight; end);
    for x in set do
        if not x in newset then
            Add(newset, x);
        fi;
    od;
    return newset;
end;


AutGroupOps.Add := Add;


AutGroupOps.DerivationAuts := function (E)
    #
    # Compute derivations for the extension, and turn these into
    # automorphisms of the extension group.

    local dertoaut, ocb, D, auts, x, a;

    dertoaut := function (deriv)
        local aut, i;
        aut := [];
        for i in [1..E.n] do
            Add(aut, E.extension.generators[i] * 
                AutGroupOps.vec2ext(E,deriv{(i-1)*E.d+[1..E.d]}));
        od;
        Append(aut, List([1..E.d], j -> E.extension.generators[E.n+j]));

        return AutOps.MakeAut(E.extension, aut,
                       AutGroupOps.AutWeight(E.layer,2));
    end;

    ocb := OneCocycles(E.extension, E.subgroup);
    D := Base(ocb.oneCocycles);

    auts := [];
    for x in D do
        a := dertoaut(x);
        Add(auts, a);
    od;

    return rec(auts := auts, length := Length(D));
end;



AutGroupOps.AutWeight := function (layer, type)
    #
    # Return an integer representing the weight of an automorphism.  The
    # weight determines where in the LG-series it first starts acting. If
    # G/G_layer+1 is the first quotient of G by a term of the LG-series on
    # which an aut acts nontrivially, then the aut has weight either
    # AutWeight(layer,1) or AutWeight(layer,2) --- the former iff it acts
    # non-trivially on the quotient G_layer/G_layer+1.
    
    # level=1 means nontrivial action on the layer of the LG-series
    #
    # level=2 means nontrivial action only on the factor group by the
    #         layer of the LG-series.
    
    return 2 * (layer - 1) + type;
end;



AutGroupOps.AdjustChain := function (E, orbstab)
    #
    # take an orbit-stabiliser record and use the transversal weights
    # field to adjust the orders of the chain of normal subgroups in
    # the automorphism group.

    local i, len;
    
    InfoAutgroup("#I autChain: from ", E.autChain);

    for i in [1..Length(E.autChain)] do
        len := Number(orbstab.tweight, j -> j >= i);
        E.autChain[i] := E.autChain[i] / len;
    od;

    InfoAutgroup(" to ", E.autChain, "\n");
    
end;



AutGroupOps.AutGroupFactor := function (A, i)
    #
    # Given aut-group <A>, and weight <i>, return a description of the
    # weight <i> factor in the normal series of <A>. This will either be
    # simply a list [p,e] indicating an elementary abelian factor of
    # size p^e, or a list of matrices generating a matrix group isomorphic
    # to the factor.

    local G, chain, layer, first, p, size, gens, mats, last, dim, a, mat, 
          ret;

    G := A.group;
    chain := A.ordersChain;

    layer := Int((i-1)/2) + 1;  # invert AutOps.Weight !
    
    first := G.first[layer];
    p := G.weights[first][3];

    size := chain[i];
    if i < Length(chain) then size := size / chain[i+1]; fi;

    if i mod 2 = 1 then
        gens := Filtered(A.generators, a -> a.weight = i);
        mats := [];
        first := G.first[layer];
        last := G.first[layer+1] - 1;
        p := G.weights[first][3];
        dim := G.first[layer+1] - G.first[layer];
        for a in gens do
            mat := List(a.images{[first..last]},
                        x -> ExponentsAgWord(x, first, last, Z(p)));
            Add(mats, mat);
        od;
        ret := rec(matrices := mats, 
                   identity := IdentityMat(dim, GF(p)),
                   dim := dim,
                   field := GF(p),
                   size := size
                   );
    else
        ret := rec(prime := p, 
                   exponent := LogInt(size, p),
                   size := size);
    fi;
    
    return ret;
end;


AutGroupFactors := function (arg)
    #
    # Parse the information in the record returned by AutGroupSagGroup
    # and give details of factors in a normal series for the aut group.
    #
    # AutGroupFactors(A [, flag ])
    #
    # if flag=true, then print the information to screen, otherwise
    # return a list describing all of the factors.

    local A, print, chain, factors, i, factor, size;

    A := arg[1];
    if Length(arg) > 1 then
        print := arg[2];
    else
        print := false;
    fi;

    # this is the list of the orders of subgroups of the aut group
    # in the normal series
    chain := A.ordersChain;

    if print then
        Print("\n Order of full automorphism group is ", A.size);
        if not (IsPrime(A.size) or A.size = 1) then
            Print(" = ", StringPP(A.size)); 
        fi;
        Print("\n");
    else
        factors := [];
    fi;

    # loop over the weights of subgroups
    #
    for i in [1..Length(chain)] do

        # find the description
        factor := AutGroupOps.AutGroupFactor(A, i);
        
        # check the order of the factor at this weight (may be trivial)
        size := factor.size;
        if size <> 1 then
            # 
            # non-trivial factor

            if print then
                Print("\n Factor of size ", size);
                if not IsPrime(size) then 
                    Print(" = ", StringPP(size)); 
                fi;
            fi;

            if i mod 2 = 1 then
                # 
                # we have a matrix group representation of this factor

                if print then
                    Print(" (matrix group, weight ", i, ")\n");
                    if Length(factor.matrices) = 0 then
                        Print(" (no matrices...)\n");
                    else
                        Nice(factor.matrices);
                    fi;
                else
                    Add(factors, Group(factor.matrices, factor.identity));
                fi;

            else
                #
                # we know that this factor is elementary abelian

                if print then
                    Print(" (elementary abelian, weight ", i, ")\n");
                else
                    Add(factors, [factor.prime, factor.exponent]);
                fi;
            fi;
        fi;
    od;
    if print then 
        Print("\n"); 
    else
        return factors;
    fi;
end;


AutGroupStructure := function(A)
    #
    # Display aut-group structure

    AutGroupFactors(A, true);

end;


AutGroupOps.TrimGenerators := function (E, gens)
    #
    # Remove redundant generators from generating set <A> for the aut group
    # of <E.extension>.
    #
    # For each <weight> we have either a matrix group representation for
    # the corresponding quotient (when <weight> is odd), or we know that it
    # is elementary abelian of a given order (when <weight> is even).
    #
    # In the former case, we just make sure that each generator of that
    # <weight> produces a different action in the matrix group
    # representation, and throw away ones that produce duplicates.
    #
    # In the latter case, we currently do nothing.

    local rt, A, newgens, weights, maxweight, w, i, factor, size, id, mats, 
          j;

    rt := Runtime();

    # this is the form that AutGroupFactor expects:
    A := rec(group := E.parent, 
             ordersChain := E.autChain,
             generators := gens);

    # we will accumulate the irredundant generators in <newgens>
    newgens := [];

    # get a list of the generator weights
    weights := List(A.generators, a -> a.weight);
 
    # the maximum weight is given by the length of the aut-chain order list
    maxweight := Length(E.autChain);

    # for each weight, find a list of generator indices for that weight
    w := List([1..maxweight], i -> PositionsProperty(weights, j -> j=i));

    for i in [1 .. Length(E.autChain)] do

        # compute the representation we have for weight <i>
        factor := AutGroupOps.AutGroupFactor(A, i);
        
        # if trivial size factor, we don't do anything (ie throw away any
        # "generators" that we may have for that weight.
        size := factor.size;

        if size > 1 and Length(w[i]) > 0 then

            if i mod 2 = 1 then
                #
                # matrix group representation

                id := factor.identity;
                # <mats> will store matrices without repeats
                mats := [];
                for j in [1..Length(factor.matrices)] do
                    if not factor.matrices[j] in mats 
                       and factor.matrices[j] <> id then
                        Add(newgens, A.generators[w[i][j]]);
                        Add(mats, factor.matrices[j]);
                    fi;
                od;
            else
                #
                # elementary abelian factor, just add the gens

                Append(newgens, A.generators{w[i]});
            fi;
        fi;
    od;

    InfoAutgroup("#I Trimmed generators from ", Length(gens),
            " to ", Length(newgens), " (", Time(Runtime()-rt), ")\n");

    return newgens;

end;


AutGroupSeries := function (A)
    #
    # return a sequence of subgroups of <A> which reflect the subnormal
    # series of <A> (the weighted subgroups).

    local chain, series, i, size, gens, sub;

    chain := A.ordersChain;
    series := [ ];

    # loop over the weights of subgroups
    #
    for i in [1..Length(chain)] do
        size := chain[i];
        if i < Length(chain) then size := size / chain[i+1]; fi;
        if size <> 1 then
            gens := Filtered(A.generators,
                            a -> a.weight >= i);
            sub := Subgroup(A, gens);
            sub.weight := i;
            sub.isFullWeight := true;
            Add(series, sub);
        fi;
    od;
    return series;
end;


AutGroupConverted := function (A)
    #
    # convert the aut group <A> into a GAP group generated by GrpHomByIms
    # records.
    
    local ans;

    ans := Group(List(A.generators, a -> AutOps.hom(a)),
                 AutOps.hom(A.identity));
    ans.size := A.size;

    # Normalize strips off the name of the group:
    if IsBound(A.group.name) then
        ans.identity.source.name := A.group.name;
    fi;
    return ans;
end;
