###########################################################################
##
#A  split.g                  autag package                 Michael J Smith
##
##  November 1996
##
##  This file forms part of a package for computing automorphism groups of
##  finite soluble groups which are given in terms of special soluble group
##  presentations.
##
###########################################################################



AutGroupOps.LiftSplit := function (E, auts)
    #
    # Take an extension and the auts of the acting group return a list of
    # generators for the aut group of the extension
    #
    # The subgroup of automorphisms that lift to the extension group is
    # computed by finding the stabiliser of the module under the action of
    # automorphisms of the factor group. If $M$ is the module with
    # associated representation $\xi$, then the automorphism $\alpha$ act
    # on $M$ to give the module with associated representation
    # $\xi'=\alpha^{-1} \xi$.

    local n, d, p, rt1, rt, C, x, aut, y, r, D, cgens, cmat, inners, size;

    n := E.n;  d := E.d;
    p := E.module.field.size;

    GModOps.ComputeDecomp(E.module);

    rt1 := Runtime();

    InfoAutgroup("#I LiftSplit: find compatible pairs\n");
    rt := Runtime();

    C := AutGroupOps.CompatiblePairs(E, auts, true);

    InfoAutgroup("#I LiftSplit: compatible pairs found");
    InfoAutgroup(" (", Time(Runtime()-rt), ")\n");

    # check compatible pairs --- debugging
    if AutGroupOps.Debugging then
        if false in List(C, c -> PairOps.IsCompatible(E, c.kappa, c.nu)) then
            Error("\n\nNB: problem with compatible pairs\n\n");
        fi;
    fi;

    # Take generators of the subgroup of compatible pairs, and return
    # generators for the automorphism group of the split extension

    auts := [];

    # first the automorphisms from compatible pairs.
    for x in C do
        # images of the first n generators
        if IsBound(x.inner) then
            aut := AutOps.MakeAut(E.extension,
                           AutGroupOps.Transversal(E, x.inner), 
                           x.weight);
            AutGroupOps.Add(auts, aut);
        else
            aut := List(x.kappa.images, g -> AutGroupOps.Transversal(E,g)); 
            y := List(x.nu, i -> List(i, j -> Int(j)));
            for r in y do           # for each row of the companion matrix
                Add(aut, AutGroupOps.vec2ext(E, r));
            od;
            aut := AutOps.MakeAut(E.extension, aut, x.weight);
            AutGroupOps.Add(auts, aut);
        fi;
    od;

    # All automorphisms associated with derivations are inner. Since
    # we are effectively working modulo inner automorphisms, we need
    # only adjust the order of the automorphism group appropriately.
    
    InfoAutgroup("#I order: (Derivations) changing from ", E.autOrder);
    D := Centralizer(E.subgroup, E.extension);
    cgens := D.generators;
    cmat := List(cgens, c -> ExponentsAgWord(c,n+1,n+d,Z(p)));
    TriangulizeMat(cmat);
    inners := List(Difference([1..d], List(cmat, 
                      x -> PositionProperty(x, IsNonZeroElt))),
                   i -> AutOps.MakeAut(E.extension, 
                           E.extension.generators[n+i], 
                           AutGroupOps.AutWeight(E.layer, 2)));
    Append(auts, inners);

    size := Size(E.subgroup) / Size(D);
    E.autOrder := E.autOrder * size;
    E.autChain := E.autChain * size;
    E.autChain[AutGroupOps.AutWeight(E.layer,2)] := size;
    InfoAutgroup(" to ",E.autOrder,"\n");

    
    InfoAutgroup("#I LiftSplit: total time = ", Time(Runtime()-rt1), "\n");

    return auts;

end;



AutGroupOps.ActOnModule := function (M, a)
    #
    # <M> a G-module, <a> in Aut(G).  Return the module <M>^<a>.  It is
    # the module whose representation is g->g^{<a>^-1*<M>.rep}

    local a1, N, i, indec;

    a1 := a^-1;

    # Make a copy of the module M
    #
    N := ShallowCopy(M);

    # Now change those components that are acted upon by the automorphism
    N.genimages := List(a1.images,
                        g -> AutGroupOps.Map(N.genimages, g));
    N.generators := N.genimages{N.groupgens};
    N.matrices := N.generators;
    # N.rep := AutOps.hom(a1)*M.rep;

    # Note that the decomposition we have computed remains valid. We just
    # need to adjust the generators in each of them.
    if IsBound(M.indecomposables) then
        N.indecomposables := [];
        for i in [1..Length(M.indecomposables)] do
            indec := ShallowCopy(M.indecomposables[i]);
            indec.genimages := List(a1.images,
                                    g -> AutGroupOps.Map(indec.genimages, g));
            indec.generators := indec.genimages{N.groupgens};
            indec.matrices := indec.generators;
            N.indecomposables[i] := indec;
        od;
    fi;

    return N;
end;



AutGroupOps.CompatiblePairs := function (E, auts, flag)
    #
    # Compute the generating set for the group of compatible pairs.
    #
    # NB: if flag=false then do not compute module automorphisms
    # This is useful in the nonsplit lifting section, where no liftings
    # of the identity act on the module, so they are not needed.

    local fnAutModule, ans, g, pair, n, d, F, G, gens, M, Compat, rt, 
          trivaut, orbstab, autags, triv, x;
    
    fnAutModule := rec();

    fnAutModule.action := AutGroupOps.ActOnModule;

    fnAutModule.equivalent :=  function (stabiliser, new, trans, Mnew, Morb, 
                                       weight)
        # <new> and <trans> lie in the same coset of the stabiliser if and
        # only if the two modules <Mnew> and <Morb> are isomorphic.  When
        # they are, the isomorphism is the mate for the Schreier generator
        # (to make compatible).
                              
        local ans, g, pair;
        ans := GModOps.IsomModules(Mnew, Morb);  

        if not ans = false then
            g := new*trans^-1;
            # marry the G-aut and module-aut
            pair := PairOps.MakePair(E, g, ans);
            pair.weight := weight;
            AutGroupOps.Add(stabiliser, pair);
        fi;
        return not ans = false;
    end;
    
    n := E.n;  d := E.d;
    F := E.module.field;

    G := E.group;
    gens := G.generators;

    M := E.module;

    Compat := [];
    if Length(auts) > 0 then

        # find the stabilizer in A of the isomorphism type of M
        rt := Runtime();
        trivaut := AutOps.MakeAut(E.group, E.group.identity,
                           AutGroupOps.AutWeight(E.layer, 2));
        orbstab := OrbitOps.OrbitStabiliser(auts, trivaut,
                           M, fnAutModule);
        InfoAutgroup("#I Found compatible pairs stabilizer");
        InfoAutgroup("  (", Time(Runtime()-rt), ")\n");

        AutGroupOps.AdjustChain(E, orbstab);
        InfoAutgroup ("#I order: changing from ", E.autOrder);
        E.autOrder := E.autOrder / Length(orbstab.transversal);
        InfoAutgroup (" to ",E.autOrder,"\n");

        # extract the compatible pairs
        Compat := orbstab.stabiliser;

    fi;

    # return early if asked not to compute module automorphisms - 
    # i.e. compatible pairs with trivial action on the group and
    # non-trivial action on the module.
    #
    if flag = false then
        return Compat;
    fi;

    # now the module automorphisms that commute with action of the group
    if d = 1 and Size(F)=2 then
        # trivial general linear group
        Ignore();
    else
        rt := Runtime();
        autags := GModOps.ModuleAutGens(M);

        InfoAutgroup("#I Found generating set for module automorphisms");
        InfoAutgroup(" (", Time(Runtime()-rt), ")\n");
        
        InfoAutgroup ("#I order: changing from ", E.autOrder);
        E.autOrder := E.autOrder * GModOps.ModuleAutSize(M);
        InfoAutgroup (" to ",E.autOrder,"\n");

        # now blow up the orders of the aut subgroups by this latest
        # factor...
        E.autChain := E.autChain * GModOps.ModuleAutSize(M);

        # ...and add its order at the bottom of the list.
        E.autChain[AutGroupOps.AutWeight(E.layer,1)] := 
          GModOps.ModuleAutSize(M);

        triv := AutOps.MakeAut(E.group, E.group.identity);

        for x in autags do
            pair := PairOps.MakePair(E, triv, x);
            pair.weight := AutGroupOps.AutWeight(E.layer, 1);
            AutGroupOps.Add(Compat, pair);
        od;

    fi;

    return Compat;

end;


