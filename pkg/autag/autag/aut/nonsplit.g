###########################################################################
##
#A  nonsplit.g               autag package                 Michael J Smith
##
##  November 1996
##
##  This file forms part of a package for computing automorphism groups of
##  finite soluble groups which are given in terms of special soluble group
##  presentations.
##
###########################################################################

## Other References:

# Derek J. Robinson's paper ``Applications of cohomology to the theory of
# groups'', C.M. Campbell and E.F. Robertson editors, Groups St Andrews
# 1981, LMS Lecture Note Series volume 71, pp 46-80.

## Changes:

# 23:16 Thu 15 Aug 1996 - update and improve code documentation.

# 22:14 Wed 10 Apr 1996 - modified to use relations from Sylow p-subgroup
# in order to vastly speed.

# 17:03 Fri 2 Aug 1996 - keep track of information about which relations
# provide equations in LiftAut. Try to solve the system early if the rank
# hasn't changed in a while.




AutGroupOps.LiftNonsplit := function (E, auts)
    # 
    # Take an extension and the auts of the factor group and return a list
    # of generators for the aut group of the extension
    # 
    # Two orbit-stabiliser calculations are performed. The first stabiliser
    # determines the compatible pairs. The second has the compatible pairs
    # modulo module automorphism pairs acting on the cohomology class of
    # the extension. This is a modification of the action of the compatible
    # pairs on the cohomology as described in Robinson.

    local rt1, n, d, p, a, g, rt, compatpairs, fnRec, s, nukap, nu, newpz, 
          i, elem, cocycle, trivpair, trivpoint, orbstab, inducibles, D;

    rt1 := Runtime();

    n := E.n;  d := E.d;
    p := E.module.field.size;

    # compute definitions for the generators from the module.
    AutGroupOps.ComputeDefinitions(E);

    # for computing cohomology of a Sylow p-subgroup
    E.pgens := PositionsProperty([1..n], i -> E.parent.weights[i][3] = p);

    # define the public Sylow subgroup
    E.sylow := Subgroup(E.group, E.group.generators{E.pgens});

    # Make sure that auts normalise the public Sylow p-subgroup
    for a in auts do
        if false in List(E.pgens, i -> a.images[i] in E.sylow) then
            # adjust with inner automorphism to normalise E.sylow
            g := RepresentativeOperation(E.group, 
                         Subgroup(E.group, a.images{E.pgens}), E.sylow);
            a := a * AutOps.MakeAut(E.group, g);
        fi;
    od;

    InfoAutgroup("#I LiftNonsplit: find subgroup of compatible auts\n");
    rt := Runtime();

    # first find compatible pairs --- but we do not want the module
    # automorphism ones, so pass "false" as third arg
    compatpairs := AutGroupOps.CompatiblePairs(E, auts, false);

    InfoAutgroup("#I LiftNonsplit: subgroup of compatible pairs found");
    InfoAutgroup("  (", Time(Runtime()-rt), ")\n");

    # now compute the inducible pairs
    rt := Runtime();

    # The following orbit-stabiliser calculation evolved out of the action
    # of compatible pairs on cohomology classes of extension as described
    # in Robinson (1981). However, the problem with using this to lift aut
    # groups is that the compatible pair subgroup often contains an
    # extremely large subgroup of liftings of the identity --- and *none*
    # of these will lift to the extension here. Hence, all they do is blow
    # out the size of the transversal, which is extremely damaging to the
    # run-time.

    # Instead (as will be described in a paper soon to be written), I
    # project the compatible pairs onto the component corresponding to
    # Aut(G) in order to get the subgroup K of Aut(G) which consists of all
    # auts which can be in a compatible pair. This subgroup can act
    # directly on the cohomology class of the extension, and its stabiliser
    # is the subgroup of auts that lift to the extension.

    # In order to acheive this, compatible pairs are used but the action on
    # the module is adjusted on the fly in the equivalence routine so that
    # a modified pair that is inducible to the extension is obtained iff
    # the kappa component can lift to the extension.


    # The points being passed around below contain two components --- the
    # field <.g> gives the group element used to obtain the point, and <.z>
    # gives the cocycle for the corresponding extension. The <.g> component
    # is needed so that the coycle routine can always work on the original
    # point <z0> (the cocycle of the extension we start with).  The <.z>
    # components are used to check equivalence of extension (differencing
    # and comparing to the coboundaries stored in E.2cobounds).

    fnRec := rec();

    fnRec.equivalent := function (stabiliser, new, trans, newp, orbp, weight)
        #
        # This is the equivalence check for the cohomological stabiliser.
        # The element <s> is constructed from the two compatible pairs
        # <new> and <trans>. If the aut <s.kappa> of <E.group> lifts to
        # <E.ext>, then the action of the lifting on the module can be
        # deduced by its action of the Sylow p-subgroup in this nilpotent
        # layer --- say, <nukap>. The cocycle corresponding to a
        # "corrected" compatible pair <new> can be deduced, and this
        # compared against the cocycle for the transversal representative.

        local s, nukap, nu, newpz, i;

        s := new * trans^-1;

        # now find the module-aut that must be used if the new point
        # is in the same equivalence class as the orbit point.

        nukap := AutGroupOps.ImpliedAction(E, s.kappa);
        if RankMat(nukap) < d then
            #
            # the action on the module is singluar, indicating that <s>
            # cannot possibly be an element of the stabiliser.
            return false;
        fi;

        # if <new> could be adjusted by a compatible pair with trival G
        # action, so that the result is equivalent to <trans>, then this is
        # the action on the module we need to adjust <new> by.
        nu := new.nu^-1 * nukap * trans.nu;

        # now use this on the cocycle associated with the new point
        newpz := ShallowCopy(newp.z);
        for i in [0..Length(newpz)/d-1] do
            newpz{i*d + [1..d]} := newpz{i*d + [1..d]} * nu;
        od;

        # check whether the adjusted cocycle is equivalent to the given one
        if (newpz - orbp.z) in E.2cobounds.rowspace then
            s.weight := weight;
            s.nu := nukap;
            AutGroupOps.Add(stabiliser, s);
            return true;
        else
            return false;
        fi;
    end;

    fnRec.action := function (point, g)
        local elem, cocycle;
        elem := point.g * g;
        cocycle := AutGroupOps.PairCocycle(E, elem);
        return rec( g := elem, z := cocycle);
    end;

    fnRec.weight := function (pair)
        return pair.weight;
    end;

    if Length(compatpairs) > 0 then

        # first we need the 2-coboundaries of E
        AutGroupOps.Compute2Cobounds(E);

        trivpair := PairOps.TrivialPair(E);
        trivpair.weight := AutGroupOps.AutWeight(E.layer, 2);

        # now compute the 2-cocycle corresponding to our extension
        E.z0 := AutGroupOps.PairCocycle(E, trivpair);

        # now do orbit-stabiliser calculation to find stabiliser of
        # the cohomology class of E under the action of the projection
        # onto Aut(G) of the compatible pairs.
        trivpoint := rec(g := trivpair, z := E.z0);
        orbstab := OrbitOps.OrbitStabiliser(compatpairs, 
                           trivpair, trivpoint, fnRec);
        inducibles := AutGroupOps.Set(orbstab.stabiliser);

        AutGroupOps.AdjustChain(E, orbstab);
        InfoAutgroup("#I changing order from ", E.autOrder);
        E.autOrder := E.autOrder / Length(orbstab.transversal);
        InfoAutgroup (" to ",E.autOrder,"\n");

    else
        inducibles := [];
    fi;

    InfoAutgroup("#I LiftNonsplit: inducible auts found");
    InfoAutgroup("  (", Time(Runtime()-rt), ")\n");

    # now AutGroupOps.LiftAut these inducible pairs to the extension group
    rt := Runtime();

    # now lift the compatible pairs to the extension
    auts := List(inducibles, pair -> AutGroupOps.LiftAut(E, pair));

    InfoAutgroup("#I LiftNonsplit: lifted automorphisms");
    InfoAutgroup("  (", Time(Runtime()-rt), ")\n");


    # there are no automorphisms that act on this layer nontrivially
    # for the first time
    E.autChain[AutGroupOps.AutWeight(E.layer,1)] := 1;


    # now the derivations between G and the module M give rise to aut's
    rt := Runtime();
    D := AutGroupOps.DerivationAuts(E);

    InfoAutgroup("#I Calculated derivations from group to module");
    InfoAutgroup("  (", Time(Runtime()-rt), ")\n");

    InfoAutgroup("#I Derivations, length D = ",D.length,"\n");

    InfoAutgroup("#I order: Derivations, changing from ", E.autOrder);
    E.autOrder := E.autOrder * p^D.length;
    InfoAutgroup(" to ",E.autOrder,"\n");

    E.autChain := E.autChain * p^D.length;
    E.autChain[AutGroupOps.AutWeight(E.layer,2)] := p^D.length;

    Append(auts, D.auts);
    
    InfoAutgroup("#I LiftNonsplit: total time = ", Time(Runtime()-rt1), "\n");

    return auts;
end;



AutGroupOps.ComputeRHS := function (E)
    #
    # Store the exponent vectors of the right-hand sides of relations to be
    # used later in collection.

    local n, d, p, gens, i, q, j;

    n := E.n;  d := E.d;
    p := E.module.field.size;
    gens := E.group.generators;

    # store the exponent vectors of the rhs of relations
    E.rhs := List([1..n], x -> []);
    for i in [1..n] do
        q := E.parent.weights[i][3];
        E.rhs[i][i] := ExponentsAgWord(gens[i]^q);
        for j in [i+1..n] do
            E.rhs[j][i] := ExponentsAgWord(gens[j]*gens[i]);
        od;
    od;
end;



AutGroupOps.PairCocycle := function (E, pair)
    #
    # Takes a compatible pair for E and returns a cocycle corresponding to
    # the action of the pair on the extension.  Note that only the
    # relations of the Sylow p-subgroup are used.
    #
    # The cocycle is represented as a vector of tails that would occur in
    # the pcp of the extension that it determines --- a la Alex Wegner's
    # thesis. This representation of cocycles is not faithful, but the one
    # it induces on cohomology classes is faithful.

    local phi, collect, n, d, p, gens, weights, kappa1, nu, all, i, q, m, 
          k, j;

    if not IsBound(E.rhs) then
        AutGroupOps.ComputeRHS(E);
    fi;

    # This is the action of a pair on a cocycle. This action comes from
    # Robinson (1981).
    phi := function (x,y)
        return AutGroupOps.FactorSet(E, 
                       x^kappa1, y^kappa1) * nu;
    end;

    collect := function(j, i)
        local right, m, u, k, e;
        right := gens[1]^0;
        m := [1..d]*0*Z(p);
        u := E.rhs[j][i];
        for k in [n, n-1 .. 1] do
            for e in [1 .. u[k]] do
                m := m + phi(gens[k], right);
                right := gens[k] * right;
            od;
        od;
        return m;
    end;

    n := E.n;  d := E.d;
    p := E.module.field.size;
    gens := E.group.generators;
    weights := E.parent.weights;
    
    kappa1 := pair.kappa^-1;

    nu := pair.nu;

    all := [];
    
    for i in [1..n] do
        
        # NB: only work in the Sylow p-subgroup
        if i in E.pgens then

            # power relation
            q := weights[i][3];
            m := [1..d]*0*Z(p);
            for k in [0..q-1] do
                m := m + phi(gens[i], gens[i]^k);
            od;
            m := m - collect(i,i);

            Append(all,m);

            # conjugation relations
            for j in [i+1..n] do
            
                if j in E.pgens then
                    m := + phi(gens[j], gens[i]) 
                         - collect(j,i);
                    
                    Append(all, m);
                fi;
            od;
        fi;
    od;
    
    return all;
    
end;





AutGroupOps.Compute2Cobounds := function (E)
    #
    # This function calls the GAP SQ functions to compute the coboundaries
    # for the module associated with E. Note that it uses the Sylow
    # p-subgroup instead of the whole group.

    local n, d, p, zero, F, M, C, b;

    n := Length(E.sylow.generators);    # n is number of Sylow p-gens.
    d := E.d;
    p := Size(E.module.field);
    zero := [1..d*n*(n+1)/2] * E.module.field.zero;

    F := FpGroup(E.sylow);
    M := rec (  isMatGroup := true,
                dimension := d,
                field := E.module.field,
                generators := E.module.genimages{E.pgens} );
    C := CollectorSQ(F, M);
    b := AutGroupOps.TwoCoboundaries( C, F, M );
    
    IsMat(b);

    E.2cobounds := rec();
    E.2cobounds.rowspace := RowSpace(GF(p), b, zero);

end;



InfoLift := Ignore;

AutGroupOps.LiftAut := function (E, pair)
    #
    # This nonsplit lifting code computes a system of linear equations
    # (inhomogeneous) that determine the values of a map psi:G->M which
    # allows a compatible pair to be lifted.
    #
    # We assume that a module element m_i is associated with each generator
    # g_i of G, and we compute an inhomogeneous system of linear equations
    # that the m_i satisfy if they define the correct coboundary to lift
    # the automorphism with.
    #
    # Generally not all relations are required to complete the system of
    # equations, so we will keep an eye out for when the rank of the system
    # has remained constant for a while, and try to solve the system as
    # early as possible. If we succeed we return the lifting early,
    # otherwise we continue processing more relations.
    #
    # It is also true that there is a small set of relations which produce
    # equations for each pair that needs to be lifted. We keep track of those
    # relations which produce equations in <E.liftAutInfo.goodrelations>. We
    # get an enormous saving this way.
    
    local makeaut, n, d, info, count, eqns, processrel, x, i, j, s, aut;

    makeaut := function (s)
        #
        # take a solution to the system and construct the images of the
        # lifted automorphism. Returns an aut record.

        local images, kappa, nu, i, aut;

        s := List([1..n], i -> s{(i-1)*d + [1..d]});
        images := [];
        kappa := pair.kappa; 
        nu := pair.nu;
        for i in [1..n] do
            Add(images, AutGroupOps.Map(E.extension, kappa.images[i])
                * AutGroupOps.vec2ext(E, s[i]));
        od;
        Append(images, List(nu, row -> AutGroupOps.vec2ext(E,row)));

        aut := AutOps.MakeAut(E.extension, images, pair.weight);
        return aut;
    end;


    n := E.n;  d := E.d;

    if not IsBound(E.liftAutInfo) then
        #
        # set up info that speeds up subsequent calls to this routine
        E.liftAutInfo := rec();
        E.liftAutInfo.tooshort := 0;
        E.liftAutInfo.goodrelations := [];
        E.liftAutInfo.enough := false;
    fi;
    info := E.liftAutInfo;
    count := 1;

    # set up a new system of equations
    #
    eqns := EqnOps.NewEqns(n*d, E.module.field);

    processrel := function (i,j,addgood)
        # this function does all the work

        # External: count, aut
        local ret, rows, vec, tmp, s;

        ret := AutGroupOps.LiftingEquations(E, pair, j, i);
        rows := Cat(ret.psi);
        vec := ret.vec;

        tmp := Length(eqns.mat);

        EqnOps.AddEqns(eqns, TransposedMat(rows), -vec);

        if Length(eqns.mat) > tmp then
            # reset counter for this new rank
            if addgood then
                Add(info.goodrelations, [i,j]);
            fi;
            count := 1;
        else
            count := count + 1;
        fi;

        if ( (info.enough <> false and Length(eqns.mat) >= info.enough) 
             or count > eqns.dim/10)
           and Length(eqns.mat) > info.tooshort then
            #
            # We possible have enough equations so far, or the equations
            # have remained at current rank for a reasonably long time (10%
            # of max possible rank) --- try to solve system now and check
            # whether we get an automorphism

            if (info.enough <> false and Length(eqns.mat) >= info.enough) then
                InfoLift("#L enough = ", info.enough, " try solving...\n");
            else
                InfoLift("#L had ", Length(eqns.mat)," for a while, ",
                        "try solving...\n");
            fi;

            s := EqnOps.SolveEqns(eqns);
            aut := makeaut(s);
            if IsAutomorphism(AutOps.hom(aut,false)) then
                info.enough := Length(eqns.mat);
                return true;
            fi;
            # otherwise, keep going
            info.tooshort := Length(eqns.mat);
            info.enough := false;
        fi;


        if eqns.fail = true then
            return true;
        fi;
        return false;
    end;


    InfoLift("#L LiftAut: Processing ", Length(info.goodrelations),
            " good relations first...\n");
    for x in info.goodrelations do
        if processrel(x[1], x[2], false) then
            # processrel returns true when it has the answer
            if eqns.fail = true then
                return false;
            else
                return aut;
            fi;
        fi;
    od;

    InfoLift("#L LiftAut: Processing other relations...\n");
    for i in [1..n] do
        for j in [i..n] do
            if not [i,j] in info.goodrelations then
                if processrel(i, j, true) then
                    # processrel returns true when it has the answer
                    if eqns.fail = true then
                        return false;
                    else
                        return aut;
                    fi;
                fi;
            fi;
        od;
    od;
    
    s := EqnOps.SolveEqns(eqns);
    if s = false then return false; fi;

    aut := makeaut(s);
    info.enough := Length(eqns.mat);

    if AutGroupOps.Debugging then
        if not IsAutomorphism(AutOps.hom(aut,false)) then
            Error("\n\nAutGroupOps.LiftAut: constructed aut is bad\n\n");
        fi;
    fi;
        
    return aut;

end;




AutGroupOps.LiftingEquations := function (E, pair, j, i)
    #
    # This function computes part of an inhomogeneous system of linear
    # equations to determine a 2-coboundary which must be used to lift a
    # compatible pair. It is called by AutGroupOps.LiftAut.
    #
    # We need to compute elements z_i of the module so that the pair
    #
    #  ( g_1 -> w_1,  g_2 -> w_2, ...  g_n -> w_n ;  nu )
    #
    # lifts to an automorphism of the extension:
    #
    #  ( g_1 -> w_1 z_1,  g_2 -> w_2 z_2 ... g_n -> w_n z_n,
    #    g_n+1 -> nu_1 ... )

    local n, d, F, null, id, gens, ggens, q, nu, kappa, psi, w, tail, 
          right, vec, k, e;

    n := E.n;  d := E.d;
    F := E.module.field;
    null := NullMat(d, d, F);
    id := IdentityMat(d, F);
    gens := E.extension.generators;
    ggens := E.group.generators;
    q := List(E.parent.weights{[1..n]}, x -> x[3]);
    
    nu := pair.nu;
    kappa := pair.kappa.images;

    psi := List([1..n], x -> null);

    w := E.rhs[j][i];

    if i = j then
        tail := AutGroupOps.ext2vec(E, gens[i]^q[i]);
    else
        tail := AutGroupOps.ext2vec(E, gens[j]*gens[i]);
    fi;

    right := ggens[1]^0;
    vec := tail * nu;

    for k in [n, n-1 .. 1] do
        for e in [1..w[k]] do
            psi[k] := psi[k] + AutGroupOps.Map(E.module.genimages, right);
            vec := vec + AutGroupOps.FactorSet(E, kappa[k], right);
            right := kappa[k]*right;
        od;
    od;

    if i = j then
        right := ggens[1]^0;
        for e in [1..q[i]] do
            psi[i] := psi[i] - AutGroupOps.Map(E.module.genimages, right);
            vec := vec - AutGroupOps.FactorSet(E, kappa[i], right);
            right := kappa[i] * right;
        od;
    else
        psi[i] := psi[i] - id;
        psi[j] := psi[j] - AutGroupOps.Map(E.module.genimages, kappa[i]);
        vec := vec - AutGroupOps.FactorSet(E, kappa[j], kappa[i]);
    fi;

    return rec(psi := psi, vec := vec);
end;




# the following function has been taken from the soluble quotient section
# of the standard GAP library and modified slightly to reorder the
# relations.  Note in particular that the relations are taken to be
#
#                           $g_j g_i = g_i ...$ 
#
# in order to eliminate the inversion of the matrix generators.

#############################################################################
##
#F  TwoCoboundariesSQ( <C>, <G>, <M> )  . . . .  compute the two coboundaries
##
AutGroupOps.TwoCoboundaries := function( C, G, M )
    local n, R, r, i, j, d, x, m, e, k;

    # InfoSQ2( "#I  computing two coboundaries\n" );

    # start with zero matrix
    n := Length(G.generators);
    R := [];
    r := n*(n+1)/2;
    for i  in [ 1 .. n ]  do
        R[i] := [];
        for j in [ 1 .. r ] do
            R[i][j] := C.mzero;
        od;
    od;

    # compute inverse generators
    M  := M.generators;
    d  := Length(M[1]);

    x := 0;
    # loop over all relators
    for i  in [ 1 .. n ]  do
        for j in  [ i .. n ]  do

            x := x + 1;
            #x := (j^2-j)/2 + i;

            # power relator
            if i = j  then
                m := C.mone;
                for e  in [ 1 .. C.orders[j] ]  do
                    R[j][x] := R[j][x] - m;  m := M[j] * m;
                od;

            # conjugate
            else
                R[j][x] := R[j][x] - M[i];
                R[i][x] := R[i][x] - C.mone;
            fi;

            # compute fox derivatives
            m := C.mone;
            r := C.relators[j][i];
            if (j <> i) then
                r := Cat([i,1],r);   # <-- don't need inverses anymore
            fi;
            if r <> 0  then
                for k  in [ Length(r)-1, Length(r)-3 .. 1 ]  do
                    for e  in [ 1 .. r[k+1] ]  do
                        R[r[k]][x] := R[r[k]][x] + m;
                        m := M[r[k]] * m;
                    od;
                od;
            fi;
        od;
    od;

    # make one list
    m := [];
    r := n*(n+1)/2;
    for i  in [ 1 .. n ]  do
        for k  in [ 1 .. d ]  do
            e := [];
            for j  in [ 1 .. r ]  do
                Append( e, R[i][j][k] );
            od;
            Add( m, e );
        od;
    od;

    # compute a base for <m>
    return BaseMat(m);

end;



AutGroupOps.ComputeDefinitions := function (E)
    #
    # Compute definitions for generators in the module. We take the p-group
    # in the current nilpotent layer and compute definitions that only
    # involve generators that centralise the module --- these definitions
    # are used later to compute actions on the module that are implied by
    # actions on the p-subgroup of the factor group.
    #
    # Let E be the extension with l-th layer of the LG-series as the module.
    #
    # This function find definitions for the generators in the l-th layer
    # of the LG-series.
    #
    # Returns a list. The elements of the list give definitions of each
    # generator in the l-th layer of the LG-series, as products of p-th
    # powers and commutators of generators in the public Sylow p-group of
    # the nilpotent layer containing the l-th layer.
    #
    # Definitions hold modulo the (l+1)-st layer of the LG-series.

    local G, l, first, next, d, weights, np, p, c, F, gens, top, bot, mat, 
          def, j, i, I, modulegendefs, n, defns, t, comm;

    G := E.parent;
    l := E.layer;

    first := G.first[l];
    next := G.first[l+1];
    d := next - first;

    weights := G.weights;

    np := weights[first][1];    # the nilpotent layer containing the l-th layer
    p := weights[first][3];
    c := weights[first][2];
    F := GF(p);

    gens := G.generators;

    # these are the generators at the top of the p-group (ie class 1 gens)
    top := PositionsProperty(weights, x -> x[1] = np and x[3] = p and x[2] = 1);

    # these are the generators from the class c-1 layer of the p-group
    bot := PositionsProperty(weights, x -> x[1] = np and x[3] = p and x[2] = c-1);

    mat := [];
    def := [];
    for j in bot do
        Add(mat, ExponentsAgWord(gens[j]^p,first,next-1,F.root));
        Add(def, [j]);
        for i in top do
            if i < j then
                Add(mat, ExponentsAgWord(
                        Comm(gens[j],gens[i]),first,next-1,F.root));
                Add(def, [j,i]);
            fi;
        od;
    od;

    I := IdentityMat(Length(mat),F);
    for i in [1..Length(mat)] do
        Append(mat[i], I[i]);
    od;

    TriangulizeMat(mat);

    modulegendefs := [];
    for i in [1..d] do
        if mat[i][i] <> F.one then
            Error("problem with finding definitions of module gens\n");
        fi;
        modulegendefs[i] := [];
        for j in [1 .. Length(def)] do
            if mat[i][d+j] <> F.zero then
                Add(modulegendefs[i], [ def[j], Int(mat[i][d+j]) ]);
            fi;
        od;
    od;

    # now change these to a different form

    n := E.n;  d := E.d;
    p := E.module.field.size;

    defns := [];
    for i in [1..d] do
        def := [];
        for t in modulegendefs[i] do
            if Length(t[1]) = 1 then
                # power of gen
                Append(def, List([1..t[2]*p], i -> t[1][1]));
            else
                # commutator of gens
                comm := [-t[1][1], -t[1][2], t[1][1], t[1][2]];
                Append(def, Cat(List([1..t[2]], i -> comm)));
            fi;
        od;
        defns[n+i] := def;
    od;
    E.defns := defns;
end;


AutGroupOps.ImpliedAction := function ( E, aut )
    #
    # Take a list of definitions for elements of a module layer and apply
    # the aut to them to find the action on the module implied by the
    # higher action.
    #
    # Note: assumes that the definition of each generator is in terms of
    # elements that centralise the module generators --- for example, only
    # involving generators from the p-group in the same nilpotent section
    # as the module.
    
    local n, d, p, defns, a, ans, i, m, x, e;
    
    n := E.n;  d := E.d;

    if not IsBound(E.defns) then
        AutGroupOps.ComputeDefinitions(E);
    fi;
    defns := E.defns;

    a := aut.images;

    ans := [];
    for i in [1..d] do

        m := E.extension.generators[1]^0;

        for x in defns[n+i] do
            e := SignInt(x);
            x := AbsInt(x);
            m := m * AutGroupOps.Transversal(E, a[x])^e;
        od;
        Add(ans, AutGroupOps.ext2vec(E, m));

    od;
    IsMat(ans);
    return ans;

end;
