###########################################################################
##
#A  extend.g                 autag package                 Michael J Smith
##
##  November 1996
##
##  This file forms part of a package for computing automorphism groups of
##  finite soluble groups which are given in terms of special soluble group
##  presentations.
##
###########################################################################


## Changes:

# 17:12 Fri 2 Aug 1996 - instead taking a shallow copy of the current
# extension record, define a new one and copy across the automorphism group
# order information.

# 14:04 Sun 28 Jul 1996 - improved debugging print information.


AutGroupOps.MakeExtension := function (Grp, l)
    #
    # <Grp> is a sag-group, <l> is the index of a layer of the LG-series of
    # <Grp>. Returns a record describing the extension of (l-1)-st quotient
    # of <Grp> by the l-th layer of the LG-series.
    
    local r, gens, first, next, n, d, S, H, A, G;
    
    r := rec();
    
    # Keep the parent group around
    r.parent := Grp;
    gens := Grp.generators;
    r.layer := l;

    # find the position of the l-th term of the LG-series
    first := Grp.first[l];
    next := Grp.first[l+1];
    n := first - 1;
    d := next - first;
    r.n := n;
    r.d := d;
    
    # make the extension
    if l < Maximum(Grp.layers) then
        # Kernel of map from whole group onto the extension
        S := AgSubgroup(Grp, 
                     gens{[n+d+1 .. Length(gens)]}, true);
        H := FactorGroup(Grp,S);
        H.name := Cat("G/G", String(l));
    else
        H := Grp;
    fi;
    r.extension := H;
    
    # make the group we are extending by
    A := AgSubgroup(H, H.generators{[n+1..n+d]},true);
    G := FactorGroup(H, A);
    G.name := Cat("G/G", String(l-1));
    r.subgroup := A;
    r.group := G;
    
    # define the module
    r.module := AutGroupOps.MakeModule(Grp, l);
    r.module.group := G;

    return r;
end;




AutGroupOps.NextExtension := function (E)
    # 
    # Take extension record <E> and return the next extension.

    local r, Grp, gens, l, first, next, n, d, S, H, A, G;

    r := rec();

    # copy over the parent group (ie our target)
    r.parent := E.parent;

    # copy over the order of the aut group of the current factor
    r.autOrder := E.autOrder;

    # copy over the list of orders of the subgroups in the 
    # aut group normal series
    r.autChain := E.autChain;
    
    Grp := r.parent;
    gens := Grp.generators;
    l := E.layer + 1;

    # next layer
    r.layer := l;

    # find the position of the l-th term of the LG-series
    first := Grp.first[l];
    next := Grp.first[l+1];
    n := first - 1;
    d := next - first;
    r.n := n;
    r.d := d;

    # Kernel of map from whole group onto the extension
    S := AgSubgroup(Grp, 
                 gens{[n+d+1 .. Length(gens)]}, true);
    
    # make the extension
    if l < Maximum(Grp.layers) then
        # Kernel of map from whole group onto the extension
        S := AgSubgroup(Grp, 
                     gens{[n+d+1 .. Length(gens)]}, true);
        H := FactorGroup(Grp,S);
        H.name := Cat("G/G", String(l));
    else
        H := Grp;
    fi;
    r.extension := H;
    
    # make the group we are extending by, and keep its relations
    A := AgSubgroup(H, H.generators{[n+1..n+d]},true);
    r.subgroup := A;

    G := E.extension;
    r.group := G;
    
    # define the module
    r.module := AutGroupOps.MakeModule(Grp, l);
    r.module.group := G;

    return r;
    
end;


AutGroupOps.ext2vec := function (E, g)
    # 
    # Take a element <g> of the extension group in <E> and return
    # a vector lying in the module.

    return ExponentsAgWord(g, E.n+1, E.n+E.d, E.module.field.root);
end;


AutGroupOps.vec2ext := function (E, m)
    #
    # Take a vector in the module of extension <E> and return the
    # corresponding element of the extension group.

    return Product([1..E.d], i -> E.extension.generators[E.n+i]^Int(m[i]));
end;


AutGroupOps.Transversal := function (E, g)
    return AutGroupOps.Map(E.extension, g);
end;


AutGroupOps.FactorSet := function (E, x, y)
    #
    # Return the evaluation of the factor set associated with extension <E>
    # on the pair of group elements <x>, <y> in E.group.

    return ExponentsAgWord(
                   AutGroupOps.Map(E.extension, x * y)^-1
                   * AutGroupOps.Map(E.extension, x)
                   * AutGroupOps.Map(E.extension, y),
                   E.n + 1, E.n + E.d, E.module.field.root);
end;



AutGroupOps.MakeModule := function( G, l )
    # 
    # This is a modified version of a routine from the sag-group library
    # files. It returns the action of the sag-group <G> on the <l>-th layer
    # of the LG-series of <G>. This version stores the action of all
    # generators of <G>, not just the ones in nilpotent-heads.

    local first, next, gens, N, M, NM, field, matgens, allmatgens, 
          groupgens, NMgens, i, mat, I;

    first := G.first[l];
    next := G.first[l+1];
    gens := G.generators;
    N := G.operations.AgSubgroup( G, gens{[first .. Length(gens)]}, true );
    M := G.operations.AgSubgroup( G, gens{[next  .. Length(gens)]}, true );
    NM := N mod M;
    field := GF ( G.weights[ first ][ 3 ] );
    matgens := [ ];
    allmatgens := [ ];
    groupgens := [ ];

    # take all necessary generators of G and calculate matrix representation
    NMgens := gens{[ first .. next-1 ]};
    for i  in [ 1 .. first-1 ]  do
        mat := List ( NMgens , x -> Exponents ( NM, x^gens[i], field ) );
        if     G.weights[i][2] = 1
               and G.weights[i][1] <> G.weights[first][1]
               then	
            Add ( matgens, mat );
            Add ( groupgens, i);
        fi;
        Add ( allmatgens, mat );
    od;

    if 0 = Length(matgens)  then
        I := IdentityMat( next-first, field );
        matgens := [ I ];
        groupgens := [ 1 ];
    fi;

    M := Gmodule(matgens);
    M.genimages := allmatgens;
    M.groupgens := groupgens;

    return M;
end;

