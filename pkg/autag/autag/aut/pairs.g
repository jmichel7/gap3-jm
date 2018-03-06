###########################################################################
##
#A  pairs.g                  autag package                 Michael J Smith
##
##  November 1996
##
##  This file forms part of a package for computing automorphism groups of
##  finite soluble groups which are given in terms of special soluble group
##  presentations.
##
###########################################################################


# First, functions for dealing with compatible pairs, ie elements
# of the direct product of Aut(G)xAut(M).


PairOps.IsCompatible := function (arg)
    local E, a, nu;
    E := arg[1];
    if Length(arg) > 2 then
        a := arg[2]; nu := arg[3];
    else
        a := arg[2].kappa; nu := arg[2].nu;
    fi;
    return not (false in List(E.group.generators{E.module.groupgens},
                   g -> AutGroupOps.Map(E.module.genimages,g)^nu
                   = AutGroupOps.Map(E.module.genimages, (g^a))));
end;


PairOps.MakePair := function (E, kappa, nu)
    # pair.kappa is aut record
    local r;
    r := rec();
    r.kappa := kappa;
    if AutOps.IsKnownInner(r.kappa) and AutOps.IsNontrivial(r.kappa) then
        r.inner := r.kappa.inner;
    fi;
    r.nu := nu;
    r.operations := PairOps;
    r.extRec := E;
    return r;
end;


PairOps.TrivialPair := function (E)
    local kappa, nu;
    kappa := AutOps.MakeAut(E.group, E.group.generators);
    nu := IdentityMat(E.module.dim, E.module.field);
    return PairOps.MakePair(E, kappa, nu);
end;


PairOps.Print := function (pair)
    Print (" ( ", pair.kappa," , ", pair.nu, ") ");
end;

PairOps.\* := function (pair1, pair2)
    return PairOps.MakePair(pair1.extRec, 
                   pair1.kappa * pair2.kappa,
                   pair1.nu * pair2.nu);
end;

PairOps.\^ := function (pair, n)
    local r;
    if n = 0 then
        r := PairOps.TrivialPair(pair.extRec);
    elif n = -1 and IsBound(pair.inverse) then
        r := pair.inverse;
    else
        r := PairOps.MakePair(pair.extRec, pair.kappa^n, pair.nu^n);
        if n = -1 then
            pair.inverse := r;
        fi;
    fi;
    return r;
end;

PairOps.\= := function (pair1, pair2)
    return (pair1.extRec = pair2.extRec) and
           (pair1.kappa = pair2.kappa) and 
           (pair1.nu = pair2.nu);
end;


PairOps.inverse := function (pair)
    if not IsBound(pair.inverse) then
        pair.inverse := PairOps.MakePair(pair.extRec,
                                pair.kappa^-1, pair.nu^-1);
    fi;
    return pair.inverse;
end;
