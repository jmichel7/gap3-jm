#############################################################################
##
#W  random.g		   	Matrix Packages                  Frank Celler
#W                                                              Eamon O'Brien
##
#H  @(#)$Id: random.g,v 1.1 1997/03/10 13:52:07 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains functions for pseudo random elements.
##
RevisionMatrix.random_g :=
    "@(#)$Id: random.g,v 1.1 1997/03/10 13:52:07 gap Exp $";


#############################################################################
##

#F  RandomSeedFlag( <grp> )
##
RandomSeedFlag := function( grp )
    if IsBound(grp.randomSeed)  then
        return grp.randomSeed;
    else
        return "unknown";
    fi;
end;


#############################################################################
##
#F  SetRandomSeedFlag( <grp>, <seed> )
##
SetRandomSeedFlag := function( grp, seed )
    grp.randomSeed := seed;
end;


#############################################################################
##
#F  GeneratorsFlag( <grp> )
##
GeneratorsFlag := function( grp )
    return grp.generators;
end;


#############################################################################
##
#F  IdentityFlag( <grp> )
##
IdentityFlag := function( grp )
    return grp.identity;
end;


#############################################################################
##

#F  InitPseudoRandom( <grp>, <seed_len>, <initial_scramble>
##
InitPseudoRandom := function( grp, len, scramble )
    local   gens,  seed,  i;

    # we need at least as many seeds as generators
    gens := GeneratorsFlag(grp);
    if 0 = Length(gens)  then
        SetRandomSeedFlag( grp, [] );
        return;
    fi;
    len := Maximum( len, Length(gens), 2 );

    # add random generators
    seed := ShallowCopy(gens);
    for i  in [ Length(gens)+1 .. len ]  do
        seed[i] := Random(gens);
    od;
    SetRandomSeedFlag( grp, seed );

    # scramble seed
    for i  in [ 1 .. scramble ]  do
        PseudoRandom(grp);
    od;

end;


#############################################################################
##
#F  PseudoRandom( <grp> )
##
PseudoRandom := function( grp )
    local   seed,  i,  j;

    # set up the seed
    seed := RandomSeedFlag(grp);
    if seed = "unknown"  then
        i := Length(GeneratorsFlag(grp));
        InitPseudoRandom( grp, i+10, Maximum( i*10, 100 ) );
        seed := RandomSeedFlag(grp);
    fi;
    if 0 = Length(seed)  then
        return IdentityFlag(grp);
    fi;

    # construct the next element
    i := Random([ 1 .. Length(seed) ]);
    repeat
        j := Random([ 1 .. Length(seed) ]);
    until i <> j;
    if Random([true,false])  then
        seed[j] := seed[i] * seed[j];
    else
        seed[j] := seed[j] * seed[i];
    fi;
    SetRandomSeedFlag( grp, seed );
    return seed[j];
end;


#############################################################################
##

#E  random.g  . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
##
