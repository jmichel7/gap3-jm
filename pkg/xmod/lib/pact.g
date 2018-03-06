##  pact.g,  version 29/ 6/96   (should be part of fpsgpres.g)
##           edited out 4 lines  #WN....  as per VF email (28/6/96)

#############################################################################
##
#F  PresentationAugmentedCosetTable( <aug>  [,<print level>] ) . . . create a
#F                                                              Tietze record
##
##  'PresentationAugmentedCosetTable'  creates a presentation,  i.e. a Tietze
##  record, from the given augmented coset table.
##
PresentationAugmentedCosetTable := function ( arg )

    local aug, coFacTable, comps, convert, gens, i, invs, lengths, numgens,
          numrels, pointers, printlevel, rel, rels, T, tietze, total, tree,
          treelength, treeNums;

    # check the first argument to be an augmented coset table.
    aug := arg[1];
    if not ( IsRec( aug ) and IsBound( aug.isAugmentedCosetTable ) and
        aug.isAugmentedCosetTable ) then
        Error( "first argument must be an augmented coset table" );
    fi;

    # check the second argument to be an integer.
    printlevel := 1;
    if Length( arg ) = 2 then  printlevel := arg[2];  fi;
    if not IsInt( printlevel ) then
        Error (" second argument must be an integer" );
    fi;

    # initialize some local variables.
    rels := Copy( aug.subgroupRelators );
    gens := Copy( aug.subgroupGenerators );
    coFacTable := aug.cosetFactorTable;
    tree := ShallowCopy( aug.tree );
    treeNums := Copy( aug.treeNumbers );
    treelength := Length( tree[1] );

    # create the Tietze record.
    T := rec( );
    T.isTietze := true;
    T.operations := PresentationOps;

    # construct the relator lengths list.
    numrels := Length( rels );
    lengths := 0 * [ 1 .. numrels ];
    total := 0;
    for i in [ 1 .. numrels ] do
        lengths[i] := Length( rels[i] );
        total := total + lengths[i];
    od;

    # initialize the Tietze stack.
    tietze := 0 * [ 1 .. TZ_LENGTHTIETZE ];
    tietze[TZ_NUMRELS] := numrels;
    tietze[TZ_RELATORS] := rels;
    tietze[TZ_LENGTHS] := lengths;
    tietze[TZ_FLAGS] := 1 + 0 * [ 1 .. numrels ];
    tietze[TZ_TOTAL] := total;

    # renumber the generators in the relators, if necessary.
    numgens := Length( gens );
    if numgens < treelength then
        convert := 0 * [ 1 .. treelength ];
        for i in [ 1 .. numgens ] do
            convert[treeNums[i]] := i;
        od;
        for rel in rels do
            for i in [ 1 .. Length( rel ) ] do
                if rel[i] > 0 then
                    rel[i] := convert[rel[i]];
                else
                    rel[i] := - convert[-rel[i]];
                fi;
            od;
        od;
    fi;

    # construct the generators and the inverses list, and save the generators
    # as components of the Tietze record.
    invs := 0 * [ 1 .. 2 * numgens + 1 ];
    comps := 0 * [ 1 .. numgens ];
    pointers := [ 1 .. treelength ];
    for i in [ 1 .. numgens ] do
        invs[numgens+1-i] := i;
        invs[numgens+1+i] := - i;
        T.(String( i )) := gens[i];
        comps[i] := i;
        pointers[treeNums[i]] := treelength + i;
    od;

    # define the remaining Tietze stack entries.
    tietze[TZ_NUMGENS] := numgens;
    tietze[TZ_GENERATORS] := gens;
    tietze[TZ_INVERSES] := invs;
    tietze[TZ_NUMREDUNDS] := 0;
    tietze[TZ_STATUS] := [ 0, 0, -1 ];
    tietze[TZ_MODIFIED] := false;

    # define some Tietze record components.
    T.generators := tietze[TZ_GENERATORS];
    T.tietze := tietze;
    T.components := comps;
    T.nextFree := numgens + 1;
    T.identity := IdWord;

    # initialize the Tietze options by their default values.
    T.eliminationsLimit := 100;
    T.expandLimit := 150;
    T.generatorsLimit := 0;
    T.lengthLimit := "infinity";
    T.loopLimit := "infinity";
    T.printLevel := 0;
    T.saveLimit := 10;
    T.searchSimultaneous := 20;

    # save the tree as component of the Tietze record.
    tree[TR_TREENUMS] := treeNums;
    tree[TR_TREEPOINTERS] := pointers;
    tree[TR_TREELAST] := treelength;
    T.tree := tree;

    # save the definitions of the primary generators as words in the original
    # group generators.
    T.primaryGeneratorWords := aug.primaryGeneratorWords;

    # handle relators of length 1 or 2, but do not eliminate any primary
    # generators.
    T.protected := tree[TR_PRIMARY];
    TzHandleLength1Or2Relators( T );
    T.protected := 0;

    # sort the relators.
    TzSort( T );
    T.printLevel := printlevel;

    # return the Tietze record.
    return T;
end;

