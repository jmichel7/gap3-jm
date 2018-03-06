#############################################################################
##
#A  cryst2.g                CrystGap library                     Bettina Eick
#A                                                              Franz G"ahler
#A                                                              Werner Nickel
##
#Y  Copyright 1990-1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  CrystGap - the crystallographic groups package for GAP (part 2)
##  

if not IsBound( InfoOperation1 ) then InfoOperation1 := Ignore; fi;
if not IsBound( InfoOperation2 ) then InfoOperation2 := Ignore; fi;

#############################################################################
##
#F  StabilizerInfiniteGroup( G, d, opr ) . . . . stabilizer of infinite group
#F                              without membership test, acting on finite set
##
StabilizerInfiniteGroup := function ( G, d, opr )
    local   stb,        # stabilizer, result
            orb,        # orbit
            rep,        # representatives for the points in the orbit <orb>
            set,        # orbit <orb> as set for faster membership test
            gen,        # generator of the group <G>
            pnt,        # point in the orbit <orb>
            img,        # image of the point <pnt> under the generator <gen>
            sch,        # schreier generator of the stabilizer
            lst;        # list of schreier generators

    # standard operation
    if   opr = OnPoints  then
        InfoOperation1("#I  Stabilizer |<gens>|=0\c");
        orb := [ d ];
        set := [ d ];
        rep := [ G.identity ];
        lst := [];
        for pnt  in orb  do
            for gen  in G.generators  do
                img := pnt ^ gen;
                if not img in set  then
                    Add( orb, img );
                    AddSet( set, img );
                    Add( rep, rep[Position(orb,pnt)]*gen );
                else
                    sch := rep[Position(orb,pnt)]*gen
                           / rep[Position(orb,img)];
                    AddSet( lst, sch );
                    InfoOperation2("\r#I  Stabilizer |<gens>|=",
                                   Length(stb.generators), "\c" );
                fi;
            od;
        od;
        InfoOperation1("\r#I  Stabilizer |<gens>|=",
                       Length(stb.generators),"\n");

    # compute iterated stabilizers for the operation on pairs or on tuples
    elif opr = OnPairs  or opr = OnTuples  then
        InfoOperation1("#I  Stabilizer |<fixed points>|=0\n");
        stb := G;
        for pnt in d  do
            stb := StabilizerInfiniteGroup( stb, pnt, OnPoints );
            InfoOperation2("#I  Stabilizer |<fixed points>|=",
                           Position( d, pnt ), "\n" );
        od;
        InfoOperation1("#I  Stabilizer |<fixed points>|=",
                       Length(d),"\n");

    # other operation
    else
        InfoOperation1("#I  Stabilizer |<gens>|=0\c");
        orb := [ d ];
        set := [ d ];
        rep := [ G.identity ];
        lst := [];
        for pnt  in orb  do
            for gen  in G.generators  do
                img := opr( pnt, gen );
                if not img in set  then
                    Add( orb, img );
                    AddSet( set, img );
                    Add( rep, rep[Position(orb,pnt)]*gen );
                else
                    sch := rep[Position(orb,pnt)]*gen
                           / rep[Position(orb,img)];
                    AddSet( lst, sch );
                    InfoOperation2("\r#I  Stabilizer |<gens>|=",
                                   Length(stb.generators), "\c" );
                fi;
            od;
        od;
        InfoOperation1("\r#I  Stabilizer |<gens>|=",
                       Length(stb.generators),"\n");

    fi;

    # return the stabilizer <stb>
    stb := G.operations.Subgroup( Parent( G ), lst );
    return stb;
end;


#############################################################################
##
#F  CrystGroupOps.Orbit( G, d, opr ) . . . . . orbit of a point under a group
##
CrystGroupOps.Orbit := function ( G, d, opr )
    local   orb,        # orbit, result
            set,        # orbit <orb> as set for faster membership test
            gen,        # generator of the group <G>
            pnt,        # point in the orbit <orb>
            img,        # image of the point <pnt> under the generator <gen>
            HasComparison;

    # if d is an infinite domain we cannot compare orbit elements
    HasComparison := not IsDomain( d ) or IsFinite( d );

    orb := [ d ];
    if HasComparison then
        set := [ d ];
    else
        set := orb;
    fi;

    # standard operation
    if   opr = OnPoints  then
        InfoOperation1("#I  Orbit |<d>^<G>|=\c");
        for pnt  in orb  do
            for gen  in G.generators  do
                img := pnt ^ gen;
                if not img in set  then
                    Add( orb, img );
                    if HasComparison then 
                        AddSet( set, img ); 
                    fi;
                fi;
            od;
        od;
        InfoOperation1("\r#I  Orbit |<d>^<G>|=",Length(orb),"\n");

    # special case for operation on pairs
    elif opr = OnPairs  then
        InfoOperation1("#I  Orbit |<pair>^<G>|=\c");
        for pnt  in orb  do
            for gen  in G.generators  do
                img := [ pnt[1]^gen, pnt[2]^gen ];
                if not img in set  then
                    Add( orb, img );
                    if HasComparison then 
                        AddSet( set, img ); 
                    fi;
                fi;
            od;
        od;
        InfoOperation1("\r#I  Orbit |<pair>^<G>|=",Length(orb),"\n");

    # other operation
    else
        InfoOperation1("#I  Orbit |opr(<d>,<G>)|=\c");
        for pnt  in orb  do
            for gen  in G.generators  do
                img := opr( pnt, gen );
                if not img in set  then
                    Add( orb, img );
                    if HasComparison then 
                        AddSet( set, img ); 
                    fi;
                fi;
            od;
        od;
        InfoOperation1("\r#I  Orbit |opr(<d>,<G>)|=",Length(orb),"\n");

    fi;

    # return the orbit <orb>
    return orb;
end;


#############################################################################
##
#F  CrystGroupOps.Orbits( G, D, opr )  . . . orbits of a domain under a group
##
CrystGroupOps.Orbits := function ( G, D, opr )

    local orbs, orb, HasComparison;

    # if D contains infinite domains we cannot compare orbit elements
    HasComparison := not IsDomain( D[1] ) or IsFinite( D[1] );

    if HasComparison then
        D := Set( D );
    fi;
    orbs := [];

    while D <> []  do
        orb := G.operations.Orbit( G, D[1], opr );
        Add( orbs, orb );
        if HasComparison then
            SubtractSet( D, orb );
        else
            D := Filtered( D, x -> not x in orb );
        fi;
    od;

    return orbs;

end;


CrystGroupOps.OrbitsOld := function ( G, D, opr )
    local   orbs,       # orbits, result
            orb,        # orbit
            set,        # orbit <orb> as set for faster membership test
            gen,        # generator of the group <G>
            pnt,        # point in the orbit <orb>
            img,        # image of the point <pnt> under the generator <gen>
            HasComparison;

    # if D contains infinite domains we cannot compare orbit elements
    HasComparison := not IsDomain( D[1] ) or IsFinite( D[1] );

    # standard operation
    if   opr = OnPoints  then
        InfoOperation1("#I  Orbits |orbs|=\c");
        if HasComparison then
            D := Set( D );
        fi;
        orbs := [];
        while D <> []  do
            orb := [ D[1] ];
            if HasComparison then
                set := [ D[1] ];
            else
                set := orb;
            fi;
            for pnt  in orb  do
                for gen  in G.generators  do
                    img := pnt ^ gen;
                    if not img in set  then
                        Add( orb, img );
                        if HasComparison then
                            AddSet( set, img );
                        fi;
                    fi;
                od;
            od;
            Add( orbs, orb );
            if HasComparison then
                SubtractSet( D, orb );
            else
                D := Filtered( D, x -> not x in orb );
            fi;
            InfoOperation2("\r#I  Orbits |orbs|=",Length(orbs),"\c");
        od;
        InfoOperation1("\r#I  Orbits |orbs|=",Length(orbs),"\n");

    # special case for operation on pairs
    elif opr = OnPairs  then
        InfoOperation1("#I  Orbits |orbs|=\c");
        if HasComparison then
            D := Set( D );
        fi;
        orbs := [];
        while D <> []  do
            orb := [ D[1] ];
            if HasComparison then
                set := [ D[1] ];
            else
                set := orb;
            fi;
            for pnt  in orb  do
                for gen  in G.generators  do
                    img := [ pnt[1]^gen, pnt[2]^gen ];
                    if not img in set  then
                        Add( orb, img );
                        if HasComparison then
                            AddSet( set, img );
                        fi;
                    fi;
                od;
            od;
            Add( orbs, orb );
            if HasComparison then
                SubtractSet( D, orb );
            else
                D := Filtered( D, x -> not x in orb );
            fi;
            InfoOperation2("\r#I  Orbits |orbs|=",Length(orbs),"\c");
        od;
        InfoOperation1("\r#I  Orbits |orbs|=",Length(orbs),"\n");

    # other operation
    else
        InfoOperation1("#I  Orbits |orbs|=\c");
        if HasComparison then
            D := Set( D );
        fi;
        orbs := [];
        while D <> []  do
            orb := [ D[1] ];
            if HasComparison then
                set := [ D[1] ];
            else
                set := orb;
            fi;
            for pnt  in orb  do
                for gen  in G.generators  do
                    img := opr( pnt, gen );
                    if not img in set  then
                        Add( orb, img );
                        if HasComparison then
                            AddSet( set, img );
                        fi;
                    fi;
                od;
            od;
            Add( orbs, orb );
            if HasComparison then
                SubtractSet( D, orb );
            else
                D := Filtered( D, x -> not x in orb );
            fi;
            InfoOperation2("\r#I  Orbits |orbs|=",Length(orbs),"\c");
        od;
        InfoOperation1("\r#I  Orbits |orbs|=",Length(orbs),"\n");

    fi;

    # return the orbits <orbs>
    return orbs;
end;


#############################################################################
##
#F  CrystGroupOps.OrbitLengths( G, D, opr )  lengths of the orbits of a group
##
CrystGroupOps.OrbitLengths := function ( G, D, opr )

    local lens, orb, HasComparison;

    # if D contains infinite domains we cannot compare orbit elements
    HasComparison := not IsDomain( D[1] ) or IsFinite( D[1] );

    if HasComparison then
        D := Set( D );
    fi;
    lens := [];

    while D <> []  do
        orb := G.operations.Orbit( G, D[1], opr );
        Add( lens, Length( orb ) );
        if HasComparison then
            SubtractSet( D, orb );
        else
            D := Filtered( D, x -> not x in orb );
        fi;
    od;

    return lens;

end;


CrystGroupOps.OrbitLengthsOld := function ( G, D, opr )
    local   lens,       # orbit lengths, result
            orb,        # orbit
            set,        # orbit <orb> as set for faster membership test
            gen,        # generator of the group <G>
            pnt,        # point in the orbit <orb>
            img,        # image of the point <pnt> under the generator <gen>
            HasComparison;

    # if D contains infinite domains we cannot compare orbit elements
    HasComparison := not IsDomain( D[1] ) or IsFinite( D[1] );

    InfoOperation1("#I  OrbitLengths |orbs|=\c");
    if HasComparison then
        D := Set( D );
    fi;
    lens := [];

    # standard operation
    if   opr = OnPoints  then
        while D <> []  do
            orb := [ D[1] ];
            if HasComparison then
                set := [ D[1] ];
            else
                set := orb;
            fi;
            for pnt  in orb  do
                for gen  in G.generators  do
                    img := pnt ^ gen;
                    if not img in set  then
                        Add( orb, img );
                        if HasComparison then
                            AddSet( set, img );
                        fi;
                    fi;
                od;
            od;
            Add( lens, Length(orb) );
            if HasComparison then
                SubtractSet( D, orb );
            else
                D := Filtered( D, x -> not x in orb );
            fi;
            InfoOperation2("\r#I  OrbitLengths |orbs|=",
                           Length(lens),"\c");
        od;
        InfoOperation1("\r#I  OrbitLengths |orbs|=",
                       Length(lens),"\n");

    # special case for operation on pairs
    elif opr = OnPairs  then
        while D <> []  do
            orb := [ D[1] ];
            if HasComparison then
                set := [ D[1] ];
            else
                set := orb;
            fi;
            for pnt  in orb  do
                for gen  in G.generators  do
                    img := [ pnt[1]^gen, pnt[2]^gen ];
                    if not img in set  then
                        Add( orb, img );
                        if HasComparison then
                            AddSet( set, img );
                        fi;
                    fi;
                od;
            od;
            Add( lens, Length(orb) );
            if HasComparison then
                SubtractSet( D, orb );
            else
                D := Filtered( D, x -> not x in orb );
            fi;
            InfoOperation2("\r#I  OrbitLengths |orbs|=",
                           Length(lens),"\c");
        od;
        InfoOperation1("\r#I  OrbitLengths |orbs|=",
                       Length(lens),"\n");

    # other operation
    else
        while D <> []  do
            orb := [ D[1] ];
            if HasComparison then
                set := [ D[1] ];
            else
                set := orb;
            fi;
            for pnt  in orb  do
                for gen  in G.generators  do
                    img := opr( pnt, gen );
                    if not img in set  then
                        Add( orb, img );
                        if HasComparison then
                            AddSet( set, img );
                        fi;
                    fi;
                od;
            od;
            Add( lens, Length(orb) );
            if HasComparison then
                SubtractSet( D, orb );
            else
                D := Filtered( D, x -> not x in orb );
            fi;
            InfoOperation2("\r#I  OrbitLengths |orbs|=",
                           Length(lens),"\c");
        od;
        InfoOperation1("\r#I  OrbitLengths |orbs|=",
                       Length(lens),"\n");

    fi;

    # return the orbit lengths <lens>
    return lens;
end;


#############################################################################
##
#F  CrystGroupOps.Stabilizer( G, d, opr ) stabilizer of a point under a group
##
CrystGroupOps.Stabilizer := function ( G, d, opr )
    local   stb,        # stabilizer, result
            orb,        # orbit
            rep,        # representatives for the points in the orbit <orb>
            set,        # orbit <orb> as set for faster membership test
            gen,        # generator of the group <G>
            pnt,        # point in the orbit <orb>
            img,        # image of the point <pnt> under the generator <gen>
            sch,        # schreier generator of the stabilizer
            HasComparison;

    # if d is an infinite domain, we cannot compare orbit elements
    HasComparison := not IsDomain( d ) or IsFinite( d );

    # standard operation
    if   opr = OnPoints  then
        InfoOperation1("#I  Stabilizer |<gens>|=0\c");
        orb := [ d ];
        if HasComparison then
            set := [ d ];
        else
            set := orb;
        fi;
        rep := [ G.identity ];
        stb := Subgroup( Parent( G ), [] );
        for pnt  in orb  do
            for gen  in G.generators  do
                img := pnt ^ gen;
                if not img in set  then
                    Add( orb, img );
                    if HasComparison then
                        AddSet( set, img );
                    fi;
                    Add( rep, rep[Position(orb,pnt)]*gen );
                else
                    sch := rep[Position(orb,pnt)]*gen
                           / rep[Position(orb,img)];
                    if not sch in stb  then
                        stb := stb.operations.Closure( stb, sch );
                        InfoOperation2("\r#I  Stabilizer |<gens>|=",
                                       Length(stb.generators), "\c" );
                    fi;
                fi;
            od;
        od;
        InfoOperation1("\r#I  Stabilizer |<gens>|=",
                       Length(stb.generators),"\n");

    # compute iterated stabilizers for the operation on pairs or on tuples
    elif opr = OnPairs  or opr = OnTuples  then
        InfoOperation1("#I  Stabilizer |<fixed points>|=0\n");
        stb := G;
        for pnt in d  do
            stb := stb.operations.Stabilizer( stb, pnt, OnPoints );
            InfoOperation2("#I  Stabilizer |<fixed points>|=",
                           Position( d, pnt ), "\n" );
        od;
        InfoOperation1("#I  Stabilizer |<fixed points>|=",
                       Length(d),"\n");

    # other operation
    else
        InfoOperation1("#I  Stabilizer |<gens>|=0\c");
        orb := [ d ];
        if HasComparison then
            set := [ d ];
        else
            set := orb;
        fi;
        rep := [ G.identity ];
        stb := Subgroup( Parent(G), [] );
        for pnt  in orb  do
            for gen  in G.generators  do
                img := opr( pnt, gen );
                if not img in set  then
                    Add( orb, img );
                    if HasComparison then
                        AddSet( set, img );
                    fi;
                    Add( rep, rep[Position(orb,pnt)]*gen );
                else
                    sch := rep[Position(orb,pnt)]*gen
                           / rep[Position(orb,img)];
                    if not sch in stb  then
                        stb := stb.operations.Closure( stb, sch );
                        InfoOperation2("\r#I  Stabilizer |<gens>|=",
                                       Length(stb.generators), "\c" );
                    fi;
                fi;
            od;
        od;
        InfoOperation1("\r#I  Stabilizer |<gens>|=",
                       Length(stb.generators),"\n");

    fi;

    # return the stabilizer <stb>
    return stb;
end;


#############################################################################
##
#F  CrystGroupOps.RepresentativeOperation(G,d,e,opr) . . repr. op. of a point
##
CrystGroupOps.RepresentativeOperation := function ( G, d, e, opr )
    local   rep,        # representative, result
            orb,        # orbit
            set,        # orbit <orb> as set for faster membership test
            gen,        # generator of the group <G>
            pnt,        # point in the orbit <orb>
            img,        # image of the point <pnt> under the generator <gen>
            by,         # <by>[<pnt>] is a gen taking <frm>[<pnt>] to <pnt>
            frm,        # where <frm>[<pnt>] lies earlier in <orb> than <pnt>
            HasComparison;

    # if d is infinite domain, we cannot compare orbit elements
    HasComparison := not IsDomain( d ) or IsFinite( d );

    orb := [ d ];
    if HasComparison then
        set := [ d ];
    else
        set := orb;
    fi;
    by  := [ G.identity ];
    frm := [ 1 ];

    # standard operation
    if   opr = OnPoints  then
        InfoOperation1("#I  RepresentativeOperation |<d>^<G>|=\c");
        if d = e  then return G.identity;  fi;
        for pnt  in orb  do
            for gen  in G.generators  do
                img := pnt ^ gen;
                if img = e  then
                    rep := gen;
                    while pnt <> d  do
                        rep := by[ Position(orb,pnt) ] * rep;
                        pnt := frm[ Position(orb,pnt) ];
                    od;
                    InfoOperation1("\r#I  RepresentativeOperation returns ",
                                   rep, "\n" );
                    return rep;
                elif not img in set  then
                    Add( orb, img );
                    if HasComparison then
                        AddSet( set, img );
                    fi;
                    Add( frm, pnt );
                    Add( by,  gen );
                fi;
            od;
        od;
        InfoOperation1("\r#I  RepresentativeOperation returns 'false'\n");
        return false;

    # special case for operation on pairs
    elif opr = OnPairs  then
        InfoOperation1("#I  RepresentativeOperation |<pair>^<G>|=\c");
        if d = e  then return G.identity;  fi;
        for pnt  in orb  do
            for gen  in G.generators  do
                img := [ pnt[1]^gen, pnt[2]^gen ];
                if img = e  then
                    rep := gen;
                    while pnt <> d  do
                        rep := by[ Position(orb,pnt) ] * rep;
                        pnt := frm[ Position(orb,pnt) ];
                    od;
                    InfoOperation1("\r#I  RepresentativeOperation returns ",
                                   rep, "\n" );
                    return rep;
                elif not img in set  then
                    Add( orb, img );
                    if HasComparison then
                        AddSet( set, img );
                    fi;
                    Add( frm, pnt );
                    Add( by,  gen );
                fi;
            od;
        od;
        InfoOperation1("\r#I  RepresentativeOperation returns 'false'\n");
        return false;

    # other operation
    else
        InfoOperation1("#I  RepresentativeOperation |opr(<d>,<G>)|=\c");
        if d = e  then return G.identity;  fi;
        for pnt  in orb  do
            for gen  in G.generators  do
                img := opr( pnt, gen );
                if img = e  then
                    rep := gen;
                    while pnt <> d  do
                        rep := by[ Position(orb,pnt) ] * rep;
                        pnt := frm[ Position(orb,pnt) ];
                    od;
                    InfoOperation1("\r#I  RepresentativeOperation returns ",
                                   rep, "\n" );
                    return rep;
                elif not img in set  then
                    Add( orb, img );
                    if HasComparison then
                        AddSet( set, img );
                    fi;
                    Add( frm, pnt );
                    Add( by,  gen );
                fi;
            od;
        od;
        InfoOperation1("\r#I  RepresentativeOperation returns 'false'\n");
        return false;

    fi;

end;


#############################################################################
##
#F  CrystGroupOps.RepresentativesOperation( G, d, opr ) rep. ops. of an orbit
##
CrystGroupOps.RepresentativesOperation := function ( G, d, opr )
    local   reps,       # representatives for the points in the orbit, result
            orb,        # orbit
            set,        # orbit <orb> as set for faster membership test
            gen,        # generator of the group <G>
            pnt,        # point in the orbit <orb>
            img,        # image of the point <pnt> under the generator <gen>
            HasComparison;

    # if d is an infinite domain we cannot compare orbit elements
    HasComparison := not IsDomain( d ) or IsFinite( d );

    orb := [ d ];
    if HasComparison then
        set := [ d ];
    else
        set := orb;
    fi;
    reps := [ G.identity ];

    InfoOperation1("#I  RepresentativesOperation |<orb>|=\c");

    # standard operation
    if   opr = OnPoints  then
        for pnt  in orb  do
            for gen  in G.generators  do
                img := pnt ^ gen;
                if not img in set  then
                    Add( orb, img );
                    if HasComparison then
                        AddSet( set, img );
                    fi;
                    Add( reps, reps[Position(orb,pnt)]*gen );
                fi;
            od;
        od;

    # other operation
    else
        for pnt  in orb  do
            for gen  in G.generators  do
                img := opr( pnt, gen );
                if not img in set  then
                    Add( orb, img );
                    if HasComparison then
                        AddSet( set, img );
                    fi;
                    Add( reps, reps[Position(orb,pnt)]*gen );
                fi;
            od;
        od;

    fi;

    InfoOperation1("\r#I  RepresentativesOperation |<orb>|=",
                   Length(orb),"\n");

    # return the representatives <reps>
    return reps;
end;


#############################################################################
##
#F  CrystGroupOps.Operation( G, D, opr ) . . . . . . . . . . . . . .Operation
##
CrystGroupOps.Operation := function ( G, D, opr )
    local   grp,        # operation group, result
            prms,       # generators of the operation group <grp>
            prm,        # generator of the operation group <grp>
            gen,        # generator of the group <G>
            i;          # loop variable

    # standard operation
    if   opr = OnPoints  then
        InfoOperation1("#I  Operation(<G>,<D>)              \c");

        # make a sorted copy <set> of the domain <D> and remember <pos>
        InfoOperation2("\r#I  Operation(<G>,<D>) sorting <D>\c");

        # make the permutations <prms> and the operation group <grp>
        prms := [];
        for gen  in G.generators  do
            InfoOperation2("\r#I  Operation(<G>,<D>) make perm ",
                           Position( G.generators, gen ), "\c");
            prm := [];
            for i  in [1..Length(D)]  do
                prm[i] := Position( D, D[i]^gen );
            od;
            Add( prms, PermList( prm ) );
        od;
        grp := Group( prms, () );
        grp.operation := rec( genimages := prms );
        InfoOperation1("\r#I  Operation(<G>,<D>) returns       \n");

    # special case for operation on pairs
    elif opr = OnPairs  then
        InfoOperation1("#I  Operation(<G>,<D>)              \c");

        # make a sorted copy <set> of the domain <D> and remember <pos>
        InfoOperation2("\r#I  Operation(<G>,<D>) sorting <D>\c");

        # make the permutations <prms> and the operation group <grp>
        prms := [];
        for gen  in G.generators  do
            InfoOperation2("\r#I  Operation(<G>,<D>) make perm ",
                           Position( G.generators, gen ), "\c");
            prm := [];
            for i  in [1..Length(D)]  do
                prm[i] := Position( D,[ D[i][1]^gen, D[i][2]^gen ] );
            od;
            Add( prms, PermList( prm ) );
        od;
        grp := Group( prms, () );
        grp.operation := rec( genimages := prms );
        InfoOperation1("\r#I  Operation(<G>,<D>) returns       \n");

    # other operation
    else
        InfoOperation1("#I  Operation(<G>,<D>)              \c");

        # make a sorted copy <set> of the domain <D> and remember <pos>
        InfoOperation2("\r#I  Operation(<G>,<D>) sorting <D>\c");

        # make the permutations <prms> and the operation group <grp>
        prms := [];
        for gen  in G.generators  do
            InfoOperation2("\r#I  Operation(<G>,<D>) make perm ",
                           Position( G.generators, gen ), "\c");
            prm := [];
            for i  in [1..Length(D)]  do
                prm[i] := Position( D, opr( D[i], gen ) );
            od;
            Add( prms, PermList( prm ) );
        od;
        grp := Group( prms, () );
        grp.operation := rec( genimages := prms );
        InfoOperation1("\r#I  Operation(<G>,<D>) returns       \n");

    fi;

    # return the permutation group <grp>
    return grp;
end;


#############################################################################
##
#F  CrystGroupOpHomOps . . . . . . . . . . . . . . . . . . CrystGroupOpHomOps
##
CrystGroupOpHomOps := ShallowCopy( OperationHomomorphismOps );
CrystGroupOpHomOps.name := "CrystGroupOpHomOps";


#############################################################################
##
#F  CrystGroupOpHomOps.ImageElm( hom, elm ) . . . . . . . . . . . . .ImageElm
##
CrystGroupOpHomOps.ImageElm := function( hom, elm )
    local D, opr, prm, i;
    D   := hom.range.operation.domain;
    opr := hom.range.operation.operation;
    prm := [1..Length( D )];
    for i in [1..Length( D )] do
        prm[i] := Position( D, opr( D[i], elm ) );
    od;
    return PermList( prm );
end;


#############################################################################
##
#F  CrystGroupOpHomOps.ImagesElm( hom, elm ) . . . . . . . . . . . .ImagesElm
##
CrystGroupOpHomOps.ImagesElm := function( hom, elm )
    return [ CrystGroupOpHomOps.ImagesElm( hom, elm ) ];
end;


#############################################################################
##
#F  CrystGroupOps.IsFinite( G ) . . . . . . . . . . . . . . . . . . .IsFinite
##
CrystGroupOps.IsFinite := function( G )
    return TranslationsCrystGroup( G ) = [];
end;


#############################################################################
##
#F  CrystGroupOps.IsSolvable( S ) . . . . . . . . . . . . . . . . .IsSolvable
##
CrystGroupOps.IsSolvable := function( S )
    return IsSolvable( PermGroup( PointGroup( S ) ) );
end;


#############################################################################
##
#F  CrystGroupOps.Index( G, H ) . . . . . . . . . . . . . . . . . . . . Index
##
CrystGroupOps.Index := function( G, H )
    if IsFinite( G ) then
        return Size( G ) / Size( H );
    elif Length( TranslationsCrystGroup( G ) ) <>
             Length( TranslationsCrystGroup( H ) ) then
        return "infinity";
    else
        return Length( RightCosets( G, H ) );
    fi;
end;


#############################################################################
##
#F  CrystGroupOps.Group( D ) . . . . . . . . .convert a subgroup into a group
##
CrystGroupOps.Group := function ( D )
    local G;
    G := CrystGroup( GroupOps.Group( D ) );
    if IsBound( D.translations ) then
        AddTranslationsCrystGroup( G, D.translations );
    fi;
    return G;
end;


#############################################################################
##
#F  CrystGroupOps.RightCoset( H ) . . . . . . . . . . . . . . . . .RightCoset
##
CrystGroupOps.RightCoset := GroupOps.RightCoset;


#############################################################################
##
#F  CrystGroupOps.RightCosets( G, H ) . . . . . . . . . . . . . . RightCosets
##
CrystGroupOps.RightCosets := function( G, H )

    local orb, pnt, img, gen, oprec;

    # first some simple checks
    if IsFinite( G ) then
        return GroupOps.RightCosets( G, H );
    elif Length( TranslationsCrystGroup( G ) ) <>
             Length( TranslationsCrystGroup( H ) ) then
        Error("sorry, there are infinitely many cosets");
    fi;

    # we cannot compare infinite cosets, so we use a simplified algorithm
    orb := [ RightCoset( H ) ];
    for pnt in orb do
        for gen in G.generators do
            img := OnRight( pnt, gen );
            if not img in orb  then
                Add( orb, img );
            fi;
        od;
    od;

    # make the orbit a sorted list
    oprec := ShallowCopy( orb[1].operations );
    oprec.\< := function( a, b )
                    return Position( orb, a ) < Position( orb, b );
                end;
    for pnt in orb do
        pnt.operations := oprec;
    od;
    return orb;

end;    


#############################################################################
##
#F  CrystGroupOps.Normalizer( G, H ) . . . . . . . . Normalizer in CrystGroup
##
CrystGroupOps.Normalizer := GroupOps.Normalizer;


#############################################################################
##
#F  ElementsConjugacyClassInfiniteSubgroups( C ) .Elements of conjugacy class 
#F  . . . . . . . . . . . . . . . . . . . . . . . . . . of infinite subgroups
##
ElementsConjugacyClassInfiniteSubgroups := function ( C )

    local e, elems, oprec;

    if not IsBound( C.normalizer )  then
        C.normalizer := Normalizer( C.group, C.representative );
    fi;

    # get the list of elements
    elems := List( RightTransversal( C.group, C.normalizer ),
                   t -> C.representative^t );

    # make it a sorted list
    oprec := ShallowCopy( elems[1].operations );
    oprec.\< := function( a, b )
                    return Position( elems, a ) < Position( elems, b );
                end;
    for e in elems do
        e.operations := oprec;
    od;
    return elems;
     
end;


#############################################################################
##
#F  CrystGroupOps.ConjugacyClassSubgroups( G, H ) . . ConjugacyClassSubgroups
##
CrystGroupOps.ConjugacyClassSubgroups := function( G, H )

    local C;

    C := GroupOps.ConjugacyClassSubgroups( G, H );
    if not IsFinite( H ) then
        C.operations := ShallowCopy( C.operations );
        C.operations.Elements := ElementsConjugacyClassInfiniteSubgroups;
    fi;
    return C;

end;


#############################################################################
##
#F  CrystGroupOps.Closure( G, obj ) .  closure of group with group or element
##
CrystGroupOps.Closure := function ( G, obj )
    local   C,          # closure '\< <G>, <obj> \>', result
            gen,        # generator of <G> or <C>
            reps,       # representatives of cosets of <G> in <C>
            rep;        # representative of coset of <G> in <C>

    # handle the closure of a group with a subgroup
    if IsGroup( obj )  then
        if IsParent( G )  then
            C := G;
        elif IsParent( obj )  then
            C := obj;
        else
            C := G;
            for gen in obj.generators  do
                C := G.operations.Closure( C, gen );
            od;
        fi;

    # closure of a group with a single element
    else

        # try to avoid adding an element to a group that already contains it
        if   IsParent( G )
          or obj in G.generators
          or obj^-1 in G.generators
          or obj in G
        then
            return G;
        fi;

        # make the closure group
        C := G.operations.Subgroup( Parent(G),
                                    Concatenation( G.generators, [obj] ) );

        # if <G> is nonabelian then so is <C>
        if IsBound( G.isAbelian ) and not G.isAbelian  then
            C.isAbelian := false;
        elif IsBound( G.isAbelian )  then
            C.isAbelian := ForAll( G.generators,
                                   gen -> Comm(gen,obj)=G.identity );
        fi;

        # if <G> is infinite then so is <C>
        if IsBound( G.isFinite ) and not G.isFinite  then
            C.isFinite := false;
            C.size     := "infinity";
        fi;

        # if the elements of <G> are known then extend this list
        if IsBound( G.elements ) and IsFinite( Parent( G ) ) then
            # if <G>^<obj> = <G> then <C> = <G> * <obj>
            if ForAll( G.generators, gen -> gen ^ obj in G.elements )  then
                InfoGroup2( "#I   new generator normalizes\n" );
                C.elements := ShallowCopy( G.elements );
                rep := obj;
                while not rep in G.elements  do
                    Append( C.elements, G.elements * rep );
                    rep := rep * obj;
                od;
                C.elements := Set( C.elements );
                C.isFinite := true;
                C.size     := Length( C.elements );

            # otherwise use a Dimino step
            else
                InfoGroup2( "#I   new generator normalizes not\n" );
                C.elements := ShallowCopy( G.elements );
                reps := [ G.identity ];
                InfoGroup2( "\r#I   |<cosets>| = ",Length(reps), "\c" );
                for rep  in reps  do
                    for gen  in C.generators  do
                        if not rep * gen in C.elements  then
                            Append( C.elements, G.elements * (rep * gen) );
                            Add( reps, rep * gen );
                            InfoGroup2( "\r#I   |<cosets>| = ",
                                        Length(reps),"\c");
                        fi;
                    od;
                od;
                InfoGroup2( "\n" );
                C.elements := Set( C.elements );
                C.isFinite := true;
                C.size := Length( C.elements );

            fi;
        fi;

    fi;

    # return the closure
    return C;

end;


#############################################################################
##
#F  CrystGroupOps.DerivedSubgroup . . . . . . . . . . . . . . DerivedSubgroup
##
CrystGroupOps.DerivedSubgroup := GroupOps.DerivedSubgroup;


#############################################################################
##
#F  CrystGroupOps.CommutatorSubgroup . . . . . . . . . . . CommutatorSubgroup 
##
CrystGroupOps.CommutatorSubgroup := GroupOps.CommutatorSubgroup;


#############################################################################
##
#F  CrystGroupOps.CompositionSeries . . . . . . . . . . . . CompositionSeries
##
CrystGroupOps.CompositionSeries := GroupOps.CompositionSeries;


#############################################################################
##
#F  CrystGroupOps.ConjugacyClasses( G ) . . . . . . . . . . .ConjugacyClasses
##
CrystGroupOps.ConjugacyClasses := function( G )
    if IsFinite( G ) and IsFinite( Parent( G ) ) then
        return MatGroupOps.ConjugacyClasses( G );
    elif IsFinite( G ) then
        return GroupOps.ConjugacyClasses( G );
    else
        Error("sorry, there are infinitely many conjugacy classes");
    fi;
end;


#############################################################################
##
#F  CrystGroupOps.PermGroup( G ) . . . . . . . . . . . . . . . . . .PermGroup
##
CrystGroupOps.PermGroup := function( G )
    if IsFinite( G ) then
        return MatGroupOps.PermGroup( G );
    else
        Error("sorry, this group is not finite");
    fi;
end;


#############################################################################
##
#F  CrystGroupOps.SylowSubgroup( G, p ) . . . . . . . . . . . . SylowSubgroup
##
CrystGroupOps.SylowSubroup := function( G, p )
    if IsFinite( G ) then
        return MatGroupOps.SylowSubgroup( G, p );
    else
        Error("sorry, this group is not finite");
    fi;
end;


#############################################################################
##
#F  Unbind functions which make no sense or don't work 
##
Unbind(CrystGroupOps.KroneckerProduct);
Unbind(CrystGroupOps.InvariantForm);
Unbind(CrystGroupOps.Random);


#############################################################################
##
#F  NormalizerGL( G ) . . . . .Normalizer of integral matrix group in GL(d,Z)
##
NormalizerGL := function( G )
    if not IsBound( G.normalizerGL ) then
        G.normalizerGL := G.operations.NormalizerGL( G );
    fi;
    return G.normalizerGL;
end;


#############################################################################
##
#F  MatGroupOps.NormalizerGL( G ) . . . . . . . . . . . Normalizer in GL(d,Z)
##
MatGroupOps.NormalizerGL := function( G )

    local S, t, M, N, B;

    # G is the point group of a space group from the library
    if IsBound( G.crystGroup ) then
        S := G.crystGroup;
        if IsBound( S.crTransposedSpaceGroupTyp ) then
            t := S.crTransposedSpaceGroupTyp;
            M := NormalizerZClass( t[1], t[2], t[3], t[4] );
            N := Group( List( M.generators, TransposedMat ), M.identity );
            N.size := M.size;
            return N;
        fi;
    else
        Error("Sorry, NormalizerGL currently is implemented only for\n",
              "point groups of space groups from the library");
    fi;

    # G knowns its Bravais group B, and the NormalizerGL of B
    if IsBound( G.bravaisGroup ) then
        B := G.bravaisGroup;
        if IsBound( B.normalizerGL ) then
            N := B.normalizerGL;
            if IsBound( N.isFinite ) and N.isFinite or
                     IsBound( N.size ) and N.size <> "infinity" then
                N := N.operations.Stabilizer( N, G, OnPoints );
            else
                N := StabilizerInfiniteGroup( N, G, OnPoints );
            fi;
        else
            Error("NormalizerGL of the Bravais group of G is unknown");
        fi;
    else
        Error("Bravais group of G and its NormalizerGL are unknown");
    fi;

    # return the result
    return N;

end;


#############################################################################
##
#F  PointGroupOps.NormalizerGL( G ) . . . . . . . . . . Normalizer in GL(d,Z)
##
PointGroupOps.NormalizerGL := MatGroupOps.NormalizerGL;


#############################################################################
##
#F  CentralizerGL( G ) . . . .Centralizer of integral matrix group in GL(d,Z)
##
CentralizerGL := function( G )
    if not IsBound( G.centralizerGL ) then
        G.centralizerGL := G.operations.CentralizerGL( G );
    fi;
    return G.centralizerGL;
end;


#############################################################################
##
#F  MatGroupOps.CentralizerGL( G ) . . . . . . . . . . Centralizer in GL(d,Z)
##
MatGroupOps.CentralizerGL := function( G )

    local N, S, t, B, C;

    # get the NormalizerGL of G, or that of its Bravais group
    if IsBound( G.normalizerGL ) then
        N := G.normalizerGL;

    # G is the point group of a space group from the library
    elif IsBound( G.crystGroup ) then
        S := G.crystGroup;
        if IsBound( S.crTransposedSpaceGroupTyp ) then
            t := S.crTransposedSpaceGroupTyp;
            N := NormalizerZClass( t[1], t[2], t[3], t[4] );
            N.generators := List( N.generators, x -> TransposedMat( x ) );
        fi;
    else
        Error("Sorry, CentralizerGL currently is implemented only for\n",
              "point groups of space groups from the library");
    fi;

#    # G knowns its Bravais group B, and the NormalizerGL of B
#    elif IsBound( G.bravaisGroup ) then
#        B := G.bravaisGroup;
#        if IsBound( B.normalizerGL ) then
#            N := B.normalizerGL;
#        else
#            Error("NormalizerGL of Bravais group of G is unknown");
#        fi;
#    else
#        Error("Bravais group of G is unknown");
#    fi;

    # now get the CentralizerGL
    if IsBound( N.isFinite ) and N.isFinite or
          IsBound( N.size ) and N.size <> "infinity" then
        C := N.operations.Stabilizer( N, G.generators, OnTuples );
    else
        C := StabilizerInfiniteGroup( N, G.generators, OnTuples );
    fi;
    return C;

end;


#############################################################################
##
#F  PointGroupOps.CentralizerGL( G ) . . . . . . . . . Centralizer in GL(d,Z)
##
PointGroupOps.CentralizerGL := MatGroupOps.CentralizerGL;


#############################################################################
##
#F  QuadFormEquations( m ) . . . .solutions of m * X * TransposedMat( m ) = X
##
##  Defining quations for the entries of all matrices
##  X with m^-1 * X - X * TransposedMat( m ) = 0
##
QuadFormEquations := function( m )

local d, t, e, n, r, i, j;

   d := Length( m );
   t := m^-1;
   e := IdentityMat( d );
   n := NullMat( d, d );
   r := Copy( n );
   for i in [1..d] do
       for j in [1..d] do
           if i = j then  
               r[i][j] := e*t[i][j] - m; 
           else 
               r[i][j] := e*t[i][j]; 
           fi; 
       od;
   od;
   return FlattenMatMat( r );

end;


#############################################################################
##
#F  QuadFormSpace( G ) . . . . . . space of quadratic forms invariant under G
##
QuadFormSpace := function( G )

   local d, M, S, i, j, k, Q, N;

   d := Parent( G ).dimension;

   # equations determining the quadratic forms
   M := Concatenation( List( G.generators, QuadFormEquations ) );

   # we want only symmetric solutions
   S := NullMat( d*(d-1)/2, d*d );
   k := 1;
   for i in [1..d-1] do
       for j in [i+1..d] do
           S[k][(i-1)*d+j] := 1;
           S[k][(j-1)*d+i] :=-1;
           k := k+1;
       od;
   od;
   M := Concatenation( M, S );

   Q := IdentityMat( Length(M[1]) );
    
   # first diagonalize M
   M := RowEchelonForm( M );
   while not IsDiagonalMat( M ) do
       M := TransposedMat( M );
       M := RowEchelonFormT( M, Q );
       if not IsDiagonalMat( M ) then
           M := TransposedMat( M );
           M := RowEchelonForm( M );
       fi;
   od;

   # and then determine its kernel (that of the original M)
   if Length( M ) < Length( Q ) then
       N := Q{[ Length( M )+1..Length( Q ) ]};
   else
       N := [];
   fi; 

   return List( N, x -> List( [1..d], i -> x{(i-1)*d+[1..d]} ) );

end;


#############################################################################
##
#F  PointGroupsBravaisClass( <B> [, <norm> ] ) . reps of subgroups in Bravais
#F                                                   class of Bravais group B
##
PointGroupsBravaisClass := function( arg )

    local B, norm, P, A, gen, h, r, d, g, N, opr, new, orb;

    B := arg[1];
    if Length(arg)>1 then 
       norm := Filtered( arg[2], x -> not x in B ); 
    else 
       norm := []; 
    fi;

    # get subgroup conjugacy classes
    P := PermGroup( B );
    if IsSolvable( P ) then
        A := AgGroup( P );
        gen := List( A.generators, 
                     x -> Image( P.bijection, Image( A.bijection, x ) ) );
        h := GroupHomomorphismByImages( B, A, gen, A.generators );
    else
        A := P;
        h := GroupHomomorphismByImages( B, A, A.generators,
                     List( A.generators, x -> Image( A.bijection, x ) ) );
    fi;
    r := List( ConjugacyClassesSubgroups( A ), Representative );
    r := List( r, x -> PreImage( h, x ) );
    r := List( r, x -> x.operations.Group( x ) );

    # filter those with the same space of quadratic forms
    d := Length( QuadFormSpace( B ) );
    r := Filtered( r, x -> Length( QuadFormSpace( x ) ) = d );

    # eliminate those equivalent under the normalizer
    if norm <> [] then
        N := Group( Concatenation( B.generators, norm ), B.identity );
        opr := function( G, g ) 
            return Group( List( G.generators, x -> x^g ), G.identity ); 
        end;
        new := [];
        while r <> [] do
            orb := Orbit( N, r[1], opr );
            Add( new, r[1] );
            r := Filtered( r, x -> not x in orb );
        od;
        r := new;
    fi;

    for g in r do
        g.bravaisGroup := B;
    od;

    return r;

end;


#############################################################################
##
#F  UnionModule( M1, M2 ) . . . . . . . . . . . . union of two free Z-modules
##
UnionModule := function( M1, M2 )
    return ReducedTranslationBasis( Concatenation( M1, M2 ) );
end;


#############################################################################
##
#F  IntersectionModule( M1, M2 ) . . . . . intersection of two free Z-modules
##
IntersectionModule := function( M1, M2 )

    local M, Q, T;

    if M1 = [] or M2 = [] then
        return [];
    fi;
    M := Concatenation( M1, M2 );
    M := M*Lcm( List( Flat( M ), Denominator ) );
    Q := IdentityMat( Length( M ) );
    M := RowEchelonFormT( M, Q );
    T := Q{[Length(M)+1..Length(Q)]}{[1..Length(M1)]}*M1;
    return ReducedTranslationBasis( T );

end;


#############################################################################
##
#F  VectorModL . . . . . . . . . . . . . . . . .vector modulo a free Z-module
##
VectorModL := function( v, L )

    local l, i, x, j;

    for l in L do
        i := PositionProperty( l, x -> x<>0 );
        x := v[i]/l[i];
        j := Int( x );
        if x < 0 and not IsInt( x ) then j := j-1; fi;
        v := v - j*l;
    od;
    return v;

end;


#############################################################################
##
#F  IntSolutionMat( M, b ) . . integer solution for inhom system of equations
##
IntSolutionMat := function( M, b )

    local b, Q, M, den, sol, i, x;

    if M = [] then
        return false;
    fi;

    b := ShallowCopy(b);
    Q := IdentityMat( Length(M) );

    den := Lcm( List( Flat( M ), x -> Denominator( x ) ) );
    if den <> 1 then
        M := den*M;
        b := den*b;
    fi;

    M := TransposedMat(M);
    M := RowEchelonFormVector( M,b );
    while not IsDiagonalMat(M) do
        M := TransposedMat(M);
        M := RowEchelonFormT(M,Q);
        if not IsDiagonalMat(M) then
            M := TransposedMat(M);
            M := RowEchelonFormVector(M,b);
        fi;
    od;

    # are there integer solutions?
    sol:=[];
    for i in [1..Length(M)] do
        x := b[i]/M[i][i];
        if IsInt( x ) then
            Add( sol, x );
        else 
            return false;
        fi;
    od;

    # are there solutions at all?
    for i in [Length(M)+1..Length(b)] do
        if b[i]<>0 then
            return false;
        fi;
    od;

    return sol*Q{[1..Length(sol)]};

end;
    

#############################################################################
##
#F  CrystGroupOps.Intersection( G, obj ) . . .intersection of two CrystGroups
##
CrystGroupOps.Intersection := function( G1, G2 )

    local G, d, P1, P2, P, T1, T2, T, L, gen, gen1, gen2, orb, set, 
          rep, stb, pnt, i, img, sch, new, g, g1, g2, t1, t2, s, t, R;

    # get the intersections of the point groups and the translation groups
    G  := Parent( G1 );
    d  := G.dimension - 1;
    P1 := Subgroup( PointGroup( G ), PointGroup( G1 ).generators );    
    P2 := Subgroup( PointGroup( G ), PointGroup( G2 ).generators );
    P  := Intersection( P1, P2 );
    T1 := TranslationsCrystGroup( G1 );
    T2 := TranslationsCrystGroup( G2 );
    T  := IntersectionModule( T1, T2 );
    L  := UnionModule( T1, T2 );

    gen  := P.generators;
    gen1 := List( gen, x -> PreImagesRepresentative( G1.pointHomom, x ) );
    gen2 := List( gen, x -> PreImagesRepresentative( G2.pointHomom, x )^-1 );

    orb := [ G.identity ];
    set := [ G.identity ];
    rep := [ P.identity ];
    stb := Subgroup( P, [] );

    # get the subgroup of P that can be lifted to the intersection
    for pnt  in orb  do
        for i in [1..Length( gen )] do
            img := gen2[i]*pnt*gen1[i];
            img[d+1]{[1..d]} := VectorModL( img[d+1]{[1..d]}, L ); 
            if not img in set  then
                Add( orb, img );
                AddSet( set, img );
                Add( rep, rep[Position(orb,pnt)]*gen[i] );
            else
                sch := rep[Position(orb,pnt)]*gen[i]
                       / rep[Position(orb,img)];
                if not sch in stb  then
                    stb := stb.operations.Closure( stb, sch );
                fi;
            fi;
        od;
    od;

    # determine the lift of stb
    new := [];
    for g in stb.generators do
        g1 := PreImagesRepresentative( G1.pointHomom, g );
        g2 := PreImagesRepresentative( G2.pointHomom, g );
        t1 := g1[d+1]{[1..d]};
        t2 := g2[d+1]{[1..d]};
        s  := IntSolutionMat( Concatenation( T1, -T2 ), t2-t1 ); 
        g1[d+1]{[1..d]} := t1+s{[1..Length(T1)]}*T1;
        Add( new, g1 );
    od;

    # add the translations
    for t in T do
        g1 := IdentityMat( d+1 );
        g1[d+1]{[1..d]} := t;
        Add( new, g1 );
    od;

    R := Subgroup( G, new );
    AddTranslationsCrystGroup( R, T );
    return R;

end;


#############################################################################
##
#F  CentralizerElement
##
CentralizerElement := function( G, u, TT )

    local d, P, T, I, L, U, orb, set, rep, stb, pnt, gen, img, sch, v;

    d := G.dimension - 1;
    P := PointGroup( G );
    T := TranslationsCrystGroup( G );
    I := IdentityMat( d );

    L := Concatenation( List( P.generators, x -> T*(x - I) ) );
    if Length( L ) > 0 then
        L := ReducedTranslationBasis( L );
    fi;

    U := TT*(u{[1..d]}{[1..d]} - I);
   
    orb := [ u ];
    set := [ u ];
    rep := [ G.identity ];
    stb := Subgroup( Parent( G ), [] );
    for pnt  in orb  do
        for gen  in G.generators  do
            img := pnt^gen;
            # reduce image mod L
            img[d+1]{[1..d]} := VectorModL( img[d+1]{[1..d]}, L );
            if not img in set  then
                Add( orb, img );
                AddSet( set, img );
                Add( rep, rep[Position(orb,pnt)]*gen );
            else
                sch := rep[Position(orb,pnt)]*gen
                       / rep[Position(orb,img)];
                # check if a translation conjugate of sch is in stabilizer
                v := u^sch - u;
                v := v[d+1]{[1..d]};
                v := IntSolutionMat( U, v );
                if v <> false then
                    sch[d+1]{[1..d]} :=sch[d+1]{[1..d]} - v*U; 
                    if not sch in stb  then
                        stb := stb.operations.Closure( stb, sch );
                    fi;
                fi;
            fi;
        od;
    od;
    return stb;
end;


#############################################################################
##
#F  CrystGroupOps.Centralizer( G, obj ) . . . centralizer of subgroup/element
##
CrystGroupOps.Centralizer := function ( G, obj )

    local d, P, T, e, I, M, m, L, i, U, o, gen, Q, C, u;

    d := G.dimension - 1;
    P := PointGroup( G );
    T := TranslationsCrystGroup( G );
    e := Length( T );
    I := IdentityMat( d );

    # we first determine the subgroup of G that centralizes the 
    # point group and the translation group of obj or its span

    if IsGroup( obj ) then
        M := PointGroup( obj );
        L := List( [1..e], x -> [] );
        for i in [ 1..Length( M.generators ) ] do
            L{[1..e]}{[1..d]+(i-1)*d} := T*(M.generators[i] - I);
        od;
        P := Centralizer( P, Subgroup( P, M.generators ) );
        P := Stabilizer( P, TranslationsCrystGroup( obj ), OnRight );
        U := Filtered( obj.generators, x -> x{[1..d]}{[1..d]} <> I );
    else
        M := obj{[1..d]}{[1..d]};
        L := T*(M - I);
        P := Centralizer( P, M );
        o := OrderMat( M );
        m := obj^o;
        P := Stabilizer( P, m[d+1]{[1..d]} );
        if o > 1 then U := [ obj ]; else U := []; fi;
    fi; 

    gen := List( P.generators, 
                 x -> PreImagesRepresentative( G.pointHomom, x ) );

    # if G is finite
    if e = 0 then
        return G.operations.Subgroup( Parent( G ), gen );
    fi;
    
    # otherwise, keep only translation generators which centralize obj
    Q := IdentityMat( e );
    if L<>[] then
        L := RowEchelonFormT( L, Q );
    fi;
    for i in [ Length( L )+1..e ] do
        M := IdentityMat( d+1 );
        M[d+1]{[1..d]} := Q[i]*T;
        Add( gen, M );
    od;

    # C centralizes the point group and the translation group of obj
    C := G.operations.Subgroup( Parent( G ), gen );     

    # now find the centralizer for each u in U
    for u in U do
        C := CentralizerElement( C, u, T );
    od;

    return C;

end;


#############################################################################
##
#F  TranslationNormalizer( S ) . . . . . . . . . . . translational normalizer
##
TranslationNormalizer := function( S )
    if not IsCrystGroup( S ) then
        Error("S must be a CrystGroup");
    fi;
    if not IsBound( S.translationNormalizer ) then
        S.translationNormalizer :=
            S.operations.TranslationNormalizer( S );
    fi;
    return S.translationNormalizer;
end;


#############################################################################
##
#F  CrystGroupOps.TranslationNormalizer( S ) . . . . . translation normalizer
##
CrystGroupOps.TranslationNormalizer := function( S )

    local P, T, d, N, M, I, g, i, Q, L, K, k, l, j, mul, gen;

    # do we have a space group?
    P := PointGroup( S );
    T := TranslationsCrystGroup( S );
    if not IsSpaceGroup( S ) then
        Error("S must be a space group");
    fi;

    d := S.dimension - 1;

    if Length( P.internalGenerators ) = 0 then
        N := Group( [], IdentityMat( d+1 ) );
        N.continuousTranslations := IdentityMat( d );
        return N;
    fi;

    M := List( [1..d], i->[] ); i := 0;
    I := IdentityMat( d );
    for g in P.internalGenerators do
        g := g - I;
        M{[1..d]}{[1..d]+i*d} := g;
        i := i+1;
    od;
    
    # first diagonalize M
    Q := IdentityMat( Length(M) );
    M := TransposedMat(M);
    M := RowEchelonForm( M );
    while not IsDiagonalMat(M) do
        M := TransposedMat(M);
        M := RowEchelonFormT(M,Q);
        if not IsDiagonalMat(M) then
            M := TransposedMat(M);
            M := RowEchelonForm(M);
        fi;
    od;

    # and then determine the solutions of x*M=0 mod Z
    if Length(M)>0 then
        L := List( [1..Length(M)], i -> [ 0 .. M[i][i]-1 ] / M[i][i] );
        L := List( Cartesian( L ), l -> l * Q{[1..Length(M)]} );
    else
        L := NullMat( 1, Length(Q) );
    fi;

    # get the kernel
    if Length(M) < Length(Q) then
        K := Q{[Length(M)+1..Length(Q)]};
        TriangulizeMat( K );
    else
        K := [];
    fi; 

    # reduce to basis modulo kernel
    Append( L, IdentityMat( d ) );
    for k in K do
        j := PositionProperty( k, x -> x=1 );
        for l in L do
            l := l-l[j]*k;
        od;
    od;
    mul := Lcm( Integers, List( Flat(L), Denominator ) );
    L   := HermiteNormalForm( L*mul ) / mul;

    # conjugate if not standard
    if not S.isStandard then
        L := L*T;
    fi;

    # get generators
    gen := List( L, x -> IdentityMat( d+1 ) );
    for i in [1..Length(L)] do
        gen[i][d+1]{[1..d]} := L[i];
    od;

    N := Group( gen, IdentityMat( d+1 ) );
    N.continuousTranslations := K;

    return N;

end;


#############################################################################
##
#F  AffineLift( pnt, d )
##
AffineLift := function( pnt, d )

    local M, b, i, I, p, m, Q, j, s; 

    M := List( [1..d], i->[] ); b := []; i := 0;
    I := IdentityMat( d );
    for p in pnt do
        m := p[1]{[1..d]}{[1..d]} - I;
        M{[1..d]}{[1..d]+i*d} := m;
        Append( b, p[2] - p[1][d+1]{[1..d]} );
        i := i+1;
    od;

    Q := IdentityMat( d );
    
    M := TransposedMat(M);
    M := RowEchelonFormVector( M,b );
    while not IsDiagonalMat(M) do
        M := TransposedMat(M);
        M := RowEchelonFormT(M,Q);
        if not IsDiagonalMat(M) then
            M := TransposedMat(M);
            M := RowEchelonFormVector(M,b);
        fi;
    od;

    ##  Check if we have any solutions modulo Z.
    for j in [Length(M)+1..Length(b)] do
        if not IsInt( b[j] ) then
            return [];
        fi;
    od;
    s := List( [1..Length(M)], i -> b[i]/M[i][i] );
    for i in [Length(M)+1..d] do 
        Add( s, 0);
    od;
    return s*Q;

end;


#############################################################################
##
#F  AffineNormalizer( S ) . . . . . . . . . . . . . . . . . affine normalizer
##
AffineNormalizer := function( S )
    if not IsCrystGroup( S ) then
        Error("S must be a CrystGroup");
    fi;
    if not IsBound( S.affineNormalizer ) then
        S.affineNormalizer := S.operations.AffineNormalizer( S );
    fi;
    return S.affineNormalizer;
end;


#############################################################################
##
#F  CrystGroupOps.AffineNormalizer( S ) . . . . . . . . . . affine normalizer
##
CrystGroupOps.AffineNormalizer := function( S )

    local d, P, T, N, gens, Pi, Si, hom, opr, orb, g, m, set, rep, 
          lst, pnt, img, t, sch, n, nn, normgens, TN, AN;

    d := S.dimension - 1;
    P := PointGroup( S );
    T := TranslationsCrystGroup( S );
    if not IsSpaceGroup( S ) then
        Error("S must be a space group");
    fi;
    N := NormalizerGL( P ); 

    # we work in a standard representation
    if not S.isStandard then
        gens := List( N.generators, x -> x^(T^-1) );
    else
        gens := N.generators;
    fi;
    Pi  := Group( P.internalGenerators, P.identity );
    Si  := Group( S.internalGenerators, S.identity );
    hom := GroupHomomorphismByImages( Si, Pi, 
                       Si.generators, Pi.generators );

    # the operation we shall need in the stabilizer algorithm
    opr := function( data, g )
        local m, mm, res; 
        m  := data[1]{[1..d]}{[1..d]};
        mm := m^g;
        if m = mm then
            res := [ data[1], List( data[2]*g, FractionModOne ) ];
        else
           m := PreImagesRepresentative( hom, mm );
           m[d+1]{[1..d]} := List( m[d+1]{[1..d]}, FractionModOne );
           res := [ m, List( data[2]*g, FractionModOne ) ];
        fi;
        return res;
    end;

    orb := [];
    for g in Si.generators do
        m := Copy( g );
        m[d+1]{[1..d]} := List( m[d+1]{[1..d]}, FractionModOne );
        Add( orb, [ m, m[d+1]{[1..d]} ] );
    od;
    orb := [ orb ];
    set := ShallowCopy( orb );

    rep := [ N.identity ];
    lst := [];
    for pnt  in orb  do
        for g  in gens  do
            img := List( pnt, x -> opr( x, g ) );
            if not img in set  then
                Add( orb, img );
                AddSet( set, img );
                Add( rep, rep[Position(orb,pnt)]*g );
            else
                t := AffineLift( img, d );
                if t<>[] then
                    sch := rep[Position(orb,pnt)]*g;
                    n := IdentityMat( d+1 );
                    n{[1..d]}{[1..d]} := sch;
                    n[d+1]{[1..d]} := t;
                    AddSet( lst, n );
                fi;
            fi;
        od;
    od;

    if IsBound( N.isFinite ) and N.isFinite or
       IsBound( N.size ) and N.size <> "infinity" then
        nn := Subgroup( N, [] );
        normgens := [];
        for g in lst do
            m := g{[1..d]}{[1..d]};
            if not m in nn then
                Add( normgens, g );
                nn := nn.operations.Closure( nn, m );
            fi;
        od;
    else
        normgens := lst;
    fi;

    m := IdentityMat( d+1 );
    m{[1..d]}{[1..d]} := T;
    if not S.isStandard then
        normgens := List( normgens, x -> x^m );
    fi;
    
    TN := TranslationNormalizer( S );
    Append( normgens, TN.generators );
    AN := Group( normgens, S.identity );
    AN.continuousTranslations := TN.continuousTranslations;

    return AN;

end;


#############################################################################
##
#F  AffineInequivalentSubgroups( sub )  reps of affine inequivalent subgroups
##
AffineInequivalentSubgroups := function( sub )

    local S, C, A, opr, reps, orb, pnt, gen, img;

    if sub = [] then
        return sub;
    fi;
    S := Parent( sub[1] );
    if not IsSpaceGroup( S ) then
        Error("parent of groups in sub must be a space group");
    fi;
    C := ShallowCopy( sub );
    A := AffineNormalizer( S );

    reps := [];
    while C <> []  do
        orb := [ C[1] ];
        for pnt  in orb  do
            for gen  in A.generators  do
                img := S.operations.ConjugateSubgroup( pnt, gen);
                if not img in orb then
                    Add( orb, img );
                fi;
            od;
        od;
        Add( reps, orb[1] );
        C := Filtered( C, x -> not x in orb );
    od;

    return reps;

end;






























