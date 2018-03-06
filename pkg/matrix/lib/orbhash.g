#############################################################################
##
#W  orbhash.g	 	   	Matrix Packages                  Frank Celler
##
#H  @(#)$Id: orbhash.g,v 1.1 1997/03/10 13:52:06 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This  file contains high   level  dispatcher functions for the  classical
##  group recognition.
##
RevisionMatrix.orbhash_g :=
    "@(#)$Id: orbhash.g,v 1.1 1997/03/10 13:52:06 gap Exp $";


#############################################################################
##

#F  InfoMatOrbitN . . . . . . . . . . . . . . . . . . . information functions
##
if not IsBound(InfoMatOrbit1)  then InfoMatOrbit1 := Ignore;  fi;
if not IsBound(InfoMatOrbit2)  then InfoMatOrbit2 := Ignore;  fi;


#############################################################################
##

#F  VectorHashTable . . . . . . . . . . . .  hash table functions for vectors
##
VectorHashTable := rec();


VectorHashTable.NumberVector1 := function( tab, v )
    local   n,  i,  z;

    n := 0;
    for i  in [ 1 .. tab.dimension ]  do
        z := v[i];
        if z <> tab.zero  then
            n := n + tab.powers[i]*(LogFFE(z,tab.root)+1);
        fi;
    od;
    return n;
end;


VectorHashTable.NumberVector2 := function( tab, v )
    return IntVecFFE(v) * tab.powers;
end;


VectorHashTable.ResizeTable := function( tab )
    local   elements,  i;

    #Print( "#I  resizing hash table: mod = ", tab.modulus, ", size = ",
    #       tab.size, ", missses = ", tab.misses, "\n" );
    elements     := tab.elements;
    tab.elements := [];
    tab.size     := 0;
    tab.modulus  := NextPrimeInt( 3*tab.modulus );
    for i  in elements  do
        VectorHashTable.InsertVector( tab, i );
    od;
end;


#############################################################################
##
#F  VectorHashTable.HashTable( <d>, <q>, <init> )
##
VectorHashTable.HashTable := function( d, q, init )
    local   tab;

    tab := rec();
    tab.modulus     := NextPrimeInt( Minimum( d^q, init ) );
    tab.root        := Z(q);
    tab.powers      := List( [ 0 .. d-1 ],  x -> q^x );
    tab.size        := 0;
    tab.elements    := [];
    tab.misses      := 0;
    tab.dimension   := d;
    tab.zero        := 0*tab.root;

    if IsPrimeInt(q)  then
        tab.hash := VectorHashTable.NumberVector2;
    else
        tab.hash := VectorHashTable.NumberVector1;
    fi;

    return tab;
end;


#############################################################################
##
#F  VectorHashTable.InsertVector( <tab>, <v> )
##
VectorHashTable.InsertVector := function( tab, v )
    local   val,  key;

    if not IsInt(v)  then
        val := tab.hash( tab, v );
    else
        val := v;
    fi;
    key := ( val mod tab.modulus ) + 1;
    if not IsBound(tab.elements[key])  then
        tab.elements[key] := val;
        tab.size := tab.size + 1;
        if tab.modulus < 3*tab.size  then
            VectorHashTable.ResizeTable(tab);
        fi;
        return true;
    fi;
    while IsBound(tab.elements[key]) and tab.elements[key] <> val  do
        key := ( key mod tab.modulus ) + 1;
        tab.misses := tab.misses + 1;
    od;
    if IsBound(tab.elements[key])  then
        return false;
    else
        tab.elements[key] := val;
        tab.size := tab.size + 1;
        if tab.modulus < 3*tab.size  then
            VectorHashTable.ResizeTable(tab);
        fi;
        return true;
    fi;
end;


#############################################################################
##

#F  OrbitMat( <G>, <v>, <base>, <max> )	. . . . . . . orbit for matrix groups
##
OrbitMat := function ( G, v, base, max )

    local   orb,        # orbit, result
            tab,        # orbit <orb> as hash table for fast membership test
            gen,        # generator of the group <G>
            pos,        # current position in the orbit <orb>
            img,        # image of the point <pnt> under the generator <gen>
            dim,        # dimension
            fld,        # field
            i,          # loop me
            p;          # position to enter new point

    orb := [ v ];
    dim := Length(G.identity);
    fld := G.field;
    tab := VectorHashTable.HashTable( dim, Size(fld), 2000 );
    VectorHashTable.InsertVector( tab, v );
    pos := 1;
    while pos <= Length(orb)  do
        for gen  in G.generators  do
            img := orb[pos]*gen;
            if VectorHashTable.InsertVector( tab, img )  then
                Add( orb, img );
                if max <= tab.size  then
                    return false;
                fi;
            fi;
        od;
        pos := pos + 1;
    od;

    # update the base (can be done faster!)
    for i  in orb  do
        if Length(base) < dim  then
            if Length(base) < RankMat(Concatenation(base,[i]))  then
                Add( base, i );
            fi;
        fi;
    od;

    # return the orbit <orb>
    return orb;

end;


#############################################################################
##
#F  NiceEigenVectorsMatGroup( <G> )
##
NiceEigenVectorsMatGroup := function( G )
    local   j,  i,  elms,  cpls,  ns,  int,  found;

    found := [];
    for  j  in [ 1 .. 5 ]  do

        # get random elements and their characteristic polynomials
        elms := List( [1..2], x -> PseudoRandom(G) );
        cpls := [];
        for i  in [ 1 .. Length(elms) ]  do
            cpls[i] := CharacteristicPolynomial(elms[i]);
            cpls[i].baseRing := G.field;
            cpls[i] := Filtered( Set(Factors(cpls[i])), x -> Degree(x)=1 );
        od;

        # compute the nullspace
        ns := [];
        for i  in [ 1 .. Length(cpls) ]  do
            if 0 < Length(cpls[i])  then
                ns[i] := List( cpls[i], x -> NullspaceMat(Value(x,elms[i])) );
                ns[i] := Concatenation(ns[i]);
            else
                ns[i] := [];
            fi;
        od;

        # check if there is a common nullspace
        int := IntersectionMat( ns[1], ns[2] );
        if 0 < Length(int) and Length(int) < 10  then
            InfoMatOrbit2( "#I  found common eigenspace of dimension ",
                           Length(int), "\n" );
            return int;

        elif 0 < Length(ns[1]) and Length(ns[1]) < 10  then
            InfoMatOrbit2( "#I  found single eigenspace of dimension ",
                           Length(ns[1]), "\n" );
            return ns[1];

        elif 0 < Length(ns[2]) and Length(ns[2]) < 10  then
            InfoMatOrbit2( "#I  found single eigenspace of dimension ",
                           Length(ns[2]), "\n" );
            return ns[2];

        elif 0 < Length(int)  then
            Add( found, Random(int) );

        elif 0 < Length(ns[1])  then
            Add( found, Random(ns[1]) );

        elif 0 < Length(ns[2])  then
            Add( found, Random(ns[2]) );
        fi;

    od;
    if 0 < Length(found)  then
        InfoMatOrbit2("#I  found ",Length(found)," random eigenvectors\n");
        return found;
    else
        return false;
    fi;

end;


#############################################################################
##
#F  PermGroupRepresentation( <G>, <max> )
##
PermGroupRepresentation := function( G, max )
    local   dim,  base,  fld,  orbs,  tries,  space,  ns,  v,  
            orb,  m,  f;

    if IsBound(G.permGroupP)  then
        return G.permGroupP;
    fi;

    # start with no basis vectors
    dim   := Length(G.identity);
    base  := [];
    fld   := G.field;
    orbs  := [];
    tries := 0;
    space := fld^dim;

    # try to find a base with small orbits
    while Length(base) < dim  do

        # use nullspaces
        ns := NiceEigenVectorsMatGroup(G);
        if ns = false  then
            InfoMatOrbit2( "#I  using random vector\n" );
            ns := [ Random(space) ]; 
        fi;
        ns := Difference( ns, orbs );
        InfoMatOrbit2( "#I  try next ", Length(ns), " vectors\n" );
        if 0 = Length(ns)  then
            tries := tries + 1;
        fi;
        if 5 < Length(ns)  then
            ns := ns{[1..5]};
        fi;

        # try all vectors
        for v  in ns  do
            if Length(base) < dim  then
                orb := OrbitMat( G, v, base, max );
                if orb <> false  then
                    InfoMatOrbit2( "#I  found orbit of length ", Length(orb),
                                   " and base ", Length(base), "\n" );
                    UniteSet( orbs, orb );
                else
                    InfoMatOrbit2( "#I  failed to find small orbit\n" );
                    tries := tries + 1;
                fi;
            fi;
        od;
        if 10 < tries  then
            InfoMatOrbit2( "#I  giving up after ", tries, " tries\n" );
            return false;
        fi;
    od;
    InfoMatOrbit2( "#I  found basis\n" );

    # construct the permutations
    m := [];
    for v  in G.generators  do
        Add( m, PermList(List(orbs,x->PositionSorted(orbs,x*v))) );
    od;
    f := Group( m, () );
    base := List( base, x -> PositionSorted( orbs, x ) );
    #This value was originally 200, but it failed for 2m11d5.gap -- EOB
    StabChain( f, rec( knownBase := base, random := 1000 ) );

    # and return
    G.permDomain := orbs;
    G.permGroupP := f;
    G.permGroupP.operation := rec();
    G.permGroupP.operation.genimages := m;
    return f;

end;


#############################################################################
##

#E  orbhash.g	. . . . . . . . . . . . . . . . . . . . . . . . . . ends here
##
