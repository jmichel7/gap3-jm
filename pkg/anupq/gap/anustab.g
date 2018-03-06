#############################################################################
##
#A  stabcalc.g                  GAP share library                Frank Celler
#A                                                           & Benedikt Rothe
##
#A  @(#)$Id$
##
#Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
#H  $Log$
##
if not IsBound(ANUPQ_field)  then ANUPQ_field := false;  fi;
if not IsBound(ANUPQ_s)      then ANUPQ_s     := false;  fi;
if not IsBound(ANUPQ_r)      then ANUPQ_r     := false;  fi;
if not IsBound(ANUPQ_d)      then ANUPQ_d     := false;  fi;
if not IsBound(ANUPQ_q)      then ANUPQ_q     := false;  fi;


#############################################################################
##
#F  InfoANUPQ1	. . . . . . . . . . . . . . . . . . . . . . debug information
#F  InfoANUPQ2	. . . . . . . . . . . . . . . . . . . . . . debug information
##
if not IsBound(InfoANUPQ1)  then InfoANUPQ1 := Ignore;   fi;
if not IsBound(InfoANUPQ2)  then InfoANUPQ2 := Ignore;   fi;


#############################################################################
##
#F  ANUPQpolycyclicGenerators( <G> ) . . . . .  polycyclic generators for <G>
##
ANUPQpolycyclicGenerators := function( G )
    local   P,  L;

    # convert <G> into a permutation group
    P := PermGroup(G);

    # if <P> is trivial return a empty set
    if 1 = Size(P)  then
    	return [];

    # use 'CompositionSeries' if <G> is solvable
    elif IsSolvable(P)  then
        if not IsBound(P.bsgs)  then
            Unbind(P.orbit);
            Unbind(P.transversal);
            Unbind(P.stabilizer);
            P.operations.TryElementaryAbelianSeries(P);
        fi;
    	L := List( P.bssgs, x -> Image(P.bijection,x) );
    	return Reversed(L);

    # otherwise give up
    else
    	return false;
    fi;

end;


#############################################################################
##
#F  ANUPQoutputResult( <G>, <name> ) . . . . . . . . . . . create File <name>
##
ANUPQoutputResult := function( G, name )
    local   dim,  gens,  mat,  row,  i,  log,  tab,  isStdout;

    # remove <name> if it exists
    if name <> "*stdout*"  then
        Exec( ConcatenationString( "rm -f ", name ) );
    fi;

    # if <G> is solvable get polycyclic generators
    isStdout := name = "*stdout*";
    gens     := ANUPQpolycyclicGenerators( G );
    if gens = false  then
    	if isStdout  then
    	    Print( "0\n" );
    	else
            PrintTo( name, "0\n" );
    	fi;
    	gens := G.generators;
    else
    	if isStdout  then
    	    Print( "-1\n" );
    	else
            PrintTo( name, "-1\n" );
    	fi;
    fi;
    if isStdout  then
    	Print( Length(gens), "\n" );
    else
        AppendTo( name, Length(gens), "\n" );
    fi;

    # print the generators in PQ format
    if 0 < Length(gens)  then
        dim := Length(gens[1][1]);
        tab := IntegerTable( ANUPQ_field );
        for mat  in gens  do
            for row  in mat  do
                for i  in [ 1 .. dim ]  do
                    log := LogVecFFE( row, i );
                    if log = false  then
    	    	    	if isStdout  then
    	    	    	    Print( " 0\n" );
    	    	    	else
                            AppendTo( name, " 0 " );
    	    	    	fi;
                    else
    	    	    	if isStdout  then
    	    	    	    Print( " ", tab[log+1], "\n" );
    	    	    	else
                            AppendTo( name, " ", tab[log+1], " " );
    	    	    	fi;
                    fi;
                od;
    	    	if not isStdout  then
    	            AppendTo( name, "\n" );
    	    	fi;
    	    od;
        od;
    fi;
end;


#############################################################################
##
#F  ANUPQmatrixLabel( <label>, <offs> )  . . . . compute matrix from <label>
##
ANUPQmatrixLabel := function( label, offs )
    local   pos,  defs,  mat,  row,  col,  start,  next,  char;

    # adjust <label> and find <defs>
    label := label - 1;
    pos   := PositionSorted( offs[2], label );
    if Length(offs[1]) < pos or offs[2][pos] <> label  then
        pos := pos - 1;
    fi;
    defs := offs[1][pos];

    # initialize matrix
    mat  := List( [ 1 .. ANUPQ_s ], x -> [1..ANUPQ_q] * ANUPQ_field.zero );
    next := 1;
    for col  in [ 1 .. ANUPQ_r ]  do
        if defs[col]  then
      	    mat[next][col] := ANUPQ_field.one;
      	    next := next + 1;
    	fi;
    od;

    # construct the matrix
    label := label - offs[2][pos];
    row   := 1;
    start := 1;
    char  := ANUPQ_field.char;
    while row <= ANUPQ_s  do
        while not defs[start]  do
            start := start + 1;
    	od;
    	start := start + 1;
    	for col  in [ start .. ANUPQ_q ]  do
      	    if col > ANUPQ_r or not defs[col]  then
            	mat[row][col] := ANUPQ_field.one * RemInt( label, char );
            	label := QuoInt( label, char );
      	    fi;
    	od;
    	row := row + 1;
    od;
    return mat;

end;

#############################################################################
##
#F  ANUPQlabelMatrix( <mat>, <offs> )  . . . . . .  compute label from <mat>
##
ANUPQlabelMatrix := function( mat, offs )
  local defs,label,row,startcol, col, pf, lst, tab, char, log;

    # construct <defs>
    defs := BlistList( [ 1 .. ANUPQ_r ], [] );
    col  := 1;
    for row  in [ 1 .. Length(mat) ]  do
    	while mat[row][col] = ANUPQ_field.zero  do
      	    col := col + 1;
    	od;
    	defs[col] := true;
    	col := col + 1;
    od;

    # construct <pf>
    pf   := 0;
    lst  := ANUPQ_r;
    row  := ANUPQ_s;
    tab  := IntegerTable( ANUPQ_field );
    char := ANUPQ_field.char;
    repeat
    	while not defs[lst]  do
      	    lst := lst - 1;
    	od;
    	col := ANUPQ_q;
    	while col > lst  do
      	    if col > ANUPQ_r or not defs[col]  then
    	    	log := LogVecFFE( mat[row], col );
    	    	if log = false  then
    	    	    pf := pf * char;
    	    	else
    	    	    pf := pf * char + tab[log+1];
    	    	fi;
    	    fi;
    	    col := col - 1;
    	od;
    	lst := lst - 1;
    	row := row - 1;
    until row < 1;

    # return the label
    return offs[2][PositionSorted(offs[1],defs)] + pf + 1;

end;


#############################################################################
##
#F  ANUPQonLabels( <mat>, <label>, <offs> )  . .  action of <mat> on <label>
##
ANUPQonLabels := function( mat, label, offs)
    mat := ANUPQmatrixLabel( label, offs ) * Transposed(mat);
    TriangulizeMat(mat);
    return ANUPQlabelMatrix( mat, offs );
end;


#############################################################################
##
#F  ANUfastProduct( <to>, <prev>, <mist> ) . . . . . . . . . fast subproduct
##
ANUfastProduct := function( to, prev, mlist )
    local   IDFU, DIDFU;

    IDFU  := RootInt( to, 3 );
    DIDFU := IDFU^3;
    if prev >= DIDFU  then
    	return false;
    else
    	return [ DIDFU+1, mlist[IDFU] ];
    fi;
end;


#############################################################################
##
#F  ANUPQstabilizer( <label>, <olen>, <glb> )  . . . . . . . . . . . .  main
##
ANUPQstabilizer := function( label, olen, glb )

    local   grpDp,  sizeStab,  cnst,  base,  offs,  i,  sum,  idQp,
    	    ntQp,  stab,  orbitTree,  word,  wordQp,  wordDp, mlist,
    	    count,  next, previous, earlier,  newStab,  newLabel, orb,
    	    speedUp,  invDp,  newStabQp,  mlistQp,  invQp,  genQp,  genDp;

    # set the global variables
    ANUPQ_field := glb.F;
    ANUPQ_s := glb.s;
    ANUPQ_r := glb.r;
    ANUPQ_d := glb.d;
    ANUPQ_q := glb.q;
    genDp   := glb.genD;
    genQp   := glb.genQ;

    # calculate stab size using the orbit length <olen>
    grpDp    := MatGroup( genDp, ANUPQ_field );
    sizeStab := Size( grpDp ) / olen;
    
    # calculate <offs>
    cnst := ANUPQ_q * ANUPQ_s - ANUPQ_s * (ANUPQ_s-1) / 2;
    base := BlistList( [1 .. ANUPQ_r ], [ 1 .. ANUPQ_s ] );
    offs := [ Arrangements( base, ANUPQ_r ), [0] ];
    for i  in [ 2 .. Length(offs[1]) ]  do
    	sum := Sum( ListBlist( [ 1 .. ANUPQ_r ], offs[1][i-1] ) );
    	offs[2][i] := offs[2][i-1] + ANUPQ_field.char^(cnst-sum);
    od;

    # construct a list of position of non-trivial <genQp>
    idQp := genQp[1]^0;
    ntQp := Filtered( [ 1 .. Length(genQp) ], x -> genQp[x] <> idQp );

    # check if any of the generators already fixes <label>
    stab := Subgroup( grpDp, List( Filtered( [ 1 .. Length(genDp) ],
    	    	x -> label = ANUPQonLabels( genQp[x], label, offs ) ),
    	    	y -> genDp[y] ) );

    # <orbitTree> contains the orbit as tree
    orbitTree := rec();

    # <word> is the current word as list of generator numbers
    word := [];

    # <wordQp> is the current word in <genQp>
    wordQp := IdentityMat( ANUPQ_q, ANUPQ_field.one );

    # <wordDp> is the current word in <genDp>
    wordDp := IdentityMat( ANUPQ_d, ANUPQ_field.one );
    
    # <mlist> will hold a few subproducts
    mlist   := [];
    #D mlistQp := [];

    # we are done as soon as we have a stabilizer of size <sizeStab>
    count := 1;
    InfoANUPQ1( "#I  |<stab>| = ",Size(stab),", needed = ",sizeStab,"\n" );

    #change by EOB -- September 1996 
    #while Size(stab) <> sizeStab  do

    while Size(stab) < sizeStab  do
    	next   := RandomList(ntQp);
        wordDp := wordDp * genDp[next];
        wordQp := wordQp * genQp[next];
    	if RootInt( Length(word)+1, 3 ) > Length(mlist)  then
      	    Add( mlist, wordDp );
    	    #D Add( mlistQp, wordQp );
    	fi;
    	newLabel  := ANUPQonLabels( wordQp, label, offs );
    	invDp     := wordDp ^ -1;
    	newStab   := invDp;
    	#D invQp     := wordQp ^ -1;
    	#D newStabQp := invQp;

    	# find earlier occurences of <newLabel>
    	orb := orbitTree;
    	while IsBound(orb.L) and orb.L <> newLabel  do
    	    if newLabel < orb.L  then
      	    	orb := orb.s;
    	    else
    	    	orb := orb.b;
    	    fi;
    	od;
    	if not IsBound(orb.L)  then
    	    orb.L := newLabel;
    	    orb.F := [];
    	    orb.s := rec();
    	    orb.b := rec();
    	fi;
    	orb := orb.F;

    	# add new stabilizer elements
    	previous := 1;
    	earlier  := 1;
        while Size(stab) <> sizeStab and earlier <= Length(orb)  do
    	    InfoANUPQ2( "#I  <word> = ", word, ", ",
    	    	    	"<orb>[<earlier>] = ", orb[earlier], ", ",
    	    	        "<previous> = ", previous, "\n" );
      	    speedUp := ANUfastProduct( orb[earlier], previous, mlist );
      	    if speedUp <> false  then
            	newStab := invDp
    	    	    	   * speedUp[2]
    	    	    	   * Product( [ speedUp[1] .. orb[earlier] ],
                                      i -> genDp[word[i]] );
      	    else
           	newStab := newStab
    	    	    	   * Product( [ previous .. orb[earlier] ],
                                      i -> genDp[word[i]] );
      	    fi;

            #D speedUp := ANUfastProduct( orb[earlier], previous, mlistQp );
            #D if speedUp <> false  then
            #D     newStabQp := invQp
            #D                  * speedUp[2]
            #D                  * Product( [ speedUp[1] .. orb[earlier] ],
            #D                             i -> genQp[word[i]] );
            #D else
            #D     newStabQp := newStabQp
            #D                  * Product( [ previous .. orb[earlier] ],
            #D                             i -> genQp[word[i]] );
            #D fi;
            #D if label <> ANUPQonLabels(newStabQp,label,offs)  then
            #D     Error( "wrong stabilizer element" );
            #D else
            #D     InfoANUPQ1( "#I  stabilizer calculation OK\n" );
            #D fi;

            if not newStab in stab  then
                stab := Closure( stab, newStab );
                InfoANUPQ1( "#I  |<stab>| = ", Size(stab), ", needed = ",
                            sizeStab, "\n" );
            fi;
      	    previous := orb[earlier] + 1;
      	    earlier  := earlier + 1;
    	od;
    	Add( orb, count );
    	count := count+1;
    	Add( word, next );
    od;
    return stab;

end;
