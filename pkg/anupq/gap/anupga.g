#############################################################################
##
#A  anupga.g                    GAP share library                Frank Celler
#A                                                           & Eamonn O'Brien
#A                                                           & Benedikt Rothe
##
#A  @(#)$Id$
##
#Y  Copyright 1992-1994,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
#Y  Copyright 1992-1994,  School of Mathematical Sciences, ANU,     Australia
##
#H  $Log$
##
if IsBound(ANUPQgroups)  then
    ThisIsAHack := ANUPQgroups;
else
    ThisIsAHack := [];
fi;
ANUPQgroups := ThisIsAHack;

if IsBound(ANUPQautos)  then
    ThisIsAHack := ANUPQautos;
else
    ThisIsAHack := [];
fi;
ANUPQautos := ThisIsAHack;

if IsBound(ANUPQtmpDir)  then
    ThisIsAHack := ANUPQtmpDir;
else
    ThisIsAHack := "ThisIsAHack";
fi;
ANUPQtmpDir := ThisIsAHack;

if IsBound(ANUPQmagic)  then
    ThisIsAHack := ANUPQmagic;
else
    ThisIsAHack := "ThisIsAHack";
fi;
ANUPQmagic := ThisIsAHack;


#############################################################################
##
#F  InfoANUPQ1  . . . . . . . . . . . . . . . . . . . . . . debug information
#F  InfoANUPQ2  . . . . . . . . . . . . . . . . . . . . . . debug information
##
if not IsBound(InfoANUPQ1)  then InfoANUPQ1 := Print;    fi;
if not IsBound(InfoANUPQ2)  then InfoANUPQ2 := Ignore;   fi;


#############################################################################
##
#F  DimensionFrattiniFactor( <G>, <p> ) . . . . . . . rank of <G> / <G>^p<G>'
##
DimensionFrattiniFactor := function( G, p )
    if not IsBound(G.rank)  then
        G.rank := LogInt( Index( G, PRump(G,p) ), p );
    fi;
    return G.rank;
end;


#############################################################################
##
#F  ANUPQauto( <G>, <gens>, <imgs> )  . . . . . . . .  construct automorphism
##
# ANUPQauto := function( G, g, i )
#     local   f;
# 
#     # if <frattineImages> is unbound construct this images
#     if not IsBound(G.frattiniGenerators) or G.frattiniGenerators <> g  then
#         G.abstractGens := WordList( Length(g), "a" );
#         G.abstractIgs  := AbstractIgs( G, g, G.abstractGens ).abstractIgs;
#         G.frattiniGenerators := g;
#     fi;
#     g := Igs(G);
#     i := List( G.abstractIgs, x -> MappedWord( x, G.abstractGens, i ) );
#     f := GroupHomomorphismByImages(G,G,g,i);
#     f.isHomomorphism := true;
#     f.isGroupHomomorphism := true;
#     f.isMapping := true;
#     f.isAutomorphism := true;
#     return f;
# end;

ANUPQauto := function( G, g, i )
   local   f;
   f := GroupHomomorphismByImages(G,G,g,i);
   f.isHomomorphism := true;
   f.isGroupHomomorphism := true;
   f.isMapping := true;
   f.isAutomorphism := true;
   return f;
end;


#############################################################################
##
#F  ANUPQautoList( <G>, <gens>, <L> ) . . . . . . . construct a list of autos
##
ANUPQautoList := function( G, gens, L )
    local   D,  i,  d,  l,  igs;

    # construct direct product elements
    D := [];
    for i  in [ 1 .. Length(gens) ]  do
	d := [];
	for l  in L  do
	    Add( d, l[i] );
	od;
	Add( D, DirectProductElement( d ) );
    od;

    # and compute the abstract igs simultaneously
    igs := AbstractIgs( G, gens, D );

    # construct the images list
    D := [];
    for i  in [ 1 .. Length(L) ]  do
	D[i] := [];
	for l  in igs.abstractIgs  do
	    Add( D[i], l.element[i] );
	od;
    od;

    # and then the automorphisms
    return List( D, x -> ANUPQauto( G, igs.igs, x ) );

end;


#############################################################################
##
#F  SavePqList( <file>, <lst> ) . . . . . . . . .  save a list of descendants
##
SavePqList := function( arg )
    local   file,  list,  appendExp,  entries,  entry,  l,  G,  p,
            i,  w,  str,  word,  j,  r,  gens;

    # check arguments
    if 3 < Length(arg) or Length(arg) < 2  then
        Error( "usage: SavePqList( <file>, <list> )" );
    fi;
    file := arg[1];
    list := arg[2];

    # function to add exponent vector
    appendExp := function( str, word )
        local   first, s, oldLen, i, w;

        first  := true;
        s      := str;
        oldLen := 0;
        for i  in [ 1 .. Length (word) ]  do
            if word[i] <> 0 then
                w := ConcatenationString( "G.", String (i) );
                if word[i] <> 1  then
                    w := ConcatenationString( w, "^", String(word[i]) );
                fi;
                if not first  then
                    w := ConcatenationString( "*", w );
                fi;
                if Length(s)+Length(w)-oldLen >= 77  then
                    s := ConcatenationString( s, "\n" );
                    oldLen := Length(s);
                fi;
                s := ConcatenationString( s, w );
                first := false;
            fi;
        od;
        if first  then
            s := ConcatenationString( s, "G.1^0" );
        fi;
        return s;
    end;



    # <entries> hold a list of record entries which will be saved
    entries := [ "pqIdent", "exponentLaw", "isAgAutomorphisms",
                 "isCapable", "rank", "nuclearRank" ];
    if Length(arg) = 3  then
        if IsList(arg[3])  then
            Append( entries, arg[3] );
       else
            Add( entries, arg[3] );
       fi;
    fi;

    # print head of file
    PrintTo(  file, "ANUPQgroups := [];\n"    );
    AppendTo( file, "Unbind(ANUPQautos);\n\n" );

    # run through all groups in <list>
    for l  in [ 1 .. Length(list) ]  do
        G    := list[l];
        gens := G.generators;
        p    := RelativeOrder(gens[1]);
        AppendTo( file, "## group number: ", l, "\n"                     );
        AppendTo( file, "ANUPQgroups[", l, "] := function( L )\n"        );
        AppendTo( file, "local   G,  H,  A,  B;\n"                       );
        AppendTo( file, "G := FreeGroup( ", Length(gens), ", \"G\" );\n" );
        AppendTo( file, "G.relators := [\n"                              );

        # at first the power relators
        for i in [ 1 .. Length(gens) ]  do
            if 1 < i  then
                AppendTo( file, ",\n" );
            fi;
            w   := gens[i]^p;
            str := ConcatenationString( "G.", String(i), "^", String(p) );
            if w <> G.identity  then
                word := Exponents( G, w );
                str  := ConcatenationString( str, "/(" );
                str  := appendExp( str,word );
                str  := ConcatenationString( str, ")" );
            fi;
            AppendTo( file, str );
        od;

        # and now the commutator relators
        for i  in [ 1 .. Length(gens)-1 ]  do
            for j  in [ i+1 .. Length(gens) ]  do
                w := Comm( gens[j], gens[i] );
                if w <> G.identity  then
                    word := Exponents( G, w );
                    str  := ConcatenationString(
                                ",\nComm( G.", String(j),
                                ", G.", String(i), " )/(" );
                    str := appendExp( str, word );
                    AppendTo( file, str, ")" );
                fi;
            od;
        od;
        AppendTo( file, "];\n" );

        # convert group into an ag group, save presentation
        AppendTo( file, "H := G;\n"                                );
        AppendTo( file, "G := AgGroupFpGroup(G);\n"                );
        AppendTo( file, "G.relators := H.relators;\n"              );
        AppendTo( file, "G.abstractGenerators := H.generators;\n"  );

        # add automorphisms
        if IsBound(G.automorphisms)  then
            DimensionFrattiniFactor( G, p );
            AppendTo( file, "A := [];\nB := [" );
    	    for r  in [ 1 .. G.rank ]  do
                AppendTo( file, "G.", r );
                if r <> G.rank  then
                    AppendTo( file, ", " );
    	    	else
    	    	    AppendTo( file, "];\n" );
                fi;
            od;
            for j  in [ 1 .. Length (G.automorphisms) ]  do
                AppendTo( file, "A[", j, "] := [");
                for r  in [ 1 .. G.rank ]  do
                    word := Image( G.automorphisms[j], gens[r] );
                    word := ExponentsAgWord( word );
                    AppendTo( file, appendExp( "", word ) );
                    if r <> G.rank  then
                        AppendTo (file, ", \n");
                    fi;
                od;
                AppendTo( file, "]; \n");
            od;
    	    AppendTo( file, "G.automorphisms := ANUPQautoList( G, B, A );\n" );
        fi;

        # add entries stored in <entries>
        for entry  in entries  do
            if IsBound( G.(entry) )  then
                AppendTo( file, "G.", entry, " := " );
                if IsString(G.(entry))  then
                    AppendTo( file, "\"", G.(entry), "\"" );
                else
                    AppendTo( file, G.(entry) );
                fi;
                AppendTo( file, ";\n" );
            fi;
        od;
        AppendTo( file, "Add( L, G );\n" );
        AppendTo( file, "end;\n\n\n"     );
    od;

    # write a magic string to the files
    AppendTo( file, "ANUPQmagic := \"groups saved to file\";\n" );
end;


#############################################################################
##
#F  PqList( <file> ) . . . . . . . . . . . . . . .  get a list of descendants
##
PqList := function( arg )
    local   lst,  help1,  help2,  sublist,  func;

    # check arguments
    if 2 < Length(arg) or Length(arg) < 1  then
        Error( "usage: PqList( <file> )" );
    fi;

    # save <ANUPQgroups> in case somebody has defined this variable
    if IsBound(ANUPQgroups)  then help1 := ANUPQgroups;  fi;
    Unbind( ANUPQgroups );
   
    # save <ANUPQautos> in case somebody has defined this variable
    if IsBound(ANUPQautos)  then help2 := ANUPQautos;  fi;
    Unbind( ANUPQautos );
   
    # try to read <file>
    Unbind(ANUPQmagic);
    if not READ(arg[1]) or not IsBound(ANUPQmagic)  then
        Unbind(ANUPQgroups);
        Unbind(ANUPQautos);
        if IsBound(help1)  then ANUPQgroups := help1;  fi;
        if IsBound(help2)  then ANUPQautos  := help2;  fi;
        return false;
    fi;

    # <lst> will hold the groups
    lst := [];
    if IsBound(ANUPQgroups) then
        if Length(arg) = 2  then
            if IsList(arg[2])  then
                sublist := arg[2];
            else
                sublist := [arg[2]];
            fi;
        else
            sublist := [ 1 .. Length(ANUPQgroups) ];
        fi;
        for func  in sublist  do
            ANUPQgroups[func](lst);
            if IsBound(ANUPQautos) and IsBound(ANUPQautos[func])  then
                ANUPQautos[func](lst[Length(lst)]);
            fi;
        od;
    fi;

    # restore <ANUPQgroups>
    Unbind(ANUPQgroups);
    Unbind(ANUPQautos);
    if IsBound(help1)  then ANUPQgroups := help1;  fi;
    if IsBound(help2)  then ANUPQautos  := help2;  fi;

    # return the groups
    return lst;

end;


#############################################################################
##
#F  ANUPQprintExps( <pqi>, <lst> ) . . . . . . . . . . .  print exponent list
##
ANUPQprintExps := function( pqi, lst )
    local   first,  l,  j;

    l := Length(lst);
    first := true;
    for j  in [1 .. l]  do
        if lst[j] <> 0  then
          if not first  then
              AppendTo( pqi, "*" );
          fi;
          first := false;
          AppendTo( pqi, "g", j, "^", lst[j] );
        fi;
    od;
end;


#############################################################################
##
#F  ANUPQerror( <param> ) . . . . . . . . . . . . . .report illegal parameter
##
ANUPQerror := function( param )
    Error(
    "Valid Options:\n",
    "    \"ClassBound\", <bound>\n",
    "    \"OrderBound\", <order>\n",
    "    \"StepSize\", <size>\n",
    "    \"AgAutomorphisms\"\n",
    "    \"RankInitialSegmentSubgroups\", <rank>\n",
    "    \"SpaceEfficient\"\n",
    "    \"AllDescendants\"\n",
    "    \"Exponent\", <exponent>\n",
    "    \"Metabelian\"\n",
    "    \"TmpDir\"\n",
    "    \"Verbose\"\n",
    "    \"SetupFile\", <file>\n",
    "Illegal Parameter: \"", param, "\"" );
end;


#############################################################################
##
#F  ANUPQextractArgs( <args>) . . . . . . . . . . . . . . parse argument list
##
ANUPQextractArgs := function( args )
    local   CR,  i,  act,  G,  match;

    # allow to give only a prefix
    match := function( g, w )
    	return 1 < Length(g) and SubString(w,1,Length(g)) = g;
    end;

    # extract arguments
    G  := args[1];
    CR := rec( group := G );
    i  := 2;
    while i <= Length(args)  do
        act := args[i];

        # "ClassBound", <class>
        if match( act, "ClassBound" )  then
            i := i + 1;
            CR.ClassBound := args[i];
            if CR.ClassBound <= G.pClass  then
                Error( "\"ClassBound\" must be at least ", G.pClass+1 );
            fi;

        # "OrderBound", <order>
        elif match( act, "OrderBound" )  then
            i := i + 1;
            CR.OrderBound := args[i];

        # "StepSize", <size>
        elif match( act, "StepSize" )  then
            i := i + 1;
            CR.StepSize := args[i];

        # "AgAutomorphisms"
        elif match( act, "AgAutomorphisms" )  then
            CR.AgAutomorphisms := true;

        # "RankInitialSegmentSubgroups", <rank>
        elif match( act, "RankInitialSegmentSubgroups" )  then
            i := i + 1;
            CR.RankInitialSegmentSubgroups := args[i];

        # "SpaceEfficient"
        elif match( act, "SpaceEfficient" ) then
            CR.SpaceEfficient := true;

        # "AllDescendants"
        elif match( act, "AllDescendants" )  then
            CR.AllDescendants := true;

        # "Exponent", <exp>
        elif match( act, "Exponent" )  then
            i := i + 1;
            CR.Exponent := args[i];

        # "Metabelian"
        elif match( act, "Metabelian" ) then
            CR.Metabelian := true;

        # "Verbose"
        elif match( act, "Verbose" )  then
            CR.Verbose := true;

        # "SubList"
        elif match( act, "SubList" )  then
            i := i + 1;
            CR.SubList := args[i];

        # temporary directory
        elif match( act, "TmpDir" )  then
            i := i + 1;
            CR.TmpDir := args[i];

        # "SetupFile", <file>
        elif match( act, "SetupFile" )  then
            i := i + 1;
            CR.SetupFile := args[i];

        # signal an error
        else
            return act;
        fi;
        i := i + 1;
    od;
    return CR;

end;


#############################################################################
##
#F  ANUPQinstructions( <pqi>, <param>, <p> )  . . . . . .  construct PQ input
##
ANUPQinstructions := function( pqi, param, p )
    local   G, firstStep,  ANUPQbool,  CR,  f,  i,  RF1,  RF2;

    G := param.group;
    firstStep := 0;

    # function to print boolean value
    ANUPQbool := function(b)
        if b  then
            return 1;
        else
            return 0;
        fi;
    end;

    # <CR> will hold the parameters
    CR := rec ();
    CR.ClassBound                  := G.pClass + 1;
    CR.StepSize                    := -1;
    CR.OrderBound                  := -1;
    CR.AgAutomorphisms             := IsBound(G.isAgAutomorphisms)
                                      and G.isAgAutomorphisms;
    CR.Verbose                     := false;
    CR.RankInitialSegmentSubgroups := 0;
    CR.SpaceEfficient              := false;
    CR.AllDescendants              := false;
    CR.ExpSet                      := false;
    CR.Exponent                    := 0;
    CR.Metabelian                  := false;
    CR.group                       := G;
    CR.SubList                     := -1; 

    # merge arguments with default parameters
    RF1 := RecFields(CR);
    RF2 := RecFields(param);
    for f  in RF2  do
        if not f in RF1  then
            ANUPQerror(f);
        else
            CR.(f) := param.(f);
        fi;
    od;

    # sanity check
    if CR.OrderBound <> -1 and CR.OrderBound <= LogInt(Size(G), p)  then
        return [false];
    fi;
    if CR.SpaceEfficient and not CR.AgAutomorphisms  then
        f := "\"SpaceEfficient\" is only allowed in conjunction with ";
        f := ConcatenationString( f, "\"AgAutomorphisms\"" );
        return f;
    fi;
    if CR.StepSize <> -1 and CR.OrderBound <> -1  then
        f := "\"StepSize\" and \"OrderBound\" must not be set ";
        f := ConcatenationString( f, "simultaneously" );
        return f;
    fi;

    # generate instructions
    AppendTo( pqi, "5\n", CR.ClassBound, "\n" );
    if CR.StepSize <> -1  then
        AppendTo( pqi, "0\n" );
        if CR.ClassBound = G.pClass + 1  then
            if IsList(CR.StepSize)  then
                if Length(CR.StepSize) <> 1  then
                    return "Only one \"StepSize\" must be given";
                else
                    CR.StepSize := CR.StepSize[1];
                fi;
            fi;
            AppendTo( pqi, CR.StepSize, "\n" );
            firstStep := CR.StepSize;
        else
            if IsList(CR.StepSize)  then
                if Length (CR.StepSize) <> CR.ClassBound - G.pClass  then
                    f := "The difference between maximal class and class ";
                    f := ConcatenationString(f, 
                      "of the starting group is ", 
                      String(CR.ClassBound-G.pClass),
                      ".\nTherefore you must supply ",
                      String(CR.ClassBound-G.pClass),
                      " step-sizes in the \"StepSize\" list \n" );
                    return f;
                fi;
                AppendTo( pqi, "0\n" );
                for i  in CR.StepSize  do
                    AppendTo( pqi, i, " " );
                od;
                AppendTo( pqi, "\n" );
                firstStep := CR.StepSize[1];
            else
                AppendTo( pqi, "1\n", CR.StepSize, "\n" );
            fi;
        fi;
    elif CR.OrderBound <> -1  then
        AppendTo( pqi, "1\n1\n", CR.OrderBound, "\n" );
    else
        AppendTo( pqi, "1\n0\n" );
    fi;    
    AppendTo( pqi, ANUPQbool(CR.AgAutomorphisms), "\n0\n",
                   CR.RankInitialSegmentSubgroups, "\n" );
    if CR.AgAutomorphisms  then
        AppendTo( pqi, ANUPQbool(CR.SpaceEfficient), "\n" );
    fi;
    if IsBound(G.nuclearRank) and firstStep <> 0  and
        firstStep > G.nuclearRank then
            f := ConcatenationString( "\"StepSize\" (=", String(firstStep),
                   ") must be smaller or equal the \"Nuclear Rank\" (=",
                   String(G.nuclearRank), ")" );
            return f;
    fi;
    AppendTo( pqi, ANUPQbool(CR.AllDescendants), "\n" );
    AppendTo( pqi, CR.Exponent, "\n" );
    AppendTo( pqi, ANUPQbool(CR.Metabelian), "\n", "1\n" );

    # return success
    return [true, CR];

end;


#############################################################################
##
#F  PqDescendants( <G>, ... ) . . . . . . . . .  construct descendants of <G>
##
PqDescendants := function( arg )
    local   G,  lst,  dir,  pqi,  p,  l,  i,  j,  rank,  aut,  gens,
            res,  cmd,  desc,  CR, x, r;

    # check arguments, return if <G> is not capable
    if 0 = Length (arg)  then
        Error( "usage: PqDescendants( <G>, ... )" );
    fi;
    G := arg[1];
    if not IsAgGroup(G)  then
        Error( "<G> must be an Ag-Group" );
    fi;

    # get exponent-p class
    p := Order( G, G.generators[Length(G.generators)] );
    if not IsBound(G.pClass)  then
        G.pClass := Length( PCentralSeries( G, p ) ) - 1;
    fi;

    # extract arguments in case the second arg is not a argument record
    if Length(arg) < 2 or not IsRec(arg[2])  then
        CR := ANUPQextractArgs( arg );
        if IsString(CR)  then
            ANUPQerror(CR);
        fi;
    else
        CR := ShallowCopy(arg[2]);
        CR.group := G;
    fi;

    rank := DimensionFrattiniFactor( G, p );

    # if automorphisms are not supplied and group has p-class 1, 
    # construct automorphisms, else signal Error 
    if (not IsBound (G.automorphisms)) then 
        if (G.pClass = 1) then 
           G.automorphisms := [];
           for x  in GeneralLinearGroup(rank,p).generators  do
               aut := [];
               for i  in [ 1 .. rank ]  do
                   r := G.identity;
                   for j  in [ 1 .. rank ]  do
                       r := r * G.generators[j]^Int(x[i][j]);
                   od;
                   aut[i] := r;
               od;
               Add( G.automorphisms, GroupHomomorphismByImages( G, G,
                       G.generators, aut ) );
           od;
        else 
           Error ("<G> must have class 1 or <G>.automorphisms must be bound\n");
        fi;
    fi;

    # if <G> is not capable and we want to compute something return
    if     IsBound(G.isCapable)
       and not G.isCapable 
       and not IsBound(CR.SetupFile)
    then
        return [];
    fi;

    # we only want to set up an input file for ANU pq
    if IsBound(CR.SetupFile)  then
        pqi := CR.SetupFile;
        Unbind(CR.SetupFile);

    # otherwise construct a temporary directory
    elif not IsBound(CR.TmpDir) and ANUPQtmpDir = "ThisIsAHack"  then
        dir := TmpName();

        # create the directory
        Exec( ConcatenationString( "mkdir ", dir ) );
        pqi := ConcatenationString( dir, "/PQ_INPUT" );

    # use a given directory and try to construct a random subdir
    else
        if IsBound(CR.TmpDir)  then
            dir := CR.TmpDir;
            Unbind(CR.TmpDir);
        else
            dir := ANUPQtmpDir;
        fi;

        # try to get a random number
        i := Runtime();
        i := i + RandomList( [ 1 .. 2^16 ] ) * RandomList( [ 1 .. 2^16 ] );
        i := i * Runtime();
        i := i mod 19^8;
        dir := ConcatenationString( dir, "/", LetterInt(i), ".apq" );

        # create the directory
        Exec( ConcatenationString( "mkdir ", dir ) );
        pqi := ConcatenationString( dir, "/PQ_INPUT" );
    fi;

    # write first instruction to start ANU pq p-group generation
    PrintTo(  pqi, "1\n"                     );
    AppendTo( pqi, "prime ", p, " \n"        );
    AppendTo( pqi, "class ", G.pClass, " \n" );

    # print generators of <G>
    l := Length(G.generators);
    AppendTo( pqi, "generators {" );
    for i  in [ 1 .. l ]  do
        AppendTo( pqi, "g", i );
        if i <> l  then
            AppendTo( pqi, ", " );
        fi;
    od;
    AppendTo( pqi, " }\n" );

    # print relators of <G>
    AppendTo( pqi, "relations {" );
    for i  in [ 1 .. l ]  do
        if i <> 1  then
            AppendTo( pqi, ", " );
        fi;
        lst := ExponentsAgWord( G.generators[i]^p );
        AppendTo( pqi, "g", i, "^", p );
        if ForAny( lst, x -> x<>0 )  then
            AppendTo( pqi, "=" );
            ANUPQprintExps( pqi, lst );
        fi;
    od;
    for i  in [ 1 .. l ]  do
        for j  in [ 1 .. i-1 ]  do
            lst := ExponentsAgWord( Comm( G.generators[i],
                                          G.generators[j] ) );
            AppendTo( pqi, ", [g", i, ", g", j, "]" );
            if ForAny( lst, x -> x<>0 )  then
                AppendTo( pqi, "=" );
                ANUPQprintExps( pqi, lst );
            fi;
        od;
    od;
    AppendTo( pqi, "} \n" );

   
    # enter p-group generation
    AppendTo( pqi, "; \n7\n9\n1\n" );

    # print automorphisms of <G>
    AppendTo( pqi, Length(G.automorphisms), "\n" );
    for aut  in G.automorphisms  do
        for gens  in [ 1 .. rank ]  do
            for i in ExponentsAgWord(Image(aut, G.generators[gens]))  do
                AppendTo( pqi, i, " " );
            od;
            AppendTo( pqi, "\n" );
        od;
    od;

    # now construct the instruction from the args
    res := ANUPQinstructions( pqi, CR, p );
    if IsString(res)  then
        if IsBound(dir)  then
            Exec( ConcatenationString( "rm -rf ", dir ) );
        fi;
        Error(res);
    elif not res[1]  then
        if IsBound(dir)  then
            Exec( ConcatenationString( "rm -rf ", dir ) );
            return [];
        fi;
        return;
    fi;

    #the next two lines were added by EOB 
    CR := res[2];
    res := CR.Verbose;

    AppendTo( pqi, "0\n0\n" );

    # return if we only want to set up a input file
    if not IsBound(dir)  then
    	Print( "#I  input file '", pqi, "' written, ",
    	       "run 'pq' with '-k' flag\n" );
        return true;
    fi;

    # start PQ (verbose if <res> is true)
    cmd := ConcatenationString( "-v ",EXENAME,"gap.sh -k -g < ", pqi );
    if not res  then
        cmd := ConcatenationString( cmd, " > PQ_LOG" );
    fi;
    ExecPkg( "anupq", "bin/pq", cmd, dir );

    # read in the library file written by pq
    if CR.SubList <> -1 then 
       desc := PqList( ConcatenationString(dir,"/GAP_library"), CR.SubList );
    else 
       desc := PqList( ConcatenationString(dir,"/GAP_library") );
    fi;
    if desc = false  then
        Exec( ConcatenationString( "rm -rf ", dir ) );
        Error( "cannot execute ANU pq,  please check installation" );
    fi;

    # add 'isCapable'
    for G  in desc  do
        if not IsBound(G.isCapable)  then
            G.isCapable := false;
        fi;
    od;

    # remove temporary directory and return
    Exec( ConcatenationString( "rm -rf ", dir ) );
    return desc;

end;
