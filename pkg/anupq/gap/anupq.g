#############################################################################
##
#A  anupq.g                     GAP share library              Eamonn O'Brien
#A                                                             & Frank Celler
##
#A  @(#)$Id$
##
#Y  Copyright 1992-1994,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
#Y  Copyright 1992-1994,  School of Mathematical Sciences, ANU,     Australia
##
#H  $Log$
##
if not IsBound(G)  then Unbind(G);  fi;


#############################################################################
##
#F  ANUPQerrorPq( <param> ) . . . . . . . . . . . . . . . . . report an error
##
ANUPQerrorPq := function( param )
    Error(
    "Valid Options:\n",
    "    \"ClassBound\", <bound>\n",
    "    \"Prime\", <prime>\n",
    "    \"Exponent\", <exponent>\n",
    "    \"Metabelian\"\n",
    "    \"OutputLevel\", <level>\n",
    "    \"Verbose\"\n",
    "    \"SetupFile\", <file>\n",
    "Illegal Parameter: \"", param, "\"" );
end;


#############################################################################
##
#F  ANUPQextractPqArgs( <args> )  . . . . . . . . . . . . . extract arguments
##
ANUPQextractPqArgs := function( args )
    local   CR,  i,  act,  match;

    # allow to give only a prefix
    match := function( g, w )
    	return 1 < Length(g) and SubString(w,1,Length(g)) = g;
    end;

    # extract arguments
    CR := rec();
    i  := 2;
    while i <= Length(args)  do
        act := args[i];

    	# "ClassBound", <class>
        if match( act, "ClassBound" ) then
            i := i + 1;
            CR.ClassBound := args[i];

    	# "Prime", <prime>
        elif match( act, "Prime" )  then
            i := i + 1;
            CR.Prime := args[i];

    	# "Exponent", <exp>
        elif match( act, "Exponent" )  then
            i := i + 1;
            CR.Exponent := args[i];

        # "Metabelian"
        elif match( act, "Metabelian" ) then
            CR.Metabelian := true;

    	# "Output", <level>
        elif match( act, "OutputLevel" )  then
            i := i + 1;
            CR.OutputLevel := args[i];
    	    CR.Verbose     := true;

    	# "SetupFile", <file>
        elif match( act, "SetupFile" )  then
    	    i := i + 1;
            CR.SetupFile := args[i];

    	# "Verbose"
        elif match( act, "Verbose" ) then
            CR.Verbose := true;

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
#F  Pq( <G>, ... )  . . . . . . . . . . . . . . . . . . . . .  prime quotient
##
if IsBound(G)  then ThisIsAHack := G;  else ThisIsAHack := "G";  fi;
G := ThisIsAHack;

Pq := function( arg )
    local   CR,  F,  file,  out,  oldG,  newG,  x,  gens,  cmd,  r;

    # check arguments
    if Length(arg) < 1  then
    	Error( "usage: Pq( <F>, <control-args>, ... )" );
    fi;
    F := arg[1];
    if not IsFpGroup(F)  then
    	Error( "<F> must be a finitely presented group" );
    fi;
    if Length(arg) < 2 or not IsRec(arg[2])  then
    	CR := ANUPQextractPqArgs( arg );
    	if IsString(CR)  then
    	    ANUPQerrorPq(CR);
    	fi;
    else
    	CR := ShallowCopy(arg[2]);
    	x := Set( RecFields(CR) );
    	SubtractSet( x, Set( [ "Exponent", "ClassBound", "Verbose",
    	    	    	       "Prime", "OutputLevel", "SetupFile" ] ) );
    	if 0 < Length(x)  then
    	    ANUPQerrorPq(x);
    	fi;
    fi;

    # at least "Prime" and "Class" must be given
    if not IsBound(CR.Prime)  then
    	Error( "you must supply a prime" );
    fi;
    if not IsBound(CR.ClassBound)  then
    	Error( "you must supply a class bound" );
    fi;

    # set default values
    if not IsBound(CR.Exponent)  then 
    	CR.Exponent := 0;
    fi;
    if not IsBound(CR.SetupFile)  then 
    	file := TmpName();
    else
    	file := CR.SetupFile;
    fi;
    if not IsBound(CR.Verbose) then 
    	CR.Verbose := false;
    fi;

    PrintTo( file, "#input file for pq\n" );

    # setup input file
    AppendTo( file, "1\n" );
    AppendTo( file, "prime ", CR.Prime, " \n" );
    AppendTo( file, "class ", CR.ClassBound, " \n" );
    if CR.Exponent <> 0  then
        AppendTo( file, "exponent ", CR.Exponent, "\n" );
    fi;
    if IsBound(CR.Metabelian)  then 
        AppendTo( file, "metabelian\n" );
    fi;
    if IsBound(CR.OutputLevel)  then
    	AppendTo( file, "output ", CR.OutputLevel, " \n" );
    fi;

    # create generic generators "g1" ... "gn"
    gens := WordList( Length(F.generators), "g" );
    AppendTo( file, "generators {" );
    for x  in gens  do
    	AppendTo( file, x, ", " );
    od;
    AppendTo( file, " }\n" );

    # write the presentation using these generators
    AppendTo( file, "relations {" );
    if IsBound(F.relators)  then
    	for r  in F.relators  do
    	    AppendTo( file, MappedWord(r,F.generators,gens), ",\n" );
    	od;
    fi;
    AppendTo( file, "}\n;\n" );
    
    # if we only want to setup the file we are ready now
    if IsBound(CR.SetupFile)  then
        AppendTo( file, "8\n25\nPQ_OUTPUT\n2\n0\n0\n" );
    	Print( "#I  input file '", CR.SetupFile, "' written,\n",
    	       "#I    run 'pq' with '-k' flag, the result will be saved in ",
    	       "'PQ_OUTPUT'\n" );
    	return true;
    fi;

    # otherwise append code to save the output in a temporary file
    out := TmpName();
    AppendTo( file, ConcatenationString( "8\n25\n", out, "\n2\n0\n" ) );
    AppendTo( file, "0\n" );

    # and finally start the pq
    cmd:=ConcatenationString("-v ",EXENAME,"gap.sh ");
    if CR.Verbose  then 
        #cmd := ConcatenationString( "pq -k < ", file );
        cmd := ConcatenationString( cmd,"-k < ", file );
    else 
        #cmd := ConcatenationString( "pq -k < ", file, " > /dev/null" );
        cmd := ConcatenationString( cmd, "-k < ", file, " > /dev/null" );
    fi;
    ExecPkg( "anupq", "bin/pq", cmd, "." );
    #Exec(cmd);
    
    # read group from file
    if IsBound(G)  then oldG := G;  fi;
    if not READ(out)  then
        Exec( ConcatenationString( "rm -f ", file, " ", out ) );
    	Error( "cannot execute ANU pq,  please check installation" );
    fi;
    G.pqFpGroup := F;
    newG := G;
    Unbind(G);
    if IsBound(oldG)  then G := oldG;  fi;

    # remove intermediate files and return
    Exec( ConcatenationString( "rm -f ", file, " ", out ) );
    return newG;

end;


#############################################################################
##
#F  PqHomomorphism( <G>, <img> )  . . . . transfer a homomorphism to quotient
##
PqHomomorphism := function( G, img )
    local   gi;

    # map the images into the quotient
    gi := List( img, x -> MappedWord(x,G.pqFpGroup.generators,G.pqImages) );

    # construct a homomorphism and return
    return GroupHomomorphismByImages( G, G, G.pqImages, gi );

end;
