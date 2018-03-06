##############################################################################
##
#A  sq.g                        29 Mar 1994                   Alice C Niemeyer
#A                                   and functions of            Werner Nickel
##
##  This file contains the interface to the ANU SQ program. 
##    

#############################################################################
##
#V  The level of information printed by the ANU Sq is determined by
##  InfoSq*. It has the following effect
##
##
##  InfoSq1 :   the program will indicate which prime is currently used
##              and the dimension of the computed module
##  InfoSq2 :   the program will also mark the begin of the execution
##                 of the following basic steps in the algorithm
##                    - Call to AddDefinitions() which adds new generators
##                    - Call to Consistency(), the consistency check function
##                    - Call to LiftEpimorphism(), to lift the epimorphism
##                    - Call to the vector enumerator
##                    - Call to UpdatePresentation, the function which updates
##                            the presentation according to VE output
##               and the presentation at the end of a completion 
##                 of pStep() is printed
##
InfoSq1 := Print;
InfoSq2 := Ignore;


#############################################################################
##
#V  SqPresentation  . . . . . . . . . . . . . . . . . . . .  globale variable
##
SqPresentation := [];


#############################################################################
##
#F  SqUsage() . . . . . . . . . . . . . . . show usage of 'NilpotentQuotient'
##
SqUsage := function()
    return Error("usage: Sq( <file>|<fpgroup> [<Lseries>] )");
end;
	

#############################################################################
##
#F  SqPresPrintToFile( <file>, <fp> ) . . . . . .  print presentation to file
##
##  Print a finite presentation in Sq format. 
##  (modified from Werner Nickel's function NqPresPrintToFile)
##
SqPresPrintToFile := function( file, fp )
    local   i,  gens,  append,  size;

    # append a relator (this is a hack)
    append := function( rel )
	local   pos,  len,  max;

	max := 10;
	pos := 1;
	rel := MappedWord( rel, fp.generators, gens );
	len := LengthWord(rel);
	while 0 < len  do
	   if len <= max  then
	        AppendTo( file, Subword(rel,pos,pos+len-1) );
		pos := pos+len;
		len := 0;
	    else
	    	AppendTo( file, Subword(rel,pos,pos+max-1), "*\n      " );
		pos := pos+max;
		len := len-max;
	    fi;
	od;
    end;

    # raise screen size
    size := SizeScreen();
    SizeScreen( [ 100, 24 ] );

    # print presentation to file using generators "x1" ... "xn"
    PrintTo( file, "< " );
    if 0 < Length(fp.generators)  then
        gens := WordList( Length(fp.generators), "x" );
	AppendTo( file, gens[1] );
    fi;
    for i  in [2..Length(fp.generators)]  do
	AppendTo( file, ", ", gens[i] );
    od;
    AppendTo( file, " |\n    " );
    if IsBound(fp.relators)  then
        if 0 < Length(fp.relators)  then
	    append( fp.relators[1] );
	fi;
	for i  in [2..Length(fp.relators)]  do
	    AppendTo( file, ",\n    " );
	    append( fp.relators[i] );
	od;
    fi;
    AppendTo( file, "\n>\n" );

    # restore screen size
    SizeScreen( size );

end;

IsLseries := function( ls )
	local i, l;

	l := Length(ls);
        for i in [ 1 .. l ] do
	    if not IsPrime(ls[i][1]) then return false; fi;
        od;
	return true;
end;

#############################################################################
##
#F  Sq( <F>, <Lseries> )  . . . . . . . . . . . . . . soluble quotient of <F>
##
##  The interface to the Solublequotient standalone.
##
Sq := function( arg )
    local   ll,  lseries,  dir,  name,  res,  cmd;

    if not Length(arg) in [1,2]  then SqUsage();  fi;

    # Check if a list is an L-series

    lseries := [];
    if Length(arg) = 2 then
	lseries := arg[2];
	if not IsLseries(lseries)  then return SqUsage();  fi;
    fi;

    # create a tmp directory
    dir := TmpName();
    Exec(Concatenation( "mkdir ", dir ));
    name := Concatenation( dir, "/SQ_INPUT" );
    res  := Concatenation( dir, "/SQ_OUTPUT" );

    # set up SQ input file SQ_INPUT in <dir>
    if IsRec( arg[1] ) then
	SqPresPrintToFile( name, arg[1] );
    elif IsString( arg[1] ) then
        Exec(Concatenation( "cp ", arg[1], " ", name ));
    else
        Exec(Concatenation( "rm -rf ", dir ));
	return SqUsage();
    fi;

    cmd:="";
#   if InfoSq1 = Print then
#       cmd := ConcatenationString( cmd, " -p 1" );
#   elif InfoSq2 = Print then
#       cmd := ConcatenationString( cmd, " -p 2" );
#   else
#       cmd := ConcatenationString( cmd, " -p 0" );
#   fi;

    # add lower central series if known
    if lseries <> [] then
	for ll in lseries do 
	    AppendTo( name, ll[1], " ", ll[2], "\n" );
	od;
    fi;

    # call the sq
    cmd := ConcatenationString( cmd, " -v", LOADED_PACKAGES.ve,"bin/");
    cmd := ConcatenationString( cmd, " < ", name, " " );
    cmd := ConcatenationString( cmd, " > ", res );
    ExecPkg( "anusq", "bin/Sq", cmd, dir );

    # read in the result
    AppendTo( res, ";\n" );
    Read(res);

    # remove the tmp dir
    Exec(Concatenation( "rm -rf ", dir ));

    # and return
    return SqPresentation;

end;
