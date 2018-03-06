PrintPkgInit(rec(name:="ve",version:="3.04"));
#############################################################################
##
##  Notify the functions to be defined for {\VE}.
##
AUTO( ReadPkg( "ve", "gap", "ve" ),
      VE, PrintVEInput, CallVE, VEOutput );

AUTO( ReadPkg( "ve", "gap", "mtx" ),
      VEMeatAxe, VEMeatAxeMat );

# the following forces to read lib/algfp.g but this is unavoidable
#############################################################################
##
#F  FpAlgebraOps.OperationQuotientModule( <A>, <Q>, <opr> )
#F  FpAlgebraOps.OperationQuotientModule( <A>, <Q>, "mtx" )
##
##  takes a finitely presented algebra <A> and a quotient module <Q> of a
##  free module, and returns the matrix representation computed by {\VE}.
##
##  If the third argument is the string '\"mtx\"', the output is an algebra
##  of {\MeatAxe} matrices.
##
FpAlgebraOps.OperationQuotientModule := function( A, M, opr )

  local file,       # stem of the presentation file name
	i,          #
	alpha,      # alphabet over which names of generators are written
	lalpha,     # length of 'alpha'
	names,      # list of new names for algebra generators
	filepres,   # name of the presentation file
	outfile,    # stem of the output file name
	commandstr, # command string for {\VE}
	output,     # output of {\VE}
	result;     # returned matrix algebra

  # Check the arguments.
  if not IsFpAlgebra( A ) or not IsList( M ) then
    Error( "<A> must be f.p. algebra, <M> list of submodule generators" );
  fi;

  # Check that the input can be processed by 'me', 'qme', or 'zme'.
  if A.field.char = 0 then

    if   IsIdentical( A.field, Integers ) then
#T no '=' for Integers ...
      commandstr:= "zme";
    elif IsIdentical( A.field, Rationals ) then
      commandstr:= "qme";
    else
      Error( "characteristic 0: 'Integers' and 'Rationals' only" );
    fi;
    
    if opr = "mtx" then
      Error( "MeatAxe output only for nonzero characteristic" );
    fi;

  elif A.field.char > 255 then
    Error( "'me' allows finite characteristic up to 255 only" );
  else
    commandstr:= "me";
  fi;

  # Construct a list 'names' of names for the generators
  # consisting only of letters.
  # Provide that no generator name is initial part of another
  # by choosing nonzero minimal length, namely
  # if we have $k$ generators and an alphabet of length $n$ then
  # choose words of length $i+1$ where $n^i \< k \leq n^{i+1}$,
  # the first word having number $n + n^2 + \cdots + n^i$.
  alpha:= "abcdefghijklmnopqrstuvwxyz";
  lalpha:= Length( alpha );
  i:= 1;
  while lalpha ^ i < Length( A.generators ) do i:= i+1; od;
  i:= Sum( [ 1 .. i-1 ], x -> lalpha^x );
  names:= List( [ 1 .. Length( A.generators ) ],
		x -> WordAlp( alpha, x+i ) );
  
  # Produce the input file for {\VE}.
  file:= TmpName();
  filepres:= Concatenation( file, ".pres" );
  PrintTo( filepres, PrintVEInput( A, M, names ) );

  # Prepare the output file name.
  outfile:= TmpName();

  # Call 'me', 'qme', or 'zme' with standard options.
  if opr <> "mtx" then

    # Choose {\GAP} output.
    CallVE( commandstr, file, outfile, Concatenation( VE.options, " -G" ) );

    # Get the output.
    output:= VEOutput( A, M, names, outfile );

  else

    # Choose {\MeatAxe} output.
    CallVE(commandstr, file, outfile, Concatenation( VE.options, " -m -H" ) );

    # Get the output.
    output:= VEOutput( A, M, names, outfile, "mtx" );

  fi;

  # Make the algebra.
  # Check for the special case of total collapse.
  if output.gens[1] = [] or Dimensions( output.gens[1] )[1] = 0 then

    output.operation.genimages:= List( output.gens, x -> EmptyMat );
    result:= NullAlgebra( A.field );

  else

    result:= UnitalAlgebra( A.field, output.gens );

  fi;

  result.operation:= output.operation;
  result.operation.genpreimages:=
	A.generators{ List( result.generators,
		      x -> Position( output.operation.genimages, x ) ) };

  # Make clean.
  EXEC( "rm ", filepres );

  # Return the algebra.
  return result;
end;
