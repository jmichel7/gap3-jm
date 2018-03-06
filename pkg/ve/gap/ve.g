#############################################################################
##
#A  ve.g                    GAP Share Library                    Steve Linton
##
#A  @(#)$Id: ve.g,v 1.1.1.1 1996/12/11 12:39:41 werner Exp $
##
#Y  Copyright 1994-1995,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains the interface between {\GAP} and {\VE}.
##
#H  $Log: ve.g,v $
#H  Revision 1.1.1.1  1996/12/11 12:39:41  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 3.0  1994/03/24  17:14:35  sam
#H  Initial Revision under RCS
#H
##

#############################################################################
##
#V  VE . . . . global variable used for the interface to {\VE}
##
##  The {\VE} share library uses the global variable 'VE' for
##  storing necessary information to call the standalones; intermediate
##  results that are read from files are stored in the component 'out'.
##
##  'VE' is a record with components
##
##  'Path': \\ the full path name of the directory that contains the
##             executables of 'me', 'qme', 'zme',
##
##  'options': \\ a string with the options for {\VE}, will be appended to
##                the command string, so may overwrite chosen options,
##                thus you should *not* use '-o' in this variable.
##                The options should be among 'a', 'e', 'l', 'v'.
##  the '-a' and '-e' options of {\VE}; 'false' means the
##                         default value, otherwise the values must be
##                         strings (see Chapter 4 of the {\VE} manual for
##                         details).
#T options!!
##
VE := rec( options := "",
           Path := Concatenation( PKGNAME[1], "ve/bin/" ) );

#############################################################################
##
#F  PrintVEInput( <A>, <M>, <names> )  . . .  input for the Vector Enumerator
##
##  takes a finitely presented algebra <A> and submodule generators <M>
##  (a list of vectors over <A>), and a list <names> of names the generators
##  have in the presentation for {\VE},
##  and prints the presentation to be input into the vector enumeration
##  program.
##
PrintVEInput := function( A, M, names )

    local i, j, k,            # loop variables
          char,               # characteristic of the algebra
          ncols,              # line length
          free,               # space in actual line
          Pr,                 # print function taking care of line length
          printalgelmstring,  # local function printing algebra elements
          rels,               # relators of the algebra
          grouptyperels,      # group type relators
          otherrels,          # non-group type relators
          id,
          one,
          invertible,
          involutions,
          dim,
          entry,
          Agenerators;
    
    ncols:= SizeScreen()[1] - 2;
    free:= ncols;

    # 1. Characteristic of the algebra
    char:= A.field.char;
    if char <> 0 then
      if char > 255 or A.field <> GF( char ) then
        Error( "calculations only over prime fields of char < 256" );
      fi;
    fi;
    Print( char, ".\n" );
    
    Agenerators:= List( A.generators, x -> x.words[1] );

    # print <string> without producing '\' at the end of the lines
    Pr:= function( string )
    local len;
    len:= Length( string );
    if len <= free then
      Print( string );
      free:= free - len;
    else
      Print( "\n", string );
      free:= ncols - len;
    fi;
    end;

    # local function printing an algebra element,
    # in terms of the new generator names
    printalgelmstring:= function( elm )
    local j, k, coeff;
    
    # 'elm' is either an algebra element or a word
    if not IsRec( elm ) then
      elm:= rec( words:= [ elm ], coeff:= [ A.field.one ] );
    fi;
    
    # zero element
    if Length( elm.words ) = 0 then
      Pr( "0" );
    fi;
    
    for j in [ 1 .. Length( elm.words ) ] do
      if char = 0 then
        coeff:= elm.coeff[j];
      else
        coeff:= Int( elm.coeff[j] );
      fi;
      if j > 1 and coeff > 0 then
        Pr( "+" );
      fi;
      if coeff <> 0 and coeff <> 1 then
        Pr( String( coeff ) );
        if elm.words[j] <> id then
          Pr( "*" );
        fi;
      fi;
      if coeff = 1 and elm.words[j] = id then
        Pr( "1" );
      fi;
      for k in [ 1 .. LengthWord( elm.words[j] ) ] do
        if k > 1 then
          Pr( "*" );
        fi;
        Pr( names[ Position( Agenerators,
                             Subword( elm.words[j], k, k ) ) ] );
      od;
    od;
    end;
    
    # Distribute algebra relators to
    # group type relators and other relators.
    # 'grouptyperels' will contain all <w> for relators of form
    # '\pm(<w>-1)' in '<A>.relators'.
    # 'otherrels' will contain the other relators *themselves*,
    # i.e. as records.
    # *Note* that at this stage also those relators are considered
    # to be group type relators that may contain non-invertible
    # generators.
    rels:= A.relators;
    grouptyperels:= [];
    otherrels:= [];
    id:= A.generators[1].words[1] ^ 0;
    one:= A.field.one;
    for entry in rels do
      if     Length( entry.words ) = 2
         and entry.words[1] = id
         and ( ( entry.coeff[1] = one and entry.coeff[2] = - one )
            or ( entry.coeff[2] = one and entry.coeff[1] = - one ) )
        then
        Add( grouptyperels, entry.words[2] );
      else
        Add( otherrels, entry );
      fi;
    od;
    
    # 2. Generators of the algebra:
    #    new names, separated by ' '
    for i in names do
      Pr( i );
      Pr( " " );
    od;
    Pr( ".\n" );
    free:= ncols;
    
    # 2a. Generators not known to be invertible:
    #     Take the difference of all generators and those
    #     that have a left and right inverse.
    invertible:= Intersection( 
                  List( grouptyperels,
                        x -> Subword( x, 1, 1 ) ),
                  List( grouptyperels,
                        x -> Subword( x, LengthWord(x), LengthWord(x) ) ) );
    
    if Length( invertible ) = 0 then
    
      # no invertible generators
      Pr( "*..\n" );
      free:= ncols;
    
    else
    
      for i in Difference( Agenerators, invertible ) do
        Pr( names[ Position( Agenerators, i ) ] );
        Pr( " " );
      od;
      Pr( ".\n" );
      free:= ncols;
    
    # 2b. Generators not known to be involutions:
    #     Take those invertible generators for that
    #     the square is not a group type relator.
    
      involutions:= Filtered( invertible, x -> x^2 in grouptyperels );
      if Length( involutions ) = 0 then
        Pr( "*.\n" );
        free:= ncols;
      else
        for i in Difference( invertible, involutions ) do
          Pr( names[ Position( Agenerators, i ) ] );
          Pr( " " );
        od;
        Pr( ".\n" );
      fi;
      free:= ncols;
    fi;
    
    # Move those relators from 'grouptyperels' to 'otherrels'
    # that contain non-invertible generators.
    for i in [ 1 .. Length( grouptyperels ) ] do
      if ForAny( [ 1 .. LengthWord( grouptyperels[i] ) ],
                 x -> not Subword( grouptyperels[i], x, x ) in invertible ) then
        AddSet( otherrels,
                FpAlgebraElement( A, A.field.one * [ 1, -1 ],
                                          [ grouptyperels[i], id ] ) );
        Unbind( grouptyperels[i] );
      fi;
    od;
    
    # 4. Submodule Generators
    if IsRec( M ) then M:= M.generators; fi;
    
    # rank of the module
    dim:= Length( M[1] );
    Pr( "{" );
    Pr( String( dim ) );
    Pr( "}" );
    
    # generators of the submodule, separated by ',',
    # each a enclosed in '( )',
    for i in [ 1 .. Length( M ) ] do
      Pr( "(" );
      for j in [ 1 .. dim ] do
        printalgelmstring( M[i][j] );
        if j < dim then
          Pr( "," );
        fi;
      od;
      Pr( ")" );
      if i < Length( M ) then
        Pr( ",\n" );
        free:= ncols;
      fi;
    od;
    Pr( ".\n" );
    free:= ncols;
    
    # 3. Algebra Relators, first the group type relators ...
    for i in [ 1 .. Length( grouptyperels ) ] do
      if IsBound( grouptyperels[i] ) then
        printalgelmstring( grouptyperels[i] );
        Pr( ", " );
      fi;
    od;
    Pr( ":\n" );
    free:= ncols;
    
    #    ... and then the others
    for i in otherrels do
      printalgelmstring( i );
      Pr( " = 0, " );
    od;
    
    # The presentation is complete.
    Print( ".\n" );
    end;

#############################################################################
##
#F  CallVE( <commandstr>, <infile>, <outfile>, <options> )
##
CallVE := function( commandstr, infile, outfile, options )

#T ExecPkg?
    EXEC( VE.Path,
          commandstr,
          " -p ", infile,          # presentation file
          " -o ", outfile,         # output file
          " -i",                   # images
          " -P",                   # preimages
          " -v0",                  # suppress comments
          " -L'#I  '",             # but set prefix if someone wants comments
          " -Y VE.out ",           # assign result to 'VE.out'
          options,                 # and add the other options
          "" );
    end;

#############################################################################
##
#F  VEOutput( <A>, <M>, <names>, <outfile> [, "mtx"] )
##
VEOutput := function( arg )

    local A,          # 1st argument
          M,          # 2nd argument
          names,      # 3rd argument
          outfile,    # 4th argument
          result,     # matrix algebra, result
          preim,      # preimages of module generators
          zero,       # zero of the algebra
          gens,       # list of generators
          dir,        # temporary directory name (for {\MeatAxe} output)
          i,          # loop over generators list
          name,       # one generator name (for {\MeatAxe} output)
          dimensions; # dimensions of matrices (for {\MeatAxe} output)

    A:=       arg[1];
    M:=       arg[2];
    names:=   arg[3];
    outfile:= arg[4];

    # Read the output, and make clean.
    VE.out:= false;
    Read( Concatenation( outfile, ".g" ) );
    EXEC( "rm ", outfile, ".g " );

    # Check whether the output file contained the result.
    if VE.out = false then
      Error( "output file was not readable" );
    fi;
    result:= VE.out;
    Unbind( VE.out );

    if Length( arg ) = 4 then

      # {\GAP} output.
      # 1. Decode the 'preImages' component.
      preim:= [];
      zero:= List( [ 1 .. Length( M[1] ) ], x -> A.zero );
      gens:= List( names, i -> result.( Concatenation( "g", i ) ) );
      for i in [ 1 .. Length( result.preImages ) ] do
        preim[i]:= ShallowCopy( zero );
        preim[i][ result.preImages[i].modGen ]:=
          MappedWord( result.preImages[i].word, gens, A.generators );
      od;

      # 2. Add the information needed to construct
      #    the operation (algebra) homomorphism.
      result.operation:= rec( genimages  := result.gens,
                              moduleinfo := rec(
                                  preimagesBasis := preim,
                                  genimages := result.images_mat ) );

    else

      # {\MeatAxe} output (constructed with option '-H')
      # 1. Get the dimensions of the matrices.
      dimensions:= [ result.dim, result.dim ];
      
#T collapse?
      # 2. Convert to internal format.
      #    Change the file names in order to have a common stem, extended
      #    by numbers only.
      dir:= VEMeatAxe.TmpName();
      EXEC( "mkdir ", dir );
      AddSet( VEMeatAxe.dirs, dir );

      gens:= [];
      for i in [ 1 .. Length( names ) ] do

        # The new file name
        name:= Concatenation( dir, "/g.", String( i ) );

#T field info file?
        # Convert to internal format.
        EXEC( VEMeatAxe.PATH, "zcv ", outfile, ".", names[i], " ", name );

        # Make clean.
        EXEC( "rm ", outfile, ".", names[i] );

        # Notify the {\MeatAxe} matrix.
        gens[i]:= VEMeatAxeMat( name, A.field, dimensions );
        gens[i].gensname:= Concatenation( dir, "/g" );

      od;

      # 3. Store the generators.
      result.gens:= gens;

      # 4. Add the information needed to construct
      # the operation (algebra) homomorphism.
      result.operation:= rec( genimages := gens,
                              gensname  := Concatenation( dir, "/g" ) );

    fi;

    return result;
    end;
