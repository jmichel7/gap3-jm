#############################################################################
##
#A  sisgprin.g               GAP Share Library               Martin Wursthorn
##
#A  @(#)$Id: sisgprin.g,v 3.0 1994/05/19 14:09:41 sam Exp $
##
#Y  Copyright 1994-1995,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains those functions of the interface between {\SISYPHOS}
##  and {\GAP} that deal with group rings of $p$-groups.
##
#H  $Log: sisgprin.g,v $
#H  Revision 3.0  1994/05/19  14:09:41  sam
#H  Initial Revision under RCS
#H
##

#############################################################################
##
#F  NormalizedUnitsGroupRing( <P> )
#F  NormalizedUnitsGroupRing( <P>, <n> )
##
##  returns the group of normalized units of the group ring $FP$ of the
##  $p$-group <P> over the field $F$ with $p$ elements.
##
##  If a second argument <n> is given, the group of normalized units of
##  $FP / I^n$ is returned, where $I$ denotes the augmentation ideal of
##  $FP$.
##
##  The returned group is represented as polycyclicly presented group.
##
NormalizedUnitsGroupRing := function( arg )

    local P,       # $p$-group, argument
          n,       # power of the augmentation ideal, argument
          f1, f2,  # files for input and output
          globp,   # save global variable 'p'
          isoP,    # isomorphism from 'P' onto group that is compatible
                   # with the Jennings series of 'P'
          jser,    # Jennings series of 'p'
          weights, # list of Jennings weights of the generators of 'p'
          mtable,  # string containing SISYPHOS command for computing
                   # a multiplication table for 'p' (or empty)
          tsize,   # size of the multiplication table
          j,       # loop variable
          i;       # loop variable

    # Check the arguments.
    if Length( arg ) < 1 or Length( arg ) > 2
       or not IsGroup( arg[1] )
       or ( Length( arg ) = 2 and not IsInt( arg[2] ) ) then
      Error( "usage: NormalizedUnitsGroupRing( <P> [, <n>] )" );
    fi;

    P:= arg[1];

    f1:= SISYPHOS.TmpName1(); PrintTo( f1, " " );
    f2:= SISYPHOS.TmpName2();

    if IsBound( p ) then
      globp:= p;
    fi;

    # Compute an isomorphic group with presentation compatible
    # with the Jennings series.
    isoP:= IsomorphismAgGroup( JenningsSeries( P ) );
    p:= isoP.range;
    if not IsBound( p.1 ) then
      for i in [ 1 .. Length( p.generators ) ] do
        p.(i):= p.generators[i];
      od;
    fi;

    # Compute the weights of the generators w.r. to the Jennings series.
    jser:= JenningsSeries( p );
    weights:= [];
    for i in p.generators do
      j:= 2;
      while i in jser[j] do
        j:= j+1;
      od;
      Add( weights, j-1 );
    od;

    if Length( arg ) = 2 then
      n:= Maximum ( arg[2], 1 );
    else
      n:= Length( DimensionsLoewyFactors( p ) );
    fi;

    # compute amount of memory needed
    SISYPHOS.SISPMEM := "300000";
    tsize := Size ( p ) * Size ( p ) * 4;
    if tsize < 5400000 then
        tsize := tsize + 600000;
        mtable := "use (multiplication table);\n";
        Print ("#D use multiplication table\n");    
    else
        tsize := 600000;
        mtable := "";
    fi;
    SISYPHOS.SISTMEM := String ( tsize );

    # Prepare the input file for {\SISYPHOS}
    PrintTo( f1, PrintSisyphosInputPGroup( p, "p", "pcgroup", weights ),
                 "setdomain( groupring( p ) );\n",
                 mtable,
                 "unitgroup( ", n, ",\"SISYPHOS.SISISO\", 1 );\n",
                 "quit;\n" );

    SISYPHOS.SISISO:= 0;

    # Call {\SISYPHOS}, read the output, make clean.
    ExecPkg( "sisyphos", SISYPHOS.SISCALL,
             Concatenation( " -t ", SISYPHOS.SISTMEM,
                            " -m ", SISYPHOS.SISPMEM,
                            " <", f1, " >", f2 ), "." );
    Read( f2 );
    SISYPHOS.Exec( Concatenation( "rm ", f1, " ", f2 ) );

    # Check whether the output file contained the result.
    if SISYPHOS.SISISO = 0 then
      Error( "output file was not readable" );
    fi;

    # Restore global variable 'p'.
    if IsBound( globp ) then p:= globp; fi;

    P:= SISYPHOS.SISISO;
    Unbind( SISYPHOS.SISISO );

    # Return the result.
    return P;
    end;

#############################################################################
##
#E  Emacs . . . . . . . . . . . . . . . . . . . . . . . local emacs variables
##
##  Local Variables:
##  mode:               outline
##  outline-regexp:     "#F\\|#V\\|#E"
##  fill-column:        73
##  fill-prefix:        "##  "
##  eval:               (hide-body)
##  End:
##
