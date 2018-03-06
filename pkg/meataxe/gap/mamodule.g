if not IsBound(MetaAxe) then MeatAxe:=MeatAxe;fi;
#############################################################################
##
#A  mamodule.g                  GAP Library                     Thomas Breuer
##
#A  @(#)$Id: mamodule.g,v 1.1.1.1 1996/12/11 12:40:13 werner Exp $
##
#Y  Copyright 1994-1995,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains the functions for {\Meataxe} modules.
##
#H  $Log: mamodule.g,v $
#H  Revision 1.1.1.1  1996/12/11 12:40:13  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 3.2  1996/11/13 15:43:18  sam
#H  replaced 'mkpeak' and 'gkond' by 'pwkond' twice
#H
#H  Revision 3.1  1994/05/19 14:05:19  sam
#H  changes for GAP-3.4
#H
#H  Revision 3.0  1994/03/28  15:31:13  sam
#H  initial revision under RCS
##

#############################################################################
##
#F  InfoMeatAxe( ... ) . . . . . . . info function for the {\MeatAxe} package
##
    if not IsBound( InfoMeatAxe ) then InfoMeatAxe:= Ignore; fi;

#############################################################################
##
#V  MeatAxeModuleOps . . . . . . . . operations record for {\MeatAxe} modules
##
MeatAxeModuleOps := OperationsRecord( "MeatAxeModuleOps", ModuleOps );

#############################################################################
##
#F  MeatAxeModuleOps.Sorry( arg )
##
MeatAxeModuleOps.Sorry := function( arg )
    Error( "this function is not provided for MeatAxe modules" );
    end;

MeatAxeModuleOps.Elements := MeatAxeModuleOps.Sorry;
MeatAxeModuleOps.Random :=   MeatAxeModuleOps.Sorry;
MeatAxeModuleOps.Zero :=     MeatAxeModuleOps.Sorry;

#############################################################################
##
#F  MeatAxeModuleOps.GapObject( <M> )
##
MeatAxeModuleOps.GapObject := function( M )
    if   IsBound( M.factorNum ) then
      return GapObject( M.factorNum ) / GapObject( M.factorDen );
    elif IsBound( M.parent ) then
      return Submodule( GapObject( Parent( M ) ),
                        GapObject( Generators( M ) ) );
    else
      return Module( GapObject( M.ring ), GapObject( Generators( M ) ) );
    fi;
    end;

#############################################################################
##
#F  MeatAxeModuleOps.Submodule( <M>, <gens> )
##
MeatAxeModuleOps.Submodule := function( M, gens )
    local S;
    S:= ModuleOps.Submodule( M, gens );
    S.operations:= MeatAxeModuleOps;
    return S;
    end;

#############################################################################
##
#F  MeatAxeModuleOps.\=( M1, M2 )
##
MeatAxeModuleOps.\= := function( M1, M2 )
    if M1.ring = M2.ring and M1.generators = M2.generators then
      return true;
    else
      return false;
    fi;
    end;

#############################################################################
##
#F  MeatAxeModuleOps.\<( M1, M2 )
##
MeatAxeModuleOps.\< := function( M1, M2 )
    return    M1.ring < M2.ring
           or ( M1.ring = M2.ring and M1.generators < M2.generators );
    end;

#############################################################################
##
#F  MeatAxeModuleOps.Intersection( <M1>, <M2> ) . . . intersection of modules
##
MeatAxeModuleOps.Intersection := function( M1, M2 )
    local intersect;           # intersection of modules, result
    if Dimensions( M1.generators )[2] = Dimensions( M2.generators )[2] 
       and M1.ring = M2.ring
       and M1.generators.field = M2.generators.field  then
      intersect:= MeatAxeMatOps.SumIntersectionSpaces( M1.generators,
                                                       M2.generators );
      MeatAxe.Unbind( intersect[1] );
      intersect:= Module( M1.ring, intersect[2] );
      intersect.basis:= Basis( intersect, intersect.generators );
    else
      Error( "need same field, same number of columns, same acting ring" );
    fi;
    return intersect;
    end; 

#############################################################################
##
#F  MeatAxeModuleOps.CompositionFactors( <M> )
##
##  returns a list of {\MeatAxe} modules which are the factors of a
##  composition series of the {\MeatAxe} module <M>.
##
##  Each factor occurs in this list with its multiplicity in <M>.
##
MeatAxeModuleOps.CompositionFactors := function( M )

    local G,        # the algebra acting on 'M'
          f,        # finite field 'G' is written over
          name,     # the initial part of names of generators of the group
          name2,    # temporary file name
          cf,       # list of composition factors, result
          i,        # loops over composition factors
          j,        # loops over generators of composition factor
          gens,     # generators of composition factor
          dims;     # dimensions of composition factors
 

    G:= Parent( M ).ring;
    f:= G.field;

    # zero module
    if Dimension( M ) = 0 then return []; fi;

    # first get the equivalent natural module if necessary
    if Dimension( M ) <> Dimensions( G.generators[1] )[1] then
      M:= NaturalModule( Operation( G, M ) );
      G:= M.ring;
    fi;

    name := G.generators[1].gensname;
    name2:= TmpName();

    MeatAxe.EXEC( "chop -G -g ", Length( G.generators ), " ", name,
          " > ", name2 );

    Unbind( MeatAxe.CompositionFactors );
    Unbind( MeatAxe.CompositionSeries );
    Read( name2 );
    EXEC( "rm ", name2 );

    # get the list of file names of actions on composition factors
    dims:= List( MeatAxe.CompositionFactors, x -> MeatAxe.NUMBER( x[1] ) );
    cf:= List( MeatAxe.CompositionFactors,
               x -> Concatenation( name, x[1] ) );

    if InfoMeatAxe = Print then
      Print( "#I   Name Mult  SF\n" );
      for i in [ 1 .. Length( MeatAxe.CompositionFactors ) ] do
        Print( "#I ", String( MeatAxe.CompositionFactors[i][1], 6 ),
                      String( MeatAxe.CompositionFactors[i][2], 5 ),
                      String( MeatAxe.CompositionFactors[i][3], 4 ), "\n" );
      od;
      Print( "#I\n#I  Ascending composition series:\n#I  " );
      for i in MeatAxe.CompositionSeries do
        Print( MeatAxe.CompositionFactors[i][1], "  " );
      od;
      Print( "\n" );
    fi;

    # notify the {\MeatAxe} matrices to {\GAP}
    for i in [ 1 .. Length( cf ) ] do
      gens:= [];
      for j in [ 1 .. Length( G.generators ) ] do
        gens[j]:= MeatAxeMat( Concatenation( cf[i], ".", String( j ) ),
                              f, [ dims[i], dims[i] ] );
        gens[j].gensname:= cf[i];
        if IsBound( G.generators[j].abstract ) then
          gens[j].abstract:= G.generators[j].abstract;
        fi;
      od;
      cf[i]:= Algebra( f, gens );
      if IsBound( G.freeAlgebra ) then
        cf[i].freeAlgebra:= G.freeAlgebra;
      fi;
    od;

    cf:= List( cf, NaturalModule );
    for i in [ 1 .. Length( cf ) ] do
      cf[i].isIrreducible:= true;
      cf[i].isAbsolutelyIrreducible:= MeatAxe.CompositionFactors[i][3] = 1;
      cf[i].sf:= MeatAxe.CompositionFactors[i][3];
      cf[i].compositionFactors:= [ cf[i] ];
    od;

    cf:= cf{ MeatAxe.CompositionSeries };

    Unbind( MeatAxe.CompositionFactors );
    Unbind( MeatAxe.CompositionSeries );

    return cf;
    end;

#############################################################################
##
#F  MeatAxeModuleOps.Constituents( <M> )
##
MeatAxeModuleOps.Constituents := MeatAxeModuleOps.CompositionFactors;

#############################################################################
##
#F  MeatAxeModuleOps.KroneckerProduct( <M1>, <M2> )
##
MeatAxeModuleOps.KroneckerProduct := function( M1, M2 )

    local gens,      # list of generators
          j,         # loop over 'gens'
          kron,      # the algebra of the kronecker product
          result;    # the Kronecker product of modules, result
    
    if Length( M1.ring.generators ) <> Length( M2.ring.generators )  then
      Error( "representations not compatible" );
    elif M1.ring.field <> M2.ring.field then
      Error( "representations are over different fields" );
    fi;
    
    gens:= List( [ 1 .. Length( M1.ring.generators ) ],
                 i -> KroneckerProduct( M1.ring.generators[i],
                                        M2.ring.generators[i] ) );

    # Attach abstract generators.
    for j in [ 1 .. Length( gens ) ] do
      if IsBound( M1.ring.generators[j].abstract ) then
        gens[j].abstract:= M1.ring.generators[j].abstract;
      fi;
    od;

    # Construct the algebra.
    kron:= Algebra( M1.ring.field, gens );
    
#T avoid to make this twice!
#T (Use 'mv'!)
    MeatAxe.Unbind( gens );
    
    if Dimension( M1 ) = Dimensions( M1.ring.generators[1] )[1] and
       Dimension( M2 ) = Dimensions( M2.ring.generators[1] )[1] then
    
      result:= NaturalModule( kron );
    
    else
      result:= Module( kron,
                       KroneckerProduct( M1.generators, M2.generators ) );
    fi;
    
    return result;
    end;

#############################################################################
##
#V  BasisMeatAxeModuleOps
##
BasisMeatAxeModuleOps := OperationsRecord( "BasisMeatAxeModuleOps" );
BasisMeatAxeModuleOps.Print := BasisRowSpaceOps.Print;
BasisMeatAxeModuleOps.Coefficients := MeatAxeModuleOps.Sorry;
BasisMeatAxeModuleOps.SiftedVector := MeatAxeModuleOps.Sorry;
BasisMeatAxeModuleOps.IsSemiEchelonBasis := MeatAxeModuleOps.Sorry;

#############################################################################
##
#V  SemiEchelonBasisMeatAxeModuleOps
##
SemiEchelonBasisMeatAxeModuleOps := OperationsRecord(
                                "SemiEchelonBasisMeatAxeModuleOps",
                                BasisMeatAxeModuleOps );
SemiEchelonBasisMeatAxeModuleOps.Print := SemiEchelonBasisRowSpaceOps.Print;

#############################################################################
##
#V  StandardBasisMeatAxeModuleOps
##
StandardBasisMeatAxeModuleOps := OperationsRecord(
                                "StandardBasisMeatAxeModuleOps",
                                BasisMeatAxeModuleOps );
StandardBasisMeatAxeModuleOps.Print := StandardBasisModuleOps.Print;

#############################################################################
##
#F  MeatAxeModuleOps.CanonicalBasis( <M> )
##
MeatAxeModuleOps.CanonicalBasis := MeatAxeModuleOps.Sorry;

#############################################################################
##
#F  MeatAxeModuleOps.SemiEchelonBasis( <M> [,<vectors>] )
##
##  Since 'zsb' does not give access to the semi-echelonized basis, we cannot
##  use 'zsb' here.
##  Instead we call 'zsp' for spinning up the space and computing a
##  semi-echelonized basis.
##
MeatAxeModuleOps.SemiEchelonBasis := function( arg )

    local M,         # the module
          vectors,   # basis vectors
          gens,      # generators of the acting ring
          name,      # file for echelonized module generators
          name2,     # file for 'vectors'
          name3;     # intermediate file for output message

    # Check and get the arguments
    if Length( arg ) < 1 or Length( arg ) > 2 or not IsModule( arg[1] ) then
      Error( "usage: SemiEchelonBasis( <M> [,<vectors>] )" );
    fi;

    M:= arg[1];

    if Length( arg ) = 2 then
      vectors:= arg[2];
    else

      # Make sure that the generators are linearly independent.
      name:= TmpName();
      MeatAxe.EXEC( "zef -G -Q ", M.generators.file, " ", name,
                    " > /dev/null" );

      name2:= MeatAxe.TmpName( MeatAxe.files );
      name3:= TmpName();

      gens:= Parent( M ).ring.generators;

      if Length( gens ) = 2 then

        MeatAxe.EXEC( "zsp -G -d ", M.generators.dimensions[2], " -n -ss ",
                      gens[1].file, " ", gens[2].file,
                      " ", name, " ", name2, " > ", name3 );

      else

        # the generators have appropriate names ...
        MeatAxe.EXEC( "zsp -G -d ", M.generators.dimensions[2], " -n -ss -g ",
                      Length( gens ), " ",
                      gens[1].gensname, " ", name, " ",
                      name2, " > ", name3 );

      fi;

      Unbind( MeatAxe.SplitDimensions );
      Read( name3 );

      # 'name2' contains the basis
      vectors:= MeatAxeMat( name2, M.generators.field,
                            [ MeatAxe.SplitDimensions[1],
                              Dimensions( M.generators )[2] ] );

      # Make clean.
      EXEC( "rm ", name, " ", name3 );
      Unbind( MeatAxe.SplitDimensions );

    fi;

    # Return the basis record.
    return rec( vectors             := vectors,
                isBasis             := true,
                isSemiEchelonBasis  := true,
                structure           := M,
                operations          := SemiEchelonBasisMeatAxeModuleOps );
    end;

#############################################################################
##
#F  MeatAxeModuleOps.Basis( <M> [,<vectors>] )
##
##  We do *not* allow to input basis vectors explicitly (as second argument),
##  so we may assume that every basis except the standard basis is
##  semi-echelonized.
##
MeatAxeModuleOps.Basis := MeatAxeModuleOps.SemiEchelonBasis;

#############################################################################
##
#F  MeatAxeModuleOps.StandardBasis( <M> )
#F  MeatAxeModuleOps.StandardBasis( <M>, <seed> )
##
##  returns the standard basis of the module <M>.
##  When the vector <v> is processed, the generators of '<M>.ring' are
##  applied to <v> in turn, and if the image is linear independent of the
##  basis vectors found up to this time, it is added to the basis.
##
##  If no list <seedvectors> is given, '<M>.generators' is taken instead.
##
MeatAxeModuleOps.StandardBasis := function( arg )

    local M,         # the module
          seed,      # list of seed vectors
          vectors,   # file for the standard basis vectors
          gens,      # generators of the acting ring
          file;      # intermediate file for output message

    # Check and get the arguments
    if Length( arg ) < 1 or Length( arg ) > 2 or not IsModule( arg[1] ) then
      Error( "usage: StandardBasis(<M>) resp. StandardBasis(<M>,<seedvecs>)" );
    fi;

    M:= arg[1];

    if Length( arg ) = 2 then
      seed:= arg[2];
    else
      if not IsBound( M.generators ) then
        Error( "no generators present" );
      fi;
      seed:= M.generators;
    fi;

    # Compute the standard basis.
    gens:= Parent( M ).ring.generators;
    vectors:= MeatAxe.TmpName( MeatAxe.files );
    file:= TmpName();
    MeatAxe.EXEC( "zsb -G -Q -g ", String( Length( gens ) ), " ",
          gens[1].gensname, " ", seed.file, " ", vectors, " > ", file );

    Read( file );
    if not IsBound( MeatAxe.Dimension ) then
      Error( "output not readable" );
    fi;

    vectors:= MeatAxeMat( vectors, seed.field,
                          [ MeatAxe.Dimension, seed.dimensions[2] ] );

    # Make clean.
    EXEC( "rm ", file );
    Unbind( MeatAxe.Dimension );
    Unbind( MeatAxe.CyclicDimensions );

    # Return the standard basis record.
    return rec( vectors          := vectors,
                seed             := seed,
                isBasis          := true,
                isStandardBasis  := true,
                structure        := M,
                operations       := StandardBasisMeatAxeModuleOps );
    end;

#############################################################################
##
#F  MeatAxeModuleOps.IsIrreducible( <M> )
##
MeatAxeModuleOps.IsIrreducible := function( M )
    local store;   # store the value of 'InfoMeatAxe'
    store:= InfoMeatAxe;
    InfoMeatAxe:= Ignore;
    CompositionFactors( M );
    InfoMeatAxe:= store;
    return Length( M.compositionFactors ) = 1;
    end;

#############################################################################
##
#F  MeatAxeModuleOps.IsAbsolutelyIrreducible( <M> )
##
MeatAxeModuleOps.IsAbsolutelyIrreducible := function( M )
    return IsIrreducible( M ) and
           M.compositionFactors[1].sf = 1;
    end;

#############################################################################
##
#F  MeatAxeModuleOps.IsEquivalent( <M1>, <M2> )
##
MeatAxeModuleOps.IsEquivalent := function( M1, M2 )

    local dim1, dim2,     # dimensions of modules
          eq,             # are modules equal, result
          f,
          dirsum, 
          name,           # name of the directory for the matrices
          d,
          gens,
          module,
          comp,
          i,
          store;          # store value of 'InfoMeatAxe'
          
    if    Length( Parent( M1 ).ring.generators )
       <> Length( Parent( M2 ).ring.generators ) then
      Error( "acting rings have different numbers of generators" );
    fi;

    dim1:= Dimension( M1 );
    dim2:= Dimension( M2 );

    if dim1 <> dim2 then
      eq:= false;
    else

      # first get the equivalent natural modules if necessary
      if dim1 <> Dimensions( Parent( M1 ).ring.generators[1] )[1] then
        M1:= NaturalModule( Operation( Parent( M1 ).ring, M1 ) );
      fi;
      if dim2 <> Dimensions( Parent( M2 ).ring.generators[1] )[1] then
        M2:= NaturalModule( Operation( Parent( M2 ).ring, M2 ) );
      fi;

      # construct a new directory
      name:= MeatAxe.TmpName( MeatAxe.dirs );

      # form the direct sum of matrices
      f:= Parent( M1 ).ring.field;
      dirsum:= [];
      for i in [ 1 .. Length( Parent( M1 ).ring.generators ) ] do
        Add( dirsum, Concatenation( name, "/g.", String(i) ) );
        MeatAxe.EXEC( "zpt -r 2 -c 2 ", name, "/g.", i, " ",
              M1.ring.generators[i].file, " - - ",
              M2.ring.generators[i].file );
      od;
      d:= [ dim1+dim2, dim1+dim2 ];
      gens:= List( dirsum, x -> MeatAxeMat( x, f, d ) );
      module:= NaturalModule( Algebra( f, gens ) );

      # compute the composition factors of the direct sum
      store:= InfoMeatAxe;
      InfoMeatAxe:= Ignore;
      comp:= Set( CompositionFactors( module ) );
      InfoMeatAxe:= store;

      EXEC( "rm -r ", name );

      if Length( comp ) = 1 and Dimension( comp[1] ) = dim1 then
        eq:= true;
      elif ForAny( comp, x -> Dimension(x) = dim1 ) then
        eq:= false;
      else
        Error( "sorry, modules are both reducible" );
      fi;

    fi;

    return eq; 
    end;

#############################################################################
##
#F  MeatAxeModuleOps.SplittingField( <M> )
##
MeatAxeModuleOps.SplittingField := function( M )
    local index,
          F;
    index:= Lcm( List( CompositionFactors( M ), x -> x.sf ) );
    if index = 1 then
      F:= M.ring.field;
    else
      F:= GF( M.ring.field.char ^ ( M.ring.field.dimension * index ) );
    fi;
    return F;
    end;

#############################################################################
##
#F  MeatAxeModuleOps.Lattice( <M> )
##
##  returns the submodule lattice of the {\MeatAxe} module <M>.
##  *Note* that it is *not* decomposed into direct summands, in order to
##  allow 'genmod' to get submodule generators.
##  The function 'LatticeSummands' computes the list of direct summands of
##  the lattice.
##
##  (uses 'vbd' script of Klaus Lux)
##
MeatAxeModuleOps.Lattice := function( M )

    local i,
          G,
          name,
          store,
          latt;

    if Dimension(M) <> Dimensions( Parent(M).ring.generators[1] )[1] then
      # submodule, first construct the equivalent natural module
      M:= NaturalModule( Operation( Parent( M ).ring, M ) );
    fi;

    G:= Parent( M ).ring;
    name:= G.generators[1].gensname;
    store:= InfoMeatAxe;
    InfoMeatAxe:= Ignore;
    CompositionFactors( M );
    InfoMeatAxe:= store;

    MeatAxe.EXEC( "pwkond ",   name, " >> ", name, ".log" );
    MeatAxe.EXEC( "mkcycl ",   name, " >> ", name, ".log" );
    MeatAxe.EXEC( "mkinc ",    name, " >> ", name, ".log" );
    MeatAxe.EXEC( "mkdotl ",   name, " >> ", name, ".log" );
    MeatAxe.EXEC( "mksub ",    name, " >> ", name, ".log" );

    Read( ConcatenationString( name, ".lat" ) );
    latt:= MeatAxe.Lattice;
    Unbind( MeatAxe.Lattice );

    latt:= rec( dimensions:= List( latt, z -> z[1] ),
                maxes:= List( latt, z -> List( z[2], y -> y[1] ) ),
                weights:= List( latt, z -> List( z[2], y -> y[2] ) ),
                module:= M );
    latt.XGAP:= List( [ 1 .. Length( latt.maxes ) ],
                      y -> [ y, latt.dimensions[y], latt.maxes[y] ] );

    return latt;

    end;

#############################################################################
##
#F  LatticeSummands( <M> )
##
##  returns the list of direct summands of the submodule lattice of the
##  {\MeatAxe} module <M>.
##
##  (uses 'vbd' script of Klaus Lux)
##
LatticeSummands := function( M )

    local i,
          G,
          name,
          store,
          latt;

    if Dimension(M) <> Dimensions( Parent(M).ring.generators[1] )[1] then
      # submodule, first construct the equivalent natural module
      M:= NaturalModule( Operation( Parent( M ).ring, M ) );
    fi;

    G:= Parent( M ).ring;
    name:= G.generators[1].gensname;
    store:= InfoMeatAxe;
    InfoMeatAxe:= Ignore;
    CompositionFactors( M );
    InfoMeatAxe:= store;

    MeatAxe.EXEC( "pwkond ",   name, " >> ", name, ".log" );
    MeatAxe.EXEC( "mkcycl ",   name, " >> ", name, ".log" );
    MeatAxe.EXEC( "mkinc ",    name, " >> ", name, ".log" );
    MeatAxe.EXEC( "mkdotl ",   name, " >> ", name, ".log" );
    MeatAxe.EXEC( "mksub -b ", name, " >> ", name, ".log" );

    latt:= [];
    Unbind( MeatAxe.Lattice );
    i:= 1;
    while READ( ConcatenationString( name, ".lat.", String( i ) ) ) do
      Add( latt, MeatAxe.Lattice );
      i:= i+1;
    od;
    Unbind( MeatAxe.Lattice );

    latt:= List( latt,
       x -> rec( dimensions:= List( x, z -> z[1] ),
                 maxes:= List( x, z -> List( z[2], y -> y[1] ) ),
                 weights:= List( x, z -> List( z[2], y -> y[2] ) ),
                 module:= M ) );
    for i in latt do
      i.XGAP:= List( [ 1 .. Length( i.maxes ) ],
                     y -> [ y, i.dimensions[y], i.maxes[y] ] );
    od;

    return latt;

    end;

#############################################################################
##
#F  GeneratorsSubmodule( <lattice>, <nr> ) . . . list of submodule generators
##
GeneratorsSubmodule := function( lattice, nr )

    local module,
          ring,
          name,
          result;

    module:= lattice.module;
    ring:= module.ring;
    name:= ring.generators[1].gensname;

    if nr = 1 then

      # trivial submodule will not be constructed by 'genmod'
      result:= [];

    elif nr <= Length( lattice.dimensions ) then

      # 'genmod' will put the module generators
      # into the file 'name'.s'$<nr>-1$'
      MeatAxe.EXEC( "genmod -Q ", name, " ", nr-1 );
      name:= Concatenation( name, ".s", String( nr-1 ) );
      result:= MeatAxeMat( name, ring.field,
                           [ lattice.dimensions[ nr ],
                             Dimensions( module.generators )[2] ] );

    else
      Error( "not so many submodules" );
    fi;

    return result;

    end;

#############################################################################
##
#F  GeneratorsSubmodules( <M> ) . . . . . .  list of submodules of module <M>
##
##  <M> must be a list of matrices or a *natural* module.
##
GeneratorsSubmodules := function( M )

    local F,
          dim,
          gensname,
          i,
          submodulegens,
          latt;

    if IsList( M ) then
#T das soll gehen ?

      if IsList( M[1] ) then

        F:= Field( Flat( M ) );
        dim:= Dimensions( M[1] );
        gensname:= MeatAxe.TmpName();
        M:= List( [ 1 .. Length( M ) ],
                  x -> MeatAxeMat( M[x], F,
                              Concatenation( gensname, ".", String(x) ) ) );
        for i in M do
          i.gensname:= gensname;
        od;
        
        submodulegens:= List( GeneratorsSubmodules( M ), GapObject );
        MeatAxe.Unbind( M );

      else

        # list of {\MeatAxe} matrices, compute the submodule lattice
        latt:= Lattice( NaturalModule( M ) );

        # compute the list of submodule generators
        submodulegens:= GeneratorsSubmodules(
                             NaturalModule( Algebra( M[1].field, M ) ) );

      fi;

    elif IsModule( M ) then

      latt:= Lattice( M );

      # compute the list of submodule generators
      submodulegens:= List( [ 1 .. Length( latt.dimensions ) ],
                            x -> GeneratorsSubmodule( latt, x ) );

    else
      Error( "Submodules(<M>) where <M> is a module or list of generators" );
    fi;

    return submodulegens;
    end;

#############################################################################
##
#V  MeatAxeFactorModuleOps
##
MeatAxeFactorModuleOps := OperationsRecord( "MeatAxeFactorModuleOps",
                                            FactorModuleOps );

#############################################################################
##
#F  MeatAxeFactorModuleOps.GapObject( <M> )
##
MeatAxeFactorModuleOps.GapObject := function( M )
    if IsBound( M.parent ) then
      return Submodule( GapObject( Parent( M ) ),
                        GapObject( Generators( M ) ) );
    else
      return GapObject( M.factorNum ) / GapObject( M.factorDen );
    fi;
    end;

#############################################################################
##
#F  MeatAxeModuleOps.FactorModule( <M>, <U> )
#F  MeatAxeModuleOps.\/( <M>, <U> )
##
MeatAxeModuleOps.\/ := function( M, U )
#T warum spezialis.?
    local F;
    F:= ModuleOps.\/( M, U );
    F.operations:= MeatAxeFactorModuleOps;
    return F;
    end;

#############################################################################
##
#F  MeatAxeModuleOps.FixedSubmodule( <M> ) . submodule of fixed points in <M>
##
##  calculate intersection of nullspaces of <g>-1 for all generators <g>
##
MeatAxeModuleOps.FixedSubmodule := function( M )

    local i,           # loops over generators
          mat,         # one generator minus identity
          fixed,       # basis of fixed space
          nullsp,      # fixed space of one generator
          newfixed,    # fixed space also for next generator
          S;           # the submodule, result

    mat:= M.ring.generators[1] - One( M.ring );
    fixed:= BaseNullspace( mat );
    MeatAxe.Unbind( mat );

    if not IsNaturalModule( M ) then
      newfixed:= MeatAxeMatOps.SumIntersectionSpaces( fixed, Base( M ) );
      MeatAxe.Unbind( newfixed[1], fixed );
      fixed:= newfixed[2];
    fi;

    i:= 2;
    while     i <= Length( M.ring.generators )
          and Dimensions( fixed )[1] > 0 do

      mat:= M.ring.generators[i] - One( M.ring ); 
      nullsp:= BaseNullspace( mat );
      newfixed:= MeatAxeMatOps.SumIntersectionSpaces( fixed, nullsp );
      MeatAxe.Unbind( newfixed[1], fixed, nullsp, mat );
      fixed:= newfixed[2];
      i:= i+1;

    od;
    S:= Submodule( Parent( M ), fixed );
    S.basis:= Basis( S, fixed );
    return S;
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
