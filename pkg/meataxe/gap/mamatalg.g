if not IsBound(MetaAxe) then MeatAxe:=MeatAxe;fi;
#############################################################################
##
#A  mamatalg.g                  GAP Library                     Thomas Breuer
##
#A  @(#)$Id: mamatalg.g,v 1.1.1.1 1996/12/11 12:40:13 werner Exp $
##
#Y  Copyright 1994-1995,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains the functions for {\Meataxe} matrices.
##
##  A {\MeatAxe} matrix algebra is a record with components
##
#H  $Log: mamatalg.g,v $
#H  Revision 1.1.1.1  1996/12/11 12:40:13  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 3.2  1996/11/13 15:38:41  sam
#H  fixed call of 'Algebra' in 'Operation'
#H
#H  Revision 3.1  1994/05/19 14:05:19  sam
#H  changes for GAP-3.4
#H
#H  Revision 3.0  1994/03/28  15:31:13  sam
#H  initial revision under RCS
#H
#H
##

#############################################################################
##
#F  InfoMeatAxe( ... ) . . . . . . . info function for the {\MeatAxe} package
##
    if not IsBound( InfoMeatAxe ) then InfoMeatAxe:= Ignore; fi;

#############################################################################
##
#V  MeatAxeMatAlgebraOps
##
MeatAxeMatAlgebraOps := OperationsRecord( "MeatAxeMatAlgebraOps",
                                          UnitalMatAlgebraOps );

#############################################################################
##
#V  MeatAxeMatGroupOps
##
MeatAxeMatGroupOps := OperationsRecord( "MeatAxeMatGroupOps", MatGroupOps );

#############################################################################
##
#F  MeatAxeMatAlgebraOps.GapObject( <A> )
##
MeatAxeMatAlgebraOps.GapObject := function( A )
    return Algebra( A.field, List( A.generators, GapObject ) );
    end;

#############################################################################
##
#F  MeatAxeMatGroupOps.GapObject( <G> )
##
MeatAxeMatGroupOps.GapObject := function( G )
    return Group( List( G.generators, GapObject ), GapObject( G.identity ) );
    end;

#############################################################################
##
#F  MeatAxeMatGroupOps.Order
##
MeatAxeMatGroupOps.Order := function( MeatAxeMatrices, mtxmat )
    return OrderMeatAxeMat( mtxmat );
    end;

#############################################################################
##
#F  MeatAxeMatAlgebraOps.\in
##
MeatAxeMatAlgebraOps.\in := function( A, a )
    if IsMeatAxeMat( a ) and a.structure = A then
      return true;
    else
      return DomainOps.\in( A, a );
    fi;
    end;

#############################################################################
##
#F  MeatAxeMatAlgebraOps.Order
##
MeatAxeMatAlgebraOps.Order := function( MeatAxeMatrices, mtxmat )
    return OrderMeatAxeMat( mtxmat );
    end;

#############################################################################

MeatAxeMatAlgebraOps.Algebra := MeatAxeMatricesOps.Algebra;


MeatAxeMatGroupOps.\= := function( G1, G2 )
    if Group( G1 ) then
      if IsGroup( G2 ) then
        return   List( G1.generators, x -> x.file )
               = List( G2.generators, x -> x.file );
      else
        return false;
      fi;
    else
      return false;
    fi;
    end;

MeatAxeMatGroupOps.\< := function( G1, G2 )
    if IsGroup( G1 ) then
      if IsGroup( G2 ) then
        return   List( G1.generators, x -> x.file )
               < List( G2.generators, x -> x.file );
      else
        return false;
      fi;
    else
      return true;
    fi;
    end;

MeatAxeMatAlgebraOps.\= := function( G1, G2 )
    if IsAlgebra( G1 ) then
      if IsAlgebra( G2 ) then
        return   List( G1.generators, x -> x.file )
               = List( G2.generators, x -> x.file );
      else
        return false;
      fi;
    else
      return false;
    fi;
    end;

MeatAxeMatAlgebraOps.\< := function( G1, G2 )
    if IsAlgebra( G1 ) then
      if IsAlgebra( G2 ) then
        return   List( G1.generators, x -> x.file )
               < List( G2.generators, x -> x.file );
      else
        return false;
      fi;
    else
      return true;
    fi;
    end;

#############################################################################
##
#F  MeatAxeMatGroupOps.Operation( <H>, <M>, <opr> )
##
##  for a group <H> of {\MeatAxe} matrices, a {\MeatAxe} module <M>, and an
##  operation <opr>, returns a permutation group with action on the points
##  equivalent to that of <H> on the orbits of the generators of <M>.
##
##  <H> must be a subgroup of the algebra '<M>.ring'.
##  (This is *not* checked.)
##
##  If <M> is a quotient module, it must be a quotient of the entire space.
##
MeatAxeMatGroupOps.Operation := function( H, M, opr )

    local  gens,     # generators of the result
           dir,      # directory containing the generators of the result
           i,        # loops over generators of 'H'
           name,     # temporary file name
           name2,    # temporary file name
           matrix,   # intermediate {\MeatAxe} matrix
           g;        # one generator of the result

    if IsFactorModule( M ) then
      H:= Algebra( H.generators[1].field, H.generators );
      M:= NaturalModule( Operation( H, M ) );
      H:= Group( M.ring.generators, M.generators );
    fi;

    # make a new directory where the permutations are stored
    dir:= MeatAxe.TmpName( MeatAxe.dirs );
    name:= Concatenation( dir, "/g" );

    # compute the permutation actions of the generators
    MeatAxe.EXEC( "zvp -g ", Length( H.generators ), " ",
          H.generators[1].gensname, " ", M.generators.file, " ",
          name, " > /dev/null" );

    # read the permutation generators into {\GAP}
    gens:= [];
    name2:= TmpName();
    for i in [ 1 .. Length( H.generators ) ] do
      Unbind( MeatAxe.Perms );
      MeatAxe.EXEC( "zpr -G ", name, ".", String(i), " > ", name2 );
      Read( name2 );
      Add( gens, MeatAxe.Perms[1] );
    od;

    # make clean
    Unbind( MeatAxe.Perms );
    EXEC( "rm ", name2 );
    MeatAxe.Unbind( dir );

    # return the permutation group
    return Group( gens, () );
    end;

#############################################################################
##
#F  MeatAxeMatAlgebraOps.Operation( <H>, <M>, <opr> )
##
##  for an algebra <H> of {\MeatAxe} matrices, a {\MeatAxe} module <M>, and
##  an operation <opr>, returns a {\MeatAxe} matrix algebra of same dimension
##  as <M> with action on its natural module equivalent to that of <H> on
##  <M>.
##
##  <H> must be a subalgebra of '<M>.ring'.  (This is *not* checked.)
##
#T  If <M> is a quotient module, it must be a quotient of the entire space.
##
MeatAxeMatAlgebraOps.Operation := function( H, M, opr )

    local  gens,     # generators of the result
           dim,      # dimension of the module
           f,        # field
           dir,      # directory containing the generators of the result
           i,        # loops over generators of 'H'
           name,     # temporary file name
           matrix,   # intermediate {\MeatAxe} matrix
           result;   # result algebra

    gens:= [];
    dim:= Dimension( M );

    if dim = 0 then Error( "no 0x0 matrices" ); fi;

    dir:= MeatAxe.TmpName( MeatAxe.dirs );

    if IsFactorModule( M ) then

      f:= M.factorNum.ring.field;
      for i in [ 1 .. Length( H.generators ) ] do
        name:= Concatenation( dir, "/g.", String( i ) );
        MeatAxe.EXEC( "zqt -i ", Base( M.factorDen ).file , " ",
              H.generators[i].file, " ", name );
        Add( gens, MeatAxeMat( name, f, [ dim, dim ] ) );
      od;

    else

      f:= H.field;
      for i in [ 1 .. Length( H.generators ) ] do
        name:= Concatenation( dir, "/g.", String( i ) );
        matrix:= Base( M ) * H.generators[i];
        MeatAxe.EXEC( "zcl ", Base( M ).file, " ", matrix.file,
              " /dev/null ", name );
        MeatAxe.Unbind( matrix );
        Add( gens, MeatAxeMat( name, f, [ dim, dim ] ) );
      od;

    fi;

    for i in [ 1 .. Length( gens ) ] do
      gens[i].gensname:= Concatenation( dir, "/g" );
      gens[i].abstract:= H.generators[i].abstract;
    od;

    result:= Algebra( f, gens, 0*gens[1] );

    # Add the information needed to setup the operation homomorphism.
    result.operation:= rec( genimages    := gens,
                            genpreimages := H.generators );
#T module info?

    # Return the result.
    return result;
    end;

#############################################################################
##
#F  MeatAxeMatAlgebraOps.Fingerprint( <A>, <list> ) . . . . . . . fingerprint
##
MeatAxeMatAlgebraOps.Fingerprint := function( A, range )
    local name,     # temporary file name
          fp;       # fingerprint

    name:= TmpName();
    if range = "standard" then
    
      MeatAxe.EXEC( "zsm -G fp ", A.generators[1].file, " ",
            A.generators[2].file, " > ", name );
 
    else

      MeatAxe.EXEC( "zsm -G -g ", String( Length( A.generators ) ),
                          " fp", range[1], "-", range[ Length(range) ],
                          " ", A.generators[1].gensname, " > ", name );

    fi;

    Unbind( MeatAxe.Fingerprint ); 
    Read( name );
    if not IsBound( MeatAxe.Fingerprint ) then
      Error( "not a group for fingerprint" );
    fi;
    fp:= MeatAxe.Fingerprint;
  
    # Make clean.
    Unbind( MeatAxe.Fingerprint ); 
    EXEC( "rm ", name );

    # Return the result.
    return fp;
    end;

#############################################################################
##
#F  MeatAxeMatGroupOps.Dual( <matgrp> )
##
MeatAxeMatGroupOps.Dual := function( matgrp )

    local dir,     # directory for generators of the dual group
          gens,    # generators for the dual group
          names,   # store files of inverse matrices
          i,       # loops over generators
          name;    # one generator file

    # new directory for the generators of the dual group
    dir:= MeatAxe.TmpName( MeatAxe.dirs );

    # inverses of the generators; store names
    gens:= List( matgrp.generators, x -> x^-1 );
    names:= List( gens, x -> x.file );

    # transposed matrices of the inverses, in new directory
    for i in [ 1 .. Length( gens ) ] do
      name:= Concatenation( dir, "/g.", String( i ) );
      MeatAxe.EXEC( "ztr ", gens[i].file, " ", name );
      gens[i]:= MeatAxeMat( name, gens[i].field, gens[i].dimensions );
      gens[i].gensname:= Concatenation( dir, "/g" );
    od;

    # remove inverses of original generators
    MeatAxe.Unbind( names );

    # return the dual group
    return Group( gens, matgrp.identity );
    end;
   
#############################################################################
##
#F  MeatAxeMatAlgebraOps.Dual( <matalg> )
##
MeatAxeMatAlgebraOps.Dual := function( matalg )
    return Algebra( matalg.field,
             Dual( Group( matalg.generators, matalg.identity ) ).generators,
             matalg.one );
    end;

#############################################################################
##
#F  MeatAxeMatAlgebraOps.Module( <matalg>, <gens>, <dummy> )
##
##  is a module generated by the list <gens> (which must be a {\MeatAxe}
##  matrix) and acted on by the {\MeatAxe} matrix algebra <matalg>
##
MeatAxeMatAlgebraOps.Module := function( R, gens, dummy )
    local M;

    # let the default function do the main work
    M:= RingOps.Module( R, gens, dummy );

    # overlay the right operations record
    M.operations:= MeatAxeModuleOps;

    # return the module
    return M;
    end;

#############################################################################
##
#F  MeatAxeMatAlgebraOps.NaturalModule( <matalg> )
##
MeatAxeMatAlgebraOps.NaturalModule := function( matalg )
    local gens, M;

    gens:= ShallowCopy( One( matalg ) );
    Unbind( gens.abstract );
    M:= Module( matalg, gens );
#T  M:= Module( matalg, gens, "basis" );
#T mit "basis" 
    M.basis:= SemiEchelonBasis( M, gens );
    M.isNaturalModule:= true;
    M.operations:= MeatAxeModuleOps;
    return M;
    end;

#############################################################################
##
#F  RandomOrders(<G>)  . . .  random element orders in group (or algebra) <G>
##
RandomOrders := function( G )
    local name,     # temporary file name
          ro;       # list of random orders

    if not ( IsGroup( G ) or IsAlgebra( G ) )
       or not ForAll( G.generators, IsMeatAxeMat ) then
      Error( "works for MeatAxe matrix groups only" );
    fi;
    name := TmpName();
    MeatAxe.EXEC( "zsm -G -g ", Length( G.generators ), " ro ",
          G.generators[1].gensname, " > ", name );

    Unbind( MeatAxe.RandomOrders ); 
    Read( name );
    ro:= MeatAxe.RandomOrders;
    Unbind( MeatAxe.RandomOrders ); 
  
    EXEC( "rm ", name );

    return ro;
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
