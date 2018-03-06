if not IsBound(MetaAxe) then MeatAxe:=MeatAxe;fi;
#############################################################################
##
#A  mamatrix.g                  GAP Library                     Thomas Breuer
##
#A  @(#)$Id: mamatrix.g,v 1.3 1998/02/27 10:31:21 gap Exp $
##
#Y  Copyright 1994-1995,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains the functions for {\Meataxe} matrices.
##
##  A {\MeatAxe} matrix is a record with components
##
##  'isMeatAxeMat'  : always true
##  'isMatrix'      : always true
##  'domain'        : 'MeatAxeMatrices'
##  'file'          : the name of the file that contains the matrix in
##                    {\MeatAxe} format
##  'field'         : the (finite) field the matrix is viewed over
##  'abstract'      : (optional) an element of a free algebra representing
##                    the construction of the matrix in terms of generators
##  'dimensions'    : list containing the numbers of rows and columns
##  'operations'    : 'MeatAxeMatOps'
##
#H  $Log: mamatrix.g,v $
#H  Revision 1.3  1998/02/27 10:31:21  gap
#H  Fix by Thomas. AH
#H
#H  Revision 1.2  1997/01/21 13:08:34  gap
#H  replaced 'Add' by 'AddSet' twice
#H
#H  Revision 1.1.1.1  1996/12/11 12:40:13  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 3.6  1996/12/06 16:33:45  sam
#H  fixed that a file name is added to 'MeatAxe.files'
#H
#H  Revision 3.5  1996/12/03 15:46:28  sam
#H  fixed *SILLY* bug in conjugation of matrices
#H
#H  Revision 3.4  1996/11/13 15:37:11  sam
#H  fixed 'Trace' call, added '-Q' option for 'zcv'
#H
#H  Revision 3.3  1994/08/23 11:06:47  sam
#H  replaced 'PrintZEVPols' by 'StringZEVPols',
#H  improved 'OrderMeatAxeMat',
#H  changed 'zsm mw' call
#H
#H  Revision 3.2  1994/06/14  10:35:00  sam
#H  improved 'InvariantForm' and 'Display'
#H
#H  Revision 3.1  1994/05/19  14:05:19  sam
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

############################################################################
##
#V  MeatAxeMatricesOps . . . . . . operations record for {\MeatAxe} matrices
##
MeatAxeMatricesOps := OperationsRecord( "MeatAxeMatricesOps",
                                        MatricesOps );

#############################################################################
##
#V  MeatAxeMatOps . . . . . . . . . operations record for a {\MeatAxe} matrix
##
MeatAxeMatOps := OperationsRecord( "MeatAxeMatOps" );

############################################################################
##
#V  MeatAxeMatrices  . . . . . . . . . . . domain of all {\MeatAxe} matrices
##
MeatAxeMatrices := rec( isDomain  := true,
                        name      := "MeatAxeMatrices",
                        isFinite  := false,
                        size      := "infinity",
                        operations:= MeatAxeMatricesOps );

##############################################################################
##
#F  PrintMeatAxeInput( <obj>, <Fsize>, <dim> )
##
##  prints the {\MeatAxe} ASCII input of the matrix or permutation <obj>,
##  over the field of size <Fsize>, with dimensions <dim> (a list with numbers
##  of rows and columns).
##
PrintMeatAxeInput := function( obj, Fsize, dim )

    local fflist,    # 'FFList' list for field of size 'Fsize'
          i,         # loop over rows
          j,         # loop over columns
          zero,      # string corresponding to zero row of a block
          iblocks,   # nonzero blocks in i-th row block
          pos,
          cols,
          ii,        # loop over row blocks
          jj;        # loop over column blocks

    # step 1.
    # initialize the matrix with
    # mode, field size, number of rows, number of columns
    # (free format, i.e., separated by ' ')

    fflist:= FFLISTS[ Fsize ];

    if IsPerm( obj ) then
      Print( "2 ", Fsize, " ", dim[1], " ", dim[2], "\n" );
    else
      Print( "6 ", Fsize, " ", dim[1], " ", dim[2], "\n" );
    fi;

    # step 2.
    # print the entries of the matrix
    if IsMat( obj ) then

      # ordinary matrix

      if IsInt( obj[1][1] ) then

        for i in [ 1 .. dim[1] ] do
          for j in [ 1 .. dim[2] ] do
            Print( obj[i][j], "\n" );
          od;
        od;

      else

        for i in [ 1 .. dim[1] ] do
          for j in [ 1 .. dim[2] ] do
            Print( Position( fflist, obj[i][j] ) - 1, "\n" );
          od;
        od;

      fi;

    elif not IsPerm( obj ) then

      # block matrix

      zero:= Concatenation( List( [ 1 .. obj.shorz ], x -> "0\n" ) );

      for i in [ 1 .. obj.vert ] do

        # for all row blocks
        iblocks:= Filtered( obj.mat, x -> x[1] = i );
        cols:= List( iblocks, x -> x[2] );

        for ii in [ 1 .. obj.svert ] do
          # for all rows in the block
          for j in [ 1 .. obj.horz ] do
            # for all column blocks
            pos:= Position( cols, j );
            if pos <> false then
              for jj in iblocks[ pos ][3][ ii ] do
                Print( Position( fflist, jj ) - 1, "\n" );
              od;
            else
              Print( zero );
            fi;
          od;
        od;
      od;

    else

      # permutation
      for i in [ 1 .. dim[1] ] do Print( i^obj, "\n" ); od;

    fi;
    end;

############################################################################
##
#F  IsMeatAxeMat( obj )  . . . . . . . . . . .  is <obj> a {\MeatAxe} matrix
##
IsMeatAxeMat := function( obj )
    return IsRec( obj ) 
           and IsBound( obj.isMeatAxeMat )
           and obj.isMeatAxeMat = true;
    end;

############################################################################
##
#F  MeatAxeMat( <mat>  [, <F> ] [, <abstract> ] [, <filename> ] )
##
##  returns a {\MeatAxe} matrix corresponding to the matrix (or block
##  matrix) <mat>, viewed over the finite field <F>, or over the field
##  containing all elements of <mat>.
##  (Of course <F> must contain all entries of <mat>,
##  or <mat> must be a matrix over the integers.)
##
##  If <mat> is already a {\MeatAxe} matrix then the call means that it
##  shall now be viewed over the field <F> which may be smaller or larger
##  than the field <mat> was viewed over.
##
##  The optional argument <abstract> is an element of a free algebra that
##  represents the matrix in terms of generators.
##
##  If the optional argument <filename> is given, the {\MeatAxe} matrix is
##  written to the file with this name.  Otherwise {\GAP} creates a
##  temporary file.
##
#F  MeatAxeMat( <perm>, <F>, <dim> [,<abstract>][,<filename>] )
##
##  does the same for a permutation <perm> that shall be converted into
##  a permutation matrix over the field <F>, with '<dim>[1]' rows and
##  '<dim>[2]' columns.
##
#F  MeatAxeMat( <file>, <F>, <dim> [, <abstract> ] )
##
##  is the {\MeatAxe} matrix stored on file <file>, viewed over the field
##  <F>, with dimensions <dim>, and representation <abstract>.
##  This may be used to make a shallow copy of a {\MeatAxe} matrix, or to
##  notify {\MeatAxe} matrices that were not produced by {\GAP}.
##
##  *Note*:  No field change is allowed here.
##
MeatAxeMat := function( arg )

    local obj,          # object the {\MeatAxe} matrix is constructed for
          F,            # finite field
          name,
          name2,
          dim,          # dimensions of the matrix
          abstract,
          mtxmat;

    # check the arguments
    if Length( arg ) < 1 or Length( arg ) > 6 then
      Error( "usage: MeatAxeMat(<mat>,[<F>][,<abstract>][,<file>])\n",
        " resp. MeatAxeMat(<perm>,<F>,<dim>[,<abstract>][,<file>])\n",
        " resp. MeatAxeMat(<file>,<F>,<dim>[,<abstract>])" );
    fi;

    # the object
    obj:= arg[1];

    # the field
    if Length( arg ) >= 2 and IsField( arg[2] ) then
      F:= arg[2]; 
    elif IsList( arg[1] ) then
      F:= Field( Flat( arg[1] ) );
    else
      Error( "field must be given in case of file name, perm or record" );
    fi;

    # the dimensions
    if    IsMat( obj )
       or IsMeatAxeMat( obj )
       or IsBlockMat( obj ) then
      dim:= Dimensions( obj );
    elif ( IsString( obj ) or IsPerm( obj ) ) and Length( arg ) >= 3 then
      dim:= arg[3];
    else
      Error( "dimensions of the matrix not known" );
    fi;

    # the file name
    if Length( arg ) > 1 and IsString( arg[ Length( arg ) ] ) then
      name:= arg[ Length( arg ) ];
      Unbind( arg[ Length( arg ) ] );
    fi;

    # abstract
    if IsFpAlgebraElement( arg[ Length( arg ) ] ) then
      abstract:= arg[ Length( arg ) ];
    fi;

    if   IsString( obj ) then

      # 'obj' is the file name of a {\MeatAxe} matrix
      mtxmat:= rec( isMeatAxeMat := true,
                    isMatrix     := true,
                    domain       := MeatAxeMatrices,
                    file         := obj,
                    field        := F,
                    dimensions   := dim,
                    operations   := MeatAxeMatOps );

    elif IsMeatAxeMat( obj ) then

      # 'obj' is already a {\MeatAxe} matrix
      # (this means we must perhaps apply 'zcf'.)

      if   obj.field = F and
           not IsBound( name ) then

        # no change of the field necessary, no file name
        # explicitly given, nothing to do

        mtxmat:= obj;

      elif Size( F ) mod Size( obj.field ) = 0 or
           Size( obj.field ) mod Size( F ) = 0 then

        # embed 'obj' into a larger  field or
        # view  'obj' over a smaller field

        if not IsBound( name ) then
          name := MeatAxe.TmpName( MeatAxe.files );
        fi;

        if F = obj.field then
          EXEC( "cp ", obj.file, " ", name );
        else

          # provide the field information file if necessary
          if not IsBound( MeatAxe.fields[ Size( F ) ] ) then
            MeatAxe.Maketab( Size( F ) );
          fi;
          MeatAxe.EXEC( "zcf ", String( Size( F ) ),
                " ", obj.file, " ", name, " > /dev/null" );

        fi;
        mtxmat:= MeatAxeMat( name, F, obj.dimensions );

      else

        Error( "fields are not compatible" );

      fi;

    else

      # we have to work, i.e., to set up the {\MeatAxe} matrix

      # provide the field information file if necessary
      if not IsBound( MeatAxe.fields[ Size( F ) ] ) then
        MeatAxe.Maketab( Size( F ) );
      fi;

      # print to the file and convert to internal format
      name2:= TmpName();
      PrintTo( name2, PrintMeatAxeInput( obj, Size( F ), dim ) );

      if   not IsBound( name ) then
        name:= MeatAxe.TmpName( MeatAxe.files );
      elif name[1] <> '/' then
        name:= Concatenation( MeatAxe.direc, "/", name );
        AddSet( MeatAxe.files, name );
      fi;
      MeatAxe.EXEC( "zcv -Q ", name2, " ", name );
      EXEC( "rm ", name2 );
      mtxmat:= MeatAxeMat( name, F, dim );

    fi;

    # if the optional argument <abstract> was given, add this component
    if IsBound( abstract ) then
      mtxmat.abstract:= abstract;
    fi;

    # return the {\MeatAxe} matrix
    return mtxmat;

    end;

#############################################################################
##
#F  MeatAxeMatOps.GapObject( <mtxmat> )
##
##  returns the {\GAP} matrix corresponding to the {\MeatAxe} matrix <mtxmat>
##
MeatAxeMatOps.GapObject := function( mtxmat )

    local name,   # temporary file name
          result; # the {\GAP} matrix, result
  
    name:= TmpName();
    MeatAxe.EXEC( "zpr -G ", mtxmat.file, " ", name );
    Read( name );
    if not IsBound( MeatAxe.Matrix ) then
      Error( "file '", mtxmat.file, "'does not contain a MeatAxe matrix" );
    fi;
    result:= MeatAxe.Matrix;
  
    EXEC( "rm ", name );
    Unbind( MeatAxe.Matrix );
  
    return result;
    end;

#############################################################################
##
#F  MeatAxeMatOps.Display( <obj>, <dummyrec> )
##
MeatAxeMatOps.Display := function( obj, dummyrec )
    if not IsBound( MeatAxe.fields[ Size( obj.field ) ] ) then
      MeatAxe.Maketab( Size( obj.field ) );
    fi;
    MeatAxe.EXEC( "zpr -G ", obj.file );
    end;

#############################################################################
##
#F  MeatAxeMatOps.\=
##
MeatAxeMatOps.\= := function( m1, m2 )

    local iseq;

    if IsMeatAxeMat( m1 ) then

      if IsMeatAxeMat( m2 ) then

        if m1.file = m2.file then
          return true;
        else
  
          # Construct a shell script if necessary
          if not IsBound( MeatAxe.IsEqual ) then
            MeatAxe.IsEqual:= MeatAxe.TmpName( MeatAxe.files );
            EXEC( "echo '#!/bin/csh -f' > ", MeatAxe.IsEqual, " ; " );
            EXEC( "echo 'cmp $1 $2 >/dev/null ; echo $status' >> ",
                                                          MeatAxe.IsEqual );
            EXEC( "chmod 777 ", MeatAxe.IsEqual );
          fi;
  
          # Compare the files without reading into {\GAP}
          iseq:= TmpName();
          EXEC( "echo 'MeatAxe.iseq:=' > ", iseq );
          EXEC( MeatAxe.IsEqual, " '", m1.file, "' '", m2.file, "' >> ", iseq );
          EXEC( "echo ';' >> ", iseq );
          Unbind( MeatAxe.iseq );
          Read( iseq );
          if MeatAxe.iseq = 1 then
            EXEC( "rm ", iseq );
            Unbind( MeatAxe.iseq );
            return false;
          else
            EXEC( "rm ", iseq );
            Unbind( MeatAxe.iseq );
            return true;
          fi;
  
        fi;

      elif IsMatrix( m2 ) and Dimensions( m1 ) = Dimensions( m2 ) then
        return GapObject( m1 ) = GapObject( m2 );
      else
        return false;
      fi;

    elif IsMeatAxeMat( m2 ) then

      if IsMatrix( m1 ) and Dimensions( m1 ) = Dimensions( m2 ) then
        return GapObject( m1 ) = GapObject( m2 );
      else
        return false;
      fi;

    else
      return GapObject( m1 ) = GapObject( m2 );
    fi;
    end;

#############################################################################
##
#F  MeatAxeMatOps.\<( <m1>, <m2> ) . . . . .  ordering of {\MeatAxe} matrices
##
MeatAxeMatOps.\< := function( m1, m2 )
    return m1.file < m2.file;
    end;

#############################################################################
##
#F  MeatAxeMatOps.\+
##
MeatAxeMatOps.\+ := function( m1, m2 )
    local name,     # file name of the result
          result,   # result
          field;    # common field of the two matrices

    if IsMeatAxeMat( m1 ) and IsMeatAxeMat( m2 ) then

      if m1.field = m2.field then

        name:= MeatAxe.TmpName( MeatAxe.files );
        MeatAxe.EXEC( "zad ", m1.file, " ", m2.file, " ", name );
        result:= MeatAxeMat( name, m1.field, m1.dimensions );

      elif m1.field.char = m2.field.char then

        field:= GF( m1.field.char ^
                           LcmInt( m1.field.degree, m2.field.degree ) );
        result:= MeatAxeMat( m1, field ) + MeatAxeMat( m2, field );

      else
        Error( "characteristic not compatible" );
      fi;

      if IsBound( m1.abstract ) and IsBound( m2.abstract ) then
        result.abstract:= m1.abstract + m2.abstract;
      fi;

    else
      result:= GapObject( m1 ) + GapObject( m2 );
    fi;
    return result;
    end;

#############################################################################
##
#F  MeatAxeMatOps.\-
##
MeatAxeMatOps.\- := function( m1, m2 )
    local name,     # file name of the result
          minusm2,  # negative of second matrix
          scalar,   # {\MeatAxe} matrix containing '-1'
          result,   # result
          field;    # common field of the two matrices

    if IsMeatAxeMat( m1 ) and IsMeatAxeMat( m2 ) then

      name:= MeatAxe.TmpName( MeatAxe.files );

      if m1.field = m2.field then

        if m1.field.char = 2 then
          MeatAxe.EXEC( "zad ", m1.file, " ", m2.file, " ", name );
        else
          scalar:= MeatAxeMat( [[-m2.field.one]], m2.field );
          minusm2:= scalar * m2;
          MeatAxe.EXEC( "zad ", m1.file, " ", minusm2.file, " ",
                name );
          MeatAxe.Unbind( scalar, minusm2 );
        fi;
        result:= MeatAxeMat( name, m1.field, m1.dimensions );

      elif m1.field.char = m2.field.char then

        field:= GF( m1.field.char ^
                           LcmInt( m1.field.degree, m2.field.degree ) );
        m1:= MeatAxeMat( m1, field );
        m2:= MeatAxeMat( m2, field );
        result:= m1 - m2;
        MeatAxe.Unbind( m1, m2 );

      else

        Error( "characteristic not compatible" );

      fi;

      if IsBound( m1.abstract ) and IsBound( m2.abstract ) then
        result.abstract:= m1.abstract - m2.abstract;
      fi;

    else
      result:= GapObject( m1 ) - GapObject( m2 );
    fi;

    return result;
    end;

#############################################################################
##
#F  MeatAxeMatOps.\*
##
MeatAxeMatOps.\* := function( m1, m2 )
    local name,     # file name of the result
          mtxmat,   # result
          field;    # common field of the two matrices

    if IsMeatAxeMat( m1 ) and IsMeatAxeMat( m2 ) then

      if m1.field = m2.field then

        name:= MeatAxe.TmpName( MeatAxe.files );
        MeatAxe.EXEC( "zmu ", m1.file, " ", m2.file, " ", name );
        mtxmat:= MeatAxeMat( name, m1.field, m1.dimensions );

        if IsBound( m1.abstract ) and IsBound( m2.abstract ) then
          mtxmat.abstract:= m1.abstract * m2.abstract;
        fi;

      elif m1.field.char = m2.field.char then

        field:= GF( m1.field.char ^
                           LcmInt( m1.field.degree, m2.field.degree ) );
        m1:= MeatAxeMat( m1, field );
        m2:= MeatAxeMat( m2, field );
        mtxmat:= m1 * m2;
        MeatAxe.Unbind( m1, m2 );

      else
        Error( "characteristic not compatible" );
      fi;

    elif ( IsInt( m1 ) or IsFFE( m1 ) ) and IsMeatAxeMat( m2 ) then

      if IsBound( m2.abstract ) then
        mtxmat:= m1 * m2.abstract;
      fi;
      m1:= MeatAxeMat( [ [ m1 * m2.field.one ] ], m2.field );
      name:= MeatAxe.TmpName( MeatAxe.files );
      MeatAxe.EXEC( "zmu ", m1.file, " ", m2.file, " ", name );
      if IsBound( m2.abstract ) then
        mtxmat:= MeatAxeMat( name, m2.field, m2.dimensions, mtxmat );
      else
        mtxmat:= MeatAxeMat( name, m2.field, m2.dimensions );
      fi;
      MeatAxe.Unbind( m1 );

    elif IsFFE( m2 ) and IsMeatAxeMat( m1 ) then

      if IsBound( m1.abstract ) then
        mtxmat:= m1.abstract * m2;
      fi;
      m2:= MeatAxeMat( [ [ m2 * m1.field.one ] ], m1.field );
      name:= MeatAxe.TmpName( MeatAxe.files );
      MeatAxe.EXEC( "zmu ", m1.file, " ", m2.file, " ", name );
      if IsBound( m1.abstract ) then
        mtxmat:= MeatAxeMat( name, m1.field, m1.dimensions, mtxmat );
      else
        mtxmat:= MeatAxeMat( name, m1.field, m1.dimensions );
      fi;
      MeatAxe.Unbind( m2 );

    else
      mtxmat:= GapObject( m1 ) * GapObject( m2 );
    fi;
    return mtxmat;
    end;

#############################################################################
##
#F  MeatAxeMatOps.\^( <m1>, <m2> )
##
##  <mtxmat> ^ <mtxmat>  conjugation
##  <mtxmat> ^ <int>     powering
##  <mtxmat> ^ <frobaut> apply Frobenius automorphism
##
MeatAxeMatOps.\^ := function( m1, m2 )
    local name,          # file name
          mtxmat;        # result

    if IsMeatAxeMat( m1 ) then

      if IsMeatAxeMat( m2 ) then

        # conjugate matrices, return m2^-1 * m1 * m2
        mtxmat:= m2^-1 * m1 * m2;

      elif IsInt( m2 ) then

        # <mat> ^ <int>
        if m2 > 1 then

          name:= MeatAxe.TmpName( MeatAxe.files );
          MeatAxe.EXEC( "zsm pwr", m2, " ", m1.file, " ", name );
          mtxmat:= MeatAxeMat( name, m1.field, m1.dimensions );

          if IsBound( m1.abstract ) then
            mtxmat.abstract:= m1.abstract ^ m2;
          fi;

        elif m2 = 1 then

          mtxmat:= m1;

        elif m2 = 0 then

          if IsBound( m1.structure ) and IsBound( m1.structure.one ) then

            # the identity matrix of the algebra/group is already stored.
            mtxmat:= m1.structure.one;

          else

            # Compute the identity matrix.
            name:= MeatAxe.TmpName( MeatAxe.files );
            MeatAxe.EXEC( "zsm -Q mw1id ", m1.file, " ", m1.file,
                  " ", name );
            mtxmat:= MeatAxeMat( name, m1.field, m1.dimensions );
            if IsBound( m1.abstract ) then
              mtxmat.abstract:= m1.abstract ^ 0;
            fi;

            # Store the identity matrix.
            if IsBound( m1.structure ) then
              m1.structure.one:= mtxmat;
            fi;

          fi;

        else

          if not IsBound( m1.inverse ) then

            # write the inverse of 'm1' to 'name'
            name:= MeatAxe.TmpName( MeatAxe.files );
            MeatAxe.EXEC( "ziv ", m1.file, " ", name );
            m1.inverse:= MeatAxeMat( name, m1.field, m1.dimensions );
            m1.inverse.inverse:= m1;

            # no chance to invert arbitrary element of a f.p. algebra,
            # but in case of a word with coefficient 1 it is allowed
            if     IsBound( m1.abstract )
               and Length( m1.abstract.coeff ) = 1
               and m1.abstract.coeff[1] = m1.field.one then
              m1.inverse.abstract:=
                   FpAlgebraElement( m1.abstract.domain,
                                     m1.abstract.coeff,
                                     [ m1.abstract.words[1] ^ -1 ] );
            fi;

          fi;

          if m2 = -1 then
            mtxmat:= m1.inverse;
          else
            name:= MeatAxe.TmpName( MeatAxe.files );
            MeatAxe.EXEC( "zsm pwr", -m2, " ", m1.inverse.file, " ", name );
            mtxmat:= MeatAxeMat( name, m1.field, m1.dimensions );
          fi;

        fi;

      elif IsFrobeniusAutomorphism( m2 )
           and m2.source.char = m1.field.char then

        name:= MeatAxe.TmpName( MeatAxe.files );
        MeatAxe.EXEC( "zfr ", m1.file, " ", name, " >/dev/null" );
        mtxmat:= MeatAxeMat( name, m1.field, m1.dimensions );

      else
        mtxmat:= GapObject( m1 ) ^ GapObject( m2 );
      fi;
    else
      mtxmat:= GapObject( m1 ) ^ GapObject( m2 );
    fi;

    return mtxmat;
    end;

#############################################################################
##
#F  MeatAxeMatOps.Base( <mtxmat> )
##
MeatAxeMatOps.Base := function( mtxmat )

    local name,       # file that will contain the basis in semi-echelon form
          name2;      # file that will contain the rank

    if not IsBound( mtxmat.base ) then

      name := MeatAxe.TmpName( MeatAxe.files );
      name2:= TmpName();
      MeatAxe.EXEC( "zef -G ", mtxmat.file, " ", name, " > ", name2 );
      Unbind( MeatAxe.Rank ); 
      Read( name2 );
      if not IsBound( MeatAxe.Rank ) then
        Error( "no MeatAxe matrix" );
      fi;
      EXEC( "rm ", name2 );
      mtxmat.base:= MeatAxeMat( name, mtxmat.field,
                                [ MeatAxe.Rank, mtxmat.dimensions[2] ] );
      Unbind( MeatAxe.Rank ); 

    fi;
    return mtxmat.base;
    end;

#############################################################################
##
#F  MeatAxeMatOps.Print . . . . . . . . . . . . .   print a {\MeatAxe} matrix
##
MeatAxeMatOps.Print := function( obj )
    Print( "MeatAxeMat( \"", obj.file, "\", ", obj.field, ", ",
           obj.dimensions );
    if IsBound( obj.abstract ) then
      Print( ", ", obj.abstract );
    fi;
    Print( " )" );
    end;

############################################################################
##
#F  MeatAxeMatricesOps.\in
##
MeatAxeMatricesOps.\in := function( mat, MeatAxeMatrices )
    return IsMeatAxeMat( mat );
    end;

#############################################################################
##
#F  MeatAxeMatOps.Rank( <mtxmat> ) . . . . . . .  rank of a {\MeatAxe} matrix
##
MeatAxeMatricesOps.Rank := function( mtxmat )
    if IsBound( mtxmat.baseNullspace ) then
      return   Dimensions( mtxmat )[1]
             - Dimensions( mtxmat.baseNullspace )[1];
    else
      return Dimensions( Base( mtxmat ) )[1];
    fi;
    end;

#############################################################################
##
#F  MeatAxeMatOps.Dimensions( <mtxmat> )
##
MeatAxeMatOps.Dimensions := function( mtxmat )
    return mtxmat.dimensions;
    end;

#############################################################################
##
#F  OrderMeatAxeMat( <mtxmat> )
##
OrderMeatAxeMat := function( mtxmat )
    local name;

    if not IsBound( mtxmat.order ) then

      name:= TmpName();
      MeatAxe.EXEC( "zor -G ", mtxmat.file, " >", name );
  
      Unbind( MeatAxe.Order );
      Read( name );
      if not IsBound( MeatAxe.Order ) then
        Error( "object is no MeatAxe matrix or permutation" );
      fi;
      EXEC( "rm ", name );
      mtxmat.order:= MeatAxe.Order;
      Unbind( MeatAxe.Order );

    fi;
    return mtxmat.order;
    end;

#############################################################################
##
#F  MeatAxeMatricesOps.Order
##
MeatAxeMatricesOps.Order := function( MeatAxeMatrices, mtxmat )
    return OrderMeatAxeMat( mtxmat );
    end;

#############################################################################
##
#F  MeatAxeMatOps.Trace( MeatAxeMatrices, <mtxmat> )
##
MeatAxeMatOps.Trace := function( MeatAxeMatrices, mtxmat )
    local name;

    name:= TmpName();
    MeatAxe.EXEC( "ztc -G ", mtxmat.file, " >", name );

    Unbind( MeatAxe.Trace );
    Read( name );
    EXEC( "rm ", name );
    if not IsBound( MeatAxe.Trace ) then
      Error( "no MeatAxe matrix" );
    fi;
    name:= MeatAxe.Trace;
    Unbind( MeatAxe.Trace );
    return name;
    end;

#############################################################################
##
#F  MeatAxeMatOps.Transposed( <mtxmat> )
##
MeatAxeMatOps.Transposed := function( mtxmat )
    local name;
    name:= MeatAxe.TmpName( MeatAxe.files );
    MeatAxe.EXEC( "ztr ", mtxmat.file, " ", name );
    return MeatAxeMat( name, mtxmat.field, mtxmat.dimensions );
    end;

#############################################################################
##
#F  MeatAxeMatOps.KroneckerProduct( <m1>, <m2> )
##
MeatAxeMatOps.KroneckerProduct := function( m1, m2 )
    local name,         # file name
          dim1, dim2;   # dimensions of the matrices
    if m1.field = m2.field then
      name:= MeatAxe.TmpName( MeatAxe.files );
      dim1:= Dimensions( m1 );
      dim2:= Dimensions( m2 );
      MeatAxe.EXEC( "zte ", m1.file, " ", m2.file, " ", name );
      return MeatAxeMat( name, m1.field,
                         [ dim1[1]*dim2[1], dim1[2]*dim2[2] ] );
    else
      Error( "<m1> and <m2> must have the same field" );
    fi;
    end;

#############################################################################
##
#F  MeatAxeMatOps.BaseNullspace( <mtxmat> )
##
MeatAxeMatOps.BaseNullspace := function( mtxmat )

    local name,       # file that will contain the basis in semi-echelon form
          name2;      # file that will contain the nullity

    name := MeatAxe.TmpName( MeatAxe.files );
    name2:= TmpName();
    MeatAxe.EXEC( "znu -G ", mtxmat.file, " ", name, " > ", name2 );

    Unbind( MeatAxe.Nullity );
    Read( name2 );
    EXEC( "rm ", name2 );
    if not IsBound( MeatAxe.Nullity ) then
      Error( "no MeatAxe matrix" );
    fi;
    name2:= MeatAxe.Nullity;

    Unbind( MeatAxe.Nullity );
    return MeatAxeMat( name, mtxmat.field, [ name2, mtxmat.dimensions[2] ] );

    end;

############################################################################
##
#F  MeatAxeMatricesOps.CharacteristicPolynomial( D, mtxmat )
##
MeatAxeMatricesOps.CharacteristicPolynomial := function( D, mtxmat )
    local name;
    if not IsBound( mtxmat.characteristicPolynomial ) then
      name:= TmpName();
      MeatAxe.EXEC( "zcp -G ", mtxmat.file, " >", name );
      Unbind( MeatAxe.CharPol );
      Read( name );
      MeatAxe.CharPol:= List( MeatAxe.CharPol,
                              x -> Polynomial( mtxmat.field, x ) );
      mtxmat.characteristicPolynomial:= Product( MeatAxe.CharPol );
      mtxmat.characteristicPolynomial.factors:= MeatAxe.CharPol;
      Unbind( MeatAxe.CharPol );
      EXEC( "rm ", name );
    fi;
    return mtxmat.characteristicPolynomial;
    end;

#############################################################################
##
#F  MeatAxeMatOps.SumIntersectionSpaces(<m1>,<m2>) . .
#F                                         sum and intersection of row spaces
##
MeatAxeMatOps.SumIntersectionSpaces := function( m1, m2 )

    local namesum,    # file of the sum
          nameint,    # file of the intersection
          name2,      # file of the output
          result;     # list with sum and intersection of spaces, result

    if   Dimensions( m1 )[1] = 0 then
      result:= [ m2, m1 ];
    elif Dimensions( m2 )[1] = 0 then
      result:= [ m1, m2 ];

    elif Dimensions( m1 )[2] = Dimensions( m2 )[2] 
       and m1.field = m2.field then

      namesum:= MeatAxe.TmpName( MeatAxe.files );
      nameint:= MeatAxe.TmpName( MeatAxe.files );
      name2  := TmpName();
      MeatAxe.EXEC( "zsi -G ", m1.file, " ", m2.file, " ", namesum,
            " ", nameint, " > ", name2 );
      Unbind( MeatAxe.SumDimension );
      Unbind( MeatAxe.IntDimension );
      Read( name2 );
      EXEC( "rm ", name2 );
      result:= [ MeatAxeMat( namesum, m1.field,
                       [ MeatAxe.SumDimension, Dimensions( m1 )[2] ] ),
                 MeatAxeMat( nameint, m1.field,
                       [ MeatAxe.IntDimension, Dimensions( m1 )[2] ] ) ];
      Unbind( MeatAxe.SumDimension );
      Unbind( MeatAxe.IntDimension );

    else
      Error( "need same field and same number of columns" );
    fi;
    return result;
    end; 

#############################################################################
##
#F  MeatAxeMatricesOps.Group( MeatAxeMatrices, <gens>,<id>)
##
##  should work by call of 'MatricesOps.Group', but there is a call of
##  'Field' ... (will also cause problems for groups of block matrices ... )
##
#T  Better get rid of the crazy 'identity' !!
##
MeatAxeMatricesOps.Group := function( MeatAxeMatrices, gens, id )
    local   F,           # field over which all generators are written
            degree,      # degree of the field
            i,
            G,
            d,           # dimensions of the generators
            g,           # one generator
            name,
            try,
            gensname;

    # check that the generators are all of the same size
    d := Dimensions( id );
    for g in gens do
      if Dimensions( g ) <> d then
        Error("<gens> must be a list of (invertible) square matrices");
      fi;
    od;

    # make sure that the generators are all written over the same field
    F:= gens[1].field;
    if not ForAll( gens, x -> x.field.char = F.char ) then
      Error( "characteristic not compatible" );
    elif not ForAll( gens, x -> x.field = F ) then

      # change the field (and put the new generators into a new directory)
      degree:= 1;
      for g in gens do
        degree:= Lcm( degree, g.field.degree );
      od;
      F:= GF( F.char ^ degree );

      name:= MeatAxe.TmpName( MeatAxe.dirs );
      gens:= Copy( gens );

      for i in [ 1 .. Length( gens ) ] do
        gens[i]:= MeatAxeMat( gens[i], F,
                              Concatenation( name, "/g.", String(i) ) );

        gens[i].gensname:= Concatenation( name, "/g" );
      od;

    fi;

    if not ForAll( gens, x -> IsBound( x.gensname ) )
       or Length( Set( List( gens, x -> x.gensname ) ) ) > 1 then

      # The generators must have names appropriate for 'zsp', 'chop'.
      # If they do not then make a new directory and copy the files there.

      try:= Position( Reversed( gens[1].file ), '.' );
      if try <> false then
        try:= Length( gens[1].file ) - try;
      fi;
      if try <> false then
        gensname:= gens[1].file{ [ 1 .. try ] };
      fi;

      i:= 2;
      while try <> false and i <= Length( gens ) do
        if    Length( gens[i].file ) < try
           or gens[i].file{ [ 1 .. try ] } <> gensname
           or gens[i].file{ [ try+2 .. Length( gens[i].file ) ] }
                       <> String(i) then
          try:= false;
        fi;
        i:= i+1;
      od;

      if try = false then

        name:= MeatAxe.TmpName( MeatAxe.dirs );
        gensname:= Concatenation( name, "/g" );

        gens:= Copy( gens );

        for i in [ 1 .. Length( gens ) ] do
          EXEC( "cp ", gens[i].file, " ", gensname, ".", i );
          gens[i].file:= Concatenation( gensname, ".", String( i ) );
          gens[i].gensname:= gensname;
          AddSet( MeatAxe.files, gens[i].file );
        od;

      else
        for i in [ 1 .. Length( gens ) ] do
          gens[i].gensname:= gensname;
        od;
      fi;

    fi;

    # make the group record
    G            := rec();
    G.isDomain   := true;
    G.isGroup    := true;
    G.identity   := id;
    G.generators := ShallowCopy( gens );

    # add the matrix group tag
    G.isMatGroup     := true;

    # add the known information
    G.dimension         := d[1];
    G.field             := id.field;

    # add the operations record
    G.operations        := MeatAxeMatGroupOps;

    # add abstract generators if not yet bound for group generators
    if ForAll( gens, x -> not IsBound( x.abstract ) ) then

      # make a new name of a generator
      name:= WordAlp( MeatAxe.alpha, Length( MeatAxe.gennames ) + 1 );
      Add( MeatAxe.gennames, name );

      G.freeAlgebra:= FreeAlgebra( id.field, Length( gens ), name );
      for i in [ 1 .. Length( gens ) ] do
        gens[i].abstract:= G.freeAlgebra.generators[i];
      od;
      G.identity.abstract:= G.freeAlgebra.one;
    fi;

    # return the group record
    return G;
    end;

#############################################################################
##
#F  MeatAxeMatricesOps.Algebra( <R>, MeatAxeMatrices, <gens>, <zero> )
##
##  *Note* that <zero> may be a dummy (equal to <R>) if <gens> is nonempty.
##
MeatAxeMatricesOps.Algebra := function( R, MeatAxeMatrices, gens, zero )

    local   i,
            G,
            d,
            g,
            name,
            try,
            gensname;

    # Check that the generators are all of the same size.
    if Length( gens ) > 0 then
      d:= Dimensions( gens[1] );
    else
      d:= Dimensions( zero );
    fi;

    for g in gens do
      if Dimensions( g ) <> d then
        Error("<gens> must be a list of square matrices");
      fi;
    od;

    # Make sure that the generators are all written over 'R'.
    if   ForAny( gens, x -> x.field.char <> R.char ) then
      Error( "characteristic not compatible" );
    elif ForAny( gens, x -> x.field <> R ) then

      # change the field (and put the new generators into a new directory)

      name:= MeatAxe.TmpName( MeatAxe.dirs );
      gens:= Copy( gens );

      for i in [ 1 .. Length( gens ) ] do
        gens[i]:= MeatAxeMat( gens[i], R,
                              Concatenation( name, "/g.", String(i) ) );

        gens[i].gensname:= Concatenation( name, "/g" );
      od;

    fi;

    if not ForAll( gens, x -> IsBound( x.gensname ) )
       or Length( Set( List( gens, x -> x.gensname ) ) ) > 1 then

      # The generators must have names appropriate for 'zsp', 'chop'.
      # If they do not then make a new directory and copy the files there.

      # The generators do already have appropriate names if and only if 
      # a '.' occurs in the file name, followed only by digits.
      try:= Position( Reversed( gens[1].file ), '.' );
      if try <> false then
        try:= Length( gens[1].file ) - try;
        if ForAny( [ try+1 .. Length( gens[1].file ) ],
                   x -> not gens[1].file[x] in "0123456789" ) then
          try:= false;
        fi;
      fi;
      if try <> false then
        gensname:= gens[1].file{ [ 1 .. try ] };
      fi;

      i:= 2;
      while try <> false and i <= Length( gens ) do
        if    Length( gens[i].file ) < try
           or gens[i].file{ [ 1 .. try ] } <> gensname
           or gens[i].file{ [ try+2 .. Length( gens[i].file ) ] }
                      <> String(i) then
          try:= false;
        fi;
        i:= i+1;
      od;

      if try = false then

        name:= MeatAxe.TmpName( MeatAxe.dirs );
        gensname:= Concatenation( name, "/g" );

        gens:= Copy( gens );

        for i in [ 1 .. Length( gens ) ] do
          EXEC( "cp ", gens[i].file, " ", gensname, ".", i );
          gens[i].file:= Concatenation( gensname, ".", String( i ) );
          gens[i].gensname:= gensname;
          AddSet( MeatAxe.files, gens[i].file );
        od;

      else
        for i in [ 1 .. Length( gens ) ] do
          gens[i].gensname:= gensname;
        od;
      fi;

    fi;

    # Make the algebra record.
    G              := rec();
    G.isDomain     := true;
    G.isAlgebra    := true;
    G.isUnitalAlgebra := true;
    G.isMatAlgebra := true;
    G.field        := R;
    G.generators   := ShallowCopy( gens );

    # Add the operations record.
    G.operations   := MeatAxeMatAlgebraOps;

    # add domain info for elements
    for i in gens do
      i.structure:= G;
    od;

    # add abstract generators if not yet bound for group generators
    if ForAll( gens, x -> not IsBound( x.abstract ) ) then

      # make a new name of a generator
      name:= WordAlp( MeatAxe.alpha, Length( MeatAxe.gennames ) + 1 );
      Add( MeatAxe.gennames, name );

      G.freeAlgebra:= FreeAlgebra( R, Length( gens ), name );
      for i in [ 1 .. Length( gens ) ] do
        gens[i].abstract:= G.freeAlgebra.generators[i];
      od;

    else
      G.freeAlgebra:= gens[1].abstract.domain;
    fi;

    # Check whether the zero element is given explicitly.
    if zero <> R then
      G.zero:= ShallowCopy( zero );
      G.zero.abstract:= G.freeAlgebra.zero;
    fi;

    # return the algebra record
    return G;
    end;

#############################################################################
##
#F  MeatAxeMatricesOps.UnitalAlgebra( <R>, MeatAxeMatrices, <gens>, <zero> )
##
MeatAxeMatricesOps.UnitalAlgebra := MeatAxeMatricesOps.Algebra;

#############################################################################
##
#F  MeatAxeMatricesOps.InvariantForm( <matrices> )
##
MeatAxeMatricesOps.InvariantForm := function( D )

    local i,          # loop variable
          F,          # field of matrix entries
          dim,        # dimensions of matrices
          pos,        #
          fp,         # fingerprints of natural module
          gensname,   # 'gensname' entry of all generators
          name,       # temporary file name
          ns,         # file for a nullspace
          exec,       # command string
          sb1,        # standard basis of 'D'
          sb2,        # standard basis of contragredient representation 'T'
          pr1, pr2,
          M,
          result;

    if ForAll( D, x -> IsBound( x.gensname ) ) and
       Length( Set( List( D, x -> x.gensname ) ) ) = 1 then
      gensname:= D[1].gensname;
    else
      gensname:= MeatAxe.TmpName();
      for i in [ 1 .. Length( D ) ] do
        EXEC( "cp ", D[i].file, " ", gensname, ".", i );
      od;
      D:= List( [ 1 .. Length( D ) ],
                x -> MeatAxeMat( Concatenation( gensname, ".", String( x ) ),
                                 D[1].field, D[1].dimensions ) );
    fi;

    # compute fingerprints of the natural module
    name:= TmpName();
    pos:= 2;
    fp:= [];
    while ForAll( fp, x -> x[2] <> 1 ) do
      InfoMeatAxe( "#I InvariantForm: looking for word of nullity 1\n" );
      MeatAxe.EXEC( "zsm -G -g ", String( Length( D ) ),
                          " fp", pos, "-", pos+9, " ",
                          gensname, " > ", name );
      Unbind( MeatAxe.Fingerprint );
      Read( name );
      fp:= MeatAxe.Fingerprint;
      Unbind( MeatAxe.Fingerprint );
      pos:= pos+10;
    od;

    # write the 'D' nullspace of the nullity 1 word to 'ns'
    ns:= TmpName();
    pos:= First( fp, x -> x[2] = 1 );
    MeatAxe.EXEC( "zsm -Q -g ", String( Length( D ) ), " mw",
               String( pos[1] ), " ", gensname, " ", name, " ", ns );

    # compute the standard basis for 'D',
    # starting with word of nullity 1, write output to 'sb1'

    sb1:= MeatAxe.TmpName();
    MeatAxe.EXEC( "zsb -G -Q -g ", String( Length( D ) ), " ",
          gensname, " ", ns, " ", sb1, " > ", name );

    # check whether the whole space is spanned
    Read( name );
    if MeatAxe.Dimension < D[1].dimensions[1] then
      Error( "matrices do not act irreducibly" );
    fi;
    Unbind( MeatAxe.Dimension );

    # the contragredient representation 'T'
    gensname:= MeatAxe.TmpName();
    for i in [ 1 .. Length( D ) ] do
      MeatAxe.EXEC( "ziv ", D[i].file, " ", name );
      MeatAxe.EXEC( "ztr ", name, " ", gensname, ".", i );
    od;

    # write the 'T' nullspace of the nullity 1 word to 'ns'
    MeatAxe.EXEC( "zsm -Q -g ", String( Length( D ) ), " mw",
               String( pos[1] ), " ", gensname, " ", name, " ", ns );

    # compute the standard basis for 'T',
    # starting with word of nullity 1, write output to 'sb2'

    sb2:= TmpName();
    MeatAxe.EXEC( "zsb -Q -g ", String( Length( D ) ), " ",
          gensname, " ", ns, " ", sb2 );

    # if 'D' and 'T' are equivalent then
    # the invariant matrix is 'M = sb1^(-1) * sb2',
    # since 'sb1 * D[i] * sb1^(-1) = sb2 * T[i] * sb2^(-1)' implies
    # that D[i] * M * D^{tr} = M'.

    MeatAxe.EXEC( "ziv ", sb1, " ", name );
    MeatAxe.EXEC( "zmu ", name, " ", sb2, " ", sb1 );

    # make clean
    EXEC( "rm ", name, " ", sb2, " ", ns );

    # check for equality
    M:= MeatAxeMat( sb1, D[1].field, D[1].dimensions );
    for i in [ 1 .. Length( D ) ] do
      pr1:= D[i] * M;
      pr2:= M * MeatAxeMat( Concatenation( gensname, ".", String( i ) ),
                            D[1].field, D[1].dimensions );
      if pr1 <> pr2 then

        # make clean
        MeatAxe.Unbind( pr1, pr2, sb1 );
        return false;

      else
        MeatAxe.Unbind( pr1, pr2 );
      fi;
    od;  

    # return the result
    return M;
    end;

###############################################################################
##
#F  ZEVPols( <F>, <ord> )
##
##  Given a finite field <F> and an element of order <ord>, this program
##  returns the possible polynomials in {\MeatAxe} notation and the
##  corresponding Brauer character values.
##
##  If the necessary finite fields are not available
##  in {\GAP}, a suitable Conway polynomial must be available in the
##  function 'ConwayPol'.
##
ZEVPols := function( F, ord )
    local ffl,       # 'FFList( F )'
          Ee,
          pf,        # prime field
          factors,l,dimext,extfieldsize,xpol,list,i,extfield,fac,
          conwaypol,fieldroot,elmroot,exp,power,res,j;
    
    if ord mod F.char = 0 then
      Error( "<ord> must be relatively prime to the characteristic of <F>" );
    fi;
    
    ffl:= FFList( F );
    
    # irreducible factors of cyclotomic polynomial over 'F'
    factors:= Factors( X(F)^ord - X(F)^0 );
    
    # splitting field for element of order 'ord'
    dimext:= F.degree * OrderMod( Size( F ), ord );
    extfieldsize:= F.char ^ dimext; 
    
    list:=[];
    
    if extfieldsize <= 65536 then
    
      # {\GAP} supports the field
      extfield:= GF( extfieldsize );
      Ee:= E( extfieldsize - 1 );
    
      for i in [1..Length(factors)] do
        list[i]:=[ 0,
                   Reversed( List( factors[i].coefficients,
                                   j -> Position( ffl, j ) - 1 ) ) ];
    
        factors[i].baseRing:= extfield;
        for fac in Factors( factors[i] ) do
          if F.char = 2 then
            list[i][1]:=   list[i][1]
                         + Ee^LogFFE(fac.coefficients[1],extfield.root);
          else
            list[i][1]:=   list[i][1]
                         - Ee^LogFFE(fac.coefficients[1],extfield.root);
          fi;
        od;
      od;
    
    else
    
      # splitting field is not available in {\GAP}
      Error( "splitting field is too large for GAP" );
    
  #    Ee:= E( ord );
  #  
  #    for i in [1..Length(factors)] do
  #      list[i]:=[ 0,
  #                 Reversed( List( factors[i].coefficients,
  #                                 j->Position(ffl,j)-1)) ];
  #  
  #      pf:= GF( F.char );
  #      conwaypol:= Polynomial( pf,
  #                    ConwayPol( F.char, dimext ) * pf.one );
  #      xpol:= X( pf );
  #      fieldroot:= PowerMod( xpol, (extfieldsize-1) / (Size(F)-1), conwaypol);
  #      elmroot:=   PowerMod( xpol, (extfieldsize-1) / ord,         conwaypol);
  #      l:=Length(factors[i].coefficients);
  #      for exp in [0..ord-1] do
  #        power:=PowerMod(elmroot,exp,conwaypol);
  #        if factors[i].coefficients[l] = F.zero then
  #          res:=0;
  #        else
  #          res:=PowerMod(fieldroot,
  #               LogFFE(factors[i].coefficients[l], F.root),conwaypol);
  #        fi;
  #        for j in [1..l-1] do
  #          if factors[i].coefficients[l-j] = F.zero then
  #                 res:=res*power;
  #          else
  #                 res:=res*power+PowerMod(fieldroot,
  #                      LogFFE(factors[i].coefficients[l-j], F.root),
  #                      conwaypol);
  #          fi;
  #          res:=res mod conwaypol;
  #        od;
  #        if res=0*xpol then
  #          list[i][1]:=list[i][1] + Ee ^ exp;
  #        fi;
  #      od;
  #    od;
    fi;
    
    # return the list of polynomials
    return list;
    end;
    
###############################################################################
##
#F  StringZEVPols( <F>, <ord> )
##
##  Given a finite field <F> and an element of order <ord>, this program
##  returns the string of an input file for the {\MeatAxe} program 'zev'.
##
##  If the necessary finite fields are not available in {\GAP},
##  a suitable Conway polynomial must be supplied.
##
StringZEVPols := function( F, ord )

    local str,     # string to be returned
          entry,   # loop over output of 'ZEVPols'
          i,       # loop over 'entry[2]'
	  j,
	  pols,
    	  names;
    
    # header
    str:= String( Size( F ) );
    Add( str, 'p' );
    Append( str, String( ord ) );
    Add( str, '\n' );
    
    # get the possible polynomials
    pols:= ZEVPols( F, ord );
    names:= rec();
    for j in [ 1 .. Length( pols ) ] do
      entry:= pols[j];

      Add( str, ' ' );
      Append( str, "MeatAxe.ZEV." );
      Append( str, String( j ) );
      names.( j ):= entry[1];
      Add( str, ' ' );
      for i in entry[2] do
        Append( str, String( i ) );
        Add( str, ' ' );
      od;
      Add( str, '\n' );
    od;
    IsString( str );
    return [str,names];
    end;

#############################################################################
##
#F  BrauerCharacterValue( <mat> )
##
BrauerCharacterValue := function( mat )

    local name, name2,   # temporary file names
          str,
          val;

    name:= TmpName();
    val:= OrderMeatAxeMat( mat );
    if val = 1 then
      val:= mat.dimensions[1];
    else

      # Compute the possible polynomials
      # (This is so complicated because we have to avoid outputs
      # of 'InfoRead' etc. inside the 'PrintTo' call.)
      str:= StringZEVPols( mat.field, val );
      PrintTo( name, str[1], "\n" );
  
      # Run 'zev' to get the eigenvalues.
      name2:= TmpName();
      MeatAxe.EXEC( "zev -G ", mat.file, " ", name, " > ", name2 );
  
      # get the output
      Unbind( MeatAxe.BrauerChar );
      MeatAxe.ZEV:=str[2];
      Read( name2 );
      val:= MeatAxe.BrauerChar;
      Unbind( MeatAxe.BrauerChar );
      Unbind(MeatAxe.ZEV);
  
      # Make clean.
      EXEC( "rm ", name, " ", name2 );

    fi;
    return val;
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
