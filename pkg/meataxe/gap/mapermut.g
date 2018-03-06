if not IsBound(MetaAxe) then MeatAxe:=MeatAxe;fi;
#############################################################################
##
#A  mapermut.g                  GAP Library                     Thomas Breuer
##
#A  @(#)$Id: mapermut.g,v 1.2 1997/01/21 13:08:36 gap Exp $
##
#Y  Copyright 1994-1995,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains the functions for {\Meataxe} permutations.
##
##  A {\MeatAxe} permutation is a record with components
##
##  'isMeatAxePerm' : always true
##  'isPermutation' : always true
##  'domain'        : 'MeatAxePermutations'
##  'file'          : the name of the file that contains the permutation in
##                    {\MeatAxe} format
##  'operations'    : 'MeatAxePermOps'
##
#H  $Log: mapermut.g,v $
#H  Revision 1.2  1997/01/21 13:08:36  gap
#H  replaced 'Add' by 'AddSet' twice
#H
#H  Revision 1.1.1.1  1996/12/11 12:40:13  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 3.2  1996/12/06 16:34:10  sam
#H  fixed that a file name is added to 'MeatAxe.files'
#H
#H  Revision 3.1  1996/11/13 16:01:07  sam
#H  changed format of MeatAxe permutations (now need largest moved point)
#H
#H  Revision 3.0  1994/06/17 16:36:18  sam
#H  initial revision under RCS
#H
##

#############################################################################
##
#F  InfoMeatAxe( ... ) . . . . . . . info function for the {\MeatAxe} package
##
    if not IsBound( InfoMeatAxe ) then InfoMeatAxe:= Ignore; fi;

#############################################################################
##
#V  MeatAxePermOps . . . . . . operations record for a {\MeatAxe} permutation
##
MeatAxePermOps := OperationsRecord( "MeatAxePermOps" );

############################################################################
##
#F  IsMeatAxePerm( obj )  . . . . . . . .  is <obj> a {\MeatAxe} permutation
##
IsMeatAxePerm := function( obj )
    return IsRec( obj ) 
           and IsBound( obj.isMeatAxePerm )
           and obj.isMeatAxePerm = true;
    end;

############################################################################
##
#V  MeatAxePermutationsOps  .  operations record for {\MeatAxe} permutations
##
MeatAxePermutationsOps := OperationsRecord( "MeatAxePermutationsOps" );

############################################################################
##
#V  MeatAxePermutations  . . . . . . . domain of all {\MeatAxe} permutations
##
MeatAxePermutations := rec( isDomain  := true,
                            name      := "MeatAxePermutations",
                            isFinite  := false,
                            size      := "infinity",
                            operations:= MeatAxePermutationsOps );

#############################################################################
##
#F  OrderMeatAxePerm( <mtxperm> )
##
OrderMeatAxePerm := function( mtxperm )

    local name;

    if not IsBound( mtxperm.order ) then

      name:= TmpName();
      MeatAxe.EXEC( "zor -G ", mtxperm.file, " >", name );
  
      Unbind( MeatAxe.Orders );
      Read( name );
      EXEC( "rm ", name );
      if not IsBound( MeatAxe.Orders ) then
        Error( "object is no no MeatAxe matrix or permutation" );
      fi;
      mtxperm.order:= MeatAxe.Orders[1];
      Unbind( MeatAxe.Orders );

    fi;
    return mtxperm.order;
    end;

#############################################################################
##
#F  MeatAxePermutationsOps.Order( MeatAxePermutations, <mtxperm> )
##
MeatAxePermutationsOps.Order := function( MeatAxePermutations, mtxperm )
    return OrderMeatAxePerm( mtxperm );
    end;

#############################################################################
##
#F  MeatAxePermOps.GapObject( <mtxperm> )
##
##  returns the {\GAP} permutation corresponding to the {\MeatAxe}
##  permutation <mtxperm>.
##
MeatAxePermOps.GapObject := function( mtxperm )

    local name,   # temporary file name
          result; # the {\GAP} permutation, result
  
    name:= TmpName();
    MeatAxe.EXEC( "zpr -G ", mtxperm.file, " ", name );
    Read( name );
    if not IsBound( MeatAxe.Perms ) then
      Error( "file '",mtxperm.file,"'does not contain MeatAxe permutations");
    fi;
#T zpr should be improved!
    if MeatAxe.Perms = [] then
      result:= ();
    else
      result:= MeatAxe.Perms[1];
    fi;
  
    EXEC( "rm ", name );
    Unbind( MeatAxe.Perms );
  
    return result;
    end;

#############################################################################
##
#F  MeatAxePermOps.Print . . . . . . . . . .   print a {\MeatAxe} permutation
##
MeatAxePermOps.Print := function( obj )
    Print( "MeatAxePerm( \"", obj.file, "\", ", obj.maxpoint, " )" );
    end;

#############################################################################
##
#F  MeatAxePermOps.\*( <p1>, <p2> ) . . .  product of {\MeatAxe} permutations
##
MeatAxePermOps.\* := function( p1, p2 )
    local name,     # file name of the result
          mtxperm;  # result

    if IsMeatAxePerm( p1 ) and IsMeatAxePerm( p2 ) then

      name:= MeatAxe.TmpName( MeatAxe.files );
      MeatAxe.EXEC( "zmu ", p1.file, " ", p2.file, " ", name );
      mtxperm:= MeatAxePerm( name, p1.maxpoint );

    fi;
    return mtxperm;
    end;

#############################################################################
##
#F  MeatAxePermOps.\^( <p1>, <p2> )
##
##  <mtxperm> ^ <mtxperm>  conjugation
##  <mtxperm> ^ <int>      powering
##
MeatAxePermOps.\^ := function( p1, p2 )

    local name, name2,
          mtxperm;        # result

    if IsMeatAxePerm( p1 ) then

      if IsMeatAxePerm( p2 ) then

        # conjugate permutations, return p1^-1 * p2 * p1
        mtxperm:= p1^-1 * p2 * p1;

      elif IsInt( p2 ) then

        # <perm> ^ <int>
        if p2 > 1 then

          name:= MeatAxe.TmpName( MeatAxe.files );
          MeatAxe.EXEC( "zsm pwr", p2, " ", p1.file, " ", name );
          mtxperm:= MeatAxePerm( name, p1.maxpoint );

          if IsBound( p1.abstract ) then
            mtxperm.abstract:= p1.abstract ^ p2;
          fi;

        elif p2 = 1 then

          mtxperm:= p1;

        elif p2 = 0 then

          if     IsBound( p1.structure )
             and IsBound( p1.structure.identity ) then
            mtxperm:= p1.structure.identity;
          else

            mtxperm:= MeatAxePerm( (), p1.maxpoint );
            if IsBound( p1.abstract ) then
              mtxperm.abstract:= p1.abstract ^ 0;
            fi;
  
          fi;

        else

          if not IsBound( p1.inverse ) then

            # write the inverse of 'p1' to 'name'
            name:= MeatAxe.TmpName( MeatAxe.files );
            MeatAxe.EXEC( "ziv ", p1.file, " ", name );
            p1.inverse:= MeatAxePerm( name, p1.maxpoint );
            p1.inverse.inverse:= p1;

            if IsBound( p1.abstract ) then
              p1.inverse.abstract:= p1.abstract ^ -1;
            fi;

          fi;

          if p2 = -1 then
            mtxperm:= p1.inverse;
          else
            name2:= MeatAxe.TmpName( MeatAxe.files );
            MeatAxe.EXEC( "zsm pwr", -p2, " ", p1.inverse.file,
                  " ", name2 );
            mtxperm:= MeatAxePerm( name2, p1.maxpoint );

            if IsBound( p1.inverse.abstract ) then
              mtxperm.abstract:= p1.inverse.abstract ^ -p2;
            fi;

          fi;

        fi;

      else
        mtxperm:= GapObject( p1 ) ^ GapObject( p2 );
      fi;

    else
      mtxperm:= GapObject( p1 ) ^ GapObject( p2 );
    fi;

    return mtxperm;
    end;

#############################################################################
##
#F  MeatAxePermOps.Display( <obj>, <dummyrec> )
##
MeatAxePermOps.Display := function( obj, dummyrec )
    MeatAxe.EXEC( "zpr -G ", obj.file );
    end;

#############################################################################
##
#F  MeatAxePermOps.\=
##
MeatAxePermOps.\= := MeatAxeMatOps.\=;

#############################################################################
##
#F  MeatAxePermOps.\<( <p1>, <p2> ) . . . ordering of {\MeatAxe} permutations
##
MeatAxePermOps.\< := function( p1, p2 )
    return p1.file < p2.file;
    end;

############################################################################
##
#F  MeatAxePerm( <perm>, <maxpoint> [,<file>] ) .  construct {\MeatAxe} perm
#F  MeatAxePerm( <file>, <maxpoint> )  . . . . . . .  notify {\MeatAxe} perm
##
MeatAxePerm := function( arg )

    local linelen,  # line length
          maxpoint, # largest moved point
          npl,      # no. of points per line
          nol,      # no. of lines
          func,     # function that prints the permutation
          i, j,     # loop over points
          perm,     # permutation, first argument
          file,     # file that will contain the {\MeatAxe} permutation
          tmp;      # temporary file with ASCII permutation

    if   Length( arg ) = 2 and IsString( arg[1] ) and IsInt( arg[2] ) then

      # Just notify a {\MeatAxe} permutation.
      file:= arg[1];
      maxpoint:= arg[2];

    elif Length( arg ) >= 2 and IsPerm( arg[1] ) and IsInt( arg[2] ) and
         (   Length( arg ) = 2 or
           ( Length( arg ) = 3 and IsString( arg[3] ) ) ) then

      # Set the print parameters.
      perm:= arg[1];
      maxpoint:= arg[2];
      linelen:= 80;
      npl:= Int( linelen / ( 1 + Length( String( maxpoint ) ) ) );
      nol:= Int( maxpoint / npl );
  
      # Define a function that prints the permutation.
      func:= function()
        Print( "12 1 ", maxpoint, " 1\n" );
        for i in [ 0 .. nol-1 ] do
          for j in [ 1 .. npl ] do
            Print( ( i * npl + j ) ^ perm, " " );
          od;
          Print( "\n" );
        od;
        for i in [ npl * nol + 1 .. maxpoint ] do
          Print( i ^ perm, " " );
        od;
        Print( "\n" );
        end;
  
      # Print the permutation to the intermediate file.
      tmp:= TmpName();
      PrintTo( tmp, func() );
  
      # Get the file name.
      if Length( arg ) = 3 then
        file:= arg[3];
        if file[1] <> '/' then
          file:= Concatenation( MeatAxe.direc, "/", file );
          AddSet( MeatAxe.files, file );
        fi;
      else
        file:= MeatAxe.TmpName( MeatAxe.files );
      fi;
  
      # Convert the permutation to internal format.
      MeatAxe.EXEC( "zcv ", tmp, " ", file );
  
      # Make clean.
      EXEC( "rm ", tmp );

    else

      Error( "usage: MeatAxePerm( <perm>, <maxpoint> [,<file>] )\n",
             " or MeatAxePerm( <file>, <maxpoint> )" );

    fi;

    # Return the {\MeatAxe} permutation.
    return rec(
                isPermutation := true,
                isMeatAxePerm := true,
                domain        := MeatAxePermutations,
                file          := file,
                maxpoint      := maxpoint,
                operations    := MeatAxePermOps
               );
    end;

