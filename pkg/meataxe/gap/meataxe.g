
#############################################################################
##
#A  meataxe.g                   GAP Library                     Thomas Breuer
##
#A  @(#)$Id: meataxe.g,v 1.1.1.1 1996/12/11 12:40:13 werner Exp $
##
#Y  Copyright 1991-1993,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains the interface between {\GAP} and the {\MeatAxe}.
##
##  Functions that deal with {\MeatAxe} matrices, {\MeatAxe} permutations,
##  {\MeatAxe} matrix algebras, and {\MeatAxe} modules can be found in the
##  files 'mamatrix.g', 'mapermut.g', 'mamatalg.g', 'mamodule.g',
##  respectively.
##
#H  $Log: meataxe.g,v $
#H  Revision 1.1.1.1  1996/12/11 12:40:13  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 3.2  1996/11/13 16:05:21  sam
#H  use 'ExecPkg' and 'LOADED_PACKAGES'
#H
#H  Revision 3.1  1994/05/19 14:05:19  sam
#H  changes for GAP-3.4
#H
#H  Revision 3.0  1994/03/28  15:31:13  sam
#H  initial revision under RCS
#H
##

############################################################################
##
#V  MeatAxe
##
##  is a record containing relevant information about the usage of MeatAxe
##  under {\GAP}.
##  Currently there are the following components.
##
##  'PATH':     the path name of the {\MeatAxe} directory
##
##  'fields':   a list where position <i> is bound if and only if
##              the field of order <i> has already been constructed by
##              the 'maketab' command; in this case it contains the
##              name of the 'pxxx.zzz' file
##  'files':    a list of all file names that were constructed by calls
##              to {\MeatAxe} (for allowing to make clean \ldots)
##  'nrnames':  the number of names already produced
##  'dirs':     a list of all directory names that were constructed by calls
##              to {\MeatAxe} (for allowing to make clean \ldots)
##  'gennames': list of strings that are used as generator names
##              in 'abstract' components of {\MeatAxe} matrices
##  'alpha':    alphabet ober which 'gennames' entries are formed
##  'Unbind':   function to delete superfluous files (see below)
##  'Maketab':  function that computes the field information file,
##  'TmpName':  function to create new file names in the {\MeatAxe}
##              directory
##  'direc':    the {\MeatAxe} directory
##  'NUMBER':   ...
##
##  Furthermore there are some components which are intermediately bound
##  when {\MeatAxe} output files are read.
##
MeatAxe := rec(
                PATH     := Concatenation(LOADED_PACKAGES.meataxe,"bin/"),
                lib     := Concatenation(LOADED_PACKAGES.meataxe,"lib/"),
                fields   := [],
                files    := [],
                nrnames  := 0,
                dirs     := [],
                gennames := [],
                alpha    := "abcdefghijklmnopqrstuvwxyz",
                direc    := TmpName(),
                istmpdir := true
               );

EXEC( "mkdir ", MeatAxe.direc );
Print("#I  All MeatAxe files will be placed in the directory  '",MeatAxe.direc,"'\n",
      "#I  Use 'MeatAxe.SetDirectory( <path> )' if you want to change.\n" );

#############################################################################
##
#F  MeatAxe.Unbind( <obj1>, <obj2>, ..., <objn> )
#F  MeatAxe.Unbind( <listofobjects> )
##
##  Called without arguments, this removes all files and directories
##  constructed by calls of 'MeatAxeMat' and 'Group', provided they are still
##  notified in 'MeatAxe.files', 'MeatAxe.dirs' and 'MeatAxe.fields'.
##
##  Otherwise all those {\MeatAxe} matrices and field information files in
##  'MeatAxe.files', 'MeatAxe.dirs' and 'MeatAxe.fields' are removed that are
##  specified in the argument list.
##
MeatAxe.Unbind := function( arg )
    local i, list, obj, pos;

    if Length( arg ) = 0 then
      for i in MeatAxe.files  do EXEC( "rm "   , i ); od;
      for i in MeatAxe.dirs   do EXEC( "rm -r ", i ); od;
      for i in MeatAxe.fields do EXEC( "rm "   , i ); od;
      EXEC( "rm -r ", MeatAxe.direc );
      MeatAxe.files   := [];
      MeatAxe.dirs    := [];
      MeatAxe.fields  := [];
    else

      if Length( arg ) = 1 and IsList( arg[1] ) then
        list:= arg[1];
      else
        list:= arg;
      fi;

      for i in list do
        if IsRec( i ) and IsBound( i.file ) then
          obj:= i.file;
        else
          obj:= i;
        fi;
        if   obj in MeatAxe.files then
          EXEC( "rm ", obj );
          RemoveSet( MeatAxe.files, obj );
        elif obj in MeatAxe.dirs then
          EXEC( "rm -r ", obj );
          RemoveSet( MeatAxe.dirs, obj );
        fi;
        pos:= Position( MeatAxe.fields, obj );
        if pos <> false then
          EXEC( "rm ", MeatAxe.fields[ pos ] );
          Unbind( MeatAxe.fields[ pos ] );
        fi;
      od;
    fi;
    end;

#############################################################################
##
#F  MeatAxe.Maketab( <n> )
##
MeatAxe.Maketab := function( n )

    InfoMeatAxe( "#I  calling 'maketab' for field of size ", n, "\n" );

    if   n <= 256 then
      MeatAxe.EXEC( "maketab -Q ", String( n ) );
      MeatAxe.fields[n]:=
          Concatenation( MeatAxe.direc, "/p",
                         List( [ 1 .. 2-LogInt( n, 10 ) ],
                               x -> '0' ),
                         String( n ),
                         ".zzz" );
    elif n <= 65536 then
      MeatAxe.EXEC( "bigmktab -Q ", n );
      MeatAxe.fields[n]:=
          Concatenation( MeatAxe.direc, "/p",
                         List( [ 1 .. 4-LogInt( n ) ],
                               x -> '0' ),
                         String( n ),
                         ".zzz" );
    else
      Error( "sorry, this field is too large" );
    fi;

    # compute conversion information and store it
    FFList( GF(n) );
    end;

#############################################################################
##
#F  MeatAxe.TmpName( arg )
##
##  filename in directory 'MeatAxe.direc'.
##  This should be used *only* for those files that are not removed in the
##  same function call that causes the creation, because the variable
##  'MeatAxe' keeps track of the files in order to be able to remove them
##  in the end.
##
MeatAxe.TmpName := function( arg )
    local list,
          name;
    MeatAxe.nrnames:= MeatAxe.nrnames + 1;
    name:= Concatenation( MeatAxe.direc, "/",
                  WordAlp( MeatAxe.alpha, MeatAxe.nrnames ) );
    if Length( arg ) > 0 then
      list:= arg[1];
      AddSet( list, name );
      if list = MeatAxe.dirs then
        EXEC( "mkdir ", name );
      fi;
    fi;
    return name;
    end;

#############################################################################
##
#F  MeatAxe.SetDirectory( <path> )
##
MeatAxe.SetDirectory := function( path )
    if not IsString( path ) then
      Error( "<path> must be a string that denotes the directory name" );
    elif MeatAxe.istmpdir then
      MeatAxe.istmpdir:= false;
      if Length( MeatAxe.files ) = 0 then
        if BANNER then
          Print( "#I  Removing empty directory '", MeatAxe.direc, "'\n" );
        fi;
        EXEC( "rm -r ", MeatAxe.direc );
      elif BANNER then
        Print( "#I  There are still files in the old directory.\n" );
      fi;
    fi;
    MeatAxe.direc:= path;
    end;

#############################################################################
##
#F  MeatAxe.EXEC( arg )
##
MeatAxe.EXEC := function( arg )local cmd;
 #  ExecPkg( "meataxe", Concatenation( List( arg, String ) ),
 #           "", MeatAxe.direc );
    cmd:=Concatenation("cd ",MeatAxe.direc,"; export MTXLIB=", MeatAxe.lib,";", 
      MeatAxe.PATH, Concatenation( List( arg, String ) ) );
    InfoRead1( "#I  MeatAxe: executing ",cmd,"\n");
    Exec(cmd);
end;

#############################################################################
##
#F  MeatAxe.NUMBER( <string> )
##
##  takes a string whose initial part consists of digits, and returns the
##  corresponding number
##
MeatAxe.NUMBER := function( string )
    local i, digits, number;
    digits:= "0123456789";
    i:= 1;
    number:= 0;
    while i <= Length( string ) and string[i] in digits do
      number:= 10 * number + Position( digits, string[i] ) - 1;
      i:= i+1;
    od;
    return number;
    end;

#############################################################################
##
#F  GapObject( <obj> ) . . . . . . . . . {\GAP} object corresponding to <obj>
##
GapObject := function( obj )
    if     IsRec( obj )
       and IsBound( obj.gapObject ) then
      return obj.gapObject;
    elif   IsRec( obj )
       and IsBound( obj.operations )
       and IsBound( obj.operations.GapObject ) then
      obj.gapObject:= obj.operations.GapObject( obj );
      if IsRec( obj.gapObject ) then
        obj.gapObject.meatAxeObject:= obj.gapObject;
      fi;
      return obj.gapObject;
    else
      return obj;
    fi;
    end;

#############################################################################
##
#F  SplittingField( <obj> )
##
SplittingField := function( obj )
    if not IsBound( obj.splittingField ) then
      obj.splittingField:= obj.operations.SplittingField( obj );
    fi;
    return obj.splittingField;
    end;

#############################################################################
##
#F  MeatAxe.InverseMatMod( <intmat>, <p> )
##
##  does the same as 'InverseMatMod', but uses the {\MeatAxe}
##
MeatAxe.InverseMatMod := function( intmat, p )

    local tmp, tmp2;  # temporary file names

    # Make the field information file if necessary.
    if not IsBound( MeatAxe.fields[p] ) then
      MeatAxe.Maketab( p );
    fi;
    
    # Print the matrix to a temporary file, invert it, and decode the result.
    tmp:= TmpName();
    PrintTo( tmp, PrintMeatAxeInput( intmat, p, Dimensions( intmat ) ) );
    tmp2:= TmpName();
    MeatAxe.EXEC( "zcv "   , tmp , " "  , tmp2 );
    MeatAxe.EXEC( "ziv "   , tmp2, " "  , tmp  );
    MeatAxe.EXEC( "zpr -G ", tmp , " > ", tmp2 );

    # Get the inverse.
    Read( tmp2 );
    intmat:= MeatAxe.Matrix;

    # Make clean.
    EXEC( "rm ", tmp, " ", tmp2 );

    # Return the result.
    return intmat ;
    end;

