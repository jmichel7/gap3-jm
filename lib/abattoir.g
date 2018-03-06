#############################################################################
##
#A  abattoir.g                  GAP library                  Martin Schoenert
##
#A  @(#)$Id: abattoir.g,v 1.3 1997/10/12 10:54:02 gap Exp $
##
#Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file is only there to catch some loose ends.
##
#H  $Log: abattoir.g,v $
#H  JM 2015/02/25 introduced OpsOps.String from chevie
#H  Revision 1.3  1997/10/12 10:54:02  gap
#H  DOS/MAC
#H
#H  Revision 1.2  1997/03/27 12:21:37  gap
#H  Small change to ReadPkg. Now tries bare filename first and then .g file
#H  if that fails. Enables Specht to use ReadPkg.
#H
#H  	SL
#H
#H  Revision 1.1.1.1  1996/12/11 12:37:42  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 3.16.1.2  1994/11/17  16:41:38  sam
#H  added functions 'True' and 'False'
#H
#H  Revision 3.16.1.1  1995/06/26  10:44:28  mschoene
#H  added OpsOps.( "=" )
#H
#H  Revision 3.16  1994/06/23  11:53:51  sam
#H  added 'IsOddInt', 'IsEvenInt' (Frank needs ...)
#H
#H  Revision 3.15  1994/04/21  11:56:12  sam
#H  added 'EXEC', and operations records share their operations records now
#H
#H  Revision 3.14  1994/02/02  10:26:45  sam
#H  moved some functions to 'dispatch.g',
#H  added 'OperationsRecord'
#H
#H  Revision 3.13  1993/07/30  07:27:58  fceller
#H  added 'ChiefSeries'
#H
#H  Revision 3.12  1993/05/25  08:34:09  fceller
#H  moved 'IsIdentical' into the kernel
#H
#H  Revision 3.11  1993/02/10  18:00:21  martin
#H  added nondestructive blist functions
#H
#H  Revision 3.10  1993/02/09  14:12:55  fceller
#H  changed 'PKGNAME'
#H
#H  Revision 3.9  1993/02/04  11:22:41  martin
#H  added 'SubString', and 'ConcatenationString'
#H
#H  Revision 3.8  1992/12/07  07:42:58  fceller
#H  added 'ExecPkg'
#H
#H  Revision 3.7  1992/12/04  13:41:12  fceller
#H  changed 'ReadPkg' to allow variable number of args
#H
#H  Revision 3.6  1992/12/03  12:23:01  fceller
#H  renamed 'IsEquivalent' to 'IsIdentical'
#H
#H  Revision 3.5  1992/11/19  14:34:43  fceller
#H  added package functions
#H
#H  Revision 3.4  1992/08/19  10:00:24  fceller
#H  added 'Save'
#H
#H  Revision 3.3  1992/08/18  09:30:20  fceller
#H  added Juergen's lattice print functions
#H
#H  Revision 3.2  1992/08/14  16:09:31  fceller
#H  added a primitive 'IsEquivalent' for lists
#H
#H  Revision 3.1  1992/04/07  20:07:51  martin
#H  initial revision under RCS
##

#############################################################################
##
#F  SubString(<args>) . . . . . . . . . . . . . . . . . substring of a string
##
##  'SubString' is obsolete since strings are lists.
##
SubString := function ( arg )
    local    string, from, to;
    string := arg[1];
    if not IsString( string )  then
        Error("usage: SubString( <string>, <from> [, <to>] )");
    fi;
    from := arg[2];
    if not IsInt( from )  then
        Error("usage: SubString( <string>, <from> [, <to>] )");
    fi;
    if Length(arg) = 3  then
        to := arg[3];
    else
        to := Length(string);
    fi;
    if to < from  then return "";  fi;
    if from <= 0  then  from := 1;  fi;
    if Length(string) < from   then from := Length(string);  fi;
    if to   <= 0  then  to   := 1;  fi;
    if Length(string) < to     then to   := Length(string);  fi;
    return string{[from..to]};
end;

#############################################################################
##
#F  ConcatenationString(<string>,..)  . . . . . . .  concatenation of strings
##
##  'ConcatenationString' is obsolete since strings are lists.
##
ConcatenationString := function ( arg ) local   res,  str;
  res := "";
  for str  in arg  do Append( res, str ); od;
  IsString( res );
  return res;
end;

#############################################################################
##
#F  Edit(<name>)  . . . . . . . . . . . . . . . . . . . . . . . . edit a file
##
if not IsBound( EDITOR )  then EDITOR := "vi";  fi;

Edit := function ( name )
    Exec( ConcatenationString( EDITOR, " ", name ) );
    Read( name );
end;

#############################################################################
##
#F  ProductPol( <f>, <g> )  . . . . . . . . . . .  product of two polynomials
##
ProductPol := function ( f, g )
    local  prod,  q,  m,  n,  i,  k;
    m := Length(f);  while 1 < m  and f[m] = 0  do m := m-1;  od;
    n := Length(g);  while 1 < n  and g[n] = 0  do n := n-1;  od;
    prod := [];
    for i  in [ 2 .. m+n ]  do
        q := 0;
        for k  in [ Maximum(1,i-n) .. Minimum(m,i-1) ]  do
            q := q + f[k] * g[i-k];
        od;
        prod[i-1] := q;
    od;
    return prod;
end;


#############################################################################
##
#F  ValuePol( <f>, <x> )  . . . . . . . . . . . evaluate a polynom at a point
##
ValuePol := function ( f, x )
    local  value, i, id;
    id := x ^ 0;
    value := 0 * id;
    i := Length(f);
    while 0 < i  do
        value := x * value + id * f[i];
        i := i-1;
    od;
    return value;
end;


#############################################################################
##
#F  MergedRecord(<rec1>, .., <recn>)   . . . . . merge the fields of records
##
MergedRecord := function ( arg )
    local   res,        # merged record, result
            record,     # one of the arguments
            name;       # name of one component of <record>
    res := rec();
    for record  in  arg do
        for name  in RecFields( record )  do
            if IsBound( res.(name) )  then
                Unbind( res.(name) );
            else
                res.(name) := record.(name);
            fi;
        od;
    od;
    return res;
end;

#############################################################################
##
#F  UnionBlist( <blist1>, .. , <blistn> ) . . . . . . . . . . union of blists
##
UnionBlist := function ( arg ) local  U, i;
    if Length( arg ) = 1  then arg := arg[1]; fi;
    U := Copy( arg[1] );
    for i  in [2..Length(arg)]  do
        UniteBlist( U, arg[i] );
    od;
    return U;
end;

#############################################################################
##
#F  IntersectionBlist( <blist1>, <blist2> ) . . . . .  intersection of blists
##
IntersectionBlist := function ( arg )
    local  I, i;
    if Length( arg ) = 1  then
        arg := arg[1];
    fi;
    I := Copy( arg[1] );
    for i  in [2..Length(arg)]  do
        IntersectBlist( I, arg[i] );
    od;
    return I;
end;


#############################################################################
##
#F  DifferenceBlist( <blist1>, <blist2> ) . . . . . . .  difference of blists
##
DifferenceBlist := function ( blist1, blist2 )
    local  D;
    D := Copy( blist1 );
    SubtractBlist( D, blist2 );
    return D;
end;


#############################################################################
##
#F  SetPrintLevel( <L>, <lev> ) . . . . . . . . . . .  set print level of <L>
##
SetPrintLevel := function( L, lev )
   L.operations.SetPrintLevel( L, lev );
end;


#############################################################################
##
#F  Save( <file>, <obj>, <name> ) . . . . . . . . . save some strange objects
##
Save := function( F, G, N )
    if not IsRec(G) or not IsBound(G.operations.Save)  then
        Error( "sorry, I do not know how to save <G>" );
    fi;
    G.operations.Save( F, G, N );
end;


#############################################################################
##
#V  PKGNAME . . . . . . . . . . . . . . . . . . . location of share libraries
##
SetPkgname := function( path )
    local   i,  l,  p;

    # copy old path
    path := Copy(path);

    # append final ';'
    if path[Length(path)] <> ';'  then
        Add( path, ';' );
    fi;

    # replace "lib/;" by "pkg/;"
    for i  in [ 1 .. Length(path)-4 ]  do
    	if path{[i..i+4]} = "lib/;"  then
    	    path{[i..i+4]} := "pkg/;";
    	elif path{[i..i+4]} = "lib\\;"  then  # DOS
    	    path{[i..i+4]} := "pkg\\;";
    	elif path{[i..i+4]} = "lib:;"  then   # MacOS
    	    path{[i..i+4]} := "pkg:;";
    	fi;
    od;

    # now split paths
    p := [];
    l := 1;
    for i  in [ 1 .. Length(path) ]  do
    	if path[i] = ';'  then
    	    Add( p, path{[l..i-1]} );
    	    IsString( p[Length(p)] );
    	    l := i+1;
    	fi;
    od;

    # and return
    return p;

end;

PKGNAME := SetPkgname(LIBNAME);


#############################################################################
##
#F  ReadPkg( <lib>, <name> )  . . . . . . . . . .   read a share library file
##
LOADED_PACKAGES := rec();

ReadPkg := function( arg )
    local   ind,  fln,  i, fln1;

    # store old indent value, add two spaces
    ind := ReadIndent;
    ReadIndent := ConcatenationString( ReadIndent, "  " );

    # construct complete path
    fln := Copy( LOADED_PACKAGES.(arg[1]) );
    for i  in [ 2 .. Length(arg)-1 ]  do
    	Append( fln, arg[i] );
    	Add( fln, '/' );
    od;
    Append( fln, arg[Length(arg)] );
    IsString(fln);
            
    # read in file -- try bare name first. If that fails, try with .g
    if not READ(fln)  then
        Append( fln, ".g" );
        InfoRead1( "#I", ReadIndent, "ReadPkg( \"", fln, "\" )\n" );
                
        # read in file
        if not READ(fln)  then
            Error("share library file \"",fln,"\" must exist and be readable");
         fi;
    else
        InfoRead1( "#I", ReadIndent, "ReadPkg( \"", fln, "\" )\n" );
    fi;

    # restore old indentation
    ReadIndent := ind;

end;


#############################################################################
##
#F  ExecPkg( <lib>, <cmd>, <ags>, <dir> ) . . . . .  execute a package binary
##
##  Change to the directory <dir> and execute <cmd> with arguments <ags>.
##
ExecPkg := function( lib, cmd, ags, dir )
    local   del,  new,  i,  sub;

    # prefix <cmd> with path
    new := Copy( LOADED_PACKAGES.(lib) );
    Append( new, cmd );

    # construct the command line
    cmd := ConcatenationString( "cd ", dir, "; ", new, " ", ags );
    InfoRead1( "#I  ExecPkg: executing ", cmd, "\n" );
    Exec(cmd);

end;


#############################################################################
##
#F  LoadPackage( <name> ) . . . . . . . . . . .  load a share library package
##
LoadPackage := function( name )
    local   path,  init,  ind;

    # store old indent value, add two spaces
    ind := ReadIndent;
    ReadIndent := ConcatenationString( ReadIndent, "  " );

    # find the share library <name>
    for path  in PKGNAME  do

    	# check next <path>
    	init := Copy(path);
    	Append( init, name );
    	Append( init, "/" );
    	IsString(init);
    	LOADED_PACKAGES.(name) := Copy(init);
    	Append( init, "init.g" );
    	IsString(init);

    	# give read info
    	InfoRead1( "#I  LoadPackage tries \"", init, "\"\n" );

    	# try to read the init file
    	if READ(init)  then
    	    ReadIndent := ind;
    	    return init;
    	fi;
    od;

    # signal an error
    Unbind( LOADED_PACKAGES.(name) );
    ReadIndent := ind;
    Error( "share library \"", name, "\" is not installed" );

end;


#############################################################################
##
#F  RequirePackage( <name> )  . . . . . . . . . .  make sure <name> is loaded
##
RequirePackage := function( arg )local name;
  for name in arg do # check if <name> is already loaded
    if not IsBound(LOADED_PACKAGES.(name)) then LoadPackage( name );fi;
  od;
end;

PkgDate:=function(arg)local months,res;
  res:=String(arg[1]);
  if Length(arg)>=2 then Append(res," ");Append(res,NameMonth[arg[2]]);fi;
  if Length(arg)>=3 then Append(res," ");Append(res,String(arg[3]));fi;
  return res;
end;

#############################################################################
##
#F  PrintPkgInit( <Name>, <Version> )  . . Normalized greeting message of Pkg
##
PrintPkgInit:=function(opt)local fill,ver,i;
  if QUIET or not BANNER then return;fi;
  fill:=function(s,n)local res,i;res:=" ";Append(res,s);Append(res," ");
    for i in [1..n-Length(res)] do Add(res,'-');od;
    return res;
  end;
  ver:="version";
  if IsBound(opt.version) then Append(ver," ");Append(ver,opt.version);fi;
  if IsBound(opt.date) then Append(ver," of ");
     Append(ver,ApplyFunc(PkgDate,opt.date));
  fi;
  Print("--- Loading package",fill(opt.name,14),fill(ver,44),"\n");
  if IsBound(opt.copyright) then 
    if IsString(opt.copyright) then Print("    ",opt.copyright,"\n");
    else for i in [1..Length(opt.copyright)] do
      Print("    ",opt.copyright[i],"\n");od;
    fi;
  fi;
end;

#############################################################################
##
#F  IsOperationsRecord( <obj> ) . . . . . . . . . . .  category test function
##
IsOperationsRecord := function( obj )
    return     IsRec( obj )
           and IsBound( obj.name )
           and IsBound( obj.operations )
           and IsRec( obj.operations )
           and IsBound( obj.operations.name )
           and obj.operations.name = "OpsOps";
end;

##############################################################################
##
#V  OpsOps
##
## JM 2015/02/25 introduced String from chevie
OpsOps := rec( name:="OpsOps",
               String:=function(obj)return obj.name;end,
	       Print:=function(obj)Print(String(obj));end);
               
OpsOps.( "=" ) := function( oprec1, oprec2 )
    if IsOperationsRecord( oprec1 )  then
        if IsOperationsRecord( oprec2 )  then
            return oprec1.name = oprec2.name;
        else
            return false;
        fi;
    elif not IsOperationsRecord( oprec2 )  then
        Error( "panic, neither argument is an operations record" );
    fi;
end;

OpsOps.operations := OpsOps;

##############################################################################
##
#F  OperationsRecord( <name> )
#F  OperationsRecord( <name>, <parent> )
##
OperationsRecord := function( arg ) local oprec;

    if Length( arg ) = 1 then oprec:= rec();
    else oprec:= Copy( arg[2] );
    fi;
    oprec.name:= arg[1];
    oprec.operations:= OpsOps;
    return oprec;
end;  

##############################################################################
##
#F  EXEC( <str1>, <str2>, ... )
##
##  This should become the standard of 'Exec'
##
EXEC := function( arg )
    Exec( Concatenation( List( arg, String ) ) );
end;

IsOddInt := n->n mod 2 = 1;

IsEvenInt := n->n mod 2 = 0;

False := function( arg )
    return false;
end;

True := function( arg )
    return true;
end;
