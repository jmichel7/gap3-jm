#######################################################################
##  SPECHT 3.1 - init.g : Specht initialization                      ##
##                                                                   ##
##     A GAP package for calculating the decomposition numbers of    ##
##     Hecke algebras of type A (over fields of characteristic       ##
##     zero). The functions provided are primarily combinatorial in  ##
##     nature. Many of the combinatorial tools from the (modular)    ##
##     representation theory of the symmetric groups appear in the   ##
##     package.                                                      ##
##                                                                   ##
##     These programs, and the enclosed libraries, are distributed   ##
##     under the usual licensing agreements and conditions of GAP.   ##
##                                                                   ##
##     Andrew Mathas          mathas@maths.usyd.edu.au               ##
##     University of Sydney   Sydney, 1997                           ##
##                                                                   ##
##        Version 3.1:  ???????   ????  University of Sydney         ##
##        Version 2.4:  October   1997  University of Sydney         ##
##        Version 2.3:  December  1996  Imperial College             ##
##        Version 2.2:  June      1996  Imperial College             ##
##        Version 2.1:  April     1996  Imperial College             ##
##        Version 2.0:  March     1996  Imperial College             ##
##        Version 1.0:  September 1995  Imperial College             ##
##                                                                   ##
#######################################################################

## Change log

## 2.4:
##   - replaced SPECHTHOME, SPECHTVERSION etc with SPECHT.Home,
##     SPECHT.Version etc.
##
## 2.2: June 1996: various changes requested by the referee.
##   - moved handling functions into induce.g and deleted SpechtParts.
##   - changed various function names
##
## 2.1: April 1996:
##  - Changed LabelPartition so that it no longer generates all
##    partitions of n when the given a 'new' partition.
##  - Fixed typo in AUTO command for IsSpecht in init.g.
##
## 2.0: March 1996: 
##   - broke specht.g up into smaller files and added AUTO commands.
##     The component files of Specht are now:
##       init.g:     general handling functions (this file).
##       specht.g:   functions which are heavily dependent upon
##                   Specht() and Schur() records, or decomposition
##                   matrices. In particular this file contains the
##                   functions Specht, Schur, and all decomposition
##                   matrix routines.
##       symmcomb.g: functions for combinatorics on Young diagrams and
##                   partitions (combinatorics of ths symmetric groups).
##       tableaux.g: operations on (semistandard) tableaux.
##       tex.g:      TeX() functions (this really should be moved into
##                   the GAP library).
##   - added Specialized() function.
##   - added SpechtParts() giving streamlined access to partitions and
##     associated strings.
##   - default decomposition matrix libraries reduced to those that 
##     Specht is unable to calculate using LLT. 
##   - added tex.g library to allow easy TeXing of (crystallized)
##     decomposition matrices with polynomial entries.
##   - put tableaux functions into tableaux.g
##   - added opening message
##
## 1.0: December 1995: initial release.

###########################################################################

PrintPkgInit(rec(
  name:="specht",
  version:="3.1 - unofficial (development)",
  copyright:=[
"Decomposition numbers of Hecke algebras of type A and q-Schur algebras.",
#"    This release is unofficial and contains *undocumented* functions",
"(C) Andrew Mathas  mathas@maths.usyd.edu.au     Sydney"]));

###########################################################################
SPECHT:=OperationsRecord("Specht data");
SPECHT.Home:=Concatenation(LOADED_PACKAGES.specht,"gap/");
SPECHT.Library:=Concatenation(LOADED_PACKAGES.specht,"lib/");
SPECHT.Version:="version 3.1, (C) Andrew Mathas 2004.";
SPECHT.Email:="mathas@maths.usyd.edu.au";
SPECHT.DecompositionNumberOK:=true;
SpechtInfo:=Ignore;

# Caching of basis elements of the Fock space..............................
# We cache the cannonical bases elements Pq(mu) in SPECHT.CACHE.Pq. For a  
# given e and a partition mu of n the cannonical bases element Pq(mu) is   
# cached as                                                                
#    SPECHT.Cache[H.e][n].Pq.b[Position(SPECHT.Cache[H.e][n].Pq.mu,mu)]    
# and smilarly for the Aq() basis. The cache is accessed via the functions 
# SPECHT.PutCached(H,item,b) and SPECHT.GetCached(H.item,mu) where, in both
# cases item is either "Pq" and "Aq" and the other arguments have the      
# obvious meaning.                                                         
SPECHT.Cache:=[];
SPECHT.PutCached:=function(H,item,b) local mu, muStr, n, Mu;
  mu:=b.parts[Length(b.parts)];
  n:=Sum(mu);
 if IsBound(SPECHT.Cache[H.e][n]) then Mu:=Position(SPECHT.Cache[H.e][n].(item).mu,mu);
  else 
    SPECHT.Cache[H.e][n]:=rec(Aq:=rec(b:=[],mu:=()),Pq:=rec(b:=[],mu:=[]));
    Mu:=false;
  fi;
  if Mu=false then
    Add(SPECHT.Cache[H.e][n].(item).mu,mu);
    Add(SPECHT.Cache[H.e][n].(item).b,b);
  fi;
end;
SPECHT.GetCached:=function(H,item,mu) local n, muStr, Mu, Pq;
  n:=Sum(mu);
  if IsBound(SPECHT.Cache[H.e][n]) then
    Mu:=Position(SPECHT.Cache[H.e][n].(item).mu,mu);
    if Mu<>false then return ShallowCopy(SPECHT.Cache[H.e][n].(item).b[Mu]); fi;
  fi;
  return false;
end;

#######################################################################

ReadSpecht:=function(name)
   if not ReadPath(SPECHT.Home, name, ".g", "ReadSpecht") then
      Error("SPECHT library file '", name, "' must exist and be readable");
  fi;
end;

###########################################################################

ReadSpecht( "auto" );  
AUTO(ReadSpecht("induce"), SpechtPrintFn );  

## The following directory is searched by ReadDecompositionMatrix()
## when it is looking for decomposition matrices. By default, it points
## to the current directory (if set, the current directory is not
## searched).
if not IsBound(SpechtDirectory) then SpechtDirectory:=""; fi;

## This variable is what is used in the decomposition matrices files saved
## by SaveDecompositionMatrix() (and also the variable which contains them
## when they are read back in).
A_Specht_Decomposition_Matrix:=false;
