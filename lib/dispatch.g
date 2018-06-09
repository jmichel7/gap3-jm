#############################################################################
##
#A  dispatch.g                  GAP library                     Thomas Breuer
##
#A  @(#)$Id: dispatch.g,v 1.1.1.1 1996/12/11 12:37:51 werner Exp $
##
#Y  Copyright 1990-2900,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains dispatchers of the following two types.
##
##  The first expect exactly one argument, a record, and return
##  the corresponding record component.
##  If necessary this component is computed using the
##  corresponding function in the operations record of the argument, and
##  stored in the argument.  They perform no additional checks.
##
##  The second expect either a matrix or a record.
##  In case of a record the same holds as for the dispatchers of first type,
##  in case of a matrix the function with name extension 'Mat' is called.
##  
#H  JM 2015/02/25 ported from chevie improvement to Eigenvalues
#H
#H  JM 2016/03/15 ported from chevie Dispatcher, AttributeDispatcher
#H
#H  Revision 1.1.1.1  1996/12/11 12:37:51  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 3.5.1.2  1995/12/05  10:31:27  sam
#H  'SylowSystem' does no longer store its result in the group record
#H
#H  Revision 3.5.1.1  1994/10/26  16:18:25  ahulpke
#H  Added 'AutomorphismGroup'
#H
#H  Revision 3.5  1994/06/18  13:02:48  sam
#H  added bad hack in 'Base'
#H
#H  Revision 3.4  1994/06/17  11:12:14  sam
#H  fixed 'Centralizer'
#H
#H  Revision 3.3  1994/06/15  12:49:59  sam
#H  fixed 'Components'
#H
#H  Revision 3.0  1993/12/09  11:05:31  sam
#H  Initial revision under RCS
#H
##

Dispatcher:=function(name)
  return function(arg)local obj; obj:=arg[1];
    if IsRec(obj) and IsBound(obj.operations) and IsBound(obj.operations.(name))
    then return ApplyFunc(obj.operations.(name),arg);
    else Error(obj," has no method for ",name);
    fi;
  end;
end;
  
AttributeDispatcher:=function(arg)local name,fieldname,upper,lower;
  name:=arg[1];
  if Length(arg)=2 then fieldname:=arg[2];
  else
    fieldname:=ShallowCopy(name);
    upper:="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    lower:="abcdefghijklmnopqrstuvwxyz";
    fieldname[1]:=lower[Position(upper,fieldname[1])];
  fi;
  return function(arg)local obj;obj:=arg[1];
    if not IsRec(obj) then Error(obj," has no method for ",name);fi;
    if not IsBound(obj.(fieldname)) then
      if IsBound(obj.operations) and IsBound(obj.operations.(name)) then
	obj.(fieldname):=ApplyFunc(obj.operations.(name),arg);
      else Error(obj," has no method for ",name);
      fi;
    fi;
    return obj.(fieldname);
  end;
end;

AddGenerator:=Dispatcher("AddGenerator");
AddGenerators:=Dispatcher("AddGenerators");
Coefficient:=Dispatcher("Coefficient");
Degree:=Dispatcher("Degree");
Dual:=Dispatcher("Dual");
Expand:=Dispatcher("Expand");
LeadingCoefficient:=Dispatcher("LeadingCoefficient");
Representation:=Dispatcher("Representation");
SylowSystem:=Dispatcher("SylowSystem");

AbelianInvariants:=AttributeDispatcher("AbelianInvariants");
AbsoluteIrreducibilityTest:=AttributeDispatcher("AbsoluteIrreducibilityTest");
AsRing:=AttributeDispatcher("AsRing");
AsVectorSpace:=AttributeDispatcher("AsVectorSpace");
Automorphisms:=AttributeDispatcher("Automorphisms");
CanonicalBasis:=AttributeDispatcher("CanonicalBasis");
CanonicalRepresentative:=AttributeDispatcher("CanonicalRepresentative");
Centre:=AttributeDispatcher("Centre");
CharacterDegrees:=AttributeDispatcher("CharacterDegrees");
ChiefSeries:=AttributeDispatcher("ChiefSeries");
ClassInvariants:=AttributeDispatcher("ClassInvariants");
CommutatorFactorGroup:=AttributeDispatcher("CommutatorFactorGroup");
CompositionFactors:=AttributeDispatcher("CompositionFactors");
CompositionSeries:=AttributeDispatcher("CompositionSeries");
ConjugacyClasses:=AttributeDispatcher("ConjugacyClasses");
ConjugacyClassesPerfectSubgroups:=
  AttributeDispatcher("ConjugacyClassesPerfectSubgroups");
ConjugacyClassesSubgroups:=AttributeDispatcher("ConjugacyClassesSubgroups");
Constituents:=AttributeDispatcher("Constituents");
DerivedSeries:=AttributeDispatcher("DerivedSeries");
DerivedSubgroup:=AttributeDispatcher("DerivedSubgroup");
Dimension:=AttributeDispatcher("Dimension");
Discriminant:=AttributeDispatcher("Discriminant");
Exponent:=AttributeDispatcher("Exponent");
FittingSubgroup:=AttributeDispatcher("FittingSubgroup");
FpGroup:=AttributeDispatcher("FpGroup");
FrattiniSubgroup:=AttributeDispatcher("FrattiniSubgroup");
GaloisGroup:=AttributeDispatcher("GaloisGroup");
Generators:=AttributeDispatcher("Generators");
Idempotents:=AttributeDispatcher("Idempotents");
Identity:=AttributeDispatcher("Identity");
IdentityMapping:=AttributeDispatcher("IdentityMapping");
InvariantSubspace:=AttributeDispatcher("InvariantSubspace");
InverseMapping:=AttributeDispatcher("InverseMapping");
IrreducibilityTest:=AttributeDispatcher("IrreducibilityTest");
IsAbelian:=AttributeDispatcher("IsAbelian");
IsAssociative:=AttributeDispatcher("IsAssociative");
IsAutomorphism:=AttributeDispatcher("IsAutomorphism");
IsBijection:=AttributeDispatcher("IsBijection");
IsBijective:=AttributeDispatcher("IsBijective");
IsCyclic:=AttributeDispatcher("IsCyclic");
IsElementaryAbelian:=AttributeDispatcher("IsElementaryAbelian");
IsEndomorphism:=AttributeDispatcher("IsEndomorphism");
IsEpimorphism:=AttributeDispatcher("IsEpimorphism");
IsFaithful:=AttributeDispatcher("IsFaithful");
IsHomomorphism:=AttributeDispatcher("IsHomomorphism");
IsInjective:=AttributeDispatcher("IsInjective");
IsIsomorphism:=AttributeDispatcher("IsIsomorphism");
IsMonomorphism:=AttributeDispatcher("IsMonomorphism");
IsNilpotent:=AttributeDispatcher("IsNilpotent");
IsNormalized:=AttributeDispatcher("IsNormalized");
IsParent:=AttributeDispatcher("IsParent");
IsPerfect:=AttributeDispatcher("IsPerfect");
IsSimple:=AttributeDispatcher("IsSimple");
IsSolvable:=AttributeDispatcher("IsSolvable");
IsSurjective:=AttributeDispatcher("IsSurjective");
IsTrivial:=AttributeDispatcher("IsTrivial");
KernelGroupHomomorphism:=AttributeDispatcher("KernelGroupHomomorphism");
Lattice:=AttributeDispatcher("Lattice");
LatticeSubgroups:=AttributeDispatcher("LatticeSubgroups");
LowerCentralSeries:=AttributeDispatcher("LowerCentralSeries");
MaximalElement:=AttributeDispatcher("MaximalElement");
MaximalNormalSubgroups:=AttributeDispatcher("MaximalNormalSubgroups");
MinimalGeneratingSet:=AttributeDispatcher("MinimalGeneratingSet");
Normalized:=AttributeDispatcher("Normalized");
NormalSubgroups:=AttributeDispatcher("NormalSubgroups");
Omega:=AttributeDispatcher("Omega");
One:=AttributeDispatcher("One");
Radical:=AttributeDispatcher("Radical");
RationalClasses:=AttributeDispatcher("RationalClasses");
Representative:=AttributeDispatcher("Representative");
RepresentativesPerfectSubgroups:=
  AttributeDispatcher("RepresentativesPerfectSubgroups");
Resultant:=AttributeDispatcher("Resultant");
SizesConjugacyClasses:=AttributeDispatcher("SizesConjugacyClasses");
SmallestGenerators:=AttributeDispatcher("SmallestGenerators");
SupersolvableResiduum:=AttributeDispatcher("SupersolvableResiduum");
SylowComplements:=AttributeDispatcher("SylowComplements");
TrivialSubgroup:=AttributeDispatcher("TrivialSubgroup");
UpperCentralSeries:=AttributeDispatcher("UpperCentralSeries");
Valuation:=Dispatcher("Valuation");
Zero:=AttributeDispatcher("Zero");

##############################################################################
##
##  From here on the functions are dispatchers of the second type.
##

##############################################################################
##
#F  Endomorphisms( <arg> )
##
Endomorphisms:=function(arg)
  if Length(arg)=1 then
       return ApplyFunc(AttributeDispatcher("Endomorphisms"),arg);
  else return ApplyFunc(Dispatcher("Endomorphisms"),arg);
  fi;
end;

##############################################################################
##
#F  Dimensions( <struct> )
##
Dimensions := function( struct )
    if struct = [] then
      return [ 0, 0 ];
    elif     IsList( struct )
         and IsBound( struct[1] )
         and IsList( struct[1] ) then

      # Allow not only matrices but also lists of lists of matrices ...
      return [ Length( struct ), Length( struct[1] ) ];
    else
      return AttributeDispatcher("Dimensions")(struct);
    fi;
end;

##############################################################################
##
#F  Determinant( <struct> )
##
Determinant := function( struct )
    if   IsMat( struct ) then return DeterminantMat( struct );fi;
    return AttributeDispatcher("Determinant")(struct);
end;

##############################################################################
##
#F  Rank( <struct> )
##
Rank := function( struct )
    if   IsMat( struct ) then return RankMat( struct );fi;
    return AttributeDispatcher("Rank")(struct);
end;

##############################################################################
##
#F  Transposed( <struct> )
##
Transposed := function( struct )
    if   IsMat( struct ) then return TransposedMat( struct ); fi;
    return AttributeDispatcher("Transposed")(struct);
end;

#############################################################################
##
#F  AutomorphismGroup( <struct> )
##
AutomorphismGroup := function ( struct ) local d;
    if IsRec( struct ) then
      return AttributeDispatcher("AutomorphismGroup")(struct);
    else
      d:=Domain([struct]);
      if IsBound(d.operations.AutomorphismGroup) then
        return d.operations.AutomorphismGroup(struct);
      else
	Error("sorry, don't know how to compute the automorphism group");
      fi;
    fi;
end;

##############################################################################
##
#F  IsMonomial( <struct> )
##
IsMonomial := function( struct )
    if IsRec( struct ) then
      return AttributeDispatcher("IsMonomial")(struct);
    else
      return Domain( [ struct ] ).operations.IsMonomial( struct );
    fi;
end;

##############################################################################
##
##  From here on the functions are individual dispatchers.
##
##############################################################################
##
#F  Components( <struct> )
##
Components := function( arg )
    if 1 = Length(arg)  then
         return AttributeDispatcher("Components")(arg[1]);
    else return ApplyFunc(Dispatcher("Components"),arg);
    fi;
end;

##############################################################################
##
#F  Basis( <struct> )
#F  Basis( <struct>, <vectors> )
##
Basis := function( arg )
    local struct;

    if Length( arg ) < 1 or not IsRec( arg[1] ) then
      Error( "usage: Basis( <struct> [, <vectors>] )" );
    fi;

    struct:= arg[1];
    if   Length( arg ) = 1 then
      if not IsBound( struct.basis ) then
        if not IsBound( struct.operations ) then
          Error( "first argument must have 'operations' entry" );
        fi;
        struct.basis:=
            struct.operations.Basis( struct );
      fi;
      return struct.basis;
    elif Length( arg )>=2 then
      if not IsBound( struct.operations ) then
        Error( "first argument must have 'operations' entry" );
      fi;
      return ApplyFunc(struct.operations.Basis,
          Concatenation([struct], arg{[2..Length(arg)]}));
    fi;
    end;

##############################################################################
##
#F  StandardBasis( <struct> )
#F  StandardBasis( <struct>, <seed> )
##
StandardBasis := function( arg )
    local struct;

    if Length( arg ) < 1 or Length( arg ) > 2 or not IsRec( arg[1] ) then
      Error( "usage: StandardBasis( <struct> [, <seed>] )" );
    fi;

    struct:= arg[1];
    if   Length( arg ) = 1 then
      if not IsBound( struct.standardBasis ) then
        if not IsBound( struct.operations ) then
          Error( "first argument must have 'operations' entry" );
        fi;
        struct.standardBasis:=
            struct.operations.StandardBasis( struct );
      fi;
      return struct.standardBasis;
    elif Length( arg ) = 2 then
      if not IsBound( struct.operations ) then
        Error( "first argument must have 'operations' entry" );
      fi;
      return struct.operations.StandardBasis( struct, arg[2] );
    fi;
    end;

##############################################################################
##
#F  Display( <obj> [,<arec>] )
##
Display := function( arg )

    if Length( arg ) = 1 then

      if not IsRec( arg[1] ) or not IsBound( arg[1].operations ) then
        Error( "<obj> must be a record with 'operations' entry" );
      fi;
      arg[1].operations.Display( arg[1], rec() );

    elif Length( arg ) = 2 then

      if not IsRec( arg[1] ) or not IsBound( arg[1].operations ) then
        Error( "<obj> must be a record with 'operations' entry" );
      fi;
      arg[1].operations.Display( arg[1], arg[2] );

    else
      Error( "usage: Display( <obj> [, <record>] )" );
    fi;
    end;

#############################################################################
##
#F  IsIrreducible( <obj> )
#F  IsIrreducible( [<R>,] <obj> )
##
##  irreducibility of <obj>: character, module, representation, ring element
##
IsIrreducible := function( arg )
    local R, r;

    # get and check the arguments
    if   Length( arg ) = 1 then

        # group module from smash
        # gap records
        if IsRec( arg[1] ) then
           return ApplyFunc(AttributeDispatcher("IsIrreducible"),arg);

        # anything else
        else
            R := DefaultRing( arg[1] );
            r := arg[1];
        fi;

    elif Length(arg) = 2 then
        R := arg[1];
        if not IsRing(R)  then
            Error( "<R> must be a ring" );
        fi;
        r := arg[2];
        if not r in R  then
            Error( "<r> must be an element of <R>" );
        fi;
    else
        Error( "usage: IsIrreducible( <struct> )",
               " resp. IsIrreducible( [<R>,] <r> )" );
    fi;

    # permform the test
    return R.operations.IsIrreducible( R, r );

end;

#############################################################################
##
#F  IsEquivalent( <obj1>, <obj2> )
##
IsEquivalent := function( obj1, obj2 )
    return IsIdentical( obj1, obj2 ) or
           obj1.operations.IsEquivalent( obj1, obj2 );
    end;

#############################################################################
##
#F  Kernel( <obj> ) . . . . . . . kernel of a homomorphism, character, fusion
##
##  very preliminary .... (hopefully)
##
Kernel := function( obj )
    if not IsBound( obj.kernel )  then
      if IsBound( obj.operations.Kernel ) then
        obj.kernel:= obj.operations.Kernel( obj );
      else
        obj.kernel:= obj.source.operations.Kernel( obj );
      fi;
    fi;
    return obj.kernel;
    end;

#############################################################################
##
#F  FusionConjugacyClasses( <S>, <R> )
##
FusionConjugacyClasses := function( S, R )
    return S.operations.FusionConjugacyClasses( S, R );
    end;

#############################################################################
##
#F  KroneckerProduct(<M1>,...,<Mn>) . . . .  Kronecker product of n objects
##
KroneckerProduct:=function(arg)local res, i, Kronecker2;
  Kronecker2:=function(obj1,obj2)
    if IsMat( obj1 ) then
      if IsMat( obj2 ) then return MatricesOps.KroneckerProduct( obj1, obj2 );
      else return obj2.operations.KroneckerProduct( obj1, obj2 );
      fi;
    else return obj2.operations.KroneckerProduct( obj1, obj2 );
    fi;
  end;
  res:=arg[1];
  for i in [2..Length(arg)] do res:=Kronecker2(res,arg[i]);od;
  return res;
end;

#############################################################################
##
#F  Closure(<G>,<obj>)  . closure of struct with another struct or an element
##
#T  maybe the 'IsRec' calls should be replaced.
##
Closure := function ( G, obj )
    local   C;          # closure of <G> with <obj>, result

    # check the arguments
    if not IsRec( G )  then
      Error( "<G> must be a record" );
    elif not IsRec( obj ) and not obj in Parent( G )  then
      Error( "<obj> must be an element of the parent of <G>" );
    elif IsRec( obj ) and not obj in Parent( G )
                      and Parent( obj ) <> Parent( G ) then
      Error( "<obj> must be a substructure of the parent of <G>" );
    fi;

    # compute the closure
    C := G.operations.Closure( G, obj );

    # return the closure
    return C;
    end;

#############################################################################
##
#F  Centralizer( <G>, <obj> ) . . centralizer of a substructure or an element
##
Centralizer := function ( G, obj )
    local   C;          # centralizer of <obj> in <G>, result

    # check the arguments
    if not IsRec( G )  then
      Error( "<G> must be at least a record" );

#T What shall I do with this hack (invented by Alexander)?

    elif IsGroup( G ) and not IsGroup( obj ) and not obj in Parent( G )  then

	# centralizer for elements in not defined inifinite supergroup
	if IsBound(Domain([obj]).operations.Group)  # test for possible
	                                            # group object
	   and obj^0 in Parent(G) then
          C:=Stabilizer(G,obj,function(e,g) return e^g;end);
	  return C;
	elif IsBound(G.operations.Centralizer) then
          return G.operations.Centralizer( G, obj );
        else
	  Error( "<obj> must be an element of the parent of <G>" );
	fi;
    elif IsGroup( G ) and IsGroup( obj ) and Parent( obj ) <> Parent( G )  then
      if obj.identity=G.identity then
	# centralizer for subgroup of not defined inifinite supergroup
	return Intersection(List(obj.generators,i->Centralizer(G,i)));
      else
        Error( "<obj> must be a subgroup of the parent of <G>" );
      fi;
    fi;

#T end of the hack ...

    # compute the centralizer
    if IsGroup( obj ) and G = Parent( obj ) then
      if not IsBound( obj.centralizer ) then
        obj.centralizer:= G.operations.Centralizer( G, obj );
      fi;
      C:= obj.centralizer;
    else
      C:= G.operations.Centralizer( G, obj );
    fi;

    # return the centralizer
    return C;
    end;

#############################################################################
##
#F  IsCentral( <G>, <U> ) . . . test if a structure is centralized by another
##
IsCentral := function ( G, U )
    local   isCen;      # 'true' if <U> is central in <G>, result

    # check that the arguments are groups with a common parent
    if not IsRec( G )  then
        Error( "<G> must be a group" );
    elif not IsRec( U )  then
        Error( "<U> must be a group" );
    fi;

    # <U> and <G> must have a common parent
    Parent( G, U );

    # if <G> is the parent, use the entry '<U>.isCentral'
    if IsParent( G )  then
        if not IsBound( U.isCentral )  then
            U.isCentral := G.operations.IsCentral( G, U );
        fi;
        isCen := U.isCentral;

    # otherwise
    else
        if IsBound( U.isCentral ) and U.isCentral  then
            isCen := true;
        else
            isCen := G.operations.IsCentral( G, U );
        fi;
    fi;

    # return the result
    return isCen;
end;

#############################################################################
##
#F  Base( <object> )  . . . . . . . . . . . . . . . . . . . base of an object
##
##  Determines  a base of an object by either returning a corresponding entry
##  in  the  domain record  or  by  calculating it, using the function in the
##  operations  record.  In the latter case the record field 'isComputedBase'
##  is set to (possibly) indicate special properties of the base.
##
Base := function( obj )
    local   base;

#T remove this bad hack as soon as 'Base' is for permutation groups only!
    if   IsPermGroup( obj ) then
      base:= obj.operations.Base( obj );
    elif IsRec( obj ) then
      if   IsBound( obj.base ) then
        base:= obj.base;
      elif IsBound( obj.operations.Base ) then
        obj.base:= obj.operations.Base( obj );
        obj.isComputedBase:= true;
        base:= obj.base;
      else
        Error( "sorry, can't compute a base for <object>" );
      fi;
    else
      base:= Domain( [ obj ] ).operations.Base( obj );
    fi;
    return base;
    end;

##############################################################################
##
#F  IsZero( <struct>, <obj> ) . . . . . . .  is an object the zero of a domain
##
IsZero := function( struct, obj )
    if not IsRec( struct ) then
      Error( "<struct> must be a record" );
    fi;
    if   IsBound( struct.zero ) then
      return obj = struct.zero;
    elif IsBound( struct.operations ) then
      return struct.operations.IsZero( struct, obj );
    else
      Error( "<struct> must have 'operations' entry" );
    fi;
    end;

##############################################################################
##
#F  Eigenvalues( <mat> )
#F  Eigenvalues( <classfun>, <class> )
#F  Eigenvalues( <tbl>, <charvals>, <class> )
#F  Eigenvalues( ... )
##
## JM 2015/02/25 ported from chevie patch for arbitrary arg length
Eigenvalues := function( arg )
  if Length(arg)=0 then Error( "need at least one argument" );fi;
  if Length(arg)=1 and IsMat(arg[1]) then
    return Domain(arg).operations.Eigenvalues(arg[1]);fi;
  return ApplyFunc(Dispatcher("Eigenvalues"),arg);
end;

#############################################################################
##
#F  TransitiveClosure( <M> ) . . . . . Transitive closure
##  
##  if <M> is a square boolean matrix representing a relation, returns
##  a  boolean matrix representing the transitive closure of this relation.
##  Uses Floyd-Warshall algorithm so is quite fast.
##
##  Otherwise is a dispatcher for TransitiveClosure.
#
TransitiveClosure:=function(m)local i,k;
  if IsList(m) then
    m:=Copy(m);
    for k in [1..Length(m)] do for i in [1..Length(m)] do 
      if m[i][k] then UniteBlist(m[i],m[k]);fi;
    od;od;
    return m;
  else 
    return Dispatcher("TransitiveClosure")(m);
  fi;
end;
