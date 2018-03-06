###############################################################################
##
##  t.g  GLISSANDO ver 1.0  Christof Noebauer   1995, 1996
##                                                       
##
#############################################################################
##
#V  TRANSFORMATION_PRINT_LEVEL: this global variable determines how
##                              detailed transformations are displayed.
##
TRANSFORMATION_PRINT_LEVEL := -1; # set default: display detailed

#############################################################################
##
#F  TransformationPrintLevel( <arg> ). . set print detail for transformations
##
TransformationPrintLevel := function( arg )
  if Length( arg ) > 1 or 
     ( Length( arg ) = 1 and not IsString( arg[1] ) ) then
    Error( "Usage: TransformationPrintLevel() or \n", 
           "TransformationPrintLevel( <string> ) where <string> must be\n",
           "\"long\" or \"short\" or \"l\" or \"s\"" );
  fi;

  if Length( arg ) = 0 then
    if TRANSFORMATION_PRINT_LEVEL = 0 then
      return "short";
    else
      return "long";
    fi;
  else
    if arg[1][1] = 'l' then
      TRANSFORMATION_PRINT_LEVEL := -1;
    else
      TRANSFORMATION_PRINT_LEVEL := 0;
    fi;
  fi;
  return;
end;

TPL := TransformationPrintLevel;

#############################################################################
##
#F  IsTransformation( <obj> ) . . . . . . . test if <obj> is a transformation
##
IsTransformation := function( obj )
  return IsRec( obj )
    and IsBound( obj.isTransformation ) and obj.isTransformation;
end;

#############################################################################
##
#F  IsSetTransformation( <obj> ). . . . test if <obj> is a set transformation
##
IsSetTransformation := function( obj )
  return IsRec( obj )
    and IsBound( obj.isSetTransformation ) and obj.isSetTransformation;
end;

#############################################################################
##
#F  IsGroupTransformation( <obj> ). . test if <obj> is a group transformation
##
IsGroupTransformation := function( obj )
  return IsRec( obj )
    and IsBound( obj.isGroupTransformation ) and obj.isGroupTransformation;
end;

#############################################################################
##
#V  TransformationOps . . . . . . . . . operations record for transformations
##
##  'TransformationOps' is the operations record for transformations.  This 
##  is initially a copy of 'MappingOps'.  This way all the default methods 
##  for mappings are inherited.                                                        ##                                                                     
SetTransformationOps   := Copy( MappingOps );
GroupTransformationOps := Copy( MappingOps ); 

SetTransformationOps.name   := "SetTransformationOps";
GroupTransformationOps.name := "GroupTransformationOps";

#############################################################################
##
#F  AsTransformation( <T> ). . . . . . . make a mapping into a transformation
##
##  Function for making any mapping with the same range and source into a 
##  Transformation.
##
AsTransformation := function( T )
  local t,l;
  if not ( IsRec( T ) and IsMapping( T ) and ( T.source = T.range ) and
           ( IsGroup( T.source ) or IsSet( T.range ) ) ) then
    Error( "Usage: AsTransformation( <T> ) where <T> must be a mapping with",
           "\nthe equal source and range which must be a set or a ",
           "group" );
  fi;

  # make the transformation
  t := Copy( T );
  # enter the tag components
  t.isTransformation        := true;
  if IsSet( T.source ) then
    t.isSetTransformation   := true;
    t.elements              := T.source;
  elif IsGroup( T.source ) then
    t.isGroupTransformation := true;
    t.isGroupElement        := true;
    t.elements              := Elements( T.source );
  fi;
  # make the transformation list
  l := List( t.elements, e -> Position( t.elements, Image( T, e ) ) );
  t.tfl                     := l; # a list which represents the transformation
                                  # in the obvious way. ('TransFormationList')
  # enter the operations record
  if IsSet( T.source ) then
    t.operations            := SetTransformationOps;
  elif IsGroup( T.source ) then
    t.operations            := GroupTransformationOps;
  fi;
  
  return t;
end;

#############################################################################
##
#F  IdentityTransformation( <D> ). . . . . create the identity transformation
##
##  Constructor function for the identity Transformation on <D> where <D> 
##  is either a set or a group.
##
IdentityTransformation_gliss := function( D )
  local t,l;
  # check on the validity of the input
  if not ( IsSet( D ) or IsGroup( D ) ) then
    Error( "Usage: IdentityTransformation( <D> ) where <D> must be a set\n",    
           "or a group" );
  fi;
  
  l := List( [1..Size( D )], i -> i );
  # make the transformation
  t := rec();
  # enter the tag components
  t.isGeneralMapping        := true;
  t.domain                  := Mappings;
  t.isMapping               := true;
  t.isTransformation        := true;
  if IsSet( D ) then
    t.isSetTransformation   := true;
    t.elements              := D;
  elif IsGroup( D ) then
    t.isGroupTransformation := true;
    t.isGroupElement        := true;
    t.elements              := Elements( D );
  fi;
  # enter source and range
  t.source                  := D;
  t.range                   := D;
  t.tfl                     := l; # a list which represents the transformation
                                  # in the obvious way. ('TransFormationList')
  # enter the operations record
  if IsSet( D ) then
    t.operations            := SetTransformationOps;
  elif IsGroup( D ) then
    t.operations            := GroupTransformationOps;
  fi;
  # return the transformation
  return t;
end;

#############################################################################
##
#F  Transformation( <D>, <l> ). . . . . . . . . . . . create a transformation
##
##  Constructor function for a Transformation on <D> where <D> is either a
##  set or a group
##
Transformation_gliss := function( D, l )
  local t;
  # check on the validity of the input
  if not ( ( IsSet( D ) or IsGroup( D ) ) and IsVector( l ) and 
           ForAll( l, x -> x in [1..Size( D )] ) and 
           Length( l ) = Size( D )
         ) then
    Error( "Usage: Transformation( <obj>, <tfl> ) where <obj> must be a ", 
           "set\nor a group and <tfl> must be a valid ", 
           "transformation list for <obj>" );
  fi;
  
  # make the transformation
  t := rec();
  # enter the tag components
  t.isGeneralMapping        := true;
  t.domain                  := Mappings;
  t.isMapping               := true;
  t.isTransformation        := true;
  if IsSet( D ) then
    t.isSetTransformation   := true;
    t.elements              := D;
  elif IsGroup( D ) then
    t.isGroupTransformation := true;
    t.isGroupElement        := true;
    t.elements              := Elements( D );
  fi;
  # enter source and range and the transformation list
  t.source                  := D;
  t.range                   := D;
  t.tfl                     := l; # a list which represents the transformation
                                  # in the obvious way. ('TransFormationList')
  # enter the operations record
  if IsSet( D ) then
    t.operations            := SetTransformationOps;
  elif IsGroup( D ) then
    t.operations            := GroupTransformationOps;
#    t.operations.ImageRepresentative := 
#                 ConjugationGroupHomomorphismOps.ImageRepresentative;
#    t.operations.PreImageRepresentative := 
#                 ConjugationGroupHomomorphismOps.PreImageRepresentative;
  fi;
  
  # return the transformation
  return t;
end;

#############################################################################
##
#F  DisplayTransformation( <t> )  . . . . . . . . . display a transformation
##
##  Display a transformation <t>.
##
DisplayTransformation := function( t )
  local elm;
  if TRANSFORMATION_PRINT_LEVEL <> 0 then
    if IsSetTransformation( t ) then
      Print( "Transformation on ", t.elements, ":\n" );
      for elm in t.elements do
        Print( "  ", elm, " -> ", t.tfl[elm], "\n" );
      od;
    else
      Print( "Transformation on ", t.source, ":\n" );
      for elm in t.elements do
        Print( "  ", elm, " -> ", Image( t, elm ), "\n" );
      od;
    fi;
  else
    Print( t.tfl, " " );
  fi;
end;

SetTransformationOps.ImageElm := function( t, elm )
  return t.elements[ t.tfl[ Position( t.elements, elm ) ] ];
end;

SetTransformationOps.ImagesElm := function( t, elm )
  return [ t.elements[ t.tfl[ Position( t.elements, elm ) ] ] ];
end;

SetTransformationOps.ImagesRepresentative := function( t, elm )
  return t.elements[ t.tfl[ Position( t.elements, elm ) ] ];
end;

SetTransformationOps.ImagesSource := function( t )
  return Set( Filtered( t.elements, 
              elm -> Position( t.elements, elm ) in t.tfl ) );
end;

SetTransformationOps.Print := function( t )
  Print( "Transformation( ", t.source, ", ", t.tfl, " )" );
end;

SetTransformationOps.\* := function( t1, t2 )
  local prd;
  # product of two transformations ( IF both arguments are transformations )
  if IsRec( t1 ) and IsRec( t2 ) and
    IsTransformation( t1 ) and IsTransformation( t2 ) and 
    t1.source = t2.range then
      prd := Transformation( t1.source, t2.tfl{ t1.tfl } );
  else 
    # delegate to MappingOps.\*
    prd := MappingOps.\* ( t1, t2 );
  fi;
  return prd;
end;

SetTransformationOps.Kernel := function( t )
  local ker, # help variable: a list in which the eq.classes of the 
             # relation Ker are put together
        elm; # help variable
  ker := [];
  for elm in Image( t ) do
    AddSet( ker, Filtered( t.elements, e -> Image( t, e ) = elm ) ); 
  od;
  return ker;
end;

SetTransformationOps.Rank := function( t )
  return Length( Image( t ) );
end;

GroupTransformationOps := Copy( SetTransformationOps );

GroupTransformationOps.\+ := function( t1, t2 )
  local tfl,  # the tfl of the resulting transformation 
        elms, # a list of the elements involved
        l;    # number of elements
  # sum of two group transformations ( IF both arguments are group tf's )
  if IsRec( t1 ) and IsRec( t2 ) and
    IsGroupTransformation( t1 ) and IsGroupTransformation( t2 ) and 
    t1.source = t2.source and t1.range = t2.range then
      elms := t1.elements;
      l    := Length( elms );
      tfl  := List( [1..l], i -> 
                    Position( elms, elms[ t1.tfl[i] ] * elms[ t2.tfl[i] ] )
                  );
  else 
    Error( "Usage: <t1> + <t2> where <t1> and <t2> must both be group\n", 
           "transformations on the same group" );
  fi;
  return Transformation( t1.source, tfl );
end;

GroupTransformationOps.\- := function( t1, t2 )
  local tfl,  # the tfl of the resulting transformation 
        elms, # a list of the elements involved
        l;    # number of elements
  # difference of two group tf's ( IF both arguments are group tf's )
  if IsRec( t1 ) and IsRec( t2 ) and
    IsGroupTransformation( t1 ) and IsGroupTransformation( t2 ) and 
    t1.source = t2.source and t1.range = t2.range then
      elms := t1.elements;
      l    := Length( elms );
      tfl  := List( [1..l], i -> 
                    Position( elms, elms[ t1.tfl[i] ] * elms[ t2.tfl[i] ]^-1 )
                  );
  else 
    Error( "Usage: <t1> - <t2> where <t1> and <t2> must both be group\n", 
           "transformations on the same group" );
  fi;
  return Transformation( t1.source, tfl );
end;

