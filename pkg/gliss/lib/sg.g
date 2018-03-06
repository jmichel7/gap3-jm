###############################################################################
##
##  sg.g  GLISSANDO ver 1.0  Christof Noebauer   1995, 1996
##                                                       
##
#############################################################################
##
#V  SgOps. . . . . . . . . . . . . . . . . . operations record for semigroups
##
##  'SgOps' is the operations record for semigroups.  This is initially a 
##  copy of 'DomainOps'.  This way all the default methods for domains are 
##  inherited.                                                        
##                                                                     
SgOps := Copy( DomainOps );

#############################################################################
##
#F  IsSemigroup( <S> ) . . . . . . . . . . . . . . test if <S> is a semigroup
##
IsSemigroup := function( S )
  return IsRec( S ) and IsBound( S.isSemigroup ) and S.isSemigroup; 
end;

#############################################################################
##
#F  IsTransformationSemigroup( <S> ). . . . . . test if <S> is a tf semigroup
##
IsTransformationSemigroup := function( S )
  return IsRec( S ) and 
    IsBound( S.isTransformationSemigroup ) and S.isTransformationSemigroup; 
end;

#############################################################################
##
#F  AllFunctions( <n> ) . create all f's from the set {1,2,...,n} into itself 
##
##  Those functions are returned as lex. ordered list with [1,...,1] first, 
##  and with the constant function [n,n,...,n] last.
##
AllFunctions := function( n )
  local af,i,j,k,done;

  i := []; 
  for j in [1..n] do i[j] := 1; od;
  k := -1;
  af := [];
  
  while k <= n do
    AddSet( af, Reversed( i ) );
    k := 1; done := false;
    while not done and k <= n do
      if i[k] = n then  i[k] := 1; k := k+1; 
      else  i[k] := i[k]+1; done := true;
      fi;
    od;
  od;

  return af;
end;

#############################################################################
##
#F  TransformationSemigroup(<gen>,...). . . . . . . . . . .create a semigroup
##
##  Constructor function for transformation semigroups
##
TransformationSemigroup := function( arg )
  
  local S,    # the semigroup to be returned
        t,    # transformations
        set,  # help variable: a set
        gens; # generators of the semigroup
  
  arg := Flat( arg );
  S := rec();

  if Length( arg ) = 1 and IsInt( arg[1] ) and arg[1] > 0 then

    gens := []; set := List( [1..arg[1]], i -> i );
    for t in AllFunctions( arg[1] ) do
      AddSet( gens, Transformation( set, t ) );
    od;
    S.elements := gens;
   
  else

    # check the validity of the input
    if not ForAll( arg, t -> IsTransformation( t ) and 
                        t.source = arg[1].source and t.range = arg[1].range
                 ) then
      Error( "Usage: TransformationSemigroup( <t1>, <t2>, ... ) where all ",
             "arguments\nmust be transformations on the same set or \n", 
             "TransformationSemigroup( <n> ) where <n> must be a positive ", 
             "integer" );
    fi;
  
    # make sure that all sources and ranges of the generators are not only
    # equal but identical.
    gens := [];
    for t in arg do
      t.source := arg[1].source; t.range := arg[1].source;
      AddSet( gens, t );
    od;

  fi;

  # make the domain
  S.isDomain                  := true;
  S.isSemigroup               := true;
  S.isTransformationSemigroup := true;
  S.generators                := gens;
  S.multiplication            := function( x, y ) return x * y; end;
  S.operations                := SgOps;

  return S;
end;

#############################################################################
##
##  SgOps.\=( <S1>, <S2> ). . . . . . . compare two transformation semigroups 
##
SgOps.\= := function( S1, S2 )
  if IsTransformationSemigroup( S1 ) and IsTransformationSemigroup( S2 ) then
    return ForAll( S1.generators, x -> x in S2 ) and
           ForAll( S2.generators, x -> x in S1 );
  else
    return DomainOps.\= ( S1, S2 );
  fi;
end;

#############################################################################
##
##  SgOps.Print( <S> ) . . . . . . . . . . . print a transformation semigroup
##
SgOps.Print := function( S )
  local i;
  Print( "TransformationSemigroup( " );
  for i in [ 1..Length( S.generators ) - 1 ] do 
    Print( S.generators[i], ", " ); 
  od;
  Print( S.generators[ Length( S.generators )], " ) " );
end;

#############################################################################
##
#F  SgOps.Elements( <S> ). compute the elements of a transformation semigroup
##
SgOps.Elements := function ( S )
  local tset,             # set of transformation lists
        previously_added, # set of tf's added in the previous loop execution
        newly_added,      # set of tf's added in the actual loop execution
        t1,t2,t;          # help variables: transformation lists

  if IsBound( S.elements ) then 
    return S.elements;
  else
    # perform a simple closure algorithm
    tset := Set( List( S.generators, t -> t.tfl ) ); 
    previously_added := Copy( tset ); 
    
    repeat
      newly_added := []; 
      for t1 in tset do
        for t2 in previously_added do
          t := t1{ t2 };
          if not t in tset and not t in previously_added 
            and not t in newly_added
          then 
            AddSet( newly_added, t );      # Print( "Adding ", t, "\n" ); 
          fi;
          if t1 <> t2 then
            t := t2{ t1 };
            if not t in tset and not t in previously_added 
              and not t in newly_added
            then 
              AddSet( newly_added, t );    # Print( "Adding ", t, "\n" ); 
            fi;
          fi;
        od;
      od;
      for t1 in previously_added do
        for t2 in previously_added do
          t := t1{ t2 };
          if not t in tset and not t in previously_added 
            and not t in newly_added
          then 
            AddSet( newly_added, t );      # Print( "Adding ", t, "\n" ); 
          fi;
        od;
      od;
      UniteSet( tset, previously_added ); 
      previously_added := Copy( newly_added ); 

#      Print("size: ", Length( tset ) + Length( newly_added ), "\n");
    
    until newly_added = [];

    if not IsBound( S.size ) then S.size := Length( tset ); fi;
    
    S.elements := [];
    for t in tset do
      AddSet( S.elements, Transformation( S.generators[1].source, t ) );
    od;
    
    if not IsBound( S.rank ) then S.rank := S.operations.Rank( S ); fi;
    
    return S.elements;
  fi;
end;

#############################################################################
##
#F  SgOps.Rank( <S> ). . . . . . . . . . . Rank of a transformation semigroup
##
##  This function returns the rank of a transformation semigroup S,i.e. 
##  min(rank(t)) for all transformations t in S (see Lallement p. 20)
##
SgOps.Rank := function ( S )
  return Minimum( List( Elements( S ), t -> Length( Set( t.tfl ) ) ) );
end;

#############################################################################
##
#F  SmallestIdeal( <S> ) . . . . smallest ideal of a transformation semigroup
##  Dispatcher function  
##
SmallestIdeal := function( S )
  
  if not IsTransformationSemigroup( S ) then
    Error( "Usage: SmallestIdeal( <S> ) where <S> must be a ",
           "\ntransformation semigroup" );
  fi;
  
  if not IsBound( S.smallestIdeal ) then
    S.smallestIdeal := S.operations.SmallestIdeal( S );
  fi;
  
  return S.smallestIdeal;
end;

#############################################################################
##
#F  SgOps.IsSimple( <S> ) . check if a transformation semigroup <S> is simple
##
SgOps.IsSimple := function( S )  

  return S = SmallestIdeal( S );

end;

#############################################################################
##
#F  SgOps.SmallestIdeal( <S> ) . smallest ideal of a transformation semigroup
##
##  this function computes the smallest ideal of a transformation semigroup S.
##  Take all the transformations of minimal rank.
##
SgOps.SmallestIdeal := function ( S )
  local r,     # help variable: rank of the semigroup S
        minid; # help variable: the smallest ideal

  r     := S.operations.Rank( S ); 
  minid := Filtered( Elements( S ), t -> Length( Set( t.tfl ) ) = r );
  return minid;
end;

#############################################################################
##
#F  SgOps.PrincipalLeftIdeal( <S>, <t> ) . . . . . . . . . . . . . . . . . .
##
##  Compute the principal left ideal of a transformation semigroup generated 
##  by the element <t>.
##
SgOps.PrincipalLeftIdeal := function( S, t )
  local elms,elm,pli,PLI;
  
  elms := Elements( S ); 
  pli  := [ t.tfl ];
  
  for elm in elms do
    AddSet( pli, elm.tfl{ t.tfl } );
  od;

  PLI := [];
  for elm in pli do
    AddSet( PLI, Transformation( S.generators[1].source, elm ) );
  od;
  return PLI;
end;

#############################################################################
##
#F  SgOps.PrincipalRightIdeal( <S>, <t> ) . . . . . . . . . . . . . . . . . .
##
##  Compute the principal right ideal of a transformation semigroup generated 
##  by the element <t>.
##
SgOps.PrincipalRightIdeal := function( S, t )
  local elms,elm,pri,PRI;
  
  elms := Elements( S ); 
  pri  := [ t.tfl ];
  
  for elm in elms do
    AddSet( pri, t.tfl{ elm.tfl } );
  od;

  PRI := [];
  for elm in pri do
    AddSet( PRI, Transformation( S.generators[1].source, elm ) );
  od;
  return PRI;
end;

#############################################################################
##
#F  SgOps.PrincipalIdeal( <S>, <t> ). . . . . . . . . . . . . . . . . . . . .
##
##  Compute the principal ideal of a transformation semigroup generated 
##  by the element <t>.
##
SgOps.PrincipalIdeal := function( S, t )
  local elm,PI,tl;
  
  PI  := [ t ];
  
  tl := S.operations.PrincipalLeftIdeal( S, t ); 
  for elm in tl do
    UniteSet( PI, S.operations.PrincipalRightIdeal( S, elm ) );
  od;

  return PI;
end;

#############################################################################
##
#F  Green( <S>, <type> ) . . . . . compute Green's relations for a tf sg <S>
##
##  Dispatcher function for all of Green's relations
##
Green := function( S, type )

  if not ( IsTransformationSemigroup( S ) and 
           type in [ "L","l","R","r","J","j","D","d","H","h" ] ) then
    Error( "Usage: Green( <S>, <type> ) where <S> must be a transformation",
           "\nsemigroup and <type> must be one of the following: ", 
           "\"l\",\"r\",\"j\",\"d\", or \"h\" "); 
  fi;

  if type in [ "l","L" ] then
    if not IsBound( S.greenL ) then 
      S.greenL := S.operations.GreenL( S );
    fi;
    return S.greenL;
  elif type in [ "r","R" ] then
    if not IsBound( S.greenR ) then 
      S.greenR := S.operations.GreenR( S );
    fi;  
    return S.greenR;
  elif type in [ "j","J" ] then
    if not IsBound( S.greenJ ) then 
      S.greenJ := S.operations.GreenJ( S );
    fi;  
    return S.greenJ;
  elif type in [ "d","D" ] then
    if not IsBound( S.greenD ) then 
      S.greenD := S.operations.GreenD( S );
    fi;  
    return S.greenD;
  elif type in [ "h","H" ] then
    if not IsBound ( S.greenH ) then 
      S.greenH := S.operations.GreenH( S );
    fi;  
    return S.greenH;
  else
    Error( "panic, should never come here" );
  fi;
  
  return -1;
end;

#############################################################################
##
#F  SgOps.GreenL( <S> ). . . . . Green's left relation of a transformation Sg
##
SgOps.GreenL := function( S )
  local l,       # this holds the equivalence classes of Greens left relation  
        tf,      # help variable: a transformation in S
        classes, # set of classes
        class,cl,# one equivalence class of Green's left relation
        elms;    # holds the elements of S
  
  elms := Elements( S );
  
  # make a pre - classification
  classes := [];
  for tf in elms do
    AddSet( classes, Filtered( elms, t -> Kernel( t ) = Kernel( tf ) ) ); 
  od;
  l := Copy( classes );
  
  # refine the equivalence relation if necessary
  cl := Length( S.generators[1].source );
  if Size( S ) < cl^cl then 
    l := [];
    for class in classes do
      elms := Copy( class );
      repeat
        tf := elms[1];
        cl := Filtered( elms, t -> S.operations.PrincipalLeftIdeal( S, t ) = 
                                   S.operations.PrincipalLeftIdeal( S, tf ) );
        AddSet( l, cl ); 
        SubtractSet( elms, cl );
      until elms = [];
    od;
  fi;

  return l;
end;

#############################################################################
##
#F  SgOps.GreenR( <S> ) . . . . Green's right relation of a transformation Sg
##
SgOps.GreenR := function( S )
  local r,       # this holds the equivalence classes of Greens right relation  
        tf,      # help variable: a transformation in S
        classes, # set of classes
        class,cl,# one equivalence class of Green's right relation
        elms;    # holds the elements of S
  
  elms := Elements( S );
  
  # make a pre - classification
  classes := [];
  for tf in elms do
    AddSet( classes, Filtered( elms, t -> Image( t ) = Image( tf ) ) ); 
  od;
  r := Copy( classes );
  
  # refine the equivalence relation if necessary
  cl := Length( S.generators[1].source );
  if Size( S ) < cl^cl then 
    r := [];
    for class in classes do
      elms := Copy( class );
      repeat
        tf := elms[1];
        cl := Filtered( elms, t -> S.operations.PrincipalRightIdeal( S, t ) = 
                                   S.operations.PrincipalRightIdeal( S, tf ) 
                      );
        AddSet( r, cl ); 
        SubtractSet( elms, cl );
      until elms = [];
    od;
  fi;

  return r;
end;

#############################################################################
##
#F  SgOps.GreenJ( <S> ). . . . . . . join of Green's left and right relations 
##
SgOps.GreenJ := function( S )
  local j,       # this holds the equivalence classes of Greens D relation  
        l,r,     # Green's L resp. R relations
        class,   # a class of Green's left relation
        cl,      # a class of Green's J relation
        elm;     # help variable, an element
  
  l := Green( S, "L" );
  r := Green( S, "R" );
  
  j := [];
  for class in l do 
    cl := [];
    for elm in class do
      UniteSet( cl, First( r, c -> elm in c ) );
    od;
    AddSet( j, cl );
  od;

  return j;
end;

#############################################################################
##
#F  SgOps.GreenD( <S> ). . . . . . . . . . . . . . . . . . Green's D relation 
##  D = J in a finite semigroup!
##
SgOps.GreenD := SgOps.GreenJ;

#############################################################################
##
#F  SgOps.GreenH( <S> ) . . . . Intersection of Green's left and right rel's.
##
SgOps.GreenH := function( S )
  local h,        # this holds the equivalence classes of Greens H relation  
        l,r,      # Green's L resp. R relations    
        cl,       # a class of Green's h relation
        cl_left,  # a class of Green's L relation
        cl_right, # a class of Green's right relation
        elm,      # help variable, an element
        elms;     # holds the elements of S
  
  elms := Copy( Elements( S ) );
  l := Green( S, "L" );
  r := Green( S, "R" );
  
  h := [];
  repeat
    elm := elms[1]; 
    cl_left  := First( l, cl -> elm in cl );
    cl_right := First( r, cl -> elm in cl );
    cl := IntersectionSet( cl_left, cl_right );
    AddSet( h, cl );
    SubtractSet( elms, cl );
  until elms = [];

  return h;
end;

#############################################################################
##
#F  SgOps.DisplayTable( <S> ) . . . print a Cayley table of the semigroup <S>
##
SgOps.DisplayTable := function( S )

  local elms,    # the elements of the semigroup
        n,       # the size of the semigroup
        symbols, # a list of the symbols for the elements of the semigroup
        tw,      # the width of a table
        spc,     # local function which prints the right number of spaces
        bar,     # local function for printing the right length of the bar
        ind,     # help variable, an index
        i,j;     # loop variables
  
  elms    := Elements( S );
  n       := Size( elms );
  symbols := [ AbstractGenerator( "s0" ) ];
  if n > 1 then
    symbols := Concatenation( [ AbstractGenerator( "s0" ) ], 
                              AbstractGenerators( "s" , n-1 ) );
  fi;
  # compute the number of characters per line required for the table
  if n < 11 then tw := 3*n + 6; else tw := 4*n + 8; fi;
  spc     := function( i, max )
               if max < 11 then return " "; fi;
               if i < 11 then return "  "; else return " "; fi;
             end;
  bar     := function( max )
               if max < 11 then return "---"; else return "----"; fi;
             end;
  
# if SizeScreen()[1] - 3 < tw then
#   Print( "The table of a semigroup of order ", n, " will not ",
#          "look\ngood on a screen with ", SizeScreen()[1], " characters per ", 
#          "line.\nHowever, you may want to set your line length to a ",
#          "greater\nvalue by using the GAP function 'SizeScreen'.\n" );
#   return;
# fi;
  Print( "Let:\n" );
  for i in [1..n] do Print( symbols[i], " := ", elms[i], "\n" ); od;
  
  # print the multiplication table
  Print( "\n\n  * ", spc( 0, n ), "| " ); 
    for i in [1..n] do Print( symbols[i], spc( i, n ) ); od;
  Print( "\n ---", bar( n ) ); for i in [1..n] do Print( bar( n ) ); od;
  for i in [1..n] do
    Print( "\n  ", symbols[i], spc( i, n ), "| " );
    for j in [1..n] do
      ind := Position( elms, S.multiplication( elms[i], elms[j] ) );
      Print( symbols[ ind ], spc( ind, n ) );  
    od;
  od;
  
  Print( "\n\n" );
end;

#############################################################################
##
#F  SgOps.Identity( <S> ) . . . . . . . . compute identity of a semigroup <S>
##
SgOps.Identity := function( S )
  local id, elms;
  
  id := IdentityTransformation( S.generators[1].source );
  if id in Elements( S ) then
    return [ id ];
  else
    elms := Elements( S );
    for id in elms do
      if ForAll( elms, e -> e*id = e and id*e = e ) then
        return [ id ];
      fi;
    od;
  fi;
  
  return [ ];
end;

#############################################################################
##
#F  SgOps.IsCommutative( <S> ). . . . . . . . . . . test commutativity of <S>
##
SgOps.IsCommutative := function( S )
  local mul, elms;
  
  mul  := S.multiplication;
  elms := Elements( S );
  
  return ForAll( elms, c -> ForAll( elms, e -> mul( c, e ) = mul( e, c ) ) );
end;

#############################################################################
##
#F  SgOps.IdempotentElements( <S> ). . . . . . . . . . idempotent elms of <S> 
##
SgOps.IdempotentElements := function( S )
  local mul;

  mul  := S.multiplication;

  return Filtered( Elements( S ), i -> mul( i, i ) = i );

end;

#############################################################################
##
#F  LibrarySemigroup: function extracts a semigroup from the semigroups 
##                    library file SG.LIB
##  V1.1 1.3.95.
##
LibrarySemigroup := function( order, type )
  local i,j,          # loop variables
        clmax,        # the total number of isomorphism classes in a sg
        gens,         # the generators (=elements) of the semigroup
        f,            # the function representing the semigroup
        rho,          # a transformation being constructed
        af,           # help variable: containing the relevant functions
        sg,           # help variable: containing semigroup data
        num2fun;      # help variable: function to convert library shorthand
                      #                of a function into "the real thing".
  
  num2fun := function( num, b )
    local f, remainder, i;
    f := []; num := num - 1 ;
    for i in [1..b] do
      remainder := num mod b;
      num := ( num - remainder ) / b;
      Add( f, remainder + 1 );
    od;
    return Reversed( f );
  end;

  if not ( IsInt( order ) and IsInt( type ) and order > 0 and type > 0 ) then
    Error( "Usage: LibrarySemigroup( <order>, <type> ) where <order> and\n",
    "      <type> must both be positive integers" );
  fi;
  
  if   ( order = 1 ) then
    sg := SG1; 
  elif ( order = 2 ) then
    sg := SG2;
  elif ( order = 3 ) then
    sg := SG3;
  elif ( order = 4 ) then
    sg := SG4;
  elif ( order = 5 ) then
    sg := SG5;
  else
    Print( "There are only orders 1 to 5 in the semigroups library.\n" ); 
    return;
  fi;
  
  clmax := Length( RecFields( sg ) );
  if type > clmax then
    Print( "There are only ", clmax, 
           " isomorphism classes of semigroups of order ", order, ".\n" );
    return;
  fi;

  f    := sg.(type).phi;
  gens := [];
  af   := [];

  for i in [1..order] do
    if not IsBound( af[f[i]] ) then
      af[f[i]] := num2fun( f[i], order );
    fi;
  od;

  # check if the semigroup has an identity
  i := First( [1..order+1], j -> j = order+1 or
              ( af[f[j]] = [1..order] and
                ForAll( [1..order], x -> af[f[x]][j] = x ) )
            );

  if i <= order then  
    
    # the case if the sg has an identity
    
    for i in [1..order] do
      rho := [];
      for j in [1..order] do
        Add( rho, af[f[j]][i] );
      od;
      Add( gens, Transformation( [1..order], rho ) );
    od;
  
  else
  
    # the case if the sg has no identity
    
    for i in [1..order] do
      rho := [];
      for j in [1..order] do
        Add( rho, af[f[j]][i] );
      od;
      Add( rho, i );
      Add( gens, Transformation( [1..order+1], rho ) );
    od;
  
  fi;

  return TransformationSemigroup( gens );
end;

#############################################################################
##
#F  AllLibrarySemigroups: this function extracts all semigroups of a
##  specified class from the semigroups library file SG.LIB
##  V1.0 14.2.95.
##
AllLibrarySemigroups := function( order, type )
  local i,j,          # loop variables
        clmax,        # the total number of isomorphism classes in a sg
        gens,         # the generators (=elements) of the semigroup
        f,            # the function representing the semigroup
        autos,a,      # automorphisms
        rho,          # a transformation being constructed
        af,           # help variable: containing the relevant functions
        sgps,         # return value: a list of semigroups
        sg,           # help variable: containing semigroup data
        num2fun;      # help variable: function to convert library shorthand
                      #                of a function into "the real thing".
  
  num2fun := function( num, b )
    local f, remainder, i;
    f := []; num := num - 1 ;
    for i in [1..b] do
      remainder := num mod b;
      num := ( num - remainder ) / b;
      Add( f, remainder + 1 );
    od;
    return Reversed( f );
  end;
 
  if not ( IsInt( order ) and IsInt( type ) and order > 0 and type > 0 ) then
    Error( "Usage: AllLibrarySemigroups( <order>, <type> ) where <order>", 
           "and\n      <type> must both be positive integers" );
  fi;
  
  if   ( order = 1 ) then
    sg := SG1; 
  elif ( order = 2 ) then
    sg := SG2;
  elif ( order = 3 ) then
    sg := SG3;
  elif ( order = 4 ) then
    sg := SG4;
  elif ( order = 5 ) then
    sg := SG5;
  else
    Print( "There are only orders 1 to 5 in the semigroups library.\n" ); 
    return;
  fi;
  
  clmax := Length( RecFields( sg ) );
  if type > clmax then
    Print( "There are only ", clmax, 
           " isomorphism classes of semigroups of order ", order, ".\n" );
    return;
  fi;

  f     := sg.(type).phi;
  af    := [];
  sgps  := [];
  autos := [];
  for a in sg.(type).bijs_yielding_iso_sgps do
    Add( autos, num2fun( a, order ) );
  od;
  
  for i in [1..order] do
    if not IsBound( af[f[i]] ) then
      af[f[i]] := num2fun( f[i], order );
    fi;
  od;

  # check if the semigroup has an identity
  i := First( [1..order+1], j -> j = order+1 or
              ( af[f[j]] = [1..order] and
                ForAll( [1..order], x -> af[f[x]][j] = x ) ) );

  if i <= order then  
    
    # the case if the sg has an identity
    
    for a in autos do
      gens := [];
      for i in [1..order] do
        rho := [];
        for j in [1..order] do
          rho[a[j]] := a[ af[f[j]][i] ];
        od;
        Add( gens, Transformation( [1..order], rho ) );
      od;
      AddSet( sgps, TransformationSemigroup( gens ) );  
    od;

  else
  
    # the case if the sg has no identity
    
    for a in autos do
      gens := [];
      for i in [1..order] do
        rho := [];
        for j in [1..order] do
          rho[a[j]] := a[ af[f[j]][i] ];
        od;
        Add( rho, a[i] );
        Add( gens, Transformation( [1..order+1], rho ) );
      od;
      AddSet( sgps, TransformationSemigroup( gens ) );  
    od;
  
  fi;

  return sgps;
end;

