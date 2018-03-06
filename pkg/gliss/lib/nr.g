###############################################################################
##
##  nr.g  GLISSANDO ver 1.0  Christof Noebauer   1995, 1996
##                                                       
##
#############################################################################
##  Provide some GAP library functions for near-rings
#############################################################################
##
#############################################################################
##
#V  NearringOps . . . . . . . . . . . . . . . operations record for nearrings
##
##  'NearringOps' is the operation record for nearrings. This is initially
##  a copy of 'DomainOps'. This way all the default methods for domains are
##  inherited.
##
NearringOps := Copy( DomainOps );

#############################################################################
##
##  IsNearring( <obj> ) . . . . . . . . . . . . . test if <obj> is a nearring
##
IsNearring := function( obj )
  return IsRec( obj ) and IsBound( obj.isNearring ) and obj.isNearring; 
end;

#############################################################################
##
##  IsTransformationNearring( <obj> ). . . . . test if <obj> is a tf nearring
##
IsTransformationNearring := function( obj )
  return IsRec( obj ) and 
    IsBound( obj.isTransformationNearring ) and obj.isTransformationNearring; 
end;

#############################################################################
##
#F  Endomorphisms( <D> ). . . . . . . . . . create all endo's on a domain <D> 
##  V1.0  20.2.95
##  Dispatcher function for computing all endomorphisms on <D>.
##  Works for <D> = a group or a nearring.
##  JM 15-3-2016 moved to lib/dispatch.g since conflicted with chevie
##

#############################################################################
##
#F  Automorphisms( <D> ). . . . . . . . . . create all auto's on a domain <D> 
##  V1.0  20.2.95
##  Dispatcher function for computing all automorphisms on <D>.
##  Works for <D> = a group or a nearring.
##  JM 15-3-2016 moved to lib/dispatch.g since conflicted with sisyphos
##
#############################################################################  
##
##  IsNrMultiplication( <G>, <mul> ). . . . . . . test if <mul> is a nearring
##                                                multiplication on <G>
##
##  This function tests if a specified multiplication function <mul> is a 
##  nearring multiplication on the group <G>.
##
IsNrMultiplication := function( G, mul )
  local elms;
  elms := Elements( G );
  
  # check that the arguments are really a group and a WELL-DEFINED function
  if not ( IsGroup( G ) and IsFunc( mul ) and 
           ForAll( elms, n1 -> ForAll( elms, n2 -> mul( n1, n2 ) in G ) )
         ) then
    Error( "Usage: IsNrMultiplication( <G>, <mul> ) where <G> must be a ", 
           "group\nand <mul> must be a function <G> x <G> -> <G>" );
  fi;
  
  # check if mul is ASSOCIATIVE
  if not ForAll( elms, n1 -> ForAll( elms, n2 -> ForAll( elms, n3 ->
    mul( n1, mul( n2, n3 ) ) = mul( mul( n1, n2 ), n3 ) ) ) ) 
  then
    Print( "specified multiplication is not associative.\n" );
    return false;
  fi;

  # check if mul is RIGHT DISTRIBUTIVE 
  # (note that the addition is denoted by '*' )
  if not ForAll( elms, n1 -> ForAll( elms, n2 -> ForAll( elms, n3 ->
    mul( n1 * n2 , n3 ) =  mul( n1, n3 ) * mul( n2, n3 )  ) ) )
  then
    Print( "specified multiplication is not right distributive.\n" );
    return false;
  fi;

  return true;
end;

#############################################################################
##
#F  Nearring( <arg> ) . . . . . . . . . . . . . . . . . . . create a nearring
##
##  Constructor function for a nearring
##  So far there are two possibilities to construct a nearring:
##  1.) enter a group and a nearring multiplication on this group.
##  2.) enter a few group transformations and consider the generated
##      nearring. ( in this case Elements(.) is the most important function )
##
Nearring := function( arg )
  local G,    # the additive group of a nearring to be defined
        mul,  # a multiplication which makes G into a nearring
        gens, # generators of a nearring
        gt,   # a group transformation
        NR;   # the nearring to be returned
  
  arg := Flat( arg );
  if  Length( arg ) in [ 2, 3 ]  and IsGroup( arg[1] ) then 
    
    if Length( arg ) = 2 then 
      if not IsNrMultiplication( arg[1], arg[2] ) then return; fi; 
    fi;
    G   := arg[1];
    mul := arg[2];
    NR  := rec();
    # enter category components
    NR.isDomain       := true;
    NR.isNearring     := true; 
    # enter identification components  
    NR.group          := G;
    NR.addition       := function( x, y ) return x * y;    end;
    NR.subtraction    := function( x, y ) return x * y^-1; end;
    NR.multiplication := mul;
    # enter operations record
    NR.operations     := NearringOps;
    
  elif ForAll( arg, gt -> IsGroupTransformation( gt ) and    
               gt.source = arg[1].source and gt.range = arg[1].range
             ) then
    # make sure that all sources and ranges of the generators are not only
    # equal but identical.
    gens := [];
    for gt in arg do
      gt.source := arg[1].source; gt.range := arg[1].source;
      AddSet( gens, gt );
    od;
    NR  := rec();
    # enter category components
    NR.isDomain                 := true;
    NR.isNearring               := true; 
    NR.isTransformationNearring := true;
    # enter identification components  
    NR.generators               := gens;
    NR.group                    := "?";
    NR.addition                 := function( x, y ) return x + y; end;
    NR.subtraction              := function( x, y ) return x - y; end;
    NR.multiplication           := function( x, y ) return y * x; end; # N.B.!
    # enter operations record
    NR.operations     := NearringOps;

  else
    Error( "Usage: Nearring( <G>, <mul> ) where <G> must be a group\n", 
           "and <mul> must be a valid function <G> x <G> -> <G> or",
           "\nNearring( <t1>, <t2>, ... ) where all arguments must ",
           "be\ntransformations on the same group" );
  fi;
  return NR;
end;

#############################################################################
##
#F  NearringOps.Elements( <N> ). . . compute the elements of the nearring <N>
##
NearringOps.Elements := function( N )
  local closure,        # the constructed closure
        elms,           # the elements of the group on which the tf's operate
        l,              # the number of the elements
        firstrun,       # help var: indicates if first loop run
        tfl,tfl1,tfl2,  # help vars: transformation lists
        elmset,         # help var: contains the elms curr. in the closure 
        changed_elmset, # help var: indicates if elmset has changed
        done;           # help var: indicates, when it is time to stop
        
  if IsTransformationNearring( N ) then
    
    # perform a simple closure algorithm
    closure  := Set( List( N.generators, g -> g.tfl ) );
    elms     := N.generators[1].elements;
    l        := Length( elms );
    firstrun := true;
#    Print( closure, "\n" );

    repeat
      # step 1: add all sums  
#      Print( "Step 1: Building sums of full transformations\n" );
      elmset         := closure;
      changed_elmset := true;
      done           := true;
      while changed_elmset do
        closure        := Copy( elmset );
        changed_elmset := false;
        for tfl1 in closure do
          for tfl2 in closure do  
            tfl := List( [1..l], i -> 
                     Position( elms, elms[ tfl1[i] ] * elms[ tfl2[i] ] ) );
            if not tfl in elmset then 
              AddSet( elmset, tfl ); 
#              Print( "Adding", tfl, "\n" ); 
              changed_elmset := true; 
              done           := false;
            fi;
          od;
        od;
      od; 
      if not done or firstrun then 
        # step 2 build the multiplicative closure    
#        Print( "Step 2: Multiplying full transformations\n" );
        firstrun       := false;
        elmset         := closure;
        changed_elmset := true;
        done           := true;
        while changed_elmset do
          closure        := Copy( elmset );
          changed_elmset := false;
          for tfl1 in closure do
            for tfl2 in closure do
              tfl := tfl1{ tfl2 };
              if not tfl in elmset then 
                AddSet( elmset, tfl ); 
#                Print( "Adding", tfl, "\n" ); 
                changed_elmset := true; 
                done           := false;
              fi;
            od;
          od;
#          Print( "size:", Length( elmset ), "\n" );
        od;
      fi;
    until done;
    N.elements := [];
    for tfl in closure do
      AddSet( N.elements, Transformation( N.generators[1].source, tfl ) );
    od;
  
  else
    
    N.elements := Elements( N.group );
  
  fi;
  
  return N.elements;
end;

#############################################################################
##
#F  NearringOps.Print( <N> ) . . . . . . . . . . . . . . . . print a nearring
##
NearringOps.Print := function( N )
  local i;
  if IsTransformationNearring( N ) then
    Print( "Nearring( " );
    for i in [ 1..Length( N.generators ) - 1 ] do 
      Print( N.generators[i], ", " ); 
    od;
    Print( N.generators[ Length( N.generators )], " ) " );
  elif
    IsBound( N.isLibraryNearring ) and N.isLibraryNearring then
      Print( "LibraryNearring( \"", N.group, "\", ", N.idNum, " )" );
  else
    Print( "Nearring( ", N.group, ", ", N.multiplication, " )" );
  fi;
end;

#############################################################################
##
#F  FindGroup( <N> ). . . determine the additive group of a transformation nr 
##                        <N> as a GAP permutation group.
##  
##  Find the embedding which maps the additive group of transformations of
##  a transformation nearring <N> into the symmetric group Sn of all 
##  permutations on {1..n} where n is the size of the nearring.
##  The record fields 'group' and 'PHI' will be added to N.
##  The list PHI is a list of the permutations of the subgroup of Sn
##  such that there is a 1-1 correspondence between the list 
##  elms := Elements( N ) ( which is identical to N.elements ) and the 
##  list PHI s.t. PHI[1] = elms[1], PHI[2] = elms[2], ... ,PHI[n] = elms[n]. 
##
FindGroup := function( N )
  local TFL,   # the list of the transformation lists of the tf's of N
        telms, # the elms of the group a tf of N works on ( "transf'd elms" )
        s,     # the size of the group a tf in N works on
        PHI,   # the list of all permutations, return value
        tfl1,  # a transformation list of a transformation in N
        phi;   # the permutation derived for one fixed transformation in N

  if not IsTransformationNearring( N ) then
    Error( "Usage: FindGroup( <N> ) where <N> must be a transformation ", 
           "nearring" );
  fi;
  if not ( IsBound( N.group ) and IsBound( N.PHI ) ) then
    
    TFL   := List( Elements( N ), e -> e.tfl );
    telms := Elements( N.generators[1].source );
    s     := Size( telms );
    PHI   := [];
  
    for tfl1 in TFL do
      phi := PermList( List( TFL, tfl2 -> 
              Position( TFL, List( [1..s], i -> 
                Position( telms, telms[ tfl1[i] ] * telms[ tfl2[i] ] ) 
                     ))    )     );
      Add( PHI, phi );
    od;
    
    N.group := Group( SmallestGeneratingSystem( Group( PHI, () ) ), () );
    N.PHI   := PHI;
  fi;
  return N.group;
end;

#############################################################################
##
#F  NearringOps.Identity( <N> ) . . . . . . . . . . . compute identity of <N>
##
NearringOps.Identity := function( N )
  local mul, elms;
  
  mul  := N.multiplication;
  elms := Elements( N );
  
  return Filtered( elms, i -> ForAll( elms, e -> mul( i, e ) = e and 
                                                 mul( e, i ) = e ) );
end;

#############################################################################
##
#F  NearringOps.IsCommutative( <N> ). . . . . . . . test commutativity of <N>
##
NearringOps.IsCommutative := function( N )
  local mul, elms;
  
  mul  := N.multiplication;
  elms := Elements( N );
  
  return ForAll( elms, c -> ForAll( elms, e -> mul( c, e ) = mul( e, c ) ) );
end;

#############################################################################
##
#F  NearringOps.Endomorphisms( <N> ). . . create all endo's on a nearring <N> 
##
NearringOps.Endomorphisms := function( N )
  local elms,endos;
  if IsTransformationNearring( N ) then
    FindGroup( N );
  fi;
  elms  := Elements( N );
  endos := Filtered( Endomorphisms( N.group ), e ->
    ForAll( elms, x -> ForAll( elms, y ->
      elms[ e.tfl[ Position( elms, N.multiplication( x, y ) ) ]  ] =
      N.multiplication( elms[ e.tfl[ Position( elms, x ) ] ], 
                        elms[ e.tfl[ Position( elms, y ) ] ] ) ) ) );
  return endos;
end;

#############################################################################
##
#F  NearringOps.Automorphisms( <N> ). . . create all auto's on a nearring <N> 
##
NearringOps.Automorphisms := function( N )
  local elms,autos;
  if IsTransformationNearring( N ) then
    FindGroup( N );
  fi;
  elms  := Elements( N );
  autos := Filtered( Automorphisms( N.group ), e ->
    ForAll( elms, x -> ForAll( elms, y ->
      elms[ e.tfl[ Position( elms, N.multiplication( x, y ) ) ]  ] =
      N.multiplication( elms[ e.tfl[ Position( elms, x ) ] ], 
                        elms[ e.tfl[ Position( elms, y ) ] ] ) ) ) );
  return autos;
end;

#############################################################################
##
#F  NearringOps.DisplayTable( <arg> )  . . . . . print a Cayley table of a nr 
##
NearringOps.DisplayTable := function( arg )

  local N,       # the nearring
        elms,    # the elements of the nearring
        n,       # the size of the nearring
        symbols, # a list of the symbols for the elements of the nearring
        tw,      # the width of a table
        spc,     # local function which prints the right number of spaces
        bar,     # local function for printing the right length of the bar
        ind,     # help variable, an index
        print_addition, print_multiplication, # status variables
        i,j;     # loop variables
  
  print_addition := true; print_multiplication := true;
  if Length( arg ) = 2 then
    if arg[ 2 ] = "a" or arg[ 2 ] = "A" then
      print_multiplication := false;
    else
      print_addition := false;
    fi;
  fi;
  N       := arg[ 1 ];
  elms    := Elements( N );
  n       := Size( elms );
  symbols := [ AbstractGenerator( "n0" ) ];
  if n > 1 then
    symbols := Concatenation( [ AbstractGenerator( "n0" ) ], 
                              AbstractGenerators( "n" , n-1 ) );
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
#   Print( "The table of a nearring of order ", n, " will not look\n",
#          "good on a screen with ", SizeScreen()[1], " characters per ", 
#          "line.\nHowever, you may want to set your line length to a ",
#          "greater\nvalue by using the GAP function 'SizeScreen'.\n" );
#   return;
# fi;
  if IsTransformationNearring( N ) then FindGroup( N ); fi;
  if print_addition and print_multiplication then
    Print( "Let:\n" );
    for i in [1..n] do 
      Print( symbols[i], " := ", elms[i] ); 
      if IsTransformationNearring( N ) and TRANSFORMATION_PRINT_LEVEL = 0 then 
        Print( " ( <-> ", Elements( N.group )[i], " )" ); 
      fi;
      Print( "\n" );
    od;
  fi;
  
  if print_addition then
    # print the addition table
    Print( "\n  + ", spc( 0, n ), "| " ); 
      for i in [1..n] do Print( symbols[i], spc( i, n ) ); od;
    Print( "\n ---", bar( n ) ); for i in [1..n] do Print( bar( n ) ); od;
    for i in [1..n] do
      Print( "\n  ", symbols[i], spc( i, n ), "| " );
      for j in [1..n] do
        ind := Position( elms, N.addition( elms[i], elms[j] ) );
        Print( symbols[ ind ], spc( ind, n ) );  
      od;
    od;
  fi;
  
  if print_multiplication then
    # print the multiplication table
    Print( "\n\n  * ", spc( 0, n ), "| " ); 
      for i in [1..n] do Print( symbols[i], spc( i, n ) ); od;
    Print( "\n ---", bar( n ) ); for i in [1..n] do Print( bar( n ) ); od;
    for i in [1..n] do
      Print( "\n  ", symbols[i], spc( i, n ), "| " );
      for j in [1..n] do
        ind := Position( elms, N.multiplication( elms[i], elms[j] ) );
        Print( symbols[ ind ], spc( ind, n ) );  
      od;
    od;
  fi;
  
  Print( "\n\n" );
end;

#############################################################################
##
##  NearringOps.IsNrIdeal( <N>, <I> ) . . . . check if <I> is an ideal of <N>
##
NearringOps.IsNrIdeal := function( N, I )

  local isid,  # return value: a record with the boolean record fields: 
               #               isLeftIdeal, isRightIdeal, isIdeal.
        ielms; # help var: the elements of  the subgroup being considered

  if not IsNearring( N ) then
    Error( N, " must be a nearring" );
  fi;
  if IsTransformationNearring( N ) then
    FindGroup( N );
  fi;
  if not IsSubgroup( N.group, I ) then
    Error( I, " must be a subgroup of ", N, ".group" );
  fi;
  if not IsNormal( N.group, I ) then
    Error( I, " must be a normal subgroup of ", N, ".group" );
  fi;

  # in case of a transformation nearring take those elements 
  # (=transformations) of N which form the subgroup
  if IsTransformationNearring( N ) then
    ielms := Filtered( Elements( N ), n -> 
      Elements( N.group )[ Position( Elements( N ), n ) ] in Elements( I ) );
  else
    ielms := Elements( I );
  fi;
#  Print( "ielms: ", ielms, "\n" );
  
  isid := rec(
    isLeftIdeal := false,
    isRightIdeal := false,
    isIdeal := false );

  if ForAll( Elements( N ), n -> ForAll( ielms, i ->
             N.multiplication( i, n ) in ielms ) )
  then
    isid.isRightIdeal := true;
  fi;

  if ForAll( Elements( N ), m -> ForAll( Elements( N ), n ->
     ForAll( ielms, i ->
             N.subtraction(
             N.multiplication( m, N.addition( n, i) ),
             N.multiplication( m, n ) ) in ielms ) ) )
  then
    isid.isLeftIdeal := true;
  fi;

  if isid.isLeftIdeal and isid.isRightIdeal then
    isid.isIdeal := true;
  fi;

#  Print("I =   ",I,"\n");
#  Print("isid= ",isid,"\n");
  
  return isid;
end;

#############################################################################
##
##  NearringIdeals( <arg> ). . . . . . . . . . . . compute all ideals of a nr
##
NearringIdeals := function( arg )

  local N,          # the nearring
        L,          # the lattice of subgroups of the add. group 
        Rep,        # representative of the <i>-th class
        normalizer, # normalizer of <I> in <N.group>
        reps,       # transversal of <normalizer> in <N.group>
        I,          # the subgroup being considered
        elms_I,     # the elements of the subgroup M
        elms_N,     # the elements of the (group of the) nearring
        ideals,     # the list of ideals
        add,
        sub,
        mul,        # the nearring multiplication
        right_only,
        left_only,
        i,k;        # loop variables
  
  if not ( 
           Length( arg ) = 1 and IsNearring( arg[1] ) or ( 
           Length( arg ) = 2 and
           ( arg[2] = "l" or arg[2] = "L" or arg[2] = "r" or arg[2] = "R" ) ) 
         ) 
           then
    Error( "Usage: NearringIdeals( <N> ) or NearringIdeals( <N>, <\"r\"> )\n",
           "or NearringIdeals( <N>, <\"l\"> ) where <N> must be a nearring" );
  fi;

  N := arg[ 1 ]; right_only := false; left_only := false;
  if Length( arg ) = 2 then 
    if arg[2] = "r" or arg[2] = "R" then
      right_only := true;
    else
      left_only := true;
    fi;
  fi;

  if right_only then
    if IsBound( N.rightIdeals ) then return N.rightIdeals; fi;
  elif left_only then
    if IsBound( N.leftIdeals ) then return N.leftIdeals; fi;
  else
    if IsBound( N.ideals ) then return N.ideals; fi;
  fi;

  if IsTransformationNearring( N ) then FindGroup( N ); fi;
  
  add    := N.addition;
  sub    := N.subtraction;
  mul    := N.multiplication;
  ideals := [ Subgroup( N.group, [ ] ) ];
  elms_N := Elements( N );
  L      := Lattice( N.group );
  
  for i in [ 2..Length( L.classes )-1 ] do
    
    Rep := L.classes[ i ].representative;
    # get the transversal
    normalizer := Normalizer( L.group, ShallowCopyNoSC(Rep) );
    reps := RightTransversal( L.group, ShallowCopyNoSC(normalizer) );
    
    # consider all normal subgroups of N.group
#   for k in [ 1..Length( reps ) ] do
    if Length( reps ) = 1 then  
#     I := Rep^reps[ k ];
      I := Rep;
      if IsNormal( N.group, I ) then
      
        elms_I := List( Elements( I ), e -> 
                  elms_N[ Position( Elements( N.group ), e ) ] );

        # this is the check for the (right) (left) ideal condition
        if right_only then
          if ForAll( elms_I, i -> ForAll( elms_N, n ->
             mul( i, n ) in elms_I ) ) then 
            Add( ideals, I );
          fi;
        elif left_only then  
          if ForAll( elms_N, n -> ForAll( elms_N, m -> ForAll( elms_I, i ->
             sub( mul( n, add( m, i ) ), mul( n, m ) ) in elms_I ) ) ) then  
            Add( ideals, I );
          fi;
        else
          if ForAll( elms_I, i -> ForAll( elms_N, n ->
                     mul( i, n ) in elms_I ) ) and 
            ForAll( elms_N, n -> ForAll( elms_N, m -> ForAll( elms_I, i ->
            sub( mul( n, add( m, i ) ), mul( n, m ) ) in elms_I ) ) ) then  
            Add( ideals, I );
          fi;
        fi;
      fi;
    
#   od;
    fi;

  od;
  
  Add( ideals, N.group );

  if right_only then
    N.rightIdeals := ideals;
  elif left_only then
    N.leftIdeals := ideals;
  else
    N.ideals := ideals;
  fi;
 
  return ideals;
end;

#############################################################################
##
##  InvariantSubnearrings( <N> ). . . . . . . compute all inv subnr's of <N>.
##
InvariantSubnearrings := function( N )

  local L,          # the lattice of subgroups of the add. group 
        Rep,        # representative of the <i>-th class
        normalizer, # normalizer of <I> in <N.group>
        reps,       # transversal of <normalizer> in <N.group>
        M,          # the subgroup being considered
        elms_M,     # the elements of the subgroup M
        elms_N,     # the elements of the (group of the) nearring
        inv_sub_nrs,# the list of invariant subnearrings
        mul,        # the nearring multiplication
        i,k;        # loop variables
  
  if not IsNearring( N ) then
    Error( "Usage: InvariantSubnearrings( <N> ) where <N> must be a ",
           "nearring" );
  fi;

  if IsTransformationNearring( N ) then FindGroup( N ); fi;
  
  mul         := N.multiplication;
  inv_sub_nrs := [];
  elms_N      := Elements( N );
  L           := Lattice( N.group );
  
  for i in [1..Length(L.classes)] do
    
    Rep := L.classes[i].representative;
    # get the transversal
    normalizer := Normalizer( L.group, ShallowCopyNoSC(Rep) );
    reps := RightTransversal( L.group, ShallowCopyNoSC(normalizer) );
    
    # consider all subgroups of N.group
    for k in [1..Length(reps)] do
      
      M      := Rep^reps[ k ];
      elms_M := List( Elements( M ), e -> 
                elms_N[ Position( Elements( N.group ), e ) ] );

      # this is the check for the invariant subnr condition
      if ForAll( elms_N, n -> ForAll( elms_M, m ->
        mul( m, n ) in elms_M and mul( n, m ) in elms_M ) ) then
        Add( inv_sub_nrs, M );
      fi;
    
    od;
  
  od;
  
  if IsTransformationNearring( N ) then
    return List( inv_sub_nrs, i -> Nearring( List( Elements( i ), e ->
      Elements( N )[ Position( Elements( N.group ), e ) ] ) ) );
  else
    return List( inv_sub_nrs, i -> Nearring( i, mul, "n" ) );
  fi;
end;

#############################################################################
##
##  Subnearrings( <N> ) . . . . . . . . . . . . . compute all subnr's of <N>.
##
Subnearrings := function( N )

  local L,          # the lattice of subgroups of the add. group 
        Rep,        # representative of the <i>-th class
        normalizer, # normalizer of <I> in <N.group>
        reps,       # transversal of <normalizer> in <N.group>
        M,          # the subgroup being considered
        elms_M,     # the elements of the subgroup M
        sub_nrs,    # the list of subnearrings
        mul,        # the nearring multiplication
        i,k;        # loop variables
  
  if not IsNearring( N ) then
    Error( "Usage: Subnearrings( <N> ) where <N> must be a nearring" );
  fi;

  if IsTransformationNearring( N ) then FindGroup( N ); fi;
  
  mul     := N.multiplication;
  sub_nrs := [];
  L       := Lattice( N.group );
  
  for i in [1..Length(L.classes)] do
    
    Rep := L.classes[i].representative;
    # get the transversal
    normalizer := Normalizer( L.group, ShallowCopyNoSC(Rep) );
    reps := RightTransversal( L.group, ShallowCopyNoSC(normalizer) );
    
    # consider all subgroups of N.group
    for k in [1..Length(reps)] do
      
      M      := Rep^reps[ k ];
      elms_M := List( Elements( M ), e -> 
                Elements( N )[ Position( Elements( N.group ), e ) ] );

      # this is the check for the subnr condition
      if ForAll( elms_M, n -> ForAll( elms_M, m -> mul( m, n ) in elms_M ) ) 
        then Add( sub_nrs, M );
      fi;
    
    od;
  
  od;
  
  if IsTransformationNearring( N ) then
    return List( sub_nrs, i -> Nearring( List( Elements( i ), e ->
      Elements( N )[ Position( Elements( N.group ), e ) ] ) ) );
  else
    return List( sub_nrs, i -> Nearring( i, mul, "n" ) );
  fi;
end;

#############################################################################
##
#F  LibraryNearring( <name>, <num> ). . . . . . . . get a nr from the library
##
##  This function 'extracts' a nearring from the nearring library files.
##
LibraryNearring := function( name, num )

  local n,     # help var: a nearring
        clmax, # the maximal number of equivalence classes of nearrings
        G,     # the additive group of the nr to be returned
        NR,    # the nearring to be returned
        elms,  # help var: the elements of G
        i,     # help var: a loop variable
        tfle,  # help var: the record that holds the tfl's of the group endos
        f,     # help var: a valid function that represents a class of nr's
        vf,endos,g,a,a_inv,h,compute_all,
        mul;   # local function: the multiplication of the nearring

  # check the arguments 
  if not ( IsString( name ) and IsInt( num ) and num > 0 ) then
    Error( "Usage: LibraryNearring( <name>, <num> ) where <name> must be ", 
           "the\n       name of a group and <num> must be a positive ",
           "integer which\ndetermines an isomorphism class" );
  fi;
  
  if   ( name = "C2" ) then  
    n := NR_C2; clmax := 3;
  elif ( name = "C3" ) then
    n := NR_C3; clmax := 5;
  elif ( name = "C4" ) then
    n := NR_C4; clmax := 12;
  elif ( name = "V4" ) then
    n := NR_V4; clmax := 23;
  elif ( name = "C5" ) then
    n := NR_C5; clmax := 10;
  elif ( name = "C6" ) then
    n := NR_C6; clmax := 60;
  elif ( name = "S3" ) then
    n := NR_S3; clmax := 39;
  elif ( name = "C7" ) then
    n := NR_C7; clmax := 24; 
  elif ( name = "C8" ) then
    n := NR_C8; clmax := 135;
  elif ( name = "C2xC4" ) then
    n := NR_C2xC4; clmax := 1159;
  elif ( name = "C2xC2xC2" ) then
    n := NR_C2xC2xC2; clmax := 834;
  elif ( name = "D8" ) then
    n := NR_D8; clmax := 1447;
  elif ( name = "Q8" ) then
    n := NR_Q8; clmax := 281;
  elif ( name = "C9" ) then
    n := NR_C9; clmax := 222;
  elif ( name = "C3xC3" ) then
    n := NR_C3xC3; clmax := 264;
  elif ( name = "C10" ) then
    n := NR_C10; clmax := 329;
  elif ( name = "D10" ) then
    n := NR_D10; clmax := 206;
  elif ( name = "C11" ) then
    n := NR_C11; clmax := 139;
  elif ( name = "C12" ) then
    n := NR_C12; clmax := 1749;
  elif ( name = "C2xC6" ) then
    n := NR_C2xC6; clmax := 3501;
  elif ( name = "D12" ) then
    
    if   num in [1..5000] then
      n := NR_D12_1;
    elif num in [5001..10000] then
      n := NR_D12_2;
    elif num in [10001..15000] then
      n := NR_D12_3;
    elif num in [15001..20000] then
      n := NR_D12_4;
    elif num in [20001..25000] then
      n := NR_D12_5;
    elif num in [25001..30000] then
      n := NR_D12_6;
    elif num in [30001..35000] then
      n := NR_D12_7;
    elif num in [35001..40000] then
      n := NR_D12_8;
    elif num in [40001..45000] then
      n := NR_D12_9;
    else 
      n := NR_D12_10;
    fi; clmax := 48137;

  elif ( name = "A4" ) then
    n := NR_A4; clmax := 483;
  elif ( name = "T" ) then
    n := NR_T; clmax := 824;
  elif ( name = "C13" ) then
    n := NR_C13; clmax := 454;
  elif ( name = "C14" ) then
    n := NR_C14; clmax := 2716;
  elif ( name = "D14" ) then
    n := NR_D14; clmax := 1821;
  elif ( name = "C15" ) then
    n := NR_C15; clmax := 3817;
  else
    Print( "There is no group name '", name, 
           "' in the nearrings library.\n" );
    return;
  fi;
  
  if num > clmax then
    Print( "There are only ", clmax, " isomorphism classes of nearrings ",
           "on the group ", name, ".\n" );
    return;
  fi;
    
  # put the group of the nearring together and define a few help variables
  G      := Group( n.group_generators, () );
  G.name := n.group_name;
  elms   := Elements( G );
  tfle   := n.group_endomorphisms;
  f      := n.classes.(num).phi;
  G.phi  := f;
  G.a_y_i_nrs := n.classes.(num).autos_yielding_iso_nrs;
  # retrieve the group endomorphisms from the Nearrings record
  if not IsBound( G.endomorphisms ) then
    i := 1; G.endomorphisms := [];         # convert the endomorphism record
    while IsBound( tfle.(i) ) do           # into a list of endomorphisms
      Add( G.endomorphisms, Transformation( G, tfle.(i) ) );   
      i := i + 1;
    od;
  fi;
  
  compute_all := false;
  vf := [ f ]; 
  if compute_all then
    endos := []; g:= [];
    for i in [1..Length( RecFields( tfle ) )] do
      Add( endos, tfle.(i) );
    od;
    for a in n.classes.(num).autos_yielding_iso_nrs do
      a := endos[a]; a_inv := [];
      for i in [1..Length( a )] do a_inv[a[i]] := i; od;
      for i in [1..Length( a )] do
        h := a{endos[f[i]]};
        g[ a[i] ] := Position( endos, h{a_inv} );
      od;
      AddSet( vf, Copy( g ) );
    od;
  fi;  

  # define a RIGHT distributive multiplication  
  mul := function( x, y )
    return elms[ tfle.(f[ Position( elms, y ) ]) [ Position( elms, x ) ] ];
  end;
# define a LEFT distributive multiplication  
#  mul := function( x, y )
#    return elms[ tfle.(f[ Position( elms, x ) ]) [ Position( elms, y ) ] ];
#  end;
  
  # put the nearring together
  NR  := rec();
  # enter category components
  NR.isDomain       := true;
  NR.isNearring     := true; 
  NR.isLibraryNearring := true;
  # enter identification components  
  NR.group          := G;
  NR.idNum          := num;
  NR.addition       := function( x, y ) return x * y;    end;
  NR.subtraction    := function( x, y ) return x * y^-1; end;
  NR.multiplication := mul;
  # enter operations record
  NR.operations     := NearringOps;
  
  if compute_all then
    NR.vf_of_iso_nrs := vf;
  fi;
  
  return NR;
end;

#############################################################################
##
#F  Distributors( <N> ) . . compute the set of distributors on a nearring <N> 
##
##  Dispatcher function for computing the distributors of <N>.
##
Distributors := function( N )
  if IsNearring( N ) then
    if not IsBound( N.distributors ) then
      N.distributors := N.operations.Distributors( N ); 
    fi;
  else
    Error( "Usage: Distributors( <N> ) where <N> must be a nearring" );
  fi;

  return N.distributors;
end;
 
#############################################################################
##
#F  DistributiveElements( <N> ). . . compute the distributive elements of <N> 
##
##  Dispatcher function for computing the distributive elements of <N>.
##
DistributiveElements := function( N )
  if IsNearring( N ) then
    if not IsBound( N.distributiveElements ) then
      N.distributiveElements := N.operations.DistributiveElements( N ); 
    fi;
  else
    Error( "Usage: DistributiveElements( <N> ) where <N> must be a nearring");
  fi;

  return N.distributiveElements;
end;
 
#############################################################################
##
#F  ZeroSymmetricElements( <N> ) . compute the zero-symmetric elements of <N> 
##
##  Dispatcher function for computing the zero-symmetric elements of <N>,
##  i.e. all elements n s.t. n0 = 0 (note: 0n = 0 is always true) 
##
ZeroSymmetricElements := function( N )
  if IsNearring( N ) then
    if not IsBound( N.zeroSymmetricElements ) then
      N.zeroSymmetricElements := N.operations.ZeroSymmetricElements( N ); 
    fi;
  else
    Error("Usage: ZeroSymmetricElements( <N> ) where <N> must be a nearring");
  fi;

  return N.zeroSymmetricElements;
end;
 
#############################################################################
##
#F  IdempotentElements( <D> ). . . . . compute the idempotent elements of <D> 
##
##  Dispatcher function for computing the idempotent elements of <D>.
##
IdempotentElements := function( D )
  if IsNearring( D ) or IsSemigroup( D ) then
    if not IsBound( D.idempotentElements ) then
      D.idempotentElements := D.operations.IdempotentElements( D ); 
    fi;
  else
    Error( "Usage: IdempotentElements( <D> ) where <D> must be a ", 
           "nearring or a semigroup" );
  fi;

  return D.idempotentElements;
end;
 
#############################################################################
##
#F  NilpotentElements( <N> ). . . . . . . . compute the nilpotent elms of <N>
##
##  Dispatcher function to compute the nilpotent elements of a nearring <N> 
##
NilpotentElements := function( N )
  
  if not IsNearring( N )  then
    Error( "Usage: NilpotentElements( <N> ) where <N> must be a nearring" );
  fi;

  if not IsBound( N.nilpotentElements ) then
    N.nilpotentElements := N.operations.NilpotentElements( N );
  fi;

  return N.nilpotentElements;
end;

#############################################################################
##
#F  QuasiregularElements( <N> ). . . . . compute the quasiregular elms of <N>
##
##  Dispatcher function to compute the quasiregular elements of a nr <N>. 
##
QuasiregularElements := function( N )
  
  if not IsNearring( N )  then
    Error( "Usage: QuasiregularElements( <N> ) where <N> must be a ",
           "nearring" );
  fi;

  if not IsBound( N.quasiregularElements ) then
    N.quasiregularElements := N.operations.QuasiregularElements( N );
  fi;

  return N.quasiregularElements;
end;

#############################################################################
##
#F  RegularElements( <N> ). . . . . . . . . . compute the regular elms of <N>
##
##  Dispatcher function to compute the regular elements of a nr <N>. 
##
RegularElements := function( N )
  
  if not IsNearring( N )  then
    Error( "Usage: RegularElements( <N> ) where <N> must be a nearring" );
  fi;

  if not IsBound( N.regularElements ) then
    N.regularElements := N.operations.RegularElements( N );
  fi;

  return N.regularElements;
end;

#############################################################################
##
#F  IsAbstractAffineNearring( <N> ) . . . . . . . . . . . test if <N> is a.a. 
##
IsAbstractAffineNearring := function( N )
  if IsNearring( N ) then
    if not IsBound( N.isAbstractAffine ) then
      if IsTransformationNearring( N ) then FindGroup( N ); fi;
      N.isAbstractAffine := IsAbelian( N.group ) and
                 ( ZeroSymmetricElements( N ) = DistributiveElements( N ) );
    fi;
  else
    Error( "Usage: IsAbstractAffineNearring( <N> ) where <N> must be a ", 
           "nearring" );
  fi;

  return N.isAbstractAffine;
end;

#############################################################################
##
#F  IsDistributiveNearring( <N> ) . . . . . . . . test if <N> is distributive 
##
IsDistributiveNearring := function( N )
  if IsNearring( N ) then
    if not IsBound( N.isDistributive ) then
      N.isDistributive := 
        Size( DistributiveElements( N ) ) = Size( Elements( N ) ); 
    fi;
  else
    Error( "Usage: IsDistributiveNearring( <N> ) where <N> must be a ",
           "nearring" );
  fi;

  return N.isDistributive;
end;

#############################################################################
##
#F  IsBooleanNearring( <N> ) . . . . . . . . . . . . . test if <N> is boolean 
##
IsBooleanNearring := function( N )
  if IsNearring( N ) then
    if not IsBound( N.isBoolean ) then
      N.isBoolean := 
        Size( IdempotentElements( N ) ) = Size( Elements( N ) ); 
    fi;
  else
    Error( "Usage: IsBooleanNearring( <N> ) where <N> must be a nearring" );
  fi;

  return N.isBoolean;
end;

#############################################################################
##
#F  IsDgNearring( <N> ) . . . . . . . test if <N> is distributively generated 
##
IsDgNearring := function( N )
  local Nd, elms;
  
  if IsNearring( N ) then
    if not IsBound( N.isDg ) then
      elms := Elements( N );
      if IsTransformationNearring( N ) then 
        FindGroup( N ); 
        Nd := List( DistributiveElements( N ), de -> 
                    Elements( N.group)[ Position( elms, de ) ] );
      else
        Nd := DistributiveElements( N );
      fi;
      N.isDg := Group( Nd, () ) = N.group;
    fi;
  else
    Error( "Usage: IsDgNearring( <N> ) where <N> must be a nearring" );
  fi;

  return N.isDg;
end;

#############################################################################
##
#F  IsIntegralNearring( <N> ) . . . . . . . . . . . . test if <N> is integral 
##
##  A nr is called integral if it has no zero divisors
##
IsIntegralNearring := function( N )
  local mul, elms, zero, non_zero_elms;
  
  if IsNearring( N ) then
    if not IsBound( N.isIntegral ) then
      mul  := N.multiplication;
      elms := Elements( N );
      zero := elms[ 1 ];  ## the first element is always the zero!
      non_zero_elms := Copy( elms ); RemoveSet( non_zero_elms, zero );
      N.isIntegral := ForAll( non_zero_elms, x -> 
                      ForAll( non_zero_elms, y -> mul( x, y ) <> zero ) );
    fi;
  else
    Error( "Usage: IsIntegralNearring( <N> ) where <N> must be a nearring" );
  fi;

  return N.isIntegral;
end;

#############################################################################
##
#F  IsNilNearring( <N> ) . . . . . . . . . . . . . . . . . test if <N> is nil 
##
IsNilNearring := function( N )
  if IsNearring( N ) then
    if not IsBound( N.isNil ) then
      N.isNil := 
        Length( NilpotentElements( N ) ) = Size( Elements( N ) ); 
    fi;
  else
    Error( "Usage: IsNilNearring( <N> ) where <N> must be a nearring" );
  fi;

  return N.isNil;
end;

#############################################################################
##
#F  IsNilpotentNearring( <N> ) . . . . . . . . . . . test if <N> is nilpotent 
##
IsNilpotentNearring := function( N )
  local mul, elms, prod, previous_prod, m, n;
  
  if IsNearring( N ) then
    if not IsBound( N.isNilpotent ) then
      if ( Size( IdempotentElements( N ) ) > 1 ) or 
         ( not IsNilNearring( N ) ) then
          N.isNilpotent := false;
      else
        
        mul := N.multiplication;
        elms := Elements( N );
        previous_prod := Copy( elms );
        
        repeat
          
          prod := []; 
          for m in previous_prod do
            for n in elms do
              AddSet( prod, mul( m, n ) );
            od;
          od;
          if prod = previous_prod then 
            N.isNilpotent := false;
          elif Size( prod ) = 1 then 
            N.isNilpotent := true;
          else 
            previous_prod := Copy( prod );
          fi;
        until IsBound( N.isNilpotent );

      fi;
    fi;
  else
    Error( "Usage: IsNilpotentNearring( <N> ) where <N> must be a nearring" );
  fi;

  return N.isNilpotent;
end;

#############################################################################
##
#F  IsPrimeNearring( <N> ) . . . . . . . . . . . . . . . test if <N> is prime 
##
##  A nearring N is called prime if { 0 } is a prime ideal
##
IsPrimeNearring := function( N )
  local mul, ideals, elms, elms_G, zero;
  
  if IsNearring( N ) then
    if not IsBound( N.isPrime ) then
      
      if IsIntegralNearring( N ) then
        N.isPrime := true;
      else
        ideals := Copy( NearringIdeals( N ) ); Unbind( ideals[ 1 ] );
        mul    := N.multiplication;
        if IsTransformationNearring( N ) then 
          FindGroup( N ); 
          elms   := Elements( N );
          elms_G := Elements( N.group );
          zero   := elms[ 1 ];
          N.isPrime := ForAll( ideals, I -> ForAll( ideals, J ->
                     ForAny( Elements( I ), i -> ForAny( Elements( J ), j ->
                     mul( elms[ Position( elms_G, i ) ], 
                          elms[ Position( elms_G, j ) ] ) <> zero  ) ) ) );
        else
          N.isPrime := ForAll( ideals, I -> ForAll( ideals, J ->
                       ForAny( Elements( I ), i -> ForAny( Elements( J ), j ->
                       mul( i, j ) <> ()  ) ) ) );
        fi;
      fi;
    
    fi;
  
  else
    Error( "Usage: IsPrimeNearring( <N> ) where <N> must be a nearring" );
  fi;

  return N.isPrime;
end;

#############################################################################
##
#F  IsQuasiregularNearring( <N> ) . . . . . . . . . . . . . test if <N> is qr 
##
IsQuasiregularNearring := function( N )
  if IsNearring( N ) then
    if not IsBound( N.isQuasiregular ) then
      N.isQuasiregular := 
        Size( QuasiregularElements( N ) ) = Size( Elements( N ) ); 
    fi;
  else
    Error( "Usage: IsQuasiregularNearring( <N> ) where <N> must be a ",
           "nearring" );
  fi;

  return N.isQuasiregular;
end;

#############################################################################
##
#F  IsRegularNearring( <N> ) . . . . . . . . . . . . . test if <N> is regular 
##
IsRegularNearring := function( N )
  if IsNearring( N ) then
    if not IsBound( N.isRegular ) then
      N.isRegular := 
        Size( RegularElements( N ) ) = Size( Elements( N ) ); 
    fi;
  else
    Error( "Usage: IsRegularNearring( <N> ) where <N> must be a nearring" );
  fi;

  return N.isRegular;
end;

#############################################################################
##
#F  IsNilpotentFreeNearring( <N> ) . . test if <N> is w/o non-zero nilpotents 
##
IsNilpotentFreeNearring := function( N )
  if IsNearring( N ) then
    if not IsBound( N.isNilpotentFree ) then
      N.isNilpotentFree := Length( NilpotentElements( N ) ) = 1; 
    fi;
  else
    Error( "Usage: IsNilpotentFreeNearring( <N> ) where <N> must be a ",
           "nearring" );
  fi;

  return N.isNilpotentFree;
end;

#############################################################################
##
#F  CouplingFunctions( <N> ). . compute the coupling functions <N> is made of 
##
##  this is a support function for the function IsPlanarNearring.
##                                              =-=-=-=-=-=-=-=-=
CouplingFunctions := function( N )

  local cf, elms, mul, a, f_a;

  if not IsNearring( N ) then 
    Print( "Usage: CouplingFunctions( <N> ) where <N> must be a nearring" );
  fi;

  if IsTransformationNearring( N ) then FindGroup( N ); fi;
  elms := Elements( N );
  mul  := N.multiplication;
  cf   := [ ];

  for a in elms do
    f_a := List( elms, b -> Position( elms, mul( b, a ) ) );
    Add( cf, Transformation( N.group, f_a ) );
  od;

  return cf;
end;

#############################################################################
##
#F  IsPlanarNearring( <N> ) . . . . . . . . . . . . . . test if <N> is planar 
##
IsPlanarNearring := function( N )
  local phi, size, endos, CF, zero, identity;

  if IsBound( N.isLibraryNearring ) then
    if not IsBound( N.isPlanar ) then
      phi   := Set( N.group.phi );
      size  := Size( N.group );
      endos := Endomorphisms( N.group );
      N.isPlanar := Size( phi ) >= 3 and
        ForAll( phi, p -> p = 1 or p = Length( endos ) or 
          ( Size( Set( endos[p].tfl ) ) = size and
            ForAll( [2..size], i -> endos[p].tfl[i] <> i ) ) );
    fi;
 
  elif IsNearring( N ) then
    if not IsBound( N.isPlanar ) then
      if IsTransformationNearring( N ) then FindGroup( N ); fi;
      CF       := Set( CouplingFunctions( N ) );
      size     := Size( N );
      zero     := Transformation( N.group, List( [1..size], i -> 1 ) );
      identity := IdentityTransformation( N.group );

      N.isPlanar := Size( CF ) >= 3 and 
        ForAll( CF, cf -> cf = zero or cf = identity or
                          ( 
                            Size( Set( cf.tfl ) ) = size  # cf is auto 
                            and ForAll( [2..size], i -> cf.tfl[i] <> i ) 
                          )
              );
    fi;

  else
    Error( "Usage: IsPlanarNearring( <N> ) where <N> must be a nearring" );
  fi;

  return N.isPlanar;
end;

#############################################################################
##
#F  IsNearfield( <N> ) . . . . . . . . . . . . . . test if <N> is a nearfield 
##
IsNearfield := function( N )
  local mul, id, elms;

  if IsNearring( N ) then
    if not IsBound( N.isNearfield ) then
      mul  := N.multiplication;
      elms := Elements( N );
      id   := Identity( N );
      N.isNearfield := id <> [ ] and
        Size( 
              Filtered( elms, e -> ForAny( elms, x -> mul( e, x ) = id[1] ) ) 
            ) = Size( elms ) - 1; 
    fi;
  else
    Error( "Usage: IsNearfield( <N> ) where <N> must be a nearring" );
  fi;

  return N.isNearfield;
end;

#############################################################################
##
#F  NearringOps.Distributors( <N> ). . compute distributors of a nearring <N> 
##
NearringOps.Distributors := function( N )
  local elms, a, b, c, add, sub, mul, dbs;

  add := N.addition;
  sub := N.subtraction;
  mul := N.multiplication;
  dbs := [];
  elms := Elements( N );
  for a in elms do
    for b in elms do
      for c in elms do
        AddSet( dbs, sub( mul( a, add( b, c ) ), 
                          add( mul( a, b ), mul( a, c ) ) ) ); 
      od;
    od;
  od;

  return dbs;
end;

#############################################################################
##
#F  NearringOps.DistributiveElements( <N> ). . . . . distributive elms of <N> 
##
##  Note: this function works only for RIGHT nearrings, i.e. it computes
##        the LEFT distributive elements.
##
NearringOps.DistributiveElements := function( N )
  local add, mul, elms;

  add  := N.addition;
  mul  := N.multiplication;
  elms := Elements( N );  

  return Filtered( elms, d -> ForAll( elms, a -> ForAll( elms, b ->
              mul( d, add( a, b ) ) = add( mul( d, a ), mul( d, b ) ) ) ) );

end;

#############################################################################
##
#F  NearringOps.ZeroSymmetricElements( <N> ) . . . zero-symmetric elms of <N> 
##
##  Note: this function works only for RIGHT nearrings, i.e. it computes
##        all elements n s.t n0 = 0. (Note that in a RIGHT nearring 
##        0n = 0 is always true).
##
NearringOps.ZeroSymmetricElements := function( N )
  local mul, elms, zero;

  mul  := N.multiplication;
  elms := Elements( N );  
  zero := elms[ 1 ];  # the first element is the zero element!
  
  return Filtered( elms, n -> mul( n, zero ) = zero );

end;

#############################################################################
##
#F  NearringOps.IdempotentElements( <N> ). . . . . . . idempotent elms of <N> 
##
NearringOps.IdempotentElements := function( N )
  local mul;

  mul  := N.multiplication;

  return Filtered( Elements( N ), i -> mul( i, i ) = i );

end;

#############################################################################
##
#F  NearringOps.NilpotentElements( <N> ). . compute nilpotent elements of <N>
##
NearringOps.NilpotentElements := function( N )
  local mul, elms, elm, e, size, zero, npelms, k, old_e;
  
  mul    := N.multiplication;
  elms   := Elements( N );
  size   := Size( elms );
  zero   := elms[ 1 ];   ## the first element is the zero!
  npelms := [ ];
  
  for elm in elms do
    k := 1; e := Copy( elm );
    old_e := zero;
    while e <> zero and e <> old_e and k < size do  
      k := k + 1;
      old_e := Copy( e );
      e := mul( e, elm );
    od;
    if e = zero then Add( npelms, [ elm, k ] ); fi;
  od;
  
  return npelms;
end;

#############################################################################
##
#F  NearringOps.QuasiregularElements( <N> ). . . . compute qr elements of <N>
##
NearringOps.QuasiregularElements := function( N )
  local mul, sub, elms, z, elms_to_test, qr_elms, A, li, Lz;
  
  mul     := N.multiplication;
  sub     := N.subtraction;
  elms    := Elements( N );
  qr_elms := List( NilpotentElements( N ), n -> n[1] );
  elms_to_test := Copy( elms );
  SubtractSet( elms_to_test, IdempotentElements( N ) );
  SubtractSet( elms_to_test, qr_elms );
  
  for z in elms_to_test do
    A := Set( List( elms, n -> sub( n, mul( n, z ) ) ) );   
    if IsTransformationNearring( N ) then
      FindGroup( N );
      li := List( NearringIdeals( N, "l" ), i -> 
        List( Elements(i), e -> elms[ Position( Elements(N.group), e ) ] ) );
    else
      li := List( NearringIdeals( N, "l" ), i -> Elements( i ) );
    fi;
  
    Lz := First( li, i -> IsSubset( i, A ) );
    if z in Lz then AddSet( qr_elms, z ); fi;
  od;
  
  return qr_elms;
end;

#############################################################################
##
#F  NearringOps.RegularElements( <N> ). . . . . . . . . . regular elms of <N> 
##
NearringOps.RegularElements := function( N )
  local mul, elms;

  mul  := N.multiplication;
  elms := Elements( N );

  return 
    Filtered( elms, x -> ForAny( elms, y -> mul( x, mul( y, x ) ) = x ) );

end;

#############################################################################
##
#F  LibraryNearringInfo( <name>, <list> ). . . . info about library nearrings
##
LibraryNearringInfo := function( arg )
  local N, elms, n, symbols, help, i, k, name, list, string, letters;

  if not ( Length( arg ) in [ 2, 3 ] and IsString( arg[1] ) and 
           IsList( arg[2] ) and ForAll( arg[2], l -> IsInt( l ) ) ) then
    Error( "Usage: LibraryNearringInfo( <name>, <list> ) where <name>",
           " must be a\ngroup name and <list> must be a list of numbers",
           "of classes" );
  fi;

  name := arg[1]; list := arg[2];
  if IsBound( arg[3] ) then letters := arg[3]; else letters := ""; fi;
  if 'C' in letters or 'c' in letters then
    Print( "A ... abstract affine\n" );
    Print( "B ... boolean\n" );
    Print( "C ... commutative\n" );
    Print( "D ... distributive\n" );
    Print( "F ... nearfield\n" );
    Print( "G ... distributively generated\n" );
    Print( "I ... integral\n" );
    Print( "N ... nilpotent\n" );
    Print( "O ... planar\n" );
    Print( "P ... prime\n" );
    Print( "Q ... quasiregular\n" );
    Print( "R ... regular\n" );
    Print( "W ... without non-zero nilpotent elements\n" );
    Print( "---------------------------------------",
           "---------------------------------------\n" );
  fi;
  
  N := LibraryNearring( name, list[ 1 ] );
  elms    := Elements( N );
  n       := Size( elms );
  symbols := [ AbstractGenerator( "n0" ) ];
  if n > 1 then
    symbols := Concatenation( [ AbstractGenerator( "n0" ) ], 
                              AbstractGenerators( "n" , n-1 ) );
  fi;
  Print( "---------------------------------------",
         "---------------------------------------" );
  Print( "\n>>> GROUP: ", N.group.name, "\nelements: " );
  Print( symbols, "\n" );
  
  Print( "addition table:\n" );
  NearringOps.DisplayTable( N, "a" );
  
  Print( "group endomorphisms:\n" );
  for i in [1..Length(Endomorphisms( N.group ))] do
    if i < 10 then 
      Print( i, ":   " ); 
    else
      Print( i, ":  " ); 
    fi;
    Print( List( Endomorphisms(N.group)[i].tfl, e -> symbols[ e ] ), "\n" );
  od;

  Print( "\nNEARRINGS:\n" );
  Print( "---------------------------------------",
         "---------------------------------------" );
  
  for k in list do  
    N := LibraryNearring( name, k );
    Print( "\n", k, ")  phi: ", N.group.phi, ";  " );
    for i in N.group.a_y_i_nrs do Print( i, ";" ); od;
    string := [];
    if IsAbstractAffineNearring( N ) then Add( string, 'A' ); 
    else Add( string, '-' ); fi;
    if IsBooleanNearring( N ) then Add( string, 'B' ); 
    else Add( string, '-' ); fi;
    if IsCommutative( N ) then Add( string, 'C' ); 
    else Add( string, '-' ); fi;
    if IsDistributiveNearring( N ) then Add( string, 'D' ); 
    else Add( string, '-' ); fi;
    if IsNearfield( N ) then Add( string, 'F' ); 
    else Add( string, '-' ); fi;
    if IsDgNearring( N ) then Add( string, 'G' ); 
    else Add( string, '-' ); fi;
    if IsIntegralNearring( N ) then Add( string, 'I' ); 
    else Add( string, '-' ); fi;
    if IsNilpotentNearring( N ) then Add( string, 'N' ); 
    else Add( string, '-' ); fi;
    if IsPlanarNearring( N ) then Add( string, 'O' ); 
    else Add( string, '-' ); fi;
    if IsPrimeNearring( N ) then Add( string, 'P' ); 
    else Add( string, '-' ); fi;
    if IsQuasiregularNearring( N ) then Add( string, 'Q' ); 
    else Add( string, '-' ); fi;
    if IsRegularNearring( N ) then Add( string, 'R' ); 
    else Add( string, '-' ); fi;
    if IsNilpotentFreeNearring( N ) then Add( string, 'W' ); 
    else Add( string, '-' ); fi;
    Print( "  ", string );
    if Identity( N ) <> [] then Print( ";  I = ", 
      symbols[ Position( Elements( N.group ), Identity(N)[1] ) ], "\n" );  
    else
      Print( "\n" );
    fi;
    
    if 'M' in letters or 'm' in letters then
      Print( "multiplication table:" );
      NearringOps.DisplayTable( N, "m" );
    fi;
    if 'I' in letters or 'i' in letters then
      Print( "ideals:\n" );
      help := NearringIdeals( N );
      for i in [1..Length(help)] do
        Print( i, ".  ", List( Elements(help[i]), 
               elm -> symbols[Position(elms,elm)] ), "\n" );
      od;
    fi;
    if 'L' in letters or 'l' in letters then
      Print( "left ideals:\n" );
      help := NearringIdeals( N, "l" );
      for i in [1..Length(help)] do
        Print( i, ".  ", List( Elements(help[i]), 
               elm -> symbols[Position(elms,elm)] ), "\n" );
      od;
    fi;
    if 'R' in letters or 'r' in letters then
      Print( "right ideals:\n" );
      help := NearringIdeals( N, "r" );
      for i in [1..Length(help)] do
        Print( i, ".  ", List( Elements(help[i]), 
               elm -> symbols[Position(elms,elm)] ), "\n" );
      od;
    fi;
    if 'V' in letters or 'v' in letters then
      Print( "invariant subnearrings:\n" );
      help := InvariantSubnearrings( N );
      for i in [1..Length(help)] do
        Print( i, ".  ", List( Elements(help[i]), 
               elm -> symbols[Position(elms,elm)] ), "\n" );
      od;
    fi;
    if 'S' in letters or 's' in letters then
      Print( "subnearrings:\n" );
      help := Subnearrings( N );
      for i in [1..Length(help)] do
        Print( i, ".  ", List( Elements(help[i]), 
               elm -> symbols[Position(elms,elm)] ), "\n" );
      od;
    fi;
    if 'E' in letters or 'e' in letters then
      Print( "nearring endomorphisms: " );
      help := Endomorphisms( N );
      for i in List( help, e -> Position( Endomorphisms(N.group), e ) ) do 
        Print( i, "; " );
      od; Print( "\n" );
    fi;
    if 'A' in letters or 'a' in letters then
      Print( "nearring automorphisms: " );
      help := Automorphisms( N );
      for i in List( help, e -> Position( Endomorphisms(N.group), e ) ) do 
        Print( i, "; " );
      od; Print( "\n" );
    fi;
    Print( "---------------------------------------",
           "---------------------------------------" );
  
  od;
  Print( "\n" );
  return;
end;

#############################################################################
##
#F  NearringInfo( <N>, <string> ). . . . . . . . . . . . info about nearrings
##
NearringInfo := function( arg )
  local N, elms, elms_g, n, symbols, help, i, k, name, list, string, letters;

  if not Length( arg ) in [ 1, 2 ] and IsNearring( arg[1] ) and 
         ( not IsBound( arg[2] ) or IsString( arg[2] )  ) then
    Error( "Usage: NearringInfo( <N>, {string} ) where <N>",
           " must be a\nNearring and {string} must be a string" );
  fi;

  N := arg[1]; 
  if IsBound( arg[2] ) then letters := arg[2]; else letters := ""; fi;
  if 'C' in letters or 'c' in letters then
    Print( "A ... abstract affine\n" );
    Print( "B ... boolean\n" );
    Print( "C ... commutative\n" );
    Print( "D ... distributive\n" );
    Print( "F ... nearfield\n" );
    Print( "G ... distributively generated\n" );
    Print( "I ... integral\n" );
    Print( "N ... nilpotent\n" );
    Print( "O ... planar\n" );
    Print( "P ... prime\n" );
    Print( "Q ... quasiregular\n" );
    Print( "R ... regular\n" );
    Print( "W ... without non-zero nilpotent elements\n" );
    Print( "---------------------------------------",
           "---------------------------------------\n" );
  fi;
  
  elms := Elements( N ); elms_g := Elements( N );
  if IsTransformationNearring( N ) then 
    FindGroup( N ); elms_g := N.PHI;
  fi;
  n       := Size( elms );
  symbols := [ AbstractGenerator( "n0" ) ];
  if n > 1 then
    symbols := Concatenation( [ AbstractGenerator( "n0" ) ], 
                              AbstractGenerators( "n" , n-1 ) );
  fi;
  Print( "---------------------------------------",
         "---------------------------------------\n" );
  DisplayCayleyTable( N );
 
  Print( "PROPERTIES:  " );
 
  string := [];
  if IsAbstractAffineNearring( N ) then Add( string, 'A' ); 
  else Add( string, '-' ); fi;
  if IsBooleanNearring( N ) then Add( string, 'B' ); 
  else Add( string, '-' ); fi;
  if IsCommutative( N ) then Add( string, 'C' ); 
  else Add( string, '-' ); fi;
  if IsDistributiveNearring( N ) then Add( string, 'D' ); 
  else Add( string, '-' ); fi;
  if IsNearfield( N ) then Add( string, 'F' ); 
  else Add( string, '-' ); fi;
  if IsDgNearring( N ) then Add( string, 'G' ); 
  else Add( string, '-' ); fi;
  if IsIntegralNearring( N ) then Add( string, 'I' ); 
  else Add( string, '-' ); fi;
  if IsNilpotentNearring( N ) then Add( string, 'N' ); 
  else Add( string, '-' ); fi;
  if IsPlanarNearring( N ) then Add( string, 'O' ); 
  else Add( string, '-' ); fi;
  if IsPrimeNearring( N ) then Add( string, 'P' ); 
  else Add( string, '-' ); fi;
  if IsQuasiregularNearring( N ) then Add( string, 'Q' ); 
  else Add( string, '-' ); fi;
  if IsRegularNearring( N ) then Add( string, 'R' ); 
  else Add( string, '-' ); fi;
  if IsNilpotentFreeNearring( N ) then Add( string, 'W' ); 
  else Add( string, '-' ); fi;
  Print( string );
  if Identity( N ) <> [] then Print( ";  I = ", 
    symbols[ Position( Elements( N ), Identity(N)[1] ) ], "\n" );  
  else
    Print( "\n" );
  fi;
  
  Print( "\nidempotent elements:\n" );
  Print( List( IdempotentElements( N ), elm ->
    symbols[ Position( elms, elm ) ] ), "\n" );
  Print( "quasiregular elements:\n" );
  Print( List( QuasiregularElements( N ), elm ->
    symbols[ Position( elms, elm ) ] ), "\n" );
  Print( "nilpotent elements:\n" );
  Print( List( NilpotentElements( N ), elm ->
    symbols[ Position( elms, elm[1] ) ] ), "\n" );
  
  if 'I' in letters or 'i' in letters then
    Print( "\nideals:\n" );
    help := NearringIdeals( N );
    for i in [1..Length(help)] do
      Print( i, ".  ", List( Elements(help[i]), 
             elm -> symbols[Position(elms_g,elm)] ), "\n" );
    od;
  fi;
  if 'L' in letters or 'l' in letters then
    Print( "left ideals:\n" );
    help := NearringIdeals( N, "l" );
    for i in [1..Length(help)] do
      Print( i, ".  ", List( Elements(help[i]), 
             elm -> symbols[Position(elms_g,elm)] ), "\n" );
    od;
  fi;
  if 'R' in letters or 'r' in letters then
    Print( "right ideals:\n" );
    help := NearringIdeals( N, "r" );
    for i in [1..Length(help)] do
      Print( i, ".  ", List( Elements(help[i]), 
             elm -> symbols[Position(elms_g,elm)] ), "\n" );
    od;
  fi;
  if 'V' in letters or 'v' in letters then
    Print( "invariant subnearrings:\n" );
    help := InvariantSubnearrings( N );
    for i in [1..Length(help)] do
      Print( i, ".  ", List( Elements(help[i]), 
             elm -> symbols[Position(elms,elm)] ), "\n" );
    od;
  fi;
  if 'S' in letters or 's' in letters then
    Print( "subnearrings:\n" );
    help := Subnearrings( N );
    for i in [1..Length(help)] do
      Print( i, ".  ", List( Elements(help[i]), 
             elm -> symbols[Position(elms,elm)] ), "\n" );
    od;
  fi;
  if 'E' in letters or 'e' in letters then
    Print( "\nnearring endomorphisms:\n" );
    for i in Endomorphisms( N ) do
      Print( List( i.tfl, x -> symbols[x] ), "\n" );
    od; Print( "\n" );
  fi;
  if 'A' in letters or 'a' in letters then
    Print( "nearring automorphisms:\n" );
    for i in Automorphisms( N ) do
      Print( List( i.tfl, x -> symbols[x] ), "\n" );
    od; 
  fi;
  Print( "---------------------------------------",
         "---------------------------------------\n" );
  
  return;
end;
