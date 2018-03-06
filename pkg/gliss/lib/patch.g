#############################################################################
##  Add some support functions for Groups.
#############################################################################
#############################################################################
##
##  GroupOps.Endomorphisms( <G> ). . . . compute the endomorphisms of a group  
##
GroupOps.Endomorphisms := function( G )

  local elms,        # the elements of G
        gens,        # the generators of G
        m,           # the number of generators of G
        n,           # the number of elements of G
        k,           # a loop variable
        TFL,         # a list of transformation lists ( tfl's )
        im,          # the image of an endomorphism candidate
        orders_gens, # a list of the orders of the generators of G
        orders_elms, # a list of the orders of the elements of G
        tfl,         # a transformation list
        E,           # the list of endomorphisms on G
        h,           # an endomorphism candidate
        done;        # a flag variable
        
  elms        := Elements( G );
  gens        := G.generators;
  m           := Length( gens ); 
  n           := Size( elms );
  k           := -1;
  TFL         := [];
  im          := List( [1..m], j -> 1 );
  orders_gens := List( gens, gen -> Order( G, gen ) );
  orders_elms := List( elms, elm -> Order( G, elm ) );
  
  # consider all functions: f: gens -> elms 
  #                         s.t. Order( f(gen) ) divides Order( gen )
  # ( f represented as Length(gens)-tuples "im" of elements in elms )
  while k <= m do
    
    if ForAll( [1..m], j -> RemInt( orders_gens[j], orders_elms[im[j]]) = 0 )
    then
      
      h := GroupHomomorphismByImages( G, G, gens, List( im, j -> elms[j] ) );
      if IsGroupHomomorphism( h ) then
      # NOTE: this additional "if" may indeed seem a little awkward, but 
      # "IsGroupHomomorphism( h )" alone won't work in GAP 3.2. 
        if MappingOps.IsGroupHomomorphism( h ) then
          AddSet( TFL, 
                  List( elms, elm -> Position( elms, Image( h, elm ) ) ) );
        fi;
      fi;
    fi;
    k := 1; done := false;
    while not done and k <= m do
      if im[k] = n then im[k] := 1; k := k+1; 
      else im[k] := im[k]+1; done := true;
      fi;
    od;
  od;
  
  # Put I in the last position  
  TFL := Filtered( TFL, e -> e <> [1..n] );
  Add( TFL, List( [1..n], j -> j ) );
  
  # return the result as group transformations
  E := [];
  for tfl in TFL do
    h := Transformation( G, tfl );
    h.isGroupEndomorphism := true; 
    if Size( Set( tfl ) ) = n then h.isGroupAutomorphism := true; fi;
    Add( E, h );
  od;

  return E;
end;
 
#############################################################################
##
##  GroupOps.Automorphisms( <G> ). . . . compute the automorphisms of a group  
##
GroupOps.Automorphisms := function( G )

  local elms,        # the elements of G
        gens,        # the generators of G
        m,           # the number of generators of G
        n,           # the number of elements of G
        k,           # a loop variable
        TFL,         # a list of transformation lists ( tfl's )
        im,          # the image of an endomorphism candidate
        orders_gens, # a list of the orders of the generators of G
        orders_elms, # a list of the orders of the elements of G
        tfl,         # a transformation list
        E,           # the list of endomorphisms on G
        h,           # an endomorphism candidate
        done;        # a flag variable
        
  elms        := Set(Elements( G ));
  gens        := G.generators;
  m           := Length( gens ); 
  n           := Size( elms );
  k           := -1;
  TFL         := [];
  im          := List( [1..m], j -> 2 );
  orders_gens := List( gens, gen -> Order( G, gen ) );
  orders_elms := List( elms, elm -> Order( G, elm ) );
  # consider all functions: f: gens -> elms 
  #                         s.t. Order( f(gen) ) divides Order( gen )
  # ( f represented as Length(gens)-tuples "im" of elements in elms )
  while k <= m do
    if ( Size( Set( im ) ) = m and
         ForAll( [1..m], j -> orders_gens[j] = orders_elms[im[j]] ) ) then
      h := GroupHomomorphismByImages( G, G, gens, List( im, j -> elms[j] ) );
      if IsBijection( h ) then
      # NOTE: this additional "if" may indeed seem a little awkward, but 
      # "IsGroupHomomorphism( h )" alone won't work in GAP 3.2. 
        if MappingOps.IsGroupHomomorphism( h ) then
          AddSet( TFL, 
                  List( elms, elm -> Position( elms, Image( h, elm ) ) ) );
        fi;
      fi;
    fi;
    k := 1; done := false;
    while not done and k <= m do
      if im[k] = n then im[k] := 2; k := k+1; 
      else im[k] := im[k]+1; done := true;
      fi;
    od;
  od;
  
  # Put I in the last position  
  TFL := Filtered( TFL, e -> e <> [1..n] );
  Add( TFL, List( [1..n], j -> j ) );
  
  # return the result as group transformations
  E := [];
  for tfl in TFL do
    h := Transformation( G, tfl );
    h.isGroupEndomorphism := true; 
    h.isGroupAutomorphism := true; 
    Add( E, h );
  od;

  return E;
end;

#############################################################################
##
##  GroupOps.InnerAutomorphisms( <G> ). . . . compute inner auto's of a group  
##
GroupOps.InnerAutomorphisms := function( G )
  local I,E,g,id,i;

  I := []; E := [];
  for g in Elements( G ) do
    AddSet( I, InnerAutomorphism( G, g ) );
  od;

  for i in I do
    i := AsTransformation( i );
    i.isGroupEndomorphism := true; 
    i.isGroupAutomorphism := true; 
    i.isInnerAutomorphism := true;
    if i <> IdentityTransformation( G ) then
      AddSet( E, i );
    else 
      id := Copy( i );
    fi;
  od;
  Add( E, id );
  return E;
end;

#############################################################################
##
##  GroupOps.DisplayTable( <G> ). . . . . print a Cayley table of a group <G> 
##
GroupOps.DisplayTable := function( G )

  local elms,    # the elements of the group
        n,       # the size of the group
        symbols, # a list of the symbols for the elements of the group
        tw,      # the width of a table
        spc,     # local function which prints the right number of spaces
        bar,     # local function for printing the right length of the bar
        ind,     # help variable, an index
        i,j;     # loop variables
  
  elms    := Elements( G );
  n       := Length( elms );
  symbols := [ AbstractGenerator( "g0" ) ];
  if n > 1 then
    symbols := Concatenation( [ AbstractGenerator( "g0" ) ], 
                              AbstractGenerators( "g" , n-1 ) );
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
#   Print( "The table of a group of order ", n, " will not ",
#          "look\ngood on a screen with ", SizeScreen()[1], " characters per ", 
#          "line.\nHowever, you may want to set your line length to a ",
#          "greater\nvalue by using the GAP function 'SizeScreen'.\n" );
#   return;
# fi;
  Print( "Let:\n" );
  for i in [1..n] do Print( symbols[i], " := ", elms[i], "\n" ); od;
  
  # print the addition table
  Print( "\n  + ", spc( 0, n ), "| " ); 
    for i in [1..n] do Print( symbols[i], spc( i, n ) ); od;
  Print( "\n ---", bar( n ) ); for i in [1..n] do Print( bar( n ) ); od;
  for i in [1..n] do
    Print( "\n  ", symbols[i], spc( i, n ), "| " );
    for j in [1..n] do
      ind := Position( elms, elms[i] * elms[j] );
      Print( symbols[ ind ], spc( ind, n ) );  
    od;
  od;

  Print( "\n\n" );
end;
