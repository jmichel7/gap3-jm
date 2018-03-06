###############################################################################
##
##  csg.g  GLISSANDO ver 1.0  Christof Noebauer   1995, 1996
## 
#############################################################################
##  This source file contains the stuff to compute and classify semigroups. 
#############################################################################
##
#############################################################################
##
#F AllFunctions: create all functions from the set [1..n] into the set [1..n]
##               s.t. they are lex. ordered with [1,...,1] first, and with 
##               the constant function [n,n,...,n] last.
##
AllFunctions := function( n )
  local af,i,j,k,success;

  i := []; 
  for j in [1..n] do i[j] := 1; od;
  k := -1;
  af := [];
  
  while k <= n do
    AddSet( af, Reversed( i ) );
    k := 1; success := false;
    while not success and k <= n do
      if i[k] = n then i[k] := 1; k := k+1; 
      else i[k] := i[k]+1; success := true;
      fi;
    od;
  od;

  return af;
end;

#############################################################################
##
#F  ValidFunctionsSg( m ) . . . . . Determine all functions f: S -> T (= S^S) 
##                                  with f(f(a)(b)) = f(a) o f(b).
##  V1.0  15.2.95
##  input parameter: m.........a positive integer defining a size 
##  return value:    valid_f...a list of all functions with the above property
##                  
##  the format of valid_f is a list of lists s.t. each of those lists 
##  represents a function f: S -> T
##
##  This is an implementation of algorithm 1 on page 10 of the GLISSANDO
##  overview paper.
##
ValidFunctionsSg := function( m )

  local n,       # help var: the size of T
        i,       # help var: a list of indices
        tuples,  # help var: all tuples out of [1..m]
        valid_f, # a list of the functions with the above property
        total,   # total number of functions to be considered
        k,       # a loop variable
        count,   # the number of the function currently being considered
        done,    # help var: indicates when to stop the loop
        f;       # help var: a function to be considered

  Print( "computing all functions...\n" );
  T       := AllFunctions( m );
  n       := m^m;
  i       := List( [1..m], k -> 1 );
  tuples  := Tuples( [1..m], 2 );
  valid_f := [];
  # compute the number of loop executions 
  total     := n^m; 
  Print( "computing valid functions...\n" );

  k     := -1;
  count :=  0; 
  while k <= m and i[m] <= n do
    count := count + 1;
    Print( "  considering function ", count, " of ", total, "...\r" );    
    f := Sublist( T, Reversed( i ) );
    if ForAll( tuples, t -> f[ f[t[1]][t[2]] ] = f[t[1]] { f[t[2]] } ) then  
      Add( valid_f, Reversed( i ) );
    fi;
    k := 1; 
    done := false;
    while not done and k <= m do
      if i[k] = n 
        then i[k] := 1; k := k+1; 
      else 
        i[k] := i[k]+1; done := true;
      fi;
    od;
  od;

  Print( "\n" );
  return valid_f;

end;

#############################################################################
##
#F  ClassifySg( valid_f ). . . . . . . . . . . .determine iso classes of sgps
##  V1.0  15.2.95
##  input parameters: valid_f...a list of lists representing 
##                              valid functions G -> E
##  return value:     classes...a record of iso classes
##                  
##
##  This is an implementation of algorithm 2 on page 10 of the GLISSANDO
##  overview paper.
##
ClassifySg := function( valid_f )

  local R,       # a list, the entries stand for the elements of the set
        AF,      # the list of all functions
        A,       # a list of the bijections: R -> R
        vf,      # help var: the set of lists representing valid functions, 
                 # will be dynamically reduced in each loop execution
        noc,     # contains the number of the current class
        classes, # return value: record which contains the computed classes
        vfcount, # help var: counts how many functions remain
        f,g,     # lists, representing functions in vf
        a,       # a list, representing a bijection
        loa;     # a list of automorphisms ( represented as lists )
  
  Print( "computing all functions...\n" );
  AF      := AllFunctions( Length( valid_f[1] ) );
  R       := [1..Length( AF[1] )];
  Print( "computing bijections...\n" );
  A       := Filtered( AF, e -> Size( Set( e ) ) = Size( R ) ); 
             Add( A, [] );
  vf      := Set( valid_f ); 
  noc     := 0;
  classes := rec();
  vfcount := Length( vf );

  Print( "classifying...\n" );
  for f in vf do
    noc := noc + 1;
    classes.(noc) := rec( phi := f );
    # initialize loa with the identity 
    loa := [ Position( AF, R ) ]; # this adds the identity!
    
    Print( "functions to go: ", vfcount, "    \r" );
    
    Unbind( vf[ Position( vf, f ) ] ); vfcount := vfcount - 1;
    for g in vf do
      a := First( A, a -> a = [] 
               # this is isomorphisms the right way
               or ForAll( R, x -> a{ AF[ f[x] ] } = AF[ g[a[x]] ]{ a } ) 
#               # this is anti-isomorphisms
#               or ForAll( R, x -> ForAll( R, y -> 
#                                    a[AF[f[y]][x]] = AF[g[a[x]]][a[y]] )) 
                );
      if a <> [] then
        AddSet( loa, Position( AF, a ) );
        Unbind( vf[Position( vf, g )] ); vfcount := vfcount - 1;
      fi;
    od;
    classes.(noc).bijs_yielding_iso_sgps := loa;
  od;

  Print( "\n" );
  return classes;
end;

