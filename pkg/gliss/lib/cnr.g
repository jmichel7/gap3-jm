###############################################################################
##
##  cnr.g  GLISSANDO ver 1.0  Christof Noebauer   1995, 1996
##                                                       
#############################################################################
##  This source file contains the stuff to compute and classify nearrings. 
#############################################################################
##
#############################################################################
##
#F  ValidFunctionsNr( <G> ). . Determine all functions f: G -> E ( = End(G) ) 
##                             with f(f(a)(b)) = f(a) o f(b)
##  V1.0 23.2.95  V0.9 6.8.94
##  input parameter: G........a group 
##  return value:    valid_f..a list of all functions with the above property
##                  
##  the format of valid_f is a list of lists s.t. each of those lists 
##  represents a function f: G -> E
##
##  This is an implementation of algorithm 3 on page 16 of the GLISSANDO
##  overview paper. It follows Clay's method to generate 
##  nearrings as Yearby describes it on pp. 7 - 10 of his dissertation.
##
ValidFunctionsNr := function( G )

  local TFL,     # the list of the tf lists of the tf's of the group
        m,       # help var: the number of transformed elements
        n,       # help var: total number of endomorphisms
        i,       # help var: a list of indices
        tuples,  # help var: all tuples out of [1..m]
        valid_f, # a list of the functions with the above property
        idpot,   # number of idempotent endomorphisms
        total,   # total number of functions to be considered
        k,       # a loop variable
        count,   # the number of the function currently being considered
        done,    # help var: indicates when to stop the loop
        f;       # help var: a function to be considered

  Print( "computing endomorphisms...\n" );
  TFL     := List( Endomorphisms( G ), e -> e.tfl );
  m       := Length( TFL[1] );
  n       := Length( TFL );
  i       := List( [1..m], k -> 1 );
  tuples  := Tuples( [1..m], 2 );
  valid_f := [];
  idpot   := 0; 
  # determine the number of idempotent endomorphisms in TFL.
  for k in TFL do if k{ k } = k then idpot := idpot + 1; fi; od;
  # compute the number of loop executions 
  total     := (idpot-1) * n^(m-1); 
  Print( "total number of endomorphisms on ", G, ": ", n, "\n" );
  Print( "total number of idempotent endomorphisms on ", G, ": ", 
         idpot, "\n" );
  Print( "computing valid functions...\n" );

  k     := -1;
  count :=  0; 
  while k <= m and i[m] < n do
    count := count + 1;
    Print( "  considering function ", count, " of ", total, "...\r" );    
    f := Sublist( TFL, Reversed( i ) );
    if ForAll( tuples, t -> f[ f[t[1]][t[2]] ] = f[t[1]] { f[t[2]] } ) then  
      Add( valid_f, Reversed( i ) );
    fi;
    k := 1; 
    done := false;
    while not done and k <= m do
      if i[k] = n then i[k] := 1; k := k+1; 
      else 
        if k = m then 
          repeat
            i[k] := i[k]+1;
          # until the next idempotent or I is found  
          until TFL[i[k]]{ TFL[i[k]] } = TFL[i[k]] or TFL[i[k]] = [1..m]; 
          done := true;
        else
          i[k] := i[k]+1; 
          done := true;
        fi;
      fi;
    od;
  od;

  Print( "\n" );
  # Add [I,...,I] to the list of valid functions
  Add( valid_f, List( [1..m], k -> Position( TFL, [1..m] ) ) );
  return valid_f;

end;

#############################################################################
##
#F  ClassifyNr( <valid_f>, <E> ) . . . . . . . . determine iso classes of nrs
##  V1.0 23.2.95  V0.15 7.10.94
##  input parameters: valid_f...a list of lists representing 
##                              valid functions G -> E
##                    E.........a list of endomorphisms on G
##  return value:     classes...a record of iso classes
##                  
##  This is an implementation of algorithm 4 on page 16 of the GLISSANDO
##  overview paper. 
##
ClassifyNr := function( valid_f, E )

  local TFL,     # the list of the tf lists of the tf's of the group
        R,       # a list, the entries stand for the elements of the group
        A,       # a list of the tfl's that represent automorphisms
        vf,      # help var: the set of lists representing valid functions, 
                 # will be dynamically reduced in each loop execution
        noc,     # contains the number of the current class
        classes, # return value: record which contains the computed classes
        vfcount, # help var: counts how many functions remain
        f,g,     # lists, representing functions in vf
        a,       # a tfl, representing an automorphism in A
        loa;     # a list of automorphisms ( represented as lists )
  
  TFL     := List( E, e -> e.tfl );
  R       := [1..Length( TFL[1] )];
  Print( "computing automorphisms...\n" );
  A       := Filtered( TFL, e -> Size( Set( e ) ) = Size( R ) ); 
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
    loa := [ Length( E ) ]; # ( note that I must be the last entry in E )
    Print( "functions to go: ", vfcount, "   \r" );
    Unbind( vf[ Position( vf, f ) ] ); vfcount := vfcount - 1;
    for g in vf do
      a := First( A, a -> a = [] 
               # this is isomorphisms the right way
               or ForAll( R, x -> a{ TFL[ f[x] ] } = TFL[ g[a[x]] ]{ a } ) 
#               # this is anti-isomorphisms
#               or ForAll( R, x -> ForAll( R, y -> 
#                                    a[TFL[f[y]][x]] = TFL[g[a[x]]][a[y]] )) 
                );
      if a <> [] then
        AddSet( loa, Position( TFL, a ) );
        Unbind( vf[Position( vf, g )] ); vfcount := vfcount - 1;
      fi;
    od;
    classes.(noc).autos_yielding_iso_nrs := loa;
  od;

  Print( "\n" );
  return classes;
end;

#############################################################################
##
#F  ConstructNr( <G>, <classes> ) . . . . . . . . . construct nearring record
##  V1.0 23.2.95  
##  input parameters: G.........a group 
##                    classes...a record of nr classes on G as 
##                              constructed by ClassifyNr
##  return value:     rec.......a nearring record with all neccessary info
##
##  This function just constructs a nearring record ready for storage as
##  a nearring library file.
##
ConstructNr := function( G, classes )
  
  local E, elms, endos, i, NR;

  E := Endomorphisms( G );
  elms  := rec();
  endos := rec();
  for i in [1..Size( G )] do elms.(i) := Elements( G )[i]; od;
  for i in [1..Length( E )] do endos.(i) := E[i].tfl; od;

  for i in [1..Length( RecFields( classes ) )] do
    Sort( classes.(i).autos_yielding_iso_nrs );
  od;

  NR := rec( 
    group_name          := G.name,
    group_generators    := G.generators,
    elements            := elms,
    group_endomorphisms := endos,
    classes             := classes 
  );
  
  return NR;
end;

#############################################################################
##  The two following functions are an implementation of Yearby's method
##  to construct nearrings as he describes it on pp 74 - 91 of his 
##  dissertation.
#############################################################################
#############################################################################
##
##  Extend( <S>, <f> ) . . . . . . . . . . . . extend partial multiplications
##
##  input parameters: S: a list of the numbers of elements that are mapped
##                    f: a list of tfl's ( maybe with holes )
##  Example: S = [ 3, 5 ], f = [ ,,[...],,[...] ]
##
##                  
##  This is an implementation of algorithm 5 on page 24 of the GLISSANDO
##  overview paper. 
##
Extend := function( S, f )
  
  local f1, L1, L2, S2, ST, s, t, c, f2, C, pairs;
  
  f1 := Copy( f ); L1 := Copy( S ); C := [];
  while true do
    # generate L2 and S2
    L2 := Copy( L1 ); S2 := []; ST := []; 
    for s in L1 do
      for t in L1 do
        c := f1[s][t];
        AddSet( L2, c );
        # construct the set C for the compatibility check
        if c in L1 then
          if [ s, t ] in C then 
            RemoveSet( C, [s,t] ); 
          else 
            AddSet( C, [s,t] );
          fi;
        fi;
        # store the pair [ s, t ] for later use
        if not ( c in L1 ) and not ( c in S2 ) then
          AddSet( S2, c );
          ST[c] := [ s, t ];
        fi;
      od;
    od;
    # perform compatibility check
    if not ForAll( C, pairs -> f1[ f1[pairs[1]][pairs[2]] ] = 
                               f1[pairs[1]]{ f1[pairs[2]] } ) then 
      return [ false, [] ];
    fi;
    if S2 = [] then return [ f1, L1 ]; fi;
    # compute f2
    f2 := [];
    for c in L2 do
      if c in L1 then
        f2[c] := f1[c];
      else
        s := ST[c][1]; t := ST[c][2];
        f2[c] := f1[s]{ f1[t] };
      fi;
    od;
    # reassign f1 and L1 and restart the loop
    f1 := f2; L1 := L2; 
  od;
  return;
end;

#############################################################################
##
#F  ValidFunctionsYearby( <E> ). . . . . . Determine all functions f: S -> E 
##                                         with f(f(a)(b)) = f(a) o f(b)
##  V1.0 28.2.95  V0.9 9.11.94
##  input parameter: E....a set s.t. ( E, gamma ) is a Yearby pair w.r.t. S
##  return value:    valid_f...a list of all valid functions 
##                  
##  Examples for Yearby pairs: G a group: ( End(G), Aut(G) ),
##                             S a set :  ( T(G), Sym(G) ).
##
##  the format of valid_f is a list of lists s.t. each of those lists 
##  represents a function f:S -> E
##
##  This is an implementation of the method Yearby introduced in chapter IV
##  of his dissertation.
##  It is based on extending and backtracking partial nearring
##  multiplications.
##
ValidFunctionsYearby := function( E )
  local TFL, m, n, G, rec_set, valid_f, done, l, f, S, ext, elm, i,
        f1, L, diff, count, f2;

  TFL := List( E, e -> e.tfl );
  m   := Length( TFL[1] );
  n   := Length( TFL );
  G   := [1..m];
  rec_set := [ [ [ TFL[1] ], 1 ] ];
  valid_f := [];
  done    := false;
  count   := 0;
  
  repeat
    
    l    := Length( rec_set );
    f    := rec_set[l][1];
    S    := Filtered( [1..Length(f)], i -> IsBound(f[i]) );
    ext  := Extend( S, f );
    elm  := rec_set[l][2];
    i    := Position( TFL, f[elm] );
    f1   := ext[1];
    L    := ext[2];
    diff := Difference( G, L );

    if f1 = false or diff = [] then
      if diff = [] then
        f2 := List( f1, k -> Position( TFL, k ) );
        Add( valid_f, f2 );
      fi;
      if l = 1 then
        if i < n then
          rec_set := [ [ [ TFL[i+1] ], 1 ] ];
        else
          done := true;
        fi;
      else  # if l > 1
        if i < n then
          rec_set[l][1][elm] := TFL[i+1]; # assign a new image
        else
          repeat
            Unbind( rec_set[l] );
            l := l-1;
            elm := rec_set[l][2];
            f := rec_set[l][1];
            i := Position( TFL, f[elm] );
          until i < n or Length( rec_set ) = 1;
          if i < n then
            rec_set[l][1][elm] := TFL[i+1]; # assign a new image
          else
            done := true;
          fi;
        fi;
      fi;
   else  # extension succeeded only partially
      elm := diff[1];
      f1[elm] := TFL[1];
      Add( rec_set, [ f1, elm ] );
    fi;

  count := count + 1;
  
  until done;

  Print( "\nloop executions :", count, "\n\n" );
  return valid_f;
end;
