RequirePackage("chevie");
CoxeterDirectory := "coxeter_1.01";
KazhdanLusztigPolsFromCoxeter := false;
InfoChevie := Print;
firstavalues := function(w)
  # Let w be a Coxeter group w
  # calculates all the a-values of all elements of w
  # The result is a vector with one integer for each element of w
  # the order is that of `CoxeterWords(w)' (not that of `Elements(w)'!)
  local a,c,elnr,h,l,lc,lp,ls,name,s,tabl,u,x,y;
  u := X(Cyclotomics);
  u.name := "u";
  h := Hecke(w,u^2,u);
  c := Basis(h,"C");
  l := CoxeterWords(w);
  lp := List(l,x->PermCoxeterWord(w,x));
  ls := ShallowCopy(lp);
  tabl := [1..Length(l)];
  SortParallel(ls,tabl);
  lc := List(l,c);
  s := List(lc,x->List(lc,y->x*y));
  Print("Have structure constants...\n");
  a := 0*[1..Length(l)];
  for x in [1..Length(l)] do
    for y in [1..Length(l)] do
      for c in [1..Length(s[x][y].coeff)] do
        elnr := tabl[Position(ls,s[x][y].elm[c])];
        if s[x][y].coeff[c].valuation < a[elnr] then
          a[elnr] := s[x][y].coeff[c].valuation;
        fi;
      od;
    od;
    Print(".\c");
  od;
  Print("\n");
  return a;
end;
  
SortLengthLex := function(a,b)
  if Length(a) < Length(b) then return true;
  elif Length(a) > Length(b) then return false;
  else return a < b; fi;
end;

CoxeterElementsLengthLex := function(w)
  local i,l,ll;
  l := [];
  for i in [0..w.N] do
    Append(l,CoxeterElementsLength(w,i));
  od;
  ll := List(l,x->CoxeterWord(w,x));
  SortParallel(ll,l,SortLengthLex);
  return l;
end;

TuneCoxeterGroup := function(w)
  # w must be a CoxeterGroup record. Some additional information is calculated
  # and stored in the record to improve performance. 
  local i,j,k,l,len,ll,lll,r,s;

  # enumerate all elements, then they are in w.elements:
  Elements(w);

  # first enumerate all elements by length, sort words lexicographically
  # within one length
  w.elslenlex := [];
  w.lennr := [];
  w.lenfirst := [];
  w.lens := [];
  w.words := [];
  for i in [0..w.N] do   # all possible lengths
    ll := CoxeterElementsLength(w,i);
    lll := List(ll,x->CoxeterWord(w,x));
    SortParallel(lll,ll);   # sort lexicographically
    Add(w.lennr,Length(ll));
    Add(w.lenfirst,Length(w.elslenlex)+1);
    Append(w.elslenlex,ll);
    Append(w.lens,0*[1..Length(ll)]+i);
    Append(w.words,lll);
  od;
  w.elslenlex := CoxeterElementsLengthLex(w);

  # Store descent sets as bit lists:
  r := [1..w.semisimpleRank];
  w.ldescent := List(w.elslenlex,x->BlistList(r,LeftDescentSet(w,x)));
  w.rdescent := List(w.elslenlex,x->BlistList(r,RightDescentSet(w,x)));

  # create a lookup from sorted elements into elslenlex:
  l := ShallowCopy(w.elslenlex);
  w.elslookup := [1..Size(w)];
  SortParallel(l,w.elslookup);

  # calculate mult table from left by left-descent-generators:
  w.lmulttab := [];
  for i in [1..Size(w)] do
    l := [];
    for j in [1..w.semisimpleRank] do
      Add(l,w.elslookup[Position(w.elements,w.generators[j]*w.elslenlex[i])]);
    od;
    Add(w.lmulttab,l);
  od;

  # calculate inv table:
  w.invtab := [];
  for i in [1..Size(w)] do
    Add(w.invtab,w.elslookup[Position(w.elements,w.elslenlex[i]^-1)]);
  od;

  # calculate a bruhat table (triangular bitlist):
  # w.bruhat[i][j] iff Bruhat(w,w.elslenlex[j],w.elslenlex[i])
  w.bruhat := List([1..Size(w)],i->BlistList([1..Size(w)],[]));
  w.bruhat[1][1] := true;
  for len in [1..w.N] do
    InfoChevie("Bruhat order for elements of length ",len,"...\n");
    for i in [w.lenfirst[len+1]..w.lenfirst[len+1]+w.lennr[len+1]-1] do
      for j in [w.lenfirst[len]..w.lenfirst[len]+w.lennr[len]-1] do
        # Slowly:
        #if Bruhat(w,w.elslenlex[j],w.elslenlex[i],len-1,len) then
        #  UniteBlist(w.bruhat[i],w.bruhat[j]);
        #fi;
        # We go faster and use:
        # let sy < y and (z := sx if sx<x and z:=x otherwise)
        # then x <= y iff z <= sy
        # Proof: Think about subwords!
        s := First([1..w.semisimpleRank],k->w.ldescent[i][k]);
        if w.ldescent[j][s] then
          k := w.lmulttab[j][s];
        else
          k := j;
        fi;
        if w.bruhat[w.lmulttab[i][s]][k] then
          UniteBlist(w.bruhat[i],w.bruhat[j]);
        fi;
      od;
      w.bruhat[i][i] := true;
    od;
  od;

  # prepare for Kazhdan-Lusztig-Polynomials:
  w.fastklpols := [];
  w.extrpairs := List([1..Size(w)],x->[]);
  w.extrpols := List([1..Size(w)],x->[]);

  # mark group as tuned:
  w.isTuned := true;
end;

ExtremalPair := function(W,y,w)
  # Calculates an extremal pair for y,w, given as numbers in W.elslenlex.
  # `TuneCoxeterGroup' must have been called. Only feasible for y <= w.
  local s;
  if not(W.bruhat[w][y]) then
    return false;
  fi;
  if W.invtab[w] < w then
    w := W.invtab[w];
    y := W.invtab[y];
  fi;
  if y = w then return [w,w]; fi;
  s := 1;
  while s <= W.semisimpleRank do
    if W.ldescent[w][s] and not(W.ldescent[y][s]) then
      y := W.lmulttab[y][s];
      if y = w then return [w,w]; fi;
      s := 0;   # start again with first generator
    elif W.rdescent[w][s] and not(W.rdescent[y][s]) then
      # the following would be   W.rmulttab[y][s]
      # if we had an rmulttab:
      y := W.invtab[W.lmulttab[W.invtab[y]][s]];
      if y = w then return [w,w]; fi;
      s := 0;   # start again with first generator
    fi;
    s := s + 1;
  od;
  return [y,w];
end;

StoreKazhdanLusztigPolsFromCoxeter := function(W)
  # The output of `coxtogap.py' for the correct type must have been read.
  # Uses the globale variable `KazhdanLusztigPolsFromCoxeter'.
  # Restores it to its original value.
  local d2,i,j,l,len,p,pos,r;
  r := KazhdanLusztigPolsFromCoxeter;
  if not(IsBound(W.isTuned)) then
    TuneCoxeterGroup(W);
  fi;
  InfoChevie("Calculating mue-not-zero table...\n");
  W.fastklpols := r.pols;
  W.mueneq := List([1..Size(W)],x->BlistList([1..Size(W)],[]));
  for p in r.pairs do
    if p[1] <> p[2] then
      pos := PositionSorted(W.extrpairs[p[1]],p[2]);
      l := Length(W.extrpairs[p[1]]);
      W.extrpairs[p[1]]{[pos+1..l+1]} := W.extrpairs[p[1]]{[pos..l]};
      W.extrpairs[p[1]][pos] := p[2];
      W.extrpols[p[1]]{[pos+1..l+1]} := W.extrpols[p[1]]{[pos..l]};
      W.extrpols[p[1]][pos] := p[3];
      d2 := W.lens[p[1]] - W.lens[p[2]] - 1;
      if IsEvenInt(d2) and d2/2 = Length(W.fastklpols[p[3]])-1 then
        W.mueneq[p[1]][p[2]] := true;
        W.mueneq[p[2]][p[1]] := true;
        W.mueneq[W.invtab[p[1]]][W.invtab[p[2]]] := true;
        W.mueneq[W.invtab[p[2]]][W.invtab[p[1]]] := true;
        if W.fastklpols[p[3]][d2/2+1] <> 1 then
          InfoChevie("Mue not equal to 1 for pair: ",p,"\n");
        fi;
      fi;
    fi;
  od;
  # Now mueneq is not correct!
  # We have to add all pairs (y,w) with l(y) = l(w)-1 and y <= w
  for len in [1..W.N] do
    for i in [W.lenfirst[len]..W.lenfirst[len]+W.lennr[len]-1] do
      for j in [W.lenfirst[len+1]..W.lenfirst[len+1]+W.lennr[len+1]-1] do
        if W.bruhat[j][i] then
          W.mueneq[j][i] := true;
          W.mueneq[i][j] := true;
        fi;
      od;
    od;
  od;
  # Just for the sake of completeness: add diagonal
  for i in [1..Size(W)] do
    W.mueneq[i][i] := true;
  od;
  if IsBound(KazhdanLusztigPolsFromCoxeter.orig) then
    KazhdanLusztigPolsFromCoxeter := KazhdanLusztigPolsFromCoxeter.orig;
  else
    Unbind(KazhdanLusztigPolsFromCoxeter);
  fi;
end;

#############################################################################
##
#F  FastKazhdanLusztigPolynomial( <W>, <y>, <w> )  . . . . . . . . . . . . .
##  . . . . . . . . . . . . . . . . . . . . . . . Kazhdan-Lusztig polynomial
##
##  'FastKazhdanLusztigPolynomial'  returns  the Kazhdan - Lusztig polynomial 
##  as coefficient list over the rationals corresponding to the elments <y> 
##  and <w>  (given as numbers in elslenlex) in the Coxeter group <W>. 
##  `TuneCoxeterGroup' must have been called beforehand!
##
FastKazhdanLusztigPolynomial:=function(W,y,w)
  local ep,l,p,po,pos,pos2;
  if not W.bruhat[w][y] then
    return [];
  fi;
  ep:=ExtremalPair(W,y,w);
  if ep[1] = ep[2] then
    return [1];
  fi;
  po := Position(W.extrpairs[ep[2]],ep[1]);
  if po = false then
    Print("Did not find Kazhdan-Lusztig-Polynomial of extremal pair ",ep,"\n");
    p := KazhdanLusztigPolynomial(W,W.elslenlex[ep[1]],W.elslenlex[ep[2]]);
    p := p.coefficients;
    pos := PositionSorted(W.extrpairs[ep[2]],ep[1]);
    l := Length(W.extrpairs[ep[2]]);
    W.extrpairs[ep[2]]{[pos+1..l+1]} := W.extrpairs[ep[2]]{[pos..l]};
    W.extrpairs[ep[2]][pos] := ep[1];
    W.extrpols[ep[2]]{[pos+1..l+1]} := W.extrpols[ep[2]]{[pos..l]};
    pos2 := Position(W.fastklpols,p);
    if pos2 = false then
      AddSet(W.fastklpols,p);
      pos2 := Length(W.fastklpols);
    fi;
    W.extrpols[ep[2]][pos] := pos2;
    return p;
  fi;
  return W.fastklpols[W.extrpols[ep[2]][po]];
end;

PrintBoolArray := function(M)
  local i,r;
  for r in M do
    for i in r do
      if i = true then Print("+"); else Print("."); fi;
    od;
    Print("\n");
  od;
end;

RightCells := function(w)
  # Calculates the right cells of w. Uses `LeftCells' and inversion.
  # Facts: x <= y  <=>  x^-1 <= y^-1   (Bruhat-Order!)
  #        P_{x,y} = P_{x^-1,y^-1}     --> Brenti
  #        L(x) = R(x^-1)
  # ==> right cells are "inverses" of left cells.
  return List(LeftCells(w),c->[List(c[1],x->ReducedCoxeterWord(w,Reversed(x))),
                               ShallowCopy(c[2])]);
end;

ClosureEquivalenceRelation := function(ll,rr)
  # ll and rr are both partitions of the set [1..n] for some n (lists of sets).
  # Returns the partition of [1..n] corresponding to the closure of both
  # equivalence relations.
  local c,cc,g,g1,g2,grouptab,mergetab,newkey,postab,result;
  postab := [];
  for c in [1..Length(ll)] do
    for cc in ll[c] do
      postab[cc] := c;
    od;
  od;
  mergetab := [1..Length(ll)];
  grouptab := List([1..Length(ll)],i->[i]);
  newkey := Length(ll)+1;
  # the first entry is the key
  for c in rr do
    for cc in [1..Length(c)-1] do
      g1 := mergetab[postab[c[cc]]];   # first group
      g2 := mergetab[postab[c[cc+1]]]; # second group
      if not(IsIdentical(grouptab[g1],grouptab[g2])) then
        # merge groups g1 and g2:
        Append(grouptab[g1],grouptab[g2]);
        grouptab[g2] := grouptab[g1];
      fi;
    od;
  od;
  result := rec(partition := [],mergetab := []);
  for g in [1..Length(grouptab)] do
    if grouptab[g][1] <> false then  # not yet done
      Add(result.partition,Concatenation(ll{grouptab[g]}));
      Add(result.mergetab,ShallowCopy(grouptab[g]));
      grouptab[g][1] := false;  # done
    fi;
  od;
  return result;
end;

DufloElements := function(W)
  local a,a1,c,d,i,inv,lc,pc;
  lc:=W.leftcells;
  W.duflos := [];
  W.avalues := [];
  InfoChevie("#I Searching for Duflo involutions...\n");
  for c in lc do
    InfoChevie(".\c");
    pc:=W.elslenlex{c};
    inv:=Filtered([1..Length(c)],x->pc[x]^2=());
    d:=inv[1];
    a:=W.lens[c[d]]-2*(Length(FastKazhdanLusztigPolynomial(W,1,c[d]))-1);
    for i in [2..Length(inv)] do
      a1:=W.lens[c[inv[i]]]-
           2*(Length(FastKazhdanLusztigPolynomial(W,1,c[inv[i]]))-1);
      if a1<a then
        d:=inv[i];
        a:=a1;
      fi;
    od;
    Add(W.duflos,c[d]);
    Add(W.avalues,a);
  od;
  InfoChevie(" done.\n");
  W.rlcellsavalues := List([1..Length(W.rlcells)],
                           i->W.avalues[W.leftcellfusion[i][1]]);
  return W;
end;

Lless := function( W,x,y )
  return W.mueneq[x][y] and not(IsSubsetBlist(W.ldescent[y],W.ldescent[x]));
end;

Rless := function( W,x,y )
  return W.mueneq[x][y] and not(IsSubsetBlist(W.rdescent[y],W.rdescent[x]));
end;

# From chevie (improved for performance):
MyDecomposedLeftCells := function ( W, c)
  # W is a CoxeterGroup already prepared by `TuneCoxeterGroup' and with
  # loaded data from "coxeter" via `StoreKazhdanLusztigPolsFromCoxeter'.
  # c is a list of numbers of group elements in the ordering given by
  # `W.elslenlex'. They are all group elements having the same right
  # descent set and the lengths of the elements are non-decreasing.
  # Searches the distribution into left cells by looking at mue values
  # and left descent sets.
  #from now on it seems no longer necessary, that the elements
  #are in non-decreasing length.
  local ce,d,i,j,l,l1,l2,rest,s1,s2;
  ce := [  ];
  rest := Set(c);
  while rest <> [  ]  do
      l := Length(rest);
      s1 := [ rest[1] ];
      s2 := ShallowCopy(s1);
      l1 := BlistList(rest,[rest[1]]);
      for j  in s1  do
          for i  in [2..Length(rest)] do
              if not(l1[i]) and Lless( W, j, rest[i] )  then
                  Add( s1, rest[i] );
                  l1[i] := true;
              fi;
          od;
      od;
      l2 := BlistList(rest,[rest[1]]);
      for j  in s2  do
          for i  in [ 2 .. Length( rest ) ]  do
              if not(l2[i]) and Lless( W, rest[i], j )  then
                  Add( s2, rest[i] );
                  l2[i] := true;
              fi;
          od;
      od;
      IntersectBlist( l1, l2 );
      d := ListBlist( rest, l1 );
      Add( ce, d );
      SubtractSet( rest, d );
  od;
  return ce;
end;

#############################################################################
##
#F  LeftCells( <W> )  . . . . . . . . . . . . . . . . . . the left cells of W
##
##  'LeftCells' returns a list of pairs. The first  component  of  each  pair
##  consists of the reduced words in the Coxeter group <W> which lie in 
## one left cell C,  the  second  component  consists  of the corresponding 
##  matrix of  highest coefficients mue(y,w), where y,w are in C.
##  Tuned by Max for efficiency.
##
MyLeftCells := function ( W )
  local i,res,rw,t,x;
  if not(IsBound(W.isTuned)) then
    TuneCoxeterGroup(W);
  fi;
  t := List(Combinations( [ 1 .. W.semisimpleRank ] ),
            x->BlistList([1..W.semisimpleRank],x));
  Sort(t);
  rw := List([ 1 .. Length( t ) ],j->[  ]);
  for i in [1..Size(W)] do
    Add( rw[Position( t, W.rdescent[i] )],i );
  od;
  # Important: In every list in rw the length of the elements is non-decreasing
  res := [  ];
  for i  in [1..Length(rw)]  do
      InfoChevie( "#I  R(w) = ", Filtered([1..W.semisimpleRank],x->t[i][x]),
		   " :  Size = ", Length( rw[i] ) , "\c" );
      x := MyDecomposedLeftCells( W, rw[i]);
      if Length( x ) = 1  then
          InfoChevie( ",  ", Length( x ), " new cell \n" );
      else
          InfoChevie( ",  ", Length( x ), " new cells \n" );
      fi;
      Append( res, x );
  od;
  W.leftcells := res;
  return res;
end;

# again own code:

RLCells := function(w)
  # Calculates RL cells by calculating left and right cells and
  # merging. Changes the group object by adding some fields:
  # leftcells: list of left cells, each cell is a list of numbers of 
  #            elements of w in the order as in elslenlex
  # rightcells: list of right cells, analogously
  # rlcells: list of rl cells, given as list of numbers of
  #          elements of w each
  # leftcellfusion: for each left cell the number of the rl cell
  #                 in which it is contained
  # 
  # `DufloElements' is called afterwards, such that even more fields
  # are filled, see documentation there. Especially all a invariants
  # are calculated.
  local c,cc,i,t,r;
  if not(IsBound(w.isTuned)) then
    TuneCoxeterGroup(w);
  fi;
  InfoChevie("#I Calculating LeftCells...\n");
  MyLeftCells(w);
  InfoChevie("#I ",Length(w.leftcells)," left cells. ",
             "Calculating RightCells from LeftCells...\n");
  w.rightcells := [];
  for c in w.leftcells do
    cc := [];
    for i in c do
      Add(cc,w.invtab[i]);
    od;
    Sort(cc);
    Add(w.rightcells,cc);
  od;
  # Now we merge:
  InfoChevie("#I Calculating RLCells...\n");
  t := ClosureEquivalenceRelation(w.leftcells,w.rightcells);
  w.rlcells := t.partition;
  w.leftcellfusion := t.mergetab;
  InfoChevie("#I ",Length(w.rlcells)," rl cells. Done.\n");
  return DufloElements(w);
end;

Lleq := function ( W, x, y )
  local  i;
  if x = y  then
    return false;
  elif CoxeterLength(W,x) <= CoxeterLength(W,y) and 
       KazhdanLusztigMue(W,x,y) = 0  then
    return false;
  elif CoxeterLength(W,x) >= CoxeterLength(W,y) and 
       KazhdanLusztigMue(W,y,x) = 0  then
    return false;
  else
    for i  in [ 1 .. W.semisimpleRank ]  do
      if W.rootInclusion[i]^x>W.parentN
                        and W.rootInclusion[i]^y<=W.parentN  then
        return true;
      fi;
    od;
    return false;
  fi;
end;

PosetByPairs := function(l,n)
  # l is a list of pairs of integers in [1..n]. A pair [a,b] means "a<=b".
  # This function constructs the poset with elements [1..n] generated by
  # these pairs.
  # Ignores pairs with identical entries.
  local FindDown,FindUp,a,known,li,next,nrcomps,p,todo;
  p := rec(domain := [1..n]);     # the list of points in the poset
  p.n := n;
  p.maxes := List([1..n],x->[]);  # for each point the list of maximals
  p.mins := List([1..n],x->[]);   # for each point list of min. superpoints
  p.gt := [];         # all points greater
  p.leq := [];        # all points less than or equal
  # in the beginning, those two may contain to many points!
  p.doneup := List([1..n],x->false);   # did we already go up?
  p.donedo := List([1..n],x->false);   # did we already go down?
  for a in l do    # first collect the information  
    if a[1] <> a[2] then
      AddSet(p.maxes[a[2]],a[1]);
      AddSet(p.mins[a[1]],a[2]);
    fi;
  od;
  p.comps := [];     # list of connected compontents
  p.compnr := [];    # for each point nr. of component it is in
  known := 0;        # we do not know any point yet
  nrcomps := 0;      # number of connected components known

  FindUp := function(p,a)
    # calculates the set of all points x > a and for all those x the
    # set of all points y > x. Returns list of points.
    # `todo' is changed: all points newly found are added.
    local b,greater;
    greater := [];
    for b in p.mins[a] do   # those are greater
      if not(IsBound(p.compnr[b])) then
        p.compnr[b] := nrcomps;
        AddSet(p.comps[nrcomps],b);
        Add(todo,b);
      fi;
      if not(p.doneup[b]) then   # not yet done!
        UniteSet(greater,FindUp(p,b));
      else
        UniteSet(greater,p.gt[b]);
      fi;
    od;
    # now "greater" is the set of points reachable by 2 or more steps down
    SubtractSet(p.mins[a],greater);
    #Print("SubtractUp: ",a," ",p.mins[a],greater,"\n");
    UniteSet(greater,p.mins[a]);
    p.doneup[a] := true;
    p.gt[a] := greater;
    return greater;
  end;
  
  FindDown := function(p,a)
    # calculates the set of all points x < a and for all those x the
    # set of all points y < x. Returns list of points.
    # `todo' is changed: all points newly found are added.
    local b,lessthan;
    lessthan := [];
    for b in p.maxes[a] do   # those are greater
      if not(IsBound(p.compnr[b])) then
        p.compnr[b] := nrcomps;
        AddSet(p.comps[nrcomps],b);
        Add(todo,b);
      fi;
      if not(p.donedo[b]) then   # not yet done!
        UniteSet(lessthan,FindDown(p,b));
      else
        UniteSet(lessthan,p.leq[b]);
      fi;
    od;
    # now "lessthan" is the set of points reachable by 2 or more steps down
    SubtractSet(p.maxes[a],lessthan);
    #Print("SubtractDown: ",a," ",p.maxes[a],lessthan,"\n");
    UniteSet(lessthan,p.maxes[a]);
    p.donedo[a] := true;
    p.leq[a] := lessthan;
    return lessthan;   # Note that a is included in there!
  end;

  # enumerate connected components by starting at new points
  # and going up and down:
  next := 0;
  while known < n do
    # find next point that we do not know yet:
    repeat
      next := next + 1;
    until not(IsBound(p.compnr[next]));
    # begin a new connected component:
    nrcomps := nrcomps + 1;
    p.compnr[next] := nrcomps;
    p.comps[nrcomps] := [next];   # our current connected component 
    li := [next];     # current new list
    while li <> [] do   # as long as there are new points
      todo := [];
      for a in li do
        if not(p.doneup[a]) then FindUp(p,a); fi;
        if not(p.donedo[a]) then FindDown(p,a); fi;
      od;
      # note that todo is modified in `FindUp' and `FindDown'
      li := todo;
    od;
    known := known + Length(p.comps[nrcomps]);
  od;

  for a in [1..n] do
    AddSet(p.leq[a],a);
  od;
  Unbind(p.doneup);
  Unbind(p.donedo);
  return p;
end;

CalcLevels := function(p)
  # Distributes points into levels. All minimals are in Level 0.
  # In level 1 are all points whose maximals are all in level 0.
  # In level 2 are all points whose maximals are all in levels < 2.
  # ...
  local c,done,i,j,le,news;
  p.level := [];
  for c in [1..Length(p.comps)] do
    # First the minimals:
    news := Filtered(p.comps[c],i->p.maxes[i] = []);
    for i in news do
      p.level[i] := 0;
    od;
    done := Length(news);
    le := 0;
    while done < Length(p.comps[c]) do
      le := le + 1;   # a new level
      news := [];
      for i in p.comps[c] do
        if not(IsBound(p.level[i])) and 
           ForAll(p.maxes[i],j->IsBound(p.level[j])) then
          Add(news,i);
        fi;
      od;
      for j in news do
        p.level[j] := le;
      od;
      done := done + Length(news);
    od; 
  od;
  return p;
end;

PrintToDaVinciPoset := function(f,p,c)
  local i,ii,j,jj;
  PrintTo(f,"[\n");   # Initialize file
  for i in [1..Length(p.comps[c])] do
    ii := p.comps[c][i];
    # First the node with its attribute itself:
    AppendTo(f," l(\"",String(ii),"\",n(\"nodetype\",[a(\"OBJECT\",\"",
               String(ii),"\")],\n  [\n");
    # Now the list of children which is the list of edges down:
    if p.maxes[ii] = [] then
      AppendTo(f,"  ]\n ))");
    else
      for j in [1..Length(p.maxes[ii])] do
        jj := p.maxes[ii][j];
        AppendTo(f,"  e(\"edgetype\",[],r(\"",String(jj),"\"))");
        if j <> Length(p.maxes[ii]) then
          AppendTo(f,",\n");
        else
          AppendTo(f,"\n  ]\n ))");
        fi;
      od;
    fi;

    if i <> Length(p.comps[c]) then
      AppendTo(f,",\n");
    else
      AppendTo(f,"\n]\n");
    fi;
  od;
end;

PrintToStringPoset := function(f,p,c)
  local i,j;
  AppendTo(f,"(\n");
  for i in p.comps[c] do
    AppendTo(f," (");
    AppendTo(f,String(i));
    AppendTo(f," (");
    if p.mins[i] <> [] then
      AppendTo(f,String(p.mins[i][1]));
    fi;
    for j in p.mins[i]{[2..Length(p.mins[i])]} do
      AppendTo(f," ");
      AppendTo(f,String(j));
    od;
    AppendTo(f,"))\n");
  od;
  AppendTo(f,")\n");
end;

PrintToStringPoset2 := function(f,p,c)
  local i,j;
  PrintTo(f,"<title>Display Diamond</title>\n");
  AppendTo(f,"<hr>\n");
  AppendTo(f,"<applet codebase=\"http://www.math.vanderbilt.edu/~pjipsen/gap/\"\n code=\"Poset.class\" width=800 height=600>\n");
  AppendTo(f,"<param name=maxX value=\"25\">\n"); 
  AppendTo(f,"<param name=maxY value=\"");
  AppendTo(f,String(Maximum(p.level)));
  AppendTo(f,"\">\n"); 
  AppendTo(f,"<param name=poset value='"); 
  for i in p.comps[c] do
    AppendTo(f,"\n[");
    AppendTo(f,String(i));
    AppendTo(f,",");
    AppendTo(f,String(Random([0..25])));
    AppendTo(f,",");
    AppendTo(f,String(p.level[i]));
    AppendTo(f,",\"");
    AppendTo(f,String(i));
    AppendTo(f,"\",[");
    if p.mins[i] <> [] then
      AppendTo(f,String(p.mins[i][1]));
    fi;
    for j in p.mins[i]{[2..Length(p.mins[i])]} do
      AppendTo(f,",");
      AppendTo(f,String(j));
    od;
    AppendTo(f,"]],");
  od;
  AppendTo(f,"'>\n</applet>\n"); 
  AppendTo(f,"<hr>\n"); 
  AppendTo(f,"<a href=\"http://www.math.vanderbilt.edu/~pjipsen/gap/",
             "Poset.java\">The source.</a>\n");
end;

PrintToXGAPPoset := function(f,p,c)
  local i,j,l,le,s;
  PrintTo(f,"s := GraphicPoset(\"Poset\",1024,768);\n");
  l := p.comps[c];
  le := p.level{l};
  for i in Set(le) do
    AppendTo(f,"CreateLevel(s,",-i,",\"",i,"\");\n");
  od;
  AppendTo(f,"PosetVertices := [];\n");
  for i in l do
    AppendTo(f,"Add(PosetVertices,Vertex(s,",i,",rec(levelparam := ",
               -p.level[i],", label := \"",i,"\")));\n");
  od;
  for i in [1..Length(l)] do
    for j in p.maxes[l[i]] do
      AppendTo(f,"Edge(s,PosetVertices[",i,"],PosetVertices[",
                 Position(l,j),"]);\n");
    od;
  od;
end;

PosetLeftCells := function(w)
  # Takes a CoxeterGroup after "RLCells" and builds a poset of the 
  # left cells with <=_L as partial order.
  local i,j,k,l,len,li,pairs,rd,yes;
  len := Length(w.leftcells);
  li := w.leftcells;
  rd := List(w.leftcells,x->w.rdescent[x[1]]);
  pairs := [];
  for i in [1..len] do
    for j in [1..len] do
      if i <> j and IsSubsetBlist(rd[i],rd[j]) then   
        # it might be that leftcell[i] <=_L leftcell[j] !
        yes := false;
        k := 1;
        while not(yes) and k <= Length(li[i]) do
          l := 1;
          while not(yes) and l <= Length(li[j]) do
            if Lless(w,li[i][k],li[j][l]) then
              yes := true;
            fi;
            l := l + 1;
          od;
          k := k + 1;
        od;
        if yes then Add(pairs,[i,j]); fi;
      fi;
    od;
  od;
  w.pairsforLlt := pairs;
  return PosetByPairs(pairs,Length(w.leftcells));
end;

PosetRightCells := function(w)
  # Takes a CoxeterGroup after "RLCells" and builds a poset of the 
  # right cells with <=_R as partial order.
  # Note: Gives nothing new because exactly the same as poset of left cells.
  local i,j,k,l,len,re,pairs,ld,yes;
  len := Length(w.rightcells);
  re := w.rightcells;
  ld := List(w.rightcells,x->w.ldescent[x[1]]);
  pairs := [];
  for i in [1..len] do
    for j in [1..len] do
      if i <> j and IsSubsetBlist(ld[i],ld[j]) then   
        # it might be that rightcell[i] <=_R rightcell[j] !
        yes := false;
        k := 1;
        while not(yes) and k <= Length(re[i]) do
          l := 1;
          while not(yes) and l <= Length(re[j]) do
            if Rless(w,re[i][k],re[j][l]) then
              yes := true;
            fi;
            l := l + 1;
          od;
          k := k + 1;
        od;
        if yes then Add(pairs,[i,j]); fi;
      fi;
    od;
  od;
  w.pairsforRlt := pairs;
  return PosetByPairs(pairs,Length(w.rightcells));
end;

IsSuperExtremal := function(W,y,w)
  local s;
  for s in [1..W.semisimpleRank] do
    if W.ldescent[w][s] and not(W.bruhat[W.lmulttab[w][s]][y]) then
      return false;
    fi;
    if W.rdescent[w][s] and 
       not(W.bruhat[W.invtab[W.lmulttab[W.invtab[w]][s]]][y]) then
      return false;
    fi;
  od;
  return true;
end;

#############################################################################
##
#F  MyLeftCellRepresentation( <H> , <cell> ) . . . . . . . . . representation
#F  of the Hecke algebra associated to one or more left cells of a CoxeterGroup
##
##  'MyLeftCellRepresentation' returns a list of matrices giving the left  cell 
##  representation of the  Hecke algebra <H> associated to <cell>, given as
##  a list of numbers of group elements in the CoxeterGroup record
##  underlying H.
##
MyLeftCellRepresentation := function ( H, cell )
    local W,c,f1,f2,i,j,p,rep,s,t,v;
    if ForAny([1..Length(H.parameter)],i->not IsBound(H.sqrtParameter[i]))
    then
      Error("sqrtParameters must be bound to compute cell representations");
    elif Length( Set( H.parameter ) ) > 1 then
      Error("cell representations for unequal parameters not yet implemented");
    else 
      v := H.sqrtParameter[1];
    fi;
    c := Set(cell);
    W:=CoxeterGroup(H);
    ##pc := List( c, function ( i ) return PermCoxeterWord ( W, i ); end );
    rep := [  ];
    for s  in [1..Length(W.generators)]  do
        f1 := [  ];
        for i  in c  do
            if not(W.ldescent[i][s]) then
                Add( f1, true );
            else
                Add( f1, false );
            fi;
        od;
        f2 := List( c, i->W.lmulttab[i][s] );
        t := [  ];
        for i  in [ 1 .. Length( c ) ]  do
            t[i] := [  ];
            for j  in [ 1 .. Length( c ) ]  do
                if i = j  then
                    if f1[i]  then
                        t[i][i] := v ^ 2;
                    else
                        t[i][i] := -v^0;
                    fi;
                elif i > j  then
                    if f2[j] = c[i]  then
                        t[i][j] := v;
                    else
                        t[i][j] := 0*v;
                    fi;
                else
                    if f1[j] and W.mueneq[c[j]][c[i]] and not f1[i]  then
                        p := FastKazhdanLusztigPolynomial(W,c[i],c[j]);
                        t[i][j] := p[Length(p)] * v;
                    else
                        t[i][j] := 0*v;
                    fi;
                fi;
            od;
        od;
        Add( rep, t );
    od;
    return rep;
end;

#############################################################################
##
#F  MyRightCellRepresentation( <H> , <cell> ) . . . . . . . . . representation
#F  of the Hecke algebra associated to one or more right cells of a CoxeterGroup
##
##  'MyRightCellRepresentation' returns a list of matrices giving the rightcell 
##  representation of the  Hecke algebra <H> associated to <cell>, given as
##  a list of numbers of group elements in the CoxeterGroup record
##  underlying H.
##
MyRightCellRepresentation := function ( H, cell )
    local W,c,f1,f2,i,j,p,rep,s,t,v;
    if ForAny([1..Length(H.parameter)],i->not IsBound(H.sqrtParameter[i]))
    then
      Error("sqrtParameters must be bound to compute cell representations");
    elif Length( Set( H.parameter ) ) > 1 then
      Error("cell representations for unequal parameters not yet implemented");
    else 
      v := H.sqrtParameter[1];
    fi;
    c := Set(cell);
    W:=CoxeterGroup(H);
    ##pc := List( c, function ( i ) return PermCoxeterWord ( W, i ); end );
    rep := [  ];
    for s  in [1..Length(W.generators)]  do
        f1 := [  ];
        for i  in c  do
            if not(W.rdescent[i][s]) then
                Add( f1, true );
            else
                Add( f1, false );
            fi;
        od;
        f2 := List( c, i->W.invtab[W.lmulttab[W.invtab[i]][s]] );
        t := [  ];
        for i  in [ 1 .. Length( c ) ]  do
            t[i] := [  ];
            for j  in [ 1 .. Length( c ) ]  do
                if i = j  then
                    if f1[i]  then
                        t[i][i] := v ^ 2;
                    else
                        t[i][i] := -v^0;
                    fi;
                elif i > j  then
                    if f2[j] = c[i]  then
                        t[i][j] := v;
                    else
                        t[i][j] := 0*v;
                    fi;
                else
                    if f1[j] and W.mueneq[c[j]][c[i]] and not f1[i]  then
                        p := FastKazhdanLusztigPolynomial(W,c[i],c[j]);
                        t[i][j] := p[Length(p)] * v;
                    else
                        t[i][j] := 0*v;
                    fi;
                fi;
            od;
        od;
        Add( rep, t );
    od;
    return rep;
end;

LeftCell := function(W,el)
  return First([1..Length(W.leftcells)],i->el in W.leftcells[i]);
end;

NrWord := function(W,w)
  return Position(W.words,ReducedCoxeterWord(W,w));
end;

Cox := function(type,rank)
  # Does all preparations to work:
  local W;
  W := CoxeterGroup(type,rank);
  Read(Concatenation(CoxeterDirectory,"/",type,String(rank),".gap"));
  StoreKazhdanLusztigPolsFromCoxeter(W);
  RLCells(W);
  W.plc := PosetLeftCells(W);
  CalcLevels(W.plc);
  return W;
end;

IdealsOfPoset := function(p)
  # Calculates all ideals of a poset p as delivered by `PosetFromLeftCells'.
  # A list of generator lists is returned.
  local allones,allzeros,els,elssort,g,gens,i,j,l,ll,new,princs;
  princs := List([1..p.n],i->BlistList([1..p.n],p.leq[i]));
  gens := [];
  els := [];
  elssort := [];
  allzeros := BlistList([1..p.n],[]);
  allones := BlistList([1..p.n],[1..p.n]);
  for i in [1..p.n] do
    Add(gens,[i]);
    Add(els,princs[i]);
    AddSet(elssort,princs[i]);
    for j in [i+1..p.n] do
      new := UnionBlist(princs[i],princs[j]);
      if not(new in elssort) then
        Add(gens,[i,j]);
        Add(els,new);
        AddSet(elssort,new);
      fi;
    od;
    Print(".\c");
  od;
  Print("\n");
  i := 1;
  while i <= Length(gens) do
    ll := ShallowCopy(allzeros);
    for g in ListBlist([1..p.n],els[i]) do
      UniteBlist(ll,BlistList([1..p.n],p.mins[g]));
    od;
    SubtractBlist(ll,els[i]);
    for j in ListBlist([1..p.n],ll) do
      new := UnionBlist(els[i],princs[j]);
      if not(new in elssort) then
        l := Filtered(gens[i],x->not(princs[j][x]));
        Add(l,j);
        Add(gens,l);
        Add(els,new);
        AddSet(elssort,new);
      fi;
    od;
    i := i + 1;
    if i mod 100 = 0 then
      Print(i,"(",Length(gens),") \c");
    fi;
  od;
  Print("\nDone.\n");
  return [gens,els];
end;

IdealsOfPosetOld := function(p)
  # Calculates all ideals of a poset p as delivered by `PosetFromLeftCells'.
  # A list of generator lists is returned.
  local allones,els,elssort,g,gens,i,j,l,new,princs;
  princs := List([1..p.n],i->BlistList([1..p.n],p.leq[i]));
  gens := [];
  els := [];
  elssort := [];
  allones := BlistList([1..p.n],[1..p.n]);
  for i in [1..p.n] do
    Add(gens,[i]);
    Add(els,princs[i]);
    AddSet(elssort,princs[i]);
    for j in [i+1..p.n] do
      new := UnionBlist(princs[i],princs[j]);
      if not(new in elssort) then
        Add(gens,[i,j]);
        Add(els,new);
        AddSet(elssort,new);
      fi;
    od;
    Print(".\c");
  od;
  Print("\n");
  i := 1;
  while i <= Length(gens) do
    for j in ListBlist([1..p.n],DifferenceBlist(allones,els[i])) do
        new := UnionBlist(els[i],princs[j]);
        if not(new in elssort) then
          l := Filtered(gens[i],x->not(princs[j][x]));
          Add(l,j);
          Add(gens,l);
          Add(els,new);
          AddSet(elssort,new);
        fi;
    od;
    i := i + 1;
    if i mod 100 = 0 then
      Print(i,"(",Length(gens),") \c");
    fi;
  od;
  Print("\nDone.\n");
  return [gens,els];
end;

IntervalsOfPoset := function(p)
  # Calculates all intervals of a poset as delivered by `PosetFromLeftCells'.
  # Returns a pair, first element is a list of pairs [i,j] and the second
  # element is a list of elements in the interval [i,j] respectively.
  local els,gens,i,j,l;
  gens := [];
  els := [];
  for i in [1..p.n] do
    for j in [1..p.n] do
      if i in p.leq[j] then
        Add(gens,[i,j]);
        l := Intersection(p.leq[j],p.gt[i]);
        AddSet(l,i);
        Add(els,l);
      fi;
    od;
  od;
  return [gens,els];
end;

InduceElements := function(W,w,V,els)
  # V must be a reflection subgroup of W, typically type A_{n-1}. 
  # w must be an isomorphic group with full tuning infrastructure.
  # els is a list of numbers of elements of V as in V.elslenlex.
  # Returns a list of left cells, such that the inclusion of els and
  # their translates is a union of those left cells.
  local allwords,i,j,l,lclist,ll,poss,reps,words;

  # Calculate the inclusion of elements in els:
  words := List(els,x->CoxeterWord(V,V.elslenlex[x]));
  allwords := List(W.elslenlex,x->CoxeterWord(W,x));
  poss := List(words,x->Position(allwords,x));

  # Now the distinguished left coset representatives:
  reps := List(ReducedRightCosetRepresentatives(W,w),x->x^-1);

  # Now we calculate all products:
  l := [];
  for i in reps do
    for j in W.elslenlex{poss} do
      AddSet(l,W.elslookup[Position(W.elements,i*j)]);
    od;
  od;

  # Now calculate intersections with left cells:
  lclist := [];
  for i in [1..Length(W.leftcells)] do
    ll := Intersection(l,W.leftcells[i]);
    if ll <> [] then
      Add(lclist,i);
      if ll <> W.leftcells[i] then
        Error("!!! Alert !!! Left cell no. ",i," only partially contained!\n");
      fi;
    fi;
  od;
      
  return lclist;
end;

RobinsonSchensted := function(n,pi)
  local P,Q,dummy,i,j,pi,pos,z;
  if IsPerm(pi) then
    pi := ListPerm(pi);
  fi;
  while Length(pi) < n do
    Add(pi,Length(pi)+1);
  od;
  P := [];
  Q := [];
  for i in [1..Length(pi)] do
    j := pi[i];
    z := 1;
    while j <> 0 do
      if not(IsBound(P[z])) then
        P[z] := [j];
        Q[z] := [i];
        j := 0;
      else
        pos := PositionSorted(P[z],j);
        if pos > Length(P[z]) then
          Add(P[z],j);
          Add(Q[z],i);
          j := 0;
        else
          dummy := P[z][pos];
          P[z][pos] := j;
          j := dummy;
          z := z + 1;
        fi;
      fi;
    od;
  od;
  return [P,Q];
end;

CompositionDominates := function(p,q)
  local i,m,sp,sq;
  if Sum(p) <> Sum(q) then 
    Error("Comparing compositions of different sum!");
  fi;
  sp := 0;
  sq := 0;
  i := 1;
  m := Minimum(Length(p),Length(q));
  while i <= m do
    sp := sp + p[i];
    sq := sq + q[i];
    if sp < sq then
      return false;
    fi;
    i := i + 1;
  od;
  return true;
end;
      

##############################################################################
#F StandardTableaux 
## Robinson-Schensted.

# Liste standard Tableaux
StandardTableaux := function( n )
    local i, j, z, x, t, T, p, P;
    P := Partitions( n );
    T := [];
    for p in P do
        t := [];
        z := 1;
        for i in p do
            x := [];
            for j in [ 1 .. i ] do
              Add( x, z );
              z := z + 1;
            od;
            Add( t, x);
        od;
        Add( T, t );
    od;
    return T;
end;

#########################################################################
#F PermutationsTableau 
# Berechne zu gegebenem Tableau Urbild unter Robinson-Schensted-Abbildung

PermutationsTableau := function ( tbl )
    local   prms, prm, len, n, m, r, l, t, s;

    # remembering the lengths makes everything easier
    n := 0;
    len := [];
    for r  in [1..Length(tbl)]  do
        n := n + Length(tbl[r]);
        len[r] := Length(tbl[r]);
    od;
    for r  in [Length(tbl)+1..n+1]  do
        len[r] := 0;
    od;

    # output
    prm := [];
    prms := [];

    # loop variables
    m := n;
    r := m;

    # loop until done
    while 0 < r  or m < n  do

        # if there is a hook in line <r>
        if 0 < r  and len[r+1] < len[r]  then

            # take this hook from the tableau
            l := len[r];
            t := tbl[r][l];
            len[r] := l - 1;
            r := r - 1;

            # pop <t> to the correct position in the previous row
            while 0 < r  do
                while l < len[r]  and tbl[r][l+1] < t  do
                    l := l + 1;
                od;
                s := tbl[r][l];
                tbl[r][l] := t;
                t := s;
                r := r - 1;
            od;

            # insert <t> in the permutation
            prm[m] := t;
            m := m - 1;

            # start with new choices
            r := m;

        # if there is no hook but further choices
        elif 0 < r  then

            # go for the next choice
            r := r - 1;

        # else push
        elif 0 < m  then

            # take <t> from the permutation
            m := m + 1;
            t := prm[m];
            r := 1;
            l := m;

            # push <t> to the correct position in the next row
            while 0 < len[r]  and t < tbl[r][len[r]]  do
                if len[r] < l  then l := len[r];  fi;
                while 1 < l  and t < tbl[r][l-1]  do
                    l := l - 1;
                od;
                s := tbl[r][l];
                tbl[r][l] := t;
                t := s;
                r := r + 1;
            od;

            # insert <t> as hook into the tableau
            len[r] := len[r] + 1;
            tbl[r][len[r]] := t;
            r := r - 1;

        # no more choices
        else

            # we have a permutation
            Add( prms, ShallowCopy(prm) );
            m := m + 1;
            len[1] := 1;

        fi;

    od;

    return prms;
end;

MakePSchlange := function(t)
  local n,p;
  p := PermutationsTableau(t)[1];
  n := Maximum(p);
  p := (n+1)-p;
  return RobinsonSchensted(n,p)[1];
end;

STableauOurDominates := function(l,m)
  local i,j,l,ll,m,mm,n,p,q,tabl,tabm;
  n := Sum(List(l,Length));
  if Sum(List(m,Length)) <> n then
    Error("Comparing Tableau of different sizes!");
  fi;
  p := List(l,Length);
  q := List(m,Length);
  if not(CompositionDominates(List(l,Length),List(m,Length))) then
    return false;
  fi;
  l := Copy(l);
  m := Copy(m);
  ll := MakePSchlange(l);
  mm := MakePSchlange(m);
  tabl := [];
  tabm := [];
  for i in [1..Length(l)] do 
    for j in l[i] do
      tabl[j] := i;
    od;
  od;
  for i in [1..Length(m)] do
    for j in m[i] do
      tabm[j] := i;
    od;
  od;
  Unbind(l[tabl[n]][Length(l[tabl[n]])]);
  Unbind(m[tabm[n]][Length(m[tabm[n]])]);
  if not(STableauOurDominates(l,m)) then
    return false;
  fi;
  tabl := [];
  tabm := [];
  for i in [1..Length(ll)] do 
    for j in ll[i] do
      tabl[j] := i;
    od;
  od;
  for i in [1..Length(mm)] do
    for j in mm[i] do
      tabm[j] := i;
    od;
  od;
  Unbind(ll[tabl[n]][Length(ll[tabl[n]])]);
  Unbind(mm[tabm[n]][Length(mm[tabm[n]])]);
  if not(STableauOurDominates(ll,mm)) then
    return false;
  fi;
  return true;
end;

STableauDominates := function(l,m)
  local i,j,n,p,q,tabl,tabm;
  n := Sum(List(l,Length));
  if Sum(List(m,Length)) <> n then
    Error("Comparing Tableau of different sizes!");
  fi;
  p := List(l,Length);
  q := List(m,Length);
  if not(CompositionDominates(p,q)) then
    return false;
  fi;
  tabl := [];
  tabm := [];
  for i in [1..Length(l)] do 
    for j in l[i] do
      tabl[j] := i;
    od;
  od;
  for i in [1..Length(m)] do
    for j in m[i] do
      tabm[j] := i;
    od;
  od;
  for i in [n,n-1..2] do
    p[tabl[i]] := p[tabl[i]]-1;
    q[tabm[i]] := q[tabm[i]]-1;
    if not(CompositionDominates(p,q)) then
      return false;
    fi;
  od;
  return true;
end;

Theory := function(l,m)
  local i,j,n,tabl,tabm;
  n := Sum(List(l,Length));
  if Sum(List(m,Length)) <> n then
    Error("Comparing Tableau of different sizes!");
  fi;
  tabl := [];
  tabm := [];
  for i in [1..Length(l)] do 
    for j in [1..Length(l[i])] do
      tabl[l[i][j]] := j;
    od;
  od;
  for i in [1..Length(m)] do
    for j in [1..Length(m[i])] do
      tabm[m[i][j]] := j;
    od;
  od;
  for i in [1..n] do
    if tabl[i] > tabm[i] then
      return false;
    fi;
  od;
  return true;
end;

LeftCellInclusion := function(W,w)
  # w must be a parabolic subgroup of W, both from "Cox", calculates
  # the inclusion map of left cells.
  local c,d,incl,res,x,y;
  incl := List(w.words,x->Position(W.words,x));
  res := [];
  for c in [1..Length(w.leftcells)] do
    x := w.leftcells[c][1];
    y := incl[x];
    d := First([1..Length(W.leftcells)],i->y in W.leftcells[i]);
    Add(res,d);
  od;
  return res;
end;

VersuchLCP := function(w)
  local cells,celltab,i,j,k,l,ll,n,s;
  l := [];
  n := Size(w);
  celltab := List([1..n],i->[i]);
  for i in [1..n] do
    for s in [1..w.semisimpleRank] do
      if not(w.ldescent[i][s]) then
        j := w.lmulttab[i][s];
        if not(IsSubsetBlist(w.ldescent[j],w.ldescent[i])) then
          if celltab[i] <> celltab[j] then
            UniteSet(celltab[i],celltab[j]);
            ll := celltab[j];
            for k in ll do
              celltab[k] := celltab[i];
            od;
          fi;
        else
          Add(l,[j,i]);
        fi;
      fi;
    od;
  od;
  cells := Set(celltab);
  celltab := [1..n];
  for i in [1..Length(cells)] do
    for j in cells[i] do
      celltab[j] := i;
    od;
  od;
  # Now we have a list of pairs of elements, induce pairs of cells:
  ll := [];
  for i in l do
    Add(ll,[celltab[i[1]],celltab[i[2]]]);
  od;
  return rec(cells := cells,pairs := ll,l := l);
end;

DoRightStar := function(g,y,i)
  local ys,j;
  j := i+1;
  if g.type[1][1] = "D" and i = 1 then
    j := 3;
  fi;
  if (g.rdescent[y][i] and g.rdescent[y][j]) or
     (not(g.rdescent[y][i]) and not(g.rdescent[y][j])) then
    return false;   # not defined on y
  fi;
  ys := g.invtab[g.lmulttab[g.invtab[y]][i]];
  if (g.rdescent[ys][i] and g.rdescent[ys][j]) or
     (not(g.rdescent[ys][i]) and not(g.rdescent[ys][j])) then
    ys := g.invtab[g.lmulttab[g.invtab[y]][j]];
  fi;
  return ys;
end; 

ApplyAllRightStarsOnPairs := function(W,l)
  # l must be a list of pairs [y,w] with y <=_L w. Calculates all pairs
  # which can be achieved by applying right star operations.
  local i,k,ll,p,q;
  ll := Set(l);
  k := 1;
  while k <= Length(l) do
    p := ll[k];
    for i in [1..W.semisimpleRank-1] do
      q := [DoRightStar(W,p[1],i),DoRightStar(W,p[2],i)];
      if IsInt(q[1]) and IsInt(q[2]) and not(q in ll) then
        Add(l,q);
        AddSet(ll,q);
      fi;
    od;
    if k mod 1000 = 0 then
      Print(k,"(",Length(l),")\n");
    fi;
    k := k + 1;
  od;
  return ll;
end;

CellTab := function(W)
  local c,celltab,i;
  celltab := [];
  for c in [1..Length(W.leftcells)] do
    for i in W.leftcells[c] do
      celltab[i] := c;
    od;
  od;
  return celltab;
end;

MakeLeftAscentInCell := function(W,cellnr)
  local celltab,el,i,j,n;
  n := W.semisimpleRank;
  W.celltab := CellTab(W);
  for el in W.leftcells[cellnr] do
    Print("\nWord: ",W.words[el],"\n");
    for i in [1..n] do
      if not(W.ldescent[el][i]) then
        Print(" apply gen. ",i," from the left:\n");
        Print("  ",W.words[W.lmulttab[el][i]]," (cell ",
              W.celltab[W.lmulttab[el][i]],")\n");
        for j in ListBlist([1..Size(W)],W.bruhat[el]) do
          if W.mueneq[el][j] and W.ldescent[j][i] then
            Print("  ",W.words[j]," (cell ",W.celltab[j],")\n");
          fi;
        od;
        Print("\n");
      fi;
    od;
  od;
  return;
end;

TestTheory101 := function(W)
  local c,cbl,cells,co,i,j;
  co := Combinations([1..W.semisimpleRank]);
  for c in co do
    cbl := BlistList([1..W.semisimpleRank],c);
    cells := Filtered([1..Length(W.leftcells)],
                      i-> cbl = W.rdescent[W.leftcells[i][1]]);
    for i in cells do
      for j in cells do
        if (j in W.plc.leq[i]) <>
           (W.bruhat[W.duflos[j]][W.duflos[i]]) then
          Print("Passt nicht bei Zellen ",i," und ",j,"!\n");
        fi;
        if not(j in W.plc.leq[i] or i in W.plc.leq[j]) then
          Print("Alarm: Keine totale Ordnung bei ",i," und ",j,"!\n");
        fi;
        if not(W.bruhat[W.duflos[j]][W.duflos[i]] or
               W.bruhat[W.duflos[i]][W.duflos[j]]) then
          Print("Alarm: Keine totale Ordnung bei Duflos bei ",i," und ",j,
                "!\n");
        fi;
      od;
    od;
  od;
end;

##########################################
# from here on old stuff:

PutCoxeterInfoIntoCoxeterGroup := function(W,p)
  # p is a record as coming from "coxeter" and "coxtogap.py", containing
  # all KL-Polynomials of "extremal" pairs. w is a CoxeterGroup of the
  # correct type.
  local 2N,l,lw,otc1,pw,py,tmp,tup,w,wcr,y,z;

  # Set up the infra-structure in W:
  if not IsBound(W.criticalPairs) then
    W.criticalPairs:=List([1..W.parentN+1],x->[]);
    W.klpol:=[];
  fi;
  
  # Now go through pairs:
  l := CoxeterElementsLengthLex(W);
  for tup in p.pairs do
    y := l[tup[2]];
    w := l[tup[1]];
    lw := CoxeterLength(W,w);
    # We now know P_{y,w} = p.pols[tup[3]] as list of coefficients
    # and y,w is a critical pair
    pw:=Position(CoxeterElementsLength(W,lw),w);
    wcr:=W.criticalPairs[lw+1];

    # encode y in some way to save memory (an element is uniquely determined
    # by its images of the simple roots)
    2N:=2*W.parentN;
    tmp:=OnTuples(W.rootInclusion{[1..W.semisimpleRank]},y);
    otc1:=tmp[1];
    for z in [2..Length(tmp)] do
      otc1:=otc1*2N+tmp[z];
    od;
    if not(IsBound(wcr[pw])) then
      wcr[pw]:=[[],[]];
    fi;
    py:=Position(wcr[pw][1],otc1);
    if py=false then   # otherwise nothing to do!
      Add(wcr[pw][1],otc1);
      py:=Position(W.klpol,p.pols[tup[3]]);
      if py=false then
        Add(W.klpol,p.pols[tup[3]]);
        py:=Length(W.klpol);
      fi;
      Add(W.criticalPairs[lw+1][pw][2],py);
    else # do a check:
      if W.klpol[wcr[pw][2][py]] <> p.pols[tup[3]] then
        Print("Alarm: Polynomial wrong! ",tup,"\n");
      fi;
    fi;
  od;
end;

ExtremalPairPlus := function(W,y,w)
  # Calculates an extremal pair for y,w, given as numbers in W.elslenlex.
  # `TuneCoxeterGroup' must have been called. Only feasible for y <= w.
  local s;
  if not(W.bruhat[w][y]) then
    return false;
  fi;
  if W.invtab[w] < w then
    w := W.invtab[w];
    y := W.invtab[y];
  fi;
  if y = w then return [w,w]; fi;
  s := 1;
  while s <= W.semisimpleRank do
    if s = 1 and W.invtab[w] < w then
      w := W.invtab[w];
      y := W.invtab[y];
    fi;
    if W.ldescent[w][s] and not(W.ldescent[y][s]) then
      if not(W.bruhat[W.lmulttab[w][s]][y]) then
        w := W.lmulttab[w][s];
      fi;
      y := W.lmulttab[y][s];
      if y = w then return [w,w]; fi;
      s := 0;   # start again with first generator
    elif W.rdescent[w][s] and not(W.rdescent[y][s]) then
      if not(W.bruhat[W.invtab[W.lmulttab[W.invtab[w]][s]]][y]) then
        # this either shortens the length difference by 2 or does a shift!
        w := W.invtab[W.lmulttab[W.invtab[w]][s]];  # see below
      fi;
      # the following would be   W.rmulttab[y][s]
      # if we had an rmulttab:
      y := W.invtab[W.lmulttab[W.invtab[y]][s]];
      if y = w then return [w,w]; fi;
      s := 0;   # start again with first generator
    fi;
    s := s + 1;
  od;
  return [y,w];
end;

OldExtremalPair := function(W,y,w)
  # Calculates an extremal pair for y,w, given as numbers in W.elslenlex.
  # `TuneCoxeterGroup' must have been called.
  local s;
  if not(W.bruhat[w][y]) then
    return false;
  fi;
  if W.invtab[w] < w then
    w := W.invtab[w];
    y := W.invtab[y];
  fi;
      #elif W.invtab[w] = w then
      #  if W.invtab[y] > y then
      #    y := W.invtab[y];
      #  fi;
  s := 1;
  while s <= W.semisimpleRank do
    # first try on the left hand side:
    if W.ldescent[w][s] and
      #then   # only if s makes w shorter
      #if not(W.bruhat[W.lmulttab[w][s]][y]) then
      #  # this either shortens the length difference by 2 or does a shift!
      #  w := W.lmulttab[w][s];
      #  y := W.lmulttab[y][s];
      #  s := 0;  # begin anew, s is incremented at the end of the loop!
      #  if w=y then return [w,w]; fi;
      #elif 
      not(W.ldescent[y][s]) then
        # we can make y longer to shorten the interval
        y := W.lmulttab[y][s];
        if w=y then return [w,w]; fi;
        s := 0;   # begin anew, s is incremented at the end of the loop!
        # otherwise s does nothing good for us
      #fi;
    elif W.rdescent[w][s] and
      # then  # only if s makes w shorter
      # the following would be
      #  W.bruhat[W.rmulttab[w][s]]
      # if we had an rmulttab!
      #if not(W.bruhat[W.invtab[W.lmulttab[W.invtab[w]][s]]][y]) then
      #  # this either shortens the length difference by 2 or does a shift!
      #  w := W.invtab[W.lmulttab[W.invtab[w]][s]];  # see above
      #  y := W.invtab[W.lmulttab[W.invtab[y]][s]];  # see above
      #  s := 0;  # begin anew, s is incremented at the end of the loop!
      #  if w=y then return [w,w]; fi;
      #elif 
      not(W.rdescent[y][s]) then
        # we can make y longer to shorten the interval
        y := W.invtab[W.lmulttab[W.invtab[y]][s]];
        if w=y then return [w,w]; fi;
        s := 0;   # begin anew, s is incremented at the end of the loop!
        # otherwise s does nothing good for us
      #fi;
    fi;
    s := s + 1;
  od;
  # now the pair is reduced
  return [y,w];
end;

MyMue := function(W,x,y)
  # Returns mue(x,y) by first trying to lookup and then calling
  # x and y are numbers of elements in W
  local lx,ly,w0,mue;
  if W.muetab[x][y] then
    return W.muevals[x][Position(W.mueels[x],y)];
  fi;
  lx := Length(W.words[x]);
  ly := Length(W.words[y]);
  if lx + ly > W.N  then
     w0 := LongestCoxeterElement(W);
     mue := KazhdanLusztigMue(W,w0*W.elements[y],w0*W.elements[x],
                                W.N-ly,W.N-lx);
  else 
     mue := KazhdanLusztigMue(W,W.elements[x],W.elements[y],lx,ly);
  fi;
  W.muetab[x][y] := true;
  Add(W.mueels[x],y);
  Add(W.muevals[x],mue);
  SortParallel(W.mueels[x],W.muevals[x]);
  return mue;
end;
  
# Merken: http://www.math.hawaii.edu/~ralph/LatDraw/
