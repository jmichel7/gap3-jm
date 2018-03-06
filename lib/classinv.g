#############################################################################
##
#A  classinv.g               CHEVIE library                     Frank Luebeck
##
#A  $Id: classinv.g,v 1.1 1997/01/21 13:46:21 gap Exp $
##  
#Y  Copyright (C) 1992 - 1996  Lehrstuhl D fuer Mathematik, RWTH Aachen, IWR
#Y  der Universitat Heidelberg, University of St. Andrews, and   University 
#Y  Paris VII.
##
##  This file contains  general  function   for computing invariants    of
##  conjugacy classes  in   permutation   groups.  
##
#H  08/12/97  Frank Luebeck
#H  made PositionClass work for one element cosets.

#############################################################################
##
#F  ClassInvariants( <g>[, <func1>[, ....]..] ) . . . . . . . class invariants
#F  for the group <g> (dispatcher function)
#F  PositionClass( <g>, <x> ) . . . . . . . . . . . . . class of given element
#F  (dispatcher function)
##  
##  <g> should be a domain such that 'ConjugacyClasses(<g>)' is defined.
##  
##  'PositionClass( <g>, <x> )' with <x> in <g> should return:  
##          Position(ConjugacyClasses(<g>), ConjugacyClass(<g>, <x>))
##  
##  'ClassInvariants'  returns a record,  which is used by 'PositionClass'
##  to identify the conjugacy class of a given element in <g>.
##  
##  The  optional functions    <func1>,      ... in  the    argument    of
##  'ClassInvariants' should have the following properties:
##  
##    - their argument is of form <invrec>, <l> 
##      where <invrec> is a partial result record of 'ClassInvariants'
##      and <l> is a sublist of [1..Length(ConjugacyClasses(<g>)]
##    - they return a function f of the following form:
##      - the argument of f is an element x of <g>
##      - f returns an invariant of the conjugacy class of x
##  
##  The  easiest way to understand  this may be  to  look at the functions
##  'PermGroupClassInvariants' and 'CentreMultFunction' below.
##  
PositionClass:=function(g,x)
  if IsBound(g.operations.PositionClass) then
    return g.operations.PositionClass(g,x);
  else return PositionProperty(ConjugacyClasses(g),y->x in y);
  fi;
end;

#############################################################################
##
#F  CycleStructurePermOrbits( <perm>, <orbs> ) . . . . . . . cycle structures
#F  of <perm> on orbits <orbs>
##  
##  <orbs> must be a list of  lists of points containing the orbits of the
##  permutation <perm>. Then CycleStructurePermOrbits( <perm>, <orbs> ) is
##  faster code for:
##        List(orbs,orb->CycleStructurePerm(RestrictedPerm(perm,orb)))
##  
CycleStructurePermOrbits:=function(perm,orb)local cysl,cys,mark,i,j,len,p;
  if perm = ()  then return [];fi;
  mark := BlistList([],[]); cysl := [  ];
  for j in orb do
    cys:=[];
    for i  in j  do
      if not IsBound(mark[i]) then
	len:=0; mark[i]:=true; p:=i^perm;
	while not IsBound(mark[p]) do len:=len+1;mark[p]:=true;p:=p^perm;od;
	if 0 < len  then
	  if IsBound( cys[len] )  then cys[len] := cys[len] + 1;
	  else cys[len] := 1;
	  fi;
	fi;
      fi;
    od;
    Add(cysl,cys);
  od;  
  return cysl;
end;

#############################################################################
##
#F  ShortClassListFunction( <r>, <l> ) . . . . . . returns function checking
#F  membership in classes via lists of elements
##  
##  This        function     can    be     used     as    argument      of
##  'PermGroupOps.ClassInvariants'.   It   checks if   the    classes   of
##  <r>.reps{<l>} have size at  most MAXSHORTCLASSLIST (default 200) or if
##  their elements  are already   computed. In this    case it  returns  a
##  function f which returns f(x)=p if x  is conjugate to <r>.reps{<l>}[p]
##  (provided x is conjugate to one of <r>.reps{<l>}).  Here f  checks for 
##  membership in the lists of elements (using binary search).
##  
##  If the classes corresponding  to <l> are  too large and their elements
##  are not known, then a constant function is returned (which immediately
##  is found not to be able to distinguish classes).
##  
MAXSHORTCLASSLIST:=200;

# better, because no new list of elements is stored in .help:
ShortClassListFunction:=function(r,l)local h, a, b, x;
  
  # checking sizes of classes or if elements already computed:
  if ForAll(l,x->Size(ConjugacyClasses(r.g)[x])<=MAXSHORTCLASSLIST
            or IsBound(ConjugacyClasses(r.g)[x].elements)) then
    h:=Length(r.help)+1;
    if Length(l)>2 then
      # don't need to check the last one:
      r.help[h]:=List(l{[1..Length(l)-1]},
                      x->Elements(ConjugacyClasses(r.g)[x]));

      return function(el)local p; p:=1;
	  while p<=Length(r.help[h]) and not el in r.help[h][p] do p:=p+1;od;
	  return p;
	end;
    else
      r.help[h]:=Elements(ConjugacyClasses(r.g)[l[1]]);
      return el->el in r.help[h];
    fi;
  else
    # trivial function, if classes too large:
    return el->0;
  fi;
end;

#############################################################################
##
#F  CentreMultFunction( <r>, <l> ) . . . . . . . . . . returns function
#F  which returns cycle types of elements multiplied by center elements
##  
##  This   can be  used   as argument  for 'PermGroupOps.ClassInvariants'.
##  'CentreMultFunction' returns a function f.
##  
##  If it  can be used  to distinguish classes  corresponding to  <l> then
##  f(<x>) returns a list of cycle types  of <x>*z where  z runs through a
##  set Z' of center elements of <r>.g. The function tries to minimize the
##  set Z'. 
##  
##  If no such cycle type distinguishes classes corresponding to  <l> then
##  f is a constant function.
##  
CentreMultFunction:=function(r,l)
  local c, tmp, pos, h, i, len, len1;
  
  # cycle structures for elements multiplied by non trivial center
  # elements:
  c:=Elements(Centre(r.g));
  c:=c{[2..Length(c)]};
  if IsBound(r.orbits) then
    tmp:=List(r.reps{l},x->List(c,y->CycleStructurePermOrbits(x*y,r.orbits)));
  else
    tmp:=List(r.reps{l},x->List(c,y->CycleStructurePerm(x*y)));
  fi;
  
  len:=Length(Set(tmp));
  # check,  if useful to distinguish classes:
  if len>1 then
    # throwing out unnecessary center elements 
    len1:=1;
    pos:=[];
    for i in [1..Length(c)] do
      if len1<len and 
         Length(Set(tmp{[1..Length(tmp)]}{Concatenation(pos,[i])}))>len1 then
        Add(pos,i);
      fi;
    od;
    h:=Length(r.help)+1;
    r.help[h]:=c{pos};
    
    # returning the function:
    if IsBound(r.orbits) then
      return function(x) return List(r.help[h],y->
                             CycleStructurePermOrbits(y*x,r.orbits));end;
    else
      return function(x) return List(r.help[h],y->
                                     CycleStructurePerm(x*y));end;
    fi;
  else
    # trivial function if not useful
    return el->0;
  fi;
end;


#############################################################################
##
#F  PointsRepOrb( <g> ) . . . . . orbits of points under permutation group <g>
#F  and transversal
##  
##  <g>  must be  a  permutation  group.   'PointsRepOrb'  returns a  list
##  [orb,rep].   Here    orb is    a     list   of the    <g>-orbits    on
##  [1..LargestMovedPoint(<g>)]. rep[i][j]  is    an x in  <g>   such that
##  orb[i][j]^x=orb[i][1].
##  
# args: G[, m] (m maximal  point)
PointsRepOrb:=function ( arg )
    local  G, orbs, orb, max, g, gs, new, gen, Ggen, p, pnt, fst, img;
    G:=arg[1];

    if Length(arg)>1 then
      max:=arg[2];
    else
      max := PermGroupOps.LargestMovedPoint(G);
    fi;

    Ggen:=G.generators;
    new := BlistList( [ 1 .. max ], [ 1 .. max ] );
    orbs := [  ];
    gs:=[];
    fst := 1;
    while fst <> false  do
	orb := [ fst ];
	g:=[()];
	new[fst] := false;
	p:=1;
	while p<=Length(orb)  do
	    for gen  in Ggen  do
		img := orb[p] ^ gen;
		if new[img]  then
		    Add( orb, img );
		    Add(g,gen mod g[p]);
		    new[img] := false;
		fi;
	    od;
            p:=p+1;
	od;
	Add( orbs, orb );
	Add( gs, g );
	fst := Position( new, true, fst );
    od;
    return [orbs,gs];
end;


#############################################################################
##
#F  FingerPrintFunction( <r>, <l> ) . . . . . . returns function which computes
#F  finger print of elements 
##  
##  If  <r>.g  moves   at  most MAXNUMBERPOINTSFINGERPRINT (default   200)
##  points, then a function f is  returned which computes the fingerprints
##  of cycles of fixed length of an element.
##  
##  (Else a constant function is returned.)
##  
##  This is a refinement of the *fingerprint* used in the programs for the
##  Dixon-Schneider algorithm:
##  
##  Number the orbits   of <r>.g on the set   of pairs of moved  points of
##  <r>.g by 1..k. For a given element x  in <r>.g the fingerprint f(x) is
##  a list   whose  i-th  entry tells   for   each l in [1..k]    how many
##  (i+1)-cycles [p0, p1, .., pi] of x have the following property:
##  
##    [p0,p1] is contained in orbit l of pairs of points
##  
##  (Obviously l does not depend on the choice  of the pair [p0,p1] in the
##  cycle.)
##  
##  This is clearly a class  invariant, since any y in  <r>.g maps a cycle
##  via conjugation onto
##  
##    [p0, p1, .., pi]^y = [p0^y, p1^y, .., pi^y]
##  
MAXNUMBERPOINTSFINGERPRINT:=200;
FingerPrintFunction:=function( r, l)
  local po, m, k, h, sp, p, a, b, i, j, S, orbs, tmp;

  m:=PermGroupOps.LargestMovedPoint(r.g);

  # orbits of points and transversals:
  po:=PointsRepOrb(r.g);
  po[1]:=Filtered(po[1],a->Length(a)>1);
  po[2]:=Filtered(po[2],a->Length(a)>1);

  # if there are too many moved points:
  if Sum(List(po[1],Length))>MAXNUMBERPOINTSFINGERPRINT then
    return x->0;
  fi;

  # elements mapping p to smallest point in orbit
  sp:=[];
  for a in [1..Length(po[1])] do
    for i in [1..Length(po[1][a])] do
      sp[po[1][a][i]]:=po[2][a][i];
    od;
  od;

  # determine different orbits of pairs of points:
  orbs:=[];
  k:=1;
  for a in [1..Length(po[1])] do
    p:=po[1][a][1];
    S:=Stabilizer(r.g,p);
    tmp:=Orbits(S,[1..m]);
    orbs[p]:=[];
    for b in tmp do
      if IsBound(sp[b[1]]) then
        for j in b do 
          orbs[p][j]:=k;
        od;
        k:=k+1;
      fi;
    od;
  od;

  # not usable if r.g is two times transitive on moved points: 
  if k=2 then
    return x->0;
  fi;

  h:=Length(r.help)+1;
  r.help[h]:=sp;
  r.help[h+1]:=orbs;
  
  # return a function which gives list of fingerprints 
  # on cycles with fixed length:
  return 
      function(x)
	local k, y, cys, degree, mark, i, j, len, cyc;
        if x=() then
          return [];
        fi;
	degree := LargestMovedPointPerm( x );
	mark := BlistList( [ 1 .. degree ], [  ] );
	cys := [  ];
	for i  in [ 1 .. degree ]  do
	    if not mark[i]  then
		cyc := CyclePermInt( x, i );
		len := Length( cyc ) - 1;
		if 0 < len  then
		    y:=r.help[h][i];
		    k:=r.help[h+1][i^y][(i^x)^y];
		    if IsBound( cys[len] )  then
			Add(cys[len],k);
		    else
			cys[len] := [k];
		    fi;
		fi;
		for j  in cyc  do
		    mark[j] := true;
		od;
	    fi;
	od;
	for i in [1..Length(cys)] do
	  if IsBound(cys[i]) then
	    cys[i]:=Collected(cys[i]);
	  fi;
	od;
	return cys;
      end;
end;

#############################################################################
##
#F  ConjugacyTestFunction( <r>, <l> ) . . . . . . . . . returns a function 
#F  computing an invariant using backtrack
##  
##  This   can  be used    as argument for 'PermGroupOps.ClassInvariants'.
##  'ConjugacyTestFunction' returns a function f. If x  is conjugate to an
##  element of <r>.reps{<l>} then f(x) returns  p such that x is conjugate
##  to <r>.reps{<l>}[p].  p  is computed  by using  the 'in' operator  for
##  conjugacy classes.
##  
ConjugacyTestFunction:=function(r,l)
  local a, h, inClass;
  
  h:=Length(r.help)+1;
  r.help[h]:=ConjugacyClasses(r.g){l{[1..Length(l)-1]}};
  # compute centralizers, if necessary:
  for a in r.help[h] do
    if not IsBound(a.elements) and not IsBound(a.centralizer) then
      a.centralizer:=Centralizer(a.group,a.representative);
    fi;
  od;
  
  # The usual 'in' function for conjugacy classes checks, if
  # x is an element of a.group. We want to avoid this for using
  # it in cosets:
  inClass:=function(a,x)
    if IsBound(a.elements) then return x in a.elements; fi;
    return a.group.operations.RepresentativeConjugationElements(a.group,
                   x,a.representative,a.centralizer)<>false;
  end;
  
  # the function increases a counter until the class is found:
  return function(x) local n; n:=1; while n<=Length(r.help[h]) and 
    not inClass(r.help[h][n],x) do n:=n+1; od; return n;end;
end;




#############################################################################
##
#F  PermGroupOps.ClassInvariants( <g>,  ... ) . . . . . 'ClassInvariants' for 
#F  permutation groups
#F  PermGroupOps.PositionClass( <g>, <x> ) . . . . . . . 'PositionClass' for 
#F  permutation groups
##  
##  
##  This implements  'ClassInvariants' and 'PositionClass' for permutation
##  groups.
##  
##  If  'PermGroupOps.ClassInvariants'   is called  without  the  optional
##  functions then several  default   functions are  used.  These   return
##  functions which compute similar invariants as used in the programs for
##  the Dixon-Schneider algorithm.
##  
##  Let <funlist> be  the list of such functions.
##  
##  'PermGroupOps.ClassInvariants' works  as follows: It always  starts to
##  compute the cycle types of elements on the orbits of <g> on the points
##  moved by <g>.
##  
##  This creates a record r with components:
##  
##    .g       <g>
##    .help    just an empty list which may be used by the functions
##             computing invariants
##    .reps    representatives of the conjugacy classes
##    .inv     cycle types corresponding to the elements in .reps
##    .invset  sorted list of invariants (now Set(.inv))
##    .ret     a list of Length(.invset) such that .ret[i] is
##               - either the unique position of .invset[i] in .inv
##               - or the list of positions of .invset[i] in .inv
##    .orbits  non trivial orbits of <g> on its moved points, if there
##             is more than one such orbit 
##  Now, as long as there are still lists of positions in .ret it proceeds
##  as follows:
##   - take the next function F from <funlist>
##   - for each list of positions l in .ret do
##       - f:=F( r, l);
##       - apply f to elements in r.reps{l}
##       - if this is not a constant function then
##          - substitute l in .ret by the function f
##          - refine invariants of corresponding classes with
##            the values of f.
##          - add these new invariants to the set .invset
##          - update .ret to contain in .ret[i]
##              - either the unique position of .invset[i] in .inv 
##              - or the list of positions of .invset[i] in .inv
##              - or a function which computes an invariant
##  
##  As final function in <funlist> there is always 'ConjugacyTestFunction'
##  which returns a function checking conjugacy  by a backtrack algorithm.
##  This ensures that  finally there are only  positions and functions  in
##  r.ret.
##  
##  Now 'PermGroupOps.PositionClass( <g>, <x> )' works as follows:
##   - it computes the cycle type of <x> as first invariant and determines
##     its position p in ClassInvariants( <g> ).invset
##   - while ClassInvariants( <g> ).ret[p] is a function this is used to
##     refine the invariant and the position p of the new invariant is 
##     determined
##   - finally ClassInvariants( <g> ).ret[p] is returned as result
##  
PermGroupOps.ClassInvariants:=function(arg)
  local g, fs, ff, f, invset, res, ccl, a, b, x, y, i, j, tmp;
  
  g:=arg[1];
  if Length(arg)>1 then
    fs:=arg{[2..Length(arg)]};
  else
    # the default:
    fs:=[CentreMultFunction,ShortClassListFunction,FingerPrintFunction];
  fi;
  # explicit conjugacy tests are always the last chance:
  Add(fs,ConjugacyTestFunction);
  
  # trivial group
  if Size(g)=1 then
    if IsBound(g.generators[1]) then
      x:=g.generators[1];
    else
      x:=Elements(g)[1];
    fi;
    a:=CycleStructurePerm(x);
    return rec(inv:=[a],invset:=[a],ret:=[1],reps:=[x]);
  fi;
  
  res:=rec(g:=g, help:=[]);
  
  res.reps:=List(ConjugacyClasses(g),Representative);
  
  # to work with cosets moving "outer" points:
  if IsBound(g.phi) and g.phi<>() then
    tmp:=[1..Maximum(List(Concatenation(g.generators,[g.phi]),
                                         LargestMovedPointPerm))];
  else
    tmp:=[1..Maximum(List(g.generators,LargestMovedPointPerm))];
  fi;
  res.orbits:=Filtered(PermGroupOps.Orbits(g,tmp,OnPoints),a->Length(a)>1);
  
  if Length(res.orbits)=1 then
    Unbind(res.orbits);
    res.inv:=List(res.reps,x->CycleStructurePerm(x));
  else
    res.inv:=List(res.reps,x->CycleStructurePermOrbits(x,res.orbits));
  fi;  
  res.invset:=Set(res.inv);
  res.ret:=[];
  for x in res.invset do 
    a:=Position(res.inv,x);
    b:=Position(res.inv,x,a);
    if b<>false then
      a:=[a];
      repeat
        Add(a,b);
        b:=Position(res.inv,x,b);
      until b=false;
    fi;
    Add(res.ret,a);
  od;
  
  for ff in fs do
    invset:=[];
    for i in [1..Length(res.ret)] do
      x:=res.ret[i];
      # if there are classes not yet distinguished
      if IsList(x) then
        f:=ff(res,x);
        tmp:=List(res.reps{x},f);
        # if f distinguishes more classes:
        if Length(Set(tmp))>1 then
          res.ret[i]:=f;
          for j in [1..Length(x)] do
            res.inv[x[j]]:=[res.inv[x[j]],tmp[j]];
          od;
          Append(invset,res.inv{x});
        fi;
      fi;
    od;
    invset:=Set(invset);
    Append(res.invset,invset);
    for x in invset do
      a:=Position(res.inv,x);
      b:=Position(res.inv,x,a);
      if b<>false then
        a:=[a];
        repeat
          Add(a,b);
          b:=Position(res.inv,x,b);
        until b=false;
      fi;
      Add(res.ret,a);
    od;
  od;
  
  # sorting res.invset such that we can do binary searches 
  # in 'PositionClass'
  SortParallel(res.invset,res.ret);
    
  return res;
end;

PermGroupOps.PositionClass:=function(g,x)local ci, inv, p;
  
  ci:=ClassInvariants(g);
  
  if IsBound(ci.orbits) then inv:=CycleStructurePermOrbits(x,ci.orbits);
  else                       inv:=CycleStructurePerm(x);
  fi;
  
  p:=PositionSet(ci.invset,inv);
  while IsFunc(ci.ret[p]) do
    inv:=[inv,ci.ret[p](x)];
    p:=PositionSet(ci.invset,inv);
  od;
  return ci.ret[p];
end;


#############################################################################
##
##  In examples   it seems   to be   (sometimes much)   more efficient  to
##  overwrite 'GroupOps.FusionConjugacyClasses'  for   subgroups  by a
##  function which uses 'PositionClass'. But is this true in general?
##  

PermGroupOps.FusionConjugacyClasses:=function(u,g)
  if Parent(u)=Parent(g) then
    return List(List(ConjugacyClasses(u),Representative),x->PositionClass(g,x));
  else
    return GroupOps.FusionConjugacyClasses(u,g);
  fi;
end;
