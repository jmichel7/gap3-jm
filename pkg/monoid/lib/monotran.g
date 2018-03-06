#############################################################################
##
#A  GAP                                                 Goetz.Pfeiffer@UCG.IE
##
#A  $Id: monotran.g,v 2.6 1997/11/05 13:21:34 goetz Exp $
##
#Y  Copyright (C) 1997, Mathematics Dept, University College Galway, Ireland.
##
##  This file defines the functions for  transformation semigroups & monoids.
##
##  The theory behind these algorithms is developed in 
##  
##  [LPRR1] S. A.   Linton, G.  Pfeiffer, E.  F.  Robertson, and N.   Ruskuc,
##  Groups  and actions in  transformation semigroups,  to appear  in Math Z.
##  (1997).
##
##  The algorithms themselves are described in
##
##  [LPRR2] S. A.   Linton, G.  Pfeiffer, E.  F.  Robertson, and N.   Ruskuc,
##  Computing transformation semigroups, (1996), in preparation.
##  
##  Another reference is
##
##  [LM]  G.  Lallement and R.    McFadden, On the   determination of Green's
##  relations in finite transformation semigroups, J. Symbolic Computation 10
##  (1990), 481--489.
##

#############################################################################
##
#F  InfoMono? . . . . . . . . . . . . . . . . . . . . . . . . info functions.
##
if not IsBound(InfoMono1) then InfoMono1:= Ignore; fi;
if not IsBound(InfoMono2) then InfoMono2:= Ignore; fi;

#############################################################################
##
#V  TransSemiGroupOps  . . . . . . . . . . . . . . . . . . operations record.
##
##  This   is    the  operations   record   for   transformation  semigroups.
##  Transformation   semigroups are   a   special kind   of  semigroups, thus
##  'TransSemiGroupOps' inherits from 'SemiGroupOps'.
##
TransSemiGroupOps:= OperationsRecord("TransSemiGroupOps", SemiGroupOps);

#############################################################################
##
#V  TransMonoidOps . . . . . . . . . . . . . . . . . . . . operations record.
##
##  This  is    the   operations record   for    transformation   monoids.  A
##  Transformation monoid is a special kind  of monoid, thus 'TransMonoidOps'
##  inherits from 'MonoidOps'.
##
TransMonoidOps:= OperationsRecord("TransMonoidOps", MonoidOps);

#############################################################################
##
#F  IsTransMonoid( <obj> ) . . . . . . . . . . . . . . . . . . .  type check.
##
IsTransMonoid:= function(obj)

   return IsRec(obj) and IsBound(obj.isTransMonoid) and obj.isTransMonoid;

end;

#############################################################################
##
#F  TransSemiGroupOps.Degree( <S> )  . . . . . . . . . . . . . . . .  degree.
##
##  The *degree* of  a transformation semigroup  is the  number  of points it
##  acts upon.
##
#C  this does not work for the empty semigroup!
##
TransSemiGroupOps.Degree:= S -> Degree(S.1);

#############################################################################
##
#F  TransMonoidOps.Degree( <M> ) . . . . . . . . . . . . . . . . . .  degree.
##
##  The *degree* of a transformation monoid is the  number  of points it acts
##  upon.
##
TransMonoidOps.Degree:= M -> Degree(M.identity);

#############################################################################
##
#F  TransSemiGroupOps.Elements( <S> )  . . . . . . . . . . . . . .  elements.
##
##  The set of elements of a transformation  semigroup can be determined by a
##  simle orbit algorithm.
##
TransSemiGroupOps.Elements:= function(S)

   local gens, elts, orbt, n, w, s;

   gens:= Generators(S);
   orbt:= Set(gens);  elts:= ShallowCopy(orbt);
   for w in orbt do
      for s in gens do
         n:= w * s; 
         if not n in elts then
            Add(orbt, n);
            AddSet(elts, n);
         fi;
      od;
   od;

   return elts;

end;

#############################################################################
##
#F  TransMonoidOps.Elements( <M> ) . . . . . . . . . . . . . . . .  elements.
##
##  The set of elements of a transformation monoid is determined as the union
##  of the element sets of its R classes.
##
TransMonoidOps.Elements:= function(M)

   return Set(Concatenation(List(RClasses(M), Elements)));

end;

#############################################################################
##
#F  <x> in <M> . . . . . . . . . . . . . . . . . . . . . . . membership test.
##
##  A transformation <x>  lies in a transformation monoid  <M> if it  has the
##  same degree as <M> and if it is contained in one of the R classes of <M>.
##
TransMonoidOps.\in:= function(x, M)

   local i, j, k, pos, R, ker;

   # check degree.
   if not IsTransformation(x) or Degree(x) <> Degree(M)  then
      return false;
   fi;

   # unfold monoid, if necessary.
   RClasses(M);

   # check image.
   k:= Rank(x);
   pos:= Position(M.images[k], Image(x));
   if pos = false then
      return false;
   fi;

   # locate representing R class.
   j:= M.imagePos[k][pos];
   R:= M.orbitClasses[j];

   # check kernels.
   ker:= Kernel(x);
   for i in [1..Length(M.kernels[j])] do
       if M.kernels[j][i] = ker then
          if M.lTrans[j][i] * x in R then
             return true;
          fi;
       fi;
   od;

   # if that fails.
   return false;

end;

#############################################################################
##
#F  IsRegularTrans( <M>, <x> ) . . . . . . . . . . . . . .  regularity check.
##
##  A transformation <x> is regular  inside the transformation monoid <M>, if
##  <x> lies insied a regular D class of <M>.
##
##  The transformation <x>  is regular if  the orbit of  $img <x>$ contains a
##  cross section of $ker <x>$.  If $r$ is the rank of  <x> then a set $n$ of
##  size $r$ is a  cross section of $ker <x>$  iff the set of images  $n^<x>$
##  has size $r$. [LM 3.2]
##
IsRegularTrans:= function(M, x)

   local r, orb, gens, s, p, n;

   # initialize.
   n:= Set(x.images);   r:= Size(n);

   # maybe <x> lies inside a group.
   if Size(Set(x.images{n})) = r then
      return true;
   fi;

   # otherwise form the orbit of img x.
   orb:= [n];   gens:= Generators(M);
   for p in orb do
      for s in gens do
         n:= Set(s.images{p});
         if Size(n) = r and not n in orb then

            # did we find a cross section?
            if Size(Set(x.images{n})) = r then
               return true;
            fi;

            Add(orb, n);
         fi;
      od;
   od;

   # if we arrive here, no cross section has been found.
   return false;

end;

#############################################################################
##
#F  OnKernelsAntiAction( <ker>, <s> )  . . . . . . . . . . . . . . .  action.
##
##  This implementation relies heavily on the data structure of
##  transformations!!!
##
OnKernelsAntiAction:= function(ker, s)

   local n, pos, new, loc, i;

   n:= Length(s.images);  # should be 'Degree(s)', but thats too slow :-(
   pos:= [];  new:= [];  loc:= [];

   # construct transformation 'pos' with kernel 'ker'.
   for i in [1..Length(ker)] do
      pos{ker[i]}:= List(ker[i], x-> i);
   od;

   # apply 's' from the left.
   pos:= pos{s.images};

   # determine kernel.
   for i in [1..n] do 
      if IsBound(loc[pos[i]]) then
         Add(new[loc[pos[i]]], i);
      else
         Add(new, [i]);
         loc[pos[i]]:= Length(new);
      fi;
   od;

   # return the kernel.
   return new;

end;

#############################################################################
##
#F  OnTuplesOfSetsAntiAction( <tup>, <s> ) . . . . . . . . . . . . .  action.
##
##  We definitively need something more efficient here!
##
OnTuplesOfSetsAntiAction:= function(tup, s)

   s:= s^-1;
   return List(tup, x-> x^s);
 
end;

#############################################################################
##
#V  TransRClassOps . . . . . . . . . . . . . . . . . . . . operations record.
##
##  This is the  operations record for  R classes  of transformation monoids.
##  It inherits from the general 'RClassOps'.
##
TransRClassOps:= OperationsRecord("TransRClassOps", RClassOps);

#############################################################################
##
#F  TransMonoidOps.RClass( <M>, <x> )  . . . . . . . . . . construct R class.
##  
##  The initial construction of an R class is delegated to 'MonoidOps.RClass'
##  where the R class of <x> in <M> is represented by a record with basically
##  two components, 'representative' for <x>  and 'monoid' for <M>.  Then  we
##  overlay the operations record with  the more efficient (or only existing)
##  functions in 'TransRClassOps'.
##
TransMonoidOps.RClass:= function(M, x)

   local R;

   R:= MonoidOps.RClass(M, x);
   R.operations:= TransRClassOps;

   return R;

end;

#############################################################################
##
#F  TransRClassOps.SchutzenbergerGroup( <R> )  . . . . . . .  unfold R class.
##
##  This function unfolds  the R class  <R>.  It determines (and returns) the
##  right Schutzenberger group   of the representative  of <R>.   Moreover it
##  determines (and installs) the components 'images' and 'rMults'.
##
TransRClassOps.SchutzenbergerGroup:= function(R)

   local img, orbit, i, j, back, rep, s, pnt, new, a, n, set, sets, 
      perms, gens;

   # determine starting point.
   rep:= R.representative;  img:= Image(rep);

   # form the (weak, but graded) orbit.
   orbit:= [img];  sets:= [img]; n:= Size(img);  i:= 0;  back:= [[]];
   gens:= R.monoid.generators;

   for pnt in orbit do

      # keep track of position of 'pnt'.
      i:= i+1;

      # loop over the generators.
      for s in gens do
         new:= OnTuples(pnt, s);  set:= Set(new);

         # discard points of lower grading.
         if Size(set) = n then

            j:= Position(sets, set);
   
            # install new point, if necessary.
            if j = false then
               Add(orbit, new);  Add(sets, set);  Add(back, []);
               j:= Length(orbit);
            fi;
   
            # remember predecessor.
            AddSet(back[j], i);

         fi;
      od;
   od;

   # form the transitive closure.
   n:= Length(orbit);
   for j in [1..n] do
      for i in [1..n] do
         if j in back[i] then 
            UniteSet(back[i], back[j]);
         fi;
      od;
   od;

   # select strong orbit of point 1.
   AddSet(back[1], 1);
   orbit:= orbit{back[1]};  sets:= sets{back[1]};

   # find multipliers. 
   perms:= [];
   for pnt in orbit do
      Add(perms, MappingPermListList(pnt, img));
   od;

   # determine Schutz grp.
   new:= [];
   for i in [1..Length(sets)] do
      pnt:= sets[i];
      for s in gens do
         j:= Position(sets, OnSets(pnt, s));
         if j <> false then
            Add(new, PermLeftQuoTrans(rep, rep/perms[i] * (s*perms[j])));
         fi;
      od;
   od;

   # for the records...
   R.images:= sets;
   R.rMults:= perms;

   # return the group.
   return Group(Set(new), ());

end;

#############################################################################
##
#F  TransRClassOps.Size( <R> ) . . . . . . . . . . . . . . . . . . . .  size.
##
TransRClassOps.Size:= R -> Size(SchutzenbergerGroup(R)) * Length(R.images);

#############################################################################
##
#F  TransRClassOps.Elements( <R> ) . . . . . . . . . . . . . . . .  elements.
##
TransRClassOps.Elements:= function(R)

   local m, x, elts, grp;

   grp:= SchutzenbergerGroup(R);
   x:= Representative(R);

   elts:= [];
   for m in R.rMults do
      Append(elts, x * (grp * m^-1));
   od;

   return Set(elts);

end; 

#############################################################################
##
#F  TransRClassOps.Rank( <R> ) . . . . . . . . . . . . . . . . . . . .  rank.
##
TransRClassOps.Rank:= R -> Rank(R.representative);

#############################################################################
##
#F  <x> in <rClass>  . . . . . . . . . . . . . . . . . . . . membership test.
##
TransRClassOps.\in:= function(x, R)

   local i, rep, grp;

   # what were we  talkin about?
   if not IsTransformation(x) then
      return false;
   fi;

   # check degree, rank, and kernel.
   rep:= Representative(R); 
   if Degree(x) <> Degree(rep) or Rank(x) <> Rank(rep) 
        or Kernel(x) <> Kernel(rep) then
      return false;
   fi;

   # unfold class, if necessary.
   grp:= SchutzenbergerGroup(R);

   # check image.
   i:= Position(R.images, Image(x));
   if i = false then 
      return false;
   fi;

   # check the group.
   return PermLeftQuoTrans(rep, x * R.rMults[i]) in grp;

end;

#############################################################################
##
#F  TransMonoidOps.RClasses( <M> ) . . . . . . . . . . . . . . . . R classes.
##
TransMonoidOps.RClasses:= function(M)

   local rClasses, lTrans, rReps, images, class, n, one, gens, orb, img,
         x, s, r, k, i, j, pos, ker, kernels, new;

   # initialize.
   n:= Degree(M);  one:= M.identity;  gens:= M.generators;
   images:= List([1..n], x-> []);  class:= List([1..n], x-> []);
   rClasses:= [];  lTrans:= [];  rReps:= []; kernels:= [];
   orb:= [one];

   # loop over the orbit.
   for x in orb do

      # locate image.
      img:= Image(x);  k:= Size(img);  j:= Position(images[k], img);

      # new image means new bigger class.
      if j = false then
         r:= RClass(M, x);  Size(r);
         Add(rClasses, r);  Add(lTrans, [one]);  
         Add(rReps, [x]);  Add(kernels, [Kernel(x)]);
         Append(images[k], r.images);  j:= Length(rClasses);
         Append(class[k], List(r.images, x-> j)); 

InfoMono1("n\c");

         # install descendants in the queue.
         for s in gens do 
            Add(orb, s * x);
         od;

      else
         pos:= class[k][j];  r:= rClasses[pos];

         # adjust transformation st. img x = img r.representative.
         #C should get rid of 'Position' call here!
         x:= x * r.rMults[Position(r.images, img)];

         # locate R class of 'x'.
         ker:= Kernel(x);  new:= true;
         for i in [1..Length(kernels[pos])] do
            if new and ker = kernels[pos][i] and lTrans[pos][i] * x in r then
               new:= false;
            fi;
         od;

         if new then               

            # remember representative.
            Add(rReps[pos], x);

            # construct new l!
            Add(lTrans[pos], TransformationNC(List([1..n], 
                 i-> Position(x.images, i^r.representative))));  #  :-)

            Add(kernels[pos], ker);
InfoMono1(k, "\c");

            # install descendants in the queue.
            for s in gens do 
               Add(orb, s * x);
            od;

         fi;
      fi;
   od;

   # install things in M.  What do we need exactly??
   M.images:= images;
   M.imagePos:= class;
   M.kernels:= kernels;
   M.rClassReps:= rReps;
   M.lTrans:= lTrans;
   M.orbitClasses:= rClasses;

   # make actual list of R classes.
   orb:= [];
   for j in [1..Length(rClasses)] do
      for class in rReps[j] do
         class:= RClass(M, class);
         class.images:= rClasses[j].images;
         class.schutzenbergerGroup:= rClasses[j].schutzenbergerGroup;
         class.rMults:= rClasses[j].rMults;
         Add(orb, class);
      od;
   od;
     
   # return list of R classes.
   return orb;

end;

#############################################################################
##
#V  TransLClassOps . . . . . . . . . . . . . . . . . . . . operations record.
##
##  This is the  operations record for  L classes  of transformation monoids.
##  It inherits from the general 'LClassOps'.
##
TransLClassOps:= OperationsRecord("TransLClassOps", LClassOps);

#############################################################################
##
#F  TransMonoidOps.LClass( <M>, <x> )  . . . . . . . . . . construct L class.
##
TransMonoidOps.LClass:= function(M, x)

   local L;

   L:= MonoidOps.LClass(M, x);
   L.operations:= TransLClassOps;

   return L;

end;

#############################################################################
##
#F  TransLClassOps.SchutzenbergerGroup( <L> )  . . . . . . .  unfold L class.
##
##  this is very much the exact dual of 'TransRClassOps.SchutzenbergerGroup'.
##
TransLClassOps.SchutzenbergerGroup:= function(L)

   local ker, orbit, i, j, back, s, pnt, new, a, n, set, sets, 
      relts, gens, img, rep;

   # determine starting point.
   rep:= L.representative;  ker:= Kernel(rep);

   # form the (weak, but graded) orbit.
   orbit:= [ker];  sets:= [ker]; n:= Size(ker);  i:= 0;  back:= [[]];
   gens:= L.monoid.generators;

   for pnt in orbit do

      # keep track of position of 'pnt'.
      i:= i+1;

      # loop over the generators.
      for s in gens do
         new:= OnTuplesOfSetsAntiAction(pnt, s);  
         set:= Set(new);

         # discard points of lower grading.
         if not [] in set then

            j:= Position(sets, set);
   
            # install new point, if necessary.
            if j = false then
               Add(orbit, new);  Add(sets, set);  Add(back, []);
               j:= Length(orbit);
            fi;
   
            # remember predecessor.
            AddSet(back[j], i);

         fi;
      od;
   od;

   # form the transitive closure.
   n:= Length(orbit);
   for j in [1..n] do
      for i in [1..n] do
         if j in back[i] then 
            UniteSet(back[i], back[j]);
         fi;
      od;
   od;

   # select strong orbit of point 1.
   AddSet(back[1], 1);
   orbit:= orbit{back[1]};  sets:= sets{back[1]};

   # find multipliers.
   relts:= [];
   for pnt in orbit do
      new:= [];
      for i in [1..Length(ker)] do
         new{ker[i]}:= List(ker[i], x-> pnt[i]);
      od;
      Add(relts, Relation(new));
   od;

   # determine Schutz grp.
   new:= [];
   for i in [1..Length(sets)] do
      pnt:= sets[i];
      for s in gens do
         j:= Position(sets, OnKernelsAntiAction(pnt, s));
         if j <> false then
            Add(new, PermLeftQuoTrans(rep, 
                TransRel(relts[j] * (s * TransRel(relts[i]^-1 * rep)))));
         fi;
      od;
   od;

   # for the records...
   L.kernels:= sets;
   L.lMults:= relts;

   # return the group.
   return Group(Set(new), ());

end;


#############################################################################
##
#F  TransLClassOps.Size( <L> ) . . . . . . . . . . . . . . . . . . . .  size.
##
TransLClassOps.Size:= function(L)

   return Size(SchutzenbergerGroup(L)) * Length(L.kernels);

end;

#############################################################################
##
#F  TransLClassOps.Elements( <L> ) . . . . . . . . . . . . . . . .  elements.
##
TransLClassOps.Elements:= function(L)

   local m, x, elts, grp;

   grp:= SchutzenbergerGroup(L);
   x:= Representative(L);

   elts:= [];
   for m in L.lMults do
      Append(elts, TransRel(m^-1 * x) * grp);
   od;

   return Set(elts);

end;

#############################################################################
##
#F  TransLClassOps.Rank( <L> ) . . . . . . . . . . . . . . . . . . . .  rank.
##
TransLClassOps.Rank:= L -> Rank(L.representative);

#############################################################################
##
#F  <x> in <lClass>  . . . . . . . . . . . . . . . . . . . . membership test.
##
TransLClassOps.\in:= function(x, L)

   local i, rep, grp;

   # what were we  talkin about?
   if not IsTransformation(x) then
      return false;
   fi;

   # check degree, rank, and image.
   rep:= L.representative; 
   if Degree(x) <> Degree(rep) or Rank(x) <> Rank(rep) or 
        Image(x) <> Image(rep) then
      return false;
   fi;

   # unfold class, if necessary.
   grp:= SchutzenbergerGroup(L);

   # check kernel.
   i:= Position(L.kernels, Kernel(x));
   if i = false then 
      return false;
   fi;

   # check the group.
   return PermLeftQuoTrans(rep, TransRel(L.lMults[i] * x)) in grp;

end;

#############################################################################
##
#V  TransHClassOps . . . . . . . . . . . . . . . . . . . . operations record.
##
##  This is the  operations record for  H classes  of transformation monoids.
##  It inherits from the general 'HClassOps'.
##
TransHClassOps:= OperationsRecord("TransHClassOps", HClassOps);

#############################################################################
##
#F  TransMonoidOps.HClass( <M>, <x> )  . . . . . . . . . . construct H class.
##
TransMonoidOps.HClass:= function(M, x)

   local H;

   H:= MonoidOps.HClass(M, x);
   H.operations:= TransHClassOps;

   return H;

end;


#############################################################################
##
#F  TransHClassOps.Size( <H> ) . . . . . . . . . . . . . . . . . . . .  size.
##
TransHClassOps.Size:= function(H)

   return Size(SchutzenbergerGroup(H));

end;

#############################################################################
##
#F  TransHClassOps.Elements( <H> ) . . . . . . . . . . . . . . . .  elements.
##
TransHClassOps.Elements:= function(H)

   return Set(Representative(H) * SchutzenbergerGroup(H));

end; 

#############################################################################
##
#F  TransHClassOps.Rank( <H> ) . . . . . . . . . . . . . . . . . . . .  rank.
##
TransHClassOps.Rank:= H -> Rank(H.representative);

#############################################################################
##
#F  <x> in <hClass>  . . . . . . . . . . . . . . . . . . . . membership test.
##
TransHClassOps.\in:= function(x, H)

   local rep, img, grp;

   # what were we  talkin about?
   if not IsTransformation(x) then
      return false;
   fi;

   # check degree, rank, image and kernel.
   rep:= H.representative; 
   if Degree(x) <> Degree(rep) or Rank(x) <> Rank(rep) 
        or Image(x) <> Image(rep) or Kernel(x) <> Kernel(rep) then
      return false;
   fi;

   # check the group.
   return PermLeftQuoTrans(rep, x) in SchutzenbergerGroup(H);  

end;

#############################################################################
##
#F  TransHClassOps.SchutzenbergerGroup( <H> )  . . . .  Schutzenberger group.
##
##  This  function unfolds an H class  <H>.   Additionally the components 'L'
##  and 'R' are installed, containing  the L class and  the <R> class of  the
##  representative of <H>, respectively.
##
TransHClassOps.SchutzenbergerGroup:= function(H)

   local M, x;

   M:= H.monoid;  x:= H.representative;
   H.R:= RClass(M, x);  H.L:= LClass(M, x);
   return Intersection(SchutzenbergerGroup(H.L), SchutzenbergerGroup(H.R));

end;


#############################################################################
##
#V  TransDClassOps . . . . . . . . . . . . . . . . . . . . operations record.
##
##  This is the  operations record for  D classes  of transformation monoids.
##  It inherits from the general 'DClassOps'.
##
TransDClassOps:= OperationsRecord("TransDClassOps", DClassOps);

#############################################################################
##
#F  TransMonoidOps.DClass( <M>, <x> )  . . . . . . . . . . construct D class.
##
TransMonoidOps.DClass:= function(M, x)

   local D;

   D:= MonoidOps.DClass(M, x);
   D.operations:= TransDClassOps;

   return D;

end;

#############################################################################
##
#F  TransDClassOps.SchutzenbergerGroup( <D> )  . . . .  Schutzenberger group.
##
##  This function unfolds    a D  class.   It   determines,  for   the  given
##  representative 'x', its H class, L class and R class,  and stores them in
##  components 'H', 'L', and 'R', respectively.
##
##  Additionally, the component 'rCosets' contains a right transversal of the
##  small group in the right group.  (our favorite point of view!)
##
TransDClassOps.SchutzenbergerGroup:= function(D)

   local x, M, rGrp, grp;

   x:= D.representative;  M:= D.monoid;
   D.H:= HClass(M, x);  grp:= SchutzenbergerGroup(D.H);
   D.L:= D.H.L;  D.R:= D.H.R; 

   rGrp:= SchutzenbergerGroup(D.R);
   D.rCosets:= List(Cosets(rGrp, AsSubgroup(rGrp, grp)), Representative);

   return grp;

end;

#############################################################################
##
#F  TransDClassOps.Size( <D> ) . . . . . . . . . . . . . . . . . . . .  size.
##
TransDClassOps.Size:= function(D)

   return Size(SchutzenbergerGroup(D))^-1 * Size(D.R) * Size(D.L);

end;

#############################################################################
##
#F  TransDClassOps.Elements( <D> ) . . . . . . . . . . . . . . . .  elements.
##
TransDClassOps.Elements:= function(D)

   local c, e, m, elts;

   # unfold class.
   SchutzenbergerGroup(D);

   # list elements.
   elts:= [];
   for c in D.rCosets do 
      for m in D.R.rMults do
         for e in Elements(D.L) do
            Add(elts, e * c / m);
         od;
      od;
   od;

   return Set(elts);

end;

#############################################################################
##
#F  TransDClassOps.Rank( <D> ) . . . . . . . . . . . . . . . . . . . .  rank.
##
TransDClassOps.Rank:= D -> Rank(D.representative);

#############################################################################
##
#F  <x> in <dClass>  . . . . . . . . . . . . . . . . . . . . membership test.
##
TransDClassOps.\in:= function(x, D)

   local i, c, rep, ker, img, grp, quo;

   # what were we  talkin about?
   if not IsTransformation(x) then
      return false;
   fi;

   # check degree and rank.
   rep:= D.representative; 
   if Degree(x) <> Degree(rep) or Rank(x) <> Rank(rep) then
      return false;
   fi;

   # unfold class, if necessary.
   SchutzenbergerGroup(D); 

   # check image, and adjust.
   img:= Image(x);  i:= Position(D.R.images, img);
   if i = false then
      return false;
   fi;
   x:= x * D.R.rMults[i];
   img:= Image(x);

   # check kernel, and adjust.
   ker:= Kernel(x);  i:= Position(D.L.kernels, ker);
   if i = false then 
      return false;
   fi;
   x:= TransRel(D.L.lMults[i] * x);

   # check the (cosets of the) group.
   grp:= SchutzenbergerGroup(D.L);
   quo:= PermLeftQuoTrans(rep, x);
   for c in D.rCosets do
      if quo/c in grp then
         return true;
      fi;
   od;

   # if that fails.
   return false;

end;

#############################################################################
##
#F  TransMonoidOps.DClasses( <M> ) . . . . . . . . . . . . . . . . D classes.
##
##  relies on the component 'rClassReps' of <M> stored by 'RClasses'.
##
TransMonoidOps.DClasses:= function(M)

   local classes, reps, d;

   # start with R classes.
   RClasses(M);

   # split big classes into D classes.
   classes:= [];
   for reps in M.rClassReps do
      repeat
         d:= DClass(M, reps[1]);
         Add(classes, d);
         reps:= Filtered(reps, x-> not x in d);
      until reps = [];
   od;

   return classes;

end;

#############################################################################
##
#F  IsRegularDClass( <D> ) . . . . . . . . . . . . . . . .  check regularity.
##
##  There must be a better way to find this if the D class is unfolded!
##
IsRegularDClass:= function(D)

   if not IsBound(D.isRegular) then   
      D.isRegular:= IsRegularTrans(D.monoid, D.representative);
   fi;

   return D.isRegular;

end;

#############################################################################
##
#F  TransMonoidOps.Size( <M> ) . . . . . . . . . . . . . . . . . . . .  size.
##
TransMonoidOps.Size:= function(M)

   # trigger structure calculations.
   RClasses(M);

   # calculate size.
   return List(M.orbitClasses, Size) * List(M.rClassReps, Length);

end;

#############################################################################
##
#F  TransMonoidOps.Display( <M> )  . . . . . . . . . . . . . . . . . display.
##
TransDClassOps.Display:= function(D, opt)

   local sh;

   if IsRegularDClass(D) then         
      Print("[", Size(D.H), ".", Length(D.R.rMults), ".",
            Length(D.L.lMults), "] \c");
   else
      sh:= Size(D.H);
      Print("{", sh,
            ".", Length(D.R.rMults), "x", 
               Size(SchutzenbergerGroup(D.R))/sh,   
            ".", Length(D.L.lMults), "x", 
            Size(SchutzenbergerGroup(D.L))/sh,
            "} \c");
   fi;
end;

TransMonoidOps.Display:= function(M, opt)

   local dc, i, len, sh, D, layer;

   # determine D classes and sort according to rank.
   layer:= List([1..Degree(M)], x-> []);
   for D in DClasses(M) do
      Add(layer[Rank(D)], D);
   od;

   # loop over the layers.
   len:= Length(layer);
   for i in [len, len-1 .. 1] do
      if layer[i] <> [] then
         Print("Rank ", i, ": \c");

         # loop over D classes.
         for D in layer[i] do
            Display(D);
         od;
         Print("\n");
      fi;
   od;

end;

#############################################################################
##
#F  TransDClassOps.LClasses( <D> ) . . . . . . . . . . . . . . . . L classes.
##
TransDClassOps.LClasses:= function(D)

   local M, x, c, d, m, classes, new, gens;

   # initialize.
   classes:= [];  M:= D.monoid;  x:= D.representative;  Size(D);
   gens:= Generators(SchutzenbergerGroup(D.L));

   # loop over the cosets.
   for c in D.rCosets do

      # loop over the L class representatives.
      for m in D.R.rMults do
         d:= c/m;
         new:= LClass(M, x * d);
         new.kernels:= D.L.kernels;
         new.lMults:= D.L.lMults;
         new.schutzenbergerGroup:= Group(List(gens, x-> x^d), ());
         Add(classes, new);
      od;
   od;

   # return the list of L classes.
   return classes;

end;

#############################################################################
##
#F  TransDClassOps.RClasses( <D> ) . . . . . . . . . . . . . . . . R classes.
##
TransDClassOps.RClasses:= function(D)

   local M, x, l, c, d, classes, new, grp, sch, cos;

   # initialize. 
   classes:= [];  M:= D.monoid;  x:= D.representative;  Size(D);
   grp:= SchutzenbergerGroup(D.L);
   sch:= AsSubgroup(grp, SchutzenbergerGroup(D.H));
   cos:= List(LeftCosets(grp, sch), Representative);
   grp:= SchutzenbergerGroup(D.R);

   # loop over R class reps.
   for l in D.L.lMults do
      d:= TransRel(l^-1 * x);

      # loop over cosets.
      for c in cos do
         new:= RClass(M, d * c);
         new.images:= D.R.images;
         new.rMults:= D.R.rMults;
         new.schutzenbergerGroup:= grp;
         Add(classes, new);
      od;
   od;

   # return list of R classes.
   return classes;

end;

#############################################################################
##
#F  TransMonoidOps.LClasses( <M> ) . . . . . . . . . . . . . . . . L classes.
##
TransMonoidOps.LClasses:= function(M)

   local classes, D;

   classes:= [];
   for D in DClasses(M) do
      Append(classes, LClasses(D));
   od;

   return classes;

end;

#############################################################################
##
#F  TransMonoidOps.HClasses( <M> ) . . . . . . . . . . . . . . . . H classes.
##
TransMonoidOps.HClasses:= M -> Concatenation(List(DClasses(M), HClasses));

#############################################################################
##
#F  TransDClassOps.HClasses( <D> ) . . . . . . . . . . . . . . . . H classes
##
TransDClassOps.HClasses:= D -> Concatenation(List(LClasses(D), HClasses));

#############################################################################
##
#F  TransRClassOps.HClasses( <R> ) . . . . . . . . . . . . . . . . H classes.
##
TransRClassOps.HClasses:= function(R)

   local M, D, x, c, d, m, classes, new;

   # initialize.
   classes:= [];  M:= R.monoid;  x:= R.representative;
   D:= DClass(M, x);  Size(D);

   # loop over the cosets.
   for c in D.rCosets do

      # loop over the class representatives.
      for m in R.rMults do
         d:= c/m;
         new:= HClass(M, x * d);
         Add(classes, new);
      od;
   od;

   # return the list of H classes.
   return classes;

end;

#############################################################################
##
#F  TransLClassOps.HClasses( <L> ) . . . . . . . . . . . . . . . . H classes
##
##  determines the list of H classes inside a  given L class.
##
TransLClassOps.HClasses:= function(L)

   local M, D, x, l, c, d, classes, new, grp, sch, cos;

   # initialize. 
   classes:= [];  M:= L.monoid;  x:= L.representative;
   D:= DClass(M, x);  Size(D);

   # determine groups and coset reps.
   grp:= SchutzenbergerGroup(D.L);
   sch:= AsSubgroup(grp, SchutzenbergerGroup(D.H));
   cos:= List(LeftCosets(grp, sch), Representative);

   # loop over R class reps.
   for l in D.L.lMults do
      d:= TransRel(l^-1 * x);

      # loop over cosets.
      for c in cos do
         new:= HClass(M, d * c);
         new.schutzenbergerGroup:= sch;  # remins the same throughout <L>.
         Add(classes, new);
      od;
   od;

   # return list of H classes.
   return classes;

end;

#############################################################################
##
#F  KernelsTransMonoid( <M> )  . . . . . . . . . . . . . . . . . . . kernels.
##
KernelsTransMonoid:= function(M)

   return Union(
     GradedOrbit(M, Kernel(M.identity), Size, OnKernelsAntiAction));

end; 

#############################################################################
##
#F  ImagesTransMonoid( <M> ) . . . . . . . . . . . . . . . . . . . .  images.
##
ImagesTransMonoid:= function(M)

   return Union(GradedOrbit(M, Image(M.identity), Size, OnSets));

end;

#############################################################################
##
#F  TransMonoidOps.Idempotents( <M> )  . . . . . . . . . . . . . idempotents.
##
TransMonoidOps.Idempotents:= function(M)

   local isCrossSection, idempotent, pt, ker, img, kers, imgs, i, n, idm;

   # how to determine whether 'img' is a cross section of 'ker'.
   isCrossSection:= function(ker, img)
      return ForAll(ker, k-> Number(k, i-> i in img) = 1);
   end;  # there might be a more efficient way...

   # how to construct an idempotent with kernel 'ker' and image 'img'.
   idempotent:= function(ker, img)
      local e, l;
      e:= [];
      for l in ker do  
         e{l}:= 0*l + Intersection(l, img)[1];  # :-)
      od;
      return TransformationNC(e);
   end;

   n:= Degree(M);  idm:= [];

   # determine kernels with grading.
   kers:= List([1..n], x-> []);
   for pt in KernelsTransMonoid(M) do
      Add(kers[Size(pt)], pt);
   od;

   # determine images with grading.
   imgs:= List([1..n], x-> []);
   for pt in ImagesTransMonoid(M) do
      Add(imgs[Size(pt)], pt);
   od;

   # loop over all ranks.
   for i in [1..n] do

      # loop over the kernels.
      for ker in kers[i] do

         # loop over the images.
         for img in imgs[i] do

            # check for cross section.
            if isCrossSection(ker, img) then
               Add(idm, idempotent(ker, img));
            fi;
         od;
      od;
   od;

   # return the set of idempotents.
   return Set(idm);

end;

#############################################################################
##
#F  TransDClassOps.Idempotents( <D> ) . . . . . . . . . . . . .  idempotents.
##
TransDClassOps.Idempotents:= function(D)

   local isCrossSection, idempotent, ker, img, idm;

   # how to determine whether 'img' is a cross section of 'ker'.
   isCrossSection:= function(ker, img)
      return ForAll(ker, k-> Number(k, i-> i in img) = 1);
   end;  # there might be a more efficient way...

   # how to construct an idempotent with kernel 'ker' and image 'img'.
   idempotent:= function(ker, img)
      local e, l;
      e:= [];
      for l in ker do  
         e{l}:= 0*l + Intersection(l, img)[1];  # :-)
      od;
      return TransformationNC(e);
   end;

   idm:= [];  Size(D);  # unpack D class.

   # loop over the kernels.
   for ker in D.L.kernels do

      # loop over the images.
      for img in D.R.images do

         # check for cross section.
         if isCrossSection(ker, img) then
            Add(idm, idempotent(ker, img));
         fi;
      od;
   od;

   # return the list of idempotents.
   return idm;

end;

#############################################################################
##
#F  TransRClassOps.Idempotents( <R> ) . . . . . . . . . . . . .  idempotents.
##
TransRClassOps.Idempotents:= function(R)

   local isCrossSection, idempotent, ker, img, idm;

   # how to determine whether 'img' is a cross section of 'ker'.
   isCrossSection:= function(ker, img)
      return ForAll(ker, k-> Number(k, i-> i in img) = 1);
   end;  # there might be a more efficient way...

   # how to construct an idempotent with kernel 'ker' and image 'img'.
   idempotent:= function(ker, img)
      local e, l;
      e:= [];
      for l in ker do  
         e{l}:= 0*l + Intersection(l, img)[1];  # :-)
      od;
      return TransformationNC(e);
   end;

   idm:= [];  Size(R);  # unpack R class.

   ker:= Kernel(R.representative);

   # loop over the images.
   for img in R.images do

      # check for cross section.
      if isCrossSection(ker, img) then
         Add(idm, idempotent(ker, img));
      fi;
   od;

   # return the list of idempotents.
   return idm;

end;

#############################################################################
##
#F  TransLClassOps.Idempotents( <L> ) . . . . . . . . . . . . .  idempotents.
##
TransLClassOps.Idempotents:= function(L)

   local isCrossSection, idempotent, ker, img, idm;

   # how to determine whether 'img' is a cross section of 'ker'.
   isCrossSection:= function(ker, img)
      return ForAll(ker, k-> Number(k, i-> i in img) = 1);
   end;  # there might be a more efficient way...

   # how to construct an idempotent with kernel 'ker' and image 'img'.
   idempotent:= function(ker, img)
      local e, l;
      e:= [];
      for l in ker do  
         e{l}:= 0*l + Intersection(l, img)[1];  # :-)
      od;
      return TransformationNC(e);
   end;

   idm:= [];  Size(L);  # unpack L class.

   img:= Image(L.representative);

   # loop over the kernels.
   for ker in L.kernels do

      # check for cross section.
      if isCrossSection(ker, img) then
         Add(idm, idempotent(ker, img));
      fi;
   od;

   # return the list of idempotents.
   return idm;

end;

#############################################################################
##
#F  TransHClassOps.Idempotents( <H> ) . . . . . . . . . . . . .  idempotents.
##
TransHClassOps.Idempotents:= function(H)

   local isCrossSection, idempotent, ker, img, idm;

   # how to determine whether 'img' is a cross section of 'ker'.
   isCrossSection:= function(ker, img)
      return ForAll(ker, k-> Number(k, i-> i in img) = 1);
   end;  # there might be a more efficient way...

   # how to construct an idempotent with kernel 'ker' and image 'img'.
   idempotent:= function(ker, img)
      local e, l;
      e:= [];
      for l in ker do  
         e{l}:= 0*l + Intersection(l, img)[1];  # :-)
      od;
      return TransformationNC(e);
   end;

   idm:= [];  Size(H);  # unpack L class.

   img:= Image(H.representative);
   ker:= Kernel(H.representative);

   # check for cross section.
   if isCrossSection(ker, img) then
      Add(idm, idempotent(ker, img));
   fi;

   # return the list of idempotents.
   return idm;

end;

#############################################################################
##
#V  FullTransMonoidOps . . . . . . . . . . . . . . . . . . operations record.
##
FullTransMonoidOps:= OperationsRecord("FullTransMonoidOps", TransMonoidOps);

#############################################################################
##
#F  FullTransMonoid( <n> ) . . . . . . . . . . .  full transformation monoid.
##
##  The full transformation monoid $T_n$ on $n$  points is generated by three
##  transformations
##  
##     a:= Transformation([2, 1, 3, 4, ..., n]),
##     b:= Transformation([n, 1, 2, 3, ..., n-1]), and
##     c:= Transformation([2, 2, 3, 4, ..., n]).
##  
##  Note that 'a' and 'b' generate the symmetric  group on $n$ points and 'c'
##  is  the   transformation  that  identifies  points 1   and 2   and leaves
##  everything else fixed.
##
FullTransMonoid:= function(n)

   local T, a, b, c;

   # check argument.
   if n < 1 then
      Error("<n> must be a positive integer");
   fi;

   # trivial cases.
   if n = 1 then 
      T:= Monoid(TransformationNC([1]));

   elif n = 2 then
      T:= Monoid(List([[2, 1], [2, 2]], TransformationNC));

   else

      a:= [1..n];  a[1]:= 2;  a[2]:= 1;
      b:= [0..n-1]; b[1]:= n;
      c:= [1..n];  c[1]:= 2;

      T:= Monoid(List([a, b, c], TransformationNC));
   fi;

   T.operations:= FullTransMonoidOps;

   return T;

end;

############################################################################# 
##
#F  FullTransMonoidOps.Size  . . . . . . . . . . . . . . . . . . . . .  size.
##
##  The size of the full transformation monoid $T_n$ on  $n$ points is $n^n$,
##  there are $n$ choices for the images of $n$ places.
##
FullTransMonoidOps.Size:=  M -> Degree(M)^Degree(M);

#############################################################################
##
#F  FullTransMonoidOps.Random( <D> ) . . . . . . . . . . . .  random element.
##
FullTransMonoidOps.Random:= function(D)

   local N;

   N:= [1..Degree(D)];

   return TransformationNC(List(N, x-> Random(N)));

end;

#############################################################################
##
#V  PartialTransMonoidOps  . . . . . . . . . . . . . . . . operations record.
##
PartialTransMonoidOps:= OperationsRecord("PartialTransMonoidOps", 
   TransMonoidOps);

#############################################################################
##
#F  PartialTransMonoid( <n> )  . . . . . . .  partial transformations monoid.
##
##  The undefined point is 'n+1'.
##
PartialTransMonoid:= function(n)

   local T, a, b, c, d;

   # define the generators.
   a:= [1..n+1];  a[1]:= 2;  a[2]:= 1;
   b:= [0..n];  b[1]:= n;  b[n+1]:= n+1;
   c:= [1..n+1];  c[1]:= n+1;
   d:= [1..n+1];  d[1]:= 2;

   # define the monoid and install operations.
   T:= Monoid(List([a, b, c, d], Transformation));
   T.operations:= PartialTransMonoidOps;

   # return the monoid.
   return T;

end;

#############################################################################
##
#F  PartialTransMonoidOps.Size( <M> )  . . . . . . . . . . . . . . . .  size.
##
PartialTransMonoidOps.Size:= function(M)

   local n;

   n:= Degree(M);

   return (n+1)^n;

end;

############################################################################# 
##
#F  PartialTransMonoidOps.Random( <D> )  . . . . . . . . . .  random element.
##
PartialTransMonoidOps.Random:= function(D)

   local N, lst;

   N:= [1..Degree(D)+1];
   lst:= List(N, x-> Random(N));  lst[Length(lst)]:= Length(lst);

   return TransformationNC(lst);

end;

#############################################################################
##
#E  Emacs  . . . . . . . . . . . . . . . . . . . . . . local emacs variables.
##
##  Local Variables:
##  mode:               outline
##  outline-regexp:     "#F\\|#V\\|#A\\|#E"
##  fill-column:        77
##  fill-prefix:        "##  "
##  eval:               (hide-body)
##  End:
##
