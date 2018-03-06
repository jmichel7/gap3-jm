#############################################################################
##
#A  GAP                                                 Goetz.Pfeiffer@UCG.IE
##
#A  $Id: action.g,v 2.4 1997/11/24 17:37:01 goetz Exp $
##
#Y  Copyright (C) 1997, Mathematics Dept, University College Galway, Ireland.
##
##  This file contains the functions for actions of monoids.
##

#############################################################################
##
#F  MonoidOps.Orbit( <M>, <d>, <op> )  . . . . . . . . . . . . . . . . orbit.
##
##  The *orbit* of point <d> under the action of the monoid <M> via operation
##  <op> is the set of all points that can be reached.
##
MonoidOps.Orbit:= function(M, d, op)

   local pnt, orbit, gens, s, new;

   # initialize.
   orbit:= [d];  gens:= M.generators;

   # loop over the orbit.
   for pnt in orbit do
      for s in gens do
         new:= op(pnt, s);
         if not new in orbit then
            Add(orbit, new);
         fi;
      od;
   od;

   # return the orbit.
   return orbit;

end;


#############################################################################
##
#F  MonoidOps.Orbits( <M>, <D>, <op> ) . . . . . . . . . . . . . . .  orbits.
##
MonoidOps.Orbits:= function(M, D, op)

   local pnt, orbits, orb, gens, s, new, pos;

   # initialize.
   D:= Set(D);  orbits:= [];  gens:= M.generators;

   # keep constructing orbits.
   while D <> [] do
      orb:= [D[1]];
      for pnt in orb do
         for s in gens do
            new:= op(pnt, s);
            if not new in orb then
               Add(orb, new);
            fi;
         od;
      od;

      Sort(orb);

      pos:= PositionProperty(orbits, x-> IsSubset(orb, x));
      if pos = false then
         Add(orbits, orb);
         pos:= Length(orbits);
      else
         orbits[pos]:= orb;
      fi;
      SubtractSet(D, orb);

   od;

   # return the orbits.
   return orbits;

end;


#############################################################################
##
#F  ShortOrbit( <M>, <d>, <op>, <grad> ) . . . . . . . . . . . . short orbit.
##
##  Let <M> be  a monoid acting  on the  set <D>.  A  *grading*  is a map  $g
##  \colon D \to \Z$ such that $g(d) \geq (d^m)$ for all $d \in D$ and all $m
##  \in M$.
##  
##  The *short orbit* of  the point <d> in  $D$ under the  action of <M> with
##  respect to   the grading  <grad>  is  the  set  $\{d\^m   \mid  m  \in M,
##  <grad>(d\^m) = <grad>(d)\}$.
##
ShortOrbit:= function(M, d, op, grad)
   return M.operations.ShortOrbit(M, d, op, grad);
end;

#############################################################################
##
#F  MonoidOps.ShortOrbit( <M>, <d>, <grad>, <op> ) . . . . . . . short orbit.
##
MonoidOps.ShortOrbit:= function(M, d, op, grad)

   local orb, n, pnt, new, s;

   # initialize.
   orb:= [d];  n:= grad(d);

   # from the orbit.
   for pnt in orb do
      for s in M.generators do
         new:= op(pnt, s);
         if grad(new) = n and not new in orb then
            Add(orb, new);
         fi;
      od;
   od;

   # return the orbit.
   return orb;

end;

#############################################################################
##
#F  GradedOrbit( <M>, <d>, <op>, <grad> )  . . . . . . . . . .  graded orbit.
##
##  The *graded orbit* of the point  <d> in $D$ under  the action of <M> with
##  respect to the grading <grad> is the list '[$O_1$, $O_2$, ...  ]' of sets
##  $O_i = \{d\^m \mid m  \in M, <grad>(d\^m) =  i\}$.  Thus the orbit of <d>
##  is simply the union of the sets $O_i$.
##
GradedOrbit:= function(M, d, op, grad)
   return M.operations.GradedOrbit(M, d, op, grad); 
end;

#############################################################################
##
#F  MonoidOps.GradedOrbit( <M>, <d>, <op>, <grad> )  . . . . .  graded orbit.
##
MonoidOps.GradedOrbit:= function(M, d, op, grad)

   local i, orb, n, pnt, new, s, g;

   # initialize.
   n:= grad(d); 
   orb:= List([1..n], x-> []);  Add(orb[n], d);

   # form the orbit.
   for i in [n, n-1..1] do
      for pnt in orb[i] do
         for s in M.generators do
            new:= op(pnt, s);
            g:= grad(new);
            if not new in orb[g] then
               Add(orb[g], new);
            fi;
         od;
      od;
   od;

   # return the orbit.
   return orb;

end;

#############################################################################
##
#F  StrongOrbit( <M>, <d>, <op>[, <grad> ] ) . . . . . . . . .  strong orbit.
##
##  The *strong orbit* of a point <d> under the action of <M> consists of all
##  those points of the orbit of <d> from where  <d> can be reached, i.e. all
##  predecessors of <d>.  One can take advantage of a grading.
##
StrongOrbit:= function(arg)

   # default grading is 'x -> 1'.
   if Length(arg) = 3 then
      return arg[1].operations.StrongOrbit(arg[1], arg[2], arg[3], x-> 1);

   elif Length(arg) = 4 then
      return arg[1].operations.StrongOrbit(arg[1], arg[2], arg[3], arg[4]);

   else
      Error("usage: StrongOrbit( <M>, <d>, <op>[, <grad> ] )");

   fi;

end;

#############################################################################
##
#F  MonoidOps.StrongOrbit( <M>, <d>, <opr>, <grad> ) . . . . .  strong orbit.
##
#C  This algorithm computes connected components and can be improved!
##
MonoidOps.StrongOrbit:= function(M, d, opr, grad)

   local orbit, i, j, back, s, pnt, new, a, n;

   # initialize.
   orbit:= [d];  i:= 0;  back:= [[]];  n:= grad(d);

   # form the (weak, but graded) orbit.
   for pnt in orbit do

      # keep track of position of 'pnt'.
      i:= i+1;

      # loop over the generators.
      for s in M.generators do
         new:= opr(pnt, s);

         # discard point of lower grading.
         if grad(new) = n then

            j:= Position(orbit, new);
   
            # install new point, if necessary.
            if j = false then
               Add(orbit, new);  Add(back, []);
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

   # return predecessors of point 1.
   AddSet(back[1], 1);
   return orbit{back[1]};

end;
      
#############################################################################
##
#F  StrongOrbits( <M>, <D>, <opr>[, <grad> ] ) . . . . . . . . strong orbits.
##
StrongOrbits:= function(arg)

   # default grading is 'x-> 1'.
   if Length(arg) = 3 then
      return arg[1].operations.StrongOrbits(arg[1], arg[2], arg[3], x-> 1);

   elif Length(arg) = 4 then
      return arg[1].operations.StrongOrbits(arg[1], arg[2], arg[3], arg[4]);

   else
      Error("usage: StrongOrbits( <M>, <D>, <opr>[, <grad> ] )");

   fi;

end;

#############################################################################
##
#F  MonoidOps.StrongOrbits( <M>, <D>, <opr>, <grad> )  . . . . strong orbits.
##
MonoidOps.StrongOrbits:= function(M, D, opr, grad)

   local orbit, orbits, i, j, back, s, pnt, new, a, n;

   D:= Set(D);  orbits:= [];

   while D <> [] do

      # initialize.
      orbit:= [D[1]];  i:= 0;  back:= [[]];  n:= grad(orbit[1]);
   
      # form the (weak, but graded) orbit.
      for pnt in orbit do
   
         # keep track of position of 'pnt'.
         i:= i+1;
   
         # loop over the generators.
         for s in M.generators do
            new:= opr(pnt, s);
   
            # discard point of lower grading.
            if grad(new) = n then
   
               j:= Position(orbit, new);
      
               # install new point, if necessary.
               if j = false then
                  Add(orbit, new);  Add(back, []);
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
   
      # select predecessors of point 1.
      AddSet(back[1], 1);
      orbit:= orbit{back[1]};
   
      Add(orbits, orbit);
      SubtractSet(D, orbit);

   od;

   # return list of orbits.
   return orbits;

end;
      
#############################################################################
##
#F  Action( <M>, <D> ) . . . . . . . . . . .  action of a monoid on a domain.
##
##  'Action'  returns a  transformation  monoid  like 'Operation' returns   a
##  permutation group.
##
Action:= function(arg)

   local mon;

   # default operation is 'OnPoints'.
   if Length(arg) = 2  then
      mon:= arg[1].operations.Action(arg[1], arg[2], OnPoints);
      if not IsBound(mon.action) then
         mon.action:= rec();
      fi;
      mon.action.structure := arg[1];
      mon.action.domain    := arg[2];
      mon.action.operation := OnPoints;
   elif Length(arg) = 3  then
      mon:= arg[1].operations.Action(arg[1], arg[2], arg[3]);
      if not IsBound(mon.action) then
         mon.action:= rec();
      fi;
      mon.action.structure := arg[1];
      mon.action.domain    := arg[2];
      mon.action.operation := arg[3];
   else
      Error("usage: Action( <M>, <D> [, <action> ] ) ");
   fi;

   # return the action.
   return mon;

end;

#############################################################################
##
#F  MonoidOps.Action( <M>, <D>, <opr> )  . . . . . . . . . .  generic action.
##
MonoidOps.Action:= function (M, D, opr)

   local mon, new, trn, x, i;

   # make the transformations.
   new:= [];
   for x in M.generators  do
      trn:= [];
      for i in [1..Length(D)] do
         trn[i]:=  Position(D, opr(D[i], x));
      od;
      Add(new, Transformation(trn));    # check arguments!
   od;
   mon:= Monoid(new);
   mon.action:= rec(images:= new);

   # return the transformation monoid.
   return mon;

end;

#############################################################################
##
#F  ActionWithZero( <M>, <D> ) . . . . . . .  action of a monoid on a domain.
##
##  'ActionWithZero' returns a transformation monoid like 'Operation' returns
##  a permutation group.
##
ActionWithZero:= function(arg)

   local mon;

   # default operation is 'OnPoints'.
   if Length(arg) = 2  then
      mon:= arg[1].operations.ActionWithZero(arg[1], arg[2], OnPoints);
      if not IsBound(mon.action) then
         mon.action:= rec();
      fi;
      mon.action.structure := arg[1];
      mon.action.domain    := arg[2];
      mon.action.operation := OnPoints;
   elif Length(arg) = 3  then
      mon:= arg[1].operations.ActionWithZero(arg[1], arg[2], arg[3]);
      if not IsBound(mon.action) then
         mon.action:= rec();
      fi;
      mon.action.structure := arg[1];
      mon.action.domain    := arg[2];
      mon.action.operation := arg[3];
   else
      Error("usage: ActionWithZero( <M>, <D> [, <action> ] )");
   fi;

   # return the action.
   return mon;

end;

#############################################################################
##
#F  MonoidOps.ActionWithZero( <M>, <D>, <opr> )  . . . . . .  generic action.
##
MonoidOps.ActionWithZero:= function (M, D, opr)

   local mon, new, trn, pos, x, i, zero;

   # make the transformations.
   new:= [];  zero:= Length(D) + 1;
   for x in M.generators  do
      trn:= [];
      for i in [1..Length(D)] do
         pos:=  Position(D, opr(D[i], x));
         if pos = false then
            trn[i]:= zero;
         else
            trn[i]:= pos;
         fi;
      od;
      trn[zero]:= zero;
      Add(new, TransformationNC(trn));
   od;
   mon:= Monoid(new);
   mon.action:= rec(images:= new);

   # return the transformation monoid.
   return mon;

end;

#############################################################################
##
#F  OnLeftAntiAction( <l>, <r> ) . . . . . . . . . . . . . . . . left action.
##
OnLeftAntiAction:= OnLeftAntiOperation;

#############################################################################
##
#E  Emacs  . . . . . . . . . . . . . . . . . . . . . . local emacs variables.
##
##  Local Variables:
##  mode:               outline
##  outline-regexp:     "#F\\|#V\\|#E\\|#A"
##  fill-column:        77
##  fill-prefix:        "##  "
##  eval:               (hide-body)
##  End:
##
