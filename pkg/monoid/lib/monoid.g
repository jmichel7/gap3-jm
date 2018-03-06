#############################################################################
##
#A  GAP                                                 Goetz.Pfeiffer@UCG.IE
##
#A  $Id: monoid.g,v 2.1 1997/10/29 17:52:25 goetz Exp $
##
#Y  Copyright (C) 1997, Mathematics Dept, University College Galway, Ireland.
##
##  This file contains the definition of the domains of semigroups + monoids.
##

#############################################################################
##
#F  InfoMono? . . . . . . . . . . . . . . . . . . . . . . . . info functions.
##
if not IsBound(InfoMono1) then InfoMono1:= Ignore; fi;
if not IsBound(InfoMono2) then InfoMono2:= Ignore; fi;

#############################################################################
##
#V  MonoidElementsOps  . . . . . . . . . . . . . . . . . . operations record.
##
MonoidElementsOps:= OperationsRecord("MonoidElementsOps", DomainOps);  

#############################################################################
##
#V  MonoidElements . . . . . . . . . . . . . . . . . . . . . . . . .  domain.
##
MonoidElements:= rec();
MonoidElements.name:= "MonoidElements";
MonoidElements.isDomain:= true;
MonoidElements.isFinite:= false;
MonoidElements.size:= "infinity";
MonoidElements.operations:= MonoidElementsOps;

#############################################################################
##
#F  IsMonoidElement( <obj> ) . . . . . . . . . . . . . . . . . .  type check.
##
IsMonoidElement:= function(obj)

   return IsRec(obj) and IsBound(obj.isMonoidElement)
                     and obj.isMonoidElement;

end;

#############################################################################
##
#F  <x> in MonoidElements  . . . . . . . . . . . . . . . . . membership test.
##
MonoidElementsOps.\in:= function(x, MonoidElements)

   return IsMonoidElement(x);

end;

#############################################################################
##
#V  SemiGroupOps . . . . . . . . . . . . . . . . . . . . . operations record.
##
SemiGroupOps:= OperationsRecord("SemiGroupOps", DomainOps);  

#############################################################################
##
#F  SemiGroup( MonoidElements, <gens> )  . . . . . . . construct a semigroup.
##
MonoidElementsOps.SemiGroup:= function(MonoidElements, gens)

   local S;

   S:= rec();
   S.isDomain:= true;
   S.isSemiGroup:= true;
   S.generators:= ShallowCopy(gens);
   S.operations:= SemiGroupOps;

   return S;

end;

#############################################################################
##
#V  MonoidOps  . . . . . . . . . . . . . . . . . . . . . . operations record.
##
MonoidOps:= OperationsRecord("MonoidOps", SemiGroupOps);  

#############################################################################
##
#F  Monoid( MonoidElements, <gens>, <id> ) . . . . . . .  construct a monoid.
##
MonoidElementsOps.Monoid:= function(MonoidElements, gens, id)

   local M;

   M:= rec();
   M.isDomain:= true;
   M.isMonoid:= true;
   M.identity:= id;
   M.generators:= ShallowCopy(gens);
   M.operations:= MonoidOps;

   return M;

end;

#############################################################################
##
#F  SemiGroup( <gen1>, ... ) . . . . . . . . . . . . . construct a semigroup.
##
SemiGroup:= function(arg)

   local  S, D, gens, i;

   if Length(arg) = 1 and IsDomain(arg[1])  then
      S:= arg[1].operations.SemiGroup(arg[1]);
   elif Length(arg) = 1 and IsList(arg[1]) then
      gens:= arg[1];
      D:= Domain(gens);
      S:= D.operations.SemiGroup(D, gens);
   elif Length(arg) = 2 and IsMat(arg[1])  then
      gens:= arg;
      D:= Domain(gens);
      S:= D.operations.SemiGroup(D, gens);
   elif Length(arg) = 2 and IsList(arg[1])  then
      gens:= arg[1];
      D:= Domain(gens);
      S:= D.operations.SemiGroup(D, gens);
   elif 0 < Length(arg)  then
      gens:= arg;
      D:= Domain(gens);
      S:= D.operations.SemiGroup(D, gens);
   else
      Error("usage: SemiGroup( <gen1>, ... ) or SemiGroup( <D> )");
   fi;

   for i  in [1..Length(S.generators)]  do
      S.(i):= S.generators[i];
   od;

   return S;

end;

#############################################################################
##
#F  SemiGroupOps.Elements( <S> ) . . . . . . . . . . . . . . . . .  elements.
##
SemiGroupOps.Elements:= function(S)

   local gens, elts, n, w, s;

   gens:= Generators(S);

   #C  Is this the place to make sure there are no doubles in gens?
   elts:= Set(gens);
   for w in elts do
      for s in gens do
         n:= w * s;
         if not n in elts then
            Add(elts, n);
         fi;
      od;
   od;

   return Set(elts);

end;

#############################################################################
##
#F  SemiGroupOps.Print( <S> )  . . . . . . . . . . . . . . . . . . . . print.
##
SemiGroupOps.Print:= function(S)

   if IsBound(S.name) then
      Print(S.name);
   else
      Print("SemiGroup( ", Generators(S), " )");
   fi;

end;

#############################################################################
##
#F  IsSemiGroup( <obj> ) . . . . . . . . . . . . . . . . . . . .  type check.
##
IsSemiGroup:= function(obj)

   return IsRec(obj) and IsBound(obj.isSemiGroup) and obj.isSemiGroup;

end;

#############################################################################
##
#F  Monoid( <gen1>, ... )  . . . . . . . . . . . . . . .  construct a monoid.
##
Monoid:= function(arg)

   local  M, D, gens, id, i;

   if Length(arg) = 1 and IsDomain(arg[1])  then
      M:= arg[1].operations.Monoid(arg[1]);
   elif Length(arg) = 1 and IsList(arg[1]) then
      gens:= arg[1];
      id:= gens[1]^0;
      D:= Domain(gens);
      M:= D.operations.Monoid(D, gens, id);
   elif Length(arg) = 2 and IsMat(arg[1])  then
      gens:= arg;
      id:= arg[1]^0;
      D:= Domain(gens);
      M:= D.operations.Monoid(D, gens, id);
   elif Length(arg) = 2 and IsList(arg[1])  then
      gens:= arg[1];
      id:= arg[2];
      D:= Domain(Concatenation(gens, [ id ]));
      M:= D.operations.Monoid(D, gens, id);
   elif 0 < Length(arg)  then
      gens:= arg;
      id:= arg[1]^0;
      D:= Domain(gens);
      M:= D.operations.Monoid(D, gens, id);
   else
      Error("usage: Monoid(<gen1>,...) or Monoid(<gens>,<id>) or Monoid(<D>)");
   fi;
   for i  in [1..Length(M.generators)]  do
      M.(i):= M.generators[i];
   od;

   return M;

end;

#############################################################################
##
#F  MonoidOps.Elements( <M> )  . . . . . . . . . . . . . . . . . .  elements.
##
MonoidOps.Elements:= function(M)

   local gens, elts, n, w, s;

   gens:= Generators(M);

   #C  Is this the place to make sure there are no doubles in gens?
   elts:= [M.identity];
   for w in elts do
      for s in gens do
         n:= w * s;
         if not n in elts then
            Add(elts, n);
         fi;
      od;
   od;

   return Set(elts);

end;

#############################################################################
##
#F  MonoidOps.Print( <M> ) . . . . . . . . . . . . . . . . . . . . . . print.
##
MonoidOps.Print:= function(M)

   if IsBound(M.name) then
      Print(M.name);
   else
      Print("Monoid( ", Generators(M), " )");
   fi;

end;

#############################################################################
##
#F  IsMonoid( <obj> )  . . . . . . . . . . . . . . . . . . . . .  type check.
##
IsMonoid:= function(obj)

   return IsRec(obj) and IsBound(obj.isMonoid) and obj.isMonoid;

end;

#############################################################################
##
#V  RClassOps  . . . . . . . . . . . . . . . . . . . . . . operations record.
##
RClassOps:= OperationsRecord("RClassOps", DomainOps);

#############################################################################
##
#F  RClass( <M>, <x> ) . . . . . . . . . . . . . . . . . .  define a R class.
##
RClass:= function(M, x)

   return M.operations.RClass(M, x);

end;

MonoidOps.RClass:= function(M, x)

   local rClass;

   rClass:= rec();
   rClass.isRClass:= true;
   rClass.isDomain:= true;
   rClass.monoid:= M;
   rClass.representative:= x;
   rClass.operations:= RClassOps;

   return rClass;

end;

#############################################################################
##
#F  IsRClass( <obj> )  . . . . . . . . . . . . . . . . . . . . .  type check.
##
IsRClass:= function(obj)

   return IsRec(obj) and IsBound(obj.isRClass) and obj.isRClass;

end;

#############################################################################
##
#F  RClassOps.Print( <rClass> )  . . . . . . . . . . . . . . . . . . . print.
##
RClassOps.Print:= function(rClass)

   Print("RClass( ", rClass.monoid, ", ", rClass.representative, " )");

end;
   
#############################################################################
##
#F  <l> = <r>  . . . . . . . . . . . . . . . . . . . . . . .  equality check.
##
RClassOps.\=:= function(l, r)

   if IsRClass(l) and IsRClass(r) then
      return Representative(r) in l;
   else
      return false;
   fi;

end;

#############################################################################
##
#F  OnRClassesAntiAction( <R>, <x> ) . . . . . . .  left action on R classes.
##
OnRClassesAntiAction:= function(R, x)

   return RClass(R.monoid, x * R.representative);

end;

#############################################################################
##
#V  LClassOps  . . . . . . . . . . . . . . . . . . . . . . operations record.
##
LClassOps:= OperationsRecord("LClassOps", DomainOps);

#############################################################################
##
#F  LClass( <M>, <x> ) . . . . . . . . . . . . . . . . . .  define a L class.
##
LClass:= function(M, x)

   return M.operations.LClass(M, x);

end;

MonoidOps.LClass:= function(M, x)

   local lClass;

   lClass:= rec();
   lClass.isLClass:= true;
   lClass.isDomain:= true;
   lClass.monoid:= M;
   lClass.representative:= x;
   lClass.operations:= LClassOps;

   return lClass;

end;

#############################################################################
##
#F  IsLClass( <obj> )  . . . . . . . . . . . . . . . . . . . . .  type check.
##
IsLClass:= function(obj)

   return IsRec(obj) and IsBound(obj.isLClass) and obj.isLClass;

end;

#############################################################################
##
#F  LClassOps.Print( <lClass> )  . . . . . . . . . . . . . . . . . . . print.
##
LClassOps.Print:= function(lClass)

   Print("LClass( ", lClass.monoid, ", ", lClass.representative, " )");

end;
   
#############################################################################
##
#F  <l> = <r>  . . . . . . . . . . . . . . . . . . . . . . .  equality check.
##
LClassOps.\=:= function(l, r)

   if IsLClass(l) and IsLClass(r) then
      return Representative(l) in r;
   else
      return false;
   fi;

end;

#############################################################################
##
#F  OnLClasses( <L>, <x> ) . . . . . . . . . . . . . . . action on L classes.
##
##  Since the L classes of a monoid $M$ form a right (?) congruence, $M$ acts
##  from the right on the set of all L classes.
##
OnLClasses:= function(L, x)

   return LClass(L.monoid, L.representative * x);

end;

#############################################################################
##
#V  HClassOps  . . . . . . . . . . . . . . . . . . . . . . operations record.
##
HClassOps:= OperationsRecord("HClassOps", DomainOps);

#############################################################################
##
#F  HClass( <M>, <x> ) . . . . . . . . . . . . . . . . . .  define a H class.
##
HClass:= function(M, x)

   return M.operations.HClass(M, x);

end;

MonoidOps.HClass:= function(M, x)

   local hClass;

   hClass:= rec();
   hClass.isHClass:= true;
   hClass.isDomain:= true;
   hClass.monoid:= M;
   hClass.representative:= x;
   hClass.operations:= HClassOps;

   return hClass;

end;

#############################################################################
##
#F  IsHClass( <obj> )  . . . . . . . . . . . . . . . . . . . . .  type check.
##
IsHClass:= function(obj)

   return IsRec(obj) and IsBound(obj.isHClass) and obj.isHClass;

end;

#############################################################################
##
#F  HClassOps.Print( <hClass> )  . . . . . . . . . . . . . . . . . . . print.
##
HClassOps.Print:= function(hClass)

   Print("HClass( ", hClass.monoid, ", ", hClass.representative, " )");

end;
   
#############################################################################
##
#F  <l> = <r>  . . . . . . . . . . . . . . . . . . . . . . .  equality check.
##
HClassOps.\=:= function(l, r)

   if IsHClass(l) and IsHClass(r) then
      return Representative(r) in l;
   else
      return false;
   fi;

end;

#############################################################################
##
#V  DClassOps  . . . . . . . . . . . . . . . . . . . . . . operations record.
##
DClassOps:= OperationsRecord("DClassOps", DomainOps);

#############################################################################
##
#F  DClass( <M>, <x> ) . . . . . . . . . . . . . . . . . .  define a D class.
##
DClass:= function(M, x)

   return M.operations.DClass(M, x);

end;

MonoidOps.DClass:= function(M, x)

   local dClass;

   dClass:= rec();
   dClass.isDClass:= true;
   dClass.isDomain:= true;
   dClass.monoid:= M;
   dClass.representative:= x;
   dClass.operations:= DClassOps;

   return dClass;

end;

#############################################################################
##
#F  IsDClass( <obj> )  . . . . . . . . . . . . . . . . . . . . .  type check.
##
IsDClass:= function(obj)

   return IsRec(obj) and IsBound(obj.isDClass) and obj.isDClass;

end;

#############################################################################
##
#F  DClassOps.Print( <dClass> )  . . . . . . . . . . . . . . . . . . . print.
##
DClassOps.Print:= function(dClass)

   Print("DClass( ", dClass.monoid, ", ", dClass.representative, " )");

end;
   
#############################################################################
##
#F  <l> = <r>  . . . . . . . . . . . . . . . . . . . . . . .  equality check.
##
DClassOps.\=:= function(l, r)

   if IsDClass(l) and IsDClass(r) then
      return Representative(r) in l;
   else
      return false;
   fi;

end;

#############################################################################
##
#F  DClasses( <obj> )  . . . . . . . . . . . . . . . . . . . . . . D classes.
##
DClasses:= function(obj)

   if IsRec(obj) then
      if IsBound(obj.dClasses) then
         return obj.dClasses;
      elif IsBound(obj.operations) and IsBound(obj.operations.DClasses) then
         obj.dClasses:= obj.operations.DClasses(obj);
         return obj.dClasses;
      else
         Error("don't know how to construct D Classes of <obj>");
      fi;
   else
      Error("<obj> must be a record");
   fi;

end;

#############################################################################
##
#F  MonoidOps.DClasses( <M> )  . . . . . . . . . . . . . . . . . . D classes.
##
MonoidOps.DClasses:= function(M)

   Error("no generic method for D classes known");

end;

#############################################################################
##
#F  LClasses( <obj> )  . . . . . . . . . . . . . . . . . . . . . . L classes.
##
LClasses:= function(obj)

   if IsRec(obj) then
      if IsBound(obj.lClasses) then
         return obj.lClasses;
      elif IsBound(obj.operations) and IsBound(obj.operations.LClasses) then
         obj.lClasses:= obj.operations.LClasses(obj);
         return obj.lClasses;
      else
         Error("don't know how to construct L Classes of <obj>");
      fi;
   else
      Error("<obj> must be a record");
   fi;

end;

#############################################################################
##
#F  MonoidOps.LClasses( <M> )  . . . . . . . . . . . . . . . . . . L Classes.
##
MonoidOps.LClasses:= function(M)

   local D, lcl;

   # initialize.
   lcl:= [];

   # loop over the D classes.
   for D in DClasses(M) do 
      Append(lcl, LClasses(D));
   od;

   # return the list of L classes.
   return lcl;

end;

#############################################################################
##
#F  RClasses( <obj> )  . . . . . . . . . . . . . . . . . . . . . . R classes.
##
RClasses:= function(obj)

   if IsRec(obj) then
      if IsBound(obj.rClasses) then
         return obj.rClasses;
      elif IsBound(obj.operations) and IsBound(obj.operations.RClasses) then
         obj.rClasses:= obj.operations.RClasses(obj);
         return obj.rClasses;
      else
         Error("don't know how to construct R Classes of <obj>");
      fi;
   else
      Error("<obj> must be a record");
   fi;

end;

#############################################################################
##
#F  MonoidOps.RClasses( <M> )  . . . . . . . . . . . . . . . . . . R Classes.
##
MonoidOps.RClasses:= function(M)

   local D, rcl;

   # initialize.
   rcl:= [];

   # loop over the D classes.
   for D in DClasses(M) do
      Append(rcl, RClasses(D));
   od;

   #  return the list of R classes.
   return rcl;

end;

#############################################################################
##
#F  HClasses( <obj> )  . . . . . . . . . . . . . . . . . . . . . . H classes.
##
HClasses:= function(obj)

   if IsRec(obj) then
      if IsBound(obj.hClasses) then
         return obj.hClasses;
      elif IsBound(obj.operations) and IsBound(obj.operations.HClasses) then
         obj.hClasses:= obj.operations.HClasses(obj);
         return obj.hClasses;
      else
         Error("don't know how to construct H Classes of <obj>");
      fi;
   else
      Error("<obj> must be a record");
   fi;

end;

#############################################################################
##
#F  MonoidOps.HClasses( <M> )  . . . . . . . . . . . . . . . . . . H Classes.
##
MonoidOps.HClasses:= function(M)

   local D, hcl;

   # initialize.
   hcl:= [];

   # loop over the D classes.
   for D in DClasses(M) do 
      Append(hcl, HClasses(D));
   od;

   # return the list of H classes.
   return hcl;

end;

#############################################################################
##
#F  SchutzenbergerGroup( <xClass> )  . . . . . . . .  Sch\"utzenberger group.
#F  SchutzenbergerGroup( <M>, <x> )  . . . . . . . .  Sch\"utzenberger group.
##
SchutzenbergerGroup:= function(arg)

   local  class, schGrp;

   # the argument might be an H class, or...
   if Length(arg) = 1 then
      class:= arg[1];
      if IsRec(class) and IsBound(class.schutzenbergerGroup) then
         schGrp:= class.schutzenbergerGroup;
      elif IsBound(class.operations) and 
           IsBound(class.operations.SchutzenbergerGroup) then
         schGrp:= class.operations.SchutzenbergerGroup(class);
         class.schutzenbergerGroup:= schGrp;
      else
         Error("don't know how to compute Schutzenberger group of <class>");
      fi;

   # the arguments are a monoid and one of its elements...
   elif Length(arg) = 2 and IsMonoid(arg[1]) then
      schGrp:= SchutzenbergerGroup(HClass(arg[1], arg[2]));

   # ...and nothing else.
   else
      Error("usage: SchutzenbergerGroup( <hClass> ) or\n",
     "              SchutzenbergerGroup( <M>, <x> )");
   fi;

   # return the group.
   return schGrp;

end;                     

#############################################################################
ReadPkg("monoid", "lib", "action");

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
