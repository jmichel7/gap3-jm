#############################################################################
##
#A  GAP                                                 Goetz.Pfeiffer@UCG.IE
##                                  & Werner Nickel <werner@dcs.st-and.ac.uk>
##
#A  $Id: relation.g,v 2.1 1997/09/05 11:44:11 goetz Exp $
##
#Y  Copyright (C) 1997, Mathematics Dept, University College Galway, Ireland.
##
##  This  file  defines the  domain of and   the functions for  finite binary
##  relations.
##
##  This should be the only file that  knows about (and changes) the internal
##  representation of  a binary relation,  everyone  else can  only   ask for
##  properties, multiply or invert them.
##
##  A finite binary relation on $n$ points is a multivalued  map from the set
##  $[1..n]$ to itself.  It is represented by a list of succssors sets.  Each
##  successors set is  internally represented as  a  subset of $[1..n]$ by  a
##  boolean list  of  length $n$.  All binary   relations on $[1..n]$ form  a
##  monoid.
##

#############################################################################
##
#V  RelationsOps . . . . . . . . . . . . . . . . . . . . . operations record.
##
##  The domain 'Relations' inherits from 'MonoidElements'.
##
RelationsOps:= OperationsRecord("RelationsOps", MonoidElementsOps);

#############################################################################
##
#V  Relations  . . . . . . . . . . . . . . . . . . . . . . . . . . .  domain.
##
Relations:= rec(); 
Relations.isDomain:= true;
Relations.name:= "Relations";
Relations.isFinite:= false;
Relations.size:= "infinity";
Relations.operations:= RelationsOps;

#############################################################################
##
#F  IsRelation( <obj> )  . . . . . . . . . . . . . . . . . . . .  type check.
##
IsRelation:= obj-> IsRec(obj) and IsBound(obj.isRelation) and obj.isRelation;

#############################################################################
##
#F  <obj> in Relations . . . . . . . . . . . . . . . . . . . membership test.
##
RelationsOps.\in:= function(obj, Relations) return IsRelation(obj);  end;

#############################################################################
##
#V  RelationOps  . . . . . . . . . . . . . . . . . . . . . operations record.
##
##  Every 'Relation' inherits from 'Mapping'.
##
RelationOps:= OperationsRecord("RelationOps", MappingOps);

#############################################################################
##
#F  Relation( <lst> )  . . . . . . . . . . . . . . construct binary relation.
##
##  A binary  relation  is defined by  its  list of  image sets  <lst>.   The
##  argument <lst> must be a list of sets (i.e. sorted!) of integers, each of
##  which must lie in the range '[1..Length(<lst>)]'.
##
##  Internally, each  image set is  stored as a  boolean  list of length $n$.
##  The  range is not  checked by 'BlistList',  but unwanted  input is simply
##  ignored.
##
Relation:= function(lst)

   local rel, len, i, rng;

   # initialize.
   rel:= rec();

   # declare domain.
   rel.isRelation:= true;
   rel.domain:= Relations;

   # immediate case.
   if IsBlist(lst[1]) and Length(lst[1]) = Length(lst) then
      rel.successors:= lst;

   # conversion case.
   else
      len:= Length(lst);  rng:= [1..len];  rel.successors:= [];
      for i in rng do
         rel.successors[i]:= BlistList(rng, lst[i]);
      od;
   fi;

   # install operations.
   rel.operations:= RelationOps;

   # return binary relation.
   return rel;

end;

#############################################################################
##
#F  RandomRelation( <n> )  . . . . . . . . . . . . .  random binary relation
##
##  returns a random binary relation on <n> points.
##
RandomRelation:= function(n)

   local i, j, booleans, rel;

   if not IsInt(n) then
      return Error("usage: RandomRelation( <number of points> )");
   fi;

   # initialize.
   booleans:= [false, false, true, false, false];
   rel:= List( [1..n], x->[] );

   # loop over pairs.
   for i in [1..n] do
      for j in [1..n] do
          # throw the coin.
          rel[i][j] := Random(booleans);
      od;
      IsBlist( rel[i] );
   od;

   # return the random relation.
   return Relation(rel);

end;

#############################################################################
##
#F  RelationOps.Print( <rel> ) . . . . . . . . . . . . . . . . . . . . print.
##
RelationOps.Print:= function(rel)

   local img, len, rng, i;

   img:= rel.successors;  len:= Length(img);  rng:= [1..len];
   Print("Relation( [ ");
   for i in [1..len-1] do
      Print(ListBlist(rng, img[i]), ", ");
   od;
   Print(ListBlist(rng, img[len]), " ] )");   

end;

#############################################################################
##
#F  RelationOps.Kernel( <rel> )  . . . . . . . . . . . . . . . . . .  kernel.
##
RelationOps.Kernel:= function(rel)

   Error("not yet implemented");

end;

#############################################################################
##
#F  RelationOps.Degree( <rel> )  . . . . . . . . . . . . . . . . . .  degree.
##
RelationOps.Degree:= rel -> Length(rel.successors);

#############################################################################
##
#F  <l> = <r>  . . . . . . . . . . . . . . . . . . . . . . .  equality check.
##
##  Two binary relations are equal if their image lists are equal.
##
RelationOps.\=:= function(l, r)

   if IsRelation(l) and IsRelation(r) then
      return l.successors = r.successors;
   else
      return false;
   fi;

end;

#############################################################################
##
#F  <l> < <r>  . . . . . . . . . . . . . . . . . . . . . . . . .  comparison.
##
##  makes Relations bigger than everything else.
##
RelationOps.\<:= function(l, r)

   if IsRelation(l) then
      if IsRelation(r) then
         return l.successors < r.successors;
      else
         return false;
      fi;
   else
      if IsRelation(r) then
         return true;
      else
         Error("Panic: either <l> or <r> should be a Relation");
      fi;
   fi;

end;

#############################################################################
##
#F  RelTrans( <trans> )  . . . . . . . . . . . . . . . . . . type conversion.
##
##  returns the binary relation defined by the transformation <trans>.
##
RelTrans:= function(trans)

   local i, img, rng, new;

   # initialize.
   img:= trans.images;  rng:= [1..Length(img)];  new:= [];

   # loop over the range.
   for i in rng do
      new[i]:= BlistList(rng, []);
      new[i][img[i]]:= true;
   od;

   # return the binary relation.
   return Relation(new);

end;

#############################################################################
##
#F  TransRel( <rel> )  . . . . . . . . . . . . . . . . . . . type conversion.
##
##  returns the  transformation defined by  the binary relation <rel>.  This
##  should only be applied if every image of <rel> has size 1.
##
TransRel:= function(rel)

   local img, new;

   # initialize.
   new:= [];

   # loop over successors.
   for img in rel.successors do 
      if SizeBlist(img) = 1 then
         Add(new, Position(img, true));
      else
         Error("<rel> is not a transformation");
      fi;
   od;

   # return the transformation
   return TransformationNC(new);

end;

#############################################################################
##
#F  <l> * <r>  . . . . . . . . . . . . . . . . . . . . . . . . . . . product.
##
##  Note that by  this definition of  a product Relations act  from the
##  right.
##
RelationOps.ProductRelRel:= function(l, r)

   local i, j, rng, img, new;

   # initialize.
   img:= l.successors;  rng:= [1..Length(img)];  new:= [];

   # loop over the range.
   for i in rng do
      new[i]:= BlistList(rng, []);
      for j in rng do
         if img[i][j] then
            UniteBlist(new[i], r.successors[j]);
         fi;
      od;
   od;

   # return the product.
   return Relation(new);

end;

RelationOps.ProductRelTrans:= function(l, r)

   local i, j, rng, img, new;

   # initialize.
   img:= l.successors;  rng:= [1..Length(img)];  new:= [];

   # loop over the range.
   for i in rng do
      new[i]:= BlistList(rng, []);
      for j in rng do
         if img[i][j] then
            new[i][r.images[j]]:= true;
         fi;
      od;
   od;

   # return the product.
   return Relation(new);

end;

RelationOps.ProductRelPerm:= function(l, r)
  
   local rng, img, new;

   # initialize.
   img:= l.successors;  rng:= [1..Length(img)]; 

   # loop over successors. ;-)
   new:= img{rng}{OnTuples(rng, r^-1)};

   # return the product.
   return Relation(new);

end;

RelationOps.ProductPermRel:= function(l, r)

   local rng, img, new;

   # initialize.
   img:= r.successors;  rng:= [1..Length(img)]; 

   # loop over successors.
   new:= img{OnTuples(rng, l)};

   # return the product.
   return Relation(new);

end;

RelationOps.\*:= function(l, r)

   if IsRelation(l) then

      # Relation by Relation.
      if IsRelation(r) then
         return RelationOps.ProductRelRel(l, r);

      # Relation by perm.
      elif IsPerm(r) then
         return RelationOps.ProductRelPerm(l, r);

      # Relation by list.
      elif IsList(r) then
         return List(r, x-> l * x);

      else
         Error("don't know how to multiply Relation <l> by <r>");
      fi;

   elif IsRelation(r) then

      # perm by Relation.
      if IsPerm(l) then
         return RelationOps.ProductPermRel(l, r);

      # trans by Relation.
      elif IsTransformation(l) then

         #C There must be a more efficient way...
         return RelTrans(l) * r;

      # list by Relation.
      elif IsList(l) then
         return List(l, x-> x * r);

      else
         Error("don't know how to multiply <l> by Relation");
      fi;

  else
     Error("Panic: neither <l> nor <r> is a Relation!");
  fi;

end;

#############################################################################
##
#F  IdentityRelation( <n> )  . . . . . . . . . . . . . . . . . . .  identity.
##
IdentityRelation:= function(n)

   local i, new, rng;

   # initialize.
   rng:= [1..n];  new:= [];

   # loop over the range.
   for i in rng do
      new[i]:= BlistList(rng, []);
      new[i][i]:= true;
   od;

   # return the identity.
   return Relation(new);

end;

#############################################################################
##
#F  EmptyRelation( <n> ) . . . . . . . . . . . . . . . . . .  empty relation.
##
EmptyRelation:= function(n)

   local i, new, rng;

   # initialize.
   rng:= [1..n];  new:= [];

   # loop over the range.
   for i in rng do
      new[i]:= BlistList(rng, []);
   od;

   # return the relation.
   return Relation(new);

end;

#############################################################################
##
#F  InverseRelation( <rel> ) . . . . . . . . . . . . . . . . . . the inverse.
##
##  The inverse of a binary relation <rel> arises from <rel> by ``reverting
##  arrows''.   Note that, in  general, the product  of a binary relation and
##  its inverse does  not equal the  identity.   Neither is  it equal to  the
##  product of the inverse and the binary relation.
##
InverseRelation:= function(rel)

   local i, j, img, new, rng;

   # initialize.
   img:= rel.successors;  rng:= [1..Length(img)];  new:= [];
   for i in rng do
      new[i]:= BlistList(rng, []);
   od;

   # loop over image.
   for i in rng do
      for j in rng do 
         if img[i][j] then
            new[j][i]:= true;
         fi;
      od;
   od;

   # return the inverse.
   return Relation(new);

end;

#############################################################################
##
#F  InverseTransformation( <trans> ) . . . . . . . . . . . . . . the inverse.
##
##  is a shorthand for 'InverseRelation(RelationTrans( <trans> )'
##  and  returns  the inverse  of   the  transformation <trans> as  a  binary
##  relation.
##
InverseTransformation:= function(trans)

   local i, j, img, new, rng;

   # initialize.
   img:= trans.images;  rng:= [1..Length(img)];  new:= [];
   for i in rng do
      new[i]:= BlistList(rng, []);
   od;

   # loop over image.
   for i in rng do
      new[img[i]][i]:= true;
   od;

   # return the inverse.
   return Relation(new);

end;

#############################################################################
##
#F  <l> ^ <r>  . . . . . . . . . . . . . . . . . . . . . . .  exponentiation.
##
RelationOps.\^:= function(l, r)

   #  power of Relation.
   if IsRelation(l) then
      if IsInt(r) then

         # the one.
         if r = 0 then
            return IdentityRelation(Length(l.successors));
    
         # identity.
         elif r = 1 then
            return l;
    
         # the square.
         elif r = 2 then
            return l * l;
   
         #C this should better be done iteratively:
         # higher powers.
         elif r > 2 then
            if r mod 2 = 1 then
               return (l^((r-1)/2))^2 * l;
            else
               return (l^(r/2))^2;
            fi;
    
         # inverse.
         elif r = -1 then
            return InverseRelation(l);
   
         else  # r < -1
            return InverseRelation(l^-r);
    
         fi;

      else
         Error("don't know how to form power of binary relation by <r>");
      fi;

   elif IsRelation(r) then
 
      # image of <int> under <rel>.
      if IsInt(l) then
         return ListBlist([1..Length(r.successors)], r.successors[l]);

      # image of <set> under <rel>.
      elif IsSet(l) then
         return Union(List(l, x-> x^r));
 
      # image of <list> under <rel>.
      elif IsList(l) then
         return List(l, x-> x^r);
 
      else
         Error("don't know how to form power of <l> by transformation");
      fi;

   # this should never happen! 
   else
      Error("PANIC: neither <l> nor <r> is a transformation!");
   fi;

end;

#############################################################################
##
#F  IsReflexive( <obj> ) . . . . . . . . . . . . . . . . . . . property test.
##
IsReflexive:= function(obj)

   if IsRec(obj) then
      if not IsBound(obj.isReflexive) then
         if IsBound(obj.operations) 
             and IsBound(obj.operations.IsReflexive) then
            obj.isReflexive:= obj.operations.IsReflexive(obj);
         else
            Error("don't know how to determine property");
         fi;
      fi;
      return obj.isReflexive;
   else
      return false;
   fi;

end;

#############################################################################
##
#F  RelationOps.IsReflexive( <rel> ) . . . . . . . . . . . . . property test.
##
RelationOps.IsReflexive:= function(rel)

   local img, i;

   img:= rel.successors;
   for i in [1..Length(img)] do
      if not img[i][i] then
         return false;
      fi;
   od;

   return true;

end;

#############################################################################
##
#F  ReflexiveClosure( <rel> )  . . . . . . . . . . . . . . reflexive closure.
##
ReflexiveClosure:= function(rel)

   local img, new, rng, i;

   # initialize.
   img:= rel.successors;  rng:= [1..Length(img)];  new:= [];

   # loop over the successors.
   for i in rng do
      new[i]:= Copy(img[i]);
      new[i][i]:= true;
   od;

   # make a binary relation.
   new:= Relation(new);
   new.isReflexive:= true;

   # return the result.
   return new;

end;

#############################################################################
##
#F  IsSymmetric( <obj> ) . . . . . . . . . . . . . . . . . . . property test.
##
IsSymmetric:= function(obj)

   if IsRec(obj) then
      if not IsBound(obj.isSymmetric) then
         if IsBound(obj.operations) 
             and IsBound(obj.operations.IsSymmetric) then
            obj.isSymmetric:= obj.operations.IsSymmetric(obj);
         else
            Error("don't know how to determine property");
         fi;
      fi;
      return obj.isSymmetric;
   else
      return false;
   fi;

end;

#############################################################################
##
#F  RelationOps.IsSymmetric( <rel> ) . . . . . . . . . . . . . property test.
##
RelationOps.IsSymmetric:= function(rel)

   local rng, img, i, j;

   img:= rel.successors;  rng:= [1..Length(img)];
   for i in rng do
      for j in rng do
         if img[i][j] and not img[j][i] then
            return false;
         fi;
      od;
   od;

   return true;

end;

#############################################################################
##
#F  SymmetricClosure( <rel> )  . . . . . . . . . . . . . . symmetric closure.
##
SymmetricClosure:= function(rel)

   local img, new, rng, i, j;

   # initialize.
   img:= rel.successors;  rng:= [1..Length(img)];  new:= Copy(img);

   # loop over successors.
   for i in rng do
      for j in rng do
         if img[i][j] then
            new[j][i]:= true;
         fi;
      od;
   od;

   # make binary relation.
   new:= Relation(new);
   new.isSymmetric:= true;

   # return the closure.
   return new;

end;

#############################################################################
##
#F  IsTransitiveRel( <obj> ) . . . . . . . . . . . . . . . . . property test.
##
#C  should be called 'IsTransitive'.  But this name is already occupied by
#C  perm groups!!
##
IsTransitiveRel:= function(obj)

   if IsRec(obj) then
      if not IsBound(obj.isTransitive) then
         if IsBound(obj.operations) 
             and IsBound(obj.operations.IsTransitive) then
            obj.isTransitive:= obj.operations.IsTransitive(obj);
         else
            Error("don't know how to determine property");
         fi;
      fi;
      return obj.isTransitive;
   else
      return false;
   fi;

end;

#############################################################################
##
#F  RelationOps.IsTransitive( <rel> )  . . . . . . . . . . . . property test.
##
RelationOps.IsTransitive:= function(rel)

   local img, rng, i, j;

   # initialize.
   img:= rel.successors;  rng:= [1..Length(img)];

   # loop over range.
   for i in rng do
      for j in rng do
         if img[i][j] and not IsSubsetBlist(img[i], img[j]) then
            return false;
         fi;
      od;
   od;

   # if nothing fails.
   return true;

end;

#############################################################################
##
#F  TransitiveClosure( <rel> ) . . . . . . . . . . . . .  transitive closure.
##
##  Warshall's algorithm.
##
RelationOps.TransitiveClosure:= function(rel)

   local new;

   new:= Relation(TransitiveClosure(rel.successors));
   new.IsTransitive:= true;

   return new;

end;

#############################################################################
##
#F  IsAntisymmetric( <obj> ) . . . . . . . . . . . . . . . .  property check.
##
IsAntisymmetric:= function(obj)

   if IsRec(obj) then
      if not IsBound(obj.isAntisymmetric) then
         if IsBound(obj.operations) 
             and IsBound(obj.operations.IsAntisymmetric) then
            obj.isAntisymmetric:= obj.operations.IsAntisymmetric(obj);
         else
            Error("don't know how to determine property");
         fi;
      fi;
      return obj.isAntisymmetric;
   else
      return false;
   fi;

end;

#############################################################################
##
#F  RelationOps.IsAntisymmetric( <rel> ) . . . . . . . . . .  property check.
##
RelationOps.IsAntisymmetric := function( rel )
    local   n,  i,  j;

    rel:= rel.successors;
    
    n := Length( rel );
    for i in [1..n] do
        for j in [i+1..n] do
            if rel[i][j] and rel[j][i] then
                return false;
            fi;
        od;
    od;
    
    return true;
end;

#############################################################################
##
#F  IsPreOrder( <rel> )  . . . . . . . . . . . . . . . . . . . property test.
##
IsPreOrder:= function(rel)

   if not IsBound(rel.isPreOrder) then
      rel.isPreOrder:= IsReflexive(rel) and IsTransitiveRel(rel);
   fi;

   return rel.isPreOrder;

end;

#############################################################################
##
#F  IsPartialOrder( <rel> )  . . . . . . . . . . . . . . . . . property test.
##
IsPartialOrder:= function(rel)

   if not IsBound(rel.isPartialOrder) then
      rel.isPartialOrder:= IsPreOrder(rel) and IsAntisymmetric(rel);
   fi;

   return rel.isPartialOrder;

end;

#############################################################################
##
#F  IsEquivalence( <rel> ) . . . . . . . . . . . . . . . . . . property test.
##
IsEquivalence:= function(rel)

   if not IsBound(rel.isEquivalence) then
      rel.isEquivalence:= IsPreOrder(rel) and IsSymmetric(rel);
   fi;

   return rel.isEquivalence;

end;

#############################################################################
##
#F  EquivalenceClasses( <rel> )  . . . . . . . . . . . . equivalence classes.
##
EquivalenceClasses:= function(rel)

   local img, rng, cls;

   if not IsEquivalence(rel) then
      Error("<rel> must be an equivalence relation");
   fi;

   img:= rel.successors;  rng:= [1..Length(img)];  cls:= Set(rel.successors);

   return List(cls, x-> ListBlist(rng, x));

end;

#############################################################################
##
#F  HasseDiagram( <rel> )  . . . . . . . . . . . . . . . . . . Hasse diagram.
##
HasseDiagram := function( rel )
    local   n,  Hr,  i,  j;
    
    if not IsPartialOrder( rel ) then
        return Error( "usage: HasseDiagram( <partial ord> )" );
    fi;
    
    rel:= rel.successors;
    n := Length( rel );
    Hr := Copy( rel );
    
    # make Hr non-reflexive
    for i in [1..n] do
        Hr[i][i] := false;
    od;
    
    for i in [1..n] do
        # for each point keep only those neighbours that can't be reached
        # through any other neighbour.
        for j in [1..n] do
            if Hr[i][j] then
                SubtractBlist( Hr[i], Hr[j] );
            fi;
        od;
    od;
    
    return Relation(Hr);
end;

#############################################################################
##
#F  Monoid( Relations, <gens> )  . . . . . . . . . . . .  construct a monoid.
##
RelationsOps.Monoid:= function(Relations, gens, id)

   local s, M, r, t;

   # determine and check degree.
   r:= Degree(id);
   if not ForAll(gens, x-> Degree(x) = r) then
      Error("<gens> must have the same degree");
   fi;

   M:= MonoidElementsOps.Monoid(Relations, gens, id);
   M.isRelMonoid:= true;
   M.isFinite:= true;
   M.operations:= RelMonoidOps;

   # return the monoid.
   return M;

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

