#############################################################################
##
#A  GAP                                                 Goetz.Pfeiffer@UCG.IE
##
#A  $Id: transfor.g,v 2.2 1997/09/05 11:53:44 goetz Exp $
##
#Y  Copyright (C) 1997, Mathematics Dept, University College Galway, Ireland.
##
##  This file defines the domain of and functions for transformations.
##  
##  Let $n >  0$.  A *transformation*  of degree  $n$ is a  map  from the set
##  $[1..n]$ into itself.  The transformations of degree $n$ form monoid: the
##  full transformation monoid of degree $n$.
##

#############################################################################
##
#V  TransformationsOps . . . . . . . . . . . . . . . . . . operations record.
##
##  The domain 'Transformations' inherits from 'MonoidElements'.
##
TransformationsOps:= OperationsRecord("TransformationsOps", 
                                      MonoidElementsOps);

#############################################################################
##
#F  IsTransformation( <obj> )  . . . . . . . . . . . . . . . . .  type check.
##
IsTransformation:= function(obj)
   return IsRec(obj) and IsBound(obj.isTransformation) 
                     and obj.isTransformation;
end;

#############################################################################
##
#V  Transformations  . . . . . . . . . . . . . . . . . . . . . . . .  domain.
##
Transformations:= rec(); 
Transformations.isDomain:= true;
Transformations.name:= "Transformations";
Transformations.isFinite:= false;
Transformations.size:= "infinity";
Transformations.operations:= TransformationsOps;

#############################################################################
##
#F  <obj> in Transformations . . . . . . . . . . . . . . . . membership test.
##
TransformationsOps.\in:= function(obj, Transformations)
   return IsTransformation(obj);
end;

#############################################################################
##
#V  TransformationOps  . . . . . . . . . . . . . . . . . . operations record.
##
TransformationOps:= OperationsRecord("TransformationOps", MappingOps);

#############################################################################
##
#F  Transformation( <lst> )  . . . . . . . . . . .  construct transformation.
##
##  A transformation is defined by  its list of images.  This  list must be a
##  vector  of positive integers,  where no entry  exceeds  the length of the
##  list.  
##
##  'TransformationNC' is the version without having the arguments checked.
##
TransformationNC:= function(lst)

   local trans;

   # initialize.
   trans:= rec();

   # install components
   trans.isTransformation:= true;
   trans.domain:= Transformations;
   trans.images:= lst;
   trans.operations:= TransformationOps;

   # return the transformation.
   return trans;

end;

Transformation_mon:= function(lst)

   local i, deg;

   # check argument.
   if not IsList(lst) then 
      Error("<lst> must be a list");
   fi;

   deg:= Length(lst);
   for i in lst do
      if not (IsInt(i) and i > 0 and i <= deg) then
         Error("image must not be > ", deg);
      fi;
   od;

   # construct transformation.
   return TransformationNC(lst);

end;
      

#############################################################################
##
#F  TransPerm( <n>, <perm> ) . . . . . . . . . . . . . . . . type conversion.
##
##  turns the  permutation <perm> into  a transformation  of  degree <n>.  
##
#W  Warning:   No   error    is  issued     if      <n>   is  smaller    that
#W  'LargestMovedPoint(<perm>)'.   This   can  cause  the   result  to be  no
#W  transformation at all!!!
##
TransPerm:= function(n, perm)

   return TransformationNC(OnTuples([1..n], perm));

end;

#############################################################################
##
#F  PermTrans( <trans> ) . . . . . . . . . . . . . . . . . . type conversion.
##
##  returns the permutation defined by the transformation <trans>.
##
##  'PermList' *will* check that <trans> is a bijection.
##
PermTrans:= trans-> PermList(trans.images);

#############################################################################
##
#F  IdentityTransformation ( <n> ) . . . . . . . . . identity transformation.
##
IdentityTransformation_mon:= n-> TransformationNC([1..n]);

#############################################################################
##
#F  TransformationOps.Print( <trans> ) . . . . . . . . . . . . . . . . print.
##
TransformationOps.Print:= function(trn)

   Print("Transformation( ", trn.images, " )");

end;

#############################################################################
##
#F  <l> = <r>  . . . . . . . . . . . . . . . . . . . . . . .  equality check.
##
TransformationOps.\=:= function(l, r)

   if IsTransformation(l) and IsTransformation(r) then
      return l.images = r.images;
   else
      return false;
   fi;

end;

#############################################################################
##
#F  <l> < <r>  . . . . . . . . . . . . . . . . . . . . . . . . .  comparison.
##
##  makes transformations bigger than everything else.
##
TransformationOps.\<:= function(l, r)

   if IsTransformation(l) then
      if IsTransformation(r) then
         return l.images < r.images;
      else
         return false;
      fi;
   else
      if IsTransformation(r) then
         return true;
      else
         Error("either <l> or <r> should be a Transformation");
      fi;
   fi;

end;


#############################################################################
##
#F  <l> * <r>  . . . . . . . . . . . . . . . . . . . . . . . . . . . product.
##
##  Note that by  this definition of  a product transformations act  from the
##  right.
##
TransformationOps.\*:= function(l, r)

   if IsTransformation(l) then

      # transformation by transformation.
      if IsTransformation(r) then
         return TransformationNC(r.images{l.images});              # :-)

      # transformation by perm.
      elif IsPerm(r) then
         return TransformationNC(OnTuples(l.images, r));

      # transformation by list.
      elif IsList(r) then
         return List(r, x-> l * x);

      else
         Error("don't know how to multiply transformation by <r>");
      fi;

   elif IsTransformation(r) then

      # perm by transformation.
      if IsPerm(l) then
        return TransformationNC(r.images{OnTuples([1..Length(r.images)],l)});

      # Relation by Transformation.
      elif IsRelation(l) then  
         return RelationOps.ProductRelTrans(l, r);

      # list by transformation.
      elif IsList(l) then
         return List(l, x-> x * r);

      else
         Error("don't know how to multiply <l> by transformation");
      fi;

  # this should never happen.
  else
     Error("Panic: neither <l> nor <r> is a transformation!");
  fi;

end;

#############################################################################
##
#F  <l> ^ <r>  . . . . . . . . . . . . . . . . . . . . . . .  exponentiation.
##
TransformationOps.\^:= function(l, r)

   #  power of transformation
   if IsTransformation(l) then
      if IsInt(r) then

         # the one.
         if r = 0 then
            return IdentityTransformation(Length(l.images));
    
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
            return InverseTransformation(l);
   
         else  # r < -1
            return InverseTransformation(l^-r);
    
         fi;
 
      else
         Error("don't know how to form power of transformation by <r>");
      fi;
 
   elif IsTransformation(r) then
 
      #  image of <int> under <trans>.
      if IsInt(l) then
         return r.images[l];
    
      #  image of <vec> under <trans>.
      elif IsVector(l) then
         return r.images{l};
    
      #  image of <list> under <trans>.
      elif IsList(l) then
         return List(l, x-> x^r);

      else
         Error("don't know how to form power of <l> by transformation");
      fi;
 
   else
      Error("PANIC: neither <l> nor <r> is a transformation!");
   fi;

end;

#############################################################################
##
#F  TransformationOps.ImagesSource( <trans> )  . . . . . . . . . . .  images.
##
TransformationOps.ImagesSource:= function(trans)

   return  Set(trans.images);

end;

#############################################################################
##
#F  TransformationOps.PreImagesElm( <trans>, <elt> ) . . . . . . . preimages.
##
##  returns the  set of preimages of an element <elt>  under a  map <trans> which is  not necessarily
##  bijective.
##
TransformationOps.PreImagesElm:= function(trans, elt)

   local i, map, pre;

   # initialize.
   map:= trans.images; pre:= [];

   # construct preimage.
   for i in [1..Length(map)] do
      if map[i] = elt then
	 Add(pre, i);
      fi;
   od;

   # return preimage.
   return pre;

end;

#############################################################################
##
#F  TransformationOps.PreImagesRange( <trans> )  . . . . . . . . . preimages.
##
##  List of single preimages. (amounts to the inverse of trans!)
##
TransformationOps.PreImagesRange:= function(trans)

   local map;

   map:= trans.images;
   return List([1..Length(map)], x-> PreImage(trans, x));

end;

#############################################################################
##
#F  TransformationOps.Kernel( <trans> )  . . . . . . . . . . . . . .  kernel.
##
TransformationOps.Kernel:= function(trans)

   local ker, imgs, i;

   # initialize.
   ker:= []; imgs:= trans.images;
   for i in imgs do 
      ker[i]:= [];
   od;

   # compute preimages.
   for i in [1..Length(imgs)] do
      Add(ker[imgs[i]], i);
   od;

   # return kernel.
   return Set(ker);

end;

#############################################################################
##
#F  TransformationOps.Degree( <trans> )  . . . . . . . . . . . . . .  degree.
##
##  The *degree* of a transformation is the number of points it acts upon.
##
#N  will someone please move 'Degree' dispatcher out of "polynom.g"?
##
TransformationOps.Degree:= function(trans)

   return Length(trans.images);

end;

#############################################################################
##
#F  TransformationOps.Rank( <trans> )  . . . . . . . . . . . . . . . .  rank.
##
##  The *rank* of a transformation is the size of its set of images.
##
TransformationOps.Rank:= function(trans)

  return Size(Set(trans.images));

end;

#############################################################################
##
#F  TransformationOps.Restricted( <trans>, <alpha> ) . . . . . . .  restrict.
##
##  restricts <trans> to subset <alpha>.  Do we  need to check arguments, put
##  restrictions on them?
##
#C  this should rather return a (binary) relation defined *only* on <alpha> 
#C  rather than supposing the identity outside of <alpha>.
##
TransformationOps.Restricted:= function(argv)

   local new, trans, alpha;

   trans:= argv[1]; alpha:= argv[2];

   new:= [1..Degree(trans)];
   new{alpha}:= alpha^trans;

   return TransformationNC(new);

end;

#############################################################################
##
#F  TransformationOps.IsInjective( <trans> ) . . . . . . . injectivity check.
##
##  Maybe there should be more of this: after all, a transformation is a map!
##
TransformationOps.IsInjective:= function(trans)

   return Size(Image(trans)) = Degree(trans);

end;

TransformationOps.IsSurjective:= TransformationOps.IsInjective;
TransformationOps.IsBijective:= TransformationOps.IsInjective;
TransformationOps.IsBijection:= TransformationOps.IsInjective;

#############################################################################
##
#F  PermLeftQuoTrans( <tra1>, <tra2> ) . . . . . . . . . . . . . .  quotient.
##
##  Given transformations <tra1> and <tra2> with  equal kernel and image, the
##  permutation induced by <tra1>^-1  * <tra2> on  the set Image( <tra1> ) is
##  computed.
##
##  This high speed implementation is used in time critical parts, it heavily
##  relies on the representation of transformations as lists of images.
##
PermLeftQuoTrans:= function(l, r)

   local n, alpha, perm;

   # construct cross section of kernel.
   n:= Length(l.images);  alpha:= [];
   alpha{l.images}:= [1..n];  # :-) is this a bug or a feature???
   alpha:= Set(alpha);

   # calculate induced permutation.  :-))
   perm:= [1..n];  perm{l.images{alpha}}:= r.images{alpha};  

   # return permutation.
   return PermList(perm);

end;

#############################################################################
##
#F  MappingTransformation( <trans> ) . . . . . . . . . . . . type conversion.
##
##  seems to be useless :-(
##
MappingTransformation:= function(trans)

   local n;

   n:= [1..Length(trans.images)];

   return MappingByFunction( n, n, i-> trans.images[i] );

end;

#############################################################################
##
#F  SemiGroup( Transformations, <gens> ) . . . . . . . construct a semigroup.
##
TransformationsOps.SemiGroup:= function(Transformations, gens)

   local S;

   S:= MonoidElementsOps.SemiGroup(Transformations, gens);
   S.isTransSemiGroup:= true;
   S.isFinite:= true;
   S.operations:= TransSemiGroupOps;

   return S;

end;

#############################################################################
##
#F  Monoid( Transformations, <gens> )  . . . . . . . . .  construct a monoid.
##
TransformationsOps.Monoid:= function(Transformations, gens, id)

   local s, M, r, t;

   # determine and check rank.
   r:= Length(id.images);
   if not ForAll(gens, x-> Length(x.images) = r) then
      Error("<gens> must have the same rank");
   fi;

   M:= MonoidElementsOps.Monoid(Transformations, gens, id);
   M.isTransMonoid:= true;
   M.isFinite:= true;
   M.operations:= TransMonoidOps;

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
