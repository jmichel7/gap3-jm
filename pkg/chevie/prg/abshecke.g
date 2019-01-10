#############################################################################
##
#A  abshecke.g                CHEVIE library                      Jean Michel 
##
#Y  Copyright (C) 1992 - 1999  Lehrstuhl D fuer Mathematik, RWTH Aachen, and
#Y  Universite Paris VII.
##
##  This file contains  GAP functions for working with Hecke algebras for
##  arbitrary reflection groups or cosets.
##

##############################################################################
##
## AbsHeckeOps: operations for Hecke algebras.
##

AbsHeckeOps:=OperationsRecord("AbsHeckeOps");

AbsHeckeOps.CompactPara:=function(p) p:=ShallowCopy(p);
  if Length(Set(p))=1 then p:=[p[1]];fi;
  p:=Zip(p,List(p,Length),function(p,e)
    if IsList(p) and p{[2..Length(p)]}=List([1..e-1],j->E(e)^j) 
    then return p[1];else return p;fi;end);
  if Length(p)=1 and not IsList(p[1]) then p:=p[1];fi;
  return p;
end;

AbsHeckeOps.Format:=function(H,opt)local W,i,n,para,root,res,format;
  format:=function(l)
    if IsList(l) then 
      if IsBound(opt.TeX) then 
         if ForAny(l,IsList) then return Join(List(l,format),";");
         else return Join(List(l,format),",\\allowbreak ");
	 fi;
      else return SPrint("[",Join(List(l,format)),"]");
      fi;
    else return Format(l,opt);
    fi;
  end;
  if IsBound(H.spets) then W:=Spets(H);H:=Hecke(H);
  else W:=Group(H);
  fi;
  para:=H.parameter;
  if not IsBound(opt.expand) then para:=AbsHeckeOps.CompactPara(para);fi;
  if ForAny([1..Length(H.parameter)],i->IsBound(H.rootParameter[i])) then
    root:=ShallowCopy(H.rootParameter);
    if Length(root)>0 and ForAll(root,x->x=root[1]) then root:=root[1];fi;
  fi;
  if IsBound(opt.GAP) then res:=SPrint("Hecke(",Format(W,opt));
  else
    if IsBound(W.operations.ReflectionName) then n:=ReflectionName(W,opt);
    else n:=String(W);
    fi;
    if IsBound(opt.TeX) then res:=SPrint("{\\mathcal H}_{",n,"}(");
    else res:=SPrint("Hecke(",n);
    fi;
  fi;
  if para<>1 or (IsBound(root) and root<>1) then 
     if not IsBound(opt.TeX) then PrintToString(res,",");fi;
     PrintToString(res,format(para));
  fi;
  if IsBound(root) and root<>1 then PrintToString(res,",",format(root));fi;
  PrintToString(res,")");
  return res;
end;

AbsHeckeOps.String:=Format;

AbsHeckeOps.Print:=function(H)Print(String(H));end;

AbsHeckeOps.Group:=H->H.reflectionGroup;

#############################################################################
##
#F  Hecke(<W>[,<parameters>, [<rootParameter>]] )  Make Hecke algebra of W
##
##  This  function  can  be  applied  to  any  (complex)  reflection group.
##  <parameters>  is a list of  lengths W.nbGeneratingReflections whose ith
##  element  is of length  W.OrdersGeneratingReflections[i]. The generating
##  elements  T(s) of the algebra are indexed by the generating reflections
##  of  W. If we set q:=<parameters>  and take s in W.generatingReflections
##  they satisfy:
##     (T(s)-q[s][1])..(T(s)-q[s][W.OrdersGeneratingReflections[s]])=0  
##  and  the same braid relations as the generators of W.
##
##  <parameters> can be entered in some abbreviated forms.
##  If it is omitted, it is interpreted  as the ring element 1 (then, by
##  further conventions  given below,  the Hecke algebra  constructed is
##  just the algebra of the group W).
##  If one Ring element q is given, it is interpreted as [q].
##  If <parameters>  is a  list [x]  of length 1,  it is  interpreted as
##  [x,...,x] (replicated W.nbGeneratingReflections times)
##  If parameters[s] is a single ring element q it is interpreted as
##   [q,E(e),..,E(e)^(e-1)] where e=W.OrdersGeneratingReflections[s]
##
##  Parameters  corresponding  to  conjugate  reflections  must  be  equal.
##  Parameters  can  be  left  unbound  for  some  reflections  in an orbit
##  provided they are given at least for the reflexion of smallest index in
##  that orbit.
##
##  For certain operations (e.g. character tables, Kazhdan-Lusztig bases of
##  Hecke  algebras of Coxeter groups) some roots of some of the parameters
##  are  necessary. In  the case  of Coxeter  groups only  square roots are
##  needed.  There  is  a  general  mechanism  to  specify a square root of
##  -q[s][0]/q[s][1]  in  the  Coxeter  case  by  giving  an extra argument
##  rootParameter.   If  not  given,  rootParameter  will  be  filled  with
##  GetRoot(-q[s][0]/q[s][1],2) when needed.
##
##  Example: Hecke(CoxeterGroup("B",3),[q^2,t^2],[-q,t]);
##    parameters interpreted as [[q^2,-1],[t^2,-1],[t^2,-1]]. 
##    Square root of q^2  taken to be -q and of t^2 to be t;
##
##  The code  works for any reflection  group, so in particular  for the
##  two most general  classes which are AbsCoxOps  et PermRootOps. Since
##  we don't have a  type which is union of both, we  stuff the code for
##  the  moment in  PermRootOps (where  no further  operations on  Hecke
##  algebras are defined).
##
CHEVIE.Cache.HeckeAlgebras:=true;
PermRootOps.Hecke:=function(arg)local W,H,i,j,e,res,o;
  W:=arg[1];
  H:=rec(reflectionGroup:=W,
 	 isDomain:=true, # for "Group" to work
	 operations:=AbsHeckeOps);
  if IsBound(arg[2]) then H.parameter:=ShallowCopy(arg[2]);
                     else H.parameter:=1;fi;
  if not IsList(H.parameter) then
    H.parameter:=List(W.generatingReflections,x->H.parameter);
  fi;
  if Length(H.parameter)=1 then
    H.parameter:=List(W.generatingReflections,x->H.parameter[1]);
  fi;
  if IsBound(arg[3]) then 
    H.rootParameter:=arg[3];
    if not IsList(H.rootParameter) then
      H.rootParameter:=List(W.generatingReflections,i->H.rootParameter);
    fi;
  else H.rootParameter:=[];
  fi;
  if IsFinite(W) then
    o:=List(W.reflections{W.generatingReflections},s->PositionClass(W,s));
    # better for CRG where there are too many orbits of roots
  else
    o:=W.orbitRepresentative{W.generatingReflections}; # for infinite coxeter
  fi;
  for i in W.generatingReflections do
    j:=Position(o,o[i]);
    if IsBound(H.parameter[i]) then
      e:=W.OrdersGeneratingReflections[i];
      if not IsList(H.parameter[i]) then
       H.parameter[i]:=Concatenation([H.parameter[i]],List([1..e-1],i->E(e)^i));
      elif e<>Length(H.parameter[i]) then
        Error("the ",i,"-th parameter list should be of length ",e);
      fi;
      if j<i then
        if H.parameter[j]<>H.parameter[i] then
	   Error("parameters should be equal for conjugate reflections ",
	         i," and ", j);
	fi;
      fi;
    elif j<i then H.parameter[i]:=H.parameter[j];
    else 
Error("parameters should be defined for at least one reflection in each orbit");
    fi;
    if not IsBound(H.rootParameter[i]) and IsBound(H.rootParameter[j]) then 
      H.rootParameter[i]:=H.rootParameter[j];
    fi;
    if IsBound(H.rootParameter[i]) and W.OrdersGeneratingReflections[i]=2 and 
      H.rootParameter[i]^2<>-Product(H.parameter[i])
    then Error("rootParameter[",i,"]^2 is not equal to ",-Product(H.parameter[i]));
    fi;
  od;
  if W.nbGeneratingReflections=0 then H.unit:=1;
  else H.unit:=Product(Concatenation(H.parameter))^0;
    # .unit introduced to try to deal with problems with ring embeddings in GAP3
  fi;
  return CHEVIE.GetCached(W,"HeckeAlgebras",H,
                                         x->[x.parameter,x.rootParameter]);
end; 

#############################################################################
##
#F  HeckeSubAlgebra( <H>, <roots> ) . . . . . . . . . . . . . . . . . . . or
#F  HeckeSubAlgebra( <H>, <subgroup> ) . . . . . . . .   Hecke Sub-Algebra
## 
##  Given a Hecke Algebra H and either a set of roots of Group(H) 
##  given as their index in the roots of W, or a reflection subgroup of W, 
##  return the Hecke sub-algebra generated by the T_s corresponding to 
##  these roots. The roots must be simple roots if the parameters are not
##  those of the group algebra of W.
##  As for Subgroup, a subalgebra of a subalgebra is given as a subalgebra
##  of the largest algebra.
##
AbsHeckeOps.HeckeSubAlgebra:=function(H,subW)local W,subroots,res,s,p;
  W:=Group(H);
  if IsList(subW) then subW:=ReflectionSubgroup(W,subW);fi;
  subroots:=W.rootRestriction{subW.rootInclusion{subW.generatingReflections}};
  s:=W.rootRestriction{W.orbitRepresentative{subroots}};
  p:=List(W.OrdersGeneratingReflections{s},i->List([1..i],j->E(i)^(j-1)));
  if H.parameter{s}=p then return Hecke(subW,p);
  elif IsSubset(W.generatingReflections,subroots) then
  res:=Hecke(subW,H.parameter{subroots},Sublist(H.rootParameter,subroots));
   res.parent:=H;
   return res;
  else Error("Generators of a sub-Hecke algebra should be simple reflections");
  fi;
end;

########################################################################
##
#F  CheckHeckeDefiningRelations( <H> , <t> )  . . . . . . . .  . . . . .
#F  check  the  defining  relations  of  a Hecke  algebra  for  a  given
#F  representation
##  
##  'CheckHeckeDefiningRelations'  returns true  or false,  according to
##  whether a  given set <t>  corresponding to  the standard
##  generators of the Hecke algebra <H> defines a representation.
## 
CheckHeckeDefiningRelations:=function(H,t)local W,r,id,res,e;
  W:=Group(H);res:=true;
  id:=t[1]^0;
  for r in W.generatingReflections do
    e:=Product(H.parameter[r],q->(t[r]-q*id));
    if e<>0*e then
      InfoChevie("#I  Error in ",Ordinal(r)," parameter relation\n");
      res:=false;
    fi;
  od;
  for r in BraidRelations(W) do
    e:=Product(t{W.rootRestriction{r[1]}})-Product(t{W.rootRestriction{r[2]}});
    if e<>0*e then
      InfoChevie("#I Error in relation ",r[1],"=",r[2],"\n");
      res:=false;
    fi;
  od;
  return res;
end;

IsHeckeAlgebra:=H->IsRec(H) and IsBound(H.parameter) 
                            and IsBound(H.reflectionGroup);

#############################################################################
##
## An exemple of CreateHeckeBasis: create the 'T' basis for
## AbsHeckeOps (i.e. all Hecke algebras for reflection groups...)
##
CreateHeckeBasis("T",rec(T:=x->x),     # method to convert to T
  AbsHeckeOps);
