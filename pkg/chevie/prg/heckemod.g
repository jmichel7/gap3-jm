#############################################################################
##
#A  heckemod.g             CHEVIE library       Jean Michel, Raphael Rouquier
##
#Y  Copyright (C) 1997 - 1999  Universite Paris VII.
##
##  This file contains  GAP functions for working with elements of module of
##  an Hecke algebra of a Coxeter group, in various bases.
##

HeckeModuleOps:=OperationsRecord("HeckeModuleOps",HeckeEltOps);
   # operations for Hecke modules

#############################################################################
##
#F  HeckeModuleOps.MakeRec(H,basis,reflections,chi,elm,coeff) Internal service
#F    routine: makes an Hecke module element out of its fields
##
##  The code has been gathered in one place for cleaner reference; also if we
##  want to modify something this is the only place needing a change.
##
HeckeModuleOps.MakeRec:=function(H,basis,reflections,chi,elm,coeff)
  return rec(coeff:=coeff, elm:=elm, basis:=basis,
	     operations:=H.operations.(basis),hecke:=H,
	     reflections:=reflections,chi:=chi);
end;

##     ModuleBasis(<H>,<basis> [, <reflections> [,<chi>]])
## or  ModuleBasis(<module element>)
##
##  if not given by default reflections:=[1..W.nbGeneratingReflections-1]
##                      and chi:=-1 (meaning sign character)
##
ModuleBasis:=function(arg)local H,basis,reflections,chi,W;
  if IsHeckeAlgebra(arg[1]) then
     H:=arg[1];basis:=arg[2];
     W:=Group(H);
     if Length(arg)=2 then reflections:=[1..W.nbGeneratingReflections-1];
     else reflections:=arg[3];
     fi;
     if Length(arg)<=3 then chi:=List(reflections,x->-1);
     else chi:=arg[4];
     fi;
     if not IsList(chi) then chi:=List(reflections,x->chi);fi;
  else
     H:=Hecke(arg[1]);
     reflections:=arg[1].reflections;
     chi:=arg[1].chi;
     basis:=arg[1].basis;
  fi;
  if not IsBound(H.operations.(basis)) or not IsBound(H.operations.(basis).MT)
  then Error("basis ",basis," unknown");fi;
  if IsBound(H.operations.(basis).init) then 
     H.operations.(basis).init(H,reflections,chi); fi;
  return function(arg) return 
   ApplyFunc(H.operations.MakeModuleElt,
     Concatenation([H,basis,reflections,chi],[arg]));end;
end;

HeckeModuleOps.\+:=function(x,y)
  if x=0 then return y;
  elif y=0 then return x;fi; # to get Sum to work...
  if not IsIdentical(Hecke(x),Hecke(y)) then 
    Error("elements do not have the same Hecke algebra");
  fi;
  if x.chi<>y.chi or x.reflections<>y.reflections then
    Error("not elements of the same module");
  fi;
  if x.basis <> y.basis then 
    return ModuleBasis(Hecke(x),"MT",x.reflections,x.chi)(x)+
           ModuleBasis(Hecke(x),"MT",x.reflections,x.chi)(y);
  else
    x:=ShallowCopy(x);
    x.elm:=Concatenation(x.elm,y.elm);
    x.coeff:=Concatenation(x.coeff,y.coeff);
    CollectCoefficients(x);return x;
  fi;
end;

HeckeModuleOps.\*:=function(x,y)
  local q,s,W,res,H,M,MT,T,i,j,tmp,tempcoeff,tempelm,basis,sj,elm;
  if not IsRec(x) or not IsBound(x.hecke) or not IsBound(x.elm) then 
  # assume x is a scalar by which to multiply y
    res:=ShallowCopy(y);res.coeff:=res.coeff*x;
    return res;
  fi;
  H:=Hecke(y);
  W:=Group(H);
  T:=Basis(H,"T");
  basis:=y.basis;
  M:=ModuleBasis(y);
  y:=y.operations.MT(y);
  MT:=ModuleBasis(y);
  res:=MT([],[]);
  x:=T(x);
  for i in [1..Length(x.coeff)] do
    tmp:=y;
    elm:=x.elm[i]^-1;
    while elm<>W.identity do
      s:=FirstLeftDescending(W,elm);elm:=W.reflections[s]*elm;
      tempcoeff:=[];tempelm:=[];
      for j in [1..Length(tmp.coeff)] do
	sj:=W.reflections[s]*tmp.elm[j];
	if ForAny(RightDescentSet(W,sj),x->x in y.reflections) then
	    Add(tempcoeff,
         y.chi[Position(W.reflections{W.rootRestriction{y.reflections}},
	      tmp.elm[j]^-1*sj)]*tmp.coeff[j]);
	  Add(tempelm,tmp.elm[j]);
	elif IsLeftDescending(W,tmp.elm[j],s) then
	  q:=H.parameter[s];
	  Add(tempcoeff,-q[1]*q[2]*tmp.coeff[j]);Add(tempelm,sj);
	  Add(tempcoeff,(q[1]+q[2])*tmp.coeff[j]);Add(tempelm,tmp.elm[j]);
	else
	  Add(tempcoeff,tmp.coeff[j]);Add(tempelm,sj);
	fi;
      od;
      tmp:=MT(tempelm,tempcoeff);
      CollectCoefficients(tmp);
    od;
    res:=res+x.coeff[i]*tmp;
  od;
  return M(res);
end;

HeckeModuleOps.Specialization:=function(t,H2,f)local res;
  res:=ModuleBasis(H2,t.basis,t.reflections,t.chi)(t.elm,List(t.coeff,f));
  CollectCoefficients(res);return res;
end;

## makes an element of the module Wa/WI over a polynomial ring 
## returns a record consisting of two lists: elm=reduced expressions
##                                           coeff=polynomials
##
#############################################################################
##
#F algebraops.(basis).MakeModuleElt(H,arg)  after setting X:=ModuleBasis(H,"X")
#F   the hecke-module-element-of-the-'X'-basis-making-function X will contain
#F   the code found in H.operations.X.MakeBasisElt. If the default version given
#F   below of MakeBasisElt is used it will accept the forms (here W=Group(H)):
## 
##1.  X(<group element g>) Interpreted as X([g],[1])
##2.  X([s_1,...,s_n])     where s_i is the name of generator S_i of W:
##                         Interpreted as Product_i X(S_i)
##3.  X(s_1,...,s_n)       Interpreted as X([s_1,...,s_n])
##4.  X([elts],[coeffs])   Makes the element of basis "X" with elements [elts]
##                         and coefficients [coeffs]
##5.  X(h)              In this form h is a module element and the function
##                  tries to convert h to basis "X". It first looks if h has a
##                  method .("X"), and if not converts h to MT and then to X.
CoxeterHeckeAlgebraOps.MakeModuleElt:=function(H,basis,reflections,chi,arg)
   local w,res,s,h,W,ops,basisops;
  ops:=H.operations;basisops:=H.operations.(basis);
  if Length(arg)=1 and IsRec(arg[1]) and IsBound(arg[1].hecke) and 
     IsBound(arg[1].elm) then  # case 5. assume arg is some hecke element h
    h:=arg[1];
    if h.basis=basis then return h;
    elif IsBound(h.operations.(basis)) then return h.operations.(basis)(h);
    else return basisops.(basis)(h.operations.MT(h));
    fi;
  fi;
  W:=Group(H);
  if arg=[] then # case 2.
    return HeckeModuleOps.MakeRec(H,basis,reflections,chi,[W.identity],[H.unit]);
  fi;
  if IsWordFor(W,arg) then  # case 3.
    h:=arg;arg:=[arg];
  fi;
  h:=arg[1];
  if IsBound(arg[2]) then # case 4.
    return HeckeModuleOps.MakeRec(H,basis,reflections,chi,h,arg[2]);
  elif not IsWordFor(W,h)
  then  # assume h is a group element --- case 1.
    if ForAny(RightDescentSet(W,h),x->x in reflections) then
      Error("argument should be ",reflections,"-reduced");
    fi;
    return HeckeModuleOps.MakeRec(H,basis,reflections,chi,[h],[H.unit]);
  fi;
  # now case 2. Check [s_1,...,s_n] is reduced -- else assume basis is "T"
  w:=EltWord(W,h);
  if CoxeterLength(W,w)=Length(h) then 
    if ForAny(RightDescentSet(W,w),x->x in reflections) then
      Error("argument should be ",reflections,"-reduced");
    fi;
    return HeckeModuleOps.MakeRec(H,basis,reflections,chi,[w],[H.unit]);
  else
   Error("reduced words needed to construct module elements");
  fi;
end;

CreateHeckeModuleBasis:=function(basis,basisops,algebraops)
  if not IsRec(basisops) or not IsBound(basisops.MT) or not IsBound(basisops.(basis))
  then Error("The operations record must contain methods .MT and .",basis);
  fi;
  algebraops.(basis):=OperationsRecord(
    Concatenation("HeckeModule",basis,"Ops"),HeckeModuleOps);
  Inherit(algebraops.(basis),basisops);
end;

CreateHeckeModuleBasis("MT",rec(MT:=x->x),CoxeterHeckeAlgebraOps);

CreateHeckeModuleBasis("MC'",rec(
   init:=function(H,reflections,chi)
     if not Set(chi)=[-1] then 
        Error("MC' only implemented for the sign character");fi;
     H.operations.InitKL(H,"MC' basis");
   end,
   MT:=function(x)local H;
     H:=Hecke(x);
     return Sum([1..Length(x.elm)],i->x.coeff[i]*
	Basis(H,"C'")(x.elm[i])*ModuleBasis(H,"MT",x.reflections,x.chi)());
     end,
   ("MC'"):=function(x)local H,t,c,W;
     H:=Hecke(x);W:=Group(H);
     t:=Basis(H,"C'")(Basis(H,"T")(x.elm,x.coeff));
     c:=List(t.elm,y->not ForAny(RightDescentSet(W,y),z->z in x.reflections));
     return ModuleBasis(H,"MC'",x.reflections,x.chi)
                            (ListBlist(t.elm,c),ListBlist(t.coeff,c));
  end),CoxeterHeckeAlgebraOps);
