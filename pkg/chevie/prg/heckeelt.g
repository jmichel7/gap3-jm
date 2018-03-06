#############################################################################
##
#A  heckeelt.g   CHEVIE library    Meinolf Geck, Andrew Mathas and Jean Michel 
##
#Y  Copyright (C) 1992 - 1999  Lehrstuhl D fur Mathematik, RWTH Aachen, and
#Y  Universite Paris VII.
##
##  This file contains  GAP functions for working with elements of an
##  Hecke algebra of a Coxeter group, in various bases.
##

HeckeEltOps:=OperationsRecord("HeckeEltOps");

IsHeckeElt:=h->IsRec(h) and IsBound(h.hecke) and IsBound(h.elm);

#############################################################################
##
#F  HeckeElt(H,basis,elm,coeff[,rec of extra args])  .... 
#F  or HeckeElt(h,elm,coeff)  .... 
#F  Internal service routine making an Hecke element out of its fields
##
##  The code has been gathered in one place for cleaner reference; also if we
##  want to modify something this is the only place needing a change.
##
HeckeElt:=function(arg)local res;
  if Length(arg)=3 then res:=rec(coeff:=arg[3],elm:=arg[2],
    basis:=arg[1].basis,hecke:=arg[1].hecke,
      operations:=arg[1].hecke.operations.(arg[1].basis));
    if IsBound(res.operations.extraFields) then
      Inherit(res,arg[1],res.operations.extraFields);fi;
  else
    res:=rec(coeff:=arg[4],elm:=arg[3],basis:=arg[2],
               operations:=arg[1].operations.(arg[2]),hecke:=arg[1]);
    if IsBound(arg[1].spets) then res.coset:=arg[1];res.hecke:=arg[1].hecke;fi;
  fi;
  if Length(arg)=5 then Inherit(res,arg[5]);fi;
  return res;
end;

HeckeEltOps.GetExtra:=function(h)local res;res:=rec();
  if IsBound(h.operations.extraFields) then 
    Inherit(res,h,h.operations.extraFields);fi;
  return res;
end;

##############################################################
##
#F  HeckeEltOps.Hecke(h)  returns Hecke algebra of h
##
HeckeEltOps.Hecke:=h->h.hecke;

#############################################################################
##
#F  HeckeEltOps.\+     . . . . . . . . . . . . . . Addition of Hecke elements
#F   If they are not of the same basis, both are converted to T.
##
HeckeEltOps.\+:=function(x,y)
  if x=0 then return y;
  elif y=0 then return x;fi; # to get Sum to work...
  if not IsIdentical(Hecke(x),Hecke(y)) then 
    Error("not elements of the same algebra");
  fi;
  if x.basis <> y.basis then 
    return Basis(Hecke(x),"T")(x)+Basis(Hecke(x),"T")(y);
  else
    x:=ShallowCopy(x);
    x.elm:=Concatenation(x.elm,y.elm);
    x.coeff:=Concatenation(x.coeff,y.coeff);
    CollectCoefficients(x);
    return x;
  fi;
end;

#############################################################################
##
#F  HeckeEltOps.\+     . . . . . . . . . . . . . Subtraction of Hecke elements
##
HeckeEltOps.\-:=function(h1,h2) return h1+(-1)*h2; end;

#############################################################################
##
#F  HeckeEltOps.\*(x,y) . . . . . . . . . . Multiplication of Hecke elements
#F   If x is a scalar, just multiply the coefficients of y by x.
#F   If they are not of the same basis return T(x)*T(y).
#F   If they are of T basis, use quadratic relations.
#F   If they are both of Basis X, return X(T(x)*T(y)).
##
HeckeEltOps.\*:=function(x,y)local H,W,ILD,multsi,res;
  if IsList(x) then return List(x,z->z*y);
  elif IsList(y) then return List(y,z->x*z);
  fi;
  if IsHeckeElt(y) then H:=Hecke(y);
    if IsHeckeElt(x)=false then # assume x is a scalar by which to multiply y
      if x=0*x then return HeckeElt(H,y.basis,[],[]);fi;
      y:=ShallowCopy(y);y.coeff:=y.coeff*(x*Hecke(y).unit);return y;
    fi;
  else # assume y is a scalar by which to multiply x
    if y=0*y then return HeckeElt(Hecke(x),x.basis,[],[]);fi;
    x:=ShallowCopy(x);x.coeff:=x.coeff*(y*Hecke(x).unit);return x;
  fi;
  if not IsIdentical(H,Hecke(x)) then 
    Error("not elements of the same algebra");
  fi;
  if x.basis<>y.basis then 
    if x.operations.\* <> y.operations.\* then return x.operations.\*(x,y);
    else return Basis(H,"T")(x)*Basis(H,"T")(y);
    fi;
  elif x.basis="T" then 
    W:=Group(H);
    if not IsCoxeterGroup(W) then
      Error("cannot multiply cyclotomic Hecke algebras elements");
      return false;
    fi;
    ILD:=W.operations.IsLeftDescending;
  multsi:=function(i,h)local q,res,up,down,s,j; #  T_{s_i}*h
    q:=H.parameter[i];s:=W.reflections[i]; res:=HeckeElt(H,h.basis,[],[]);
    up:=[];down:=[];
    for j in [1..Length(h.elm)] do 
      if ILD(W,h.elm[j],i) then Add(down,j);else Add(up,j);fi;
    od;
    Append(res.elm,h.elm{down}); Append(res.coeff,(q[1]+q[2])*h.coeff{down});
    Append(res.elm,s*h.elm{down}); Append(res.coeff,-q[1]*q[2]*h.coeff{down});
    Append(res.elm,s*h.elm{up}); Append(res.coeff,h.coeff{up});
    return res;
  end;
    if Length(x.elm)=0 then return HeckeElt(H,x.basis,[],[]);fi;
    res:=Sum([1..Length(x.elm)],function(i)local res,w,j; res:=x.coeff[i]*y;
      w:=W.rootRestriction{Reversed(CoxeterWord(W,x.elm[i]))};
      for j in [1..Length(w)] do res:=multsi(w[j],res);
        if j mod 2=0 then CollectCoefficients(res);fi; # 2 experimentally best
      od;
      return res;
    end);
    if Length(x.elm)=1 then CollectCoefficients(res);fi;
    return res;
  elif IsBound(x.operations.extraFields) then
    return Basis(H,x.basis)(Basis(H,"T")(x)*Basis(H,"T")(y),
      HeckeEltOps.GetExtra(x));
  else return Basis(H,x.basis)(Basis(H,"T")(x)*Basis(H,"T")(y));
  fi;
end;

#############################################################################
##
#F  HeckeEltOps.\^(x,y) . . . . . . . . . . Exponentiation for Hecke elements
#F   If y is an integer, return x^y using multiplication. y is allowed to be
#F    negative only for bases which have an .inverse method.
#F   Otherwise return the conjugate y^-1*x*y
##
HeckeEltOps.\^:=function(h,n) local i,p,H,W,T,q,e,res;
  H:=Hecke(h);W:=Group(H);
  if not IsInt(n) then return n^-1*h*n;fi; # assume n is a Hecke element
  if n<0 then
    if n=-1 then
      if IsBound(h.operations.inverse) then return h.operations.inverse(h);
      else Error(h," has no method for inverse");
      fi;
    else return (h^-1)^(-n);
    fi;
  fi;
  if n=0 then if IsBound(h.operations.extraFields) then
   return Basis(H,h.basis)(Basis(H,"T")(),HeckeEltOps.GetExtra(h));
  else return Basis(H,h.basis)(Basis(H,"T")());fi;
  fi; 
  # note that the above possibly expensive expression is called only for h^0
  p:=false;
  while true do
   if n mod 2 <> 0 then 
     if p=false then p:=h;
     else p:=p*h;
     fi;
   fi;
   n:=QuoInt(n,2);
   if n=0 then return p;fi;
   h:=h*h;
  od;
end;

#############################################################################
##
#F  HeckeEltOps.\/(x,y) . . . . . . . . . . Division for Hecke elements
##
HeckeEltOps.\/:=function(a,b) return a*b^-1;end;

#############################################################################
##
#F AlphaInvolution(h)    
## The involution on Hecke Elements defined by T_w->T_{w^{-1}}
## (and same in other bases)
HeckeEltOps.AlphaInvolution:=h->HeckeElt(h,List(h.elm,x->x^-1),h.coeff);

CHEVIE.PrintHecke:=rec();

#############################################################################
##
#F  HeckeEltOps.Format(h) . . . . . . . . . . Format method for Hecke elements
##
HeckeEltOps.Format:=function(h,option)local res,b,W,opt;
  opt:=ShallowCopy(CHEVIE.PrintHecke);
  Inherit(opt,option);
  if h.elm=[] then return "0";fi;
  if IsBound(h.operations.FormatBasis) then b:=h.operations.FormatBasis(h);
  else b:=h.basis; fi;
  W:=Group(Hecke(h));
  if IsCoxeterGroup(W) then h:=[List(h.elm,x->CoxeterWord(W,x)),h.coeff];
  else                      h:=[h.elm,h.coeff];
  fi;
  h:=TransposedMat(h);
  SortBy(h,x->[Length(x[1]),x[1]]);
  res:=Concatenation(List(h,function(i)local s;
    s:=FormatCoefficient(i[2],SPrint(b,"(",Join(i[1]),")"),opt);
    if s[1]<>'-' then return SPrint("+",s);else return s;fi;end));
  if res[1]='+' then res:=res{[2..Length(res)]};fi;
  return res;
end;

HeckeEltOps.String:=h->Format(h,rec());

#############################################################################
##
#F  HeckeEltOps.Print(h) . . . . . . . . . . Print method for Hecke elements
##  .Print calls .String as is general policy in CHEVIE
##
HeckeEltOps.Print:=function(h)Print(String(h));end;

#############################################################################
##
#F  HeckeEltOps.Coefficient(h,e) . . . . . . . . . . . Coefficient of e in h
#F     e is an element of Group(h) or a word in the generators of Group(h)
##  Uses that Coefficient is generic.
##
HeckeEltOps.Coefficient:= function(T,elm)local p,W;
  W:=Group(Hecke(T));
  if IsCoxeterGroup(W) and IsWordFor(W,elm) then elm:=EltWord(W,elm);fi;
  p:=PositionSet(T.elm,elm);
  if p=false then return 0*Hecke(T).unit;else return T.coeff[p];fi;
end;

#############################################################################
##
#F  HeckeClassPolynomials( <h> )  . . . . . . class polynomials of <h>
##
##  returns  the  class  polynomials  of  the  element  <h> with respect to
##  representatives  of minimal length in  the (F-)conjugacy classes of the
##  Coxeter  group (or coset). <h> is an element of <H> given in any basis.

HeckeEltOps.HeckeClassPolynomials:=function(h)
  local H,minl,min,l,maxl,new,i,o,orb,W,WF,p;
  H:=Hecke(h);
  if IsBound(h.coset) then WF:=Spets(h.coset); W:=Group(WF);
  else W:=Group(H);WF:=W;fi;
  minl:=List(ChevieClassInfo(WF).classtext,Length);
  h:=Basis(H,"T")(h);
# Since  vF is not of minimal length in its class there exists wF conjugate
# by   cyclic  shift  to  vF  and  a  generating  reflection  s  such  that
# l(swFs)=l(vF)-2. Return T_sws.T_s^2
  orb:=function(orbit)local q,w,sw,s,sws,ILD;
    ILD:=W.operations.IsLeftDescending; # avoid dispatching overhead
    for w in orbit do
      for s in W.generatingReflections do
        if ILD(W,w,s) then
          sw:=W.reflections[s]*w;
          sws:=sw*W.reflections[s];
          if ILD(W,sw^-1,s) then q:=H.parameter[s];
            return rec(elm:=[sws,sw],coeff:=[-q[1]*q[2],q[1]+q[2]]);
          elif not sws in orbit then Add(orbit,sws);
          fi;
        fi; 
      od;
    od;
    Error("Geck-Kim-Pfeiffer theory");
  end;

  min:=minl*0*H.unit;
  while Length(h.elm)>0 do
    new:=rec(elm:=[],coeff:=[]);
    l:=List(h.elm,x->W.operations.CoxeterLength(W,x));
    maxl:=Maximum(l); 
    for i in [1..Length(h.elm)] do
      if l[i]<maxl then Add(new.elm,h.elm[i]);Add(new.coeff,h.coeff[i]);
      else
        p:=PositionClass(WF,h.elm[i]);
        if minl[p]=maxl then min[p]:=min[p]+h.coeff[i];
        else o:=orb([h.elm[i]]);
          Append(new.elm,o.elm);Append(new.coeff,o.coeff*h.coeff[i]);
        fi;
      fi;
    od;
    CollectCoefficients(new);h:=new;
  od;
  return min;
end;

#############################################################################
##
#F  HeckeCharValues( <h> [,<irreds>] )  . . character values of <irreds> on <h>
##
##  h  is an element of an Iwahori-Hecke algebra H (belonging to any basis)
##  and  <irreds> is a set of  irreducible characters of H. HeckeCharValues
##  returns the values of irreds on T, using HeckeClassPolynomials.
##   If absent, <irreds> is taken to be all irreducible characters of <H>.

HeckeEltOps.HeckeCharValues:=function(arg)local irrs,H,cox;
  if IsBound(arg[1].coset) then H:=arg[1].coset;cox:=IsCoxeterCoset(Spets(H));
  else H:=Hecke(arg[1]);cox:=IsCoxeterGroup(Group(H));
  fi;
  if cox then
    if IsBound(arg[2]) then irrs:=arg[2]*H.unit^0;
    else irrs:=CharTable(H).irreducibles;fi;
    return irrs*HeckeClassPolynomials(arg[1]);
  else return Sum(Zip(arg[1].elm,arg[1].coeff,function(elm,coeff)
    return HeckeCharValues(H,elm)*coeff;end));
  fi;
end;

#############################################################################
##
#F  Representation( <h> , <n>)  . . value in nth representation
##                                  or in representation given by list l
##
##  h  is an element of an Iwahori-Hecke algebra H (belonging to any basis)
##  returns the values of repr. n on T.

HeckeEltOps.Representation:=function(h,n)local H,W,r;
  if IsBound(h.coset) then H:=h.coset; else H:=Hecke(h); fi;
  W:=Group(H);h:=Basis(H,"T")(h);
  if IsInt(n) then r:=Representation(H,n);else r:=n;fi;
  if Length(h.coeff)=0 then return 0*r[1];fi;
  return Sum(Zip(h.coeff,h.elm,function(c,w)
    if w=() then return r[1]^0*c;else return Product(r{CoxeterWord(W,w)})*c;
    fi;end));
end;

#############################################################################
##
#F  HeckeEltOps.Frobenius(W,h,i) . . . . . . . . . . . Frobenius(W)(h[,i])
#F     W is a CoxeterCoset, and h an element of an HeckeCoset of W or of
#F     an hecke algebra of the group of which W is a coset.
##
HeckeEltOps.Frobenius:=function(W,x,i)
  x:=ShallowCopy(x);x.elm:=List(x.elm,y->Frobenius(W)(y,i));
  return x;
end;

#############################################################################
##
#F CreateHeckeBasis(basis,basisops[,algebraops]) . . . . . . . Create a new
#F   basis of Hecke algebras which have algebraops as operations.
##
## basis, a string, is  the name of the new basis  and basisops should be
## a  record  containing the  operations  in  HeckeEltOps which  will  be
## overridden for that particular basis. basisops should contain at least
## methods .T  and .(basis)  to convert  to and from  the T-basis.  If it
## contains a  method init(H) this  is called when  using Basis(H,basis).
## Most  other methods  of HeckeEltOps  can  be overriden  if needed,  in
## particular MakeBasisElt called  by 'Basis' is overriden  by the Murphy
## basis code (see the package 'Specht').
##
CreateHeckeBasis:=function(arg)local basis,basisops,algebraops;
# default 3rd argument for compatibility with previous versions of CHEVIE
  basis:=arg[1];basisops:=arg[2];
  if Length(arg)=2 then algebraops:=HeckeAlgebraOps;
  else algebraops:=arg[3];
  fi;
# Print("Create(",arg[1],",",RecFields(arg[2]),",",algebraops,")\n");
  if not IsBound(algebraops.Basis) then # first basis ever created for algebra
    algebraops.Basis:=function(H,basis)local basisops;
    # we use that Basis is generic
     if not IsBound(H.operations.(basis)) or 
        not IsBound(H.operations.(basis).T)  # test it is not a module basis
     then 
        Error("basis ",basis," unknown");
     fi;
     basisops:=H.operations.(basis);
     if IsBound(basisops.init) then basisops.init(H);fi;
     return function(arg) return basisops.MakeBasisElt(H,basis,arg);end;
    end;
  fi;

  if not IsRec(basisops) or not IsBound(basisops.T) 
                         or not IsBound(basisops.(basis)) then 
    Error("The operations record must contain methods .T and .",basis);
  fi;
  algebraops.(basis):=OperationsRecord(Concatenation("Hecke",basis,"Ops"),
                                       HeckeEltOps);
  Inherit(algebraops.(basis),basisops);
end;

#############################################################################
##
#F HeckeEltOps.MakeBasisElt(H,basis,arg)   . . after setting X:=Basis(H,"X")
#F   the hecke-element of the 'X' basis making function X will contain
#F    ``function(arg) return H.operations.X.MakeBasisElt(H,basis,arg);end;''
#F   If the default version given
#F   below of MakeBasisElt is used it will accept the forms (here W=Group(H)):
## 
##1.  X(<group element g>[,extra]) Interpreted as X([g],[1])
##2.  X([s_1,...,s_n][,extra])     where s_i is the name of generator S_i of W:
##        Interpreted as X(Product_i S_i) if s_1..s_n reduced. Otherwise
##        basis must be T and returns Product_i T(S_i) otherwise. As a 
##        special convention, a negative s_i represents the inverse of the
##        corresponding generator.
##3.  X(s_1,...,s_n[,extra])       Interpreted as X([s_1,...,s_n])
##4.  X([elts],[coeffs][,extra])   Makes the element of basis "X" with
##                           elements [elts] and coefficients [coeffs]
##5.  X(h[,extra])  In this form h is a hecke element and the function
##                  tries to convert h to basis "X". It first looks if h has a
##                  method .("X"), and if not converts h to T and then to X.
HeckeEltOps.MakeBasisElt:=function(H,basis,arg)
  local w,res,s,h,basisops,W,isRefl,isExtra,extra,id;
  basisops:=H.operations.(basis);
  isExtra:=x->IsBound(basisops.extraFields) and IsRec(x)
    and Set(RecFields(x))=Set(basisops.extraFields);
  W:=Group(H);
  if IsCoxeterGroup(W) then id:=W.identity;else id:=[];fi;
  if arg=[] then return HeckeElt(H,basis,[id],[H.unit]); # case 2.
  elif Length(arg)=1 and isExtra(arg[1]) then
    return HeckeElt(H,basis,[id],[H.unit],arg[1]); # case 2.
  fi;
  if IsHeckeElt(arg[1])then # case 5. 
    h:=arg[1];
    if h.basis=basis then return h;
    elif IsBound(h.operations.(basis)) then 
      return ApplyFunc(h.operations.(basis),arg);
    elif Length(arg)=1 then return basisops.(basis)(h.operations.T(h));
    elif Length(arg)=2 and isExtra(arg[2]) then
      return basisops.(basis)(h.operations.T(h),arg[2]);
    else Error("expecting ",basis,"(<Hecke elt>[,<extra>])\n");
    fi;
  fi;
  isRefl:=x->x in W.reflectionsLabels or(IsInt(x)and -x in W.reflectionsLabels);
  if (Length(arg)>=2 and IsList(arg[1]) and not isRefl(arg[1]) and 
     IsList(arg[2]) and not isRefl(arg[2])) and 
     (Length(arg)=2 or (Length(arg)=3 and isExtra(arg[3]))) then # case 4.
    res:=ApplyFunc(HeckeElt,Concatenation([H,basis],List(arg,ShallowCopy)));
    CollectCoefficients(res);
    return res;
  fi;
  if not isRefl(arg[1]) and (Length(arg)=1 or Length(arg)=2 and isExtra(arg[2]))
  then
    h:=arg[1];
    if not(IsList(h) and ForAll(h,isRefl)) then
      return ApplyFunc(HeckeElt,Concatenation([H,basis,[h],[H.unit]],
         arg{[2..Length(arg)]}));
      # assume h is a group element --- case 1.
    fi;
    if Length(arg)=2 then extra:=arg[2];fi;
  else h:=arg; # case 3
  fi;
  # now case 3. Check [s_1,...,s_n] is reduced -- else assume basis is "T"
  if Length(h)>0 and isExtra(h[Length(h)]) then 
    extra:=h[Length(h)];h:=h{[1..Length(h)-1]};fi;
  if ForAll(h,x->not IsInt(x) or x>0) then
    if IsCoxeterGroup(W) then
      w:=EltWord(W,h);
      if CoxeterLength(W,w)=Length(h) then 
        if IsBound(extra) then return HeckeElt(H,basis,[w],[H.unit],extra);
        else return HeckeElt(H,basis,[w],[H.unit]);
        fi;
      fi;
    else return HeckeElt(H,basis,[h],[H.unit]);
    fi;
  fi;
  if basis<>"T" then
    Error("Construction with non-reduced words implemented only for basis T");
  fi;
  res:=HeckeElt(H,basis,[id],[H.unit]);
  for s in h do 
    if s>0 then res:=res*HeckeElt(H,basis,
        [W.reflections[W.operations.ReflectionFromName(W,s)]],[H.unit]);
    else res:=res/HeckeElt(H,basis,
        [W.reflections[W.operations.ReflectionFromName(W,-s)]],[H.unit]);
    fi;
  od;
  return res;
end;

HeckeEltOps.Specialization:=function(t,H2,f)
  return Basis(H2,t.basis)(t.elm,List(t.coeff,f));
end;
