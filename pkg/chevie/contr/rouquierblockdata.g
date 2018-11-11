###########################################################################
# Contribution to the Chevie Package
#
# rouquierblockdata.g
#
# (C) July 2015 --- Maria Chlouveraki and Jean Michel
#
# The  main  function  of  this  file  is  'RouquierBlocks'  which computes
# Rouquier  blocks  of  1-cyclotomic  Hecke  algebras of complex reflection
# groups.
#
# For each type of irreducible algebras, 'RouquierBlockData' is computed
# which could be stored and make the computation very fast...
###########################################################################
#
# A  1-cyclotomic Hecke algebra H is  an Hecke algebra whose j_th parameter
# for  the i-th reflection (of order e) is of the form E(e)^j x^m_{i,j} for
# some  rational numbers m_{i,j};  thus such an  algebra specializes to the
# group algebra for x->1. Here in CHEVIE x must be Mvp("x").
# 
# A  tool  to  determine  the  Rouquier  blocks  of  H  are  the "essential
# hyperplanes"  which are linear  forms in the  $m_{i,j}$ determined by the
# Schur elements of the generic algebra associated to H. For each essential
# hyperplane  h  there  is  an  associated  1-cyclotomic  algebra A_h whose
# $m_{i,j}$  annihilate h  and no  other essential  hyperplane. Let us call
# h-blocks  the Rouquier blocks of  A_h. Then the Rouquier  blocks of H are
# the  Lcm of the h-blocks  for h annihilating the  $m_{i,j}$ of H. In case
# the  $m_{i,j}$ annihilate  no hyperplane  we get  the 0-blocks  which are
# common to all 1-cyclotomic algebras with that property.

# PValuation(c,p) returns the p-adic valuation v of c where p is a prime
#  number and c is:
#  -a cyclotomic number: v is the largest integer (positive or negative)
#   such that cp^-v is an algebraic integer.
#  -a list: v is the minimum of the valuations of its entries.
#  -an Mvp: v is the minimum of the valuations of the coefficiants.
#
PValuation:=function(c,p)local e;
  if c=0*c then return 0;Error("first argument should not be zero");
  elif IsInt(c) then e:=0; while c mod p=0 do e:=e+1;c:=c/p;od; return e;
  elif IsRat(c) then 
    return PValuation(Numerator(c),p)-PValuation(Denominator(c),p);
  elif IsCyc(c) then e:=Conjugates(c);
    return PValuation(Product(e),p)/Length(e);
  elif IsList(c) then return Minimum(List(c,x->PValuation(x,p)));
  elif IsMvp(c) then return PValuation(c.coeff,p);
  else Error("PValuation not implemented for ",c);
  fi;
end;

# Given  a  complex  reflection  group  W  , the function RouquierBlockData
# returns a list of [essential hyperplane h, corresponding h-blocks]
# h  is represented  as a  list of  integers of  same length as the list of
#    parameters for the Hecke algebra of W.
# h-blocks is a partition of [1..NrConjugacyClasses(W)]
#
# The  first entry in the result list has h=[0,...,0] and the corresponding
# h-blocks are the 0-blocks.
#
RouquierBlockData:=function(W)
  local hplanes,bl0,sch,NRPARA,Generic1SchurElements;
  NRPARA:=5; # how many random algebras A_h to consider

# Generic1SchurElements returns the Schur elements for the Hecke algebra of
# W  with parameters  \zeta_e^i x_{j,i},  a variant  of the generic algebra
# which specializes to the group algebra for x_{i,j}->1.
#    The  result is  a record  with fields  .coeff,.mon  representing the 
# leading monomial, and a field .vcyc which is a list of records with fields
# .coeff, .mon representing a monomial m and a field .pol holding a cycpol
# such that the record represents pol(m)
  Generic1SchurElements:=function(W)local c,o,v,H,vnames,vars;
    c:=Set(W.orbitRepresentative);
    o:=List(W.EigenvaluesGeneratingReflections{c},x->1/x);
    vars:="xyz";
    v:=[];v{c}:=List([1..Length(o)],i->List([1..o[i]]-1,
		     j->Mvp(SPrint([vars[i]],j))*E(o[i])^j));
    H:=ApplyFunc(Hecke,[W,v]);
    vnames:=Variables(H.parameter);
    return List(FactorizedSchurElements(H),function(s)local res,montovec,f;
      montovec:=function(mon)local res;
	res:=List(vnames,x->0);
	res{List(mon.elm[1].elm,x->Position(vnames,x))}:=mon.elm[1].coeff;
	return res;
      end;
      if IsInt(s.factor) then f:=Mvp(s.factor);else f:=s.factor;fi;
      res:=rec(vars:=List(vnames,Mvp),
	vcyc:=List(s.vcyc,p->rec(pol:=p.pol,mon:=montovec(p.monomial),
	  coeff:=p.monomial.coeff[1])),
        coeff:=f.coeff[1],
	mon:=montovec(f));
  #   res.operations:=Generic1SchurElementOps;
      SortBy(res.vcyc,x->[x.pol,x.mon,x.coeff]);
      return res;
    end);
  end;
  sch:=Generic1SchurElements(W);
  hplanes:=Concatenation(List(sch,x->List(x.vcyc,y->y.mon)));
  hplanes:=List(hplanes,v->v*Lcm(List(v,Denominator)));
  hplanes:=Set(List(hplanes,v->v/Gcd(List(v,Numerator))));
  Sort(hplanes); 
  hplanes:=Concatenation([hplanes[1]*0],hplanes); # [0,..,0]+essential hplanes
  InfoChevie("#I ",Length(hplanes)," hplanes\n");

  return List(hplanes, # for each hplane h return [h,Rouquier blocks of A_h]
    function(h)local aA,res,para,findpara,c,m,p,hh; 
#    Print("\nhplane=",h);
    hh:=Filtered(hplanes,k->k<>0*k and k<>h);
    m:=NullspaceIntMat(TransposedMat([h]));
    para:=[];
    while Length(para)<NRPARA do
      if Length(hh)=0 and h=[1,-1] then Add(para,[1,1]);
      else p:=List([1..Length(m)],i->Random([-2*Length(m)..2*Length(m)]))*m;
	if not 0 in hh*p then Add(para,p/Gcd(p));fi;
      fi;
    od;
    SortBy(para,x->x*x); # increasing "complexity"
# para holds NRPARA random lists of m_{i,j} defining each a possible A_h

    aA:=ApplyFunc(GcdPartitions,List(para,p->CollectBy([1..Length(sch)],
       i->2*sch[i].mon*p+Sum(sch[i].vcyc,x->x.mon*p*Degree(x.pol)))));
# aA holds the (a+A)-blocks common to all NRPARA possible A_h

    para:=para[1]; # choose now the first A_h with "simplest" parameters
#    Print(Stime()," para=",para);
    c:=List(sch,s->s.coeff*Product(s.vcyc,function(r)
      if para*r.mon=0 then return Value(r.pol,r.coeff);else return 1;fi;end));
# c holds the leading coefficients of the Schur elements of A_h

# compute now the p-blocks of A_h for all primes p dividing |W|
    res:=List(Set(FactorsInt(Size(W))),function(p)local bl,vp,j,x,i,cut; 
#    Print(" p=",p);
      bl:=GcdPartitions(PBlocks(W,p),aA);
# here bl holds the coarsest partition 
# - finer than the p-blocks of W
# - finer than the (a+A)-blocks
      if h=[1,-1] then return bl;fi; # for 1-parameter algebras bl is h-blocks
      vp:=List(c,x->PValuation(x,p));
      i:=1;
      while i<=Length(bl) do 
# we examine each "pseudo-block" obtained and apply various rules to refine it
	if Length(bl[i])=1 then 
	  if vp[bl[i][1]]<>0 then 
	     Error("Schur elt of v_",p,"=",vp[bl[i][1]]," alone in block");
	  else i:=i+1;
	  fi;
	else 
# a Schur element s is alone in its block iff v_p[c[s]]=0
	  for j in bl[i] do
	    if vp[j]=0 and Length(bl[i])>1 then 
	      Add(bl,[j]);bl[i]:=Difference(bl[i],[j]);
	    fi;
	  od;
# a pseudo-block of size <=3 is a block
	  if Length(bl[i])>3 then 
	    if h<>h*0 then
# if B is a pseudo-block of size>3 and C\subset B is a 0-block
# and for each s in C we have vp(c(s))=v_p(c_0(s)) then C is a block
	      for j in Filtered(bl0,x->IsSubset(bl[i],x) and x<>bl[i]) do
		if vp{j}=List(sch{j},s->PValuation(s.coeff,p)) then
		  bl[i]:=Difference(bl[i],j);Add(bl,j);
		fi;
	      od;
	    fi;
	    if Length(bl[i])>3 then 
# We cut the remaining pseudo-blocks in p-blocks by ultimate test 
# \sum_{\chi\in bl}\chi(T)/s_\chi p-integral
cut:=function(bl,para)local csch,lcm,lsch,p0,ct,ch,msch,l,Ah,getH;
  InfoChevie("#I p=",p," h",Position(hplanes,h),":",FormatGAP(h),
     " cut",FormatGAP(bl));
    getH:=function(para)local c,o,v;
      c:=Set(W.orbitRepresentative);
      o:=List(W.EigenvaluesGeneratingReflections{c},x->1/x);
      o:=List([1..Length(o)],i->Sum(o{[1..i-1]})+[1..o[i]]);
      v:=[];v{c}:=List(o,i->Zip(i,i-Minimum(i),
	function(x,y)return Mvp("x")^para[x]*E(Length(i))^y;end));
      return Hecke(W,v); # algebra A_h
    end;
  # replace para by smallest multiple such that schur elements rational
  para:=para*Lcm(List(Set(Concatenation(List(sch{bl},
           x->List(x.vcyc,y->y.mon))))*para,Denominator));
  Ah:=getH(para);
  csch:=List(sch{bl},s->s.coeff*CycPol(Mvp("x")^(para*s.mon))*
     Product(s.vcyc,x->CycPolOps.DescentOfScalars(x.pol,para*x.mon)));
  # csch holds the Schur elements for characters of Ah in bl
  lcm:=ApplyFunc(LcmCycPol,csch);
  lsch:=List(csch,x->lcm/x);
  InfoChevie(" Schur:",Stime(),"\c");
  lsch:=List(lsch,x->Value(x,Mvp("x")));
  InfoChevie(" Value:",Stime()," \c");
  ct:=TransposedMat(CharTable(Ah).irreducibles{bl});
  SortParallel(List(ChevieClassInfo(W).classtext,x->-Length(x)),ct);
  ct:=Filtered(ct,r->ForAll(r,x->IsCyc(x) or (not ForAny(x.coeff,IsUnknown) and
      Length(Variables(x))<=1)));
  if Length(ct)<>NrConjugacyClasses(W) then
    Print("\n!! Unreliable computation: ",CharTable(Ah)," partially unknown\n");
  fi;
  p0:=List([1..Length(bl)],x->[x]);
  for ch in ct do
    InfoChevie(".\c");
    msch:=Zip(lsch,ch,function(x,y)return x*y;end);
    l:=Filtered(Filtered(List(Combinations(p0),x->Set(Concatenation(x))),
       x->Length(x)>0),x->PValuation(Sum(msch{x}),p)>=0);
    if not [1..Length(bl)] in l then Error("theory");fi;
    l:=Filtered(l,x->Number(l,y->IsSubset(x,y))=1);
    p0:=LcmPartitions(l,p0);
    if Length(p0)=1 then InfoChevie(Stime()," ok\n");return [bl];fi;
  od;
  # here bl has been non-trivially cut 
  InfoChevie(Stime(),"\n  ->",FormatGAP(List(p0,x->bl{x})),"\n");
  return List(p0,x->bl{x});
end;
	      j:=cut(bl[i],para);
	      bl:=Concatenation(bl{[1..i-1]},j,bl{[i+1..Length(bl)]});
	      i:=i+Length(j);
	    else i:=i+1;
	    fi;
	  else i:=i+1;
	  fi;
	fi;
      od;
      bl:=Filtered(bl,x->x<>[]);Sort(bl);return bl;
    end);
# The h-blocks is the finest partition coarser than all p-blocks
    res:=ApplyFunc(LcmPartitions,res);
    if h=0*h then bl0:=res;fi;
    return [h,res];
  end);
end;

# The Rouquier blocks of a 1-cyclotomic algebra H is the finest partition 
# coarser than h-blocks for all h annihilated by H's parameters.
RouquierBlocks:=function(H)local W,d,p;
  W:=Group(H);
  d:=RouquierBlockData(W);
  p:=Concatenation(H.parameter{Set(W.orbitRepresentative)});
  d:=Filtered(d,function(x)local h;
    h:=x[1]*Lcm(List(x[1],Denominator));
    h:=Product(Zip(p,h,function(a,b)return a^b;end));
    return ScalMvp(h)<>false;
    end);
  return ApplyFunc(LcmPartitions,List(d,x->x[2]));
end;
