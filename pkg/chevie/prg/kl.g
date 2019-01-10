#############################################################################
##
#A  kl.g           CHEVIE library    Meinolf Geck, Jean Michel, Andrew Mathas
##
#Y  Copyright (C) 1992 - 2012  Lehrstuhl D f\"ur Mathematik, RWTH Aachen, IWR
#Y  der Universit\"at Heidelberg, University of St. Andrews, and   University 
#Y  Paris VII.
##
##  This file contains GAP functions to compute Kazhdan - Lusztig polynomials
##  and left cells, and  for  working with various bases of Hecke algebras. 
##
#############################################################################

#############################################################################
##
## An exemple of CreateHeckeBasis: create the 'T' basis for
## CoxeterHeckeAlgebraOps (i.e. all Hecke algebras for Coxeter groups...)
##
CreateHeckeBasis("T",rec(T:=x->x,     # method to convert to T
  inverse:=function(h)local getTinv,H;H:=Hecke(h);
    getTinv:=function(w)local i,q,W,res;
      if not IsBound(H.("T^-1")) then H.("T^-1"):=Dictionary();fi;
      res:=H.("T^-1").Get(w);
      if res<>false then return res;fi;
      W:=Group(H);
      if w=W.identity then res:=HeckeElt(H,"T",[W.identity],[1]);
      else i:=FirstLeftDescending(W,w); q:=H.parameter[i];
	res:=HeckeElt(H,"T",[W.identity,W.reflections[i]],
	  [Sum(q),-1]/Product(q))*getTinv(W.reflections[i]*w);
      fi;
      return H.("T^-1").Insert(w,res);
    end;
    if Length(h.elm)<>1 then Error("inverse implemented only for single T_w");
    else return h.coeff[1]^-1*getTinv(h.elm[1]^-1);fi;
  end),CoxeterHeckeAlgebraOps);

# RootParameter(H[,i][,msg])
# i is an integer in W.generatingReflections or an element of W.
# for a generator s, if (T_s-p1)(T_s-p2)=0 then H.rootParameter is sqrt(-p1p2)
# so for (T_s-q)(T_s+1)=0 it is sqrt(q) 
# and for (T_s-v)(T_s+v^-1)=0 it is 1
# Anyway the 2 parameters of T_s/RootParameter(W,s) satisfy p1*p2=-1
# For w \in W returns the product of RootParameter along a reduced expression.
RootParameter:=function(arg)local H,W,i,j,l;
  H:=arg[1]; W:=Group(H);
  if Length(arg)=1 then
    H.equal:=Length(Set(List([1..Length(H.parameter)],
                      i->RootParameter(H,i))))<=1;
    if not H.equal then Error("H should have equal parameters");fi;
    return RootParameter(H,1);
  fi;
  i:=arg[2];
  if IsInt(i) then
    if not IsBound(H.rootParameter[i]) then
      j:=W.rootRestriction[W.orbitRepresentative[i]];
      if Product(H.parameter[j])=-1 then 
        H.rootParameter[j]:=H.parameter[j][1];
      else
        l:=[-H.parameter[j][1]*H.parameter[j][2],2];
        if Length(arg)=3 then Add(l,arg[3]);fi;
        H.rootParameter[j]:=ApplyFunc(GetRoot,l);
      fi;
      for i in Filtered(W.generatingReflections,i->
         W.rootRestriction[W.orbitRepresentative[i]]=j) do
        H.rootParameter[i]:=H.rootParameter[j];
      od;
    fi;
    return H.rootParameter[i];
  elif IsBound(H.equal) and H.equal then 
    return H.rootParameter[1]^CoxeterLength(W,i);
  else return Product(CoxeterWord(W,i),
    y->H.rootParameter[Position(W.reflectionsLabels,y)]);
  fi;
end;

# The next routine  checks that the Hecke  algebra H is one of  a kind for
# which we know how to compute the KL basis:
#  According to  Lusztig, they are  defined when all  rootParameters are
# (or can  be) bound and are  elements of some totally  ordered group G.
# Then coefficients of Hecke elements are  elements of the ring Z[G]. We
# need to be able to perform the following operations in this ring:
#
#   H.Bar: the involution \sum_g a_g g -> \sum_g a_g g^-1
#   H.PositivePart:      \sum_g a_g g -> \sum_{g>=1} a_g g
#   H.NegativePart:      \sum_g a_g g -> \sum_{g<=1} a_g g
#
# We  know  how  to  construct  automatically  such  functions  when  all
# rootParameters  are powers  (GAP monomials)  of the  same variable  q.
# Otherwise they need  to be bound already (Positive for equal parameters
# and Negative in the general case: there is special code faster than the
# general case when all parameters are equal).

PrepareForPolynomials:=function(H)local q;
# prepare to work with 1-variable polynomials

  for q in H.rootParameter do # need parameters positive monomials
    if not IsPolynomial(q) or Length(q.coefficients)<>1 or
      q.coefficients[1]<>1 or Degree(q)<=0 then 
      return false;
    fi;
  od;

  H.Bar:=p->Polynomial(p.baseRing,Reversed(p.coefficients),
              1-p.valuation-Length(p.coefficients));

  # returns terms of positive degree of polynomial p
  H.PositivePart:=function(p)local v;
    v:=Maximum(0,-p.valuation);if v=0 then return p;fi;
    return Polynomial(p.baseRing,p.coefficients{[1+v..Length(p.coefficients)]},
      p.valuation+v);
  end;

  # returns terms of negative degree of polynomial p
  H.NegativePart:=p->Polynomial(p.baseRing,p.coefficients{[1..Minimum(
           Length(p.coefficients),1-p.valuation)]},p.valuation);

  return true;
end;

# to prepare Hecke algebras to use Mvps for KL theory
PrepareForMvp:=function(H)local q;
# prepare to work with Mvp -- lexicographic order

if VKCURVE.mvp=2 then
  for q in H.rootParameter do # need parameters positive monomials
    if not IsMvp(q) or Length(q.coeff)<>1 or q.coeff[1]<>1 or 
      Length(RecFields(q.elm[1]))=0 or 
      q.elm[1].(RecFields(q.elm[1])[1])<0 then return false;
    fi;
  od;

  H.Bar:=p->Mvp(List(p.elm,function(x)local res,f;res:=rec();
    for f in RecFields(x) do res.(f):=-x.(f);od;return res;end),p.coeff);
  
  H.PositivePart:=function(p)local v;
    v:=Filtered([1..Length(p.elm)],function(i)local f;i:=p.elm[i];
      f:=RecFields(i);return ForAll(f,ff->i.(ff)>0);end);
    return Mvp(p.elm{v},p.coeff{v});
  end;

  H.NegativePart:=function(p)local v;
    v:=Filtered([1..Length(p.elm)],function(i)local f;i:=p.elm[i];
      f:=RecFields(i);return ForAll(f,ff->i.(ff)<0);end);
    return Mvp(p.elm{v},p.coeff{v});
  end;
else
  for q in H.rootParameter do # need parameters positive monomials
    if not IsMvp(q) or Length(q.coeff)<>1 or q.coeff[1]<>1 or 
      Length(q.elm[1].elm)=0 or q.elm[1].coeff[1]<0 then 
      return false;
    fi;
  od;

  H.Bar:=p->Mvp(List(p.elm,x->rec(elm:=x.elm,coeff:=-x.coeff)),p.coeff);
  
  H.PositivePart:=function(p)local v;
    v:=Filtered([1..Length(p.elm)],x->p.elm[x].coeff=[] or p.elm[x].coeff[1]>0);
    return Mvp(p.elm{v},p.coeff{v});
  end;

  H.NegativePart:=function(p)local v;
    v:=Filtered([1..Length(p.elm)],x->p.elm[x].coeff=[] or p.elm[x].coeff[1]<0);
    return Mvp(p.elm{v},p.coeff{v});
  end;
fi;

  return true;
end;

CoxeterHeckeAlgebraOps.InitKL:=function(H,msg)
  if IsBound(H.Bar) then return;fi; # assume already called
  if not IsCoxeterGroup(Group(H))then Error(msg,": only for Coxeter groups");fi;
  H.equal:=Length(Set(List([1..Length(H.parameter)],
                      i->RootParameter(H,i,msg))))<=1;
  # this takes care rootParameter is filled
  if not IsBound(H.Bar) and not PrepareForPolynomials(H) 
    and not PrepareForMvp(H) then
      Error(msg,": parameters should be positive monomials",
                " or the function '.Bar' should be defined\n");
  fi;
end;

CoxeterHeckeAlgebraOps.InitD:=function(H,msg)
  CoxeterHeckeAlgebraOps.InitKL(H,msg);
  if not H.equal or not IsFinite(Group(H)) then 
    Error(msg,": only for equal parameters algebras of finite Coxeter groups");
  fi;
end;

##  Alt is The involution of H defined by v -> v^-1 and
##  T_w -> (-1)^{l(w)} RootParameter(w)^-2 T_w.
##  It swaps C_w and C'_w, and D_w and D'_w.
##  Essentially it corresponds to tensoring with the sign representation.

# The function  below is made  member of CoxeterHeckeAlgebraOps  just to
# hide it. It applies v->v^-1 to the coefficients of h. This is not Alt!
# One needs also to do some normalisation depending on the basis.
CoxeterHeckeAlgebraOps.Alt:=function(h,newbasis)local H;H:=Hecke(h);
  return HeckeElt(H,newbasis,h.elm,List(h.coeff,H.Bar));
end;

##  Beta is the involution on H defined by v -> v^-1 and
##  T_w -> v^{-l(w_0)}T_{w_0w}.
## It swaps C'_w and D_{w_0w}, and C_w and D'_{w_0w}.

# The function  below is made  member of CoxeterHeckeAlgebraOps  just to
# hide it. It is not the full Beta!
CoxeterHeckeAlgebraOps.Beta := function(h,newbasis)local H;H:=Hecke(h);
  return HeckeElt(H,newbasis,LongestCoxeterElement(Group(H))*h.elm,
     List(h.coeff,H.Bar));
end;

# now add operations to basis "T"

# The complete Alt for T basis
CoxeterHeckeAlgebraOps.T.AltInvolution:= function(h)local H;H:=Hecke(h); 
  H.operations.InitKL(H,"alt involution");
  return HeckeElt(H,"T",h.elm, Zip(h.elm,h.coeff,function(e,c)
      return RootParameter(H,e)^-2*(-1)^CoxeterLength(Group(H),e)*H.Bar(c);end));
end;

# The complete Beta for T basis
CoxeterHeckeAlgebraOps.T.BetaInvolution:=function(h)local H,W,w0;H:=Hecke(h);
  H.operations.InitKL(H,"beta involution");
  W:=Group(H);
  if not IsBound(W.isFinite) or not W.isFinite then 
    Error("BetaInvolution only for finite Coxeter groups");
  fi;
  w0:=LongestCoxeterElement(W);
  return HeckeElt(H,"T",w0*h.elm,RootParameter(H,w0)^-1*List(h.coeff,H.Bar));
end;

#############################################################################
#F  CriticalPair( <W>, <y>, <w> ) . . . . the critical pair (z,w)
#F   associated with the pair (y,w).
##
## Let  L  (resp.  R)  be  the  left  (resp.  right) descent set. A pair of
## elements  (y,w) of W  is called critical  if L(w) is  contained L(y) and
## R(w)  is  contained  in  R(y).  If  $(y,w)$  is not critical, <y> can be
## multiplied  from the left (resp. the right) by an element of L(w) (resp.
## R(w))  which is not in L(y) (resp. R(y)) until we get a new pair $(z,w)$
## critical. The function returns z. If y<=w then y<=z<=w.
##
## The  significance of  this  construction is  that the  Kahzdan-Lusztig
## polynomials P_{y,w} and P_{z,w} are equal.
##
CriticalPair:=function(W,y,w)local cr,Rw,Lw,IL;
  IL:=W.operations.IsLeftDescending; # avoid dispatching overhead.
  Lw:=W.rootRestriction{LeftDescentSet(W,w)};
  Rw:=W.rootRestriction{RightDescentSet(W,w)};
  cr:=function(y) local s;
    for s in Lw do if not IL(W,y,s) then return cr(W.reflections[s]*y);fi;od;
    for s in Rw do if not IL(W,y^-1,s) then return cr(y*W.reflections[s]);fi;od;
    return y;
  end;
  return cr(y);
end;

# The function below is made member of CoxeterHeckeAlgebraOps just to hide it
CoxeterHeckeAlgebraOps.getCp:=function(H,w)local W,iw,i,qx,x,z,s,res;
  if not IsBound(H.("C'->T")) then H.("C'->T"):=Dictionary();fi;
  res:=H.("C'->T").Get(w);
  if res<>false then return ShallowCopy(res);fi;
  W:=Group(H);
  if w=W.identity then res:=HeckeElt(H,"T",[W.identity],[H.unit]); # C'_1=1
# InfoChevie("# computed ",Length(H.("C'->T").keys)," C':",Stime());
  elif H.equal then
    if w in W.reflections{W.generatingReflections} then return
      HeckeElt(H,"T",[W.identity,w],[1/H.rootParameter[1],RootParameter(H,w)^-1]);
    fi;
    i:=FirstLeftDescending(W,w);
    s:=W.reflections[i];
    res:=H.operations.getCp(H,s)*H.operations.getCp(H,s*w);
	
    # new implementation JM and FD march 2000 replacing A.Mathas dec. 1994
    # We use the following formulae:
    # $$ C'_w=\sum_{y\le w}P_{y,w}(q)q^{-l(w)/2}T_y$$
    # and if sw<w then
    # $$ C'_s C'_{sw}=C'w+\sum_{y<sw}\mu(y,sw)C'y=\sum_{v\le w}\mu_v T_v$$
    # where
    # $$\mu_v=P_{v,w}(q)q^{-l(w)/2}+\sum_{v\le y\le sw}
    # \mu(y,sw)P_{v,y}(q)q^{-l(y)/2}$$
    #
    # It follows  that if  deg(\mu_v)>=-l(v)  then deg(\mu_v)=-l(v)  
    # with leading coefficient \mu(v,sw) 
    # (this happens exactly for y=v in the sum which occurs in
    # the formula for \mu_v).
    
    res:=res-Sum([1..Length(res.elm)],function(i)local c,e;
	e:=res.elm[i];if e=w then return 0;fi;
	c:=H.PositivePart(res.coeff[i]*RootParameter(H,e));
	if c<>0*c then return c*H.operations.getCp(H,e); else return 0; fi;
      end);
  else
    # we follow formula 2.2 in Lusztig's 'Left cells in Weyl groups'
    #
    # \bar P_{x,w}-P_{x,w}=\sum_{x<y\le w} R_{x,y} P_{y,w}
    #
    #  where R_{x,y}=\bar(T_{y^-1}^{-1}|T_x)
    #
    # thus we compute P_{x,w} by induction on l(w)-l(x) by
    # P_{x,w}=\neg \sum_{x<y\le w} R_{x,y} P_{y,w}
    res:=HeckeElt(H,"T",[w],[RootParameter(H,w)^-1]);
    res.elm:=Concatenation(Reversed(BruhatSmaller(W,w)));
    for i in [2..Length(res.elm)] do
      x:=res.elm[i];qx:=RootParameter(H,x);z:=CriticalPair(W,x,w);
      if x<>z then res.coeff[i]:=RootParameter(H,x)*
        res.coeff[Position(res.elm,z)]/qx;
      else res.coeff[i]:=-H.NegativePart(Sum([1..i-1],y->H.Bar(Coefficient(
        HeckeElt(H,"T",[res.elm[y]^-1],[1])^-1,x))*res.coeff[y])/qx)/qx;
      fi;
    od;
    CollectCoefficients(res);
  fi;
  return ShallowCopy(H.("C'->T").Insert(w,res));
end;

# The function below is made member of CoxeterHeckeAlgebraOps just to hide it
CoxeterHeckeAlgebraOps.ToKL:=function(h,target,index)local x,res,lens,H,coeff;
  H:=Hecke(h);res:=HeckeElt(H,target,[],[]); 

# To convert from "T", we use the fact that the transition matrix M from
# any  KL  bases to  the  standard  basis  is triangular  with  diagonal
# coefficient on  T_w equal  to RootParameter(w)^-1.  The transition  matrix is
# lower triangular for the C and  C' bases, and upper triangular for the
# D and D' bases which is what index(maximum or minimum) is for.

  while h.elm <> [] do
    lens:=List(h.elm,w->CoxeterLength(Group(H),w));
    x:=index(lens); lens:=Filtered([1..Length(h.elm)],i->lens[i]=x);
    coeff:=List(lens,i->h.coeff[i]*RootParameter(H,h.elm[i]));
    x:=HeckeElt(H,target,h.elm{lens},coeff);
    res:=res+x; h:=h-Basis(H,"T")(x);
  od;
  return res;
end;

CreateHeckeBasis("C'",rec(
  init:=function(H)H.operations.InitKL(H,"C' basis");end,
  BetaInvolution:=h->Hecke(h).operations.Beta(h,"D"),
  AltInvolution:=h->Hecke(h).operations.Alt(h,"C"),
  T:=function(h)local H;H:=Hecke(h);
  if Length(h.coeff)=0 then return Basis(H,"T")([],[]);
  else return h.coeff*List(h.elm,x->H.operations.getCp(H,x));fi;end,
  ("C'"):=h->Hecke(h).operations.ToKL(h,"C'",Maximum)
  ),CoxeterHeckeAlgebraOps);

CreateHeckeBasis("C",rec(
  init:=function(H)H.operations.InitKL(H,"C basis");end,
  BetaInvolution:=h->Hecke(h).operations.Beta(h,"D'"),
  AltInvolution:=h->Hecke(h).operations.Alt(h,"C'"),
  T:=function(h)local H;H:=Hecke(h);
  if Length(h.coeff)=0 then return Basis(H,"T")([],[]);
  else return h.coeff*List(h.elm,x->(-1)^CoxeterLength(Group(H),x)*
        H.operations.T.AltInvolution(H.operations.getCp(H,x)));fi;end,
  C:=h->Hecke(h).operations.ToKL(h,"C",Maximum)
  ),CoxeterHeckeAlgebraOps);

CreateHeckeBasis("D",rec(
  init:=function(H)H.operations.InitD(H,"D basis");end,
  BetaInvolution:=h->Hecke(h).operations.Beta(h,"C'"),
  AltInvolution:=h->Hecke(h).operations.Alt(h,"D'"),
  T:=function(h)local H;H:=Hecke(h);
  if Length(h.coeff)=0 then return Basis(H,"T")([],[]);
  else return h.coeff*List(h.elm,x->H.operations.T.BetaInvolution(
     H.operations.getCp(H,LongestCoxeterElement(Group(H))*x)));fi;end,
  D:=h->Hecke(h).operations.ToKL(h,"D",Minimum)
  ),CoxeterHeckeAlgebraOps);

CreateHeckeBasis("D'",rec(
  init:=function(H)H.operations.InitD(H,"D' basis");end,
  BetaInvolution:=h->Hecke(h).operations.Beta(h,"C"),
  AltInvolution:=h->Hecke(h).operations.Alt(h,"D"),
  T:=function(h)local H;H:=Hecke(h);
  if Length(h.coeff)=0 then return Basis(H,"T")([],[]);
  else return h.coeff*List(h.elm,x->(-1)^CoxeterLength(Group(H),x)*
      H.operations.T.AltInvolution(
        H.operations.T.BetaInvolution(H.operations.getCp(
	H,LongestCoxeterElement(Group(H))*x))));fi;end,
  ("D'"):=h->Hecke(h).operations.ToKL(h,"D'",Minimum)
  ),CoxeterHeckeAlgebraOps);

#############################################################################
##
#F  KazhdanLusztigPolynomial( <W>, <y>, <w>)  . . . Kazhdan-Lusztig polynomial
##
##  returns the list of coefficients of the Kazhdan-Lusztig polynomial P_{y,w}
##
KazhdanLusztigPolynomial:=function(W,y,w)local lw,s,v,z,lz,m,pol,ild,add;
  if not Bruhat(W,y,w) then return [];fi;
  y:=CriticalPair(W,y,w);
  lw:=W.operations.CoxeterLength(W,w);
  if lw-W.operations.CoxeterLength(W,y)<=2 then return [1];fi; 
  if not IsBound(W.klpol) then W.klpol:=Dictionary();fi;
  pol:=W.klpol.Get([w,y]);
  if pol<>false then return pol;fi;
  s:=FirstLeftDescending(W,w);v:=W.reflections[s]*w;
  add:=function(a,b)a:=ShallowCopy(a);AddCoeffs(a,b);
    NormalizeCoeffs(a);return a;end;
  pol:=add(KazhdanLusztigPolynomial(W,W.reflections[s]*y,v),
           ShiftedCoeffs(KazhdanLusztigPolynomial(W,y,v),1));
  lz:=lw-2;
  ild:=W.operations.IsLeftDescending;
  while (lw-lz)/2+1<=Length(pol) do
    for z in CoxeterElements(W,lz) do 
      if (lw-lz)/2+1<=Length(pol) and  pol[(lw-lz)/2+1]>0 and ild(W,z,s) and 
        Bruhat(W,y,z) then
        m:=KazhdanLusztigMue(W,z,v);
        if m<>0 then pol:=add(pol,
         -m*ShiftedCoeffs(KazhdanLusztigPolynomial(W,y,z),(lw-lz)/2));fi;
      fi;
    od;
    lz:=lz-2;
  od;
  return W.klpol.Insert([w,y],pol);
end;

#############################################################################
##
#F  KazhdanLusztigMue( <W>, <y>, <w>)  . . . . . . . . . . .
#F  . . . . . . . . . the highest coefficient of a Kazhdan-Lusztig polynomial
##
##  'KazhdanLusztigMue'  returns the coefficient of highest possible degree
##  of   the   Kazhdan-Lusztig   polynomial   of   the  Coxeter  group  <W>
##  corresponding to the elements <y> and <w>. That highest possible degree
##  is $(lw- ly -1)/2$. If <y> and <w> are not related by the Bruhat order,
##  the value 0 is returned.
##
KazhdanLusztigMue:=function(W,y,w) local ly,lw,pol,ild;
  ly:=W.operations.CoxeterLength(W,y); lw:=W.operations.CoxeterLength(W,w);
  if ly=lw or not Bruhat(W,y,w) then return 0; fi;
  if lw=ly+1 then return 1; fi; 
  ild:=W.operations.IsLeftDescending;
  if ForAny(W.generatingReflections,s->(ild(W,w,s) and not ild(W,y,s)) 
    or (ild(W,w^-1,s) and not ild(W,y^-1,s))) then return 0;
  fi;
  pol:=KazhdanLusztigPolynomial(W,y,w);
  if Length(pol)=(lw-ly+1)/2 then return pol[(lw-ly+1)/2];else return 0;fi;
end;

#############################################################################
##
#F  KazhdanLusztigCoefficient( <W>, <y>, <w>, <k> ) . . . . .
#F  . . . . . . . . . . . . . . . coefficient of a Kazhdan-Lusztig polynomial
##
##  returns the <k> th coefficient of the Kazhdan-Lusztig polynomial P_{y,w}.
##  return 0 unless y<=w for the Bruhat order.
##
KazhdanLusztigCoefficient:=function(W,y,w,k) local ly,lw,s,Lw,Rw,pol;
  ly:=W.operations.CoxeterLength(W,y); lw:=W.operations.CoxeterLength(W,w);
  if k<0 or not IsInt(k) or ly>lw or (ly<lw and 2*k>lw-ly-1) 
         or not Bruhat(W,y,w) then return 0; fi;
  if k=0 then return 1; fi;
  pol:=KazhdanLusztigPolynomial(W,y,w);
  if k+1<=Length(pol) then return pol[k+1];else return 0;fi;
end;

# all pairs [\mu(y,x),y] with non-zero mu
#muelist:=function(W,x)local v,H,t,l; v:=X(Rationals); H:=Hecke(W,[[v,-v^-1]]);
#  H.operations.InitKL(H,"muelist"); t:=H.operations.getCp(H,EltWord(W,x));
#  l:=Zip(t.coeff,t.elm,function(a,b)return [Coefficient(a,-1),Braid(W)(b)];end);
#  return Filtered(l,x->x[1]<>0);
#end;
#############################################################################
##
#F  KLMueMat( <W>, <list> )  . . . . (symmetrized) matrix of leading 
#F coefficients of Kazhdan-Lusztig polynomials of elements in a given list
##
KLMueMat:=function(W,c)local m,i,j,k,lc,n,w0c;
  w0c:=LongestCoxeterElement(W)*c; 
  n:=[1..Length(c)]; 
  lc:=List(c,i->CoxeterLength(W,i));
  m:=List(n,i->[]);
  for k in [0..Maximum(lc)-Minimum(lc)] do
    for i in n do
      for j in Filtered(n,j->lc[i]-lc[j]=k) do
	if lc[i]+lc[j]>W.N then m[i][j]:=KazhdanLusztigMue(W,w0c[i],w0c[j]);
	else m[i][j]:=KazhdanLusztigMue(W,c[j],c[i]);
	fi;
	m[j][i]:=m[i][j];
      od;
    od;
  od;
  return m;
end;

LeftCellOps:=OperationsRecord("LeftCellOps");

# A left cell must have the fields
# .reps representatives such that the left cell is the *-orbit of these.
# .group of which group it is a Left cell
# Optional (computed) fields are:
# .duflo the Duflo involution of the cell
# .character the character of the cell
# .a the a-function of the cell
# .elements the elements of the cell
LeftCellOps.Size:=function(c)
  if IsBound(c.character) then return
    Sum(List(CharTable(c.group).irreducibles,x->x[1]){c.character});
  else return Length(Elements(c));
  fi;
end;

# returns character as list of irred. numbers repeated if multiplicity
LeftCellOps.Character:=function(c)local r,ct,cc;
  if not IsBound(c.character) then
    r:=Representation(c,Hecke(c.group));
    cc:=CharRepresentationWords(r,ChevieClassInfo(c.group).classtext);
    ct:=CharTable(c.group);
    cc:=MatScalarProducts(ct,ct.irreducibles,[cc])[1];
    c.character:=Concatenation(List([1..Length(cc)],i->[1..cc[i]]*0+i));
    c.a:=ChevieCharInfo(c.group).a{c.character};
    if Length(Set(c.a))>1 then Error();else c.a:=c.a[1];fi;
  fi;
  return c.character;
end;

LeftCellOps.Print:=function(c)local f,p,uc;
  Print("LeftCell<",ReflectionName(c.group),": ");
  if IsBound(c.duflo) then 
    Print("duflo=",Join(DescribeInvolution(c.group,c.duflo)));
  fi;
  if IsBound(c.character) then
    uc:=UnipotentCharacters(c.group);
    f:=First(uc.families,f->uc.harishChandra[1].charNumbers[c.character[1]]
       in f.charNumbers); 
    f:=Position(uc.harishChandra[1].charNumbers,f.charNumbers[f.special]);
    p:=Position(c.character,f);
    p:=Concatenation([[f,1]],Collected(Drop(c.character,p)));
    Print(" character=",Join(List(p,function(v)
      local res; if v[2]<>1 then res:=String(v[2]);else res:="";fi;
      Append(res,CharNames(c.group)[v[1]]);return res;end),"+"));
  fi;
  Print(">");
end;

# r is one of BraidRelations(W); returns corresponding * op on w if applicable
# otherwise returns w
LeftStar:=function(W,r,w)local s,LeftStarNC;
  # st is a chain ststs of Length m_{s,t}
  # LeftDescentSet(W,w) contains st[1] and not st[2]
  # returns corresponding * operation applied to w
  LeftStarNC:=function(W,st,w)local w0,i,rst;rst:=W.reflections{st};
    w0:=Product(rst);i:=1;
    repeat w:=rst[i]*w;i:=i+1;w0:=w0*rst[i];
    until not W.operations.IsLeftDescending(W,w,st[i]);
    return w0*w;
  end;
  if W.operations.IsLeftDescending(W,w,W.rootRestriction[r[1][1]]) then
    if W.operations.IsLeftDescending(W,w,W.rootRestriction[r[2][1]]) then 
         return w;
    else return LeftStarNC(W,W.rootRestriction{r[1]},w);
    fi;
  elif W.operations.IsLeftDescending(W,w,W.rootRestriction[r[2][1]]) then
    return LeftStarNC(W,W.rootRestriction{r[2]},w);
  else return w;
  fi;
end;

# List of functions giving all possible left * images of w
LeftStars:=W->List(Filtered(BraidRelations(W),r->Length(r[1])>2),
  st->(w->LeftStar(W,st,w)));

LeftCellOps.CoxeterElements:=c->List(Elements(c),w->CoxeterWord(c.group,w));

LeftCellOps.Elements:=function(c)local w;
  if not IsBound(c.elements) then 
    c.elements:=FOrbit(LeftStars(c.group),c.duflo);
    for w in c.reps do Append(c.elements,FOrbit(LeftStars(c.group),w));od;
    c.elements:=Set(c.elements);
  fi;
  return c.elements;
end;

LeftCellOps.\=:=function(a,b)return a.duflo=b.duflo;end;

LeftCellOps.\in:=function(w,c)
  return w in Elements(c);end;

LeftCellOps.Mu:=function(c)
  if not IsBound(c.mu) then
    c.mu:=KLMueMat(c.group,Elements(c));
  fi;
  return c.mu;
end;

LeftCellOps.Representation:=function(c,H)local v,W,u,e,mu,w,res,l,k,s,n,value;
  W:=Group(H);
  if Length(Set(List(W.generatingReflections,
                i->RootParameter(H,i,"Left Cell Representation"))))>1 then
  # this checks one can compute rootParameters
    Error("cell representations for unequal parameters not yet implemented");
  else v := H.rootParameter[1];
  fi;
  return WGraphToRepresentation(c.group.semisimpleRank,WGraph(c),v);
end;

# returns right star operation st (a BraidRelation) of LeftCell c
RightStar:=function(st,c)local res,W,n;
  res:=ShallowCopy(c);Unbind(res.elements);
  W:=c.group;
  if IsBound(c.duflo) then
     res.duflo:=LeftStar(W,st,LeftStar(W,st,c.duflo^-1)^-1)^-1;
  fi;
  if IsBound(c.reps) then
     res.reps:=List(c.reps,w->LeftStar(W,st,w^-1)^-1);
  fi;
  if IsBound(c.elements) then
    res.elements:=List(c.elements,w->LeftStar(W,st,w^-1)^-1);
    n:=[1..Length(c.elements)];
    SortParallel(res.elements,n);
    if IsBound(c.mu) then res.mu:=c.mu{n}{n}; fi;
    if IsBound(c.graph) then res.orderGraph:=c.orderGraph{n};fi;
  fi;
  return res;
end;

# LeftCell containing w
LeftCell:=function(W,w)local word,v,g,sst,cell,l;
  l:=KLeftCellRepresentatives(W);
  sst:=Filtered(BraidRelations(W),r->Length(r[1])>2);
  word:=MinimalWordProperty(w,List(sst,st->(w->LeftStar(W,st,w^-1)^-1)),
     w->ForAny(l,c->w in c));
  v:=w;
  for g in Reversed(word) do v:=LeftStar(W,sst[g],v^-1)^-1;od;
  cell:=First(l,c->v in c);
  for g in word do cell:=RightStar(sst[g],cell);od;
  return cell;
end;

OldKLeftCellRepresentatives:=function(W)
  local st,rw,c,mu,n,Lleq,m,x,i,e,rd,d,ild;
  ild:=W.operations.IsLeftDescending; 
  st:=List(Filtered(BraidRelations(W),r->Length(r[1])>2),
     st->(c->RightStar(st,c)));
  rw:=CollectBy(Elements(W),x->RightDescentSet(W,x));
  rw:=List(rw,x->rec(rd:=RightDescentSet(W,x[1]),elements:=Set(x)));
  SortBy(rw,x->Length(x.elements));
  W.cells0:=[];
  while Length(rw)>0 do
    c:=rw[1].elements;
    InfoChevie("#I R(w)=",rw[1].rd," : #Elts=",Length(c),"\c");
    mu:=KLMueMat(W,c); n:=[1..Length(c)];
    Lleq:=List(n,x->List(n,y->x=y or (mu[x][y]<>0 and 
      ForAny(W.generatingReflections,i->ild(W,c[x],i) and not ild(W,c[y],i)))));
    Lleq:=TransitiveClosure(Lleq);
    m:=TransposedMat(Lleq);
    m:=Set(List(n,i->ListBlist(n,IntersectionBlist(Lleq[i],m[i]))));
    x:=List(m,d->rec(elements:=c{d},mu:=mu{d}{d},isDomain:=true,
      group:=W,operations:=LeftCellOps));
    while Length(x)>0 do
      c:=x[1];n:=[1..Length(c.elements)];SortParallel(c.elements,n);
      c.elements:=Set(c.elements); c.mu:=c.mu{n}{n};
      i:=Filtered(c.elements,x->x^2=W.identity);
      if Length(i)=1 then c.duflo:=i[1];
      else m:=List(i,x->CoxeterLength(W,x)
                           -2*Length(KazhdanLusztigPolynomial(W,(),x)));
        c.a:=Minimum(m);
        c.duflo:=i[Position(m,c.a)];# Duflo involutions minimize Delta
      fi;
      i:=Filtered(FOrbits(LeftStars(W),c.elements),x->not c.duflo in x);
      c.reps:=List(i,x->x[1]);
      Add(W.cells0,c);
      n:=FOrbit(st,c);
      InfoChevie(", ",Length(n)," new cell" );
      if Length(n)>1 then InfoChevie("s");fi;
      for e in n do
        rd:=RightDescentSet(W,e.duflo);
        i:=PositionProperty(rw,x->x.rd=rd);
        if i=1 then x:=Filtered(x,c->not c.elements[1] in Elements(e));
	elif i<>false then 
	  rw[i].elements:=Difference(rw[i].elements,Elements(e));
        fi;
      od;
    od;
    InfoChevie(" \n");
    rw:=Filtered(rw,x->Length(x.elements)>0);
    rw:=rw{[2..Length(rw)]};
    SortBy(rw,x->Length(x.elements));
  od;
  return W.cells0;
end;
  
#############################################################################
##
#F  LeftCells( <W> [,i] ) . . . left cells of W [in i-th 2-sided cell]
#   for the 1-parameter Hecke algebra
##
## 'LeftCells'  returns a list  of pairs. The  first component of each pair
## consists  of the reduced words in the Coxeter group <W> which lie in one
## left  cell C, the second component  consists of the corresponding matrix
## of highest coefficients mu(y,w), where y,w are in C.
## options: family= only leftcells in that uc.family
##
LeftCells:=function(arg)local W,ch,cc,opt,uc,st;
  W:=arg[1]; cc:=KLeftCellRepresentatives(W);
  if cc=false then cc:=OldKLeftCellRepresentatives(W);fi;
  if Length(arg)=2 then
    uc:=UnipotentCharacters(W);
    cc:=Filtered(cc,c->uc.harishChandra[1].charNumbers[Character(c)[1]]
       in uc.families[arg[2]].charNumbers);
  fi;
  st:=List(Filtered(BraidRelations(W),r->Length(r[1])>2),
     st->(c->RightStar(st,c)));
  return Union(List(cc,c->FOrbit(st,c)));
end;

# WGraph of LeftCell c
LeftCellOps.WGraph:=function(c)local e,mu,n,p,l,u,w,s,k,value,nodes;
  if not IsBound(c.graph) then 
    e:=Elements(c);mu:=LeftCellOps.Mu(c);n:=Length(e);
    nodes:=List(e,x->LeftDescentSet(c.group,x));
    p:=[1..n]; SortParallel(nodes,p);mu:=mu{p}{p};
    c.orderGraph:=p;
    nodes:=Concatenation(List(Collected(nodes),
      function(p)if p[2]=1 then return [p[1]];else return [p[1],p[2]-1];fi;
      end));
    c.graph:=[nodes,[]];
    l:=Concatenation(List([1..n],i->List([1..i-1],j->[mu[i][j],mu[j][i],i,j])));
    l:=Filtered(l,x->x[1]<>0 or x[2]<>0);
    l:=CollectBy(l,x->x{[1,2]});
    for u in l do
      if u[1][1]=u[1][2] then value:=u[1][1];else value:=u[1]{[1,2]};fi;
      w:=[value,[]];
      s:=CollectBy(List(u,x->x{[3,4]}),y->y[1]);
      for k in s do Add(w[2],Concatenation([k[1][1]],List(k,x->x[2])));od;
      Add(c.graph[2],w);
    od;
  fi;
  return c.graph;
end;

########################################################################
#               Functions for W-graphs 
# (Jean Michel june/december 2003 from  code/data of Geck, Marin, Alvis,
# Naruse, Howlett,Yin)
#    WGraphToRepresentation(semisimpleRank,graph,v)
# or WGraphToRepresentation(H,graph)
# Chevie stores some representations of some equal-parameter Hecke algebras
# as  $W$-graphs. For a Coxeter system $(W,S)$  a $W$-graph is defined by a
# set  of vertices  $C$; to  $x\in C$  is attached  $I(x)\subset S$  and to
# $(x,y)\in  C^2$  is  attached  an  ``edge''  $\mu(x,y)$  in  the field of
# definition  of $W$;  this defines  a representation  of the Hecke algebra
# with single rootparameter $v$ on a space with basis $e_y_{y \in C}$ by:
#
#  $$ T_s(e_y)=\cases{-e_y&                            if $s\in I(y)$\cr
#              v^2 e_y+\sum_{x\mid s\in I(x)} v\mu(x,y)e_x&otherwise\cr}$$
#
# The W-graphs  are stored in a  compact format to save  space. They are
# represented  as a  pair. 
# -The  first element is a list describing C; its elements are either a set
# I(x),  or an integer n  specifying to repeat the  previous element n more
# times.
# -The  second element is a list which  specifies mu. We first describe the
# mu-list  for symmetric  W-graphs (when  \mu(x,y)=\mu(y,x)). There  is one
# element  of the  mu-list for  each non-zero  value m  taken by \mu, which
# consists of a pair whose first element is m and whose second element is a
# list  of  lists;  if  l  is  one  of  these  lists  each pair [l[1],l[i]]
# represents  an  edge  (x=l[1],y=l[i])  such that \mu(x,y)=\my(y,x)=m. For
# non-symmetric  W-graphs, the first element of each pair in the mu-list is
# a  pair [m1,m2] and each edge [x,y] obtained from the lists in the second
# element has to be interpreted as mu(x,y)=m1 and mu(y,x)=m2.

# The next function given a W-graph gr for some Hecke algebra of rank rk
# with rootparameter v constructs the rk matrices it specifies
WGraphToRepresentation:=function(arg)local l,x,y,i,j,n,S,V,mu,H,rk,gr,v;
  gr:=arg[2];
  if IsInt(arg[1]) then
    rk:=arg[1];v:=arg[3];
    V:=[];
    for S in gr[1] do 
      if  IsInt(S) then Append(V,List([1..S],i->V[Length(V)]));
      else Add(V,S);
      fi;
    od;
    n:=Length(V);
    S:=List([1..rk],i->IdentityMat(n)*v^2);
    for j in [1..n] do for i in V[j] do S[i][j][j]:=-v^0;od;od;
    for i in gr[2] do 
      if IsList(i[1]) then mu:=i[1];else mu:=[i[1],i[1]];fi;
      for l in i[2] do 
        x:=l[1];
        for y in l{[2..Length(l)]} do
          for j in Difference(V[y],V[x]) do S[j][y][x]:=mu[2]*v;od;
          for j in Difference(V[x],V[y]) do S[j][x][y]:=mu[1]*v;od;
        od;
      od;
    od;
    return S;
  fi;
  H:=arg[1];
  S:=-H.parameter[1][2]*WGraphToRepresentation(Length(H.parameter),gr,
     RootParameter(H)/H.parameter[1][2]);
  CheckHeckeDefiningRelations(H,S);
  return S;
end;

############################################################################
# How to interpret W-graphs for complex reflection groups with one orbit of
# reflections, for Hecke(W,[vars]).

WGraph2Representation:=function(a,vars)local pos,nodes,n,dim,R,j,r,k;
  nodes:=a[1];
  pos:=function(n,j)local p;
    if IsList(n[1]) then p:=PositionProperty(n,x->j in x);
      if p=false then p:=Length(vars);fi;
    elif j in n then p:=1; else p:=2; fi;
    return p;
  end;
  n:=Maximum(Flat(nodes));# number of generators
  dim:=Length(nodes);
  R:=List([1..n],j->DiagonalMat(List([1..dim],k->vars[pos(nodes[k],j)])));
  for r in a[2] do for k in [3,4] do
    if IsList(r[k]) then
      for j in [2,4..Length(r[k])] do R[r[k][j-1]][r[k-2]][r[5-k]]:=r[k][j];od;
    else
      j:=Filtered([1..n],i->pos(nodes[r[k-2]],i)<pos(nodes[r[5-k]],i));
      R{j}[r[k-2]][r[5-k]]:=List(j,x->r[k]);
    fi;
  od;od;
  return R;
end;

# the next function returns the dual W-graph of gr (for an Hecke algebra of
# rank rk). A dual W-graph corresponds to a Curtis Dual representation.
DualWGraph:=function(rk,gr)
  return [List(gr[1],function(x)if IsInt(x) then return x;
                                else return Difference([1..rk],x);fi;end),
          List(gr[2],function(x)if IsList(x[1]) then return 
	  [-Reversed(x[1]),x[2]]; else return [-x[1],x[2]];fi;end)];
end;

############################################################################
# AsymptoticAlgebra(W,i) returns the asymptotic algebra with support the
# i-th two-sided cell of W
#
AsymptoticAlgebra:=function(W,i)local f,l,v,H,C,Cp,T,e,a,A,t,j,w0;
  l:=LeftCells(W,i);f:=Union(List(l,Character));a:=l[1].a;e:=List(l,Elements);
  for j in [1..Length(l)] do SortBy(e[j],x->CoxeterLength(W,x));od;
  e:=Concatenation(e);
  t:=List(l,x->Position(e,x.duflo));
  v:=X(Rationals);
  H:=Hecke(W,v^2,v);C:=Basis(H,"C");T:=Basis(H,"T");Cp:=Basis(H,"C'");
  A:=rec(field:=Rationals,
   operations:=OperationsRecord("AsympAlgebraOps",FDAlgebraOps),
   type:="Asymptotic algebra",
   parameters:=List(e,x->IntListToString(CoxeterWord(W,x))),
   basisname:="t");
  A.identification:=[A.type,i,W];
  A.zero:=AlgebraElement(A,[]);
  A.dimension:=Length(e);
  A.operations.underlyingspace(A);
  w0:=LongestCoxeterElement(W);
  # The algorithm below follows D. Alvis, "Subrings of the asymptotic Hecke
  # algebra of type H4" Experimental Math. 17 (2008) 375--383
  A.structureconstants:=List(e,x->List(e,function(y)local F,lx,ly,sc;
#   InfoChevie(".\c");
    F:=T(x)*T(y); lx:=CoxeterLength(W,x); ly:=CoxeterLength(W,y);
    sc:=List(e,function(z)local c,ez,lz,s; z:=z^-1;c:=T(Cp(w0*z));
      lz:=CoxeterLength(W,z);ez:=(-1)^lz;
      s:=ez*Sum([1..Length(F.elm)],function(i)local w,lw;
        w:=F.elm[i];lw:=CoxeterLength(W,w);
        return Coefficient(c,w0*w)*Coefficient(F,w)*(-1)^lw;end);
#      Print("deg=",[Valuation(s),Degree(s)]," pdeg=",a+lx+ly-W.N,"\n");
      return [Coefficient(s,a+lx+ly-W.N),Position(e,z^-1)];end);
    sc:=Filtered(sc,x->x[1]<>0);
    return sc;
  end));
  A.multiplication:=function(i,j)return A.structureconstants[i][j];end;
  A.operations.Print:=function(A)Print(A.type," dim.",A.dimension);end;
  A.operations.underlyingspace(A);
  A.one:=Sum(A.basis{t});
  A.operations.CharTable:=function(A)local tbl;
    tbl:=rec(domain:=A,field:=A.field,operations:=rec(Print:=TablePrint));
    tbl.irreducibles:=TransposedMat(List(e,x->List(HeckeCharValues(C(x)){f},
     p->(-1)^a*Coefficient(p,-a))));
    tbl.basistraces:=tbl.irreducibles;
    tbl.matrix:=tbl.irreducibles;
    tbl.characterDegrees:=List(tbl.matrix,x->Sum(x{t}));
    tbl.columns:=A.parameters;
    tbl.rows:=CharNames(W){f};
    return tbl;
  end;
  return A;
end;
############################################################################
# Virtual representations \cA_w and \ca_w of [Lus85] 5.10.2 and 5.11.6
Lusztigaw:=function(W,w)local v,l;
  v:=Indeterminate(Rationals);
  l:=HeckeCharValues(Basis(Hecke(W,v^2,v),"T")(w))*(-v)^-CoxeterLength(W,w);
  return Zip(l,ChevieCharInfo(W).a,function(c,a)return Coefficient(c,-a);end);
end;
LusztigAw:=function(W,w)local v,l;
  v:=Indeterminate(Rationals);
  l:=HeckeCharValues(Basis(Hecke(W,v^2,v),"T")(w))*v^-CoxeterLength(W,w);
  return Zip(l,ChevieCharInfo(W).A,function(c,A)return Coefficient(c,W.N-A);end);
end;
