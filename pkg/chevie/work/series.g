# series.g J.Michel 4/2003
# This file contains code to work with d-Harish-Chandra series of Spetses
#------------------------- Utility functions ---------------------------

# CuspidalUnipotentCharacters(WF[,d]) indices of the unipotent characters
# of the Spets WF which are d-cuspidal (d=0 if not specified)
CuspidalUnipotentCharacters:=function(arg)local WF,d,ad,ud;
  WF:=arg[1];
  if Size(WF)=1 then return [1];fi;
  if Length(arg)>1 then d:=arg[2];else d:=0;fi;
  ad:=Number(RelativeDegrees(WF,d),x->x<>1);
# if ad=0 then Error(d," should divide one of the degrees");fi;
  ud:=CycPolUnipotentDegrees(WF);
  return Filtered([1..Length(ud)],i->ad=Valuation(ud[i],d));
end;

# CuspidalPairs(W[,d[,ad]]) returns the pairs (LF,\lambda) where LF is a d-split
# Levi [with d-center of dimension ad] and lambda a d-cuspidal character of LF 
CuspidalPairs:=function(arg)local W,d,ad,WF;
  W:=arg[1];
  if IsSpets(W) then WF:=W;W:=Group(WF); else WF:=Spets(W); fi;
  if Length(arg)=1 then arg[2]:=0; fi;
  d:=arg[2];
  if Length(arg)=2 then 
    return Concatenation(List([0.. Length(RelativeDegrees(WF,d))],
      ad->CuspidalPairs(WF,d,ad)));
  fi;
  ad:=arg[3];
  if IsInt(d) and d<>0 then d:=1/d;fi;
  return Concatenation(List(SplitLevis(WF,d,ad),HF->
    List(CuspidalUnipotentCharacters(HF,d),char->[HF,char])));
end;

# Permutation of the conjugacy classes induced by an automorphism of W
PermutationOnClasses:=function(W,aut)
  return PermList(List(ConjugacyClasses(W),
     c->PositionClass(W,Representative(c)^aut)));
end;

# Permutation of the characters induced by a perm. automorphism of W
# PermutationOnCharacters(W,aut [,charlist])
PermutationOnCharacters:=function(arg)local W,p,ct;
  W:=arg[1]; p:=PermutationOnClasses(W,arg[2]);
  ct:=CharTable(W).irreducibles;if Length(arg)=3 then ct:=ct{arg[3]};fi;
  return PermListList(ct,List(ct,r->Permuted(r,p)));
end;

# Permutation of the unipotent characters induced by an automorphism of W
# PermutationOnUnipotents(W,aut [,uniplist])
PermutationOnUnipotents:=function(arg)local aut,W,p,t,l,uc;
  aut:=arg[2];W:=arg[1];
  uc:=UnipotentCharacters(W);
  if Length(arg)=3 then l:=arg[3];else l:=[1..Size(uc)];fi;
  t:=List(DeligneLusztigCharacterTable(W),x->x{l});
  Add(t,Eigenvalues(uc,l));t:=TransposedMat(t);
  if Length(Set(t))<Length(t) then 
    t:=List(l,x->Position(uc.harishChandra[1].charNumbers,x));
    if ForAll(t,x->x<>false) then return PermutationOnCharacters(W,aut,t);
    else Error("Rw + eigen cannot disambiguate\n");
    fi;
  fi;
  return PermListList(t,List(t,r->Permuted(r,PermutationOnClasses(W,aut))));
end;

# s is a Set of tuples. Return E_1,...,E_n such that
# s=List(Cartesian(E_1,...,E_n),Concatenation)
# Assumes all E_i but one are of size 2
FactorsSet:=function(s)local i,j,s1,s2,r,a;
  if s=[[]] then return [];fi;
  for i in s[1] do
    s1:=Filtered(s,y->i in y);
    if Length(s1)=Length(s)/2 then
      i:=Intersection(s1);
      r:=Difference(s[1],i);
      s2:=Filtered(s,y->Length(Intersection(y,i))=0);
      a:=PositionProperty(s2,y->IsSubset(y,r));
      if Length(s2)=Length(s)/2 and a<>false then j:=Difference(s2[a],r);
	s1:=Set(List(s1,x->Difference(x,i)));
	s2:=Set(List(s2,x->Difference(x,j)));
	if Length(i)=1 then i:=i[1];fi;
	if Length(j)=1 then j:=j[1];fi;
	if s1=s2 then return Concatenation([[i,j]],FactorsSet(s1));fi;
      fi;
    fi;
  od;
  return [s];
end;

LIMSubsetsSum:=10;
# l is a matrix and S a list of same length as a row of l.
# Find subsets P of [1..Length(l)] such that Sum(l{P})=S.
# in  addition, lv is a vector of same  length as l, v is a sub-multiset of
# lv and the chosen subsets should satisfy lv{P}=v as multisets.
SubsetsSum:=function(S,l,v,lv)local c,inner,found,sievev;
  sievev:=function(good,v)local i,p;
    for i in good do 
      p:=Position(v,lv[i]); if p=false then return [];fi;
      v:=v{Filtered([1..Length(v)],i->i<>p)};
    od;
    return v;
  end;
  c:=0;found:=0;
  # s assumed to be in P at this stage
  # S= initial S minus Sum(l{s})
  # t= remaining elements of [1..Length(l)] which could be in P
  # nonsolved= indices of nonsolved entries of S
  # v= remaining v to match
  inner:=function(S,s,t,nonsolved,v,factor)
    local bad,good,p,sols,res,i,sol,f,ll,solved;
  # Print("#solved=",Length(l[1])-Length(nonsolved)," ",Join(s),"=>",Join(t),
  #       " v=",Join(List(Collected(v),x->Join(x,":"))," "),"\n");
    c:=c+1;
    if c mod 1000=0 then
      InfoChevie("# ",factor,": xcols:",Length(nonsolved)-1," xrows:",Length(t),
            " found:",found,"\n");
    fi;
    t:=Filtered(t,i->lv[i] in v);
    if Length(t)=0 then 
      if S=S*0 then found:=found+1; return [[]];
      else return [];
      fi;
    fi;
    ll:=List(nonsolved,i->rec(pos:=i,cand:=Filtered(t,j->l[j][i]<>0)));
    SortBy(ll,x->Length(x.cand));
    if Length(ll[1].cand)>LIMSubsetsSum then ll:=[ll[1]];
    else ll:=Filtered(ll,x->Length(x.cand)<=LIMSubsetsSum);fi;
    solved:=[];good:=[];bad:=[];
    for p in ll do
      p.sols:=Filtered(Combinations(p.cand),e->Sum(l{e}[p.pos])=S[p.pos]);
      if Length(p.sols)=0 then return [];
      elif Length(p.sols)=1 then Add(solved,p.pos);
      fi;
      good:=Union(good,Intersection(p.sols)); #lines part of any solution
      bad:=Union(bad,Difference(p.cand,Union(p.sols)));#part of no solution
    od;
    nonsolved:=Difference(nonsolved,solved);
    if Length(good)+Length(bad)>0 then # progress
      return List(inner(S-Sum(l{good}),Union(s,good),
        Difference(t,Union(good,bad)),
        nonsolved,sievev(good,v),factor),r->Concatenation(good,r));
    else
      res:=[];
#     p:=Minimum(List(ll,x->Length(x.sols)));p:=First(ll,x->Length(x.sols)=p);
      p:=Maximum(List(ll,x->Length(x.cand)));p:=First(ll,x->Length(x.cand)=p);
      f:=Length(p.sols);
      nonsolved:=Difference(nonsolved,[p.pos]);
      InfoChevie("# ",factor,": xcols:",Length(nonsolved)," xrows:",Length(t),
	 " in comb(",Length(p.cand),")=>",Length(p.sols),"\n");
#     if Length(p.sols)>15 then Error();fi;
      for sol in p.sols do
	good:=sol;
	bad:=Difference(p.cand,sol);
	if Intersection(good,bad)=[] then 
         Append(res,List(inner(S-Sum(l{good}),Union(s,good),
	   Difference(t,Union(good,bad)),nonsolved,sievev(good,v),
	   SPrint(factor,":",f)),r->Concatenation(good,r)));
	fi;
	f:=f-1;
      od;
      return res;
    fi;
  end;
  return inner(S,[],[1..Length(l)],[1..Length(S)],v,"");
end;
  
#test:=function(S,l,v,lv)local groups,fuse,sure;
#  groups:=List([1..Length(l)],x->[x]);lv:=List(lv,x->[x]);sure:=[];
#  fuse:=function()local try,ll,p,part,good,bad,i,j;
#  try:=Filtered([1..Length(l[1])],i->Number([1..Length(l)],j->l[j][i]<>0) in [1..10]);
#  ll:=List(try,i->rec(pos:=i,cand:=Filtered([1..Length(l)],j->l[j][i]<>0)));
#  SortBy(ll,x->Length(x.cand));
#  for p in ll do
#    p.sols:=Filtered(Combinations(p.cand),e->Sum(l{e}[p.pos])=S[p.pos]);
#  od;
#  for i in [Length(ll),Length(ll)-1..1] do for j in [i-1,i-2..1] do
#    if IsSubset(ll[i].cand,ll[j].cand) then
#      ll[i].sols:=Filtered(ll[i].sols,s->Intersection(s,ll[j].cand) in ll[j].sols);
#    fi;
#  od;od;
#  for p in ll do
#    good:=Intersection(p.sols);
#    if Length(good)>0 then
#      S:=S-Sum(l{good});
#      UniteSet(sure,Union(groups{good}));
#      v:=DifferenceMultiSet(v,Concatenation(lv{good}));
#      Print("*",Length(v)+Length(sure),"*");
#      l:=Drop(l,good);groups:=Drop(groups,good);lv:=Drop(lv,good);
#      return true;
#    fi;
#    bad:=Difference(p.cand,Union(p.sols));
#    if Length(bad)>0 then
#      l:=Drop(l,bad);groups:=Drop(groups,bad);lv:=Drop(lv,bad);
#      return true;
#    fi;
#    part:=ApplyFunc(GcdPartitions,List(p.sols,x->[x,Difference(p.cand,x)]));
#    part:=Filtered(part,x->Length(x)>1);
#    if Length(part)>0 then
#      part:=part[1];i:=part{[2..Length(part)]};
#      l[part[1]]:=Sum(l{part});l:=Drop(l,i);
#      groups[part[1]]:=Union(groups{part});groups:=Drop(groups,i);
#      lv[part[1]]:=Concatenation(lv{part});lv:=Drop(lv,i);
#      return true;
#    fi;
#  od;
#  return false;
#  end;
#  while fuse() do Print(Length(l),"\n");od;
#  Error();
#end;

# FitParameter(sch,m,d) given:
# - A list sch of length e of schur elements for H(Z/e) given as cycpols 
# - A list m of length e of rationals 
# - A root of unity given as d in Q/Z
#
# find all permutations s of [1..e] such that the list of parameters 
#   p_k:=E(e)^(s(k)-1)(q/E(d))^m_k
# gives the given sch (by  the formula sch[k]:=\prod_{j\ne k}(1-p_k/p_j). 
# Each permutation is represented modulo a translation (can't be detected). 
#
# The  result is a list of pairs [v1,v2] telling that the elements of sch
# of indices in v1 are to be mapped to indices v2. The list of v1 is the
# same as CollectBy([1..Length(sch)],List(sch,x->Position(sch,x)));
#
FitParameter:=function(sch,m,d)
  local e,res,a,v,p,den,poss,pp,G,para,positive,bad,good;
  sch:=List(sch,x->CycPolOps.EnnolaTwist(x,E(Denominator(d))^Numerator(d)));
  positive:=p->ForAll(p.vcyc,x->x[2]>0); # tells if a CycPol is really a pol
  den:=Lcm(List(m,Denominator)); 
  e:=Length(m); a:=Collected(m);
  sch:=List(sch,x->CycPolOps.DescentOfScalars(x,den));
  poss:=function(i)local v,k,avail,good,p,term;
  # for each element of Collected(m) tries to find a minimal
  # corresponding possible set of s(i)-s(k) for k such that m_k=m
    v:=sch[i];
    avail:=[1..Length(a)];p:=List(avail,k->[0..e-1]);good:=[];
    term:=function(j,k)
      return CycPol(1-E(e)^j*X(Cyclotomics)^(den*(m[i]-a[k][1])));end;
    while true do
      for k in avail do 
        p[k]:=Filtered(Difference(p[k],Concatenation(p{good})),
                                   j->m[i]=a[k][1] or positive(v/term(j,k)));
      od;
      good:=Filtered(avail,i->Length(p[i])=a[i][2]);
      if Length(good)=0 then return p;fi;
      avail:=Difference(avail,good);
#     Print("good=",good,"\n");
      v:=v/Product(Filtered(good,k->m[i]<>a[k][1]),k->Product(p[k],j->term(j,k)));
    od;
  end;
  pp:=List([1..e],poss);
  p:=PositionProperty([1..e],i->List(a,x->x[2])=List(pp[i],Length));
  if p=false then return false;fi;
  v:=List([0..e-1],i->List(i+pp[p],x->Set(List(x,y->y mod e))));
  v:=List([1..e],k->Filtered([1..e],i->ForAll([1..Length(a)],j->IsSubset(
    pp[k][j],v[i][j]))));
  G:=List(sch,x->Position(sch,x));
  v:=List(Set(G),x->v[x]);
  G:=List(G,x->Position(Set(G),x));
  bad:=[1..Length(v)];
  repeat
    good:=Filtered(bad,i->Length(v[i])=Number([1..e],j->G[j]=i));
    bad:=Difference(bad,good);
    for p in bad do for a in good do v[p]:=Difference(v[p],v[a]);od;od;
  until Length(good)=0;
  if Sum(v,Length)>e then ChevieErr("non-unique solution\n");return false;fi;
  res:=List([1..Length(v)],x->[]);
  for p in [1..e] do Add(res[G[p]],p);od;
  pp:=[];for p in [1..Length(res)] do pp{res[p]}:=v[p];od;
  para:=List([1..e],i->E(e)^pp[i]*X(Cyclotomics)^(den*m[i]));
  if List([1..e],k->Product(Filtered([1..e],j->j<>k),
   j->CycPol(1-para[k]/para[j])))<>sch then 
    ChevieErr("schur elms don't match\n");return false;fi;
  return TransposedMat([CollectBy([1..Length(G)],G),v]);
end;

#---------------------- SeriesOps -----------------------------------------
# A d-Harish-Chandra series \CE(\BG,(\BL,\lambda) is a record with fields:
#  .d            AsRootOfUnity(\zeta) such that L=C_G(V_\zeta)
#  .spets        \BG
#  .levi         \BL
#  .classno      class of s.levi with \zeta-eigenspace V_\zeta
#  .element      representative of classno
#  .projector    .element-equivariant projector on V_\zeta
#  .cuspidal     \lambda (index in UnipotentCharacters(.levi))
#  .principal    true iff \lambda=\Id
#  .WGL          W_\BG(\BL,\lambda)
#                as a relgroup, contains parentMap (refs->elts of W)
#                and reflists (generators->gens of parab of W)
#  .WGLdims      irr dims of WGL
#  .charNumbers  \CE(\BG,(\BL,\lambda) (indices in UnipotentCharacters(.spets))
#  .degree       The degree of R_\BL^\BG(\lambda)=|G/L|_{q'} deg\lambda
#  .eps          For each \chi in .charNumbers sign of <\chi,R_\BL^\BG(\lambda)>
#  .dims         \chi(1) for \gamma_\chi
#  .Hecke        H_G(L,\lambda)
#
# Fields only present when WGL is cyclic
# .e         |WGL|

SeriesOps:=OperationsRecord("SeriesOps");

IsSeries:=s->IsRec(s) and IsBound(s.operations) and s.operations=SeriesOps;

Series:=function(spets,levi,cuspidal,d)local s,q,eig,e;
  if IsInt(d) and d>0 then d:=Mod1(1/d);fi;
  if not IsSpets(spets) then spets:=Spets(spets);fi;
  s:=rec(spets:=spets,levi:=levi,cuspidal:=cuspidal,d:=d,operations:=SeriesOps);
  s.principal:=UnipotentCharacters(s.levi).a[s.cuspidal]=0
               and UnipotentCharacters(s.levi).A[s.cuspidal]=0;
# find simplest regular eigenvalue q of s.levi
  eig:=Union(List(RegularEigenvalues(s.levi),function(x) 
    if IsInt(x) then return PrimeResidues(x)/x;else return [x];fi;end));
  e:=Minimum(List(eig,Denominator));
  e:=Minimum(Filtered(eig,x->Denominator(x)=e));
  eig:=ReflectionEigenvalues(s.levi);
  q:=Maximum(List(eig,x->Number(x,y->y=e)));
  s.classno:=PositionsProperty(eig,x->Number(x,y->y=e)=q);
  if Length(s.classno)>1 then Error("classno=",s.classno);fi;
  s.classno:=s.classno[1];
  s.element:=EltWord(s.levi,ChevieClassInfo(s.levi).classtext[s.classno]);
  q:=Minimum(List(ReflectionEigenvalues(s.levi),x->Number(x,y->y=d)));
  q:=PositionsProperty(ReflectionEigenvalues(s.levi),x->Number(x,y->y=d)=q);
  s.projector:=EigenspaceProjector(s.spets,
      Representative(ConjugacyClasses(s.levi)[q[1]]),d);
  q:=X(Cyclotomics);
  s.degree:=GenericSign(s.spets)^-1*GenericOrder(s.spets,q)/
           (GenericSign(s.levi)^-1* GenericOrder(s.levi,q));
  s.degree.valuation:=0;
  s.degree:=CycPol(s.degree)*CycPolUnipotentDegrees(s.levi)[s.cuspidal];
  return s;
end;

CuspidalSeries:=function(arg)local d;
  if Length(arg)=1 then d:=1;else d:=arg[2];fi;
  return List(ApplyFunc(CuspidalPairs,arg),x->Series(arg[1],x[1],x[2],d));
end;

SeriesOps.String:=s->Format(s,rec());

SeriesOps.Print:=function(s)Print(String(s));end;

SeriesOps.Display:=function(s,opt)
  opt.screenColumns:=SizeScreen()[1];
  Print(Format(s,opt));end;

SeriesOps.Format:=function(s,opt)local m,f,uw,e,t,rowlab,res,n,cname;
  cname:=CharNames(UnipotentCharacters(s.levi),opt)[s.cuspidal];
  if IsBound(opt.TeX) then n:="\\lambda";else n:="c";fi;
  if s.spets=s.levi then 
    res:=SPrint(s.d,"-cuspidal ",cname," of ",ReflectionName(s.spets,opt));
  else
  res:=SPrint(s.d,"-series ");
  if IsBound(opt.TeX) then PrintToString(res,"$");fi;
  PrintToString(res,"R^{",ReflectionName(s.spets,opt),"}_",
   "{",ReflectionName(s.levi,opt),"}(");
  if IsBound(opt.TeX) then PrintToString(res,"\\lambda=");fi;
  PrintToString(res,cname,")");
  fi;
  if IsBound(s.e) then 
    if IsBound(opt.TeX) then 
      PrintToString(res,"\\quad W_G(L,\\lambda)=Z_{",s.e,"}");
    else PrintToString(res," W_G(L,",n,")=Z",s.e);fi;
  elif IsBound(s.WGL)then 
    if IsBound(opt.TeX) then PrintToString(res,"\\quad ");fi;
    if IsBound(s.WGL.reflections) then
      PrintToString(res," W_G(L,",n,")=",ReflectionName(s.WGL,opt));
    else
      PrintToString(res," |W_G(L,",n,")|=",Size(s.WGL));
    fi;
  elif IsBound(s.WGLdims) then 
    if IsBound(opt.TeX) then PrintToString(res,"\\quad ");fi;
      PrintToString(res," |W_G(L,",n,")|=",Sum(s.WGLdims,x->x^2));
  fi;
  if IsBound(s.translation) then
    if IsBound(opt.TeX) then PrintToString(res,"\\quad ");fi;
    PrintToString(res," translation=",s.translation);
  fi;
  if IsBound(opt.TeX) then PrintToString(res,"$");fi;
  if not IsBound(opt.screenColumns) then return res;fi;
  uw:=UnipotentCharacters(s.spets);
  e:=Length(s.charNumbers);
  f:=function(arg)local vals;
     if IsBound(opt.TeX) then Add(rowlab,arg[1]);vals:=arg[Length(arg)];
     else Add(rowlab,arg[2]); vals:=arg[3];fi;
     Add(m,List(vals,x->Format(x,opt)));
  end;
  m:=[];rowlab:=[];
  f("\\hbox{Character}","Name",CharNames(uw,opt){s.charNumbers});
# f("\\hbox{degree}","degree",CycPolUnipotentDegrees(s.spets){s.charNumbers});
  f("\\hbox{specializes}","specializes",CharNames(s.WGL,opt));
  f("\\varepsilon","eps",s.eps);
  if IsBound(s.eigen) then f("\\hbox{eigen}","eigen",s.eigen);fi;
  if IsBound(s.e) then
    if IsBound(s.Hecke) then f("\\hbox{parameter}","param",s.Hecke.parameter[1]);
    elif IsBound(s.mC)  then f("\\hbox{mC}","degparam",s.mC);fi;
  fi;
  f("\\hbox{family \\#}","family #",
   List(s.charNumbers,j->PositionProperty(uw.families,k->j in k.charNumbers)));
  if IsBound(s.permutable) and ForAny(s.permutable,x->x<>false) then
    f("\\hbox{permutable}","permutable",s.permutable);
  fi;
  m:=TransposedMat(m);
  opt.rowLabels:=s.charNumbers;
  opt.columnLabels:=rowlab;
  PrintToString(res,"\n",FormatTable(m,opt));
  return res;
end;

# Degree in q of the parameters (normalized so the smallest is 0)
SeriesOps.mC:=function(s)local D0,aA,uc,pG,pL,xiL,xiG,e,cn,lpi;
  e:=HyperplaneOrbits(s.WGL);
  if Length(e)>1 then return;
  else e:=e[1].e_s;fi;
  uc:=UnipotentCharacters(s.spets);
  cn:=s.charNumbers{Filtered([1..Length(s.dims)],i->s.dims[i]=1)};
  aA:=uc.a{cn}+uc.A{cn};
  lpi:=W->Sum(ReflectionDegrees(W)+ReflectionCoDegrees(W));
  if s.principal then # Id in L
    if Minimum(aA)<>0 then Error("id not in RLG(1)");fi;
    pG:=lpi(Group(s.spets)); pL:=lpi(Group(s.levi));
    D0:=pG-pL;
    # compute xi_L
    xiL:=AsRootOfUnity(PhiOnDiscriminant(s.levi))*Denominator(s.d);
    xiG:=AsRootOfUnity(PhiOnDiscriminant(s.spets))*Denominator(s.d);
    if xiL<>xiG then 
      ChevieErr("fixing dimension of variety by xiL-xiG=",xiL-xiG,"\n");
    D0:=D0+xiL-xiG;
    fi;
    # Id in chars
    if not IsInt(D0*s.d) then 
      ChevieErr(s,"\n    => (l(\pi_G)=",pG,")-(l(\pi_L)=",pL,")+(xiL*d=",
      xiL,")-(xiG*d=",xiG,") not =0(mod d=",s.d,
      ") Max(aA)-Min(aA)=",Maximum(aA)-Minimum(aA),"\n");
    fi;
    if not IsInt(D0/e) then 
      ChevieErr("fixing dimension of variety to be divisible by e=",e,"\n");
      D0:=D0-(D0 mod e);
    fi;
    s.mC:=(D0-aA)/e;return s.mC;
  else
    # (JM+Gunter 18-3-2004) in any case we normalize so that
    # the smallest mC is 0 since the above choice gives unpleasant
    # results for G27
    D0:=Maximum(aA)-Minimum(aA); # just so that mC are positive
    s.mC:=(D0+Minimum(aA)-aA)/e;
    return s.mC;
  fi;
end;


if false then
# Find W_G(L,\lambda)
SeriesOps.RelativeGroup:=function(s) 
  local refs,rH,r,c,H,rel,W,ud,eig,WF,L,rel,NLF,NLFb,N,hom,phi,V,m,func,
        mats,reflists,id;
  if IsBound(s.WGL) then return s.WGL;fi;
  W:=Group(s.spets);
  refs:=Set(List(Reflections(W),x->Position(W.reflections,x)));
  mats:=[];
  mats{refs}:=List(refs,r->MatXPerm(W,Reflection(W,r)));
  id:=s.projector^0;
  if s.projector=id then #central series
    s.WGL:=Copy(W); s.WGL.parentMap:=s.WGL.generators;
    s.WGL.reflists:=List(W.rootInclusion{W.generatingReflections},x->[x]);
  else
    V:=NullspaceMat(s.projector-id);  # The E(d)-eigenspace
    rH:=List(refs,r->SumIntersectionMat(NullspaceMat(mats[r]-id),V)[2]);
    rH:=Set(Filtered(rH,H->Length(H)<Length(V)));
  # Print("rH=",rH,"\n"); # Hyperplanes of W_G(L)
    rel:=[];
    for H in rH do
      if Length(H)=0 then # cyclic WGL
	ConjugacyClasses(Group(s.levi)); # hack to fix a (bug in GAP?)
	H:=W;
      else H:=Filtered(refs,r->H*mats[r]=H);
	   H:=ReflectionSubgroup(W,W.rootInclusion{H});
      fi;
      Add(rel,rec(refs:=H.rootInclusion{H.generatingReflections},
	   WH:=Normalizer(H,Group(s.levi))));
    od;
    if s.spets.phi=() then WF:=W;L:=Group(s.levi);
    else WF:=Group(Concatenation(W.generators,[s.spets.phi]),W.identity);
      L:=Subgroup(WF,Group(s.levi).generators);
    fi;
    if Size(L)=1 then
      s.WGL:=Centralizer(W,s.levi.phi);
      for H in rel do H.WH:=Centralizer(H.WH,s.levi.phi);od;
      func:=x->x;
    else
      NLF:=Normalizer(WF,L);
      if NLF=L then 
        if s.levi<>s.spets then Error(s," N=L\n");return false;fi;
	s.WGL:=CoxeterGroup();s.WGL.reflists:=[];s.WGLdims:=[1];
	return s.WGL;
      fi;
      NLFb:=FactorGroup(NLF,L);
      hom:=NaturalHomomorphism(NLF,NLFb);
      phi:=Image(hom,s.levi.phi);
      s.WGLdegrees:=RelativeDegrees(s.spets,s.d);
      if Number(RelativeDegrees(s.levi,s.d),x->x<>1)=0 and 
         Number(s.WGLdegrees,x->x<>1)=1 
      then
	if  s.levi.phi in W and Order(NLFb,phi)=Product(s.WGLdegrees) then 
	  s.WGL:=Subgroup(NLFb,[phi]);
	else N:=Centralizer(Subgroup(WF,W.generators),s.levi.phi);
	  if Size(N)=Size(L)*Product(s.WGLdegrees)
	  then s.WGL:=FactorGroup(N,L);
	  fi;
	fi;
      fi;
      if not IsBound(s.WGL) then  N:=Subgroup(WF,W.generators);
	N:=Subgroup(N,Intersection(N,NLF).generators);
	s.WGL:=Centralizer(FactorGroup(N,L),phi);
      fi;
      if Length(rel)=1 and rel[1].WH=NLF then rel[1].WH:=s.WGL;
      else rel:=List(rel,function(x)local N,WH;
	N:=Intersection(s.WGL,FactorGroup(Subgroup(NLF,x.WH.generators),L));
	if IsGroup(N) then x.WH:=N;
	else x.WH:=Subgroup(s.WGL,N);
	fi;
	return x;end);
      fi;
      func:=x->x.element.representative;
    fi;
    ud:=CycPolUnipotentDegrees(s.levi);
    eig:=Eigenvalues(UnipotentCharacters(s.levi));
    ud:=Filtered([1..Length(ud)],
        i->ud[i]=ud[s.cuspidal] and eig[i]=eig[s.cuspidal]);
    if Length(ud)>1 then
      c:=Size(s.WGL);
      s.WGL:=Stabilizer(s.WGL,Position(ud,s.cuspidal),function(c,g)
	return c^PermutationOnUnipotents(s.levi,func(g),ud);end);
      if c<> Size(s.WGL) then
	ChevieErr("# WGL:",c,"/",Size(s.WGL)," fix c\n");
      fi;
      for H in rel do H.WH:=Intersection(H.WH,s.WGL);od;
      rel:=Filtered(rel,x->Size(x.WH)>1);
    fi;
    if not ForAll(rel,x->IsCyclic(x.WH)) then Error("a");fi;
    hom:=List(rel,x->First(Elements(x.WH),y->Order(x.WH,y)=Size(x.WH)));
    if Subgroup(s.WGL,hom)<>s.WGL then Error("b");fi;
    reflists:=List(rel,x->x.refs);
    m:=Concatenation(V,NullspaceMat(s.projector))^-1;
    rel:=List(hom,function(x)local M;M:=MatXPerm(W,func(x))^m;
      return M{[1..Length(V)]}{[1..Length(V)]};end);
    rel:=List(rel,AsReflection);
    for c in [1..Length(rel)] do
      r:=AsRootOfUnity(rel[c].eigenvalue);
      m:=Numerator(r);r:=Denominator(r);
      if m<>1 then
	rel[c]:=AsReflection(Reflection(rel[c].root,rel[c].coroot)^(1/m mod r));
      fi;
    od;
    s.WGL:=PermRootGroup(List(rel,x->x.root),List(rel,x->x.coroot));
    s.WGL.parentMap:=hom;
    s.WGL.reflists:=reflists{List(s.WGL.roots{s.WGL.generatingReflections},
      x->PositionProperty(rel,y->ProportionalityCoefficient(y.root,x)<>false))};
  fi;
  s.WGLdims:=List(CharTable(s.WGL).irreducibles,x->x[1]);
  if ForAll(s.WGLdims,x->x=1) then s.e:=Length(s.WGLdims);fi;
  return s.WGL;
end;
else
SeriesOps.RelativeGroup:=function(s)local id,W,L,N,ud,eig,c,V,refs,r,m,v,
  getreflection,reflist,func,hplane,smalltobig,bigtosmall,rrefs,NF,LF,NFQ;
  if IsBound(s.WGL) then return s.WGL;fi;
  id:=s.projector^0;
  W:=Group(s.spets);L:=Group(s.levi);
  if s.projector=id then #central series
    s.WGL:=Copy(W); s.WGL.parentMap:=s.WGL.generators;
    s.WGL.reflists:=List(W.rootInclusion{W.generatingReflections},x->[x]);
    s.WGLdims:=List(CharTable(s.WGL).irreducibles,x->x[1]);
    if IsCyclic(s.WGL) then s.e:=Size(s.WGL);fi;
    return s.WGL;
  fi;
  N:=Normalizer(W,L);
  if s.levi.phi<>() then 
    if Size(L)=1 then N:=Centralizer(N,s.levi.phi);
    elif IsCoxeterCoset(s.spets) then
      N:=Group(List(N.generators,x->ReducedInRightCoset(Group(s.levi),x)),());
      N:=Centralizer(N,s.levi.phi);
      N:=Subgroup(W,Concatenation(N.generators,Group(s.levi).generators));
    else #   N:=Stabilizer(N,s.levi); # is shorter but slower...
      NF:=Centralizer(N,s.levi.phi);
      if Size(NF)=Size(L)*Product(RelativeDegrees(s.spets,s.d)) then
        N:=NF;
      else
	NF:=Group(Concatenation(N.generators,[s.levi.phi]),());
	LF:=Subgroup(NF,L.generators);
	NFQ:=NF/LF;
	N:=Subgroup(NFQ,List(N.generators,x->x^NaturalHomomorphism(NF,NFQ)));
	NF:=Centralizer(NFQ,s.levi.phi^NaturalHomomorphism(NF,NFQ));
	NFQ:=Intersection(NF,N);
	N:=Subgroup(W,Concatenation(L.generators,List(NFQ.generators,x->
	  x.element.representative)));
      fi;
    fi;
  fi;
  eig:=Eigenvalues(UnipotentCharacters(s.levi));
  eig:=Filtered([1..Length(eig)],i->eig[i]=eig[s.cuspidal]);
  if Length(eig)>1 then
    ud:=CycPolUnipotentDegrees(s.levi);
    ud:=Filtered([1..Length(eig)],i->ud[i]=ud[s.cuspidal]);
    if Length(ud)>1 then
      c:=Size(N);
      N:=Stabilizer(N,Position(ud,s.cuspidal),function(c,g)
	return c^PermutationOnUnipotents(s.levi,g,ud);end);
      if c<> Size(N) then
	ChevieErr("# WGL:",Size(N)/Size(L),"/",c/Size(L)," fix c\n");
      fi;
    fi;
  fi;
  refs:=Set(List(Reflections(W),x->Position(W.reflections,x)));
  if Size(L)=1 then s.WGL:=N;func:=x->x;
  elif Size(N)=Size(L) then 
    if s.levi<>s.spets then Error(s," N=L\n");return false;fi;
    s.WGL:=CoxeterGroup();s.WGL.reflists:=[];s.WGLdims:=[1];
    return s.WGL;
  else s.WGL:=N/L; func:=x->x.element.representative;
  fi;
  V:=NullspaceMat(s.projector-id);  # The E(d)-eigenspace
  m:=Concatenation(V,NullspaceMat(s.projector))^-1;
  hplane:=r->SumIntersectionMat(NullspaceMat(MatXPerm(W,Reflection(W,r))-id),
     V)[2];
  smalltobig:=h->List(h,x->Concatenation(x,[1..W.rank-Length(V)]*0))*m^-1;
  bigtosmall:=function(x)x:=x^m;return x{[1..Length(V)]}{[1..Length(V)]};end;
  getreflection:=function(rr)local rH,H,r,res,n;
    rH:=hplane(rr[1]); # Print("rH=",rH,"\n"); # Hyperplane of W_G(L)
    if Length(rH)=Length(V) then return false;fi;
    if Length(rH)=0 then # cyclic WGL
      ConjugacyClasses(L); # hack to fix a (bug in GAP?)
      H:=W;
    else H:=ReflectionSubgroup(W,Concatenation(W.rootInclusion{rr},
                              L.rootInclusion{L.generatingReflections}));
    fi;
    res:=rec(refs:=H.rootInclusion{H.generatingReflections});
    H:=Intersection(H,N);
    if Size(L)<>1 then H:=H/L;fi;
    if Size(H)=1 then return false;fi;
    if not IsCyclic(H) then Error("a");fi;
    res.hom:=First(Elements(H),y->Order(H,y)=Size(H));
    r:=bigtosmall(MatXPerm(W,func(res.hom)));
    Inherit(res,AsReflection(r));
    n:=AsRootOfUnity(res.eigenvalue);n:=1/Numerator(n) mod Denominator(n);
    r:=r^n;
    Inherit(res,AsReflection(r)); res.WH:=H; return res;
  end;
  v:=function(h)if Length(h)=0 then return VectorSpace([[0]],Cyclotomics);
    else return VectorSpace(h,Cyclotomics);fi;end;
  rrefs:=CollectBy(refs,x->v(hplane(x)));
  rrefs:=Filtered(rrefs,x->not W.rootInclusion[x[1]] in L.rootInclusion);
  reflist:=[];
  for r in rrefs do
    Add(reflist,getreflection(r));
    if Subgroup(s.WGL,List(reflist,x->x.hom))=s.WGL then 
      s.WGL:=PermRootGroup(List(reflist,x->x.root),List(reflist,x->x.coroot));
      reflist:=List(s.WGL.matgens,x->smalltobig(NullspaceMat(x-x^0)));
      reflist:=List(reflist,h->First(rrefs,rr->v(hplane(rr[1]))=v(h)));
      s.WGL.reflists:=List(reflist,getreflection);
      s.WGL.parentMap:=List(s.WGL.reflists,x->x.hom);
      s.WGL.reflists:=List(s.WGL.reflists,x->x.refs);
      s.WGLdims:=List(CharTable(s.WGL).irreducibles,x->x[1]);
      if IsCyclic(s.WGL) then s.e:=Size(s.WGL);fi;
      return s.WGL;
    fi;
  od;
  Error("b");
end;
fi;

# takes a d-series s with s.spets split; fills in s.charNumbers, s.eps, s.dims
SeriesOps.CharNumbers:=function(s)local ad,q,ud,v,eig,f,c,t,cand,i,check,Ed;
  if IsBound(s.charNumbers) then return s.charNumbers;fi;
  if not IsBound(s.WGL) and RelativeGroup(s)=false then return false;fi;
  ud:=CycPolUnipotentDegrees(s.spets);
  ad:=Number(RelativeDegrees(s.levi,s.d),x->x<>1);
  cand:=Filtered([1..Length(ud)],i->ad=Valuation(ud[i],s.d));
  s.RLG:=LusztigInduction(s.spets,UnipotentCharacter(s.levi,s.cuspidal));
  if s.RLG=false and s.levi.phi=() then 
   s.RLG:=HarishChandraInduction(s.spets,UnipotentCharacter(s.levi,s.cuspidal));
  fi;
  if s.RLG=false then ChevieErr(s,":RLG failed\n");
  elif s.degree<>CycPol(Degree(s.RLG))then 
    ChevieErr(s,":Deg RLG<>Sum(eps[i]*ud[i])\n");
  fi;
# now use that   S_\phi(q)=\eps_\phi Deg(RLG(\lambda))/Deg(\gamma_\phi)
  cand:=List(cand,c->rec(charNumbers:=c,sch:=s.degree/ud[c]));
  cand:=Filtered(cand,c->ForAll(c.sch.vcyc,x->x[2]>0)); # positive
# Print("after deg cands=",List(cand,c->c.charNumbers),"\n");
  Ed:=E(Denominator(s.d))^Numerator(s.d);
  q:=Mvp("q");
  ad:=CycPol(q-Ed)^ad;
  f:=s.degree/ad;
  if ForAny(f.vcyc,x->x[2]<0) then
    ChevieErr(s," cuspidal is not\n");return false;
  fi;
  v:=Value(f,Ed);
  for c in cand do 
    c.span:=Degree(c.sch)-c.sch.valuation;
# now use that 
    f:=Sum(s.WGLdims,x->x^2)*Value(ud[c.charNumbers]/ad,Ed)/v;
    c.dims:=AbsInt(f);c.eps:=SignInt(f);
  od;
  cand:=Filtered(cand,c->c.dims in s.WGLdims);
  eig:=Eigenvalues(UnipotentCharacters(s.levi),[s.cuspidal])[1];
  eig:=eig*List([1..Denominator(s.d)^2],i->E(Denominator(s.d)^2)^i);
  cand:=Filtered(cand,c->Eigenvalues(UnipotentCharacters(s.spets),
      [c.charNumbers])[1] in eig);
# Print("after eig cands=",List(cand,c->c.charNumbers),"\n");
  if Length(cand)<Length(s.WGLdims) then 
    ChevieErr(s,":not enough left with predicted eigenvalues in ",
        FormatGAP(List(eig,AsRootOfUnity)),"\n");return false;
  fi;
  SortBy(cand,x->x.dims);Sort(s.WGLdims);
  f:=function(arg)  # f("field"[,indices])
    if Length(arg)=1 then return List(cand,x->x.(arg[1]));fi;
    if IsList(arg[2]) then return List(cand{arg[2]},x->x.(arg[1]));
                     else return cand[arg[2]].(arg[1]);fi;
  end;
  check:=function()local n;
    for n in ["charNumbers","eps","dims","span"] do s.(n):=f(n);od;
    if s.RLG<>false and (Filtered([1..Length(s.RLG.v)],
       i->s.RLG.v[i]<>0)<>Set(s.charNumbers)
       or s.RLG.v{s.charNumbers}<>Zip(s.dims,s.eps,function(x,y)return x*y;end))
    then ChevieErr(s,":RLG does not match");
    fi;
    return s.charNumbers;
  end;
  if Length(cand)=Length(s.WGLdims) then return check();fi;
  ud:=Zip(CycPolUnipotentDegrees(s.spets){f("charNumbers")},f("eps"),f("dims"),
     function(x,y,z)return x*y*z;end);
  t:=Maximum(List(ud,Degree));
  c:=p->Concatenation(Coefficients(Value(p,q),"q"),[1..t-Degree(p)]*0);
  v:=SubsetsSum(c(s.degree),List(ud,c),s.WGLdims,f("dims"));
  InfoChevie("# ",Length(v)," found\n");
  if Length(v)>10000 then
    InfoChevie("# ",Length(v)," combinations sum to dimRLG\n");
  elif Length(v)=0 then
    ChevieErr(s," no combination sums to dimRLG\n");return false;
  fi;
  
  if IsBound(s.e) then
    s.charNumbers:=f("charNumbers"); s.dims:=f("dims"); SeriesOps.mC(s);
    if Length(v)>1 then
      v:=Filtered(v,a->f("span",a)=
        List(a,i->Sum(Difference(a,[i]),j->AbsInt(s.mC[i]-s.mC[j]))));
      if Length(v)>10000 then
        InfoChevie("# ",Length(v)," combinations have right span\n");
      fi;
    elif Length(v)=0 then
      ChevieErr(s," no combination has right span\n");return false;
    fi;
    Unbind(s.charNumbers);Unbind(s.dims);Unbind(s.mC);
  fi;

  if Length(v)>1 then
    InfoChevie("# after span ",Length(v)," combinations\n");
    InfoChevie("# Warning: using Mackey with tori for ",s,"\n");
    i:=FusionConjugacyClasses(s.levi,s.spets);
    c:=Zip(CharTable(s.spets).centralizers{i},CharTable(s.levi).centralizers,
      function(a,b)return a/b;end);
    t:=Zip(TransposedMat(DeligneLusztigCharacterTable(s.levi))[s.cuspidal],c,
      function(a,b)return a*b;end);
    t:=List([1..Length(ConjugacyClasses(s.spets))],
            k->Sum(t{Filtered([1..Length(i)],j->i[j]=k)}));
    c:=TransposedMat(DeligneLusztigCharacterTable(s.spets));
    v:=Filtered(v,a->Zip(f("dims",a),f("eps",a),
           function(x,y)return x*y;end)*c{f("charNumbers",a)}=t);
  fi;

#  Print(" after Mackey with tori ",Length(v)," combinations\n");
  if Length(v)>1 then
     ChevieErr(s," ",
     Join(List(FactorsSet(List(v,x->f("charNumbers",x))),x->FormatGAP(x)),"x"),
       " chars candidates: using RLG\n");
     if s.RLG=false then return false; fi;
     v:=Filtered(v,l->ForAll(l,i->s.RLG.v[f("charNumbers",i)]<>0));
     if Length(v)<>1 then Error();fi;
  elif Length(v)=0 then ChevieErr(s," no candidates left\n");return false;
  fi;
  cand:=cand{v[1]};return check();
end;

COMPACTCOHOMOLOGY:=true;
SeriesOps.fill:=function(s)local uc,Schur,r,p,ratio,LFrob,predeigen,map,
   series,noncus,unique,i,m,quality,rr,FractionToRoot,a,param,j,o,u,nid;
  if not IsBound(s.e) then Error("fill assumes .e bound\n");fi;
  if not IsBound(s.charNumbers) and CharNumbers(s)=false then return false;fi;
  if not IsBound(s.mC) then SeriesOps.mC(s);fi;
  uc:=UnipotentCharacters(s.spets);
  Schur:=CycPolUnipotentDegrees(s.spets){s.charNumbers};
  Schur:=List([1..s.e],x->s.degree/Schur[x]*s.eps[x]);
  s.eigen:=Eigenvalues(uc,s.charNumbers);
  LFrob:=AsRootOfUnity(Eigenvalues(UnipotentCharacters(s.levi))[s.cuspidal]);
  m:=ReflectionDegrees(Group(s.spets));
  s.delta:=Lcm(List(Filtered(ReflectionDegrees(s.spets),x->x[1]<>1),
        x->Denominator(AsRootOfUnity(x[2]))));
  rr:=function(j,i)return (i-1)/s.e-s.mC[j]*s.d;end;
  FractionToRoot:=x->E(Denominator(x))^Numerator(x);
  param:=function(j,i)return Mvp("q")^s.mC[j]*FractionToRoot(rr(j,i));end;
# parameters of Hecke algebra are List([1..s.e],i->param(i,i))
  predeigen:=function(j,i) #eigenvalue for mC[j] and \zeta_e^{i-1}
    return FractionToRoot(s.delta*rr(j,i)*s.e*s.d+LFrob);
  end;
  # as fraction predeigen_i:=delta di -delta m_i d^2e 
  if COMPACTCOHOMOLOGY then map:=FitParameter(Schur,s.mC,s.d);
  else map:=FitParameter(Schur,-s.mC,s.d);
  fi;
  if map=false then ChevieErr(s," FitParameter failed\n");return false;fi;
  # possible perms encoded by map[i][1]->map[i][2]
  nid:=uc.almostHarishChandra[1].charNumbers[PositionId(s.spets)];
  if nid in s.charNumbers then # id should correspond to id(WGL)
     predeigen:=function(j,i) #eigenvalue for mC[j] and \zeta_e^{i-1}
       return FractionToRoot(s.d*s.e*s.delta*(rr(j,i)
                                +s.d*s.mC[Position(s.charNumbers,nid)]));
     end;
  fi; 
  series:=List(s.charNumbers,
     x->PositionProperty(uc.harishChandra,y->x in y.charNumbers));
# unique:=Filtered([1..s.e],
#      i->Length(uc.harishChandra[series[i]].charNumbers)>1);
  # non-cuspidal thus assumed known eigenvalue
# if Length(unique)=0 then unique:=[1..s.e];fi;
# unique:=Filtered(map,p->Length(p[1])=1 and p[1][1] in unique);
  unique:=Filtered(map,p->Length(p[1])=1);
  ratio:=List(unique,p->s.eigen[p[1][1]]/predeigen(p[1][1],p[2][1]));
# Print("ratio=",ratio,"\n");
  if Length(Set(ratio))>1 then 
    ChevieErr(s," eigenvalue ratios=",ratio,"\n");return false;fi;
  ratio:=AsRootOfUnity(ratio[1]);
# now find integer translation t such that mod 1. we have t delta d=ratio
  ratio:=ratio*Denominator(s.d*s.delta);
  if not IsInt(ratio) then 
    ChevieErr(s,"non-integral ratio=",ratio,"\n");return false;fi;
  if ratio=0 then r:=0;
  else r:=QuotientMod(ratio,Numerator(s.d*s.delta),Denominator(s.d*s.delta));
  fi;
  if r=false then Error(); fi;# should not happen
  map:=List(map,x->[x[1],List(x[2]+r-1,y->(y mod s.e)+1)]);
  r:=[];
  s.permutable:=List(s.charNumbers,x->false);j:=1;
  for i in map do
    a:=Arrangements(i[2],Length(i[2]));
    p:=PositionsProperty(a,A->s.eigen{i[1]}=
      List(A,j->predeigen(i[1][1],j)));
    if p=[] then
      ChevieErr(s,"predicted eigenvalues cannot match actual\n");return false;
    else 
      if Length(p)>1 then 
        o:=Orbits(ApplyFunc(Group,List(a{p},x->PermListList(a[p[1]],x))),
	  [1..Length(i[2])]);
        if Length(p)<>Product(o,x->Factorial(Length(x))) then Error();fi;
	for u in o do s.permutable{i[1]{u}}:=j+u*0;j:=j+1;od;
      fi;
      r{i[1]}:=a[p[1]];
    fi;
  od;
  p:=SortingPerm(r);
  for i in ["mC","charNumbers","eigen","span","eps","dims","permutable"] do 
    s.(i):=Permuted(s.(i),p);
  od;
  s.translation:=Filtered([0..s.e-1],t->s.eigen=
        List([1..s.e],i->predeigen(i,(i+t)mod s.e)));
  if Length(s.translation)=1 then Unbind(s.translation);
  else p:=s.translation[2];
    if ForAny(s.translation,x->x mod p<>0) or Length(s.translation)<>s.e/p
    then Error();
    fi;
#   Error();
    quality:=List(s.translation,t->List([1..s.e],i->param((t+i-1)mod s.e+1,i)));
    quality:=List(quality,x->NF(Set(Product(x,y->Mvp("x")-y).coeff)));
    quality:=1+ListBlist(s.translation,
      List(quality,x->ForAll(quality,j->IsSubset(j,x))));
    if quality=[] then quality:=[1];fi;
    m:=Rotations(s.mC);
    m:=Position(m,Maximum(m{quality}));
    m:=Rotation([1..s.e],m-1);
    for i in ["mC","charNumbers","eigen","span","eps","dims","permutable"] do 
      s.(i):=s.(i){m};od;
    s.translation:=p;
  fi;
  s.Hecke:=Hecke(s.WGL,[List([1..s.e],i->param(i,i))]);
  return s.Hecke;
end;

SeriesOps.Hecke:=function(s)local H;
  if IsBound(s.Hecke) then return s.Hecke;fi;
  CharNumbers(s);
  if IsBound(s.e) then SeriesOps.fill(s);
  elif CHEVIE.relativeSeries then 
    InfoChevie("\n      # Relative: ",String(s));
    s.relativeSeries:=SeriesOps.RelativeSeries(s);
  fi;
  if IsBound(s.Hecke) then return s.Hecke;else return false;fi;
end;

AllProperSeries:=function(W)local l;
  l:=Set(List(Flat(ReflectionEigenvalues(W)),Denominator));
  return Concatenation(List(l,d->Concatenation(List(
      [1..Length(RelativeDegrees(W,d))],i->CuspidalSeries(W,d,i)))));
end;


findfractionalpowers:=function(W)local uc,h,L,cusp,sers,s,p,fix,i,n,
  reasons,frac,UNBOUND;
  if W.nbGeneratingReflections=0 then return [0];fi;
  uc:=UnipotentCharacters(W);
  UNBOUND:=99;uc.fractions:=[1..Size(uc)]*0+UNBOUND;
  reasons:=List([1..Size(uc)],x->[]);
  for h in uc.harishChandra do
    if Length(h.levi)<>W.nbGeneratingReflections then
      L:=ReflectionSubgroup(W,h.levi);
      cusp:=FindCuspidalInLevi(h.cuspidalName,L);
      n:=findfractionalpowers(L);
      if n[cusp]<>UNBOUND then
        uc.fractions{h.charNumbers}:=h.charNumbers*0+n[cusp];
	for i in h.charNumbers do Add(reasons[i],IntListToString(h.levi));od;
      fi;
    fi;
  od;
  sers:=AllProperSeries(W);
  List(sers,Hecke);
  sers:=Set(Filtered(sers,x->IsBound(x.mC) and IsBound(x.e)));
  while Length(sers)<>0 do
  for s in sers do
    p:=PositionProperty(s.charNumbers,i->IsBound(uc.fractions[i]));
    if p=false then Print("reticent series",s,"\n");
    else
      frac:=s.mC[p]*s.e*s.d;
      fix:=Mod1(uc.fractions[s.charNumbers[p]]-frac);
      if fix<>0 then
        Print(" **** Badly normalized series ", s," adding ",fix,"\n");
      fi;
      Print("\n",s,"=>");
      for i in [1..s.e] do
        Print(s.charNumbers[i],".\c");
        frac:=Mod1(s.mC[i]*s.e*s.d+fix);
        if uc.fractions[s.charNumbers[i]]=UNBOUND then
	  uc.fractions[s.charNumbers[i]]:=frac;
	  Add(reasons[s.charNumbers[i]],s.d);
        elif uc.fractions[s.charNumbers[i]]<>frac then
	  Print("Failed! ",s.charNumbers[i],"=",
	    TeXStrip(uc.TeXCharNames[s.charNumbers[i]])," in ",s,
	    "\n     where mC mod1=",frac,"\n conflicts with ",
	     reasons[s.charNumbers[i]],"\n");
	else Add(reasons[s.charNumbers[i]],s.d);
	fi;
      od;
      sers:=Difference(sers,[s]);
    fi;
  od;
  od;
  for i in uc.harishChandra do 
   if Length(i.charNumbers)>1 then
     Print(Format([uc.fractions{i.charNumbers},
                                 List(reasons{i.charNumbers},Length)]),"\n");
   else Print(uc.fractions[i.charNumbers[1]],reasons[i.charNumbers[1]],"\n");
   fi;
   p:=Set(uc.fractions{i.charNumbers});
   if IsBound(i.qEigen) then
     if p<>[i.qEigen] then 
       if p=[UNBOUND] then
       Print("!!!!!!!!!!!! for ",i.charNumbers," qEigen should be UNBOUND is ",
         i.qEigen,"\n");
       else
       Error("for HCseries ",i," of ",ReflectionName(W)," qEigen should be ",p);
       uc.fractions{i.charNumbers}:=i.qEigen+0*i.charNumbers;
       fi;
     fi;
   else 
     if p<>[0] then Print("!!!!!!!!!!!!!! qEigen unbound should be ",p,"\n");fi;
   fi;
  od;
  return uc.fractions;
end;

# series of d-regular element w
PrincipalSeries:=function(W,d)local s;
  s:=SplitLevis(W,d,Length(RelativeDegrees(W,d)));
  if Length(s)<>1 then Error("not one ",d,"-Sylow");fi;
  s:=s[1];
  return Series(W,s,Position(UnipotentCharacters(s).A,0),d);
end;

CheckMaximalPairs:=function(W)local divs,l,s,d,c,e,res;
  divs:=Set(List(Flat(ReflectionEigenvalues(W)),Denominator));
  res:=List(divs,d->PrincipalSeries(W,d));
  for s in res do
    e:=Elements(s.levi);
    ForAny(e,function(x)local Z;
      Z:=Centralizer(Group(s.spets),x);
      Z:=Group(Concatenation(Z.generators,Group(s.levi).generators),());
      if Size(Z)<>Size(RelativeGroup(s))*Size(s.levi) then Error();fi;
      return true;
      end);
  od;
  return res;
end;

ReflectionCosetFromType:=function(arg)local g,t,g1,i;
  if Length(arg)=0 then return CoxeterGroup();fi;
  if Length(arg)>1 then Error(ReflectionName(t)," non-irred not yet");fi;
  for t in arg do
    g:=ReflectionGroup(t.orbit[1]);
    if Length(t.orbit)>1 then 
      g1:=g;
      for i in [1..Length(t.orbit)-1] do g:=g*g1;od;
      i:=MatXPerm(g,t.twist*PermList(
        Concatenation(Rotations(List(g.type,x->x.indices))[2])));
      g:=Spets(g,i);
    else g:=Spets(g,t.twist);
    fi;
  od;
  return g;
end;

SpetsFromType:=function(t)local W;
  W:=ApplyFunc(ReflectionGroup,t.orbit);
  return Spets(W,t.twist);
end;

# get Hecke of series s. If scalar-twisted untwist first
getHecke:=function(s)local t,scal,g,l,s1,p,e,c;
  t:=ReflectionType(s.spets);
  if not IsBound(s.charNumbers) then CharNumbers(s);fi;
 if Length(t)=1 and IsBound(t[1].scalar) and not ForAll(t[1].scalar,x->x=1) then
    t:=t[1];
    scal:=t.scalar;t:=ShallowCopy(t);
    InfoChevie("   # removing scal=",scal,"\n");
    t.scalar:=List(t.scalar,x->1);
    if t.orbit[1].series="B" and t.orbit[1].rank=2 and t.twist=(1,2) then
      t.orbit[1].cartanType:=E(8)-E(8)^3;
    fi;
    g:=SpetsFromType(t);
    if g=false then return false;fi;
    e:=WordEnumerator(Group(s.spets));
    p:=EltWord(Group(g),e.Get(s.levi.phi/s.spets.phi));
#   p:=OnTuples(Group(s.spets).rootInclusion,s.levi.phi/s.spets.phi);
#   p:=PermListList(Group(s.spets).rootInclusion,p);
#   if p<>() then Error();fi;
    Reflections(Group(g));
    l:=SubSpets(g,List(Group(s.levi).generators,
      x->Position(Group(g).reflections,EltWord(Group(g),e.Get(x)))),p);
    if Length(scal)>1 then ChevieErr("scal=",scal," unimplemented\n");
      return false;
    fi;
    scal:=AsRootOfUnity(scal[1]);
    # one should do an Ennola of the Levi's cuspidal by the absorbed part
    # of scal by the center of the levi
    if not s.cuspidal in CuspidalUnipotentCharacters(l,Mod1(s.d-scal)) then
      e:=Ennola(Group(l))[1];c:=AbsInt(e.ls[s.cuspidal]);
      if not c in CuspidalUnipotentCharacters(l,Mod1(s.d-scal)) then
        c:=AbsInt(LsPuiss(e.ls,-1)[s.cuspidal]);
      fi;
    else c:=s.cuspidal;
    fi;
    s1:=Series(g,l,c,Mod1(s.d-scal));
    p:=getHecke(s1);
    if p<>false then 
      p:=Value(p,["q",E(Denominator(scal))^-Numerator(scal)*Mvp("q")]);
    fi;
    return p;
  else
    if ForAny(t,u->IsBound(u.scalar) and not ForAll(u.scalar,x->x=1)) then
      ChevieErr("scals=",List(t,u->u.scalar)," unimplemented\n");
    fi;
#   CHECK.title:="getHecke";
    SeriesOps.fill(s);
    if IsBound(s.Hecke) then return s.Hecke.parameter[1];
    else return false;fi;
  fi;
end;

SeriesOps.RelativeSeries:=function(s)local res,p,ud,u1,f,o,aA,tt;
  if not IsBound(s.charNumbers) then CharNumbers(s);fi;
  res:=List(s.WGL.reflists,function(r)local R,w,i,l;
    i:=Position(s.WGL.reflists,r);
    R:=SubSpets(s.spets,r,s.element/s.spets.phi);
    l:=Group(s.levi).rootInclusion;
    if not IsSubset(Group(R).rootInclusion,l) then
      r:=List(r,function(w)local p,rr;rr:=Group(s.spets).roots;
        p:=PositionProperty(l,
          y->ProportionalityCoefficient(rr[w],rr[y])<>false); 
        if p<>false then return l[p];fi;
        return w;
        end);
      R:=SubSpets(s.spets,r,s.element/s.spets.phi);
      if not IsSubset(Group(R).rootInclusion,r) then
        l:=Set(List(l,function(w)local p,rr;rr:=Group(s.spets).roots;
          p:=PositionProperty(Group(R).rootInclusion,
            y->ProportionalityCoefficient(rr[w],rr[y])<>false); 
          if p<>false then return Group(R).rootInclusion[p];fi;
          return w;
          end));
        if not IsSubset(Group(R).rootInclusion,l) then
         Error("could not change r=",r," so that Subspets(r) contains",
           l,"\n");
        fi;
      fi;
    fi;
    p:=Series(R,SubSpets(R,l,s.element/R.phi),s.cuspidal,s.d);
    p.orbit:=s.WGL.orbitRepresentative[i];
    r:=List(Set(s.WGL.reflections),x->Position(s.WGL.reflections,x));
    if s.WGL.generators=s.WGL.parentMap then
      r:=Filtered(r,x->s.WGL.reflections[x] in Group(p.spets));
    else
      w:=List(r,x->GetWord(s.WGL,s.WGL.reflections[x]));
      w:=List(w,x->Product(s.WGL.parentMap{x}));
      r:=ListBlist(r,List(w,
	function(x)local g;g:=Group(p.spets);
	  if IsRec(x) then return Representative(x.element) in g;
	  else return x in g;fi;end));
    fi;
    p.WGL:=ReflectionSubgroup(s.WGL,r);
    p.WGLdims:=List(CharTable(p.WGL).irreducibles,x->x[1]);
    if ForAll(p.WGLdims,x->x=1) then p.e:=Length(p.WGLdims);fi;
    return p;end);
  s.relativeSpets:=List(res,x->x.spets);
# tt:=CHECK.title;
  p:=List(res,getHecke);
# CHECK.title:=tt;
  if false in p then return res;fi;
  s.Hecke:=Hecke(s.WGL,p);
  u1:=List(SchurElements(s.Hecke),Mvp);
  if ForAny(u1,x->ForAny(x.elm,y->not ForAll(y.coeff,IsInt))) then 
     ChevieErr(s.Hecke," wrong set of SchurElements");return res;fi;
  u1:=List(u1,x->s.degree/CycPol(Mvp(x)));
  ud:=List(CycPolUnipotentDegrees(s.spets){s.charNumbers},
   x->x*SignInt(Value(x/CycPolUnipotentDegrees(s.levi)[s.cuspidal],
      E(Denominator(s.d))^Numerator(s.d))));
  p:=PermListList(u1,ud);
  # the permutation should also take in account eigenvalues
  if p=false then ChevieErr(s.Hecke," wrong set of SchurElements");return res;fi;
  s.charNumbers:=Permuted(s.charNumbers,p);
  if IsBound(s.span) then s.span:=Permuted(s.span,p);fi;
  aA:=List(u1,x->E(Denominator(s.d)^2)^(x.valuation+Degree(x)));
  p:=PositionRegularClass(s.WGL,s.d);
  if p=false and Size(s.WGL)=1 then p:=1;fi;
  o:=List(CharTable(s.WGL).irreducibles,x->x[p]/x[1]);
  s.predictedEigen:=List([1..Length(o)],i->aA[i]*o[i])*
    Eigenvalues(UnipotentCharacters(s.levi))[s.cuspidal];
  return res;
end;

Checkzegen:=function(W)local l,uc,aA,e,z,s,ucL,aAL;
  l:=AllProperSeries(W);
# for s in l do SeriesOps.RelativeSeries(s);od;
  Print("with  no Hecke:",Filtered(l,s->Hecke(s)=false),"\n");
  Print("    center=",OrderCenter(W),"\n");
  l:=Filtered(l,s->IsBound(s.Hecke));
  uc:=UnipotentCharacters(W);
  aA:=uc.a+uc.A;
  for s in l do
    e:=gete(s.Hecke);
    z:=OrderCenter(s.WGL);
    ucL:=UnipotentCharacters(s.levi);
    aAL:=ucL.A+ucL.a;aAL:=aAL[s.cuspidal];
    Print(s,"\n");
    Print(FormatTable(
       [List((aA{s.charNumbers}-aAL)/z,Mod1),List(aA{s.charNumbers}/z,Mod1),
	  List(e,x->Mod1(1/x))],
      rec(rowLabels:=["(aA-aAL)/z","aA/z","local Hecke irr"])));
  od;
end;

# checks minimal polynomials for Cyclotomic hecke algebras
CheckRatCyc:=function(s)local l,fact,F;
  if not IsBound(s.Hecke) then 
    return SPrint(FormatTeX(s)," Hecke not bound");fi;
  l:=CollectBy(s.Hecke.parameter[1],Degree);
  fact:=List(l,x->Product(x,y->Mvp("T")-y));
  F:=FieldOfDefinition(Group(s.spets));
# Append(fact,[ForAll(fact,x->ForAll(x.coeff,y->y in F)),
#              ForAll(fact,x->ForAll(x.elm,y->ForAll(y.coeff,IsInt)))]);
  return SPrint(FormatTeX(s)," $",
    Join(List(fact,x->SPrint("(",FormatTeX(x),")")),""),"$");
end;

CheckRatSer:=function(arg)local W,s,c,b;
  W:=ApplyFunc(ComplexReflectionGroup,arg);
  s:=AllProperSeries(W);
  s:=Filtered(s,x->x.spets<>x.levi);
  c:=Join(List(s,CheckRatCyc),"\\hfill\\break\n");
  return c;
# b:=Filtered([1..Length(s)],i->not IsString(c[i])
#   and (not c[i][Length(c[i])] or not c[i][Length(c[i])-1]));
# return TransposedMat([s{b},c{b}]);
end;

CheckPiGPiL:=function(n)local W,s,D0;
  W:=ComplexReflectionGroup(n);
  s:=AllProperSeries(W);List(s,Hecke);
  s:=Filtered(s,ser->IsBound(ser.e) and ser.principal);
  D0:=s->
    Sum(ReflectionDegrees(Group(s.spets))+ReflectionCoDegrees(Group(s.spets)))-
    Sum(ReflectionDegrees(Group(s.levi))+ReflectionCoDegrees(Group(s.levi)));
  return List(s,x->D0(x)/x.e);
end;

Checkdovere:=function(n)local W,s;
  W:=ComplexReflectionGroup(n);
  s:=AllProperSeries(W);List(s,Hecke);
  s:=Filtered(s,ser->IsBound(ser.e) and ser.principal);
  return Filtered(s,x->x.e/x.d>2);
end;

CheckxiL:=function(n)local W,s;
  W:=ComplexReflectionGroup(n);
  s:=AllProperSeries(W);List(s,Hecke);
  s:=Filtered(s,ser->IsBound(ser.e) and ser.principal);
  return Filtered(s,x->AsRootOfUnity(PhiOnDiscriminant(x.levi))*x.d<>0);
end;

# check formula for product parameters
CheckCoN:=function(i) local W,s;
  W:=ComplexReflectionGroup(i);
  s:=AllProperSeries(W);List(s,Hecke);
  s:=Filtered(s,x->Size(x.levi)=1 and IsBound(x.Hecke));
  return List(s,function(ser)local m,sg;
    m:=DeterminantMat(MatXPerm(W,ser.element));
    sg:=(-1)^Sum(ReflectionDegrees(ser.WGL)-1);
    if Length(HyperplaneOrbits(ser.WGL))>1 then return false;
    else return m*sg*
    Product(ser.Hecke.parameter[1])^HyperplaneOrbits(ser.WGL)[1].N_s;
    fi;end);
end;

CheckLuCox:=function(s)local p,p1;
  p:=Product(Filtered(s.Hecke.parameter[1],x->x<>1),x->1-x);
  p1:=Value(s.degree,Mvp("q"));
#  return CycPol(p);
  return [p,p1];
end;

# data for linear char:
# returns d, e_W/d (mod d), \omega_c (mod d)
datafor:=function(W)local l,e;
  l:=List(RegularEigenvalues(W),d->PrincipalSeries(W,d));
  e:=Sum(ReflectionDegrees(W)+ReflectionCoDegrees(W));
  return List(l,s->
    [s.d,Gcd(e*s.d,Denominator(s.d))/
    Gcd(Gcd(List(HyperplaneOrbits(s.WGL),x->x.N_s)),Denominator(s.d))]);
end;

# description of centralizers of regular elements
cent:=function(W)local l,s,res,e;
  l:=List(RegularEigenvalues(W),d->PrincipalSeries(W,d));
  res:=List(l,s->[s.d,ReflectionName(s.WGL)]);
  l:=Set(List(res,x->x[2]));
  res:=List(l,x->[List(Filtered(res,y->y[2]=x),z->z[1]),x]);
  Sort(res);
  return List(res,x->SPrint(Join(x[1]),":",x[2]));
end;

EigenspaceNumbers:=function(W)local d,e;
  d:=ReflectionDegrees(W);
  if not IsSpets(W) then return Union(List(d,DivisorsInt));fi;
  e:=List(d,p->Lcm(p[1],Denominator(AsRootOfUnity(p[2]))));
  e:=Union(List(e,DivisorsInt));
  return Filtered(e,x->ForAny(d,p->E(x)^p[1]=p[2]));
end;

# nrconjecture: si kd est le multiple maximum de d tq
# RelativeDegrees(W,kd)<>[], alors
# Length(RelativeDegrees(W,d))=Length(RelativeDegrees(W,kd))*Phi(kd)/Phi(d)
nrconjecture:=function(W)local d,e,p,r,f;
  e:=Difference(EigenspaceNumbers(W),RegularEigenvalues(W));
  e:=Filtered(e,x->not x[1] in r);
  for d in e do
    f:=Filtered(r,x->x mod d=0);
    if Length(f)>0 then
      Print("**** failed: ",ReflectionName(W),d,f,"\n");
    fi;
    p:=First(Reversed(e),i->i[1] mod d=0);
    if p<>d then
      if Length(RelativeDegrees(W,d))*Phi(d)<>p[2]*Phi(p[1]) then 
        Print("**** failed: ",ReflectionName(W),p,d,"\n");
      fi;
    fi;
  od;
end;

# An element which has a maximal \zeta_d-eigenspace lives in the minimal
# parabolic subgroup such that the d-rank is the same
minimaldparabolic:=function(W,d)local l;
  l:=List(W.generatingReflections,i->Difference(W.generatingReflections,[i]));
  l:=Filtered(l,x->Length(RelativeDegrees(ReflectionSubgroup(W,x),d))=
                   Length(RelativeDegrees(W,d)));
  l:=List(l,x->ReflectionName(ReflectionSubgroup(W,x)));
  return l;
end;
