#############################################################################
##
#A  series.g            CHEVIE library             Jean Michel 
##
#Y  Copyright (C) 2003--2020   University Paris VII.
##
# This file contains code to work with d-Harish-Chandra series of Spetses
#------------------------- Utility functions ---------------------------

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
# res:=List([1..Length(v)],x->[]);
# for p in [1..e] do Add(res[G[p]],p);od;
  res:=CollectBy([1..Length(G)],G);
# pp:=[];for p in [1..Length(res)] do pp{res[p]}:=v[p];od;
  pp:=Concatenation(v);SortParallel(Concatenation(res),pp);
  para:=List([1..e],i->E(e)^pp[i]*X(Cyclotomics)^(den*m[i]));
  if List([1..e],k->Product(Filtered([1..e],j->j<>k),
   j->CycPol(1-para[k]/para[j])))<>sch then 
    ChevieErr("schur elms don't match\n");return false;fi;
  return TransposedMat([res,v]);
end;

#---------------------- SeriesOps -----------------------------------------
SeriesOps:=OperationsRecord("SeriesOps");

IsSeries:=s->IsRec(s) and IsBound(s.operations) and s.operations=SeriesOps;

# A d-Harish-Chandra series \CE(\BG,(\BL,\lambda)) is a record with fields:
#  .spets        \BG
#  .levi         \BL
#  .d            AsRootOfUnity(\zeta) such that L=C_G(V_\zeta)
#  .cuspidal     \lambda (index in UnipotentCharacters(.levi))
#  .principal    true iff \lambda=\Id
SeriesNC:=function(spets,levi,cuspidal,d)local s,q,e;
  if IsInt(d) and d>0 then d:=Mod1(1/d);fi;
  if not IsSpets(spets) then spets:=Spets(spets);fi;
  s:=rec(spets:=spets,levi:=levi,cuspidal:=cuspidal,d:=d,operations:=SeriesOps);
  s.principal:=UnipotentCharacters(s.levi).a[s.cuspidal]=0
               and UnipotentCharacters(s.levi).A[s.cuspidal]=0;
  return s;
end;

SeriesOps.String:=s->Format(s,rec());

SeriesOps.Print:=function(s)Print(String(s));end;

SeriesOps.Display:=function(s,opt)
  opt.screenColumns:=SizeScreen()[1];
  Print(Format(s,opt));end;

SeriesOps.Format:=function(s,opt)local m,f,uw,e,t,rowlab,res,n,cname;
  cname:=CharNames(UnipotentCharacters(s.levi),opt)[s.cuspidal];
  if IsBound(opt.TeX) then n:="\\lambda";else n:="c";fi;
  e:=String(E(Denominator(s.d))^Numerator(s.d));
  if s.spets=s.levi then 
    return SPrint(e,"-cuspidal ",cname," of ",ReflectionName(s.spets,opt));
  else
  res:=SPrint(e,"-series ");
  if IsBound(opt.TeX) then PrintToString(res,"$");fi;
  PrintToString(res,"R^{",ReflectionName(s.spets,opt),"}_",
   "{",ReflectionName(s.levi,opt),"}(");
  if IsBound(opt.TeX) then PrintToString(res,"\\lambda=");fi;
  PrintToString(res,cname,")");
  fi;
  if IsBound(s.cyclic) and s.cyclic and (s.e>3 or not IsBound(s.Hecke)) then 
    if IsBound(opt.TeX) then 
      PrintToString(res,"\\quad W_G(L,\\lambda)=Z_{",s.e,"}");
    else PrintToString(res," W_G(L,",n,")=Z",s.e);fi;
  elif IsBound(s.WGL)then 
    if IsBound(opt.TeX) then PrintToString(res,"\\quad ");fi;
    if IsBound(s.Hecke) then PrintToString(res," H_G(L,",n,")=",s.Hecke);
    elif IsBound(s.WGL.reflections) then
      PrintToString(res," W_G(L,",n,")=",ReflectionName(s.WGL,opt));
    else
      PrintToString(res," |W_G(L,",n,")|=",Size(s.WGL));
    fi;
  elif IsBound(s.WGLdims) then 
    if IsBound(opt.TeX) then PrintToString(res,"\\quad ");fi;
      PrintToString(res," |W_G(L,",n,")|=",Sum(s.WGLdims,x->x^2));
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
  if IsBound(s.translation) then
    e:=SPrint("WGLname(mod ",s.translation,")");
    f(Concatenation("\\hbox{",e,"}"),e,CharNames(s.WGL,opt));
  else f("\\hbox{WGLname}","WGLname",CharNames(s.WGL,opt));
  fi;
  f("\\varepsilon","eps",s.eps);
  if s.cyclic then
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

SeriesOps.projector:=function(s)local q;
  if IsBound(s.projector) then return s.projector;fi;
  s.ad:=Number(ReflectionDegrees(s.levi),x->x[1]=1 and
          x[2]=E(Denominator(s.d))^Numerator(s.d));
  q:=Minimum(List(ReflectionEigenvalues(s.levi),x->Number(x,y->y=s.d)));
  if q<>s.ad then Error("bad start"); fi;
  q:=PositionsProperty(ReflectionEigenvalues(s.levi),x->Number(x,y->y=s.d)=q);
# found an element w\in Levi such that V_\zeta=\zeta-eigenspace of w
# .projector= w-equivariant projector on V_\zeta
  s.projector:=EigenspaceProjector(s.spets,
      Representative(ConjugacyClasses(s.levi)[q[1]]),s.d);
  return s.projector;
end;
  
# Find W_G(s.levi,s.cuspidal)
#   as a relgroup, contains parentMap (refs->elts of W)
#   and reflists (gens->gens of parab of W)
SeriesOps.RelativeGroup:=function(s)local id,W,L,N,ud,eig,c,V,refs,r,m,v,
  getreflection,reflist,hplane,smalltobig,GetRelativeAction,rrefs,NF,LF,NFQ,q;
  if IsBound(s.WGL) then return s.WGL;fi;
  id:=IdentityMat(Rank(s.spets));
  W:=Group(s.spets);L:=Group(s.levi);
  if SeriesOps.projector(s)=id then #central series L=1
    s.WGL:=Copy(W); s.WGL.parentMap:=s.WGL.generators;
    s.WGL.reflists:=List(W.rootInclusion{W.generatingReflections},x->[x]);
    s.WGLdims:=List(CharTable(s.WGL).irreducibles,x->x[1]);
    s.cyclic:=IsCyclic(s.WGL);
    s.e:=Size(s.WGL);
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
  if Size(L)=1 then s.WGL:=N;
  elif Size(N)=Size(L) then 
    if s.levi<>s.spets then Error(s," N=L\n");return false;fi;
    s.WGL:=CoxeterGroup();s.WGL.reflists:=[];s.WGLdims:=[1];s.cyclic:=true;s.e:=1;
    return s.WGL;
  else s.WGL:=N/L; 
  fi;
  V:=NullspaceMat(SeriesOps.projector(s)-id);  # The E(d)-eigenspace
  # Ker(r-1)\cap V
  hplane:=r->SumIntersectionMat(NullspaceMat(MatXPerm(W,Reflection(W,r))-id),
     V)[2];
  m:=Concatenation(V,NullspaceMat(SeriesOps.projector(s)));
  GetRelativeAction:=function(x) # restriction of matrix x to V
    x:=m*x*m^-1;return x{[1..Length(V)]}{[1..Length(V)]};end;
  # Hplane of V to global HPlane
  smalltobig:=h->List(h,x->Concatenation(x,[1..W.rank-Length(V)]*0))*m;
  # for reflections of W with same Hplane in V return rec
  # .hom element of W above reflection of WGL
  # .refs generators of parabolic of W fixing Hplane
  # .WH  cyclic subgroup of WGL fixing Hplane
  # .eigenvalue of refl. of WGL
  # .root       of refl. of WGL
  # .coroot     of refl. of WGL
  getreflection:=function(rr)local rH,H,r,res,n;
    rH:=hplane(rr[1]); # Print("rH=",rH,"\n"); # Hyperplane of W_G(L)
    if Length(rH)=0 then # cyclic WGL
      ConjugacyClasses(L); # hack to fix a (bug in GAP?)
      H:=W;
    else H:=ReflectionSubgroup(W,Concatenation(W.rootInclusion{rr},
                              L.rootInclusion{L.generatingReflections}));
    fi;
    res:=rec(refs:=H.rootInclusion{H.generatingReflections});
    H:=Intersection(H,N);
    if Size(L)<>1 then H:=H/L;fi;
    if Size(H)=1 then Error("H=L");fi;
    if not IsCyclic(H) then Error("H not cyclic");fi;
    res.hom:=First(Elements(H),y->Order(H,y)=Size(H));
    if Size(L)=1 then r:=res.hom; else r:=res.hom.element.representative;fi;
    r:=GetRelativeAction(MatXPerm(W,r));
    Inherit(res,AsReflection(r));
    n:=AsRootOfUnity(res.eigenvalue);n:=1/Numerator(n) mod Denominator(n);
    r:=r^n; # distinguished reflection
    res.hom:=res.hom^n;
    Inherit(res,AsReflection(r)); res.WH:=H; return res;
  end;
  v:=function(h)if Length(h)=0 then return VectorSpace([[0]],Cyclotomics);
    else return VectorSpace(h,Cyclotomics);fi;end;
  rrefs:=CollectBy(refs,x->v(hplane(x)));
  rrefs:=Filtered(rrefs,x->not hplane(x[1])=V);
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
      s.cyclic:=IsCyclic(s.WGL);
      s.e:=Size(s.WGL);
      return s.WGL;
    fi;
  od;
  Error("b");
end;

# Degree in q of the parameters (normalized so the smallest is 0)
SeriesOps.mC:=function(s)local D0,aA,uc,pG,pL,xiL,xiG,e,cn,lpi;
  e:=HyperplaneOrbits(s.WGL);
  if Length(e)<>1 then return;
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

# Degree(R_s.levi^s.spets(s.cuspidal))
#  .degree       The degree of R_\BL^\BG(\lambda)=|G/L|_{q'} deg\lambda
SeriesOps.Degree:=function(s)local q;
  if IsBound(s.degree) then return s.degree;fi;
  q:=X(Cyclotomics);
  s.degree:=GenericSign(s.spets)^-1*GenericOrder(s.spets,q)/
           (GenericSign(s.levi)^-1* GenericOrder(s.levi,q));
  s.degree.valuation:=0;
  s.degree:=CycPol(s.degree)*CycPolUnipotentDegrees(s.levi)[s.cuspidal];
  return s.degree;
end;

# candidates for s.charNumbers
# list of rec(charNumbers,sch,dims,eps,span)
SeriesOps.candfromdeg:=function(s)local ad,cand,Ed,q,f,v,c,ud;
  ad:=Number(RelativeDegrees(s.levi,s.d),x->x<>1);
  ud:=CycPolUnipotentDegrees(s.spets);
  cand:=Filtered([1..Length(ud)],i->ad=Valuation(ud[i],s.d));
# now use that   S_\phi(q)=\eps_\phi Deg(RLG(\lambda))/Deg(\gamma_\phi)
  cand:=List(cand,c->rec(charNumbers:=c,sch:=Degree(s)/ud[c]));
  cand:=Filtered(cand,c->ForAll(c.sch.vcyc,x->x[2]>0)); # positive
  Ed:=E(Denominator(s.d))^Numerator(s.d);
  q:=Mvp("q");
  ad:=CycPol(q-Ed)^ad;
  f:=Degree(s)/ad;
  if ForAny(f.vcyc,x->x[2]<0) then Error(s," cuspidal is not\n"); fi;
  v:=Value(f,Ed);
  for c in cand do 
    c.span:=Degree(c.sch)-c.sch.valuation;
# now use that 
    f:=Sum(s.WGLdims,x->x^2)*Value(ud[c.charNumbers]/ad,Ed)/v;
    c.dims:=AbsInt(f);c.eps:=SignInt(f);
  od;
  cand:=Filtered(cand,c->c.dims in s.WGLdims);
  Print("after deg cands=",Length(cand),"(",Length(s.WGLdims),")\n");
  return cand;
end;

# takes a d-series s with s.spets split; fills in s.charNumbers, s.eps, s.dims
SeriesOps.CharNumbers:=function(s)
  local ad,q,ud,v,eig,f,c,t,cand,i,check,Ed,uc,g;
  if IsBound(s.charNumbers) then return s.charNumbers;fi;
  if not IsBound(s.WGL) and RelativeGroup(s)=false then return false;fi;
  s.RLG:=LusztigInduction(s.spets,UnipotentCharacter(s.levi,s.cuspidal));
  if s.RLG=false and s.levi.phi=() then 
   s.RLG:=HarishChandraInduction(s.spets,UnipotentCharacter(s.levi,s.cuspidal));
  fi;
  uc:=UnipotentCharacters(s.spets);
  if s.RLG=false then ChevieErr(s,":RLG failed\n");
  elif Degree(s)<>CycPol(Degree(s.RLG))then 
    ChevieErr(s,":Deg RLG<>Sum(eps[i]*ud[i])\n");
  else
    s.charNumbers:=Filtered([1..Length(s.RLG.v)],i->s.RLG.v[i]<>0);
    s.eps:=List(s.RLG.v{s.charNumbers},SignInt);
    s.dims:=List(s.RLG.v{s.charNumbers},AbsInt);
    s.span:=Degree(Degree(s))-Valuation(Degree(s))+uc.a-uc.A;
    return s.charNumbers;
  fi;
  # here RLG failed; we try to carry on
  cand:=SeriesOps.candfromdeg(s);
  f:=field->List(cand,x->x.(field));
  check:=function()local n;
    for n in ["charNumbers","eps","dims","span"] do s.(n):=f(n);od;
    return s.charNumbers;
  end;
  if Length(cand)=Length(s.WGLdims) then return check();fi;
  if false then
    eig:=Eigenvalues(UnipotentCharacters(s.levi),[s.cuspidal])[1];
    eig:=eig*List([1..Denominator(s.d)^2],i->E(Denominator(s.d)^2)^i);
    cand:=Filtered(cand,c->Eigenvalues(UnipotentCharacters(s.spets),
        [c.charNumbers])[1] in eig);
  # Print("after eig cands=",List(cand,c->c.charNumbers),"\n");
    if Length(cand)<Length(s.WGLdims) then 
      ChevieErr(s,":not enough left with predicted eigenvalues in ",
          FormatGAP(List(eig,AsRootOfUnity)),"\n");
      return false;
    fi;
  fi;
  SortBy(cand,x->x.dims);Sort(s.WGLdims);
  ud:=Zip(CycPolUnipotentDegrees(s.spets){f("charNumbers")},f("eps"),f("dims"),
     function(x,y,z)return x*y*z;end);
  t:=Maximum(List(ud,Degree));
  c:=p->Concatenation(Coefficients(Value(p,Mvp("q")),"q"),[1..t-Degree(p)]*0);
  v:=SubsetsSum(c(Degree(s)),List(ud,c),s.WGLdims,f("dims"));
  InfoChevie("# ",Length(v)," found\n");
  if Length(v)>10000 then
    InfoChevie("# ",Length(v)," combinations sum to dimRLG\n");
  elif Length(v)=0 then
    ChevieErr(s," no combination sums to dimRLG\n");return false;
  fi;
  
  if Length(v)=1 then cand:=cand{v[1]};return check();fi;

  InfoChevie("# after span ",Length(v)," combinations\n");
  InfoChevie("# Warning: using Mackey with tori for ",s,"\n");
  i:=FusionConjugacyClasses(s.levi,s.spets);
  c:=Zip(CharTable(s.spets).centralizers{i},CharTable(s.levi).centralizers,
    function(a,b)return a/b;end);
  t:=Zip(TransposedMat(DeligneLusztigCharacterTable(s.levi))[s.cuspidal],c,
    function(a,b)return a*b;end);
  t:=List([1..NrConjugacyClasses(s.spets)],
          k->Sum(t{Filtered([1..Length(i)],j->i[j]=k)}));
  c:=TransposedMat(DeligneLusztigCharacterTable(s.spets));
  g:=function(field,indices)return List(cand{indices},x->x.(field));end;
  v:=Filtered(v,a->Zip(g("dims",a),g("eps",a),
         function(x,y)return x*y;end)*c{g("charNumbers",a)}=t);

  if Length(v)=1 then cand:=cand{v[1]};return check();fi;

  if s.cyclic then
    InfoChevie("# Warning: using span for",s,"\n");
    s.charNumbers:=f("charNumbers"); s.dims:=f("dims"); SeriesOps.mC(s);
    if Length(v)>1 then
      v:=Filtered(v,a->g("span",a)=
        List(a,i->Sum(Difference(a,[i]),j->AbsInt(s.mC[i]-s.mC[j]))));
      if Length(v)>10000 then
        InfoChevie("# ",Length(v)," combinations have right span\n");
      fi;
    elif Length(v)=0 then
      ChevieErr(s," no combination has right span\n");return false;
    fi;
    Unbind(s.charNumbers);Unbind(s.dims);Unbind(s.mC);
  fi;

  if Length(v)=1 then cand:=cand{v[1]};return check();fi;
  if Length(v)>1 then ChevieErr(s," ",
     Join(List(FactorsSet(List(v,x->g("charNumbers",x))),x->FormatGAP(x)),"x"),
       " chars candidates\n");
  elif Length(v)=0 then ChevieErr(s," no candidates left\n");
  fi;
  return false;
end;

COMPACTCOHOMOLOGY:=true;
SeriesOps.fill:=function(s)local uc,Schur,r,p,ratio,LFrob,predeigen,map,
   series,noncus,unique,i,m,quality,rr,a,param,j,o,u,nid,FractionToRoot;
  FractionToRoot:=x->E(Denominator(x))^Numerator(x);
  if not s.cyclic then Error("fill assumes WGL cyclic\n");fi;
  if not IsBound(s.charNumbers) and CharNumbers(s)=false then return false;fi;
  if s.e=1 then s.Hecke:=Hecke(s.WGL);return s.Hecke;fi;
  if not IsBound(s.mC) then SeriesOps.mC(s);fi;
  uc:=UnipotentCharacters(s.spets);
  Schur:=CycPolUnipotentDegrees(s.spets){s.charNumbers};
  Schur:=List([1..s.e],x->Degree(s)/Schur[x]*s.eps[x]);
  s.eigen:=Eigenvalues(uc,s.charNumbers);
  LFrob:=AsRootOfUnity(Eigenvalues(UnipotentCharacters(s.levi))[s.cuspidal]);
  m:=ReflectionDegrees(Group(s.spets));
  s.delta:=Lcm(List(Filtered(ReflectionDegrees(s.spets),x->x[1]<>1),
        x->Denominator(AsRootOfUnity(x[2]))));
  rr:=function(j,i)return (i-1)/s.e-s.mC[j]*s.d;end;
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

CHEVIE.relativeSeries:=true;
SeriesOps.Hecke:=function(s)local H;
  if IsBound(s.Hecke) then return s.Hecke;fi;
  if CharNumbers(s)=false then return false;fi;
  if s.cyclic then SeriesOps.fill(s);
  elif CHEVIE.relativeSeries then 
    InfoChevie("# Relative: ",String(s),"\n");
    s.relativeSeries:=SeriesOps.RelativeSeries(s);
  fi;
  if IsBound(s.Hecke) then return s.Hecke;else return false;fi;
end;

# get Hecke of series s. If scalar-twisted untwist first
getHecke:=function(s)local t,scal,g,l,s1,p,e,c;
  t:=ReflectionType(s.spets);
  if not IsBound(s.charNumbers) then CharNumbers(s);fi;
 if Length(t)=1 and IsBound(t[1].scalar) and not ForAll(t[1].scalar,x->x=1) then
    t:=t[1];
    scal:=t.scalar;t:=ShallowCopy(t);
    InfoChevie("# removing scal=",scal,"\n");
    t.scalar:=List(t.scalar,x->1);
    if t.orbit[1].series="B" and t.orbit[1].rank=2 and t.twist=(1,2) then
      t.orbit[1].cartanType:=E(8)-E(8)^3;
    fi;
    g:=ReflectionGroup(t);
    if g=false then return false;fi;
    p:=EltWord(Group(g),GetWord(Group(s.spets),s.levi.phi/s.spets.phi));
#   p:=OnTuples(Group(s.spets).rootInclusion,s.levi.phi/s.spets.phi);
#   p:=PermListList(Group(s.spets).rootInclusion,p);
#   if p<>() then Error();fi;
    Reflections(Group(g));
    l:=SubSpets(g,List(Group(s.levi).generators,
      x->Position(Group(g).reflections,
        EltWord(Group(g),GetWord(Group(s.spets),x)))),p);
    if Length(scal)>1 then ChevieErr("scal=",scal," unimplemented\n");
      return false;
    fi;
    scal:=AsRootOfUnity(scal[1]);
    # one should do an Ennola of the Levi's cuspidal by the absorbed part
    # of scal by the center of the levi
    if not s.cuspidal in CuspidalUnipotentCharacters(l,Mod1(s.d-scal)) then
      InfoChevie("# using Ennola\n");
      e:=Ennola(Group(l));c:=AbsInt(s.cuspidal^e);
      if not c in CuspidalUnipotentCharacters(l,Mod1(s.d-scal)) then
        c:=AbsInt(s.cuspidal^(e^-1));
      fi;
    else c:=s.cuspidal;
    fi;
    s1:=SeriesNC(g,l,c,Mod1(s.d-scal));
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

# .classno  class of s.levi with \zeta-eigenspace V_\zeta
# .element  representative of classno
SeriesOps.RelativeSeries:=function(s)local res,p,ud,u1,f,o,aA,tt,eig,e,q;
  if not IsBound(s.charNumbers) then CharNumbers(s);fi;
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
  res:=List(s.WGL.reflists,function(r)local R,w,i,l;
    i:=Position(s.WGL.reflists,r);
#   Print("r:=",r,";\n");
#   Print("phi:=",s.element/s.spets.phi,";\n");
    R:=SubSpets(s.spets,r,s.element/s.spets.phi);
#   Print("R:=",FormatGAP(R),";\n");
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
    p:=SeriesNC(R,SubSpets(R,l,s.element/R.phi),s.cuspidal,s.d);
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
    p.cyclic:=ForAll(p.WGLdims,x->x=1);
    p.e:=Sum(p.WGLdims,x->x^2);
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
  u1:=List(u1,x->Degree(s)/CycPol(Mvp(x)));
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

# Series(spets,levi,cuspidal,d) or Series(spets[,d[,ad]])
Series:=function(arg)local s,d;
  if Length(arg)=4 then s:=ApplyFunc(SeriesNC,arg); Hecke(s); return s;
  elif Length(arg)=1 then d:=1;else d:=arg[2];fi;
  return List(ApplyFunc(CuspidalPairs,arg),x->Series(arg[1],x[1],x[2],d));
end;

#--------------------------------------------------------------------
ProperSeries:=function(W)local l;
  l:=Set(List(Flat(ReflectionEigenvalues(W)),Denominator));
  return Concatenation(List(l,d->Concatenation(List(
      [1..Length(RelativeDegrees(W,d))],i->Series(W,d,i)))));
end;
ProperSeriesNC:=function(W)local l;
  l:=Set(List(Flat(ReflectionEigenvalues(W)),Denominator));
  return Concatenation(List(l,d->Concatenation(List(
      [1..Length(RelativeDegrees(W,d))],i->
   List(CuspidalPairs(W,d,i),p->SeriesNC(W,p[1],p[2],i))))));
end;
