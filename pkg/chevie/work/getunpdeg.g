# Programs to get the list of ud (Unipotent Degrees) + Frobenius eigenvalues
# as outlined in ``Split Spetses...''

# returns string description as list of ranges of sorted list of integers
FormatAsRanges:=function(l)local res,er,br,i;
  if Length(l)=0 then return "[]";fi;
  res:=[]; br:=l[1];er:=br;
  for i in [2..Length(l)] do
    if l[i]<>er+1 then Add(res,[br,er]);br:=l[i];er:=br;
    else er:=er+1;
    fi;
  od;
  Add(res,[br,er]);
  return Join(List(res,function(x)
    if x[1]+1<x[2] then return Join(x,"..");
    elif x[1]+1=x[2] then return Join(x);
    else return String(x[1]);
    fi;end));
end;

FractionToRoot:=function(d) # computes exp(2i\pi d) for Rational d
  if IsInt(d) and d<>0 then d:=Mod1(1/d);fi;
  return E(Denominator(d))^Numerator(d);
end;

# i^n where n fraction i AsRootOfUnity computed in 'right' way
# so it agrees with GetRoot
PowerRoot:=function(i,n)local d,j,n1,k;
  d:=Denominator(i);j:=1;
  n1:=Denominator(n);
  if n1=1 then return Mod1(n*i);fi;
  repeat k:=Gcd(n1,d);n1:=n1/k;j:=j*k;until k=1;
  return Mod1((Numerator(i)*Numerator(n)*GcdRepresentation(n1,d)[1])/j/d);
end;

# returns the list of [a,A] not fully accounted by the unipotent
# characters of indices l
# Test for completion is \sum_l|ud|^2=\sum_l|fd|^2
# The fds are read from harishChandra[1] but could be from the knowledge
# of Rouquier families
UnaccountedAa:=function(W,l)local ud,uw,aA,h,q,fd,remain;
  uw:=UnipotentCharacters(W);ud:=CycPolUnipotentDegrees(W);
  aA:=TransposedMat([uw.a,uw.A]);
  h:=uw.harishChandra[1].charNumbers;
  aA:=List(Set(aA{h}),x->rec(a:=x[1],A:=x[2],
     fd:=Filtered([1..Length(h)],i->aA[h[i]]=x),ud:=Filtered(l,i->aA[i]=x)));
  q:=X(Cyclotomics);q.name:="q";fd:=FakeDegrees(W,q);
  remain:=Filtered(aA,function(r)local p;
    p:=Sum(fd{r.fd},x->x^2)-Sum(ud{r.ud},x->Value(x*ComplexConjugate(x),q));
    return p<>0*p;end);
  return List(remain,x->[x.a,x.A]);
end;

eW:=function(W) 
  if IsSpets(W) then return eW(Group(W));fi;
  return Sum(ReflectionDegrees(W)+ReflectionCoDegrees(W));
end;

# A *parameter list* is a list of pairs [f,d] of rationals representing
# a list of parameters FractionToRoot(f)*q^d (for an Hecke algebra)

# tests if a list of parameters can give for prod(x-param) a polynomial in x,q 
# That is for each non integral d we must have corresponding Sum(E(f))=0.
RationalMinPolParam:=p->ForAll(Filtered(Set(List(p,x->x[2])),x->not IsInt(x)),
  v->Sum(Filtered(p,x->x[2]=v),x->FractionToRoot(x[1]))=0);

# returns CycPol(i-th schur elem) of Hecke(Z/Length(p),p) for parameter
# list p
CycPolCyclicSchur:=function(p,i)
  p:=CollectBy(List(Drop(p,i),x->p[i]-x),x->x[2]);
  p:=List(p,x->[x[1][2],List(x,y->y[1])]);
  return Product(p,function(l)local v;
    if IsInt(l[1]) then return 
      Product(l[2],x->CycPol(1-FractionToRoot(x)*Mvp("x")^l[1]));
    fi;
    v:=Product(l[2],x->1-FractionToRoot(x)*Mvp("x")^l[1]);
    if ForAny(v.elm,y->not ForAll(y.coeff,IsInt))
    then return CycPol(0);
    else return CycPol(v);fi;end);
end;

# arrangements of val completed by an arrangement of other with in the
# resulting arrangement val[i] positioned at an element of pos[i].
# this is recursive and positions which serves just for recursion is 
# initialized to [0..Length(val)+Length(other)]
solveit:=function(positions,val,pos,other)local p,pp,v,v1,res,a,head,tail;
  if Length(val)=0 then return List(Arrangements(other,Length(other)),
      a->TransposedMat([a,positions]));fi;
  p:=pos[1];
  pp:=Filtered([1..Length(val)],i->pos[i]=p);
  v:=val{pp};
  res:=[];
  for a in Arrangements(other,Length(p)-Length(pp)) do
    v1:=Concatenation(v,a);
    head:=List(Arrangements(v1,Length(v1)),x->TransposedMat([x,p]));
    tail:=solveit(Difference(positions,p),
         Drop(val,pp),Drop(pos,pp),DifferenceMultiSet(other,a));
    Append(res,List(Cartesian(head,tail),Concatenation));
  od;
  return res;
end;

# Returns [[partition of [1..|Irr(W)|] in Rouquier families],
# [list of corresponding a+A]]
# (this could be obtained in theory from Maria's program with no
# knowledge of UnipotentCharacters)
RouquierFamilies:=function(W)local uc,ff,fd;
  uc:=UnipotentCharacters(W);
  ff:=List(uc.families,x->x.charNumbers);
  ff:=List(ff,x->Filtered(x,y->y in uc.harishChandra[1].charNumbers));
  ff:=List(ff,x->List(x,y->Position(uc.harishChandra[1].charNumbers,y)));
  fd:=FakeDegrees(W,X(Cyclotomics));
  return [ff,List(ff,x->Minimum(List(fd{x},y->y.valuation)))+
       List(ff,x->Maximum(List(fd{x},Degree)))];
end;

# determine d-principal series of W: d should be a regular number such that 
# W_d has <=2 hyperplane orbits 
# known: list of indices of already known (deg,eig) of unipotent characters
# in case W_d has 2 hplane orbits we need known to include indices of all 
# unipotent chars corresponding to linear chars of W_d
FindRelativeHecke:=function(W,d,known)local ud,uw,s,ser,par,tbl,i,j,poss,a,p,H,
  res,sch,cnt,c,bad,ho,ss,tocyc,pp,pars,omegachi,Ed,degparam;
  s:=PrincipalSeries(W,d);
  if Size(s.levi)>1 then Error(d," is not regular\n");fi;
  Print("For ",d,"-series W_G(L)=",Format(RelativeGroup(s)),
    " find Hecke parameters\n");
  ho:=HyperplaneOrbits(s.WGL);
  ud:=CycPolUnipotentDegrees(W);uw:=UnipotentCharacters(W);
  Ed:=FractionToRoot(d);
  sch:=List(ud,x->SignInt(Value(x,Ed))*s.degree/x);
  # known parameters <-> linear chars of WGL corresp to known => ser.ind
  ser:=rec(ind:=Filtered(known,i->AbsInt(Value(ud[i],Ed))=1));
  ser.aA:=uw.a{ser.ind}+uw.A{ser.ind};
  ser.degparam:=(eW(W)-ser.aA)/(ho[1].e_s*ho[1].N_s);
  ser.eigen:=List(Eigenvalues(uw,ser.ind),AsRootOfUnity);
# degparam(W,d[,l]) 
#  d in Q/Z describes \zeta regular eigenvalue
#  l: known unipotent degrees given as CycPol list (assumed empty if not given)
# returns degrees of parameters for H_W(\zeta) for each Hplane orbit in
# W(\zeta), each list encoded as [[d,n],..]
#   where n is a majoration of the multiplicity of degree d
#
# A better majoration of the multiplicity is obtained if l is nonempty.
# We use formula 5.22 in "Split Spetses":
# \sum_{\chi\in\CF}|\Feg(\chi)(\zeta)|^2=
# \sum_{\psi\in W_d|\rho_\psi\in\CF}\psi(1)^2
# and that \rho(\zeta)=0 unless \rho=\rho_\psi and \psi(1) else.
#
# We need that W(\zeta) has <=2 orbits. An exact list is returned if
# 2 orbits Hplanes (then assume all \rho_\psi for \psi linear known): 
# or if W(\zeta) is cyclic.
degparam:=function(arg)
  local W,d,ff,fd,faA,n,s,aA,ho,v,max,Ed,decomposeassum,lin;
decomposeassum:=function(l,card,max)local t;Sort(l);
# Let  l  be  a  multiset  (list  with  repetitions) of length
# Product(card) of rationals in [0,Sum(max)].
# If  it is  possible to  write l=List(Cartesian(S1,S2),Sum)  where Si is a
# list  of length card[i] of rationals  in [0,max[i]] then returns S_1, S_2
# else gives Error.
  t:=Cartesian(List([1..2],i->Filtered(Combinations(l,card[i]),
    x->Minimum(x)=0 and Maximum(x)=max[i])));
  t:=Filtered(t,i->Collected(List(Cartesian(i),Sum))=Collected(l));
  if Length(t)<>1 then Error(Length(t)," decompositions as sum");fi;
  return t[1];
end;
  W:=arg[1];d:=arg[2];Ed:=FractionToRoot(d);
  if Length(arg)=3 then aA:=TransposedMat(
    [eW(W)-List(arg[3],x->Degree(x)+Valuation(x)),
           List(arg[3],x->Value(x,Ed)^2)]);
    lin:=Filtered(aA,x->x[2]=1);
  else aA:=[];lin:=[];fi;
  s:=PrincipalSeries(W,d);
  ho:=HyperplaneOrbits(RelativeGroup(s));
  if Length(ho)>1 then
    if Length(lin)<>Product(ho,x->x.e_s) then 
      Error("unimplemented: only ",Length(lin)," out of ",
      Product(ho,x->x.e_s)," linear characters known and >=2 Hplane orbits\n");
      return false;
    fi; # else all linear chars known
    max:=List(ho,o->eW(ReflectionSubgroup(W,s.WGL.reflists[o.s]))*o.N_s);
    max:=decomposeassum(List(lin,x->x[1]),List(ho,x->x.e_s),max);
    return List([1..Length(ho)],i->max[i]/ho[i].e_s/ho[i].N_s);
  fi;
# Print("maxdeg=",eW(s.WGL)/ho[1].e_s,"\n");
  ff:=RouquierFamilies(W);faA:=eW(W)-ff[2];ff:=ff[1];
  fd:=FakeDegrees(W,X(Cyclotomics));
  n:=List(ff,x->Sum(fd{x},y->Value(y,Ed)*Value(y,Ed^-1)));
  n:=TransposedMat([faA/eW(s.WGL),n]);
  n:=Filtered(n,x->x[2]<>0 and x[1]>=0);
  n:=List(CollectBy(n,x->x[1]),y->[y[1][1],Sum(y,x->x[2])]);
  if IsBound(s.e) then return n;fi; #cyclic case
  aA:=Filtered(aA,x->x[2]>1 and x[1]>=0);
  aA:=List(CollectBy(aA,x->x[1]),y->[y[1][1],Sum(y,x->x[2])]);
  aA:=List(aA,x->[x[1]/eW(s.WGL),x[2]]);
  for v in aA do 
   if v[2]>1 then p:=PositionProperty(n,x->x[1]=v[1]); n[p][2]:=n[p][2]-v[2];fi;
  od;
  n:=Filtered(n,x->x[2]<>0 and x[1]>=0);
  Sort(n);
  return n;
end;
  p:=degparam(W,d,ud{known}); # list of [degparam, majoration of multiplicity]
  if p=false then return false;fi;
  if IsBound(s.e) then # cyclic case
    Print( "\t",Length(ser.ind)," out of ",s.e," schur elements known\n");
  # In the cyclic case, for a parameter x=\zeta_e^i(\zeta\inv q)^m_i,
  # where AsRootOfUnity(\zeta)=d
  # the Frobenius eigenvalue Fr is \zeta^{(e_W-em_i)*d+i}, 
  # so Fr+em_i*d^2-e_W*d^2=i*d (mod 1) or if d=num/den then 
  # Fr*num'*den+(em_i-e_W)*d=i (mod den) where num'=num^(-1) (mod den)
  # ser.specialization contains the possible values of i for a param
    if d=0 then ser.specialization:=
      List(ser.degparam*s.e-eW(W),x->[0..s.e-1]+(x mod 1));
    else
    ser.specialization:=List(ser.eigen*Denominator(d)*
    QuotientMod(1,Numerator(d),Denominator(d))+ser.degparam*s.e*d
    -eW(W)*d,x->[0,Denominator(d)..s.e-Denominator(d)]+(x mod Denominator(d)));
    fi;
    poss:=Position(ser.ind,Position(ud,CycPol(1)));
    ser.specialization[poss]:=[0];
    poss:=Difference([1..Length(ser.specialization)],[poss]);
    ser.specialization{poss}:=
      List(ser.specialization{poss},x->Filtered(x,y->y<>0));
#   tbl:=TransposedMat([ser.eigen,ser.degparam,ser.specialization]);
#   SortBy(tbl,x->x[3]);
#   Print("Known part of ",d,"-series\n",FormatTable(tbl,
#      rec(rowLabels:=ser.ind,columnLabels:=["eig","degparam","specialize"])));
    p:=Concatenation(List(p,x->[1..x[2]]*0+x[1]));
    p:=DifferenceMultiSet(p,ser.degparam);
    if Length(p)=0 then # all are known
      p:=List(sch{ser.ind},i->Position(sch{ser.ind},i));
      poss:=FitParameter(sch{ser.ind},ser.degparam,d);
      poss:=Concatenation(List(poss,TransposedMat));
      SortBy(poss,x->x[2]); poss:=List(poss,y->ser.degparam[y[1]]);
      p:=Position(poss,Maximum(poss));
      if p<>1 then SortParallel(List(E(s.e)^(1-p)*List([0..s.e-1],x->E(s.e)^x),
        AsRootOfUnity),poss);
      fi;
      poss:=[poss];
    else poss:=solveit([0..s.e-1],ser.degparam,ser.specialization,p);
      poss:=List(poss,function(v)local res,x;res:=[];
	for x in v do res[x[2]+1]:=x[1];od;return res;end);
    fi;
   # get pairs [root of unity part of ith param, m_i]
    poss:=Filtered(poss,x->Maximum(x)=x[1]);
    if Length(poss)>1 then 
      Print( "\t",Length(poss)," parameter lists to test\n");
    fi;
    cnt:=0;
    poss:=List(poss,function(l)
      cnt:=cnt+1;if cnt mod 1000=0 then Print(".\c");fi;
     return List([1..Length(l)],
      i->[Mod1((i-1)/s.e-PowerRoot(d,l[i])),l[i]]);end);
    poss:=Filtered(poss,RationalMinPolParam);
    if Length(poss)>1 then 
      Print( "\t",Length(poss)," give an algebra defined over C[q]\n");
    fi;
    cnt:=0;
    poss:=Filtered(poss,function(x)local k;
      cnt:=cnt+1;if cnt mod 100=0 then Print(".\c");fi;
      if not CycPolCyclicSchur(x,1) in sch then return false;fi;
      return true;end);
    if Length(poss)>1 then Error("more than one poss");fi;
    H:=Hecke(ComplexReflectionGroup(s.e,1,1),
      [List(poss[1],x->FractionToRoot(x[1])*Mvp("q")^x[2])]);
    Print("found \n");Cut(Format(H));
    p:=List(SchurElements(H),CycPol);
    if Length(p)=2 then p:=Reversed(p);fi;
    if ForAny([1..Length(ser.ind)],i->not sch[ser.ind[i]]in p{ser.specialization[i]+1})
    then CHEVIE.Check.EqLists(sch{ser.ind},p{List(ser.specialization,x->x[1])+1});
      Error("expected not in schur");fi;
  elif Length(ho)=1 then # one orbit of Hplanes
    p:=Concatenation(List(p,x->List([1..
      Minimum(x[2],ho[1].e_s-Length(ser.eigen))],i->x[1])));
    p:=Filtered(Combinations(p,ho[1].e_s-Length(ser.eigen)),x->
        Sum(x)+Sum(ser.degparam)=W.N/ho[1].N_s);
    p:=List(p,x->Concatenation(ser.degparam,x));
  # Print("p=",p,"\n");
  # here each line of p is possible degparam for other parameters
    res:=[];
    pp:=CharParams(s.WGL){Concatenation([PositionId(s.WGL)],ho[1].det_s)};
    for c in p do
      pars:=List(Arrangements(c,ho[1].e_s),a->
        List([0..ho[1].e_s-1],i->FractionToRoot(i/ho[1].e_s-PowerRoot(d,a[i+1]))*Mvp("q")^a[i+1]));
      Print( "\t",Length(pars)," to try:\c");
      for par in pars do
        Print(".\c");
        H:=Hecke(s.WGL,[par]);
        ss:=List(pp,x->SchurElement(H,x));
        if ForAll(ss,x->ForAll(x.elm,y->ForAll(y.coeff,IsInt))) then
          ss:=List(ss,CycPol);
          if ForAll(ss,i->i in sch) then Add(res,par);fi;
        fi;
      od;
    od;
    Print(Length(res)," selected\n");
    res:=Filtered(res,function(par) H:=Hecke(s.WGL,[par]);
      ss:=List(SchurElements(H),CycPol); 
      return ForAll(sch,i->i=0*i or i in ss);
     end);
    # Print("possible params:",res,"\n");
    # test that res[1] has highest degree param monic
    res:=Filtered(res,function(p)
      p:=Filtered(p,x->Degree(x)=Maximum(List(p,Degree)));
      if Length(p)>1 then Error("more that one param of max. degree");fi;
      return Value(p[1],["q",Ed])=1;end);
    if Length(res)>1 then Error("possibilities:",res);fi;
    H:=Hecke(s.WGL,res);
  else # for Length(ho)>1 we assume p contains correct list of degparam
    p:=List(p,Reversed);
    pars:=Cartesian(List([1..Length(ho)],i->List(Arrangements(p[i]
     {[2..ho[i].e_s]},ho[i].e_s-1),a->Concatenation([p[i][1]],a))));
          Print( "\t",Length(pars)," to try:\c");
    pars:=List(pars,x->List([1..Length(ho)],j->List([0..ho[j].e_s-1],
       i->FractionToRoot(i/ho[j].e_s-PowerRoot(d,x[j][i+1]))*Mvp("q")^x[j][i+1])));
      Print( "\t",Length(pars)," to try:\c");
    res:=[];
    pp:=[CharParams(s.WGL)[PositionId(s.WGL)]];
    for par in pars do
      Print(".\c");
      i:=[];
      for j in [1..Length(par)] do i[ho[j].s]:=par[j];od;par:=i;
      H:=Hecke(s.WGL,i);
      ss:=List(pp,x->SchurElement(H,x));
      if ForAll(ss,x->ForAll(x.elm,y->ForAll(y.coeff,IsInt))) then
	ss:=List(ss,CycPol);
	if ForAll(ss,i->i in sch) then Add(res,par);fi;
      fi;
    od;
    if Length(res)>1 then Error("possibilities:",res);fi;
    H:=Hecke(s.WGL,res[1]);
  fi;
  return H;
end;

# EigenAndDegHecke(s) 
# s is a series. Uses only s.Hecke, s.d, 
# (and s.degree, s.levi and s.cuspidal for non-principal series)
# 
# Assume H=s.Hecke is a spetsial exp(2i pi d)-cyclotomic algebra, 
# for d in Q/Z describing a central element of Group(H).
# Returns a list of records for each character chi of H with fields
#   deg:=Generic degree = S_Id/S_chi
#   eig:=rootUnity part of Eigenvalue of Frobenius
#   frac:= fractional power in q-part of Eigenvalue of Frobenius
EigenAndDegHecke:=function(s)
  local W,H,d,i,omegachi,om2,ss,frac,ct,ct1,ct2,p,n,good,zeta,d1;
  H:=s.Hecke;d:=s.d;
  if IsInt(d) and d<>0 then d:=1/d;fi;
  W:=Group(H);ct:=CharTable(H).irreducibles;
  ct1:=CharTable(W).irreducibles;
  zeta:=FractionToRoot(d);
  ct2:=ScalMvp(Value(ct*Mvp("q")^0,["q",zeta]));
  n:=[1..Length(ct2)];
  good:=Filtered(n,i->not ForAny(ct2{n}[i],IsUnknown));
  p:=PermListList(ct1{n}{good},ct2{n}{good}); #Permuted(ct,p) specializes
  if p<>() then Print("***** perm=",p,"\n"); 
    if IsCyclic(W) then ChevieErr("should not have perm");fi;fi;
  d1:=d*Denominator(d)/Gcd(Denominator(d),OrderCenter(W));
  i:=PositionRegularClass(W,d1); # this class represents F
  omegachi:=List(Permuted(ct,p),x->Mvp(x[i]/x[1]).coeff[1]);
  frac:=List(HeckeCentralMonomials(H),Degree)*d1;
  om2:=Zip(List(ct1,x->x[i]/x[1]),frac,function(o,p)return o/
      FractionToRoot(PowerRoot(d,p));end);
  if omegachi<>om2 then omegachi:=om2;
    ChevieErr(ChevieClassInfo(W).classtext[i],"^",1/d1,
      " not equal to pi(",ReflectionName(W),")\n");
  fi;
  ss:=List(Permuted(SchurElements(H),p),CycPol);
  ss:=List(ss,x->s.degree/x);
# omegachi:=omegachi*Eigenvalues(UnipotentCharacters(s.levi))[s.cuspidal];
  zeta:=FractionToRoot(PowerRoot(d1,frac[PositionId(W)]));
  if IsBound(s.delta) and  s.delta<>1 and IsCyclic(W) then
    omegachi:=List([1..s.e],
       i->FractionToRoot(s.d*s.e*s.delta*((i-1)/s.e-s.mC[i]*s.d)));
    zeta:=FractionToRoot(s.d*s.e*s.delta*s.d*s.mC[1]);
  fi;
  return Zip(ss,omegachi,frac,function(deg,eig,frac)
     return rec(deg:=deg,eig:=eig*zeta,frac:=Mod1(frac));end);
end;

# Given a exp(2i\pi d)-cyclotomic algebra H with parameters monomials in q,
# Ennola-twist it by exp(2i\pi z) getting an exp(2i\pi (d+z))-cyclotomic algebra
# here d and z are Rationals.
EnnolaTwistHecke:=function(H,z,d)local zeta;
  zeta:=FractionToRoot(Mod1(d+z));
  z:=Value(H.parameter*Mvp("q")^0,["q",Mvp("q")/FractionToRoot(z)]);
  z:=List(z,function(p)
    SortBy(p,x->AsRootOfUnity(ScalMvp(Value(x,["q",zeta]))));return p;end);
  H:=Hecke(Group(H),z);
# Print("by ",z,"-twisting:",Format(H),"\n");
  return H;
end;

# Find (unipotent degrees, eigenvalues) for W
# see "split spetses" chapter 5
# getunpdeg(W[,opt]) 
# if given opt=rec(principal:=list of indices of 1-series to use)
#                             ie indices in uc.harishChandra
getunpdeg:=function(arg)
  local W,H,chars,d,n,uc,h,z,k,l,s,process,opt,ser,sers,finddeinW,cyclic,todo,e;
  # find positions in Uch(W) of elts of list il of rec(deg:=+-xx[,eig:=xx])
  finddeinW:=function(s,il)local l,ud,res,p,m,eigs,uw,pos,r;
    res:=[];ud:=CycPolUnipotentDegrees(W);l:=ShallowCopy(il);
    uw:=UnipotentCharacters(W);eigs:=Eigenvalues(uw);
    pos:=[1..Length(l)];
    while Length(pos)>0 do
      r:=l[pos[1]];
      m:=Filtered(pos,i->l[i].deg in [r.deg,-r.deg]);
      if IsBound(r.eig) then m:=Filtered(m,i->l[i].eig=r.eig);fi;
      p:=List(PositionsSgn(ud,r.deg),AbsInt);
      if IsBound(r.eig) then p:=Filtered(p,y->eigs[y]=r.eig);fi;
      if Length(p)<>Length(m) then 
	Print("W has ",Length(p)," with deg=",r.deg," eig=",r.eig,
	      " when ",Length(m)," occur in given list\n");
	p:=List(PositionsSgn(ud,r.deg),AbsInt);
	Print("they are in positions ",p," with eig ",eigs{p},"\n");
	Print("this error is for ",Position(il,r),"-th in given list\n");
	Error();
      else
        if Length(m)>1 then 
          if not IsBound(s.ambig) then s.ambig:=[];fi;
          Add(s.ambig,[m,p]);
        fi;
	res{m}:=p;
      fi;
      pos:=Difference(pos,m);
    od;
    s.charNumbers:=res;
  end;
  ser:=function(d,H)local s,e;
    s:=Series(Spets(W),
      SubSpets(Spets(W),[],Representative(ConjugacyClasses(W)[
        PositionRegularClass(W,d)])),1,d);
    RelativeGroup(s);
    s.Hecke:=H;
    e:=EigenAndDegHecke(s);
    s.eigen:=List(e,x->x.eig);
    s.frac:=List(e,x->x.frac);
    finddeinW(s,e);
    return s;
  end;
  W:=arg[1];
  if Length(arg)=2 then opt:=arg[2];else opt:=rec(ennola:=true);fi;
  chars:=[];
  uc:=UnipotentCharacters(W);
  process:=function(s)local d,n,nbu; nbu:=Size(uc);
    n:=Difference(s.charNumbers,chars);
    d:=Difference(n,chars);
    if Length(d)<>0 then 
      Print("By ",s.d,"-series found ",FormatAsRanges(d),"\n");
    fi;
    p:=PositionProperty(sers,ss->ss.d=s.d and ss.cuspidal=s.cuspidal);
    if p=false then Add(sers,s);
    else if sers[p].Hecke.parameter<>s.Hecke.parameter then Error();fi;
    fi;
    UniteSet(chars,n);
    if Length(chars)=nbu then 
      if Length(d)<>0 then Print("****** Found all!!!\n");fi;
      return true;
    elif Length(d)>0 then Print("  --- still to find: ",
      FormatAsRanges(Difference([1..nbu],chars)),"\n");
      return false;
    fi;
  end;
  z:=OrderCenter(W); H:=Hecke(W,Mvp("q"));
  sers:=[];
  if IsBound(opt.ennola) then
    for k in [0..z-1]/z do process(ser(k,EnnolaTwistHecke(H,k,0)));od;
  else process(ser(0,H));
  fi;
  if IsBound(opt.principal) then
    for k in opt.principal do
      h:=uc.harishChandra[k];
      if h.relativeType<>[] and h.levi<>[] then
        H:=SubSpets(Spets(W),h.levi);
        s:=Series(Spets(W),H,FindCuspidalInLevi(h.cuspidalName,H),0);
        RelativeGroup(s);
        s.Hecke:=UnipotentCharactersOps.RelativeHecke(uc,k,Mvp("q"));
        e:=EigenAndDegHecke(s);
        s.eigen:=List(e,x->x.eig);
        s.frac:=List(e,x->x.frac);
        finddeinW(s,e);
        Print("Using ",s,"\n");
        process(s);
      fi;
    od;
  fi;
  cyclic:=Filtered(RegularEigenvalues(W),x->Length(RelativeDegrees(W,x))=1);
  todo:=Concatenation(List(Reversed(cyclic),d->List(PrimeResidues(d),k->k/d)));
  for d in todo do
    H:=FindRelativeHecke(W,d,chars);
    for k in [0..z-1] do 
      process(ser(Mod1(d+k/z),EnnolaTwistHecke(H,k/z,d)));
      todo:=Difference(todo,[Mod1(d+k/z)]);
    od;
  od;
  todo:=Reversed(Difference(RegularEigenvalues(W),
    Union(cyclic,DivisorsInt(z))));
  todo:=Union(List(todo,d->List(PrimeResidues(d),k->k/d)));
  for d in todo do
    H:=FindRelativeHecke(W,d,chars); 
    if H=false then Print("**** Could not compute ",d,"-series!!!!\n");
    else process(ser(d,H));todo:=Difference(todo,[d]);fi;
  od;
  return [sers,chars];
end;

# find pairs of families with same a and A
dbfams:=function(W)local uc,aA,db,l,f;
  uc:=UnipotentCharacters(W);
  aA:=List(uc.families,
    f->[uc.a[f.charNumbers[1]],uc.a[f.charNumbers[1]]]);
  db:=List(Filtered(Collected(aA),x->x[2]>1),x->
    Filtered([1..Length(aA)],i->aA[i]=x[1]));
  for l in db do
    Print("families ",l," [a,A,#]:");
    for f in l do
      f:=uc.families[f].charNumbers;Print([uc.a[f[1]],uc.A[f[1]],Length(f)]);
    od;
    Print("\n");
  od;
end;

# given a Hecke algebra [where parameters are monomials in q], 
# find the minimal field containing coeffs of prod(x-p_i) where p_i 
# are the parameters normalized so that highest degree has coeff 1.
fh:=function(H)local l;
  l:=List(H.parameter,x->x/x[1].coeff[1]);
  l:=Set(List(l,x->Product(Mvp("x")-x)));
  l:=List(l,x->NF(Set(x.coeff)));
  l:=List(l,String);
  if Length(Set(l))=1 then return l[1];else return Join(l);fi;
end;

testg:=function(arg)local l;
  l:=getunpdeg(ApplyFunc(ComplexReflectionGroup,arg))[1];
  l:=Filtered(l,x->x<>false);
  Print("-----------------------------------------\n");
  Print(Format(List(l,x->[x,fh(Hecke(x))])),"\n");
end;

# s is a series returned by getunpdeg; show qEigen
CheckFrac:=function(W)local ct,r,s,uc,i,v;
  r:=getunpdeg(W)[1];uc:=UnipotentCharacters(W);
  for s in r do
    ct:=CharTable(s.Hecke).irreducibles*Mvp("q")^0;
    ct:=List(ct,x->Set(Flat(List(x,y->List(y.elm,z->z.coeff)))));
    ct:=List(ct,x->Union(Set(List(x,Mod1)),[0]));
    ct:=List(ct,x->Mod1(1/Lcm(List(x,Denominator))));
    s.irch:=ct;
  od;
  s:=[];
  for i in [1..Size(uc)] do
    v:=Filtered(r,s->i in s.charNumbers);
    v:=List(v,function(s)local p;
      p:=Position(s.charNumbers,i);
      return [s.d,s.irch[p],s.frac[p]];end);
    if Length(Set(List(v,x->x[2])))>1 or Length(Set(List(v,x->x[3])))>1 then
      Error();fi;
    if Length(v)>0 then
    Add(s,rec(no:=i,sers:=Set(List(v,x->x[1])),irch:=v[1][2],frac:=v[1][3],
      fam:=PositionProperty(uc.families,x->i in x.charNumbers)));
    fi;
  od;
  if ForAny(s,function(x)local h;
    h:=First(uc.harishChandra,h->x.no in h.charNumbers);
    if x.frac<>h.qEigen then Error("for ",x.no," computed qEig=",
      x.frac," while stored qEig=",h.qEigen,"\n");fi;
    return x.frac<>h.qEigen;end) then Error("oopos");fi;
  s:=Filtered(s,x->x.irch<>0 or x.frac<>0);
  Print(FormatTable(List(s,x->[x.irch,x.frac,x.fam,x.sers]),
     rec(rowLabels:=List(s,x->x.no), rowsLabel:="no",
    columnLabels:=["irr.ch.","frac.eig.","family","series"])),"\n");
end;

findzeta:=function(W)local ct,r,s,uc,i,v;
  r:=getunpdeg(W)[1];uc:=UnipotentCharacters(W);
  s:=[];
  for i in [1..Size(uc)] do
    v:=Filtered(r,s->i in s.charNumbers);
    s[i]:=List(v,s->FractionToRoot(s.d)^(uc.a[i]+uc.A[i]
         -Sum(ReflectionDegrees(W)+ReflectionCoDegrees(W))));
  od;
  return s;
end;

checkseries:=function(s)local e,n,W;
  W:=s.spets;
  if IsBound(s.eigen) then
    e:=Eigenvalues(UnipotentCharacters(W)){s.charNumbers};
    if s.eigen<>e then
      n:=CharNames(UnipotentCharacters(W)){s.charNumbers};
      ChevieErr(s," actual eigen differ from predicted eigen");
      CHEVIE.Check.EqTables(
                 rec(rowLabels:=["actual"],columnLabels:=n,scalar:=[e]),
                 rec(rowLabels:=["pres"],columnLabels:=n,scalar:=[s.eigen]));
    fi;
  else Print("no eigen");
  fi;
end;
