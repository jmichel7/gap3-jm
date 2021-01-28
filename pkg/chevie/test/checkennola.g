
# According to Michel (principal series) and Gunter (other cases)
# the eigenvalue of Ennola_E(z)^k(rho_chi) is  
#  EigenEnnola(W,i,k)/Ennolascalar(rho,E(z)^k)
# where z=|ZW| and rho_chi is i-th unip. char.
# check that Ennola is multiply by Ennola(special) in fusion algebra
CHEVIE.AddTest("Ennola",
function(W)local uc,en,f,i,B,r,es,stored,implemented,rr,t;
  if not IsSpets(W) then W:=Spets(W);fi;
  t:=W.type[1];
  if t.orbit[1].series in ["G","I","F"] and t.twist<>() then
    Print("not applicable\n");return;
  fi;
  implemented:=not IsBound(t.orbit[1].p) and not t.orbit[1].series in ["B","D"];
  en:=SpetsEnnola(t,false);
  if ForAny(en,x->Length(x)>1) then ChevieErr("Ambiguity in Ennola:",en,"\n");fi;
  en:=Cartesian(en);
  if en=[] then ChevieErr("Ennola failed for ",W,"\n");return;fi;
  en:=en[1];
  if not implemented then Print("Ennola=",en,"\n");return; fi;
  stored:=not t.orbit[1].series in ["A"];
  if not stored then
    if SpetsEnnola(t)<>Ennola(W) then Error("value of Ennola");fi;
    return;
  fi;
  uc:=UnipotentCharacters(W);
  for i in [1..Length(uc.families)] do
    f:=uc.families[i];
    if IsBound(f.ennola) then stored:=f.ennola;else stored:=f.special;fi;
    if implemented and not stored=en[i] then 
      Error("for family ",i,"=",TeXStrip(f.name)," ennola stored=",stored,
            " ennola=",en[i],"\n");
    fi;
  od;
end,
W->IsSpetsial(W) and Length(W.type)=1);

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
  local W,H,d,i,omegachi,om2,ss,frac,ct,ct1,ct2,p,n,good,zeta,d1,FractionToRoot;
  FractionToRoot:=x->E(Denominator(x))^Numerator(x);
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
  d1:=d*Denominator(d)/Gcd(Denominator(d),OrderCenter(W.type[1]));
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

CheckSerie:=function(s)local W,l,s,c,e,Ed,pred,n;
  W:=s.spets;
  InfoChevie("\n   # ",String(s));
  RelativeGroup(s);
  if s.principal then
    Ed:=E(Denominator(s.d))^Numerator(s.d);
    c:=Value(GenericOrder(s.spets,X(Cyclotomics))/Sum(s.WGLdims,x->x^2)/
    GenericOrder(s.levi,X(Cyclotomics)),Ed);
    if c<>1 then  ChevieErr(s," => |G|/(|L||WGL|) mod Phi=",c,"\n");fi;
    e:=GenericSign(s.spets)*Ed^(Sum(ReflectionCoDegrees(Group(s.spets)))-
       Sum(ReflectionCoDegrees(Group(s.levi))))/GenericSign(s.levi);
    if e<>1 then  ChevieErr(s," => epsL/c mod Phi=",e,"\n");fi;
  fi;
  if CharNumbers(s)=false then return;fi;
  Hecke(s);
  e:=Eigenvalues(UnipotentCharacters(W)){s.charNumbers};
  if IsBound(s.Hecke) then 
    pred:=List(EigenAndDegHecke(s),x->x.eig)*
      Eigenvalues(UnipotentCharacters(s.levi))[s.cuspidal];
    if e<>pred and (not IsBound(s.delta) or s.delta=1 or IsCyclic(s.WGL)) then
      n:=CharNames(UnipotentCharacters(W)){s.charNumbers};
      ChevieErr(s," actual eigen differ from predicted eigen");
      c:=Set(Zip(pred,e,function(x,y)return x/y;end));
      if Length(c)=1 then ChevieErr("predicted=",c[1],"*actual");
      else
        CHEVIE.Check.EqTables(
                 rec(rowLabels:=["actual   "],columnLabels:=n,scalar:=[e]),
                 rec(rowLabels:=["predicted"],columnLabels:=n,scalar:=[pred]));
      fi;
    fi;
  fi;
end;

# test(WF[,d[,ad]]) or test(series)
CHEVIE.AddTest("Series",function(arg)local W,l,s,c,e,Ed,pred,n;
  W:=arg[1];
  if Length(arg)>=2 then l:=Filtered(ApplyFunc(Series,arg),s->s.levi<>s.spets);
  else l:=ProperSeries(W);
  fi;
  for s in l do CheckSerie(s);od;
  InfoChevie("\n   ");
  return l;
end,
IsSpetsial);

# ParameterExponents(W[,HC series no])
# check that the parameter exponents of the relative hecke algebras
# agree with unipotent degrees corresponding to cyclic
# relative groups of sub-series
CHEVIE.AddTest("ParameterExponents",
function(arg)local h,i,I,L,hh,ud,t,H,s,exp,W,ser,G,m,e;
  W:=arg[1];
  ser:=UnipotentCharacters(W).harishChandra;
  if Length(arg)=1 then
    for i in [1..Length(ser)] do
      CHEVIE.Test("ParameterExponents",W,i);
    od;
    InfoChevie("\n   ");
    return;
  fi;
  h:=ser[arg[2]];
  I:=Concatenation(List(h.relativeType,x->x.indices));
  if IsSpets(W) then I:=List(I,x->Cycle(W.phi,x));
  else 
#     InfoChevie("# ",IntListToString(h.levi),"\n");
#     CHEVIE.Testing(IntListToString(h.levi));
    G:=RelativeGroup(W,h.levi);
    if G.relativeIndices{Concatenation(List(G.type,t->t.indices))}<>I 
    then ChevieErr("indices computed=",
      IntListToString(G.relativeIndices{Concatenation(List(G.type,t->t.indices))}),
      " stored=",IntListToString(I),"\n");
    fi;
    I:=List(I,x->[x]);
  fi;
  for i in [1..Length(I)] do
    L:=ReflectionSubgroup(W,Concatenation(h.levi,I[i]));
    t:=ReflectionType(L);
    H:=ReflectionSubgroup(L,h.levi);
    InfoChevie("\n   # ParameterExponents from ",ReflectionName(H),":",I[i]);
    CHEVIE.Testing("from ",ReflectionName(H),":",I[i]);
    if t=false or (not IsBound(t[1].indices) and not IsBound(t[1].orbit))then
      ChevieErr("Levi ",Concatenation(h.levi,I[i])," could not be identified\n");
    else
      if not IsSpets(L) then L:=Spets(L);H:=SubSpets(L,h.levi);fi;
      hh:=FindSeriesInParent(h,H,L,UnipotentCharacters(L).harishChandra).ser;
      if Length(hh.charNumbers)<>2 then 
        s:=SeriesNC(L,H,Position(UnipotentCharacters(H).TeXCharNames,
               h.cuspidalName),1);
        SeriesOps.RelativeGroup(s);
        if CharNumbers(s)=false or SeriesOps.fill(s)=false then 
          Error("could not fill ",s);
        else 
          if IsInt(hh.parameterExponents[1]) then 
            exp:=[1..s.e]*0;exp[1]:=hh.parameterExponents[1];
          else exp:=hh.parameterExponents[1];
          fi;
          if IsBound(s.translation) then
          t:=Filtered([0,s.translation..s.e],function(d)local v;
            v:=1+List([1..s.e]+d,x->x mod s.e);
            return s.mC{v}=exp and s.charNumbers{v}=hh.charNumbers;end);
          else t:=[1];
          fi;
          if Length(t)=0 then Error("unexpected");fi;
#    Ok("decs=",t);
        fi;
      else
        ud:=UnipotentDegrees(L,X(Cyclotomics));
        ud:=ud[hh.charNumbers[1]]/ud[hh.charNumbers[2]];
        if Length(ud.coefficients)<>1 or ud.coefficients[1]<>1 then
          ChevieErr("not monomial");
        elif h.parameterExponents[i]<>ud.valuation then
          ChevieErr(ReflectionName(L),": wrong parameter ",
            h.parameterExponents[i],
            " instead of ",ud.valuation,"\n");
        fi;
#         Ok("=",ud.valuation);
      fi;
    fi;
  od;
end,
IsSpetsial);

# this was never called. Suppress?
# Check that the stored parameters of 1-HC series are correct
CHEVIE.AddTest("HCSeries",
function(WF)local sers,uc,ss,s,ul,p,para,stored_para,i;
  if not IsSpets(WF) then WF:=Spets(WF);fi;
  uc:=UnipotentCharacters(WF);
  sers:=uc.harishChandra;
  ss:=Filtered(CuspidalPairs(WF,1),x->SemisimpleRank(x[1])<SemisimpleRank(WF));
  for s in ss do
    s:=Series(WF,s[1],s[2],1);
    ul:=UnipotentCharacters(s.levi);
    p:=PositionProperty(sers,h->ul.TeXCharNames[s.cuspidal]=h.cuspidalName
      and Eigenvalues(ul)[s.cuspidal]=h.eigenvalue);
    para:=Hecke(s).parameter;
    para:=List(para,x->List(x,Mvp)); # should be unnecessary
    stored_para:=sers[p].parameterExponents;
    for i in s.WGL.generatingReflections do
      stored_para[i]:=stored_para[s.WGL.orbitRepresentative[i]];od;
    stored_para:=List([1..Length(para)],function(i)
      if IsRat(stored_para[i]) then return 
        Concatenation([stored_para[i]],[1..Length(para[i])-1]*0);
      else return stored_para[i];
      fi;end);
    para:=List(para,v->List(v,Degree));
    if para<>stored_para and para<>Reversed(stored_para) then 
      CHEVIE.Check.EqLists(s.charNumbers,sers[p].charNumbers,"charnum","stored");
      CHEVIE.Check.EqLists(para,stored_para,"para","stored para");
    fi;
  od;
end,
IsSpetsial);

# for a reflection group of rank r: Discriminant(G)
#          returns a list of linear factors as Mvps in x1,x2,...,xr
ReflectionDiscriminant:=function(W)local h,res,coroot,w,vars;
  res:=[];
  vars:=List([1..Length(ReflectionDegrees(W))],
		i->Mvp(Concatenation("x",String(i))));
  for h in HyperplaneOrbits(W) do
    coroot:=W.simpleCoroots[h.s];
    for w in List(Elements(ConjugacyClasses(W)[h.classno[1]]),x->
      RepresentativeOperation(W,W.generators[h.s],x)) do
      Append(res,List([1..h.e_s],i->MatXPerm(W,w^-1)*coroot));
    od;
  od;
  res:=List(res,x->Sum([1..Length(x)],i->Mvp(vars[i])*x[i]));
  return res;
end;

CHEVIE.AddTest("Discriminant",
function(W)local j,r,ii;
  r:=Product(ReflectionDiscriminant(W))*Mvp("x")^0;
  j:=Discriminant(W);
  if j=false then ChevieErr("not implemented\n");return;fi;
  ii:=Invariants(W);
  ii:=List(ii,x->ApplyFunc(x,List([1..W.rank],i->Mvp(SPrint("x",i)))));
  j:=ApplyFunc(j,ii);
  if r<>j*r.coeff[1]/j.coeff[1] then ChevieErr("disagrees\n");fi;
end,
x->not IsSpets(x) and Size(x)<1152); # F4 first painful client
