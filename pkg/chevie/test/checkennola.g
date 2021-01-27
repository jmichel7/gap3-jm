
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

# i^n where n fraction i AsRootOfUnity computed in 'right' way
# so it agrees with GetRoot
PowerRoot:=function(i,n)local d,j,n1,k;
  d:=Denominator(i);j:=1;
  n1:=Denominator(n);
  if n1=1 then return Mod1(n*i);fi;
  repeat k:=Gcd(n1,d);n1:=n1/k;j:=j*k;until k=1;
  return Mod1((Numerator(i)*Numerator(n)*GcdRepresentation(n1,d)[1])/j/d);
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
