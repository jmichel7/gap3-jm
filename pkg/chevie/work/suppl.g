# tests which cannot be in chevie/test since they use undocumented
# functions

CHEVIE.AddTest("Discriminant",
function(W)local j,r,ii;
  r:=Product(ReflectionDiscriminant(W))*Mvp("x")^0;
  j:=InvariantsDiscriminant(W);
  if j=false then ChevieErr("not implemented\n");return;fi;
  ii:=Invariants(W);
  ii:=List(ii,x->ApplyFunc(x,List([1..W.rank],i->Mvp(SPrint("x",i)))));
  j:=ApplyFunc(j,ii);
  if r<>j*r.coeff[1]/j.coeff[1] then ChevieErr("disagrees\n");fi;
end,
x->not IsSpets(x) and Size(x)<1152); # F4 first painful client

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

# According to Michel (principal series) and Gunter (other cases)
# the eigenvalue of Ennola_E(z)^k(rho_chi) is  
#  EigenEnnola(W,i,k)/Ennolascalar(rho,E(z)^k)
# where z=|ZW| and rho_chi is i-th unip. char.
CHEVIE.AddTest("Ennola",
function(W)
  if Length(Ennola(W))>1 then ChevieErr("Ambiguity in Ennola\n");fi;
end,
W->IsSpetsial(W) and not IsSpets(W) and OrderCenter(W)>1);

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
        s:=Series(L,H,Position(UnipotentCharacters(H).TeXCharNames,
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

IsClassical:=W->ForAll(ReflectionType(W),function(n)
  if IsBound(n.series) then return n.series in ["A","B","C","D","G","I"] 
    or n.series="ST" and IsBound(n.p);
  elif IsBound(n.orbit) then return n.orbit[1].series in ["A","B","C","D","I"]
    and OrderPerm(n.twist)<>3;
  fi;end);

CHEVIE.AddTest("UniDegClassical",
function(W)local cs,ud,i,vud;
  ud:=CycPolUnipotentDegrees(W);
  cs:=List(UnipotentCharacters(W).charSymbols,x->x[1]);
  vud:=List(cs,CycPolGenericDegreeSymbol);
  for i in [1..Length(cs)] do cmpcycpol(vud[i],ud[i],
     SPrint("Deg",StringSymbol(cs[i])),"ud");
  od;
end,
W->IsClassical(W) and (Length(ReflectionName(W))=1 or 
   ReflectionName(W){[1,2]}<>"2A"));

CHEVIE.AddTest("FakeDegClassical",
function(W)local uc,fd,i,vfd,n,cs;
  n:=ReflectionName(W);
  uc:=UnipotentCharacters(W);
  cs:=List(uc.almostCharSymbols{uc.almostHarishChandra[1].charNumbers},x->x[1]);
  fd:=List(FakeDegrees(W,X(Cyclotomics)),CycPol);
  if Length(n)>1 and n{[1,2]} in ["2D","2I","2B","2G"] then 
       vfd:=List(cs,x->CycPolFakeDegreeSymbol(x,1));
  else vfd:=List(cs,CycPolFakeDegreeSymbol);
  fi;
  for i in [1..Length(cs)] do cmpcycpol(vfd[i],fd[i],
    SPrint("Deg",StringSymbol(cs[i])),"fd");
  od;
end,
W->IsClassical(W) and (Length(ReflectionName(W))=1 or 
   ReflectionName(W){[1,2]}<>"2A"));

CHEVIE.relativeSeries:=true;
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
  if Length(arg)>=2 then 
       l:=Filtered(ApplyFunc(CuspidalSeries,arg),s->s.levi<>s.spets);
  else l:=AllProperSeries(W);
  fi;
  for s in l do CheckSerie(s);od;
  InfoChevie("\n   ");
  return l;
end,
IsSpetsial);
