#############################################################################
##
#A  sscoset.g           CHEVIE library      Fran\c cois Digne and Jean Michel
##
#Y  Copyright (C) 2017  University  Paris VII and Universite de Picardie
##
##  This file contains functions dealing with semisimple elements of algebraic
##  cosets.
##
## An algebraic coset G.\sigma where \sigma quasi-central is represented
## by a CoxeterCoset WF, where WF.F0Mat is the action of \sigma on X(T).
##
## A finite order quasi-semisimple element t\sigma is represented as t, an 
## element of Y(T)\otimes Q/Z=(Q/Z)^n, that is a list of length n of rationals r
## such that 0<=r<1.
##
#####################################################################

# IsSpecial(WF,c) c is an orbit of WF.phi on WF.rootInclusion
# return true iff c is special in the sense of Digne-Michel
IsSpecial:=function(WF,c)local p,W;
  if Length(c) mod 2<>0 then return false;fi;
  W:=Group(WF);
  return W.roots[W.rootRestriction[c[1]]]+
         W.roots[W.rootRestriction[c[1+Length(c)/2]]] in W.roots;
end;

# Computes X_\sigma X^\sigma, Y_\sigma, Y^\sigma, R(\sigma) of ss
# see 1.1 to 1.7 of \cite{ss}
RelativeDatum:=function(WF)local n,W,Phis,cPhis,cc;
  if not IsBound(WF.Rs) then
    W:=Group(WF);
    n:=OrderMat(WF.F0Mat); # matrix of sigma on X
    WF.pi:=Sum([0..n-1],i->WF.F0Mat^i)/n;
    WF.X_s:=BaseIntMat(n*WF.pi)/n; # basis of X_\sigma
    WF.Y_s:=BaseIntMat(n*TransposedMat(WF.pi))/n; # basis of Y_\sigma
    WF.Xs:=(WF.X_s*TransposedMat(WF.Y_s))^-1*WF.X_s; # basis of X^\sigma
    WF.Ys:=(WF.Y_s*TransposedMat(WF.X_s))^-1*WF.Y_s; # basis of Y^\sigma
    cc:=Cycles(WF.phi,W.rootInclusion{W.generatingReflections});
    Phis:=List(cc,function(c)local res;
      res:=Sum(W.roots{W.rootRestriction{c}})*W.simpleRoots;
      if IsSpecial(WF,c) then res:=2*res;fi; return res;end);
    cPhis:=List(cc,c->Sum(W.coroots{W.rootRestriction{c}}*W.simpleCoroots)/
      Length(c));
    WF.Rs:=CoxeterGroup(List(Phis,x->SolutionMat(WF.Xs,x)),
                       List(cPhis,x->SolutionMat(WF.Y_s,x)));
  fi;
  return WF.Rs; # rootdatum R(\sigma)
end;

# computes the centralizer of t\sigma as an ExtendedReflectionGroup
CoxeterCosetOps.Centralizer:=function(WF,t)
  local W,Rs,cRs,o,t,C,i,refC,p,good,labels;
  W:=Group(WF);
  if not IsBound(WF.Cso) then # compute constants C_\sigma,\alpha
    WF.Cso:=[1..W.N]*0+1;
    for o in ReflectionType(WF) do
      if OrderPerm(o.twist)=2 and 
         o.orbit[1].series="A" and o.orbit[1].rank mod 2=0 then
        for p in o.orbit do WF.Cso{p.indices{p.rank/2+[0,1]}}:=[-1,-1];od;
      fi;
    od;
    C:=function(i)local p,j;
      for j in [1..W.semisimpleRank] do
        if W.roots[i][j]>0 then p:=Position(W.roots,W.roots[i]-W.roots[j]);
          if p<>false then return WF.Cso[j]*WF.Cso[p];fi;
        fi;
      od;
    end;
    for i in [W.semisimpleRank+1..W.N] do WF.Cso[i]:=C(i);od;
    Append(WF.Cso,WF.Cso);
  fi;
  RelativeDatum(WF);
  refC:=Centralizer(WF.Rs,
    SemisimpleElement(WF.Rs,SolutionMat(WF.Y_s,t.v*TransposedMat(WF.pi))));
  Rs:=List(Cycles(WF.phi,W.rootInclusion{[1..W.N]}),function(c)local res;
    res:=[c,Sum(W.roots{W.rootRestriction{c}}*W.simpleRoots)/Length(c),
          Sum(W.coroots{W.rootRestriction{c}})*W.simpleCoroots];
    if IsSpecial(WF,c) then res[3]:=2*res[3];fi;
    return res;end);
  labels:=List(Cycles(WF.phi,W.rootInclusion{[1..W.N]}),IntListToString);
  if t.additive  then good:=List(Rs,p->Mod1(Sum(p[1],
      i->t^Parent(W).roots[i])+AsRootOfUnity(WF.Cso[p[1][1]]))=0);
  else good:=List(Rs,p->Product(p[1],
      i->t^Parent(W).roots[i])*WF.Cso[p[1][1]]=t.v[1]^0);
  fi;
  Rs:=ListBlist(Rs,good); labels:=ListBlist(labels,good);
  cRs:=List(Rs,x->x[3]); cRs:=List(cRs,x->SolutionMat(WF.Ys,x));
  cRs:=Filtered(cRs,x->not x in List(Cartesian(cRs,cRs),Sum));
  Rs:=List(Rs,x->x[2]); Rs:=List(Rs,x->SolutionMat(WF.X_s,x));
  good:=List(Rs,x->not x in List(Cartesian(Rs,Rs),Sum));
  Rs:=ListBlist(Rs,good);
  labels:=ListBlist(labels,good);
  if Length(Rs)>0 then C:=CoxeterGroup(Rs,cRs);
  else C:=Torus(Length(WF.Xs));
  fi;
# C.reflectionsLabels:=labels;
  C.operations.ReflectionFromName:=function(W,x)
    return Position(W.rootInclusion,x);end;
  p:=List(WF.Xs,x->SolutionMat(WF.X_s,x));
  # transfer matrix on X^\sigma to X_\sigma
  if refC.F0s=[] then return ExtendedReflectionGroup(C);fi;
  return ExtendedReflectionGroup(C,ApplyFunc(Group,List(refC.F0s,x->x^p)));
end;

CoxeterCosetOps.Centralizer:=function(WF,t)
  local W,Rs,cRs,o,t,C,i,refC,p,good,labels;
  if  not IsSemisimpleElement(t) then
    Error(t, " must be a element semisimple element");
  fi;
  W:=Group(WF);
  if not IsBound(WF.Cso) then # compute constants C_\sigma,\alpha
    WF.Cso:=[1..W.N]*0+1;
    for o in ReflectionType(WF) do
      if OrderPerm(o.twist)=2 and 
         o.orbit[1].series="A" and o.orbit[1].rank mod 2=0 then
        for p in o.orbit do WF.Cso{p.indices{p.rank/2+[0,1]}}:=[-1,-1];od;
      fi;
    od;
    C:=function(i)local p,j;
      for j in [1..W.semisimpleRank] do
        if W.roots[i][j]>0 then p:=Position(W.roots,W.roots[i]-W.roots[j]);
          if p<>false then return WF.Cso[j]*WF.Cso[p];fi;
        fi;
      od;
    end;
    for i in [W.semisimpleRank+1..W.N] do WF.Cso[i]:=C(i);od;
    Append(WF.Cso,WF.Cso);
  fi;
  RelativeDatum(WF);
  refC:=Centralizer(WF.Rs,
    SemisimpleElement(WF.Rs,SolutionMat(WF.Y_s,t.v*TransposedMat(WF.pi))));
  Rs:=List(Cycles(WF.phi,W.rootInclusion{[1..W.N]}),function(c)local res;
    res:=[c,Sum(W.roots{W.rootRestriction{c}}*W.simpleRoots)/Length(c),
          Sum(W.coroots{W.rootRestriction{c}})*W.simpleCoroots];
    if IsSpecial(WF,c) then res[3]:=2*res[3];fi;
    return res;end);
  labels:=List(Cycles(WF.phi,W.rootInclusion{[1..W.N]}),IntListToString);
  if t.additive  then good:=List(Rs,p->Mod1(Sum(p[1],
      i->t^Parent(W).roots[i])+AsRootOfUnity(WF.Cso[p[1][1]]))=0);
  else good:=List(Rs,p->Product(p[1],
      i->t^Parent(W).roots[i])*WF.Cso[p[1][1]]=t.v[1]^0);
  fi;
  Rs:=ListBlist(Rs,good); labels:=ListBlist(labels,good);
  cRs:=List(Rs,x->x[3]); cRs:=List(cRs,x->SolutionMat(WF.Ys,x));
  cRs:=Filtered(cRs,x->not x in List(Cartesian(cRs,cRs),Sum));
  Rs:=List(Rs,x->x[2]); Rs:=List(Rs,x->SolutionMat(WF.X_s,x));
  good:=List(Rs,x->not x in List(Cartesian(Rs,Rs),Sum));
  Rs:=ListBlist(Rs,good);
  labels:=ListBlist(labels,good);
  if Length(Rs)>0 then C:=CoxeterGroup(Rs,cRs);
  else C:=Torus(Length(WF.Xs));
  fi;
# C.reflectionsLabels:=labels;
  C.operations.ReflectionFromName:=function(W,x)
    return Position(W.rootInclusion,x);end;
  p:=List(WF.Xs,x->SolutionMat(WF.X_s,x));
  # transfer matrix on X^\sigma to X_\sigma
  if refC.F0s=[] then return ExtendedReflectionGroup(C);fi;
  return ExtendedReflectionGroup(C,ApplyFunc(Group,List(refC.F0s,x->x^p)));
end;
# returns representatives of quasi-isolated classes of G.\sigma
CoxeterCosetOps.QuasiIsolatedRepresentatives:=function(arg)local WF,p;
  WF:=arg[1]; if Length(arg)=2 then p:=arg[2];else p:=0;fi;
  return List(QuasiIsolatedRepresentatives(RelativeDatum(WF),p),
    x->SemisimpleElement(Group(WF),x.v*WF.Y_s));
end;

# whether t\sigma is isolated
CoxeterCosetOps.IsIsolated:=function(WF,t)
  t:=SemisimpleElement(WF.Rs,SolutionMat(WF.Y_s,t.v*TransposedMat(WF.pi)));
  return IsIsolated(WF.Rs,t);
end;
