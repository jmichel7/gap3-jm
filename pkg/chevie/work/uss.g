########################################################################
# Until the end algos for USS
########################################################################
# returns the iterated transports of s along the cycling orbit of b
# (returns only the looping part)
CTransportOrbit:=function(b,s)local t,tl,b0,p,CTransport0,CTransport,
InitialFactor;
# CTransport(x,s)
# Given x in SSS(x), s simple  such that x^s in SSS(x) and c(x) the 
# element such that Cycling(x)=x^c(x), CTransport returns c(x)^-1*s*c(x^s)
CTransport0:=function(x,s) local M,R,b;
  M:=x.monoid;
  if Length(x.elm)=0 then return s;fi;
  b:=M.alpha2(M.RightComplementToDelta(x.elm[1]),M.DeltaAction(s,-1))[1];
  R:=M.LeftGcdSimples(M.AddToNormal(x.elm{[2..Length(x.elm)]},s)[1],b)[1];
  if R=M.identity then return M.delta;
  else return M.DeltaAction(R,-x.pd);
  fi;
end;

InitialFactor:=function(b)local M;M:=b.monoid;
 if Length(b.elm)=0 then return M.identity;fi;
 return M.DeltaAction(b.elm[1],-b.pd);
end;

# CTransport0 is fast but buggy...
CTransport:=function(x,s)
  return InitialFactor(x)^-1*s*InitialFactor(x^x.monoid.Elt([s]));
end;

  t:=s;tl:=[s];b0:=b;
  while true do
    repeat
#     if CTransport0(b,t)<>CTransport(b,t) then Error();fi;
      t:=CTransport(b,t);
      b:=Cycling(b)[1]; until b=b0;
    p:=Position(tl,t);
    if p<>false then return tl{[p..Length(tl)]};
    else Add(tl,t);
    fi;
  od;
end;

MainPullback:=function(y,s)local co,po,cycle,i,l,k,CPullback,CyclingOrbit;
# pullback is an almost-inverse of the transport
CPullback:=function(y,s)local M,b0,br,i,b;
  M:=y.monoid;
  if Length(y.elm)=0 then return s;fi;
  b0:=M.DeltaAction(M.RightLcmSimples(M.DeltaAction(
     M.RightComplementToDelta(y.elm[1]),-y.pd),s)[2],-1);
  br:=M.DeltaAction(s,y.pd);
  for i in [2..Length(y.elm)] do
    br:=M.RightLcmSimples(br,y.elm[i])[3];
  od;
  b:=M.RightLcmSimples(b0,br)[1];
  return MinConjugating.SS(y,b,function(arg)return arg[1];end);
end;

# Assumes b in SSS(b), otherwise may loop infinitely
CyclingOrbit:=function(b)local l;
  l:=[b];
  while true do b:=Cycling(b)[1];
    if b=l[1] then return l;else Add(l,b);fi;
  od;
end;

  co:=CyclingOrbit(y);
  po:=[];
  cycle:=function(s)
    for i in [Length(co),Length(co)-1..1] do s:=CPullback(co[i],s);od;
    return s;
  end;
  repeat Add(po,s); s:=cycle(s); until s in po;
  k:=Position(po,s)-1;l:=Length(po)-k;
  return po[1+k+((l-k) mod l)];
end;

########################################################################
# MinConjugating.USS(a,x)                                                #
#   This  program finds the minimal right multiple  m of the simple x  #
#   such that  a^m is in USS(a).                                        #
#   This m is a simple.                                                #
########################################################################
MinConjugating.USS:=function(y,b,F)local a,to,p,M;
  a:=MinConjugating.SS(y,b,F);
  to:=CTransportOrbit(y,a);
  M:=y.monoid;
  p:=PositionProperty(to,x->M.RightLcmSimples(b,x)[3]=M.identity);
  if p<>false then return to[p];fi;
  a:=MainPullback(y,b);
  to:=CTransportOrbit(y,a);
  p:=PositionProperty(to,x->M.RightLcmSimples(b,x)[3]=M.identity);
  if p<>false then return to[p];
  else Error("theory");
  fi;
end;
