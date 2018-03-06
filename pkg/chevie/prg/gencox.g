#############################################################################
##
#A  gencox.g              CHEVIE library                          Jean Michel
##
##
#Y  Copyright (C) 1992 - 1999  Lehrstuhl D f\"ur Mathematik, RWTH Aachen, IWR
#Y  der Universit\"at Heidelberg, University of St. Andrews, and   University 
#Y  Paris VII.
##
##  This file contains functions for various ways of defining some
##  abstract Coxeter groups.
##

GenCoxOps:=OperationsRecord("GenCoxOps",MatGroupOps);
Inherit(GenCoxOps,AbsCoxOps);
# inherit from MatGroupOps and AbsCoxOps

GenCoxOps.String:=W->W.name;

GenCoxOps.Print:=function(W)Print(String(W));end;

GenCoxOps.IsLeftDescending:=function(W,w,i)local s; s:=Sum(w[i]);
  if IsInt(s) then return s<0;else return evalf(s)<0;fi;
end;

#############################################################################
##
#F  CoxeterGroupByCoxeterMatrix(<M>) . . . . Makes the Coxeter group defined
#F    by Coxeter matrix M. An infinite entry is represented by 0.
##
CoxeterGroupByCoxeterMatrix:=function(M)local tmp,W,i,j;
  tmp:=rec();
  tmp.nbGeneratingReflections:=Length(M);
  tmp.generatingReflections:=[1..tmp.nbGeneratingReflections];
  tmp.identity:=M^0;
  tmp.coxeterMat:=M;
  tmp.cartan:=CartanMatFromCoxeterMatrix(M);
  tmp.reflections:=List(tmp.generatingReflections,i->M^0);
  M:=M^0-tmp.cartan;
  for i in tmp.generatingReflections do 
    tmp.reflections[i]{tmp.generatingReflections}[i]:=M[i]; 
  od;
  W:=ApplyFunc(Group,tmp.reflections);
  Inherit(W,tmp);
  W.operations:=GenCoxOps;
  W.name:=SPrint("CoxeterGroupByCoxeterMatrix(",FormatGAP(W.coxeterMat),")");
  AbsCoxOps.CompleteCoxeterGroupRecord(W);
  return W;
end;

#############################################################################
##
#F  CoxeterGroupByCartanMatrix(<M>) . . . . Makes the Coxeter group defined
#F    by Cartan matrix M
##
CoxeterGroupByCartanMatrix:=function(M)local tmp,W,k;
  if ForAll(FiniteCoxeterTypeFromCartanMat(M),x->x<>false) then 
     return CoxeterGroup(M);fi;
  tmp:=rec();
  tmp.nbGeneratingReflections:=Length(M);
  tmp.generatingReflections:=[1..tmp.nbGeneratingReflections];
  tmp.identity:=M^0;
  tmp.cartan:=M;
  tmp.reflections:=List(tmp.generatingReflections,i->M^0);
  M:=M^0-M;
  for k in tmp.generatingReflections do 
    tmp.reflections[k]{tmp.generatingReflections}[k]:=M[k]; 
  od;
  W:=ApplyFunc(Group,tmp.reflections);
  Inherit(W,tmp);
  W.operations:=GenCoxOps;
  W.name:=SPrint("CoxeterGroupByCartanMatrix(",FormatGAP(W.cartan),")");
  AbsCoxOps.CompleteCoxeterGroupRecord(W);
  return W;
end;

#############################################################################
#
# Example I: Symmetric groups as Coxeter groups
#
CoxeterGroupSymmetricGroup:=function(n)local W;
  W:=SymmetricGroup(n);
  W.nbGeneratingReflections:=n-1;
  W.reflections:=Concatenation(List([1..n],i->List([i+1..n],j->(j-i,j))));
  W.operations.IsLeftDescending:=function(W,w,i) return i^w>(i+1)^w;end;
  AbsCoxOps.CompleteCoxeterGroupRecord(W);
  W.operations.FirstLeftDescending:=function(W,w)local i;
    for i in [1..n-1] do if i^w>(i+1)^w then return i;fi;od;return false;end;
  W.operations.ReflectionDegrees:=W->[2..n];
  W.operations.Reflections:=W->W.reflections;
  W.operations.ReflectionCharValue:=
                         function(W,w)return Number([1..n],i->i=i^w)-1;end;
  W.operations.ReflectionLength:=function(W,w)local l,i,cyc,mark;l:=n;mark:=[];
    for i in [1..n] do if not IsBound(mark[i]) then 
      cyc:=CyclePermInt(w,i); mark{cyc}:=cyc; l:=l-1;
    fi; od;
    return l;
  end;
  W.operations.String:=W->SPrint("CoxeterGroupSymmetricGroup(",n,")");
  W.name:=String(W);
  # existence of next saves a factor 2 in Braid multiplications
  W.operations.LeftGcdSimples:=function(arg)local res,i,s,try;res:=();
    try:=function(i)local j,b;j:=i+1;
      for b in arg do if i^b<j^b then return;fi;od;
      s:=s*(i,j);
    end;
    while true do
      s:=();for i in [1..n-1] do try(i);od;
      if s=() then return Concatenation([res],arg);fi;
      res:=res*s;arg:=s^-1*arg;
    od;
  end;
  W.operations.DualBraidMonoid:=PermRootOps.DualBraidMonoid;
  W.operations.CoxeterLength:=function(W,w)local i ,j, length;
    length:=0;
    for j in [1..n] do for i in [1..j-1] do
        if i^w>j^w then length:=length+1; fi;
    od; od;
    return length;
  end;
  return W;
end;

#############################################################################
#
# Example II: HyperOctaedral groups as Coxeter groups
#
CoxeterGroupHyperoctaedralGroup:=function(n)local W,i,f,conj;
  conj:=i->2*n+1-i;
  f:=[(n,conj(n))];for i in [n,n-1..2] do Add(f,(i,i-1)(conj(i),conj(i-1)));od;
  W:=Group(f,());
  W.reflections:=f;
  W.operations.IsLeftDescending:=function(W,w,i) return (n+1-i)^w>(n+2-i)^w;end;
  AbsCoxOps.CompleteCoxeterGroupRecord(W);
  W.operations.Reflections:=function(W)local i;
    if Length(W.reflections)>=n^2 then return W.reflections{[1..n^2]};fi;
    W.reflections:=Concatenation(W.reflections{[1..n]},
                                List([1..n-1],i->(i,conj(i))));
    for i in [1..n] do  Append(W.reflections,
      List(Concatenation([i+2..n],List([1..i-1],conj)),j->(i,j)(conj(i),conj(j))));
    od;
    return W.reflections{[1..n^2]};
  end;
  W.operations.ReflectionLength:=function(W,w)local c,c1;c:=Cycles(w);
    c1:=Filtered(c,x->Set(x)=Set(List(x,conj)));
    return Sum(c1,Length)/2+Sum(Difference(c,c1),x->Length(x)-1)/2;
  end;
  W.operations.DualBraidMonoid:=PermRootOps.DualBraidMonoid;
  return W;
end;
