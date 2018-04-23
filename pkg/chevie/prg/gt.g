# gt.g (C) Jean Michel 2016
#
# This file deals with class types for reductive grouops.
#
RationalUnipotentClasses:=function(WF,p)local t,u;
# Print("WF=",WF," p=",p,"\n");
  u:=UnipotentClasses(WF,p);
  t:=GreenTable(u,rec(classes:=true));
  return List([1..Length(t.locsys)],i->
      rec(card:=CycPol(t.cardClass[i]),
          class:=u.classes[t.locsys[i][1]],
          classno:=t.locsys[i][1],
          AuNo:=t.locsys[i][2]));
end;

ClassTypesOps:=OperationsRecord("ClassTypesOps");

ClassTypesOps.Format:=function(r,opts)local c,res,t,o,x,u,v,nc;
  res:=String(r);Append(res,"\n");
  nc:=function(p)local d;p:=Mvp(p);
    d:=Lcm(List(p.coeff,Denominator));
    if d=1 then return Format(p,opts);fi;
    return SPrint(BracketIfNeeded(Format(p*d,opts),"+-"),"/",d);
  end;
  if IsBound(opts.nrClasses) then NrConjugacyClasses(r);fi;
  if IsBound(opts.unip) then
    o:=rec(rowLabels:=[],rowsLabel:="Type",columnLabels:=[]);
    if IsBound(opts.nrClasses) then Add(o.columnLabels,"nrClasses");fi;
    Add(o.columnLabels,"u");
    Add(o.columnLabels,"Centralizer");
    t:=[]; 
    for x in r.ss do
      u:=RationalUnipotentClasses(x.CGs,r.p);
      for c in u do
        v:=[];
        if c.card=CycPol(1) then 
          if IsBound(opts.nrClasses) then Add(v,nc(x.nrClasses));fi;
          Add(o.rowLabels,ReflectionName(x.CGs,opts));
        else Add(o.rowLabels," ");if IsBound(opts.nrClasses) then Add(v,"");fi;
        fi;
        Add(v,Name(c.class,Inherit(rec(class:=c.AuNo),opts)));
        Add(v,Format(x.cent/c.card,Inherit(rec(vname:="q"),opts)));
        Add(t,v);
      od;
    od;
  else
    o:=rec(rowLabels:=List(r.ss,x->ReflectionName(x.CGs,opts)));
    o.rowsLabel:="Type";
    t:=[]; o.columnLabels:=[];
    if IsBound(opts.nrClasses) then 
      Add(t,List(NrConjugacyClasses(r),nc));
      Add(o.columnLabels,"nrClasses");
    fi;
    Add(t,List(r.ss,x->Format(x.cent,Inherit(rec(vname:="q"),opts)))); 
    Add(o.columnLabels,"Centralizer");
    t:=TransposedMat(t);
  fi;
  Append(res,FormatTable(t,Inherit(o,opts)));
  return res;
end;

ClassTypesOps.Value:=function(arg)local r;
  r:=ShallowCopy(arg[1]);
  r.specialized:=List([2,4..Length(arg[2])],i->arg[2]{[i-1,i]});
  r.ss:=List(r.ss,function(x)local res;
    res:=ShallowCopy(x);
    NrConjugacyClasses(r);
    res.nrClasses:=Value(x.nrClasses,arg[2]);
    return res;end);
  return r;
end;

ClassTypesOps.Display:=function(r,opts)
  Print(Format(r,opts));
end;

ClassTypesOps.String:=function(r)local res;
 res:=SPrint("ClassTypes(",r.spets);
 if r.p=0 then Append(res,",good characteristic)");
 else PrintToString(res,",char. ",r.p,")");
 fi;
 if IsBound(r.specialized) then
  PrintToString(res," ",Join(List(r.specialized,x->SPrint(x[1],"=",x[2]))," "));
 fi;
 return res;
end;

ClassTypesOps.Print:=function(r)Print(String(r));end;

# ClassTypes(W[,p])
ClassTypes:=function(arg)local W,WF,l,p;
  W:=arg[1];
  if Length(arg)=2 then p:=arg[2];else p:=0;fi;
  if IsSpets(W) then WF:=W;W:=Group(W);else WF:=Spets(W);fi;
  if Length(arg)=2 then p:=arg[2];else p:=0;fi;
  l:=SemisimpleCentralizerRepresentatives(Group(WF),p);
  l:=Concatenation(List(l,x->Twistings(WF,x)));
  return rec(p:=p,spets:=WF,ss:=List(l,x->rec(CGs:=x,
    cent:=CycPol(GenericOrder(x,Mvp("q"))),
    unip:=RationalUnipotentClasses(x,p))),
    operations:=ClassTypesOps);
end;

# returns the Poset of closed subsystems of the root system of W
ClosedSubsets:=function(W)local sum,closure,l,w,new,n,p,P,covered,f;
  if IsBound(W.closedsubsets) then return W.closedsubsets;fi;
  sum:=List([1..2*W.N],i->List([1..2*W.N],function(j)local p;
    p:=Position(W.roots,W.roots[i]+W.roots[j]);
    if p<>false then return p;else return 0;fi;end));
  closure:=function(l,new)local i,j,nnew;
    nnew:=new;
    repeat 
      new:=nnew;nnew:=[];
      for i in new do for j in l do
        if sum[i][j]<>0 then AddSet(nnew,sum[i][j]);fi;od;od;
      l:=Union(l,new);
      nnew:=Difference(nnew,l);
    until nnew=[];
    return l;
  end;
  l:=[[]];new:=[1];covered:=[[]];
  for w in new do
    for f in Difference([1..W.N],l[w]) do
      n:=closure(l[w],[f,f+W.N]);
      p:=Position(l,n);
      if p=false then Add(l,n);Add(covered,[w]);Add(new,Length(l));
      else AddSet(covered[p],w);fi;
    od;
  od;
  P:=Poset(covered);
  P.elements:=l;
  P.label:=function(P,i,opt)return Format(P.elements[i],opt);end;
  W.closedsubsets:=P;return P;
end;

# See Fleischmann-Janiszczak AAECC 1996 definition 2.1
ClassTypesOps.NrConjugacyClasses:=function(C)local HF,W,H,o,P,l,less,mu,n,i,r,b;
  W:=Group(C.spets);b:=Set(Factors(PermRootOps.BadNumber(W)));
  if Size(FundamentalGroup(W))>1 then
   Print("# Nr classes each type implemented only for simply connected groups");
   return;
  fi;
  for r in C.ss do
    if not IsBound(r.nrClasses) then 
    HF:=r.CGs;H:=Group(HF); P:=Copy(ClosedSubsets(W));
    o:=Filtered([1..Size(P)],i->OnSets(P.elements[i],HF.phi)=P.elements[i]);
    o:=Filtered(o,i->ForAll(H.rootInclusion{H.generatingReflections},j->
      j in P.elements[i]));
    P:=Restricted(P,o);P.elements:=P.elements{o};
    # here P poset of HF.phi-stable closed subsets containing roots(H)
    InfoChevie("# ",HF,"=>",P,"\c");
    l:=List(P.elements,x->Spets(ReflectionSubgroup(W,x),HF.phi));
    l:=List(l,function(RF)local res,d,R;R:=Group(RF);
      res:=Product(Filtered(ReflectionDegrees(RF),y->y[1]=1),
         p->Mvp("q")-p[2]);
      if R.semisimpleRank=0 then return res;fi;
      d:=SmithNormalFormMat(R.roots{R.generatingReflections}*R.simpleRoots);
      d:=Filtered(Flat(d),x->not x in [0,1,C.p]);
      return res*Product(d,function(i)
         if i=2 and 2 in b and C.p<>2 then return 2;fi;
         if C.p mod i=1 then return i;fi;
         return Mvp(SPrint("q_",i));
         end);
    end);
    less:=i->Difference(ListBlist([1..Length(l)],Incidence(P)[i]),[i]);
    o:=LinearExtension(P);mu:=[];mu[o[Size(P)]]:=1;
    for i in o{[Size(P)-1,Size(P)-2..1]} do mu[i]:=-Sum(mu{less(i)});od;
    n:=Stabilizer(W,Set(H.rootInclusion{H.generatingReflections}),OnSets);
    n:=mu*l/Size(Centralizer(n,HF.phi));
    InfoChevie("=>",n,":",Stime(),"\n");
    r.nrClasses:=n;
    fi;
  od;
  C.ss:=Filtered(C.ss,x->x.nrClasses<>0);
  return List(C.ss,x->x.nrClasses);
end;
