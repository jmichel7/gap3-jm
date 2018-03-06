clean:=function(l)local res;
  if not IsString(l[1]) then 
    if IsList(l[1][1]) then 
       if l[1][2]=0 then res:= List(l,x->IntListToString(x[1]));
       else res:= List(l,PartitionTupleToString);
       fi;
    else res:= List(l,IntListToString);fi;
  else res:= l;
  fi;
  return List(res,function(s)if s="1." then return "11";
    elif s=".1" then return "2";else return s;fi;end);
end;

cleanlab:=function(l)
  if IsString(l) then
    if l{[1..5]}<>"cusp(" then Error("lab");fi;
    l:=l{[6..Length(l)-1]};l:=Split(l);
  else l:=[IntListToString(l[1]),l[2]];
  fi;
  return l;
end;

findclass:=function(uc,n)local p;
  p:=PositionProperty(uc.classes,u->n=u.name
    or IsBound(u.mizuno) and u.mizuno=n);
  return p;
end;

# check a Luebeck record g again UnipotentClasses(g.group)
tr:=function(g)local n,r,f,ser,makep,ftop,q,W,uc,s,dsq;
  q:=X(Rationals);W:=g.group;
  ftop:=function(f)local tr,digits; digits:="01235"; tr:=[0,1,2,3,5];
    return tr[Position(digits,f[1])];end;
  makep:=function(y)if IsList(y) then 
    return q^y[2]*ValuePol(y[1],q);else return y*q^0;fi;end;
  ser:=function(r,ind)local p,cp,cp1,t,m,ci,cc,t,R;
   if r.levi=[] then
     p:=PermList(List(r.classtext,x->PositionClass(W,EltWord(W,x))));
     if p=false then Error("classtext");else Unbind(r.classtext);fi;
     cc:=ChevieClassInfo(W);
     if clean(r.classnames)<>Permuted(cc.classnames,p)
     then Error("classnames");else Unbind(r.classnames);fi;
     if r.classes<>Permuted(cc.classes,p)
     then Error("classes");else Unbind(r.classes);fi;
     cp:=PermListList(r.irreducibles,
	List(CharTable(W).irreducibles,x->Permuted(x,p)));
     if cp=false then Error("irred");else Unbind(r.irreducibles);fi;
     if clean(r.GreenLabels)<>Permuted(cc.classnames,p)
     then Error("GreenLabels");else Unbind(r.GreenLabels);fi;
     if not r.GreenLabel in[["A",0],["B",0]] then Error("GreenLabel");
     else Unbind(r.GreenLabel);fi;
     ci:=ChevieCharInfo(W);
     if clean(r.PtLabels)<>Permuted(ci.charnames,cp) 
     then CHEVIE.Check.EqLists(clean(r.PtLabels),Permuted(ci.charnames,cp),
       rec(na:="r.PtLabels",nb:= "R.charnames"));
        cp1:=PermListList(clean(r.PtLabels),ci.charnames);
     else Unbind(r.PtLabels);cp1:=cp;fi;
     t:=ICCTable(uc,1,q);
     if r.sqpowers<>2*Permuted(t.dimBu,cp1) 
     then CHEVIE.Check.EqLists(r.sqpowers,2*Permuted(t.dimBu,cp1),
           rec(na:="r.sqpowers",nb:="t.dimBu"));
     else Unbind(r.sqpowers);fi;
     if r.Ptransposed=TransposedMat(OnMatrices(t.scalar,cp1))
     then Unbind(r.Ptransposed);else 
      CHEVIE.Check.EqLists(r.Ptransposed,TransposedMat(OnMatrices(t.scalar,cp1)),
         rec(na:="r.Ptransposed",nb:="t.scalar"));
     fi;
     Unbind(r.levi);
     r.Lambda:=[];
     for m in r.Lambdas do
       if IsDiagonalMat(m) then Append(r.Lambda,List(DiagonalOfMat(m),x->[[x]]));
       else Add(r.Lambda,m);
       fi;
     od;
     Unbind(r.Lambdas);
     m:=OnMatrices(t.L,cp);
     m:=List(DecomposedMat(m),x->m{x}{x});
     if List(m,Length)<>List(r.Lambda,Length) then Error("length");fi;
     m:=Value(m,Mvp("x"));r.Lambda:=Value(r.Lambda,Mvp("x"));
     m:=Zip(m,r.Lambda,function(x,y)return y/x;end);
     if not ForAll(m,IsDiagonalMat) then Error("diagonal");fi;
     m:=Concatenation(List(m,DiagonalOfMat));
 #   p:=Set(Zip(m,Permuted(t.dimBu,cp),function(x,z)return
 #   x*Mvp("x")^(2*z);end));
     p:=Set(m);
     if Length(p)<>1 then Error("non constant");fi;
#    if p[1]<>Mvp("x")^W.N then Error("mal compris");fi;
     if p[1]<>1 then Error("mal compris");fi;
     Unbind(r.Lambda);
     if clean(r.YiLabels)<>Permuted(ci.charnames,cp) 
      and (not IsBound(ci.spaltenstein) or 
	   r.YiLabels<>Permuted(ci.spaltenstein,cp)) then Error("YiLabels");
     else Unbind(r.YiLabels);fi;
   elif r.levi=[1..W.semisimpleRank] then
     if r.classtext<>[[]] then Error("class");else Unbind(r.classtext);fi;
     if clean(r.classnames)<>["1"] then Error("cname");else Unbind(r.classnames);fi;
     if r.classes<>[1] then Error();else Unbind(r.classes);fi;
     if r.irreducibles<>[[1]] then Error();else Unbind(r.irreducibles);fi;
     if Length(r.GreenLabels)<>1 or clean(r.GreenLabels)<>["1"] 
     then Error("GreenLabels");else Unbind(r.GreenLabels);fi;
     t:=ReflectionType(W)[1];
     if r.GreenLabel<>[t.series,t.rank] then Error("GreenLabel");
     else Unbind(r.GreenLabel);fi;
     if r.PtLabels<>[[1]] then Error("charnames");else Unbind(r.PtLabels);fi;
     if r.Ptransposed<>[[q^0]] then Error("Ptransposed ",f);else Unbind(r.Ptransposed);fi;
     if GenericOrder(W,q)/r.Lambdas[1][1][1]<>q^r.sqpowers[1] then
       Error("cuspLambda");else Unbind(r.Lambdas);fi;
     r.dsq:=(r.sqpowers[1]-Length(r.levi))/2;
     Unbind(r.levi); Unbind(r.sqpowers);
     if Length(r.YiLabels)<>1 then Error("YiLabels");fi;
     r.lab:=cleanlab(r.YiLabels[1]);Unbind(r.YiLabels);
     p:=findclass(uc,r.lab[1]);
     if p=false then Error("class ",r.lab[1]," not found");fi;
     if uc.classes[p].dimBu=r.dsq then Unbind(r.dsq);else Error("dsq");fi;
     p:=PositionsProperty(uc.springerSeries,x->Length(x.locsys)=1 and x.locsys[1][1]=p);
     if Length(p)=1 and Length(r.lab)=1 then Unbind(r.lab);
     elif Length(p)=0 then Error("cuspidal series ",r.lab," not found");fi;
   else
     R:=RelativeGroup(W,r.levi);
     cc:=Copy(ChevieClassInfo(R));
     if ReflectionName(R)="G2" then cc.classtext:=3-cc.classtext;fi;
     p:=PermListList(List(r.classtext,x->PositionClass(W,EltWord(W,x))),
       List(cc.classtext,function(x)if x=[] then return PositionClass(W,());
       else return PositionClass(W,Product(R.parentMap{x}));fi;end));
     if p=false then Error("classtext"); else Unbind(r.classtext); fi;
     if clean(r.classnames)<>Permuted(cc.classnames,p) then 
       CHEVIE.Check.EqLists(clean(r.classnames),Permuted(cc.classnames,p),
         rec(na:="r.classnames",nb:="R.classnames"));
     else Unbind(r.classnames);fi;
     if clean(r.GreenLabels)<>Permuted(cc.classnames,p)
     then CHEVIE.Check.EqLists(clean(r.GreenLabels),Permuted(cc.classnames,p),
       rec(na:="r.GreenLabels",nb:= "R.classnames"));
     else Unbind(r.GreenLabels);fi;
     t:=List(ReflectionType(ReflectionSubgroup(W,r.levi)),x->[x.series,x.rank]);
     if not r.GreenLabel in [t,t[1]] then Error("GreenLabel");
     else Unbind(r.GreenLabel);fi;
     if r.classes<>Permuted(CharTable(R).classes,p)
     then Error("classes");else Unbind(r.classes);fi;
     cp:=PermListList(r.irreducibles,
	List(CharTable(R).irreducibles,x->Permuted(x,p)));
     if cp=false then Error("irred");else Unbind(r.irreducibles);fi;
     ci:=CharNames(R);
     if ReflectionName(R)="G2" then ci:=Permuted(ci,(3,4));fi;
     if clean(r.PtLabels)<>Permuted(ci,cp) 
     then CHEVIE.Check.EqLists(clean(r.PtLabels),Permuted(ci,cp),
        rec(na:="r.PtLabels",nb:="R.charnames"));
       cp1:=PermListList(clean(r.PtLabels),ci);
     else Unbind(r.PtLabels);cp1:=cp;fi;
     t:=ICCTable(uc,ind,q);
     if r.Ptransposed=TransposedMat(OnMatrices(t.scalar,cp1)) then Unbind(r.Ptransposed);else 
     CHEVIE.Check.EqLists(r.Ptransposed,TransposedMat(OnMatrices(t.scalar,cp1)),
       rec(na:="r.Ptransposed",nb:="t.scalar"));
     fi;
     dsq:=r.sqpowers-2*Permuted(t.dimBu,cp1);
     dsq:=Set(dsq);
     if Length(dsq)<>1 then Error("dsq");fi;
     dsq:=dsq[1];Print("dsq=",dsq,"\n");
     r.Lambda:=[];
     for m in r.Lambdas do
       if IsDiagonalMat(m) then Append(r.Lambda,List(DiagonalOfMat(m),x->[[x]]));
       else Add(r.Lambda,m);
       fi;
     od;
     Unbind(r.Lambdas);
     m:=OnMatrices(t.L,cp1);
     m:=List(DecomposedMat(m),x->m{x}{x});
     if List(m,Length)<>List(r.Lambda,Length) then Error("length");fi;
     m:=Value(m,Mvp("x"));r.Lambda:=Value(r.Lambda,Mvp("x"));
     m:=Zip(m,r.Lambda,function(x,y)return y/x;end);
     if not ForAll(m,IsDiagonalMat) then Error("diagonal");fi;
     m:=Concatenation(List(m,DiagonalOfMat));
#    p:=Set(Zip(m,Permuted(t.dimBu,cp1),function(x,z)return
#    x*Mvp("x")^(2*z);end));
     p:=Set(m);
     if Length(p)<>1 then Error("non constant");fi;
 #   p:=p[1]*GenericOrder(R,Mvp("x"))/
 #     GenericOrder(W,Mvp("x"))/Mvp("x")^R.N*Mvp("x")^dsq;
     if p[1]<>1 then Error("mal compris");fi;
     Unbind(r.Lambda);
     Unbind(r.sqpowers);
     r.dsq:=(dsq-Length(r.levi))/2;
   fi;
  end;
  for f in Difference(RecFields(g),["group"]) do
    r:=g.(f);
    uc:=UnipotentClasses(W,ftop(f));
    for n in uc.classes do 
      if not IsBound(n.dimBu) then n.dimBu:=W.N-DimUnipotentClass(W,n.dynkin)/2;fi;
    od;
    n:=List([1..Length(r.Ls)],i->rec(
      classtext:=List(r.NLLRepresentatives[i],x->x[1]),
      classnames:=List(r.NLLRepresentatives[i],x->x[2]),
      classes:=r.NLLclasses[i],
      irreducibles:=r.NLLirreducibles[i],
      levi:=r.Ls[i],
      GreenLabels:=List(r.GreenLabels[i],x->x[2]),
      GreenLabel:=Set(List(r.GreenLabels[i],x->x[1])),
      Lambdas:=List(r.Lambdas[i],m->List(m,l->List(l,makep))),
      PtLabels:=r.PtLabels[i],
      Ptransposed:=List(r.Ptransposed[i],l->List(l,makep)),
      YiLabels:=r.YiLabels[i],
      sqpowers:=r.sqpowers[i]
    ));
    Unbind(r.NLLRepresentatives); Unbind(r.NLLclasses);
    Unbind(r.NLLirreducibles); Unbind(r.Ls);
    Unbind(r.GreenLabels); Unbind(r.Lambdas);
    Unbind(r.PtLabels); Unbind(r.Ptransposed);
    Unbind(r.YiLabels); Unbind(r.sqpowers);
    if Length(RecFields(r))>0 then Error("fields remaining");fi;
    g.(f):=n;
    Error();
    for r in n do
      if Length(r.GreenLabel)>1 then Error("GreenLabel too long");
      else r.GreenLabel:=r.GreenLabel[1];fi;
      s:=PositionsProperty(uc.springerSeries,x->x.levi=r.levi);
      Print("series number=",s,"\n");
      ser(r,s[1]);
    od;
    Print("done ",f,"\n");
  od;
end;
