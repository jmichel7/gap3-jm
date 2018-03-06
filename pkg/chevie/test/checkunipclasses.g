# Check various consistencies for UnipotentClasses(W[,p])
CHEVIE.AddTest("UnipotentClasses",
function(arg)
  local W,uc,b,u,a,ad,bu,s,nid,du,i,j,k,pl,cl,order,t,l,p,L,w,bc,name;
  pl:=function(u,i)return SPrint(Position(uc.classes,u),".",i,"=",u.name,".",
       CharNames(u.Au)[i]);end;
  W:=arg[1];
  if Length(arg)=1 then uc:=UnipotentClasses(W);
  else uc:=UnipotentClasses(W,arg[2]);fi;
  s:=uc.springerSeries[1];
  b:=LowestPowerFakeDegrees(W);du:=[];
  for u in uc.classes do
    nid:=PositionId(u.Au);
    if IsBound(u.dimBu) then a:=u.dimBu;fi;
    if IsBound(u.dynkin) and W.N<>0 then 
      ad:=W.N-Number(W.roots*u.dynkin,x->not x in [0,1])/2;
      if IsBound(u.dimBu) and a<>ad then
       ChevieErr(pl(u,nid)," dimBu=>",a," and from Dynkin=>",ad,"\n");
      else a:=ad;
      fi;
    fi;
    Add(du,a);
    i:=PositionProperty(s.locsys,y->y[1]=Position(uc.classes,u) and y[2]=nid);
    bu:=b[i];
    if a=bu then InfoChevie(".\c");
    else ChevieErr(pl(u,nid)," dimBu=>",a," and Springer[1.",i,"]=>",bu,"\n");
    fi;
    if IsBound(u.dimBu) and a<>u.dimBu then
      ChevieErr(pl(u,nid)," dynkin=>",a," and dimBu=",u.dimBu,"\n");
    fi;
    for i in [1..NrConjugacyClasses(u.Au)] do a:=[];
      for j in [1..Length(uc.springerSeries)] do
        for k in [1..Length(uc.springerSeries[j].locsys)] do
	  if [Position(uc.classes,u),i]=uc.springerSeries[j].locsys[k] then
	    Add(a,[j,k]);fi;
	od;
      od;
      if Length(a)=0 then ChevieErr(pl(u,i)," not accounted for\n");
      elif Length(a)>1 then ChevieErr(pl(u,i)," at ",a,"\n");
      fi;
    od;
  od;
  for i in [1..Length(uc.orderClasses)] do
    for a in uc.orderClasses[i] do
      if du[i]<=du[a] then 
        ChevieErr("bu(",i,":",uc.classes[i].name,")=",du[i]," <=bu(",
	  a,":",uc.classes[a].name,")=",du[a],"\n");
      fi;
    od;
  od;
  for s in uc.springerSeries do
    cl:=List(s.locsys,x->x[1]);
    order:=Incidence(Poset(uc)){cl}{cl};
    t:=ICCTable(uc,Position(uc.springerSeries,s));
    for i in [1..Length(cl)] do for j in [1..Length(cl)] do
      if t.scalar[i][j]<>0*t.scalar[i][j] and not order[i][j] then
        ChevieErr(Position(uc.springerSeries,s),"[",i,",",j,"]<>0 and not ",
	uc.classes[cl[i]].name,"<",uc.classes[cl[j]].name,"\n");
      fi;
    od;od;
  od;
  Print("\n");
  bc:=BalaCarterLabels(W);
  if ForAll(uc.classes,cl->IsBound(cl.dynkin)) then
  SortBy(bc,x->PositionProperty(uc.classes,cl->x[1]=cl.dynkin));
  bc:=List(bc,x->x[2]);
  for i in [1..Size(uc)] do
    cl:=uc.classes[i];
    if IsBound(cl.balacarter) then
      if cl.balacarter<>bc[i] then Error("balacarter");fi;
    else Print("bala-carter[",cl.name,"] unbound, should be",bc[i],"\n");
    fi;
  od;
  if ForAll(W.type,x->x.series in ["E","F","G"]) then
    for i in [1..Size(uc)] do
      name:=IsomorphismType(ReflectionSubgroup(W,List(bc[i],AbsInt)),
        rec(TeX:=true));
      if name="" then name:="1";fi;
      l:=Number(bc[i],x->x<0);
      if l<>0 then
        if '+' in name then name:=Replace(name,"+",SPrint("(a_",l,")+"));
        else Append(name,SPrint("(a_",l,")"));fi;
      fi;
      name:=String(Replace(name,"+","{+}"));
      cl:=uc.classes[i];
      if name<>cl.name then 
        Print("name=<",cl.name,"> but from Bala-Carter <",name,">\n");
      fi;
    od;
  fi;
  fi;
  for u in uc.classes do
    if IsBound(u.dimred) and IsBound(u.red) and
      u.dimred<>Dimension(u.red) then 
       ChevieErr("for ",u.name," dimred=",u.dimred,"<> Dim(",
          IsomorphismType(u.red,rec(torus:=true)),")=",Dimension(u.red),"\n");
    fi;
  od;
end,
#W->(IsCoxeterGroup(W) or IsCoxeterCoset(W)) and ForAll(Flat(CartanMat(W)),IsInt)];
W->IsCoxeterGroup(W) and ForAll(Flat(CartanMat(W)),IsInt));
