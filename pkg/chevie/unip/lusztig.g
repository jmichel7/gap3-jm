# find series h of levi HF in WF (sers may be .harish or .almostHarish)
# return a record with 2 fields: ser= found series
#   op: an element of W which conjugates h.levi to ser.levi
#       (never needed for Coxeter groups and standard levis H)
FindSeriesInParent:=function(h,HF,WF,sers)local n,p,W,y,split;
  split:=function(s)local res,p,l; #JM: I should replace by more accurate code
    res:=[];
    repeat p:=PositionSublist(s,"\\otimes ");l:=Length("\\otimes ");
      if p<>false then Add(res,s{[1..p-1]});s:=s{[p+l..Length(s)]};fi;
    until p=false;
    Add(res,s); res:=Filtered(res,x->x<>"");Sort(res); return res;
  end;
  if IsSpets(WF) then W:=Group(WF); else W:=WF;fi;
  n:=split(h.cuspidalName);
  for y in sers do
    if h.eigenvalue=y.eigenvalue and Length(h.levi)=Length(y.levi) then
      p:=split(y.cuspidalName);
      if p=n then 
         p:=RepresentativeOperation(W,Set(h.levi),Set(y.levi),OnSets);
        if p<>false then  return rec(ser:=y,op:=p);fi;
         p:=RepresentativeOperation(W,
          Set(Reflections(ReflectionSubgroup(W,h.levi))),
          Set(Reflections(ReflectionSubgroup(W,y.levi))),OnSets);
        if p<>false then  return rec(ser:=y,op:=p);fi;
      fi;
    fi;
  od;
  Error("series not found");
end;

# find the number of cuspidal of name n in UnipotentCharacters(HF)
FindCuspidalInLevi:=function(n,HF)local strip,cusp;
  strip:=function(n,s)local l;l:=Length(s);
    while Length(n)>l-1 and n{[1..l]}=s do n:=n{[l+1..Length(n)]};od;
    while Length(n)>l-1 and n{Length(n)-[l-1,l-2..0]}=s
    do n:=n{[1..Length(n)-l]};od;
    return n;
  end;
  cusp:=Position(List(UnipotentCharacters(HF).TeXCharNames,
     x->strip(x,"\\otimes ")),strip(n,"\\otimes "));
  if cusp=false then
    Error("cuspidal ",n," not found in ",HF,"\n");
  fi;
  return cusp;
end;

# l is a list of vectors each of length n. FindIntSol returns roots of unity
# x_i such that l[i]*[1,x2,..xn] is an integer for each i.
FindIntSol:=function(l)
  local apply,vals,try,smallnorm,stuck,simple,v,vars,simplify;
  vars:=Concatenation([1],List([2..Length(l[1])],i->SPrint("x",i)));
  simplify:=function()
    l:=List(l,function(p)
      if p=0 or Length(p.elm[1].elm)>0 or not IsRat(p.coeff[1]) then return p;fi;
      if p.coeff[1]<0 then p:=-p;fi;
      return p-Int(p.coeff[1]);end);
    l:=Set(Filtered(l,x->x<>0));
  end;

  l:=l*List(vars,Mvp);
  simplify();

# Print("FindIntSol",l,"\n");
# variable v is val: mvp or list-of-rootsofunity-possibilities.
  apply:=function(v,val)local i,f;
    i:=Position(vars,v);
    if IsBound(vals[i]) then # vals[i] is a list of possibilities
      if IsMvp(val) then return true;fi;
      vals[i]:=Intersection(vals[i],val);
      if Length(vals[i])=0 then return false;
      elif Length(vals[i])>1 then return true;fi;
      val:=vals[i];
    fi;
    vals[i]:=val;
    if IsList(vals[i]) then 
       if Length(vals[i])>1 then return true;
       else val:=vals[i][1];
       fi;
    fi;
    for i in Filtered([1..Length(vals)],j->IsBound(vals[j]) and IsMvp(vals[j]))
    do vals[i]:=Value(vals[i],[v,val]);
       f:=ScalMvp(vals[i]);if f<>false then vals[i]:=[f];fi;
    od;
    l:=List(l,x->Value(x,[v,val]));
    simplify();
#   Print(v,"->",val," vals=",vals,"\n");
    return true;
  end;
      
# if p has small norm must evaluate to 0
  smallnorm:=p->Sum(p.coeff,x->GetRoot(evalf(x*ComplexConjugate(x))))<15/16;

  simple:=function(p)local N,i;
    if Length(p.coeff)>2 or 
       (Length(p.coeff)=2 and Length(p.elm[1].elm)>0) then return false;fi;
    if Length(p.coeff)=1 then
      if Length(p.elm[1].elm)=0 then 
        InfoChevie("#I WARNING! FindIntSol cannot make:",p," integral\n");
        return 0;
      fi;
      v:=[0,p.coeff[1],p.elm[1].elm[1]];
    else v:=[p.coeff[1],p.coeff[2],p.elm[2].elm[1]];
    fi;
    N:=NofCyc(v{[1,2]});if N mod 2<>0 then N:=2*N;fi;
    i:=Filtered(List([0..N-1],i->E(N)^i),j->IsInt(v[1]+v[2]*j));
    if Length(i)=0 then
      InfoChevie("#I WARNING! FindIntSol cannot make:",p," integral\n");
      return 0;
    fi;
    return [v[3],i];
  end;

  try:=function()local v,i;
    if Length(l)=0 then return true;fi;
    for i in [1..Length(l)] do
      v:=simple(l[i]);
      if v<>false then
        if v=0 then return false;
        else l:=Drop(l,i);
          return apply(v[1],v[2]);
        fi;
      fi;
    od;
    i:=PositionProperty(l,smallnorm);
    if i<>false then v:=l[i];
      i:=PositionProperty(v.elm,x->Length(x.elm)>0);
      return apply(v.elm[i].elm[1],Mvp(v.elm[i].elm[1])-v/v.coeff[i]);
    fi;
    stuck:=true;return false;
  end;

  vals:=[[1]];stuck:=false;
  while try() do if Length(l)=0 then return Cartesian(vals);fi;od;
  if stuck then InfoChevie("#I WARNING! FindIntSol: stuck ",l,"\n");fi;
  return false;
end;

LusztigInductionPieces:=function(res)
  local LF,WF,W,L,h,uL,uW,hw,p,WFGL,LFGL,w,WGL,res,Jb,rh,L,ser,piece,i,j;
  WF:=res.g;LF:=res.u;
  L:=Group(LF); W:=Group(WF);
  if IsCoxeterGroup(L) and not 
    IsSubset(W.rootInclusion{W.generatingReflections},
           L.rootInclusion{L.generatingReflections}) then
    w:=StandardParabolic(W,L);L:=L^w;LF:=Spets(L,LF.phi^w);
  fi;
  uL:=UnipotentCharacters(LF); uW:=UnipotentCharacters(WF);
  if uL=false or uW=false then return false;fi;
  res.gNames:=function(t,option)return CharNames(uW,option);end;
  res.uNames:=function(t,option)return CharNames(uL,option);end;
  res.pieces:=[];
  hw:=uW.almostHarishChandra;
  for h in uL.almostHarishChandra do
    p:=FindSeriesInParent(h,LF,WF,hw);ser:=p.ser;
    L:=ReflectionSubgroup(W,ser.levi); # L^op contained in LF
    if IsBound(WF.isCoxeterCoset) then
      WFGL:=RelativeCoset(WF,h.levi);
      LFGL:=RelativeCoset(LF,h.levi);
      LFGL:=CoxeterSubCoset(WFGL,
        Group(WFGL).rootInclusion{List(Group(LFGL).relativeIndices,
          x->Position(Group(WFGL).relativeIndices,x))},
            Group(WFGL).MappingFromNormalizer(LF.phi*WF.phi^-1));
    else
      Jb:=Concatenation(List(ser.relativeType,
         function(x)if IsBound(x.orbit) then return
          Concatenation(List(x.orbit,y->y.indices));
          else return x.indices;fi;end));
      WFGL:=RelativeCoset(WF,ser.levi,Jb);
      if WFGL=false then return false;fi;
      WGL:=Group(WFGL);
      rh:=OnTuples(Concatenation(List(h.relativeType,
       function(x)if IsBound(x.orbit) then return
        Concatenation(List(x.orbit,y->y.indices));
        else return x.indices;fi;end)),p.op);
      w:=WGL.MappingFromNormalizer((LF.phi^p.op)*WF.phi^-1);
      if w=false then Error("Could not compute MappingFromNormalizer\n");fi;
      rh:=List(rh,function(x)local r;
         r:=GetRelativeRoot(W,L,x);
         return First(WGL.rootInclusion,
           i->Reflection(WGL,WGL.rootRestriction[i])=
            PermMatX(WGL,Reflection(r.root,r.coroot)));end);
      LFGL:=SubSpets(WFGL,rh,w);
      ReflectionName(LFGL); ReflectionName(WFGL); # to shut up CharTable
    fi;
    piece:=rec(wnum:=ser.charNumbers,
               hnum:=h.charNumbers,
               uNames:=function(t,option)
                 return AlmostCharNames(uL,option){t.hnum};end,
               gNames:=function(t,option)
                 return AlmostCharNames(uW,option){t.wnum};end,
                scalar:=ComplexConjugate(InductionTable(LFGL,WFGL).scalar),
#              scalar:=InductionTable(LFGL,WFGL).scalar,
               u:=LFGL,
               g:=WFGL,
               what:="piece",
               operations:=InductionTableOps);
    if h.levi=[] then piece.head:=function(t,option)
      if IsBound(option.TeX) then
       return SPrint("from $",ReflectionName(LF,option),
                     "$ to $",ReflectionName(WF,option),"$");
      else return SPrint("from ",ReflectionName(LF,option),
                     " to ",ReflectionName(WF,option));
      fi;end;
    else piece.head:=function(t,option)local sub,a,b;
      if IsBound(option.TeX) then sub:=L->SPrint("_{",L,"}");
                             else sub:=L->SPrint("_[",L,"]");fi;
      a:=SPrint("W",sub(ReflectionName(LF,option)),"(",sub(Join(h.levi)),",",
       TeXStrip(h.cuspidalName,option),")=",ReflectionName(t.u,option));
      b:=SPrint("W",sub(ReflectionName(WF,option)),"(",sub(Join(ser.levi)),",",
       TeXStrip(ser.cuspidalName,option),")=",ReflectionName(t.g,option));
      if IsBound(option.TeX) then return SPrint("from $",a,"$\nto $",b,"$");
      else return SPrint("from ",a,"\nto   ",b);
      fi;end;
    fi;
    Add(res.pieces,piece);
  od;
  return res;
end;

CHEVIE.Cache.LusztigInductionMaps:=false;
# arguments (LF,WF[,force])
# if there is a third argument then return res (with the pieces)
# even in case of failure
LusztigInductionTable:=function(arg)
  local LF,WF,res,uW,uL,hh,scalars,fL,fWinv,map,maps,ret,p;
  LF:=arg[1];WF:=arg[2];
  if not IsSpets(WF) then WF:=Spets(WF);fi;
  if not IsSpets(LF) then LF:=Spets(LF);fi;
# if not IsBound(WF.isCoxeterCoset) and WF.phi<>() then 
#   ChevieErr("Non-trivial parent Spets not implemented");
#   return false;
# fi;
  res:=rec(u:=LF, g:=WF, what:="LusztigInduction",
    head:=function(t,option)local math;
      if IsBound(option.TeX) then math:=x->SPrint("$",x,"$");else math:=x->x;fi;
      return SPrint("Lusztig Induction from ",math(ReflectionName(t.u,option)),
        " to ",math(ReflectionName(t.g,option)));end, 
    operations:=InductionTableOps);
  uW:=UnipotentCharacters(WF);
  res:=CHEVIE.GetCached(uW,"LusztigInductionMaps",res,x->
    [Group(x.u).rootInclusion{Group(x.u).generatingReflections},
      x.u.phi*x.g.phi^-1]);
  if IsBound(res.scalar) then return res;fi;
  res:=LusztigInductionPieces(res);if res=false then return false;fi;
  uL:=UnipotentCharacters(LF);
  fL:=UnipotentCharactersOps.Fourier(uL);
  hh:=uL.almostHarishChandra;
  fWinv:=UnipotentCharactersOps.FourierInverse(uW);
  maps:=List([1..Length(hh)],function(i)local piece;
    map:=List([1..Size(uW)],y->[1..Size(uL)]*0);
    piece:=res.pieces[i];
    map{piece.wnum}{piece.hnum}:=piece.scalar;
    return fWinv*map*fL;
  end);
  ret:=function(map)
    if map=false then 
      ChevieErr("Failed\n");
      if Length(arg)=2 then return false; else return res; fi;
    fi;
    res.scalar:=map;
    return res;
  end;
  map:=Sum(maps);
  if ForAll(map,v->ForAll(v,IsInt)) then 
    if LF.phi=WF.phi and not ForAll(Flat(map),x->x>=0) then 
      ChevieErr("non-positive RLG for untwisted L");
    fi;
    return ret(map);
  fi;
  scalars:=FindIntSol(Set(TransposedMat(List(maps,Concatenation))));
  if scalars=false then return ret(scalars);fi;
  if Length(scalars)>1 then
    ChevieErr("#I WARNING: ambiguity in scalars:",FormatGAP(scalars),"\n");
  fi;
  scalars:=scalars[1];
  if ForAny(scalars,IsMvp) then Error();fi;
  if ForAny(scalars,x->x<>1) then
    res.scalars:=scalars;
    if ForAll(scalars,IsInt) then p:="#I signs are ";
    else InfoChevie("#I non-sign scalars needed:",FormatGAP(scalars),"\n");
    fi;
#   ChevieErr(p,Join(List([1..Length(scalars)],function(i)local n;
#     n:=hh[i].cuspidalName;if n="" then n:="Id";fi;
#     return SPrint(n,"->",scalars[i]);end),";"),"\n");
  fi;
  scalars:=Sum([1..Length(scalars)],i->maps[i]*scalars[i]);
  if LF.phi=WF.phi and not ForAll(Flat(scalars),x->x>=0) then 
    ChevieErr("non-positive RLG for untwisted L");
  fi;
  return ret(scalars);
end;

CHEVIE.Cache.HCInductionMaps:=true;
HarishChandraInductionTable:=function(HF,WF)
  local W,H,h,uh,uw,name,p,Wi,Hi,rh,Jb,t,N,piece,L,res,st,getHi;
# st:=Runtime();
  if not IsSpets(WF) then WF:=Spets(WF); fi;
  uw:=UnipotentCharacters(WF); W:=Group(WF);
  if not IsSpets(HF) then HF:=Spets(HF); fi;
  uh:=UnipotentCharacters(HF); H:=Group(HF);
  res:=rec(u:=HF, g:=WF, what:="HCInductionTable",
    uNames:=function(t,option)return CharNames(uh,option);end,
    gNames:=function(t,option)return CharNames(uw,option);end,
    head:=function(t,option)local math;
      if IsBound(option.TeX) then math:=x->SPrint("$",x,"$");else math:=x->x;fi;
      return SPrint("Harish-Chandra Induction from ",
      math(ReflectionName(t.u,option))," to ",math(ReflectionName(t.g,option)));
    end, 
    operations:=InductionTableOps
  );
# Print("HF=",HF," WF=",WF,"\n");
  res:=CHEVIE.GetCached(uw,"HCInductionMaps",res,
    x->Group(x.u).rootInclusion{Group(x.u).generatingReflections});
  if IsBound(res.scalar) then return res;fi;
  res.pieces:=[];res.scalar:=List([1..Size(uw)],y->[1..Size(uh)]*0);
  for h in uh.harishChandra do
    p:=FindSeriesInParent(h,HF,WF,uw.harishChandra);
    Jb:=Concatenation(List(p.ser.relativeType,x->x.indices));
    if IsCoxeterCoset(WF) then
      Wi:=ApplyFunc(ReflectionGroup,p.ser.relativeType);
#        InfoChevie("W_",ReflectionName(WF),"(",IntListToString(p.ser.levi),")=",
#           ReflectionName(Wi),"\n");
      if p.op<>() then
        rh:=RelativeGroup(H,h.levi).generators;
        rh:=Filtered(Arrangements([1..Length(Wi.generators)],Length(rh)),
          a->ForAll([1..Length(a)],i->OrderPerm(Wi.generators[a[i]])=
                                                OrderPerm(rh[i])));
        if Length(rh)>1 then
          ChevieErr("WARNING: embedding ambiguous:",rh,"\n");
        fi;
        Hi:=ReflectionSubgroup(Wi,rh[1]);
      else Hi:=ReflectionSubgroup(Wi,
             List(Concatenation(List(h.relativeType,x->x.indices)),
               x->PositionProperty(Jb,y->x in Cycle(WF.phi,y))));
      fi;
    else
      L:=ReflectionSubgroup(W,p.ser.levi);
      Wi:=RelativeGroup(W,p.ser.levi,Jb);
#     if p.op=() then 
      if false then
        if h.levi=[] then rh:=H.rootInclusion{H.generatingReflections};
        else rh:=Filtered(Jb,i->Reflection(Parent(W),i) in H);
        fi;
      else 
        rh:=ReflectionSubgroup(W,
              OnTuples(H.rootInclusion{H.generatingReflections},p.op));
        if not IsSubset(rh.rootInclusion,L.rootInclusion) then Error();fi;
        rh:=Filtered(rh.rootInclusion,i->not Reflection(Parent(W),i) in L);
      fi;
      getHi:=function()local rr,r,sH,x,p,Hi;
        sH:=Size(RelativeGroup(H,h.levi));
        rr:=[];
        for x in rh do
          r:=GetRelativeRoot(W,L,x);
          p:=PositionProperty(Wi.rootInclusion,
           i->Reflection(Wi,Wi.rootRestriction[i])=
            PermMatX(Wi,Reflection(r.root,r.coroot)));
          if p<>false then AddSet(rr,p);fi;
          Hi:=ReflectionSubgroup(Wi,rr);
          if Size(Hi)=sH then return Hi;fi;
        od;
        return ReflectionSubgroup(Wi,[]);
      end;
      Hi:=getHi();
    fi;
#    ReflectionName(Wi);ReflectionName(Hi); # to shut up CharTable
    piece:=rec(
     uNames:=function(t,option)return CharNames(uh,option){h.charNumbers};end,
     gNames:=function(t,option)return CharNames(uw,option){p.ser.charNumbers};end,
     scalar:=InductionTable(Hi,Wi).scalar,
     u:=Hi,
     g:=Wi,
     what:="piece",
     operations:=InductionTableOps);
    Add(res.pieces,piece);
    res.scalar{p.ser.charNumbers}{h.charNumbers}:=piece.scalar;
  od;
# Print("HCInd^",WF,"_",HF,":",Elapsed(st),"\n");
  return res;
end;
