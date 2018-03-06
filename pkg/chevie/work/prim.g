ParamChars:=function(p,q,r)local t,i,tt,res;
  res:=[];
  for t in PartitionTuples(r,p) do
    tt:=List([1..q]*p/q,i->t{Concatenation([i+1..Length(t)],[1..i])});
    if t=Minimum(tt) then
      r:=Position(tt,t);
      t:=t{[1..r*p/q]};
      if r=q then Add(res,t);
      else r:=q/r;for i in [0..r-1] do Add(res,Concatenation(t,[E(r)^i]));od;
      fi;
    fi;
  od;
  return res;
end;

embed:=function(w,p,q)
  if p>q then return Replace(w-1,[1],Concatenation([1..p-1]*0+1,[2,1]),
                             [0],[1..q]*0+1);
  else return Replace(w,[1],Concatenation([1..p-1]*0+1,[2,1]));
  fi;
end;

CharMap:=function(p,q,r)local W,W1,ct1,p,cl,v,ct,c,cc,k,z,pos,w,vw,good,cv;
  W:=ComplexReflectionGroup(p,q,r);
  cc:=ParamChars(p,q,r);
  W1:=ComplexReflectionGroup(p,1,r);
  ct1:=CharTable(W1).irreducibles;
  ct:=CharTable(W).irreducibles;
  cl:=List(ChevieClassInfo(W).classtext,x->PositionClass(W1,
    EltWord(W1,embed(x,p,q))));
  good:=Filtered([1..Length(cl)],i->Number(cl,j->j=cl[i])=1);
  for c in cc do
    if IsCyc(c[Length(c)]) then
      k:=p/(Length(c)-1);z:=c[Length(c)];
      pos:=Position(CharParams(W1),
        [Concatenation(List([1..k],i->c{[1..Length(c)-1]}))]);
      cv:=List([0..k-1],l->c[Length(c)]^l);
      vw:=List(ChevieClassInfo(W).classtext,w->
        List([0..k-1],i->PositionClass(W1,W1.1^i*EltWord(W1,embed(w,p,q)))));
      vw:=TransposedMat(vw);
      v:=List(vw,x->ct1[pos]{x});
      v:=cv*v/k;
      pos:=Filtered([1..Length(ct)],i->ct[i]{good}=v{good});
      Print(CHEVIE.imp.CharName(c),"=>",CharNames(W){pos},"\n");
    else 
      pos:=Position(CharParams(W1),[c]);
      v:=ct1[pos]{cl}; pos:=Position(ct,v);
      Print(CHEVIE.imp.CharName(c),"=>",CharNames(W)[pos],"\n");
    fi;
  od;
end;

reps:=function(p,q,r,H)local W,t,para;
  W:=Group(H);
  return List(ParamChars(p,q,r),S->CHEVIE.GpqrRep(p,q,r,H.parameter,S));
end;

c2:=function(p,q,r,para)local H,i,W,cl,ct,c,n,pos,p,rep,par,res;
  W:=ComplexReflectionGroup(p,q,r);H:=Hecke(W,para);
  cl:=ChevieClassInfo(W).classtext;
  n:=CharNames(W);
  if ForAll(H.parameter,x->x=List([1..Length(x)],i->E(Length(x))^(i-1))) then
    ct:=CharTable(W).irreducibles;
  else ct:=CharTable(H).irreducibles;
  fi;
  par:=List(ParamChars(p,q,r),CHEVIE.imp.CharName);
  rep:=reps(p,q,r,H);
  res:=[];
  for i in [1..Length(rep)] do
     c:=CharRepresentationWords(rep[i],cl)*Product(H.parameter,Product)^0;
     pos:=Position(ct,c);
     if pos=false then c:=SPrint(c,":not found");else c:=n[pos];fi;
     Print(i,":",par[i]," => ",pos,":",c,"\n");
     Add(res,pos);
  od;
  return res;
end;
