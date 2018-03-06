# RightTauCells(W[,l]) cuts in right-tau-cells the list l
RightTauCells:=function(arg)local W,d,sts,l,i,st,p,changed,prev,s,t,cand,ILD;
  W:=arg[1];
  ILD:=W.operations.IsLeftDescending;
  sts:=Filtered(BraidRelations(W),r->Length(r[1])>2);
  if Length(arg)=1 then l:=Elements(W);else l:=arg[2];fi;
  l:=List(CollectBy(Elements(W),x->RightDescentSet(W,x)),Set);
  changed:=[1..Length(l)];
  Print("#I \c");
  repeat
    prev:=changed;changed:=[]; i:=1;
    while i<=Length(l) do
      d:=LeftDescentSet(W,l[i][1]^-1);
      for st in sts do s:=Intersection(List(st,x->x[1]),d);
        if Length(s)=1 then 
	  t:=Difference(List(st,x->x[1]),s)[1];s:=s[1];
	  st:=First(st,x->x[1]=s);
	  cand:=l{Filtered(prev,function(i)local w;w:=l[i][1]^-1;
	    return ILD(W,w,t) and not ILD(W,w,s);end)};
	  p:=CollectBy(l[i],function(w)local ss; ss:=LeftStarNC(W,st,w^-1)^-1;
	    return PositionProperty(cand,y->ss in y); end);
	  if Length(p)>1 then
	    p:=List(p,Set);
	    Add(changed,i);Append(changed,Length(l)+[1..Length(p)-1]);
	    l[i]:=p[1];Append(l,p{[2..Length(p)]});
	  fi;
	fi;
      od;
      i:=i+1;
    od;
    Print(Length(l)," \c");
  until changed=[];
  Print("\n");
  return l;
end;

RightStarsOnSet:=W->List(Filtered(BraidRelations(W),r->Length(r[1])>2),
  st->function(l)local s,rst;
  s:=Intersection(RightDescentSet(W,l[1]),List(st,x->x[1]));
  rst:=First(st,x->x[1]=s[1]);
  return Set(List(l,w->LeftStarNC(W,rst,w^-1)^-1));
end);

condense:=function(W,l)
  l:=List(FOrbits(RightStarsOnSet(W),l),x->x[1]);
  return List(l,y->List(FOrbits(LeftStars(W),y),z->z[1]));
end;

expand:=function(W,l)
  l:=List(l,x->Union(List(x,w->FOrbit(LeftStars(W),w))));
  return Union(List(l,x->FOrbit(RightStarsOnSet(W),x)));
end;

RightStarOnList:=function(W,st)return function(l)local s,rst;
  s:=Intersection(RightDescentSet(W,l[1]),List(st,x->x[1]));
  if Length(s)<>1 then return l;fi;
  rst:=First(st,x->x[1]=s[1]);
  return List(l,w->LeftStarNC(W,rst,w^-1)^-1);
end;end;

CheckKL:=function(W,l)local sts,i,st,e,p,n;
  sts:=Filtered(BraidRelations(W),x->Length(x[1])>2);
  for i in [1..Length(l)] do
    Print(i,"->");
    for st in sts do 
      e:=RightStarOnList(W,st)(Elements(l[i]));
      p:=PositionProperty(l,x->Set(e)=Elements(x));
      if List(Elements(l[i]),x->LeftDescentSet(W,x))<>
         List(e,x->LeftDescentSet(W,x)) then Error("LD");fi;
      Print(p);
      n:=[1..Length(e)];SortParallel(e,n);
      if l[p].mu=l[i].mu{n}{n} then Print("\c* ");
      else Error();
      fi;
    od;
    Print("\n");
  od;
end;

rep2graph:=function(r,v)local dim,nodes,mu,i,j,k,l,n;
  dim:=Length(r[1]);
  nodes:=List([1..dim],i->Filtered([1..Length(r)],j->r[j][i][i]=-1));
  n:=[1..dim];
  SortParallel(nodes,n);
  r:=List(r,x->x{n}{n});
  mu:=List([1..dim],i->[1..dim]*0+99);
  for j in [1..Length(r)] do
    for i in [1..dim] do
      for k in Difference([1..dim],[i]) do
        if r[j][i][k]<>0 then
	  if not j in nodes[k] then Error("a");fi;
	  if mu[k][i]<>99 and mu[k][i]<>r[j][i][k]/v then Error("b");fi;
	  mu[k][i]:=r[j][i][k]/v;
	fi;
  od;od;od;
  l:=Concatenation(List([1..dim],i->List([i+1..dim],j->[mu[i][j],mu[j][i],i,j])));
  l:=Filtered(l,x->x[1]<>99 or x[2]<>99);
  l:=List(l,function(j)
    if j[1]=99 or j[1]=j[2] then return j{[2,3,4]};
    elif j[2]=99 then return j{[1,3,4]};
    else return [j{[1,2]},j[3],j[4]];
    fi;end);
  l:=CollectBy(l,x->x[1]);
  l:=List(l,function(v)
    v:=CollectBy(v,x->x[2]);
    return [v[1][1][1],List(v,x->Concatenation([x[1][2]],List(x,y->y[3])))];end);
  nodes:=Concatenation(List(Collected(nodes),
    function(p)if p[2]=1 then return [p[1]];else return [p[1],p[2]-1];fi;
    end));
  return [nodes,l];
end;

CheckCellReps:=function(W)local H,tbl,cl,c,ch,char,i,r,l;
  H:=Hecke(W);tbl:=CharTable(W);cl:=ChevieClassInfo(W).classtext;
  l:=LeftCells(W);
  for c in l do
    r:=Representation(c,H);ch:=CharRepresentationWords(r,cl);
    ch:=List(tbl.irreducibles,x->ScalarProduct(tbl,x,ch));
    char:=[];
    for i in [1..Length(ch)] do Append(char,[1..ch[i]]*0+i);od;
    if PermListList(char,c.character)=false then 
     Error("in cell no. ",Position(l,c)," should be ",char," not ",c.character);fi;
    Print(".\c");
  od;
end;
