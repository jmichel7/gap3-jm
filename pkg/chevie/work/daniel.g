# Compare:
# order on chars coming from Hecke/Cherednik algebras (formule)
# order coming from generalized Springer correspondence
#   (partial order on unipotents, or vanishing of IC coefficients)
#
formule:=function(W)local p,cl;
  p:=List(HyperplaneOrbits(W),x->x.classno[1]);cl:=ChevieClassInfo(W).classes;
  return List(CharTable(W).irreducibles,ch->List(p,j->(1-ch[j]/ch[1])*cl[j]));
end;

Printformule:=function(W)local l,d,m;
  l:=formule(W);
  d:=1+List([1,2],i->Maximum(List(l,p->p[i])));
  m:=List([1..d[1]],i->List([1..d[2]],j->"."));
  for i in [1..Length(l)] do m[1+l[i][1]][1+l[i][2]]:=i;od;
  Print(Format(Reversed(m)),"\n");
end;

# order(W,coeffs on Hplanes[, use indices instead of charnames])
order:=function(arg)local W,c,n;
  W:=arg[1];c:=arg[2];
  if Length(arg)=2 then n:=CharNames(W);
  else n:=[1..NrConjugacyClasses(W)];
  fi;
  Print(Join(List(CollectBy(n,formule(W)*c),Join),"<"),"\n");
end;

# spr(W,seriesno[,opt])
# opt.ic:  use ICCtable instead of po on unip classes
# opt.no:  use indices instead of charnames
spr:=function(arg)local W,i,uc,s,p,opt;
  W:=arg[1];i:=arg[2];uc:=UnipotentClasses(W);
  if Length(arg)=2 then opt:=rec();else opt:=arg[3];fi;
  s:=uc.springerSeries[i];
  if IsBound(opt.ic) then
    p:=Poset(List(ICCTable(uc,i).scalar,x->List(x,y->y<>0*y)));
  else p:=Restricted(Poset(uc),List(s.locsys,x->x[1]));
  fi;
  p.group:=s.relgroup;
  if IsBound(opt.no) then
    p.label:=function(p,n,opt)return n;end;
  else
    p.label:=function(p,n,opt)
      return CharName(p.group,CharParams(p.group)[n],opt);end;
  fi;
  Display(p);
end;

# Printformule(CoxeterGroup("C",3));
# spr(RootDatum("spin",18),4,rec(no:=true));
# spr(RootDatum("spin",18),4,rec(no:=true,ic:=true));
# order(CoxeterGroup("C",3),-[2,1],true);
# spr(RootDatum("spin",12),3,rec(no:=true));
# spr(RootDatum("spin",12),3,rec(no:=true,ic:=true));
# order(CoxeterGroup("C",3),-[1,2],true);
# spr(RootDatum("spin",22),3,rec(no:=true));
# spr(RootDatum("spin",22),3,rec(no:=true,ic:=true));
# order(CoxeterGroup("C",3),-[4,1],true);
# spr(CoxeterGroup("C",3),1,rec(no:=true));
# spr(CoxeterGroup("C",3),1,rec(no:=true,ic:=true));
# order(CoxeterGroup("C",3),-[1,1],true);
# spr(RootDatum("spin",13),3,rec(no:=true));
# spr(RootDatum("spin",13),3,rec(no:=true,ic:=true));
# spr(RootDatum("spin",15),2,rec(no:=true));
# spr(RootDatum("spin",15),2,rec(no:=true,ic:=true));
# spr(RootDatum("spin",15),3,rec(no:=true));
# spr(RootDatum("spin",15),3,rec(no:=true,ic:=true));
