#############################################################################
##
#A checkcharpara.g   CHEVIE library      Jean Michel and Ulrich Thiel
##
#Y  Copyright (C) 2017  University  Paris VII and Universitat Stuttgart
##
## This file is destined to work even with old versions of CHEVIE/GAP3
## It checks that the parameters/names for characters of Complex reflection
## groups agree with the description (since 2016) in the CHEVIE manual.
## If not, it tells what permutation of the data is needed.

if not IsBound(ChevieErr) then
  ChevieErr:=function(arg)
    Print("****** WARNING!: CheckCharParams:");
    ApplyFunc(Print,arg);
  end;
fi;

# check the charparams are consistent with actual data in various ways
CheckCharParams:=function(W) local ct,q,fd,db,n,l,nm,t,f,fd,ok,ddb,occursin,
  issign,fakdeghighterms,fakdegsmaller,isconjugate,inducedfromid,value,property,
  rules,mintensor,lexpol,perm,DecomposeTensor,Galois,isgal;
  ct:=CharTable(W).irreducibles;
  q:=X(Rationals);q.name:="q";fd:=FakeDegrees(W,q);
  db:=List(fd,x->[Value(x,1),x.valuation]);
  n:=ReflectionName(W);
  if IsBound(ChevieCharInfo(W).malle) then
    l:=ChevieCharInfo(W).malle;
    nm:=List(l,x->CHEVIE.R("CharName","G2")(x,rec()));
  elif n[1] in "EFGH" and Length(n)<4 then
    l:=Filtered(Collected(CharParams(W)),x->x[2]>1);
    if Length(l)>0 then ChevieErr("Duplicated:",l,"\n");fi;
    l:=List(CharParams(W),x->x[1]);
    nm:=CharNames(W);
  elif n[1] in ".ABCD" or Length(n)>=4 then return;
  else
    l:=List(CharParams(W),x->x[1]);
    nm:=CharNames(W);
  fi;
  ddb:=Filtered(List(Set(db),x->Filtered([1..Length(fd)],i->db[i]=x)),
      x->Length(x)>1);
  f:=x->Position(nm,x);
  # lexicographic polynomial order starting with highest coeff
  lexpol:=function(p,q)local a,b;
    a:=Concatenation([1..p.valuation]*0,p.coefficients);
    b:=Concatenation([1..q.valuation]*0,q.coefficients);
    if Length(a)<Length(b) then Append(a,[1..Length(b)-Length(a)]*0);fi;
    if Length(b)<Length(a) then Append(b,[1..Length(a)-Length(b)]*0);fi;
    return Reversed(a)<=Reversed(b);
  end;
  # property(n,f)
  # chars of name n+qualifiers satisfy f
  property:=function(n,fn)local n,qq,a,i,j;
    n:=PositionsProperty(nm,x->Length(x)>=Length(n) and x{[1..Length(n)]}=n);
    qq:=List(nm{n},i->Number(i,u->u='\'')+2*Number(i,u->u='"'));
    SortParallel(qq,n);
    a:=First(Arrangements([1..Length(n)],Length(n)),x->fn(nm{n{x}}));
    perm:=perm*MappingPermListList(n,n{a});
    for i in n do j:=PositionProperty(ddb,x->i in x);
      if j=false then if i=n[1] then Print("useless test: ",nm{n},"\n");fi;
      else ddb[j]:=Difference(ddb[j],[i]);
       if Length(ddb[j])<=1 then ddb:=Drop(ddb,j);fi;
      fi;
    od;
  end;
  ok:=function(arg)
    if not arg[1] then ApplyFunc(ChevieErr,arg{[2..Length(arg)]});fi;
    return arg[1];
  end;
  DecomposeTensor:=function(arg)local t;
    t:=List(TransposedMat(ct{arg}),Product);
    return MatScalarProducts(CharTable(W),ct,[t])[1];
  end;
  Galois:=function(x,e)
    if IsList(x) then return List(x,y->Galois(y,e));else return x^e;fi;
  end;
  occursin:=function(arg)# arg[1] occurs in arg[2] tensor .. tensor arg[n]
   return ok(0<>ApplyFunc(DecomposeTensor,List(arg{[2..Length(arg)]},f))
    [f(arg[1])],arg[1]," should occur in tensor",arg{[2..Length(arg)]});
  end;
  fakdeghighterms:=function(a,p) return ok(Degree(fd[f(a)]-p)<p.valuation,
    "FakeDegree(",a,") should have high terms ",p);
  end;
  fakdegsmaller:=function(a,b) return ok(Degree(fd[f(a)])<Degree(fd[f(b)]),
    "FakeDeg(",a,") should have a smaller degree than FakeDeg(",b,")");
  end;
  issign:=function(a,b)local det;
    det:=List(W.generators,x->PositionClass(W,x));
    det:=PositionProperty(ct,l->ForAll(l{det},y->y=-1));
    return ok(1=DecomposeTensor(det,f(b))[f(a)],
       a," should be tensored by sign of ",b);
  end;
  isconjugate:=function(a,b) return ok(ct[f(a)]=ComplexConjugate(ct[f(b)]),
       a," should be complex conjugate to ",b);
  end;
  isgal:=function(a,b,g) return ok(ct[f(a)]=Galois(ct[f(b)],g),
       a," should be conjugate by ",g," to ",b);
  end;
  inducedfromid:=function(a,J)local L,t;
    L:=ReflectionSubgroup(W,J);
    t:=InductionTable(L,W);
    return ok(t.scalar[f(a)][PositionId(L)]<>0,
       a," should occur in the induced of Id from ",ReflectionName(L));
  end;
  value:=function(a,c,v)c:=PositionClass(W,W.(c));
    return ok(CharTable(W).irreducibles[f(a)][c]=v,
      a," should take value ",v," at class ",c);
  end;
  # lexicographically minimum tensor of (d,b)-distinguishable chars
  # in which char no. i occurs
  mintensor:=function(W,i)local bydb,b,a,ch,gb,m;
    bydb:=Filtered([2..Length(ct)],i->not '\'' in nm[i]);
    b:=ChevieCharInfo(W).b;
    SortParallel(List(bydb,i->[ct[i][1],b[i]]),bydb);
    for a in bydb do
      ch:=List([1..Length(ct[i])],j->ComplexConjugate(ct[a][j])*ct[i][j]);
      gb:=Filtered(bydb,j->ct[j][1]<=ct[i][1]*ct[a][1]);
      m:=MatScalarProducts(CharTable(W),ct{gb},[ch])[1];
      b:=PositionProperty(m,x->x<>0);
      if b<>false then return [a,gb[b],m[b]];fi;
    od;
  end;
  rules:=function()local b,p,q,bad,pp,g,u,bb,o,rule;
    q:=i->Number(nm[i],u->u='\'')+2*Number(nm[i],u->u='"');
    g:=Field(Flat(CartanMat(W)));
    if g<>Rationals then
      g:=List(GaloisGroup(g).generators,x->PermListList(ct,Galois(ct,x)));
      g:=Orbits(Group(g,()),[1..Length(nm)]);
      o:=i->First(g,o->i in o);
    fi;
    rule:=[];
    bad:=List(ddb,function(p)local fin;
      fin:=function(arg)
        if not arg[2] then Add(rule,-arg[1]);
          ApplyFunc(ChevieErr,Concatenation(["Rule ",arg[1]," failed:"],
            arg{[3..Length(arg)]}));
          return (p[1],p[2]);
        else Add(rule,arg[1]);return ();
        fi;
      end;
      if Length(p)>2 then return true;fi;
      if fd[p[1]]<>fd[p[2]] then
	return fin(1,lexpol(fd[p[1]],fd[p[2]])=(q(p[1])<q(p[2])),
	  nm[p[1]],"->",fd[p[1]]," and ",nm[p[2]],"->",fd[p[2]],"\n");
      fi;
      if IsBound(o) then
	pp:=List(p,i->Filtered(o(i),j->not j in Flat(ddb)));
	if Number(pp,x->x<>[])=1 then
	  u:=PositionProperty(pp,x->x<>[]);
	  return fin(2,q(p[u])=1,nm[p[u]],"~",Join(nm{pp[u]}),"\n");
	fi;
	pp:=List(p,o);
	bb:=List(pp,i->List(i,y->db[y][2]));for u in bb do Sort(u);od;
	if bb[1]<>bb[2] then
	  return fin(3,(bb[1]<bb[2])=(q(p[1])<q(p[2])),
	   nm[p[1]],"~",Join(nm{pp[1]})," and ",nm[p[2]],"~",Join(nm{pp[2]}),"\n");
	fi;
      fi;
      bb:=List(p,i->mintensor(W,i));
      pp:=List(bb,a->[db[a[1]][1],db[a[1]][2],db[a[2]][1],db[a[2]][2],a[3]]);
      if pp[1]<>pp[2] then
	return fin(4,(pp[1]<pp[2])=(q(p[1])<q(p[2])),
	      nm[p[1]],"|",Join(nm{bb[1]{[1,2]}},"o "),"=",pp[1][5]," & ", 
	      nm[p[2]],"|",Join(nm{bb[2]{[1,2]}},"o "),"=",pp[2][5],"\n");
      fi;
      return true;
    end);
    if Length(rule)<>0 then InfoChevie("#I applied  rules ",
      Join(List(Collected(rule),x->Join(x,"x"))),"\n");fi;
    perm:=Product(Filtered(bad,x->x<>true));
    if perm=1 then perm:=();fi;
    ddb:=ListBlist(ddb,List(bad,x->x=true));
  end;
  if IsBound(l) and db<>List(l,x->x{[1,2]}) then
    ChevieErr("disagree with (\\chi(1),b_\\chi)\n");
    PrintArray(Filtered(TransposedMat([db,List(l,x->x{[1,2]})]),
      x->x[1]<>x[2]));
    Error();
  fi;
  perm:=();
  if n in ["G334","G3,3,4"] then rules();
    property("phi{6,5}",v->occursin(v[1],"phi{4,1}","phi{4,1}"));
  elif n in ["G335","G3,3,5"] then rules();
    property("phi{30,7}",
    v->occursin(v[2],"phi{5,1}","phi{5,1}","phi{5,1}","phi{5,1}"));
  elif n in ["G443","G4,4,3"] then rules();
    property("phi{3,2}",v->isconjugate(v[2],"phi{3,1}"));
    property("phi{3,6}",v->isconjugate(v[2],"phi{3,5}"));
  elif n="G2" then rules();
    property("phi{1,3}",v->value(v[1],1,-1)); # W.rootLengths=[3,1]
  elif n="G5" then rules();
    property("phi{1,4}'",v->value(v[1],1,1));
    property("phi{1,8}",v->isconjugate(v[2],"phi{1,16}")
      and isconjugate(v[1],"phi{1,4}'")
      and isconjugate(v[3],"phi{1,4}''"));
    property("phi{1,12}",v->value(v[1],1,E(3)) and isconjugate(v[1],v[2]));
    property("phi{2,5}",v->isconjugate(v[2],"phi{2,1}") and value(v[1],1,-1));
    property("phi{2,7}",v->isconjugate(v[1],"phi{2,5}'") 
      and isconjugate(v[2],"phi{2,5}'''"));
    property("phi{2,3}",v->value(v[1],1,-E(3)) and isconjugate(v[2],v[1]));
  elif n="G7" then rules();
    property("phi{1,4}",v->value(v[1],2,1));
    property("phi{1,8}",v->isconjugate(v[2],"phi{1,16}") and
                           isconjugate(v[1],"phi{1,4}'"));
    property("phi{1,10}",v->value(v[1],2,1));
    property("phi{1,12}",v->value(v[1],2,E(3)));
    property("phi{1,14}",v->isconjugate(v[2],"phi{1,22}") and
                           isconjugate(v[1],"phi{1,10}'"));
    property("phi{1,18}",v->value(v[1],2,E(3)));
    property("phi{2,3}",v->value(v[1],2,-E(3)));
    property("phi{2,5}",v->isgal(v[2],"phi{2,1}",NFAutomorphism(CF(12),5)) and
                    isgal(v[1],"phi{2,13}'",NFAutomorphism(CF(12),5)));
    property("phi{2,7}",v->isgal(v[3],"phi{2,1}",NFAutomorphism(CF(12),7)) and
                    isgal(v[1],"phi{2,13}'",NFAutomorphism(CF(12),7)));
    property("phi{2,9}",v->isconjugate(v[1],"phi{2,15}") and
                           isconjugate(v[3],"phi{2,3}'"));
    property("phi{2,11}",v->isconjugate(v[2],"phi{2,1}") and
                           isconjugate(v[1],"phi{2,13}'"));
    property("phi{2,13}",v->value(v[1],2,-1));
  elif n="G27" then
    property("phi{3,5}",v->fakdegsmaller(v[1],v[2]));
    property("phi{3,20}",v->fakdegsmaller(v[1],v[2]));
    property("phi{8,9}",v->fakdegsmaller(v[2],v[1]));
    property("phi{5,15}",v->inducedfromid(v[1],[1,3]));
    property("phi{5,6}",function(v)local L,t;
      L:=ReflectionSubgroup(W,[1,3]); t:=InductionTable(L,W);
      return ok(t.scalar[f(v[1])][PositionId(L)]=1,
         "phi{5,6}' should occur once in the induced of Id from ",
          ReflectionName(L));end);
  elif n="F4" then # W.rootLengths=[3,1]
    property("phi{1,12}",v->inducedfromid(v[1],[3,4]));
    property("phi{2,4}",v->inducedfromid(v[1],[3,4]));
    property("phi{2,16}",v->issign(v[1],"phi{2,4}''"));
    property("phi{4,7}",v->inducedfromid(v[1],[3,4]));
    property("phi{6,6}",v->occursin(v[2],"phi{4,1}","phi{4,1}"));
    property("phi{8,9}",v->inducedfromid(v[1],[3,4]));
    property("phi{8,3}",v->issign(v[1],"phi{8,9}''"));
    property("phi{9,6}",v->inducedfromid(v[1],[3,4]));
  elif n="G29" then
    property("phi{15,4}",v->occursin(v[2],"phi{4,1}","phi{4,3}"));
    property("phi{15,12}",v->issign(v[2],"phi{15,4}''"));
    property("phi{6,10}",v->occursin(v[3],"phi{4,1}","phi{4,1}") 
      and isconjugate(v[4],v[3]) and inducedfromid(v[1],[1,2,4]));
  elif n="H4" then 
    property("phi{30,10}",v->occursin(v[1],"phi{9,2}","phi{6,12}"));
  elif n="G31" then
    property("phi{15,8}",v->occursin(v[1],"phi{4,1}","phi{20,7}"));
    property("phi{15,20}",v->occursin(v[1],"phi{4,1}","phi{20,7}"));
    property("phi{20,13}",v->isconjugate(v[1],"phi{20,7}"));
    property("phi{30,10}",v->fakdeghighterms(v[1],q^50+q^46));
    property("phi{45,8}",v->occursin(v[2],"phi{4,1}","phi{20,7}"));
    property("phi{45,12}",v->issign(v[1],"phi{45,8}'"));
  elif n="G33" then 
    property("phi{10,8}",v->fakdeghighterms(v[1],q^28+q^26));
    property("phi{10,17}",v->issign(v[1],"phi{10,8}'"));
    property("phi{40,5}",v->fakdeghighterms(v[1],q^31+q^29+2*q^27));
    property("phi{40,14}",v->issign(v[1],"phi{40,5}'"));
  elif n="G34" then
    property("phi{20,33}",v->occursin(v[1],"phi{6,1}","phi{15,14}"));
    property("phi{70,9}",v->occursin(v[2],"phi{6,1}","phi{15,14}")
      and isconjugate(v[1],v[1]));
    property("phi{70,45}",v->issign(v[2],"phi{70,9}''") 
      and isconjugate(v[1],v[1]));
    property("phi{105,8}",v->isconjugate(v[1],"phi{105,4}"));
    property("phi{120,21}",v->inducedfromid(v[1],[1,2,4,5,6]));
    property("phi{280,12}",v->occursin(v[1],"phi{6,1}","phi{336,17}"));
    property("phi{280,30}",v->occursin(v[2],"phi{6,1}","phi{336,17}"));
    property("phi{540,21}",v->occursin(v[1],"phi{6,1}","phi{105,20}"));
    property("phi{560,18}",v->occursin(v[3],"phi{6,1}","phi{336,17}")
      and isconjugate(v[1],v[1]));
    property("phi{840,13}",v->isconjugate(v[1],"phi{840,11}"));
    property("phi{840,23}",v->isconjugate(v[1],"phi{840,19}"));
  elif n in ["G6","G8","G9","G10","G11","G13","G14","G15","G16","G17",
    "G18","G19","G20","G21","G25","G26","G32"] then
    rules();
  else return;
  fi;
  if ddb<>[] then Print("not separated:",ddb,"\n");fi;
  if perm<>() then 
    Print("Permutation to do:\n");
    perm:=Cycles(perm);
    if ForAll(perm,x->Length(x)=2) then
       for t in perm do Print(nm[t[1]], " <=> ",nm[t[2]],"\n");od;
    else Print(perm,"\n");
    fi;
  fi;
end;
