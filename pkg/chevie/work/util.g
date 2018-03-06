# precautions so things work well..
CHEVIE.CheckIndexChars:=true;

############################################################################

# Computes an Hecke algebra with parameters depending on type:
GenericHecke:=function(arg)local W,type,power,r,o,vars,v,i,j,q,p,newp,para;
  vars:="xyztuvwabcdefghijklmnopqrs";
  if Length(arg)=0 then Print(
" GenericHecke(W,type[,power])\n",
" type= 0: [a,b,c]  10: [x0,x1,x2]\n",
"       1: [1,b,c]  11: [1,x1,x2]\n",
"       2: [1,q,q^2] (if several Hplanes take Lcm)\n",
"       3: [q,E3,E3^2] (Spets) 13: [X(Cyclotomics),E3,E3^2]\n",
"       4: [1,E3,E3^2]\n",
"       5: [1,2,3] (primes)\n",
"       6: [a,E(3)*b,E(3)^2*c] 16: [x0,E(3)*x1,E(3)^2*x2]\n",
" if power given the variables are raised to that power\n");
  return;
  fi;
  W:=arg[1];type:=arg[2];
  if Length(arg)=2 then power:=1;else power:=arg[3];fi;
  p:=1;
  newp:=function()
    p:=p+1;
    while not IsPrimeInt(p) do p:=p+1;od;
    return p;
  end;
  r:=Set(W.orbitRepresentative{W.generatingReflections});
  o:=List(r,i->OrderPerm(W.(i)));
  v:=0;
  para:=[];
  for i in [1..Length(r)] do para[r[i]]:=List([0..o[i]-1],function(j)
    if type in [0,10] then 
      if type=10 then return Mvp(SPrint([vars[i]],j))^power;
      else  v:=v+1;   return Mvp([vars[v]])^power;
      fi;
    elif type in [1,11] then
      if j=0 then return Mvp("q")^0;
      elif type=11 then return Mvp(SPrint([vars[i]],j))^power;
      else  v:=v+1;return Mvp([vars[v]])^power;
      fi;
    elif type=2 then return Mvp("q")^(power*j*Lcm(o)/o[i]);
    elif type in [3,13] then
      if j=0 then if type=13 then return X(Cyclotomics)^power;
      else return Mvp("q")^power;fi;
      else return E(o[i])^j;
      fi;
    elif type=4 then return E(o[i])^j;
    elif type=5 then
      if j=0 then return 1;
      else return newp()^power;
      fi;
    elif type in [6,16] then 
      if type=16 then return E(o[i])^j*Mvp(SPrint([vars[i]],j))^power;
      else  v:=v+1;return E(o[i])^j*Mvp([vars[v]])^power;
      fi;
    fi;
  end);
  od;
  return Hecke(W,para);
end;

###############################################################################
#  Various tests of integrity of the data                                     #
###############################################################################

# FindRepresentation(H,r [,true])
# find which char corresponds; false if none; check all char values if 3rd arg
FindRepresentation:=function(arg)local gr,check,W,i,ct,pos,O,l,iszero,r,t,gens;
  W:=arg[1];gr:=arg[2]; iszero:=x->x=0*x;
  check:=Length(arg)>2; O:=W;
  if IsHeckeAlgebra(O) then W:=Group(O);
    t:=List(ChevieClassInfo(W).classtext,x->W.rootRestriction{x});
  elif IsSpets(O) then W:=Group(O);
    t:=List(ChevieClassInfo(O).classtext,x->W.rootRestriction{x});
  elif O.operations=HeckeCosetOps then W:=Spets(O);
    t:=List(ChevieClassInfo(W).classtext,x->Group(W).rootRestriction{x});
  else t:=List(ChevieClassInfo(W).classtext,x->W.rootRestriction{x});
  fi;
  ct:=CharTable(O).irreducibles;
  l:=[1..Length(t)];SortBy(l,i->Length(t[i]));t:=t{l};pos:=[1..Length(l)];
  if IsRec(gr) then gens:=ShallowCopy(gr.gens);Add(gens,gr.F);
    for i in [1..Length(t)] do Add(t[i],Length(gens));od;
  else gens:=gr;
  fi;
  for i in [1..Length(t)] do
    r:=CharRepresentationWords(gens,[t[i]])[1];
    if IsHeckeAlgebra(O) then r:=r*O.unit;fi;
    pos:=Filtered(pos,j->iszero(ct[j][l[i]]-r));
    if not check and Length(pos)=1 then return pos[1];
    elif Length(pos)=0 then return false;
    fi;
  od;
  if Length(pos)=1 then return pos[1];
  elif Length(pos)=0 then return false;
  else Error("characters ",pos," are equal");
  fi;
end;

# CheckRepresentations(H [,l])
# check representations in list l or of dim. l (default: all)
# Check that the representations of the group, coset, Hecke algebra  or
# Hecke coset W:
# 1- satisfy the braid and quadratic relations [and action of F]
# 2- agree with chartable
CheckRepresentations:=function(arg)local W,i,gr,r,ct,pos,H,O,l,f,coset,cl,WF;
  W:=arg[1];O:=W;
  if IsHeckeAlgebra(W) then H:=W;W:=Group(H);coset:=false;
    cl:=ChevieClassInfo(W).classtext;
  elif IsSpets(W) then H:=Hecke(Group(W));coset:=true;WF:=W;
    W:=Group(WF);cl:=ChevieClassInfo(W).classtext;
  elif W.operations=HeckeCosetOps then H:=W.hecke;coset:=true;
    WF:=Spets(W);cl:=ChevieClassInfo(WF).classtext;
  else H:=Hecke(W);coset:=false;
    cl:=ChevieClassInfo(W).classtext;
  fi;
  ct:=CharTable(O).irreducibles;
  if Length(arg)=1 then l:=[1..Length(ct)];
  elif IsInt(arg[2]) then l:=Filtered([1..Length(ct)],i->ct[i][1]=arg[2]);
  else l:=arg[2];
  fi;
  for i in l do
    Print("Representation #",i);
    gr:=Representations(O,i);
    if gr=false then Print("=false\n");
    else
      if coset then r:=gr.gens;
        if gr.F^OrderPerm(WF.phi)<>gr.F^0 then 
	  ChevieErr("F^",OrderPerm(WF.phi),"<>1");
	fi;
	if Permuted(r,WF.phi^-1)<>List(r,x->x^gr.F) then
	  ChevieErr("F does not act as ",
	    RestrictedPerm(WF.phi,W.rootInclusion{W.generatingReflections}),
	     " on generators");
	fi;
      else r:=gr;fi;
      Print(" dim:",Length(r[1]),"...\c");
      CheckHeckeDefiningRelations(H,r);
      pos:=FindRepresentation(O,gr,true);
      if pos=i then Print("\n");
      elif pos<>false then ChevieErr("character found at ",pos,"\n");
      else ChevieErr("character does not match\n");
	pos:=TransposedMat([ct[i],CharRepresentationWords(gr,cl)]);
	f:=List(pos,x->x[1]<>x[2]);
	pos:=List(pos,x->List(x,FormatGAP));
	Print(FormatTable(ListBlist(pos,f),
	  rec(rowLabels:=ListBlist([1..Length(pos)],f),
	      columnLabels:=["is","should be"])));
      fi;
    fi;
  od;
end;

# CheckG4_22Hecke(H[,dimension  of chars to explore])
CheckG4_22Hecke:=function(arg)local W,H,rows,t,v,u,g,name,ix,para,tp;
  H:=arg[1];W:=Group(H);tp:=ReflectionType(W);para:=H.parameter;
  if Length(tp)>1 or not IsBound(tp[1].ST) or not tp[1].ST in [4..22] then
    Error("H should be a Hecke algebra of G4..G22");
  fi;
  tp:=tp[1];para:=CHEVIE.Data("GetParams",tp,para);
  g:=CHEVIE.Data("Generic",tp);
  u:=CHEVIE.Data("EigenvaluesGeneratingReflections",rec(series:="ST",ST:=g));
  u:=Concatenation(List(u,x->List([0..1/x-1],j->E(1/x)^j)));
  # assume non-cyc params are monomials in one Mvp
  t:=TransposedMat([Concatenation(para),u]);
  t:=List(t,function(x)local x,v,p;
    v:=x[2];x:=x[1];
    if not IsMvp(x) or Length(x.elm)<>1 or Length(x.elm[1].elm)<>1 or 
      AsRootOfUnity(x.coeff[1])=false then return false;fi;
    p:=(AsRootOfUnity(v/x.coeff[1])+
        [0..AbsInt(x.elm[1].coeff[1])-1])/x.elm[1].coeff[1];
    return [x.elm[1].elm[1],List(p,x->E(Denominator(x))^Numerator(x))];
  end);
  t:=Filtered(Set(t),x->x<>false);
  v:=Set(List(t,x->x[1]));
  t:=List(v,x->[x,Intersection(List(Filtered(t,y->y[1]=x),y->y[2]))[1]]);
  Print("specialization=",Join(List(t,x->SPrint(x[1],"->",x[2]))),"\n");
  rows:=List(CharTable(Hecke(ComplexReflectionGroup(g),para)).irreducibles,
    x->x{ChevieClassInfo(W).indexclasses});
  u:=List(rows,x->List(x,function(y)if IsCyc(y) then return y;
                 else return ScalMvp(Value(y,Concatenation(t)));fi;end));
  t:=CHEVIE.Data("CharTable",tp).irreducibles;
  Print("how chars of G",g," for these params specialize to chartable\n");
  v:=List(u,x->Position(u,x));
  v:=List(Set(v),y->[y,Filtered([1..Length(v)],i->v[i]=y)]);
  v:=Set(List(v,x->[Position(t,u[x[1]]),x[2]]));
  u:=Maximum(List(v,x->Length(x[2])));
  for ix in v do Append(ix[2],List([1..u-Length(ix[2])],i->""));od;
  ix:=CHEVIE.Data("CharInfo",tp).indexchars;
  v:=List(v,function(x)local ok;ok:="no";
    if x[1]<>false then if ix[x[1]] in x[2] then ok:="";fi;
         return Concatenation([x[1],ix[x[1]],ok],x[2]);
    else return Concatenation([false,false,ok],x[2]);fi;end);
  if Length(arg)=2 then v:=Filtered(v,x->t[x[1]][1]=arg[2]);fi;
  Print(FormatTable(TransposedMat(v),
     rec(rowLabels:=Concatenation(["char","chosen","ok"],List(v[1],x->"")))));
  t:=Representations(H);
  u:=List(t,r->FindRepresentation(H,r));
  if u<>[1..Length(u)] then 
    u:=TransposedMat([[1..Length(u)],u]);
    u:=Filtered(u,x->x[1]<>x[2]);
    u:=TransposedMat(u);
    Print("representations do not match:\n");
    Print(FormatTable([u[2]],rec(rowLabels:=["found at"],columnLabels:=u[1])));
  fi;
  return v; # list of vectors [n,chosen,n in poss,poss1,..,possr]
end;

# shows info about which chars cause problems for which algebras for G_n
# GuessG4_22(n[,types])
GuessG4_22:=function(arg)local n,t,l, ll,opt;
  n:=arg[1];if Length(arg)=2 then t:=arg[2];else t:=[0..4];fi;
  l:=List(t,i->CheckG4_22Hecke(GenericHecke(ComplexReflectionGroup(n),i)));
  for ll in l do SortBy(ll,x->x[1]);od;
  l:=TransposedMat(l);
  l:=List(l,function(v)local res,p,j;
    if Length(Set(v))=1 and v[1][3]<>"no" then return true;fi;
    res:=[v[1][1],v[1][2]];
    for p in Set(v) do
      j:=Filtered([1..Length(v)],i->v[i]=p);
      Add(res,v[j[1]]{[4..Length(v[j[1]])]});
      Add(res,t{j});
    od;
    return res;
  end);
  l:=Filtered(l,x->x<>true);
  for t in l do
    Print(t[1],":chosen ",t[2]);
    if ForAll([3,5..Length(t)-1],i->t[2] in t[i]) then Print(" OK!!");
    else for n in [3,5..Length(t)-1]do 
      Print(" ",IntListToString(t[n]),"for types",IntListToString(t[n+1]));od;
    fi;
    Print("\n");
  od;
  return l; 
end;

# Check SchurElements(H) satisfy Schur relations
CheckSchurRelations:=function(H)local i,s,p,ct,lcm,n,un;
  un:=Product(Flat(H.parameter),x->x^0);
  n:=NrConjugacyClasses(Group(H));
  ct:=CharTable(H);
  if ct=false then return;fi;
  if IsMvp(un) then
    s:=FactorizedSchurElements(H);
    lcm:=FactorizedSchurElementsOps.Lcm(s);
    s:=List(s,x->lcm/x);
    Print("expanding Lcm(S_\\chi)/S_\\chi quotients..\c");
    s:=List(s,Expand)*un;
    lcm:=Expand(lcm);
    Print("done in ", Stime(),"\n");
  else
    s:=SchurElements(H);
    lcm:=s[PositionId(Group(H))];
    s:=List(s,x->lcm/x);
  fi;
  for i in [1..n] do
    p:=s*ct.irreducibles{[1..n]}[i];
    if i=1 and p<>lcm then
      ChevieErr("Sum_chi chi(1)/S_chi not 1\n");
    elif i<>1 and p<>0*p then
      ChevieErr("Sum_chi chi(",Ordinal(i)," class)/S_chi not 0\n");
    else Print(".\c");
    fi;
  od;
  Print(" done in ", Stime(),"\n");
# Ok("generic*CharTable(H)");
end;

CheckHeckeSpecializes:=function(W)local H,i;
  Unbind(W.HeckeAlgebras);
  CHEVIE.checkroots:=true;
  for i in [0..4] do
    H:=GenericHecke(W,i);
    Print("parameters=",H.parameter,"\n");
    CharTable(H);
  od;
end;
  
CheckOpdam:=function(W)local ct,complex,x,fd,ci,c,i,f,opdam;
  ct:=CharTable(W).irreducibles;
  complex:=PermListList(ct,ComplexConjugate(ct));
  x:=Mvp("x");
  fd:=FakeDegrees(W,x);
  ci:=ChevieCharInfo(W);
  if IsBound(ci.opdam) then opdam:=ci.opdam;else opdam:=();fi;
  # c_\chi in Gordon-Griffeth
  c:=function(W,i)local cl;
    cl:=ChevieClassInfo(W).classes;
    return Sum(ReflectionDegrees(W)-1)-Sum(HyperplaneOrbits(W),
      h->Sum(h.classno,c->cl[c]*ct[i][c]/ct[i][1]));
  end;
  for i in [1..NrConjugacyClasses(W)] do
    f:=fd[i^complex];
    if fd[i^opdam]<>x^c(W,i)*Value(f,["x",x^-1]) then Print("**** ",i,"\n");fi;
  od;
end;

goodlines:=function(f)local m,ind,special,ss;
  m:=f.fourierMat;ind:=[1..Length(m)];
  if IsBound(f.signs) then m:=m^DiagonalMat(f.signs);fi;
  special:=Filtered(ind,x->not 0 in m[x]);
  m:=TransposedMat(m);
  if Length(m)>15 then Print("lines to check:",special,"\n");fi;
  ss:=Filtered(special,function(ss)local s;
    s:=List(m,x->x/x[ss]);Print(ss," \c");
    return ForAll(ind,i->ForAll(ind,j->ForAll(List(ind,k->m[k][i]*m[k][j])*s,
       res->IsInt(res) and res>=0)));end);
  if ss<>special then 
   Print("******* of non-0 lines:",special," only ",ss," are positive\n");fi;
  return ss;
end;

allfams:=function()
 return Concatenation(List(CHEVIE.TestObjs,function(W)local n,f,ud,i;
   W:=ApplyFunc(W[1],W{[2..Length(W)]});
   f:=UnipotentCharacters(W);
   if f=false then return [];fi;
   f:=f.families; n:=ReflectionName(W);
#  if IsCoxeterGroup(W) or IsCoxeterCoset(W) then return [];fi;
   ud:=CycPolUnipotentDegrees(W);
   for i in [1..Length(f)] do
     f[i].label:=SPrint(n,".",i); f[i].ud:=ud{f[i].charNumbers};
     f[i].operations:=ShallowCopy(f[i].operations);
     f[i].operations.Print:=function(f)
       if not IsBound(f.name) then f.name:="?";fi;
       Print(f.label,"#",Length(f.charNumbers),"=",f.name);end;
   od;
   return Filtered(f,x->Length(x.charNumbers)>1);
   end));
end;
   
testmu:=function(x)local l,r;
  if IsList(x) then return ForAll(x,testmu);fi;
  if x=0 then return true;fi;
  l:=Conjugates(x);
  r:=RootInt(1/Product(l),Length(l));
  if not IsInt(r) then return false;fi;
  if 1/Product(l)<>r^Length(l) then return false;fi;
  return AsRootOfUnity(x*r)<>false;
end;
