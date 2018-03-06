ReadChv("work/init");
ReadChv("work/getunpdeg");
ReadChv("work/mvptools");
# programmes accompagnant le bilan de Luminy 2017
CHEVIE.CheckIndexChars:=true;

KnownSeries:=function(W)
  if not IsBound(W.series) then W.series:=getunpdeg(W)[1];fi; 
  return W.series;
end;

# Calcule la (q,t)-Trace(D_x | H^*_c(X(y))) pour y = \pi^ny et x =\pi^nx
# suivant la formule
# \sum_{\chi in \Irr(H_y)]\chi(x) exp(2i\pi n_x n_y \alpha_\chi)
#    q^{n_x(N+N^*-\alpha_\chi)} t^{n_y(N+N^*-\alpha_\chi)}
# ou \alpha_\chi=a_\chi+A_\chi
hxy:=function(arg) # hxy(W,nx,ny[,opt])
# W groupe de ref complexe, nx et nx sont des fractions
# opt: si absent retourne un UnipotentCharacter
#      si opt.var présent retourne un polynome en variables u_{y,\chi}
#      si opt.sh présent retourne un polynome en variables su_{y,\chi}
  local W,nx,ny,s,cx,q,t,signes,chix,aA,NN,asroot,u,opt,cy,oy,eig,i,ud,ct;
  W:=arg[1];nx:=arg[2];ny:=arg[3];
  if Length(arg)=4 then opt:=arg[4];else opt:=rec();fi;
  asroot:=function(eta)
    return E(Denominator(eta))^(Numerator(eta));
  end;
  s:=First(KnownSeries(W),x->x.d=Mod1(ny)); # utilise H^*(X(y))=H^*(X(y\pi))
  if not (Denominator(nx) in RegularEigenvalues(s.WGL)) then 
    #Error(nx," is not a regular number for ",s.WGL);
    return false;
  fi;
  q:=Mvp("q");t:=Mvp("t");
  ud:=UnipotentDegrees(W,q){s.charNumbers};
  signes:=List(ud,x->SignInt(Value(x,["q",asroot(ny)])));
  ct:=CharTable(s.WGL).irreducibles;
  cx:=PositionRegularClass(s.WGL,Mod1(nx));
  chix:=List(ct,chi->chi[cx]);
  aA:=List(ud,P->Degree(P)+Valuation(P));
  NN:=Sum(ReflectionDegrees(W))+Sum(ReflectionCoDegrees(W));
  if IsBound(opt.var) then
    if not IsBound(W.hxy) then W.hxy:=rec(eigen:=rec());fi;
    u:=rho->Mvp(SPrint("u",Mod1(ny),".",rho)); 
    cy:=PositionRegularClass(s.WGL,Mod1(ny));
    oy:=List(ct,chi->chi[cy]/chi[1]);
    eig:=Zip(oy,aA,function(o,a)return o*asroot(ny^2*a);end);
    for i in [1..Length(ud)] do
      W.hxy.eigen.(SPrint("u",Mod1(ny),".",i)):=eig[i];
    od;
  elif IsBound(opt.sh) then u:=rho->Mvp(SPrint("su",ny,".",rho));
  else u:=rho->UnipotentCharacter(W,s.charNumbers[rho]);
  fi;
  return Sum([1..Length(ud)],i->chix[i]*asroot(PowerRoot(ny,aA[i]*nx))*
    q^((NN-aA[i])*nx)*t^((NN-aA[i])*ny)*signes[i]*u(i));
end;

Shintani:=function(W,u) #calcule l'effet de shintani sur un  UnipotentCharacter
# à coefficients dans C[q,q^{-1},t,t^{-1}] (change q en qt)
  local uc,S,O,Sh;
  uc:=UnipotentCharacters(W);
  S:=Fourier(uc);
  O:=DiagonalMat(Eigenvalues(uc));
  if ForAny(uc.families,f->IsBound(f.lusztig)) then
    Sh:=O*S^(-1)*O^(-1);
  else
    Sh:=O*S^(-1)*O;
  fi;
  return UnipotentCharacter(W,Value(u.v,["q",Mvp("q")*Mvp("t")])*Sh);
end; 
  
UnipCharOps.Projection:=function(u,L)
  return UnipotentCharacter(u.group,List([1..Length(u.v)],
     function(i) if i in L then return u.v[i]; else return 0; fi; end));
end;

UnambigFamilies:=function(W)local N,F;
  N:=Filtered(KnownSeries(W),s->IsBound(s.ambig));
  N:=Union(List(N,s->Union(List(s.ambig,x->x[2]))));
  F:=UnipotentCharacters(W).families;
  F:=Filtered(F,f->Intersection(N,f.charNumbers)=[]);
  return Union(List(F,f->f.charNumbers));
end;

AdmissiblePairs:=function(W)
  local D,res;
  D:=List(KnownSeries(W),s->s.d);
  res:=List(Cartesian(D,D),i->[Mod1(i[2]-i[1]),i[1]]);
  return res;
  return Filtered(res,x->x[2] in D);
end;

g:=function(W,i,j)
  local Un,a,b;
  Un:=UnambigFamilies(W);
  a:=hxy(W,i,j);
  if a=false then return false;fi;
  b:=hxy(W,i+j,j);
  if b=false then return false;fi;
  return Projection(b-Shintani(W,a),Un);
end;
  
# Coefficients of a polynomial for each monomial in q,t
explode2:=function(p)
  p:=Zip(p.elm,p.coeff,function(e,c)local pt,ct,pq,cq,left;
    pt:=Position(e.elm,"t");if pt=false then ct:=0;else ct:=e.coeff[pt];fi;
    pq:=Position(e.elm,"q");if pq=false then cq:=0;else cq:=e.coeff[pq];fi;
    left:=Filtered([1..Length(e.elm)],i->i<>pt and i<>pq);
    return [c,rec(elm:=e.elm{left},coeff:=e.coeff{left}),[ct,cq]];
    end);
  p:=CollectBy(p,x->x[3]);
  return List(p,x->Mvp(List(x,x->x[2]),List(x,x->x[1])));
end;

# Cut a unipotent character according to each monomial in q,t
explode:=function(u)local mm,l;
  l:=u.v*Mvp("x")^0;
  mm:=Union(List(l,x->x.elm));
  return List(mm,m->List(l,function(P)local p;
    p:=Position(P.elm,m);
    if p=false then return 0;
    else return P.coeff[p];fi;end));
end;

# Cut a polynomial in u_{y,\chi} according to each eigenvalue of u_{y,\chi}
explodeEigen:=function(W,p)local e,eig;
  e:=Mvp(p[1])-p[2];
  eig:=List(Variables(e),v->W.hxy.eigen.(v));
  eig:=CollectBy(Variables(e),eig);
  return List(eig,function(c)local s;
     s:=Filtered([1..Length(e.elm)],i->e.elm[i].elm[1] in c);
     return Mvp(e.elm{s},e.coeff{s});
     end);
end;

# utilise omega^-1*Sh*Omega^-1=Sh*Omega^-1*Sh
usebraid:=function(W,p)local og,rm,lm,a,b,vv,v,l;
  l:=W.hxy.sols;
  og:=W.hxy.eigen.(p[1]{[2..Length(p[1])]});
  rm:=Mvp(p[2].elm,Zip(p[2].coeff,p[2].elm,
    function(c,e)return c*og^-1*W.hxy.eigen.(e.elm[1])^-1;end));
  a:=List(p[2].elm,e->rec(coeff:=[1],elm:=[ConcatenationString("s",e.elm[1])]));
  b:=Zip(p[2].coeff,p[2].elm,
    function(c,e)return c*W.hxy.eigen.(e.elm[1])^-1;end);
  lm:=Mvp(a,b);
  rm:=rm-lm;
  for v in Filtered(Variables(rm),x->x[1]='s') do
    vv:=First(l,x->x[1]=v);
    rm:=Value(rm,vv);
  od;
  return rm;
end;

poltovec:=function(W,p)local v,res,i;
  v:=RecFields(W.hxy.eigen);
  res:=List(v,x->0);
  for i in [1..Length(p.elm)] do
    res[Position(v,p.elm[i].elm[1])]:=p.coeff[i];
  od;
  return res;
end;
  
shmatrix:=function(W)local v,m,l;
  l:=W.hxy.sols;l:=Filtered(l,x->x[1][1]='s');
  v:=RecFields(W.hxy.eigen);
  if Length(l)<>Length(v) then Error("not all variables s... found");fi;
  m:=List(l,e->poltovec(W,e[2]));
  SortParallel(List(l,x->Position(v,x[1]{[2..Length(x[1])]})),m);
  return m;
end;

omegamatrix:=function(W)
  return DiagonalMat(List(RecFields(W.hxy.eigen),x->W.hxy.eigen.(x)));
end;

testsh:=function(W)local l,a,b,p,v,addsol;
  addsol:=function(v)local p;
    Add(W.hxy.sols,v);for p in W.hxy.sols do p[2]:=Value(p[2],v);od;
  end;
  l:=[];
  for p in AdmissiblePairs(W) do
    a:=hxy(W,p[1],p[2],rec(sh:=true));
    if a<>false then
      a:=Value(a,["q",Mvp("q")*Mvp("t")]);
      b:=hxy(W,p[1],p[1]+p[2],rec(var:=true));
      if b<>false then Add(l,b-a);fi;
    fi;
  od;
  W.hxy.sols:=[];
  l:=Concatenation(List(l,explode2));
  while Length(l)>0 do
    v:=solve(l[1]);addsol(v);l:=Value(l,v);l:=sortbylg(l);
  od;
  l:=Filtered(W.hxy.sols,x->x[1] in RecFields(W.hxy.eigen) and 
    Length(x[2].coeff)>1);
  l:=Concatenation(List(l,p->explodeEigen(W,p)));
  while Length(l)>0 do
    v:=solve(l[1]);addsol(v);l:=Value(l,v);l:=sortbylg(l);
  od;
  p:=PositionProperty(W.hxy.sols,x->ForAny(Variables(x[2]),
    v->v[1]='s'));
  if p<>false then 
    Error("one of the s.. variables not resolved",W.hxy.sols[p],"\n");
  fi;
  for p in Filtered(W.hxy.sols,x->x[1][1]='s') do
    v:=usebraid(W,p); if v<>0 then v:=solve(v);addsol(v);fi;
  od;
  W.hxy.sols:=Set(W.hxy.sols);
  return W.hxy.sols;
end;
