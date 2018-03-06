# This file contains code to find Ennola-permutation,
# and code to find Fourier matrices (detfam)
#------------------ EnnolaOps -----------------------------------------
# Ennola-theory axioms:
# Let z=|ZW|. We call xi-Ennola the map q -> E(z)^xi q.
# (This corresponds to Ennola_E(z)^-1 of "towards Spetses I").
# - For any character sheaf A there exists e_A such that xi-Ennola(A)=e_A^xi A
# Thus  Ennola stabilizes  families and almost-HC-series; 
# -  If  A  has  cuspidal  data  (\lambda,\chi) where \lambda is a cuspidal
# character  sheaf  of  L  and  \chi\in\Irr(W_G(L,\lambda)),  there  exists
# E_\lambda such that e_A=E_\lambda*\omega_\chi(E(z))
#
# For  a split Spetses,  cuspidal character sheaves  correspond to cuspidal
# unipotent characters; if U\in\Uch(G) has cuspidal data
# \chi\in\Irr(W_G(L,\lambda)) then we have
# Frob(Ennola_xi(U))/Frob(U)=omegachi(E(z))^-xi E(z^2)^(xi^2*(a+A))E_\lambda^xi

EnnolaOps:=OperationsRecord("EnnolaOps");

# possible action of xi-Ennola on i-th family (here xi in Q/Z)
# (list of possible destinations for each char, taking just degree in account)
EnnolaBete:=function(W,i,xi)local ud;
  xi:=E(Denominator(xi))^Numerator(xi);
  ud:=CycPolUnipotentDegrees(W){UnipotentCharacters(W).families[i].charNumbers};
  return List(ud,p->PositionsSgn(ud,CycPolOps.EnnolaTwist(p,xi)));
end;

EnnolaOps.init:=function(W)local e,uc;
  e:=rec(W:=W,scalars:=[[1]]);
  uc:=UnipotentCharacters(W);
  e.z:=OrderCenter(W);
  e.families:=List([1..Length(uc.families)],i->rec(poss:=
    List([1..e.z-1],j->EnnolaBete(W,i,j/e.z))));
  e.hc:=List([1..Size(uc)],n->PositionProperty(uc.harishChandra,
            h->n in h.charNumbers));
  e.operations:=EnnolaOps;
  return e;
end;

EnnolaOps.Display:=function(e,o)
  o.rowLabels:=List(UnipotentCharacters(e.W).harishChandra,
    h->TeXStrip(h.cuspidalName,o));
  o.rowsLabel:="Name";
  o.columnLabels:=[1/e.z];
  o.screenColumns:=SizeScreen()[1];
  Print("Ennola scalars for cuspidals of ",ReflectionName(e.W,o),"\n",
    FormatTable(List(e.scalars,x->[x]),o));
end;

# Ennola for i-th family as an Ls
EnnolaOps.PossToLs:=function(e,i,xi)local r,gg;
  r:=Copy(e.families[i].poss[xi]);
  if Product(r,Length)>1000000 then return;fi;
  gg:=function(l,poss)local x,res;
    if Length(poss)=0 then return [[]];fi;
    res:=[];
    for x in poss[1] do
      if not AbsInt(x) in l then
	Append(res,List(gg(Union(l,[AbsInt(x)]),poss{[2..Length(poss)]}),
	  v->Concatenation([x],v)));
      fi;
    od;
    return res;
  end;
  return gg([],r);
end;

# if Ennola_xi(U_i)=U_j returns deduced En(\lambda)^xi for cuspidal of U_i
EigenEnnola:=function(W,i,j,xi)local uw,z;
  z:=OrderCenter(W);
  uw:=UnipotentCharacters(W);
  return OmegaChi(uw)[i]^-xi*E(z^2)^(xi^2*(uw.a[i]+uw.A[i]))*
    Eigenvalues(uw,[i])[1]/Eigenvalues(uw,[j])[1];
end;

# in e set xi-ennola scalar of cuspidal corresp. to charno to be in scal
EnnolaOps.SetScal:=function(e,scal,charno,xi)local k,n;
  k:=e.hc[charno];
  if xi=1 and not IsBound(e.scalars[k]) then e.scalars[k]:=scal;
  else 
    n:=Filtered(e.scalars[k],x->x^xi in scal);
    if e.scalars[k]<>n then 
      InfoChevie("\n   # ",xi,"-Ennola scalars[",k,"] ",FormatGAP(e.scalars[k]),
         "->",FormatGAP(n));
    fi;
    e.scalars[k]:=n;
  fi;
end;
  
EnnolaOps.GetScal:=function(e,charno,xi)local k,n;
  return Set(List(e.scalars[e.hc[charno]],x->x^xi));
end;

# restrict ennola-scalars from frobenius eigenvalues using
# knowledge on ennola-ps and
# En(\lambda)^k=EigenEnnola(W,i,j,k) if Ennola_xi(U_i)=U_j
EnnolaOps.ScalarsFromFrob:=function(e)local W,uc,i,f,j,xi,z;
  W:=e.W;uc:=UnipotentCharacters(W).families;z:=OrderCenter(W);
  for i in [1..Length(e.families)] do f:=uc[i].charNumbers;
    for xi in [1..z-1] do 
      for j in [1..Length(f)] do
	EnnolaOps.SetScal(e,Set(List(e.families[i].poss[xi][j],
	  x->EigenEnnola(W,f[j],f[AbsInt(x)],xi))),f[j],xi);
      od;
    od;
  od;
end;

# get info on ennola-ps from ennola-scalars and frobenius-eigenvalues
EnnolaOps.PossFromScal:=function(e)local W,i,r,f,j,l,s,xi;W:=e.W;
  for i in [1..Length(e.families)] do 
    f:=UnipotentCharacters(W).families[i].charNumbers;
    for xi in [1..e.z-1] do r:=e.families[i].poss[xi];
      for j in [1..Length(f)] do
	s:=EnnolaOps.GetScal(e,f[j],xi);
	l:=Filtered([1..Length(r[j])],
	  n->EigenEnnola(W,f[j],f[AbsInt(r[j][n])],xi)in s);
	if Length(l)<Length(r[j]) then
	  InfoChevie("   # Family#",i," En",xi/e.z,"[",j,"]:restricted ",
	    Join(r[j]),"->",Join(r[j]{l}),"\n");
	fi;
	r[j]:=r[j]{l};
      od;
    od;
  od;
end;

# EnnolaOps.ScalarsFromFourier(e,fno[,fourierMat)
EnnolaOps.ScalarsFromFourier:=function(arg)
  local e,fno,f,F,n,xi,j,k,l,m,poss,es,o,F1,uc;
  e:=arg[1];fno:=arg[2];
  uc:=UnipotentCharacters(e.W);
  f:=uc.families[fno];
  o:=OmegaChi(uc){f.charNumbers};
  if Length(arg)=2 then
  F:=ComplexConjugate(f.fourierMat);
  if IsBound(f.lusztig) and f.lusztig then 
    F:=Permuted(F,PermListList(ComplexConjugate(f.eigenvalues),f.eigenvalues));
  fi;
  else F:=ComplexConjugate(arg[3]);
  fi;
  n:=Length(F);
  for xi in [1..e.z-1] do
    poss:=e.families[fno].poss[xi];
    es:=List(f.charNumbers,n->EnnolaOps.GetScal(e,n,xi));
    F1:=TransposedMat(Zip(TransposedMat(F),o,function(x,y)return x*y^xi;end));
    for m in [1..n] do
      poss[m]:=Filtered(poss[m],j->ForAll([1..n],
        function(k)local a,b;
	  a:=F[AbsInt(j)][k]; b:=F1[m][k];
	  if b=0 then a:=ScalMvp(a);return a=false or a=0;fi;
	  a:=ScalMvp(a/b);
	  if a=false then return true;fi;
          return SignInt(j)*a in es[k];end));
      if Length(poss[m])=0 then return false;fi;
      if Length(poss[m])=1 then
        j:=poss[m][1];
        l:=Filtered([1..n],k->F[m][k]<>0);
        for k in l do
          EnnolaOps.SetScal(e,[SignInt(j)*F[AbsInt(j)][k]/F1[m][k]],
           f.charNumbers[k],xi);
        od;
      fi;
    od;
  od;
  return true;
end;

# EnnolaFromFusion(W,i[,poss])
# possibilities of Ennola from fusion algebra. If poss given
# tries to fit with poss.
EnnolafromFusion:=function(arg)local W,i,poss,A,b,res;
  W:=arg[1];i:=arg[2];
  if Length(arg)=3 then poss:=arg[3];fi;
  A:=FusionAlgebra(UnipotentCharacters(W).families[i]);
  b:=A.basis;
  res:=List(Concatenation([1..Length(b)],-[1..Length(b)]),function(i)local l;
    if i<0 then l:=-b[-i]*b;else l:=b[i]*b;fi;
    if ForAny(l,x->Length(x.coefficients)<>1 or not
      x.coefficients[1][1] in [1,-1]) then return false;fi;
    return [i,List(l,x->SignInt(x.coefficients[1][1])*x.coefficients[1][2])];
    end);
  res:=Filtered(res,x->x<>false);
  if not IsBound(poss) then return res;fi;
  A:=Copy(poss);
  poss:=List(Set(poss),function(p)local r;
    r:=Filtered(res,x->ForAll([1..Length(x[2])],j->x[2][j] in p[j]));
    if Length(r)=1 then return r[1];fi;
    if Length(r)=0 then return false;
    else Error("too many sols");fi;end);
  if ForAll(poss,x->x=false) then poss:=A;
    ChevieErr("no solution for family ",i,"\n");
    poss:=List(Set(poss),function(p)local r;
      r:=Filtered(res,x->Number([1..Length(x[2])],j->x[2][j] in p[j])
        >=Length(x[2])-2);
      if Length(r)=1 then return r[1];fi;
      if Length(r)=0 then return false;
      else Error("too many sols");fi;end);
    poss:=Filtered(poss,x->x<>false);
  fi;
  return Filtered(poss,x->x<>false);
end;
    
Ennola:=function(W)local e,i,ff,le,fam;
  if OrderCenter(W)=1 then Print("W has trivial center!\n");return false;fi;
  e:=EnnolaOps.init(W);
  EnnolaOps.ScalarsFromFrob(e);
  EnnolaOps.PossFromScal(e);
  fam:=UnipotentCharacters(W).families;
  for i in [1..Length(fam)] do EnnolaOps.ScalarsFromFourier(e,i);od;
  EnnolaOps.PossFromScal(e);
  ff:=Filtered([1..Length(fam)],i->
    ForAny(e.scalars{e.hc{fam[i].charNumbers}},j->Length(j)>1));
  le:=List(Cartesian(e.scalars),s->rec(W:=e.W,z:=e.z,
    families:=Copy(e.families),scalars:=List(s,x->[x]),hc:=e.hc));
  le:=Filtered(le,e->ForAll(ff,i->EnnolaOps.ScalarsFromFourier(e,i)));
  return List(le,function(e)local res,f;
    res:=rec(scalars:=List(e.scalars,x->x[1]));
    res.fls:=Zip([1..Length(fam)],e.families,
      function(i,r)local ls,gls;
        r.poss:=EnnolafromFusion(W,i,r.poss);
        if Length(r.poss)=0 then Error("no solution for family ",i);fi;
        if Length(r.poss)>1 then 
          InfoChevie("   # ",Length(r.poss)," sols in family ",i,"\n");fi;
        return r.poss[1][2];
       end);
    res.ls:=[];for i in [1..Length(fam)] do res.ls{fam[i].charNumbers}:=
      List(res.fls[i],function(x)
        if x<0 then return -fam[i].charNumbers[-x];
	       else return fam[i].charNumbers[x];fi;
	end);od;
    return res;end);
end;

isz:=x->x=0*x;
# omega_chi(E(z)) where z=OrderCenter(W) for all character sheaves (rho,chi)
OmegaChi:=function(uc)local relativeFake,xi0,h;
  if not IsBound(uc.omegachi) then
    relativeFake:=[];
    for h in uc.harishChandra do
      relativeFake{h.charNumbers}:=FakeDegrees(ApplyFunc(ReflectionGroup,
       List(h.relativeType,function(t)t.operations:=ReflTypeOps;return t;end)),
         X(Cyclotomics));
    od;
    xi0:=E(OrderCenter(Group(uc.group)));
    uc.omegachi:=List(relativeFake,f->Value(f,xi0)/Value(f,1));
  fi;
  return uc.omegachi;
end;

# convert family-wise ennola ps to global perm
famEn2En:=function(W,e)local ff,res,i,j;
  ff:=UnipotentCharacters(W).families;res:=[];
  for i in [1..Length(ff)] do
    res{ff[i].charNumbers}:=SignPermuted(ff[i].charNumbers,e[i]);
  od;
  return res;
end;

# Twist a subspets by xi, an element of the center of the (untwisted) parent
EnnolaTwist:=function(HF,xi)local H,W,w;
  H:=Group(HF);W:=Parent(H);
  w:=Representative(ConjugacyClasses(W)[PositionRegularClass(W,AsRootOfUnity(xi))]);
  return SubSpets(HF.parent,H.rootInclusion{H.generatingReflections},w*HF.phi);
end;

# Predict what should be LusztigInductionTable(HF,WF) by Ennola-twisting that
# of EnnolaTwist(HF,xi^-1);
PredictRLGByEnnola:=function(HF,WF,xi)local H,t,ps,GlobalEnnolaLs;
  GlobalEnnolaLs:=function(W,xi)local o;o:=AsRootOfUnity(xi)*OrderCenter(W);
    return famEn2En(W,List(Ennola(W)[1].ps,x->LsPuiss(x,o)));
  end;
  HF:=EnnolaTwist(HF,xi^-1);
  H:=Group(HF);	
  t:=LusztigInductionTable(HF,WF).scalar;
  if xi^OrderCenter(H)=1 then
    ps:=GlobalEnnolaLs(H,xi);
    t:=Permuted(TransposedMat(t),PermList(List(ps,AbsInt))^-1);
    t:=TransposedMat(List([1..Length(t)],x->SignInt(ps[x])*t[x]));
  fi;
  ps:=GlobalEnnolaLs(Group(WF),xi);
  t:=Permuted(t,PermList(List(ps,AbsInt))^-1);
  t:=List([1..Length(t)],x->SignInt(ps[x])*t[x]);
  return t;
end;

CheckRLGByEnnola:=function(HF,WF,xi)local H,t,fW,fH,pieces,p,s,i,r,j,k;
  H:=Group(HF);	
  t:=PredictRLGByEnnola(HF,WF,xi);
  fW:=UnipotentCharactersOps.Fourier(UnipotentCharacters(WF));
  fH:=UnipotentCharactersOps.FourierInverse(UnipotentCharacters(HF));
  t:=fW*t*fH;
  pieces:=LusztigInductionTable(HF,WF,true);
  Print("Checking ",pieces,"\n");
  pieces:=pieces.pieces;
  for i in [1..Length(pieces)] do
    p:=pieces[i];
    p.ns:=t{p.wnum}{p.hnum};
    t{p.wnum}{p.hnum}:=0*t{p.wnum}{p.hnum};
    s:=ProportionalityCoefficient(Concatenation(p.scalar),Concatenation(p.ns));
    if s<>false then
      if s<>1 then Print("got piece[",i,"] upto ",s,"\n");fi;
    else
      r:=TransposedMat(p.scalar);
      p.ns:=TransposedMat(p.ns);
      r:=List([1..Length(r)],i->ProportionalityCoefficient(r[i],p.ns[i]));
      if ForAny(r,x->x=false) then
        Error("extensions not proportional",r);fi;
      Print(" *** piece no.",i," of ",Length(pieces)," ***\n");
      Print(p.head(p,rec()));
#     Display(p);
      j:=Filtered([1..Length(r)],i->r[i]<>1);
      for k in j do
        Print("char. ",k,"=",p.charSubgroup[k]," of ",p.u," extension differs by ",r[k],"\n");
      od;
    fi;
  od;
  if Set(Concatenation(t))<>[0] then Error();fi;
end;

# Check that an inner Ennola E_xi sends RLG to RE_xi(L)E_xi(G)
CheckEnnolaRLG:=function(W)
  local HF,WF,p,w,rW,HFxi,rH,t,txi,l,xi,c,n,nxi,twist,e,uc,ps,psH,i;
  twist:=function(HF,w)local r;
    r:=Group(HF).rootInclusion{Group(HF).generatingReflections};
    return SubSpets(HF.parent,r,w*HF.phi);
  end;
  WF:=Spets(W);
  c:=OrderCenter(W);
  l:=Concatenation(List(ProperParabolics(W),x->Twistings(W,x)));
  e:=Ennola(W)[1].ps;
  uc:=List(UnipotentCharacters(W).families,x->x.charNumbers);
  ps:=[];
  for p in [1..Length(uc)] do
    ps{uc[p]}:=List(e[p],x->SignInt(x)*uc[p][AbsInt(x)]);
  od;
  for xi in [1..c-1] do
    for p in Filtered(List(l,x->[x,EnnolaTwist(x,E(c)^xi)]),x->x[1]>=x[2]) do
    HF:=p[1];HFxi:=p[2];
    n:=ReflectionName(HF);nxi:=ReflectionName(HFxi);
# if '[' in  n or '[' in nxi or n{[1,2]}="2A" or nxi{[1,2]}="2A" then
#   Print(E(c)^xi,":",ReflectionName(HF),"=>",ReflectionName(HFxi)," ***Cannot test***!\n");
# else
    e:=Ennola(Group(HF))[1].ps;
    uc:=List(UnipotentCharacters(Group(HF)).families,x->x.charNumbers);
    psH:=[];
    for i in [1..Length(uc)] do
      psH{uc[i]}:=List(e[i],x->SignInt(x)*uc[i][AbsInt(x)]);
    od;
    t:=TransposedMat(LusztigInductionTable(HF,WF).scalar);
    t:=Permuted(t,PermList(List(psH,AbsInt))^-1);
    t:=TransposedMat(List([1..Length(t)],x->SignInt(psH[x])*t[x]));
    t:=Permuted(t,PermList(List(ps,AbsInt))^-1);
    t:=List([1..Length(t)],x->SignInt(ps[x])*t[x]);
    txi:=LusztigInductionTable(HFxi,WF);
    if txi.scalar<>t then
      p:=Copy(txi);p.scalar:=t;
      CHEVIE.Check.EqTables(p,txi);
      Error();
    else
      Print(E(c)^xi,":",ReflectionName(HF),"=>",ReflectionName(HFxi)," OK!\n");
    fi;
  # fi;
    od;
  od;
end;

# The Fourier matrix F times unipotent  degrees are the fake degrees. If
# p=Ennola-Permutation  and e=Ennola-scalars  on character  sheaves then
# one has Fp=eF
# This program,  on i-th family  of W, using partially  computed fourier
# matrix  M, and  a  given permutation-with-signs  p,  tries to  fill-in
# Ennola scalars e. (start with some known scalars in e).
# Return false if some incompatibility discovered
scalpred:=function(xi,i,W,p,e,M)local uc,f,n,e,j,i,scal;
  uc:=UnipotentCharacters(W);
  f:=uc.families[i];
  n:=Length(f.charNumbers);
  e:=ShallowCopy(e);
  for j in [1..n] do
    i:=PositionProperty([1..n],i->not M.Valueat(i,j) in [false,0] and
       M.Valueat(AbsInt(p[i]),j)<>false);
    if i<>false then
      scal:=ComplexConjugate(M.Valueat(AbsInt(p[i]),j)*SignInt(p[i])/
                             M.Valueat(i,j));
      if not IsBound(e[j]) then e[j]:=scal;
      elif e[j]<>scal then return false;
      fi;
    fi;
  od;
  return e;
end;

# Add  to M  (linear  relations for  fourier)  resulting from  xi-Ennola
# diagonal on CS on i-th family, using a given sgnperm p for xi-Ennola
relennolaxi:=function(xi,M,p)local om,res,f,n,p,i,j,rel,uc,hc,val,i1,ijknown,W;
  W:=M.W; uc:=UnipotentCharacters(W); f:=uc.families[M.fno]; 
  n:=Length(f.charNumbers); res:=[];
  for j in [1..n] do
    ijknown:=Filtered([1..n],i->not M.Valueat(i,j) in [false,0]and p[i]<>false);
    om:=PositionProperty([1..n],i->i in ijknown and M.Valueat(AbsInt(p[i]),j)<>false);
    val:=EnnolaOps.GetScal(M.ennola,f.charNumbers[j],xi);
    if Length(val)>1 then
      if om<>false then
	om:=M.Valueat(AbsInt(p[om]),j)*SignInt(p[om])/M.Valueat(om,j);
	EnnolaOps.SetScal(M.ennola,
         [ComplexConjugate(om)/OmegaChi(uc)[f.charNumbers[j]]^xi],
           f.charNumbers[j],xi);
      else Print(" cannot determine ",Format(E(OrderCenter(W))^xi),
        "-scalar for ",TeXStrip(uc.TeXCharNames[f.charNumbers[j]]),"\n");
      fi;
    else om:=ComplexConjugate(val[1]*OmegaChi(uc)[f.charNumbers[j]]^xi);
    fi;
    if om<>false then
      for i in Filtered([1..n],i->p[i]<>false) do
        if not Makerel(M,[[M.at(i,j),om],[M.at(AbsInt(p[i]),j),-SignInt(p[i])]],0)
	then return false;fi;
      od;
    else
      for i in ijknown do
        for i1 in Filtered(ijknown,k->k>i) do
	  if not Makerel(M,[[M.at(AbsInt(p[i]),j),SignInt(p[i])/M.Valueat(i,j)],
	             [M.at(AbsInt(p[i1]),j),-SignInt(p[i1])/M.Valueat(i1,j)]],0)
	  then return false;fi;
      od;od;
    fi;
  od;
  Print(Length(res)," relations ");
  return true;
end;

# sq(M) get Mvp matrix from linear relations M
sq:=function(M)local v;v:=Values(M);
  return List([1..M.n],i->List([1..M.n],j->v[M.at(i,j)]))*Mvp("a")^0;
end;

# relconj(W,i,M[,p]) add to M for i-th family of W relations from S^2=complex 
# conjugation; if given p=sgnperm for Ennola
relconj:=function(arg)local W,i,M,p,uc,f,n,i,j,k,v,m,ud;
  W:=arg[1];i:=arg[2];M:=arg[3];
  Print("conj... \c");
  uc:=UnipotentCharacters(W);f:=uc.families[i];n:=Length(f.charNumbers);
  if Length(arg)<4 then
# possibilities for the effect of complex conjugation on ith family
# (list of possible destinations for each char)
    ud:=CycPolUnipotentDegrees(W){f.charNumbers};
    p:=List(ud,p->PositionsSgn(ud,ComplexConjugate(p)));
    m:=sq(M)^2;
    for i in [1..n] do for j in [1..n] do
      if isz(m[i][j]-1) then
	for k in [1..n] do
	  if k=i then p[i]:=[j]; else p[k]:=Difference(p[k],[j]); fi;
	od;
      fi;
    od;od;
  else p:=arg[4];
  fi;
  Print("possconj=",p,"\n");
  for i in [1..n] do
    for j in Difference([1..n],List(p[i],AbsInt)) do
      v:=List([1..n],k->M.Valueat(k,j));
      if ForAll(v,x->x<>false) then
        Makerel(M,List([1..n],k->[M.at(i,k),v[k]]),0);
      fi;
    od;
  od;
end;

# add to M relations expressing M=TransposedMat(m)
relsym:=function(M)local i,j;
  for i in [1..M.n] do for j in [1..i-1] do
    Makerel(M,[[M.at(i,j),1],[M.at(j,i),-1]],0);od;od;
end;
  
# add to M relations expressing M[i]*v=c
relprod:=function(M,v,i,c)
  return Makerel(M,List([1..Length(v)],j->[M.at(i,j),v[j]]),-c);end;

# add to M relations expressing M*vp=vq where vp, vq are vectors of polynomials
relpol:=function(M,vp,vq)local i,k,coeff,max;
  coeff:=function(p,i)
    if i<p.valuation or i>Degree(p) then return 0;fi;
    return p.coefficients[i-p.valuation+1];
  end;
  max:=Maximum(List(Concatenation(vp,vq),Degree));
  for i in [1..M.n] 
    do for k in [0..max] do 
      if not relprod(M,List(vp,p->coeff(p,k)),i,coeff(vq[i],k)) then 
        return false; fi;
    od; 
  od;
  return applynew(M);
end;

# add linear relations coming from Sbar*S=Id and known lines of S
relunitary:=function(M)local m,i,j,known,v;
  m:=List([1..M.n],i->List([1..M.n],j->M.Valueat(i,j)));
  for i in [1..M.n] do 
    for j in [i..M.n] do
      if ForAll([1..M.n],k->m[j][k]<>false or 
        (m[i][k]<>false and isz(m[i][k]))) then
	v:=List(m[j],function(x) if x=false then return 0;else return
	  ComplexConjugate(x);fi;end);
	if i<>j then v:=relprod(M,v,i,0);else v:=relprod(M,v,i,1);fi;
	if not v then return false;fi;
      fi;
    od;
  od;
  applynew(M);
end;

fromunitary:=function(M)local n;
  n:=sq(M);
  return Set(Concatenation(n*ComplexConjugate(n)-IdentityMat(Length(n))));
end;

CheckMagrees:=function(M)local f,four;
  if M.error<>true then return;fi;
  f:=UnipotentCharacters(M.W).families[M.fno];
  four:=f.fourierMat;
  checkMagrees(M,Concatenation(List([1..M.n],i->four[i]{[1..i]})),
  #inverse of at: [i,j] from linear position
  function(p)local i,j;
    i:=First([1..p],i->i*(i+1)/2>=p);
    j:=p-i*(i-1)/2;
    return [i,j];
  end);
end;
  
initdetfam:=function(W,fno)local uc,q,f,M,n;
  uc:=UnipotentCharacters(W);q:=X(Cyclotomics);
  f:=uc.families[fno]; 
  n:=Length(f.charNumbers);
  M:=LinSys(n*(n+1)/2);M.n:=n;
  # linear position of element [i][j] in symmetric matrix represented linearly
  M.varnames:=Concatenation(List([1..n],i->List([1..i],j->SPrint("x",i,"_",j))));
  M.at:=function(i,j)
    if i>=j then return i*(i-1)/2+j;else return j*(j-1)/2+i;fi;end;
  M.Valueat:=function(i,j)return Value(M,M.at(i,j));end;
  M.W:=W;M.fno:=fno;M.error:=true;
  M.sqdisp:=function(M)local z;
    if Length(M.relations)<>0 then z:=TransposedMat(sq(M));
     z:=List(z,x->List(x,Format));
     Add(z,CharNames(uc){f.charNumbers});
    fi;
    if IsBound(z) and Length(z)<15 then Print(Format(TransposedMat(z)),"\n");
    else Print(Format(CharNames(uc){f.charNumbers}),"\n");
    fi;
  end;
  return M;
end;

# ennoladetfam(M [,ps]) if ps for Ennola not given uses what holds from degrees
ennoladetfam:=function(arg)local M,e,oldknown,newknown,z,i,p,v;
  M:=arg[1];newknown:=Length(M.known)+Length(M.relations);
  if not IsBound(M.ennola) then M.ennola:=EnnolaOps.init(M.W);fi;
  EnnolaOps.ScalarsFromFrob(M.ennola);
  EnnolaOps.PossFromScal(M.ennola);
  if Length(arg)>1 then e:=arg[2];fi;
  repeat
    oldknown:=newknown;z:=OrderCenter(M.W);
    v:=[1..z-1];SortParallel(List(v,x->OrderMod(x,z)),v);
    for i in Reversed(v) do
      if IsBound(e) then p:=LsPuiss(e,i);
      else p:=List(M.ennola.families[M.fno].poss[i],function(p)
              if Length(p)=1 then return p[1];else return false;fi;end);
      fi;
      Print("Ennola ",Format(E(z)^i)," diagonal on CS\c=>");
      if not relennolaxi(i,M,p) then return false;fi;
      applynew(M);
      CheckMagrees(M);
    od;
    relunitary(M);
    newknown:=applynew(M);
  until newknown=oldknown;
  CheckMagrees(M);
end;

# detfam(W,fno[,do not compare with stored data]) 
# relations, where CC=ComplexConjugate: 
# S=Transposed(S), S*ud=fd, S*fd=CC(ud) (since S^-1=CC(S))
# S*CC(O)*ud=O*ud [since sh=CC(O)SCC(O) preserves ud]
# S*O*ud=CC(O)*S*fd [since S^2CC(O)ud=CC(O)S^2ud =>SOud=CC(O)S^2ud]
detfam:=function(arg)local W,fno,M,p,uc,q,ud,fd,f,O;
  W:=arg[1];fno:=arg[2];
  M:=initdetfam(W,fno);
  M.error:=Length(arg)<3;
  uc:=UnipotentCharacters(W);q:=X(Cyclotomics);
  f:=uc.families[fno];
  O:=Eigenvalues(f);
  fd:=FakeDegrees(uc,q){f.charNumbers};
  ud:=UnipotentDegrees(W,q){f.charNumbers}; 
  Print("S*fd=udbar\c=>");relpol(M,fd,ComplexConjugate(ud));
  CheckMagrees(M);
  ennoladetfam(M);
  Print("S*ud=fd\c=>");relpol(M,ud,fd);
  Print("S*O^-1*ud=O*ud\c=>");
  relpol(M,Zip(ComplexConjugate(O),ud,function(x,y)return x*y;end),
           Zip(O,ud,function(x,y)return x*y;end));
  Print("S*O*ud=O^-1*ComplexConjugate(ud)\c=>");
  relpol(M,Zip(O,ud,function(x,y)return x*y;end),
    ComplexConjugate(Zip(O,ud,function(x,y)return x*y;end)));
  ennoladetfam(M);
# relconj(W,fno,M);applynew(M);
  M.sqdisp(M);
  return M;
end;

# the Ennola theory says that if E:=DiagonalMat(Ennola-scalars)
# and P:=matrix of permutation-with signs that Ennola does and M=fourier,
# then one has M*P=E*M
# here p is a Ls
# result is a list of triples [i,e,-e] meaning E[i][i]=e or -e
newratios:=function(M,p)local m,g,pos1,pos2,letter,val,j,res,IsMonomial,IsScal;
  IsMonomial:=x->Length(x.coeff)=1;IsScal:=x->ScalMvp(x)<>false;
  m:=sq(M);
  res:=[];
  for j in [1..Length(m)] do
    g:=List([1..Length(m)],i->[SignInt(p[i])*m[AbsInt(p[i])][j],m[i][j]]);
    if not ForAny(g,x->ForAll(x,IsScal)) then
      pos1:=List(Filtered(g,x->IsMonomial(x[2]) and IsScal(x[1])),x->x[2]/x[1]);
      pos2:=List(Filtered(g,x->IsMonomial(x[1]) and IsScal(x[2])),x->x[1]/x[2]);
      letter:=Intersection(List(pos1,x->x.elm[1].elm[1]),
                           List(pos2,x->x.elm[1].elm[1]));
      if Length(letter)>0 then
        pos1:=First(pos1,x->x.elm[1].elm[1]=letter[1]);
        pos2:=First(pos2,x->x.elm[1].elm[1]=letter[1]);
	val:=GetRoot(pos2.coeff[1]/pos1.coeff[1],2);
        Add(res,[j,val,-val]);
      fi;
    fi;
  od;
  return res;
end;

ratios:=function(M,p)local j,res,f,v;
  res:=[];
  f:=function(j)local i,v1,v2;
    for i in [1..M.n] do
      v1:=M.Valueat(AbsInt(p[i]),j);v2:=M.Valueat(i,j);
      if  v1<>false and v2<>false and not v1=0*v1 then 
        return SignInt(p[i])*v1/v2;
      fi;
    od;
    return false;
  end;
  for j in [1..M.n] do v:=f(j);if v<>false then res[j]:=v;fi; od;
  return res;
end;

tryratio:=function(W,fno,M,p,pos,ratio)local f,i,oldratios,newratios,inds;
  f:=UnipotentCharacters(W).families[fno];
  oldratios:=ratios(M,p);
  newratios:=ShallowCopy(oldratios);newratios[pos]:=ratio;
  while newratios<>oldratios do
    inds:=Filtered([1..M.n],i->IsBound(newratios[i]) and
			      not IsBound(oldratios[i]));
    for pos in inds do
      for i in [1..M.n] do
        Makerel(M,[[M.at(AbsInt(p[i]),pos),1],
	          [M.at(i,pos),-SignInt(p[i])*newratios[pos]]],0);
      od;
    od;
    applynew(M);
    oldratios:=newratios;
    newratios:=ratios(M,p);
  od;
end;

ratiosmonomialmat:=function(m,p)local A,n;
  n:=Length(m);
  A:=List([1..n],function(j)local v;
    v:=Filtered([1..n],i->not isz(m[i][j]));
    return Set(List(v,i->[m[AbsInt(p[i])][j]*SignInt(p[i]),m[i][j]]));
    end);
  A:=List(A,function(v)local ratio,c;
    ratio:=false;
    for c in v do
      if IsCyc(c[2]) or c[2].elm[1].elm=[] then ratio:=c[1]/c[2];fi;
    od;
    if ratio<>false then  v:=List(v,c->c[1]-ratio*c[2]);fi;
    return Set(v);end);
  return A;
  #return Filtered(A,x->Length(x)>1 or x[1].elm[1].elm<>[]);
end;

Makerelandconj:=function(M,v,rhs)local n;n:=M.total/2;
 if not Makerel(M,v,rhs) then return false;fi;
 if not Makerel(M,List(v,function(x)
   if x[1]>n then return [x[1]-n,ComplexConjugate(x[2])];
   else return [x[1]+n,ComplexConjugate(x[2])];fi;end),ComplexConjugate(rhs))
 then return false;fi;
 return true;
end;

detfrob2:=function(M)local oldknown,newknown,i,j,n;
  newknown:=applynew(M);n:=M.total/2;
  Print("S*T*S=Tbar*Sbar*Tbar=>\c");
  repeat
    oldknown:=newknown;
    for i in Filtered(M.known,x->x<=n) do
      for j in Filtered(M.known,x->x<=n and x<=i) do
        if not Makerelandconj(M,List([1..n],k->[k,M.S[i][k]*M.S[k][j]]),
         -ComplexConjugate(Value(M,i)*Value(M,j)*M.S[i][j]))
	 then return false;fi;
    od; od;
    for i in Filtered(M.known,x->x<=n) do
      for j in Difference([1..n],M.known) do
      if not Makerelandconj(M,Concatenation(List([1..n],k->[k,M.S[i][k]*M.S[k][j]]),
         [[j+n,-ComplexConjugate(Value(M,i))*ComplexConjugate(M.S[i][j])]]),0)
       then return false;fi;
    od; od;
    newknown:=applynew(M);
  until newknown=oldknown;
end;

# detfrob (W,fno[,fouriermat])
# (uses fourierMat for fno-th family of W if arg[3] not given)
detfrob:=function(arg)local W,fno,uc,f,S,n,M,hc,i,j,h,
       oldknown,newknown,cdeg,deg,m,ud,q;
  W:=arg[1];fno:=arg[2];
  uc:=UnipotentCharacters(W);
  cdeg:=function(deg,p)
    if deg<p.valuation or deg>Degree(p) then return 0;
    else return p.coefficients[1+deg-p.valuation];fi;
  end;
  f:=uc.families[fno];
  if Length(arg)=2 then S:=f.fourierMat;else S:=arg[3];fi;
  n:=Length(f.charNumbers);
  M:=LinSys(2*n);
  hc:=List(f.charNumbers,
        x->PositionProperty(uc.harishChandra,h->x in h.charNumbers));
  M.known:=Filtered([1..n],i->uc.harishChandra[hc[i]].relativeType<>[]);
  M.values:=List(M.known,i->uc.harishChandra[hc[i]].eigenvalue);
  M.known:=Concatenation(M.known,M.known+n);
  M.values:=Concatenation(M.values,ComplexConjugate(M.values));
  # relations coming from Shintani principal series
  # si p= indices de la serie principale alors (Sbar*T*ud){p}=ud{p}
  Print("Shintani principal series=>\n");
  q:=X(Cyclotomics);
  ud:=UnipotentDegrees(W,q){f.charNumbers};
  m:=Maximum(List(ud,Degree));
  for j in Filtered([1..n],
       i->f.charNumbers[i] in uc.harishChandra[1].charNumbers) do
    for deg in [0..m] do
      if not Makerelandconj(M,List([1..n],
       i->[i,ComplexConjugate(S[j][i])*cdeg(deg,ud[i])]),-cdeg(deg,ud[j]))
      then return false;fi;
    od;
  od;
  M.S:=S;
# detfrob2(M);
  if Length(arg)=2 then
  checkMagrees(M,Concatenation(f.eigenvalues,ComplexConjugate(f.eigenvalues)));
  fi;
  return M;
end;

# get Frobenius from result of detfrob, trying to get more values by
# using relation x*xbar=1
Frob:=function(M)local n,res,t,optimize;
  # input: 2 vectors x, xbar  of Mvp. Try to use whenever x.xbar is monomial
  # to get more values as Laurent monomials
  optimize:=function(t)local e,v,i;
    while true do
      e:=List([1..n],j->t[1][j]*t[2][j]);
      v:=List(e,ScalMvp);
      if ForAll(v,x->x=1) then return t;fi;
      if ForAny(v,x->x<>false and x<>1) then return false;fi;
      i:=PositionProperty(e,x->Length(x.elm)=1 and Length(x.elm[1].elm)<>0);
      if i=false then return t;fi;
      e:=e[i];
      v:=e.elm[1].elm[1];
      t:=List(t,r->List(r,c->Value(c,[v,Mvp(v)/e])));
    od;
  end;
  res:=Values(M); n:=Length(res)/2;
  res:=[res{[1..n]},res{n+[1..n]}]*Mvp("a")^0;
  t:=List(Set(List(TransposedMat(res),x->Product(x)-1)),quadratMvp);
  t:=Filtered(t,x->x<>false);
  t:=List(Cartesian(t),Set);
  Print("t=",t,"\n");
  t:=Filtered(t,x->Length(Set(List(x,y->y[1])))=Length(x));
  Print("t=",t,"\n");
  if Length(t)=0 then return res;
  else res:=List(t,x->List(res,a->List(a,b->Value(b,Concatenation(x)))));
  fi;
  return Filtered(List(res,optimize),x->x<>false);
end;

# ambiguity group of fno-th family of W (permutations fixing unipotent degrees)
ambig:=function(W,fno)local uc,f,q,ud,p,pp,d,res,n;
  uc:=UnipotentCharacters(W);
  f:=uc.families[fno];
  ud:=CycPolUnipotentDegrees(W);
  p:=List(f.charNumbers,x->Position(ud,ud[x]));
  pp:=List(Filtered(Collected(p),x->x[2]>1),x->x[1]);
  res:=[];
  n:=Length(f.charNumbers);
  for d in pp do
    q:=Filtered([1..n],i->p[i]=d);
    Print(q,"=",CharNames(uc){f.charNumbers{q}},"\n");
    Add(res,q);
  od;
  return Group(Concatenation(List(res,q->List([1..Length(q)-1],
     i->(q[i],q[i+1])(n+q[i],n+q[i+1])))),());
end;

# List of representatives of each class of permutations-with-signs for
# xi-Ennola of fno-th family f of W  under ambiguity group of f.
# the representative in the orbit which occurs in CHEVIE data is taken
# to be compatible with CHEVIE data.
# The third argument is Ennola(W)
# [detfam is to be called with W.ennola unbound if testing nonoccuring orbit]
repsEnnola:=function(W,fno,e)local p,n,orbits,res,e;
  p:=EnnolaBete(W,fno,1/OrderCenter(W));
  n:=Length(p);
  p:=DistinctCartesian(p);
  p:=Filtered(p,x->Length(Intersection(x,-x))=0);
  p:=Orbits(ambig(W,fno),List(p,LsToPs));
  orbits:=List(p,o->List(o,ps->PsPermuted([1..n],ps)));
  p:=List(e,x->x.ls[fno]);
  n:=Set(List(p,ps->PositionProperty(orbits,y->ps in y)));
  if Length(n)<>1 then Error("theory");fi;
  Print("All possibilities occuring in W.ennola fall in orbit ",n[1],"\n");
  res:=List(orbits,x->x[1]);
  for n in p do res[PositionProperty(orbits,y->n in y)]:=n;od;
  return res;
end;

# set the fno-th family of W to f
# [at least fourierMat and eigenvalues should have been filled in f]
tryfam:=function(W,fno,f)local uc,i,h;
  uc:=UnipotentCharacters(W);
  uc.families[fno].fourierMat:=f.fourierMat;
  uc.families[fno].eigenvalues:=f.eigenvalues;
  Unbind(uc.eigenvalues);
  for i in [1..Length(uc.families[fno].charNumbers)] do
    h:=PositionProperty(uc.harishChandra,h->uc.families[fno].charNumbers[i]
       in h.charNumbers);
    uc.harishChandra[h].eigenvalue:=f.eigenvalues[i];
  od;
  return W;
end;

# check that Ennola is multiply by Ennola(special) in fusion algebra
checkEnnolaBonnafe:=function(W)local uc,en,f,i,B,r,es;
  uc:=UnipotentCharacters(W);
  en:=Ennola(W);
  for i in [1..Length(uc.families)] do
    f:=uc.families[i];
    B:=Basis(FusionAlgebra(f));
    for r in en do
      es:=r.fls[i][f.special];
      if es<0 then es:=-B[-es]*B;else es:=B[es]*B;fi;
      es:=List(es,function(b)
        if Length(b.coefficients)<>1 then Error();fi;
	if not b.coefficients[1][1] in [-1,1] then Error();fi;
	return Product(b.coefficients[1]);
	end);
      if es<>r.fls[i] then Error();fi;
    od;
  od;
end;
