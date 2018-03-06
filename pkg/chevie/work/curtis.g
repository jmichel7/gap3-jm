# Curtis sign of HC series h
CSign:=function(W,h)
  return (-1)^(W.semisimpleRank-Sum(h.relativeType,t->t.rank));
end;

# Effect of Curtis on a generic degree
DegCur:=function(W,p)local x;
  x:=Variables(p);
  if Length(x)=0 then x:="x";else x:=x[1];fi;
  return ComplexConjugate(Value(p,[x,Mvp(x)^-1]))*Mvp(x)^W.N;
end;

# Check that Curtis of chars a are chars b
CheckDeg:=function(W,a,b)local uc,ud,nud,psgn,sgn;
  Print(a,"->",b,"\n");
  uc:=UnipotentCharacters(W); ud:=UnipotentDegrees(W,Mvp("x"));
  nud:=List(ud{a},p->DegCur(W,p));
  psgn:=List(a,x->PositionProperty(uc.harishChandra,h->x in h.charNumbers));
  psgn:=List(psgn,x->CSign(W,uc.harishChandra[x]));
  if psgn<>List(b,SignInt) then 
    CHEVIE.Check.EqLists(psgn,List(b,SignInt));
    Error("signs ");
  fi;
  if List(b,i->SignInt(i)*ud[AbsInt(i)])<>nud then Error("degrees");fi;
end;

# Which part of curtisperm can be detected by degrees
CurtisByDeg:=function(W)local ud,x,nud,res,i,p;
  x:=Mvp("x"); ud:=UnipotentDegrees(W,x); nud:=List(ud,p->DegCur(W,p));
  res:=[];
  for i in [1..Length(ud)] do
    p:=PositionsSgn(ud,nud[i]);
    if Length(p)=1 then res[i]:=p[1];fi;
  od;
  return res;
end;

# alvis-curtis duality on unipotent characters
# in each HC series check the symmetry of parameters and flip them
CurtisPerm:=function(W)
  local uc,res,h,a,hh,G,h1,x,p1,c,cc,gO,cG,pp,ud,cb,sgn,L,uL;
  uc:=UnipotentCharacters(W);x:=Mvp("x");
  res:=[];
  ud:=UnipotentDegrees(W,x);
  hh:=Copy(uc.harishChandra);
  for h in hh do
    if Length(h.relativeType)>0 then
      G:=ReflectionGroup(h.relativeType[1]);
      gO:=G.OrdersGeneratingReflections;
    fi;
    h.parameterExponents:=List([1..Length(h.parameterExponents)],
      function(i)local p;
        p:=h.parameterExponents[i];
        if IsList(p) then 
	  return List([1..Length(p)],j->x^p[j]*E(Length(p))^(j-1));
	else return Concatenation([x^p],List([1..gO[i]-1],j->E(gO[i])^j)*x^0);
	fi;end);
  od;
  cb:=CurtisByDeg(W);
  for h in hh do
    if Length(h.relativeType)>1 then Error();fi;
    if Length(h.relativeType)=0 then 
      if IsBound(cb[h.charNumbers[1]]) then
         res{h.charNumbers}:=cb{h.charNumbers};
	Print("cusp\n");
         CheckDeg(W,h.charNumbers,res{h.charNumbers});
      fi;
    else 
      G:=ReflectionGroup(h.relativeType[1]);
      sgn:=CSign(W,h);
      if Set(G.OrdersGeneratingReflections)=[2] then
        cG:=CharTable(G).irreducibles;
	a:=PermListList(cG,List(cG,ch->ComplexConjugate(
	 Zip(ch,cG[PositionDet(G)],function(x,y)return x*y;end))));
	a:=Permuted(h.charNumbers,a);
	res{h.charNumbers}:=a*sgn;
	Print("rel=",ReflectionName(G),"\n");
        CheckDeg(W,h.charNumbers,res{h.charNumbers});
      else L:=ReflectionSubgroup(W,h.levi);
        c:=FindCuspidalInLevi(h.cuspidalName,L);
	uL:=UnipotentCharacters(L);
        p:=h.parameterExponents;
	p1:=List(p,l->List(l,m->x^Maximum(List(l,Degree))*m^-1));
	if not IsBound(uL.curtis) then 
	  if Number(hh,H->H.levi=h.levi)=1 then h1:=Position(hh,h);
	  else Error("no curtis for ",ReflectionName(L),"\n");fi;
	else
	  c:=TeXStrip(uL.TeXCharNames[AbsInt(uL.curtis[c][1])]);
	  h1:=PositionProperty(hh,h->h.cuspidalName=c);
	  if h1=false then Error();fi;
	fi;
	if hh[h1]=h then #inner
	  c:=List([1..Length(p)],i->PositionProperty([0..Length(p[i])-1],j->
	    Collected(p[i]*E(Length(p[i]))^j)=Collected(p1[i])));
	  if false in c then Error(); else
	    c:=List([1..Length(p)],i->E(Length(p[i]))^(c[i]-1));fi;
	  cG:=CharTable(G).irreducibles;
	  cc:=PositionProperty(cG,ch->ch[1]=1 and ForAll([1..Length(c)],
	    j->ch[PositionClass(G,G.(j))]=c[j]));
	  pp:=PermListList(cG,ComplexConjugate(List(cG,
	    ch->Zip(ch,cG[cc],function(x,y)return x*y;end))));
	  res{h.charNumbers}:=Permuted(h.charNumbers,pp)*sgn;
	Print("rel=",ReflectionName(G),"\n");
	  CheckDeg(W,h.charNumbers,res{h.charNumbers});
	else 
	   Print(Position(hh,h),"->",h1,"\n");
	   h1:=hh[h1];p:=h1.parameterExponents;
	   c:=List([1..Length(p)],i->PositionProperty([0..Length(p[i])-1],j->
	    Collected(p[i]*E(Length(p[i]))^j)=Collected(p1[i])));
	   if false in c then Error(); else
	    c:=List([1..Length(p)],i->E(Length(p[i]))^(c[i]-1));fi;
	   Print("c=",c,"\n");
	   cG:=CharTable(G).irreducibles;
	   cc:=PositionProperty(cG,ch->ch[1]=1 and ForAll([1..Length(c)],
	     j->ch[PositionClass(G,G.(j))]=c[j]));
	   pp:=PermListList(cG,ComplexConjugate(List(cG,
	     ch->Zip(ch,cG[cc],function(x,y)return x*y;end))));
	   res{h.charNumbers}:=Permuted(h1.charNumbers,pp)*sgn;
	Print("rel=",ReflectionName(G),"\n");
	   CheckDeg(W,h.charNumbers,res{h.charNumbers});
	fi;
      fi;
    fi;
  od;
  return res;
end;

# Check duality of exponents in principal d-series agrees with given c
CurtisSeries:=function(W,d,c)local s,l,l1,c,G,cG,pp,i,e,ud,nud,uc,sgn,x;
  s:=PrincipalSeries(W,d);x:=Mvp("x");
  l:=getHecke(s);l:=Value(l,["q",x]);e:=Length(l);
  l1:=List(l,m->m^-1*x^Maximum(List(l,Degree)));
  pp:=PositionProperty([1..e],j->Collected(l*E(e)^j)=Collected(l1));
  pp:=PermListList(l1,l*E(e)^pp);
  pp:=Permuted(s.charNumbers,pp);
  ud:=UnipotentDegrees(W,x); nud:=List(ud,p->DegCur(W,p));
  uc:=UnipotentCharacters(W);
  for i in [1..Length(s.charNumbers)] do
    if not ud[s.charNumbers[i]] in [nud[pp[i]],-nud[pp[i]]] then
      Error("grave");fi;
    if IsBound(c[s.charNumbers[i]]) then
      if AbsInt(c[s.charNumbers[i]])<>pp[i] then 
        Error("found ",s.charNumbers[i],"=>",pp[i]," instead of ",
	  s.charNumbers[i],"=>",c[s.charNumbers[i]],"\n");fi;
    else 
      sgn:=nud[s.charNumbers[i]]/ud[pp[i]];
      sgn:=ScalMvp(sgn);
      c[s.charNumbers[i]]:=pp[i]*sgn;
    fi;
  od;
  CheckDeg(W,s.charNumbers,c{s.charNumbers});
  return c;
end;

# Find some possible (probably wrong) Curtis from action on degrees
CurtisBete:=function(W)local ud,x;
  x:=Mvp("x"); ud:=UnipotentDegrees(W,x);
  return SignedPermListList(ud,List(ud,p->DegCur(W,p)));
end;

# Is Steinberg alternate sum of Ind^W_{W_I} (no!)
CurtisByInd:=function(W)
  return Sum(Combinations(W.generatingReflections),function(I)local L;
    L:=ReflectionSubgroup(W,I);
    return (-1)^Length(I)*TransposedMat(InductionTable(L,W).scalar)
      [PositionId(L)];
    end);
end; 

# complexconjugation on ith family of unipotent characters
ConjUnip:=function(W,i) local uc,res,ud,eig,p,i,f;
  uc:=UnipotentCharacters(W);
  f:=uc.families[i].charNumbers;
  res:=[];
  ud:=UnipotentDegrees(W,Mvp("x")){f};
  eig:=Eigenvalues(uc){f};
  p:=TransposedMat([ud,eig]);
  res:=List([1..Length(ud)],i->PositionsProperty(p,
    c->ComplexConjugate(c[1]) in [ud[i],-ud[i]] and
       ComplexConjugate(c[2])=eig[i]));
  for i in [1..Length(res)] do
    if Length(res[i])<>1 then Unbind(res[i]);
    elif ComplexConjugate(ud[i])=-ud[res[i][1]] then res[i]:=-res[i][1];
    else res[i]:=res[i][1];
    fi;
  od;
  return res;
end;

ConjUnip2:=function(W) local ud;
  ud:=UnipotentDegrees(W,Mvp("x"));
  return List([1..Length(ud)],function(i)
    if ComplexConjugate(ud[i]) in ud then 
      return Position(ud,ComplexConjugate(ud[i]));
    else
      return -Position(ud,-ComplexConjugate(ud[i]));
    fi;
    end);
end;

CheckConjUnip:=function(W) local uc,a,G,ct,b,f,s2,i,h;
  uc:=UnipotentCharacters(W);
  for h in uc.harishChandra do
    a:=h.charNumbers;
    if Length(h.relativeType)>0 then
    G:=ReflectionGroup(h.relativeType[1]);
    ct:=CharTable(G).irreducibles;
    b:=Permuted(a,PermListList(ComplexConjugate(ct),ct));
    for f in uc.families do
      s2:=MatToLs(Fourier(f)^2);
      for i in Intersection(f.charNumbers,a) do
	if b[Position(a,i)]<>f.charNumbers[AbsInt(s2[Position(f.charNumbers,i)])]
	then Error(i,h);
	fi;
      od;
    od;
    fi;
  od;
end;
