# compares lists or matrices up to sign
showsgn:=function(a,b)
  if IsList(a[1]) then return 
    Concatenation(Zip(a,b,function(a,b)return SPrint(showsgn(a,b),"\n");end));
  fi;
  return Zip(a,b,function(x,y)
    if x=y then if x=0 then return '0'; else return '+'; fi;
    elif x=-y then return '-'; else return '*'; fi;end);
end;

# returns list of line discrepancies up to sign
cmpsgn:=function(m,m1)return Zip(m,m1,function(l,l1)
  return Number(TransposedMat(l,l1), x->x[1]<>x[2] and x[1]<>-x[2]);end);
end;

# try to make family i-th family of W fit family record g 
fitfam:=function(W,i,g)local n,v,p,f,r,uc;
  uc:=UnipotentCharacters(W);
  f:=uc.families[i];
  p:=PsMatMat(g.fourierMat,f.fourierMat,g.eigenvalues,f.eigenvalues);
  p:=Elements(PsMatStab(g.fourierMat,g.eigenvalues))*p;
  v:=List(p,x->LsToPermAndSigns(PsToLs(x,Length(g.fourierMat))));
  SortParallel(List(v,x->Number(x[2],y->y<0)),v);
  Print("perm=",v[1][1],"\n");
  Print(FormatTable(List(v,r->Permuted(r[2],r[1]^-1)),
        rec(rowLabels:=List(v,r->Number(r[2],x->x<0)),
            columnLabels:=CharNames(uc){f.charNumbers})));
  if IsBound(g.name) then n:=g.name;else n:="Z3";fi;
  for r in v do
    Cut(SPrint("Family(\"",n,"\",",FormatGAP(Permuted(f.charNumbers,r[1])),
	  ",\n\trec(signs:=",FormatGAP(r[2]),"))"));
  od;
  return v;
end;

# EnnolaTest(W,a) show Ennola of families of W isomorphic to
# family a in the CS basis.
EnnolaTest:=function(W,a)
  local e,f,p,g,m,W,i,j,Z,A,find,uc,sgn,tbl,l,VecChi,ind,En,ud,fd,special;
# vector of (x,\psi)->\chi(x)\omega_\psi(z)(x,psi) où l'élément central z 
# est donné par son numéro de classe et le caractère central chi par son 
# numéro de caractère.
# One of these should describe Ennola in the basis of Character sheaves
  VecChi:=function(f,z,chi)local g,basis,z,fchi; g:=f.group;
    fchi:=f.operations.chi();
    basis:=List([1..Length(fchi)],i->rec(x:=f.operations.x()[i],chi:=fchi[i]));
    z:=ConjugacyClasses(g)[z].representative;
    return List(basis,function(b)local cent; cent:=Centralizer(g,b.x);
      return fchi[chi][PositionClass(g,b.x)]*b.chi[PositionClass(cent,z)]/
	b.chi[PositionClass(cent,())];end);
  end;
  uc:=UnipotentCharacters(W);
  ud:=CycPolUnipotentDegrees(W);
  fd:=uc.operations.FakeDegrees(uc,X(Rationals));
  En:=Ennola(W);
  if Length(En)>1 then Error("ambiguity in Ennola");fi;
  En:=En[1];
  ind:=Filtered([1..Length(uc.families)],i->Length(uc.families[i].charNumbers)=
    Length(a.fourierMat));
  for i in ind do
    f:=uc.families[i];
    sgn:=(-1)^uc.A[f.charNumbers[1]];
    special:=PositionProperty(f.charNumbers,j->ud[j].valuation=fd[j].valuation);
    Print("\nFamily ",i," sgn=",sgn," a.special=",a.special,"\n");
    if special<>f.special then 
      Print("*** special=",special," but f.special=",f.special,"\n");fi;
    e:=En.ls[i];
    m:=f.fourierMat;
    if IsBound(f.signs) then m:=m^DiagonalMat(f.signs);fi;
    p:=PermMatMat(a.fourierMat,m,a.eigenvalues,f.eigenvalues);
    Print("perm from M(Gamma) to F:",p,"\n");
    if p=false then
      if IsBound(f.name) then e:=f.name;else e:="??";fi;
      p:=PsMatMat(a.fourierMat,m,a.eigenvalues,f.eigenvalues);
      if p=false then Print("!!!!!!!!!!! no match for family ",e,"\n");
      else Print("!!!!!!! family ",e," could be matched by ",
        PsToLs(p,Length(m)),"\n");
      fi;
    else
    e:=SignedPermutationMat(e);
    if IsBound(f.signs) then e:=e^DiagonalMat(f.signs);fi;
    e:=OnMatrices(e,p^-1);special:=special/p;
    if IsCoxeterGroup(W) then g:=MatStab(a.fourierMat,a.eigenvalues);
    else g:=PsMatStab(a.fourierMat,a.eigenvalues);
    fi;
    e:=List(Elements(g),x->rec(mat:=e^SignedPermutationMat(PsToLs(x,Length(e))),
      perm:=x,special:=1+((special^x-1)mod Length(e))));
    for p in e do p.mat:=p.mat^(a.fourierMat^-1);od;
    if not ForAny(e,x->IsDiagonalMat(x.mat))then 
      Error("Ennola should be diagonal");fi;
    if not ForAll(e,x->IsDiagonalMat(x.mat))then 
      Z:=Filtered(e,x->not IsDiagonalMat(x.mat));
      e:=Filtered(e,x->IsDiagonalMat(x.mat));
      Z:=List(Z,x->x.perm);
      Print("!!!! bad permutations (giving non diagonal Ennola):\n");
      Print(Z,"\n");
    fi;
    for p in e do p.mat:=DiagonalOfMat(p.mat);od;
    Sort(e);
    g:=a.group;
    Z:=List(Elements(Centre(g)),x->PositionClass(g,x));
    A:=Filtered([1..Length(ConjugacyClasses(g))],
       i->a.operations.chi()[i][PositionClass(g,())]=1);
    find:=function(v)local i,j;
      for i in Z do for j in A do
	if v=sgn*VecChi(a,i,j) then return rec(z:=i,chi:=j,sgn:="");
        elif v=-sgn*VecChi(a,i,j) then return rec(z:=i,chi:=j,sgn:="-");
        fi;
      od;od;
      Error("not found");
      return false;
    end;
    l:=[];
    for j in [1..Length(e)] do p:=find(e[j].mat);
    if false then
      if j=1 or l[Length(l)].mat<>p then 
        Add(l,rec(perm:=e[j].perm,mat:=p));
      else Add(l[Length(l)].perm,e[j].perm);
      fi;
    else e[j].mat:=p; Add(l,e[j]);
    fi;
    od;
    Print(Join(CharNames(uc){f.charNumbers}," "),"\n");
    Print(FormatTable(List(l,m->[Representative(ConjugacyClasses(g)[m.mat.z]),
        a.operations.chi()[m.mat.chi],m.mat.sgn,m.special]),
       rec(columnLabels:=["z","chi","sgn","special"],
	      rowLabels:=List(l,m->OrderPerm(m.perm)))));
    fi;
  od;
end;

# display one of Gunter's records
displayfam:=function(f)local mat;
  mat:=TransposedMat(Concatenation([
  List(f.unpdeg,x->String(CycPol(x))),
  List(f.fakdeg,x->String(CycPol(x))),
  List(f.eigenvalues,Format)],List(f.fourierMat,x->List(x,Format))));
  Print(FormatTable(mat,
   rec(rowLabels:=List([1..Length(f.fourierMat)],String),
       columnLabels:=Concatenation(["UnipDeg","Fake Deg","Eigen"],
       List([1..Length(f.fourierMat)],x->"")))));
end;

transf:=function(r)local i,res,f,q;
  res:=[];q:=X(Cyclotomics);
  for i in [1..Length(r.frob)] do
    res[i]:=rec(fourierMat:=r.foumat[i],eigenvalues:=r.frob[i],
      unpdeg:=r.unpdeg[i],fakdeg:=r.fakdeg[i]);
  od;
  for i in [1..Length(res)] do
    for f in [Length(res[i].fakdeg)+1..Length(res[i].unpdeg)] do
      res[i].fakdeg[f]:=res[i].fakdeg[1]*0;
    od;
    res[i].fakdeg:=res[i].fakdeg*q^0;
    res[i].unpdeg:=res[i].unpdeg*q^0;
  od;
  return res;
end;

cmpfam:=function(arg) # cmpfam(transf(Gunter's record), W [, family numbers])
  local i,j,gun,gfam,myfam,W,indices,perm,uc,n,p,q,ud,fd,a,b,p1;
  gun:=arg[1];W:=arg[2];
  if Length(arg)>=3 then 
    if IsList(arg[3]) then indices:=arg[3];else indices:=[arg[3]];fi;
  else indices:=[1..Length(gun)];
  fi;
  q:=X(Cyclotomics);q.name:="q";
  uc:=UnipotentCharacters(W);
  ud:=UnipotentDegrees(W,q);
  fd:=FakeDegrees(uc,q);
  for i in indices do
    gfam:=gun[i];
    myfam:=uc.families[i];n:=myfam.charNumbers;
    myfam.unpdeg:=ud{n};myfam.fakdeg:=fd{n};
    if ForAll(["fakdeg","unpdeg","fourierMat","eigenvalues"],
      f->gfam.(f)=myfam.(f)) then Print("family #",i," agrees\n");
    else 
      p:=function(f,i)return [f.fakdeg[i],Set([f.unpdeg[i],-f.unpdeg[i]])];end;
      p:=PermListList(List([1..Length(myfam.eigenvalues)],i->p(myfam,i)),
                      List([1..Length(gfam.eigenvalues)],i->p(gfam,i)));
      if p<>false then
        p1:=function(f,i)return [f.eigenvalues[i],f.fakdeg[i],
	   Set([f.unpdeg[i],-f.unpdeg[i]])];end;
        p1:=PermListList(List([1..Length(myfam.eigenvalues)],i->p1(myfam,i)),
                         List([1..Length(gfam.eigenvalues)],i->p1(gfam,i)));
        if p1<>false then p:=p1;fi;
        Print("family #",i,": p=",p,"\n");
        myfam:=myfam^(p^-1);
    #   if gfam.fourierMat<>myfam.fourierMat then
    #     myfam.fourierMat:=per(myfam.fourierMat,cv(myfam.fourierMat,gfam.fourierMat,cmpeq));
    #   fi;
	fi;
      fi;
      if gfam.fourierMat<>myfam.fourierMat then
         if ComplexConjugate(myfam.fourierMat)=gfam.fourierMat then
           Print("## Gunter's fourierMat is conjugate to mine\n");
         elif Set(cmpsgn(gfam.fourierMat,myfam.fourierMat))<>[0] then
         Print("## gunter's family #",i,"\n");displayfam(gfam);
         Print("##       my family #",i,"\n");displayfam(myfam);
         fi;
	 Print(showsgn(gfam.fourierMat,myfam.fourierMat));
      fi;
      if gfam.eigenvalues<>myfam.eigenvalues then
         for j in [1..Length(gfam.eigenvalues)] do
	   if gfam.eigenvalues[j]<>myfam.eigenvalues[j] then
       Print("## gunter's frob#",i,".",j,":",Format(gfam.eigenvalues[j]),"\n");
       Print("##       my frob#",i,".",j,":",Format(myfam.eigenvalues[j]),"\n");
	   fi;
	 od;
      fi;
      if gfam.fakdeg<>myfam.fakdeg then
         for j in [1..Length(gfam.fakdeg)] do
	   if gfam.fakdeg[j]<>myfam.fakdeg[j] then
         Print("## gunter's fakdeg#",i,".",j,":",Format(gfam.fakdeg[j]),"\n");
         Print("##       my fakdeg#",i,".",j,"=",uc.families[i].charNumbers[j],
                         ":",Format(myfam.fakdeg[j]),"\n");
	   fi;
	 od;
      fi;
      if gfam.unpdeg<>myfam.unpdeg then
         for j in [1..Length(gfam.unpdeg)] do
	   if gfam.unpdeg[j]<>myfam.unpdeg[j] then
	     b:=[];
	     for a in [1..Length(myfam.unpdeg)] do
	       if myfam.unpdeg[a]=gfam.unpdeg[j] then Add(b,a);
	       elif myfam.unpdeg[a]=-gfam.unpdeg[j] then Add(b,-a);
	       fi;
	     od;
	     if Length(b)<>0 then
               Print("## gunter's unpdeg#",i,".",j," found at",b,
                 "=",myfam.charNumbers{List(b,AbsInt)},"\n");
	     else
         Print("## gunter's unpdeg#",i,".",j,":",CycPol(gfam.unpdeg[j]),"\n");
         Print("##       my unpdeg#",i,".",j,":",CycPol(myfam.unpdeg[j]),"\n");
	     fi;
	   fi;
	 od;
      fi;   
  od;
end;
