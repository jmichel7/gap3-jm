Checkgd:=function(W)local ud,uc,gd,q,i,pr,l,j;
  pr:=function(i)
    Print("[",i,"=",uc.charNames[i],"=",StringSymbol(uc.charSymbols[i][1]),"]");
  end;
  q:=X(Cyclotomics);
  ud:=UnipotentDegrees(W,q);
  uc:=UnipotentCharacters(W);
  gd:=List(uc.charSymbols,x->Value(CycPolGenericDegreeSymbol(x[1]),q));
  if ud=gd then Print("OK!\n");return;fi;
  for i in [1..Length(ud)] do
    if ud[i]<>gd[i] then 
      pr(i);Print("\n   ",
          CycPol(ud[i])," stored while\n   ",CycPol(gd[i])," computed\n");
      l:=Filtered([1..Length(ud)],j->ud[j]=gd[i]);
      if Length(l)>0 then 
        Print("   occurs for ");for j in l do pr(j);Print("\n     ");od;
      else
      l:=Filtered([1..Length(ud)],j->ud[j]=-gd[i]);
        Print("   opposite occurs for ");for j in l do pr(j);Print("\n     ");od;
      fi;
      Print("\n");
    fi;
  od;
end;

#MakeFourier(2,[0,1]);
#MakeFourier(3,[0,0,1]);
#MakeFourier(3,[0,1,2]);        ## G(3,3,3)
#MakeFourier(3,[0,0,1,1,2,2]);  ## G(3,3,3)
#MakeFourier(3,[0,0,1,1,2,3]);  ## G(3,3,4)
#MakeFourier(4,[0,0,1,1]);      ## G(4,4,2)
#MakeFourier(4,[0,0,1,2]);      ## G(4,4,3)
#MakeFourier(4,[0,0,0,1,1,1,2,2]);  ## G(4,4,3)
#MakeFourier(4,[0,1,2,3]);
#MakeFourier(5,[0,0,0,1,1]);      ## G(5,5,2)
#MakeFourier(6,[0,0,0,0,1,1]);    ## G(6,6,2)
#MakeFourier(6,[0,0,0,1,1,1]);    ## G(6,6,3)

checkcyc:=function(e)local uc,W,r;
  W:=ComplexReflectionGroup(e,1,1);
  r:=FourierImprimitive(e,Concatenation(List([1..e-1],x->0),[1,1]));
  uc:=UnipotentCharacters(W);
  return MatPermPerm(r.mat,uc.families[2].fourierMat);
end;

Principal:=function(S)local e,d,m;
  e:=Length(S);m:=QuoInt(Sum(S,Length),e);d:=Sum(S,Length)-e*m;
  return (d=0 and ForAll(S,x->Length(x)=m)) or
         (d=1 and Number(S,x->Length(x)<>m)=1);
end;

mf:=function(S)local e,m,d,q,res,r,ud,fd;
  e:=Length(S);m:=QuoInt(Sum(S,Length),e);d:=Sum(S,Length)-e*m;
  r:=FourierImprimitive(e,Collected(Concatenation(S)));
  q:=X(Cyclotomics);
  res:=rec(symbs:=r[1],mat:=r[2],frob:=r[3]);
  if d=0 then res.symbs:=List(res.symbs,x->x[1]);fi;
  ud:=List(res.symbs,x->Value(CycPolGenericDegreeSymbol(x),q));
  fd:=List(res.symbs,function(x)if Principal(x) then return
                       Value(CycPolFakeDegreeSymbol(x),q);
                else return 0*q;fi;end);
  if fd*List(res.mat,x->List(x,y->GaloisCyc(y,-1)))<>ud then Error();fi;
  if res.mat=IdentityMat(Length(res.mat)) then
    res.symbs:=res.symbs{[1]};
    res.frob:=res.frob{[1]};
    res.mat:=[[1]];
  fi;
  return res;
end;

PossSymbol:=function(s)local e;
  e:=Length(s);
 return Filtered(List([0..e-1],i->s{1+List([i..i+e-1],x->x mod e)}),
    x->DefShape(List(x,Length))=0);
end;

checkfam:=function(W)local uc,r,f,s,p,v,v2,cs,i,j,sgn;
  uc:=UnipotentCharacters(W);
  for i in [1..Length(uc.families)] do
  f:=uc.families[i];
  s:=uc.charSymbols[f.charNumbers[1]];
  r:=mf(s[1]);
  cs:=List(uc.charSymbols,x->x[1]);
  v2:=List(r.symbs,x->PositionProperty(cs,y->x in PossSymbol(y)));
  sgn:=v2*0+1;
  for j in [1..Length(v2)] do
    if r.symbs[j]<>cs[v2[j]] then
      sgn[j]:=Value(CycPolGenericDegreeSymbol(r.symbs[j])/CycPolGenericDegreeSymbol(cs[v2[j]]),1);
      r.mat[j]:=r.mat[j]*sgn[j];
      r.mat{[1..Length(v2)]}[j]:=r.mat{[1..Length(v2)]}[j]*sgn[j];
      Print(i,": should be ",StringSymbol(r.symbs[j])," instead of ",
        StringSymbol(cs[v2[j]])," rel. sgn.",sgn[j],"\n");
    fi;
  od;
  p:=List(f.charNumbers,x->Position(v2,x));
  if p=[false] then p:=[1];fi;
  s:=MatPermPerm(r.mat{p}{p},f.fourierMat);
  if s=false then 
      Print(f.charNumbers,"\n");
      showsgn(r.mat{p}{p},f.fourierMat);
      Error();fi;
  if s<>[1..Length(s)] then Print(i,":", s,"\n");fi;
  if r.frob{p}<>f.eigenvalues then
    Print("eigen computed ",r.frob{p}," instead of:",f.eigenvalues ,"\n");
  fi;
  od;
end;

# frob for d=1
frob:=function(s) local e,m;e:=Length(s);m:=QuoInt(Sum(s,Length),e);
  return E(24)^(-2*(e^2-1)*m)*E(2*e)^Sum([0..e-1],i->-(i^2+e*i)*Length(s[i+1]));
end;

hc:=function(p,r)local l,pp,cusp,ser;
ser:=S->[ReducedSymbol(List(S,x->[0..Length(x)-1])),List(S,PartBeta)];
l:=CHEVIE.R("CharSymbols","imp")(p,1,r);
pp:=List(l,ser);
cusp:=Set(List(pp,x->x[1]));
SortParallel(List(cusp,RankSymbol),cusp);
return List(cusp,function(c)local cr,res,rp;
  cr:=RankSymbol(c);
  rp:=List(CHEVIE.R("CharSymbols","imp")(p,1,r-cr)
   {[1..Length(PartitionTuples(r-cr,p))]},x->List(x,PartBeta));
  res:=rec(levi:=[1..cr]);
  if cr<r then res.parameterExponents:=[List(c,Length)];
  else res.parameterExponents:=[];
  fi;
  Append(res.parameterExponents,[2+cr..r]*0+1);
  res.relativeType:=rec(series:="ST",indices:=[1+cr..r],rank:=r-cr,p:=p,q:=1);
  res.cuspidalName:=IntListToString(List(c,Length));
  res.eigenvalue:=frob(c);
  res.charNumbers:=
                List(rp,x->Position(l,SymbolPartitionTuple(x,List(c,Length))));
  return res;end);
end;

DefShape:=function(s)local e;e:=Length(s);
 return (Binomial(e,2)*QuoInt(Sum(s),e)-s*[0..e-1]) mod e;
end;

Rotations:=function(s)local e;e:=Length(s);
  return Filtered(List([1..e-1],i->s{1+List([i..i+e-1],x->x mod e)}),
     x->DefShape(x)=0);
end;

# possible shapes for cuspidal symbols of rank<=r, length e, Inhalt=d mod e
ShapesSymbols:=function(r,e,d)local f,res,m,new;
  f:=function(lim2,sum,nb,max)local res,a;
    if nb=1 then 
      if sum=0 then return [[sum]];
      else return [];
      fi;
    fi;
    res:=[];a:=QuoInt(sum,nb-1);
    while a<=max and Binomial(a,2)<=lim2 and a<=sum do
      Append(res,List(f(lim2-Binomial(a,2),sum-a,nb-1,a),
	   x->Concatenation([a],x)));
      a:=a+1;
    od;
    return res;
  end;
  res:=[];m:=0;
  repeat new:=f(r+QuoInt((m*e+d-1)*(m*e+d-e+1),2*e),d+m*e,e,d+m*e);
	 Append(res,new); m:=m+1;
  until Length(new)=0;
  res:=Concatenation(List(res,x->Arrangements(x,e)));
  return Filtered(res,s->DefShape(s)=0 and ForAll(Rotations(s),x->not x>s));
end;

# SymbolsFor(e,e,r)=unipotent symbols for G(e,e,r)
# SymbolsFor(e,1,r)=unipotent symbols for G(e,1,r)
SymbolsFor:=function(p,q,r)
  return Concatenation(List(ShapesSymbols(r,p,q mod p),s->
    Filtered(List(PartitionTuples(r-RankSymbol(List(s,x->[0..x-1])),p),
       x->SymbolPartitionTuple(x,s)),IsReducedSymbol)));
end;

SchurElem:=function(e,q,r,phi,para)local res,m,v,cusp,S;
  cusp:=[1..e]*0;#cusp[1]:=1;
  S:=SymbolPartitionTuple(phi,cusp);
  m:=QuoInt(Sum(S,Length),e);
  q:=-para[2][1]/para[2][2];v:=Copy(para[1]);#v{[2..e]}:=q*v{[2..e]};
  res:=Product([1..e],i->Product(S[i],l->Product([1..l],k->
      Product(v,w->q^k*v[i]-w))));
  res:=res*Product([1..e],i->Product([i+1..e],j->(v[i]-v[j])^m));
  res:=res*q^Sum([1..m-1],i->Binomial(i*e+1,2));
  res:=res/Product([1..e],i->Product([i+1..e],j->Product(S[i],
   l->Product(S[j],m->q^l*v[i]-q^m*v[j]))));
  res:=res/(-1)^(Binomial(e,2)*Binomial(m,2)+r*(e-1));
  res:=res/(q-1)^r;
  res:=res/q^(Sum(S,Sum)-r);
  res:=res/Product(v)^r;
  res:=res/Product([1..e],i->Product(Filtered(Cartesian(S[i],S[i]),
       x->x[1]>x[2]),x->v[i]*(q^x[1]-q^x[2])));
  return res;
end;

sch1:=function(p,q,r)local q,para;
  q:=X(Cyclotomics);para:=[List([0..p-1],i->E(p)^i),[q,-1]];
  para[1][1]:=q;
  return List(PartitionTuples(r,p),x->SchurElem(p,1,r,x,para));
end;

sch2:=function(p,q,r)return SchurElements(Hecke(ComplexReflectionGroup(p,q,r),
   X(Cyclotomics)));
end;

sch:=function(p,q,r,phi,para)local Z,w,m,S;
    S:=SymbolPartitionTuple(phi,[1..p]*0);
    m:=QuoInt(Sum(S,Length),p);
    q:=-para[2][1]/para[2][2];Z:=List([0..p-1],i->E(p)^i);
    return Product([1..p],i->Product(S[i],l->Product([1..l],k->
	Product(Z,w->q^k*Z[i]-w))))*
	Product([1..p],i->Product([i+1..p],j->(Z[i]-Z[j])^m))*
	q^Sum([1..m-1],i->Binomial(i*p+1,2))/
	Product([1..p],i->Product([i+1..p],j->Product(S[i],
     l->Product(S[j],m->q^l*Z[i]-q^m*Z[j]))))/
     (-1)^(Binomial(p,2)*Binomial(m,2)+r*(p-1))/(q-1)^r/q^(Sum(S,Sum)-r)/
     Product(Z)^r/Product([1..p],i->Product(Filtered(Cartesian(S[i],S[i]),
	 x->x[1]>x[2]),x->Z[i]*(q^x[1]-q^x[2])));
end;

ct:=function(p,r)local W,H,s,ss,q;
  W:=ComplexReflectionGroup(p,p,r);
  q:=X(Cyclotomics);H:=Hecke(W,q);
  s:=SchurElements(H);
  ss:=List(CHEVIE.imp.CharSymbols(p,p,r),phi->sch(p,q,r,List(phi,PartBeta),H.parameter));
  return List([1..Length(s)],i->s[i]/ss[i]);
end;

classreps:=function(p,r)local res,w;
 res:=List(PartitionTuples(r,p),s->Concatenation(List([1..Length(s)],
      i->List(s[i],t->[t,i-1]))));
 for w in res do Sort(w,function(a,b)return a[1]<b[1] or a[1]=b[1] and
 a[2]>b[2];end);od;
    w:=function(S)local l,res,i,j;
      l:=0;res:=[];
      for i in [1..Length(S)] do
	for j in [1..S[i][2]] do Append(res,[l+1,l..1]);Append(res,[2..l+1]);od;
	Append(res,[l+2..l+S[i][1]]);l:=l+S[i][1];
      od;
      return res;
    end;
    return List(res,w);
end;

cch:=function(p,r,q)local H,T;
  if q=1 then return CharTable(ComplexReflectionGroup(p,1,r)).irreducibles;fi;
  H:=Hecke(ComplexReflectionGroup(p,1,r),q);
  if p=2 then
    T:=Basis(H,"T");
    return TransposedMat(List(classreps(2,r),x->HeckeCharValues(T(x))));
  else return CharTable(H).irreducibles;
  fi;
end;

cch2:=function(p,r,q)local W,H,r,cl;
  W:=ComplexReflectionGroup(p,1,r);
  H:=Hecke(W,q);
  r:=Representations(H);
  cl:=ChevieClassInfo(W).classtext;
  if IsList(q) then q:=Product(q,Product);fi;
  r:=[CharTable(H).irreducibles, 
         List(r,chi->q^0*CharRepresentationWords(chi,cl))];
  CHEVIE.Check/EqObj(r[1],r[2]);
  return r;
end;

chv:=function(W)local q,v,tv;
  q:=X(Cyclotomics);
  v:=List(UnipotentDegrees(W,q),Degree);
  tv:=List(UnipotentCharacters(W).charSymbols,x->HighestPowerGenericDegreeSymbol(x[1]));
  if v<>tv then return [v-tv];fi;
end;
