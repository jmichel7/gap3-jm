###########################################################################
##
#A  hecbloc.g               The CHEVIE package                 Meinolf Geck
##
#Y  Copyright (C)  03/1999   Equipe des groupes finis, Universite Paris VII
##
##  This file contains  functions  for  computing  blocks  and  defects  of
##  characters of Iwahori--Hecke algebras  specialized  at  roots  of unity 
##  over the rational nunmbers. 
##
##  We assume that  W  is any finite Coxeter group and H is a corresponding
##  generic Iwahori--Hecke algebra  where  all parameters are even positive
##  powers of an indeterminate v  (over  the complex numbers).  Throughout, 
##  we only consider specializations where  v  is mapped to  the  primitive 
##  2e-th root E(2*e), i.e., the parameters of H are e-th roots of unity. 
##
##  We also include functions that try to compute Brauer trees from partial
##  information  provided  by a certain 'restricted' character  table,  and
##  functions for computing the generic degrees of H (or, rather, the Schur 
##  elements in the case of unequal parameters).
##
##  The functions in this file also  work  for  W  of  non-crystallographic 
##  type. At some places, we use the  observation  that  X_1^{tr} X_v  is a 
##  matrix of  polynomials with  rational coefficients,  where X_1  is  the
##  character table of W and X_v is the character table of H.
##
##  The main functions are:
##
##  'AlgoSchurElements' . . . . . .  compute Schur elements/generic degrees
##  'HeckeNonSemisimple'  . . . . . find all non-semisimple specializations
##  'HeckeCentralTable' . . . . . . . . . . the table of central characters
##  'HeckeBlocks' . . . . . . . . . . . . . . . . . .  blocks of characters
##  'HeckeModularRank'  . . . . . . . . . . .  the number of simple modules
##  'HeckeDefect1Blocks'  . . . . . . . . . .  computing blocks of defect 1
##
v:=X(Cyclotomics); v.name:="v";
###########################################################################
##
#F  LinEqGenDegs( <H> ) . . . . . . . . . . . .  system of linear equations 
##  . . . . . . . . .  for the generic degrees of an Iwahori--Hecke algebra
##
##  'LinEqGenDegs' return the matrix  of  a  system  of linear equations of 
##  which the generic degrees for  the  Iwahori--Hecke  algebra  <H>  are a 
##  solution. The equations are determined using (1)  induction from proper
##  parabolic subalgebras  (and  the  knowledge  of the generic degrees for 
##  them),  (2) the  duality  operation and  (3) character values on 1-good 
##  elements. It is assumed  that  all  parameters of <H> are powers of one 
##  variable v^2.
## 
LinEqGenDegs:=function(H)
  local W,P,J,ct,system,lcm,r,r1,mon,i,j,iw0,row,sh1,ind,sgn,neu,x,z;
  W:=CoxeterGroup(H);
  P:=PoincarePolynomial(H);
  ct:=CharTable(W);
  r:=CoxeterConjugacyClasses(W);
  system:=[];

  # 1. Step: Add all equations coming from maximal parabolic subgroups
  InfoChevie("#I parabolic subgroups : \c");
  for i in [1..W.semisimpleRank] do
    J:=Difference([1..W.semisimpleRank],[i]);
    ind:=TransposedMat(InductionTable(ReflectionSubgroup(W,J),W).scalar);
    sh1:=SchurElements(HeckeSubAlgebra(H,J));
    lcm:=Lcm(Concatenation(sh1,[P]));
    for j in [1..Length(ind)] do
      if lcm=P then 
        row:=Copy(ind[j]);
      else
        row:=(lcm/P)*ind[j];
      fi;
      Add(row,lcm/sh1[j]);
      if not row in system then 
        Add(system,row);
        InfoChevie(".\c");
      fi;
    od;
  od;
  InfoChevie(" ",Length(system)," equations \n");
  
  # 2. Step: dual characters
  InfoChevie("#I Dual characters     : \c");
  sgn:=ct.irreducibles[PositionSgn(W)];
  neu:=[];
  for i in [1..Length(ct.irreducibles)] do
    x:=List([1..Length(ct.irreducibles)],j->sgn[j]*ct.irreducibles[i][j]);
    j:=Position(ct.irreducibles,x);
    if not Set([i,j]) in neu then 
      Add(neu,Set([i,j]));
    fi;
  od;
  iw0:=H.parameter[1][1]^0;
  for i in LongestCoxeterWord(W) do
    iw0:=iw0*H.parameter[i][1];
  od;
  z:=HeckeCentralMonomials(H);
  for i in neu do 
    if Length(i)=2 then
      row:=0*[1..Length(ct.irreducibles)];
      row[i[1]]:=1;
      row[i[2]]:=-iw0*z[i[1]]^(-1);
      Add(row,0*H.parameter[1][1]);
      if not row in system then 
        Add(system,row);
      fi;
    fi;
  od;
  InfoChevie(Length(system)," equations \n");
  
  # 3. Step: Add all equations coming from values on regular elements
  InfoChevie("#I 1-good elements     : \n");
  r1:=Filtered([1..Length(r)],x->Set(r[x])=[1..W.semisimpleRank] 
                                  and Length(GoodCoxeterWord(W,r[x]))=1);
  mon:=List(r1,x->HeckeCharValuesGood(H,r[x]));
  for j in [1..Length(r1)] do
    row:=[];
    for i in [1..Length(ct.irreducibles)] do
      if ct.irreducibles[i][r1[j]]=0 then
        row[i]:=0;
      else
        row[i]:=ct.irreducibles[i][r1[j]]
                    *X(Cyclotomics)^(mon[j][i].valuation/ct.orders[r1[j]]);
      fi;
    od;
    Add(row,H.parameter[1][1]*0);
    if not row in system then 
      Add(system,row);
    fi;
    Add(system,row);
  od;
  InfoChevie("#I ",Length(system)," equations \n");

  return system;
end;

###########################################################################
##
#F  AlgoSchurElements( <H> ) . . . . . . . . . . . . compute Schur elements  
##
##  'AlgoSchurElements' applies 'LinEqGenDegs' to get a  system  of  linear 
##  equations for the generic degrees of <H>.  That  system  is  solved for
##  various specializations of the parameters of <H> and then the solutions
##  are interpolated  to get the Schur  elements  of <H>. If all parameters
##  of <H> then some  simplifications  can be made  and the generic degrees 
##  are returned.
##
AlgoSchurElements:=function(H)
  local P,Pv,v,vs,max,i,j,iw0,gd,res,n,pols,gmat,mat;
  P:=PoincarePolynomial(H);
  iw0:=H.parameter[1][1]^0;
  for i in LongestCoxeterWord(CoxeterGroup(H)) do
    iw0:=iw0*H.parameter[i][1];
  od;
  if Length(Set(List(H.parameter,x->x[1])))=1 then 
    Print("#I All parameters equal : computing generic degrees\n");
    max:=Valuation(iw0);
  else
    Print("#I Unequal parameters : computing Schur elements\n");
    max:=2*Valuation(iw0);
  fi;
  gmat:=LinEqGenDegs(H);
  res:=[];
  vs:=[];
  v:=0;
  InfoChevie("#I Specializing v -> 1 to ",max+1," \c");
  while Length(res)<=max do
    v:=v+1;
    if IsInt(v/10) then 
      InfoChevie(".\c");
    fi;
    mat:=[];
    for i in [1..Length(gmat)] do
      mat[i]:=[];
      for j in [1..Length(gmat[1])] do
        if IsCyc(gmat[i][j]) then
          mat[i][j]:=gmat[i][j];
        else
          mat[i][j]:=Value(gmat[i][j],v);
        fi;
      od;
    od;
    TriangulizeMat(mat);
    mat:=Filtered(mat,x->x<>0*x);
    Pv:=Value(P,v);
    if Length(mat)=Length(mat[1])-1 then 
      gd:=List(mat,i->i[Length(mat[1])]);
      if ForAll(gd,x->x<>0) then 
        Add(vs,v);
        if Length(Set(H.parameter))=1 then 
          Add(res,gd);
        else
          Add(res,List(gd,x->v^iw0.valuation*Pv/x));
        fi;
      fi;
    fi;
  od;
  res:=TransposedMat(res);
  pols:=[];
  InfoChevie("\n#I Interpolating    \c");
  for i in [1..Length(res)] do
    Print(".\c");
    if Length(Set(H.parameter))=1 then 
      Add(pols,InterpolatedPolynomial(Cyclotomics,vs,res[i]));
    else
      Add(pols,InterpolatedPolynomial(Cyclotomics,vs,res[i])/iw0);
    fi;
  od;
  return pols;
end;

###########################################################################
##
#F  DefectPolynomial( <f>, <phi> ) . . . . . the phi-defect of a polynomial
## 
##  'DefectPolynomial'  returns  the  largest $m\geq 0$ such that $<phi>^m$
##  divides <f>.
##
DefectPolynomial:=function(f,phi)
  local m,qr;
  if Degree(f)=0 then return 0;fi;
  m:=0;
  qr:=QuotientRemainder(f,phi);
  if qr[2]<>0*qr[2] then return 0; fi;
  while qr[2]=0*qr[2] do
    m:=m+1;
    qr:=QuotientRemainder(qr[1],phi);
  od;
  return m;
end;

###########################################################################
##
#F  HeckeNonSemisimple( <H> ) . . . . . . . non-semisimple specializations
##
##  'HeckeNonSemisimple' returns the list of all $e$ such that <H> becomes 
##  non-semisimple under the specialization v -> E(2*e)  (primitive  2e-th  
##  root of unity over the rationals). 
##
HeckeNonSemisimple:=function(H)
  local lcm,m,she,d,i,phid,ind,sh,ds,e;
  ind:=H.parameter[1][1]^0;
  for i in LongestCoxeterWord(CoxeterGroup(H)) do
    ind:=ind*H.parameter[i][1];
  od;
  sh:=ind*SchurElements(H);
  lcm:=Lcm(sh);
  ds:=[];
  e:=1;
  while Length(lcm.coefficients)>1 do 
    if Value(lcm,E(2*e))=0 then
      phid:=Value(CyclotomicPolynomial(Cyclotomics,e),X(Cyclotomics)^2);
      m:=DefectPolynomial(lcm,phid);
      Add(ds,[e,m]);
      lcm:=Quotient(lcm,phid^m);
    fi;
    e:=e+1;
  od;
  return ds;
end;

###########################################################################
##
#F  'HeckeCharDefect( <H>, <d> ) . . . . . . . .  defects of the characters
##
##  'HeckeCharDefect' returns the list of <d>-defects of the Schur elements 
##  of  <H>,  i.e., for each  irreducible  character  $\chi$, the order  of 
##  E(2*d) as a zero of  the Schur element of $\chi$.
##
HeckeCharDefect:=function(H,d)
  local ind,f,i,x,phid;
  ind:=H.parameter[1][1]^0;
  for i in LongestCoxeterWord(CoxeterGroup(H)) do
    ind:=ind*H.parameter[i][1];
  od;
  phid:=[];
  for x in SchurElements(H) do 
    f:=ind*x; i:=0;
    while Value(f,E(2*d))=0 do 
      i:=i+1;
      f:=Quotient(f,X(Cyclotomics)-E(2*d));
    od;
    Add(phid,i);
  od;
  return phid;
end;

###########################################################################
##
#F  'HeckeCentralTable( <H> ) . . . . . . . the table of central characters
##
##  'HeckeCentralTable'  returns  the  matrix  $(\omega_\chi(z_C))$,  where 
##  $\chi$ runs over the irreducible characters  of <H> and $z_C$ runs over
##  the basis  of the  center of <H>, as constructed by Geck-Rouquier. This 
##  is done by inverting  the  character  table  of  <H>. The  computations
##  are performed using finite fields, interpolation and chinese remainder.
##  
HeckeCentralTable:=function(H)
  local ListPolFFE,Cetv,ApplyChinRem,d,ct,ctw,R,s,lcm,ind,tbl,tblv,sh,
                                            i,j,k,p2,intp,pr,x,res,fertig;
  if IsBound(H.centralTable) then
    return H.centralTable;
  fi;
  # Specialize v to integers mod pr and return the central tables.
  Cetv:=function(R,lcm,pr,d)
    local valFFE,v,vs,lcmv,tblv,z;
    valFFE:=function(pr,f,x)
      local val,i,id;
      id:=x^0; val:=0*id; i:=Length(f.coefficients);
      while 0<i do val:=val*x+id*(f.coefficients[i] mod pr); i:=i-1; od;
      if 0<>f.valuation then val:=val*x^f.valuation; fi;
      return val;
    end;
    if pr=1 then z:=1; else z:=Z(pr)^0; fi;
    vs:=[]; v:=1; tblv:=[];
    while Length(Set(vs))<d do
      if pr>1 and v>=pr then Print("  not enough values \n");return false;fi;
      if IsInt(v/10) then InfoChevie("+\c"); fi;
      lcmv:=valFFE(pr,lcm,z*v);
      if lcmv<>0*z then
        Add(vs,z*v); Add(tblv,lcmv*List(R,x->List(x,y->valFFE(pr,y,z*v)))^-1);
      fi;
      v:=v+1;
    od;
    return [vs,tblv];
  end;
  # return the coefficients of a polynomial over GF(pr) as a list of integers. 
  ListPolFFE:=function(d,pr,l)
    local i,p,coeff;
    coeff:=0*[1..d];
    for i in [1..Length(l)] do
      if l[i]<>0*l[i] then
        p:=IntFFE(l[i]);
        if p>(pr-1)/2 then coeff[i]:=p-pr; else coeff[i]:=p; fi;
      fi;
    od;
    return coeff;
  end;
  # Apply Chinese Remainder to a matrix of lists.
  ApplyChinRem:=function(mat1,mat2,m1,m2)
    local i,j,k,g,x,m;
    m:=m1*m2;
    g:=Gcdex(m1,m2);
    for i in [1..Length(mat1)] do 
      for j in [1..Length(mat1)] do
        for k in [1..Length(mat1[i][j])] do 
          x:=(mat1[i][j][k]*g.coeff2*m2+mat2[i][j][k]*g.coeff1*m1) mod m;
          if x>Int((m-1)/2) then mat1[i][j][k]:=x-m; else mat1[i][j][k]:=x; fi;
        od;
      od;
    od;
  end;
  ind:=H.parameter[1][1]^0;
  for i in LongestCoxeterWord(CoxeterGroup(H)) do ind:=ind*H.parameter[i][1]; od;
  sh:=ind*SchurElements(H); 
  ct:=TransposedMat(CharTable(H).irreducibles*X(Cyclotomics)^0);
  #if ForAll(Concatenation(CoxeterGroup(H).cartan),IsRat) then
  #  ctw:=ct^0;
  #  R:=ct;
  #  lcm:=Lcm(sh);
  #  d:=2*Degree(lcm)+1;
  #else
    ctw:=CharTable(CoxeterGroup(H)).irreducibles*X(Cyclotomics)^0;
    IsMat(ctw); R:=ct*ctw;
    lcm:=Size(CoxeterGroup(H))*Lcm(sh);
    d:=2*Degree(lcm)+1;
  #fi;
  # initialize matrix where we store the coefficients of polynomials.
  fertig:=false;
  pr:=1; p2:=NextPrimeInt(3*Degree(lcm)+1);
  while not fertig do
    InfoChevie("#I Trying prime : ",p2,"  \c");
    tblv:=Cetv(R,lcm,p2,d);
    if tblv<>false then 
      tbl:=[]; InfoChevie("\n#I Interpolating    : \c ");
      intp:=List([1..Length(tblv[1])],i->List(tblv[1],x->x^(i-1)));
      IsMat(intp); intp:=intp^-1;
      for i in [1..Length(ctw)] do
        tbl[i]:=[]; InfoChevie(".\c");
        for j in [1..Length(ctw)] do
          tbl[i][j]:=ListPolFFE(d,p2,List([1..d],x->tblv[2][x][i][j])*intp);
        od;
      od;
      if pr=1 then
        res:=Copy(tbl);
      else
        InfoChevie("\n#I      chinese remainder (",pr,",",p2,")   \c");
        ApplyChinRem(res,tbl,pr,p2);
      fi;
      pr:=pr*p2;
      fertig:=true; InfoChevie("     (check inverse \c");
      i:=1;
      while fertig and i<=Length(sh) do  
        j:=1; InfoChevie("\^\c");
        while fertig and j<=Length(sh) do
          x:=0*X(Cyclotomics);
          for k in [1..Length(sh)] do
            x:=x+Polynomial(Cyclotomics,res[i][k])*R[k][j];
          od;
          if (i=j and x<>lcm) or (i<>j and x<>0*x) then fertig:=false; fi;
          j:=j+1;
        od;
        i:=i+1;
      od;
      InfoChevie(" ",fertig,") \n");
      p2:=NextPrimeInt(p2);
    fi;
  od;
  if ctw<>ctw^0 then
    tbl:=ctw*List(res,x->List(x,y->Polynomial(Cyclotomics,y)));
  else
    tbl:=List(res,x->List(x,y->Polynomial(Cyclotomics,y)));
  fi;
  for i in [1..Length(sh)] do tbl[i]:=(sh[i]*tbl[i])/lcm/ind; od;
  InfoChevie("     final check \c");
  i:=1;
  while fertig and i<=Length(sh) do  
    j:=1; InfoChevie("\^\c");
    while fertig and j<=Length(sh) do
      x:=0*X(Cyclotomics);
      for k in [1..Length(sh)] do x:=x+tbl[i][k]*ct[k][j]; od;
      if (i=j and x<>sh[i]/ind) or (i<>j and x<>0*x) then fertig:=false; fi;
      j:=j+1;
    od;
    i:=i+1;
  od;
  InfoChevie(" ",fertig," \n");
  H.centralTable:=Copy(tbl);
  return tbl;
end;

###########################################################################
##
#F  'HeckeBlocks( <H>, <centbl>, <d> ) . . . . . . . . blocks of characters
##
##  'HeckeBlocks'  returns  the  block  distribution  of   the  irreducible 
##  characters of <H> for the specialization v->E(2*d),   using the central
##  table <centbl> of <H> (see 'HeckeCentralTable').
##
HeckeBlocks:=function(H,centbl,d)local c,p,bls,cetd,ds,bl,bld,phi,d,def,tbl;
  def:=HeckeCharDefect(H,d); 
  InfoChevie("#I Defects : ",def,"\n");
  cetd:=List(centbl,x->List(x,y->Value(y,E(2*d))));
  bls:=[];
  for c in Set(cetd) do 
    bl:=Filtered([1..Length(cetd)],x->cetd[x]=c);
    Add(bls,[bl,def{bl}]);
  od;
  SortBy(bls,x->x[2][1]);return bls;
end;

###########################################################################
##
#F  'HeckeModularRank( <H>, <d> ) . . . . . .  the number of simple modules
##  . . . . . . . . . . . . . . . . . . . . . . . for a specialized algebra
##
##  'HeckeModularRank'  computes the  number of  simple  modules  after the
##  specialization of $v$ to E(2 * <d>). This is done by computing the rank 
##  of the specialized character table.
##
HeckeModularRank:=function(H,d)local Koerper,v,base,rels,spec,Hcts,f,b,j;
  Koerper:=function(e,min)local i,q;
    i:=1;
    while true do
      q:=1+2*e*i;
      if q>min and IsPrime(q) then return [i,q]; fi;
      i:=i+1;
    od;
  end;
  if not ForAll(Flat(CharTable(CoxeterGroup(H)).irreducibles),IsRat) then 
    Hcts:=(TransposedMat(CharTable(CoxeterGroup(H)).irreducibles)
                           *X(Cyclotomics)^0)*CharTable(H).irreducibles;
    spec:=List(Hcts,x->List(x,y->Value(y,E(2*d))));
    return RankMat(spec);
  else
    Hcts:=CharTable(H).irreducibles;
    spec:=List(Hcts,x->List(x,y->Value(y,E(2*d))));
    v:=Koerper(d,2*Maximum(Factors(Size(CoxeterGroup(H)))));
    f:=[[],[]];
    for j in [1..v[2]] do
      f[1][j]:=j-(v[2]+1)/2;
      f[2][j]:=f[1][j]*Z(v[2])^0;
    od;
    rels:=[];
    base:=NullSpaceMat(List(Hcts,x->List(x,y->Value(y,Z(v[2])^v[1]))));
    for b in base do
      for j in [1..Length(b)] do
        b[j]:=f[1][Position(f[2],b[j])];
      od;
      if b*spec=0*b then
        Add(rels,b);
      else
        InfoChevie("wrong relation \n");
      fi;
    od;
    if Length(rels)=Length(base) then 
      return Length(Hcts)-Length(rels);
    else
      Print("#W something wrong");
    fi;
  fi;
end;

###########################################################################
##
#F  'HeckeRegularTable( <H> ) . . . . . .  character values on non-cuspidal 
##  . . . . . . . . . . . . . . . . . . . . . . . . . . and regular classes
##
##  'HeckeRegularTable'  returns  those  columns  of the character table of
##  <H>  which correspond to  non-cuspidal  classes  (i.e.,  classes  which
##  have  empty  intersection  with  every  proper  parabolic subgroup)  or  
##  classes  containing regular (or `1-good') elements.   The  point  about  
##  this table is that it is computed using the character tables of  proper  
##  parabolic  subgroups  and the character table of the underlying Coxeter 
##  group.
##
HeckeRegularTable:=function(H)local W,c,ct,r,reg;
  if not IsBound(H.regularTable) then 
    W:=CoxeterGroup(H);
    r:=CoxeterConjugacyClasses(W);
    c:=Filtered([1..Length(r)],x->Length(Set(r[x]))<W.semisimpleRank or 
                                       Length(GoodCoxeterWord(W,r[x]))=1);
    H.regularTable:=CharTable(H).irreducibles{[1..Length(r)]}{c};
  fi;
  return H.regularTable;
end;
    

###########################################################################
##
#F  'HeckeDefect1Blocks( <H>, <d> ) . . . . .  computing blocks of defect 1
##
##  'HeckeDefect1Blocks' tries to compute blocks of defect 1 of <H> without
##  using the central table of <H>. We only use  the  regular table and the
##  monomials by which T_{w_0}^2 acts (see 'HeckeCentralMonomials'). 
##
HeckeDefect1Blocks:=function(H,d)
  local p,ct,bl,W,w0,reg,i,ind,sh,x,y,res,null,null1,b,mon,sinv,phid,inv,ls;
  inv:=List(HeckeCharDefect(H,d),x->[x]);
  W:=CoxeterGroup(H);
  w0:=LongestCoxeterElement(W);
  #if  ForAll(W.generators,s->s*w0=w0*s) then 
  #  mon:=List(HeckeCentralMonomials(H),x->E(2*d)^(x.valuation/2));
  #else
    mon:=List(HeckeCentralMonomials(H),x->Value(x,E(2*d)));
  #fi;
  InfoChevie("#I Using central monomials ... \c");
  ind:=H.parameter[1][1]^0;
  for i in LongestCoxeterWord(W) do
    ind:=ind*H.parameter[i][1];
  od;
  phid:=Value(CyclotomicPolynomial(Cyclotomics,d),X(Cyclotomics)^2);
  phid:=CyclotomicPolynomial(Cyclotomics,2*d);
  sh:=List(SchurElements(H),x->Value(QuotientRemainder(x,phid)[1],E(2*d)));
  InfoChevie("Schur elements ... \c");
  for x in [1..Length(inv)] do 
    Add(inv[x],mon[x]);
    Add(inv[x],Set([sh[x],-sh[x]]));
  od;
  p:=Collected(Factors(d));
  if Length(p)=1 then
    InfoChevie("p-blocks of W ... \c");
    ct:=CharTable(W);
    ct.irreducibles:=List(CharTable(H).irreducibles,x->List(x,y->Value(y,-1)));
    bl:=PrimeBlocks(ct,p[1][1]);
    for i in [1..Length(inv)] do
      Add(inv[i],bl.block[i]);
    od;
  fi;
  InfoChevie("\n");
  ls:=Filtered([1..Length(inv)],x->inv[x][1]<>0);
  sinv:=List(Set(inv{ls}),x->Filtered(ls,y->inv[y]=x));
  reg:=List(HeckeRegularTable(H),x->List(x,y->Value(y,E(2*d))));
  res:=[];
  for x in sinv do
    null:=[];
    for y in NullSpaceMat(reg{x}) do 
      b:=0*[1..Length(mon)];
      for i in [1..Length(x)] do 
        b[x[i]]:=y[i];
      od;
      if b*reg<>0*reg[1] then return false;fi;
      Add(null,b);
    od;
    Add(res,[x,null]);
  od;
  return res;
end;

###########################################################################
##
#F  'NicePrintBlocks( <H>, <lbl>, <chnames> ) . . . . . .  print nicely the 
##  . . . . . . . . . . . . . . . . . . . . . . . . . output of HeckeBlocks 
##
NicePrintBlocks:=function(W,lbl,chnames)local a,i,j,b,d,e,l,z,bls,defs;
  Print("\\begin\{table\} \\caption\{Blocks of positive defect of \$H( )\$");
  Print(" at roots of unity\}\n");
  Print("\\label\{tab\:blocks-?? \} \\begin\{center\}\n");
  Print("\$\\begin\{array\}\{ccl\}\\hline e\&\\mbox\{defect\}\&\n");
  Print("\\mbox\{blocks of characters of positive defect\}\\\\");
  a:=ChevieCharInfo(W).a;
  for e in lbl do for l in e[2] do SortBy(l[1],x->a[x]); od; od;
  for e in lbl do
    defs:=Set(List(e[2],x->Maximum(x[2])));
    z:=1;
    Print("\n\\hline ",e[1]);
    for d in Filtered(defs,i->i>0) do 
      bls:=Filtered(e[2],x->Maximum(x[2])=d);
      Print("\& ",d," \& ");
      for b in [1..Length(bls)] do
        Print("\\{");
        for i in [1..Length(bls[b][1])] do
          Print(chnames[bls[b][1][i]]);
          if i<Length(bls[b][1]) then Print(",");fi;
        od;
        Print("\\}");
        if b<Length(bls) then Print(",\\;");else Print("\\\\");fi;
      od;
    od;
  od;
  Print("\n\\hline \\end\{array\}\$ \\end\{center\}\n\\end\{table\}\n");
end;
