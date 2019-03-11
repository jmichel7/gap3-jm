#############################################################################
##
#A  tbl/coxi.g               CHEVIE library      Meinolf Geck, Jean Michel
##
#Y  Copyright (C) 1992 - 2010  The CHEVIE Team
##
##  This file contains the character table of the Hecke algebras of type I.
##
##  The  characters of the groups I_2(m) are parametrized by [d,b] where
##  d  is their degree and b is  their b-function, i.e. the valuation of
##  the  fake degrees. This identifies uniquely the characters when m is
##  odd.  When m is  even, there are  two characters with [d,b]=[1,m/2].
##  The  one  which  maps  the  generators  [1,2]  to  [1,-1] is denoted
##  [1,m/2,"'"]  and  the  one  which  maps  [1,2]  to [-1,1] is denoted
##  [1,m/2,"''"].
##
CHEVIE.AddData("CartanMat","I",function(arg)local bond,type,m;
  m:=[[2,0],[0,2]];
  bond:=arg[1];
  if bond=2 then return m;fi;
  if Length(arg)=2 then type:=arg[2];
  elif bond mod 2=0 then type:=1;
  else type:=E(2*bond)+E(2*bond)^-1;
  fi;
  m[1][2]:=-type;m[2][1]:=(2+E(bond)+E(bond)^-1)/m[1][2];return m;
end);

CHEVIE.AddData("PrintDiagram","I",function(arg)local bond,indices,type;
  Print(arg[3]," ");bond:=arg[1];indices:=arg[2];
  if Length(arg)=4 then type:=arg[4];else type:=E(2*bond)+E(2*bond)^-1;fi;
  if type=E(2*bond)+E(2*bond)^-1 then 
    Print(indices[1]," -",String(bond),"- ",indices[2],"\n");
  elif type=1 then  Print(indices[1]," >",String(bond),"> ",indices[2],"\n");
  else Print(indices[1]," ?",String(bond),"? ",indices[2],"\n");
  fi;
end);

CHEVIE.AddData("ReflectionName","I",function(arg)local bond,type,opt;
  bond:=arg[1];opt:=arg[2];
  if Length(arg)=3 then type:=arg[3];
  elif bond mod 2=0 then type:=1;
  else type:=E(2*bond)+E(2*bond)^-1;
  fi;
  if type=1 then 
    if IsBound(opt.TeX) then return SPrint("I_2(",bond,")");
    elif IsBound(opt.arg) then return SPrint("\"I\",2,",bond);
    else return SPrint("I2(",bond,")");fi;
  elif type=E(2*bond)+E(2*bond)^-1 then 
    if bond mod 2=1 then 
      if IsBound(opt.TeX) then return SPrint("I_2(",bond,")");
      elif IsBound(opt.arg) then return SPrint("\"I\",2,",bond);
      else return SPrint("I2(",bond,")");fi;
    else 
      if IsBound(opt.TeX) then return SPrint("I_{\\hbox{sym}2}(",bond,")");
      elif IsBound(opt.arg) then return SPrint("\"Isym\",2,",bond);
      else return SPrint("Isym2(",bond,")");fi;
    fi;
  elif IsBound(opt.TeX) then
       return SPrint("I_?(",Format(type^2/(2+E(bond)+E(bond)^-1),opt),
                     ")(",bond,")");
  elif IsBound(opt.arg) then return SPrint("\"Isym\",2,",bond,",",
     Format(type^2/(2+E(bond)+E(bond)^-1),opt));
  else return SPrint("I?(",type^2/(2+E(bond)+E(bond)^-1),")(",bond,")");
  fi;
end);

CHEVIE.AddData("SemisimpleRank","I",2);

CHEVIE.AddData("GeneratingRoots", "I", function(m)local a, b, r;
  a := E(2*m)^(m-1);
  b := ComplexConjugate(a);
  if m mod 2 = 0 then r := ER(m/2); else r := 1; fi;
  return [[1, 0], [r*(a+b)/2, r*(a-b)/2/E(4)]];
end);

CHEVIE.AddData("EigenvaluesGeneratingReflections","I",m->[-1,-1]);

CHEVIE.AddData("Size", "I", function(arg) return 2*arg[1]; end);

CHEVIE.AddData("ReflectionDegrees","I",m->[2,m]);

CHEVIE.AddData("NrConjugacyClasses","I", m->QuoInt(m+3,2)+((m+1)mod 2)*2);

CHEVIE.AddData("ParabolicRepresentatives", "I", function(m,s)
  return CHEVIE.R("ParabolicRepresentatives","imp")(m,m,2,s);end);

#how to make a .charname from a .charparam
CHEVIE.AddData("CharName","I",function(m,x,option)local s;
  if IsList(x[1]) then return PartitionTupleToString(x);
  else
    if IsBound(option.TeX) then s:="\\phi";else s:="phi";fi;
    s:=SPrint(s,"{",x[1],",",x[2],"}");
    if Length(x)=3 then Append(s,x[3]);fi;
    return String(s);
  fi;
end);

CHEVIE.AddData("CharInfo","I",function(m)local res,applyf,v;
  res:=rec(charparams:=[[1,0]]);
  if m mod 2=0 then res.extRefl:=[1,5,4]; 
     Append(res.charparams,[[1,m/2,"'"],[1,m/2,"''"]]);
  else res.extRefl:=[1,3,2];fi;
  Add(res.charparams,[1,m]);
  Append(res.charparams,List([1..QuoInt(m-1,2)],i->[2,i]));
  res.b:=List(res.charparams,x->x[2]);
  res.B:=List(res.charparams,function(phi)
    if phi[1]=1 then return phi[2]; else return m-phi[2]; fi;
  end);
  res.a:=List(res.charparams,function(phi)
     if phi[1]<>1 or phi[2]=m/2 then return 1;else return phi[2]; fi;
  end);
  res.A:=List(res.charparams,function(phi)
     if phi[1]=1 or phi[2]=m/2 then return m-1;else return phi[2];fi;
  end);
  res.charSymbols:=List([1..QuoInt(m-1,2)],function(l)local S,k; 
      S:=List([1..m],i->[0]);k:=0;
      if k<>0 then S[1]:=[0,1];S[1+((k+l)mod m)]:=[0,1];S[k+1]:=[];S[l+1]:=[];
      else S[1]:=[1];S[l+1]:=[1];
      fi;
      return S;end);
  v:=List([1..m],x->[0]);v[m]:=[1,2];
  res.charSymbols:=Concatenation([v],res.charSymbols);
  if m mod 2=0 then
    v:=List([1..m],x->[0]);v[m]:=[1];v[m/2]:=[1];
    res.charSymbols:=Concatenation([v],res.charSymbols);
    v:=List([1..m],x->[0]);v[m]:=[1];v[m/2]:=[1];
    res.charSymbols:=Concatenation([v],res.charSymbols);
  fi;
  v:=List([1..m],x->[0,1]);v[m]:=[2];
  res.charSymbols:=Concatenation([v],res.charSymbols);
  res.malleParams:=List(res.charSymbols,x->List(x,PartBeta));
  if m mod 2=0 then
    res.malleParams[2]:=Concatenation(res.malleParams[2]{[1..m/2]},[1]);
    res.malleParams[3]:=Concatenation(res.malleParams[3]{[1..m/2]},[-1]);
  fi;
  return res;
end);

CHEVIE.AddData("WordsClassRepresentatives", "I", function(m) local r, x, i;
  if IsInt(m/2) then r:=[[],[1],[2]]; else r:=[[],[1]]; fi;
  x:=[1,2];
  for i in [1..QuoInt(m,2)] do Add(r, ShallowCopy(x)); Append(x,[1,2]); od;
  return r;
end);

#############################################################################
##
#F  CHEVIE.R("ClassInfo","I")( <m> ) . .. conjugacy classes for type I.
##
CHEVIE.AddData("ClassInfo","I",function(m)local r,i,clnp,cl,g1,g2,gen,perm;
  r := CHEVIE.R("WordsClassRepresentatives", "I")(m);
  clnp := List(r,IntListToString);
  g1:=();i:=2; while 2*i<=m+1 do g1:=g1*(i,m-i+2); i:=i+1;od;
  g2:=();i:=1; while 2*i<=m do g2:=g2*(i,m-i+1); i:=i+1;od;
  gen := [g1,g2];
  perm := function(l) if Length(l)=0 then return ();else return
                      Product(gen{l});fi;end;
  if m mod 2=0 then cl:=[1,m/2,m/2]; Append(cl,[1..m/2-1]*0+2);Add(cl,1);
  else cl:=[1,m];Append(cl,[1..(m-1)/2]*0+2);
  fi;
  return rec(classtext:=r, classnames:=clnp, classparams:=clnp,
    orders:=List(r,i->OrderPerm(perm(i))), classes:=cl);
end);

##########################################################################
##
#F  HeckeCharTable( <m>, param, rootparam)  . . . . . . . . .
##  . . . . . . . .  character table of the Hecke algebra of dihedral type
##
##
CHEVIE.AddData("HeckeCharTable","I",function(m,param,rootparam)
  local u, v, squv, cl, r, ct, tbl;
  u:=-param[1][1]/param[1][2]; v:=-param[2][1]/param[2][2];
  if m mod 2<>0 then squv:=u;
  elif IsBound(rootparam[1]) and IsBound(rootparam[2]) then
       squv:=rootparam[1]*rootparam[2];
  else squv:=GetRoot(u*v,2,"CharTable(Hecke(I2(",m,")))");
  fi;
  ct:=[[u,v]];if m mod 2=0 then Append(ct,[[u,-u^0],[-v^0,v]]);fi;
  Add(ct,[-v^0,-v^0]);
  cl:=CHEVIE.R("ClassInfo","I")(m); r:=cl.classtext;
  ct:=List(ct,i->List(r,x->Product(i{x})));
  Append(ct,List([1..QuoInt(m-1,2)],function(j)local l,i,k;
    l:=[];
    for i in [1..Length(r)] do
      k:=Length(r[i])/2;
      if r[i]=[] then l[i]:=2*v^0;
      elif r[i]=[1] then l[i]:=u-1;
      elif r[i]=[2] then l[i]:=v-1;
      else l[i]:=squv^k*(E(m)^(k*j)+E(m)^(-k*j));
      fi;
    od;
  return l;end));
  tbl:=rec(identifier:=SPrint("H(I2(",m,"))"),cartan:=CartanMat("I",2,m), 
    size:=2*m,
    irredinfo:=List(CHEVIE.R("CharInfo","I")(m).charparams,x->rec(
       charparam:=x,charname:=CHEVIE.R("CharName","I")(m,x,rec(TeX:=true)))),
    parameter:=[u,v], powermap:=[], irreducibles:=ct*v^0);
  Inherit(tbl,cl);
  tbl.centralizers:=List(tbl.classes,i->tbl.size/i);
  tbl:=CHEVIE.compat.MakeCharacterTable(tbl);
  CHEVIE.compat.AdjustHeckeCharTable(tbl,param);
  return tbl;
end);

CHEVIE.AddData("Representation","I",function(m,i)return 
  CHEVIE.R("HeckeRepresentation","I")(m,[[1,-1],[1,-1]],[1,1],i);
end);

CHEVIE.AddData("HeckeRepresentation","I",function(m,param,rootparam,i)
 local u,v,squv;
  if i=1 then return [[[param[1][1]]],[[param[2][1]]]];fi;
  if m mod 2=0 then i:=i-2;fi;
  if i=0 then return [[[param[1][1]]],[[param[2][2]]]];
  elif i=1 then return [[[param[1][2]]],[[param[2][1]]]];
  elif i=2 then return [[[param[1][2]]],[[param[2][2]]]];
  else 
    u:=-param[1][1]/param[1][2]; v:=-param[2][1]/param[2][2];
    if m mod 2<>0 then squv:=u;
    elif IsBound(rootparam[1]) and IsBound(rootparam[2]) then
         squv:=rootparam[1]*rootparam[2];
    else squv:=GetRoot(u*v,2,"Representation(Hecke(I2(",m,")),[",i,"])");
    fi;
    return [-[[-u^0,u^0],[0*u,u]]*param[1][2],
     -[[v,0*v],[u+v+squv*(E(m)^(i-2)+E(m)^(2-i)),-v^0]]*param[2][2]];
  fi;
end);

CHEVIE.AddData("Frobenius","I",function(m,sqrtu,j)
  return [[0,1/sqrtu/(E(2*m)^j+E(2*m)^-j)],
                           [sqrtu*(E(2*m)^j+E(2*m)^-j),0]]*sqrtu^0;
end);

CHEVIE.AddData("PoincarePolynomial","I",function(m,param)local u,v;
  u:=-param[1][1]/param[1][2];v:=-param[2][1]/param[2][2];
  if IsInt(m/2) then return Sum([1..m/2],i->(u*v)^(i-1))*(u+1)*(v+1);
  else return Sum([1..m],i->u^(i-1))*(u+1);
  fi;
end);

#############################################################################
##
#F  SchurElement( <m>, <phi>, <u> [, <v> [,<squv>]] )   . . . . . . .
#F   . . . . . . . . . . . . . . . . . . Schur element for dihedral type
##
##  'CHEVIE.R("SchurElement","I")'  returns  the  constant arising
##  from  the  Schur  relations  for  the  irreducible  characters  with
##  charparam phi of the Hecke algebra of type I with parameters u, v.
##
##  The  constant corresponding to phi  is defined to be  P/D where P is
##  the Poincare polynomial and D is the generic degree of phi.
##
CHEVIE.AddData("SchurElement","I",function(m,phi,para,rootpara)
  local u,v,ruv,e,ci;
  if m mod 2=1 then 
    ci:=CHEVIE.R("CharInfo","I")(m);
    ci:=ci.malleParams[Position(ci.charparams,phi)];
    return CHEVIE.R("SchurElement","imp")(m,1,2,ci,
             [List([0..m-1],i->E(m)^i),para[2]],[])/m; 
  fi;
  u:=-para[1][1]/para[1][2];v:=-para[2][1]/para[2][2];
  if phi[1]=1 then
    if phi[2]=m/2 then e:=Sum([0..m/2-1],i->(u/v)^i)*(u+1)*(v+1)/v;
      if phi[3]="'" then return e;else return (v/u)^(m/2)*e;fi;
    else e:=Sum([0..m/2-1],i->(u*v)^i)*(u+1)*(v+1);
      if phi[2]=0 then return e;else return (u*v)^(-m/2)*e;fi;
    fi;
  else e:=E(m)^phi[2]+E(m)^-phi[2];
    if ForAll([1,2],i->IsBound(rootpara[i])) then ruv:=Product(rootpara);
    else ruv:=GetRoot(u*v,2,"SchurElement(Hecke(I2(",m,"),",phi,"))");
    fi;
    return -m*(u*v+1-ruv*e)*(u+v+e*ruv)/(u*v*(e^2-4));
  fi;
end);

CHEVIE.AddData("FakeDegree","I",function(m,phi,q)
  if phi[1]=1 then return q^phi[2];else return q^phi[2]+q^(m-phi[2]);fi;
end);

CHEVIE.AddData("CharTable","I",function(m)local res;
  res:=CHEVIE.R("HeckeCharTable","I")(m,[[1,-1],[1,-1]],[1,1]);
  res.identifier:=SPrint("W(I2(",m,"))");
  return res;
end);

CHEVIE.AddData("DecompositionMatrix","I",function(n,p)local T, m; 
  T:=CHEVIE.R("CharTable","I")(n); T.name:=T.identifier; 
  m:=DecompositionMatrix(T mod p);
  return List(BlocksMat(m),c->[c[1],m{c[1]}{c[2]}]);
end);

CHEVIE.AddData("FactorizedSchurElement","I",function(arg)local ci;
  if arg[1] mod 2=0 and arg[3][1]<>arg[3][2] then Error("not implemented");fi;
  ci:=CHEVIE.R("CharInfo","I")(arg[1]);
  ci:=ci.malleParams[Position(ci.charparams,arg[2])];
  return CHEVIE.R("FactorizedSchurElement","imp")
    (arg[1],arg[1],2,ci,arg[3],1);
end);

CHEVIE.AddData("Invariants","I",function(arg)local e,type,m;
  e:=arg[1];
  if Length(arg)=2 then type:=arg[2];
  elif e mod 2=0 then type:=1;
  else type:=-E(e)^((e+1)/2)-E(e)^((e+3)/2);
  fi;
  m:=DiagonalMat(1+E(e)^-1,-type)*
    CHEVIE.R("GeneratingRoots","imp")(e,e,2);
  return List(CHEVIE.R("Invariants","imp")(e,e,2),f->function(arg)
    return ApplyFunc(f,arg*m);end);
end);

# the next function is only for cuspidal symbols
CHEVIE.AddData("SymbolToParameter","I",function(S)
  if S[1]<>[0,1] or not [] in S then return false;fi;
  if Length(S) mod 2=1 then S:=Reversed(S);
    return [Position(S,[]), Position(S,[0,1])-Position(S,[])];
  else return Position(S,[])+[-Position(S{[2..Length(S)]},[0,1]),0]-1;
  fi;
end);

# The symbols returned are rotations of those given by Gunter.
# They are reduced in the sense of SymbolsDefect(e,2,0,0)
CHEVIE.AddData("ParameterToSymbol","I",function(e,p)local S;
  if p=[0] then S:=List([1..e],x->[0]);S[e]:=[2];
  elif p=[1] then S:=List([1..e],x->[0,1]);S[e]:=[1,2];
  elif Length(p)=3 then S:=List([1..e/2-1],x->[0]);Append(S,[[1],2,(p[3]+1)/2]);
  elif e mod 2=0 then S:=List([1..e],x->[0]);
    if p[1]=0 then S{[e,e-p[2]]}:=[[1],[1]];
    else S{1+[0,(p[2]-p[1])mod e]}:=[[0,1],[0,1]];
         S{1+[-p[1]mod e,p[2]]}:=[[],[]];
    fi;
  else S:=List([1..e],i->[0]);
    if p[1]<>0 then S{1+[0,-Sum(p)mod e]}:=[[0,1],[0,1]];
                    S{1+List(-p,x->x mod e)}:=[[],[]];
    else S{e+[-((p[2]-p[1]) mod e),0]}:=[[1],[1]];
    fi;
  fi;
  return S;
end);

CHEVIE.AddData("UnipotentCharacters","I",function(e)local cusp,uc,f;
  f:=QuoInt(e,2); uc:=rec();
  uc.harishChandra:=[rec(
    relativeType:=rec(series:="I",indices:=[1,2],rank:=2,bond:=e),
    parameterExponents:=[1,1], levi:=[], eigenvalue:=1,
    cuspidalName:="")];
  if e mod 2<>0 then uc.harishChandra[1].charNumbers:=[1..f+2];
  else uc.harishChandra[1].charNumbers:=Concatenation([1,3,4,2],4+[1..f-1]);
  fi;
# For I2(e) there are 3 families: Id, St and a big one.
# in the big one the cuspidal chars are S(k,l) where 0<k<l<e-k
  cusp:=Concatenation(List([1..f-1],k->List([k+1..e-k-1],l->[k,l])));
  f:=f+1-(e mod 2); # number principal series chars in big family
  Append(uc.harishChandra,List([1..Length(cusp)],x->
     rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
         parameterExponents:=[], levi:=[1,2],
         eigenvalue:=E(e)^-Product(cusp[x]),
         cuspidalName:=SPrint("I_2(",e,")",FormatGAP(cusp[x])),
         charNumbers:=[x+f+2])));
  uc.families:=[Family(CHEVIE.families.Dihedral(e),[1..Length(cusp)+f]+2),
     Family("C1",[1]),Family("C1",[2])];
  uc.parameters:=Concatenation([[0],[1]],uc.families[1].parameters);
  uc.charSymbols:=List(uc.parameters,p->CHEVIE.I.ParameterToSymbol(e,p));
# for S(k,l) the b is min(k+l,e-k-l)
  uc.a:=Concatenation([0,e],List(uc.families[1].parameters,x->1));
  uc.A:=Concatenation([0,e],List(uc.families[1].parameters,x->e-1));
  if e=5 then uc.curtis:=[2,1,3,4,6,5];fi;
  return uc;
end);
