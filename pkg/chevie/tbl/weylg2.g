#############################################################################
##
#A  tbl/weylg2.g             CHEVIE library         Meinolf Geck, Jean Michel
##
#Y  Copyright (C) 1994 - 2001  The CHEVIE Team
##
##  This file contains the character table of the Hecke algebra of type G2.
##
CHEVIE.AddData("CartanMat","G2",function(arg)local a,type;
  if Length(arg)>0 then type:=arg[1];else type:=1;fi;
  a:=[[2,-1],[-3,2]];
  a[1][2]:=-type;a[2][1]:=3/a[1][2];
  return a;
end);

CHEVIE.AddData("PrintDiagram","G2",function(indices,title,type)
  Print(title," ",indices[1]);
  if type=1 then Print(" >>> ");
  elif type=ER(3) then Print(" =6= ");
  else Print(" ?6? ");
  fi;
  Print(indices[2]," \n");
end);

CHEVIE.AddData("ReflectionName","G2",function(arg)local i,opt,type;
  if Length(arg)=1 then return "G2(?)"; fi;
  type:=arg[2]; opt:=arg[1];
  if type=1 then 
    if IsBound(opt.TeX) then return "G_2";
    elif IsBound(opt.arg) then return "\"G\",2";
    else return "G2";fi;
  elif type=ER(3) then 
    if IsBound(opt.TeX) then return "G_{\\hbox{sym}2}";
    elif IsBound(opt.arg) then return "\"Gsym\",2";
    else return "Gsym2";fi;
  elif IsBound(opt.TeX) then return SPrint("G_2(",Format(type^2/3,opt),")");
  elif IsBound(opt.arg) then return SPrint("\"G\",",2,",",Format(type^2/3,opt));
  else return SPrint("G2(",Format(type^2/3,opt),")");
  fi;
end);

CHEVIE.AddData("ParabolicRepresentatives", "G2", 
  s->CHEVIE.R("ParabolicRepresentatives","imp")(6,6,2,s));

# in dim 3, as in Bourbaki
CHEVIE.AddData("GeneratingRoots", "G2", [ [ 1, -1, 0 ], [ -2, 1, 1 ] ]);

CHEVIE.AddData("HyperplaneRepresentatives","G2",[1,2]);

CHEVIE.AddData("Size", "G2", 12);

CHEVIE.AddData("ReflectionDegrees","G2",[2,6]);

CHEVIE.AddData("NrConjugacyClasses","G2",6);

CHEVIE.AddData("CharInfo","G2",function()local res;
  res:=rec(charparams:=[[1,0],[1,6],[1,3,1],[1,3,2],[2,1],[2,2]],
  extRefl:=[1,5,2],
  a:=[0,6,1,1,1,1],
  A:=[0,6,5,5,5,5]);
  res.b:=List(res.charparams,x->x[2]);
  res.B:=[0,6,3,3,5,4];
  # charnames in Spaltenstein's "Sous-groupes de Borel et classes unipotentes"
  res.spaltenstein:=["1","\\varepsilon","\\varepsilon_l","\\varepsilon_c",
                     "\\theta'","\\theta''"];
  return res;
  end
);

CHEVIE.AddData("ClassNames","G2",
  ["A_0","\\tilde A_1","A_1","G_2","A_2","A_1+\\tilde A_1"]);

CHEVIE.AddData("ClassInfo","G2",
  rec(classtext:=[[],[2],[1],[1,2],[1,2,1,2],[1,2,1,2,1,2]],
    classnames:=CHEVIE.R("ClassNames","G2"),
    classparams:=CHEVIE.R("ClassNames","G2"),
    orders:=[1,2,2,6,3,2],
    classes:=[1,3,3,2,2,1]));

CHEVIE.AddData("PowerMaps","G2",[,[1,1,1,5,5,1], [1,2,3,6,1,6]]);

# description of fake degrees by powers of q
CHEVIE.AddData("sparseFakeDegrees","G2",
  [[1,0],[1,6],[1,3],[1,3],[1,1,1,5],[1,2,1,4]]);

#############################################################################
##
#F  CHEVIE.R("ClassParameter","G2")( <w> )    . . . class parameter of w
##
##  given a Coxeter word for w in G_2,  returns the name of w's conjugacy class.
##
CHEVIE.AddData("ClassParameter","G2",w->CHEVIE.R("ClassNames","G2")
  [PositionProperty([[[]],[[2],[1,2,1],[2,1,2,1,2]],[[1],[2,1,2],[1,2,1,2,1]],
    [[2,1],[1,2]],[[2,1,2,1],[1,2,1,2]],[[1,2,1,2,1,2]]],x->w in x)]);

##  SchurElements, and CharTable and Representations of Hecke algebra
##  need the square root of the product of the parameters for 2-dimensional
##  representations.
##  If  sqrtpara are not bound, the function squv knows how to compute this
##  root only if u and v are equal or one is the cube of the other one.
##
CHEVIE.AddData("squv","G2",function(para,sqrtpara)
  local u,v;
  u:=Product(para[1]); v:=Product(para[2]);
  if u=v then return u;
  elif u=v^3 then return -v^2;
  elif v=u^3 then return -u^2;
  elif IsBound(sqrtpara[1]) and IsBound(sqrtpara[2]) then
       return sqrtpara[1]*sqrtpara[2];
  else return GetRoot(u*v,2,"Hecke(G2)");
  fi;
end);

CHEVIE.AddData("HeckeCharTable","G2",function(para,sqrtpara)
  local x,y,z,t,tbl,f1,f2,one;
  x:=para[1][1];y:=para[1][2];z:=para[2][1];t:=para[2][2];
  one:=(x*y*z*t)^0;
  f1:=function(u,v)return [1,v,u,v*u,v^2*u^2,v^3*u^3]*one;end;
  f2:=function(x,y,z,t,eps)local squv;
    squv:=eps*CHEVIE.R("squv","G2")(para,sqrtpara);
    return [2,z+t,x+y,-squv,-x*y*z*t,2*squv^3]*one;end;
  tbl:=rec( identifier:="H(G2)", parameter:=[[x,y],[z,t]], size:=12,
    powermap:=CHEVIE.R("PowerMaps","G2"),
    irreducibles:=
     [f1(x,z),f1(y,t),f1(y,z),f1(x,t),f2(x,y,z,t,1),f2(x,y,z,t,-1)],
    irredinfo:=CHEVIE.R("IrredInfo","G2"));
  Inherit(tbl,CHEVIE.R("ClassInfo","G2"));
  tbl.centralizers:=List(tbl.classes,x->tbl.size/x);
  tbl := CHEVIE.compat.MakeCharacterTable(tbl);
  return tbl;
end);

CHEVIE.AddData("HeckeRepresentation","G2",function(para,sqrtpara,i)
  local one,squv,x,y,z,t; one:=Product(para[1])^0*Product(para[2])^0;
  x:=para[1][1];y:=para[1][2];z:=para[2][1];t:=para[2][2];
  if i=1 then return [[[x]],[[z]]]*one;
  elif i=2 then return [[[y]],[[t]]]*one;
  elif i=3 then return [[[y]],[[z]]]*one;
  elif i=4 then return [[[x]],[[t]]]*one;
  else
    squv:=CHEVIE.R("squv","G2")(para,sqrtpara);
    if i=6 then squv:=-squv;fi;
    return [[[y,-1], [0,x]],[[z,0],[squv+y*z+x*t,t]]]*one;
  fi;
end);

CHEVIE.AddData("Representation","G2",function(i)local para;
  return CHEVIE.R("HeckeRepresentation","G2")([[1,-1],[1,-1]],[1,1],i);
end);

#############################################################################
##
#F  PoincarePolynomial(u,v) . . . . Poincare Polynomial for G2
##
##  Same as Schur element for the identity character.
##
CHEVIE.AddData("PoincarePolynomial","G2",function(param)local u,v;
  u:=-param[1][1]/param[1][2]; v:=-param[2][1]/param[2][2];
  return (1+u)*(v+1)*(1+u*v+u^2*v^2);
end);

CHEVIE.AddData("SchurModels","G2",rec(
  f1:=rec(vcyc:=[[[1,-1,0,0],1],[[0,0,1,-1],1],[[1,-1,1,-1],3]]),
  f2:=rec(coeff:=-2,root:=[1,-1,1,-1]/2,factor:=[-1,1,0,0],
    vcyc:=[[[0,0,0,0,1],3],[[0,0,-1,1,1],3]])));

CHEVIE.AddData("SchurData","G2",[
  rec(name:="f1",order:=[1,2,3,4]),rec(name:="f1",order:=[2,1,4,3]),
  rec(name:="f1",order:=[2,1,3,4]),rec(name:="f1",order:=[1,2,4,3]),
  rec(name:="f2",order:=[1,2,3,4],rootPower:=-1),
  rec(name:="f2",order:=[1,2,3,4],rootPower:=1)]);

#############################################################################
##
#F  SchurElement( <phi>, <para>, <sqrtpara> ) . .  Schur elements for type G2
##
CHEVIE.AddData("SchurElement","G2",function(phi,para,sqrtpara)local u,v,squv,p;
  u:=-para[1][1]/para[1][2]; v:=-para[2][1]/para[2][2];
  p:=Position(CHEVIE.R("CharInfo","G2")().charparams,phi);
  if p=1 then return (1+u)*(v+1)*(u^2*v^2+u*v+1);
  elif p=2 then return (1+u)*(v+1)*(u^2*v^2+u*v+1)/u^3/v^3;
  elif p=3 then return (u^2+v^2+u*v)*(1+u)*(v+1)/u^3;
  elif p=4 then return (u^2+v^2+u*v)*(1+u)*(v+1)/v^3;
  fi;
  squv:=CHEVIE.R("squv","G2")(para,sqrtpara)/para[1][2]/para[2][2];
  if p=6 then squv:=-squv;fi;
  return 2*(u*v)^-1*(u*v+1+squv)*(u+v-squv);
end);

CHEVIE.AddData("UnipotentCharacters","G2",function()
  return rec(
  harishChandra:=[
    rec(relativeType:=rec(series:="G",indices:=[1..2],rank:=2),
 levi:=[], parameterExponents:=[1,1],
        charNumbers:=[1..6], 
  eigenvalue:=1, cuspidalName:=""),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
levi:=[1..2], parameterExponents:=[],
	charNumbers:=[10], eigenvalue:=E(3)^2, cuspidalName:="G_2[\\zeta_3^2]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
levi:=[1..2], parameterExponents:=[],
	charNumbers:=[7], eigenvalue:=-1, cuspidalName:="G_2[-1]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
levi:=[1..2], parameterExponents:=[],
	charNumbers:=[9], eigenvalue:=E(3), cuspidalName:="G_2[\\zeta_3]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
levi:=[1..2], parameterExponents:=[],
	charNumbers:=[8], eigenvalue:=1, cuspidalName:="G_2[1]")],
  families:=[
    Family("S3",[5,6,4,3,8,7,9,10],rec(ennola:=-5)),
    Family("C1",[1]),
    Family("C1",[2])],
  a:=[0,6,1,1,1,1,1,1,1,1],
  A:=[0,6,5,5,5,5,5,5,5,5],
  charSymbols:=[[[0],[0],[0],[0],[0],[2]],[[0,1],[0,1],[0,1],[0,1],[0,1],[1,2]],
    [[0],[0],[1],2,0],[[0],[0],[1],2,1],
    [[0],[0],[0],[0],[1],[1]],[[0],[0],[0],[1],[0],[1]],
    [[0,1],[0],[0,1],[],[0],[]],[[0,1],[0,1],[0],[],[],[0]],
    [[0,1],[0],[0],[0,1],[],[]],[[0,1],[0,1],[],[0],[0],[]]]);
end);

CHEVIE.AddData("Invariants","G2",
 [function(x,y)return -3*x*y+3*x^2+y^2;end,
  function(x,y)return x^2*y^4-6*x^3*y^3+13*x^4*y^2-12*x^5*y+4*x^6;end]);

CHEVIE.AddData("Discriminant","G2",function()return
  function(x,y)return 4*x^3*y-27*y^2;end;end);

# c is CartanType -- does not matter
CHEVIE.AddData("UnipotentClasses","G2",function(c,p)local uc,Z;
  if p=0 then p:=1;fi; Z:=n->ComplexReflectionGroup(n,1,1);
  uc:=rec(classes:=[
   rec(name:="1",succ:=["A1"],dynkin:=[0,0],balacarter:=[],
     red:=CoxeterGroup("G",2)),
   rec(name:="A_1",succ:=["~A1"],dynkin:=[1,0],balacarter:=[1],red:=Z(2)),
   rec(name:="\\tilde A_1",succ:=["G2(a1)"],dynkin:=[0,1],balacarter:=[2],
     red:=Z(2-(Gcd(p,3)-1)/2)),
   rec(name:="G_2(a_1)",succ:=["G2"],dynkin:=[2,0],balacarter:=[1,-2],
     Au:=CoxeterGroup("A",2-(Gcd(p,3)-1)/2)),
   rec(name:="G_2",succ:=[],dynkin:=[2,2],Au:=Z(Gcd(p,6)),balacarter:=[1,2])],
springerSeries:=[rec(relgroup:=CoxeterGroup("G",2),levi:="",Z:=[],
  locsys:=[[5,1],[1,1],[4,2],[2,1],[4,3],[3,1]]),
  rec(relgroup:=CoxeterGroup(),levi:=[1,2],Z:=[],locsys:=[[4,1]],
   parameter:=[8])]); # Fourier transform of 8th unip. character
  if p=2 then uc.springerSeries[1].locsys[1]:=[5,2];
    Add(uc.springerSeries,
      rec(relgroup:=CoxeterGroup(),levi:=[1,2],Z:=[],locsys:=[[5,1]]));
  elif p=3 then 
    Add(uc.classes,rec(name:="(\\tilde A_1)_3",
      succ:=["~A1"],dimBu:=3,red:=Z(2),Au:=CoxeterGroup()));
    Add(uc.classes[1].succ,"(~A1)3");
    uc.classes[3].dimBu:=2;Unbind(uc.classes[3].dynkin);
    uc.springerSeries[1].locsys{[3,5]}:=[[6,1],[4,2]];
    for c in [2,3] do Add(uc.springerSeries,
      rec(relgroup:=CoxeterGroup(),levi:=[1,2],Z:=[],locsys:=[[5,c]]));od;
  fi;
  uc.orderClasses:=List(uc.classes,c->List(c.succ,
    n->PositionProperty(uc.classes,c->UnipotentClassOps.Name(c)=n)));
  for c in uc.classes do Unbind(c.succ);
    if not IsBound(c.red) then c.red:=Z(1);fi;
    if not IsBound(c.Au) then c.Au:=Z(1);fi;
    c.AuAction:=ExtendedReflectionGroup(c.red,List(c.Au.generators,
        x->IdentityMat(c.red.rank)));
  od;
  return uc;
end);

CHEVIE.AddData("KLeftCellRepresentatives","G2",
[rec(character:=[1], duflo:=[1,2], reps:=""),
 rec(character:=[2], duflo:=[7,8], reps:=""),
 rec(character:=[3,5,6], duflo:=[5,8], reps:=[[6,10],[12,3]]),
 rec(character:=[4,5,6], duflo:=[7,3], reps:=[[5,10],[12,4]])]);
