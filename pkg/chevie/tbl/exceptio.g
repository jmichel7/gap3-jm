#############################################################################
##
#A  tbl/exceptio.g             CHEVIE library                     Jean Michel
##
#Y  Copyright (C) 1999 - 2017  Lehrstuhl D fur Mathematik, RWTH Aachen,
#Y  University of St. Andrews, and  University Paris VII.
##
##  This file contains various code which is common to several (but not all)
##  types of reflection groups.
##

CHEVIE.IndirectAddData("CharName",["2E6","E6","E7","E8","2F4","F4","G2",
  "H3","H4","2G5","G24","G25","G26","G27","G29","G31","G32","G33","G34"],
  t->function(x,option)local s,f;
    for f in ["frame","kondo","spaltenstein","gp"] do
      if IsBound(option.(f)) then s:=CHEVIE.R("CharInfo",t)();
        if IsBound(s.(f)) then 
          s:=s.(f)[Position(s.charparams,x)];
          if IsBound(option.TeX) then return s; else return TeXStrip(s);fi;
        fi;
      fi;
    od;
    if IsBound(option.TeX) then s:="\\phi_";else s:="phi";fi;
    Append(s,SPrint("{",x[1],",",x[2],"}"));
    if Length(x)=3 then Append(s,List([1..x[3]],y->'\''));fi;
    return String(s);
end);

CHEVIE.IndirectAddData("IrredInfo",["G24","G25","G26","G27","G29","G31","G32",
  "G33","G34","H3","H4","2E6","2F4","3D4","E6","E7","E8","F4","G2"],
  t->List(CHEVIE.R("CharInfo",t)().charparams,x->
    rec(charparam:=x,charname:=CHEVIE.R("CharName",t)(x,rec(TeX:=true)))));

#CHEVIE.AddData("ClassName",["H3","H4","G4_22",
#  "G24","G25","G26","G27","G29","G31","G32","G33","G34"],x->x);

CHEVIE.IndirectAddData("CartanMat",["G25","G26","G29","G31","G32","G34"],
  function(t)local r,rbar,e;
  r:=CHEVIE.R("GeneratingRoots",t);
  rbar:=ComplexConjugate(r);
  e:=CHEVIE.R("EigenvaluesGeneratingReflections",t);
  e:=1-List(e,x->E(Denominator(x))^Numerator(x));
  e:=List([1..Length(e)],i->e[i]*rbar[i]/(rbar[i]*r[i]));
  return List(e,x->List(r,y->x*y));
end);

CHEVIE.IndirectAddData("ReflectionName",
  ["G24","G25","G26","G27","G29","G31","G32","G33","G34","E6","E7","E8","2E6",
   "2F4","3D4","H3","H4"],
  t->function(option)local i,o;
  i:=["G24","G25","G26","G27","G29","G31","G32","G33","G34","E6","E7","E8",
   "2E6","2F4","3D4","H3","H4"];
  o:=["G_{24}","G_{25}","G_{26}","G_{27}","G_{29}","G_{31}","G_{32}","G_{33}",
   "G_{34}","E_6","E_7","E_8","{}^2E_6","{}^2F_4","{}^3D_4","H_3","H_4"];
  if IsBound(option.TeX) then return o[Position(i,t)];else return t;fi;end);

CHEVIE.IndirectAddData("ReflectionName",["A","D","2A","2D"],
t->function(r,option)local i,o;
  i:=["A","D","2A","2D"];o:=["A","D","{}^2A","{}^2D"];
  if IsBound(option.arg) then return SPrint(FormatGAP(t),",",r);
  elif IsBound(option.TeX) then 
       return SPrint(o[Position(i,t)],"_",TeXBracket(r));
  else return SPrint(t,r);fi;
  end);

CHEVIE.IndirectAddData("CharTable",["3D4","E6","2E6","E7","E8",
 "F4","2F4","G2","H3","H4"],
 t->function()local res,rank;
  rank:=Position("12345678",t[Length(t)]);
  res:=CHEVIE.R("HeckeCharTable",t)(List([1..rank],x->[1,-1]),
                                          List([1,rank],x->1));
  CHEVIE.compat.ChangeIdentifier(res, SPrint("W(",t,")"));
  return res;
end);

CHEVIE.IndirectAddData("PoincarePolynomial",["G24","G27","G29","G33","G34",
  "H3","H4","E6","E7","E8"],
t->function(q)return Product(CHEVIE.R("ReflectionDegrees",t),
  x->Sum([0..x-1],y->(-q[1][1]/q[1][2])^y));
end);

CHEVIE.IndirectAddData("Representation",["G24","G25","G26","G27","G29"],
  t->function(i)local para;
  para:=CHEVIE.R("EigenvaluesGeneratingReflections",t);
  para:=List(para,x->List([0..1/x-1],j->E(1/x)^j));
  return CHEVIE.R("HeckeRepresentation",t)(para,[],i);
end);

CHEVIE.IndirectAddData("SemisimpleRank",
["G2","F4","H3","E6","G24","G25","G26","G27","G29","G31","G32","G33","G34"], 
function(t)local r;r:=CHEVIE.R("GeneratingRoots",t);
  if IsFunc(r) then r:=r();fi;
  return Length(r[1]);end);

CHEVIE.IndirectAddData("SemisimpleRank",["A","B","D"],t->(r->r));

CHEVIE.IndirectAddData("FakeDegree",
["G2","F4","H3","E6","G24","G25","G26","G27","G29","G32","G33","G34"], 
  t->function(phi,q)local f;f:=CHEVIE.R("sparseFakeDegrees",t)
   [Position(CHEVIE.R("CharInfo",t)().charparams,phi)];
   return Sum([1,3..Length(f)-1],i->f[i]*q^f[i+1]);
end);

CHEVIE.IndirectAddData("FakeDegree",["H4","E7","E8","G31"],
  t->function(phi,q)local f,res;f:=CHEVIE.R("cycpolfakedegrees",t)
   [Position(CHEVIE.R("CharInfo",t)().charparams,phi)];
   if IsList(f[1]) then res:=ValuePol(f[1],q^2);else res:=f[1];fi;
   f:=ShallowCopy(f);f[1]:=1;
   return res*Value(CycPol(f),q);
end);
  
CHEVIE.IndirectAddData("HighestPowerFakeDegrees",["H4","E7","E8","G31"],
t->function()return List(CHEVIE.R("cycpolfakedegrees",t),
  function(f)local res;
  if IsList(f[1])then res:=2*Length(f[1])+f[2]-2;else res:=f[2];fi;
  return res+Sum(f{[3..Length(f)]},Phi);end);
end);

CHEVIE.IndirectAddData("HighestPowerFakeDegrees",
["E6","G32","G33","G34","G2","F4","H3","G24","G25","G26","G27","G29"],
t->function()return List(CHEVIE.R("sparseFakeDegrees",t),x->x[Length(x)]);end);

CHEVIE.IndirectAddData("LowestPowerFakeDegrees",
["G2","F4","H3","H4","G24","G25","G26","G27","G29","E6","E7","E8",
 "G31","G32","G33","G34"], 
t->function()return List(CHEVIE.R("sparseFakeDegrees",t),x->x[2]);end);

CHEVIE.IndirectAddData("PrintDiagram",["E6","E7","E8"],
  t->function(indices,title)local i,r,digits,l;digits:="678";
  Print(title," ");
  r:=Position(digits,t[2])+5;
  l:=Length(String(indices[1]))+Length(String(indices[3]))+4;
  Print(String("",l-1),indices[2],"\n");
  Print(String("",Length(title)+l),"|\n");
  Print(SPrint(String("",Length(title)-2),indices[1]));
  for i in [3..r] do Print(" - ",indices[i]);od;
  Print("\n");
end);

CHEVIE.IndirectAddData("PrintDiagram",["H3","H4"],t->function(indices,title)
  local i;
  Print(title," ");
  Print(SPrint(String("",Length(String(indices[1]))-1),"5 \n"));
  Print(String("",Length(title)-1),indices[1]," - ",indices[2]," - ",indices[3]);
  if t="H4" then Print(" - ",indices[4]);fi;
  Print("\n");
end);

CHEVIE.IndirectAddData("HighestPowerGenericDegrees",
 ["G24","G27","G29","G33","G34","H3","H4","E6","E7","E8"],
 t->function()local N;
  N:=Sum(CHEVIE.R("ReflectionDegrees",t),x->x-1);
  return List(CHEVIE.R("CycPolSchurElements",t),
                          x->N-Degree(CycPol(x)));end);

CHEVIE.IndirectAddData("LowestPowerGenericDegrees",
["G24","G27","G29","G33","G34","H3","H4","E6","E7","E8"],
 t->function()return List(CHEVIE.R("CycPolSchurElements",t),x->-x[2]);
end);

CHEVIE.IndirectAddData("DecompositionMatrix",["F4","G2","G25","G26"],
  t->function(p)local T, m; 
  T:=CHEVIE.R("CharTable",t)(); T.name:=T.identifier; 
  m:=DecompositionMatrix(T mod p);
  return List(BlocksMat(m),c->[c[1],m{c[1]}{c[2]}]);
  end);

CHEVIE.IndirectAddData("SchurElement",["G24","G27","G29","G33","G34",
 "E6","E7","E8","H3","H4"],
t->function(arg)return Value(CycPol(CHEVIE.R("CycPolSchurElements",t)
  [Position(CHEVIE.R("CharInfo",t)().charparams,arg[1])]),
   -arg[2][1][1]/arg[2][1][2]);
end);

CHEVIE.IndirectAddData("FactorizedSchurElement",["G24","G27","G29","G33",
 "G34","E6","E7","E8","H3","H4"], t->function(arg)local c,q,res,v,e;
#arg= [phi,q] for G24--G34, [phi,q,rootparam] for E6--H4
  c:=CHEVIE.R("CycPolSchurElements",t)
      [Position(CHEVIE.R("CharInfo",t)().charparams,arg[1])];
  q:=-arg[2][1][1]/arg[2][1][2];
  res:=rec(factor:=Mvp(c[1]*q^c[2]),operations:=FactorizedSchurElementsOps);
  res.vcyc:=List(c{[3..Length(c)]},v->rec(monomial:=q,pol:=CycPol([1,0,v])));
  return FactorizedSchurElementsOps.Simplify(res);
end);

CHEVIE.IndirectAddData("FactorizedSchurElement",["G2","F4","G25","G26","G32"],
  t->function(arg)local Y,ci;
  Y:=Concatenation(arg[2]{CHEVIE.R("HyperplaneRepresentatives",t)});
  ci:=CHEVIE.R("SchurData",t)[
    Position(CHEVIE.R("CharInfo",t)().charparams,arg[1])];
  return ApplyFunc(VFactorSchurElement,
    Concatenation([Y,CHEVIE.R("SchurModels",t).(ci.name),ci],
    arg{[3..Length(arg)]}));
end);

CHEVIE.IndirectAddData("SchurElement",["F4","G25","G26","G32"],
  t->function(arg)local Y,ci;
  Y:=Concatenation(arg[2]{CHEVIE.R("HyperplaneRepresentatives",t)});
  ci:=CHEVIE.R("SchurData",t)[
    Position(CHEVIE.R("CharInfo",t)().charparams,arg[1])];
  return VcycSchurElement(Y,CHEVIE.R("SchurModels",t).(ci.name),ci);
end);

############################################################################
#  VcycSchurElement(Y,schur model[,schur data])
#
#  This function computes the Schur elements for G4-22,  G25-26, G28, G32
#  according to the data computed by M. Chlouveraki.
#  Y is the list of parameters of the algebra.
#  schur model describes the shape of the Schur element: it has the fields
#   .factor=(possibly fractional) vecmonomial
#   .coeff= a constant
#   [nothing] or [.root=vecmonomial] or [.rootUnity]  
#   vcyc= a list of pairs [vecmonomial, cyclotomic polynomial index]
#   rootCoeff=  a constant by which multiply .root before taking root
#  vecmonomial=vector of powers for elts of Y (plus possibly
#     the power to which to raise root or rootUnity)
#  schur data describes the Schur element in its Galois orbit : it has fields
#   order: in which order to take the variables
#   rootPower: by which E(root)^i multiply .root
VcycSchurElement:=function(arg)local r,data,i,para,res,n,monomial,den,root;
  n:=Length(arg[1]);
  if Length(arg)=3 then data:=arg[3];para:=arg[1]{data.order};
                   else para:=ShallowCopy(arg[1]);fi;
  monomial:=v->Product([1..Length(v)],i->para[i]^v[i]);
  r:=arg[2];
  if IsBound(r.coeff) then res:=r.coeff;else res:=1;fi;
  if IsBound(r.factor) then res:=res*monomial(r.factor);fi;
  if IsBound(r.root) then
    para:=para+0*Product(para);para[n+1]:=ChevieIndeterminate(para);
  elif IsBound(r.rootUnity) then para[n+1]:=r.rootUnity^data.rootUnityPower;fi;
  res:=res*Product(r.vcyc,
    x->Value(CyclotomicPolynomial(Cyclotomics,x[2]),monomial(x[1])));
  if IsBound(r.root) then
    den:=Lcm(List(r.root,Denominator));
    root:=monomial(den*r.root);
    if IsBound(r.rootCoeff) then root:=root*r.rootCoeff;fi;
    return EvalPolRoot(res,root,den,data.rootPower);
  else return res;
  fi;
end;

############################################################################
#  FactorizedSchurElement(parameters Y, schurModel r [,schur data])
#
#  This function computes the Schur elements for G4-22, G25-26, G28, G32
#  according to the data computed by M. Chlouveraki.
#
VFactorSchurElement:=function(arg)local para,r,data,res,n,monomial,den,root;
  n:=Length(arg[1]);
  if Length(arg)>=3 then data:=arg[3];para:=arg[1]{data.order};
                    else para:=ShallowCopy(arg[1]);fi;
  monomial:=v->Product([1..Length(v)],i->para[i]^v[i]);
  r:=arg[2]; res:=rec();
  if IsBound(r.coeff) then res.factor:=r.coeff;else res.factor:=1;fi;
  if IsBound(r.factor) then res.factor:=res.factor*monomial(r.factor);fi;
  if IsBound(r.root) then
    den:=Lcm(List(r.root,Denominator));
    root:=monomial(r.root*den);
    if IsBound(r.rootCoeff) then root:=root*r.rootCoeff;fi;
    para[n+1]:=GetRoot(root,den);
    if IsBound(data) then para[n+1]:=para[n+1]*data.rootPower;fi;
#   Print("root=",r.root,"\n");
#   Print(den,"-th root.",data[n+2]," of ","f=",f,"=>",para[n+1],"\n");
  elif IsBound(r.rootUnity) then para[n+1]:=r.rootUnity^data.rootUnityPower;
  fi;
  res.vcyc:=List(r.vcyc,
     v->rec(monomial:=monomial(v[1]),pol:=CycPol([1,0,v[2]])));
  if res.factor=0 or res.vcyc=[] then return res.factor;fi;
  res.operations:=FactorizedSchurElementsOps;
# res.factor:=Mvp(res.factor);
  return FactorizedSchurElementsOps.Simplify(res);
end;
