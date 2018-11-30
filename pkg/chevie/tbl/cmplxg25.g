#############################################################################
##
#A  tbl/cmplxg25.g       CHEVIE library          Gunter Malle and Jean Michel
##
#Y  Copyright (C) 1998 - 2001  The CHEVIE Team
##
##  This file contains data about the complex reflection group
##  of type G25 in the Shephard-Todd classification.
##
CHEVIE.AddData("PrintDiagram","G25",function(indices,title)
  Print(title," ",indices[1],"(3)--(3)",indices[2],"--(3)",indices[3],"\n");
end);

CHEVIE.AddData("GeneratingRoots","G25",
 [[0,0,-1],(-(2*E(3)^2+1)/3)*[1,1,1],[0,1,0]]);

CHEVIE.AddData("EigenvaluesGeneratingReflections","G25",[1/3,1/3,1/3]);

CHEVIE.AddData("HyperplaneRepresentatives","G25",[1]);

CHEVIE.AddData("BraidRelations","G25",[[[1,2,1],[2,1,2]],
  [[1,3],[3,1]],[[2,3,2],[3,2,3]]]);

CHEVIE.AddData("Size","G25", 648);

CHEVIE.AddData("ReflectionDegrees","G25",[6,9,12]);

CHEVIE.AddData("NrConjugacyClasses","G25", 24);

CHEVIE.AddData("ParabolicRepresentatives","G25",# repr. of conj. classes
function(s)local t;t:=[[[]],[[1]],[[1,2],[1,3]],[[1..3]]];
  return t[s+1];end);

# Position in classes of G26
# [1,4,7,8,11,13,15,17,19,21,24,27,28,29,32,34,36,37,39,40,42,43,45,47]

CHEVIE.AddData("ClassNames","G25",
 [".","cc","31","3131","12231223","1223","d","dd","z","zz","2231223",
  "d1","1","131","3221223221","11","1122","12","12z","122312231223",
  "332112","212","c","cz"]);

CHEVIE.AddData("WordsClassRepresentatives","G25",
 List(CHEVIE.R("ClassNames","G25"),
  x->Replace(x,".",[],"1",[1],"2",[2],"3",[3],"c",[1,2,3],
    "z",[1,2,3,1,2,3,1,2,3,1,2,3],"d",[1,2,3,2])));

CHEVIE.AddData("PowerMaps","G25",
 [,[1,9,4,3,15,5,8,7,10,9,3,15,16,14,5,13,16,13,4,1,10,20,
  2,21],[1,20,1,1,1,20,9,10,1,1,20,20,1,1,1,1,20,20,20,20,20,22,22,22],,[1,21,
  4,3,15,12,8,7,10,9,19,6,16,14,5,13,18,17,11,20,2,22,24,23],,[1,2,3,4,5,6,7,8,
  9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],,,,[1,21,4,3,15,12,8,7,10,9,
  19,6,16,14,5,13,18,17,11,20,2,22,24,23]]);

CHEVIE.AddData("ClassInfo","G25",
  rec(classtext:=CHEVIE.R("WordsClassRepresentatives","G25"),
  classnames:=CHEVIE.R("ClassNames","G25"),
  classparams:=CHEVIE.R("ClassNames","G25"),
  orders:=[1,6,3,3,3,6,9,9,3,3,6,6,3,3,3,3,6,6,6,2,6,4,12,12],
  classes:=[1,9,12,12,12,36,72,72,1,1,36,36,12,24,12,12,36,36,36,9,9,54,54,54]
));

CHEVIE.AddData("CharInfo","G25",function()local res;
  res:=rec(charparams:=
 [[1,0],[1,24],[1,12],[2,15],[2,3],[2,9],[3,6],[3,5,2],
  [3,5,1],[3,17],[3,13,2],[3,1],[3,13,1],[6,8,2], [6,8,1],[6,2],
  [6,4,2],[6,10],[6,4,1],[8,3],[8,9],[8,6],[9,5],[9,7]],
# The labelling is determined as follows:
# phi{3,5}' is complexconjugate of phi{3,1}
# phi{3,13}' is complexconjugate of phi{3,17}
# phi{6,8}' is complexconjugate of phi{6,10}
# phi{6,4}' is complexconjugate of phi{6,2}
  extRefl:=[1,12,8,3]);
  res.b:=List(res.charparams,x->x[2]);
  return res;
end);

CHEVIE.AddData("HeckeCharTable","G25",function(para,root)
  local u,v,w,f10,f23,f31,f62,f83,f97,res,c;
  u:=para[1][1];  v:=para[1][2];  w:=para[1][3];c:=(u*v*w)^0;
  res:=rec( name:="H(G25)", identifier:="H(G25)",
   parameter:=para, size:=648,order:=648,
  dim:=3, degrees:=[6,9,12], reflclasses:=[13],
  powermap:=CHEVIE.R("PowerMaps","G25"),
  irredinfo:=CHEVIE.R("IrredInfo","G25"));
  f10:=y->List(res.classtext,w->y^Length(w));
  f23:=function(u,v,w)return [2,-2*(u*v)^3,u^2+v^2,u^4+v^4,(u*v)^2*(u^4+v^4),
  -u*v*(u^2+v^2),-u^2*v^2,-u^4*v^4,2*u^6*v^6,2*u^12*v^12,-v^3*u^3*(u+v),
  -v^2*u^2*(u+v),u+v,(u+v)*(u^2-u*v+v^2),u^4*v^4*(u^2+v^2),u^2+v^2,
  -u*v*(u^2+v^2),u*v,u^7*v^7,-u^3*v^3*(u^2+v^2)*(v^4-u^2*v^2+u^4),-2*u^3*v^3,
  0,0,0];end;
  f31:=function(u,v,w)return [3,-u^4*v^2,2*u*v+u^2,2*u^2*v^2+u^4,
  u^4*v^4+2*u^6*v^2,-u^2*v^2,0,0,3*u^8*v^4,3*u^16*v^8,-u^4*v^3,-u^4*v,2*u+v,
  u*v^2+u^2*v+u^3,2*u^6*v^4+u^8*v^2,2*u^2+v^2,-u*v^3-u^3*v+u^4,u*v+u^2,
  u^9*v^5+u^10*v^4,-u^6*v^6,u^2*v^4-2*u^5*v,u^3,u^2*v,u^10*v^5];end;
  f62:=function(u,v,w)return[6, 2*u^3*v^2*w, 2*u*v+2*u*w+u^2+v^2,
   2*u^2*v^2+2*u^2*w^2+u^4+v^4, u^2*(v^4*w^2+2*u^2*v^2*w^2+2*v^4*u^2+u^4*w^2),
   u*w*(u^2+v^2), 0, 0, 6*u^6*v^4*w^2, 6*u^12*v^8*w^4, u^3*v^2*w*(w+u),
   u^2*v^2*(w+u), 3*u+2*v+w, v^2*u+u*w^2+v*u^2+u^2*w+u^3+v^3, 
   v^2*u^4*(3*v^2*w^2+2*u^2*w^2+u^2*v^2), 3*u^2+2*v^2+w^2,
   -(u^2+v^2)*(u*v-w^2-u^2), u*(u+v), v^4*u^7*w^2*(u+v),
   u^3*w^3*(u^2+v^2)*(v^4-u^2*v^2+u^4), v*(-2*u^3*w^2+3*u^4*v+v^3*w^2),
   -u*(-u^2+v*w), 0, 0];end;
  f83:=function(u,v,w) return [8, 0, 2*(w+u)*(u+v), 2*(w^2+u^2)*(u^2+v^2),
   2*u^2*v*w*(w^2+u^2)*(u^2+v^2), 0,-u^2*v*w, -u^4*v^2*w^2, 8*u^6*v^3*w^3,
   8*u^12*v^6*w^6, 0, 0, 4*u+2*v+2*w, v^2*u+u*w^2+v*w^2+v*u^2+u^2*w+v^2*w+2*u^3,
   2*u^4*v*w*(2*v^2*w^2+u^2*w^2+u^2*v^2),4*u^2+2*v^2+2*w^2,
   -u*v^3-u*w^3+u^2*v^2+u^2*w^2+v^2*w^2-u^3*v-u^3*w+u^4, u*(u+v+w),
   v^3*u^7*w^3*(u+v+w), 0, -2*u^3*v^3-2*u^3*w^3+v^3*w^3+3*u^4*v*w,
   -u*(-u^2+v*w), 0, 0];end;
  f97:=function(u,v,w,J) return [9,-3*J^2*u^2*v^2*w^2,(u+v+w)^2,(u^2+w^2+v^2)^2,
   J*(u^2*v^2+u^2*w^2+v^2*w^2)^2, -J^2*(u^2*v^2+u^2*w^2+v^2*w^2), 0, 0,
   9*J*u^4*v^4*w^4, 9*J^2*u^8*v^8*w^8, -v^2*J^2*u^2*w^2*(u+v+w),
   -v*u*w*J^2*(u*v+u*w+v*w), 3*u+3*v+3*w, (u+v+w)*(u^2+w^2+v^2),
   3*J*u^2*v^2*w^2*(u^2*v^2+u^2*w^2+v^2*w^2), 3*u^2+3*v^2+3*w^2,
   -u*v^3-u*w^3-v*w^3+u^2*v^2+u^2*w^2+v^2*w^2-u^3*v-u^3*w-v^3*w, u*v+u*w+v*w,
   v^4*J*u^4*w^4*(u*v+u*w+v*w), -u^6*v^6-u^6*w^6-v^6*w^6,
   -J*u*v*w*(2*w^3+2*v^3-3*u*v*w+2*u^3),-u*v*w,-J*u*v*w, -J^2*u^5*v^5*w^5];end;
  Inherit(res,CHEVIE.R("ClassInfo","G25"));
  res.centralizers:=List(res.classes,x->res.order/x);
# Position in chars of G26
# [1,3,4,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47]
  res.irreducibles:=[f10(u),f10(w),f10(v),f23(v,w,u),f23(u,v,w),f23(u,w,v),
[3,3*u^2*v^2*w^2,u^2+v^2+w^2,u^4+v^4+w^4,u^4*v^4+u^4*w^4+v^4*w^4,
u^2*v^2+u^2*w^2+v^2*w^2,0,0,3*u^4*v^4*w^4,3*u^8*v^8*w^8,
u^2*v^2*w^3+u^2*v^3*w^2+u^3*v^2*w^2,u*v^2*w^2+u^2*v*w^2+u^2*v^2*w,u+v+w,
u^3+v^3+w^3,u^2*v^4*w^4+u^4*v^2*w^4+u^4*v^4*w^2,u^2+v^2+w^2,
u^2*v^2+u^2*w^2+v^2*w^2,0,0,u^6*v^6+u^6*w^6+v^6*w^6,3*u^2*v^2*w^2,-u*v*w,
-u*v*w,-u^5*v^5*w^5],
   f31(v,u,w),f31(u,w,v),f31(w,v,u),f31(w,u,v),f31(u,v,w),
   f31(v,w,u),f62(w,u,v),f62(v,w,u),f62(u,v,w),f62(v,u,w),
   f62(w,v,u),f62(u,w,v),f83(u,v,w),f83(w,v,u),f83(v,u,w),
   f97(u,v,w,E(3)^2),f97(u,v,w,E(3))]*c;
  res := CHEVIE.compat.MakeCharacterTable(res);
  return res;
end);

CHEVIE.AddData("CharTable","G25",function()
  return CHEVIE.R("HeckeCharTable","G25")([[1,E(3),E(3)^2]],[]);end);

CHEVIE.AddData("sparseFakeDegrees","G25",
[[1,0],[1,24],[1,12],[1,15,1,21],[1,3,1,9],[1,9,1,15],[1,6,1,12,1,18],[1,5,1,
8,1,11],[1,5,1,8,1,11],[1,17,1,20,1,23],[1,13,1,16,1,19],[1,1,1,4,1,7],[1,13,
1,16,1,19],[1,8,1,11,2,14,1,17,1,20],[1,8,1,11,2,14,1,17,1,20],[1,2,1,5,2,8,1,
11,1,14],[1,4,1,7,2,10,1,13,1,16],[1,10,1,13,2,16,1,19,1,22],[1,4,1,7,2,10,1,
13,1,16],[1,3,2,6,2,9,2,12,1,15],[1,9,2,12,2,15,2,18,1,21],[1,6,2,9,2,12,2,15,
1,18],[1,5,1,8,3,11,2,14,2,17],[2,7,2,10,3,13,1,16,1,19]]);

# The data below was computed by Maria Chlouveraki
CHEVIE.AddData("SchurModels","G25",rec(
f1_0:=rec(
 vcyc:=[[[1,-1,0],1],[[1,-1,0],1],[[1,0,-1],1],[[1,0,-1],1],[[1,-1,0],4],
 [[1,0,-1],4],[[3,-2,-1],1],[[3,-1,-2],1],[[2,-1,-1],3],[[2,-1,-1],2],
 [[1,-1,0],6],[[1,0,-1],6]]),
f2_3:=rec(
 vcyc:=[[[1,0,-1],1],[[1,0,-1],1],[[0,1,-1],1],[[0,1,-1],1],[[1,0,-1],2],
[[0,1,-1],2],[[1,-1,0],1],[[1,-1,0],1],[[1,1,-2],3],[[1,1,-2],2],[[-1,1,0],6]]),
f3_1:=rec(
 vcyc:=[[[1,-1,0],1],[[-1,1,0],1],[[1,0,-1],1],[[1,0,-1],1],[[1,0,-1],2],
  [[0,1,-1],1],[[1,1,-2],2],[[2,-1,-1],2],[[1,0,-1],6],[[1,-1,0],4],
  [[2,1,-3],1]]),
f3_6:=rec(
 vcyc:=[[[1,-1,0],1],[[1,-1,0],1],[[1,0,-1],1],[[1,0,-1],1],[[0,1,-1],1],
  [[0,-1,1],1],[[-1,-1,2],2],[[-1,2,-1],2],[[-2,1,1],2]]),
f6_2:=rec(
 vcyc:=[[[-1,1,0],1],[[1,0,-1],1],[[-1,0,1],1],[[0,1,-1],1],[[0,1,-1],1],
  [[1,0,-1],2],[[1,0,-1],6],[[1,-2,1],2],[[0,1,-1],2],[[3,-2,-1],1]]),
f8_3:=rec(
 vcyc:=[[[0,1,-1],1],[[0,-1,1],1],[[-1,0,1],1],[[-1,1,0],1],[[2,-3,1],1],
  [[2,-1,-1],3],[[2,1,-3],1]]),
f9_7:=rec(rootUnity:=E(3),
  vcyc:=[[[0,0,0,1],1],[[0,0,0,2],2],[[-1,1,0],6],[[1,0,-1],6],
  [[0,-1,1],6],[[2,-1,-1,1],1],[[-1,2,-1,1],1],[[-1,-1,2,1],1]])));

CHEVIE.AddData("SchurData","G25",[
  rec(name:="f1_0",order:=[1,2,3]), rec(name:="f1_0",order:=[3,2,1]), 
  rec(name:="f1_0",order:=[2,1,3]), rec(name:="f2_3",order:=[2,3,1]),
  rec(name:="f2_3",order:=[1,2,3]), rec(name:="f2_3",order:=[1,3,2]),
  rec(name:="f3_6",order:=[1,3,2]), rec(name:="f3_1",order:=[2,1,3]),
  rec(name:="f3_1",order:=[1,3,2]), rec(name:="f3_1",order:=[3,2,1]),
  rec(name:="f3_1",order:=[3,1,2]), rec(name:="f3_1",order:=[1,2,3]),
  rec(name:="f3_1",order:=[2,3,1]), rec(name:="f6_2",order:=[3,1,2]),
  rec(name:="f6_2",order:=[2,3,1]), rec(name:="f6_2",order:=[1,2,3]),
  rec(name:="f6_2",order:=[2,1,3]), rec(name:="f6_2",order:=[3,2,1]),
  rec(name:="f6_2",order:=[1,3,2]), rec(name:="f8_3",order:=[1,2,3]),
  rec(name:="f8_3",order:=[3,2,1]), rec(name:="f8_3",order:=[2,1,3]),
  rec(name:="f9_7",order:=[1,2,3],rootUnityPower:=1), 
  rec(name:="f9_7",order:=[1,2,3],rootUnityPower:=2)]);

CHEVIE.AddData("HeckeRepresentation","G25",function(para,root,i)
  local u,v,w,f1,f2,f31,f32,f6,f8,f9,rep;
  u:=para[1][1];v:=para[1][2];w:=para[1][3];
  f1:=u->[[[u]],[[u]],[[u]]];
  f2:=function(v,w) return
    WGraph2Representation([[[1,3],[2]],[[1,2,-1,v*w]]],[w,v]);end;
  f31:=function(u,v)return WGraph2Representation(
    [[[1],[2],[3]],[[1,2,u,-v],[2,3,-v,u]]],[u,v]);end;
  f32:=function(u,v,w) return WGraph2Representation(
   [[[[2],[]],[[],[1,2,3]],[[1,3],[]]],
   [[1,2,-1,u*w+v^2],[1,3,v,v],[2,3,-u*w-v^2,1]]],[u,v,w]);end;
  f6:=function(v,u,w)return WGraph2Representation(
  [[[[2],[]],[[],[1,2]],[[1],[]],[[],[2,3]],[[3],[]],[[],[1,3]]],
  [[1,2,-1,v*w+u^2],[1,3,u,u],[1,4,-1,v*w+u^2],[1,5,-u,-u],[1,6,w,0],
   [2,3,-v*w-u^2,1],[2,6,-u*w,1],[4,5,v*w+u^2,-1],[4,6,-u*w,1]]],[v,u,w]);end;
  f8:=function(u,w,v)return WGraph2Representation(
   [[[[2,3],[]],[[3],[1,2]],[[1,3],[]],[[2],[3]],[[1,3],[]],[[2],[1]],
    [[1],[2,3]],[[1,2],[]]],
   [[1,2,-u*v-w^2,1],[1,3,w,w],[1,4,v*w-w^2,0],[1,5,0,-1],[2,3,-1,u*v+w^2],
    [2,4,[1,0,3,w],-u],[2,5,0,-w],[2,6,-1,0],[3,6,[1,0,3,v-w],-u],
    [3,7,u*w+w^2,-1],[3,8,-w,-w],[4,5,-u,[1,v,3,0]],[4,7,0,v],
    [5,6,[1,0,3,1],-u*w],[5,7,-u,v-w],[5,8,0,v*w-w^2],[6,7,u*w,[1,-1,3,0]],
    [6,8,0,v-w],[7,8,-1,u*v+w^2]]], [u,w,v]);end;
  f9:=function(u,v,w,a) return WGraph2Representation(
[[[[2],[]],[[],[1,2,3]],[[1],[3]],[[1,3],[]],[[2],[1]],[[1],[2]],[[2],[3]],
[[3],[2]],[[3],[1]]],[[1,2,-1,u*w+v^2],[1,3,v,[1,v,3,0]],[1,4,-a*v,0],
[1,5,0,a^2*u-v],[1,6,0,a^2*u],[1,7,0,a^2*u-v],[1,8,0,-a^2*u],[1,9,v,[1,0,3,v]],
[2,3,-u*w-v^2,1],[2,4,-u*w+a*v^2,0],[2,5,-a^2*v*w,0],[2,7,-a^2*v*w,0],
[2,9,-u*w-v^2,1],[3,4,0,u+a^2*v],[3,5,0,u],[3,6,-a^2*w,a*v],[3,7,w,0],
[4,5,[1,0,3,-w],u],[4,6,-a*w,0],[4,7,[1,-w,3,0],u],[4,8,a*w,0],[4,9,u+a^2*v,0],
[5,6,-u,v],[5,9,0,w],[7,8,u,-v],[7,9,u,0],[8,9,-a*v,a^2*w]]],[u,v,w]);end;
  rep:=[[f1,u],[f1,w],[f1,v],[f2,v,w],[f2,u,v],[f2,u,w],[f32,u,v,w],
   [f31,u,v],[f31,w,u],[f31,v,w],[f31,u,w],[f31,v,u],[f31,w,v],
   [f6,v,u,w], [f6,u,w,v], [f6,w,v,u], [f6,w,u,v], [f6,u,v,w], [f6,v,w,u], 
   [f8,u,v,w],[f8,w,u,v],[f8,v,w,u],[f9,u,v,w,E(3)],[f9,u,v,w,E(3)^2]];
  return ApplyFunc(rep[i][1],rep[i]{[2..Length(rep[i])]})*Product(para[1])^0;
end);

CHEVIE.AddData("UnipotentCharacters","G25",function()local J;J:=E(3);
  return rec(
   harishChandra:=[
    rec(relativeType:=rec(series:="ST",indices:=[1..3],rank:=3,ST:=25),
 levi:=[], parameterExponents:=[1,1,1],
        charNumbers:=[1..24], eigenvalue:=1, cuspidalName:=""),
    rec(relativeType:=rec(series:="ST",indices:=[3,2],rank:=2,p:=3,q:=1),
        levi:=[1],parameterExponents:=[1,3],
	charNumbers:=[39,31,30,41,38,40,25,27,26], eigenvalue:=J^2,
	cuspidalName:=ImprimitiveCuspidalName([[],[0,1],[0,1]])),
    rec(relativeType:=rec(series:="ST",indices:=[2],rank:=1,p:=6,q:=1),
        levi:=[1,3],
        parameterExponents:=[[ 3, 3, 2, 0, 0, 2 ]],
	charNumbers:=[ 29, 28, 32, 44, 43, 33 ], eigenvalue:=J, 
	cuspidalName:=Concatenation(ImprimitiveCuspidalName([[],[0,1],[0,1]]),
	 "\\otimes ",ImprimitiveCuspidalName([[],[0,1],[0,1]]))),
    rec(relativeType:=rec(series:="ST",indices:=[3],rank:=1,p:=3,q:=1),
        levi:=[1..2],parameterExponents:=[[0,4,4]],
	charNumbers:=[42,34,35], eigenvalue:=-1,
        cuspidalName:="G_4"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
        levi:=[1..3], parameterExponents:=[],
	charNumbers:=[36], eigenvalue:=-J,
        cuspidalName:="G_{25}[-\\zeta_3]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
        levi:=[1..3], parameterExponents:=[],
	charNumbers:=[37], eigenvalue:=J,
        cuspidalName:="G_{25}[\\zeta_3]")],
   families:=[
     Family("C1",[1]),
     Family(CHEVIE.families.X(3),[12,9,25],rec(signs:=[1,1,-1],ennola:=-2)),
     Family(CHEVIE.families.QZ(3),[20,16,19,6,28,26,5,27,29],
       rec(signs:=[1,1,1,1,1,-1,1,1,1],special:=2,cospecial:=3,ennola:=5)),
     Family(CHEVIE.families.X(6),[17,23,7,24,14,32,34,30,36,8,37,31,11,35,33],
       rec(signs:=[1,1,1,1,1,1,-1,-1,1,-1,1,-1,1,1,-1],ennola:=-10)),
     Family(CHEVIE.families.X(3),[22,21,38],rec(signs:=[1,1,-1],ennola:=1)),
     Family(CHEVIE.families.X(3),[15,18,39],rec(signs:=[1,1,-1],ennola:=-2)),
     Family(SubFamilyij(CHEVIE.families.ExtPowCyclic(6,3),1,2,-ER(2)/ER(-1)),
      [3,13,40,10,41,2,43,42,4,44],rec(signs:=[1,1,1,1,-1,1,-1,1,-1,-1],
       cospecial:=6,ennola:=6))],
  a:=[0,12,12,12,2,2,4,4,1,12,4,1,12,4,8,2,4,8,2,2,6,6,4,4,1,2,2,2,2,4,4,4,4,
      4,4,4,4,6,8,12,12,12,12,12],
  A:=[0,24,24,24,16,16,20,20,11,24,20,11,24,20,22,16,20,22,16,16,21,21,20,20,
      11,16,16,16,16,20,20,20,20,20,20,20,20,21,22,24,24,24,24,24]);
end);

CHEVIE.AddData("Invariants","G25",
  [function(x1,x2,x3)
    return -10*x1^3*x2^3-10*x1^3*x3^3-10*x2^3*x3^3+x1^6+x2^6+x3^6;end,
   function(x1,x2,x3)
    return -x1^3*x2^6+x1^3*x3^6-x2^3*x3^6+x1^6*x2^3-x1^6*x3^3+x2^6*x3^3;end,
    function(x1,x2,x3) return 2*x1^3*x2^3*x3^6+2*x1^3*x2^6*x3^3+x1^3*x2^9
     +x1^3*x3^9+x2^3*x3^9+2*x1^6*x2^3*x3^3-4*x1^6*x2^6-4*x1^6*x3^6-4*x2^6*x3^6
     +x1^9*x2^3+x1^9*x3^3+x2^9*x3^3;end]
);

# the discriminant as a polynomial in the invariants
CHEVIE.AddData("Discriminant","G25",function()return function(t1,t2,t3)
  return  36*t1*t2^2*t3-t1^2*t3^2-32*t3^3+t1^3*t2^2+108*t2^4;
end;end);
