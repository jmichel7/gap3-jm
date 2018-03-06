#############################################################################
##
#A  tbl/cmpxtimp.g        CHEVIE library        Jean Michel
##
#Y  Copyright (C) 2011  The CHEVIE Team
##
##  Separated by JM on 7-8-2011.
##  This file contains data about imprimitive spetses attached to G333
##  and 3G422
##
CHEVIE.AddData("PhiFactors","timp",function(p,q,r,phi)local res,o;
  o:=OrderPerm(phi);
  if p=q then
    if p mod o=0 then
      if phi=(1,2,4) then return [1,1,E(3)^2];
      else res:=[1..r]*0+1;res[Length(res)]:=E(o); return res;
      fi;
    elif [p,q,r,o]=[3,3,3,4] then return [E(4),1,-E(4)];
    fi;
  elif p=2*q then return [1,E(3)];
  fi;
  Error("wrong arguments");
end);

CHEVIE.AddData("ReducedInRightCoset","timp",function(W,phi)
  local e,y,sets,sets2,i,v;
  # Case of G(3,3,3)
  # quads of roots which have the same CartanMat and are representatives of 
  # W-orbits of quads of reflections satisfying Corran-Picantin relations
  sets:=[[1,2,3,44],[21,3,1,32],[3,11,2,36],[22,3,2,16]]; # 3G333
  sets2:=[[1,50,3,12,2], [3,52,2,23,11], [1,16,3,43,38],  # 4G333
          [2,37,3,15,14], [50,3,52,38,53], [1,23,3,22,45]];
  for i in sets do
    y:=Set(List(i,j->Reflection(W,j)));
    e:=RepresentativeOperation(W,y,OnSets(y,phi),OnSets);
    if e<>false then 
      return rec(gen:=i{[1..3]},phi:=phi/e);
    fi;
  od;
  for i in sets2 do
    y:=W.rootInclusion{i{[1..4]}};
    v:=Set(y);
    e:=RepresentativeOperation(W,v,OnSets(v,phi),OnSets);
    if e<>false and PermListList(y,OnTuples(y,phi/e))=(1,2,3,4) then 
      return rec(gen:=i{[1,5,3]},phi:=phi/e);
    fi;
  od;
  return false;
  end
);

CHEVIE.AddData("ClassInfo","timp",function(p,q,r,phi)
  if [p,q,r]=[3,3,3] then
    if phi=(1,4,2) then
      return rec( classtext:=[[],[3],[1],[2,1],[3,1],[3,2,1],[2,3,1,2,3,1]],
        classparams:=[[],[3],[1],[2,1],[3,1],[3,2,1],[2,3,1,2,3,1]],
        classnames:=["Id","3","1","21","31","321","231231"]);
    elif phi=(1,2,4) then
      return rec( classtext:=[[],[3],[1],[1,2],[3,1],[3,1,2],[1,3,2,1,3,2]],
        classparams:=[[],[3],[1],[1,2],[3,1],[3,1,2],[1,3,2,1,3,2]],
        classnames:=["Id","3","1","12","31","312","132132"]);
    elif phi=(1,2,3,4) then
      return rec( classtext:=[[],[1],[1,2],[1,2,3],[1,2,1],[1,2,1,3]],
        classparams:=[[],[1],[1,2],[1,2,3],[1,2,1],[1,2,1,3]],
        classnames:=["Id","1","12","123","121","1213"]);
    else Error("should not happen");
    fi;
  elif [p,q,r]=[4,2,2] then
    return rec(classtext:=[[],[1],[1,2,3,1,2,3],[1,2,3,1,2,3,1,2,3]],
    classes:=[4,4,4,4],
    classnames:=["Id","1","cc","z"]);
  else ChevieErr("ClassInfo not implemented");return false;
  fi;
end);

CHEVIE.AddData("NrConjugacyClasses","timp",function(p,q,r,phi)
  return Length(CHEVIE.R("ClassInfo","timp")(p,q,r,phi).classtext);
end);

CHEVIE.AddData("CharInfo","timp",function(p,q,r,phi)
  if [p,q,r]=[3,3,3] then
    if phi=(1,4,2) or phi=(1,2,4) then
      return rec(charparams:=[[[],[],[3]],[[],[],[1,1,1]],[[],[],[2,1]],[[],[1,1],[1]],
       [[],[2],[1]],[[],[1],[1,1]],[[],[1],[2]]],
       extRefl:=[1,5,6,2]);
    elif phi=(1,2,3,4) then return rec(charparams:=[[[],[],[3]],[[],[],[1,1,1]],
     [[],[1,1],[1]],[[],[2],[1]],[[],[1],[1,1]],[[],[1],[2]]],extRefl:=[1,4,5,2]);
    else Error("should not happen");
    fi;
  elif [p,q,r]=[4,2,2] then return rec(charparams:=[[[],[],[2],[]],[[],[],[],[1,1]],
    [[],[1],[1],[]],[[],[],[1],[1]]],extRefl:=[1,4,2]);
  else ChevieErr("CharInfo not implemented");return false;
  fi;
end);

CHEVIE.AddData("CharName","timp",
  function(arg)return ApplyFunc(CHEVIE.imp.CharName,Drop(arg,4));end);

CHEVIE.AddData("CharTable","timp",function(p,q,r,phi)local res;
  if [p,q,r]=[3,3,3] then
    if phi=(1,4,2) then res:=
    rec(size:=54, order:=54, centralizers:=[ 18, 6, 6, 18, 3, 6, 18 ],
    identifier:="3'G(3,3,3)", name:="3'G(3,3,3)",
    classes:=[ 3, 9, 9, 3, 18, 9, 3 ],
    irreducibles:=[[1,1,1,1,1,1,1],[1,-1,-1,1,1,-1,1],[2,0,0,2,-1,0,2],
     E(3)*[-ER(-3),            -1,-E(3)^2,2*E(3)+E(3)^2,0,-E(3),-E(3)-2*E(3)^2],
     E(3)*[-ER(-3),             1, E(3)^2,2*E(3)+E(3)^2,0, E(3),-E(3)-2*E(3)^2],
      [-2*E(3)-E(3)^2,-E(3)^2,     -1, ER(-3),      0,-E(3), E(3)+2*E(3)^2],
      [-2*E(3)-E(3)^2, E(3)^2,      1, ER(-3),      0, E(3), E(3)+2*E(3)^2]]);
    elif phi=(1,2,4) then res:=
    rec(size:=54, order:=54, centralizers:=[ 18, 6, 6, 18, 3, 6, 18 ],
    identifier:="3G(3,3,3)", name:="3G(3,3,3)",
    classes:=[ 3, 9, 9, 3, 18, 9, 3 ],
    irreducibles:=[[1,1,1,1,1,1,1],[1,-1,-1,1,1,-1,1],[2,0,0,2,-1,0,2],
      E(3)^2*[ER(-3),-1,-E(3),(-3-ER(-3))/2,0,-E(3)^2,(3-ER(-3))/2],
      E(3)^2*[ER(-3),1,E(3),(-3-ER(-3))/2,0,E(3)^2,(3-ER(-3))/2],
      [(3+ER(-3))/2,-E(3),-1,-ER(-3),0,-E(3)^2,(-3+ER(-3))/2],
      [(3+ER(-3))/2,E(3),1,-ER(-3),0,E(3)^2,(-3+ER(-3))/2]]);
    elif phi=(1,2,3,4) then res:=
      rec(size:=54, order:=54, centralizers:=[6,6,6,6,6,6],
      identifier:="4G(3,3,3)", name:="4G(3,3,3)",
      classes:=[ 9,9,9,9,9,9],
      irreducibles:=[[1,1,1,1,1,1],[1,-1,1,-1,-1,1],
       [1,E(3),E(3)^2,1,E(3)^2,E(3)],[1,-E(3),E(3)^2,-1,-E(3)^2,E(3)],
       [1,E(3)^2,E(3),1,E(3),E(3)^2],[1,-E(3)^2,E(3),-1,-E(3),E(3)^2]]);
    else Error("should not happen");
    fi;
  elif [p,q,r]=[4,2,2] then 
    res:=rec(size:=16, order:=16, centralizers:=[ 4, 4, 4, 4 ],
    classes:=[ 4, 4, 4, 4 ],
    identifier:="3G(4,2,2)", name:="3G(4,2,2)",
    irreducibles:=[[1,1,1,1],[1,-1,1,-1],[-1,E(4),1,-E(4)],[-1,-E(4),1,E(4)]]);
  else ChevieErr("CharTable not implemented");return false;
  fi;
  res.text:="origin: Dixon's Algorithm";
  return CHEVIE.compat.MakeCharacterTable(res);
end);

CHEVIE.AddData("UnipotentCharacters","timp",function(p,q,r,phi)local res,a;
  if [p,q,r]=[3,3,3] then
    if phi=(1,4,2) then 
  return rec(
  harishChandra:=[
 #  rec(relativeType:=rec(series:="ST",indices:=[1,2,44,3],p:=3,q:=1,rank:=2),
    rec(relativeType:=rec(series:="ST",indices:=[1,2],p:=3,q:=1,rank:=2),
      levi:=[  ], eigenvalue:=1,
      parameterExponents:=[[2,0,1],1], cuspidalName:="",
      charNumbers:=[7,3,5,2,4,9,1,6,8])],
  almostHarishChandra:=[
    rec(relativeType:=rec(orbit:=[rec(series:="ST",indices:=[1..3],
     p:=3,q:=3,rank:=3)],
      twist:=(1,2,4)),
      levi:=[  ], eigenvalue:=1, cuspidalName:="",
      charNumbers:=[1..7]),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0), levi:=[1..3],
	charNumbers:=[8], eigenvalue:=E(3),cuspidalName:="G_{3,3,3}[\\zeta_3]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0), levi:=[1..3],
	charNumbers:=[9], eigenvalue:=E(3)^2, 
	cuspidalName:="G_{3,3,3}[\\zeta_3^2]")],
  families:=[
    Family("C1",[1]),
    Family("C1",[2]),
    Family("C1",[3]),
    Family(ComplexConjugate(CHEVIE.families.X(3)),[7,5,8]),
    Family(ComplexConjugate(CHEVIE.families.X(3)),[6,4,9])],
  a:=[0,9,3,4,1,4,1,1,4],
  A:=[0,9,6,8,5,8,5,5,8]);
    elif phi=(1,2,4) then 
  return rec(
  harishChandra:=[
  # rec(relativeType:=rec(series:="ST",indices:=[1,2,44,3],p:=3,q:=1,rank:=2),
    rec(relativeType:=rec(series:="ST",indices:=[1,2],p:=3,q:=1,rank:=2),
      levi:=[  ], eigenvalue:=1,
      parameterExponents:=[[2,1,0],1], cuspidalName:="",
      charNumbers:=[7,5,3,9,4,2,1,8,6])],
  almostHarishChandra:=[
    rec(relativeType:=rec(orbit:=[rec(series:="ST",indices:=[1..3],
     p:=3,q:=3,rank:=3)],
      twist:=(1,2,4)),
      levi:=[  ], eigenvalue:=1, cuspidalName:="",
      charNumbers:=[1..7]),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0), levi:=[1..3],
	charNumbers:=[8], eigenvalue:=E(3), 
	cuspidalName:="G_{3,3,3}[\\zeta_3]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0), levi:=[1..3],
	charNumbers:=[9], eigenvalue:=E(3)^2, 
	cuspidalName:="G_{3,3,3}[\\zeta_3^2]")],
  families:=[
    Family("C1",[1]),
    Family("C1",[2]),
    Family("C1",[3]),
    Family(CHEVIE.families.X(3),[7,5,8]),
    Family(CHEVIE.families.X(3),[6,4,9])],
  a:=[0,9,3,4,1,4,1,1,4],
  A:=[0,9,6,8,5,8,5,5,8]);
    elif phi=(1,2,3,4) then 
  res:=rec(
  harishChandra:=[
    rec(relativeType:=rec(series:="ST",indices:=[1],p:=6,q:=1,rank:=1),
      levi:=[  ], eigenvalue:=1,
      parameterExponents:=[[3,1,2,0,2,1]], cuspidalName:="",
      charNumbers:=[1,5,4,2,6,3]),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0), levi:=[1..3],
	parameterExponents:=[],charNumbers:=[7], eigenvalue:=E(3)^2,
	cuspidalName:="{}^4G_{3,3,3}[\\zeta_3^2]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0), levi:=[1..3],
	parameterExponents:=[],charNumbers:=[8], eigenvalue:=E(3),
	cuspidalName:="{}^4G_{3,3,3}[\\zeta_3]")],
  almostHarishChandra:=[
    rec(relativeType:=rec(orbit:=[rec(series:="ST",indices:=[1..3],
     p:=3,q:=3,rank:=3)],
      twist:=(1,2,3,4)),
      levi:=[  ], eigenvalue:=1, cuspidalName:="",
      charNumbers:=[1..6]),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0), levi:=[1..3],
	charNumbers:=[7], eigenvalue:=E(3), 
	cuspidalName:="G_{3,3,3}[\\zeta_3]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0), levi:=[1..3],
	charNumbers:=[8], eigenvalue:=E(3)^2, 
	cuspidalName:="G_{3,3,3}[\\zeta_3^2]")],
  families:=[
    Family("C1",[1]),
    Family("C1",[2]),
    Family(ComplexConjugate(CHEVIE.families.X(3)),[3,5,7],rec(signs:=[1,1,-1])),
    Family(CHEVIE.families.X(3),[4,6,8],rec(signs:=[1,1,-1]))],
  a:=[0,9,4,1,4,1,4,1],
  A:=[0,9,8,5,8,5,8,5]);
    res.families[3].eigenvalues[3]:=E(3)^2;
    res.families[4].eigenvalues[3]:=E(3);
    a:=1;
    res.families[3].fourierMat[3]:=a*res.families[3].fourierMat[3];
    res.families[4].fourierMat[3]:=GaloisCyc(a,-1)*res.families[4].fourierMat[3];
    return res;
    else Error("should not happen");return false;
    fi;
  else ChevieErr("UnipotentCharacters not implemented");return false;
  fi;
end);
