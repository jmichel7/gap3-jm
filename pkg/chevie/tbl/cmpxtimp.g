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
      elif phi^-1=(1,2,4) then return [1,1,E(3)];
      else res:=List([1..r],x->E(1));res[Length(res)]:=E(o); return res;
      fi;
    elif [p,q,r,o]=[3,3,3,4] then return [E(4),1,-E(4)];
    fi;
  elif p=2*q then return [1,E(3)];
  fi;
  Error("wrong arguments");
end);

#Group( ( 1, 2,44)( 4, 8,49)( 5,48,24)( 6, 9,26)( 7,25,47)(10,11,50)(14,17,38)
#(15,36,43)(16,37,32)(18,19,39)(21,22,52)(41,42,54), ( 2,12,37)( 3,10,15)
#( 5,26,40)( 6,24,27)( 7,42,13)( 8,33,19)(14,21,46)(16,53,44)(18,51,49)
#(22,38,31)(28,54,47)(30,36,50)); # a section G of the Diagram automorphisms

# here W should have type[1] be G333 with indices [1,2,3]
# the representatives chosen are for the 24 elements of G: 12 in G and 12 in -G
CHEVIE.AddData("ReducedInRightCoset","timp",function(W,phi)
  local e,sets,sets2,i,m1perm,o,g;
  # Case of G(3,3,3)
  m1perm:=(1,4)(2,8)(3,13)(5,22)(6,14)(7,10)(9,17)(11,25)(12,33)(15,42)(16,18)
  (19,37)(20,35)(21,24)(23,29)(26,38)(27,46)(28,30)(31,40)(32,39)(34,45)(36,54)
  (41,43)(44,49)(47,50)(48,52)(51,53); # -1
  m1perm:=MappingPermListList(W.rootInclusion,Permuted(W.rootInclusion,m1perm));
  if phi in W then return rec(gen:=W.rootInclusion{[1..3]},phi:=());fi;
  g:=PermMatX(W,MatXPerm(W,phi)*MatXPerm(W,m1perm));
  if g in W then return rec(gen:=W.rootInclusion{[1..3]},phi:=g^-1*phi);fi;
  # quads of roots which have the same CartanMat and are representatives of 
  # <-1,W>-orbits of quads of reflections satisfying Corran-Picantin relations
  sets:=[[1,2,3,44],[2,12,11,37],[3,11,2,36],[1,12,10,16]]; # 3G333
  for g in [(),m1perm] do
    for i in List(sets,x->W.rootInclusion{x}) do
      for o in [[4,1,3,2],[2,4,3,1]] do
        e:=RepresentativeOperation(W,i{o},OnTuples(i,g*phi),OnTuples);
        if e<>false then return rec(gen:=i{[1..3]},phi:=phi/e); fi;
      od;
    od;
  od;
  # sextuples of roots which have the same CartanMat and are representatives of 
  # W-orbits of each element of order 4 of G
  sets2:=[[1,2,32,16,3,36,30,10],[3,10,30,36,17,21,53,38],
          [2,12,16,53,11,10,43,36],[2,44,16,37,3,43,30,11],
          [1,37,32,44,15,30,50,3],[1,12,32,53,10,50,36,15]]; # 4G333
  for i in List(sets2,x->W.rootInclusion{x}) do
    e:=RepresentativeOperation(W,i{[2,3,4,1]},OnTuples(i{[1..4]},phi),OnTuples);
    if e<>false then return rec(gen:=i{[1,2,5]},phi:=phi/e);fi;
  od;
  return false;
  end
);

CHEVIE.AddData("ClassInfo","timp",function(p,q,r,phi)
  if [p,q,r]=[3,3,3] then
    if phi=(1,4,2) then
      return rec(classes:=[3,9,9,3,18,9,3],
        classtext:=[[],[3],[1],[2,1],[3,1],[3,2,1],[2,3,1,2,3,1]],
        classparams:=[[],[3],[1],[2,1],[3,1],[3,2,1],[2,3,1,2,3,1]],
        classnames:=["Id","3","1","21","31","321","231231"]);
    elif phi=(1,2,4) then
      return rec(classes:=[3,9,9,3,18,9,3],
        classtext:=[[],[3],[1],[1,2],[3,1],[3,1,2],[1,3,2,1,3,2]],
        classparams:=[[],[3],[1],[1,2],[3,1],[3,1,2],[1,3,2,1,3,2]],
        classnames:=["Id","3","1","12","31","312","132132"]);
    elif OrderPerm(phi)=4 then
      return rec(classes:=[9,9,9,9,9,9],
        classtext:=[[],[1],[1,2],[1,2,3],[1,2,1],[1,2,1,3]],
        classparams:=[[],[1],[1,2],[1,2,3],[1,2,1],[1,2,1,3]],
        classnames:=["Id","1","12","123","121","1213"]);
    elif OrderPerm(phi)=2 then return CHEVIE.imp.ClassInfo(3,3,3);
    else Error("phi=",phi," not implemented");
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

CHEVIE.AddData("CharInfo","timp",function(p,q,r,phi)local res;
  if [p,q,r]=[3,3,3] then
    if OrderPerm(phi)=3 then
      res:=rec(charparams:=[[[],[],[3]],[[],[],[1,1,1]],[[],[],[2,1]],
       [[],[1,1],[1]],[[],[2],[1]],[[],[1],[1,1]],[[],[1],[2]]],
       extRefl:=[1,5,6,2]);
    elif OrderPerm(phi)=4 then res:=rec(charparams:=[[[],[],[3]],[[],[],[1,1,1]],
   [[],[1,1],[1]],[[],[2],[1]],[[],[1],[1,1]],[[],[1],[2]]],extRefl:=[1,4,5,2]);
    else Error("phi=",phi," not implemented");
    fi;
  elif [p,q,r]=[4,2,2] then res:=rec(charparams:=[[[],[],[2],[]],[[],[],[],
    [1,1]],[[],[1],[1],[]],[[],[],[1],[1]]],extRefl:=[1,4,2]);
  else ChevieErr("CharInfo not implemented");return false;
  fi;
  res.charnames:=List(res.charparams,PartitionTupleToString);
  return res;
end);

CHEVIE.AddData("CharTable","timp",function(p,q,r,phi)local res;
  if [p,q,r]=[3,3,3] then
    if phi=(1,4,2) then res:=
    rec(size:=54, order:=54, centralizers:=[ 18, 6, 6, 18, 3, 6, 18 ],
    identifier:="3'G(3,3,3)", name:="3'G(3,3,3)",
    classes:=[ 3, 9, 9, 3, 18, 9, 3 ],
    irreducibles:=[[1,1,1,1,1,1,1],[1,-1,-1,1,1,-1,1],[2,0,0,2,-1,0,2],
     E(3)*[-ER(-3),        -1,-E(3)^2,2*E(3)+E(3)^2,0,-E(3),-E(3)-2*E(3)^2],
     E(3)*[-ER(-3),         1, E(3)^2,2*E(3)+E(3)^2,0, E(3),-E(3)-2*E(3)^2],
      [-2*E(3)-E(3)^2,-E(3)^2, -1, ER(-3),      0,-E(3), E(3)+2*E(3)^2],
      [-2*E(3)-E(3)^2, E(3)^2,  1, ER(-3),      0, E(3), E(3)+2*E(3)^2]]);
    elif phi=(1,2,4) then res:=
    rec(size:=54, order:=54, centralizers:=[ 18, 6, 6, 18, 3, 6, 18 ],
    identifier:="3G(3,3,3)", name:="3G(3,3,3)",
    classes:=[ 3, 9, 9, 3, 18, 9, 3 ],
    irreducibles:=[[1,1,1,1,1,1,1],[1,-1,-1,1,1,-1,1],[2,0,0,2,-1,0,2],
      E(3)^2*[ER(-3),-1,-E(3),(-3-ER(-3))/2,0,-E(3)^2,(3-ER(-3))/2],
      E(3)^2*[ER(-3),1,E(3),(-3-ER(-3))/2,0,E(3)^2,(3-ER(-3))/2],
      [(3+ER(-3))/2,-E(3),-1,-ER(-3),0,-E(3)^2,(-3+ER(-3))/2],
      [(3+ER(-3))/2,E(3),1,-ER(-3),0,E(3)^2,(-3+ER(-3))/2]]);
    elif OrderPerm(phi)=4 then res:=
      rec(size:=54, order:=54, centralizers:=[6,6,6,6,6,6],
      identifier:="4G(3,3,3)", name:="4G(3,3,3)",
      classes:=[ 9,9,9,9,9,9],
      irreducibles:=[[1,1,1,1,1,1],[1,-1,1,-1,-1,1],
       [1,E(3),E(3)^2,1,E(3)^2,E(3)],[1,-E(3),E(3)^2,-1,-E(3)^2,E(3)],
       [1,E(3)^2,E(3),1,E(3),E(3)^2],[1,-E(3)^2,E(3),-1,-E(3),E(3)^2]]);
    else Error("phi=",phi," not implemented");
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
    elif OrderPerm(phi)=4 then 
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
    else Error("phi=",phi," not implemented");
    fi;
  else 
    if q=1 or q=p then ChevieErr("UnipotentCharacters not implemented");fi;
    return false;
  fi;
end);
