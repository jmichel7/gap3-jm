#############################################################################
##
#A  tbl/weyl2f4.g            CHEVIE library                      Meinolf Geck
##
#Y  Copyright (C) 1998 - 2001  The CHEVIE Team
##
##  This file contains information about the coset W.F where W is
##  the group F4:   1-2=3-42      and F is (1,4)(2,3).
##
##  The HeckeCharTable was obtained by explicit computation using
##  representations of the Hecke algebra of type F_4.
##
CHEVIE.AddData("NrConjugacyClasses", "2F4", 11);

CHEVIE.AddData("ClassInfo","2F4",function()local res;
  res:=rec(classtext:=[[],[2,3,2],[1],[1,2,1,3,2,1,4,3,2,1],[1,2],[2],
    [1,2,1,3,2,1,3,2],[1,2,3,2],[1,2,1,3,2,1,3,2,4,3,2,1],
    [1,2,1,3,2,1,3,2,4,3,2,1,3,2,4,3,2,1],[1,2,1,3,2,1]],
  orders:=[ 2, 8, 4, 24, 24, 8, 8, 12, 4, 8, 8],
  classes:=[ 72, 144, 288, 96, 96, 144, 72, 192, 24, 12, 12 ],
     classnames:=["2a","8a","4a","24a","24b","8a","8b","12a","4b","8c","8d"]);
  res.classparams:=List(res.classtext,IntListToString);
  return res;
end);

CHEVIE.AddData("CharInfo","2F4",function()local res,resparams;
  res:=rec(extRefl:=[1,9,7,10,2],
  charparams:=[[1,0],[1,24],[4,8],[9,2],[9,10],[6,6,1],[6,6,2],[12,4],
   [4,1],[4,13],[16,5]],
  kondo:=["1_1","1_4","4_1","9_1","9_4","6_1","6_2","12","4_2","4_5","16"]);
  resparams:=CHEVIE.R("CharInfo","F4")().charparams;
  res.charRestrictions:=List(res.charparams,x->Position(resparams,x));
  res.nrGroupClasses:=Length(resparams);
  res.b:=List(res.charparams,x->x[2]);
  return res;
end);

CHEVIE.AddData("cyclestructure","2F4",[ [ 24 ], [ 4,,,,,, 5 ], [ 2,, 11 ],
  [ ,,,,,,,,,,,,,,,,,,,,,, 2 ], [,,,,,,,,,,,,,,,,,,,,,, 2 ], [4,,,,,, 5 ],
  [ ,,,,,, 6 ], [,,,,,,,,,, 4 ], [ ,, 12 ], [ ,,,,,, 6 ], [ ,,,,,, 6 ] ]);

CHEVIE.AddData("generators","2F4",
  [ ( 1,25)( 2, 5)( 6, 9)( 7,11)(10,14)(12,15)(16,18)(21,23)(26,29)(30,33)
    (31,35)(34,38)(36,39)(40,42)(45,47), ( 1, 5)( 2,26)( 3, 7)( 8,12)( 9,13)
    (14,17)(18,20)(19,21)(25,29)(27,31)(32,36)(33,37)(38,41)(42,44)(43,45),
  ( 2, 6)( 3,27)( 4, 8)( 5, 9)(12,16)(15,18)(17,19)(20,22)(26,30)(28,32)
    (29,33)(36,40)(39,42)(41,43)(44,46), ( 3, 8)( 4,28)( 6,10)( 7,12)( 9,14)
    (11,15)(13,17)(22,24)(27,32)(30,34)(31,36)(33,38)(35,39)(37,41)(46,48) ]);

CHEVIE.AddData("phi","2F4",
( 1, 4)( 2, 3)( 5, 8)( 6, 7)( 9,12)(10,11)(13,16)(14,15)(17,18)(19,20)(21,22)
(23,24)(25,28)(26,27)(29,32)(30,31)(33,36)(34,35)(37,40)(38,39)(41,42)(43,44)
(45,46)(47,48));

CHEVIE.AddData("CartanMat","2F4",
  [[2,-1,0,0],[-1,2,-ER(2),0],[0,-ER(2),2,-1],[0,0,-1,2 ]]);

#############################################################################
##
#F  CHEVIE.R("ClassParameter","2F4")( <w> )  . .  class parameter of w
##
##  given an element w  of a Coxeter group W of type F_4 as word in standard
##  generators, returns the classparam of its F-conjugacy class under the
##  nontrivial F action permuting the generators 1<->4 2<->3.
##
CHEVIE.AddData("ClassParameter","2F4",function(w)local x,l;
  if Length(w)=1 then
    if CycleStructurePerm(x*CHEVIE.R("phi","2F4"))=[2,,11] then
         return CHEVIE.R("classparams","2F4")[3];
    else return CHEVIE.R("classparams","2F4")[6];
    fi;
  else
    l:=List(CHEVIE.R("WordsClassRepresentatives","2F4"),Length);
    return CHEVIE.R("classparams","2F4")[Position(l,Length(w))];
  fi;
end);

###########################################################################
##
#F  CHEVIE.R("HeckeCharTable","2F4")( <v> )  . . . . . . . . . .
##  . . character table of the Hecke algebra of type 2F4 with parameter v^2.
##
CHEVIE.AddData("HeckeCharTable","2F4",function(param,sqrtparam)local v,q,tbl;
 q:=-param[1][1]/param[1][2];
 if not IsBound(sqrtparam[1]) then v:=GetRoot(q,2,"CharTable(Hecke(2F4))");
 else v:=-sqrtparam[1]/param[1][2];
 fi;
 tbl:=rec(
  identifier:="H(2F4)",
  parameter:=[q,q,q,q],
  sqrtparameter:=[v,v,v,v],
  cartan:=CHEVIE.R("CartanMat","2F4"),
  size:=1152,
  irreducibles:=[
   [ 1, v^6, v^2, v^20, v^4, v^2, v^16, v^8, v^24, v^36, v^12 ],
   [ 1, -1, -1, 1, 1, -1, 1, 1, 1, 1, 1 ],
   [ 2,v^6-1,v^2-1,-v^10,-v^2,v^2-1,v^10+v^6,-v^4,2*v^12,2*v^18,2*v^6 ],
   [1,v^6,-1,0,0,v^2,v^12-2*v^10,0,-3*v^16,3*v^24,3*v^8 ],
   [-1,1,-v^2,0,0,1,2*v^6-v^4,0, 3*v^8, -3*v^12, -3*v^4 ],
   [0,0,0,v^10,v^2,0,2*v^8,-v^4,-4*v^12,-2*v^18,-2*v^6 ],
   [-2,-v^6+1,-v^2+1,v^10,v^2,-v^2+1,v^10-2*v^8+v^6,-v^4,2*v^12,4*v^18,4*v^6],
   [ 2,v^6-1,v^2-1,v^10,v^2,v^2-1,-v^10-v^6,-v^4,2*v^12,-2*v^18,-2*v^6 ],
   [ 0, -ER(2)*v^3,0,-ER(2)*v^15, ER(2)*v^3, ER(2)*v, ER(2)*(v^13-v^11), 0,0,
     -2*ER(2)*v^27, 2*ER(2)*v^9 ],
   [ 0, -ER(2)*v^3, 0, ER(2)*v^5, -ER(2)*v,ER(2)*v, ER(2)*(v^5-v^3), 0, 0,
       2*ER(2)*v^9, -2*ER(2)*v^3 ],
   [ 0, 0, 0, ER(2)*v^10, -ER(2)*v^2, 0, ER(2)*(v^10-2*v^8+v^6), 0,
       0, -4*ER(2)*v^18, 4*ER(2)*v^6 ] ]*v^0,
  irredinfo:=CHEVIE.R("IrredInfo","2F4"));
  Inherit(tbl,CHEVIE.R("ClassInfo","2F4")());
  tbl.centralizers:=List(tbl.classes,x->tbl.size/x);
  tbl := CHEVIE.compat.MakeCharacterTable(tbl);
  CHEVIE.compat.AdjustHeckeCharTable(tbl,param);
  return tbl;
end);

CHEVIE.AddData("PhiFactors","2F4",[1,-1,1,-1]);

CHEVIE.AddData("Representation","2F4",function(i)return
  CHEVIE.R("HeckeRepresentation","2F4")
  (List([1..4],x->[1,-1]),[1,1,1,1],i);end);

CHEVIE.AddData("HeckeRepresentation","2F4",function(para,rootpara,i)local v,res;
  if IsBound(rootpara[1]) then v:=rootpara[1]*para[1][2];
  else v:=GetRoot(-para[1][1]/para[1][2],2,
                  "Representation(Hecke(2F4),[",i,"])");
  fi;
  res:=[rec(gens:=1,F:=[[1]]), rec(gens:=4,F:=[[1]]),
   rec(gens:=9,F:=Permuted(IdentityMat(4),(1,4))),
   rec(gens:=10,F:=Permuted(IdentityMat(9),(1,9)(2,6)(4,8)(5,7))),
   rec(gens:=13,F:=-Permuted(IdentityMat(9),(1,9)(2,6)(4,8)(5,7))),
   rec(gens:=14,F:=Permuted(IdentityMat(6),(1,3)(2,5)(4,6))),
   rec(gens:=15,F:=-[[0,0,0,0,0,2],[0,0,0,0,1,0],[0,0,1,0,0,0],
		[0,0,0,1,0,0],[0,1,0,0,0,0],[1/2,0,0,0,0,0]]),
   rec(gens:=16,F:=Permuted(IdentityMat(12),(1,12)(2,5)(3,10)(4,9)(7,11))),
   rec(gens:=17,F:=-[[0,0,0,1],[0,0,1,0],[0,2,0,0],[2,0,0,0]]/ER(2)),
   rec(gens:=20,F:=[[0,0,0,2],[0,0,2,0],[0,1,0,0],[1,0,0,0]]/ER(2)),
   rec(gens:=25,
     F:=-[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
   [0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0],
   [0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0],[0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
   [0,0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0],
   [0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0],
   [0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0],[0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0],
   [0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0],[0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0],
   [0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0],[2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]/ER(2))];
  res:=res[i];
  res.gens:=-para[1][2]*WGraphToRepresentation(4,
                      CHEVIE.R("WGraph","F4")(res.gens),v)*v^0;
  return res;
end);

CHEVIE.families.X2:=rec(name:="X_2",
   fourierMat:=ER(2)/2*[[-1,-1],[-1,1]], eigenvalues:=[E(8)^3,-E(8)], 
   charLabels:=["1","2"], special:=1, sh:=[1,-1]);

CHEVIE.AddData("UnipotentCharacters","2F4",function()local res;
  return rec(
  harishChandra:=[
    rec(relativeType:=rec(series:="I",indices:=[1,2],bond:=8,rank:=2),
      parameterExponents:=[2,4],
      levi:=[], eigenvalue:=1,
      cuspidalName:="",
      charNumbers:=[1,5,4,2,3,20,19]),
    rec(relativeType:=rec(series:="A",indices:=[1],rank:=1),
      parameterExponents:=[12],
      levi:=[2,3], eigenvalue:=E(8)^3,
      cuspidalName:="{}^2B_2[1,3]",
      charNumbers:=[10,9]),
    rec(relativeType:=rec(series:="A",indices:=[1],rank:=1),
      parameterExponents:=[12],
      levi:=[2,3], eigenvalue:=-E(8),
      cuspidalName:="{}^2B_2[1,5]",
      charNumbers:=[13,12]),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
        levi:=[1..4], parameterExponents:=[],
  	charNumbers:=[6], eigenvalue:=-E(3),
        cuspidalName:="{}^2F_4[-\\zeta_3]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
        levi:=[1..4], parameterExponents:=[],
  	charNumbers:=[7], eigenvalue:=-1,
        cuspidalName:="{}^2F_4[-1]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
        levi:=[1..4], parameterExponents:=[],
  	charNumbers:=[8], eigenvalue:=-1,
        cuspidalName:="{}^2F_4^2[-1]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
        levi:=[1..4], parameterExponents:=[],
  	charNumbers:=[11], eigenvalue:=-1,
        cuspidalName:="{}^2F_4^3[-1]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
        levi:=[1..4], parameterExponents:=[],
  	charNumbers:=[14], eigenvalue:=-E(4),
        cuspidalName:="{}^2F_4[-i]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
        levi:=[1..4], parameterExponents:=[],
  	charNumbers:=[15], eigenvalue:=E(4),
        cuspidalName:="{}^2F_4[i]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
        levi:=[1..4], parameterExponents:=[],
  	charNumbers:=[16], eigenvalue:=E(4),
        cuspidalName:="{}^2F_4^2[i]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
        levi:=[1..4], parameterExponents:=[],
  	charNumbers:=[17], eigenvalue:=-E(4),
        cuspidalName:="{}^2F_4^2[-i]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
        levi:=[1..4], parameterExponents:=[],
  	charNumbers:=[18], eigenvalue:=-E(3)^2,
        cuspidalName:="{}^2F_4[-\\zeta_3^2]"),
    rec(relativeType:=rec(series:="A",indices:=[],rank:=0),
        levi:=[1..4], parameterExponents:=[],
  	charNumbers:=[21], eigenvalue:=-1,
        cuspidalName:="{}^2F_4^4[-1]")],
  families:=[
    Family("C1",[ 1 ]),
    Family("C1",[ 2 ]),
    Family("C1",[ 4 ]),
    Family(rec(name:="C''_1",group:="C1", charLabels:=[""], 
      fourierMat:=[[-1]], eigenvalues:=[1], sh:=[1]),[5]),
    Family("X2", [9,12]),
    Family("X2", [10,13]),
    Family(rec(name:="X13",
	fourierMat:=Permuted(TransposedMat(
         [[ 3, 0,-6,3, ER(2)*3, ER(2)*3, ER(2)*3, ER(2)*3, 3,-3, 0, 0, 0],
          [ 3, 0,-6,3,-ER(2)*3,-ER(2)*3,-ER(2)*3,-ER(2)*3, 3,-3, 0, 0, 0],
          [ 6, 0, 0,6,       0,       0,       0,       0,-6, 6, 0, 0, 0],
          [-6,-4,-4,2,       0,       0,       0,       0,-6,-2,-4,-4, 0],
          [-3, 4,-2,1,-ER(2)*3, ER(2)*3, ER(2)*3,-ER(2)*3,-3,-1, 4, 4, 0],
          [-3, 4,-2,1, ER(2)*3,-ER(2)*3,-ER(2)*3, ER(2)*3,-3,-1, 4, 4, 0],
          [-3, 0, 0,3, ER(2)*3,-ER(2)*3, ER(2)*3,-ER(2)*3, 3, 3, 0, 0, 6],
          [-3, 0, 0,3, ER(2)*3, ER(2)*3,-ER(2)*3,-ER(2)*3, 3, 3, 0, 0,-6],
          [-3, 0, 0,3,-ER(2)*3,-ER(2)*3, ER(2)*3, ER(2)*3, 3, 3, 0, 0,-6],
          [-3, 0, 0,3,-ER(2)*3, ER(2)*3,-ER(2)*3, ER(2)*3, 3, 3, 0, 0, 6],
          [ 0, 4, 4,4,       0,       0,       0,       0, 0,-4, 4,-8, 0],
          [ 0, 4, 4,4,       0,       0,       0,       0, 0,-4,-8, 4, 0],
          [ 0,-8, 4,4,       0,       0,       0,       0, 0,-4, 4, 4, 0]]/12),
          (8,13)),
        special:=4,
        eigenvalues:=[1,1,1,-1,-1,-1,-E(4),E(4),E(4),-E(4),-E(3),-E(3)^2,-1],
        sh:=[1,1,1,1,1,E(4),-E(4),-1,1,1,E(3),E(3)^2,-1]),
        [3,19,20,7,8,11,14,15,16,17,6,18,21])],
  almostHarishChandra:=[
    rec(relativeType:=rec(
         orbit:=[rec(series:="F",cartanType:=ER(2),indices:=[1..4],rank:=4)],
         twist:=(1,4)(2,3)),
         levi:=[], eigenvalue:=1, cuspidalName:="",
      charNumbers:= [1,2,3,4,5,19,20,7,9,10,8]),
    rec(relativeType:=rec(orbit:=[rec(series:="B",indices:=[1,4],
         cartanType:=ER(2),rank:=2)],twist:=(1,2)),levi:=[2,3],eigenvalue:=-1,
	 cuspidalName:="B_2",charNumbers:=[12,13,15]),
    rec(relativeType:=rec(series:="A",indices:=[  ],rank:=0),
 levi:=[ 1 .. 4 ], eigenvalue:=-1, 
   cuspidalName:="F_4[-1]", charNumbers:=[ 21 ]),
    rec(relativeType:=rec(series:="A",indices:=[  ],rank:=0),
 levi:=[ 1 .. 4 ], eigenvalue:=-E(4), 
   cuspidalName:="F_4[-i]", charNumbers:=[ 14 ]),
    rec(relativeType:=rec(series:="A",indices:=[  ],rank:=0),
 levi:=[ 1 .. 4 ], eigenvalue:=E(4), 
   cuspidalName:="F_4[i]", charNumbers:=[ 11 ]),
    rec(relativeType:=rec(series:="A",indices:=[  ],rank:=0),
 levi:=[ 1 .. 4 ], eigenvalue:=E(3), 
   cuspidalName:="F_4[\\zeta_3]", charNumbers:=[ 6 ]),
    rec(relativeType:=rec(series:="A",indices:=[  ],rank:=0),
 levi:=[ 1 .. 4 ], eigenvalue:=E(3)^2, 
   cuspidalName:="F_4[\\zeta_3^2]", charNumbers:=[ 18 ]),
    rec(relativeType:=rec(series:="A",indices:=[  ],rank:=0),
 levi:=[ 1 .. 4 ], eigenvalue:=1, 
   cuspidalName:="F_4[1]", charNumbers:=[ 16 ]),
    rec(relativeType:=rec(series:="A",indices:=[  ],rank:=0),
 levi:=[ 1 .. 4 ], eigenvalue:=1, 
   cuspidalName:="F_4^2[1]", charNumbers:=[ 17 ])],
 a:=[ 0, 24, 4, 2, 10, 4, 4, 4, 1, 13, 4, 1, 13, 4, 4, 4, 4, 4, 4, 4, 4 ],
 A:=[ 0,24,20,14,22,20,20,20,11,23,20,11,23,20,20,20,20,20,20, 20, 20 ]);
end);
