#############################################################################
##
#A  tbl/weyl3d4.g            CHEVIE library                       Jean Michel
##
#Y  Copyright (C) 1994 - 2001  The CHEVIE Team
##
##  This file contains information about the coset W.F where W is
##  the group D4:   2      and F is (1,2,4).
##                1-3-4
##
##  The class names come from F.Luebeck's tables.
##
CHEVIE.AddData("cyclestructure","3D4",[[2,,,,3],[,6],[3,,,,3],[ 3, 2,,, 2 ],
  [ ,,,,,,,,,, 2 ],[ , 8 ], [ ,,,, 4 ] ]);

CHEVIE.AddData("generators","3D4",
  [(1,13)(3,5)(6,8)(7,9)(10,11)(15,17)(18,20)(19,21)(22,23),
    ( 2,14)( 3, 6)( 5, 8)( 7,10)( 9,11)(15,18)(17,20)(19,22)(21,23),
    ( 1, 5)( 2, 6)( 3,15)( 4, 7)(11,12)(13,17)(14,18)(16,19)(23,24),
    ( 3, 7)( 4,16)( 5, 9)( 6,10)( 8,11)(15,19)(17,21)(18,22)(20,23) ]);

CHEVIE.AddData("phi","3D4",(1,2,4)(5,6,7)(8,10,9)
  (13,14,16)(17,18,19)(20,22,21));

CHEVIE.AddData("CharName","3D4",function(arg)return
  PartitionTupleToString(arg[1]);end);

CHEVIE.AddData("CartanMat","3D4",[[2,0,-1,0],[0,2,-1,0],[-1,-1,2,-1],
  [0,0,-1,2]]);

#############################################################################
##
#F  CHEVIE.R("ClassParameter","3D4")( <w> ) . . . . class parameter of w
##
##  given an element w  of a Coxeter group W of type D_4 as word in standard
##  generators, returns the classparam of its F-conjugacy class under the
##  nontrivial F action permuting the generators 1->2->4->1.
##
CHEVIE.AddData("ClassParameter","3D4",function(w)local x;
  if w=[] then x:=();
  else x:=Product(CHEVIE.R("generators","3D4"){w});
  fi;
  return CHEVIE.R("classparams","3D4")
   [Position(CHEVIE.R("cyclestructure","3D4"),
                          CycleStructurePerm(x*CHEVIE.R("phi","3D4")))];
end);

CHEVIE.AddData("ClassInfo","3D4",function()local res;
  res:=rec(
  classtext :=[[1],[],[1,2,3,1,2,3],[3],[1,3],[1,2,3,1,2,4,3,2],[1,2,3,2]],
  classnames:=["C_3","\\tilde A_2","C_3+A_1","\\tilde A_2+A_1","F_4",
     "\\tilde A_2+A_2","F_4(a_1)"],
  orders:=[ 6, 3, 6, 6, 12, 3, 6 ],
  classes:=[ 48, 16, 16, 48, 48, 8, 8 ]);
  res.classparams:=res.classnames;
  return res;
end);

CHEVIE.AddData("NrConjugacyClasses","3D4",7);

CHEVIE.AddData("CharInfo","3D4",function()return
  rec(extRefl:=[1,5,4,6,2],
      charparams:=[[[],[4]],[[],[1,1,1,1]],[[],[2,2]],
          [[1,1],[2]],[[1],[3]],[[1],[1,1,1]],[[1],[2,1]]],
      charRestrictions:=[13,4,10,5,11,3,6],
      nrGroupClasses:=13,
      b:=[ 0, 12, 4, 4, 1, 7, 3 ],
      B:=[ 0, 12, 8, 8, 5, 11, 9 ]);end);

###########################################################################
##
#F  HeckeCharTable( <q> )  . . . . . . . . . . . . . . . . . . . .
##  . . . . . . . . . . .  character table of the Hecke algebra of type 3D4
##
##  HeckeCharTable returns the character table of the Hecke algebra
##  of type  3D4  with  parameter  q.
##
CHEVIE.AddData("HeckeCharTable","3D4",function(param,sqrtparam)local q,tbl;
  q:=-param[1][1]/param[1][2];
  tbl:=rec(
  identifier:="H(3D4)",
  parameter:=[q,q,q,q],
  sqrtparameter:=[],
  cartan:=CHEVIE.R("CartanMat","3D4"),
  size:=192,
  irreducibles:=[ [   q, 1,           q^6,   q, q^2,    q^8,    q^4 ],
		  [  -1, 1,             1,  -1,   1,      1,      1 ],
		  [ q-1, 2,         2*q^3, q-1,  -q,   -q^4,   -q^2 ],
		  [   0, 0, q^4-2*q^3+q^2,   0,  -q,  3*q^4,  3*q^2 ],
		  [   q, 1,     q^5-2*q^4,  -1,   0, -2*q^6,  2*q^3 ],
		  [  -1, 1,      -2*q^2+q,   q,   0, -2*q^2,    2*q ],
		  [ q-1, 2,      -q^4-q^2, q-1,   0,  2*q^4, -2*q^2 ] ]*q^0,
  irredinfo:=CHEVIE.R("IrredInfo","3D4"));
  Inherit(tbl,CHEVIE.R("ClassInfo","3D4")());
  tbl.centralizers:=List(tbl.classes,x->tbl.size/x);
  tbl := CHEVIE.compat.MakeCharacterTable(tbl);
  CHEVIE.compat.AdjustHeckeCharTable(tbl,param);
  return tbl;
end);

CHEVIE.AddData("PhiFactors","3D4",[1,E(3),1,E(3)^2]);

CHEVIE.AddData("Representation","3D4",function(i)
  return CHEVIE.R("HeckeRepresentation","3D4")
    ([[1,-1],[1,-1],[1,-1],[1,-1]],[1,1,1,1],i);end);

CHEVIE.AddData("HeckeRepresentation","3D4",
  function(param,sqrtparam,i)local q,v,res,x;
  q:=-param[1][1]/param[1][2];
  if not IsBound(sqrtparam[1]) then
    v:=GetRoot(q,2,"Representation(Hecke(3D4)[",i,"])");
  else v:=-sqrtparam[1]/param[1][2];
  fi;
  res:=[ rec(gens:=[[[ v^2 ]],[[ v^2 ]],[[ v^2 ]],[[ v^2 ]]],F :=[[ 1 ]]),
  rec(gens:=[[[ -1 ] ],[[ -1 ] ],[[ -1 ] ],[[ -1 ] ] ],F :=[[ 1 ] ] ),
  rec(gens:=[[[-1,v],[0,v^2]],[[-1,v],[0,v^2]],[[v^2,0],[v,-1]],
             [[-1,v],[0,v^2]]],F :=[[ 1,0 ],[ 0,1 ] ] ),
  rec(gens:= [[[ v^2,0,0,0,0,0 ],[ 0,-1,0,0,0,0 ],[ 0,0,v^2,0,0,0 ],
	       [ 0,0,0,v^2,0,0 ],[ 0,0,v,0,-1,0 ],[ 0,0,0,v,0,-1 ] ],
	      [[ -1,0,0,v,0,0 ],[ 0,v^2,0,0,0,0 ],[ 0,0,-1,0,0,0 ],
	       [ 0,0,0,v^2,0,0 ],[ 0,v,0,0,-1,0 ],[ 0,0,0,0,0,v^2 ] ],
	      [[ v^2,0,0,0,0,0 ],[ 0,-1,0,0,v,v ],[ -v,0,-1,0,v,0 ],
	       [ v,0,0,-1,0,v ], [ 0,0,0,0,v^2,0 ],[ 0,0,0,0,0,v^2 ] ],
	      [[ -1,0,-v,0,0,0 ],[ 0,v^2,0,0,0,0 ],[ 0,0,v^2,0,0,0 ],
	       [ 0,0,0,-1,0,0 ],[ 0,0,0,0,v^2,0 ],[ 0,v,0,0,0,-1 ] ] ],
	  F :=[[ 0,0,0,0,0,-1 ],[ 0,0,-1,0,0,0 ],[ 0,0,0,1,0,0 ],
	       [ 0,-1,0,0,0,0 ],[ 1,0,0,0,0,0 ],[ 0,0,0,0,-1,0 ] ] ),
  rec( gens := [[[ -1,v,0,0 ],[ 0,v^2,0,0 ],[ 0,0,v^2,0 ],[ 0,0,0,v^2 ] ],
               [[ v^2,0,0,0 ],[ 0,v^2,0,0 ],[ 0,v,-1,0 ],[ 0,0,0,v^2 ] ],
	       [[ v^2,0,0,0 ],[ v,-1,v,v ], [ 0,0,v^2,0 ],[ 0,0,0,v^2 ] ],
	       [[ v^2,0,0,0 ],[ 0,v^2,0,0 ], [ 0,0,v^2,0 ],[ 0,v,0,-1 ] ] ],
           F :=[[ 0,0,1,0 ],[ 0,1,0,0 ],[ 0,0,0,1 ],[ 1,0,0,0 ] ]),
  rec( gens := [[[ -1,0,0,0 ],[ 0,v^2,0,0 ], [ 0,v,-1,0 ],[ 0,0,0,-1 ] ],
		[[ v^2,0,0,0 ],[ 0,-1,0,0 ], [ v,0,-1,0 ],[ 0,0,0,-1 ] ],
		[[ -1,0,v,0 ],[ 0,-1,v,0 ], [ 0,0,v^2,0 ],[ 0,0,v,-1 ] ],
		[[ -1,0,0,0 ],[ 0,-1,0,0 ], [ 0,0,-1,v ],[ 0,0,0,v^2 ] ] ],
            F :=[[ 0,0,0,1 ],[ 1,0,0,0 ],[ 0,0,1,0 ],[ 0,1,0,0 ] ] ),
  rec( gens :=[[[-1,v,0,0,0,0,0,2*v],[0,v^2,0,0,0,0,0,0],[0,v,-1,0,0,0,0,0],
                [0,0,0,v^2,0,0,0,0],[0,0,0,0,v^2,0,0,0],[0,0,0,v,0,-1,0,0],
                [0,0,0,0,v,0,-1,0],[0,0,0,0,0,0,0,v^2]],
               [[-1,v,0,0,0,0,2*v,0],[0,v^2,0,0,0,0,0,0],[0,0,v^2,0,0,0,0,0],
                [0,v,0,-1,0,0,0,0],[0,0,0,0,v^2,0,0,0],[0,0,v,0,0,-1,0,0],
                [0,0,0,0,0,0,v^2,0],[0,0,0,0,v,0,0,-1]],
               [[v^2,0,0,0,0,0,0,0],[v,-1,0,0,0,0,0,0],[0,0,-1,0,0,v,v,0],
                [0,0,0,-1,0,v,0,v],[0,0,0,0,-1,0,v,v],[0,0,0,0,0,v^2,0,0],
                [0,0,0,0,0,0,v^2,0],[0,0,0,0,0,0,0,v^2]],
               [[-1,v,0,0,0,2*v,0,0],[0,v^2,0,0,0,0,0,0],[0,0,v^2,0,0,0,0,0],
	        [0,0,0,v^2,0,0,0,0],[0,v,0,0,-1,0,0,0],[0,0,0,0,0,v^2,0,0],
                [0,0,v,0,0,0,-1,0],[0,0,0,v,0,0,0,-1]]],
            F:=[[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,0,1,0,0,0,0],
                [0,0,0,0,1,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,1],
                [0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0]])];
  res:=res[i];res.gens:=res.gens*v^0;
  return res;
end);

CHEVIE.AddData("UnipotentCharacters","3D4",function()local res;
  return rec(
  harishChandra:=[
    rec(relativeType:=rec(series:="G",indices:=[1,3],rank:=2),
      levi:=[  ], eigenvalue:=1,
      parameterExponents:=[3,1], cuspidalName:="",
      charNumbers:=[1,2,6,5,7,3]),
    rec(relativeType:=rec(series:="A",indices:=[  ],rank:=0),
      levi:=[ 1 .. 4 ], eigenvalue:=1,
      parameterExponents:=[], cuspidalName:="{}^3D_4[1]",
      charNumbers:=[4]),
    rec(relativeType:=rec(series:="A",indices:=[  ],rank:=0),
      levi:=[ 1 .. 4 ], eigenvalue:=-1,
      parameterExponents:=[], cuspidalName:="{}^3D_4[-1]",
      charNumbers:=[8])],
  almostHarishChandra:=[
    rec(relativeType:=rec(orbit:=[rec(series:="D",indices:=[1..4],rank:=4)],
      twist:=(1,2,4)),
      levi:=[], eigenvalue:=1,
      cuspidalName:="",
      charNumbers:=[1..7]),
    rec(relativeType:=rec(orbit:=[rec(series:="A",indices:=[  ],rank:=0)],
          twist:=()),
      levi:=[ 1 .. 4 ], eigenvalue:=-1,
      cuspidalName:="D_4",
      charNumbers:=[ 8 ])],
  families:=[
    Family("C1",[1]),
    Family("C1",[2]),
    Family("C1",[5]),
    Family("C1",[6]),
    Family("C2",[7, 4, 3, 8 ])],
  a:=[ 0, 12, 3, 3, 1, 7, 3, 3 ],
  A:=[ 0, 12, 9, 9, 5, 11, 9, 9 ]);
end);

CHEVIE.AddData("UnipotentClasses","3D4",function(p)local uc,p,c,g,class;
  class:=n->First(uc.classes,x->x.name=n);
  uc:=Copy(CHEVIE.R("UnipotentClasses","D")(4,p));
  for p in [["11111111",(1,2,4)],["221111",(1,2,3)]] do class(p[1]).F:=p[2];od;
  c:=class("3311");
  g:=CoxeterGroup("G",2);
  c.red:=ReflectionSubgroup(g,[]);
  c.F:=(g.1*g.2)^2;
  c.AuAction:=ExtendedReflectionGroup(c.red,(g.1*g.2)^3);
  return uc;
end);
