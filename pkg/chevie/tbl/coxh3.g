#############################################################################
##
#A  tbl/coxh3.g                 CHEVIE library      Meinolf Geck, Jean Michel
##
#Y  Copyright (C) 1992 - 2004  The CHEVIE Team
##
##  This file contains data for Coxeter groups and Hecke algebras of type H3.
##
CHEVIE.AddData("ReflectionDegrees","H3",[2,6,10]);

CHEVIE.AddData("Size", "H3", 120);

# from Humphreys, "Reflection Groups and Coxeter Groups"
CHEVIE.AddData("GeneratingRoots", "H3", function()local a;a:=(1+ER(5))/2;
  return [[a,-1,a-1],[-a,1,a-1],[1,a-1,-a]]/2;
end);

# consider using the following Galois-stable model w.r.t. the automorphism
# Aut(H3,[[1,2,1,2,3,2,1,2,1],[3],[2]]), equal to the 
# usual model conjugate by [[(5-ER(5))/2,-ER(5),-1+ER(5)],
#   [-2*ER(5),1,0],[2*ER(5),1,0]]
#H3:=CoxeterGroup([[(5-ER(5))/2,-ER(5),-1+ER(5)],[-2*ER(5),1,0],[2*ER(5),1,0]],
#[[(5+ER(5))/40,(-1-ER(5))/4,(-3+ER(5))/16],[(-3*ER(5))/20,1/2,(2-ER(5))/8],
# [(3*ER(5))/20,1/2,(2+ER(5))/8]]);

CHEVIE.AddData("NrConjugacyClasses","H3", 10);

CHEVIE.AddData("cyclestructure","H3",[]);
CHEVIE.AddData("generators","H3",[]);

CHEVIE.AddData("CartanMat","H3",[[2,-(1+ER(5))/2,0],[-(1+ER(5))/2,2,-1],[0,-1,2]]);

CHEVIE.AddData("PowerMaps","H3",[,[1,1,7,1,5,3,3,5,7,1],[1,2,7,4,1,9,3,10,6,10],
,[ 1, 2, 1, 4, 5, 10, 1, 8, 10, 10 ],, [ 1, 2, 7, 4, 5, 9, 3, 8, 6, 10 ] ]);

CHEVIE.AddData("WordsClassRepresentatives","H3",
 [[],[1],[1,2],[1,3],[2,3],[1,2,3],[1,2,1,2],[1,2,1,2,3],
  [1,2,1,2,3,2,1,2,3],[1,2,1,2,1,3,2,1,2,1,3,2,1,2,3]]);

CHEVIE.AddData("ParabolicRepresentatives","H3",
function(s)local t;t:=[[[]],[[1]],[[1,2],[1,3],[2,3]],[[1..3]]];
  return t[s+1];end);

CHEVIE.AddData("ClassInfo","H3",function()local res;
  res:=rec(
  classtext:=CHEVIE.R("WordsClassRepresentatives","H3"),
  orders:=[1,2,5,2,3,10,5,6,10,2],
  classes:=[ 1, 15, 12, 15, 20, 12, 12, 20, 12, 1 ]);
  res.classnames:=List(res.classtext,IntListToString);res.classnames[1]:=".";
  res.classparams:=res.classnames;
  return res;
end);

CHEVIE.AddData("CharInfo","H3",function()local res;
  res:=rec(
  charparams:=[[1,15],[1,0],[5,5],[5,2],[3,6],[3,8],[3,1],[3,3],[4,3],[4,4]],
  gp:=["1_r'","1_r","5_r'","5_r","3_s","overline{3}_s","3_s'","overline{3}_s'",
     "4_r'","4_r"],
  opdam:=(9,10),
  extRefl:=[2,7,5,1]);
  res.b:=List(res.charparams,x->x[2]);
  return res;
end);

CHEVIE.AddData("vpolheckeirreducibles","H3",
[[[[1],0],[[-1],0],[[1],0],[[1],0],[[1],0],[[-1],0],[[1],0],[[-1],0],[[-1],0],
[[-1],0]],[[[1],0],[[1],1],[[1],2],[[1],2],[[1],2],[[1],3],[[1],4],[[1],5],
[[1],9],[[1],15]],[[[5],0],[[-3,2],0],[[1,-1],0],[[2,-2,1],0],[[1,-2],0],[[],
0],[[1,0,-1],0],[[1],2],[[],0],[[-5],6]],[[[5],0],[[-2,3],0],[[-1,1],1],[[1,
-2,2],0],[[-2,1],1],[[],0],[[-1,0,1],2],[[-1],3],[[],0],[[5],9]],[[[3],0],
[[-2,1],0],[[1,(-1+ER(5))/2],0],[[1,-2],0],[[1,-1],0],[[(1-ER(5))/2],1],[[1,0,
(-1-ER(5))/2],0],[[],0],[[(1+ER(5))/2],3],[[3],5]],[[[3],0],[[-2,1],0],[[1,
(-1-ER(5))/2],0],[[1,-2],0],[[1,-1],0],[[(1+ER(5))/2],1],[[1,0,(-1+ER(5))/2],
0],[[],0],[[(1-ER(5))/2],3],[[3],5]],[[[3],0],[[-1,2],0],[[(-1+ER(5))/2,1],1],
[[-2,1],1],[[-1,1],1],[[(-1+ER(5))/2],2],[[(-1-ER(5))/2,0,1],2],[[],0],
[[(-1-ER(5))/2],6],[[-3],10]],[[[3],0],[[-1,2],0],[[(-1-ER(5))/2,1],1],[[-2,
1],1],[[-1,1],1],[[(-1-ER(5))/2],2],[[(-1+ER(5))/2,0,1],2],[[],0],
[[(-1+ER(5))/2],6],[[-3],10]],[[[4],0],[[-2,2],0],[[-1],1],[[1,-2,1],0],[[1,
-1,1],0],[[1],3/2],[[-1],2],[[-1],5/2],[[1],9/2],[[-4],15/2]],[[[4],0],[[-2,
2],0],[[-1],1],[[1,-2,1],0],[[1,-1,1],0],[[-1],3/2],[[-1],2],[[1],5/2],[[-1],
9/2],[[4],15/2]]]);

CHEVIE.AddData("CycPolSchurElements","H3",
[ [1,-15,2,2,2,3,5,6,10], [1,0,2,2,2,3,5,6,10], [1,-5,2,2,2,3,6],
  [1,-2,2,2,2,3,6], [(5+ER(5))/2,-6,2,2,2,2/5,3/5,1/10,9/10],
  [(5-ER(5))/2,-6,2,2,2,1/5,4/5,3/10,7/10],
  [(5+ER(5))/2,-1,2,2,2,2/5,3/5,1/10,9/10],
  [(5-ER(5))/2,-1,2,2,2,1/5,4/5,3/10,7/10], [2,-3,3,5], [2,-3,3,5] ]
);

# description of fake degrees by powers of q
CHEVIE.AddData("sparseFakeDegrees","H3",[[1,15],[1,0],[1,5,1,7,1,9,1,11,1,13],
[1,2,1,4,1,6,1,8,1,10],[1,6,1,10,1,14],[1,8,1,10,1,12],[1,1,1,5,1,9],
[1,3,1,5,1,7],[1,3,1,7,1,9,1,11],[1,4,1,6,1,8,1,12]]);

###########################################################################
##
#F HeckeCharTable( <v> )  character table of the Hecke algebra of type H3
##

CHEVIE.AddData("HeckeCharTable","H3",function(param,sqrtparam)
 local a, q, v, ci, tbl;
 a:=(1+ER(5))/2;
 q:=-param[1][1]/param[1][2];
 if not IsBound(sqrtparam[1]) then v:=GetRoot(q,2,"CharTable(Hecke(H3))");
 else v:=-sqrtparam[1]/param[1][2];
 fi;
 ci := CHEVIE.R("ClassInfo", "H3")();
 tbl:=rec(
  identifier:="H(H3)",
  text:="the representing matrices are those of Lusztig (1981)",
  parameter := [ q, q, q ],
  cartan:=CHEVIE.R("CartanMat","H3"),
  size := 120, order := 120,
  powermap:=CHEVIE.R("PowerMaps","H3"),
  irreducibles:=List(CHEVIE.R("vpolheckeirreducibles","H3"),
      i->List(i,function(j)local res;res:=ValuePol(j[1],q);
       if IsInt(j[2]) then res:=res*q^j[2];else res:=res*v^(2*j[2]);fi;
       return res;end)),
  irredinfo:=CHEVIE.R("IrredInfo","H3"));
  Inherit(tbl,ci);
  tbl.centralizers:=List(tbl.classes,x->tbl.size/x);
  tbl := CHEVIE.compat.MakeCharacterTable(tbl);
  CHEVIE.compat.AdjustHeckeCharTable(tbl,param);
  return tbl;
end);

CHEVIE.AddData("Representation","H3",function(i)
  return CHEVIE.R("HeckeRepresentation","H3")
  ([[1,-1],[1,-1],[1,-1]],[1,1,1],i);end);

# W-graphs given by Ivan Marin.
# see the comments before WGraphToRepresentation for the format.
# An integer instead of a graph says that the graph is dual to that
# of index specified by the integer.
CHEVIE.AddData("WGraphs","H3",[[[[1,2,3]],[]],1,
 [[[2],[1,2],[1,3],[1,3],[2,3]],[[-1,[[1,3],[2,4],[3,5],[4,5]]]]],3,
 [[[1,2],[1,3],[2,3]],[[-1,[[1,2]]],[(-1-ER(5))/2,[[2,3]]]]],
 [[[1,2],[1,3],[2,3]],[[-1,[[1,2]]],[(-1+ER(5))/2,[[2,3]]]]],
 5,6,[[[1],[2],[1,3],[2,3]],[[1,[[1,2,3],[2,3,4],[3,4]]]]],9]);

CHEVIE.AddData("WGraph","H3",function(i)local gr;
 gr:=CHEVIE.R("WGraphs","H3");
 if IsInt(gr[i]) then return DualWGraph(3,gr[gr[i]]); else return gr[i]; fi;
end);

CHEVIE.AddData("HeckeRepresentation","H3",function(param,sqrtparam,i)local v;
 if not IsBound(sqrtparam[1]) then 
     v:=GetRoot(-param[1][1]/param[1][2],2,"Representation(Hecke(H3),[",i,"])");
 else v:=-sqrtparam[1]/param[1][2];
 fi;
 return -param[1][2]*WGraphToRepresentation(3,CHEVIE.R("WGraph","H3")(i),v);
end);

CHEVIE.AddData("UnipotentCharacters","H3",function()local res;
  res:=rec(
  harishChandra:=[
    rec(relativeType:=rec(series:="H",indices:=[ 1 .. 3 ],rank:=3),
 levi:=[  ], eigenvalue:=1,
      parameterExponents:=[ 1, 1, 1 ], cuspidalName:="",
      charNumbers:=[ 1 .. 10 ]),
    rec(relativeType:=rec(series:="A",indices:=[ 3 ],rank:=1),
 levi:=[ 1 .. 2 ], eigenvalue:=E(5)^2,
      parameterExponents:=[ 5 ], cuspidalName:="I_2(5)[1,3]",
      charNumbers:=[ 11, 13 ]),
    rec(relativeType:=rec(series:="A",indices:=[ 3 ],rank:=1),
 levi:=[ 1 .. 2 ], eigenvalue:=E(5)^3,
      parameterExponents:=[ 5 ], cuspidalName:="I_2(5)[1,2]",
      charNumbers:=[ 12, 14 ]),
    rec(relativeType:=rec(series:="A",indices:=[  ],rank:=0),
 levi:=[ 1 .. 3 ], eigenvalue:=E(4), qEigen:=1/2,
      parameterExponents:=[  ], cuspidalName:="H_3[i]",
      charNumbers:=[ 15 ]),
    rec(relativeType:=rec(series:="A",indices:=[  ],rank:=0),
 levi:=[ 1 .. 3 ], eigenvalue:=-E(4), qEigen:=1/2,
      parameterExponents:=[  ], cuspidalName:="H_3[-i]",
      charNumbers:=[ 16 ])],
  families:=[
    Family("C1",[ 2 ]),
    Family(CHEVIE.families.Dihedral(5),[ 7, 8, 14, 13 ],rec(ennola:=-1)),
    Family("C1",[ 4 ]),
    Family("C'\"2",[ 9, 10, 15, 16 ],rec(ennola:=3)),
    Family("C1",[ 3 ],rec(ennola:=-1)),
    Family(CHEVIE.families.Dihedral(5),[ 5, 6, 12, 11 ],rec(ennola:=1)),
    Family("C1",[ 1 ],rec(ennola:=-1))],
  a:=[ 15, 0, 5, 2, 6, 6, 1, 1, 3, 3, 6, 6, 1, 1, 3, 3 ],
  A:=[ 15, 0, 13, 10, 14, 14, 9, 9, 12, 12, 14, 14, 9, 9, 12, 12 ]);
  return res;
end);

CHEVIE.AddData("Invariants","H3",function()local r,C;
  C:=CHEVIE.R("CartanMat","H3");r:=RootsCartan(C)*C; 
 return List(CHEVIE.R("ReflectionDegrees","H3"),
   d-> function(arg) return Sum(r,a->(arg*a)^d);end);
end);

# the discriminant as a polynomial in the invariants
CHEVIE.AddData("Discriminant","H3",function()return function(a,b,c)
return 131835937500*a*b^3*c-100195312500*a^2*b*c^2+
395507812500*c^3-28369140625*a^3*b^4+1371093750*a^4*b^2*c+175781250000*b^5+
1191796875*a^5*c^2+1162187500*a^6*b^3-74250000*a^7*b*c-22233750*a^9*b^2+
438750*a^10*c+213700*a^12*b-829*a^15;
end;end);

CHEVIE.AddData("KLeftCellRepresentatives","H3",
[rec(character:=[2], duflo:=[1,2,3], reps:=""),
 rec(character:=[1], duflo:=[16,17,18], reps:=""),
 rec(character:=[3], duflo:=[1,24,3], reps:=""),
 rec(character:=[4], duflo:=[2,1,28], reps:=""),
 rec(character:=[6,5], duflo:=[1,20,18], reps:=[[7,19,24]]),
 rec(character:=[8,7], duflo:=[1,6,18], reps:=[[9,2,27]]),
 rec(character:=[10,9], duflo:=[8,18,17], 
   reps:=[[11,17,25],[11,27,10],[14,30,4]]),
 rec(character:=[10,9], duflo:=[13,30,8],
   reps:=[[10,29,5],[12,21,22],[13,22,23]])]);
