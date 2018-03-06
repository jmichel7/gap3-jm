#############################################################################
##
#A  affine.g              CHEVIE library        Jean Michel, Raphael Rouquier
##
##
#Y  Copyright (C) 1992 - 2000  Lehrstuhl D f\"ur Mathematik, RWTH Aachen, IWR
#Y  der Universit\"at Heidelberg, University of St. Andrews, and   University 
#Y  Paris VII.
##
##  This file contains functions for affine Weyl groups.
##  Since they are Coxeter groups, all functions for general Coxeter groups
##  like Hecke algebras, hecke modules etc.. are supported by AbsCoxOps
##
##  Following   Kac  "Infinite  dimensional  Lie  algebras",  we  represent
##  matrices for an Affine group corresponding to a Weyl group of rank n by
##  n+2  by n+2 matrices. The upper (n+1)x(n+1) left corner of the matrices
##  form  another representation, but in which  it is harder to compute the
##  AffineRootAction.

AffineCoxeterGroupOps:=OperationsRecord("AffineCoxeterGroupOps",GenCoxOps);

AffineCoxeterGroupOps.String:=W->SPrint("Affine(",W.linear,")");

AffineCoxeterGroupOps.Print:=function(W)Print(String(W));end;

AffineCoxeterGroupOps.ReflectionName:=function(W,opt)
  return SPrint("~",ReflectionName(W.linear,opt));
end;

#############################################################################
##
#F  Affine(<W>) . . . . . . . . . . . Makes the affine group corresponding
#F    to the irreducible Weyl group W
##
Affine:=function(W)local tmp,Wa,i,j,n;
  if Length(ReflectionType(W))>1 then
    Error("Affine(W) implemented only for irreducible W");
  fi;
# if not Forall(W.cartan,IsInt) then
#   Error("Affine(W) implemented only for crystallographic W");
# fi;
  tmp:=rec();
  tmp.operations:=AffineCoxeterGroupOps;
  tmp.nbGeneratingReflections:=W.nbGeneratingReflections+1;
  n:=tmp.nbGeneratingReflections;
  if
  Set(W.rootInclusion{[1..n-1]})<>[1..n-1]
  then Error("Affine Coxeter groups implemented only for W such that\n",
  "Set(W.rootInclusion{W.generatingReflections})=[1..W.nbGeneratingReflections]");
  fi;
  tmp.reflectionsLabels:=[1..n];
##
## Extended Cartan matrix
##
  tmp.cartan:=List(W.cartan,ShallowCopy);
  for i in [1..W.nbGeneratingReflections] do 
               Add(tmp.cartan[i],-W.cartan[i]*W.roots[W.N]);od;
  Add(tmp.cartan,Concatenation(-W.coroots[W.N]*W.cartan,[2]));
##
## Action of the affine Weyl tmp group on the dual Cartan algebra
## (a n+1 dimensional vector space) following Kac sections 1.1 and 3.7
##
  tmp.identity:=IdentityMat(n+1);
  tmp.rank:=n+1;
  tmp.reflections:=List([1..n],i->List(tmp.identity,ShallowCopy));
  for i in [1..n] do for j in [1..n] do 
    tmp.reflections[i][j][i]:=tmp.reflections[i][j][i]-tmp.cartan[i][j];
  od;od;
  tmp.reflections[n][n+1][n]:=-1;
  Wa:=ApplyFunc(Group,tmp.reflections);
  Inherit(Wa,tmp);
  Wa.linear:=W;
  Wa.isFinite:=false;
  if W.semisimpleRank=1 then Wa.coxeterMat:=[[1,0],[0,1]];fi;#~A1 infinite bond!
  AbsCoxOps.CompleteCoxeterGroupRecord(Wa);
  return Wa;
end;

## Given an affine Weyl group W, x a vector in the basis of 
## simple roots of W.linear and w in W, returns the image of x under w.
##
AffineRootAction:=function(W,w,x) local y;
  y:=Concatenation(x,[0,1])*w;
  return y{[1..Length(x)]}-y[Length(x)+1]*W.linear.roots[W.linear.N];
end;

AffineCoxeterGroupOps.IsLeftDescending:=function(W,w,i)
  return Sum(w[i])<0;
end;

CHEVIE.AddData("PrintDiagram","AffineA",function(v)local s,n,r,o;
  r:=Length(v)-1;
  if r=1 then Print("A1~  ",v[1]," oo ",v[2],"\n");
  else 
    n:=SPrint("A",r,"~   ");
    s:=SPrint(Join(List(v{[1..r]},Format)," - "));o:=Length(s)-4;
    Print(String("",Length(n)),"  ",
       Join(List([1..QuoInt(o,4)],i->"- "),""),Format(v[r+1]),
       Join(List([1..QuoInt(o,4)],i->" -"),""),"\n");
    Print(String("",Length(n))," /",String("",o),"\\\n");
    Print(n,s,"\n");
  fi;
end);

CHEVIE.AddData("PrintDiagram","AffineD",function(v)local i,r;
  r:=Length(v)-1;
  Print("D",r,"~  ",v[1],String("",4*(r-3)-1),v[r+1],"\n");
  Print(String("\\",Length(String(r))+6),String("",4*(r-3)-3),"/");
  Print(String("",Length(String(r))+6),"\n",
		String("",Length(String(r))+6),v[3]);
  for i in [4..r-1] do Print(" - ",v[i]);od;Print("\n");
  Print(String("/",Length(String(r))+6),String("",4*(r-3)-3),"\\\n");
  Print(String("",Length(String(r))+4),v[2],String("",4*(r-3)-1),v[r],"\n");
end);

CHEVIE.AddData("PrintDiagram","AffineI",function(v,bond)
  if bond mod 2 = 1 then
    Print(String("",7+Length(String(bond)))," ",v[3],"\n");
    Print("      ",bond," /",bond,"\\ ",bond,"\n");
    Print("I2(",bond,")~ ",v[1]," - ",v[2],"\n");
  elif bond mod 4 =2 then 
    Print(String("",7+Length(String(bond))),bond,"\n");
    Print("I2(",bond,")  ",v[1]," < ",v[2]," - ",v[3],"\n");
  else 
    Print(String("",8+Length(String(bond))),bond,
	  String("",4-Length(String(bond))),bond,"\n");
    Print("I2(",bond,")  ",v[3]," > ",v[1]," < ",v[2],"\n");
  fi;
end);

#############################################################################
#F  PrintDiagram( <rec> ) . . prints Dynkin diagram of affine Coxeter
##
AffineCoxeterGroupOps.PrintDiagram:=function(W)local i, v, r, n, t, a, o, s;
  a:=W.linear.type[1];
  v:=W.reflectionsLabels;r:=a.rank;t:=a.series;
  if t="I" then CHEVIE.R("PrintDiagram","AffineI")(W.reflectionsLabels,a.bond);
  elif t="B" then
    if (IsBound(a.cartanType) and a.cartanType=1) or r=2 then
      Print("C",r,"~   ",v[1]," > ",Join(v{[2..r]}," - ")," < ",v[r+1],"\n");
    else
      s:=SPrint("B",r,"~   ",v[1]," < ",Join(v{[2..r-1]}," - "));
      Print(String("",Length(s)-1),Format(v[r+1]),"\n");
      Print(String("",Length(s)-1),"|\n");
      Print(s," - ",Format(v[r]),"\n");
    fi;
  elif t="A" then CHEVIE.R("PrintDiagram","AffineA")(W.reflectionsLabels);
  elif t="D" then CHEVIE.R("PrintDiagram","AffineD")(W.reflectionsLabels);
  elif t="E" then
    if r=6 then Print("             ",v[7],"\n             |\n",
	  "             ",v[2],"\n             |\n",
          "E6~  ",Join(v{[1,3,4,5,6]}," - "),"\n");
    elif r=7 then Print("                 ",v[2],"\n                 |\nE7~  ",
	    Join(v{[8,1,3,4,5,6,7]}," - "),"\n");
    else Print("             ",v[2],"\n             |\nE8~  ",
	   Join(v{[1,3,4,5,6,7,8,9]}," - "),"\n");
    fi;
  elif t="G" then Print("G2~  ",v[3]," - ",v[1]," > ",v[2]," \n");
  elif t="F" then
    Print("F4~  ",v[5]," - ",v[1]," - ",v[2]," > ",v[3]," - ",v[4],"\n");
  elif t="H" then
    if r=3 then
      Print("     ",v[4],"\n    5 \\\nH3~    ",v[2]," - ",v[3],"\n",
            "    5 /\n     ",v[1],"\n");
    else Print("        5           5\nH",r,"    ",Join(v," - "),"\n");
    fi;
  fi;
end;

##############################################################################
# The following function follows Lewis-McCammond-Petersen-Schwer,
# "computing reflection length in an affine Coxeter group"
#
AffineCoxeterGroupOps.ReflectionLength:=function(W,w)
  local Id,mov,W0,plat,dimw,l,p;
  W0:=W.linear;
  Id:=Concatenation(IdentityMat(W0.semisimpleRank),[0*[1..W0.semisimpleRank]]);
  mov:=List(Id,v->AffineRootAction(W,w,v)-v);
  l:=Concatenation(List(W0.generatingReflections,i->Reflection(W0,i)),
     [Reflection(W0,W0.N)]);
  p:=ReflectionLength(W0,Product(CoxeterWord(W,w),i->l[i]));
  dimw:=Minimum(List(Filtered(ParabolicSubgroups(W0),
    x->RankMat(Concatenation(mov,W0.roots{x}))=Length(x)),Length));
  return 2*dimw-p;
end;
