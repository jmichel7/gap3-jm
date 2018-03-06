###########################################################################
##
#A  xy.g     The CHEVIE package        Jean Michel and Raphael Rouquier
##
#Y  Copyright (C)  1997   Equipe des groupes finis, Universite Paris VII
##
##  Display graphically elements of Hecke modules for affine Weyl groups
##  of rank 2.
##
# An example of use:
# W:=Affine(CoxeterGroup("G",2));
# l:=Filtered(CoxeterWords(W,15),x->x[1]=3);
# q:=X(Rationals);q.name:="q";
# H:=Hecke(W,q^2);
# MT:=ModuleBasis(H,"MT");
# MC:=ModuleBasis(H,"MC'");
# XYDvi(W,MT(MC(l[1])));
###########################################################################
##
#F  XYPrint( <W> ,<a>) 
##   
##  <W> should be an affine Weyl group of rank 2.
##
## Given  a=[elm,coeff], returns an  XYPic program to  print in the alcoves
## parametrized by the elements of elm the corresponding coeff.
## elm should be a list of elements of W,
## all of them being of minimal length in the left class modulo
## the finite Weyl group (ie, correspond to dominant alcoves)
##
XYPrint:=function(W,a) local pt,ori,res,i,texpol,w1,w2,type;
  if W.linear.semisimpleRank<>2 then Error(W," should have rank 2"); fi;
  if IsRec(a) then a:=[a.elm,a.coeff];fi;
  a[2]:=List(a[2],i->FormatTeX(i));
  ori:=[1/3,1/3]*TransposedMat(W.linear.cartan)^-1;
## [1/3,1/3] may be more specifically adapted to the different types ?
## Set the dominant weight basis
  res:="$$\\xy\n<2cm,0cm>:";
  type:=W.linear.type[1].series;
  if type="A" then Append(res,"<1cm,1.732cm>::\n");
  elif type="C" then Append(res,"<1cm,1cm>::\n");
  elif type="G" then Append(res,"<1cm,0.577cm>::\n");
  fi;
  w1:=TransposedMat(W.linear.cartan^-1)[1];
  w2:=TransposedMat(W.linear.cartan^-1)[2];
## Determines the wall where to stop drawing
  pt:=Maximum(List(a[1],x->Maximum(
  [Sum(AffineRootAction(W,x,w1)*TransposedMat(W.linear.cartan)),
   Sum(AffineRootAction(W,x,[0,0])*TransposedMat(W.linear.cartan)),
   Sum(AffineRootAction(W,x,w2)*TransposedMat(W.linear.cartan))])));
## We draw the walls 
  if type="A" then
    for i in [1..pt] do
      PrintToString(res,"(0,",i,");(",i,",0)**\\dir{-},\n");
      PrintToString(res,"(",i-1,",0);(",i-1,",",pt-i+1,")**\\dir{-},\n");
      PrintToString(res,"(0,",i-1,");(",pt-i+1,",",i-1,")**\\dir{-},\n");
    od;
  elif type="C" then
    for i in [1..pt] do
      PrintToString(res,"(0,",i,");(",i,",0)**\\dir{-},\n");
      PrintToString(res,"(",i-1,",0);(",i-1,",",pt-i+1,")**\\dir{-},\n");
      PrintToString(res,"(",i,",0);(",Maximum(0,2*i-pt),",",Minimum(2*pt-2*i,2*i),")**\\dir{-},\n");
      PrintToString(res,"(0,",2*(i-1),");(",Maximum(0,pt-2*i+2),",",2*i-2,")**\\dir{-},\n");
    od;
  elif type="G" then
    for i in [1..pt] do
      PrintToString(res,"(0,",i,");(",i,",0)**\\dir{-},\n");
      PrintToString(res,"(",i-1,",0);(",i-1,",",pt-i+1,")**\\dir{-},\n");
    od;
  fi;
## We write the polynomials
  for i in [1..Length(a[1])] do
    pt:=AffineRootAction(W,a[1][i],ori)*TransposedMat(W.linear.cartan);
    Print(pt,"\n");
    if (type="A" or type="G") then
     PrintToString(res,"(",evalf(pt[1],4),",",evalf(pt[2],4),")*{",a[2][i],"},\n");
## Should be programmed more nicely
    elif type="C" then
     PrintToString(res,"(",evalf(pt[2],4),",",evalf(pt[1],4),")*{",a[2][i],"},\n");
    fi;
    od;
  PrintToString(res,"\\endxy$$\n"); 
  return res;
end;
## Runs automatically TeX + xdvi on the file produced by XYPrint
##
XYDvi:=function(W,a) local str;
  str:=Concatenation("\\nopagenumbers\n\\input xypic\n",XYPrint(W,a),"\\end\n");
  PrintTo("outXY.tex",str);
  Exec("tex outXY > /dev/null");
  Exec("xdvi -paper a4 -expert -keep outXY.dvi ");
  Exec("rm outXY.*");
end;
