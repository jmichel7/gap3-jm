ReadVK("contrib/common");
W:=ComplexReflectionGroup(27);
res:=init(W);

if IsBound(VKCURVE.useChevieInvars) then
invars:=List(Invariants(W),f->ApplyFunc(f,vars));
res.H:=Hessian(invars[1],varnames);
#invars[2]:=myDet(res.H)/(25*27*(1+ER(-15)));
#invars[3]:=bord(invars[1],invars[2])/3240/(1+ER(-15));
elif false then
invars:=[ 10*x^3*y^3+9*z*(x^5+y^5)-45*x^2*y^2*z^2-135*x*y*z^4+27*z^6];
res.H:=Hessian(invars[1],varnames);
invars[2]:=myDet(res.H)/20250;
invars[3]:=bord(invars[1],invars[2])/24300;
else
invars:=[-90*x^2*y^2*z^2+180*x^2*y^3*z+30*x^2*y^4-135*x^2*z^4+135*y^2*z^4
+90*x^4*y*z-30*x^4*y^2+45*x^4*z^2+45*y^4*z^2+18*y^5*z+10*x^6-10*y^6-27*z^6];
res.H:=Hessian(invars[1],varnames);
invars[2]:=myDet(res.H)/81000;
invars[3]:=bord(invars[1],invars[2])/388800;
fi;

common24_33(res);

#gap> M;                 
#[ [ 6x, 12y^2, 30z+52/3x^3y ], 
#  [ 12y, -6z, -34/3xz-908/9x^2y^2+104/3y^3+200/9x^4y-20/9x^6 ], 
#  [ 30z, -34/3xyz-908/9x^2y^3+26/3x^3z+104/3y^4+200/9x^4y^2-20/9x^6y, 
#  6448/27xy^4+220/3y^2z-6698/27x^2yz-23636/81x^3y^3+1342/27x^4z+3160/81x^5\
#	  y^2+20/81x^7y-40/81x^9 ] ]
#  gap> infty(MinorMat(M));
#  [ [ -400/81x^12y, -160/27x^9y, -80/3x^6y^2 ], 
#   [ 880/27x^9y^2, -80/27x^10, -40/3x^7y ], [ -80/3x^6y^2, -40/3x^7,
#    -144y^3 ]  ]
# Donc (x,y,z)=(0,0,1) est le seul point possible a l'infini dans l'adherence
# du lieu singulier du discriminant.
# En particulier, la 2-section ci-dessous verifie les conditions
# de tranversalite a l'infini:

res.twoplane:=["x",1,"y",y,"z",x];
res.curve:=Value(res.disc,res.twoplane);

# VKcurve en exact donne au bout de 14h sur grobner2 (220 segments):
#
# 1: stst=tsts
# 2: tutut=ututu
# 3: susus=ususu
# 4: ustSus=stSust
# 5: tuTstu=stuTst
# 6: tustust=stustus
# 7: tusUtusU=usUtusUt
#
# t->stS gives length 56
#
# 1: utu=tut
# 2: stst=tsts
# 3: susus=ususu
# 4: tsTuts=utsTut
# 5: tsuStsuSts=suStsuStsu
#
# s->tsT gives length 40
#
# 1: sus=usu
# 2: utu=tut
# 3: stst=tsts
# 4: tusUtSustu=SustusUtus <=> tstustustust=stustustusts
#
# (stu)^5 is central
#
# s->usU gives length 38
#
# 1: sus=usu
# 2: utu=tut
# 3: tstst=ststs
# 4: suStUsut=tsuStsuS <=> tstsutsuts=stsutsutsu
#
# (sut)^5 is central
#
# u->suS gives length 36
#
# 1: sus=usu
# 2: utut=tutu
# 3: tstst=ststs
# 4: ustSus=stusUt <=> ustustu=stustut
#
# (stu)^5 is central
