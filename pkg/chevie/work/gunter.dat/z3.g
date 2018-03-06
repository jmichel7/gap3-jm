###################################
##  Unipotente Grade  Z3         ##
###################################
x:=X(Cyclotomics); x.name:="q";
#x:=X(NF([E(3)])); x.name:="x";
## pols:
p1:=x-1; p3:=x^2+x+1; 
## irrat:
z:=E(3); zb:=z^2; i3:=E(3)-E(3)^2;

Z3:=rec(
unpdeg:=[
 [1],
 [1/3*i3*zb*x*(x-zb),
  -1/3*i3*z*x*(x-z),
  -1/3*i3*x*p1]],
fakdeg:=[
 [1],
 [x,x^2]],
foumat:=[
 [[1]],
 1/3*i3*[[zb,-z,-1],[-z,zb,1],[-1,1,1]]
],
frob:=[
 [1],
 [1,1,zb]]
);
