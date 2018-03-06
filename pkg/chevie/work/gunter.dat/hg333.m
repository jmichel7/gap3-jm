h:=transpose(matrix(10,10,[
 [1,1,1,1,1,1,1,1,1,1],
 [1,x,x^2,x^2,x^2,x^3,x^3,x^4,x^6,x^6],
 [2,x+1,x^2+1,x,x,0,0,-x^2,-2*x^3,-2*x^3],
 [2,x+1,x,x^2+1,x,0,0,-x^2,-2*x^3,-2*x^3],
 [2,x+1,x,x,x^2+1,0,0,-x^2,-2*x^3,-2*x^3],
 [2,x+1,x,x,x,0,-(x+1)*(x^2-x+1),-x*(x^2+1),-2*x^3,(x^2+1)*(x^4-x^2+1)],
 [3,2+x,x+1,x+1,x+1,1/2*(1+sqrt(-3))*x,1/2*(-1+sqrt(-3))*x^2,1/2*(-1+sqrt(-3))*(x+1)*x,3/2*(-1+sqrt(-3))*x^2,-1/2*(1+sqrt(-3))*(x^3-2)*x],
 [3,2+x,x+1,x+1,x+1,1/2*(1-sqrt(-3))*x,1/2*(-1-sqrt(-3))*x^2,1/2*(-1-sqrt(-3))*(x+1)*x,3/2*(-1-sqrt(-3))*x^2,-1/2*(1-sqrt(-3))*(x^3-2)*x],
 [3,2*x+1,x*(x+1),x*(x+1),x*(x+1),-1/2*x^2*(-1+sqrt(-3)),-1/2*(1+sqrt(-3))*x,-1/2*x^2*(1+sqrt(-3))*(x+1),-3/2*(1+sqrt(-3))*x^4,-1/2*(-1+sqrt(-3))*(2*x^3-1)*x^2],
  [3,2*x+1,x*(x+1),x*(x+1),x*(x+1),-1/2*x^2*(-1-sqrt(-3)),-1/2*(1-sqrt(-3))*x,-1/2*x^2*(1-sqrt(-3))*(x+1),-3/2*(1-sqrt(-3))*x^4,-1/2*(-1-sqrt(-3))*(2*x^3-1)*x^2]
])):

#rank(h);

g:=vector(10,[]);
erg:=multiply(h,g):

lprint(`***`);
solve({erg[1]-1,erg[2],erg[3],erg[4],erg[5],erg[6],erg[7],erg[8],erg[9],erg[10]},
{g[1],g[2],g[3],g[4],g[5],g[6],g[7],g[8],g[9],g[10]});
lprint(");
quit;

Damit generische Grade:
g[1]:=x^9/(x-1)/(x^2+x+1)/(x^2-x+1)^3;
g[2]:=-1/(x-1)/(x^2+x+1)/(x^2-x+1)^3;
g[3]:=-x^3/(x^2-x+1)^3;
g[4]:=-x^3/(x^2-x+1)^3;
g[5]:=-x^3/(x^2-x+1)^3;
g[6]:=-x^3/(x^2-x+1)^3;
g[7]:=-1/3*sqrt(-3)*x^4*zb/(-1+x)/(x-z)/(x+z)^3;
g[8]:=1/3*sqrt(-3)*x^4*z/(-1+x)/(x-zb)/(x+zb)^3;
g[9]:=1/3*sqrt(-3)*x*zb/(-1+x)/(x-z)/(x+z)^3;
g[10]:=-1/3*sqrt(-3)*x*z/(-1+x)/(x-zb)/(x+zb)^3;


ind:=(q^2-1/(q^2+1)^2/(q^2-1)^3;
x:=-q^6; u:=q; v:=-q; w:=-q^4;


with(linalg):

T[1]:=diag(1):      #  IDENTITY
T[2]:=T1:             #  S1
T[3]:=multiply(T1,T2):        #  S1 S2
T[4]:=multiply(T1,T3):        #  S1 S3
T[5]:=multiply(T2,T3):        #  S2 S3
T[6]:=multiply(T[3],T3):      #  S1 S2 S3
T[7]:=multiply(T[4],T2):      #  S1 S3 S2
T[8]:=multiply(T[3],T[4]):    #  S1 S2 S1 S3
T[9]:=multiply(T[6],T[6]):    #  S1 S2 S3 S1 S2 S3
T[10]:=multiply(T[7],T[7]):   #  S1 S3 S2 S1 S3 S2

for i from 1 to 10 do lprint(factor(trace(T[i]))) od;

T1:=matrix(1,1,[[1]]); T2:=matrix(1,1,[[1]]); T3:=matrix(1,1,[[1]]);
T1:=matrix(1,1,[[x]]); T2:=matrix(1,1,[[x]]); T3:=matrix(1,1,[[x]]);

T1:=matrix(2,2,[[1,0],[a,x]]); T2:=matrix(2,2,[[1,0],[a,x]]);
T3:=matrix(2,2,[[x,-x/a],[0,1]]);

T1:=matrix(2,2,[[1,0],[a,x]]); T2:=matrix(2,2,[[x,-x/a],[0,1]]);
T3:=matrix(2,2,[[1,0],[a,x]]);

T1:=matrix(2,2,[[1,0],[a,x]]); T2:=matrix(2,2,[[x,-x/a],[0,1]]);
T3:=matrix(2,2,[[x,-x/a],[0,1]]);

T1:=matrix(2,2,[[1,0],[a,x]]); T2:=matrix(2,2,[[x,-x/a],[0,1]]);
T3:=matrix(2,2,[[1+x,-1/a],[a*x,0]]);

T1:=matrix(3,3,[[1,0,0],[0,1,0],[0,c,x]]); T2:=matrix(3,3,[[1,0,0],[0,x,-x/c],[0,0,1]]);
T3:=matrix(3,3,[[(x+1)*x/(x+zb),z*c*a3,a3],[z*(x-z)*x*(x+z)/c/a3/(x+zb)^2,-(x-z)*z/(x+zb),z*x/c/(x+zb)],[-(x+z)*(x-z)*x/a3/(x+zb)^2,-z*c*x/(x+zb),zb/(x+zb)]]);
z:=(-1+sqrt(-3))/2; zb:=z^2;

T1:=matrix(3,3,[[x,0,0],[0,1,0],[0,c,x]]); T2:=matrix(3,3,[[x,0,0],[0,x,-x/c],[0,0,1]]);
T3:=matrix(3,3,[[zb*(x+1)/(x+zb),zb*c*a3,a3],[-z*x*(x-z)*(x+z)/c/a3/(x+zb)^2,x*(x-z)/(x+zb),x*z/(x+zb)/c],[zb*(x+z)*(x-z)*x/a3/(x+zb)^2,-z*c*x/(x+zb),x^2/(x+zb)]]);


wg:FREE(s1,s2,s3);
wg.relations:s1^2=s2^2=s3^2=1, s1*s2*s1=s2*s1*s2, s1*s3*s1=s3*s1*s3,
  s2*s3*s2=s3*s2*s3, s1*s2*s3*s1*s2*s3=s3*s1*s2*s3*s1*s2;
print order(wg); print(classes(wg));
print character table(wg);    #  G(3,3,3)


with(linalg):

R1:=multiply(T1,T2):
S1:=multiply(R1,T1):
R2:=multiply(T2,T1):
S2:=multiply(R2,T2):
S3:=scalarmul(S2,-1):
S4:=add(S1,S3);

R1:=multiply(T1,T3):
S1:=multiply(R1,T1):
R2:=multiply(T3,T1):
S2:=multiply(R2,T3):
S3:=scalarmul(S2,-1):
S4:=add(S1,S3);

R1:=multiply(T1,T2):
R2:=multiply(T3,T1):
R3:=multiply(T2,T3):
S1:=multiply(multiply(R1,R2),R3):
S2:=multiply(multiply(R2,R3),R1):
S3:=scalarmul(S2,-1):
S4:=add(S1,S3);

for i to 3 do for j to 3 do lprint(factor(S4[i,j])) od; od;

R1:=multiply(T2,T3):
S1:=multiply(R1,T2):
R2:=multiply(T3,T2):
S2:=multiply(R2,T3):
S3:=scalarmul(S2,-1):
S4:=add(S1,S3);

#####################################################
#  Ist Normalteiler N erzeugt von 
  S3,  S2*S3*S2^-1,  S1*S2*S3*S2^-1*S1^-1
von WG_26. (Klasse der Spiegelungen)
Faktorgruppe  SL_2(3).
  mit Elementen der Ordnung 1, 2, 3, 4, 6.

  N + s1*s2^-1*s1*s2^-1 --> Gruppe der Ordnung 108 = 2*54,     20 Klassen (zentral)
  N + s1 --> Gruppe der Ordnung 162 = 3*54, n"amlich B_3^3,    22 Klassen 
  N + s1*s2*s1 --> Gruppe der Ordnung 216 = 4*54,              28 Klassen
  N + s1*s2^-1*s1*s2^-1, s1 --> Gruppe der Ordnung 324 = 6*54, 44 Klassen

Also ^2G = Ennola(G)
     ^3G durch Einbettung in B_3^3    H = B_2^3  (q^3-1)*(q^6-1)*(q^3-z^2)
     ^4G                              H = B_1^6  (q^6-1)*(q^6+1)
     ^6G = Ennola(^3G)



