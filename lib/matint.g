#############################################################################
##
#W  matint.gi                GAP library                        A. Storjohann
#W                                                              R. Wainwright
#W                                                                 F. Gaehler
#W                                                                    D. Holt
#W                                                                  A. Hulpke
##
#H  $Id: matint.gi,v 4.28.2.3 2004/03/17 11:49:50 gap Exp $
##
#Y  Copyright (C) 2003 The GAP Group
##
##  This file contains functions that compute Hermite and Smith normal forms 
##  of integer matrices, with or without the HNF/SNF  expressed as the linear 
##  combination of the input.
##
##  Backported to GAP3 by Jean Michel on 10th Feb. 2005
########################################################
##
##  	auxiliary + main code for all in one function
##
##  MATINTsplit 
##  MATINTrgcd
##  MATINTmgcdex
##  MATINTbezout
##  SNFofREF
##  NormalFormIntMat
##

Is_Zero:=x->x=0*x;

########################################################
#
# MATINTsplit(<N>,<a>) - returns product of prime factors of N which are not factors of a.
#
MATINTsplit:=function(N,a) local x,t;

  x:=a;
  t:=N;
  while x<>1 do
    x:=GcdInt(x,t);
    t:=QuoInt(t,x);
  od;
  return t;
end;

################################################
#
#   MATINTrgcd(<N>,<a>) - Returns smallest nonnegative c such that
#		   gcd(N,a+c) = 1
#
MATINTrgcd:=function(N,a)local k,r,d,i,c,g,q;

  if N=1 then return 0; fi;
  k := 1;
  r:=[(a-1) mod N];
  d:=[N];
  c:=0;
  while true do
    for i in [1..k] do r[i]:=(r[i]+1) mod d[i]; od;
    i:=PositionProperty(r,x->x<=0);
    if i=false then
      g:=1;i:=0; 
      while g=1 and i<k do
	i:=i+1;
	g:=GcdInt(r[i],d[i]);
      od;
      if g=1 then return c; fi;
      q:=MATINTsplit(QuoInt(d[i],g),g);
      if q>1 then
	k:=k+1;
	r[k]:=r[i] mod q;
	d[k]:=q;
      fi;
      r[i]:=0;
      d[i]:=g;    
    fi;
    c:=c+1;
  od;

end;

#######################################################
# 
#  MATINTmgcdex(<N>,<a>,<v>) - Returns c[1],c[2],...c[k] such that
#
#   gcd(N,a+c[1]*b[1]+...+c[n]*b[k]) = gcd(N,a,b[1],b[2],...,b[k])
#
MATINTmgcdex:=function(N,a,v) local h,g,M,c,i,d,b,l;
  l:=Length(v);
  c:=[]; M:=[];
  h := N;
  for i in [1..l] do
    g := h;
    h:=GcdInt(g,v[i]);  
    M[i]:=QuoInt(g,h);
  od;
  h:=GcdInt(a,h);
  g:=QuoInt(a,h);

  for i in [l,l-1..1] do
    b:=QuoInt(v[i],h);
    d:=MATINTsplit(M[i],b);
    if d=1 then
      c[i]:=0;
    else
      c[i]:=MATINTrgcd(d,g/b mod d);
      g:=g+c[i]*b;
    fi;
  od;

  return c;

end;




#####################################################
#
#  MATINTbezout(a,b,c,d) - returns row transform , P, to transform, A, to hnf :
#
#  PA=H;
#
#  [ s  t ] [ a  b ]   [ e  f ] 
#  [      ] [      ] = [       ]   
#  [ u  v ] [ c  d ]   [    g ]
#
MATINTbezout:=function(a,b,c,d) local e,f,g,q;

  e := Gcdex(a,c);
  f := e.coeff1*b+e.coeff2*d;
  g := e.coeff3*b+e.coeff4*d;
  if g<0 then
    e.coeff3 := -e.coeff3;
    e.coeff4 := -e.coeff4;
    g := -g; 
  fi;
  if g>0 then
    q := QuoInt(f-(f mod g),g);
    e.coeff1 := e.coeff1-q*e.coeff3;
    e.coeff2 := e.coeff2-q*e.coeff4;
  fi;
  return e;

end;

#####################################################
##
## SNFofREF - fast SNF of REF matrix
##
##
SNFofREF:=function(R,destroy)
local k,g,b,ii,m1,m2,t,tt,si,n,m,i,j,r,jj,piv,d,gt,tmp,A,T,TT,kk;

  n := Length(R);
  m := Length(R[1]);

  piv := List(R,x->PositionProperty(x,y->y<>0));
  r := PositionProperty(piv,x->x=false);
  if r=false then
    r := Length(piv);
  else
    r := r-1;  
    piv := piv{[1..r]}; 
  fi;
  Append(piv,Difference([1..m],piv));

  if destroy then
      T:=R;
      ##  Need to be careful: we are trying to permute the cols in place
      for i in [1..r] do
          T[i]{[1..m]} := T[i]{piv};
      od;
  else
      T := NullMat(n,m);
      for j in [1..m] do
          for i in [1..Minimum(r,j)] do T[i][j]:=R[i][piv[j]]; od;
      od;
  fi;

  si := 1;
  A := [];
  d := 2;
  for k in [1..m] do
    if k<=r then
	d := d*AbsInt(T[k][k]);
	Apply(T[k],x->x mod (2*d));
    fi;

   t := Minimum(k,r);
   for i in [t-1,t-2..si] do
      t := MATINTmgcdex(A[i],T[i][k],[T[i+1][k]])[1];
      if t<>0 then
         T[i]:=T[i]+T[i+1]*t; 
         Apply(T[i],x->x mod A[i]); 
      fi;
   od;

   for i in [si..Minimum(k-1,r)] do
      g := Gcdex(A[i],T[i][k]);
      T[i][k] := 0;
      if g.gcd<>A[i] then
         b := QuoInt(A[i],g.gcd);
         A[i] := g.gcd;
         for ii in [i+1..Minimum(k-1,r)] do
            T[ii]:=T[ii]+List(T[i]*-g.coeff2*QuoInt(T[ii][k],A[i]),x->x mod A[ii]);
            T[ii][k] := b*T[ii][k];

            Apply(T[ii],x->x mod A[ii]);
         od;
         if k<=r then 
            t := g.coeff2*QuoInt(T[k][k],g.gcd);
            T[k]:=T[k]+T[i]*-t;
            T[k][k]:=b*T[k][k];
         fi;
         Apply(T[i],x->x mod A[i]);
         if A[i]=1 then si := i+1; fi;
      fi;
   od;

   if k<=r then 
      A[k] := AbsInt(T[k][k]);
      Apply(T[k],x->x mod A[k]);
   fi;

  od;

  for i in [1..r] do T[i][i] := A[i]; od;

  return T;

end;

BITLISTS_NFIM:=
  [ [ false, false, false, false, false ], [ true, false, false, false, false ],
    [ false, true, false, false, false ], [ true, true, false, false, false ], 
    [ false, false, true, false, false ], [ true, false, true, false, false ], 
    [ false, true, true, false, false ], [ true, true, true, false, false ], 
    [ false, false, false, true, false ], [ true, false, false, true, false ], 
    [ false, true, false, true, false ], [ true, true, false, true, false ], 
    [ false, false, true, true, false ], [ true, false, true, true, false ], 
    [ false, true, true, true, false ], [ true, true, true, true, false ], 
    [ false, false, false, false, true ], [ true, false, false, false, true ], 
    [ false, true, false, false, true ], [ true, true, false, false, true ], 
    [ false, false, true, false, true ], [ true, false, true, false, true ], 
    [ false, true, true, false, true ], [ true, true, true, false, true ], 
    [ false, false, false, true, true ], [ true, false, false, true, true ], 
    [ false, true, false, true, true ], [ true, true, false, true, true ], 
    [ false, false, true, true, true ], [ true, false, true, true, true ], 
    [ false, true, true, true, true ], [ true, true, true, true, true ] ];


###########################################################
#
# DoNFIM(<mat>,<options>)
#
# Options bit values:
#
# 1  - Triangular / Smith
# 2  - No / Yes  Reduce off diag entries
# 4  - No / Yes  Row Transforms 
# 8  - No / Yes  Col Transforms
# 16 - change original matrix in place (The rows still change) -- save memory
#
# Compute a Triangular, Hermite or Smith form of the n x m 
# integer input matrix A.  Optionally, compute n x n / m x m
# unimodular transforming matrices which satisfy Q C A = H 
# or  Q C A B P = S.
#
# Triangular / Hermite :
#
# Let I be the min(r+1,n) x min(r+1,n) identity matrix with r = rank(A).
# Then Q and C can be written using a block decomposition as
#
#             [ Q1 |   ]  [ C1 | C2 ]
#             [----+---]  [----+----]  A  =  H.
#             [ Q2 | I ]  [    | I  ]
#
# Smith :
#
#  [ Q1 |   ]  [ C1 | C2 ]     [ B1 |   ]  [ P1 | P2 ]
#  [----+---]  [----+----]  A  [----+---]  [----+----] = S.
#  [ Q2 | I ]  [    | I  ]     [ B2 | I ]  [ *  | I  ]
#
# * - possible non-zero entry in upper right corner...
#				
#
DoNFIM:=function(mat,opt)
local sig, n, m, A, C, Q, B, P, r, c2, rp, c1, j, k, N, L, b, a, g, c,
      t, tmp, i, q, R, rank, signdet;

  if mat=[] then mat:=[[]];fi;
  if not (IsMatrix(mat) or (IsList(mat) and Length(mat)=1 
    and IsList(mat[1]) and Length(mat[1])=0)) or not IsInt(opt) then 
    Error("syntax is DoNFIM(<matrix>,<options>)"); 
  fi;

  opt := BITLISTS_NFIM[opt+1];
  sig:=1;

  #Embed mat in 2 larger "id" matrix
  n := Length(mat)+2;
  m := Length(mat[1])+2;
  k:=[1..m]*0;
  if opt[5] then
    # change the matrix
    A:=mat;
    for i in [n-1,n-2..2] do
      A[i]:=ShallowCopy(k);
      A[i]{[2..m-1]}:=A[i-1];
    od;
  else
    A := [];
    for i in [2..n-1] do
      #A[i] := [0];
      #Append(A[i],mat[i-1]);
      #A[i][m] := 0;
      A[i]:=ShallowCopy(k);
      A[i]{[2..m-1]}:=mat[i-1];
    od;
  fi;
  A[1]:=ShallowCopy(k);
  A[n]:=k;
  A[1][1] := 1;
  A[n][m] := 1;

  if opt[3] then 
    C := IdentityMat(n); 
    Q := NullMat(n,n);
    Q[1][1] := 1; 
  fi;

  if opt[1] and opt[4] then 
    B := IdentityMat(m);
    P := IdentityMat(m);
  fi;

  r := 0;
  c2 := 1;
  rp := [];
  while m>c2 do
    r := r+1;
    c1 := c2;
    rp[r] := c1;
    if opt[3] then Q[r+1][r+1] := 1; fi;

    j := c1+1;
    while j<=m do
      k := r+1;
      while k<=n and A[r][c1]*A[k][j]=A[k][c1]*A[r][j] do k := k+1; od;
      if k<=n then c2 := j; j := m; fi;
      j := j+1;
    od;
    #Smith with some transforms..
    if opt[1] and (opt[4] or opt[3]) and c2<m then
      N := Gcd(Flat(A{[r..n]}[c2]));
      L := [c1+1..c2-1];
      Append(L,[c2+1..m-1]);
      Add(L,c2);
      for j in L do
	if j=c2 then
	  b:=A[r][c2];a:=A[r][c1];
	  for i in [r+1..n] do
	    if b<>1 then
	      g:=Gcdex(b,A[i][c2]);
	      b:=g.gcd;
	      a:=g.coeff1*a+g.coeff2*A[i][c1];
	    fi; 
	  od;
	  N:=0;
	  for i in [r..n] do  
	    if N<>1 then N:=GcdInt(N,A[i][c1]-QuoInt(A[i][c2],b)*a);fi;
	  od;
	else
	  c := MATINTmgcdex(N,A[r][j],A{[r+1..n]}[j]);
	  b := A[r][j]+c*A{[r+1..n]}[j];
	  a := A[r][c1]+c*A{[r+1..n]}[c1];
	fi;
	t := MATINTmgcdex(N,a,[b])[1];
	tmp := A[r][c1]+t*A[r][j];
	if tmp=0 or tmp*A[k][c2]=(A[k][c1]+t*A[k][j])*A[r][c2] then
	  t := t+1+MATINTmgcdex(N,a+t*b+b,[b])[1];
	fi;
	if t>0 then
	  for i in [1..n] do A[i][c1] := A[i][c1]+t*A[i][j]; od;
	  if opt[4] then B[j][c1] := B[j][c1]+t; fi;
	fi;
      od;
      if A[r][c1]*A[k][c1+1]=A[k][c1]*A[r][c1+1] then
	for i in [1..n] do A[i][c1+1] := A[i][c1+1]+A[i][c2]; od;
	if opt[4] then B[c2][c1+1] := 1; fi; 
      fi;
      c2 := c1+1;
    fi;

    c := MATINTmgcdex(AbsInt(A[r][c1]),A[r+1][c1],A{[r+2..n]}[c1]);
    for i in [r+2..n] do 
      if c[i-r-1]<>0 then
	A[r+1]:=A[r+1]+A[i]*c[i-r-1];
	if opt[3] then 
	  C[r+1][i] := c[i-r-1];  
	  Q[r+1]:=Q[r+1]+Q[i]*c[i-r-1]; 
	fi;
      fi;
    od;

    i := r+1;
    while A[r][c1]*A[i][c2]=A[i][c1]*A[r][c2] do i := i+1; od;
    if i>r+1 then
      c := MATINTmgcdex(AbsInt(A[r][c1]),A[r+1][c1]+A[i][c1],[A[i][c1]])[1]+1;
      A[r+1]:=A[r+1]+A[i]*c;
      if opt[3] then 
	C[r+1][i] := C[r+1][i]+c; 
	Q[r+1]:=Q[r+1]+Q[i]*c; 
      fi;
    fi;
    
    g := MATINTbezout(A[r][c1],A[r][c2],A[r+1][c1],A[r+1][c2]);
    sig:=sig*SignInt(A[r][c1]*A[r+1][c2]-A[r][c2]*A[r+1][c1]);
    A{[r,r+1]} := [[g.coeff1,g.coeff2],[g.coeff3,g.coeff4]]*A{[r,r+1]};
    if opt[3] then 
      Q{[r,r+1]} := [[g.coeff1,g.coeff2],[g.coeff3,g.coeff4]]*Q{[r,r+1]};
    fi;

    for i in [r+2..n] do
      q := QuoInt(A[i][c1],A[r][c1]);
      A[i]:=A[i]+A[r]*-q;
      if opt[3] then Q[i]:=Q[i]+Q[r]*-q; fi;
      q := QuoInt(A[i][c2],A[r+1][c2]);
      A[i]:=A[i]+A[r+1]*-q;
      if opt[3] then Q[i]:=Q[i]+Q[r+1]*-q; fi;
    od;

  od; 
  rp[r+1] := m;
  if n=m and r+1<n then sig:=0;fi;

  #smith w/ NO transforms - farm the work out...
  if opt[1] and not (opt[3] or opt[4]) then
    #R:=rec(normal:=SNFofREF(A{[2..n-1]}{[2..m-1]}),rank:=r-1);
    for i in [2..n-1] do
      A[i-1]:=A[i]{[2..m-1]};
    od;
    Unbind(A[n-1]);
    Unbind(A[n]);
    R:=rec(normal:=SNFofREF(A,opt[5]),rank:=r-1);
    if n=m then R. signdet:=sig;fi;
    return R;
  fi;

  # hermite or (smith w/ column transforms)
  if (not opt[1] and opt[2]) or (opt[1] and opt[4]) then
    for i in [r, r-1 .. 1] do
      for j in [i+1 .. r+1] do
	q := QuoInt(A[i][rp[j]]-(A[i][rp[j]] mod A[j][rp[j]]),A[j][rp[j]]);
	A[i]:=A[i]+A[j]*-q;
	if opt[3] then Q[i]:=Q[i]+Q[j]*-q; fi;
      od;
      if opt[1] and i<r then
	for j in [i+1..m] do 
	  q := QuoInt(A[i][j],A[i][i]);
	  for k in [1..i] do A[k][j] := A[k][j]-q*A[k][i]; od;
	  if opt[4] then P[i][j] := -q; fi;      
	od;
      fi;
    od;
  fi;

  #Smith w/ row but not col transforms
  if opt[1] and opt[3] and not opt[4] then
    for i in [1..r-1] do
      t := A[i][i];
      A[i] := List([1..m],x->0);
      A[i][i] := t;
    od;
    for j in [r+1..m-1] do
      A[r][r] := GcdInt(A[r][r],A[r][j]);
      A[r][j] := 0;
    od;
  fi;

  #smith w/ col transforms
  if opt[1] and opt[4] and r<m-1 then
    c := MATINTmgcdex(A[r][r],A[r][r+1],A[r]{[r+2..m-1]});
    for j in [r+2..m-1] do
      A[r][r+1] := A[r][r+1]+c[j-r-1]*A[r][j];
      B[j][r+1] := c[j-r-1];
      for i in [1..r] do P[i][r+1] := P[i][r+1]+c[j-r-1]*P[i][j]; od;
    od;
    P[r+1] := List([1..m],x->0);
    P[r+1][r+1] := 1;
    g := Gcdex(A[r][r],A[r][r+1]);
    A[r][r] := g.gcd;
    A[r][r+1] := 0;
    for i in [1..r+1] do
      t := P[i][r];
      P[i][r] := P[i][r]*g.coeff1+P[i][r+1]*g.coeff2;
      P[i][r+1] := t*g.coeff3+P[i][r+1]*g.coeff4;
    od;
    for j in [r+2..m-1] do  
      q := QuoInt(A[r][j],A[r][r]);
      for i in [1..r+1] do P[i][j] := P[i][j]-q*P[i][r]; od;
      A[r][j] := 0;
    od;
    for i in [r+2..m-1] do
      P[i] := List([1..m],x->0);
      P[i][i] := 1;
    od;
  fi;

  #row transforms finisher
  if opt[3] then for i in [r+2..n] do Q[i][i]:= 1; od; fi;

  for i in [2..n-1] do
    A[i-1]:=A[i]{[2..m-1]};
  od;
  Unbind(A[n-1]);
  Unbind(A[n]);
  R:=rec(normal:=A);

  if opt[3] then 
    R.rowC:=C{[2..n-1]}{[2..n-1]}; 
    R.rowQ:=Q{[2..n-1]}{[2..n-1]}; 
  fi;

  if opt[1] and opt[4] then
    R.colC:=B{[2..m-1]}{[2..m-1]}; 
    R.colQ:=P{[2..m-1]}{[2..m-1]}; 
  fi;

  R.rank:=r-1;
  if n=m then R.signdet:=sig;fi;
  return R;

end;

#############################################################################
##
#F  NormalFormIntMat(<mat>,<options>)
##
NormalFormIntMat:=function(mat,options)
  local r,opt;
  r:=DoNFIM(mat,options);
  opt := BITLISTS_NFIM[options+1];

  if opt[3] then r.rowtrans:=r.rowQ*r.rowC; fi;

  if opt[1] and opt[4] then  r.coltrans:=r.colC*r.colQ; fi;
  return r;
end;

TriangulizedIntegerMat:=mat->DoNFIM(mat,0).normal;
TriangulizedIntegerMatTransform:=mat->NormalFormIntMat(mat,4);
TriangulizeIntegerMat:=mat->DoNFIM(mat,16);
HermiteNormalFormIntegerMat:=mat->DoNFIM(mat,2).normal;
HermiteNormalFormIntegerMatTransform:=mat->NormalFormIntMat(mat,6);
SmithNormalFormIntegerMat:=mat->DoNFIM(mat,1).normal;
SmithNormalFormIntegerMatTransforms:=mat->NormalFormIntMat(mat,13);
DiagonalizeIntMat:=mat->DoNFIM(mat,17);
BaseIntMat:=function(mat)local norm;
  norm:=NormalFormIntMat(mat,2);
  return norm.normal{[1..norm.rank]};
end;

#############################################################################
##
#M  BaseIntersectionIntMats(<m1>,<m2>)
##
BaseIntersectionIntMats:=function(M1, M2 )
local M, Q, r, T;
  M := Concatenation( M1, M2 );
  r := NormalFormIntMat( M, 4 );
  T := r.rowtrans{[r.rank+1..Length(M)]}{[1..Length(M1)]};
  if not Length(T)=0 then T := T * M1; fi;
  return BaseIntMat( T );
end;

#############################################################################
##
#M  ComplementIntMat(<m1>,<m2>)
##
ComplementIntMat:=function( full,sub ) local F, S, M, r, T, R;
  F := BaseIntMat( full );
  if Length(sub)=0 or Is_Zero( sub ) then 
    return rec( complement := F, sub := [], moduli := [] );
  fi;
  S := BaseIntersectionIntMats( F, sub );
  if S <> BaseIntMat( sub ) then
    Error( "sub must be submodule of full" );
  fi;
  # find T such that BaseIntMat(T*F) = S
  M := Concatenation( F, S );
  T := NormalFormIntMat( M, 4 ).rowtrans^-1;
  T := T{[Length(F)+1..Length(T)]}{[1..Length(F)]};

  # r.rowtrans * T * F = r.normal * r.coltrans^-1 * F
  r := NormalFormIntMat( T, 13 );
  M := r.coltrans^-1 * F;
  R := rec( complement := BaseIntMat( M{[1+r.rank..Length(M)]} ),
	    sub := r.rowtrans * T * F,
	    moduli    := List( [1..r.rank], i -> r.normal[i][i] ) );
  return R;
end;

NullspaceIntMat:=function(mat)local norm, kern;
  norm := NormalFormIntMat( mat, 4 );
  kern := norm.rowtrans{[norm.rank+1..Length(mat)]};
  return BaseIntMat( kern );
end;

#############################################################################
##
#M  SolutionIntMat(<mat>,<vec>)
##
SolutionIntMat:= function( mat,v ) local norm, rs, t, M, r;
  if Is_Zero(mat) then
    if Is_Zero(v) then return [1..Length(mat)]*0;
    else return false;
    fi;
  fi;
  norm := NormalFormIntMat( mat, 4 );
  t := norm.rowtrans;
  rs :=  norm.normal{[1..norm.rank]};
  M := Concatenation( rs, [v] );
  r := NormalFormIntMat( M, 4 );
  if r.rank = Length(r.normal) or
     r.rowtrans[Length(M)][Length(M)] <> 1 then
    return false;
  fi;
  return -r.rowtrans[Length(M)]{[1..r.rank]} * t{[1..r.rank]};
end;

#############################################################################
##
#M  SolutionNullspaceIntMat(<mat>,<vec>)
##
SolutionNullspaceIntMat:=function( mat,v )
local norm, rs, t, M, r, kern, len;
  if Is_Zero(mat) then
    len := Length(mat);
    if Is_Zero(v) then
      return [[1..len]*0, IdentityMat(len)];
    else
      return [false, IdentityMat(len)];
    fi;
  fi;
  norm := NormalFormIntMat( mat, 4 );
  kern := norm.rowtrans{[norm.rank+1..Length(mat)]};
  kern := BaseIntMat( kern );
  t := norm.rowtrans;
  rs :=  norm.normal{[1..norm.rank]};
  M := Concatenation( rs, [v] );
  r := NormalFormIntMat( M, 4 );
  if r.rank = Length(r.normal) or
     r.rowtrans[Length(M)][Length(M)] <> 1 then
    return [false,kern];
  fi;
  return [-r.rowtrans[Length(M)]{[1..r.rank]} * t{[1..r.rank]}, kern];
end;


#############################################################################
##
#F  DeterminantIntMat(<mat>)
##
DeterminantIntMat:=function(mat)
local sig, n, m, A, r, c2, c1, j, k, c, i, g, q;

  sig:=1;

  #Embed mat in 2 larger "id" matrix
  n := Length(mat)+2;
  # Crossover point roughly 20x20 matrices, so farm the work if smaller..
  if n<22 then return DeterminantMat(mat);fi;
  m := Length(mat[1])+2;
  
  if not n=m then 
    Error( "DeterminantIntMat: <mat> must be a square matrix" );
  fi;

  A := [List([1..m],x->0)];
  for i in [2..n-1] do
    A[i] := [0];
    Append(A[i],mat[i-1]);
    A[i][m] := 0;
  od;
  A[n] := List([1..m],x->0);
  A[1][1] := 1;      A[n][m] := 1;

  r := 0;    c2 := 1;
  while m>c2 do
    r := r+1;
    c1 := c2;

    j := c1+1;
    while j<=m do
      k := r+1;
      while k<=n and A[r][c1]*A[k][j]=A[k][c1]*A[r][j] do k := k+1; od;
      if k<=n then c2 := j; j := m; fi;
      j := j+1;
    od;

    c := MATINTmgcdex(AbsInt(A[r][c1]),A[r+1][c1],A{[r+2..n]}[c1]);
    for i in [r+2..n] do 
      if c[i-r-1]<>0 then
	A[r+1]:=A[r+1]+A[i]*c[i-r-1];
      fi;
    od;

    i := r+1;
    while A[r][c1]*A[i][c2]=A[i][c1]*A[r][c2] do 
    i := i+1; 
    od;

    if i>r+1 then
      c := MATINTmgcdex(AbsInt(A[r][c1]),A[r+1][c1]+A[i][c1],[A[i][c1]])[1]+1;
      A[r+1]:=A[r+1]+A[i]*c;
    fi;
    
    g := MATINTbezout(A[r][c1],A[r][c2],A[r+1][c1],A[r+1][c2]);
    sig:=sig*SignInt(A[r][c1]*A[r+1][c2]-A[r][c2]*A[r+1][c1]);
    if sig=0 then return 0;fi;
    A{[r,r+1]} := [[g.coeff1,g.coeff2],[g.coeff3,g.coeff4]]*A{[r,r+1]};

    for i in [r+2..n] do
      q := QuoInt(A[i][c1],A[r][c1]);
      A[i]:=A[i]+A[r]*-q;
      q := QuoInt(A[i][c2],A[r+1][c2]);
      A[i]:=A[i]+A[r+1]*-q;
    od;
  od; 

  for i in [2..r+1] do
    sig:=sig*A[i][i];
  od;

  return sig;

end;

IntersectionLatticeSubspace:=function(m)local r,i;
  m:=m*Lcm(List(Concatenation(m),Denominator));
  r:= SmithNormalFormIntegerMatTransforms(m);
  for i in [1..Length(r.normal)] do
    if r.normal[i]<>0*r.normal[i] then
      r.normal[i]:=r.normal[i]/Maximum(List(r.normal[i],AbsInt));
    fi;
  od;
  return r.rowtrans^-1*r.normal*r.coltrans;
end;

#############################################################################
##
#F  DiaconisGraham(m,moduli)
#
# m  should be  an integral  matrix with  n columns where n=Length(moduli),
# such  that each  line (with  the i-th  element taken  mod moduli[i]) of m
# represents an element of the group A=Z/moduli[1] x ... x Z/moduli[n], and
# such  that the set of  lines of m generates  A. The moduli should be such
# that moduli[i+1] divides moduli[i] for all i.
#
# The function returns  a record  r with fields
# r.normal          the Diaconis-Graham normal form, a matrix of same shape
#    as  m where either the  first n lines are  the identity matrix and the
#    remaining  lines are  0, or  Length(m)=n and  .normal differs from the
#    identity  matrix only  in the  entry .normal[n][n],  which is prime to
#    moduli[n]
# r.rowtrans        a unimodular matrix such that  
#   r.normal=List(r.rowtrans*m,v->Zip(v,moduli,function(x,y)return x mod y;end))
#
DiaconisGraham:=function(m,moduli)local r,res,i,l,l1,e,n,ZipMod;
  if moduli=[] then return rec(rowtrans:=[],normal:=m);fi;
  if ForAny([1..Length(moduli)-1],i->moduli[i] mod moduli[i+1]<>0) then
    Error("DiaconisGraham(m,moduli): moduli[i+1] should divide moduli[i] for all i");
  fi;
  r:=HermiteNormalFormIntegerMatTransform(m);
  res:=r.rowtrans;
  m:=r.normal;
  n:=Length(moduli);
  if Length(m)>0 and n<>Length(m[1]) then
    Error("DiaconisGraham(m,moduli): moduli and m[1] should have same length");
  fi;
  ZipMod := function(m, moduli)
      return List(m,v->Zip(v,moduli,function(x,y)return x mod y;end));
  end;
  for i in [1..Minimum(n,Length(m)-1)] do
    l:=m[i][i];
    if m[i][i]<>1 then
      if Gcd(m[i][i],moduli[i])<>1 then return false;fi;
      l1:=InverseMod(l,moduli[i]);
      e:=IdentityMat(Length(m));
      e{[i,i+1]}{[i,i+1]}:=[[l1,l1-1],[1-l*l1,l+1-l*l1]];
      res:=e*res;
      m:=ZipMod(e*m,moduli);
    fi;
  od;
  r:=HermiteNormalFormIntegerMatTransform(m);
  m:=ZipMod(r.normal,moduli);
  res:=r.rowtrans*res;
  if Length(m)=n then
    if m[n][n] > QuoInt(moduli[n], 2) then
      m[n][n] := -m[n][n] mod moduli[n];
      res[n] := res[n]*-1;
    fi;
    l:=m[n][n];
    if Gcd(l,moduli[n])<>1 then return false;fi;
    l1:=InverseMod(l,moduli[n]);
    for i in [1..n-1] do
      if m[i][n]<>0 then
        e:=IdentityMat(Length(m));
        e{[i,n]}{[i,n]}:=[[1,-m[i][n]*l1],[0,1]];
        res:=e*res;
        m:=ZipMod(e*m,moduli);
      fi;
    od;
  fi;
  return rec(rowtrans:=res,normal:=m);
end;
