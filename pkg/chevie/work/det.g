# variations on DeterminantMat trying to handle more matrices

myDet:=function(m)local n,f,i,j,k,det,compl,v;
  compl:=function(m,i,j)local v;v:=[1..Length(m)];
    return m{Difference(v,[i])}{Difference(v,[j])};end;
  n:=Length(m);
  if n=1 then return m[1][1];
  elif n=2 then return m[1][1]*m[2][2]-m[1][2]*m[2][1];
  fi;
  det:=1;
  repeat
    i:=PositionProperty(m,v->Number(v,x->x<>0*x)<=1);
    if i<>false then 
      if Number(m[i],x->x<>0*x)=0 then return 0;fi;
      j:=PositionProperty(m[i],x->x<>0*x);
      det:=det*(-1)^(i+j)*m[i][j];
      m:=compl(m,i,j);
      if Length(m)>71 then Print("=>",Length(m),"\c");fi;
    fi;
  until i=false;
  if m=[] then return det;fi;
  v:=[1..Length(m)];
  for j in v do
    i:=PositionProperty(m{v}[j],x->IsRec(x) and Length(x.elm)=1);
    if i<>false then Print(Length(m),":",[j,i],"\n");
      m:=ShallowCopy(m);
      f:=m[i][j];
      for k in Difference(v,[i]) do m[k]:=m[k]-m[k][j]*f^-1*m[i];od;
      return det*(-1)^(i+j)*m[i][j]*myDet(compl(m,i,j));
    fi;
  od;
  v:=List(m,x->Number(x,y->y<>0));
  j:=Position(v,Minimum(v));
  if Length(m)>71 then Print("\n",Length(m),":",v[j],"\c");fi;
  if v[j]>5 then return det*DeterminantMat(m);fi;
  return det*Sum([1..Length(m)],function(i)
#   if Length(m) mod 5=0 then Print(Length(m),":",i," \c");fi;
    if m[i][j]=0*m[i][j] then return 0;
    else return (-1)^(i+1)*m[i][j]*myDet(compl(m,i,j));
    fi;end);
end;

myDet2:=function(m)local n,f,i,j,k,det,compl,v,mvp;
  compl:=function(m,i,j)local v;v:=[1..Length(m)];
    return m{Difference(v,[i])}{Difference(v,[j])};end;
  mvp:=IsMvp(m[1][1]);
  n:=Length(m);
  if n=1 then return m[1][1];
  elif n=2 then return m[1][1]*m[2][2]-m[1][2]*m[2][1];
  fi;
  det:=1;
  repeat
    i:=PositionProperty(m,v->Number(v,x->x<>0*x)<=1);
    if i<>false then 
      if Number(m[i],x->x<>0*x)=0 then return 0;fi;
      j:=PositionProperty(m[i],x->x<>0*x);
      det:=det*(-1)^(i+j)*m[i][j];
      m:=compl(m,i,j);
      Print("=>",Length(m),"\c");
    fi;
  until i=false;
  repeat
    if m=[] then return det;fi;
    v:=List(m,x->Number(x,y->y<>0));
    i:=Position(v,Minimum(v));
    if mvp then
       j:=List(m[i],function(x)
         if x=0 then return 1000;
	 elif Mvp(x)=false then return 900; 
	 else return Length(Mvp(x).elm);fi;end);
       j:=Position(j,Minimum(j));
    else j:=PositionProperty(m[i],x->x<>0*x);
    fi;
    Print("\n",Length(m),":",v[i]," \c");
    Print(m[i][j],"\c");
    for k in Difference([1..Length(m)],[i]) do
      m[k]:=m[k]-m[k][j]/m[i][j]*m[i];od;
    det:=det*(-1)^(i+j)*m[i][j];
    m:=compl(m,i,j);
  until m=[];
  return det;
end;

# next best thing to be able to invert m
myCoF:=function(m)local v,j,res,i;
  if Length(m)=1 then return [[1]];fi;
  v:=[1..Length(m)];
  res:=m*0;
  for i in v do for j in v do Print("i=",i," j=",j,"\n");
    res[i][j]:=(-1)^(i+j)*myDet(m{Filtered(v,x->x<>i)}{Filtered(v,x->x<>j)});
  od;od;
  return TransposedMat(res);
end;

# determinant of a matrix of Polynomials(Rationals)
PolynomialDet:=function(m)local n,v,i,j,k,M,l,d,dircoeff,dnorm,norm,coeff;
  n:=Length(m);
  if Length(m)=1 then return m[1][1];
  elif Length(m)=2 then return m[1][1]*m[2][2]-m[1][2]*m[2][1];
  fi;
  dircoeff:=function(p)
    if Length(p.coefficients)=0 then return 0;fi;
    return p.coefficients[Length(p.coefficients)];
  end;
  dnorm:=function(p)if Length(p.coefficients)=0 then return 1;fi;
    return Lcm(List(p.coefficients,Denominator));
  end;
  norm:=function(p)if Length(p.coefficients)=0 then return 0;fi;
    return Gcd(p.coefficients);
  end;
  m:=ShallowCopy(m);
  coeff:=1;
  repeat
    l:=List(m,function(x)local d;
      d:=Lcm(List(x,dnorm));return d/Gcd(List(x*d,norm));
    end);
    coeff:=coeff/Product(l);
#   Print("coeff=",coeff,"\n");
    for i in [1..n] do m[i]:=m[i]*l[i];od;
    M:=Maximum(List(m,x->Maximum(List(x,y->Length(y.coefficients)))));
    M:=List(m,x->MinPos(List(x,function(y)local l;l:=Length(y.coefficients);
      if l=0 then return [M+1,0,0];fi;
      return [l,AbsInt(y.valuation),
        Length(String(Maximum(List(y.coefficients*dnorm(y),AbsInt))))];
    end)));
    i:=MinPos(List(M,x->x[1]))[2];M:=M[i];j:=M[2];M:=M[1];
    v:=dircoeff(m[i][j]);
    for k in Filtered([1..n],j->j<>i) do
      l:=dircoeff(m[k][j])/v;
      if l<>0 then
	if M[1]=1 then m[k]:=m[k]-EuclideanQuotient(m[k][j],m[i][j])*m[i];
	else d:=Degree(m[k][j])-Degree(m[i][j]);
	     m[k]:=m[k]-List(m[i],function(p)
		p:=ShallowCopy(p);
		p.coefficients:=p.coefficients*l;
		p.valuation:=p.valuation+d;
		return p;
		end);
	fi;
      fi;
    od;
    Print(Length(m),":",M," \c");
  until M[1]=1;
  return coeff*(-1)^(i+j)*m[i][j]*PolynomialDet(
    m{Concatenation([1..i-1],[i+1..n])}{Concatenation([1..j-1],[j+1..n])});
end;
