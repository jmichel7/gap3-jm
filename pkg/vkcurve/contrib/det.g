# variation on DeterminantMat trying to handle more matrices
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
