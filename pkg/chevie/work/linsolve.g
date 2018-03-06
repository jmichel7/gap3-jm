# linsys represent the status of partially solved linear equations
# for variables. One wants to solve for x_i a system of equations 
#  $\{\sum_i=1^n a_{i,j} x_i=c_j\}_j$
#
#  Usage:
#
#  to start:
#  M:=LinSys(n); # n=number of unknowns
#
#  the value n is then stored in the field M.total
#
#  to add an equation:
#
#  Makerel(M,v,a)
#   each elm of v is [index in [1..M.total],coeff]
#   specifies the equation $a+\sum_{p\in v} p[2] x_{p[1]}=0$
#  
#  Knowledge  is deduced  about the  x_i by  applying TriangulizeMat.  This is
#  applied  automatically when M.freq  (default value 100)  new equations have
#  been added. It can be applied explicitely by calling
#  applynew(M);
#
# If there is a chance that contradictory equations are given one should check
# the  return value of Makerel and applynew: it is true for Makerel and equal
# to the number of independent relations for makerel normally, but is false if
# a  call  to  TriangulizeMat  has  occured  and  has  detected  contradictory
# relations.
#
#  Information can be extracted form the system by calling
#   Values(M)
#  
#  This returns a vector of length M.total, containing in the ith slot the
#  value of x_i if known, or otherwise a linear form in some Mvp's expressing
#  the relationship with other x_j.
#
#  Here is a complete example of use:
#
#  gap> M:=LinSys(6);   # there are 6 unknowns
#  known:0 unknown:6
#  gap> Makerel(M,[[1,1],[2,1]],3); # specify x_1+x_2=3
#  true
#  gap> Makerel(M,[[2,1],[3,1]],4); # specify x_2+x_3=4
#  true
#  gap> Makerel(M,[[3,1],[4,1]],5); # specify x_3+x_4=5
#  true
#  gap> Makerel(M,[[4,1],[5,1]],6); # specify x_4+x_5=6
#  true
#  gap> Makerel(M,[[6,1]],0);       # specify x_6=0
#  true
#  gap> applynew(M);      # one variable is completely known.
#  known:1 unknown:1=5-4  # for the others we have 4 equations for 5 unknowns
#  5                      # 5 independent relations exist
#  gap> Values(M);
#  [ 2+a, -5-a, 1+a, -6-a, a, 0 ]  # known information on M
#  gap> Values(M,["x"]);
#  [ 2+x, -5-x, 1+x, -6-x, x, 0 ]  # a second argument tells variables names
                                   # to use

LinSysOps:=OperationsRecord("LinSysOps");
LinSysOps.Treshold:=2;

# At start only the total number n of x_i must be specified
LinSys:=n->rec(newrel:=[],known:=[],values:=[],relations:=[],total:=n,
    freq:=100,operations:=LinSysOps);

# value of x_i indicated by M
LinSysOps.Value:=function(M,i)local p;
  p:=Position(M.known,i);
  if p=false then return false;
  else return M.values[p];
  fi;
end;

LinSysOps.String:=function(M)local res;
  res:=SPrint("known:",Length(M.known),
              " unknown:",M.total-Length(M.known)-Length(M.relations));
  if Length(M.relations)>0 then 
    PrintToString(res,"=",Length(M.relations[1])-1,"-",Length(M.relations));
  fi;
  return res;
end;

LinSysOps.Print:=function(M) Print(String(M));end;

# To save space we apply each new set of relations as soon as possible.
# We represent our current knowledge as a record |M| with 3 fields:
# \begin{itemize}
#  \item |M.known| is the indices i of the variables $x_i$ that we know already.
#  \item |M.values| is the values of the variables $x_i$ that we know already.
#  \item |M.relations| is a set of vectors representing the relations on
#     the remaining variables.
# \end{itemize}
# The next  function takes  a bunch  of new  relations and  modifies |M|
# accordingly.  It uses  the routine  |BaseMat| to  find all  completely
# known basis  vectors resulting from |M.relations|,  and then supresses
# the corresponding  columns from |M.relations|, adding  entries instead
# to |M.known| and |M.values|.
# It returns |true| iff at the end of the process all values are known.

apply:=function(M,newrels)local compl,i,newval,newind,IsZero,tmp;
  IsZero:=x->x=0*x;
  if Length(newrels)=0 then return true;fi;
  compl:=Filtered([1..Length(newrels[1])-1],i->not i in M.known);
  if Length(M.relations)>0 and Length(M.relations[1])>LinSysOps.Treshold then
    InfoChevie(Length(newrels)," new relations\c");
  fi;
  if Length(M.known)>0 then
    newrels:=List(newrels,
                x->Concatenation(x{compl},[x[Length(x)]+x{M.known}*M.values]));
  fi;
  tmp:=BaseMat(newrels);
  if Length(M.relations)>0 and Length(M.relations[1])>LinSysOps.Treshold then
    InfoChevie("(",Length(tmp)," independent) \c");
  fi;
  if IsBound(M.stopindep) and Length(tmp)=M.stopindep then return false;fi;
  newrels:=tmp;
  Append(M.relations,newrels);
  M.relations:=BaseMat(M.relations);
  if Length(M.relations)>0 and Length(M.relations)>0 and 
     PositionProperty(M.relations[Length(M.relations)],
                    y->not IsZero(y))=Length(M.relations[1])
  then Print("!!!contradictory relations");
    if IsBound(M.error) then Error();fi;
    return false;fi;
  i:=List(M.relations,x->Number(x{[1..Length(x)-1]},y->not IsZero(y))=1);
  newval:=ListBlist(M.relations,i);
  M.relations:=ListBlist(M.relations,List(i,x->not x));
  newind:=List(newval,x->PositionProperty(x,y->not IsZero(y)));
  M.relations:=List(M.relations,x->Concatenation(
	x{Filtered([1..Length(x)-1],y->not y in newind)},[x[Length(x)]]));
  Append(M.known,compl{newind});
  Append(M.values,List(newval,x->-x[Length(x)]));
  SortParallel(M.known,M.values);
  if Length(M.relations)>0 and Length(M.relations[1])>LinSysOps.Treshold then 
    InfoChevie(String(M),"\n");
  fi;
  return true;
end;

# checkMagrees(M,values [,transl]) check M agrees with values  
# if given transl=translation of coordinates (like ij) in error message
checkMagrees:=function(arg)local M,values,transl,i,unknown,r;
  M:=arg[1];values:=arg[2];
  if Length(arg)=3 then transl:=arg[3];else transl:=x->x;fi;
  for i in [1..Length(M.known)] do
    if M.values[i]<>values[M.known[i]] then 
      Error("f",transl(i)," computed  as ",M.values[i],
                          " instead of ",values[M.known[i]]);fi;
  od;
  if Length(M.relations)>0 then
    unknown:=Difference([1..M.total],M.known);
    for r in M.relations do
      if Sum([1..Length(r)-1],i->r[i]*values[unknown[i]])+r[Length(r)]<>0 then
	Error("relation ",Filtered(List([1..Length(unknown)],
	  x->[transl(unknown[x]),r[x]]),x->not x[2]=0*x[2]),"+",
	  r[Length(r)],"=0 not satisfied");
      fi;
    od;
  fi;
  Print("known values and relations agree!\n");
end;

applynew:=function(M)local r;
  if not apply(M,M.newrel) then return false;fi;
  M.newrel:=[];return Length(M.known)+Length(M.relations);
end;

# each elm of v is [index as in M.known,coeff]
# rhs is -right hand side of relation
# so this descrives the equation sum_{e in v} e[2] x_e[1] +rhs=0
Makerel:=function(M,v,rhs) local rel,i,IsZero;
  IsZero:=x->x=0*x;
  v:=Filtered(v,x->x[2]<>0);
  if ForAll(v,x->x[1] in M.known) then # check v, rhs for correctness
   if not IsZero(Sum(v,x->M.values[Position(M.known,x[1])]*x[2])+rhs) then 
     Print("!!!!!! new relation ",v,"+",rhs,"=0 is contradictory\n");
     if IsBound(M.error) then Error();fi;
     return false;
   fi;
  else rel:=[1..M.total]*0*rhs;
    for i in v do rel[i[1]]:=rel[i[1]]+i[2];od;
    Add(rel,rhs);Add(M.newrel,rel);
    if Length(M.newrel)>M.freq then 
      if false=applynew(M) then return false;fi;
    fi;
  fi;
  return true;
end;

# get values in M, giving linear combs of Mvp's for unknowns
# Values(M[,varnames])
Values:=function(arg)local i,p,unknown,q,res,M,NotZero,unknownpos;
  NotZero:=x->x<>0*x;
  M:=arg[1];
  unknownpos:=Difference([1..M.total],M.known);
  res:=[]; res{M.known}:=M.values;
  if not IsBound(M.varnames) then 
    M.varnames:=List([1..M.total],i->SPrint("x",i));
  fi;
  if Length(M.relations)=0 then 
    res{unknownpos}:=M.varnames{unknownpos};return res;
  fi;
  unknown:=[1..Length(M.relations[1])-1];
  q:=Difference(unknown,List(M.relations,p->PositionProperty(p,NotZero)));
# Print("unknownpos=",unknownpos," unknown=",unknown," q=",q,"\n");
  unknown:=List(unknown,function(p)local i;
    i:=Position(q,p);
    if i<>false then return Mvp(M.varnames[unknownpos[q[i]]]);
    else return 0*Mvp("a");fi;end);
  Add(unknown,Mvp("a")^0);
# Print("M.varnames=",M.varnames," unknown=",unknown,"\n");
  for p in M.relations do 
    i:=PositionProperty(p,NotZero);
    if NotZero(p[i]-1) then Error("relation does not start with 1");fi;
    p:=ShallowCopy(p);p[i]:=p[i]*0;unknown[i]:=-p*unknown;
  od;
  res{unknownpos}:=unknown{[1..Length(unknown)-1]};
  return res;
end;

# return rec if linear form, false otherwise
MvpLinForm:=function(p)local c,res,i;
  if Length(p.coeff)<>Length(p.elm) then Error();fi;
  c:=PositionProperty(p.elm,x->Length(x.coeff)=0);
  res:=rec();
  if c<>false then
    res.c:=p.coeff[c];
    i:=Filtered([1..Length(p.coeff)],i->i<>c);
  else if Length(p.coeff)>0 then res.c:=0*p.coeff[1]; else res.c:=0;fi;
    i:=[1..Length(p.coeff)];
  fi;
  if ForAny(p.elm{i},x->Length(x.coeff)<>1 or x.coeff[1]<>1) then
    Error("not a linear form:",p);
    return false;fi;
  res.vars:=List(p.elm{i},x->x.elm[1]);
  res.coeff:=p.coeff{i};
  return res;
end;

# make a relation from an Mvp (obtained by Values(M) e.g.)
Addrel:=function(M,p)
  p:=MvpLinForm(p);
  return Makerel(M,List([1..Length(p.coeff)],
    i->[Position(M.varnames,p.vars[i]),p.coeff[i]]),p.c);
end;

# below ll is a bunch of linear forms in some Mvps.
# it is made into a  LinSys to be solved; track is made of the
# variable names
MvpSolve:=function(arg)local ll,varnames,l,M,v;
  # collect varnames
  ll:=List(arg[1],MvpLinForm);
  varnames:=Union(List(ll,x->x.vars));
  M:=LinSys(Length(varnames));M.varnames:=varnames;
  if Length(arg)=2 then M.stopindep:=arg[2];fi;
  for l in ll do
    v:=List([1..Length(l.coeff)],i->[Position(varnames,l.vars[i]),l.coeff[i]]);
    Makerel(M,v,l.c);
  od;
  if false=applynew(M) then return false;fi;
# return M;
  return Concatenation(Filtered(TransposedMat([varnames,Values(M)]),x->
    x[2]<>Mvp(x[1])));
end;
