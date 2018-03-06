ReadChv("contr/rouquierblockdata");
Read("RouquierBlocks.g");
# tester vs. Maria rouquierblockdata.
# faire cmp(n) pour n\in[4..37]
#-------------------------------------------------------------------------------
#Generic1SchurElementOps:=OperationsRecord("Generic1SchurElementOps");
#
#Generic1SchurElementOps.Print:=function(s)local pmon;
#  pmon:=function(v)local pos,neg,mon,res;
#    pos:=List(v,function(x)if x>0 then return x;else return 0;fi;end);
#    neg:=pos-v;
#    mon:=function(v)local den,res;
#      den:=Lcm(List(v,Denominator));
#      res:=Product(Zip(v*den,s.vars,function(a,b)return b^a;end));
#      if den=1 then return Format(res);
#      else return SPrint("(",res,")^(1/",den,")");fi;
#    end;
#    res:=String(mon(pos));
#    if neg<>0*neg then PrintToString(res,"/",mon(neg));fi;
#    return res;
#  end;
#  Print(FormatCoefficient(s.coeff,pmon(s.mon),rec()),"*");
#  Print(Join(List(s.vcyc,v->SPrint(v.pol,"(",
#     FormatCoefficient(v.coeff,pmon(v.mon),rec()),")")),""));
#end;
#
#Generic1SchurElementOps.Value:=function(s,v)local mon;
#  mon:=v->Product([1..Length(v)],i->s.vars[i]^v[i]);
#  s:=s.coeff*mon(s.mon)*Product(s.vcyc,x->Value(x.pol,x.coeff*mon(x.mon)));
#  return Value(s,v);
#end;
 
permbl:=function(a,p)
  a:=List(a,x->OnTuples(x,p));
  for p in a do Sort(p);od;
  Sort(a);
  return a;
end;

cmpbl:=function(a,b)local dif,res,n,vdif,i,a1,difif,m,pairs;
  a:=Copy(a);n:=Maximum(Flat(a));
  dif:=function(v1,v2)return Number(TransposedMat(List([v1,v2],Concatenation)),
    x->x[1]<>x[2]);end;
  vdif:=function(v1,v2)return Sum(Zip(v1,v2,dif));end;
  difif:=function(perm)local a1;
     a1:=Copy(a);a1:=List(a1,u->permbl(u,perm));return vdif(a1,b);
  end;
  res:=();
  i:=Filtered([1..Length(a)],i->a[i]<>b[i]);
  if i=[] then return res;fi;
  if Collected(List(a[i[1]],Length))<>Collected(List(b[i[1]],Length)) then
    a1:=Intersection(a[i[1]],b[i[1]]);
    a:=Difference(a[i[1]],a1);b:=Difference(b[i[1]],a1);
    Print("different blocks",a,b,"\n");
    return false;
  fi;
  while vdif(a,b)>0 do
    Print("vdif=",vdif(a,b),"=>");
    pairs:=Zip(a,b,function(v1,v2)
      return TransposedMat(List([v1,v2],Concatenation));end);
    pairs:=Set(Concatenation(pairs));
    pairs:=Filtered(pairs,x->x[1]<>x[2]);
    pairs:=List(pairs,x->(x[1],x[2]));
#   pairs:=List(Combinations([1..n],2),i->(i[1],i[2]));
    m:=Minimum(List(pairs,difif));
    m:=First(pairs,i->difif(i)=m);
    a1:=Copy(a);a1:=List(a1,u->permbl(u,m));
    Print(vdif(a1,b),".\c");
    if vdif(a1,b)<vdif(a,b) then a:=a1;res:=res*m;fi;
  od;
  Print(res,"\n");
  return res;
end;

cmp:=function(n)local nm,a,b,pp,d,p,rouq,maria,W,i,bad,perm,a1,b1,int;
  W:=ComplexReflectionGroup(n);
  rouq:=RouquierBlockData(W);
  rouq:=List(rouq,x->rec(cond:=x[1],block:=x[2]));
  rouq[1].cond:=[];
  rouq:=Concatenation([rouq[1]],
    Filtered(rouq{[2..Length(rouq)]},x->x.block<>rouq[1].block));
  nm:=SPrint("G",n);
  maria:=BlockResults.(nm);SortBy(maria,x->x.cond);
  # first exclude line with different blocks
  a:=Copy(maria);b:=Copy(rouq);
  if Length(a)<>Length(b) then Error("length");fi;
  if List(a,x->x.cond)<>List(b,x->x.cond) then Error("cond");fi;
  bad:=[];
  for i in [1..Length(a)] do
    if Collected(List(a[i].block,Length))<>
       Collected(List(a[i].block,Length))
    then Add(bad,i);
      a1:=a[i].block;b1:=b[i].block;
      int:=Intersection(a1,b1);
      a1:=Difference(a1,int);b1:=Difference(b1,int);
      Print("different blocks: p=",p," ",a1,b1,"\n");
    fi;
  od;
  a:=Drop(a,bad); b:=Drop(b,bad);
  # then find perm
  cmpbl(List(a,x->x.block),List(a,x->x.block));
end;
