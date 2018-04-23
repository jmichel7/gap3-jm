#############################################################################
##
#A  check.g           CHEVIE library                            Jean Michel
#Y  Copyright (C) august 2017 -  University  Paris Diderot.
#
# This defines utility functions for comparing various data
# 
CHEVIE.Check:=rec();

# EqLists(a,b [,na,nb])
# compares lists a and b (whose descriptions are strings na and nb)
CHEVIE.Check.EqLists:=function(arg)local a,b,na,j,t,pa,p,nb;
  a:=arg[1];b:=arg[2];
  if a=b then return;fi;
  if IsBound(arg[3]) then na:=arg[3];else na:="a";fi;
  if IsBound(arg[4]) then nb:=arg[4];else nb:="b";fi;
# ChevieErr(na,"<>",nb,"\n");
  if Length(a)<>Length(b) then 
    ChevieErr("Length(",na,")=",Length(a),
      " while Length(",nb,")=",Length(b),"\n");
  fi;
  if -a=b then ChevieErr(na,"=-",nb,"\n");return;fi;
  pa:=PermListList(a,b);
  if pa<>false then ChevieErr(na,"=Permuted(",nb,"),",pa,"\n");return;fi;
  pa:=PermListList(a,-b);
  if pa<>false then ChevieErr(na,"=-Permuted(",nb,"),",pa,"\n");return;fi;
  p:=j->SPrint(na,"[",j,"]");
  for j in [1..Length(a)] do
    if a[j]<>b[j] then
      if a[j]=-b[j] then ChevieErr(p(j),"=-",nb,"[",j,"]\n");
      else t:=PositionsProperty(b,k->k=a[j]);
        if Length(t)>0 then ChevieErr(p(j)," found at ",t,"\n");
        else t:=PositionsProperty(b,k->k=-a[j]);
          if Length(t)>0 then ChevieErr(p(j)," found at -",t,"\n");
          else ChevieErr(p(j)," not found\n");
          fi;
        fi;
      fi;
    fi;
  od;
end;

# compare two objects recursively
# CHEVIE.Check.EqObj(a,b[,option record[,history]])
# options: 
#   na:=description of a
#   nb:=description of b
#   depth:=how deep recursion level to go [default 20]
#   diffs:=how many diffs to report [default 5]
#
# return 0=unknown 1=equal 2=difer
CHEVIE.Check.EqObj:=function(arg)local a,b,d,fa,fb,nopt,depth,pr,diffs,
  compared,res,history,nh,p,opt,nbdiffs;
  a:=arg[1];b:=arg[2];arg:=arg{[3..Length(arg)]};
  if Length(arg)>0 then opt:=arg[1];arg:=arg{[2..Length(arg)]};
  else opt:=rec();fi;
  if not IsBound(opt.na) then opt.na:="a";fi;
  if not IsBound(opt.nb) then opt.nb:="b";fi;
  pr:=function(arg)if depth<0 then return;fi;
    Cut(Replace(Concatenation(List(arg,Format)),"{a}",opt.na,"{b}",opt.nb),
      rec(after:=", ="));
  end;
  if Length(arg)>0 then history:=arg[1]; else history:=[[],[]];fi;
  nh:=function()return[Concatenation(history[1],[a]),
                       Concatenation(history[2],[b])];end;
  if IsBound(opt.depth) then depth:=opt.depth;else depth:=20;fi;
  if IsBound(opt.diffs) then nbdiffs:=opt.diffs-1;else nbdiffs:=5;fi;
  if depth<=-5 then return 0;fi;
  if a=b then return 1;fi;
  fa:=PositionProperty(history[1],x->IsIdentical(x,a));
  fb:=PositionProperty(history[2],x->IsIdentical(x,b));
  if fa=fb and fa<>false then pr("{a} already seen");return 1;fi;
  if IsRec(a) then
    if not IsRec(b) then pr("{a} is a record but in {b} is not");return 2;fi;
    diffs:=0;compared:=true;
    fa:=RecFields(a);fb:=RecFields(b);
    d:=Difference(fa,fb);
    if d<>[] then 
      if Length(d)=1 then p:=" is a field";else p:=" are fields";fi;
      pr(Join(d),p," in {a} but not in {b}");diffs:=diffs+1;fi;
    d:=Difference(fb,fa);
    if d<>[] then 
      if Length(d)=1 then p:=" is a field";else p:=" are fields";fi;
      pr(Join(d),p," in {b} but not in {a}");diffs:=diffs+1;fi;
    for d in Intersection(fa,fb) do 
      if d in ["operations","parent", "group","domain", "baseRing",
        "classInvariants","stabChain","groundRing", "relativeGroups",
	"transversal","Spetss", "SubCosets","ReflectionSubgroups"]
      then ;# pr("will not compare {a}.",d," and {b}.",d);
      else
        nopt:=ShallowCopy(opt);nopt.depth:=depth-1;
        nopt.na:=SPrint(opt.na,".",d); nopt.nb:=SPrint(opt.nb,".",d);
	res:=CHEVIE.Check.EqObj(a.(d),b.(d),nopt,nh());
	if res=0 then pr(nopt.na," not compared with ",nopt.nb);compared:=false;
	elif res=2 then 
	  if depth=0 then pr(nopt.na," differ with {b}");fi;diffs:=diffs+1;
	fi;
      fi;
    od;
    if diffs>0 then return 2;elif compared then return 1;else return 0;fi;
  elif IsString(a) then
    if not IsString(b) then pr("{a} is a string but {b} is not");return 2;
    elif a<>b then pr("{a}=",a," but {b}=",b);return 2;
    else return 1;
    fi;
  elif IsList(a) then
    if not IsList(b) then pr("{a} is a list but {b} is not");return 2;fi;
    if Length(a)<>Length(b) then pr("Length({a})=",Length(a),
         " but Length({b})=",Length(b));return 2;fi;
    if IsMat(a) and IsMat(b) and Length(a)=Length(a[1])
       and Length(b)=Length(b[1]) then
      p:=PermMatMat(a,b);
      if p<>false then pr("{a}=OnMatrices({b},",p,")");return 2;fi;  
    fi;
    if ForAll([1..Length(a)],i->IsBound(a[i])) and
       ForAll([1..Length(b)],i->IsBound(b[i])) then
      p:=PermListList(a,b);
      if p<>false then pr("{a}=Permuted({b},",p,")");return 2;fi;  
    fi;
    diffs:=0;compared:=0;
    for d in  [1..Length(a)] do 
      nopt:=ShallowCopy(opt);nopt.depth:=depth-1;
      nopt.na:=SPrint(opt.na,"[",d,"]"); nopt.nb:=SPrint(opt.nb,"[",d,"]");
      if IsBound(a[d]) then
	if IsBound(b[d]) then 
	  res:=CHEVIE.Check.EqObj(a[d],b[d],nopt,nh());
	  if res=0 then pr(nopt.na," not compared with ",nopt.nb);
	    compared:=compared+1;
            if compared>nbdiffs then 
	      pr("...");if diffs>0 then return 2;else return 0;fi;
	    fi;
	  elif res=2 then 
	    if depth=0 then pr(nopt.na," differ with {b}");fi;diffs:=diffs+1;
	    if diffs>nbdiffs then pr("...");return 2;fi;
	  fi;
	else pr(nopt.na," is bound but ",nopt.nb," is not");diffs:=diffs+1;
	fi;
      elif IsBound(b[d]) then
        pr(nopt.na," is not bound but ",nopt.nb," is");diffs:=diffs+1;
      fi;
    od;
    if diffs>0 then return 2;elif compared=0 then return 1;else return 0;fi;
  elif IsFunc(a) then 
    if not IsFunc(b) then pr("{a} is a function but {b} is not");return 2;
    elif a<>b then pr("Function {a} is not the same as {b}");return 2;
    else return 1;
    fi;
  elif a<>b then pr("{a}=",a," but {b}=",b);return 2;
  else return 1;
  fi;
end;

# cmptable(t1,t2[,nz]) nz: show all non-zero columns anyway
# compare two chevie tables or inductiontables
CHEVIE.Check.EqTables:=function(arg)local r,c,msg,t,p,opt,m,i;
  t:=[Copy(arg[1]),Copy(arg[2])];opt:=[];m:=[];
  for i in [1,2] do
    if IsRec(t[i])then
      if IsBound(t[i].operations) and t[i].operations=InductionTableOps then
        opt[i]:=rec(rowLabels:=t[i].gNames(t[i],rec()),
                  columnLabels:=t[i].uNames(t[i],rec()));
        m[i]:=t[i].what; t[i]:=t[i].scalar;
      else opt[i]:=rec(rowLabels:=t[i].rowLabels,
                       columnLabels:=t[i].columnLabels);
           m[i]:=""; t[i]:=t[i].scalar;
      fi;
    else opt[i]:=rec(rowLabels:=[1..Length(t[i])],
      columnLabels:=[1..Length(t[i][1])]);
      m[i]:="";
    fi;
  od;
  if Length(arg)<3 or arg[3]=1 then
    r:=Filtered([1..Length(t[1])],i->t[1][i]<>t[2][i]);
  else r:=Filtered([1..Length(t[1])],i->ForAny(t[1][i],x->x<>0*x)
    or ForAny(t[2][i],x->x<>0*x));
  fi;
  msg:=SPrint("[",Length(r),"/",Length(t[1]),",");
  if Length(r)=0 then InfoChevie("Tables agree!\n");return;fi;
  p:=PermListList(t[1]{r},t[2]{r});
  if p<>false then
    ChevieErr("Permuted lines:\n",Join(List(TransposedMat(
     [opt[1].rowLabels{r},Permuted(opt[1].rowLabels{r},p)]),
      x->Join(x,"->")),"\n"),"\n");
    return;
  fi;
  p:=PermListList(TransposedMat(t[1]),TransposedMat(t[2]));
  if p<>false then
    ChevieErr("Permuted columns:\n",
      Join(List(Cycles(p),c->Join(opt[1].columnLabels{c},"->"))," "),"\n");
    return;
  fi;
  if Length(arg)<3 or arg[3]=2 then
    c:=Filtered([1..Length(t[1][1])],i->t[1]{r}[i]<>t[2]{r}[i]);
  else c:=[1..Length(t[1][1])];
  fi;
  PrintToString(msg,Length(c),"/",Length(t[1][1]),"] of ");
  ChevieErr(msg,"\n",Join(List([1,2],i->SPrint(m[i],"\n",FormatTable(t[i],
     Inherit(opt[i],rec(rows:=r,columns:=c,screenColumns:=70)))))));
end;

CHEVIE.Check.EqCycPol:=function(arg)local a,b,na,nb,q;
  a:=arg[1];b:=arg[2];
  if a=b then return;fi;
  if IsBound(arg[3]) then na:=arg[3];else na:="a";fi;
  if IsBound(arg[4]) then nb:=arg[4];else nb:="b";fi;
  if b.coeff=0 then 
    if a.coeff<>0 then ChevieErr(na,"=",a," but ",nb,"=",b,"\n");fi;
    return;
  fi;
  q:=a/b;
  if Length(q.vcyc)=0 then ChevieErr(na,"=",
 #  FormatCoefficient(q,nb,rec(GAP:=true))," where ",nb,"=",b,"\n");
    FormatCoefficient(q,nb,rec())," where ",nb,"=",b,"\n");
  else ChevieErr(na,"=",a," but ",nb,"=",b,"\n");
  fi;
end;
