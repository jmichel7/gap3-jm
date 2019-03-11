#############################################################################
##
#A  tbl/cmplximp.g        CHEVIE library        Gunter Malle and  Jean Michel
##
#Y  Copyright (C) 1998 - 2011  The CHEVIE Team
##
##  This file contains data about imprimitive complex reflection groups
##
CHEVIE.AddData("PrintDiagram","imp",function(arg)
  local p,q,r,indices,j,indent,title,g;
  p:=arg[1];q:=arg[2];r:=arg[3];indices:=arg[4];title:=arg[5];
  Print(title," ");
  indent:=Length(title)+1;g:=i->String("",indent-i);
  if q=1 then
    Print(indices[1],"(",p,")");
    if Length(indices)>1 then Print("==");fi;
    Print(Join(indices{[2..Length(indices)]},"--"),"\n");
  elif p=q then
    Print(indices[1],"\n",g(0),"|\\\n");
    if p<>3 then Print(String(p,indent));else Print(g(0));fi;
    Print("|=",indices[3]);
    for j in [4..r] do Print(" - ",indices[j]);od;Print("\n");
    Print(g(0),"|/\n",g(0),indices[2],"\n");
  elif q=2 then
    Print(indices[2],"\n",g(2),"/3|");
    if r>=3 then Print("\\");fi;
    Print("\n");
    if p/q>2 then Print(g(Length(String(p/q))+5),"(",p/q,")");
    else Print(g(3));
    fi;
    Print(indices[1],"  | ");
    for j in [3..r] do Print(indices[j+1]);if j<>r then Print("-");fi;od;
    Print("\n",g(2),"\\ |");
    if r>=3 then Print("/");fi;
    Print("\n",String(indices[3],indent+1),
      "   ",IntListToString(indices{[1,2,3]}),
      "=",IntListToString(indices{[2,3,1]}),
      "=",IntListToString(indices{[3,1,2]}),"\n");
  else
    Print(indices[2],"\n",g(2),"/",q+1," ");
    if r>=3 then Print("\\");fi;
    Print("\n");
    if p/q>2 then Print(g(Length(SPrint(p/q))+5),"(",p/q,")");
    else Print(g(3));
    fi;
    Print(indices[1],"   ");if r>=3 then Print("=");fi;
    for j in [3..r] do Print(indices[j+1]);if j<>r then Print("-");fi;od;
    Print("\n",g(2),"\\  ");
    if r>=3 then Print("/");fi;
    Print("\n",String(indices[3],indent+1));
    j:=CHEVIE.R("BraidRelations","imp")(p,q,r);
    for g in [1..Minimum(3,r)] do
      Print("   ",IntListToString(indices{j[g][1]}),
              "=",IntListToString(indices{j[g][2]}));
    od;
    Print("\n");
  fi;
end);

CHEVIE.AddData("SemisimpleRank","imp",function(p,q,r)return r;end);

CHEVIE.AddData("BraidRelations","imp",function(p,q,r)local i,b,res;
  b:=function(i,j,o)local p;
    p:=function(i,j)return List([1..o],k->i*(k mod 2)+j*((1-k)mod 2));end;
    return [p(i,j),p(j,i)];
  end;
  res:=[];
  if q=1 then
    if r>=2 then if p=1 then Add(res,b(1,2,3)); else Add(res,b(1,2,4));fi;fi;
    Append(res,List([3..r],i->b(i,i-1,3)));
    for i in [3..r] do Append(res,List([1..i-2],j->b(i,j,2)));od;
  elif p=q then
    Add(res,b(1,2,p));
    if r>=3 then
      Append(res,[[[1,2,3,1,2,3],[3,1,2,3,1,2]],b(1,3,3),b(2,3,3)]);fi;
    Append(res,List([4..r],i->b(i,i-1,3)));
    for i in [4..r] do Append(res,List([1..i-2],j->b(i,j,2)));od;
  else Add(res,[[1,2,3],[2,3,1]]);
    i:=b(2,3,q-1);
    Add(res,[Concatenation([1,2],i[2]),Concatenation([3,1],i[1])]);
    if r>=3 then
      if q<>2 then Add(res,[[2,3,4,2,3,4],[4,2,3,4,2,3]]);fi;
      Append(res,[b(2,4,3),b(3,4,3),b(1,4,2)]);
    fi;
    Append(res,List([5..r+1],i->b(i,i-1,3)));
    for i in [5..r+1] do Append(res,List([1..i-2],j->b(i,j,2)));od;
  fi;
  return res;
end);

CHEVIE.AddData("Size","imp",function(p,q,r)return p^r*Factorial(r)/q;end);

CHEVIE.AddData("ReflectionName","imp",function(arg)local n,option;
  option:=arg[4];
  if arg[3]=1 and arg[2]=1 then 
    if IsBound(option.TeX) then return SPrint("Z_{",arg[1],"}");
    else return SPrint("Z",arg[1]);fi;
  fi;
  if IsBound(option.TeX) then n:=SPrint("G_{",Join(arg{[1..3]}),"}");
  else                     n:=SPrint("G",IntListToString(arg{[1..3]}));fi;
  if Length(arg)=5 then PrintToString(n,"(",Format(arg[4],option),")");fi;
  return n;
end);

CHEVIE.AddData("GeneratingRoots","imp",function(p,q,r)local roots,v,i;
  if q=1 then roots:=[Concatenation([1],[2..r]*0)];
  else
    if q<>p then roots:=[Concatenation([1],[2..r]*0)];
    fi;
    v:=Concatenation([-E(p),1],[3..r]*0);
    if r=2 and q>1 and q mod 2=1 then v:=v*E(p);fi; # so only 2 orbits
    if q=p then roots:=[v];else Add(roots,v);fi;
  fi;
  for i in [2..r] do v:=[1..r]*0;v[i]:=1;v[i-1]:=-1;Add(roots,v);od;
  return roots;
end);

CHEVIE.AddData("EigenvaluesGeneratingReflections","imp",
 function(p,q,r)local res;res:=[1..r]*0+1/2;
  if q=1 then res[1]:=1/p;
  elif q<>p then res:=Concatenation([q/p],res);
  fi;
  return res;
end);

CHEVIE.AddData("CartanMat", "imp",function(p,q,r)local rt,rbar,e;
  rt:=CHEVIE.R("GeneratingRoots","imp")(p,q,r);
  rbar:=ComplexConjugate(rt);
  e:=CHEVIE.R("EigenvaluesGeneratingReflections","imp")(p,q,r);
  e:=1-List(e,x->E(Denominator(x))^Numerator(x));
  e:=List([1..Length(e)],i->e[i]*rbar[i]/(rbar[i]*rt[i]));
  return List(e,x->List(rt,y->x*y));
end);

CHEVIE.AddData("ReflectionDegrees","imp",function(p,q,r)
  return Concatenation(p*[1..r-1],[r*p/q]);end);

CHEVIE.AddData("ReflectionCoDegrees","imp",function(p,q,r)local res;
  res:=p*[0..r-1];
  if p=q and p>=2 and r>2 then res[r]:=res[r]-r;fi;
  return res;
end);

CHEVIE.AddData("ParabolicRepresentatives","imp",function(p,q,r,s)local t;
  if q=1 then 
    if p=1 then 
      if s=0 then return [[]];fi;
      return List(Concatenation(List([1..r+1-s],i->Partitions(s,i))),j->
        Concatenation(List([1..Length(j)],k->Sum(j{[1..k-1]})+k-1+[1..j[k]])));
    else return Concatenation(List([0..s],i->List(
            CHEVIE.imp.ParabolicRepresentatives(1,1,r-i-1,s-i),j->
          Concatenation([1..i],i+1+j))));
    fi;
  elif r=2 then 
    if q=2 then t:=[[[]],[[1],[2],[3]],[[1..3]]];return t[s+1];
    elif p=q then 
      if p mod 2=0 then t:=[[[]],[[1],[2]],[[1,2]]];return t[s+1];
                   else t:=[[],[1],[1,2]];return t[s+1];fi;
    else return false;
    fi;
  else return false;
  fi;
end);

CHEVIE.AddData("NrConjugacyClasses","imp",function(p,q,r)
  if [q,r]=[2,2] then return p*(p+6)/4;
  elif q=1 then return NrPartitionTuples(r,p);
  else return Length(CHEVIE.R("ClassInfo","imp")(p,q,r).classtext);
  fi;
end);

CHEVIE.AddData("ClassInfo","imp",function(p,q,r)local res,times,trans,I,i,j,a,S;
  times:=function(e,o)return Concatenation(List([1..e],x->o));end;
  if [q,r]=[2,2] and not IsBound(CHEVIE.othermethod) then 
    res:=rec(classtext:=[],classparams:=[],classnames:=[]);
    for i in [0..p-1] do for j in [0..QuoInt(p-i-1,2)] do
      Add(res.classparams,Concatenation([1..j]*0+1,[1..i]*0));
      Add(res.classtext,Concatenation([1..j]*0+1,times(i,[1,2,3])));
      Add(res.classnames,String(Concatenation(times(j,"1"),times(i,"z"))));
    od;od;
    for j in [2,3] do
      for i in [0..p/2-1] do Add(res.classparams,Concatenation([j],[1..i]*0));
	Add(res.classtext,Concatenation([j],times(i,[1,2,3])));
	Add(res.classnames,String(Concatenation(String(j),times(i,"z"))));
      od;
    od;
    res.malle:=[];
    for a in [0..p-1] do 
       Append(res.malle,List([0..QuoInt(p-a-1,2)],m->[3,a,m]));od;
    Append(res.malle,List([0..p/2-1],m->[1,m]));
    Append(res.malle,List([0..p/2-1],m->[2,m]));
    res.orders:=List(res.classparams,function(c)
     if Length(c)>0 and c[1] in [2,3] then return Lcm(2,p/Gcd(Number(c,x->x=0),p));
     else return Lcm(p/Gcd(Number(c,x->x=0),p),(p/2)/Gcd(Number(c,x->x=1),p/2));
     fi;end);
    res.classes:=List(res.classparams,function(c)
      if Length(c)>0 and c[1] in [2,3] then return p/q;
      elif 1 in c then return 2;else return 1;fi;end);
    return res;
  elif q=1 then
    res:=rec(classparams:=PartitionTuples(r,p));
    res.classtext:=List(res.classparams,function(S)local l,w,d;
      S:=Concatenation(List([1..p], i->List(S[i],t->[t,i-1])));
      SortBy(S,a->[a[1],-a[2]]);
      l:=0;w:=[];
      for d in S do
        Append(w,times(d[2],Concatenation([l+1,l..2],[1..l+1])));
	# non-reduced word because this is the one used by Halverson-Ram
	# for characters of the Hecke algebra (see below).
	Append(w,[l+2..l+d[1]]);
	l:=l+d[1];
      od;
      return w;
    end);
    res.classnames:=List(res.classparams,CHEVIE.R("ClassName","imp"));
    res.orders:=List(res.classparams,m->Lcm(List([1..Length(m)],function(i)
       if Length(m[i])=0 then return 1;
       else return Lcm(m[i]*p/Gcd(i-1,p));fi;end)));
    res.centralizers:=List(res.classparams,m->p^Sum(m,Length)*
     Product(List(m,pp->Product(Collected(pp),y->Factorial(y[2])*y[1]^y[2]))));
    res.classes:=List(res.centralizers,x->p^r*Factorial(r)/x);
    return res;
  else
  # According  to Hugues  ``On  decompositions  in complex  imprimitive 
  # reflection groups'' Indagationes 88 (1985) 207--219:                
  #
  # Let l=(S_0,..,S_{p-1}) be  a p-partition of r specifying  a class C 
  # of G(p,1,r) as  in the above code;  C is in G(p,q,r)  iff q divides 
  # sum_i i*|S_i|;  C splits  in d  classes for  the largest  d|q which 
  # divides all parts  of all S_i and  such that |S_i|=0 if  d does not 
  # divide i;  if w is in  C and t  is the first generator  of G(p,1,r) 
  # then t^i w t^-i for i in [0..d-1] are representatives of classes of 
  # G(p,q,r) which meet C.                                              

    trans:=function(w)local d,res,l,i,add,word;
    # translate words  in G(p,1,r) into  words of G(p,q,r); use  that if
    # t,s2 (resp. s1,s2)  are the first 2 generators  of G(p,1,r) (resp.
    # G(p,p,r)) then  s1=s2^t thus  s2^(t^i)= (s1s2)^i  s2 [the  first 3
    # generators of G(p,q,r) are t^q,s1,s2].
    #Print(IntListToString(w),"=>");
     d:=0;res:=[];
     word:=function(l,i)return List(i+[l,l-1..1],j->1+(j mod 2));end;
     add:=function(a)local l; # here we try to reduce words
       l:=Length(res);
       if l>0 and res[l]=a then res:=res{[1..l-1]};
       elif p=q and a in [1,2] and l>=q and res{[l-q+1..l]}=word(q,3-a)
       then res:=Concatenation(res{[1..l-q]},word(q-1,3-a));
       else Add(res,a);
       fi;
     end;
     for l in w do
       if l=1 then d:=d+1;
       elif l<>2 then add(l);
       else d:=d mod p;
	 if d=0 then add(2);
	 else for i in [1..p-d-1] do add(1);add(2);od;add(1);
	 fi;
       fi;
     od;
     d:=d mod p;
     if d mod q<>0 then Error();
     elif d<>0 then res:=Concatenation(1+res,[1..d/q]*0+1);
     elif p<>q then res:=1+res;
     fi;
    #Print(IntListToString(res),"\n");
     return res;
    end;

    I:=CHEVIE.R("ClassInfo","imp")(p,1,r);
    res:=rec(classtext:=[],classparams:=[],classnames:=[],orders:=[],
      centralizers:=[]);
    for i in Filtered([1..Length(I.classparams)],i->
      List(I.classparams[i],Length)*[0..p-1] mod q=0)
    do
      S:=I.classparams[i];
      a:=Concatenation(S);Add(a,q);
      Append(a,Filtered([1..p],j->Length(S[j])<>0)-1);
      a:=ApplyFunc(Gcd,a); # number of pieces the class splits
      for j in [0..a-1] do
    #   Print("i=",i," j=",j," text=<",IntListToString(I.classtext[i]),"> ");
	Add(res.classtext,
	  trans(Concatenation([1..j]*0+1,I.classtext[i],[1..p-j]*0+1)));
	if a>1 then Add(res.classparams,Concatenation(S,[p*j/a]));
	else Add(res.classparams,S);
	fi;
	Add(res.orders,I.orders[i]);
	Add(res.centralizers,I.centralizers[i]*a/q);
      od;
    od;
    res.classes:=List(res.centralizers,x->res.centralizers[1]/x);
    res.classnames:=List(res.classparams,CHEVIE.R("ClassName","imp"));
    return res;
  fi;
end);

CHEVIE.AddData("ClassName", "imp", function(p)local j,p1;
  if IsList(p) and ForAll(p, IsList) then 
    if Sum(p,Sum)=1 then return Format(E(Length(p))^(Position(p,[1])-1));
    else return PartitionTupleToString(p);
    fi;
  elif IsList(p) and ForAll(p, IsInt) then return IntListToString(p);
  elif IsList(p) and ForAll(p{[1..Length(p)-1]},IsList) and IsInt(p[Length(p)])
  then p1:=p{[1..Length(p)-1]};Append(p1,[Length(p1),p[Length(p)]]);
    return PartitionTupleToString(p1);
  else Error(); # should not happen
  fi;
end);

CHEVIE.AddData("PowerMaps","imp",function(p,q,r)local pow,pp,pw,res;
  if q=1 then
    pow:=function(p,n)local e,res,k,l,g,j;
      e:=Length(p);
      res:=List([1..e],x->[]);
      for k in [1..e] do
       for l in p[k] do
	 g:=Gcd(n,l);
	 for j in [1..g] do Add(res[1+(QuoInt(n*(k-1),g)mod e)],l/g);od;
       od;
      od;
      for k in [1..e] do
	Sort(res[k]);res[k]:=Reversed(res[k]);
      od;
      return res;
    end;

    pp:=CHEVIE.R("ClassInfo","imp")(p,q,r).classparams;
    res:=[];
    for pw in Set(Factors(Factorial(r)*p)) do
       res[pw]:=List(pp,x->Position(pp,pow(x,pw)));
    od;
    return res;
  else
    InfoChevie("# PowerMaps not implemented for G(",p,",",q,",",r,")\n");
    return false;
  fi;
end);

CHEVIE.AddData("CharInfo","imp",function(de,e,r)local d,ct,res,t,tt,s,fd;
  res:=rec();d:=QuoInt(de,e);
  if e=1 then res.charparams:=PartitionTuples(r,de); s:=[1..d]*0;s[1]:=1;
    res.charSymbols:=List(res.charparams,x->SymbolPartitionTuple(x,s));
  else
    res.charparams:=[];
    for t in PartitionTuples(r,de) do
      tt:=List([1..e]*d,i->Rotation(t,i));
      if t=Minimum(tt) then
	s:=Position(tt,t);
	if s=e then Add(res.charparams,t);
	else t:=t{[1..s*d]}; s:=e/s;
	  Append(res.charparams,List([0..s-1],i->Concatenation(t,[s,i])));
	fi;
      fi;
    od;
    if d=1 then
      res.charSymbols:=List(res.charparams,x->SymbolPartitionTuple(x,0));
    fi;
    if d>1 and e mod 2=0 and r=2 then
# .malle: indexing of chars as in Malle's paper on rank 2 cyclotomic algebras.
      res.malle:=List(res.charparams,function(t)local pos,de;
	  if IsInt(t[Length(t)]) then
	    if t[Length(t)]=0 then return [1,2,1,Position(t,[1])];
	    else                   return [1,1,2,Position(t,[1])];
	    fi;
	  else de:=Length(t)/2;
	    pos:=Filtered([1..Length(t)],i->Length(t[i])>0);
	    if Length(pos)=1 then
	      if t[pos[1]]=[2] then return [1,1,1,pos[1]-de];
	      else                  return [1,2,2,pos[1]-de];
	      fi;
	    elif pos[1]<=de then return [2,-1,pos[1],pos[2]-de];
	    else                return [2,1,pos[2]-de,pos[1]-de];
	    fi;
	  fi;
	end);
    elif [de,e,r]=[3,3,3] then res.malle:=
     [[2,3,2],[2,3,3],[2,3,1],[3,4],[3,5],[1,9],[3,2],[3,1],[2,3,4],[1,0]];
    elif [de,e,r]=[3,3,4] then res.malle:=
     [[12,6],[4,10],[6,8],[4,11],[1,18],[12,3],[6,5,2],[8,4],[8,5],
      [6,5,1],[3,9],[6,2],[2,6],[4,2],[4,1],[3,3],[1,0]];
# here the labeling is defined by phi_{6,5}' being the one which appears
# in the tensor square of the reflection representation phi_{4,1}
    elif [de,e,r]=[3,3,5] then res.malle:=
     [[30,10],[20,12],[5,19],[10,14],[10,15],[5,20],[1,30],[30,7,1],[40,6],
      [30,7,2],[10,11],[15,10],[20,9],[20,8],[15,11],[10,12],[4,18],
      [30,4],[20,5],[10,8],[10,7],[20,6],[5,12],[20,3],[10,6],[15,4],[15,5],
      [10,5],[6,9],[10,3],[10,2],[5,6],[5,2],[5,1],[4,3],[1,0]];
  # here the labeling is defined by phi_{30,7}'' being the one which appears
  # in the tensor 4th power of the reflection representation phi_{5,1}
    elif [de,e,r]=[4,4,3] then res.malle:=[[6,3],[3,6,1],[3,5],[3,6,2],
      [1,12],[3,2,1],[3,2,2],[3,1],[2,4],[1,0]];
  # here the labeling is defined by phi_{3,2}'' being the complex
  # conjugate of phi_{3,1} and phi_{3,6}'' the complex conjugate of phi_{3,5}
    fi;
  fi;
  t:=List([r,r-1..0],function(i)local v;
    v:=List([1..de],x->[]);if i>0 then v[1]:=[i];fi;
    v[2]:=[1..r-i]*0+1;return v;end);
  if e>1 then t:=List(t,v->Minimum(List([1..e],i->Rotation(v,i*d))));fi;
  res.extRefl:=List(t,v->Position(res.charparams,v));
  if e=1 or d=1 then
    res.A:=List(res.charSymbols,HighestPowerGenericDegreeSymbol);
    res.a:=List(res.charSymbols,LowestPowerGenericDegreeSymbol);
    res.B:=List(res.charSymbols,HighestPowerFakeDegreeSymbol);
    res.b:=List(res.charSymbols,LowestPowerFakeDegreeSymbol);
  fi;
  if e>1 and d>1 then
    res.opdam:=PermListList(res.charparams,List(res.charparams,
      function(s)
        if not IsList(s[Length(s)]) then
          s:=Copy(s);t:=(Length(s)-2)/d;
          s{[0..t-1]*d+1}:=Rotation(s{[0..t-1]*d+1},1);
          s{[1..Length(s)-2]}:=Minimum(List([1..t],
              i->Rotation(s{[1..Length(s)-2]},i*d)));
          return s;
        fi;
        s:=ShallowCopy(s);s{[0..e-1]*d+1}:=Rotation(s{[0..e-1]*d+1},1);
        return Minimum(List([1..e],i->Rotation(s,i*d)));
        end));
  fi;
  return res;
end);

CHEVIE.AddData("LowestPowerFakeDegrees","imp",function(p,q,r)local ci;
  if q=1 or p=q then Error("should not be called");fi;
  return false;end);

CHEVIE.AddData("HighestPowerFakeDegrees","imp",function(p,q,r)local ci;
  if q=1 or p=q then Error("should not be called");fi;
  return false;end);

CHEVIE.AddData("CharSymbols","imp",function(p,q,r)local s,ss,res;
  if q=1 then return SymbolsDefect(p,r,0,1);
  elif q=p then ss:=SymbolsDefect(p,r,0,0);
    res:=[];
    for s in ss do p:=Position(Rotations(s){[2..Length(s)]},s);
      if p=false then Add(res,s);
      else Append(res,List([0..Length(s)/p-1],
        i->Concatenation(List(s{[1..p]},ShallowCopy),[Length(s)/p,i])));
      fi;
    od;
    return res;
  else return false;
  fi;
end);

CHEVIE.AddData("FakeDegree","imp",function(p,q,r,c,v)
  if q=1 then c:=CycPolFakeDegreeSymbol(SymbolPartitionTuple(c,1));
  elif q=p then c:=CycPolFakeDegreeSymbol(SymbolPartitionTuple(c,[1..p]*0));
  else return false;
  fi;
  return Value(c,v);
end);

CHEVIE.AddData("CharName","imp",function(p,q,r,s,option)
  if RankSymbol(s)=1 then 
       return Format(E(Length(s))^(Position(s,[1])-1),option);
  else return PartitionTupleToString(s,option);
  fi;
end);

CHEVIE.AddData("SchurModel","imp",function(p,q,r,phi)
  local l,i,j,res,s,t,ci,GenHooks,v,h,d;
  if q=1 then # cf. Chlouveraki, arxiv 1101.1465
    GenHooks:=function(l,m)if Length(l)=0 then return [];fi;
      m:=AssociatedPartition(m);Append(m,[1..l[1]-Length(m)]*0);
      m:=1+m-[1..Length(m)];
      return Concatenation(List([1..Length(l)],i->l[i]-i+m{[1..l[i]]}));
    end;
    res:=rec(coeff:=(-1)^(r*(p-1)),factor:=[1..p]*0,vcyc:=[]);
    l:=Concatenation(phi);Sort(l);Add(res.factor,([1..Length(l)]-Length(l))*l);
    for s in [1..p] do for t in [1..p] do for h in GenHooks(phi[s],phi[t]) do
      v:=[1..p]*0;
      if s<>t then v{[s,t]}:=[1,-1];Add(v,h);Add(res.vcyc,[v,1]);
      else Add(v,1);
        for d in DivisorsInt(h) do if d>1 then Add(res.vcyc,[v,d]);fi;od;
      fi;
    od;od;od;
    return res;
  elif [q,r]=[2,2] then
    ci:=CHEVIE.imp.CharInfo(p,q,r);
    phi:=ci.malle[Position(ci.charparams,phi)];
    if phi[1]=1 then
      res:=rec(coeff:=1,factor:=[1..4+p/2]*0, vcyc:=[]);
      for l in [[1,-1,0,0],[0,0,1,-1]] do
        Append(l,[1..p/2]*0);Add(res.vcyc,[l,1]);
      od;
      for i in [2..p/2] do for l in [[0,0,0,0,1],[1,-1,1,-1,1]] do
        Append(l,[1..p/2-1]*0);l[4+i]:=-1;Add(res.vcyc,[l,1]);
      od;od;
    else
      res:=rec(coeff:=-2,factor:=[1..4+p/2]*0, vcyc:=[],root:=[1..4+p/2]*0);
      res.rootCoeff:=E(p/2)^(2-phi[3]-phi[4]);
      res.root{[1..6]}:=[1,1,1,1,1,1]/2;
      for i in [3..p/2] do for j in [1,2] do
	l:=[1..4+p/2]*0;l{4+[j,i]}:=[1,-1];Add(res.vcyc,[l,1]);
      od;od;
      if IsBound(CHEVIE.old) then
      for l in [[0,-1,0,-1,-1,0],[0,-1,-1,0,-1,0],
                [-1,0,-1,0,-1,0],[-1,0,0,-1,-1,0]] do
        Append(l,[1..p/2-2]*0);Add(l,1);Add(res.vcyc,[l,1]);
      od;
      else
      for l in [[0,-1,0,-1,-1,0],[0,-1,-1,0,0,-1],
                [-1,0,-1,0,-1,0],[-1,0,0,-1,0,-1]] do
        Append(l,[1..p/2-2]*0);Add(l,1);Add(res.vcyc,[l,1]);
      od;
      fi;
    fi;
    return res;
  else Error("not implemented");
  fi;
end);

CHEVIE.AddData("SchurData","imp",function(p,q,r,phi)local ci,res;
  if [q,r]=[2,2] then 
    ci:=CHEVIE.imp.CharInfo(p,q,r);
    phi:=ci.malle[Position(ci.charparams,phi)];
    if phi[1]=1 then
      res:=rec(order:=[phi[2],3-phi[2],2+phi[3],5-phi[3],4+phi[4]]);
      Append(res.order,4+Difference([1..p/2],[phi[4]]));
      return res;
    else
      res:=rec(order:=[1,2,3,4,4+phi[3],4+phi[4]]);
      Append(res.order,4+Difference([1..p/2],phi{[3,4]}));
      res.rootPower:=phi[2]*E(p)^(phi[3]+phi[4]-2);
      return res;
    fi;
  else Error("not implemented");
  fi;
end);

CHEVIE.AddData("SchurElement","imp",function(p,q,r,phi,para,root)local m;
  if r=1 then return VcycSchurElement(Concatenation(para[1],[0]),
    CHEVIE.imp.SchurModel(p,q,r,phi));
  elif p=1 then return VcycSchurElement([0,-para[1][1]/para[1][2]],
    CHEVIE.imp.SchurModel(p,q,r,phi));
  elif q=1 then return VcycSchurElement(Concatenation(para[1],
    [-para[2][1]/para[2][2]]),CHEVIE.imp.SchurModel(p,q,r,phi));
  elif [q,r]=[2,2] then return VcycSchurElement(Concatenation(para{[2,3,1]}),
      CHEVIE.imp.SchurModel(p,q,r,phi),CHEVIE.imp.SchurData(p,q,r,phi));
  elif p=q then
    if IsInt(phi[Length(phi)]) then m:=Length(phi)-2;phi:=FullSymbol(phi);
    else m:=p;
    fi;
    return CHEVIE.imp.SchurElement(p,1,r,phi,
      Concatenation([List([0..p-1],i->E(p)^i)],para{[2..Length(para)]}),[])/m;
  elif para[2]=para[3] then
    if IsInt(phi[Length(phi)]) then m:=Length(phi)-2;phi:=FullSymbol(phi);
    else m:=p;
    fi;
    if para[1]=List([1..p/q],i->E(p/q)^(i-1)) then
        para:=[List([0..p-1],i->E(p)^i),para[2]];
    else para:=[Concatenation(TransposedMat(List(para[1],i->List([0..q-1],
      j->E(q)^j)*GetRoot(i,q)))),para[2]];
    fi;
    return p/q*CHEVIE.imp.SchurElement(p,1,r,phi,para,[])/m;
  else CHEVIE.compat.InfoChevie("# SchurElements(H(G(",
                     p,",",q,",",r,"),",para,") not implemented\n");
    return false;
  fi;
end);

CHEVIE.AddData("FactorizedSchurElement","imp",function(p,q,r,phi,para,root)
  local m,F;
  if r=1 then return VFactorSchurElement(Concatenation(para[1],[0]),
    CHEVIE.imp.SchurModel(p,q,r,phi));
  elif p=1 then return VFactorSchurElement([0,-para[1][1]/para[1][2]],
    CHEVIE.imp.SchurModel(p,q,r,phi));
  elif q=1 then return VFactorSchurElement(Concatenation(para[1],
    [-para[2][1]/para[2][2]]),CHEVIE.imp.SchurModel(p,q,r,phi));
  elif [q,r]=[2,2] then return VFactorSchurElement(Concatenation(para{[2,3,1]}),
      CHEVIE.imp.SchurModel(p,q,r,phi),CHEVIE.imp.SchurData(p,q,r,phi));
  elif p=q then
    if IsInt(phi[Length(phi)]) then m:=Length(phi)-2;phi:=FullSymbol(phi);
    else m:=p;
    fi;
    F:=CHEVIE.imp.FactorizedSchurElement(p,1,r,phi,
      Concatenation([List([0..p-1],i->E(p)^i)],para{[2..Length(para)]}),[]);
    F.factor:=F.factor/m;
    return F;
  elif para[2]=para[3] then
    if IsInt(phi[Length(phi)]) then m:=Length(phi)-2;phi:=FullSymbol(phi);
    else m:=p;
    fi;
    if para[1]=List([1..p/q],i->E(p/q)^(i-1)) then
        para:=[List([0..p-1],i->E(p)^i),para[2]];
    else para:=[Concatenation(TransposedMat(List(para[1],i->List([0..q-1],
      j->E(q)^j)*GetRoot(i,q)))),para[2]];
    fi;
    F:=CHEVIE.imp.FactorizedSchurElement(p,1,r,phi,para,[]);
    F.factor:=p/(q*m)*F.factor;
    return F;
    #return p/q*CHEVIE.imp.FactorizedSchurElement(p,1,r,phi,para)/m;
   else CHEVIE.compat.InfoChevie("# FactorizedSchurElements(H(G(",
                     p,",",q,",",r,"),",para,") not implemented\n");
    return false;
  fi;
end);

CHEVIE.AddData("HeckeCharTable","imp",function(p,q,r,para,root)
  local X,Y,Z,res,cl,GenericEntry,pow,d,I,LIM,HooksBeta,StripsBeta,Strips,
                                Delta,StripsCache,chiCache,code,j,ci;
  res:=rec();
  res.name:=SPrint("H(G(",p,",",q,",",r,"))"); res.identifier:=res.name;
  res.degrees:=CHEVIE.R("ReflectionDegrees","imp")(p,q,r);
  res.size:=Product(res.degrees); res.order:=res.size;
  res.dim:=r;
  cl:=CHEVIE.R("ClassInfo","imp")(p,q,r);
  if r=1 then
    Inherit(res,cl,["classes","orders"]);
    res.irreducibles:=List([1..p],i->List([0..p-1],j->para[1][i]^j));
    res.powermap:=CHEVIE.R("PowerMaps","imp")(p,q,r);
  elif q=1 then
    Inherit(res,cl);
    res.powermap:=CHEVIE.R("PowerMaps","imp")(p,q,r);

    # character table of the Hecke algebra of G(p,1,r) parameters v and [q,-1]
    # according to "Characters of Iwahori-Hecke algebras of G(r,p,n)"
    # Halverson and Ram, Canadian Journal of math. 50 (1998) 167--192

    # The next function returns the list of all hooks of area less than or equal
    # to s in the Ferrers-Young diagram of the partition with beta-numbers S.
    # Each hook is returned as a record with the following fields:
    #  .area
    #  .hooklength number of rows -1
    #  .start      the beta-number which is decreased by removing the hook
    #  .startpos   the position of that beta-number in the list
    #  .stoppos    the position it will occupy after being decreased
    #  .DC         the list of dull corners, given by their axial position
    #  .SC         the list of sharp corners, given by their axial position
    HooksBeta:=function(S,s)local res,i,j,e,k,z,zi;
      res:=[];e:=Length(S);
      if e=0 then return res;fi;
      j:=e;
      for i in [S[e]-1,S[e]-2..0] do
	if not(i in S) then
	while j>0 and S[j]>i do j:=j-1;od;
	k:=j+1;
	while k<=e and S[k]-i<=s do
	  z:=[i];Append(z,S{[j+1..k-1]});
	  zi:=Filtered([2..Length(z)],i->z[i]-z[i-1]>1);
	  Add(res,rec(area:=S[k]-i,hooklength:=k-j-1,start:=S[k],startpos:=k,
	    stoppos:=j+1,DC:=z{zi}-e,
	    SC:=z{Concatenation(zi-1,[Length(z)])}+1-e));
	  k:=k+1;
	od;
	fi;
      od;
      return res;
    end;

    # The next function returns as a list of lists all broken border strips
    # (union of disjoint hooks) of area less than or equal to s in the diagram
    # of the partition  with beta-numbers S.
    StripsBeta:=function(S,s)local res,j,hook,hs,h;
      res:=[[]];
      for hook in HooksBeta(S,s) do
	if s=hook.area then Add(res,[hook]);
	else j:=hook.stoppos-1-Length(S);
	  for hs in StripsBeta(S{[1..hook.stoppos-1]},s-hook.area) do
	    for h in hs do h.SC:=h.SC+j;h.DC:=h.DC+j;od;
	    Add(hs,hook);Add(res,hs);
	  od;
	fi;
      od;
      return res;
    end;

    StripsCache:=rec();

    code:=function(arg)local S,res,p;res:=[];
      for S in arg do for p in S do Append(res,p);Add(res,-1);od;od;
      p:=".0123456789abcdefghijklmnopqrstuvwxyz";return p{2+res};
    end;

    # The next function returns as a GAP list all collections of broken border
    # strips of total area equal to s coming from the various beta-lists in the
    # symbol S. Each collection is represented by a GAP record containing the
    # statistical information about it necessary to compute the function Delta
    # in Ram-Halverson 2.17. These records have the following fields:
    #   area       total area
    #   cc         number of connected components (hooks)
    #   hooklength Sum of hooklengths of individual hooks
    #   DC         the list of all dull corners, represented by a pair:
    #                which beta list they come from, axial position
    #   SC         the list of all sharp corners, represented by a pair:
    #                which beta list they come from, axial position
    #   remainder  the symbol left after removing the strip collection
    Strips:=function(S,s)
      local apply, e, name, res, hs, ss, a, r;

      apply:=function(S,hs)local h;
	S:=ShallowCopy(S);
	for h in hs do S{[h.stoppos..h.startpos]}:=
	  Concatenation([h.start-h.area],S{[h.stoppos..h.startpos-1]});
        od;
	while Length(S)>0 and S[1]=0 do S:=S{[2..Length(S)]}-1;od;
	return S;
      end;
      e:=Length(S);
      if e=0 then if s=0 then
       return [rec(SC:=[],DC:=[],cc:=0,hooklength:=0,area:=0,remainder:=[])];
      else return [];fi;fi;
      name:=code(S,[[s]]);
      if IsBound(StripsCache.(name)) then return StripsCache.(name);fi;
      res:=[];
      for hs in StripsBeta(S[e],s) do
	hs:=rec(area:=Sum(hs,x->x.area),
		cc:=Length(hs),
		hooklength:=Sum(hs,x->x.hooklength),
		SC:=Concatenation(List(hs,x->List(x.SC,y->[e,y]))),
		DC:=Concatenation(List(hs,x->List(x.DC,y->[e,y]))),
		remainder:=apply(S[e],hs)
		);
	for a in Strips(S{[1..e-1]},s-hs.area) do
          # since there is no `Copy' in GAP 4
          ss := rec();
          for r in RecFields(a) do
            ss.(r) := ShallowCopy(a.(r));
          od;
	  Append(ss.SC,hs.SC);
	  Append(ss.DC,hs.DC);
	  Add(ss.remainder,hs.remainder);
	  ss.cc:=ss.cc+hs.cc;
	  ss.hooklength:=ss.hooklength+hs.hooklength;
	  ss.area:=ss.area+hs.area;
	  Add(res,ss);
	od;
      od;
      StripsCache.(name):=res;
      return res;
    end;

    # the function Delta of Ram-Halverson 2.17, modified to take in account that
    # our eigenvalues for T_2..T_r are Q[1] and Q[2] instead of q and -q^-1
    Delta:=function(k,hs,Q,v)local res,ctSC,ctDC,q,
      ElementarySymmetricFunction,HomogeneousSymmetricFunction;
      res:=1;
      if hs.cc>1 then
	  if k=1 or Q[1]=-Q[2] then return 0;
	  else res:=res*(Q[1]+Q[2])^(hs.cc-1);
	  fi;
      fi;
      q:=-Q[1]/Q[2];
      res:=res*(-1)^hs.hooklength*Q[1]^(hs.area-hs.cc)*q^(-hs.hooklength);
      if k=0 then return res;fi;
      ctSC:=List(hs.SC,x->v[x[1]]*q^x[2]);
      ctDC:=List(hs.DC,x->v[x[1]]*q^x[2]);
      res:=res*Product(ctSC)*Product(ctDC)^-1;
      if k=1 then return res;fi;

      ElementarySymmetricFunction:=function(t,v)
	return Sum(Combinations([1..Length(v)],t),x->Product(v{x}));end;

      HomogeneousSymmetricFunction:=function(t,v)
	return Sum(Combinations(Concatenation(List([1..t],
	          x->[1..Length(v)])),t),x->Product(v{x}));end;

      return res*(-1)^(hs.cc-1)*Sum(List([0..Minimum(Length(ctDC),k-hs.cc)],
	t->(-1)^t*ElementarySymmetricFunction(t,ctDC)*
	HomogeneousSymmetricFunction(k-t-hs.cc,ctSC)));
    end;

    chiCache:=rec(); LIM:=r;

    GenericEntry:=function(lambda,mu)local bp,i,rest,res,name,n;
      n:=Sum(lambda,Sum);
      if n=0 then return 1;fi;
      if n<LIM then
	name:=code(lambda,mu);
	if IsBound(chiCache.(name)) then return chiCache.(name);fi;
      fi;
      bp:=Maximum(Concatenation(lambda));i:=PositionProperty(lambda,x->bp in x);
    # here choice of bp and i corresponds to choice (Sort) in classtext
      rest:=ShallowCopy(lambda);rest[i]:=rest[i]{[2..Length(rest[i])]};
      res:=(-Product(para[2]))^((i-1)*(n-bp))*Sum(Strips(mu,bp),
	function(x)local d;d:=Delta(i-1,x,para[2],para[1]);
	 if d=0 then return d;
	 else return d*GenericEntry(rest,x.remainder);
	 fi;
	 end);
      if n<LIM then chiCache.(name):=res;fi;
      return res;
    end;

    res.irreducibles:=List(List(cl.classparams,x->List(x,BetaSet)),
        x->List(cl.classparams,y->GenericEntry(y,x)));
  elif [q,r]=[2,2] and not IsBound(CHEVIE.othermethod) then
    Inherit(res,cl,["classes","orders"]);
    X:=para[2];Y:=para[3];Z:=para[1];
    ci:=CHEVIE.R("CharInfo","imp")(p,q,r);
    GenericEntry:=function(char,class)local w;
      char:=ci.malle[Position(ci.charparams,char)];
      if char[1]=1 then
        w:=[Z[char[4]],X[char[2]],Y[char[3]]];
        return Product(class,function(i)
         if i=0 then return Product(w); else return w[i]; fi;end);
      else
        w:=char[2]*GetRoot(X[1]*X[2]*Y[1]*Y[2]*Z[char[3]]*Z[char[4]]*
          E(p/q)^(2-char[3]-char[4]),2)*E(p)^(char[3]+char[4]-2);
        class:=List([0..3],i->Number(class,j->i=j));
        if class[2]>0 then char:=Sum(Z{char{[3,4]}},x->x^class[2]);
        elif class[3]>0 then char:=Sum(X);
        elif class[4]>0 then char:=Sum(Y);
        else char:=2;
        fi;
        return w^class[1]*char;
      fi;
    end;
    res.irreducibles:=List(ci.charparams,
      char->List(cl.classparams,class->GenericEntry(char,class)));
  else 
    Inherit(res,cl,["centralizers","orders","classnames"]);
    res.classes:=List(res.centralizers,x->res.size/x);
    res.irreducibles:=List([1..Length(res.classes)],i->CharRepresentationWords(
      CHEVIE.R("HeckeRepresentation","imp")(p,q,r,para,[],i),cl.classtext));
  fi;
  res.centralizers:=List(res.classes,x->res.size/x);
  res.parameter:=para;
  res.irreducibles:=res.irreducibles*Product(para,Product)^0;
  return CHEVIE.compat.MakeCharacterTable(res);
end);

CHEVIE.AddData("HeckeRepresentation","imp",function(p,q,r,para,root,i)
  local X,Y,t,x,a,v,d,T,S,m,extra,l,m1,p1rRep,f;
  if not IsList(para) then para:=[para];fi;
  if [q,r]=[1,2] then X:=para[2];Y:=para[1];#integral matrices in this case
    t:=PartitionTuples(2,p)[i];
    if Number(t,x->x<>[])=1 then
      p:=PositionProperty(t,x->x<>[]);
      if t[p]=[2] then return X[1]^0*[[[Y[p]]],[[X[1]]]];
      else return X[1]^0*[[[Y[p]]],[[X[2]]]];
      fi;
    else  p:=Filtered([1..Length(t)],i->t[i]<>[]);
      return X[1]^0*[[[Y[p[1]],0],[-1,Y[p[2]]]],
		     [[X[1],X[1]*Y[p[1]]+X[2]*Y[p[2]]],[0,X[2]]]];
    fi;
  elif [p,q,r]=[3,3,3] then 
    x:=-para[2][1]/para[2][2];
    f:=function(x,j)return [[[-1,0,0],[0,0,1],[0,x,-1+x]],[[-1,0,0],
     [x-x^2,-1+x,j^2],[j*x-j*x^2,j*x,0]],[[0,1,0],[x,-1+x,0],[0,0,-1]]];
    end;
    r:=x^0*[[[[-1,0],[-1,x]],[[x,-x],[0,-1]],[[x,-x],[0,-1]]],
      [[[-1,0],[-1,x]],[[x,-x],[0,-1]],[[-1,0],[-1,x]]],
      [[[-1,0],[-1,x]],[[x,-x],[0,-1]],[[-1+x,1],[x,0]]],
      f(x,E(3)),f(x,E(3)^2),
      [[[-1]],[[-1]],[[-1]]],
      -x*f(x^-1,E(3)^2),-x*f(x^-1,E(3)),
      [[[-1,0],[-1,x]],[[-1,0],[-1,x]],[[x,-x],[0,-1]]],
      [[[x]],[[x]],[[x]]]];
    return r[i];
  elif [p,q,r]=[2,2,4] then 
    x:=-para[1][1]/para[1][2];r:=[
 x->[[[-1+x,-1,0],[-x,0,0],[x-x^2,-1+x,-1]],[[0,1,0],[x,-1+x,0],[0,0,-1]],
     [[-1,0,0],[0,0,1],[0,x,-1+x]],[[0,1,0],[x,-1+x,0],[0,0,-1]]],
 x->[[[0,1,0],[x,-1+x,0],[0,0,-1]],[[-1+x,-1,0],[-x,0,0],[x-x^2,-1+x,-1]],
     [[-1,0,0],[0,0,1],[0,x,-1+x]],[[0,1,0],[x,-1+x,0],[0,0,-1]]],
 x->[[[-1,0,0,0],[0,-1+x,-1,0],[0,-x,0,0],[0,0,0,-1]],[[-1,1-x,1-x,0],[0,0,1,0],
      [0,x,-1+x,0],[0,-1+x,-1+x,-1]],[[-1+x,-x,0,0],[-1,0,0,0],[0,0,-1,0],
      [0,0,0,-1]],[[0,0,0,1],[0,-1,0,0],[0,0,-1,0],[x,0,0,-1+x]]],
 x->[[[-1]],[[-1]],[[-1]],[[-1]]],
 x->[[[x,1-x,-1+x,-x+x^2,x-x^2,0],[0,-1+x,0,0,-x,x-x^2],[0,0,-1+x,-x,0,x-x^2],
      [0,0,-1,0,0,-1+x],[0,-1,0,0,0,-1+x],[0,0,0,0,0,-1]],[[x,0,0,0,0,0],
      [0,0,0,0,x,0],[0,0,0,x,0,0],[0,0,1,-1+x,0,0],[0,1,0,0,-1+x,0],
      [0,0,0,0,0,-1]],[[0,0,x,0,0,0],[0,-1,0,0,0,0],[1,0,-1+x,0,0,0],
      [0,0,0,x,0,0],[0,0,0,0,0,x],[0,0,0,0,1,-1+x]],[[-1,0,0,0,0,0],
      [0,-1+x,1,0,0,0],[0,x,0,0,0,0],[0,0,0,0,x,0],[0,0,0,1,-1+x,0],
      [0,0,0,0,0,x]]],
 x->[[[-1+x,0,-1,0,0,0,0,0],[0,0,0,0,1,0,0,0],[-x,0,0,0,0,0,0,0],
      [0,0,0,0,0,0,1,0],[0,x,0,0,-1+x,0,0,0],[0,0,0,0,0,-1+x,0,x],
      [0,0,0,x,0,0,-1+x,0],[0,0,0,0,0,1,0,0]],[[0,0,1,0,0,0,0,0],
      [0,0,-1+x,0,1,0,(-1+x)/x,0],[x,0,-1+x,0,0,0,0,0],[0,0,0,-1+x,0,0,-1,0],
      [x-x^2,x,0,-1+x,-1+x,0,(1-2*x+x^2)/x,0],[-x+x^2,0,0,-x+x^2,0,-1+x,1-x,x],
      [0,0,0,-x,0,0,0,0],[0,0,1-x,x-x^2,0,1,-1+x,0]],[[0,1,0,0,0,0,0,0],
      [x,-1+x,0,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,x,-1+x,0,0,0,0],
      [0,0,0,0,-1+x,0,-1,0],[0,0,0,0,0,x,0,0],[0,0,0,0,-x,0,0,0],
      [0,0,0,0,0,-x,0,-1]],[[-1,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],
      [0,0,-1,-x,0,0,0,0],[0,0,0,x,0,0,0,0],[0,0,0,0,0,1,0,0],
      [0,0,0,0,x,-1+x,0,0],[0,0,0,0,0,0,x,0],[0,x,0,0,0,0,0,-1+x]]],
 x->[[[-1,-1,0],[0,x,0],[0,1,-1]],[[-1,-1,0],[0,x,0],[0,1,-1]],[[-1+x,x,0],
      [1,0,0],[0,0,-1]],[[0,0,1],[0,-1,0],[x,0,-1+x]]],
    1,2,
 x->[[[x,0],[-1,-1]],[[x,0],[-1,-1]],[[0,1],[x,-1+x]],[[x,0],[-1,-1]]],
    3,7,4];
    if IsInt(r[i]) then return -x*r[r[i]](x^-1);else return r[i](x)*x^0;fi;
  elif [p,q,r]=[3,3,4] then 
    x:=-para[2][1]/para[2][2];
    m:=function(i)local f1,f2,f3,f5,f7,f8,f11,f13;
     f1:=x->x^0*[[[x,-1,0,0,0,0,0,0,0,0,1-x-x^2+x^3,0],
     [0,-1,0,0,0,0,0,0,0,0,0,0],[0,0,-1+x,0,x,0,-x,0,0,0,x-x^2,0],
     [0,0,0,-1+x,0,0,-x,0,0,0,x-x^2,0],[0,0,1,-1,0,0,0,0,0,0,
0,0],[0,0,0,0,0,0,0,0,0,0,0,-x],[0,0,0,-1,0,0,0,0,0,0,-1+x,0],[0,0,0,0,0,0,0,
-1,0,0,0,0],[0,0,0,0,0,0,0,0,-1+x,1,-1+x,0],[0,0,0,0,0,0,0,0,x,0,-1+x,0],[0,0,
0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,-1,0,0,0,0,0,-1+x]],[[0,x,0,0,0,0,0,0,0,0,0,
0],[1,-1+x,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,x,0,0,0,0,0],[0,0,0,0,x,0,0,1-x,
0,0,0,0],[0,0,0,1,-1+x,0,0,1-x,0,0,0,0],[0,0,0,0,0,0,0,1-x,-1+x,1,0,1-x+x^2],
[0,0,1,0,0,0,-1+x,0,0,0,0,0],[0,0,0,0,0,0,0,-1,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,
0,x],[0,0,0,0,0,x,0,x-x^2,-x,-1+x,0,x-x^2],[0,0,0,0,0,0,0,0,0,0,-1,0],[0,0,0,
0,0,0,0,0,1,0,0,-1+x]],[[0,-1+2*x-x^2,1-x,x,-x+x^2,0,0,0,-1+2*x-x^2,0,0,0],[0,
-1+x,1,0,0,0,0,0,-1+x,0,0,0],[0,x,0,0,0,0,0,0,-1+x,0,0,0],[1,-1+x,0,-1+x,0,0,
1-x,0,0,0,0,0],[0,0,0,0,-1+x,0,1,0,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,0,0],[0,0,0,
0,x,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,-1+x,0,0,0,-1],[0,0,0,0,0,0,0,0,-1,0,0,0],
[0,0,0,0,0,0,0,0,0,0,x,0],[0,0,0,0,0,0,0,0,0,1,-1+x,0],[0,0,0,0,0,0,0,-x,0,0,
0,0]],[[-1,0,0,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,
x,0,0,0],[0,0,0,0,0,x,0,0,0,0,1-x,x-x^2],[0,0,0,0,0,0,0,0,0,1,0,x],[0,0,0,1,0,
-1+x,-1+x,0,0,0,1-x,0],[0,0,0,0,0,0,0,0,0,0,0,x],[0,0,0,0,0,0,0,x,0,0,-1,0],
[0,0,1,0,0,0,0,0,-1+x,0,0,0],[0,0,0,0,x,0,-x,0,0,-1+x,0,0],[0,0,0,0,0,0,0,0,0,
0,-1,0],[0,0,0,0,0,0,1,0,0,0,0,-1+x]]];
    f2:=function(x,j)return [[[-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,x,x,x]],
  [[-1,0,0,0],[0,-1,0,0],[0,-j^2,x,1],[0,0,0,-1]],[[-1,0,0,0],[x,x,-j*x,1],
   [0,0,-1,0],[0,0,0,-1]],[[x,1,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]]];
   end;
    f3:=x->[[[x,-1,0,-x,0,0],[0,-1,0,0,0,0],[0,0,-1+x,0,1,x],[0,0,0,-1,0,0],
   [0,0,x,0,0,x],[0,0,0,0,0,-1]],[[-1,0,0,0,0,0],[0,0,0,0,x,1],[0,0,-1,0,0,0],
   [-1,0,-1,x,0,-1+x],[0,1,0,0,-1+x,1],[0,0,0,0,0,-1]],[[0,x,1,-1,-1,0],
   [1,-1+x,1,-1,-1,0],[0,0,-1,0,0,0],[0,0,0,-1,0,0],[0,0,0,0,-1,0],
   [0,0,0,1,1,x]],[[x,-1,0,0,1,x],[0,-1,0,0,0,0],[0,0,-1,0,0,0],[0,0,-1,x,1,x],
   [0,0,0,0,-1,0],[0,0,0,0,0,-1]]];
    f5:=x->[[[-1]],[[-1]],[[-1]],[[-1]]];
    f7:=function(x,j)return [[[-1,0,0,0,0,0],[x,x,0,0,0,0],[x,0,x,0,0,0],[0,0,0,
-1,0,0],[0,0,0,0,-1,0],[0,0,0,-j*x^2,x,x]],[[x,1,0,0,0,0],[0,-1,0,0,0,0],
[0,j^2,x,0,0,0],[0,0,0,-1,0,0],[0,0,0,x,x,1],[0,0,0,0,0,-1]],[[x,0,1,0,1,
0],[0,x,j*x,0,0,1],[0,0,-1,0,0,0],[0,0,0,x,1,-j^2*x^-1],[0,0,0,0,-1,0],
[0,0,0,0,0,-1]],[[-1,0,0,0,0,0],[0,-1,0,0,0,0],[0,0,-1,0,0,0],[0,0,0,x,0,0],
[x,0,0,0,x,0],[0,x,0,0,0,x]]];end;
    f8:=function(x,j)return [[[-1,0,0,0,0,0,0,0],[1,x,0,0,0,0,1,0],[1,0,x,0,
0,0,0,0],[0,0,0,-1,0,0,0,0],[0,0,0,0,-1,0,0,0],[0,0,0,-j*x,x,x,0,0],[0,0,0,
0,0,0,-1,0],[0,0,0,0,(j^2-j)*x,0,1,x]],[[x,x,0,0,0,0,-j^2,0],[0,-1,0,
0,0,0,0,0],[0,j,x,0,0,0,0,0],[0,0,0,-1,0,0,0,0],[0,0,0,1,x,1,0,0],[0,0,0,0,
0,-1,0,0],[0,0,0,0,0,0,-1,0],[0,0,0,0,0,-2*j^2-j,1,x]],[[x,0,x,0,x,0,0,
0],[0,x,j^2*x,0,0,1,0,0],[0,0,-1,0,0,0,0,0],[0,0,0,x,x,-j^2,0,0],[0,0,0,
0,-1,0,0,0],[0,0,0,0,0,-1,0,0],[0,0,0,0,0,0,x,x],[0,0,0,0,0,0,0,-1]],[[-1,0,0,
0,0,0,0,0],[0,-1,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0],[0,0,0,x,0,0,-j^2,0],[1,0,
0,0,x,0,0,0],[0,x,0,0,0,x,0,0],[0,0,0,0,0,0,-1,0],[0,0,(j^2-j)*x,0,0,0,1,x]]];
    end;
    f11:=x->[[[x,1,0],[0,-1,
0],[0,0,-1]],[[x,1,0],[0,-1,0],[0,0,-1]],[[-1,0,0],[x,x,1],[0,0,-1]],[[-1,0,
0],[0,-1,0],[0,x,x]]];
    f13:=x->[[[-1,0],[x,x]],[[-1,0],[x,x]],[[x,1],[0,-1]],[[-1,0],[x,x]]];
     r:=[f1(x),f2(x,E(3)),f3(x),f2(x,E(3)^2),f5(x),-x*f1(x^-1),f7(x,E(3)),
       f8(x,E(3)),f8(x,E(3)^2),-x*f7(x^-1,E(3)),f11(x),-x*f3(x^-1),f13(x),
       -x*f2(x^-1,E(3)^2),-x*f2(x^-1,E(3)),-x*f11(x^-1),-x*f5(x^-1)];
     return x^0*r[i];end;
     return m(i);
  elif [p,q,r]=[3,3,5] then 
     x:=-para[2][1]/para[2][2];
     m:=function(i)local r,f1,f2,f3,f4,f8,f9,f11,f12,f13,f17,f20,f23,f29;
     f1:=function(x)return
  x^0*[[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[E(3)^2-E(3)^2*x,E(3)^2-E(3)^2*x,0,0,0,E(3)^2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,
-E(3)^2*x+E(3)^2*x^2,0,E(3)^2-E(3)^2*x,0,0,0,E(3)^2,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0],[-E(3)^2-ER(-3)*x+E(3)*x^2,E(3)-E(3)*x+E(3)*x^2,0,0,0,
E(3)-E(3)*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,
E(3)^2*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,E(3),-1+x,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,
E(3)^2*x+ER(-3)*x^2-E(3)*x^3,0,E(3)-E(3)*x+E(3)*x^2,0,0,0,E(3)-E(3)*x,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,E(3),0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0],[0,0,E(3)^2-E(3)^2*x,0,0,0,0,0,0,0,E(3)^2-E(3)^2*x,0,0,
0,0,E(3)^2,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,E(3)^2*x,0,0,-1+x,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,E(3)*x,0,-1+x,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,-E(3)^2-ER(-3)*x+E(3)*x^2,0,0,0,0,0,0,0,E(3)-E(3)*x+E(3)*x^2,0,0,0,0,
E(3)-E(3)*x,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
E(3)*x,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2*x,
0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2,0,-1+x,0,0,0,0,0,
0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3),0,-1+x,0,0,0,0,0,0,0,0,0,
0],[0,0,E(3)^2-E(3)^2*x-E(3)^2*x^2+E(3)^2*x^3,0,E(3)-ER(-3)*x-E(3)^2*x^2,0,0,
0,-E(3)^2+E(3)^2*x,0,-E(3)+ER(-3)*x+E(3)^2*x^2,0,0,0,0,E(3)^2-E(3)^2*x,0,0,0,
0,-1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2,
0,0,0,0,0,0],[E(3)-E(3)*x,0,0,2*E(3)-E(3)*x^-1-E(3)*x,0,0,0,0,0,0,
-E(3)+2*E(3)*x-E(3)*x^2,0,0,0,0,-E(3)+E(3)*x,0,0,0,E(3)-E(3)*x,0,0,
E(3)-E(3)*x,0,E(3),0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
E(3)*x,0,-1+x,0,0,0,0,0,0],[-E(3)+ER(-3)*x+E(3)^2*x^2,0,x-2*x^2+x^3,
(1-3*ER(-3))/2+E(3)*x^-1+((1+3*ER(-3))/2)*x+E(3)^2*x^2,0,0,0,0,0,0,
E(3)-2*E(3)*x+2*E(3)*x^2-E(3)*x^3,0,0,0,0,E(3)-2*E(3)*x+E(3)*x^2,0,
-E(3)+E(3)*x,0,-E(3)+2*E(3)*x-E(3)*x^2,0,0,E(3)^2-E(3)^2*x+E(3)^2*x^2,0,
E(3)^2-E(3)^2*x,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-1+x,0,-E(3)^2+E(3)^2*x,0,0,
0,-1+x^-1,0,E(3)^2-E(3)^2*x,0,0,0,0,0,-1,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2*x,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,E(3),-1+x,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,-1+x,E(3)*x],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,E(3)^2,0]],[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0],[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0],[0,x,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0],[1-x,0,0,0,0,0,-1+x,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[-1+x^-1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,
x,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-2*x+x^2,1-2*x+x^2,
1-2*x+x^2,1-x,2-x^-1-x,1-x,0,0,1-x^-1,-1+x,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,
1-x,0,0,0,0,0,0,0,-1+x,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-2*x+x^2,-x+x^2,
1-2*x+x^2,1-x,1-x,0,0,0,0,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,
1-x,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[2-x^-1-x,1-x,0,0,0,
1-x,-1+x,-1+x,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,
0,0,x,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,x-x^2,-1+2*x-x^2,0,0,0,
x-x^2,0,0,0,0,0,-1+x,0,0,-1+x,0,x,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1+x,0,0,0,0,
0,0,0,0,0,0,0,0,0,-1+x,0,x,0,0,0,0,0,0,0,0,0,0],[2-x^-1-x,0,1-x,0,0,0,-1+x,0,
0,0,0,1-x,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1-x^-1,0,0,0,0,0,0,0,0,
0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
-1,0,0,0,0,0,0,0,0,0],[-1+3*x-3*x^2+x^3,-1+3*x-3*x^2+x^3,x-2*x^2+x^3,2-x^-1-x,
-2+x^-1+2*x-x^2,-1+2*x-x^2,0,1-2*x+x^2,-2+x^-1+x,0,0,-1+x,-1+x,2-x^-1-x,1-x,0,
0,0,1-x,0,0,-1+x,0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,1,0,0,0,0,0],[-3+x^-1+3*x-x^2,-1+2*x-x^2,0,1-2*x+x^2,0,-1+2*x-x^2,
1-2*x+x^2,1-2*x+x^2,-1+x,x-x^2,0,-1+2*x-x^2,0,1-x,1-x,0,-1+x,0,0,0,0,x,0,0,0,
0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,-1+x,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0],
[1-3*x+3*x^2-x^3,1-2*x+2*x^2-x^3,1-3*x+3*x^2-x^3,1-2*x+x^2,3-x^-1-3*x+x^2,
1-2*x+x^2,0,0,2-x^-1-x,0,-1+2*x-x^2,0,1-x,0,0,0,0,-1+x,0,0,1-x,0,1-2*x+x^2,0,
1-x,0,-1+x,x,0,0],[0,0,0,2-x^-1-x,0,-1+x,0,0,0,-1+x,0,0,0,0,0,-2+x^-1+x,0,
2-x^-1-x,0,1-x,-1+x^-1,0,-1+x,0,0,0,1,0,0,0],[-1+2*x-x^2,0,0,0,0,0,1-2*x+x^2,
-x+x^2,0,0,0,0,0,0,0,1-x,-1+x,0,0,0,0,0,1-x,0,0,1-x,0,0,0,x],[0,0,1-2*x+x^2,
3-x^-1-3*x+x^2,0,0,1-x,1-2*x+x^2,0,0,-1+x,0,0,2-x^-1-x,0,0,0,0,1-x,0,0,0,
2-x^-1-x,0,1-x^-1,1-x,0,0,1,-1+x]],[[0,-x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0],[-1,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,
0,0,0,0,0,0,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,x,0,-1+x,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,
0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,-1,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0,0],[0,0,0,0,0,
0,0,x,0,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,-1+x,
0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,
0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],[0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,x,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
-1+x,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,-1+x,0,0,0,
0,0,0],[0,x-x^2,0,0,0,0,0,0,0,1-2*x+x^2,-1+2*x-x^2,0,0,0,0,0,0,-1+x,0,0,1-x,0,
1-2*x+x^2,0,0,0,-1+x,x,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+x,0,0,0,0,0,
1-x,0,0,0,0,0,0,x],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,-1+x,
0,0,0],[1-x,0,0,2-x^-1-x,0,0,0,0,0,0,-1+2*x-x^2,0,0,0,0,-1+x,0,0,0,1-x,0,0,
1-x,0,1,0,0,-1+x,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,-1,0],[0,0,0,0,0,0,0,0,0,0,0,1-x,0,0,0,0,0,1-x^-1,0,0,0,0,0,0,0,1,0,0,0,
-1+x]],[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,
0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,x,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0],[1-x,0,0,0,0,0,-1+x,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,x,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-x,1-x,0,0,
0,1,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,x-x^2,0,0,0,0,
0,0,0,-x+x^2,0,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0,0,0,0,0],[0,0,
-x+x^2,0,1-x,0,0,0,1,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,
0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1+x,0,0,0,0,0,0,0,
0,0,0,0,0,0,-1+x,0,x,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,
-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,-1+x,0,0,0,
0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,
0,0,0],[0,0,1-x,0,0,0,0,0,0,0,1-x,0,0,0,0,1,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0],[0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,-1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,
0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1+x,0,0,0,0],[0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x],[0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,x,-1+x,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,1,0,0,-1+x]],[[0,0,-x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0],[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1,0,-1+x,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,x,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,
1-x,0,0,0,0,0,0,0,1-x,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,x,0,0,
-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,-1,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1-x,0,0,0,0,0,
-1+x,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],[1-x,0,0,0,0,0,-1+x,x,0,0,0,0,0,0,
0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,x,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,
0,0,0,0],[1-x,0,0,2-x^-1-x,0,0,0,0,0,0,-1+2*x-x^2,0,0,0,0,-1+x,0,0,0,1-x,0,0,
1-x,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,
0,0,0],[0,-x+x^2,0,0,0,0,1-x,0,0,0,0,0,0,0,x,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,
0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,-1+x,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0],[0,0,x-x^2,-1+2*x-x^2,0,0,0,x-x^2,0,0,0,
0,0,-1+x,0,0,-1+x,0,x,0,0,0,0,0,-1+x,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
x,0,0,0,0,-1+x,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
-1+x,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0],[0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1]]];end;
    f2:=function(q)return
    q^0*[[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-E(3)+ER(-3)*q+E(3)^2*q^2,
E(3)^2-E(3)^2*q,0,0,-E(3)+E(3)*q,0,ER(-3)-E(3)*q^-1+E(3)^2*q,
-E(3)^2+E(3)^2*q^-1+E(3)^2*q,0,0,2*E(3)-E(3)*q^-1-E(3)*q,0,0,0,0,0,0,0,0,0],
[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0],[0,0,0,0,0,0,0,0,0,0,E(3)^2,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,
0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],
[E(3)*q-E(3)*q^2,E(3)*q,0,0,0,0,E(3)-E(3)*q,E(3)-E(3)*q,0,0,E(3)-E(3)*q,0,0,0,
0,0,0,0,0,0],[0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,
0,0,E(3),0,0,0,0,0,0,0],[0,0,0,0,E(3)*q,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0],[0,
0,0,0,0,0,E(3)^2-E(3)^2*q,0,E(3)^2-E(3)^2*q,0,0,E(3)^2-E(3)^2*q,0,0,0,E(3)^2,
0,0,0,0],[0,0,0,0,0,0,0,0,0,E(3)^2*q,0,0,-1+q,0,0,0,0,0,0,0],[-E(3)+E(3)*q,0,
0,0,0,E(3)-E(3)*q,0,0,0,0,0,0,0,-1+q,-E(3),0,0,0,0,0],[1-q,0,0,0,0,-1+q,0,0,0,
0,0,0,0,-E(3)^2*q,0,0,0,0,0,0],[0,0,0,0,0,0,-E(3)^2-ER(-3)*q+E(3)*q^2,0,
-E(3)^2-ER(-3)*q+E(3)*q^2,0,0,E(3)-E(3)*q+E(3)*q^2,0,0,0,E(3)-E(3)*q,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3),0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
E(3)^2*q,-1+q,0,0],[E(3)*q-2*E(3)*q^2+E(3)*q^3,E(3)*q-E(3)*q^2,0,0,0,0,
E(3)-2*E(3)*q+E(3)*q^2,E(3)-2*E(3)*q+E(3)*q^2,0,0,E(3)-2*E(3)*q+E(3)*q^2,0,
-E(3)+E(3)*q,0,0,0,0,-E(3)+E(3)*q,0,E(3)^2*q],[0,0,0,0,0,0,0,-E(3)+E(3)*q,0,
E(3)-E(3)*q,0,0,0,0,0,0,E(3)-E(3)*q,0,E(3),-1+q]],[[-1,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0],[0,-1+q,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,
q^3-q^4,1-q,-1+q,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,0,0,0,0,
0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,q,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0],[1-2*q+q^2,0,
q^3-q^4,0,0,0,0,0,1-q,-1+q,0,0,1,1-q,0,0,0,0,0,0],[0,0,q^3-q^4,1-q,q,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],[1-2*q+q^2,0,
q^3-q^4,0,0,0,0,0,1-q,q,0,0,0,0,1-q,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,
-1,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,-q,-1+q,0,0,0,0,0],[0,0,0,0,0,0,0,0,
0,0,0,q,0,0,0,-1+q,0,0,0,0],[0,0,0,1-q,0,-1+q,0,0,0,0,0,2-q^-1-q,0,0,0,1-q^-1,
-1+q,1,0,0],[0,0,0,1-q,0,-1+q,0,0,0,0,0,1-q,0,0,0,0,q,0,0,0],[-1+2*q-q^2,0,0,
0,q-q^2,1-2*q+q^2,0,0,0,0,0,-1+q,0,0,-1+q,0,0,0,-1+q,q],[2-q^-1-q,0,
-q^3+2*q^4-q^5,-1+2*q-q^2,1-2*q+q^2,-2+q^-1+q,0,0,0,0,-1+q,0,0,-1+q,2-q^-1-q,
-1+q^-1,0,0,1,0]],[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,
0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0],[1-q,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,-q,0,0,0,0,0,
0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,
0,q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,1,-1+q,0,0,0,0,0,0,0,0,0,0],[0,0,0,
0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,
0,0,0,0],[1-q,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,q,0],[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1+q,0,0,0],[0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,q],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1+q,0],[0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1+q]],[[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,
0,0,0],[0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-q^-2,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0],[0,0,-q^3,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,
-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
[-q,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,
0,0,0,0,0],[0,0,0,0,0,q,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,-q,0,0,0],[0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0],[-1+q,0,0,0,
0,1-q,0,0,0,0,0,-1+q,0,0,-1,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,
0,0],[0,0,0,0,0,0,1-q,0,1-q,0,0,1-q,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1-q,0,1-q,0,
0,-q,0,0,0,0,0,0,0,0],[1-q,0,0,0,0,-1+q,0,0,0,0,0,0,0,q,1-q,-1+q,0,0,0,0],[0,
0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,-1+q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,
0,-1+q,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,-1]],[[-1+q,0,-q^3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[1-2*q+q^2,0,q^3-q^4,0,0,0,0,0,1-q,-1+q,0,0,1,1-q,0,0,0,0,0,0],[-q^-2,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-q,
0,0,0,0,-1+q,0,0,0,0,0,0,0,0,1,0,0,0,0,0],[0,0,0,1,0,-1+q,0,0,0,0,0,0,0,0,0,0,
0,0,0,0],[0,0,0,0,0,0,0,0,-q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,q,0,0,
0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,-1+q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,
1,0,-1+q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0,0,0,0,0],[0,0,
0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0],[q-q^2,q,0,0,0,0,1-q,1-q,0,0,1-q,0,-1+q,
0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,-1,0,0,-1+q,0,0,0,0,0,0],[0,0,q^3-q^4,1-q,
q,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,
0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,-1,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,-1]]];end;
  f3:=function(q,j)return
  q^0*[[[-1,0,0,0,0],[0,-1,0,0,0],[1,0,0,0,-1],[1+j*q,0,1+j*q,-1,-1-j*q],
[-q,0,-q,0,-1+q]],[[-1,0,0,0,0],[0,-1,0,0,0],[1,0,0,-1,0],[-q,0,-q,-1+q,0],
[-j,0,-j,j,-1]],[[-1,0,-1,0,0],[0,-1,1,0,0],[0,0,q,0,0],[0,0,0,-1,0],
[0,0,0,0,-1]],[[q,0,0,0,0],[-1,-1,0,0,0],[-q,0,-1,0,0],[1,0,0,-1,0],[1,0,0,0,
-1]],[[0,1,0,0,0],[q,-1+q,0,0,0],[0,0,-1,0,0],[1,1,0,-1,0],[1,1,0,0,-1]]];end;
   f4:=function(q,j)return 
   q^0*[[[-1+q,0,0,q,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,1],[0,0,-1,0,0,0,0,0,0,0],[1,0,
0,0,0,0,0,0,0,0],[0,0,j^2*q+(-j^2+j)*q^2-j*q^3,0,j-j*q,0,
-j+j*q-j*q^2,0,0,0],[0,0,0,0,0,-1,0,0,0,0],[0,0,j^2*q-j^2*q^2,
0,-j^2,0,j^2-j^2*q,0,0,0],[-q+q^2,0,j^2*q-2*j^2*q^2+j^2*q^3,
-q+q^2,-j^2+j^2*q,0,-j+(-j^2+j)*q+j^2*q^2,-1,0,0],[0,q-q^2,
j^2*q-2*j^2*q^2+j^2*q^3,0,-j^2+j^2*q,0,
-j+(-j^2+j)*q+j^2*q^2,0,-1,q-q^2],[0,q,0,0,0,0,0,0,0,-1+q]],[[0,0,
j^2*q-j^2*q^2,j^2*q,0,0,0,0,0,0],[0,-1+q,-q+q^2,0,0,0,0,0,0,j],[0,
0,-1,0,0,0,0,0,0,0],[j,0,q-q^2,-1+q,0,0,0,0,0,0],[0,0,0,0,-1+q,0,-q,0,0,0],
[q-q^2,-j^2*q+j^2*q^2,j^2*q-j^2*q^2-j^2*q^3+j^2*q^4,
j^2*q^2-j^2*q^3,0,-1,0,0,0,-1+q],[0,0,0,0,-1,0,0,0,0,0],[0,0,0,0,0,0,0,
-1,0,0],[0,0,0,0,0,0,0,0,-1,0],[0,j^2*q,-j^2*q+j^2*q^2,0,0,0,0,0,0,
0]],[[-1,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,
q,-1+q,0,0,0,0,0,0],[0,0,0,0,0,0,0,-1,0,0],[0,0,0,0,0,-1+q,0,0,0,q],[0,0,0,0,
0,0,-1,0,0,0],[0,0,0,0,-q,0,0,-1+q,0,0],[0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,1,0,
0,0,0]],[[0,0,0,0,0,0,0,0,0,1],[0,-1+q,0,q,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,
0],[0,1,0,0,0,0,0,0,0,0],[0,0,0,0,-1,0,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0],[0,0,0,
0,0,0,-1,0,0,0],[0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,0,-q,-1+q,0],[q,0,0,0,0,0,
0,0,0,-1+q]],[[-1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0],[0,0,-1,0,0,0,0,0,
0,0],[0,0,0,-1,0,0,0,0,0,0],[0,-q+q^2,-q^2+q^3,0,-1+q,0,0,0,0,j*q],[0,0,0,
-j^2*q+j^2*q^2,0,0,j^2-j^2*q,-j^2,0,0],[0,-q,0,0,0,0,-1+q,0,0,
0],[0,-q+q^2,0,q^2-q^3,0,-j*q,0,-1+q,0,0],[0,0,0,0,0,0,0,0,-1,0],[0,0,
-j^2*q+j^2*q^2,0,j^2,0,-j^2+j^2*q,0,0,0]]];end;
   f11:=function(q,j)return 
   q^0*[[[q,0,0,0,0,0,0,0,0,0],[-j^2+j^2*q,j^2-j^2*q,j^2,0,0,0,0,0,0,
0],[-j+(-j^2+j)*q+j^2*q^2,j-j*q+j*q^2,j-j*q,0,0,0,
0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,q,0,-1+q,0,0,0,0],
[0,0,0,0,q,0,-1+q,0,0,0],[0,0,0,0,0,0,0,-1,0,0],[0,0,0,0,0,0,0,0,-1,0],[0,0,0,
0,0,0,0,0,0,-1]],[[q,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,q,-1+q,0,0,0,
0,0,0,0],[1-q,0,0,-1+q,0,j,0,0,0,0],[1-q,0,0,0,-1+q,0,j,0,0,0],
[-j^2*q+j^2*q^2,0,0,j^2*q,0,0,0,0,0,0],[-j^2*q+j^2*q^2,0,0,0,
j^2*q,0,0,0,0,0],[-j^2+2*j^2*q-j^2*q^2,j^2-j^2*q,
j^2-j^2*q,-j^2*q+j^2*q^2,0,-1+q,0,-1,0,0],[0,0,0,
j^2*q-j^2*q^2,-j^2*q+j^2*q^2,1-q,-1+q,0,-1,0],
[-j^2+2*j^2*q-j^2*q^2,j^2-j^2*q,j^2-j^2*q,0,
-j^2*q+j^2*q^2,0,-1+q,0,0,-1]],[[0,-1,0,0,0,0,0,0,0,0],[-q,-1+q,0,0,0,0,
0,0,0,0],[0,0,q,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0],[0,0,0,0,-1,0,0,0,0,0],
[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,q,0,-1+q,0,0],[0,0,0,0,
0,0,0,0,-1,0],[0,0,0,0,0,0,q,0,0,-1+q]],[[-1,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,
0,0,0,0],[-1+q,0,0,1-q,0,-j,0,0,0,0],[0,-q,0,-1+q,0,0,0,0,0,0],[0,0,0,0,-1,
0,0,0,0,0],[j^2*q-j^2*q^2,-j^2*q+j^2*q^2,-j^2*q,0,0,-1+q,0,0,0,
0],[0,0,0,0,0,0,-1,0,0,0],[0,0,0,0,0,0,0,q,0,0],[0,0,0,0,0,0,0,0,-1+q,-q],[0,
0,0,0,0,0,0,0,-1,0]],[[-1,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0],[0,0,-1,0,
0,0,0,0,0,0],[0,0,0,0,-1,0,0,0,0,0],[0,0,0,-q,-1+q,0,0,0,0,0],[0,0,0,0,0,0,-1,
0,0,0],[0,0,0,0,0,-q,-1+q,0,0,0],[0,0,0,0,0,0,0,0,0,-1],[0,0,0,0,0,0,0,0,q,0],
[0,0,0,0,0,0,0,-q,0,-1+q]]];end;
   f12:=function(q,j)return
q^0*[[[0,0,0,0,-j^2,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,-j^2*q,0,0,0,0,0,0,0,0,
0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-j^2,0,0,0,0,0,0,0,0],
[-j*q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0],[0,-j,0,0,0,-1+q,0,0,0,0,0,0,0,0,
0],[0,0,0,-j*q,0,0,-1+q,0,0,0,0,0,0,0,0],[-j^2+j^2*q,0,0,0,
j-j*q,0,0,-1,0,0,0,0,0,0,0],[0,0,0,j^2-j^2*q,0,0,-j+j*q,0,
-1,0,0,0,0,0,0],[0,2*j^2-j^2*q^-1-j^2*q,0,0,0,j-2*j*q+j*q^2,
0,0,0,-1,0,0,0,0,0],[0,-1+2*q-q^2,-q,0,0,1-2*q+q^2,0,0,0,q,q,0,0,0,0],[0,0,0,
0,0,0,0,0,0,0,0,-1,0,0,0],[j-j*q,0,0,-q+q^2,j^2+2*j+q^-1-j*q,0,
1-q,-q,-q,0,0,0,q,0,0],[0,-2*j^2+j^2*q^-1+j^2*q,0,0,0,
-j+2*j*q-j*q^2,0,0,0,0,0,0,0,-1,0],[0,1-2*q+q^2,0,0,0,-1+2*q-q^2,0,0,
0,0,0,q,0,q,q]],[[-1+q,0,0,0,-1,0,0,0,0,0,0,0,0,0,0],[0,-1+q,0,0,0,-q,0,0,0,0,
0,0,0,0,0],[-1+q,0,-1,0,-1+q^-1,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1+q,0,0,-1,0,0,0,
0,0,0,0,0],[-q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,-q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0],[0,0,0,
-q+q^2,0,0,1-q,0,-1,0,0,0,0,0,0],[0,-q+q^2,0,0,0,q-q^2,0,0,0,-1,0,0,0,0,0],[0,
-1+2*q-q^2,-q,1-q,1-q,-q+q^2,1-q,0,0,q,q,0,0,0,0],[0,0,0,-1+q,0,0,-1+q^-1,0,0,
0,0,-1,0,0,0],[-j^2+j^2*q,1-q^-1,0,1-2*q+q^2,-j^2+j^2*q,-1+q,1-q,
-q,-q,0,0,0,q,0,0],[0,q-q^2,0,0,0,-q+q^2,0,0,0,0,0,0,0,-1,0],[1-q^-1,
2+q^-2-2*q^-1-2*q+q^2,0,2-q^-1-q,1-q^-1,-2+q^-1+2*q-q^2,1-q^-1,0,0,0,0,q,0,q,
q]],[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,
-1+q,0,-1,0,0,0,0,0,0,0,0,0,0],[q,-1,0,q,0,0,0,0,0,0,0,0,0,0,0],[0,0,-q,0,0,0,
0,0,0,0,0,0,0,0,0],[0,1-2*q+q^2,0,0,0,-1+2*q-q^2,0,0,0,-q,0,0,0,0,0],[0,0,0,0,
0,0,0,0,0,0,q,0,0,0,0],[0,0,-j^2-2*j-q^-1+j*q,0,
-j+q^-2+(j^2+2*j)*q^-1,0,0,-1,0,0,0,0,0,0,0],
[-j+(j^2+2*j)*q+q^2,1-q,0,0,-j^2-2*j-q^-1+j*q,0,-1+q,q,q,0,
1-q,0,-1,0,0],[0,3-q^-1-4*q+3*q^2-q^3,0,0,0,-2+3*q-3*q^2+q^3,0,0,0,-q+q^2,0,0,
0,0,0],[0,0,0,0,0,0,1,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q],[0,
0,0,0,0,0,0,0,0,0,0,0,-1,0,0],[0,-4+q^-1+6*q-4*q^2+q^3,0,0,0,2-5*q+4*q^2-q^3,
0,0,0,-1+2*q-q^2,0,0,0,-1,0],[0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1+q]],[[-1,0,0,0,
0,0,0,0,0,0,0,0,0,0,0],[0,0,0,q^2,0,0,0,0,0,0,0,0,0,0,0],[-j+j*q,0,0,0,
-j^2-2*j-q^-1+j*q,0,0,q,0,0,0,0,0,0,0],[0,q^-1,0,-1+q,0,0,0,0,0,0,0,
0,0,0,0],[0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,q,0,0,0,0,0,0,0,0],[0,
0,0,0,0,1,-1+q,0,0,0,0,0,0,0,0],[-j+j*q,0,1,0,
-j^2-2*j-q^-1+j*q,0,0,-1+q,0,0,0,0,0,0,0],[0,-1+q,0,0,0,1-q,0,0,-1+q,
-1,0,0,0,0,0],[0,0,0,-q^2+q^3,0,0,q-q^2,0,-q,0,0,0,0,0,0],[-j+j*q,0,0,
q-q^2,-j^2-2*j-q^-1+j*q,0,-1+q,q,q,0,0,0,-1,0,0],[0,2-q^-1-q,0,0,0,
-2+q^-1+q,0,0,0,0,0,-1+q,0,-1,0],[0,1-2*q+q^2,q,0,0,-1+2*q-q^2,0,0,0,-q,-q,0,
-1+q,0,0],[0,0,0,-q+2*q^2-q^3,0,0,1-2*q+q^2,0,0,0,0,-q,0,0,0],[0,0,0,0,0,0,0,
0,0,0,0,0,0,0,-1]],[[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],[q,-1+q,0,0,0,0,0,0,0,0,0,
0,0,0,0],[0,-1+2*q-q^2,0,0,0,1-2*q+q^2,0,0,0,q,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,
0,0,0,0,0,0],[0,0,0,0,0,q,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,-1+q,0,0,0,0,0,0,0,0,
0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0],[0,2*j^2+j-j^2*q^-1+q,0,0,0,
-j^2+j^2*q,0,0,0,0,0,0,0,1,0],[0,0,0,1-q,0,0,1-q^-1,0,0,0,0,1,0,0,0],
[1-2*q+q^2,0,1,0,2-q^-1-q,0,0,0,0,-1+q,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,-1,0,0,
0,0],[0,0,0,q-q^2,0,0,-1+q,0,q,0,0,-1+q,0,0,0],[0,1-2*q+q^2,0,0,0,-1+2*q-q^2,
0,0,0,0,0,q,0,q,q],[j^2+(-2*j^2-j)*q-q^2,0,0,0,j^2-j^2*q,0,0,q,
0,0,0,0,0,-1+q,0],[j-j*q,0,0,-q+q^2,j^2+2*j+q^-1-j*q,0,1-q,-q,
-q,0,0,0,1,0,-1+q]]];end;
f13:=function(q,j) return 
q^0*[[[j-j*q,-1+q-q^2,-j^2+2*j^2*q-j^2*q^2,
-j+(-j^2+j)*q+j^2*q^2,-j^2+(j^2-j)*q+j*q^2,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0],[-1,j^2-j^2*q,-1+q,1-q,-1+q,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-q,q,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0],[0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,-q,q,0,0,
0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[-q,0,
-q+q^2,0,-q+q^2,1-q,0,j^2-j^2*q,0,0,0,0,0,0,0,q,0,0,0,0],[0,0,1-q,-1+q,
0,0,-1+q,0,j^2-j^2*q,0,1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,q,0,0,0,
0,0,0,0,0,0,0],[0,0,-j^2+2*j^2*q-j^2*q^2,
-j+(-j^2+j)*q+j^2*q^2,0,0,j^2+(-j^2+j)*q-j*q^2,0,
1-q+q^2,0,j-j*q,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,j^2,0,
0,0,0,0,0],[0,0,0,0,1-q,-1+q,1-q,0,0,0,0,0,j^2-j^2*q,0,j^2*q,0,0,0,0,
0],[0,0,0,0,0,0,0,0,0,0,0,j*q,0,-1+q,0,0,0,0,0,0],[0,0,0,0,2-q^-1-q,
2*j^2+j-j^2*q^-1+q,-2*j^2-j-q^-1+j^2*q,0,0,0,0,0,
-j+j*q^-1+j*q,0,j-j*q,0,0,0,0,0],[0,-1+q-q^2,-q+q^2,
-j+(-j^2+j)*q+j^2*q^2,1-q,j^2-j+j*q^-1-j^2*q,0,
-1+q^-1+q,0,0,0,0,0,0,0,j-j*q,0,0,0,0],[0,0,0,0,
j^2-j-j^2*q^-1+j*q,0,-2*j+j*q^-1+j*q,0,0,
-j^2+j-j*q^-1+j^2*q,0,0,0,0,0,0,j-j*q,1-q^-1-q,0,0],[0,0,0,
0,-1+q,0,-1+q,0,0,1-q,0,0,0,0,0,0,-q,j^2-j^2*q,0,0],[0,
j*q+(j^2-j)*q^2-j^2*q^3,-q+2*q^2-q^3,0,1-2*q+q^2,0,1-q-q^2+q^3,
-j+(-j^2+j)*q+j^2*q^2,j*q+(j^2-j)*q^2-j^2*q^3,0,
-q+q^2,0,-j+(-j^2+j)*q+j^2*q^2,0,j^2*q-j^2*q^2,q-q^2,0,0,-1,
0],[1-q,j+(j^2-j)*q-j^2*q^2,0,0,-1+q^-1-q+q^2,0,-1+q^-1-q+q^2,0,
j+(j^2-j)*q-j^2*q^2,2-q^-1-q,-1+q,-j+j*q^-1,0,-1+q^-1,0,0,
1-q,j^2-j+j*q^-1-j^2*q,0,-1]],[[-1+q,-j^2*q,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0],[-j,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-q,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,
0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,-q,q,0,0,0,0,0,0,0,0,0,0,0,0,
0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[-j*q,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,j*q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,j,0,0,0,0,0,0,0,0,0],[0,0,0,0,
0,0,-1,0,0,q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,j^2*q,0,-1+q,0,0,0,0,0,
0,0,0,0],[0,0,-1+q,0,0,0,1-q,0,0,-1+q,0,-1+q,0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,
0,0,0,0,0,0,0,q,0,0,0,0,0],[0,0,-1+q,0,0,0,0,0,0,q-q^2,0,q,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1+q,0,0,0,0,0],[0,-j^2*q,0,0,0,0,0,j^2,0,
0,0,0,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,-j^2,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-j*q,0,0,0],[0,j-j*q,0,0,0,0,0,
j-j*q^-1,j-j*q,0,j^2-j^2*q,0,j-j*q^-1,0,-j+j*q,
-j^2+j^2*q,0,0,-1,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1]],[[-1,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,-1+q,0,0,j,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0],[0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,1,-1+q,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0],[0,j^2*q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,
0,0,j^2,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,j^2*q,0,0,0,0,0,0,0,0,
0,0,0],[0,0,0,0,0,j*q,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,j,0,
-1+q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],[0,0,0,
0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,q,0,-1+q,0,0,0,0,0,0,0,
0],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0],[0,0,-q+q^2,q-q^2,0,0,q-q^2,0,
-j^2*q+j^2*q^2,0,-q,0,0,q,0,0,0,0,0,0],[0,
j+(j^2-j)*q-j^2*q^2,-1+2*q-q^2,0,-2+q^-1+q,0,-1+q^-1-q+q^2,
-j^2+j-j*q^-1+j^2*q,j+(j^2-j)*q-j^2*q^2,0,-1+q,0,
-j^2+j-j*q^-1+j^2*q,0,j^2-j^2*q,1-q,0,0,-1,0],[-1,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,j^2],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0],[-1+q,
-j+3*j*q-3*j*q^2+j*q^3,j-3*j*q+3*j*q^2-j*q^3,
-q+2*q^2-q^3,j^2+3*j-j*q^-1+(-2*j^2-4*j)*q+(j^2+2*j)*q^2,
-j^2+(2*j^2+j)*q+q^2,
j^2+2*j-j*q^-1+(-2*j^2-j)*q+(j^2-j)*q^2+j*q^3,
-2*j^2-3*j-q^-1+(j^2+3*j)*q-j*q^2,
-j+(j^2+3*j)*q+(-2*j^2-3*j)*q^2-q^3,0,
-j^2+(j^2-j)*q+j*q^2,0,2-q^-1-2*q+q^2,0,-1+q-q^2,
-j+(-j^2+j)*q+j^2*q^2,0,0,j-j*q,0],[0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,j*q,0,0,-1+q]],[[q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0,0,0],[0,0,0,
0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,q,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0],[0,0,0,q,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,
0,0,0,0,0,0,0,0,0],[0,q,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,
0,0,0,0,0,0,1,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0],[0,0,
0,0,0,0,0,0,0,0,0,0,0,0,j^2*q,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,-j^2,0,0],[0,0,0,0,0,0,0,0,q,0,0,0,-1+q,0,0,0,0,0,0,0],[0,0,0,0,1-q,0,
1-q,0,0,-1+q,0,0,0,0,0,0,q,-j^2+j^2*q,0,0],[0,0,0,0,0,0,0,0,0,0,j,0,
0,0,-1+q,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0],[0,0,-1+q,0,0,
0,1-q,0,0,-1+q,0,-1+q,0,1,0,0,-1+q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-j*q,0,0,
0,0,0,-1+q,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0],[0,
-j*q+(-j^2+j)*q^2+j^2*q^3,q-2*q^2+q^3,0,-1+2*q-q^2,0,-1+q+q^2-q^3,
j+(j^2-j)*q-j^2*q^2,-j*q+(-j^2+j)*q^2+j^2*q^3,0,q-q^2,
0,j+(j^2-j)*q-j^2*q^2,0,-j^2*q+j^2*q^2,-q+q^2,0,0,q,q]],[[0,
0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,q,0,
0,0,0,0,0,0,0,0,0],[0,0,0,0,q,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,
0,0,0,0,0,j*q,0,0,0,0,0,0,0,0],[0,q,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,
0],[0,0,0,0,0,1,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0],[-q,0,0,0,0,0,0,0,0,0,-1+q,0,
0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,j^2,0,0,0,-1+q,0,0,0,0,0,0,0,0],[0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0],[-q,0,-q+q^2,0,-q+q^2,1-q,0,j^2-j^2*q,0,
0,0,0,0,-1+q,0,q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,j*q,0,0,0],[0,0,
-1+q,0,0,0,1-q,0,0,-1+q,-1,-1+q,0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,
j^2,0,-1+q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,-1+q,0,0],[q-q^2,
j*q+(j^2-j)*q^2-j^2*q^3,0,0,1-q-q^2+q^3,0,1-q-q^2+q^3,0,
j*q+(j^2-j)*q^2-j^2*q^3,-1+2*q-q^2,-q+q^2,j-j*q,0,1-q,0,0,
q-q^2,j+(j^2-j)*q-j^2*q^2,0,-q],[0,
j+(j^2-j)*q-j^2*q^2,-1+2*q-q^2,0,-2+q^-1+q,0,-1+q^-1-q+q^2,
-j^2+j-j*q^-1+j^2*q,j+(j^2-j)*q-j^2*q^2,0,-1+q,0,
-j^2+j-j*q^-1+j^2*q,0,j^2-j^2*q,1-q,0,0,-1,-1+q]]];end;
  f17:=function(q)return
q^0*[[[-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[-q^4,q^3,-q^2,q]],[[-1,0,0,0],[0,-1,0,0],
[0,0,-1,0],[-q^4,q^3,-q^2,q]],[[-1,0,0,0],[0,-1,0,0],[0,0,0,1],[0,0,q,-1+q]],
[[-1,0,0,0],[0,0,1,0],[0,q,-1+q,0],[0,0,0,-1]],[[0,1,0,0],[q,-1+q,0,0],[0,0,
-1,0],[0,0,0,-1]]];end;
  f20:=function(q,j)return
q^0*[[[q,0,-q^2,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0],[0,0,
0,-1,0,0,0,0,0,0],[0,0,0,-q^2,q,0,0,0,0,0],[0,-q^2,0,0,0,q,0,0,0,0],
[-j*q+(-j^2+j)*q^2+j^2*q^3,j^2*q+(-j^2+j)*q^2-j*q^3,
-j^2*q^2+2*j^2*q^3-j^2*q^4,j^2*q-2*j^2*q^2+j^2*q^3,
j+(j^2-j)*q-j^2*q^2,0,j-j*q,-j+j*q-j*q^2,0,0],
[j^2*q-j^2*q^2,j^2*q-j^2*q^2,-j^2*q^2+j^2*q^3,
j^2*q-j^2*q^2,-j^2+j^2*q,0,-j^2,j^2-j^2*q,0,0],
[-j*q^2+(-j^2+j)*q^3+j^2*q^4,-j*q^2+2*j*q^3-j*q^4,
j^2*q^2+(-2*j^2+j)*q^3+(2*j^2-j)*q^4-j^2*q^5,
j^2*q^2+(-j^2+j)*q^3-j*q^4,0,j+(j^2-j)*q-j^2*q^2,0,
-j*q+j*q^2-j*q^3,j-j*q,(j-j*q+j*q^2)/q],[0,0,
-j^2*q^3+j^2*q^4,0,-j^2*q^2+j^2*q^3,j^2*q-j^2*q^2,
-j^2*q^2,0,j^2*q,j^2-j^2*q]],[[q,0,-q^2,0,0,0,0,0,0,0],[0,-1,0,0,
0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0],[0,0,0,-q^2,q,0,0,
0,0,0],[0,-q^2,0,0,0,q,0,0,0,0],[0,0,0,0,0,0,-1+q,-q,0,0],[0,0,0,0,0,0,-1,0,0,
0],[0,0,0,0,0,0,0,-q^2,-1+q,1],[0,0,0,0,0,0,-q^2,0,q,0]],[[-1+q,0,q,0,0,0,0,0,
0,0],[0,0,0,0,0,0,0,-1,0,0],[1,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,
0,q,-1+q,0,0,0,0,0],[0,0,0,0,0,0,0,-q,0,1/q],[0,0,0,0,0,0,-1,0,0,0],[0,-q,0,0,
0,0,0,-1+q,0,0],[0,0,0,0,0,0,-q^2,0,q,0],[0,-q^3,0,0,0,q^2,0,0,0,-1+q]],[[q,0,
-q^2,0,0,0,0,0,0,0],[0,-1+q,0,q,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0],[0,1,0,0,
0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,q,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,
1,0],[0,0,0,0,0,0,0,q,0,-1/q],[0,0,0,0,0,0,q,0,-1+q,0],[0,0,0,0,0,0,0,0,0,
-1]],[[0,0,0,0,1,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,
q,-1+q,0,0,0,0,0,0],[q,0,0,0,-1+q,0,0,0,0,0],[0,-q^2,0,0,0,q,0,0,0,0],[0,0,0,
0,0,0,-1,0,0,0],[0,0,0,0,0,0,0,-1,0,0],[0,0,0,0,0,0,-q^2,0,q,0],[0,0,0,0,0,0,
0,-q^2,0,q]]];end;
   f23:=function(q)return
q^0*[[[-1,0,0,0,0],[0,-1,0,0,0],[q^3,-q^2,q,0,0],[0,0,0,-1,0],[-q^3,0,0,-q^2,q]],
[[-1,0,0,0,0],[0,-1,0,0,0],[q^3,-q^2,q,0,0],[0,0,0,-1,0],[-q^3,0,0,-q^2,q]],
[[-1,0,0,0,0],[0,0,1,0,0],[0,q,-1+q,0,0],[0,0,0,0,1],[0,0,0,q,-1+q]],[[0,1,0,
0,0],[q,-1+q,0,0,0],[0,0,-1,0,0],[0,0,0,-1,0],[-q^4,q^3,-q^2,-q^2,q]],[[-1,0,
0,0,0],[0,0,0,1,0],[0,0,0,0,1],[0,q,0,-1+q,0],[0,0,q,0,-1+q]]];end;
   f29:=function(q)return
q^0*[[[q,q^4,0,0,-q^2,0],[0,-1,0,0,0,0],[0,0,-1,0,0,0],[0,q^3,-q^2,q,0,0],[0,0,0,
0,-1,0],[0,0,q^4,0,-q^3,q]],[[q,q^4,0,0,-q^2,0],[0,-1,0,0,0,0],[0,0,-1,0,0,0],
[0,q^3,-q^2,q,0,0],[0,0,0,0,-1,0],[0,0,q^4,0,-q^3,q]],[[-1+q,0,0,0,q,0],[0,-1,
0,0,0,0],[0,0,0,1,0,0],[0,0,q,-1+q,0,0],[1,0,0,0,0,0],[0,0,0,0,0,q]],[[0,0,0,
0,0,1],[0,0,1,0,0,0],[0,q,-1+q,0,0,0],[0,0,0,-1,0,0],[0,0,0,0,q,0],[q,0,0,0,0,
-1+q]],[[-1+q,0,0,q,0,0],[0,q,0,0,0,0],[0,0,0,0,1,0],[1,0,0,0,0,0],[0,0,q,0,
-1+q,0],[0,0,0,0,0,-1]]];end;
    f9:=function(q) return
q^0*[[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0],[-q,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-q,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[E(3)^2-2*E(3)^2*q+E(3)^2*q^2,E(3)+
(E(3)^2-E(3))*q-E(3)^2*q^2,E(3)^2+(-E(3)^2+E(3))*q-E(3)*q^2,0,0,0,0,
E(3)-E(3)*q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
E(3)-E(3)*q+E(3)*q^2,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,E(3)^2,0,0,0,0,0],[-1+2*q-q^2,1-2*q+q^2,-1+2*q-q^2,-1+q,0,
-1+q,0,1-q,0,E(3)-E(3)*q,0,0,0,E(3)*q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,1-q,0,1-2*q+q^2,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2*q,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,q,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-3+q^-1+3*q-q^2,3-q^-1-3*q+
q^2,-3+q^-1+3*q-q^2,-2+q^-1+q,0,E(3)^2+2*E(3)+q^-1-E(3)*q,0,2-q^-1-q,0,
-E(3)^2+E(3)^2*q^-1+E(3)^2*q,0,0,0,E(3)^2-E(3)^2*q,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,-E(3)^2-2*E(3)+E(3)*q^-1-q,0,2-q^-1-2*q+q^2,0,0],[-E(3)+
E(3)*q,E(3)-E(3)*q,0,-2*E(3)+E(3)*q^-1+E(3)*q,0,0,0,0,0,0,0,-E(3)+E(3)*q,0,0,
E(3)-E(3)*q,0,0,0,E(3)^2,0,0,0,E(3)-E(3)*q,0,0,0,0,0,0,0,0,0,0,0,0,
2*E(3)-E(3)*q^-1-E(3)*q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,-E(3)^2+E(3)^2*q,0,
E(3)^2-E(3)^2*q,0,0,E(3)^2-E(3)^2*q,0,0,0,0,0,0,E(3)^2*q,0,0,0,0,0,E(3)^2,0,0,
0,0,0,0,-E(3)^2+E(3)^2*q,0,0,0,0],[0,0,0,0,0,E(3)^2-E(3)^2*q^-1,0,0,0,0,0,
-E(3)^2+E(3)^2*q,0,0,0,0,E(3)^2-E(3)^2*q,0,0,E(3)^2,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,E(3)^2-E(3)^2*q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,q,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-2*q+q^2,E(3)+(-E(3)^2-2*E(3))*q-q^2,0,
3-q^-1-3*q+q^2,0,0,0,0,0,0,-1+q,1+(E(3)^2+2*E(3))*q-E(3)*q^2,0,0,E(3)-E(3)*q+
E(3)*q^2,0,0,0,E(3)^2-E(3)^2*q,0,0,0,-1+2*q-q^2,0,0,0,0,0,0,0,0,0,0,0,0,
E(3)^2+3*E(3)-E(3)*q^-1+(-2*E(3)^2-3*E(3))*q-q^2,0,0,0,0],[0,0,0,0,0,
-2*E(3)^2+E(3)^2*q^-1+E(3)^2*q,0,0,0,0,0,E(3)^2+(-E(3)^2+E(3))*q-E(3)*q^2,0,0,
0,0,E(3)-E(3)*q+E(3)*q^2,0,0,E(3)-E(3)*q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)+
(E(3)^2-E(3))*q-E(3)^2*q^2,0,0,0],[0,0,-E(3)+E(3)*q,0,E(3)-E(3)*q,0,-E(3)+
E(3)*q,0,0,0,0,0,0,0,0,0,E(3)-2*E(3)*q+E(3)*q^2,0,0,E(3)-E(3)*q,E(3)-E(3)*q,0,
0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1+q,0,0,-2+q^-1+q,0,0,0,0,0,0,2-q^-1-q,
-1+q,-2+q^-1+q,0,0,-2+q^-1+q,0,1-q,0,0,0,E(3)-E(3)*q,1-q,0,0,E(3),0,0,-1+q^-1,
0,0,0,0,0,0,2-q^-1-q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,E(3),0,0,0,0,0,0,0,0,0,0,0,
-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,1-q^-1,0,0,0,0,0,-1+q,0,0,
0,0,0,1-q,0,1,0,0,0,E(3)^2-E(3)^2*q,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],[0,0,
-2*E(3)^2+E(3)^2*q^-1+E(3)^2*q,0,2*E(3)^2+E(3)+q^-1-E(3)^2*q,0,-2*E(3)^2-E(3)+
E(3)^2*q^-1-q,0,0,0,0,E(3)-2*E(3)*q+E(3)*q^2,0,0,0,0,
2*E(3)^2-E(3)^2*q^-1-2*E(3)^2*q+E(3)^2*q^2,0,0,2*E(3)^2-E(3)^2*q^-1-E(3)^2*q,
-1+q^-1+q,0,0,0,E(3)^2-E(3)^2*q,0,0,0,0,0,0,0,0,0,0,0,-2*E(3)+E(3)*q^-1+
E(3)*q,0,0,0],[1+(E(3)^2+2*E(3))*q-E(3)*q^2,0,0,-2*E(3)^2-3*E(3)-q^-1+(E(3)^2+
3*E(3))*q-E(3)*q^2,0,0,0,0,0,0,-3+q^-1+3*q-q^2,-E(3)+2*E(3)*q-E(3)*q^2,
3-q^-1-3*q+q^2,0,0,2-q^-1-2*q+q^2,0,E(3)+(-E(3)^2-2*E(3))*q-q^2,0,0,0,
E(3)^2-E(3)^2*q+E(3)^2*q^2,-1+2*q-q^2,0,0,E(3)^2-E(3)^2*q,0,0,2-q^-1-q,0,0,0,
0,0,0,-3+q^-1+3*q-q^2,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,-E(3)^2,0,0,0,0,0,0],[0,0,-1+q,0,0,0,-1+q,0,0,0,0,0,1-q,0,
0,0,0,0,0,1-q,0,0,0,E(3)^2-2*E(3)^2*q+E(3)^2*q^2,0,0,0,E(3)-E(3)*q,0,q,1-q,0,
0,0,0,0,0,0,0,0],[0,0,0,-q+q^2,0,0,1-q,0,0,0,-E(3)+E(3)*q-E(3)*q^2,0,E(3)+
(E(3)^2-E(3))*q-E(3)^2*q^2,0,0,E(3)-E(3)*q+E(3)*q^2,0,0,0,0,0,0,0,0,0,0,0,0,
E(3)-E(3)*q,0,0,0,0,0,0,-E(3)+(-E(3)^2+E(3))*q+E(3)^2*q^2,0,0,0,0],[0,0,
E(3)^2-E(3)+E(3)*q^-1-E(3)^2*q,0,0,-2+q^-1+q,2*E(3)^2-E(3)^2*q^-1-E(3)^2*q,0,
0,0,0,3-q^-1-3*q+q^2,-E(3)^2+E(3)+E(3)^2*q^-1-E(3)*q,0,0,0,0,-2+q^-1+q,0,
2*E(3)-E(3)*q^-1-E(3)*q,0,0,0,2-q^-1-2*q+q^2,0,0,0,-1+q^-1+q,0,
E(3)^2-E(3)^2*q,2*E(3)-E(3)*q^-1-E(3)*q,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,-1+q,0,
0,0,0,0,1-q,0,0,0,0,-E(3)+E(3)*q-E(3)*q^2,E(3)+(E(3)^2-E(3))*q-E(3)^2*q^2,0,0,
0,0,0,1-q+q^2,0,0,0,0,0,0,E(3)-E(3)*q,0,0,0,0,0,-E(3)+(-E(3)^2+E(3))*q+
E(3)^2*q^2,0,0,0],[-1+2*q-q^2,2*E(3)^2-E(3)-E(3)^2*q^-1+(-E(3)^2+
2*E(3))*q-E(3)*q^2,E(3)+(E(3)^2-2*E(3))*q+(-2*E(3)^2+E(3))*q^2+E(3)^2*q^3,
7*E(3)^2+E(3)^2*q^-2-4*E(3)^2*q^-1-7*E(3)^2*q+4*E(3)^2*q^2-E(3)^2*q^3,0,
-E(3)^2+2*E(3)-E(3)*q^-1+(2*E(3)^2-E(3))*q-E(3)^2*q^2,
3*E(3)^2-E(3)^2*q^-1-3*E(3)^2*q+E(3)^2*q^2,2*E(3)^2-E(3)-E(3)^2*q^-1+
(-2*E(3)^2+E(3))*q+E(3)^2*q^2,0,2-q^-1-2*q+q^2,-5*E(3)^2+2*E(3)-E(3)^2*q^-2+
(3*E(3)^2-E(3))*q^-1+(5*E(3)^2-E(3))*q-2*E(3)^2*q^2,3*E(3)^2-E(3)-E(3)^2*q^-1+
(-4*E(3)^2+2*E(3))*q+(3*E(3)^2-E(3))*q^2-E(3)^2*q^3,E(3)^2-2*E(3)+E(3)^2*q^-2+
(-2*E(3)^2+E(3))*q^-1+E(3)*q,-1+2*q-q^2,-E(3)^2+2*E(3)^2*q-2*E(3)^2*q^2+
E(3)^2*q^3,2*E(3)^2-2*E(3)+E(3)^2*q^-2+(-2*E(3)^2+E(3))*q^-1+(-2*E(3)^2+
E(3))*q+E(3)^2*q^2,0,-2*E(3)^2+E(3)+E(3)^2*q^-1+(E(3)^2-2*E(3))*q+E(3)*q^2,-1+
2*q-q^2,-2*E(3)^2+E(3)+E(3)^2*q^-1+(2*E(3)^2-E(3))*q-E(3)^2*q^2,0,-2+q^-1+
2*q-q^2,-3*E(3)^2+E(3)+E(3)^2*q^-1+(4*E(3)^2-2*E(3))*q+(-3*E(3)^2+E(3))*q^2+
E(3)^2*q^3,-E(3)^2-3*E(3)+E(3)*q^-1+(2*E(3)^2+4*E(3))*q+(-E(3)^2-3*E(3))*q^2+
E(3)*q^3,0,-2+q^-1+q,E(3)-E(3)*q+E(3)*q^2,-2+q^-1+2*q-q^2,2*E(3)^2-E(3)+
E(3)^2*q^-2+(-2*E(3)^2+E(3))*q^-1-E(3)^2*q,E(3)^2-2*E(3)^2*q+E(3)^2*q^2,
-2*E(3)^2+E(3)+E(3)^2*q^-1+(2*E(3)^2-E(3))*q-E(3)^2*q^2,E(3)-E(3)*q,1-q+q^2,0,
0,-5*E(3)^2+2*E(3)-E(3)^2*q^-2+(3*E(3)^2-E(3))*q^-1+
(6*E(3)^2-E(3))*q-4*E(3)^2*q^2+E(3)^2*q^3,0,2*E(3)^2-E(3)-E(3)^2*q^-1+
(-2*E(3)^2+2*E(3))*q+(2*E(3)^2-E(3))*q^2-E(3)^2*q^3,0,0],[0,-2+q^-1+q,1-2*q+
q^2,-4-q^-2+3*q^-1+3*q-q^2,0,2-q^-1-q,-2+q^-1+q,-2+q^-1+q,0,-2*E(3)+E(3)*q^-1+
E(3)*q,-3*E(3)^2-4*E(3)+q^-2-3*q^-1+(E(3)^2+2*E(3))*q,-3+q^-1+3*q-q^2,-1-q^-2+
2*q^-1,E(3)-E(3)*q,1-2*q+q^2,3*E(3)^2+2*E(3)-q^-2+(-3*E(3)^2-2*E(3))*q^-1+q,0,
2-q^-1-q,E(3)-E(3)*q,2-q^-1-q,0,2*E(3)-E(3)*q^-1-E(3)*q,3-q^-1-3*q+q^2,
2*E(3)^2-E(3)-E(3)^2*q^-1+(-2*E(3)^2+E(3))*q+E(3)^2*q^2,0,E(3)-E(3)*q^-1,1-q,
2*E(3)-E(3)*q^-1-E(3)*q,-1-q^-2+2*q^-1,-1+q,2-q^-1-q,1,E(3)^2-E(3)^2*q,1,0,4+
q^-2-3*q^-1-3*q+q^2,0,3*E(3)^2+2*E(3)+q^-1+(-3*E(3)^2-2*E(3))*q-q^2,0,0],[0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-E(3)*q,0,0,0,0,0,0,-1+q,0,0,
0,0,0,0],[0,0,0,0,0,0,0,0,E(3)*q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,-1+q,0,0,0,0,0],[0,0,0,-q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,q,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0],[-E(3)^2+E(3)^2*q,E(3)^2-E(3)^2*q,-E(3)^2+
E(3)^2*q,0,0,0,0,E(3)^2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,E(3)^2-E(3)^2*q,0,0],[0,0,0,-1+q,0,-1+q,0,0,0,0,0,0,0,E(3)*q,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,1-q,0,1-q,0,E(3)-E(3)*q,E(3)*q],[3-q^-1-3*q+q^2,
-3+q^-1+3*q-q^2,3-q^-1-3*q+q^2,-E(3)^2+E(3)^2*q,0,-E(3)^2+E(3)^2*q^-1,0,-2+
q^-1+q,-1+q,E(3)^2-E(3)^2*q^-1-E(3)^2*q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,2-q^-1-q,E(3)^2+2*E(3)-E(3)*q^-1+q,-E(3)^2-2*E(3)+E(3)*q^-1-q,-2+
q^-1+2*q-q^2,-E(3)^2+E(3)^2*q^-1+E(3)^2*q,E(3)^2-E(3)^2*q]],[[-1,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-q,q,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,
-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0],[0,0,-q,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,q,0,0],[1-q,0,1-q,0,1-q^-1,0,0,0,-1+q,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1-q,0,0,1-q,0,0,0,
-1+q,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,
0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,
0,0,0,0,-1,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,1,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0],[0,0,0,1-q,0,0,1-q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,
0,1,0,0,0,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,q,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,
E(3)^2*q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,q,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2*q,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
0,0,0,0],[0,0,0,-2+q^-1+q,0,0,-1+q^-1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,1-q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3),0,0,
0,0,0,0,0,0,0,0,E(3),0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,E(3),0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-2+q^-1+q,2-q^-1-q,
-1+2*q-q^2,3-q^-1-3*q+q^2,-2+q^-1+q,0,2-q^-1-q,2-q^-1-q,-E(3)+E(3)*q,
2*E(3)-E(3)*q^-1-E(3)*q,E(3)-E(3)*q,1-2*q+q^2,0,0,-1+2*q-q^2,0,3*E(3)^2+
2*E(3)+q^-1+(-3*E(3)^2-2*E(3))*q-q^2,0,-E(3)+E(3)*q,-2+q^-1+q,-2+q^-1+q,0,-1+
2*q-q^2,0,E(3)^2-E(3)^2*q,0,-1+q,0,0,0,0,0,0,-1,-2+q^-1+q,-1+2*q-q^2,-2+q^-1+
q,-3*E(3)^2-2*E(3)-q^-1+(3*E(3)^2+2*E(3))*q+q^2,-2*E(3)+E(3)*q^-1+E(3)*q,
E(3)-E(3)*q],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
E(3)^2*q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,-q,0,0,0,0,q,0,0,0,0,0,0,0,
0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,E(3),0,-1+q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,-q,0,0,0,0,0,0,E(3)^2*q,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,-1+q,E(3)^2*q,0,0,0,0,0,0,
0],[2*E(3)-E(3)*q^-1-E(3)*q,-2*E(3)+E(3)*q^-1+E(3)*q,E(3)-2*E(3)*q+E(3)*q^2,
-3*E(3)+E(3)*q^-1+3*E(3)*q-E(3)*q^2,2*E(3)-E(3)*q^-1-E(3)*q,0,-2*E(3)+
E(3)*q^-1+E(3)*q,-2*E(3)+E(3)*q^-1+E(3)*q,E(3)^2-E(3)^2*q,-2*E(3)^2+
E(3)^2*q^-1+E(3)^2*q,-E(3)^2+E(3)^2*q,-E(3)+2*E(3)*q-E(3)*q^2,0,0,
E(3)-2*E(3)*q+E(3)*q^2,0,E(3)^2+3*E(3)-E(3)*q^-1+(-E(3)^2-3*E(3))*q+E(3)*q^2,
0,E(3)^2-E(3)^2*q,2*E(3)-E(3)*q^-1-E(3)*q,2*E(3)-E(3)*q^-1-E(3)*q,0,
E(3)-2*E(3)*q+E(3)*q^2,0,-1+q,0,E(3)-E(3)*q,0,0,0,0,E(3),0,E(3),
2*E(3)-E(3)*q^-1-E(3)*q,E(3)-2*E(3)*q+E(3)*q^2,2*E(3)-E(3)*q^-1-E(3)*q,
-E(3)^2-3*E(3)+E(3)*q^-1+(E(3)^2+3*E(3))*q-E(3)*q^2,
2*E(3)^2-E(3)^2*q^-1-E(3)^2*q,-E(3)^2+E(3)^2*q],[-E(3)+2*E(3)*q-E(3)*q^2,-1+
2*q-q^2,-E(3)+2*E(3)*q-E(3)*q^2,-E(3)^2+2*E(3)^2*q-E(3)^2*q^2,
2*E(3)^2-E(3)^2*q^-1-E(3)^2*q,2-q^-1-q,-E(3)^2+2*E(3)^2*q-E(3)^2*q^2,
E(3)-E(3)*q,-E(3)^2+2*E(3)^2*q-E(3)^2*q^2,-E(3)+E(3)*q,E(3)^2-2*E(3)^2*q+
E(3)^2*q^2,-1+2*q-q^2,0,0,-q+q^2,0,1-2*q+q^2,0,0,-E(3)+E(3)*q,1-q,0,
E(3)*q-E(3)*q^2,0,0,0,-q,0,0,0,0,0,0,0,-E(3)+E(3)*q,E(3)^2-2*E(3)^2*q+
E(3)^2*q^2,1-2*q+q^2,-1+2*q-q^2,E(3)-E(3)*q,0],[1-q,0,0,0,1-q,0,0,0,q,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-q,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0],[0,0,0,0,0,-1,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0],[0,0,0,0,0,
0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,
0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q],[0,
0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
-1+q]],[[0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0],[1,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,q,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,q,0,0,0,0],[0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,E(3)*q,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,-1,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,-1+q,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,
E(3)^2,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0],[0,0,0,0,0,0,1,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1+q,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,-1+q,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,E(3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2*q,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2*q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
E(3)^2,0,0,0,0,0,0,0,0,0,0,0,0],[-2+q^-1+q,2-q^-1-q,-1+2*q-q^2,3-q^-1-3*q+q^2,
-2+q^-1+q,0,2-q^-1-q,2-q^-1-q,-E(3)+E(3)*q,2*E(3)-E(3)*q^-1-E(3)*q,
E(3)-E(3)*q,1-2*q+q^2,0,0,-1+2*q-q^2,0,3*E(3)^2+2*E(3)+q^-1+
(-3*E(3)^2-2*E(3))*q-q^2,0,-E(3)+E(3)*q,-2+q^-1+q,-2+q^-1+q,0,-1+2*q-q^2,0,
E(3)^2-E(3)^2*q,0,-1+q,0,0,0,0,0,0,-1,-2+q^-1+q,-1+2*q-q^2,-2+q^-1+q,
-3*E(3)^2-2*E(3)-q^-1+(3*E(3)^2+2*E(3))*q+q^2,-2*E(3)+E(3)*q^-1+E(3)*q,
E(3)-E(3)*q],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0,0,0,0,0,q,0,0,0,0,0,
0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3),0,0,0,0,0,0,-1+
q,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
E(3)*q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1-q,0,0,1-q,0,0,0,-1+q,0,0,
0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,q,0,0,0,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2,0,0,0,0,0,0,0,0],[0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)*q,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)*q,0,-1+q,0,0,
0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
E(3)^2,0,-1+q,0,0,0,0,0,0,0],[-3*E(3)^2+E(3)^2*q^-1+(4*E(3)^2+E(3))*q+
(-3*E(3)^2-2*E(3))*q^2-q^3,3*E(3)^2-E(3)^2*q^-1+(-4*E(3)^2-E(3))*q+(3*E(3)^2+
2*E(3))*q^2+q^3,-4*E(3)^2-2*E(3)-q^-1+(6*E(3)^2+2*E(3))*q+
(-4*E(3)^2-E(3))*q^2+E(3)^2*q^3,E(3)^2-3*E(3)+E(3)*q^-1+(-2*E(3)^2+4*E(3))*q+
(E(3)^2-3*E(3))*q^2+E(3)*q^3,2*E(3)-E(3)*q^-1-2*E(3)*q+E(3)*q^2,2-q^-1-q,
-2*E(3)+E(3)*q^-1+2*E(3)*q-E(3)*q^2,3*E(3)^2-E(3)^2*q^-1-3*E(3)^2*q+
E(3)^2*q^2,0,3-q^-1-3*q+q^2,0,E(3)^2+(-2*E(3)^2+E(3))*q+(E(3)^2-2*E(3))*q^2+
E(3)*q^3,0,E(3)-E(3)*q,-E(3)^2+(2*E(3)^2-E(3))*q+(-E(3)^2+
2*E(3))*q^2-E(3)*q^3,0,3*E(3)-E(3)*q^-1-4*E(3)*q+3*E(3)*q^2-E(3)*q^3,0,-1+
(-2*E(3)^2-E(3))*q+E(3)^2*q^2,2*E(3)-E(3)*q^-1-2*E(3)*q+E(3)*q^2,
2*E(3)-E(3)*q^-1-2*E(3)*q+E(3)*q^2,0,-E(3)^2+(2*E(3)^2-E(3))*q+(-E(3)^2+
2*E(3))*q^2-E(3)*q^3,0,-1+q-q^2,0,0,0,0,0,0,0,0,E(3)-E(3)*q,-E(3)^2+
2*E(3)-E(3)*q^-1+(2*E(3)^2-E(3))*q-E(3)^2*q^2,E(3)-3*E(3)*q+
3*E(3)*q^2-E(3)*q^3,-E(3)^2+2*E(3)-E(3)*q^-1+(2*E(3)^2-E(3))*q-E(3)^2*q^2,
4*E(3)^2-E(3)^2*q^-1-6*E(3)^2*q+4*E(3)^2*q^2-E(3)^2*q^3,3*E(3)^2+
E(3)-E(3)^2*q^-1+(-3*E(3)^2-2*E(3))*q-q^2,-E(3)^2+(2*E(3)^2+E(3))*q+q^2],
[q-q^2,-q+q^2,q-q^2,0,0,0,0,-q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,q,0,0,-q+q^2,0,0],[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3),0],[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2*q,0,-1+q,0],[0,0,0,0,0,0,0,
0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q]],[[0,0,1,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,
0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[q,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0],[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0],[0,q,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0],[0,0,0,q,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[q-q^2,0,q-q^2,0,-1+q,0,0,0,-q+q^2,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0],[0,0,0,0,0,0,0,0,-1+q,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0],[0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,
0,1,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0],[0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,q],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2*q,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1+q,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,E(3)*q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1-q,0,0,1-q,0,0,0,-1+q,
0,0,0,0,0,0,0,0,-1+q,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,
0,0,0,0,0,0,0,0,q,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,1-q^-1,0,0,0,0,0,-1+q,0,0,0,0,1-q,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,1-q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3),0,0,0,0,0,0,0,-1+q,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2,0,0,0,
0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,E(3)^2*q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0],[0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3),0,0,0,-1+q,0,0,0,0,0,0,0,0,0,
0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1+q,0,0,0,0,
0,0,0,0,0],[-1+2*q-q^2,1-2*q+q^2,q-2*q^2+q^3,1-3*q+3*q^2-q^3,-1+2*q-q^2,0,
1-2*q+q^2,1-2*q+q^2,E(3)*q-E(3)*q^2,E(3)-2*E(3)*q+E(3)*q^2,-E(3)*q+E(3)*q^2,
-q+2*q^2-q^3,0,0,q-2*q^2+q^3,0,-1+(-3*E(3)^2-2*E(3))*q+(3*E(3)^2+2*E(3))*q^2+
q^3,0,E(3)*q-E(3)*q^2,-1+2*q-q^2,-1+2*q-q^2,0,q-2*q^2+q^3,0,-E(3)^2*q+
E(3)^2*q^2,0,q-q^2,0,0,0,0,q,0,q,-1+2*q-q^2,q-2*q^2+q^3,-1+2*q-q^2,1+
(3*E(3)^2+2*E(3))*q+(-3*E(3)^2-2*E(3))*q^2-q^3,-E(3)+2*E(3)*q-E(3)*q^2,
-E(3)*q+E(3)*q^2],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3),0,
0,0,0,0,q,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,-1,0,0,0,0,0,0],[-1+q,1-q,-1+q,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,1-q,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1+q,0,0,0],[0,0,0,0,0,0,0,0,
q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,
0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0],[0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
-1]],[[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0],[q,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1-q,0,0,1-q,0,0,0,-1+q,0,0,
0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,
0,0,0,0,0,0,0,E(3)^2*q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,
0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,
0,0],[0,0,0,0,0,1,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,E(3),0,0,0,0,-1+q,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,-1+q,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,E(3)^2*q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,
0,0,0,0,0,0,q,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0],[0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)*q,0,0,0,0,0,0,0,0,0,
0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,1,0],[-1+q,1-q,-1+q,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,1-q,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3),0,0,0,
0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)*q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,
0,0,q,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,q],[0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,E(3)^2*q,0,0,0,0,0,0,0],[0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,
0,0,0,0],[q-q^2,0,q-q^2,0,-1+q,0,0,-q,-q+q^2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,-1+q,0,0,0,0,0,q,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,E(3)^2,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0],[-2+q^-1+q,2-q^-1-q,-1+
2*q-q^2,3-q^-1-3*q+q^2,-2+q^-1+q,0,2-q^-1-q,2-q^-1-q,-E(3)+E(3)*q,
2*E(3)-E(3)*q^-1-E(3)*q,E(3)-E(3)*q,1-2*q+q^2,0,0,-1+2*q-q^2,0,3*E(3)^2+
2*E(3)+q^-1+(-3*E(3)^2-2*E(3))*q-q^2,0,-E(3)+E(3)*q,-2+q^-1+q,-2+q^-1+q,0,-1+
2*q-q^2,0,E(3)^2-E(3)^2*q,0,-1+q,0,0,0,0,-1+q,0,-1,-2+q^-1+q,-1+2*q-q^2,-2+
q^-1+q,-3*E(3)^2-2*E(3)-q^-1+(3*E(3)^2+2*E(3))*q+q^2,-2*E(3)+E(3)*q^-1+E(3)*q,
E(3)-E(3)*q],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
-1,0,0,0,0,0,0,0],[0,-1+2*q-q^2,-q+2*q^2-q^3,-3+q^-1+4*q-3*q^2+q^3,0,1-2*q+
q^2,-1+2*q-q^2,-1+2*q-q^2,0,-E(3)+2*E(3)*q-E(3)*q^2,3-q^-1+(3*E(3)^2+
4*E(3))*q+(-E(3)^2-2*E(3))*q^2,-1+3*q-3*q^2+q^3,-2+q^-1+q,-E(3)*q+E(3)*q^2,-q+
2*q^2-q^3,3*E(3)^2+2*E(3)+q^-1+(-3*E(3)^2-2*E(3))*q-q^2,0,1-2*q+q^2,-E(3)*q+
E(3)*q^2,1-2*q+q^2,0,E(3)-2*E(3)*q+E(3)*q^2,1-3*q+3*q^2-q^3,E(3)^2+(-2*E(3)^2+
E(3))*q+(2*E(3)^2-E(3))*q^2-E(3)^2*q^3,0,E(3)-E(3)*q,-q+q^2,E(3)-2*E(3)*q+
E(3)*q^2,-2+q^-1+q,q-q^2,1-2*q+q^2,-q,-E(3)^2*q+E(3)^2*q^2,0,0,3-q^-1-4*q+
3*q^2-q^3,0,-1+(-3*E(3)^2-2*E(3))*q+(3*E(3)^2+2*E(3))*q^2+q^3,0,0],[0,0,0,0,0,
0,0,0,0,0,-1+q,0,1-q,0,0,1-q,0,0,0,0,0,0,q,0,0,0,0,0,1,0,0,0,0,0,0,-1+q,0,0,0,
0],[0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+
q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0]]];end;
  f8:=function(q) return
q^0*[[[E(3)-E(3)*q,E(3)*q,0,-q+q^2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0],[-E(3)^2+E(3)^2*q^-1+E(3)^2*q,E(3)^2-E(3)^2*q,0,-E(3)+((-3+
ER(-3))/2)*q+q^2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,q,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,q,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,E(3)^2*q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0],[0,0,0,0,0,E(3),-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0],[-2+q^-1+q,1-q,-E(3)^2+E(3)^2*q^-1,0,0,0,0,E(3)^2-E(3)^2*q,0,1,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,E(3)^2-E(3)^2*q,(5+ER(-3))/2-q^-1+
(-2-ER(-3))*q-E(3)^2*q^2,0,0,0,E(3)^2-E(3)^2*q,0,E(3)^2,0,0,0,0,0,0,0,0,0,0,0,
-E(3)^2+E(3)^2*q,0,0,0,0,0,0,0],[2*E(3)^2-E(3)^2*q^-1-2*E(3)^2*q+E(3)^2*q^2,
-E(3)^2+2*E(3)^2*q-E(3)^2*q^2,(-3+ER(-3))/2+q^-1-E(3)*q,-E(3)^2+
2*E(3)^2*q-E(3)^2*q^2,0,0,0,1-q+q^2,0,E(3)-E(3)*q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0],[0,0,0,E(3)-ER(-3)*q-E(3)^2*q^2,3*E(3)^2-E(3)^2*q^-1-4*E(3)^2*q+
3*E(3)^2*q^2-E(3)^2*q^3,0,0,0,E(3)-E(3)*q+E(3)*q^2,0,E(3)-E(3)*q,0,0,0,0,0,0,
0,0,E(3)^2-E(3)^2*q,0,0,E(3)^2-2*E(3)^2*q+E(3)^2*q^2,0,0,0,0,0,0,0],[0,0,
E(3)^2-E(3)^2*q^-1,0,1-2*q+q^2,0,0,0,0,0,0,E(3)^2-E(3)^2*q,0,0,E(3)^2,0,0,0,0,
0,0,0,1-q,0,0,0,0,0,0,0],[0,0,-2+ER(-3)-E(3)*q^-1+((5-ER(-3))/2)*q-q^2,0,1-q,
0,0,0,0,0,0,0,E(3)^2-E(3)^2*q,0,0,0,0,E(3)*q,0,0,0,0,0,0,0,0,0,1-q,0,0],[0,0,
-2+ER(-3)-E(3)*q^-1+((5-ER(-3))/2)*q-q^2,-E(3)^2+E(3)^2*q,0,0,0,0,0,0,0,0,0,
E(3)^2-E(3)^2*q,0,E(3)^2*q,0,0,0,0,0,0,0,0,0,0,0,1-q,0,0],[0,0,
ER(-3)-E(3)*q^-1+E(3)^2*q,0,E(3)^2+(2+ER(-3))*q+((-5-ER(-3))/2)*q^2+q^3,0,0,0,
0,0,0,E(3)-E(3)*q+E(3)*q^2,0,0,E(3)-E(3)*q,0,0,0,0,-1+q,0,0,-1+2*q-q^2,0,0,0,
0,0,0,0],[0,0,-4-q^-2+3*q^-1+3*q-q^2,ER(-3)-E(3)*q^-1+E(3)^2*q,0,0,0,0,0,0,0,
0,0,-E(3)+E(3)*q^-1+E(3)*q,0,E(3)-E(3)*q,0,0,0,0,0,0,0,0,0,0,-E(3)^2+E(3)^2*q,
2-q^-1-q,0,0],[3-q^-1-3*q+q^2,-1+2*q-q^2,0,0,-1+3*q-3*q^2+q^3,1-2*q+q^2,
E(3)^2*q-2*E(3)^2*q^2+E(3)^2*q^3,E(3)-ER(-3)*q-E(3)^2*q^2,0,-1+q,0,
E(3)-ER(-3)*q-E(3)^2*q^2,0,0,-E(3)^2+E(3)^2*q,0,-1,0,0,0,0,0,-1+2*q-q^2,0,1-q,
E(3)^2*q-E(3)^2*q^2,0,0,0,0],[0,0,-4*E(3)-E(3)*q^-2+3*E(3)*q^-1+
3*E(3)*q-E(3)*q^2,0,(-3+ER(-3))/2+q^-1-E(3)*q,0,0,0,0,0,0,0,-E(3)^2+
E(3)^2*q^-1+E(3)^2*q,0,0,0,0,E(3)-E(3)*q,0,0,0,0,0,0,0,0,-1+q,
2*E(3)-E(3)*q^-1-E(3)*q,0,0],[-2+q^-1+q,0,0,-3+q^-1+3*q-q^2,1-3*ER(-3)+
ER(-3)*q^-1+((-5+7*ER(-3))/2)*q-4*E(3)*q^2+E(3)*q^3,-1+
2*ER(-3)-ER(-3)*q^-1-3*E(3)*q+E(3)*q^2,(-3+ER(-3))/2+((7-ER(-3))/2)*q-3*q^2+
q^3,4*E(3)^2-E(3)^2*q^-1+(2+3*ER(-3))*q-2*ER(-3)*q^2+E(3)*q^3,0,(5+
ER(-3))/2-q^-1+(-2-ER(-3))*q-E(3)^2*q^2,0,-2-ER(-3)-E(3)^2*q^-1+((5+
ER(-3))/2)*q-q^2,(-3+ER(-3))/2+q^-1-E(3)*q,3-q^-1-3*q+q^2,-2+q^-1+q,-1+
2*q-q^2,0,-1+q,-1,0,-E(3)^2-ER(-3)*q+E(3)*q^2,0,-3*E(3)+E(3)*q^-1+
3*E(3)*q-E(3)*q^2,0,2*E(3)-E(3)*q^-1-E(3)*q,-1+2*q-q^2,1-2*q+q^2,
2*E(3)^2-E(3)^2*q^-1-E(3)^2*q,0,-E(3)+E(3)*q],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,E(3)*q,0,0,0,0,0,0,0],[-4+2*q^-1+3*q-q^2,1-2*q+q^2,0,-2+q^-1+q,
0,0,0,-1-2*ER(-3)-E(3)^2*q^-1+((1+5*ER(-3))/2)*q-ER(-3)*q^2,0,(5+
ER(-3))/2-q^-1+((-3-ER(-3))/2)*q,0,0,0,2-q^-1-q,0,-1+q,0,0,0,0,E(3)-E(3)*q,0,
0,0,0,0,1-q,0,0,-E(3)],[0,0,0,0,2-ER(-3)+E(3)*q^-1+((-7+ER(-3))/2)*q+
3*q^2-q^3,0,0,0,0,0,0,-1-2*ER(-3)-E(3)^2*q^-1+((1+5*ER(-3))/2)*q-ER(-3)*q^2,
-2+q^-1+q,0,2*E(3)^2-E(3)^2*q^-1-E(3)^2*q,0,0,-1+q,0,-1+q^-1,0,E(3)-E(3)*q,
1-2*q+q^2,0,0,0,-1+q,0,-E(3)^2,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
E(3)^2,0,0,-1+q,0,0,0,0,0,0,0],[5+2*ER(-3)+((-3-ER(-3))/2)*q^-1+
(-6-3*ER(-3))*q+(3+2*ER(-3))*q^2+E(3)^2*q^3,E(3)^2-3*E(3)^2*q+
3*E(3)^2*q^2-E(3)^2*q^3,0,4-q^-1-6*q+4*q^2-q^3,(-3+7*ER(-3))/2-ER(-3)*q^-1+
((11-9*ER(-3))/2)*q+((-15+5*ER(-3))/2)*q^2+((9-ER(-3))/2)*q^3-q^4,
(3-5*ER(-3))/2+ER(-3)*q^-1+(-4+2*ER(-3))*q+((7-ER(-3))/2)*q^2-q^3,
(3-ER(-3))/2+((-9+ER(-3))/2)*q+(5+ER(-3))*q^2+
((-5-3*ER(-3))/2)*q^3-E(3)^2*q^4,(3+5*ER(-3))/2+E(3)^2*q^-1+
((-3-11*ER(-3))/2)*q+6*ER(-3)*q^2+(1-3*ER(-3))*q^3+E(3)*q^4,1+((-3+
ER(-3))/2)*q-E(3)*q^2,-4-ER(-3)+q^-1+((11+5*ER(-3))/2)*q+
(-3-2*ER(-3))*q^2-E(3)^2*q^3,-E(3)+E(3)*q,(3+5*ER(-3))/2+E(3)^2*q^-1+
((-3-9*ER(-3))/2)*q+((1+7*ER(-3))/2)*q^2-ER(-3)*q^3,3-q^-1-3*q+q^2,-3+q^-1+
((9+ER(-3))/2)*q+((-7-ER(-3))/2)*q^2+q^3,-3*E(3)^2+E(3)^2*q^-1+
3*E(3)^2*q-E(3)^2*q^2,1+((-5-ER(-3))/2)*q+((5+ER(-3))/2)*q^2-q^3,0,1-2*q+q^2,
0,2-q^-1-q,E(3)^2+((1+3*ER(-3))/2)*q+((1-3*ER(-3))/2)*q^2+E(3)*q^3,E(3)^2+
ER(-3)*q-E(3)*q^2,(-3+ER(-3))/2+(4-ER(-3))*q+((-7+ER(-3))/2)*q^2+q^3,-1,
2-ER(-3)+E(3)*q^-1+((-5+ER(-3))/2)*q+q^2,1+((-5-ER(-3))/2)*q+(2+ER(-3))*q^2+
E(3)^2*q^3,q-2*q^2+q^3,-E(3)^2+2*E(3)^2*q-E(3)^2*q^2,E(3)^2-E(3)^2*q,
E(3)-2*E(3)*q+E(3)*q^2],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
E(3)^2*q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3),-1+q,
0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)*q,-1+q,0,0],[0,0,
-4-q^-2+3*q^-1+4*q-3*q^2+q^3,0,-E(3)+2*E(3)*q-E(3)*q^2,0,0,0,0,0,0,
3*E(3)-E(3)*q^-1-4*E(3)*q+3*E(3)*q^2-E(3)*q^3,2-q^-1-2*q+q^2,0,
2*E(3)-E(3)*q^-1-E(3)*q,0,0,1-2*q+q^2,0,0,0,-E(3)+E(3)*q-E(3)*q^2,
E(3)^2-E(3)^2*q,0,0,0,-E(3)^2+2*E(3)^2*q-E(3)^2*q^2,2-q^-1-2*q+q^2,
E(3)^2-E(3)^2*q,0],[-5+2*q^-1+5*q-3*q^2+q^3,2-4*q+3*q^2-q^3,1-2*q+q^2,
-ER(-3)-E(3)^2*q^-1+((-5+ER(-3))/2)*q+3*q^2-q^3,0,0,0,-1-2*ER(-3)-E(3)^2*q^-1+
(1+3*ER(-3))*q+((-1-5*ER(-3))/2)*q^2+ER(-3)*q^3,0,(7+ER(-3))/2-q^-1+
(-4-ER(-3))*q+((3+ER(-3))/2)*q^2,0,0,0,2-q^-1-2*q+q^2,0,-1+2*q-q^2,0,0,0,0,
-E(3)^2+E(3)^2*q-E(3)^2*q^2,0,0,0,0,0,0,E(3)-E(3)*q,0,E(3)^2-E(3)^2*q]],[[0,q,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1,-1+q,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0],[0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
[1-2*q+q^2,q-q^2,0,0,-E(3)+E(3)*q,-1+q,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0],[-1+q,0,0,0,E(3)-E(3)*q,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0],[0,0,0,0,0,0,0,0,0,E(3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,
0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,
E(3)^2*q,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,q,0,
-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2*q,0,0,
0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,
0,0],[0,0,0,0,0,0,0,0,0,0,0,q,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,
0,0,0,0,0,0,0,0,0,0,1,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,E(3),0,0,
0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,
0,0,0,0,0,0,0,0,0],[0,0,0,0,q-q^2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,q,0,0,
0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1],[0,
0,0,0,0,0,0,0,0,0,0,0,0,0,1-q,0,0,1-q,0,0,0,0,0,0,0,0,0,0,-E(3),0],[0,0,0,0,
-q+q^2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0],[-2+ER(-3)-E(3)*q^-1+(4-ER(-3))*q+
((-7+ER(-3))/2)*q^2+q^3,-E(3)+(-2+ER(-3))*q+((5-ER(-3))/2)*q^2-q^3,0,-E(3)+
2*E(3)*q-E(3)*q^2,1+2*ER(-3)+E(3)^2*q^-1+((-1-5*ER(-3))/2)*q+ER(-3)*q^2,0,
-E(3)+(-2+ER(-3))*q+((3-ER(-3))/2)*q^2,0,-E(3)+2*E(3)*q-E(3)*q^2,0,-E(3)+
E(3)*q,0,0,0,0,0,0,0,0,-E(3)+E(3)*q,0,0,0,0,-1+q,E(3)*q,0,0,0,0],[1-2*q+q^2,
1-q,0,1-2*q+q^2,1-2*q+q^2,(-5-ER(-3))/2+q^-1+((3+ER(-3))/2)*q,0,0,1-q,0,0,0,0,
0,0,0,0,0,0,2-q^-1-q,0,0,1-q,0,E(3)^2,0,0,0,0,0],[0,0,-E(3)+E(3)*q,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,E(3),0,0],[0,0,q-q^2,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)^2*q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,
E(3)^2*q-E(3)^2*q^2,1-q,0,0,0,0,0,0,0,0,-E(3)^2*q,0,0,0,0,0,0,-1+q,0],[0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0,0,0,0,0,0,0,-1+q]],[[0,0,0,-q,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-q,q,0,1-q,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0],[-1,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0],[0,0,0,0,0,E(3)^2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,
0,0,0,E(3)*q,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,
0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,-1+q,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,-E(3)+E(3)*q,0,0,0,-1+q,0,0,0,
0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0],[1-q,0,0,1-2*q+q^2,E(3)-2*E(3)*q+E(3)*q^2,
(-5-ER(-3))/2+q^-1+(2+ER(-3))*q+E(3)^2*q^2,-q+q^2,0,1-2*q+q^2,0,0,0,0,0,0,0,0,
0,0,1-q,0,0,0,0,E(3)^2-E(3)^2*q,-q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0],[0,0,-E(3)^2+E(3)^2*q,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,
0,0,0,0,q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,
0],[0,0,0,0,0,0,0,0,0,(3+ER(-3))/2-q^-1+E(3)^2*q,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,
0,0,0,0,0,0,-E(3)],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,-1+q,0,0,0,0,0,0,0,0,0,0,
0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1-q^-1,0,-E(3)^2,0,0,0,0,0,0,0,0,0,0,
0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,E(3)-E(3)*q,0,0,-E(3)*q,-1+q,0,0,0,0,0,0,0,0,
0,0,0],[0,0,0,0,0,-1+q,0,0,-q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,E(3)^2,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q+q^2,0,0,
-q+q^2,0,0,0,0,0,-1+q,0,0,0,0,E(3)*q,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,E(3)*q,0,-1+q,0,0,0,0,0],[0,0,0,-1+q,2-ER(-3)+E(3)*q^-1+((-5+
ER(-3))/2)*q+q^2,0,-1+q,0,-1+q,0,-1,0,0,0,0,0,0,0,0,0,0,0,1-q,0,0,-1+q,0,0,0,
0],[0,0,0,0,0,0,0,E(3)^2-E(3)^2*q,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,E(3)^2-E(3)^2*q^-1,0,-E(3)+E(3)*q,0,0,0,0,E(3)^2,0,0,0,0,
0,0],[0,0,0,0,0,0,0,0,0,-E(3)^2-ER(-3)*q+E(3)*q^2,0,0,0,0,0,-E(3)^2*q,0,0,0,0,
0,0,0,0,0,0,0,0,0,0]],[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
0,0,0],[0,0,0,0,1-q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1-q^-1,0,0,1,0,0,0,0,0,0,0],
[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,E(3),
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,E(3)^2*q,-1+q,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1-q,-2+ER(-3)-E(3)*q^-1+
((5-ER(-3))/2)*q-q^2,0,0,0,1-q,0,1,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,
0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,
0,E(3)^2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-2*q+q^2,q-q^2,0,0,-E(3)+E(3)*q,-1+
q,q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,q,0,0,0,
-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1+q,-q,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,E(3)*q,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-E(3)^2,0,0,0,0,0,0,0,0,0,0,0,0],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,-E(3)*q,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,1-q,0,-1+q,0,0,0,0,-1,0,0,0,0,0,0],[q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,-1+q,0,0,0,0,0,
0,0,0],[1-q,q,0,-E(3)^2*q+E(3)^2*q^2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,
0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q-q^2,0,-q,0,0,0,0,0,0,0,0,0,
0,0],[0,0,0,0,0,q-q^2,0,0,-q+q^2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0],
[0,0,0,-1+2*q-q^2,(5-3*ER(-3))/2+E(3)*q^-1+((-9+3*ER(-3))/2)*q+
((7-ER(-3))/2)*q^2-q^3,0,q-q^2,0,-1+2*q-q^2,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,
1-2*q+q^2,0,0,q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
-1,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0],[0,0,
0,0,0,0,0,0,0,q-q^2,0,0,0,0,0,-q+q^2,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,E(3)^2*q],
[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,E(3),0]],[[-1,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1+q,0,E(3)^2*q,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0],[0,0,E(3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0],[0,0,0,0,0,0,0,E(3)^2*q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0],[-2+q^-1+q,1-q,-E(3)^2+E(3)^2*q^-1,0,0,0,0,E(3)^2-E(3)^2*q,0,1,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,E(3),0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0],[1-2*q+q^2,q-q^2,0,0,-E(3)+E(3)*q,-1+q,q,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,
0],[0,0,0,0,0,0,0,0,0,0,0,0,-E(3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,
0,0,0,0,0,0,0,-E(3)^2*q,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,
0,0,q,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,-q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1+q,0,0,
0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,-E(3)^2*q,0,
0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1+q,0,0,0,0,0,0,0,0,
0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-E(3),0,0,0,0,0,0,0,0,0,0,0,0,0],[0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0,0],[0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,-E(3),0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-q+
q^2,-E(3)+E(3)*q,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-E(3),0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+
q,0,-E(3)^2*q+E(3)^2*q^2,0,0,0,0,q,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,-E(3)^2*q,0,0,0,0,0,0,0,0,0],[4-2*q^-1-3*q+q^2,-1+2*q-q^2,0,
2-q^-1-q,0,0,0,1+2*ER(-3)+E(3)^2*q^-1+((-1-5*ER(-3))/2)*q+ER(-3)*q^2,0,
(-5-ER(-3))/2+q^-1+((3+ER(-3))/2)*q,0,0,0,-2+q^-1+q,0,1-q,0,0,0,0,-E(3)+
E(3)*q,0,0,0,0,0,-1+q,0,0,E(3)],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,
0,0,0,0,-1+q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-E(3)^2*q,0,
0,0,0,-1+q,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-E(3)^2+E(3)^2*q,0,0,-E(3)^2+
E(3)^2*q,0,0,0,0,0,0,0,0,0,0,q,0],[(-1-3*ER(-3))/2-E(3)^2*q^-1+((-1+
5*ER(-3))/2)*q+(1-2*ER(-3))*q^2+E(3)*q^3,-E(3)^2+((-1-3*ER(-3))/2)*q+((-1+
3*ER(-3))/2)*q^2-E(3)*q^3,0,-E(3)^2+2*E(3)^2*q-E(3)^2*q^2,(-7-ER(-3))/2+q^-1+
(4+ER(-3))*q+((-3-ER(-3))/2)*q^2,0,-E(3)^2+((-1-3*ER(-3))/2)*q+ER(-3)*q^2,0,
-E(3)^2+2*E(3)^2*q-E(3)^2*q^2,0,-E(3)^2+E(3)^2*q,0,0,0,0,0,0,0,0,-E(3)^2+
E(3)^2*q,0,0,0,0,-E(3)+E(3)*q,E(3)^2*q,0,0,0,-1+q]]];end;
    r:=[f1(x),f2(x),f3(x,E(3)),f4(x,E(3)),f4(x,E(3)^2),f3(x,E(3)^2),
     [[[-1]],[[-1]],[[-1]],[[-1]],[[-1]]],f8(x),f9(x),f8(1/x)*-x,f11(x,E(3)),
     f12(x,E(3)),f13(x,E(3)),f13(x,E(3)^2),f12(x,E(3)^2),f11(x,E(3)^2),f17(x),
     f1(1/x)*-x,f13(1/x,E(3)^2)*-x,f20(x,E(3)),f20(x,E(3)^2),f13(1/x,E(3))*-x,
     f23(x),f2(1/x)*-x,f11(1/x,E(3)^2)*-x,f12(1/x,E(3))*-x,f12(1/x,E(3)^2)*-x,
     f11(1/x,E(3))*-x,f29(x),f4(1/x,E(3)^2)*-x,f4(1/x,E(3))*-x,
     f23(1/x)*-x,f3(1/x,E(3)^2)*-x,f3(1/x,E(3))*-x,
     f17(1/x)*-x,[[[x]],[[x]],[[x]],[[x]],[[x]]]];
     return r[i];
     end;
     return m(i);
  elif [p,q,r]=[4,4,3] then 
     x:=-para[2][1]/para[2][2];
     r:=x^0*
 [[[[x,-1,-1,0,0,0],[0,-1,0,0,0,0],[0,0,-1,0,0,0],[0,0,0,-1,0,0],[0,0,-1+x,-1,
x,0],[0,1,-1,-1,0,x]],[[-1,0,0,0,0,0],[-x,x,0,0,0,x],[0,0,0,0,-x,0],[0,0,0,x,
0,0],[0,0,-1,0,-1+x,0],[0,0,0,0,0,-1]],[[x,-1,0,0,0,-1],[0,-1,0,0,0,0],[0,0,x,
0,1,-1],[0,-x,0,x,-1,1-x],[0,0,0,0,-1,0],[0,0,0,0,0,-1]]],[[[-1,0,0],[0,0,1],
[0,x,-1+x]],[[x,0,0],[-1,-1,0],[1,0,-1]],[[0,1,0],[x,-1+x,0],[0,0,-1]]],[[[x,
0,0],[x,-1,0],[x,0,-1]],[[-1,2,0],[0,x,0],[0,(-E(4)+1)*x,-1]],[[-1,0,1],[0,-1,
(E(4)+1)/2],[0,0,x]]],[[[x,0,0],[x,-1,0],[x,0,-1]],[[-1,2,0],[0,x,0],[0,
(E(4)+1)*x,-1]],[[-1,0,1],[0,-1,(-E(4)+1)/2],[0,0,x]]],[[[-1]],[[-1]],[[-1]]],
[[[x,-1,0],[0,-1,0],[0,-1,x]],[[-1+x,0,1],[0,x,0],[x,0,0]],[[0,x,0],[1,-1+x,
0],[0,0,x]]],[[[-1,0,0],[-1,x,0],[-1,0,x]],[[x,-2*x,0],[0,-1,0],[0,
-E(4)-1,x]],[[x,0,-x],[0,x,((E(4)-1)/2)*x],[0,0,-1]]],[[[-1,0,0],[-1,x,0],[-1,
0,x]],[[x,-2*x,0],[0,-1,0],[0,E(4)-1,x]],[[x,0,-x],[0,x,((-E(4)-1)/2)*x],[0,0,
-1]]],[[[-1,0],[-1,x]],[[-1,0],[-1,x]],[[x,-x],[0,-1]]],[[[x]],[[x]],[[x]]]];
    return r[i];
  else
    S:=CHEVIE.imp.CharInfo(p,q,r).charparams[i];
    # S is a p-tuple of partitions of area r
    p1rRep:=function()local Q,pos,ct;
    # Model of Ariki,  Halverson-Ram for reps of G(p,1,r): 
    # needs rational fractions if the parameters are indeterminates
      if r>1 then Q:=-para[2][1]/para[2][2];else Q:=0;fi;
      pos:=function(t,i)local j,k,l;# return [j,k,p] if i at [k,p] in t[j]
	for j in [1..Length(t)] do for k in [1..Length(t[j])] do
	  l:=Position(t[j][k],i); if l<>false then return [j,k,l];fi;
	od; od;
      end;
      ct:=p->para[1][p[1]]*Q^(p[3]-p[2]);
      T:=Tableaux(S);
      return Concatenation([DiagonalMat(List(T,S->ct(pos(S,1))))],
	List([2..r],i->List([1..Length(T)],function(j)local S,v,a,b,p,tll;
	  S:=T[j];  a:=pos(S,i); b:=pos(S,i-1);
	  S:=List(S, a-> List(a, ShallowCopy)); # no `Copy' in GAP 4
	  S[a[1]][a[2]][a[3]]:=i-1; S[b[1]][b[2]][b[3]]:=i; 
	  if para[2][1]=-para[2][2] then
	    if a[1]=b[1] then tll:=para[2][1]/(a[3]+b[2]-a[2]-b[3]);
	    else tll:=0; 
	    fi;
	  else tll:=Sum(para[2])/(1-ct(b)/ct(a));
	  fi;
	  v:=[1..Length(T)]*0; v[j]:=tll; 
	  p:=Position(T,S);if p<>false then v[p]:=tll-para[2][2];fi;
	  return v;end)))*Product(para,Product)^0;
    end;
    if q=1 then return p1rRep();
    elif p=q then para:=[List([0..p-1],i->E(p)^i),para[1]];
    else
      if para[2]<>para[3] then
        if q mod 2=0 and r=2 then
          S:=CHEVIE.imp.CharInfo(p,q,r).malle[i];
	  if S[1]=1 then 
	    return [[[para[1][1+((S[4]-1) mod (p/q))]]],
	     [[para[2][S[2]]]],[[para[3][S[3]]]]];
	  else Y:=para[2];T:=para[3];
	    if q>2 then X:=List(para[1],y->GetRoot(y,q/2));
	      X:=Concatenation(List([1..q/2],i->E(q/2)^i*X));
	    else X:=para[1];fi;
	    X:=X{S{[3,4]}};
	    v:=S[2]*GetRoot(Product(X)*Product(Y)*Product(T)*
	      E(p/q)^(2-S[3]-S[4]),2)*E(p)^(S[3]+S[4]-2);
            d:=1+Sum(X)*0+Sum(Y)*0+Sum(T)*0;
	    return [(d*[[X[1],Sum(Y,y->1/y)-X[2]/v*Sum(T)],[0,X[2]]])^(q/2),
	      [[Sum(Y),1/X[1]],[-Product(Y)*X[1],0]],
	      [[0,-Product(T)/v],[v,Sum(T)]]];
	  fi;
	else Error("should not happen");
	fi;
      elif para[1]=List([1..p/q],i->E(p/q)^(i-1)) then
	  para:=[List([0..p-1],i->E(p)^i),para[2]];
      else para:=[Concatenation(TransposedMat(List(para[1],i->List([0..q-1],
	j->E(q)^j)*GetRoot(i,q)))),para[2]];
      fi;
    fi;
    extra:=false;
    if IsInt(S[Length(S)]) then extra:=E(S[Length(S)-1])^S[Length(S)]; 
      d:=Length(S)-2;
      S:=FullSymbol(S);
    fi;
    v:=p1rRep();
    if p=q then v:=Concatenation([v[2]^v[1]],v{[2..Length(v)]});
    elif q>1 then v:=Concatenation([v[1]^q,v[2]^v[1]],v{[2..Length(v)]});
    fi;
    if extra<>false then
      m:=PermListList(T,List(T,S->S{Concatenation([d+1..p],[1..d])}));
      m:=Cycles(m,[1..Length(T)]);
      l:=List([0,-1..1-p/d],i->extra^i);m1:=List(m,x->x[1]);
      return List(v,x->List(m,c->l*x{c}{m1}));
    else return v;
    fi;
  fi;
end);

CHEVIE.AddData("Representation","imp",function(p,q,r,i)local o;
  o:=CHEVIE.R("EigenvaluesGeneratingReflections","imp")(p,q,r);
  o:=List(o,Denominator);
  return CHEVIE.R("HeckeRepresentation","imp")
         (p,q,r,List(o,x->List([0..x-1],i->E(x)^i)),[],i);
end);

CHEVIE.AddData("CharTable","imp",function(p,q,r)local o;
  o:=CHEVIE.R("EigenvaluesGeneratingReflections","imp")(p,q,r);
  o:=List(o,Denominator);
  return CHEVIE.R("HeckeCharTable","imp")
	 (p,q,r,List(o,x->List([0..x-1],i->E(x)^i)),[]);
end);

##############################################################################
##  FamilyImprimitive(S) . . . . . . . . .  returns symbols, Fourier matrix, 
##    Frobenius eigenvalues for the family of the symbol S, following
##  G. Malle, "Unipotente Grade...", J. Algebra 177 (1995), Sect. 4  and 6 
##  for G(e,1,n) or G(e,e,n).
##  Initial writing  GM 26.10.2000 modified JM 10.08.2011
##############################################################################
FamilyImprimitive:=function(S)# fourier for family of symbol S
  local e,Scoll,ct,d,m,ll,eps,equiv,nrSymbols,epsreps,trace,roots,i,j,mat,
   frobs,symbs,newsigns,schon,orb,mult,res,IsReducedSymbol;
  e:=Length(S);Scoll:=Collected(Flat(S));
  ct:=Concatenation(List(Scoll,x->[1..x[2]]*0+x[1])); 
# Fourier matrix of the family of symbols with content ct:
# Let F be the set of functions ct->[0..e-1] which are injective restricted
# to  a given value in  ct, and satisfy a  certain condition for the sum of
# their values (see below). Then for f\in F the list of preimages of f is a
# symbol  S(f). Conversely for a  symbol S there is  a 'canonical' map f(S)
# which is increasing on entries of ct of given value. Then the formula is
#
# mat[S][T]=\sum_{f\mid S(f)=S}\eps(f)\eps(f(T))\zeta_e^{f*f(T)}
#
# where for f in F with ordered image im(f):=List(ct,f)
# where \eps(f)=(-1)^{number of non-inversions in the list im(f)}
# and f*f(T) is the scalar product of vectors im(f) and im(f(T))
# 
# To  compute this reasonably  fast, it can  be decomposed as  a product of
# sums, each relative to a set of consecutive equal entries in ct.
  d:=Length(ct) mod e;
  if not(d in [0,1]) then
    Error("Length(",IntListToString(ct),") should be 0 or 1 mod ",e," !\n");
  fi;
  m:=(Length(ct)-d)/e; 
  j:=(m*Binomial(e,2)) mod e; # for f in F we must have Sum(ct,f) mod e = j
  ll:=Cartesian(List(Scoll,i->[0..e-1]));
  ll:=Filtered(ll,x->Sum(x)mod e=j);
  ll:=List(ll,c->Zip(Scoll,c,function(x,y)
    return Filtered(Combinations([0..e-1],x[2]),c->Sum(c)mod e=y);end));
  nrSymbols:=Sum(ll,x->Product(x,Length));
# trace:=nrSymbols>10;
# if trace then InfoChevie("e=",e," ct=",IntListToString(ct),
# " nrSymbols=",nrSymbols,"\c");fi;
  ll:=Concatenation(List(ll,Cartesian)); # equiv classes of F
  eps:=l->(-1)^Sum([1..Length(l)],i->Number(l{[i+1..Length(l)]},j->l[i]<j));
  equiv:=List(ll,x->rec(globaleps:=(-1)^Sum([1..Length(x)],i->
    Sum([i+1..Length(x)],j->Sum(x[i],y->Number(x[j],k->y<k)))),
   aa:=List(x,function(y)local a;
    a:=Arrangements(y,Length(y));
    return List(a,x->rec(l:=x,eps:=eps(x)));end)));
  # a record in equiv describes f in F with a given f(S)
  # .aa is the list of possible im(f)
  epsreps:=List(ll,x->eps(Concatenation(x)));
# if trace then InfoChevie(" so far:",Stime(),
#   " next ",nrSymbols,"x",Join(List(equiv[1].aa,Length),"x"),":\c");
# fi;
  roots:=List([0..e-1],i->E(e)^i);
  mat:=List(equiv,i->i.globaleps*List([1..nrSymbols],k->epsreps[k]*
    Product([1..Length(i.aa)],
      l->Sum(i.aa[l],j->j.eps*roots[1+((-j.l*ll[k][l])mod e)]))));
  mat:=(-1)^(m*(e-1))*mat/(E(4)^Binomial(e-1,2)*ER(e)^e)^m;
# if trace then InfoChevie(Stime(),"\n");fi;
  frobs:=E(12)^(-(e^2-1)*m)*List(ll,i->E(2*e)^(-Sum(i,j->j*j)-e*Sum(i,Sum)));
  symbs:=List(ll,function(l)local sy,j;sy:=List([1..e],j->[]);
    Zip(Concatenation(l),ct,function(v,c)Add(sy[v+1],c);return 1;end);
    return sy;
  end);
  newsigns:=(-1)^(Binomial(e,2)*Binomial(m,2))*List(symbs,
    i->(-1)^([0..e-1]*List(i,x->Binomial(Length(x),2))));
  mat:=Zip(newsigns,mat,function(s,l)
    return s*Zip(newsigns,l,function(x,y)return x*y;end);end);
# if trace then InfoChevie("finally:",Stime(),"\n");fi;
  if d=0 then # compact entries...
    IsReducedSymbol:=s->ForAll(Rotations(s){[2..Length(s)]},
      x->s=x or LessSymbols(x,s));
    schon:=List(symbs,IsReducedSymbol);mult:=[];
    for i in [1..nrSymbols] do if schon[i] then
      orb:=Set(Rotations(symbs[i])); Add(mult,e/Length(orb)); # Symmetriegruppe
      for j in Filtered([i+1..nrSymbols],j->symbs[j] in orb)
      do schon[j]:=false;od;
    fi; od;
    frobs:=Concatenation(Zip(mult,ListBlist(frobs,schon),
      function(m,f)return [1..m]*0+f;end));
    symbs:=Concatenation(Zip(mult,ListBlist(symbs,schon),
      function(m,s)if m=1 then return [s];
        else return List([0..m-1],j->Concatenation(s{[1..e/m]},[m,j]));
	fi;end));
    mat:=Concatenation(Zip(mult,ListBlist(mat,schon),function(m,l) 
     return List([1..m],i->Concatenation(Zip(mult,ListBlist(l,schon),
       function(n,c)return [1..n]*0+e*c/m/n;end)));end));
    mult:=Concatenation(List(mult,m->[1..m]*0+m));
    nrSymbols:=Length(symbs);
    for i in [1..nrSymbols] do
      for j in [1..nrSymbols] do 
	if FullSymbol(symbs[i])=FullSymbol(symbs[j]) then 
	  mat[i][j]:=mat[i][j]-1/mult[i];
	  if symbs[i]=symbs[j] then mat[i][j]:=mat[i][j]+1; fi;
	fi; 
      od;
    od;
#   if trace then InfoChevie("reduction to ",nrSymbols," symbols:",Stime(),"\n");fi;
    if (mat*DiagonalMat(frobs))^3<>mat^0 then Print("** WARNING: (S*T)^3<>1\n");
    fi;
  fi;
  res:=rec(symbols:=symbs,fourierMat:=mat,eigenvalues:=frobs);
  res.name:=IntListToString(ct);
  res.explanation:="classical family";
# res.name:=String(Concatenation(List(Scoll,function(x)
#     if x[2]=1 then return String(x[1]);
#     else return SPrint(x[1],"^{",x[2],"}");fi;end)));
  res.special:=1;
  # the next should be improved
  res.charLabels:=List([1..Length(res.symbols)],String);
  res.size:=Length(res.symbols);
  res.operations:=FamilyOps;
  return res;
end;

# makes  family of symbol S in spets G(e,e,r) or G(e,1,r)
MakeFamilyImprimitive:=function(S,uc)local r,f;
  f:=x->Position(uc.charSymbols,x);
  if Length(S)=1 then return Family("C1",List(S,f));fi;
  r:=FamilyImprimitive(FullSymbol(S[1]));
  r.charNumbers:=List(r.symbols,f);
  r.special:=PositionProperty(r.charNumbers,x->uc.a[x]=uc.b[x]);
  r.cospecial:=PositionProperty(r.charNumbers,x->uc.A[x]=uc.B[x]);
  if Length(DecomposedMat(r.fourierMat))>1 then Error();fi;
  return r;
end;

# JM 30/11/2000 we have to decide how to represent cuspidals of imprimitive
# groups -- the function below is an ad-hoc solution for now
#  ImprimitiveCuspidalName(<symbol>) returns the TeX name
ImprimitiveCuspidalName:=function(S) local d,r,s,p;
  r:=RankSymbol(S);d:=Length(S);s:=IntListToString(List(S,Length));
  if r=0 then return "";fi;
  if Sum(S,Length) mod d=1 then # G(d,1,r)
    if r=1 then
      if d=3 then return "Z_3";else return SPrint("Z_{",d,"}^{",s,"}");fi;
    else return SPrint("G_{",d,",1,",r,"}^{",s,"}");
    fi;
  else # G(d,d,r)
    if r=2 then
      if d=4 then return "B_2";
      elif d=6 then p:=rec(212010:=-1,221001:=1,211200:=E(3)^2,220110:=E(3));
	return SPrint("G_2[",FormatTeX(p.(s)),"]");
      else p:=CHEVIE.R("SymbolToParameter","I")(S);
	return SPrint("I_2(",d,")",FormatGAP(p));
      fi;
    elif r=3 and d=3 then p:=rec(300:=E(3),330:=E(3)^2);
      return SPrint("G_{3,3,3}[",FormatTeX(p.(s)),"]");
    elif r=3 and d=4 then p:=rec(3010:=E(4),3230:=-E(4));
      return SPrint("G_{4,4,3}[",FormatTeX(p.(s)),"]");
    else return SPrint("G_{",d,",",d,",",r,"}^{",s,"}");
    fi;
  fi;
end;

CHEVIE.AddData("UnipotentCharacters","imp",function(p,q,r)
 local uc,cusp,f,l,ci,seteig,s,extra,addextra;
 if not (q in [1,p]) then return false;fi;
 uc:=rec(charSymbols:=CHEVIE.R("CharSymbols","imp")(p,q,r));
 uc.a:=List(uc.charSymbols,LowestPowerGenericDegreeSymbol);
 uc.A:=List(uc.charSymbols,HighestPowerGenericDegreeSymbol);
 ci:=CHEVIE.R("CharInfo","imp")(p,q,r);
 if q=1 then
   cusp:=Set(List(uc.charSymbols,S->List(S,Length)-Minimum(List(S,Length))));
   cusp:=List(cusp,x->List(x,y->[0..y-1]));
   SortBy(cusp,RankSymbol);
   uc.harishChandra:=List(cusp,function(c)local cr,res;
    cr:=RankSymbol(c);
    res:=rec(levi:=[1..cr]);
    if cr<r then res.parameterExponents:=[List(c,Length)];
    else res.parameterExponents:=[];
    fi;
    Append(res.parameterExponents,[2+cr..r]*0+1);
    if r=cr then res.relativeType:=rec(series:="A",indices:=[],rank:=0);
    else         res.relativeType:=rec(series:="ST",indices:=[1+cr..r],
                                       rank:=r-cr,p:=p,q:=1);
    fi;
    res.eigenvalue:=E(24)^(-2*(p^2-1)*QuoInt(Sum(c,Length),p))*
	E(2*p)^Sum([0..p-1],i->-(i^2+p*i)*Length(c[i+1]));
    res.charNumbers:=List(List(CHEVIE.R("CharSymbols","imp")(p,1,r-cr)
      {[1..Length(PartitionTuples(r-cr,p))]},x->List(x,PartBeta)),
         x->Position(uc.charSymbols,SymbolPartitionTuple(x,List(c,Length))));
    res.cuspidalName:=ImprimitiveCuspidalName(c);
    return res;end);
    uc.b:=[];uc.B:=[];
    uc.b{uc.harishChandra[1].charNumbers}:=ci.b;
    uc.B{uc.harishChandra[1].charNumbers}:=ci.B;
    for f in uc.harishChandra{[2..Length(uc.harishChandra)]} do
     uc.b{f.charNumbers}:=f.charNumbers*0; uc.B{f.charNumbers}:=f.charNumbers*0;
    od;
    uc.families:=List(CollectBy(uc.charSymbols,x->Collected(Concatenation(x))),
      y->MakeFamilyImprimitive(y,uc));
    SortBy(uc.families,x->x.charNumbers);
    if r=1 then l:=List(uc.charSymbols{uc.families[2].charNumbers},
      function(S)local p;p:=Position(S,[]);if p=false then return 1;
                                           else return (-1)^p;fi;end);
      uc.families[2].fourierMat:=uc.families[2].fourierMat^DiagonalMat(l);
      uc.cyclicparam:=List(uc.charSymbols,function(s)
        if Number(Flat(s),x->x=1)=1 then return [1];
        else s:=Copy(s);l:=PositionProperty(s,p->1 in p);s[l]:=[];
          return [PositionProperty(s,p->1 in p)-1,l-1];fi;end);
    elif r=2 and p=3 then
      uc.families[4].fourierMat:=uc.families[4].fourierMat^DiagonalMat(-1,1,1);
      uc.families[1].fourierMat:=uc.families[1].fourierMat^
        DiagonalMat(1,-1,-1,1,1,1,1,1,1);
    fi; #Dudas' sign change
    return uc;
  elif p=q then
   uc.families:=[];
   for f in CollectBy([1..Length(uc.charSymbols)],i->
     Collected(Concatenation(FullSymbol(uc.charSymbols[i])))) do
     if Length(Set(List(uc.charSymbols{f},FullSymbol)))>1 
     then Add(uc.families,rec(charNumbers:=f));
     else Append(uc.families,List(f,x->Family("C1",[x])));
     fi;
   od;
   SortBy(uc.families,x->x.charNumbers);
   uc.harishChandra:=List(CollectBy([1..Length(uc.charSymbols)],
    function(i)local s,l;
     s:=FullSymbol(uc.charSymbols[i]);l:=List(s,Length);
     return [Sum(s,x->Sum(PartBeta(x))),l-Minimum(l)];end),
     l->rec(charNumbers:=l));
   Sort(uc.harishChandra);
   extra:=[];
   for f in uc.harishChandra do
     addextra:=false;
     s:=FullSymbol(uc.charSymbols[f.charNumbers[1]]);
     l:=r-Sum(s,x->Sum(PartBeta(x))); f.levi:=[1..l];
     s:=List(s,Length);s:=s-Minimum(s);
     f.eigenvalue:=E(24)^(-2*(p^2-1)*QuoInt(Sum(s),p))*
	E(2*p)^Sum([0..p-1],i->-(i^2+p*i)*s[i+1]);
     if l=r then f.relativeType:=rec(series:="A",indices:=[],rank:=0);
       f.parameterExponents:=[];
       if Length(f.charNumbers)=2 then addextra:=true;fi;
     elif l=0 then 
       f.relativeType:=rec(series:="ST",indices:=[1..r],rank:=r,p:=p,q:=q);
       f.parameterExponents:=[1..r]*0+1;
     else 
       f.relativeType:=rec(series:="ST",indices:=[l+1..r],
	 rank:=r-l,p:=p,q:=1);
       f.parameterExponents:=[1..r-l]*0+1;f.parameterExponents[1]:=s;
     fi;
     s:=List(s,x->[0..x-1]);
     f.cuspidalName:=ImprimitiveCuspidalName(s);
     if addextra then 
       s:=Copy(f.charNumbers);f.charNumbers:=s{[1]};
       f:=Copy(f);f.charNumbers:=s{[2]};
       Add(f.cuspidalName,'2');
       Add(extra,f);
     fi;
   od;
   Append(uc.harishChandra,extra);
   for f in uc.families do f.eigenvalues:=List(f.charNumbers,i->
     First(uc.harishChandra,s->i in s.charNumbers).eigenvalue);
   od;
   uc.b:=[];uc.B:=[];
   uc.b{uc.harishChandra[1].charNumbers}:=ci.b;
   uc.B{uc.harishChandra[1].charNumbers}:=ci.B;
   for f in uc.harishChandra{[2..Length(uc.harishChandra)]} do
    uc.b{f.charNumbers}:=f.charNumbers*0; uc.B{f.charNumbers}:=f.charNumbers*0;
   od;
   if [p,q,r]=[3,3,3] then
     uc.families[6]:=Family(ComplexConjugate(CHEVIE.families.X(3)),[8,7,11],
       rec(signs:=[1,1,-1]));
     uc.families[4]:=Family(CHEVIE.families.X(3),[4,5,12],
       rec(signs:=[1,1,-1]));
     uc.curtis:=[ 1, 2, 3, 7, 8, 10, 4, 5, 9, 6, -12, -11 ];
   elif [p,q,r]=[3,3,4] then
     uc.families[2]:=Family(CHEVIE.families.X(3),[2,4,23],
       rec(signs:=[1,1,-1]));
     uc.families[6]:=Family(CHEVIE.families.QZ(3),[13,9,8,10,19,22,7,21,20],
      rec(signs:=[1,1,1,-1,-1,1,-1,-1,1],special:=3,cospecial:=2));
     uc.families[9]:=Family(ComplexConjugate(CHEVIE.families.X(3)),[15,14,18],
       rec(signs:=[1,1,-1]));
   elif [p,q,r]=[3,3,5] then
     uc.families[3]:=Family(CHEVIE.families.X(3),[3,6,51],
       rec(signs:=[1,1,-1]));
     uc.families[4]:=Family(CHEVIE.families.X(3),[4,5,54],
       rec(signs:=[1,1,-1]));
     uc.families[6]:=Family(CHEVIE.families.QZ(3),[9,10,8,21,44,46,20,49,45],
       rec(signs:=[1,1,1,1,1,1,1,-1,-1]));
     uc.families[7]:=Family(CHEVIE.families.QZ(3),[23,11,16,12,42,50,15,48,40],
       rec(signs:=[1,-1,-1,1,1,1,1,-1,-1],special:=4,cospecial:=7));
     uc.families[8]:=Family(ComplexConjugate(CHEVIE.families.X(3)),[14,13,41],
       rec(signs:=[1,1,-1]));
     uc.families[11]:=Family(CHEVIE.families.X(3),[19,22,47],
       rec(signs:=[1,1,-1]));
     uc.families[13]:=Family(CHEVIE.families.QZ(3),[32,27,26,28,38,53,25,52,39],
       rec(signs:=[1,1,1,-1,-1,1,-1,-1,1],special:=3,cospecial:=2));
     uc.families[15]:=Family(ComplexConjugate(CHEVIE.families.X(3)),[31,30,37],
       rec(signs:=[1,1,-1]));
     uc.families[16]:=Family(ComplexConjugate(CHEVIE.families.X(3)),[34,33,43],
       rec(signs:=[1,1,-1]));
   elif [p,q,r]=[4,4,3] then
     uc.families[2]:=Family(CHEVIE.families.X(4),[3,2,4,14,16,13],
	 rec(signs:=[1,1,1,1,-1,-1]));
     uc.families[4]:=Family(ComplexConjugate(CHEVIE.families.X(4)),
       [8,6,7,12,15,11],rec(signs:=[1,1,1,1,-1,-1]));
     uc.curtis:=[1,6,7,8,10,2,3,4,9,5,14,13,12,11,-16,-15];
   elif [p,q,r]=[4,4,4] then
     uc.families[5]:=Family(CHEVIE.families.X(4),[5,8,9,46,53,47],
 	rec(signs:=[1,1,1,-1,-1,1]));
     uc.families[6]:=Family("C2",[12,7,6,42]);
     uc.families[7]:=Family(CHEVIE.families.X(4),[13,10,11,41,55,43],
      rec(signs:=[1,1,1,1,1,-1],special:=3,cospecial:=1));
     uc.families[9]:=Family(CHEVIE.families.QZ(4),
     [18,21,28,22,23,49,39,54,56,40,15,36,19,52,37,51],
     rec(signs:=[1,1,1,1,1,1,-1,-1,1,-1,-1,-1,1,1,-1,-1],
     special:=2,cospecial:=4));
    uc.families[10]:=Family(ComplexConjugate(CHEVIE.families.X(4)),
     [16,17,20,38,50,34],rec(signs:=[1,1,1,-1,1,1],special:=3,cospecial:=1));
     uc.families[12]:=Family("C2",[27,26,25,35]);
     uc.families[13]:=Family(ComplexConjugate(CHEVIE.families.X(4)),
       [30,29,31,44,48,45],rec(signs:=[1,1,1,1,1,-1],special:=3,cospecial:=1));
   else uc.families:=List(uc.families,x->MakeFamilyImprimitive(
     uc.charSymbols{x.charNumbers},uc));
   fi; 
   return uc;
  fi;
end);

CHEVIE.AddData("Invariants","imp",function(p,q,r)local v;
  v:=List([1..r-1],i->function(arg)
    return Sum(Arrangements([1..r],i),a->Product(arg{a})^p);end);
  Add(v,function(arg)return Product(arg)^(p/q);end);
  return v;
end);
