##############################################################################
##
#A  pres.g       VKCURVE package         Jean Michel
##
#Y  Copyright (C) 2001 - 2002  University Paris VII, France.
##
##  This file holds various functions to deal with presentations.
## 
#############################################################################
ShrinkPresentation:=function(arg)local g,lim,rot,count,test;
  g:=arg[1];
  if Length(arg)=2 then lim:=arg[2];else lim:=1000;fi;
  rot:=function(i)local l;l:=g.tietze[6][i];
    l:=l{Concatenation([2..Length(l)],[1])};g.tietze[6][i]:=l;
  end;
  test:=function()local v,i,j,before,tt,t;
    g.tietze[6]:=Filtered(g.tietze[6],x->Length(x)>0);
    tt:=g.tietze[6];
    if Product(tt,Length)<lim then
      v:=List(tt,x->0);
      if Length(v)=0 then return false;fi;
      while true do
	before:=[Length(tt),Sum(tt,Length)];
	j:=Length(v);
	while v[j]=Length(tt[j])-1 do 
	  rot(j);v[j]:=0;j:=j-1;
	  if j=0 then return false;fi;
        od;
	rot(j);v[j]:=v[j]+1;
	TzGoGo(g);
	tt:=g.tietze[6];
	if Length(tt)<before[1] or Sum(tt,Length)<before[2] then return true;fi;
      od;
    else 
      for i in [1..lim] do
	t:=Sum(tt,Length);
	i:=Random([1..t]);
	j:=1;
	while i>Length(tt[j]) do i:=i-Length(tt[j]);j:=j+1;od;
	before:=[Length(tt),t];rot(j);TzGoGo(g);
	tt:=g.tietze[6];
	if Length(tt)<before[1] or Sum(tt,Length)<before[2] then return true;fi;
      od;
    fi;
    return false;
  end;
  count:=0; if g.tietze[6]<>[] then while test() do count:=count+1;od; fi;
end;
    
DisplayPresentation:=function(arg)local g,min,maj,l,i,w,n,used,f,lw,m;
  f:=function(i)
    i:=(i mod lw)+1;
    return Concatenation(w{[i..Minimum(i+lw/2-1,lw)]},w{[1..i-lw/2-1]});
  end;
  min:="abcdefghijklmnopqrstuvwxyz";
  maj:="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  used:=[];
  g:=arg[1];
  l:=List(g.tietze[6],x->List(x,function(y)
    AddSet(used,AbsInt(y));
    if y<0 then return maj[-y];else return min[y];fi;
    end));
  if g.tietze[1]>Length(used) then
    Print("There are ",g.tietze[1]-Length(used)," free generators\n");
  fi;
  for i in [1..Length(l)] do
    w:=l[i];lw:=Length(w);
    if lw mod 2=1 then Print(i,": ",w,"=1\n");
    elif lw>0 then
      m:=List([1..lw],i->Number(f(i),x->x in min));
      n:=Maximum(m);
      m:=PositionProperty(m,i->i=n);
      Print(i,": ",f(m),"=",List(Reversed(f(m+lw/2)),
	function(x) if x in maj then return min[Position(maj,x)];
		    else             return maj[Position(min,x)];
		    fi;
		    end),"\n");
    fi;
  od;
end;

# <presentation>, "abA"
ConjugatePresentation:=function(p,s)local f,x,l,n,minmaj;
 minmaj:="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
 l:=List(s,l->Position(minmaj,l));
 f:=FpGroupPresentation(p);
 x:=AbstractGenerator("x");
 if l[1]-l[3]>0 then n:=f.(l[3])^-1*x*f.(l[3]);
 else n:=f.(l[1])*x*f.(l[1])^-1;
 fi;
 f.relators:=List(f.relators,y->EliminatedWord(y,f.(l[2]),n));
 f.relators:=List(f.relators,y->EliminatedWord(y,x,f.(l[2])));
 f:=SimplifiedFpGroup(f);
 p:=PresentationFpGroup(f);
 ShrinkPresentation(p,100);
 return p;
end;

TryConjugatePresentation:=function(arg)
  local expand,p1,p2,m,n,Tally,applicable,minmaj,tp,c,p,toletter;
  minmaj:="abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
  toletter:=function(n)
    if n>0 then return minmaj[n];
    else return minmaj[-n+26];
    fi;
  end;
# s: "Cac" means take a->Cac

  Tally:=p->Sum( p.tietze[6], Length);

  applicable:=function(p)local res,v,i,r;
   res:=[];
   for v in p.tietze[6] do
     for i in [1..Length(v)-2] do
       if v[i]=-v[i+2] then
         if v[i+1]>0 then r:=List(v{[i..i+2]},toletter);
         else r:=List(-v{[i..i+2]},toletter);
	 fi;
	 if not r in res then Add(res,r);fi;
       fi;
     od;
   od;
   return List(res,String);
  end;

  expand:=function(p)local res;
    res:=List(applicable(p),x->[x,ConjugatePresentation(p,x)]);
    SortParallel(List(res,x->Tally(x[2])),res);
    return res;
  end;

  p:=Copy(arg[1]);
  TzGoGo(p);
  if Length(arg)>=2 then tp:=arg[2];else tp:=Tally(p);fi;
  p1:=[];
  for c in applicable(p) do
    Print(c);
    n:=ConjugatePresentation(p,c);
    if Length(arg)=3 and Tally(n)<=arg[3] then
      InfoChevie2("# ",[c[2]],"->",c," gives length ",Tally(n),"\n");
      DisplayPresentation(n);
    fi;
    if Tally(n)<tp then
      InfoChevie("# ",[c[2]],"->",c," gives length ",Tally(n),"\n");
      return n;
    fi;
    if Tally(n)<tp+10 then Add(p1,[c,n]);fi;
  od;
  SortParallel(List(p1,x->Tally(x[2])),p1);
  for p2 in p1 do 
    for c in applicable(p2[2]) do
      n:=ConjugatePresentation(p2[2],c);
      Print(p2[1],",",c);
      if Length(arg)=3 and Tally(n)<=arg[3] then
	InfoChevie2("# ",[p2[1][2]],"->",p2[1]," then ",[c[2]],"->",c,
	  " gives length ",Tally(n),"\n");
	DisplayPresentation(n);
      fi;
      if Tally(n)<tp then
        InfoChevie("# ",[p2[1][2]],"->",p2[1]," then ",[c[2]],"->",c,
         " gives length ",Tally(n),"\n");
        return n;
      fi;
    od;
  od;
  Print("# could not shrink presentation");
  return p;
end;
