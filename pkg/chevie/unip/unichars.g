#############################################################################
##
#A  unichar.g    Unipotent functions                           Jean Michel
##
#Y  Copyright 1996-2010, Univ. Paris VII
##  
##  This file contains functions for linear combinations of unipotent 
##  characters cotrresponding to Coxeter groups, Coxeter cosets and Spets.
##  
UnipCharOps:=OperationsRecord("UnipCharOps");
IsUnipotentCharacter:=x->IsRec(x) and 
  IsBound(x.operations) and x.operations=UnipCharOps;

CHEVIE.PrintUniChars:=rec(short:=1);

UnipCharOps.Format:=function(r,option)local m,c,i,s,n,res;
  res:="";
  option:=Inherit(ShallowCopy(CHEVIE.PrintUniChars),option);
  s:=CharNames(UnipotentCharacters(r.group),option);
  m:=Maximum(List(s,Length))+3;
  for i in [1..Length(r.v)] do
    n:=ConcatenationString("<",s[i],">");
    c:=Format(r.v[i]);
    if IsBound(option.short) then
      if c<>"0" then
	if c="1" then Append(res,"+");
	elif c="-1" then Append(res,"-");
	else  if '+' in c{[2..Length(c)]} or '-' in c{[2..Length(c)]} then
		  c:=Concatenation("(",c,")");
	      fi;
	      if not c[1] in "+-" then Append(res,"+");fi;
	      Append(res,c);
	fi;
	Append(res,n);
      fi;
    elif c<>"0" or not IsBound(option.nozero) then
      PrintToString(res,"\n",n,String("",m-Length(n)),c);
    fi;
  od;
  if Length(res)=0 then res:="0";fi;
  if res[1]='+' then res:=res{[2..Length(res)]};fi;
  if IsBound(r.name) then 
    res:=SPrint("DLvar[",ReflectionName(r.group),",",r.name,"]=",res);
  else
    res:=SPrint("[",ReflectionName(r.group),"]=",res);
  fi;
  return String(res);
end;

UnipCharOps.Display:=function(r,opt)
  Print(Format(r,opt),"\n");
end;

UnipCharOps.String:=r->Format(r,CHEVIE.PrintUniChars);

UnipCharOps.Print:=function(r)Print(String(r));end;

UnipCharOps.\-:=function(a,b)return UnipotentCharacter(a.group,a.v-b.v);end;

UnipCharOps.\+:=function(a,b)return UnipotentCharacter(a.group,a.v+b.v);end;

UnipCharOps.\*:=function(a,b)local res;
  if IsUnipotentCharacter(b) then
     if IsUnipotentCharacter(a) then return a.v*b.v;
     else return UnipotentCharacter(b.group,a*b.v);
     fi;
  else return UnipotentCharacter(a.group,b*a.v);
  fi;
end;

UnipCharOps.Degree:=function(a)
   return a.v*UnipotentDegrees(a.group,X(Cyclotomics));
end;

UnipotentCharacter:=function(W,v)local r,n;
  n:=CharNames(UnipotentCharacters(W));
  if IsInt(v) then r:=List(n,x->0);r[v]:=1;v:=r;
  elif IsString(v) then  r:=List(n,x->0);r[Position(n,v)]:=1;v:=r;
  fi;
  return rec(v:=v,group:=W,operations:=UnipCharOps);
end;

# induce unip char u to WF
LusztigInduction:=function(WF,u)local t;
  t:=LusztigInductionTable(u.group,WF);
  if t=false then return false; 
  else return UnipotentCharacter(WF,t.scalar*u.v);
  fi;
end;

# Restrict unip char u to HF
LusztigRestriction:=function(HF,u)
  return UnipotentCharacter(HF,u.v*LusztigInductionTable(HF,u.group).scalar);
end;

# induce unip char u to WF
HarishChandraInduction:=function(WF,u)return 
  UnipotentCharacter(WF,HarishChandraInductionTable(u.group,WF).scalar*u.v);
end;

# Restrict unip char u to HF
HarishChandraRestriction:=function(HF,u)return 
  UnipotentCharacter(HF,u.v*HarishChandraInductionTable(HF,u.group).scalar);
end;

DeligneLusztigCharacterTable:=function(W)local uc;
  if not IsBound(W.rwTable) then
    uc:=UnipotentCharacters(W);
    W.rwTable:=TransposedMat(ComplexConjugate(CharTable(W).irreducibles))*
      Fourier(uc){uc.almostHarishChandra[1].charNumbers};
  fi;
  return W.rwTable;
end;
  
#      DeligneLusztigCharacter(W,<cox. word>)
#   or DeligneLusztigCharacter(W,<perm>)
#   or DeligneLusztigCharacter(W,integer i) (for ith class)
#   W may be a group or a Spets.
DeligneLusztigCharacter:=function(W,w)local ct,uc;
  if IsList(w) then w:=EltWord(W,w);fi;
  if IsPerm(w) then w:=PositionClass(W,w);fi;
  return UnipotentCharacter(W,DeligneLusztigCharacterTable(W)[w]);
end;

AlmostCharacter:=function(W,i)local ct,dl;
  dl:=List([1..NrConjugacyClasses(W)],i->DeligneLusztigCharacter(W,i));
  ct:=CharTable(W);
  return Zip(ct.irreducibles[i],ct.classes,
            function(x,y)return x*y;end)/Size(W)*dl;
end;

# usage DeligneLusztigLefschetz(hecke element[,i])
# if not i given |X_T^{F^m}| for divisible h
# if i given use eigenvalues^i
DeligneLusztigLefschetz:=function(arg)local W,h,ct,uc,i;
  h:=arg[1];
  if IsBound(h.coset) 
  then W:=ReflectionCoset(h.coset);
  else W:=Group(Hecke(h));
  fi;
  uc:=UnipotentCharacters(W);
  if Length(arg)=2 then i:=arg[2];else i:=0;fi;
  return UnipotentCharacter(W,Zip(ComplexConjugate(HeckeCharValues(h))*
    uc.operations.Fourier(uc){uc.almostHarishChandra[1].charNumbers},
      List(uc.operations.Eigenvalues(uc,[1..Size(uc)]),x->x^i),
      function(a,b)return a*b;end));
end;

DeligneLusztigLefschetzTable:=function(H)local WF,t,uc;
  if IsBound(H.spets) then WF:=ReflectionCoset(H);
  else WF:=Group(H);
  fi;
  t:=CharTable(H).irreducibles;
# t:=t*(1+0*Sum(t,Sum));
  uc:=UnipotentCharacters(WF);
  return ComplexConjugate(TransposedMat(t))*
    uc.operations.Fourier(uc){uc.almostHarishChandra[1].charNumbers};
end;

UnipCharOps.Frobenius:=function(WF,x,i)local W,p,uc,t,pt;
  W:=x.group;
  p:=List(ConjugacyClasses(W),x->PositionClass(W,Representative(x)^WF.phi));
  p:=PermList(p);
  uc:=UnipotentCharacters(W);
  t:=DeligneLusztigCharacterTable(W);Add(t,Eigenvalues(uc));
  pt:=Permuted(t,p);
  t:=TransposedMat(t);
  pt:=TransposedMat(pt);
  p:=List(pt,x->Positions(t,x));
  if ForAny(p,x->Length(x)>1) then 
    Error("Rw + eigen cannot disambiguate\n");
  fi;
  p:=PermList(List(p,x->x[1]));
  x:=ShallowCopy(x);
  x.v:=Permuted(x.v,p^-i);
  return x;
end;
