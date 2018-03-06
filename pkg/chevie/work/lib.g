# JM's library of useful functions
2A:=n->RootDatum("psu",n+1); 2D:=n->RootDatum("2D",n); 2E6:=RootDatum("2E6");
3D4:=RootDatum("3D4"); 2B2:=RootDatum("2B2"); 2G2:=RootDatum("2G2");
2F4:=RootDatum("2F4"); 2I:=n->RootDatum("2I",n); 2G5:=Spets("2G5");
3G333:=Spets("3G333"); 3pG333:=Spets("3pG333"); 4G333:=Spets("4G333");
3G422:=Spets("3G422");

cmpvec:=CHEVIE.Check.EqLists;
cmpobj:=CHEVIE.Check.EqObj;
cmptbl:=CHEVIE.Check.EqTables;

#IsZero:=x->x=0*x;

# order of the center of the complex reflection group W
OrderCenter:=W->Product(List(W.type,t->Gcd(ReflectionDegrees(t))));

# elements of Cartesian([l1,l2,...,ln]) with all coordinates distinct
DistinctCartesian:=function(l)
  if Length(l)=1 then return List(l[1],x->[x]);fi;
  return Concatenation(List(DistinctCartesian(l{[2..Length(l)]}),
    r->List(Difference(l[1],r),e->Concatenation([e],r))));
end;

# decompose a character in term of irreducibles
DecomposeChar:=function(W,char)local ct;ct:=CharTable(W);
  return List(ct.irreducibles,x->ScalarProduct(ct,x,char));
end;

ProperStableParabolics:=function(WF)local W,l;
  if IsSpets(WF) then W:=Group(WF);else W:=WF;fi;
  l:=Combinations(W.rootInclusion{W.generatingReflections});
  if IsSpets(WF) then l:=Filtered(l,x->OnSets(x,WF.phi)=x);fi;
  SortBy(l,Length);
  return Difference(l,[W.rootInclusion{W.generatingReflections}]);
end;

ProperParabolics:=function(WF)local W,l;
  if IsSpets(WF) then W:=Group(WF);else W:=WF;fi;
  l:=Combinations(W.rootInclusion{W.generatingReflections});
  SortBy(l,Length);
  return l{[1..Length(l)-1]};
end;

braidmatrix:=function(l)
  Print(Format(List(l,i->List(l,j->BraidRelation(i,j)))),"\n");
end;
  
# Etend GaloisCyc (et l'action d'un element d'un groupe de Galois
# qui est un HomomorphismByImages) a marcher avec des listes, matrices, etc...
GaloisAction:=function(x,e)
  if IsList(x) then return List(x,y->GaloisAction(y,e));
  elif IsPolynomial(x) then
    x:=ShallowCopy(x);x.coefficients:=GaloisAction(x.coefficients,e);return x;
  elif IsMvp(x) then x:=ShallowCopy(x);
    x.coeff:=GaloisAction(x.coeff,e);return x;
  elif IsRatFrac(x) then x:=ShallowCopy(x);
    x.num:=GaloisAction(x.num,e);x.den:=GaloisAction(x.den,e);
    return x;
  elif IsInt(e) then return GaloisCyc(x,e);
  else return x^e;
  fi;
end;

# Given lists l1 and l2 of vectors for each vector in l2 find
# which vectors of l1 define the same line and proportionality coefficients.
permcoeff:=function(l1,l2)local res,v,c,d;
  res:=[];
  for v in l2 do
    c:=List(l1,x->ProportionalityCoefficient(v,x));
    d:=Filtered([1..Length(l1)],i->c[i]<>false);
    Add(res,[d,c{d}]);
  od;
  return res;
end;

####################  TeX #################################
TeXs:=function(arg)local name,s;
  s:=arg[1];
  if Length(arg)=1 then name:="tmp.tex";
  else name:=Concatenation(arg[2],".tex");fi;
  PrintTo(name,"\\tabskip=0.1cm\n");
  AppendTo(name,"\\input amssym.def\n\\input amssym.tex\n");
  AppendTo(name,"\\def\\BZ{{\\Bbb Z}}\n");
  Cut(s,rec(before:="_{+-<",after:="}=&",file:=name));
  AppendTo(name,"\\bye\n");
  Exec(Concatenation("tex ",name));
  Exec(Concatenation("xdvi -expert -s 5 ",Replace(name,".tex",".dvi"),"&"));
end;

TeX:=function(arg)local name,s,o;
  o:=rec(TeX:=1);
  if Length(arg)>1 and IsRec(arg[2]) then
    Inherit(o,arg[2]);arg:=Concatenation([arg[1]],arg{[3..Length(arg)]});
  fi;
  arg[1]:=Format(arg[1],o);
  ApplyFunc(TeXs,arg);
end;

LaTeXs:=function(arg)local name,s;
  s:=arg[1];
  if Length(arg)=1 then name:="tmp.tex";
  else name:=Concatenation(arg[2],".tex");fi;
  PrintTo(name,"\\documentclass{article}\n");
  AppendTo(name,"\\usepackage{amsmath}\n");
  AppendTo(name,"\\usepackage{amssymb}\\newcommand{\\BZ}{{\\mathbb Z}}\n");
  AppendTo(name,"\\usepackage{array}\\newcolumntype{R}{>{$}r<{$}}\n");
  AppendTo(name,"\\usepackage{supertabular}\n");
  AppendTo(name,"\\begin{document}\n");
  Cut(s,rec(before:="_{+-<",after:="}=&",file:=name));
  AppendTo(name,"\\end{document}\n");
  Exec(Concatenation("latex ",name));
  Exec(Concatenation("xdvi -expert -s 5 ",Replace(name,".tex",".dvi"),"&"));
end;

LaTeX:=function(arg)local name,s,o;
  o:=rec(TeX:=1,LaTeX:=1);
  if Length(arg)>1 and IsRec(arg[2]) then
    Inherit(o,arg[2]);arg:=Concatenation([arg[1]],arg{[3..Length(arg)]});
  fi;
  arg[1]:=Format(arg[1],o);
  ApplyFunc(LaTeXs,arg);
end;

####################  Structural functions #################################
# dump record
dr:=function(r)local p,f;
  f:=RecFields(r);
  Print("rec(");
  for p in f do 
    Print(p,":=");
    if IsFunc(r.(p)) then Print("<<function>>"); 
    else Print(r.(p)); 
    fi;
    if p<>f[Length(f)] then Print(",\n");fi;
  od;
  Print(")\n");
end;

# another way to dump record
dp:=function(r)local p;
  Print("rec(");
  for p in RecFields(r) do
    Print(p,":=",Format(r.(p)),"\n");
  od;
  Print(")\n");
end;

ReadChv("test/check");

# compare two cycpols.
cmpcycpol:=function(a,b,na,nb)local q;
 if a=b then return;fi;
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
