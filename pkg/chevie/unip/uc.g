#############################################################################
##
#A  uc.g      Unipotent Characters    Frank Luebeck and Jean Michel
##
#Y  Copyright 1996-1998, Rwth Aachen and Univ. Paris VII
##  
##  This file contains basic functions for unipotent characters of
##  Coxeter groups, Coxeter cosets and Spets.
##  

UnipotentCharactersOps:=OperationsRecord("UnipotentCharOps");

UnipotentCharactersOps.Print:=function(r)
  Print("UnipotentCharacters( ", ReflectionName(r.group), " )");
end;

# Fourier times the vector of unip. degrees is the vector of fake degrees
UnipotentCharactersOps.Fourier:=function(uc)local f;
  if not IsBound(uc.fourier) then
    uc.fourier:=IdentityMat(Size(uc));
    for f in uc.families do 
      uc.fourier{f.charNumbers}{f.charNumbers}:=f.fourierMat;
    od;
  fi;
  return uc.fourier;
end;

# FourierInverse times the vector of fake degrees is the vector of unip degrees
UnipotentCharactersOps.FourierInverse:=function(uc)local f;
  if not IsBound(uc.fourierinverse) then
    uc.fourierinverse:=IdentityMat(Size(uc));
    for f in uc.families do 
      uc.fourierinverse{f.charNumbers}{f.charNumbers}:=
	  TransposedMat(ComplexConjugate(f.fourierMat));
    od;
  fi;
  return uc.fourierinverse;
end;

UnipotentCharactersOps.Eigenvalues:=function(arg)local uc,f;
  uc:=arg[1];
  if not IsBound(uc.eigenvalues) then
    uc.eigenvalues:=[];
    for f in uc.families do uc.eigenvalues{f.charNumbers}:=f.eigenvalues;od;
  fi;
  if Length(arg)=2 then return uc.eigenvalues{arg[2]};fi;
  return uc.eigenvalues;
end;

UnipotentCharactersOps.qEigen:=function(uc)local f;
  if not IsBound(uc.qEigen) then
    uc.qEigen:=[];
    for f in uc.harishChandra do
      if IsBound(f.qEigen) 
      then uc.qEigen{f.charNumbers}:=f.charNumbers*0+f.qEigen;
      else uc.qEigen{f.charNumbers}:=f.charNumbers*0;
      fi;
    od;
  fi;
  return uc.qEigen;
end;

UnipotentCharactersOps.RelativeHecke:=function(uc,i,q)local hw,t;
  hw:=uc.harishChandra[i];
# for t in hw.relativeType do t.operations:=ReflTypeOps;od;
  return Hecke(ApplyFunc(ReflectionGroup,hw.relativeType),
    List(hw.parameterExponents,function(i)
    if IsList(i) then return List([1..Length(i)],j->E(Length(i))^(j-1)*q^i[j]);
#   elif i<0  then return -q^-i; # JM 14/2/2018 I think obsolete
    else return q^i;
    fi; end));
end;

# next is a field to hide its name
UnipotentCharactersOps.ParamsAndNames:=function(sers)local res,ser,t,n,params;
  res:=rec(charParams:=[],TeXCharNames:=[]);
  for ser in sers do
    t:=ser.relativeType;n:=ser.cuspidalName;
    params:=CHEVIE.Data("CharInfo",t).charparams;
    res.charParams{ser.charNumbers}:=List(params,x->[n,x]);
    res.TeXCharNames{ser.charNumbers}:=List(params,function(x)local s;
      s:=Copy(n);
      if IsBound(t.orbit) then t.rank:=t.orbit[1].rank;fi;
      if Length(s)>0 and t.rank>0 then Add(s,':');fi;
      if t.rank>0 then  Append(s,CHEVIE.Data("CharName",t,x,rec(TeX:=1)));fi;
      return String(s);
      end);
  od;
  return res;
end;

UnipotentCharactersOps.CharNames:=function(uc,opt)
  if IsBound(opt.cyclicparam) and IsBound(uc.cyclicparam) then
    return List(uc.cyclicparam,function(x)
     if Length(x[1])=1 then return "Id";
     else return TeXStrip(SPrint("\\rho_{",x[1][1],",",x[1][2],"}"),opt);
     fi;end);
  else return List(uc.TeXCharNames,x->TeXStrip(x,opt));
  fi;
end;

FixRelativeType:=function(t)
# fix illegal relativeTypes B1 and C2 which appear in HC or almost HC
# series of classical groups
  if t.relativeType.series="B" then
    if t.relativeType.rank=1 then
      t.relativeType.series:="A";
      t.charNumbers{[1,2]}:=t.charNumbers{[2,1]}; # map B1->A1
    elif t.relativeType.rank=2 and IsBound(t.relativeType.cartanType) and
      t.relativeType.cartanType=1 then
      t.relativeType.cartanType:=2;
      t.relativeType.indices:=Reversed(t.relativeType.indices);
      t.charNumbers{[1,5]}:=t.charNumbers{[5,1]}; # map C2->B2
      if IsBound(t.parameterExponents) then
	t.parameterExponents:=Reversed(t.parameterExponents);
      fi;
    fi;
  fi;
end;

#############################################################################
##
#F  UnipotentCharacters( <rdr> ) . . . . . information about the
#F  unipotent characters of G(q) given by its complete root datum.
##  Works for Spets.
##  
#   The fields are as follows:
#
#   harishChandra: the Harish-Chandra series
#    each is a record with fields
#     levi:  The indices in the  diagram of W of  the Levi L^F where is the
#        cuspidal (Given in terms of W's parent)
#     cuspidalName: The name of that cuspidal
#     eigenvalue:   Its Frobenius Eigenvalue
#     relativeType: A record describing W_G^F(L). In addition to series and
#        rank,  indices  is  the  index  in  W's parent of reflexions whose
#        F-orbits generate W_G^F(L)
#     parameterExponents:  In the  same order  as relativeType.indices, the
#        exponents  of  the  powers  of  q  which are the parameters of the
#        relative Hecke algebra.
#     charNumbers:  the  numbers  of  the  unipotent  characters  in the HC
#        series, in the same order as the characters of W_G^F(L)
#
#   almostHarishChandra: the almost Harish-Chandra series
#    each is a record with fields
#     levi:  The indices  in the  diagram of  W of  the Levi L where is the
#        cuspidal character sheaf (Given in terms of W's parent)
#     cuspidalName: The name of that cuspidal
#     eigenvalue:   Its Frobenius Eigenvalue
#     relativeType:  A record describing (W_G(L),F).  In addition to series
#        and  rank, indices is the index  in W's parent of reflexions which
#        generate W_G(L), and the action of Frobenius is described.
#     charNumbers:  the numbers of the characters  in the almost HC series,
#        in the same order as the characters of W_G(L)
#
#   The argument must be the reflectiontype  of a Spets

ReflTypeOpsUnipotentCharacters:=function(t) local uc,a,s,i,indices,r;
  if t.orbit[1].series="B" and IsBound(t.orbit[1].cartanType) 
     and t.orbit[1].cartanType=1 then 
       uc:=CHEVIE.Data("UnipotentCharacters",t,t.orbit[1].cartanType);
  else uc:=CHEVIE.Data("UnipotentCharacters",t);
  fi;
  if uc=false then 
    Print("Warning: ",ReflectionName(t)," Is not a Spets!!");
    return uc;
  fi;
  Inherit(uc,UnipotentCharactersOps.ParamsAndNames(uc.harishChandra));
  if not IsBound(uc.charSymbols) then
    uc.charSymbols:=uc.charParams;
  fi;
  if false then
  # JM experimental code 4-2007: adjust principal series for scalar-twisted
  # should work only if order of scalar prime to ZW
  # ChevieErr(ReflectionName(t)," Ennola-twisted: eigenvalues may be wrong\n");
  if IsBound(t.scalar) and t.scalar[1]<>1 then
    if Length(t.orbit)>1 then Error("pas envisage");fi;
    uc.families:=Copy(uc.families);
    for s in uc.families do
      a:=List(s.charNumbers,i->t.scalar[1]^(-uc.a[i]-uc.A[i]));
      s.eigenvalues:=Zip(s.eigenvalues,a,function(x,y)return x*y;end);
      s.fourierMat:=TransposedMat(Zip(TransposedMat(s.fourierMat),a,
            function(x,y)return x*y;end));
    od;
  fi;
  fi;
  # adjust things for descent of scalars
  # we would like to adjust indices so they fit with those stored in t
  # but we cannot when indices mention non-generating reflections!
  a:=Length(t.orbit);
  if a>1 then
    if IsBound(uc.a) then uc.a:=a*uc.a;fi;
    if IsBound(uc.A) then uc.A:=a*uc.A;fi;
    for s in uc.harishChandra do
      s.parameterExponents:=a*s.parameterExponents;
      s.eigenvalue:=s.eigenvalue^a;
      s.cuspidalName:=Join(List([1..a],i->s.cuspidalName),"\\otimes ");
      s.relativeType.operations:=ReflTypeOps;
    od;
  else
    for s in uc.harishChandra do
      s.relativeType.operations:=ReflTypeOps;
    od;
  fi;
# indices:=Concatenation(List(t.orbit,x->x.indices));
# for s in uc.harishChandra do
#   s.levi:=indices{s.levi};
#   s.relativeType.indices:=indices{s.relativeType.indices};
# od;
  if not IsBound(uc.almostHarishChandra) then
    uc.almostHarishChandra:=List(uc.harishChandra, function(s)local res,a;
      res:=rec();
      Inherit(res,s,["levi", "cuspidalName", "eigenvalue", "charNumbers"]);
      res.relativeType:=rec(orbit:=[ShallowCopy(s.relativeType)],twist:=());
      if t.twist<>() then
	a:=t.orbit[1].indices{s.relativeType.indices};
	res.relativeType.twist:=MappingPermListList(a,OnTuples(a,t.twist));
      fi;
      return res;
      end);
  else
    for s in uc.almostHarishChandra do
      s.relativeType.operations:=ReflTypeOps;
      if not IsBound(s.relativeType.orbit) then 
	s.relativeType:=rec(orbit:=[s.relativeType],twist:=());
      fi;
    od;
  fi;
  if not IsBound(uc.almostCharSymbols) then
    uc.almostCharSymbols:=uc.charSymbols;
  fi;
  a:=UnipotentCharactersOps.ParamsAndNames(uc.almostHarishChandra);
  for s in [["almostCharParams","charParams"],["almostTeXCharNames","TeXCharNames"]]
  do uc.(s[1]):=a.(s[2]);od;
  return uc;
end;

HasTypeOpsUnipotentCharacters:=function(WF) 
  local CartesianSeries, type, simp, t, H, s, a, r, f, tmp, uc, i, res, W;

  CartesianSeries:=function(sers)local res;
    res:=rec();
    res.levi:=Concatenation(List(sers,x->x.levi));
    res.relativeType:=Filtered(List(sers,x->x.relativeType),x->x.rank<>0);
    if IsBound(sers[1].eigenvalue) then
      res.eigenvalue:=Product(List(sers,x->x.eigenvalue));
    fi;
    if ForAny(sers,x->IsBound(x.qEigen)) then
      res.qEigen:=Sum(List(sers,function(x)
       if not IsBound(x.qEigen) then return 0;
       elif x.qEigen=false then return false;
       else return x.qEigen;
       fi;end));
    else 
      res.qEigen:=0;
    fi;
    if ForAll(sers,x->IsBound(x.parameterExponents)) then
      res.parameterExponents:=Concatenation(List(sers,x->x.parameterExponents));
    fi;
    res.charNumbers:=Cartesian(List(sers,x->x.charNumbers));
    res.cuspidalName:=Join(List(sers,x->x.cuspidalName),"\\otimes ");
    return res;
  end;

  type:=ReflectionType(WF);
  if Length(type)=0 then
  # UnipotentCharacters(CoxeterGroup());
    return rec( 
      harishChandra:=[
	rec(relativeType:=[], 
	    levi:=[], parameterExponents:=[],
	    cuspidalName:="", eigenvalue:=1, charNumbers :=[ 1 ])],
      families := [Family("C1",[1])],
      charParams := [ [ "", [ 1 ] ] ],
      TeXCharNames := [ "" ],
      charSymbols := [ [ "", [ 1 ] ] ],
      size:=1,
      almostHarishChandra:=[
        rec(relativeType:=[], levi:=[], 
            cuspidalName := "", eigenvalue:=1, charNumbers := [1])],
      almostCharSymbols := [ [ "", [ 1 ] ] ],
      almostCharParams := [ [ "", [ 1 ] ] ],
      almostTeXCharNames := [ "" ],
      a := [ 0 ],
      A := [ 0 ],
      group:=WF,
      operations:=UnipotentCharactersOps);
  fi;

  simp:=[]; W:=Group(WF);
  for t in ReflectionType(WF) do
# adjust indices of Levis, almostLevis, relativetypes so they agree with
# Parent(Group(WF))
    uc:=ReflTypeOps.UnipotentCharacters(t);
    if uc=false then return false;fi;
    H:=List(t.orbit,x->ReflectionSubgroup(W,W.rootInclusion{x.indices}));
    for s in uc.harishChandra do
      s.levi:=Concatenation(List(H,x->x.rootInclusion{s.levi}));
      s.relativeType.indices:=H[1].rootInclusion{s.relativeType.indices};
    od;
    for s in uc.almostHarishChandra do
      s.levi:=Concatenation(List(H,x->x.rootInclusion{s.levi}));
      s.relativeType.orbit:=Concatenation(List(H,x->
        List(s.relativeType.orbit,function(r)
	  r:=ShallowCopy(r);
	  r.indices:=x.rootInclusion{r.indices};
	  return r;end)));
      s.relativeType.twist:=s.relativeType.twist^MappingPermListList(
         [1..Length(H[1].roots)],H[1].rootInclusion);
    od;

    Add(simp,uc);
  od;

  # "Kronecker product" of records in simp:
  r:=simp[1];f:=RecFields(r);res:=rec();
  for a in f do
    if Length(simp)=1 then res.(a):=List(r.(a),x->[x]);
    elif ForAll(simp,x->IsBound(x.(a))) then
      res.(a):=Cartesian(List(simp,x->x.(a)));
    fi;
  od;
  
  for a in ["TeXCharNames","almostTeXCharNames"]
  do res.(a):=List(res.(a),x->Join(x,"\\otimes "));od;

  res.size:=Length(res.TeXCharNames);
  
  for a in [ "harishChandra", "almostHarishChandra" ] 
  do res.(a):=List(res.(a),CartesianSeries);od;

  if Length(res.families[1])=1 then
    res.families:=List(res.families,function(f)
      f:=f[1];f.charNumbers:=List(f.charNumbers,x->[x]);return f;end);
  else res.families:=List(res.families,x->ApplyFunc(FamilyOps.\*,x));
  fi;
  
  for a in ["a", "A"] do 
    if IsBound(res.(a)) then res.(a):=List(res.(a),Sum);fi;
  od;

  # finally the new 'charNumbers' lists
  tmp:=Cartesian(List(simp,a->[1..Length(a.TeXCharNames)]));
  for a in [ "harishChandra", "almostHarishChandra", "families"] do
    for s in res.(a) do
      s.charNumbers:=List(s.charNumbers,y->Position(tmp,y));
    od;
  od;

  # trying to save some memory
  for a in ["charSymbols","almostCharParams"] do
    if res.(a)=res.charParams then res.(a):=res.charParams; fi;
  od;
  if res.almostCharSymbols=res.almostCharParams then
    res.almostCharSymbols:=res.almostCharParams;
  fi;
  if res.almostTeXCharNames=res.TeXCharNames then
    res.almostTeXCharNames:=res.TeXCharNames;
  fi;
  if Length(res.almostHarishChandra)=Length(res.harishChandra) then
  for i in [1..Length(res.almostHarishChandra)] do
    for a in [ "levi", "cuspidalName", "charNumbers"] do
      if res.almostHarishChandra[i].(a)=res.harishChandra[i].(a) then
	res.almostHarishChandra[i].(a):=res.harishChandra[i].(a);
      fi;
    od;
  od;
  fi;
  res.group:=WF;
  res.operations:=UnipotentCharactersOps;
  return res;
end;

UnipotentCharactersOps.AlmostCharNames:=function(uc,option)
  return List(uc.almostTeXCharNames,x->TeXStrip(x,option));
end;

UnipotentCharacters:=function(W)
  if not IsSpets(W) then W:=Spets(W);fi;
  return AttributeDispatcher("UnipotentCharacters")(W);
end;

UnipotentDegrees:=function(W,q)local uc,v;
  uc:=UnipotentCharacters(W);v:=X(Cyclotomics);
  if not IsBound(uc.degrees) then 
    uc.degrees:=UnipotentCharactersOps.FourierInverse(uc)*FakeDegrees(uc,v);
  fi;
  if IsIdentical(q,v) then return uc.degrees;
  else return List(uc.degrees,x->Value(x,q));
  fi;
end;

CycPolUnipotentDegrees:=function(W)local uc;
  uc:=UnipotentCharacters(W);
  if not IsBound(uc.CycPolDegrees) then 
    uc.CycPolDegrees:=List(UnipotentDegrees(W,X(Cyclotomics)),
      function(p)p:=CycPol(p);p.vname:="q";return p;end);
  fi;
  return uc.CycPolDegrees;
end;

UnipotentCharactersOps.FakeDegrees:=function(uc,q)local fd;
  fd:=[1..Size(uc)]*0*q;
  fd{uc.almostHarishChandra[1].charNumbers}:=FakeDegrees(uc.group,q);
  return fd;
end;

UnipotentCharactersOps.items:=
     ["Name","Degree","FakeDegree","Eigenvalue","Family"];
#  ["n0","Name","Symbol","Degree","FakeDegree","Eigenvalue","Family","Signs"];

UnipotentCharactersOps.possitems:=rec(
  n0:=function(uc,opt)return ["n^0",[1..Size(uc)]];end,
  Degree:=function(uc,opt)
    return ["Deg($\\gamma$)",CycPolUnipotentDegrees(uc.group)];end,
  FakeDegree:=function(uc,opt)local q;
    if IsBound(opt.vname) then q:=opt.vname;else q:="q";fi;
    return ["FakeDegree",
    List(FakeDegrees(uc,X(Cyclotomics)),function(p)
	p:=CycPol(p);p.vname:=q;return p;end)];end,
  Eigenvalue:=function(uc,opt)local q;
    if IsBound(opt.vname) then q:=opt.vname;else q:="q";fi;
    return ["Fr($\\gamma$)", Zip(Eigenvalues(uc,[1..Size(uc)]),
	  uc.operations.qEigen(uc),function(x,y)return x*Mvp(q)^y;end)];end,
  Name:=function(uc,opt)return["$\\gamma$",ShallowCopy(CharNames(uc,opt))];end,
  Family:=function(uc,opt)local n,l,f;
    n:="Label";l:=[];
    for f in uc.families do
      if IsBound(opt.TeX) then l{f.charNumbers}:=f.charLabels;
      else l{f.charNumbers}:=List(f.charLabels,TeXStrip);fi;
    od;
    return [n,l];end,
  Symbol:=function(uc,opt)local l;
    if uc.charSymbols=uc.charParams then l:=false;
    else l:=List(uc.charSymbols,x->Join(List(x,StringSymbol)));
    fi;
    return ["Symbol",l];end,
  Signs:=function(uc,opt)local l,f,n;
    if IsBound(opt.TeX) then n:="$\\varepsilon$";fi;
    l:=[];
    for f in uc.families do
      if IsBound(f.signs) then l{f.charNumbers}:=f.signs;
      else  l{f.charNumbers}:=f.charNumbers*0+1;
      fi;
    od;
    return [n,l];end);

UnipotentCharactersOps.Format:=function(uc,opt)local items,fields,p,i,W,q,res,
  l,n,f,TeX,LaTeX,head,tbl,start,center,fams;
  TeX:=IsBound(opt.TeX); LaTeX:=IsBound(opt.LaTeX);
  W:=uc.group;
  if IsBound(opt.items) then items:=opt.items;
  else items:=UnipotentCharactersOps.items;
  fi;
  opt:=ShallowCopy(opt);
  center:=function(s)
    if LaTeX then return SPrint("\\begin{center}",s,"\\end{center}");
    elif TeX then return SPrint("\\centerline{",s,"}");
    else return s;
    fi;
  end;
  if IsBound(opt.vname) then q:=opt.vname;else q:="q";fi;
  n:=ReflectionName(W,opt);if TeX then n:=SPrint("$",n,"$");fi;
  res:=SPrint(center(SPrint("Unipotent characters for ",n)),"\n");
  head:=[];tbl:=[];
  for n in items do
    if IsBound(UnipotentCharactersOps.possitems.(n)) then
      l:=UnipotentCharactersOps.possitems.(n)(uc,opt);
      if l[2]<>false then
        if not ForAll(l[2],IsString) then l[2]:=List(l[2],x->Format(x,opt));fi;
        Add(tbl,l[2]);
        if LaTeX then l[1]:=SPrint("\\mbox{",l[1],"}");
        elif TeX then l[1]:=SPrint("\\hbox{",l[1],"}");
        else l[1]:=TeXStrip(l[1]);
        fi;
        Add(head,l[1]);
      fi;
    else Error("unknown item:",n,"\tPossibilities are:\n",
      Join(RecFields(UnipotentCharactersOps.possitems)));
    fi;
  od;
  if IsBound(opt.byFamily) then 
    fields:=[2..Length(head)];
    p:=Position(items,"Family");
    fams:=uc.families;
    opt.rowsLabel:=head[1];
    if p<>false then fields:=Difference(fields,[p-1]);fi;
    if opt.byFamily="long" then
      if TeX then Append(res,"\\medskip");fi;
      PrintToString(res,center("Trivial families"),"\n");
      if LaTeX then n:="\\mbox{Family n$^0$}";
      elif TeX then n:="\\hbox{Family n$^0$}";
      else n:="Family n0";
      fi;
      f:=Filtered(uc.families,x->Length(x.charNumbers)=1);
      if ForAll(f,g->not IsBound(g.signs) or ForAll(g.signs,x->x=1)) then
	p:=Position(items,"Signs");
	if p<>false then fields:=Difference(fields,[p-1]);fi;
      fi;
      f:=List(f,g->g.charNumbers[1]);
      l:=List(f,x->PositionProperty(uc.families,g->x in g.charNumbers));
      opt.rowLabels:=tbl[1]{f};
      opt.columnLabels:=Concatenation(head{fields},[n]);
      Append(res,FormatTable(TransposedMat(Concatenation(tbl{fields}{f},[l])),
        opt));
      fams:=Filtered(uc.families,x->Length(x.charNumbers)>1);
    else opt.tbl:=[]; opt.rowLabels:=[]; opt.columnLabels:=[];
      opt.separators:=[];
    fi;
    for f in fams do
      if opt.byFamily="long" then
	if TeX then PrintToString(res,"\\medskip",
	    center(SPrint("Family $n^0.",Position(uc.families,f),
	       "$. of type $",f.name,"$")),"\n");
	  if IsBound(f.comment) then 
	    PrintToString(res,center(SPrint("\\it ",f.comment)),"\n");fi;
	else Append(res,SPrint("\nFamily ",TeXStrip(f.name),"\n"));
	fi;
      fi;
      p:=f.charNumbers[f.special];
      if TeX then n:="*\\hfill ";else n:="*";fi;
      tbl[1][p]:=Concatenation(n,tbl[1][p]);
      if IsBound(f.cospecial) and f.cospecial<>f.special then
	p:=f.charNumbers[f.cospecial];
	if TeX then n:="\\#\\hfill ";else n:="#";fi;
	tbl[1][p]:=Concatenation(n,tbl[1][p]);
      fi;
      fields:=[2..Length(head)];
      if not IsBound(f.signs) or ForAll(f.signs,x->x=1) then
	p:=Position(items,"Signs");
	if p<>false then fields:=Difference(fields,[p-1]);fi;
      fi;
      if opt.byFamily="long" then
	opt.rowLabels:=tbl[1]{f.charNumbers};
	opt.columnLabels:=head{fields};
	Append(res,FormatTable(TransposedMat(tbl{fields}){f.charNumbers},
	   opt));
      else
        Add(opt.separators,Length(opt.tbl));
	Append(opt.rowLabels,tbl[1]{f.charNumbers});
	Append(opt.tbl,TransposedMat(tbl{fields}){f.charNumbers});
      fi;
    od;
    if opt.byFamily<>"long" then
      opt.columnLabels:=head{fields};
      Append(res,FormatTable(opt.tbl,opt));
    fi;
  else 
    if IsBound(opt.chars) then f:=opt.chars;else f:=[1..Size(uc)];fi;
    opt.rowLabels:=tbl[1]{f};
    opt.columnLabels:=head{[2..Length(head)]};
    opt.rowsLabel:=head[1];
    Append(res,FormatTable(TransposedMat(tbl{[2..Length(tbl)]}{f}),opt));
  fi;
  return res;
end;

UnipotentCharactersOps.Display:=function(uc,opt)
  opt:=ShallowCopy(opt);opt.screenColumns:=SizeScreen()[1];
  Print(Format(uc,opt));
end;

UnipotentCharactersOps.LowestPowerGenericDegrees:=function(uc)local uc,ud;
  ud:=CycPolUnipotentDegrees(Group(uc));
  return List(ud{uc.harishChandra[1].charNumbers},x->x.valuation);
end;

UnipotentCharactersOps.HighestPowerGenericDegrees:=function(uc)local uc,ud;
  ud:=CycPolUnipotentDegrees(Group(uc));
  return List(ud{uc.harishChandra[1].charNumbers},Degree);
end;

if not IsBound(CHEVIE.families) then ReadChv("unip/families");fi;
