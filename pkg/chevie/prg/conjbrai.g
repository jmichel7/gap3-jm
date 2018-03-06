###############################################################################
##                                                                            #
#A  conjbraid.g             CHEVIE library         Nuno Franco, Jean Michel.  #
##                                                                            #
#Y  Copyright (C) 2002-  The CHEVIE Team                                      #
##                                                                            #
##  This file contains the  Franco & Gonzalez-Meneses and Gebhardt Algorithms #
##  for conjugacy in Garside monoids. Some references:                        #
##                                                                            #
## 1: Franco-Gonzalez Conjugacy problem for braid groups and Garside groups,  #
##    J. Algebra 266 (2003) 112--132                                          #
## 2: Franco-Gonzalez Computation of normalizers in braid groups and Garside  #
##    groups,  Rev. Mat. Iberoamericana 19 (2003) 367--384                    #
## 3: Gonzalez-Gebhardt Solving the conjugacy problem in Garside groups by    #
##    cyclic sliding, J. Symb. Computation 45 (2010) 629--656                 #
###############################################################################

###############################################################################
# Code for finite categories
# They are records with fields:
#  .obj: list of objects
#  .atoms:  same length as .obj, atoms[i]  is a list of records (.map,.tgt)
#  where .map represents a morphism from obj[i] to obj[tgt] (in a conjugacy
#  category obj[i]^map=obj[tgt])

CategoryOps:=OperationsRecord("CategoryOps");

CategoryOps.String:=C->SPrint("category with ",Length(C.obj)," objects and ",
      Sum(C.atoms,Length)," atoms");

CategoryOps.Print:=function(C)Print(String(C));end;

CategoryOps.Display:=function(C,opt)local isolated;
  isolated:=Filtered([1..Length(C.obj)],i->Length(C.atoms[i])=0
    and not ForAny(C.atoms,x->ForAny(x,m->i=m.tgt)));
  if Length(isolated)<>0 then
    Print("Isolated objects: ",Join(C.obj{isolated}),"\n");
  fi;
  ShowMaps(Concatenation(List([1..Length(C.obj)],i->
    List(C.atoms[i],m->rec(src:=C.obj[i],map:=m.map,tgt:=C.obj[m.tgt])))));
end;

##############################################################################
# CategoryByAtoms(b,atoms):   Construct  a  category  from an object b
# and a function atoms(c) giving the atoms of source c as a list of records
# with fields .map=the atom and .tgt the target object.
###############################################################################
CategoryByAtoms:=function(b,atoms)local i,m,p,res;
  res:=rec(atoms:=[[]],obj:=[b]); i:=1;
  while i<=Length(res.obj) do b:=res.obj[i];
    for m in atoms(b) do p:=Position(res.obj,m.tgt);
#	Print(b,"^",m.map,"->",m.tgt,"\n");
      if p=false then 
        Add(res.obj,m.tgt);Add(res.atoms,[]);p:=Length(res.obj);
      fi;
      Add(res.atoms[i],rec(map:=m.map,tgt:=p));
    od;
    i:=i+1;
    if i mod 100=0 then InfoChevie(".\c");fi;
  od;
  res.operations:=CategoryOps;
  return res;
end;

###############################################
# Endomorphisms(C,o) for category C, returns  #
# generators of the endomorphisms of C.obj[o] #
###############################################
CategoryOps.Endomorphisms:=function(C,o)local i,j,paths,maps,gens,nmap,t,foo;
  paths:=[];paths[o]:=[];maps:=[];
  foo:=function()local reached,m;reached:=[o];
    for i in reached do t:=C.atoms[i];
      for j in [1..Length(t)] do m:=t[j];
	if not m.tgt in reached then
	   paths[m.tgt]:=Concatenation(paths[i],[[i,j]]);
	   if i<>o then maps[m.tgt]:=maps[i]*m.map;else maps[m.tgt]:=m.map;fi;
	   Add(reached,m.tgt);
	   if Length(reached)=Length(C.obj) then return; fi;
	fi;
      od;
    od;
  end;
  foo();
  # here paths[p] describes a path to get from obj o to  obj p
  gens:=[];
  for i in [1..Length(C.obj)] do
     t:=C.atoms[i];
     for j in [1..Length(t)] do
       if Concatenation(paths[i],[[i,j]])<>paths[t[j].tgt] then
	 if i=o then nmap:=t[j].map;else nmap:=maps[i]*t[j].map;fi;
	 if t[j].tgt=o then AddSet(gens,nmap);
	 elif nmap<>maps[t[j].tgt] then AddSet(gens,nmap/maps[t[j].tgt]);fi;
       fi;
     od;
  od;
  return gens;
end;
    
#########################################################
# Finds and prints a covering set of maps for a category
ShowMaps:=function(maps)local new,m,p,f,l,l1,l2,i,a,composed;
  maps:=List(maps,x->[x.src,x.map,x.tgt]);
  for f in [1..3] do
    new:=[];
    for m in maps do
      p:=PositionProperty(new,x->x[Length(x)]=m[1]);
      if p=false then 
	p:=PositionProperty(new,x->x[1]=m[Length(m)]);
	if p=false then Add(new,m);
	else new[p]:=Concatenation(m{[1..Length(m)-1]},new[p]);
	fi;
      else new[p]:=Concatenation(new[p],m{[2..Length(m)]});
      fi;
    od;
    maps:=new;
  od;
  for f in maps do
    l1:="";l2:="";
    for i in [1,3..Length(f)-2] do 
      a:=String(f[i]);l:=Maximum(2,Length(String(f[i+1])));
      if Length(l1)+Length(a)+l+1>SizeScreen()[1] then
       Print(l1,"\n",l2,"\n");l1:="";l2:="";fi;
      PrintToString(l2,a,List([1..l],x->'-'),">");
      PrintToString(l1,String("",Length(a)),String(f[i+1],2)," ");
    od;
    Print(l1,"\n",l2,f[Length(f)],"\n");
  od;
end;

############# End of category code ##############################

# faster than the naive code:
# 
# InitialFactor:=function(b)local M;M:=b.monoid;
#   if Length(b.elm)=0 then return M.identity;fi;
#   return M.DeltaAction(b.elm[1],-b.pd);
# end;
#
# PreferredPrefix:=b->b.monoid.LeftGcdSimples(InitialFactor(b),
#                                             InitialFactor(b^-1))[1];
PreferredPrefix:=function(arg)local b,F,M,o; 
  b:=arg[1];M:=b.monoid;
  if Length(b.elm)=0 then return M.identity;fi;
  if Length(arg)=1 then F:=function(x,i)return x;end;else F:=arg[2];fi;
  o:=b.elm[Length(b.elm)];
  return F(M.LeftQuotient(o,M.alpha2(o,F(M.DeltaAction(b.elm[1],-b.pd),1))[1]),-1);
end;

#PositiveSimpleConjugation(b,r[,F])
# conjugation of b by simple r such that inf(b^r)>=inf(b)
# About 2 times faster than (b,r)->b^b.monoid.Elt([r]) for long words
PositiveSimpleConjugation:=function(arg)local b,r,F,M,l;
  b:=arg[1];r:=arg[2];
  if Length(arg)=3 then F:=arg[3];else F:=x->x;fi;
  M:=b.monoid; if r=M.identity then return b;fi;
  l:=M.AddToNormal(ShallowCopy(b.elm),F(r));
  l:=M.PrefixToNormal(M.LeftQuotient(M.DeltaAction(r,b.pd),l[1]),l{[2..Length(l)]});
  return GarsideEltOps.Normalize(M.Elt(l,b.pd));
end;

# RepresentativeSC(b[,F]) returns
# rec(conj:=minimal r such that r^-1*b*F(r) in sliding circuit, 
#     circuit:=sliding circuit)
RepresentativeSC:=function(arg)local b,seen,l,t,r,F;
  b:=arg[1];seen:=[];l:=[[b,b^0]];
  if Length(arg)=1 then F:=function(arg)return arg[1];end;else F:=arg[2];fi;
  repeat AddSet(seen,b); 
    r:=PreferredPrefix(b,F);
    b:=PositiveSimpleConjugation(b,r,F);
    Add(l,[b,l[Length(l)][2]*b.monoid.Elt([r])]);
  until b in seen;
  t:=PositionProperty(l,x->x[1]=b);
  return rec(conj:=l[t][2],circuit:=List(l{[t..Length(l)-1]},x->x[1]));
end;

###############################################################################
# in  the following functions (fields of MinConjugating) a is an element of
# some  Garside monoid M and F is a Garside automorphism of M. x is usually
# a  simple, sometimes an atom  of M. The function returns the minimal (for
# divisibility)  simple  m  such  that  x  divides  m  and m belongs to the
# appropriate  conjugacy  category  --  that  is  m^-1a^F(m) belongs to the
# appropriate conjugacy set.
# The function returns false if no such m exists or a divisor of m would be
# returned for another atom x.
###############################################################################
MinConjugating:=rec();

# Cyc(a,x,F)  Returns atom x if x left-divides a and false otherwise.
MinConjugating.Cyc:=function(a,x,F)local M;
  if a.pd>0 then return x;fi;
  if Length(a.elm)=0 then return false;fi;
  M:=a.monoid;
  if not M.IsLeftDescending(a.elm[1],Position(M.atoms,x)) then return false;fi;
  return x;
end;

# Inf(a,x,F)  minimal m such that x<m and inf(m^-1*a^F(m))>=inf(a). 
# m is simple; see algorithm 2 in Franco-Gonzales 1.
MinConjugating.Inf:=function(a,x,F)local min,s,M;M:=a.monoid;min:=x;
  repeat
    x:=M.DeltaAction(min,a.pd);
    for s in a.elm do x:=M.RightLcmSimples(x,s)[3];od;
    s:=M.RightLcmSimples(x,F(min,1)); min:=F(s[1],-1);
  until s[3]=M.identity;
  return min;
end;

# Pos(a,x,F) minimal m such that x<m and m^-1*a*F(m) is positive. m is simple.
MinConjugating.Pos:=function(a,x,F)
  if a.pd>0 then return x;else return MinConjugating.Inf(a,x,F);fi;
end;

# SS(a,x,F)  Assumes a in SSS(a). Returns minimal m such that x<m and
# m^-1*a*F(m) is in SSS(a). m is simple. See algorithm 5 in Franco-Gonzales 1
MinConjugating.SS:=function(a,x,F)local min,ai;ai:=a^-1;
  repeat min:=x;
    x:=MinConjugating.Inf(a,x,F);
    x:=F(MinConjugating.Inf(ai,F(x),function(x,i)return F(x,-i);end),-1);
  until x=min;
  return min;
end;

# MinConjugating.SC(a,x,F) minimal m such that x<m and m^-1*a*F(m) is in SC(a).
# m is a simple.
MinConjugating.SC:=function(a,x,F)local M,ggF,f,p,l;M:=a.monoid;
# Gebhart-Gonzalez function F for a in SC such that a^x in SSS
  ggF:=function(a,x,F)local f,y,M,r;M:=a.monoid;
    f:=[];y:=a;
    repeat # find the history under sliding of the pair [a,x]
      Add(f,[y,x]);
      r:=PreferredPrefix(y,F);
      x:=M.LeftQuotient(r,
        M.Product(x,PreferredPrefix(PositiveSimpleConjugation(y,x,F),F)));
      if x=M.identity then return [M.identity];fi;
      y:=PositiveSimpleConjugation(y,r,F);
      p:=Position(f,[y,x]);
    until p<>false;
    return List(Filtered(f{[p..Length(f)]},x->x[1]=a),x->x[2]);
  end;
  x:=MinConjugating.SS(a,x,F);
  f:=ggF(a,x,F);
  if f<>[M.identity] then
    p:=PositionProperty(f,s->M.LeftGcdSimples(x,s)[2]=M.identity);
    if p<>false then return f[p];
    else return false;
    fi;
  else p:=PreferredPrefix(a,F);
    if M.LeftGcdSimples(x,p)[2]<>M.identity then return false;fi;
#   l:=Filtered(Concatenation(LeftDivisorsSimple(M,p)),
#            s->s<>x and M.LeftGcdSimples(x,s)[2]=M.identity);
    l:=LeftDivisorsSimple(M,M.LeftQuotient(x,p));
    l:=M.Product(x,Concatenation(l{[2..Length(l)]}));
#   Print("Warning: for b=",a," F=1 & x=",x," divides p=",p," ",Length(l),"\n");
    return First(l,x->x=p or x in ggF(a,x,F));
  fi;
end;

AtomicMaps:=function(arg)local a,F,type,res,i,M,minc,pos,tgt;
  a:=arg[1];arg:=arg{[2..Length(arg)]};
  if Length(arg)>0 and IsFunc(arg[1])then F:=arg[1];arg:=arg{[2..Length(arg)]};
  else F:=function(arg)return arg[1];end;
  fi;
  if Length(arg)=1 then type:=arg[1];else type:="SC";fi;
  M:=a.monoid;res:=[];
  if not IsBound(MinConjugating.(type)) then
    Error(type," should be one of ",RecFields(MinConjugating),"\n");
  fi;
  for i in [1..M.nrAtoms] do
    minc:=MinConjugating.(type)(a,M.atoms[i],F);
    if minc<>false and not 
      ForAny([i+1..M.nrAtoms],k->M.IsLeftDescending(minc,k)) then
      if not type in ["Pos","Cyc"] then # inf does not decrease
           tgt:=PositiveSimpleConjugation(a,minc,F); 
      else tgt:=M.Elt([minc])^-1*a*M.Elt([F(minc)]);
      fi;
      Add(res,rec(map:=M.Elt([minc]),tgt:=tgt)); 
    fi;
  od;
  pos:=x->x.pd>=0;
  return Filtered(res,x->Number(res,y->pos(y.map^-1*x.map))=1);
end;

############################################################################
# RepresentativeConjugation(b,c[,F][,type]) Given b,c in a Garside         #
# group returns false if b and c are not [F-]conjugate, and e if b^e=c.    #
############################################################################
RepresentativeConjugation:=function(arg)local a,class,m,b,bconj,c,cconj,res,e,F,type;
  b:=arg[1];c:=arg[2];arg:=arg{[3..Length(arg)]};
  if not IsBound(b.pd) then Error("Only for Garside monoids");fi;
  if Length(arg)>0 and IsFunc(arg[1])then F:=arg[1];arg:=arg{[2..Length(arg)]};
  else F:=function(arg)return arg[1];end;
  fi;
  if Length(arg)=1 then type:=arg[1];else type:="SC";fi;
  if type in ["SC","SS","USS"] then 
    bconj:=RepresentativeSC(b,F);cconj:=RepresentativeSC(c,F);
#   Print("b->",bconj," c->",cconj,"\n");
    b:=bconj.circuit[1]; c:=cconj.circuit[1]; 
    if b.pd<>c.pd or Length(b.elm)<>Length(c.elm) then return false;fi;
    bconj:=bconj.conj;cconj:=cconj.conj;
  else bconj:=b^0;cconj:=b^0;
  fi;
  if b=c then return bconj*cconj^-1;fi;
  res:=[bconj]; class:=[b];
  for a in class do for m in AtomicMaps(a,F,type) do
    if not m.tgt in class then
      e:=res[Position(class,a)]*m.map;
      if m.tgt=c then return e*cconj^-1;fi;
      Add(class,m.tgt); Add(res,e);
    fi;
  od; od;
  return false;
end;

# ConjugacyCategory(w[,F][,type])
ConjugacyCategory:=function(arg)local b,type,F;
  b:=arg[1];arg:=arg{[2..Length(arg)]};
  if not IsBound(b.pd) then Error("Only for Garside monoids");fi;
  if Length(arg)>0 and IsFunc(arg[1])then F:=arg[1];arg:=arg{[2..Length(arg)]};
  else F:=function(arg)return arg[1];end;
  fi;
  if Length(arg)=1 then type:=arg[1];else type:="SC";fi;
  if type in ["SS","SC","USS"] then b:=RepresentativeSC(b,F).circuit[1];fi;
  return CategoryByAtoms(b,a->AtomicMaps(a,F,type));
end;

ConjugacySet:=function(arg)
  return ApplyFunc(ConjugacyCategory,arg).obj;
end;

#############################################################################
# CentralizerGenerators(a[,F][,type]): for a in a Garside group G, returns  #
# the generating set for the C_G(a) obtained as a generating set for End(a) #
# in the 'type' conjugacy category of a.                                    #
#############################################################################
CentralizerGenerators:=function(arg)local a,b,type,F;
  b:=arg[1];arg:=arg{[2..Length(arg)]};
  if not IsBound(b.pd) then Error("Only for Garside monoids");fi;
  if Length(arg)=0 or not IsFunc(arg[1]) then
    F:=function(arg)return arg[1];end;
  else F:=arg[1];arg:=arg{[2..Length(arg)]};
  fi;
  if Length(arg)=1 then type:=arg[1];else type:="SC";fi;
  if type in ["SS","SC","USS"] then
    b:=RepresentativeSC(b,F);a:=b.conj;b:=b.circuit[1];
  else a:=b^0;
  fi;
  return a*Endomorphisms(ConjugacyCategory(b,F,type),1)*a^-1;
end;

# initially on my laptop with profile the following in B7 takes 14.15 sec
# ConjugacyCategory(B(1,3,4,5,4,3,6,5,1,2,3,2,4,5,4,6,5,4,3,2,1,3,5,6,5,4,
#   3,4,6,5),"USS")
