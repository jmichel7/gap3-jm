# positions(-with-sign) in l where o or -o appears
PositionsSgn:=function(l,o)
  return Concatenation(Positions(l,o),-Positions(l,-o));
end;

# A permutation-with signs has 3 forms:
# list with signs (Ls)
# Element of the hyperoctaedral group (as perm of 1..2n) (Ps)
# Element of the hyperoctaedral group (as monomial matrix with signs)

# Ls to Ps
LsToPs:=function(ls)local n;
  n:=Length(ls);
  ls:=List(ls,function(i)if i<0 then return -i+n;else return i;fi;end);
  ls:=Concatenation(ls,List(ls+n,x->1+(x-1) mod (2*n)));
  return PermList(ls)^-1;
end;

# if N=PsOnMatrices(M,LsToPs(l)) gives [p,v] such that
#    M=OnMatrices(N,p)^DiagonalMat(v)
LsToPermAndSigns:=l->[PermList(List(l,AbsInt)),
  Permuted(List(l,SignInt),PermList(List(l,AbsInt)))];

PermAndSignsToLs:=function(p,s)
  return Permuted(Zip([1..Length(s)],s,function(x,y)return x*y;end),p^-1);
end;

#  Ps to Ls. Must specify n
PsToLs:=function(ps,n)
  return Permuted(Concatenation([1..n],-[1..n]),ps){[1..n]};
end;

# Ls to the i-th power
LsPuiss:=function(ls,i)
  return PsToLs(LsToPs(ls)^i,Length(ls));
end;

# matrix to Ls
MatToLs:=function(m)local n,ps,i,nz;
  n:=Length(m);
  ps:=[];
  for i in [1..n] do
    nz:=Filtered([1..n],x->m[i][x]<>0);
    if Length(nz)<>1 then Print("****** not monomial\n");return false;fi;
    nz:=nz[1];
    if m[i][nz]=1 then ps[i]:=nz;
    elif m[i][nz]=-1 then ps[i]:=-nz;
    else Print("****** not +1,-1-monomial\n");return false;
    fi;
  od;
  return ps;
end;

PsPermuted:=function(l,ps)
  return Permuted(Concatenation(l,-l),ps){[1..Length(l)]};
end;

DisplayLs:=ls->DisplayPs(LsToPs(ls),Length(ls));

DisplayPs:=function(ps,n)local l;
  l:=Filtered(Cycles(ps),x->ForAny(x,x->x<=n));
  l:=List(l,x->List(x,function(y)if y>n then return n-y;else return y;fi;end));
  l:=l{Filtered([1..Length(l)],i->ForAny(l[i],y->not -y in
           Concatenation(l{[1..i-1]})))};
  return String(Concatenation(List(l,x->SPrint("(",Join(x),")"))));
end;

PsOnMatrices:=function(M,p)return PsPermuted(List(M,y->PsPermuted(y,p)),p);end;

PsGroup:=n->Group(Concatenation(List([1..n-1],i->(i,i+1)(n+i,n+i+1)),
List([1..n],i->(i,n+i))),());

# PsMatStab(M [,extra]) find permutations with signs stabilizing M
# (and such that the associated permutations in addition stabilizes extra
#  -- which could be for instance a list of eigenvalues)
PsMatStab:=function(arg)local ind,stab,g,r,l,I,p,ss,dup,M,extra;
  M:=arg[1];
  ss:=x->Set([x,-x]);
  if M<>TransposedMat(M) then Error("M should be symmetric");fi;
  if Length(arg)>1 then extra:=arg[2];else extra:=List(M,x->1);fi;
  ind:=function(I)local iM;
    iM:=List(I,i->[Collected(List(M[i]{I},ss)),M[i][i],extra[i]]);
    return List(Set(iM),x->I{Filtered([1..Length(I)],j->iM[j]=x)});
  end;
  dup:=function(M)
    M:=TransposedMat(Concatenation(M,-M));
    return Concatenation(M,-M);
  end;
  stab:=function(I)local inds,g,p,v,ss;
    inds:=ind(I);
    if Length(inds)>1 then return Concatenation(List(inds,J->stab(J)));
    else
      p:=MappingPermListList([1..2*Length(I)],Concatenation(I,I+Length(M)));
      if Length(I)>5 then 
        InfoChevie("# Block of Length ",Length(I),":",I,"\n");
      fi;
      return [rec(gens:=List(MatStab(PsGroup(Length(I)),
	                     dup(M{I}{I})).generators,x->x^p),ind:=I)];
    fi;
  end;
  l:=stab([1..Length(M)]);
  SortBy(l,x->-Size(Group(x.gens,())));
  g:=Group(());I:=[];
  for r in l do
    if Length(I)>0 then
    p:=MappingPermListList([Length(I)+1..2*Length(I)],
                           [Length(I)+1..2*Length(I)]+Length(r.ind));
    g:=Group(OnTuples(g.generators,p),());
#   Print("prev=",g.generators,"\n");
    fi;
    Append(I,r.ind);
    p:=MappingPermListList(Concatenation(I,I+Length(M)),[1..2*Length(I)]);
#   Print("newind=",Concatenation(r.ind,r.ind+Length(M)),
#         "newgens=",r.gens,"\n");
#   Print("tind=",OnTuples(Concatenation(r.ind,r.ind+Length(M)),p),
#         "tgens=",OnTuples(r.gens,p),"\n");
    g:=Group(Concatenation(g.generators,OnTuples(r.gens,p)),());
#   Print("g=",Size(g),"\c");
    g:=MatStab(g,dup(M{I}{I}));
#   Print(" stab=",Size(g),"\n");
  od;
  return Group(OnTuples(g.generators,p^-1),());
end;

# PsMatMat(M, N[, extra1, extra2]) find p such that PsOnMatrices(M,p)=N
# [and such that Permuted(extra1,p)=extra2]
PsMatMat:=function(arg)
  local ind,l,I,J,r,p,e,g,n,h,trans,psg,tr,q,dup,ss,M,N,extra1,extra2;
  M:=arg[1];N:=arg[2];
  ss:=x->Set([x,-x]);
  if M<>TransposedMat(M) then Error("M should be symmetric");fi;
  if N<>TransposedMat(N) then Error("N should be symmetric");fi;
  if Length(arg)=2 then extra1:=List(M,x->1);extra2:=List(N,x->1);
  else extra1:=arg[3];extra2:=arg[4];
  fi;
  psg:=function(n)local gens,i,j;gens:=[];
    for i in [1..n-1] do 
      for j in [i+1..n] do Add(gens,(i,j)(i+n,j+n));od;
    od;
    Append(gens,List([1..n],i->(i,i+n)));
    return Group(gens,());
  end;
  dup:=function(M)
    M:=TransposedMat(Concatenation(M,-M));
    return Concatenation(M,-M);
  end;
  ind:=function(I,J)local iM,iN,p,n;
    iM:=List(I,i->[Collected(List(M[i]{I},ss)),M[i][i],extra1[i]]);
    iN:=List(J,i->[Collected(List(N[i]{J},ss)),N[i][i],extra2[i]]);
    if Set(iM)<>Set(iN) then InfoChevie("content differs");return false;fi;
    iM:=List(Set(iM),x->I{Filtered([1..Length(I)],j->iM[j]=x)});
    iN:=List(Set(iN),x->J{Filtered([1..Length(J)],j->iN[j]=x)});
    if List(iM,Length)<>List(iN,Length) then 
       InfoChevie("collected content differs");return false;fi;
    if Length(iM)=1 then
      if Length(I)>6 then InfoChevie("large block:",Length(I),"\n");
        p:=DistHelpedRepresentativeOperation(psg(Length(I)),
           M{I}{I},N{J}{J},PsOnMatrices,
	   function(M,N)return Sum(M-N,x->Number(x,y->y<>0));end);
      else p:=RepresentativeOperation(psg(Length(I)),
         M{I}{I},N{J}{J},PsOnMatrices);
      fi;
      if p=false then InfoChevie("could not match block");return false;fi;
      return [[I,J,p]];
    else p:=Zip(iM,iN,ind);
      if false in p then return false;
      else return Concatenation(p);
      fi;
    fi;
  end;
  l:=ind([1..Length(M)],[1..Length(N)]);
  if l=false then return false;fi;
  I:=[];J:=[];g:=Group(());tr:=();
  for r in l do
    n:=Length(r[1]);
    q:=MappingPermListList([1..Length(I)]+Length(I),[1..Length(I)]+Length(I)+n);
    p:=MappingPermListList([1..2*n],
          Concatenation([1..n]+Length(I),[n+1..2*n]+2*Length(I)));
    Append(I,r[1]);Append(J,r[2]);
    Print("#I=",Length(I),"\c");
    if Comm(r[3]^p,tr^q)<>() then Error("noncomm");fi;
    tr:=tr^q*r[3]^p;
    h:=OnTuples(PsMatStab(M{r[1]}{r[1]}).generators,p);
    g:=Group(Concatenation(OnTuples(g.generators,q),h),());
    Print(" #g=",Size(g),"\c");
    e:=RepresentativeOperation(g,M{I}{I},
      PsOnMatrices(N{J}{J},tr^-1),PsOnMatrices);
    if e=false then return false;
    else if Comm(e,tr)<>() then Print("*** tr does not commute to e\n");fi;
         tr:=e*tr;
    fi;
    g:=MatStab(g,dup(M{I}{I}));
    Print(" #stab=",Size(g),"\n");
  od;
  # transporter of a ps from [1..Length(I)] to I
  trans:=I->MappingPermListList([1..2*Length(I)],Concatenation(I,I+Length(M)));
  return trans(I)^-1*tr*trans(J);
end;
