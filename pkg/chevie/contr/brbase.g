#########################################################################
##
#A  brbase.g                                Meinolf Geck and Sungsoon Kim
##
#Y  Copyright (C) 1996     Equipe des groupes finis, Universite Paris VII
##
##  This  file contains  functions for  computing bi-grassmannians  and the
##  base  of finite Coxeter groups.  You need CHEVIE to  work with them. We
##  also redefine the function 'Bruhat' in the case where the base has been
##  already computed.
##
##  The  base is defined to be the set  of all elements w in W which cannot
##  be  obtained as the supremum  of a subset not  containing w. Denote the
##  base by B. Then every element w in W can be coded by the boolean vector
##  (e_b)_{b  in B} where e_b=1  if b <= w  and e_b=0 otherwise. The Bruhat
##  order  on W is given by the  operation of `boolean OR' on such vectors.
##  The  bases for all finite Coxeter  groups are determined in the article
##  'Bases for the Bruhat--Chevalley order on all finite Coxeter groups'.
## 
##  Example:  If W:=CoxeterGroup(\"type\",rank) is the  record of a Coxeter
##  group  then the command  `BaseBruhat(W)' returns the  base elements and
##  adds  a component `incidence' to the record `W' which contains for each
##  element  in `W.elts' the  associated boolean vector  defined above. The
##  function 'RectricesBruhat(W,w)' returns the maximal elements in B which
##  are smaller than or equal to w. If one is only interested in the set of
##  bi-grassmannians  and the  set of  base elements  (especially for large
##  groups), one should use the commands:
##
##              gap> W := CoxeterGroup( \"F\", 4 );       # for example
##              gap> bg := BiGrassmannians( W );
##              gap> base := FindBaseBruhat( W, bg );
##
#H  written MG, SK, Oct 1996
##
#F BiGrassmannians( <W> ) . . . . . . . . . . . .  Bi-Grassmannians of W 
#F FindBaseBruhat( <W>, <bg> ) . . . . . .  the base of the Bruhat order
#F BaseBruhat( <W> ) . . . . . . . .  the base and the coding of elements
#F RectricesBruhat( <W>, <w> ) . compute the rectrices in the Bruhat order
##
###########################################################################
##
#F BiGrassmannians( <W> ) . . . . . . . . . . . . . . Bi-Grassmannians of W 
##
## 'BiGrassmannians' computes the set of bi-grassmannians in the Coxeter
## group <W>. By definition, an element of <W> is a bi-grassmannian if both
## the left descent set and the right descent set of the element contain
## just one generator of <W>
##
BiGrassmannians:=function(W)local bi,I,s,a,r,i,k,j,ea,l,en,x,n;
  bi:=[];
  for k in [1..W.semisimpleRank] do
    bi[k]:=List([1..W.semisimpleRank],i->[]);
    I:=Difference([1..W.semisimpleRank],[k]);
    n:=1;
    ea:=[()];
    repeat 
      en:=[ ];
      for a in ea  do
        for r in [1..W.semisimpleRank]  do
          if W.rootInclusion[r]/a<=W.parentN then
            x:=a*W.generators[r];
            if not x in en then
              j := 1;
              while j<=Length(I) and W.rootInclusion[I[j]]^x<=W.parentN do
                j:=j+1;
              od;
              if j>=W.semisimpleRank then Add(en,x); n:=n+1; fi;
            fi;
         fi;
       od;
     od;
     for i in en do
       l:=RightDescentSet(W,i);
       if Length(l)=1 then Add(bi[k][W.rootRestriction[l[1]]],i); fi;
     od;
     ea:=en;
    until en=[];
    InfoChevie("#I ",n,": ",List(bi[k],Length),"\n");
  od;
  InfoChevie("#I No of bigrassmannians = ",Sum(Concatenation(bi),Length),"\n");
  return bi;
end;

###########################################################################
##
#F FindBaseBruhat( <W>, <bg> ) . . . . . . . . the base of the Bruhat order
##
## 'FindBaseBruhat' computes the base of the Bruhat order on the Coxeter
## group <W>, where <bg> is the set of bi-grassmannians of <W> as computed
## by the function 'Bigrassmannians'.
## 
FindBaseBruhat:=function(W,bg)
  local leqz,z,z1,mins,l1,l,r,p,s,b,inz,i,j;
  mins:=[]; 
  InfoChevie("#I \c");
  for r in [1..Length(bg)] do
    for s in [1..Length(bg[1])] do
      InfoChevie(" +\c");
      l1:=bg[r][s];
      if l1<>[] then
        inz:=[];
        for i in [1..Length(l1)] do
          inz[i]:=0*[1..Length(l1)];
          for j in [1..i] do
            if Bruhat(W,l1[j],l1[i],CoxeterLength(W,l1[j]),
                                           CoxeterLength(W,l1[i])) then
              inz[i][j]:=1;
            fi;
          od;
        od;
        for z in [1..Length(l1)] do
          leqz:=Filtered([1..z],j->inz[z][j]=1);
          leqz:=Difference([1..Length(l1)],leqz);
          for i in leqz do
            z1:=Filtered(leqz,j->j<=i);
            if Sum(inz[i]{z1})=1 and not bg[r][s][i] in mins then
              Add(mins,bg[r][s][i]);
            fi;
          od;
        od;
        if ForAll(inz,i->i[1]=1) then
          Add(mins,bg[r][s][1]);
        fi;
      fi;
    od;
    InfoChevie("\n#I ");
  od;
  InfoChevie("No of base elements = ",Length(mins),"\n");
  p:=Sortex(List(mins,i->CoxeterLength(W,i)));
  W.base:=Permuted(mins,p);
  return W.base;
end;

###########################################################################
##
#F BaseBruhat( <W> ) . . . . . . . . . . the base and the coding of elements
##
## 'BaseBruhat' returns the base of the Bruhat order on the Coxeter group 
## <W>. Furthermore, the function adds a component 'incidence' to the record
## of <W> which contains the coding of the elements of <W> by boolean vectors
## defined by the base. The ordering is the same as in the component 'elts'.
## 
BaseBruhat:=function(W)
  local l,i,j,mat,z,lb;
  if not IsBound(W.base) then
    W.base:=FindBaseBruhat(W,BiGrassmannians(W));
  fi;
  if IsBound(W.incidence) then
    return W.base;
  fi;
  lb:=List(W.base,i->CoxeterLength(W,i));
  mat:=[];
  InfoChevie("#I Calculating incidence matrix...");
  for l in [1..W.N+1] do
    InfoChevie("+\c");
    mat[l]:=[];
    z:=0;
    for i in CoxeterElements(W,l-1) do
      z:=z+1;
      mat[l][z]:=[];
      for j in [1..Length(W.base)] do
        mat[l][z][j]:=Bruhat(W,W.base[j],i,lb[j],l-1); 
      od;
      IsBlist(mat[l][z]);
    od;
  od;
  InfoChevie("\n");
  W.incidence:=mat;
  return W.base;
end;

###########################################################################
##
#F RectricesBruhat( <W>, <w> ) . . compute the rectrices in the Bruhat order
##
## 'RectricesBruhat' returns the maximal elements in the base which are 
## smaller than or equal to the element <w> in the Coxeter group <W>.
##
RectricesBruhat:=function(W,w)
  local b,z,e,list,mat,inz,pw,sum,i,j;
  if w=[] then
    return [];
  fi;
  mat:=Concatenation(W.incidence);
  e:=Concatenation(W.elts);
  pw:=Position(e,EltWord(W,w));
  b:=List(BaseBruhat(W),i->Position(e,i));
  list:=Filtered([1..Length(mat[1])],i->mat[pw][i]);
  inz:=mat{List(list,i->b[i])}{list};
  sum:=0*[1..Length(inz[1])];
  for i in [1..Length(inz)] do
    for j in [1..Length(inz[1])] do
      if inz[i][j] then
        sum[j]:=sum[j]+1;
      fi;
    od;
  od;
  return List(Filtered([1..Length(sum)],i->sum[i]=1),j->
                                   CoxeterWord(W,e[b[list[j]]]));
end;

BruhatOld:=Bruhat;

Bruhat:=function(arg)local c,py,pw,W,y,w,ly,lw;
  W:=arg[1];y:=arg[2];w:=arg[3];
  if not IsBound(W.incidence) then return BruhatOld(W,y,w);fi;
  if Length(arg)>3 then ly:=arg[4];lw:=arg[5];
  else ly:=W.operations.CoxeterLength(W,y);lw:=W.operations.CoxeterLength(W,w);
  fi;
  py:=Position(W.elts[ly+1],y);# PositionSet works, but Position is faster!
  pw:=Position(W.elts[lw+1],w);
  return IsSubsetBlist(W.incidence[lw+1][pw],W.incidence[ly+1][py]);
end;
