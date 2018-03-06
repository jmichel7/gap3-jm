if not IsBound(MetaAxe) then MeatAxe:=MeatAxe;fi;
#############################################################################
##
##  induce.g                                                Alexander Hulpke
##
## Induction of representations and Block Matrix Operations:
##

#############################################################################
##
#F  IsBlockMat( <mat> )
##
IsBlockMat := function( mat )
    return IsRec( mat ) and IsBound( mat.isBlockMat ) and mat.isBlockMat;
    end;

#############################################################################
##
#V  BlockMatrices  . . . . . . . . . . . . . . domain of all block matrices
#V  BlockMatricesOps . . .  operations for the domain of all block matrices
##
BlockMatricesOps := ShallowCopy(GroupElementsOps);

BlockMatricesOps.\in := function( m, D )
    return IsBlockMat( m );
    end;

BlockMatrices:= rec(
    isDomain   := true,
    name       := "BlockMatrices",
    isFinite   := false,
    size       := "infinity",
    operations := BlockMatricesOps );

#############################################################################
##
#V  BlockMatOps . . . . . . . . . . . . . . . operations for block matrices
##
BlockMatOps := Copy( GroupElementOps );


#############################################################################
##
#F  BlockMat(<m>,<v>,<h>[,<sv>,<sh>][,<field>]) . . . .create a block matrix
##
##  A block matrix is a matrix, that contains entries only at selected
##  positions: [i,j,e] is a subblock e at i,j. e may
##  again be a block matrix (for double induction etc).
##  m is the list of entries (possibly empty), v the vertical, h the
##  horizontal number of blocks. sv and sh are the vertical/horizontal sizes of
##  the small matrices (needed, if m is empty). sv,sh and field are
##  calculated by the routine, if possible.
##  The list m must be sorted and no subblock may be zero, since this would
##  confuse the comparison routine.
##  Warning: Submatrices of block matrices may be IDENTICAL!! (to save memory!)

BlockMat:=function(arg)
local f,m,n,v,h,sv,sh,i;
  m:=arg[1];
  v:=arg[2];
  h:=arg[3];
  if Length(arg)<5 then
    sv:=Length(m[1][3]);
    sh:=Length(m[1][3][1]);
  else
    sv:=arg[4];
    sh:=arg[5];
  fi;
  if Length(m)>0 then
    f:=Field(Flat(m[1][3]));
  else
    f:=Field(arg[Length(arg)]);
  fi;
  # remove zero blocks, that would only confuse!
  n:=[];
  for i in m do
    if IsBlockMat(i[3]) then
      Add(n,i);
    elif i[3]<>0*i[3] then
      Add(n,i);
    fi;
  od;
  Sort(n);
  return rec(mat:=n,
             isBlockMat:=true,
             domain:=BlockMatrices,
             operations:=BlockMatOps,
             field:=f,
             vert:=v,
             horz:=h,
             svert:=sv,
             shorz:=sh );
    end;

#############################################################################
##
#F  BlockMatOps.Transposed( <m> ) . . . . . . . . . . transposed block matrix
##
BlockMatOps.Transposed := function( m )
    return BlockMat( List( m.mat, i-> [ i[2], i[1], Transposed( i[3] ) ] ),
                     m.vert, m.horz, m.svert, m.shorz, m.field ); 
    end;

#############################################################################
##
#F  MatBlockMat(<m>) . . . . . . . . . create matrix from (block) matrix
##
MatBlockMat:=function(m)
local mat,nu,i,j,k,l,li,z,ko;
  if IsMat(m) then
    return m;
  else
    mat:=[];
    nu:=List([1..m.shorz],i->m.field.zero);
    nu:=List([1..m.svert],i->nu);
    for i in [1..m.vert] do
      li:=List([1..m.horz],i->nu);
      for j in List(Filtered(m.mat,j->j[1]=i),i->[i[1],i[2],
                                  MatBlockMat(i[3])]) do
        li[j[2]]:=j[3];
      od;
      for j in [1..m.svert] do
        z:=[];
        for k in [1..m.horz] do
          ko:=(k-1)*m.shorz+1;
          for l in [1..m.shorz] do
            z[ko+l-1]:=li[k][j][l];
          od;
        od;
        Add(mat,z);
      od;
    od;
    return mat;
  fi;
end;

#############################################################################
##
#F  BlockMatMat(<m>,<v>,<h>) . . . . . . . . create block matrix from matrix
##  The resulting block matrix has v resp. h blocks vertical resp.
##  horizontal.
BlockMatMat:=function(m,v,h)
local c,li,i,j,io,jo,sv,sh,ma;
  sv:=Length(m)/v;
  sh:=Length(m[1])/h;
  if not(IsInt(sv) and IsInt(sh)) then
    Error("does not fit");
  fi;
  c:=[];
  for i in [1..v] do
    io:=(i-1)*sv+1;
    li:=Sublist(m,[io..io+sv-1]);
    for j in [1..h] do
      jo:=(j-1)*sh+1;
      ma:=List(li,i->Sublist(i,[jo..jo+sh-1]));
      if ma<>0*ma then
        Add(c,[i,j,ma]);
      fi;
    od;
  od;
  return BlockMat(c,v,h,sv,sh);
end;


#############################################################################
##
#F  Arithmetics for block matrices
##

BlockMatOps.("="):=function(a,b)
  return IsBlockMat(a) and IsBlockMat(b) and
   a.mat=b.mat and a.vert=b.vert and a.horz=b.horz
   and a.svert=b.svert and a.shorz=b.shorz;
end;

BlockMatOps.("<"):=function(a,b)
local i,j,ap,bp,av,bv,n,zai,zbi;
  if not(IsRec(a) and IsBlockMat(a)) then
    return true;
  elif a.field<>b.field then
    return a.field<b.field;
  else
    zai:=[a.vert,a.horz,a.svert,a.shorz];
    zbi:=[b.vert,b.horz,b.svert,b.shorz];
    if zai=zbi then
      n:=List([1..a.shorz],i->a.field.zero);
      for i in [1..a.vert] do
  zai:=Filtered(a.mat,j->j[1]=i);
  zbi:=Filtered(b.mat,j->j[1]=i);

  # get rid of iterated block matrices
  for j in zai do
    j[3]:=MatBlockMat(j[3]);
  od;
  for j in zbi do
    j[3]:=MatBlockMat(j[3]);
  od;
  for j in [1..a.svert] do
    ap:=1;
    bp:=1;
    av:=zai[ap][2];
    if ap<=Length(zai) then
      av:=zai[ap][2];
    else
      av:=a.horz+1;
    fi;
    if bp<=Length(zbi) then
      bv:=zbi[bp][2];
    else
      bv:=a.horz+1;
    fi;
    bv:=zbi[bp][2];
    while av<=a.horz or bv<=a.horz do
      if av=bv then 
        if zai[ap][3][j]<>zbi[bp][3][j] then
    return zai[ap][3][j]<zbi[bp][3][j];
        else
    ap:=ap+1;
    if ap<=Length(zai) then
      av:=zai[ap][2];
    else
      av:=a.horz+1;
    fi;
    bp:=bp+1;
    if bp<=Length(zbi) then
      bv:=zbi[bp][2];
    else
      bv:=a.horz+1;
    fi;
        fi;
      else
        while av<=a.horz and av<bv do
    if zai[ap][3][j]<n then
      return true;
    elif zai[ap][3][j]>n then 
      return false;
    fi;
    ap:=ap+1;
    if ap<=Length(zai) then
      av:=zai[ap][2];
    else
      av:=a.horz+1;
    fi;
        od;
        while bv<=a.horz and bv<av do
    if zbi[bp][3][j]<n then
      return false;
    elif zbi[bp][3][j]>n then 
      return true;
    fi;
    bp:=bp+1;
    if bp<=Length(zbi) then
      bv:=zbi[bp][2];
    else
      bv:=a.horz+1;
    fi;
        od;
      fi;
    od;
  od;

      od;
      # if equal
      return false;
    else 
      return zai<zbi;
    fi;
  fi;
end;

BlockMatOps.("+"):=function(x,y)
local a,b,c,e,ap,bp,av,bv,k,m;
  if y=0 then 
    return x;
  elif x=0 then
    return y;
  elif IsBlockMat(x) and IsBlockMat(y)
     and x.vert=y.vert and x.horz=y.horz
     and x.svert=y.svert and x.shorz=y.shorz then
    e:=Maximum(x.vert,x.horz);
    a:=x.mat;
    b:=y.mat;
    if Length(a)=0 then 
      return y;
    elif Length(b)=0 then
      return x;
    else
      c:=[];
      ap:=1;
      av:=Sublist(a[ap],[1,2]);
      bp:=1;
      bv:=Sublist(b[bp],[1,2]);
      while ap<=Length(a) or bp<=Length(b) do
  if av=bv then
          m:=a[ap][3]+b[bp][3];
    if m<>0*m then
      Add(c,[a[ap][1],a[ap][2],m]);
          fi;
    ap:=ap+1;
    if ap<=Length(a) then
      av:=Sublist(a[ap],[1,2]);
    else
      av:=[e,e];
    fi;         
    bp:=bp+1;
    if bp<=Length(b) then
      bv:=Sublist(b[bp],[1,2]);
    else
      bv:=[e,e];
    fi;         
  else
    while bv<av do
      Add(c,b[bp]);
      bp:=bp+1;
      if bp<=Length(b) then
        bv:=Sublist(b[bp],[1,2]);
      else
        bv:=[e,e];
      fi;         
    od;
    while av<bv do
      Add(c,a[ap]);  
      ap:=ap+1;
      if ap<=Length(a) then
        av:=Sublist(a[ap],[1,2]);
      else
        av:=[e,e];
      fi;         
    od;
  fi;
      od;
      return BlockMat(c,x.vert,x.horz,x.svert,x.shorz);
    fi;
  else
    return MatBlockMat(a)+MatBlockMat(b);
  fi;
end;

BlockMatOps.("-"):=function(a,b)
  return a+((-1)*b);
end;

BlockMatOps.("*"):=function(a,b)
local c,i,j,k,li,m,el;
  if IsInt(b) then
    c:=b;
    b:=a;
    a:=c;
  fi;
  if IsInt(a) and IsBlockMat(b) then
  return BlockMat( List(b.mat,i -> [ i[1],i[2],a*i[3] ] ),
            b.vert,b.horz,b.svert,b.shorz);
  elif IsBlockMat(a) and IsBlockMat(b)
       and a.horz=b.vert
       and a.shorz=b.svert then
    c:=[];
    for i in [1..a.vert] do
      li:=Filtered(a.mat,j->j[1]=i);
      if Length(li)>0 then
        for j in [1..b.horz] do
          m:=0;
          for k in li do
            el:=Filtered(b.mat,i->i[1]=k[2] and i[2]=j);
            if Length(el)>0 then
              m:=m+k[3]*el[1][3];
            fi;
          od;
          if m<>0*m then 
            Add(c,[i,j,m]);
    fi;
        od;
      fi;
    od;
    return BlockMat(c,a.vert,b.horz,a.svert,b.shorz);
  elif IsVector(a) then
    return MatBlockMat(BlockMatMat([a],1,b.vert)*b)[1];
  elif IsVector(b) then
    return b*TransposedBlockMat(a);
  else
    return MatBlockMat(a)*MatBlockMat(b);
  fi;
end;

BlockMatOps.("^"):=function(a,b)
local c,e,i;
  if IsInt(b) and a.vert=a.horz and a.svert=a.shorz then
    if b<0 then
      # Inverse berechnen
      a:=BlockMatMat(MatBlockMat(a)^(-1),a.vert,a.vert);
      b:=-b;
    fi;
    if b>0 then 
      c:=a;
      b:=b-1;
      while b>0 do
        c:=c*a;
        b:=b-1;
      od;
      return c;
    else # b=0
      c:=[];
      if Length(a.mat)=0 then
        e:=IdentityMat(a.svert);
      else
        e:=a.mat[1][3]^0; # to allow for recursive block matrices !
      fi;
      for i in [1..a.vert] do
        Add(c,[i,i,e]);
      od;
      return BlockMat(c,a.vert,a.horz,a.svert,a.shorz);
    fi;
  else
    return b^(-1)*a*b;
  fi;
end;

BlockMatOps.Print:=function(m)
    Print("BlockMat(",m.mat,",",m.vert,",",m.horz,
          ",",m.svert,",",m.shorz,")");
    end;

BlockMatOps.Field := function(m)
    return m.field;
    end;


#############################################################################
## 
#F  BlockMatricesOps.DimensionsMat( BlockMatrices, <m> )
##
BlockMatricesOps.DimensionsMat := function( BlockMatrices, m )
    return [ m.vert * m.svert, m.horz * m.shorz ];
    end;

# induction part

#############################################################################
##
#F  GroupRepresentation(<G>,<gi>)  representation of <G>. 
##  <gi> are the images of the generators of <G>.
##
GroupRepresentation:=function(G,gi)
local hom;
  hom := GroupHomomorphismByImages(G,Group(gi,gi[1]^0),G.generators,gi);
  hom.isHomomorphism:=true;
  hom.isMapping:=true;
  return hom;
end;

InducedRepresentationGenerators:=function(H,r,d,g)
local dim,i,j,m,s;
  m:=[];
  for i in [1..Length(r)] do
    s:=RightCoset(H,r[i]*g);
    if IsBound(s.smallest) then
      j:=Position(r,s.smallest);
    else
      j:=1;
      while not(r[j] in s) do
  j:=j+1;
      od;
    fi;
    Add(m,[i,j,(r[i]*g/r[j])^d]);
  od;
  return BlockMat(m,Length(r),Length(r),
        Length(d.range.1),Length(d.range.1[1]));
end;

#############################################################################
##
#F  InducedRepresentation(<G>,<H>,<d>,<r>)  induced reperesentation
##
##  induce representation <d> from <H> to <G>. <r> is right transversal of
##  G/H
InducedRepresentation:=function(G,H,d,r)
local l;
  if IsBound(RightCoset(H,r[1]).smallest) then
    # Standardize Representatives !
    r:=List(r,i->RightCoset(H,i).smallest);
  fi;
  l:=List(G.generators,i->InducedRepresentationGenerators(H,r,d,i));
  return GroupRepresentation(G,l);
end;

