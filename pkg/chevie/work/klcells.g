#########################################################################
##
#A  klcells.g                                                Meinolf Geck
##
#Y  Copyright (C) 2003       Institut Girard Desargues, Universite Lyon 1
##
##  This file provides the function  'KazhdanLusztigCells'  for computing 
##  the  Kazhdan-Lustig pre-order on a finite Coxeter group.  The program 
##  works in the general  multi-parameter  case.  In  the  course  of the 
##  computations,  some of   Lusztig's  conjectures  on  left cells  with  
##  unequal parameters  are checked.  The parameters may be powers of one
##  fixed indeterminate  (in which case they are  simply specified by the
##  exponents) or monomials in two variables.  In  the  latter  case  the  
##  programs   require  an  arithmetic  for  (sparse)  polynomials in two 
##  variables and a monomial order on them. (We provide these programs as 
##  well.) 
##
##  Examples:
##     gap> W := CoxeterGroup( "B", 4 );
##     gap> Q := Indeterminate1Sparse2( ["Q", "q"], LexOrder2 );
##     gap> q := Indeterminate2Sparse2( ["Q", "q"], LexOrder2 );
##     gap> KazhdanLusztigCells( W, [ Q, q, q, q ] );    
##
##  In the case where all parameters are powers of one indeterminate (and
##  where we take the standard monomial order), the function call is:
##
##     gap> KazhdanLusztigCells( W, [ 5, 3, 3, 3 ] ); 
##     gap> KazhdanLusztigCells( W, [ 1, 1, 1, 1 ] ); (equal parameters)    
##
##  The output is a record with several components.  Example:
##     gap> W := CoxeterGroup( "B", 2 );
##     gap> KazhdanLusztigCells( W, [ 4, 1 ] );
##     #I Initializing (Bruhat-Chevalley order etc.) ...
##     #I Computing KL polynomials for elements of length:
##     #I        1 2 3 4 
##     #I Computing left preorder for elements of length:
##     #I        1 2 3 4 
##     #I ---> Checking existence of distinguished involutions and 
##     #I       leading coefficients of KL-pols (should be +1,-1) ...   (OK)
##     #I ---> Checking compatibility of a-invariant and preorder ...   (OK)
##     #I Computing 2-sided preorder: ++++
##     #I ---> Checking Lusztig's property (no relations in 2-cell) ... (OK)
##     #I ---> Checking if w and w^-1 are in the same 2-sided cell ...  (OK)
##     #I #############################################################
##     #I -------> There are 6 left cells and 5 twosided cells <-------
##     #I #############################################################
##     #I The characters carried by the left cells (with the a-invariants):
##     #I (with number in the list of 2-cells and a-invariants):
##     #I ------------------------------------------------------
##     #I    1     0    2.
##     #I ---------------------------------------------------------
##     #I    2     1    .2
##     #I ---------------------------------------------------------
##     #I    3     4    1.1
##     #I ---------------------------------------------------------
##     #I    4     7    11.
##     #I ---------------------------------------------------------
##     #I    5    10    .11
##     #I ---------------------------------------------------------
##     #I (Lines separate two-sided cells)
##     rec( group := CoxeterGroup("B", 2), parameter := [ v^4, v ],
##        words := [ [  ], [ 2 ], [ 1 ], [ 2, 1 ], [ 1, 2 ], 
##                      [ 2, 1, 2 ], [ 1, 2, 1 ], [ 1, 2, 1, 2 ] ],
##        involutions := [ 1, 2, 6, 3, 7, 8 ],  a := [ 0, 1, 4, 4, 7, 10 ],
##        twocells := [ [ 1 ], [ 2 ], [ 6, 3, 5, 4 ], [ 7 ], [ 8 ] ],
##        hasse2 := [ [ 2, 1 ], [ 3, 2 ], [ 4, 3 ], [ 5, 4 ] ],
##        leftcells := [ [ 1 ], [ 2 ], [ 6, 5 ], [ 3, 4 ], [ 7 ], [ 8 ] ],
##        incidence := [ [ true,  true,  true,  true,  true,  true ], 
##                       [ false, true,  true,  false, false, true ], 
##                       [ false, false, true,  false, false, true ],
##                       [ false, false, false, true,  true,  true ],
##                       [ false, false, false, false, true,  true ],
##                       [ false, false, false, false, false, true ] ],
##        characters := [ [ 0, 0, 0, 1, 0 ], [ 0, 0, 0, 0, 1 ], 
##                        [ 0, 1, 0, 0, 0 ], [ 0, 1, 0, 0, 0 ], 
##                        [ 1, 0, 0, 0, 0 ], [ 0, 0, 1, 0, 0 ] ] )
##
##  Comments on the record components:
##    group      : the group under consideration
##    parameter  : the chosen parameters
##    words      : reduced expressions of all elements, ordered by length
##    involutions: the distinguished involutions in the left cells
##                 (specified by their position in the list of words)
##    a          : the a-invariants of the distinguished involutions
##                 (the involutions are ordered by increasing a-invariant)
##    twocells   : the partition of the elements into two-sided cells
##    hasse2     : Hasse relations (=immediate neighbour relations) for 
##                   the LR-preorder
##    leftcells  : the partition of the elements into left cells
##    incidence  : the incidence matrix of the left preorder 
##    characters : the characters of the left cell representations 
##                 (scalar products with the irreducible characters)
##    
##  Some additional information is available  in  record  components of
##  the group 'W' which are created during the computation.
##
##  W.klpol      : the list of all the different Kazhdan-Lusztig 
##                 polynomials that exist for W 
##  W.mpol       : the list of all the different M-polynomials
##  W.klmat      : the matrix of all  Kazhdan-Lusztig and  M-polynomials.
##                 The coefficient  W.klmat[i][j] is bound if  1<=i<=|W|,
##                 1<=j<=i  and the j-th element is less than or equal to
##                 the i-th  element in  the Bruhat-Chevalley order. Then 
##                 W.klmat[i][j]  itself  is a list whose first component
##                 is the  corresponding  Kazhdan-Lusztig element and the
##                 following  components  contain the  M-polynomials (for 
##                 those generators for which they are defined).
##                 
#H  written MG, Nov 2003
##

#########################################################################
##
#A  sparse2.g                                                Meinolf Geck
##
#Y  Copyright (C) 2003       Institut Girard Desargues, Universite Lyon 1
##
##  Here are the functions for dealing with sparse Laurent polynomials in 
##  two variables, where  a  monomial order has been fixed.  
##
##  A polynomial  in  two  variables is  represented as a list of triples 
##  [c1,c2,c3] which  correpond to the various monomials involved in that 
##  polynomial,  and  where  'c1'  is the exponent of the first variable, 
##  'c2'  is  the  exponent  of  the  second  variable, and 'c3'  is  the 
##  coefficient (a rational number). 
##
##  The  monomial order is  specified by a function which takes two lists 
##  as   arguments,  and  which  defines  an ordering  on the  first  two 
##  components  of  these lists. (The   first  two  components  of a list  
##  specify  the exponents in a monomial in two variables.)  For example,
##  given  the  way  in  which   GAP  defines   the  ordering  of  lists,  
##  the lexicographic  order is simply given by  '<'  (see  the  function 
##  'LexOrder2' below).   Another  example  is  provided  by the function 
##  called 'AlmostEqual2' below.  
##
##  Example: create two variables Q and q, with lexicographic order Q>q.
##     gap> LexOrder2 := function(c1,c2) 
##     >      return c1{[1,2]}<c2{[1,2]}; 
##     >    end;
##     gap> Q := Indeterminate1Sparse2( ["Q", "q"], LexOrder2 );
##     gap> q := Indeterminate2Sparse2( ["Q", "q"], LexOrder2 );
##     gap> ( Q - q ^ -1 ) ^ 3;
##     -q^-3+3*Q*q^-2-3*Q^2*q^-1+Q^3              # terms are printed in 
##                                                # increasing order
## 
#H  written MG, Jun 2003
##

#########################################################################
##
##  Some monomial orders
##
LexOrder2:=function(c1,c2) 
  return c1{[1,2]}<c2{[1,2]}; 
end;

AlmostEqual2:=function(c1,c2)
  local s1,s2;
  s1:=c1[1]+c1[2];
  s2:=c2[1]+c2[2];
  if s1<s2 then
    return true;
  elif s1=s2 then
    return c1[1]<c2[1];
  else
    return false;
  fi;
end;

WeightOrder2:=function(c1,c2)
  local s1,s2;
  s1:=4*c1[1]+3*c1[2];
  s2:=4*c2[1]+3*c2[2];
  if s1<s2 then
    return true;
  elif s1=s2 then
    return c1[1]<c2[1];
  else
    return false;
  fi;
end;

###########################################################################
##
#F  IsSparse2( <obj> )  . . . . . . . . . . . . . . . . . . . . type check.
##
IsSparse2:=function(obj)
  return IsRec(obj) and IsBound(obj.isSparse2) and obj.isSparse2=true;
end;

Sparse2Ops:=OperationsRecord("Sparse2Ops");

###########################################################################
##
#V  Sparse2Ops . . . . . . . . . . . . . . . .  operations for polynomials.
##
Sparse2Ops.\=:=function(f,g)
  if IsRec(f) and IsRec(g) then 
    return f.lcoeff=g.lcoeff;
  else
    return false;
  fi;
end;

Sparse2Ops.\+:=function(f,g)
  local res,i,nam,un; 
  if f=0 then return g;
  elif g=0 then return f;
  fi;
  if IsRat(f) then 
    res:=Copy(g.lcoeff);
    if res[Length(res)]{[1,2]}<[0,0] then 
      Add(res,[0,0,f]);
      return rec(isSparse2:=true,names:=g.names,
                 lcoeff:=res,monOrder2:=g.monOrder2,operations:=Sparse2Ops);
    elif [0,0]<res[Length(res)][1] then 
      return rec(isSparse2:=true,names:=g.names,
                 lcoeff:=Concatenation([[0,0,f]],res),
		 monOrder2:=g.monOrder2,operations:=Sparse2Ops);
    else
      Add(res,[0,0,f]);
      nam:=g.names;
    fi;
  elif IsRat(g) then 
    res:=Copy(f.lcoeff);
    if res[Length(res)]{[1,2]}<[0,0] then 
      Add(res,[0,0,g]);
      return rec(isSparse2:=true,names:=f.names,lcoeff:=res,
                    monOrder2:=f.monOrder2,operations:=Sparse2Ops);
    elif [0,0]<res[Length(res)][1] then 
      return rec(isSparse2:=true,names:=f.names,
                 lcoeff:=Concatenation([[0,0,g]],res),
	         monOrder2:=f.monOrder2,operations:=Sparse2Ops);
    else
      Add(res,[0,0,g]);
      nam:=f.names;
    fi;
  else
    if f.lcoeff[Length(f.lcoeff)]{[1,2]}<g.lcoeff[1]{[1,2]} then 
      return rec(isSparse2:=true,names:=f.names,
                 lcoeff:=Concatenation(f.lcoeff,g.lcoeff),
		 monOrder2:=f.monOrder2,operations:=Sparse2Ops);
    elif g.lcoeff[Length(g.lcoeff)]{[1,2]}<f.lcoeff[1]{[1,2]} then 
      return rec(isSparse2:=true,names:=f.names,
                 lcoeff:=Concatenation(g.lcoeff,f.lcoeff),
		 monOrder2:=f.monOrder2,operations:=Sparse2Ops);
    else
      res:=Concatenation(Copy(f.lcoeff),Copy(g.lcoeff));
      nam:=f.names;
    fi;
  fi;
  Sort(res);
  un:=false;
  i:=1;
  while i<Length(res) do 
    if res[i]{[1,2]}=res[i+1]{[1,2]} then 
      res[i][3]:=res[i][3]+res[i+1][3];
      if res[i][3]=0 then Unbind(res[i]); fi;
      Unbind(res[i+1]);
      i:=i+2;
      un:=true;
    else
      i:=i+1;
    fi;
  od;
  if un then res:=Set(res); fi;
  if res=[] then return 0;
  else
    return rec(isSparse2:=true,names:=nam,lcoeff:=res,
                 monOrder2:=f.monOrder2,operations:=Sparse2Ops);
  fi;
end;

Sparse2Ops.\-:=function(f,g)
  return f+(-g);
end;

Sparse2Ops.\*:=function(f,g)
  local res,x,c;
  if f=0 or g=0 then
    return 0;
  elif IsList(f) then return List(f,x->x*g);
  elif IsRat(f) then
    return rec(isSparse2:=true,names:=g.names,
               lcoeff:=List(g.lcoeff,x->[x[1],x[2],x[3]*f]),
               monOrder2:=g.monOrder2,operations:=Sparse2Ops);
  elif IsRat(g) then
    return rec(isSparse2:=true,names:=f.names,
               lcoeff:=List(f.lcoeff,x->[x[1],x[2],x[3]*g]),
               monOrder2:=f.monOrder2,operations:=Sparse2Ops);
  fi;
  if Length(f.lcoeff)<Length(g.lcoeff) then 
    res:=0;
    for c in f.lcoeff do 
      res:=res+rec(isSparse2:=true,names:=f.names,
                   lcoeff:=List(g.lcoeff,x->[x[1]+c[1],x[2]+c[2],x[3]*c[3]]),
                   monOrder2:=f.monOrder2,operations:=Sparse2Ops);
    od;
    return res;
  else
    res:=0;
    for c in g.lcoeff do 
      res:=res+rec(isSparse2:=true,names:=g.names,
                   lcoeff:=List(f.lcoeff,x->[x[1]+c[1],x[2]+c[2],x[3]*c[3]]),
                   monOrder2:=g.monOrder2,operations:=Sparse2Ops);
    od;
    return res;
  fi;
end;
  
Sparse2Ops.\^:=function(f,p)
  local i,g;
  if p=0 then
    return rec(isSparse2:=true,names:=f.names,lcoeff:=[[0,0,1]],
                 monOrder2:=f.monOrder2,operations:=Sparse2Ops);
  fi;
  if Length(f.lcoeff)=1 then 
    return rec(isSparse2:=true,names:=f.names,
          lcoeff:=[[p*f.lcoeff[1][1],p*f.lcoeff[1][2],f.lcoeff[1][3]^p]],
	                   monOrder2:=f.monOrder2,operations:=Sparse2Ops);
  fi;
  if p>0 then
    g:=f;
    for i in [1..p-1] do
      g:=g*f;
    od;
    return g;
  fi;
end;
  
Sparse2Ops.\/:=function(f,g)return f*g^-1;end;

Sparse2Ops.String:=function(r)local i,j,c,n,res;
  Sort(r.lcoeff,r.monOrder2);res:="";
  for n in [1..Length(r.lcoeff)] do
    i:=r.lcoeff[n][1];
    j:=r.lcoeff[n][2];
    c:=r.lcoeff[n][3];
    if i=0 and j=0 then
      if c>0 and n>1 then PrintToString(res,"+"); fi;
      PrintToString(res,c);
    else
      if c<0 then 
        if c=-1 then PrintToString(res,"-"); else PrintToString(res,c,"\*"); fi;
      else 
        if n>1 then PrintToString(res,"+"); fi;
        if c<>1 then PrintToString(res,c,"\*"); fi;
      fi;
      if i<>0 then 
        PrintToString(res,r.names[1]);
	if i<>1 then PrintToString(res,"\^(",i,")"); fi;
      fi;
      if j<>0 then
        if i<>0 then PrintToString(res,"\*"); fi;
	PrintToString(res,r.names[2]);
        if j<>1 then PrintToString(res,"\^(",j,")"); fi;
      fi;
    fi;
  od;
  return res;
end;
  
Sparse2Ops.Print:=function(r)Print(String(r));end;

PosPartSparse2:=function(f)
  local nf;
  if f=0 then
    return 0;
  fi;
  nf:=Filtered(f.lcoeff,cf->f.monOrder2([0,0,1],cf));
  if nf=[] then
    return 0;
  else
    return rec(isSparse2:=true,names:=f.names,lcoeff:=nf,
      monOrder2:=f.monOrder2,operations:=Sparse2Ops);
  fi;
end;

NonNegPartSparse2:=function(f)
  local nf;
  if f=0 then
    return 0;
  fi;
  nf:=Filtered(f.lcoeff,cf->cf{[1,2]}=[0,0] or f.monOrder2([0,0,1],cf));
  if nf=[] then
    return 0;
  else
    return rec(isSparse2:=true,names:=f.names,lcoeff:=nf,
                 monOrder2:=f.monOrder2,operations:=Sparse2Ops);
  fi;
end;

NegPartSparse2:=function(f)
  local nf;
  if f=0 then
    return 0;
  fi;
  nf:=Filtered(f.lcoeff,cf->f.monOrder2(cf,[0,0,1]));
  if nf=[] then
    return 0;
  else
    return rec(isSparse2:=true,names:=f.names,lcoeff:=nf,
                  monOrder2:=f.monOrder2,operations:=Sparse2Ops);
  fi;
end;

ConstPartSparse2:=function(f)
  local p;
  if f=0 then
    return 0;
  fi;
  p:=Position(List(f.lcoeff,x->x{[1,2]}),[0,0]);
  if p<>false then 
    return f.lcoeff[p][3];
  else
    return 0;
  fi;
end;

HighestTermSparse2:=function(f)
  if f=0 then
    return 0;
  elif IsRat(f) then
    return 1;
  else
    Sort(f.lcoeff,f.monOrder2);
    #l:=ShallowCopy(f.lcoeff[Length(f.lcoeff)]);
    #l[3]:=1;
    return rec(isSparse2:=true,names:=f.names,
                lcoeff:=[f.lcoeff[Length(f.lcoeff)]],
                monOrder2:=f.monOrder2,operations:=Sparse2Ops);
  fi;
end;
    
HighestCoeffSparse2:=function(f)
  if f=0 then
    return 0;
  elif IsRat(f) then
    return f;
  else
    Sort(f.lcoeff,f.monOrder2);
    return f.lcoeff[Length(f.lcoeff)][3];
  fi;
end;

LowestTermSparse2:=function(f)
  if f=0 then
    return 0;
  elif IsRat(f) then
    return 1;
  else
    Sort(f.lcoeff,f.monOrder2);
    #l:=ShallowCopy(f.lcoeff[1]);
    #l[3]:=1;
    return rec(isSparse2:=true,names:=f.names,lcoeff:=[f.lcoeff[1]],
                 monOrder2:=f.monOrder2,operations:=Sparse2Ops);
  fi;
end;

BarSparse2:=function(f)
  if f=0 then 
    return 0;
  else
    return rec(isSparse2:=true,names:=f.names,
               lcoeff:=Reversed(List(f.lcoeff,x->[-x[1],-x[2],x[3]])),
               monOrder2:=f.monOrder2,operations:=Sparse2Ops);
  fi;
end;

###########################################################################
##
#F  Indeterminate1Sparse2( <names>, <order> ) . . .  produce indeterminates
#F  Indeterminate2Sparse2( <names>, <order> ) . . .  produce indeterminates
##
Indeterminate1Sparse2:=function(name,order)
  return rec(isSparse2:=true,names:=name,lcoeff:=[[1,0,1]],
                       monOrder2:=order,operations:=Sparse2Ops);
end;

Indeterminate2Sparse2:=function(name,order)
  return rec(isSparse2:=true,names:=name,lcoeff:=[[0,1,1]],
                       monOrder2:=order,operations:=Sparse2Ops);
end;

Q:=Indeterminate1Sparse2(["Q","q"],LexOrder2);
q:=Indeterminate2Sparse2(["Q","q"],LexOrder2);

###########################################################################
##
#F  EvalSparse2  . . . . . . . . . . . . . . . . . . . .  evaluate polynom.
##
EvalSparse2:=function(f,x,y)
  local cf,z;
  if IsRat(f) then 
    return f;
  fi;
  z:=0;
  for cf in f.lcoeff do 
    z:=z+cf[3]*x^cf[1]*y^cf[2];
  od;
  return z;
end;

##########################################################################
##
#F KazhdanLusztigCells( <W>, <para> ) . . . . . . . . KL and M-polynomials  
##
## 'KazhdanLusztigCells' computes the Kazhdan-Lusztig cells (both the left 
## and the two-sided cells) of the finite Coxeter  group  <W> with respect 
## to the parameters <para>. This involves, in particular, the computation 
## of all Kazhdan-Lusztig polynomials and M-polynomials.  A third optional
## argument is permitted, where  one specifies  names for  the irreducible
## characters of <W>.
##
KazhdanLusztigCells:=function(arg)
  local W,MonOrder2,elms,left,right,ld,bruhat,b,w0,pw0,py0,ind,
           delta,i,j,mat,mats,lleft,cmpl,richtig,v,inv,vs,duflo,rest,
	      incmat2,len,h,chn,p,y,w,z,tc,s,sy,sw,m,mp,klmat,klpol,mpol;
  W:=arg[1];
  if ForAll(arg[2],IsInt) then
    MonOrder2:=function(c1,c2) 
      return c1[1]<c2[1]; 
    end; 
    v:=Indeterminate1Sparse2(["v","v1"],MonOrder2);
    vs:=List(arg[2],s->v^s);
  else
    vs:=Copy(arg[2]);
    MonOrder2:=vs[1].monOrder2;
  fi;
  W.vs:=vs;
  if Length(arg)=3 then
    chn:=arg[3];
  else
    chn:=ChevieCharInfo(W).charnames;
  fi;
  # check input
  for i in [1..Length(vs)] do
    for j in [i+1..Length(vs)] do 
      if OrderPerm(W.generators[i]*W.generators[j]) mod 2=1 
           and vs[i]<>vs[j] then
	     Error("#E  parameters are not all right !!!\n");
      fi;
    od;
  od;
  elms:=Flat(List([1..W.N+1],i->CoxeterElementsLength(W,i-1)));
  W.elms:=Copy(elms);
  w0:=LongestCoxeterElement(W);
  len:=List(elms,w->CoxeterLength(W,w));
  inv:=List(elms,w->Position(elms,w^-1));
  Print("#I Initializing (Bruhat-Chevalley order etc.) \c");
  left:=List(elms,w->List(W.generators,s->Position(elms,s*w)));
  right:=List(elms,w->List(W.generators,s->Position(elms,w*s)));
  ld:=[];
  for w in [1..Size(W)] do
    ld[w]:=LeftDescentSet(W,elms[w]);
    Sort(ld[w],function(x,y) 
                 return MonOrder2(vs[x].lcoeff[1],vs[y].lcoeff[1]); 
               end);
  od;
  bruhat:=[];
  for w in [1..Size(W)] do 
    if w>1 and len[w]>len[w-1] then 
      Print(".\c"); 
    fi;
    b:=List([1..w],x->false); 
    IsBlist(b);
    for y in [1..w] do
      if y=1 or y=w then
        b[y]:=true;
      elif y<w then 
        s:=ld[w][1];
        if (left[y][s]<y and left[y][s]<=left[w][s] and 
	                           bruhat[left[w][s]][left[y][s]]) or
            (left[y][s]>y and y<=left[w][s] and bruhat[left[w][s]][y]) then
          b[y]:=true;
        fi;
      fi;
    od;
    bruhat[w]:=ShallowCopy(b);
  od;
  Print("\n");
  ind:=List(elms,w->vs[1]^0*Product(List(CoxeterWord(W,w),s->vs[s])));
  Print("#I Computing KL polynomials for elements of length:\n#I        ");
  klmat:=[]; 
  klpol:=[vs[1]^0]; 
  mpol:=List(W.generators,x->[]);
  for w in [1..Size(W)] do
    if w>1 and len[w]>len[w-1] then 
      Print(len[w]," \c"); 
    fi;
    klmat[w]:=[];
    for y in Reversed([1..w]) do
      if y=w then
        klmat[w][w]:=[1];
      elif bruhat[w][y] then     
        s:=1;
        while s<=W.semisimpleRank and (left[y][s]<y or left[w][s]>w) do
          s:=s+1;
        od;
        if s<=W.semisimpleRank then
          h:=klpol[klmat[w][left[y][s]][1]];
        else
          s:=1;
          while s<=W.semisimpleRank and (right[y][s]<y or right[w][s]>w) do 
            s:=s+1;
          od;
          if s<=W.semisimpleRank then 
            h:=klpol[klmat[w][right[y][s]][1]];
	  else 
	    # now we have to go into the recursion
            s:=ld[w][1]; 
            sw:=left[w][s];
	    sy:=left[y][s];
            h:=ind[w]^-1*klpol[klmat[sw][sy][1]];
	    if y<=sw and bruhat[sw][y] then 
	      h:=h+vs[s]*ind[sw]^-1*klpol[klmat[sw][y][1]];
            fi;
            for z in Reversed([y..sw-1]) do
              if left[z][s]<z and bruhat[z][y] and bruhat[sw][z] then 
                m:=mpol[s][klmat[sw][z][s+1]];
                if m<>0 then 
	          h:=h-m*ind[z]^-1*klpol[klmat[z][y][1]]; 
                fi;
              fi;
            od;
	    h:=ind[w]*h;
	  fi;
	fi;
	# now we have the KL-polynomial for (y,w)
	p:=Position(klpol,h);
	if IsInt(p) then 
	  klmat[w][y]:=[p];
        else
	  Add(klpol,h);
	  klmat[w][y]:=[Length(klpol)];
        fi;
	# compute the corresponding M-polynomials
	for s in [1..W.semisimpleRank] do
	  if left[y][s]<y and left[w][s]>w then 
	    if len[y]+len[w]>W.N then
	      pw0:=Position(elms,elms[w]*w0);
	      py0:=Position(elms,elms[y]*w0);
	      m:=-(-1)^(len[w]-len[y])*mpol[s][klmat[py0][pw0][s+1]];
	    else
              m:=NonNegPartSparse2(vs[s]*ind[y]*ind[w]^-1*
	                                   klpol[klmat[w][y][1]]);
              if Set(arg[2])<>[1] then 
	        for z in Reversed([y+1..w-1]) do
                  if left[z][s]<z and bruhat[z][y] and bruhat[w][z] then
                    mp:=PosPartSparse2(mpol[s][klmat[w][z][s+1]]);
                    if mp<>0 then
                      m:=m-NonNegPartSparse2(mp*ind[y]*ind[z]^-1*
	                                       klpol[klmat[z][y][1]]);
                    fi;
                  fi;
                od;
                if m<>0 then
                  m:=BarSparse2(m)+m-ConstPartSparse2(m);
                fi;
              fi;
	    fi;
	    p:=Position(mpol[s],m);
	    if IsInt(p) then
	      klmat[w][y][s+1]:=p;
	    else
	      Add(mpol[s],m);
	      klmat[w][y][s+1]:=Length(mpol[s]);
	    fi;
	  fi;
        od;
      fi;
    od;
  od;
  delta:=[];
  for w in [1..Size(W)] do
    if ForAll(arg[2],IsInt) and ForAny(arg[2],s->s<=0) then
      b:=ind[w];
    else
      b:=ind[w]*HighestTermSparse2(klpol[klmat[w][1][1]])^-1;
    fi;
    Add(delta,b.lcoeff[1]{[1,2]});
  od;
  Print("\n");
  Print("#I Computing left preorder for elements of length:\n#I        ");
  cmpl:=[];
  for w in [1..Size(W)] do
    if w>1 and len[w]>len[w-1] then 
      Print(len[w]," \c"); 
    fi;
    b:=[w];
    for s in [1..W.semisimpleRank] do
      if left[w][s]>w and not left[w][s] in b then
        Add(b,left[w][s]);
      fi;
    od;
    for z in [1..w-1] do
      if bruhat[w][z] and not z in b then
        if ForAny([1..W.semisimpleRank],s->left[z][s]<z and left[w][s]>w
	    and mpol[s][klmat[w][z][s+1]]<>0*mpol[s][klmat[w][z][s+1]]) then
          Add(b,z);
        fi;
      fi;
    od;
    Add(cmpl,b);
  od;
  # lleft contains the simple left preorder relations
  lleft:=Copy(cmpl);
  # Make transitive closure: cmpl[w] contains all y <= w (left preorder). 
  for z in cmpl do
    for w in z do
      for y in lleft[w] do
        if (not y in z) then
          Add(z,y);
        fi;
      od;
    od;
    Sort(z);
  od;
  Print("\n");
  duflo:=rec(group:=W, parameter:=vs, 
	     words:=List(elms,w->CoxeterWord(W,w)),
             involutions:=[], a:=[],
	     twocells:=[], hasse2:=[],
	     leftcells:=[], incidence:=[]);
  duflo.leftcells:=[];
  Print("#I ---> Checking existence of distinguished involutions and\n");
  Print("#I        leading coefficients of KL-pols (should be +1,-1) ...\c");
  richtig:=true;
  rest:=[1..Size(W)];
  while rest<>[] do
    w:=rest[1];
    for y in rest do
      if MonOrder2(delta[y],delta[w]) then 
        w:=y; 
      fi;
    od;
    b:=Filtered(cmpl[w],x->w in cmpl[x]);
    Sort(b,function(x,y) return MonOrder2(delta[x],delta[y]); end);
    # Check for distinguished involutions
    if Length(b)>1 and (delta[b[1]]=delta[b[2]] or elms[b[1]]^2<>()) then
      Print("\n#  WARNING: cell without distinguished involution !!!\n");
      richtig:=false;
    else
      Add(duflo.involutions,b[1]);
      if not HighestCoeffSparse2(klpol[klmat[b[1]][1][1]]) in [1,-1] then
        Print("\n#  WARNING: leading coefficient not +1 or -1 !!!\n");
        richtig:=false;
      fi;
      if ForAll(arg[2],IsInt) then 
        Add(duflo.a,delta[b[1]][1]);
      else
        Add(duflo.a,delta[b[1]]);
      fi;
    fi;
    rest:=Difference(rest,b);
    Add(duflo.leftcells,b);
  od;
  if richtig then 
    Print("   (OK)\n"); 
  else 
    Print("   (not OK!)\n");
  fi;
  Print("#I ---> Checking compatibility of a-invariant and preorder ...\c");
  richtig:=true;
  for w in duflo.involutions do 
    Add(duflo.incidence,List(duflo.involutions,x->x in cmpl[w]));
  od;
  for y in [1..Length(duflo.involutions)] do
    for z in [1..Length(duflo.involutions)] do 
      if y<z and duflo.incidence[z][y] then 
        richtig:=false;
      elif y<>z and duflo.incidence[y][z] and 
         (delta[duflo.involutions[y]]=delta[duflo.involutions[z]]
            or not MonOrder2(delta[duflo.involutions[y]],
                             delta[duflo.involutions[z]])) then 
        richtig:=false;
      fi;
    od;
  od;      
  if richtig then 
    Print("    (OK)\n"); 
  else 
    Print("    (not OK!)\n");
  fi;
  Print("#I Computing 2-sided preorder: \c");
  # Make transitive closure of cmpl for two-sided preorder 
  for z in duflo.involutions do
    if z>1 and len[z]>len[z-1] then 
      Print("+\c"); 
    fi;
    for w in cmpl[z] do
      for y in lleft[w] do 
        if (not y in cmpl[z]) then
          Add(cmpl[z],y);
        fi;
      od;
      for y in lleft[inv[w]] do
        if (not inv[y] in cmpl[z]) then
          Add(cmpl[z],inv[y]);
        fi;
      od;
    od;
    Sort(cmpl[z]);
  od; 
  for z in Difference([1..Size(W)],duflo.involutions) do
    cmpl[z]:=cmpl[Filtered(duflo.leftcells,x->z in x)[1][1]];
  od;
  Print("\n");
  Print("#I ---> Checking Lusztig's property (no relations in 2-cells) ...\c");
  richtig:=true;
  rest:=Copy(duflo.involutions);
  while rest<>[] do
    w:=rest[1];
    for y in rest do
      if MonOrder2(delta[y],delta[w]) then 
        w:=y; 
      fi;
    od;
    b:=Filtered(cmpl[w],x->w in cmpl[x]);
    Sort(b,function(x,y) return MonOrder2(delta[x],delta[y]); end);
    # Check Lusztig's property (no relations inside a 2-sided cell)
    p:=Filtered([1..Length(duflo.involutions)],x->duflo.involutions[x] in b);
    if Length(Set(delta{duflo.involutions{p}}))>1 then 
      richtig:=false;
    fi;
    for y in p do 
      for z in p do 
        if y<>z and duflo.incidence[y][z] then
	  Print("\n#  WARNING: 2-sided cell with non-trivial relation !!!\n");
	  richtig:=false;
	fi;
      od;
    od;
    rest:=Difference(rest,b);
    Add(duflo.twocells,b);
  od;
  if richtig then 
    Print(" (OK)\n"); 
  else 
    Print(" (not OK!)\n");
  fi;
  Print("#I ---> Checking if w and w^-1 are in the same 2-sided cell ...\c");
  if Set(List(duflo.twocells,z->Set(z)=Set(List(z,y->inv[y]))))=[true] then
    Print("   (OK)\n"); 
  else 
    Print("   (not OK!)\n");
  fi;
  incmat2:=[];
  for w in duflo.twocells do 
    Add(incmat2,List(duflo.twocells,x->x[1] in cmpl[w[1]]));
  od;
  duflo.hasse2:=[];
  for y in [1..Length(incmat2)] do
    for z in [1..Length(incmat2)] do
      if y<>z and incmat2[y][z] then
        if not ForAny([1..Length(incmat2)],x->x<>y and x<>z 
                    and incmat2[y][x] and incmat2[x][z]) then 
          Add(duflo.hasse2,[z,y]);
        fi;
      fi;
    od;
  od;
  Print("#I ######################################################");
  if Length(duflo.leftcells)>=10 then 
    Print("#");
  fi;
  if Length(duflo.leftcells)>=100 then 
    Print("#");
  fi;
  if Length(duflo.twocells)>=10 then 
    Print("#");
  fi;
  if Length(duflo.twocells)>=100 then 
    Print("#");
  fi;
  Print("\n#I ---> There are ",Length(duflo.leftcells));
  Print(" left cells and ",Length(duflo.twocells));
  Print(" two-sided cells <---\n");
  Print("#I ######################################################");
  if Length(duflo.leftcells)>=10 then 
    Print("#");
  fi;
  if Length(duflo.leftcells)>=100 then 
    Print("#");
  fi;
  if Length(duflo.twocells)>=10 then 
    Print("#");
  fi;
  if Length(duflo.twocells)>=100 then 
    Print("#");
  fi;
  Print("\n");
  W.klpol:=klpol;
  W.mpol:=mpol;
  W.klmat:=klmat;
  W.leftcells:=duflo.leftcells;
  duflo.characters:=[];
  for b in duflo.leftcells do 
    mats:=[];
    for s in [1..W.semisimpleRank] do
      mat:=0*IdentityMat(Length(b));
      for j in [1..Length(b)] do
        if left[b[j]][s]<b[j] then 
          mat[j][j]:=-1;
        else
          for i in [1..Length(b)] do 
            if i=j then 
              mat[j][j]:=1;
            elif b[i]=left[b[j]][s] and b[i]>b[j] then
  	      mat[j][i]:=1;
  	    elif b[i]<b[j] and bruhat[b[j]][b[i]] and left[b[i]][s]<b[i] then
  	      mat[j][i]:=-(-1)^(len[b[j]]-len[b[i]])
  	              *EvalSparse2(W.mpol[s][W.klmat[b[j]][b[i]][s+1]],1,1);
            fi;
          od;
        fi;
      od;
      Add(mats,mat);
    od;
    if not CheckHeckeDefiningRelations(Hecke(W,1),mats) then
      Print("\n#  WARNING: representation doesn't satisfy relations !!!\n");
    fi;
    Add(duflo.characters,
         MatScalarProducts(CharTable(W),CharTable(W).irreducibles,
           [CharRepresentationWords(mats,CoxeterConjugacyClasses(W))])[1]);
  od;
  Print("#I The characters carried by the left cells\n");
  Print("#I (with number in the list of 2-cells and a-invariants): \n");
  Print("#I ------------------------------------------------------\n");
  b:=[];
  for tc in [1..Length(duflo.twocells)] do
    for y in duflo.characters{List(Intersection(duflo.twocells[tc],
               duflo.involutions),x->Position(duflo.involutions,x))} do 
      if not y in b then 
        Add(b,y); 
	Print("#I");
        if tc<100 then 
          Print("   ");
	fi;
	if tc<10 then   
          Print(" ");
	fi;
	Print(tc);
        i:=duflo.a[Position(duflo.characters,y)];
        if IsInt(i) then 
	  Print("   ");
          if i<100 then 
            Print(" ");
	  fi;
          if i<10 then
            Print(" ");
      	  fi;
        elif IsList(i) then
	  if i[1]>=0 then
	    Print(" ");
	  fi;
	  if i[2]>=0 then
	    Print(" ");
	  fi;
	  if AbsInt(i[1])<100 then 
	    Print(" ");
	  fi;
	  if AbsInt(i[1])<10 then 
	    Print(" ");
	  fi;
	  if AbsInt(i[2])<100 then 
	    Print(" ");
	  fi;
	  if AbsInt(i[2])<10 then 
	    Print(" ");
	  fi;
        fi;
        Print(i,"    ");
        z:=0;
        for i in [1..Length(y)] do
          if y[i]<>0 then
            if z=1 then 
	      Print(" + ");
	    fi;
            if y[i]<>1 then 
	      Print(y[i],"*");
	    fi;
            Print(chn[i]);
            z:=1;
          fi;
        od;
        Print("\n");
      fi;
    od;
    Print("#I ---------------------------------------------------------\n");
  od;
  Print("#I (Lines separate two-sided cells) \n");
  return Copy(duflo);
end;
