########################################################################
##
#A murphy.g                    Andrew Mathas   mathas@maths.usyd.edu.au
## 
#Y Copyright (C) July 1998     University of Sydney
## 
## These programs allow calculations with the Murphy basis of the Hecke 
## algebra of type A. They require CHEVIE 3.4 or better (together with 
## some functions in this file which are included from SPECHT).
##
## Multiplication of Murphy basis elements is done using the Garnir 
## tableaux as described in Murphy's paper 'The representations of Hecke 
## algebras of type A', J. Algebra  173 (1995), 97-121. This also lets us 
## convert from the T-basis to the Murphy basis since
##     T_w = M([[1],...,[n]], [[1],...,[n]]) * T_w.
## (We use "M" for the Murphy basis since X is a shorthand for the
## Indeterminate() function.)
##
## As with the T-basis, Murphy basis records have parallel list components
## .elm and .coeff which in our case point to the standard tableaux pairs 
## and the corresponding coefficients. The .elm list is an ordered list of 
## lists of the form [ mu, s, t ] where mu, s and t are all integers and
##   H.partition[ mu ] is the associated partition
##   H.tableaux[ mu][ s and t ] are records describing s and t (these 
## records are described in the function init() below).
##
## Throughout memory considerations are thrown to the wind as we cache
## many of the more horrible expansions as we go along in order to save
## time when we next need them.
##
## There are a number of utility functions below, mainly for converting
## between (standard) tableaux and right coset representatives; however
## the bulk of the functions for working with the Murphy basis are pushed
## into the HeckeAlgebraOps record using CHEVIE's CreateHeckeBasis.
##
########################################################################
RequirePackage("chevie");
RequirePackage("specht");
if not IsBound(HeckeAlgebraOps) then  ## we'd better read hecke.g
  ReadChv("prg/hecke");
fi;

if not IsBound( CHEVIE.MurphySpechtModules ) then
  CHEVIE.MurphySpechtModules:=false;
  Print("\n# Type M:=Basis(H,\"Murphy\"); to access the Murphy basis.\n",
        "\n# Set CHEVIE.MurphySpechtModules:=true to work just inside",
        "\n# the Specht modules. This will make calculations with the",
        "\n# Murphy basis inside Specht modules much faster but will also",
        "\n# mean that T-basis to to Murphy-basis conversions will not work,",
        "\n# even if CHEVIE.MurphySpechtModules is set to false later.\n\n");
fi;

########################################################################

## Given a permuation <x> in Sym(n) return the corresponding CoxeterWord
## (we assume that <x> does indeed belong to W=Sym(n) (certainly it belongs
## in some symmetric group).
CoxeterWordSymmPerm:=function(x) local w, i;
  w:=[];
  while x<>() do
    i:=1; while i^x<(i+1)^x do i:=i+1; od;
    Add(w,i);
    x:=(i,i+1)*x;
  od;
  return w;
end;

## Return the CoxeterGroup permutation corresponding to the symmetric group
## permutation <x>.
CoxeterPermSymmPerm:=function(W,x)
  return EltWord(W,CoxeterWordSymmPerm(x));
end;

## Given two standard tableau <s> and <t> - which are assumed to be of
## the same shape - return the permuation w such that <t>=<s>w.
MappingPermTableaux:=function(s,t)
  return MappingPermListList(Flat(s),Flat(t));
end;

########################################################################

## Finally the Murphy basis code proper:

## We pass a record to CreateHeckeBasis() to add the Murphy basis operations
## to HeckeAlgebraOps.Murphy (we could do this directly ourselves but this
## way we have access to any operations which are added to HeckeEltOps in
## the future).
## The functions passed to CreateHeckeBasis can be accessed via
##   HeckeAlgebraOps.Murphy.<function name>()
CreateHeckeBasis("Murphy", rec(

## Adds various record components to <H> which are used by 
## the Murphy basis functions. init() is called the first time that 
## the Murphy basis is used.
init:=function(H) local mu;
  if IsBound(H.partitions) then return; fi;

  if ReflectionType(CoxeterGroup(H))[1].series<>"A" then
    Error("the Murphy basis is implemented only for type A");
  fi;
  Append(ReadIndent, "  ");
  InfoRead1("#I", ReadIndent, "Initializing Murphy basis...\n");
  
  ##  H.partitions = [ partitions ];
  H.partitions:=[ [1..CoxeterGroup(H).semisimpleRank+1]*0+1 ]; # just 1^n initially

  ## H.tableaux  = [ # partitions ][ rec( ind = index of tableau in H.Tableaux[mu]
  ##                                       mu = partition index in H.partition
  ##                                       wd = associated word in S_n (as a
  ##                                            CoxeterGroup permutation)
  ##                                        s = Xst ( s an index=integer )     
  ##                                   Garnir = [..] will hold the Garnir
  ##                                              expansions          ) ]
  ## Creating a list of all of the standard tableaux is a big overhead for
  ## Specht modules of large dimension so we create these tableaux records
  ## only as needed, storing them t.tableax[mu][...] with a parallel list
  ## h.Tableaux[mu][...] containing the list of the corresponding tableaux.
  ## The bookkeeping to maintain and access this list is done by the
  ## function Tableaux().
  H.tableaux:=[ [rec( ind:=1, mu:=1, wd:=(), Garnir:=[], xst:=[Basis(H,"T")()])] ];
  H.Tableaux:=[ [FirstStandardTableau(H.partitions[1])] ];

  ## T-basis to Murphy-basis conversions; store as calculated, elements
  ## cross-indexed in H.TMBasisElts. Start with identity element.
  H.TwToMurphy:=[HeckeAlgebraOps.Murphy.MakeBasisElt(H,"Murphy",[[[1,1,1]],[H.unit]])];
  H.TMBasisElts:=[()]; # stored as Coxeter words

  InfoRead1("#I",ReadIndent,"...done\n");
  ReadIndent:=ReadIndent{[1..Length(ReadIndent)-2]};
end,

## rewrite Murphy basis in terms of the T-basis
T:=function(M) local m, H;
  H:=Hecke(M);
  m:=Sum([1..Length(M.elm)], 
           m->M.coeff[m]*HeckeAlgebraOps.Murphy.Xst(H, 
                           H.tableaux[M.elm[m][1]][M.elm[m][2]],
                           H.tableaux[M.elm[m][1]][M.elm[m][3]]) );
  ## since HeckeEltOps.\+ normalizes I don't see why it is necessary 
  ## to normalize here but it is...
  CollectCoefficients(m);
  return m;
end,

## rewrite T-basis in terms of the Murphy-basis
Murphy:=function(h) local w, H;
  H:=Hecke(h);
  if CHEVIE.MurphySpechtModules then
    Print("\n# WARNING: because CHEVIE.MurphySpechtModules=true the answer",
          "\n# this function returns will almost certainly be incorrect.\n\n");
  fi;
  return Sum([1..Length(h.elm)], 
               w->h.coeff[w]*HeckeAlgebraOps.Murphy.TwToMurphy(H,h.elm[w]) );
end,

## This function recursively expands T_w into a linear combination of
## Murphy basis elements (using Garnir expansions). Note that we know
## how to write 1 in terms of the Murphy basis. 
## As we go along we cache the expansion of T_w in the list H.TwToMurphy
## which is indexed by the parallel list H.TMBasisElts.
TwToMurphy:=function(H,w) local wpos, W, r;
  wpos:=Position(H.TMBasisElts,w);
  if wpos=false then   
    W:=Group(H);
    r:=FirstLeftDescending(W,w^-1);
    Add( H.TwToMurphy, HeckeAlgebraOps.Murphy.TwToMurphy(H,w*W.(r))*Basis(H,"T")(r) );
    Add(H.TMBasisElts, w);
    wpos:=Length(H.TMBasisElts);
    if w<>w^-1 then
      Add(H.TMBasisElts,w^-1);
      Add(H.TwToMurphy, AlphaInvolution(H.TwToMurphy[wpos]));
    fi;
  fi;
  return ShallowCopy(H.TwToMurphy[wpos]);
end,

# The command Basis(H,"Murphy")(<s>,<t>) gets redirected to this
# function which will return a record which corresponds to the Murphy
# basis element M_{<s>,<t>}. Note that we use HeckeElt() to
# create the actual basis element.
MakeBasisElt:=function(H,basis,tabs)local s, t, mu;
  if Length(tabs)=1 then
    if Basis(tabs[1])="Murphy" then return tabs[1];
    else return HeckeAlgebraOps.Murphy.Murphy(Basis(H,"T")(tabs[1]));
    fi;
  elif Length(tabs)<>2 then Error("usage, M(<s>, <t>)"); 
  elif tabs[2]=[] or not IsList(tabs[2][1]) then
    ## assume tabs[1]=<elm> and tabs[2]=<coeffs>
    return HeckeElt(H,"Murphy",tabs[1],tabs[2]);
  fi;
  s:=tabs[1]; t:=tabs[2];
  if not ( IsStandardTableau(s) and IsStandardTableau(t) ) then
    Error("<s> and <t> must standard tableaux"); 
  fi;
  mu:=List(s,Length);
  if mu<>List(t,Length) then
    Error("<s> and <t> must be standard tableaux of the same shape"); 
  fi;
  mu:=HeckeAlgebraOps.Murphy.Partition(H,mu);
  if IsList(s) and IsList(t) then
    s:=HeckeAlgebraOps.Murphy.Tableau(H,mu,s);
    t:=HeckeAlgebraOps.Murphy.Tableau(H,mu,t);
  fi;
  return HeckeElt(H,"Murphy", [[s.mu, s.ind, t.ind]], [H.unit]);
end,

#F returns a `tight' string for the tableau <t>
StringTableau:=function(H,t)
  if CHEVIE.PrintHecke="GAP" then return StringTableau(H.Tableaux[t.mu][t.ind]);
  else return VeryTightTableau(H.Tableaux[t.mu][t.ind]);
  fi;
end,

## printing is done by calling String()
String:=function(h) local H,i,j,coeff,s,needsbrackets,res;
  H:=Hecke(h);
  if h.elm=[] then return "0"; fi;
  res:="";
  for i in [1..Length(h.elm)] do
    coeff:=String(h.coeff[i]);
    needsbrackets:=false;
    for j in [2..Length(coeff)] do
      if (coeff[j]='+' or coeff[j]='-') and coeff[j-1]<>'^' then
	      needsbrackets:=true;fi;
    od;
    if needsbrackets then coeff:=Concatenation("(",coeff,")");
    elif coeff="1" then coeff:="";
    elif coeff="-1" then coeff:="-";
    fi;
    if CHEVIE.PrintHecke<>"GAP" and i>1 then Append(res,"\n"); fi;
    if Position(coeff,'-')<>1 and i>1 then Append(res,"+");fi;
    Append(res,String(coeff));
    if CHEVIE.PrintHecke="GAP" then Append(res,"*");fi;
    if CHEVIE.MurphySpechtModules then
      res:=Concatenation(res,"S(", 
             HeckeAlgebraOps.Murphy.StringTableau(H,H.tableaux[h.elm[i][1]][h.elm[i][3]] ), ")");
     else
      res:=Concatenation(res,"M(", 
             HeckeAlgebraOps.Murphy.StringTableau(H,H.tableaux[h.elm[i][1]][h.elm[i][2]] ), ", ",  
             HeckeAlgebraOps.Murphy.StringTableau(H,H.tableaux[h.elm[i][1]][h.elm[i][3]] ), ")");
     fi;
  od;
  return String(res);
end,

## returns a latex string for h
TeX:=function(h) local H, res, coeff, needsbrackets, i, j;
  if h.elm=[] then return "0"; fi;
  H:=Hecke(h);
  res:="";
  for i in [1..Length(h.elm)] do
    coeff:=String(CycPol(h.coeff[i]));
    needsbrackets:=false;
    for j in [2..Length(coeff)] do
      if (coeff[j]='+' or coeff[j]='-') and coeff[j-1]<>'^' then
	      needsbrackets:=true;
      fi;
    od;
    if needsbrackets then coeff:=Concatenation("(",coeff,")");
    elif coeff="1" then coeff:="";
    elif coeff="-1" then coeff:="-";
    fi;
    if Position(coeff,'-')<>1 and i>1 then PrintToString(res,"+");fi;
    PrintToString(res,coeff);
    if CHEVIE.MurphySpechtModules then
      PrintToString(res, TeXTableau(H.Tableaux[h.elm[i][1]][h.elm[i][3]]),"\n");
     else
      PrintToString(res,"M(", 
             TeXTableau(H.Tableaux[h.elm[i][1]][h.elm[i][2]]), ", ",  
             TeXTableau(H.Tableaux[h.elm[i][1]][h.elm[i][3]]), ")\n");
     fi;
  od;
  return res;
end,

## returns a latex string for h
Display:=function(h,junk) local H, res, coeff, needsbrackets, i, j;
  if h.elm=[] then return "0"; fi;
  H:=Hecke(h);
  res:="";
  for i in [1..Length(h.elm)] do
    coeff:=String(CycPol(h.coeff[i]));
    needsbrackets:=false;
    for j in [2..Length(coeff)] do
      if (coeff[j]='+' or coeff[j]='-') and coeff[j-1]<>'^' then
	      needsbrackets:=true;
      fi;
    od;
    if needsbrackets then coeff:=Concatenation("(",coeff,")");
    elif coeff="1" then coeff:="";
    elif coeff="-1" then coeff:="-";
    fi;
    if Position(coeff,'-')<>1 and i>1 then PrintToString(res,"+");fi;
    PrintToString(res,coeff);
    if CHEVIE.MurphySpechtModules then
      PrintToString(res," S(", 
             HeckeAlgebraOps.Murphy.StringTableau(H,H.tableaux[h.elm[i][1]][h.elm[i][3]]), ")");
    else
      PrintToString(res," M(", 
             HeckeAlgebraOps.Murphy.StringTableau(H,H.tableaux[h.elm[i][1]][h.elm[i][2]]), ", ",  
             HeckeAlgebraOps.Murphy.StringTableau(H,H.tableaux[h.elm[i][1]][h.elm[i][3]]), ")");
    fi;
    PrintToString(res,"\n");
  od;
  Print(res);
end,


## given a tableau <t> return the basis element x_mu*T_<t> with respect to
## the T-basis
Xt:=function(H, t) local J, WJ, xt, row;
  if not IsBound(H.tableaux[t.mu][t.ind].xst[1]) then
    if t.ind=1 then 
      J:=[];  # -> generating set for Sym(<t.mu>)
      for row in H.Tableaux[t.mu][1] do
        Append(J, row{[1..Length(row)-1]});
      od;
      WJ:=List(CoxeterWords(ReflectionSubgroup(CoxeterGroup(H),J)),
                                i->EltWord(CoxeterGroup(H),i));
      xt:=HeckeElt(H,"T",WJ,[1..Length(WJ)]*0+H.unit);
    else
      xt:=ShallowCopy(HeckeAlgebraOps.Murphy.Xt(H, H.tableaux[t.mu][1]));
      xt.elm:=xt.elm*HeckeAlgebraOps.Murphy.CoxeterPermTableau(H,t); 
    fi;
    CollectCoefficients(xt);
    H.tableaux[t.mu][t.ind].xst[1]:=xt;
  fi;
  return H.tableaux[t.mu][t.ind].xst[1];
end,

## return the basis element <t>.(<s>.ind) = T_<s>^* x_mu T_<t>  with
## respect to the T-basis.
Xst:=function(H,s,t) local xst;
  if not IsBound(H.tableaux[t.mu][t.ind].xst[s.ind]) then
    if s.ind=1 then return ShallowCopy(HeckeAlgebraOps.Murphy.Xt(H,t));
    elif t.ind=1 then 
      xst:=AlphaInvolution(HeckeAlgebraOps.Murphy.Xt(H,s));
    elif IsBound(s.(t.ind)) then 
      xst:=AlphaInvolution(s.(t.ind));
    else
      xst:=Basis(H,"T")(HeckeAlgebraOps.Murphy.CoxeterPermTableau(H,s)^-1)
                 * HeckeAlgebraOps.Murphy.Xt(H,t);
    fi;
    CollectCoefficients(xst);
    H.tableaux[t.mu][t.ind].xst[s.ind]:=xst;
  fi;
  return ShallowCopy(H.tableaux[t.mu][t.ind].xst[s.ind]);
end,

## Given a partition <mu> , which is assumed to be a partition of n, return the 
## index of <mu> in the list H.partitions, or add <mu> to this list if it is not
## already there AND intilize  the list H.tableaux[mu] and H.Tableaux[...].
Partition:=function(H,mu) local t, s;
  t:=Position(H.partitions, mu);
  if t=false then
    t:=Length(H.partitions)+1;
    H.partitions[t]:=Copy(mu);
    H.tableaux[t]:=[ rec( ind:=1, mu:=t, wd:=(), Garnir:=[], xst:=[] ) ];
    H.Tableaux[t]:=[ FirstStandardTableau(H.partitions[t]) ];
  fi;
  return t;
end, 

## given a tableau <tab> , which is assumed to be standard, return the corresponding 
## element of <H>.tableau (a "CoxeterGroup tableau"). Here, <mu> is the index of the 
## shape of <tab> in H.partitions[].
Tableau:=function(H,mu,tab) local t, s;
  t:=Position(H.Tableaux[mu],tab);
  if t=false then
    t:=Length(H.tableaux[mu])+1;
    H.tableaux[mu][t]:=rec(mu:=mu, ind:=t, Garnir:=[], xst:=[]);
    H.Tableaux[mu][t]:=Copy(tab);
  fi;
  return H.tableaux[mu][t];
end, 

## returns the position of the integer <i> in the tableau record <t>.
PositionInTableau:=function(H,t,i) local row, col;
  for row in [1..Length(H.Tableaux[t.mu][t.ind])] do
    col:=Position(H.Tableaux[t.mu][t.ind][row], i);
    if col<>false then return rec(row:=row,col:=col); fi;
  od;
  Error(i, " is not contained in the tableau <t>");
end,

## HeckeEltOps gives control of \* to basisOps.\* if it exists
## As we know how to multiply the Murphy basis we seize the reins here
\*:=function(m, h) 
  local H, W, q, mh, w, mw, mr, t, nodeR, nodeS, tr, s, a, r, i;
  H:=Hecke(h);
  if not IsRec(m) or not IsBound(m.hecke) or not IsBound(m.elm) then
  # assume m is a scalar by which to multiply h
    return HeckeElt(H,h.basis,h.elm,h.coeff*(m*H.unit));
  fi;
  if not IsIdentical(H,Hecke(m)) then
    Error("not elements of the same algebra");
  fi;
  if m.basis<>"Murphy" then  # it must be true that h.basis="Murphy"
    return AlphaInvolution( AlphaInvolution(h)*AlphaInvolution(m) );
  fi;
  if h.basis<>"T" then h:=Basis(H,"T")(h); fi; #(and m.basis="Murphy")
  W:=CoxeterGroup(H);
  q:=H.parameter[1][1];
  # mh=return value, initially zero with respect to the Murphy basis
  mh:=HeckeElt(H,"Murphy",[],[]);  
  for a in [1..Length(h.elm)] do 
    w:=CoxeterWord(W,h.elm[a]); # multiply one simple reflecton at a time
    mw:=ShallowCopy(m);
    for r in w do
      mr:=HeckeElt(H,"Murphy",[],[]);
      for i in [1..Length(mw.elm)] do
        t:=H.tableaux[ mw.elm[i][1] ][ mw.elm[i][3] ];
        ## we are interested in the two nodes, <nodeR> and <nodeS> 
        ## which are swapped by the transposition r=(r,r+1). Thus,
        ## these are the nodeRs such that t<nodeR>=r and t<nodeS>=r+1.
        nodeR:=HeckeAlgebraOps.Murphy.PositionInTableau(H,t, r);
        nodeS:=HeckeAlgebraOps.Murphy.PositionInTableau(H,t, r+1);
        if nodeR.row=nodeS.row then       ## same row
          Add(mr.coeff, q*mw.coeff[i]);
          Add(mr.elm, mw.elm[i]);
        elif nodeR.col<>nodeS.col then    ## different row and column => t*r is still standard
          tr:=Copy(H.Tableaux[t.mu][t.ind]);
          tr[nodeR.row][nodeR.col]:=r+1;
          tr[nodeS.row][nodeS.col]:=r;
          tr:=HeckeAlgebraOps.Murphy.Tableau(H,t.mu,tr);
          if nodeR.row<nodeS.row then     ## up in the Bruhat order
            Add(mr.coeff, mw.coeff[i]); 
            Add(mr.elm, [tr.mu, mw.elm[i][2], tr.ind]);
          else
            Append(mr.coeff, [q, q-1]*mw.coeff[i] );
            Append(mr.elm, [ [tr.mu, mw.elm[i][2], tr.ind], mw.elm[i] ] );
          fi;
        else  ## The hard part: here in the tableau t, r+1 occupies the nodeR 
              ## below r; so nodeS.row=nodeR.row+1 and nodeS.col=nodeR.col and
              ## interchanging them gives a non-standard tableau.
          s:=H.tableaux[ mw.elm[i][1] ][ mw.elm[i][2] ];
          mr:=mr+mw.coeff[i]*HeckeAlgebraOps.Murphy.GarnirExpansion(H,nodeR,s,t);
        fi;
      od;
      mw:=mr;CollectCoefficients(mw);
    od;
    mh:=mh+h.coeff[a]*mw;
  od;
  CollectCoefficients(mh);
  return mh;
end,

## The real work: here <H> is a Hecke algebra record
##    <node>=(r,c) is the coordinate where the tableau <t> becomes
## non-standard; ie. if <t>(r,c)=x then <t>(r+1,c)=x+1 and we want
## to expand <t>*(x,x+1). We are actually expanding 
##         M_{<s>,<t>} * T_{(x,x+1)}
## whexe we know that <t>*(x,x+1) is not standard. Once we know how
## to expxess M_{1,<t>}*T_x in the Murphy basis we store this as
## t.Garnir[x] for future reference.
GarnirExpansion:=function(H,node,s,t)
  local rt, W, gtab, a, b, g, w, rg, tab, coeff, tnu, d, J, K, h;
  rt:=H.Tableaux[t.mu][t.ind][node.row][node.col]; # acting on t by (rt,rt+1)
  if not IsBound( t.Garnir[rt] ) then 
    W:=CoxeterGroup(H);
    ## a typical situation here is
    ##         1 2 7
    ##   <t> = 3 5 8     with <node>=[2,2];
    ##         4 6 9
    ## thus, we want to expand the tableau
    ##         1 2 7 
    ##    t' = 3 6 8 = <t>*(5,6) [note that 5 occupies position [2,2] in <t>]
    ##         4 5 9
    ## into a linear combination of standard tableaux. To do this we first
    ## pretend that we started with
    ##         1 2 3
    ##    g  = 4 6 8
    ##         5 7 9
    ## (that is we put the numbers in order upto, but not including [node.row,node.col]
    ## and then enter them in order starting from the next row down, filling
    ## up the nodes around <node> and then continuing on. ie. what almost Murphy
    ## called a Garnir tableau).

    gtab:=Copy(H.Tableaux[t.mu][1]);
    a:=gtab[node.row][node.col];   ## first number being moved; above a=5 and b=8
    b:=gtab[node.row+1][node.col]; ## last number being moved
    gtab[node.row]{[node.col+1..Length(gtab[node.row])]}:=[a+node.col+1..b];
    gtab[node.row+1]{[1..node.col-1]}:=[a..a+node.col-2];
    gtab[node.row][node.col]:=a+node.col-1;
    gtab[node.row+1][node.col]:=a+node.col;
    g:=HeckeAlgebraOps.Murphy.Tableau(H,t.mu,gtab);

    ## w is the permutation such that t=g*w => T_t = T_g*T_w
    w:=CoxeterWordSymmPerm(MappingPermTableaux(gtab,H.Tableaux[t.mu][t.ind]));

    # we will compute the Garnir expansion for the Garnir tableaux first,
    # cache this and then use it to cmpute the expansion of m_st.
    rg:=H.Tableaux[g.mu][g.ind][node.row][node.col]; # acging on g by (rg,rg+1)
    if not IsBound(H.tableaux[g.mu][g.ind].Garnir[rg]) then
      ## first note that, by an astute look at right coset sums,
      ##      1 2 3   1 2 3   1 2 3   1 2 3   1 2 3         1 2 3       1 2 3
      ## (*)  4 5 6 + 4 5 7 + 4 5 8 + 4 6 7 + 4 6 8 + ... + 4 7 8 = h * 4
      ##      7 8 9   6 8 9   6 7 9   5 8 9   5 7 9         5 6 9       5 6 7 8
      ##                                                                9
      ## Because of our choice of g all of the LH tableaux are standard, except
      ## the last  and the term on the RHS. We'll worry about the RHS later.
      ## First we spin out the tableaux on the left hand side.
      tab:=[]; coeff:=[];
      for J in Combinations([a..b],node.col) do
        if J<>[a..a+node.col-1] then
          gtab[node.row+1]{[1..node.col]}:=J;
          gtab[node.row]{[node.col..Length(gtab[node.row])]}:=Difference([a..b],J);
          ## note that we set <s>=t^mu below; this is because we later
          ## have to multiply by T_s^*
          Add(tab,[g.mu,1,HeckeAlgebraOps.Murphy.Tableau(H,g.mu,gtab).ind]); 
          Add(coeff, -1);
        fi;
      od;
      H.tableaux[g.mu][g.ind].Garnir[rg]:=HeckeElt(H,"Murphy",tab,coeff);
      CollectCoefficients(H.tableaux[g.mu][g.ind].Garnir[rg]);

      ## Next, if CHEVIE.MurphySpechtModules=false (in which case we 
      ## just work in the Specht module), we look after the right hand term
      ## in (*) above. In general it won't correspond to a partition but we
      ## can find a partition <nu> and a <d> in <W> such that 
      ## T_d<RHS>=<x_nu>T_<d>
      if CHEVIE.MurphySpechtModules=false then
        ## First we set <tab> equal to the the new tableau on the right in (*);
        ## above, <tab>:=[[1,2,3],[4],[5,6,7,8],[9]].
        tab:=gtab{[1..node.row-1]};
        if node.col>1 then Add(tab,gtab[node.row]{[1..node.col-1]}); fi;
        Add(tab, [a..b]);
        if node.col<Length(gtab[node.row+1]) then 
          Add(tab, gtab[node.row+1]{[node.col+1..Length(gtab[node.row+1])]}); 
        fi;
        if IsBound(gtab[node.row+2]) then Append(tab, gtab{[node.row+2..Length(gtab)]}); fi;
        ## Now we reorder <tab> so that the diagram has the shape of a
        ## partition. Our <tab> above becomes [[5,6,7,8],[1,2,3],[4],[9]].
        SortBy(tab, x->[-Length(x),x[1]]);
        ## which gives us the (shape of the) new tableau
        tnu:=H.tableaux[ HeckeAlgebraOps.Murphy.Partition(H,List(tab,Length)) ][1];

        ## and finally we have <d>. The point is that tab = T_d^-1*tnu*T_d.
        d:=CoxeterPermSymmPerm(W, MappingPermTableaux(H.Tableaux[tnu.mu][tnu.ind],tab));

        ## <tab> is now under control, but we still need to compute <h>
        ## from (*). The point here is that we are essentially writing
        ##             $H.I = \bigcup H.d = \bigcup d' I$
        ## where H and I are two sugroups and d and d' run over coset
        ## representatives of H and I.
        gtab:=H.Tableaux[g.mu][1];
        J:=gtab[node.row]{[1..Length(gtab[node.row])-1]};
        Append(J, gtab[node.row+1]{[1..Length(gtab[node.row+1])-1]});
        K:=Difference(J,[a-1,b]);
        h:=List(ReducedRightCosetRepresentatives(
                             ReflectionSubgroup(W,W.rootInclusion{J}),
                             ReflectionSubgroup(W,W.rootInclusion{K}) ),
                 coeff->coeff^-1);
        h:=HeckeElt(H,"T",h,[1..Length(h)]*0+1);

        ## the multiplication below is quite costly as it is recursive; 
        ## but it is only done once as we store the result in g.Garnir.
        H.tableaux[g.mu][g.ind].Garnir[rg] 
           := H.tableaux[g.mu][g.ind].Garnir[rg] 
                + h * Basis(H,"T")(d)^-1 
                    * Basis(H,"Murphy")(H.Tableaux[tnu.mu][1],H.Tableaux[tnu.mu][1])
                    * Basis(H,"T")(d);
      fi;  ## end of if CHEVIE.MurphySpechtModules=false
    fi;    ## end of Garnir expansion for the Garnir tableau

    ## Next we worry about the element <w> above (remember t=g*w).
    ## This multiplication is usually recursive.
    if w<>[] then
      H.tableaux[t.mu][t.ind].Garnir[rt]:=H.tableaux[g.mu][g.ind].Garnir[rg]*Basis(H,"T")(w);
    fi;
  fi;  ## end of if IsBound( t.Garnir[t.tab[node.row][node.col]] )

  ## Finally we have to put <s> back into the equation. If we are working
  ## in just the Specht module <s> is almost irrelevant; but in general it 
  ## affects tnu in strange ways (hence it might be better to cache the
  ## full expansion rather than just the right hand side). 
  if s.ind=1 then return ShallowCopy(H.tableaux[t.mu][t.ind].Garnir[rt]);
  else
    return AlphaInvolution(AlphaInvolution(H.tableaux[t.mu][t.ind].Garnir[rt])
           * Basis(H,"T")(HeckeAlgebraOps.Murphy.CoxeterPermTableau(H,s)) );
  fi;
end,

## This is the anti-isomorphism of the Hecke algebra given by T_i -> T_i;
## on the Murphy basis we have M_{s,t} -> M_{t,s} (also called * by many)
AlphaInvolution:=function(h) local hstar, t;
  hstar:=HeckeElt(Hecke(h),"Murphy",[],[]);  # zero for M-basis
  hstar.elm:=List(h.elm, t->[t[1],t[3],t[2]]);
  hstar.coeff:=Copy(h.coeff);
  CollectCoefficients(hstar);
  return hstar;
end,

## Return the coefficient of [s,t] in h                                 
Coefficient:=function(h,tab) local s, t, c;
  if IsStandardTableau(tab) then
    s:=HeckeAlgebraOps.Murphy.Tableau(Hecke(h),List(tab,Length),tab);
    c:=Position(h.elm,[s.mu,1,s.ind]);
  else
    s:=HeckeAlgebraOps.Murphy.Tableau(Hecke(h),List(tab,Length),tab[1]);
    t:=HeckeAlgebraOps.Murphy.Tableau(Hecke(h),s.mu,tab[2]);
    c:=Position(h.elm,[s.mu,s.ind,t.ind]);
  fi;
  if c=false then return 0*h.coeff[1];
  else return h.coeff[c];
  fi;
end,

## Given a standard tableau <t> for the Hecke algebra <H> return
## the permutation w, as a CoxeterGroup permutation, such that t=t^mu*w.
CoxeterPermTableau:=function(H,t)
  if not IsBound(t.wd) then
    t.wd:=EltWord(CoxeterGroup(H), CoxeterWordSymmPerm(
            MappingPermTableaux(H.Tableaux[t.mu][1],H.Tableaux[t.mu][t.ind])));
  fi;
  return t.wd;
end,

## Given a standard tableau <t> for the Hecke algebra <H> return
## the permatution w, as a CoxeterGroup word, such that t=t^mu*w.
## (we never use this function; it's included for completeness).
CoxeterWordTableau:=function(H,t)
  if not IsBound(t.wd) then
    t.wd:=EltWord(CoxeterGroup(H),CoxeterWordSymmPerm(
            MappingPermTableaux(H.Tableaux[t.mu][1],H.Tableaux[t.mu][t.ind])));
  fi;
  return CoxeterWord(CoxeterGroup(H), t.wd);
end

));  ## end of CreateHeckeBasis

########################################################################

## Compute the Gram matrix of a Specht module w.r.t. its Murphy basis.
GramMatrix:=function(H,mu) local M, s, t, g, h, tab;
  if not CHEVIE.MurphySpechtModules then
    Print("# WARNING: in the interests of speed, this function has just \n",
          "#          disabled T-basis to Murphy basis convertions.\n");
    CHEVIE.MurphySpechtModules:=true;
  fi;
  tab:=StandardTableaux(mu);
  M:=Basis(H,"Murphy");
  g:=List([1..Length(tab)], s->[]);
  for s in [1..Length(tab)] do
    for t in [s..Length(tab)] do
      h:=M(tab[1],tab[s])*M(tab[t],tab[1]);
      if h=0*h then g[s][t]:=0; else g[s][t]:=h.coeff[1]; fi;
      if s<>t then g[t][s]:=Copy(g[s][t]); fi;
    od;
  od;
  return g;
end;

#U MurphyOperators(H)....................................................
#M Returns a function, L say, for producing the Jucys-Murphy operators of
#M the Hecke algbebra H. These functions return the elements             
#M  L(i) := q^-1T(i-1)+q^-2T(i-1,i-2,i-1)+...+q^(1-i)T(i-1,...,1,...,i-1)
#M In fact, the function L() has a second optional argument so that      
#M  L(i,m) := L(i)-[m]_q*T().                                            
MurphyOperators:=function(H)
  H.MurphyOperators:=[0*Basis(H,"T")()];
  return function(arg)
    local i,T, q, Li, j;
    i :=arg[1];
    T:=Basis(H,"T");
    q:=H.parameter[1][1];
    if not IsBound(H.MurphyOperators[i]) then
      H.MurphyOperators[i]:=0*T();
      for j in [1..i-1] do
        H.MurphyOperators[i]:=q^-1*T(j)+q^-1*T(j)*H.MurphyOperators[i]*T(j);
      od;
    fi;
    if Length(arg)=1 then return H.MurphyOperators[i];
    elif IsInt(arg[2]) then return H.MurphyOperators[i]-QuantumInteger(q,arg[2])*T();
    else return H.MurphyOperators[i]-arg[2]*T();
    fi;
  end;
end;
MurphyOperatorsAffine:=function(H)
  H.MurphyOperatorsAffine:=[0*Basis(H,"T")()];
  return function(arg)
    local i,T, q, Li, j;
    i :=arg[1];
    T:=Basis(H,"T");
    q:=H.parameter[1][1];
    if not IsBound(H.MurphyOperatorsAffine[i]) then
      H.MurphyOperatorsAffine[i]:=(1-q^-1)*T(1)+T();
      for j in [2..i-1] do
        H.MurphyOperatorsAffine[i]:=q^-1*T(j)*H.MurphyOperatorsAffine[i]*T(j);
      od;
    fi;
    if Length(arg)=1 then return H.MurphyOperatorsAffine[i];
    else return H.MurphyOperatorsAffine[i]-QuantumInteger(q,arg[2])*T();
    fi;
  end;
end;

#U HeckeAlgebra(n).......................................................
#M Fast definition of the Hecke algebra for the symmetric group Sym(n).  
#HeckeAlgebra:=n->Hecke(CoxeterGroup("A",n-1),q);

#U SpechtModule(H,mu)....................................................
#M A function whic returns a function for working with the Murphy basis  
#M of M the Specht module S(<mu>).                                       
SpechtModule:=function(H,mu)
  local tmu, l, i;
  if Sum(mu)<>Length(H.parameter)+1 then 
    Error("<mu> must be a partition of ",Length(H.parameter)+1);
  fi;
  CHEVIE.MurphySpechtModules:=true;
  tmu:=[];
  l:=0;
  for i in mu do
    Add(tmu,[l+1..l+i]);
    l:=l+i;
  od;
  return t->Basis(H,"Murphy")(tmu,t);
end;

#U MurphyTest([n]).......................................................
#M Test the Murphy basis code by having it check that T(M(T(w)))=T(w) for
#M all w in Sym(n), for some <n>. If <n> is not supplied then it defaults
#M to the case n=4.                                                      
MurphyTest:=function(arg) local nTest, q, W, H, n, T, M, elts, t, w;
  if Length(arg)=0 then nTest:=[4];
  elif Length(arg)=1 and IsInt(arg[1]) then nTest:=[arg[1]];
  else
    nTest:=Flat(arg);
  fi;
  q:=Indeterminate(Rationals); q.name:="q";
  for t in nTest do
    if IsInt(t) then
      W:=CoxeterGroup("A",t-1);
      H:=Hecke(W,q);
    elif IsHeckeAlgebra(t) then
      H:=t;
      W:=Group(H);
    else Error("MurphyTest( n/H )");
    fi;
    n:=W.rank+1;
    T:=Basis(H,"T");
    M:=Basis(H,"Murphy");
    Print("Testing the Murphy basis functions for the Hecke algebra of Sym(",n,")\n");
    elts:=CoxeterWords(W);SortBy(elts, a->[Length(a),a]);
    for w in elts do
      Print(" - checking w=",w,"\n");
      if T(w)<>T(M(T(w))) then 
        Print("w = ", w,"\n");
        Print("Difference = ", T(w)-T(M(T(w))),".\n");
        Print("T(w) =",T(w),".\n");
        Print("M(T(w))=",M(T(w)),"\n");
        Print("T(M(T(w)))=",T(M(T(w))),"\n");
        Error(" murphy.g FAILED for Sym(",n,") at T(",Join(w),")  :(\n");
      fi;
    od;
    Print("\n ** murphy.g PASSED the test for Sym(",n,") !!\n");
  od;
end;

#U MurphyTest2([n])......................................................
#M Test the Murphy basis code by having it check that M(T(M(s,t)))=M(s,t)
#M all pairs of standard tabelaux (s,t) of the same shape.               
MurphyTest2:=function(arg) local nTest, q, W, H, T, M, std, n, mu, s, t;
  if Length(arg)=0 then nTest:=[4];
  elif Length(arg)=1 and IsInt(arg[1]) then nTest:=[arg[1]];
  else
    nTest:=Flat(arg);
    if not ForAll(nTest, IsInt) then Error("MurphyTest([n]"); fi;
  fi;
  q:=Indeterminate(Rationals); q.name:="q";
  for n in nTest do
    Print("Testing the Murphy basis functions for the Hecke algebra of Sym(",n,")\n");
    W:=CoxeterGroup("A",n-1);
    H:=Hecke(W,q);
    T:=Basis(H,"T");
    M:=Basis(H,"Murphy");
    #Sort(elts, function(a,b) return Length(a)<=Length(b); end);
    for mu in Reversed(Partitions(n)) do
      std:=StandardTableaux(mu);
      for s in std do
        for t in std do
          if M(s,t)<>M(T(M(s,t))) then
            Error(" murphy.g FAILED for Sym(",n,") at (s,t)=(",
                  VeryTightTableau(s),", ",VeryTightTableau(t),"):\n  ",
                  M(s,t)," --> ", M(T(M(s,t))),"\n");
          else
            Print("OK for (s,t)=(",
                  VeryTightTableau(s),", ",VeryTightTableau(t),")\n");
          fi;
        od;
      od;
    od;
    Print("\n ** murphy.g PASSED the test for Sym(",n,") !!\n");
  od;
end;

#U SpechtModel(H,mu).....................................................
#M Returns the representation of the Hecke algebra H of type A indexed by
#M partition mu, that is the list of matrices of the T_i.
SpechtModel:=function(H,mu) local T, S, tabs, n, zero, mats, St, ti, t,i;
  T:=Basis(H,"T");
  S:=SpechtModule(H,mu);
  tabs:=StandardTableaux(mu);
  for t in tabs do S(t); od; # initialize tableaux in H.Tableaux
  n:=Sum(mu)-1;
  zero:=Sum(H.parameter,Sum)*0;
  mats:=List([1..n],i->List(tabs, j->[1..Length(tabs)]*zero));
  for t in tabs do
    St:=S(t);
    for i in [1..n] do
      ti:=St*T(i);
      mats[i][St.elm[1][3]]{List(ti.elm,x->x[3])}:=ti.coeff;
    od;
  od;
  return mats;
end;
