#############################################################################
##
#A  cycpol.g           CHEVIE library     Jean Michel 13-9-96 and later
##
#Y  Copyright (C)  1992-2012 Lehrstuhl D f\"ur  Mathematik, RWTH Aachen,
#Y  and  University Paris VII.
##
##############################################################################
## Useful function provided by Thomas Breuer in may 1999
##
#F  AsRootOfUnity( <c> )
##
##  Given  a cyclotomic c  `AsRootOfUnity' returns the  rational e/n if c is
##  equal to E(n)^e, and false otherwise. In particular returns 0 for 1.
##
##############################################################################
##  Suppose that we know the coefficients of $\pm E(n)^i$ w.r.t. the $n$-th
##  Zumbroich  basis (see~"ZumbroichBasis"). These  coefficients are either
##  all  $1$ or all $-1$. More precisely, they arise in the base conversion
##  from  (formally)  successively  multiplying  $\pm  E(n)^i$  by  $1  = -
##  \sum_{j=1}^{p-1}  E(n)^{jn/p}$, for suitable prime divisors $p$ of $n$,
##  and  then treating the summands  $\pm E(n)^{i + jn/p}$  in the same way
##  until roots in the basis are reached. It should be noted that all roots
##  obtained this way are distinct.
##
##  Suppose  the above procedure must be  applied for the primes $p_1, p_2,
##  \ldots,  p_r$. Then $ E(n)^i$  is equal to $(-1)^r \sum_{j_1=1}^{p_1-1}
##  \cdots  \sum_{j_r=1}^{p_r-1}  E(n)^{i  +  \sum_{k=1}^r j_k n/p_k}$. The
##  number of summands is $m = \prod_{k=1}^r (p_k-1)$. Since these root are
##  all  distinct, we can compute the sum $s$ of their exponents modulo $n$
##  from  the known coefficients of $ E(n)^i$, and we get $s \equiv m ( i +
##  r n/2 ) \pmod{n}$. Either $m = 1$ or $m$ is even, hence this congruence
##  determines $ E(n)^i$ at most up to its sign.
##
##  Now  suppose  that  $g  =  \gcd(  m,  n  )$  is nontrivial. Then $i$ is
##  determined  only modulo $n/g$, and  we have to decide  which of the $g$
##  possible  values $i$  is. This  could be  done by  computing the values
##  $j_{0,p}$ for one candidate $i$, where $p$ runs over the prime divisors
##  $p$ of $n$ for which $m$ is the product of $(p-1)$.
##
##  (Recall  that each $n$-th root of unity  is of the form $\prod_{p\in P}
##  \prod_{k_p=1}^{\nu_p-1}  E(n)^{j_{k,p}n/p^{k_p+1}}$, with  $j_{0,p} \in
##  \{  0, 1, \ldots, p-1  \}$, $j_{k,2} \in \{  0, 1 \}$ for  $k > 0$, and
##  $j_{k,p}  \in \{ -(p-1)/2, \ldots, (p-1)/2 \}$ for $p > 2$. The root is
##  in the $n$-th Zumbroich basis if and only if $j_{0,2} = 0$ and $j_{0,p}
##  \not= 0$ for $p > 2$.)
##
##  But note that we do not have an easy access to the decomposition of $m$
##  into factors $(p-1)$, although in fact this decomposition is unique for
##  $n \leq 65000$.
##
##  So the exponent is identified by dividing off one of the candidates and
##  then  identifying the quotient, which  is a $g$-th or  $2 g$-th root of
##  unity. (Note that $g$ is small compared to $n$. $m$ divides the product
##  of  $(p-1)$  for  prime  divisors  $p$  that  divide  $n$ at least with
##  exponent  $2$. The maximum  of $g$ for  $n \leq 65000$  is $40$, so $2$
##  steps suffice for these values of $n$.)

AsRootOfUnity := function( root )

    local coeffs,   # Zumbroich basis coefficients of `root'
          n,        # conductor of `n'
          sum,      # sum of exponents with nonzero coefficients
          num,      # number of nonzero coefficients
          i,        # loop variable
          val,      # one coefficient
          coeff,    # one nonzero coefficient (either `1' or `-1')
          exp,      # candidate for the exponent
          g,        # `Gcd( n, num )'
          groot;    # root in recursion

    # Added JM 19/5/2003: reject non-roots
    if root*GaloisCyc(root,-1)<>1 then return false;fi;

    # Handle the trivial cases that `root' is an integer.
    if root = 1 then return 0;
    elif root = -1 then return 1/2;
    fi;

    # Compute the Zumbroich basis coefficients,
    # and number and sum of exponents with nonzero coefficient (mod `n').
    coeffs:= COEFFSCYC( root );
    if ForAny(coeffs,x->not x in [0,1,-1]) then return false;fi;
    n:= Length( coeffs );
    sum:= 0;
    num:= 0;
    for i in [ 1 .. n ] do
      val:= coeffs[i];
      if val <> 0 then
        sum:= sum + i; num:= num + 1; coeff:= val;
      fi;
    od;
    sum:= sum - num;

    # `num' is equal to `1' if and only if `root' or its negative
    # belongs to the basis.
    # (The coefficient is equal to `-1' if and only if either
    # `n' is a power of $2$ or
    # `n' is odd and `root' is a primitive $2 `n'$-th root of unity.)
    if num = 1 then
      if coeff < 0 then
        if n mod 2 = 0 then sum:= sum + n/2;
        else sum:= 2*sum + n; n:= 2*n; sum:= sum mod n;
        fi;
      fi;
      if root=E(n)^sum then return (sum mod n)/n; else return false;fi;
    fi;

    # Let $N$ be `n' if `n' is even, and equal to $2 `n'$ otherwise.
    # The exponent is determined modulo $N / \gcd( N, `num' )$.
    g:= GcdInt( n, num );
    if g = 1 then

      # If `n' and `num' are coprime then `root' is identified up to its sign.
      exp:= ( sum / num ) mod n;
      if root <> E(n)^exp then
        exp:= 2*exp + n; n:= 2*n; exp:= exp mod n;
      fi;

    elif g = 2 then

      # `n' is even, and again `root' is determined up to its sign.
      exp:= ( sum / num ) mod ( n / 2 );
      if root <> E(n)^exp then exp:= ( exp + n / 2 ) mod n; fi;

    else

      # Divide off one of the candidates.
      # The quotient is a `g'-th or $2 `g'$-th root of unity,
      # which can be identified by recursion.
      exp:= ( sum / num ) mod ( n / g );
      groot:= AsRootOfUnity( root * E(n)^(n-exp) );
      if groot=false then return false;fi;
      if n mod 2 = 1 and Denominator(groot) mod 2 = 0 then
        exp:= 2*exp; n:= 2*n;
      fi;
      exp:= ( exp + Numerator(groot) * n / Denominator(groot) ) mod n;

    fi;

    if root=E(n)^exp then return (exp mod n)/n; else return false;fi;
end; 

#############################################################################
##  The  rest  of  this  file  contains  functions  dealing  with  rational
##  fractions  with all poles and zeroes equal to 0 or roots of unity. They
##  are  represented as  'CycPol', which  are used  in various files of the
##  CHEVIE package whenever possible.
##
##  A CycPol is a record of the form:
##  (coeff:=c,valuation:=v,vcyc:=vec,operations:=CycPolOps)
##
##  it  represents a rational fraction in one variable (which by default is
##  printed  as 'q'  but this  can be  changed, see  'Format') of  the form
##  c*q^v*product(terms for the entries in the list vec)
##
##  where:
##  c   can be anything (an element of a suitable ring...)
##      This  gives  a  great  flexibility  of  use  (e.g.  one  can accept
##      arbitrary   polynomials  by  putting  the  non-cyclotomic  part  in
##      coeff)...
##  vec is a  list of entries. An  entry is of the form  [k/d, m] where
##      k<d is a rational in [0,1[ and m is an integer; it represents
##      (q-E(d)^k)^m
##
##  The  entries  [k/d,m]  are  sorted  by  [Denominator(k/d), k/d, m], and
##  CycPolOps.normalize  makes  sure  all  entries  have  k/d  different by
##  gathering factors.
##
##  The  advantage  of  representing  rational  fractions  as  CycPol  when
##  possible  should  be  obvious:  less  storage,  faster  multiplication,
##  division  and evaluation, and not  least nice display (factorized). The
##  big  disadvantage is that addition  and subtraction are not implemented
##  :)
#############################################################################

CycPolOps:= OperationsRecord("CycPolOps");

IsCycPol:=x->IsRec(x) and IsBound(x.operations) and x.operations=CycPolOps;

# makes sure a in in normal form
CycPolOps.normalize:=function(a)
  if a.coeff=0 then a.vcyc:=[];a.valuation:=0;fi;
  if a.vcyc=[] then return a;fi;
  a.vcyc:=CollectBy(a.vcyc,x->[Denominator(x[1]),Mod1(x[1])]);
  a.vcyc:=List(a.vcyc,x->[Mod1(x[1][1]),Sum(x,y->y[2])]);
  a.vcyc:=Filtered(a.vcyc,x->x[2]<>0);
  return a;
end;

CycPolOps.\=:=function(a,b)
  if not IsRec(a) or not IsRec(b) or not IsBound(a.operations) or
   not IsBound(b.operations) or a.operations<>b.operations then
   return false;fi;
  return a.coeff=b.coeff and  a.valuation=b.valuation and a.vcyc=b.vcyc;
end;

CycPolOps.\*:=function(a,b)
  if IsCyc(a) then b:=ShallowCopy(b);b.coeff:=a*b.coeff;
    b.vcyc:=ShallowCopy(b.vcyc);return b; fi;
  if IsCyc(b) then a:=ShallowCopy(a);a.coeff:=b*a.coeff;
    a.vcyc:=ShallowCopy(a.vcyc);return a; fi;
  return CycPolOps.normalize(rec(coeff:=a.coeff*b.coeff,
                        valuation:=a.valuation+b.valuation,
   vcyc:=Concatenation(a.vcyc,b.vcyc),operations:=CycPolOps));
end;

CycPolOps.\/:=function(a,b)local res;
  if IsList(a) then return List(a,x->x/b);fi;
  if IsCyc(a) then return rec(coeff:=a/b.coeff,valuation:=-b.valuation,
               vcyc:=List(b.vcyc,v->[v[1],-v[2]]),operations:=CycPolOps);
  fi;
  if IsCyc(b) then return rec(coeff:=a.coeff/b,valuation:=a.valuation,
                  vcyc:=ShallowCopy(a.vcyc),operations:=CycPolOps);
  fi;
  return CycPolOps.normalize(rec(
    coeff:=a.coeff/b.coeff,
    valuation:=a.valuation-b.valuation,
    vcyc:=Concatenation(a.vcyc,List(b.vcyc,v->[v[1],-v[2]])),
    operations:=CycPolOps));
end;

CycPolOps.\^:=function(a,n)local res;
  res:=ShallowCopy(a);
  res.coeff:=res.coeff^n;res.valuation:=n*res.valuation;
  res.vcyc:=List(res.vcyc,x->[x[1],n*x[2]]);
  return res;
end;

CycPolOps.\<:=function(a,b)
  if not IsCycPol(a) then return a<b.vcyc;
  elif not IsCycPol(b) then return a.vcyc<b;
  elif a.coeff<>b.coeff then return a.coeff<b.coeff;
  elif a.valuation<>b.valuation then return a.valuation<b.valuation;
  else return a.vcyc<b.vcyc;
  fi;
end;

CycPolOps.Value:=function(v,q)local l,res,foo;
  foo:=function(t,q)return (q-E(Denominator(t[1]))^Numerator(t[1]))^t[2];end;
  l:=CollectBy(v.vcyc,x->Denominator(x[1]));
  res:=q^v.valuation*Product(l,function(p)local r;
    if IsPolynomial(q) and ForAll(p,t->t[2]>=0) then
      r:=Product(p,t->foo(t,Indeterminate(Cyclotomics)));
      if q<>Indeterminate(Cyclotomics) then r:=Value(r,q);fi;
    else r:=List(p,t->foo(t,q));
      if false in r then
	Error("Cannot evaluate the non-Laurent polynomial CycPol ",v);
      else r:=Product(r);
      fi;
    fi;
    return r;end);
  if IsRec(v.coeff) and IsBound(v.coeff.operations) and
    IsBound(v.operations.Value) then return Value(v.coeff,q)*res;
  else return v.coeff*res;
  fi;
end;

# Fast routine for CycPol(Value(p,q*e)) for a root of unity e
CycPolOps.EnnolaTwist:=function(v,e)local res;res:=ShallowCopy(v);
  res.coeff:=res.coeff*e^Degree(res);
  e:=AsRootOfUnity(e); res.vcyc:=List(res.vcyc,p->[Mod1(p[1]-e),p[2]]);
  return CycPolOps.normalize(res);
end;

# Fast routine for  CycPol(Value(p,q^n))
CycPolOps.DescentOfScalars:=function(p,n)local res,c;
  res:=ShallowCopy(p);res.valuation:=n*p.valuation;
  if n>0 then
    res.vcyc:=Concatenation(List(p.vcyc,x->List(x[1]/n+[0..n-1]/n,r->[r,x[2]])));
  elif n=0 then return CycPol(Value(p,1));
  else
    res.valuation:=res.valuation+n*Sum(p.vcyc,x->x[2]);
    res.coeff:=res.coeff*(-1)^Sum(p.vcyc,x->x[2]);
    c:=Sum(p.vcyc,y->y[1]*y[2]);
    res.coeff:=res.coeff*E(Denominator(c))^Numerator(c);
    res.vcyc:=Concatenation(List(p.vcyc,x->List(x[1]/n+[0..-n-1]/-n,r->[r,x[2]])));
  fi;
  return CycPolOps.normalize(res);
end;
    
CycPolOps.Degree:=function(a)local deg;
  deg:=a.valuation+Sum(a.vcyc,t->t[2]);
  if IsPolynomial(a.coeff) then deg:=deg+Degree(a.coeff);fi;
  return deg;
end;

# Valuation(p [,d]) if d absent return valuation
# else  return d-valuation i.e. returns power of (q-\zeta)-factor of CycPol
# p where if IsInt(d) and d<>0 then \zeta=E(d) else d=AsRootOfUnity(zeta)
CycPolOps.Valuation:=function(arg)local p,d;
  p:=arg[1];
  if Length(arg)=1 then return p.valuation;fi;
  d:=arg[2];if IsInt(d) and d<>0 then d:=Mod1(1/d);fi;
  return Sum(Filtered(p.vcyc,x->x[1]=d),x->x[2]);
end;

CycPolOps.decs:=[,,,,,,,
  [[1,5],[3,7],[1,7],[3,5],[1,3],[5,7]],,,,  #8
  [[1,5],[7,11],[1,7],[5,11],[1,11],[5,7],[1],[5],[7],[11]],,,  #12
  [[1,4,11,14],[2,7,8,13],[1,4,7,13],[2,8,11,14],
   [1,4],[7,13],[11,14],[2,8]],  #15
  [[1,7,9,15],[3,5,11,13]],,,,  #16
  [[1,9,11,19],[3,7,13,17],[1,9,13,17],[3,7,11,19]],  #20
  [[1,4,10,13,16,19],[2,5,8,11,17,20]],,,  #21
  [[1,7,13,19],[5,11,17,23],[1,7,17,23],[5,11,13,19],  #24
   [1,5,19,23],[7,11,13,17],[1,11,17,19],[5,7,13,23],
   [7,13],[1,19],[5,23],[11,17]],,,,,,  
  [[1,11,19,29],[7,13,17,23],[1,7,13,19],[11,17,23,29],
   [11,29],[17,23],[1,19],[7,13]],,,,,,,,,,,,  #30
  [[1,13,19,25,31,37],[5,11,17,23,29,41]]];  #42

# returns  the  list  of  factors  of  Phi_d(q)  which  are  recognized  when
# formatting CycPols.
# Each factor is returned as the list of i such that its roots are E(d)^i.
# The results are cached in CycPolOps.decompositions to go faster. 
# The first factor is Phi_d(q). Then if d=p^e or d=2p^e then Z/d^\times is 
# cyclic and there are 2 factors Phi' and Phi'' corresponding  to the product 
# of the even (resp odd) powers of a primitive d-th root of unity.
# For  small values  of d  (<=42) there are more
# factors given by the explicit list 'decs' below.
CycPolOps.PhiDecompositions:=function(d)local dec,r,phi;
  if not IsBound(CycPolOps.decompositions) then CycPolOps.decompositions:=[];fi;
  if not IsBound(CycPolOps.decompositions[d]) then
    dec:=[PrimeResidues(d)];
    if d=1 then dec:=[[1]];
    elif d>2 then r:=PrimitiveRootMod(d);
      if r<>false then phi:=[0,2..Phi(d)-2];
	Append(dec,[List(phi,i->r^i mod d),List(phi+1,i->r^i mod d)]);
      fi;
    fi;
    if IsBound(CycPolOps.decs[d]) then Append(dec,CycPolOps.decs[d]);fi;
    CycPolOps.decompositions[d]:=dec;
  fi;
  return CycPolOps.decompositions[d];
end;

# dumps a.vcyc to a list on the way to formatting CycPol a. Returns list of:
#  [n,d]     specifying Phi_d^n
#  [n,d,i]   specifying (i-th factor of Phi_d)^n
#  [n,i/d]   specifying (q-E(d)^i)^n
CycPolOps.decompose:=a->Concatenation(List(
  CollectBy(a.vcyc,t->Denominator(t[1])),
  function(t)local d,v,res,r,i,dec,n;
    d:=Denominator(t[1][1]);
    if d=1 then return [[t[1][2],d]];fi;
    res:=[];
    v:=[1..d]*0;v{List(t,p->Numerator(p[1]))}:=List(t,p->p[2]);
    dec:=CycPolOps.PhiDecompositions(d);
    for i in [1..Length(dec)] do
      r:=dec[i]; n:=Minimum(v{r});
      if n>0 then v{r}:=v{r}-n;
      else n:=Maximum(v{r});if n<0 then v{r}:=v{r}-n; else n:=0;fi;
      fi;
      if n<>0 then if i=1 then Add(res,[n,d]);else Add(res,[n,d,i-1]);fi;fi;
    od;
    for i in [1..d] do if v[i]<>0 then Add(res,[v[i],i/d]);fi;od;
    return res;
  end));

# options:
# by defaut print coeff.q^valuation.product of  Phi-factors
#            (see PhiDecompositions for recognized factors)
# TeX     same but in TeX
# expand  print in extenso of each cyclotomic polynomial
# vname   use for variable name
CycPolOps.Format:=function(a,option)local res,var,m;
  if IsBound(option.vname) then var:=option.vname;
  elif IsBound(a.vname) then var:=a.vname;
  elif IsBound(Indeterminate(Cyclotomics).name) then 
    var:=Indeterminate(Cyclotomics).name;
  else var:="q";
  fi;
  res:=List(CycPolOps.decompose(a),function(e)local p,n,c;
    n:=e[1];p:=e[2];
    if IsBound(option.GAP) and not IsBound(option.expand) then
      if Length(e)=3 then 
        return Concatenation(List([1..n],
	    i->CycPolOps.PhiDecompositions(p)[e[3]]/p));
      else return List([1..n],x->p);
      fi;
    elif Length(e)=3 then
#     AddSet(CycPolOps.Used,e{[2,3]});
      if IsBound(option.expand) then
	c:=Product(CycPolOps.PhiDecompositions(p)[e[3]+1],
	   l->Indeterminate(Cyclotomics)-E(p)^l).coefficients;
	c:=Concatenation("(",FormatPolynomial(c,0,var,option),")");
      elif IsBound(option.TeX) then 
	if e[3]>4 then c:=SPrint("{\\Phi^{(",e[3],")}_",TeXBracket(p),"}");
	else c:=SPrint("{\\Phi",List([1..e[3]],y->'\''),"_",TeXBracket(p),"}");
	fi;
      elif e[3]>6 then c:=SPrint("P(",e[3],")",p);
      else c:=List([1..QuoInt(e[3],2)],y->'"');
	if e[3] mod 2=1 then Add(c,'\'');fi;
	c:=SPrint("P",c,p);
      fi;
    elif IsInt(p) then
      if IsBound(option.expand) then
	c:=SPrint("(",FormatPolynomial(CyclotomicPol(p),0,var,option),")");
      elif IsBound(option.TeX) then c:=SPrint("\\Phi_",TeXBracket(p));
      else c:=SPrint("P",p);
      fi;
    else 
      c:=Concatenation("(",
       FormatPolynomial([-E(Denominator(p))^Numerator(p),1],0,var,option),")");
    fi;
    return FormatMonomial(c,n,option);
  end);
  if IsBound(option.expand) and IsBound(option.GAP) then res:=Join(res,"*");
  else res:=Concatenation(res);fi;
  if IsBound(option.GAP) and not IsBound(option.expand) then return
    SPrint("CycPol(",FormatGAP(Flat([a.coeff,a.valuation,res])),")");fi;
  m:=FormatMonomial(var,a.valuation,option);
  if m<>"1" then res:=Concatenation(m,res);fi;
  if IsRec(a.coeff)then 
       return FormatCoefficient(Value(a.coeff,Mvp(var)),res,option);
  else return FormatCoefficient(a.coeff,res,option);
  fi;
end;

CycPolOps.String:=a->Format(a);

CycPolOps.Print:=function(a)Print(String(a));end;

CycPolOps.ComplexConjugate:=function(a)
  a:=ShallowCopy(a);a.coeff:=ComplexConjugate(a.coeff);
  a.vcyc:=List(a.vcyc,v->[Mod1(1-v[1]),v[2]]);
  return CycPolOps.normalize(a);
end;

###########################################################################
#  CycPol "constructor" accepts the following forms
#  CycPol(vector v) 
#     returns a CycPol with  coeff:=v[1] and valuation:=v[2],
#     vcyc:=v{3..Length(v)} if v[3] is a list
#     otherwise  v{3..Length(v)} should be a  list of integers or rationals
#     mod1.   A  rational  p/q  represents  (X-E(q)^p)  and  an  integer  p
#     represents Phi_p(X).
#  CycPol(polynomial over the cyclotomics or a subfield) 
#    returns the CycPol representing the polynomial
#  CycPol(cyclotomic)
#    for the special case of a constant polynomial.
#  CycPol(obj)
#    calls obj.operations.CycPol if it exists otherwise.
##
CycPol:=function(p)
  local res,d,a,l,i,r,e,n,q,tested,try,testcyc,testall,conductor,bounds,found;
  if IsList(p) then
    res:=rec(coeff:=p[1],valuation:=p[2],vcyc:=[],operations:=CycPolOps);
    if Length(p)>2 and IsList(p[3]) then res.vcyc:=p{[3..Length(p)]};
    else res.vcyc:=[];
      for a in p{[3..Length(p)]} do
        if IsInt(a) and a<>0 then Append(res.vcyc,PrimeResidues(a)/a);
	else Add(res.vcyc,a);
	fi;
      od;
      res.vcyc:=Collected(res.vcyc);
    fi;
  elif IsCyc(p) then
    res:=rec(coeff:=p,valuation:=0,vcyc:=[],operations:=CycPolOps);
  elif IsPolynomial(p) then
    # lot of code to be as efficient as possible in all cases
    res:=rec(valuation:=p.valuation,vcyc:=[],operations:=CycPolOps);
    q:=Indeterminate(Cyclotomics);
    if IsBound(q.name) then res.vname:=q.name;fi;
    if Length(p.coefficients)=0 then    # p=0
      res.coeff:=0;
    elif Length(p.coefficients)=1 then  # p=ax^s
      res.coeff:=p.coefficients[1];
    elif 2=Number(p.coefficients,x->x<>x*0) then # p=ax^s+bx^t
      d:=Length(p.coefficients)-1;
      res.coeff:=p.coefficients[d+1];
      r:=-p.coefficients[1]/res.coeff;
      a:=AsRootOfUnity(r);
      if a=false then res.coeff:=Value(p,q)/q^p.valuation;
      else res.vcyc:=List([0..d-1],i->[(a+i)/d,1]);
      fi;
    else
      res.coeff:=p.coefficients[Length(p.coefficients)];
      p:=Polynomial(Cyclotomics,p.coefficients/res.coeff); # now p is monic

      # returns the list of i such that Phi(i)/Phi(Gcd(i,conductor))<=deg pol
      bounds:=function(pol)local f,t,p,tp,pw,pp,t1,prod,conductor,d;
	d:=Degree(pol);
        if d=0 then return [];fi;
        conductor:=Lcm(List(pol.coefficients,NofCyc));
	if conductor=1 then f:=[];else f:=Collected(Factors(conductor));fi;
	t:=[];t1:=[];
	for p in f do
	  tp:=[1];pw:=p[1];while pw<=d do Add(tp,pw);pw:=pw*p[1];od;
	  Add(t,tp);Add(t1,tp*p[1]^p[2]);
	od;
	pp:=Filtered(Primes,x->x<=d+1);pp:=Difference(pp,List(f,x->x[1]));
	for p in pp do
	  tp:=[1,p-1];pw:=p*(p-1);while pw<=d do Add(tp,pw);pw:=pw*p;od;
	  Add(t,tp);Add(t1,List([1..Length(tp)]-1,i->p^i));
	od;
	prod:=function(l,d)local p;
	  p:=Filtered([1..Length(l[1])],x->l[1][x]<=d);
	  if Length(l)=1 then return List(p,x->[x]);fi;
	  return Concatenation(List(p,i->List(prod(l{[2..Length(l)]},d/l[1][i]),
	    x->Concatenation([i],x))));
	end;
	p:=List(prod(t,d),l->Product(Zip(l,t1,function(i,j)return j[i];end)));
	p:=Union(List(p,DivisorsInt));
	p:=Difference(p,tested);
	SortBy(p,x->x/Length(DivisorsInt(x)));
	return p;
      end;

      # find factors Phi_i
      testcyc:=function(i)local qr,found;
        found:=false;
	while true do
	  qr:=QuotientRemainder(p,CyclotomicPolynomial(Cyclotomics,i));
	  if qr[2]=0*qr[2] then Append(res.vcyc,PrimeResidues(i)/i);
	    p:=qr[1];found:=true;
	  else return found;
	  fi;
	od;
      end;

      # find other primitive i-th roots of unity
      testall:=function(i)local r,found;found:=false;
        for r in PrimeResidues(i) do
	  while ValuePol(p.coefficients,E(i)^r)=0 do # faster than QuotientRemainder
	    found:=true;
	    p:=p/(q-E(i)^r);Add(res.vcyc,r/i);
	    conductor:=Lcm(List(p.coefficients,NofCyc));
	    if Degree(p)<Phi(i)/Phi(Gcd(i,conductor)) then return found;fi;
	  od;
	od;
	return found;
      end;

      # first try commonly occuring fields
      tested:=[1,2,4,3,6,8,12,5,10,9,18,24,16,20,7,14,15,30,36,28,21,42];
      for i in tested do testcyc(i); testall(i);od;
      
      # if not finished do a general search.
      # p is in Q(zeta_conductor)[x] so can only have a root in mu_i for i below
      try:=bounds(p);
#     Print("tested=",tested,"\n");
      conductor:=Lcm(List(p.coefficients,NofCyc));
      i:=1;
      while i<=Length(try) do 
	Add(tested,try[i]);
	if conductor=1 then # All factors are Phi_i
	     found:=testcyc(try[i]);
	else found:=testall(try[i]);
	fi;
	if found  then 
	  try:=bounds(p);i:=1;
#	  Print("tested=",tested,"\n");
	else i:=i+1;
	fi;
      od;
      res.vcyc:=Collected(res.vcyc);
      if Degree(p)=0 then res.coeff:=res.coeff*p.coefficients[1];
      else res.coeff:=res.coeff*p;
      fi;
    fi;
  elif IsRec(p) and IsBound(p.operations) and IsBound(p.operations.CycPol) then
    return p.operations.CycPol(p);
  else Error("CycPol(",p,") not implemented");
  fi;
  return CycPolOps.normalize(res);
end;

LcmCycPol:=function(arg)local res,lcm2,p;
  lcm2:=function(a,b)
   a:=rec(coeff:=1,valuation:=Maximum(a.valuation,b.valuation),
    vcyc:=Concatenation(a.vcyc,b.vcyc),operations:=CycPolOps);
   a.vcyc:=CollectBy(a.vcyc,x->[Denominator(x[1]),Mod1(x[1])]);
   a.vcyc:=List(a.vcyc,x->[Mod1(x[1][1]),Maximum(List(x,y->y[2]))]);
   a.vcyc:=Filtered(a.vcyc,x->x[2]<>0);
   return a;
  end;
  res:=arg[1];
  for p in arg{[2..Length(arg)]} do res:=lcm2(res,p);od;
  return res;
end;
