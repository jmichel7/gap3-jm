###########################################################################
# Contribution to the Chevie Package ---- cp.g
#
# (C) July 2017 --- Jean Michel and Georges Neaime
#
# This  file defines the function 'CorranPicantinMonoid' which computes the
# interval  monoid generalizing  the Corran-Picantin  monoid, as defined in
# [1] G. Neaime, http://arxiv.org/abs/1707.06864.
#
# CorranPicantinMonoid(e,n[,k])  returns the monoid with delta the diagonal
# matrix  whose diagonal entries except the first are equal to E(e)^k; this
# monoid is the Corran-Picantin monoid for G(e,e,n) when Gcd(k,e)=1
###########################################################################

CorranPicantinMonoid:=function(arg)local e,n,k,M,W,i,deltamat,r,cr;
  e:=arg[1];n:=arg[2];
  if Length(arg)=2 then k:=1;
  else k:=arg[3] mod e;
    if k=0 then Error("k must not be divisible by e");fi;
  fi;
  r:=[Concatenation([-E(e),1],[1..n-2]*0)];
  for i in [1..n-1] do Add(r,Concatenation([1..i-1]*0,[-1,1],[1..n-i-1]*0));od;
  cr:=Copy(r);cr[1][1]:=-E(e)^-1;
  W:=PermRootGroup(r,cr);
  M:=rec(group:=W, operations:=rec(
    Print:=function(M)Print("CorranPicantinMonoid(",e,",",e,",",n,")");
    if k<>1 then Print("[",k,"]");fi;
  end));
  M.atoms:=[W.2,W.1];
  for i in [2..e-1] do Add(M.atoms,M.atoms[i-1]^M.atoms[i]); od;
  Append(M.atoms,W.generators{[3..n]});
  M.mat:=s->MatXPerm(W,s);
  # for the next function, see Proposition 3.11 in [1]
  M.IsLeftDescending:=function(simp,i)local m,c1,c2;
    m:=M.mat(simp);
    if i<=e then
      c1:=PositionProperty(m[1],x->x<>0);
      c2:=PositionProperty(m[2],x->x<>0);
      if c1<c2 then return m[2][c2]<>1;
      else return m[1][c1]=E(e)^-(i-1);
      fi;
    fi;
    i:=i-e+2;
    c1:=PositionProperty(m[i-1],x->x<>0);
    c2:=PositionProperty(m[i],x->x<>0);
    if c1<c2 then return m[i][c2]<>1;
    else return  m[i-1][c1]=1;
    fi;
  end;
  M.IsRightDescending:=function(simp,i)
    return M.IsLeftDescending(simp^-1,i);
  end;
  M.IsRightAscending:=function(simp,i) 
    return M.IsLeftDescending(simp^-1*M.delta,i);
  end;
  M.IsLeftAscending:=function(simp,i) 
    return M.IsRightDescending(M.delta*simp^-1,i);
  end;
  M.identity:=W.identity;
  deltamat:=IdentityMat(n)*E(e)^k;deltamat[1][1]:=E(e)^(k*(1-n));
  M.delta:=PermMatX(W,deltamat);
  M.stringDelta:="d";
  M.orderDelta:=Order(W,M.delta);
  M.showmat:=function(s)local t;
    if IsPerm(s) then PrintArray(M.mat(s));
    else for t in s.elm do M.showmat(t);od;
    fi;
  end;
  if k<>e-k then
    if Length(arg)=4 then M.revMonoid:=arg[4];
    else M.revMonoid:=CorranPicantinMonoid(e,n,e-k,M);
    fi;
    M.Reverse:=function(b)local res,s;
      if Length(b.elm)=0 then return M.revMonoid.Elt([],b.pd);fi;
      res:=[];
      for s in List(Reversed(b.elm),y->M.revMonoid.DeltaAction(y^-1,b.pd))
      do res:=M.revMonoid.AddToNormal(res,s);od;
      return GarsideEltOps.Normalize(M.revMonoid.Elt(res,b.pd));
    end;
  fi;
  CompleteGarsideRecord(M,rec(interval:=true));
  return M;
end;
