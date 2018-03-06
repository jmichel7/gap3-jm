RevPrefix:=b->PreferredPrefix(ReversedWord(b))^-1;

# minimal element conjugating  x to its sss
MinimalToSSS:=function(b)local sb,sinf,ssup,min,delta,star,M;M:=b.monoid;
  star:=x->delta/x;
  sb:=RepresentativeSC(b).circuit[1];delta:=b.monoid.Elt([],1);
  sinf:=sb.pd;ssup:=sb.pd+Length(sb.elm);
  min:=b^0;
  while sinf>b.pd or ssup<b.pd+Length(b.elm) do
  if sinf>b.pd and ssup<b.pd+Length(b.elm) then
    b:=Sliding(b);min:=star(M.Elt([b[2]]))*min;b:=b[1]^delta;
  elif sinf>b.pd then min:=star(iota(b))*min;b:=(b^iota(b))^delta;
  else min:=star(iota(b^-1))*min;b:=(b^iota(b^-1))^delta;
  fi;
  od;
  return min;
end;

# testconj(b[,nb repetitions=]) test time computing various conjugacy sets
testconj:=function(arg)local b,v,n; b:=arg[1];v:=[];
  if Length(arg)=1 then n:=10;else n:=arg[2];fi;
  Add(v,Sum([1..n],i->[Length(ConjugacySet(b,"SS")),Dtime()])/n);
  Add(v,Sum([1..n],i->[Length(ConjugacySet(b,"USS")),Dtime()])/n);
  Add(v,Sum([1..n],i->[Length(ConjugacySet(b,"SC")),Dtime()])/n);
  Print("\n",FormatTable(v,rec(rowLabels:=["SS","USS","SC"],
    columnLabels:=["Size","Time"])));
end;

STransport:=function(x,s)local p,M;M:=x.monoid;
  if IsList(s) then return List(s,y->STransport(x,y));fi;
  p:=PreferredPrefix(x);
  return M.Elt([p^-1*GarsideAlpha(s)*PreferredPrefix(x^s)]);
end;

STransportN:=function(N,x,s)local i,p,M;M:=x.monoid;
  s:=GarsideAlpha(s);
  for i in [1..N] do
    p:=PreferredPrefix(x);
    s:=p^-1*s*PreferredPrefix(PositiveSimpleConjugation(x,s));
    x:=PositiveSimpleConjugation(x,p);
  od;
  return M.Elt(s);
end;

# orbit of braid b under function f
History:=function(b,f)local h;
  h:=[[b,b^0]];
  while true do
    b:=f(b);
    if b[1] in List(h,x->x[1]) then return Concatenation(h,[b]);fi;
    Add(h,b);b:=b[1];
  od;
end;

# canonical length of braid b
CanonicalLength:=b->Length(b.elm);

foo:=function(x)local h;h:=History(x,Sliding);
# return CanonicalLength(h[Length(h)][1])<CanonicalLength(x) and CanonicalLength(h[2][1])= CanonicalLength(x);
  return h[Length(h)][1].pd>x.pd and CanonicalLength(h[2][1])= CanonicalLength(x);
end;

maxss:=function(x)local h;h:=History(x,Sliding);
  return PositionProperty(h,x->CanonicalLength(x[1])=CanonicalLength(h[Length(h)][1]));
end;

iota:=function(b)local M;
 if Length(b.elm)=0 then return b^0;fi;
 M:=b.monoid; return M.Elt([M.DeltaAction(b.elm[1],-b.pd)]);
end;

minimal2:=function(b)local sb,sinf,ssup,min,delta,star;
  star:=x->delta/x;
  sb:=RepresentativeSC(b).circuit[1];delta:=b.monoid.Elt([],1);
  sinf:=sb.pd;ssup:=sb.pd+Length(sb.elm);
  min:=b^0;
  while sinf>b.pd do min:=star(iota(b))*min;b:=(b^iota(b))^delta;od;
  while ssup<b.pd+Length(b.elm) do
     min:=star(iota(b^-1))*min;b:=(b^iota(b^-1))^delta;od;
  return min;
end;

minimal3:=function(b)local sb,sinf,ssup,min,delta,star,M;M:=b.monoid;
  star:=x->delta/x;
  sb:=RepresentativeSC(b).circuit[1];delta:=b.monoid.Elt([],1);
  sinf:=sb.pd;ssup:=sb.pd+Length(sb.elm);
  min:=b^0;
  while sinf>b.pd or ssup<b.pd+Length(b.elm) do
    b:=Sliding(b);min:=star(M.Elt([b[2]]))*min;b:=b[1]^delta;od;
  return min;
end;

minimalr:=function(b)local sb,sinf,ssup,min,delta,star;
  star:=x->delta/x;
  sb:=RepresentativeSC(b).circuit[1];delta:=b.monoid.Elt([],1);
  sinf:=sb.pd;ssup:=sb.pd+Length(sb.elm);
  min:=b^0;b:=ReversedWord(b);
  while sinf>b.pd do min:=star(iota(b))*min;b:=(b^iota(b))^delta;od;
  while ssup<b.pd+Length(b.elm) do
     min:=star(iota(b^-1))*min;b:=(b^iota(b^-1))^delta;od;
  return ReversedWord(min);
end;

# History(b,x->Cycling(x)^Delta)
CSHistory:=function(b)local h;
  h:=[[b,b^0]];
  while true do
    b:=Cycling(b);
    b[1].elm:=List(b[1].elm,x->b[1].monoid.DeltaAction(x,1));
    b[2]:=b[1].monoid.Elt([b[1].monoid.delta/h[Length(h)][1].elm[1]]);
    if b[1] in List(h,x->x[1]) then return Concatenation(h,[b]);fi;
    Add(h,b);b:=b[1];
  od;
end;

# pullback starting at z of s starting a sliding(z)
# assumes z is in SSS
SPullback:=function(z,s)local y,M;
  y:=Sliding(z);M:=z.monoid;
  return RightGcd(y[2]*s,ReversedWord(Sliding(ReversedWord(y[1]^s))[2]))[2];
end;

SPullback2:=function(z,s)local y,M;
  y:=Sliding(z);M:=z.monoid;
  return RightGcd(y[2]*s,ReversedWord(Sliding(ReversedWord(y[1]^s))[2])){[2,3]};
end;

Hom:=function(x,y)local M,s;
  M:=x.monoid;
  s:=LeftDivisorsSimple(M,M.delta);
  s:=List(Concatenation(s),x->M.B(x));
  return Filtered(s,z->x^z=y);
end;

PreimageSliding:=function(b)local d,M;M:=b.monoid;
  d:=List(Concatenation(LeftDivisorsSimple(M,GarsideAlpha(ReversedWord(b)))),
    x->x^-1);
  d:=List(d,function(x)local y;y:=b.monoid.Elt([x]);return y*b*y^-1;end);
  return Set(Filtered(d,x->Sliding(x)[1]=b));
end;

PreimageSlidingSSS:=b->Filtered(PreimageSliding(b),
  x->Length(x.elm)=Length(b.elm) and x.pd=b.pd);

PreimageCategory:=function(b)local c,f,m,l,d,i,j;
  c:=PreimageSliding(b);l:=[1..Length(c)];
  if l=[] then Print(b," is not in the image of the sliding\n");return;fi;
  f:=function(a,b)a:=Sliding(a)[2]*Sliding(b)[2]^-1;
    if a.pd>=0 then return a;else return "-";fi;
  end;
  d:=PreimageSlidingSSS(b);
  d:=List(l,function(i)
    if c[i]in d then return SPrint(i,"*");else return i;fi;end);
  m:=List(c,a->List(c,b->f(a,b)));
  SortBy(l,x->[-Number(l,i->m[i][x]="-"),Number(m[x],y->y="-")]);
  Print("b=",b,"\n",FormatTable(m{l}{l},d,l,rec()));
  ShowMaps(Concatenation(List(l,i->List(Filtered(l,j->j<>i and m[i][j]<>"-"),
    j->[i,m[i][j],j]))));
end;

SimpleEndomorphisms:=function(b)local M,s;
  M:=b.monoid;
  s:=List(Elements(Centralizer(M.group,EltBraid(b))),M.B);
  s:=Filtered(s,x->EltBraid(x)<>M.group.identity);
  return Filtered(s,x->b=b^x);
end;

SimpleMorphisms:=function(b)local M,s;
  M:=b.monoid;
  s:=Concatenation(SimpleLeftDivisors(M,M.delta));
  s:=Filtered(s,x->x<>M.group.identity);
  s:=List(s,M.B);
  return Filtered(s,function(x)x:=b^x;return x.pd>=0;end);
end;

imtrans:=function(b)local l,s;
  l:=PreimageSliding(b);
  l:=List(l,x->[x,SimpleEndomorphisms(x)]);
  for s in l do s[2]:=List(s[2],y->[y,STransport(s[1],y)]);od;
  return [SimpleEndomorphisms(b),
    Set(Concatenation(List(l,x->List(x[2],y->y[2]))))];
end;

# preimages de f partant de a
PreimageSTransport:=function(a,f)local b,c,cp;
  b:=a^f;
  if b.pd<0 then Error(f," not a morphism from ",a);fi;
  c:=Cartesian(PreimageSliding(a),PreimageSliding(b));
  for cp in c do Add(cp,Sliding(cp[1])[2]*f*Sliding(cp[2])[2]^-1);od;
  return Filtered(c,x->x[3].pd>=0);
end;

PreimageSTransportSSS:=function(a,f)local b,c,cp;
  b:=a^f;
  if b.pd<0 then Error(f," not a morphism from ",a);fi;
  c:=Cartesian(PreimageSlidingSSS(a),PreimageSlidingSSS(b));
  for cp in c do Add(cp,Sliding(cp[1])[2]*f*Sliding(cp[2])[2]^-1);od;
  return List(Filtered(c,x->x[3].pd>=0),x->x{[1,3]});
end;
  
# preimages de f partant de a
PreimageSTransportSSS2:=function(a,f)local c;
  c:=a^f;
  if c.pd<0 then Error(f," not a morphism from ",a);fi;
  c:=PreimageSlidingSSS(a);
  return List(c,function(z)local s; s:=SPullback(z,f);
     return [z,s,STransport(z,s)];end);
end;

PreimageSTransportSSS3:=function(a,f)local c;
  c:=a^f;
  if c.pd<0 then Error(f," not a morphism from ",a);fi;
  c:=PreimageSlidingSSS(a);
  return List(c,function(z)local s; s:=SPullback2(z,f);
     return [z,s[1],STransport(z,s[1]),s[2]];end);
end;

IsInSSS:=function(b)local l;
  l:=RepresentativeSC(b).circuit[1];
  return  b.pd=l.pd and Length(b.elm)=Length(l.elm);
end;

IsInSC:=b->b in RepresentativeSC(b).circuit;

test3:=function(b)local l;
  if not IsInSSS(b) then return true;fi;
  l:=RepresentativeSC(b).circuit;
  if b in l then return true;fi;
  if Sliding(b)[1] in l then return true;fi;
  return false;
end;

RSliding:=function(y) y:=Sliding(ReversedWord(y));
  return [ReversedWord(y[1]),y[2]^-1];
end;

sr:=function(y)local z; y:=RSliding(y);z:=Sliding(y[1]);
  return [z[1],y[2]^-1*z[2]];
end;

test2:=b->IsInSSS(b) and  Length(PreimageSlidingSSS(b))>0 and b<>sr(b);

rs:=function(y)local z; y:=Sliding(y); z:=RSliding(y[1]);
  return [z[1],y[2]/z[2]];
end;

srcat:=function(b)local cat,arr,i,c;
  cat:=[b];arr:=[];i:=1;
  repeat
    c:=sr(cat[i])[1];p:=Position(cat,c);
    if p=false then Add(cat,c);Add(arr,[i,2,Length(cat)]);
    else Add(arr,[i,2,p]);
    fi;
    c:=rs(cat[i])[1];p:=Position(cat,c);
    if p=false then Add(cat,c);Add(arr,[i,3,Length(cat)]);
    else Add(arr,[i,3,p]);
    fi;
    i:=i+1;
  until i>Length(cat);
# ShowMaps(Filtered(arr,x->x[1]<>x[3]));
  return [cat,arr];
end;

# n braids of length l such that f
braidssuch:=function(n,f,M)local res,b,c;
  res:=[];c:=0;
  repeat b:=RandomBraid(M,Random([5..15]));
    if f(b) then Add(res,b);Print("#\c");
    else Print(".\c");c:=c+1;
    fi;
  until Length(res)=n;
  Print("after trying ",c," braids...\n");
  return res;
end;

tt:=function(b,f)local sb,sf;
  sb:=List(Sliding(ReversedWord(b)),ReversedWord);
  sf:=sb[2]*f*ReversedWord(Sliding(ReversedWord(b^f))[2])^-1;
  return STransport(sb[1],sf);
end;

imtrans2:=function(b)local l;
  l:=SimpleMorphisms(b);
  return [l,Filtered(l,f->Length(PreimageSTransport(b,f))>0)];
end;

imtrans2SSS:=function(b)local l;
  l:=SimpleMorphisms(b);
  return [l,Filtered(l,f->Length(PreimageSTransportSSS(b,f))>0)];
end;

imtrans3SSS:=function(b)local l,c;
  l:=SimpleMorphisms(b); c:=PreimageSlidingSSS(b);
  return [l,Filtered(l,f->ForAny(c,z->STransport(z,SPullback(z,f))=f))];
end;

comp:=function(b)local l0,l;
  l0:=imtrans2(b)[2];
  l:=List(l0,x->x*imtrans2(b^x)[2]);
  l:=Set(Concatenation(l));
  l:=Filtered(l,x->(x.pd=0 and Length(x.elm)=1)or(x.pd=1 and Length(x.elm)=0));
  return [Difference(l0,l),Difference(l,l0)];
end;

compSSS:=function(b)local l0,l;
  l0:=imtrans2SSS(b)[2];
  l:=List(l0,x->x*imtrans2SSS(b^x)[2]);
  l:=Set(Concatenation(l));
  l:=Filtered(l,x->(x.pd=0 and Length(x.elm)=1)or(x.pd=1 and Length(x.elm)=0));
  return [Difference(l0,l),Difference(l,l0)];
end;

compSSS2:=function(b)local l0,l;
  l0:=imtrans3SSS(b)[2];
  l:=List(l0,x->x*imtrans3SSS(b^x)[2]);
  l:=Set(Concatenation(l));
  l:=Filtered(l,x->(x.pd=0 and Length(x.elm)=1)or(x.pd=1 and Length(x.elm)=0));
  return [Difference(l0,l),Difference(l,l0)];
end;

PreimageCategorySSS:=function(b)local c,f,m,l;
  c:=PreimageSlidingSSS(b);
  if Length(c)=0 then Print(b," not image of sliding in SSS\n");return;fi;
  f:=function(a,b)a:=Sliding(a)[2]*Sliding(b)[2]^-1;
    if a.pd>=0 then return a;else return "-";fi;end;
  m:=List(c,a->List(c,b->f(a,b)));l:=[1..Length(c)];
  SortBy(l,x->Number(m[x],y->y="-"));
  Print("b=",b,"\n");
  Print(FormatTable(m{l}{l},l,l,rec()));
  ShowMaps(Concatenation(List(l,i->List(Filtered(l,j->j<>i and m[i][j]<>"-"),
    j->[i,m[i][j],j]))));
end;

findcl:=function(l,cat)local res,l,cl;
  res:=[];
  while Length(l)>0 do
  cl:=ConjugacySet(l[1],"Pos");
  Add(res,[Length(cl),Length(cat(cl[1]))]);
  l:=Difference(l,cl);
  od;
  Sort(res);
  return res;
end;
