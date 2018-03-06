CheckUnipotentCentralizers:=function(arg)local W,p,c,t,cl,cc,q,w,sum,f,u,g,N;
  W:=arg[1];
  if Length(arg)=1 then p:=0;else p:=arg[2];fi;
  t:=RationalUnipotentClasses(W,p);q:=Mvp("q");sum:=0;
  u:=UnipotentClasses(W,p);
  g:=CycPol(GenericOrder(W,q));
  for c in t do
    cl:=u.classes[c.classno];
    cc:=q^cl.dimunip;
    if NrConjugacyClasses(cl.Au)>1 then
      if cl.red.rank>0 then
        if IsBound(cl.F) then 
          if IsPerm(cl.F) then f:=MatXPerm(cl.red,cl.F);else f:=cl.F;fi;
        else f:=MatXPerm(cl.red,());
        fi;
        if IsBound(cl.AuAction) and ForAny(cl.AuAction.F0s,x->x<>x^0) then
          w:=ChevieClassInfo(cl.Au).classtext[c.AuNo];
          w:=Product(cl.AuAction.F0s{w})*f;
          if w=1 then w:=cl.AuAction.F0s[1]^0;fi;
          cc:=cc*GenericOrder(Spets(cl.red,w),q);
        else
#   Print("cl=",cl," red=",cl.red," F=",cl.F,"\n");
           cc:=cc*GenericOrder(Spets(cl.red,f),q);
        fi;
      fi;
      cc:=cc*Size(cl.Au)/ChevieClassInfo(cl.Au).classes[c.AuNo];
    elif IsBound(cl.F) then 
      cc:=cc*GenericOrder(Spets(cl.red,cl.F),q);
    else cc:=cc*GenericOrder(cl.red,q);
    fi;
    sum:=sum+Value(c.card,q);
    cc:=CycPol(cc);
    if cc<>g/c.card then 
     Print(TeXStrip(cl.name));
     Print(".",c.AuNo);
     Print(":",UnipotentClassOps.FormatCentralizer(cl,rec()),
       " => ",cc," (is ",g/c.card,")\n");
    fi;
  od;
  if IsSpets(W) then N:=Group(W).N;else N:=W.N;fi;
  if sum<>q^(2*N) then
    Print("found nr. unip=",sum," instead of ",q^(2*N),"\n");
  fi;
end;

CheckNrSemisimple:=function(WF)local l,g,q;
  q:=Mvp("q");g:=GenericOrder(WF,q);
  l:=ClassTypes(WF);NrConjugacyClasses(l);
  l:=List(l.ss,x->x.nrClasses*g/Value(x.cent,q)*q^(2*Group(x.CGs).N));
  if Sum(l) <> g then
    Print("found nr elem=",Sum(l)," instead of ",g,"\n");
    Error();
  fi;
end;

CheckNrSemisimple2:=function(WF)local l,a,g,q;
  q:=X(Cyclotomics);g:=CycPol(GenericOrder(WF,q));
  l:=ClassTypes(WF);NrConjugacyClasses(l);
  a:=Sum(l,x->x.nrClasses*Value(g/x.cent,q));
  return CycPol(a);
end;
