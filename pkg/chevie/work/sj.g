SimplifyCycDecompositionJ:=function(r)local i,p,p1,v,tr;
  tr:=function(arg)local sterm;sterm:=v->SPrint("P",v.cyc,"(",v.monomial,")");
    Print(sterm(r.vcyc[i]),Join(List(r.vcyc{[3..Length(arg)]},sterm),""));
    r.vcyc[i]:=rec(monomial:=arg[1],cyc:=arg[2]);
    Print("=>",sterm(r.vcyc[i]),"\n");
    r.vcyc:=Drop(r.vcyc,arg{[3..Length(arg)]});
    i:=1;
  end;
  i:=1;
#   SortBy(r.vcyc,v->[v.cyc,v.monomial.elm[1]]);
  while i<=Length(r.vcyc) do
    v:=r.vcyc[i];
    p:=Position(r.vcyc,rec(monomial:=E(3)*v.monomial,cyc:=v.cyc));
    if p<>false then
      p1:=Position(r.vcyc,rec(monomial:=E(3)^2*v.monomial,cyc:=v.cyc));
      if p1<>false and Gcd(v.cyc,3)=1 then 
        # Product([0..2],i->Pd(E3^ix))=Pd(x^3) if Gcd(d,3)=1
        tr(v.monomial^3,v.cyc,p,p1);
      elif v.cyc in [1,2] then # Pd(x)Pd(E3x)=P3d(x) if d in [1,2]
        tr(v.monomial*E(3)^2,3*v.cyc,p);
      fi;
    fi;
    p:=Position(r.vcyc,rec(monomial:=-v.monomial,cyc:=v.cyc));
    if p<>false and Gcd(v.cyc,2)=1 then #Pd(x)Pd(-x)=(-1)^phi(d)Pd(x^2) d odd
      r.factor:=(-1)^Phi(v.cyc)*r.factor; tr(v.monomial^2,v.cyc,p);
    elif p<>false and v.cyc=2 then #P2(x)P2(-x)=-P1(x^2)
      r.factor:=-r.factor; tr(v.monomial^2,1,p);
    fi;
    p:=Position(r.vcyc,rec(monomial:=v.monomial,cyc:=2*v.cyc));
    if p<>false and Gcd(v.cyc,2)=1 then #Pd(x)P2d(x)=Pd(x^2)
      tr(v.monomial^2,v.cyc,p);
    fi;
    p:=Position(r.vcyc,rec(monomial:=E(3)/v.monomial,cyc:=v.cyc));
    if p<>false and v.cyc=1 then
      r.factor:=r.factor*-E(3)/v.monomial; tr(E(3)*v.monomial,3,p);
    fi;
    if IsMvp(v.monomial) and v.monomial.coeff[1]=-1 then
      if Gcd(v.cyc,2)=1 then
        if v.cyc=1 then r.factor:=-r.factor;fi; tr(-v.monomial,2*v.cyc);
      elif Gcd(v.cyc,4)=2 then
        if v.cyc=2 then r.factor:=-r.factor;fi; tr(-v.monomial,v.cyc/2);
      else v.monomial:=-v.monomial;
      fi;
    fi;
    i:=i+1;
  od;
  return r;
end;

