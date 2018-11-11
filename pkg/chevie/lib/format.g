#############################################################################
##
#A  format.g              CHEVIE library                          Jean Michel
##
#Y  Copyright (C) 1999 - 2006  University Paris VII.
##
##  This file contains routines to format various objects in various ways.
##  In particular in provides general support for outputting objects in
##  TeX or LaTeX.
##

#############################################################################
##
#F  IntListToString( <part>, [<brackets>] ) . . . . . . . . as the name says
##  
##  <part> must be a list of integers. If all of them are in [0..9] then
##  a  string  of  digits  corresponding  to  the  entries  of <part> is
##  returned. Otherwise the elements of <part> are converted to strings,
##  concatenated  with separating  commas and  the result  surrounded by
##  brackets (default brackets are parentheses):
##  
##     IntListToString([4,2,2,1,1])      ="42211"
##     IntListToString([14,2,2,1,1])     ="(14,2,2,1,1)"
##     IntListToString([14,2,2,1,1],"[]")="[14,2,2,1,1]"
##  
IntListToString:=function(arg) local digits,brackets,p;
  p:=arg[1];
  if Length(arg)=2 then brackets:=arg[2];else brackets:="()";fi;
  if ForAll(p,x->x<10 and x>=0) then digits:="0123456789"; 
     return String(digits{1+p});
  else return ConcatenationString([brackets[1]],Join(p),[brackets[2]]);
  fi;
end;

#############################################################################
##
#F  FormatTable( <object>, [<option>] ) . . . . . . . .  as the name says
##
##  General routine to format a table.
##  Options:
##     rowLabels           Labels for rows
##     columnLabels        Labels for columns
##     rowsLabel           Label for column of rowLabels
##     separators          line numbers after which put a separator
##     columnRepartition   cut in vertical pieces by these column numbers
##     lines               show only these lines
##     columns             show only these columns
##     TeX
##     LaTeX
##
FormatTable:=function(t,option)
  local rowlab,collab,colrep,rowslab,sep,res,pad,i,j,cols,lastcol,fsep,
        t,colwidth,labwidth,cut,ci,res,l,hline,k;
  rowlab:=List(option.rowLabels,x->Format(x,option));
  t:=List(t,x->List(x,y->Format(y,option))); 
  if IsBound(option.rowsLabel) then rowslab:=option.rowsLabel;else rowslab:="";
  fi;
  if IsBound(option.separators) then sep:=option.separators;else sep:=[0];fi;
  if IsBound(option.columnLabels) then 
    collab:=List(option.columnLabels,x->Format(x,option));fi;
  if IsBound(option.rows) then 
    rowlab:=rowlab{option.rows};t:=t{option.rows};
  fi;
  if IsBound(option.columns) then 
    collab:=collab{option.columns};t:=List(t,x->x{option.columns});
  fi;
  if IsBound(option.TeX) then
    if IsBound(option.columnRepartition) then colrep:=option.columnRepartition;
    else colrep:=[Length(t[1])];
    fi;
    if IsBound(option.LaTeX) then fsep:="\\\\";else fsep:="\\cr";fi;
    lastcol:=0;res:="";
    for j in colrep do
      cols:=lastcol+[1..j];lastcol:=lastcol+j;
      if IsBound(option.LaTeX) then
	Append(res,"\\tablehead{\\hline ");
      else
	Append(res,"{\\tabskip 0mm");pad:="\\kern 1mm";
	PrintToString(res,"\\halign{\\strut\\vrule",pad,"\\hfill$#$",pad,
	  "\\vrule&");
	PrintToString(res,Join(List(cols,i->SPrint(pad,"\\hfill$#$",pad)),"&"),
	  "\\vrule",fsep,"\n");
	Append(res,"\\noalign{\\hrule}\n");
      fi;
      PrintToString(res,rowslab);
      PrintToString(res,"&",Join(collab{cols},"&"),fsep);
      if IsBound(option.LaTeX) then
        if 0 in sep then Append(res,"\\hline}\n");fi;
	Append(res,"\\tabletail{\\hline}\n");
	PrintToString(res,"\\begin{supertabular}{|R|",List(cols,i->'R'),"|}\n");
	Append(res,"\\shrinkheight{30pt}\n");
      else
	Append(res,"\n");
	if 0 in sep then Append(res,"\\noalign{\\hrule}\n");fi;
      fi;
      for i in [1..Length(t)] do 
	PrintToString(res,rowlab[i],"&",Join(t[i]{cols},"&"),fsep,"\n");
	if i in sep then 
	  if IsBound(option.LaTeX) then Append(res,"\\hline\n");
	                           else Append(res,"\\noalign{\\hrule}\n");
	  fi;
	fi;
      od;
      if IsBound(option.LaTeX) then Append(res,"\\end{supertabular}\n");
      else Append(res,"\\noalign{\\hrule}}\n}\n");
      fi;
    od;
  else 
    colwidth:=List(TransposedMat(t),x->Maximum(List(x,Length)));
    if IsBound(option.columnLabels) then 
      colwidth:=Zip(colwidth,collab,function(x,y)return Maximum(x,Length(y));end);
      collab:=Zip(collab,colwidth,function(x,y)return String(x,y);end);
    fi;
    labwidth:=Maximum(Length(rowslab),Maximum(List(rowlab,Length)));
    rowslab:=String(rowslab,labwidth);
    rowlab:=List(rowlab,x->String(x,-labwidth));
    hline:=function(ci)
      PrintToString(res,List([1..Sum(colwidth{ci}+1)+labwidth+2],i->'_'),"\n");
    end;
    cut:=function(l,max)local res,i,len,n;# cut l in parts of sum<max
      res:=[];len:=0;n:=0;
      for i in l do len:=len+i;
	if len>=max then 
	  if n=0 then Add(res,1);n:=0;len:=0;
	  else Add(res,n);n:=1;len:=i; 
	  fi;
        else n:=n+1; fi;
      od;
      Add(res,n);return res;
    end;
    if IsBound(option.columnRepartition) then colrep:=option.columnRepartition;
    elif IsBound(option.screenColumns) then
      colrep:=cut(1+colwidth,option.screenColumns-labwidth-3);
    else colrep:=[Length(t[1])];
    fi;
    ci:=[0];
    res:=[];
    for  k in colrep do
      ci:=ci[Length(ci)]+[1..k];
      if IsBound(option.columnLabels) then
	PrintToString(res,rowslab," |",Join(collab{ci}," "),"\n");
	if 0 in sep then hline(ci);fi;
      fi;
      for l in [1 .. Length(t)]  do
	PrintToString(res,rowlab[l]," |",
	  Join(List(ci,i->String(t[l][i],colwidth[i]))," "),"\n");
	if l in sep then hline(ci);fi;
      od;
      if ci[Length(ci)]<>Length(colwidth) then Append(res,"\n");fi;
    od;
  fi;
  return res;
end;

# compute a descriptive string from a TeX string (for names, etc...)
TeXStrip:=function(arg)local s,opt,f;
  f:=function(s)local p,i;
    p:=PositionSublist(s,"ph");
    if p=false then return s;fi;
    i:=p+2; while i<=Length(s) and s[i] in "0123456789," do i:=i+1;od;
    return Concatenation(s{[1..p-1]},"phi{",s{[p+2..i-1]},"}",
         f(s{[i..Length(s)]}));
  end;
  s:=arg[1];if Length(arg)=1 then opt:=rec();else opt:=arg[2];fi;
  if not IsBound(opt.TeX) then
    s:=Replace(s,
      "_","",
      "$","",
      "\\tilde ","~",
      "{","",
      "}","",
      "\\phi","ph",
      "\\otimes",",",
      "\\cdot","",
      "\\zeta","E",
      "\\varepsilon","eps",
      "\\gamma","gamma",
      "\\lambda","lambda",
      "\\theta","theta",
      "\\nu","nu",
      "\\rho","r",
      "\\sigma","sigma",
      "\\BZ","Z",
      "\\!\\wedge\\!",".",
      "\\wedge","",
      "\\overline","'",
      " ","");
    s:=List([1..Length(s)],function(i)if s[i]='i' and 
      (i=Length(s) or not s[i+1] in "abcdefghijklmnopqrstuvwxyz")
      then return 'I';else return s[i];fi;end);
    if Length(s)>0 and s[1]='^' then s:=s{[2..Length(s)]};fi;
    s:=String(f(s));
  fi;
  return Copy(s);
end;

TeXBracket:=function(s)
  s:=String(s);
  if Length(s)=1 then 
    if s[1] in "abcdefghijklmnopqrstuvwxyABCDEFGHIJKLMNOPQRSTUVWXYZ" then
         return Concatenation(" ",s);
    else return s;
    fi;
  else return Concatenation("{",s,"}");
  fi;
end;

# Formats name^p
FormatMonomial:=function(name,p,option)local res;
  p:=String(p);
  if p="0" then return "1";fi;
  if IsBound(option.TeX) then p:=TeXBracket(p);
  elif IsBound(option.Maple) and p[1]='-' then p:=SPrint("(",p,")");fi;
  if p<>"1" then return SPrint(name,"^",BracketIfNeeded(p,"+-/*"));
  else return String(name);fi;
end;

# if in an expression like 2*a brackets are needed around a returns (a) 
# Brackets are needed if <a> contains one of <ops> at bracket level 0
# Typically <ops> are "+-" . In some contexts could be "*/"
BracketIfNeeded:=function(a,ops)local par,bracket,i;
  par:=0;bracket:=false;
  for i in [1..Length(a)] do
    if a[i] in "({" then par:=par+1;
    elif a[i] in ")}" then par:=par-1;
    fi;
    if i>1 and a[i] in ops and a[i-1]<>'^' and par=0 then bracket:=true;fi;
  od;
  if bracket then return Concatenation("(",a,")");
  else return a;
  fi;
end;

# formats coeff*p
FormatCoefficient:=function(coeff,p,option)
  if not IsString(coeff) then coeff:=Format(coeff,option);fi;
  if coeff="0" then return "0";fi;
  if coeff="1" and Length(p)<>0 then return p;fi;
  if coeff="-1" and Length(p)<>0 then return Concatenation("-",p);fi;
  if p="1" or p="" then return coeff;fi;
  coeff:=BracketIfNeeded(coeff,"+-");
  if (IsBound(option.GAP) or IsBound(option.Maple)) and not IsBound(option.TeX)
  then Append(coeff,"*");fi;
  return ConcatenationString(coeff,p);
end;

FormatPolynomial:=function(coeffs,valuation,name,option)local res,x,p;
  res:="";
  for x in [Length(coeffs),Length(coeffs)-1..1] do
    p:=FormatCoefficient(coeffs[x],FormatMonomial(name,x+valuation-1,option),option);
    if p<>"0" then
      if p[1]<>'-' then Add(res,'+');fi;
      Append(res,p);
    fi;
  od;
  if Length(res)=0 then res:="0";fi;
  if res[1]='+' then res:=res{[2..Length(res)]};fi;
  return String(res);
end;

# q/d for a quadratic q
FormatQuadratic:=function(q,d,option)local res,gcd,c,p;
  c:=[q.a,q.b];gcd:=Gcd(c);c:=c/gcd;
  res:=rec(s:="",d:=gcd/(q.d*d));
  c:=c*Numerator(res.d);res.d:=Denominator(res.d);
  if c[1]<>0 then Append(res.s,String(c[1])); fi;
  if IsBound(option.TeX) then p:=SPrint("\\sqrt ",TeXBracket(q.root));
  elif IsBound(option.Maple) then p:=SPrint("sqrt(",q.root,")");
  else p:=SPrint("ER(",q.root,")");
  fi;
  p:=FormatCoefficient(c[2],p,option);
  if p[1]<>'-' and res.s<>"" then Add(res.s,'+');fi;
  Append(res.s,p);
  return res;
end;

# a/b
FormatQuotient:=function(a,b,option)
  b:=Format(b,option);
  if b="1" then return a;fi;
  if IsBound(option.LaTeX) then return SPrint("\\frac",TeXBracket(a),TeXBracket(b));
  elif IsBound(option.TeX) then return SPrint("{",a,"\\over",b,"}");
  else return SPrint(BracketIfNeeded(a,"+-/*"),Concatenation("/",b));
  fi;
end;

FormatCyclotomic:=function(c,option)
  local N,root,v,d,q,res,x,try,FormatRootUnity;
  FormatRootUnity:=function(x)local res;res:="";
    if (Denominator(x) mod 4)=2 or x=3/4 then x:=Mod1(x+1/2);res:="-";fi;
    return SPrint(res,FormatMonomial(root(Denominator(x)),Numerator(x),option));
  end;
  N:=NofCyc(c); res:="";
  root:=function(N)
    if IsBound(option.TeX) then 
      if N=4 then return "i";fi;
      return Concatenation("\\zeta_",TeXBracket(N));
    elif IsBound(option.GAP) then return Concatenation("E(",String(N),")");
    else if N=4 then return "I";fi;
      return Concatenation("E",String(N));
    fi;
  end;
  v:=COEFFSCYC(c);
  d:=Lcm(List(v,Denominator));
  if not N in [1,4] and Number(v,x->x<>0)>1 and not IsBound(option.noQuadrat) then
    q:=Quadratic(c*d);
    if q<>false then res:=FormatQuadratic(q,d,option);
      return FormatQuotient(res.s,res.d,option);
    fi;
    q:=AsRootOfUnity(c*d);
    if q<>false then
      return FormatQuotient(FormatRootUnity(q),d,option);
    fi;
    try:=[(1+E(4))/2,(1-E(4))/2,E(3)/2,E(3)^2/2,1+E(3),1-E(3),
           1+E(3)^2,1-E(3)^2,ER(-3)/2,E(3)^2/ER(-3),E(3)^2*(1+E(4))/2];
    for x in try do # hack for unip. deg. coeffs in G_6, G_14, G_27
      q:=Quadratic(c*d/x);
      if q<>false then 
        res:=BracketIfNeeded(FormatCyclotomic(2*x,option),"+-");
	q:=FormatQuadratic(q,2*d,option);
	return FormatQuotient(FormatCoefficient(q.s,res,option),q.d,option);
      fi;
    od;
  fi;
  x:=AsRootOfUnity(c); if x<>false then return FormatRootUnity(x); fi;
  return FormatQuotient(FormatPolynomial(d*v,0,root(N),option),d,option);
end;

# List function, but accepts lists with unbound entries
ListUnbnd:=function(l,f)local res,i;
  res:=[];
  for i in [1..Length(l)] do
    if IsBound(l[i]) then res[i]:=f(l[i]);fi;
  od;
  return res;
end;

# display an object of certain types I recognize:
# scalars,vectors,matrices of cyclotomics or multi-variable
# polynomials over the cyclotomics
Format:=function(arg)local obj,v,res,s,i,m,option,l;
  obj:=arg[1];
  if Length(arg)=1 then option:=rec();else option:=arg[2];fi;
  if IsString(obj) then 
    if IsBound(option.GAP) then return SPrint("\"",Replace(obj,
	  "\\","\\\\", "\"","\\\"", "\n","\\n"),"\"");
    else return obj;
    fi;
  elif IsList(obj) then
    if IsBound(option.GAP) or IsBound(option.Maple) then 
         return SPrint("[",Join(List([1..Length(obj)],
	   function(i)if IsBound(obj[i]) then return Format(obj[i],option);
	     else return "";fi;end)),"]");
    elif ForAll(obj,IsList) and Length(Set(ListUnbnd(obj,Length)))<=1 then 
    # IsMat
      m:=ListUnbnd(obj,v->List(v,x->Format(x,option)));
      if IsBound(option.LaTeX) then 
        return Join(ListUnbnd(m,v->SPrint(Join(v,"&"),"\\\\\n")),"");
      elif IsBound(option.TeX) then 
        return Join(ListUnbnd(m,v->SPrint(Join(v,"&"),"\\cr\n")),"");
      else l:=Set(ListUnbnd(obj,Length))[1];
	s:=List([1..l],i->Maximum(Set(ListUnbnd(m,x->Length(x[i])))));
        return Join(List([1..Length(m)],function(v)
	  if IsBound(m[v]) then
	    return Join(List([1..l],i->String(m[v][i],s[i]))," ");
	  else return "";
	  fi;end),"\n");
      fi;
    else return SPrint(Join(List(obj,x->Format(x,option))," "));
    fi;
  elif IsInt(obj) then return String(obj);
  elif IsRat(obj) then return FormatQuotient(Format(Numerator(obj),option),
                                              Denominator(obj),option);
  elif IsCyc(obj) then return FormatCyclotomic(obj,option);
  elif IsPolynomial(obj) then 
    if IsBound(obj.name) then
      return FormatPolynomial(obj.coefficients,obj.valuation,obj.name,option);
    elif IsBound(obj.baseRing.indeterminate) then
        return FormatPolynomial(obj.coefficients,obj.valuation,
        String(obj.baseRing.indeterminate),option);
    else return FormatPolynomial(obj.coefficients,obj.valuation,"x",option);
    fi;
  elif IsBool(obj) then 
    if IsBound(option.GAP) then 
         if obj then return "true"; else return "false"; fi;
    else if obj then return "T"; else return "F"; fi;
    fi;
  elif IsUnknown(obj) then if IsBound(option.GAP) then return("Unknown()");
                           else return "?";fi;
  elif IsRec(obj) then
    if IsBound(obj.operations) then
      if IsBound(obj.operations.Format)
      then return obj.operations.Format(obj,option);
      elif IsBound(obj.operations.String)
      then return obj.operations.String(obj);
      fi;
    fi;
    return SPrint("rec(",Join(List(RecFields(obj),
        f->SPrint(f,":=",Format(obj.(f),option))),",\n"),")");
  elif IsString([obj]) then return [obj];
  else return String(obj);
  fi;
end;

FormatGAP:=function(arg)local option;
  if Length(arg)=2 then option:=arg[2];else option:=rec();fi;
  option.GAP:=true;
  return Format(arg[1],option);
end;

FormatMaple:=function(arg)local option;
  if Length(arg)=2 then option:=arg[2];else option:=rec();fi;
  option.Maple:=true;
  return Format(arg[1],option);
end;

FormatTeX:=function(arg)local option;
  if Length(arg)=2 then option:=arg[2];else option:=rec();fi;
  option.TeX:=true;
  return Format(arg[1],option);
end;

FormatLaTeX:=function(arg)local option;
  if Length(arg)=2 then option:=arg[2];else option:=rec();fi;
  option.TeX:=true;option.LaTeX:=true;
  return Format(arg[1],option);
end;
