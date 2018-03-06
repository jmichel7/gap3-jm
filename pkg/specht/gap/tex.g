########################################################################
## TeX()  GAP routines for TeXing (polynomials) GAP objects           ##
##                                                                    ##
##        In the case of records, the general handler TeX() checks    ##
##        for an <X>.operations.TeX() command and calls this function ##
##        if possible. Otherwise special functions are provided here  ##
##        only for polynomials and (trivially) matrices (actually     ## 
##        LaTeX code).                                                ##
##                                                                    ##
##        These routines are included with the SPECHT package in      ##
##        order to allow "crystalized decomposition matrices" (which  ##
##        have polynomial entries), to be TeXed.                      ##
##                                                                    ##
##        (So not much in the way of documentation...)                ##
##                                                                    ##
##        Andrew Mathas                                               ##
##                                                                    ##
########################################################################

#D Generic functions for producing TeX code for Gap elements

# July 2003
#  o changed the function so that it returns a string
## July 1997
##   o fixed minor bug TeXing lists.

## We define TeX() below; for now we need to know that it exists.
TeX:=Ignore;

# return the catenation of the String()s of all of the arguments
AddToString:=function(arg)
  local tex, strings;
  strings:=List(arg, String);
  return ApplyFunc(ConcatenationString, strings);
end;

#P This function does the real work in TeXing polynomials. The functions 
## may not work well for polynomials with non-rational coefficients; I've 
## never checked. 
TeXPolynomial:=function(poly)
  local tex, i, val, coeffs, polyname, TeXCoeff;
 
  # TeXing the term a_n x^n where
  #   a_n = coeff
  #   n = val[uation]
  #   x = polyname
  # and len is the number of terms in the polynomials being TeXed.
  TeXCoeff:=function(coeff, polyname, val, len)local tex;
    tex:="";
    if not IsRat(coeff) then 
      tex:=TeXPolynomial(coeff);
      if val<> 0 then
        if val = 1 then val:="";
        else val:=Concatenation("^{", String(val), "}");
        fi;
        PrintToString(tex,polyname, val);
      fi;
    elif coeff <> 0 then
      if val = 0 then
        if (coeff <> 1 and coeff <> -1) or len > 1 then
          PrintToString(tex,coeff);
        elif coeff = -1 then PrintToString(tex,"-1");
        elif coeff = 1 then PrintToString(tex,"1");
        fi;
      else
        if val = 1 then val:="";
        else val:=Concatenation("^{", String(val), "}");
        fi;
        if coeff = -1 then PrintToString(tex,"-", polyname, val);
        elif coeff = 1 then PrintToString(tex,polyname, val);
        else PrintToString(tex,coeff, polyname, val);
        fi;
      fi;
    fi;
    return tex;
  end; # TeXCoeff
  
  if not IsPolynomial(poly) then return TeX(poly);
  elif Length(poly.coefficients) = 0 then return "0"; # poly = 0
  else
    if not IsBound(poly.baseRing.indeterminate.name) then polyname:="x";
    else polyname:=poly.baseRing.indeterminate.name;  
    fi;
    coeffs:=Copy(poly.coefficients);
    coeffs:=coeffs{[Length(coeffs),Length(coeffs)-1..1]};
    val:=poly.valuation+Length(coeffs)-1;
    tex:=TeXCoeff(coeffs[1], polyname, val, Length(coeffs));
    if Length(coeffs) > 1 then
      for i in [2..Length(coeffs)] do
        if IsRat(coeffs[1]) and coeffs[i] > 0 then PrintToString(tex,"+"); fi;
        val:=val - 1;
	PrintToString(tex,TeXCoeff(coeffs[i], polyname, val, Length(coeffs)));
      od;
    fi;
    return tex;
  fi;
end;      # __TeXPolynomial

#P LaTeX code for a matrix
TeXMatrix:=function(M)
  local tex, i, j;

  tex:=SPrint("\\left(\\begin{array}{*{", Length(M[1]),"}{l}}\n");
  for i in [1..Length(M)] do
    for j in [1..(Length(M[i])-1)] do
      PrintToString(tex,TeX(M[i][j]),"& ");
      if ( j mod 10 = 0 ) then PrintToString(tex,"\n"); fi;
    od;
    PrintToString(tex,TeX(M[i][Length(M[i])]),"\\\\\n");
  od;
  PrintToString(tex,"\\end{array}\\right)\n");
  return tex;
end;

#P Usage: TeXWideMatrix(m1 [,m2, ...], C)
##   C is the number of columns per page. All of the matrices given to 
## TeXWideMatrix() are printed on the same page (if this is possible).
## *** undocumented
TeXWideMatrix:=function(arg)
  local tex, mats, m, C, r, c, col, cmax, arraytop;

  mats:=arg{[1..Length(arg)-1]};  ## matrices
  if mats=[] or not IsInt(arg[Length(arg)])  then
    Error("usage: TeXWideMatrix(m1 [,m2...], col)\n");
    return;
  fi;
  C:=arg[Length(arg)];            ## columns per page

  cmax:=Maximum(List(mats, m->Maximum(List(m, r->Length(r)))));
  arraytop:=Concatenation("\\begin{array}{*{",String(cmax),"}{l}}\n");
  col:=[1-C..0];
  tex:="";
  while col[C]<cmax do
    col:=col+C;
    for m in mats do
      if col[1]=1 then PrintToString(tex,"$\\left(");
      else PrintToString(tex,"$\\left.");
      fi;
      PrintToString(tex,arraytop);
      for r in [1..Length(m)] do
        if IsBound(m[r][col[1]]) then 
          TeX(m[r][col[1]]); 
          for c in col{[2..C]} do
            if c<=Length(m[r]) then TeX("& ",m[r][c]); fi;
          od;
        fi;
        PrintToString(tex,"\\\\\n");
      od;
      if col[C]>=cmax then PrintToString(tex,"\\end{array}\\right)$\n\n\n");
      else PrintToString(tex,"\\end{array}\\right.$\n\n\\medskip\n\n");
      fi;
    od;
    if col[C]<cmax then PrintToString(tex,"\\newpage\n\n");fi;
  od;
  return tex;
end;


#P The TeX() handler; calls the above routines and checks operations record
TeX:=function(arg) local tex, a, i;
  tex:="";
  for a in arg do
    if IsRat(a) or IsString(a) or IsPerm(a) then PrintToString(tex, a);
    elif IsMatrix(a) then PrintToString(tex, TeXMatrix(a));
    elif IsList(a) then 
      if a=[] then PrintToString(tex, "[]");
      else
	tex:=SPrint("[ ", TeX(a[1]));
        for i in [2..Length(a)] do 
	  PrintToString(tex,", ", TeX(a[i])); 
        od;
        PrintToString(tex," ]");
      fi;
    elif IsPolynomial(a) then PrintToString(tex, TeXPolynomial(a));
    elif IsRec(a) and IsBound(a.operations) and IsBound(a.operations.TeX)
      then PrintToString(tex, a.operations.TeX(a));
    elif a=true then tex:="T";
    elif a=false then tex:="F";
    else Error("TeX(<a>), don't know how to TeX <a>.\n");
    fi;
  od;
  return tex;
end;

#U TeXTop(heading)                                                            
#U TeXTop(heading[,preamble1,preamble2,...])                                  
#M Prints a reasonable preamable for a tex file.                              
TeXTop:=function(arg) local heading, str, p;
  heading:=arg[1];
  str:="\\documentclass{amsart}\n\\usepackage{fullpage}\n";
  for p in arg{[2..Length(arg)]} do
    PrintToString(str,p,"\n");
  od;
  PrintToString(str,"\\nofiles\n\\markboth{",heading,"}{",heading,"}\n");
  PrintToString(str,"\\begin{document}\n\n");
  return str;
end;
