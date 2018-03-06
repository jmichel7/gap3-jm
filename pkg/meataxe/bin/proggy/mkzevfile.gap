###############################################################################
##
#F FFList := function(field)
##
## This program was originally written by M. Geck.
## For a given finite field it returns a list giving the correspondence
## between the MeatAxe numbering and the Gap numbering of the finite field
## elements.
##
FFList := function(field)
local char,deg,elms,list;
char:=Factors(Size(field))[1];
deg:=Length(Factors(Size(field)));
elms:=List(Cartesian(List([1..deg],i->[0 .. char-1])),Reversed);
list:=List(elms,i->Sum(List([1..deg],j->i[j]*field.root^(j-1))));
return list;
end;

###############################################################################
##
#F ConvPol:=function(polynomial)
##
## Given a Gap polynomial this program returns this polynomial in reversed 
## order with coefficients in the MeatAxe numbering.
##
ConvPol:=function(polynomial)
local ffl,coeffs;
ffl:=FFList(polynomial.baseRing);
coeffs:=Reversed(List(polynomial.coefficients,i->Position(ffl,i)-1));
return coeffs;
end;    

###############################################################################
##
#F NumberPol:=function(polynomial)
##
## Given a Gap polynomial, this program returns the number for the ZSM MW option
## to compute this polynomial.
##
NumberPol:=function(polynomial)
local coeffs,number;
coeffs:=Reversed(ConvPol(polynomial));
number:=Sum(List([0 .. Length(polynomial.coefficients)-1],
                 i->Size(polynomial.baseRing)^(2^i-1)*coeffs[i+1]));
return number;
end;    

###############################################################################
##
#F BrauerValueSmall:=function(polynomial)
##
## Given an irreducible Gap polynomial, this program return the corresponding
## Brauer character value using a splitting field of the given polynomial.
##
BrauerValueSmall:=function(polynomial)
local fac,extfield,pol,charvalue;
extfield:=GF(Size(polynomial.baseRing)^(Length(polynomial.coefficients)-1));
pol:=Copy(polynomial);
pol.baseRing:=extfield;
charvalue:=0;
for fac in Factors(pol) do
    charvalue:=charvalue+(-1)^Size(pol.baseRing)
              *E(Size(extfield)-1)^LogFFE(fac.coefficients[1]);
od;
return charvalue;
end;

###############################################################################
##
#F BrauerValueBig:=function(polynomial,elmorder,conwaypolynomial)
##
## Given an irreducible Gap polynomial, this program return the corresponding
## Brauer character value using a suitable Conway polynomial.
##
BrauerValueBig:=function(polynomial,elmorder,conwaypolynomial)
local field,nullpolynomial,xpolynomial,root,charvalue,exp,power,result,j;
field:=polynomial.baseRing;
nullpolynomial:=Polynomial(field,[0]*field.one);
xpolynomial:=Polynomial(field,[0,1]*field.one);
root:=PowerMod(xpolynomial,
        (Size(field)^(Length(conwaypolynomial.coefficients)-1)-1)/elmorder,
         conwaypolynomial);
charvalue:=0;
for exp in [0 .. elmorder-1] do
    power:=PowerMod(root,exp,conwaypolynomial);
    result:=nullpolynomial;
    for j in [1 .. Length(polynomial.coefficients)] do
        result:=result+PowerMod(power,j-1,conwaypolynomial)
                      *polynomial.coefficients[j];
    od;
    result:=result mod conwaypolynomial;
    if result = nullpolynomial then
       charvalue:=charvalue+E(elmorder)^exp;
    fi;
od;
return charvalue;
end;

###############################################################################
##
#F BrauerChar:=function(field,elmorder[,conwaypolynomial])
##
## Given a finite field and an element order, this program returns the
## possible polynomials in MeatAxe notation and the corresponding
## Brauer character values. If the necessary finite fields are not available
## in GAP, a suitable Conway polynomial must be supplied.
##
BrauerChar:=function(arg)
local xpolynomial,list,factors,i;
xpolynomial:=Polynomial(arg[1],arg[1].one*[0,1]);
factors:=Factors(xpolynomial^arg[2]-xpolynomial^0); 
list:=[];
for i in [1..Length(factors)] do
    list[i]:=[];
    list[i][2]:=ConvPol(factors[i]);
    if Length(arg)=2 then
       list[i][1]:=BrauerValueSmall(factors[i]);
    fi;
    if Length(arg)=3 then
       list[i][1]:=BrauerValueBig(factors[i],arg[2],arg[3]);
    fi;
od;
return list;
end;

###############################################################################
##
#F MakeZEVFile:=function(field,elmorder[,conwaypolynomial])
##
## Given a finite field and an element order, this program writes an input
## file for ZEV. If the necessary finite fields are not available in GAP, 
## a suitable Conway polynomial must be supplied.
##
MakeZEVFile:=function(arg)
local list,name,entry,i;
name:=ConcatenationString(String(Size(arg[1])),"p",String(arg[2]));
PrintTo(name,name,"\n");
if Length(arg)=2 then
   list:=BrauerChar(arg[1],arg[2]);
fi;
if Length(arg)=3 then
   list:=BrauerChar(arg[1],arg[2],arg[3]);
fi;
for entry in list do
    AppendTo(name," ",entry[1]," ");
    for i in [1..Length(entry[2])] do
        AppendTo(name,entry[2][i]," ");
    od;
    AppendTo(name,"\n");
od;
return;
end;

