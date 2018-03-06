#############################################################################
##
#A  global.g       VKCURVE package         Jean Michel
##
#Y  Copyright (C) 2001 - 2002  University Paris VII, France.
##
##  This file holds the controlling functions which calls the various
##  steps of the Van Kampen algorithm.
## 
#############################################################################

# VKCURVE.Loops(r)
#  r should be a record with the fields
#      .roots   -- roots of the curve discriminant
#      .ismonic -- tells if the discriminant is monic in x
#  the function computes the following fields describing loops around the
#  .roots around from a basepoint:
#
#  .loops --  a list of loops, each described by a list of indices in r.segments
#             (a negative index tells to follow the reverse segment)
#  .segments -- oriented segments represented as a pair of indices in r.points.
#  .points -- a list of points stored a complex decimal numbers.
#  .basepoint -- holds the chosen basepoint
#
VKCURVE.Loops:=function(r)local m,i,P,l,segmentNumbers,uniquePoints;
  Inherit(r,LoopsAroundPunctures(r.roots));
  # here we  have loops around the  'true' roots and around  the 'extra'
  # roots Difference(r.roots,r.trueroots). We get rid of the extra loops
  # and the associated segments and points, first saving the basepoint.
  # (its location is known now, and maybe not later? J.M.) 
  if r.loops[1][1]<0 then r.basepoint:=r.segments[-r.loops[1][1]][2];
                     else r.basepoint:=r.segments[r.loops[1][1]][1];
  fi;
  if r.ismonic
      then r.loops:=r.loops{List(r.roots,z->Position(r.roots,z))};
      else r.loops:=r.loops{List(r.trueroots,z->Position(r.roots,z))};
  fi;
  segmentNumbers:=Union(List(r.loops,x->List(x,AbsInt)));
  r.loops:=List(r.loops,x->List(x,function(y)
    if y<0 then return -Position(segmentNumbers,-y);
           else return Position(segmentNumbers,y);
    fi;
    end));
  r.segments:=r.segments{segmentNumbers};
  uniquePoints:=Union(r.segments);
  r.segments:=List(r.segments,x->List(x,y->Position(uniquePoints,y)));
  r.points:=r.points{uniquePoints};
  r.basepoint:=Position(uniquePoints,r.basepoint);
  if VKCURVE.showSegments then
    Print("# There are ",Length(r.segments),
          " segments in ",Length(r.loops)," loops\n");
  fi;
  if VKCURVE.showWorst then
    l:=List([1..Length(r.segments)],function(i)local m,s;
     s:=r.segments[i];
     m:=MinPos(List(r.roots,z->DistSeg(z,r.points[s[1]],r.points[s[2]])));
     return [m[1],i,m[2]];
    end);
    Sort(l);
    Print("worst segments:\n");
    for i in [1..Minimum(5,Length(l))] do
      P:=l[i];
      Print("segment ",P[2],"=",r.segments[P[2]]," dist to ",P[3],
	     "-th root is ",P[1],"\n");
    od;
  fi;
# find the minimum distance m between two roots
  if Length(r.roots)>1 then
    m:=Dispersal(r.roots);
    if VKCURVE.showRoots then
      Print("\nMinimum distance=",evalf(m[1])," between roots ",m[2][1]," and ",
	 m[2][2]," of discriminant\n");
    fi;
    r.dispersal:=m[1];
  else
    r.dispersal:=1/1000;
  fi;
# and round points to m/100
  m:=-DecimalLog(Rational(r.dispersal)/100);
  r.points:=List(r.points,y->evalf(y,m));
end;

# VKCURVE.Zeros(r)
#   r should be a record with fields
#   .points
#   .curve
#
# It computes r.zeros as the zeros of the .curve at each of the .points
VKCURVE.Zeros:=function(r)local i,min,m;
  if VKCURVE.showRoots then
    Print("Computing zeros of curve at the ",Length(r.points),
          " segment extremities...\n"); 
  fi;
  r.zeros:=[];min:=[];
  for i in [1..Length(r.points)] do
    if VKCURVE.showZeros then Print("<",i,"/",Length(r.points),">");fi;
    r.zeros[i]:=SeparateRoots(Value(r.curve,
			 ["y",ComplexRational(r.points[i])]),1000);
    if Length(r.zeros[i])>1 then
      m:=Dispersal(r.zeros[i]);
      m[1]:=evalf(m[1],-DecimalLog(m[1]/1000));
      if VKCURVE.showZeros then Print(" d=",m,"\n");fi;
      min[i]:=[m[1],i];
    fi;
  od;
  if VKCURVE.showWorst and Length(r.zeros[1])>1 then
    Sort(min);
    Print("worst points:\n");
    for i in [1..Minimum(5,Length(min))] do
      Print(min[i][2],": mindist(zeros)=",min[i][1],"\n");
    od;
  fi;
end;

# VKCURVE.PrepareCurve(curve)
#   curve should be either a polynomial in y with coeffs polynomials in x
#   or an Mvp in x and y.
#   This  function  makes  sure  the  curve  is  quadratfrei  and makes its
#   coefficients  Rational or in  CF(4)=GaussianRationals (the coeffs could
#   be  decimals or complex  decimals) returns a  record with fields .curve
#   and .ismonic if the curve is monic in x
VKCURVE.PrepareCurve:=function(curve)local d;
  if IsPolynomial(curve) then 
    curve.coefficients:=List(curve.coefficients,z->Value(z,Mvp("x")));
    curve:=Value(curve,Mvp("y"));
  fi;
  if ForAny(curve.coeff,IsComplex) then
    curve.coeff:=List(curve.coeff,x->Cyclotomic(Complex(x)));
  fi;
  d:=MvpGcd(curve,Derivative(curve,"x"));
  if Length(Coefficients(d,"x"))>1 then 
    Print("**** Warning: curve is not quadratfrei: dividing by ",d,"\n");
    curve:=curve/d;
  fi;
  d:=Coefficients(curve,"x");d:=d[Length(d)];
  return rec(curve:=curve,ismonic:=Degree(d)=0);
end;

# VKCURVE.Discy(r)
#  r should be a record with field r.curve, a quadratfrei Mvp in x,y.
#  The discriminant of this curve with respect to x (a polynomial in y)
#  is computed.
#  First, the curve is split in
#    r.curveVerticalPart  -- the Gcd of the coeffs in x (an Mvp in y).
#    r.nonVerticalPart    -- curve/curveVerticalPart
#  Then disc=discriminant of r.nonVerticalPart is computed. Its quadratfrei
#  part  is computed, stripped  of factors common  with d and then factored
#  (if possible which in GAP3 means it is a polynomial over the rationals),
#  disc  is  stored  in  r.discy  as  a  GAP polynomial, and its factors in
#  r.discyFactored as a list of Mvp in x.
#  Some of  these computations may be  too costly for GAP,  in which case
#  one has a better hope to complete them in MAPLE.
VKCURVE.Discy:=function(r)local d,l,F,i,common;
  r.curveVerticalPart:=ApplyFunc(MvpGcd,Coefficients(r.curve,"x"));
  if VKCURVE.showRoots and Degree(r.curveVerticalPart)>0 then
    Print("Curve has ",Degree(r.curveVerticalPart)," linear factors in y\n"); 
  fi;
  r.nonVerticalPart:=r.curve/r.curveVerticalPart;
  d:=Discy(r.nonVerticalPart); # an Mvp in y
  if d=0 then # d=Discriminant of the curve minus vertical components
    Error("Discriminant is 0 but ",r.curve," should be quadratfrei"); 
  fi;
  if VKCURVE.showRoots then
    Print("Discriminant has ",Degree(d)," roots, "); 
  fi;
  d:=d/MvpGcd(d,Derivative(d,"y")); # make d quadratfrei
  if VKCURVE.showRoots then
    Print(" of which ",Degree(d)," are distinct\n"); 
  fi;
  common:=MvpGcd(d,r.curveVerticalPart);
  if VKCURVE.showRoots and Degree(common)>0 then
    Print(" and of which ",Degree(common)," are roots of linear factors\n"); 
  fi;
  d:=d/common; # make sure roots contain no vertical component
  d:=d/d.coeff[Length(d.coeff)];
  l:=ScalMvp(Coefficients(d,"y"));r.discy:=Polynomial(DefaultField(l),l);
  if r.discy.baseRing=GaussianRationals then
    d.coeff:=List(d.coeff,Complex);
    r.discyFactored:=[d];
  else r.discyFactored:=List(Factors(r.discy),x->Value(x,Mvp("x")));
  fi;
end;

VKCURVE.GetDiscyRoots:=function(r)local m;
  if VKCURVE.showRoots then
    Print("Computing roots of discriminant...\n"); 
  fi;
  if Degree(r.curveVerticalPart)=0 then r.roots:=[];
  else   r.roots:=SeparateRoots(r.curveVerticalPart,1000);
  fi;
  Append(r.roots,Concatenation(List(r.discyFactored,p->SeparateRoots(p,1000))));
end;

VKCURVE.Braids:=function(r)local i,l,bb,pr;
  if VKCURVE.showgetbraid then pr:=Print;else pr:=Ignore;fi;
  pr("# Computing monodromy braids\n");
  r.braids:=[];
  for i in [1..Length(r.loops)] do
    l:=Filtered(List(r.loops[i],AbsInt),s->not IsBound(r.monodromy[s]));
    if Length(l)>0 then 
      pr("# loop[",i,"] missing segments ",l,"\n");
    else
      bb:=Product(r.loops[i],function(s)if s<0 then return r.monodromy[-s]^-1;
	                                    else return r.monodromy[s];fi;end);
      pr("# loop[",i,"]=",bb,"\n");
      Add(r.braids,bb);
     fi;
  od;
  if VKCURVE.shrinkBraid then
    r.rawBraids:=r.braids;
    r.braids:=ShrinkGarsideGeneratingSet(r.braids);
  fi;
end;

VKCURVE.Segment:=function(r,segno)local pr,tm;
  if IsBound(r.name) then 
    PrintTo(SPrint(r.name,".",segno),"");
    pr:=function(arg)
      ApplyFunc(AppendTo,Concatenation([SPrint(r.name,".",segno)],arg));
      ApplyFunc(Print,arg);
    end;
  elif VKCURVE.showSegments then pr:=Print;
  else pr:=Ignore;
  fi;
  tm:=Runtime();
  if VKCURVE.monodromyApprox then
    r.monodromy[segno]:=ApproxFollowMonodromy(r,segno,pr);
  else
    r.monodromy[segno]:=FollowMonodromy(r,segno,pr);
  fi;
  tm:=Runtime()-tm;
  if IsBound(r.name) then pr(r.name,".");fi;
  pr("monodromy[",segno,"]:=");
  pr(AsWord(r.monodromy[segno]),";\n");
  pr("# segment ",segno,"/",Length(r.segments),
     " Time=",evalf(tm/1000,1),"sec\n");
end;

VKCURVE.SetPrintLevel:=function(printlevel)
  VKCURVE.showSingularProj:=printlevel>=2;
  VKCURVE.showBraiding:=printlevel>=2;
  VKCURVE.showLoops:=printlevel>=2;
  VKCURVE.showAction:=printlevel>=2;
  VKCURVE.showSegments:=printlevel>=1;
  VKCURVE.showInsideSegments:=printlevel>=2;
  VKCURVE.showWorst:=printlevel>=2;
  VKCURVE.showZeros:=printlevel>=2;
  VKCURVE.showNewton:=printlevel>=2;
  VKCURVE.showgetbraid:=printlevel>=1;
  VKCURVE.showRoots:=printlevel>=2;
end;

VKCURVE.TrivialCase:=function(r)
  r.presentation:=PresentationFpGroup(FreeGroup(
    Length(Coefficients(r.curve,"x"))-1));
  r.operations:=rec(Print:=function(r)DisplayPresentation(r.presentation);end);
  return r;
end;

# VKCURVE.Segments(name[,range])
# builds  files  describing  monodromy  along  segments  of r.segments (all
# segments by default, those in range if given).
# If  r  has  a  .name  field  then  computed segments are written on files
# beginning with r.name, otherwise are added to r.monodromy
# if r does not have the right fields then one tries to read (r.name).tmp
VKCURVE.Segments:=function(arg)local r,range,segno,tm,pr;
  r:=arg[1];
  Read(Concatenation(r.name,".tmp"));
  r.B:=Braid(CoxeterGroupSymmetricGroup(Length(r.zeros[1])));
  if not IsBound(r.monodromy) then r.monodromy:=[];fi;
  if Length(arg)>1 and IsList(arg[2]) then range:=arg[2];
  else range:=[1..Length(r.segments)];
  fi;
  for segno in range do VKCURVE.Segment(r,segno);od;
end;

VKCURVE.Finish:=function(r)
  VKCURVE.Braids(r);
  if r.ismonic then r.rawPresentation:=VKQuotient(r.braids);
               else r.rawPresentation:=DBVKQuotient(r);
  fi;
  r.presentation:=PresentationFpGroup(r.rawPresentation);
  ShrinkPresentation(r.presentation);
  r.rawPresentation:=PresentationFpGroup(r.rawPresentation);
  r.operations:=rec(Print:=function(r)DisplayPresentation(r.presentation);end);
  return r;
end;

VKCURVE.SearchHorizontal:=function(r0)local height,section,r;
  # Searching for a good horizontal
  height:=9;
  repeat height:=height+1;
    section:=Value(r0.curve,["x",height]);
    section:=section/MvpGcd(Derivative(section,"y"),section);
  until  Degree(section)=Length(Coefficients(r0.curve,"y"))-1
     and Degree(MvpGcd(Value(r0.discy,Mvp("y")),section))=0;
  section:=section/MvpGcd(section,r0.curveVerticalPart);
  Print("Curve is not monic in x -- Trivializing along horizontal line x = ",
          height,"\n");
  r:=rec(height:=height,input:=r0.curve,curve:=r0.curve*(x-height));
  if IsBound(r0.name) then r.name:=r0.name;fi;
  VKCURVE.Discy(r);
  # set trueroots  to roots  of Discy  which do  not only  correspond to
  # intersections of the curve with the chosen horizontal
  r.trueroots:=r0.roots;
  r.verticallines:=r0.roots{[1..Degree(r0.curveVerticalPart)]};
  r.roots:=ShallowCopy(r.trueroots);
  r.ismonic:=false;
  Append(r.roots,SeparateRoots(section,1000));
  return r;
end;

# curve --  an Mvp in x and y describing a curve in complex^2
FundamentalGroup:=function(arg)local c,r,i;
  c:=arg[1];
  if IsPolynomial(c) or IsMvp(c) then
    if Length(arg)=2 then VKCURVE.SetPrintLevel(arg[2]);
		     else VKCURVE.SetPrintLevel(0);fi;
    r:=VKCURVE.PrepareCurve(c);
    VKCURVE.Discy(r);VKCURVE.GetDiscyRoots(r);
    if Length(r.roots)=0 then return VKCURVE.TrivialCase(r);fi;
    if not r.ismonic then r:=VKCURVE.SearchHorizontal(r); fi;
    VKCURVE.Loops(r); VKCURVE.Zeros(r);
    if Length(r.zeros[1])=0 then return VKCURVE.TrivialCase(r);fi;
    r.B:=Braid(CoxeterGroupSymmetricGroup(Length(r.zeros[1])));
    r.monodromy:=[];
    for i in [1..Length(r.segments)] do VKCURVE.Segment(r,i);od;
    return VKCURVE.Finish(r);
  elif IsRec(c) and IsBound(c.operations) and 
    IsBound(c.operations.FundamentalGroup) then
    return ApplyFunc(c.operations.FundamentalGroup,arg);
  else
    Error(c," has no method for FundamentalGroup");
  fi;
end;

# Guarantees on LoopsAroundPunctures:
# For a set Z of zeroes and z in Z, let R(z):=1/2 dist(z,Z-z).
# The  input of  LoopsAroundPunctures is  a set  Z of approximate zeroes of
# r.discy such that for any z one of the zeroes is closer than R(z)/S where
# S is a global constant of the program (in practice we may take S=100).
# Let  d=inf_{z in  Z}(R(z)); we  return points  with denominator  10^-k or
# 10^-k<d/S'  (in practive we take S'=100) and  such that the distance of a
# segment to a zero of r.discy is guaranteed >= d-d/S'-d/S

# curve --  an Mvp in x and y "monic in x" describing a curve P(x,y)
PrepareFundamentalGroup:=function(curve,name)local r,i,f,fname;
  r:=VKCURVE.PrepareCurve(curve);r.name:=name;
  VKCURVE.Discy(r);VKCURVE.GetDiscyRoots(r);
  if Length(r.roots)=0 then return VKCURVE.TrivialCase(r);fi;
  if not r.ismonic then r:=VKCURVE.SearchHorizontal(r); fi;
  VKCURVE.Loops(r);VKCURVE.Zeros(r);
  if Length(r.zeros[1])=0 then return VKCURVE.TrivialCase(r);fi;
  fname:=Cat(r.name,".tmp");PrintTo(fname,"");
  for f in ["curve","discy","ismonic","roots","loops",
            "segments","points","zeros"] do
    Cut(SPrint(r.name,".",f,":=",FormatGAP(r.(f)),";\n"),rec(places:=",+*",
      file:=fname));
  od;
  Print("     ----------------------------------\n");
  Print("Data saved in ",fname,"\n");
  Print("You can now compute segments 1 to ",Length(r.segments),"\n");
  Print("in different GAP sessions by doing in each of them:\n");
  Print("    ",r.name,":=rec(name:=\"",r.name,"\");\n");
  Print("    VKCURVE.Segments(",r.name,",[1..",Length(r.segments),"]);\n");
  Print("(or some other range depending on the session)\n");
  Print("Then when all files ",r.name,".xx have been computed finish by\n");
  Print("    ",r.name,":=rec(name:=\"",r.name,"\");\n");
  Print("    FinishFundamentalGroup(",r.name,");\n");
end;

FinishFundamentalGroup:=function(r)local i;
  Read(Concatenation(r.name,".tmp"));
  if not IsBound(r.monodromy) then r.monodromy:=[];fi;
  r.B:=Braid(CoxeterGroupSymmetricGroup(Length(r.zeros[1])));
  for i in [1..Length(r.segments)] do
    Read(SPrint(r.name,".",i));
    if not IsBound(r.monodromy[i]) then Print("***** ",i," missing ****\n");fi;
    r.monodromy[i]:=ApplyFunc(r.B,r.monodromy[i]);
  od;
  return VKCURVE.Finish(r);
end;
