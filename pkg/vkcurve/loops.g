#############################################################################
##
#A  loops.g       VKCURVE package         David Bessis
##
#Y  Copyright (C) 2001 - 2002  University Paris VII, France.
##
##  This file holds the implementation of LoopsAroundPunctures
## 
#############################################################################
#########################################################################
# Computes a list of generators of the fundamental group of
#       C - {y1,...,yn}
# Input: list of complex numbers, [y1,...,yn].
# Output: list of PL-loops encoded as a record with 3 fields
#         .points : list of all endpoints of all segments
#         .segments : list of all segments, each of them being
#                     encoded as the list of the position of its 2
#		      endpoints in .points
# 	  .loops  : list of the generators, each of them being encoded
#	             as a list of the positions in .segments of its 
#		     successive linear pieces 
########################################################################
# the debugging printing is controlled by VKCURVE.showLoops
#######################################################################
LoopsAroundPunctures:=function(originalroots) 
   local roots,res,y,ys,h,k,sy,z,t,newfriends,n,
         cut,newcirc,distneighbours,circleorigin,
	          rs,is,maxr,minr,maxi,mini,box,average,
	tan,cmplxscal,cnorm,cycorder,lineq,crossing,boundpaths,
	mediatrix,convert,neighbours,detectsleftcrossing,shrink;
		  
tan:= x -> x.i/x.r;
cmplxscal:=function(x,y) return x.r*y.r+x.i*y.i; end;
cnorm:= x->cmplxscal(x,x);

# Ordonne une liste de points trigonometriquement autour d'un centre
cycorder:=function(list,center) local y,right,left,top,bottom;
  right:=[]; left:=[]; top:=[]; bottom:=[];
  for y in list-center do
    if y.r > 0 then Add(right,y);
    elif y.r < 0 then Add(left,y);
    elif y.i > 0 then Add(top,y);
    else Add(bottom,y);
    fi;
  od;
  SortParallel(List(right,tan),right);
  SortParallel(List(left,tan),left);
  return Concatenation(right,top,left,bottom)+center;
end;

# Input: (list of complex numbers,complex number)
# Output: sublist of "neighbours" of the second input,
#         x and y are neighbours iff no z is in the disk of diameter
#                  [x,y]

neighbours:=function(list,center) local l,isneighbour;
  isneighbour:=function(y)local d,z;
    d:=cnorm(y-center);
    for z in l do
      if z<>y and cnorm(y-z)+cnorm(z-center)<= d then return false; fi;
    od;
    return true;
  end;
  l:=Filtered(list,y->y <> center);
  l:=Filtered(l,isneighbour);
  return cycorder(l,center);
end;

# value at z of an equation of the line (x,y)
lineq:=function(x,y,z)
  if x.r=y.r then  
    if x.i=y.i then Print("Undefined line\n"); return false;
    else return z.r-x.r; fi;
  else return (y.i-x.i)*(z.r-x.r)/(y.r-x.r)+x.i-z.i;
  fi;
end;

mediatrix:=function(x,y) local mid;
  if x = y then Print("Undefined mediatrix"); return false; fi;
  return (x+y)/2+[E(4),-E(4)]*(x-y);
end;

# Computes the intersecting point of two lines, each given by either a pair
# of points or a vector; returns false if the lines are parallel or a pair
# is a single
crossing:=function(arg) local x1,x2,y1,y2,lambdax,mux,lambday,muy,res,resr,resi;
if Length(arg)=4 then x1:=arg[1]; x2:=arg[2]; y1:=arg[3]; y2:=arg[4];
    else  x1:=arg[1][1]; x2:=arg[1][2]; y1:=arg[2][1]; y2:=arg[2][2];
fi;
if (x1=x2) or (y1=y2) then return false; fi; #Undefined line
if x1.r <> x2.r then
       lambdax:=(x1.i-x2.i)/(x1.r-x2.r);
       mux:=-lambdax*x1.r+x1.i;
       if y1.r <> y2.r then
           lambday:=(y1.i-y2.i)/(y1.r-y2.r);
           muy:=-lambday*y1.r+y1.i;
	   if lambdax=lambday then return false; fi;
	   resr:=(muy-mux)/(lambdax-lambday);
	   resi:=lambdax*resr+mux;
	   res:=Complex(resr,resi);
       else res:=crossing(E(3)*x1,E(3)*x2,E(3)*y1,E(3)*y2);
            if res = false then return false; fi;
            res:=res/E(3);
       fi;
else res:=crossing(E(4)*x1,E(4)*x2,E(4)*y1,E(4)*y2);
     if res = false then return false; fi;
     res:=res/E(4);
fi;
return res;
end;

detectsleftcrossing:=function(c,w,y,z) local res,med,a,b,x,k,xx;
res:=[1..Length(c)-1]; med:=mediatrix(y,z); a:=med[1]; b:=med[2];
for k in [1..Length(c)-1] do
    if lineq(a,b,c[k])*lineq(a,b,c[k+1])<= 0
       then
       x:=crossing(a,b,c[k],c[k+1]);
       if x = false then res[k]:=false;
          else xx:=(z-y)/(w[k]-y);
          if xx.i >= 0 then res[k]:=true;
                           else res[k]:=false;
          fi;
       fi;
       else res[k]:=false;
    fi;
od;
return res;
end;

# y must be an element of ys
boundpaths:=function(ys,sy,path,y) local z;
if not(IsBound(y.path)) then 
   y.path:=Concatenation(path,[y.y]);
   for z in y.lovers do boundpaths(ys,sy,y.path,sy(z)); od;
fi;
end;

# eliminates trivial segments and contracts pairs [a,b],[b,a]
shrink:=function(l)local k;
  k:=PositionProperty([1..Length(l)-1],i->l[i]=l[i+1]);
  if k<>false then
    return shrink(l{Cat([1..k],[k+2..Length(l)])});
  else
    k:=PositionProperty([1..Length(l)-2],i->l[i]=l[i+2]);
    if k<>false then return shrink(l{Cat([1..k],[k+3..Length(l)])});
    else return l;
    fi;
  fi;
end;

# converts old result format (list of loops) to the new format :
#  .points   : set of all endpoints of all segments
#  .segments : set of all segments used, where endpoints are indexed
#              as in points
#     .loops : list of sequence of numbers of used segments
convert:=function(oldres) local res,segments,loops,points;
res:=ShallowCopy(oldres);
segments:=Set(Concatenation(
        List(res,l->List([2..Length(l)],i->Set([l[i-1],l[i]])) )));
loops:=List(res,l->List([2..Length(l)],
	    function(i) local seg;
	      seg:=[l[i-1],l[i]];
	      if Position(segments,seg) = false then
	         return -Position(segments,Reversed(seg));
	      else
	         return Position(segments,seg);
	      fi;
	    end ));
points:=Set(Concatenation(segments));
segments:=List(segments,s->[Position(points,s[1]),Position(points,s[2])]);
return rec(points:=points,segments:=segments,loops:=loops);
end;

roots:=originalroots;
n:=Length(roots);
average:=Sum(roots)/n;
Sort(roots,function(x,y) local dx,dy;
    dx:=cnorm(x-average); dy:=cnorm(y-average);
    return dx < dy; end );
if n=1 then return convert([roots[1]+[Complex(1,0),Complex(0,1),
               Complex(-1,0),Complex(0,-1),Complex(1,0)]]);
fi;
ys:=List(roots,x->rec(y:=x));
sy:= y-> ys[Position(roots,y)];
for y in ys do
    y.neighbours:=neighbours(roots,y.y);
    y.friends:=[y.y];
    y.lovers:=[];
od;
if VKCURVE.showLoops then Print("neighbours computed\n");fi;
for y in ys do
    for z in y.neighbours do
        if not(z in y.friends) then
	   Add(y.lovers,z);
	   Add(sy(z).lovers,y.y);
	   newfriends:=Concatenation(y.friends,sy(z).friends);
           for t in y.friends do sy(t).friends:=newfriends; od;
	   for t in sy(z).friends do sy(t).friends:=newfriends; od;
	fi;
    od;
od;
for y in ys do
    distneighbours:=List(y.neighbours,z->(y.y.r-z.r)^2+(y.y.i-z.i)^2);
    SortParallel(distneighbours,y.neighbours);
od;
# To avoid trouble with points on the border of the convex hull,
# we make a box around all the points;
# The evalf trick is just in case we are dealing with cyclotomics
#rs:=List(ys,y->evalf(y.y.r)); is:=List(ys,y->evalf(y.y.i));
rs:=List(ys,y->y.y.r); is:=List(ys,y->y.y.i);
minr:=ys[Position(rs,Minimum(rs))].y.r;
maxr:=ys[Position(rs,Maximum(rs))].y.r;
mini:=ys[Position(is,Minimum(is))].y.i;
maxi:=ys[Position(is,Maximum(is))].y.i;
box:=[Complex(minr-2,mini-2),Complex(minr-2,maxi+2),
      Complex(maxr+2,mini-2),Complex(maxr+2,maxi+2),
      Complex((maxr+minr)/2,mini-(maxr-minr)/2-2),
      Complex((maxr+minr)/2,maxi+(maxr-minr)/2+2),
      Complex(minr-(maxi-mini)/2-2,(maxi+mini)/2),
      Complex(maxr+(maxi-mini)/2+2,(maxi+mini)/2)];
for y in ys do
    y.cycorder:=cycorder(
                  Concatenation(Filtered(roots,z-> (z<>y.y)) , box),
		  y.y);
    k:=Position(y.cycorder,y.neighbours[1]);
    y.cycorder:=Concatenation(
             Sublist(y.cycorder,[k..Length(y.cycorder)]),
	     Sublist(y.cycorder,[1..k-1]));
    Add(y.cycorder,y.cycorder[1]);
    y.circle:=[(y.y+y.neighbours[1])/2];
    y.witness:=[y.neighbours[1]];
    for z in Sublist(y.cycorder,[2..Length(y.cycorder)]) do
        cut:=detectsleftcrossing(y.circle,y.witness,y.y,z);
        if true in cut
        then
	    k:=Position(cut,true);
     	    y.circle:=Sublist(y.circle,[1..k]);
            y.witness:=Sublist(y.witness,[1..k]);
        fi;
        k:=Length(y.circle);
        newcirc:=crossing(mediatrix(y.y,y.witness[k]),mediatrix(y.y,z));
        if newcirc <> false then
	   Add(y.circle,newcirc);
           Add(y.witness,z);
	fi;
	if z in y.lovers then
           Add(y.circle,(y.y+z)/2); Add(y.witness,z);
	fi;
    od;
od;
if VKCURVE.showLoops then Print("circles computed\n"); fi;
boundpaths(ys,sy,[],ys[1]);
for y in ys do
    k:=Length(y.path);
    if k > 1 then circleorigin:=(y.y+y.path[k-1])/2;
      k:=Position(y.circle,circleorigin);
      y.circle:=Concatenation(
        Sublist(y.circle,[k..Length(y.circle)]),Sublist(y.circle,[1..k-1]));
    fi;
od;
for y in ys do
    k:=Length(y.path);
    y.handle:=Concatenation(List([1..k-1],i->
          Sublist( sy(y.path[i]).circle,
	[1..Position(sy(y.path[i]).circle,(y.path[i]+y.path[i+1])/2)])));
    y.loop:=Concatenation([y.handle,y.circle,Reversed(y.handle)]);
od;
for y in ys do y.loop:=shrink(y.loop); od;
Sort(ys,
       function(y1,y2)
       return Position(originalroots,y1.y)< Position(originalroots,y2.y);
       end);
return convert(List(ys,y->y.loop));
end;
