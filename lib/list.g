#############################################################################
##
#A  list.g                      GAP library                  Martin Schoenert
#A                                                           &  Werner Nickel
##
#A  @(#)$Id: list.g,v 1.5 1997/04/13 13:20:22 gap Exp $
##
#Y  Copyright 1990-1992,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains the library functions that  deal  mainly  with  lists.
##  It does not contain functions that deal with sets, vectors  or  matrices.
##
#H  $Log: list.g,v $
#H  Revision 1.5  1997/04/13 13:20:22  gap
#H  Again please makeinit
#H
#H  Revision 1.4  1997/04/13 12:58:25  gap
#H  Formatting to please makeinit
#H
#H  Revision 1.3  1997/04/13 12:49:05  gap
#H  Minimal formatting
#H
#H  Revision 1.2  1997/03/26 15:49:58  gap
#H  Applied Jean Michel's improvements to Sort, Sortex, SortParallel,
#H  Permuted, Reversed and BetaSet. Also added new functions PositionSet,
#H  PermListList, SortingPerm.
#H
#H      SL
#H
#H  Revision 1.1.1.1  1996/12/11 12:37:46  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 3.18.1.3  1995/12/05  14:52:11  mschoene
#H  added 'Compacted' and 'PositionBound'
#H
#H  Revision 3.18.1.2  1995/03/27  08:59:01  sam
#H  removed unnesseccary calls of 'Add' in 'List'
#H
#H  Revision 3.18.1.1  1994/10/24  16:25:16  htheisse
#H  improved 'Reversed' for ranges
#H
#H  Revision 3.18  1994/06/27  09:52:29  sam
#H  improved 'Maximum' and 'Minimum' for ranges
#H
#H  Revision 3.17  1994/06/12  13:47:57  mschoene
#H  changed 'RandomList', so that 'R_X' makes it into the init file
#H
#H  Revision 3.16  1993/02/04  18:38:32  martin
#H  changed 'List' to explode strings again (damn backward compatibility)
#H
#H  Revision 3.15  1992/01/02  15:27:51  jmnich
#H  added function Equivalenceclasses
#H
#H  Revision 3.14  1991/12/23  12:18:24  martin
#H  improved 'Maximum' and 'Minimum'
#H
#H  Revision 3.13  1991/12/19  13:02:20  martin
#H  renamed 'SupportPerm' to 'LargestMovedPointPerm'
#H
#H  Revision 3.12  1991/10/14  11:03:53  martin
#H  added 'Collected'
#H
#H  Revision 3.11  1991/10/09  12:36:32  martin
#H  added 'SortParallel', 'Sortex', and 'Permuted'
#H
#H  Revision 3.10  1991/10/09  08:33:24  martin
#H  added 'Flat'
#H
#H  Revision 3.9  1991/09/26  12:40:36  fceller
#H  changed 'PositionProperty' to return 'false'
#H
#H  Revision 3.8  1991/09/25  19:40:03  fceller
#H  added 'PositionProperty'
#H
#H  Revision 3.7  1991/09/25  13:35:23  fceller
#H  added 'Apply'
#H
#H  Revision 3.6  1991/06/28  15:00:00  martin
#H  moved 'ForAll', 'ForAny', and 'First' to the list package
#H
#H  Revision 3.5  1991/06/28  14:00:00  martin
#H  moved 'Maximum' and 'Minimum' to the list package
#H
#H  Revision 3.4  1991/06/28  13:00:00  martin
#H  added 'RandomList'
#H
#H  Revision 3.3  1991/06/28  12:00:00  martin
#H  removed 'GcdList' and 'LcmList'
#H
#H  Revision 3.2  1991/06/01  12:00:00  martin
#H  changed 'Quo' to 'QuoInt'
#H
#H  Revision 3.1  1990/11/12  12:00:00  martin
#H  added 'Number' and 'Cartesian'
#H
#H  Revision 3.0  1990/01/01  12:00:00  martin
#H  initial revision under RCS
#H
##


#############################################################################
##
#F  List( <obj> ) . . . . . . . . . . . . . . . . . . . . . convert to a list
##
List := function ( arg )
    local  lst,  i,  fun;

    if Length(arg) = 1  then
        if IsList(arg[1])  then
            lst := arg[1];
        elif IsPerm(arg[1])  then
            lst := [];
            for i  in [1..LargestMovedPointPerm(arg[1])]  do
                lst[i] := i ^ arg[1];
            od;
        elif IsWord(arg[1])  then
            lst := [];
            for i  in [1..LengthWord(arg[1])]  do
                lst[i] := Subword( arg[1], i, i );
            od;
        else
            Error("can't convert ",arg[1]," into a list");
        fi;

    elif Length(arg) = 2  and IsList(arg[1])  and IsFunc(arg[2])  then
        lst := [];  fun := arg[2];
        for i  in [ 1 .. Length( arg[1] ) ] do
            lst[i]:= fun( arg[1][i] );
        od;

    else
        Error("usage: List( <obj> ) or List( <list>, <func> )");

    fi;

    return lst;
end;


#############################################################################
##
#F  Apply( <list>, <func> ) . . . . . . . .  apply a function to list entries
##
##  'Apply' will  apply <func>  to  every  member  of <list> and replaces  an
##  entry by the corresponding return value.  Warning:  The previous contents
##  of <list> will be lost.
##
Apply := function ( list, func )
    local i;

    if not IsList( list ) or not IsFunc( func ) then
        Error( "usage: Apply( <list>, <func> )" );
    fi;

    for i in [1..Length( list )] do
        list[i] := func( list[i] );
    od;
end;


#############################################################################
##
#F  Concatenation( <list>, <list> ) . . . . . . . . . concatentation of lists
##
Concatenation := function ( arg )
    local  cat, lst;

    if Length(arg) = 1  and IsList(arg[1])  then
        cat := [];
        for lst  in arg[1]  do
            Append( cat, lst );
        od;

    else
        cat := [];
        for lst  in arg  do
            Append( cat, lst );
        od;

    fi;

    return cat;
end;

Cat:=Concatenation;

#############################################################################
##
#F  Flat( <list> )  . . . . . . . list of elements of a nested list structure
##
Flat := function ( lst )
    local   flt,        # list <lst> flattened, result
            elm;        # one element of <lst>

    # make the flattened list
    flt := [];
    for elm  in lst  do
        if not IsList(elm)  then
            Add( flt, elm );
        else
            Append( flt, Flat(elm) );
        fi;
    od;

    # and return it
    return flt;
end;

#############################################################################
##
#F  Reversed( <obj> )  . . . . . . . . . . .  reverse the elements in a list
##

Reversed := function( obj ) local  len, t;
    t:=TYPE(obj);
    if t="range" then
      len := Length( obj );
      return [ obj[ len ], obj[ len - 1 ] .. obj[ 1 ] ];
    elif t="record" then
      if IsBound(obj.operations) and IsBound(obj.operations.Reversed) then
        return  obj.operations.Reversed(obj);
      fi;
    else
       len := Length( obj );
       return obj{[len,len-1..1]};
    fi;
end;

#############################################################################
##
#F  Sublist( <list>, <list> ) . . . . . . . . . . .  extract a part of a list
##
Sublist:=function(list,ind)local i,sub;
  sub:=[];
  for i in [1..Length(ind)] do
    if IsBound(list[ind[i]]) then sub[i]:=list[ind[i]];fi;
  od;
  return sub;
end;

#############################################################################
##
#F  Filtered( <list>, <func> )  . . . . extract elements that have a property
##
Filtered := function ( list, func )
    local  flt,  elm;

    flt := [];
    for elm  in list  do
        if func( elm )  then
            Add( flt, elm );
        fi;
    od;

    return flt;
end;


#############################################################################
##
#F  Number( <list> [, <func>] ) . . . . . count elements that have a property
##
Number := function ( arg )
    local  nr,  elm;

    if Length(arg) = 1  then
        nr := 0;
        for elm  in arg[1]  do
            nr := nr + 1;
        od;

    elif Length(arg) = 2  then
        nr := 0;
        for elm  in arg[1]  do
            if arg[2]( elm )  then
                nr := nr + 1;
            fi;
        od;

    else
        Error("usage: Number( <list> ) or Number( <list>, <func> )");
    fi;

    return nr;
end;


#############################################################################
##
#F  Compacted( <list> ) . . . . . . . . . . . . . .  remove holes from a list
##
Compacted := function ( list )
    local    res,       # compacted of <list>, result
             elm;       # element of <list>
    res := [];
    for elm  in list  do
        Add( res, elm );
    od;
    return res;
end;


#############################################################################
##
#F  Collected( <list> ) . . . . . 
##
Collected := function ( list )
    local   col,        # collected of <list>, result
            elm,        # element of <list>
            nr,         # number of elements of <list> equal to <elm>
            i;          # loop variable

    col := [];
    for elm  in Set( list )  do
        nr := 0;
        for i  in list  do
            if i = elm  then
                nr := nr + 1;
            fi;
        od;
        Add( col, [ elm, nr ] );
    od;
    return col;
end;


#############################################################################
##
#F   CollectBy(<l>, <f>) . . . . . Collect by value of <f>
##
##   Let  v1,..,vn be  the distinct  values the  function <f>  takes on the
##   elements  of  <l>.  Returns  a  list  whose  i-th  item is the list of
##   elements  of <l> where <f> takes the value  vi. <f> may also be a list
##   representing List(l,f).

CollectBy:=function(l,vals)local res,i;
  if IsFunc(vals) then vals:=List(l,vals);else vals:=ShallowCopy(vals);fi;
  if Length(l)=0 then return [];fi;
  l:=ShallowCopy(l);SortParallel(vals,l);
  res:=[[l[1]]];
  for i in [2..Length(l)] do 
    if vals[i]=vals[i-1] then Add(res[Length(res)],l[i]);
    else Add(res,[l[i]]);
    fi;
  od;
  return res;
end;

#############################################################################
##
#F  Equivalenceclasses( <list>, <function> )  . calculate equivalence classes
##
##
##  returns
##
##      rec(
##          classes := <list>,
##          indices := <list>
##      )
##
Equivalenceclasses := function( list, isequal )
    local ecl, idx, len, new, i, j;

    if not IsList( list ) or not IsFunc( isequal ) then
        Error( "usage: Equivalenceclasses( <list>, <function> )" );
    fi;

    len := 0;
    ecl := [];
    idx := [];
    for i in [1..Length( list )] do
        new := true;
        j   := 1;
        while new and j <= len do
            if isequal( list[i], ecl[j][1] ) then
                Add( ecl[j], list[i] );
                Add( idx[j], i );
                new := false;
            fi;
            j := j + 1;
        od;
        if new then
            len := len + 1;
            ecl[len] := [ list[i] ];
            idx[len] := [ i ];
        fi;
    od;
    return rec( classes := ecl, indices := idx );
end;

#############################################################################
##
#F   Zip(<a1>,...,<an>,<f>) . . . Zips with function <f> the lists a1,..,an
##
##  Given the lists a1,..,an returns the list of f(a1[i],a2[i],..,an[i])
##
Zip:=function(arg)local f;f:=arg[Length(arg)];
  return List(TransposedMat(arg{[1..Length(arg)-1]}),v->ApplyFunc(f,v));
end;

#############################################################################
##
#F  ForAll( <list>, <func> )  . .  test a property for all elements of a list
##
ForAll := function ( list, func )
    local  l;

    for l  in list  do
        if not func( l )  then
            return false;
        fi;
    od;

    return true;
end;

#############################################################################
##
#F  ForAny( <list>, <func> )  . . . test a property for any element of a list
##
ForAny := function ( list, func )
    local  l;

    for l  in list  do
        if func( l )  then
            return true;
        fi;
    od;

    return false;
end;

#############################################################################
##
#F  First( <list>, <func> ) . .  find first element in a list with a property
##
First := function ( list, func )
    local  l;

    for l  in list  do
        if func( l )  then
            return l;
        fi;
    od;

    Error("at least one element of <list> must fulfill <func>");
end;

#############################################################################
##
#F  PositionProperty( <list>, <func> ) position of an element with a property
##
PositionProperty := function ( list, func )
    local   i;

    for i  in [ 1 .. Length( list ) ]  do
        if func( list[ i ] )  then
            return i;
        fi;
    od;
    return false;

end;

#############################################################################
##
#F  PositionBound( <list> ) . . . . . . . . . . position of first bound entry
##
PositionBound := function ( list )
    local   i;

    # look for the first bound element
    for i  in [1..Length(list)]  do
        if IsBound( list[i] )  then
            return i;
        fi;
    od;

    # no bound element found
    return false;
end;

#############################################################################
##
#F  Cartesian( <list>, <list>.. ) . . . . . . . .  cartesian product of lists
##
Cartesian2 := function ( list, n, tup, i )
    local  tups,  l;
    if i = n+1  then
        tup := ShallowCopy(tup);
        tups := [ tup ];
    else
        tups := [];
        for l  in list[i]  do
            tup[i] := l;
            Append( tups, Cartesian2( list, n, tup, i+1 ) );
        od;
    fi;
    return tups;
end;

Cartesian := function ( arg )
    if Length(arg) = 1  then
        return Cartesian2( arg[1], Length(arg[1]), [], 1 );
    else
        return Cartesian2( arg, Length(arg), [], 1 );
    fi;
end;

#############################################################################
##
#F  Sort( <list> )  . . . . . . . . . . . . . . . . . . . . . . . sort a list
##
##  Sort() uses Shell's diminishing increment sort, which extends bubblesort.
##  The bubble sort works by  running  through  the  list  again  and  again,
##  each time exchanging pairs of adjacent elements which are out  of  order.
##  Thus large elements "bubble" to the top, hence the name  of  the  method.
##  However elements need many moves to come close to their  final  position.
##  In shellsort the first passes do not compare element j with its  neighbor
##  but with the element j+h, where h is larger than one.  Thus elements that
##  aren't at their final position make large moves towards the  destination.
##  This increment h is diminished, until during the last  pass  it  is  one.
##  A good sequence of incremements is given by Knuth:  (3^k-1)/2,... 13,4,1.
##  For this sequence shellsort uses on average  approximatly  N^1.25  moves.
##
##  Shellsort is the method of choice to  sort  lists  for  various  reasons:
##  Shellsort is quite easy to get right, much easier than,  say,  quicksort.
##  It runs as fast as quicksort for lists with  less  than  ~5000  elements.
##  It handles both  almost  sorted  and  reverse  sorted  lists  very  good.
##  It works well  in  the  presence  of  duplicate  elements  in  the  list.
##  Says Sedgewick: "In short, if you have a sorting problem,  use the  above
##  program, then determine whether the extra effort required to  replace  it
##  with a sophisticated method will be worthwile."
##
##  Donald Knuth, The Art of Computer Programming, Vol.3, AddWes 1973, 84-95
##  Donald Shell, CACM 2, July 1959, 30-32
##  Robert Sedgewick, Algorithms 2nd ed., AddWes 1988, 107-123
##
##  In the case where theer is just one argument, we use the code supplied by Jean
##  Michel, using the internal function Set
##
Sort := function ( arg )
    local   list,  isLess,  i,  k,  h,  v, l, both;

    if Length(arg) = 1 and IsList(arg[1])  then
        list := arg[1];
        l:=Length(list);
        both := [  ];
        for i  in [ 1 .. l ]  do
            both[i] := [ list[i], i ];
        od;
        both:=Set( both );
        list{[1..l]} := both{[1..l]}[1];
    elif Length(arg) = 2  and IsList(arg[1])  and IsFunc(arg[2])  then
        list := arg[1];  isLess := arg[2];
        h := 1;  while 9 * h + 4 < Length(list)  do h := 3 * h + 1;  od;
        while 0 < h  do
            for i  in [ h+1 .. Length(list) ]  do
                v := list[i];  k := i;
                while h < k  and isLess( v, list[k-h] )  do
                    list[k] := list[k-h];   k := k - h;
                od;
                list[k] := v;
            od;
            h := QuoInt( h, 3 );
        od;

    else
        Error("usage: Sort( <list> ) or Sort( <list>, <func> )");
    fi;

end;

#############################################################################
##
#F  SortParallel(<list>,<list2> [,<func>])  . . .  sort two lists in parallel
##
SortParallel := function ( arg )
    local   lst,        # list <lst> to be sorted, first argument
            par,        # list <par> to be sorted parallel, second argument
            isLess,     # comparison function, optional third argument
            gap,        # gap width
            l, p,       # elements from <lst> and <par>
            i, k,       # loop variables
            both;       # special list to pass to Set for fast sorting

    if Length(arg) = 2  and IsList(arg[1])  then
        lst := arg[1];
        par := arg[2];
        l:=Length(lst);
        both := [  ];
        for i  in [ 1 .. l ]  do
            both[i] := [ lst[i], i , par[i]];
        od;
        both:=Set( both );
        for i  in [ 1 .. l ]  do
            lst[i] := both[i][1];
            par[i] := both[i][3];
        od;
        
    elif Length(arg) = 3  and IsList(arg[1])  and IsFunc(arg[3])  then
        lst := arg[1];
        par := arg[2];
        isLess := arg[3];
        gap := 1;  while 9*gap+4 < Length(lst)  do gap := 3*gap+1;  od;
        while 0 < gap  do
            for i  in [ gap+1 .. Length(lst) ]  do
                l := lst[i];  p := par[i];  k := i;
                while gap < k  and isLess( l, lst[k-gap] )  do
                    lst[k] := lst[k-gap];  par[k] := par[k-gap];  k := k-gap;
                od;
                lst[k] := l;  par[k] := p;
            od;
            gap := QuoInt( gap, 3 );
        od;

    else
        Error("usage: SortParallel(<lst>,<par>[,<func>])");
    fi;

end;

#############################################################################
##
#F   SortBy(<list>, <func>) . . . . . Sort <list> by value of <func>
##
SortBy:=function(list,func)
  SortParallel(List(list,func),list);
end;

#############################################################################
##
#F  Sortex(<list>) . . . sort a list (stable), return the applied permutation
##
Sortex := function ( list )
    local   both, perm, i;

    # make a new list that contains the elements of <list> and their indices
    both := [];
    for i  in [1..Length(list)]  do
        both[i] := [ list[i], i ];
    od;

    # sort the new list according to the first item (stable)
    both := Set( both );

    # copy back and remember the permutation
    perm := [];
    for i  in [1..Length(list)]  do
        list[i] := both[i][1];
        perm[i] := both[i][2];
    od;

    # return the permutation mapping old <list> onto the sorted list
    return PermList( perm )^(-1);
end;

############################################################################
##
#F  Permuted( <list>, <perm> ) . . .apply permutation <perm> to list <list>
##  
##  make maximum use of kernel functions
##
Permuted := function ( list, perm )
  return list{OnTuples([1..Length(list)],perm^-1)}; 
end;

#############################################################################
##
#F  PositionSorted( <list>, <elm> ) . . . .  find an element in a sorted list
##
##  'PositionSorted' uses a binary search instead of the linear  search  used
##  'Position'.  This takes log to base 2 of  'Length( <list> )' comparisons.
##  The list <list> must be  sorted  however  for  'PositionSorted'  to work.
##
##  Jon Bentley, Programming Pearls, AddWes 1986, 85-88
##
PositionSorted := function ( arg )
    local   list,  elm,  isLess,  l,  m,  h;

    if Length(arg) = 2  and IsList(arg[1])  then
        list := arg[1];  elm := arg[2];
        l := 0;  h := Length(list)+1;
        while l+1 < h  do               # list[l]<elm & elm<=list[h] & l+1<h
            m := QuoInt( l + h, 2 );    # l < m < h
            if list[m] < elm  then l := m;
            else                   h := m;
            fi;
        od;
        return h;                       # list[l]<elm & elm<=list[h] & l+1=h

    elif Length(arg) = 3  and IsList(arg[1])  and IsFunc(arg[3])  then
        list := arg[1];  elm := arg[2];  isLess := arg[3];
        l := 0;  h := Length(list)+1;
        while l+1 < h  do               # list[l]<elm & elm<=list[h] & l+1<h
            m := QuoInt( l + h, 2 );    # l < m < h
            if isLess( list[m], elm )  then l := m;
            else                            h := m;
            fi;
        od;
        return h;                       # list[l]<elm & elm<=list[h] & l+1=h

    else
        Error("usage: PositionSorted( <list>, <elm> [, <func>] )");
    fi;
end;

#############################################################################
##
#F  Product( <list> ) . . . . . . . . . . . product of the elements in a list
##
##  'Product( <list> )' \\
##  'Product( <list>, <func> )'
##
##  When used in the first way 'Product' returns the product of the  elements
##  of the list <list>.  When used in the second way  'Product'  applies  the
##  function <func>, which must  be  a  function  taking  one  argument,  and
##  returns the product of the results.  In either case if  <list>  is  empty
##  'Product' returns 1.
##
Product := function ( arg )
    local  list,  func,  prod,  i;

    if Length(arg) = 1  then
        list := arg[1];
        if Length(list) = 0  then
            prod := 1;
        else
            prod := list[1];
            for i  in [ 2 .. Length(list) ]  do
                prod := prod * list[i];
            od;
        fi;

    elif Length(arg) = 2  and IsList(arg[1])  and IsFunc(arg[2])  then
        list := arg[1];  func := arg[2];
        if Length(list) = 0  then
            prod := 1;
        else
            prod := func( list[1] );
            for i  in [ 2 .. Length(list) ]  do
                prod := prod * func( list[i] );
            od;
        fi;

    else
        Error("usage: Product( <list> ) or Product( <list>, <func> )");

    fi;

    return prod;
end;

#############################################################################
##
#F  Sum( <list> ) . . . . . . . . . . . . . . . sum of the elements of a list
##
Sum := function ( arg )
    local  list,  func,  sum,  i;

    if Length(arg) = 1  then
        list := arg[1];
        if Length(list) = 0  then
            sum := 0;
        else
            sum := list[1];
            for i  in [ 2 .. Length(list) ]  do
                sum := sum + list[i];
            od;
        fi;

    elif Length(arg) = 2  and IsList(arg[1])  and IsFunc(arg[2])  then
        list := arg[1];  func := arg[2];
        if Length(list) = 0  then
            sum := 0;
        else
            sum := func( list[1] );
            for i  in [ 2 .. Length(list) ]  do
                sum := sum + func( list[i] );
            od;
        fi;

    else
        Error("usage: Sum( <list> ) or Sum( <list>, <func> )");

    fi;

    return sum;
end;


#############################################################################
##
#F  Iterated( <list>, <func> )  . . . . . . .  iterate a function over a list
##
Iterated := function ( list, func )
    local  res,  i;
    if Length(list) = 0  then
        Error("Iterated: <list> must contain at least one element");
    fi;
    res := list[1];
    for i  in [ 2 .. Length(list) ]  do
        res := func( res, list[i] );
    od;
    return res;
end;

#############################################################################
##
#F  Maximum( <obj>, <obj>... )  . . . . . . . . . . . . . maximum of integers
##
Maximum := function ( arg )
    local   max, elm;
    if   Length(arg) = 1  and IsRange(arg[1])  then
        if Length(arg[1]) = 0  then
            Error("Maximum: <list> must contain at least one element");
        fi;
        max := arg[1][Length(arg[1])];
        if max < arg[1][1] then max:= arg[1][1]; fi;
    elif Length(arg) = 1  and IsList(arg[1])  then
        if Length(arg[1]) = 0  then
            Error("Maximum: <list> must contain at least one element");
        fi;
        max := arg[1][Length(arg[1])];
        for elm  in arg[1]  do
            if max < elm  then
                max := elm;
            fi;
        od;
    elif  Length(arg) = 2  then
        if arg[1] > arg[2]  then return arg[1];
        else                     return arg[2];
        fi;
    elif Length(arg) > 2  then
        max := arg[Length(arg)];
        for elm  in arg  do
            if max < elm  then
                max := elm;
            fi;
        od;
    else
        Error("usage: Maximum( <obj>, <obj>... ) or Maximum( <list> )");
    fi;
    return max;
end;

#############################################################################
##
#F  Minimum( <obj>, <obj>... )  . . . . . . . . . . . . . minimum of integers
##
Minimum := function ( arg )
    local   min, elm;
    if   Length(arg) = 1  and IsRange(arg[1])  then
        if Length(arg[1]) = 0  then
            Error("Minimum: <list> must contain at least one element");
        fi;
        min := arg[1][Length(arg[1])];
        if min > arg[1][1] then min:= arg[1][1]; fi;
    elif Length(arg) = 1  and IsList(arg[1])  then
        if Length(arg[1]) = 0  then
            Error("Minimum: <list> must contain at least one element");
        fi;
        min := arg[1][Length(arg[1])];
        for elm  in arg[1]  do
            if min > elm  then
                min := elm;
            fi;
        od;
    elif  Length(arg) = 2  then
        if arg[1] < arg[2]  then return arg[1];
        else                     return arg[2];
        fi;
    elif Length(arg) > 2  then
        min := arg[Length(arg)];
        for elm  in arg  do
            if min > elm  then
                min := elm;
            fi;
        od;
    else
        Error("usage: Minimum( <obj>, <obj>... ) or Minimum( <list> )");
    fi;
    return min;
end;

#############################################################################
##
#F  RandomList( <list> )  . . . . . . . . return a random element from a list
##
#N  31-May-91 martin 'RandomList' should be internal
##
R_N := 1;
R_X := [];

RandomList := function ( list )
    R_N := R_N mod 55 + 1;
    R_X[R_N] := (R_X[R_N] + R_X[(R_N+30) mod 55+1]) mod 2^28;
    return list[ QuoInt( R_X[R_N] * Length(list), 2^28 ) + 1 ];
end;

RandomSeed := function ( n )
    local  i;
    R_N := 1;  R_X := [ n ];
    for i  in [2..55]  do
        R_X[i] := (1664525 * R_X[i-1] + 1) mod 2^28;
    od;
    for i  in [1..99]  do
        R_N := R_N mod 55 + 1;
        R_X[R_N] := (R_X[R_N] + R_X[(R_N+30) mod 55+1]) mod 2^28;
    od;
end;

if R_X = []  then RandomSeed( 1 );  fi;


#############################################################################
##
#F  PositionSet( <l>, <x>[, <less> ) . . . . like 'Position', but the user
#F  is responsible for <l> beeing sorted
##  
##  This is only a slight variation of 'PositionSorted'. Due to Chevie team
##  The difference from PositionSorted is that it returns false if the object is
##  not present, rather than returning the position where it would be inserted
##  

PositionSet := function(arg)
    local  list, elm, isLess, l, m, h;
    if Length( arg ) = 2 and IsList( arg[1] )  then
        list := arg[1];
        elm := arg[2];
        l := 0;
        h := Length( list ) + 1;
        while l + 1 < h  do
            m := QuoInt( l + h, 2 );
            if list[m] < elm  then
                l := m;
            else
                h := m;
            fi;
        od;
    elif Length( arg ) = 3 and IsList( arg[1] ) and IsFunc( arg[3] )  then
        list := arg[1];
        elm := arg[2];
        isLess := arg[3];
        l := 0;
        h := Length( list ) + 1;
        while l + 1 < h  do
            m := QuoInt( l + h, 2 );
            if isLess( list[m], elm )  then
                l := m;
            else
                h := m;
            fi;
        od;
    else
        Error( "usage: PositionSet( <list>, <elm> [, <func>] )" );
    fi;
    if IsBound(list[h]) and list[h]=elm then
      return h;
    else
      return false;
    fi;
end;

#############################################################################
##
#F  Positions( <list> , <obj>) . . list of positions in <list> holding <obj>
##  
Positions:=function(l,o)local i,res;res:=[];
  for i in [1..Length(l)] do if l[i]=o then Add(res,i);fi;od;return res;
end;

#############################################################################
##
#F  PositionsProperty( <list> , <f>) . . list of positions i in <list> such
##    f(l[i]) holfd.
PositionsProperty:=function(l,f)local i,res;res:=[];
  for i in [1..Length(l)] do if f(l[i]) then Add(res,i);fi;od;return res;
end;

############################################################################
# PositionSublist(s,sub)
# returns the position of the first occurence of sub in s, or false
# if there is no such occurence
PositionSublist:=function(s,sub)local i,j;
  j:=1;
  for i in [1..Length(s)] do
    if s[i]=sub[j] then
      if j=Length(sub) then return i-Length(sub)+1;
      else j:=j+1;
      fi;
    else j:=1;
    fi;
  od;
  return false;
end;

#############################################################################
##
#F  SortingPerm( <list> ) . . . . . . . returns the same as 'Sortex( <list> )'
#F  but does *not* change the argument
##  
##    

SortingPerm := function ( list )
  local  both, perm, i, l;
  l:=Length(list);
  both:=[];
  for i in [1..l]  do
    both[i]:=[list[i],i];
  od;
  both:=Set(both);
  perm := [];
  perm{[1..l]} := both{[1..l]}[2];
  return PermList(perm)^-1;
end;

#############################################################################
##
#F  PermListList( <lst>, <lst2> ) . . . . what permutation of <lst> is <lst2>
##
##  PermListList finds which permutation p of [1..Length(lst)] is such 
##  that lst[i^p]=lst2[i] 
##  It returns false if there is no such permutation.
##
PermListList := function(l1,l2) local res;
  l1:=ShallowCopy(l1);l2:=ShallowCopy(l2); # to not destroy l1 and l2
  res:=Sortex(l2)* Sortex(l1)^-1;
  if l1<>l2 then return false; else return res;fi;
end;
