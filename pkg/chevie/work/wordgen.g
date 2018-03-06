# WordEnumerator(G[,options])
#  returns an enumerator for the words in the generators of the group G.
#
#  Here is an example:
#
#  W:=ComplexReflectionGroup(4);
#  e:=WordEnumerator(W);
#  w:=Elements(Centre(W))[2];;
#  e.Get(w);
#  [-2,1,-2,1]
#  e.Get(w,"all");
#  [ [ -2, 1, -2, 1 ], [ -1, 2, -1, 2 ], [ 2, -1, 2, -1 ], [ 1, -2, 1, -2 ] ]
#
#  e.Get(w) returns a minimal word in the generators of G representing w in G.
#  In this word, the inverse of a generator is given as a negative index.
#
#  e.Get(w,"all") returns all minimal words for w.
#
#  e.UpTo(l[,"all"]) returns a decomposition (all decompositions if argument
#  "all" given) of all words of length <=l
#
#  e.First(n[,"all"]) returns a decomposition (all decompositions if argument
#  "all" given) of the first n words (in order of increasing length)
#
# options is a record of options allowing some variations in the performed task:
#  
#  .monoid  if this option is set, minimal positive words are returned.
#  .generators  if set, this is used instead of the generators of G
#  .genNames    if set, this is used as the alphabet for the returned words.
#
#  Examples with the same group and w as above:
#
#  e:=WordEnumerator(W,rec(monoid:=true));e.Get(w);
#  [2,1,1,2,1,1]
#  e:=WordEnumerator(W,rec(generators:=Reflections(W){[1..6]}));e.Get(w);
#  [ 5, 2, 1 ]
#
#  The data fields of a WordEnumerator are:
#
#  .elements     A list where the computed elements are accumulated
#  .cayleyGraph  A list of same length as .elements representing the Cayley
#                graph. The k-th element of .cayleyGraph is a list of pairs 
#                (i,j) such that .elements[k]=.elements[i]*.generators[j]
#  .nbLength     records how many words of each length have been found.
#  .group        W
#  .permGroup    true if W is a permutation group. For such groups there is
#                special code to represent compactly the elements as the image
#                of the base.
#  .generators   the generators used.
#  .genNames     the alphabet for the returned words.
#
WordEnumerator:=function(arg)local G,opt,e,g;
  G:=arg[1];
  if Length(arg)=2 then opt:=arg[2];else opt:=rec();fi;
  e:=rec(cayleyGraph:=[0],nbLength:=[1],permGroup:=IsPermGroup(G),group:=G);
  if e.permGroup then e.elements:=[Base(G)];else e.elements:=[G.identity];fi;
  if IsBound(opt.generators) then e.generators:=opt.generators;
    if IsBound(opt.genNames) then e.genNames:=opt.genNames;
    else e.genNames:=[1..Length(e.generators)];
    fi;
  else
    e.generators:=ShallowCopy(G.generators);
    e.genNames:=[1..Length(e.generators)];
    if not IsBound(opt.monoid) then
      for g in G.generators do
	if g<>g^-1 then Add(e.generators,g^-1);
	                Add(e.genNames,-Position(e.generators,g));fi;
      od;
    fi;
  fi;
  e.spin:=function(w,cond)local bag,new,nb; bag:=Set(e.elements);
    while not cond(w,bag) do
      if e.permGroup and Length(e.elements)=Size(e.group) then return false;fi;
      new:=List([1+Sum(e.nbLength{[1..Length(e.nbLength)-1]})..Length(e.elements)],
        h->List([1..Length(e.generators)],function(g)local n;
	 if e.permGroup then n:=OnTuples(e.elements[h],e.generators[g]);
	 else n:=e.elements[h]*e.generators[g];fi;
	 return [n,[h,e.genNames[g]]];end));
      new:=Concatenation(new);
      new:=CollectBy(new,x->x[1]);
      new:=List(new,x->[x[1][1],List(x,y->y[2])]);
      nb:=Difference(List(new,y->y[1]),bag);
      new:=Filtered(new,x->x[1] in nb);
      Append(e.cayleyGraph,List(new,x->x[2]));
      new:=List(new,x->x[1]);
      Append(e.elements,new);
      UniteSet(bag,new);
      if Length(new)>10 then
      InfoChevie("#I ",Length(new),
          " elements of length ",Length(e.nbLength),"\n");
      fi;
      Add(e.nbLength,Length(new));
    od;
    return true;
  end;
  e.decode:=function(i)
    if i=1 then return [[]];else return Concatenation(List(e.cayleyGraph[i],
    x->List(e.decode(x[1]),y->Concatenation(y,[x[2]]))));fi;
  end;
  e.Get:=function(arg)local w,i;w:=arg[1];
    if e.permGroup then w:=OnTuples(e.elements[1],w);fi;
    i:=Position(e.elements,w);
    if i=false and not e.spin(w,function(w,bag)return w in bag;end) 
    then Error("Element given not in ",e.group);return false;fi;
    i:=Position(e.elements,w);
    if Length(arg)>1 then return e.decode(i);fi;
    w:=[];
    while i<>1 do i:=e.cayleyGraph[i][1];Add(w,i[2]); i:=i[1];od;
    return Reversed(w);
  end;
  e.decodeAll:=function(arg)local n,i,decoded;
    n:=arg[1];
    if Length(arg)=1 then decoded:=[[]];
      for i in [2..n] do
	decoded[i]:=Concatenation(decoded[e.cayleyGraph[i][1][1]],
				    [e.cayleyGraph[i][1][2]]);
      od;
    else decoded:=[[[]]];
      for i in [2..n] do decoded[i]:=Concatenation(List(e.cayleyGraph[i],c->
	  List(decoded[c[1]],x->Concatenation(x,[c[2]]))));
      od;
    fi;
    return decoded;
  end;
  e.First:=function(arg)local n;n:=arg[1];
    e.spin(n,function(n,bag)return Sum(e.nbLength)>=n;end);
    return ApplyFunc(e.decodeAll,arg);
  end;
  e.UpTo:=function(arg)local l;l:=arg[1];
    e.spin(l,function(l,bag)return Length(e.nbLength)>l;end);
    arg[1]:=Sum(e.nbLength{[1..l+1]});
    return ApplyFunc(e.decodeAll,arg);
  end;
  e.operations:=rec(Print:=function(e)
    Print("<enumerated ",Length(e.elements)," words for ",e.group,">");end);
  return e;
end;
