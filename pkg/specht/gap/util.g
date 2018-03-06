#D Generic utility functions

#U FactoredInteger(n)
#M Prints a nice factorization of the integer <n>.
FactoredInteger:=function(n) local f;
  for f in Collected(Factors(n)) do
    Print(f[1]);
    if f[2]>1 then Print("^",f[2]," "); else Print(" ");fi;
  od;
  Print("\n");
end;

#U pAdicExpansionInt(n,p)
#M Return the <p>-adic expansion of <n>.
pAdicExpansionInt:=function(n,p)
  local expansion, pk, r;
  if not IsInt(n) or n<1 then 
    Error("pAdicExpansionInt: <n> must be a positive integer"); 
  fi;
  expansion:=[];
  while n>=0 do
    r:= n mod p;
    Add(expansion, r);
    n:=(n-r)/p;
  od;
  return expansion;
end;

#U remove the <r>th entry from the list <list>. 
DeleteListEntry:=function(list, r)
  if IsBound(list[r]) then
    list{[r..Length(list)-1]}:=list{[r+1..Length(list)]};
    Unbind(list[Length(list)]);
  fi;
end;

#U SwapListEntries(list, i,j )                                           
#M Swap the <i>th and <j>th entries of the list <list>.                  
SwapListEntries:=function(list, i,j)
  local mi;
  mi:=list[i];
  list[i]:=list[j];
  list[j]:=mi;
end;

#U IsProperSubset(A,B)
#M Return true if B is a proper subset of A
IsProperSubset:=function(A,B)
  return A<>B and IsSubset(A,B);
end;

#U ReplaceAllStrings(str,old,new)
#M Replace ALL occurrences of <old> in <str> with <new>
ReplaceAllStrings:=function ( str, old, new )
  local i, res;
  i:=1;
  res := "";
  while i <= Length(str)-Length(old)+1 do
    if str{[i ..i+Length(old)-1]} = old  then
      Append(res, new);
      i:=i+Length(old);
    else
      Add(res, str[i]);
      i:=i+1;
    fi;
  od;
  if i>Length(str)-Length(old)+1 and i<=Length(str) then
    Append(res, str{[i..Length(str)]});
  fi;
  return res;
end;

## return a string which will reconstruct the record <rec>
SaveRecord:=function(savor)
  local str, sep, s;
  if IsList(savor) then
    str:=""; sep:="[";
    for s in savor do
      str:=AddToString(str, sep, SaveRecord(s));
      sep:=",";
    od;
    return AddToString(str, "]");
  else
    if IsRec(savor) and IsBound(savor.operations) and IsBound(savor.operations.Save) 
    then return savor.operations.Save(savor);
    else Error("no save operation for <savor>");
    fi;
  fi;
end;

#U GoodString
#M A hack to get around a silly bug in String() for polynomials
GoodString:=function(str)
  str:=String(str);
  str:=ReplaceAllStrings(str, "(","");
  str:=ReplaceAllStrings(str, ")","");
  str:=ReplaceAllStrings(str, "+-","-");
  return str;
end;

#U PPower(p,n)                                                           
#M Find the largest power of <p> dividing <n>.                           
PPower:=function(p, n) local r;
  r:=0;
  while n mod p =0 do
    n:=n/p;
    r:=r+1;
  od;
  return r;
end;
