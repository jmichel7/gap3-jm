#M#################################################################
#M  induce.g: mostly handling functions                          ##
#M                                                               ##
#M     In a future release of SPECHT this file will contain all  ##
#M     of the functions for inducing and restricting and their   ##
#M     q-analogues, together with other funcitons needed for     ##
#M     calculations in the Fock space(s) etc.                    ##
#M                                                               ##
#M     Andrew Mathas                                             ##
#M                                                               ##
#M#################################################################

#D Functions for inducing and restricting modules.

#M Change log
#M 2.2: June 1996 
#M        o  moved these functions out of specht.g

#M General handling functions 

#M Usage: Specialized(x) or Specialized(x,a); defaults to a=1 (**undocumented)
Specialized:=function(arg) local x;
   x:=arg[1];
   if IsRec(x) and IsBound(x.operations) 
   and IsBound(x.operations.Specialized) then
     if Length(arg)=1 then return x.operations.Specialized(x,1);
     else return x.operations.Specialized(x, arg[2]);
     fi;
   else Error("Specialized(<x>), don't know how to specialize <x>.\n");
   fi;
end;

InducedModule:=function(arg)
  if arg=[] or arg[1] = false then  return false;
  elif IsBound(arg[1]) and IsRec(arg[1]) and IsBound(arg[1].operations)
  and IsBound(arg[1].operations.InducedModule) then
    return arg[1].operations.InducedModule(arg[1],arg{[2..Length(arg)]});
  else Error("InducedModule(<arg>), don't know how to induce <arg>");
  fi;
end;

SInducedModule:=function(arg)
  if arg=[] or arg[1] = false then  return false;
  elif IsBound(arg[1]) and IsRec(arg[1]) and IsBound(arg[1].operations)
  and IsBound(arg[1].operations.InducedModule) then
    return arg[1].operations.SInducedModule(arg[1],arg{[2..Length(arg)]});
  else Error("SInducedModule(), don't know how to S-induce ", arg);
  fi;
end;

RestrictedModule:=function(arg)
  if arg=[] or arg[1] = false then  return false;
  elif IsBound(arg[1]) and IsRec(arg[1]) and IsBound(arg[1].operations)
  and IsBound(arg[1].operations.RestrictedModule) then
    return arg[1].operations.RestrictedModule(arg[1],arg{[2..Length(arg)]});
  else Error("RestrictedModule(), don't know how to restrict ", arg);
  fi;
end;

SRestrictedModule:=function(arg)
  if arg=[] or arg[1] = false then  return false;
  elif IsBound(arg[1]) and IsRec(arg[1]) and IsBound(arg[1].operations)
  and IsBound(arg[1].operations.SRestrictedModule) then
    return arg[1].operations.SRestrictedModule(arg[1],arg{[2..Length(arg)]});
  else Error("SRestrictedModule(), don't know how to s-restrict ", arg);
  fi;
end;

InnerProduct:=function(a,b)
  if IsRec(b) and IsBound(b.operations) and IsBound(b.operations.InnerProduct) 
  then return b.operations.InnerProduct(a,b);
  else Error("InnerProduct(<a>,<b>), don't know how to take the inner product",
             " of <a> or <b>");
  fi;
end;

PositiveCoefficients:=function(x)
  if IsRec(x) and IsBound(x.operations) 
  and IsBound(x.operations.PositiveCoefficients) then
    return x.operations.PositiveCoefficients(x);
  else
    Error("PositiveCoefficients(<x>), don't know how to test <x> for ",
          "positive coefficients");
  fi;
end;

IntegralCoefficients:=function(x)
  if IsRec(x) and IsBound(x.operations) 
  and IsBound(x.operations.IntegralCoefficients) then
    return x.operations.IntegralCoefficients(x);
  else
    Error("IntegralCoefficients(<x>), don't know how to test <x> for ",
          "positive coefficients");
  fi;
end;

#M################################################################

#M Returns a string of the form "l1,l2,...,lk" for the list [l1,l2,..,lk]"
#M which contains no spaces, and has first element 1.
TightStringList:=function(list) local s, l;
  if list=[] then return ""; fi;
  s:=String(list[1]);
  for l in [2..Length(list)] do 
    PrintToString(s,",",list[l]); 
  od;
  return s;
end;

#M################################################################

#U RealLabelPartition(mu).................................................
#M Given a partition <mu> return a string for it in exponential notation. 
RealLabelPartition:=function(mu) local label, sep, tuple;
  if mu=[] then return "0"; fi;
  label:="";
  sep:="";
  for tuple in Reversed(Collected(mu)) do
    PrintToString(label,sep,String(tuple[1]));
    if tuple[2]>1 then PrintToString(label,"^",String(tuple[2])); fi;
    sep:=",";
  od;
  return label;
end;

#U TightLabelPartition(mu).................................................
#M Given a partition <mu> return a string for it in exponential notation. 
TightLabelPartition:=function(mu) local label, tuple;
  label:="";
  for tuple in Reversed(Collected(mu)) do
    PrintToString(label,String(tuple[1]));
    if tuple[2]>1 then PrintToString(label,"^",String(tuple[2])); fi;
  od;
  return label;
end;


#U LabelQuotient(quot)....................................................
#M Given a quotient return a string for printing it.                      
LabelQuotient:=function(quot) local label, sep, mu;
  label:="<";
  sep:="";
  for mu in quot do
    if mu=[] or mu=[0] then PrintToString(label,sep,".");
    else PrintToString(label,sep,LabelPartition(mu));
    fi;
    sep:="|";
  od;
  PrintToString(label,">");
  return label;
end;

#U TightTeXQuotient(q)...................................................
#M Given a quotient return a compact string for TeXing it. This string is
#M of the form <i_mu_1}j_mu_j}...> for the non-empty parts of q.         
TightTeXQuotient:=function(q) local str, sep, mu;
  str:="\\<";
  sep:="";
  for mu in [1..Length(q)] do
    if q[mu]<>[] then 
      PrintToString(str,sep,mu-1);
      sep:="|";
      if q[mu]<>[1] then PrintToString(str,"_{",TightLabelPartition(q[mu]),"}"); fi;
    fi;
  od;
  PrintToString(str,"\\>");
  return str;
end;

#U LabelPartition(mu)
#M Return a compact string for the partition <mu>.
LabelPartition:=function(mu) 
  return RealLabelPartition(mu); 
end;

#U LabelMultiPartition(mu).....................................
#M Returna a string for a multipartition.
LabelMultiPartition:=function(mu) local str, m;
  str:=LabelPartition(mu[1]);
  for m in mu{[2..Length(mu)]} do
    PrintToString(str,"|",LabelPartition(m));
  od;
  return str;
end;

#U SmallLabelPartition(mu)
#M Return the concatentation of the parts of mu which are assume 
#M "small" in the sense that they are less than 10.
SmallLabelPartition:=function(mu)
  local str, m;
  str:="";
  for m in mu do
    PrintToString(str,m);
  od;
  return str;
end;

#U LabelByEQuotient(e,mu)
#M Return a compact string for the e-quotient of the partition <mu>.
LabelByEQuotient:=function(e,mu) 
  local core, w, beads, qmu, s, i;
  core:=ECore(e,mu);
  w:=(Sum(mu)-Sum(core))/e;       # the e-weight of the block
  beads:=Length(core)+w*e;        # need at least this many beads
  #beads:=beads+((-beads-1) mod e);# we want most dominant in block to
                                  # have it e-quotient in the e^th component
  qmu:=List(EAbacusRunners(e,mu,beads), s->PartitionBetaSet(Set(s)));
  s:=""; 
  for i in qmu do
    PrintToString(s, String(SmallLabelPartition(i),-w), "|");
  od; 
  return s;
end;

#M Returns a string for PrintModule() from SpechtParts.labels, adding this 
#M string if it does not already exist.
StringPartition:=function(mu) local m, string, p;
  if mu=[] or mu=[0] then return "0";
  else return TightStringList(mu);
  fi;
end;

#M The function used by Specht to decide the format of partitions when
#M printing; see SpechtPrettyPrint).
SpechtPrintFn:=StringPartition;

#P Toggles te way in which Specht prints partition labels.
SpechtPrettyPrint:=function(arg)
  if Length(arg)=0 then
    if SpechtPrintFn=LabelPartition then
      SpechtPrintFn:=StringPartition;
    else SpechtPrintFn:=LabelPartition;
    fi;
  elif Length(arg)=1 and IsBool(arg[1]) then
    if arg[1] then SpechtPrintFn:=LabelPartition;
    else SpechtPrintFn:=StringPartition;
    fi;
  else Error("usage: SpechtPrettyPrint() or SpechtPrettyPrint(<bool>)\n");
  fi;
end;
