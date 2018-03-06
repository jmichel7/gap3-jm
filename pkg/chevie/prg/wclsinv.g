#############################################################################
##
#A  wclassinv.g               CHEVIE library                     Frank Luebeck
##
#Y  Copyright (C) 1992 - 1996  Lehrstuhl D fuer Mathematik, RWTH Aachen, IWR
#Y  der Universitat Heidelberg, University of St. Andrews, and   University 
#Y  Paris VII.
##  
##  This  file contains   special  functions for  computing invariants  of
##  conjugacy classes in Coxeter groups and Coxeter cosets.
##  
##  The  main  difference  to  the  general functions 'ClassInvariants' and
##  'PositionClass'  for permutation groups is  the use of the 'ClassParam'
##  functions defined for the different types of irreducible groups. If the
##  class  of an element  cannot be identified  cheaply then the element is
##  transformed  into a  word and  the CHEVIE  label for  its corresponding
##  conjugacy class is computed by quite efficient functions.
##  
#############################################################################
##
#F  ComponentWordsPerm( <W>, <w> ) . . . . . . .  similar to 'CoxeterWord'
#F  but giving lists of  words in standard generators of simple components
##    
CoxeterGroupOps.ComponentWords:=function(W,w)
  w:=W.rootRestriction{CoxeterWord(W,w)};
  return List(W.type,t->List(Filtered(w,y->y in
    t.indices),i->Position(t.indices,i)));
end;

#F  ComponentWordsPerm( <WF>, <w> ). . . . . similar to 'CoxeterWord'
#F  but giving lists of  words in standard generators of simple components
CoxeterCosetOps.ComponentWords:=function(W,w)
  w:=W.reflectionGroup.rootRestriction{CoxeterWord(W,w)};
  return List(W.type,function(o)
    w:=List(o.orbit,t->List(Filtered(w,y->y in t.indices),
      i->Position(t.indices,i)));
    # change to other representative in same F-conjugacy class,
    # which is of form (x,1,...,1) on F-orbits of simple components:
    return Concatenation(w[1],
      Concatenation(List([Length(w),Length(w)-1..2],j->w[j])));
    end);
end;

#############################################################################
##
#F  ClassParamCheckFunction( <r>, <l> ) . . . . . . . . . returns function
#F  which computes a label for the (F-)conjugacy class of a given element 
##  
##  This function can be used as argument for 'ClassInvariants' for Coxeter
##  groups or Coxeter cosets, respectively. The returned function is always
##  able  to distinguish  all classes  (and is  quite effective compared to
##  backtrack conjugacy tests).
##  
##  It uses the classification.
##  
CoxeterGroupOps.ClassParamCheckFunction:=function(r,l)
  return x->Zip(r.g.type,CoxeterGroupOps.ComponentWords(r.g,x),
     function(t,w)return CHEVIE.Data("ClassParameter",t,w);end);
end;
       
CoxeterCosetOps.ClassParamCheckFunction:=function(r,l)
  return x->Zip(r.g.coset.type,CoxeterCosetOps.ComponentWords(r.g.coset,x),
      function(t,w)return CHEVIE.Data("ClassParameter",t,w);end);
end;

#############################################################################
##
#F  CoxeterGroupOps.ClassInvariants(  ... )  . . . . . 'ClassInvariants' for
#F  Coxeter groups
##
##  Here we use ClassParameter for irreducible types to  distinguish classes
##  which are not distinguished by cycle  type (and cycle type of elements
##  multiplied by center elements and not lying in small classes).
##  
CoxeterGroupOps.ClassInvariants:=function(arg)
  if IsList(arg[1]) then arg:=arg[1]; fi;
  if Length(arg)=1 then
    return PermGroupOps.ClassInvariants(arg[1],
                   CentreMultFunction,
                   ShortClassListFunction,
                   CoxeterGroupOps.ClassParamCheckFunction);
  else
    return PermGroupOps.ClassInvariants(arg);
  fi;
end;

#############################################################################
##
#F  CoxeterCosetOps.ClassInvariants(  ... )  . . . . . 'ClassInvariants' for
#F  Coxeter cosets
##
##  Here we use the  ClassParameter function to  find the class of elements
##  which are not lying in small classes.
##  
CoxeterCosetOps.ClassInvariants:=function(WF)local tmp;
  if IsList(WF) then WF:=WF[1];fi;
  
  # in non twisted case use the function for Coxeter groups:
  if WF.phi=() then return ClassInvariants(Group(WF));fi;
  
  tmp:=Copy(Group(WF).generators);Add(tmp,WF.phi);
  tmp:=Group(tmp,());
  tmp.conjugacyClasses:=ConjugacyClasses(WF);
  tmp.size:=Size(WF);
  tmp.coset:=WF;
  
  return PermGroupOps.ClassInvariants(tmp,
    ShortClassListFunction,
    CoxeterCosetOps.ClassParamCheckFunction);
end;
