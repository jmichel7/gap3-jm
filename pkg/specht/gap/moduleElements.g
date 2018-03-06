#######################################################################
##  moduleElements.g - defines a "module" record structure           ##
##                                                                   ##
##     Andrew Mathas          mathas@maths.usyd.edu.au               ##
##     University of Sydney   Sydney, 2006                           ##
##                                                                   ##
#######################################################################

#D Defines generic records for working with modules and with algebras

# Next we define the various operations for the elements of a module,
# placing them in the generic ModuleOps operations record. We follow
# Jean Michel's mvp code.
ModuleOperations:=OperationsRecord("ModuleOps");

#U MakeModuleElement(module,elts,coeffs)                               
#M This is function is used to create all elements. It is not meant to 
#M be used "externally", so there is no error checking.                
#M An element of the space is a record with the following components:  
#M   module     = <module> = module record                             
#M   elts       = a set of elements indexing the basis                 
#M   coeffs     = a parallell list of coefficents to elts              
#M   operations = <module>.Operations                                  
#C To maximize speed, this function does not check for syntax errors   
#C and it does not normalize the module element that is returned.      
ModuleOperations.MakeElement:=function(parent,elts,coeffs) 
  return rec(parent:=parent, elts:=elts, coeffs:=coeffs,
             operations:=parent.Operations);
end;

#U MakeModuleElement(module,elt)..........................................
#U MakeModuleElement(module,elts,coeffs)..................................
#M This is function is the front end to the function                      
#M ModuleOperations.MakeElement(). It is a little more flexible in        
#M its usage and it does some error checking and it returns a normalized  
#M basis element.                                                         
MakeModuleElement:=function(arg) local parent, elts, coeffs;
  parent:=arg[1];
  if Length(arg)=3 then     # parallel lists of elts and coefficients
    elts:=arg[2];
    coeffs:=arg[3];
  elif Length(arg)=2 then   # a single index
    elts:=[arg[2]];
    coeffs:=[parent.baseRing.one];
  else Error("MakeModuleElement(<parent>, elts [, coeffs])");
  fi;
  return parent.Operations.Normalize(
           parent.Operations.MakeElement(parent,elts,coeffs));
end;

#U IsModuleElement(x).....................................................
#M A crude test to see if a record is an element of some module.          
IsModuleElement:=function(x)
  return IsRec(x) and IsBound(x.parent);
end;

#U Normalize(b)                                                           
#M Collect together all common elts in a module record.                   
ModuleOperations.Normalize:=function(b)
local elts, coeffs, i;
  if b.elts=[] then return b; fi;
  SortParallel(b.elts,b.coeffs);
  elts:=[b.elts[1]];coeffs:=[b.coeffs[1]];
  for i in [2..Length(b.elts)] do
    if b.elts[i]<>b.elts[i-1] then
      Add(elts,b.elts[i]);
      Add(coeffs,b.coeffs[i]);
    else 
      coeffs[Length(coeffs)]:=coeffs[Length(coeffs)]+b.coeffs[i];
    fi;
  od;
  i:=List(coeffs,x->x<>x*0);
  b.coeffs:=ListBlist(coeffs,i); b.elts:=ListBlist(elts,i);
  return b;
end;

# Zero element of module
ModuleOperations.Zero:=a->a.operations.MakeElement(a.parent,[],[]);

# Addition of module elements
ModuleOperations.\+:=function(a,b) local B;
  if a=0 then return b;
  elif b=0 then return a;           # to get Sum() to work...
  elif not IsIdentical(a.operations,b.operations)  then
    Error("<a> and <b> belong to different algebras");
  fi;
  return
  a.operations.Normalize(a.operations.MakeElement(a.parent,
      Concatenation(a.elts,b.elts),Concatenation(a.coeffs,b.coeffs)) );
end;

# Subtraction of module elements
ModuleOperations.\-:=function(a,b) local B;
  if a=0 then return -1*b;
  elif b=0 then return a;           # to get Sum() to work...
  elif not IsIdentical(a.operations,b.operations) then
    Error("<a> and <b> belong to different algebras");
  fi;
  return
  a.operations.Normalize(a.operations.MakeElement(a.parent,
      Concatenation(a.elts,b.elts),Concatenation(a.coeffs,-b.coeffs)) );
end;

# Scalar multiplication of module elements
ModuleOperations.\*:=function(a,b) local c;
  if IsList(b) then return List(b, c->a*c);
  elif IsModuleElement(a) and b in a.parent.baseRing then
    return a.operations.MakeElement(a.parent, a.elts, a.coeffs*b);
  elif IsModuleElement(b) and a in b.parent.baseRing then
    return b.operations.MakeElement(b.parent, b.elts, b.coeffs*a);
  else Error("I do not know how to compute the product a*b");
  fi;
end;

#M This function returns a string for the module element, where the way  
#M that the indexes are printed is controled by the function indexString.
#M The intention is that this function can be called with different      
#M indexString() functions to create functions like String(), Print(),   
#M Display() and TeX() for the module elements.                          
ModuleOperations.GeneralString:=function(a,indexString,coeffString)
  local non0, elts, coeffs, res, coeff, i;
  if IsBound(a.name) then return a.name; fi;
  # we now use a hack so that we can print <a> as a linear
  # combination of basis elements if it has .basis and .basiscoeffs
  # records; if not, we just print in the usual way. Each element of
  # <basis> should have a <name> record
  if IsBound(a.basis) and IsBound(a.basiscoeffs) then 
    non0:=Filtered([1..Length(a.basis)], i->a.basiscoeffs[i]<>0*a.basiscoeffs[i]);
    elts:=a.basis{non0};
    coeffs:=a.basiscoeffs{non0};
    indexString:=x->x.name;
  else
    elts:=a.elts;
    coeffs:=a.coeffs;
  fi;
  if coeffs=[] and elts=[] then return coeffString(0); fi;
  res:="";
  for i in [1..Length(elts)] do
    coeff:=coeffString(coeffs[i]);
    if coeff="1" then coeff:="";
    elif coeff="-1" then coeff:="-";
    elif not IsInt(coeffs[i]) and ('+' in coeff or '-' in coeff) then 
	coeff:=Concatenation("(",coeff,")");
    fi;
    if Position(coeff,'-')<>1 and i<>1 then Append(res,"+");fi;
    PrintToString(res,String(coeff));
    if a.parent.PrintOption="GAP" and coeff<>"" and coeff<>"-" then 
      Append(res,"*");
    fi;
    PrintToString(res,indexString(elts[i]));
  od;
  return String(res);
end;
# Display does "normal" printing
ModuleOperations.Display:=function(a,dummy) local b;
  b:=ShallowCopy(a);
  if IsBound(b.name) then Unbind(b.name); fi;
  if IsBound(b.basis) then Unbind(b.basis); fi;
  Print(b.operations.GeneralString(b,b.operations.indexString,String));
end;
ModuleOperations.Print:=function(a)
  Print(a.operations.GeneralString(a,a.operations.indexString,String));
end;
ModuleOperations.String:=function(a)
  return a.operations.GeneralString(a,a.operations.indexString,String);
end;

##########################################################################
  
#U ModuleRecord(name,indexString,baseRing[, attributes])..................
#M This function returns a single "module" record which can be used to    
#M work within the a module. The elements in this module are created via  
#M function calls of the form MakeModuleElement(module, *). The elements  
#M of this "module" are records with the following components:            
#M   elts = list of indexing elements for the additive structure          
#M   coeffs = parallel list of coefficents                                
#M   operations = <name>Ops = operations record which control printing    
#M            (via the string function above), addition etc.              
#M The function <indexstring> is used by operations.GeneralString to make 
#M functions Print(), String(), TeX() etc for the module.                 
ModuleRecord:=function(arg) local name, indexString, baseRing, module, att;
  if not (Length(arg) in [3,4] and IsString(arg[1]) and IsFunc(arg[2]) 
          and IsRing(arg[3]) ) then
    Error("usage, ModuleRecord(name,indexString,baseRing[, attributes])");
  fi;
  name:=arg[1];
  indexString:=arg[2];
  baseRing:=arg[3];
  module:=rec(name:=name, 
	            baseRing:=baseRing,
              PrintOption:="Pretty",
	            operations:=OperationsRecord((SPrint(name,"Ops")))
  );
  module.operations.Print:=M->Print("ModuleRecord(",M.name,")");
  # membership test
  module.operations.\in:=function(x, M)
    return IsRec(x) and IsBound(x.parent) and x.parent=M;
  end;
  if Length(arg)=4 and IsRec(arg[4]) then 
    for att in RecFields(arg[4]) do
      module.(att):=arg[4].(att);
    od;
  fi;
  # finally define the module operations record which gives operations
  # for the elements of the module. We inherit operations from the
  # generic module operations record ModuleOperations
  module.Operations:=OperationsRecord((SPrint(name,"Ops")),ModuleOperations);
  module.Operations.indexString:=indexString;
  return module;
end;

##########################################################################

# Operations for algebras inherit the module operations
AlgebraOperations:=OperationsRecord("AlgebraOps",ModuleOperations);

#C We should really define different functions but...                    

#F IsAlgebraElement
# IsAlgebraElement:=IsModuleElement; # comment JM 2/2015: conflict with gap3
#F MakeAlgebraElement
MakeAlgebraElement:=MakeModuleElement;
AlgebraOperations.MakeElement:=AlgebraOperations.MakeElement;

# Multiplication of algebra elements - we deal with scalars and lists      
# generically and pass other multiplications to operations.multiplication()
AlgebraOperations.\*:=function(a,b) local elts, coeffs, c, i, j;
  if IsList(b) then return List(b, c->a*c);
  elif IsModuleElement(a) and b in a.parent.baseRing then
    return a.operations.MakeElement(a.parent, a.elts, a.coeffs*b);
  elif IsModuleElement(b) and a in b.parent.baseRing then
    return b.operations.MakeElement(b.parent, b.elts, b.coeffs*a);
  else
    elts:=[];coeffs:=[];
    for i in [1..Length(a.elts)] do
      for j in [1..Length(b.elts)] do
	c:=a.operations.multiplication(a.parent,a.elts[i],b.elts[j]);
	Append(elts,c.elts);
	Append(coeffs,a.coeffs[i]*b.coeffs[j]*c.coeffs);
      od;
    od;
    return a.operations.Normalize(a.operations.MakeElement(a.parent,elts,coeffs) );
  fi;
end;

# exponentiation
AlgebraOperations.\^:=function(a,k) local ak, b, r;
  if k=0 then return a.parent.identity;
  elif k>0 then
    # use repeated squaring to evaluate a^k
    ak:=a.parent.identity;       # this will become a^k
    b:=ShallowCopy(a);               # we keep squaring b
    while k>0 do
      r:=k mod 2;
      k:=(k-r)/2;
      if r=1 then ak:=ak*b; fi;
      b:=b*b;
    od;
    return ak;
  elif k<0 then
    if IsBound(a.operations.inverse) then
      b:=a.operations.inverse(a);
#     if IsAlgebraElement(b) then return b^k;
      if IsModuleElement(b) then return b^k;
      else Error("a^k: do not know how to compute negative powers of <a>");
      fi;
    fi;
  fi;
end;

#U AlgebraRecord(name,indexString,baseRing,identity,ops)..................
#M This is a general function for creating records for algebras. it looks 
#M and behaves very much like the function ModuleRecord() above, the only 
#M difference being that elements of an algebra can be multiplied. All of 
#M the arguments to this function have the same meaning as above and, in  
#M addition, <identity> is the identity element of the algebra and <ops>  
#M is a record which contains additional functions which are to be added  
#M to the operations records of of the algebras elements. One of these    
#M functions must be <ops>.multiplication=function(A,a,b) which specifies 
#M how to compute the product a*b for two elements <a> and <b> of the     
#M algebra <A>. Note that the function <ops>.multiplication should only   
#M multiply two elements of the algebra as multiplications by scalars,    
#M lists etc are handled generically. Other functions which can be given  
#M include an inverse() funciton,MakeElements() etc.                      
AlgebraRecord:=function(name,indexString,baseRing,identity,ops) 
  local algebra, o;
  if not ( IsString(name) and IsFunc(indexString) 
     and (IsRing(baseRing) or IsField(baseRing))
     and IsRec(ops) and IsBound(ops.multiplication) )
  then 
    Error("usage: AlgebraRecord(name,indexString,baseRing,identity,ops)"); 
  fi;
  algebra:=rec(name:=name,
	       baseRing:=baseRing,
               isDomain:=true,
               PrintOption:="Pretty",
	       operations:=OperationsRecord((SPrint(name,"Ops")))
  );
  algebra.operations.Print:=M->Print("AlgebraRecord(",M.name,")");
  # membership test
  algebra.operations.\in:=function(x, M)
    return IsRec(x) and IsBound(x.parent) and x.parent=M;
  end;
  # finally define the module operations record which gives operations
  # for the elements of the module. We inherit operations from the
  # generic module operations record ModuleOperations
  algebra.Operations:=OperationsRecord((SPrint(algebra.name,"Ops")),
                                       AlgebraOperations);
  algebra.identity:=MakeAlgebraElement(algebra,[identity],[One(algebra.baseRing)]);
  algebra.Operations.indexString:=indexString;
  algebra.operations.One:=a->a.identity;
  algebra.operations.Zero:=a->algebra.Operations.MakeElement(a,[],[]);
  for o in RecFields(ops) do
      algebra.Operations.(o):=ops.(o);
  od;
  return algebra;
end;

##########################################################################

#U TriangularBasis(basis)                                                 
#M Given a basis <basis> return, if possible, a reordering of <basis>,    
#M together with a set of distinguished 'words', that lets us use a       
#M triangular algorithm to rewrite an arbitrary element in the span of    
#M <basis> as a linear combination of <basis>.                            
TriangularBasis:=function(basis)
  local tribasis, goodwords, bad, allwords, count, pw, nonefound, w, b;
  # now we try and find unique identifers for each element of basis
  allwords:=[];
  for b in basis do
    Append(allwords, b.elts);
  od;
  allwords:=Set(allwords);
  # this list counts the number of different basis elements that that 
  # a given word appears in.
  count:=[1..Length(allwords)]*0;
  for b in basis do
    for w in [1..Length(b.elts)] do
      pw:=Position(allwords,b.elts[w]);
      count[pw]:=count[pw]+1;
    od;
  od;
  # to start with every basis element is bad
  bad:=[1..Length(basis)];
  # slowly, however, they become good with markers in goodwords
  tribasis:=[];
  goodwords:=[];
  repeat
    nonefound:=true;
    for b in bad do
      # If a given word w appears only once in the remaining basis
      # elements then count[pw]=1 and this word can be used as a
      # marker for the corresponding basis element.  We should,
      # perhaps, try to find such a basis element which has
      # coefficient 1 in basis[b], but...
      w:=Length(basis[b].elts)+1; # 'common' elements tend to come first
      repeat
	w:=w-1;
        pw:=Position(allwords,basis[b].elts[w]);
      until w=1 or count[pw]=1;
      if count[pw]=1 then
	Add(tribasis, basis[b]);
	Add(goodwords, allwords[pw]);
	RemoveSet(bad, b);
	nonefound:=false;
	# now adjust the multplicities in count...
	for w in [1..Length(basis[b].elts)] do
          pw:=Position(allwords,basis[b].elts[w]);
	  count[pw]:=count[pw]-1;
	od;
      fi;
    od;
  until bad=[] or nonefound;
  #if Length(bad)<>0 then
    #Print("Bad basis elements:\n");
    #for b in bad do Print("  ", b, "\n"); od;
    #Print("\n");
  #fi;
  if bad=[] then return [goodwords,tribasis];
  else return [goodwords,tribasis,bad];
  fi;
end;

#U RewriteWithRespectToBasis(elt, basis)                                  
#U RewriteWithRespectToBasis(elt, basis, tbasis)                          
#M This function tries to write <elt> as a linear combination of the      
#M elements in the (alleged) basis <basis>. It does this by using         
#M TriangularBasis(), unless <tbasis> is supplied,  to find a set of      
#M words which uniquely determine the coefficient of a given element of   
#M <basis> in <elt>, modulo higher terms.                                 
RewriteWithRespectToBasis:=function(arg)
  local elt, basis, tribasis, goodwords, nelt, coeffs, b, w, wb, coeff;
  elt:=arg[1];
  basis:=arg[2];
  if Length(arg)=3 then 
    goodwords:=arg[3][1];
    tribasis:=arg[3][2];
  elif Length(arg)=2 then 
    tribasis:=TriangularBasis(basis);
    goodwords:=tribasis[1];
    tribasis:=tribasis[2];
  else Error("usage: RewriteWithRespectToBasis(elt, basis [, tribasis])");
  fi;
  # hopefully we now have a triangular basis...
  nelt:=ShallowCopy(elt);
  coeffs:=[1..Length(basis)]*0;
  b:=1;
  while b<=Length(tribasis) and nelt<>0*nelt do
    w:=Position(nelt.elts, goodwords[b]);
    if IsInt(w) then 
      wb:=Position(tribasis[b].elts, goodwords[b]);
      coeff:=nelt.coeffs[w]/tribasis[b].coeffs[wb];
      coeffs[Position(basis, tribasis[b])]:=coeff;
      nelt:=nelt-coeff*tribasis[b];
    fi;
    b:=b+1;
  od;
  if nelt=0*nelt then 
    elt.basis:=basis;
    elt.basiscoeffs:=coeffs;
    return elt;
  else return nelt; fi;
end;

##########################################################################

#T Generic code for computing structure constants for a particular basis

#U MultiplicationTableBasis(basis)                                        
#M Work out, if possible, the multiplication table of the subalgebra of   
#M the group algebra RG which ha basis <basis>.                           
MultiplicationTableBasis:=function(basis) local mtable, tbasis, b, c;
  mtable:=[];
  tbasis:=TriangularBasis(basis);
  for b in [1..Length(basis)] do
    Add(mtable,[]);
    for c in [1..Length(basis)] do
      mtable[b][c]:=RewriteWithRespectToBasis(basis[b]*basis[c],basis,tbasis);
    od;
    Print(".\c");
  od;
  Print("\n");
  return mtable;
end;

#U DimensionOfRadicalFromStructureConstants(structure)                    
#M Given the matrix <structure> of the structure constants for some       
#M algebra compute the rank of the radical; by [Garsia-Reutenauer; 1.4]   
#M this is equal to the codimension of the matrix with (i,j)th entry equal
#M to the the trace of a_i*a_j where a_k=<structure>[k] is the matrix     
#M giving the action of the basis element b_k with respect to the basis   
#M <b_1,...,b_n>.                                                         
DimensionOfRadicalFromStructureConstants:=function(structure)
  local traces, i, j;
  traces:=List([1..Length(structure)], i->[]);
  for i in [1..Length(structure)] do
    for j in [i..Length(structure)] do
      traces[i][j]:=TraceMat(structure[i]*structure[j]);
      traces[j][i]:=traces[i][j];
    od;
  od;
  return Length(structure)-RankMat(traces);
end;

#U IsCommutativeFromBasis(basis)                                          
#M Use <basis> to check whether or not the corresponding algebra is       
#M commutative and return true if it is.                                  
IsCommutativeFromBasis:=function(basis)
  local comm, i, j;
  comm:=true;
  for i in [1..Length(basis)-1] do
    for j in [i+1..Length(basis)] do
      if basis[i]*basis[j]<>basis[j]*basis[i] then
	Print("(i,j) = (",i,",",j,")\n");
	comm:=false;
      fi;
    od;
    Print(".\c");
  od;
  Print("\n");
  return comm;
end;

