#############################################################################
##
#A  Matrix package                                               Anthony Pye 
##                                                      
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##
BaseRing := function(G)
    if IsBound(G.field)=false then return "unknown"; fi;
    return G.field;
end;

RowSpaceOps.\^ := function(space,mat)
    return RowSpace(space.field,Base(space)*mat);
end;

GetFpHomomorphism := function(g)
    local f, h, op_gens;
    
    op_gens := PermGroupPFlag(g).operation.genimages;
    f := FreeGroup(Length(op_gens));
    h := GroupHomomorphismByImages(PermGroupPFlag(g),f,
                 op_gens,Generators(f));
    h.isHomomorphism := true;
    h.isMapping := true;
    h.isBijection := true;
    SetFpHomomorphismFlag(g,h);
    SetFpGroupFlag(g,f);
    SetAbstractGeneratorsFlag(g,Generators(FpGroupFlag(g)));
    if Length(GeneratorsFlag(g)) <> Length(AbstractGeneratorsFlag(g)) then
        Error("Something has gone wrong in 'GetFpHomomorphism' ");
    fi;
    
end;

###############################################################################
##
#F  RandomRelsPerm( <g> ) . . . . . . . . . . construct random relators
##
##  'RandomRelsPerm' will return a "random" relator for <g> as word in
##  <g.abstractGenerators> using a permutation representation.
##
RandomRelsPerm := function(g)
    local w, p;
    
    # construct a random relation
    SetupPermRep(g);
    w := Product([1..Random([1..100])],x->Random(Generators(FpGroupFlag(g))));
    p := MappedWord(w,Generators(FpGroupFlag(g)),PermGroupPFlag(g).operation.genimages);
    return w/Image(FpHomomorphismFlag(g),p);
    
end;
