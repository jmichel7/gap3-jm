#############################################################################
##
#A  Matrix package                                               Anthony Pye 
##                                                      
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##
###############################################################################
##
##  RandomRelation(<G>,<P>,<maps>) . .  Given a matrix group <G>,
##
##  with permutation representation <P> return a random relation.
##
RandomRelation := function(arg)
    local i, G, P, genimages, rel;
    
    # check for the correct number of arguments
    if Length(arg) <> 2 and Length(arg) <> 3 then
        Error(" Inproper use of RandomRelation ");
    fi;
    G := arg[1];
    P := arg[2];
    
    # set the permutation group to be P
    SetPermGroupPFlag(G,P);
    SetPermDomainFlag(G,[]);
    
    # deal with maps if they are given as an argument
    if Length(arg) = 3 then
        genimages := [];
        for i in arg[3] do
            if i = 0 then
                Add(genimages,P.1^0);
            else
                Add(genimages,Generators(P)[i]);
            fi;
        od;
    else
        genimages := Generators(P);
    fi;
    G.permGroupP.operation := rec(genimages := genimages);
    
    # get the relation using the function 'RandomRelsPerm'
    rel := RandomRelsPerm(G);
    
    # unbind information we don't need anymore from <G>
    UndoPermGroupPFlag(G);
    UndoPermDomainFlag(G);
    UndoFpHomomorphismFlag(G);
    UndoAbstractGeneratorsFlag(G);
    
    # and return the relation
    return rel;
    
end;

###############################################################################
##
##  EvaluateRelation(<reln>,<G>) . . . . . . . . evaluate the relation
##
##  <reln> on the generators of <G>.
##
EvaluateRelation := function(reln,G)
    
    return MappedWord(reln,Generators(FpGroupFlag(G)),GeneratorsFlag(G));
    
end;
