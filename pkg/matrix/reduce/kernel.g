#############################################################################
##
#A  Matrix package                                               Anthony Pye 
##                                                      
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
## 
ReduceLayerNumberByOne := function(step)
    
    SetLayerNumberFlag(step,LayerNumberFlag(step)-1);
    if KernelFlag(step)="unknown" then
        return;
    fi;
    ReduceLayerNumberByOne(KernelFlag(step));
    
end;

##########################################################################
##
##  ApproximateKernel( <G>, <P>, <m>, <n>, <maps> ) . .  approximate the kernel
##
##  of the mapping from the matrix group <G> to the permutation group <P>.
##  <m> is the maximum number of relations to be stripped and <n> is the
##  successive number of relations to be stripped during the course of the
##  computation.
##  <maps> is an additional parameter to indicate the correspondence 
##  between the generators of <G> and <P>, if not given we assume the
##  obvious one.
##
ApproximateKernel := function(arg)
    local G, P, step, genimages, i;
    
    # check for number of arguments
    if Length(arg) <> 4 and Length(arg) <> 5 then
        Error(" Bad use of ApproximateKernel: Wrong number of arguments");
    fi;
    G := arg[1];
    P := arg[2];
    
    # check if the group is irreducible
    if IsMatGroup(G) then
        G := GModule(G);
    fi;
    if IsIrreducible(G)=false then
        Error("Bad use of ApproximateKernel: G must be irreducible");
    fi;
    
    # check that <m> and <n> are set sensibly
    if arg[4] > arg[3] then 
        Error("Bad use of Approximate Kernel: n > m");
    fi;
    
    # initialise <step.quotient>
    step := InitSplit(G,arg[3],arg[4]);
    SetTypeFlag(step,"Perm");
    SetQuotientFlag(step,G);
    SetIdentityQuotientFlag(step,IdentityFlag(step));
    Print ("#W  Applying Size to permutation group\n");
    SetSizeQuotientFlag(step,Size(P));
    SetDimensionQuotientFlag(step,DimensionFlag(step));
    SetGeneratorsFlag(step,GeneratorsFlag(G));
    SetBasisSubmoduleFlag(step,[]);
    SetBasisFlag(step,IdentityFlag(step));
    
    # set kernel
    InitialiseKernel(step);
    
    # assuming the quotient is irreducible
    SetReducibleFlag(QuotientFlag(step),false);
    
    # store the permutation group information in the record
    SetPermGroupPFlag(QuotientFlag(step),P);
    SetPermDomainFlag(QuotientFlag(step),[]);
    SetIsFaithfulFlag(QuotientFlag(step),true);
    
    # deal with correspondence between generators
    if Length(arg) = 5 then
        genimages := [];
        for i in arg[5] do
            if i = 0 then
                Add(genimages,P.1^0);
            else
                Add(genimages,Generators(P)[i]);
            fi;
        od;
    else
        genimages := Generators(P);
    fi;
    step.quotient.permGroupP.operation := rec(genimages := genimages);
    
    # and into the code we go ..
    GoDownGeneratorsRels(step,GeneratorsFlag(G),GeneratorsFlag(G));
    GoDownChain(step,[]);
    
    # subtract 1 from each of the levels so that we start from level 0
    ReduceLayerNumberByOne(step);
    
    # return what we have found
    return step.kernel;
    
end;
