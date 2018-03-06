#############################################################################
##
#F  WeightDistribution( <C> ) . . . returns the weight distribution of a code
##
WeightDistribution := function (C)
    if not IsBound(C.weightDistribution) then
        C.weightDistribution := C.operations.WeightDistribution(C);
    fi;
    return C.weightDistribution;
end;

CodeOps.WeightDistribution := function (C)
    local El, nl, newwd;
    if IsLinearCode(C) then
        return WeightDistribution(C);
    fi;
    El := VectorCodeword(Elements(C));
    nl := VectorCodeword(NullWord(C));
    newwd := DistancesDistributionVecFFEsVecFFE(El, nl);
    if Length( newwd ) > WordLength(C)+1  then
        newwd := newwd{[1..WordLength(C)+1]};
    fi;
    return newwd;
end;

LinCodeOps.WeightDistribution := function(C)
    local G, nl, k, n, q, wd, newwd, oldrow, newrow, i, j;
    n := WordLength(C);
    k := Dimension(C);
    q := Size(Field(C));
    nl := VectorCodeword(NullWord(C));
    if k = 0 then
        G := NullVector(n+1);
        G[1] := 1;
        newwd := G;
    elif k = n then
        newwd := List([0..n], i->Binomial(n, i));
    elif k <= Int(n/2) then
	G := Copy(GeneratorMat(C));
        newwd := DistancesDistributionMatFFEVecFFE(G, q, nl);
    else
        G := Copy(CheckMat(C));
        wd := DistancesDistributionMatFFEVecFFE(G, q, nl);
        Print(wd);
        newwd := [Sum(wd)];
        oldrow := List([1..n+1], i->1);
        newrow := [];
        for i in [2..n+1] do
            newrow[1] := Binomial(n, i-1) * (q-1)^(i-1);
            for j in [2..n+1] do
                newrow[j] := newrow[j-1] - (q-1) * oldrow[j] - oldrow[j-1];
            od;
            newwd[i] := newrow * wd;
            oldrow := Copy(newrow);
        od;
        newwd:= newwd / (q ^ Redundancy(C));
    fi;
    Print(newwd);
    
    if Length(newwd)>n+1 then
        newwd:=newwd{[1..n+1]};
    fi;
    return newwd;
end;

CycCodeOps.WeightDistribution := LinCodeOps.WeightDistribution;
