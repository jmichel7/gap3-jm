###########################################################################
##
#A  misc.g                   autag package                 Michael J Smith
##
##  November 1996
##
##  This file is part of the autag package, which contains GAP code for
##  computing automorphism groups of finite soluble groups given by
##  special soluble group presentations.
## 
###########################################################################


## miscellaneous function definitions

IsSagGroup := function (G)
    return IsBound(G.isSagGroup) and G.isSagGroup = true;
end;

Select := function (cond, iftrue, iffalse)
    # if <cond> true, return <iftrue>, else return <iffalse>
    if cond then return iftrue; else return iffalse; fi;
end;


PrintVec := function (vec)
    Print("    ",List(vec, x -> Int(x)));
    if IsFFE(vec[1]) then Print(" * ", vec[1]^0); fi;
end;


PrintMat := function (mat)
    local i;
    Print(" [ ");
    if Length(mat) > 0 then
        Print(List(mat[1], x -> Int(x)));
        for i in [2..Length(mat)] do
            Print(",\n   ");
            Print(List(mat[i], x -> Int(x)));
        od;
        Print(" ]");
        if IsFFE(mat[1][1]) then Print(" * ",mat[1][1]^0); fi;
    else
        Print(" [ [ ] ]");
    fi;
end;

NiceMatList := function (mats)
    #
    # <mats> should be a list of matrices of the same dimensions over
    # either the same prime field, or over the integers. This function
    # displays them out in a compact form. Note that <mats> must be
    # nonempty.

    local elem_width, flat, mat_width, screen_width, per_row, first, last, 
          i, k, j, done;

    # work out <elem_width>, the max width of the elements as strings
    # (NB includes the space between columns)
    if IsFFE(mats[1][1][1]) then
        Print("  Field: ", Field(mats[1][1][1]), "\n");
        elem_width := Length(String(Size(Field(mats[1][1][1]))));
    else
        flat := Flat(mats);
        elem_width := Maximum(Length(String(Maximum(flat))),
                              Length(String(Minimum(flat))));
    fi;
    mat_width := 6 + Length(mats[1][1])*(elem_width+1);
    screen_width := SizeScreen()[1];
    per_row := Maximum(1, Int((screen_width-2)/mat_width));

    first := 1; last := Minimum(per_row, Length(mats));
    repeat 
        for i in [1..Length(mats[1])] do    # rows
            for k in [first..last] do       # matrices in this block 
                Print("    [");
                for j in [1..Length(mats[k][i])] do
                    Print(String(Int(mats[k][i][j]), elem_width));
                    if j < Length(mats[k][i]) then Print(" "); fi;
                od;
                Print("]");
            od;
            Print("\n");
        od;
        first := last + 1;
        done := last = Length(mats);
        last := Minimum(first - 1 + per_row, Length(mats));
        if not done then Print("\n"); fi;
    until done;
end;


Nice := function (arg)
    local i, j;

    # figure out whether to call NiceMatList -- we need to know if there
    # is a single argument which is a list of matrices all of the same type
    if Length(arg) = 1 and
       IsList(arg[1]) and
       Length(arg[1]) > 0 and
       IsMat(arg[1][1]) and
       ForAll(arg[1], a -> IsMat(a)) and
       Length(Set(List(arg[1], a -> DimensionsMat(a)))) = 1 and
       IsVector(Flat(arg[1])) then
        NiceMatList(arg[1]);

    else
        for i in [1..Length(arg)] do
            if IsMat(arg[i]) then 
                NiceMatList([arg[i]]);
            elif IsVector(arg[i]) then
                PrintVec(arg[i]);
            elif IsList(arg[i]) then
                Print("[\n");
                for j in [1..Length(arg[i])] do
                    if IsBound(arg[i][j]) = false then
                        Print("\n");
                    else
                        Nice(arg[i][j]);
                    fi;
                    if j < Length(arg[i]) then Print(","); fi;
                    Print("\n");
                od;
                Print("]\n");
            else
                Print(arg[i]);
            fi;
            if i < Length(arg) then Print(","); fi;
        od;
    fi;
end;

Plural := function (arg)
    if arg[1] = 1 then return ""; else return "s"; fi;
end;

Seconds := function (rt)
    local ans;
    ans := Cat(String(Int(rt/1000)),".");
    if rt mod 1000 = 0 then
        ans := Cat(ans, "000");
    else
        ans := Cat(ans, List([1..3 - Length(String(rt mod 1000))], x -> '0'),
                   String(rt mod 1000));
    fi;
    return ans;
end;

Time := function (rt)
    return ConcatenationString( StringTime(rt), " seconds" );
end;

IsNonZeroElt := function (x)
    return x <> x * 0;
end;

IsZeroElt := function (x)
    return x = x * 0;
end;

RenamedGensSagGroup := function (G, name)
    # the following is needed to get better generator names...
    return SpecialAgGroup(AgGroupFpGroup(G.operations.FpGroup(G,name)));
end;
