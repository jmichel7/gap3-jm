PrintPkgInit(rec(name:="pcqa",date:=[1995]));

PCQAexec := Concatenation( LOADED_PACKAGES.pcqa, "bin/pcqa" );
PCQAcpp := 0;
PCQAhom := 0;
PCQAsts := 0;

IsDerivedSeries := function( quot )
    return IsBound( quot.isDerivedSeries ) and quot.isDerivedSeries;
end;

ConvDOUBLE := function( seq )
local pos, string;

    string := "";
    for pos in [1..Length(seq)] do
        Append( string, String(seq[pos][1]));
        Append( string, " ");
        Append( string, String(seq[pos][2]));
        Append( string, " ");
    od;
    return string;
end;

ConvELM := function( seq )
local pos, string;

    string := "";
    for pos in [1..Length(seq)] do
        Append( string, String(seq[pos]));
        Append( string, " ");
    od;
    Append( string, "\n");
    return string;
end;

ConvFPG := function ( fpg )
    local string, g, r, rl, pos;

    string := "< ";

    for pos in [1..Length(fpg.generators)] do
        Append( string, fpg.namesGenerators[pos] );
        Append( string, ", " );
    od;
    Unbind( string[ Length(string) ] );
    Unbind( string[ Length(string) ] );
    
    Append( string, " | " );
    
    for r in fpg.relators do
        rl := Reversed(List( r ));
        for g in rl do
            pos := Position( fpg.generators, g );
            if pos <> false then
                Append( string, fpg.namesGenerators[pos] );
                Append( string, " " );
            else
                pos := Position( fpg.generators, g^-1 );
                Append( string, fpg.namesGenerators[pos] );
                Append( string, "^-1 " );
            fi;
        od;
        Unbind( string[ Length(string) ] );
        Append( string, ", " );
    od;
    if Length(fpg.relators) <> 0 then
        Unbind( string[ Length(string) ] );
    fi;
    Unbind( string[ Length(string) ] );
    Append( string, " >\n" );

    return string;
end;

ConvCPP := function ( cpp )
    local string, number, pos1, pos2, elm, id;

    string := "";
    number := cpp.Generators;
    id := 0*[1..number];
    Append( string, String(number));
    Append( string, "\n");
    for pos1 in [1..number] do
        Append( string, String(cpp.ExponentList[pos1]));
        Append( string, " ");
    od;
    Append( string, "\n");
    for pos1 in [1..number-1] do
        for pos2 in [pos1+1..number] do
            if IsBound(cpp.npRelations[pos1][pos2-pos1]) then
              elm := cpp.npRelations[pos1][pos2-pos1];
            else elm := id; fi;
            Append( string, ConvELM(Reversed(elm)));
        od;
    od;
    for pos1 in [1..number-1] do
        for pos2 in [pos1+1..number] do
            if IsBound(cpp.nnRelations[pos1][pos2-pos1]) then
              elm := cpp.nnRelations[pos1][pos2-pos1];
            else elm := id; fi;
            Append( string, ConvELM(Reversed(elm)));
        od;
    od;
    for pos1 in [1..number-1] do
        for pos2 in [pos1+1..number] do
            if IsBound(cpp.ppRelations[pos1][pos2-pos1]) then
              elm := cpp.ppRelations[pos1][pos2-pos1];
            else elm := id; fi;
            Append( string, ConvELM(Reversed(elm)));
        od;
    od;
    for pos1 in [1..number-1] do
        for pos2 in [pos1+1..number] do
            if IsBound(cpp.pnRelations[pos1][pos2-pos1]) then
              elm := cpp.pnRelations[pos1][pos2-pos1];
            else elm := id; fi;
            Append( string, ConvELM(Reversed(elm)));
        od;
    od;
    for pos1 in [1..Length(cpp.PowerRelations)] do
        if IsBound(cpp.PowerRelations[pos1]) then
          elm := cpp.PowerRelations[pos1];
          Append( string, ConvELM(Reversed(elm)));
        fi;
    od;
    return string;
end;

ConvHOM := function( homom )
    local pos, string, size;

  string := "";
  for pos in [1..Length(homom.Epimorphism)] do
      Append( string, ConvELM(Reversed(homom.Epimorphism[pos])));
  od;
  size := 0;
  for pos in [1..Length(homom.InverseMap)] do
      size := size+Length(homom.InverseMap[pos])+1;
  od;
  Append( string, String(size));
  Append( string, "\n");
  for pos in [1..Length(homom.InverseMap)] do
      Append( string, ConvDOUBLE(Reversed(homom.InverseMap[pos])));
      Append( string, "0\n");
  od;
  return string;
end;

OpCPP := function(cpp)
    local pos1, pos2, number, temp;

    temp := cpp.npRelations;
    cpp.npRelations := cpp.ppRelations;
    cpp.ppRelations := temp;
    temp := cpp.nnRelations;
    cpp.nnRelations := cpp.pnRelations;
    cpp.pnRelations := temp;
    number := cpp.Generators;
    for pos1 in [1..number-1] do
        for pos2 in [pos1+1..number] do
            if IsBound(cpp.ppRelations[pos1][pos2-pos1]) then
              cpp.ppRelations[pos1][pos2-pos1] :=
                Reversed(cpp.ppRelations[pos1][pos2-pos1]);
            fi;
            if IsBound(cpp.pnRelations[pos1][pos2-pos1]) then
              cpp.pnRelations[pos1][pos2-pos1] :=
                Reversed(cpp.pnRelations[pos1][pos2-pos1]);
            fi;
            if IsBound(cpp.npRelations[pos1][pos2-pos1]) then
              cpp.npRelations[pos1][pos2-pos1] :=
                Reversed(cpp.npRelations[pos1][pos2-pos1]);
            fi;
            if IsBound(cpp.nnRelations[pos1][pos2-pos1]) then
              cpp.nnRelations[pos1][pos2-pos1] :=
                Reversed(cpp.nnRelations[pos1][pos2-pos1]);
            fi;
        od;
    od;
    for pos1 in [1..Length(cpp.PowerRelations)] do
        if IsBound(cpp.PowerRelations[pos1]) then
          cpp.PowerRelations[pos1] := Reversed(cpp.PowerRelations[pos1]);
        fi;
    od;
end;

OpHOM := function(homom)
    local pos;

    for pos in [1..Length(homom.Epimorphism)] do
        homom.Epimorphism[pos] := Reversed(homom.Epimorphism[pos]);
    od;
    for pos in [2..Length(homom.InverseMap)] do
        homom.InverseMap[pos] := Reversed(homom.InverseMap[pos]);
    od;
end;

ExtendPCQA := function ( arg )
    local temp1, temp2, temp3, temp4, temp5, temp6, temp7, string, result, next,
          fpg, cpp, homom, flag, cent;

    if Length(arg) > 5 or Length(arg) < 3 then return;
    fi;
    fpg := arg[1];
    cpp := arg[2];
    homom := arg[3];
    if Length(arg) >= 4 then flag := arg[4];
    else flag := 0; fi;
    if Length(arg) = 5 then cent := arg[5];
    else cent := 0; fi;
    temp1 := TmpName();
    PrintTo( temp1, ConvFPG(fpg) );
    temp2 := TmpName();
    PrintTo( temp2, ConvCPP(cpp) );
    temp3 := TmpName();
    PrintTo( temp3, ConvHOM(homom) );
    temp4 := TmpName();
    PrintTo( temp4, "" );
    temp5 := TmpName();
    PrintTo( temp5, "" );
    temp6 := TmpName();
    PrintTo( temp6, "" );
    temp7 := TmpName();
    PrintTo( temp7, "" );

    string := "r\n1\n";
    Append( string, temp1);
    Append( string, "\n3\n");
    Append( string, temp2);
    Append( string, "\n4\n");
    Append( string, temp3);
    Append( string, "\nq\n2\n1\nq\n3\n");
    if flag > 0 then
        Append( string, "9\n1\n");
        Append( string, String(flag));
        Append( string, "\nq\n");
    fi;
    if cent <> 0 then Append( string, "9\n4\nq\n");
    fi;
    Append( string, "g\n0\nq\n4\nc\nG\n1\nq\n");
    if flag < 0 then Append( string, "7\n3\ng\n0\nq\n4\n"); fi;
    Append( string, "2\n");
    Append( string, temp4);
    Append( string, "\ny\nq\ny\n4\n3\ny\n5\n3\nq\ns\n3\n");
    Append( string, temp5);
    Append( string, "\ny\n4\n");
    Append( string, temp6);
    Append( string, "\ny\nq\nq\ny\n");
    PrintTo( temp7, string);
    Exec(Concatenation(PCQAexec, " < ", temp7, " > /dev/null"));
    Read( temp4 );
    if PCQAsts = -1 then
        result := rec ( QuotientStatus := PCQAsts );
    elif PCQAsts > 0 then
        result := rec ( QuotientStatus := PCQAsts );
    else
        Read( temp5 );
        Read( temp6 );
        Exec( Concatenation("rm ", temp5));
        Exec( Concatenation("rm ", temp6));
        OpCPP(PCQAcpp);
        OpHOM(PCQAhom);
        result := rec (
            QuotientStatus := PCQAsts,
            PolycyclicPresentation := PCQAcpp,
            Homomorphisms := PCQAhom,
            Next := cpp.Generators+1
        );
    fi;
    Exec( Concatenation("rm ", temp1));
    Exec( Concatenation("rm ", temp2));
    Exec( Concatenation("rm ", temp3));
    Exec( Concatenation("rm ", temp4));
    Exec( Concatenation("rm ", temp7));
    return result;
end;

CallPCQA := function(fpg , n)
    local i, temp1, temp2, temp3, temp4, string, result, arr;
 
    if n >= 1 then
        temp1 := TmpName();
        PrintTo( temp1, ConvFPG(fpg) );
        temp2 := TmpName();
        PrintTo( temp2, "" );
        temp3 := TmpName();
        PrintTo( temp3, "" );
        temp4 := TmpName();
        PrintTo( temp4, "" );
        string := "1\n";
        i := 1;
        Append( string, temp1);
        Append( string, "\nc\nsmith\n1\nG\n1\nq\n5\n3\nq\ns\n3\n");
        Append( string, temp3);
        Append( string, "\ny\n4\n");
        Append( string, temp4);
        Append( string, "\ny\nq\nq\ny\n");
        PrintTo( temp2, string);
        Exec( Concatenation(PCQAexec, " < ", temp2," > /dev/null"));
        Read( temp3 );
        Read( temp4 );
        PCQAsts := 0;
        OpCPP(PCQAcpp);
        OpHOM(PCQAhom);
        arr := [PCQAcpp.Generators];
        i := 2;
        while i <= n and PCQAsts = 0 do
            result := ExtendPCQA( fpg, PCQAcpp, PCQAhom );
            if result.QuotientStatus = 0 then
                Add(arr, PCQAcpp.Generators);
                PCQAcpp := result.PolycyclicPresentation;
                PCQAhom := result.Homomorphisms;
                i := i+1;
            fi;
        od;
        Exec( Concatenation("rm ", temp1));
        Exec( Concatenation("rm ", temp2));
        Exec( Concatenation("rm ", temp3));
        Exec( Concatenation("rm ", temp4));
        result := rec (
            isDerivedSeries := true,
            DerivedLength := i-1,
            QuotientStatus := PCQAsts,
            PolycyclicPresentation := PCQAcpp,
            Homomorphisms := PCQAhom,
            MembershipArray := arr
        );
        return result;
    fi;
end;

AbelianComponent := function( quot )
    local i, j, cmp, mat, cpp;

    mat := [];
    cpp := quot.PolycyclicPresentation;
    if IsDerivedSeries(quot) then j := quot.MembershipArray[1];
    else return; fi;
    for i in [1..j] do
        if IsBound(cpp.PowerRelations[2*i-1]) then
            Append(mat,[Sublist(cpp.PowerRelations[2*i-1],[1..j])]);
        else Append(mat,[0*[1..j]]);
        fi;
        mat[i][i] := -cpp.ExponentList[i];
    od;
    cmp := [ElementaryDivisorsMat(mat)];
    j := j+1;
    for i in [2..quot.DerivedLength] do
        Append( cmp,
            [Sublist(cpp.ExponentList,
                [j..quot.MembershipArray[i]])]);
        j := quot.MembershipArray[i]+1;
    od;
#    for i in [1..Length(cmp)] do
#        cmp[i] := ElementaryDivisorsOfList(cmp[i]);
#    od;
    return cmp;
end;

HirschLength := function( cpp )
    local i, leng;

    leng := 0;
    for i in [1..cpp.Generators] do
        if cpp.ExponentList[i] = 0 then leng := leng+1; fi;
    od;
    return leng;
end;

ModuleAction := function( quot )
    local i, j, size, dl, ft, lt, matlist, mat, cpp, flag;

    cpp := quot.PolycyclicPresentation;
    flag := false;
    matlist := [];
    if IsDerivedSeries(quot) then
        dl := quot.DerivedLength;
        if dl > 1 then
            ft := quot.MembershipArray[dl-1];
            lt := quot.MembershipArray[dl];
            flag := true;
        fi;
    else
        ft := quot.Next-1;
        lt := cpp.Generators;
      flag := true;
    fi;
    if flag then
        size := lt-ft;
        for i in [1..ft] do
            mat := [];
            for j in [ft-i+1..lt-i] do
                Add( mat, Sublist(cpp.ppRelations[i][j],[ft+1..lt]));
            od;
            Add(matlist, Copy(mat));
        od;
        return matlist;
    fi;
end;



