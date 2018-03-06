#############################################################################
##
#A  fpgtorws.g                  GAP library                  Derek Holt
##
##  This file contains a function to calculate a rewriting-system from
##  a finitely presented group.
##  It is derived from such a function for group presentations
##  written by Werner Nickel.
##
##
#############################################################################
##
#F  FpGroupToRWS(<G>[,<case_change>]). calculate a rewriting-system for group G
##
##  <G> must be a finitely presented group. This function calculates and
##  returns a rewriting-system for <G>.
##  By default, the string used for inverse of generator with printing string
##  "g" is "g^-1". If <case_change> is set to be true, then inverse strings
##  are derived by changing case.
##
##  Public function.
FpGroupToRWS := function ( arg )
    local    G, case_change, rws, gens, names, invgens, ggens, gnames, grels,
             lhs, rhs, g, l, hl, ng, r, i, j, k, ct, istr;

    if Number(arg)<1 or Number(arg)>2 then
       Error("Number of arguments must be 1 or 2");
    fi;
    G := arg[1];
    if Number(arg)=2 and arg[2]=true then case_change:=true;
    else case_change := false;
    fi;
    # check the given argument to be an fp-group.
    if not IsFpGroup(G) then
       Error("Argument must be a finitely presented group");
    fi;

    rws := rec();
    rws.isRWS  := true;
    rws.isInternalRWS  := true;
    rws.operations := RWSOps;
    rws.gpMon := G;

    rws.generators := [];
    gens := rws.generators;
    rws.namesGenerators := [];
    names := rws.namesGenerators;
    rws.invGenerators := [];
    invgens := rws.invGenerators;
    rws.ordering := "shortlex";
    
    ggens := G.generators;
    if IsBound(G.namesGenerators)  then
        gnames := G.namesGenerators;
    else
        gnames := [];
        for i  in [ 1 .. Length(G.generators) ]  do
            Add( gnames, Concatenation( "_x", String(i) ) );
        od;
    fi;
    ng := Length(ggens);
    if IsBound(G.relators) then grels := ShallowCopy(G.relators);
    else grels := [];
    fi;

    ct := 0;
    for i in [1..ng] do
       Add(gens,ggens[i]);
       Add(names,gnames[i]);
       ct := ct+1;
       #Check to see if this generator is involutory.
       #If so, remove all occurrences of its inverse from the relators.
       #If not, adjoin an inverse generator to gens.
       if Position(grels,ggens[i]^2)<>false or
          Position(grels,ggens[i]^-2)<>false  then
          for j in [1..Length(grels)] do
             k := 1;
             while k <= LengthWord(grels[j]) do
                if Subword(grels[j],k,k)=ggens[i]^-1 then
                   grels[j] := SubstitutedWord(grels[j],k,k,ggens[i]);
                   k := 1;
                else k := k+1;
                fi;
             od;
          od;
          Add(rws.invGenerators,ct);
       else
         Add(gens,ggens[i]^-1);
         if case_change then
            istr := CaseChange(gnames[i]);
            if Position(gnames,istr)<>false then
               Error("Case-change for inverses causes repetition of name");
            fi;
            Add(names,istr);
         else
            Add(names,Concatenation(gnames[i],"^-1"));
         fi;
         Append(rws.invGenerators,[ct+1,ct]);
         ct := ct+1;
       fi;
    od;

    ng := ct;
    rws.numGenerators := ng;
    rws.generatorOrder := [1..ng];
    for i in [1..ng] do rws.(i) := gens[i]; od;
    
    rws.equations := [];
    for r in grels do
       if r <> IdWord then
         #balance the equation
         l := LengthWord(r); hl := Int((l+1)/2);
         lhs := Subword(r,1,hl); rhs := Subword(r,hl+1,l)^-1;
         #rhs may have acquired some redundant inverses again!
         for k in [1.. LengthWord(rhs)] do
            g := Subword(rhs,k,k);
            if Position(gens,g) = false then
               rhs := SubstitutedWord(rhs,k,k,g^-1);
            fi;
         od;
         if lhs <> rhs then
             Add(rws.equations,[rws.operations.WordToList(lhs,gens),
                            rws.operations.WordToList(rhs,gens)]);
         fi;
       fi;
    od;

    return rws;
end;

#############################################################################
##
#F  CaseChange(<string>) . . . . . . . change case of all letters in string
##
CaseChange := function ( string )
    local  lc, uc, newstring, pl, pu, i;
    lc := "abcdefghijklmnopqrstuvwxyz";
    uc := "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    newstring := [];
    for i in [1..Length(string)] do
      pl := Position(lc,string[i]);
      pu := Position(uc,string[i]);
      if pl<>false then newstring[i]:=uc[pl];
      elif pu<>false then newstring[i]:=lc[pu];
      else newstring[i]:=string[i];
      fi;
    od;

    return (String(newstring));
end;
