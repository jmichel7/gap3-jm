#############################################################################
##
#A  rws.g                  GAP library                  Derek Holt
##
##  This file contains those functions that deal with rewriting systems.
##
##  15.3.95.
##  Each (internal) rewriting-system now has components "gpMon" (for the
##  associated group or monoid), "generators" (for the generators, which
##  will include those of gpMon, but may also include inverses
##  in the group case), and "namesGenerators", which again include
##  those of gpMon, but will have names with "^-1" adjoined for inverses. 
##
##  When an externally created file containing a rewriting-system is read in
##  to GAP, a preprocessing external program called "ppgap" is run, which
##  creates a file called "file.gap", which includes necessary declarations
##  of a suitable underlying monoid.
##
##  23. 2. 95. The internal storage of a rewriting system was changed so
##  that generators are simply numbers in the range [1..ng] for some ng,
##  and words are lists of generator numbers.
##
if not IsBound(InfoRWS)  then InfoRWS := Ignore;  fi;

#############################################################################
#V  _RWS                external variable - the name of the rewriting system
#V  _RWS.gpMon 		external variable - name of underlying group or monoid
#V  _KBExtDir 		external variable - directory of external executables
#V  _KBTmpFileName      external variable - name of temporary file.
##
_RWS   := rec();
_RWS.gpMon := rec();

#_KBExtDir  :=  "/maths/geom/dfh/kbmag/bin/";
_KBTmpFileName := TmpName();

#############################################################################
##
#F  RWSOps  . . . . . . . operations record for RWS
##
RWSOps := rec();
#############################################################################
##
#F  IsRWS(<x>) . . . . . . . test whether x is an RWS
##
##  Public function.
IsRWS := function ( x )
    return IsRec(x) and IsBound(x.isRWS) and x.isRWS;
end;

#############################################################################
##
#F  IsInternalRWS(<x>) . . . . . . . test whether x is an InternalRWS
##
##  Public function.
IsInternalRWS := function ( x )
    return IsRec(x) and IsBound(x.isInternalRWS) and x.isInternalRWS;
end;

#############################################################################
##
#F  IsConfluentRWS(<x>) . . . . . . . test whether x is a confluent rws
##
##  Public function.
IsConfluentRWS := function ( x )
    return IsRec(x) and IsBound(x.isInternalRWS) and IsBound(x.isConfluent)
	   and x.isConfluent=true;
end;

#############################################################################
##
#F  IsAvailableNormalFormRWS(<x>) . . . . test whether x has a normal form
##
##  Public function.
IsAvailableNormalFormRWS := function ( x )
    return IsRec(x) and IsBound(x.isInternalRWS) and
           IsBound(x.isAvailableNormalForm) and x.isAvailableNormalForm=true;
end;

#############################################################################
##
#F  IsAvailableReductionRWS(<x>) . . test whether x has a reduction algorithm
##
##  Public function.
IsAvailableReductionRWS := function ( x )
    return IsRec(x) and IsBound(x.isInternalRWS) and
           IsBound(x.isAvailableReduction)
                 and x.isAvailableReduction=true;
end;

#############################################################################
##
#F  SetOrderingRWS(<rws>,<ord>[,list]) 
##                          . . .  specify the ordering of a rewriting system
##
##  <rws> should be a rewriting system, and <ord> one of the strings that
##  defines an ordering on the words in the generators of <rws>.
##  These are "shortlex", "recursive", rt_recursive", "wtlex" and "wreathprod".
##  In the case of "wtlex" and "wreathprod", the extra parameter <list> is
##  required, and it should be a list of ng (= number of generators of <rws>)
##  non-negative integers, specifying the weights or the levels of the
##  generators, respectively, for this ordering.
##  Public function.
SetOrderingRWS := function ( arg )
    local rws, ord, list, ng, go, i;
    if Length(arg)<2 or Length(arg)>3 then
       Error("SetOrderingRWS has 2 or 3 arguments");
    fi;
    rws:=arg[1];
    ord:=arg[2];
    if Length(arg)=3 then
      list:=arg[3];
    else
      list:=[];
    fi;
    if not IsInternalRWS(rws)  then
       Error("First argument is not an internal rewriting system.");
    fi;
    if not IsInternalRWS(rws)  then
       Error("Second argument is not a string.");
    fi;
    ng := rws.numGenerators;
    if Length(arg)=3 then
      if not IsList(list) or Length(list)<>ng then
         Error("Third argument is not a list of length <ng>.");
      fi;
      for i in [1..ng] do
        if not IsInt(list[i]) or list[i]<0 then
          Error("Third argument is not a list of non-negative integers.");
        fi;
      od;
    fi;
    
    if ord="shortlex" or ord="recursive" or ord="rt_recursive" or
       ord="wtlex" or ord="wreathprod" then
      rws.ordering:=ord;
    else
      Error("Unknown ordering",ord);
    fi;
    if (ord="wtlex" or ord="wreathprod") and list=[] then
      Error("Third argument required for ordering",ord);
    fi;
    if ord="wtlex" then
      go:=rws.generatorOrder;
      rws.weight:=[];
      for i in [1..ng] do
        rws.weight[go[i]]:=list[i];
      od;
    fi;
    if ord="wreathprod" then
      go:=rws.generatorOrder;
      rws.level:=[];
      for i in [1..ng] do
        rws.level[go[i]]:=list[i];
      od;
    fi;
end;

#############################################################################
##
#F  ReorderGeneratorsRWS(<rws>,<p>) . reorder generators of a rewriting system
##
##  <rws> should be a rewriting system, and <p>  a permutation of the set
##  [1..ng], where <rws> has <ng> generators.
##  The generators of <rws> are reordered by applying the permutation <p> to
##  their existing order.
##  Public function.
ReorderGeneratorsRWS := function ( rws, p )
    local ng, ngo, i;
    if not IsInternalRWS(rws)  then
       Error("First argument is not an internal rewriting system.");
    fi;
    if not IsPerm(p) then
       Error("Second argument is not a permutation.");
    fi;
    ng := rws.numGenerators;
    if LargestMovedPointPerm(p) > ng then
       Error("Permutation is on more points than there are generators!");
    fi;
    
    ngo:=[];
    for i in [1..ng] do
      ngo[i^p]:=rws.generatorOrder[i];
    od;
    rws.generatorOrder:=ngo;
end;

#############################################################################
##
#F  ReadRWS(<filename>) . . . . . . . .read and convert a rewriting system
##
##  ReadRWS reads a rewriting system, which must be declared with the
##  external variable name "_RWS" from the file <filename>, and converts it
##  to internal format. First it creates and reads the GAP preproceesor file
##  <filename>.gap, for the declarations of variable names.
##  This is created using the external program "ppgap".
##  The rewriting system is returned.
##  Public function.
ReadRWS := function ( filename )
    local filename, gens, names, i, invgens, nig, eqn, ng, go, ri, ig, rgm, rws;
    
    Exec(Concatenation(_KBExtDir,"ppgap ",filename));
    Read(Concatenation(filename,".gap"));
    Exec(Concatenation("/bin/rm -f ",filename,".gap")); 
 
    rgm := _RWS.gpMon; #This is about to get overwritten, so we remember it! 

    if not READ(filename) then Error("Could not open file for reading"); fi;

    rws := _RWS; _RWS := rec(); #reset _RWS for next time
    rws.isInternalRWS := true;
    rws.operations := RWSOps;
    
    rws.gpMon := rgm;
    rws.generators := Copy(rgm.generators);
    gens := rws.generators;
    rws.namesGenerators := Copy(rgm.namesGenerators);
    names := rws.namesGenerators;
    ng := Length(gens);
    go := rws.generatorOrder;
    nig := 0;
    if IsFpGroup(rgm) then
      invgens := List( [1..ng], i->gens[i]^-1 );
      for i in [1..ng] do
        if Position(go,invgens[i])<>false then
           Add(gens,invgens[i]);
           nig:=nig+1;
           Add(names,Concatenation(names[i],"^-1"));
        fi;
      od;
      ng := ng+nig;
    fi;
    
    if ng<>Length(go) then
       Error("Length of generator order field wrong");
    fi;
    rws.numGenerators :=  ng;

    for i in [1..ng] do rws.(i) := gens[i]; od;

    for i in [1..ng] do go[i] := Position(gens,go[i]); od;

    ri := rws.inverses;
    ig := []; rws.invGenerators := ig;
    for i in [1..ng] do
      if IsBound(ri[i]) then ig[go[i]] := Position(gens,ri[i]);
      else ig[go[i]] := false;
      fi;
    od;
    Unbind(rws.inverses);

    for eqn in rws.equations do
        eqn[1] := rws.operations.WordToList(eqn[1],gens);
        eqn[2] := rws.operations.WordToList(eqn[2],gens);
    od;

    return rws;
end;

#############################################################################
##
#F  RWSOps.UpdateRWS(<rws>, <filename>, <kb>)
##                              . . update rws, after run of external program
##
## <rws> should be a rewriting system on which an external program has been
## run (using KB or Autgroup). <filename> should be the file in which the
## original rewriting-system was stored. <kb> should be true or false,
## according to whether the function is being called from a Knuth-Bendix
## application (KB) or automatic groups (Automata).
## This function updates the rewriting system <rws> after the call of the
## external program.
## In the Knuth-Bendix case, <filename>.kbprog is read in, for the updated
## version of the equations. Then <filename>.reduce is read in
## for the reduction machine.
## In the automatic groups case, <filename>.wa is read in for the word-acceptor,
## and then <filename>.diff1c for the word-difference machine used in word
## reduction.
## Unfortunately, the generator order may be different internally and
## externally, so the states of the relevant finite state automata may need
## to be permuted, to comply with the internal ordering.
## 

RWSOps.UpdateRWS := function(rws, filename, kb)
	local gens, x, i, j, k, l, eqn, ng, go, fsa, fsa2, twovar, row, newrow,
              ig;

    gens := rws.generators;
    go := rws.generatorOrder;
    ng := Length(gens);

    #Make preprocessing file
    rws.operations.ExtVarHandler(rws,filename);
    Read(Concatenation(filename,".gap1"));
    rws.oldNames := _RWS.oldNames;
    for i in [1..ng] do
      if IsBound(_RWS.(i+ng)) then rws.(i+ng):=_RWS.(i+ng); fi;
    od;

    if kb then
      #Read in updated version of equations.
      if not READ(Concatenation(filename,".kbprog")) then
         Error("Could not open output of external Knuth-Bendix program.");
      fi;
      rws.equations := _RWS.equations;
      for eqn in rws.equations do
        eqn[1] := rws.operations.WordToList(eqn[1],gens);
        eqn[2] := rws.operations.WordToList(eqn[2],gens);
      od;
      rws.isConfluent := _RWS.isConfluent;
    fi;

    # read reduction machine - unfortunately, the generators
    # might be in the wrong order, so we have to permute the transition table.
    if kb then
      if not READ(Concatenation(filename,".reduce")) then
         Error("Could not open reduction machine file");
      fi;
      rws.reductionFSA := _RWS.reductionFSA;
      fsa := rws.reductionFSA;
    else
      if not READ(Concatenation(filename,".wa")) then
         Error("Could not open word-acceptor file");
      fi;
      rws.wa := _RWS.wa;
      if not READ(Concatenation(filename,".diff1c")) then
         Error("Could not open word-difference file");
      fi;
      rws.diff1c := _RWS.diff1c;
      fsa := rws.wa; fsa2 := rws.diff1c;
    fi;

    #Reset any lost existing externals
    if rws.oldNames then
      _RWS := rws;
      Read(Concatenation(filename,".gap2"));
      for i in [1..ng+2] do Unbind(rws.(i+ng)); od;
    fi;

    _RWS := rec();
    
    InitializeFSA(fsa);
    #Make sure the table is stored densely
    DenseDTableFSA(fsa);
    Unbind(fsa.sparseTable);
    Unbind(fsa.alphabet.names); fsa.alphabet.type:="simple";
    #Now permute the edge labels.
    for i in [1.. fsa.states.size] do
       row := fsa.denseDTable[i];
       newrow := [];
       for j in [1..ng] do newrow[go[j]] := row[j]; od;
       fsa.denseDTable[i] := newrow;
    od;
    if not kb then
       InitializeFSA(fsa2);
       #Make sure the table is stored densely
       DenseDTableFSA(fsa2);
       Unbind(fsa2.sparseTable);
       fsa2.alphabet.base.printingStrings:=rws.namesGenerators;
       fsa2.states.printingStrings:=rws.namesGenerators;
       #Now permute the edge labels.
       InitializeSR(fsa2.alphabet);
       for i in [1.. fsa2.states.size] do
         row := fsa2.denseDTable[i];
         newrow := [];
         go[ng+1] := ng+1; #for convenience
         for j in [1..ng+1] do for k in [1..ng+1] do
           if j<=ng or k <= ng then
             newrow[(go[j]-1)*(ng+1)+go[k]] := row[(j-1)*(ng+1)+k];
           fi;
         od; od;
         fsa2.denseDTable[i] := newrow;
       od;
    fi;

    #Finally tidy up an irritating problem in the KB case!
    #If some of the generators have names like "a^-1", then any equations of
    #form [a*a^-1,IdWord] will have had their left-hand-sides reduced to
    #IdWord on reading in.  We correct this now.
    if kb then
       ig := rws.invGenerators;
       for i in [1..ng] do
          j := ig[i];
          if IsInt(j) then 
            k := fsa.denseDTable[1][i];
            if k>0 and fsa.denseDTable[1][j]>0 then
               k := -fsa.denseDTable[k][j];
               #k is the number of the equation which should have LHS [i,j].
               rws.equations[k][1] := [i,j];
            fi;
          fi;
       od;
    fi;

end;

#############################################################################
##
#F  RWSOps.WriteRWS(<rws>, [<filename>], [<endsymbol>])
##           . . . . . . . . . . . .write an rws to a file in external format
##
##  WriteRWS prints the rws <rws> to the  file <filename> formatting nicely.
##  It works by building up the material to be printed line by line as strings,
##  and calling LinePrintRWS to print each individual line.
##  If <filename> is not present, or empty, then writing is to the terminal
##  and is simply of form rec(..).
##  Otherwise, printing takes form _RWS := rec(...)<endsymbol>
##  where <endsymbol> is a string which is ";" by default.
##  (_RWS is a global variable.)
##
##  Public function.
RWSOps.WriteRWS := function ( arg )
    local rws, name, filename, gapfilename, line, i, eqn, endsymbol,
          ng, en, gn, ls, go, ig;

    if Length(arg)<1 then
       Error("WriteRWS has 1, 2 or 3 arguments");
    fi;
    rws := arg[1];
    filename := "";
    if Length(arg)>=2 then filename := arg[2]; fi;
    if filename="" then endsymbol := ""; else endsymbol := ";"; fi;
    if Length(arg)>=3 then endsymbol := arg[3]; fi;
    
    if not IsInternalRWS(rws) then
      Error("First argument is not an internal rewriting system.");
    fi;

    ng := rws.numGenerators;
    en := rws.namesGenerators;

    #Now print main file
    if filename="" then Print("rec(\n");
    else PrintTo(filename,"_RWS := \nrec (\n");
    fi;

    line := String("isRWS",16);
    line := Concatenation(line," := true,");
    rws.operations.LinePrintRWS(line,filename);

    if IsBound(rws.isConfluent) then
      line := String("isConfluent",16);
      line := Concatenation(line," := ",String(rws.isConfluent),",");
      rws.operations.LinePrintRWS(line,filename);
    fi;

#Now come all of the optional parameters
    if IsBound(rws.tidyint) then
      line := String("tidyint",16);
      line := Concatenation(line," := ",String(rws.tidyint),",");
      rws.operations.LinePrintRWS(line,filename);
    fi;
    if IsBound(rws.maxeqns) then
      line := String("maxeqns",16);
      line := Concatenation(line," := ",String(rws.maxeqns),",");
      rws.operations.LinePrintRWS(line,filename);
    fi;
    if IsBound(rws.maxstates) then
      line := String("maxstates",16);
      line := Concatenation(line," := ",String(rws.maxstates),",");
      rws.operations.LinePrintRWS(line,filename);
    fi;
    if IsBound(rws.maxreducelen) then
      line := String("maxreducelen",16);
      line := Concatenation(line," := ",String(rws.maxreducelen),",");
      rws.operations.LinePrintRWS(line,filename);
    fi;
    if IsBound(rws.confnum) then
      line := String("confnum",16);
      line := Concatenation(line," := ",String(rws.confnum),",");
      rws.operations.LinePrintRWS(line,filename);
    fi;
    if IsBound(rws.maxwdiffs) then
      line := String("maxwdiffs",16);
      line := Concatenation(line," := ",String(rws.maxwdiffs),",");
      rws.operations.LinePrintRWS(line,filename);
    fi;
    if IsBound(rws.maxstoredlen) then
      line := String("maxstoredlen",16);
      line := Concatenation(line, " := [",
              String(rws.maxstoredlen[1]),",",String(rws.maxstoredlen[2]),"],");
      rws.operations.LinePrintRWS(line,filename);
    fi;
    if IsBound(rws.sorteqns) then
      line := String("sorteqns",16);
      line := Concatenation(line," := ",String(rws.sorteqns),",");
      rws.operations.LinePrintRWS(line,filename);
    fi;
    if IsBound(rws.maxoplen) then
      line := String("maxoplen",16);
      line := Concatenation(line," := ",String(rws.maxoplen),",");
      rws.operations.LinePrintRWS(line,filename);
    fi;
    if IsBound(rws.silent) then
      line := String("silent",16);
      line := Concatenation(line," := ",String(rws.silent),",");
      rws.operations.LinePrintRWS(line,filename);
    fi;
    if IsBound(rws.verbose) then
      line := String("verbose",16);
      line := Concatenation(line," := ",String(rws.verbose),",");
      rws.operations.LinePrintRWS(line,filename);
    fi;
    if IsBound(rws.veryVerbose) then
      line := String("veryVerbose",16);
      line := Concatenation(line," := ",String(rws.veryVerbose),",");
      rws.operations.LinePrintRWS(line,filename);
    fi;

    go := rws.generatorOrder;
    line := Concatenation(String("generatorOrder",16)," := [");
    for i in [1..ng] do
      if i > 1 then
        line := Concatenation(line,",");
      fi;
      if i=1 or Length(line)+Length(en[go[i]]) <= 76 then
        line := Concatenation(line,en[go[i]]);
      else
        rws.operations.LinePrintRWS(line,filename);
        line := String("",21);
        line := Concatenation(line,en[go[i]]);
      fi;
    od;
    line := Concatenation(line,"],");
    rws.operations.LinePrintRWS(line,filename);

    ig := rws.invGenerators;
    line := Concatenation(String("inverses",16)," := [");
    for i in [1..ng] do
      if i > 1 then line := Concatenation(line,","); fi;
      if IsInt(ig[go[i]]) then
        if i=1 or Length(line)+Length(en[ig[go[i]]]) <= 76 then
          line := Concatenation(line,en[ig[go[i]]]);
        else
          rws.operations.LinePrintRWS(line,filename);
          line := String("",21);
          line := Concatenation(line,en[ig[go[i]]]);
        fi;
      fi;
    od;
    line := Concatenation(line,"],");
    rws.operations.LinePrintRWS(line,filename);

    if not IsString(rws.ordering) then
       Error("Can only output orderings that are strings");
    fi;
    line := String("ordering",16);
    line := Concatenation(line," := \"",rws.ordering,"\",");
    rws.operations.LinePrintRWS(line,filename);

    if rws.ordering="wtlex" and IsBound(rws.weight) then
      line := Concatenation(String("weight",16)," := [");
      for i in [1..ng] do
        if i > 1 then
          line := Concatenation(line,",");
        fi;
        line := Concatenation(line,String(rws.weight[go[i]]));
      od;
      line := Concatenation(line,"],");
      rws.operations.LinePrintRWS(line,filename);
    fi;

    if rws.ordering="wreathprod" and IsBound(rws.level) then
      line := Concatenation(String("level",16)," := [");
      for i in [1..ng] do
        if i > 1 then
          line := Concatenation(line,",");
        fi;
        line := Concatenation(line,String(rws.level[go[i]]));
      od;
      line := Concatenation(line,"],");
      rws.operations.LinePrintRWS(line,filename);
    fi;

    line := Concatenation(String("equations",16)," := [");
    for i in [1..Length(rws.equations)] do
      if i > 1 then line := Concatenation(line,","); fi;
      rws.operations.LinePrintRWS(line,filename);
      eqn := rws.equations[i];
      line := Concatenation(String("[",10),
                              rws.operations.ListToWordString(eqn[1],en),",");
      if Length(line)>=40 then
        rws.operations.LinePrintRWS(line,filename);
        line := String("",10);
      fi;
      line :=Concatenation(line,rws.operations.ListToWordString(eqn[2],en),"]");
    od;
    rws.operations.LinePrintRWS(line,filename);
    line := String("]",8);
    rws.operations.LinePrintRWS(line,filename);
    line := Concatenation(")",endsymbol);
    rws.operations.LinePrintRWS(line,filename);
end;

#############################################################################
##
#F  RWSOps.Print(<rws>) . . . print the rws "rws"
##
##  The display function for an rws.
##  It merely calls the WriteRWS function
RWSOps.Print := function ( rws )
    RWSOps.WriteRWS(rws);
end;

#############################################################################
##
#F  IsReducedWordRWS(<rws>,<w>) . . . . tests if a word is reduced
##
##  IsReducedWordRWS tests whether the word <w>
##  is reduced, using the  rewriting system <rws>.
##  <w> can be given either as a list of integers (internal format) or as
##  a word in the generators of the underlying group or monoid.
##  Either the word-acceptor (automatic group case) or the reduction FSA
##  must be defined.
##  It merely calls the corresponding FSA function.
##
##  Public function.
IsReducedWordRWS := function ( rws, w )
    local i, ng, gens;
    if not IsInternalRWS(rws)  then
       Error("First argument is not an internal rewriting system.");
    fi;
    if not IsList(w) and not IsWord(w) then
       Error("Second argument is not a word or list.");
    fi;
    ng := rws.numGenerators;
    gens := rws.generators;
    if IsWord(w) then w:=rws.operations.WordToList(w,gens); fi;
    for i in w do
      if not i in [1..ng] then Error("Invalid entry in word or list");fi;
    od;
    if IsBound(rws.wa) then
    # Automatic group case - use word-acceptor
      return IsAcceptedWordDFA( rws.wa,w );
    fi;
    if not IsBound(rws.reductionFSA)  then
       Error("First argument does not have initialized dfa as field.");
    fi;
    return IsAcceptedWordDFA( rws.reductionFSA,w );
end;
      
#############################################################################
##
#F  ReduceWordRWS(<rws>,<w>) . . . . reduces a word using rewriting-system
##
##  ReduceWordRWS reduces the word <w>, using the rewriting-system <rws>.
##  <w> can be given either as a list of integers (internal format) or as
##  a word in the generators of the underlying group or monoid.
##  Either the reduction FSA, or one of the word-difference automata (in the
##  automatic group case) must be defined for <rws>.
##  In the latter case, the separate function ReducedWordWD is called.
##
##  Public function.
ReduceWordRWS := function ( rws, w )
    local fsa, pos, label, state, history, eqn, sublen, table, ng, go, i,
	  gens, word;
    if not IsInternalRWS(rws)  then
       Error("First argument is not an internal rewriting system.");
    fi;
    if not IsAvailableReductionRWS(rws) then
       Error(
        "Reduction algorithm unavailable. You must run KB or Automata first.");
    fi;
    if not IsList(w) and not IsWord(w) then
       Error("Second argument is not a word or list.");
    fi;
    ng := rws.numGenerators;
    gens := rws.generators;
    if IsWord(w) then word :=true; w:=rws.operations.WordToList(w,gens);
    else word := false;
    fi;
    for i in w do
      if not i in [1..ng] then Error("Invalid entry in word or list");fi;
    od;
    for i in w do if not i in [1..ng] then Error("Invalid entry in list");fi;od;
    if IsBound(rws.warningOn) and rws.warningOn then
      Print(
   "#WARNING: system is not confluent, so reductions may not be to normal form.\n");
      rws.warningOn:=false;
           # only give the warning once, or it could become irritating!
    fi;
    go := rws.generatorOrder;
    if IsBound(rws.diff2) then
     # automatic group case
       w := rws.operations.ReducedWordWD(rws.diff2,w,go);
    elif IsBound(rws.diff1c) then
     # automatic group case
       w := rws.operations.ReducedWordWD(rws.diff1c,w,go);
    elif IsBound(rws.diff1) then
     # automatic group case
       w := rws.operations.ReducedWordWD(rws.diff1,w,go);
    elif IsBound(rws.reductionFSA)  then
       fsa := rws.reductionFSA;
       if not IsInitializedFSA(fsa) or IsDeterministicFSA(fsa)=false  then
          Error("First argument does not have initialized dfa as field.");
       fi;
   
       state := fsa.initial[1];
       pos := 1;
       history:= [];
       history[1] := state; # history[i] = state before reading i-th symbol.
       table := DenseDTableFSA(fsa);
       while pos <= Length(w) do
         state := table[state][w[pos]];
         if state>0 then
           pos := pos+1;
           history[pos] := state;
         else
           state := -state;
           eqn := rws.equations[state];
           sublen := Length(eqn[1]);
           fsa.operations.SubstitutedList(w,pos-sublen+1,pos,eqn[2]);
           pos := pos-sublen+1;
           state := history[pos];
         fi;
       od;
    else
       Error("First argument does not have initialized dfa as field.");
    fi;

    if word then w := rws.operations.ListToWord(w,gens); fi;
    return w;
end;
      
#############################################################################
##
#F  RWSOps.ReducedWordWD(<wd>,<w>,<go>)
##                   . . . . .reduces a word using word-difference automaton
##
##  ReducedWordWD calculates the reduction of the word <w> (list of integers)
##  using the word-difference automaton <wd>.
##  <wd> should be two-variable, where <w> is a list of integers in the range
##  [1..ng], where ng is the size of the base alphabet.
##  go is the order of generators in the underlying rewriting system.
##  WARNING: No validity checks are carried out!
##
##  Private function.
RWSOps.ReducedWordWD := function ( wd, w, go)
    local  ndiff, ngens, ng1, identity, level, cf, gpref, gct, gptr,
           diff, newdiff, deqi, gen1, gen2, donesub, donediffs, subvert,dosub,
           g2ltg1, diffct, t, nlen, olen, i, l, table, gpos;
    if not IsInitializedFSA(wd)  or IsDeterministicFSA(wd)=false  then
       Error("First argument is not an initialized dfa.");
    fi;

    ndiff := wd.states.size;
    ngens := wd.alphabet.base.size;
    #The inverse of the generator order list is more useful - so work it out!
    gpos:=[]; for i in [1..ngens] do gpos[go[i]] := i; od;
    ng1 := ngens+1;            
    identity := wd.initial[1];
    if Length(w) <= 0 then
      return w;
    fi;
    cf := [];
    # cf is used as a characteristic function, when constructing a subset of the
    # set  D  of word differences.
    gpref := []; gct := 0; gpref[1] := 0; gptr := [];
    # gpref[n]  is the number of "vertices" that have been defined after
    # reading the first n-1 elements of the word.
    # These vertices are gptr[1],...,gptr[gpref[n]].
    # A vertex is a record with 4 components, backptr, genno, diffno, sublen,
    # It represents a vertex in the graph of strings that may eventually
    # be used as substituted strings in the word w.
    # backptr is either undefined or another vertex.
    # gen is the generator at the vertex.
    # diffno is the word-difference number of the string at which the vertex
    #        is at the end - this string is reconstructed using backptr.
    # sublen is plus or minus the length of this string. sublen is positive
    #        iff the string lexicographically precedes the corresponding
    #        prefix of the word being reduced.

    # Now we start scanning the word.
    table := DenseDTableFSA(wd);
    level := 1;
    while level <= Length(w) do
      for i in [1..ndiff] do cf[i] := false; od;
      gen1 := w[level];
      # The next loop is over the identity, and the subset of the set of
      # word-differences (states of wd) defined at the previous level (level-1)

      diff := identity;
      donesub := false;
      donediffs := false;
      while not donesub and not donediffs do
        deqi := diff = identity;
      # First look for a possible substitution of a shorter string
        newdiff := table[diff][ng1*gen1];
        if newdiff=identity then
          #Make substitution  reducing length of word by 1
          wd.operations.SubstitutedList(w,level,level,[]);
          i := level-1;
          if not deqi then
            subvert := gptr[diffct];
	    dosub := true;
            while dosub do
              w[i] := subvert.gen;
              i := i-1;
              if IsBound(subvert.backptr) then
	        subvert := subvert.backptr;
              else
                dosub := false;
              fi;
            od;
          fi;
          #Whenever we make a substitution, we have to go back one level more
          #than expected, because of our policy of looking ahead for
          #substitutions that reduce the length by 2.
          if i>0 then level:=i-1; else level:=i; fi;
          gct := gpref[level+1];
          donesub := true;
        elif newdiff>0 and level<Length(w) then
          #See if there is a substitution reducing length by 2.
          gen2 := w[level+1];
          t := table[newdiff][ng1*gen2];
          if t=identity then
            #Make substitution  reducing length of word by 2
            wd.operations.SubstitutedList(w,level,level+1,[]);
            i := level-1;
            if not deqi then
              subvert := gptr[diffct];
  	      dosub := true;
              while dosub do
                w[i] := subvert.gen;
                i := i-1;
                if IsBound(subvert.backptr) then
  	          subvert := subvert.backptr;
                else
                  dosub := false;
                fi;
              od;
            fi;
            if i>0 then level:=i-1; else level:=i; fi;
            gct := gpref[level+1];
            donesub := true;
          fi;
        fi;

        if not donesub then
          #Now we loop over the generator that is a candidate for
          #substitution at this point.
          for gen2 in [1..ngens] do
            g2ltg1 := gpos[gen2] < gpos[gen1];
            newdiff := table[diff][ng1*(gen1-1)+gen2];
            if donesub then
              donesub := true;
              #i.e. do nothing - we really want to break from the for loop here.
            elif newdiff=identity then
              if deqi then #only occurs when gen2 and gen1 are equal in group
                if g2ltg1 then
                  w[level] := gen2;
                  if level>1 then level:=level-2; else level:=level-1; fi;
                  gct := gpref[level+1];
                  donesub := true;
                fi;
              elif gptr[diffct].sublen>0 then
                #Make a substitution by a string of equal length.
                w[level] := gen2;
                i := level-1;
                subvert := gptr[diffct];
    	        dosub := true;
                while dosub do
                  w[i] := subvert.gen;
                  i := i-1;
                  if IsBound(subvert.backptr) then
    	            subvert := subvert.backptr;
                  else
                    dosub := false;
                  fi;
                od;
                if i>0 then level:=i-1; else level:=i; fi;
                gct := gpref[level+1];
                donesub := true;
              fi;
            elif newdiff>0 then
              if cf[newdiff] then
                #We have this word difference stored already, but we will check
                #to see if the current string precedes the existing one.
                i := gpref[level];
                repeat
                  i := i+1;
                  subvert := gptr[i];
                until subvert.diffno=newdiff;
                olen := subvert.sublen;
                if deqi then 
                  if g2ltg1 then nlen:=1; else nlen:= -1; fi;
                else
                  l := gptr[diffct].sublen;
                  if l>0 then nlen:=l+1; else nlen:=l-1; fi;
                fi;
                if nlen > olen then # new string is better than existing one
                  subvert.gen := gen2;
                  subvert.sublen := nlen;
                  if deqi then Unbind(subvert.backptr);
                  else subvert.backptr := gptr[diffct];
                  fi;
                fi;
              else
               # this is a new word-difference at this level, so
               # we define a new vertex.
                gct := gct+1;
                gptr[gct] := rec();
                if deqi then 
                  if g2ltg1 then nlen:=1; else nlen:= -1; fi;
                else
                  l := gptr[diffct].sublen;
                  if l>0 then nlen:=l+1; else nlen:=l-1; fi;
                fi;
                subvert := gptr[gct];
                subvert.gen := gen2;
                subvert.diffno := newdiff;
                subvert.sublen := nlen;
                if not deqi then subvert.backptr := gptr[diffct]; fi;
                cf[newdiff] := true;
              fi;
            fi;
          od; # End of loop over gen2

          if not donesub then
            #Go on to next word-difference from the previous level
            if diff=identity then
              if level=1 then
                donediffs := true;
              else
                diffct := gpref[level-1] + 1;
              fi;
            else
              diffct := diffct+1;
            fi;
            if not donesub and not donediffs then
              if diffct > gpref[level] then
                donediffs := true;
              else
                diff := gptr[diffct].diffno;
              fi;
            fi;
          fi;
        fi;
      od; #end of loop over word-differences at previous level

      level := level+1;
      gpref[level] := gct;
    od;
    return w;
end;

#############################################################################
##
#F  SizeRWS(<rws>>) . . . . . number of reduced words in a rewriting system
##
##  This merely calls the corresponding FSA function.
##
##  Public function.
SizeRWS := function ( rws )
    if not IsInternalRWS(rws)  then
       Error("First argument is not a rewriting system.");
    fi;
    if not IsAvailableReductionRWS(rws) then
       Error(
        "Size algorithm unavailable. You must run KB or Automata first.");
    fi;
    if IsBound(rws.wa) then
     # automatic group case
      return LSizeDFA( rws.wa );
    fi;
    return LSizeDFA( rws.reductionFSA );
end;

#############################################################################
##
#F  EnumerateRWS(<rws>, <min>, <max>) . . . enumerate reduced words in a rws
##
##  This merely calls the corresponding FSA function.
##  Words are converted to words in generators of underlying group or monoid
##  before returning.
##
##  Public function.
EnumerateRWS := function ( rws, min, max )
    local ret, x, i;
    if not IsInternalRWS(rws)  then
       Error("First argument is not a rewriting system.");
    fi;
    if not IsAvailableReductionRWS(rws) then
       Error(
        "Enumeration algorithm unavailable. You must run KB or Automata first.");
    fi;
    if IsBound(rws.wa) then
     # automatic group case
      ret := LEnumerateDFA( rws.wa,min,max );
    else
      ret := LEnumerateDFA( rws.reductionFSA,min,max );
    fi;
    for i in [1..Length(ret)] do
      x := ret[i]; ret[i] := rws.operations.ListToWord(x,rws.generators);
    od;
    return ret;
end;

#############################################################################
##
#F  SortEnumerateRWS(<rws>, <min>, <max>)  . . enumerate reduced words and sort
##
##  This merely calls the corresponding FSA function.
##  Words are converted to words in generators of underlying group or monoid
##  before returning.
##
##  Public function.
SortEnumerateRWS := function ( rws, min, max )
    local ret, x, i;
    if not IsInternalRWS(rws)  then
       Error("First argument is not a rewriting system.");
    fi;
    if not IsAvailableReductionRWS(rws) then
       Error(
        "Enumeration algorithm unavailable. You must run KB or Automata first.");
    fi;
    if IsBound(rws.wa) then
     # automatic group case
      ret := SortLEnumerateDFA( rws.wa,min,max );
    else
      ret := SortLEnumerateDFA( rws.reductionFSA,min,max );
    fi;
    for i in [1..Length(ret)] do
      x := ret[i]; ret[i] := rws.operations.ListToWord(x,rws.generators);
    od;
    return ret;
end;

#############################################################################
##
#F  SizeEnumerateRWS(<rws>, <min>, <max>)  . . . . number of reduced words 
##
##  This merely calls the corresponding FSA function.
##
##  Public function.
SizeEnumerateRWS := function ( rws, min, max )
    if not IsInternalRWS(rws)  then
       Error("First argument is not a rewriting system.");
    fi;
    if not IsAvailableReductionRWS(rws) then
       Error(
        "Enumeration algorithm unavailable. You must run KB or Automata first.");
    fi;
    if IsBound(rws.wa) then
     # automatic group case
      return SizeLEnumerateDFA( rws.wa,min,max );
    fi;
    return SizeLEnumerateDFA( rws.reductionFSA,min,max );
end;

#############################################################################
##
#F  KB(<rws>)  . . . . call external Knuth-Bendix program on rws.
##
##  This returns true if a confluent rewriting system results - otherwise
##  false. In the latter case, words can still be rewritten with respect to
##  the current equations, but they are not guaranteed to reduce to the unique
##  representative of the group element.
##  An error message results if the external program aborts without outputting.
##  Public function.
KB := function ( rws )
    if not IsInternalRWS(rws)  then
       Error("First argument is not a rewriting system.");
    fi;
    if IsConfluentRWS(rws) then
       Print("#The rewriting system is already confluent.\n");
       Print("#Call - ResetRWS(<rws>) to restart.\n");
       return;
    fi;
    rws.operations.WriteRWS(rws,_KBTmpFileName);
    #Keep the original equations, in case we want them again.
    if not IsBound(rws.originalEquations) then
      rws.originalEquations := Copy(rws.equations);
    fi;
    InfoRWS("#I Calling external Knuth-Bendix program.\n");
    Exec(Concatenation(_KBExtDir,"kbprog ",_KBTmpFileName));
    rws.operations.UpdateRWS(rws,_KBTmpFileName,true);
    Exec(Concatenation("/bin/rm -f ",_KBTmpFileName,"*"));
    InfoRWS("#I External Knuth-Bendix program complete.\n");
    if rws.isConfluent then
      InfoRWS("#I System computed is confluent.\n");
      rws.isAvailableNormalForm := true;
      rws.isAvailableReduction := true;
      rws.warningOn := false;
    else
      InfoRWS("#I System computed is NOT confluent.\n");
      rws.isAvailableNormalForm := false;
      rws.isAvailableReduction := true;
      rws.warningOn := true;
    fi;
    return rws.isConfluent;
end;

#############################################################################
##
#F  ResetRWS(<rws>)  . . . . . . . . . . . reset rws after a call of KB.
##
##  Public function.
##  This resets a rewriting system back to the original equations, after a
##  call of KB or Automata.
##  Perhaps useful if order of generators is to be changed.
ResetRWS := function ( rws )
    if not IsInternalRWS(rws)  then
       Error("First argument is not a rewriting system.");
    fi;
    Unbind(rws.isConfluent);
    Unbind(rws.isAvailableNormalForm);
    Unbind(rws.isAvailableReduction);
    Unbind(rws.warningOn);
    Unbind(rws.reductionFSA);
    Unbind(rws.wa);
    Unbind(rws.diff1c);
    if IsBound(rws.originalEquations) then
       Unbind(rws.equations);
       rws.equations := rws.originalEquations;
       Unbind(rws.originalEquations);
    fi;
end;

#############################################################################
##
#F  AddOriginalEqnsRWS(<rws>). add original equations to rws after a call of KB.
##
##  This appends the original equations to the list of equations, after a
##  call of KB. Useful for a re-check, if the external program may have
##  deleted some equations.
##  After this function, rewrtiting is no longer possible.
##  Public function.
AddOriginalEqnsRWS := function ( rws )
    if not IsInternalRWS(rws)  then
       Error("First argument is not a rewriting system.");
    fi;
    Unbind(rws.reductionFSA);
    Unbind(rws.isConfluent);
    if IsBound(rws.originalEquations) then
      Append(rws.equations,rws.originalEquations);
    fi;
end;

#############################################################################
##
#F  Automata(<rws>, [<large>], [<filestore>], [<diff1>]) 
##                      . . . . call external automatic group program on rws.
##
##  This returns true if the automatic group programs succeed - otherwise
##  false.
##  The optional parameters are all boolean, and false by default.
##  <large> means problem is large - the external programs use bigger tables.
##  <filestore> means external programs use less core memory and more external
##         files - they run a little slower.
##  <diff1> is necessary on some examples - see manual for information.
##  Public function.
Automata := function ( arg )
    local  narg, rws, large, filestore, diff1, callstring, optstring;
    narg := Number(arg);
    if narg<1  or  not IsInternalRWS(arg[1]) then
       Error("First argument is not a rewriting system.");
    fi;
    rws := arg[1];
    if not rws.ordering = "shortlex" then
       Error("Ordering must be shortlex for automatic group programs");
    fi;
    if IsConfluentRWS(rws) then
       Print("#The rewriting system is already confluent.\n");
       Print("#Call - ResetRWS(<rws>) to restart.\n");
       return;
    fi;
    large:=false; filestore:=false; diff1:=false;
    if narg>=2 and arg[2]=true then large:=true; fi;
    if narg>=3 and arg[3]=true then filestore:=true; fi;
    if narg>=4 and arg[4]=true then diff1:=true; fi;
    rws.operations.WriteRWS(rws,_KBTmpFileName);
    callstring := Concatenation(_KBExtDir,"autgroup ");
    optstring := "";
    if large then Append(optstring," -l "); fi;
    if filestore then Append(optstring," -f "); fi;
    if diff1 then Append(optstring," -d "); fi;
    if IsBound(rws.silent) and rws.silent then Append(optstring," -s "); fi;
    if IsBound(rws.verbose) and rws.verbose then Append(optstring," -v "); fi;
    if IsBound(rws.veryVerbose) and rws.veryVerbose then
                      Append(optstring," -vv "); fi;
    callstring := Concatenation(callstring,optstring,_KBTmpFileName);
    InfoRWS("#I Calling external automatic groups program.\n");
    Exec(callstring);
    callstring := Concatenation(_KBExtDir,"gpminkb ");
    optstring := "";
    if IsBound(rws.silent) and rws.silent then
                      optstring := Concatenation(optstring," -s "); fi;
    if IsBound(rws.verbose) and rws.verbose then
                      optstring := Concatenation(optstring," -v "); fi;
    if IsBound(rws.veryVerbose) and rws.veryVerbose then
                      optstring := Concatenation(optstring," -vv "); fi;
    callstring := Concatenation(callstring,optstring,_KBTmpFileName);
    if READ(Concatenation(_KBTmpFileName,".success")) then
     InfoRWS("#I Computation was successful - automatic structure computed.\n");
      Exec(callstring);
      rws.operations.UpdateRWS(rws,_KBTmpFileName,false);
      Exec(Concatenation("/bin/rm -f ",_KBTmpFileName,"*"));
      rws.isAvailableNormalForm := true;
      rws.isAvailableReduction := true;
      rws.warningOn := false;
      return true;
    else
      Exec(Concatenation("/bin/rm -f ",_KBTmpFileName,"*"));
      InfoRWS("#I Computation was not successful.\n");
      return false;
    fi;
end;

#############################################################################
##
#F  RWSOps.WordToList( <w>, <genlist> ) .  converts word to list of integers
##
## <w> should be a word in the generators contained in the list <genlist>.
## The corresponding list of integers is computed an retruned.
RWSOps.WordToList := function ( w, genlist )
   local i, l;
   l := [];
   for i in [1..LengthWord(w)] do
       Add(l,Position(genlist,Subword(w,i,i)));
   od;

   return l;
end;

#############################################################################
##
#F  RWSOps.ListToWord( <l>, <genlist> ) .  converts list of integers to word
##
## <l> should be a list of integers in range [1..<ng>], where <genlist>
## is a list of generators of length <ng>.
## The corresponding word in the generators in <genlist> is computed.
RWSOps.ListToWord := function ( l, genlist )
   local i, w;
   w := IdWord;
   for i in [1..Length(l)] do
       w := w*genlist[l[i]];
   od;

   return w;
end;

#############################################################################
##
#F  RWSOps.ListToWordString( <list>, <names> )
##                                     . . . . . converts <list> to a string
##
## <list> is a list of integers in the range [1..ng] for some ng>=0.
## <names> is a list of strings of length ng.
## The list is converted into a string representing the word for printing.
##  Private function.
RWSOps.ListToWordString := function ( list, names )
   local string, i, l, ls, ng, g, ct, lg;
   l := Length(list);
   if l=0 then
      return "IdWord";
   fi;
   string := "";
   lg := 0; ct := 0;
   for i in [1..l] do
     g := list[i];
     if g <> lg and lg <> 0 then
        if string<>"" then string := Concatenation(string,"*"); fi;
        ng := names[lg];
        if ct > 1 then
           #Check to see if names[lg] ends "^-1"
           ls := Length(ng);
           if ls>3 and ng[ls-2]='^' and ng[ls-1]='-' and ng[ls]='1' then
             string := Concatenation(string,ng{[1..ls-1]},String(ct));
           else
             string := Concatenation(string,ng,"^",String(ct));
           fi;
        else string := Concatenation(string,ng);
        fi;
        ct := 0;
     fi;
     lg := g; ct := ct+1;
   od;
   if string<>"" then string := Concatenation(string,"*"); fi;
   ng := names[lg];
   if ct > 1 then
       #Check to see if names[lg] ends "^-1"
       ls := Length(ng);
       if ls>3 and ng[ls-2]='^' and ng[ls-1]='-' and ng[ls]='1' then
         string := Concatenation(string,ng{[1..ls-1]},String(ct));
       else
         string := Concatenation(string,ng,"^",String(ct));
       fi;
   else string := Concatenation(string,ng);
   fi;
 
   return string;
end;

WAppendTo:=function(arg)local s; # "wide" AppendTo
  s:=SizeScreen();SizeScreen([250,70]);
  ApplyFunc(AppendTo,arg);
  SizeScreen(s);
end;

#############################################################################
##
#F  RWSOps.LinePrintRWS(<line> [,<filename>]) . . . . . . . print the line x
##
##  LinePrintRWS prints the line (a string) to the terminal (default)
##  or to file filename  if specified, formatting nicely.
##  It works by building up the material to be printed line by line as strings,
##  and calling LinePrintRWS to print each individual line.
RWSOps.LinePrintRWS := function ( arg )
    local line, filename;
    line := arg[1];
    if Length(arg) = 1 then
      filename := "";
    else
      filename := arg[2];
    fi;
    if filename = "" then
      Print(line,"\n");
    else
      WAppendTo(filename,line,"\n");
    fi;
end;

#############################################################################
##
#F  RWSOps.ExtVarHandler(<rws>, <filename>) . . write file to handle externals
##
## This is hopefully a temporary hack, for use until GAP V4, where should be
## a better solution. A GAP file is written to preserve any existing values
## of external variables corresponding to the generator names of the
## rewriting system <rws>, and then to declare these variables to be equal
## to the corresponding generators of <rws>. This is necessary for reading
## back in the output of KB or Automata, which uses these names.
## This first file is called <rws>.gap1.
## A second file <rws>.gap2 is written for reading afterwards, which restores
## the values of any previously existing externals with those names.
##
RWSOps.ExtVarHandler  := function(rws, filename)
    local file1, file2, ng, names, line, i, ni, l;
    file1 := Concatenation(filename,".gap1");
    file2 := Concatenation(filename,".gap2");
    _RWS := rws;
    PrintTo(file1,"_RWS.oldNames:=false;\n");
    PrintTo(file2,"");

    ng := rws.numGenerators;
    names := rws.namesGenerators;
    for i in [1..ng+1] do Unbind(_RWS.(i+ng)); od;

    if Length(names[1])>2 and
       (names[1]{[1,2]}="f." or names[1]{[1,2]}="F.") then
     # This is the default printing string for generators of fp-groups.
     # We need to define a record variable "f" (and possibly "F"),
     # or we will not be able to # read back in the output of the external
     # program.
        line := Concatenation("if IsBound(f) then _RWS.",String(2*ng+1),":=f;",
                "_RWS.oldNames:=true; fi;\n");
        line := Concatenation(line,"f:=rec();\n");
        WAppendTo(file1,line);
        line := Concatenation("if IsBound(_RWS.",String(2*ng+1),") then ",
                "f:=_RWS.",String(2*ng+1),"; fi;\n");
        WAppendTo(file2,line);
        line := Concatenation("if IsBound(F) then _RWS.",String(2*ng+2),":=F;",
                "_RWS.oldNames:=true; fi;\n");
        line := Concatenation(line,"F:=rec();\n");
        WAppendTo(file1,line);
        line := Concatenation("if IsBound(_RWS.",String(2*ng+2),") then ",
                "F:=_RWS.",String(2*ng+2),"; fi;\n");
        WAppendTo(file2,line);
    fi;
    for i in [1..ng] do
      ni := names[i]; l := Length(ni);
      if l <= 3 or ni{[l-2..l]} <> "^-1" then
        line := Concatenation("if IsBound(",ni,") and ",ni,
                " <> _RWS.",String(i)," then _RWS.",String(i+ng),":=",
                ni,"; _RWS.oldNames:=true; fi;\n");
        line := Concatenation(line,ni,":=_RWS.",String(i),";\n");
        WAppendTo(file1,line);
        line := Concatenation("if IsBound(_RWS.",String(i+ng),") then ",
                ni,":=_RWS.",String(i+ng),"; fi;\n");
        WAppendTo(file2,line);
      fi;
    od;
    WAppendTo(file1,"_:=IdWord;\n");
end;
