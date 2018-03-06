#H##########################################################################
##A  DCE.G      Steve Linton    1993, 1994, 1995 
##
##
##  GAP double coset enumerator.
##
## Documentation of the data structures, invariant assumptions
## and the like is in a separate file, for sanity
##
#############################################################################
##
#F IsDCEUniverse(u)   --  test whether an object is a DCE universe
##

IsDCEUniverse := function(u)
    return IsRec(u) and IsBound(u.isDCEUniverse) and u.isDCEUniverse;
end;

#############################################################################
##
#F DCEInfoPrint(u[,stuff]*) - print stuff tidily
#F InfoDCE<n> - should be DCEInfoPrint or Ignore
## Various functions and definitions used in testing and logging
##   all pretty self explanatory


OldDCEInfoPrint := function(arg)
    local a,u;
    u := arg[1];
    if not IsDCEUniverse(u) then
        Error("Info printing bad first arg");
    fi;
    Print("#I ",String(Runtime()-u.starttime,-12));
    for a in arg{[2..Length(arg)]} do
        Print(a," ");
    od;
    Print("\n");
end;

DCEInfoPrint := function(arg)
    local a,u;
    u := arg[1];
    if not IsDCEUniverse(u) then
        Error("Info printing bad first arg");
    fi;
    Print("#I ");
    for a in arg{[2..Length(arg)]} do
        Print(a," ");
    od;
    Print("\n");
end;


if not IsBound(InfoDCE1) then
    InfoDCE1 := DCEInfoPrint;
#   InfoDCE1 := Ignore;
fi;

if not IsBound(InfoDCE2) then
#    InfoDCE2 := DCEInfoPrint;
   InfoDCE2 := Ignore;
fi;

if not IsBound(InfoDCE3) then
    InfoDCE3 := Ignore;
#   InfoDCE3 := DCEInfoPrint;
fi;

if not IsBound(InfoDCE4) then
    InfoDCE4 := Ignore;
#   InfoDCE4 := DCEInfoPrint;
fi;




#############################################################################
##
#F RightCosetReps(G,U)  Streamlined permgroups only functions produce
#F LeftCosetReps(G,U)     random set of reps
##
## RightTransversal makes this largely unnecessary but never mind.
## we still want to avoid unnecessary checking, but use cached values
## 
##

DCERightCosetReps := function(G,U)
    if IsParent(G) then 
        if not IsBound(U.rightTransversal) then
            if IsBound(U.rightCosets) then
                U.rightTransversal :=  List(U.rightCosets,c->c.representative);
            else
                U.rightTransversal := G.operations.RightTransversal(G,U);
            fi;
        fi;
        return U.rightTransversal;
    fi;    
    return G.operations.RightTransversal(G,U);
end;

DCELeftCosetReps := function(G,U)
    if IsParent(G) then
        if not IsBound(U.leftTransversal) then
            if IsBound(U.leftCosets) then
                U.leftTransversal := List(U.leftCosets,c->c.representative);
            else
                U.leftTransversal := List(DCERightCosetReps(G,U),c->c^-1);
            fi;
        fi;
        return U.leftTransversal;
    fi;
    return List(G.operations.RightTransversal(G,U),c->c^-1);
end;



#############################################################################
##
#F DoubleCosetReps(g,h,k) Does not check arguments and checks special cases
##
##
 
DCEDoubleCosetReps := function(g,h,k)
     local sg;
    sg := Size(g);
    if sg = Size(h) or sg = Size(k) then
        return [g.identity];
    fi;
    if Length(h.generators) = 0 then
        if Length(k.generators) = 0 then
            return Elements(g);
        else
            return DCELeftCosetReps(g,k);
        fi;
    else
        if Length(k.generators) = 0 then
            return DCERightCosetReps(g,h);
        else
            return List(g.operations.DoubleCosets(g,h,k),d->d.representative);
        fi;
    fi;
end;

#############################################################################
##
#F Faster Subgroup, avoids checking
##

DCESubgroup :=  function(G,gens)   
    return  G.operations.Subgroup( G, gens );
end;


#############################################################################
##
#F DCEPermString(p)   Avoid printing very large permutations
##

DCEPermString := function(p)
    if IsPerm(p) and p <> () and LargestMovedPointPerm(p) > 15 then
        return ConcatenationString("<",String(OrderPerm(p)),">");
    else
        return String(p);
    fi;
end;


#############################################################################
##
#F DCEDKToString Convert a dk pair as used in various places to a string
##  for printing. Can't be done with an operations record because the dk pair
##   is a list not a record
##


DCEDKToString := function(dk)
    return ConcatenationString(String(dk[1].name),".",DCEPermString(dk[2]));
end;



############################################################################
##
#V  DCERowOps  - operations record for a row of the DCT
##

DCERowOps := OperationsRecord("DCERowOps");

#############################################################################
##
##   Some Operations for DCERows
##

DCERowOps.Print := function(d)
    local col,x,cohort;
    Print("<< ");
    Print("Row ", String(d.name));
    if d.deleted then
        Print(" Deleted, replaced by ",
              DCEDKToString(d.replace));
    else
        if not IsPermGroup(d.muddle) or
           Size(d.muddle.operations.MovedPoints(d.muddle)) <= 15  then
            Print(" Muddle Group ",d.muddle);
        else 
            Print(" Muddle Group Order ",Size(d.muddle));
        fi;
        for x in [1..Length(d.cohorts)] do
            Print(" ",x,": ");
            cohort := d.cohorts[x];
            for col in [1..Length(cohort)] do
                if IsBound(cohort[col]) then
                    if IsList(cohort[col]) then
                        Print(" ",col," ",DCEDKToString(cohort[col]));
                    else
                        Print(" ",col," -> ",cohort[col].col);
                    fi;
                else
                    Print(" Unknown");
                fi;
            od;
        od;
    fi;
    Print(" >>");
end;


DCEPrintRow := function(row)
    Print(row,"\n");
end;

#############################################################################
##
#F DCEPrintTab(u) PRint the double coset table
##

DCEPrintTab := function(u)
    local d;
    for d in u.table do
        Print(d,"\n\n");
    od;
end;



#############################################################################
##
#F DCEUniverseOps An operations record for DCE universes. Mostly to
##  override print which gets very tangled otherwise.
##
##

DCEUniverseOps := OperationsRecord("DCEUniverseOps");

DCEUniverseOps.Print := function(u)
    Print("<< Double coset table \"",u.name,"\" ",u.status);
    if u.status in ["closed", "early-closed"] then
        Print(" ",Number(u.table)," double ",u.degree," single");
    fi;
   Print(" >>");
end;
            
    

#############################################################################
##
#F DCEXonEl(u,l,x) l should be in Lx. Then returns l^x
##
    

DCEXonEl := function(u,l,x)
    local act;
    act := u.Xrecs[x].action;
    if act = false then
        return l;
    else
        return l^act;
    fi;
end;
    


#############################################################################
##
#F DCEXonSG(u,h,x)  h should be a subgroup of Lx. Then returns h^x
##


DCEXonSG := function(u,h,x)
    local act;
    act := u.Xrecs[x].action;
    if act = false then
        return h;
    elif act in u.K then
        return h^act;
    else 
        return DCESubgroup(u.K,List(h.generators,l->l^act));
    fi;
end;

#############################################################################
##
#F DCELin(u,m,x,k) Compute M intersect Lx^k
##

DCELin := function(u,m,x,k)
    local xrec,gg,l;
    if Length(m.generators) = 0 then
        return u.triv;
    fi;
    xrec := u.Xrecs[x];
    gg := u.GGrecs[xrec.wgg];
    if Length(gg.dom) = 1 then
        return m;
    fi;
    if k = u.id then
        if IsIdentical(m,u.K) then
            return AsSubgroup(u.K,xrec.L);
        fi;
        l := gg.stab(m,xrec.pt,gg.op); 
    else
        if IsIdentical(m,u.K) then
            return AsSubgroup(u.K,xrec.L)^k;
        fi;
        l := gg.stab(m,gg.op(xrec.pt,k),gg.op); 
    fi;
    if Length(l.generators) = 0 then
        return u.triv;
    else  
        return l;
    fi;
end;

#############################################################################
##
#F DCEColumnRegular(u,k,x)    Which column of the cohort do we look up the
## action of k.x, assuming trivial muddle group.
##

DCEColumnRegular := function(u,k,x)
    local xrec,ggrec,im;
    xrec := u.Xrecs[x];
    ggrec := u.GGrecs[xrec.wgg];
    if Length(ggrec.dom) = 1 then
        return 1;
    fi;
    return PositionSorted(ggrec.dom,ggrec.op(xrec.pt,k^-1));
end;



#############################################################################
##
#F DCEColumn(u,dk,x)       Which column of the x cohort of 
##                          row d to actually look up
##                          dk.x   - may alter dk
##


DCEColumn := function(u,dk,x)
    local col,d,xrec,ggrec,cohort,b;
    d := dk[1];
    # col := DCEColumnRegular(u,dk[2],x); Inline this as 
    xrec := u.Xrecs[x];
    ggrec := u.GGrecs[xrec.wgg];
    if Length(ggrec.dom) = 1 then
        return 1;
    fi;
    col := PositionSorted(ggrec.dom,ggrec.op(xrec.pt,dk[2]^-1));
    #
    # End of in-lining
    #
    cohort := d.cohorts[x];
    if IsBound(cohort[col]) and not IsList(cohort[col]) then
        b := cohort[col];
        dk[2] := b.elt*dk[2];
        return b.col;
    fi;
    return col;
end;


#############################################################################
##
#F DCEFollow( u, dk, r)          Compute the image of d.k under r which
##  is either an element of K or an Integer (element of X). dk is passed
##  as a two-element list which is modified and returned. If the image
##  cannot be computed (the table entry is empty) then dk is not modified
##  and false is returned.
##

DCEFollow := function(u,dk,r)
    local col,d,xrec,ggrec,cohort,op,p2s,i,im,score;
    d := dk[1];
    if IsInt(r) then
        #        col := DCEColumn(u,dk,r); Inline this as 
      # col := DCEColumnRegular(u,dk[2],r); Inline this as 
        xrec := u.Xrecs[r];
        ggrec := u.GGrecs[xrec.wgg];
        cohort := d.cohorts[r];
        #
        # First deal with the super-easy case
        #
        if Length(ggrec.dom) = 1 then
            if not IsBound(cohort[1]) then
                return false;
            fi;
            dk[1] := cohort[1][1];
            if xrec.action = false then
                dk[2] := cohort[1][2]*dk[2];
            else
                dk[2] := cohort[1][2]*DCEXonEl(u,dk[2],r);
            fi;
            return dk;
        fi;
        #
        # Otherwise we need to do some work
        #
        op := ggrec.op;
        col := PositionSorted(ggrec.dom,ggrec.op(xrec.pt,dk[2]^-1));
        #
        # End of in-lining
        #
        if IsBound(cohort[col]) and IsRec(cohort[col]) then
            dk[2] := cohort[col].elt*dk[2];
            col := cohort[col].col;
        fi;
         # 
         # End of outer in-lining
         #
        if not IsBound(cohort[col]) then
            return false;
        fi;
        dk[1] := cohort[col][1];
#
# Avoid calling DCEXonEl in the fast case
#
        if false = xrec.action then
            dk[2] := cohort[col][2] *(xrec.Reps[col] mod dk[2]);
        else
            dk[2] := cohort[col][2]*DCEXonEl(u,xrec.Reps[col] mod dk[2],r);
        fi;
    else
        dk[2] := dk[2]*r;
    fi;
    return dk;
end;
       
#############################################################################
##
#F DCEMultiFollow( u, dk, rel, f, b, backwards)    
##
## Follow dk through relator rel, starting with letter f and ending
##  before letter b, going backwards if indicated. Return the first letter
##  not used (ie b for success, f for total failure, etc.
##

DCEMultiFollow := function(u,dk,rel,f, b, backwards)
    local col,d,xrec,ggrec,cohort,r,op,im,i,score,p2s;
    while f <> b do
        d := dk[1];
        r := rel[f];
        if IsInt(r) then
            xrec := u.Xrecs[r];
            ggrec := u.GGrecs[xrec.wgg];
            cohort := d.cohorts[r];
         #
         # First deal with the super-easy case
         #
            if Length(ggrec.dom) = 1 then
                if not IsBound(cohort[1]) then
                    return f;
                fi;
                dk[1] := cohort[1][1];
                if xrec.action = false then
                    dk[2] := cohort[1][2]*dk[2];
                else
                    dk[2] := cohort[1][2]*DCEXonEl(u,dk[2],r);
                fi;
            else
         #
         # Otherwise we need to do some work
         #
                op := ggrec.op;
                col := PositionSorted(ggrec.dom,op(xrec.pt,dk[2]^-1));
            #
            # End of in-lining
            #
                if IsBound(cohort[col]) and IsRec(cohort[col]) then
                    dk[2] := cohort[col].elt*dk[2];
                    col := cohort[col].col;
                fi;
               # 
               # End of outer in-lining
               #
                if not IsBound(cohort[col]) then
                    return f;
                fi;
                dk[1] := cohort[col][1];
#
# Avoid calling DCEXonEl in the fast case
#
                if false = xrec.action then
                    dk[2] := cohort[col][2] *(xrec.Reps[col] mod dk[2]);
                else
                    dk[2] := cohort[col][2]*DCEXonEl(u,xrec.Reps[col] mod dk[2],r);
                fi;
            fi;
        else
            dk[2] := dk[2]*r;
        fi;
        if backwards then
            f := f-1;
        else
            f := f+1;
        fi;
    od;
    return f;
end;
       


#############################################################################
##
#F DCEDegree(u)  The number of single cosets in the table of u
##

DCEDegree := function(u)
    return Sum(List(Filtered(u.table,d->not d.deleted),d->Index(u.K,d.muddle)));
end;


############################################################################
##
#V  DCEWordOps 
##

DCEWordOps := OperationsRecord("DCEWordOps",GroupElementOps);


#############################################################################
##
#F IsDCEWord(x)
##

IsDCEWord := function(x)
    return IsRec(x) and IsBound(x.isDCEWord) and x.isDCEWord;
end;

#############################################################################
##
#F DCEWord(K,x) Convert x into a DCE relation. x may be an element of K
##           or a word in abstract generators
##

DCEWord := function(K,x)
    local r,y,i;
    if not IsGroup(K) then
        Error("Usage DCEWord(K,l)");
    fi;
    if not IsList(x) then
        x := [x];
    fi;
    if not ForAll(x,y->IsWord(y) or y in K) then
        Error("Cannot convert ",x," into a DCE relation");
    fi;
    r := rec();
    r.relext := [];
    for y in x do
        if IsWord(y) and LengthWord(y) > 1 then
            for i in [1..LengthWord(y)] do
                Add(r.relext,Subword(y,i,i));
            od;
        else
            Add(r.relext,y);
        fi;
    od;
    r.isDCEWord := true;
    r.operations := DCEWordOps;
    r.exponent := 1;
    r.K := K;
    return r;
end;


#############################################################################
##
#F DCEWordOps functions
##

DCEWordOps.Print := function(r)
    local rl,i;
    if IsBound(r.name) then 
        Print(r.name);
        return;
    fi;
    Print("DCEWord(",r.K,",[");
    
    for i in [1..Length(r.relext)] do
        rl := r.relext[i];
        if i <> 1 then 
            Print(", ");
        fi;
        if IsPerm(rl) then
            Print(DCEPermString(rl));
        else
            Print(rl);
        fi;
    od;
    Print("])");
    if r.exponent <> 1 then
        Print("^",r.exponent);
    fi;
end;

DCEWordOps.\* := function(r,r1)
    local prod,re,i,lhs;
    if not (IsDCEWord(r) and IsDCEWord(r1)) then
        Error("can only multiply DCE Words with one another");
    fi;
    if r.K <> r1.K then
        Error("Can only multiply DCW Words with the same K");
    fi;
    re := [];
    for i in [1..r.exponent] do
        Append(re,r.relext);
    od;
    for i in [1..r1.exponent] do
        Append(re,r1.relext);
    od;
    prod := DCEWord(r.K,re);
    return prod;
end;

DCEWordOps.\^ := function(r,r1)
    local pow;
    if IsInt(r1) then
        pow := ShallowCopy(r);
        if r1 < 0 then
            pow.relext := List(Reversed(pow.relext),rl->rl^-1);
            r1 := -r1;
        fi;
        pow.exponent := pow.exponent*r1;
        return pow;
    fi;
    return r1^-1*r*r1;
end;


#############################################################################
##
#F DCEPerm(u,word) Make the permutation of the single cosets indicated by word
##



DCEOPerm := function(u,word)
    local s,d,rcs,ims,dk,w,f,d1,cos;
    s := 0;
    dk := [false,false];
    if not IsBound(u.haveSCs) then
        InfoDCE1(u,"Starting To Add Cosets");
        for d in u.table do
            if not d.deleted then
                d.offset := s;
                rcs := Set(RightCosets(u.K,d.muddle));
                d.cosets := rcs;
                InfoDCE2(u,d.name,s);
                s := s+ Size(rcs);
            fi;
        od;
        u.haveSCs := true;
    fi;
    InfoDCE1(u,"Computing action");
    ims := [];
    for d in u.table do
        if not d.deleted then
            for cos in d.cosets do
                dk[1] := d;
                dk[2] := cos.representative;
                DCEMultiFollow(u,dk,word,1,Length(word)+1,false);
                d1 := dk[1];
                Add(ims,d1.offset+Position(d1.cosets,d1.muddle*dk[2]));
            od;
            InfoDCE2(u,d.name);
        fi;
    od;
    return PermList(ims);
end;

            
DCEIPerm := function(u,word)
    local s,d,rcs,ims,dk,w,f,d1,cos,CCE;
    if not IsPermGroup(u.K) then
        return DCEOPerm(u,word);
    fi;
    s := 0;
    dk := [false,false];
    if not IsBound(u.haveSCs) then
        InfoDCE1(u,"Starting To Add Cosets");
        for d in u.table do
            if not d.deleted then
                PermGroupOps.MinimalStabChain(d.muddle);
                CCE := x -> MainEntryCCEPermGroup(u.K,d.muddle,x);
                d.offset := s;
                rcs := Set(List(DCERightCosetReps(u.K,d.muddle),CCE));
                d.cosets := rcs;
                s := s+ Size(rcs);
                InfoDCE2(u,d.name,s);
            fi;
        od;
        u.haveSCs := true;
    fi;
    InfoDCE1(u,"Done cosets, starting image");
    ims := [];
    CCE := function(d,k) return MainEntryCCEPermGroup(u.K,d.muddle,k); end;
    for d in u.table do
        if not d.deleted then
            for cos in d.cosets do
                dk[1] := d;
                dk[2] := cos;
                DCEMultiFollow(u,dk,word,1,Length(word)+1,false);
                d1 := dk[1];
                Add(ims,d1.offset+PositionSorted(d1.cosets,CCE(d1,dk[2])));
            od;
            InfoDCE2(u,d.name);
        fi;
    od;
    return PermList(ims);
end;
            

DCEPerm := function(u,dword)
    local ri,l,p,i;
    if not (IsDCEUniverse(u) and IsDCEWord(dword)) then
        Error("Usage: DCEPerm(u,word)");
    fi;
    if not u.status in ["closed","early-closed"] then
        Error("You can't do this until the enumeration is complete");
    fi;
    if dword.K <> u.K then
        Error("Word and Universe must be over the same K");
    fi;
    ri := [];
    for i in [1..dword.exponent] do
        for l in dword.relext do
            if l in u.K then
                Add(ri,l);
            else
                p := Position(u.xnames,l);
                if p = false then
                    Error("Word contains unknown generator");
                fi;
                Add(ri,u.xnums[p]);
            fi;
        od;
    od;
    return DCEIPerm(u,ri);
end;

        

#############################################################################
##
#F DCEPerms(u) return the permutations corresponding to all the
##  generators of G.
##

DCEPerms := function(u)
    local pgens,x;
    if not IsDCEUniverse(u) or not u.status in ["closed","early-closed"] then
        Error("Argument must be [early-]closed DCE universe");
    fi;
    pgens := [];
    for x in [1..u.nX] do
        if u.Xrecs[x].inverse >= x then
            Add(pgens,DCEIPerm(u,[x]));
        fi;
    od;
    for x in u.K.generators do
        Add(pgens,DCEIPerm(u,[x]));
    od;
    return pgens;
end;


#############################################################################
##
#F DCEPack(u)  pack u, deleteing delerted rows. Should be fast, so do it often
##

DCEPack := function(u)
    local d;
    for d in [1..Length(u.table)] do
        if IsBound(u.table[d]) and u.table[d].deleted then
            Unbind(u.table[d]);
        fi;
    od;
end;


#############################################################################
##
#F DCESuperPack(u)  pack u, recovering blank spaces. Only worth doing
##  this when we have finished, as we are only recovering 4bytes per row.
## 
##

DCESuperPack := function(u)
    local d,ntab;
    ntab := [];
    for d in u.table do
        if not d.deleted then
            Add(ntab,d);
        fi;
    od;
    u.table := ntab;
end;



#############################################################################
##
#F DCERin(u,r) invert relator component r.
##

DCERin := function(u,r)
    if IsInt(r) then
        return u.Xrecs[r].inverse;
    else
        return r^-1;
    fi;
end;
    

#############################################################################
##
#F DCERoot(u,dk) modifies dk to its undeleted image
##

DCERoot := function(u,dk)
    local rep;
    if dk[1].deleted then
        InfoDCE3(u,"Root: ",DCEDKToString(dk));
        rep := DCERoot(u,dk[1].replace); # path compression happens automatically 
        dk[1] := rep[1];
        dk[2] := rep[2]*dk[2];
        InfoDCE3(u," -> ",DCEDKToString(dk));
    fi;
    return dk;
end;

#############################################################################
##
#F DCEInitStacks(u)  Add the stack records to u
##

DCEInitStacks := function(u)
    u.AStack := [];
    u.BStack := [];
    u.CStack := [];
    if u.strategy.Felsch then
        u.Deductions := [];
        u.PDL := [];
        u.nPDLs := 0;
        u.LastPDL := 0;
    fi;
end;



#############################################################################
##
#F DCEStackPD(u,dk,x) Stack a preferred definition
##

DCEStackPD := function(u,dk,x)
    local n;
    n := (u.LastPDL) mod 20 + 1;
    if n = 21 then
        n := 1;
    fi;
    if u.nPDLs < 20 then
        u.nPDLs := u.nPDLs+1;
    fi;
    u.PDL[n] := [dk,x];
    u.LastPDL := n;
    return;
end;


#############################################################################
##
#F DCEClearPDs(u) Clear the preferred definition Stack
##

DCEClearPDs := function(u)
    u.PDL := [];
    u.nPDLs := 0;
    u.LastPDL := 0;
end;
    


#############################################################################
##
#F DCEDedComp(d1,d2) -- the list of rows with outstanding deductions
##                      is ordered by increasing muddle group size (and
##  then used from the tail of the list). The hope is that biug double
##  cosets with small muddle groups will be deleted (or at least
##  contracted) before we get to them.
##

DCEDedComp := function(d1,d2) 
    local s1,s2;
    s1 := Size(d1.muddle);
    s2 := Size(d2.muddle);
    if s1 = s2 then
        return d1.name > d2.name;
    else
        return s1 > s2;
    fi;
end;


#############################################################################
##
#F DCEStackDed(u,d,x,col,back) Record that the entry at d.cohorts[x][col] 
##  has been filled in
##  and is due to be scanned. m is the associated group
##

DCEStackDed := function(u,d,x,col,back)
    local i,olded;
    olded := u.Deductions;
    i := PositionSorted(olded,d,DCEDedComp);
    if i > Length(olded) or not IsIdentical(olded[i],d) then
        u.Deductions := olded{[1..i-1]};
        Add(u.Deductions,d);
        Append(u.Deductions,olded{[i..Length(olded)]});
    fi;
    d.deductions[x][col] := back;
end;

#############################################################################
##
#F DCEStackC(u,dk,mgens)  Stack the coincidence mgens fix d. 
##

DCEStackC := function(u,dk,mgens)
    if mgens = u.id or IsList(mgens) and ForAll(mgens,x->x=u.id) then
        return;
    fi;
    Add(u.CStack,dk);
    Add(u.CStack,mgens);
end;



#############################################################################
##
#F DCEStackA(u,dk1,dk2)  Stack the coincidence dk1 = dk2
##

DCEStackA := function(u,dk1,dk2)
    local k,s;
    if IsIdentical(dk1[1],dk2[1]) then
        if dk1[2] <> dk2[2] then
            k := dk1[2]/dk2[2];
            InfoDCE3(u,"Type A ",DCEDKToString(dk1),DCEDKToString(dk2),
                    "recast as C",DCEPermString(k));
            dk1[2] := u.id;
            DCEStackC(u,dk1,k);
        fi;
        return;
    fi;
    Add(u.AStack,dk1);
    Add(u.AStack,dk2);
end;

#############################################################################
##
#F DCEStackB(u,dk1,x,dk2)  Stack the coincidence dk1.x = dk2
##

DCEStackB := function(u,dk1,x,dk2)
    Add(u.BStack,dk1);
    Add(u.BStack,x);
    Add(u.BStack,dk2);
end;



#############################################################################
##
#F DCEEmptyOne(u,d,x,col) empty the col entry in the x cohort of the 
## d row onto the B stack
##
## NB this does not increment u.blanks for the entry it was told to
##  delete as that is always about to be destroyed, but does for the
##  other end if relevant
##
DCEEmptyOne := function(u,d,x,col)
    local entry,col1,xi,d1,cohort,xrec,cohort1;
    cohort := d.cohorts[x];
    entry := cohort[col];
    Unbind(cohort[col]);
    Unbind(d.deductions[x][col]);
    if IsList(entry) then
        xrec := u.Xrecs[x];
        DCEStackB(u,[d,xrec.Reps[col]],x,entry);
        xi := xrec.inverse;
        col1 := DCEColumn(u,entry,xi);
        d1 := entry[1];
        cohort1 := d1.cohorts[xi];
        if IsBound(cohort1[col1]) then
            Unbind(cohort1[col1]);
            Unbind(d1.deductions[xi][col1]);
            u.blanks := u.blanks+1;
        fi;
    fi;
end;


#############################################################################
##
#F DCEProcessA(u,dk1,dk2)  Process the coincidence dk1 = dk2
##

DCEProcessA := function(u,dk1,dk2)
    local d1,d2,k,m,row,dk,i,j,l,x,cohort,xrec;
    InfoDCE2(u,"Processing A",DCEDKToString(dk1),DCEDKToString(dk2));
    dk1 := DCERoot(u,dk1);
    dk2 := DCERoot(u,dk2);
    if IsIdentical(dk1[1],dk2[1]) then
        DCEStackC(u,[dk1[1],u.id],dk1[2]/dk2[2]);
        return;
    elif dk1[1].name> dk2[1].name then
        d1 := dk2[1];
        dk2[2] := dk2[2]/dk1[2];
        d2 := dk1[1];
        dk := dk2;
    else
        d1 := dk1[1];
        dk1[2] := dk1[2]/dk2[2];
        d2 := dk2[1];
        dk := dk1;
    fi;
## 
## Here d1 < d2 and both are undeleted and the
## coincidence is d2 = dk    
##    
    for x in [1..u.nX] do
        cohort := d2.cohorts[x];
        xrec := u.Xrecs[x];
        for i in [1..xrec.CohortSize] do
            if IsBound(cohort[i]) then
                DCEEmptyOne(u,d2,x,i);
            else
                u.blanks := u.blanks-1; 
            fi;
        od;
    od;
    m := d2.muddle;
    u.degree := u.degree - Index(u.K,m);
    u.dcct := u.dcct - 1;
    DCEStackC(u,dk,m.generators);
    if u.strategy.Felsch then
        Unbind(d2.deductions);
        i := PositionSorted(u.Deductions,d2,DCEDedComp);
        if i <= Length(u.Deductions) and IsIdentical(u.Deductions[i],d2) then
            l :=Length(u.Deductions); 
            for j in [i+1..l] do
                u.Deductions[j-1] := u.Deductions[j];
            od;
            Unbind(u.Deductions[l]);
        fi;
        DCEClearPDs(u);
    fi;
    d2.deleted := true;
    d2.replace  := dk;
    Unbind(d2.cohorts);
    Unbind(d2.muddle);
end;




#############################################################################
##
#F DCEProcessB(u,dk1,x,dk2)  Process the coincidence d1.k1.x = d2.k2
##

DCEProcessB := function(u,dk1,x,dk2)
    local d1,d2,k1,k2,m1,m2,flag,mm,xi,col1,col2,im1,im2,l1,h,
          k01,k02,l2,ki1,ki2,mg,mg1,cohort1,cohort2;
    InfoDCE2(u,"Processing B",DCEDKToString(dk1),x,DCEDKToString(dk2));
## First we replace everything by undeleted images
    dk1 := DCERoot(u,dk1);
    dk2 := DCERoot(u,dk2);
## Now load up a few locals    
    d1 := dk1[1];
    d2 := dk2[1];
    xi := u.Xrecs[x].inverse;
    col1 := DCEColumn(u,dk1,x);
    col2 := DCEColumn(u,dk2,xi);
    k01 := u.Xrecs[x].Reps[col1];
    k02 := u.Xrecs[xi].Reps[col2];
    k1 := dk1[2];
    k2 := dk2[2];
    l1 := LeftQuotient(k01,k1);
    l2 := LeftQuotient(k02,k2);
    ki1 := k2/DCEXonEl(u,l1,x);
    ki2 := k1/DCEXonEl(u,l2,xi);
    
# The input should imply that d1.k01.x = d2.ki1 and d2.k02.xi = d1.ki2    
# col1 and col2 are the columns relevant to these two facts   
    
#
# First we check for the entries being alreadu full.
#
    
    flag := false;
    cohort1 := d1.cohorts[x];
    cohort2 := d2.cohorts[xi];
    if IsBound (cohort1[col1]) then
        im1 := ShallowCopy(cohort1[col1]);
        dk2[2] := ki1;
        InfoDCE3(u,"Stacked A in first place",DCEDKToString(im1),DCEDKToString(dk2));
        DCEStackA(u,im1,dk2);
        flag := true;
    fi;
    if IsBound (cohort2[col2]) then
        im2 := ShallowCopy(cohort2[col2]);
        dk1[2] := ki2;
        InfoDCE3(u,"Stacked A in second place",DCEDKToString(im2),
                DCEDKToString(dk1));
        DCEStackA(u,im2,dk1);
        flag := true;
    fi;
    if flag then 
        return;
    fi;
    
## Now we check the consistency conditions    
    m1 := d1.muddle;
    m2 := d2.muddle;
    
    if not IsIdentical(m1,u.triv) then
        mm := DCELin(u,m1,x,k01^-1);
        mg := List(mm.generators,g -> DCEXonEl(u,g^k01,x)^(ki1^-1));
        mg := Filtered(mg,g-> not g in m2);
        if Length(mg) <> 0 then
            flag := true;
            DCEStackC(u,[d2,u.id],mg);
        fi;
    fi;
    if not IsIdentical(m2,u.triv) then
        mm := DCELin(u,m2,xi,k02^-1);
        mg := List(mm.generators,g -> DCEXonEl(u,g^k02,xi)^(ki2^-1));
        mg := Filtered(mg,g-> not g in m1);
        if Length(mg) <> 0 then
            flag := true;
            DCEStackC(u,[d1,u.id],mg);
        fi;
    fi;
    
    if flag then
        dk1[2] := k01;
        dk2[2] := ki1;
        DCEStackB(u,dk1,x,dk2);
        InfoDCE3(u,"Restacked");
        return;
    fi;
    
    
    #
    # Here we have really got new information to put into the table
    #
    #
    # If we are working Felsch style then we have to scan the relators for
    # any possible conclusions to be derived from this definition
    #
    
    if IsIdentical(d1,d2) and x = xi and col1 = col2 then
        h := k1*DCEXonEl(u,LeftQuotient(l2,l1),x)/k2;
        DCEStackC(u,[d1,u.id],h); 
        dk2[2] := ki1;
        cohort1[col1] := dk2;
        InfoDCE3(u,"Entries in same place ",d1.name,col1,h,DCEDKToString(dk2));
        if u.strategy.Felsch then
            DCEStackDed(u,d1,x,col1,false);
        fi;
        u.blanks := u.blanks-1;
    else
        dk2[2] := ki1;
        cohort1[col1] := dk2;
        dk1[2] := ki2;
        cohort2[col2] := dk1;
        InfoDCE3(u,"Two new entries",d1.name,col1,DCEDKToString(dk2),"\n#I  ",
                d2.name,col2,DCEDKToString(dk1));
        u.blanks := u.blanks - 2;
        if u.strategy.Felsch then
            DCEStackDed(u,d1,x,col1,false);
            DCEStackDed(u,d2,xi,col2,true);
        fi;
    fi;
end;



#############################################################################
##
#F DCECompBetas(u,m,ggr)  Compute the beta entries needed for
##                        muddle group m, gain group encoded by
##                        ggr
##   returns a list with unbound positions for the alpha entries
##
#F DCEExtendBetas avoids recalculating things that were known already
##

DCEExtendBetas := function(u,d,x,newgens)
    local op,xrec,ggr,q,cohort,ix,i,nb,alphas,isalpha,unused,seed,pt,dom,
          im,gen,neworb,seed1,j,el,nleft,elt,orbs,score,p2s,p2list;
    xrec := u.Xrecs[x];
    ggr := u.GGrecs[xrec.wgg];
    dom := ggr.dom;
    op := ggr.op;
    cohort := d.cohorts[x];
    ix := [1..xrec.CohortSize];
    alphas := Filtered(ix,i->not IsBound(cohort[i]) or IsList(cohort[i]));
    IsSet(alphas);
    isalpha := BlistList(ix,alphas);
    unused := BlistList(ix,ix);
    nb := [];
    nleft := xrec.CohortSize;
    orbs := [];
    for i in ix do
        if isalpha[i] then
            j := i;
        else
            j := cohort[i].col;
        fi;
        if IsBound(orbs[j]) then
            Add(orbs[j],i);
        else
            orbs[j] := [i];
        fi;
    od;
    for seed in alphas do
        if unused[seed] then
            q := orbs[seed];
            for j in q do
                if j <> seed then
                    nb[j] := cohort[j];
                fi;
                unused[j] := false;
            od;
            nleft := nleft - Length(q);
            for i in q do                
                if nleft = 0 then
                    return nb;
                fi;
                pt := dom[i];
                if i = seed then
                    elt := u.id;
                else
                    elt := nb[i].elt;
                fi;
                for gen in newgens do
                    im := PositionSorted(dom,op(pt,gen));
                    if unused[im] then
                        if isalpha[im] then
                            seed1 := im;
                            el := u.id;
                        else
                            seed1 := cohort[im].col;
                            el := cohort[im].elt^-1;
                        fi;
                        el := elt*gen*el;
                        neworb := orbs[seed1];
                        for j in neworb do
                            if j <> seed1 then
                                nb[j] := cohort[j];
                                nb[j].col := seed;
                                nb[j].elt := el*nb[j].elt;
                            else
                                nb[j] := rec(col := seed, elt := el);
                            fi;
                            unused[j] := false;
                        od;
                        nleft := nleft-Length(neworb);
                        Append(q,neworb);
                    fi;
                od;
            od;
        fi;
    od;
end;
                        

DCECompBetas := function(u,m,ggr)
    local op,q,l,pt,gen,gens,seed,dom,im,rep,i,inds,pt0,pt1,
          nleft,unused,score,p2s,p2list;
    op := ggr.op;
    dom := ggr.dom;
    pt0 := PositionSorted(dom,ggr.pt0);
    inds := [1..Length(dom)];
    l := [];
    gens := m.generators;
    nleft := Length(dom);
    unused := BlistList([1..nleft],[1..nleft]);
    for seed in inds do
        if unused[seed] then
            q := [seed];
            nleft := nleft-1;
            unused[seed] := false;
            for pt in q do
                if nleft = 0 then
                    return l;
                fi;
                if pt = seed then
                    rep := u.id;
                else
                    rep := l[pt].elt;
                fi;
                pt1 := dom[pt];
                for gen in gens do
                    im := PositionSorted(dom,op(pt1,gen));
                    if unused[im] then
                        l[im] := rec(elt := rep*gen, col := seed);
                        Add(q,im);
                        nleft := nleft-1;
                        unused[im] := false;
                    fi;
                od;
            od;
        fi;
    od;
end;
                

#############################################################################
##
#F DCEProcessC(u,d,k)  Process the coincidence k fixes dk
##

DCEProcessC := function(u,dk,k)
    local newm,d,s,newbetas,nb,xrec,co,cs,i,wgg,col,k1,dedp,l,x,m,cohort;
    dk := DCERoot(u,dk);
    k1 := dk[2];
    k := k1*k/k1;
    d := dk[1];
    m := d.muddle;
    if not IsList(k) then
        k := [k];
    fi;
    k := Filtered(k, x-> not x in m);
    if Length(k) = 0 then
        return;
    fi;
    if u.strategy.Felsch then
        dedp := PositionSorted(u.Deductions,d,DCEDedComp);
        l := Length(u.Deductions);
        if dedp <= l and IsIdentical(u.Deductions[dedp],d) then
            for i in [dedp+1..Length(u.Deductions)] do
                u.Deductions[i-1] := u.Deductions[i];
            od;
            Unbind(u.Deductions[l]);
        fi;
    fi;
    u.degree := u.degree - Index(u.K,m);
    Unbind(m.elements);
    newm := m;
    for i in k do
        newm := Closure(newm,i);
    od;
    InfoDCE2(u,"Processing C",DCEDKToString(dk),List(k,x->DCEPermString(x)));    
    if newm = u.K then
        newm := u.K;
    fi;
    d.muddle := newm;
    u.degree := u.degree + Index(u.K,newm);
    InfoDCE3(u,"   new muddle group",newm);
    #
    # Now sort out the beta-entries
    #
    newbetas := [];
    for x in [1..u.nX] do
        xrec := u.Xrecs[x];
        cs := xrec.CohortSize;
        wgg := xrec.wgg;
        cohort := d.cohorts[x];
        #
        # Cohorts with only one alpha entry left can be ignored
        #
        if ForAny([2..cs],col->(not IsBound(cohort[col])) or IsList(cohort[col])) then
            if not IsBound(newbetas[wgg]) then
                if IsIdentical(m,u.triv) then
                    newbetas[wgg] := DCECompBetas(u,newm,u.GGrecs[wgg]);
                else
                    newbetas[wgg] := DCEExtendBetas(u,d,x,k);
                fi;
            fi;
            nb := newbetas[wgg];
            for i in [2..cs] do
                if IsBound(nb[i]) then
                    if IsBound(cohort[i]) then
                        if IsList(cohort[i]) then
                            DCEEmptyOne(u,d,x,i);
                        fi;
                    else
                        u.blanks := u.blanks-1;
                    fi;
                    cohort[i] := nb[i];
                fi;
            od;
        fi;
        #
        # We have to propagate the muddle group info along any
        # edges that we haven't removed.
        #
        for col in [1..cs] do
            if IsBound(cohort[col]) and IsList(cohort[col]) then
                l := xrec.Reps[col];
                if xrec.action = false then
                    DCEStackC(u,ShallowCopy(cohort[col]),
                            OnTuples(DCELin(u,newm,x,l^-1).generators,l));
                else
                    DCEStackC(u,ShallowCopy(cohort[col]),
                            List(DCELin(u,newm,x,l^-1).generators,
                                 g->DCEXonEl(u,g^l,x)));
                fi;
            fi;
        od;
    od;
    if u.strategy.Felsch then
        DCEClearPDs(u);
        if ForAny(d.deductions,x->Length(x) <> 0) then
            dedp := PositionSorted(u.Deductions,d,DCEDedComp);
            for i in [Length(u.Deductions),Length(u.Deductions)-1..dedp] do
                u.Deductions[i+1] := u.Deductions[i];
            od;
            u.Deductions[dedp] := d;
        fi;
    fi;
end;



#############################################################################
##
#F DCEDoStacks(u)  Clear the coincidence stacks
##

DCEDoStacks := function(u)
    local dk,dk1,x,dk2,m,l;
    while true do
        if Length(u.CStack) <> 0 then
            l := Length(u.CStack);
            dk := u.CStack[l-1];
            m := u.CStack[l];
            Unbind(u.CStack[l-1]);
            Unbind(u.CStack[l]);
            DCEProcessC(u,dk,m);
        elif Length(u.AStack) <> 0 then
            l := Length(u.AStack);
            dk1 := u.AStack[l-1];
            dk2 := u.AStack[l];
            Unbind(u.AStack[l]);
            Unbind(u.AStack[l-1]);
            DCEProcessA(u,dk1,dk2);
        elif Length(u.BStack) <> 0 then
            l := Length(u.BStack);
            dk1 := u.BStack[l-2];
            x := u.BStack[l-1];
            dk2 := u.BStack[l];
            Unbind(u.BStack[l]);
            Unbind(u.BStack[l-1]);
            Unbind(u.BStack[l-2]);
            DCEProcessB(u,dk1,x,dk2);
        else
            return true;
        fi;
    od;
end;


#########a####################################################################
##
#F DCEScan(u,d,x,col,edpgp,m)  This does the Felsch style scanning
##  for the positions in edpgp
##   when the col entry in row d has just been filled in.
##
##  m is the group DCELin(u,d.muddle^u.Xrecs[x].Reps[col],x) which has
## always been pre-calculated by the time we get called
##
## flipped is true when this is the 'back half of a deduction that
## created two new entries. In this case, certain relators can be ignored
##

DCEScan := function(u,d,x,col,eg,m)
    local k0,e,todo,lx,l,dkf,dkb,f,b,r,ri,n,iri,rt,s,i,
          flag,ct,es,m1,xrec,cohort,flag,kx;
    rt := false;
    cohort := d.cohorts[x];
    xrec := u.Xrecs[x];
    k0 := xrec.Reps[col];
    lx := xrec.L;
    n := eg[1];
    m := AsSubgroup(Parent(lx),m);
    m1 := d.muddle;
    if Index(lx,m)/RootInt(Size(n)) >= 200 then
        es := Filtered(eg[2],e->not e.ismulti);
    else
        es := eg[2];
    fi;
    if Length(es)= 0 then
        return false;
    fi;
    todo := DCEDoubleCosetReps(lx,m,n);
    for e in es do
        r := e.rel;
        InfoDCE3(u,"Scanning in ",r,"Starting at ",e.pos,
                Length(todo)," cases");
        dkf := []; dkb := [];
        for l in [1..Length(todo)] do
            dkf[1] := d;
            dkf[2] := k0*todo[l];
            dkb[1] := d;
            dkb[2] := dkf[2];
            f := e.pos;
            b := e.pos+r.intlen;
            ri := r.relint3;
            iri := r.irelint3;
            flag := true;
            f := DCEMultiFollow(u,dkf,ri,f,b,false);
            if f < b then
                b := DCEMultiFollow(u,dkb,iri,b-1,f-1,true)+1;
            fi;
            if f + 3 >= b then
                if f = b then
                    InfoDCE3(u,"Got round, got coincidence");
                    DCEStackA(u,ShallowCopy(dkf),ShallowCopy(dkb));
                elif f + 1 = b then
                    InfoDCE3(u,"Got deduction");
                    DCEStackB(u,ShallowCopy(dkf),ri[f],ShallowCopy(dkb));
                elif f + 2 = b or not IsInt(ri[f+1]) then
                    DCEStackPD(u,ShallowCopy(dkf),ri[f]); 
                fi;
            fi;
            DCEDoStacks(u);
            if d.deleted or not IsList(cohort[col]) then 
                return false;
            fi;
            if not IsIdentical(d.muddle,m1) and Length(todo)-l > 3 then
               return  true;
            fi;
        od;
    od;
    return false;
end;
   

#############################################################################
##
#F DCEDoubleCosetRepstabs(u,g,h,k,r) gets double cosets reps and pt stabilisers
##  in h. r is a set of Canonical right coset reps for k in g. 
## 
 
DCEDoubleCosetRepStabs := function(u,g,h,k,r,ir)
    local sg,MyCCE,results,bl,rep,st,seed,q,gens,gen,el,pt,im,results,reps,iseed,
          sch,sgf,uk,stg,sh,stfull,ofull,i,todo,norm,rt;
    sg := Size(g);
    sh := Size(h);
    if sg = sh then
        return [[g.identity,k]];
    fi;
    if sg = Size(k) then
        return [[g.identity,h]];
    fi;
    if sh = 1 then
        if Length(k.generators) = 0 then
            return List(Elements(g),x->[x,h]);
        else
            return List(ir,x->[x,h]);
        fi;
    else
        if Length(k.generators) = 0 then
            return List(DCERightCosetReps(g,h),x->[x,k]);
        else
            #
            # Here we go seriously to work, inspired by the double coset code
            #
            #
            # First we must make sure that k is has the canonical base
            #
            if not IsBound( k.smallestBase )  or k.smallestBase <> Base( k)  then
                MakeStabChain( k, k.operations.MovedPoints( k ) );
                k.smallestBase := Base( k );
            fi;
            MyCCE:=g.operations.MainEntryCCE;
            bl := BlistList(r,[]);
            gens := h.generators;
            results := [];
            uk := u.K;
            sgf := uk.operations.Subgroup;
            norm :=  IsNormal(g,h);
            if norm then 
                todo := [1];
            else
                todo := [1..Length(r)];
            fi;
            for seed in todo do
                if not bl[seed] then
                    iseed := ir[seed];
                    st := u.triv;
                    stg := [];
                    reps := [];
                    q := [seed];
                    bl[seed] := true;
                    reps[seed] := h.identity;
                    stfull := false;
                    ofull := false;
                    i := 1;
                    while not ofull do
                        pt := q[i];
                        i := i+1;
                        rep := reps[pt];
                        el := r[pt];
                        for gen in gens do
                            im := PositionSorted(r,MyCCE(g,k,el*gen));
                            if not bl[im] then
                                Add(q,im);
                                bl[im] := true;
                                reps[im] := rep*gen;
                            elif not stfull then
                                sch := (rep*gen/reps[im])^iseed;
                                if sch <> u.id and not sch in st then 
                                    Add(stg,sch);
                                    st := sgf(uk,stg);
                                    stfull := 2*Size(st)*Length(q) > sh;
                                fi;
                            fi;
                        od;
                        #
                        # Check if we have everything
                        #
                        ofull := Length(q)*Size(st) = sh;
                    od;
                    Add(results,[iseed,st]);
                fi;
            od;
            if norm then
                rt := DCERightCosetReps(g,h);
                for seed in [2..Length(r)] do
                    if not bl[seed] then
                        el := r[seed];
                        Add(results,[ir[seed],st]);
                        for pt in rt do
                            bl[PositionSorted(r,MyCCE(g,k,el*pt))] := true;
                        od;
                    fi;
                od;
            fi;
                    
            return results;
        fi;
    fi;
end;



#############################################################################
##
#F DCEDoMultiScan(u,d,x,col,multi,m) Use MultiScan to see whether any information
##  can be gained from a relator and process it accordingly
##

DCEDoMultiScan := function(u,d,x,col,multi,m)
    local l,k0,r,rf,g,results,status,i,cos,rl,h,newh,split,m1,dk,ll,cos1,dk1,
          j,sg,cos2,mm,nstatus,dkf,dkb,f,b,ri,iri,ci,cj,bgps,oldm,xrec,cohort,
          reps,ireps,spl,breps,ibreps,flag,kx,cases;
    xrec := u.Xrecs[x];
    cohort := d.cohorts[x];
    l := xrec.L;
    m := AsSubgroup(Parent(l),m);
    g := multi.gp;
    if Index(l,m)/RootInt(Size(g)) < 200 then
        return false;
    fi;
    r := multi.rel;
    ri := r.relint3;
    iri := r.irelint3;
    mm := AsSubgroup(u.K,m);
    k0 := xrec.Reps[col];
    oldm := d.muddle;
    InfoDCE2(u,"Multi-scanning",d.name,u.Xrecs[x].name,col);
    InfoDCE3(u,"Would have been",Length(DCEDoubleCosetReps(l,m,g)));
    rf := multi.pos;
    results := [];
    status := [[u.id,mm,[d,k0]]];
    i := 1;
    while i < multi.limit and Length(status) <> 0 do
        rl := ri[rf+i-1];
        if IsInt(rl) then
            h := multi.fgps[i];
            nstatus := [];
            newh := multi.fgps[i+1];
            reps := multi.freps[i];
            ireps := multi.ifreps[i];
            for cos in status do
                m1 := cos[2];
                dk := cos[3];
                split := DCEDoubleCosetRepStabs(u,h,m1,newh,reps,ireps);
                if Length(split) = 1 then
                    cos[2]:= split[1][2];
                    if false = DCEFollow(u,dk,rl) then
                        Add(cos,i);
                        Add(results,cos);
                    else
                        Add(nstatus,cos);
                    fi;
                else
                    for spl in split do
                        ll := spl[1];
                        dk1 := ShallowCopy(dk);
                        dk1[2] := dk1[2]*ll;
                        cos1 := [cos[1]*ll,spl[2],dk1];
                        if false = DCEFollow(u,dk1,rl) then
                            Add(cos1,i);
                            Add(results,cos1);
                        else
                            Add(nstatus,cos1);
                        fi;
                    od;
                fi;
            od;
            status := nstatus;
        else
            for cos in status do
                cos[1] := cos[1]^rl;
                cos[2] := cos[2]^rl;
                cos[3][2] := cos[3][2]*rl;
            od;
        fi;
        i := i+1;
    od;
    InfoDCE2(u,"Got forward to",i,Length(status),"still going,",
          Length(results),"blocked");
    #
    # At this point, results contains the stuff that stopped along the way
    # while status contains the stuff that got past the point where multi-scanning
    # does any good. 
    #
    # First we deal with the stuff that's still going
    #
    dkb := [];
    ci := multi.conjs[i];
    for cos in status do
        dkf := cos[3];
        DCERoot(u,dkf);
        dkb[1] := d;
        dkb[2] := k0*cos[1]^ci;
        f := rf+i-1;
        b := multi.pos+r.intlen;
        f := DCEMultiFollow(u,dkf,ri,f,b,false);
        if f < b then
            b :=  DCEMultiFollow(u,dkb,iri,b-1,f-1,true)+1;
        fi;
        if f + 3 >= b then
            if f +1 >= b then
                if f = b then
                    InfoDCE3(u,"Got round, got coincidence");
                    DCEStackA(u,ShallowCopy(dkf),ShallowCopy(dkb));
                else
                    InfoDCE3(u,"Got deduction");
                    DCEStackB(u,ShallowCopy(dkf),ri[f],ShallowCopy(dkb));
                fi;
                DCEDoStacks(u);
                if d.deleted or not IsList(cohort[col]) then
                    return false;
                fi;
                if not IsIdentical(oldm,d.muddle) then
                    return true;
                fi;
            elif f + 2 = b or not IsInt(ri[f+1]) then
                DCEStackPD(u,ShallowCopy(dkf),ri[f]); 
            fi;
        fi;
    od;
    InfoDCE2(u,"Done the ones that continued");
    #
    # Finally things that got blocked up
    #
    sg := Size(g);
    cases := 0;
    for cos in results do
        i := cos[4];
        ci := multi.conjs[i];
        j := r.intlen;
        ll := cos[1]^ci;
        bgps := multi.bgps[i];
        breps := multi.breps[i];
        ibreps := multi.ibreps[i];
        h := bgps[j];
        status := [[u.id,cos[2]^ci,[d,k0*ll]]];
        while Length(status) <> 0 and Size(h) > sg do
            rl := iri[multi.pos + j-1];
            if IsInt(rl) then
                nstatus := [];
                newh := bgps[j-1];
                reps := breps[j-1];
                ireps := ibreps[j-1];
                for cos1 in status do
                    m1 := cos1[2];
                    split := DCEDoubleCosetRepStabs(u,h,m1,newh,reps,ireps);
                    if Length(split) = 1 then
                         if false <> DCEFollow(u,cos1[3],rl) then
                            cos1[2] := split[1][2];
                            Add(nstatus,cos1);
                        else
                            cases := cases+1;
                        fi;
                    else
                        for spl in split do
                            ll := spl[1];
                            dk1 := ShallowCopy(cos1[3]);
                            dk1[2] := dk1[2]*ll;
                            if false <> DCEFollow(u,dk1,rl) then
                                cos2 := [cos1[1]*ll,spl[2],dk1];
                                Add(nstatus,cos2);
                            else
                                cases := cases+1;
                            fi;
                        od;
                    fi;
                od;
                status := nstatus;
                h := newh;
            else
                for cos1 in status do
                    cos1[1] := cos1[1]^rl;
                    cos1[2] := cos1[2]^rl;
                    cos1[3][2] := cos1[3][2]*rl;
                od;
                h := bgps[j-1];  
            fi;
            j := j-1;
        od;
        #
        # Finally, anything left in status is still going when the group
        # at the back died, so we need to scan it normally
        #
        cases := cases+Length(status);
        cj := (multi.conjs[j+1]/multi.conjs[r.intlen+1])/multi.conjs[i];
        for cos1 in status do
            DCERoot(u,cos[3]);
            dkf := ShallowCopy(cos[3]);
            DCERoot(u,cos1[3]);
            dkb := cos1[3]; 
            dkf[2] := dkf[2]*cos1[1]^cj;
            f := rf+i-1;
            b := rf+j;
            f := DCEMultiFollow(u,dkf,ri,f,b,false);
            if f < b then
                b := 1+ DCEMultiFollow(u,dkb,iri,b-1,f-1,true);
            fi;
            if f + 3 >= b then
                if f +1 >= b then
                    if f = b then
                        InfoDCE3(u,"Got round, got coincidence");
                        DCEStackA(u,dkf,dkb);
                    else
                        InfoDCE3(u,"Got deduction");
                        DCEStackB(u,dkf,ri[f],dkb);
                    fi;
                    DCEDoStacks(u);
                    if d.deleted or not IsList(cohort[col]) then
                        return false;
                    fi;
                    if not IsIdentical(d.muddle,oldm) then
                        return true;
                    fi;
                elif f + 2 = b or not IsInt(ri[f+1]) then
                    DCEStackPD(u,dkf,ri[f]); 
                fi;
            fi;
        od;
    od;
    InfoDCE2(u,"Total of",cases,"cases");
    return false;
end;
                    
         
#############################################################################
##
#F DCEProcessDed(u,d,col,back) This does all the scanning corresponding to
##  the col entry in row d
##


DCEProcessDed := function(u,d,x,col,back)
    local k0,m,edps,e,eg,multi,flag,xrec,cohort;
    xrec := u.Xrecs[x];
    cohort := d.cohorts[x];
    if d.deleted or not IsList(cohort[col]) then
        return;
    fi;
    InfoDCE1(u,"Processing deduction:",d.name,u.Xrecs[x].name,"column",col,back);
    k0 := xrec.Reps[col];
    edps := u.Edps[x];
    if back then
        edps := Filtered(List(edps,eg-> [eg[1], Filtered(eg[2],e->not e.flip)]),
                        eg-> Length(eg[2]) > 0);
        if Length(edps) = 0 then 
            return;
        fi;
    fi;
    m := DCELin(u,d.muddle,x,k0^-1)^k0;
    for eg in edps do
        flag := true;
        while flag do
            flag := DCEScan(u,d,x,col,eg,m);
        od;
        if d.deleted or not IsList(cohort[col]) then
            return;
        fi;
    od;
    for multi in u.Multis[x] do
        if not back or not multi.flip then
            flag := true;
            while flag do
                flag := DCEDoMultiScan(u,d,x,col,multi,m);
                if d.deleted or not IsList(cohort[col]) then
                    return;
                fi;
            od;
        fi;
    od;
end;



#############################################################################
##
#F DCEDefine(u,dk,x)   Make a definition to be dk.x. Stacks an
##  appropriate type B and returns the value of dkx
##

DCEDefine := function(u,dk,x)
    local nrow,ndk;
    if not dk[1].deleted then
        nrow := rec( muddle := u.triv, name := u.nextD,
                     cohorts := List([1..u.nX],x->[]),
                     deductions := List([1..u.nX],x->[]),
                     deleted := false,
                     operations := DCERowOps);
        u.table[u.nextD] := nrow;
        u.nextD := u.nextD+1;
        u.blanks := u.blanks + u.nCols;
        u.degree := u.degree+Size(u.K);
        u.dcct := u.dcct +1;
        ndk := [nrow,u.id];
        DCEStackB(u,dk,x,ShallowCopy(ndk));
        return ndk;
    fi;
end;



#############################################################################
##
#F DCEPush(u,d,k,r) Push the relator r from the single coset d.k
##

DCEPush := function(u,d,k,r)
    local dkf,dkb,f,b,lr,rl,ndk,nrow,m,ri,iri;
    InfoDCE2(u,d.name,DCEPermString(k));
    dkf := [d,k];
    dkb := ShallowCopy(dkf);
    f := 0;
    ri := r.relint;
    iri := r.irelint;
    lr := r.intlen;
    #
    # First go forwards as far as we can
    #
    b := lr;
    f := DCEMultiFollow(u,dkf,ri,f+1,b+1,false)-1;
    if f < b then
        b := DCEMultiFollow(u,dkb,iri,b,f,true);
    fi;
    InfoDCE2(u,f,b);
    if b = f then
        DCEStackA(u,dkf,dkb);
    else
        for f in [f+1..b-1] do
            rl := ri[f];
            if IsInt(rl) then
                dkf := DCEDefine(u,dkf,rl);
            else
                dkf[2] := dkf[2]*rl;
            fi;
        od;
        DCEStackB(u,dkf,ri[b],dkb);
    fi;
    DCEDoStacks(u);
end;


#############################################################################
##
#F DCEDPush(u,ds,r) Push relator number r from a double coset
##
## This dunction is split into two to get around the lack of break 
## in GAP. The first function returns true if it finished or
## false if it didn't.
##


DCEDPush1 := function(u,d,r)
    local n,d,m,todo,i,l;
    n := r.group;
    if d.deleted then
        return true;
    fi;
    m := d.muddle;
    InfoDCE1(u,d.name,r);
    todo := DCEDoubleCosetReps(u.K,m,n);
    #
    # So when we get here we know what we have to do.
    #
    i := 1;
    l := Length(todo);
    InfoDCE1(u," ",l,"cases");
    while true do
        DCEPush(u,d,todo[i],r);
        if i = l or d.deleted then
            return true;
        fi;
        if not IsIdentical(d.muddle,m) then
            return false;
        fi;
        if  u.blanks = 0 and u.degree in u.strategy.EC then
            u.EarlyClosed := true;
            u.status := "early-closed";
            return true;
        fi;
        i := i+1;
    od;
end;    



DCEDPush := function(u,d,r)
    local flag;
    flag := false;
    while not flag do
        flag := DCEDPush1(u,d,r);
    od;
    if u.dcct > 10 and u.blanks < u.dcct and u.degree in u.strategy.EC then
        u.PreEarlyClosed := true;
        u.status := "in end game";
    fi;
end;




#############################################################################
##
#F DCEHandlePEC(u)  Handle Pre-early closure that is a nearly full table.
##                  do this by pushing all relators in order of inclreasing
##                  weight from all rows with a blank entry.
##                  We hope this produces early-closing.

DCEHandlePEC := function(u)
    local brows, sortedrels, cohortsizes, xs,r,d;
    InfoDCE1(u,"Entering Pre-early closing", u.dcct, u.degree, u.blanks);
    xs := [1..u.nX];
    cohortsizes := List(xs,x->u.Xrecs[x].CohortSize);
    sortedrels := ShallowCopy(u.relators);
    Sort(sortedrels, function(r1,r2) return r1.weight < r2.weight; end);
    brows := 
      Filtered(u.table, d-> not(d.deleted) and ForAny(xs,x-> Number(d.cohorts[x]) <> cohortsizes[x]));    
    for d in brows do
        for r in sortedrels do
            if d.name >= r.next then
                DCEDPush(u,d,r);
                if u.EarlyClosed then
                    u.PreEarlyClosed := false;
                    return;
                fi;
            fi;
        od;
    od;
    u.PreEarlyClosed := false;
    InfoDCE1(u,"Leaving Pre-early closing", u.dcct, u.degree, u.blanks);
end;



#############################################################################
##
#F DCEPushAll(u,w)   Push everything at weight w
##

DCEPushAll := function(u,w)
    local r,d,flag,top;
    if not ForAny(u.relators, r-> w > r.weight and 
               r.next < u.FirstOfWeight[w-r.weight+1]) then
        InfoDCE2(u,"Nothing to do at weight",w);
        u.FirstOfWeight[w+1] := Length(u.table)+1;
        return;
    fi;
    InfoDCE1(u,"Pushing at weight",w,"\n#I     ",
            Number(u.table,d->not d.deleted),"double",u.degree,"single",
            u.blanks,"blanks");
    flag := true;
    while flag do
        flag := false;
        for r in u.relators do
            if w >= r.weight then
                d := r.next;
                top := u.FirstOfWeight[w-r.weight+1];
                while d < top and 
                  (not IsBound(u.table[d]) or u.table[d].deleted) do
                    d := d+1;
                od;
                if d < top then
                    DCEDPush(u,u.table[d],r);
                    if u.PreEarlyClosed then
                        DCEHandlePEC(u);
                    fi;
                    if u.EarlyClosed then return; fi;
                    d := d+1;
                    if d < top then
                        flag := true;
                    fi;
                fi;
                r.next := d;    
            fi;
        od;
    od;
    DCEPack(u);
    u.FirstOfWeight[w+1] := Length(u.table)+1;
end;



#############################################################################
##
#F DCEClearDeds(u)
##

DCEClearDeds := function(u)
    local l,d,deds,back,col,i,j,x,xrec,deds1,flag;
    while Length(u.Deductions) > 0 do
        d := u.Deductions[Length(u.Deductions)];
        deds := d.deductions;
        flag := true;
        EC := false;
        while not d.deleted and flag do
            flag := false;
            for x in [1..u.nX] do
                xrec := u.Xrecs[x];
                deds1 := deds[x];
                if Length(deds1) <> 0 then
                    for col in [1..xrec.CohortSize] do
                        if not EC and not d.deleted and  IsBound(deds1[col]) then
                            back := deds1[col];
                            Unbind(deds1[col]);
                            DCEProcessDed(u,d,x,col,back);
                            flag := true;
                            if u.blanks = 0 and u.degree in u.strategy.EC then
                                u.EarlyClosed := true;
                               u.status := "early-closed";
                                return;
                            fi;
                        fi;
                    od;
                fi;
            od;
        od;
        if not d.deleted then
            i := PositionSorted(u.Deductions,d,DCEDedComp);
            l := Length(u.Deductions);
            if i <= l and IsIdentical(u.Deductions[i],d) then
                for j in [i+1..l] do
                    u.Deductions[j-1] := u.Deductions[j];
                od;
                Unbind(u.Deductions[l]);
            fi;
        fi;
    od;
end;    

#############################################################################
##
#F DCERunF(u) Do the double coset enumeration in u pure felsch style
##

DCERunF := function(u)
    local s,d,col,r,defd,defcol,l,flag,back,deds,defx,cs,x,cohort;
    u.status := "running";
    for s in u.subgens do
        DCEPush(u,u.table[1],u.id, s); # problem
    od;
    InfoDCE1(u,"Done subgroup generators");
    for r in u.sgrels  do
        DCEPush(u,u.table[1],u.id,r);
    od;
    InfoDCE1(u,"Done relators in subgroup");
    defd := 1;
    defx := 1;
    defcol := 1;
    while true do
        #
        # First clear outstanding deductions
        #
        DCEClearDeds(u);
        if u.EarlyClosed then
            return;
        fi;
        #
        # Now we will need to make a definition
        #
        flag := true;
        while flag and defd <= Length(u.table) do
            while defd <= Length(u.table) and 
                    (not IsBound(u.table[defd]) or u.table[defd].deleted) do
                defd := defd+1;
                defx := 1;
                defcol := 1;
            od;
            if defd > Length(u.table) then
                u.status := "closed";
                return;
            fi;
            d := u.table[defd];
            while flag and defx <= u.nX do
                cs := u.Xrecs[defx].CohortSize;
                cohort := d.cohorts[defx];
                while defcol <= cs and IsBound(cohort[defcol]) do
                    defcol := defcol+1;
                od;
                if defcol >cs then
                    defx := defx + 1;
                    defcol := 1;
                else
                    flag := false;
                fi;
            od;
            if defx > u.nX then
                defd := defd+1;
                defx := 1;
                defcol := 1;
            else
                flag := false;
            fi;
        od;
        if flag then
            u.status := "closed";
            return;
        fi;
        InfoDCE1(u,"Defining: ",u.nextD,"= row",d.name,u.Xrecs[defx].name,"column",defcol);
        DCEDefine(u,[d,u.Xrecs[defx].Reps[defcol]],defx);
        DCEDoStacks(u);
        DCEPack(u);
        defcol := defcol+1;
    od;
end;

#############################################################################
##
#F DCERunHavas(u) Do the double coset enumeration in hybrid style
##

DCERunHavas := function(u)
    local s,d,col,r,defd,defcol,l,flag,back,nfull,hrow,defx,cs,cohort,i;
    u.status := "running";
    for s in u.subgens do
        DCEPush(u,u.table[1],u.id, s);
    od;
    InfoDCE1(u,"Done subgroup generators");
    for r in u.sgrels  do
        DCEPush(u,u.table[1],u.id,r);
    od;
    InfoDCE1(u,"Done relators in subgroup");
    defd := 1;
    defx := 1;
    defcol := 1;
    nfull := 0;
    hrow := 1;
    while true do
        DCEClearDeds(u);
        #
        # Now we will need to make a definition
        #
        # Make a Felsch Style definition
        #   see what the next definition to make in the usual
        #  sequence is
        #
        flag := true;
        while flag and defd <= Length(u.table) do
            while defd <= Length(u.table) and 
              (not IsBound(u.table[defd]) or u.table[defd].deleted) do
                defd := defd+1;
                defx := 1;
                defcol := 1;
            od;
            if defd > Length(u.table) then
                u.status := "closed";
                return;
            fi;
            d := u.table[defd];
            while flag and defx <= u.nX do
                cs := u.Xrecs[defx].CohortSize;
                cohort := d.cohorts[defx];
                while defcol <= cs and IsBound(cohort[defcol]) do
                    defcol := defcol+1;
                od;
                if defcol >cs then
                    defx := defx + 1;
                    defcol := 1;
                else
                    flag := false;
                fi;
            od;
            if defx > u.nX then
                defd := defd+1;
                defx := 1;
                defcol := 1;
                nfull := nfull+1;
            else
                flag := false;
            fi;
        od;
        if flag then
            u.status := "closed";
            return;
        fi;
        
        #
        # Now decide whether we actually go with this or whether we
        # try something from our preferred definition list
        # 
        #
        
        if nfull < u.strategy.HavN then
            if u.nPDLs = 0 or defd*u.strategy.FF
               <= Length(u.table) then
                InfoDCE1(u,"Defining",u.nextD,"= row",d.name,
                        u.Xrecs[defx].name,"column",defcol);
                DCEDefine(u,[d,u.Xrecs[defx].Reps[defcol]],defx);
                defcol := defcol+1;
            else
                s := (u.LastPDL - u.nPDLs) mod 20 +1;
                l := u.PDL[s];
                Unbind(u.PDL[s]);
                u.nPDLs := u.nPDLs-1;
                InfoDCE1(u,"Defining",u.nextD," =",DCEDKToString(l[1]),l[2],"From PDL");
                DCEDefine(u,l[1],l[2]);
            fi;
            DCEDoStacks(u);
            DCEPack(u);
        else
            for i in [1..u.strategy.HavK] do
                while hrow < Length(u.table) and 
                  (not IsBound(u.table[hrow]) or u.table[hrow].deleted) do
                    hrow := hrow+1;
                od;
                if hrow > Length(u.table) then
                    u.status := "closed";
                    return;
                fi;
                for r in u.relators do
                    DCEDPush(u,u.table[hrow],r);
                    if u.EarlyClosed then return; fi;
                od;
                hrow := hrow+1;
                nfull := 0;
            od;
        fi;
    od;
end;


#############################################################################
##
#F DCERun(u)    Do the double coset enumeration set up in u in HLT style
##

DCERun := function(u)
    local w,s,r;
    u.status := "running";
    for s in u.subgens do
        DCEPush(u,u.table[1],u.id, s);
    od;
    InfoDCE1(u,"Done subgroup generators");
    for r in u.sgrels do
        DCEPush(u,u.table[1],u.id, r);
    od;
    InfoDCE1(u,"Also done relators in subgroup");
    DCEPack(u);
    u.FirstOfWeight := [1,Length(u.table)+1];
    w := 2;
    while not ForAll(u.relators,r->r.next > Length(u.table)) and 
      not u.EarlyClosed do
        DCEPushAll(u,w);
        w := w+1;
    od;
    if u.EarlyClosed then
        u.status := "early-closed";
    else
        u.status := "closed";
    fi;
end;
    

#############################################################################
##
#F DCESetupCols(u)  This sets up the column structure,
##
##  It is assumed that DCESetupGens has already been called on u so
##  a certain amount is already in u.


DCESetupCols := function(u)
    local o,x,xrec,gg,i;
    for x in [1..u.nX] do # Go through assigning columns
        xrec := u.Xrecs[x];
        gg := u.GGrecs[xrec.wgg];
        xrec.CohortSize := Length(gg.dom);
        xrec.Reps := [];
        for i in [1..Length(gg.dom)] do
            Add(xrec.Reps,RepresentativeOperation(u.K,gg.dom[i],xrec.pt,gg.op));
        od;
    od;
    u.nCols := Sum(u.Xrecs,xr->xr.CohortSize);
    if not IsBound(u.K.name) then
        u.K.name := "K";
    fi;
    u.id := u.K.identity;
    u.triv := DCESubgroup(u.K,[]);
    u.triv.name := "Trivial";
    u.table := [rec( muddle := u.triv, name := 1,
                     cohorts := List([1..u.nX],x->[]),
                     deductions := List([1..u.nX],x->[]),
                     operations := DCERowOps,
                     deleted := false)];
    u.blanks := u.nCols;
    u.dcct := 1;
    u.degree := Size(u.K);
    u.nextD := 2;
    u.EarlyClosed := false;
    u.PreEarlyClosed := false;
    if IsPermGroup(u.K) then
        if IsTrivial(u.K) then
            u.Kdeg := 0;
        else
            u.Kdeg := Maximum(List(u.K.generators,LargestMovedPointPerm));
        fi;
    fi;
    return u;
end;
    


#############################################################################
##
#F DCEProcessRelator(u,r)  Convert a relator r from input format (a list
##  of which each entry is either in K (as input) or a symbol
##  inverse of a symbol)
##  to internal format (a list of generator numbers or elemments of u.K)
##

DCEProcessRelator := function(u,r)
    local i,j,nr,rl,r1;
    
    #
    # First lets combine multiple permutations
    #
    
    if not IsList(r) then 
        r1 := [r];
    else
        r1 := ShallowCopy(r);
    fi;
    for i in [2..Length(r1)] do
        if not IsWord(r1[i])  then
            j := i-1;
            while j >= 1 and not IsBound(r1[j]) do 
                j := j-1;
            od;
            if j > 0 and not IsWord(r1[j]) then
                r1[j] := r1[j]*r1[i];
                Unbind(r1[i]);
            fi;
        fi;
    od;
    
    #
    # Now copy the list to get rid of holes and expand words
    # 
    nr := [];
    for rl in r1 do
        if IsWord(rl) and LengthWord(rl) > 1 then
            Append(nr,List([1..LengthWord(rl)],x->Subword(rl,x,x)));
        else
            Add(nr,rl);
        fi;
    od;
    if Length(nr) > 1 and not IsWord(nr[1])  and not IsWord(nr[Length(nr)]) then
        nr[1] := nr[Length(nr)]*nr[1];
        Unbind (nr[Length(nr)]);
    fi;
    
    #    
    # Now translate everything to internal format
    #

            
    return List(nr, 
                function(rl)
        local p;
        if not IsWord(rl) then
            return rl;
        else 
            p := Position(u.xnames,rl);
            if false <> p then
                return u.xnums[p];
            else
                p := Position(u.xnames,rl^-1);
                if false <> p and u.Xrecs[u.xnums[p]].invol = true then
                    return u.xnums[p];
                else
                    Error("Unhandleable stuff in relator ",rl);
                fi;
            fi;
        fi;
    end);
end;


#############################################################################
##
#F DCERelGroups(u,r) Construct the 'muddle' group of the relator and the
##  approriate list of conjugates
##

DCERelGroups := function(u,r)
    local n,rl,rgs,i;
    n := u.K;
    for rl in r do
        if IsInt(rl) then
            n := DCEXonSG(u,DCELin(u,n,rl,u.id),rl);
        else
            n := n^rl;
        fi;
    od;
    if IsTrivial(n) then
        return List(r,rl->u.triv);
    fi;
    if n = u.K then
        return List(r,rl->u.K);
    fi;
    rgs := [];
    rgs[Length(r)+1] := n;
    for i in [Length(r),Length(r)-1..1] do
        rl := r[i];
        if not IsInt(rl) then
            n := n^(rl^-1);
        else
            n := DCEXonSG(u,n,u.Xrecs[rl].inverse);
        fi;
        rgs[i] := n;
    od;
#    if rgs[1] <> rgs[Length(r)+1] then
#        Print("Relator has non-trivial action in K - K can be reduced");
#    fi;
    rgs := List(rgs,function(g) 
        if Length(g.generators) = 0 then 
            return u.triv; 
        else return g; fi; end);
    return rgs;
end;



#############################################################################
##
#F DCESetupGens(u) This is the first stage of the setup. In particular
##  it copes with generator inverses. Note that at this stage we only
##  have external notation for the elements of K. The only field of u which
## is assumed here is pres.
##

DCESetupGens := function(u)
    local i,act,xrec,irecs,irec,gg,id,j,score;
    u.Xrecs := [];
    irecs := [];
    u.nX := Length(u.pres.gens);
    u.K := u.pres.groupK;
    u.GGrecs := List(u.pres.gainGroups,ShallowCopy);
    for gg in u.GGrecs do
        if not IsRec(gg) then
            Error("bad gain group entry");
        fi;
        if not IsBound(gg.dom) then # with this field absent we have
            #  no action
            gg.dom := 0;
            gg.op := function(pt,gen) return pt; end;
            gg.stab := function(g,pt,op) return g; end;
        fi;
        gg.pt0 := gg.dom;
        if not IsBound(gg.op) then
            if IsInt(gg.dom) then
                gg.op := OnPoints;
            elif IsSet(gg.dom) then
                gg.op := OnSets;
            else
                Error("Can't guess an action");
            fi; 
        else
            if not IsFunc(gg.op) then
                Error("Operation must be a function");
            fi;
        fi;
        gg.dom := Set(Orbit(u.K,gg.dom,gg.op));
        if gg.op = OnPoints then
            IsRange(gg.dom); # If this orbit is a range lets find out
        fi;
        if not IsBound(gg.stab) then
            gg.stab := Stabilizer;
        elif not IsFunc(gg.stab) then
            Error("Stabilizer function must be a function");
        fi;
        if not IsBound(gg.L) then
            gg.L := gg.stab(u.K,gg.pt0,gg.op);
            if gg.L = u.K then
                gg.L := u.K;
            fi;
        elif not IsGroup(gg.L) or not IsSubgroup(u.K,gg.L) then
            Error("Gain group must be a subgroup of K");
        fi;
    od;
    for i in [1..u.nX] do
        xrec := ShallowCopy(u.pres.gens[i]);
        if not IsRec(xrec) or 
           not (IsBound(xrec.wgg) and xrec.wgg in [1..Length(u.GGrecs)]) or
           not (IsBound(xrec.name) and IsWord(xrec.name)) then
            Error("Bad generator record ",i);
        fi;
        gg := u.GGrecs[xrec.wgg];
        if not IsBound(xrec.invol) then
            xrec.invol := false;
        elif not IsBool(xrec.invol) then
            Error("invol must be boolean");
        fi;
        if not IsBound(xrec.ggconj) then
            xrec.ggconj := u.K.identity;
        elif not xrec.ggconj in u.K then
            Error("Gain group conjugator must lie in K");
        fi;
        if not IsBound(xrec.action) then
            xrec.action := false;
        elif not (xrec.action = false or xrec.action in u.K or 
                IsIsomorphism(xrec.action)) then
            Error("Don't understand this action");
        fi;
        if xrec.invol then
            xrec.inverse := i;
        elif IsBound(xrec.inverse) then
            xrec.inverse :=
                  PositionProperty(u.pres.gens,gen -> gen.name = xrec.inverse);
        else
            u.nX := u.nX+1;
            xrec.inverse := u.nX;
            irec := rec(name := xrec.name^-1,
                        invol := false,
                        inverse := i);
            act := xrec.action;
            if act = false then
                irec.wgg := xrec.wgg;
                irec.ggconj := xrec.ggconj;
                irec.action := false;
            elif act in u.pres.groupK then
                irec.action := act^-1;
                irec.wgg := xrec.wgg;
                irec.ggconj := xrec.ggconj*act;
            else
                Error("Can't invert this sort of action");
            fi;
            irec.pt := gg.op(gg.pt0,irec.ggconj);
            irec.L :=gg.L^irec.ggconj;
            if IsPermGroup(u.K) and irec.L <> u.K then
                irec.L := Group(irec.L);
            fi;
            Add(irecs,irec);
        fi;
        xrec.pt := gg.op(gg.pt0,xrec.ggconj);
        xrec.L := gg.L^xrec.ggconj; 
        if IsPermGroup(u.K) and xrec.L <> u.K then
            xrec.L := Group(xrec.L);
        fi;
        Add(u.Xrecs,xrec);
    od;
    Append(u.Xrecs,irecs);
    u.xnums := [1..u.nX];
    u.xnames := List(u.Xrecs,x->x.name);
    SortParallel(u.xnames,u.xnums);
    return u;
end;


#############################################################################
##
#F DCEIsFlippable(u,r,i) Tests whether the relator r is flippable about
##  the ith position
##

DCEIsFlippable := function(u,r,i)
    local b,f;
    if i = r.intlen then
        b := i+1;
    else
        b := i + r.intlen+1;
    fi;
    f := i-1;
    return ForAll([1..r.intlen],x -> r.relint3[f+x] = r.irelint3[b-x]);
end;

#############################################################################
##
#F DCESetupEDPs(u) Find all the essentially different positions where a
##    generator occurrs in a relator for Felsch Scanning
##

DCESetupEDPs := function(u)
    local x,edps,i,xi,r,edp,g,eg,j,l,h,rl,h1,hg,k,rl1,multis,multi,c,bgps,s,h,
          breps,ibreps,h2;
    u.Edps := [];
    u.Multis := [];
    for x in [1..u.nX] do
        l := u.Xrecs[x].L;
        xi := u.Xrecs[x].inverse;
        edps := [];
        multis := [];
        for r in u.relators do
            for i in [1..r.intlen/r.exponent] do
                if r.relint[i] = x then
                    edp := rec(rel := r, pos := i);
                    edp.flip := DCEIsFlippable(u,r,i);
                    g := AsSubgroup(Parent(l),r.groups[i]);
                    if not IsBound(r.domulti) then
                        r.domulti := r.intlen > 10 and Index(u.K,r.groups[i]) > 1000 and
                                     ForAll(r.relint, rl -> not IsInt(rl) or
                                            u.Xrecs[rl].action = false);
                    fi;
                    if r.domulti then
                        edp.ismulti := true;
                        multi := rec(rel := r, pos := i, gp := g, flip := edp.flip);
                        h := AsSubgroup(u.K,l);
                        c := u.id;
                        multi.fgps := [h];
                        multi.freps := [];
                        multi.ifreps := [];
                        multi.ogps := [h];
                        multi.conjs := [c];
                        multi.bgps := [];
                        multi.breps := [];
                        multi.ibreps := [];
                        j := 0;
                        while j < r.intlen and Size(h) > Size(g) do
                            rl := r.relint3[i+j];
                            if IsInt(rl) then
                                h1 := DCELin(u,h,rl,());
                                multi.freps[j+1] := 
                                  Set(List(RightCosets(h,h1),c->c.representative));
                                multi.ifreps[j+1] := List(multi.freps[j+1],x->x^-1);
                                h := h1;
                                h1 := h1^c;
                                Add(multi.ogps,h1);
                                bgps:= [];
                                breps := [];
                                ibreps := [];
                                s := r.intlen;
                                bgps[s] := h1;
                                while s > j and Size(h1) > Size(g) do
                                    rl1 := r.irelint3[i+s-1];
                                    s := s-1;
                                    if IsInt(rl1) then
                                        h2 := DCELin(u,h1,rl1,u.id);
                                        breps[s] := Set(List(RightCosets(h1,h2),
                                                            c->c.representative));
                                        ibreps[s] := List(breps[s],x->x^-1);
                                        h1 := h2;
                                    else
                                        h1 := h1^rl1;
                                    fi;
                                    bgps[s] := h1;
                                od;
                                multi.bgps[j+1] := bgps;
                                multi.breps[j+1] := breps;
                                multi.ibreps[j+1] := ibreps;
                            else
                                h := h^rl;
                                c := rl mod c;
                                Add(multi.ogps,multi.ogps[Length(multi.ogps)]);
                            fi;
                            Add(multi.conjs,c);
                            Add(multi.fgps,h);
                            j := j+1;
                        od;
                        multi.limit := j;
                        for j in [j..r.intlen-1] do
                            if not IsInt(r.relint3[i+j]) then
                                c := r.relint3[i+j] mod c;
                            fi;
                            Add(multi.conjs,c);
                        od;
                        Add(multis,multi);
                    else
                        edp.ismulti := false;
                    fi;
                    j :=1; 
                    while j <= Length(edps) and g <> edps[j][1] do
                        j := j+1;
                    od;
                    if j > Length(edps) then
                        Add(edps,[g,[edp]]);
                    else
                        Add(edps[j][2],edp);
                    fi;
                fi;
            od;
        od;
        Add(u.Multis,multis);
        Add(u.Edps,edps);
    od;
end;
                
 
#############################################################################
##
#F DCEAddExtraRels(u)
##

DCEAddExtraRels := function(u)
    local xrec,rel;
    for xrec in u.Xrecs do
        rel := DCEWord(u.K,[xrec.name,u.Xrecs[xrec.inverse].name]);
        rel.weight := 100;
        Add(u.relators,rel);
    od;
end;
      
      
#############################################################################
##
#F DCESetupRels(u) Process the relators and subgroup generators from u.pres
## ready for use.
## 
DCESetupRels := function(u)
    local r,i,ct;
    u.relators := List(u.pres.relators,
                        function(r) if IsDCEWord(r)
        then return(ShallowCopy(r)); 
    else return DCEWord(u.K,r); 
    fi; end);
    if u.strategy.whichStrategy = "HLT" then
        DCEAddExtraRels(u);
    fi;
    u.sgrels := [];
    ct := 1;
    for r in u.relators do
        r.num := ct;
        ct := ct+1;
        if IsBound(r.insg) and r.insg = true then
            Add(u.sgrels,r);
        fi;
        r.relintbase := DCEProcessRelator(u,r.relext);
        if IsBound(r.exponent) then
            r.relint := [];
            for i in [1..r.exponent] do
                Append(r.relint,r.relintbase);
            od;
        else
            r.relint := r.relintbase;
            r.exponent := 1;
        fi;
        Unbind(r.relintbase);
        r.intlen := Length(r.relint);
        r.groups := DCERelGroups(u,r.relint);
        r.group := r.groups[1];
        r.irelint := List(r.relint,r->DCERin(u,r));
        if u.strategy.Felsch then
            r.relint3 := Concatenation(r.relint,r.relint);
            r.irelint3 := Concatenation(r.irelint,r.irelint);
        fi;
        r.next := 1;
        if not IsBound(r.weight) then
            r.weight := LogInt(r.intlen,2);
        fi;
    od;
    u.subgens := List(u.pres.subgens,function(r) if IsDCEWord(r)
        then return(ShallowCopy(r)); else return DCEWord(u.K,r); fi; end);
        for r in u.subgens do
        r.relintbase := DCEProcessRelator(u,r.relext);
        if IsBound(r.exponent) then
            r.relint := [];
            for i in [1..r.exponent] do
                Append(r.relint,r.relintbase);
            od;
        else
            r.relint := r.relintbase;
            r.exponent := 1;
        fi;
        Unbind(r.relintbase);
        r.intlen := Length(r.relint);
        r.irelint := List(r.relint,rl->DCERin(u,rl));
    od;
    u.nRels := ct;
    if u.strategy.Felsch then
        DCESetupEDPs(u);
    fi;
    
end;



############################################################################
##
#V  DCEDefaultStrategy
##

DCEDefaultStrategy := rec(
                          whichStrategy := "HLT",
                          HavN := 10,
                          HavK := 1,
                          FF := 5,
                         EC := []);

#############################################################################
##
#F DCESetupStrat(u) Copy the strategy from u.pres and fill in defaults
##

DCESetupStrat := function(u)
    if IsBound(u.pres.strategy) then
        if not IsRec(u.pres.strategy) then
            Error("strategy must be a record");
        fi;
        u.strategy := MergedRecord(MergedRecord(DCEDefaultStrategy,u.pres.strategy),
                              u.pres.strategy);
    else
        u.strategy := DCEDefaultStrategy;
    fi;
    u.strategy.Felsch := u.strategy.whichStrategy <> "HLT";
end;
        

#############################################################################
##
#F DCESetup(pres)   Split off from DCE. Does everything except the calculation
##
            
DCESetup := function(pres)
    local u,r,rep,a;
    if not IsRec(pres) or
       not (IsBound(pres.groupK) and IsGroup(pres.groupK)) or
       not (IsBound(pres.gainGroups) and IsList(pres.gainGroups)) or
       not (IsBound(pres.gens) and IsList(pres.gens)) or
       not (IsBound(pres.relators) and IsList(pres.relators)) or
       not (IsBound(pres.subgens) and IsList(pres.subgens)) then
        Error("Bad presentation, not a record or missing compulsary components");
    fi;
    u := rec(pres := pres);
    u.operations := DCEUniverseOps;
    u.starttime := Runtime();
    u.isDCEUniverse := true;
    if IsBound(pres.name) then
        u.name := pres.name;
    else
        u.name := "No name";
    fi;
    u.status := "Setting up";
    DCESetupStrat(u);
    DCESetupGens(u);
    InfoDCE1(u,"Set up generators and inverses");
    DCESetupCols(u);
    InfoDCE1(u,"Set up column structure:",u.nCols,"columns");
    DCEInitStacks(u);
    DCESetupRels(u);
    InfoDCE1(u,"Pre-processed relators");
    u.status := "Set up";
    return u;
end;

#############################################################################
##
#F DCE(pres)  sets up and runs the Double Coset Enumeration described in pres
##            and returns the completed universe
##


DCE := function(pres)
    local u;
    u := DCESetup(pres);
    if u.strategy.whichStrategy = "Havas" then
        DCERunHavas(u);
    elif u.strategy.whichStrategy = "Felsch" then
        DCERunF(u);
    elif  u.strategy.whichStrategy = "HLT" then
        DCERun(u);
    else
        Error("Unrecognised strategy ",u.strategy.whichStrategy);
    fi;
    DCESuperPack(u);
    return u;
end;

                      


#############################################################################
##
#F DCECheckDK(u,dk)  check validity of a dk pair
##

DCECheckDK := function(u,dk)
    return  IsList(dk) and
            Length(dk) = 2 and
            ForAny(u.table, x-> IsIdentical(dk[1],x)) and
            dk[2] in u.K;
end;
 

#############################################################################
##
#F DCECheckTab  - check the table for internal consistency
##

DCECheckTab := function(u)
    local d,row,col,beta,col1,dk,xrec,cohort,x;
    if not IsDCEUniverse(u) or
       not IsBound(u.table) or
       not IsList(u.table) or
       Length(u.table) = 0 then
        Print("Badly structured, can't start\n");
        return false;
    fi;
    row := 0;
    for d in u.table
      do
        row := row+1;
        if not IsRec(d) or not IsBound(d.operations) or
           not (IsBound(d.deleted) and IsBound(d.name)) then
            Print("Badly structured row " ,row,"\n");
            return false;
        fi;
        Print("Checking ",row," name ",d.name,"\n");
        if d.deleted then
            if not IsBound(d.replace) or not DCECheckDK(u,d.replace) then
                Error("Bad replacement pair in row ",d.name,"\n");
                return false;
            fi;
        else
            if not (IsBound(d.muddle) and IsBound(d.cohorts)) then
                Error("Mandatory components missing",d.name);
                return false;
            fi;
            if not IsSubgroup(u.K,d.muddle) then
                Error("Bad muddle group in row ",d.name,"\n");
                return false;
            fi;
            for x in [1..u.nX] do
                xrec := u.Xrecs[x];
                cohort := d.cohorts[x];
                for col in [1..xrec.CohortSize] do
                    if IsBound(cohort[col]) then
                        if IsRec(cohort[col]) then
                            beta := cohort[col];
                            if not beta.elt in d.muddle or 
                               not col in [1..xrec.CohortSize] then
                                Error("Bad beta ",d.name," ",x," ",col," ",beta);
                                return false;
                            fi;
                            if beta.col <>  
                               DCEColumnRegular(u,beta.elt*xrec.Reps[col],x) then
                                Error("Two parts of beta entry do not match ",
                                      d.name," ",x," ",col);
                                return false;
                            fi;
                            if IsBound(cohort[beta.col]) and IsRec(cohort[beta.col]) then
                                Error("Beta entry points to beta entry at row ",
                                      d.name," col ",col,"\n");
                                return false; 
                            fi;
                        else   
                            if not DCECheckDK(u,cohort[col]) then
                                Error("Bad alpha entry at row ",
                                      d.name," ",x," ",col,"\n");
                                return false;
                            fi;
                            dk := DCEFollow(u,ShallowCopy(cohort[col]), xrec.inverse);
                            if not IsIdentical(dk[1],d) or 
                               not dk[2]/xrec.Reps[col] in d.muddle then
                                Error("Inverse failure row ",d.name," ",x," ",col,"\n");
                                return false;
                            fi;
                        fi;
                    fi;
                od;
            od;
        fi;
    od;
    return true;
end;
                           
#############################################################################
##
#F DCEWrite(u,fn)  Write the information from the dct that is computed
##                 by the program out to file "fn" in s format to be read back
##                 by DCERead. This can only legitimately be called when
##                 the stacks are empty
##

DCEWrite := function(u,fn)
    local i,d,oldKname,oldTrivname,l;
    oldKname := u.K.name;
    Unbind(u.K.name);
    oldTrivname := u.triv.name;
    Unbind(u.triv.name);
    PrintTo(fn,"DCEReadin.table := [\n");
    for d in u.table do
        if not d.deleted then
            l := [d.name,d.muddle,
                  List(d.cohorts, c -> List(c,
                          function(e)
                          if IsList(e) then return [e[1].name,e[2]];  
                      else return e; fi; end))];
            AppendTo(fn,l,",\n");
        fi;
    od;
    AppendTo(fn,"];\n");
    AppendTo(fn,"DCEReadin.status :=\"",u.status,"\";\n");
    u.K.name := oldKname;
    u.triv.name := oldTrivname;
end;
                            
   
#############################################################################
##
#F DCERead(u,fn)  Read the information from the dct that is computed
##                 by the program back from file "fn" in the format written
##                 by DCEWrite. 

DCEReadin := rec();

DCERead := function(u,fn)
    local row,d,mud,c,e,names;
    Read(fn);
    InfoDCE1(u,"Read the file");
    u.table := [];
    names := [];
    for row in DCEReadin.table do
        d := rec(name := row[1]);
        mud := row[2];
        mud := DCESubgroup(u.K,mud.generators);
        if Length(mud.generators) = 0 then
            mud := u.triv;
        fi;
        d.muddle := mud;
        d.cohorts := row[3];
        Unbind(row[3]);
        d.deleted := false;
        Add(u.table,d);
        Add(names,d.name);
        IsSet(names);
        InfoDCE2(u,d.name);
    od;
    Unbind(DCEReadin.table);
    for d in u.table do
        for c in d.cohorts do
            for e in c do
                if IsList(e) then
                    e[1] := u.table[PositionSorted(names,e[1])];
                fi;
            od;
        od;
    od;
    u.dcct := Length(u.table);
    u.degree := Sum(u.table,d->Index(u.K,d.muddle));
    if IsBound(DCEReadin.status) then
        u.status := DCEReadin.status;
        Unbind(DCEReadin.status);
    else
        u.status := "closed";
    fi;
end;
                            
   
#
#
# Collapsed Adjancency stuff for computing on a 
# closed DCT
#
#

#############################################################################
##
#F DCEWordGroup(u,w) w a word in internal form - returns the 'muddle group
## appropriate for computing H-orbits and collapsed adjacency
##

DCEWordGroup := function(arg)
    local n,rl,rgs,i,n0,u,w;
    if Length(arg) < 2 or Length(arg) > 3 then
        Error("Usage DCEWordGroup(u,w[,M)");
    fi;
    u := arg[1];
    w := arg[2];
    if Length(arg) = 2 then
        if IsBound(u.M) then
            n0 := u.M;
        else
            n0 := u.K;
        fi;
    else
        n0 := arg[3];
    fi;
    n := n0;
    for rl in w do
        if IsInt(rl) then
            n := DCEXonSG(u,DCELin(u,n,rl,u.id),rl);
        else
            n := n^rl;
        fi;
    od;
    n := Intersection(n,n0);
    if IsTrivial(n) or n = n0 then 
        return n;
    fi;
    for i in [Length(w),Length(w)-1..1] do
        rl := w[i];
        if not IsInt(rl) then
            n := n^(rl^-1);
        else
            n := DCEXonSG(u,n,u.Xrecs[rl].inverse);
        fi;
    od;
    return n;
end;


#############################################################################
##
#F DCEInitMBetas -- sets up the data structures for recording which M cosets
##                  lie in what (M_d,M) double cosets and which H orbits
##                   those belong to.
##


DCEInitMBetas := function(u,d,allones)
    local mud,orbs,orb,rep,pt,s,mudim;
    mud := d.muddle;
    d.horbs := [];
    if mud = u.triv then
        d.hsizes := allones;
        return;
    fi;
    if mud = u.K then
        d.hsizes := [1];
        for pt in [2..Length(allones)] do
            d.horbs[pt] := -1;
        od;
        return;
    fi;
    d.hsizes := [];
    mudim := Image(u.Mhom,mud);
    s := Size(mud);
    orbs := Orbits(mudim,[1..u.Mindex]);
    for orb in orbs do
        rep := orb[1];
        d.hsizes[rep] := Length(orb)*Size(u.M)/s;
        for pt in orb do
            if pt <> rep then
                d.horbs[pt] := -rep;
            fi;
        od;
    od;
    return;
end;

    
#############################################################################
##
#F DCEHOrbitsRegular(u)   Add the H orbit information to the table, returns a
##  list of reps - Assumes that K <= H!!
##

DCEHOrbitsRegular := function(u)
    local gens,h,orbreps,orbdcts,orbsizes,d,q,ct,size,d1,d2,
          dcr,k,toadd,dk,orb;
    if IsBound(u.orbreps) then
        return;
    fi;
    gens := Filtered(u.subgens, r -> not ForAll(r.relint,l->l in u.K));
    for h in gens do
        if h.intlen > 1 then
            h.group := DCEWordGroup(u,h.relint);
        fi;
    od;
    orbreps := [];
    orbdcts := [];
    orbsizes := [];
    orb := 1;
    for d in u.table do
        if not (d.deleted or IsBound(d.orbit)) then
            q := [d];
            ct := 1;
            size := Index(u.K,d.muddle);
            Add(orbreps,[d,1]);
            d.orbit := orb;
            for d1 in q do
                for h in gens do
                    if h.intlen = 1 then
                        toadd := List(Filtered(d1.cohorts[h.relint[1]],IsList),dk->dk[1]);
                    else
                        dcr := DCEDoubleCosetReps(u.K,d1.muddle,h.group);
                        toadd := [];
                        for k in dcr do
                            dk := [d1,k];
                            DCEMultiFollow(u,dk,h.relint,1,h.intlen+1,false);
                            Add(toadd,dk[1]);
                        od;
                    fi;
                    for d2 in toadd do
                        if not IsBound(d2.orbit) then
                            d2.orbit := orb;
                            ct := ct+1;
                            size := size+Index(u.K,d2.muddle);
                            Add(q,d2);
                        fi;
                    od;
                od;
            od;
            Add(orbdcts,ct);
            Add(orbsizes,size);
            orb := orb+1;
        fi;
    od;
    u.orbreps := orbreps;
    u.orbdcts := orbdcts;
    u.orbsizes := orbsizes;
end;
                           

#############################################################################
##
#F DCEHOrbits -- determined the H orbits
##

DCEHOrbits := function(arg)
    local M,cos,gens,allones,orbreps,orbsizes,orb,rep,i,genims,
          d,gens,orbreps,orbsizes,orb,rep,dc,q,size,dc1,dk2,dcr,k,i2,d2,h,m1,d1,i1,
          k1,u,thresh;
    if Length(arg) = 1 then
        u := arg[1];
        thresh := 10;
    elif Length(arg) = 2 then
        u := arg[1];
        thresh := arg[2];
    fi;
    if IsBound(u.orbreps) then
        return;
    fi;
    if not (IsDCEUniverse(u) and u.status in ["closed","early-closed"]) then
        Error("The argument must be the DCE universe of a completed enumeration");
    fi;
    if not IsBound(u.M) then
        M := u.table[1].muddle; # H ^ K
        u.M := M;
        if M = u.K then
            return DCEHOrbitsRegular(u);
        fi;
        if M = u.triv then
            Print("You might as well extract permutations\n");
            return;
        fi;
      #
      # So now we have to figure out the orbits of M
      # on single cosets.
      #
        cos := Set(LeftCosets(u.K,M));
        u.Mcos := cos;
        if Size(u.K)*thresh < u.degree then
            u.mcPrecalc := true;
            u.whichMcos := [];
            for i in [1..Length(cos)] do
                InfoDCE3(u,"Coset",i);
                for k in Elements(cos[i]) do
                    u.whichMcos[Position(Elements(u.K),k)] := i;
                od;
                Unbind(cos[i].elements);
            od;
        else
            u.mcPrecalc := false;
        fi;
        u.Mreps := List(cos,c->c.representative);
        u.Mindex := Index(u.K,u.M);
        gens := u.K.generators;
        genims := List(gens,g->Permutation(g,cos,OnLeftAntiOperation));
        u.KonM := Group(genims,u.id);
        u.Mhom := GroupHomomorphismByImages(u.K,u.KonM,List(gens,x->x^-1),genims);
        if not IsHomomorphism(u.Mhom) then
            Error("Construction problems with the hom.");
        fi;
    else
        M := u.M; 
        cos := u.Mcos;
        if u.M = u.K then
            return DCEHOrbitsRegular(u);
        fi;
    fi;
    InfoDCE1(u,"Completed preliminaries, index of M is ",u.Mindex);
        allones := List(cos,x->Size(M));
    for d in u.table do
        if not d.deleted then
            DCEInitMBetas(u,d,allones);
        fi;
    od;
    InfoDCE1(u,"Annotated table");
    gens := Filtered(u.subgens, 
                    r -> not (IsBound(r.central) and  r.central) and 
                    not ForAll(r.relint,l->l in u.K));
    M := Group(M);
    for h in gens do
        h.group := AsSubgroup(M,DCEWordGroup(u,h.relint));
    od;
    orbreps := [];
    orbsizes := [];
    orb := 1;
    for d in u.table do
        if not d.deleted then
            for i in [1..u.Mindex] do
                if not IsBound(d.horbs[i]) then
                    dc := [d,i];
                    q := [dc];
                    size := d.hsizes[i];
                    Add(orbreps,dc);
                    d.horbs[i] := orb;
                    for dc1 in q do
                        d1 := dc1[1];
                        i1 := dc1[2];
                        k1 := u.Mreps[i1];
                        if d1.muddle = u.K then
                            m1 := M;
                        elif d1.muddle = u.triv then
                            m1 := TrivialSubgroup(M);
                        else
                            m1 := Intersection(d1.muddle^k1,M);
                            m1 := M.operations.AsSubgroup(M,m1);
                        fi;
                        for h in gens do
                            dcr := DCEDoubleCosetReps(M,m1,h.group);
                            for k in dcr do
                                dk2 := [d1,k1*k];
                                DCEMultiFollow(u,dk2,h.relint,1,h.intlen+1,false);
                                if u.mcPrecalc then
                                    i2 := u.whichMcos[Position(Elements(u.K),dk2[2])];
                                else
                                    i2 := 1^Image(u.Mhom,dk2[2]^-1);
                                fi;
                                d2 := dk2[1];
                                if IsBound(d2.horbs[i2]) and d2.horbs[i2] < 0 then
                                    i2 := -d2.horbs[i2];
                                fi;
                                if not IsBound(d2.horbs[i2]) then
                                    d2.horbs[i2] := orb;
                                    Add(q,[d2,i2]);
                                    size := size+d2.hsizes[i2];
                                    InfoDCE3(u,d2.name,i2,size);
                                elif d2.horbs[i2] <> orb then
                                    Error("Orbit malfunction");
                                fi;
                            od;
                        od;
                    od;
                    Add(orbsizes,size);
                    InfoDCE1(u,"Completed orbit ",orb," size ",size);
                    orb := orb+1;
                fi;
            od;
        fi;
    od;
    u.orbsizes := orbsizes;
    u.orbreps := orbreps;
end;

#############################################################################
##
#F DCEAddSchreier(u)   Add Schreier vector information to the double coset entries
##

DCEAddSchreier := function(u)
    local d,d1,q,col,x,nreps,doorbs,i,d2,hard,o;
    if IsBound(u.table[1].sv) then
        return;
    fi;
    q := [u.table[1]];
    doorbs := IsBound(u.orbreps);
    hard := IsBound(u.M) and u.M <> u.K;
    d := u.table[1];
    if doorbs then
        if hard then
            nreps := [];
            for i in [1..u.Mindex] do
                o := d.horbs[i];
                if o > 0 and not  IsBound(nreps[o]) then
                    nreps[o] := [d,i];
                fi;
            od;
        else
            nreps := [d];
        fi;
    fi;
    u.table[1].sv := 0;
    for d in q do
        for x in [1..u.nX] do
            for col in d.cohorts[x] do
                if IsList(col) and not IsBound(col[1].sv) then
                    col[1].sv := [col[2],u.Xrecs[x].inverse];
                    Add(q,col[1]);
                    if doorbs then
                        if hard then
                            d2 := col[1];
                            for i in [1..u.Mindex] do
                                o := d2.horbs[i];
                                if o > 0 and not  IsBound(nreps[o]) then
                                    nreps[o] := [d2,i];
                                fi;
                            od;
                        else
                            if not IsBound(nreps[col[1].orbit]) then
                                nreps[col[1].orbit] := col[1];
                            fi;
                        fi;
                    fi;
                fi;
            od;
        od;
    od;
    if doorbs then
        u.orbreps := nreps;
    fi;
end;




#############################################################################
##
#F DCERepWord(u,dk)
##

DCERepWord := function(u,dk)
    local w,k,dk1;
    w := [];
    dk1 := ShallowCopy(dk);
    while dk1[1].sv <> 0 do
        k := dk1[2] mod dk1[1].sv[1];
        if k <> u.id then
            Add(w,k);
            dk1[2] := dk1[1].sv[1];
        fi;
        Add(w,dk1[1].sv[2]);
        DCEFollow(u,dk1,dk1[1].sv[2]);
    od;
    if dk1[2] <> u.id then 
        Add(w,dk1[2]^-1);
    fi;
    return Reversed(List(w,r->DCERin(u,r)));
end;
    

#############################################################################
##
#F DCEColAdj(u)
##
##  
DCEColAdjContRegular := function(u,ds)
    local mat,reps,gps,sizes,dcs,d,dc,dk,i,ms,ireps,m;
    ireps := [1..Length(u.orbreps)];
    mat := List(ireps, y->List(ireps,r -> List(ireps, x->0)));
    reps := List(ireps,o->DCERepWord(u,[u.orbreps[o],u.id])); 
    gps := List(ireps,i->DCEWordGroup(u,reps[i]));
    sizes := List(gps,Size);
    InfoDCE1(u,"Pre-calculation completed, starting main loop");
    for d in ds do
        if not d.deleted then
            ms := Size(d.muddle);
            for i in ireps do
                m := mat[d.orbit][i];
                if ms = 1 then
                    dcs := DCELeftCosetReps(u.K,gps[i]);
                    for dc in dcs do
                        dk := [d,dc];
                        DCEMultiFollow(u,dk,reps[i],1,Length(reps[i])+1,false);
                        m[dk[1].orbit] := m[dk[1].orbit] + sizes[i];
                    od;
                elif sizes[i] = 1 then
                    dcs := DCERightCosetReps(u.K,d.muddle);
                    for dc in dcs do
                        dk := [d,dc];
                        DCEMultiFollow(u,dk,reps[i],1,Length(reps[i])+1,false);
                        m[dk[1].orbit] := m[dk[1].orbit] + 1;
                    od;
                else
                    dcs := DoubleCosets(u.K,d.muddle,gps[i]);
                    for dc in dcs do
                        dk := [d,dc.representative];
                        DCEMultiFollow(u,dk,reps[i],1,Length(reps[i])+1,false);
                        m[dk[1].orbit] := m[dk[1].orbit] + Size(dc)/ms;
                    od;
                fi;
            od;
        fi;
        InfoDCE1(u,"Added contribution from",d.name);
    od;
    return mat;
end;


DCEColAdjCont := function(u,dis)
    local mat,reps,gps,sizes,dcs,d,dc,dk,i,ms,ireps,m,M,di,i1,
          k,d2,i2,m1,k,d2ho,fact;
    DCEHOrbits(u);
    DCEAddSchreier(u);
    if not IsBound(u.M) or u.M = u.K then
        return DCEColAdjContRegular(u,List(dis,di->di[1]));
    fi;
    ireps := [1..Length(u.orbreps)];
    mat := List(ireps,x->0);
    mat := List(ireps,x->ShallowCopy(mat));
    mat := List(ireps,x->Copy(mat));
    reps := List(u.orbreps,o->DCERepWord(u,[o[1],u.Mreps[o[2]]])); 
    M := Group(u.M);
    gps := List(reps,r->AsSubgroup(M,DCEWordGroup(u,r)));
    sizes := List(gps,Size);

    for di in dis do
        d := di[1];
        i1 := di[2];
        k := u.Mreps[i1];
        if d.muddle = u.K then
            m1 := M;
        elif d.muddle = u.triv then
            m1 := TrivialSubgroup(M);
        else
            m1 := Intersection(d.muddle^k,u.M);
            m1 := M.operations.AsSubgroup(M,m1);
        fi;
        ms := Size(m1);
        for i in ireps do
            m := mat[d.horbs[i1]][i];
            if ms = 1 then
                dcs := DCELeftCosetReps(M,gps[i]);
                for dc in dcs do
                    dk := [d,k*dc];
                    DCEMultiFollow(u,dk,reps[i],1,Length(reps[i])+1,false);
                    if IsBound(u.whichMcos) then
                        i2 := u.whichMcos[PositionSorted(Elements(u.K),dk[2])];
                    else
                        i2 := 1^Image(u.Mhom,dk[2]^-1);
                    fi;
                    d2ho := dk[1].horbs;
                    if d2ho[i2] < 0 then
                        i2 := -d2ho[i2];
                    fi;
                    m[d2ho[i2]] := m[d2ho[i2]] + sizes[i];
                    InfoDCE4(u,d.horbs[i1],d2ho[i2],sizes[i],"a");
                od;
            elif sizes[i] = 1 then
                dcs := DCERightCosetReps(M,m1);
                for dc in dcs do
                    dk := [d,k*dc];
                    DCEMultiFollow(u,dk,reps[i],1,Length(reps[i])+1,false);
                    if IsBound(u.whichMcos) then
                        i2 := u.whichMcos[PositionSorted(Elements(u.K),dk[2])];
                    else
                        i2 := 1^Image(u.Mhom,dk[2]^-1);
                    fi;
                    d2ho := dk[1].horbs;
                    if d2ho[i2] < 0 then
                        i2 := -d2ho[i2];
                    fi;
                    m[d2ho[i2]] := m[d2ho[i2]] +  1;
                    InfoDCE4(u,d.horbs[i1],d2ho[i2],1,"b");
                od;
            else
                dcs := DoubleCosets(M,m1,gps[i]);
                for dc in dcs do
                    dk := [d,k*dc.representative];
                    DCEMultiFollow(u,dk,reps[i],1,Length(reps[i])+1,false);
                    if IsBound(u.whichMcos) then
                        i2 := u.whichMcos[PositionSorted(Elements(u.K),dk[2])];
                    else
                        i2 := 1^Image(u.Mhom,dk[2]^-1);
                    fi;
                    d2ho := dk[1].horbs;
                    if d2ho[i2] < 0 then
                        i2 := -d2ho[i2];
                    fi;
                    m[d2ho[i2]] := m[d2ho[i2]] + Size(dc)/ms;
                    InfoDCE4(u,d.horbs[i1],d2.horbs[i2],Size(dc)/ms,i,sizes[i],ms,d.name,i1,"d");
                od;
            fi;
        od;
        InfoDCE1(u,"Added contribution from ", d.name,"part",i1);
    od;
    return mat;
end;
    
    

    
DCEColAdj := function(u)
    local ix;
    if not (IsDCEUniverse(u) and u.status in ["closed","early-closed"]) then
        Error("Argument must be an [early-]closed DCE Universe");
    fi;
    if IsBound(u.Mindex) then
        ix := [1..u.Mindex];
        return DCEColAdjCont(u,
                       Concatenation(List(u.table,d->List(Filtered(ix,x->d.horbs[x]>0),x->[d,x]))));
    else
        return DCEColAdjCont(u,List(u.table,d->[d,1]));
    fi;
end;
    


#############################################################################
##
#F DCEColAdjSingle(u,o)
##

DCEColAdjSingle := function(u,o)
    local m,ix;
    if not (IsDCEUniverse(u) and u.status in ["closed","early-closed"]) then
        Error("Argument must be an [early-]closed DCE Universe");
    fi;
    if not (IsInt(o) and o > 0 and o < u.degree) then
        Error("Second argument cannot possibly be an orbit number");
    fi;
    if not IsBound(u.orbsizes) then
        DCEHOrbits(u);
    fi;
    if o > Length(u.orbsizes) then
        Error("Second argument is larger than the number of orbits");
    fi;
    if IsBound(u.Mindex) then
        ix := [1..u.Mindex];
        m := DCEColAdjCont(u,Concatenation(List(u.table,d->List(Filtered(ix,x->d.horbs[x]=o),x->[d,x]))));
    else
        m := DCEColAdjCont(u,List(Filtered(u.table,d->d.orbit = o),d->[d,1]));
    fi;
    return m[o];
end;
                           

#############################################################################
##
#F SetupSymmetricPresentation(K,name[,base[,op]])  
##    setup to use DCE on a symmetric presentation in the style of Curtis
##    K is the group permuting the symmetric generators
##    name is an AbstractGenerator which names the generators
##    base is the point in the orbit of K corresponding to the generator (rather
##         than to some conjugate) (default 1)
##    op   is the operation by which K acts on the generators
##
##    The value returned is a record with two fields
##     skeleton is a partial DCE presentation. The fields K, gainGroups
##              and gens are present, relators and subgens must still be added
##     makeGen is a function which converts points in Orbit(k,base,op) into
##              DCEWords (conjugates of DCEWord(name)).
##

SetupSymmetricPresentation := function(arg)
    local K,genname,basept,op,usage,skel,func;
    usage := "Usage SymGenFunc(K,genname[,basept[,op]])";
    if Length(arg) < 2 then
        Error(usage);
    fi;
    K := arg[1]; 
    genname := arg[2];
    if not (IsGroup(K) and IsWord(genname) and LengthWord(genname) = 1) then
        Error(usage);
    fi;
    if Length(arg) = 2 then
        basept := 1;
        op := OnPoints;
    elif Length(arg) = 3 then
        basept := arg[3];
        op := OnPoints;
    else
        basept := arg[3];
        op := arg[4];
    fi;
    skel := rec(groupK := K,
                gainGroups := [rec(dom := basept, 
                        op := op)],
                gens := [rec(name := genname, wgg := 1)] ); 
    func := function(pt)
        local rep;
        if pt = basept then
            return DCEWord(K,genname);
        fi;
        rep := RepresentativeOperation(K,basept,pt,op);
        if rep = false then
            Error("point not in orbit");
        fi;
        return DCEWord(K,[rep^-1,genname,rep]);
    end;
    return rec(skeleton := skel,makeGen := func);
end;
    


#############################################################################
##
#F DCEEqual(u,dk1,dk2)   Test if two names refer to the same single coset
##                       

DCEEqual := function(u,dk1,dk2)
    return IsIdentical(dk1[1],dk2[1]) and dk1[2]/dk2[2] in dk1[1];
end;


##M###############################################################################
#E  Emacs . . . . . . . . . . . . . . . . . . . . . . . local emacs variables##
##  Local Variables:
##  mode:               gap
##  fill-column:        73
##  fill-prefix:        "##  "
##  tab-width:	 			3
##  tab-stop-list:		(3 6 9 12 15 18 21 24 27 30 33 36 39 42 45 48 51 54 57 60 63 66 69 72 75 78)
##  register-alist:		((118 . "\n############################################################################\n##\n#V  Variable declaration\n##\n") (102 . "\n#############################################################################\n##\n#F Function Header \n##\n"))
##  End:
##

