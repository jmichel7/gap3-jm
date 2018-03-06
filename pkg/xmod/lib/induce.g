###############################################################################
##
##  induce.g                    for GAP 3.4                    version  9/ 1/97
##
###############################################################################
##
#A  induce.g                    GAP library                       Chris Wensley
##                                                                    Murat Alp
#Y  Copyright
##
##  This file contains functions that manipulate crossed modules,
##  pre-crossed modules, group graphs and associated constructions.
##
#H  $Log: induce.g,v $
#H  Revision 1.1  1997/03/27 13:35:50  gap
#H      Added xmod 1.31
#H
#H  	    SL
#H
#H

###############################################################################
##  Layout of groups for induced crossed module calculations:
## 
##          Fin. Pres. Groups               Perm. Groups
##          =================               ============
##
##          N=M -------> I              RN=RM -------> RI
##            |          |                  |          |
##            |          |                  |          |
##            V          V                  V          V
##            P -------> Q                 RP -------> RQ
##
##############################################################################
##  Layout of groups for induced cat1-group calculations:
##
##                                       iotastar
##         G ------------> I       PG ------------> PI
##         ||             ||       ||               || 
##         ||             ||       ||               ||  
##         th            t*h*      th              t*h* 
##         ||             ||       ||               || 
##         ||             ||       ||               || 
##         VV             VV       VV               VV 
##         R ------------> Q       PR ------------> PQ
##                                        iota
##############################################################################


##############################################################################
##
#F  InclusionInducedXModData     set up data for function InclusionInducedXMod
##

InclusionInducedXModData := function( arg )

    local  Qdata, Pdata, Mdata, Tdata, IXMdata, nargs, X, iota, imgeniota,
           T, t, qP, qT, pos, posp2f, posf2p, i, diff, ok,
           Q, oQ, q, Qrec, RQ, elRQ, genRQ, degRQ, indQP,
           P, oP, p, Prec, RP, genRP,
           iRP, eliRP, iRN, geniRN,
           M, oM, m, Mrec, RM, elRM, genRM, ngRM, kept,
           genM, ngM, Mgen, relM, nrM, rm,
           N, Nrec, RN, genRN, ngRN, genN, ngN, l1, l2, m1, m2, elM, presM,
           fN, genfN, subfN, defN, relN, nrN, n1, n2, posn, posm,
           presM, f2pM, f2pN;

    nargs := Length( arg );
    if ( nargs <> 3 ) then
        Print( "\nUsage: InclusionInducedXModData( X, iota, T );\n" );
        return false;
    fi;
    X := arg[1];
    iota := arg[2];
    Mdata := X.source;
    Pdata := X.range;
    Qdata := iota.range;
    Tdata := arg[3];

    IXMdata := rec( );
    
    # process Qdata
    if ( IsBound( Qdata.isFpPair ) and Qdata.isFpPair ) then
        Qrec := Qdata;
    elif ( IsBound( Qdata.isPermGroup ) and Qdata.isPermGroup ) then
        Qrec := rec( );
        Qrec.perm := Qdata;
    else
        Error( " Qdata not of correct type" );
    fi;
    RQ := Qrec.perm;
    oQ := Size( RQ );
    elRQ := Elements( RQ );
    genRQ := RQ.generators;
    degRQ := PermGroupOps.NrMovedPoints( RQ );
    
    # process Pdata
    if ( IsBound( Pdata.isFpPair ) and Pdata.isFpPair ) then
        Prec := Pdata;
    elif ( IsBound( Pdata.isPermGroup ) and Pdata.isPermGroup ) then
        Prec := rec( );
        Prec.perm := Pdata;
    else
        Error( "Pdata not of correct type" );
    fi;
    RP := Prec.perm;
    genRP := RP.generators;
    if ( Parent( RQ ) = Parent( RP ) and IsSubgroup( RQ, RP ) ) then
        iRP := RP;
    else
        imgeniota := List( genRP, p -> Image( iota, p ) );
        iRP := Subgroup( RQ, imgeniota );
    fi;
    Prec.image := iRP;
    oP := Size( RP );
    eliRP := Elements( iRP );
    indQP := oQ/oP;
    
    # process Mdata
    if ( IsBound( Mdata.isFpPair ) and Mdata.isFpPair ) then
        Mrec := Mdata;
    elif ( IsBound( Mdata.isPermGroup ) and Mdata.isPermGroup ) then
        Mrec := rec( );
        Mrec.perm := Mdata;
    else
        Error( "Mdata not of correct type" );
    fi;
    RM := Mrec.perm;
    if not IsNormal( RP, RM ) then
        Error( "RM not a normal subgroup of RP" );
    fi;
    oM := Size( RM );
    elRM := Elements( RM );
    genRM := RM.generators;
    ngRM := Length( genRM );
    
    # find a presentation for RM
    if ( IsBound( RM.fpPair ) and IsFpPair( RM.fpPair ) ) then
        Mrec := RM.fpPair;
        M := Mrec.fp;
        genM := M.generators;
        ngM := Length( genM );
        f2pM := Mrec.f2p;
    else
        presM := PresentationViaCosetTable( RM );
        TzInitGeneratorImages( presM );
        M := FpGroupPresentation( presM );
        if ( ngRM <> Length( presM.generators ) ) then
            Print( "\n\nWARNING: extra generators defined - \n" );
            Print( "             the isomorphism f2pM may be incorrect!\n\n" );
        fi;
        Mrec := rec( );
        Mrec.perm := RM;
        Mrec.fp := M;
        genM := M.generators;
        ngM := Length( genM );
        f2pM := GroupHomomorphismByImages( M, RM, genM, genRM );
        if not IsIsomorphism( f2pM ) then
            Error( "f2pM not an isomorphism" );
        fi;
        Mrec.f2p := f2pM;
        Mrec.isFpPair := IsFpPair( Mrec );
    fi;
    
    presM := PresentationFpGroup( M );
    TzInitGeneratorImages( presM );
    if ( XModPrintLevel > 2 ) then
        Print( "initial presentation for M :-\n\n" );
        TzPrint( presM );
    fi;
    relM := M.relators;
    nrM := Length( relM );
    elM := Elements( M );   

    # Determine the positions of isomorphic images M <-> RM
    posf2p := 0 * [1..oM];
    posp2f := 0 * [1..oM];
    for i in [1..oM] do
        m := elM[i];
        rm := Image( f2pM, m );
        pos := Position( elRM, rm );
        posf2p[i] := pos;
        posp2f[pos] := i;
    od;
    if ( XModPrintLevel > 1 ) then
        Print( "posf2p = ", posf2p, "\n" );
        Print( "posp2f = ", posp2f, "\n" );
    fi;

    # determine  genRN = closure of  genRM  under conjugation by  RP
    genRN := Copy( genRM );
    Mgen := Copy( genM );
    ngN := 0;
    ngRN := Length( genRN );
    if ( XModPrintLevel > 1 ) then
        Print( "checking that M.generators is closed under P-action\n" );
    fi;
    while ngN < ngRN do
        ngN := ngN + 1;
        n1 := genRN[ngN];
        for p in genRP do
            n2 := n1^p;
            posn := Position( genRN, n2 );
            if ( posn = false ) then
                Add( genRN, n2 );
                posm := Position( elRM, n2 );
                if ( posm = false ) then
                    Error( "M is not a normal subgroup of P" );
                else
                    m2 := elM[ posp2f[ posm ] ];
                    Add( Mgen, m2 );
                    ngRN := ngRN + 1;
                fi;    
            fi;
        od;
    od;
    if ( XModPrintLevel > 1 ) then
        Print( "Closure of RM.generators under action of P :-\n" );
        Print( genRN, "\n\n" );
    fi;
    # prepare copy N of M with more generators
    fN := FreeGroup( ngRN, "fN" );
    genfN := fN.generators;
    subfN := Sublist( genfN, [1..ngM] );
    defN := List( Mgen, g -> MappedWord( g, genM, subfN ) );
    relN := List( relM, r -> MappedWord( r, genM, subfN ) );
    for i in [ (ngM+1)..ngRN ] do
        n1 := defN[i];
        n2 := genfN[i]^(-1);
        Add( relN, n1*n2 );
    od;
    if ( XModPrintLevel > 1 ) then
        Print( "Extended set of relators for N :- \n", relN, "\n\n" );
    fi;
    nrN := Length( relN );
    
    genN :=  [ ];
    for n1 in genRN do
        n2 := PreImagesRepresentative( f2pM, n1 );
        posn := Position( elRM, n1 );
        posm := posp2f[ posn ];
        m1 := elM[ posm ];
        l1 := LengthWord( m1 );
        l2 := LengthWord( n2 );
        if ( l2 < l1 ) then
            Add( genN, n2 );
        else
            Add( genN, m1 );
        fi;
    od;
    diff := ngRN - ngM;
    RN := Subgroup( RP, genRN );
    geniRN := List( genRN, r -> Image( iota, r ) );
    iRN := Subgroup( RQ, geniRN );
    if ( ( XModPrintLevel > 1 ) and ( diff > 0 ) ) then
        Print( "Closed generating set for N :-\n", genN, "\n" );
    fi;
    
    N := Subgroup( M, genN );
    
    # process Tdata
    if ( Tdata <> [ ] ) then
        T := Tdata;
    else
        T := CommonTransversal( RQ, iRP );
    fi;
    ok := IsCommonTransversal( RQ, iRP, T );
    if not ok then
        Error( "T fails to be a common transversal" );
    fi;
    if ( XModPrintLevel > 1 ) then
        Print( "Using transversal :- \n", T, "\n\n" );
    fi;
    # express each  q in RQ  as  p.t with  p in iRP, t in T
    qP := 0 * [1..oQ];
    qT := 0 * [1..oQ];
    for t in T do
        for p in eliRP do
            q := p*t;
            pos := Position( elRQ, q );
            qP[pos] := p;
            qT[pos] := t;
        od;
    od;
    if ( XModPrintLevel > 1 ) then
        Print( "qO = ", qP, "\nqT = ", qT, "\n\n" );
    fi;
    RQ.qP := qP;
    RQ.qT := qT;
    
    # set up Nrec
    Nrec := rec( );
    Nrec.perm := RN;
    Nrec.image := iRN;
    Nrec.fp := N;
    Nrec.free := fN;
    N.relators := relN;
    f2pN := GroupHomomorphismByImages( N, RN, genN, genRN );
    
    IXMdata.Qrec := Qrec;
    IXMdata.Prec := Prec;
    IXMdata.Mrec := Mrec;
    IXMdata.Nrec := Nrec;
    IXMdata.trans := T;
    IXMdata.xmod := X;
    IXMdata.iota := iota;
    IXMdata.isInducedXModData := true;
    return IXMdata;
end;            

###############################################################################
##
#F  InclusionInducedXModByCopower        standard construction for induced XMod
##

InclusionInducedXModByCopower := function( data )

    local  iota, imgeniota,
           Q, RQ, Qrec, oQ, genRQ, ngRQ, actQ, elRQ, p2fQ,
           P, RP, Prec, oP, indQP, qP, qT, t, T, words, total, tietze,
           M, RM, Mrec, oM, genRM, ngRM, xgM, imRM,
           N, RN, Nrec, genN, relN, nrN,
           fN, genfN, genRN, ngRN,
           iRP, iRN, geniRN,
           I, RI, Irec, oI, fI, ngI, xgI, nxI, genI, ordI, relI, degRI, imIQ,
           newI, newgenI, f2pI, presI, newngI, imRI, idRI, imact, genRI, 
           imRIRQ, homIRQ, kerIRQ, genK, genIK, IK, big,
           partition, numpart, factors, partI, ordI, grpI, IdI,
           ck, cl, cm, fk, fl, gk, gl, gm, hk, hl, hm,
           i, ik, il, im, j, jk, jl, jm, k, m, pm, rm,
           tk, tl, tm, zk, zl, zm, u, v, x, y, z, ix, iy, 
           imD, gimD, elimD, genD, ind, gpos, vpos, imrem, genpos, genim,
           X, IX, bdy, act, aut, mor, morsrc, Iname, Xname, GpId;
    
    if not (    IsBound( data.isInducedXModData )
            and data.isInducedXModData ) then
        return false;
    fi;
    Qrec := data.Qrec;
    Prec := data.Prec;
    Mrec := data.Mrec;
    Nrec := data.Nrec;
    RQ := Qrec.perm;
    RP := Prec.perm;
    iRP := Prec.image;
    RM := Mrec.perm;
    RN := Nrec.perm;
    iRN := Nrec.image;
    M := Mrec.fp;
    N := Nrec.fp;
    T := data.trans;
    X := data.xmod;
    iota := data.iota;
    qP := RQ.qP;
    qT := RQ.qT;
    genN := N.generators;
    relN := N.relators;
    nrN := Length( relN );
    fN := Nrec.free;
    genfN := fN.generators;
    genRN := RN.generators;
    geniRN := iRN.generators;
    genRM := RM.generators;
    genRQ := RQ.generators;
    oQ := Size( RQ );
    oP := Size( RP );
    indQP := oQ/oP;
    elRQ := RQ.elements;
    ngRN := Length( genRN );
    ngRM := Length( genRM );
    ngRQ := Length( genRQ );
    xgM := ngRN-ngRM;
    ngI := ngRN*indQP;
    nxI :=  xgM*indQP;

    fI := FreeGroup( ngI, "fI" );
    genI := [ ];
    for i in [1..indQP] do
        j := (i-1) * ngRN + 1;
        k := i*ngRN;
        Add( genI, Sublist( fI.generators, [j..k]) );
    od;
    ordI := 0 * [1..ngI];
    actQ := 0 * [1..ngRQ];
    for i in [1..ngRQ] do
        actQ[i] := Copy( ordI );
    od;
    for i in [1..ngRN] do
        ordI[i] := Order( RM, genRN[i] );
    od;
    for i in [(ngRN+1)..ngI] do
        ordI[i] := ordI[i-ngRN];
    od;
    if ( XModPrintLevel > 1 ) then
        Print( "Orders of the generators of I :- \n", ordI, "\n\n" );
    fi;
    # Images of the generators of I in Q
    imIQ := 0 * [1..ngI];
    for i in [1..indQP] do
        t := T[i];
        for j in [1..ngRN] do
            x := geniRN[j];
            y := x^t;
            k := (i-1)*ngRN + j;
            imIQ[k] := y;
        od;
    od;
    if ( XModPrintLevel > 1 ) then
        Print( "Images in RQ of the generators of I :- \n", imIQ, "\n" );
    fi;
    # Action of the generators of Q on the generators of I
    for jk in [1..ngRQ] do
        fk := genRQ[jk];
        for il in [1..indQP] do
            for jl in [1..ngRN] do
                tl := T[il];
                fl := geniRN[jl];
                cl := (il-1)*ngRN + jl;
                zl := tl*fk;
                pm := Position( elRQ, zl );
                tm := qT[pm];
                im := Position( T, tm );
                rm := qP[pm];
                zm := fl^rm;
                jm := Position( geniRN, zm );
                if ( jm = false ) then
                    Print( "\n\n !!! Position Error !!! \n\n" );
                fi;
                cm := (im-1)*ngRN + jm;
                actQ[jk][cl] := cm;
            od;
        od;
    od;
    Print( "\nAction of RQ on generators of I :- \n" );
    for i in [1..ngRQ] do
        Print( "  ", genRQ[i], " : ", PermList( actQ[i] ), "\n" );
    od;
    Print( "\n" );
    relI := [ ];
    for i in [1..indQP] do
        for j in [1..nrN] do
            u := relN[j];
            v := MappedWord( u, genfN, genI[i] );
            Add( relI, v );
        od;
    od;

    gimD := [];
    for ik in [1..indQP] do
        for jk in [1..ngRM] do                 # no longer [1..ngRN]
            tk := T[ik];
            fk := geniRN[jk];
            gk := genI[ik][jk];
            ck := (ik-1)*ngRN + jk;
            if ( ordI[ck] > 2 ) then
                hk := gk^-1;
            else
                hk := gk;
            fi;
            zk := fk^tk;
            gpos := Position( gimD, zk );
            if ( gpos = false ) then
                Add( gimD, zk );
            fi;
            for il in [1..indQP] do
                for jl in [1..ngRM] do         # no longer [1..ngRN]
                    tl := T[il];
                    fl := geniRN[jl];
                    gl := genI[il][jl];
                    cl := (il-1)*ngRN + jl;
                    if ( ordI[cl] > 2 ) then
                        hl := gl^-1;
                    else 
                        hl := gl;
                    fi;
                    zl := tl*zk;
                    m := Position( elRQ, zl );
                    tm := qT[m];
                    im := Position( T, tm );
                    rm := qP[m];
                    zm := fl^rm;
                    jm := Position( geniRN, zm );
                    if ( jm = false ) then
                        Print("\n\n !!! Position Error !!! \n\n");
                    fi;
                    gm := genI[im][jm];
                    cm := (im-1)*ngRN + jm;
                    if ( ordI[cm] > 2 ) then
                        hm := gm^-1;
                    else
                        hm := gm; 
                    fi;
                    if (ik <> il) then
                        big := cm;
                        if ( big < cl ) then 
                            big := cl; 
                        fi;
                        if ( big < ck ) then 
                            big := ck;
                        fi;
                        if ( big = cm ) then 
                            v := gm * hk * hl * gk;
                        elif ( big = cl ) then 
                            v := gl * gk * hm * hk;
                        elif ( ( hk = gk ) and ( cl < cm ) ) then
                            v := hk * hl * gk * gm;
                        else
                            v := gk * gm * hk * hl;
                        fi;
                        vpos := Position( relI, v );
                        if ( vpos = false ) then 
                            Add( relI, v );            # new relator!
                        fi;
                    fi;
                od;
            od;
        od;
    od;

    I := fI / relI;
    presI := PresentationFpGroup( I );
    TzInitGeneratorImages( presI );
    presI.protected := ngRM;
    Print( "#I  Protecting the first ", ngRM, " generators.\n" );
    presI.oldGenerators := Copy( presI.generators );
    if ( XModPrintLevel > 1 ) then
        Print( "\nFull set of relators for I :- \n", relI, "\n" );
        Print( "\nApplying PresentationFpGroup to I :- \n" );
        TzPrint( presI );
        Print( "\nApplying TzPartition & TzGo to I :- \n" );
    fi;

    i := 0;
    repeat
        i := i + 1;
        tietze := presI.tietze;
        total := tietze[TZ_TOTAL];
        TzPartition( presI );
        TzGo( presI );
    until ( ( total = tietze[TZ_TOTAL] ) or ( i > 9 ) );

    partition := TzPartition( presI );
    numpart := Length( partition );
    if ( numpart > 1 ) then
        factors := FactorsPresentation( presI );
        Print( "partitioning the generators: ", presI.partition, "\n" );
        if ( XModPrintLevel > 1 ) then
            Print( "factors = ", factors, "\n" );
        fi;
    else
        factors := [ presI ];
    fi;

    Print( "\nSimplified presentation for I :- \n");
    TzPrint( presI );
    newI := FpGroupPresentation( presI );
    oI := Size( newI );
    newgenI := newI.generators;
    newngI := Length( newgenI );
    Print( "\n I has Size: ", oI, "\n");
    Print( "**************** \n\n ");
    GpId := [ ];
    for i in [1..numpart] do
        partI := factors[i];
        TzGo( partI );
        grpI := FpGroupPresentation( partI );
        ordI := Size( grpI );
        if ( ordI > 1 ) then
            Irec := FpPair( grpI );
            RI := Irec.perm;
            degRI := PermGroupOps.NrMovedPoints( RI );
            f2pI := Irec.f2p;
        fi;
        # now identify grpI (if possible)
        if ( ordI > 1 ) then
            if  IsAbelian( I ) then
                Print( "factor ", i, " is abelian" );
                IdI := GroupId( RI );
                Print( "  with invariants: ", IdI.abelianInvariants, "\n" );
                Add( GpId, IdI );
            elif ( ordI <= 100 ) then
                Print( "Searching Solvable Groups Library:\n" );
                IdI := GroupId( RI );
                Add( GpId, IdI );
            fi;
        fi;
    od;
    if not IsAbelian( I ) then
        Print( "\n GroupId = \n" );
        for i in [1..numpart] do
            Print( GpId[i], "\n" );
        od;
    fi;
    if ( oI > 1 ) then
        Irec := FpPair( newI );
        RI := Irec.perm;
        if IsBound( RM.name ) then
            RI.name := Concatenation( "i*(", RM.name, ")" );
        else
            RI.name := "i*RM";
        fi;
        RI.fpPair := Irec;
        degRI := PermGroupOps.NrMovedPoints( RI );
        f2pI := Irec.f2p;
        imD := Subgroup( RQ, gimD );
        ind := Index( RQ, imD );
        Print( "\nImage of I has index ", ind, " in RQ " );
        Print( "and is generated by : \n", gimD, "\n" );
        genD := 0 * [1..newngI];
        for i in [1..newngI] do
            x := newgenI[i];
            j := Position( presI.oldGenerators, x );
            genD[i] := imIQ[j];
        od;
        # Print( "\n\npresI.oldGenerators, genD = \n" );
        # Print( presI.oldGenerators,"\n", genD,"\n\n" );
        homIRQ := GroupHomomorphismByImages( newI, RQ, newgenI, genD );
        kerIRQ := Kernel( homIRQ );
        # Print( "kerIRQ = \n", kerIRQ, "\n" );
        genK := kerIRQ.generators;
        genIK := List( genK, k -> Image( f2pI, k ) );
        # Print( "genK, genIK = \n", genK, "\n", genIK, "\n" );
        IK := Subgroup( RI, genIK );
        # Print( "IK =\n", IK, "\n" );
        words := List( presI.imagesOldGens,
                         w -> AbstractWordTietzeWord( w, presI.generators ) );
        # Print( "\n     words = ", words, "\n" ); 
        imrem := List( words, w -> Image( f2pI, w ) );
        if ( XModPrintLevel > 1 ) then
            Print( "Initial generators in terms of final generators :- \n" );
            for i in [ 1..Length( presI.imagesOldGens ) ] do
                Print( presI.oldGenerators[i], " : " );
                Print( presI.imagesOldGens[i], " --> ", imrem[ i ], "\n" );
            od;
        fi;
    fi;

    idRI := IdentityMapping( RI );
    genRI := RI.generators;
    genpos := List( genRI, g -> Position( imrem, g ) );
    if ( XModPrintLevel > 2 ) then
        Print( "genpos = ", genpos, "\n" );
    fi;
    imRI := 0 * [1..ngRQ];
    imact := 0 * [1..ngRQ];
    for i in [1..ngRQ] do
        imRI[i] := List( actQ[i], j -> imrem[j] );
        genim := List( genpos, p -> imRI[i][p] );
        imact[i] := GroupHomomorphismByImages( RI, RI, genRI, genim );
    od;
    imRIRQ := List( genpos, p -> imIQ[p] );
    if ( XModPrintLevel > 2 ) then
        Print( "imIQ = ", imIQ, "\n" );
        Print( "imRIRQ = ", imRIRQ, "\n" );
    fi;
    bdy := GroupHomomorphismByImages( RI, RQ, genRI, imRIRQ );
    imRM := Sublist( imrem, [1..ngRM] );
    morsrc := GroupHomomorphismByImages( RM, RI, genRM, imRM );
    aut := Group( imact, idRI );
    act := GroupHomomorphismByImages( RQ, aut, genRQ, imact );

    IX := rec( );
    IX.source := RI;
    IX.range := RQ;
    IX.boundary := bdy;
    IX.action := act;
    IX.aut := aut;
    IX.kernel := IK;
    IX.isDomain := true;
    IX.operations := XModOps;
    IX.isInducedXMod := true;
    IX.xmod := X;
    IX.name := Concatenation( "i*(", XModName( X ), ")" );
    IX.isXMod := IsXMod( IX );
    mor := XModMorphism( X, IX, [morsrc, iota] );
    if not IsXModMorphism( mor ) then
        Print( "mor: X -> IX  not an xmod morphism!\n" );
    fi;
    IX.morphism := mor;
    if IsBound( RI.groupId ) then
        IX.sourceId := RI.groupId;
    fi;
    return IX;
end;

##############################################################################
##
#F  SurjectiveInducedXMod                   constructs  M/[M,K], K = ker(iota)
##

SurjectiveInducedXMod := function( X, iota )

    local  S, genS, R, bdy, act, K, genK, s, r, a, x, H, genH, rcos, reps,
           Q, lenQ, genQ, preQ, I, genI, imi, istar, ophom, imb, bdystar,
           i, autgen, imI, imS, aut, actstar, autstar, idI, IX, mor;

    R := X.range;
    S := X.source;
    genS := S.generators;
    bdy := X.boundary;
    act := X.action;
    K := Kernel( iota );
    genK := K.generators;
    genH := [ ];
    H := Subgroup( S, genH );
    for r in genK do
        a := Image( act, r );
        for s in genS do
            x := s^(-1) * Image( a, s );
            if not ( x in H ) then
                Add( genH, x );
                H := Subgroup( S, genH );
            fi;
        od;
    od;
    Q := iota.range;
    genQ := Q.generators;
    preQ := List( genQ, q -> PreImagesRepresentative( iota, q ) );
    rcos := RightCosets( S, H );
    reps := List( rcos, r -> r.representative );
    imb := List( genS, r -> Image( iota, Image( bdy, r ) ) );
    I := Operation( S, rcos, OnRight );
    genI := I.generators;
    if IsBound( S.name ) then
        I.name := Concatenation( S.name, "/ker" );
    fi;
    ophom := OperationHomomorphism( S, I );
    imi := List( genS, s -> Image( ophom, s ) );
    if ( genI <> imi ) then
        Error( "unequal images: genI <> imi" );
    fi;
    istar := GroupHomomorphismByImages( S, I, genS, imi );
    bdystar := GroupHomomorphismByImages( I, Q, genI, imb );
    lenQ := Length( genQ );
    autgen := 0 * [1..lenQ];
    for i in [1..lenQ] do
        r := preQ[i];
        a := Image( act, r );
        imS := List( genS, s -> Image( a, s ) );
        imI := List( imS, s -> Image( ophom, s ) );
        aut := GroupHomomorphismByImages( I, I, genI, imI );
        autgen[i] := aut;
    od;
    idI := InclusionMorphism( I, I );
    autstar := Group( autgen, idI );
    actstar := GroupHomomorphismByImages( Q, autstar, genQ, autgen );
    IX := XMod( bdystar, actstar );
    IX.isInducedXMod := true;
    IX.xmod := X;
    IX.name := Concatenation( "i*(", XModName( X ), ")" );
    mor := XModMorphism( X, IX, [istar,iota] );
    IX.morphism := mor;
    return IX;
end;

##############################################################################
##
#F  InducedXMod           calls various functions to construct an induced xmod
##

InducedXMod := function( arg )

    local  usage1, usage2, usage3, usage4, nargs, data, mono, surj,
           Qdata, Pdata, Mdata, Tdata, X, iota, IX, Rdata, genP, imQ,
           ires, JX, JXsrc, JXrng, Jsrcgen, Ksrcgen, KXsrc, KX,
           J2Ksrc, J2Krng, inc;

    usage1 := "\nUsage: InducedXModData( Qdata, Pdata, Mdata [, Tdata] );";
    usage2 := "\n where  Q >= P >= M  and  T  is a transversal for Q/P";
    usage3 := "\n   or: InducedXModData( X, iota [, Tdata] );";
    usage4 := "\n where  X  is a conjugation XMod and iota a monomorphism.\n";
    nargs := Length( arg );
    if ( nargs > 4 ) then
        Print( usage1, usage2, usage3, usage4 );
        return false;
    fi;
    if not IsRec( arg[1] ) then
        Print( usage1, usage2, usage3, usage4 );
        return false;
    fi;
    Tdata := [ ];
    if IsGroup( arg[1] ) then
        if ( ( nargs < 3 ) or not IsNormal( arg[2], arg[1] ) ) then
            Print( usage1, usage2, usage3, usage4 );
            return false;
        fi;
        Qdata := arg[1];
        Pdata := arg[2];
        Mdata := arg[3];
        X := ConjugationXMod( Pdata, Mdata );
        iota := InclusionMorphism( Pdata, Qdata );
        if ( nargs = 4 ) then
            Tdata := arg[4];
        fi;
    elif IsXMod( arg[1] ) then
        X := arg[1];
        Mdata := X.source;
        Pdata := X.range;
        iota := arg[2];
        if ( ( nargs > 3 ) or ( iota.source <> Pdata ) ) then
            Print( usage1, usage2, usage3, usage4 );
            return false;
        fi;
        Qdata := iota.range;
        if ( nargs = 3 ) then
            Tdata := arg[3];
        fi;
    fi;

    mono := IsMonomorphism( iota );
    surj := IsSurjective( iota );
    if mono then
        data := InclusionInducedXModData( X, iota, Tdata );
        IX := InclusionInducedXModByCopower( data );
    elif surj then
        IX := SurjectiveInducedXMod( X, iota );
    else  # split in two
        Rdata := Image( iota );
        genP := Pdata.generators;
        imQ := iota.genimages;
        ires := GroupHomomorphismByImages( Pdata, Rdata, genP, imQ );
        JX := SurjectiveInducedXMod( X, ires );
        if ( XModPrintLevel > 1 ) then
            Print( "Calculating SurjectionIndXMod then InclusionIndXMod\n" );
            XModPrint( JX );
        fi;
        JXsrc := JX.source;
        JXrng := JX.range;
        Jsrcgen := JXsrc.generators;
        Ksrcgen := List( Jsrcgen, j -> Image( JX.boundary, j ) );
        KXsrc := Subgroup( JXrng, Ksrcgen );
        J2Ksrc := GroupHomomorphismByImages( JXsrc, KXsrc, Jsrcgen, Ksrcgen );
        J2Krng := InclusionMorphism( JXrng, JXrng );
        KX := IsomorphicXMod( JX, [ J2Ksrc, J2Krng ] );
        inc := InclusionMorphism( Rdata, Qdata );
        data := InclusionInducedXModData( KX, inc, [ ] );
        IX := InclusionInducedXModByCopower( data );
        if IsBound( Mdata.name ) then
            IX.source.name := Concatenation( "i*(", Mdata.name, ")" );
        else
            IX.source.name := "i*RM";
        fi;
        XModName( IX );
    fi;
    return IX;
end;

##############################################################################
##
#F  AllInducedXMods                      given Q, finds all XMods  M <= P <= Q
##

AllInducedXMods := function( arg )

    local  nargs, usage, Q, L, lenL, reps, nreps, r, i, j, k, a, b,
           norm, nnorm, n, sizes, keep, coll, P, M, id, X, num, line;

    usage := "\nUsage:  AllInducedXMods( Q [, P ] );\n";
    nargs := Length( arg );
    if ( ( nargs < 1 ) or ( nargs > 2 ) ) then
        Print( usage );
        return false;
    fi;
    line := "--------------------------------------";
    Q := arg[1];
    if ( nargs = 2 ) then
        P := arg[2];
        if not IsSubgroup( Q, P ) then
            Print( usage );
            Error( "P is not a subgroup of Q" );
        fi;
        id := Subgroup( P, [ ] );
        reps := [ Q, P, id ];
        nreps := 3;
    else
        L := Lattice( Q );
        reps := Reversed( List( L.classes, c -> c.representative ) );
        nreps := Length( reps );
        Print( "non-trivial reps = \n" );
        for r in [ 2 .. nreps-1 ] do
            Print( reps[r], "\n" );
        od;
    fi;
    num := 0;
    Print( "\nAll induced crossed modules  M --> IM" );
    Print( "\n                             |     | " );
    Print( "\n                             P --> Q\n");
    Print( "\ngenQ = ", Q.generators, "\n" );
    Print( "\n", line, line, "\n\n" );
    for r in [ 2 .. nreps-1 ] do
        P := reps[r];
        norm := NormalSubgroups( P );
        # find representatives of conjugacy classes in Q
        sizes := List( norm, n -> Size( n ) );
        coll := Collected( sizes );
        keep := List( norm, n -> true );
        k := 1;
        for i in [ 2 .. ( Length( coll ) - 1 ) ] do
            j := k + 1;
            k := k + coll[i][2];
            for a in [ j .. k-1 ] do
                if keep[a] then
                    for b in [ a+1 ..k ] do
                        if IsConjugate( Q, norm[a], norm[b] ) then
                            keep[b] := false;
                        fi;
                    od;
                fi;
            od;
        od;
        nnorm := Length( norm );
        norm[ nnorm ] := P;
        for n in [2..nnorm] do
            if keep[n] then
                M := norm[n];
                Print( "genM = ", M.generators, "\n" );
                Print( "genP = ", P.generators, "\n" );
                X := InducedXMod( Q, P, M );
                XModPrint( X );
                num := num + 1;
                Print( line, line, "\n\n" );
            fi;
        od;
    od;
    Print( "Number of induced crossed modules calculated = ", num, "\n" );
    return num;
end;

############################################################################
##
#F  InducedCat1Data      set up the data for function InducedCat1
##

InducedCat1Data := function( arg )

local  Qdata, Rdata, Gdata,        # 3 permutation groups
       C, iota,                    # C = [Gdata ==iota==> Rdata]
       ICGdata,                    # InducedCat1Data
       nargs,                      # number of arg
       usage,                      # usage message
       Q, R, G,                    # Fin. Pres. Group
       PQ, PR, PG,                 # Perm groups
       elPQ, elPR, elPG,           # elements of perm groups
       genPQ, genPR, genPG,        # generating set of perm groups
       oQ, oR, oG,                 # size of perm groups
       Qrec, Rrec, Grec,           # record fields 
       ngPG,                       # number of generating set for PG
       indPQ,                      # oQ/oR
       degPQ,                      # PermGroupOps.NrMovedPoints( PQ )
       genG,                       # generating set for Fin.Pres. Group G
       ngG,                        # number of generating set for G
       relG,                       # relators for G
       nrG,                        # number of all relators
       elG,                        # elements           
       presG,                      # PresentationViaCosetTable for G
       f2pG,                       # record field for G
       qR, qT,                     # record field for PQ
       posp2f, posf2p,             # positions of isomorphic images G <-> PG
       pos,                        # position variable
       T,                          # Tdata
       ok,                         # checking variable
       t, i, q, p, m, rm;          # variables

    usage := "\nUsage: InducedCat1Data( C, iota );\n";
    nargs := Length( arg );
    if ( nargs <> 2 ) then
        Print( usage );
        return false;
    fi;
    C := arg[1];
    Gdata := C.source;
    Rdata := C.range;
    iota := arg[2];
    Qdata := iota.range;

    ICGdata := rec( );
    
    # process Qdata
    if ( IsBound( Qdata.isFpPair ) and Qdata.isFpPair ) then
        Qrec := Qdata;
    elif ( IsBound( Qdata.isPermGroup ) and Qdata.isPermGroup ) then
        Qrec := rec( );
        Qrec.perm := Qdata;
    else
        Error( " Qdata not of correct type" );
    fi;
    PQ := Qrec.perm;
    oQ := Size( PQ );
    elPQ := Elements( PQ );
    genPQ := PQ.generators;
    degPQ := PermGroupOps.NrMovedPoints( PQ );
    
    # process Rdata
    if ( IsBound( Rdata.isFpPair ) and Rdata.isFpPair ) then
        Rrec := Rdata;
    elif ( IsBound( Rdata.isPermGroup ) and Rdata.isPermGroup ) then
        Rrec := rec( );
        Rrec.perm := Rdata;
    else
        Error( " Rdata not of correct type" );
    fi;
    PR := Rrec.perm;
   # if not IsSubgroup( PQ, PR ) then
   #     Error( " PR not a subgroup of PQ" );
   # fi;
    oR := Size( PR );
    elPR := Elements( PR );
    genPR := PR.generators;
    indPQ := oQ/oR;
    
    # process Gdata
    if ( IsBound( Gdata.isFpPair ) and Gdata.isFpPair ) then
        Grec := Gdata;
    elif ( IsBound( Gdata.isPermGroup ) and Gdata.isPermGroup ) then
        Grec := rec( );
        Grec.perm := Gdata;
    fi;
    PG := Grec.perm;
    #if not IsNormal( PR, PG ) then
    #    Error( " PR not a normal subgroup of PG" );
    #fi;
    oG := Size( PG );
    elPG := Elements( PG );
    genPG := PG.generators;
    ngPG := Length( genPG );
    
    # find a presentation for PG
    if ( IsBound( PG.fpPair ) and IsFpPair( PG.fpPair ) ) then
        Grec := PG.fpPair;
        G := Grec.fp;
        genG := G.generators;
        ngG := Length( genG );
        f2pG := Grec.f2p;
    else
        presG := PresentationViaCosetTable( PG );
        TzInitGeneratorImages( presG );
        G := FpGroupPresentation( presG );
        if ( ngPG <> Length( presG.generators ) ) then
            Print( "\n\nWARNING: extra generators defined - \n" );
            Print( "             the isomorphism f2pG may be incorrect!\n\n" );
        fi;
        Grec := rec( );
        Grec.perm := PG;
        Grec.fp := G;
        genG := G.generators;
        ngG := Length( genG );
        f2pG := GroupHomomorphismByImages( G, PG, genG, genPG );
        if not IsIsomorphism( f2pG ) then
            Error( "f2pG not an isomorphism" );
        fi;
        Grec.f2p := f2pG;
        Grec.isFpPair := IsFpPair( Grec );
    fi;
    
    presG := PresentationFpGroup( G );
    TzInitGeneratorImages( presG );
    if ( XModPrintLevel > 2 ) then
        Print( "initial presentation for G :-\n\n" );
        TzPrint( presG );
    fi;
    relG := G.relators;
    nrG := Length( relG );
    elG := Elements( G );   

    # Determine the positions of isomorphic images G <-> PG
    posf2p := 0 * [1..oG];
    posp2f := 0 * [1..oG];
    for i in [1..oG] do
        m := elG[i];
        rm := Image( f2pG, m );
        pos := Position( elPG, rm );
        posf2p[i] := pos;
        posp2f[pos] := i;
    od;
    if ( XModPrintLevel > 1 ) then
        Print( "posf2p = ", posf2p, "\n" );
        Print( "posp2f = ", posp2f, "\n" );
    fi;

    ICGdata.Qrec := Qrec;
    ICGdata.Rrec := Rrec;
    ICGdata.Grec := Grec;
    ICGdata.cat1 := C; 
    ICGdata.iota := iota;
    ICGdata.isInducedCat1Data := true;
    return ICGdata;
end;            

##############################################################################
##
#F  InducedCat1ByFreeProduct              construct an induced cat1-group
##

InducedCat1ByFreeProduct := function( data )

    local  Q,              # Fin. Pres. Group
           PQ,             # Perm Group 
           Qrec,           # Record field 
           oQ,             # Size of Q
           genPQ,          # generating set for Perm group Q 
           ngPQ,           # number of generating set of Perm Group Q
           genPG,          # generating set perm group G
           ngPG,           # number of generating set of perm G
           ngI,            # total length of ngPG+ngPQ
           C,              # Cat1 
           Csrc,           # Cat1 source 
           Crng,           # Cat1 range 
           t, h, e,        # tail, head, embedding 
           genCsrc,        # generating set of sourxe group
           genCrng,        # generating set of range group
           fI,             # Free group 
           genfI,          # generating set of free group 
           imPG, imPQ, 
           relI,           # all relators
           Gfp,            # FpPair 
           genGfp,         # generating set  
           Grel,           # relators
           len,            # Length of relators
           Qfp,            # FpPair for Q
           genQfp,         # generating set 
           Qrel,           # relators
           iota,            # inclusion map from Crng to PQ
           imembed,        # relations produced by embedding
           imiota,         # relations produced by iota
           uuQ,            # List variable
           wQ, uQ, wG, uG, uQG,
           kert, kerh,     # kernel of tail and head
           genkert,        # generating set of kert
           genkerh,        # generating set of kerh
           imt, imh, 
           tG, hG, 
           com,            # Commutator subgroup
           Yh, Yt,         # conjugations for tail and head
           YYt, YYh,       # List of conjugations
           I, genI,        # new free group and its generating set
           presI,          # Presentation
           newIfp,         # FpPair
           PI,             # new permutational group
           oI, genPI,      # Size and generating set of new perm group
           iotastar,       # homomorphism from Csrc to nep perm group
           imh1, imh2, 
           hstar,          # new head homomorphism for Induced Cat1 group
           imt1, imt2, 
           tstar,          # new tail homomorphism for Induced Cat1-Group
           imm, imag, images, 
           estar,          # new embed homomorphism for Ind.cat1
           IC,              # Induced Cat1 variable
           mor,            # Cat1Morphism from C to IC 
           u, v, j, x, i, g;# using variables

    PQ := data.Qrec.perm;
    C := data.cat1;
    iota := data.iota;
    Csrc := C.source;
    Crng := C.range;
    t := C.tail;
    h := C.head;
    e := C.embedRange;
    genPQ := PQ.generators;
    genCsrc := Csrc.generators;
    genCrng := Crng.generators;
Print( "genCrng = ", genCrng, "\n" );
    ngPG := Length( genCsrc );
    ngPQ := Length( genPQ );
    ngI := ngPG+ngPQ;
    fI := FreeGroup( ngI, "fI" );
    genfI := fI.generators;
    imPQ := Sublist( genfI, [ngPG+1..ngI] );
    imPG := Sublist( genfI, [1..ngPG] ); 
    relI := [ ];
    # Creating the relations of G
    Gfp := FpPair( Csrc );
    genGfp := Gfp.fp.generators;
    Grel := Gfp.fp.relators;
    len := Length( Grel );
    for j in [1..len] do
        u := Grel[j];
        v := MappedWord( u, genGfp, imPG );
        Add( relI, v );
    od;
    # Adding extra relations from Q
    Qfp := FpPair( PQ );
    genQfp := Qfp.fp.generators;
    Qrel := Qfp.fp.relators;
    len := Length( Qrel );
    for j in [1..len] do
        u := Qrel[j];
        v := MappedWord( u, genQfp, imPQ );
        Add( relI, v );
    od;
    # Adding extra relations from embedding and iota
    uuQ := [ ];
    imembed := List( genCrng, x -> Image( e, x ) );
    wG := List( imembed, x -> Image( Gfp.p2f, x ) ); 
    uG := List( wG, g -> MappedWord( g, genGfp, imPG ) );
    imiota := List( genCrng, x -> Image( iota, x ) ); 
    wQ := List( imiota, x -> Image( Qfp.p2f, x ) );
    uQ := List( wQ, u -> MappedWord( u, genQfp, imPQ ) );
    for i in [1..Length(uG)] do     
        uQG := uG[i]*uQ[i]^-1;
        Add( uuQ, uQG );
    od;
    relI := Concatenation( relI, uuQ );
    # Finding the Peiffer subgroup
    YYt := [ ];
    YYh := [ ];
    kert := Kernel( t );
    kerh := Kernel( h );
    genkert := kert.generators;
    genkerh := kerh.generators;
    imt := List( genkert, x -> Image( Gfp.p2f, x ) );
    tG := List( imt, i -> MappedWord( i, genGfp, imPG ) );
    imh := List( genkerh, x -> Image( Gfp.p2f, x ) );
    hG := List( imh, i -> MappedWord( i, genGfp, imPG ) );
    for u in genfI do
    Yt := List( tG, x -> x^u );
    YYt := Concatenation( YYt, Yt );
    Yh := List( hG, x -> x^u );
    YYh := Concatenation( YYh, Yh );
    od;
    for i in YYt do
        for j in YYh do
            com := Comm( i, j );
            Add( relI, com );
        od;
    od;
    I := fI / relI;
    genI := I.generators;
    presI := PresentationFpGroup( I );
    TzInitGeneratorImages( presI );
    presI.protected := Length( genI );
    Print( "#I  Protecting the first ", ngPG, " generators.\n" );
    presI.oldGenerators := Copy( presI.generators );
     if ( XModPrintLevel > 1 ) then
        Print( "\nFull set of relators for I :- \n", relI, "\n" );
        Print( "\nApplying PresentationFpGroup to I :- \n" );
        TzPrint( presI );
        Print( "\nApplying TzPartition & TzGo to I :- \n" );
    fi;
    TzPrint( presI ); 
    TzGoGo( presI );
    TzPrint( presI );
    imPQ := Sublist( genI, [ngPG+1..ngI] );
    imPG := Sublist( genI, [1..ngPG] ); 
    newIfp := FpPair( I );
    PI := newIfp.perm;
    oI := Size( PI );
    genPI := PI.generators;
    Print("new perm group size ", oI, "\n");
    Print("******************** \n");
    imPG := Sublist( genPI, [1..ngPG] );
    iotastar := GroupHomomorphismByImages( Csrc, PI, genCsrc, imPG );
    imh1 := List( genCsrc, x -> Image( h, x ) );
    imh2 := List( imh1, x -> Image( iota, x ) );
    imh := Concatenation( imh2, genPQ );
    hstar := GroupHomomorphismByImages( PI, PQ, genPI, imh );
    imt1 := List( genCsrc, x -> Image( t, x ) );
    imt2 := List( imt1, x -> Image( iota, x ) );
    imt := Concatenation( imt2, genPQ );
    tstar := GroupHomomorphismByImages( PI, PQ, genPI, imt );
    imm := List( genPQ, x -> Image( Qfp.p2f, x ) );
    imag := List( imm, x -> MappedWord( x, genQfp, imPQ ) );
    images := List( imag, x -> Image( newIfp.f2p, x ) );
    estar := GroupHomomorphismByImages( PQ, PI, genPQ, images );
    IC := Cat1( PI, tstar, hstar, estar );
    IC.isCat1 := IsCat1( IC );
    mor := Cat1Morphism( C, IC, [ iotastar, iota ] );
    if not ( IsCat1Morphism( mor ) ) then
        Print( " mor : C --> IC not a cat1-group morphism \n" );
    fi;
    IC.morphism := mor;
    IC.name := Concatenation( "<ICG(", Cat1Name( C ), ")>" );
    IC.isDomain := true;
    IC.operations := Cat1Ops;
    IC.cat1 := C;
    IC.isInducedCat1 := true;
    return IC;
end;

##############################################################################
##
#F  InducedCat1      calls various functions to construct an induced xmod
##

InducedCat1 := function( arg )

    local  usage1, usage2, usage3, usage4, nargs, data,
           Qdata, Pdata, Gdata, C, iota, IC;

    usage1 := "\nUsage: InducedCat1( Qdata, Pdata, Gdata );";
    usage2 := "\n where  Q >= P  and  C = [ G ==> P ]";
    usage3 := "\n   or: InducedCat1( C, iota );";
    usage4 := "\n where  C=[G==>P] is a Cat1, iota:P->Q an inclusion.\n";
    nargs := Length( arg );
    if ( nargs > 2 ) then
        Print( usage1, usage2, usage3, usage4 );
        return false;
    fi;
    if not IsRec( arg[1] ) then
        Print( usage1, usage2, usage3, usage4 );
        return false;
    fi;
    if IsGroup( arg[1] ) then
        if ( ( nargs < 3 ) or not IsNormal( arg[2], arg[1] ) ) then
            Print( usage1, usage2, usage3, usage4 );
            return false;
        fi;
        Qdata := arg[1];
        Pdata := arg[2];
        Gdata := arg[3];
        C := ConjugationCat1( Pdata, Gdata );
        Gdata := C.source;
        iota := InclusionMorphism( Pdata, Qdata );
    elif IsCat1( arg[1] ) then
        C := arg[1];
        Gdata := C.source;
        Pdata := C.range;
        iota := arg[2];
        #if ( ( nargs > 2 ) or ( iota.source <> Pdata )
        #                   or not IsMonomorphism( iota ) ) then
        #    Print( usage1, usage2, usage3, usage4 );
        #    return false;
        #fi;
        Qdata := iota.range;
    fi;

    data := InducedCat1Data( C, iota );
    IC := InducedCat1ByFreeProduct( data );
    return IC;
end;

##############################################################################
##
#F  AllInducedCat1s                     
##

AllInducedCat1s := function( arg )

local  nargs, usage, Q, L, lenL, reps, nreps, r, i, j, k, a, b,
           norm, nnorm, n, sizes, keep, coll, P, M, id, data, 
           IC, num, line, C, iota;

    usage := "\nUsage:  AllInducedCat1s( Q [, P ] );\n";
    nargs := Length( arg );
    if ( ( nargs < 1 ) or ( nargs > 2 ) ) then
        Print( usage );
        return false;
    fi;
    line := "--------------------------------------";
    Q := arg[1];
    if ( nargs = 2 ) then
        P := arg[2];
        if not IsSubgroup( Q, P ) then
            Print( usage );
            Error( "P is not a subgroup of Q" );
        fi;
        id := Subgroup( P, [ ] );
        reps := [ Q, P, id ];
        nreps := 3;
    else
        L := Lattice( Q );
        reps := Reversed( List( L.classes, c -> c.representative ) );
        nreps := Length( reps );
        Print( "non-trivial reps = \n" );
        for r in [ 2 .. nreps-1 ] do
            Print( reps[r], "\n" );
        od;
    fi;
    num := 0;
    Print( "\n All induced cat1-groups       M --> IM" );
    Print( "\n                              ||    || " );
    Print( "\n                               P --> Q\n");
    Print( "\n genQ = ", Q.generators, "\n" );
    Print( "\n", line, line, "\n\n" );
    for r in [ 2 .. nreps-1 ] do
        P := reps[r];
        norm := NormalSubgroups( P );
        # find representatives of conjugacy classes in Q
        sizes := List( norm, n -> Size( n ) );
        coll := Collected( sizes );
        keep := List( norm, n -> true );
        k := 1;
        for i in [ 2 .. ( Length( coll ) - 1 ) ] do
            j := k + 1;
            k := k + coll[i][2];
            for a in [ j .. k-1 ] do
                if keep[a] then
                    for b in [ a+1 ..k ] do
                        if IsConjugate( Q, norm[a], norm[b] ) then
                            keep[b] := false;
                        fi;
                    od;
                fi;
            od;
        od;
        nnorm := Length( norm );
        norm[ nnorm ] := P;
        for n in [2..nnorm] do
            if keep[n] then
                M := norm[n];
                Print( "genM = ", M.generators, "\n" );
                Print( "genP = ", P.generators, "\n" );
                C := ConjugationCat1(P, M);
                iota := InclusionMorphism( P, Q );
                IC := InducedCat1( C, iota );
                Cat1Print( IC );
                num := num + 1;
                Print( line, line, "\n\n" );
            fi;
        od;
    od;
    Print( "Number of induced cat1-groups calculated = " );
    return num;
end;


# end of file:  induce.g
