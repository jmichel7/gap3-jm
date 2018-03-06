#############################################################################
##
#A  wylat.g                 CrystGap library                     Bettina Eick
#A                                                              Franz G"ahler
#A                                                              Werner Nickel
##
#Y  Copyright 1990-1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  CrystGap - the crystallographic groups package for GAP (Wyckoff lattice)
##  

#############################################################################
##
## Routines for the determination of a Wyckoff lattice
##
#############################################################################


#############################################################################
##
#F  IsSubspaceAffineSubspace( <U>, <V> ) . . . . . . . . .  V subspace of U ?
##
IsSubspaceAffineSubspace := function( U, V )

    local W, w, h, u;

    # span(V.basis) contained in span(U.basis)?
    W := ShallowCopy( V.basis );
    for w in W do
        h:=1;
        for u in U.basis do
            while u[h] = 0 do h := h+1; od;
            if w[h] <> 0 then
                w := w - w[h] * u / u[h];
            fi;
        od;
        if w<>0*w then
            return false;
        fi;
    od;

    # are translations compatible?
    w := U.translation - V.translation;
    h := 1;
    for u in U.basis do
        while u[h] = 0 do h := h+1; od;
        if w[h] <> 0 then
            w := w - w[h] * u / u[h];
        fi;
    od;
    w := List( w, FractionModOne );

    return w=0*w;

end;


#############################################################################
##
#F  WyckoffPosRelations( <S>, <W> ) . . . Incidence rels of Wyckoff positions
##
WyckoffPosRelations := function( S, W )
    local m, O, o, i, j, k, index;

    if not S.isStandard then
        W := Copy( W );
        ConjugateWyckoffPositions( W, S.translations^-1 );
    fi;

    m := NullMat( Length(W), Length(W) );
    for i in [1..Length(W)] do
        O := OrbitAffineSpaces( S, W[i] );
        for j in [1..Length(W)] do
            index := Size(W[j].stabilizer) / Size(W[i].stabilizer);
            if Length(W[j].basis) < Length(W[i].basis) and IsInt(index) then
                if ForAny( O, o -> IsSubspaceAffineSubspace( o, W[j] ) ) then
                    m[j][i] := index;
                fi;
            fi;
        od;
    od;

    for i in Reversed([1..Length(W)]) do
        for j in Reversed([1..i-1]) do
            if m[j][i]<>0 then
                for k in [1..j-1] do
                    if m[k][j]<>0 then m[k][i]:=0; fi;
                od;
            fi;
        od;
    od;

    return m;

end;


#############################################################################
##
#F  WyckoffLatticeRecord( <S> ) . . . . . . Create record for Wyckoff lattice
##
WyckoffLatticeRecord := function( S )

    local W, m, lat, i, j, max, w, classes, classlen, classpos;

    W := ShallowCopy( WyckoffPositions( S ) );
    Sort( W, function(a,b) return a.class > b.class; end );

    m := WyckoffPosRelations( S, W );
    lat := rec( vertices := W,
                sizes    := List( W, w -> Size(w.stabilizer) ) );
    max := List( W, w -> [] );
    for i in [1..Length(W)] do
        for j in [1..i-1] do
            if m[j][i]<>0 then Add(max[j],i); fi;
        od;
    od;
    lat.maximals := max;
    classes  := Set( List(W, w -> w.class) );
    classlen := List( classes, w -> 0 );
    classpos := [1..Length(W)];
    for i in [1..Length(W)] do
        j := Position( classes, W[i].class );
        classlen[j] := classlen[j]+1;
        classpos[i] := [ j, classlen[j] ];
    od;
    lat.classPositions := classpos;

    return lat;
end;


#############################################################################
##
#F  WyckoffLatticeOps . . . . . . . .  operations record for Wyckoff lattices
##
WyckoffLatticeOps := Copy(GraphicLatticeRecordOps);


#############################################################################
##
#F  WyckoffLatticeOps.MakeX( <sheet> )  . . . . . .  make x coordinates
##
WyckoffLatticeOps.MakeX := function( sheet )
    
    local   i,  j,  k,  l,	# loop variable
            ri,                 # local rows
            cl,                 # class length
            coor, coors,        # local x coordinate list
            x,                  # x coordinate
            n,                  # number of branches
            pos1,  pos2,        # positions
            noClasses,          # number of classes of a given order
            noElements,         # number of vertices of a given order
            tmp;
    
    # if the lattice is trivial (one or two vertices) return
    if Length(sheet.lattice.vertices) <= 2 then
        sheet.init.x := [ QuoInt(sheet.width, 2), QuoInt(sheet.width, 2) ];
        return;
    fi;

    # count the number of classes of a given order
    noClasses := List( sheet.init.orders, x -> 
                       [ Number( sheet.init.orderReps, y->y=x ), x ] );

    # count the number of vertices of a given order
    noElements := List( sheet.init.orders, x ->
                        Number( sheet.init.vertices, y -> y[2] = x ) );

    # compute the layers <sheet>.l
    sheet.init.l := [];
    for i  in [ 1 .. Length(FactorsInt(Maximum(sheet.init.orders))) ]  do
        sheet.init.l[i] := [];
    od;
    for i  in [ 2 .. Length(sheet.init.orderReps) ]  do
        AddSet( sheet.init.l[Length(FactorsInt(sheet.init.orderReps[i])) ],
                sheet.init.orderReps[i] );
    od;
    sheet.init.l := Filtered( sheet.init.l, i -> 0 < Length(i) ); 
    
    # compute the branches <sheet>.b
    sheet.init.b := [];
    tmp := Maximum(List(sheet.init.l, x->Length(x)));
    for i in [ 1..tmp]  do
        sheet.init.b[i] := [];
        for j in [1..Length(sheet.init.l)] do
             if IsBound(sheet.init.l[j][i]) then
                 Add(sheet.init.b[i], sheet.init.l[j][i]);
             fi;
        od;
    od;
    sheet.init.b := Filtered( sheet.init.b, i -> 0 < Length(i) );
    
    # compute the row <sheet>.r of <sheet>.b and <sheet>.l
    sheet.init.r := [];
    ri  := [];
    for i  in [ 1 .. Length(sheet.init.l) ]  do
        for j  in [ 1 .. Length(sheet.init.b) ]  do
            tmp := Intersection( sheet.init.l[i], sheet.init.b[j] );
            if 0 < Length(tmp) then
                if Length(tmp) = 1 then
                    Add( ri, tmp );
                    Add( sheet.init.r, [ [i,j], tmp ] );
                else
                    for k  in [ 1 .. Length(tmp) ]  do
                        Add( ri, [tmp[k]] );
                        Add( sheet.init.r, [ [i,j], [tmp[k]] ] );
                    od;
                fi;
            fi;
        od;
    od;   
  
    # compute the x coordinates
    sheet.init.x := List( sheet.init.vertices, x -> 0 );
    sheet.init.x[1] := QuoInt( sheet.width, 2 );
#    sheet.init.x[Length(sheet.init.vertices)] := QuoInt( sheet.width, 2 );
   
    # divide x axis
    n := Length(sheet.init.b);
    if n = 2  then
        coors := [ QuoInt(sheet.width,3), 2*QuoInt(sheet.width,3) ];
    else
        coors := [ QuoInt(sheet.width,2) ];
        for j  in [ 1 .. QuoInt(n+1,2)-1 ]  do
            Add( coors, QuoInt( (2*j-1)*sheet.width, 2*n ) );
        od;
        for j  in [ QuoInt(n+1,2)+1 .. n ]  do
            Add( coors, QuoInt( (2*j-1)*sheet.width, 2*n ) );
        od;
    fi;
    coor := [];
    for i  in [ 1 .. Length(sheet.init.b) ]  do
        for j in [ 1 .. Length(sheet.init.b[i]) ]  do
          pos1:=noElements[Position(sheet.init.orders,sheet.init.b[i][j])];
          pos2:=noClasses[Position(sheet.init.orders,sheet.init.b[i][j])][1];
          if 1 = Length(coors)
               and 1 < pos1
               and 2*sheet.circle*pos1 < sheet.width
          then
              coor[ Position(ri,[sheet.init.b[i][j]]) ] :=
                QuoInt(sheet.width,2) - QuoInt(pos1+pos2-1,2)*sheet.circle;
          elif   1 = Length(coors)
             and 1 < pos1
             and pos1 > sheet.width/(2*sheet.circle) 
          then
              coor[ Position(ri,[sheet.init.b[i][j]]) ] :=
                QuoInt(sheet.width,2) - QuoInt(pos1,2) * sheet.circle;
#          elif 1 = pos1
#             and j in [ 2..Length(sheet.init.b[i])-1 ]
#             and noElements[Position(sheet.init.orders,sheet.init.b[i][j+1])]
#                 =  1 
#          then
#              tmp := QuoInt(Maximum(noElements) + 1, 2);
#              coor[ Position(ri,[sheet.init.b[i][j]]) ] := 
#                coors[i] + (-1)^(j+1) * tmp * sheet.circle;
          else
              coor[ Position(ri,[sheet.init.b[i][j]]) ] := coors[i];
          fi;
        od;
    od;
    for i  in [ 1.. Length(sheet.init.r) ]  do
        sheet.operations.MakeXFirstClass( sheet, coor, i );
    od;
    
    # set x coordinates of all other elements
    i := 1;
    while i <= Length(sheet.init.vertices)  do
        if sheet.init.x[i] = 0  then
            if sheet.init.vertices[i-1][2] <= sheet.init.vertices[i][2]  then
                j := 1;
                while noClasses[j][2] <> sheet.init.vertices[i][2]  do
                    j := j + 1;
                od;
                if 2 * noElements[j] * sheet.circle > sheet.width  then
                    sheet.init.x[i] := sheet.init.x[i-1] + sheet.circle;  
                else
                    sheet.init.x[i] := sheet.init.x[i-1] + 2*sheet.circle;
                fi;
                if     sheet.init.x[i] > sheet.width
                   and sheet.init.vertices[i][1] in sheet.init.reps
                then
                    sheet.init.x[i] := sheet.circle;
                fi;
                cl := sheet.init.classLengths[sheet.init.vertices[i][4]];
                i  := i + 1;
                j  := 1;
                x  := sheet.init.x[i-1]; 
                while j < cl  do
                    sheet.init.x[i] := x + sheet.circle * j;
                    j := j + 1;
                    i := i + 1;
                od;
            else
                Error();
            fi;
        else
            i := i + 1;
        fi;
    od;
end;    


#############################################################################
##
#F  WyckoffLatticeOps.MakeMenus( <sheet> )  . . . . menus for a lattice sheet
##
WyckoffLatticeOps.MakeMenus := function( sheet )
    local   tmp,  menu;

    # <updateMenus0> menus not available if no selection is made
    sheet.updateMenus0 := [];

    # <updateMenus1> menus not available if less than 2 selections are made
    sheet.updateMenus1 := [];

    # <updateMenus2> menus not available if less than 3 selections are made
    sheet.updateMenus2 := [];

    # <updateMenusEq1> menus available if 1 selection is made
    sheet.updateMenusEq1 := [];

    # create the resize menu
    Menu( sheet, "Resize",
        [ "Double Lattice",       sheet.operations.RMDoubleLattice,
          "Halve Lattice",        sheet.operations.RMHalveLattice,
          "Resize Lattice",       sheet.operations.RMResizeLattice,
          ,                       Ignore,
          "Resize Graphic Sheet", sheet.operations.RMResizeGraphicSheet
        ] );

    # create the clean up menu
    if sheet.color.model = "color"  then
        tmp := [
            "Average Y Levels",       	sheet.operations.CMAverageYLevels,
            "Average X Levels",         sheet.operations.CMAverageXLevels,
            "Rotate Conjugates",        sheet.operations.CMRotateConjugates,
            ,			        Ignore,
            "Deselect All",             sheet.operations.CMDeselectAll,
#            "Select Representatives",   sheet.operations.CMSelectReps,
            "Relabel Vertices",         sheet.operations.CMRelabelVertices,
            ,                           Ignore,
            "Use Black&White",          sheet.operations.CMUseBlackWhite ];
    else
        tmp := [
            "Average Y Levels",         sheet.operations.CMAverageYLevels,
            "Average X Levels",         sheet.operations.CMAverageXLevels,
            "Rotate Conjugates",        sheet.operations.CMRotateConjugates,
            ,			        Ignore,
            "Deselect All",             sheet.operations.CMDeselectAll,
#            "Select Representatives",   sheet.operations.CMSelectReps ];
            "Relabel Vertices",         sheet.operations.CMRelabelVertices ];
    fi;
    sheet.cleanUpMenu := Menu( sheet, "CleanUp", tmp );
    Add( sheet.updateMenus0, [ sheet.cleanUpMenu, [
         "Deselect All" ] ] );
    Add( sheet.updateMenus1, [ sheet.cleanUpMenu, [
         "Average X Levels", "Rotate Conjugates" ] ] );

end;


#############################################################################
##
#F  WyckoffLatticeOps.SubgroupVertex( <sheet>, <ver> )   subgroup of a vertex
##
WyckoffLatticeOps.SubgroupVertex := function( sheet, ver )
    local p;

    if not IsBound(ver.group)  then
        p := Position( sheet.lattice.classPositions, ver.ident );
        ver.group := sheet.lattice.vertices[p];
    fi;
    return ver.group;   # this isn't a group, but a Wyckoff position
end;


#############################################################################
##
#F  WyckoffLatticeOps.SIStabDim( <sheet>, <ver>, <grp> )    dim. of stabspace
##
WyckoffLatticeOps.SIStabDim := function( sheet, ver, grp )

    local p;

    if not IsBound(ver.info.stabdim) then
        p := Position( sheet.lattice.classPositions, ver.ident );
        ver.info.stabdim := Length(sheet.lattice.vertices[p].basis);
    fi;

    if ver.info.stabdim < 20  then
        return String(ver.info.stabdim);
    else
        return StringPP(ver.info.stabdim);
    fi;
end;


#############################################################################
##
#F  WyckoffLatticeOps.SIClassSize( <sheet>, <ver>, <grp> ) size of Wyckoff cl.
##
WyckoffLatticeOps.SIClassSize := function( sheet, ver, grp )

    if not IsBound(ver.info.index) then
        ver.info.index := sheet.pointGroupSize / Size( grp );
    fi;

    if ver.info.index < 20  then
        return String(ver.info.index);
    else
        return StringPP(ver.info.index);
    fi;
end;


#############################################################################
##
#F  WyckoffLatticeOps.SIConjugacyClassInfo( <sheet>, <ver>, <grp> )   
##
WyckoffLatticeOps.SIConjugacyClassInfo := function( sheet, ver, grp )

    local c,r,l;

    if not IsBound(ver.info.ccInfo) then
        c := ConjugacyClasses(grp);
        r := List(c,Representative);
        l := [ [1..Length(c)],
               List(r,x->Order(grp,x)),
               List(r,x->TraceMat(x)-1),
               List(r,DeterminantMat),
               List(c,Size)
             ];
        l:=TransposedMat(l);
        l:=Concatenation([["class","order","trace","det","size"]],l);
        ver.info.ccInfo := l;
    fi;
    Print("#I Conjugacy Class Information About Vertex ",
              ver.label.text, " \n");
    PrintArray(ver.info.ccInfo);

    return "done";

end;


#############################################################################
##
#F  WyckoffLatticeOps.SITranslation( <sheet>, <ver>, <grp> )   
##
WyckoffLatticeOps.SITranslation := function( sheet, ver, grp )

    local p;

    if not IsBound(ver.info.translation)  then
        p := Position( sheet.lattice.classPositions, ver.ident );
        ver.info.translation := sheet.lattice.vertices[p].translation;
    fi;

    Print("#I Translation of Vertex ", ver.label.text, " \n");
    PrintArray(ver.info.translation);

    return "done";

end;


#############################################################################
##
#F  WyckoffLatticeOps.SIBasis( <sheet>, <ver>, <grp> )   
##
WyckoffLatticeOps.SIBasis := function( sheet, ver, grp )

    local p;

    if not IsBound(ver.info.basis)  then
        p := Position( sheet.lattice.classPositions, ver.ident );
        ver.info.basis := sheet.lattice.vertices[p].basis;
    fi;

    Print("#I Basis of Vertex ", ver.label.text, " \n");
    PrintArray(ver.info.basis);

    return "done";

end;


#############################################################################
##
#F  WyckoffLatticeOps.PMInformation( <sheet>, <obj> ) . . . . show group info
##
WyckoffLatticeOps.PMInformation := function( sheet, obj )
    local   g,  grp,  text,  info,  str,  i,  func1,  func2;
    
    # destroy other text selectors flying around
    if IsBound(sheet.selector)  then
        Close(sheet.selector);
    fi;
   
    # get the group of <obj>
    if not IsBound(obj.group) then
        sheet.operations.SubgroupVertex( sheet, obj );
    fi;
    if not IsBound(obj.grp) then
        g := obj.group.stabilizer;
        obj.grp := Group( g.generators, g.identity );
    fi;
    grp := obj.grp;

    # construct info texts (text, record component, function)
    info := [
      [ "StabDim",     "stabdim",     sheet.operations.SIStabDim     ],
      [ "StabSize",    "size",        sheet.operations.SISize        ],
      [ "ClassSize",   "index",       sheet.operations.SIClassSize   ],
#      [ "Index",       "index",       sheet.operations.SIIndex       ],
      [ "IsAbelian",   "isAbelian",   sheet.operations.SIIsAbelian   ],
#      [ "IsCentral",   "isCentral",   sheet.operations.SIIsCentral   ],
      [ "IsCyclic",    "isCyclic",    sheet.operations.SIIsCyclic    ],
      [ "IsNilpotent", "isNilpotent", sheet.operations.SIIsNilpotent ],
#      [ "IsNormal",    "isNormal",    sheet.operations.SIIsNormal    ],
      [ "IsPerfect",   "isPerfect",   sheet.operations.SIIsPerfect   ],
      [ "IsSimple",    "isSimple",    sheet.operations.SIIsSimple    ],
      [ "IsSolvable",  "isSolvable",  sheet.operations.SIIsSolvable  ],
      [ "Isomorphism", "groupId",     sheet.operations.SIIsomorphism ],
      [ "ConjClassInfo","ccInfo",sheet.operations.SIConjugacyClassInfo],
      [ "Translation", "translation", sheet.operations.SITranslation ],
      [ "Basis",       "basis",       sheet.operations.SIBasis       ]
    ];
   
    # text select function
    func1 := function( sel, tid )
        tid  := sel.selected;
        text := Copy(sel.labels);
        str  := String( info[tid][1], -14 );
        Append( str, info[tid][3]( sheet, obj, grp ) );
        text[tid] := str;
        Relabel( sel, text );
    end;

    # construct the string
    text := [];
    if not IsBound(obj.info)  then
        obj.info := rec();
    fi;
    for i  in info  do
        str := String( i[1], -14 );
        if IsBound(obj.info.(i[2]))  then
            Append( str, i[3]( sheet, obj, grp ) );
        else
            Append( str, "unknown" );
        fi;
        Add( text, str );
        Add( text, func1 );
    od;

    # button select functions
    func1 := function( sel, bt )
        Close(sel);
        Unbind(sheet.selector);
    end;
    func2 := function( sel, bt )
        local   i;
        for i  in [ 1 .. Length(sel.labels) ]  do
            sel.selected := i;
            sel.textFuncs[i]( sel, sel.labels[i] );
        od;
        Enable( sel, "all", false );
    end;

    # construct text selector
    sheet.selector := TextSelector(
        Concatenation( " Information about ", obj.label.text ),
        text,
        [ "all", func2, "close", func1 ] );
               
end;


#############################################################################
##
#F  WyckoffLattice( <S>, <x>, <y> ) . . . . . . . . display a Wyckoff lattice
##
WyckoffLattice := function( arg )
    local   R,  S,  i,  j,  tmp,  def,  match,  permlist,  d;
    
    # we need at least one argument: the space group
    if Length(arg) < 1  then
        Error("usage: WyckoffLattice( <S>[, <x>, <y>][, \"prime\"] )");
    fi;

    if not IsSpaceGroup( arg[1] ) then
        Error("S must be a space group");
    fi;

    # match a substring
    match := function( b, a )
        return a{[1..Length(b)]} = b;
    end;

    # parse the other argument: <x>, <y> or "prime"
    def := rec( prime := false, normal := false, x := 800, y := 600 );
    i := 2;
    while i <= Length(arg)  do
        if IsInt(arg[i])  then
            def.x := arg[i];
            def.y := arg[i+1];
            i := i+2;
        elif IsString(arg[i])  then
            if match( arg[i], "prime" )  then
                def.prime := true;
            elif match( arg[i], "Prime" )  then
                def.prime := true;
            elif match( arg[i], "title" )  then
                def.title := arg[i];
            elif match( arg[i], "Title" )  then
                def.title := arg[i];
            else
                Error( "unkown option \"", arg[i], "\",\n",
                       "options are: \"prime\" or \"title\"" );
            fi;
            i := i+1;
        else
            Error( "unkown option ", arg[i] );
        fi;
    od;

    # construct a nice title
    if not IsBound(def.title)  then
        def.title := "Wyckoff Lattice";
    fi;

    # open a graphic sheet
    S := GraphicSheet( def.title, def.x, def.y );
    S.defaultTitle := def.title;
    S.close := function(S)
        if IsBound(S.selector)  then
            Close(S.selector);
        fi;
    end;

    # <S>.circle is the diameter of a vertex
    S.circle := 2*QuoInt(FONTS.tiny[3]+4*(FONTS.tiny[1]+FONTS.tiny[2])+5,3);

    # <S>.init holds valueable information for the initial setup
    S.init := rec();
    S.init.primeOrdering  := def.prime;

    # change ops for lattice sheet
    S.operations := WyckoffLatticeOps;

    # select a color model
    S.color := rec();
    if COLORS.red <> false or COLORS.lightGray <> false  then
        S.color.model := "color";
    else
        S.color.model := "monochrome";
    fi;
    S.operations.MakeColors(S);
    
    # no objects are selected at first
    S.selected := [];
        
    # compute the lattice
    SetTitle( S, "Computing Wyckoff Lattice" );
    S.lattice := WyckoffLatticeRecord( arg[1] );
    S.operations.MakeLattice(S);

    # store size of the full point group
    d:=Length(arg[1].identity)-1;
    S.pointGroupSize := Size( Group( List( arg[1].generators, 
                          g -> g{[1..d]}{[1..d]} ), IdentityMat(d) ) );
    
    # sort orders according to the number of factors
    if S.init.primeOrdering  then
        Sort( S.init.orders, function(a,b)
              local la, lb;
              la := Length(Factors(a));
              lb := Length(Factors(b));
              if la = lb  then
                  return a < b;
              else
                  return la < lb;
              fi;
        end );
    fi;

    # sort the maximal subgroups
    SetTitle( S, "Computing Maximal Subgroups" );
    S.operations.SortMaximals(S);

    # compute the x-coordinates
    SetTitle( S, "Computing Coordinates" );
    S.operations.MakeX(S);
      
    # compute the y-coordinates
    S.operations.MakeY(S);
    
    # draw lattice
    SetTitle( S, "Drawing" );
    S.operations.MakeVertices(S);

# the following is an unwanted permutation, and moreover leads
# to a problem with some unbound string - so why not cancel it!
#
#    permlist := List(S.init.vertices, x -> x[1]);
#    permlist := Permuted(S.lattice.vertices, Sortex(permlist));
#    for i in [ 1..Length(permlist) ] do
#        Relabel( S.vertices[i], String(permlist[i]));
#    od;
    
    S.operations.MakeConnections(S);

    # and unbind unused information
    # Unbind(S.init);
    
    # add pointer actions to <S>
    InstallGSMethod( S, "LeftPBDown",      S.operations.LeftButton      );
    InstallGSMethod( S, "RightPBDown",     S.operations.RightButton     );
    InstallGSMethod( S, "ShiftLeftPBDown", S.operations.ShiftLeftButton );
    InstallGSMethod( S, "CtrlLeftPBDown",  S.operations.ShiftLeftButton );
    
    # create menus
    S.operations.MakeMenus(S);

    # and enable/disable entries
    S.operations.UpdateMenus(S);

    # reset the title
    SetTitle( S, S.defaultTitle );

    # that's it
    return S;

end;
