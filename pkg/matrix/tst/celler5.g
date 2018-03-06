# timings

if not IsBound(NOS)  then NOS := 1;      fi;

#############################################################################
ReadDataPkg( "matrix", "tst", "nsimple.g" );
r := RecogniseClassical( a7d4q2, "clg", 100 );
if not 7 in PossibleAlternatingGroupsFlag(r)  then
    Error( "failed on a7d4q2" );
fi;
Print( "#I  CLG a7d4q2:       OK\n" );

r := RecogniseClassical( m11d5q3, "clg", 100 );
if not "M11" in PossibleAlmostSimpleFlag(r)  then
    Error( "failed on m11d5q3" );
fi;
Print( "#I  CLG m11d5q3:      OK\n" );

r := RecogniseClassical( 2m12d6q3, "clg", 100 );
if not "M12" in PossibleAlmostSimpleFlag(r)  then
    Error( "failed on 2m12d6q3" );
fi;
Print( "#I  CLG 2m12d6q3:     OK\n" );

r := RecogniseClassical( m23d11q2_1, "clg", 100 );
if not "M23" in PossibleAlmostSimpleFlag(r)  then
    Error( "failed on m23d11q2_1" );
fi;
Print( "#I  CLG m23d11q2_1:   OK\n" );

r := RecogniseClassical( m23d11q2_2, "clg", 100 );
if not "M23" in PossibleAlmostSimpleFlag(r)  then
    Error( "failed on m23d11q2_2" );
fi;
Print( "#I  CLG m23d11q2_2:   OK\n" );

r := RecogniseClassical( m24d11q2_1, "clg", 100 );
if not "M24" in PossibleAlmostSimpleFlag(r)  then
    Error( "failed on m24d11q2_1" );
fi;
Print( "#I  CLG m24d11q2_1:   OK\n" );

r := RecogniseClassical( m24d11q2_2, "clg", 100 );
if not "M24" in PossibleAlmostSimpleFlag(r)  then
    Error( "failed on m24d11q2_2" );
fi;
Print( "#I  CLG m24d11q2_2:   OK\n" );

r := RecogniseClassical( psl27d3q11_1, "clg", 100 );
if not ["A",1,7,1] in PossibleChevalleyGroupsFlag(r)  then
    Error( "failed on psl27d3q11_1" );
fi;
Print( "#I  CLG psl27d3q11_1: OK\n" );

r := RecogniseClassical( psl27d3q11_2, "clg", 100 );
if not ["A",1,7,1] in PossibleChevalleyGroupsFlag(r)  then
    Error( "failed on psl27d3q11_2" );
fi;
Print( "#I  CLG psl27d3q11_2: OK\n" );


#############################################################################
for q  in [ 2, 9 ]  do
    for dd  in [ 1 .. 2 ]  do


#############################################################################
    d := [ 3, 9, 20 ];
    d := d[dd];

        # compute a set of test groups
        grp := GL( d, q );
        gls := [];
        for i  in [ 1 .. NOS ]  do
            repeat
                new := Group( PseudoRandom(grp),
                              PseudoRandom(grp),
                              PseudoRandom(grp) );
            until IsSLContainedFlag(RecogniseClassical(new)) = true;
            Add( gls, new.generators );
        od;
            
        # CLG
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            grp := ApplyFunc( Group, gls[i] );
            x := RecogniseClassical( grp, "clg" );
            if IsSLContainedFlag(x) <> true  then fail := fail + 1;  fi;
        od;
        Print( "#I  CLG GL(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );
            

#############################################################################
    d := [ 3, 9, 21 ];
    d := d[dd];

        # compute a set of test groups
        grp := GU( d, q );
        gls := [];
        for i  in [ 1 .. NOS ]  do
            repeat
                new := Group( PseudoRandom(grp),
                              PseudoRandom(grp),
                              PseudoRandom(grp) );
            until IsUnitaryGroupFlag(RecogniseClassical(new)) = true;
            Add( gls, new.generators );
        od;
            

        # CLG
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            grp := ApplyFunc( Group, gls[i] );
            x := RecogniseClassical( grp, "clg" );
            if IsUnitaryGroupFlag(x) <> true  then fail := fail + 1;  fi;
        od;
        Print( "#I  CLG SU(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );
            
            
#############################################################################
    d := [ 4, 8, 20 ];
    d := d[dd];

        # compute a set of test groups
        grp := SP( d, q );
        gls := [];
        for i  in [ 1 .. NOS ]  do
            repeat
                new := Group( PseudoRandom(grp),
                              PseudoRandom(grp),
                              PseudoRandom(grp) );
            until IsSymplecticGroupFlag(RecogniseClassical(new)) = true;
            Add( gls, new.generators );
        od;
            

        # CLG
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            grp := ApplyFunc( Group, gls[i] );
            x := RecogniseClassical( grp, "clg" );
            if IsSymplecticGroupFlag(x) <> true  then fail := fail + 1;  fi;
        od;
        Print( "#I  CLG SP(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );
            
            
#############################################################################
    d := [ 10, 12, 20 ];
    d := d[dd];

        # compute a set of test groups
        grp := O( +1, d, q );
        gls := [];
        for i  in [ 1 .. NOS ]  do
            repeat
                new := Group( PseudoRandom(grp),
                              PseudoRandom(grp),
                              PseudoRandom(grp) );
            until IsOrthogonalGroupFlag(RecogniseClassical(new)) = true;
            Add( gls, new.generators );
        od;
            
        # CLG
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            grp := ApplyFunc( Group, gls[i] );
            x := RecogniseClassical( grp, "clg" );
            if IsOrthogonalGroupFlag(x) <> true  then
                fail := fail + 1;
            elif ClassicalTypeFlag(x) <> "orthogonalplus"  then
                Error( "in 'ClassicalForms'" );
            fi;
        od;
        Print( "#I  CLG O+(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );
            
            
#############################################################################
    d := [ 10, 12, 30 ];
    d := d[dd];

        # compute a set of test groups
        grp := O( -1, d, q );
        gls := [];
        for i  in [ 1 .. NOS ]  do
            repeat
                new := Group( PseudoRandom(grp),
                              PseudoRandom(grp),
                              PseudoRandom(grp) );
            until IsOrthogonalGroupFlag(RecogniseClassical(new)) = true;
            Add( gls, new.generators );
        od;
            

        # CLG
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            grp := ApplyFunc( Group, gls[i] );
            x := RecogniseClassical( grp, "clg" );
            if IsOrthogonalGroupFlag(x) <> true  then
                fail := fail + 1;
            elif ClassicalTypeFlag(x) <> "orthogonalminus"  then
                Error( "in 'ClassicalForms'" );
            fi;
        od;
        Print( "#I  CLG O-(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );
            
            
#############################################################################
    d := [ 9, 11, 31 ];
    d := d[dd];

    if not q in [ 2, 4, 8, 16 ]  then

        # compute a set of test groups
        grp := O( 0, d, q );
        gls := [];
        for i  in [ 1 .. NOS ]  do
            repeat
                new := Group( PseudoRandom(grp),
                              PseudoRandom(grp),
                              PseudoRandom(grp) );
            until IsOrthogonalGroupFlag(RecogniseClassical(new)) = true;
            Add( gls, new.generators );
        od;
            

        # CLG
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            grp := ApplyFunc( Group, gls[i] );
            x := RecogniseClassical( grp, "clg" );
            if IsOrthogonalGroupFlag(x) <> true  then
                fail := fail + 1;
            elif ClassicalTypeFlag(x) <> "orthogonalcircle"  then
                Error( "in 'ClassicalForms'" );
            fi;
        od;
        Print( "#I  CLG O0(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );

    fi;

    od;
od;
