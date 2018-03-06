# timings
RecogniseClassicalNP;
InfoRecog1 := Ignore;
InfoRecog2 := Ignore;

if not IsBound(NOS)  then NOS := 2;      fi;
if not IsBound(ALL)  then ALL := false;  fi;

Print( "#I  number of test runs: ", NOS, "\n" );
for q  in [ 2, 3, 5, 4, 9 ]  do
    for dd  in [ 1 .. 10 ]  do

#############################################################################
    d := [ 3, 4, 9, 10, 11, 19, 30, 31, 60, 80 ];
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
            

        # NP
if ALL then
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            grp := ApplyFunc( Group, gls[i] );
            x := RecogniseClassicalNP( grp );
            if x <> true  then fail := fail + 1;  fi;
        od;
        Print( "#I  NP  GL(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );
fi;

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
if ALL  then
        Print( "\n" );
fi;
            

#############################################################################
    d := [ 3, 4, 9, 10, 11, 20, 29, 31, 60, 80 ];
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
            

        # NP
if ALL then
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            grp := ApplyFunc( Group, gls[i] );
            x := RecogniseClassicalNP( grp );
            if x <> true  then fail := fail + 1;  fi;
        od;
        Print( "#I  NP  SU(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );
fi;

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
if ALL  then
        Print( "\n" );
fi;
            
            
#############################################################################
    d := [ 4, 8, 10, 12, 20, 28, 30, 32, 60, 80 ];
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
            

        # NP
if ALL then
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            grp := ApplyFunc( Group, gls[i] );
            x := RecogniseClassicalNP( grp );
            if x <> true  then fail := fail + 1;  fi;
        od;
        Print( "#I  NP  SP(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );
fi;

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
if ALL  then
        Print( "\n" );
fi;
            
            
#############################################################################
    d := [ 10, 12, 14, 16, 18, 20, 22, 26, 28, 30 ];
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
            

        # NP
if ALL then
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            grp := ApplyFunc( Group, gls[i] );
            x := RecogniseClassicalNP( grp );
            if x <> true  then fail := fail + 1;  fi;
        od;
        Print( "#I  NP  O+(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );
fi;

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
if ALL  then
        Print( "\n" );
fi;
            
            
#############################################################################
    d := [ 10, 12, 14, 16, 18, 20, 22, 26, 28, 30 ];
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
            

        # NP
if ALL then
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            grp := ApplyFunc( Group, gls[i] );
            x := RecogniseClassicalNP( grp );
            if x <> true  then fail := fail + 1;  fi;
        od;
        Print( "#I  NP  O-(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );
fi;

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
if ALL  then
        Print( "\n" );
fi;
            
            
#############################################################################
    d := [ 9, 11, 13, 15, 17, 19, 21, 27, 29, 31 ];
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
            

        # NP
if ALL then
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            grp := ApplyFunc( Group, gls[i] );
            x := RecogniseClassicalNP( grp );
            if x <> true  then fail := fail + 1;  fi;
        od;
        Print( "#I  NP  O0(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );
fi;

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

if ALL  then
        Print( "\n" );
fi;
    fi;

    od;
od;
