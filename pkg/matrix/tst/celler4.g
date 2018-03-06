if not IsBound(NOS)  then NOS := 2;      fi;

Print( "#I  number of test runs: ", NOS, "\n" );
for q  in [ 2, 3, 4, 5, 8, 9 ]  do
    for dd  in [ 1 .. 5 ]  do

#############################################################################
    d := [ 4, 6, 8, 10, 20 ];
    d := d[dd];

        # compute a set of test groups
        grp := SP( d, q );
        gls := [];
        for i  in [ 1 .. NOS ]  do
            Add(gls,[PseudoRandom(grp),PseudoRandom(grp),PseudoRandom(grp)]);
        od;

        # CLG
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            tmp  := RandomInvertibleMat( d, GF(q) );
            grp  := ApplyFunc( Group, OnTuples( gls[i], tmp ) );
            form := ClassicalForms(grp);
            if 1 <> Length(form)  then
                fail := fail+1;
            elif form[1][1] <> "symplectic"  then
                if form[1][1] <> "unknown"  then
                    Error( "illigal form found" );
                fi;
                fail := fail+1;
            fi;
        od;
        Print( "#I  SP(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );

#############################################################################
    d := [ 4, 5, 7, 10, 19 ];
    d := d[dd];

        # compute a set of test groups
        grp := SU( d, q );
        gls := [];
        for i  in [ 1 .. NOS ]  do
            Add(gls,[PseudoRandom(grp),PseudoRandom(grp),PseudoRandom(grp)]);
        od;

        # CLG
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            tmp  := RandomInvertibleMat( d, GF(q) );
            grp  := ApplyFunc( Group, OnTuples( gls[i], tmp ) );
            form := ClassicalForms(grp);
            if 1 <> Length(form)  then
                fail := fail+1;
            elif form[1][1] <> "unitary"  then
                if form[1][1] <> "unknown"  then
                    Error( "illigal form found" );
                fi;
                fail := fail+1;
            fi;
        od;
        Print( "#I  SU(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );

#############################################################################
    d := [ 4, 6, 8, 10, 12 ];
    d := d[dd];

        # compute a set of test groups
        grp := O( +1, d, q );
        gls := [];
        for i  in [ 1 .. NOS ]  do
            Add(gls,[PseudoRandom(grp),PseudoRandom(grp),PseudoRandom(grp)]);
        od;

        # CLG
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            tmp  := RandomInvertibleMat( d, GF(q) );
            grp  := ApplyFunc( Group, OnTuples( gls[i], tmp ) );
            form := ClassicalForms(grp);
            if 1 <> Length(form)  then
                fail := fail+1;
            elif form[1][1] <> "orthogonalplus"  then
                if form[1][1] <> "unknown"  then
                    Error( "illigal form found" );
                fi;
                fail := fail+1;
            fi;
        od;
        Print( "#I  O+(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );

#############################################################################
    d := [ 4, 6, 8, 10, 12 ];
    d := d[dd];

        # compute a set of test groups
        grp := O( -1, d, q );
        gls := [];
        for i  in [ 1 .. NOS ]  do
            Add(gls,[PseudoRandom(grp),PseudoRandom(grp),PseudoRandom(grp)]);
        od;

        # CLG
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            tmp  := RandomInvertibleMat( d, GF(q) );
            grp  := ApplyFunc( Group, OnTuples( gls[i], tmp ) );
            form := ClassicalForms(grp);
            if 1 <> Length(form)  then
                fail := fail+1;
            elif form[1][1] <> "orthogonalminus"  then
                if form[1][1] <> "unknown"  then
                    Error( "illigal form found" );
                fi;
                fail := fail+1;
            fi;
        od;
        Print( "#I  O-(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );

#############################################################################
    d := [ 4, 6, 8, 10, 12 ] + 1;
    d := d[dd];
    if not q in [ 2, 4, 8, 16, 32 ]  then

        # compute a set of test groups
        grp := O( 0, d, q );
        gls := [];
        for i  in [ 1 .. NOS ]  do
            Add(gls,[PseudoRandom(grp),PseudoRandom(grp),PseudoRandom(grp)]);
        od;

        # CLG
        start := Runtime();
        fail  := 0;
        for n  in [ 1 .. NOS ]  do
            tmp  := RandomInvertibleMat( d, GF(q) );
            grp  := ApplyFunc( Group, OnTuples( gls[i], tmp ) );
            form := ClassicalForms(grp);
            if 1 <> Length(form)  then
                fail := fail+1;
            elif form[1][1] <> "orthogonalcircle"  then
                if form[1][1] <> "unknown"  then
                    Error( "illigal form found" );
                fi;
                fail := fail+1;
            fi;
        od;
        Print( "#I  O0(", String(d,2), ",", q, ") : ",
               String(QuoInt(Runtime()-start,NOS),11), " msec with ",
               String(fail,3), " failures\n" );
    fi;

    od;
od;
