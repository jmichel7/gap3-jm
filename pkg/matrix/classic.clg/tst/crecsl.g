for d  in [ 2 .. 30 ]  do
    for q  in [ 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23 ]  do
        gl := GL( d, q );
        Print( "#I  SL(", d, ",", q, "): " );
        for x  in [ 1 .. 5 ]  do
            repeat
                g  := Group( PseudoRandom(gl), PseudoRandom(gl) );
                cl := RecogniseClassical( g );
            until IsSLContainedFlag(cl) = true;
            repeat
                cl := CRecognizeSL( g, g.generators );
            until cl.containsSL;
            for t  in [ 1 .. 20 ]  do
                x := PseudoRandom(g);
                w := Rewrite( cl, x );
                y := Value( w, g.generators );
                if x <> y  then
                    Error( "mismatch" );
                fi;
            od;
            Print( ".\c" );
        od;
        Print( " PASSED\n" );
    od;
od;
