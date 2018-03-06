for d  in 2 * [ 2 .. 15 ]  do
    for q  in [ 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 17, 19, 23 ]  do
        sp := SP( d, q );
        Print( "#I  SP(", d, ",", q, "): " );
        for x  in [ 1 .. 5 ]  do
            repeat
                g  := Group( PseudoRandom(sp), PseudoRandom(sp) );
                cl := RecogniseClassical( g );
            until IsSymplecticGroupFlag(cl) = true;
            form := InvariantFormFlag(cl);
            repeat
                cl := CRecognizeSP( g, g.generators, form );
            until cl.isSP;
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
