#############################################################################
##
#W  matrix.g	 	   	Matrix Packages                  Frank Celler
##
#H  @(#)$Id: matrix.g,v 1.1 1997/03/10 13:52:04 gap Exp $
##
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This file contains utility functions for matrices.
##
RevisionMatrix.matrix_g :=
    "@(#)$Id: matrix.g,v 1.1 1997/03/10 13:52:04 gap Exp $";


#############################################################################
##

#F  DisplayMat( <field>, <mat> )  . . . . . . . . . . . . pretty print matrix
##
DisplayMat := function( arg )
    local   m,  f,  p,  w,  t,  x,  v,  z,  y;

    # get field
    if Length(arg) = 1  then
        m := arg[1];
        if IsMatrix(m)  then
            f := Field( List( m, x -> Field(x).root ) );
        else
            f := Set( Flat( m ) );
            f := Field( f );
        fi;
    else
        m := arg[2];
        f := arg[1];
    fi;
    if IsList(m) and 0 < Length(m) and IsMatrix(m[1])  then
        for x  in m  do
            DisplayMat( f, x );
            Print( "\n" );
        od;
        return;
    fi;

    # if it is a finite prime field,  use integers for display
    if IsFinite(f) and 1 = f.degree  then

        # compute maximal width
        p := f.char;
        w := LogInt( p, 10 ) + 2;

        # create strings
        t := [];
        for x  in [ 2 .. p ]  do
            t[x] := String( x-1, w );
        od;
        t[1] := String( ".", w );

        # print matrix
        for v  in m  do
            for x  in IntVecFFE(v)  do
                Print( t[x+1] );
            od;
            Print( "\n" );
        od;

    # if it a finite,  use mixed integers/z notation
    elif  IsFinite(f)  then

        # compute maximal width
        w := LogInt( Size(f)-1, 10 ) + 4;

        # create strings
        t := [];
        p := Elements(GF(f.char));
        z := f.root;
        for x  in [ 0 .. Size(f)-2 ]  do
            y := z^x;
            if y in p  then
                t[x+2] := String( IntFFE(y), w );
            else
                t[x+2] := String( Concatenation( "z^", String(x) ), w );
            fi;
        od;
        t[1] := String( ".", w );

        # print matrix
        for v  in m  do
            for x  in v  do
                if x = f.zero  then
                    Print( t[1] );
                else
                    Print( t[LogFFE(x,f.root)+2] );
                fi;
            od;
            Print( "\n" );
        od;
        #Print( "\nz := Z(", Size(f), ")\n" );

    # sorry not implement
    else
        Error( "sorry,  pretty print for field ", f, " not implemented" );
    fi;

end;


#############################################################################
##

#E  matrix.g  . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
##
