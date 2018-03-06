#############################################################################
##
#A  lines.g			XGAP Test File	                 Frank Celler
##
#H  $Log: lines.g,v $
#H  Revision 1.1.1.1  1996/12/11 12:40:11  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 1.3  1995/08/10  18:08:33  fceller
#H  changed 'Box', 'Rectangle', 'Diamond' and 'Line'
#H  to expect the start position and a width and height
#H
#H  Revision 1.2  1995/07/24  09:52:14  fceller
#H  changed 'Delete' to use only one argument
#H
#H  Revision 1.1  1995/07/24  09:40:08  fceller
#H  Initial revision
##


# create a new graphic sheet
s := GraphicSheet( "LINE TEST", 400, 400 );


# give some info to the user
PrintGS( s, 10, 20, Concatenation(
    "This test will create some lines which will move around the graphic ",
    "sheet changing color.  Press CTRL-C to stop." ) );


# check which colors are available
col := [];
for name  in [ "black", "dimGrey", "lightGrey", "red", "blue", "green" ]  do
    if COLORS.(name) <> false  then
        Add( col, COLORS.(name) );
    fi;
od;
m   := Length(col);
cc  := 0;

# start with one a from (50,0) to (0,50)
x := 50;  y := 0;  w := 50;  h := 50;

# the general direction is down/right
rx := 1;  ry := 1;  dx := 5;  dy := 5;


# compute the next line
NextLine := function()

    # change direction
    dx := Random( [ 2 ..  8 ] );
    dy := Random( [ 4 .. 12 ] );

    # change the length of our line
    w := w + Random( [-1..1] );
    h := h + Random( [-1..1] );

    # check that it doesn't get to long or short
    if w < 40  then w := 40;  fi;
    if h < 40  then h := 40;  fi;
    if 90 < w  then w := 90;  fi;
    if 90 < h  then h := 90;  fi;

    # make sure we don't leave our sheet
    if rx*dx+x < 0          then rx :=  1;  fi;
    if ry*dy+y < 0          then ry :=  1;  fi;
    if s.width <= rx*dx+x   then rx := -1;  fi;
    if s.height <= ry*dy+y  then ry := -1;  fi;

    # ok, draw it
    x := x + rx*dx;  y := y + ry*dx;
    return Line( s, x, y, -rx*w, ry*h, rec(color:=col[1+QuoInt(cc,10)]) );
end;


l := List( [ 1 .. 20 ], x -> NextLine() );
i := Length(l);
while true  do
    i  := (i mod Length(l))+1;
    cc := (cc+1) mod (10*m);
    Delete( l[i] );
    l[i] := NextLine();
od;
