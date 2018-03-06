#############################################################################
##
#A  color.g			XGAP Test File	                 Frank Celler
##
#H  $Log: color.g,v $
#H  Revision 1.1.1.1  1996/12/11 12:40:11  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 1.2  1995/08/10  18:06:34  fceller
#H  changed 'Box', 'Rectangle', 'Diamond' and 'Line'
#H  to expect the start position and a width and height
#H
#H  Revision 1.1  1995/07/24  09:40:08  fceller
#H  Initial revision
##

# create a new graphic sheet
G := GraphicSheet( "COLOR ME", 400, 370 );

# give some info to the user
PrintGS( G, 5, 20, Concatenation( 
    "This test will create various colored vertices, boxes, ",
    "lines. The color model used is '", COLORS.model, "'." ) );

# start position
x := 200;
y :=  30;

# loop over the various colors
for color  in [ "black", "dimGray", "lightGray", "red", "blue", "green" ]  do

    y := y + 50;
    if COLORS.(color) <> false  then
        Vertex( G, x, y, rec( label := color, color := COLORS.(color) ) );
        Box( G, 10, y-10, x-40, 20, rec( color := COLORS.(color) ) );
	Line( G, x+30, y-10, 360-x, 0, rec( color := COLORS.(color) ) );
	Line( G, x+30, y,    360-x, 0, rec( color := COLORS.(color) ) );
	Line( G, x+30, y+10, 360-x, 0, rec( color := COLORS.(color) ) );
    #else
    #    Vertex( G, x, y, rec( label := "BLCK", color := COLORS.black ) );
    fi;

od;

