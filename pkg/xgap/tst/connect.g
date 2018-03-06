#############################################################################
##
#A  connect.g			XGAP Test File	                 Frank Celler
##
#H  $Log: connect.g,v $
#H  Revision 1.1.1.1  1996/12/11 12:40:11  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 1.1  1995/07/24  09:40:08  fceller
#H  Initial revision
#H
##


# create a new graphic sheet
G := GraphicSheet( "MOVE ME", 400, 400 );

# give some info to the user
PrintGS( G, 10, 20, Concatenation(
    "This test will create two connected vertices,  the 'MOV' vertex ",
    "can be moved around using the mouse." ) );

# create two connected vertices
V1 := Vertex( G, 200,  50 );  Relabel( V1, "MOV" );
V2 := Vertex( G, 200, 250 );  Relabel( V2, "FIX" );
C  := Connection( V1, V2 );
Relabel( C, "Test Me" );

# cycle through the shapes
shape := 1;
Reshape( V1, shape );

# install a drag as the pointer down action
InstallGSMethod( G, "LeftPBDown", function( G, x, y )
    shape := shape mod 7 + 1;
    Reshape( V1, shape );
    Drag( G, x, y, BUTTONS.left, function( x, y )
        FastUpdate(G); Move( V1, x, y ); FastUpdate(G,false); end );
end );

    
