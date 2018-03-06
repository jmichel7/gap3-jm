#############################################################################
##
#A  objects.g			XGAP Test File	                 Frank Celler
##
#H  $Log: objects.g,v $
#H  Revision 1.1.1.1  1996/12/11 12:40:12  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 1.1  1995/08/10  18:06:34  fceller
#H  Initial revision
#H
##


# create a new graphic sheet
sheet := GraphicSheet( "GRAPHIC OBJECTS", 500, 500 );
sheet.new := [];


# give some info to the user
l := PrintGS( sheet, 10, 20, Concatenation(
    "Press any of the buttons below.  Press the left mouse button on the ",
    "empty space to create a such a graphic object.  Use the right ",
    "mouse button to move it" ) );


# highlight a button
HighlightButton := function( sheet, but )
    local   b;

    # unhighlight the other buttons
    for b  in sheet.buttons  do
        SetWidth( b.rectangle, 1 );
    od;

    # highlight <but>
    SetWidth( but.rectangle, 2 );

end;


# create a rectangle
CreateRectangle := function( sheet, x, y )
    local   r;

    r := Rectangle( sheet, x, y, 1, 1 );
    Drag( sheet, x, y, BUTTONS.left, function( nx, ny )
        Reshape( r, nx-x, ny-y );
    end );
    Add( sheet.new, r );

end;


# rectangle button: highlight buttton
RectangleButton := function( sheet, but )
    HighlightButton( sheet, but );
    sheet.create := CreateRectangle;
end;


# create a circle
CreateCircle := function( sheet, x, y )
    local   r;

    r := Circle( sheet, x, y, 1 );
    Drag( sheet, x, y, BUTTONS.left, function( nx, ny )
        Reshape( r, Maximum( nx-x, ny-y ) );
    end );
    Add( sheet.new, r );

end;


# circle button: highlight buttton
CircleButton := function( sheet, but )
    HighlightButton( sheet, but );
    sheet.create := CreateCircle;
end;


# create a disc
CreateDisc := function( sheet, x, y )
    local   r;

    r := Disc( sheet, x, y, 1 );
    Drag( sheet, x, y, BUTTONS.left, function( nx, ny )
        Reshape( r, Maximum( nx-x, ny-y ) );
    end );
    Add( sheet.new, r );

end;


# disc button: highlight buttton
DiscButton := function( sheet, but )
    HighlightButton( sheet, but );
    sheet.create := CreateDisc;
end;


# create a box
CreateBox := function( sheet, x, y )
    local   r;

    r := Box( sheet, x, y, 1, 1 );
    Drag( sheet, x, y, BUTTONS.left, function( nx, ny )
        Reshape( r, nx-x, ny-y );
    end );
    Add( sheet.new, r );

end;


# box button: highlight buttton
BoxButton := function( sheet, but )
    HighlightButton( sheet, but );
    sheet.create := CreateBox;
end;


# create a line
CreateLine := function( sheet, x, y )
    local   r;

    r := Line( sheet, x, y, 1, 1 );
    Drag( sheet, x, y, BUTTONS.left, function( nx, ny )
        Reshape( r, nx-x, ny-y );
    end );
    Add( sheet.new, r );

end;


# line button: highlight buttton
LineButton := function( sheet, but )
    HighlightButton( sheet, but );
    sheet.create := CreateLine;
end;


# create a diamond
CreateDiamond := function( sheet, x, y )
    local   r;

    r := Diamond( sheet, x, y, 1, 1 );
    Drag( sheet, x, y, BUTTONS.left, function( nx, ny )
        Reshape( r, nx-x, ny-y );
    end );
    Add( sheet.new, r );

end;


# diamond button: highlight buttton
DiamondButton := function( sheet, but )
    HighlightButton( sheet, but );
    sheet.create := CreateDiamond;
end;


# create a text
CreateText := function( sheet, x, y )
    local   r;

    r := Text( sheet, FONTS.large, x, y, "HALLO" );
    Add( sheet.new, r );

end;


# text button: highlight buttton
TextButton := function( sheet, but )
    HighlightButton( sheet, but );
    sheet.create := CreateText;
end;


# create the buttons
x := 10;
y := 20 + (FONTS.normal[1]+FONTS.normal[2]+2) * (l+1);

b := Button( sheet, x, y, "rectangle", RectangleButton );
x := x + b.rectangle.w + 5;

b := Button( sheet, x, y, "circle", CircleButton );
x := x + b.rectangle.w + 5;

b := Button( sheet, x, y, "disc", DiscButton );
x := x + b.rectangle.w + 5;

b := Button( sheet, x, y, "box", BoxButton );
x := x + b.rectangle.w + 5;

b := Button( sheet, x, y, "line", LineButton );
x := x + b.rectangle.w + 5;

b := Button( sheet, x, y, "diamond", DiamondButton );
x := x + b.rectangle.w + 5;

b := Button( sheet, x, y, "text", TextButton );


# left mouse button pressed: check if a button is pressed
low := y + 20;
LPBD := function( sheet, x, y )
    local   pos,  b;

    # check for a button
    pos := [x,y];
    for b  in sheet.buttons  do
        if pos in b.rectangle  then
            return b.func( sheet, b );
        fi;
    od;

    # if <y> is low enough call the create function
    if low < y and IsBound(sheet.create)  then
        return sheet.create( sheet, x, y );
    fi;

end;

# delta move a object
RPBD := function( sheet, x, y )
    local   pos,  o;

    # check if the click lies inside an object
    pos := [x,y];
    for o  in sheet.new  do
        if pos in o  then
            Drag( sheet, x, y, BUTTONS.right, function( nx, ny )
                MoveDelta( o, nx-x, ny-y );
                x := nx;
                y := ny;
            end );
            return;
        fi;
    od;

end;


# move a object
SRPBD := function( sheet, x, y )
    local   pos,  o;

    # check if the click lies inside an object
    pos := [x,y];
    for o  in sheet.new  do
        if pos in o  then
            Drag( sheet, x, y, BUTTONS.right, function( nx, ny )
                Move( o, nx, ny );
            end );
            return;
        fi;
    od;

end;


# change the line width and color
colors := Filtered( ["black","red","blue","green","lightGrey","dimGrey"],
                    x -> COLORS.(x) <> false );

SLPBD := function( sheet, x, y )
    local   pos,  o;

    # check if the click lies inside an object
    pos := [x,y];
    for o  in sheet.new  do
        if pos in o  then
            if IsBound(o.width)  then
                SetWidth( o, 3-o.width );
            fi;
            Recolor( o, COLORS.(Random(colors)) );
            return;
        fi;
    od;

end;


# print the info text
CLPBD := function( sheet, x, y )
    local   pos,  o;

    # check if the click lies inside an object
    pos := [x,y];
    for o  in sheet.new  do
        if pos in o  then
            Print( o, "\n" );
            o.operations.PrintInfo(o);
            return;
        fi;
    od;

end;


# install a function for "mouse button down" events
InstallGSMethod( sheet, "RightPBDown", RPBD );
InstallGSMethod( sheet, "ShiftRightPBDown", SRPBD );
InstallGSMethod( sheet, "LeftPBDown", LPBD );
InstallGSMethod( sheet, "ShiftLeftPBDown", SLPBD );
InstallGSMethod( sheet, "CtrlLeftPBDown", CLPBD );

