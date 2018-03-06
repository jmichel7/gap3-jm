#############################################################################
##
#A  fonts.g			XGAP Test File	                 Frank Celler
##
#H  $Log: fonts.g,v $
#H  Revision 1.1.1.1  1996/12/11 12:40:11  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 1.4  1995/08/16  12:53:39  fceller
#H  changed parameters of 'Line' slightly
#H
#H  Revision 1.3  1995/08/10  18:14:19  fceller
#H  changed 'Box', 'Rectangle', 'Diamond' and 'Line'
#H  to expect the start position and a width and height
#H
#H  Revision 1.2  1995/08/09  10:58:35  fceller
#H  created a box-text to check if the font is monospaced
#H
#H  Revision 1.1  1995/07/24  09:40:08  fceller
#H  Initial revision
##


# create a new graphic sheet
G := GraphicSheet( "FONT TEST", 450, 450 );


# give some info to the user
l := PrintGS( G, 10, 20, Concatenation(
    "This test will create a new window showing the five ",
    "possible font sizes which can be used in graphic sheets." ) );


# start position
x := 10;
y := 20 + (FONTS.normal[1]+FONTS.normal[2]+2) * (l+1);


# tiny font
h := FONTS.tiny[1]+FONTS.tiny[2]+2;
Text( G, FONTS.tiny, x,   y, "Tiny Font ---- TINY FONT" );
Text( G, FONTS.tiny, x, h+y, "Test Me  TEST ME Test Me" );
y := y + 2*h + 10;


# small font
h := FONTS.small[1]+FONTS.small[2]+2;
Text( G, FONTS.small, x,   y, "Small Font ---- SMALL FONT" );
Text( G, FONTS.small, x, h+y, "Test Me TEST ME Test ME TE" );
y := y + 2*h + 10;


# normal font
h := FONTS.normal[1]+FONTS.normal[2]+2;
Text( G, FONTS.normal, x,   y, "Normal Font ---- NORMAL FONT" );
Text( G, FONTS.normal, x, h+y, "Test Me TEST ME Test ME TEST" );
y := y + 2*h + 10;


# large font
h := FONTS.large[1]+FONTS.large[2]+2;
Text( G, FONTS.large,  x,   y, "Large Font ---- LARGE FONT" );
Text( G, FONTS.large,  x, h+y, "Test Me TEST ME Test ME TE" );
y := y + 2*h + 10;


# huge font
h := FONTS.huge[1]+FONTS.huge[2]+2;
Text( G, FONTS.huge,   x,   y, "Huge Font ---- HUGE FONT" );
Text( G, FONTS.huge,   x, h+y, "Test Me TEST  ME Test ME" );
y := y + 2*h + 10;

# font dimension test
def := rec();
if COLORS.green <> false  then
  def.color := COLORS.green;
elif COLORS.lightGray <> false  then
  def.color := COLORS.lightGray;
fi;

for font  in [ "tiny", "small", "normal", "large", "huge" ]  do
  t := "MI!*aghjypqwm";
  w := Length(t)*FONTS.(font)[3];
  z := FONTS.(font)[2];
  c := -FONTS.(font)[1];
  Line( G, x, y,   w, 0, def );
  Line( G, x, y+z, w, 0, def );
  Line( G, x, y+c, w, 0, def );
  for i  in [ 0 .. Length(t) ]  do
      Line( G, x+i*FONTS.(font)[3], y+z, 0, c-z, def );
  od;
  Text( G, FONTS.(font), x, y, t );
  y := y + FONTS.(font)[1]+FONTS.(font)[2] + 12;
od;

