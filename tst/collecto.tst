#############################################################################
##
#A  collecto.tst                GAP Test                         Frank Celler
##
#A  @(#)$Id: collecto.tst,v 1.1.1.1 1996/12/11 12:44:37 werner Exp $
##
#Y  Copyright 1990-1991,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##  This files tests   the different   collectors in  soluble  and  nilpotent
##  groups.  It does not test divisions,  commutator or inversions. Note that
##  this test may run very long. At needs at least 3 MByte.
##
#H  $Log: collecto.tst,v $
#H  Revision 1.1.1.1  1996/12/11 12:44:37  werner
#H  Preparing 3.4.4 for release
#H
#H  Revision 3.2  1993/06/01  09:50:31  fceller
#H  fixed a few things for GAP 3.2
#H
#H  Revision 3.1  1991/09/13  12:44:26  fceller
#H  Initial Release under RCS
#H
##

# Define a few abstract generators.
g1 := AbstractGenerator( "g1" );;
g2 := AbstractGenerator( "g2" );;
g3 := AbstractGenerator( "g3" );;
g4 := AbstractGenerator( "g4" );;
g5 := AbstractGenerator( "g5" );;
g6 := AbstractGenerator( "g6" );;
g7 := AbstractGenerator( "g7" );;
g8 := AbstractGenerator( "g8" );;
g9 := AbstractGenerator( "g9" );;
g10 := AbstractGenerator( "g10" );;
g11 := AbstractGenerator( "g11" );;
g12 := AbstractGenerator( "g12" );;
g13 := AbstractGenerator( "g13" );;
g14 := AbstractGenerator( "g14" );;
g15 := AbstractGenerator( "g15" );;
g16 := AbstractGenerator( "g16" );;
g17 := AbstractGenerator( "g17" );;
g18 := AbstractGenerator( "g18" );;
g19 := AbstractGenerator( "g19" );;
g20 := AbstractGenerator( "g20" );;
g21 := AbstractGenerator( "g21" );;
g22 := AbstractGenerator( "g22" );;
g23 := AbstractGenerator( "g23" );;
g24 := AbstractGenerator( "g24" );;
g25 := AbstractGenerator( "g25" );;
g26 := AbstractGenerator( "g26" );;
g27 := AbstractGenerator( "g27" );;
g28 := AbstractGenerator( "g28" );;
g29 := AbstractGenerator( "g29" );;
g30 := AbstractGenerator( "g30" );;
g31 := AbstractGenerator( "g31" );;
g32 := AbstractGenerator( "g32" );;
g33 := AbstractGenerator( "g33" );;
g34 := AbstractGenerator( "g34" );;
g35 := AbstractGenerator( "g35" );;
g36 := AbstractGenerator( "g36" );;
g37 := AbstractGenerator( "g37" );;
g38 := AbstractGenerator( "g38" );;
g39 := AbstractGenerator( "g39" );;
g40 := AbstractGenerator( "g40" );;
g41 := AbstractGenerator( "g41" );;
g42 := AbstractGenerator( "g42" );;
g43 := AbstractGenerator( "g43" );;
g44 := AbstractGenerator( "g44" );;
g45 := AbstractGenerator( "g45" );;
g46 := AbstractGenerator( "g46" );;
g47 := AbstractGenerator( "g47" );;
g48 := AbstractGenerator( "g48" );;
g49 := AbstractGenerator( "g49" );;
g50 := AbstractGenerator( "g50" );;
g51 := AbstractGenerator( "g51" );;
g52 := AbstractGenerator( "g52" );;
g53 := AbstractGenerator( "g53" );;
g54 := AbstractGenerator( "g54" );;
g55 := AbstractGenerator( "g55" );;
g56 := AbstractGenerator( "g56" );;
g57 := AbstractGenerator( "g57" );;
g58 := AbstractGenerator( "g58" );;
g59 := AbstractGenerator( "g59" );;
g60 := AbstractGenerator( "g60" );;
g61 := AbstractGenerator( "g61" );;
g62 := AbstractGenerator( "g62" );;
g63 := AbstractGenerator( "g63" );;
g64 := AbstractGenerator( "g64" );;
g65 := AbstractGenerator( "g65" );;
g66 := AbstractGenerator( "g66" );;
g67 := AbstractGenerator( "g67" );;
g68 := AbstractGenerator( "g68" );;
g69 := AbstractGenerator( "g69" );;
g70 := AbstractGenerator( "g70" );;
g71 := AbstractGenerator( "g71" );;
g72 := AbstractGenerator( "g72" );;
g73 := AbstractGenerator( "g73" );;
g74 := AbstractGenerator( "g74" );;
g75 := AbstractGenerator( "g75" );;
g76 := AbstractGenerator( "g76" );;
g77 := AbstractGenerator( "g77" );;
g78 := AbstractGenerator( "g78" );;
g79 := AbstractGenerator( "g79" );;

# The following groups were computed using Alice's PQ.
GROUP := rec( generators := [
g1,
g2,
g3,
g4,
g5,
g6,
g7,
g8,
g9,
g10,
g11,
g12,
g13,
g14,
g15,
g16,
g17,
g18,
g19,
g20,
g21,
g22,
g23,
g24,
g25,
g26,
g27,
g28,
g29,
g30,
g31,
g32,
g33,
g34,
g35,
g36,
g37,
g38,
g39,
g40,
g41,
g42,
g43,
g44,
g45,
g46,
g47,
g48,
g49,
g50,
g51,
g52,
g53,
g54,
g55,
g56,
g57,
g58,
g59,
g60,
g61,
g62,
g63,
g64,
g65,
g66,
g67,
g68,
g69,
g70,
g71,
g72,
g73,
g74,
g75,
g76,
g77,
g78,
g79],
relations := [
g1^3 / ( g7 ),
g2^3 / ( g8 ),
g3^3 / ( g9 ),
g4^3 / ( g18^2*g24^2*g42*g50^2*g56*g57*g59^2 ),
g5^3 / ( g19^2*g32^2*g45*g51^2 ),
g6^3 / ( g20^2*g39^2*g48*g52^2 ),
g7^3 / ( g21 ),
g8^3 / ( g22 ),
g9^3 / ( g23 ),
g10^3 / ( g42^2*g56^2 ),
g11^3 / ( g43^2*g57^2 ),
g12^3 / ( g44^2*g58^2 ),
g13^3 / ( g45^2 ),
g14^3 / ( g46^2 ),
g15^3 / ( g47^2 ),
g16^3 / ( g48^2 ),
g17^3 / ( g49^2 ),
g18^3 / ( g50 ),
g19^3 / ( g51 ),
g20^3 / ( g52 ),
g21^3 / ( g53 ),
g22^3 / ( g54 ),
g23^3 / ( g55 ),
g24^3 / ( IdWord ),
g25^3 / ( IdWord ),
g26^3 / ( IdWord ),
g27^3 / ( IdWord ),
g28^3 / ( IdWord ),
g29^3 / ( IdWord ),
g30^3 / ( IdWord ),
g31^3 / ( IdWord ),
g32^3 / ( IdWord ),
g33^3 / ( IdWord ),
g34^3 / ( IdWord ),
g35^3 / ( IdWord ),
g36^3 / ( IdWord ),
g37^3 / ( IdWord ),
g38^3 / ( IdWord ),
g39^3 / ( IdWord ),
g40^3 / ( IdWord ),
g41^3 / ( IdWord ),
g42^3 / ( IdWord ),
g43^3 / ( IdWord ),
g44^3 / ( IdWord ),
g45^3 / ( IdWord ),
g46^3 / ( IdWord ),
g47^3 / ( IdWord ),
g48^3 / ( IdWord ),
g49^3 / ( IdWord ),
g50^3 / ( IdWord ),
g51^3 / ( IdWord ),
g52^3 / ( IdWord ),
g53^3 / ( IdWord ),
g54^3 / ( IdWord ),
g55^3 / ( IdWord ),
g56^3 / ( IdWord ),
g57^3 / ( IdWord ),
g58^3 / ( IdWord ),
g59^3 / ( IdWord ),
g60^3 / ( IdWord ),
g61^3 / ( IdWord ),
g62^3 / ( IdWord ),
g63^3 / ( IdWord ),
g64^3 / ( IdWord ),
g65^3 / ( IdWord ),
g66^3 / ( IdWord ),
g67^3 / ( IdWord ),
g68^3 / ( IdWord ),
g69^3 / ( IdWord ),
g70^3 / ( IdWord ),
g71^3 / ( IdWord ),
g72^3 / ( IdWord ),
g73^3 / ( IdWord ),
g74^3 / ( IdWord ),
g75^3 / ( IdWord ),
g76^3 / ( IdWord ),
g77^3 / ( IdWord ),
g78^3 / ( IdWord ),
g79^3 / ( IdWord ),
Comm( g2, g1 ) / ( g4 ),
Comm( g3, g1 ) / ( g5 ),
Comm( g4, g1 ) / ( g10 ),
Comm( g5, g1 ) / ( g13 ),
Comm( g6, g1 ) / ( g12^2*g14*g26*g28*g29^2*g30^2*g36*g37^2*g44*g61^2
    *g62*g63*g67^2*g68*g69*g72^2 ),
Comm( g8, g1 ) / ( g18^2*g24^2*g27*g42*g43^2*g50^2*g56*g59^2*g60*g65^2
     ),
Comm( g9, g1 ) / ( g19^2*g32^2*g38*g45*g47^2*g51^2 ),
Comm( g10, g1 ) / ( g24 ),
Comm( g11, g1 ) / ( g25*g57^2*g59*g60^2*g65 ),
Comm( g12, g1 ) / ( g29 ),
Comm( g13, g1 ) / ( g32 ),
Comm( g14, g1 ) / ( g26*g29^2*g33*g58*g61^2*g63*g68*g71^2*g72^2*g73^2
    *g77 ),
Comm( g15, g1 ) / ( g34 ),
Comm( g16, g1 ) / ( g28*g30*g35*g61^2*g63^2*g67*g68*g69*g72^2*g74^2
    *g75*g76^2*g78 ),
Comm( g17, g1 ) / ( g31^2*g36^2*g37^2*g70*g73*g76^2*g77^2 ),
Comm( g18, g1 ) / ( g42 ),
Comm( g19, g1 ) / ( g45 ),
Comm( g20, g1 ) / ( g44^2*g46*g58^2*g67 ),
Comm( g22, g1 ) / ( g50^2 ),
Comm( g23, g1 ) / ( g51^2 ),
Comm( g24, g1 ) / ( g56 ),
Comm( g25, g1 ) / ( g59 ),
Comm( g26, g1 ) / ( g62 ),
Comm( g27, g1 ) / ( g65 ),
Comm( g28, g1 ) / ( g68 ),
Comm( g29, g1 ) / ( g71 ),
Comm( g30, g1 ) / ( g74 ),
Comm( g31, g1 ) / ( g77 ),
Comm( g33, g1 ) / ( g58^2*g62^2*g71^2 ),
Comm( g35, g1 ) / ( g61^2*g63^2*g68*g72 ),
Comm( g37, g1 ) / ( g64*g73*g77 ),
Comm( g39, g1 ) / ( g67^2 ),
Comm( g40, g1 ) / ( g70*g76 ),
Comm( g41, g1 ) / ( g79^2 ),
Comm( g3, g2 ) / ( g6 ),
Comm( g4, g2 ) / ( g11 ),
Comm( g5, g2 ) / ( g14 ),
Comm( g6, g2 ) / ( g16 ),
Comm( g7, g2 ) / ( g18 ),
Comm( g9, g2 ) / ( g20^2*g39^2*g41*g48*g49^2*g52^2 ),
Comm( g10, g2 ) / ( g25 ),
Comm( g11, g2 ) / ( g27 ),
Comm( g12, g2 ) / ( g30 ),
Comm( g13, g2 ) / ( g33 ),
Comm( g14, g2 ) / ( g35 ),
Comm( g15, g2 ) / ( g37 ),
Comm( g16, g2 ) / ( g39 ),
Comm( g17, g2 ) / ( g40 ),
Comm( g18, g2 ) / ( g43 ),
Comm( g19, g2 ) / ( g46 ),
Comm( g20, g2 ) / ( g48 ),
Comm( g21, g2 ) / ( g50 ),
Comm( g23, g2 ) / ( g52^2 ),
Comm( g24, g2 ) / ( g57 ),
Comm( g25, g2 ) / ( g60 ),
Comm( g26, g2 ) / ( g63 ),
Comm( g27, g2 ) / ( g66 ),
Comm( g28, g2 ) / ( g69 ),
Comm( g29, g2 ) / ( g72 ),
Comm( g30, g2 ) / ( g75 ),
Comm( g31, g2 ) / ( g78 ),
Comm( g4, g3 ) / ( g12 ),
Comm( g5, g3 ) / ( g15 ),
Comm( g6, g3 ) / ( g17 ),
Comm( g7, g3 ) / ( g19 ),
Comm( g8, g3 ) / ( g20 ),
Comm( g10, g3 ) / ( g26 ),
Comm( g11, g3 ) / ( g28 ),
Comm( g12, g3 ) / ( g31 ),
Comm( g13, g3 ) / ( g34 ),
Comm( g14, g3 ) / ( g36 ),
Comm( g15, g3 ) / ( g38 ),
Comm( g16, g3 ) / ( g40 ),
Comm( g17, g3 ) / ( g41 ),
Comm( g18, g3 ) / ( g44 ),
Comm( g19, g3 ) / ( g47 ),
Comm( g20, g3 ) / ( g49 ),
Comm( g21, g3 ) / ( g51 ),
Comm( g22, g3 ) / ( g52 ),
Comm( g24, g3 ) / ( g58 ),
Comm( g25, g3 ) / ( g61 ),
Comm( g26, g3 ) / ( g64 ),
Comm( g27, g3 ) / ( g67 ),
Comm( g28, g3 ) / ( g70 ),
Comm( g29, g3 ) / ( g73 ),
Comm( g30, g3 ) / ( g76 ),
Comm( g31, g3 ) / ( g79 ),
Comm( g5, g4 ) / ( g26*g29^2*g58^2*g62*g73^2*g77 ),
Comm( g6, g4 ) / ( g28*g30^2*g67^2*g69*g76^2*g78 ),
Comm( g7, g4 ) / ( g42 ),
Comm( g8, g4 ) / ( g43*g57*g66^2 ),
Comm( g9, g4 ) / ( g44*g58*g79^2 ),
Comm( g10, g4 ) / ( g57^2*g59 ),
Comm( g11, g4 ) / ( g60^2*g65 ),
Comm( g12, g4 ) / ( g72^2*g74 ),
Comm( g13, g4 ) / ( g58^2*g62^2*g71^2 ),
Comm( g14, g4 ) / ( g61^2*g63*g68*g72^2 ),
Comm( g15, g4 ) / ( g64*g73*g77 ),
Comm( g16, g4 ) / ( g67^2*g69^2*g75^2 ),
Comm( g17, g4 ) / ( g70*g76*g78 ),
Comm( g6, g5 ) / ( g36*g37^2 ),
Comm( g7, g5 ) / ( g45 ),
Comm( g8, g5 ) / ( g46 ),
Comm( g9, g5 ) / ( g47 ),
Comm( g10, g5 ) / ( g58^2*g62 ),
Comm( g11, g5 ) / ( g61^2*g68 ),
Comm( g12, g5 ) / ( g73^2*g77 ),
Comm( g14, g5 ) / ( g64^2*g73 ),
Comm( g7, g6 ) / ( g44^2*g46 ),
Comm( g8, g6 ) / ( g48 ),
Comm( g9, g6 ) / ( g49 ),
Comm( g10, g6 ) / ( g61^2*g63 ),
Comm( g11, g6 ) / ( g67^2*g69 ),
Comm( g12, g6 ) / ( g76^2*g78 ),
Comm( g8, g7 ) / ( g50^2 ),
Comm( g9, g7 ) / ( g51^2 ),
Comm( g9, g8 ) / ( g52^2 ),
] );;
G3 := AgGroupFpGroup( GROUP );;
w  := Product( List( Cgs(G3), x -> x^(RelativeOrder(x)-1) ) );;
G3 := Group( List( G3.generators, x -> x ^ w ), G3.identity );;

GROUP := rec( generators := [
g1,
g2,
g3,
g4,
g5,
g6,
g7,
g8,
g9,
g10,
g11,
g12,
g13,
g14],
relations := [
g1^7 / ( g8^6*g9^3*g10^4*g11^5*g12^6*g13^4*g14^2 ),
g2^7 / ( IdWord ),
g3^7 / ( IdWord ),
g4^7 / ( IdWord ),
g5^7 / ( IdWord ),
g6^7 / ( IdWord ),
g7^7 / ( IdWord ),
g8^7 / ( IdWord ),
g9^7 / ( IdWord ),
g10^7 / ( IdWord ),
g11^7 / ( IdWord ),
g12^7 / ( IdWord ),
g13^7 / ( IdWord ),
g14^7 / ( IdWord ),
Comm( g2, g1 ) / ( g3 ),
Comm( g3, g1 ) / ( g4 ),
Comm( g4, g1 ) / ( g5 ),
Comm( g5, g1 ) / ( g6 ),
Comm( g6, g1 ) / ( g7 ),
Comm( g7, g1 ) / ( g8 ),
Comm( g8, g1 ) / ( g10 ),
Comm( g9, g1 ) / ( g10^5*g11*g12^3*g13^5*g14 ),
Comm( g10, g1 ) / ( g11 ),
Comm( g11, g1 ) / ( g12 ),
Comm( g12, g1 ) / ( g13 ),
Comm( g13, g1 ) / ( g14 ),
Comm( g3, g2 ) / ( g5^6*g6^3*g7^3*g8*g9^6*g11^6*g12^3*g13^2*g14^4 ),
Comm( g4, g2 ) / ( g6^6*g7*g8^5*g9^6*g10^5*g12^6*g13^4*g14 ),
Comm( g5, g2 ) / ( g7^4*g8^6*g9^4*g10^5*g12^4*g13^2*g14^2 ),
Comm( g6, g2 ) / ( g8^2*g9^2*g10^6*g11^3*g12^5*g13^2*g14^6 ),
Comm( g7, g2 ) / ( g9 ),
Comm( g8, g2 ) / ( g10*g12^4*g13^3*g14^5 ),
Comm( g9, g2 ) / ( g12^5*g13*g14^6 ),
Comm( g10, g2 ) / ( g12^6*g13^3*g14 ),
Comm( g11, g2 ) / ( g13^6*g14^2 ),
Comm( g12, g2 ) / ( g14^4 ),
Comm( g4, g3 ) / ( g7^2*g9^3*g10*g11^2*g12*g13 ),
Comm( g5, g3 ) / ( g8^2*g9^6*g10^3*g11^6*g12^5*g13*g14 ),
Comm( g6, g3 ) / ( g9^6*g10*g11^6*g12^6*g13^2*g14^5 ),
Comm( g7, g3 ) / ( g10^4*g13^3*g14 ),
Comm( g8, g3 ) / ( g11*g12*g13*g14^4 ),
Comm( g9, g3 ) / ( g12^5*g13^5*g14^6 ),
Comm( g10, g3 ) / ( g14^6 ),
Comm( g11, g3 ) / ( g14^2 ),
Comm( g5, g4 ) / ( g9*g10^5*g11^6*g12^4*g13*g14^4 ),
Comm( g6, g4 ) / ( g10^5*g11^4*g12^5*g13^5*g14^4 ),
Comm( g7, g4 ) / ( g11^3*g12^5*g13^5*g14^3 ),
Comm( g8, g4 ) / ( g12*g13*g14^4 ),
Comm( g9, g4 ) / ( g13^5*g14^4 ),
Comm( g10, g4 ) / ( g14^5 ),
Comm( g6, g5 ) / ( g11^2*g12^4*g13^4*g14^4 ),
Comm( g7, g5 ) / ( g12^2*g13^3*g14^5 ),
Comm( g8, g5 ) / ( g13*g14^3 ),
Comm( g9, g5 ) / ( g14 ),
Comm( g7, g6 ) / ( g13*g14^6 ),
Comm( g8, g6 ) / ( g14 )
] );;
G7 := AgGroupFpGroup( GROUP );;

GROUP := rec( generators := [
g1,
g2,
g3,
g4,
g5,
g6,
g7,
g8,
g9,
g10,
g11,
g12,
g13,
g14,
g15,
g16,
g17,
g18,
g19,
g20,
g21,
g22,
g23,
g24,
g25,
g26,
g27,
g28,
g29,
g30,
g31,
g32,
g33,
g34,
g35,
g36,
g37,
g38,
g39,
g40,
g41,
g42,
g43,
g44,
g45,
g46,
g47,
g48,
g49,
g50,
g51,
g52,
g53,
g54,
g55],
relations := [
g1^2 / ( g4 ),
g2^2 / ( g5 ),
g3^2 / ( g6*g8*g11*g14*g16*g20*g21*g25*g28*g30*g34*g37*g38
    *g42*g44*g48*g53 ),
g4^2 / ( g9 ),
g5^2 / ( g10 ),
g6^2 / ( g11*g14*g19*g20*g21*g25*g28*g35*g37*g38*g42*g44
    *g48*g51 ),
g7^2 / ( g12*g15*g20*g22*g23*g26*g29*g37*g43*g46*g49*g52
     ),
g8^2 / ( g16*g25*g28*g42*g51 ),
g9^2 / ( g17 ),
g10^2 / ( g18 ),
g11^2 / ( g19*g25*g33*g34*g35*g42*g48 ),
g12^2 / ( g20*g26*g34*g36*g37*g43*g49 ),
g13^2 / ( g22*g27*g36*g39*g40*g45*g50 ),
g14^2 / ( g28*g42*g48 ),
g15^2 / ( g29*g43*g49 ),
g16^2 / ( g30 ),
g17^2 / ( g31 ),
g18^2 / ( g32 ),
g19^2 / ( g33*g42 ),
g20^2 / ( g34*g43 ),
g21^2 / ( g35*g44 ),
g22^2 / ( g36*g45 ),
g23^2 / ( g38*g46 ),
g24^2 / ( g39*g47 ),
g25^2 / ( g48 ),
g26^2 / ( g49 ),
g27^2 / ( g50 ),
g28^2 / ( g51 ),
g29^2 / ( g52 ),
g30^2 / ( g53 ),
g31^2 / ( g54 ),
g32^2 / ( g55 ),
g33^2 / ( IdWord ),
g34^2 / ( IdWord ),
g35^2 / ( IdWord ),
g36^2 / ( IdWord ),
g37^2 / ( IdWord ),
g38^2 / ( IdWord ),
g39^2 / ( IdWord ),
g40^2 / ( IdWord ),
g41^2 / ( IdWord ),
g42^2 / ( IdWord ),
g43^2 / ( IdWord ),
g44^2 / ( IdWord ),
g45^2 / ( IdWord ),
g46^2 / ( IdWord ),
g47^2 / ( IdWord ),
g48^2 / ( IdWord ),
g49^2 / ( IdWord ),
g50^2 / ( IdWord ),
g51^2 / ( IdWord ),
g52^2 / ( IdWord ),
g53^2 / ( IdWord ),
g54^2 / ( IdWord ),
g55^2 / ( IdWord ),
Comm( g2, g1 ) / ( g3 ),
Comm( g3, g1 ) / ( g6 ),
Comm( g5, g1 ) / ( g6*g7*g8*g11*g14*g16*g20*g21*g22*g23*g25*g28*g30
    *g34*g37*g42*g44*g46*g48*g53 ),
Comm( g6, g1 ) / ( g11 ),
Comm( g7, g1 ) / ( g12*g20*g21*g22*g23*g35*g36*g43*g45 ),
Comm( g8, g1 ) / ( g14 ),
Comm( g10, g1 ) / ( g11*g12*g14*g15*g16*g20*g21*g22*g23*g24*g25*g28
    *g30*g37*g38*g39*g40*g42*g43*g44*g45*g46*g48*g49*g50*g52*g53
     ),
Comm( g11, g1 ) / ( g19 ),
Comm( g12, g1 ) / ( g21 ),
Comm( g13, g1 ) / ( g23 ),
Comm( g14, g1 ) / ( g25 ),
Comm( g15, g1 ) / ( g20*g21*g26*g34*g36*g38*g43*g45*g46 ),
Comm( g16, g1 ) / ( g28 ),
Comm( g18, g1 ) / ( g19*g20*g25*g26*g28*g29*g30*g33*g34*g37*g38*g39
    *g43*g44*g45*g46*g47*g48*g51*g53 ),
Comm( g19, g1 ) / ( g33 ),
Comm( g20, g1 ) / ( g35 ),
Comm( g21, g1 ) / ( g34 ),
Comm( g22, g1 ) / ( g38 ),
Comm( g23, g1 ) / ( g36*g37*g38 ),
Comm( g24, g1 ) / ( g39 ),
Comm( g25, g1 ) / ( g42 ),
Comm( g26, g1 ) / ( g44 ),
Comm( g27, g1 ) / ( g46 ),
Comm( g28, g1 ) / ( g48 ),
Comm( g29, g1 ) / ( g34*g35*g43*g44*g49 ),
Comm( g30, g1 ) / ( g51 ),
Comm( g32, g1 ) / ( g33*g34*g42*g43*g48*g49*g51*g52*g53 ),
Comm( g3, g2 ) / ( g7 ),
Comm( g4, g2 ) / ( g8 ),
Comm( g6, g2 ) / ( g12 ),
Comm( g7, g2 ) / ( g13 ),
Comm( g8, g2 ) / ( g15 ),
Comm( g9, g2 ) / ( g16 ),
Comm( g11, g2 ) / ( g20 ),
Comm( g12, g2 ) / ( g22 ),
Comm( g13, g2 ) / ( g24 ),
Comm( g14, g2 ) / ( g26 ),
Comm( g15, g2 ) / ( g27 ),
Comm( g16, g2 ) / ( g29 ),
Comm( g17, g2 ) / ( g30 ),
Comm( g19, g2 ) / ( g34 ),
Comm( g20, g2 ) / ( g36 ),
Comm( g21, g2 ) / ( g37 ),
Comm( g22, g2 ) / ( g39 ),
Comm( g23, g2 ) / ( g40 ),
Comm( g24, g2 ) / ( g41 ),
Comm( g25, g2 ) / ( g43 ),
Comm( g26, g2 ) / ( g45 ),
Comm( g27, g2 ) / ( g47 ),
Comm( g28, g2 ) / ( g49 ),
Comm( g29, g2 ) / ( g50 ),
Comm( g30, g2 ) / ( g52 ),
Comm( g31, g2 ) / ( g53 ),
Comm( g4, g3 ) / ( g14*g20*g21*g35*g37*g38*g43 ),
Comm( g5, g3 ) / ( g12*g13*g15*g23*g27*g36*g38*g39*g40 ),
Comm( g6, g3 ) / ( g20*g21*g35*g37*g38*g43 ),
Comm( g7, g3 ) / ( g22*g23*g38*g45 ),
Comm( g8, g3 ) / ( g20*g21*g34*g37*g38*g44 ),
Comm( g9, g3 ) / ( g28*g34*g35*g43*g44 ),
Comm( g10, g3 ) / ( g20*g22*g26*g27*g29*g34*g36*g37*g38*g39*g40*g41
    *g46*g50 ),
Comm( g11, g3 ) / ( g34*g35 ),
Comm( g12, g3 ) / ( g37*g38 ),
Comm( g13, g3 ) / ( g39*g40 ),
Comm( g14, g3 ) / ( g43*g44 ),
Comm( g15, g3 ) / ( g36*g37*g45*g46 ),
Comm( g16, g3 ) / ( g34*g35*g43*g44 ),
Comm( g17, g3 ) / ( g51 ),
Comm( g18, g3 ) / ( g34*g36*g43*g45*g49*g50*g52 ),
Comm( g5, g4 ) / ( g15*g16*g25*g28*g29*g30*g42*g43*g48*g49*g52*g53
     ),
Comm( g6, g4 ) / ( g25*g34*g35*g48 ),
Comm( g7, g4 ) / ( g20*g21*g26*g34*g36*g38*g44*g45*g46*g49 ),
Comm( g8, g4 ) / ( g25*g28*g42*g48 ),
Comm( g10, g4 ) / ( g29*g30*g43*g47*g48*g49*g51*g53 ),
Comm( g11, g4 ) / ( g42 ),
Comm( g12, g4 ) / ( g34*g35*g44 ),
Comm( g13, g4 ) / ( g36*g37*g46 ),
Comm( g14, g4 ) / ( g42*g48 ),
Comm( g15, g4 ) / ( g43*g49 ),
Comm( g16, g4 ) / ( g48*g51 ),
Comm( g18, g4 ) / ( g52*g53 ),
Comm( g6, g5 ) / ( g20*g22*g26*g34*g36*g37*g43*g49 ),
Comm( g7, g5 ) / ( g22*g24*g27*g36*g39*g40*g45*g50 ),
Comm( g8, g5 ) / ( g27*g29*g43*g49 ),
Comm( g9, g5 ) / ( g29*g30 ),
Comm( g11, g5 ) / ( g34*g36*g43 ),
Comm( g12, g5 ) / ( g36*g39*g45 ),
Comm( g13, g5 ) / ( g39*g41*g47 ),
Comm( g14, g5 ) / ( g45*g49 ),
Comm( g15, g5 ) / ( g47*g50 ),
Comm( g16, g5 ) / ( g50*g52 ),
Comm( g17, g5 ) / ( g52*g53 ),
Comm( g7, g6 ) / ( g36*g38 ),
Comm( g8, g6 ) / ( g34*g35*g43*g44 ),
Comm( g9, g6 ) / ( g48 ),
Comm( g10, g6 ) / ( g34*g36*g43*g45*g49 ),
Comm( g8, g7 ) / ( g45*g46 ),
Comm( g9, g7 ) / ( g34*g35*g43*g44*g49 ),
Comm( g10, g7 ) / ( g36*g39*g45*g47*g50 ),
Comm( g9, g8 ) / ( g48*g51 ),
Comm( g10, g8 ) / ( g50*g52 ),
Comm( g10, g9 ) / ( g52*g53 )
] );;
G2 := AgGroupFpGroup( GROUP );;

# Now use 'IsConsistent' to check the collectors
InfoAgGroup1 := Print;;
InfoAgGroup2 := Print;;

# single collector
SetCollectorAgWord( G2.identity, "single" );
SetCollectorAgWord( G3.identity, "single" );
SetCollectorAgWord( G7.identity, "single" );
Unbind( G2.isConsistent );
Unbind( G7.isConsistent );
IsConsistent( G2, true );
#>true
IsConsistent( G3, true );
#>true
IsConsistent( G7, true );
#>true

# combinatorial collector
SetCollectorAgWord( G2.identity, "combinatorial" );
SetCollectorAgWord( G3.identity, "combinatorial" );
SetCollectorAgWord( G7.identity, "combinatorial" );
Unbind( G2.isConsistent );
Unbind( G3.isConsistent );
Unbind( G7.isConsistent );
IsConsistent( G2, true );
#>true
IsConsistent( G3, true );
#>true
IsConsistent( G7, true );
#>true

# triple collector
SetCollectorAgWord( G3.identity, "triple", 3 );
SetCollectorAgWord( G7.identity, "triple", 7 );
Unbind( G3.isConsistent );
Unbind( G7.isConsistent );
IsConsistent( G3, true );
#>true
IsConsistent( G7, true );
#>true

# quadruple collector
SetCollectorAgWord( G3.identity, "quadruple", 3 );
SetCollectorAgWord( G7.identity, "quadruple", 7 );
Unbind( G3.isConsistent );
Unbind( G7.isConsistent );
IsConsistent( G3, true );
#>true
IsConsistent( G7, true );
#>true


#############################################################################
##
##  That is it for collecting with agwords.
##
Print( "$Id: collecto.tst,v 1.1.1.1 1996/12/11 12:44:37 werner Exp $\n",
       QuoInt( 6423720000, time ), " GAPstones\n" );
