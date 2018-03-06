#############################################################################
##
#A  pq.tst    		        GAP Test                       Alice Niemeyer
##
#A  @(#)$Id: pq.tst,v 1.1.1.1 1996/12/11 12:44:37 werner Exp $
##
#Y  Copyright 1991, Mathematics Research Section, SMS, ANU Canberra, Australia
##
##  This files contains test examples for the p-quotient algorithm. 
##
#H  @(#)$Log: pq.tst,v $
#H  @(#)Revision 1.1.1.1  1996/12/11 12:44:37  werner
#H  @(#)Preparing 3.4.4 for release
#H  @(#)
#H  Revision 3.2  1993/01/22  15:35:47  martin
#H  changed 'Word' to 'AbstractGenerator', and '.relations' to '.relators'
#H
#H  Revision 3.1  1991/09/13  14:31:57  fceller
#H  Initial Release under RCS
#H
##
InfoPQ1 := Ignore;;
InfoPQ2 := Ignore;;

a := AbstractGenerator("a");;
b := AbstractGenerator("b");;
c := AbstractGenerator("c");;
d := AbstractGenerator("d");;
e := AbstractGenerator("e");;
f := AbstractGenerator("f");;
g := AbstractGenerator("g");;

G := rec( generators := [ a, b ], relators := [ a*b ], exponents := [ 5 ] );;
pQuotient(G,5,5).dimensions;
#>[ 2, 2, 3, 4, 7 ]

G := rec( generators := [ a, b, c ], relators := [ a, b, c, c*a*b^-1 ],
   exponents := [ 4, 2, 2, 1 ] );;
pQuotient(G,2,0).dimensions;
#>[ 2, 1 ]

G := rec( generators := [ a, b ], relators := [  ] );;
pQuotient(G,13,5).dimensions;
#>[ 2, 3, 5, 8, 14 ]

G := rec( generators := [a,b,c,d,e,f],
relators := [a^5/(f ),b^5,c^5,d^5,e^5,f^5,
Comm( b, a )/(c *d *f^4),
Comm( c, a )/(d^2*e^2),Comm( c, b )/(e ),
Comm( d, a )/(e^3*f^4),Comm( d, b )/(f^2),
Comm( e, a )/(f^4) ],
exponents := [1,1,1,1,1,1,1,1,1,1,1,1] );;
pQuotient(G,5,0).dimensions;
#>[ 2, 1, 1, 1, 1 ]

G := rec( generators := [ a, b, c, d, e, f ],
   relators := [ b^-1*a^-1*b*a*c^-1, c^-1*a^-1*c*a,
      c^-1*b^-1*c*b, d^-1*a^-1*d*a, d^-1*b^-1*d*b*f^-1,
      e^-1*a^-1*e*a*f^-1, e^-1*b^-1*e*b, f^-1*a^-1*f*a,
           f^-1*b^-1*f*b, a^2*d^-1, b^2*e^-1, c^2*f^-1, d^2, e^2,
      f^2 ], exponents := [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ] );;
pQuotient( G, 2, 0 ).dimensions;
#>[ 2, 3, 1 ]

G := rec( generators := [ a, b, c, d, e, f ],
   relators := [ b^-1*a^-1*b*a*c^-1, c^-1*a^-1*c*a*f^-1,
     c^-1*b^-1*c*b, d^-1*a^-1*d*a, d^-1*b^-1*d*b,
     e^-1*a^-1*e*a*f^-1, e^-1*b^-1*e*b, f^-1*a^-1*f*a,
      f^-1*b^-1*f*b, a^2*d^-1, b^2*e^-1, c^2*f^-1, d^2, e^2,
      f^2 ], exponents := [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ] );;
pQuotient( G, 2, 0 ).dimensions;
#>[ 2, 3, 1 ]

G := rec( generators := [ a, b, c, d, e, f, g ],
   relators := [ b^-1*a^-1*b*a*d^-1, c^-1*a^-1*c*a, 
      c^-1*b^-1*c*b*g^-1*f^-1, d^-1*a^-1*d*a*f^-1, 
      d^-1*b^-1*d*b, d^-1*c^-1*d*c, e^-1*a^-1*e*a, 
      e^-1*b^-1*e*b, e^-1*c^-1*e*c, f^-1*a^-1*f*a, 
      f^-1*b^-1*f*b, f^-1*c^-1*f*c, g^-1*a^-1*g*a, 
      g^-1*b^-1*g*b, g^-1*c^-1*g*c, a^2*e^-1, b^2*d^-1, c^2, 
      d^2*f^-1, e^2*g^-1, f^2, g^2 ],
   exponents := [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
       1, 1 ] );;
pQuotient( G, 2, 0 ).dimensions;
#>[ 3, 2, 2 ]


#############################################################################
##
##  That is it for collecting with agwords.
##
Print( "$Id: pq.tst,v 1.1.1.1 1996/12/11 12:44:37 werner Exp $\n",
       QuoInt( 32900000, time+1 ), " GAPstones\n" );
