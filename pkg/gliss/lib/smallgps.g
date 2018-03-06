###############################################################################
##
##  smallgps.g  GLISSANDO ver 1.0  Christof Noebauer   1995, 1996
##                                                       
##
#############################################################################
##  This defines all groups of low order without using any library.
#############################################################################
##
C2            := Group( (1,2) );      
C2.name       := "C2";  # cyclic group order 2

C3            := Group( (1,2,3) );       
C3.name       := "C3";  # cyclic group order 3  

C4            := Group( (1,2,3,4) );        
C4.name       := "C4";  # cyclic group order 4 
V4            := DirectProduct( C2, C2 ); 
V4.name       := "V4";  # Klein's four group

C5            := Group( (1,2,3,4,5) );         
C5.name       := "C5";  # cyclic group order 5

C6            := Group( (1,2,3,4,5,6) );        
C6.name       := "C6";  # cyclic group order 6 
S3            := Group( (1,2), (1,2,3) );     
S3.name       := "S3";  # S3 

C7            := Group( (1,2,3,4,5,6,7) );        
C7.name       := "C7";  # cyclic group order 7   

C8            := Group( (1,2,3,4,5,6,7,8) );        
C8.name       := "C8";  # cyclic group order 8      
C2xC4         := DirectProduct( C2, C4 );
C2xC4.name    := "C2xC4";
C4xC2         := DirectProduct( C4, C2 );
C4xC2.name    := "C4xC2";
C2xC2xC2      := DirectProduct( C2, C2, C2 );
C2xC2xC2.name := "C2xC2xC2";
D8            := Group( (1,2,3,4), (2,4) );      
D8.name       := "D8";  # dihedral group of order 8 
Q8            := Group( (1,2,3,4)(5,6,7,8), (1,5,3,7)(2,8,4,6) ); 
Q8.name       := "Q8";  # quaternion group of order 8

C9            := Group( (1,2,3,4,5,6,7,8,9) );        
C9.name       := "C9";  # cyclic group of order 9
C3xC3         := DirectProduct( C3, C3 );
C3xC3.name    := "C3xC3";

C10           := Group( (1,2,3,4,5,6,7,8,9,10) );         
C10.name      := "C10"; # cyclic group of order 10 
D10           := Group( (1,2,3,4,5), (2,5)(3,4) );     
D10.name      := "D10"; # dihedral group of order 10

C11           := Group( (1,2,3,4,5,6,7,8,9,10,11) );       
C11.name      := "C11"; # cyclic group of order 11

C12           := Group( (1,2,3,4,5,6,7,8,9,10,11,12) );       
C12.name      := "C12"; # cyclic group of order 12 
C2xC6         := DirectProduct( C2, C6 );
C2xC6.name    := "C2xC6";
D12           := Group( (1,2,3,4,5,6), (2,6)(3,5) );
D12.name      := "D12";
A4            := Group( (1,2,4), (2,3,4) );
A4.name       := "A4";
# replace T (bad name) with ZS
ZS            := Group( ( 1, 2, 4, 7, 9,12)( 3, 6,10,11, 8, 5), 
                        ( 1, 3, 7,11)( 2, 5, 9,10)( 4, 8,12, 6) );
ZS.name       := "T";  # ( ZS-metacyclic )

C13           := Group( (1,2,3,4,5,6,7,8,9,10,11,12,13) );
C13.name      := "C13";

C14           := Group( (1,2,3,4,5,6,7,8,9,10,11,12,13,14) );
C14.name      := "C14";
D14           := Group( (1,2,3,4,5,6,7), (2,7)(3,6)(4,5) );
D14.name      := "D14";

C15           := Group( (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15) );
C15.name      := "C15";
