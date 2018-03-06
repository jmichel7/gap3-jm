/*****************************************************************************
**
**    word.c                          NQ                       Werner Nickel
**                                         Werner.Nickel@math.rwth-aachen.de
*/


#include "nq.h"

void	printGen( g, c )
gen	g;
char	c;

{	putchar( c + (g-1) % 26 );
	if( (g-1) / 26 != 0 )
	    printf( "%d", (g-1) / 26 );
}

void	printWord( w, c )
word	w;
char	c;

{	if( w == (word)0 || w->g == EOW ) {
	    printf( "Id" );
	    return;
	}
	
	while( w->g != EOW ) {
	    if( w->g > 0 ) {
		printGen( w->g, c );
		if( w->e != 1 )
		    printf( "^%d", w->e );
	    }
	    else {
		printGen( -w->g, c );
		printf( "^%d", -w->e );
	    }
	    w++;
	    if( w->g != EOW ) putchar( '*' );
	}
}
