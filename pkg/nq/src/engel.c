/****************************************************************************
**
**    engel.c                         NQ                       Werner Nickel
**                                         Werner.Nickel@math.rwth-aachen.de
*/


#include "nq.h"
#include "engel.h"

static	int	LeftEngel = 0, RightEngel = 0, Engel = 0;
static  int     NrEngelGens = 0;
static	int	NrWords;
static  int     Needed;
static	word	A;
        int     SemigroupOnly  = 0;
        int     SemigroupFirst = 0;
        int     CheckFewInstances = 0;

static
void	Error( v, w, type )
word	v, w;
char	type;

{	printf( "Overflow in collector computing [ " );
	printWord( v, 'a' );
	if( type == 'e' ) printf( " , %d ", Engel );
	if( type == 'l' ) printf( " , %d ", LeftEngel );
	if( type == 'r' ) printf( " , %d ", RightEngel );
	printWord( w, 'a' );
	printf( " ]\n" );
}

static
void	evalEngelRel( v, w )
word	v, w;

{	word	vs=v, ws=w, v1, vv = v;
	long	n, needed;

/*	printf( "evalEngelRel() called with : " );
	printWord( v, 'A' ); printf( "    " );
	printWord( w, 'A' ); putchar( '\n' ); */

	NrWords++;
	/* Calculate [ v, w, .., w ] */
	if( (v = Commutator( v, w )) == (word)0 ) {
	    Error( vv, w, 'e' );
	    return;
	}
	n = Engel-1;
	while( n-- > 0 ) { 
	    if( (v1 = Commutator( v, w )) == (word)0 ) {
		Error( vv, w, 'e' );
		free( v );
		return;
	    }
	    free( v );
	    v = v1;
	}

	needed = addRow( ExpVecWord( v ) );
	if( needed ) {
   	    printf( "#    [ " );
	    printWord( vs, 'a' );
	    for( n = Engel-1; n >= 0; n-- ) {
		printf( ", " );
	    	printWord( ws, 'a' );
	    }
	    printf( " ]\n" );
	}
        if( CheckFewInstances ) Needed |= needed;
        else                    Needed = 1;

	free( v );
}

static
void	buildPairs( u, i, g, v, wt, which )
word	u, v;
long	i;
gen	g;
long	wt, which;

{	long	save_wt;


        /* First we check if the Engel condition is trivially
           satisfied for weight reasons. The commutator
           [u, n v] is 1 if w(u) + n*w(v) > Class+1. */
        if( which == 1 && i == 1 &&
           Wt(abs(u[0].g)) + Engel*Wt(abs(v[0].g)) > Class+1 )
            return;

	if( wt == 0 && which == 1 && i > 0 ) {
	    evalEngelRel( u, v );
	    return;
	}

	/* Keep u and start to build v. */
	if( i > 0 && which == 2 ) buildPairs( v, 0, 1, u, wt, 1 );

	if( g > NrPcGens ) return;

	save_wt = wt;
	while( !EarlyStop && g <= NrPcGens && Wt(g) <= wt ) {
	    u[i].g   = g;
	    u[i].e   = 0;
 	    u[i+1].g = EOW;
	    while( !EarlyStop && Wt(g) <= wt ) {
		u[i].e++;
		if( Exponent[g] > 0 && Exponent[g] == u[i].e ) break;
		wt -= Wt(g);
		buildPairs( u, i+1, g+1, v, wt, which );
                /* now build the same word with negative exponent */
                if( !EarlyStop && !SemigroupOnly && Exponent[g] == 0 ) {
                    u[i].g *= -1;
                    buildPairs( u, i+1, g+1, v, wt, which );
                    u[i].g *= -1;
                }
	    }
	    wt = save_wt;
	    g++;
	}
	u[i].g = EOW;
	u[i].e = 0;
	if( EarlyStop || SemigroupOnly || !SemigroupFirst ) return;

        while( !EarlyStop && g <= NrPcGens && Wt(g) <= wt ) {
	    u[i].g   = -g;
	    u[i].e   = 0;
 	    u[i+1].g = EOW;
	    while( !EarlyStop && Wt(g) <= wt ) {
		u[i].e++;
		if( Exponent[g] > 0 && Exponent[g] == u[i].e ) break;
		wt -= Wt(g);
		buildPairs( u, i+1, g+1, v, wt, which );
                if( EarlyStop ) return;
                /* now build the same word with negative exponent */
                if( !EarlyStop && !SemigroupOnly && Exponent[g] == 0 ) {
                    u[i].g *= -1;
                    buildPairs( u, i+1, g+1, v, wt, which );
                    u[i].g *= -1;
                }
	    }
	    wt = save_wt;
	    g++;
	}
	u[i].g = EOW;
	u[i].e = 0;
}

static
void	evalEngel()

{	word	u, v;
	long	c;

	u = (word)Allocate( (NrPcGens+NrCenGens+1) * sizeof(gpower) );
	v = (word)Allocate( (NrPcGens+NrCenGens+1) * sizeof(gpower) );

        /* For `production purposes' I don't want tot run through      */
        /* those classes that don't yield non-trivial instances of the */
        /* Engel law. Therefore, we stop as soon as we ran through a   */
        /* class that didn't yield any non-trivial instances. This is  */
        /* done through the static variable Needed which is set by     */
        /* evalEngelRel() as soon as a non-trivial instance has been   */
        /* found if the flag CheckFewInstances (option -c) is set.     */
        Needed = 1;
	for( c = 2; !EarlyStop && Needed && c <= Class+1; c++ ) {
            Needed = 0;
	    u[0].g = EOW; u[0].e = 0;
	    v[0].g = EOW; v[0].e = 0;
	    NrWords = 0;
	    if(Verbose)
                printf("#    Checking pairs of words of weight %d\n",c);
	    buildPairs( u, 0, 1, v, c, 2 );
	    if(Verbose) printf( "#    Checked %d words.\n", NrWords );
	}

	for( ; !EarlyStop && c <= Class+1; c++ )
	    printf("#    NOT checking pairs of words of weight %d\n",c);
            
	free( u ); free( v );
}

static
void	evalRightEngelRel( w )
word	w;

{	word	ws = w, v, v1;
	long	n, needed;

/*	printf( "evalRightEngelRel() called with : " );*/
/*	printWord( w, 'A' );*/
/*	putchar( '\n' );*/

	NrWords++;
	/* Calculate [ a, w, .., w ] */
	if( (v = Commutator( A, w )) == (word)0 ) {
	    Error( A, w, 'r' );
	    return;
	}
	n = RightEngel-1;
	while( n-- > 0 ) {
	    if( (v1 = Commutator( v, w )) == (word)0 ) {
		Error( A, w, 'r' );
		free( v );
		return;
	    }
	    free( v );
	    v = v1;
	}

	needed = addRow( ExpVecWord( v ) );
	if( needed ) {
   	    printf( "#    [ " );
	    printWord( A, 'a' );
	    for( n = RightEngel-1; n >= 0; n-- ) {
		printf( ", " );
	    	printWord( ws, 'a' );
	    }
	    printf( " ]\n" );
	}

	free( v );
}

static
void	evalLeftEngelRel( w )
word	w;

{	word	ws = w, v, v1;
	long	n, needed;

/*	printf( "evalLeftEngelRel() called with : " );*/
/*	printWord( w, 'A' );*/
/*	putchar( '\n' );*/

	NrWords++;
	/* Calculate [ w, a, .., a ] */
	if( (v = Commutator( w, A )) == (word)0 ) {
	    Error( w, A, 'l' );
	    return;
	}
	n = LeftEngel-1;
	while( n-- > 0 ) {
	    if( (v1 = Commutator( v, A )) == (word)0 ) {
		Error( w, A, 'l' );
		free( v );
		return;
	    }
	    free( v );
	    v = v1;
	}

	needed = addRow( ExpVecWord( v ) );
	if( needed ) {
   	    printf( "#    [ " );
	    printWord( ws, 'a' );
	    for( n = LeftEngel-1; n >= 0; n-- ) {
		printf( ", " );
	    	printWord( A, 'a' );
	    }
	    printf( " ]\n" );
	}

	free( v );
}

static
void	buildWord( u, i, g, wt )
word	u;
long	i, wt;
gen	g;

{	long	save_wt;

	if( wt == 0 && i > 0 ) { 
	    if( RightEngel ) evalRightEngelRel( u );
	    if( LeftEngel )  evalLeftEngelRel( u );
	    return;
	}

	if( g > NrPcGens ) return;

	save_wt = wt;
	while( !EarlyStop && g <= NrPcGens && Wt(g) <= wt ) {
	    u[i].g   = g;
	    u[i].e   = 0;
	    u[i+1].g = EOW;
	    while( !EarlyStop && Wt(g) <= wt ) {
		u[i].e++;
		if( Exponent[g] > 0 && Exponent[g] == u[i].e ) break;
		wt -= Wt(g);
		buildWord( u, i+1, g+1, wt );
                /* now build the same word with negative exponent */
                if( !EarlyStop && !SemigroupOnly &&
                    !SemigroupFirst && Exponent[g] == 0 ) {
                    u[i].g *= -1;
                    buildWord( u, i+1, g+1, wt );
                    u[i].g *= -1;
                }
	    }
	    wt = save_wt;
	    g++;
	}
	u[i].g = EOW;
	u[i].e = 0;
        if( EarlyStop || SemigroupOnly || !SemigroupFirst ) return;
	while( !EarlyStop && g <= NrPcGens && Wt(g) <= wt ) {
	    u[i].g   = -g;
	    u[i].e   = 0;
	    u[i+1].g = EOW;
	    while( !EarlyStop && Wt(g) <= wt ) {
		u[i].e++;
		if( Exponent[g] > 0 && Exponent[g] == u[i].e ) break;
		wt -= Wt(g);
		buildWord( u, i+1, g+1, wt );
	    }
	    wt = save_wt;
	    g++;
	}
	u[i].g = EOW;
	u[i].e = 0;
}

static
void	evalLREngel() {

	word	u;
	int	n;
	long	cl;

	A = (word)Allocate( 2*sizeof(gpower) );
        u = (word)Allocate( (NrPcGens+NrCenGens+1) * sizeof(gpower) );

        for( n = 1; !EarlyStop && n <= NrEngelGens; n++ ) {
	    A[0].g = n;   A[0].e = 1;
	    A[1].g = EOW; A[1].e = 0;

	    for( cl = 2; !EarlyStop && cl <= Class+1; cl++ ) {
	        u[0].g = EOW; u[0].e = 0;
	        NrWords = 0;
	        if(Verbose) printf("#    Checking words of weight %d\n",cl-1);
	        buildWord( u, 0, 1, cl-1 );
	        if(Verbose) printf( "#    Checked %d words.\n", NrWords );
	    }
        }
	free( u ); free( A );
}

void	EvalEngel() {	

	long	t;

	if( Verbose ) t = RunTime();

	if( LeftEngel || RightEngel ) evalLREngel();
	if( Engel ) evalEngel();

	if( Verbose )
	    printf("#    Evaluated Engel condition (%d msec).\n",RunTime()-t);
}

void	InitEngel( l, r, e, n )
int	l, r, e;

{	LeftEngel = l;
	RightEngel = r;
	Engel = e;
        NrEngelGens = n;
}
