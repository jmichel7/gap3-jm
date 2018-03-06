/*****************************************************************************
**
**    pc.c                            NQ                       Werner Nickel
**                                         Werner.Nickel@math.rwth-aachen.de
*/


#include "nq.h"

int	 Class = 0;
int      NrPcGens = 0;
int      NrCenGens = 0;
int      *Dimension;
 
word     *Generators;
int	 *Weight;
 
gen      *Commute;
exp      *Exponent;
word     *Power;
word     **Conjugate;
def	 *Definition;
char	 **PcGenName;

/*
**    InitPcPres() initializes those parts of pc-presentation which are
**    not initialized in elimMap().
*/
void	InitPcPres() {

	long	i, j, t;

	if( Verbose ) t = RunTime();

	Class = 1;

	Generators = (word*)malloc( (2*NrCenGens+1)*sizeof(word) );
	if( Generators == NULL ) {
	    perror( "InitPcPres(), Generators" );
	    exit( 2 );
	}
	Generators += NrCenGens;
	Weight = (int *)malloc( (NrCenGens+1)*sizeof(int) );
	if( Weight == NULL ) {
	    perror( "InitPcPres, Weight" );
	    exit( 2 );
	}
	for( i = -NrCenGens; i <= NrCenGens; i++ ) {
	    if( i == 0 ) continue;
	    if( i > 0  ) Weight[i] = Class;
	    Generators[i] = (word)malloc( 2*sizeof(gpower) );
	    if( Generators[i] == NULL ) {
		perror( "InitPcPres(), Generators[]" );
		exit( 2 );
	    }
	    Generators[i][0].g = i;   Generators[i][0].e = 1;
	    Generators[i][1].g = EOW; Generators[i][1].e = 0;
	}
	PcGenName = (char **)Allocate( (NrCenGens+1)*sizeof(char *) );

	Dimension = (int*)malloc( (Class+1)*sizeof(int) );
	if( Dimension == NULL ) {
	    perror( "InitPcPres(), Dimension" );
	    exit( 2 );
	}
	Dimension[Class] = NrCenGens;

	Conjugate = (word**)malloc( (2*NrCenGens+1)*sizeof(word*) );
	if( Conjugate == NULL ) {
	    perror( "InitPcPres(), Conjugate" );
	    exit( 2 );
	}
	Conjugate += NrCenGens;
	for( j = 1; j <= NrCenGens; j++ ) {
	    /* the length of Conjugate[j] is 2*(j-1)+1 */
	    Conjugate[j] = (word*)malloc( (2*j-1)*sizeof(word) );
	    if( Conjugate[j] == NULL ) {
		perror( "InitPcPres(), Conjugate[]" );
		exit( 2 );
	    }
	    Conjugate[j] += j-1;
	    for( i = -(j-1); i <= j-1; i++ )
		Conjugate[j][i] = Generators[j];
	    if( Exponent[j] == 0 ) {
		Conjugate[-j] = (word*)malloc( (2*j-1)*sizeof(word) );
		if( Conjugate[-j] == NULL ) {
		    perror( "InitPcPres(), Conjugate[]" );
		    exit( 2 );
		}
		Conjugate[-j] += j-1;
		for( i = -(j-1); i <= j-1; i++ )
		    Conjugate[-j][i] = Generators[-j];
	    }
	}

	/* Here central generators change their status to pc-generators. */
	NrPcGens += NrCenGens;
	NrCenGens = 0;

	if( Verbose )
	  printf("#    Initialized pc-presentation (%d msec).\n",RunTime()-t);
}

void	ExtPcPres() {

	long	i, j, c, N, oldsize, newsize, t;
	word	*tmp, **ttmp;

	if( Verbose ) t = RunTime();

	Class++;

	Weight = (int *)realloc( Weight, (NrPcGens+NrCenGens+1)*sizeof(int) );
	if( Weight == NULL ) {
	    perror( "InitPcPres, Weight" );
	    exit( 2 );
	}

	tmp = (word*)malloc( (2*(NrPcGens+NrCenGens)+1)*sizeof(word) );
	if( tmp == NULL ) {
	    perror( "InitPcPres(), tmp" );
	    exit( 2 );
	}
	tmp += NrPcGens+NrCenGens;

	for( i = -(NrPcGens+NrCenGens); i <= (NrPcGens+NrCenGens); i++ ) {
	    if( i == 0 ) continue;
	    if( i > NrPcGens ) Weight[i] = Class;
	    if( i < -NrPcGens || i > NrPcGens  ) {
		tmp[i] = (word)malloc( 2*sizeof(gpower) );
		if( tmp[i] == NULL ) {
		    perror( "InitPcPres(), tmp[]" );
		    exit( 2 );
		}
		tmp[i][0].g = i;   tmp[i][0].e = 1;
		tmp[i][1].g = EOW; tmp[i][1].e = 0;
	    }
	    else
		tmp[i] = Generators[i];
	}
	free( Generators - NrPcGens );
	Generators = tmp;

	Dimension = (int*)realloc( Dimension, (Class+1)*sizeof(int) );
	if( Dimension == NULL ) {
	    perror( "InitPcPres(), Dimension" );
	    exit( 2 );
	}
	Dimension[Class] = NrCenGens;


	/* Now Conjugate[] has to be enlarged. */
	ttmp = (word**)malloc( (2*(NrPcGens+NrCenGens)+1)*sizeof(word*) );
	if( ttmp == NULL ) {
	    perror( "extPcPres(), tmp" );
	    exit( 2 );
	}
	ttmp += NrPcGens + NrCenGens;
	/* The contents of Conjugate[] must be copied to the new array. */
	for( i = -NrPcGens; i <= NrPcGens; i++ )
	    ttmp[i] = Conjugate[i];

	free( Conjugate-NrPcGens );
	Conjugate = ttmp;

	/*
	** The next nilpotency class to be calculated is Class+1. Therefore
	** commutators of weight Class+1, which are currently trivial, will
	** get new generators and tails. For the corresponding conjugates
	** space must be created in the array Conjugate[].
	**
	** Only those entries in Conjugate[] which do not have exceeded their
	** maximal length yet must be enlarged. This business is a little
	** bit tricky because the amount by which Conjugate[N], for a
	** generator N, has to be enlarged depends on the class of N.
	** The generators of highest class do not yet have any conjugates. They
	** will get a conjugate relation for each generator of weight 1,
	** therefore the size of the array for those generators is
	** 2*Dimension[1]+1. The array for generators of Class-1 has to be
	** enlarged by 2*Dimension[2] and so on.
	*/
	N = NrPcGens + NrCenGens;
	newsize = 0;
	for( c = Class; c > Class - c; c-- ) {
	    oldsize = newsize;
	    /* Compute the new size of the array for generators of class c.
	    ** Those generators get a new conjugate relation for each generator
	    ** of weight Class-c+1. */
	    newsize += Dimension[ Class - c + 1 ];
	    for( i = 1; i<= Dimension[c]; i++ ) {
		tmp = (word*)malloc( (2*min(N-1,newsize)+1)*sizeof(word) );
		if( tmp == NULL ) {
		    perror( "extPcPres(), tmp" );
		    exit( 2 );
		}
		tmp += min(N-1,newsize);
		if( c < Class ) {
		    /* Copy the contents to the new array. */
		    for( j = -oldsize; j <= oldsize; j++ )
			tmp[j] = Conjugate[N][j];
		    free( Conjugate[N] - oldsize );
		}
		Conjugate[N] = tmp;
		/* Initialise the new space. */
		for( j = oldsize+1; j <= min(N-1,newsize); j++ )
		    Conjugate[N][j] = Conjugate[N][-j] = Generators[N];

		if( Exponent[ N ] != 0 ) { N--; continue; }
		/* If the generator N is of infinite order, it also has
		** conjugate relations `on the other side'. All that has to
		** be done is exactly the same as before just for negative N.
		*/
		tmp = (word*)malloc( (2*min(N-1,newsize)+1)*sizeof(word) );
		if( tmp == NULL ) {
		    perror( "extPcPres(), tmp" );
		    exit( 2 );
		}
		tmp += min(N-1,newsize);
		if( c < Class ) {
		    for( j = -oldsize; j <= oldsize; j++ )
			tmp[j] = Conjugate[-N][j];
		    free( Conjugate[-N] - oldsize );
		}
		Conjugate[-N] = tmp;
		for( j = oldsize+1; j <= min(N-1,newsize); j++ )
		    Conjugate[-N][j] = Conjugate[-N][-j] = Generators[-N];
		N--;
	    }
	}

	/* Now the central generators have conjugate relations and so they
	** change theit status to pc-generators. */
	NrPcGens += NrCenGens;
	NrCenGens = 0;

	if( Verbose )
	    printf("#    Extended pc-presentation (%d msec).\n",RunTime()-t);
}

void	PrintPcPres() {

	gen	g;
	long	i, j, first = 1;

	if( Gap ) putchar( '#' );
	printf( "    <" );
	for( i = 1; i <= NrPcGens; i++ ) {
	    printGen( i, 'A' );
	    if( i < NrPcGens+NrCenGens ) putchar( ',' );
	}
	printf( "\n" );
	if( Gap ) putchar( '#' );
	printf( "     " );
	for( ; i <= NrPcGens+NrCenGens; i++ ) {
	    printGen( i, 'A' );
	    if( i < NrPcGens+NrCenGens ) putchar( ',' );
	}
	printf( " |" );
	for( i = 1; i <= NrPcGens+NrCenGens; i++ )
	    if( Exponent[i] != 0 ) {
		if( first ) { putchar( '\n' ); first = 0; }
		else	    printf( ",\n" );
		if( Gap ) putchar( '#' );
		printf( "        " );
		printGen( i, 'A' );
		printf( "^%d", Exponent[i] );
		if( Power[i] != (word)0 && Power[i]->g != EOW ) {
		    printf( " = " );
		    printWord( Power[i], 'A' );
		}
	    }

	for( j = 1; j <= NrPcGens; j++ ) {
	    i = 1;
	    while( i < j && Wt(i) + Wt(j) <= Class + (NrCenGens==0?0:1) ) {
		/* print Conjugate[j][i] */
		if( first ) { putchar( '\n' ); first = 0; }
		else	    printf( ",\n" );
		if( Gap ) putchar( '#' );
		printf( "        " );
		printGen( j, 'A' );
		putchar( '^' );
		printGen( i, 'A' );
		if( (g = Conjugate[j][i][1].g) != EOW
		     && Definition[g].h == j && Definition[g].g == i )
			printf( "           =: " );
		else	printf( "           =  " );
		printWord( Conjugate[j][i], 'A' );
		if( Exponent[i] == 0 ) {
		    if( first ) { putchar( '\n' ); first = 0; }
		    else          printf( ",\n" );
		    if( Gap ) putchar( '#' );
		    printf( "        " );
		    printGen( j, 'A' );
		    putchar( '^' );
		    putchar( '(' );
		    printGen( i, 'A' );
		    printf( "^-1)      =  " );
		    printWord( Conjugate[j][-i], 'A' );
		}
		if( 0 && Exponent[j] == 0 ) {
		    if( first ) { putchar( '\n' ); first = 0; }
		    else          printf( ",\n" );
		    if( Gap ) putchar( '#' );
		    printf( "        " );
		    putchar( '(' );
		    printGen( j, 'A' );
		    printf( "^-1)^" );
		    printGen( i, 'A' );
		    printf( "      =  " );
		    printWord( Conjugate[-j][i], 'A' );
		}
		if( 0 && Exponent[i] + Exponent[j] == 0 ) {
		    if( first ) { putchar( '\n' ); first = 0; }
		    else          printf( ",\n" );
		    if( Gap ) putchar( '#' );
		    printf( "        " );
		    putchar( '(' );
		    printGen( j, 'A' );
		    printf( "^-1)^" );
		    putchar( '(' );
		    printGen( i, 'A' );
		    printf( "^-1) =  " );
		    printWord( Conjugate[-j][-i], 'A' );
		}
		i++;
	    }
	}
	printf( " >\n" );

	printf( "\n#    Class : %d\n", Class );
	printf( "#    Nr of generators of each class :" );
	for( i = 1; i <= Class; i++ ) printf( " %d", Dimension[i] );
	printf( "\n" );
}

void	sizePcPres() {

	long	size = 0, nrPt = 0;
	gen	g, h;

	nrPt += 1;
	size += sizeof(gpower);                /* Identty */

	/* First calculate the size of all those arrays. */
	nrPt += 1;
	size += Class * sizeof(int);           /* Dimension[]. */
	nrPt += 1+2*NrPcGens;
	size += (2*NrPcGens+1) * sizeof(word); /* Generators[]. */
	size += 2*NrPcGens * 2*sizeof(gpower); /* Generators. */
	nrPt += 1;
	size += NrPcGens * sizeof(int);        /* Weight[]. */
	nrPt += 1;
	size += NrPcGens * sizeof(gen);        /* Commute[]. */
	nrPt += 1;
	size += NrPcGens * sizeof(exp);        /* Exponent[]. */
	nrPt += 1;
	size += NrPcGens * sizeof(Definition); /* Definition[]. */
	nrPt += 1;
	size += NrPcGens * sizeof(word);       /* Power[]. */
	for( g = 1; g <= NrPcGens; g++ )
	    if( Exponent[g] != 0 ) {
		nrPt += 1;
		size += sizeof(gpower)*(WordLength(Power[g])+1);
	    }

	nrPt += 1;
	size += (2*NrPcGens+1) * sizeof(word*);/* Conjugate[]. */
	for( h = 1; h <= NrPcGens; h++ ) {
	  nrPt += 1;
	  size += sizeof(word);
	  if( Exponent[h] == 0 ) {
	    nrPt += 1;
	    size += sizeof(word);
	  }
	  for( g = 1; g < h; g++ )
	    if( Wt(h)+Wt(g) <= Class+(NrCenGens==0?0:1) ) { 
	      nrPt += 1;
	      size += sizeof(word);
	      size += sizeof(gpower)*(WordLength(Conjugate[h][g])+1);
	      if( Exponent[h] == 0 ) {
		nrPt += 1;
	        size += sizeof(word);
	 	size += sizeof(gpower)*(WordLength(Conjugate[-h][ g])+1);
	      }
	      if( Exponent[g] == 0 ) {
		nrPt += 1;
	        size += sizeof(word);
		size += sizeof(gpower)*(WordLength(Conjugate[ h][-g])+1);
	      }
	      if( Exponent[h] + Exponent[g] == 0 ) {
		nrPt += 1;
	        size += sizeof(word);
		size += sizeof(gpower)*(WordLength(Conjugate[-h][-g])+1);
	      }
	    }
	}

	printf( "size of the presentation : %d bytes\n", size );
	printf( "pointers in the presentation : %d\n", nrPt );
}
