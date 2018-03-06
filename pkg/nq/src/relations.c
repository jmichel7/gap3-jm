/*****************************************************************************
**
**    relations.c                     NQ                       Werner Nickel
**                                         Werner.Nickel@math.rwth-aachen.de
*/


#include "nq.h"
#include "presentation.h"

static word	*Image;

void	NqEvalRelations() {

	long	t;
	expvec	ev;
	node	*r;
	word	w;

	if( Verbose ) t = RunTime();

	r = FirstRelation();
	while( !EarlyStop && r != (node *)0 ) {
	    w = (word)EvalNode( r );
	    ev = ExpVecWord( w );
	    addRow( ev );
	    r = NextRelation();
	}

	if( Verbose )
	    printf("#    Evaluated Relations (%d msec).\n",RunTime()-t);
}

/*
**    InitEpim() sets up the map from the generators of a finitely presented
**    group onto the generators of a free abelian group. It also sets up
**    the necessary data structures for collection.
*/
void	InitEpim() {

	long	i, t, nrGens;
	
	if( Verbose ) t = RunTime();

	/* Set the number of central generators to the number of generators
	** in the finite presentation. */
	nrGens = NumberOfGens();
	NrCenGens = nrGens;

	/* Initialize Exponent[]. */
	Exponent = (exp*) calloc( (NrCenGens+1), sizeof(exp) );
	if( Exponent == NULL ) {
	    perror( "initEpim(), Exponent" );
	    exit( 2 );
	}

	/* Initialize Commute[]. */
	Commute = (gen*)malloc( (NrCenGens+1) * sizeof(gen) );
	if( Commute == NULL ) {
	    perror( "initEpim(), Commute" );
	    exit( 2 );
	}
	for( i = 0; i <= NrCenGens; i++ ) Commute[i] = i;

	/* initialize the epimorphism onto the pc-presentation. */
	Image = (word*)malloc( (nrGens+1)*sizeof(word) );
	if( Image == NULL ) {
	    perror( "initEpim(), Image" );
	    exit( 2 );
	}
	for( i = 1; i <= nrGens; i++ ) {
	    Image[i] = (word)malloc( 2*sizeof(struct gpower) );
	    if( Image[i] == NULL ) {
		perror( "initEpim(), Image[]" );
		exit( 2 );
	    }
	    Image[i][0].g = i;
	    Image[i][0].e = 1;
	    Image[i][1].g = EOW;
	    Image[i][1].e = 0;
	}

	if( Verbose )
	    printf("#    Initialized epimorphism (%d msec).\n",RunTime()-t);
}

int	ExtendEpim() {

	int	j, l, G, nrGens;
	word	w;

	G = NrPcGens;
	nrGens = NumberOfGens();

	/* If there is an epimorphism, we have to add pseudo-generators
	** to the right hand side of images which are not definitions. */
	for( j = 1; j <= Dimension[1]; j++ )
	    Image[ -Definition[j].h ] =
		  (word)((unsigned long)(Image[-Definition[j].h]) | 0x1);

	for( j = 1; j <= nrGens; j++ )
	    if( !((unsigned long)(Image[j]) & 0x1) ) {
		G++;
		l = 0;
		if( Image[j] != (word)0 ) l = WordLength( Image[ j ] );
		w = (word)malloc( (l+2)*sizeof(gpower) );
		if( Image[j] != (word)0 ) WordCopy( Image[ j ], w );
		w[l].g   = G;   w[l].e   = 1;
		w[l+1].g = EOW; w[l+1].e = 0;
		if( Image[ j ] != (word)0 ) free( Image[ j ] );
		Image[ j ] = w;
	    }

	for( j = 1; j <= Dimension[1]; j++ )
	    Image[ -Definition[j].h ] =
		 (word)((unsigned long)(Image[-Definition[j].h]) & ~0x1);

	return G - NrPcGens;
}

int	ElimAllEpim( n, M, renumber  )
int	n;
expvec	*M;
int	*renumber;

{	int	i, j, l, nrGens;
	word	w;

	nrGens = NumberOfGens();

	/* first we eliminate ALL central generators that occur in the
	** epimorphism. */
	for( j = 1; j <= Dimension[1]; j++ )
	    Image[ -Definition[j].h ] =
		(word)((unsigned long)(Image[-Definition[j].h]) | 0x1);

	for( j = 1, i = 0; j <= nrGens; j++ )
	    if( !((unsigned long)(Image[j]) & 0x1) ) {
		l = WordLength( Image[j] );
		w = (word)Allocate( (l+NrCenGens+1-n)*sizeof(gpower) );
		WordCopy( Image[j], w );
		l--;
		l += appendExpVector( w[l].g+1-NrPcGens, M[i], w+l, renumber );

		if( Image[j] != (word)0 ) free( Image[j] );
		if( l == 1 ) {
		    Image[j] = (word)0;
		    free( w );
		}
		else
		    Image[j] = (word)realloc( w, l*sizeof(gpower) );
		i++;
	    }

	for( j = 1; j <= Dimension[1]; j++ )
	    Image[ -Definition[j].h ] =
		(word)((unsigned long)(Image[-Definition[j].h]) & ~0x1);

	return i;
}

void	ElimEpim() {

	long	i, j, h, l, n = 0, t;
	gen	*renumber;
	expvec	*M;
	word	w;
	extern  expvec	*MatrixToExpVecs();

	if( Verbose ) t = RunTime();

	M = MatrixToExpVecs();

	renumber = (gen*) malloc( (NrCenGens+1) * sizeof(gen) );
	if( renumber == NULL ) {
	    perror( "ElimEpim(), renumber" );
	    exit( 2 );
	}

	/* first assign a new number to each generator which is
	   not to be eliminated. */
	for( h = 1, i = 0; h <= NrCenGens; h++ )
	    if( i >= NrRows || h != Heads[i] )
	        renumber[ h ] = h - n;
	    else if( M[i][h] != 1 ) { /* h will become a torsion element */
		renumber[ h ] = h - n;
		Exponent[ renumber[h] ] = M[i][h];
		i++;
	    }
	    else {                    /* h will be eliminated. */
		n++;
		i++;
	    }

	/* allocate memory for Power[], note that n is the number of
	   generators to be eliminated. */
	Power = (word*) malloc( (NrCenGens-n+1) * sizeof(word) );
	if( Power == NULL ) {
	    perror( "ElimEpim(), Power" );
	    exit( 2 );
	}

	/* allocate memory for Definition[]. */
	Definition = (def*)malloc( (NrCenGens-n+1)*sizeof(def) );
	if( Definition == NULL ) {
	    perror( "ElimEpim(), Definition" );
	    exit( 2 );
	}

	/* Now eliminate and renumber generators. */
	for( h = 1, i = 0; h <= NrCenGens; h++ ) {
	    /* h runs through all generators. Only if a generator is
	    ** encountered that occurs as the i-th head we have to work. */
	    if( i >= NrRows || h != Heads[i] ) {
		/* generator i survives and does not get a power relation */
		Image[h][0].g = renumber[ h ];
		Definition[ renumber[h] ].h = -h;
		Definition[ renumber[h] ].g = 0;
		continue;
	    }

	    /* From here on we have that  h = Heads[i]. */
	    w = (word)malloc( (NrCenGens+1-h)*sizeof(gpower) );
	    if( w == NULL ) {
		perror( "ElimEpim(), w" );
		exit( 2 );
	    }

	    /* Copy the exponent vector M[i] into w. */
	    for( l = 0, j = h+1; j <= NrCols; j++ )
		if( M[i][j] > 0 ) {
		    w[l].g = -renumber[j];
		    w[l].e = M[i][j];
		    l++;
		}
		else if( M[i][j] < 0 ) {
		    w[l].g = renumber[j];
		    w[l].e = -M[i][j];
		    l++;
		}
	    w[l].g = EOW;
	    w[l].e = 0;
	    l++;

	    if( M[i][h] == 1 ) {
		/* generator h has to be eliminated. */
		free( Image[h] );
		Image[h] = (word)realloc( w, l*sizeof(gpower) );
	    }
	    else {
		/* generator h survives and gets a power relation. */
		Image[h][0].g = renumber[ h ];
		Definition[ renumber[h] ].h = -h;
		Definition[ renumber[h] ].g = 0;
		Power[ renumber[h]]   = (word)realloc( w, l*sizeof(gpower) );
	    }
	    i++;
	}

	/* Now adjust the sizes of the arrays */
	Commute = (gen*)realloc( Commute, (NrCenGens+1-n)*sizeof(gen) );
	Exponent = (exp*)realloc( Exponent, (NrCenGens+1-n)*sizeof(exp) );

	free( renumber );

	freeExpVecs( M );
	NrCenGens -= n;

	if( Verbose )
	    printf("#    Eliminated generators (%d msec).\n",RunTime()-t);
}

void	PrintEpim() {

	long	i, nrGens;

	if( Image == NULL ) {
	    printf( "#    No map set.\n" );
	    return;
	}

	nrGens = NumberOfGens();
	for( i = 1; i <= nrGens; i++ ) {
	    printf( "#    " ); printGen( i, 'a' ); printf( "|---> " );
	    printWord( Image[i], 'A' );
	    putchar( '\n' );
	}
}

word	Epimorphism( g )
gen	g;

{	return Image[g];	}
