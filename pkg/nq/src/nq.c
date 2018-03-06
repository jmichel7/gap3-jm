/*****************************************************************************
**
**    nq.c                            NQ                       Werner Nickel
**                                         Werner.Nickel@math.rwth-aachen.de
*/


#include "nq.h"
#include "engel.h"

int	Debug = 0;
int	Gap = 0;
int     AbelianInv = 0;
int	Verbose = 0;

char	*InputFile;

static char	*ProgramName;
static int	Cl;

static
void	usage( error )
char	*error;

{	int	i;

	if( error != (char *)0 ) fprintf( stderr, "%s\n", error );
	fprintf( stderr, "usage: %s", ProgramName );
	fprintf( stderr, " [-a] [-d] [-g] [-v] [-s] [-f] [-c]\n" );
	for( i = strlen(ProgramName)+7; i > 0; i-- ) fputc( ' ', stderr );
        fprintf( stderr, " [-t <n>] [-l <n>] [-r <n>] [-n <n>] [-e <n>]\n" );
	for( i = strlen(ProgramName)+7; i > 0; i-- ) fputc( ' ', stderr );
	fprintf( stderr, " <presentation>  [<class>]\n" );
    	exit( 1 );
}

static int  leftEngel  = 0,
	    rightEngel = 0,
	 	 engel = 0,
           nrEngelGens = 1;

static
void	printHeader() {

	char	*s, hostname[128];

	gethostname( hostname, 128 );

	printf( "#\n" );
	printf(
	  "#    The ANU Nilpotent Quotient Program (Version %s)\n", VERSION );
	printf( "#    Calculating a nilpotent quotient\n" );
	printf( "#    Input: %s", InputFile );
	if( leftEngel ) {
            if( nrEngelGens > 1 )
                printf( " & the first %d generators are", nrEngelGens );
            else
		printf( " &" );
	    s = "th";
	    if( leftEngel == 1 ) s = "st";
	    if( leftEngel == 2 ) s = "nd";
	    if( leftEngel == 3 ) s = "rd";
	    printf( " %d%s left Engel", leftEngel, s );
	}
	if( rightEngel ) {
            if( nrEngelGens > 1 )
                printf( " & the first %d generators are", nrEngelGens );
            else
		printf( " &" );
	    s = "th";
	    if( rightEngel == 1 ) s = "st";
	    if( rightEngel == 2 ) s = "nd";
	    if( rightEngel == 3 ) s = "rd";
	    printf( " %d%s right Engel", rightEngel, s );
	}
	if( engel ) {
	    s = "th";
	    if( engel == 1 ) s = "st";
	    if( engel == 2 ) s = "nd";
	    if( engel == 3 ) s = "rd";
	    printf( " %d%s Engel", engel, s );
	}
	printf( "\n" );
	if( Cl != 666 ) printf( "#    Nilpotency class: %d\n", Cl );
	printf( "#    Program: %s", ProgramName );
	printf( "     Machine: %s\n#\n", hostname );
}

main( argc, argv )
int	argc;
char	*argv[];

{	FILE	*fp;
	long	t, time, start, begin;
	gen	g;
	extern	int	NrGens;
	extern	word	Epimorphism();

	CatchSignals();
	start = sbrk(0);
	begin = RunTime();

	ProgramName = argv[0]; argc--; argv++;

	setbuf( stdout, NULL );

	while( argc > 0 && argv[0][0] == '-' ) {
	    if( argv[0][2] != '\0' ) {
		fprintf( stderr, "unknown option: %s\n", argv[0] );
		usage( NULL );
	    }
	    switch( argv[0][1] ) {
		case 'r': if( --argc <= 1 ) usage("-r requires an argument");
			  argv++;
			  if( (rightEngel = atoi(argv[0])) <= 0 ) {
			      fprintf( stderr, "%s\n", argv[0] );
			      usage( "<n> must be positive." );
			  }
			break;
		case 'l': if( --argc <= 1 ) usage("-l requires an argument.");
			  argv++;
			  if( (leftEngel = atoi(argv[0])) <= 0 ) {
			      fprintf( stderr, "%s\n", argv[0] );
			      usage( "<n> must be positive." );
			  }
			break;
		case 'n': if( --argc <= 1 ) usage("-n requires an argument.");
			  argv++;
			  if( (nrEngelGens = atoi(argv[0])) <= 0 ) {
			      fprintf( stderr, "%s\n", argv[0] );
			      usage( "<n> must be positive." );
			  }
			break;
		case 'e': if( --argc <= 1 ) usage("-e requires an argument");
			  argv++;
			  if( (engel = atoi(argv[0])) <= 0 ) {
			      fprintf( stderr, "%s\n", argv[0] );
			      usage( "<n> must be positive." );
			  }
			break;
		case 't': if( --argc <= 1 ) usage("-t requires an argument");
			  argv++;
			  if( (t = atoi(argv[0])) <= 0 ) {
			      fprintf( stderr, "%s\n", argv[0] );
			      usage( "<n> must be positive." );
			  }
			  switch( argv[0][strlen(argv[0])-1] ) {
			    case 'd' : t *= 24;
			    case 'h' : t *= 60;
			    case 'm' : t *= 60;
			  }
			  SetTimeOut( t );
			break;
		case 'g': Gap = !Gap;
			break;
		case 'a': AbelianInv = !AbelianInv;
		        break;
		case 'v': Verbose = !Verbose;
			break;
		case 'd': Debug = !Debug;
			break;
		case 's': SemigroupOnly = !SemigroupOnly;
                        break;
                case 'c': CheckFewInstances = !CheckFewInstances;
                        break;
                case 'f': SemigroupFirst = !SemigroupFirst;
                        break;
                default : fprintf( stderr, "unknown option: %s\n", argv[0] );
			  usage( NULL );
			break;
	    }
	    argc--; argv++;
	}

	switch( argc ) {
	    case 1: InputFile = argv[0];
		     Cl = 666;
		 break;
	    case 2: InputFile = argv[0];
		     Cl = atoi( argv[1] );
		     if( Cl <= 0 ) usage( "<class> must be positive." );
		 break;
	    default: usage( NULL );
		 break;
	}
	
	if( (fp = fopen( InputFile, "r" )) == NULL ) {
	    perror( InputFile );
	    exit( 1 );
	}

	/* Read in the finite presentation. */
	Presentation( fp, InputFile );
	/* Set the number of generators. */
	WordInit( Epimorphism );
	InitEngel( leftEngel, rightEngel, engel, nrEngelGens );
	InitPrint( stdout );

	printHeader();

	time = RunTime();

	if( Gap ) printf( "NqLowerCentralSeries := [\n" );

	printf( "#    Calculating the abelian quotient ...\n" );
	InitEpim();

	NqEvalRelations();
	EvalEngel();

	ElimEpim();

	if( NrCenGens == 0 ) {
	    printf( "#    trivial abelian quotient\n" );
	    goto end;
	}

	if( Cl == 1 ) goto end;

	InitPcPres();

	printf( "#    The abelian quotient" );
	printf( " has %d generators\n", Dimension[Class] );
	printf( "#        with the following exponents:" );
	for( g = NrPcGens-Dimension[Class]+1; g <= NrPcGens; g++ )
	    printf( " %d", Exponent[g] );
	printf( "\n" );
	if(Verbose) {
	    printf("#    runtime       : %d msec\n",RunTime()-time);
	    printf("#    total runtime : %d msec\n",RunTime()-begin);
	    printf("#    total size    : %d byte\n",sbrk(0)-start);
	}
	printf( "#\n" );

	while( Class < Cl ) {
	    time = RunTime();
	    printf( "#    Calculating the class %d quotient ...\n", Class+1 );

	    AddGenerators();
	    Tails();
	    Consistency();
	    NqEvalRelations();
	    EvalEngel();
	    ElimGenerators();
	    if( NrCenGens == 0 ) goto end;
	    ExtPcPres();

	    printf( "#    Layer %d of the lower central series", Class );
	    printf( " has %d generators\n", Dimension[Class] );
	    printf( "#          with the following exponents:" );
	    for( g = NrPcGens-Dimension[Class]+1; g <= NrPcGens; g++ )
		printf( " %d", Exponent[g] );
	    printf( "\n" );
	    if(Verbose) {
		printf("#    runtime       : %d msec\n",RunTime()-time);
	        printf("#    total runtime : %d msec\n",RunTime()-begin);
	        printf("#    total size    : %d byte\n",sbrk(0)-start);
	    }
	    printf( "#\n" );

	}

end:
	TimeOutOff();
	printf( "\n\n#    The epimorphism :\n");
	PrintEpim();
	printf( "\n\n#    The nilpotent quotient :\n");
	PrintPcPres();
	printf( "\n" );
	printf( "#    total runtime : %d msec\n", RunTime() - begin );
	printf( "#    total size    : %d byte\n", sbrk(0) - start );

	if( Gap ) printf( "];\n" );
	TimeOutOn();

	return 0;
}
