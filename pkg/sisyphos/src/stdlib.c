/* supplement for standard library */

/* 	$Id: stdlib.c,v 1.3 1995/03/07 12:09:12 pluto Exp $	 */
/* 	$Log: stdlib.c,v $
 * Revision 1.3  1995/03/07  12:09:12  pluto
 * Changed name of 'to_(upper|lower)' routines to avoid conflicts with
 * GNU readline.
 *
 * Revision 1.2  1995/01/05  17:12:20  pluto
 * Initial version under RCS control.
 *	 */

#ifndef lint
static char vcid[] = "$Id: stdlib.c,v 1.3 1995/03/07 12:09:12 pluto Exp $";
#endif /* lint */

#include <ctype.h>

char *stpblk (register char *p)
{
	while ( (*p < 33) && (*p != 0) ) ++p;
	return ( p );
}

char my_to_upper (register char c)
{
	if ( islower ( c ) )
		return ( toupper ( c ) );
	else
		return ( c );
}

char my_to_lower (register char c)
{
	if ( isupper ( c ) )
		return ( tolower ( c ) );
	else
		return ( c );
}
