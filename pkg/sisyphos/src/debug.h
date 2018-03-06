/* 	$Id: debug.h,v 1.2 1995/01/05 17:29:30 pluto Exp $	 */
/* 	$Log: debug.h,v $
 * Revision 1.2  1995/01/05  17:29:30  pluto
 * Initial version under RCS control.
 *	 */

#define _IS(t,x)	(((t)1 << (x)) != 0)
#define NBITS(t)	(4 * (1 + _IS(t, 4) + _IS(t, 8) + _IS(t,12) + _IS(t,16) \
					   + _IS(t,20) + _IS(t,24) + _IS(t,28) + _IS(t,32) ) )
					   
