/* file rws.h 
 * 18/12/95 change for V1.4.2 name of wtshortlex ordering changed to wtlex
 * throughout.
 * 1.4. 95. components rkminlen, rkmineqns added
 * 6/10/94.
 */

typedef enum {SHORTLEX, RECURSIVE, RT_RECURSIVE, WTLEX, WREATHPROD}
	orderings;
/* Default ordering is SHORTLEX, the other is RECURSIVE path ordering.
 * WTLEX is where the generators have weights, and the length i
 * determined by adding up the weights of the generators in the strings.
 * (so SHORTLEX is WTLEX with all weights 1).
 * WREATHPROD is as defined in Sims' book, where generators have levels.
 * Note RECURSIVE is a special case of this, where levels are 1,2,3,...
 * If more orderings are required, add more constants here and alter
 * the compare routine.
 */

typedef struct {
  char *lhs;		/* left hand side of equation stored as string */
  char *rhs;		/* right hand side of equation stored as string */
  boolean done;		/* true if equation has been handled on a earlier run */
  boolean changed;	/* true if changed during tidying up */
  boolean eliminated;	/* true if eliminated during tidying up */
} reduction_equation;

typedef struct {
  char name[256];		/* name of system for printing */
  orderings ordering;		/* ordering used on strings */
  int *weight;			/* for wtlex ordering */
  int *level;			/* for wreathprod ordering */
  boolean confluent;		/* true if system has been proved confluent */
  int num_gens;
  int num_eqns;
  int num_inveqns;		/* number of equations at beginning of list
				 * of length 2 and form g*inv_of[g] = 1 */
  int  *inv_of;			/* list of inverse numbers of generators */
  char **gen_name;		/* names of generators for printing */
  reduction_equation *eqns;	/* list of equations */
  fsa *reduction_machine;	/* finite state automaton for word reduction */
  int maxlenleft;		/* max stored length of lhs of equation */
  int maxlenright;		/* max stored length of rhs of equation */
  int maxoverlaplen;		/* max overlap length processed */
  int rkminlen;			/* min length for Rabin-Karp word reduction */
  int rkmineqns;		/* min no of eqns for Rabin-Karp reduction */
  int *history;			/* state history on reduction */
  int **slowhistory;
  int *slowhistorysp;		/* for use in slow-reduction */
  int *preflen;
  int *prefno;			/* for recording prefixes of left-hand-sides
                                 * when building reduction fsa. */
  char *testword1;
  char *testword2;		/* for storing words during reduction */
  boolean sorteqns;		/* true if equations are to be sorted 
				 * by length before output. */
  int maxoplen;			/* max length of equations output -
				 * only used if sorteqns is true */
} rewriting_system;
