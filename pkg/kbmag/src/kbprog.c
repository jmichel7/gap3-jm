/* kbprog.c

29/5/95 - introduced -mo option to limit length of overlap length in
          processed pairs.

27/3/95 - introduced Rabin-Karp as alternative reduction procedure for
	longer left hand sides, + associated "-rk minlen mineqns" option.

14/3/95 - introduced -ro option to resume, but re-including the original
          equations

3/3/95 - started introduction of wtlex and wreathprod orderings.

31/1/95 - bug corrected due to building full table when one lhs is a prefix
of another. Patch by doing repeated tidyings before building full table.
Also introduced variable lostinfo - this is set true if an equation is dicarded
because it is too long during tidying - when this haapens, the monoid may have
changed.
Warning message printed at end when this occurs.

18/1/95 various reorganisation - extra components to rws structure
 (including test1, test2 -> rws.testword1,rws.testword2).
 automatic doubling of maxslowhistoryspace when required
 option -mrl n  to set maxreducelen
 file rwsio.c reorganised and split into rwsio.c and rwsio2.c.

13/1/95 introduced exit code 0 for "successful" completeion.

24/12/94 introduced variable current_maxstates -
unlike maxstates, this may be increased (doubled) if necessary.
If maxstates is set explicitly, then current_maxstates and maxstates are
equal. Otherwise, the initial value of current_maxstates is determined
by allotting INIT_FSASPACE space to the reduction automaton.

21/12/94 introduced wd_record array to remember progress of calculation and
try to decide when to stop; also options -mt and -hf (see below).

9/11/94 introduced options -v, -silent for verbosity of output - same in
all other programs.

6/10/94 - reduction algorithms moved to file reduce.c
Also rws structure introduced.

27/9/94 -
(i) Changed field-name generators to generatorOrder in input/output files.
(ii) Introduced -op option to output prefixes corresponding to states of fsa.

22/9/94 -
(i) introduced boolean array done.
For i an equation number, done[i] true means that equation number i has
been processed. So if done[i] and done[j] are both true, then there is no
need to test for overlaps. This is mainly useful when re-running. 
(ii) The option -r (resume) to take input from group.kbprog, presumed
to be output of a previous run.
(iii) -cn confnum option introduced. If more than this number of pairs
of equations.eqns are tested for overlaps with no new equations.eqns being
discovered, then a fast check for confluence is made using routine conf_check.
confnum=0 switches this test off.

This is the version with GAP input/output by dfh starting 13/9/94
(using some ideas and code by Peter Proehle of 9/5/94)
A new feature - not all generators need have inverses listed, but all
listed inverses must be 2-sided, and the inv_of function must be involutory.
Only generators with inverses will be cancelled in the routine "insert".

!! options added by dfh starting 30/3/94. 
   -mlr maxlenleft maxlenright
   Only equations with lhs having length at most maxlenleft,
   rhs having length at most maxlenright are stored.

   -sort maxoplen
   Output rules up to length of lhs up to maxoplen in increasing order of length
   of lhs. (maxoplen=0 means no maximum)

   Some changes were made to this code by Jamie P. Curmi during
   February/March 1992.

NEW SYNOPSIS (for GAP version)
	  kbprog  [-r] [-ro] [-t tidyint] [-me maxeqns] [-ms maxstates]
		[-sort maxoplen] [-mlr maxlenleft maxlenright]
		[-mo maxoverlaplen] [-mrl maxreducelen]
		[-lex] [-rec] [-rtrec] [-wtlex] [-wreath]
		[-cn confnum] [-op] [-mt min_time] [-hf halting_factor]
		[-rk minlen mineqns]
		[-wd] [-mwd maxwdiffs] [-v] [-silent] groupname

OPTIONS
     -r		Resume. Take input from groupname.kbprog - output also
		goes to groupname.kbprog and overwrites input.

     -ro	Resume. Take input from groupname.kbprog, and also re-read
		the original equations from groupname and append them
		to the end - output also goes to groupname.kbprog and
		overwrites input.

     -t	    tidyint
	  	Set the value of tidyint equal to tidyint.

     -me    maxeqns
		At most maxeqns equations are allowed. If more are
		generated, tidying up and output occur.

     -ms    maxstates
		At most maxstates states are allowed in the fsa.
		If more are generated, tidying up and output occur.

     -sort  maxoplen
		The equations are sorted into increasing length order
		before output. Only those with lhs up to length maxoplen
		are output. maxoplen=0 means no restriction.

     -mlr   maxlenleft maxlenright
		Only equations with lhs and rhs of lengths at most
		maxlenleft and maxlenright are stored.

     -mo    maxoverlaplen
		Only overlaps of total length at most maxoverlaplen are
		processed.

     -lex	Short-lexicographical ordering
		** THE DEFAULT SETTING **
	  	This is the original ordering used in kbprog.c.
		The ordering is lexicographic, with longer words having
		higher order than shorter words.

     -rec	Recursive Path ordering
		ordering based on that described in the book "Confluent String
		Rewriting" by Matthias Jantzen, Defn 1.2.14, page 24.
		The algorithm is described at the header of 'rec_compare'

     -rtrec	Recursive Path ordering
		ordering based on that described in the book "Confluent String
		Rewriting" by Matthias Jantzen, Defn 1.2.14, page 24.
		The algorithm is described at the header of 'rt_rec_compare'
		rtrec orders from the right of the words, whereas
		rec orders from the left - sometimes one is more
		efficient than the other in particular examples.

    -wtlex	Weighted short-lexicographical ordering
		A generalisation of shortlex. All genertors have a weight
		(specified in the input file), and the weightsof the generators
		in a word are added up to determine the "length" of the word.
		So shortlex is the special case where all weights are 1.

    -wreath	Wreath product ordering, as defined in Sims' book.
		A generalisation of recursive, where the generators have
		levels, and words in genertators at a fixed level are
		compared using shortlex. Recursive is the special case
		where the levels are 1,2,3,...

    -cn	    confnum
		If confnum pairs of equations are processed for overlaps
		and no new equations are discovered, then a fast
		confluence test is performed. Setting confnum=0 turns this
		off.
    -mt     min_time
    -hf	    halting_factor
		Only applies when word-differences are being computed.
		(Roughly) if halting_factor is set to be > 0, and
		the number of equations and the number of states
		has increased by more than halting_factor percent since the
		number of word-differences was last less than what it now is,
		and at least min_time seconds have elapsed, then halt.
		Alternatively, halt irrespective of elapsed time if number of
		equations has increased by 2*halting_factor percent.
		These conditions are highly experimental.
    -op		The states of the fsa correspond to prefixes of the left
		hand sides of equations When this option is set, these
		prefixes are output as state-names.
    -rk	   minlen mineqns
		When there are at least mineqns equations, start using the
		Rabin-Karp reduction method on left-hand-sides having length
		at least minlen
    -wd		After each tidying, calculate word-differences and re-order
		equations so that those which result in the word-difference
		fsa changing get higher priority in the Knuth-Bendix.
		On completeion, the word-difference machines are output
		instead of the updated rewriting system.

    -mwd    maxwdiffs
		Set the maximum number of worddiffs to be maxwdiffs.

    -silent	No output to stdout.

    -v		Verbose output.


EXIT STATUS:
  0 if completion is successful, in the sense that either a confluent set
  is produced, or the halting-factor conditions are satsified when word-
  differences are being computed. This is for use in the "automata" shellscript.
  1 for badusage, failure of malloc, etc.,
  2 if there is output, but completion is not successful in the sense above.

CHANGES TO KBPROG
      The orderings options -red, -rtred above.  I have tried to make the code
      as easy as possible to modify in case new ordering options are to be
      included (such as weighted lexicographical ordering) in the future.

      The code has been 'beautified' using a C-Beautifier on our system.
      This helped with my understanding of the algorithms, and
      modifications.		- J.P.Curmi (1992)

*/


#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <sys/times.h>

#include "defs.h"
#include "fsa.h"
#include "rws.h"
#include "definitions.h"

#ifdef SYSTEMV
#define HZTIME 		100
#else
#define HZTIME		60
#endif
#define PROGNAME 	"kbprog"
#define TIDYINT		100
#define MAXEQNS		32767
#define MAXREDUCELEN	32767
#define INIT_FSASPACE	262144
#define MAXWDIFFS	512
#define MAXSLOWHISTORYSPACE 	131072
#define CONFNUM		500
#define MAXCYCLES	16384

#define	USAGE\
	"Usage: kbprog [-r] [-ro] [-t tidyint] [-me maxeqns] [-ms maxstates]\n\
	[-sort maxoplen] [-mlr maxlenleft maxlenright] [-mo maxoverlaplen]\n\
	[-mrl maxreducelen] [-v] [-silent]\n\
	[-lex] [-rec] [-rtrec] [-wtlex] [-wreath]\n\
	[-wd] [-mwd maxwdiffs] [-mt min_time] [-hf halting_factor]\n\
	[-rk minlen mineqns] [-cn confnum] [-op] [groupname]\n"


int 	        tidyint = TIDYINT,
                maxeqns = MAXEQNS,
                maxstates = 0,
		init_fsaspace = INIT_FSASPACE,
		current_maxstates,
		confnum = CONFNUM,
		oldconfnum,
		maxslowhistoryspace = MAXSLOWHISTORYSPACE,
                maxreducelen = MAXREDUCELEN,
                maxwdiffs = MAXWDIFFS,
		num_states,
                **table,
		exit_status = 0;
boolean         tidyon,
		tidyon_now,
		worddiffs = FALSE,
		onintr = FALSE,
                halting = FALSE,
                do_conf_test = FALSE,
                lostinfo = FALSE,
                tidyintset = FALSE,
		maxeqnsset = FALSE,
		hadmaxeqns = FALSE,
		maxstatesset = FALSE,
		orderingset = FALSE,
		silentset = FALSE,
		verboseset = FALSE,
		confnumset = FALSE,
                maxreducelenset = FALSE,
                maxwdiffsset = FALSE,
		resume = FALSE,
		resume_with_orig = FALSE,
		outputprefixes = FALSE,
		wd_reorder = TRUE,
		double_states = FALSE,
		rk_on = FALSE,
		hadlongoverlap = FALSE,
		*new_wd;
static int	maxpreflen,
                nneweqns,
		tot_eqns, /* like rws.num_eqns, but not decreased if we
                          * throw equations away when maxeqns exceeded. */
		hadct,
		maxhad = 0,
                lasttid;
static char	inf[100],   /* name of input file */
          	outf[100], /* name of output file for rewriting system 
                              = inf.kbprog   */
		outfr[100], /* name of output file for reduction_fsa
			     for rewriting system - inf.reduce */
          	outf1[100], /* name of output file for first wd-machine 
                              = inf.diff1   */
          	outf2[100], /* name of output file for second wd-machine 
                              = inf.diff2   */
		outflog[100]; /* name of logfile for automata run */
static	FILE	*rfile,
		*wfile;
struct tms  	time_buffer;
fsa		reduction_fsa;		
fsa		wd_fsa; /* The word-difference machine if required. */		
rewriting_system  rws;

void (*reduce_word)(); /* The word-reduction function - it can be any of
			* rws_reduce, slow_rws_reduce, slow_rws_reduce_rk.
                        */
boolean (*check_reduce_word)();/* The functionthat checks if a word is reduced.
			     * slow_check_rws_reduce or slow_check_rws_reduce_rk
			     */
/* The following items are used to record progress of calculation,
 * and for trying to decide when to stop - this is a difficult problem!
 * They are only used when we are calculating word-differences.
 */
struct wdr {
  int num_eqns;
  int num_states;
  int num_diff;
} *wd_record; /* one entry for each calculation of word-differences, after
	       * tidying */
int num_cycles = 0; /* number of times we have been through the tidying loop */
int eqn_factor;   /* Percentage by which number of equations has increased since
	           * number of word-differences increased. No sense in
		   * halting unless > 0.  */
int states_factor; /* Percentage by which number of states has increased since
	            * number of word-differences increased.
		    */
int min_time = 0;
int halting_factor = 0;
/* If halting_factor is set nonzero (by command-line option -hf) then we halt if
 * the program has been running for at least min_time seconds, and if
 * both eqn_factor and states_factor exceed halting_factor.
 * OR if eqn_factor exceeds 2*halting_factor. 
 */


/* Functions defined in this file: */
void interrupt();
int insert();
void consider();
int tidyup();
void build_quicktable();
int modify_table();
void build_fulltable();
int lex_compare();
int rec_compare();
int rt_rec_compare();
int compare();
int conf_check();
void sort_eqns();
int wd_sort_eqns();
void should_we_halt();
void output_and_exit();
void badusage();

/* Functions used in this file defined in other files: */
void read_kbinput();
void read_extra_kbinput();
void print_kboutput();
void print_wdoutput();
void rws_reduce();
void slow_rws_reduce();
boolean slow_check_rws_reduce();
void initialise_wd_fsa();
void build_wd_fsa();
void make_full_wd_fsa();
void fsa_table_init();
void fsa_table_dptr_init();
void clear_wd_fsa();
void rws_clear();
void fsa_clear();
void rk_init();
void rk_reset();
void rk_clear();
void add_lhs();
void slow_rws_reduce_rk();
boolean slow_check_rws_reduce_rk();

void
interrupt()
/* When the program receives an interrupt  signal, it continues until the
 * next tidying, and then stops and outputs.
 */
{
	onintr = TRUE;
	signal(SIGINT, SIG_DFL);
	signal(SIGKILL, SIG_DFL);
	signal(SIGQUIT, SIG_DFL);
}

main(argc, argv)
	int             argc;
	char           *argv[];
{
	int             i,
	                j,
	                k,
	                l,
	                onum_eqns,
	                onum_eqnsi,
	                cct,
	                arg;

	setbuf(stdout, (char *) 0);
	setbuf(stderr, (char *) 0);

/* First set some default values of the rewriting-system fields. They may be
 * overridden by command-line or input in file.
 */
	rws.maxlenleft = 0;
	rws.maxlenright = 0;
	rws.rkminlen = 0;
	rws.rkmineqns = 0;
	rws.weight = 0;
	rws.level = 0;
	rws.sorteqns = FALSE;
	rws.maxoplen = 0;
	rws.maxoverlaplen = 0;
	rws.ordering = SHORTLEX;

	arg = 1;
	while (argc > arg) {
		if (strcmp(argv[arg], "-r") == 0) {
			resume = TRUE;
		}
		else if (strcmp(argv[arg], "-ro") == 0) {
			resume = resume_with_orig = TRUE;
		}
		else if (strcmp(argv[arg], "-t") == 0) {
			tidyintset = TRUE;
			arg++;
			if (arg >= argc)
				badusage();
			tidyint = atoi(argv[arg]);
                }
		else if (strcmp(argv[arg], "-me") == 0) {
			maxeqnsset = TRUE;
			arg++;
			if (arg >= argc)
				badusage();
			maxeqns = atoi(argv[arg]);
		}
		else if (strcmp(argv[arg], "-ms") == 0) {
			maxstatesset = TRUE;
			arg++;
			if (arg >= argc)
				badusage();
			maxstates = atoi(argv[arg]);
		}
		else if (strcmp(argv[arg], "-sort") == 0) {
			rws.sorteqns = TRUE;
			arg++;
			if (arg >= argc)
				badusage();
			rws.maxoplen = atoi(argv[arg]);
		}
		else if (strcmp(argv[arg], "-mlr") == 0) {
			arg++;
			if (arg >= argc)
				badusage();
			rws.maxlenleft = atoi(argv[arg]);
			arg++;
			if (arg >= argc)
				badusage();
			rws.maxlenright = atoi(argv[arg]);
			if (rws.maxlenleft<=0 || rws.maxlenright<=0)
                        /* invalid setting */
				rws.maxlenleft=rws.maxlenright=0;
		}
		else if (strcmp(argv[arg], "-mo") == 0) {
			arg++;
			if (arg >= argc)
				badusage();
			rws.maxoverlaplen = atoi(argv[arg]);
		}
		else if (strcmp(argv[arg], "-mrl") == 0) {
			maxreducelenset = TRUE;
			arg++;
			if (arg >= argc)
				badusage();
			maxreducelen = atoi(argv[arg]);
		}
		else if (strcmp(argv[arg], "-rec") == 0) {
			orderingset = TRUE;
			rws.ordering = RECURSIVE;
		}
		else if (strcmp(argv[arg], "-rtrec") == 0) {
			orderingset = TRUE;
			rws.ordering = RT_RECURSIVE;
		}
		else if (strcmp(argv[arg], "-lex") == 0) {
			orderingset = TRUE;
			rws.ordering = SHORTLEX;	/* DEFAULT! */
		}
		else if (strcmp(argv[arg], "-wtlex") == 0) {
			orderingset = TRUE;
			rws.ordering = WTLEX;
		}
		else if (strcmp(argv[arg], "-wreath") == 0) {
			orderingset = TRUE;
			rws.ordering = WREATHPROD;
		}
		else if (strcmp(argv[arg], "-cn") == 0) {
			confnumset = TRUE;
			arg++;
			if (arg >= argc)
				badusage();
			confnum = atoi(argv[arg]);
                }
		else if (strcmp(argv[arg], "-op") == 0) {
			outputprefixes = TRUE;
                }
		else if (strcmp(argv[arg], "-rk") == 0) {
			arg++;
			if (arg >= argc)
				badusage();
			rws.rkminlen = atoi(argv[arg]);
			arg++;
			if (arg >= argc)
				badusage();
			rws.rkmineqns = atoi(argv[arg]);
			if (rws.rkminlen<=0 || rws.rkmineqns<0)
                        /* invalid setting */
				rws.rkminlen=rws.rkmineqns=0;
		}
		else if (strcmp(argv[arg], "-wd") == 0) {
			worddiffs = TRUE;
		}
		else if (strcmp(argv[arg], "-mwd") == 0) {
			maxwdiffsset = TRUE;
			arg++;
			if (arg >= argc)
				badusage();
			maxwdiffs = atoi(argv[arg]);
                }
		else if (strcmp(argv[arg], "-silent") == 0) {
			silentset = TRUE;
			print_level = 0;
		}
		else if (strcmp(argv[arg], "-v") == 0) {
			print_level = 2;
			verboseset = TRUE;
		}
		else if (strcmp(argv[arg], "-vv") == 0) {
			print_level = 3;
			verboseset = TRUE;
		}
		else if (strcmp(argv[arg], "-mh") == 0) {
			arg++;
			if (arg >= argc)
				badusage();
			maxhad = atoi(argv[arg]);
		}
		else if (strcmp(argv[arg], "-mt") == 0) {
			arg++;
			if (arg >= argc)
				badusage();
			min_time = atoi(argv[arg]);
                }
		else if (strcmp(argv[arg], "-hf") == 0) {
			arg++;
			if (arg >= argc)
				badusage();
			halting_factor = atoi(argv[arg]);
                }
		else {
			if (argv[arg][0] == '-')
				badusage();
			strcpy(inf, argv[arg]);
			strcpy(outf, inf);
			strcat(outf, ".");
			strcat(outf, PROGNAME);
			strcpy(outfr, inf);
			strcat(outfr,".reduce");
			if (worddiffs){
				strcpy(outf1, inf);
				strcat(outf1, ".");
				strcat(outf1, "diff1");
				strcpy(outf2, inf);
				strcat(outf2, ".");
				strcat(outf2, "diff2");
				strcpy(outflog, inf);
				strcat(outflog, ".");
				strcat(outflog, "log");
			}
		}
		arg++;
	}
	if (strlen(inf)==0)
		badusage();

	if (resume) {
	    if ((rfile=fopen(outf,"r"))==0){
	  	fprintf(stderr,"Error: cannot open file %s\n",inf);
    		exit(1);
	    }
	    read_kbinput(rfile,FALSE);
/* If we are re-running from previous output, validity checks on input
 * should not be necessary.
 */
	}
	else {
	    if ((rfile=fopen(inf,"r"))==0){
	  	fprintf(stderr,"Error: cannot open file %s\n",inf);
    		exit(1);
	    }
	    read_kbinput(rfile,TRUE);
	}
	fclose(rfile);

	if (resume_with_orig) {
	    if ((rfile=fopen(inf,"r"))==0){
	  	fprintf(stderr,"Error: cannot open file %s\n",inf);
    		exit(1);
	    }
	    read_extra_kbinput(rfile,FALSE);
	    fclose(rfile);
	}

	if (worddiffs){
		if (rws.ordering != SHORTLEX) {
        		fprintf(stderr,
  "Word-difference calculation only makes sense for the shortlex ordering.\n");
                        exit(1);
                }
		if (!confnumset)
			confnum *=8;
/* If we are calculating word-differences, we are not interested in
 * lots of confluence tests.
 */
		for (i=1;i<=rws.num_gens;i++) if (rws.inv_of[i]==0){
			fprintf(stderr,
   "Word-difference calculation requires that all generators have inverses.\n");
			exit(1);
		}
		initialise_wd_fsa(&wd_fsa,reduction_fsa.alphabet,maxwdiffs);
		tmalloc(wd_record,struct wdr,MAXCYCLES+1);
		hadct = 0;
		tot_eqns = 0;
	}

	tidyon = tidyon_now = FALSE;

	/* Now we are ready to start the algorithm. The policy is that
	 * whenever a new relation is found, it is immediately compared with
	 * all of the inverse relations. This seems to be most efficient. It
	 * might also be a good idea to immediately compare new relations
	 * with any other very short relations. Otherwise, a new relation has
	 * to wait in the queue before being compared with the others.
         */


	signal(SIGINT, interrupt);
	signal(SIGKILL, interrupt);
	signal(SIGQUIT, interrupt);
        tmalloc(rws.history,int,maxreducelen);
        tmalloc(rws.slowhistory,int *,maxreducelen);
        tmalloc(rws.slowhistorysp,int,maxslowhistoryspace);

	reduce_word = slow_rws_reduce;
	check_reduce_word = slow_check_rws_reduce;
	for (k = rws.num_inveqns + 1; k <= rws.num_eqns; k++) {
		for (l = 1; l <= rws.num_inveqns; l++) {
			consider(k, l);
			if (tidyon_now)
				break;
			consider(l, k);
			if (tidyon_now)
				break;
		}
		if (tidyon_now)
			break;
	}

	if (onintr)
		tidyon = tidyon_now = TRUE;
/* When tidy_on becomes true, we are due to tidy up, but we deal with
 * overlaps with inverse equations first - this can avoid an excessive
 * number of word-differences which later trun out to be redundant.
 * However, if tidy_on now becomes true, we have too many new equations,
 * and we tidy up immediately.
 */

restart:
	nneweqns = 0;
	cct = 0;
	for (i = rws.num_inveqns + 1; i <= rws.num_eqns; i++) {
		onum_eqnsi = rws.num_eqns;
		for (j = 1; j <= i; j++) {
			if (rws.eqns[i].done && rws.eqns[j].done)
				continue;
			if (tidyon)
				break; /* only there for when i = 1 */
			cct++;
			onum_eqns = rws.num_eqns;
			consider(i, j);
			if (tidyon_now)
				break;
			if (i != j)
				consider(j, i);
			if (tidyon_now)
				break;
			if (rws.num_eqns > onum_eqns)
				/* Compare new relations with inverse
				 * relations */
				for (k = onum_eqns + 1;k <= rws.num_eqns; k++) {
					for (l = 1; l <= rws.num_inveqns; l++) {
						consider(k, l);
						if (tidyon_now)
							break;
						consider(l, k);
						if (tidyon_now)
							break;
					}
					if (tidyon_now)
						break;
				}
			if (tidyon)
				break;
		}
		/* If we have more than tidyint new relations since last
		 * tidying up, then tidy; that is, go through all relations,
		 * remove the redundant ones, and then redefine the FSA. This
		 * will usually reduce its number of states. */
		if (onintr) {
			halting = tidyon = TRUE;
			exit_status = 2;
		}
		if (tidyon) {
			i = tidyup(i) - 1;
			if (halting) {
				if (!worddiffs) {
					build_fulltable();
					reduce_word = rws_reduce;
				}
				output_and_exit(TRUE);
			}
			if (worddiffs){
				if (wd_reorder)
					tmalloc(new_wd,boolean,rws.num_eqns+1)
				else new_wd = 0;
				build_wd_fsa(&wd_fsa,rws,new_wd);
				should_we_halt();
				if (halting)
					output_and_exit(FALSE);
				if (wd_fsa.states->size > maxwdiffs/2) {
				    char **newnames;
				    if (print_level>=3)
					printf("    #Doubling maxwdiffs.\n");
				    maxwdiffs*=2;
				    tmalloc(newnames,char *,maxwdiffs+1);
				    for (j=1;j<=wd_fsa.states->size;j++)
					newnames[j] = wd_fsa.states->names[j];
				    tfree(wd_fsa.states->names);
				    wd_fsa.states->names = newnames;
				    tfree(wd_fsa.table->table_data_ptr[0]);
				    tfree(wd_fsa.table->table_data_ptr);
  				    fsa_table_init(wd_fsa.table,maxwdiffs,
					wd_fsa.alphabet->size);
				    tfree(wd_fsa.table->table_data_dptr);
				    fsa_table_dptr_init(&wd_fsa);
				}
				clear_wd_fsa(&wd_fsa);
				if (wd_reorder){
					i = wd_sort_eqns(i) - 1;
					tfree(new_wd);
				}
			}
			tidyon = tidyon_now = FALSE;
                        nneweqns = 0;
		}
		else if (onum_eqnsi==rws.num_eqns) {
/* Number of equations has not increased for this i.
*/
			if (confnum>0 && cct >=confnum) {
				do_conf_test = TRUE;
				i = tidyup(i);
        if (print_level>=2)
	  printf("  #No new eqns for some time - testing for confluence\n");
				k = conf_check();
				if (k==1) { 
					if (hadlongoverlap && print_level>0) {
						printf(
 "#System is not guaranteed confluent, since not all overlaps were processed.\n"
						);
					}
					else {
        				    if (print_level>0)
					      printf("#System is confluent.\n");
					    rws.confluent = TRUE;
					}
					output_and_exit(TRUE);
				}
				else if (k== -1) { 
        				if (print_level>0)
	  					printf(
    "#System is not confluent - halting because new equations are too long.\n");
					output_and_exit(TRUE);
				}
				else if (halting)
					output_and_exit(TRUE);
				do_conf_test = FALSE;
				cct = 0;
			}
		}
		else
			cct = 0;
		rws.eqns[i].done = TRUE;
	}

	/* Algorithm has now completed. */
	if (print_level>=2)
		printf("  #Search for overlaps is complete.\n");
          
	halting = TRUE;
        tidyup(rws.num_eqns);
	k = conf_check();
	if (k==1){
		if (hadlongoverlap && print_level>0) {
			printf(
 "#System is not guaranteed confluent, since not all overlaps were processed.\n"
			);
		}
		else {
       		    if (print_level>0)
		      printf("#System is confluent.\n");
		    rws.confluent = TRUE;
		}
        }
	else if (k== -1 && print_level>0) 
	  printf(
   "#System is not confluent - halting because new equations are too long.\n");
	else if (k==0) {
	  printf(
           "#Search for overlaps is complete, but system is not confluent.\n");
	  if (!onintr) {
                printf("#Let's try again then!\n");
		halting = FALSE;
		build_quicktable();
		reduce_word=slow_rws_reduce;
		goto restart;
	  }
	}
	output_and_exit(TRUE);
}

int
insert(lhs,rhs)
	char             **lhs, **rhs;

/* Look at the words in rws.testword1 and rws.testword2, and remove any
 * common prefixes and suffixes of generators with inverses.
 * If they are not both empty at that stage, then use procedure
 * "compare" to see which word comes first in the ordering.
 * If using shortlex order, and length(bigger) - length(smaller) > 2, then
 * use inverses of generators to transfer generators from testword1 to
 * testword2.  Finally, copy testword1 and testword2 into *lhs, and *rhs
 * (after assigning space).
 * Return value is 0 if lhs and rhs are equal, 1 if the new equation is
 * inserted, and -1 if it is not inserted because it is too long
 * according to maxlenleft and maxlenright.
 */
{
	int             bigger;
	char           *ptr1,
	               *ptr2,
	               *t1,
	               *t2,
	               *ts,
	               *ptre1,
	               *ptre2;

	ptr1 = rws.testword1;
	ptr2 = rws.testword2;
	if (strcmp(ptr1,ptr2)==0)
		return 0;
	ptre1 = ptr1 + strlen(ptr1);
	ptre2 = ptr2 + strlen(ptr2);
	while (ptr1 < ptre1)
		if (*ptr1 == *ptr2 && rws.inv_of[*ptr1])
			*(ptr1++) = *(ptr2++) = 0;
		else
			break;
	if (ptr1 == ptre1 && ptr2 == ptre2)
		return 0;
	if (*ptr1 && *ptr2)
		while (*(ptre1-1) == *(ptre2-1)  && rws.inv_of[*(ptre1-1)] &&
			ptre1 > ptr1 && ptre2 > ptr2) {
			ptre1--;
			ptre2--;
		}
	t1 = ptr1;
	t2 = ptr2;
	*ptre1 = *ptre2 = '\0';

	/* Common prefixes and suffixes have now been stripped.
	 * Now make t1 the larger.  */
	bigger = compare(t1, t2);
	if (bigger == 2) {
		ts = t1;
		t1 = t2;
		t2 = ts;
		ts = ptre1;
		ptre1 = ptre2;
		ptre2 = ts;
	}
	/* Now transfer terms from lhs to rhs if possible */
	if (rws.ordering == SHORTLEX){
	    while (strlen(t1) > strlen(t2)+2) {
                if (rws.inv_of[*(ptre1-1)]==0)
                  break;
		ptre1--;
		*ptre2 = rws.inv_of[*ptre1];
		*ptre1 = '\0';
		*(++ptre2) = '\0';
	    }
	    if (strlen(t1)==strlen(t2)+2 && t2[0] && t1[0]>t2[0] &&
			rws.inv_of[*(ptre1-1)]){
                ptre1--;
                *ptre2 = rws.inv_of[*ptre1];
                *ptre1 = '\0';
                *(++ptre2) = '\0';
            }
	}
	if (rws.maxlenleft>0 &&
              (strlen(t1)>rws.maxlenleft || strlen(t2)>rws.maxlenright)) 
		return -1;

	/* t1 is biggest according to the ordering, so it becomes the LHS of
	 * the new relation.
         */
	tmalloc(*lhs, char, strlen(t1) + 1)
		ptr1 = t1;
	ptr2 = *lhs;
	while (ptr1 < ptre1)
		*(ptr2++) = *(ptr1++);
	*(ptr2++) = '\0';
	tmalloc(*rhs, char, strlen(t2) + 1)
		ptr1 = t2;
	ptr2 = *rhs;
	while (ptr1 < ptre2)
		*(ptr2++) = *(ptr1++);
	*(ptr2++) = '\0';

	return 1;
}


void
consider(k, l)
	int             k,
	                l;

/* The left hand sides of the relations k and l are considered for common
 * parts, according to the Knuth-Bendix procedure, and any new relations
 * produced are processed.
 * We look for overlaps of either a suffix of equation k with a prefix of
 * equation l, or of a subword of eqn. k with the whole of eqn. l.
 */
{
	char           *midwd,
	               *ptr1,
	               *ptr2,
	               *ptr,
	               *testwd1,
	               *testwd2;
	boolean		ok;
	int		r;

	midwd = testwd1 = rws.eqns[k].lhs;
	testwd2 = rws.eqns[l].lhs;
	while (*midwd != '\0') {
		if (rws.maxoverlaplen>0 &&
			(midwd-testwd1)+strlen(testwd2) > rws.maxoverlaplen) {
		/* Any further overlap would exceed the permitted length */
			break;
		}
		ptr1 = midwd;
		ptr2 = testwd2;
		ok = TRUE;
		while (*ptr1 != '\0' && *ptr2 != '\0') {
			if (*ptr1 != *ptr2) {
				ok = FALSE;
				break;
			}
			ptr1++; ptr2++;
		}
		if (ok)
			/* An overlap has been found. The new equation
			 * for consideration is inserted into words 
			 * rws.testword1 and rws.testword2. */
		{
			if (*ptr1 == '\0') {
				ptr1 = testwd1 - 1;
				ptr = rws.testword1 - 1;
				while (++ptr1 < midwd)
					*(++ptr) = *ptr1;
				ptr1 = rws.eqns[l].rhs - 1;
				while (*(++ptr1) != '\0')
					*(++ptr) = *ptr1;
				*(++ptr) = '\0';
				ptr1 = rws.eqns[k].rhs - 1;
				ptr = rws.testword2 - 1;
				while (*(++ptr1) != '\0')
					*(++ptr) = *ptr1;
				ptr1 = ptr2 - 1;
				while (*(++ptr1) != '\0')
					*(++ptr) = *ptr1;
				*(++ptr) = '\0';
			} else {
				ptr2 = testwd1 - 1;
				ptr = rws.testword1 - 1;
				while (++ptr2 < midwd)
					*(++ptr) = *ptr2;
				ptr2 = rws.eqns[l].rhs - 1;
				while (*(++ptr2) != '\0')
					*(++ptr) = *ptr2;
				ptr2 = ptr1 - 1;
				while (*(++ptr2) != '\0')
					*(++ptr) = *ptr2;
				*(++ptr) = '\0';
				ptr2 = rws.eqns[k].rhs - 1;
				ptr = rws.testword2 - 1;
				while (*(++ptr2) != '\0')
					*(++ptr) = *ptr2;
				*(++ptr) = '\0';
			}
			/* Now reduce rws.testword1 and rws.testword2,
			 * using the current relations, and then install
			 * them if they are different. */
			reduce_word(rws.testword1);
			reduce_word(rws.testword2);
			if (insert( &(rws.eqns[rws.num_eqns+1].lhs),
                                    &(rws.eqns[rws.num_eqns+1].rhs) )>0){
/* We have a new equation */
			    r = modify_table(rws.num_eqns + 1);
			    if (r == -1) {
				tfree(rws.eqns[rws.num_eqns + 1].lhs)
				tfree(rws.eqns[rws.num_eqns + 1].rhs)
                                return;
			    }
			    if (r==1) {
				    rws.num_eqns++;
				    tot_eqns++;
				    rws.eqns[rws.num_eqns].done = FALSE;
				    if (rws.num_eqns >= maxeqns){
		if (print_level>0)
			printf("#Maximum number of equations exceeded.\n");
					hadmaxeqns = TRUE;
					hadct++;
					tidyon = tidyon_now = TRUE;
					if (!worddiffs || hadct>maxhad){
					    halting = TRUE;
					    exit_status = 2;
					}
			/* If we are calculating word-differences we do not
			 * always need to stop.
			 */
					return;
				    }
				/* Decide if it is time to tidy up. */
				    if (++nneweqns == tidyint)
					    tidyon = TRUE;
				    if (nneweqns == 2*tidyint)
					    tidyon_now = TRUE;
			    }
			}
		}
		midwd++;
	}
}

int
tidyup(crelno)
	int             crelno;

/* Remove redundant relations. "crelno" is the current relation being
 * processed.
 * Return value is the new value of crelano, which may have changed as
 * a result of earlier equations being eliminated.
 */
{
	int             i,
			iv,
	                nnum_eqns,
			totlenl,
			totlenr,
	                ret;
	boolean         moving,
	                retain,
			red,
	                some_changed;
	char           **newlhs,
	               **newrhs,
			*testword1 = rws.testword1,
			*testword2 = rws.testword2;

	lasttid = rws.num_eqns;
	ret = crelno;
repeat:
	tmalloc(newlhs, char *, rws.num_eqns + 1)
	tmalloc(newrhs, char *, rws.num_eqns + 1)
	for (i = 1; i <= rws.num_eqns; i++)
		rws.eqns[i].changed = rws.eqns[i].eliminated = FALSE;
	nnum_eqns = 0;
	moving = FALSE;
	some_changed = FALSE;
	totlenl = totlenr = 0;
	for (i = 1; i <= rws.num_eqns; i++) {
		strcpy(testword1, rws.eqns[i].lhs);
		strcpy(testword2, rws.eqns[i].rhs);
		red = check_reduce_word(testword1, i);
		if (red) {
		/* LHS is irreducible by other rws.eqns */
			retain = TRUE;
			reduce_word(testword2);
			if (strcmp(testword2, rws.eqns[i].rhs))
		/* RHS is reducible */
				some_changed = rws.eqns[i].changed = TRUE;
			totlenl += strlen(testword1);
			totlenr += strlen(testword2);
		} else {
		/* LHS can be reduced using other equations */
			reduce_word(testword1);
			reduce_word(testword2);
			retain = strcmp(testword1, testword2);
			if (retain) {
				some_changed = rws.eqns[i].changed = TRUE;
				totlenl += strlen(testword1);
				totlenr += strlen(testword2);
			}
		}
		if (rws.eqns[i].changed){
			iv = insert(newlhs+i,newrhs+i);
			retain = iv > 0;
			if (iv < 0)
				lostinfo = TRUE;
			if (retain)
				rws.eqns[i].done = FALSE;
		}
	
		if (!retain) {
			moving = TRUE;
			rws.eqns[i].eliminated = TRUE;
		}
	}
	if (moving || some_changed) {
		moving = FALSE;
		for (i = 1; i <= rws.num_eqns; i++)
			if (!rws.eqns[i].eliminated) {
				nnum_eqns++;
				if (rws.eqns[i].changed) {
					fflush(stdout);
					tfree(rws.eqns[i].lhs)
                                        tfree(rws.eqns[i].rhs)
					rws.eqns[i].lhs = newlhs[i];
					rws.eqns[i].rhs = newrhs[i];
				}
				if (moving)
					/* "moving" means that some previous
					 * relations have not been retained,
					 * and so we have to renumber words. */
				{
				  fflush(stdout);
				  rws.eqns[nnum_eqns].lhs = rws.eqns[i].lhs;
				  rws.eqns[nnum_eqns].rhs = rws.eqns[i].rhs;
				  rws.eqns[nnum_eqns].done = rws.eqns[i].done;
				  rws.eqns[i].lhs = rws.eqns[i].rhs = 0;
				  rws.eqns[i].done = FALSE;
				}
				if (i <= crelno)
					ret = nnum_eqns;
			} else {
				moving = TRUE;
				if (i <= rws.num_inveqns)
					rws.num_inveqns--;
				fflush(stdout);
				tfree(rws.eqns[i].lhs)
				tfree(rws.eqns[i].rhs)
			}
	} else
		nnum_eqns = rws.num_eqns;
	/* If any relations have been changed or eliminated, then we
	 * completely reconstruct the FSA from scratch.
	 * This is not too time consuming provided that we do not
	 * tidy up too often.
         */
	if (num_states > current_maxstates) { /* shouldn't happen */
		if (print_level>0)
			printf(
			"#Maximum number of states exceeded (why?).\n");
		halting = TRUE;
	        exit_status = 2;
	}
	tfree(newlhs) tfree(newrhs)
	if (moving || some_changed || rk_on) {
		tot_eqns -= (rws.num_eqns - nnum_eqns);
		rws.num_eqns = nnum_eqns;
		if (rws.rkmineqns>0) {
			if (!rk_on && !halting && rws.num_eqns>=rws.rkmineqns) {
				rk_on = TRUE;
				reduce_word = slow_rws_reduce_rk;
				check_reduce_word=slow_check_rws_reduce_rk;
				oldconfnum = confnum; confnum = 0;
				rk_init(maxeqns);
			}
			else if (rk_on && (rws.num_eqns<rws.rkmineqns ||
						     halting && !worddiffs)) {
                                rk_on = FALSE;
                                reduce_word = slow_rws_reduce;
				check_reduce_word = slow_check_rws_reduce;
				confnum = oldconfnum;
				rk_clear();
                        }
		}
		build_quicktable();
/* On the last time, we need to be certain that the set of equations is reduced,
 * so we repeat.
 */
		if ((halting || do_conf_test) && !worddiffs)
			goto repeat;
	}
	times(&time_buffer);
	i = time_buffer.tms_utime / HZTIME;
        if (print_level>=2)
	  printf(
	      "  #%d eqns; total len: lhs, rhs = %d, %d; %d states; %d secs.\n",
		rws.num_eqns, totlenl, totlenr, num_states, i);

	return ret;
}

void
build_quicktable()
{
	int	i;
    restart:
	if (double_states) {
	    current_maxstates *= 2;
	    tfree(reduction_fsa.table->table_data_ptr[0]);
	    tfree(reduction_fsa.table->table_data_ptr);
	    fsa_table_init(reduction_fsa.table,current_maxstates,rws.num_gens);
            table = reduction_fsa.table->table_data_ptr;
	    tfree(rws.preflen);
	    tfree(rws.prefno);
	    tmalloc(rws.preflen,int,current_maxstates+1);
	    tmalloc(rws.prefno,int,current_maxstates+1);
	    double_states = FALSE;
	}
	if (rk_on)
		rk_reset(maxeqns);
	num_states = 1;
	maxpreflen=0;
	rws.preflen[1] = rws.prefno[1] = 0;
	for (i = 1; i <= rws.num_gens; i++)
		set_dense_target(table,i,1,0);
	for (i = 1; i <= rws.num_eqns; i++){
          if (modify_table(i) == -1){
		double_states = TRUE;
		goto restart;
          }
	}
}

int
modify_table(relno)
	int             relno;

/* This version returns a table which only rejects the left hand sides
 * of the rws.eqns themselves.
 * Return value is -1 if current_maxstates is exceeded -  otherwise 1.
 */
{
	int             genno,
			nextgenno,
			state,
			x,
			l,
			len,
	                i;
	char           *ptr;

	ptr = rws.eqns[relno].lhs;
	len = strlen(ptr);
	if (rk_on && len>=rws.rkminlen) {
		rk_add_lhs(relno);
		return 1;
	}
        if (len + num_states > current_maxstates){
/* The number of states of the table could exceed the maximum, so we give up 
 * immediately - if maxstates is set we stop completely, and otherwise we
 * double current_maxstates after next tidying.
 */
		if (current_maxstates==maxstates) {
                    if (print_level>0)
		  	printf(
                    	    "#Maximum number of states exceeded.\n");
			halting = TRUE;
			exit_status = 2;
		}
		else {
			if (print_level >= 3)
		  		printf(
		   "    #current_maxstates exceeded - will tidy and double.\n");
			double_states = TRUE;
		}
		tidyon = tidyon_now = TRUE;
		return -1;
	}
	state = 1;
	genno = *(ptr++);
	l=0;
	while ((nextgenno = *(ptr++))) {
		l++;
		x = dense_target(table,genno,state);
		if (x<0)
			return 1;
/* If that happens, the equation must have a prefix which is also a lhs.
 * We keep the equation for now - it will be reduced on the next tidying.
 */
		if (x==0){
			num_states++;
			for (i = 1; i <= rws.num_gens; i++)
				set_dense_target(table,i,num_states,0);
			rws.preflen[num_states]=l;
			if (l>maxpreflen)
				maxpreflen=l;
			rws.prefno[num_states]=relno;
                        x = num_states;
			set_dense_target(table,genno,state,num_states);
		}
		state=x;
		genno=nextgenno;
	}
	set_dense_target(table,genno,state,-relno);
	return  1;
}

void
build_fulltable()
/* Extends the table to one rejecting all strings containing a lhs,
 * builds edges in E_2 - E_3 as in Sims' book, p. 118-9.
*/
{
	int	i,
		j,
		k,
		l;
	char	*w,
		*we;
	int *preflen = rws.preflen;
	int *prefno = rws.prefno;
	for (i=1;i<=rws.num_gens;i++)
	    if (dense_target(table,i,1)==0)
		 set_dense_target(table,i,1,1);
	for (l=1;l<=maxpreflen;l++) {
		for (j=2;j<=num_states;j++) if (preflen[j]==l) {
			w=rws.eqns[prefno[j]].lhs;
			we = w+l;
			k = 1;
			while (++w < we)
				k = dense_target(table,*w,k);
		/* k is the image of 1 under the prefix representing state  j,
		 * but with its first letter removed. */
			for (i=1;i<=rws.num_gens;i++)
			    if (dense_target(table,i,j)==0)
				set_dense_target(table,i,j,
                                                 dense_target(table,i,k));
		}
	}
}

int
lex_compare(w1,w2)
	char		*w1,
			*w2;

/* Compare words w1 and w2 to see which is bigger according to the ordering.
 * The ordering used here is longer words are bigger, and amongst equal
 * length words, lexicographical ordering according to the order of the
 * generators as input is used. (I use 'bigger' to mean greater in the
 * ordering.)  Returns  0 if w1=w2, 1 if w1 bigger than w2, otherwise 2.
 */
{
	int len1,len2,bigger=0;
	char *we1,*we2;

	len1=strlen(w1);
	len2=strlen(w2);
	we1=w1+len1;
	we2=w2+len2;
	if (len1>len2)
		bigger=1;
	else if (len2>len1)
		bigger=2;
	else
		while (w1<we1) {
			if (*w1>*w2) {
				bigger=1;
				break;
			} else if (*w2>*w1) {
				bigger=2;
				break;
			}
			w1++;
			w2++;
		}
	return bigger;
} 

int
wtlex_compare(w1,w2)
	char		*w1,
			*w2;

/* Compare words w1 and w2 to see which is bigger according to the ordering.
 * The ordering used here is longer words are bigger, where length is computed
 * by adding up the weights of the generators in the words, and amongst equal
 * length words, lexicographical ordering according to the order of the
 * generators as input is used. (I use 'bigger' to mean greater in the
 * ordering.)  Returns  0 if w1=w2, 1 if w1 bigger than w2, otherwise 2.
 */
{
	int len1,len2, wtlen1, wtlen2, bigger=0;
	char *wb1,*wb2,*we1,*we2;

	len1=strlen(w1); len2=strlen(w2);
	wb1=w1; wb2=w2;
	we1=w1+len1; we2=w2+len2;
	wtlen1 = wtlen2 = 0;
	while (wb1 < we1) {
		wtlen1 += rws.weight[*wb1];
		wb1++;
	}
	while (wb2 < we2) {
		wtlen2 += rws.weight[*wb2];
		wb2++;
	}
	if (wtlen1>wtlen2)
		bigger=1;
	else if (wtlen2>wtlen1)
		bigger=2;
	else
		while (w1<=we1) {
			if (*w1>*w2) {
				bigger=1;
				break;
			} else if (*w2>*w1) {
				bigger=2;
				break;
			}
			w1++; w2++;
		}
	return bigger;
} 

int
rec_compare(w1, w2)
	char           *w1,
	               *w2;

/* Compare words w1 and w2 to see which is 'bigger' according to the
   ordering. The ordering used here is recursive path ordering (based on
   that described in the book "Confluent String Rewriting" by Matthias Jantzen,
   Defn 1.2.14, page 24).
   Returns 1 if w1 is bigger than w2, 2 if w2 is bigger than w1,
   0 if equal.
   ----------------------------------------------------------------------
   The ordering is as follows:

   let u, v be elements of X*
   u >= v iff one of the following conditions is fulfilled;

   1) u = v			OR
   u = u'a, v = v'b for some a,b elements of X, u',v' elements of X*
   and then:
   2) a = b and u' >= v'		OR
   3) a > b and u  > v'		OR
   4) b > a and u'> v
   ----------------------------------------------------------------------
   Written by : Jamie P. Curmi (February 1992) 
   Corrected by : dfh (December 1993) (must return 0 if words are equal).
   Altered by dfh (4.6.94.) to reverse word direction in order,
   and to eliminate recursion.
*/
{
        int     lastmoved=0;
        char    *p1,
                *p2;
        p1 = w1 + strlen(w1) - 1;
        p2 = w2 + strlen(w2) - 1;
        while (1) {
                if (p1 < w1){
                        if (p2 < w2)
                                return lastmoved;
                        return 2;
                }
                if (p2 < w2)
                        return 1;
                if (*p1 == *p2){
                        p1--;
                        p2--;
                }
                else if (*p1 < *p2){
                        p1--;
                        lastmoved = 1;
                }
                else if (*p2 < *p1){
                        p2--;
                        lastmoved = 2;
                }
        }
}

int
rt_rec_compare(w1, w2)
	char           *w1,
	               *w2;

/* Compare words w1 and w2 to see which is 'bigger' according to the
   ordering. The ordering used here is recursive path ordering (based on
   that described in the book "Confluent String Rewriting" by Matthias Jantzen,
   Defn 1.2.14, page 24).
   Returns 1 if w1 is bigger than w2, 2 if w2 is bigger than w1,
   0 if equal.
   ----------------------------------------------------------------------
   The ordering is as follows:

   let u, v be elements of X*
   u >= v iff one of the following conditions is fulfilled;

   1) u = v			OR
   u = a'u, v = bv' for some a,b elements of X, u',v' elements of X*
   and then:
   2) a = b and u' >= v'		OR
   3) a > b and u  > v'		OR
   4) b > a and u'> v
   ----------------------------------------------------------------------
   Written by : Jamie P. Curmi (February 1992) 
   Corrected by : dfh (December 1993) (must return 0 if words are equal).
*/
{
        int     lastmoved=0;
        char    *p1,
                *p2;
        p1 = w1 + strlen(w1) - 1;
        p2 = w2 + strlen(w2) - 1;
        while (1) {
                if (p1 < w1){
                        if (p2 < w2)
                                return lastmoved;
                        return 2;
                }
                if (p2 < w2)
                        return 1;
                if (*w1 == *w2){
                        w1++;
                        w2++;
                }
                else if (*w1 < *w2){
                        w1++;
                        lastmoved = 1;
                }
                else if (*w2 < *w1){
                        w2++;
                        lastmoved = 2;
                }
        }
}

int
wreath_compare(w1,w2)
        char *w1, *w2;
/* Compare w1 and w2 to see which comes first in the wreath-product ordering
 * (as defined in Sims' book), using the level function rws.level.
 * Note that the recursive ordering is the special case of this with
 * levels 1,2,3,...
 * Returns  0 if w1=w2, 1 if w1 is bigger than w2, otherwise 2.
 */
{
        int     winning=0,
		waiting=0,
		winning_level=0,
		level1, level2;
        char    *p1,
                *p2;
        p1 = w1 + strlen(w1) - 1;
        p2 = w2 + strlen(w2) - 1;
        while (1) {
                if (p1 < w1){
                        if (p2 < w2)
                                return winning;
			if (winning!=1 || rws.level[*p2]>=winning_level)
                        	return 2;
			p2--;
                }
                if (p2 < w2) {
			if (winning!=2 || rws.level[*p1]>=winning_level)
                        	return 1;
			p1--;
		}
		level1 = rws.level[*p1];
		level2 = rws.level[*p2];
		if (waiting==1) {
			if (level2>level1) {
				waiting=2;
				winning=1; winning_level=0;
				p1--;
			} else if (level2==level1) {
				waiting=0;
				if (*p1>*p2) {
					winning=1; winning_level=level1;
				}
				else if (*p2>*p1) {
					winning=2; winning_level=level1;
				}
				else if (level1>winning_level)
					winning_level=0;
				p1--; p2--;
			}
			else {
				if (winning==1 && level2>=winning_level) {
					winning=2; winning_level=0;
				}
				p2--;
			}
                }
		else if (waiting==2) {
			if (level1>level2) {
				waiting=1;
				winning=2; winning_level=0;
				p2--;
			} else if (level2==level1) {
				waiting=0;
				if (*p2>*p1) {
					winning=2; winning_level=level2;
				}
				else if (*p1>*p2) {
					winning=1; winning_level=level2;
				}
				else if (level2>winning_level)
					winning_level=0;
				p1--; p2--;
			}
			else {
				if (winning==2 && level1>=winning_level) {
					winning=1; winning_level=0;
				}
				p1--;
			}
                }
		else if (level2>level1 && level2>=winning_level) {
			waiting=2;
			if (winning==2 && level1>=winning_level)  {
				winning=1; winning_level=0;
			}
			p1--;
		}
		else if (level1>level2 && level1>=winning_level) {
			waiting=1;
			if (winning==1 && level2>=winning_level)  {
				winning=2; winning_level=0;
			}
			p2--;
		}
		else {
			if (*p1>*p2 && level1>=winning_level) {
				winning=1; winning_level=level1;
			}
			else if (*p2>*p1 && level1>=winning_level) {
				winning=2; winning_level=level1;
			}
			else if (*p1==*p2 && level1>winning_level)
				winning_level=0;
			p1--; p2--;
		}
        }
}

int
compare(w1,w2)
/* COMPARE: Compares two words 'w1' and 'w2' to see which is 'bigger'
 * according to the ordering to be used (Set by the global 'ordering').
 * If more ordering options are to be provided, the new function should
 * be given a different name and called from a 'case' statement in this
 * function.
 * Returns  0 if w1=w2, 1 if w1 is bigger than w2, otherwise 2.
 */
	char *w1, *w2;
{
	switch (rws.ordering) {

		case RECURSIVE:		/* Recursive path ordering */
			return (rec_compare(w1,w2));

		case RT_RECURSIVE:		/* Recursive path ordering */
			return (rt_rec_compare(w1,w2));

		case WTLEX:		/* Weighted shortlex ordering */
			return (wtlex_compare(w1,w2));

		case WREATHPROD:		/* Wreath product ordering */
			return (wreath_compare(w1,w2));

		default: 			/* shortlex ordering */
			return (lex_compare(w1, w2));
	}
}

int
conf_check()
/* Performs fast check for confluence.
 * Return value is 1 if confluent, 0 if not confluent and new equations added,
 * and -1 of not confluent but no new equations added because they are longer
 * than allowed by the values of maxlenleft and maxlenright.
 * Algorithm used is index_confluent, from p.117 of Sims' book.
 */
{
	int	i,
		st,
		la,
		la2,
		lb,
		r,
		iv,
		rv,
		beta[1000],
		betar,
		y[1000],
		num_neweqns=0;
	char	*w,
		*wa,
		*wb,
		*wx,
		*we,
		*testword1 = rws.testword1,
		*testword2 = rws.testword2;
	boolean discarded=FALSE,
		bt;
	int *preflen = rws.preflen,
	    *prefno = rws.prefno;
	build_fulltable();
	reduce_word = rws_reduce;
	hadlongoverlap = FALSE;

	for (i=1;i<=rws.num_eqns;i++) {
		wa = rws.eqns[i].lhs;
		la = strlen(rws.eqns[i].lhs);
		la2 = strlen(rws.eqns[i].rhs);
		w = wa;
		we = w + la;
		st=1;
		while (++w < we)
			st = dense_target(table,*w,st);
		beta[0] = st;
		r = 0;
		bt = FALSE;
		while (r>=0 && !bt) {
			if (rws.maxoverlaplen>0 && r+la>rws.maxoverlaplen) {
		/* Any further overlap would exceed the permitted length */
				hadlongoverlap = TRUE;
				break;
			}
			betar = beta[r];
			bt = ((betar<0 && strlen(rws.eqns[-betar].lhs)<=r) ||
				(betar>=0 && preflen[betar]<=r));
			if (betar<0 && !bt) {
				strcpy(testword1,rws.eqns[i].rhs);
				wx = testword1 + la2;
				wb = rws.eqns[-betar].lhs;
				lb = strlen(wb);
				w = wb + lb - r;
				we = wb + lb;
				while (w < we)
					*(wx++) = *(w++);
				*(wx) = 0;
				strcpy(testword2,rws.eqns[i].lhs);
				wx = testword2 + la - lb + r;
				strcpy(wx,rws.eqns[-betar].rhs);
				reduce_word(testword1);
				reduce_word(testword2);
				if ((iv=insert(&(rws.eqns[rws.num_eqns+1].lhs),
					   &(rws.eqns[rws.num_eqns+1].rhs)))>0){
					rws.eqns[++rws.num_eqns].done=FALSE;
					num_neweqns++;
					if (num_neweqns>=tidyint ||
							rws.num_eqns==maxeqns) {
                		if (print_level>=2)
				  printf("  #System is not confluent.\n");
					  tidyon = TRUE;
					  if (rws.num_eqns==maxeqns) {
                				if (print_level>0)
                        printf("#Maximum number of equations exceeded.\n");
                                        	hadmaxeqns = TRUE;
                                        	hadct++;
                                        	tidyon_now = TRUE;
                                        	if (!worddiffs || hadct>maxhad){
                                            	  halting = TRUE;
                                            	  exit_status = 2;
                                        	}
					  }
					  if (!halting) {
					      build_quicktable();
					      reduce_word = slow_rws_reduce;
					  }
					  return 0;
					}
				}
				if (iv== -1)
/* We have discarded the equation because it is too long.
 * So we must remember that the system is not confluent.
 */
					discarded = TRUE;
				bt = TRUE;
			}
			if (!bt) {
				beta[r+1] = dense_target(table,1,betar);
				y[r+1] = 2;
				r++;
			}
			else {
				while (bt && r>0) {
					if (y[r] <= rws.num_gens) {
					    beta[r] = dense_target(table,y[r],
                                                                     beta[r-1]);
					    y[r]++;
					    bt = FALSE;
					}
					else
						r--;
				}
			}
		}
	}
	if (num_neweqns>0) {
 		if (print_level>=2)
               	  printf("  #System is not confluent (finished conf_check).\n");
		build_quicktable();
		reduce_word = slow_rws_reduce;
		return 0;
	}

	halting = TRUE;
	if (discarded) {
/* We will mark all equations as not done, since we may want to resume with
 * higher limits on equation length.
 */
		for (i=1;i<=rws.num_eqns;i++)
			rws.eqns[i].done = FALSE;
		exit_status = 2;
		return -1;
	}
 
	exit_status = 0;
	return 1;
}

void
sort_eqns(n)
	int n;
/* The equations are re-ordered into order of increasing LHS.
 * Only equations up to length n are output.
 */
{	int i, l, ct, maxl;

	char	 **newlhs,
	         **newrhs;
	boolean   *newdone;

	maxl=0;
        for (i=1;i<=rws.num_eqns;i++) if ((l=strlen(rws.eqns[i].lhs)) > maxl)
          maxl = l;
	if (n>0 && maxl>n && print_level>0){
		printf("#Warning - not all equations will be output.\n");
		if (rws.confluent)
			printf("#The output system will not be confluent.\n");
	}
	if (n==0 || n>maxl)
		n=maxl;
	tmalloc(newlhs, char *, rws.num_eqns + 1)
	tmalloc(newrhs, char *, rws.num_eqns + 1)
	tmalloc(newdone, boolean, rws.num_eqns + 1)
	ct = 0;
	for (l=1;l<=n;l++)
/* Find equations of which LHS has length i */
		for (i=1;i<=rws.num_eqns;i++) if (strlen(rws.eqns[i].lhs)==l){
			newlhs[++ct] = rws.eqns[i].lhs;
			newrhs[ct] = rws.eqns[i].rhs;
			newdone[ct] = rws.eqns[i].done;
		}
/* and free the rest */
	for (i=1;i<=rws.num_eqns;i++) if (strlen(rws.eqns[i].lhs)>n){
		tfree(rws.eqns[i].lhs);
		tfree(rws.eqns[i].rhs);
	}

        for (i=1;i<=rws.num_eqns;i++){
          rws.eqns[i].lhs = newlhs[i];
          rws.eqns[i].rhs = newrhs[i];
          rws.eqns[i].done = newdone[i];
        }
	tfree(newlhs);
	tfree(newrhs);
	tfree(newdone);
	rws.num_eqns = ct;
        build_quicktable();
        build_fulltable();
        if (print_level>0)
	  printf(
		"#Output: %d eqns; table has %d states.\n",
		rws.num_eqns, num_states);
}

int
wd_sort_eqns(x)
	int x;
/* The equations are re-ordered so that those that resulted in new entries
 * in the word-difference table come first.
 * Those are recorded in the boolean array new_wd.
 * x is the current equation number in the KB.
 */
{	int i,  j,  ct,  ret, newmaxeqns;

	char	 **newlhs,
	         **newrhs;
	boolean   *newdone;
	double r;

	if (hadmaxeqns){
/* In this case we will not keep all equations, to allow room for more */
		r = (double)x/rws.num_eqns;
		if (r > 0.75){
			if (print_level>=2)
				printf(" #Halting - r=%f, hadct=%d.\n",r,hadct);
			halting = TRUE;
			exit_status = 2;
			hadmaxeqns = FALSE;
		}
		else
			newmaxeqns = rws.num_eqns*(1+r)/2;
	}
	tmalloc(newlhs, char *, rws.num_eqns + 1)
	tmalloc(newrhs, char *, rws.num_eqns + 1)
	tmalloc(newdone, boolean, rws.num_eqns + 1)
	ct = 0;
	ret = 0;
/* Find equations that are already handled or gave new word-diffs */
	for (i=1;i<=rws.num_eqns;i++) if (rws.eqns[i].done || new_wd[i]){
			newlhs[++ct] = rws.eqns[i].lhs;
			newrhs[ct] = rws.eqns[i].rhs;
			newdone[ct] = rws.eqns[i].done;
			if (i==x) ret=ct;
		}
/* Now take the rest. */
	for (i=1;i<=rws.num_eqns;i++) if (!rws.eqns[i].done && !new_wd[i]){
			newlhs[++ct] = rws.eqns[i].lhs;
			newrhs[ct] = rws.eqns[i].rhs;
			newdone[ct] = rws.eqns[i].done;
			if (i==x) ret=ct;
			if (hadmaxeqns && ct>=newmaxeqns && ret>0) {
				if (print_level>0)
				  printf("Keeping only %d equations.\n",ct);
				/* free the rest */
				for (j=i+1;j<=rws.num_eqns;j++)
				  if (!rws.eqns[j].done && !new_wd[j]){
					tfree(rws.eqns[j].lhs);
					tfree(rws.eqns[j].rhs);
				  }
				rws.num_eqns = ct;
				hadmaxeqns = FALSE;
				break;
			}
		}

        for (i=1;i<=rws.num_eqns;i++){
          rws.eqns[i].lhs = newlhs[i];
          rws.eqns[i].rhs = newrhs[i];
          rws.eqns[i].done = newdone[i];
        }
	tfree(newlhs);
	tfree(newrhs);
	tfree(newdone);
        build_quicktable();
	return ret;
}

void
should_we_halt()
/* Try to decide whether we should halt, using number of word-differences. */
{ int	i, ndiff, t;
	times(&time_buffer);
	t = time_buffer.tms_utime / HZTIME;
	num_cycles++;
	if (num_cycles >= MAXCYCLES) {
		halting = TRUE;
		exit_status = 2;
	}
	wd_record[num_cycles].num_eqns = tot_eqns;
	wd_record[num_cycles].num_states = num_states;
	ndiff = wd_record[num_cycles].num_diff = wd_fsa.states->size;
	for (i=num_cycles-1;i>=1;i--) if (wd_record[i].num_diff < ndiff)
		break;
	i++;
	eqn_factor = wd_record[i].num_eqns==0 ? 0 :
				100*tot_eqns/wd_record[i].num_eqns - 100;
	states_factor =  hadct>0 ? 0 :
				100*num_states/wd_record[i].num_states - 100;
	if (halting_factor>0 && ((t>min_time && eqn_factor>=halting_factor
			 && (states_factor==0 || states_factor>=halting_factor))
	|| (eqn_factor>=2*halting_factor&&states_factor>=2*halting_factor))) {
		if (print_level>1) 
		    printf("  #eqn_factor=%d, states_factor=%d - halting.\n",
			eqn_factor,states_factor);
		halting = TRUE;
		exit_status = 0;
	}
}

void
output_and_exit(make_wd)
	boolean make_wd;
{
	int ndiff1, ndiff2;
	if (worddiffs){
		if (make_wd) {
			new_wd = 0;
			build_wd_fsa(&wd_fsa,rws,new_wd);
			should_we_halt(); /* to calculate factors */
		}
		ndiff1 = wd_fsa.states->size;
		wfile = fopen(outf1,"w");
		print_wdoutput(wfile,".diff1");
		fclose(wfile);
	  	make_full_wd_fsa(&wd_fsa,rws.inv_of,1);
		ndiff2 = wd_fsa.states->size;
		wfile = fopen(outf2,"w");
		print_wdoutput(wfile,".diff2");
		fclose(wfile);
		fsa_clear(&wd_fsa);
/* Writing info to a log-file - not sure if we want this at present.
		wfile = fopen(outflog,"w");
		fprintf(wfile,"eqn_factor:\t\t%d\n",eqn_factor);
		fprintf(wfile,"states_factor:\t\t%d\n",states_factor);
		fprintf(wfile,"diff1:\t\t%d\n",ndiff1);
		fprintf(wfile,"diff2:\t\t%d\n",ndiff2);
		fclose(wfile);
*/
		if (print_level>0) {
		  printf("#Halting with %d equations.\n",rws.num_eqns);
		  printf(
                   "#First word-difference machine with %d states computed.\n",
									ndiff1);
		  printf(
                   "#Second word-difference machine with %d states computed.\n",
									ndiff2);
		}
		if (print_level>=2) {
		  printf("  #eqn_factor=%d, states_factor=%d\n",
					eqn_factor,states_factor);
		}
		if (exit_status==0 && print_level>0)
	  printf("#System is confluent, or halting factor condition holds.\n");
		tfree(wd_record);
		if (rk_on)
			rk_clear();
	}
	else {
		if (rws.sorteqns)
			sort_eqns(rws.maxoplen);
		wfile = fopen(outf,"w");
		print_kboutput(wfile);
		fclose(wfile);
		wfile = fopen(outfr,"w");
		fsa_print(wfile,&reduction_fsa,rws.name);
		fclose(wfile);
		if (print_level>0)
		  printf("#Halting with %d equations.\n",rws.num_eqns);
	}
	rws_clear(&rws);
	fsa_clear(rws.reduction_machine);
	if (lostinfo && print_level>0) {
		printf(
"WARNING: The monoid defined by the presentation may have changed, since\n\
         equations have been discarded.\n\
         If you re-run, include the original equations.\n");
	}
	exit(exit_status);
}

void
badusage()
{
	fprintf(stderr, USAGE);
	exit(1);
}
