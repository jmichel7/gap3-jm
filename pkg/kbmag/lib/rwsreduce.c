/* file reduce.c - created 6/10/94 
 * This file contains the various routines for reducing words,
 * using a rewriting system.
 */

#include "defs.h"
#include "fsa.h"
#include "rws.h"
#include "externals.h"

extern int	maxslowhistoryspace,
		maxreducelen;
extern rewriting_system	rws;

/* Functions defined in this file: */
void rws_reduce();
void slow_rws_reduce();
boolean slow_check_rws_reduce();

void
rws_reduce(w)
	char     	      	*w;

/* Reduce "w", by replacing any occurrences of the LHS of the current
 * equations in the rewriting system by their RHS.
 * The reduction fsa for the rewriting system rws is used for this,and it is
 * assumed to be stored in dense format. The complete sequence of
 * states that it goes through on reading "w" is remembered in the array
 * "rws.history".
 */
{
	int             len,
	                st,
			**table,
	                longer,
			*history = rws.history;
	char           *midwd,
	               *ptr1,
	               *ptr2,
	               *ptr;

	table = rws.reduction_machine->table->table_data_ptr;
	midwd = w;
	len = 0;
	history[0] = 1;
	st = 1;
	while (*midwd != '\0') {
		st = history[++len] = dense_target(table,*midwd,st);
		if (st < 0)
			/* st= -n means that we have just read the LHS or
			 * relation  n. Replace it by RHS, and go back to the
			 * beginning of that subword. */
		{
			st = -st;
			ptr1 = midwd;
			ptr2 = rws.eqns[st].rhs - 1;
			if ((longer= strlen(rws.eqns[st].rhs)  -
                                     strlen(rws.eqns[st].lhs))  > 0)
			{
                                if (strlen(w)+longer>=maxreducelen){
              fprintf(stderr,
        "#Error: word too long in reduction - try increasing maxreducelen.\n");
				 exit(1);
				}
				ptr = w + strlen(w);
				while (ptr > midwd) {
					*(ptr + longer) = *ptr;
					ptr--;
				}
				ptr1 += longer;
			}
			len -= strlen(rws.eqns[st].lhs);
			midwd = w + len - 1;
			ptr = midwd;
			while (*(++ptr2) != '\0')
				*(++ptr) = *ptr2;
			if (ptr != ptr1) {
				while (*(++ptr1) != '\0')
					*(++ptr) = *ptr1;
				*(++ptr) = '\0';
			}
			st = history[len];
			/*	midwd--;*/
		}
		midwd++;
	}
}


void
slow_rws_reduce(w)
	char     	      	*w;
/* The version of reduce for the reduction automaton that recognises
 * left hand sides of relations only.
 * This is both slower, and needs more space, since the history has to
 * be stored as a list of subsets of the states, rather than a list of states.
 * The subsets are stored contiguously in the array rws.slowhistorysp, and
 * pointed at by rws.slowhistory.
 */
{
	int             len,
	                st,
			**table,
		       *histptr,
		       *histptre,
		       *nextptr,
			subrelno,
	                longer,
			*slowhistorysp = rws.slowhistorysp,
			**slowhistory = rws.slowhistory;
	char           *midwd,
	               *ptr1,
	               *ptr2,
	               *ptr;

restart:
	midwd = w;
	len = 0;
	table = rws.reduction_machine->table->table_data_ptr;
	nextptr=slowhistory[0]=slowhistorysp;
	while (*midwd != '\0') {
		subrelno=0;
		histptr=slowhistory[len++];
		slowhistory[len]=histptre=nextptr;
		if (nextptr+len-slowhistorysp > maxslowhistoryspace){
			if (print_level>=3)
  			  printf(
          	  "    #maxslowhistoryspace too small; doubling it.\n");
		  	tfree(rws.slowhistorysp);
			maxslowhistoryspace *= 2;
			tmalloc(rws.slowhistorysp,int,maxslowhistoryspace);
			slowhistorysp = rws.slowhistorysp;
			goto restart;
		}
/* First adjoin image of start state */
		st = dense_target(table,*midwd,1);
		if (st<0)
			subrelno= -st;
		else{
			if (st>0)
				*(nextptr++) = st;
			while (histptr<histptre){
				st = dense_target(table,*midwd,*(histptr++));
				if (st>0)
					*(nextptr++) = st;
				else if (st<0){
					subrelno= -st;
					break;
				}
			}
		}
		if (subrelno>0){
			/* subrelno=n means that we have just read the LHS of
			 * relation  n. Replace it by RHS, and go back to the
			 * beginning of that subword.
                         */
			ptr1 = midwd;
			ptr2 = rws.eqns[subrelno].rhs - 1;
			if ((longer = (strlen(rws.eqns[subrelno].rhs) -
                                       strlen(rws.eqns[subrelno].lhs)))>0) {
                                if (strlen(w)+longer>=maxreducelen){
              fprintf(stderr,
       "#Error: word too long in reduction - try increasing maxreducelen.\n");
				 exit(1);
				}
				ptr = w + strlen(w);
				while (ptr > midwd) {
					*(ptr + longer) = *ptr;
					ptr--;
				}
				ptr1 += longer;
			}
			len -= strlen(rws.eqns[subrelno].lhs);
			midwd = w + len - 1;
			ptr = midwd;
			while (*(++ptr2) != '\0')
				*(++ptr) = *ptr2;
			if (ptr != ptr1) {
				while (*(++ptr1) != '\0')
					*(++ptr) = *ptr1;
				*(++ptr) = '\0';
			}
			histptr = slowhistory[len];
			nextptr = slowhistory[len+1];
		}
		midwd++;
	}
}

boolean
slow_check_rws_reduce(w,i)
	char           		*w;
	int            		i;

/* This is similar to slow_rws_reduce, but it does not change the word  w.
 * It merely checks whether it is reduced or not, and returns true or false.
 * If the second parameter i is greater than 0, then the check is for
 * reducibility ignoring equation number i. (Usually this is used on the left-
 * hand-side of equation i itself.
 */
{
	int     st,
		**table,
                *histptr,
                *histptre,
                *nextptr,
		*slowhistorysp = rws.slowhistorysp,
		**slowhistory = rws.slowhistory;
	char	*wc = w;

restart:
	table = rws.reduction_machine->table->table_data_ptr;
        histptr=nextptr=slowhistory[0]=slowhistorysp;
	while (*wc != '\0') {
                histptre=nextptr;
		if (nextptr+(wc-w)-slowhistorysp > maxslowhistoryspace){
			if (print_level>=3)
  			  printf(
          		   "    #maxslowhistoryspace too small doubling it.\n");
		  	tfree(rws.slowhistorysp);
			maxslowhistoryspace *= 2;
			tmalloc(rws.slowhistorysp,int,maxslowhistoryspace);
			slowhistorysp = rws.slowhistorysp;
			goto restart;
		}
                while (histptr<histptre){
                        st = dense_target(table,*wc,*(histptr++));
                        if (st>0)
                                *(nextptr++) = st;
			else if (st < 0 && st != -i)
				return FALSE;
		}
                st = dense_target(table,*wc,1);
                if (st>0) 
                        *(nextptr++) = st;
                else if (st<0 && st!= -i)
                        return(0);
		wc++;
	}
	return TRUE;
}
