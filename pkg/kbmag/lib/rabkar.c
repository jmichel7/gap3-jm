/* file rabkar.c   27.3.95.
 * This file contains all routines related to Rabin-Karp data structures for
 * word reduction under kbprog.
 * The word-reduction routines are modified versions of those in rwsreduce.
 * The ones here use both the reduction machine and the Rabin-Karp hash-table.
 */
#include "defs.h"
#include "fsa.h"
#include "rws.h"
#include "externals.h"

extern int      maxslowhistoryspace,
                maxreducelen;
extern rewriting_system rws;

static unsigned
prime_constants[] = {20011,60013,120011,200003,300007,400009,500009,700001,
		      1000003,1500007,2000003,5000011,10000019,32000011};
static int num_prime_constants = 14;
/* For use as modulus in hash-tables */
static unsigned bighashmod = 32000011; /* the big modulus chosen */
static unsigned hashmod; /* the modulus chosen */
static int shiftbits; /* 2^shiftbits is factor in hash-function */
static unsigned **plusinc, *minusinc;
                      /* Lookoup table of increments for hash-function */
static int *first_tail, *next_tail;
static unsigned *bighashvaleqn;
  /* All strings have a large hashvalue (taken mod bighashmod) and an
   * ordinary hash value (same thing taken mod hashmod).
   * first_tail[i] for 0<=i < hashmod, is the first equation number whose
   * tail of lhs has hash-value i. *next_tail[i] is the number of an equation
   * with the same hash-value as eqn number i - or -1 if no more.
   * bighashvaleqn holds the big hash value of equation number i.
   * Before we actually compare strings,we check that the big hash values
   * correspond.
   */

/* Functions defined in this file: */
void rk_init();
void rk_reset();
void rk_clear();
void rk_add_lhs();
void slow_rws_reduce_rk();
boolean slow_check_rws_reduce_rk();

void
rk_init(maxeqns)
	int maxeqns;
/* Initialize hash-table and lookup tables for this value of maxeqns */
{ int ng, i, j, tp;
  /* First select a suitable prime */
  for (i=0;i<num_prime_constants;i++) if (prime_constants[i] > (3*maxeqns/5)) {
    hashmod = prime_constants[i];
    break;
  }
  if (i==num_prime_constants) {
    fprintf(stderr,"Value of maxeqns is absurdly large! Cannot cope.\n");
    exit(1);
  }

  /* Now calculate shiftbits */
  ng = rws.num_gens;
  tp = 2; shiftbits=1;
  while (ng > tp) {
    tp*=2; shiftbits++;
  }
  if (print_level >= 3)
    printf(
 "    #Initializing Rabin-Karp tables. Maxeqns=%d; Modulus=%d; Shiftbits=%d.\n",
     maxeqns,hashmod,shiftbits);

  /* Now initialize lookup tables */
  tmalloc(plusinc,unsigned *,rws.rkminlen);
  tmalloc(minusinc,unsigned,ng+1);
  for (i=0;i<rws.rkminlen;i++) {
    tmalloc(plusinc[i],unsigned,ng+1);
    if (i==0) for (j=1;j<=ng;j++)
      plusinc[i][j] = j-1; /* j-1 to make values in range [0,ng-1] */
    else for (j=1;j<=ng;j++)
      plusinc[i][j] = (plusinc[i-1][j]<<shiftbits) % bighashmod;
    if (i==rws.rkminlen-1) for (j=1;j<=ng;j++)
      minusinc[j] = plusinc[i][j] ? bighashmod-plusinc[i][j] : 0;
  }

  /* Initialize hash-table */
  tmalloc(first_tail,int,hashmod);
  tmalloc(next_tail,int,maxeqns+1);
  tmalloc(bighashvaleqn,unsigned,maxeqns+1);
}

void
rk_reset(maxeqns)
	int maxeqns;
{ int i;
  for (i=0;i<hashmod;i++)
    first_tail[i] = 0;
}


void
rk_clear()
{ int i;
  if (print_level>=3)
    printf("    #Calling rk_clear.\n");
  for (i=0;i<rws.rkminlen;i++)
    tfree(plusinc[i]);
  tfree(plusinc);
  tfree(minusinc);
  tfree(first_tail);
  tfree(next_tail);
  tfree(bighashvaleqn);
}
    
void
rk_add_lhs(n)
	int n;
/* Calculate the hash-value of the tail of the lhs of equation number n
 * (it is assumed that the length of this lhs is at least rws.rkminlen),
 * and adjoin it to the hash-table.
 */
{ char *w, *we;
  int i, eqno, neweqno;
  unsigned bighashval, hashval;
  w = rws.eqns[n].lhs;
  we = w+strlen(w);
  bighashval = 0;
  for (i=0;i<rws.rkminlen;i++)
    bighashval += plusinc[i][*(--we)];
  bighashval %= bighashmod;
  bighashvaleqn[n] = bighashval;
  hashval = bighashval % hashmod;
  if ((eqno=first_tail[hashval]) == 0)
    first_tail[hashval]=n;
  else {
    while (neweqno=next_tail[eqno]) 
      eqno=neweqno;
    next_tail[eqno] = n;
  }
  next_tail[n]=0;
}

void
slow_rws_reduce_rk(w)
	char     	      	*w;
/* This is the routine slow_rws_reduce (from rwsreduce.c) adapted to deal
 * in addition with reductions arising from the Rabin-Karp table.
 */
{
	int             len,
	                st,
			j,
			l,
			rkminlen = rws.rkminlen,
			**table,
		       *histptr,
		       *histptre,
		       *nextptr,
			subrelno,
	                longer,
			*slowhistorysp = rws.slowhistorysp,
			**slowhistory = rws.slowhistory;
	unsigned	bighashval,
			hashval;
	char           *midwd,
	               *ptr1,
	               *ptr2,
	               *ptr,
		       *lhs;

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
		if (subrelno==0 && len>=rkminlen) {
	/* Look for a substitution using Rabin-Karp table  -
	 * First calculate the hash-value for the word ending at this point
	 */
			if (len==rkminlen) {
				bighashval = 0;
				ptr1=midwd;
				for (j=0;j<rkminlen;j++)
					bighashval += plusinc[j][*(ptr1--)];
			}
			else {
				bighashval += minusinc[*(midwd-rkminlen)];
				bighashval = bighashval<<shiftbits;
				bighashval += (*midwd - 1);
			}
			bighashval %= bighashmod;
			hashval = bighashval % hashmod;
			subrelno=first_tail[hashval];
			while (subrelno) {
			/* check if we have read lhs of equation no subrelno */
				lhs = rws.eqns[subrelno].lhs;
				l = strlen(lhs);
				if (l <= len &&
					  bighashval==bighashvaleqn[subrelno]) {
					ptr1 = lhs;
					ptr2 = midwd-l+1;
					ptr = ptr1+l;
					while (ptr1<ptr) {
						if (*(ptr1) != *(ptr2))
							break;
						ptr1++; ptr2++;
					}
					if (ptr1==ptr) {
						break;
					}
                                }
				subrelno = next_tail[subrelno];
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
			if (len>=rkminlen) {
				bighashval = 0;
				ptr1=midwd;
				for (j=0;j<rkminlen;j++)
					bighashval += plusinc[j][*(ptr1--)];
				bighashval %= bighashmod;
			}
		}
		midwd++;
	}
}

boolean
slow_check_rws_reduce_rk(w,i)
	char     	      	*w;
	int			i;
/* This is the routine slow_check_rws_reduce (from rwsreduce.c) adapted to deal
 * in addition with reductions arising from the Rabin-Karp table.
 */
{
	int             len,
	                st,
			j,
			l,
			rkminlen = rws.rkminlen,
			**table,
			subrelno,
		       *histptr,
		       *histptre,
		       *nextptr,
	                longer,
			*slowhistorysp = rws.slowhistorysp,
			**slowhistory = rws.slowhistory;
	unsigned	bighashval,
			hashval;
	char           *midwd,
	               *ptr1,
	               *ptr2,
	               *ptr,
		       *lhs;

restart:
	midwd = w;
	len = 0;
	table = rws.reduction_machine->table->table_data_ptr;
	nextptr=slowhistory[0]=slowhistorysp;
	while (*midwd != '\0') {
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
		if (st<0 && st != -i)
			return FALSE;
		else{
			if (st>0)
				*(nextptr++) = st;
			while (histptr<histptre){
				st = dense_target(table,*midwd,*(histptr++));
				if (st>0)
					*(nextptr++) = st;
				else if (st<0 && st!= -i)
					return FALSE;
			}
		}
		if (len>=rkminlen) {
	/* Look for a substitution using Rabin-Karp table  -
	 * First calculate the hash-value for the word ending at this point
	 */
			if (len==rkminlen) {
				bighashval = 0;
				ptr1=midwd;
				for (j=0;j<rkminlen;j++)
					bighashval += plusinc[j][*(ptr1--)];
			}
			else {
				bighashval += minusinc[*(midwd-rkminlen)];
				bighashval = bighashval<<shiftbits;
				bighashval += (*midwd - 1);
			}
			bighashval %= bighashmod;
			hashval = bighashval % hashmod;
			subrelno=first_tail[hashval];
			while (subrelno) {
			/* check if we have read lhs of equation no subrelno */
				lhs = rws.eqns[subrelno].lhs;
				l = strlen(lhs);
				if (l<=len && subrelno!=i &&
					  bighashval==bighashvaleqn[subrelno]) {
					ptr1 = lhs;
					ptr2 = midwd-l+1;
					ptr =  ptr1+l;
					while (ptr1<ptr) {
						if (*(ptr1) != *(ptr2))
							break;
						ptr1++; ptr2++;
					}
					if (ptr1==ptr)
						break;
                                }
				subrelno = next_tail[subrelno];
                        }
			if (subrelno>0)
				return FALSE;
                }
		midwd++;
	}
	return TRUE;
}
