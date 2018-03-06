/* definition file for structure used in computing 
   power-commutator presentation and for array y */

#ifndef __PCP_VARS__
#define __PCP_VARS__

#define MAXIDENT 100 

struct pcp_vars {

   int     p;			/* prime */
   int     pm1;			/* prime - 1 */
    
   int     m;                   /* number of automorphisms */

   int     cc;			/* current class */
   int     ccbeg;		/* begin current class */
   int     clend;		/* end current class */

   int     newgen;		/* number of generators of nucleus */
   int     lastg;		/* last generator of group */
   int     first_pseudo;	/* first pseudo-generator */
   int     redgen;              /* redundant generator */

   int     fronty;		/* first storage position available in y */
   int     dgen;		/* storage location for generators */
   int     relp;		/* storage location for relations */
   int     lused;		/* last position used in y from front */
   int     gspace;		/* first garbage location available in y */
   int     words;		/* storage position of words in y */
   int     submlg;		/* position subgrp - lastg */
   int     subgrp;		/* storage position of subgroup information */
   int     structure;		/* storage position of structure information */
   int     ppower;		/* base position for power information */
   int     ppcomm;		/* base position for pointers to commutators */
   int     backy;		/* last storage place available in y */

   int     extra_relations;	/* indicate whether exponent law is imposed */
   int     start_wt;		/* start weight for exponent checking */
   int     end_wt;		/* end weight for exponent checking */

   int     ndgen;		/* number of defining generators */
   int     ndrel;		/* number of defining relations */
   int     ncomm;		/* number of commutators */
   int     nwords;		/* number of words */
   int     nsubgp;		/* number of subgroups */

   int     nocset;		/* number of occurrences parameter */
   int     complete;		/* is the group complete? */
   int     ncset;		/* is next class set up? */

   char    ident[MAXIDENT];	/* identifier of group */ 

   Logical middle_of_tails;     /* middle of tails calculation? */
   Logical update;              /* update of generators performed? */
   Logical dummy1;              /* dummy variables which can be used later */
   Logical dummy2;
   Logical dummy3;
   Logical dummy4;

   Logical metabelian;          /* is the group metabelian? */
   Logical fullop;		/* indicate nature of output */
   Logical diagn;		/* indicate nature of output */
   Logical eliminate_flag;	/* indicate that generator is eliminated */
   Logical valid;		/* indicate that input is valid */
   Logical overflow;		/* indicate integer or space overflow */
   Logical multiplicator;	/* p-multiplicator is to be computed */
   Logical cover;		/* p-covering group is to be computed */

#if defined (LIE) 
#define NRELS 3                    /* number of multilinear relations */
   int     mlin_relations[NRELS];  /* array storing degree of multilinear
                                      conditions to be imposed */
#endif

   /* variables that Magma needs */
#ifdef Magma
   t_handle group;		/* the group fed in at the top */
   t_handle y_handle;		/* handle to the work space */
   t_int output_level;		/* the actual output level */
   t_int group_type;            /* type of the group */
#endif

};

#ifndef Magma
int     *y_address;     /* definition of storage for presentation */
#endif

extern int*compact_description(Logical write_to_file, struct pcp_vars *pcp);
extern int***determine_action(int format,int *nmr_of_auts,struct pcp_vars *pcp);
extern int echelon (struct pcp_vars *pcp);
extern int find_definition ( int generator, int pointer, int weight, struct pcp_vars *pcp);
extern int layer ( int generator, struct pcp_vars *pcp);
extern int pquotient(int max_class,int output,FILE *file, int format, struct pcp_vars *pcp);
extern int pretty_filter ( FILE *file, int *max_class, int *output, struct pcp_vars *pcp);
extern int pretty_read_generators (struct pcp_vars *pcp);
extern int***read_auts_from_file(FILE*file,int*nmr_of_auts,struct pcp_vars *pcp);
extern int***read_auts(int option,int*nmr_of_auts,int*nmr_of_exponents,struct pcp_vars*pcp);
extern int vector_to_string (int cp, int str, struct pcp_vars *pcp);
extern int vector_to_word ( int cp, int ptr, struct pcp_vars *pcp);
extern void Allocate_WorkSpace (int work_space, struct pcp_vars *pcp);
extern void assemble_matrix ( int **A, int t, int** auts, struct pcp_vars *pcp);
extern void bubble_sort ( int *x, int len, struct pcp_vars *pcp);
extern void calculate_commutator (int format, struct pcp_vars *pcp);
extern void calculate_jacobi (struct pcp_vars *pcp);
extern void CAYLEY_presentation (FILE_TYPE file, struct pcp_vars *pcp);
extern void class1_eliminate (struct pcp_vars *pcp);
extern void close_queue(Logical report, int list_length, int limit, int *head, int *list, int *queue, int queue_length, struct pcp_vars *pcp);
extern void close_relations(Logical report,int limit,int queue_type,int*head,int*list,int*queue,int length,int*long_queue,int*long_queue_length,struct pcp_vars *pcp);
extern void collect_def_comm (int ptr, int cp, struct pcp_vars *pcp);
extern void collect_gen_word ( int ptr, int length, int cp, struct pcp_vars *pcp);
extern void collect (int pointer, int collected_part, struct pcp_vars *pcp);
extern void collectp2 ( int pointer, int collected_part, struct pcp_vars *pcp);
extern void collect_relations (struct pcp_vars *pcp);
extern void collect_word (int ptr, int cp, struct pcp_vars *pcp);
extern void commute_defining_generators ( int format, struct pcp_vars *pcp);
extern void compact (struct pcp_vars *pcp);
extern void consistency ( int type, int *queue, int *queue_length, int wc, struct pcp_vars *pcp);
extern void copy ( int old, int length, int new, struct pcp_vars *pcp);
extern void delete_tables(int type, struct pcp_vars *pcp);
extern void down_class ( int ptr, struct pcp_vars *pcp);
extern void eliminate (Logical middle_of_tails, struct pcp_vars *pcp);
extern void evaluate_formula (int *queue, int *queue_length, struct pcp_vars *pcp);
extern void extend_automorphisms(int ***auts, int nmr_of_auts, struct pcp_vars *pcp);
extern void Extend_Auts ( int **head, int **list, int start, struct pcp_vars *pcp);
extern void find_commutator(int cp1,int cp2,int cp3,int cp4,int result,struct pcp_vars *pcp);
extern void GAP_presentation ( FILE_TYPE file, struct pcp_vars * pcp);
extern void image_to_word(int string,int*image,struct pcp_vars *pcp);
extern void initialise_pcp (int output, struct pcp_vars *pcp);
extern void interactive_pq (Logical group_present, int format, int output_level, int **head,int **list, struct pcp_vars *pcp);
extern void invert_generator ( int gen, int exp, int cp, struct pcp_vars *pcp);
extern void invert_word ( int ptr, int cp, struct pcp_vars *pcp);
extern void isom_options ( int format, struct pcp_vars *pcp);
extern void jacobi ( int c, int b, int a, int ptr, struct pcp_vars *pcp);
extern void last_class (struct pcp_vars *pcp);
extern void List_Auts ( int *head, int *list, int first, int last, struct pcp_vars *pcp);
extern void list_commutators ( int *queue, int *queue_length, struct pcp_vars *pcp);
extern void List_Commutators ( int *queue, int *queue_length, struct pcp_vars *pcp);
extern void Magma_Auts(int *head,int *list,int start,int first, int last, struct pcp_vars *pcp);
extern void Magma_presentation ( FILE_TYPE file, struct pcp_vars *pcp);
extern void next_class ( Logical report, int **head, int **list, struct pcp_vars *pcp);
extern void options ( int call, int format, struct pcp_vars *pcp);
extern void output_information(int *sequence,int nmr_of_exponents,struct pcp_vars *pcp);
extern void pgroup_generation ( Logical *group_present, struct pcp_vars *pcp);
extern void power (int exp, int cp, struct pcp_vars *pcp);
extern void pretty_read_relations ( int output, int *max_class, struct pcp_vars *pcp);
extern void pretty_read_word(FILE *file,int disp, int type, struct pcp_vars *pcp);
extern void print_auts(int nmr_auts,int nmr_gens,int***auts,struct pcp_vars*pcp);
extern void print_level ( int *output, struct pcp_vars *pcp);
extern void print_map (struct pcp_vars *pcp);
extern void print_presentation ( Logical full, struct pcp_vars *pcp);
extern void print_structure (int first, int last, struct pcp_vars *pcp);
extern void print_word ( int ptr, struct pcp_vars *pcp);
extern void read_class_bound (int *class_bound, struct pcp_vars *pcp);
extern void read_order_bound ( int *order_bound, struct pcp_vars *pcp);
extern void read_parameters (int format,int *max_class, int *output,struct pcp_vars *pcp);
extern void read_relations (struct pcp_vars *pcp);
extern void read_relator_file ( int *queue, int *queue_length, struct pcp_vars *pcp);
extern void read_word( FILE_TYPE file, int disp, int type, struct pcp_vars *pcp);
extern void restore_automorphisms (FILE_TYPE ifp, int **head, int **list, struct pcp_vars *pcp);
extern void restore_pcp ( FILE *ifp, struct pcp_vars *pcp);
extern void save_auts(FILE_TYPE ofp,int *head,int *list,struct pcp_vars *pcp);
extern void save_pcp (FILE *ofp, struct pcp_vars *pcp);
extern void set_maxoccur (struct pcp_vars *pcp);
extern void Setup_Action(int **head,int **list,int ***auts,int nmr_of_exponents,struct pcp_vars *pcp);
extern void setup_defgen_word_to_collect(FILE_TYPE file, int format, int type, int cp, struct pcp_vars *pcp);
extern void setup_echelon ( int *queue, int *queue_length, int cp, struct pcp_vars *pcp);
extern void setup (struct pcp_vars *pcp);
extern void setup_symbols (struct pcp_vars *pcp);
extern void setup_to_solve_equation ( int format, struct pcp_vars *pcp);
extern void setup_word_to_collect(FILE_TYPE file,int format,int type,int cp, struct pcp_vars *pcp);
extern void setup_word_to_print(char *type,int cp, int str,struct pcp_vars *pcp);
extern void solve_equation ( int cp1, int cp2, int result, struct pcp_vars *pcp);
extern void string_to_vector ( int str, int cp, struct pcp_vars *pcp);
extern void tails(int type,int work_class,int start_weight, int end_weight, struct pcp_vars *pcp);
extern void trace_relation(int*sequence,int*index,int ptr,int generator,struct pcp_vars*pcp);
extern void update_generators (struct pcp_vars *pcp);
extern void update ( int ptr, struct pcp_vars *pcp);
extern void write_CAYLEY_library ( FILE_TYPE file, struct pcp_vars *pcp);
extern void write_GAP_library ( FILE_TYPE file, struct pcp_vars * pcp);
extern void write_Magma_library ( FILE_TYPE file, struct pcp_vars *pcp);
extern Logical is_genlim_exceeded(struct pcp_vars *pcp);
extern Logical is_space_exhausted(int required, struct pcp_vars *pcp);
extern void calculate_power(int exp, int ptr, int cp, struct pcp_vars *pcp);
extern void calculate_tails(int final_class,int start_weight,int end_weight,struct pcp_vars *pcp);
extern void check_input(int output, int *max_class, struct pcp_vars *pcp);
extern void collect_image_of_generator(int cp,int *auts,struct pcp_vars *pcp);
extern void collect_image_of_string(int string, int cp, int **auts,struct pcp_vars *pcp);
extern void complete_echelon(Logical trivial,int redgen,struct pcp_vars *pcp);
extern void create_tail(int address, int f, int s, struct pcp_vars *pcp);
extern void extend_automorphism(int **auts, struct pcp_vars *pcp);
extern void extend_commutator(int cp1,int cp2,int u,int v,int **auts,struct pcp_vars *pcp);
extern void extend_power(int cp1,int cp2,int u,int**auts,struct pcp_vars *pcp);
extern void extend_tail(int address, int f, int s, struct pcp_vars *pcp);
extern void invalid_group (struct pcp_vars *pcp);
extern void print_pcp_relations(struct pcp_vars *pcp);
extern void traverse_list(int exponent, int head, int *list, int cp, struct pcp_vars *pcp);
#endif
