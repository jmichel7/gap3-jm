/* definition file for structure used in p-group generation */

#ifndef __PGA_VARS__
#define __PGA_VARS__

struct pga_vars {
   int     p;			/* prime */

   int     ndgen;               /* rank of Frattini quotient */

   int     multiplicator_rank;  /* rank of p-multiplicator */
   int     nuclear_rank;	/* rank of nucleus */
   int     step_size;		/* step size */

   /* values of the above parameters relative to 
      chosen characteristic subgroup */
   int     q;   
   int     r;    
   int     s;    

   Logical final_stage;		/* indicates whether in intermediate stage */
   Logical capable;		/* indicates whether group is capable */

   int    fixed;		/* number of generators of the p-multiplicator
				   which cannot be eliminated */

   int     m;			/* number of automorphisms */
   int     nmr_centrals;	/* number of central automorphisms */
   int     nmr_stabilisers;	/* number of generators for stabiliser */
   int     Degree;		/* degree of permutation group */
   int    *powers;		/* store powers of prime */
   int    *inverse_modp;	/* store inverses of 0 .. p - 1 */
   int    *list;		/* list of definition sets */
   int    *available;		/* number of available positions for each set */
   int    *offset;		/* offset for each definition set */
   int    nmr_def_sets;		/* number of definition sets */
   int    nmr_subgroups;	/* number of subgroups processed */

   int    *rep;                 /* list of orbit representatives */
   int    nmr_orbits;           /* number of orbits */
   int    nmr_of_descendants;   /* number of immediate descendants */
   int    nmr_of_capables;      /* number of capable immediate descendants */

   int    *map;                 /* map from automorphisms to permutations */
   int    nmr_of_perms;         /* number of permutations */

   /* series of print flags */

   Logical print_extensions;
   Logical print_automorphism_matrix;

   Logical print_degree;
   Logical print_permutation;

   Logical print_subgroup;
   Logical print_reduced_cover;
   Logical print_group;
   Logical print_nuclear_rank;
   Logical print_multiplicator_rank;

   Logical print_orbits;
   Logical print_orbit_summary;
   Logical print_orbit_arrays;

   Logical print_commutator_matrix;
   Logical print_automorphisms;
   Logical print_automorphism_order;
   Logical print_stabiliser_array;

   Logical trace;               /* trace details of algorithm */

   /* algorithm flags */
   Logical space_efficient; 
   Logical soluble;
   Logical combined;
   Logical terminal;        /* completely process terminal descendants */
   Logical metabelian;      /* ensure descendant is metabelian */

   int exponent_law;        /* ensure descendant satisfies exponent law */

   int orbit_size;          /* total orbit size in constructing group */

   Logical dummy1;          /* dummy variables */ 
   Logical dummy2;   

   Logical upper_bound;     /* only automorphism group order upper 
                               bound stored */

#ifdef LARGE_INT 
   MP_INT aut_order;        /* order of automorphism group */
#endif 
};

extern char*find_permutation(int *b,char *c,struct pga_vars *pga);
extern int*bitstring_to_subset(int K,struct pga_vars *pga);
extern int***central_automorphisms(struct pga_vars *pga, struct pcp_vars *pcp);
extern int close_subgroup ( int k, int ***auts, struct pga_vars *pga, struct pcp_vars *pcp);
extern int**commutator_matrix(struct pga_vars *pga, struct pcp_vars *pcp);
extern int construct ( int call_depth, struct pga_vars *flag, int option,FILE_TYPE output_file, FILE_TYPE start_file, int k, int order_bound, int group_nmr,struct pga_vars *pga,  struct pcp_vars *pcp);
extern int echelonise_matrix (int **a, int nmr_rows, int nmr_columns, int p,int *subset,   struct pga_vars *pga);
extern int**find_allowable_subgroup(int option,FILE_TYPE cover_tmp_file,FILE_TYPE group_tmp_file,int*bit_string,int**subset,struct pga_vars*pga,struct pcp_vars*pcp);
extern int find_image ( int label, int **auts, struct pga_vars *pga, struct pcp_vars *pcp);
extern int*find_orbit_reps(int *a,int *b,struct pga_vars*pga);
extern int find_soluble_auts ( int ***auts, struct pga_vars *pga, struct pcp_vars *pcp);
extern int**find_stabiliser(Logical*identity_map,int non_standard,int***auts,int**perms, int*a,int*b,char*c,int*orbit_length,struct pga_vars *pga,struct pcp_vars *pcp);
extern int**finish_pga_run(Logical*identity_map,FILE_TYPE cover_tmp_file,FILE_TYPE group_file,int***auts,struct pga_vars*pga,struct pcp_vars*pcp);
extern int**group_completed(int***auts,struct pga_vars*pga,struct pcp_vars*pcp);
extern int***immediate_descendant(FILE_TYPE descendant_file,struct pga_vars *pga,struct pcp_vars *pcp);
extern int***invert_automorphisms(int***auts,struct pga_vars *pga,struct pcp_vars *pcp);
extern int**label_to_subgroup(int*Index,int**subset,int label,struct pga_vars *pga);
extern int**permute_subgroups(FILE_TYPE LINK_input,int**a,int**b,char**c,int***auts,struct pga_vars *pga,struct pcp_vars *pcp);
extern int preimage ( int perm, struct pga_vars *pga);
extern int process_identity_perm ( int *a, int *b, char *c, struct pga_vars *pga);
extern int***read_stabiliser_gens(int nmr_of_generators,int ***soluble_generators,struct pga_vars *pga);
extern int ***read_stabiliser_gens(int nr_gens, int *** sol_gens, struct pga_vars * pga);
extern int reduced_covers ( FILE_TYPE descendant_file, FILE_TYPE covers_file, int k,int ***auts,  struct pga_vars *pga, struct pcp_vars *pcp);
extern int***restore_group(Logical rewind_flag,FILE_TYPE input_file,int group_number,struct pga_vars *pga,struct pcp_vars *pcp);
extern int***restore_pga(FILE *ifp,struct pga_vars *pga,struct pcp_vars *pcp);
extern int***setup_identity_auts(int nmr_of_generators,int***auts,struct pga_vars*pga);
extern int***stabiliser_of_rep(int**perms,int rep,int orbit_length,int*a,int*b,char*c,char*d,int***auts,struct pga_vars *pga,struct pcp_vars *pcp);
extern int**start_pga_run(Logical*identity_map,int***auts,struct pga_vars*pga,struct pcp_vars*pcp);
extern int subgroup_to_label (int **S, int K, int *subset, struct pga_vars *pga);
extern void autgp_order(struct pga_vars *pga, struct pcp_vars *pcp);
extern void combined_computation ( int ***auts, int **a, int **b, char **c, int **perms,int **orbit_length, struct pga_vars *pga, struct pcp_vars *pcp);
extern void compute_degree (struct pga_vars *pga);
extern void compute_images(int**A,int K,int depth,int *permutation,struct pga_vars *pga);
extern void compute_orbits (int **a, int **b, char **c, int **perms, struct pga_vars *pga);
extern void compute_permutation(int*permutation,int**A,struct pga_vars *pga);
extern void copy_flags ( struct pga_vars *flag, struct pga_vars *pga);
extern void defaults_pga ( int option, int *k, struct pga_vars *flag, struct pga_vars *pga,struct pcp_vars *pcp);
extern void evaluate_generators(int pointer,int nmr_of_generators,int*** stabiliser,int ***auts,struct pga_vars *pga, struct pcp_vars *pcp);
extern void factorise_subgroup ( int **S, int index, int *subset, struct pga_vars *pga, struct pcp_vars *pcp);
extern void find_available_positions(int K,int**A,int**Image,int**row,int**column,int depth,struct pga_vars *pga);
extern void find_padic(int x,int k,int p,int *expand,struct pga_vars *pga);
extern void free_space(Logical soluble_computation, int **perms, int *orbit_length, int *a, int *b, char *c, struct pga_vars *pga);
extern void GAP_auts(FILE_TYPE file,int***central,int***stabiliser,struct pga_vars*pga,struct pcp_vars * pcp);
extern void image_of_generator(int pointer,int generator,int***auts,struct pga_vars *pga,struct pcp_vars *pcp);
extern void initialise_pga (struct pga_vars *pga, struct pcp_vars *pcp);
extern void insoluble_compute_orbits(int **orbit,int **backptr, char **schreier, int **perms,struct pga_vars *pga);
extern void interactive_pga ( Logical group_present, FILE_TYPE StartFile, int group_nmr, int ***auts, struct pga_vars *pga, struct pcp_vars *pcp);
extern void iteration_information ( int *subgroup_rank, struct pga_vars *flag, int*class_bound, int *order_bound, int **step_sequence, struct pga_vars *pga,struct pcp_vars *pcp);
extern void iteration(int call_depth,int*step_sequence,int subgroup_rank,struct pga_vars*flag,FILE_TYPE input_file,int nmr_of_descendants,int class_bound,int order_bound,struct pga_vars*pga,struct pcp_vars *pcp);
extern void map_relations ( int **map, struct pga_vars *pga, struct pcp_vars *pcp);
extern void orbit_summary ( int *length, struct pga_vars *pga);
extern void print_aut_description ( int ***central, int ***stabiliser, struct pga_vars *pga,struct pcp_vars *pcp);
extern void print_orbit_information ( int *a, int *b, char *c, struct pga_vars *pga);
extern void reduce_matrix(int **a,int nmr_rows,int nmr_columns,int p,struct pga_vars *pga);
extern void report_autgp_order ( struct pga_vars *pga, struct pcp_vars *pcp);
extern void save_pga ( FILE *ofp, int ***central, int ***stabiliser, struct pga_vars *pga,struct pcp_vars *pcp);
extern void set_defaults (struct pga_vars *flag);
extern void setup_identity_perm (int *permutation, struct pga_vars *pga);
extern void setup_reps(int*reps,int nmr_of_reps,int*orbit_length,int**perms,int*a,int*b,char*c,int ***auts, FILE_TYPE descendant_file,FILE_TYPE covers_file, struct pga_vars *pga, struct pcp_vars *pcp);
extern void set_values (struct pga_vars *pga, struct pcp_vars *pcp);
extern void stabiliser_generators(int**perms,int rep,int*a,int*b,char*c,char*d,int***auts,struct pga_vars *pga, struct pcp_vars *pcp);
extern void standard_presentation ( Logical *identity_map, int standard_output, int ***auts,struct pga_vars *pga, struct pcp_vars *pcp);
extern void start_GAP_file ( int  *** auts, struct pga_vars   * pga);
extern void StartGapFile (struct pga_vars   * pga);
extern void start_group(FILE_TYPE*StartFile,int***auts,struct pga_vars*pga,struct pcp_vars*pcp);
extern void start_stage ( FILE_TYPE descendant_file, int k, int ***auts, struct pga_vars *pga, struct pcp_vars *pcp);
extern void step_range ( int k, int *lower_step, int *upper_step, int ***auts,struct pga_vars *pga, struct pcp_vars *pcp);
extern void store_definition_sets (int r, int lower_step, int upper_step, struct pga_vars *pga);
extern void strip_identities (int ***auts, struct pga_vars *pga, struct pcp_vars *pcp);
extern void update_autgp_order (int orbit_length, struct pga_vars *pga, struct pcp_vars *pcp);
extern void update_image(int**A,int column,int**Image,int row,struct pga_vars *pga);
extern void add_to_list (int *subset, struct pga_vars *pga);
extern void enforce_laws ( struct pga_vars *flag, struct pga_vars *pga, struct pcp_vars *pcp);
extern void get_definition_sets(struct pga_vars *pga);
extern void orbit_option(int option,int**perms,int**a,int**b,char**c,int**orbit_length,struct pga_vars *pga);
extern void orbits (int *permutation,int *a,int *b,char *c,struct pga_vars *pga);
extern void print_group_details(struct pga_vars *pga, struct pcp_vars *pcp);
extern void query_aut_group_information (struct pga_vars *pga);
extern void query_degree_aut_information (struct pga_vars *pga);
extern void query_exponent_law (struct pga_vars *pga);
extern void query_group_information(int p, struct pga_vars *pga);
extern void query_metabelian_law (struct pga_vars *pga);
extern void query_orbit_information (struct pga_vars *pga);
extern void query_perm_information (struct pga_vars *pga);
extern void query_solubility (struct pga_vars *pga);
extern void query_space_efficiency (struct pga_vars *pga);
extern void query_terminal (struct pga_vars *pga);
extern void read_step_size (struct pga_vars *pga, struct pcp_vars *pcp);
extern void report(int nmr_of_capables, int nmr_of_descendants, int nmr_of_covers, struct pga_vars *pga, struct pcp_vars *pcp);
extern void space_for_orbits(int **a,int **b,char **c,struct pga_vars *pga);
extern void stabiliser_option(int option,int***auts,int **perms,int *a, int *b, char *c, int *orbit_length, struct pga_vars *pga, struct pcp_vars *pcp);
#endif
