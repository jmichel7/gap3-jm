/* 	$Id: siscode.h,v 1.2 1995/08/10 11:51:29 pluto Exp $	 */
/* 	$Log: siscode.h,v $
 * 	Revision 1.2  1995/08/10 11:51:29  pluto
 * 	Added code to handle gmodule and cohomology structures.
 *	 */

void node_to_code                 _(( node p ));
node code_to_node                 _(( void ));
void pcgroup_to_code              _(( PCGRPDESC *g, int bracket ));
PCGRPDESC *code_to_pcgroup        _(( DYNLIST p ));
void hom_to_code                  _(( HOM *a ));
HOM *code_to_hom                  _(( DYNLIST p ));
GRHOM *code_to_grhom              _(( DYNLIST p ));
void grhom_to_code                _(( GRHOM *a ));
void gmodule_to_code              _(( GMODULE *gm, int bracket ));
GMODULE *code_to_gmodule          _(( DYNLIST p ));
void cohomology_to_code           _(( COHOMOLOGY *cohomol ));
COHOMOLOGY *code_to_cohomology    _(( DYNLIST p ));








