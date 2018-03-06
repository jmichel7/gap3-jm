/****************************************************************************
**
**                  Memory Allocation Header File
**
*/

extern void Error();

extern  GroupWord NewGroupWord( /* <nr> */ );
extern  GroupWord GroupWordOne();
extern  GroupWord CopyGroupWord( /* <w> */ );
extern  RingWord NewRingWord( /* <nr> */ );
extern  RingWord RingWordOne();
extern  RingWord CopyRingWord( );
extern  Vector NewVector( /* <nr> */ );
extern  Vector VectorOne();
extern  Vector NormalNewVector( /* <nr> */ );
extern  Vector CopyVector( /* <nr> */ );
extern  ExtensionElement *NewExtensionElement();
extern  ExtensionElement *CopyExtensionElement(/*<e>*/);
extern int PrintMERelation(ExtensionElement *lhs, ExtensionElement *rhs, FILE
*f );
extern  Presentation *NewPresentation();
extern  void PrintExPresentation (FILE *f);
extern  void PrintGroupWord (GroupWord g,FILE *f,int ve);
extern  void PrintGeneratorVE( FILE *fpt, int i );
extern  void PrintGenerator( FILE *fpt, int i );
extern  void PrintExtensionElement (ExtensionElement *e, FILE *f);

extern  void FreeGroupGenerator( /* <g> */ );
extern  void FreeGroupWord( /* <w> */ );
extern  void FreeRingWord( /* <w> */ );
extern  void FreeModuleGenerator( /* <m> */ );
extern  void FreeVector( /* <v> */ );
extern  void FreeExtensionElement( /* <e> */ );
