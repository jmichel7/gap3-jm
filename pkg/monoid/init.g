#############################################################################
##
#A  init.g  
##
##  This is the init file of the MONOID package.
##
PrintPkgInit(rec(name:="monoid",version:="2.2"));

AUTO(  ReadPkg( "monoid", "lib", "action" ),
  ShortOrbit,
  GradedOrbit,
  StrongOrbit,
  StrongOrbits,
  Action,
  ActionWithZero );

AUTO( ReadPkg( "monoid", "lib", "monoid" ),
  MonoidElementsOps,
  SemiGroupOps,
  MonoidOps,
  RClassOps,
  LClassOps,
  HClassOps,
  DClassOps,
  IsMonoidElement,
  SemiGroup,
  IsSemiGroup,
  Monoid,
  IsMonoid,
  RClass,
  IsRClass,
  OnRClassesAntiAction,
  LClass,
  IsLClass,
  OnLClasses,
  HClass,
  IsHClass,
  DClass,
  IsDClass,
  DClasses,
  LClasses,
  RClasses,
  HClasses,
  SchutzenbergerGroup);

AUTO( ReadPkg( "monoid", "lib", "monorela" ), RelMonoidOps );

AUTO( ReadPkg( "monoid", "lib", "monotran" ),
  TransSemiGroupOps,
  TransMonoidOps,
  TransRClassOps,
  TransLClassOps,
  TransHClassOps,
  TransDClassOps,
  FullTransMonoidOps,
  PartialTransMonoidOps,
  IsTransMonoid,
  IsRegularTrans,
  OnKernelsAntiAction,
  OnTuplesOfSetsAntiAction,
  IsRegularDClass,
  KernelsTransMonoid,
  ImagesTransMonoid,
  FullTransMonoid,
  PartialTransMonoid );

AUTO( ReadPkg( "monoid", "lib", "relation" ),
  RelationsOps,
  RelationOps,
  Relation,
  IsRelation,
  RandomRelation,
  RelTrans,
  TransRel,
  IdentityRelation,
  EmptyRelation,
  InverseRelation,
  InverseTransformation,
  IsReflexive,
  ReflexiveClosure,
  IsSymmetric,
  SymmetricClosure,
  IsTransitiveRel,
  IsAntisymmetric,
  IsPreOrder,
  IsPartialOrder,
  IsEquivalence,
  EquivalenceClasses,
  HasseDiagram );

AUTO( ReadPkg( "monoid", "lib", "transfor" ),
  TransformationsOps,
  TransformationOps,
  IsTransformation,
  TransformationNC,
  IdentityTransformation_mon, Transformation_mon,# see gap3_jm/lib/glissmon.g
  TransPerm,
  PermTrans,
  PermLeftQuoTrans,
  MappingTransformation );
