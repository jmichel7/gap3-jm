###############################################################################
##
##  init.g  GLISSANDO ver 1.0  Christof Noebauer   1995, 1996
##
##
#############################################################################
PrintPkgInit(rec(name:="gliss",version:="1.0",date:=[1996,9,27],
  copyright:=
  "(C) Christof Noebauer -- Glissando: small semigroups and near rings"));

##  make sure to read the lib file group.g first, so that Parent can be
##  overlaid.
##
#############################################################################
##
#F  Parent( <G>, ... )  . . . . . . . . . . . . . . . . . common parent group
##
##  overlay the original Parent function because sets don't have record
##  components.  
##  JM 17-3-2016 moved to gap/lib/group.g
##
#############################################################################

#############################################################################
##  make the semigroup library files AUTO readable
#############################################################################

AUTO( ReadPkg( "gliss", "data", "sg1_4" ), SG1, SG2, SG3, SG4 );
AUTO( ReadPkg( "gliss", "data", "sg5"   ), SG5 );

#############################################################################

#############################################################################
##  make the near-ring library files AUTO readable
#############################################################################

AUTO( ReadPkg( "gliss", "data", "nr_2_7" ), 
      NR_C2, NR_C3, NR_C4, NR_V4, NR_C5, NR_C6, NR_S3, NR_C7 );

AUTO( ReadPkg( "gliss", "data", "nr8_1" ), NR_C8 ); 
AUTO( ReadPkg( "gliss", "data", "nr8_2" ), NR_C2xC4 ); 
AUTO( ReadPkg( "gliss", "data", "nr8_3" ), NR_C2xC2xC2 ); 
AUTO( ReadPkg( "gliss", "data", "nr8_4" ), NR_D8 ); 
AUTO( ReadPkg( "gliss", "data", "nr8_5" ), NR_Q8 ); 

AUTO( ReadPkg( "gliss", "data", "nr9_1" ), NR_C9 ); 
AUTO( ReadPkg( "gliss", "data", "nr9_2" ), NR_C3xC3 ); 

AUTO( ReadPkg( "gliss", "data", "nr10_1" ), NR_C10 ); 
AUTO( ReadPkg( "gliss", "data", "nr10_2" ), NR_D10 ); 

AUTO( ReadPkg( "gliss", "data", "nr11_1" ), NR_C11 ); 

AUTO( ReadPkg( "gliss", "data", "nr12_1" ), NR_C12 ); 
AUTO( ReadPkg( "gliss", "data", "nr12_2" ), NR_C2xC6 ); 

AUTO( ReadPkg( "gliss", "data", "nr12_3_0" ), NR_D12_1 ); 
AUTO( ReadPkg( "gliss", "data", "nr12_3_1" ), NR_D12_2 ); 
AUTO( ReadPkg( "gliss", "data", "nr12_3_2" ), NR_D12_3 ); 
AUTO( ReadPkg( "gliss", "data", "nr12_3_3" ), NR_D12_4 ); 
AUTO( ReadPkg( "gliss", "data", "nr12_3_4" ), NR_D12_5 ); 
AUTO( ReadPkg( "gliss", "data", "nr12_3_5" ), NR_D12_6 ); 
AUTO( ReadPkg( "gliss", "data", "nr12_3_6" ), NR_D12_7 ); 
AUTO( ReadPkg( "gliss", "data", "nr12_3_7" ), NR_D12_8 ); 
AUTO( ReadPkg( "gliss", "data", "nr12_3_8" ), NR_D12_9 ); 
AUTO( ReadPkg( "gliss", "data", "nr12_3_9" ), NR_D12_10 ); 

AUTO( ReadPkg( "gliss", "data", "nr12_4" ), NR_A4 ); 
AUTO( ReadPkg( "gliss", "data", "nr12_5" ), NR_T ); 

AUTO( ReadPkg( "gliss", "data", "nr13_1" ), NR_C13 ); 

AUTO( ReadPkg( "gliss", "data", "nr14_1" ), NR_C14 ); 
AUTO( ReadPkg( "gliss", "data", "nr14_2" ), NR_D14 ); 

AUTO( ReadPkg( "gliss", "data", "nr15_1" ), NR_C15 ); 

##############################################################################
##############################################################################

AUTO( ReadPkg( "gliss", "lib", "t" ), 
TransformationPrintLevel, IsTransformation, IsSetTransformation,
IsGroupTransformation, AsTransformation, 
IdentityTransformation_gliss, Transformation_gliss,# see gap3_jm/lib/glissmon.g
DisplayTransformation, SetTransformationOps,
GroupTransformationOps);

AUTO(ReadPkg( "gliss", "lib", "g" ),
InnerAutomorphisms, SmallestGeneratingSystem, IsIsomorphicGroup);

AUTO(ReadPkg( "gliss", "lib", "sg" ),
IsSemigroup, IsTransformationSemigroup, AllFunctions, TransformationSemigroup,
SgOps, SmallestIdeal, Green, LibrarySemigroup, AllLibrarySemigroups);
  
AUTO(ReadPkg( "gliss", "lib", "nr" ),
IsNearring, IsTransformationNearring,
IsNrMultiplication, Nearring, NearringOps, FindGroup, NearringIdeals,
InvariantSubnearrings, Subnearrings, LibraryNearring, Distributors,
DistributiveElements, ZeroSymmetricElements, IdempotentElements,
NilpotentElements, QuasiregularElements, RegularElements,
IsAbstractAffineNearring, IsDistributiveNearring, IsBooleanNearring,
IsDgNearring, IsIntegralNearring, IsNilNearring, IsNilpotentNearring,
IsPrimeNearring, IsQuasiregularNearring, IsRegularNearring,
IsNilpotentFreeNearring, CouplingFunctions, IsPlanarNearring,
IsNearfield, LibraryNearringInfo, NearringInfo);

AUTO(ReadPkg( "gliss", "lib", "smallgps" ),
C2,C3,C4,V4,C5,C6,S3,C7,C8,C2xC4,C4xC2,C2xC2xC2,Q8,D8,C9,C3xC3,C10,D10,C11,
C12,C2xC6,D12,A4,ZS,C13,C14,D14,C15);

AUTO(ReadPkg( "gliss", "lib", "dispatch" ),
DisplayCayleyTable, IsCommutative);

ReadPkg( "gliss", "lib", "patch" ); 
