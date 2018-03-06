###########################################################################
##
#A  initial.g                autag package                 Michael J Smith
##
##  November 1996
##
##  This file forms part of a package for computing decompositions of
##  modules into indecomposable summands, as well as computing generating
##  sets for module automorphisms.
##
##  It is an integral part of the soluble group automorphism package.
##
###########################################################################

GModOps := rec();
GModOps.GrpAlg := rec();

GModRecOps := rec();

InfoDecompose := Ignore;
InfoHom := Ignore;
InfoStack := Ignore;
