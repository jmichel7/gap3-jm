###########################################################################
##
#A  initial.g                autag package                 Michael J Smith
##
##  November 1996
##
##  This file forms part of a package for computing automorphism groups of
##  finite soluble groups which are given in terms of special soluble group
##  presentations.
##
###########################################################################

# These are for debugging or tracing execution (or for the curious)
InfoOrbit := Ignore;
InfoOrbit2 := Ignore;
InfoAutgroup := Ignore;

# This is for helpful information
InfoAut := Print; 

AutGroupOps := OperationsRecord("AutGroupOps",GroupOps);
AutOps := OperationsRecord("AutOps",GroupElementOps);
AutGroupElements := rec();

OrbitOps := rec();
PairOps := rec();
