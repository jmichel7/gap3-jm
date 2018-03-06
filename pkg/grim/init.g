PrintPkgInit(rec(name:="grim",version:="1.0",copyright:=
  "(C) beals@math.arizona.edu -- Groups of Rational and Integer Matrices"));

# read grim's init file
AUTO(ReadPkg( "grim", "lib", "grim" ),
CycPolsWithSmallDegrees, HermiteNormalForm_grim, IntegerMatrixGroupIsFinite, 
InvariantLattice, IsFiniteDeterministic, MinPolyMonteCarlo, OrderRatMat, 
PositiveDefinite, RandomIntegerComb, SimplestRational, SymMatToVec, 
SymmetricSquareMat, VecToSymMat);

# the following forces loading lib/matgrp.g but this is unavoidable
MatGroupOps.IsFinite := function ( G )
  if IsFinite( G.field )  then return true;
  elif Field(Flat(G.generators)) = Rationals then
      if G.generators=[] or 
	 Maximum(List(Flat(G.generators), Denominator)) = 1 then
	  return IntegerMatrixGroupIsFinite(G);
      elif InvariantLattice(G) <> false and
	   IntegerMatrixGroupIsFinite(G.integerMatrixGroup) then
	  G.isFinite := true;
	  G.quadraticForm := TransposedMat(G.invariantLattice)*
			     G.integerMatrixGroup.quadraticForm*
			     G.invariantLattice;
	  return true;
      else
	  G.isFinite := false;
	  return false;
      fi;
  else
      return GroupOps.IsFinite( G );
  fi;
end;
