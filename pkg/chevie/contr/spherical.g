###########################################################################
# Contribution to the Chevie Package
#
# spherical.g
#
# (C) July 2018 --- Daniel Juteau
#
# The  main function of  this file is  'SphericalCriterion' which allows to
# determine  the support of  the spherical module  for a rational Cherednik
# algebra  associated with a complex reflection group W, with any choice of
# parameters  (using  the  calculation  of  Schur elements based on various
# assumptions, notably the symmetry of the Hecke algebra). This is based on
# S.  Griffeth  and  D.  Juteau,  "W-exponentials,  Schur elements, and the
# support  of  the  spherical  representation  of  the  rational  Cherednik
# algebra", arXiv:1707.08196.
#
# The first task is to create a generic Hecke algebra for W. The parameters
# for  the  first  hyperplane  orbit  are  1,  E(e)  x_1,  E(e)^2 x_2, ...,
# E(e)^(e-1)  x_(e-1),  and  similarly  for  the other hyperplane orbits if
# necessary,  with x replaced by  y, z... The hyperplane  orbits are in the
# order as given by the function HyperplaneOrbits.
#
# Then,  for each maximal parabolic subgroup WI, we compute the quotient of
# the  principal Schur elements (i.e. associated to the trivial characters)
# of  W and  WI. The  locus in  the    parameter space  where the spherical
# module is not a quotient of the induction of the spherical module for the
# rational  Cherednik algebra of WI is the  zero set of that polynomial. We
# use the function FactorizedPrincipalSchurElement which was implemented by
# Maria  Chlouveraki  and  Jean  Michel  (for the determination of Rouquier
# blocks).
#
# The  function  DisplaySphericalCriterion  gives  the  result  in  a  nice
# readable  form. In particular, it  first prints the diagran  of W and the
# hyperplane   orbits;  it  also  removes  the  monomial  factor  which  is
# irrelevant.
# 
###########################################################################

MaximalParabolicRepresentatives := function ( W )
    return ParabolicRepresentatives( W, W.semisimpleRank - 1 );
end;

FactorizedPrincipalSchurElement := function(H)
  local W;
  W := H.reflectionGroup;
  return FactorizedSchurElement(H, CharParams(W)[PositionId(W)]);
end;

SphericalCriterion := function( W )
  local hyp, names, params, var, e, i, H, P, L, res, WI, HI, I;
  
  hyp := HyperplaneOrbits( W );
  names := "xyz";
  params := [];
  for i in [ 1 .. Length(hyp) ] do
    e := hyp[i].e_s;
    var := names{[i]};
    params[hyp[i].s] := Concatenation( [ 1 ], 
      List( [ 1 .. e - 1 ], j -> E(e)^j * Mvp( SPrint( var, "_", j ) ) ) );
  od;
  H := Hecke( W, params );
  params := H.parameter;
  P := FactorizedPrincipalSchurElement( H );
  
  L := MaximalParabolicRepresentatives( W );
  res := [];
  for I in L do
    WI := ReflectionSubgroup( W, I );
    HI := Hecke( WI,
      params{ W.orbitRepresentative{ WI.rootInclusion{
        [1 .. WI.nbGeneratingReflections ] } } });
    Add(res, P / FactorizedPrincipalSchurElement( HI ) );
  od;
  return res; 
end;

DisplaySphericalCriterion := function( W )
  local pols, maxpar, i, rowLabels, rowsLabel, columnLabels;
  
  maxpar := MaximalParabolicRepresentatives( W );

  pols := SphericalCriterion( W );
  for i in pols do
    i.factor := 1;
  od;
 
  rowLabels := List( maxpar, l -> SPrint( ReflectionName( ReflectionType( ReflectionSubgroup( W, l ) ) ), " ", String( l ) ) ); 
 
  rowsLabel := "Maximal parabolic subgroups";

  columnLabels := [ "q-index" ];

  Print( FormatTable( TransposedMat( [ pols ] ), rec( rowLabels := rowLabels, columnLabels := columnLabels, rowsLabel := rowsLabel ) ) ); 
end;
