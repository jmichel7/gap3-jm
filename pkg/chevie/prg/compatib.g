##########################################################################
# Compability section: obsolete/superceded functions are defined here
# compatibility functions used inside other compatibility functions are first

AlgebraicFundamentalGroup:=function(W)
  Print("#W AlgebraicFundamentalGroup is obsolete; use FundamentalGroup\n");
  return FundamentalGroup(W);
end;

BraidWords:=function(W,l)
  Print("#W BraidWords is obsolete; use GarsideWords(BraidMonoid(W),l)\n");
  return GarsideWords(BraidMonoid(W),l);
end;

CartanName:=function(arg)
  Print("#W CartanName is obsolete; use ReflectionName\n");
  return ApplyFunc(ReflectionName,arg);
end;

CartanType:=function(arg)
  Print("#W CartanType is obsolete; use ReflectionType\n");
  return ApplyFunc(ReflectionType,arg);
end;

CharHeckeRepresentation:=function(arg)
  Print("#W CharHeckeRepresentation is obsolete; use CharRepresentationWords\n");
  return ApplyFunc(CharRepresentationWords,arg);
end;

CoxeterConjugacyClasses:=function(W)
  Print("#W CoxeterConjugacyClasses is obsolete; use ChevieClassInfo(W).classtext\n");
  return ChevieClassInfo(W).classtext;
end;

CoxeterElementsLength:=function(arg)
  Print("#W CoxeterElementsLength is obsolete; use CoxeterElements\n");
  return ApplyFunc(CoxeterElements,arg);
end;

WeylReflections:=W->List(Reflections(W),x->CoxeterWord(W,x));

CoxeterWordReflections:=WeylReflections;

DirectSumCartanMat:=function(arg)
  Print("#W DirectSumCartanMat is obsolete; use DiagonalMat\n");
  return ApplyFunc(DiagonalMat,arg);
end;

DirectSumMat:=function(arg)
  Print("#W DirectSumMat is obsolete; use DiagonalMat\n");
  return ApplyFunc(DiagonalMat,arg);
end;

DoublePartitions:=function(n)
  Print("#W DoublePartitions(n) is obsolete; use PartitionTuples(n,2)\n");
  return PartitionTuples(n,2);
end;

DoublePartitionToString:=function(arg)
  Print("#W DoublePartitionToString is obsolete; use PartitionTupleToString\n");
  return ApplyFunc(PartitionTupleToString,arg);
end;

HeckeCharTable:=function(arg)
 Error("#W HeckeCharTable is obsolete, use CharTable(Hecke(CoxeterGroup(...)))\n");
end;

HeckeFusion:=function(arg)
  Print("#W is obsolete, use Specialization and InductionTable\n");
end;

HeckePowermap:=function(arg)
 Error("#W is obsolete, use Specialization and Character tables\n");
end;

HeckeCharTableDirectProduct:=function(arg)
 Error("#W is obsolete, use CharTable(Hecke(CoxeterGroup(...)))\n");
end;

HeckeScalarProducts:=function(arg)
  Error("#W is obsolete, use Specialization and CharTable(Hecke(CoxeterGroup(...)))\n");
end;

HeckeReducedChar:=function(arg)
  Error("#W is obsolete, use Specialization and CharTable(Hecke(CoxeterGroup(...)))\n");
end;

KLCoefficient:=function(arg)
  Print("#W KLCoefficient is obsolete; use KazhdanLusztigCoefficient\n");
  return ApplyFunc(KazhdanLusztigCoefficient,arg);
end;

LeftCellRepresentation:=function(arg)
  Print("#W LeftCellRepresentation(H,c) is obsolete; use Representation(c,H)\n");
  return Representation(arg[2],arg[1]);
end;

LongestWeylWord:=function(arg)
  Print("#W LongestWeylWord is obsolete; use LongestCoxeterWord\n");
  return ApplyFunc(LongestCoxeterWord,arg);
end;

PermBraid:=function(arg)
  Print("#W PermBraid is obsolete; use EltBraid\n");
  return ApplyFunc(EltBraid,arg);
end;

PermCoxeterWord:=function(arg)
  Print("#W PermCoxeterWord is obsolete; use EltWord\n");
  return ApplyFunc(EltWord,arg);
end;

PermRepresentationRoots:=function(mats,roots)
  return List(roots*mats,SortingPerm)/SortingPerm(roots);
end;

PermWeylWord:=function(arg)
  Print("#W PermWeylWord is obsolete; use PermCoxeterWord\n");
  return ApplyFunc(PermCoxeterWord,arg);
end;

PhiFactors:=function(arg)
  Print("#W PhiFactors is obsolete; use List(ReflectionDegrees(WF),x->x[2])\n");
  return List(ReflectionDegrees(arg[1]),x->x[2]);
end;

PositionSgn:=function(arg)
  Print("#W PositionSgn is obsolete; use PositionDet\n");
  return ApplyFunc(PositionDet,arg);
end;

PrintDynkinDiagram:=function(arg)
  Print("#W PrintDynkinDiagram is obsolete; use PrintDiagram\n");
  ApplyFunc(PrintDiagram,arg);
end;

ReducedInCoxeterCoset:=function(arg)
  Print("#W ReducedInCoxeterCoset is obsolete; use ReducedInRightCoset\n");
  return ApplyFunc(ReducedInRightCoset,arg);
end;

ReducedWeylWord:=function(arg)
  Print("#W ReducedWeylWord is obsolete; use ReducedCoxeterWord\n");
  return ApplyFunc(ReducedCoxeterWord,arg);
end;

Rootsystem:=function(mat)
  Print("#W RootSystem(C) is obsolete; use CoxeterGroup(C).roots\n");
  return CoxeterGroup(mat).roots;
end;

SemisimpleCentralizer:=function(arg)
  Print("#W SemisimpleCentralizer is obsolete; use Centralizer\n");
  return ApplyFunc(Centralizer,arg);
end;

SimpleLeftDivisors:=function(arg)
  Print("#W SimpleLeftDivisors is obsolete; use LeftDivisorsSimple\n");
  return ApplyFunc(LeftDivisorsSimple,arg);
end;

SimpleReflexionMatrices:=function(mat)
  Print("#W SimpleReflectionMatrices(C) is obsolete; use CoxeterGroup(C).matgens\n");
  return CoxeterGroup(mat).matgens;
end;

Weyl:=function(arg)
  Print("#W Weyl is obsolete; use CoxeterGroup\n");
  return ApplyFunc(CoxeterGroup,arg);
end;

WeylClassPolynomials:=function(arg)
  Error("#W WeylClassPolynomials is obsolete, use 'HeckeClassPolynomials'\n");
end;

WeylConjugacyClasses:=function(W)
  Print("#W WeylConjugacyClasses is obsolete; use ChevieClassInfo(W).classtext\n");
  return ChevieClassInfo(W).classtext;
end;

WeylCosetPermRepresentation:=function(W,I)
  return PermCosetsSubgroup(W,ReflectionSubgroup(W,W.rootInclusion{I})); 
end;

WeylElements:=function(arg)
  Print("#W WeylElements is obsolete; use CoxeterWords\n");
  return ApplyFunc(CoxeterWords,arg);
end;

WeylLengthPerm:=function(arg)
  Print("#W WeylLengthPerm is obsolete; use CoxeterLength\n");
  return ApplyFunc(CoxeterLength,arg);
end;

WeylMueMat:=function(arg)
  Print("#W WeylMueMat is obsolete; use KLMueMat\n");
  return ApplyFunc(KLMueMat,arg);
end;

WeylWordPerm:=function(arg)
  Print("#W WeylWordPerm is obsolete; use CoxeterWord\n");
  return ApplyFunc(CoxeterWord,arg);
end;

WeylRightCosetRepresentatives:=function(W,I,J)
    return List(ReducedRightCosetRepresentatives(
                   ReflectionSubgroup(W,W.rootInclusion{I}),
                   ReflectionSubgroup(W,W.rootInclusion{J})),
                   x->CoxeterWord(W,x));
end;

WordBraid:=function(arg)
  Print("#W WordBraid is obsolete; use AsWord\n");
  return ApplyFunc(AsWord,arg);
end;

#########################################################################
#F  ParametersCentralizers(  <W> )  .  .  parameters for  the  conjugacy
##  classes of  centralizers of semisimple elements in a Chevalley group.
##
##  The  function  returns  a  list  of  pairs  which  parameterize  the
##  conjugacy classes.  The pairs are representatives  under W-conjugacy
##  of pairs of a subset I of the extended basis of the root system of W
##  (the fundamental roots union the  opposite of the highest root), and
##  an element of the stabilizer in W of I.
##
##  Such a pair [I,w] represents the twisted reductive subgroup
##  Spets(ReflectionSubgroup(W,I),EltWord(W,w))
##
ParametersCentralizers := function(W)
  Print("#W ParametersCentralizers is obsolete; use\n",
        "   SemisimpleCentralizerRepresentatives and Twistings\n");
  return Concatenation(List(SemisimpleCentralizerRepresentatives(W,0),
     I->List(Twistings(W,I),function(x)local W;W:=Group(x);
    return [W.rootInclusion{W.generatingReflections},
      CoxeterWord(Parent(W),x.phi)];end)));
end;
