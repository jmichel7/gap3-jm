AUTO(ReadChv("work/ax"),fitfam,cmpeq,cmpsgn,showsgn,displayfam,cmpfam);
AUTO(ReadChv("work/det"),myDet,myDet2,myCoF);
AUTO(ReadChv("work/ennola"),OmegaChi,possxi,Ennola,CheckEnnolaRLG,detfam,
  detfrob,EnnolaTwist,EigenEnnola);
AUTO(ReadChv("work/getunpdeg"),EigenAndDegHecke);
AUTO(ReadChv("work/hecke"),gete);
AUTO(ReadChv("work/lib"),dr,dp,braidmatrix,GaloisAction,OrderCenter,
  DistinctCartesian,DecomposeChar,cmpcycpol,
  2I,2B2,2G2,2F4,2A,2D,2E6,3D4,3G333,3pG333,4G333,2G5,permcoeff,
  ProperStableParabolics,ProperParabolics,TeX,TeXs,LaTeX,LaTeXs,
  cmpvec,cmpobj,cmptbl);
AUTO(ReadChv("work/linsolve"),LinSysOps,LinSys,apply,checkMagrees,applynew,
  Makerel,Values,MvpSolve,MvpLinForm);
AUTO(ReadChv("work/mvptools"),MvpValuation,mRatFrac,mMvp,solvefor,
  shr,shrint,clean,weedfacts);
AUTO(ReadChv("work/ps"),PositionsSgn,LsToPs,PsToLs,LsPuiss,MatToLs,
  PsPermuted,DisplayPs,DisplayLs,PsMatStab,PsMatMat,LsToPermAndSigns);
AUTO(ReadChv("work/series"),AllProperSeries,FactorsSet,CheckSeries,Series,
  SeriesOps,CuspidalPairs,CuspidalSeries,EigenspaceNumbers, FitParameter,
  getHecke, IsSeries, PrincipalSeries,
  PermutationOnClasses,PermutationOnCharacters,PermutationOnUnipotents);
AUTO(ReadChv("work/util"),CheckRepresentations,CheckSchurRelations, 
  CheckHeckeSpecializes, FindRepresentation, GenericHecke, CheckOpdam);
AUTO(ReadChv("work/wordgen"),WordEnumerator);
AUTO(ReadChv("work/checkunip"),CheckUnipotentCentralizers,CheckNrSemisimple);
