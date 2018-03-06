_KBExtDir := Concatenation(LOADED_PACKAGES.kbmag,"bin/");

AUTO( ReadPkg( "kbmag", "gap", "fpgtorws" ),
  FpGroupToRWS, CaseChange );

AUTO( ReadPkg( "kbmag", "gap", "fsa" ),
  FSAOps, SetRecOps, IsFSA, IsInitializedFSA, InitializeSR, InitializeFSA,
  IsDeterministicFSA, WriteFSA, DenseDTableFSA, SparseTableFSA, BackTableFSA,
  TargetDFA, TargetsFSA, SourcesFSA, IsAcceptedWordDFA, WordTargetDFA,
  AddStateFSA, DeleteStateFSA, PermuteStatesFSA, AddLetterFSA,
  DeleteLetterFSA, PermuteLettersFSA, AddEdgeFSA, DeleteEdgeFSA,
  SetAcceptingFSA, SetInitialFSA, IsAccessibleFSA, AccessibleFSA, IsTrimFSA,
  TrimFSA, IsBFSFSA, BFSFSA, PSizeFSA, LSizeDFA, ListWordSR, WordListSR,
  LEnumerateDFA, SortLEnumerateDFA, SizeLEnumerateDFA );

AUTO( ReadPkg( "kbmag", "gap", "rws" ),
  _RWS, _KBTmpFileName, RWSOps, IsRWS, IsInternalRWS, IsConfluentRWS,
  IsAvailableNormalFormRWS, IsAvailableReductionRWS, SetOrderingRWS,
  ReorderGeneratorsRWS, ReadRWS, IsReducedWordRWS, ReduceWordRWS, SizeRWS,
  EnumerateRWS, SortEnumerateRWS, SizeEnumerateRWS, KB, ResetRWS,
  AddOriginalEqnsRWS, Automata );
