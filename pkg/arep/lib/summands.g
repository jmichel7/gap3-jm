# Recognizing the Block Structures of Matrices
# SE, MP, 14.02.96 - , GAPv3.4

# 14.2.96 erste Version
# 16.2.96 Konvention f"ur MatPerm wurde ge"andert
# 26.2.96 neue Programmstruktur f"ur "Ahnlichkeit
# 27.2.96 DirectSummands neu gebaut

# Terminology
#   * A matrix is a non-trivial direct sum iff is of the form
#     DiagonalMat(A1, .., Ar) where r > 1.
#   * A matrix is a permuted direct sum iff it is of the form
#     PermutedMat(sL, A, sR) where A is a direct sum of size 
#     (nr x nc) and sL, sR are permutations in S_nr, S_nc 
#     respectively acting on the rows and on the columns.

# Similarity under permutation of rows and columns
# ================================================

# Idee:
#   F"ur eine Matrix A wird PermSimilarityOrbitMat(A) konstruiert,
#   was eine Datenstruktur f"ur das Orbit ist. Dann kann mit '=', 
#   '<' und 'in' verglichen werden. Dabei wird der Repr"asentant
#   modifiziert.

# Initial fields in a PermSimilarityOrbit-record
#   isPermSimilarityOrbit : always true
#   representative        : a matrix
#   original              : [x, y] such that PermutedMat(x, A, y) is
#                           the original matrix given as representative
#   rows, columns         : size of the matrix
#   operations            : always PermSimilarityOrbitOps

if not IsBound(PermSimilarityOrbitOps) then
  PermSimilarityOrbitOps := OperationsRecord("PermSimilarityOrbitOps");
fi;

PermSimilarityOrbitMat := function ( A )
  if not IsMat(A) then
    Error("usage: PermSimilarityOrbitMat( <mat> )");
  fi;
  return
    rec(
      isPermSimilarityOrbit := true,
      representative        := A,
      original              := [ (), () ],
      rows                  := Length(A),
      columns               := Length(A[1]),
      operations            := PermSimilarityOrbitOps
    );
end;

PermSimilarityOrbitOps.PreNormalForm := function ( S )
  local nr, nc, A, tA,
        Tr, sr, Br,
        br, Abr, sbr,
        Tc, sc, Bc,
        bc, tAbc, sbc,
        i;

  if IsBound(S.isPreNF) then
    return;
  fi;

  A  := S.representative;
  nr := Length(A);
  nc := Length(A[1]);

  # A sortieren bzgl. Zeilenbilanz und Zeilengleichheit
  Tr := List(A, Collected);
  sr := Sortex(Tr);
  A  := Permuted(A, sr);
  Br := PartitionIndex(Tr);
  for br in Br do
    Abr := Sublist(A, br);
    sbr := Sortex(Abr);
    for i in [1..Length(br)] do
      A[br[i]] := Abr[i];
    od;
    sr := sr*sbr;
  od;
  S.representative    := A;
  S.original[1]       := S.original[1] * sr;
  S.rowTotals         := Collected(Tr);
  S.rowTotalPartition := PartitionIndex(Tr);

  # A sortieren bzgl. Spaltenbilanz und Spaltengleichheit
  tA := TransposedMat(S.representative);
  Tc := List(tA, Collected);
  sc := Sortex(Tc);
  tA := Permuted(tA, sc);
  Bc := PartitionIndex(Tc);
  for bc in Bc do
    tAbc := Sublist(tA, bc);
    sbc  := Sortex(tAbc);
    for i in [1..Length(bc)] do
      tA[bc[i]] := tAbc[i];
    od;
    sc := sbc*sc;
  od;
  S.representative       := TransposedMat(tA);
  S.original[2]          := sc^-1 * S.original[2];
  S.columnTotals         := Collected(Tc);
  S.columnTotalPartition := PartitionIndex(Tc);


  # check if the preNF is sufficient
  S.preNFSufficient := 
    ForAll(S.rowTotalPartition,    b -> Length(b) = 1) and
    ForAll(S.columnTotalPartition, b -> Length(b) = 1);

  # compose a structure as the fingerprint
  S.isPreNF     := true;
  S.fingerprint := 
    [ S.rows, 
      S.columns,
      S.rowTotals,
      S.rowTotalPartition,
      S.columnTotals,
      S.columnTotalPartition
    ];
end;

PermSimilarityOrbitOps.NormalForm := function ( S )
  local Sr, Sc, A0, x0, y0, x, y, xAy;

  if IsBound(S.isNF) then
    return;
  fi;
  PermSimilarityOrbitOps.PreNormalForm(S);

  # brute force enumeration
  Sr := DirectProductSymmetricGroups(S.rowTotalPartition);
  Sc := DirectProductSymmetricGroups(S.columnTotalPartition);
  Print("#I enumerating group of size ", Size(Sr)*Size(Sc), "\n");

  A0 := S.representative;
  x0 := ();
  y0 := ();
  for x in Elements(Sr) do
    for y in Elements(Sc) do
      xAy := PermutedMat(x, S.representative, y);
      if xAy < A0 then
        A0 := xAy;
        x0 := x;
        y0 := y;
      fi;
    od;
  od;

  S.isNF           := true;
  S.representative := A0;
  S.original       := [ x0*S.original[1], S.original[2]*y0 ];
end;


PermSimilarityOrbitOps.\= := function ( A, B )
  if A.rows <> B.rows or A.columns <> B.columns then
    return false;
  fi;

  PermSimilarityOrbitOps.PreNormalForm(A);
  PermSimilarityOrbitOps.PreNormalForm(B);
  if A.fingerprint <> B.fingerprint then
    return false;
  fi;
  if A.representative = B.representative then
    return true;
  fi;
  if A.preNFSufficient or B.preNFSufficient then
    return false;
  fi;

  PermSimilarityOrbitOps.NormalForm(A);
  PermSimilarityOrbitOps.NormalForm(B);
  return A.representative = B.representative;
end;

PermSimilarityOrbitOps.\in := function ( A, B )
  return PermSimilarityOrbitMat(A) = B;
end;

PermSimilarityOrbitOps.\< := function ( A, B )
  if A.rows    < B.rows    then return true;  fi;
  if A.rows    > B.rows    then return false; fi;
  if A.columns < B.columns then return true;  fi;
  if A.columns > B.columns then return false; fi;

  PermSimilarityOrbitOps.PreNormalForm(A);
  PermSimilarityOrbitOps.PreNormalForm(B);
  if A.fingerprint < B.fingerprint then return true;  fi;
  if A.fingerprint > B.fingerprint then return false; fi;
  if A.representative = B.representative then return false; fi;

  PermSimilarityOrbitOps.NormalForm(A);
  PermSimilarityOrbitOps.NormalForm(B);
  return A.representative < B.representative;
end;

PermSimilarityOrbitOps.Print := function ( A )
  Print(
    "PermSimilarityOrbitMat( ", 
    PermutedMat(A.original[1], A.representative, A.original[2]), 
    " )"
  );
end;


# Decomposition of permuted direct sums
# =====================================

#F DirectSummandsPermutedMat( <mat> [, <match_blocks> ])
#F   returns [sL, [A1, .., Ar], sR] where r > 0, sL, sR are
#F   permutations, and A1, .., Ar are matrices such that
#F   the given matrix is
#F     A = PermutedMat(sL, DirectSum(A1, .., Ar), sR)
#F   and r is maximal. The sizes of the Ak are non-decreasing.
#F   If the flag match_blocks is not provided or true then
#F   the permutations sL, sR are chosen such that equivalent 
#F   matrices (under permutations of rows and columns) become 
#F   equal and occur next to each other in the list [A1, .., Ar]. 
#F   This feature allows detection of tensor products of the 
#F   form TensorProductMat(IdentityMat(r), A1) and others types 
#F   of tensor products involving identity matrices but 
#F   it on may take a lot of time (for 0/1-matrices known to be 
#F   NP-complete, although mostly quite fast).
#F


DirectSummandsPermutedMat := function ( arg )
  local A, match_blocks,  # arguments
        nr, nc, nb,       # number of rows, columns, and blocks of A,
        As,               # list of direct summands of A
        sL,  sR,          # permutations to seperate the blocks
        rbs, cbs, bs,     # row-, column- and common block structure
        sB,               # permutation to rearrange As to sort fingerprints
        k;                # counter

  if Length(arg) = 1 and IsMat(arg[1]) then
    A            := arg[1];
    match_blocks := true;
  elif Length(arg) = 2 and IsMat(arg[1]) and IsBool(arg[2]) then
    A            := arg[1];
    match_blocks := arg[2];
  else
    Error("usage: DirectSummandsPermutedMat( <mat> [, <match-blocks> ])");
  fi;
  nr := Length(A);
  nc := Length(A[1]);

  # compute the block structure
  bs  := BlockStructureMat(A);
  rbs := bs[1];
  cbs := bs[2];
  nb  := Length(rbs); # = Length(cbs)
  if nb = 1 then
    return [(), [ A ], ()];
  fi;

  # seperate the blocks
  As :=
    List(
      [1..nb], 
      k -> List(Sublist(A, rbs[k]), Ar -> Sublist(Ar, cbs[k]))
    );
  if not match_blocks then

    # just compute the seperating permutations and return
    sL := PermList( Concatenation( rbs ) )^-1;
    sR := PermList( Concatenation( cbs ) );
    return [ sL, As, sR ];
  fi;

  # match the blocks and sort them
  As := List(As, PermSimilarityOrbitMat);
  sB := Sortex(As);

  # compute sL, sR to seperate the blocks as sorted in As
  bs  := List(bs, bsi -> Permuted(bsi, sB));
  rbs := bs[1];
  cbs := bs[2];
  sL  := PermList( Concatenation( rbs ) )^-1;
  sR  := PermList( Concatenation( cbs ) );

  # sort within the blocks to get the representative
  for k in [1..nb] do
    sL :=    
      MappingPermListList(
        Permuted(rbs[k], As[k].original[1]),
        rbs[k]
      ) * 
      sL;
    sR :=
      sR *
      MappingPermListList(
        Permuted(cbs[k], As[k].original[2]),
        cbs[k]
      );
  od;
  
  return [sL, List(As, Ak -> Ak.representative), sR];
end;

# Beispiel
# A1 := # 1_3 tensor [[1,2],[3,1]]
#  [ [ 0, 0, 0, 2, 0, 1],
#    [ 3, 1, 0, 0, 0, 0],
#    [ 0, 0, 1, 0, 2, 0],
#    [ 1, 2, 0, 0, 0, 0],
#    [ 0, 0, 0, 1, 0, 3],
#    [ 0, 0, 3, 0, 1, 0] ];

# schlechtes Beispiel:
# [ [ 1, 1, 1, 0, 0, 0 ], [ 1, E(3), E(3)^2, 0, 0, 0 ], 
#   [ 1, E(3)^2, E(3), 0, 0, 0 ], [ 0, 0, 0, 1, 1, 1 ], 
#   [ 0, 0, 0, 1, E(3), E(3)^2 ], [ 0, 0, 0, 1, E(3)^2, E(3) ] ]
