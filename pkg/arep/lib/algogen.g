# Decomposition of Matrices with Symmetry
# MP, 12.12.97 - 2.8.98, GAPv3.4

# Literature:
#   T. Minkwitz: PhD. Thesis, University of Karlsruhe, 1993
#   S. Egner   : PhD. Thesis, University of Karlsruhe, 1997
#   M.Pueschel : PhD. Thesis, University of Karlsruhe, 1998

if not IsBound(InfoAlgogen) then
  InfoAlgogen := Ignore;  # switch it on if you like
fi;

#F MatrixDecompositionByPermPermSymmetry( <mat/amat> )
#F   decomposes <mat/amat> into a product of sparse
#F   matrices according to the perm-perm symmetry.
#F   An amat is returned representing the product.
#F

MatrixDecompositionByPermPermSymmetry := function ( M )
  local Rs, RL, RR, DL, DR, P, AL, AR, A;

  if not( IsMat(M) or IsAMat(M) ) then
    Error("usage: MatrixDecompositionByPermPermSymmetry( <mat/amat> )");
  fi;

  # calculate the permperm symmetry
  InfoAlgogen("#I CALCULATING THE PERMPERM SYMMETRY\n");
  Rs := PermPermSymmetry(M);
  RL := Rs[1];
  RR := Rs[2];

  # decompose the reps
  # if the reps are equivalent by 
  # a permutation, only one of 
  # them has to be decomposed
  InfoAlgogen("#I DECOMPOSING THE REPS\n");
  DL := DecompositionMonRep(RL);
  AL := DL.conjugation.element;
  if IsEquivalentARep(RL, RR) then
    P := ConjugationPermReps(RL, RR); 
    if P <> false then
      if IsIdentityMat(P) then
	AR := AL;
      else
        AR := P * AL;
      fi;
    else
      DR := DecompositionMonRep(RR);
      AR := DR.conjugation.element;
    fi;
  fi;

  # the block matrix
  InfoAlgogen("#I SPECIALIZING\n");
  if IsMat(M) then
    M := AMatMat(M);
  fi;
  A := 
    AL *
    AMatSparseMat(MatAMat(InverseAMat(AL) * M * AR), false) *
    InverseAMat(AR);

  # avoid that A is checked for monomiality
  # by SimplifyAMat
  A.isMonMat := IsMonMat(M);

  return SimplifyAMat(A);
end;


#F MatrixDecompositionByMonMonSymmetry( <mat/amat> )
#F   decomposes <mat/amat> into a product of sparse
#F   matrices according to the mon-mon symmetry.
#F   An amat is returned representing the product.
#F

MatrixDecompositionByMonMonSymmetry := function ( M )
  local Rs, RL, RR, DL, DR, AL, AR, A;

  if not( IsMat(M) or IsAMat(M) ) then
    Error("usage: MatrixDecompositionByMonMonSymmetry( <mat/amat> )");
  fi;

  # calculate the monmon symmetry
  InfoAlgogen("#I CALCULATING THE MONMON SYMMETRY\n");
  Rs := MonMonSymmetry(M);
  RL := Rs[1];
  RR := Rs[2];

  # decompose the reps
  InfoAlgogen("#I DECOMPOSING THE SYMMETRY\n");
  DL := DecompositionMonRep(RL);
  DR := DecompositionMonRep(RR);
  AL := DL.conjugation.element;
  AR := DR.conjugation.element;

  # the block matrix
  InfoAlgogen("#I SPECIALIZING\n");
  if IsMat(M) then
    M := AMatMat(M);
  fi;
  A := 
    AL *
    AMatSparseMat(MatAMat(InverseAMat(AL) * M * AR), false) *
    InverseAMat(AR);

  # avoid that A is checked for monomiality
  # by SimplifyAMat
  A.isMonMat := IsMonMat(M);

  return SimplifyAMat(A);
end;


#F MatrixDecompositionByPermIrredSymmetry( <mat/amat> [, <maxblocksize> ] )
#F   decomposes <mat/amat> into a product of sparse
#F   matrices according to the perm-irred symmetry
#F   returned by the function PermIrredSymmetry1.
#F   An amat is returned representing the product.
#F   Only those symmetries where all irreducibles 
#F   are of degree <= <maxblocksize> are considered. 
#F   The default for <maxblocksize> is 2.
#F   Among all symmetries [RL, RR] the best is chosen 
#F   according to the following measure.
#F   If
#F     RL ~= RR ~= directsum R_i^(n_i), and Q = sum n_i^2 * d_i^2,
#F   with d_i = deg(R_i), then the best symmetry has
#F   the smallest value of Q.
#F

MatrixDecompositionByPermIrredSymmetry := function ( arg )
  local M, max, dim, Rs, Qs, min, pos, R, D, A;

  # decode and check arguments
  if Length(arg) = 1 then
    M   := arg[1];
    max := 2;
  elif Length(arg) = 2 then
    M   := arg[1];
    max := arg[2];
  else
    Error(
      "usage: \n", 
      "  MatrixDecompositionByPermIrredSymmetry( ",
      "    <mat/amat> [, <maxblocksize> ]  )"
    );
  fi;
  if IsAMat(M) then
    M := MatAMat(M);
  fi;
  if not IsMat(M) then
    Error("<M> must be a matrix");
  fi;
  if not ( IsInt(max) and max >= 1 ) then
    Error("<max> must be a posiive integer");
  fi;

  dim := DimensionsMat(M);
  if dim[1] <> dim[2] then
    return AMatMat(M);
  fi;

  # calculate the symmetry
  # use PermIrredSymmetry1 for time reasons
  InfoAlgogen("#I CALCULATING THE PERMIRRED SYMMETRY\n");
  Rs := PermIrredSymmetry1(M, max);
  if Length(Rs) > 0 then
    
    # take the best pair for decomposition,
    # let R = directsum R_i^(n_i) be the decomposition
    # of the right (or left) side of the perm-irred symmetry
    # into irreducibles, then the quality is given by 
    # a small value of 
    #   sum n_i^2 * d_i^2, 
    # where d_i = deg(R_i)
    InfoAlgogen("#I CHOOSING SYMMETRY\n");
    Qs := 
      List(
        Rs, 
        p -> 
          Sum(
	    List(
	      Collected(
		List(p[2].rep.summands, r -> CharacterARep(r))
	      ),
	      cn -> Degree(cn[1])^2 * cn[2]^2
	    )
          )
      );
    min := Minimum(Qs);
    pos := PositionProperty(Qs, q -> q = min);
    R   := Rs[pos];

    InfoAlgogen("#I DECOMPOSING THE SYMMETRY\n");
    D := DecompositionMonRep(R[1]).conjugation.element;

    # the block matrix
    InfoAlgogen("#I SPECIALIZING\n");
    M := AMatMat(M);
    A := 
      D *
      AMatSparseMat(
        MatAMat(InverseAMat(D) * M * InverseAMat(R[2].conjugation)), 
        false
      ) *
      R[2].conjugation;

    # avoid that A is checked for monomiality
    # by SimplifyAMat
    A.isMonMat := IsMonMat(M);

    return SimplifyAMat(A);
  fi;

  # decomposition failed
  return AMatMat(M);

end;



