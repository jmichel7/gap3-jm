# Determination of Symmetry as a Pair of Representations
# MP, 12.12.97, GAPv3.4

# Literature:
#   T. Minkwitz: PhD. Thesis, University of Karlsruhe, 1993
#   S. Egner   : PhD. Thesis, University of Karlsruhe, 1997
#   M.Pueschel : PhD. Thesis, University of Karlsruhe, 1998

#F Computing Symmetry of Matrices
#F ==============================
#F

#F PermPermSymmetry( <mat/amat> )
#F   calculates the perm-perm symmetry of the <mat/amat> M
#F   as a pair [R1, R2] of "perm"-areps of the same group G.
#F   This means, that M is in the intertwining space 
#F   Int(R1, R2), i.e.
#F     R1(g) * M = M * R2(g), for all g in G.
#F   If G is solvable, then it is an aggroup.
#F

PermPermSymmetry := function ( A )
  local G, char, Gag, phi;

  # check argument
  if IsAMat(A) then
    A := MatAMat(A);
  fi;
  if not IsMat(A) then
    Error("<A> must be a matrix");
  fi;

  # calculate perm-perm symmetry
  G    := PermPermSym(A);
  char := Characteristic(DefaultField(A[1][1]));

  # switch to an aggroup, if possible
  if IsSolvable(G) then
    Gag := GroupWithGenerators(AgGroup(G));
    phi := Gag.bijection;

    return
      [ ARepByImages(
          Gag,
          List(
            Gag.theGenerators, 
            g -> PermPermSymL(G, Image(phi, g))
          ),
          G.dimensionsMat[1],
          char,
          "hom"
        ),
        ARepByImages(
          Gag,
          List(
            Gag.theGenerators, 
            g -> PermPermSymR(G, Image(phi, g))
          ),
          G.dimensionsMat[2],
          char,
          "hom"
        )
      ];
  else
    G := GroupWithGenerators(G);

    return
      [ ARepByImages(
          G,
          PermPermSymL(G, G.theGenerators),
          G.dimensionsMat[1],
          char,
          "hom"
        ),
        ARepByImages(
          G,
          PermPermSymR(G, G.theGenerators),
          G.dimensionsMat[2],
          char,
          "hom"
        )
      ]; 
  fi;
end;


#F MonMonSymmetry( <mat/amat> )
#F   calculates the mon-mon symmetry of the <mat/amat> A
#F   as a pair [R1, R2] of "mon"-areps of the same group G.
#F   This means, that A is in the intertwining space 
#F   Int(R1, R2), i.e.
#F     R1(g) * A = A * R2(g), for all g in G.
#F   If G is solvable, then it is an aggroup.
#F   The function only works for characteristic zero.
#F

MonMonSymmetry := function ( A )
  local G, char, Gag, phi;

  # check argument
  if IsAMat(A) then
    A := MatAMat(A);
  fi;
  if not IsMat(A) then
    Error("<A> must be a matrix");
  fi;
  char := Characteristic(DefaultField(A[1][1]));

  # require char = 0
  if char <> 0 then
    Error("characteristic of <A> must be zero");
  fi;

  # calculate perm-perm symmetry
  G    := MonMonSym(A);

  # switch to an aggroup, if possible
  if IsSolvable(G) then
    Gag := GroupWithGenerators(AgGroup(G));
    phi := Gag.bijection;

    return
      [ ARepByImages(
          Gag,
          List(
            Gag.theGenerators, 
            g -> MonMonSymL(G, Image(phi, g))
          ),
          "hom"
        ),
        ARepByImages(
          Gag,
          List(
            Gag.theGenerators, 
            g -> MonMonSymR(G, Image(phi, g))
          ),
          "hom"
        )
      ];
  else
    G := GroupWithGenerators(G);

    return
      [ ARepByImages(
          G,
          MonMonSymL(G, G.theGenerators),
          "hom"
        ),
        ARepByImages(
          G,
          MonMonSymR(G, G.theGenerators),
          "hom"
        )
      ]; 
  fi;
end;


#F PermIrredSymmetry( <mat/amat> [, <maxblocksize> ] )
#F   calculates a list containing all the non-trivial 
#F   perm-irred symmetries of <mat/amat>, where the degree
#F   of at least one irreducible is <= <maxblocksize>. 
#F   The default for <maxblocksize> is 2 because of the 
#F   expensive computation.
#F   A perm-irred symmetry is a pair [R1, R2] of areps of the same 
#F   group G, where R1 is a "perm"-arep and R2 is a direct sum 
#F   of irreducible "mat"-areps, conjugated by a "perm"-amat.
#F   M = <mat/amat> is in the intertwining space of R1 and R2, i.e.
#F     R1(g) * M = M * R2(g), for all g in G.
#F   If G is solvable, then it is an aggroup.
#F   The function only works for characteristic zero.
#F

PermIrredSymmetry := function ( arg )
  local A, k, dim, char, Gs, G, RG, Rs, irrs, nrirrs, b, perm;

  # decode and check arguments
  if Length(arg) = 1 then
    A := arg[1];
    if IsAMat(A) then
      A := MatAMat(A);
    fi;
    if not IsMat(A) then
      Error("<A> must be a matrix or an amat");
    fi;
    dim := DimensionsMat(A);
    if dim[1] <> dim[2] then
      return AMatMat(A);
    fi;
    k := 2;
  elif Length(arg) = 2 then
    A := arg[1];
    k := arg[2];
    if IsAMat(A) then
      A := MatAMat(A);
    fi;
    if not IsMat(A) then
      Error("<A> must be a matrix or an amat");
    fi;
    dim := DimensionsMat(A);
    if dim[1] <> dim[2] then
      return AMatMat(A);
    fi;
  else
    Error("usage: PermIrredSymmetry( <mat/amat> [, <maxblocksize> ] )");
  fi;
  char := Characteristic(DefaultField(A[1][1]));

  # require char = 0
  if char <> 0 then
    Error("characteristic of <A> must be zero");
  fi;

  # calculate non-trivial perm-block symmetry
  Gs := 
    Filtered(
      PermBlockSymBySubsets(A, [1..k]),
      G -> not IsTrivial(G)
    );
  Gs := List(Gs, GroupWithGenerators);

  # check each one for irreducibility
  # compute first the number of irreducibles
  # components of the right rep and compare
  # with length of the kbs
  Rs := [ ];
  for G in Gs do
    RG := NaturalARep(G, dim[1], char);
    nrirrs := 
      Sum(
        List(
          Irr(G),
          chi -> ScalarProduct(chi, CharacterARep(RG))
        )
      );
    if nrirrs = Length(G.kbsM) then
      irrs := [ ];
      for b in G.kbsM do
        Add(
          irrs,
          ARepByImages(
            G,
            List(
              G.theGenerators,
              g -> PermBlockSymR(G, b, g)
            ),
            "hom"
          )
        );
      od;
      perm := 
        MappingPermListList(
          [1..dim[1]],
          Concatenation(G.kbsM)
        );
      Add(
        Rs, 
        [ RG,
          ConjugateARep(
            DirectSumARep(irrs),
            AMatPerm(perm, dim[1], char)
          )
        ]
      );
    fi;
  od;

  return Rs;  
end;


#F PermIrredSymmetry1( <mat/amat> [, <maxblocksize> ] )
#F   calculates a list containing all the non-trivial 
#F   perm-irred symmetries of <mat/amat>, where the degrees
#F   of all irreducibles is <= <maxblocksize>
#F   (This is the difference to PermIrredSymmetry).
#F   The default for <maxblocksize> is 2 because of the 
#F   expensive computation.
#F   A perm-irred symmetry is a pair [R1, R2] of areps of the same 
#F   group G, where R1 is a "perm"-arep and R2 is a direct sum 
#F   of irreducible "mat"-areps, conjugated by a "perm"-amat.
#F   M = <mat/amat> is in the intertwining space of R1 and R2, i.e.
#F     R1(g) * M = M * R2(g), for all g in G.
#F   If G is solvable, then it is an aggroup.
#F   The function only works for characteristic zero.
#F

PermIrredSymmetry1 := function ( arg )
  local A, k, dim, char, Gs, G, RG, Rs, irrs, nrirrs, b, perm;

  # decode and check arguments
  if Length(arg) = 1 then
    A := arg[1];
    if IsAMat(A) then
      A := MatAMat(A);
    fi;
    if not IsMat(A) then
      Error("<A> must be a matrix or an amat");
    fi;
    dim := DimensionsMat(A);
    if dim[1] <> dim[2] then
      return AMatMat(A);
    fi;
    k := 2;
  elif Length(arg) = 2 then
    A := arg[1];
    k := arg[2];
    if IsAMat(A) then
      A := MatAMat(A);
    fi;
    if not IsMat(A) then
      Error("<A> must be a matrix or an amat");
    fi;
    dim := DimensionsMat(A);
    if dim[1] <> dim[2] then
      return AMatMat(A);
    fi;
  else
    Error("usage: PermIrredSymmetry1( <mat/amat> [, <maxblocksize> ] )");
  fi;
  char := Characteristic(DefaultField(A[1][1]));

  # require char = 0
  if char <> 0 then
    Error("characteristic of <A> must be zero");
  fi;

  # calculate non-trivial perm-block symmetry
  Gs := 
    Filtered(
      PermBlockSymBySubsets(A, [1..k]),
      G -> ForAll(G.kbsM, b -> Length(b) <= k) and not IsTrivial(G)
    );
  Gs := List(Gs, GroupWithGenerators);

  # check each one for irreducibility
  # compute first the number of irreducibles
  # components of the right rep and compare
  # with length of the kbs
  Rs := [ ];
  for G in Gs do
    RG := NaturalARep(G, dim[1], char);
    nrirrs := 
      Sum(
        List(
          Irr(G),
          chi -> ScalarProduct(chi, CharacterARep(RG))
        )
      );
    if nrirrs = Length(G.kbsM) then
      irrs := [ ];
      for b in G.kbsM do
        Add(
          irrs,
          ARepByImages(
            G,
            List(
              G.theGenerators,
              g -> PermBlockSymR(G, b, g)
            ),
            "hom"
          )
        );
      od;
      perm := 
        MappingPermListList(
          [1..dim[1]],
          Concatenation(G.kbsM)
        );
      Add(
        Rs, 
        [ RG,
          ConjugateARep(
            DirectSumARep(irrs),
            AMatPerm(perm, dim[1], char)
          )
        ]
      );
    fi;
  od;

  return Rs;  
end;

