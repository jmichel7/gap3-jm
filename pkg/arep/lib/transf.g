# Discrete Linear Signaltransforms
# SE, MP ab 27.02.96, GAP v3.4

# 27.02.96: Erste Version
# 02.04.96: Hartley-T.; DCT/DCT-IV
# 06.02.97: char > 0, ADFT, PADFT, RecognizeADFT
# 07.02.97: Korrektur: PADFT baut nicht unbedingt eine Basis
# 03.12.97: DCT_I
# 27.05.98: MP, Sinustransformationen

#F Discrete Signal Transforms
#F ==========================
#F

# Literature:
#   [1] S. Lang: Algebra, 2nd ed.; Addison-Wesley 1984. (HA Beth, D Lan).
#   [2] Elliott, Rao: Fast Transforms.
#   [3] Clausen, Baum: FFT.
#   [4] H. S. Malvar: Signal processing with lapped transforms.
#   [5] A. N"uckel, A. Klappenecker: On the Parametrization of
#       Algebraic Discrete Fourier Transforms. Extended Abstract
#       submitted to EUROCAST '97.
#   [6] A. Mertins: Signaltheorie

# Discrete Fourier Transform (DFT)
# ================================

#F DiscreteFourierTransform( <rootOfUnity> )
#F DiscreteFourierTransform( <N> [, <char>] )
#F InverseDiscreteFourierTransform( <rootOfUnity> )
#F InverseDiscreteFourierTransform( <N>, [, <char>] )
#F   constructs the Discrete Fourier Transform on N points
#F   or its inverse. (Elliott, Rao, 3.) If a prime int 
#F   <char> is supplied, then a DFT in that characteristic 
#F   is constructed if possible. If <rootOfUnity> is supplied
#F   then the transform of length Order(<rootOfUnity>) is
#F   constructed.
#F

DiscreteFourierTransform := function ( arg )
  local N, p, w, W, k;

  if 
    Length(arg) = 1 and IsInt(arg[1]) and arg[1] >= 1 or
    Length(arg) = 2 and IsInt(arg[1]) and arg[1] >= 1 and
    arg[2] = 0
  then
    N := arg[1];
    p := 0;
    w := ExpIPi(2/N);
  elif Length(arg) = 1 and arg[1] in FieldElements then
    w := arg[1];
    p := Characteristic(DefaultField(w));
    if p = 0 then
      N := OrderCyc(w);
    else
      N := OrderFFE(w);
    fi;
    if N = "infinity" then
      Error("<w> must be a root of unity");
    fi;
  elif 
    Length(arg) = 2 and 
    IsInt(arg[1]) and arg[1] >= 1 and
    IsInt(arg[2]) and arg[2] >= 2 and IsPrimeInt(arg[2])
  then
    N := arg[1];
    p := arg[2];
    k := 1;
    while not (p^k - 1) mod N = 0 do
      k := k + 1;
      if not p^k <= 2^16 then
        Error("cannot construct primitive <N>-th root of unity");
      fi;
    od;
    w := Z(p^k)^QuoInt(p^k - 1, N);
  else
    Error("wrong arguments");
  fi;

  W := [ w^0 ];
  for k in [1 .. N-1] do
    W[k+1] := W[k] * w;
  od;
  return 
    List([0..N-1], i -> List([0..N-1], j -> 
      W[((i*j) mod N)+1]
    ));
end;

InverseDiscreteFourierTransform := function ( arg )
  local F, N, NInv, k, Fk;

  if Length(arg) = 1 then
    F := DiscreteFourierTransform(arg[1]);
  elif Length(arg) = 2 then
    F := DiscreteFourierTransform(arg[1], arg[2]);
  else
    Error("wrong arguments");
  fi;

  N := Length(F);
  for k in [2..QuoInt(N+1, 2)] do
    Fk       := F[k];
    F[k]     := F[N-k+2];
    F[N-k+2] := Fk;
  od;

  NInv := (N * F[1][1]^0)^-1;
  for k in [1..N] do
    F[k] := NInv * F[k];
  od;
  return F;
end;

# Beispiel
# --------
#
# P23 := MatPerm((2,3), 4);
# I2  := IdentityMat(2);
# F2  := [[1,1], [E(4), -E(4)]]; 
# DFT(4) = 
#   P23 * TensorProductMat(I2, DFT(2)) * 
#   P23 * DiagonalMat(DFT(2), F2) * P23;


# Discrete Hartley Transform (DHT)
# ================================

#F DiscreteHartleyTransform( <N> )
#F InverseDiscreteHartleyTransform( <N> )
#F   constructs the Discrete Hartley Transform on N points 
#F   or its inverse. (Malvar, 1.2.3)
#F

DiscreteHartleyTransform := function ( N )
  local H, i, j;

  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: DiscreteHartleyTransform( <positive-int> )");
  fi;

  H := [];
  for i in [N-1, N-2 .. 0] do
    for j in [N-1, N-2 .. 0] do
      if not IsBound(H[1+i*j]) then
        H[1+i*j] := 1/Sqrt(N)*(CosPi(2*i*j/N) + SinPi(2*i*j/N));
      fi;
    od;
  od;
  return List([0..N-1], i -> List([0..N-1], j -> H[1+i*j]));
end;

InverseDiscreteHartleyTransform := 
  DiscreteHartleyTransform;


# Discrete Cosine/Sine Transforms 
# ===============================

# DCT, DST, type-I, type-II and type-IV
# note that type-III is the transposed of type-II

#F DiscreteCosineTransform( <N> )
#F InverseDiscreteCosineTransform( <N> )
#F   constructs the standard Discrete Cosine Transform (type-II) 
#F   on N points or its inverse. (Malvar, 1.2.5)
#F

DiscreteCosineTransform := function ( N )
  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: DiscreteCosineTransform( <positive-int> )");
  fi;
  return
    List([0..N-1], i -> List([0..N-1], function (j)
      if i = 0 then
        return 1/Sqrt(2) * Sqrt(2/N) * CosPi((j + 1/2)*i/N);
      else
        return 1         * Sqrt(2/N) * CosPi((j + 1/2)*i/N);
      fi;
    end ) ); 
end;

InverseDiscreteCosineTransform := function ( N )
  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: InverseDiscreteCosineTransform( <positive-int> )");
  fi;
  return TransposedMat( DiscreteCosineTransform( N ) );
end;


#F DiscreteSineTransform( <N> )
#F   constructs the standard Discrete Sine Transform (type-II) 
#F   on N points. Note that the matrix is rectangular and 
#F   has size (n-1) x n.
#F

DiscreteSineTransform := function ( N )
  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: DiscreteSineTransform( <positive-int> )");
  fi;
  return
    List([1..N-1], i -> List([0..N-1], function (j)
      return Sqrt(2/N) * SinPi((j + 1/2)*i/N);
    end ) ); 
end;

#F DiscreteCosineTransformIV( <N> )
#F InverseDiscreteCosineTransformIV( <N> )
#F   constructs the Type-IV Discrete Cosine Transform on N points
#F   or its inverse. (Malvar, 1.2.6)
#F

DiscreteCosineTransformIV := function ( N )
  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: DiscreteCosineTransformIV( <positive-int> )");
  fi;
  return
    List([0..N-1], i -> List([0..N-1], j ->
      Sqrt(2/N) * CosPi((i + 1/2)*(j + 1/2)/N)
    ) );
end;

InverseDiscreteCosineTransformIV := function ( N )
  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: InverseDiscreteCosineTransformIV( <positive-int> )");
  fi;
  return TransposedMat( DiscreteCosineTransformIV( N ) );
end;


#F DiscreteSineTransformIV( <N> )
#F   constructs the Type-IV Discrete Sine Transform on N points
#F   or its inverse. (Malvar, 1.2.6)
#F

DiscreteSineTransformIV := function ( N )
  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: DiscreteSineTransformIV( <positive-int> )");
  fi;
  return
    List([0..N-1], i -> List([0..N-1], j ->
      Sqrt(2/N) * SinPi((i + 1/2)*(j + 1/2)/N)
    ) );
end;


#F DiscreteCosineTransformI( <N> )
#F InverseDiscreteCosineTransformI( <N> )
#F   constructs the Type-I Discrete Cosine Transform on N + 1 
#F   points or its inverse. (Mertins, 2.4.4)
#F

DiscreteCosineTransformI := function ( N )
  local M, i;

  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: DiscreteCosineTransformI( <positive-int> )");
  fi;
  M := List([0..N], i -> List([0..N], j -> Sqrt(2/N) * CosPi(i*j/N)));
  M[1]     := 1/Sqrt(2) * M[1];
  M[N + 1] := 1/Sqrt(2) * M[N + 1];
  for i in [1..N + 1] do
    M[i][1]     := 1/Sqrt(2) * M[i][1];
    M[i][N + 1] := 1/Sqrt(2) * M[i][N + 1];
  od;

  return M;
end;

InverseDiscreteCosineTransformI := function ( N )
  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: InverseDiscreteCosineTransformI( <positive-int> )");
  fi;
  return TransposedMat( DiscreteCosineTransformI( N ) );
end;


#F DiscreteSineTransformI( <N> )
#F   constructs the Type-I Discrete Sine Transform on N + 1 
#F   points or its inverse. Note that the resulting matrix
#F   has size n x n.

DiscreteSineTransformI := function ( N )
  local M, i;

  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: DiscreteSineTransformI( <positive-int> )");
  fi;
  return List([1..N], i -> List([1..N], j -> Sqrt(2/N) * SinPi(i*j/N)));
end;


# Walsh-Hadamard Transformation (WHT)
# ===================================

#F WalshHadamardTransform( <N> )
#F InverseWalshHadamardTransform( <N> )
#F   constructs the Walsh-Hadamard Transform on N points
#F   or its inverse. For N being a power of 2 this is
#F   DFT(2)^(tensor n). (Clausen, Baum, 1.7; Elliott, Rao, 8.3) 
#F   For general N the transform is defined to be a tensor product 
#F   of Fourier transforms of prime sizes (sorted by size).
#F

WalshHadamardTransform := function ( N )
  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: WalshHadamardTransform( <positive-int> )");
  fi;
  if N = 1 then
    return [[1]];
  fi;
  return # [3], Kap. 1.7, S. 25., generalized
    TensorProductMat( 
      List(Factors(N), DiscreteFourierTransform) 
    );
end;

InverseWalshHadamardTransform := function ( N )
  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: InverseWalshHadamardTransform( <positive-int> )");
  fi;
  if N = 1 then
    return [[1]];
  fi;
  return 
    TensorProductMat( 
      List(Factors(N), InverseDiscreteFourierTransform) 
    );
end;


# Slant Transform (ST)
# ====================

#F SlantTransform( <N> )
#F InverseSlantTransform( <N> )
#F   constructs the Slant Transform of sidelength N or its 
#F   inverse. N must be a power of 2. (Elliott, Rao, 10.9)
#F

if not IsBound(SlantTransformTable) then
  SlantTransformTable := [

    # L = 1 from [2], (10.43)
    1/Sqrt(2) * [ 
      [ 1,  1 ],
      [ 1, -1 ]
    ],
    
    # L = 2 from [2], (10.45)
    1/2 * [       
      [ 1,  1,  1,  1 ],
      [ 3,  1, -1, -3 ] * 1/Sqrt(5),
      [ 1, -1, -1,  1 ],
      [ 1, -3,  3, -1 ] * 1/Sqrt(5)
    ]

  ];
fi;

SlantTransform := function ( N )
  local L,          # N = 2^L 
        I,          # identity matrix of size N/2-2
        aSqr, bSqr, # aSqr[k] = a_{2^k}^{2} as in [2]; b analogue
        aN, bN,     # a_N, b_N as in [2]
        T,          # composite factor 
        T11, T12,   # parts of T
        T21, T22,
        S,          # slant transform of sidelength 2^(L-1)
        k;          # counter

  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: SlantTransform( <positive-int> )");
  fi;
  L := LogInt(N, 2);

  # already computed before?
  if L = 0 then
    return [[1]];
  fi;
  if IsBound(SlantTransformTable[L]) then
    return List(SlantTransformTable[L], ShallowCopy);
  fi;

  # general case of sidelength N = 2^L >= 8, [2], (10.47)
  I := IdentityMat(N/2 - 2);

  aSqr := [ 1 ];
  bSqr := [ 1 ];
  for k in [2..L] do
    bSqr[k] := 1/(1 + 4*aSqr[k-1]);
    aSqr[k] := 4*bSqr[k]*aSqr[k-1];
  od;
  aN  := Sqrt( aSqr[L] );
  bN  := Sqrt( bSqr[L] );

  T11 := DiagonalMat([[ 1,  0 ], [  aN, bN ]],  I);
  T12 := DiagonalMat([[ 1,  0 ], [ -aN, bN ]],  I);
  T21 := DiagonalMat([[ 0,  1 ], [ -bN, aN ]],  I);
  T22 := DiagonalMat([[ 0, -1 ], [  bN, aN ]], -I);
  T   :=
    TensorProductMat([[1,0],[0,0]], T11) +
    TensorProductMat([[0,1],[0,0]], T12) +
    TensorProductMat([[0,0],[1,0]], T21) +
    TensorProductMat([[0,0],[0,1]], T22);

  S := SlantTransform(2^(L-1));

  SlantTransformTable[L] := 
    1/Sqrt(2) * T * DiagonalMat(S, S);

  return List(SlantTransformTable[L], ShallowCopy);
end;

InverseSlantTransform := function ( N )
  return TransposedMat( SlantTransform( N ) );
end;


# Haar Transform (HT)
# ===================

#F HaarTransform( <N> )
#F InverseHaarTransform( <N> )
#F   constructs the Haar Transform on N points or its 
#F   inverse. N must be a power of 2. (Elliott, Rao, 10.10)
#F

if true or not IsBound(HaarTransformTable) then
  HaarTransformTable := [

    # Ha(1) from [2], (10.53)
    [ [ 1,  1 ],
      [ 1, -1 ] 
    ],

    # Ha(2) from [2], (10.53)
    [ [ 1,  1,  1,  1 ],
      [ 1,  1, -1, -1 ],
      [ 1, -1,  0,  0 ] * Sqrt(2),
      [ 0,  0,  1, -1 ] * Sqrt(2) 
    ],

    # Ha(3) from [2], (10.53)
    [ [ 1,  1,  1,  1,  1,  1,  1,  1 ],
      [ 1,  1,  1,  1, -1, -1, -1, -1 ],
      [ 1,  1, -1, -1,  0,  0,  0,  0 ] * Sqrt(2),
      [ 0,  0,  0,  0,  1,  1, -1, -1 ] * Sqrt(2),
      [ 2, -2,  0,  0,  0,  0,  0,  0 ],
      [ 0,  0,  2, -2,  0,  0,  0,  0 ],
      [ 0,  0,  0,  0,  2, -2,  0,  0 ],
      [ 0,  0,  0,  0,  0,  0,  2, -2 ]
    ]

  ];
fi;

HaarTransform := function ( N )
  local L;

  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: HaarTransform( <positive-int> )");
  fi;
  L := LogInt(N, 2);
  if not 2^L = N then
    Error("<N> must be a power of 2");
  fi;
  if L = 0 then
    return [[1]];
  fi;
  if IsBound(HaarTransformTable[L]) then
    return 1/N * List(HaarTransformTable[L], ShallowCopy);
  fi;

  HaarTransformTable[L] := # [2], (10.54)
    Concatenation(
      TensorProductMat(N/2 * HaarTransform(N/2),     [[1,  1]]),
      TensorProductMat(Sqrt(N/2) * IdentityMat(N/2), [[1, -1]])
    );

  return 1/N * List(HaarTransformTable[L], ShallowCopy);
end;

InverseHaarTransform := function ( N )
  local L;

  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: InverseHaarTransform( <positive-int> )");
  fi;
  L := LogInt(N, 2);
  if not 2^L = N then
    Error("<N> must be a power of 2");
  fi;
  return N * TransposedMat( HaarTransform( N ) );
end;


# Rationalized Haar Transform (RHT)
# =================================

#F RationalizedHaarTransform( <N> )
#F InverseRationalizedHaarTransform( <N> )
#F   constructs the Tationalized Haar Transform on N points or its
#F   inverse. N must be a power of 2. (Elliott, Rao, 10.11).
#F

if true or not IsBound(RationalizedHaarTransformTable) then
  RationalizedHaarTransformTable := [

    # Rh(1) from [2], (10.58)
    [ [ 1,  1 ],
      [ 1, -1 ] 
    ]

  ];
fi;

RationalizedHaarTransform := function ( N )
  local L;

  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: RationalizedHaarTransform( <positive-int> )");
  fi;
  L := LogInt(N, 2);
  if not 2^L = N then
    Error("<N> must be a power of 2");
  fi;
  if L = 0 then
    return [[1]];
  fi;
  if IsBound(RationalizedHaarTransformTable[L]) then
    return List(RationalizedHaarTransformTable[L], ShallowCopy);
  fi;

  RationalizedHaarTransformTable[L] := # from [2] (induction)
    Concatenation(
      TensorProductMat(RationalizedHaarTransform(N/2), [[1,  1]]),
      TensorProductMat(IdentityMat(N/2),               [[1, -1]])
    );

  return List(RationalizedHaarTransformTable[L], ShallowCopy);
end;

InverseRationalizedHaarTransform := function ( N )
  local L, R, k;

  if not ( IsInt(N) and N >= 1 ) then
    Error("usage: InverseRationalizedHaarTransform( <positive-int> )");
  fi;
  L := LogInt(N, 2);
  if not 2^L = N then
    Error("<N> must be a power of 2");
  fi;

  R := RationalizedHaarTransform( N );
  for k in [1..N] do
    R[k] := R[k] / Sum(List(R[k], x -> x^2));
  od;
  return TransposedMat( R );
end;


# Abbreviations
# =============

#F Abbreviations
#F -------------
#F
#F   DFT       := DiscreteFourierTransform; 
#F   InvDFT    := InverseDiscreteFourierTransform;
#F   DHT       := DiscreteHartleyTransform;
#F   InvDHT    := InverseDiscreteHartleyTransform;
#F   DCT       := DiscreteCosineTransform;
#F   InvDCT    := InverseDiscreteCosineTransform;
#F   DCT_IV    := DiscreteCosineTransformIV;
#F   InvDCT_IV := InverseDiscreteCosineTransformIV;
#F   DCT_I     := DiscreteCosineTransformI;
#F   InvDCT_I  := InverseDiscreteCosineTransformI;
#F   WHT       := WalshHadamardTransform;
#F   InvWHT    := InverseWalshHadamardTransform;
#F   ST        := SlantTransform;
#F   InvST     := InverseSlantTransform;
#F   HT        := HaarTransform;
#F   InvHT     := InverseHaarTransform;
#F   RHT       := RationalizedHaarTransform;
#F   InvRHT    := InverseRationalizedHaarTransform;
#F

DFT       := DiscreteFourierTransform; 
InvDFT    := InverseDiscreteFourierTransform;
DHT       := DiscreteHartleyTransform;
InvDHT    := InverseDiscreteHartleyTransform;
DCT       := DiscreteCosineTransform;
InvDCT    := InverseDiscreteCosineTransform;
DCT_IV    := DiscreteCosineTransformIV;
InvDCT_IV := InverseDiscreteCosineTransformIV;
DCT_I     := DiscreteCosineTransformI;
InvDCT_I  := InverseDiscreteCosineTransformI;
DST       := DiscreteSineTransform;
DST_IV    := DiscreteSineTransformIV;
DST_I     := DiscreteSineTransformI;
WHT       := WalshHadamardTransform;
InvWHT    := InverseWalshHadamardTransform;
ST        := SlantTransform;
InvST     := InverseSlantTransform;
HT        := HaarTransform;
InvHT     := InverseHaarTransform;
RHT       := RationalizedHaarTransform;
InvRHT    := InverseRationalizedHaarTransform;

# weitere:
#   andere Hadamard
#   Fermat
#   Mersenne
#   Rader
