# Term-Algebra for Structured Matrices
# MP, SE 27.05.97 - GAPv3.4
# SE,    22.07.98 - IdentityAMat added

#F The category AMat of structured matrices
#F ========================================
#F
#F The category AMat contains elements representing matrices
#F (rectangular matrices over GAP-fields) in a structured way.
#F In fact AMat is a recursive term-algebra representing matrices. 
#F Basic building blocks are literal matrices (Lists of Lists), 
#F monomial matrices and permutation matrices (Permutations). Important 
#F compositions are products, direct sums and tensor products.
#F
#F We define the representation of AMat recursively in BNF as
#F the disjoint union of the following cases
#F
#F <AMat> ::= 
#F ; atomic cases
#F     <perm>                                ; "perm" (invertible)
#F   | <mon>                                 ; "mon"  (invertible)
#F   | <mat>                                 ; "mat"
#F 
#F ; composed cases
#F   | <scalar> * <AMat>                     ; "scalarMultiple"
#F   | <AMat> * .. * <AMat>                  ; "product"
#F   | <AMat> ^ <int>                        ; "power"
#F   | <AMat> ^ <AMat>                       ; "conjugate"
#F   | DirectSum(<AMat>, .., <AMat>)         ; "directSum"
#F   | TensorProduct(<AMat>, .., <AMat>)     ; "tensorProduct"
#F   | GaloisConjugate(<AMat>, <aut>)        ; "galoisConjugate".
#F
#F An AMat is a GAP-Rec A with the following mandatory 
#F fields common to all types of AMat
#F
#F   isAMat       := true, identifies AMats
#F   operations   := AMatOps, GAP operations record for A
#F   type         : a string identifying the type of A
#F   dimensions   : size of the matrix represented (= [rows, columns])
#F   char         : the characteristic of the field
#F
#F The following fields are mandatory to special types
#F
#F A.type = "perm"
#F   A.element      : defining permutation
#F
#F A.type = "mon"
#F   A.element      : defining mon-object
#F
#F A.type = "mat"
#F   A.element      : defining matrix
#F   A.isDFT        : optional field indicating that A is a DFT
#F   A.isSOR        : optional field indicating that A is a SOR
#F
#F A.type = "scalarMultiple"
#F   A.element      : the AMat multiplied
#F   A.scalar       : the scalar
#F
#F A.type = "product"
#F   A.factors      : List of AMats of suitable dimensions 
#F                    and characteristic
#F A.type = "power"
#F   A.element      : the square AMat to be raised to A.exponent
#F   A.exponent     : the exponent (an integer)
#F
#F A.type = "conjugate"
#F   A.element      : the square AMat to be conjugated
#F   A.conjugation  : the conjugating invertible AMat
#F
#F A.type = "directSum"
#F   A.summands     : List of AMats of the same characteristic
#F
#F A.type = "tensorProduct"
#F   A.factors      : List of AMats of the same characteristic
#F
#F A.type = "galoisConjugate"
#F   A.element      : the AMat to be Galois conjugated
#F   A.galoisAut    : the Galois automorphism
#F
#F Optional fields for all types of AMats
#F
#F   A.name         : a string for the name
#F   A.isInvertible : flag to indicate, that A is invertible
#F   A.inverse      : an AMat representing A^-1
#F   A.determinant  : the determinant of A
#F   A.trace        : the trace of A
#F   A.rank         : the rank of A
#F   A.isPermMat    : a flag indicating that A is A permutation matrix
#F   A.perm         : the perm represented by A
#F   A.isMonMat     : a flag indicating that A is A monomial matrix
#F   A.mon          : the mon represented by A
#F   A.mat          : the mat represented by A
#F   A.isSimplified : a flag indicating, that A is simplified
#F                    by the function SimplifyAMat
#F   

# Internal Functions for AMats
# ----------------------------
#

AMatOps :=
  OperationsRecord("AMatOps");

# AMatOps.CheckChar( <char/field> )
#   tests whether argument is a characteristic or 
#   a field and returns the characteristic

AMatOps.CheckChar := function ( x )
  if x = 0 then
    return 0;
  elif IsInt(x) and x > 0 and IsPrime(x) then
    return x;
  elif IsField(x) then
    return x.char;
  else
    Error("<x> must be 0, positive prime, or field");
  fi;
end;

# AMatOps.IsDegree( <degree>, <perm> )
#   tests whether <degree> is a degree consistent with <perm>

AMatOps.IsDegree := function ( p, n )
  if not IsPerm(p) then
    Error("usage: AMatOps.IsDegree( <degree>, <perm> )");
  fi;
  return 
    IsInt(n) and 
    n > 0 and 
    ( p = ( ) or n >= LargestMovedPointPerm(p) );
end;

# AMatOps.IsDimensions( <rc> )
#   tests whether <rc> is a pair of positive integers

AMatOps.IsDimensions := function ( rc )
  return 
    IsList(rc) and 
    Length(rc) = 2 and 
    ForAll(rc, x -> IsInt(x) and x > 0);
end;

# AMatOps.OneNC( <char> )
#   returns the one in the field of characteristic <char>.
#   The argument is not checked

AMatOps.OneNC := function ( char )
  if char = 0 then
    return 1;
  else
    return Z(char)^0;
  fi;
end;

# AMatOps.ZeroNC( <char> )
#   returns the zero in the field of characteristic <char>.
#   The argument is not checked

AMatOps.ZeroNC := function ( char )
  if char = 0 then
    return 0;
  else
    return 0*Z(char);
  fi;
end;

# AMatOps.IsSquare( <amat> )
#   tests if <amat> is a square matrix

if not IsBound(IsAMat) then
  IsAMat := "defined below";
fi;

AMatOps.IsSquare := function ( A )
  if not IsAMat(A) then
    Error("usage: AMatOps.IsSquare( <amat> )");
  fi;
  return A.dimensions[1] = A.dimensions[2];
end;

# AMatOps.IsGaloisAut( <int/field-aut> )
#   tests if <int/field-aut> defines a galois automorphism

AMatOps.IsGaloisAut := function ( aut )
  return
    IsInt(aut) or
    ( IsGeneralMapping(aut) and 
      IsField(aut.source) and 
      aut.range = aut.source and 
      IsAutomorphism(aut) 
    );
end;

# AMatOps.GaloisConjugation( <fieldelt>, <galoisAut> )
#   computes <fieldelt> ^ <galoisAut>.

AMatOps.GaloisConjugation := function ( x, galoisAut )
  if IsInt(galoisAut) and IsCyc(x) then
    return GaloisCyc(x, galoisAut);
  elif IsInt(galoisAut) and IsFFE(x) then
    return x^(CharFFE(x)^galoisAut);
  else
    if not x in galoisAut.source then
      Error("<x> not in <galoisAut>.source");
    fi;
    return x ^ galoisAut;
  fi;
end;


#F Fundamental Constructors and Tests for AMats
#F --------------------------------------------
#F

#F IsAMat( <obj> ) 
#F   tests whether <obj> is an AMat. The particular case of AMat
#F   is provided by A.type and the additional fields are available
#F   directly from the record.
#F

IsAMat := function ( A )
  return IsRec(A) and IsBound(A.isAMat) and A.isAMat;
end;

if not IsBound(MatAMat) then
  MatAMat := "defined below";
fi;

#F IsInvertibleMat( <amat> )
#F   tests if <amat> is invertible and sets the field
#F   <amat>.isInvertible.
#F

IsInvertibleMat := function ( A )

  if not IsAMat(A) then
    Error("usage: IsInvertibleMat( <amat> )");
  fi;
  
  # check if information is easy available
  if IsBound(A.isInvertible) then
    return A.isInvertible;
  fi;
  if IsBound(A.inverse) then
    A.isInvertible := true;
    return A.isInvertible;
  fi;
  if IsBound(A.determinant) then
    A.isInvertible := A.determinant <> 0*A.determinant;
    return A.isInvertible;
  fi;
  if IsBound(A.rank) then
    A.isInvertible := 
      AMatOps.IsSquare(A) and A.rank = A.dimensions[1];
    return A.isInvertible;
  fi;
  if not AMatOps.IsSquare(A) then
    A.isInvertible := false;
    return A.isInvertible;
  fi;

  # dispatch on the types
  if A.type = "perm" then 

    A.isInvertible := true;

  elif A.type = "mon" then

    A.isInvertible := true;

  elif A.type = "mat" then

    A.isInvertible := 
      AMatOps.IsSquare(A) and 
      RankMat(A.element) = A.dimensions[1];

  elif A.type = "scalarMultiple" then

    A.isInvertible := 
      A.scalar <> 0*A.scalar and
      IsInvertibleMat(A.element);

  elif A.type = "product" then

    if ForAll(A.factors, s -> AMatOps.IsSquare(s)) then
      A.isInvertible := ForAll(A.factors, IsInvertibleMat);
    else
      A.isInvertible := RankMat(MatAMat(A)) = A.dimensions[1];      
    fi;

  elif A.type = "power" then

    A.isInvertible := IsInvertibleMat(A.element);

  elif A.type = "conjugate" then

    A.isInvertible := IsInvertibleMat(A.element);

  elif A.type = "directSum" then

    # if any summand is not square, then the direct sum
    # is not invertible
    A.isInvertible := 
      ForAll(A.summands, s -> AMatOps.IsSquare(s)) and
      ForAll(A.summands, IsInvertibleMat);

  elif A.type = "tensorProduct" then

    A.isInvertible := ForAll(A.factors, IsInvertibleMat);   
 
  elif A.type = "galoisConjugate" then

    A.isInvertible := IsInvertibleMat(A.element);

  else

    Error("unknown type of amat <A>");

  fi;

  return A.isInvertible;
end;


#F IsIdentityMat( <amat> )
#F   checks whether <amat> is a square IdentityMat.
#F

if not IsBound(PermAMat) then
  PermAMat  := "defined below";
  IsPermMat := "defined below";
fi;

IsIdentityMat := function ( A )
  return IsPermMat(A) and PermAMat(A) = ();
end;

#F IdentityPermAMat( <size>       [, <char/field> ] )
#F IdentityMonAMat(  <size>       [, <char/field> ] )
#F IdentityMatAMat(  <size>       [, <char/field> ] )
#F IdentityMatAMat(  <dimensions> [, <char/field> ] )
#F IdentityAMat(     <dimensions> [, <char/field> ] )
#F   construct an "perm"/"mon"/"mat"-amat representing 
#F   the identity matrix of size <size>. In the case "mat"
#F   the matrix may be rectangle containing at the position
#F   (i, j) one, if i = j and zero else. The argument <size>
#F   is a positive integer, the argument <dimensions> is a
#F   pair of positive integers for nr. of rows, columns. 
#F   The function IdentityAMat either constructs a PermAMat
#F   if the matrix is square or a MatAMat if the matrix is
#F   rectangular. Use this function if you do not know whether
#F   the matrix is square and you do not care about the type. 
#F

if not IsBound(AMatPerm) then
  AMatPerm := "defined below";
  AMatMon  := "defined below";
  AMatMat  := "defined below";
fi;

IdentityPermAMat := function ( arg )
  if 
    Length(arg) = 1 and 
    AMatOps.IsDimensions([arg[1], arg[1]])
  then
    return AMatPerm( (), arg[1] );
  elif 
    Length(arg) = 2 and 
    AMatOps.IsDimensions([arg[1], arg[1]]) 
  then
    return AMatPerm( (), arg[1], AMatOps.CheckChar(arg[2]) );
  else
    Error("usage: IdentityPermAMat( <size> [, <char/field> ] )");
  fi;
end;

IdentityMonAMat := function ( arg )
  if 
    Length(arg) = 1 and 
    AMatOps.IsDimensions([arg[1], arg[1]])
  then
    return AMatMon( Mon( (), arg[1] ) );
  elif 
    Length(arg) = 2 and 
    AMatOps.IsDimensions([arg[1], arg[1]]) 
  then
    return AMatMon( Mon( (), arg[1], AMatOps.CheckChar(arg[2]) ) );
  else
    Error("usage: IdentityMonAMat( <size> [, <char/field> ] )");
  fi;
end;

IdentityMatAMat := function ( arg )
  local M, i, char;

  if 
    Length(arg) = 1 and 
    AMatOps.IsDimensions([arg[1], arg[1]])
  then
    return AMatMat( IdentityMat(arg[1]) );
  elif
    Length(arg) = 1 and 
    AMatOps.IsDimensions(arg[1]) 
  then 
    M := 
      List(
        [1..arg[1][1]], 
        i -> List([1..arg[1][2]], j -> 0)
      );
    for i in [1..Minimum([arg[1][1], arg[1][2]])] do
      M[i][i] := 1;
    od;

    return AMatMat(M);
  elif 
    Length(arg) = 2 and 
    AMatOps.IsDimensions([arg[1], arg[1]]) 
  then
    char := AMatOps.CheckChar(arg[2]);
    if char = 0 then
      return AMatMat( IdentityMat(arg[1]) );
    else
      return AMatMat( IdentityMat(arg[1], GF(char)) );
    fi;
  elif
    Length(arg) = 2 and 
    AMatOps.IsDimensions(arg[1]) 
  then 
    char := AMatOps.CheckChar(arg[2]);
    M := 
      List(
        [1..arg[1][1]], 
        i -> List([1..arg[1][2]], j -> AMatOps.ZeroNC(char))
      );
    for i in [1..Minimum([arg[1][1], arg[1][2]])] do
      M[i][i] := AMatOps.OneNC(char);
    od;

    return AMatMat(M);
  else
    Error("usage: IdentityMatAMat( <size> [, <char/field> ] )");
  fi;
end;

IdentityAMat := function ( arg )
  local dims, char;

  if 
    Length(arg) = 1 and
    AMatOps.IsDimensions(arg[1])
  then

    # IdentityAMat( <dimensions> )
    dims := arg[1];
    char := 0;

  elif 
    Length(arg) = 2 and
    AMatOps.IsDimensions(arg[1])
  then

    # IdentityAMat( <dimensions>, <char/field> )
    dims := arg[1];
    char := AMatOps.CheckChar(arg[2]);

  else
    Error("usage: IdentityAMat( <dimensions> [, <char/field> ] )");
  fi;

  if dims[1] = dims[2] then
    return IdentityPermAMat(dims[1], char);
  else
    return IdentityMatAMat(dims, char);
  fi;
end;


#F AllOneAMat( <size> [, <char/field> ] )
#F AllOneAMat( <dimensions> [, <char/field> ] )
#F   constructs the all-one matrix of type "mat" of given 
#F   size and characteristic. The default characteristic is 0.
#F

AllOneAMat := function ( arg )
  local dims, char, M;

  if Length(arg) in [1, 2] then
    if not IsList(arg[1]) then
      dims := [arg[1], arg[1]];
    else
      dims := arg[1];
    fi;
    if Length(arg) = 1 then
      char := 0;
    else
      char := AMatOps.CheckChar(arg[2]);
    fi;
  else
    Error(
      "usage: ",
      "  AllOneAMat( <size> [, <char/field> ] )",
      "  AllOneAMat( <dimensions> [, <char/field> ] )"
    );
  fi;

  # check arguments
  if not AMatOps.IsDimensions(dims) then
    Error(
      "usage: ",
      "  AllOneAMat( <size> [, <char/field> ] )",
      "  AllOneAMat( <dimensions> [, <char/field> ] )"
    );
  fi;

  # construct matrix
  M := 
    List(
      [1..dims[1]], 
      i -> List([1..dims[2]], j -> AMatOps.OneNC(char))
    );
  return
    rec(
      isAMat       := true,
      operations   := AMatOps,
      type         := "mat",
      dimensions   := dims,
      char         := char,
      element      := M,
 
      rank         := 1
    );      

end;


#F NullAMat( <size> [, <char/field> ] )
#F NullAMat( <dimensions> [, <char/field> ] )
#F   construct the all-zero matrix of type "mat" of given
#F   size and characteristic. The default characteristic is 0.
#F

NullAMat := function ( arg )
  local dims, char, M;

  if Length(arg) in [1, 2] then
    if not IsList(arg[1]) then
      dims := [arg[1], arg[1]];
    else
      dims := arg[1];
    fi;
    if Length(arg) = 1 then
      char := 0;
    else
      char := AMatOps.CheckChar(arg[2]);
    fi;
  else
    Error(
      "usage: ",
      "  AllOneAMat( <size> [, <char/field> ] )",
      "  AllOneAMat( <dimensions> [, <char/field> ] )"
    );
  fi;

  # check arguments
  if not AMatOps.IsDimensions(dims) then
    Error(
      "usage: ",
      "  AllOneAMat( <size> [, <char/field> ] )",
      "  AllOneAMat( <dimensions> [, <char/field> ] )"
    );
  fi;

  # construct matrix
  M := 
    List(
      [1..dims[1]], 
      i -> List([1..dims[2]], j -> AMatOps.ZeroNC(char))
    );
  return
    rec(
      isAMat       := true,
      operations   := AMatOps,
      type         := "mat",
      dimensions   := dims,
      char         := char,
      element      := M,
 
      trace        := AMatOps.ZeroNC(char),
      rank         := 0
    );      

end;

#F DiagonalAMat( <list> )
#F   construct the amat having <list> as diagonal. If all
#F   elements of <list> are <> 0, then the resulting amat 
#F   is of type "mon", else of type "directSum".
#F   Note that the elements of list must lie in a common
#F   field.
#F

if not IsBound(DirectSumAMat) then
  DirectSumAMat := "defined below";
fi;

DiagonalAMat := function ( diag )
  local char, degree;

  if not IsList(diag)then
    Error("usage: DiagonalAMat( <list> )");
  fi;
  char := DefaultField(diag).char;

  if ForAll(diag, x -> x <> 0*x) then

    # construct amat of type "mon"
    degree := Length(diag);
    return
      rec(
        isAMat       := true,
        operations   := AMatOps,
        type         := "mon",
        dimensions   := [degree, degree],
        char         := char,
        element      := Mon(diag)
      );      
  else
    
    # construct amat of type "directSum"
    return DirectSumAMat(List(diag, x -> AMatMat([[x]])));     
  fi;
end;


#F DFTAMat( <size> [, <char/field>] )
#F   constructs a "mat"-amat representing a DFT of 
#F   size <size>. The default characteristic is 0.
#F   If <char> is given, then a DFT is constructed
#F   iff <size> and <char> are coprime, if <field> is
#F   given, then a DFT is constructed iff <size>
#F   divides Size(<field>) - 1.
#F   The field <amat>.isDFT is set.
#F

DFTAMat := function ( arg )
  local size, char, A;

  if Length(arg) = 1 then
    size := arg[1];
    char := 0;
  elif Length(arg) = 2 then
    size := arg[1];
    char := AMatOps.CheckChar(arg[2]);

    # in the case char <> 0 and a given field arg[2]
    # check if the field is large enough
    if 
      char <> 0 and 
      IsField(arg[2]) and 
      (Size(arg[2]) - 1) mod size <> 0 
    then
      Error("the field <arg[2]> contains no <size>th root of unity");
    fi;
  else
    Error("usage: DFTAMat( <size> [, <char/field>] )");
  fi;

  A          := AMatMat(DFT(size, char), "invertible");
  A.isDFT    := true;
  A.isMonMat := (size < 2);
  
  return A;
end;


#F SORAMat( <size> [, <char/field>] )
#F   constructs a "mat"-amat representing an SOR of
#F   size <size>. SOR(n) (Split One Rep) is defined 
#F   as the following (n x n) - matrix
#F
#F               _                  _
#F              |1  1  1  1  .  .    |
#F              |1 -1  0  0          |
#F              |1  0 -1  0          |
#F     SOR(n) = |1  0  0 -1          |
#F              |.           .       |
#F              |.              .    |
#F              |                    |
#F   
#F   The SOR(n) is the sparsest matrix that splits off
#F   the onerep contained in a permrep. The number of
#F   entries is 3n - 2.        
#F   The default characteristic is 0. The field 
#F   <amat>.isSOR is set.
#F

SORAMat := function ( arg )
  local size, char, M, i;

  if Length(arg) = 1 then
    size := arg[1];
    char := 0;
  elif Length(arg) = 2 then
    size := arg[1];
    char := AMatOps.CheckChar(arg[2]);
  else
    Error("usage: SORAMat( <size> [, <char/field>] )");
  fi;

  if not( IsInt(size) and size > 0 ) then
    Error("<size> must be a positive integer");
  fi;
  if size in [1, 2] then
    Error("<size> is smaller than 3");
  fi;

  M := MatAMat(-IdentityPermAMat(size, char));
  for i in [1..size] do
    M[1][i] := M[1][i] + AMatOps.OneNC(char);
    M[i][1] := M[i][1] + AMatOps.OneNC(char);
  od;

  M          := AMatMat(M, "invertible");
  M.isSOR    := true;
  M.isMonMat := false;

  return M;
end;



#F AMatPerm( <perm>, <degree> [, <char/field> ] )
#F   constructs a "perm"-amat from <perm> of given degree 
#F   <degree> and characteristic <char/field>.
#F

AMatPerm := function ( arg )
  local perm, degree, char;

  if Length(arg) = 2 then
    perm   := arg[1];
    degree := arg[2];
    char   := 0;
  elif Length(arg) = 3 then
    perm   := arg[1];
    degree := arg[2];
    char   := AMatOps.CheckChar(arg[3]);
  else
    Error("usage: AMatPerm( <perm>, <degree> [, <char/field> ] )");
  fi;
  if perm = false then
    return false;
  fi;

  if not(
    IsPerm(perm) and
    AMatOps.IsDegree(perm, degree) 
  ) then
    Error("usage: AMatPerm( <perm>, <degree> [, <char/field> ] )");
  fi;

  return
    rec(
      isAMat       := true,
      operations   := AMatOps,
      type         := "perm",
      dimensions   := [degree, degree],
      char         := char,
      element      := perm
    );      

end;

#F AMatMon( <mon> )
#F   constructs a "mon"-amat from <mon>.
#F

AMatMon := function ( mon )
  local degree;

  if mon = false then 
    return false;
  fi;
  if not IsMon(mon) then
    Error("usage: AMatMon( <mon> )");
  fi;

  degree := Length(mon.diag);
  return
    rec(
      isAMat       := true,
      operations   := AMatOps,
      type         := "mon",
      dimensions   := [degree, degree],
      char         := DefaultField(mon.diag[1]).char,
      element      := mon
    );      
  
end;

#F AMatMat( <mat> [, <hint>] )
#F   constructs a "mat"-amat from <mat>. The <hint>
#F   "invertible" may be given.
#F

AMatMat := function ( arg )
  local mat, hint, A;

  if Length(arg) in [1,2] then
    mat := arg[1];
    if Length(arg) = 2 then
      hint := arg[2];
    else
      hint := "no hint";
    fi;
  else
    Error("usage: AMatMat( <mat> )");
  fi;
  if mat = false then
    return mat;
  fi;

  # check arguments
  if not(
    IsMat(mat) and
    hint in ["no hint", "invertible"]
  ) then
    Error("usage: AMatMat( <mat> )");
  fi;

  A := 
    rec(
      isAMat       := true,
      operations   := AMatOps,
      type         := "mat",
      dimensions   := DimensionsMat(mat),
      char         := DefaultField(mat[1][1]).char,
      element      := mat
    ); 

  if hint = "invertible" then
     A.isInvertible := true;
  fi;
     
  return A;  
end;


#F Structural Symbolic Constructors for AMats
#F ------------------------------------------
#F
#F ScalarMultipleAMat( <scalar>, <amat> )
#F   forms the amat of type "scalarMultiple" representing
#F   the scalar multiple of <scalar> and <amat>, which must
#F   have the same characteristic. In addition, the integers
#F   operate on the finite fields. Hence, -<amat> is possible.
#F

ScalarMultipleAMat := function ( scalar, A )
  if not scalar in FieldElements then
    Error("<scalar> must be a field element");
  fi;
  if not IsAMat(A) then
    Error("<A> must be an AMat");
  fi;
  if IsInt(scalar) then
    scalar := scalar * AMatOps.OneNC(A.char);
  fi;
  if not DefaultField(scalar).char = A.char then
    Error("<scalar> and <A> must have common char");
  fi;

  return
    rec(
      isAMat       := true,
      operations   := AMatOps,
      type         := "scalarMultiple",
      dimensions   := A.dimensions,
      char         := A.char,
      element      := A,
      scalar       := scalar
    );
end;


#F <amat> * <amat>
#F <scalar> * <amat>
#F   forms the product resp. the scalarMultiple. In the first
#F   case the sizes must be compatible, in both cases the
#F   characteristic must be compatible. In addition, the integers
#F   operate on the finite fields. Hence, -<AMat> is admissible.

AMatOps.\* := function ( A, B )
  local AB;

  # catch the case "scalarMultiple"
  if A in FieldElements and IsAMat(B) then
    return ScalarMultipleAMat(A, B);
  elif B in FieldElements and IsAMat(A) then
    return ScalarMultipleAMat(B, A);
  fi;
  if not ( IsAMat(A) and IsAMat(B) ) then
    Error("sorry, do not know how to compute <A> * <B>");
  fi;

  if not( 
    A.dimensions[2] = B.dimensions[1] and
    A.char = B.char
  ) then
    Error("amats <A> and <B> must be compatible in size and char");
  fi;

  # flatten products by associativity
  AB := [ ];
  if A.type = "product" then
    Append(AB, A.factors);
  else
    Add(AB, A);
  fi;
  if B.type = "product" then
    Append(AB, B.factors);
  else
    Add(AB, B);
  fi;

  return
    rec(
      isAMat       := true,
      operations   := AMatOps,
      type         := "product",
      dimensions   := [A.dimensions[1], B.dimensions[2]],
      char         := A.char,
      factors      := AB
    );
end;

#F <amat> / <amat>
#F   forms the quotient of the amats. Characteristic and
#F   size have to be compatible, the second amat must be 
#F   square and invertible.
#F

AMatOps.\/ := function ( A, B )
  return A * B^-1;
end;


#F PowerAMat( <amat>, <int> [, <hint> ] )
#F   forms the amat <amat> ^ <int> of type "power". If <int>
#F   is negative then <amat> is checked for invertibility if
#F   not the hint "invertible" is supplied.
#F

PowerAMat := function ( arg )
  local M, n, hint, A;
  if Length(arg) in [2, 3] then
    M := arg[1];
    n := arg[2];
    if Length(arg) = 2 then
      hint := "no hint";
    else
      hint := arg[3];
    fi;
  else
    Error("usage: PowerAMat( <amat1>, <int> [, <hint> ] )");
  fi;

  # check arguments
  if not(
    IsAMat(M) and
    IsInt(n) and
    hint in ["no hint", "invertible"] 
  ) then
    Error("usage: PowerAMat( <amat1>, <amat2> [, <hint> ] )");
  fi;
  if not AMatOps.IsSquare(M) then
    Error("<M> must be a square matrix");
  fi;
  if 
    n < 0 and 
    hint = "no hint" and 
    not IsInvertibleMat(M) 
  then
    Error("<M> must be invertible");
  fi;

  # set the field M.isInvertible if the hint is given
  if hint = "invertible" then
    M.isInvertible := true;
  fi;

  # construct result
  return
    rec(
      isAMat       := true,
      operations   := AMatOps,
      type         := "power",
      dimensions   := M.dimensions,
      char         := M.char,
      element      := M,
      exponent     := n
    );
end;


#F ConjugateAMat( <amat1>, <amat2> [, <hint> ] )
#F   forms the conjugate <amat1> ^ <amat2> of type "conjugate" . 
#F   The amat <amat2> is checked if it is invertible if not the 
#F   hint "invertible" is supplied
#F

ConjugateAMat := function ( arg )
  local M, T, hint;

  if Length(arg) in [2, 3] then
    M := arg[1];
    T := arg[2];
    if Length(arg) = 2 then
      hint := "no hint";
    else
      hint := arg[3];
    fi;
  else
    Error("usage: ConjugateAMat( <amat1>, <amat2> [, <hint> ] )");
  fi;

  # check arguments
  if not(
    IsAMat(M) and
    IsAMat(T) and
    hint in ["no hint", "invertible"] 
  ) then
    Error("usage: ConjugateAMat( <amat1>, <amat2> [, <hint> ] )");
  fi;
  if not( 
    AMatOps.IsSquare(M) and 
    M.dimensions = T.dimensions
  ) then
    Error("amats <M>, <T> must have suitable sizes");
  fi;
  if not M.char = T.char then
    Error("amats <M>, <T> must have suitable characteristics");
  fi;
  if hint = "no hint" and not IsInvertibleMat(T) then
    Error("conjugating amat <T> must be invertible");
  fi;

  return
    rec(
      isAMat       := true,
      operations   := AMatOps,
      type         := "conjugate",
      dimensions   := M.dimensions,
      char         := M.char,
      element      := M,
      conjugation  := T
    );
end;


#F <amat> ^ <int>
#F <amat> ^ <amat>
#F   forms the power of an amat resp. the conjugate of 
#F   an amat by another amat. The amats have to be square.
#F

if not IsBound(IsARep) then
  IsARep :=        "defined in arep.g";
  ConjugateARep := "defined in arep.g";
fi;

AMatOps.\^ := function ( A, B )

  # distinguish the cases "power" and "conjugate"
  if IsAMat(A) and IsInt(B) then
    return PowerAMat(A, B);
  elif IsAMat(A) and IsAMat(B) then
    return ConjugateAMat(A, B);
  elif IsFunc(IsARep) and IsARep(A) and IsAMat(B) then
    return ConjugateARep(A, B);
  else
    Error("sorry, can not compute <A> ^ <B>");
  fi;
end;


#F DirectSumAMat( <amat1>, .., <amatN> ) ; N >= 1
#F DirectSumAMat( <list-of-amat> )
#F   forms the direct sum of the amats given. Note that the
#F   amats does not have to be square, but must be of common 
#F   characteristic.
#F

DirectSumAMat := function ( arg )
  if IsList(arg[1]) then
    arg := arg[1];
  fi;

  # check arguments
  if not ForAll(arg, IsAMat) then
    Error(
      "usage: ",
      "  DirectSumAMat( <amat1>, .., <amatN> ) ; N >= 1",
      "  DirectSumAMat( <list-of-amat> )"
    );
  fi;
  if not ForAll(arg, x -> x.char = arg[1].char) then
    Error("amats must have common characteristic");
  fi;

  return
    rec(
      isAMat       := true,
      operations   := AMatOps,
      type         := "directSum",
      dimensions   := Sum(List(arg, x -> x.dimensions)),
      char         := arg[1].char,
      summands     := arg
    );
end;


#F TensorProductAMat( <amat1>, .., <amatN> ) ; N >= 1
#F TensorProductAMat( <list-of-amat> )
#F   forms the tensor (or kronecker) product of the 
#F   amats given. Note that the amats does not have to be 
#F   square, but must be of common characteristic.
#F

TensorProductAMat := function ( arg )
  local dims;

  if IsList(arg[1]) then
    arg := arg[1];
  fi;

  # check arguments
  if not ForAll(arg, IsAMat) then
    Error(
      "usage: ",
      "  TensorProductAMat( <amat1>, .., <amatN> ) ; N >= 1",
      "  TensorProductAMat( <list-of-amat> )"
    );
  fi;
  if not ForAll(arg, x -> x.char = arg[1].char) then
    Error("amats must have common characteristic");
  fi;

  dims := 
    [Product(List(arg, x -> x.dimensions[1])),
     Product(List(arg, x -> x.dimensions[2]))
    ];
  return
    rec(
      isAMat       := true,
      operations   := AMatOps,
      type         := "tensorProduct",
      dimensions   := dims,
      char         := arg[1].char,
      factors      := arg
    );
end;


#F GaloisConjugateAMat( <amat>, <gal-aut/int> )
#F   forms the Galois-conjugate of A under the field automorphism aut.
#F   The automorphism may be specified by an integer k: In case of a
#F   finite field aut = FrobeniusAutomorphism(A.baseField)^k. In case 
#F   of a number field (GAP: a subfield of a cyclotomic field) it is
#F   aut = NFAutomorphism(A.baseField, k) = x -> GaloisCyc(x, k).
#F

GaloisConjugateAMat := function ( A, aut )
  local ArX;

  # check arguments
  if not (
    IsAMat(A) and
    AMatOps.IsGaloisAut(aut) 
  ) then
    Error("usage: GaloisConjugateAMat( <amat>, <gal-aut/int> )");
  fi;

  return
    rec(
      isAMat       := true,
      operations   := AMatOps,
      type         := "galoisConjugate",
      dimensions   := A.dimensions,
      char         := A.char,
      element      := A,
      galoisAut    := aut
    );
end;


#F Comparison of AMats
#F -------------------
#F

#F AMatOps.\=( <amat1>, <amat2> )
#F   tests if <amat1> and <amat2> represent equal matrices
#F   in the mathematical sense which means that the entries 
#F   of the matrix represented are equal.
#F

AMatOps.\= := function ( A1, A2 )

  if not IsAMat(A1) and not IsAMat(A2) then
    Error("sorry, don't know how to compare <A1> = <A2>");
  fi;
  if IsAMat(A1) and not IsAMat(A2) then
    return false;
  fi;
  if not IsAMat(A1) and IsAMat(A2) then
    return false;
  fi;

  # cheap checks
  if A1.dimensions <> A2.dimensions then
    return false;
  fi;
  if A1.char <> A2.char then
    return false;
  fi;
  
  # expensive checks
  if IsPermMat(A1) and not IsPermMat(A2) then
    return false;
  elif not IsPermMat(A1) and IsPermMat(A2) then
    return false;
  elif IsPermMat(A1) and IsPermMat(A2) then
    return PermAMat(A1) = PermAMat(A2);
  fi;

  if IsMonMat(A1) and not IsMonMat(A2) then
    return false;
  elif not IsMonMat(A1) and IsMonMat(A2) then
    return false;
  elif IsMonMat(A1) and IsMonMat(A2) then
    return MonAMat(A1) = MonAMat(A2);
  fi;

  return MatAMat(A1) = MatAMat(A2);
end;


#F Pretty Printing of AMats
#F ------------------------
#F
#F AMatOps.Print( <amat> [, <indent> , <indentStep> , <bp> ] )
#F

# Syntactic operators
#   operator  binding-power(BP)  associative
#     none          0
#      +           10                yes
#    scalar*       15                no
#      *           20                yes 
#      ^           30                no
#    f(..)        100                no

AMatOps.PrintNC := function ( A, indent, indentStep, bp)
  local newline, i, bpScalar;

  newline := function ( )
    local i;

    Print("\n");
    for i in [1..indent] do
      Print(" ");
    od;
  end;

  bpScalar := function ( x )
    if IsRat(x) or IsFFE(x) then
      return 20;
    elif 
      IsCyc(x) and 
      Number(CoeffsCyc(x, NofCyc(x)), xi -> xi <> 0) <= 1
    then
      return 20;
    else
      return 0;
    fi;
  end;

  if IsBound(A.name) then
    Print(A.name);
    return;
  fi;

  if A.type = "perm" then

    if A.element = () then

      # use IdentityPermAMat
      Print("IdentityPermAMat(");
      Print(A.dimensions[1]);
      if A.char <> 0 then
        Print(", GF(");
        Print(A.char);
        Print(")");
      fi;
      Print(")");

    else

      # use AMatPerm
      Print("AMatPerm(");
      Print(A.element);
      Print(", ");
      Print(A.dimensions[1]);
      if A.char <> 0 then
        Print(", GF(");
        Print(A.char);
        Print(")");
      fi;
      Print(")");

    fi;
    return;

  elif A.type = "mon" then

    if A.element.perm = () then

      if ForAll(A.element.diag, x -> x = x^0) then

        # use IdentityMonAMat
        Print("IdentityMonAMat(");
        Print(Length(A.element.diag));
        if A.char <> 0 then
          Print(", GF(");
          Print(A.char);
          Print(")");
        fi;  
        Print(")");
  
      else
  
        # use DiagonalAMat
        Print("DiagonalAMat(");
        Print(A.element.diag);
        Print(")");

      fi;
      return;

    else

      # use AMatMon
      Print("AMatMon( ");
      MonOps.PrintNC(A.element, indent, indentStep);
      Print(" )");

    fi;
    return;

  elif A.type = "mat" then

    if 
      ForAll(
        A.element, 
        Ai -> ForAll(Ai, Aij -> Aij = 0*Aij)
      )
    then

      # use NullAMat
      Print("NullAMat(");
      if AMatOps.IsSquare(A) then
        Print(A.dimensions[1]);
      else
        Print(A.dimensions);
      fi;
      if A.char <> 0 then
        Print(", GF(");
        Print(A.char);
        Print(")");
      fi;
      Print(")");

    elif
      ForAll(
        [1..A.dimensions[1]],
        i -> 
          ForAll(
            [1..A.dimensions[2]], 
            function(j) 
              if i = j then
                return A.element[i][j] = A.element[i][j]^0;
              else
                return A.element[i][j] = 0*A.element[i][j];
              fi;
            end
          )
       )
    then

      # use IdentityMatAMat
      Print("IdentityMatAMat(");
      if A.dimensions[1] <> A.dimensions[2] then
        Print(A.dimensions);
      else 
        Print(A.dimensions[1]);
      fi;
      if A.char <> 0 then
        Print(", GF(");
        Print(A.char);
        Print(")");
      fi;  
      Print(")");
      return;      
 
    elif
      ForAll(
        A.element, 
        Ai -> ForAll(Ai, Aij -> Aij = Aij^0)
      )
    then

      # use AllOneAMat
      Print("AllOneAMat(");
      if AMatOps.IsSquare(A) then
        Print(A.dimensions[1]);
      else
        Print(A.dimensions);
      fi;
      if A.char <> 0 then
        Print(", GF(");
        Print(A.char);
        Print(")");
      fi;
      Print(")");

    elif IsBound(A.isDFT) and A.isDFT = true then

      # use DFTAMat
      Print("DFTAMat(");
      Print(A.dimensions[1]);
      if A.char <> 0 then
        Print(", ");
        Print(A.char);
      fi;
      Print(")");

    elif IsBound(A.isSOR) and A.isSOR = true then

      # use SORAMat
      Print("SORAMat(");
      Print(A.dimensions[1]);
      if A.char <> 0 then
        Print(", ");
        Print(A.char);
      fi;
      Print(")");

    else

      # use AMatMat
      Print("AMatMat(");
      indent := indent + indentStep;
        newline();
        Print(A.element);
        if IsBound(A.isInvertible) then
          Print(",");
          newline();
          Print("\"invertible\"");
        fi;
      indent := indent - indentStep;
      newline();
      Print(")");

    fi;
    return;

  elif A.type = "scalarMultiple" then

    if bp >= 15 then
      Print("(");
      for i in [2..indentStep] do
        Print(" ");
      od;
      indent := indent + indentStep;
    fi;   

    # print scalar
    if bpScalar(A.scalar) < 15 then
      Print("(");
      Print(A.scalar);
      Print(")");
    else
      Print(A.scalar);       
    fi;

    Print(" * ");
    AMatOps.PrintNC(A.element, indent, indentStep, 15);
    
    if bp >= 15 then
      indent := indent - indentStep;
      newline();
      Print(")");
    fi;

    return;

  elif A.type = "product" then

    if bp > 20 then
      Print("(");
      for i in [2..indentStep] do
        Print(" ");
      od;
      indent := indent + indentStep;
    fi;

    for i in [1..Length(A.factors)] do
      AMatOps.PrintNC(A.factors[i], indent, indentStep, 20);
      if i < Length(A.factors) then
        Print(" *");
        newline();
      fi;
    od;

    if bp > 20 then
      indent := indent - indentStep;
      newline();
      Print(")");
    fi;
  
    return;

  elif A.type = "power" then

    if bp >= 30 then
      Print("(");
      for i in [2..indentStep] do
        Print(" ");
      od;
      indent := indent + indentStep;
    fi;

    AMatOps.PrintNC(A.element, indent, indentStep, 30);
    Print(" ^ ");
    Print(A.exponent);

    if bp >= 30 then
      indent := indent - indentStep;
      newline();
      Print(")");
    fi;

    return;

  elif A.type = "conjugate" then

    Print("ConjugateAMat(");
    indent := indent + indentStep;
      newline();
      AMatOps.PrintNC(A.element, indent, indentStep, 0);
      Print(",");
      newline();
      AMatOps.PrintNC(A.conjugation, indent, indentStep, 0);
    indent := indent - indentStep;
    newline();
    Print(")");
    return;

  elif A.type = "directSum" then

    Print("DirectSumAMat(");
    indent := indent + indentStep;
      newline();
      for i in [1..Length(A.summands)] do
        AMatOps.PrintNC(A.summands[i], indent, indentStep, 0);
        if i < Length(A.summands) then
          Print(",");
          newline();
        fi;
      od;
    indent := indent - indentStep;
    newline();
    Print(")");
    return;

  elif A.type = "tensorProduct" then

    Print("TensorProductAMat(");
    indent := indent + indentStep;
      newline();
      for i in [1..Length(A.factors)] do
        AMatOps.PrintNC(A.factors[i], indent, indentStep, 0);
        if i < Length(A.factors) then
          Print(",");
          newline();
        fi;
      od;
    indent := indent - indentStep;
    newline();
    Print(")");
    return;

  elif A.type = "galoisConjugate" then

    Print("GaloisConjugateAMat(");
    indent := indent + indentStep;
      newline();
      AMatOps.PrintNC(A.element, indent, indentStep, 0);
      Print(",");
      newline();
      Print(A.galoisAut);
    indent := indent - indentStep;
    newline();
    Print(")");
    return;

  else
    Error("unknown type of amat <A>");
  fi;
end;

AMatOps.Print := function ( arg )
  local A, indent, indentStep, bp;

  if Length(arg) = 1 then
    A          := arg[1];
    indent     := 0;
    indentStep := 2;
    bp         := 0;
  elif Length(arg) = 4 then
    A          := arg[1];
    indent     := arg[2];
    indentStep := arg[3];
    bp         := arg[4];
  else
    Error(
      "usage: AMatOps.Print( <amat> [, <indent> , <indentStep> , <bp> ] )"
    );
  fi;

  # check arguments
  if not IsAMat(A) then
    Error("<A> must be an amat");
  fi;
  if not( 
    IsInt(indent) and indent >= 0 and
    IsInt(indentStep) and indentStep >= 0 and
    IsInt(bp) and bp >= 0
  ) then
    Error(
      "usage: AMatOps.Print( <amat> [, <indent> , <indentStep> , <bp> ] )"
    );
  fi;

  AMatOps.PrintNC(A, indent, indentStep, bp);
end;


#F Fundamental Operations with AMats
#F ---------------------------------
#F
#F InverseAMat( <amat> )
#F   returns an amat representing the inverse of
#F   <amat> if possible. The calculation uses the fact, 
#F   that inversion is compatible with most of the
#F   structures represented by a mat.
#F

InverseAMat := function ( A )
  local perm, i;

  if not IsAMat(A) then
    Error("usage: InverseAMat( <amat> )");
  fi;
  if IsBound(A.inverse) then
    return A.inverse;
  fi;
  if not IsInvertibleMat(A) then
    Error("<A> must be invertible");
  fi;

  # dispatch on the types
  if A.type = "perm" then 

    A.inverse := AMatPerm(A.element^-1, A.dimensions[1], A.char);

  elif A.type = "mon" then

    A.inverse := AMatMon(A.element^-1);

  elif A.type = "mat" then

    # catch DFT case
    if IsBound(A.isDFT) and A.isDFT = true then
      perm := [ 1 ];
      for i in [2..A.dimensions[1]] do
        perm[i] := A.dimensions[1] + 2 - i;
      od;
      A.inverse := 
        AMatOps.OneNC(A.char)/A.dimensions[1] *
        AMatPerm(PermList(perm), A.dimensions[1], A.char) *
        DFTAMat(A.dimensions[1], A.char);
    else
      A.inverse := AMatMat(A.element^-1, "invertible");
    fi;

  elif A.type = "scalarMultiple" then

    A.inverse := A.scalar^-1 * InverseAMat(A.element);

  elif A.type = "product" then

    A.inverse := Product(List(Reversed(A.factors), InverseAMat));

  elif A.type = "power" then

    if A.exponent < 0 then
      A.inverse := PowerAMat(A.element, -A.exponent);
    else
      A.inverse := PowerAMat(InverseAMat(A.element), A.exponent);
    fi;

  elif A.type = "conjugate" then

    A.inverse := 
      ConjugateAMat(
        InverseAMat(A.element), 
        A.conjugation, 
        "invertible"
      );

  elif A.type = "directSum" then

    A.inverse := DirectSumAMat(List(A.summands, InverseAMat));

  elif A.type = "tensorProduct" then

    A.inverse := TensorProductAMat(List(A.factors, InverseAMat));

  elif A.type = "galoisConjugate" then

    A.inverse := 
      GaloisConjugateAMat(InverseAMat(A.element), A.galoisAut);

  else
    Error("unrecognized type of amat <A>");
  fi;

  return A.inverse;
end;


#F TransposedAMat( <amat> )
#F   returns an amat representing the transposed of
#F   <amat>. The calculation uses the fact, 
#F   that transposition is compatible with most of the
#F   structures represented by a mat.
#F

TransposedAMat := function ( A )

  if not IsAMat(A) then
    Error("usage: TransposedAMat( <amat> )");
  fi;

  # dispatch on the types
  if A.type = "perm" then 

    return AMatPerm(A.element^-1, A.dimensions[1], A.char);

  elif A.type = "mon" then

    return AMatMon(TransposedMon(A.element));

  elif A.type = "mat" then

    # catch DFT case
    if IsBound(A.isDFT) and A.isDFT = true then
      return A;
    else
      return AMatMat(TransposedMat(A.element));
    fi;

  elif A.type = "scalarMultiple" then

    return A.scalar*TransposedAMat(A.element);

  elif A.type = "product" then

    return Product(List(Reversed(A.factors), TransposedAMat));

  elif A.type = "power" then

    return PowerAMat(TransposedAMat(A.element), A.exponent, "invertible");

  elif A.type = "conjugate" then

    return 
      ConjugateAMat(
        TransposedAMat(A.element),
        InverseAMat(TransposedAMat(A.conjugation)),
        "invertible"
      );

  elif A.type = "directSum" then

    return DirectSumAMat(List(A.summands, TransposedAMat));

  elif A.type = "tensorProduct" then

    return TensorProductAMat(List(A.factors, TransposedAMat));

  elif A.type = "galoisConjugate" then

    return
      GaloisConjugateAMat(TransposedAMat(A.element), A.galoisAut);

  else
    Error("unrecognized type of amat <A>");
  fi;

end;


#F DeterminantAMat( <amat> )
#F   calculates the determinant of <amat> if it is square
#F   and stores the result in <amat>.determinant
#F

DeterminantAMat := function ( A )
  local d, i;

  if not IsAMat(A) then
    Error("usage: DeterminantAMat( <amat> )");
  fi;
  if not AMatOps.IsSquare(A) then
    Error("<A> is not square");
  fi;

  # check if determinant is easy to calculate
  if IsBound(A.determinant) then
    return A.determinant;
  fi;  
  if IsBound(A.isInvertible) and A.isInvertible = false then
    A.determinant := AMatOps.ZeroNC(A.char);
    return A.determinant;
  fi;
  if IsBound(A.rank) and A.rank < A.dimensions[1] then
    A.determinant := AMatOps.ZeroNC(A.char);
    return A.determinant;
  fi;

  if A.type = "perm" then 

    d := SignPerm(A.element) * AMatOps.OneNC(A.char);

  elif A.type = "mon" then

    d := DeterminantMon(A.element);

  elif A.type = "mat" then

    d := DeterminantMat(A.element);
  
  elif A.type = "scalarMultiple" then

     d := (A.scalar ^ A.dimensions[1]) * DeterminantAMat(A.element);
    
  elif A.type = "product" then

    if ForAll(A.factors, AMatOps.IsSquare) then
      d := Product(List(A.factors, DeterminantAMat));
    else
      d := DeterminantMat(MatAMat(A));
    fi;

  elif A.type = "power" then

     d := DeterminantAMat(A.element) ^ A.exponent;

  elif A.type = "conjugate" then

     d := DeterminantAMat(A.element);

  elif A.type = "directSum" then

    if ForAny(A.summands, s -> not AMatOps.IsSquare(s)) then
      d := AMatOps.ZeroNC(A.char);
    else
      d := Product(List(A.summands, DeterminantAMat));
    fi;

  elif A.type = "tensorProduct" then

    d := AMatOps.OneNC(A.char);
    for i in [1..Length(A.factors)] do
      d := d * 
        DeterminantAMat(A.factors[i]) ^ 
        Product(
          List(
            Sublist(
              A.factors, 
              Concatenation([1..i-1], [i+1..Length(A.factors)])
            ),
            f -> f.dimensions[1]
          )
        );
    od;    

  elif A.type = "galoisConjugate" then
    
    d := 
      AMatOps.GaloisConjugation(
        DeterminantAMat(A.element), 
        A.galoisAut
      ); 

  else

    Error("unknown type of amat <A>");

  fi;

  # store and return result
  A.determinant := d;
  return d;
end;

#F TraceAMat( <amat> )
#F   calculates the trace of <amat> if it is square
#F   and stores the result in <amat>.trace
#F

TraceAMat := function ( A )
  local t, i;

  if not IsAMat(A) then
    Error("usage: TraceAMat( <amat> )");
  fi;
  if not AMatOps.IsSquare(A) then
    Error("<A> is not square");
  fi;

  # check if determinant is easy to calculate
  if IsBound(A.trace) then
    return A.trace;
  fi;  

  # dispatch on the type
  if A.type = "perm" then 

    t := AMatOps.ZeroNC(A.char);
    for i in [1..A.dimensions[1]] do
      if i^A.element = i then
        t := t + AMatOps.OneNC(A.char);
      fi;
    od;

  elif A.type = "mon" then

    t := TraceMon(A.element);

  elif A.type = "mat" then

    t := TraceMat(A.element);
  
  elif A.type = "scalarMultiple" then

    t := A.scalar * TraceAMat(A.element);
    
  elif A.type = "product" then

    t := TraceMat(MatAMat(A));

  elif A.type = "power" then

    t := TraceMat(MatAMat(A));

  elif A.type = "conjugate" then

    t := TraceAMat(A.element);

  elif A.type = "directSum" then

    if ForAll(A.summands, AMatOps.IsSquare) then
      t := Sum(List(A.summands, TraceAMat));
    else
      t := TraceMat(MatAMat(A));
    fi;

  elif A.type = "tensorProduct" then

    if ForAll(A.factors, AMatOps.IsSquare) then
      t := Product(List(A.factors, TraceAMat));
    else
      t := TraceMat(MatAMat(A));
    fi;

  elif A.type = "galoisConjugate" then
    
    t := 
      AMatOps.GaloisConjugation(
        TraceAMat(A.element), 
        A.galoisAut
      ); 

  else

    Error("unknown type of amat <A>");

  fi;

  # store and return result
  A.trace := t;
  return t;
end;


#F RankAMat( <amat> )
#F   calculates the rank of <amat>. The result is stored
#F   in A.rank.
#F

RankAMat := function ( A )
  local r;

  if not IsAMat(A) then
    Error("usage: TraceAMat( <amat> )");
  fi;

  # check if rank is easy to calculate
  if IsBound(A.rank) then
    return A.rank;
  fi;
  if 
    IsBound(A.determinant) and 
    A.determinant <> 0*A.determinant 
  then
    A.rank := A.dimensions[1];
    return A.rank;
  fi;  
  if IsBound(A.isInvertible) and A.isInvertible = true then
    A.rank := A.dimensions[1];
    return A.rank;
  fi;

  # dispatch on the types
  if A.type in ["perm", "mon"] then 

    r := A.dimensions[1];

  elif A.type = "mat" then

    r := RankMat(A.element);
  
  elif A.type = "scalarMultiple" then

    if A.scalar = AMatOps.ZeroNC(A.char) then
      r := 0;
    else
      r := RankAMat(A.element);
    fi;
    
  elif A.type in ["product", "power"] then

    r := RankMat(MatAMat(A));

  elif A.type = "conjugate" then

    r := RankAMat(A.element);

  elif A.type = "directSum" then

    r := Sum(List(A.summands, RankAMat));

  elif A.type = "tensorProduct" then

    r := RankMat(MatAMat(A));

  elif A.type = "galoisConjugate" then
    
    r := RankAMat(A.element);

  else

    Error("unknown type of amat <A>");

  fi;

  # store and return result
  A.rank := r;
  return r;
end;


# CharacteristicPolynomialAMat( <amat> )
#   ... noch nicht
#


#F Flattening out AMats
#F --------------------
#F
#F IsPermMat( <amat> )
#F IsMonMat( <amat> )
#F   test if <amat> can be flattened into a "perm" or "mon" amat.
#F   The result is memorized in <amat>.isPermMat/.isMonMat. 
#F   A "mon" amat is always invertible!
#F   Note that the names of the operations are not IsPermAMat, 
#F   IsMonAMat since <amat> can be of any type but represents
#F   a permutation/monomial matrix in a mathematical sense.
#F

if not IsBound(IsMonMat) then
  IsMonMat := "defined below";
  PermAMat := "defined below";
 MonAMat  := "defined below";
fi;

IsPermMat := function ( A )
  if not IsAMat(A) then 
    Error("usage: IsPermMat( <amat> )");
  fi;
  if IsBound(A.isPermMat) then
    return A.isPermMat;
  fi;
  A.isPermMat :=
    IsMonMat(A) and
    PermMon(MonAMat(A)) <> false;
  return A.isPermMat;
end;

IsMonMat := function ( A )
  if not IsAMat(A) then 
    Error("usage: IsMonMat( <amat> )");
  fi;
  if IsBound(A.isMonMat) then
    return A.isMonMat;
  fi;
  if IsBound(A.isPermMat) and A.isPermMat then
    A.isMonMat := true;
    return A.isMonMat;
  fi;
  if not AMatOps.IsSquare(A) then
    A.isMonMat := false;
    return A.isMonMat;
  fi;

  # dispatch on the types
  if A.type = "perm" then 

    A.isMonMat := true;

  elif A.type = "mon" then

    A.isMonMat := true;

  elif A.type = "mat" then

    A.mon      := MonMat(A.element);
    A.isMonMat := A.mon <> false;  
  
  elif A.type = "scalarMultiple" then

    A.isMonMat := 
      A.scalar <> 0*A.scalar and
      IsMonMat(A.element);
    
  elif A.type = "product" then

    if ForAll(A.factors, IsMonMat) then
      A.isMonMat := true;
    else
      A.isMonMat := MonAMat(A) <> false;
    fi;

  elif A.type = "power" then

    if IsMonMat(A.element) then
      A.isMonMat := true;
    else
      A.isMonMat := MonAMat(A) <> false;
    fi;

  elif A.type = "conjugate" then

    if IsMonMat(A.element) and IsMonMat(A.conjugation) then
      A.isMonMat := true;
    else
      A.isMonMat := MonAMat(A) <> false;
    fi;

  elif A.type = "directSum" then

    A.isMonMat := ForAll(A.summands, IsMonMat);

  elif A.type = "tensorProduct" then

    A.isMonMat := ForAll(A.factors, IsMonMat);

  elif A.type = "galoisConjugate" then
    
    A.isMonMat := IsMonMat(A.element);

  else
    Error("unknown type of amat <A>");
  fi;
  return A.isMonMat;
end;

#F PermAMat( <amat> )
#F MonAMat( <amat> )
#F MatAMat( <amat> )
#F   converts <amat> into a Perm/Mon/Mat-object and memorizes
#F   the result in <amat>.perm/mon/mat. If <amat> cannot be 
#F   converted then false is returned.
#F

PermAMat := function ( A )
  if A = false then
    return false;
  fi;
  if not IsAMat(A) then
    Error("usage: PermAMat( <amat> )");
  fi;
  
  if IsBound(A.isPermMat) and not A.isPermMat then
    return false;
  fi;
  if IsBound(A.isMonMat) and not A.isMonMat then
    return false;
  fi;
  if IsBound(A.perm) then
    return A.perm;
  fi;
  if not AMatOps.IsSquare(A) then
    A.isPermMat := false;
    return false;
  fi;

  A.perm      := PermMon(MonAMat(A));
  A.isPermMat := (A.perm <> false);
  return A.perm;
end;


MonAMat := function ( A )
  if A = false then
    return false;
  fi;
  if not IsAMat(A) then
    Error("usage: MonAMat( <amat> )");
  fi;

  if IsBound(A.isMonMat) and A.isMonMat = false then
    return false;
  fi;
  if IsBound(A.mon) then
    return A.mon;
  fi;
  if IsBound(A.perm) and A.perm <> false then
    return Mon(A.perm, A.dimensions[1], A.char);
  fi;
  if not AMatOps.IsSquare(A) then
    A.isPermMat := false;
    return false;
  fi;

  # dispatch on the types calculating the mon-object P
  # represented by A if possible, if not P = false
  if A.type = "perm" then

    A.mon := Mon(A.element, A.dimensions[1], A.char);

  elif A.type = "mon" then

    A.mon := A.element;

  elif A.type = "mat" then

    if 
      A.dimensions[1] > 1 and 
      IsBound(A.isDFT) and 
      A.isDFT = true 
    then
      A.mon := false;
    elif
      A.dimensions[1] > 1 and 
      IsBound(A.isSOR) and 
      A.isSOR = true 
    then
      A.mon := false;
    else
      A.mon := MonMat(A.element);
    fi;

  elif A.type = "scalarMultiple" then

    if A.scalar = 0*A.scalar then
      A.mon := false;
    else
      if IsMonMat(A.element) then
        A.mon := A.scalar * MonAMat(A.element);
      else
        A.mon := false;
      fi;
    fi;

  elif A.type = "product" then

    if ForAll(A.factors, IsMonMat) then
      A.mon := Product(List(A.factors, MonAMat));
    else
      A.mon := MonMat(MatAMat(A));
    fi;

  elif A.type = "power" then

    if IsMonMat(A.element) then
      A.mon := MonAMat(A.element) ^ A.exponent;
    else
      A.mon := MonMat(MatAMat(A));
    fi;

  elif A.type = "conjugate" then

    if IsMonMat(A.element) and IsMonMat(A.conjugation) then
      A.mon := MonAMat(A.element) ^ MonAMat(A.conjugation);
    else
      A.mon := 
        MonMat(
          MatAMat(InverseAMat(A.conjugation)) *
          MatAMat(A.element) *
          MatAMat(A.conjugation)
        );
    fi;

  elif A.type = "directSum" then

    if ForAll(A.summands, AMatOps.IsSquare) then
      if ForAll(A.summands, IsMonMat) then
        A.mon := DirectSumMon(List(A.summands, MonAMat));
      else
        A.mon := false;
      fi;
    else
      A.mon := false;
    fi;

  elif A.type = "tensorProduct" then

    if ForAll(A.factors, AMatOps.IsSquare) then
      if ForAll(A.factors, IsMonMat) then
        A.mon := TensorProductMon(List(A.factors, MonAMat));
      else
        A.mon := false;
      fi;
    else
      A.mon := false;
    fi;

  elif A.type = "galoisConjugate" then

    if IsMonMat(A.element) then
      A.mon := GaloisMon(MonAMat(A.element), A.galoisAut);
    else
      A.mon := false;
    fi;

  else
    Error("unknown type of amat <A>");
  fi;
  A.isMonMat := (A.mon <> false);
  return A.mon;
end;


MatAMat := function ( A )
  if A = false then
    return false;
  fi;
  if not IsAMat(A) then
    Error("usage: MatAMat( <amat> )");
  fi;
  
  if IsBound(A.mat) then 
    return A.mat;
  fi;
  if IsBound(A.mon) and A.mon <> false then
    A.mat := MatMon(A.mon);
    return A.mat;
  fi;
  if IsBound(A.perm) and A.perm <> false then
    A.mat := MatPerm(A.perm, A.dimensions[1], A.char);
    return A.mat;
  fi;

  # dispatch on the types
  if A.type = "perm" then

    A.mat := MatPerm(A.element, A.dimensions[1], A.char);

  elif A.type = "mon" then

    A.mat := MatMon(A.element);

  elif A.type = "mat" then

    A.mat := A.element;

  elif A.type = "scalarMultiple" then

    A.mat := A.scalar * MatAMat(A.element);

  elif A.type = "product" then

    A.mat := Product(List(A.factors, MatAMat));

  elif A.type = "power" then

# ...schlauer machen: ProductNC, InverseNC
    if A.exponent < 0 then
      A.mat :=
        MatAMat(InverseAMat(A.element)) ^ (-A.exponent);
    else
      A.mat := 
        MatAMat(A.element) ^ A.exponent;
    fi;

  elif A.type = "conjugate" then

    # use InverseAMat to use the structure of A.conjugation
    A.mat := 
      MatAMat(InverseAMat(A.conjugation)) *
      MatAMat(A.element) *
      MatAMat(A.conjugation);

  elif A.type = "directSum" then

    A.mat := ApplyFunc(DiagonalMat,List(A.summands, MatAMat));

  elif A.type = "tensorProduct" then

    A.mat := TensorProductMat(List(A.factors, MatAMat));

  elif A.type = "galoisConjugate" then

    A.mat :=  
      List(
        MatAMat(A.element), 
        Ai -> 
          List(
            Ai, 
            Aij -> AMatOps.GaloisConjugation(Aij, A.galoisAut)
          )
      );

  else
    Error("unknown type of amat <A>");
  fi;
  return A.mat;
end;

#F PermAMatAMat( <amat> )
#F MonAMatAMat( <amat> )
#F MatAMatAMat( <amat> )
#F   construct a "perm"/"mon"/"mat"-amat equal to <amat> if 
#F   possible or returns false otherwise.
#F

PermAMatAMat := function ( A )
  if not IsAMat(A) then
    Error("<A> must be an AMat");
  fi;
  return AMatPerm(PermAMat(A), A.dimensions[1], A.char);
end;

MonAMatAMat := function ( A )
  if not IsAMat(A) then
    Error("<A> must be an AMat");
  fi;
  return AMatMon(MonAMat(A));
end;

MatAMatAMat := function ( A )
  if not IsAMat(A) then
    Error("<A> must be an AMat");
  fi;
  return AMatMat(MatAMat(A));
end;


#F Simplifying AMats
#F -----------------
#F

#F SimplifyAMat( <amat> )
#F   simplifies <amat>, e.g. by removing identity matrices
#F   in products.
#F

SimplifyAMat := function ( A )
  local SA, SA1, L, L1, S,  s, A1, l, i, exp;

  if not IsAMat(A) then
    Error("usage: SimplifyAMat( <amat> )");
  fi;

  # look if A is already simplified
  if IsBound(A.isSimplified) and A.isSimplified = true then
    return A;
  fi;
  
  # easy cases
  if IsPermMat(A) then
    return PermAMatAMat(A);
  elif IsMonMat(A) then
    return MonAMatAMat(A);
  fi;

  # dispatch
  if A.type in ["perm", "mon", "mat"] then

    # do nothing
    return A;
  
  elif A.type = "scalarMultiple" then
  
    # A.scalar = 0 -> NullAMat, 
    if A.scalar = 0 * A.scalar then
      return NullAMat(A.dimensions, A.char);

    # A.scalar = 1 -> omit 
    elif A.scalar = A.scalar ^ 0 then
      SA := SimplifyAMat(A.element);

    else

      # recurse
      SA1 := SimplifyAMat(A.element);

      # if SA1 is of type "scalarMultiple", too,
      # then multiply
      if SA1.type = "scalarMultiple" then
	SA := 
	  ScalarMultipleAMat(A.scalar*SA1.scalar, SA1.element);
      else
	SA := ScalarMultipleAMat(A.scalar, SA1);
      fi;
    fi;

  elif A.type = "product" then

    # omit identity mats
    L := Filtered(A.factors, f -> not IsIdentityMat(f));
    if L = [ ] then
      return IdentityPermAMat(A.dimensions, A.char);
    fi;
    L := List(L, SimplifyAMat);

    # collect scalars
    L1 := [ ];
    s  := AMatOps.OneNC(A.char);
    for A1 in L do
      if A1.type = "scalarMultiple" then
        Add(L1, A1.element);
        s := s*A1.scalar;
      else
        Add(L1, A1);
      fi;
    od;
    L := L1;    

    if ForAll(L, AMatOps.IsSquare) then

      # multiply adjacent monomial matrices
      L1 := [ ];
      l  := IdentityPermAMat(A.dimensions[1], A.char);
      for i in [1..Length(L)] do
        if IsMonMat(L[i]) then
          l := MonAMatAMat(l * L[i]);
          if i = Length(L) then
	    if IsPermMat(l) then
	      Add(L1, PermAMatAMat(l));
	    else
	      Add(L1, l);
	    fi;
          fi;
        else
          if IsPermMat(l) then
            Add(L1, PermAMatAMat(l));
          else
            Add(L1, l);
          fi;
          Add(L1, L[i]);
          l := IdentityPermAMat(A.dimensions[1], A.char);
        fi;
      od;
      L := Filtered(L1, f -> not IsIdentityMat(f));
      if L = [ ] then
        return IdentityPermAMat(A.dimensions, A.char);
      fi;
    fi;

    SA := Product(L);

    # catch scalar computed above
    if s <> s^0 then
      SA := ScalarMultipleAMat(s, SA);
    fi; 

  elif A.type = "power" then

    # collect succeeding exponents
    exp := A.exponent;
    while A.element.type = "power" do
      A   := A.element;
      exp := exp * A.exponent;
    od;
    A := PowerAMat(A.element, exp);    

    # A.exponent = 0
    if A.exponent = 0 then
      SA := IdentityPermAMat(A.dimensions[1], A.char);
    fi;

    # A.exponent = 1 -> omit
    if A.exponent = 1 then
      SA := SimplifyAMat(A.element);
   
    # A.exponent = -1 -> InverseAMat
    elif A.exponent = -1 then
      SA := SimplifyAMat(InverseAMat(A.element));

    # A.exponent < 0 -> InverseAMat on A.element
    elif A.exponent < 0 then
      SA :=  
        PowerAMat(
          SimplifyAMat(InverseAMat(A.element)), 
          -A.exponent,
          "invertible"
        );

      # if SA.element is of type "scalarMultiple"
      # extract the scalar
      if SA.element.type = "scalarMultiple" then
        SA := 
          SimplifyAMat(
	    ScalarMultipleAMat(
	      SA.element.scalar^-A.exponent,
	      PowerAMat(
		SA.element.element,
		SA.exponent,
		"invertible"
	      )
	    )
          );

      # if SA.element is of type "power", multiply exponents
      elif SA.element.type = "power" then
        SA := 
          SimplifyAMat(
	    PowerAMat(
	      SA.element.element,
	      SA.element.exponent*SA.exponent,
	      "invertible"
	    )
          );
      fi;

    else

      SA :=  
        PowerAMat(
          SimplifyAMat(A.element), 
          A.exponent
        );  
  
      # if SA.element is of type "scalarMultiple"
      # extract the scalar
      if SA.element.type = "scalarMultiple" then
        SA := 
          SimplifyAMat(
	    ScalarMultipleAMat(
	      SA.element.scalar^A.exponent,
	      PowerAMat(
		SA.element.element,
		SA.exponent
	      )
	    )
          );

      # if SA.element is of type "power", multiply exponents
      elif SA.element.type = "power" then
        SA := 
          SimplifyAMat(
	    PowerAMat(
	      SA.element.element,
	      SA.element.exponent*SA.exponent
	    )
          );
      fi;

    fi;

  elif A.type = "conjugate" then

    # A.conjugation is identity mat -> omit it
    if IsIdentityMat(A.conjugation) then
      SA := SimplifyAMat(A.element);
    
    # A.conjugation is of type "scalarMultple"
    # -> omit the scalar
    elif A.conjugation.type = "scalarMultple" then
      SA :=  
        ConjugateAMat(
          SimplifyAMat(A.element), 
          SimplifyAMat(A.conjugation.element),
          "invertible"
        );

      # extract scalar if SA.element is of
      # type "scalarMultiple"
      if SA.element.type = "scalarMultiple" then
        SA :=
          SimplifyAMat(
            ScalarMultipleAMat(
              SA.element.scalar,
              ConjugateAMat(
                SA.element.element,
                SA.conjugation,
                "invertible"
              )
            )
          );
      fi;

    else
      SA :=  
        ConjugateAMat(
          SimplifyAMat(A.element), 
          SimplifyAMat(A.conjugation),
          "invertible"
        );

      # extract scalar if SA.element is of
      # type "scalarMultiple"
      if SA.element.type = "scalarMultiple" then
        SA :=
          SimplifyAMat(
            ScalarMultipleAMat(
              SA.element.scalar,
              ConjugateAMat(
                SA.element.element,
                SA.conjugation,
                "invertible"
              )
            )
          );
      fi;

    fi;

  elif A.type = "directSum" then
    
    if Length(A.summands) = 1 then

      SA := SimplifyAMat(A.summands[1]);

    else

      L := List(A.summands, SimplifyAMat);

      # build the direct sum from adjacent monomial matrices
      L1 := [ ];
      S  := [ ];
      for i in [1..Length(L)] do
        if IsMonMat(L[i]) then
          Add(S, L[i]);
        else
          if Length(S) >= 1 then
            Add(L1, SimplifyAMat(DirectSumAMat(S)));
            S := [ ];
          fi;
          Add(L1, L[i]);
        fi;
      od;
      if Length(S) >= 1 then
        Add(L1, SimplifyAMat(DirectSumAMat(S)));
      fi;
      
      SA := DirectSumAMat(L1);
    fi;

  elif A.type = "tensorProduct" then
  
    # omit identity mats of degree 1
    L := 
      Filtered(
        A.factors, 
        f -> not( IsIdentityMat(f) and f.dimensions[1] = 1 )
      );
    if L = [ ] then
      return IdentityPermAMat(A.dimensions, A.char);
    fi;

    if Length(L) = 1 then
      SA := SimplifyAMat(L[1]);
    else
      SA := TensorProductAMat(List(L, SimplifyAMat));
    fi;

  elif A.type = "galoisConjugate" then

    SA :=  
      GaloisConjugateAMat(
        SimplifyAMat(A.element), 
        A.galoisAut
      );

  else
    Error("unrecognized type of amat <A>");
  fi;

  # copy valuable information
  if IsBound(A.isPermMat) then
    SA.isPermMat := A.isPermMat;
  fi;
  if IsBound(A.isMonMat) then
    SA.isMonMat := A.isMonMat;
  fi;
  if IsBound(A.mat) then
    SA.mat := A.mat;
  fi;
  if IsBound(A.isInvertible) then
    SA.isInvertible := A.isInvertible;
  fi;
  if IsBound(A.inverse) then
    SA.inverse := A.inverse;
  fi;

  # set flag
  SA.isSimplified := true;

  return SA;
end;


#F Functions for AMats
#F -------------------
#F

#F kbsAMat( <amat1>, .., <amatN> ) ; N >= 1
#F kbsAMat( <list-of-amat> )
#F   calculates the joined kbs of the square amats 
#F   <amat1>, .., <amatN>, which must have common size.
#F   (cf. kbs() in permblk.g).
#F

kbsAMat := function ( arg )
  if IsList(arg[1]) then
    arg := arg[1];
  fi;
  if not ForAll(arg, IsAMat) then
    Error(
      "usage: \n",
      "  kbsAMat( <amat1>, .., <amatN> )\n",
      "  kbsAMat( <list-of-amat> )"
    );
  fi;

  return kbs(List(arg, MatAMat));
end;


#F SubmatrixAMat( <amat>, <list-of-indices> )
#F   calculates the submatrix of amat obtained by
#F   extracting all components with indices in <inds>.
#F   Note that the function first converts <amat>
#F   to a matrix.
#F   A "mat"-amat is returned.
#F

SubmatrixAMat := function ( A, inds )
  local M;

  if not( 
    IsAMat(A) and 
    IsList(inds) and 
    ForAll(inds, i -> IsInt(i) and i >= 1)
  ) then
    Error("usage: SubmatrixAMat( <amat>, <list-of-indices> )");
  fi;
  if not ForAll(inds, i -> i >= 0 and i <= Maximum(A.dimensions)) then
    Error("<inds> are to high");
  fi;

  # extract submatrix
  M := MatAMat(A);
  return
    AMatMat(
      Sublist(List(M, r -> Sublist(r, inds)), inds)
    );
end;

#F kbsDecompositionAMat( <amat> )
#F   decomposes the square <amat> into a conjugated 
#F   (by a "perm"-amat) direct sum of "mat"-amats as 
#F   far as possible.
#F

kbsDecompositionAMat := function ( A )
  local blocks, sortperm;

  if not IsAMat(A) then
    Error("usage: kbsDecompositionAMat( <amat> )");
  fi;
  if not AMatOps.IsSquare(A) then
    Error("<A> must be square");
  fi;

  blocks   := kbsAMat(A);
  sortperm := PermList(Concatenation(blocks));
  return
    ConjugateAMat(
      DirectSumAMat(
        List(blocks, b -> SubmatrixAMat(A, b))
      ),
      AMatPerm(sortperm, A.dimensions[1], A.char),
      "invertible"
    );
end;


#F UpperBoundLinearComplexityAMat( <amat> )
#F   computes an upper bound for the linear complexitiy
#F   according to the model of Clausen.
#F

UpperBoundLinearComplexityAMat := function ( A )
  local A1;

  if not IsAMat(A) then
    Error("usage: UpperBoundLinearComplexityAMat( <amat> )");
  fi;

  # dispatch on the types
  if A.type = "perm" then

    return 0;

  elif A.type = "mon" then

    return Length(Filtered(A.element.diag, x -> not x in [-1, 1]));

  elif A.type = "mat" then

    # check if DFT
    if IsBound(A.isDFT) then
      if A.dimensions[1] = 2 then
        return 2;
      elif A.dimensions[1] = 3 then
        return 8;
      elif A.dimensions[1] = 5 then
        return 22;
      fi;
    fi;

    return 2 * Product(A.dimensions);

  elif A.type = "scalarMultiple" then

    if A.scalar = 0 then
      return 0;
    elif A.scalar in [-1, 1] then 
      return UpperBoundLinearComplexityAMat(A.element);
    else
      return 
        Minimum(A.dimensions) + 
        UpperBoundLinearComplexityAMat(A.element);
    fi;

  elif A.type = "product" then
    
    return Sum(List(A.factors, UpperBoundLinearComplexityAMat));

  elif A.type = "power" then
 
    return A.exponent * UpperBoundLinearComplexityAMat(A.element);

  elif A.type = "conjugate" then
 
    return
      UpperBoundLinearComplexityAMat(A.conjugation) + 
      UpperBoundLinearComplexityAMat(A.element) + 
      UpperBoundLinearComplexityAMat(InverseAMat(A.conjugation));

  elif A.type = "directSum" then

    return Sum(List(A.summands, UpperBoundLinearComplexityAMat));

  elif A.type = "tensorProduct" then

    if Length(A.factors) = 1 then
      return UpperBoundLinearComplexityAMat(A.factors[1]);
    else
      A1 := 
        TensorProductAMat(
          Sublist(A.factors, [2..Length(A.factors)])
        );
      return
        A.factors[1].dimensions[1] * 
        UpperBoundLinearComplexityAMat(A1) +
        A1.dimensions[2] *
        UpperBoundLinearComplexityAMat(A.factors[1]);
    fi;

  elif A.type = "galoisConjugate" then

    return UpperBoundLinearComplexityAMat(A.element);

  else
  
    Error("unrecognized type of amat <A>");

  fi;
end;


#F AMatSparseMat( <mat> [, <match_blocks> ] )
#F   writes the given matrix <mat> as a permuted direct sum
#F   of copies of minimal direct summands. The output is
#F   always of the form
#F
#F     AMatPerm(permL) *
#F     IdentityAMat([nrL, ncL]) *
#F     DirectSumAMat(
#F       TensorProductAMat(IdentityAMat(n[1]), A[1]),
#F       ..
#F       TensorProductAMat(IdentityAMat(n[r]), A[r])
#F     ) *
#F     IdentityAMat([nrR, ncR]) *
#F     AMatPerm(permR)
#F
#F   The optional flag match_blocks (default is true)
#F   allows to suppress the possibly very slow matching
#F   of seperated blocks. In this case, the A[i] might be
#F   permuted copies of each other.
#F
#F   The function is very useful to write sparse matrices
#F   in a structured form. For a detailed description of
#F   the algorithm refer to Chapter 2 of the PhD thesis
#F   of S. Egner and to the file 'summands.g' in the AREP
#F   package. In particular, DirectSummandsPermutedMat does
#F   the major part of the work.
#F

AMatSparseMat := function ( arg )
  local
    mat, match_blocks,   # arguments
    nr, nc, char,        # nr. rows, columns, char of mat
    mat1,                # reduced matrix
    nr1, nc1,            # size of mat1
    R0, R1,              # all-zero row indices, complement
    C0, C1,              # all-zero column indices, complement
    permL1,              # left permutation of mat1
    DS1,                 # list of direct summands of mat1
    permR1,              # right permutation of mat1
    DS,                  # list of collected direct summands for mat1
    i, ni, Ai;           # current index, multiplicity, tensor factor

  # decode arg
  if Length(arg) = 1 then
    mat          := arg[1];
    match_blocks := true;
  elif Length(arg) = 2 then
    mat          := arg[1];
    match_blocks := arg[2];
  else
    Error("usage: AMatSparseMat( <mat> [, <match_blocks> ] )");
  fi;
  if not IsMat(mat) then
    Error("<mat> must be matrix (list of list of field elements)");
  fi;
  if not match_blocks in [false, true] then
    Error("<match_block> must be boolean");
  fi;

  # fetch size and char of mat
  nr   := Length(mat);
  nc   := Length(mat[1]);
  char := DefaultField(Concatenation(mat)).char;
 
  # split off all-zero rows and columns
  R0   := 
    Filtered(
      [1..nr], 
      r -> ForAll([1..nc], c -> mat[r][c] = 0*mat[r][c])
    );
  C0   := 
    Filtered(
      [1..nc], 
      c -> ForAll([1..nr], r -> mat[r][c] = 0*mat[r][c])
    );
  R1   := Difference([1..nr], R0);
  C1   := Difference([1..nc], C0);
  nr1  := Length(R1);
  nc1  := Length(C1);
  mat1 := List(R1, r -> Sublist(mat[r], C1));

  # handle the annoying all-zero case
  if nr1 = 0 or nc1 = 0 then
    return
      AMatPerm((), nr, char) *
      IdentityAMat([nr, Minimum(nr, nc)], char) *
      DirectSumAMat(
        TensorProductAMat(
          IdentityPermAMat(Minimum(nr, nc), char),
          AMatMat([[0]])
        )
      ) *
      IdentityAMat([Minimum(nr, nc), nc], char) *
      AMatPerm((), nc, char);
  fi;

  # extract non-zero part of mat and split it
  DS1 := 
    DirectSummandsPermutedMat(
      mat1,
      match_blocks
    );
  permL1 := DS1[1];
  permR1 := DS1[3];
  DS1    := DS1[2];

  # construct direct sum of tensor products
  DS := [ ];
  ni := 1;
  Ai := DS1[1];
  for i in [2..Length(DS1)] do
    if DS1[i] = Ai then
      ni := ni + 1;
    else
      Add( DS, 
       TensorProductAMat(
         IdentityPermAMat(ni, char),
         AMatMat(Ai)
       )
      );
      ni := 1;
      Ai := DS1[i];
    fi;
  od;
  Add( DS, 
    TensorProductAMat(
      IdentityPermAMat(ni, char),
      AMatMat(Ai)
    )
  );

  # construct result
  return
    AMatPerm( PermList(Concatenation(R1, R0))^-1 * permL1, nr, char) *
    IdentityAMat([nr, nr1], char) *
    DirectSumAMat(DS) *
    IdentityAMat([nc1, nc], char) *
    AMatPerm( permR1 * PermList(Concatenation(C1, C0)), nc, char);
end;

# Bsp:
# A1 := 
#   DirectSumAMat(
#     AMatPerm((1,2,3), 4)^IdentityMatAMat(4, 0),
#     3 * AMatMon(Mon((3,2), 3)),
#     GaloisConjugateAMat(
#       AMatMat([[1,2],[3,E(3)]]),
#       2
#     ) * 
#     AMatPerm((1,2), 2)^3
#   );
# A1 := TensorProductAMat(A1, A1);
