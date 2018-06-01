from collections import defaultdict 
from sage.rings.universal_cyclotomic_field import E

CharTableA = CharTableB = CharTableD = lambda *args: {}
CharTableSymmetric = CharTableWeylA = CharTableWeylB = CharTableWeylD = {}

ChevieData = defaultdict(dict)
ChevieData["compat"] = defaultdict(dict)
ChevieData["compat"]["CharTableB"] = CharTableB
ChevieData["compat"]["HeckeCharTableB"] = CharTableB

def ER(n):
  if n == 0:
      return 0
  elif n < 0:
      return E(4)*ER(-n)
  elif n%4==1:
      return sum( E(n,k**2) for k in range(1,n+1) )
  elif n%4==2:
      return (E(8)-E(8,3))*ER(n/2)
  elif n%4==3:
      return -E(4)*sum(E(n,k**2) for k in range(1,n+1))
  else:
      return 2*ER(n/4)

def just(st,n):
    if n >= 0:
        return st.rjust(n)
    else:
        return st.ljust(n)

def IntListToString(L):
    if any(a >= 10 for a in L):
        return ",".join(str(a) for a in L)
    else:
        return "".join(str(a) for a in L)

def GAPDiv(L,k):
    if isinstance(k,int):
        k = Integer(k)
    try:
        return L / k
    except TypeError:
        return [ GAPDiv(a,k) for a in L ]

def GAPMul(A,B):
    try:
        return A*B
    except TypeError:
        if isinstance(A,list):
            if isinstance(A[0],list):
                A = Matrix(A)
            else:
                A = vector(A)

        if isinstance(B,list):
            if isinstance(B[0],list):
                B = Matrix(B)
            else:
                B = vector(B)
        return A*B

#def GAPMin(A,B):
    #if isinstance(A,list) and isinstance(B,int):
        #return [a-B for a in A]
    #if isinstance(A,list) and isinstance(B,list):
        #return [a-b for a,b in zip(A,B)]
    #if isinstance(A,int) and isinstance(B,list):
        #return [A-b for b in B]
    #if isinstance(A,int) and isinstance(B,int):
        #return a-b

def Reflection(root, coroot):
    g = vector(root)
    r = vector(coroot)
    return identity_matrix(len(root)) - Matrix([ i*g for i in r])
    
def BuildPermutationRepresentation(roots, coroots):
    # only for complex non-real
    n = len(roots)
    roots = [ vector(root) for root in roots ]
    matgens = [ Reflection(*x) for x in zip(roots,coroots) ]

    images = []
    for i in range(n):
        images.append([])

    allroots = copy(roots)


    newroots = True
    while newroots:
        newroots = False
        for j in range(n):
            r = len(allroots)
            for y in [ root*matgens[j] for root in allroots[ len(images[j]):] ]:
                try:
                    p = allroots[:r].index(y)
                except ValueError:
                    p = len(allroots)
                    allroots.append(y)
                    newroots = True
                images[j].append(p+1)
    return [Permutation(pi) for pi in images]

def RootsCartanWeyl(C):
    # only Weyl
    Rall = []
    for L in DecomposeMatrix(C):
        CL = C.matrix_from_rows_and_columns(L,L)
        R = RootsIrredWeyl(CL)
        for r in R:
            v = [0]*C.nrows()
            for i,j in enumerate(L):
                v[j] = r[i]
            Rall.append(vector(v))
    return sorted(Rall,key=lambda x:(sum(x),-x))

def AllRootsCartan(C):
    Phiplus = RootsCartan(C)
    Phiminus = [ -beta for beta in Phiplus ]
    return Phiplus + Phiminus

def RootsIrredWeyl(C):
    n = C.nrows()
    I = identity_matrix(n)
    old = I.rows()

    R = []

    while old:
        R.append(old)
        new = []
        for r in old:
            for j in range(n):
                p = C.row(j)*r
                if p < 0 or ( len(R) - p - 1 > 0 and r[j] > p and r-(p+1)*I.row(j) in R[-p-2]):
                    if not r+I.row(j) in new:
                        new.append(r+I.row(j))
        old = new
    return sum(R,[])

def DecomposeMatrix(M):
    l = M.nrows()
    cc = range(l)
    nz = lambda x: bool(x)
    for i in range(l):
        for j in range(i,l):
            if (M[i,j] or M[j,i]) and cc[i] != cc[j]:
                cj = cc[j]
                for k in range(l):
                    if cc[k] == cj:
                        cc[k] = cc[i]
    return sorted(set(CollectBy(range(l),cc)))

def PermutationsSimpleReflections(C, Phi):
    # for Weyl
    N = len(Phi)
    n = C.nrows()
    gens = []
    for i in range(n):
        PhiNew = deepcopy(Phi)
        for j in range(N):
            PhiNew[j][i] -= Phi[j]*C.row(i)
        gens.append( Permutation([Phi.index(beta)+1 for beta in PhiNew]) )
    return gens

def RootsCartan(C):
    if isinstance(C,list):
        C = Matrix(C)
    # Cox not Weyl
    n = C.nrows()
    R = identity_matrix(n).rows()
    for a in R:
        for i in range(n):
            if a != R[i]:
                v = a - R[i]*(a*C[i])
                if not v in R:
                    R.append(v)
    return R

def CollectBy(A,B):
    ind = sorted(set(B))
    out = []
    for i in ind:
        out.append([])
        for j,x in enumerate(B):
            if x == i:
                out[-1].append(A[j])
    return [ tuple(x) for x in out ]

def Replace(*args):
    assert len(args) % 2 == 1
    args = list(args)
    X    = args[0]
    args = args[1:]
    if isinstance(X, str):
        X = list(X)
    for i,Y in enumerate(args):
        if isinstance(Y, str):
            args[i] = list(Y)

    assert isinstance(X, list)
    assert all( isinstance(Y,list) for Y in args)

    while args:
        Y1,Y2 = args[:2]
        args  = args[2:]
        X = replace_sublist(X, Y1, Y2)
    return X

def replace_sublist(seq, sublist, replacement):
    seq = copy(seq)
    length = len(replacement)
    index = 0
    for start, end in iter(lambda: find_first_sublist(seq, sublist, index), None):
        seq[start:end] = replacement
        index = start + length
    return seq

def find_first_sublist(seq, sublist, start=0):
    length = len(sublist)
    for index in range(start, len(seq)):
        if seq[index:index+length] == sublist:
            return index, index+length

def import_all():
    import os
    files = [f for f in os.listdir('.') if os.path.isfile(f)]
    errs = []
    for f in files:
        if f != "cmplxaux.py":
            #print "Importing %s"%f
            try:
                load(f)
            except StandardError as err:
                print "Importing %s failed!"%f
                print "Error message:"
                print err
                errs.append(f)
    return errs

def DiagonalMat(*args):
    M = args[0]
    M = Matrix(M)
    for y in args[1:]:
        y = Matrix(y)
        M = M.block_sum(y)
    return M

def ReadChv(*args):
    return

def ComplexConjugate(X):
    from sage.functions.other import conjugate
    if isinstance(X,list):
        return [ComplexConjugate(x) for x in X]
    else:
        return conjugate(X)

def Denominator(X):
    if hasattr(X,"denominator"):
        return X.denominator()
    else:
        return X

def Numerator(X):
    if hasattr(X,"numerator"):
        return X.numerator()
    else:
        return X

def eptio2(t):
    r=ChevieData[t]["GeneratingRoots"]
    rbar=ComplexConjugate(r)
    e=ChevieData[t]["EigenvaluesGeneratingReflections"]
    e=[1-x for x in map(lambda x: ER(Denominator(x))**Numerator(x),e)]
    e=map(lambda i: GAPDiv(GAPMul(e[i-1],rbar[i-1]),GAPMul(rbar[i-1],r[i-1])),range(1,len(e)+1))
    return map(lambda x: map(lambda y: GAPMul(x,y),r),e)

def SPrint(*args):
    return "".join(str(a) for a in args)
