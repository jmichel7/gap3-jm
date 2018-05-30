
def f297118(n):
  a=IdentityMat(n)
  for i in range(1,n+1):
    a[i][i-1]=2
    if i<n :
      a[i][i+1-1]=-1
    if i>1 :
      a[i][i-1-1]=-1
  return a

ChevieData["A"]["CartanMat"]=f297118

ChevieData["A"]["ReflectionDegrees"]=lambda n: range(2,n+1+1)

def f229529(r,indices,title):
  print title," ",Join(indices," - "),"\n",

ChevieData["A"]["PrintDiagram"]=f229529

def f491287(l):
  r=map(lambda i: GAPMul(0,range(1,l+1+1)),range(1,l+1))
  for i in range(1,l+1):
    [r[i][k-1] for k in [i,i+1]]=[1,-1]
  return r

ChevieData["A"]["GeneratingRoots"]=f491287

def f562376(l,s):
  return ChevieData["imp"]["ParabolicRepresentatives"](1,1,l,s)

ChevieData["A"]["ParabolicRepresentatives"]=f562376

def f829110(pi):
  w=[]
  i=0
  for l in pi:
    r=l%2
    Append(w,i+Concatenation([1,..(3,l-1-r)],[2,..(4,l+r-2)]))
    i=i+l
  return w

ChevieData["A"]["WordClass"]=f829110

def f18591(n):
  res={"classparams":Partitions(n+1)}
  res["classnames"]=map(IntListToString,res["classparams"])
  res["classtext"]=map(ChevieData["A"]["WordClass"],res["classparams"])
  res["classes"]=map(lambda pi: GAPDiv(Factorial(n+1),CharTableSymmetric.centralizers[1-1](n,pi)),res["classparams"])
  res["orders"]=map(Lcm,res["classparams"])
  return res

ChevieData["A"]["ClassInfo"]=f18591

ChevieData["A"]["NrConjugacyClasses"]=lambda n: NrPartitions(n+1)

ChevieData["A"]["WeightInfo"]=lambda n: {"minusculeWeights":range(1,n+1),
  "decompositions":map(lambda i: [i],range(1,n+1)),
  "moduli":[n+1]}

def f342502(n,w):
  x=Permutation("()")
  for i in w:
    x=GAPMul(x,Permutation("(%s,%s)"%(i,i+1)))
  res=[]
  mark=range(1,n+1+1)
  for i in range(1,n+1+1):
    if !=(mark[i-1],0) :
      cyc=CyclePermInt(x,i)
      Add(res,len(cyc))
      [mark[k-1] for k in cyc]=GAPMul(cyc,0)
  Sort(res)
  return Reversed(res)

ChevieData["A"]["ClassParameter"]=f342502

ChevieData["A"]["CharParams"]=lambda n: Partitions(n+1)

ChevieData["A"]["LowestPowerFakeDegree"]=lambda p: GAPMul(p,range(0,len(p)-1+1))

ChevieData["A"]["HighestPowerFakeDegree"]=lambda p: GAPDiv(GAPMul(Sum(p),Sum(p)-1),2)-ChevieData["A"]["LowestPowerFakeDegree"](AssociatedPartition(p))

def f416799(arg):
  return IntListToString(arg[2-1])

ChevieData["A"]["CharName"]=f416799

def f69333(n):
  res={"charparams":ChevieData["A"]["CharParams"](n)}
  res["extRefl"]=map(lambda i: Position(res["charparams"],Concatenation([n+1-i],GAPMul(range(1,i+1),0)+1)),range(0,n+1))
  res["b"]=map(ChevieData["A"]["LowestPowerFakeDegree"],res["charparams"])
  res["B"]=map(ChevieData["A"]["HighestPowerFakeDegree"],res["charparams"])
  res["a"]=res["b"]
  res["A"]=res["B"]
  return res

ChevieData["A"]["CharInfo"]=f69333

ChevieData["A"]["CharTable"]=ChevieData["compat"]["CharTableA"]

ChevieData["tmp"]=ShallowCopy(CharTableSymmetric)

ChevieData["tmp"]["identifier"]="HeckeA"

ChevieData["tmp"]["specializedname"]=lambda nq: SPrint("H(A",nq[1-1]-1,")")

ChevieData["tmp"]["order"]=lambda nq: Factorial(nq[1-1])

ChevieData["tmp"]["size"]=lambda nq: Factorial(nq[1-1])

ChevieData["tmp"]["domain"]=lambda nq: IsList(nq) and len(nq)==2 and IsInt(nq[1-1]) and nq[1-1]>0

ChevieData["tmp"]["text"]="generic character table of Hecke algebras of type_ A"

ChevieData["tmp"]["classparam"]=[lambda nq: Partitions(nq[1-1])]

ChevieData["tmp"]["charparam"]=[lambda nq: Partitions(nq[1-1])]

def f101237(nq,gamma,pi):
  n=nq[1-1]
  q=nq[2-1]
  if n==0 :
    return q**0
  k=pi[1-1]
  val=GAPMul(0,q)
  AHk=CHEVIE.R("Hk", "A").irreducibles[1][1-1]
  for alpha in Partitions(n-k):
    dif=DifferencePartitions(gamma,alpha)
    if !=(dif,false) :
      val=val+GAPMul(q-1**dif["cc"]-1,-1**dif["ll"])
  return val

ChevieData["tmp"]["irreducibles"]=[[f101237]]

def CHEVIE.tmp.matrix(nq):
  n=nq[1-1]
  q=nq[2-1]
  pm=[]
  scheme=[]
  def hooks(beta,m):
    prs=[]
    for i in beta:
      leg=0
      for j in [i-1,..(i-2,0)]:
        if in(j,beta) :
          leg=leg+1
        else:
          Add(prs,[{"from":i,
            "to":j,
            "leg":leg}])
    cbs=ShallowCopy(prs)
    hks=map(lambda x: [],range(1,m+1))
    for hk in cbs:
      for pr in prs:
        if pr[1-1]["to"]>hk[len(hk)-1]["from"] :
          Add(cbs,Concatenation(hk,pr))
      ll=Sum(hk,lambda x: x["from"]-x["to"])
      lg=Sum(hk,lambda x: x["leg"])
      lh=len(hk)
      new={"wgt":GAPMul(-1**lg,q**ll-lg-lh),
        "adr":1}
      if ll<m-1 :
        gamma=Difference(beta,map(lambda x: x["from"],hk))
        UniteSet(gamma,map(lambda x: x["to"],hk))
        if in(0,gamma) :
          j=0
          while gamma[j+1-1]==j:
            j=j+1
          gamma=[gamma[k-1] for k in range(j+1,len(gamma)+1)]-j
        new["adr"]=Position(pm[m-ll-1],gamma)
      Add(hks[ll-1],new)
    return hks
  
  
  for i in range(1,n+1):
    pm[i-1]=map(BetaSet,Partitions(i))
    scheme[i-1]=[]
    for beta in pm[i-1]:
      Add(scheme[i-1],hooks(beta,i))
  def charCol(n,t,k):
    col=[]
    for pi in scheme[n-1]:
      val=GAPMul(0,q)
      for hk in pi[k-1]:
        val=val+GAPMul(hk["wgt"],t[hk["adr"]-1])
      Add(col,val)
    return col
  
  
  pm=map(lambda x: [],range(1,n-1+1))
  for m in range(1,QuoInt(n,2)+1):
    Add(pm[m-1],charCol(m,[1],m))
    for k in range(m+1,n-m+1):
      for t in pm[k-m-1]:
        Add(pm[k-1],charCol(k,t,m))
  res=[]
  for k in range(1,n-1+1):
    for t in pm[n-k-1]:
      Add(res,charCol(n,t,k))
  Add(res,charCol(n,[1],n))
  return Matrix(res).transpose()



ChevieData["A"]["Hk"]=ShallowCopy(ChevieData["tmp"])

ChevieData["A"]["HeckeCharTable"]=ChevieData["compat"]["HeckeCharTableA"]

def f227941(n,param):
  return prod(range(1,n+1))

ChevieData["A"]["PoincarePolynomial"]=f227941

def f110529(n,alpha,param,sqrtparam):
  q=GAPDiv(-param[1][1-1],param[1][2-1])
  lambda=BetaSet(alpha)
  res=q**Binomial(len(lambda),3)
  for i in lambda:
    for j in range(0,i-1+1):
      if in(j,lambda) :
        res=GAPDiv(res,q**j)
      else:
        res=GAPMul(res,Sum(range(0,i-j-1+1),lambda e: q**e))
  return res

ChevieData["A"]["SchurElement"]=f110529

def f812103(arg):
  return ChevieData["imp"]["FactorizedSchurElement"](1,1,arg[1-1]+1,[arg[2-1]],arg[3-1],[])

ChevieData["A"]["FactorizedSchurElement"]=f812103

def f913600(n,param,sqrtparam,i):
  H=Hecke(CoxeterGroup("A",n),GAPDiv(-param[1][1-1],param[1][2-1]))
  return SpechtModel(H,Partitions(n + 1)[i-1])

ChevieData["A"]["HeckeRepresentation"]=f913600

def f74465(n,i):
  return [(CHEVIE.R("Representation", "imp"))(1, 1, n + 1, i)[k-1] for k in range(2,n+1+1)]

ChevieData["A"]["Representation"]=f74465

def f422080(n,p,q):
  return GAPDiv(ChevieData["A"]["PoincarePolynomial"](Sum(p)-1,[[q,-1]]),ChevieData["A"]["SchurElement"](Sum(p)-1,p,[[q,-1]],[]))

ChevieData["A"]["FakeDegree"]=f422080

def f596183(l,p):
  return [[range(1,NrPartitions(l+1)+1),MatrixDecompositionMatrix(DecompositionMatrix(Specht(p,p),l+1))]]

ChevieData["A"]["DecompositionMatrix"]=f596183

def f73313(l):
  ci=ChevieData["A"]["CharInfo"](l)
  return {"harishChandra":[{"levi":[],
    "relativeType":{"series":"A",
    "indices":range(1,l+1),
    "rank":l},
    "parameterExponents":GAPMul(0,range(1,l+1))+1,
    "cuspidalName":"",
    "eigenvalue":1,
    "charNumbers":range(1,len(ci["charparams"])+1)}],
    "families":map(lambda i: Family("C1",[i]),range(1,len(ci["charparams"])+1)),
    "charParams":ci["charparams"],
    "charSymbols":map(lambda x: [BetaSet(x)],ci["charparams"]),
    "a":ci["a"],
    "A":ci["A"]}

ChevieData["A"]["UnipotentCharacters"]=f73313

def f247986(n):
  m=ChevieData["A"]["GeneratingRoots"](n)
  Add(m,GAPMul(range(1,n+1+1),0)+1)
  return map(lambda i: f675256,range(2,n+1+1))

def f675256(arg):
  v=ShallowCopy(arg)
  Add(v,GAPMul(0,v[1-1]))
  v=GAPMul(v,m)
  return Sum(Arrangements(range(1,n+1+1),i),lambda a: prod([v[k-1] for k in a]))

ChevieData["A"]["Invariants"]=f247986

def f231630(n,p):
  uc={"classes":map(lambda p: {"parameter":p},Partitions(n+1)),
    "springerSeries":Concatenation(map(lambda d: map(lambda i: {"relgroup":CoxeterGroup("A",GAPDiv(n+1,d)-1),
    "Z":[ER(d)**i],
    "levi":Filtered(range(1,n+1+1),lambda i: !=(i%d,0)),
    "locsys":[]},PrimeResidues(d)),DivisorsInt(n+1)))}
  ss=lambda z: First(uc["springerSeries"],lambda x: x["Z"]==[z])
  for i in range(1,len(uc["classes"])+1):
    cl=uc.classes[i-1]
    p=cl["parameter"]
    d=Gcd(p)
    cl["name"]=IntListToString(p)
    cl["Au"]=ComplexReflectionGroup(d,1,1)
    cl["balacarter"]=Concatenation(map(lambda i: Sum([p[k-1] for k in range(1,i-1+1)])+range(1,p[i-1]-1+1),range(1,len(p)+1)))
    p=Concatenation(map(lambda x: [1-x,..(3-x,x-1)],p))
    Sort(p)
    cl["dynkin"]=map(lambda i: p[i+1-1]-p[i-1],range(1,len(p)-1+1))
    cl["red"]=[]
    p=1
    for j in Collected(cl["parameter"]):
      Append(cl["red"],range(p,p+j[2-1]-2+1))
      p=p+j[2-1]
    cl["red"]=ReflectionSubgroup(CoxeterGroup("A",p-2),cl["red"])
    cl["AuAction"]=ExtendedReflectionGroup(cl["red"],[IdentityMat(cl["red"]["rank"])])
    if d==2 :
      Add(ss(1)["locsys"],[i,2])
      Add(ss(-1)["locsys"],[i,1])
    else:
      for j in range(0,d-1+1):
        Add(ss(ER(d)**j)["locsys"],[i,j+1])
  uc["orderClasses"]=Hasse(Poset(map(lambda x: map(lambda y: Dominates(y["parameter"],x["parameter"]),uc["classes"]),uc["classes"])))
  return uc

ChevieData["A"]["UnipotentClasses"]=f231630

def f58629(n):
  def f(i):
    if !=(i,Permutation("()")) :
      i=prod(CoxeterWord(W,i))
    i=map(Length,RobinsonSchenstedCorrespondent(n+1,i)["P"])
    return Position(CharParams(W),[i])
  
  
  W=CoxeterGroup("A",n)
  l=Filtered(Elements(W),lambda x: x**2==Permutation("()"))
  return map(lambda x: {"duflo":OnTuples(range(1,n+1),x),
    "reps":[],
    "character":[f(x)]},l)

ChevieData["A"]["KLeftCellRepresentatives"]=f58629
