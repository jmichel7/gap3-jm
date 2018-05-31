
def weyla1(n):
    a=IdentityMat(n)
    for i in range(1,n+1):
        a[i-1][i-1]=2
        if i<n :
            a[i-1][i+1-1]=-1
        if i>1 :
            a[i-1][i-1-1]=-1
    return a

ChevieData["A"]["CartanMat"]=weyla1

ChevieData["A"]["ReflectionDegrees"]=lambda n: range(2,n+1+1)

def weyla2(r,indices,title):
    print title," ",Join(indices," - "),"\n",

ChevieData["A"]["PrintDiagram"]=weyla2

def weyla3(l):
    r=map(lambda i: GAPMul(0,range(1,l+1+1)),range(1,l+1))
    for i in range(1,l+1):
        for i,j in zip([i,i+1],[1,-1]):
            r[i-1][i-1]=j
    return r

ChevieData["A"]["GeneratingRoots"]=weyla3

def weyla4(l,s):
    return ChevieData["imp"]["ParabolicRepresentatives"](1,1,l,s)

ChevieData["A"]["ParabolicRepresentatives"]=weyla4

def weyla5(pi):
    w=[]
    i=0
    for l in pi:
        r=l%2
        w+=i+Concatenation(range(1,l-1-r+1,3-1),range(2,l+r-2+1,4-2))
        i=i+l
    return w

ChevieData["A"]["WordClass"]=weyla5

def weyla6(n):
    res={"classparams":Partitions(n+1)}
    res["classnames"]=map(IntListToString,res["classparams"])
    res["classtext"]=map(ChevieData["A"]["WordClass"],res["classparams"])
    res["classes"]=map(lambda pi: GAPDiv(Factorial(n+1),CharTableSymmetric["centralizers"][1-1](n,pi)),res["classparams"])
    res["orders"]=map(Lcm,res["classparams"])
    return res

ChevieData["A"]["ClassInfo"]=weyla6

ChevieData["A"]["NrConjugacyClasses"]=lambda n: NrPartitions(n+1)

ChevieData["A"]["WeightInfo"]=lambda n: {"minusculeWeights":range(1,n+1),
    "decompositions":map(lambda i: [i],range(1,n+1)),
    "moduli":[n+1]}

def weyla7(n,w):
    x=Permutation("()")
    for i in w:
        x=GAPMul(x,Permutation("(%s,%s)"%(i,i+1)))
    res=[]
    mark=range(1,n+1+1)
    for i in range(1,n+1+1):
        if mark[i-1]!=0 :
            cyc=CyclePermInt(x,i)
            res.append(len(cyc))
            for i,j in zip(cyc,GAPMul(cyc,0)):
                mark[i-1]=j
    Sort(res)
    return Reversed(res)

ChevieData["A"]["ClassParameter"]=weyla7

ChevieData["A"]["CharParams"]=lambda n: Partitions(n+1)

ChevieData["A"]["LowestPowerFakeDegree"]=lambda p: GAPMul(p,range(0,len(p)-1+1))

ChevieData["A"]["HighestPowerFakeDegree"]=lambda p: GAPDiv(GAPMul(Sum(p),Sum(p)-1),2)-ChevieData["A"]["LowestPowerFakeDegree"](AssociatedPartition(p))

def weyla8(arg):
    return IntListToString(arg[2-1])

ChevieData["A"]["CharName"]=weyla8

def weyla9(n):
    res={"charparams":ChevieData["A"]["CharParams"](n)}
    res["extRefl"]=map(lambda i: res["charparams"].index(Concatenation([n+1-i],GAPMul(range(1,i+1),0)+1))+1,range(0,n+1))
    res["b"]=map(ChevieData["A"]["LowestPowerFakeDegree"],res["charparams"])
    res["B"]=map(ChevieData["A"]["HighestPowerFakeDegree"],res["charparams"])
    res["a"]=res["b"]
    res["A"]=res["B"]
    return res

ChevieData["A"]["CharInfo"]=weyla9

ChevieData["A"]["CharTable"]=ChevieData["compat"]["CharTableA"]

ChevieData["tmp"]=copy(CharTableSymmetric)

ChevieData["tmp"]["identifier"]="HeckeA"

ChevieData["tmp"]["specializedname"]=lambda nq: SPrint("H(A",nq[1-1]-1,")")

ChevieData["tmp"]["order"]=lambda nq: Factorial(nq[1-1])

ChevieData["tmp"]["size"]=lambda nq: Factorial(nq[1-1])

ChevieData["tmp"]["domain"]=lambda nq: IsList(nq) and len(nq)==2 and IsInt(nq[1-1]) and nq[1-1]>0

ChevieData["tmp"]["text"]="generic character table of Hecke algebras of type_ A"

ChevieData["tmp"]["classparam"]=[lambda nq: Partitions(nq[1-1])]

ChevieData["tmp"]["charparam"]=[lambda nq: Partitions(nq[1-1])]

def weyla10(nq,gamma,pi):
    n=nq[1-1]
    q=nq[2-1]
    if n==0 :
        return q**0
    k=pi[1-1]
    val=GAPMul(0,q)
    AHk=ChevieData["A"]["Hk"]["irreducibles"][1-1][1-1]
    for alpha in Partitions(n-k):
        dif=DifferencePartitions(gamma,alpha)
        if dif!=false :
            val=val+GAPMul(q-1**dif["cc"]-1,-1**dif["ll"])
    return val

ChevieData["tmp"]["irreducibles"]=[[weyla10]]

def weyla11(nq):
    n=nq[1-1]
    q=nq[2-1]
    pm=[]
    scheme=[]
    def hooks(beta,m):
        prs=[]
        for i in beta:
            leg=0
            for j in range(i-1,0+1,i-2-i-1):
                if j in beta :
                    leg=leg+1
                else:
                    prs.append([{"from":i,
                        "to":j,
                        "leg":leg}])
        cbs=copy(prs)
        hks=map(lambda x: [],range(1,m+1))
        for hk in cbs:
            for pr in prs:
                if pr[1-1]["to"]>hk[len(hk)-1]["from"] :
                    cbs.append(Concatenation(hk,pr))
            ll=Sum(hk,lambda x: x["from"]-x["to"])
            lg=Sum(hk,lambda x: x["leg"])
            lh=len(hk)
            new={"wgt":GAPMul(-1**lg,q**ll-lg-lh),
                "adr":1}
            if ll<m-1 :
                gamma=Difference(beta,map(lambda x: x["from"],hk))
                UniteSet(gamma,map(lambda x: x["to"],hk))
                if 0 in gamma :
                    j=0
                    while gamma[j+1-1]==j:
                        j=j+1
                    gamma=[gamma[k-1] for k in range(j+1,len(gamma)+1)]-j
                new["adr"]=pm[m-ll-1].index(gamma)+1
            hks[ll-1].append(new)
        return hks
    
    
    for i in range(1,n+1):
        pm[i-1]=map(BetaSet,Partitions(i))
        scheme[i-1]=[]
        for beta in pm[i-1]:
            scheme[i-1].append(hooks(beta,i))
    def charCol(n,t,k):
        col=[]
        for pi in scheme[n-1]:
            val=GAPMul(0,q)
            for hk in pi[k-1]:
                val=val+GAPMul(hk["wgt"],t[hk["adr"]-1])
            col.append(val)
        return col
    
    
    pm=map(lambda x: [],range(1,n-1+1))
    for m in range(1,QuoInt(n,2)+1):
        pm[m-1].append(charCol(m,[1],m))
        for k in range(m+1,n-m+1):
            for t in pm[k-m-1]:
                pm[k-1].append(charCol(k,t,m))
    res=[]
    for k in range(1,n-1+1):
        for t in pm[n-k-1]:
            res.append(charCol(n,t,k))
    res.append(charCol(n,[1],n))
    return Matrix(res).transpose()

ChevieData["tmp"]["matrix"]=weyla11

ChevieData["A"]["Hk"]=copy(ChevieData["tmp"])

ChevieData["A"]["HeckeCharTable"]=ChevieData["compat"]["HeckeCharTableA"]

def weyla12(n,param):
    return prod(range(1,n+1))

ChevieData["A"]["PoincarePolynomial"]=weyla12

def weyla13(n,alpha,param,sqrtparam):
    q=GAPDiv(-param[1-1][1-1],param[1-1][2-1])
    lambda_=BetaSet(alpha)
    res=q**Binomial(len(lambda_),3)
    for i in lambda_:
        for j in range(0,i-1+1):
            if j in lambda_ :
                res=GAPDiv(res,q**j)
            else:
                res=GAPMul(res,Sum(range(0,i-j-1+1),lambda e: q**e))
    return res

ChevieData["A"]["SchurElement"]=weyla13

def weyla14(arg):
    return ChevieData["imp"]["FactorizedSchurElement"](1,1,arg[1-1]+1,[arg[2-1]],arg[3-1],[])

ChevieData["A"]["FactorizedSchurElement"]=weyla14

def weyla15(n,param,sqrtparam,i):
    H=Hecke(CoxeterGroup("A",n),GAPDiv(-param[1-1][1-1],param[1-1][2-1]))
    return SpechtModel(H,Partitions(n+1)[i-1])

ChevieData["A"]["HeckeRepresentation"]=weyla15

def weyla16(n,i):
    return [ChevieData["imp"]["Representation"](1,1,n+1,i)[k-1] for k in range(2,n+1+1)]

ChevieData["A"]["Representation"]=weyla16

def weyla17(n,p,q):
    return GAPDiv(ChevieData["A"]["PoincarePolynomial"](Sum(p)-1,[[q,-1]]),ChevieData["A"]["SchurElement"](Sum(p)-1,p,[[q,-1]],[]))

ChevieData["A"]["FakeDegree"]=weyla17

def weyla18(l,p):
    return [[range(1,NrPartitions(l+1)+1),MatrixDecompositionMatrix(DecompositionMatrix(Specht(p,p),l+1))]]

ChevieData["A"]["DecompositionMatrix"]=weyla18

def weyla19(l):
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

ChevieData["A"]["UnipotentCharacters"]=weyla19

def weyla20(n):
    m=ChevieData["A"]["GeneratingRoots"](n)
    m.append(GAPMul(range(1,n+1+1),0)+1)
    return map(lambda i: weyla21,range(2,n+1+1))

def weyla21(arg):
    v=copy(arg)
    v.append(GAPMul(0,v[1-1]))
    v=GAPMul(v,m)
    return Sum(Arrangements(range(1,n+1+1),i),lambda a: prod([v[k-1] for k in a]))

ChevieData["A"]["Invariants"]=weyla20

def weyla22(n,p):
    uc={"classes":map(lambda p: {"parameter":p},Partitions(n+1)),
        "springerSeries":Concatenation(map(lambda d: map(lambda i: {"relgroup":CoxeterGroup("A",GAPDiv(n+1,d)-1),
        "Z":[ER(d)**i],
        "levi":Filtered(range(1,n+1+1),lambda i: i%d!=0),
        "locsys":[]},PrimeResidues(d)),DivisorsInt(n+1)))}
    ss=lambda z: First(uc["springerSeries"],lambda x: x["Z"]==[z])
    for i in range(1,len(uc["classes"])+1):
        cl=uc["classes"][i-1]
        p=cl["parameter"]
        d=Gcd(p)
        cl["name"]=IntListToString(p)
        cl["Au"]=ReflectionGroup(d,1,1)
        cl["balacarter"]=Concatenation(map(lambda i: Sum([p[k-1] for k in range(1,i-1+1)])+range(1,p[i-1]-1+1),range(1,len(p)+1)))
        p=Concatenation(map(lambda x: range(1-x,x-1+1,3-x-1-x),p))
        Sort(p)
        cl["dynkin"]=map(lambda i: p[i+1-1]-p[i-1],range(1,len(p)-1+1))
        cl["red"]=[]
        p=1
        for j in Collected(cl["parameter"]):
            cl["red"]+=range(p,p+j[2-1]-2+1)
            p=p+j[2-1]
        cl["red"]=ReflectionSubgroup(CoxeterGroup("A",p-2),cl["red"])
        cl["AuAction"]=ExtendedReflectionGroup(cl["red"],[IdentityMat(cl["red"]["rank"])])
        if d==2 :
            ss(1)["locsys"].append([i,2])
            ss(-1)["locsys"].append([i,1])
        else:
            for j in range(0,d-1+1):
                ss(ER(d)**j)["locsys"].append([i,j+1])
    uc["orderClasses"]=Hasse(Poset(map(lambda x: map(lambda y: Dominates(y["parameter"],x["parameter"]),uc["classes"]),uc["classes"])))
    return uc

ChevieData["A"]["UnipotentClasses"]=weyla22

def weyla23(n):
    def f(i):
        if i!=Permutation("()") :
            i=prod(CoxeterWord(W,i))
        i=map(Length,RobinsonSchenstedCorrespondent(n+1,i)["P"])
        return CharParams(W).index([i])+1
    
    
    W=CoxeterGroup("A",n)
    l=Filtered(Elements(W),lambda x: x**2==Permutation("()"))
    return map(lambda x: {"duflo":OnTuples(range(1,n+1),x),
        "reps":[],
        "character":[f(x)]},l)

ChevieData["A"]["KLeftCellRepresentatives"]=weyla23
