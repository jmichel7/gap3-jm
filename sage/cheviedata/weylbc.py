
def eylbc1(arg):
    n=arg[1-1]
    if len(arg)==2 :
        type_=arg[2-1]
    else:
        type_=2
    a=ChevieData["A"]["CartanMat"](n)
    a[1-1][2-1]=-type_
    a[2-1][1-1]=GAPDiv(2,a[1-1][2-1])
    return a

ChevieData["B"]["CartanMat"]=eylbc1

def eylbc2(r,indices,title,type_):
    print title," ",
    if type_==1 :
        print indices[1-1]," >=> ",indices[2-1],
    else:
        if type_==2 :
            print indices[1-1]," <=< ",indices[2-1],
        else:
            if type_==ER(2) :
                print indices[1-1]," == ",indices[2-1],
            else:
                print indices[1-1]," ?==? ",indices[2-1],
    for i in range(3,r+1):
        print " - ",indices[i-1],
    print "\n",

ChevieData["B"]["PrintDiagram"]=eylbc2

def eylbc3(arg):
    r=arg[1-1]
    option=arg[len(arg)-1]
    if len(arg)==3 :
        type_=arg[2-1]
    else:
        type_=2
    if type_==2 :
        if "TeX" in option :
            return SPrint("B_",TeXBracket(r))
        else:
            if "arg" in option :
                return SPrint("\"B\",",r)
            else:
                return SPrint("B",r)
    else:
        if type_==1 :
            if "TeX" in option :
                return SPrint("C_",TeXBracket(r))
            else:
                if "arg" in option :
                    return SPrint("\"C\",",r)
                else:
                    return SPrint("C",r)
        else:
            if type_==ER(2) :
                if "TeX" in option :
                    return SPrint("B^{\\hbox{sym}}_",TeXBracket(r))
                else:
                    if "arg" in option :
                        return SPrint("\"Bsym\",",r)
                    else:
                        return SPrint("Bsym",r)
            else:
                if "TeX" in option :
                    return SPrint("B^?_",TeXBracket(r),"(",Format(type_,option),")")
                else:
                    if "arg" in option :
                        return SPrint("\"B?\",",r,",",type_)
                    else:
                        return SPrint("B?",r,"(",Format(type_),")")

ChevieData["B"]["ReflectionName"]=eylbc3

def eylbc4(l,type_):
    rts=map(lambda i: GAPMul(0,range(1,l+1)),range(1,l+1))
    for i in range(1,l-1+1):
        for i,j in zip([i,i+1],[1,-1]):
            rts[i-1][i-1]=j
    rts[l-1][l-1]=GAPDiv(2,type_)
    return [rts[k-1] for k in range(l,1+1,l-1-l)]

ChevieData["B"]["GeneratingRoots"]=eylbc4

def eylbc5(l,s):
    return ChevieData["imp"]["ParabolicRepresentatives"](2,1,l,s)

ChevieData["B"]["ParabolicRepresentatives"]=eylbc5

ChevieData["B"]["ReflectionDegrees"]=lambda n: GAPMul(2,range(1,n+1))

def eylbc6(arg):
    return GAPMul(2**arg[1-1],Factorial(arg[1-1]))

ChevieData["B"]["Size"]=eylbc6

def eylbc7(arg):
    return NrPartitionTuples(arg[1-1],2)

ChevieData["B"]["NrConjugacyClasses"]=eylbc7

def eylbc8(n,type_):
    if type_==2 :
        return {"minusculeWeights":[1],
            "minusculeCoweights":[n],
            "decompositions":[[1]],
            "moduli":[2]}
    else:
        return {"minusculeWeights":[n],
            "minusculeCoweights":[1],
            "decompositions":[[1]],
            "moduli":[2]}

ChevieData["B"]["WeightInfo"]=eylbc8

def eylbc9(pi):
    w=[]
    i=1
    for l in Reversed(pi[2-1]):
        w+=range(i,2+1,i-1-i)
        w+=range(1,i+l-1+1)
        i=i+l
    for l in pi[1-1]:
        r=l%2
        w+=i+Concatenation(range(1,l-1-r+1,3-1),range(2,l+r-2+1,4-2))
        i=i+l
    return w

ChevieData["B"]["WordClass"]=eylbc9

def eylbc10(n):
    res=ChevieData["imp"]["ClassInfo"](2,1,n)
    res["classtext"]=map(ChevieData["B"]["WordClass"],res["classparams"])
    res["classes"]=map(lambda x: GAPDiv(res["centralizers"][1-1],x),res["centralizers"])
    return res

ChevieData["B"]["ClassInfo"]=eylbc10

def eylbc11(n,w):
    x=Permutation("()")
    for i in w:
        if i==1 :
            x=GAPMul(x,Permutation("(%s,%s)"%(1,n+1)))
        else:
            x=GAPMul(x,Permutation("(%s,%s)"%(i-1,i))(i-1+n,i+n))
    res=[[],[]]
    mark=range(1,n+1)
    for i in range(1,n+1):
        if mark[i-1]!=0 :
            cyc=CyclePermInt(x,i)
            if i+n in cyc :
                res[2-1].append(GAPDiv(len(cyc),2))
            else:
                res[1-1].append(len(cyc))
            for j in cyc:
                if j>n :
                    mark[j-n-1]=0
                else:
                    mark[j-1]=0
    Sort(res[1-1])
    Sort(res[2-1])
    return [Reversed(res[1-1]),Reversed(res[2-1])]

ChevieData["B"]["ClassParameter"]=eylbc11

ChevieData["B"]["CharParams"]=lambda n: PartitionTuples(n,2)

def eylbc12(arg):
    return PartitionTupleToString(arg[2-1])

ChevieData["B"]["CharName"]=eylbc12

def eylbc13(p):
    pp=SymbolPartitionTuple(p,1)
    m=len(pp[2-1])
    res=GAPMul(pp[1-1],range(m,0+1,m-1-m))
    if pp[2-1]!=[] :
        res=res+GAPMul(pp[2-1],range(m-1,0+1,m-2-m-1))
    return GAPMul(2,res)+Sum(pp[2-1])-GAPDiv(GAPMul(m,m-1),6)

ChevieData["B"]["LowestPowerFakeDegree"]=eylbc13

def eylbc14(n):
    res={"charparams":ChevieData["B"]["CharParams"](n)}
    res["extRefl"]=Concatenation(map(lambda i: res["charparams"].index([[n-i],GAPMul(range(1,i+1),0)+1])+1,range(0,n-1+1)),[res["charparams"].index([[],GAPMul(range(1,n+1),0)+1])+1])
    res["a"]=map(lambda p: LowestPowerGenericDegreeSymbol(SymbolPartitionTuple(p,1)),res["charparams"])
    res["A"]=map(lambda p: HighestPowerGenericDegreeSymbol(SymbolPartitionTuple(p,1)),res["charparams"])
    res["b"]=map(ChevieData["B"]["LowestPowerFakeDegree"],res["charparams"])
    res["B"]=res["a"]+res["A"]-res["b"]
    return res

ChevieData["B"]["CharInfo"]=eylbc14

for t in ["B"]:
    ChevieData[t]["CharTable"]=ChevieData["compat"]["CharTableB"]

ChevieData["tmp"]=copy(CharTableWeylB)

ChevieData["tmp"]["identifier"]="HeckeB"

ChevieData["tmp"]["specializedname"]=lambda nq: SPrint("H(B",nq[1-1],")")

ChevieData["tmp"]["order"]=lambda nq: GAPMul(2**nq[1-1],Factorial(nq[1-1]))

ChevieData["tmp"]["size"]=lambda nq: GAPMul(2**nq[1-1],Factorial(nq[1-1]))

def eylbc15(nq):
    return IsList(nq) and len(nq)==3 and IsInt(nq[1-1]) and nq[1-1]>0

ChevieData["tmp"]["domain"]=eylbc15

ChevieData["tmp"]["text"]="generic character table of Hecke algebras of type_ B"

ChevieData["tmp"]["classparam"]=[lambda nq: PartitionTuples(nq[1-1],2)]

ChevieData["tmp"]["charparam"]=[lambda nq: PartitionTuples(nq[1-1],2)]

def eylbc16(nq,gamma,pi):
    n=nq[1-1]
    q=nq[3-1]
    Q=nq[2-1]
    if n==0 :
        return q**0
    val=GAPMul(0,q)
    BHk=ChevieData["B"]["Hk"]["irreducibles"][1-1][1-1]
    if pi[1-1]!=[] :
        k=pi[1-1][1-1]
        for alpha in PartitionTuples(n-k,2):
            dif=[]
            dif[1-1]=DifferencePartitions(gamma[1-1],alpha[1-1])
            dif[2-1]=DifferencePartitions(gamma[2-1],alpha[2-1])
            if dif[1-1]!=false and dif[2-1]!=false :
                dif={"cc":dif[1-1]["cc"]+dif[2-1]["cc"],
                    "ll":dif[1-1]["ll"]+dif[2-1]["ll"]}
                val=val+GAPMul(q-1**dif["cc"]-1,-1**dif["ll"])
    else:
        k=pi[2-1][1-1]
        nn=Sum(gamma[1-1])
        if nn>=k :
            for alpha in Partitions(nn-k):
                dif=DifferencePartitions(gamma[1-1],alpha)
                if dif!=false and dif["cc"]==1 :
                    val=val+GAPMul(Q,-1**dif["ll"])
        nn=Sum(gamma[2-1])
        if nn>=k :
            for alpha in Partitions(nn-k):
                dif=DifferencePartitions(gamma[2-1],alpha)
                if dif!=false and dif["cc"]==1 :
                    val=val+GAPMul(-1**dif["ll"]+1,q**n+dif["d"])
    return val

ChevieData["tmp"]["irreducibles"]=[[eylbc16]]

def eylbc17(nq):
    def DoublePartitions(n):
        if n==0 :
            return [[[],[]]]
        pm=map(lambda x: [],range(1,n+1))
        for m in range(1,n+1):
            pm[m-1].append([[],[m]])
            for k in range(m+1,n+1):
                for t in pm[k-m-1]:
                    s=[[],[m]]
                    s[2-1]+=t[2-1]
                    pm[k-1].append(s)
        for m in range(1,QuoInt(n,2)+1):
            pm[m-1].append([[m],[]])
            for k in range(m+1,n-m+1):
                for t in pm[k-m-1]:
                    s=[[m],t[2-1]]
                    s[1-1]+=t[1-1]
                    pm[k-1].append(s)
        res=[]
        for k in range(1,n-1+1):
            for t in pm[n-k-1]:
                s=[[k],t[2-1]]
                s[1-1]+=t[1-1]
                res.append(s)
        res.append([[n],[]])
        res+=pm[n-1]
        return res
    
    
    n=nq[1-1]
    x=nq[3-1]
    y=nq[2-1]
    pm=[]
    scheme=[]
    def hooks(beta,m):
        hks=map(lambda x: [],range(1,m+1))
        prs=[]
        lb=[len(beta[1-1]),len(beta[2-1])]
        for i in [1,2]:
            prs[i-1]=[]
            for j in beta[i-1]:
                leg=0
                for k in Reversed(range(0,j-1+1)):
                    if k in beta[i-1] :
                        leg=leg+1
                    else:
                        prs[i-1].append({"from":j,
                            "to":k,
                            "leg":leg,
                            "pow":m+k-lb[i-1]})
        cbs=map(lambda x: [[x],[]],prs[1-1])
        cbs+=map(lambda x: [[],[x]],prs[2-1])
        for hk in cbs:
            for pr in prs[1-1]:
                if hk[2-1]==[] and pr["to"]>hk[1-1][len(hk[1-1])-1]["from"] :
                    new=map(copy,hk)
                    new[1-1].append(pr)
                    cbs.append(new)
            for pr in prs[2-1]:
                if hk[2-1]==[] or pr["to"]>hk[2-1][len(hk[2-1])-1]["from"] :
                    new=map(copy,hk)
                    new[2-1].append(pr)
                    cbs.append(new)
            ll=Sum(hk[1-1],lambda x: x["from"]-x["to"])+Sum(hk[2-1],lambda x: x["from"]-x["to"])
            lg=Sum(hk[1-1],lambda x: x["leg"])+Sum(hk[2-1],lambda x: x["leg"])
            lh=len(hk[1-1])+len(hk[2-1])
            new={"wgt":[GAPMul(-1**lg,x**ll-lg-lh),0],
                "adr":1}
            if lh==1 :
                if hk[1-1][1]==None :
                    new["wgt"][2-1]=GAPMul(-1**lg,y)
                else:
                    new["wgt"][2-1]=GAPMul(-1**lg+1,x**hk[2-1][1-1]["pow"])
            if ll<m :
                gamma=[]
                for i in [1,2]:
                    gamma[i-1]=Difference(beta[i-1],map(lambda x: x["from"],hk[i-1]))
                    UniteSet(gamma[i-1],map(lambda x: x["to"],hk[i-1]))
                    if 0 in gamma[i-1] :
                        j=0
                        while j<len(gamma[i-1]) and gamma[i-1][j+1-1]==j:
                            j=j+1
                        gamma[i-1]=[gamma[i-1][k-1] for k in range(j+1,len(gamma[i-1])+1)]-j
                new["adr"]=pm[m-ll-1].index(gamma)+1
            hks[ll-1].append(new)
        return hks
    
    
    for i in range(1,n+1):
        pm[i-1]=map(lambda p: map(BetaSet,p),PartitionTuples(i,2))
        scheme[i-1]=[]
        for beta in pm[i-1]:
            scheme[i-1].append(hooks(beta,i))
    def charCol(n,t,k,p):
        col=[]
        for pi in scheme[n-1]:
            val=GAPMul(0,y)
            for hk in pi[k-1]:
                val=val+GAPMul(hk["wgt"][p-1],t[hk["adr"]-1])
            col.append(val)
        return col
    
    
    pm=map(lambda x: [],range(1,n+1))
    for m in range(1,n+1):
        pm[m-1].append(charCol(m,[1],m,2))
        for k in range(m+1,n+1):
            for t in pm[k-m-1]:
                pm[k-1].append(charCol(k,t,m,2))
    for m in range(1,QuoInt(n,2)+1):
        pm[m-1].append(charCol(m,[1],m,1))
        for k in range(m+1,n-m+1):
            for t in pm[k-m-1]:
                pm[k-1].append(charCol(k,t,m,1))
    res=[]
    for k in range(1,n-1+1):
        for t in pm[n-k-1]:
            res.append(charCol(n,t,k,1))
    res.append(charCol(n,[1],n,1))
    res+=pm[n-1]
    res=Permuted(res,GAPDiv(Sortex(DoublePartitions(n)),Sortex(PartitionTuples(n,2))))
    return Matrix(res).transpose()

ChevieData["tmp"]["matrix"]=eylbc17

ChevieData["B"]["Hk"]=ChevieData["tmp"]

del ChevieData["tmp"]

for t in ["B"]:
    ChevieData[t]["HeckeCharTable"]=ChevieData["compat"]["HeckeCharTableB"]

def eylbc18(n,para):
    q1=GAPDiv(-para[1-1][1-1],para[1-1][2-1])
    q2=GAPDiv(-para[2-1][1-1],para[2-1][2-1])
    return prod(range(0,n-1+1))

ChevieData["B"]["PoincarePolynomial"]=eylbc18

def eylbc19(arg):
    return ChevieData["imp"]["SchurElement"](2,1,arg[1-1],arg[2-1],arg[3-1],[])

ChevieData["B"]["SchurElement"]=eylbc19

def eylbc20(arg):
    return ChevieData["imp"]["FactorizedSchurElement"](2,1,arg[1-1],arg[2-1],arg[3-1],[])

ChevieData["B"]["FactorizedSchurElement"]=eylbc20

def eylbc21(arg):
    return ChevieData["imp"]["HeckeRepresentation"](2,1,arg[1-1],arg[2-1],[],arg[4-1])

ChevieData["B"]["HeckeRepresentation"]=eylbc21

def eylbc22(n,i):
    return ChevieData["imp"]["Representation"](2,1,n,i)

ChevieData["B"]["Representation"]=eylbc22

def eylbc23(n,c,q):
    return Value(CycPolFakeDegreeSymbol(SymbolPartitionTuple(c,1)),q)

ChevieData["B"]["FakeDegree"]=eylbc23

def eylbc24(l,p):
    decS=lambda i: MatrixDecompositionMatrix(DecompositionMatrix(Specht(p,p),i))
    pp=map(Partitions,range(0,l+1))
    pt=PartitionTuples(l,2)
    if p==2 :
        return [[range(1,len(pt)+1),GAPMul(map(eylbc25,pt),decS(l))]]
    else:
        dd=Concatenation([[[1]],[[1]]],map(decS,range(2,l+1)))
        return map(lambda i: [map(lambda x: pt.index(x)+1,Cartesian(pp[i+1-1],pp[l+1-i-1])),map(lambda x: map(Product,Cartesian(x)),Cartesian(dd[i+1-1],dd[l+1-i-1]))],range(0,l+1))

def eylbc25(p):
    p=LittlewoodRichardsonRule(p[1-1],p[2-1])
    return map(eylbc26,pp[l+1-1])

ChevieData["B"]["DecompositionMatrix"]=eylbc24

def eylbc27(arg):
    rank=arg[1-1]
    uc={"harishChandra":[],
        "charSymbols":[]}
    for d in 1+GAPMul(2,range(0,QuoInt(-1+RootInt(1+GAPMul(4,rank),2),2)+1)):
        s=GAPDiv(d**2-1,4)
        s={"relativeType":{"series":"B",
            "indices":range(1+s,rank+1),
            "rank":rank-s},
            "levi":range(1,s+1),
            "eigenvalue":-1**QuoInt(d+1,4),
            "parameterExponents":Concatenation([d],GAPMul(range(2+s,rank+1),0)+1),
            "cuspidalName":SPrint("B_",TeXBracket(s))}
        uc["harishChandra"].append(s)
        symbols=Symbols(rank,d)
        s["charNumbers"]=range(1,len(symbols)+1)+len(uc["charSymbols"])
        FixRelativeType(s)
        uc["charSymbols"]+=symbols
    uc["harishChandra"][1-1]["cuspidalName"]=""
    uc["a"]=map(LowestPowerGenericDegreeSymbol,uc["charSymbols"])
    uc["A"]=map(HighestPowerGenericDegreeSymbol,uc["charSymbols"])
    uc["families"]=FamiliesClassical(uc["charSymbols"])
    if len(arg)==2 and arg[2-1]==1 :
        uc["harishChandra"][1-1]["relativeType"]["cartanType"]=1
    return uc

ChevieData["B"]["UnipotentCharacters"]=eylbc27

def eylbc28(r,type_,char):
    def part2dynkin(part):
        p=Concatenation(map(lambda d: range(1-d,d-1+1,3-d-1-d),part))
        Sort(p)
        p=[p[k-1] for k in range(QuoInt(3+len(p),2),len(p)+1)]
        if type_==1 :
            res=[GAPMul(2,p[1-1])]
        else:
            res=[p[1-1]]
        res+=[p[k-1] for k in range(2,len(p)+1)]-[p[k-1] for k in range(1,len(p)-1+1)]
        return res
    
    
    def addSpringer(s):
        ss=First(uc["springerSeries"],lambda x: x["defect"]==DefectSymbol(s["symbol"]))
        if s["sp"]==[[],[]] :
            p=1
        else:
            if s["sp"]==[[1],[]] :
                p=2
            else:
                if s["sp"]==[[],[1]] :
                    p=1
                else:
                    p=CharParams(ss["relgroup"]).index([s["sp"]])+1
        ss["locsys"][p-1]=[len(uc["classes"]),CharParams(cc["Au"]).index(map(eylbc29,s["Au"]))+1]
    
    
    if type_==ER(2) :
        type_=2
        char=2
    if char==2 :
        ss=XSP(4,2,r)
    else:
        if type_==1 :
            ss=XSP(2,1,r)
        else:
            ss=XSP(2,0,r)
    l=Union(map(lambda c: map(lambda x: [DefectSymbol(x["symbol"]),Sum(x["sp"],Sum)],c),ss))
    SortBy(l,lambda x: [AbsInt(x[1-1]),-SignInt(x[1-1])])
    uc={"classes":[],
        "springerSeries":map(eylbc30,l)}
    if char!=2 :
        def symbol2para(S):
            c=Concatenation(S)
            Sort(c)
            i=1
            part=[]
            d=type_%2
            while i<=len(c):
                if i==len(c) or c[i+1-1]-c[i-1]>0 :
                    part.append(GAPMul(2,c[i-1]-i-1)+1-d)
                    i=i+1
                else:
                    l=GAPMul(2,c[i-1]-i-1)-d
                    part+=[l,l]
                    i=i+2
            Sort(part)
            part=Filtered(part,lambda y: y!=0)
            return Reversed(part)
        
        
    else:
        def symbol2para(S):
            c=Concatenation(S)
            Sort(c)
            i=1
            part=[]
            ex=[]
            while i<=len(c):
                if i==len(c) or c[i+1-1]-c[i-1]>1 :
                    part.append(GAPMul(2,c[i-1]-GAPMul(2,i-1)))
                    i=i+1
                else:
                    if c[i-1]==c[i+1-1] :
                        l=GAPMul(2,c[i-1]-GAPMul(2,i-1))-2
                        part+=[l,l]
                        ex.append(l)
                        i=i+2
                    else:
                        if c[i-1]+1==c[i+1-1] :
                            l=GAPMul(2,c[i-1]-GAPMul(2,i-1))-1
                            part+=[l,l]
                            i=i+2
            Sort(part)
            part=Filtered(part,lambda y: y!=0)
            return [Reversed(part),ex]
        
        
    if char==2 :
        type_=1
    for cl in ss:
        cc={"parameter":symbol2para(cl[1-1]["symbol"])}
        cc["Au"]=CoxeterGroup(*Concatenation(map(lambda x: ["A",1],cl[1-1]["Au"])))
        if char!=2 :
            cc["dynkin"]=part2dynkin(cc["parameter"])
            cc["name"]=IntListToString(cc["parameter"])
        else:
            type_=1
            cc["dimBu"]=cl[1-1]["dimBu"]
            cc["name"]=Join(map(eylbc31,Reversed(Collected(cc["parameter"][1-1]))),"")
        cc["red"]=CoxeterGroup()
        if char==2 :
            j=cc["parameter"][1-1]
        else:
            j=cc["parameter"]
        for j in Collected(j):
            if j[1-1]%2==type_%2 :
                cc["red"]=GAPMul(cc["red"],CoxeterGroup("C",GAPDiv(j[2-1],2)))
            else:
                if j[2-1]%2!=0 :
                    if j[2-1]>1 :
                        cc["red"]=GAPMul(cc["red"],CoxeterGroup("B",GAPDiv(j[2-1]-1,2)))
                else:
                    if j[2-1]>2 :
                        cc["red"]=GAPMul(cc["red"],CoxeterGroup("D",GAPDiv(j[2-1],2)))
                    else:
                        cc["red"]=GAPMul(cc["red"],Torus(1))
        uc["classes"].append(cc)
        for s in cl:
            addSpringer(s)
    uc["orderClasses"]=Hasse(Poset(map(lambda x: map(eylbc32,uc["classes"]),uc["classes"])))
    if char!=2 and type_==2 :
        def LuSpin(p):
            Sort(p)
            a=[]
            b=[]
            d=[0,1,0,-1]
            d=[d[k-1] for k in map(lambda x: 1+x%4,p)]
            i=1
            while i<=len(p):
                l=p[i-1]
                t=Sum([d[k-1] for k in range(1,i-1+1)])
                if 1==l%4 :
                    a.append(GAPDiv(l-1,4)-t)
                    i=i+1
                else:
                    if 3==l%4 :
                        b.append(GAPDiv(l-3,4)+t)
                        i=i+1
                    else:
                        j=i
                        while i<=len(p) and p[i-1]==l:
                            i=i+1
                        j=GAPMul(range(1,GAPDiv(i-j,2)+1),0)
                        a+=j+GAPDiv(l+l%4,4)-t
                        b+=j+GAPDiv(l-l%4,4)
            a=Filtered(a,lambda x: x!=0)
            a=Reversed(a)
            b=Filtered(b,lambda x: x!=0)
            b=Reversed(b)
            if Sum(d)>=1 :
                return [a,b]
            else:
                return [b,a]
        
        
        def addSpringer(f,i,s,k):
            ss=First(uc["springerSeries"],f)
            if s in [[[],[1]],[[],[]]] :
                p=1
            else:
                if s==[[1],[]] :
                    p=2
                else:
                    p=CharParams(ss["relgroup"]).index([s])+1
            ss["locsys"][p-1]=[i,k]
        
        
        def trspringer(i,old,new):
            for ss in uc["springerSeries"]:
                for c in ss["locsys"]:
                    if c[1-1]==i :
                        p=old.index(c[2-1])+1
                        if p!=false :
                            c[2-1]=new[p-1]
        
        
        d=0
        while GAPMul(4,d**2)-GAPMul(3,d)<=r:
            i=GAPMul(4,d**2)-GAPMul(3,d)
            if r-d%2==0 :
                l=Concatenation(range(1,i+1),range(i+2,r+1,i+4-i+2))
                uc["springerSeries"].append({"relgroup":CoxeterGroup("B",GAPDiv(r-i,2)),
                    "levi":l,
                    "Z":[-1],
                    "locsys":[]})
                i=GAPMul(4,d**2)+GAPMul(3,d)
                if i<=r and d!=0 :
                    l=Concatenation(range(1,i+1),range(i+2,r+1,i+4-i+2))
                    uc["springerSeries"].append({"relgroup":CoxeterGroup("B",GAPDiv(r-i,2)),
                        "levi":l,
                        "Z":[-1],
                        "locsys":[]})
            d=d+1
        l=Filtered(range(1,len(uc["classes"])+1),lambda i: ForAll(Collected(uc["classes"][i-1]["parameter"]),lambda c: c[1-1]%2==0 or c[2-1]==1))
        for i in l:
            cl=uc["classes"][i-1]
            s=LuSpin(cl["parameter"])
            if Size(cl["Au"])==1 :
                cl["Au"]=CoxeterGroup("A",1)
                trspringer(i,[1],[2])
                d=1
            else:
                if Size(cl["Au"])==4 :
                    cl["Au"]=CoxeterGroup("B",2)
                    trspringer(i,[1,2,3,4],[1,3,5,4])
                    d=2
                else:
                    Error("Au non-commutative of order ",GAPMul(Size(cl["Au"]),2),"  !  implemented")
            addSpringer(lambda ss: ss["Z"]==[-1] and ss["relgroup"]["rank"]==Sum(s,Sum),i,s,d)
    return uc

def eylbc29(x):
    if x :
        return [1,1]
    else:
        return [2]

def eylbc30(d):
    res={"relgroup":CoxeterGroup("C",d[2-1]),
        "defect":d[1-1],
        "locsys":[],
        "levi":range(1,r-d[2-1]+1)}
    if char==2 :
        res["Z"]=[1]
    else:
        if type_==1 :
            res["Z"]=[-1**r-d[2-1]]
        else:
            if IsInt(ER(GAPMul(2,r-d[2-1])+1)) :
                res["Z"]=[1]
            else:
                res["Z"]=[-1]
    return res

def eylbc31(x):
    res=IntListToString(GAPMul(range(1,x[2-1]+1),0)+x[1-1],"[]")
    if x[1-1] in cc["parameter"][2-1] :
        return SPrint("(",res,")")
    return res

def eylbc32(y):
    if char!=2 :
        return Dominates(y["parameter"],x["parameter"])
    m=Maximum(x["parameter"][1-1][1-1],y["parameter"][1-1][1-1])
    f=lambda x: map(lambda i: Sum(Filtered(x,lambda z: z<i))+GAPMul(i,Number(x,lambda z: z>=i)),range(1,m+1))
    fx=f(x["parameter"][1-1])
    fy=f(y["parameter"][1-1])
    for i in range(1,m+1):
        if fx[i-1]<fy[i-1] :
            return false
        else:
            if fx[i-1]==fy[i-1] and i in y["parameter"][2-1] :
                if i in Difference(x["parameter"][1-1],x["parameter"][2-1]) :
                    return false
                if i<m and fx[i+1-1]-fy[i+1-1]%2==1 :
                    return false
    return true

ChevieData["B"]["UnipotentClasses"]=eylbc28

def eylbc33(n,type_):
    m=GAPMul(range(1,n+1),0)+1
    m[1-1]=GAPDiv(2,type_)
    m=GAPMul(DiagonalMat(m),ChevieData["imp"]["GeneratingRoots"](2,1,n))
    return map(lambda f: eylbc34,ChevieData["imp"]["Invariants"](2,1,n))

def eylbc34(arg):
    return f(*GAPMul(arg,m))

ChevieData["B"]["Invariants"]=eylbc33
