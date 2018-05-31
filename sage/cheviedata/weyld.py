
def weyld1(n):
    if n<3 :
        m=3
    else:
        m=n
    a=ChevieData["A"]["CartanMat"](m)
    for i,j in zip(range(1,3+1),[[2,0,-1],[0,2,-1],[-1,-1,2]]):
        [a[k-1] for k in range(1,3+1)][i-1]=j
    return [[a[k-1] for k in range(1,n+1)][k-1] for k in range(1,n+1)]

ChevieData["D"]["CartanMat"]=weyld1

def weyld2(arg):
    return GAPMul(2**arg[1-1]-1,Factorial(arg[1-1]))

ChevieData["D"]["Size"]=weyld2

def weyld3(r,indices,title):
    print title," ",indices[1-1],"\n",
    s=just("",len(title)+1)
    print s," \\\n",s,"  ",indices[3-1],
    for i in range(4,r+1):
        print " - ",indices[i-1],
    print "\n",
    print s," /\n",s,indices[2-1],"\n",

ChevieData["D"]["PrintDiagram"]=weyld3

def weyld4(l):
    rts=[]
    for i in range(1,l-1+1):
        r=GAPMul(0,range(1,l+1))
        for i,j in zip([i,i+1],[1,-1]):
            r[i-1]=j
        rts.append(r)
    r=GAPMul(0,range(1,l+1))
    for i,j in zip([l-1,l],[1,1]):
        r[i-1]=j
    rts.append(r)
    return Reversed(rts)

ChevieData["D"]["GeneratingRoots"]=weyld4

def weyld5(n):
    if n%2==1 :
        return {"minusculeWeights":[1,2,n],
            "decompositions":[[1],[3],[2]],
            "moduli":[4]}
    else:
        return {"minusculeWeights":[1,2,n],
            "decompositions":[[1,0],[0,1],[1,1]],
            "moduli":[2,2]}

ChevieData["D"]["WeightInfo"]=weyld5

def weyld6(l,s):
    return ChevieData["imp"]["ParabolicRepresentatives"](2,2,l,s)

ChevieData["D"]["ParabolicRepresentatives"]=weyld6

def weyld7(arg):
    n=arg[1-1]
    if len(arg)==2 :
        param=map(lambda a: map(copy,a),arg[2-1])
    else:
        param=PartitionTuples(n,2)
    res=[]
    for pi in param:
        if pi[2-1]=='+' :
            pi[2-1]=[]
        if IsList(pi[2-1]) and len(pi[2-1])%2==0 :
            w=[]
            i=1
            for l in Reversed(pi[2-1]):
                if i==1 :
                    w+=range(2,i+l-1+1)
                else:
                    w+=range(i,3+1,i-1-i)
                    w+=range(1,i+l-1+1)
                i=i+l
            for l in pi[1-1]:
                r=l%2
                w+=i+Concatenation(range(1,l-1-r+1,3-1),range(2,l+r-2+1,4-2))
                i=i+l
            if w!=[] and w[1-1]==2 :
                w[1-1]=1
            if pi[2-1]==[] and ForAll(pi[1-1],lambda x: x%2==0) :
                res.append(w)
                w=copy(w)
                w[1-1]=2
            res.append(w)
    return res

ChevieData["D"]["WordsClassRepresentatives"]=weyld7

def weyld8(n):
    res=ChevieData["imp"]["ClassInfo"](2,2,n)
    res["classparams"]=map(weyld9,res["classparams"])
    res["classtext"]=ChevieData["D"]["WordsClassRepresentatives"](n,res["classparams"])
    return res

def weyld9(x):
    if len(x)==2 :
        return x
    if x[3-1]==0 :
        return [x[1-1],'+']
    else:
        return [x[1-1],'-']

ChevieData["D"]["ClassInfo"]=weyld8

def weyld10(n):
    if n%2==1 :
        return GAPDiv(NrPartitionTuples(n,2),2)
    else:
        return GAPDiv(NrPartitionTuples(n,2)+GAPMul(3,NrPartitions(GAPDiv(n,2))),2)

ChevieData["D"]["NrConjugacyClasses"]=weyld10

ChevieData["D"]["CharInfo"]=lambda n: ChevieData["imp"]["CharInfo"](2,2,n)

def weyld11(arg):
    return PartitionTupleToString(arg[2-1])

ChevieData["D"]["CharName"]=weyld11

ChevieData["D"]["gensMODA"]=[None,None,None,[[Permutation("(1,2)(7,8)"),Permutation("(3,4)(5,6)"),Permutation("(2,3)(6,7)"),Permutation("(3,5)(4,6)")],[[4],[None,None,2]],[[2],[1,None,1]]],None,[[Permutation("(1,2)(8,11)(12,14)(15,17)(16,18)(19,21)(22,25)(31,32)"),Permutation("(3,4)(5,6)(7,9)(10,13)(20,23)(24,26)(27,28)(29,30)"),Permutation("(2,3)(6,8)(9,12)(13,16)(17,20)(21,24)(25,27)(30,31)"),Permutation("(3,5)(4,6)(12,15)(14,17)(16,19)(18,21)(27,29)(28,30)"),Permutation("(5,7)(6,9)(8,12)(11,14)(19,22)(21,25)(24,27)(26,28)"),Permutation("(7,10)(9,13)(12,16)(14,18)(15,19)(17,21)(20,24)(23,26)")],[[16],[4,None,6],[1,None,None,None,5]],[[12],[2,None,6],[None,2,None,None,4]]],None,[[Permutation("(1,2)(8,11)(12,15)(16,20)(17,21)(22,26)(23,27)(28,33)(29,34)(30,35)(36,41)(37,42)(43,50)(44,51)(52,59)(60,68)(61,69)(70,77)(78,85)(79,86)(87,92)(88,93)(94,99)(95,100)(96,101)(102,106)(103,107)(108,112)(109,113)(114,117)(118,121)(127,128)"),Permutation("(3,4)(5,6)(7,9)(10,13)(14,18)(19,24)(25,31)(32,38)(39,45)(40,46)(47,53)(48,54)(49,55)(56,62)(57,63)(58,64)(65,71)(66,72)(67,73)(74,80)(75,81)(76,82)(83,89)(84,90)(91,97)(98,104)(105,110)(111,115)(116,119)(120,122)(123,124)(125,126)"),Permutation("(2,3)(6,8)(9,12)(13,17)(18,23)(20,25)(24,30)(26,32)(33,39)(34,40)(41,48)(42,49)(50,57)(51,58)(53,61)(59,67)(62,70)(68,76)(71,78)(72,79)(80,87)(81,88)(89,95)(90,96)(97,103)(99,105)(104,109)(106,111)(112,116)(117,120)(121,123)(126,127)"),Permutation("(3,5)(4,6)(12,16)(15,20)(17,22)(21,26)(23,29)(27,34)(30,37)(35,42)(39,47)(45,53)(48,56)(54,62)(57,65)(58,66)(63,71)(64,72)(67,75)(73,81)(76,84)(82,90)(87,94)(92,99)(95,102)(100,106)(103,108)(107,112)(109,114)(113,117)(123,125)(124,126)"),Permutation("(5,7)(6,9)(8,12)(11,15)(22,28)(26,33)(29,36)(32,39)(34,41)(37,44)(38,45)(40,48)(42,51)(46,54)(49,58)(55,64)(65,74)(71,80)(75,83)(78,87)(81,89)(84,91)(85,92)(88,95)(90,97)(93,100)(96,103)(101,107)(114,118)(117,121)(120,123)(122,124)"),Permutation("(7,10)(9,13)(12,17)(15,21)(16,22)(20,26)(25,32)(31,38)(36,43)(41,50)(44,52)(48,57)(51,59)(54,63)(56,65)(58,67)(62,71)(64,73)(66,75)(70,78)(72,81)(77,85)(79,88)(86,93)(91,98)(97,104)(103,109)(107,113)(108,114)(112,117)(116,120)(119,122)"),Permutation("(10,14)(13,18)(17,23)(21,27)(22,29)(26,34)(28,36)(32,40)(33,41)(38,46)(39,48)(45,54)(47,56)(52,60)(53,62)(59,68)(61,70)(67,76)(69,77)(73,82)(75,84)(81,90)(83,91)(88,96)(89,97)(93,101)(95,103)(100,107)(102,108)(106,112)(111,116)(115,119)"),Permutation("(14,19)(18,24)(23,30)(27,35)(29,37)(34,42)(36,44)(40,49)(41,51)(43,52)(46,55)(48,58)(50,59)(54,64)(56,66)(57,67)(62,72)(63,73)(65,75)(70,79)(71,81)(74,83)(77,86)(78,88)(80,89)(85,93)(87,95)(92,100)(94,102)(99,106)(105,111)(110,115)")],[[64],[16,None,24],[None,None,32],[4,None,None,None,20],[None,None,None,None,None,None,16]],[[56],[12,None,24],[6,None,28],[2,4,None,None,18],[1,None,3,None,None,None,14]]]]

def weyld12(n,w):
    x=Permutation("()")
    for i in w:
        if i==1 :
            x=GAPMul(x,Permutation("(%s,%s)"%(1,n+2))(2,n+1))
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
    if res[2-1]==[] and ForAll(res[1-1],lambda i: i%2==0) :
        if not ChevieData["D"]["gensMODA"][n]==None :
            tmp=CoxeterGroup("D",n)
            gens=PermCosetsSubgroup(tmp,ReflectionSubgroup(tmp,range(2,n+1)))
            tmp=ChevieData["D"]["ClassInfo"](n)
            tmp=[tmp["classtext"][k-1] for k in Filtered(range(1,len(tmp["classnames"])+1),lambda i: '+' in tmp["classnames"][i-1] or '-' in tmp["classnames"][i-1])]
            tmp=map(lambda a: CycleStructurePerm(prod([gens[k-1] for k in a])),tmp)
            ChevieData["D"]["gensMODA"][n-1]=[gens,[tmp[k-1] for k in GAPMul(2,range(1,GAPDiv(len(tmp),2)+1))-1],[tmp[k-1] for k in GAPMul(2,range(1,GAPDiv(len(tmp),2)+1))]]
        tmp=CycleStructurePerm(prod([ChevieData["D"]["gensMODA"][n-1][1-1][k-1] for k in w]))
        if tmp in ChevieData["D"]["gensMODA"][n-1][2-1] and not tmp in ChevieData["D"]["gensMODA"][n-1][3-1] :
            res[2-1]='+'
        else:
            if not tmp in ChevieData["D"]["gensMODA"][n-1][2-1] and tmp in ChevieData["D"]["gensMODA"][n-1][3-1] :
                res[2-1]='-'
    Sort(res[1-1])
    if IsList(res[2-1]) :
        Sort(res[2-1])
        return [Reversed(res[1-1]),Reversed(res[2-1])]
    else:
        return [Reversed(res[1-1]),res[2-1]]

ChevieData["D"]["ClassParameter"]=weyld12

ChevieData["D"]["CharTable"]=ChevieData["compat"]["CharTableD"]

ChevieData["tmp"]=copy(CharTableWeylD)

ChevieData["tmp"]["identifier"]="HeckeD"

ChevieData["tmp"]["specializedname"]=lambda nq: SPrint("H(D",nq[1-1],")")

ChevieData["tmp"]["size"]=lambda nq: GAPMul(2**nq[1-1]-1,Factorial(nq[1-1]))

ChevieData["tmp"]["order"]=lambda nq: GAPMul(2**nq[1-1]-1,Factorial(nq[1-1]))

def weyld13(nq):
    return IsList(nq) and len(nq)==2 and IsInt(nq[1-1]) and nq[1-1]>1

ChevieData["tmp"]["domain"]=weyld13

ChevieData["tmp"]["text"]="generic character table of Hecke algebras of type_ D"

ChevieData["tmp"]["classparam"]=[lambda nq: CharTableWeylD["classparam"][1-1](nq[1-1])]

ChevieData["tmp"]["charparam"]=[lambda nq: ChevieData["D"]["CharInfo"](nq[1-1])["charparams"]]

def weyld14(nq,alpha,pi):
    n=nq[1-1]
    q=nq[2-1]
    s="+-"
    if q==1 :
        return CharTableWeylD["irreducibles"][1-1][1-1](n,alpha,pi)
    AHk=ChevieData["A"]["Hk"]["irreducibles"][1-1][1-1]
    BHk=ChevieData["B"]["Hk"]["irreducibles"][1-1][1-1]
    if not IsList(alpha[2-1]) :
        delta=[alpha[1-1],alpha[1-1]]
        if not IsList(pi[2-1]) :
            vb=GAPDiv(BHk([n,1,q],delta,[pi[1-1],[]]),2)
            va=GAPMul(GAPDiv(q+1**len(pi[1-1]),2),AHk([GAPDiv(n,2),q**2],alpha[1-1],GAPDiv(pi[1-1],2)))
            if s[alpha[3-1]+1-1]==pi[2-1] :
                val=vb+va
            else:
                val=vb-va
        else:
            val=GAPDiv(BHk([n,1,q],delta,pi),2)
    else:
        if not IsList(pi[2-1]) :
            val=BHk([n,1,q],alpha,[pi[1-1],[]])
        else:
            val=BHk([n,1,q],alpha,pi)
    return val

ChevieData["tmp"]["irreducibles"]=[[weyld14]]

ChevieData["D"]["Hk"]=copy(ChevieData["tmp"])

ChevieData["D"]["HeckeCharTable"]=ChevieData["compat"]["HeckeCharTableD"]

def weyld15(arg):
    p=arg[2-1]
    n=arg[1-1]
    if p[2-1] in "+-" :
        p=[p[1-1],p[1-1]]
    return ChevieData["imp"]["FactorizedSchurElement"](2,2,n,p,arg[3-1],[])

ChevieData["D"]["FactorizedSchurElement"]=weyld15

def weyld16(arg):
    i=arg[4-1]
    n=arg[1-1]
    p=ChevieData["D"]["CharInfo"](n)["charparams"][i-1]
    if p[len(p)-1]==0 :
        i=i+1
    else:
        if p[len(p)-1]==1 :
            i=i-1
    return ChevieData["imp"]["HeckeRepresentation"](2,2,n,arg[2-1],[],i)

ChevieData["D"]["HeckeRepresentation"]=weyld16

def weyld17(n,i):
    p=ChevieData["D"]["CharInfo"](n)["charparams"][i-1]
    if p[len(p)-1]==0 :
        i=i+1
    else:
        if p[len(p)-1]==1 :
            i=i-1
    return ChevieData["imp"]["Representation"](2,2,n,i)

ChevieData["D"]["Representation"]=weyld17

def weyld18(n,para):
    q=GAPDiv(-para[1-1][1-1],para[1-1][2-1])
    return GAPMul(Sum(range(0,n-1+1),lambda k: q**k),prod(range(1,n-1+1)))

ChevieData["D"]["PoincarePolynomial"]=weyld18

ChevieData["D"]["symbolcharparam"]=lambda c: SymbolPartitionTuple(c,0)

def weyld19(n):
    m=ChevieData["imp"]["GeneratingRoots"](2,2,n)
    return map(lambda f: weyld20,ChevieData["imp"]["Invariants"](2,2,n))

def weyld20(arg):
    return f(*GAPMul(arg,m))

ChevieData["D"]["Invariants"]=weyld19

ChevieData["D"]["CycPolGenericDegree"]=lambda c: CycPolGenericDegreeSymbol(SymbolPartitionTuple(c,0))

def weyld21(n,phi,q,sqrtparam):
    return GAPDiv(ChevieData["D"]["PoincarePolynomial"](n,q),Value(ChevieData["D"]["CycPolGenericDegree"](phi),GAPDiv(-q[1-1][1-1],q[1-1][2-1])))

ChevieData["D"]["SchurElement"]=weyld21

def weyld22(n,c,q):
    return Value(CycPolFakeDegreeSymbol(SymbolPartitionTuple(c,0)),q)

ChevieData["D"]["FakeDegree"]=weyld22

def weyld23(rank):
    uc={"harishChandra":[],
        "charSymbols":[]}
    for d in GAPMul(4,range(0,RootInt(QuoInt(rank,4),2)+1)):
        r=GAPDiv(d**2,4)
        s={"relativeType":{"series":"B",
            "indices":range(1+r,rank+1),
            "rank":rank-r},
            "levi":range(1,r+1),
            "eigenvalue":-1**QuoInt(d+1,4),
            "parameterExponents":Concatenation([d],GAPMul(range(2+r,rank+1),0)+1)}
        if r<10 :
            s["cuspidalName"]=SPrint("D_",r,"")
        else:
            s["cuspidalName"]=SPrint("D_{",r,"}")
        if d==0 :
            s["relativeType"]["series"]="D"
            s["cuspidalName"]=""
            s["parameterExponents"][1-1]=1
        uc["harishChandra"].append(s)
        symbols=Symbols(rank,d)
        s["charNumbers"]=range(1,len(symbols)+1)+len(uc["charSymbols"])
        FixRelativeType(s)
        uc["charSymbols"]+=symbols
    uc["a"]=map(LowestPowerGenericDegreeSymbol,uc["charSymbols"])
    uc["A"]=map(HighestPowerGenericDegreeSymbol,uc["charSymbols"])
    uc["families"]=FamiliesClassical(uc["charSymbols"])
    return uc

ChevieData["D"]["UnipotentCharacters"]=weyld23

ChevieData["D"]["ReflectionDegrees"]=lambda n: Concatenation(GAPMul(2,range(1,n-1+1)),[n])

def weyld24(n,char):
    def addSpringer(s,i):
        ss=First(uc["springerSeries"],lambda x: x["defect"]==DefectSymbol(s["symbol"]))
        if s["sp"] in [[[],[1]],[[],[]]] :
            p=1
        else:
            if s["sp"]==[[1],[]] :
                p=2
            else:
                p=CharParams(ss["relgroup"]).index([s["sp"]])+1
        ss["locsys"][p-1]=[i,CharParams(cc["Au"]).index(map(weyld25,s["Au"]))+1]
    
    
    def partition2DR(part):
        p=Concatenation(map(lambda x: range(1-x,x-1+1,3-x-1-x),part))
        Sort(p)
        p=[p[k-1] for k in range(1+GAPDiv(len(p),2),len(p)+1)]
        return Concatenation([p[1-1]+p[2-1]],map(lambda i: p[i+1-1]-p[i-1],range(1,len(p)-1+1)))
    
    
    if char==2 :
        ss=XSP(4,0,n,1)
        def symbol2partition(S):
            c=Concatenation(S)
            Sort(c)
            i=1
            part=[]
            ex=[]
            while i<=len(c):
                if i==len(c) or c[i+1-1]-c[i-1]>1 :
                    part.append(GAPMul(2,c[i-1]-GAPMul(2,i-1))+2)
                    i=i+1
                else:
                    if c[i+1-1]-c[i-1]>0 :
                        l=GAPMul(2,c[i-1]-GAPMul(2,i-1))+1
                        part+=[l,l]
                        i=i+2
                    else:
                        l=GAPMul(2,c[i-1]-GAPMul(2,i-1))
                        part+=[l,l]
                        i=i+2
                        ex.append(l)
            Sort(part)
            part=Filtered(part,lambda y: y!=0)
            return [Reversed(part),ex]
        
        
    else:
        ss=XSP(2,0,n,1)
        def symbol2partition(S):
            c=Concatenation(S)
            Sort(c)
            i=1
            part=[]
            while i<=len(c):
                if i==len(c) or c[i+1-1]-c[i-1]>0 :
                    part.append(GAPMul(2,c[i-1]-i-1)+1)
                    i=i+1
                else:
                    l=GAPMul(2,c[i-1]-i-1)
                    part+=[l,l]
                    i=i+2
            Sort(part)
            part=Filtered(part,lambda y: y!=0)
            return Reversed(part)
        
        
    l=Union(map(lambda c: map(lambda x: [DefectSymbol(x["symbol"]),Sum(FullSymbol(x["sp"]),Sum)],c),ss))
    SortBy(l,lambda x: [AbsInt(x[1-1]),-SignInt(x[1-1])])
    uc={"classes":[],
        "springerSeries":map(weyld26,l)}
    for cl in ss:
        cc={"parameter":symbol2partition(cl[1-1]["symbol"])}
        if char==2 :
            cc["dimBu"]=cl[1-1]["dimBu"]
            cc["name"]=Join(map(weyld27,Reversed(Collected(cc["parameter"][1-1]))),"")
        else:
            cc["dynkin"]=partition2DR(cc["parameter"])
            cc["name"]=IntListToString(cc["parameter"])
        cc["Au"]=CoxeterGroup(*Concatenation(map(lambda x: ["A",1],cl[1-1]["Au"])))
        CharNames(cc["Au"])
        if char!=2 :
            cc["red"]=CoxeterGroup()
            j=cc["parameter"]
            for j in Collected(j):
                if j[1-1]%2==0 :
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
        if not IsList(cl[1-1]["sp"][2-1]) :
            cl[1-1]["sp"][3-1]=1-GAPDiv(n,2)%2
        uc["classes"].append(cc)
        for s in cl:
            addSpringer(s,len(uc["classes"]))
        if not IsList(cl[1-1]["sp"][2-1]) :
            cl[1-1]["sp"][3-1]=1-cl[1-1]["sp"][3-1]
            cc["name"].append('+')
            cc=Copy(cc)
            cc["name"][len(cc["name"])-1]='-'
            if "dynkin" in cc :
                for i,j in zip([1,2],[cc["dynkin"][k-1] for k in [2,1]]):
                    cc["dynkin"][i-1]=j
            uc["classes"].append(cc)
            for s in cl:
                addSpringer(s,len(uc["classes"]))
    if char==2 :
        uc["orderClasses"]=Hasse(Poset(map(lambda x: map(weyld28,uc["classes"]),uc["classes"])))
    else:
        uc["orderClasses"]=Hasse(Poset(map(lambda i: map(lambda j: Dominates(uc["classes"][j-1]["parameter"],uc["classes"][i-1]["parameter"]) and uc["classes"][j-1]["parameter"]!=uc["classes"][i-1]["parameter"] or i==j,range(1,len(uc["classes"])+1)),range(1,len(uc["classes"])+1))))
    if char!=2 :
        d=0
        while GAPMul(4,d**2)-d<=n:
            i=GAPMul(4,d**2)-d
            if n-d%2==0 :
                l=Concatenation(range(1,i+1),range(i+2,n+1,i+4-i+2))
                s={"relgroup":CoxeterGroup("B",GAPDiv(n-i,2)),
                    "levi":l,
                    "locsys":[]}
                if n%2==0 :
                    s["Z"]=[1,-1]
                else:
                    s["Z"]=[ER(4)]
                uc["springerSeries"].append(s)
                if d==0 :
                    l=Concatenation([1],range(4,n+1,6-4))
                s={"relgroup":CoxeterGroup("B",GAPDiv(n-i,2)),
                    "levi":l,
                    "locsys":[]}
                if n%2==0 :
                    s["Z"]=[-1,1]
                else:
                    s["Z"]=[-ER(4)]
                uc["springerSeries"].append(s)
                i=GAPMul(4,d**2)+d
                if d!=0 and i<=n :
                    l=Concatenation(range(1,i+1),range(i+2,n+1,i+4-i+2))
                    s={"relgroup":CoxeterGroup("B",GAPDiv(n-i,2)),
                        "levi":l,
                        "locsys":[]}
                    if n%2==0 :
                        s["Z"]=[1,-1]
                    else:
                        s["Z"]=[ER(4)]
                    uc["springerSeries"].append(s)
                    s={"relgroup":CoxeterGroup("B",GAPDiv(n-i,2)),
                        "levi":l,
                        "locsys":[]}
                    if n%2==0 :
                        s["Z"]=[1,1]
                    else:
                        s["Z"]=[-ER(4)]
                    uc["springerSeries"].append(s)
            d=d+1
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
        
        
        def trspringer(i,new):
            for ss in uc["springerSeries"]:
                for c in ss["locsys"]:
                    if c[1-1]==i :
                        c[2-1]=new[c[2-1]-1]
        
        
        l=Filtered(range(1,len(uc["classes"])+1),lambda i: ForAll(Collected(uc["classes"][i-1]["parameter"]),lambda c: c[1-1]%2==0 or c[2-1]==1))
        for i in l:
            cl=uc["classes"][i-1]
            s=LuSpin(cl["parameter"])
            if Size(cl["Au"])==1 :
                cl["Au"]=CoxeterGroup("A",1)
                trspringer(i,[2])
                k=[1,1]
            else:
                if Size(cl["Au"])==2 :
                    cl["Au"]=CoxeterGroup("A",1,"A",1)
                    trspringer(i,[2,4])
                    k=[1,3]
                else:
                    if Size(cl["Au"])==8 :
                        cl["Au"]=CoxeterGroup("A",1,"B",2)
                        trspringer(i,[1,6,8,5,10,3,4,9])
                        k=[2,7]
                    else:
                        Error("Au non-commutative of order ",GAPMul(Size(cl["Au"]),2),"  !  implemented")
            if not '-' in cl["name"] :
                addSpringer(lambda ss: ss["Z"] in [[1,-1],[ER(4)]] and ss["relgroup"]["rank"]==Sum(s,Sum),i,s,k[1-1])
            if not '+' in cl["name"] :
                addSpringer(lambda ss: ss["Z"] in [[-1,1],[-ER(4)]] and ss["relgroup"]["rank"]==Sum(s,Sum),i,s,k[2-1])
    return uc

def weyld25(x):
    if x :
        return [1,1]
    else:
        return [2]

def weyld26(d):
    res={"defect":d[1-1],
        "locsys":[],
        "levi":range(1,n-d[2-1]+1)}
    if n-d[2-1]%4==0 or char==2 :
        if n%2==0 :
            res["Z"]=[1,1]
        else:
            res["Z"]=[1]
    else:
        if n%2==0 :
            res["Z"]=[-1,-1]
        else:
            res["Z"]=[-1]
    if d[1-1]==0 :
        res["relgroup"]=CoxeterGroup("D",d[2-1])
    else:
        res["relgroup"]=CoxeterGroup("B",d[2-1])
    return res

def weyld27(x):
    res=IntListToString(GAPMul(range(1,x[2-1]+1),0)+x[1-1],"[]")
    if x[1-1] in cc["parameter"][2-1] :
        return SPrint("(",res,")")
    return res

def weyld28(y):
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
    if x["parameter"]==y["parameter"] and x!=y :
        return false
    return true

ChevieData["D"]["UnipotentClasses"]=weyld24
