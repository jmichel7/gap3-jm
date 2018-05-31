
def eptio1(x,option):
    for f in ["frame","kondo","spaltenstein","gp"]:
        if f in option :
            s=ChevieData[t]["CharInfo"]()
            if f in s :
                s=s[f][s["charparams"].index(x)+1-1]
                if "TeX" in option :
                    return s
                else:
                    return TeXStrip(s)
    if "TeX" in option :
        s="\\phi_"
    else:
        s="phi"
    s+=SPrint("{",x[1-1],",",x[2-1],"}")
    if len(x)==3 :
        s+=map(lambda y: '\'',range(1,x[3-1]+1))
    return str(s)

for t in ["2E6","E6","E7","E8","2F4","F4","G2","H3","H4","2G5","G24","G25","G26","G27","G29","G31","G32","G33","G34"]:
    ChevieData[t]["CharName"]=lambda t: eptio1

def eptio2(t):
    r=ChevieData[t]["GeneratingRoots"]
    rbar=ComplexConjugate(r)
    e=ChevieData[t]["EigenvaluesGeneratingReflections"]
    e=1-map(lambda x: ER(Denominator(x))**Numerator(x),e)
    e=map(lambda i: GAPDiv(GAPMul(e[i-1],rbar[i-1]),GAPMul(rbar[i-1],r[i-1])),range(1,len(e)+1))
    return map(lambda x: map(lambda y: GAPMul(x,y),r),e)

for t in ["G25","G26","G29","G31","G32","G34"]:
    ChevieData[t]["CartanMat"]=eptio2

def eptio3(option):
    i=["G24","G25","G26","G27","G29","G31","G32","G33","G34","E6","E7","E8","2E6","2F4","3D4","H3","H4"]
    o=["G_{24}","G_{25}","G_{26}","G_{27}","G_{29}","G_{31}","G_{32}","G_{33}","G_{34}","E_6","E_7","E_8","{}^2E_6","{}^2F_4","{}^3D_4","H_3","H_4"]
    if "TeX" in option :
        return o[i.index(t)+1-1]
    else:
        return t

for t in ["G24","G25","G26","G27","G29","G31","G32","G33","G34","E6","E7","E8","2E6","2F4","3D4","H3","H4"]:
    ChevieData[t]["ReflectionName"]=lambda t: eptio3

def eptio4(r,option):
    i=["A","D","2A","2D"]
    o=["A","D","{}^2A","{}^2D"]
    if "arg" in option :
        return SPrint(FormatGAP(t),",",r)
    else:
        if "TeX" in option :
            return SPrint(o[i.index(t)+1-1],"_",TeXBracket(r))
        else:
            return SPrint(t,r)

for t in ["A","D","2A","2D"]:
    ChevieData[t]["ReflectionName"]=lambda t: eptio4

def eptio5():
    rank="12345678".index(t[len(t)-1])+1
    res=ChevieData[t]["HeckeCharTable"](map(lambda x: [1,-1],range(1,rank+1)),map(lambda x: 1,[1,rank]))
    ChevieData["compat"]["ChangeIdentifier"](res,SPrint("W(",t,")"))
    return res

for t in ["3D4","E6","2E6","E7","E8","F4","2F4","G2","H3","H4"]:
    ChevieData[t]["CharTable"]=lambda t: eptio5

def eptio6(q):
    return prod(ChevieData[t]["ReflectionDegrees"])

for t in ["G24","G27","G29","G33","G34","H3","H4","E6","E7","E8"]:
    ChevieData[t]["PoincarePolynomial"]=lambda t: eptio6

def eptio7(i):
    para=ChevieData[t]["EigenvaluesGeneratingReflections"]
    para=map(lambda x: map(lambda j: ER(GAPDiv(1,x))**j,range(0,GAPDiv(1,x)-1+1)),para)
    return ChevieData[t]["HeckeRepresentation"](para,[],i)

for t in ["G24","G25","G26","G27","G29"]:
    ChevieData[t]["Representation"]=lambda t: eptio7

def eptio8(t):
    r=ChevieData[t]["GeneratingRoots"]
    if callable(r) :
        r=r()
    return len(r[1-1])

for t in ["G2","F4","H3","E6","G24","G25","G26","G27","G29","G31","G32","G33","G34"]:
    ChevieData[t]["SemisimpleRank"]=eptio8

for t in ["A","B","D"]:
    ChevieData[t]["SemisimpleRank"]=lambda t: lambda r: r

def eptio9(phi,q):
    f=ChevieData[t]["sparseFakeDegrees"][ChevieData[t]["CharInfo"]()["charparams"].index(phi)+1-1]
    return Sum(range(1,len(f)-1+1,3-1),lambda i: GAPMul(f[i-1],q**f[i+1-1]))

for t in ["G2","F4","H3","E6","G24","G25","G26","G27","G29","G32","G33","G34"]:
    ChevieData[t]["FakeDegree"]=lambda t: eptio9

def eptio10(phi,q):
    f=ChevieData[t]["cycpolfakedegrees"][ChevieData[t]["CharInfo"]()["charparams"].index(phi)+1-1]
    if IsList(f[1-1]) :
        res=ValuePol(f[1-1],q**2)
    else:
        res=f[1-1]
    f=copy(f)
    f[1-1]=1
    return GAPMul(res,Value(CycPol(f),q))

for t in ["H4","E7","E8","G31"]:
    ChevieData[t]["FakeDegree"]=lambda t: eptio10

def eptio11():
    return map(eptio12,ChevieData[t]["cycpolfakedegrees"])

def eptio12(f):
    if IsList(f[1-1]) :
        res=GAPMul(2,len(f[1-1]))+f[2-1]-2
    else:
        res=f[2-1]
    return res+Sum([f[k-1] for k in range(3,len(f)+1)],Phi)

for t in ["H4","E7","E8","G31"]:
    ChevieData[t]["HighestPowerFakeDegrees"]=lambda t: eptio11

def eptio13():
    return map(lambda x: x[len(x)-1],ChevieData[t]["sparseFakeDegrees"])

for t in ["E6","G32","G33","G34","G2","F4","H3","G24","G25","G26","G27","G29"]:
    ChevieData[t]["HighestPowerFakeDegrees"]=lambda t: eptio13

def eptio14():
    return map(lambda x: x[2-1],ChevieData[t]["sparseFakeDegrees"])

for t in ["G2","F4","H3","H4","G24","G25","G26","G27","G29","E6","E7","E8","G31","G32","G33","G34"]:
    ChevieData[t]["LowestPowerFakeDegrees"]=lambda t: eptio14

def eptio15(indices,title):
    digits="678"
    print title," ",
    r=digits.index(t[2-1])+1+5
    l=len(str(indices[1-1]))+len(str(indices[3-1]))
    print just("",l-1),indices[2-1],"\n",
    print just("",len(title)+l),"|\n",
    print SPrint(just("",len(title)-2),indices[1-1]),
    for i in range(3,r+1):
        print " - ",indices[i-1],
    print "\n",

for t in ["E6","E7","E8"]:
    ChevieData[t]["PrintDiagram"]=lambda t: eptio15

def eptio16(indices,title):
    print title," ",
    print SPrint(just("",len(str(indices[1-1]))-1),"5 \n"),
    print just("",len(title)-1),indices[1-1]," - ",indices[2-1]," - ",indices[3-1],
    if t=="H4" :
        print " - ",indices[4-1],
    print "\n",

for t in ["H3","H4"]:
    ChevieData[t]["PrintDiagram"]=lambda t: eptio16

def eptio17():
    N=Sum(ChevieData[t]["ReflectionDegrees"],lambda x: x-1)
    return map(lambda x: N-Degree(CycPol(x)),ChevieData[t]["CycPolSchurElements"])

for t in ["G24","G27","G29","G33","G34","H3","H4","E6","E7","E8"]:
    ChevieData[t]["HighestPowerGenericDegrees"]=lambda t: eptio17

def eptio18():
    return map(lambda x: -x[2-1],ChevieData[t]["CycPolSchurElements"])

for t in ["G24","G27","G29","G33","G34","H3","H4","E6","E7","E8"]:
    ChevieData[t]["LowestPowerGenericDegrees"]=lambda t: eptio18

def eptio19(p):
    T=ChevieData[t]["CharTable"]()
    T["name"]=T["identifier"]
    m=DecompositionMatrix(T%p)
    return map(lambda c: [c[1-1],[[m[k-1] for k in c[1-1]][k-1] for k in c[2-1]]],BlocksMat(m))

for t in ["F4","G2","G25","G26"]:
    ChevieData[t]["DecompositionMatrix"]=lambda t: eptio19

def eptio20(arg):
    return Value(CycPol(ChevieData[t]["CycPolSchurElements"][ChevieData[t]["CharInfo"]()["charparams"].index(arg[1-1])+1-1]),GAPDiv(-arg[2-1][1-1][1-1],arg[2-1][1-1][2-1]))

for t in ["G24","G27","G29","G33","G34","E6","E7","E8","H3","H4"]:
    ChevieData[t]["SchurElement"]=lambda t: eptio20

def eptio21(arg):
    c=ChevieData[t]["CycPolSchurElements"][ChevieData[t]["CharInfo"]()["charparams"].index(arg[1-1])+1-1]
    q=GAPDiv(-arg[2-1][1-1][1-1],arg[2-1][1-1][2-1])
    res={"factor":Mvp(GAPMul(c[1-1],q**c[2-1])),
        "operations":FactorizedSchurElementsOps}
    res["vcyc"]=map(lambda v: {"monomial":q,
        "pol":CycPol([1,0,v])},[c[k-1] for k in range(3,len(c)+1)])
    return FactorizedSchurElementsOps["Simplify"](res)

for t in ["G24","G27","G29","G33","G34","E6","E7","E8","H3","H4"]:
    ChevieData[t]["FactorizedSchurElement"]=lambda t: eptio21

def eptio22(arg):
    Y=Concatenation([arg[2-1][k-1] for k in ChevieData[t]["HyperplaneRepresentatives"]])
    ci=ChevieData[t]["SchurData"][ChevieData[t]["CharInfo"]()["charparams"].index(arg[1-1])+1-1]
    return VFactorSchurElement(*Concatenation([Y,ChevieData[t]["SchurModels"][ci["name"]],ci],[arg[k-1] for k in range(3,len(arg)+1)]))

for t in ["G2","F4","G25","G26","G32"]:
    ChevieData[t]["FactorizedSchurElement"]=lambda t: eptio22

def eptio23(arg):
    Y=Concatenation([arg[2-1][k-1] for k in ChevieData[t]["HyperplaneRepresentatives"]])
    ci=ChevieData[t]["SchurData"][ChevieData[t]["CharInfo"]()["charparams"].index(arg[1-1])+1-1]
    return VcycSchurElement(Y,ChevieData[t]["SchurModels"][ci["name"]],ci)

for t in ["F4","G25","G26","G32"]:
    ChevieData[t]["SchurElement"]=lambda t: eptio23

def VcycSchurElement(arg):
    n=len(arg[1-1])
    if len(arg)==3 :
        data=arg[3-1]
        para=[arg[1-1][k-1] for k in data["order"]]
    else:
        para=copy(arg[1-1])
    monomial=lambda v: prod(range(1,len(v)+1))
    r=arg[2-1]
    if "coeff" in r :
        res=r["coeff"]
    else:
        res=1
    if "factor" in r :
        res=GAPMul(res,monomial(r["factor"]))
    if "root" in r :
        para=para+GAPMul(0,prod(para))
        para[n+1-1]=ChevieIndeterminate(para)
    else:
        if "rootUnity" in r :
            para[n+1-1]=r["rootUnity"]**data["rootUnityPower"]
    res=GAPMul(res,prod(r["vcyc"]))
    if "root" in r :
        den=Lcm(map(Denominator,r["root"]))
        root=monomial(GAPMul(den,r["root"]))
        if "rootCoeff" in r :
            root=GAPMul(root,r["rootCoeff"])
        return EvalPolRoot(res,root,den,data["rootPower"])
    else:
        return res



def VFactorSchurElement(arg):
    n=len(arg[1-1])
    if len(arg)>=3 :
        data=arg[3-1]
        para=[arg[1-1][k-1] for k in data["order"]]
    else:
        para=copy(arg[1-1])
    monomial=lambda v: prod(range(1,len(v)+1))
    r=arg[2-1]
    res={}
    if "coeff" in r :
        res["factor"]=r["coeff"]
    else:
        res["factor"]=1
    if "factor" in r :
        res["factor"]=GAPMul(res["factor"],monomial(r["factor"]))
    if "root" in r :
        den=Lcm(map(Denominator,r["root"]))
        root=monomial(GAPMul(r["root"],den))
        if "rootCoeff" in r :
            root=GAPMul(root,r["rootCoeff"])
        para[n+1-1]=GetRoot(root,den)
        if IsBound(data) :
            para[n+1-1]=GAPMul(para[n+1-1],data["rootPower"])
    else:
        if "rootUnity" in r :
            para[n+1-1]=r["rootUnity"]**data["rootUnityPower"]
    res["vcyc"]=map(lambda v: {"monomial":monomial(v[1-1]),
        "pol":CycPol([1,0,v[2-1]])},r["vcyc"])
    if res["factor"]==0 or res["vcyc"]==[] :
        return res["factor"]
    res["operations"]=FactorizedSchurElementsOps
    return FactorizedSchurElementsOps["Simplify"](res)


