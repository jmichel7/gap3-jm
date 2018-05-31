
def coxi1(arg):
    m=[[2,0],[0,2]]
    bond=arg[1-1]
    if bond==2 :
        return m
    if len(arg)==2 :
        type_=arg[2-1]
    else:
        if bond%2==0 :
            type_=1
        else:
            type_=ER(GAPMul(2,bond))+ER(GAPMul(2,bond))**-1
    m[1-1][2-1]=-type_
    m[2-1][1-1]=GAPDiv(2+ER(bond),m[1-1][2-1])
    return m

ChevieData["I"]["CartanMat"]=coxi1

def coxi2(arg):
    print arg[3-1]," ",
    bond=arg[1-1]
    indices=arg[2-1]
    if len(arg)==4 :
        type_=arg[4-1]
    else:
        type_=ER(GAPMul(2,bond))+ER(GAPMul(2,bond))**-1
    if type_==ER(GAPMul(2,bond))+ER(GAPMul(2,bond))**-1 :
        print indices[1-1]," -",str(bond),"- ",indices[2-1],"\n",
    else:
        if type_==1 :
            print indices[1-1]," >",str(bond),"> ",indices[2-1],"\n",
        else:
            print indices[1-1]," ?",str(bond),"? ",indices[2-1],"\n",

ChevieData["I"]["PrintDiagram"]=coxi2

def coxi3(arg):
    bond=arg[1-1]
    opt=arg[len(arg)-1]
    if len(arg)==3 :
        type_=arg[2-1]
    else:
        if bond%2==0 :
            type_=1
        else:
            type_=ER(GAPMul(2,bond))+ER(GAPMul(2,bond))**-1
    if type_==1 :
        if "TeX" in opt :
            return SPrint("I_2(",bond,")")
        else:
            if "arg" in opt :
                return SPrint("\"I\",2,",bond)
            else:
                return SPrint("I2(",bond,")")
    else:
        if type_==ER(GAPMul(2,bond))+ER(GAPMul(2,bond))**-1 :
            if bond%2==1 :
                if "TeX" in opt :
                    return SPrint("I_2(",bond,")")
                else:
                    if "arg" in opt :
                        return SPrint("\"I\",2,",bond)
                    else:
                        return SPrint("I2(",bond,")")
            else:
                if "TeX" in opt :
                    return SPrint("I_{\\hbox{sym}2}(",bond,")")
                else:
                    if "arg" in opt :
                        return SPrint("\"Isym\",2,",bond)
                    else:
                        return SPrint("Isym2(",bond,")")
        else:
            if "TeX" in opt :
                return SPrint("I_?(",Format(GAPDiv(type_**2,2+ER(bond)),opt),")(",bond,")")
            else:
                if "arg" in opt :
                    return SPrint("\"Isym\",2,",bond,",",Format(GAPDiv(type_**2,2+ER(bond)),opt))
                else:
                    return SPrint("I?(",GAPDiv(type_**2,2+ER(bond)),")(",bond,")")

ChevieData["I"]["ReflectionName"]=coxi3

ChevieData["I"]["SemisimpleRank"]=2

def coxi4(m):
    a=ER(GAPMul(2,m))**m-1
    b=ComplexConjugate(a)
    if m%2==0 :
        r=ER(GAPDiv(m,2))
    else:
        r=1
    return [[1,0],[GAPDiv(GAPMul(r,a+b),2),GAPDiv(GAPDiv(GAPMul(r,a-b),2),ER(4))]]

ChevieData["I"]["GeneratingRoots"]=coxi4

ChevieData["I"]["EigenvaluesGeneratingReflections"]=lambda m: [-1,-1]

def coxi5(arg):
    return GAPMul(2,arg[1-1])

ChevieData["I"]["Size"]=coxi5

ChevieData["I"]["ReflectionDegrees"]=lambda m: [2,m]

ChevieData["I"]["NrConjugacyClasses"]=lambda m: QuoInt(m+3,2)+GAPMul(m+1%2,2)

def coxi6(m,s):
    return ChevieData["imp"]["ParabolicRepresentatives"](m,m,2,s)

ChevieData["I"]["ParabolicRepresentatives"]=coxi6

def coxi7(m,x,option):
    if IsList(x[1-1]) :
        return PartitionTupleToString(x)
    else:
        if "TeX" in option :
            s="\\phi"
        else:
            s="phi"
        s=SPrint(s,"{",x[1-1],",",x[2-1],"}")
        if len(x)==3 :
            s+=x[3-1]
        return str(s)

ChevieData["I"]["CharName"]=coxi7

def coxi8(m):
    res={"charparams":[[1,0]]}
    if m%2==0 :
        res["extRefl"]=[1,5,4]
        res["charparams"]+=[[1,GAPDiv(m,2),"'"],[1,GAPDiv(m,2),"''"]]
    else:
        res["extRefl"]=[1,3,2]
    res["charparams"].append([1,m])
    res["charparams"]+=map(lambda i: [2,i],range(1,QuoInt(m-1,2)+1))
    res["b"]=map(lambda x: x[2-1],res["charparams"])
    res["B"]=map(coxi9,res["charparams"])
    res["a"]=map(coxi10,res["charparams"])
    res["A"]=map(coxi11,res["charparams"])
    res["charSymbols"]=map(coxi12,range(1,QuoInt(m-1,2)+1))
    v=map(lambda x: [0],range(1,m+1))
    v[m-1]=[1,2]
    res["charSymbols"]=Concatenation([v],res["charSymbols"])
    if m%2==0 :
        v=map(lambda x: [0],range(1,m+1))
        v[m-1]=[1]
        v[GAPDiv(m,2)-1]=[1]
        res["charSymbols"]=Concatenation([v],res["charSymbols"])
        v=map(lambda x: [0],range(1,m+1))
        v[m-1]=[1]
        v[GAPDiv(m,2)-1]=[1]
        res["charSymbols"]=Concatenation([v],res["charSymbols"])
    v=map(lambda x: [0,1],range(1,m+1))
    v[m-1]=[2]
    res["charSymbols"]=Concatenation([v],res["charSymbols"])
    res["malleParams"]=map(lambda x: map(PartBeta,x),res["charSymbols"])
    if m%2==0 :
        res["malleParams"][2-1]=Concatenation([res["malleParams"][2-1][k-1] for k in range(1,GAPDiv(m,2)+1)],[1])
        res["malleParams"][3-1]=Concatenation([res["malleParams"][3-1][k-1] for k in range(1,GAPDiv(m,2)+1)],[-1])
    return res

def coxi9(phi):
    if phi[1-1]==1 :
        return phi[2-1]
    else:
        return m-phi[2-1]

def coxi10(phi):
    if phi[1-1]!=1 or phi[2-1]==GAPDiv(m,2) :
        return 1
    else:
        return phi[2-1]

def coxi11(phi):
    if phi[1-1]==1 or phi[2-1]==GAPDiv(m,2) :
        return m-1
    else:
        return phi[2-1]

def coxi12(l):
    S=map(lambda i: [0],range(1,m+1))
    k=0
    if k!=0 :
        S[1-1]=[0,1]
        S[1+k+l%m-1]=[0,1]
        S[k+1-1]=[]
        S[l+1-1]=[]
    else:
        S[1-1]=[1]
        S[l+1-1]=[1]
    return S

ChevieData["I"]["CharInfo"]=coxi8

def coxi13(m):
    if IsInt(GAPDiv(m,2)) :
        r=[[],[1],[2]]
    else:
        r=[[],[1]]
    x=[1,2]
    for i in range(1,QuoInt(m,2)+1):
        r.append(copy(x))
        x+=[1,2]
    return r

ChevieData["I"]["WordsClassRepresentatives"]=coxi13

def coxi14(m):
    r=ChevieData["I"]["WordsClassRepresentatives"](m)
    clnp=map(IntListToString,r)
    g1=Permutation("()")
    i=2
    while GAPMul(2,i)<=m+1:
        g1=GAPMul(g1,Permutation("(%s,%s)"%(i,m-i+2)))
        i=i+1
    g2=Permutation("()")
    i=1
    while GAPMul(2,i)<=m:
        g2=GAPMul(g2,Permutation("(%s,%s)"%(i,m-i+1)))
        i=i+1
    gen=[g1,g2]
    def perm(l):
        if len(l)==0 :
            return Permutation("()")
        else:
            return prod([gen[k-1] for k in l])
    
    
    if m%2==0 :
        cl=[1,GAPDiv(m,2),GAPDiv(m,2)]
        cl+=GAPMul(range(1,GAPDiv(m,2)-1+1),0)+2
        cl.append(1)
    else:
        cl=[1,m]
        cl+=GAPMul(range(1,GAPDiv(m-1,2)+1),0)+2
    return {"classtext":r,
        "classnames":clnp,
        "classparams":clnp,
        "orders":map(lambda i: OrderPerm(perm(i)),r),
        "classes":cl}

ChevieData["I"]["ClassInfo"]=coxi14

def coxi15(m,param,rootparam):
    u=GAPDiv(-param[1-1][1-1],param[1-1][2-1])
    v=GAPDiv(-param[2-1][1-1],param[2-1][2-1])
    if m%2!=0 :
        squv=u
    else:
        if rootparam[1]==None and rootparam[2]==None :
            squv=GAPMul(rootparam[1-1],rootparam[2-1])
        else:
            squv=GetRoot(GAPMul(u,v),2,"CharTable(Hecke(I2(",m,")))")
    ct=[[u,v]]
    if m%2==0 :
        ct+=[[u,-u**0],[-v**0,v]]
    ct.append([-v**0,-v**0])
    cl=ChevieData["I"]["ClassInfo"](m)
    r=cl["classtext"]
    ct=map(lambda i: map(lambda x: prod([i[k-1] for k in x]),r),ct)
    ct+=map(coxi16,range(1,QuoInt(m-1,2)+1))
    tbl={"identifier":SPrint("H(I2(",m,"))"),
        "cartan":CartanMat("I",2,m),
        "size":GAPMul(2,m),
        "irredinfo":map(lambda x: {"charparam":x,
        "charname":ChevieData["I"]["CharName"](m,x,{})},ChevieData["I"]["CharInfo"](m)["charparams"]),
        "parameter":[u,v],
        "powermap":[],
        "irreducibles":GAPMul(ct,v**0)}
    tbl.update(cl)
    tbl["centralizers"]=map(lambda i: GAPDiv(tbl["size"],i),tbl["classes"])
    tbl=ChevieData["compat"]["MakeCharacterTable"](tbl)
    ChevieData["compat"]["AdjustHeckeCharTable"](tbl,param)
    return tbl

def coxi16(j):
    l=[]
    for i in range(1,len(r)+1):
        k=GAPDiv(len(r[i-1]),2)
        if r[i-1]==[] :
            l[i-1]=GAPMul(2,v**0)
        else:
            if r[i-1]==[1] :
                l[i-1]=u-1
            else:
                if r[i-1]==[2] :
                    l[i-1]=v-1
                else:
                    l[i-1]=GAPMul(squv**k,ER(m)**GAPMul(k,j)+ER(m)**GAPMul(-k,j))
    return l

ChevieData["I"]["HeckeCharTable"]=coxi15

def coxi17(m,i):
    return ChevieData["I"]["HeckeRepresentation"](m,[[1,-1],[1,-1]],[1,1],i)

ChevieData["I"]["Representation"]=coxi17

def coxi18(m,param,rootparam,i):
    if i==1 :
        return [[[param[1-1][1-1]]],[[param[2-1][1-1]]]]
    if m%2==0 :
        i=i-2
    if i==0 :
        return [[[param[1-1][1-1]]],[[param[2-1][2-1]]]]
    else:
        if i==1 :
            return [[[param[1-1][2-1]]],[[param[2-1][1-1]]]]
        else:
            if i==2 :
                return [[[param[1-1][2-1]]],[[param[2-1][2-1]]]]
            else:
                u=GAPDiv(-param[1-1][1-1],param[1-1][2-1])
                v=GAPDiv(-param[2-1][1-1],param[2-1][2-1])
                if m%2!=0 :
                    squv=u
                else:
                    if rootparam[1]==None and rootparam[2]==None :
                        squv=GAPMul(rootparam[1-1],rootparam[2-1])
                    else:
                        squv=GetRoot(GAPMul(u,v),2,"Representation(Hecke(I2(",m,")),[",i,"])")
                return [GAPMul(-[[-u**0,u**0],[GAPMul(0,u),u]],param[1-1][2-1]),GAPMul(-[[v,GAPMul(0,v)],[u+v,-v**0]],param[2-1][2-1])]

ChevieData["I"]["HeckeRepresentation"]=coxi18

def coxi19(m,sqrtu,j):
    return GAPMul([[0,GAPDiv(GAPDiv(1,sqrtu),ER(GAPMul(2,m))**j+ER(GAPMul(2,m))**-j)],[GAPMul(sqrtu,ER(GAPMul(2,m))**j+ER(GAPMul(2,m))**-j),0]],sqrtu**0)

ChevieData["I"]["Frobenius"]=coxi19

def coxi20(m,param):
    u=GAPDiv(-param[1-1][1-1],param[1-1][2-1])
    v=GAPDiv(-param[2-1][1-1],param[2-1][2-1])
    if IsInt(GAPDiv(m,2)) :
        return GAPMul(Sum(range(1,GAPDiv(m,2)+1),lambda i: GAPMul(u,v)**i-1),u+1)
    else:
        return GAPMul(Sum(range(1,m+1),lambda i: u**i-1),u+1)

ChevieData["I"]["PoincarePolynomial"]=coxi20

def coxi21(m,phi,para,rootpara):
    if m%2==1 :
        ci=ChevieData["I"]["CharInfo"](m)
        ci=ci["malleParams"][ci["charparams"].index(phi)+1-1]
        return GAPDiv(ChevieData["imp"]["SchurElement"](m,1,2,ci,[map(lambda i: ER(m)**i,range(0,m-1+1)),para[2-1]],[]),m)
    u=GAPDiv(-para[1-1][1-1],para[1-1][2-1])
    v=GAPDiv(-para[2-1][1-1],para[2-1][2-1])
    if phi[1-1]==1 :
        if phi[2-1]==GAPDiv(m,2) :
            e=GAPDiv(GAPMul(Sum(range(0,GAPDiv(m,2)-1+1),lambda i: GAPDiv(u,v)**i),u+1),v)
            if phi[3-1]=="'" :
                return e
            else:
                return GAPMul(GAPDiv(v,u)**GAPDiv(m,2),e)
        else:
            e=GAPMul(Sum(range(0,GAPDiv(m,2)-1+1),lambda i: GAPMul(u,v)**i),u+1)
            if phi[2-1]==0 :
                return e
            else:
                return GAPMul(GAPMul(u,v)**GAPDiv(-m,2),e)
    else:
        e=ER(m)**phi[2-1]+ER(m)**-phi[2-1]
        if ForAll([1,2],lambda i: rootpara[i]==None) :
            ruv=prod(rootpara)
        else:
            ruv=GetRoot(GAPMul(u,v),2,"SchurElement(Hecke(I2(",m,"),",phi,"))")
        return GAPDiv(GAPMul(-m,GAPMul(u,v)+1-GAPMul(ruv,e)),GAPMul(u,v))

ChevieData["I"]["SchurElement"]=coxi21

def coxi22(m,phi,q):
    if phi[1-1]==1 :
        return q**phi[2-1]
    else:
        return q**phi[2-1]+q**m-phi[2-1]

ChevieData["I"]["FakeDegree"]=coxi22

def coxi23(m):
    res=ChevieData["I"]["HeckeCharTable"](m,[[1,-1],[1,-1]],[1,1])
    res["identifier"]=SPrint("W(I2(",m,"))")
    return res

ChevieData["I"]["CharTable"]=coxi23

def coxi24(n,p):
    T=ChevieData["I"]["CharTable"](n)
    T["name"]=T["identifier"]
    m=DecompositionMatrix(T%p)
    return map(lambda c: [c[1-1],[[m[k-1] for k in c[1-1]][k-1] for k in c[2-1]]],BlocksMat(m))

ChevieData["I"]["DecompositionMatrix"]=coxi24

def coxi25(arg):
    if arg[1-1]%2==0 and arg[3-1][1-1]!=arg[3-1][2-1] :
        Error(" !  implemented")
    ci=ChevieData["I"]["CharInfo"](arg[1-1])
    ci=ci["malleParams"][ci["charparams"].index(arg[2-1])+1-1]
    return ChevieData["imp"]["FactorizedSchurElement"](arg[1-1],arg[1-1],2,ci,arg[3-1],1)

ChevieData["I"]["FactorizedSchurElement"]=coxi25

def coxi26(arg):
    e=arg[1-1]
    if len(arg)==2 :
        type_=arg[2-1]
    else:
        if e%2==0 :
            type_=1
        else:
            type_=-ER(e)**GAPDiv(e+1,2)-ER(e)**GAPDiv(e+3,2)
    m=GAPMul(DiagonalMat(1+ER(e)**-1,-type_),ChevieData["imp"]["GeneratingRoots"](e,e,2))
    return map(lambda f: coxi27,ChevieData["imp"]["Invariants"](e,e,2))

def coxi27(arg):
    return f(*GAPMul(arg,m))

ChevieData["I"]["Invariants"]=coxi26

def coxi28(S):
    if S[1-1]!=[0,1] or not [] in S :
        return false
    if len(S)%2==1 :
        S=Reversed(S)
        return [S.index([])+1,S.index([0,1])+1-S.index([])+1]
    else:
        return S.index([])+1+[-[S[k-1] for k in range(2,len(S)+1)].index([0,1])+1,0]-1

ChevieData["I"]["SymbolToParameter"]=coxi28

def coxi29(e,p):
    if p==[0] :
        S=map(lambda x: [0],range(1,e+1))
        S[e-1]=[2]
    else:
        if p==[1] :
            S=map(lambda x: [0,1],range(1,e+1))
            S[e-1]=[1,2]
        else:
            if len(p)==3 :
                S=map(lambda x: [0],range(1,GAPDiv(e,2)-1+1))
                S+=[[1],2,GAPDiv(p[3-1]+1,2)]
            else:
                if e%2==0 :
                    S=map(lambda x: [0],range(1,e+1))
                    if p[1-1]==0 :
                        for i,j in zip([e,e-p[2-1]],[[1],[1]]):
                            S[i-1]=j
                    else:
                        for i,j in zip(1+[0,p[2-1]-p[1-1]%e],[[0,1],[0,1]]):
                            S[i-1]=j
                        for i,j in zip(1+[-p[1-1]%e,p[2-1]],[[],[]]):
                            S[i-1]=j
                else:
                    S=map(lambda i: [0],range(1,e+1))
                    if p[1-1]!=0 :
                        for i,j in zip(1+[0,-Sum(p)%e],[[0,1],[0,1]]):
                            S[i-1]=j
                        for i,j in zip(1+map(lambda x: x%e,-p),[[],[]]):
                            S[i-1]=j
                    else:
                        for i,j in zip(e+[-p[2-1]-p[1-1]%e,0],[[1],[1]]):
                            S[i-1]=j
    return S

ChevieData["I"]["ParameterToSymbol"]=coxi29

def coxi30(e):
    f=QuoInt(e,2)
    uc={}
    uc["harishChandra"]=[{"relativeType":{"series":"I",
        "indices":[1,2],
        "rank":2,
        "bond":e},
        "parameterExponents":[1,1],
        "levi":[],
        "eigenvalue":1,
        "cuspidalName":""}]
    if e%2!=0 :
        uc["harishChandra"][1-1]["charNumbers"]=range(1,f+2+1)
    else:
        uc["harishChandra"][1-1]["charNumbers"]=Concatenation([1,3,4,2],4+range(1,f-1+1))
    cusp=Concatenation(map(lambda k: map(lambda l: [k,l],range(k+1,e-k-1+1)),range(1,f-1+1)))
    f=f+1-e%2
    uc["harishChandra"]+=map(lambda x: {"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "parameterExponents":[],
        "levi":[1,2],
        "eigenvalue":ER(e)**-prod(cusp[x-1]),
        "cuspidalName":SPrint("I_2(",e,")",FormatGAP(cusp[x-1])),
        "charNumbers":[x+f]},range(1,len(cusp)+1))
    uc["families"]=[Family(ChevieData["families"]["Dihedral"](e),range(1,len(cusp)+f+1)+2),Family("C1",[1]),Family("C1",[2])]
    uc["parameters"]=Concatenation([[0],[1]],uc["families"][1-1]["parameters"])
    uc["charSymbols"]=map(lambda p: ChevieData["I"]["ParameterToSymbol"](e,p),uc["parameters"])
    uc["a"]=Concatenation([0,e],map(lambda x: 1,uc["families"][1-1]["parameters"]))
    uc["A"]=Concatenation([0,e],map(lambda x: e-1,uc["families"][1-1]["parameters"]))
    if e==5 :
        uc["curtis"]=[2,1,3,4,6,5]
    return uc

ChevieData["I"]["UnipotentCharacters"]=coxi30
