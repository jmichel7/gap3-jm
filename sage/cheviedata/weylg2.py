
def eylg21(arg):
    if len(arg)>0 :
        type_=arg[1-1]
    else:
        type_=1
    a=[[2,-1],[-3,2]]
    a[1-1][2-1]=-type_
    a[2-1][1-1]=GAPDiv(3,a[1-1][2-1])
    return a

ChevieData["G2"]["CartanMat"]=eylg21

def eylg22(indices,title,type_):
    print title," ",indices[1-1],
    if type_==1 :
        print " >>> ",
    else:
        if type_==ER(3) :
            print " ==6== ",
        else:
            print " ?6? ",
    print indices[2-1]," \n",

ChevieData["G2"]["PrintDiagram"]=eylg22

def eylg23(arg):
    if len(arg)==1 :
        return "G2(?)"
    type_=arg[1-1]
    opt=arg[2-1]
    if type_==1 :
        if "TeX" in opt :
            return "G_2"
        else:
            if "arg" in opt :
                return "\"G\",2"
            else:
                return "G2"
    else:
        if type_==ER(3) :
            if "TeX" in opt :
                return "G_{\\hbox{sym}2}"
            else:
                if "arg" in opt :
                    return "\"Gsym\",2"
                else:
                    return "Gsym2"
        else:
            if "TeX" in opt :
                return SPrint("G_2(",Format(GAPDiv(type_**2,3),opt),")")
            else:
                if "arg" in opt :
                    return SPrint("\"G\",",2,",",Format(GAPDiv(type_**2,3),opt))
                else:
                    return SPrint("G2(",Format(GAPDiv(type_**2,3),opt),")")

ChevieData["G2"]["ReflectionName"]=eylg23

ChevieData["G2"]["ParabolicRepresentatives"]=lambda s: ChevieData["imp"]["ParabolicRepresentatives"](6,6,2,s)

ChevieData["G2"]["GeneratingRoots"]=[[1,-1,0],[-2,1,1]]

ChevieData["G2"]["HyperplaneRepresentatives"]=[1,2]

ChevieData["G2"]["Size"]=12

ChevieData["G2"]["ReflectionDegrees"]=[2,6]

ChevieData["G2"]["NrConjugacyClasses"]=6

def eylg24():
    res={"charparams":[[1,0],[1,6],[1,3,1],[1,3,2],[2,1],[2,2]],
        "extRefl":[1,5,2],
        "a":[0,6,1,1,1,1],
        "A":[0,6,5,5,5,5]}
    res["b"]=map(lambda x: x[2-1],res["charparams"])
    res["spaltenstein"]=["1","\\varepsilon","\\varepsilon_l","\\varepsilon_c","\\theta'","\\theta''"]
    return res

ChevieData["G2"]["CharInfo"]=eylg24

ChevieData["G2"]["WordsClassRepresentatives"]=[[],[2],[1],[1,2],[1,2,1,2],[1,2,1,2,1,2]]

ChevieData["G2"]["ClassNames"]=["A_0","\\tilde A_1","A_1","G_2","A_2","A_1+\\tilde A_1"]

ChevieData["G2"]["ClassInfo"]={"classtext":ChevieData["G2"]["WordsClassRepresentatives"],
    "classnames":ChevieData["G2"]["ClassNames"],
    "classparams":ChevieData["G2"]["ClassNames"],
    "orders":[1,2,2,6,3,2],
    "classes":[1,3,3,2,2,1]}

ChevieData["G2"]["PowerMaps"]=[None,[1,1,1,5,5,1],[1,2,3,6,1,6]]

ChevieData["G2"]["sparseFakeDegrees"]=[[1,0],[1,6],[1,3],[1,3],[1,1,1,5],[1,2,1,4]]

ChevieData["G2"]["ClassParameter"]=lambda w: ChevieData["G2"]["ClassNames"][PositionProperty([[[]],[[2],[1,2,1],[2,1,2,1,2]],[[1],[2,1,2],[1,2,1,2,1]],[[2,1],[1,2]],[[2,1,2,1],[1,2,1,2]],[[1,2,1,2,1,2]]],lambda x: w in x)-1]

def eylg25(para,sqrtpara):
    u=prod(para[1-1])
    v=prod(para[2-1])
    if u==v :
        return u
    else:
        if u==v**3 :
            return -v**2
        else:
            if v==u**3 :
                return -u**2
            else:
                if sqrtpara[1]==None and sqrtpara[2]==None :
                    return GAPMul(sqrtpara[1-1],sqrtpara[2-1])
                else:
                    return GetRoot(GAPMul(u,v),2,"Hecke(G2)")

ChevieData["G2"]["squv"]=eylg25

def eylg26(para,sqrtpara):
    x=para[1-1][1-1]
    y=para[1-1][2-1]
    z=para[2-1][1-1]
    t=para[2-1][2-1]
    one=GAPMul(x,y)**0
    def f1(u,v):
        return GAPMul([1,v,u,GAPMul(v,u),GAPMul(v**2,u**2),GAPMul(v**3,u**3)],one)
    
    
    def f2(x,y,z,t,eps):
        squv=GAPMul(eps,ChevieData["G2"]["squv"](para,sqrtpara))
        return GAPMul([2,z+t,x+y,-squv,GAPMul(-x,y),GAPMul(2,squv**3)],one)
    
    
    tbl={"identifier":"H(G2)",
        "parameter":[[x,y],[z,t]],
        "size":12,
        "powermap":ChevieData["G2"]["PowerMaps"],
        "irreducibles":[f1(x,z),f1(y,t),f1(y,z),f1(x,t),f2(x,y,z,t,1),f2(x,y,z,t,-1)],
        "irredinfo":map(lambda x: {"charparam":x,
        "charname":ChevieData["G2"]["CharName"](x,{})},ChevieData["G2"]["CharInfo"]()["charparams"])}
    tbl.update(ChevieData["G2"]["ClassInfo"])
    tbl["centralizers"]=map(lambda x: GAPDiv(tbl["size"],x),tbl["classes"])
    tbl=ChevieData["compat"]["MakeCharacterTable"](tbl)
    return tbl

ChevieData["G2"]["HeckeCharTable"]=eylg26

def eylg27(para,sqrtpara,i):
    one=GAPMul(prod(para[1-1])**0,prod(para[2-1])**0)
    x=para[1-1][1-1]
    y=para[1-1][2-1]
    z=para[2-1][1-1]
    t=para[2-1][2-1]
    if i==1 :
        return GAPMul([[[x]],[[z]]],one)
    else:
        if i==2 :
            return GAPMul([[[y]],[[t]]],one)
        else:
            if i==3 :
                return GAPMul([[[y]],[[z]]],one)
            else:
                if i==4 :
                    return GAPMul([[[x]],[[t]]],one)
                else:
                    squv=ChevieData["G2"]["squv"](para,sqrtpara)
                    if i==6 :
                        squv=-squv
                    return GAPMul([[[y,-1],[0,x]],[[z,0],[squv+GAPMul(y,z),t]]],one)

ChevieData["G2"]["HeckeRepresentation"]=eylg27

def eylg28(i):
    return ChevieData["G2"]["HeckeRepresentation"]([[1,-1],[1,-1]],[1,1],i)

ChevieData["G2"]["Representation"]=eylg28

def eylg29(param):
    u=GAPDiv(-param[1-1][1-1],param[1-1][2-1])
    v=GAPDiv(-param[2-1][1-1],param[2-1][2-1])
    return GAPMul(1+u,v+1)

ChevieData["G2"]["PoincarePolynomial"]=eylg29

ChevieData["G2"]["SchurModels"]={"f1":{"vcyc":[[[1,-1,0,0],1],[[0,0,1,-1],1],[[1,-1,1,-1],3]]},
    "f2":{"coeff":-2,
    "root":GAPDiv([1,-1,1,-1],2),
    "factor":[-1,1,0,0],
    "vcyc":[[[0,0,0,0,1],3],[[0,0,-1,1,1],3]]}}

ChevieData["G2"]["SchurData"]=[{"name":"f1",
    "order":[1,2,3,4]},{"name":"f1",
    "order":[2,1,4,3]},{"name":"f1",
    "order":[2,1,3,4]},{"name":"f1",
    "order":[1,2,4,3]},{"name":"f2",
    "order":[1,2,3,4],
    "rootPower":-1},{"name":"f2",
    "order":[1,2,3,4],
    "rootPower":1}]

def eylg210(phi,para,sqrtpara):
    u=GAPDiv(-para[1-1][1-1],para[1-1][2-1])
    v=GAPDiv(-para[2-1][1-1],para[2-1][2-1])
    p=ChevieData["G2"]["CharInfo"]()["charparams"].index(phi)+1
    if p==1 :
        return GAPMul(1+u,v+1)
    else:
        if p==2 :
            return GAPDiv(GAPDiv(GAPMul(1+u,v+1),u**3),v**3)
        else:
            if p==3 :
                return GAPDiv(GAPMul(u**2+v**2,1+u),u**3)
            else:
                if p==4 :
                    return GAPDiv(GAPMul(u**2+v**2,1+u),v**3)
    squv=GAPDiv(GAPDiv(ChevieData["G2"]["squv"](para,sqrtpara),para[1-1][2-1]),para[2-1][2-1])
    if p==6 :
        squv=-squv
    return GAPMul(2,GAPMul(u,v)**-1)

ChevieData["G2"]["SchurElement"]=eylg210

def eylg211():
    return {"harishChandra":[{"relativeType":{"series":"G",
        "indices":range(1,2+1),
        "rank":2},
        "levi":[],
        "parameterExponents":[1,1],
        "charNumbers":range(1,6+1),
        "eigenvalue":1,
        "cuspidalName":""},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,2+1),
        "parameterExponents":[],
        "charNumbers":[10],
        "eigenvalue":ER(3)**2,
        "cuspidalName":"G_2[\\zeta_3^2]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,2+1),
        "parameterExponents":[],
        "charNumbers":[7],
        "eigenvalue":-1,
        "cuspidalName":"G_2[-1]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,2+1),
        "parameterExponents":[],
        "charNumbers":[9],
        "eigenvalue":ER(3),
        "cuspidalName":"G_2[\\zeta_3]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,2+1),
        "parameterExponents":[],
        "charNumbers":[8],
        "eigenvalue":1,
        "cuspidalName":"G_2[1]"}],
        "families":[Family("S3",[5,6,4,3,8,7,9,10]),Family("C1",[1]),Family("C1",[2])],
        "a":[0,6,1,1,1,1,1,1,1,1],
        "A":[0,6,5,5,5,5,5,5,5,5],
        "charSymbols":[[[0],[0],[0],[0],[0],[2]],[[0,1],[0,1],[0,1],[0,1],[0,1],[1,2]],[[0],[0],[1],2,0],[[0],[0],[1],2,1],[[0],[0],[0],[0],[1],[1]],[[0],[0],[0],[1],[0],[1]],[[0,1],[0],[0,1],[],[0],[]],[[0,1],[0,1],[0],[],[],[0]],[[0,1],[0],[0],[0,1],[],[]],[[0,1],[0,1],[],[0],[0],[]]]}

ChevieData["G2"]["UnipotentCharacters"]=eylg211

def eylg212(x,y):
    return GAPMul(-3,x)+GAPMul(3,x**2)

def eylg213(x,y):
    return GAPMul(x**2,y**4)-GAPMul(6,x**3)+GAPMul(13,x**4)-GAPMul(12,x**5)+GAPMul(4,x**6)

ChevieData["G2"]["Invariants"]=[eylg212,eylg213]

def eylg214():
    return eylg215

def eylg215(x,y):
    return GAPMul(4,x**3)-GAPMul(27,y**2)

ChevieData["G2"]["Discriminant"]=eylg214

def eylg216(c,p):
    if p==0 :
        p=1
    Z=lambda n: ReflectionGroup(n,1,1)
    uc={"classes":[{"name":"1",
        "succ":["A1"],
        "dynkin":[0,0],
        "balacarter":[],
        "red":CoxeterGroup("G",2)},{"name":"A_1",
        "succ":["~A1"],
        "dynkin":[1,0],
        "balacarter":[1],
        "red":Z(2)},{"name":"\\tilde A_1",
        "succ":["G2(a1)"],
        "dynkin":[0,1],
        "balacarter":[2],
        "red":Z(2-GAPDiv(Gcd(p,3)-1,2))},{"name":"G_2(a_1)",
        "succ":["G2"],
        "dynkin":[2,0],
        "balacarter":[1,-2],
        "Au":CoxeterGroup("A",2-GAPDiv(Gcd(p,3)-1,2))},{"name":"G_2",
        "succ":[],
        "dynkin":[2,2],
        "Au":Z(Gcd(p,6)),
        "balacarter":[1,2]}],
        "springerSeries":[{"relgroup":CoxeterGroup("G",2),
        "levi":"",
        "Z":[],
        "locsys":[[5,1],[1,1],[4,2],[2,1],[4,3],[3,1]]},{"relgroup":CoxeterGroup(),
        "levi":[1,2],
        "Z":[],
        "locsys":[[4,1]],
        "parameter":[8]}]}
    if p==2 :
        uc["springerSeries"][1-1]["locsys"][1-1]=[5,2]
        uc["springerSeries"].append({"relgroup":CoxeterGroup(),
            "levi":[1,2],
            "Z":[],
            "locsys":[[5,1]]})
    else:
        if p==3 :
            uc["classes"].append({"name":"(\\tilde A_1)_3",
                "succ":["~A1"],
                "dimBu":3,
                "red":Z(2),
                "Au":CoxeterGroup()})
            uc["classes"][1-1]["succ"].append("(~A1)3")
            uc["classes"][3-1]["dimBu"]=2
            del uc["classes"][3-1]["dynkin"]
            for i,j in zip([3,5],[[6,1],[4,2]]):
                uc["springerSeries"][1-1]["locsys"][i-1]=j
            for c in [2,3]:
                uc["springerSeries"].append({"relgroup":CoxeterGroup(),
                    "levi":[1,2],
                    "Z":[],
                    "locsys":[[5,c]]})
    uc["orderClasses"]=map(lambda c: map(lambda n: PositionProperty(uc["classes"],lambda c: UnipotentClassOps["Name"](c)==n),c["succ"]),uc["classes"])
    for c in uc["classes"]:
        del c["succ"]
        if not "red" in c :
            c["red"]=Z(1)
        if not "Au" in c :
            c["Au"]=Z(1)
        c["AuAction"]=ExtendedReflectionGroup(c["red"],map(lambda x: IdentityMat(c["red"]["rank"]),c["Au"]["generators"]))
    return uc

ChevieData["G2"]["UnipotentClasses"]=eylg216

ChevieData["G2"]["KLeftCellRepresentatives"]=[{"character":[1],
    "duflo":[1,2],
    "reps":""},{"character":[2],
    "duflo":[7,8],
    "reps":""},{"character":[3,5,6],
    "duflo":[5,8],
    "reps":[[6,10],[12,3]]},{"character":[4,5,6],
    "duflo":[7,3],
    "reps":[[5,10],[12,4]]}]
