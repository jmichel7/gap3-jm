
ChevieData["H3"]["ReflectionDegrees"]=[2,6,10]

ChevieData["H3"]["Size"]=120

def coxh31():
    a=GAPDiv(1+ER(5),2)
    return GAPDiv([[a,-1,a-1],[-a,1,a-1],[1,a-1,-a]],2)

ChevieData["H3"]["GeneratingRoots"]=coxh31

ChevieData["H3"]["NrConjugacyClasses"]=10

ChevieData["H3"]["cyclestructure"]=[]

ChevieData["H3"]["generators"]=[]

ChevieData["H3"]["CartanMat"]=[[2,GAPDiv(-1+ER(5),2),0],[GAPDiv(-1+ER(5),2),2,-1],[0,-1,2]]

ChevieData["H3"]["PowerMaps"]=[None,[1,1,7,1,5,3,3,5,7,1],[1,2,7,4,1,9,3,10,6,10],None,[1,2,1,4,5,10,1,8,10,10],None,[1,2,7,4,5,9,3,8,6,10]]

ChevieData["H3"]["WordsClassRepresentatives"]=[[],[1],[1,2],[1,3],[2,3],[1,2,3],[1,2,1,2],[1,2,1,2,3],[1,2,1,2,3,2,1,2,3],[1,2,1,2,1,3,2,1,2,1,3,2,1,2,3]]

def coxh32(s):
    t=[[[]],[[1]],[[1,2],[1,3],[2,3]],[range(1,3+1)]]
    return t[s+1-1]

ChevieData["H3"]["ParabolicRepresentatives"]=coxh32

def coxh33():
    res={"classtext":ChevieData["H3"]["WordsClassRepresentatives"],
        "orders":[1,2,5,2,3,10,5,6,10,2],
        "classes":[1,15,12,15,20,12,12,20,12,1]}
    res["classnames"]=map(IntListToString,res["classtext"])
    res["classnames"][1-1]="."
    res["classparams"]=res["classnames"]
    return res

ChevieData["H3"]["ClassInfo"]=coxh33

def coxh34():
    res={"charparams":[[1,15],[1,0],[5,5],[5,2],[3,6],[3,8],[3,1],[3,3],[4,3],[4,4]],
        "gp":["1_r'","1_r","5_r'","5_r","3_s","overline{3}_s","3_s'","overline{3}_s'","4_r'","4_r"],
        "opdam":Permutation("(9,10)"),
        "extRefl":[2,7,5,1]}
    res["b"]=map(lambda x: x[2-1],res["charparams"])
    return res

ChevieData["H3"]["CharInfo"]=coxh34

ChevieData["H3"]["vpolheckeirreducibles"]=[[[[1],0],[[-1],0],[[1],0],[[1],0],[[1],0],[[-1],0],[[1],0],[[-1],0],[[-1],0],[[-1],0]],[[[1],0],[[1],1],[[1],2],[[1],2],[[1],2],[[1],3],[[1],4],[[1],5],[[1],9],[[1],15]],[[[5],0],[[-3,2],0],[[1,-1],0],[[2,-2,1],0],[[1,-2],0],[[],0],[[1,0,-1],0],[[1],2],[[],0],[[-5],6]],[[[5],0],[[-2,3],0],[[-1,1],1],[[1,-2,2],0],[[-2,1],1],[[],0],[[-1,0,1],2],[[-1],3],[[],0],[[5],9]],[[[3],0],[[-2,1],0],[[1,GAPDiv(-1+ER(5),2)],0],[[1,-2],0],[[1,-1],0],[[GAPDiv(1-ER(5),2)],1],[[1,0,GAPDiv(-1-ER(5),2)],0],[[],0],[[GAPDiv(1+ER(5),2)],3],[[3],5]],[[[3],0],[[-2,1],0],[[1,GAPDiv(-1-ER(5),2)],0],[[1,-2],0],[[1,-1],0],[[GAPDiv(1+ER(5),2)],1],[[1,0,GAPDiv(-1+ER(5),2)],0],[[],0],[[GAPDiv(1-ER(5),2)],3],[[3],5]],[[[3],0],[[-1,2],0],[[GAPDiv(-1+ER(5),2),1],1],[[-2,1],1],[[-1,1],1],[[GAPDiv(-1+ER(5),2)],2],[[GAPDiv(-1-ER(5),2),0,1],2],[[],0],[[GAPDiv(-1-ER(5),2)],6],[[-3],10]],[[[3],0],[[-1,2],0],[[GAPDiv(-1-ER(5),2),1],1],[[-2,1],1],[[-1,1],1],[[GAPDiv(-1-ER(5),2)],2],[[GAPDiv(-1+ER(5),2),0,1],2],[[],0],[[GAPDiv(-1+ER(5),2)],6],[[-3],10]],[[[4],0],[[-2,2],0],[[-1],1],[[1,-2,1],0],[[1,-1,1],0],[[1],GAPDiv(3,2)],[[-1],2],[[-1],GAPDiv(5,2)],[[1],GAPDiv(9,2)],[[-4],GAPDiv(15,2)]],[[[4],0],[[-2,2],0],[[-1],1],[[1,-2,1],0],[[1,-1,1],0],[[-1],GAPDiv(3,2)],[[-1],2],[[1],GAPDiv(5,2)],[[-1],GAPDiv(9,2)],[[4],GAPDiv(15,2)]]]

ChevieData["H3"]["CycPolSchurElements"]=[[1,-15,2,2,2,3,5,6,10],[1,0,2,2,2,3,5,6,10],[1,-5,2,2,2,3,6],[1,-2,2,2,2,3,6],[GAPDiv(5+ER(5),2),-6,2,2,2,GAPDiv(2,5),GAPDiv(3,5),GAPDiv(1,10),GAPDiv(9,10)],[GAPDiv(5-ER(5),2),-6,2,2,2,GAPDiv(1,5),GAPDiv(4,5),GAPDiv(3,10),GAPDiv(7,10)],[GAPDiv(5+ER(5),2),-1,2,2,2,GAPDiv(2,5),GAPDiv(3,5),GAPDiv(1,10),GAPDiv(9,10)],[GAPDiv(5-ER(5),2),-1,2,2,2,GAPDiv(1,5),GAPDiv(4,5),GAPDiv(3,10),GAPDiv(7,10)],[2,-3,3,5],[2,-3,3,5]]

ChevieData["H3"]["sparseFakeDegrees"]=[[1,15],[1,0],[1,5,1,7,1,9,1,11,1,13],[1,2,1,4,1,6,1,8,1,10],[1,6,1,10,1,14],[1,8,1,10,1,12],[1,1,1,5,1,9],[1,3,1,5,1,7],[1,3,1,7,1,9,1,11],[1,4,1,6,1,8,1,12]]

def coxh35(param,sqrtparam):
    a=GAPDiv(1+ER(5),2)
    q=GAPDiv(-param[1-1][1-1],param[1-1][2-1])
    if not sqrtparam[1]==None :
        v=GetRoot(q,2,"CharTable(Hecke(H3))")
    else:
        v=sqrtparam[1-1]
    ci=ChevieData["H3"]["ClassInfo"]()
    tbl={"identifier":"H(H3)",
        "text":"the representing matrices are those of Lusztig(1981)",
        "parameter":[q,q,q],
        "cartan":ChevieData["H3"]["CartanMat"],
        "size":120,
        "order":120,
        "powermap":ChevieData["H3"]["PowerMaps"],
        "irreducibles":map(lambda i: map(coxh36,i),ChevieData["H3"]["vpolheckeirreducibles"]),
        "irredinfo":map(lambda x: {"charparam":x,
        "charname":ChevieData["H3"]["CharName"](x,{})},ChevieData["H3"]["CharInfo"]()["charparams"])}
    tbl.update(ci)
    tbl["centralizers"]=map(lambda x: GAPDiv(tbl["size"],x),tbl["classes"])
    tbl=ChevieData["compat"]["MakeCharacterTable"](tbl)
    ChevieData["compat"]["AdjustHeckeCharTable"](tbl,param)
    return tbl

def coxh36(j):
    res=ValuePol(j[1-1],q)
    if IsInt(j[2-1]) :
        res=GAPMul(res,q**j[2-1])
    else:
        res=GAPMul(res,v**GAPMul(2,j[2-1]))
    return res

ChevieData["H3"]["HeckeCharTable"]=coxh35

def coxh37(i):
    return ChevieData["H3"]["HeckeRepresentation"]([[1,-1],[1,-1],[1,-1]],[1,1,1],i)

ChevieData["H3"]["Representation"]=coxh37

ChevieData["H3"]["WGraphs"]=[[[[1,2,3]],[]],1,[[[2],[1,2],[1,3],[1,3],[2,3]],[[-1,[[1,3],[2,4],[3,5],[4,5]]]]],3,[[[1,2],[1,3],[2,3]],[[-1,[[1,2]]],[GAPDiv(-1-ER(5),2),[[2,3]]]]],[[[1,2],[1,3],[2,3]],[[-1,[[1,2]]],[GAPDiv(-1+ER(5),2),[[2,3]]]]],5,6,[[[1],[2],[1,3],[2,3]],[[1,[[1,2,3],[2,3,4],[3,4]]]]],9]

def coxh38(i):
    gr=ChevieData["H3"]["WGraphs"]
    if IsInt(gr[i-1]) :
        return DualWGraph(3,gr[gr[i-1]-1])
    else:
        return gr[i-1]

ChevieData["H3"]["WGraph"]=coxh38

def coxh39(param,sqrtparam,i):
    if not sqrtparam[1]==None :
        v=GetRoot(GAPDiv(-param[1-1][1-1],param[1-1][2-1]),2,"Representation(Hecke(H3),[",i,"])")
    else:
        v=sqrtparam[1-1]
    return GAPMul(-param[1-1][2-1],WGraphToRepresentation(3,ChevieData["H3"]["WGraph"](i),v))

ChevieData["H3"]["HeckeRepresentation"]=coxh39

def coxh310():
    res={"harishChandra":[{"relativeType":{"series":"H",
        "indices":range(1,3+1),
        "rank":3},
        "levi":[],
        "eigenvalue":1,
        "parameterExponents":[1,1,1],
        "cuspidalName":"",
        "charNumbers":range(1,10+1)},{"relativeType":{"series":"A",
        "indices":[3],
        "rank":1},
        "levi":range(1,2+1),
        "eigenvalue":ER(5)**2,
        "parameterExponents":[5],
        "cuspidalName":"I_2(5)[1,3]",
        "charNumbers":[11,13]},{"relativeType":{"series":"A",
        "indices":[3],
        "rank":1},
        "levi":range(1,2+1),
        "eigenvalue":ER(5)**3,
        "parameterExponents":[5],
        "cuspidalName":"I_2(5)[1,2]",
        "charNumbers":[12,14]},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "eigenvalue":ER(4),
        "qEigen":GAPDiv(1,2),
        "parameterExponents":[],
        "cuspidalName":"H_3[i]",
        "charNumbers":[15]},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "eigenvalue":-ER(4),
        "qEigen":GAPDiv(1,2),
        "parameterExponents":[],
        "cuspidalName":"H_3[-i]",
        "charNumbers":[16]}],
        "families":[Family("C1",[2]),Family(ChevieData["families"]["Dihedral"](5),[7,8,14,13]),Family("C1",[4]),Family("C'\"2",[9,10,15,16]),Family("C1",[3]),Family(ChevieData["families"]["Dihedral"](5),[5,6,12,11]),Family("C1",[1])],
        "a":[15,0,5,2,6,6,1,1,3,3,6,6,1,1,3,3],
        "A":[15,0,13,10,14,14,9,9,12,12,14,14,9,9,12,12]}
    return res

ChevieData["H3"]["UnipotentCharacters"]=coxh310

def coxh311():
    C=ChevieData["H3"]["CartanMat"]
    r=GAPMul(RootsCartan(C),C)
    return map(lambda d: coxh312,ChevieData["H3"]["ReflectionDegrees"])

def coxh312(arg):
    return Sum(r,lambda a: GAPMul(arg,a)**d)

ChevieData["H3"]["Invariants"]=coxh311

def coxh313():
    return coxh314

def coxh314(a,b,c):
    return GAPMul(131835937500,a)-GAPMul(100195312500,a**2)+GAPMul(395507812500,c**3)-GAPMul(28369140625,a**3)+GAPMul(1371093750,a**4)-GAPMul(74250000,a**7)-GAPMul(22233750,a**9)+GAPMul(438750,a**10)-GAPMul(829,a**15)

ChevieData["H3"]["Discriminant"]=coxh313

ChevieData["H3"]["KLeftCellRepresentatives"]=[{"character":[2],
    "duflo":[1,2,3],
    "reps":""},{"character":[1],
    "duflo":[16,17,18],
    "reps":""},{"character":[3],
    "duflo":[1,24,3],
    "reps":""},{"character":[4],
    "duflo":[2,1,28],
    "reps":""},{"character":[6,5],
    "duflo":[1,20,18],
    "reps":[[7,19,24]]},{"character":[8,7],
    "duflo":[1,6,18],
    "reps":[[9,2,27]]},{"character":[10,9],
    "duflo":[8,18,17],
    "reps":[[11,17,25],[11,27,10],[14,30,4]]},{"character":[10,9],
    "duflo":[13,30,8],
    "reps":[[10,29,5],[12,21,22],[13,22,23]]}]
