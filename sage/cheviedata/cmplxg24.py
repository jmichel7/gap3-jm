
def lxg241(indices,title):
    print title," ",indices[1-1],"\n",
    s=just("",len(title)-1)
    print s," / \\\n",s,indices[2-1],"=====",indices[3-1],"  ",IntListToString([indices[k-1] for k in [2,3,1,2,3,1,2,3,1]]),"==",IntListToString([indices[k-1] for k in [3,2,3,1,2,3,1,2,3]]),"\n",

ChevieData["G24"]["PrintDiagram"]=lxg241

ChevieData["G24"]["GeneratingRoots"]=[[1,ER(-7),0],[1,-ER(-7),0],[GAPDiv(-1-ER(-7),2),GAPDiv(-7-GAPMul(3,ER(-7)),6),GAPDiv(-4,3)]]

ChevieData["G24"]["GeneratingCoRoots"]=GAPDiv([[1,GAPDiv(GAPMul(-3,ER(-7)),7),0],[1,GAPDiv(GAPMul(3,ER(-7)),7),0],[GAPDiv(-1+ER(-7),2),GAPDiv(-7+GAPMul(3,ER(-7)),14),GAPDiv(-1,2)]],2)

def lxg242():
    return GAPMul(ChevieData["G24"]["GeneratingCoRoots"],Matrix(ChevieData["G24"]["GeneratingRoots"]).transpose())

ChevieData["G24"]["CartanMat"]=lxg242

ChevieData["G24"]["EigenvaluesGeneratingReflections"]=[GAPDiv(1,2),GAPDiv(1,2),GAPDiv(1,2)]

ChevieData["G24"]["BraidRelations"]=[[[1,2,1],[2,1,2]],[[1,3,1],[3,1,3]],[[3,2,3,2],[2,3,2,3]],[[2,3,1,2,3,1,2,3,1],[3,2,3,1,2,3,1,2,3]]]

ChevieData["G24"]["AltPres"]=[{"gens":[[1],[2,3,-2],[2]],
    "rels":[[[1,2,1,2],[2,1,2,1]],[[2,3,2,3],[3,2,3,2]],[[1,3,1],[3,1,3]],[[2,1,2,3,1,2,3],[1,2,3,1,2,3,1]]]},{"gens":[[2],[3],[-3,-2,1,2,3]],
    "rels":[[[1,3,1,3],[3,1,3,1]],[[3,2,3,2],[2,3,2,3]],[[1,2,1,2],[2,1,2,1]],[[2,3,1,2,3,1,2],[1,2,3,1,2,3,1]],[[2,3,1,2,3,1,2],[3,1,2,3,1,2,3]]]}]

ChevieData["G24"]["ReflectionDegrees"]=[4,6,14]

ChevieData["G24"]["Size"]=prod(ChevieData["G24"]["ReflectionDegrees"])

ChevieData["G24"]["NrConjugacyClasses"]=12

def lxg243(s):
    t=[[[]],[[1]],[[1,2],[2,3]],[range(1,3+1)]]
    return t[s+1-1]

ChevieData["G24"]["ParabolicRepresentatives"]=lxg243

ChevieData["G24"]["WordsClassRepresentatives"]=[[],[1],[2,3],[1,3],[1,2,3,1,2,3,1,2,3],[1,2,3],[2,3,2,3],[1,2,3,1,2,3,1],[1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3],[1,2,3,1,2,3],[1,2,3,1,2,3,1,2,3,1,2],[1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3]]

ChevieData["G24"]["ClassNames"]=[".","1","23","13","ccc","c","2323","cc1","cccccc","cc","ccc12","z"]

ChevieData["G24"]["PowerMaps"]=[None,[1,1,7,4,9,10,1,4,9,10,7,1],[1,2,3,1,6,5,7,12,10,9,11,12],None,[1,2,3,4,6,5,7,8,10,9,11,12],None,[1,2,3,4,12,12,7,8,1,1,11,12],None,None,None,[1,2,3,4,5,6,7,8,9,10,11,12],None,[1,2,3,4,6,5,7,8,10,9,11,12]]

ChevieData["G24"]["ClassInfo"]={"classtext":ChevieData["G24"]["WordsClassRepresentatives"],
    "classnames":ChevieData["G24"]["ClassNames"],
    "classparams":ChevieData["G24"]["ClassNames"],
    "orders":[1,2,4,3,14,14,2,6,7,7,4,2],
    "classes":[1,21,42,56,24,24,21,56,24,24,42,1]}

def lxg244():
    res={"charparams":[[1,0],[1,21],[3,8],[3,1],[3,10],[3,3],[6,2],[6,9],[7,6],[7,3],[8,4],[8,5]],
        "opdam":Permutation("(11,12)"),
        "extRefl":[1,4,5,2]}
    res["b"]=map(lambda x: x[2-1],res["charparams"])
    return res

ChevieData["G24"]["CharInfo"]=lxg244

ChevieData["G24"]["CycPolSchurElements"]=[[1,0,2,2,2,3,4,6,7,14],[1,-21,2,2,2,3,4,6,7,14],[GAPMul(2,ER(-7)),-8,2,2,2,GAPDiv(1,7),GAPDiv(2,7),GAPDiv(4,7),GAPDiv(3,14),GAPDiv(5,14),GAPDiv(13,14)],[GAPMul(-2,ER(-7)),-1,2,2,2,GAPDiv(3,7),GAPDiv(5,7),GAPDiv(6,7),GAPDiv(1,14),GAPDiv(9,14),GAPDiv(11,14)],[GAPMul(-2,ER(-7)),-8,2,2,2,GAPDiv(3,7),GAPDiv(5,7),GAPDiv(6,7),GAPDiv(1,14),GAPDiv(9,14),GAPDiv(11,14)],[GAPMul(2,ER(-7)),-1,2,2,2,GAPDiv(1,7),GAPDiv(2,7),GAPDiv(4,7),GAPDiv(3,14),GAPDiv(5,14),GAPDiv(13,14)],[2,-1,2,4,7],[2,-8,2,4,7],[1,-6,2,2,2,3,4,6],[1,-3,2,2,2,3,4,6],[2,-4,3,7],[2,-4,3,7]]

ChevieData["G24"]["sparseFakeDegrees"]=[[1,0],[1,21],[1,8,1,16,1,18],[1,1,1,9,1,11],[1,10,1,12,1,20],[1,3,1,5,1,13],[1,2,1,4,1,6,1,8,1,10,1,12],[1,9,1,11,1,13,1,15,1,17,1,19],[1,6,1,8,1,10,1,12,1,14,1,16,1,18],[1,3,1,5,1,7,1,9,1,11,1,13,1,15],[1,4,1,6,1,8,1,10,1,12,2,14,1,16],[1,5,2,7,1,9,1,11,1,13,1,15,1,17]]

def lxg245(param,root):
    r=param[1-1][1-1]
    p=param[1-1][2-1]
    u=GetRoot(GAPMul(-p,r),2)
    def f1(r):
        return [1,r,r**2,r**2,r**9,r**3,r**4,r**7,r**18,r**6,r**11,r**21]
    
    
    def f3(p,r,a):
        return [3,GAPMul(2,p)+r,p**2,GAPMul(p,r)+p**2,GAPMul(GAPDiv(-1-a,2),p**6),GAPMul(GAPDiv(-1+a,2),p**2),GAPMul(-2,p**2)+p**4,0,GAPMul(GAPDiv(-1-a,2),p**12),GAPMul(GAPDiv(-1+a,2),p**4),GAPMul(-p**7,r**4),GAPMul(3,p**14)]
    
    
    def f6(r,p):
        return [6,GAPMul(2,p)+GAPMul(4,r),GAPMul(2,p)+GAPMul(2,r**2),GAPMul(2,p)+GAPMul(2,r**2),GAPMul(p**3,r**6),GAPMul(p,r**2),GAPMul(2,r**4),0,GAPMul(-p**6,r**12),GAPMul(-p**2,r**4),0,GAPMul(-6,p**7)]
    
    
    def f7(p,r):
        return [7,GAPMul(4,p)+GAPMul(3,r),GAPMul(2,p)+p**2,GAPMul(2,p)+GAPMul(2,p**2),0,0,GAPMul(-2,p**2)+p**4,GAPMul(p**4,r**3),0,0,GAPMul(-p**6,r**5),GAPMul(7,p**12)]
    
    
    def f8(p,r,u):
        return [8,GAPMul(4,p)+GAPMul(4,r),GAPMul(2,p)+p**2,GAPMul(3,p)+p**2,GAPMul(u,p**4),GAPMul(-p,r),GAPMul(-2,p**2)+p**4,GAPMul(p**3,r**3),GAPMul(-p**9,r**9),GAPMul(-p**3,r**3),0,GAPMul(8,p**10)]
    
    
    tbl={"identifier":"H(G24)",
        "name":"H(G24)",
        "size":336,
        "order":336,
        "powermap":ChevieData["G24"]["PowerMaps"],
        "irreducibles":GAPMul([f1(r),f1(p),f3(p,r,ER(-7)),f3(r,p,ER(-7)),f3(p,r,-ER(-7)),f3(r,p,-ER(-7)),f6(r,p),f6(p,r),f7(p,r),f7(r,p),f8(p,r,u),f8(p,r,-u)],p**0),
        "galomorphisms":PermutationGroup([Permutation("( 5, 6)( 9,10)")]),
        "irredinfo":map(lambda x: {"charparam":x,
        "charname":ChevieData["G24"]["CharName"](x,{})},ChevieData["G24"]["CharInfo"]()["charparams"])}
    tbl.update(ChevieData["G24"]["ClassInfo"])
    tbl["centralizers"]=map(lambda x: GAPDiv(tbl["size"],x),tbl["classes"])
    return ChevieData["compat"]["MakeCharacterTable"](tbl)

ChevieData["G24"]["HeckeCharTable"]=lxg245

def lxg246():
    return ChevieData["G24"]["HeckeCharTable"]([[1,-1],[1,-1],[1,-1]],[])

ChevieData["G24"]["CharTable"]=lxg246

def lxg247(para,root,i):
    p=para[1-1][2-1]
    r=para[1-1][1-1]
    f1=lambda r: map(lambda x: [[r]],range(1,3+1))
    def f3(p,r,a):
        return GAPMul(WGraph2Representation([[[2,3],[1,2],[1,3]],[[1,2,p,-r],[1,3,p,-r],[2,3,GAPDiv(GAPMul(r,1-a),2),GAPDiv(GAPMul(-p,a+1),2)]]],[p,r]),p**0)
    
    
    def f7(p,r):
        return GAPMul(WGraph2Representation([[[2,3],[2,3],[1,3],[1,3],[1,2],[1,2]],[[1,4,r,-p],[1,5,r,-p],[2,3,r,-p],[2,6,p,-r],[3,5,-p,0],[3,6,GAPMul(-2,p),r],[4,5,r,0],[4,6,GAPMul(2,r),0]]],[r,p]),p**0)
    
    
    def f9(r,p):
        return GAPMul(WGraph2Representation([[[1],[1,2],[1,3],[2],[2],[3],[3]],[[1,2,0,-r],[1,3,0,p],[1,4,p,-r],[1,5,0,-r],[1,6,-p,r],[2,5,-p,0],[2,7,-p,r],[3,4,-p,0],[3,5,p,-r],[3,6,p,0],[3,7,p,0],[4,6,0,-p],[4,7,-r,p],[5,6,-r,p],[5,7,-r,0]]],[p,r]),p**0)
    
    
    def f11(x,y,e):
        v=GAPMul(e,GetRoot(GAPMul(-x,y)))
        return [[[0,0,0,0,0,0,0,-x],[0,x+y,0,0,y,0,0,0],[0,0,x,GAPMul(-v,y)+GAPMul(x,y),0,0,-x**2,0],[0,0,0,y,0,0,0,0],[0,-x,0,0,0,0,0,0],[0,0,0,x,0,x,-v-y,0],[0,0,0,0,0,0,y,0],[y,0,0,0,0,0,0,x+y]],[[x,0,0,v,0,0,0,-y],[0,x,0,v,x,0,0,0],[0,0,x+y,0,0,0,GAPMul(-x,y),0],[0,0,0,y,0,0,0,0],[0,0,0,0,y,0,0,0],[0,0,-1,x,-v,x,x,v],[0,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,y]],[[y,0,0,0,0,0,0,0],[0,x,0,0,x,0,-v,0],[GAPMul(-x,y),0,x,0,GAPMul(-v,y),GAPMul(v,y),GAPMul(v,y)-GAPMul(x,y)-x**2,0],[0,0,0,x,0,-y,-v-y,0],[0,0,0,0,y,0,0,0],[0,0,0,0,0,y,0,0],[0,0,0,0,0,0,y,0],[x,0,0,0,0,0,x,x]]]
    
    
    rep=[[f1,r],[f1,p],[f3,p,r,ER(-7)],[f3,r,p,ER(-7)],[f3,p,r,-ER(-7)],[f3,r,p,-ER(-7)],[f7,p,r],[f7,r,p],[f9,p,r],[f9,r,p],[f11,p,r,1],[f11,p,r,-1]]
    return rep[i-1][1-1](*[rep[i-1][k-1] for k in range(2,len(rep[i-1])+1)])+GAPMul(0,prod(para[1-1]))

ChevieData["G24"]["HeckeRepresentation"]=lxg247

ChevieData["families"]["X7"]={"name":"X7",
    "fourierMat":GAPDiv([[GAPDiv(-1,2),GAPDiv(1,2),GAPDiv(ER(-7),2),GAPDiv(ER(-7),2),-1,-1,-1],[GAPDiv(1,2),GAPDiv(-1,2),GAPDiv(ER(-7),2),GAPDiv(ER(-7),2),1,1,1],[GAPDiv(ER(-7),2),GAPDiv(ER(-7),2),GAPDiv(ER(-7),2),GAPDiv(-ER(-7),2),0,0,0],[GAPDiv(ER(-7),2),GAPDiv(ER(-7),2),GAPDiv(-ER(-7),2),GAPDiv(ER(-7),2),0,0,0],[-1,1,0,0,-ER(7)**6-ER(7),-ER(7)**5-ER(7)**2,-ER(7)**4-ER(7)**3],[-1,1,0,0,-ER(7)**5-ER(7)**2,-ER(7)**4-ER(7)**3,-ER(7)**6-ER(7)],[-1,1,0,0,-ER(7)**4-ER(7)**3,-ER(7)**6-ER(7),-ER(7)**5-ER(7)**2]],ER(-7)),
    "eigenvalues":[1,1,1,-1,ER(7)**4,ER(7)**2,ER(7)],
    "explanation":"mystery G24",
    "special":1,
    "cospecial":2}

def lxg248():
    return {"harishChandra":[{"relativeType":{"series":"ST",
        "indices":range(1,3+1),
        "rank":3,
        "ST":24},
        "levi":[],
        "parameterExponents":[1,1,1],
        "charNumbers":range(1,12+1),
        "eigenvalue":1,
        "cuspidalName":""},{"relativeType":{"series":"A",
        "indices":[1],
        "rank":1},
        "levi":[2,3],
        "parameterExponents":[7],
        "charNumbers":[19,13],
        "eigenvalue":-1,
        "cuspidalName":"B_2"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[17],
        "eigenvalue":ER(4),
        "qEigen":GAPDiv(1,2),
        "cuspidalName":"G_{24}[i]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[18],
        "eigenvalue":-ER(4),
        "qEigen":GAPDiv(1,2),
        "cuspidalName":"G_{24}[-i]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[20],
        "eigenvalue":ER(7)**3,
        "cuspidalName":"G_{24}[\\zeta_7^3]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[21],
        "eigenvalue":ER(7)**5,
        "cuspidalName":"G_{24}[\\zeta_7^5]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[22],
        "eigenvalue":ER(7)**6,
        "cuspidalName":"G_{24}[\\zeta_7^6]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[14],
        "eigenvalue":ER(7)**4,
        "cuspidalName":"G_{24}[\\zeta_7^4]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[15],
        "eigenvalue":ER(7)**2,
        "cuspidalName":"G_{24}[\\zeta_7^2]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[16],
        "eigenvalue":ER(7),
        "cuspidalName":"G_{24}[\\zeta_7]"}],
        "families":[Family("C1",[1]),Family("X7",[4,6,7,13,14,15,16]),Family("C1",[10]),Family("C'\"2",[11,12,17,18]),Family("C1",[9]),ComplexConjugate(Family("X7",[3,5,8,19,20,21,22])),Family("C1",[2])],
        "a":[0,21,8,1,8,1,1,8,6,3,4,4,1,1,1,1,4,4,8,8,8,8],
        "A":[0,21,20,13,20,13,13,20,18,15,17,17,13,13,13,13,17,17,20,20,20,20],
        "curtis":[2,1,6,5,4,3,8,7,10,9,12,11,19,-20,-21,-22,-18,-17,13,-14,-15,-16]}

ChevieData["G24"]["UnipotentCharacters"]=lxg248

def lxg249(x,y,z):
    return GAPMul(-42,x**2)-GAPMul(12,x**2)+GAPMul(GAPDiv(21,2),x**2)-GAPMul(GAPDiv(9,2),y**2)-GAPMul(6,y**3)+GAPMul(14,x**4)-GAPMul(GAPDiv(21,8),z**4)

def lxg2410(x,y,z):
    return GAPMul(-1960,x**2)+GAPMul(840,x**2)-GAPMul(1120,x**2)+GAPMul(1760,x**2)-GAPMul(1225,x**2)+GAPMul(525,y**2)-GAPMul(280,y**3)+GAPMul(3920,x**4)-GAPMul(980,x**4)-GAPMul(180,y**4)-GAPMul(240,y**5)+GAPMul(1568,x**6)-GAPMul(GAPDiv(416,7),y**6)-GAPMul(GAPDiv(49,2),z**6)

def lxg2411(x,y,z):
    return GAPMul(GAPDiv(-857157,4),x**2)-GAPMul(GAPDiv(4847619,32),x**2)+GAPMul(GAPDiv(1596665,8),x**2)-GAPMul(440118,x**2)+GAPMul(1608075,x**2)-GAPMul(633080,x**2)+GAPMul(269760,x**2)-GAPMul(GAPDiv(1327753,128),x**2)+GAPMul(GAPDiv(569037,128),y**2)-GAPMul(GAPDiv(122451,4),y**3)-GAPMul(GAPDiv(11176655,16),x**4)+GAPMul(432180,x**4)-GAPMul(2088870,x**4)-GAPMul(2922360,x**4)-GAPMul(24696,x**4)-GAPMul(5735940,x**4)-GAPMul(4210080,x**4)+GAPMul(2688840,x**4)-GAPMul(203456,x**4)+GAPMul(GAPDiv(11311111,64),x**4)-GAPMul(16696554,x**6)+GAPMul(3755850,x**6)-GAPMul(470400,x**6)+GAPMul(2546880,x**6)-GAPMul(GAPDiv(396459,8),y**6)+GAPMul(76734,y**7)-GAPMul(3457440,x**8)-GAPMul(2511936,x**8)-GAPMul(GAPDiv(8621991,4),x**8)-GAPMul(GAPDiv(424809,4),y**8)-GAPMul(114513,y**9)+GAPMul(9008552,x**10)-GAPMul(2304960,x**10)+GAPMul(7222208,x**10)-GAPMul(8978368,x**10)+GAPMul(6537923,x**10)-GAPMul(40392,y**10)+GAPMul(14928,y**11)-GAPMul(537824,x**12)-GAPMul(153664,x**12)+GAPMul(134456,x**12)

ChevieData["G24"]["Invariants"]=[lxg249,lxg2410,lxg2411]

def lxg2412():
    return lxg2413

def lxg2413(x,y,z):
    return [[x,GAPMul(3,y**2),GAPMul(7,z)-GAPMul(9,x**2)],[GAPMul(3,y),GAPMul(1792,z),GAPMul(64,x)+GAPMul(3136,x**4)],[GAPMul(7,z),GAPMul(64,x)+GAPMul(5376,x**2),GAPMul(GAPDiv(287,2),x)-GAPMul(GAPDiv(35,4),x**3)+GAPMul(GAPDiv(21,256),y**4)-GAPMul(1568,x**6)]]

ChevieData["G24"]["BasicDerivations"]=lxg2412

def lxg2414():
    return lxg2415

def lxg2415(x,y,z):
    return GAPMul(18,x)+GAPMul(5632,x**2)-GAPMul(1024,z**3)-GAPMul(67,x**3)-GAPMul(4352,x**4)-GAPMul(5504,x**6)-GAPMul(GAPDiv(27,3136),y**7)-GAPMul(229376,x**7)-GAPMul(114688,x**9)

ChevieData["G24"]["Discriminant"]=lxg2414
