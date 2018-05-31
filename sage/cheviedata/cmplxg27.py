
ChevieData["G27"]["AltPres"]=[{"gens":[[2],[2,3,-2],[1]],
    "rels":[[[1,3,1],[3,1,3]],[[1,2,1,2],[2,1,2,1]],[[2,3,2,3,2],[3,2,3,2,3]],[[3,2,1,3,2,1,2],[1,3,2,1,3,2,1]]]},{"gens":[[1],[-1,3,1],[2]],
    "rels":[[[1,2,1],[2,1,2]],[[2,3,2,3,2],[3,2,3,2,3]],[[1,3,1],[3,1,3]],[[2,3,2,3,1,2,3,1,2,3],[3,2,3,1,2,3,1,2,3,1]]]},{"gens":[[2],[3],[-3,-2,1,2,3]],
    "rels":[[[1,2,1,2],[2,1,2,1]],[[2,3,2,3,2],[3,2,3,2,3]],[[1,3,1,3,1],[3,1,3,1,3]],[[2,3,1,2,3,1,2],[1,2,3,1,2,3,1]],[[3,1,2,3,1,2,3,1],[1,2,3,1,2,3,1,3]]]},{"gens":[[3],[-3,2,3],[1]],
    "rels":[[[1,3,1],[3,1,3]],[[1,2,1,2],[2,1,2,1]],[[2,3,2,3,2],[3,2,3,2,3]],[[2,1,2,3,1,2,3],[1,2,3,1,2,3,1]]]}]

ChevieData["G27"]["BraidRelations"]=[[[2,1,2],[1,2,1]],[[3,1,3],[1,3,1]],[[3,2,3,2],[2,3,2,3]],[[3,2,3,1,2,3,1,2,3,1,2,3],[2,3,1,2,3,1,2,3,1,2,3,2]]]

def lxg271(indices,title):
    print title," ",indices[1-1],"\n",
    s=just("",len(title)-1)
    print s," / \\\n",s,indices[2-1],"=====",indices[3-1],"  ",IntListToString([indices[k-1] for k in [3,2,3,1,2,3,1,2,3,1,2,3]]),"==",IntListToString([indices[k-1] for k in [2,3,1,2,3,1,2,3,1,2,3,2]]),"\n",

ChevieData["G27"]["PrintDiagram"]=lxg271

ChevieData["G27"]["GeneratingRoots"]=[[GAPDiv(-ER(5),5),GAPDiv(-ER(15)**14+ER(15)**13-ER(15)**8-ER(15)**7-ER(15)**4+ER(15)**2,5),GAPDiv(GAPMul(4,ER(15)**14)-GAPMul(4,ER(15)**11)+GAPMul(2,ER(15)**4)-GAPMul(2,ER(15)),15)],[GAPDiv(GAPMul(5-ER(-15),ER(3)),10),GAPDiv(-ER(15)**13-ER(15)**8+ER(15)**7-ER(15),5),GAPDiv(GAPMul(-4,ER(15)**14)+GAPMul(4,ER(15)**11)-GAPMul(2,ER(15)**4)+GAPMul(2,ER(15)),15)],[GAPDiv(GAPMul(2,ER(5)),5),0,0]]

ChevieData["G27"]["GeneratingCoRoots"]=[[GAPDiv(-ER(5),2),GAPDiv(ER(15)**14+ER(15)**13-ER(15)**11-ER(15)**8-ER(15)**7+ER(15)**4-ER(15),2),GAPDiv(-ER(15)**14+ER(15)**11-GAPMul(2,ER(15)**4)+GAPMul(2,ER(15)),2)],[GAPDiv(GAPMul(5+ER(-15),ER(3)**2),4),GAPDiv(-ER(15)**14+ER(15)**13-ER(15)**7-ER(15)**2,2),GAPDiv(ER(15)**14-ER(15)**11+GAPMul(2,ER(15)**4)-GAPMul(2,ER(15)),2)],[ER(5),0,0]]

def lxg272():
    return GAPMul(ChevieData["G27"]["GeneratingCoRoots"],Matrix(ChevieData["G27"]["GeneratingRoots"]).transpose())

ChevieData["G27"]["CartanMat"]=lxg272

ChevieData["G27"]["EigenvaluesGeneratingReflections"]=[GAPDiv(1,2),GAPDiv(1,2),GAPDiv(1,2)]

ChevieData["G27"]["Size"]=2160

ChevieData["G27"]["ReflectionDegrees"]=[6,12,30]

ChevieData["G27"]["NrConjugacyClasses"]=34

def lxg273(s):
    t=[[[]],[[1]],[[1,2],[1,3],[2,3],[1,8]],[range(1,3+1)]]
    return t[s+1-1]

ChevieData["G27"]["ParabolicRepresentatives"]=lxg273

ChevieData["G27"]["ClassNames"]=[".","2","12","c2","23","132","c","2c2c","13","3232","2cc","ccc","23z","132132","cc","13z","2ccc","32zzzz","cc21323","zcc","12z","1z","3232zz","2323z","1zz","zzcccc","cccc","zcccc","32zzz","z","zzzzz","zz","zzzz","zzz"]

ChevieData["G27"]["WordsClassRepresentatives"]=map(lambda x: Replace(x,".",[],"1",[1],"2",[2],"3",[3],"c",[1,2,3],"z",[1,2,3,1,2,3,1,2,3,1,2,3,1,2,3]),ChevieData["G27"]["ClassNames"])

ChevieData["G27"]["PowerMaps"]=[None,[1,1,3,8,10,14,15,4,9,1,22,4,23,26,27,9,22,23,27,26,3,32,33,32,33,14,15,8,10,32,33,33,32,1],[1,2,1,8,5,12,12,4,1,10,29,28,29,4,4,34,5,5,28,28,34,10,10,2,2,8,8,12,29,34,34,1,1,34],None,[1,2,3,1,5,31,30,1,9,10,13,34,11,33,32,16,18,17,31,30,21,23,22,25,24,32,33,34,29,31,30,33,32,34],None,[1,2,3,8,5,19,20,4,9,10,11,28,13,27,26,16,17,18,6,7,21,22,23,24,25,15,14,12,29,30,31,32,33,34],None,None,None,[1,2,3,4,5,7,6,8,9,10,13,12,11,15,14,16,18,17,20,19,21,23,22,25,24,27,26,28,29,31,30,33,32,34],None,[1,2,3,8,5,19,20,4,9,10,11,28,13,27,26,16,17,18,6,7,21,22,23,24,25,15,14,12,29,30,31,32,33,34],None,None,None,[1,2,3,8,5,20,19,4,9,10,13,28,11,26,27,16,18,17,7,6,21,23,22,25,24,14,15,12,29,31,30,33,32,34],None,[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34],None,None,None,[1,2,3,8,5,20,19,4,9,10,13,28,11,26,27,16,18,17,7,6,21,23,22,25,24,14,15,12,29,31,30,33,32,34],None,None,None,None,None,[1,2,3,4,5,7,6,8,9,10,13,12,11,15,14,16,18,17,20,19,21,23,22,25,24,27,26,28,29,31,30,33,32,34]]

ChevieData["G27"]["ClassInfo"]={"classtext":ChevieData["G27"]["WordsClassRepresentatives"],
    "classnames":ChevieData["G27"]["ClassNames"],
    "classparams":ChevieData["G27"]["ClassNames"],
    "orders":[1,2,3,5,4,30,30,5,3,2,12,10,12,15,15,6,12,12,30,30,6,6,6,6,6,15,15,10,4,6,6,3,3,2],
    "classes":[1,45,120,72,90,72,72,72,120,45,90,72,90,72,72,120,90,90,72,72,120,45,45,45,45,72,72,72,90,1,1,1,1,1]}

def lxg274():
    res={"charparams":[[1,0],[1,45],[3,7],[3,22],[3,1],[3,16],[3,5,1],[3,20,1],[3,5,2],[3,20,2],[5,6,2],[5,6,1],[5,15,2],[5,15,1],[6,19],[6,4],[6,17],[6,2],[8,6],[8,9,1],[8,12],[8,9,2],[9,6],[9,9],[9,13],[9,4],[9,11],[9,8],[10,12],[10,3],[15,7],[15,10],[15,5],[15,8]],
        "opdam":Permutation("(19,20)(21,22)(23,26,28)(24,27,25)"),
        "extRefl":[1,5,10,2]}
    res["b"]=map(lambda x: x[2-1],res["charparams"])
    return res

ChevieData["G27"]["CharInfo"]=lxg274

ChevieData["G27"]["CycPolSchurElements"]=[[1,0,2,2,2,3,3,3,4,5,6,6,6,10,12,15,30],[1,-45,2,2,2,3,3,3,4,5,6,6,6,10,12,15,30],[GAPMul(-2,ER(-15)),-1,2,2,2,GAPDiv(1,3),GAPDiv(1,3),GAPDiv(1,3),GAPDiv(1,5),GAPDiv(4,5),GAPDiv(5,6),GAPDiv(5,6),GAPDiv(5,6),GAPDiv(1,10),GAPDiv(9,10),GAPDiv(2,15),GAPDiv(8,15),GAPDiv(7,30),GAPDiv(13,30)],[GAPMul(2,ER(-15)),-16,2,2,2,GAPDiv(2,3),GAPDiv(2,3),GAPDiv(2,3),GAPDiv(1,5),GAPDiv(4,5),GAPDiv(1,6),GAPDiv(1,6),GAPDiv(1,6),GAPDiv(1,10),GAPDiv(9,10),GAPDiv(7,15),GAPDiv(13,15),GAPDiv(17,30),GAPDiv(23,30)],[GAPMul(2,ER(-15)),-1,2,2,2,GAPDiv(1,3),GAPDiv(1,3),GAPDiv(1,3),GAPDiv(2,5),GAPDiv(3,5),GAPDiv(5,6),GAPDiv(5,6),GAPDiv(5,6),GAPDiv(3,10),GAPDiv(7,10),GAPDiv(11,15),GAPDiv(14,15),GAPDiv(1,30),GAPDiv(19,30)],[GAPMul(-2,ER(-15)),-16,2,2,2,GAPDiv(2,3),GAPDiv(2,3),GAPDiv(2,3),GAPDiv(2,5),GAPDiv(3,5),GAPDiv(1,6),GAPDiv(1,6),GAPDiv(1,6),GAPDiv(3,10),GAPDiv(7,10),GAPDiv(1,15),GAPDiv(4,15),GAPDiv(11,30),GAPDiv(29,30)],[GAPMul(2,ER(-15)),-1,2,2,2,GAPDiv(2,3),GAPDiv(2,3),GAPDiv(2,3),GAPDiv(1,5),GAPDiv(4,5),GAPDiv(1,6),GAPDiv(1,6),GAPDiv(1,6),GAPDiv(1,10),GAPDiv(9,10),GAPDiv(7,15),GAPDiv(13,15),GAPDiv(17,30),GAPDiv(23,30)],[GAPMul(-2,ER(-15)),-16,2,2,2,GAPDiv(1,3),GAPDiv(1,3),GAPDiv(1,3),GAPDiv(1,5),GAPDiv(4,5),GAPDiv(5,6),GAPDiv(5,6),GAPDiv(5,6),GAPDiv(1,10),GAPDiv(9,10),GAPDiv(2,15),GAPDiv(8,15),GAPDiv(7,30),GAPDiv(13,30)],[GAPMul(-2,ER(-15)),-1,2,2,2,GAPDiv(2,3),GAPDiv(2,3),GAPDiv(2,3),GAPDiv(2,5),GAPDiv(3,5),GAPDiv(1,6),GAPDiv(1,6),GAPDiv(1,6),GAPDiv(3,10),GAPDiv(7,10),GAPDiv(1,15),GAPDiv(4,15),GAPDiv(11,30),GAPDiv(29,30)],[GAPMul(2,ER(-15)),-16,2,2,2,GAPDiv(1,3),GAPDiv(1,3),GAPDiv(1,3),GAPDiv(2,5),GAPDiv(3,5),GAPDiv(5,6),GAPDiv(5,6),GAPDiv(5,6),GAPDiv(3,10),GAPDiv(7,10),GAPDiv(11,15),GAPDiv(14,15),GAPDiv(1,30),GAPDiv(19,30)],[2,-3,2,2,2,3,3,3,6,6,6],[2,-3,2,2,2,3,3,3,6,6,6],[2,-12,2,2,2,3,3,3,6,6,6],[2,-12,2,2,2,3,3,3,6,6,6],[3-ER(-3),-16,2,GAPDiv(2,3),GAPDiv(2,3),GAPDiv(2,3),4,5,GAPDiv(1,6),GAPDiv(5,12),GAPDiv(11,12),GAPDiv(1,15),GAPDiv(4,15),GAPDiv(7,15),GAPDiv(13,15)],[3+ER(-3),-1,2,GAPDiv(1,3),GAPDiv(1,3),GAPDiv(1,3),4,5,GAPDiv(5,6),GAPDiv(1,12),GAPDiv(7,12),GAPDiv(2,15),GAPDiv(8,15),GAPDiv(11,15),GAPDiv(14,15)],[3+ER(-3),-16,2,GAPDiv(1,3),GAPDiv(1,3),GAPDiv(1,3),4,5,GAPDiv(5,6),GAPDiv(1,12),GAPDiv(7,12),GAPDiv(2,15),GAPDiv(8,15),GAPDiv(11,15),GAPDiv(14,15)],[3-ER(-3),-1,2,GAPDiv(2,3),GAPDiv(2,3),GAPDiv(2,3),4,5,GAPDiv(1,6),GAPDiv(5,12),GAPDiv(11,12),GAPDiv(1,15),GAPDiv(4,15),GAPDiv(7,15),GAPDiv(13,15)],[5+ER(5),-6,3,3,3,GAPDiv(2,5),GAPDiv(3,5),GAPDiv(1,15),GAPDiv(4,15),GAPDiv(11,15),GAPDiv(14,15)],[5+ER(5),-6,3,3,3,GAPDiv(2,5),GAPDiv(3,5),GAPDiv(1,15),GAPDiv(4,15),GAPDiv(11,15),GAPDiv(14,15)],[5-ER(5),-6,3,3,3,GAPDiv(1,5),GAPDiv(4,5),GAPDiv(2,15),GAPDiv(7,15),GAPDiv(8,15),GAPDiv(13,15)],[5-ER(5),-6,3,3,3,GAPDiv(1,5),GAPDiv(4,5),GAPDiv(2,15),GAPDiv(7,15),GAPDiv(8,15),GAPDiv(13,15)],[3,-4,2,2,2,4,5,10],[3,-9,2,2,2,4,5,10],[3,-9,2,2,2,4,5,10],[3,-4,2,2,2,4,5,10],[3,-9,2,2,2,4,5,10],[3,-4,2,2,2,4,5,10],[2,-12,2,3,3,3,4,6,12],[2,-3,2,3,3,3,4,6,12],[GAPDiv(3+ER(-3),2),-5,2,2,2,GAPDiv(2,3),GAPDiv(2,3),GAPDiv(2,3),4,GAPDiv(1,6),GAPDiv(1,6),GAPDiv(1,6),GAPDiv(5,12),GAPDiv(11,12)],[GAPDiv(3-ER(-3),2),-8,2,2,2,GAPDiv(1,3),GAPDiv(1,3),GAPDiv(1,3),4,GAPDiv(5,6),GAPDiv(5,6),GAPDiv(5,6),GAPDiv(1,12),GAPDiv(7,12)],[GAPDiv(3-ER(-3),2),-5,2,2,2,GAPDiv(1,3),GAPDiv(1,3),GAPDiv(1,3),4,GAPDiv(5,6),GAPDiv(5,6),GAPDiv(5,6),GAPDiv(1,12),GAPDiv(7,12)],[GAPDiv(3+ER(-3),2),-8,2,2,2,GAPDiv(2,3),GAPDiv(2,3),GAPDiv(2,3),4,GAPDiv(1,6),GAPDiv(1,6),GAPDiv(1,6),GAPDiv(5,12),GAPDiv(11,12)]]

ChevieData["G27"]["sparseFakeDegrees"]=[[1,0],[1,45],[1,7,1,13,1,25],[1,22,1,28,1,40],[1,1,1,19,1,25],[1,16,1,34,1,40],[1,5,1,17,1,23],[1,20,1,32,1,38],[1,5,1,11,1,29],[1,20,1,26,1,44],[1,6,1,12,1,18,1,24,1,30],[1,6,1,12,1,18,1,24,1,30],[1,15,1,21,1,27,1,33,1,39],[1,15,1,21,1,27,1,33,1,39],[1,19,2,25,1,31,1,37,1,43],[1,4,2,10,1,16,1,22,1,28],[1,17,1,23,1,29,2,35,1,41],[1,2,1,8,1,14,2,20,1,26],[1,6,1,12,1,18,2,24,2,30,1,36],[1,9,2,15,2,21,1,27,1,33,1,39],[2,12,2,18,1,24,2,30,1,36],[1,9,2,15,1,21,2,27,2,33],[1,6,2,12,2,18,2,24,1,30,1,36],[1,9,1,15,2,21,2,27,2,33,1,39],[2,13,2,19,1,25,2,31,2,37],[1,4,1,10,2,16,2,22,2,28,1,34],[1,11,2,17,2,23,2,29,1,35,1,41],[2,8,2,14,1,20,2,26,2,32],[1,12,2,18,2,24,2,30,2,36,1,42],[1,3,2,9,2,15,2,21,2,27,1,33],[2,7,3,13,3,19,3,25,3,31,1,37],[2,10,3,16,3,22,3,28,3,34,1,40],[1,5,3,11,3,17,3,23,3,29,2,35],[1,8,3,14,3,20,3,26,3,32,2,38]]

def lxg275(para,root):
    q=GAPDiv(-para[1-1][1-1],para[1-1][2-1])
    r=para[1-1][1-1]
    p=para[1-1][2-1]
    tbl={"identifier":"H(G27)",
        "size":2160,
        "order":2160,
        "powermap":ChevieData["G27"]["PowerMaps"]}
    tbl.update(ChevieData["G27"]["ClassInfo"])
    f1=lambda r: map(lambda x: r**len(x),tbl["classtext"])
    def f3(r,p,j,d):
        return [3,p+GAPMul(2,r),GAPMul(p,r)+r**2,GAPMul(-1+GAPDiv(GAPMul(j,1-d),2),p)+GAPMul(GAPDiv(-2+GAPMul(1-d,1+j),2),p**2),r**2,GAPMul(-p,r**2)+GAPMul(GAPDiv(-1+j-d-GAPMul(j,d),2),p**2),GAPMul(GAPDiv(GAPMul(1-d,j),2),p),GAPMul(GAPDiv(3,2)-GAPMul(GAPDiv(3,2),j)+GAPMul(GAPDiv(1,2),d),p**2)+GAPMul(1-j+d,p**3),GAPMul(p,r)+r**2,GAPMul(-2,p**2)+r**4,GAPMul(-j,p**2),GAPMul(GAPDiv(1+d,2),p**3),GAPMul(j**2,p**5),GAPMul(GAPDiv(j,2)+d,r**2)+GAPMul(1-j+d,p**3),GAPMul(GAPDiv(GAPMul(1+d,j**2),2),p**2),GAPMul(j**2,p**5)+GAPMul(j**2,p**6),GAPMul(GAPDiv(1+d,2),p**3)+GAPMul(GAPDiv(1,2)+j,p**4),GAPMul(j**2,p**20),GAPMul(GAPDiv(GAPMul(-1+d,j**2),2),p**3)+GAPMul(j**2,p**5),GAPMul(GAPDiv(GAPMul(1+d,j),2),p**7),GAPMul(j**2,p**5)+GAPMul(j**2,p**6),GAPMul(2,j**2)+GAPMul(j**2,p**6),GAPMul(j,p**10)-GAPMul(2,j),GAPMul(j**2,p**5)-GAPMul(2,j**2),GAPMul(2,j)+GAPMul(j,p**11),GAPMul(GAPDiv(GAPMul(1-d,j**2),2),p**14),GAPMul(GAPDiv(GAPMul(1-d,j),2),p**4),GAPMul(GAPDiv(1-d,2),p**9),GAPMul(p**15,r**32),GAPMul(3,j**2),GAPMul(3,j),GAPMul(3,j),GAPMul(3,j**2),GAPMul(3,p**15)]
    
    
    def f11(r,p):
        return [5,GAPMul(2,p)+GAPMul(3,r),GAPMul(2,p)+r**2,GAPMul(-p,r**3)-GAPMul(p**2,r**2),GAPMul(2,p)+r**2,0,0,GAPMul(p**2,r**6)+GAPMul(2,p**3),GAPMul(p,r)+p**2,r**4,GAPMul(p**3,r**4),0,GAPMul(p**6,r**11)+GAPMul(2,p**7),0,0,GAPMul(2,p**6)+GAPMul(p**7,r**10),GAPMul(-p**4,r**6),GAPMul(p**24,r**38)+GAPMul(2,p**25),GAPMul(-p**4,r**7)+GAPMul(p**6,r**5),0,GAPMul(p**6,r**11)+GAPMul(2,p**7),GAPMul(3,p**6)+GAPMul(2,p**7),GAPMul(p**12,r**22),GAPMul(p**6,r**13),GAPMul(3,p**12)+GAPMul(2,p**13),0,0,0,GAPMul(p**18,r**29)+GAPMul(2,p**19),GAPMul(5,p**6),GAPMul(5,p**30),GAPMul(5,p**12),GAPMul(5,p**24),GAPMul(5,p**18)]
    
    
    def f12(r,p):
        return [5,GAPMul(2,p)+GAPMul(3,r),GAPMul(p,r)+p**2,GAPMul(-p,r**3)-GAPMul(p**2,r**2),GAPMul(2,p)+r**2,0,0,GAPMul(p**2,r**6)+GAPMul(2,p**3),GAPMul(2,p)+r**2,r**4,GAPMul(p**3,r**4),0,GAPMul(p**6,r**11)+GAPMul(2,p**7),0,0,GAPMul(p**6,r**11)+GAPMul(2,p**7),GAPMul(-p**4,r**6),GAPMul(p**24,r**38)+GAPMul(2,p**25),GAPMul(-p**4,r**7)+GAPMul(p**6,r**5),0,GAPMul(2,p**6)+GAPMul(p**7,r**10),GAPMul(3,p**6)+GAPMul(2,p**7),GAPMul(p**12,r**22),GAPMul(p**6,r**13),GAPMul(3,p**12)+GAPMul(2,p**13),0,0,0,GAPMul(p**18,r**29)+GAPMul(2,p**19),GAPMul(5,p**6),GAPMul(5,p**30),GAPMul(5,p**12),GAPMul(5,p**24),GAPMul(5,p**18)]
    
    
    def f15(r,p,j):
        return [6,GAPMul(4,p)+GAPMul(2,r),GAPMul(2,p)+GAPMul(2,p**2),GAPMul(j**2,p**2)+GAPMul(GAPMul(2,j**2)+j,p**3),GAPMul(2,p)+GAPMul(2,p**2),GAPMul(j**2,p),GAPMul(-j,p**2),GAPMul(j,p**4)-GAPMul(2,j**2)-GAPMul(GAPMul(3,j**2)+GAPMul(2,j),p**6),GAPMul(2,p)+GAPMul(2,p**2),GAPMul(2,p**4),0,GAPMul(-p**6,r**3),GAPMul(-2,j**2)-GAPMul(2,j**2),GAPMul(j,p**2)-GAPMul(2,j**2)+GAPMul(2,p**4)-GAPMul(2,j),GAPMul(j**2,p**4),GAPMul(-2,j**2)-GAPMul(2,j**2),GAPMul(-p**6,r**4)-GAPMul(p**7,r**3),GAPMul(2,j**2)+GAPMul(2,j**2),GAPMul(-2,j**2)+GAPMul(j**2,p**8),GAPMul(-j,p**14),GAPMul(-2,j**2)-GAPMul(2,j**2),GAPMul(-2,j**2)-GAPMul(4,j**2),GAPMul(2,j),GAPMul(-2,j**2),GAPMul(2,j)+GAPMul(4,j),GAPMul(j**2,p**28),GAPMul(j,p**8),GAPMul(-p**18,r**9),GAPMul(-2,p**31)-GAPMul(2,p**32),GAPMul(-6,j**2),GAPMul(-6,j),GAPMul(6,j),GAPMul(6,j**2),GAPMul(-6,p**30)]
    
    
    def f19(x,y,d,sgn):
        v=GAPMul(sgn,GetRoot(GAPMul(-x,y),2))
        return [8,GAPMul(4,x)+GAPMul(4,y),GAPMul(3,x)+x**2,GAPMul(GAPDiv(-1+d,2),v)-GAPMul(GAPDiv(1+d,2),v)-GAPMul(x,y**3)+GAPMul(GAPDiv(-3+d,2),x**2)-GAPMul(x**3,y),GAPMul(2,x)+x**2,GAPMul(GAPDiv(-1+d,2),v)-GAPMul(x,y**2)-GAPMul(x**2,y)-x**3-y**3,GAPMul(GAPDiv(-1+d,2),v),GAPMul(x**2,y**6)+GAPMul(GAPDiv(1-d,2),x**3)-GAPMul(GAPDiv(1+GAPMul(3,d),2),x**4)+GAPMul(GAPDiv(1-d,2),x**5),GAPMul(3,x)+x**2,GAPMul(-2,x**2)+x**4,0,GAPMul(GAPDiv(1-d,2),v),GAPMul(-v,x**7)-GAPMul(2,v)-GAPMul(v,x**9),GAPMul(3+d,v)+GAPMul(3+d,v),GAPMul(GAPDiv(-1+d,2),x**3),GAPMul(-v,x**7)-GAPMul(3,v)-GAPMul(v,x**9),GAPMul(GAPDiv(1-d,2),v)+GAPMul(GAPDiv(1-d,2),v),GAPMul(x**30,y**32)+GAPMul(2,x**31),GAPMul(-v,x**3)+GAPMul(GAPDiv(3+d,2),v)-GAPMul(v,x**7),GAPMul(GAPDiv(1-d,2),v),GAPMul(-v,x**7)-GAPMul(3,v)-GAPMul(v,x**9),GAPMul(-4,v)-GAPMul(4,v),GAPMul(-x**15,y**19)+GAPMul(2,x**17)-GAPMul(x**19,y**15),GAPMul(-v,x**7)+GAPMul(2,v)-GAPMul(v,x**11),GAPMul(-4,x**15)-GAPMul(4,x**16),GAPMul(GAPDiv(-1+d,2),x**21),GAPMul(GAPDiv(1+d,2),x**6),GAPMul(GAPDiv(-1+d,2),v),GAPMul(v,x**22)+GAPMul(2,v),GAPMul(-8,v),GAPMul(-8,v),GAPMul(-8,x**15),GAPMul(8,x**30),GAPMul(8,v)]
    
    
    def f23(r,p,j):
        q=GAPMul(GetRoot(GAPDiv(-r,p),3),j)
        return Zip([9,-4+GAPMul(5,q**3),1-GAPMul(3,q**3)+GAPMul(2,q**6),q**5-GAPMul(2,q**6)-q**8+q**9,1-GAPMul(2,q**3)+GAPMul(2,q**6),-q**3-q**4+q**6,-q**5,q**10-GAPMul(2,q**11)-GAPMul(2,q**12)-GAPMul(2,q**13)+GAPMul(2,q**14),1-GAPMul(3,q**3)+GAPMul(2,q**6),1-GAPMul(2,q**6)+GAPMul(2,q**12),q**11,-q**15,q**25-GAPMul(2,q**28)+GAPMul(2,q**31),q**6+GAPMul(2,q**7)-GAPMul(2,q**9)-GAPMul(4,q**10)+GAPMul(3,q**12)-GAPMul(2,q**14)-GAPMul(2,q**15),-q**10,q**25-GAPMul(3,q**28)+GAPMul(2,q**31),q**15+q**17-q**18,q**100-GAPMul(2,q**103)+GAPMul(2,q**106),q**13-GAPMul(2,q**19),-q**35,q**25-GAPMul(3,q**28)+GAPMul(2,q**31),GAPMul(-4,q**25)+GAPMul(5,q**28),q**50-GAPMul(2,q**56)+GAPMul(2,q**62),q**25-GAPMul(2,q**31)+GAPMul(2,q**37),GAPMul(-4,q**50)+GAPMul(5,q**53),-q**70,-q**20,-q**45,q**75-GAPMul(2,q**78)+GAPMul(2,q**81),GAPMul(9,q**25),GAPMul(9,q**125),GAPMul(9,q**50),GAPMul(9,q**100),GAPMul(9,q**75)],tbl["classtext"],lxg276)
    
    
    def f29(r,p):
        q=GAPDiv(-r,p)
        return Zip([10,-6+GAPMul(4,q),3-GAPMul(3,q)+q**2,GAPMul(2,q)-GAPMul(2,q**2),2-GAPMul(2,q),GAPMul(-2,q)+GAPMul(2,q**2),0,GAPMul(2,q**2)-GAPMul(2,q**4),3-GAPMul(3,q)+q**2,2-GAPMul(4,q**2),0,0,GAPMul(2,q**6)-GAPMul(2,q**7),GAPMul(-4,q)+GAPMul(8,q**2)-GAPMul(4,q**3),0,GAPMul(3,q**6)-GAPMul(3,q**7)+q**8,0,GAPMul(2,q**24)-GAPMul(2,q**25),GAPMul(-2,q**4)+GAPMul(2,q**6),0,GAPMul(3,q**6)-GAPMul(3,q**7)+q**8,GAPMul(-6,q**6)+GAPMul(4,q**7),GAPMul(2,q**12)-GAPMul(4,q**14),GAPMul(2,q**6)-GAPMul(4,q**8),GAPMul(-6,q**12)+GAPMul(4,q**13),0,0,0,GAPMul(2,q**18)-GAPMul(2,q**19),GAPMul(10,q**6),GAPMul(10,q**30),GAPMul(10,q**12),GAPMul(10,q**24),GAPMul(10,q**18)],tbl["classtext"],lxg277)
    
    
    def f31(p,r,j):
        q=GAPDiv(-r,p)
        return Zip([15,-7+GAPMul(8,q),2-GAPMul(5,q)+GAPMul(3,q**2),GAPMul(2,q**3)-GAPMul(3,q**2)+q,1-GAPMul(4,q)+GAPMul(2,q**2),GAPMul(2,q**2)-GAPMul(2,q)+1-q**3,0,q**2+GAPMul(GAPMul(3,j**2)+j,q**4),2-GAPMul(5,q)+GAPMul(3,q**2),1-GAPMul(4,q**2)+GAPMul(2,q**4),GAPMul(j,q**4),0,GAPMul(-j**2,q**8)+GAPMul(4,j**2)-GAPMul(2,j**2),1-GAPMul(2,q)+GAPMul(GAPMul(-4,j**2)-GAPMul(2,j),q**2)-GAPMul(4,q**5)+q**6,0,GAPMul(-2,j**2)+GAPMul(5,j**2)-GAPMul(3,j**2),GAPMul(-j,q**5),GAPMul(j**2,q**32)-GAPMul(4,j**2)+GAPMul(2,j**2),GAPMul(-2,j**2)+GAPMul(3,j**2)-GAPMul(j**2,q**8),0,GAPMul(-2,j**2)+GAPMul(5,j**2)-GAPMul(3,j**2),GAPMul(7,j**2)-GAPMul(8,j**2),GAPMul(j,q**16)-GAPMul(4,j)+GAPMul(2,j),GAPMul(-j**2,q**8)+GAPMul(4,j**2)-GAPMul(2,j**2),GAPMul(-7,j)+GAPMul(8,j),0,0,0,-q**24+GAPMul(4,q**25)-GAPMul(2,q**26),GAPMul(-15,j**2),GAPMul(-15,j),GAPMul(15,j),GAPMul(15,j**2),GAPMul(-15,q**24)],tbl["classtext"],lxg278)
    
    
    tbl["irreducibles"]=GAPMul([f1(r),f1(p),f3(r,p,ER(3),ER(5)),f3(p,r,ER(3),ER(5)),f3(r,p,ER(3),-ER(5)),f3(p,r,ER(3),-ER(5)),f3(r,p,ER(3)**2,ER(5)),f3(p,r,ER(3)**2,ER(5)),f3(r,p,ER(3)**2,-ER(5)),f3(p,r,ER(3)**2,-ER(5)),f11(r,p),f12(r,p),f12(p,r),f11(p,r),f15(r,p,ER(3)),f15(p,r,ER(3)),f15(r,p,ER(3)**2),f15(p,r,ER(3)**2),f19(r,p,ER(5),1),f19(r,p,ER(5),-1),f19(r,p,-ER(5),1),f19(r,p,-ER(5),-1),f23(r,p,1),f23(p,r,1),f23(p,r,ER(3)**2),f23(r,p,ER(3)**2),f23(p,r,ER(3)),f23(r,p,ER(3)),f29(r,p),f29(p,r),f31(p,r,ER(3)),f31(r,p,ER(3)),f31(p,r,ER(3)**2),f31(r,p,ER(3)**2)],q**0)
    tbl["centralizers"]=map(lambda x: GAPDiv(tbl["order"],x),tbl["classes"])
    tbl["irredinfo"]=map(lambda x: {"charparam":x,
        "charname":ChevieData["G27"]["CharName"](x,{})},ChevieData["G27"]["CharInfo"]()["charparams"])
    tbl["galomorphisms"]=PermutationGroup([Permutation("(4,8)(6,19)(7,20)(12,28)(14,27)(15,26)"),Permutation("(6,7)(11,13)(14,15)(17,18)(19,20)(22,23)(24,25)(26,27)(30,31)(32,33)")])
    tbl=ChevieData["compat"]["MakeCharacterTable"](tbl)
    return tbl

def lxg276(x,y):
    return GAPMul(x,-p**len(y))

def lxg277(x,y):
    return GAPMul(x,-p**len(y))

def lxg278(x,y):
    return GAPMul(x,-p**len(y))

ChevieData["G27"]["HeckeCharTable"]=lxg275

def lxg279():
    return ChevieData["G27"]["HeckeCharTable"]([[1,-1],[1,-1],[1,-1]],[])

ChevieData["G27"]["CharTable"]=lxg279

def lxg2710(para,root,i):
    p=para[1-1][2-1]
    r=para[1-1][1-1]
    f1=lambda r: map(lambda x: [[r]],range(1,3+1))
    def f3(r,p,j,d):
        c=GAPDiv(2+GAPMul(j**2,1-d),2)
        return GAPMul(WGraph2Representation([[[1,2],[1,3],[2,3]],[[1,2,-2,GAPMul(p,r)],[1,3,-c,GAPDiv(GAPMul(p,r),c)],[2,3,r,-p]]],[p,r]),p**0)
    
    
    def f11(r,p):
        return GAPMul(WGraph2Representation([[[1,2],[1,3],[1,3],[2],[2,3]],[[1,2,p,0],[1,3,-p,r],[1,5,p,-r],[2,4,0,p],[2,5,p,-r],[3,4,r,-p],[3,5,p,0]]],[p,r]),p**0)
    
    
    def f12(r,p):
        return GAPMul(WGraph2Representation([[[1,3],[1,2],[1,2],[3],[2,3]],[[1,2,0,r],[1,3,p,-r],[2,4,r,0],[3,4,-r,p],[1,5,-p,r],[2,5,-p,r],[3,5,0,r]]],[p,r]),p**0)
    
    
    def f15(x,y,j):
        return GAPMul(WGraph2Representation([[[1],[1],[2],[2],[3],[3]],[[1,4,x,-y],[1,5,x,-y],[2,3,y,-x],[2,6,GAPMul(j**2,x),GAPMul(-j,y)],[3,5,GAPMul(2,y),0],[3,6,GAPMul(2,x),-y],[4,5,GAPMul(-j,y),0],[4,6,GAPMul(-j,x),0]]],[x,y]),x**0)
    
    
    def f19(x,y,d,sgn):
        v=GAPMul(-sgn,GetRoot(GAPMul(-x,y),2))
        return GAPMul([[[x,0,-y,0,0,0,v-GAPDiv(GAPMul(y,1+d),2),0],[0,x+y,0,y,0,0,0,0],[0,0,y,0,0,0,0,0],[0,-x,0,0,0,0,0,0],[0,0,0,0,0,0,0,x],[0,0,v+GAPDiv(GAPMul(x,1+d),2),0,0,x,x,0],[0,0,0,0,0,0,y,0],[0,0,0,0,-y,0,0,x+y]],[[x,0,-y,y,0,v,v,-y],[0,x,GAPDiv(GAPMul(-x,3+d),2),x,0,0,0,0],[0,0,y,0,0,0,0,0],[0,0,0,y,0,0,0,0],[0,0,x,0,x,0,0,y],[0,0,0,0,0,x+y,y,0],[0,0,0,0,0,-x,0,0],[0,0,0,0,0,0,0,y]],[[0,0,-x,0,0,0,0,0],[0,x,0,x,0,0,GAPDiv(GAPMul(-y,3+d),2),0],[y,0,x+y,0,0,0,0,0],[0,0,0,y,0,0,0,0],[0,0,0,0,y,0,0,0],[v,0,v,-y,v,x,x,0],[0,0,0,0,0,0,y,0],[0,0,0,0,-x,0,v,x]]],x**0)
    
    
    def f23(p,r,e):
        x=GAPMul(e,GetRoot(r,3))
        y=GetRoot(p,3)
        return GAPMul(WGraph2Representation([[[1],[1,2],[1,2],[1,3],[1,3],[2],[2,3],[2,3],[3]],[[1,2,0,GAPMul(x,y**2)-GAPMul(x**2,y)+x**3-y**3],[1,3,0,GAPMul(-x,y**2)+GAPMul(x**2,y)-x**3],[1,4,0,GAPMul(x,y**2)-y**3],[1,8,GAPMul(-x**3,y**3),1],[2,5,y**2,GAPMul(-x**3,y)],[2,7,GAPMul(-x**2,y**2),GAPMul(x,y)],[3,4,GAPMul(-x,y**2),GAPMul(x**2,y)],[3,5,GAPMul(-x,y)+y**2,0],[3,6,GAPMul(x,y**2)-y**3,0],[3,9,x**3,-y**3],[4,6,GAPMul(x**2,y)-x**3,0],[4,7,GAPMul(-x**2,y**2),GAPMul(x,y)],[4,9,GAPMul(-x**2,y)+x**3,0],[5,6,GAPMul(x**3,y),-y**2],[5,7,0,x-y],[5,9,GAPMul(x**2,y**2)-GAPMul(x**3,y),0],[6,7,0,GAPMul(-2,x)+x**2],[6,8,0,1],[6,9,0,GAPMul(x,y**2)-y**3],[7,9,GAPMul(2,x)-y**2,0],[8,9,-1,0]]],[x**3,y**3]),x**0)
    
    
    def f29(x,y):
        return GAPMul(WGraph2Representation([[[1,3],[1,2],[1,3],[1,2],[1,3],[1,2],[3],[2],[2,3],[2,3]],[[1,4,0,GAPMul(x,y)+x**2],[1,6,0,GAPMul(x,y)+x**2],[1,7,y,0],[1,9,-1,GAPMul(x,y)],[1,10,x,0],[2,3,0,GAPMul(2,x**2)+y**2],[2,5,0,GAPMul(x,y)+x**2],[2,9,0,GAPMul(-x,y)-y**2],[2,10,x,-y],[3,4,x,0],[3,6,x,GAPMul(-2,y)],[3,7,GAPMul(2,x**2),0],[3,8,GAPMul(2,x),0],[3,9,0,-x],[3,10,GAPMul(x,y),-1],[4,5,-y,x],[4,7,GAPMul(x,y),-1],[4,8,x+y,0],[4,9,0,-y],[5,6,0,GAPMul(3,y)],[5,7,x**2,0],[5,8,x,-y],[5,9,0,GAPMul(2,x)],[5,10,0,2],[6,7,GAPMul(x,y),0],[6,8,x-y,0],[6,9,-x,y],[6,10,GAPMul(x,y)+x**2,0],[7,8,1,GAPMul(-x,y)],[7,9,0,y**2],[8,9,0,-x+GAPMul(2,y)],[8,10,0,-1]]],[x,y]),x**0)
    
    
    def f31(r,p,j):
        return GAPMul([[[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],[GAPMul(-p,r),p+r,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,p+r,0,0,0,GAPMul(j,p),0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0],[0,0,0,0,p,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,GAPDiv(r,p),0,0,0,0,0,0],[0,0,-j**2,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,p,0,0,0,0,0,0,0],[0,0,0,0,0,-p**2,0,0,p+r,0,0,0,0,0,0],[GAPMul(-p,r**3)-GAPMul(p**2,r**2)-r**4,GAPMul(p,r**2)+GAPMul(p**2,r),GAPMul(p**2,r),GAPMul(-p,r**3)-GAPMul(p**2,r**2)-r**4,GAPMul(j,p**2),GAPMul(p,r)+p**2,GAPMul(j,p**2),r**2,GAPDiv(GAPMul(-p,r)-p**2-r**2,p),r,0,GAPMul(-p,r**2)-GAPMul(p**2,r)-r**3,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,GAPDiv(-r,p),0],[0,0,0,GAPMul(p,r),0,0,0,0,0,0,0,p+r,0,0,0],[GAPMul(-j**2-GAPMul(2,j),p)+GAPMul(GAPMul(-2,j**2)-GAPMul(3,j),p**2),GAPMul(j**2+GAPMul(2,j),p)+GAPMul(GAPMul(2,j**2)+GAPMul(3,j),p**2)-p**3,GAPMul(-j**2,p)-GAPMul(j**2,r**3),GAPMul(-j,p)+GAPMul(2,p**2)-GAPMul(j**2,p**3),GAPMul(j**2,p)+GAPMul(p**2,r**2)-r**4,GAPMul(j,p)-p**2,GAPMul(-p,r**3)-r**4,GAPMul(j,p)+GAPMul(j,p**2),p-GAPMul(j,r),0,GAPMul(-j**2,p)+r,GAPMul(-j,p)+GAPMul(2,p**2)-GAPMul(j**2,p**3),r,GAPDiv(GAPMul(-j**2,p)+r,p),0],[0,0,0,0,0,0,0,0,0,0,p**2,0,0,p+r,0],[GAPMul(p**2,r**3)+GAPMul(GAPMul(-2,j**2)-j,p**3),GAPMul(-p**2,r**2)+GAPMul(GAPMul(2,j**2)+j,p**3),GAPMul(p,r**3)+GAPMul(j**2,p**3),GAPMul(p**2,r**3)+GAPMul(j**2,p**4),GAPMul(p**2,r**3)-GAPMul(j**2,p**5)-r**5,0,GAPMul(j,p)+GAPMul(p**3,r**2),GAPMul(p,r**2)+GAPMul(j**2,p**3),0,0,GAPMul(p,r)+r**2,GAPMul(p**2,r**2)+GAPMul(j**2,p**4),0,GAPDiv(GAPMul(p,r)+r**2,p),r]],[[0,0,0,-r,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0],[0,0,p,0,0,0,0,0,0,0,0,0,0,0,0],[p,0,0,p+r,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-r,0,0,0,0,0,0,0,0],[0,GAPMul(p,r),0,0,0,p+r,0,0,0,0,0,0,0,0,0],[0,0,0,0,p,0,p+r,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0],[0,0,0,0,0,0,0,0,p+r,0,0,GAPMul(p**2,r),0,0,0],[0,0,0,0,0,0,0,GAPMul(p,r),0,p+r,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-r,0,0],[0,0,0,0,0,0,0,0,GAPDiv(-1,p),0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,p,0,p+r,0,0],[GAPMul(j**2,p)+GAPMul(j**2,p**2)-GAPMul(p**4,r),GAPMul(-j**2,p)-GAPMul(j**2,p**2),0,GAPMul(j**2,p)+GAPMul(j**2,p**2)-GAPMul(p**4,r),GAPMul(j**2,p)-GAPMul(j,p**2)-GAPMul(p**4,r),GAPMul(-j**2,p)-GAPMul(j**2,p**2),GAPMul(j**2,p)-GAPMul(j,p**2)-GAPMul(p**4,r),GAPMul(-j**2,p)-GAPMul(j**2,p**2),GAPMul(-j**2,p**2)-r**2,GAPMul(-j**2,p)-GAPMul(j**2,p**2),GAPMul(-j**2,p),GAPMul(-p,r**3)-GAPMul(j**2,p**3),GAPMul(-j**2,p),r,0],[GAPMul(p**2,r**3)-GAPMul(j,p**3)+GAPMul(j**2,p**4)-r**5,r**4,0,GAPMul(p**2,r**3)-GAPMul(j,p**3)+GAPMul(j**2,p**4)-r**5,GAPMul(p**2,r**3)-GAPMul(j,p**3)+GAPMul(j**2,p**4)-r**5,r**3,GAPMul(p**2,r**3)-GAPMul(j,p**3)+GAPMul(j**2,p**4)-r**5,r**3,0,r**2,r**2,0,r**2,0,r]],[[p,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0],[-r**2,r,r,GAPMul(p,r),-p**2,0,GAPMul(j,r**2),1,0,0,0,p,0,0,0],[0,0,0,0,-r,0,0,0,0,0,0,0,0,0,0],[0,0,0,p,p+r,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0],[0,0,0,0,0,0,p,0,0,0,0,0,0,0,0],[0,GAPMul(p,r),0,0,0,0,0,p+r,0,0,0,0,0,0,0],[GAPMul(-p,r**3)+GAPMul(p**3,r)-r**4,GAPMul(p,r**2)+r**3,0,GAPMul(-j**2,p**2)-r**4,GAPMul(-p,r**3)+GAPMul(j**2+GAPMul(2,j),p**2)-GAPMul(p**3,r)-r**4,r**2,GAPMul(j,p)+GAPMul(j,p**2),GAPMul(p,r)+r**2,r,r,r,GAPMul(p,r**2)+GAPMul(p**2,r),0,0,1],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1],[0,0,0,0,0,GAPMul(p,r),0,0,0,0,p+r,0,0,0,0],[0,0,0,GAPMul(p,r)+p**2,GAPMul(p,r)+r**2,0,0,0,0,0,0,p,0,0,0],[GAPMul(p**2,r**2)-GAPMul(j,p**3)+GAPMul(j**2,p**4)-r**4,r**3,0,GAPMul(p**2,r**2)-GAPMul(j,p**3)+GAPMul(j**2,p**4)-r**4,GAPMul(p**2,r**2)-GAPMul(j,p**3)+GAPMul(j**2,p**4)-r**4,r**2,GAPMul(p**2,r**2)-GAPMul(j,p**3)+GAPMul(j**2,p**4)-r**4,r**2,0,r,r,0,r,0,1],[GAPMul(-2,p)-GAPMul(p**2,r**3)+GAPMul(GAPMul(-2,j**2)-j,p**3)-r**5,GAPMul(2,p)+GAPMul(p**2,r**2),0,GAPMul(-p,r**4)+GAPMul(p**2,r**3)-r**5,GAPMul(-p,r**4)+GAPMul(j,p**2)-GAPMul(p**4,r)-GAPMul(j**2,p**5)-r**5,GAPMul(p,r**2)+r**3,GAPMul(j,p)+GAPMul(j,p**2),GAPMul(2,p)+GAPMul(p**2,r),0,r**2,GAPMul(p,r)+r**2,GAPMul(-j**2-GAPMul(2,j),p**2)+GAPMul(p**3,r),0,r,r],[0,0,0,0,0,0,0,0,0,GAPMul(p,r),0,0,0,0,p+r]]],p**0)
    
    
    rep=[[f1,r],[f1,p],[f3,p,r,ER(3)**2,ER(5)],[f3,r,p,ER(3)**2,ER(5)],[f3,p,r,ER(3)**2,-ER(5)],[f3,r,p,ER(3)**2,-ER(5)],[f3,p,r,ER(3),ER(5)],[f3,r,p,ER(3),ER(5)],[f3,p,r,ER(3),-ER(5)],[f3,r,p,ER(3),-ER(5)],[f11,p,r],[f12,p,r],[f12,r,p],[f11,r,p],[f15,r,p,ER(3)],[f15,p,r,ER(3)],[f15,r,p,ER(3)**2],[f15,p,r,ER(3)**2],[f19,r,p,ER(5),1],[f19,r,p,ER(5),-1],[f19,r,p,-ER(5),1],[f19,r,p,-ER(5),-1],[f23,p,r,1],[f23,r,p,1],[f23,r,p,ER(3)**2],[f23,p,r,ER(3)**2],[f23,r,p,ER(3)],[f23,p,r,ER(3)],[f29,p,r],[f29,r,p],[f31,r,p,ER(3)],[f31,p,r,ER(3)],[f31,r,p,ER(3)**2],[f31,p,r,ER(3)**2]]
    return rep[i-1][1-1](*[rep[i-1][k-1] for k in range(2,len(rep[i-1])+1)])+GAPMul(0,prod(para[1-1]))

ChevieData["G27"]["HeckeRepresentation"]=lxg2710

ChevieData["families"]["Y6"]={"name":"Y_6",
    "explanation":"subcategory of DQ(B2).20",
    "fourierMat":GAPDiv([[-ER(5),-ER(5),GAPMul(-2,ER(5)),GAPMul(-2,ER(5)),-5,-5],[-ER(5),-ER(5),GAPMul(-2,ER(5)),GAPMul(-2,ER(5)),5,5],[GAPMul(-2,ER(5)),GAPMul(-2,ER(5)),-5+ER(5),5+ER(5),0,0],[GAPMul(-2,ER(5)),GAPMul(-2,ER(5)),5+ER(5),-5+ER(5),0,0],[-5,5,0,0,5,-5],[-5,5,0,0,-5,5]],10),
    "eigenvalues":[1,1,ER(5)**3,ER(5)**2,-1,1],
    "special":1,
    "cospecial":1}

def lxg2711():
    return {"harishChandra":[{"relativeType":{"series":"ST",
        "indices":range(1,3+1),
        "rank":3,
        "ST":27},
        "levi":[],
        "parameterExponents":[1,1,1],
        "charNumbers":range(1,34+1),
        "eigenvalue":1,
        "cuspidalName":""},{"relativeType":{"series":"ST",
        "indices":[2],
        "rank":1,
        "p":6,
        "q":1},
        "levi":[1,8],
        "parameterExponents":[[GAPDiv(5,2),5,0,GAPDiv(5,2),0,5]],
        "charNumbers":[56,37,79,55,77,39],
        "eigenvalue":ER(5)**2,
        "cuspidalName":"I_2(5)[1,3]"},{"relativeType":{"series":"ST",
        "indices":[2],
        "rank":1,
        "p":6,
        "q":1},
        "levi":[1,8],
        "parameterExponents":[[GAPDiv(5,2),5,0,GAPDiv(5,2),0,5]],
        "charNumbers":[58,38,80,57,78,40],
        "eigenvalue":ER(5)**3,
        "cuspidalName":"I_2(5)[1,2]"},{"relativeType":{"series":"ST",
        "indices":[1],
        "rank":1,
        "p":6,
        "q":1},
        "levi":[2,3],
        "parameterExponents":[[4,5,0,1,0,5]],
        "charNumbers":[47,35,76,74,75,36],
        "eigenvalue":-1,
        "cuspidalName":"B_2"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[65],
        "eigenvalue":ER(4),
        "cuspidalName":"G_{27}[i]",
        "qEigen":GAPDiv(1,2)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[63],
        "eigenvalue":ER(4),
        "cuspidalName":"G_{27}^2[i]",
        "qEigen":GAPDiv(1,2)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[64],
        "eigenvalue":-ER(4),
        "cuspidalName":"G_{27}[-i]",
        "qEigen":GAPDiv(1,2)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[66],
        "eigenvalue":-ER(4),
        "cuspidalName":"G_{27}^2[-i]",
        "qEigen":GAPDiv(1,2)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[83],
        "eigenvalue":ER(3),
        "cuspidalName":"G_{27}[\\zeta_3]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[84],
        "eigenvalue":ER(3),
        "cuspidalName":"G_{27}^2[\\zeta_3]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[85],
        "eigenvalue":ER(3),
        "cuspidalName":"G_{27}^3[\\zeta_3]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[54],
        "eigenvalue":ER(3),
        "cuspidalName":"G_{27}^4[\\zeta_3]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[67],
        "eigenvalue":ER(3)**2,
        "cuspidalName":"G_{27}[\\zeta_3^2]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[43],
        "eigenvalue":ER(3)**2,
        "cuspidalName":"G_{27}^2[\\zeta_3^2]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[44],
        "eigenvalue":ER(3)**2,
        "cuspidalName":"G_{27}^3[\\zeta_3^2]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[45],
        "eigenvalue":ER(3)**2,
        "cuspidalName":"G_{27}^4[\\zeta_3^2]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[86],
        "eigenvalue":-ER(3),
        "cuspidalName":"G_{27}[-\\zeta_3]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[46],
        "eigenvalue":-ER(3)**2,
        "cuspidalName":"G_{27}^2[-\\zeta_3]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[51],
        "eigenvalue":ER(9),
        "cuspidalName":"G_{27}[\\zeta_9]",
        "qEigen":GAPDiv(1,3)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[53],
        "eigenvalue":ER(9),
        "cuspidalName":"G_{27}^2[\\zeta_9]",
        "qEigen":GAPDiv(2,3)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[72],
        "eigenvalue":ER(9)**2,
        "cuspidalName":"G_{27}[\\zeta_9^2]",
        "qEigen":GAPDiv(2,3)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[68],
        "eigenvalue":ER(9)**2,
        "cuspidalName":"G_{27}^2[E9^2]",
        "qEigen":GAPDiv(1,3)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[50],
        "eigenvalue":ER(9)**4,
        "cuspidalName":"G_{27}[\\zeta_9^4]",
        "qEigen":GAPDiv(1,3)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[49],
        "eigenvalue":ER(9)**4,
        "cuspidalName":"G_{27}^2[\\zeta_9^4]",
        "qEigen":GAPDiv(2,3)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[70],
        "eigenvalue":ER(9)**5,
        "cuspidalName":"G_{27}[\\zeta_9^5]",
        "qEigen":GAPDiv(1,3)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[69],
        "eigenvalue":ER(9)**5,
        "cuspidalName":"G_{27}^2[\\zeta_9^5]",
        "qEigen":GAPDiv(2,3)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[52],
        "eigenvalue":ER(9)**7,
        "cuspidalName":"G_{27}[\\zeta_9^7]",
        "qEigen":GAPDiv(2,3)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[48],
        "eigenvalue":ER(9)**7,
        "cuspidalName":"G_{27}^2[\\zeta_9^7]",
        "qEigen":GAPDiv(1,3)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[73],
        "eigenvalue":ER(9)**8,
        "cuspidalName":"G_{27}[\\zeta_9^8]",
        "qEigen":GAPDiv(2,3)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[71],
        "eigenvalue":ER(9)**8,
        "cuspidalName":"G_{27}^2[\\zeta_9^8]",
        "qEigen":GAPDiv(1,3)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[41],
        "eigenvalue":ER(15),
        "cuspidalName":"G_{27}[\\zeta_{15}]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[42],
        "eigenvalue":ER(15)**4,
        "cuspidalName":"G_{27}[\\zeta_{15}^4]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[81],
        "eigenvalue":ER(15)**11,
        "cuspidalName":"G_{27}[\\zeta_{15}^{11}]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[82],
        "eigenvalue":ER(15)**14,
        "cuspidalName":"G_{27}[\\zeta_{15}^{14}]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[59],
        "eigenvalue":ER(20)**17,
        "cuspidalName":"G_{27}[\\zeta_{20}^{17}]",
        "qEigen":GAPDiv(1,2)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[61],
        "eigenvalue":ER(20)**13,
        "cuspidalName":"G_{27}[\\zeta_{20}^{13}]",
        "qEigen":GAPDiv(1,2)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0,
        "qEigen":GAPDiv(1,2)},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[60],
        "eigenvalue":ER(20)**7,
        "cuspidalName":"G_{27}[\\zeta_{20}^7]",
        "qEigen":GAPDiv(1,2)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[62],
        "eigenvalue":ER(20)**3,
        "cuspidalName":"G_{27}[\\zeta_{20}^3]",
        "qEigen":GAPDiv(1,2)}],
        "families":[Family("C1",[1]),Family(GAPMul(ChevieData["families"]["X"](3),Family("Y6")),[5,3,38,37,35,16,9,7,40,39,36,18,44,43,42,41,46,45],{"signs":[-1,1,1,1,1,1,-1,1,1,1,1,1,1,-1,1,1,1,1]}),Family("C2",[30,11,12,47]),ComplexConjugate(Family("Z9",[23,49,48,28,53,51,26,52,50],{"special":7,
        "cospecial":1})),ComplexConjugate(Family(ChevieData["families"]["X"](3),[33,31,54])),Family(GAPMul(Family("C'\"2"),ChevieData["families"]["Dihedral"](5)),[19,21,58,56,20,22,57,55,65,63,59,61,66,64,60,62]),Family(ChevieData["families"]["X"](3),[34,32,67]),Family("Z9",[24,69,68,27,73,71,25,72,70],{"cospecial":4,
        "signs":[1,1,-1,1,1,-1,1,1,-1]}),Family("C2",[29,13,14,74]),Family(ComplexConjugate(GAPMul(ChevieData["families"]["X"](3),Family("Y6"))),[6,4,77,78,75,15,10,8,79,80,76,17,84,83,81,82,86,85],{"signs":[-1,1,1,1,1,1,-1,1,1,1,1,1,1,-1,1,1,1,1]}),Family("C1",[2])],
        "a":[0,45,1,16,1,16,1,16,1,16,3,3,12,12,16,1,16,1,6,6,6,6,4,9,9,4,9,4,12,3,5,8,5,8,1,1,1,1,1,1,1,1,1,1,1,1,3,4,4,4,4,4,4,5,6,6,6,6,6,6,6,6,6,6,6,6,8,9,9,9,9,9,9,12,16,16,16,16,16,16,16,16,16,16,16,16],
        "A":[0,45,29,44,29,44,29,44,29,44,33,33,42,42,44,29,44,29,39,39,39,39,36,41,41,36,41,36,42,33,37,40,37,40,29,29,29,29,29,29,29,29,29,29,29,29,33,36,36,36,36,36,36,37,39,39,39,39,39,39,39,39,39,39,39,39,40,41,41,41,41,41,41,42,44,44,44,44,44,44,44,44,44,44,44,44]}

ChevieData["G27"]["UnipotentCharacters"]=lxg2711

def lxg2712(x,y,z):
    return GAPMul(-90,x**2)+GAPMul(180,x**2)-GAPMul(135,x**2)+GAPMul(135,y**2)-GAPMul(30,x**4)+GAPMul(45,x**4)-GAPMul(10,y**6)-GAPMul(27,z**6)

def lxg2713(x,y,z):
    return GAPMul(-2430,x**2)-GAPMul(58320,x**2)+GAPMul(3240,x**2)-GAPMul(13500,x**2)-GAPMul(17280,x**2)+GAPMul(1530,x**2)-GAPMul(1260,x**2)+GAPMul(756,x**2)-GAPMul(4374,x**2)+GAPMul(4374,y**2)-GAPMul(29160,x**4)-GAPMul(3240,x**4)-GAPMul(3240,x**4)+GAPMul(20250,x**4)-GAPMul(8460,x**4)+GAPMul(3960,x**4)-GAPMul(5832,y**5)-GAPMul(3240,x**6)-GAPMul(13500,x**6)+GAPMul(900,x**6)-GAPMul(2520,x**6)+GAPMul(760,x**6)-GAPMul(1080,y**6)+GAPMul(648,y**7)-GAPMul(10800,x**8)-GAPMul(3150,x**8)-GAPMul(900,x**8)-GAPMul(2550,x**8)+GAPMul(3375,x**8)-GAPMul(2160,y**9)+GAPMul(900,x**10)-GAPMul(300,x**10)+GAPMul(450,x**10)-GAPMul(486,y**10)-GAPMul(180,y**11)-GAPMul(50,x**12)-GAPMul(26,y**12)+GAPMul(729,z**12)

def lxg2714(x,y,z):
    return GAPMul(-60584274,x**2)+GAPMul(7110680580,x**2)-GAPMul(972537030,x**2)+GAPMul(18366600960,x**2)-GAPMul(23662312110,x**2)+GAPMul(33974432640,x**2)-GAPMul(1602852300,x**2)-GAPMul(79131761730,x**2)+GAPMul(23244704460,x**2)-GAPMul(843312060,x**2)+GAPMul(537689610,x**2)-GAPMul(12222900,x**2)+GAPMul(202639995,x**2)-GAPMul(26859060,x**2)+GAPMul(11022960,x**2)-GAPMul(23914845,x**2)+GAPMul(23914845,y**2)-GAPMul(10203667200,x**4)-GAPMul(11553527340,x**4)-GAPMul(6324147900,x**4)+GAPMul(91174608450,x**4)-GAPMul(79841333880,x**4)-GAPMul(130374287100,x**4)-GAPMul(164209757760,x**4)+GAPMul(12784633380,x**4)-GAPMul(186840221400,x**4)-GAPMul(136349849070,x**4)+GAPMul(178537275900,x**4)-GAPMul(47897530740,x**4)+GAPMul(91340084160,x**4)-GAPMul(17014048380,x**4)-GAPMul(14207274450,x**4)-GAPMul(6976672560,x**4)-GAPMul(1297518480,x**4)+GAPMul(772075800,x**4)-GAPMul(513320220,x**4)-GAPMul(303264000,x**4)+GAPMul(121971510,x**4)-GAPMul(10203667200,x**6)+GAPMul(7702351560,x**6)-GAPMul(364141404900,x**6)+GAPMul(511501617990,x**6)-GAPMul(11123417340,x**6)+GAPMul(76310262000,x**6)-GAPMul(122977974600,x**6)+GAPMul(74494157760,x**6)-GAPMul(37971611325,x**6)-GAPMul(5211188460,x**6)+GAPMul(40676316750,x**6)-GAPMul(4676283900,x**6)+GAPMul(927121590,x**6)-GAPMul(37902000,x**6)+GAPMul(16601900,x**6)-GAPMul(324179010,x**6)+GAPMul(324179010,y**6)-GAPMul(50961648960,x**8)-GAPMul(325540089450,x**8)-GAPMul(681614416800,x**8)+GAPMul(17328257100,x**8)-GAPMul(97421744790,x**8)-GAPMul(32091600600,x**8)+GAPMul(171333978300,x**8)-GAPMul(145036931400,x**8)-GAPMul(6063336000,x**8)+GAPMul(56405568240,x**8)-GAPMul(26667506475,x**8)+GAPMul(14669703450,x**8)-GAPMul(2415528000,x**8)+GAPMul(473133210,x**8)-GAPMul(178922700,x**8)-GAPMul(561700,x**8)-GAPMul(1925587890,x**8)-GAPMul(1925587890,y**8)+GAPMul(451724850,y**9)-GAPMul(18145757700,x**10)+GAPMul(36987309450,x**10)-GAPMul(56624054400,x**10)+GAPMul(216134922375,x**10)-GAPMul(284886099540,x**10)-GAPMul(539226829350,x**10)-GAPMul(561988359900,x**10)+GAPMul(348456860550,x**10)-GAPMul(369114570000,x**10)-GAPMul(135367042860,x**10)+GAPMul(139885574400,x**10)-GAPMul(40396615650,x**10)-GAPMul(27645478200,x**10)+GAPMul(69410938500,x**10)-GAPMul(20828602440,x**10)-GAPMul(14036532075,x**10)-GAPMul(3687846300,x**10)-GAPMul(4739024400,x**10)+GAPMul(56859000,x**10)-GAPMul(67904660,x**10)+GAPMul(6935895540,x**10)-GAPMul(5476499505,y**10)+GAPMul(3629151540,y**11)-GAPMul(18012569400,x**12)+GAPMul(150212454750,x**12)-GAPMul(436654889100,x**12)-GAPMul(308375820900,x**12)+GAPMul(67111740000,x**12)-GAPMul(48098059200,x**12)-GAPMul(63423243000,x**12)-GAPMul(77330627100,x**12)-GAPMul(59739638400,x**12)-GAPMul(146804427000,x**12)+GAPMul(24078124800,x**12)-GAPMul(2059479000,x**12)-GAPMul(8104579200,x**12)+GAPMul(1003240950,x**12)-GAPMul(208780500,x**12)+GAPMul(244927100,x**12)-GAPMul(421216200,x**12)-GAPMul(1799616690,y**12)+GAPMul(5662405440,y**13)-GAPMul(9714872700,x**14)+GAPMul(75171017250,x**14)-GAPMul(101928415500,x**14)+GAPMul(79496356500,x**14)-GAPMul(1354968000,x**14)+GAPMul(196238112750,x**14)-GAPMul(8402988600,x**14)+GAPMul(16518208500,x**14)-GAPMul(3895724970,y**14)+GAPMul(1565717040,y**15)-GAPMul(33570328500,x**16)+GAPMul(22659021000,x**16)-GAPMul(12712117500,x**16)-GAPMul(587628000,x**16)-GAPMul(742703250,x**16)+GAPMul(522877500,x**16)-GAPMul(693500,x**16)+GAPMul(2084484375,x**16)-GAPMul(1589391315,y**16)+GAPMul(1343955240,y**17)-GAPMul(10817570250,x**18)-GAPMul(12456180000,x**18)+GAPMul(6056926875,x**18)-GAPMul(8919679500,x**18)+GAPMul(30806831250,x**18)-GAPMul(2958846000,x**18)-GAPMul(229605000,x**18)+GAPMul(153548500,x**18)-GAPMul(1330486965,y**18)-GAPMul(424030140,y**19)+GAPMul(1953112500,x**20)-GAPMul(3811050000,x**20)-GAPMul(3666195000,x**20)-GAPMul(2054079000,x**20)-GAPMul(3189487500,x**20)+GAPMul(2598480000,x**20)-GAPMul(269568750,x**20)-GAPMul(240157500,x**20)+GAPMul(58590500,x**20)-GAPMul(40297500,x**22)+GAPMul(1609233750,x**22)-GAPMul(524306250,x**22)-GAPMul(342495000,x**22)+GAPMul(676740000,x**22)-GAPMul(13710000,x**22)+GAPMul(49442500,x**22)-GAPMul(36720540,y**22)-GAPMul(24685560,y**23)-GAPMul(85151250,x**24)+GAPMul(112050000,x**24)-GAPMul(3982500,x**24)-GAPMul(14607500,x**24)-GAPMul(64985625,x**24)-GAPMul(93184155,y**24)-GAPMul(17528238,y**25)-GAPMul(13162500,x**26)-GAPMul(1950000,x**26)-GAPMul(4875000,x**26)-GAPMul(3737500,x**26)+GAPMul(6946875,x**26)-GAPMul(37500,x**28)+GAPMul(12500,x**28)-GAPMul(18750,x**28)+GAPMul(341610,y**28)

ChevieData["G27"]["Invariants"]=[lxg2712,lxg2713,lxg2714]

def lxg2715():
    return lxg2716

def lxg2716(x,y,z):
    return [[GAPMul(6,x),GAPMul(12,y**2),GAPMul(30,z)+GAPMul(GAPDiv(13,3),x**3)],[GAPMul(12,y),GAPMul(-24,z),GAPMul(GAPDiv(34,3),x)-GAPMul(GAPDiv(227,9),x**2)-GAPMul(GAPDiv(26,3),y**3)-GAPMul(GAPDiv(50,9),x**4)-GAPMul(GAPDiv(5,9),x**6)],[GAPMul(30,z),GAPMul(GAPDiv(34,3),x)-GAPMul(GAPDiv(227,9),x**2)+GAPMul(GAPDiv(26,3),x**3)-GAPMul(GAPDiv(26,3),y**4)-GAPMul(GAPDiv(50,9),x**4)-GAPMul(GAPDiv(5,9),x**6),GAPMul(GAPDiv(403,27),x)-GAPMul(GAPDiv(55,3),y**2)-GAPMul(GAPDiv(3349,54),x**2)+GAPMul(GAPDiv(5909,324),x**3)-GAPMul(GAPDiv(671,54),x**4)+GAPMul(GAPDiv(395,162),x**5)-GAPMul(GAPDiv(5,324),x**7)-GAPMul(GAPDiv(5,162),x**9)]]

ChevieData["G27"]["BasicDerivations"]=lxg2715

def lxg2717():
    return lxg2718

def lxg2718(x,y,z):
    GAPMul(return5832,x)-GAPMul(1404,x)-GAPMul(5508,x**2)+GAPMul(11664,z**3)-GAPMul(3078,x**3)-GAPMul(198,x**4)-GAPMul(1944,y**5)+GAPMul(648,x**5)-GAPMul(3271,x**5)+GAPMul(954,x**6)-GAPMul(1094,x**7)+GAPMul(198,x**8)-GAPMul(204,x**9)+GAPMul(18,x**10)-GAPMul(20,x**11)-GAPMul(x**13,y)

ChevieData["G27"]["Discriminant"]=lxg2717
