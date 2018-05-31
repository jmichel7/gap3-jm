
def lxg261(indices,title):
    print title," ",indices[1-1],"===(3)",indices[2-1],"--(3)",indices[3-1],"\n",

ChevieData["G26"]["PrintDiagram"]=lxg261

ChevieData["G26"]["GeneratingRoots"]=[[0,1,-1],[0,0,1],GAPMul(GAPDiv(-ER(4),ER(3)),[1,1,1])]

ChevieData["G26"]["HyperplaneRepresentatives"]=[1,2]

ChevieData["G26"]["EigenvaluesGeneratingReflections"]=[GAPDiv(1,2),GAPDiv(1,3),GAPDiv(1,3)]

ChevieData["G26"]["BraidRelations"]=[[[1,2,1,2],[2,1,2,1]],[[1,3],[3,1]],[[2,3,2],[3,2,3]]]

def lxg262(arg):
    if len(arg)==1 :
        return "G26"
    type_=arg[1-1]
    TeX="TeX" in arg[2-1]
    if type_==1 :
        if TeX :
            return "G_{26}"
        else:
            return "G26"
    else:
        if TeX :
            return SPrint("G(",Format(type_),")_{26}")
        else:
            return SPrint("G(",Format(type_),")26")

ChevieData["G26"]["ReflectionName"]=lxg262

ChevieData["G26"]["Size"]=1296

ChevieData["G26"]["ReflectionDegrees"]=[6,12,18]

ChevieData["G26"]["NrConjugacyClasses"]=48

def lxg263(s):
    t=[[[]],[[1],[2]],[[1,2],[1,3],[2,3]],[range(1,3+1)]]
    return t[s+1-1]

ChevieData["G26"]["ParabolicRepresentatives"]=lxg263

ChevieData["G26"]["ClassNames"]=[".","1","212","c3c3","212c22c3","12","1212","12121212","c32","1212z","c32c32","1212zzz","12z","c","cc","z","zc","zcc","zz","zzz","zzzz","zzzzz","13","13z","13zz","133","c1223","2","21212","2323c","2z","2zz","2zzz","22","c12122","3322","23","23z","23zz","232323","232323z","232323zz","323","c121","c12","c3c3c3","323zzzz","c3"]

ChevieData["G26"]["WordsClassRepresentatives"]=map(lambda x: Replace(x,".",[],"1",[1],"2",[2],"3",[3],"c",[1,2,3],"z",[1,2,3,1,2,3,1,2,3]),ChevieData["G26"]["ClassNames"])

ChevieData["G26"]["PowerMaps"]=[None,[1,1,8,19,21,7,8,7,11,28,32,8,11,15,17,19,15,17,21,1,19,21,34,7,11,28,32,34,29,29,7,11,34,28,32,34,28,32,8,1,19,21,40,42,4,40,42,4],[1,2,2,40,2,2,1,1,20,20,1,20,40,16,19,20,21,22,1,20,1,20,2,40,2,2,40,1,1,20,20,1,20,1,20,40,40,2,40,40,2,40,43,46,43,46,43,46],None,[1,2,6,42,41,3,8,7,35,33,32,31,27,18,17,22,15,14,21,20,19,16,26,39,38,23,13,34,29,30,12,11,10,28,9,37,36,25,24,40,5,4,43,48,47,46,45,44],None,[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48],None,None,None,[1,2,6,42,41,3,8,7,35,33,32,31,27,18,17,22,15,14,21,20,19,16,26,39,38,23,13,34,29,30,12,11,10,28,9,37,36,25,24,40,5,4,43,48,47,46,45,44],None,[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48],None,None,None,[1,2,6,42,41,3,8,7,35,33,32,31,27,18,17,22,15,14,21,20,19,16,26,39,38,23,13,34,29,30,12,11,10,28,9,37,36,25,24,40,5,4,43,48,47,46,45,44]]

ChevieData["G26"]["ClassInfo"]={"classtext":ChevieData["G26"]["WordsClassRepresentatives"],
    "classnames":ChevieData["G26"]["ClassNames"],
    "classparams":ChevieData["G26"]["ClassNames"],
    "orders":[1,2,6,6,6,6,3,3,6,6,3,6,6,18,9,6,9,18,3,2,3,6,6,6,6,6,6,3,3,6,6,3,6,3,6,6,6,6,6,2,6,6,4,12,12,4,12,12],
    "classes":[1,9,36,9,9,36,12,12,12,12,12,12,36,72,72,1,72,72,1,1,1,1,36,36,36,36,36,12,24,24,12,12,12,12,12,36,36,36,36,9,9,9,54,54,54,54,54,54]}

def lxg264():
    res={"charparams":[[1,0],[1,9],[1,33],[1,21],[1,24],[1,12],[2,24],[2,15],[2,12],[2,3],[2,18],[2,9],[3,6],[3,15],[3,8,2],[3,5,2],[3,8,1],[3,5,1],[3,20],[3,17],[3,16,2],[3,13,2],[3,4],[3,1],[3,16,1],[3,13,1],[6,8,2],[6,11,2],[6,8,1],[6,11,1],[6,2],[6,5],[6,4,2],[6,7,2],[6,10],[6,13],[6,4,1],[6,7,1],[8,6,1],[8,3],[8,9,2],[8,12],[8,6,2],[8,9,1],[9,8],[9,5],[9,10],[9,7]],
        "opdam":Permutation("(39,40)"),
        "extRefl":[1,24,15,4]}
    res["b"]=map(lambda x: x[2-1],res["charparams"])
    return res

ChevieData["G26"]["CharInfo"]=lxg264

def lxg265(para,root):
    c=prod(para)**0
    res={"size":1296,
        "order":1296,
        "identifier":"G26",
        "name":"G26",
        "powermap":ChevieData["G26"]["PowerMaps"],
        "parameter":[para[k-1] for k in [1,2]],
        "dim":3,
        "irredinfo":map(lambda x: {"charparam":x,
        "charname":ChevieData["G26"]["CharName"](x,{})},ChevieData["G26"]["CharInfo"]()["charparams"])}
    res.update(ChevieData["G26"]["ClassInfo"])
    res["centralizers"]=map(lambda x: GAPDiv(res["order"],x),res["classes"])
    def f10(r,u):
        return GAPMul([1,r,GAPMul(r,u**2),GAPMul(r**2,u**6),GAPMul(r**3,u**9),GAPMul(r,u),GAPMul(r**2,u**2),GAPMul(r**4,u**4),GAPMul(r,u**4),GAPMul(r**5,u**8),GAPMul(r**2,u**8),GAPMul(r**11,u**20),GAPMul(r**4,u**7),GAPMul(r,u**2),GAPMul(r**2,u**4),GAPMul(r**3,u**6),GAPMul(r**4,u**8),GAPMul(r**5,u**10),GAPMul(r**6,u**12),GAPMul(r**9,u**18),GAPMul(r**12,u**24),GAPMul(r**15,u**30),GAPMul(r,u),GAPMul(r**4,u**7),GAPMul(r**7,u**13),GAPMul(r,u**2),GAPMul(r**2,u**5),u,GAPMul(r**2,u**3),GAPMul(r,u**6),GAPMul(r**3,u**7),GAPMul(r**6,u**13),GAPMul(r**9,u**19),u**2,GAPMul(r**3,u**5),u**4,u**2,GAPMul(r**3,u**8),GAPMul(r**6,u**14),u**6,GAPMul(r**3,u**12),GAPMul(r**6,u**18),u**3,GAPMul(r**3,u**3),GAPMul(r**2,u**3),GAPMul(r**3,u**9),GAPMul(r**12,u**27),GAPMul(r,u**3)],c)
    
    
    def f23(r,p,u,v,w):
        return GAPMul([2,GAPMul(2,r),GAPMul(r,u**2+v**2),GAPMul(-2,r**2),GAPMul(r**3,u**3),GAPMul(r,u+v),GAPMul(r**2,u**2+v**2),GAPMul(r**4,u**4+v**4),GAPMul(-r,u),GAPMul(-u**3,v**3),GAPMul(r**2,u**2),GAPMul(-u**9,v**9),GAPMul(-u**3,r**4),GAPMul(r,u),GAPMul(-r**2,u**2),GAPMul(-2,r**3),GAPMul(-r**4,u**4),GAPMul(r**5,u**5),GAPMul(2,r**6),GAPMul(-2,r**9),GAPMul(2,r**12),GAPMul(-2,r**15),GAPMul(r,u+v),GAPMul(-u**3,r**4),GAPMul(u**6,r**7),GAPMul(r,u**2+v**2),GAPMul(-u**2,r**2),u+v,GAPMul(r**2,u+v),GAPMul(-2,r),GAPMul(-u**3,v**3),GAPMul(u**6,r**6),GAPMul(-u**9,r**9),u**2+v**2,GAPMul(-u**2,r**3),GAPMul(-u,v),GAPMul(u,v),GAPMul(-r**3,u**4),GAPMul(r**6,u**7),GAPMul(-2,u**3),GAPMul(2,r**3),GAPMul(-2,r**6),0,0,0,0,0,0],c)
    
    
    def f36(r,p,u,v,w):
        return GAPMul([3,GAPMul(3,r),GAPMul(r,u**2+v**2),GAPMul(3,r**2),GAPMul(r**3,u),GAPMul(r,u+v),GAPMul(r**2,u**2+v**2),GAPMul(r**4,u**4+v**4),GAPMul(r,GAPMul(u**2,v**2)+GAPMul(u**2,w**2)),GAPMul(u**2,r**5),GAPMul(r**2,GAPMul(u**4,v**4)+GAPMul(u**4,w**4)),GAPMul(u**6,v**6),GAPMul(u**2,r**4),0,0,GAPMul(3,r**3),0,0,GAPMul(3,r**6),GAPMul(3,r**9),GAPMul(3,r**12),GAPMul(3,r**15),GAPMul(r,u+v),GAPMul(u**2,r**4),GAPMul(u**4,r**7),GAPMul(r,u**2+v**2),GAPMul(u,v),u+v,GAPMul(r**2,u**3+v**3),GAPMul(3,r),GAPMul(u**2,r**3),GAPMul(u**4,r**6),GAPMul(u**6,r**9),u**2+v**2,GAPMul(u,v),GAPMul(u**2,v**2)+GAPMul(u**2,w**2),0,0,0,GAPMul(3,u**2),GAPMul(3,r**3),GAPMul(3,r**6),GAPMul(-u,v),GAPMul(-r**3,u),GAPMul(-r**2,u),GAPMul(-r**3,u**3),GAPMul(-r**12,u**9),GAPMul(-r,u)],c)
    
    
    def f31(r,p,u,v,w):
        return GAPMul([3,p+GAPMul(2,r),GAPMul(u,GAPMul(-p,v)-GAPMul(r,v)+GAPMul(r,u)),GAPMul(u**4,v**2),GAPMul(p**2,r),GAPMul(u,r),GAPMul(u,r),GAPMul(r**2,u**2),GAPMul(-u**2,v),GAPMul(u**5,r**3),GAPMul(u**4,v**2),GAPMul(u**13,r**7),GAPMul(p,r**3),0,0,GAPMul(3,p),0,0,GAPMul(3,p**2),GAPMul(3,p**3),GAPMul(3,p**4),GAPMul(3,p**5),GAPMul(p,u)+GAPMul(r,u),GAPMul(u**4,r**2),GAPMul(u**8,r**4),GAPMul(p,u**2)+GAPMul(r,u**2),GAPMul(u**3,r),GAPMul(2,u)+v,GAPMul(r,u),GAPMul(-u**3,v),GAPMul(u**4,r**2),GAPMul(u**8,r**4),GAPMul(u**12,p**3),GAPMul(2,u**2)+v**2,GAPMul(u**3,p),GAPMul(u,-v**3-GAPMul(u**2,v)+u**3),GAPMul(u,u+v),GAPMul(u**5,r**2),GAPMul(u**9,r**4),GAPMul(u**3,GAPMul(-2,v**3)+u**3),GAPMul(u**7,v**2),GAPMul(u**11,v**4),u**3,GAPMul(p,r**2),GAPMul(-p,r),GAPMul(-r**3,u**6),GAPMul(p**4,r**8),GAPMul(-u**2,r)],c)
    
    
    def f62(r,p,u,v,w):
        return GAPMul([6,GAPMul(2,p)+GAPMul(4,r),GAPMul(-p,u)-GAPMul(p,u)-GAPMul(r,u)-GAPMul(r,u)+GAPMul(r,u**2),GAPMul(-2,u**3),GAPMul(u**3,p**2),GAPMul(r,u+v),GAPMul(r,GAPMul(-2,p)-GAPMul(2,p)+GAPMul(r,u**2)),GAPMul(r**2,GAPMul(2,p**2)+GAPMul(2,p**2)),GAPMul(-u,GAPMul(p,v**2)-GAPMul(2,r)-GAPMul(2,u)+GAPMul(u**2,p)),GAPMul(-r**3,u**3),GAPMul(u**2,GAPMul(p**2,v**4)+GAPMul(2,r**2)),GAPMul(-u**9,v**6),GAPMul(-u**3,p),0,0,GAPMul(-6,p),0,0,GAPMul(6,p**2),GAPMul(-6,p**3),GAPMul(6,p**4),GAPMul(-6,p**5),GAPMul(p,u)+GAPMul(p,v),GAPMul(-u**3,p),GAPMul(u**6,r**4),GAPMul(p,u**2)+GAPMul(p,v**2),GAPMul(-u**2,r),GAPMul(3,u)+GAPMul(2,v),GAPMul(r,GAPMul(-u,p)-GAPMul(p,u)-GAPMul(p,u**2)-GAPMul(u**2,p)+GAPMul(r,u**3)),GAPMul(u**2,v),GAPMul(-u**3,p),GAPMul(u**6,r**4),GAPMul(-u**9,p**3),GAPMul(3,u**2)+GAPMul(2,v**2),GAPMul(-u**2,p),GAPMul(u**2+v**2,u**2-GAPMul(u,v)+w**2),GAPMul(u,u+v),GAPMul(-u**4,p),GAPMul(u**7,r**4),GAPMul(u**2,GAPMul(3,v**2)-GAPMul(2,u)+u**4),GAPMul(-p,r**2),GAPMul(p**2,r**4),GAPMul(u,GAPMul(-v,w)+u**2),GAPMul(v,r**2),0,0,GAPMul(u**13,v**8),0],c)
    
    
    def f83(r,p,u,v,w,eps):
        s=GAPMul(eps,GetRoot(GAPMul(-r,p),2))
        return GAPMul([8,GAPMul(4,p)+GAPMul(4,r),GAPMul(p+r,GAPMul(-u,v)-GAPMul(u,w)-GAPMul(v,w)+u**2),GAPMul(-4,u**3),GAPMul(-u**3,v**3),GAPMul(u,p+r),GAPMul(-2,p)-GAPMul(2,p)-GAPMul(2,p)+GAPMul(p**2,u**2),GAPMul(2,p**2)+GAPMul(2,p**2),GAPDiv(GAPDiv(GAPMul(s,u),p),r),GAPMul(-u**3,v),GAPDiv(GAPDiv(GAPMul(-u**2,v),p),r),GAPMul(u**9,v**4),GAPMul(-u**4,v),GAPMul(-s,u),GAPMul(p,r),GAPMul(-8,p),GAPMul(-p**2,r**2),GAPMul(-p**2,r**2),GAPMul(-8,p**3),GAPMul(8,p**4),GAPMul(8,p**6),GAPMul(-8,p**7),GAPMul(p+r,w+v),GAPMul(-u**3,r),GAPMul(-u**6,v**3),GAPMul(w**2+v**2,p+r),GAPMul(-u**2,p+r),GAPMul(4,u)+GAPMul(2,v),GAPMul(-p,r)-GAPMul(p,r)-GAPMul(p,r)-GAPMul(p,r)-GAPMul(p,r)-GAPMul(p,r)+GAPMul(p**2,u**3),GAPDiv(GAPDiv(GAPMul(u**2,s),p),r),GAPMul(-2,u**3),GAPMul(-2,u**6),GAPMul(2,u**9),GAPMul(4,u**2)+GAPMul(2,v**2),GAPMul(-2,u**2),GAPMul(-u,v**3)-GAPMul(u,w**3)+GAPMul(u**2,v**2)-GAPMul(u**3,v)-GAPMul(u**3,w)+u**4,GAPMul(u,u+v),GAPMul(-u**4,p),GAPMul(-u**7,p**3),GAPMul(u**2,GAPMul(3,v**2)-GAPMul(2,u)-GAPMul(2,u)+u**4),GAPMul(-p,r),GAPMul(-p**3,r**3),GAPMul(u,GAPMul(-v,w)+u**2),GAPMul(p,r),0,0,GAPMul(u**13,v**6),0],c)
    
    
    def f97(r,p,u,v,w,j):
        return GAPMul([9,GAPMul(3,p)+GAPMul(6,r),GAPMul(-p,u)-GAPMul(p,u)-GAPMul(v,p)-GAPMul(r,u)-GAPMul(r,u)-GAPMul(v,r)+GAPMul(r,u**2),GAPMul(3,r),GAPMul(r,u),GAPMul(r,u+v),GAPMul(r,GAPMul(-2,p)-GAPMul(2,p)-GAPMul(2,v)+GAPMul(r,u**2)),GAPMul(r**2,GAPMul(2,p**2)+GAPMul(2,p**2)),GAPMul(j**2,GAPMul(-2,r)-GAPMul(2,r)-GAPMul(2,r)+GAPMul(p,u**2)),GAPMul(u**2,v**2),GAPMul(j,GAPMul(2,r**2)+GAPMul(2,r**2)),GAPMul(u**6,v**6),GAPMul(u**2,j**2),0,0,GAPMul(9,j**2),0,0,GAPMul(9,j),GAPMul(9,p**3),GAPMul(9,j**2),GAPMul(9,j),GAPMul(p+GAPMul(2,r),u+v),GAPMul(u**2,r**2),GAPMul(u**4,r**4),GAPMul(u**2+v**2,p+GAPMul(2,r)),GAPMul(u,GAPMul(2,p)+r),GAPMul(3,u)+GAPMul(3,v),GAPMul(r,GAPMul(-u,p)-GAPMul(p,u)-GAPMul(p,v)-GAPMul(p,u**2)-GAPMul(u**2,p)-GAPMul(p,v**2)+GAPMul(r,u**3)),GAPMul(-u,v),GAPMul(3,u**2),GAPMul(3,u**4),GAPMul(3,u**6),GAPMul(3,u**2)+GAPMul(3,v**2),GAPMul(3,u),GAPMul(-u,v**3)-GAPMul(u,w**3)-GAPMul(v,w**3)+GAPMul(u**2,v**2)-GAPMul(u**3,v)-GAPMul(u**3,w)-GAPMul(v**3,w),GAPMul(v,w)+GAPMul(u,v),GAPMul(u**2,j**2),GAPMul(u**4,j),GAPMul(3,u**2)-GAPMul(2,u**3)-GAPMul(2,w**3)-GAPMul(2,w**3),GAPMul(-j**2,p),GAPMul(-j,p**2),GAPMul(-u,v),GAPMul(-j**2,p),GAPMul(j,p),GAPMul(r**3,u**3),GAPMul(-j**2,p**4),GAPMul(j,r)],c)
    
    
    r=para[1-1][1-1]
    p=para[1-1][2-1]
    u=para[2-1][1-1]
    v=para[2-1][2-1]
    w=para[2-1][3-1]
    res["irreducibles"]=[f10(r,u),f10(p,u),f10(p,w),f10(p,v),f10(r,w),f10(r,v),f23(p,r,v,w,u),f23(r,p,v,w,u),f23(p,r,u,v,w),f23(r,p,u,v,w),f23(p,r,u,w,v),f23(r,p,u,w,v),f36(r,p,u,v,w),f36(p,r,u,v,w),f31(p,r,v,u,w),f31(r,p,v,u,w),f31(p,r,u,w,v),f31(r,p,u,w,v),f31(p,r,w,v,u),f31(r,p,w,v,u),f31(p,r,w,u,v),f31(r,p,w,u,v),f31(p,r,u,v,w),f31(r,p,u,v,w),f31(p,r,v,w,u),f31(r,p,v,w,u),f62(r,p,w,u,v),f62(p,r,w,u,v),f62(r,p,v,w,u),f62(p,r,v,w,u),f62(r,p,u,v,w),f62(p,r,u,v,w),f62(r,p,v,u,w),f62(p,r,v,u,w),f62(r,p,w,v,u),f62(p,r,w,v,u),f62(r,p,u,w,v),f62(p,r,u,w,v),f83(r,p,u,v,w,1),f83(r,p,u,v,w,-1),f83(r,p,w,v,u,-1),f83(r,p,w,v,u,1),f83(r,p,v,u,w,1),f83(r,p,v,u,w,-1),f97(p,r,u,v,w,ER(3)**2),f97(r,p,u,v,w,ER(3)**2),f97(p,r,u,v,w,ER(3)),f97(r,p,u,v,w,ER(3))]
    return ChevieData["compat"]["MakeCharacterTable"](res)

ChevieData["G26"]["HeckeCharTable"]=lxg265

def lxg266():
    return ChevieData["G26"]["HeckeCharTable"]([[1,-1],[1,ER(3),ER(3)**2],[1,ER(3),ER(3)**2]],[])

ChevieData["G26"]["CharTable"]=lxg266

ChevieData["G26"]["sparseFakeDegrees"]=[[1,0],[1,9],[1,33],[1,21],[1,24],[1,12],[1,24,1,30],[1,15,1,21],[1,12,1,18],[1,3,1,9],[1,18,1,24],[1,9,1,15],[1,6,1,12,1,18],[1,15,1,21,1,27],[1,8,1,14,1,20],[1,5,1,11,1,17],[1,8,1,14,1,20],[1,5,1,11,1,17],[1,20,1,26,1,32],[1,17,1,23,1,29],[1,16,1,22,1,28],[1,13,1,19,1,25],[1,4,1,10,1,16],[1,1,1,7,1,13],[1,16,1,22,1,28],[1,13,1,19,1,25],[1,8,2,14,2,20,1,26],[1,11,2,17,2,23,1,29],[1,8,2,14,2,20,1,26],[1,11,2,17,2,23,1,29],[1,2,2,8,2,14,1,20],[1,5,2,11,2,17,1,23],[1,4,2,10,2,16,1,22],[1,7,2,13,2,19,1,25],[1,10,2,16,2,22,1,28],[1,13,2,19,2,25,1,31],[1,4,2,10,2,16,1,22],[1,7,2,13,2,19,1,25],[2,6,3,12,2,18,1,24],[1,3,2,9,3,15,2,21],[1,9,2,15,3,21,2,27],[2,12,3,18,2,24,1,30],[1,6,2,12,3,18,2,24],[2,9,3,15,2,21,1,27],[1,8,3,14,3,20,2,26],[1,5,3,11,3,17,2,23],[2,10,3,16,3,22,1,28],[2,7,3,13,3,19,1,25]]

ChevieData["G26"]["SchurModels"]={"f1_0":{"coeff":-1,
    "vcyc":[[[1,-1,0,0,0],1],[[0,0,1,-1,0],1],[[0,0,1,0,-1],1],[[1,-1,1,-1,0],2],[[1,-1,1,0,-1],2],[[1,-1,2,-2,0],1],[[1,-1,2,0,-2],1],[[1,-1,3,-2,-1],2],[[1,-1,3,-1,-2],2],[[1,-1,2,-1,-1],6],[[0,0,2,-1,-1],2],[[0,0,1,-1,0],6],[[0,0,1,0,-1],6]]},
    "f2_3":{"factor":[0,0,-1,1,0],
    "vcyc":[[[1,-1,0,0,0],1],[[0,0,1,0,-1],1],[[0,0,0,1,-1],1],[[1,-1,1,0,-1],1],[[1,-1,0,1,-1],1],[[1,-1,1,0,-1],2],[[1,-1,0,1,-1],2],[[1,-1,1,-1,0],2],[[1,-1,-1,1,0],2],[[1,-1,1,1,-2],6],[[0,0,1,1,-2],2],[[0,0,1,-1,0],6]]},
    "f3_1":{"coeff":-1,
    "vcyc":[[[-1,1,0,0,0],1],[[0,0,1,-1,0],1],[[0,0,1,0,-1],1],[[0,0,1,0,-1],2],[[0,0,0,1,-1],1],[[0,0,1,1,-2],2],[[0,0,2,-1,-1],2],[[0,0,1,0,-1],6],[[1,-1,1,0,-1],2],[[1,-1,-1,1,0],2],[[1,-1,2,-2,0],1],[[1,-1,2,1,-3],2]]},
    "f3_6":{"coeff":-1,
    "vcyc":[[[1,-1,0,0,0],1],[[1,-1,0,0,0],3],[[1,-1,1,-1,0],2],[[1,-1,1,0,-1],2],[[1,-1,-1,1,0],2],[[1,-1,0,1,-1],2],[[1,-1,-1,0,1],2],[[1,-1,0,-1,1],2],[[0,0,1,1,-2],2],[[0,0,1,-2,1],2],[[0,0,-2,1,1],2]]},
    "f6_2":{"vcyc":[[[1,-1,0,0,0],1],[[0,0,-1,1,0],1],[[0,0,1,0,-1],1],[[0,0,0,1,-1],1],[[0,0,-1,0,1],2],[[0,0,1,0,-1],6],[[0,0,1,-2,1],2],[[1,-1,0,1,-1],1],[[-1,1,1,0,-1],2],[[1,-1,0,1,-1],2],[[1,-1,3,-2,-1],2]]},
    "f8_3":{"coeff":2,
    "root":GAPDiv([1,1,0,1,1],2),
    "rootCoeff":-1,
    "vcyc":[[[0,0,1,-1,0],1],[[0,0,1,0,-1],1],[[-1,1,0,-1,1],2],[[-1,1,0,1,-1],2],[[0,-1,1,0,-2,1],2],[[0,-1,1,-2,0,1],2],[[0,-1,-1,-1,1,1],1],[[0,-1,-1,1,-1,1],1],[[0,-1,1,-1,-1,1],3],[[-1,0,1,-1,-1,1],3]]},
    "f9_7":{"rootUnity":ER(3),
    "vcyc":[[[0,0,0,0,0,2],1],[[0,0,1,-1,0],6],[[0,0,-1,0,1],6],[[0,0,0,1,-1],6],[[1,-1,-2,1,1,1],2],[[1,-1,1,-2,1,1],2],[[1,-1,1,1,-2,1],2],[[-1,1,0,0,0],1],[[1,-1,0,0,0,1],1]]}}

ChevieData["G26"]["SchurData"]=[{"name":"f1_0",
    "order":[1,2,3,4,5]},{"name":"f1_0",
    "order":[2,1,3,4,5]},{"name":"f1_0",
    "order":[2,1,5,4,3]},{"name":"f1_0",
    "order":[2,1,4,3,5]},{"name":"f1_0",
    "order":[1,2,5,4,3]},{"name":"f1_0",
    "order":[1,2,4,3,5]},{"name":"f2_3",
    "order":[2,1,4,5,3]},{"name":"f2_3",
    "order":[1,2,4,5,3]},{"name":"f2_3",
    "order":[2,1,3,4,5]},{"name":"f2_3",
    "order":[1,2,3,4,5]},{"name":"f2_3",
    "order":[2,1,3,5,4]},{"name":"f2_3",
    "order":[1,2,3,5,4]},{"name":"f3_6",
    "order":[1,2,3,4,5]},{"name":"f3_6",
    "order":[2,1,3,4,5]},{"name":"f3_1",
    "order":[2,1,4,3,5]},{"name":"f3_1",
    "order":[1,2,4,3,5]},{"name":"f3_1",
    "order":[2,1,3,5,4]},{"name":"f3_1",
    "order":[1,2,3,5,4]},{"name":"f3_1",
    "order":[2,1,5,4,3]},{"name":"f3_1",
    "order":[1,2,5,4,3]},{"name":"f3_1",
    "order":[2,1,5,3,4]},{"name":"f3_1",
    "order":[1,2,5,3,4]},{"name":"f3_1",
    "order":[2,1,3,4,5]},{"name":"f3_1",
    "order":[1,2,3,4,5]},{"name":"f3_1",
    "order":[2,1,4,5,3]},{"name":"f3_1",
    "order":[1,2,4,5,3]},{"name":"f6_2",
    "order":[1,2,5,3,4]},{"name":"f6_2",
    "order":[2,1,5,3,4]},{"name":"f6_2",
    "order":[1,2,4,5,3]},{"name":"f6_2",
    "order":[2,1,4,5,3]},{"name":"f6_2",
    "order":[1,2,3,4,5]},{"name":"f6_2",
    "order":[2,1,3,4,5]},{"name":"f6_2",
    "order":[1,2,4,3,5]},{"name":"f6_2",
    "order":[2,1,4,3,5]},{"name":"f6_2",
    "order":[1,2,5,4,3]},{"name":"f6_2",
    "order":[2,1,5,4,3]},{"name":"f6_2",
    "order":[1,2,3,5,4]},{"name":"f6_2",
    "order":[2,1,3,5,4]},{"name":"f8_3",
    "order":[1,2,3,4,5],
    "rootPower":-1},{"name":"f8_3",
    "order":[1,2,3,4,5],
    "rootPower":1},{"name":"f8_3",
    "order":[1,2,5,4,3],
    "rootPower":1},{"name":"f8_3",
    "order":[1,2,5,4,3],
    "rootPower":-1},{"name":"f8_3",
    "order":[1,2,4,3,5],
    "rootPower":-1},{"name":"f8_3",
    "order":[1,2,4,3,5],
    "rootPower":1},{"name":"f9_7",
    "order":[2,1,3,4,5],
    "rootUnityPower":2},{"name":"f9_7",
    "order":[1,2,3,4,5],
    "rootUnityPower":2},{"name":"f9_7",
    "order":[2,1,3,4,5],
    "rootUnityPower":1},{"name":"f9_7",
    "order":[1,2,3,4,5],
    "rootUnityPower":1}]

def lxg267(para,root,i):
    x=para[1-1][1-1]
    y=para[1-1][2-1]
    u=para[2-1][1-1]
    v=para[2-1][2-1]
    w=para[2-1][3-1]
    def f10(x,u):
        return [[[x]],[[u]],[[u]]]
    
    
    def f23(x,u,v):
        return [[[x,0],[0,x]],[[u,0],[-u,v]],[[v,v],[0,u]]]
    
    
    f36=lambda x: [[[x,0,0],[0,x,0],[0,0,x]],[[w,0,0],[GAPMul(u,w)+v**2,v,0],[v,1,u]],[[u,-1,v],[0,v,GAPMul(-u,w)-v**2],[0,0,w]]]
    def f31(x,y,u,v):
        return [[[y,0,0],[0,y,GAPMul(-u,y)-GAPMul(v,x)],[0,0,x]],[[v,v,0],[0,u,0],[0,1,v]],[[u,0,0],[-u,v,0],[0,0,v]]]
    
    
    def f6(r,p,u,v,w):
        return [[[r,0,0,0,0,0],[0,r,0,0,0,GAPMul(r,u)+GAPMul(p,w**2)],[0,0,r,0,GAPMul(p,w)+GAPMul(r,u),GAPMul(-p,w)-GAPMul(r,u)],[0,0,0,r,GAPMul(p,w**2)-GAPMul(r,v**2),GAPMul(r,v)+GAPMul(r,v**2)],[0,0,0,0,p,0],[0,0,0,0,0,p]],[[w,-1,v,0,0,0],[0,v,GAPMul(-u,w)-v**2,0,0,0],[0,0,u,0,0,0],[0,0,GAPMul(-v,w)-v**2,v,0,0],[0,0,0,-w**-1,w,0],[0,0,1,-w**-1,0,w]],[[u,0,0,0,0,0],[GAPMul(u,w)+v**2,v,0,0,0,0],[v,1,w,0,0,0],[GAPMul(v,w),0,0,w,0,GAPMul(v,w**2)],[0,0,0,0,w,v],[0,0,0,0,0,v]]]
    
    
    def f8(r,p,u,v,w,sgn):
        s=GAPMul(sgn,GetRoot(GAPMul(-p,r)))
        return [[[p,0,0,0,0,0,0,0],[0,p,0,0,0,0,0,0],[0,0,p,0,0,0,0,0],[0,0,0,p,0,0,0,0],[0,0,-p,0,r,0,0,0],[0,GAPMul(-p,u**2),GAPMul(p,u)-GAPMul(p,u),0,0,r,0,0],[0,0,p+GAPMul(p,u),-p+GAPMul(s,v**-1),0,0,r,0],[GAPMul(p,u)+GAPMul(s,u),0,p+GAPMul(s,v**-1)-GAPMul(r,u)+GAPMul(s,u),0,0,0,0,r]],[[u,0,0,0,0,0,0,0],[GAPMul(u**-1,w**2)+v,w,0,0,GAPMul(p**-2,r)-GAPMul(u**-2,v)-GAPMul(p**-1,s)-GAPMul(p**-1,s)-GAPMul(p**-1,s)-GAPMul(p**-1,s)-GAPMul(p**-1,s)+GAPMul(p**-1,r)-GAPMul(u**-1,v),GAPMul(u**-2,v)+GAPMul(p**-1,r),0,GAPMul(-u**-2,v)-GAPMul(p**-1,s)],[w,u,v,0,GAPMul(-p**-1,s)-GAPMul(p**-1,s)-GAPMul(p**-1,s)+GAPMul(p**-1,r)-GAPMul(u**-1,v),GAPMul(u**-1,v),0,GAPMul(-u**-1,v)],[0,0,0,v,0,GAPMul(-p**-1,s),GAPMul(p**-1,s)+v,-v],[0,0,0,0,u,0,0,0],[0,0,0,0,0,u,0,0],[0,0,0,0,GAPMul(p**-1,s)+GAPMul(r**-1,s)-u+w-GAPMul(u**2,v**-1),-w,w,0],[0,0,0,0,0,0,0,u]],[[v,-u,w,GAPMul(p**-1,s),0,0,0,0],[0,w,GAPMul(-u**-1,w**2)-v,GAPMul(-p**-1,s)+GAPMul(u**-1,v),0,0,0,0],[0,0,u,0,0,0,0,0],[0,0,0,u,0,0,0,0],[0,0,0,0,u,0,0,0],[0,0,0,0,0,w,u,0],[0,0,0,0,0,0,u,0],[0,0,0,u,GAPMul(p**-1,s)-GAPMul(p**-1,r)-GAPMul(p**-1,s)+GAPMul(p**-1,r)-GAPMul(p**-1,s)+v,GAPMul(p**-1,s)+w,u,v]]]
    
    
    def f9(r,p,u,v,w,j):
        return [[[p,0,0,0,0,0,0,0,0],[0,p,0,0,0,0,0,0,0],[0,0,p,0,0,0,0,0,0],[0,0,0,p,0,0,0,0,0],[GAPMul(-j**2,p),GAPMul(j**2,p),GAPMul(-j**2,p)+GAPMul(j,p),GAPMul(r,u**-1),r,GAPMul(j,w),0,1,0],[0,0,0,0,0,p,0,0,0],[0,GAPMul(-p,u)-GAPMul(j**2,p),GAPMul(-j,p)-GAPMul(r,u)+GAPMul(j**2,p),GAPMul(-r,w),0,0,r,-u,0],[0,0,0,0,0,0,0,p,0],[0,0,GAPMul(-j,p)+GAPMul(r,u),r,0,0,0,0,r]],[[w,0,0,0,0,0,0,0,0],[GAPMul(u,w)+v**2,v,0,0,v,0,GAPMul(u**-1,v),0,0],[v,1,u,0,1,0,u**-1,0,j**2],[0,0,0,u,-u,0,0,0,GAPMul(-j**2,u)-v],[0,0,0,0,w,0,0,0,0],[GAPMul(-j**2,p)-GAPMul(j**2,r)+GAPMul(j,r)-GAPMul(j,p)+GAPMul(j**2,p),0,0,0,GAPMul(-j**2,r),v,GAPMul(-j**2,r),0,0],[0,0,0,0,0,0,w,0,0],[GAPMul(-j,p)-GAPMul(j**2,p)-GAPMul(j**2,r)-GAPMul(r,v)+GAPMul(j**2,p),GAPMul(p,u)+GAPMul(j**2,p)-GAPMul(r,w),0,0,0,GAPMul(v,w),p,u,GAPMul(r,u**-1)+GAPMul(j,p)-GAPMul(p,v)],[0,0,0,0,0,0,-1,0,v]],[[u,-1,v,0,0,0,0,0,0],[0,v,GAPMul(-u,w)-v**2,-w,0,0,0,0,0],[0,0,w,0,0,0,0,0,0],[0,0,0,w,0,0,0,0,0],[0,0,0,w,u,0,-j**2-GAPMul(u**-1,v),0,GAPMul(-u**-1,v)-GAPMul(j**2,w)],[0,0,GAPMul(j**2,r),0,0,u,0,GAPMul(-u,w**-1),0],[0,0,0,0,0,0,v,0,GAPMul(v,w)],[0,0,GAPMul(j**2,p)+GAPMul(p,u)-GAPMul(j**2,p)-GAPMul(r,w**2),GAPMul(j**2,p)-GAPMul(r,u**-1)+GAPMul(p,w),0,0,0,v,0],[0,0,0,0,0,0,0,0,w]]]
    
    
    rep=[[f10,x,u],[f10,y,u],[f10,y,w],[f10,y,v],[f10,x,w],[f10,x,v],[f23,y,v,w],[f23,x,v,w],[f23,y,u,v],[f23,x,u,v],[f23,y,u,w],[f23,x,u,w],[f36,x],[f36,y],[f31,x,y,u,v],[f31,y,x,u,v],[f31,x,y,w,u],[f31,y,x,w,u],[f31,x,y,v,w],[f31,y,x,v,w],[f31,x,y,u,w],[f31,y,x,u,w],[f31,x,y,v,u],[f31,y,x,v,u],[f31,x,y,w,v],[f31,y,x,w,v],[f6,x,y,v,u,w],[f6,y,x,v,u,w],[f6,x,y,u,w,v],[f6,y,x,u,w,v],[f6,x,y,w,v,u],[f6,y,x,w,v,u],[f6,x,y,w,u,v],[f6,y,x,w,u,v],[f6,x,y,u,v,w],[f6,y,x,u,v,w],[f6,x,y,v,w,u],[f6,y,x,v,w,u],[f8,x,y,u,v,w,-1],[f8,x,y,u,v,w,1],[f8,x,y,w,v,u,1],[f8,x,y,w,v,u,-1],[f8,x,y,v,u,w,-1],[f8,x,y,v,u,w,1],[f9,x,y,u,v,w,ER(3)**2],[f9,y,x,u,v,w,ER(3)**2],[f9,x,y,u,v,w,ER(3)],[f9,y,x,u,v,w,ER(3)]]
    if rep[i]==None :
        return rep[i-1][1-1](*[rep[i-1][k-1] for k in range(2,len(rep[i-1])+1)])+GAPMul(0,prod(para))
    else:
        return false

ChevieData["G26"]["HeckeRepresentation"]=lxg267

def lxg268():
    J=ER(3)
    i3=ER(-3)
    return {"harishChandra":[{"relativeType":{"series":"ST",
        "indices":range(1,3+1),
        "rank":3,
        "ST":26},
        "levi":[],
        "parameterExponents":[1,1,1],
        "charNumbers":range(1,48+1),
        "eigenvalue":1,
        "cuspidalName":""},{"relativeType":{"series":"ST",
        "indices":[1,3,13],
        "rank":2,
        "p":6,
        "q":2},
        "levi":[2],
        "parameterExponents":[[0,2,2],3,1],
        "charNumbers":[102,68,71,66,53,70,60,67,54,103,69,72,99,59,98,65,50,49],
        "eigenvalue":J**2,
        "cuspidalName":ImprimitiveCuspidalName([[],[0,1],[0,1]])},{"relativeType":{"series":"ST",
        "indices":[1],
        "rank":1,
        "p":6,
        "q":1},
        "levi":range(2,3+1),
        "parameterExponents":[[3,4,3,0,3,4]],
        "charNumbers":[73,61,74,104,75,62],
        "eigenvalue":-1,
        "cuspidalName":"G_4"},{"relativeType":{"series":"ST",
        "indices":[3],
        "rank":1,
        "p":6,
        "q":1},
        "levi":range(1,2+1),
        "parameterExponents":[[4,3,1,1,0,1]],
        "charNumbers":[51,55,76,81,100,78],
        "eigenvalue":J,
        "cuspidalName":ImprimitiveCuspidalName([[0],[],[0,1,2]])},{"relativeType":{"series":"ST",
        "indices":[3],
        "rank":1,
        "p":6,
        "q":1},
        "levi":range(1,2+1),
        "parameterExponents":[[4,1,0,1,1,3]],
        "charNumbers":[52,79,101,80,77,56],
        "eigenvalue":J,
        "cuspidalName":ImprimitiveCuspidalName([[0],[0,1,2],[]])},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[92],
        "eigenvalue":1,
        "cuspidalName":"G_{26}[1]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[93],
        "eigenvalue":1,
        "cuspidalName":"G_{26}^2[1]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[94],
        "eigenvalue":1,
        "cuspidalName":"G_{26}^3[1]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[82],
        "eigenvalue":-1,
        "cuspidalName":"G_{26}[-1]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[83],
        "eigenvalue":-1,
        "cuspidalName":"G_{26}^2[-1]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[88],
        "eigenvalue":ER(3),
        "cuspidalName":"G_{26}[\\zeta_3]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[89],
        "eigenvalue":ER(3),
        "cuspidalName":"G_{26}^2[\\zeta_3]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[64],
        "eigenvalue":ER(3)**2,
        "cuspidalName":"G_{26}[\\zeta_{3}^2]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[84],
        "eigenvalue":ER(3)**2,
        "cuspidalName":"G_{26}^2[\\zeta_3^2]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[85],
        "eigenvalue":ER(3)**2,
        "cuspidalName":"G_{26}^3[\\zeta_3^2]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[90],
        "eigenvalue":-ER(3),
        "cuspidalName":"G_{26}[-\\zeta_3]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[91],
        "eigenvalue":-ER(3),
        "cuspidalName":"G_{26}^2[-\\zeta_3]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[63],
        "eigenvalue":-ER(3)**2,
        "cuspidalName":"G_{26}[-\\zeta_{3}^2]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[86],
        "eigenvalue":-ER(3)**2,
        "cuspidalName":"G_{26}^2[-\\zeta_3^2]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[87],
        "eigenvalue":-ER(3)**2,
        "cuspidalName":"G_{26}^3[-\\zeta_3^2]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[57],
        "eigenvalue":ER(4),
        "cuspidalName":"G_{26}[i]",
        "qEigen":GAPDiv(1,2)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[58],
        "eigenvalue":-ER(4),
        "cuspidalName":"G_{26}[-i]",
        "qEigen":GAPDiv(1,2)},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[95],
        "eigenvalue":ER(9)**8,
        "cuspidalName":"G_{26}[\\zeta_9^8]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[96],
        "eigenvalue":ER(9)**5,
        "cuspidalName":"G_{26}[\\zeta_9^5]"},{"relativeType":{"series":"A",
        "indices":[],
        "rank":0},
        "levi":range(1,3+1),
        "parameterExponents":[],
        "charNumbers":[97],
        "eigenvalue":ER(9)**2,
        "cuspidalName":"G_{26}[\\zeta_9^2]"}],
        "families":[Family("C1",[1]),Family(ChevieData["families"]["QZ"](3),[2,18,24,12,51,49,10,50,52],{"signs":[1,1,1,-1,-1,1,-1,-1,-1],
        "special":3,
        "cospecial":2}),Family(ChevieData["families"]["QZ"](3),[13,17,23,37,56,53,31,54,55],{"signs":[1,-1,-1,1,-1,-1,1,1,-1],
        "special":7,
        "cospecial":4}),Family("C'\"2",[40,39,57,58]),Family(GAPMul(Family("C2"),ChevieData["families"]["X"](3)),[33,27,59,22,16,60,48,46,64,61,62,63],{"signs":[1,1,-1,-1,-1,-1,1,1,1,-1,-1,1]}),Family(ChevieData["families"]["X"](3),[32,38,65],{"signs":[1,1,-1]}),Family({"fourierMat":GAPDiv([[-ER(-3),ER(-3),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),9,9,9,GAPMul(3-ER(-3),3),GAPMul(3+ER(-3),3),GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(3+ER(-3),3),2),ER(-3),-ER(-3),GAPMul(ER(-3),2),GAPMul(-ER(-3),2),GAPMul(-3+ER(-3),3),GAPMul(3+ER(-3),3),GAPMul(ER(-3),6),GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(-3+ER(-3),3),2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),9,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPDiv(GAPMul(-3+ER(-3),3),2),GAPMul(3+ER(-3),3),GAPMul(3-ER(-3),3),-9,-9,GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(-3+ER(-3),3),2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(-3+ER(-3),3),2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(ER(-3),2),ER(-3),ER(-3),GAPMul(ER(-3),6),GAPMul(ER(-3),6),GAPMul(ER(-3),6)],[ER(-3),-ER(-3),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),9,9,9,GAPMul(3+ER(-3),3),GAPMul(3-ER(-3),3),GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),-ER(-3),ER(-3),GAPMul(-ER(-3),2),GAPMul(ER(-3),2),GAPMul(-3-ER(-3),3),GAPMul(3-ER(-3),3),GAPMul(-ER(-3),6),GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),9,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPDiv(GAPMul(-3+ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPMul(3-ER(-3),3),GAPMul(3+ER(-3),3),-9,-9,GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-ER(-3),2),-ER(-3),-ER(-3),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6)],[GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),9,9,-9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),0,0,-9,-9,9,9,-9,-9,0,0,0,0,0,GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),9,9,9,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),0,0,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),0,GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),0,0,0],[GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),9,9,-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,0,-9,-9,9,9,-9,-9,0,0,0,0,0,GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),9,9,9,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),0,0,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),0,GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),0,0,0],[9,9,-9,-9,9,9,9,0,0,9,9,-9,-9,9,9,0,0,0,0,0,9,-9,9,-9,-9,-9,-9,-9,-9,-9,-9,0,0,9,9,-9,9,-9,9,-9,9,-9,-9,0,-9,9,0,0,0],[9,9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),9,9,9,0,0,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),9,9,0,0,0,0,0,GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),-9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,0,9,9,GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,-9,9,0,0,0],[9,9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),9,9,9,0,0,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),9,9,0,0,0,0,0,GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),0,0,9,9,GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),0,-9,9,0,0,0],[GAPMul(3-ER(-3),3),GAPMul(3+ER(-3),3),0,0,0,0,0,18,18,0,0,0,0,GAPMul(-ER(-3),6),GAPMul(ER(-3),6),GAPMul(-3-ER(-3),3),GAPMul(-3+ER(-3),3),GAPMul(-18,ER(3)**2),GAPMul(18,ER(3)),0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(18,ER(3)**2),GAPMul(18,ER(3)),0,0,0,0,0,0,0,0,0,0,GAPMul(ER(-3),6),GAPMul(3+ER(-3),3),GAPMul(-3+ER(-3),3),0,0,0],[GAPMul(3+ER(-3),3),GAPMul(3-ER(-3),3),0,0,0,0,0,18,18,0,0,0,0,GAPMul(ER(-3),6),GAPMul(-ER(-3),6),GAPMul(-3+ER(-3),3),GAPMul(-3-ER(-3),3),GAPMul(-18,ER(3)),GAPMul(18,ER(3)**2),0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(18,ER(3)),GAPMul(18,ER(3)**2),0,0,0,0,0,0,0,0,0,0,GAPMul(-ER(-3),6),GAPMul(3-ER(-3),3),GAPMul(-3-ER(-3),3),0,0,0],[GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(3+ER(-3),3),2),-9,-9,9,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),0,0,9,9,9,9,GAPMul(-ER(-3),3),GAPMul(ER(-3),3),GAPMul(3+ER(-3),3),GAPMul(3-ER(-3),3),0,0,0,GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),9,9,9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,0,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-ER(-3),6),GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(-3+ER(-3),3),2),0,0,0],[GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),-9,-9,9,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),0,0,9,9,9,9,GAPMul(ER(-3),3),GAPMul(-ER(-3),3),GAPMul(3-ER(-3),3),GAPMul(3+ER(-3),3),0,0,0,GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),9,9,9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),0,0,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(ER(-3),6),GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),0,0,0],[GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(3+ER(-3),3),2),9,9,-9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),0,0,9,9,9,9,GAPMul(-ER(-3),3),GAPMul(ER(-3),3),GAPMul(3+ER(-3),3),GAPMul(3-ER(-3),3),0,0,0,GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),-9,-9,-9,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,0,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-ER(-3),6),GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(-3+ER(-3),3),2),0,0,0],[GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),9,9,-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,0,9,9,9,9,GAPMul(ER(-3),3),GAPMul(-ER(-3),3),GAPMul(3-ER(-3),3),GAPMul(3+ER(-3),3),0,0,0,GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),-9,-9,-9,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),0,0,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(ER(-3),6),GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),0,0,0],[ER(-3),-ER(-3),-9,-9,9,9,9,GAPMul(-ER(-3),6),GAPMul(ER(-3),6),GAPMul(-ER(-3),3),GAPMul(ER(-3),3),GAPMul(-ER(-3),3),GAPMul(ER(-3),3),-ER(-3),ER(-3),GAPMul(-ER(-3),2),GAPMul(ER(-3),2),GAPMul(ER(-3),6),GAPMul(ER(-3),6),GAPMul(-ER(-3),6),GAPMul(ER(-3),3),GAPMul(ER(-3),3),9,-9,9,9,9,-9,-9,GAPMul(-ER(-3),3),GAPMul(ER(-3),3),GAPMul(ER(-3),6),GAPMul(-ER(-3),6),-9,-9,GAPMul(ER(-3),3),GAPMul(ER(-3),3),9,-9,GAPMul(ER(-3),3),GAPMul(ER(-3),3),9,9,GAPMul(-ER(-3),2),-ER(-3),-ER(-3),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6)],[-ER(-3),ER(-3),-9,-9,9,9,9,GAPMul(ER(-3),6),GAPMul(-ER(-3),6),GAPMul(ER(-3),3),GAPMul(-ER(-3),3),GAPMul(ER(-3),3),GAPMul(-ER(-3),3),ER(-3),-ER(-3),GAPMul(ER(-3),2),GAPMul(-ER(-3),2),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6),GAPMul(ER(-3),6),GAPMul(-ER(-3),3),GAPMul(-ER(-3),3),9,-9,9,9,9,-9,-9,GAPMul(ER(-3),3),GAPMul(-ER(-3),3),GAPMul(-ER(-3),6),GAPMul(ER(-3),6),-9,-9,GAPMul(-ER(-3),3),GAPMul(-ER(-3),3),9,-9,GAPMul(-ER(-3),3),GAPMul(-ER(-3),3),9,9,GAPMul(ER(-3),2),ER(-3),ER(-3),GAPMul(ER(-3),6),GAPMul(ER(-3),6),GAPMul(ER(-3),6)],[GAPMul(ER(-3),2),GAPMul(-ER(-3),2),0,0,0,0,0,GAPMul(-3-ER(-3),3),GAPMul(-3+ER(-3),3),GAPMul(3+ER(-3),3),GAPMul(3-ER(-3),3),GAPMul(3+ER(-3),3),GAPMul(3-ER(-3),3),GAPMul(-ER(-3),2),GAPMul(ER(-3),2),GAPMul(-ER(-3),4),GAPMul(ER(-3),4),GAPMul(3+ER(-3),3),GAPMul(-3+ER(-3),3),GAPMul(-ER(-3),12),GAPMul(3-ER(-3),3),GAPMul(-3-ER(-3),3),0,0,0,0,0,0,0,GAPMul(-3+ER(-3),3),GAPMul(-3-ER(-3),3),GAPMul(-3+ER(-3),3),GAPMul(-3-ER(-3),3),0,0,GAPMul(3-ER(-3),3),GAPMul(-3-ER(-3),3),0,0,GAPMul(3-ER(-3),3),GAPMul(-3-ER(-3),3),0,0,GAPMul(-ER(-3),4),GAPMul(-ER(-3),2),GAPMul(-ER(-3),2),GAPMul(ER(-3),6),GAPMul(ER(-3),6),GAPMul(ER(-3),6)],[GAPMul(-ER(-3),2),GAPMul(ER(-3),2),0,0,0,0,0,GAPMul(-3+ER(-3),3),GAPMul(-3-ER(-3),3),GAPMul(3-ER(-3),3),GAPMul(3+ER(-3),3),GAPMul(3-ER(-3),3),GAPMul(3+ER(-3),3),GAPMul(ER(-3),2),GAPMul(-ER(-3),2),GAPMul(ER(-3),4),GAPMul(-ER(-3),4),GAPMul(3-ER(-3),3),GAPMul(-3-ER(-3),3),GAPMul(ER(-3),12),GAPMul(3+ER(-3),3),GAPMul(-3+ER(-3),3),0,0,0,0,0,0,0,GAPMul(-3-ER(-3),3),GAPMul(-3+ER(-3),3),GAPMul(-3-ER(-3),3),GAPMul(-3+ER(-3),3),0,0,GAPMul(3+ER(-3),3),GAPMul(-3+ER(-3),3),0,0,GAPMul(3+ER(-3),3),GAPMul(-3+ER(-3),3),0,0,GAPMul(ER(-3),4),GAPMul(ER(-3),2),GAPMul(ER(-3),2),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6)],[GAPMul(-3+ER(-3),3),GAPMul(-3-ER(-3),3),0,0,0,0,0,GAPMul(-18,ER(3)**2),GAPMul(-18,ER(3)),0,0,0,0,GAPMul(ER(-3),6),GAPMul(-ER(-3),6),GAPMul(3+ER(-3),3),GAPMul(3-ER(-3),3),GAPMul(18,ER(3)),GAPMul(-18,ER(3)**2),0,0,0,0,0,0,0,0,0,0,0,0,-18,-18,0,0,0,0,0,0,0,0,0,0,GAPMul(-ER(-3),6),GAPMul(-3-ER(-3),3),GAPMul(3-ER(-3),3),0,0,0],[GAPMul(3+ER(-3),3),GAPMul(3-ER(-3),3),0,0,0,0,0,GAPMul(18,ER(3)),GAPMul(18,ER(3)**2),0,0,0,0,GAPMul(ER(-3),6),GAPMul(-ER(-3),6),GAPMul(-3+ER(-3),3),GAPMul(-3-ER(-3),3),GAPMul(-18,ER(3)**2),GAPMul(18,ER(3)),0,0,0,0,0,0,0,0,0,0,0,0,18,18,0,0,0,0,0,0,0,0,0,0,GAPMul(-ER(-3),6),GAPMul(3-ER(-3),3),GAPMul(-3-ER(-3),3),0,0,0],[GAPMul(ER(-3),6),GAPMul(-ER(-3),6),0,0,0,0,0,0,0,0,0,0,0,GAPMul(-ER(-3),6),GAPMul(ER(-3),6),GAPMul(-ER(-3),12),GAPMul(ER(-3),12),0,0,GAPMul(ER(-3),18),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(-ER(-3),12),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6),0,0,0],[GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),9,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),0,0,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(ER(-3),3),GAPMul(-ER(-3),3),GAPMul(3-ER(-3),3),GAPMul(3+ER(-3),3),0,0,0,GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),9,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),-9,-9,-9,-9,0,0,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),9,-9,9,9,GAPMul(ER(-3),6),GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),0,0,0],[GAPDiv(GAPMul(-3+ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),-9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),0,0,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(ER(-3),3),GAPMul(-ER(-3),3),GAPMul(-3-ER(-3),3),GAPMul(-3+ER(-3),3),0,0,0,GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),9,9,9,9,0,0,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),-9,9,-9,-9,GAPMul(ER(-3),6),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),0,0,0],[GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),9,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),0,0,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),9,9,0,0,0,0,0,GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),-9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),-9,-9,-9,-9,0,0,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),-9,9,-9,-9,0,GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),0,0,0],[GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),-9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),0,0,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),-9,-9,0,0,0,0,0,GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),9,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),9,9,9,9,0,0,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),9,-9,9,9,0,GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),0,0,0],[9,9,9,9,-9,-9,-9,0,0,9,9,-9,-9,9,9,0,0,0,0,0,9,-9,-9,9,9,9,9,9,9,-9,-9,0,0,-9,-9,-9,9,9,-9,-9,9,9,9,0,-9,9,0,0,0],[GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),9,9,-9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),0,0,9,9,-9,-9,9,9,0,0,0,0,0,GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),9,9,9,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,0,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),0,GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),0,0,0],[GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),9,9,-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,0,9,9,-9,-9,9,9,0,0,0,0,0,GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),9,9,9,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),0,0,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),0,GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),0,0,0],[GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,0,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),-9,-9,0,0,0,0,0,-9,9,-9,9,9,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),0,0,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),9,-9,9,-9,GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),0,GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),0,0,0],[GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),-9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),0,0,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),-9,-9,0,0,0,0,0,-9,9,-9,9,9,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),0,0,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),9,-9,9,-9,GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),0,GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),0,0,0],[GAPDiv(GAPMul(-3-ER(-3),3),2),GAPDiv(GAPMul(-3+ER(-3),3),2),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,0,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-ER(-3),3),GAPMul(ER(-3),3),GAPMul(-3+ER(-3),3),GAPMul(-3-ER(-3),3),0,0,0,-9,9,-9,9,-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),0,0,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),-9,9,-9,9,GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-ER(-3),6),GAPDiv(GAPMul(-3+ER(-3),3),2),GAPDiv(GAPMul(3+ER(-3),3),2),0,0,0],[GAPDiv(GAPMul(-3+ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),-9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),0,0,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(ER(-3),3),GAPMul(-ER(-3),3),GAPMul(-3-ER(-3),3),GAPMul(-3+ER(-3),3),0,0,0,-9,9,-9,9,-9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),0,0,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),-9,9,-9,9,GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(ER(-3),6),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),0,0,0],[GAPMul(3+ER(-3),3),GAPMul(3-ER(-3),3),0,0,0,0,0,GAPMul(18,ER(3)**2),GAPMul(18,ER(3)),0,0,0,0,GAPMul(ER(-3),6),GAPMul(-ER(-3),6),GAPMul(-3+ER(-3),3),GAPMul(-3-ER(-3),3),-18,18,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(18,ER(3)**2),GAPMul(18,ER(3)),0,0,0,0,0,0,0,0,0,0,GAPMul(-ER(-3),6),GAPMul(3-ER(-3),3),GAPMul(-3-ER(-3),3),0,0,0],[GAPMul(3-ER(-3),3),GAPMul(3+ER(-3),3),0,0,0,0,0,GAPMul(18,ER(3)),GAPMul(18,ER(3)**2),0,0,0,0,GAPMul(-ER(-3),6),GAPMul(ER(-3),6),GAPMul(-3-ER(-3),3),GAPMul(-3+ER(-3),3),-18,18,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(18,ER(3)),GAPMul(18,ER(3)**2),0,0,0,0,0,0,0,0,0,0,GAPMul(ER(-3),6),GAPMul(3+ER(-3),3),GAPMul(-3+ER(-3),3),0,0,0],[-9,-9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),9,9,9,0,0,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),-9,-9,0,0,0,0,0,GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),-9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),0,0,9,9,GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,9,-9,0,0,0],[-9,-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),9,9,9,0,0,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),-9,-9,0,0,0,0,0,GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),0,0,9,9,GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),0,9,-9,0,0,0],[GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,0,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(ER(-3),3),GAPMul(-ER(-3),3),GAPMul(3-ER(-3),3),GAPMul(3+ER(-3),3),0,0,0,GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),-9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),9,9,-9,-9,0,0,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),9,-9,-9,-9,GAPMul(ER(-3),6),GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),0,0,0],[GAPDiv(GAPMul(-3+ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),9,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),0,0,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(ER(-3),3),GAPMul(-ER(-3),3),GAPMul(-3-ER(-3),3),GAPMul(-3+ER(-3),3),0,0,0,GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),9,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),-9,-9,9,9,0,0,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),-9,9,9,9,GAPMul(ER(-3),6),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),0,0,0],[GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,0,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),9,9,0,0,0,0,0,GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),9,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),9,9,-9,-9,0,0,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),-9,9,9,9,0,GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),0,0,0],[GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),9,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),0,0,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),-9,-9,0,0,0,0,0,GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),-9,-9,9,9,0,0,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),9,-9,-9,-9,0,GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),0,0,0],[GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,0,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(ER(-3),3),GAPMul(-ER(-3),3),GAPMul(3-ER(-3),3),GAPMul(3+ER(-3),3),0,0,0,9,-9,-9,9,-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,0,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),9,-9,-9,9,GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(ER(-3),6),GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),0,0,0],[GAPDiv(GAPMul(-3+ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),9,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),0,0,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(ER(-3),3),GAPMul(-ER(-3),3),GAPMul(-3-ER(-3),3),GAPMul(-3+ER(-3),3),0,0,0,-9,9,9,-9,9,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),0,0,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),-9,9,9,-9,GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(ER(-3),6),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),0,0,0],[GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),-9,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,0,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),9,9,0,0,0,0,0,9,-9,-9,9,9,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),0,0,GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),-9,9,9,-9,GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),0,GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),0,0,0],[GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),-9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),0,0,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),9,9,0,0,0,0,0,9,-9,-9,9,9,GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),0,0,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),-9,9,9,-9,GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),0,GAPMul(-9,ER(3)),GAPMul(9,ER(3)**2),0,0,0],[GAPMul(ER(-3),2),GAPMul(-ER(-3),2),0,0,0,0,0,GAPMul(ER(-3),6),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6),GAPMul(ER(-3),6),GAPMul(-ER(-3),6),GAPMul(ER(-3),6),GAPMul(-ER(-3),2),GAPMul(ER(-3),2),GAPMul(-ER(-3),4),GAPMul(ER(-3),4),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6),GAPMul(-ER(-3),12),GAPMul(ER(-3),6),GAPMul(ER(-3),6),0,0,0,0,0,0,0,GAPMul(-ER(-3),6),GAPMul(ER(-3),6),GAPMul(-ER(-3),6),GAPMul(ER(-3),6),0,0,GAPMul(ER(-3),6),GAPMul(ER(-3),6),0,0,GAPMul(ER(-3),6),GAPMul(ER(-3),6),0,0,GAPMul(-ER(-3),4),GAPMul(-ER(-3),2),GAPMul(-ER(-3),2),GAPMul(ER(-3),6),GAPMul(ER(-3),6),GAPMul(ER(-3),6)],[ER(-3),-ER(-3),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),-9,-9,-9,GAPMul(3+ER(-3),3),GAPMul(3-ER(-3),3),GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),-ER(-3),ER(-3),GAPMul(-ER(-3),2),GAPMul(ER(-3),2),GAPMul(-3-ER(-3),3),GAPMul(3-ER(-3),3),GAPMul(-ER(-3),6),GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),-9,GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPDiv(GAPMul(-3+ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPMul(3-ER(-3),3),GAPMul(3+ER(-3),3),9,9,GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPMul(-9,ER(3)**2),GAPMul(9,ER(3)),GAPDiv(GAPMul(3-ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),GAPMul(-ER(-3),2),-ER(-3),-ER(-3),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6)],[ER(-3),-ER(-3),GAPMul(-9,ER(3)**2),GAPMul(-9,ER(3)),9,9,9,GAPMul(-3+ER(-3),3),GAPMul(-3-ER(-3),3),GAPDiv(GAPMul(-3+ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPDiv(GAPMul(-3+ER(-3),3),2),GAPDiv(GAPMul(-3-ER(-3),3),2),-ER(-3),ER(-3),GAPMul(-ER(-3),2),GAPMul(ER(-3),2),GAPMul(3-ER(-3),3),GAPMul(-3-ER(-3),3),GAPMul(-ER(-3),6),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),9,GAPMul(9,ER(3)**2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)),GAPMul(-9,ER(3)**2),GAPDiv(GAPMul(3+ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),GAPMul(-3-ER(-3),3),GAPMul(-3+ER(-3),3),-9,-9,GAPDiv(GAPMul(-3-ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),GAPMul(9,ER(3)),GAPMul(-9,ER(3)**2),GAPDiv(GAPMul(-3-ER(-3),3),2),GAPDiv(GAPMul(3-ER(-3),3),2),GAPMul(9,ER(3)),GAPMul(9,ER(3)**2),GAPMul(-ER(-3),2),-ER(-3),-ER(-3),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6)],[GAPMul(ER(-3),6),GAPMul(-ER(-3),6),0,0,0,0,0,0,0,0,0,0,0,GAPMul(-ER(-3),6),GAPMul(ER(-3),6),GAPMul(ER(-3),6),GAPMul(-ER(-3),6),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(-3),6),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6),GAPMul(18,ER(9)**7)-GAPMul(18,ER(9)**2),GAPMul(-18,ER(9)**5)+GAPMul(18,ER(9)**4),GAPMul(-18,ER(9)**7)+GAPMul(18,ER(9)**5)-GAPMul(18,ER(9)**4)+GAPMul(18,ER(9)**2)],[GAPMul(ER(-3),6),GAPMul(-ER(-3),6),0,0,0,0,0,0,0,0,0,0,0,GAPMul(-ER(-3),6),GAPMul(ER(-3),6),GAPMul(ER(-3),6),GAPMul(-ER(-3),6),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(-3),6),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6),GAPMul(-18,ER(9)**5)+GAPMul(18,ER(9)**4),GAPMul(-18,ER(9)**7)+GAPMul(18,ER(9)**5)-GAPMul(18,ER(9)**4)+GAPMul(18,ER(9)**2),GAPMul(18,ER(9)**7)-GAPMul(18,ER(9)**2)],[GAPMul(ER(-3),6),GAPMul(-ER(-3),6),0,0,0,0,0,0,0,0,0,0,0,GAPMul(-ER(-3),6),GAPMul(ER(-3),6),GAPMul(ER(-3),6),GAPMul(-ER(-3),6),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(-3),6),GAPMul(-ER(-3),6),GAPMul(-ER(-3),6),GAPMul(-18,ER(9)**7)+GAPMul(18,ER(9)**5)-GAPMul(18,ER(9)**4)+GAPMul(18,ER(9)**2),GAPMul(18,ER(9)**7)-GAPMul(18,ER(9)**2),GAPMul(-18,ER(9)**5)+GAPMul(18,ER(9)**4)]],54),
        "eigenvalues":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,ER(3)**2,ER(3)**2,ER(3)**2,ER(3)**2,ER(3)**2,ER(3)**2,ER(3)**2,-1,-1,-1,ER(3),ER(3),ER(3),ER(3),ER(3),ER(3),-1,-1,ER(3)**2,ER(3)**2,-ER(3)**2,-ER(3)**2,ER(3),ER(3),-ER(3),-ER(3),1,1,1,ER(9)**8,ER(9)**5,ER(9)**2],
        "explanation":"mystery G26",
        "special":1,
        "cospecial":2},[43,42,28,34,8,41,44,29,35,15,21,45,47,6,5,11,9,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97],{"signs":[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,-1,1,1,1,1,1,1,1,-1,-1,-1]}),Family(ChevieData["families"]["QZ"](3),[14,36,30,20,100,98,26,99,101],{"signs":[1,1,1,-1,-1,1,-1,-1,-1],
        "special":3,
        "cospecial":2}),Family(ChevieData["families"]["X"](3),[25,19,102],{"signs":[1,1,-1]}),Family("X5",[4,7,104,103,3],{"signs":[1,1,-1,-1,1]})],
        "a":[0,1,21,21,6,6,21,6,6,1,6,1,2,11,6,4,2,1,16,11,6,4,2,1,16,11,4,6,6,11,2,5,4,6,6,11,2,5,3,3,6,6,6,6,6,4,6,4,1,1,1,1,2,2,2,2,3,3,4,4,4,4,4,4,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,11,11,11,11,16,21,21],
        "A":[0,17,33,33,30,30,33,30,30,17,30,17,22,31,30,26,22,17,32,31,30,26,22,17,32,31,26,30,30,31,22,25,26,30,30,31,22,25,24,24,30,30,30,30,30,26,30,26,17,17,17,17,22,22,22,22,24,24,26,26,26,26,26,26,25,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,31,31,31,31,32,33,33]}

ChevieData["G26"]["UnipotentCharacters"]=lxg268

def lxg269(x1,x2,x3):
    return GAPMul(-10,x1**3)-GAPMul(10,x1**3)-GAPMul(10,x2**3)+x1**6

def lxg2610(x1,x2,x3):
    return GAPMul(2,x1**3)+GAPMul(2,x1**3)-GAPMul(4,x1**6)-GAPMul(4,x1**6)-GAPMul(4,x2**6)+GAPMul(x1**9,x2**3)

def lxg2611(x1,x2,x3):
    return GAPMul(-2,x1**3)+GAPMul(2,x1**3)-GAPMul(2,x1**3)+GAPMul(2,x1**6)-GAPMul(6,x1**6)+GAPMul(2,x1**6)-GAPMul(2,x1**9)-GAPMul(2,x1**9)-GAPMul(2,x2**9)-GAPMul(2,x1**12)+GAPMul(x1**12,x2**6)

ChevieData["G26"]["Invariants"]=[lxg269,lxg2610,lxg2611]

def lxg2612():
    return lxg2613

def lxg2613(t1,t2,t3):
    return GAPMul(36,t1)-GAPMul(t1**2,t2**2)+GAPMul(108,t3**3)-GAPMul(32,t2**3)+GAPMul(t1**3,t3**2)

ChevieData["G26"]["Discriminant"]=lxg2612
