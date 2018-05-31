
def lximp1(arg):
    p=arg[1-1]
    q=arg[2-1]
    r=arg[3-1]
    indices=arg[4-1]
    title=arg[5-1]
    print title," ",
    indent=len(title)+1
    g=lambda i: just("",indent-i)
    if q==1 :
        print indices[1-1],"(",p,")",
        if len(indices)>1 :
            print "===",
        print Join([indices[k-1] for k in range(2,len(indices)+1)],"--"),"\n",
    else:
        if p==q :
            print indices[1-1],"\n",g(0),"|\\\n",
            if p!=3 :
                print just(p,indent),
            else:
                print g(0),
            print "|==",indices[3-1],
            for j in range(4,r+1):
                print " - ",indices[j-1],
            print "\n",
            print g(0),"|/\n",g(0),indices[2-1],"\n",
        else:
            if q==2 :
                print indices[2-1],"\n",g(2),"/3|",
                if r>=3 :
                    print "\\",
                print "\n",
                if GAPDiv(p,q)>2 :
                    print g(len(just(GAPDiv(p,q)))+5),"(",GAPDiv(p,q),")",
                else:
                    print g(3),
                print indices[1-1],"  | ",
                for j in range(3,r+1):
                    print indices[j+1-1],
                    if j!=r :
                        print "-",
                print "\n",g(2),"\\ |",
                if r>=3 :
                    print "/",
                print "\n",just(indices[3-1],indent+1),"   ",IntListToString([indices[k-1] for k in [1,2,3]]),"==",IntListToString([indices[k-1] for k in [2,3,1]]),"==",IntListToString([indices[k-1] for k in [3,1,2]]),"\n",
            else:
                print indices[2-1],"\n",g(2),"/",q+1," ",
                if r>=3 :
                    print "\\",
                print "\n",
                if GAPDiv(p,q)>2 :
                    print g(len(SPrint(GAPDiv(p,q)))+5),"(",GAPDiv(p,q),")",
                else:
                    print g(3),
                print indices[1-1],"   ",
                if r>=3 :
                    print "==",
                for j in range(3,r+1):
                    print indices[j+1-1],
                    if j!=r :
                        print "-",
                print "\n",g(2),"\\  ",
                if r>=3 :
                    print "/",
                print "\n",just(indices[3-1],indent+1),
                j=ChevieData["imp"]["BraidRelations"](p,q,r)
                for g in range(1,Minimum(3,r)+1):
                    print "   ",IntListToString([indices[k-1] for k in j[g-1][1-1]]),"==",IntListToString([indices[k-1] for k in j[g-1][2-1]]),
                print "\n",

ChevieData["imp"]["PrintDiagram"]=lximp1

def lximp2(p,q,r):
    return r

ChevieData["imp"]["SemisimpleRank"]=lximp2

def lximp3(p,q,r):
    def b(i,j,o):
        def p(i,j):
            return map(lambda k: GAPMul(i,k%2)+GAPMul(j,1-k%2),range(1,o+1))
        
        
        return [p(i,j),p(j,i)]
    
    
    res=[]
    if q==1 :
        if r>=2 :
            if p==1 :
                res.append(b(1,2,3))
            else:
                res.append(b(1,2,4))
        res+=map(lambda i: b(i,i-1,3),range(3,r+1))
        for i in range(3,r+1):
            res+=map(lambda j: b(i,j,2),range(1,i-2+1))
    else:
        if p==q :
            res.append(b(1,2,p))
            if r>=3 :
                res+=[[[1,2,3,1,2,3],[3,1,2,3,1,2]],b(1,3,3),b(2,3,3)]
            res+=map(lambda i: b(i,i-1,3),range(4,r+1))
            for i in range(4,r+1):
                res+=map(lambda j: b(i,j,2),range(1,i-2+1))
        else:
            res.append([[1,2,3],[2,3,1]])
            i=b(2,3,q-1)
            res.append([Concatenation([1,2],i[2-1]),Concatenation([3,1],i[1-1])])
            if r>=3 :
                if q!=2 :
                    res.append([[2,3,4,2,3,4],[4,2,3,4,2,3]])
                res+=[b(2,4,3),b(3,4,3),b(1,4,2)]
            res+=map(lambda i: b(i,i-1,3),range(5,r+1+1))
            for i in range(5,r+1+1):
                res+=map(lambda j: b(i,j,2),range(1,i-2+1))
    return res

ChevieData["imp"]["BraidRelations"]=lximp3

def lximp4(p,q,r):
    return GAPDiv(GAPMul(p**r,Factorial(r)),q)

ChevieData["imp"]["Size"]=lximp4

def lximp5(arg):
    option=arg[len(arg)-1]
    if arg[3-1]==1 and arg[2-1]==1 :
        if "TeX" in option :
            return SPrint("Z_{",arg[1-1],"}")
        else:
            return SPrint("Z",arg[1-1])
    if "TeX" in option :
        n=SPrint("G_{",Join([arg[k-1] for k in range(1,3+1)]),"}")
    else:
        n=SPrint("G",IntListToString([arg[k-1] for k in range(1,3+1)]))
    if len(arg)==5 :
        PrintToString(n,"(",Format(arg[4-1],option),")")
    return n

ChevieData["imp"]["ReflectionName"]=lximp5

def lximp6(p,q,r):
    if q==1 :
        v=GAPMul(range(1,r+1),0)
        v[1-1]=1
        roots=[v]
    else:
        if q!=p :
            v=GAPMul(range(1,r+1),0)
            v[1-1]=1
            roots=[v]
        else:
            roots=[]
        v=GAPMul(range(1,r+1),0)
        v[1-1]=-ER(p)
        v[2-1]=1
        if r==2 and q>1 and q%2==1 :
            v=GAPMul(v,ER(p))
        roots.append(v)
    for i in range(2,r+1):
        v=GAPMul(range(1,r+1),0)
        v[i-1]=1
        v[i-1-1]=-1
        roots.append(v)
    return roots

ChevieData["imp"]["GeneratingRoots"]=lximp6

def lximp7(p,q,r):
    res=GAPMul(range(1,r+1),0)+GAPDiv(1,2)
    if q==1 :
        res[1-1]=GAPDiv(1,p)
    else:
        if q!=p :
            res=Concatenation([GAPDiv(q,p)],res)
    return res

ChevieData["imp"]["EigenvaluesGeneratingReflections"]=lximp7

def lximp8(p,q,r):
    rt=ChevieData["imp"]["GeneratingRoots"](p,q,r)
    rbar=ComplexConjugate(rt)
    e=ChevieData["imp"]["EigenvaluesGeneratingReflections"](p,q,r)
    e=1-map(lambda x: ER(Denominator(x))**Numerator(x),e)
    e=map(lambda i: GAPDiv(GAPMul(e[i-1],rbar[i-1]),GAPMul(rbar[i-1],rt[i-1])),range(1,len(e)+1))
    return map(lambda x: map(lambda y: GAPMul(x,y),rt),e)

ChevieData["imp"]["CartanMat"]=lximp8

def lximp9(p,q,r):
    return Concatenation(GAPMul(p,range(1,r-1+1)),[GAPDiv(GAPMul(r,p),q)])

ChevieData["imp"]["ReflectionDegrees"]=lximp9

def lximp10(p,q,r):
    res=GAPMul(p,range(0,r-1+1))
    if p==q and p>=2 and r>2 :
        res[r-1]=res[r-1]-r
    return res

ChevieData["imp"]["ReflectionCoDegrees"]=lximp10

def lximp11(p,q,r,s):
    if q==1 :
        if p==1 :
            if s==0 :
                return [[]]
            return map(lambda j: Concatenation(map(lambda k: Sum([j[k-1] for k in range(1,k-1+1)])+k-1+range(1,j[k-1]+1),range(1,len(j)+1))),Concatenation(map(lambda i: Partitions(s,i),range(1,r+1-s+1))))
        else:
            return Concatenation(map(lambda i: map(lambda j: Concatenation(range(1,i+1),i+1),ChevieData["imp"]["ParabolicRepresentatives"](1,1,r-i-1,s-i)),range(0,s+1)))
    else:
        if r==2 :
            if q==2 :
                t=[[[]],[[1],[2],[3]],[range(1,3+1)]]
                return t[s+1-1]
            else:
                if p==q :
                    if p%2==0 :
                        t=[[[]],[[1],[2]],[[1,2]]]
                        return t[s+1-1]
                    else:
                        t=[[],[1],[1,2]]
                        return t[s+1-1]
                else:
                    return false
        else:
            return false

ChevieData["imp"]["ParabolicRepresentatives"]=lximp11

def lximp12(p,q,r):
    if [q,r]==[2,2] :
        return GAPDiv(GAPMul(p,p+6),4)
    else:
        if q==1 :
            return NrPartitionTuples(r,p)
        else:
            return len(ChevieData["imp"]["ClassInfo"](p,q,r)["classtext"])

ChevieData["imp"]["NrConjugacyClasses"]=lximp12

def lximp13(p,q,r):
    def times(e,o):
        return Concatenation(map(lambda x: o,range(1,e+1)))
    
    
    if [q,r]==[2,2] and not "othermethod" in ChevieData :
        res={"classtext":[],
            "classparams":[],
            "classnames":[]}
        for i in range(0,p-1+1):
            for j in range(0,QuoInt(p-i-1,2)+1):
                res["classparams"].append(Concatenation(GAPMul(range(1,j+1),0)+1,GAPMul(range(1,i+1),0)))
                res["classtext"].append(Concatenation(GAPMul(range(1,j+1),0)+1,times(i,[1,2,3])))
                res["classnames"].append(just(Concatenation(times(j,"1"),times(i,"z"))))
        for j in [2,3]:
            for i in range(0,GAPDiv(p,2)-1+1):
                res["classparams"].append(Concatenation([j],GAPMul(range(1,i+1),0)))
                res["classtext"].append(Concatenation([j],times(i,[1,2,3])))
                res["classnames"].append(just(Concatenation(just(j),times(i,"z"))))
        res["malle"]=[]
        for a in range(0,p-1+1):
            res["malle"]+=map(lambda m: [3,a,m],range(0,QuoInt(p-a-1,2)+1))
        res["malle"]+=map(lambda m: [1,m],range(0,GAPDiv(p,2)-1+1))
        res["malle"]+=map(lambda m: [2,m],range(0,GAPDiv(p,2)-1+1))
        res["orders"]=map(lximp14,res["classparams"])
        res["classes"]=map(lximp15,res["classparams"])
        return res
    else:
        if q==1 :
            res={"classparams":PartitionTuples(r,p)}
            res["classtext"]=map(lximp16,res["classparams"])
            res["classnames"]=map(ChevieData["imp"]["ClassName"],res["classparams"])
            res["orders"]=map(lambda m: Lcm(map(lximp17,range(1,len(m)+1))),res["classparams"])
            res["centralizers"]=map(lambda m: GAPMul(p**Sum(m,Length),prod(map(lambda pp: prod(Collected(pp)),m))),res["classparams"])
            res["classes"]=map(lambda x: GAPDiv(GAPMul(p**r,Factorial(r)),x),res["centralizers"])
            return res
        else:
            def trans(w):
                d=0
                res=[]
                def word(l,i):
                    return map(lambda j: 1+j%2,i+range(l,1+1,l-1-l))
                
                
                def add(a):
                    l=len(res)
                    if l>0 and res[l-1]==a :
                        res=[res[k-1] for k in range(1,l-1+1)]
                    else:
                        if p==q and a in [1,2] and l>=q and [res[k-1] for k in range(l-q+1,l+1)]==word(q,3-a) :
                            res=Concatenation([res[k-1] for k in range(1,l-q+1)],word(q-1,3-a))
                        else:
                            res.append(a)
                
                
                for l in w:
                    if l==1 :
                        d=d+1
                    else:
                        if l!=2 :
                            add(l)
                        else:
                            d=d%p
                            if d==0 :
                                add(2)
                            else:
                                for i in range(1,p-d-1+1):
                                    add(1)
                                    add(2)
                                add(1)
                d=d%p
                if d%q!=0 :
                    Error()
                else:
                    if d!=0 :
                        res=Concatenation(1+res,GAPMul(range(1,GAPDiv(d,q)+1),0)+1)
                    else:
                        if p!=q :
                            res=1+res
                return res
            
            
            I=ChevieData["imp"]["ClassInfo"](p,1,r)
            res={"classtext":[],
                "classparams":[],
                "classnames":[],
                "orders":[],
                "centralizers":[]}
            for i in Filtered(range(1,len(I["classparams"])+1),lambda i: GAPMul(map(Length,I["classparams"][i-1]),range(0,p-1+1))%q==0):
                S=I["classparams"][i-1]
                a=Concatenation(S)
                a.append(q)
                a+=Filtered(range(1,p+1),lambda j: len(S[j-1])!=0)-1
                a=Gcd(*a)
                for j in range(0,a-1+1):
                    res["classtext"].append(trans(Concatenation(GAPMul(range(1,j+1),0)+1,I["classtext"][i-1],GAPMul(range(1,p-j+1),0)+1)))
                    if a>1 :
                        res["classparams"].append(Concatenation(S,[GAPDiv(GAPMul(p,j),a)]))
                    else:
                        res["classparams"].append(S)
                    res["orders"].append(I["orders"][i-1])
                    res["centralizers"].append(GAPDiv(GAPMul(I["centralizers"][i-1],a),q))
            res["classes"]=map(lambda x: GAPDiv(res["centralizers"][1-1],x),res["centralizers"])
            res["classnames"]=map(ChevieData["imp"]["ClassName"],res["classparams"])
            return res

def lximp14(c):
    if len(c)>0 and c[1-1] in [2,3] :
        return Lcm(2,GAPDiv(p,Gcd(Number(c,lambda x: x==0),p)))
    else:
        return Lcm(GAPDiv(p,Gcd(Number(c,lambda x: x==0),p)),GAPDiv(GAPDiv(p,2),Gcd(Number(c,lambda x: x==1),GAPDiv(p,2))))

def lximp15(c):
    if len(c)>0 and c[1-1] in [2,3] :
        return GAPDiv(p,q)
    else:
        if 1 in c :
            return 2
        else:
            return 1

def lximp16(S):
    S=Concatenation(map(lambda i: map(lambda t: [t,i-1],S[i-1]),range(1,p+1)))
    Sort(S,lximp18)
    l=0
    w=[]
    for d in S:
        w+=times(d[2-1],Concatenation(range(l+1,2+1,l-l+1),range(1,l+1+1)))
        w+=range(l+2,l+d[1-1]+1)
        l=l+d[1-1]
    return w

def lximp17(i):
    if len(m[i-1])==0 :
        return 1
    else:
        return Lcm(GAPDiv(GAPMul(m[i-1],p),Gcd(i-1,p)))

ChevieData["imp"]["ClassInfo"]=lximp13

def lximp19(p):
    if IsList(p) and ForAll(p,IsList) :
        if Sum(p,Sum)==1 :
            return Format(ER(len(p))**p.index([1])+1-1)
        else:
            return PartitionTupleToString(p)
    else:
        if IsList(p) and ForAll(p,IsInt) :
            return IntListToString(p)
        else:
            if IsList(p) and ForAll([p[k-1] for k in range(1,len(p)-1+1)],IsList) and IsInt(p[len(p)-1]) :
                j=GAPDiv(p[len(p)-1],len(p)-1)
                j=Format(ER(Denominator(j))**Numerator(j))
                if j=="1" :
                    j="+"
                else:
                    if j=="-1" :
                        j="-"
                return SPrint(PartitionTupleToString([p[k-1] for k in range(1,len(p)-1+1)]),j)
            else:
                Error()

ChevieData["imp"]["ClassName"]=lximp19

def lximp20(p,q,r):
    if q==1 :
        def pow(p,n):
            e=len(p)
            res=map(lambda x: [],range(1,e+1))
            for k in range(1,e+1):
                for l in p[k-1]:
                    g=Gcd(n,l)
                    for j in range(1,g+1):
                        res[1+GAPDiv(GAPMul(n,k-1),g)%e-1].append(GAPDiv(l,g))
            for k in range(1,e+1):
                Sort(res[k-1])
                res[k-1]=Reversed(res[k-1])
            return res
        
        
        pp=ChevieData["imp"]["ClassInfo"](p,q,r)["classparams"]
        res=[]
        for pw in Set(Factors(GAPMul(Factorial(r),p))):
            res[pw-1]=map(lambda x: pp.index(pow(x,pw))+1,pp)
        return res
    else:
        InfoChevie("# PowerMaps  !  implemented for G(",p,",",q,",",r,")\n")
        return false

ChevieData["imp"]["PowerMaps"]=lximp20

def lximp21(de,e,r):
    res={}
    d=GAPDiv(de,e)
    if e==1 :
        res["charparams"]=PartitionTuples(r,de)
        s=GAPMul(range(1,d+1),0)
        s[1-1]=1
        res["charSymbols"]=map(lambda x: SymbolPartitionTuple(x,s),res["charparams"])
    else:
        res["charparams"]=[]
        for t in PartitionTuples(r,de):
            tt=map(lambda i: Rotation(t,i),GAPMul(range(1,e+1),d))
            if t==Minimum(tt) :
                s=tt.index(t)+1
                if s==e :
                    res["charparams"].append(t)
                else:
                    t=[t[k-1] for k in range(1,GAPMul(s,d)+1)]
                    s=GAPDiv(e,s)
                    res["charparams"]+=map(lambda i: Concatenation(t,[s,i]),range(0,s-1+1))
        if d==1 :
            res["charSymbols"]=map(lambda x: SymbolPartitionTuple(x,0),res["charparams"])
        if d>1 and e%2==0 and r==2 :
            res["malle"]=map(lximp22,res["charparams"])
        else:
            if [de,e,r]==[3,3,3] :
                res["malle"]=[[2,3,2],[2,3,3],[2,3,1],[3,4],[3,5],[1,9],[3,2],[3,1],[2,3,4],[1,0]]
            else:
                if [de,e,r]==[3,3,4] :
                    res["malle"]=[[12,6],[4,10],[6,8],[4,11],[1,18],[12,3],[6,5,2],[8,4],[8,5],[6,5,1],[3,9],[6,2],[2,6],[4,2],[4,1],[3,3],[1,0]]
                else:
                    if [de,e,r]==[3,3,5] :
                        res["malle"]=[[30,10],[20,12],[5,19],[10,14],[10,15],[5,20],[1,30],[30,7,1],[40,6],[30,7,2],[10,11],[15,10],[20,9],[20,8],[15,11],[10,12],[4,18],[30,4],[20,5],[10,8],[10,7],[20,6],[5,12],[20,3],[10,6],[15,4],[15,5],[10,5],[6,9],[10,3],[10,2],[5,6],[5,2],[5,1],[4,3],[1,0]]
                    else:
                        if [de,e,r]==[4,4,3] :
                            res["malle"]=[[6,3],[3,6,1],[3,5],[3,6,2],[1,12],[3,2,1],[3,2,2],[3,1],[2,4],[1,0]]
    t=map(lximp23,range(r,0+1,r-1-r))
    if e>1 :
        t=map(lambda v: Minimum(map(lambda i: Rotation(v,GAPMul(i,d)),range(1,e+1))),t)
    res["extRefl"]=map(lambda v: res["charparams"].index(v)+1,t)
    if e==1 or d==1 :
        res["A"]=map(HighestPowerGenericDegreeSymbol,res["charSymbols"])
        res["a"]=map(LowestPowerGenericDegreeSymbol,res["charSymbols"])
        res["B"]=map(HighestPowerFakeDegreeSymbol,res["charSymbols"])
        res["b"]=map(LowestPowerFakeDegreeSymbol,res["charSymbols"])
    if e>1 and d>1 :
        res["opdam"]=PermListList(res["charparams"],map(lximp24,res["charparams"]))
    return res

def lximp22(t):
    if IsInt(t[len(t)-1]) :
        if t[len(t)-1]==0 :
            return [1,2,1,t.index([1])+1]
        else:
            return [1,1,2,t.index([1])+1]
    else:
        de=GAPDiv(len(t),2)
        pos=Filtered(range(1,len(t)+1),lambda i: len(t[i-1])>0)
        if len(pos)==1 :
            if t[pos[1-1]-1]==[2] :
                return [1,1,1,pos[1-1]-de]
            else:
                return [1,2,2,pos[1-1]-de]
        else:
            if pos[1-1]<=de :
                return [2,-1,pos[1-1],pos[2-1]-de]
            else:
                return [2,1,pos[2-1]-de,pos[1-1]-de]

def lximp23(i):
    v=map(lambda x: [],range(1,de+1))
    if i>0 :
        v[1-1]=[i]
    v[2-1]=GAPMul(range(1,r-i+1),0)+1
    return v

def lximp24(s):
    if not IsList(s[len(s)-1]) :
        s=Copy(s)
        t=GAPDiv(len(s)-2,d)
        for i,j in zip(GAPMul(range(0,t-1+1),d)+1,Rotation([s[k-1] for k in GAPMul(range(0,t-1+1),d)+1],1)):
            s[i]=j
        for i,j in zip(range(1,len(s)-2+1),Minimum(map(lambda i: Rotation([s[k-1] for k in range(1,len(s)-2+1)],GAPMul(i,d)),range(1,t+1)))):
            s[i]=j
        return s
    s=ShallowCopy(s)
    for i,j in zip(GAPMul(range(0,e-1+1),d)+1,Rotation([s[k-1] for k in GAPMul(range(0,e-1+1),d)+1],1)):
        s[i]=j
    return Minimum(map(lambda i: Rotation(s,GAPMul(i,d)),range(1,e+1)))

ChevieData["imp"]["CharInfo"]=lximp21

def lximp25(p,q,r):
    if q==1 or p==q :
        Error("should  !  be called")
    return false

ChevieData["imp"]["LowestPowerFakeDegrees"]=lximp25

def lximp26(p,q,r):
    if q==1 or p==q :
        Error("should  !  be called")
    return false

ChevieData["imp"]["HighestPowerFakeDegrees"]=lximp26

def lximp27(p,q,r):
    if q==1 :
        return SymbolsDefect(p,r,0,1)
    else:
        if q==p :
            ss=SymbolsDefect(p,r,0,0)
            res=[]
            for s in ss:
                p=[Rotations(s)[k-1] for k in range(2,len(s)+1)].index(s)+1
                if p==false :
                    res.append(s)
                else:
                    res+=map(lambda i: Concatenation(map(ShallowCopy,[s[k-1] for k in range(1,p+1)]),[GAPDiv(len(s),p),i]),range(0,GAPDiv(len(s),p)-1+1))
            return res
        else:
            return false

ChevieData["imp"]["CharSymbols"]=lximp27

def lximp28(p,q,r,c,v):
    if q==1 :
        c=CycPolFakeDegreeSymbol(SymbolPartitionTuple(c,1))
    else:
        if q==p :
            c=CycPolFakeDegreeSymbol(SymbolPartitionTuple(c,GAPMul(range(1,p+1),0)))
        else:
            return false
    return Value(c,v)

ChevieData["imp"]["FakeDegree"]=lximp28

def lximp29(p,q,r,s,option):
    if RankSymbol(s)==1 :
        return Format(ER(len(s))**s.index([1])+1-1,option)
    else:
        return PartitionTupleToString(s,option)

ChevieData["imp"]["CharName"]=lximp29

def lximp30(p,q,r,phi):
    if q==1 :
        def GenHooks(l,m):
            if len(l)==0 :
                return []
            m=AssociatedPartition(m)
            m+=GAPMul(range(1,l[1-1]-len(m)+1),0)
            m=1+m-range(1,len(m)+1)
            return Concatenation(map(lambda i: l[i-1]-i+[m[k-1] for k in range(1,l[i-1]+1)],range(1,len(l)+1)))
        
        
        res={"coeff":-1**GAPMul(r,p-1),
            "factor":GAPMul(range(1,p+1),0),
            "vcyc":[]}
        l=Concatenation(phi)
        Sort(l)
        res["factor"].append(GAPMul(range(1,len(l)+1)-len(l),l))
        for s in range(1,p+1):
            for t in range(1,p+1):
                for h in GenHooks(phi[s-1],phi[t-1]):
                    v=GAPMul(range(1,p+1),0)
                    if s!=t :
                        for i,j in zip([s,t],[1,-1]):
                            v[i]=j
                        v.append(h)
                        res["vcyc"].append([v,1])
                    else:
                        v.append(1)
                        for d in DivisorsInt(h):
                            if d>1 :
                                res["vcyc"].append([v,d])
        return res
    else:
        if [q,r]==[2,2] :
            ci=ChevieData["imp"]["CharInfo"](p,q,r)
            phi=ci["malle"][ci["charparams"].index(phi)+1-1]
            if phi[1-1]==1 :
                res={"coeff":1,
                    "factor":GAPMul(range(1,4+GAPDiv(p,2)+1),0),
                    "vcyc":[]}
                for l in [[1,-1,0,0],[0,0,1,-1]]:
                    l+=GAPMul(range(1,GAPDiv(p,2)+1),0)
                    res["vcyc"].append([l,1])
                for i in range(2,GAPDiv(p,2)+1):
                    for l in [[0,0,0,0,1],[1,-1,1,-1,1]]:
                        l+=GAPMul(range(1,GAPDiv(p,2)-1+1),0)
                        l[4+i-1]=-1
                        res["vcyc"].append([l,1])
            else:
                res={"coeff":-2,
                    "factor":GAPMul(range(1,4+GAPDiv(p,2)+1),0),
                    "vcyc":[],
                    "root":GAPMul(range(1,4+GAPDiv(p,2)+1),0)}
                res["rootCoeff"]=ER(GAPDiv(p,2))**2-phi[3-1]-phi[4-1]
                for i,j in zip(range(1,6+1),GAPDiv([1,1,1,1,1,1],2)):
                    res["root"][i]=j
                for i in range(3,GAPDiv(p,2)+1):
                    for j in [1,2]:
                        l=GAPMul(range(1,4+GAPDiv(p,2)+1),0)
                        for i,j in zip(4+[j,i],[1,-1]):
                            l[i]=j
                        res["vcyc"].append([l,1])
                if "old" in ChevieData :
                    for l in [[0,-1,0,-1,-1,0],[0,-1,-1,0,-1,0],[-1,0,-1,0,-1,0],[-1,0,0,-1,-1,0]]:
                        l+=GAPMul(range(1,GAPDiv(p,2)-2+1),0)
                        l.append(1)
                        res["vcyc"].append([l,1])
                else:
                    for l in [[0,-1,0,-1,-1,0],[0,-1,-1,0,0,-1],[-1,0,-1,0,-1,0],[-1,0,0,-1,0,-1]]:
                        l+=GAPMul(range(1,GAPDiv(p,2)-2+1),0)
                        l.append(1)
                        res["vcyc"].append([l,1])
            return res
        else:
            Error(" !  implemented")

ChevieData["imp"]["SchurModel"]=lximp30

def lximp31(p,q,r,phi):
    if [q,r]==[2,2] :
        ci=ChevieData["imp"]["CharInfo"](p,q,r)
        phi=ci["malle"][ci["charparams"].index(phi)+1-1]
        if phi[1-1]==1 :
            res={"order":[phi[2-1],3-phi[2-1],2+phi[3-1],5-phi[3-1],4+phi[4-1]]}
            res["order"]+=4+Difference(range(1,GAPDiv(p,2)+1),[phi[4-1]])
            return res
        else:
            res={"order":[1,2,3,4,4+phi[3-1],4+phi[4-1]]}
            res["order"]+=4+Difference(range(1,GAPDiv(p,2)+1),[phi[k-1] for k in [3,4]])
            res["rootPower"]=GAPMul(phi[2-1],ER(p)**phi[3-1]+phi[4-1]-2)
            return res
    else:
        Error(" !  implemented")

ChevieData["imp"]["SchurData"]=lximp31

def lximp32(p,q,r,phi,para,root):
    if r==1 :
        return VcycSchurElement(Concatenation(para[1-1],[0]),ChevieData["imp"]["SchurModel"](p,q,r,phi))
    else:
        if p==1 :
            return VcycSchurElement([0,GAPDiv(-para[1-1][1-1],para[1-1][2-1])],ChevieData["imp"]["SchurModel"](p,q,r,phi))
        else:
            if q==1 :
                return VcycSchurElement(Concatenation(para[1-1],[GAPDiv(-para[2-1][1-1],para[2-1][2-1])]),ChevieData["imp"]["SchurModel"](p,q,r,phi))
            else:
                if [q,r]==[2,2] :
                    return VcycSchurElement(Concatenation([para[k-1] for k in [2,3,1]]),ChevieData["imp"]["SchurModel"](p,q,r,phi),ChevieData["imp"]["SchurData"](p,q,r,phi))
                else:
                    if p==q :
                        if IsInt(phi[len(phi)-1]) :
                            m=len(phi)-2
                            phi=FullSymbol(phi)
                        else:
                            m=p
                        return GAPDiv(ChevieData["imp"]["SchurElement"](p,1,r,phi,Concatenation([map(lambda i: ER(p)**i,range(0,p-1+1))],[para[k-1] for k in range(2,len(para)+1)]),[]),m)
                    else:
                        if para[2-1]==para[3-1] :
                            if IsInt(phi[len(phi)-1]) :
                                m=len(phi)-2
                                phi=FullSymbol(phi)
                            else:
                                m=p
                            if para[1-1]==map(lambda i: ER(GAPDiv(p,q))**i-1,range(1,GAPDiv(p,q)+1)) :
                                para=[map(lambda i: ER(p)**i,range(0,p-1+1)),para[2-1]]
                            else:
                                para=[Concatenation(Matrix(map(lambda i: GAPMul(map(lambda j: ER(q)**j,range(0,q-1+1)),GetRoot(i,q)),para[1-1])).transpose()),para[2-1]]
                            return GAPDiv(GAPMul(GAPDiv(p,q),ChevieData["imp"]["SchurElement"](p,1,r,phi,para,[])),m)
                        else:
                            ChevieData["compat"]["InfoChevie"]("# SchurElements(H(G(",p,",",q,",",r,"),",para,")  !  implemented\n")
                            return false

ChevieData["imp"]["SchurElement"]=lximp32

def lximp33(p,q,r,phi,para,root):
    if r==1 :
        return VFactorSchurElement(Concatenation(para[1-1],[0]),ChevieData["imp"]["SchurModel"](p,q,r,phi))
    else:
        if p==1 :
            return VFactorSchurElement([0,GAPDiv(-para[1-1][1-1],para[1-1][2-1])],ChevieData["imp"]["SchurModel"](p,q,r,phi))
        else:
            if q==1 :
                return VFactorSchurElement(Concatenation(para[1-1],[GAPDiv(-para[2-1][1-1],para[2-1][2-1])]),ChevieData["imp"]["SchurModel"](p,q,r,phi))
            else:
                if [q,r]==[2,2] :
                    return VFactorSchurElement(Concatenation([para[k-1] for k in [2,3,1]]),ChevieData["imp"]["SchurModel"](p,q,r,phi),ChevieData["imp"]["SchurData"](p,q,r,phi))
                else:
                    if p==q :
                        if IsInt(phi[len(phi)-1]) :
                            m=len(phi)-2
                            phi=FullSymbol(phi)
                        else:
                            m=p
                        F=ChevieData["imp"]["FactorizedSchurElement"](p,1,r,phi,Concatenation([map(lambda i: ER(p)**i,range(0,p-1+1))],[para[k-1] for k in range(2,len(para)+1)]),[])
                        F["factor"]=GAPDiv(F["factor"],m)
                        return F
                    else:
                        if para[2-1]==para[3-1] :
                            if IsInt(phi[len(phi)-1]) :
                                m=len(phi)-2
                                phi=FullSymbol(phi)
                            else:
                                m=p
                            if para[1-1]==map(lambda i: ER(GAPDiv(p,q))**i-1,range(1,GAPDiv(p,q)+1)) :
                                para=[map(lambda i: ER(p)**i,range(0,p-1+1)),para[2-1]]
                            else:
                                para=[Concatenation(Matrix(map(lambda i: GAPMul(map(lambda j: ER(q)**j,range(0,q-1+1)),GetRoot(i,q)),para[1-1])).transpose()),para[2-1]]
                            F=ChevieData["imp"]["FactorizedSchurElement"](p,1,r,phi,para,[])
                            F["factor"]=GAPMul(GAPDiv(p,GAPMul(q,m)),F["factor"])
                            return F
                        else:
                            ChevieData["compat"]["InfoChevie"]("# FactorizedSchurElements(H(G(",p,",",q,",",r,"),",para,")  !  implemented\n")
                            return false

ChevieData["imp"]["FactorizedSchurElement"]=lximp33

def lximp34(p,q,r,para,root):
    res={}
    res["name"]=SPrint("H(G(",p,",",q,",",r,"))")
    res["identifier"]=res["name"]
    res["degrees"]=ChevieData["imp"]["ReflectionDegrees"](p,q,r)
    res["size"]=prod(res["degrees"])
    res["order"]=prod(res["degrees"])
    res["dim"]=r
    ci=ChevieData["imp"]["CharInfo"](p,q,r)
    if r==1 :
        res["reflclasses"]=[2]
        res["orders"]=map(lambda i: GAPDiv(p,Gcd(i,p)),range(0,p-1+1))
        res["irreducibles"]=map(lambda i: map(lambda j: para[1-1][i-1]**j,range(0,p-1+1)),range(1,p+1))
        res["classes"]=GAPMul(range(1,p+1),0)+1
        res["powermap"]=ChevieData["imp"]["PowerMaps"](p,q,r)
    else:
        if q==1 :
            cl=ChevieData["imp"]["ClassInfo"](p,q,r)
            d=map(lambda x: [],range(1,p+1))
            d[1-1]=[p-1,1]
            res["reflclasses"]=[cl["classparams"].index(d)+1]
            d[1-1]=GAPMul(range(1,p-1+1),0)+1
            d[2-1]=[1]
            res["reflclasses"].append(cl["classparams"].index(d)+1)
            res.update(cl)
            res["powermap"]=ChevieData["imp"]["PowerMaps"](p,q,r)
            def HooksBeta(S,s):
                res=[]
                e=len(S)
                if e==0 :
                    return res
                j=e
                for i in range(S[e-1]-1,0+1,S[e-1]-2-S[e-1]-1):
                    if not i in S :
                        while j>0 and S[j-1]>i:
                            j=j-1
                        k=j+1
                        while k<=e and S[k-1]-i<=s:
                            z=[i]
                            z+=[S[k-1] for k in range(j+1,k-1+1)]
                            zi=Filtered(range(2,len(z)+1),lambda i: z[i-1]-z[i-1-1]>1)
                            res.append({"area":S[k-1]-i,
                                "hooklength":k-j-1,
                                "start":S[k-1],
                                "startpos":k,
                                "stoppos":j+1,
                                "DC":[z[k-1] for k in zi]-e,
                                "SC":[z[k-1] for k in Concatenation(zi-1,[len(z)])]+1-e})
                            k=k+1
                return res
            
            
            def StripsBeta(S,s):
                res=[[]]
                for hook in HooksBeta(S,s):
                    if s==hook["area"] :
                        res.append([hook])
                    else:
                        j=hook["stoppos"]-1-len(S)
                        for hs in StripsBeta([S[k-1] for k in range(1,hook["stoppos"]-1+1)],s-hook["area"]):
                            for h in hs:
                                h["SC"]=h["SC"]+j
                                h["DC"]=h["DC"]+j
                            hs.append(hook)
                            res.append(hs)
                return res
            
            
            StripsCache={}
            def code(arg):
                res=[]
                for S in arg:
                    for p in S:
                        res+=p
                        res.append(-1)
                p=".0123456789abcdefghijklmnopqrstuvwxyz"
                return [p[k-1] for k in 2+res]
            
            
            def Strips(S,s):
                def apply(S,hs):
                    S=ShallowCopy(S)
                    for h in hs:
                        for i,j in zip(range(h["stoppos"],h["startpos"]+1),Concatenation([h["start"]-h["area"]],[S[k-1] for k in range(h["stoppos"],h["startpos"]-1+1)])):
                            S[i]=j
                    while len(S)>0 and S[1-1]==0:
                        S=[S[k-1] for k in range(2,len(S)+1)]-1
                    return S
                
                
                e=len(S)
                if e==0 :
                    if s==0 :
                        return [{"SC":[],
                            "DC":[],
                            "cc":0,
                            "hooklength":0,
                            "area":0,
                            "remainder":[]}]
                    else:
                        return []
                name=code(S,[[s]])
                if name in StripsCache :
                    return StripsCache[name]
                res=[]
                for hs in StripsBeta(S[e-1],s):
                    hs={"area":Sum(hs,lambda x: x["area"]),
                        "cc":len(hs),
                        "hooklength":Sum(hs,lambda x: x["hooklength"]),
                        "SC":Concatenation(map(lambda x: map(lambda y: [e,y],x["SC"]),hs)),
                        "DC":Concatenation(map(lambda x: map(lambda y: [e,y],x["DC"]),hs)),
                        "remainder":apply(S[e-1],hs)}
                    for a in Strips([S[k-1] for k in range(1,e-1+1)],s-hs["area"]):
                        ss={}
                        for r in RecFields(a):
                            ss[r]=ShallowCopy(a[r])
                        ss["SC"]+=hs["SC"]
                        ss["DC"]+=hs["DC"]
                        ss["remainder"].append(hs["remainder"])
                        ss["cc"]=ss["cc"]+hs["cc"]
                        ss["hooklength"]=ss["hooklength"]+hs["hooklength"]
                        ss["area"]=ss["area"]+hs["area"]
                        res.append(ss)
                StripsCache[name]=res
                return res
            
            
            def Delta(k,hs,Q,v):
                res=1
                if hs["cc"]>1 :
                    if k==1 or Q[1-1]==-Q[2-1] :
                        return 0
                    else:
                        res=GAPMul(res,Q[1-1]+Q[2-1]**hs["cc"]-1)
                q=GAPDiv(-Q[1-1],Q[2-1])
                res=GAPMul(res,-1**hs["hooklength"])
                if k==0 :
                    return res
                ctSC=map(lambda x: GAPMul(v[x[1-1]-1],q**x[2-1]),hs["SC"])
                ctDC=map(lambda x: GAPMul(v[x[1-1]-1],q**x[2-1]),hs["DC"])
                res=GAPMul(res,prod(ctSC))
                if k==1 :
                    return res
                def ElementarySymmetricFunction(t,v):
                    return Sum(Combinations(range(1,len(v)+1),t),lambda x: prod([v[k-1] for k in x]))
                
                
                def HomogeneousSymmetricFunction(t,v):
                    return Sum(Combinations(Concatenation(map(lambda x: range(1,len(v)+1),range(1,t+1))),t),lambda x: prod([v[k-1] for k in x]))
                
                
                return GAPMul(res,-1**hs["cc"]-1)
            
            
            chiCache={}
            LIM=r
            def GenericEntry(lambda_,mu):
                n=Sum(lambda_,Sum)
                if n==0 :
                    return 1
                if n<LIM :
                    name=code(lambda_,mu)
                    if name in chiCache :
                        return chiCache[name]
                bp=Maximum(Concatenation(lambda_))
                i=PositionProperty(lambda_,lambda x: bp in x)
                rest=ShallowCopy(lambda_)
                rest[i-1]=[rest[i-1][k-1] for k in range(2,len(rest[i-1])+1)]
                res=GAPMul(-prod(para[2-1])**GAPMul(i-1,n-bp),Sum(Strips(mu,bp),lximp35))
                if n<LIM :
                    chiCache[name]=res
                return res
            
            
            res["irreducibles"]=map(lambda x: map(lambda y: GenericEntry(y,x),cl["classparams"]),map(lambda x: map(BetaSet,x),cl["classparams"]))
        else:
            if [q,r]==[2,2] and not "othermethod" in ChevieData :
                cl=ChevieData["imp"]["ClassInfo"](p,q,r)
                X=para[2-1]
                Y=para[3-1]
                Z=para[1-1]
                def GenericEntry(char,class_):
                    char=ci["malle"][ci["charparams"].index(char)+1-1]
                    if char[1-1]==1 :
                        w=[Z[char[4-1]-1],X[char[2-1]-1],Y[char[3-1]-1]]
                        return prod(class_)
                    else:
                        w=GAPMul(char[2-1],GetRoot(GAPMul(X[1-1],X[2-1]),2))
                        class_=map(lambda i: Number(class_,lambda j: i==j),range(0,3+1))
                        if class_[2-1]>0 :
                            char=Sum([Z[k-1] for k in [char[k-1] for k in [3,4]]],lambda x: x**class_[2-1])
                        else:
                            if class_[3-1]>0 :
                                char=Sum(X)
                            else:
                                if class_[4-1]>0 :
                                    char=Sum(Y)
                                else:
                                    char=2
                        return GAPMul(w**class_[1-1],char)
                
                
                res["classes"]=cl["classes"]
                res["orders"]=cl["orders"]
                res["irreducibles"]=map(lambda char: map(lambda class_: GenericEntry(char,class_),cl["classparams"]),ci["charparams"])
            else:
                cl=ChevieData["imp"]["ClassInfo"](p,q,r)
                res["centralizers"]=cl["centralizers"]
                res["orders"]=cl["orders"]
                res["classes"]=map(lambda x: GAPDiv(res["size"],x),res["centralizers"])
                res["irreducibles"]=map(lambda i: CharRepresentationWords(ChevieData["imp"]["HeckeRepresentation"](p,q,r,para,[],i),cl["classtext"]),range(1,len(res["classes"])+1))
    res["centralizers"]=map(lambda x: GAPDiv(res["size"],x),res["classes"])
    res["parameter"]=para
    res["irreducibles"]=GAPMul(res["irreducibles"],prod(para)**0)
    return ChevieData["compat"]["MakeCharacterTable"](res)

def lximp35(x):
    d=Delta(i-1,x,para[2-1],para[1-1])
    if d==0 :
        return d
    else:
        return GAPMul(d,GenericEntry(rest,x["remainder"]))

ChevieData["imp"]["HeckeCharTable"]=lximp34

