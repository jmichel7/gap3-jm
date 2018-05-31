
def ximp1(arg):
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
      if !=(p,3) :
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
        if >=(r,3) :
          print "\\",
        print "\n",
        if GAPDiv(p,q)>2 :
          print g(len(just(GAPDiv(p,q)))+5),"(",GAPDiv(p,q),")",
        else:
          print g(3),
        print indices[1-1],"  | ",
        for j in range(3,r+1):
          print indices[j+1-1],
          if !=(j,r) :
            print "-",
        print "\n",g(2),"\\ |",
        if >=(r,3) :
          print "/",
        print "\n",just(indices[3-1],indent+1),"   ",IntListToString([indices[k-1] for k in [1,2,3]]),"==",IntListToString([indices[k-1] for k in [2,3,1]]),"==",IntListToString([indices[k-1] for k in [3,1,2]]),"\n",
      else:
        print indices[2-1],"\n",g(2),"/",q+1," ",
        if >=(r,3) :
          print "\\",
        print "\n",
        if GAPDiv(p,q)>2 :
          print g(len(SPrint(GAPDiv(p,q)))+5),"(",GAPDiv(p,q),")",
        else:
          print g(3),
        print indices[1-1],"   ",
        if >=(r,3) :
          print "==",
        for j in range(3,r+1):
          print indices[j+1-1],
          if !=(j,r) :
            print "-",
        print "\n",g(2),"\\  ",
        if >=(r,3) :
          print "/",
        print "\n",just(indices[3-1],indent+1),
        j=ChevieData["imp"]["BraidRelations"](p,q,r)
        for g in range(1,Minimum(3,r)+1):
          print "   ",IntListToString([indices[k-1] for k in j[g][1-1]]),"==",IntListToString([indices[k-1] for k in j[g][2-1]]),
        print "\n",

ChevieData["imp"]["PrintDiagram"]=ximp1

def ximp2(p,q,r):
  return r

ChevieData["imp"]["SemisimpleRank"]=ximp2

def ximp3(p,q,r):
  def b(i,j,o):
    def p(i,j):
      return map(lambda k: GAPMul(i,k%2)+GAPMul(j,1-k%2),range(1,o+1))
    
    
    return [p(i,j),p(j,i)]
  
  
  res=[]
  if q==1 :
    if >=(r,2) :
      if p==1 :
        res.append(b(1,2,3))
      else:
        res.append(b(1,2,4))
    res.extend(map(lambda i: b(i,i-1,3),range(3,r+1)))
    for i in range(3,r+1):
      res.extend(map(lambda j: b(i,j,2),range(1,i-2+1)))
  else:
    if p==q :
      res.append(b(1,2,p))
      if >=(r,3) :
        res.extend([[[1,2,3,1,2,3],[3,1,2,3,1,2]],b(1,3,3),b(2,3,3)])
      res.extend(map(lambda i: b(i,i-1,3),range(4,r+1)))
      for i in range(4,r+1):
        res.extend(map(lambda j: b(i,j,2),range(1,i-2+1)))
    else:
      res.append([[1,2,3],[2,3,1]])
      i=b(2,3,q-1)
      res.append([Concatenation([1,2],i[2-1]),Concatenation([3,1],i[1-1])])
      if >=(r,3) :
        if !=(q,2) :
          res.append([[2,3,4,2,3,4],[4,2,3,4,2,3]])
        res.extend([b(2,4,3),b(3,4,3),b(1,4,2)])
      res.extend(map(lambda i: b(i,i-1,3),range(5,r+1+1)))
      for i in range(5,r+1+1):
        res.extend(map(lambda j: b(i,j,2),range(1,i-2+1)))
  return res

ChevieData["imp"]["BraidRelations"]=ximp3

def ximp4(p,q,r):
  return GAPDiv(GAPMul(p**r,Factorial(r)),q)

ChevieData["imp"]["Size"]=ximp4

def ximp5(arg):
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

ChevieData["imp"]["ReflectionName"]=ximp5

def ximp6(p,q,r):
  if q==1 :
    v=GAPMul(range(1,r+1),0)
    v[1-1]=1
    roots=[v]
  else:
    if !=(q,p) :
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

ChevieData["imp"]["GeneratingRoots"]=ximp6

def ximp7(p,q,r):
  res=GAPMul(range(1,r+1),0)+GAPDiv(1,2)
  if q==1 :
    res[1-1]=GAPDiv(1,p)
  else:
    if !=(q,p) :
      res=Concatenation([GAPDiv(q,p)],res)
  return res

ChevieData["imp"]["EigenvaluesGeneratingReflections"]=ximp7

def ximp8(p,q,r):
  rt=ChevieData["imp"]["GeneratingRoots"](p,q,r)
  rbar=ComplexConjugate(rt)
  e=ChevieData["imp"]["EigenvaluesGeneratingReflections"](p,q,r)
  e=1-map(lambda x: ER(Denominator(x))**Numerator(x),e)
  e=map(lambda i: GAPDiv(GAPMul(e[i-1],rbar[i-1]),GAPMul(rbar[i-1],rt[i-1])),range(1,len(e)+1))
  return map(lambda x: map(lambda y: GAPMul(x,y),rt),e)

ChevieData["imp"]["CartanMat"]=ximp8

def ximp9(p,q,r):
  return Concatenation(GAPMul(p,range(1,r-1+1)),[GAPDiv(GAPMul(r,p),q)])

ChevieData["imp"]["ReflectionDegrees"]=ximp9

def ximp10(p,q,r):
  res=GAPMul(p,range(0,r-1+1))
  if p==q and >=(p,2) and r>2 :
    res[r-1]=res[r-1]-r
  return res

ChevieData["imp"]["ReflectionCoDegrees"]=ximp10

def ximp11(p,q,r,s):
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

ChevieData["imp"]["ParabolicRepresentatives"]=ximp11

def ximp12(p,q,r):
  if [q,r]==[2,2] :
    return GAPDiv(GAPMul(p,p+6),4)
  else:
    if q==1 :
      return NrPartitionTuples(r,p)
    else:
      return len(ChevieData["imp"]["ClassInfo"](p,q,r)["classtext"])

ChevieData["imp"]["NrConjugacyClasses"]=ximp12

def ximp13(p,q,r):
  def times(e,o):
    return Concatenation(map(lambda x: o,range(1,e+1)))
  
  
  if [q,r]==[2,2] and !("othermethod" in CHEVIE) :
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
      res["malle"].extend(map(lambda m: [3,a,m],range(0,QuoInt(p-a-1,2)+1)))
    res["malle"].extend(map(lambda m: [1,m],range(0,GAPDiv(p,2)-1+1)))
    res["malle"].extend(map(lambda m: [2,m],range(0,GAPDiv(p,2)-1+1)))
    res["orders"]=map(ximp14,res["classparams"])
    res["classes"]=map(ximp15,res["classparams"])
    return res
  else:
    if q==1 :
      res={"classparams":PartitionTuples(r,p)}
      res["classtext"]=map(ximp16,res["classparams"])
      res["classnames"]=map(ChevieData["imp"]["ClassName"],res["classparams"])
      res["orders"]=map(lambda m: Lcm(map(ximp17,range(1,len(m)+1))),res["classparams"])
      res["centralizers"]=map(lambda m: GAPMul(p**Sum(m,Length),prod(map(lambda pp: prod(Collected(pp)),m))),res["classparams"])
      res["classes"]=map(lambda x: GAPDiv(GAPMul(p**r,Factorial(r)),x),res["centralizers"])
      return res
    else:
      def trans(w):
        d=0
        res=[]
        def word(l,i):
          return map(lambda j: 1+j%2,i+[l,..(l-1,1)])
        
        
        def add(a):
          l=len(res)
          if l>0 and res[l-1]==a :
            res=[res[k-1] for k in range(1,l-1+1)]
          else:
            if p==q and in(a,[1,2]) and >=(l,q) and [res[k-1] for k in range(l-q+1,l+1)]==word(q,3-a) :
              res=Concatenation([res[k-1] for k in range(1,l-q+1)],word(q-1,3-a))
            else:
              res.append(a)
        
        
        for l in w:
          if l==1 :
            d=d+1
          else:
            if !=(l,2) :
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
        if !=(d%q,0) :
          Error()
        else:
          if !=(d,0) :
            res=Concatenation(1+res,GAPMul(range(1,GAPDiv(d,q)+1),0)+1)
          else:
            if !=(p,q) :
              res=1+res
        return res
      
      
      I=ChevieData["imp"]["ClassInfo"](p,1,r)
      res={"classtext":[],
        "classparams":[],
        "classnames":[],
        "orders":[],
        "centralizers":[]}
      for i in Filtered(range(1,len(I["classparams"])+1),lambda i: GAPMul(map(Length,I.classparams[i-1]),range(0,p-1+1))%q==0):
        S=I.classparams[i-1]
        a=Concatenation(S)
        a.append(q)
        a.extend(Filtered(range(1,p+1),lambda j: !=(len(S[j-1]),0))-1)
        a=ApplyFunc(Gcd,a)
        for j in range(0,a-1+1):
          res["classtext"].append(trans(Concatenation(GAPMul(range(1,j+1),0)+1,I.classtext[i-1],GAPMul(range(1,p-j+1),0)+1)))
          if a>1 :
            res["classparams"].append(Concatenation(S,[GAPDiv(GAPMul(p,j),a)]))
          else:
            res["classparams"].append(S)
          res["orders"].append(I.orders[i-1])
          res["centralizers"].append(GAPDiv(GAPMul(I.centralizers[i-1],a),q))
      res["classes"]=map(lambda x: GAPDiv(res.centralizers[1-1],x),res["centralizers"])
      res["classnames"]=map(ChevieData["imp"]["ClassName"],res["classparams"])
      return res

def ximp14(c):
  if len(c)>0 and in(c[1-1],[2,3]) :
    return Lcm(2,GAPDiv(p,Gcd(Number(c,lambda x: x==0),p)))
  else:
    return Lcm(GAPDiv(p,Gcd(Number(c,lambda x: x==0),p)),GAPDiv(GAPDiv(p,2),Gcd(Number(c,lambda x: x==1),GAPDiv(p,2))))

def ximp15(c):
  if len(c)>0 and in(c[1-1],[2,3]) :
    return GAPDiv(p,q)
  else:
    if in(1,c) :
      return 2
    else:
      return 1

def ximp16(S):
  S=Concatenation(map(lambda i: map(lambda t: [t,i-1],S[i-1]),range(1,p+1)))
  Sort(S,ximp18)
  l=0
  w=[]
  for d in S:
    w.extend(times(d[2-1],Concatenation([l+1,..(l,2)],range(1,l+1+1))))
    w.extend(range(l+2,l+d[1-1]+1))
    l=l+d[1-1]
  return w

def ximp17(i):
  if len(m[i-1])==0 :
    return 1
  else:
    return Lcm(GAPDiv(GAPMul(m[i-1],p),Gcd(i-1,p)))

ChevieData["imp"]["ClassInfo"]=ximp13

def ximp19(p):
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

ChevieData["imp"]["ClassName"]=ximp19

def ximp20(p,q,r):
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

ChevieData["imp"]["PowerMaps"]=ximp20

def ximp21(de,e,r):
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
          res["charparams"].extend(map(lambda i: Concatenation(t,[s,i]),range(0,s-1+1)))
    if d==1 :
      res["charSymbols"]=map(lambda x: SymbolPartitionTuple(x,0),res["charparams"])
    if d>1 and e%2==0 and r==2 :
      res["malle"]=map(ximp22,res["charparams"])
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
  t=map(ximp23,[r,..(r-1,0)])
  if e>1 :
    t=map(lambda v: Minimum(map(lambda i: Rotation(v,GAPMul(i,d)),range(1,e+1))),t)
  res["extRefl"]=map(lambda v: res["charparams"].index(v)+1,t)
  if e==1 or d==1 :
    res["A"]=map(HighestPowerGenericDegreeSymbol,res["charSymbols"])
    res["a"]=map(LowestPowerGenericDegreeSymbol,res["charSymbols"])
    res["B"]=map(HighestPowerFakeDegreeSymbol,res["charSymbols"])
    res["b"]=map(LowestPowerFakeDegreeSymbol,res["charSymbols"])
  if e>1 and d>1 :
    res["opdam"]=PermListList(res["charparams"],map(ximp24,res["charparams"]))
  return res

def ximp22(t):
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
      if <=(pos[1-1],de) :
        return [2,-1,pos[1-1],pos[2-1]-de]
      else:
        return [2,1,pos[2-1]-de,pos[1-1]-de]

def ximp23(i):
  v=map(lambda x: [],range(1,de+1))
  if i>0 :
    v[1-1]=[i]
  v[2-1]=GAPMul(range(1,r-i+1),0)+1
  return v

def ximp24(s):
  if !(IsList(s[len(s)-1])) :
    s=Copy(s)
    t=GAPDiv(len(s)-2,d)
    [s[k-1] for k in GAPMul(range(0,t-1+1),d)+1]=Rotation([s[k-1] for k in GAPMul(range(0,t-1+1),d)+1],1)
    [s[k-1] for k in range(1,len(s)-2+1)]=Minimum(map(lambda i: Rotation([s[k-1] for k in range(1,len(s)-2+1)],GAPMul(i,d)),range(1,t+1)))
    return s
  s=ShallowCopy(s)
  [s[k-1] for k in GAPMul(range(0,e-1+1),d)+1]=Rotation([s[k-1] for k in GAPMul(range(0,e-1+1),d)+1],1)
  return Minimum(map(lambda i: Rotation(s,GAPMul(i,d)),range(1,e+1)))

ChevieData["imp"]["CharInfo"]=ximp21

def ximp25(p,q,r):
  if q==1 or p==q :
    Error("should  !  be called")
  return false

ChevieData["imp"]["LowestPowerFakeDegrees"]=ximp25

def ximp26(p,q,r):
  if q==1 or p==q :
    Error("should  !  be called")
  return false

ChevieData["imp"]["HighestPowerFakeDegrees"]=ximp26

def ximp27(p,q,r):
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
          res.extend(map(lambda i: Concatenation(map(ShallowCopy,[s[k-1] for k in range(1,p+1)]),[GAPDiv(len(s),p),i]),range(0,GAPDiv(len(s),p)-1+1)))
      return res
    else:
      return false

ChevieData["imp"]["CharSymbols"]=ximp27

def ximp28(p,q,r,c,v):
  if q==1 :
    c=CycPolFakeDegreeSymbol(SymbolPartitionTuple(c,1))
  else:
    if q==p :
      c=CycPolFakeDegreeSymbol(SymbolPartitionTuple(c,GAPMul(range(1,p+1),0)))
    else:
      return false
  return Value(c,v)

ChevieData["imp"]["FakeDegree"]=ximp28

def ximp29(p,q,r,s,option):
  if RankSymbol(s)==1 :
    return Format(ER(len(s))**s.index([1])+1-1,option)
  else:
    return PartitionTupleToString(s,option)

ChevieData["imp"]["CharName"]=ximp29

def ximp30(p,q,r,phi):
  if q==1 :
    def GenHooks(l,m):
      if len(l)==0 :
        return []
      m=AssociatedPartition(m)
      m.extend(GAPMul(range(1,l[1-1]-len(m)+1),0))
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
          if !=(s,t) :
            [v[k-1] for k in [s,t]]=[1,-1]
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
      phi=ci.malle[ci["charparams"].index(phi)+1-1]
      if phi[1-1]==1 :
        res={"coeff":1,
          "factor":GAPMul(range(1,4+GAPDiv(p,2)+1),0),
          "vcyc":[]}
        for l in [[1,-1,0,0],[0,0,1,-1]]:
          l.extend(GAPMul(range(1,GAPDiv(p,2)+1),0))
          res["vcyc"].append([l,1])
        for i in range(2,GAPDiv(p,2)+1):
          for l in [[0,0,0,0,1],[1,-1,1,-1,1]]:
            l.extend(GAPMul(range(1,GAPDiv(p,2)-1+1),0))
            l[4+i-1]=-1
            res["vcyc"].append([l,1])
      else:
        res={"coeff":-2,
          "factor":GAPMul(range(1,4+GAPDiv(p,2)+1),0),
          "vcyc":[],
          "root":GAPMul(range(1,4+GAPDiv(p,2)+1),0)}
        res["rootCoeff"]=ER(GAPDiv(p,2))**2-phi[3-1]-phi[4-1]
        [res.root[k-1] for k in range(1,6+1)]=GAPDiv([1,1,1,1,1,1],2)
        for i in range(3,GAPDiv(p,2)+1):
          for j in [1,2]:
            l=GAPMul(range(1,4+GAPDiv(p,2)+1),0)
            [l[k-1] for k in 4+[j,i]]=[1,-1]
            res["vcyc"].append([l,1])
        if "old" in CHEVIE :
          for l in [[0,-1,0,-1,-1,0],[0,-1,-1,0,-1,0],[-1,0,-1,0,-1,0],[-1,0,0,-1,-1,0]]:
            l.extend(GAPMul(range(1,GAPDiv(p,2)-2+1),0))
            l.append(1)
            res["vcyc"].append([l,1])
        else:
          for l in [[0,-1,0,-1,-1,0],[0,-1,-1,0,0,-1],[-1,0,-1,0,-1,0],[-1,0,0,-1,0,-1]]:
            l.extend(GAPMul(range(1,GAPDiv(p,2)-2+1),0))
            l.append(1)
            res["vcyc"].append([l,1])
      return res
    else:
      Error(" !  implemented")

ChevieData["imp"]["SchurModel"]=ximp30

def ximp31(p,q,r,phi):
  if [q,r]==[2,2] :
    ci=ChevieData["imp"]["CharInfo"](p,q,r)
    phi=ci.malle[ci["charparams"].index(phi)+1-1]
    if phi[1-1]==1 :
      res={"order":[phi[2-1],3-phi[2-1],2+phi[3-1],5-phi[3-1],4+phi[4-1]]}
      res["order"].extend(4+Difference(range(1,GAPDiv(p,2)+1),[phi[4-1]]))
      return res
    else:
      res={"order":[1,2,3,4,4+phi[3-1],4+phi[4-1]]}
      res["order"].extend(4+Difference(range(1,GAPDiv(p,2)+1),[phi[k-1] for k in [3,4]]))
      res["rootPower"]=GAPMul(phi[2-1],ER(p)**phi[3-1]+phi[4-1]-2)
      return res
  else:
    Error(" !  implemented")

ChevieData["imp"]["SchurData"]=ximp31

def ximp32(p,q,r,phi,para,root):
  if r==1 :
    return VcycSchurElement(Concatenation(para[1-1],[0]),ChevieData["imp"]["SchurModel"](p,q,r,phi))
  else:
    if p==1 :
      return VcycSchurElement([0,GAPDiv(-para[1][1-1],para[1][2-1])],ChevieData["imp"]["SchurModel"](p,q,r,phi))
    else:
      if q==1 :
        return VcycSchurElement(Concatenation(para[1-1],[GAPDiv(-para[2][1-1],para[2][2-1])]),ChevieData["imp"]["SchurModel"](p,q,r,phi))
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

ChevieData["imp"]["SchurElement"]=ximp32

def ximp33(p,q,r,phi,para,root):
  if r==1 :
    return VFactorSchurElement(Concatenation(para[1-1],[0]),ChevieData["imp"]["SchurModel"](p,q,r,phi))
  else:
    if p==1 :
      return VFactorSchurElement([0,GAPDiv(-para[1][1-1],para[1][2-1])],ChevieData["imp"]["SchurModel"](p,q,r,phi))
    else:
      if q==1 :
        return VFactorSchurElement(Concatenation(para[1-1],[GAPDiv(-para[2][1-1],para[2][2-1])]),ChevieData["imp"]["SchurModel"](p,q,r,phi))
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

ChevieData["imp"]["FactorizedSchurElement"]=ximp33

def ximp34(p,q,r,para,root):
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
    res["irreducibles"]=map(lambda i: map(lambda j: para[1][i-1]**j,range(0,p-1+1)),range(1,p+1))
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
        for i in [S[e-1]-1,..(S[e-1]-2,0)]:
          if in(!(i),S) :
            while j>0 and S[j-1]>i:
              j=j-1
            k=j+1
            while <=(k,e) and <=(S[k-1]-i,s):
              z=[i]
              z.extend([S[k-1] for k in range(j+1,k-1+1)])
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
            res.extend(p)
            res.append(-1)
        p=".0123456789abcdefghijklmnopqrstuvwxyz"
        return [p[k-1] for k in 2+res]
      
      
      def Strips(S,s):
        def apply(S,hs):
          S=ShallowCopy(S)
          for h in hs:
            [S[k-1] for k in range(h["stoppos"],h["startpos"]+1)]=Concatenation([h["start"]-h["area"]],[S[k-1] for k in range(h["stoppos"],h["startpos"]-1+1)])
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
            ss["SC"].extend(hs["SC"])
            ss["DC"].extend(hs["DC"])
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
      def GenericEntry(lambda,mu):
        n=Sum(lambda,Sum)
        if n==0 :
          return 1
        if n<LIM :
          name=code(lambda,mu)
          if name in chiCache :
            return chiCache[name]
        bp=Maximum(Concatenation(lambda))
        i=PositionProperty(lambda,lambda x: in(bp,x))
        rest=ShallowCopy(lambda)
        rest[i-1]=[rest[i][k-1] for k in range(2,len(rest[i-1])+1)]
        res=GAPMul(-prod(para[2-1])**GAPMul(i-1,n-bp),Sum(Strips(mu,bp),ximp35))
        if n<LIM :
          chiCache[name]=res
        return res
      
      
      res["irreducibles"]=map(lambda x: map(lambda y: GenericEntry(y,x),cl["classparams"]),map(lambda x: map(BetaSet,x),cl["classparams"]))
    else:
      if [q,r]==[2,2] and !("othermethod" in CHEVIE) :
        cl=ChevieData["imp"]["ClassInfo"](p,q,r)
        X=para[2-1]
        Y=para[3-1]
        Z=para[1-1]
        def GenericEntry(char,class):
          char=ci.malle[ci["charparams"].index(char)+1-1]
          if char[1-1]==1 :
            w=[Z[char[4-1]-1],X[char[2-1]-1],Y[char[3-1]-1]]
            return prod(class)
          else:
            w=GAPMul(char[2-1],GetRoot(GAPMul(X[1-1],X[2-1]),2))
            class=map(lambda i: Number(class,lambda j: i==j),range(0,3+1))
            if class[2-1]>0 :
              char=Sum([Z[k-1] for k in [char[k-1] for k in [3,4]]],lambda x: x**class[2-1])
            else:
              if class[3-1]>0 :
                char=Sum(X)
              else:
                if class[4-1]>0 :
                  char=Sum(Y)
                else:
                  char=2
            return GAPMul(w**class[1-1],char)
        
        
        res["classes"]=cl["classes"]
        res["orders"]=cl["orders"]
        res["irreducibles"]=map(lambda char: map(lambda class: GenericEntry(char,class),cl["classparams"]),ci["charparams"])
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

def ximp35(x):
  d=Delta(i-1,x,para[2-1],para[1-1])
  if d==0 :
    return d
  else:
    return GAPMul(d,GenericEntry(rest,x["remainder"]))

ChevieData["imp"]["HeckeCharTable"]=ximp34

def ximp36(p,q,r,para,root,i):
  if !(IsList(para)) :
    para=[para]
  if [q,r]==[1,2] :
    X=para[2-1]
    Y=para[1-1]
    t=PartitionTuples(2, p)[i-1]
    if Number(t,lambda x: !=(x,[]))==1 :
      p=PositionProperty(t,lambda x: !=(x,[]))
      if t[p-1]==[2] :
        return GAPMul(X[1-1]**0,[[[Y[p-1]]],[[X[1-1]]]])
      else:
        return GAPMul(X[1-1]**0,[[[Y[p-1]]],[[X[2-1]]]])
    else:
      p=Filtered(range(1,len(t)+1),lambda i: !=(t[i-1],[]))
      return GAPMul(X[1-1]**0,[[[Y[p[1-1]-1],0],[-1,Y[p[2-1]-1]]],[[X[1-1],GAPMul(X[1-1],Y[p[1-1]-1])+GAPMul(X[2-1],Y[p[2-1]-1])],[0,X[2-1]]]])
  else:
    if [p,q,r]==[3,3,3] :
      x=GAPDiv(-para[2][1-1],para[2][2-1])
      def f(x,j):
        return [[[-1,0,0],[0,0,1],[0,x,-1+x]],[[-1,0,0],[x-x**2,-1+x,j**2],[GAPMul(j,x)-GAPMul(j,x**2),GAPMul(j,x),0]],[[0,1,0],[x,-1+x,0],[0,0,-1]]]
      
      
      r=GAPMul(x**0,[[[[-1,0],[-1,x]],[[x,-x],[0,-1]],[[x,-x],[0,-1]]],[[[-1,0],[-1,x]],[[x,-x],[0,-1]],[[-1,0],[-1,x]]],[[[-1,0],[-1,x]],[[x,-x],[0,-1]],[[-1+x,1],[x,0]]],f(x,ER(3)),f(x,ER(3)**2),[[[-1]],[[-1]],[[-1]]],GAPMul(-x,f(x**-1,ER(3)**2)),GAPMul(-x,f(x**-1,ER(3))),[[[-1,0],[-1,x]],[[-1,0],[-1,x]],[[x,-x],[0,-1]]],[[[x]],[[x]],[[x]]]])
      return r[i-1]
    else:
      if [p,q,r]==[2,2,4] :
        x=GAPDiv(-para[1][1-1],para[1][2-1])
        r=[lambda x: [[[-1+x,-1,0],[-x,0,0],[x-x**2,-1+x,-1]],[[0,1,0],[x,-1+x,0],[0,0,-1]],[[-1,0,0],[0,0,1],[0,x,-1+x]],[[0,1,0],[x,-1+x,0],[0,0,-1]]],lambda x: [[[0,1,0],[x,-1+x,0],[0,0,-1]],[[-1+x,-1,0],[-x,0,0],[x-x**2,-1+x,-1]],[[-1,0,0],[0,0,1],[0,x,-1+x]],[[0,1,0],[x,-1+x,0],[0,0,-1]]],lambda x: [[[-1,0,0,0],[0,-1+x,-1,0],[0,-x,0,0],[0,0,0,-1]],[[-1,1-x,1-x,0],[0,0,1,0],[0,x,-1+x,0],[0,-1+x,-1+x,-1]],[[-1+x,-x,0,0],[-1,0,0,0],[0,0,-1,0],[0,0,0,-1]],[[0,0,0,1],[0,-1,0,0],[0,0,-1,0],[x,0,0,-1+x]]],lambda x: [[[-1]],[[-1]],[[-1]],[[-1]]],lambda x: [[[x,1-x,-1+x,-x+x**2,x-x**2,0],[0,-1+x,0,0,-x,x-x**2],[0,0,-1+x,-x,0,x-x**2],[0,0,-1,0,0,-1+x],[0,-1,0,0,0,-1+x],[0,0,0,0,0,-1]],[[x,0,0,0,0,0],[0,0,0,0,x,0],[0,0,0,x,0,0],[0,0,1,-1+x,0,0],[0,1,0,0,-1+x,0],[0,0,0,0,0,-1]],[[0,0,x,0,0,0],[0,-1,0,0,0,0],[1,0,-1+x,0,0,0],[0,0,0,x,0,0],[0,0,0,0,0,x],[0,0,0,0,1,-1+x]],[[-1,0,0,0,0,0],[0,-1+x,1,0,0,0],[0,x,0,0,0,0],[0,0,0,0,x,0],[0,0,0,1,-1+x,0],[0,0,0,0,0,x]]],lambda x: [[[-1+x,0,-1,0,0,0,0,0],[0,0,0,0,1,0,0,0],[-x,0,0,0,0,0,0,0],[0,0,0,0,0,0,1,0],[0,x,0,0,-1+x,0,0,0],[0,0,0,0,0,-1+x,0,x],[0,0,0,x,0,0,-1+x,0],[0,0,0,0,0,1,0,0]],[[0,0,1,0,0,0,0,0],[0,0,-1+x,0,1,0,GAPDiv(-1+x,x),0],[x,0,-1+x,0,0,0,0,0],[0,0,0,-1+x,0,0,-1,0],[x-x**2,x,0,-1+x,-1+x,0,GAPDiv(1-GAPMul(2,x)+x**2,x),0],[-x+x**2,0,0,-x+x**2,0,-1+x,1-x,x],[0,0,0,-x,0,0,0,0],[0,0,1-x,x-x**2,0,1,-1+x,0]],[[0,1,0,0,0,0,0,0],[x,-1+x,0,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,x,-1+x,0,0,0,0],[0,0,0,0,-1+x,0,-1,0],[0,0,0,0,0,x,0,0],[0,0,0,0,-x,0,0,0],[0,0,0,0,0,-x,0,-1]],[[-1,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,1],[0,0,-1,-x,0,0,0,0],[0,0,0,x,0,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,x,-1+x,0,0],[0,0,0,0,0,0,x,0],[0,x,0,0,0,0,0,-1+x]]],lambda x: [[[-1,-1,0],[0,x,0],[0,1,-1]],[[-1,-1,0],[0,x,0],[0,1,-1]],[[-1+x,x,0],[1,0,0],[0,0,-1]],[[0,0,1],[0,-1,0],[x,0,-1+x]]],1,2,lambda x: [[[x,0],[-1,-1]],[[x,0],[-1,-1]],[[0,1],[x,-1+x]],[[x,0],[-1,-1]]],3,7,4]
        if IsInt(r[i-1]) :
          return GAPMul(-x,r[r[i-1]-1](x**-1))
        else:
          return GAPMul(r[i-1](x),x**0)
      else:
        if [p,q,r]==[3,3,4] :
          x=GAPDiv(-para[2][1-1],para[2][2-1])
          def m(i):
            f1=lambda x: GAPMul(x**0,[[[x,-1,0,0,0,0,0,0,0,0,1-x-x**2+x**3,0],[0,-1,0,0,0,0,0,0,0,0,0,0],[0,0,-1+x,0,x,0,-x,0,0,0,x-x**2,0],[0,0,0,-1+x,0,0,-x,0,0,0,x-x**2,0],[0,0,1,-1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-x],[0,0,0,-1,0,0,0,0,0,0,-1+x,0],[0,0,0,0,0,0,0,-1,0,0,0,0],[0,0,0,0,0,0,0,0,-1+x,1,-1+x,0],[0,0,0,0,0,0,0,0,x,0,-1+x,0],[0,0,0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,-1,0,0,0,0,0,-1+x]],[[0,x,0,0,0,0,0,0,0,0,0,0],[1,-1+x,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,x,0,0,0,0,0],[0,0,0,0,x,0,0,1-x,0,0,0,0],[0,0,0,1,-1+x,0,0,1-x,0,0,0,0],[0,0,0,0,0,0,0,1-x,-1+x,1,0,1-x+x**2],[0,0,1,0,0,0,-1+x,0,0,0,0,0],[0,0,0,0,0,0,0,-1,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,x],[0,0,0,0,0,x,0,x-x**2,-x,-1+x,0,x-x**2],[0,0,0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,0,0,1,0,0,-1+x]],[[0,-1+GAPMul(2,x)-x**2,1-x,x,-x+x**2,0,0,0,-1+GAPMul(2,x)-x**2,0,0,0],[0,-1+x,1,0,0,0,0,0,-1+x,0,0,0],[0,x,0,0,0,0,0,0,-1+x,0,0,0],[1,-1+x,0,-1+x,0,0,1-x,0,0,0,0,0],[0,0,0,0,-1+x,0,1,0,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,0,0],[0,0,0,0,x,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,-1+x,0,0,0,-1],[0,0,0,0,0,0,0,0,-1,0,0,0],[0,0,0,0,0,0,0,0,0,0,x,0],[0,0,0,0,0,0,0,0,0,1,-1+x,0],[0,0,0,0,0,0,0,-x,0,0,0,0]],[[-1,0,0,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,x,0,0,0],[0,0,0,0,0,x,0,0,0,0,1-x,x-x**2],[0,0,0,0,0,0,0,0,0,1,0,x],[0,0,0,1,0,-1+x,-1+x,0,0,0,1-x,0],[0,0,0,0,0,0,0,0,0,0,0,x],[0,0,0,0,0,0,0,x,0,0,-1,0],[0,0,1,0,0,0,0,0,-1+x,0,0,0],[0,0,0,0,x,0,-x,0,0,-1+x,0,0],[0,0,0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,1,0,0,0,0,-1+x]]])
            def f2(x,j):
              return [[[-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,x,x,x]],[[-1,0,0,0],[0,-1,0,0],[0,-j**2,x,1],[0,0,0,-1]],[[-1,0,0,0],[x,x,GAPMul(-j,x),1],[0,0,-1,0],[0,0,0,-1]],[[x,1,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,-1]]]
            
            
            f3=lambda x: [[[x,-1,0,-x,0,0],[0,-1,0,0,0,0],[0,0,-1+x,0,1,x],[0,0,0,-1,0,0],[0,0,x,0,0,x],[0,0,0,0,0,-1]],[[-1,0,0,0,0,0],[0,0,0,0,x,1],[0,0,-1,0,0,0],[-1,0,-1,x,0,-1+x],[0,1,0,0,-1+x,1],[0,0,0,0,0,-1]],[[0,x,1,-1,-1,0],[1,-1+x,1,-1,-1,0],[0,0,-1,0,0,0],[0,0,0,-1,0,0],[0,0,0,0,-1,0],[0,0,0,1,1,x]],[[x,-1,0,0,1,x],[0,-1,0,0,0,0],[0,0,-1,0,0,0],[0,0,-1,x,1,x],[0,0,0,0,-1,0],[0,0,0,0,0,-1]]]
            f5=lambda x: [[[-1]],[[-1]],[[-1]],[[-1]]]
            def f7(x,j):
              return [[[-1,0,0,0,0,0],[x,x,0,0,0,0],[x,0,x,0,0,0],[0,0,0,-1,0,0],[0,0,0,0,-1,0],[0,0,0,GAPMul(-j,x**2),x,x]],[[x,1,0,0,0,0],[0,-1,0,0,0,0],[0,j**2,x,0,0,0],[0,0,0,-1,0,0],[0,0,0,x,x,1],[0,0,0,0,0,-1]],[[x,0,1,0,1,0],[0,x,GAPMul(j,x),0,0,1],[0,0,-1,0,0,0],[0,0,0,x,1,GAPMul(-j**2,x**-1)],[0,0,0,0,-1,0],[0,0,0,0,0,-1]],[[-1,0,0,0,0,0],[0,-1,0,0,0,0],[0,0,-1,0,0,0],[0,0,0,x,0,0],[x,0,0,0,x,0],[0,x,0,0,0,x]]]
            
            
            def f8(x,j):
              return [[[-1,0,0,0,0,0,0,0],[1,x,0,0,0,0,1,0],[1,0,x,0,0,0,0,0],[0,0,0,-1,0,0,0,0],[0,0,0,0,-1,0,0,0],[0,0,0,GAPMul(-j,x),x,x,0,0],[0,0,0,0,0,0,-1,0],[0,0,0,0,GAPMul(j**2-j,x),0,1,x]],[[x,x,0,0,0,0,-j**2,0],[0,-1,0,0,0,0,0,0],[0,j,x,0,0,0,0,0],[0,0,0,-1,0,0,0,0],[0,0,0,1,x,1,0,0],[0,0,0,0,0,-1,0,0],[0,0,0,0,0,0,-1,0],[0,0,0,0,0,GAPMul(-2,j**2)-j,1,x]],[[x,0,x,0,x,0,0,0],[0,x,GAPMul(j**2,x),0,0,1,0,0],[0,0,-1,0,0,0,0,0],[0,0,0,x,x,-j**2,0,0],[0,0,0,0,-1,0,0,0],[0,0,0,0,0,-1,0,0],[0,0,0,0,0,0,x,x],[0,0,0,0,0,0,0,-1]],[[-1,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0],[0,0,0,x,0,0,-j**2,0],[1,0,0,0,x,0,0,0],[0,x,0,0,0,x,0,0],[0,0,0,0,0,0,-1,0],[0,0,GAPMul(j**2-j,x),0,0,0,1,x]]]
            
            
            f11=lambda x: [[[x,1,0],[0,-1,0],[0,0,-1]],[[x,1,0],[0,-1,0],[0,0,-1]],[[-1,0,0],[x,x,1],[0,0,-1]],[[-1,0,0],[0,-1,0],[0,x,x]]]
            f13=lambda x: [[[-1,0],[x,x]],[[-1,0],[x,x]],[[x,1],[0,-1]],[[-1,0],[x,x]]]
            r=[f1(x),f2(x,ER(3)),f3(x),f2(x,ER(3)**2),f5(x),GAPMul(-x,f1(x**-1)),f7(x,ER(3)),f8(x,ER(3)),f8(x,ER(3)**2),GAPMul(-x,f7(x**-1,ER(3))),f11(x),GAPMul(-x,f3(x**-1)),f13(x),GAPMul(-x,f2(x**-1,ER(3)**2)),GAPMul(-x,f2(x**-1,ER(3))),GAPMul(-x,f11(x**-1)),GAPMul(-x,f5(x**-1))]
            return GAPMul(x**0,r[i-1])
          
          
          return m(i)
        else:
          if [p,q,r]==[3,3,5] :
            x=GAPDiv(-para[2][1-1],para[2][2-1])
            def m(i):
              def f1(x):
                return GAPMul(x**0,[[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[ER(3)**2-GAPMul(ER(3)**2,x),ER(3)**2-GAPMul(ER(3)**2,x),0,0,0,ER(3)**2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,GAPMul(-ER(3)**2,x)+GAPMul(ER(3)**2,x**2),0,ER(3)**2-GAPMul(ER(3)**2,x),0,0,0,ER(3)**2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-ER(3)**2-GAPMul(ER(-3),x)+GAPMul(ER(3),x**2),ER(3)-GAPMul(ER(3),x)+GAPMul(ER(3),x**2),0,0,0,ER(3)-GAPMul(ER(3),x),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,GAPMul(ER(3)**2,x),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,ER(3),-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,GAPMul(ER(3)**2,x)+GAPMul(ER(-3),x**2)-GAPMul(ER(3),x**3),0,ER(3)-GAPMul(ER(3),x)+GAPMul(ER(3),x**2),0,0,0,ER(3)-GAPMul(ER(3),x),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,ER(3)**2-GAPMul(ER(3)**2,x),0,0,0,0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,x),0,0,0,0,ER(3)**2,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,x),0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3),x),0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-ER(3)**2-GAPMul(ER(-3),x)+GAPMul(ER(3),x**2),0,0,0,0,0,0,0,ER(3)-GAPMul(ER(3),x)+GAPMul(ER(3),x**2),0,0,0,0,ER(3)-GAPMul(ER(3),x),0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3),x),0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,x),0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2,0,-1+x,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0,-1+x,0,0,0,0,0,0,0,0,0,0],[0,0,ER(3)**2-GAPMul(ER(3)**2,x)-GAPMul(ER(3)**2,x**2)+GAPMul(ER(3)**2,x**3),0,ER(3)-GAPMul(ER(-3),x)-GAPMul(ER(3)**2,x**2),0,0,0,-ER(3)**2+GAPMul(ER(3)**2,x),0,-ER(3)+GAPMul(ER(-3),x),0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,x),0,0,0,0,-1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2,0,0,0,0,0,0],[ER(3)-GAPMul(ER(3),x),0,0,GAPMul(2,ER(3))-GAPMul(ER(3),x**-1)-GAPMul(ER(3),x),0,0,0,0,0,0,-ER(3)+GAPMul(2,ER(3))-GAPMul(ER(3),x**2),0,0,0,0,-ER(3)+GAPMul(ER(3),x),0,0,0,ER(3)-GAPMul(ER(3),x),0,0,ER(3)-GAPMul(ER(3),x),0,ER(3),0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3),x),0,-1+x,0,0,0,0,0,0],[-ER(3)+GAPMul(ER(-3),x),0,x-GAPMul(2,x**2)+x**3,GAPDiv(1-GAPMul(3,ER(-3)),2)+GAPMul(ER(3),x**-1),0,0,0,0,0,0,ER(3)-GAPMul(2,ER(3))+GAPMul(2,ER(3))-GAPMul(ER(3),x**3),0,0,0,0,ER(3)-GAPMul(2,ER(3))+GAPMul(ER(3),x**2),0,-ER(3)+GAPMul(ER(3),x),0,-ER(3)+GAPMul(2,ER(3))-GAPMul(ER(3),x**2),0,0,ER(3)**2-GAPMul(ER(3)**2,x)+GAPMul(ER(3)**2,x**2),0,ER(3)**2-GAPMul(ER(3)**2,x),0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-1+x,0,-ER(3)**2+GAPMul(ER(3)**2,x),0,0,0,-1+x**-1,0,ER(3)**2-GAPMul(ER(3)**2,x),0,0,0,0,0,-1,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,x),0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3),-1+x,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+x,GAPMul(ER(3),x)],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2,0]],[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,x,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-x,0,0,0,0,0,-1+x,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1+x**-1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,x,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-GAPMul(2,x)+x**2,1-GAPMul(2,x)+x**2,1-GAPMul(2,x)+x**2,1-x,2-x**-1-x,1-x,0,0,1-x**-1,-1+x,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1-x,0,0,0,0,0,0,0,-1+x,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-GAPMul(2,x)+x**2,-x+x**2,1-GAPMul(2,x)+x**2,1-x,1-x,0,0,0,0,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1-x,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[2-x**-1-x,1-x,0,0,0,1-x,-1+x,-1+x,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,x-x**2,-1+GAPMul(2,x)-x**2,0,0,0,x-x**2,0,0,0,0,0,-1+x,0,0,-1+x,0,x,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+x,0,x,0,0,0,0,0,0,0,0,0,0],[2-x**-1-x,0,1-x,0,0,0,-1+x,0,0,0,0,1-x,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1-x**-1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0],[-1+GAPMul(3,x)-GAPMul(3,x**2)+x**3,-1+GAPMul(3,x)-GAPMul(3,x**2)+x**3,x-GAPMul(2,x**2)+x**3,2-x**-1-x,-2+x**-1-x**2,-1+GAPMul(2,x)-x**2,0,1-GAPMul(2,x)+x**2,-2+x**-1,0,0,-1+x,-1+x,2-x**-1-x,1-x,0,0,0,1-x,0,0,-1+x,0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],[-3+x**-1-x**2,-1+GAPMul(2,x)-x**2,0,1-GAPMul(2,x)+x**2,0,-1+GAPMul(2,x)-x**2,1-GAPMul(2,x)+x**2,1-GAPMul(2,x)+x**2,-1+x,x-x**2,0,-1+GAPMul(2,x)-x**2,0,1-x,1-x,0,-1+x,0,0,0,0,x,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,-1+x,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0],[1-GAPMul(3,x)+GAPMul(3,x**2)-x**3,1-GAPMul(2,x)+GAPMul(2,x**2)-x**3,1-GAPMul(3,x)+GAPMul(3,x**2)-x**3,1-GAPMul(2,x)+x**2,3-x**-1-GAPMul(3,x)+x**2,1-GAPMul(2,x)+x**2,0,0,2-x**-1-x,0,-1+GAPMul(2,x)-x**2,0,1-x,0,0,0,0,-1+x,0,0,1-x,0,1-GAPMul(2,x)+x**2,0,1-x,0,-1+x,x,0,0],[0,0,0,2-x**-1-x,0,-1+x,0,0,0,-1+x,0,0,0,0,0,-2+x**-1,0,2-x**-1-x,0,1-x,-1+x**-1,0,-1+x,0,0,0,1,0,0,0],[-1+GAPMul(2,x)-x**2,0,0,0,0,0,1-GAPMul(2,x)+x**2,-x+x**2,0,0,0,0,0,0,0,1-x,-1+x,0,0,0,0,0,1-x,0,0,1-x,0,0,0,x],[0,0,1-GAPMul(2,x)+x**2,3-x**-1-GAPMul(3,x)+x**2,0,0,1-x,1-GAPMul(2,x)+x**2,0,0,-1+x,0,0,2-x**-1-x,0,0,0,0,1-x,0,0,0,2-x**-1-x,0,1-x**-1,1-x,0,0,1,-1+x]],[[0,-x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,x,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,x,0,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1+x,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,-1+x,0,0,0,0,0,0],[0,x-x**2,0,0,0,0,0,0,0,1-GAPMul(2,x)+x**2,-1+GAPMul(2,x)-x**2,0,0,0,0,0,0,-1+x,0,0,1-x,0,1-GAPMul(2,x)+x**2,0,0,0,-1+x,x,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+x,0,0,0,0,0,1-x,0,0,0,0,0,0,x],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,-1+x,0,0,0],[1-x,0,0,2-x**-1-x,0,0,0,0,0,0,-1+GAPMul(2,x)-x**2,0,0,0,0,-1+x,0,0,0,1-x,0,0,1-x,0,1,0,0,-1+x,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,0,0,0,0,0,1-x,0,0,0,0,0,1-x**-1,0,0,0,0,0,0,0,1,0,0,0,-1+x]],[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,x,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-x,0,0,0,0,0,-1+x,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,x,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-x,1-x,0,0,0,1,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,x-x**2,0,0,0,0,0,0,0,-x+x**2,0,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0,0,0,0,0],[0,0,-x+x**2,0,1-x,0,0,0,1,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+x,0,x,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0],[0,0,1-x,0,0,0,0,0,0,0,1-x,0,0,0,0,1,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1+x,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,-1+x,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1+x]],[[0,0,-x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,x,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,1-x,0,0,0,0,0,0,0,1-x,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,x,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1-x,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],[1-x,0,0,0,0,0,-1+x,x,0,0,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0,0,0],[1-x,0,0,2-x**-1-x,0,0,0,0,0,0,-1+GAPMul(2,x)-x**2,0,0,0,0,-1+x,0,0,0,1-x,0,0,1-x,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0,0],[0,-x+x**2,0,0,0,0,1-x,0,0,0,0,0,0,0,x,0,0,0,0,0,-1+x,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,-1+x,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0],[0,0,x-x**2,-1+GAPMul(2,x)-x**2,0,0,0,x-x**2,0,0,0,0,0,-1+x,0,0,-1+x,0,x,0,0,0,0,0,-1+x,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x,0,0,0,0,-1+x,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,-1+x,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1]]])
              
              
              def f2(q):
                return GAPMul(q**0,[[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-ER(3)+GAPMul(ER(-3),q),ER(3)**2-GAPMul(ER(3)**2,q),0,0,-ER(3)+GAPMul(ER(3),q),0,ER(-3)-GAPMul(ER(3),q**-1)+GAPMul(ER(3)**2,q),-ER(3)**2+GAPMul(ER(3)**2,q**-1),0,0,GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,ER(3)**2,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[GAPMul(ER(3),q)-GAPMul(ER(3),q**2),GAPMul(ER(3),q),0,0,0,0,ER(3)-GAPMul(ER(3),q),ER(3)-GAPMul(ER(3),q),0,0,ER(3)-GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0,0,0,0,0,0,0],[0,0,0,0,GAPMul(ER(3),q),0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,ER(3)**2-GAPMul(ER(3)**2,q),0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,0,0,ER(3)**2,0,0,0,0],[0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,-1+q,0,0,0,0,0,0,0],[-ER(3)+GAPMul(ER(3),q),0,0,0,0,ER(3)-GAPMul(ER(3),q),0,0,0,0,0,0,0,-1+q,-ER(3),0,0,0,0,0],[1-q,0,0,0,0,-1+q,0,0,0,0,0,0,0,GAPMul(-ER(3)**2,q),0,0,0,0,0,0],[0,0,0,0,0,0,-ER(3)**2-GAPMul(ER(-3),q)+GAPMul(ER(3),q**2),0,-ER(3)**2-GAPMul(ER(-3),q)+GAPMul(ER(3),q**2),0,0,ER(3)-GAPMul(ER(3),q)+GAPMul(ER(3),q**2),0,0,0,ER(3)-GAPMul(ER(3),q),0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),-1+q,0,0],[GAPMul(ER(3),q)-GAPMul(2,ER(3))+GAPMul(ER(3),q**3),GAPMul(ER(3),q)-GAPMul(ER(3),q**2),0,0,0,0,ER(3)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),ER(3)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),0,0,ER(3)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),0,-ER(3)+GAPMul(ER(3),q),0,0,0,0,-ER(3)+GAPMul(ER(3),q),0,GAPMul(ER(3)**2,q)],[0,0,0,0,0,0,0,-ER(3)+GAPMul(ER(3),q),0,ER(3)-GAPMul(ER(3),q),0,0,0,0,0,0,ER(3)-GAPMul(ER(3),q),0,ER(3),-1+q]],[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,-1+q,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,q**3-q**4,1-q,-1+q,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0],[1-GAPMul(2,q)+q**2,0,q**3-q**4,0,0,0,0,0,1-q,-1+q,0,0,1,1-q,0,0,0,0,0,0],[0,0,q**3-q**4,1-q,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],[1-GAPMul(2,q)+q**2,0,q**3-q**4,0,0,0,0,0,1-q,q,0,0,0,0,1-q,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,-q,-1+q,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,-1+q,0,0,0,0],[0,0,0,1-q,0,-1+q,0,0,0,0,0,2-q**-1-q,0,0,0,1-q**-1,-1+q,1,0,0],[0,0,0,1-q,0,-1+q,0,0,0,0,0,1-q,0,0,0,0,q,0,0,0],[-1+GAPMul(2,q)-q**2,0,0,0,q-q**2,1-GAPMul(2,q)+q**2,0,0,0,0,0,-1+q,0,0,-1+q,0,0,0,-1+q,q],[2-q**-1-q,0,-q**3+GAPMul(2,q**4)-q**5,-1+GAPMul(2,q)-q**2,1-GAPMul(2,q)+q**2,-2+q**-1,0,0,0,0,-1+q,0,0,-1+q,2-q**-1-q,-1+q**-1,0,0,1,0]],[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-q,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,-q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,1,-1+q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0],[1-q,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0],[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1+q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1+q,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1+q]],[[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-q**-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-q**3,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],[-q,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,q,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0,0],[0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0],[-1+q,0,0,0,0,1-q,0,0,0,0,0,-1+q,0,0,-1,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0],[0,0,0,0,0,0,1-q,0,1-q,0,0,1-q,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1-q,0,1-q,0,0,-q,0,0,0,0,0,0,0,0],[1-q,0,0,0,0,-1+q,0,0,0,0,0,0,0,q,1-q,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,-1+q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,-1+q,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1]],[[-1+q,0,-q**3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-GAPMul(2,q)+q**2,0,q**3-q**4,0,0,0,0,0,1-q,-1+q,0,0,1,1-q,0,0,0,0,0,0],[-q**-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-q,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,1,0,0,0,0,0],[0,0,0,1,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,-q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,-1+q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,1,0,-1+q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0],[q-q**2,q,0,0,0,0,1-q,1-q,0,0,1-q,0,-1+q,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,-1,0,0,-1+q,0,0,0,0,0,0],[0,0,q**3-q**4,1-q,q,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1]]])
              
              
              def f3(q,j):
                return GAPMul(q**0,[[[-1,0,0,0,0],[0,-1,0,0,0],[1,0,0,0,-1],[1+GAPMul(j,q),0,1+GAPMul(j,q),-1,-1-GAPMul(j,q)],[-q,0,-q,0,-1+q]],[[-1,0,0,0,0],[0,-1,0,0,0],[1,0,0,-1,0],[-q,0,-q,-1+q,0],[-j,0,-j,j,-1]],[[-1,0,-1,0,0],[0,-1,1,0,0],[0,0,q,0,0],[0,0,0,-1,0],[0,0,0,0,-1]],[[q,0,0,0,0],[-1,-1,0,0,0],[-q,0,-1,0,0],[1,0,0,-1,0],[1,0,0,0,-1]],[[0,1,0,0,0],[q,-1+q,0,0,0],[0,0,-1,0,0],[1,1,0,-1,0],[1,1,0,0,-1]]])
              
              
              def f4(q,j):
                return GAPMul(q**0,[[[-1+q,0,0,q,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,1],[0,0,-1,0,0,0,0,0,0,0],[1,0,0,0,0,0,0,0,0,0],[0,0,GAPMul(j**2,q)+GAPMul(-j**2+j,q**2)-GAPMul(j,q**3),0,j-GAPMul(j,q),0,-j+GAPMul(j,q)-GAPMul(j,q**2),0,0,0],[0,0,0,0,0,-1,0,0,0,0],[0,0,GAPMul(j**2,q)-GAPMul(j**2,q**2),0,-j**2,0,j**2-GAPMul(j**2,q),0,0,0],[-q+q**2,0,GAPMul(j**2,q)-GAPMul(2,j**2)+GAPMul(j**2,q**3),-q+q**2,-j**2+GAPMul(j**2,q),0,-j+GAPMul(-j**2+j,q),-1,0,0],[0,q-q**2,GAPMul(j**2,q)-GAPMul(2,j**2)+GAPMul(j**2,q**3),0,-j**2+GAPMul(j**2,q),0,-j+GAPMul(-j**2+j,q),0,-1,q-q**2],[0,q,0,0,0,0,0,0,0,-1+q]],[[0,0,GAPMul(j**2,q)-GAPMul(j**2,q**2),GAPMul(j**2,q),0,0,0,0,0,0],[0,-1+q,-q+q**2,0,0,0,0,0,0,j],[0,0,-1,0,0,0,0,0,0,0],[j,0,q-q**2,-1+q,0,0,0,0,0,0],[0,0,0,0,-1+q,0,-q,0,0,0],[q-q**2,GAPMul(-j**2,q)+GAPMul(j**2,q**2),GAPMul(j**2,q)-GAPMul(j**2,q**2)-GAPMul(j**2,q**3)+GAPMul(j**2,q**4),GAPMul(j**2,q**2)-GAPMul(j**2,q**3),0,-1,0,0,0,-1+q],[0,0,0,0,-1,0,0,0,0,0],[0,0,0,0,0,0,0,-1,0,0],[0,0,0,0,0,0,0,0,-1,0],[0,GAPMul(j**2,q),GAPMul(-j**2,q)+GAPMul(j**2,q**2),0,0,0,0,0,0,0]],[[-1,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,q,-1+q,0,0,0,0,0,0],[0,0,0,0,0,0,0,-1,0,0],[0,0,0,0,0,-1+q,0,0,0,q],[0,0,0,0,0,0,-1,0,0,0],[0,0,0,0,-q,0,0,-1+q,0,0],[0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,1,0,0,0,0]],[[0,0,0,0,0,0,0,0,0,1],[0,-1+q,0,q,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,0,0,-1,0,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0],[0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,0,-q,-1+q,0],[q,0,0,0,0,0,0,0,0,-1+q]],[[-1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0],[0,0,-1,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0],[0,-q+q**2,-q**2+q**3,0,-1+q,0,0,0,0,GAPMul(j,q)],[0,0,0,GAPMul(-j**2,q)+GAPMul(j**2,q**2),0,0,j**2-GAPMul(j**2,q),-j**2,0,0],[0,-q,0,0,0,0,-1+q,0,0,0],[0,-q+q**2,0,q**2-q**3,0,GAPMul(-j,q),0,-1+q,0,0],[0,0,0,0,0,0,0,0,-1,0],[0,0,GAPMul(-j**2,q)+GAPMul(j**2,q**2),0,j**2,0,-j**2+GAPMul(j**2,q),0,0,0]]])
              
              
              def f11(q,j):
                return GAPMul(q**0,[[[q,0,0,0,0,0,0,0,0,0],[-j**2+GAPMul(j**2,q),j**2-GAPMul(j**2,q),j**2,0,0,0,0,0,0,0],[-j+GAPMul(-j**2+j,q),j-GAPMul(j,q)+GAPMul(j,q**2),j-GAPMul(j,q),0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0],[0,0,0,q,0,-1+q,0,0,0,0],[0,0,0,0,q,0,-1+q,0,0,0],[0,0,0,0,0,0,0,-1,0,0],[0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,0,0,0,-1]],[[q,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0],[0,q,-1+q,0,0,0,0,0,0,0],[1-q,0,0,-1+q,0,j,0,0,0,0],[1-q,0,0,0,-1+q,0,j,0,0,0],[GAPMul(-j**2,q)+GAPMul(j**2,q**2),0,0,GAPMul(j**2,q),0,0,0,0,0,0],[GAPMul(-j**2,q)+GAPMul(j**2,q**2),0,0,0,GAPMul(j**2,q),0,0,0,0,0],[-j**2+GAPMul(2,j**2)-GAPMul(j**2,q**2),j**2-GAPMul(j**2,q),j**2-GAPMul(j**2,q),GAPMul(-j**2,q)+GAPMul(j**2,q**2),0,-1+q,0,-1,0,0],[0,0,0,GAPMul(j**2,q)-GAPMul(j**2,q**2),GAPMul(-j**2,q)+GAPMul(j**2,q**2),1-q,-1+q,0,-1,0],[-j**2+GAPMul(2,j**2)-GAPMul(j**2,q**2),j**2-GAPMul(j**2,q),j**2-GAPMul(j**2,q),0,GAPMul(-j**2,q)+GAPMul(j**2,q**2),0,-1+q,0,0,-1]],[[0,-1,0,0,0,0,0,0,0,0],[-q,-1+q,0,0,0,0,0,0,0,0],[0,0,q,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0],[0,0,0,0,-1,0,0,0,0,0],[0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0,1],[0,0,0,0,0,q,0,-1+q,0,0],[0,0,0,0,0,0,0,0,-1,0],[0,0,0,0,0,0,q,0,0,-1+q]],[[-1,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0],[-1+q,0,0,1-q,0,-j,0,0,0,0],[0,-q,0,-1+q,0,0,0,0,0,0],[0,0,0,0,-1,0,0,0,0,0],[GAPMul(j**2,q)-GAPMul(j**2,q**2),GAPMul(-j**2,q)+GAPMul(j**2,q**2),GAPMul(-j**2,q),0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0],[0,0,0,0,0,0,0,q,0,0],[0,0,0,0,0,0,0,0,-1+q,-q],[0,0,0,0,0,0,0,0,-1,0]],[[-1,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0],[0,0,0,0,-1,0,0,0,0,0],[0,0,0,-q,-1+q,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0],[0,0,0,0,0,-q,-1+q,0,0,0],[0,0,0,0,0,0,0,0,0,-1],[0,0,0,0,0,0,0,0,q,0],[0,0,0,0,0,0,0,-q,0,-1+q]]])
              
              
              def f12(q,j):
                GAPMul(returnq**0,[[[0,0,0,0,-j**2,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,GAPMul(-j**2,q),0,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-j**2,0,0,0,0,0,0,0,0],[GAPMul(-j,q),0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0],[0,-j,0,0,0,-1+q,0,0,0,0,0,0,0,0,0],[0,0,0,GAPMul(-j,q),0,0,-1+q,0,0,0,0,0,0,0,0],[-j**2+GAPMul(j**2,q),0,0,0,j-GAPMul(j,q),0,0,-1,0,0,0,0,0,0,0],[0,0,0,j**2-GAPMul(j**2,q),0,0,-j+GAPMul(j,q),0,-1,0,0,0,0,0,0],[0,GAPMul(2,j**2)-GAPMul(j**2,q**-1)-GAPMul(j**2,q),0,0,0,j-GAPMul(2,j)+GAPMul(j,q**2),0,0,0,-1,0,0,0,0,0],[0,-1+GAPMul(2,q)-q**2,-q,0,0,1-GAPMul(2,q)+q**2,0,0,0,q,q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0],[j-GAPMul(j,q),0,0,-q+q**2,j**2+GAPMul(2,j)-GAPMul(j,q),0,1-q,-q,-q,0,0,0,q,0,0],[0,GAPMul(-2,j**2)+GAPMul(j**2,q**-1),0,0,0,-j+GAPMul(2,j)-GAPMul(j,q**2),0,0,0,0,0,0,0,-1,0],[0,1-GAPMul(2,q)+q**2,0,0,0,-1+GAPMul(2,q)-q**2,0,0,0,0,0,q,0,q,q]],[[-1+q,0,0,0,-1,0,0,0,0,0,0,0,0,0,0],[0,-1+q,0,0,0,-q,0,0,0,0,0,0,0,0,0],[-1+q,0,-1,0,-1+q**-1,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1+q,0,0,-1,0,0,0,0,0,0,0,0],[-q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0],[0,0,0,-q+q**2,0,0,1-q,0,-1,0,0,0,0,0,0],[0,-q+q**2,0,0,0,q-q**2,0,0,0,-1,0,0,0,0,0],[0,-1+GAPMul(2,q)-q**2,-q,1-q,1-q,-q+q**2,1-q,0,0,q,q,0,0,0,0],[0,0,0,-1+q,0,0,-1+q**-1,0,0,0,0,-1,0,0,0],[-j**2+GAPMul(j**2,q),1-q**-1,0,1-GAPMul(2,q)+q**2,-j**2+GAPMul(j**2,q),-1+q,1-q,-q,-q,0,0,0,q,0,0],[0,q-q**2,0,0,0,-q+q**2,0,0,0,0,0,0,0,-1,0],[1-q**-1,2+q**-2-GAPMul(2,q**-1)-GAPMul(2,q)+q**2,0,2-q**-1-q,1-q**-1,-2+q**-1-q**2,1-q**-1,0,0,0,0,q,0,q,q]],[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1+q,0,-1,0,0,0,0,0,0,0,0,0,0],[q,-1,0,q,0,0,0,0,0,0,0,0,0,0,0],[0,0,-q,0,0,0,0,0,0,0,0,0,0,0,0],[0,1-GAPMul(2,q)+q**2,0,0,0,-1+GAPMul(2,q)-q**2,0,0,0,-q,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,q,0,0,0,0],[0,0,-j**2-GAPMul(2,j)-q**-1+GAPMul(j,q),0,-j+q**-2,0,0,-1,0,0,0,0,0,0,0],[-j+GAPMul(j**2+GAPMul(2,j),q),1-q,0,0,-j**2-GAPMul(2,j)-q**-1+GAPMul(j,q),0,-1+q,q,q,0,1-q,0,-1,0,0],[0,3-q**-1-GAPMul(4,q)+GAPMul(3,q**2)-q**3,0,0,0,-2+GAPMul(3,q)-GAPMul(3,q**2)+q**3,0,0,0,-q+q**2,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0],[0,-4+q**-1-GAPMul(4,q**2)+q**3,0,0,0,2-GAPMul(5,q)+GAPMul(4,q**2)-q**3,0,0,0,-1+GAPMul(2,q)-q**2,0,0,0,-1,0],[0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1+q]],[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,q**2,0,0,0,0,0,0,0,0,0,0,0],[-j+GAPMul(j,q),0,0,0,-j**2-GAPMul(2,j)-q**-1+GAPMul(j,q),0,0,q,0,0,0,0,0,0,0],[0,q**-1,0,-1+q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,q,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,-1+q,0,0,0,0,0,0,0,0],[-j+GAPMul(j,q),0,1,0,-j**2-GAPMul(2,j)-q**-1+GAPMul(j,q),0,0,-1+q,0,0,0,0,0,0,0],[0,-1+q,0,0,0,1-q,0,0,-1+q,-1,0,0,0,0,0],[0,0,0,-q**2+q**3,0,0,q-q**2,0,-q,0,0,0,0,0,0],[-j+GAPMul(j,q),0,0,q-q**2,-j**2-GAPMul(2,j)-q**-1+GAPMul(j,q),0,-1+q,q,q,0,0,0,-1,0,0],[0,2-q**-1-q,0,0,0,-2+q**-1,0,0,0,0,0,-1+q,0,-1,0],[0,1-GAPMul(2,q)+q**2,q,0,0,-1+GAPMul(2,q)-q**2,0,0,0,-q,-q,0,-1+q,0,0],[0,0,0,-q+GAPMul(2,q**2)-q**3,0,0,1-GAPMul(2,q)+q**2,0,0,0,0,-q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1]],[[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],[q,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,-1+GAPMul(2,q)-q**2,0,0,0,1-GAPMul(2,q)+q**2,0,0,0,q,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,q,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,-1+q,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0],[0,GAPMul(2,j**2)+j-GAPMul(j**2,q**-1)+q,0,0,0,-j**2+GAPMul(j**2,q),0,0,0,0,0,0,0,1,0],[0,0,0,1-q,0,0,1-q**-1,0,0,0,0,1,0,0,0],[1-GAPMul(2,q)+q**2,0,1,0,2-q**-1-q,0,0,0,0,-1+q,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0],[0,0,0,q-q**2,0,0,-1+q,0,q,0,0,-1+q,0,0,0],[0,1-GAPMul(2,q)+q**2,0,0,0,-1+GAPMul(2,q)-q**2,0,0,0,0,0,q,0,q,q],[j**2+GAPMul(GAPMul(-2,j**2)-j,q)-q**2,0,0,0,j**2-GAPMul(j**2,q),0,0,q,0,0,0,0,0,-1+q,0],[j-GAPMul(j,q),0,0,-q+q**2,j**2+GAPMul(2,j)-GAPMul(j,q),0,1-q,-q,-q,0,0,0,1,0,-1+q]]])
              
              
              def f13(q,j):
                GAPMul(returnq**0,[[[j-GAPMul(j,q),-1+q-q**2,-j**2+GAPMul(2,j**2)-GAPMul(j**2,q**2),-j+GAPMul(-j**2+j,q),-j**2+GAPMul(j**2-j,q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1,j**2-GAPMul(j**2,q),-1+q,1-q,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-q,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,-q,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[-q,0,-q+q**2,0,-q+q**2,1-q,0,j**2-GAPMul(j**2,q),0,0,0,0,0,0,0,q,0,0,0,0],[0,0,1-q,-1+q,0,0,-1+q,0,j**2-GAPMul(j**2,q),0,1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,q,0,0,0,0,0,0,0,0,0,0],[0,0,-j**2+GAPMul(2,j**2)-GAPMul(j**2,q**2),-j+GAPMul(-j**2+j,q),0,0,j**2+GAPMul(-j**2+j,q)-GAPMul(j,q**2),0,1-q+q**2,0,j-GAPMul(j,q),0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,j**2,0,0,0,0,0,0],[0,0,0,0,1-q,-1+q,1-q,0,0,0,0,0,j**2-GAPMul(j**2,q),0,GAPMul(j**2,q),0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,GAPMul(j,q),0,-1+q,0,0,0,0,0,0],[0,0,0,0,2-q**-1-q,GAPMul(2,j**2)+j-GAPMul(j**2,q**-1)+q,GAPMul(-2,j**2)-j-q**-1+GAPMul(j**2,q),0,0,0,0,0,-j+GAPMul(j,q**-1),0,j-GAPMul(j,q),0,0,0,0,0],[0,-1+q-q**2,-q+q**2,-j+GAPMul(-j**2+j,q),1-q,j**2-j+GAPMul(j,q**-1)-GAPMul(j**2,q),0,-1+q**-1,0,0,0,0,0,0,0,j-GAPMul(j,q),0,0,0,0],[0,0,0,0,j**2-j-GAPMul(j**2,q**-1)+GAPMul(j,q),0,GAPMul(-2,j)+GAPMul(j,q**-1),0,0,-j**2+j-GAPMul(j,q**-1)+GAPMul(j**2,q),0,0,0,0,0,0,j-GAPMul(j,q),1-q**-1-q,0,0],[0,0,0,0,-1+q,0,-1+q,0,0,1-q,0,0,0,0,0,0,-q,j**2-GAPMul(j**2,q),0,0],[0,GAPMul(j,q)+GAPMul(j**2-j,q**2)-GAPMul(j**2,q**3),-q+GAPMul(2,q**2)-q**3,0,1-GAPMul(2,q)+q**2,0,1-q-q**2+q**3,-j+GAPMul(-j**2+j,q),GAPMul(j,q)+GAPMul(j**2-j,q**2)-GAPMul(j**2,q**3),0,-q+q**2,0,-j+GAPMul(-j**2+j,q),0,GAPMul(j**2,q)-GAPMul(j**2,q**2),q-q**2,0,0,-1,0],[1-q,j+GAPMul(j**2-j,q)-GAPMul(j**2,q**2),0,0,-1+q**-1-q+q**2,0,-1+q**-1-q+q**2,0,j+GAPMul(j**2-j,q)-GAPMul(j**2,q**2),2-q**-1-q,-1+q,-j+GAPMul(j,q**-1),0,-1+q**-1,0,0,1-q,j**2-j+GAPMul(j,q**-1)-GAPMul(j**2,q),0,-1]],[[-1+q,GAPMul(-j**2,q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-j,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-q,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,-q,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[GAPMul(-j,q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(j,q),0,0,0,0],[0,0,0,0,0,0,0,0,0,0,j,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,GAPMul(j**2,q),0,-1+q,0,0,0,0,0,0,0,0,0],[0,0,-1+q,0,0,0,1-q,0,0,-1+q,0,-1+q,0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0],[0,0,-1+q,0,0,0,0,0,0,q-q**2,0,q,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1+q,0,0,0,0,0],[0,GAPMul(-j**2,q),0,0,0,0,0,j**2,0,0,0,0,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,-j**2,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(-j,q),0,0,0],[0,j-GAPMul(j,q),0,0,0,0,0,j-GAPMul(j,q**-1),j-GAPMul(j,q),0,j**2-GAPMul(j**2,q),0,j-GAPMul(j,q**-1),0,-j+GAPMul(j,q),-j**2+GAPMul(j**2,q),0,0,-1,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1]],[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,-1+q,0,0,j,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,1,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,GAPMul(j**2,q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,j**2,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,GAPMul(j**2,q),0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,GAPMul(j,q),0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,j,0,-1+q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,q,0,-1+q,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0],[0,0,-q+q**2,q-q**2,0,0,q-q**2,0,GAPMul(-j**2,q)+GAPMul(j**2,q**2),0,-q,0,0,q,0,0,0,0,0,0],[0,j+GAPMul(j**2-j,q)-GAPMul(j**2,q**2),-1+GAPMul(2,q)-q**2,0,-2+q**-1,0,-1+q**-1-q+q**2,-j**2+j-GAPMul(j,q**-1)+GAPMul(j**2,q),j+GAPMul(j**2-j,q)-GAPMul(j**2,q**2),0,-1+q,0,-j**2+j-GAPMul(j,q**-1)+GAPMul(j**2,q),0,j**2-GAPMul(j**2,q),1-q,0,0,-1,0],[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,j**2],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0],[-1+q,-j+GAPMul(3,j)-GAPMul(3,j)+GAPMul(j,q**3),j-GAPMul(3,j)+GAPMul(3,j)-GAPMul(j,q**3),-q+GAPMul(2,q**2)-q**3,j**2+GAPMul(3,j)-GAPMul(j,q**-1)+GAPMul(GAPMul(-2,j**2)-GAPMul(4,j),q),-j**2+GAPMul(GAPMul(2,j**2)+j,q),j**2+GAPMul(2,j)-GAPMul(j,q**-1)+GAPMul(GAPMul(-2,j**2)-j,q),GAPMul(-2,j**2)-GAPMul(3,j)-q**-1+GAPMul(j**2+GAPMul(3,j),q)-GAPMul(j,q**2),-j+GAPMul(j**2+GAPMul(3,j),q)-q**3,0,-j**2+GAPMul(j**2-j,q),0,2-q**-1-GAPMul(2,q)+q**2,0,-1+q-q**2,-j+GAPMul(-j**2+j,q),0,0,j-GAPMul(j,q),0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(j,q),0,0,-1+q]],[[q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0,0,0],[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,q,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,q,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,q,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(j**2,q),0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-j**2,0,0],[0,0,0,0,0,0,0,0,q,0,0,0,-1+q,0,0,0,0,0,0,0],[0,0,0,0,1-q,0,1-q,0,0,-1+q,0,0,0,0,0,0,q,-j**2+GAPMul(j**2,q),0,0],[0,0,0,0,0,0,0,0,0,0,j,0,0,0,-1+q,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0],[0,0,-1+q,0,0,0,1-q,0,0,-1+q,0,-1+q,0,1,0,0,-1+q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,GAPMul(-j,q),0,0,0,0,0,-1+q,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0],[0,GAPMul(-j,q)+GAPMul(-j**2+j,q**2),q-GAPMul(2,q**2)+q**3,0,-1+GAPMul(2,q)-q**2,0,-1+q-q**3,j+GAPMul(j**2-j,q)-GAPMul(j**2,q**2),GAPMul(-j,q)+GAPMul(-j**2+j,q**2),0,q-q**2,0,j+GAPMul(j**2-j,q)-GAPMul(j**2,q**2),0,GAPMul(-j**2,q)+GAPMul(j**2,q**2),-q+q**2,0,0,q,q]],[[0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,q,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,GAPMul(j,q),0,0,0,0,0,0,0,0],[0,q,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0],[-q,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,j**2,0,0,0,-1+q,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0],[-q,0,-q+q**2,0,-q+q**2,1-q,0,j**2-GAPMul(j**2,q),0,0,0,0,0,-1+q,0,q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(j,q),0,0,0],[0,0,-1+q,0,0,0,1-q,0,0,-1+q,-1,-1+q,0,1,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,j**2,0,-1+q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,-1+q,0,0],[q-q**2,GAPMul(j,q)+GAPMul(j**2-j,q**2)-GAPMul(j**2,q**3),0,0,1-q-q**2+q**3,0,1-q-q**2+q**3,0,GAPMul(j,q)+GAPMul(j**2-j,q**2)-GAPMul(j**2,q**3),-1+GAPMul(2,q)-q**2,-q+q**2,j-GAPMul(j,q),0,1-q,0,0,q-q**2,j+GAPMul(j**2-j,q)-GAPMul(j**2,q**2),0,-q],[0,j+GAPMul(j**2-j,q)-GAPMul(j**2,q**2),-1+GAPMul(2,q)-q**2,0,-2+q**-1,0,-1+q**-1-q+q**2,-j**2+j-GAPMul(j,q**-1)+GAPMul(j**2,q),j+GAPMul(j**2-j,q)-GAPMul(j**2,q**2),0,-1+q,0,-j**2+j-GAPMul(j,q**-1)+GAPMul(j**2,q),0,j**2-GAPMul(j**2,q),1-q,0,0,-1,-1+q]]])
              
              
              def f17(q):
                GAPMul(returnq**0,[[[-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[-q**4,q**3,-q**2,q]],[[-1,0,0,0],[0,-1,0,0],[0,0,-1,0],[-q**4,q**3,-q**2,q]],[[-1,0,0,0],[0,-1,0,0],[0,0,0,1],[0,0,q,-1+q]],[[-1,0,0,0],[0,0,1,0],[0,q,-1+q,0],[0,0,0,-1]],[[0,1,0,0],[q,-1+q,0,0],[0,0,-1,0],[0,0,0,-1]]])
              
              
              def f20(q,j):
                GAPMul(returnq**0,[[[q,0,-q**2,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0],[0,0,0,-q**2,q,0,0,0,0,0],[0,-q**2,0,0,0,q,0,0,0,0],[GAPMul(-j,q)+GAPMul(-j**2+j,q**2),GAPMul(j**2,q)+GAPMul(-j**2+j,q**2)-GAPMul(j,q**3),GAPMul(-j**2,q**2)+GAPMul(2,j**2)-GAPMul(j**2,q**4),GAPMul(j**2,q)-GAPMul(2,j**2)+GAPMul(j**2,q**3),j+GAPMul(j**2-j,q)-GAPMul(j**2,q**2),0,j-GAPMul(j,q),-j+GAPMul(j,q)-GAPMul(j,q**2),0,0],[GAPMul(j**2,q)-GAPMul(j**2,q**2),GAPMul(j**2,q)-GAPMul(j**2,q**2),GAPMul(-j**2,q**2)+GAPMul(j**2,q**3),GAPMul(j**2,q)-GAPMul(j**2,q**2),-j**2+GAPMul(j**2,q),0,-j**2,j**2-GAPMul(j**2,q),0,0],[GAPMul(-j,q**2)+GAPMul(-j**2+j,q**3),GAPMul(-j,q**2)+GAPMul(2,j)-GAPMul(j,q**4),GAPMul(j**2,q**2)+GAPMul(GAPMul(-2,j**2)+j,q**3)-GAPMul(j**2,q**5),GAPMul(j**2,q**2)+GAPMul(-j**2+j,q**3)-GAPMul(j,q**4),0,j+GAPMul(j**2-j,q)-GAPMul(j**2,q**2),0,GAPMul(-j,q)+GAPMul(j,q**2)-GAPMul(j,q**3),j-GAPMul(j,q),GAPDiv(j-GAPMul(j,q)+GAPMul(j,q**2),q)],[0,0,GAPMul(-j**2,q**3)+GAPMul(j**2,q**4),0,GAPMul(-j**2,q**2)+GAPMul(j**2,q**3),GAPMul(j**2,q)-GAPMul(j**2,q**2),GAPMul(-j**2,q**2),0,GAPMul(j**2,q),j**2-GAPMul(j**2,q)]],[[q,0,-q**2,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0],[0,0,0,-q**2,q,0,0,0,0,0],[0,-q**2,0,0,0,q,0,0,0,0],[0,0,0,0,0,0,-1+q,-q,0,0],[0,0,0,0,0,0,-1,0,0,0],[0,0,0,0,0,0,0,-q**2,-1+q,1],[0,0,0,0,0,0,-q**2,0,q,0]],[[-1+q,0,q,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,-1,0,0],[1,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0],[0,0,0,q,-1+q,0,0,0,0,0],[0,0,0,0,0,0,0,-q,0,GAPDiv(1,q)],[0,0,0,0,0,0,-1,0,0,0],[0,-q,0,0,0,0,0,-1+q,0,0],[0,0,0,0,0,0,-q**2,0,q,0],[0,-q**3,0,0,0,q**2,0,0,0,-1+q]],[[q,0,-q**2,0,0,0,0,0,0,0],[0,-1+q,0,q,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0],[0,0,0,0,q,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,1,0],[0,0,0,0,0,0,0,q,0,GAPDiv(-1,q)],[0,0,0,0,0,0,q,0,-1+q,0],[0,0,0,0,0,0,0,0,0,-1]],[[0,0,0,0,1,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0],[0,0,0,1,0,0,0,0,0,0],[0,0,q,-1+q,0,0,0,0,0,0],[q,0,0,0,-1+q,0,0,0,0,0],[0,-q**2,0,0,0,q,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0],[0,0,0,0,0,0,0,-1,0,0],[0,0,0,0,0,0,-q**2,0,q,0],[0,0,0,0,0,0,0,-q**2,0,q]]])
              
              
              def f23(q):
                GAPMul(returnq**0,[[[-1,0,0,0,0],[0,-1,0,0,0],[q**3,-q**2,q,0,0],[0,0,0,-1,0],[-q**3,0,0,-q**2,q]],[[-1,0,0,0,0],[0,-1,0,0,0],[q**3,-q**2,q,0,0],[0,0,0,-1,0],[-q**3,0,0,-q**2,q]],[[-1,0,0,0,0],[0,0,1,0,0],[0,q,-1+q,0,0],[0,0,0,0,1],[0,0,0,q,-1+q]],[[0,1,0,0,0],[q,-1+q,0,0,0],[0,0,-1,0,0],[0,0,0,-1,0],[-q**4,q**3,-q**2,-q**2,q]],[[-1,0,0,0,0],[0,0,0,1,0],[0,0,0,0,1],[0,q,0,-1+q,0],[0,0,q,0,-1+q]]])
              
              
              def f29(q):
                GAPMul(returnq**0,[[[q,q**4,0,0,-q**2,0],[0,-1,0,0,0,0],[0,0,-1,0,0,0],[0,q**3,-q**2,q,0,0],[0,0,0,0,-1,0],[0,0,q**4,0,-q**3,q]],[[q,q**4,0,0,-q**2,0],[0,-1,0,0,0,0],[0,0,-1,0,0,0],[0,q**3,-q**2,q,0,0],[0,0,0,0,-1,0],[0,0,q**4,0,-q**3,q]],[[-1+q,0,0,0,q,0],[0,-1,0,0,0,0],[0,0,0,1,0,0],[0,0,q,-1+q,0,0],[1,0,0,0,0,0],[0,0,0,0,0,q]],[[0,0,0,0,0,1],[0,0,1,0,0,0],[0,q,-1+q,0,0,0],[0,0,0,-1,0,0],[0,0,0,0,q,0],[q,0,0,0,0,-1+q]],[[-1+q,0,0,q,0,0],[0,q,0,0,0,0],[0,0,0,0,1,0],[1,0,0,0,0,0],[0,0,q,0,-1+q,0],[0,0,0,0,0,-1]]])
              
              
              def f9(q):
                GAPMul(returnq**0,[[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-q,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-q,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[ER(3)**2-GAPMul(2,ER(3)**2)+GAPMul(ER(3)**2,q**2),ER(3)+GAPMul(ER(3)**2-ER(3),q)-GAPMul(ER(3)**2,q**2),ER(3)**2+GAPMul(-ER(3)**2+ER(3),q)-GAPMul(ER(3),q**2),0,0,0,0,ER(3)-GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)-GAPMul(ER(3),q)+GAPMul(ER(3),q**2),0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2,0,0,0,0,0],[-1+GAPMul(2,q)-q**2,1-GAPMul(2,q)+q**2,-1+GAPMul(2,q)-q**2,-1+q,0,-1+q,0,1-q,0,ER(3)-GAPMul(ER(3),q),0,0,0,GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1-q,0,1-GAPMul(2,q)+q**2,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-3+q**-1-q**2,3-q**-1-GAPMul(3,q)+q**2,-3+q**-1-q**2,-2+q**-1,0,ER(3)**2+GAPMul(2,ER(3))-GAPMul(ER(3),q),0,2-q**-1-q,0,-ER(3)**2+GAPMul(ER(3)**2,q**-1),0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-ER(3)**2-GAPMul(2,ER(3))+GAPMul(ER(3),q**-1)-q,0,2-q**-1-GAPMul(2,q)+q**2,0,0],[-ER(3)+GAPMul(ER(3),q),ER(3)-GAPMul(ER(3),q),0,GAPMul(-2,ER(3))+GAPMul(ER(3),q**-1),0,0,0,0,0,0,0,-ER(3)+GAPMul(ER(3),q),0,0,ER(3)-GAPMul(ER(3),q),0,0,0,ER(3)**2,0,0,0,ER(3)-GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),0,0,0,0],[0,0,0,0,0,0,0,0,0,0,-ER(3)**2+GAPMul(ER(3)**2,q),0,ER(3)**2-GAPMul(ER(3)**2,q),0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,0,0,ER(3)**2,0,0,0,0,0,0,-ER(3)**2+GAPMul(ER(3)**2,q),0,0,0,0],[0,0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q**-1),0,0,0,0,0,-ER(3)**2+GAPMul(ER(3)**2,q),0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,0,ER(3)**2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-GAPMul(2,q)+q**2,ER(3)+GAPMul(-ER(3)**2-GAPMul(2,ER(3)),q)-q**2,0,3-q**-1-GAPMul(3,q)+q**2,0,0,0,0,0,0,-1+q,1+GAPMul(ER(3)**2+GAPMul(2,ER(3)),q)-GAPMul(ER(3),q**2),0,0,ER(3)-GAPMul(ER(3),q)+GAPMul(ER(3),q**2),0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,0,0,-1+GAPMul(2,q)-q**2,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2+GAPMul(3,ER(3))-GAPMul(ER(3),q**-1)+GAPMul(GAPMul(-2,ER(3)**2)-GAPMul(3,ER(3)),q)-q**2,0,0,0,0],[0,0,0,0,0,GAPMul(-2,ER(3)**2)+GAPMul(ER(3)**2,q**-1),0,0,0,0,0,ER(3)**2+GAPMul(-ER(3)**2+ER(3),q)-GAPMul(ER(3),q**2),0,0,0,0,ER(3)-GAPMul(ER(3),q)+GAPMul(ER(3),q**2),0,0,ER(3)-GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)+GAPMul(ER(3)**2-ER(3),q)-GAPMul(ER(3)**2,q**2),0,0,0],[0,0,-ER(3)+GAPMul(ER(3),q),0,ER(3)-GAPMul(ER(3),q),0,-ER(3)+GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0,ER(3)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),0,0,ER(3)-GAPMul(ER(3),q),ER(3)-GAPMul(ER(3),q),0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1+q,0,0,-2+q**-1,0,0,0,0,0,0,2-q**-1-q,-1+q,-2+q**-1,0,0,-2+q**-1,0,1-q,0,0,0,ER(3)-GAPMul(ER(3),q),1-q,0,0,ER(3),0,0,-1+q**-1,0,0,0,0,0,0,2-q**-1-q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,ER(3),0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,1-q**-1,0,0,0,0,0,-1+q,0,0,0,0,0,1-q,0,1,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0],[0,0,GAPMul(-2,ER(3)**2)+GAPMul(ER(3)**2,q**-1),0,GAPMul(2,ER(3)**2)+ER(3)-GAPMul(ER(3)**2,q),0,GAPMul(-2,ER(3)**2)-ER(3)+GAPMul(ER(3)**2,q**-1)-q,0,0,0,0,ER(3)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),0,0,0,0,GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**-1)-GAPMul(2,ER(3)**2)+GAPMul(ER(3)**2,q**2),0,0,GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**-1)-GAPMul(ER(3)**2,q),-1+q**-1,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0,GAPMul(-2,ER(3))+GAPMul(ER(3),q**-1),0,0,0],[1+GAPMul(ER(3)**2+GAPMul(2,ER(3)),q)-GAPMul(ER(3),q**2),0,0,GAPMul(-2,ER(3)**2)-GAPMul(3,ER(3))-q**-1+GAPMul(ER(3)**2+GAPMul(3,ER(3)),q)-GAPMul(ER(3),q**2),0,0,0,0,0,0,-3+q**-1-q**2,-ER(3)+GAPMul(2,ER(3))-GAPMul(ER(3),q**2),3-q**-1-GAPMul(3,q)+q**2,0,0,2-q**-1-GAPMul(2,q)+q**2,0,ER(3)+GAPMul(-ER(3)**2-GAPMul(2,ER(3)),q)-q**2,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q)+GAPMul(ER(3)**2,q**2),-1+GAPMul(2,q)-q**2,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,0,2-q**-1-q,0,0,0,0,0,0,-3+q**-1-q**2,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-ER(3)**2,0,0,0,0,0,0],[0,0,-1+q,0,0,0,-1+q,0,0,0,0,0,1-q,0,0,0,0,0,0,1-q,0,0,0,ER(3)**2-GAPMul(2,ER(3)**2)+GAPMul(ER(3)**2,q**2),0,0,0,ER(3)-GAPMul(ER(3),q),0,q,1-q,0,0,0,0,0,0,0,0,0],[0,0,0,-q+q**2,0,0,1-q,0,0,0,-ER(3)+GAPMul(ER(3),q)-GAPMul(ER(3),q**2),0,ER(3)+GAPMul(ER(3)**2-ER(3),q)-GAPMul(ER(3)**2,q**2),0,0,ER(3)-GAPMul(ER(3),q)+GAPMul(ER(3),q**2),0,0,0,0,0,0,0,0,0,0,0,0,ER(3)-GAPMul(ER(3),q),0,0,0,0,0,0,-ER(3)+GAPMul(-ER(3)**2+ER(3),q),0,0,0,0],[0,0,ER(3)**2-ER(3)+GAPMul(ER(3),q**-1)-GAPMul(ER(3)**2,q),0,0,-2+q**-1,GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**-1)-GAPMul(ER(3)**2,q),0,0,0,0,3-q**-1-GAPMul(3,q)+q**2,-ER(3)**2+ER(3)-GAPMul(ER(3),q),0,0,0,0,-2+q**-1,0,GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),0,0,0,2-q**-1-GAPMul(2,q)+q**2,0,0,0,-1+q**-1,0,ER(3)**2-GAPMul(ER(3)**2,q),GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0],[0,0,0,0,0,-1+q,0,0,0,0,0,1-q,0,0,0,0,-ER(3)+GAPMul(ER(3),q)-GAPMul(ER(3),q**2),ER(3)+GAPMul(ER(3)**2-ER(3),q)-GAPMul(ER(3)**2,q**2),0,0,0,0,0,1-q+q**2,0,0,0,0,0,0,ER(3)-GAPMul(ER(3),q),0,0,0,0,0,-ER(3)+GAPMul(-ER(3)**2+ER(3),q),0,0,0],[-1+GAPMul(2,q)-q**2,GAPMul(2,ER(3)**2)-ER(3)-GAPMul(ER(3)**2,q**-1)+GAPMul(-ER(3)**2+GAPMul(2,ER(3)),q)-GAPMul(ER(3),q**2),ER(3)+GAPMul(ER(3)**2-GAPMul(2,ER(3)),q),GAPMul(7,ER(3)**2)+GAPMul(ER(3)**2,q**-2)-GAPMul(4,ER(3)**2)-GAPMul(7,ER(3)**2)+GAPMul(4,ER(3)**2)-GAPMul(ER(3)**2,q**3),0,-ER(3)**2+GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)+GAPMul(GAPMul(2,ER(3)**2)-ER(3),q)-GAPMul(ER(3)**2,q**2),GAPMul(3,ER(3)**2)-GAPMul(ER(3)**2,q**-1)-GAPMul(3,ER(3)**2)+GAPMul(ER(3)**2,q**2),GAPMul(2,ER(3)**2)-ER(3)-GAPMul(ER(3)**2,q**-1)+GAPMul(GAPMul(-2,ER(3)**2)+ER(3),q),0,2-q**-1-GAPMul(2,q)+q**2,GAPMul(-5,ER(3)**2)+GAPMul(2,ER(3))-GAPMul(ER(3)**2,q**-2)+GAPMul(GAPMul(3,ER(3)**2)-ER(3),q**-1)-GAPMul(2,ER(3)**2),GAPMul(3,ER(3)**2)-ER(3)-GAPMul(ER(3)**2,q**-1)+GAPMul(GAPMul(-4,ER(3)**2)+GAPMul(2,ER(3)),q)-GAPMul(ER(3)**2,q**3),ER(3)**2-GAPMul(2,ER(3))+GAPMul(ER(3)**2,q**-2),-1+GAPMul(2,q)-q**2,-ER(3)**2+GAPMul(2,ER(3)**2)-GAPMul(2,ER(3)**2)+GAPMul(ER(3)**2,q**3),GAPMul(2,ER(3)**2)-GAPMul(2,ER(3))+GAPMul(ER(3)**2,q**-2),0,GAPMul(-2,ER(3)**2)+ER(3),-1+GAPMul(2,q)-q**2,GAPMul(-2,ER(3)**2)+ER(3)-GAPMul(ER(3)**2,q**2),0,-2+q**-1-q**2,GAPMul(-3,ER(3)**2)+ER(3),-ER(3)**2-GAPMul(3,ER(3))+GAPMul(ER(3),q**-1),0,-2+q**-1,ER(3)-GAPMul(ER(3),q)+GAPMul(ER(3),q**2),-2+q**-1-q**2,GAPMul(2,ER(3)**2)-ER(3)+GAPMul(ER(3)**2,q**-2)-GAPMul(ER(3)**2,q),ER(3)**2-GAPMul(2,ER(3)**2)+GAPMul(ER(3)**2,q**2),GAPMul(-2,ER(3)**2)+ER(3)-GAPMul(ER(3)**2,q**2),ER(3)-GAPMul(ER(3),q),1-q+q**2,0,0,GAPMul(-5,ER(3)**2)+GAPMul(2,ER(3))-GAPMul(ER(3)**2,q**-2)+GAPMul(GAPMul(3,ER(3)**2)-ER(3),q**-1)-GAPMul(4,ER(3)**2)+GAPMul(ER(3)**2,q**3),0,GAPMul(2,ER(3)**2)-ER(3)-GAPMul(ER(3)**2,q**-1)+GAPMul(GAPMul(-2,ER(3)**2)+GAPMul(2,ER(3)),q)-GAPMul(ER(3)**2,q**3),0,0],[0,-2+q**-1,1-GAPMul(2,q)+q**2,-4-q**-2+GAPMul(3,q**-1)-q**2,0,2-q**-1-q,-2+q**-1,-2+q**-1,0,GAPMul(-2,ER(3))+GAPMul(ER(3),q**-1),GAPMul(-3,ER(3)**2)-GAPMul(4,ER(3))+q**-2-GAPMul(3,q**-1)+GAPMul(ER(3)**2+GAPMul(2,ER(3)),q),-3+q**-1-q**2,-1-q**-2+GAPMul(2,q**-1),ER(3)-GAPMul(ER(3),q),1-GAPMul(2,q)+q**2,GAPMul(3,ER(3)**2)+GAPMul(2,ER(3))-q**-2+GAPMul(GAPMul(-3,ER(3)**2)-GAPMul(2,ER(3)),q**-1),0,2-q**-1-q,ER(3)-GAPMul(ER(3),q),2-q**-1-q,0,GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),3-q**-1-GAPMul(3,q)+q**2,GAPMul(2,ER(3)**2)-ER(3)-GAPMul(ER(3)**2,q**-1)+GAPMul(GAPMul(-2,ER(3)**2)+ER(3),q),0,ER(3)-GAPMul(ER(3),q**-1),1-q,GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),-1-q**-2+GAPMul(2,q**-1),-1+q,2-q**-1-q,1,ER(3)**2-GAPMul(ER(3)**2,q),1,0,4+q**-2-GAPMul(3,q**-1)-GAPMul(3,q)+q**2,0,GAPMul(3,ER(3)**2)+GAPMul(2,ER(3))-q**2,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(-ER(3),q),0,0,0,0,0,0,-1+q,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0],[0,0,0,-q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0],[-ER(3)**2+GAPMul(ER(3)**2,q),ER(3)**2-GAPMul(ER(3)**2,q),-ER(3)**2+GAPMul(ER(3)**2,q),0,0,0,0,ER(3)**2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,0],[0,0,0,-1+q,0,-1+q,0,0,0,0,0,0,0,GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1-q,0,1-q,0,ER(3)-GAPMul(ER(3),q),GAPMul(ER(3),q)],[3-q**-1-GAPMul(3,q)+q**2,-3+q**-1-q**2,3-q**-1-GAPMul(3,q)+q**2,-ER(3)**2+GAPMul(ER(3)**2,q),0,-ER(3)**2+GAPMul(ER(3)**2,q**-1),0,-2+q**-1,-1+q,ER(3)**2-GAPMul(ER(3)**2,q**-1)-GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2-q**-1-q,ER(3)**2+GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)+q,-ER(3)**2-GAPMul(2,ER(3))+GAPMul(ER(3),q**-1)-q,-2+q**-1-q**2,-ER(3)**2+GAPMul(ER(3)**2,q**-1),ER(3)**2-GAPMul(ER(3)**2,q)]],[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-q,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-q,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0],[1-q,0,1-q,0,1-q**-1,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1-q,0,0,1-q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,-1,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,1,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1-q,0,0,1-q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,1,0,0,0,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-2+q**-1,0,0,-1+q**-1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1-q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0,0,0,0,0,0,0,0,0,0,ER(3),0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-2+q**-1,2-q**-1-q,-1+GAPMul(2,q)-q**2,3-q**-1-GAPMul(3,q)+q**2,-2+q**-1,0,2-q**-1-q,2-q**-1-q,-ER(3)+GAPMul(ER(3),q),GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),ER(3)-GAPMul(ER(3),q),1-GAPMul(2,q)+q**2,0,0,-1+GAPMul(2,q)-q**2,0,GAPMul(3,ER(3)**2)+GAPMul(2,ER(3))-q**2,0,-ER(3)+GAPMul(ER(3),q),-2+q**-1,-2+q**-1,0,-1+GAPMul(2,q)-q**2,0,ER(3)**2-GAPMul(ER(3)**2,q),0,-1+q,0,0,0,0,0,0,-1,-2+q**-1,-1+GAPMul(2,q)-q**2,-2+q**-1,GAPMul(-3,ER(3)**2)-GAPMul(2,ER(3))-q**-1+GAPMul(GAPMul(3,ER(3)**2)+GAPMul(2,ER(3)),q),GAPMul(-2,ER(3))+GAPMul(ER(3),q**-1),ER(3)-GAPMul(ER(3),q)],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,-q,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0,-1+q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,-1+q,GAPMul(ER(3)**2,q),0,0,0,0,0,0,0],[GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),GAPMul(-2,ER(3))+GAPMul(ER(3),q**-1),ER(3)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),GAPMul(-3,ER(3))+GAPMul(ER(3),q**-1)-GAPMul(ER(3),q**2),GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),0,GAPMul(-2,ER(3))+GAPMul(ER(3),q**-1),GAPMul(-2,ER(3))+GAPMul(ER(3),q**-1),ER(3)**2-GAPMul(ER(3)**2,q),GAPMul(-2,ER(3)**2)+GAPMul(ER(3)**2,q**-1),-ER(3)**2+GAPMul(ER(3)**2,q),-ER(3)+GAPMul(2,ER(3))-GAPMul(ER(3),q**2),0,0,ER(3)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),0,ER(3)**2+GAPMul(3,ER(3))-GAPMul(ER(3),q**-1)+GAPMul(-ER(3)**2-GAPMul(3,ER(3)),q),0,ER(3)**2-GAPMul(ER(3)**2,q),GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),0,ER(3)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),0,-1+q,0,ER(3)-GAPMul(ER(3),q),0,0,0,0,ER(3),0,ER(3),GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),ER(3)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),-ER(3)**2-GAPMul(3,ER(3))+GAPMul(ER(3),q**-1)-GAPMul(ER(3),q**2),GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**-1)-GAPMul(ER(3)**2,q),-ER(3)**2+GAPMul(ER(3)**2,q)],[-ER(3)+GAPMul(2,ER(3))-GAPMul(ER(3),q**2),-1+GAPMul(2,q)-q**2,-ER(3)+GAPMul(2,ER(3))-GAPMul(ER(3),q**2),-ER(3)**2+GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**2),GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**-1)-GAPMul(ER(3)**2,q),2-q**-1-q,-ER(3)**2+GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**2),ER(3)-GAPMul(ER(3),q),-ER(3)**2+GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**2),-ER(3)+GAPMul(ER(3),q),ER(3)**2-GAPMul(2,ER(3)**2)+GAPMul(ER(3)**2,q**2),-1+GAPMul(2,q)-q**2,0,0,-q+q**2,0,1-GAPMul(2,q)+q**2,0,0,-ER(3)+GAPMul(ER(3),q),1-q,0,GAPMul(ER(3),q)-GAPMul(ER(3),q**2),0,0,0,-q,0,0,0,0,0,0,0,-ER(3)+GAPMul(ER(3),q),ER(3)**2-GAPMul(2,ER(3)**2)+GAPMul(ER(3)**2,q**2),1-GAPMul(2,q)+q**2,-1+GAPMul(2,q)-q**2,ER(3)-GAPMul(ER(3),q),0],[1-q,0,0,0,1-q,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0],[0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0],[0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q],[0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1+q]],[[0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0],[0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,ER(3)**2,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,1,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2,0,0,0,0,0,0,0,0,0,0,0,0],[-2+q**-1,2-q**-1-q,-1+GAPMul(2,q)-q**2,3-q**-1-GAPMul(3,q)+q**2,-2+q**-1,0,2-q**-1-q,2-q**-1-q,-ER(3)+GAPMul(ER(3),q),GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),ER(3)-GAPMul(ER(3),q),1-GAPMul(2,q)+q**2,0,0,-1+GAPMul(2,q)-q**2,0,GAPMul(3,ER(3)**2)+GAPMul(2,ER(3))-q**2,0,-ER(3)+GAPMul(ER(3),q),-2+q**-1,-2+q**-1,0,-1+GAPMul(2,q)-q**2,0,ER(3)**2-GAPMul(ER(3)**2,q),0,-1+q,0,0,0,0,0,0,-1,-2+q**-1,-1+GAPMul(2,q)-q**2,-2+q**-1,GAPMul(-3,ER(3)**2)-GAPMul(2,ER(3))-q**-1+GAPMul(GAPMul(3,ER(3)**2)+GAPMul(2,ER(3)),q),GAPMul(-2,ER(3))+GAPMul(ER(3),q**-1),ER(3)-GAPMul(ER(3),q)],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3),q),0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1-q,0,0,1-q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,q,0,0,0,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3),q),0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3),q),0,-1+q,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2,0,-1+q,0,0,0,0,0,0,0],[GAPMul(-3,ER(3)**2)+GAPMul(ER(3)**2,q**-1)-q**3,GAPMul(3,ER(3)**2)-GAPMul(ER(3)**2,q**-1)+GAPMul(GAPMul(-4,ER(3)**2)-ER(3),q),GAPMul(-4,ER(3)**2)-GAPMul(2,ER(3))-q**-1+GAPMul(GAPMul(6,ER(3)**2)+GAPMul(2,ER(3)),q),ER(3)**2-GAPMul(3,ER(3))+GAPMul(ER(3),q**-1),GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),2-q**-1-q,GAPMul(-2,ER(3))+GAPMul(ER(3),q**-1)-GAPMul(ER(3),q**2),GAPMul(3,ER(3)**2)-GAPMul(ER(3)**2,q**-1)-GAPMul(3,ER(3)**2)+GAPMul(ER(3)**2,q**2),0,3-q**-1-GAPMul(3,q)+q**2,0,ER(3)**2+GAPMul(GAPMul(-2,ER(3)**2)+ER(3),q),0,ER(3)-GAPMul(ER(3),q),-ER(3)**2+GAPMul(GAPMul(2,ER(3)**2)-ER(3),q)-GAPMul(ER(3),q**3),0,GAPMul(3,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(4,ER(3))+GAPMul(3,ER(3))-GAPMul(ER(3),q**3),0,-1+GAPMul(GAPMul(-2,ER(3)**2)-ER(3),q),GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),0,-ER(3)**2+GAPMul(GAPMul(2,ER(3)**2)-ER(3),q)-GAPMul(ER(3),q**3),0,-1+q-q**2,0,0,0,0,0,0,0,0,ER(3)-GAPMul(ER(3),q),-ER(3)**2+GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)+GAPMul(GAPMul(2,ER(3)**2)-ER(3),q)-GAPMul(ER(3)**2,q**2),ER(3)-GAPMul(3,ER(3))+GAPMul(3,ER(3))-GAPMul(ER(3),q**3),-ER(3)**2+GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)+GAPMul(GAPMul(2,ER(3)**2)-ER(3),q)-GAPMul(ER(3)**2,q**2),GAPMul(4,ER(3)**2)-GAPMul(ER(3)**2,q**-1)-GAPMul(6,ER(3)**2)+GAPMul(4,ER(3)**2)-GAPMul(ER(3)**2,q**3),GAPMul(3,ER(3)**2)+ER(3)-GAPMul(ER(3)**2,q**-1)+GAPMul(GAPMul(-3,ER(3)**2)-GAPMul(2,ER(3)),q)-q**2,-ER(3)**2+GAPMul(GAPMul(2,ER(3)**2)+ER(3),q)],[q-q**2,-q+q**2,q-q**2,0,0,0,0,-q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,-q+q**2,0,0],[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0],[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,-1+q,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q]],[[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[q,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,q,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,q,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[q-q**2,0,q-q**2,0,-1+q,0,0,0,-q+q**2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0],[0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1-q,0,0,1-q,0,0,0,-1+q,0,0,0,0,0,0,0,0,-1+q,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,1-q**-1,0,0,0,0,0,-1+q,0,0,0,0,1-q,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1-q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1+q,0,0,0,0,0,0,0,0,0],[-1+GAPMul(2,q)-q**2,1-GAPMul(2,q)+q**2,q-GAPMul(2,q**2)+q**3,1-GAPMul(3,q)+GAPMul(3,q**2)-q**3,-1+GAPMul(2,q)-q**2,0,1-GAPMul(2,q)+q**2,1-GAPMul(2,q)+q**2,GAPMul(ER(3),q)-GAPMul(ER(3),q**2),ER(3)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),GAPMul(-ER(3),q)+GAPMul(ER(3),q**2),-q+GAPMul(2,q**2)-q**3,0,0,q-GAPMul(2,q**2)+q**3,0,-1+GAPMul(GAPMul(-3,ER(3)**2)-GAPMul(2,ER(3)),q),0,GAPMul(ER(3),q)-GAPMul(ER(3),q**2),-1+GAPMul(2,q)-q**2,-1+GAPMul(2,q)-q**2,0,q-GAPMul(2,q**2)+q**3,0,GAPMul(-ER(3)**2,q)+GAPMul(ER(3)**2,q**2),0,q-q**2,0,0,0,0,q,0,q,-1+GAPMul(2,q)-q**2,q-GAPMul(2,q**2)+q**3,-1+GAPMul(2,q)-q**2,1+GAPMul(GAPMul(3,ER(3)**2)+GAPMul(2,ER(3)),q)-q**3,-ER(3)+GAPMul(2,ER(3))-GAPMul(ER(3),q**2),GAPMul(-ER(3),q)+GAPMul(ER(3),q**2)],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0,0,0,0,0,q,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0],[-1+q,1-q,-1+q,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,1-q,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1+q,0,0,0],[0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1]],[[0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0],[0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[q,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1-q,0,0,1-q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0],[0,0,0,0,0,1,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,ER(3),0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0],[-1+q,1-q,-1+q,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1-q,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,q],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[q-q**2,0,q-q**2,0,-1+q,0,0,-q,-q+q**2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,q,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0],[-2+q**-1,2-q**-1-q,-1+GAPMul(2,q)-q**2,3-q**-1-GAPMul(3,q)+q**2,-2+q**-1,0,2-q**-1-q,2-q**-1-q,-ER(3)+GAPMul(ER(3),q),GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),ER(3)-GAPMul(ER(3),q),1-GAPMul(2,q)+q**2,0,0,-1+GAPMul(2,q)-q**2,0,GAPMul(3,ER(3)**2)+GAPMul(2,ER(3))-q**2,0,-ER(3)+GAPMul(ER(3),q),-2+q**-1,-2+q**-1,0,-1+GAPMul(2,q)-q**2,0,ER(3)**2-GAPMul(ER(3)**2,q),0,-1+q,0,0,0,0,-1+q,0,-1,-2+q**-1,-1+GAPMul(2,q)-q**2,-2+q**-1,GAPMul(-3,ER(3)**2)-GAPMul(2,ER(3))-q**-1+GAPMul(GAPMul(3,ER(3)**2)+GAPMul(2,ER(3)),q),GAPMul(-2,ER(3))+GAPMul(ER(3),q**-1),ER(3)-GAPMul(ER(3),q)],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0],[0,-1+GAPMul(2,q)-q**2,-q+GAPMul(2,q**2)-q**3,-3+q**-1-GAPMul(3,q**2)+q**3,0,1-GAPMul(2,q)+q**2,-1+GAPMul(2,q)-q**2,-1+GAPMul(2,q)-q**2,0,-ER(3)+GAPMul(2,ER(3))-GAPMul(ER(3),q**2),3-q**-1+GAPMul(GAPMul(3,ER(3)**2)+GAPMul(4,ER(3)),q),-1+GAPMul(3,q)-GAPMul(3,q**2)+q**3,-2+q**-1,GAPMul(-ER(3),q)+GAPMul(ER(3),q**2),-q+GAPMul(2,q**2)-q**3,GAPMul(3,ER(3)**2)+GAPMul(2,ER(3))-q**2,0,1-GAPMul(2,q)+q**2,GAPMul(-ER(3),q)+GAPMul(ER(3),q**2),1-GAPMul(2,q)+q**2,0,ER(3)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),1-GAPMul(3,q)+GAPMul(3,q**2)-q**3,ER(3)**2+GAPMul(GAPMul(-2,ER(3)**2)+ER(3),q)-GAPMul(ER(3)**2,q**3),0,ER(3)-GAPMul(ER(3),q),-q+q**2,ER(3)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),-2+q**-1,q-q**2,1-GAPMul(2,q)+q**2,-q,GAPMul(-ER(3)**2,q)+GAPMul(ER(3)**2,q**2),0,0,3-q**-1-GAPMul(4,q)+GAPMul(3,q**2)-q**3,0,-1+GAPMul(GAPMul(-3,ER(3)**2)-GAPMul(2,ER(3)),q),0,0],[0,0,0,0,0,0,0,0,0,0,-1+q,0,1-q,0,0,1-q,0,0,0,0,0,0,q,0,0,0,0,0,1,0,0,0,0,0,0,-1+q,0,0,0,0],[0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]])
              
              
              def f8(q):
                GAPMul(returnq**0,[[[ER(3)-GAPMul(ER(3),q),GAPMul(ER(3),q),0,-q+q**2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-ER(3)**2+GAPMul(ER(3)**2,q**-1),ER(3)**2-GAPMul(ER(3)**2,q),0,-ER(3)+GAPMul(GAPDiv(-3+ER(-3),2),q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,ER(3),-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-2+q**-1,1-q,-ER(3)**2+GAPMul(ER(3)**2,q**-1),0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),GAPDiv(5+ER(-3),2)-q**-1+GAPMul(-2-ER(-3),q)-GAPMul(ER(3)**2,q**2),0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,ER(3)**2,0,0,0,0,0,0,0,0,0,0,0,-ER(3)**2+GAPMul(ER(3)**2,q),0,0,0,0,0,0,0],[GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**-1)-GAPMul(2,ER(3)**2)+GAPMul(ER(3)**2,q**2),-ER(3)**2+GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**2),GAPDiv(-3+ER(-3),2)+q**-1-GAPMul(ER(3),q),-ER(3)**2+GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**2),0,0,0,1-q+q**2,0,ER(3)-GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,ER(3)-GAPMul(ER(-3),q)-GAPMul(ER(3)**2,q**2),GAPMul(3,ER(3)**2)-GAPMul(ER(3)**2,q**-1)-GAPMul(4,ER(3)**2)+GAPMul(3,ER(3)**2)-GAPMul(ER(3)**2,q**3),0,0,0,ER(3)-GAPMul(ER(3),q)+GAPMul(ER(3),q**2),0,ER(3)-GAPMul(ER(3),q),0,0,0,0,0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,0,ER(3)**2-GAPMul(2,ER(3)**2)+GAPMul(ER(3)**2,q**2),0,0,0,0,0,0,0],[0,0,ER(3)**2-GAPMul(ER(3)**2,q**-1),0,1-GAPMul(2,q)+q**2,0,0,0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,0,ER(3)**2,0,0,0,0,0,0,0,1-q,0,0,0,0,0,0,0],[0,0,-2+ER(-3)-GAPMul(ER(3),q**-1)+GAPMul(GAPDiv(5-ER(-3),2),q)-q**2,0,1-q,0,0,0,0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,0,0,0,GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0,1-q,0,0],[0,0,-2+ER(-3)-GAPMul(ER(3),q**-1)+GAPMul(GAPDiv(5-ER(-3),2),q)-q**2,-ER(3)**2+GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0,1-q,0,0],[0,0,ER(-3)-GAPMul(ER(3),q**-1)+GAPMul(ER(3)**2,q),0,ER(3)**2+GAPMul(2+ER(-3),q),0,0,0,0,0,0,ER(3)-GAPMul(ER(3),q)+GAPMul(ER(3),q**2),0,0,ER(3)-GAPMul(ER(3),q),0,0,0,0,-1+q,0,0,-1+GAPMul(2,q)-q**2,0,0,0,0,0,0,0],[0,0,-4-q**-2+GAPMul(3,q**-1)-q**2,ER(-3)-GAPMul(ER(3),q**-1)+GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,-ER(3)+GAPMul(ER(3),q**-1),0,ER(3)-GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0,0,-ER(3)**2+GAPMul(ER(3)**2,q),2-q**-1-q,0,0],[3-q**-1-GAPMul(3,q)+q**2,-1+GAPMul(2,q)-q**2,0,0,-1+GAPMul(3,q)-GAPMul(3,q**2)+q**3,1-GAPMul(2,q)+q**2,GAPMul(ER(3)**2,q)-GAPMul(2,ER(3)**2)+GAPMul(ER(3)**2,q**3),ER(3)-GAPMul(ER(-3),q)-GAPMul(ER(3)**2,q**2),0,-1+q,0,ER(3)-GAPMul(ER(-3),q)-GAPMul(ER(3)**2,q**2),0,0,-ER(3)**2+GAPMul(ER(3)**2,q),0,-1,0,0,0,0,0,-1+GAPMul(2,q)-q**2,0,1-q,GAPMul(ER(3)**2,q)-GAPMul(ER(3)**2,q**2),0,0,0,0],[0,0,GAPMul(-4,ER(3))-GAPMul(ER(3),q**-2)+GAPMul(3,ER(3))-GAPMul(ER(3),q**2),0,GAPDiv(-3+ER(-3),2)+q**-1-GAPMul(ER(3),q),0,0,0,0,0,0,0,-ER(3)**2+GAPMul(ER(3)**2,q**-1),0,0,0,0,ER(3)-GAPMul(ER(3),q),0,0,0,0,0,0,0,0,-1+q,GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),0,0],[-2+q**-1,0,0,-3+q**-1-q**2,1-GAPMul(3,ER(-3))+GAPMul(ER(-3),q**-1)-GAPMul(4,ER(3))+GAPMul(ER(3),q**3),-1+GAPMul(2,ER(-3))-GAPMul(ER(-3),q**-1)-GAPMul(3,ER(3))+GAPMul(ER(3),q**2),GAPDiv(-3+ER(-3),2)+GAPMul(GAPDiv(7-ER(-3),2),q)-GAPMul(3,q**2)+q**3,GAPMul(4,ER(3)**2)-GAPMul(ER(3)**2,q**-1)+GAPMul(2+GAPMul(3,ER(-3)),q)-GAPMul(2,ER(-3))+GAPMul(ER(3),q**3),0,GAPDiv(5+ER(-3),2)-q**-1+GAPMul(-2-ER(-3),q)-GAPMul(ER(3)**2,q**2),0,-2-ER(-3)-GAPMul(ER(3)**2,q**-1)+GAPMul(GAPDiv(5+ER(-3),2),q)-q**2,GAPDiv(-3+ER(-3),2)+q**-1-GAPMul(ER(3),q),3-q**-1-GAPMul(3,q)+q**2,-2+q**-1,-1+GAPMul(2,q)-q**2,0,-1+q,-1,0,-ER(3)**2-GAPMul(ER(-3),q)+GAPMul(ER(3),q**2),0,GAPMul(-3,ER(3))+GAPMul(ER(3),q**-1)-GAPMul(ER(3),q**2),0,GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),-1+GAPMul(2,q)-q**2,1-GAPMul(2,q)+q**2,GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**-1)-GAPMul(ER(3)**2,q),0,-ER(3)+GAPMul(ER(3),q)],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3),q),0,0,0,0,0,0,0],[-4+GAPMul(2,q**-1)-q**2,1-GAPMul(2,q)+q**2,0,-2+q**-1,0,0,0,-1-GAPMul(2,ER(-3))-GAPMul(ER(3)**2,q**-1)+GAPMul(GAPDiv(1+GAPMul(5,ER(-3)),2),q)-GAPMul(ER(-3),q**2),0,GAPDiv(5+ER(-3),2)-q**-1+GAPMul(GAPDiv(-3-ER(-3),2),q),0,0,0,2-q**-1-q,0,-1+q,0,0,0,0,ER(3)-GAPMul(ER(3),q),0,0,0,0,0,1-q,0,0,-ER(3)],[0,0,0,0,2-ER(-3)+GAPMul(ER(3),q**-1)-q**3,0,0,0,0,0,0,-1-GAPMul(2,ER(-3))-GAPMul(ER(3)**2,q**-1)+GAPMul(GAPDiv(1+GAPMul(5,ER(-3)),2),q)-GAPMul(ER(-3),q**2),-2+q**-1,0,GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**-1)-GAPMul(ER(3)**2,q),0,0,-1+q,0,-1+q**-1,0,ER(3)-GAPMul(ER(3),q),1-GAPMul(2,q)+q**2,0,0,0,-1+q,0,-ER(3)**2,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2,0,0,-1+q,0,0,0,0,0,0,0],[5+GAPMul(2,ER(-3)),ER(3)**2-GAPMul(3,ER(3)**2)+GAPMul(3,ER(3)**2)-GAPMul(ER(3)**2,q**3),0,4-q**-1-GAPMul(6,q)+GAPMul(4,q**2)-q**3,GAPDiv(-3+GAPMul(7,ER(-3)),2)-GAPMul(ER(-3),q**-1)+GAPMul(GAPDiv(11-GAPMul(9,ER(-3)),2),q)-q**4,GAPDiv(3-GAPMul(5,ER(-3)),2)+GAPMul(ER(-3),q**-1)-q**3,GAPDiv(3-ER(-3),2)+GAPMul(GAPDiv(-9+ER(-3),2),q)-GAPMul(ER(3)**2,q**4),GAPDiv(3+GAPMul(5,ER(-3)),2)+GAPMul(ER(3)**2,q**-1),1+GAPMul(GAPDiv(-3+ER(-3),2),q)-GAPMul(ER(3),q**2),-4-ER(-3)+q**-1-GAPMul(ER(3)**2,q**3),-ER(3)+GAPMul(ER(3),q),GAPDiv(3+GAPMul(5,ER(-3)),2)+GAPMul(ER(3)**2,q**-1)-GAPMul(ER(-3),q**3),3-q**-1-GAPMul(3,q)+q**2,-3+q**-1,GAPMul(-3,ER(3)**2)+GAPMul(ER(3)**2,q**-1)-GAPMul(ER(3)**2,q**2),1+GAPMul(GAPDiv(-5-ER(-3),2),q)-q**3,0,1-GAPMul(2,q)+q**2,0,2-q**-1-q,ER(3)**2+GAPMul(GAPDiv(1+GAPMul(3,ER(-3)),2),q),ER(3)**2+GAPMul(ER(-3),q)-GAPMul(ER(3),q**2),GAPDiv(-3+ER(-3),2)+GAPMul(4-ER(-3),q),-1,2-ER(-3)+GAPMul(ER(3),q**-1),1+GAPMul(GAPDiv(-5-ER(-3),2),q),q-GAPMul(2,q**2)+q**3,-ER(3)**2+GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**2),ER(3)**2-GAPMul(ER(3)**2,q),ER(3)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2)],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3),-1+q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3),q),-1+q,0,0],[0,0,-4-q**-2+GAPMul(3,q**-1)-GAPMul(3,q**2)+q**3,0,-ER(3)+GAPMul(2,ER(3))-GAPMul(ER(3),q**2),0,0,0,0,0,0,GAPMul(3,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(4,ER(3))+GAPMul(3,ER(3))-GAPMul(ER(3),q**3),2-q**-1-GAPMul(2,q)+q**2,0,GAPMul(2,ER(3))-GAPMul(ER(3),q**-1)-GAPMul(ER(3),q),0,0,1-GAPMul(2,q)+q**2,0,0,0,-ER(3)+GAPMul(ER(3),q)-GAPMul(ER(3),q**2),ER(3)**2-GAPMul(ER(3)**2,q),0,0,0,-ER(3)**2+GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**2),2-q**-1-GAPMul(2,q)+q**2,ER(3)**2-GAPMul(ER(3)**2,q),0],[-5+GAPMul(2,q**-1)-GAPMul(3,q**2)+q**3,2-GAPMul(4,q)+GAPMul(3,q**2)-q**3,1-GAPMul(2,q)+q**2,-ER(-3)-GAPMul(ER(3)**2,q**-1)+GAPMul(GAPDiv(-5+ER(-3),2),q)-q**3,0,0,0,-1-GAPMul(2,ER(-3))-GAPMul(ER(3)**2,q**-1)+GAPMul(1+GAPMul(3,ER(-3)),q),0,GAPDiv(7+ER(-3),2)-q**-1+GAPMul(-4-ER(-3),q),0,0,0,2-q**-1-GAPMul(2,q)+q**2,0,-1+GAPMul(2,q)-q**2,0,0,0,0,-ER(3)**2+GAPMul(ER(3)**2,q)-GAPMul(ER(3)**2,q**2),0,0,0,0,0,0,ER(3)-GAPMul(ER(3),q),0,ER(3)**2-GAPMul(ER(3)**2,q)]],[[0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-GAPMul(2,q)+q**2,q-q**2,0,0,-ER(3)+GAPMul(ER(3),q),-1+q,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1+q,0,0,0,ER(3)-GAPMul(ER(3),q),1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,ER(3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,q,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,q,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,ER(3),0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,q-q**2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,q,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,1-q,0,0,1-q,0,0,0,0,0,0,0,0,0,0,-ER(3),0],[0,0,0,0,-q+q**2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0],[-2+ER(-3)-GAPMul(ER(3),q**-1)+GAPMul(4-ER(-3),q),-ER(3)+GAPMul(-2+ER(-3),q)-q**3,0,-ER(3)+GAPMul(2,ER(3))-GAPMul(ER(3),q**2),1+GAPMul(2,ER(-3)),0,-ER(3)+GAPMul(-2+ER(-3),q),0,-ER(3)+GAPMul(2,ER(3))-GAPMul(ER(3),q**2),0,-ER(3)+GAPMul(ER(3),q),0,0,0,0,0,0,0,0,-ER(3)+GAPMul(ER(3),q),0,0,0,0,-1+q,GAPMul(ER(3),q),0,0,0,0],[1-GAPMul(2,q)+q**2,1-q,0,1-GAPMul(2,q)+q**2,1-GAPMul(2,q)+q**2,GAPDiv(-5-ER(-3),2)+q**-1,0,0,1-q,0,0,0,0,0,0,0,0,0,0,2-q**-1-q,0,0,1-q,0,ER(3)**2,0,0,0,0,0],[0,0,-ER(3)+GAPMul(ER(3),q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,ER(3),0,0],[0,0,q-q**2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0],[0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3)**2,q)-GAPMul(ER(3)**2,q**2),1-q,0,0,0,0,0,0,0,0,GAPMul(-ER(3)**2,q),0,0,0,0,0,0,-1+q,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0,0,0,0,0,0,0,-1+q]],[[0,0,0,-q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-q,q,0,1-q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,ER(3)**2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,GAPMul(ER(3),q),-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,-ER(3)+GAPMul(ER(3),q),0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-q,0,0,1-GAPMul(2,q)+q**2,ER(3)-GAPMul(2,ER(3))+GAPMul(ER(3),q**2),GAPDiv(-5-ER(-3),2)+q**-1,-q+q**2,0,1-GAPMul(2,q)+q**2,0,0,0,0,0,0,0,0,0,0,1-q,0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),-q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-ER(3)**2+GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,GAPDiv(3+ER(-3),2)-q**-1+GAPMul(ER(3)**2,q),0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,-ER(3)],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1-q**-1,0,-ER(3)**2,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)-GAPMul(ER(3),q),0,0,GAPMul(-ER(3),q),-1+q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,-1+q,0,0,-q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,1,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q+q**2,0,0,-q+q**2,0,0,0,0,0,-1+q,0,0,0,0,GAPMul(ER(3),q),0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(ER(3),q),0,-1+q,0,0,0,0,0],[0,0,0,-1+q,2-ER(-3)+GAPMul(ER(3),q**-1),0,-1+q,0,-1+q,0,-1,0,0,0,0,0,0,0,0,0,0,0,1-q,0,0,-1+q,0,0,0,0],[0,0,0,0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q**-1),0,-ER(3)+GAPMul(ER(3),q),0,0,0,0,ER(3)**2,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,-ER(3)**2-GAPMul(ER(-3),q)+GAPMul(ER(3),q**2),0,0,0,0,0,GAPMul(-ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0,0,0,0]],[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,1-q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1-q**-1,0,0,1,0,0,0,0,0,0,0],[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,ER(3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,GAPMul(ER(3)**2,q),-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,1-q,-2+ER(-3)-GAPMul(ER(3),q**-1)+GAPMul(GAPDiv(5-ER(-3),2),q)-q**2,0,0,0,1-q,0,1,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,1,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,ER(3)**2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-GAPMul(2,q)+q**2,q-q**2,0,0,-ER(3)+GAPMul(ER(3),q),-1+q,q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,q,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1+q,-q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,GAPMul(ER(3),q),0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-ER(3)**2,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(-ER(3),q),0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1-q,0,-1+q,0,0,0,0,-1,0,0,0,0,0,0],[q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,-1+q,0,0,0,0,0,0,0,0],[1-q,q,0,GAPMul(-ER(3)**2,q)+GAPMul(ER(3)**2,q**2),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q-q**2,0,-q,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,q-q**2,0,0,-q+q**2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0],[0,0,0,-1+GAPMul(2,q)-q**2,GAPDiv(5-GAPMul(3,ER(-3)),2)+GAPMul(ER(3),q**-1)-q**3,0,q-q**2,0,-1+GAPMul(2,q)-q**2,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,1-GAPMul(2,q)+q**2,0,0,q,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0],[0,0,0,0,0,0,0,0,0,q-q**2,0,0,0,0,0,-q+q**2,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,GAPMul(ER(3)**2,q)],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,ER(3),0]],[[-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-1+q,0,GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,ER(3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-2+q**-1,1-q,-ER(3)**2+GAPMul(ER(3)**2,q**-1),0,0,0,0,ER(3)**2-GAPMul(ER(3)**2,q),0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,ER(3),0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[1-GAPMul(2,q)+q**2,q-q**2,0,0,-ER(3)+GAPMul(ER(3),q),-1+q,q,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,-ER(3),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,GAPMul(-ER(3)**2,q),-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,q,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,GAPMul(-ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,-1+q,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-ER(3),0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,0,0,-ER(3),0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,-q+q**2,-ER(3)+GAPMul(ER(3),q),0,0,0,0,0,0,0,0,q,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-ER(3),0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1+q,0,GAPMul(-ER(3)**2,q)+GAPMul(ER(3)**2,q**2),0,0,0,0,q,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(-ER(3)**2,q),0,0,0,0,0,0,0,0,0],[4-GAPMul(2,q**-1)-GAPMul(3,q)+q**2,-1+GAPMul(2,q)-q**2,0,2-q**-1-q,0,0,0,1+GAPMul(2,ER(-3)),0,GAPDiv(-5-ER(-3),2)+q**-1,0,0,0,-2+q**-1,0,1-q,0,0,0,0,-ER(3)+GAPMul(ER(3),q),0,0,0,0,0,-1+q,0,0,ER(3)],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,-1+q,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,GAPMul(-ER(3)**2,q),0,0,0,0,-1+q,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,-ER(3)**2+GAPMul(ER(3)**2,q),0,0,-ER(3)**2+GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,0,0,q,0],[GAPDiv(-1-GAPMul(3,ER(-3)),2)-GAPMul(ER(3)**2,q**-1)+GAPMul(GAPDiv(-1+GAPMul(5,ER(-3)),2),q),-ER(3)**2+GAPMul(GAPDiv(-1-GAPMul(3,ER(-3)),2),q)-GAPMul(ER(3),q**3),0,-ER(3)**2+GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**2),GAPDiv(-7-ER(-3),2)+q**-1,0,-ER(3)**2+GAPMul(GAPDiv(-1-GAPMul(3,ER(-3)),2),q),0,-ER(3)**2+GAPMul(2,ER(3)**2)-GAPMul(ER(3)**2,q**2),0,-ER(3)**2+GAPMul(ER(3)**2,q),0,0,0,0,0,0,0,0,-ER(3)**2+GAPMul(ER(3)**2,q),0,0,0,0,-ER(3)+GAPMul(ER(3),q),GAPMul(ER(3)**2,q),0,0,0,-1+q]]])
              
              
              r=[f1(x),f2(x),f3(x,ER(3)),f4(x,ER(3)),f4(x,ER(3)**2),f3(x,ER(3)**2),[[[-1]],[[-1]],[[-1]],[[-1]],[[-1]]],f8(x),f9(x),GAPMul(f8(GAPDiv(1,x)),-x),f11(x,ER(3)),f12(x,ER(3)),f13(x,ER(3)),f13(x,ER(3)**2),f12(x,ER(3)**2),f11(x,ER(3)**2),f17(x),GAPMul(f1(GAPDiv(1,x)),-x),GAPMul(f13(GAPDiv(1,x),ER(3)**2),-x),f20(x,ER(3)),f20(x,ER(3)**2),GAPMul(f13(GAPDiv(1,x),ER(3)),-x),f23(x),GAPMul(f2(GAPDiv(1,x)),-x),GAPMul(f11(GAPDiv(1,x),ER(3)**2),-x),GAPMul(f12(GAPDiv(1,x),ER(3)),-x),GAPMul(f12(GAPDiv(1,x),ER(3)**2),-x),GAPMul(f11(GAPDiv(1,x),ER(3)),-x),f29(x),GAPMul(f4(GAPDiv(1,x),ER(3)**2),-x),GAPMul(f4(GAPDiv(1,x),ER(3)),-x),GAPMul(f23(GAPDiv(1,x)),-x),GAPMul(f3(GAPDiv(1,x),ER(3)**2),-x),GAPMul(f3(GAPDiv(1,x),ER(3)),-x),GAPMul(f17(GAPDiv(1,x)),-x),[[[x]],[[x]],[[x]],[[x]],[[x]]]]
              return r[i-1]
            
            
            return m(i)
          else:
            if [p,q,r]==[4,4,3] :
              x=GAPDiv(-para[2][1-1],para[2][2-1])
              r=GAPMul(x**0,[[[[x,-1,-1,0,0,0],[0,-1,0,0,0,0],[0,0,-1,0,0,0],[0,0,0,-1,0,0],[0,0,-1+x,-1,x,0],[0,1,-1,-1,0,x]],[[-1,0,0,0,0,0],[-x,x,0,0,0,x],[0,0,0,0,-x,0],[0,0,0,x,0,0],[0,0,-1,0,-1+x,0],[0,0,0,0,0,-1]],[[x,-1,0,0,0,-1],[0,-1,0,0,0,0],[0,0,x,0,1,-1],[0,-x,0,x,-1,1-x],[0,0,0,0,-1,0],[0,0,0,0,0,-1]]],[[[-1,0,0],[0,0,1],[0,x,-1+x]],[[x,0,0],[-1,-1,0],[1,0,-1]],[[0,1,0],[x,-1+x,0],[0,0,-1]]],[[[x,0,0],[x,-1,0],[x,0,-1]],[[-1,2,0],[0,x,0],[0,GAPMul(-ER(4)+1,x),-1]],[[-1,0,1],[0,-1,GAPDiv(ER(4)+1,2)],[0,0,x]]],[[[x,0,0],[x,-1,0],[x,0,-1]],[[-1,2,0],[0,x,0],[0,GAPMul(ER(4)+1,x),-1]],[[-1,0,1],[0,-1,GAPDiv(-ER(4)+1,2)],[0,0,x]]],[[[-1]],[[-1]],[[-1]]],[[[x,-1,0],[0,-1,0],[0,-1,x]],[[-1+x,0,1],[0,x,0],[x,0,0]],[[0,x,0],[1,-1+x,0],[0,0,x]]],[[[-1,0,0],[-1,x,0],[-1,0,x]],[[x,GAPMul(-2,x),0],[0,-1,0],[0,-ER(4)-1,x]],[[x,0,-x],[0,x,GAPMul(GAPDiv(ER(4)-1,2),x)],[0,0,-1]]],[[[-1,0,0],[-1,x,0],[-1,0,x]],[[x,GAPMul(-2,x),0],[0,-1,0],[0,ER(4)-1,x]],[[x,0,-x],[0,x,GAPMul(GAPDiv(-ER(4)-1,2),x)],[0,0,-1]]],[[[-1,0],[-1,x]],[[-1,0],[-1,x]],[[x,-x],[0,-1]]],[[[x]],[[x]],[[x]]]])
              return r[i-1]
            else:
              S=CHEVIE.imp.CharInfo(p, q, r).charparams[i-1]
              def p1rRep():
                if r>1 :
                  Q=GAPDiv(-para[2][1-1],para[2][2-1])
                else:
                  Q=0
                def pos(t,i):
                  for j in range(1,len(t)+1):
                    for k in range(1,len(t[j-1])+1):
                      l=t[j][k-1].index(i)+1
                      if !=(l,false) :
                        return [j,k,l]
                
                
                ct=lambda p: GAPMul(para[1][p[1-1]-1],Q**p[3-1]-p[2-1])
                T=Tableaux(S)
                return GAPMul(Concatenation([DiagonalMat(map(lambda S: ct(pos(S,1)),T))],map(lambda i: map(ximp37,range(1,len(T)+1)),range(2,r+1))),prod(para)**0)
              
              
              if q==1 :
                return p1rRep()
              else:
                if p==q :
                  para=[map(lambda i: ER(p)**i,range(0,p-1+1)),para[1-1]]
                else:
                  if !=(para[2-1],para[3-1]) :
                    if q%2==0 and r==2 :
                      S=CHEVIE.imp.CharInfo(p, q, r).malle[i-1]
                      if S[1-1]==1 :
                        return [[[para[1][1+S[4-1]-1%GAPDiv(p,q)-1]]],[[para[2][S[2-1]-1]]],[[para[3][S[3-1]-1]]]]
                      else:
                        Y=para[2-1]
                        T=para[3-1]
                        if q>2 :
                          X=map(lambda y: GetRoot(y,GAPDiv(q,2)),para[1-1])
                          X=Concatenation(map(lambda i: GAPMul(ER(GAPDiv(q,2))**i,X),range(1,GAPDiv(q,2)+1)))
                        else:
                          X=para[1-1]
                        X=[X[k-1] for k in [S[k-1] for k in [3,4]]]
                        v=GAPMul(S[2-1],GetRoot(GAPMul(prod(X),prod(Y)),2))
                        d=1+GAPMul(Sum(X),0)
                        return [GAPMul(d,[[X[1-1],Sum(Y,lambda y: GAPDiv(1,y))-GAPMul(GAPDiv(X[2-1],v),Sum(T))],[0,X[2-1]]])**GAPDiv(q,2),[[Sum(Y),GAPDiv(1,X[1-1])],[GAPMul(-prod(Y),X[1-1]),0]],[[0,GAPDiv(-prod(T),v)],[v,Sum(T)]]]
                    else:
                      Error("should  !  happen")
                  else:
                    if para[1-1]==map(lambda i: ER(GAPDiv(p,q))**i-1,range(1,GAPDiv(p,q)+1)) :
                      para=[map(lambda i: ER(p)**i,range(0,p-1+1)),para[2-1]]
                    else:
                      para=[Concatenation(Matrix(map(lambda i: GAPMul(map(lambda j: ER(q)**j,range(0,q-1+1)),GetRoot(i,q)),para[1-1])).transpose()),para[2-1]]
              extra=false
              if IsInt(S[len(S)-1]) :
                extra=ER(S[len(S)-1-1])**S[len(S)-1]
                d=len(S)-2
                S=FullSymbol(S)
              v=p1rRep()
              if p==q :
                v=Concatenation([v[2-1]**v[1-1]],[v[k-1] for k in range(2,len(v)+1)])
              else:
                if q>1 :
                  v=Concatenation([v[1-1]**q,v[2-1]**v[1-1]],[v[k-1] for k in range(2,len(v)+1)])
              if !=(extra,false) :
                m=PermListList(T,map(lambda S: [S[k-1] for k in Concatenation(range(d+1,p+1),range(1,d+1))],T))
                m=Cycles(m,range(1,len(T)+1))
                l=map(lambda i: extra**i,[0,..(-1,1-GAPDiv(p,d))])
                m1=map(lambda x: x[1-1],m)
                return map(lambda x: map(lambda c: GAPMul(l,[x{c}[k-1] for k in m1]),m),v)
              else:
                return v

def ximp37(j):
  S=T[j-1]
  a=pos(S,i)
  b=pos(S,i-1)
  S=map(lambda a: map(ShallowCopy,a),S)
  (S[a[1]])[a[2]][a[3-1]-1]=i-1
  (S[b[1]])[b[2]][b[3-1]-1]=i
  if para[2][1-1]==-para[2][2-1] :
    if a[1-1]==b[1-1] :
      tll=GAPDiv(para[2][1-1],a[3-1]+b[2-1]-a[2-1]-b[3-1])
    else:
      tll=0
  else:
    tll=GAPDiv(Sum(para[2-1]),1-GAPDiv(ct(b),ct(a)))
  v=GAPMul(range(1,len(T)+1),0)
  v[j-1]=tll
  p=T.index(S)+1
  if !=(p,false) :
    v[p-1]=tll-para[2][2-1]
  return v

ChevieData["imp"]["HeckeRepresentation"]=ximp36

def ximp38(p,q,r,i):
  o=ChevieData["imp"]["EigenvaluesGeneratingReflections"](p,q,r)
  o=map(Denominator,o)
  return ChevieData["imp"]["HeckeRepresentation"](p,q,r,map(lambda x: map(lambda i: ER(x)**i,range(0,x-1+1)),o),[],i)

ChevieData["imp"]["Representation"]=ximp38

def ximp39(p,q,r):
  o=ChevieData["imp"]["EigenvaluesGeneratingReflections"](p,q,r)
  o=map(Denominator,o)
  return ChevieData["imp"]["HeckeCharTable"](p,q,r,map(lambda x: map(lambda i: ER(x)**i,range(0,x-1+1)),o),[])

ChevieData["imp"]["CharTable"]=ximp39

