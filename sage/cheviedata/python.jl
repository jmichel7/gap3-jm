"""
This  file is an attempt to automatic  translate the GAP3 Chevie data files
to Python.

It  currently fails through limitations of  python. 

For  example, because it  is impossible to  generate a multi-line anonymous
function in Python, the parser treats the anonymous functions by putting on
a  stack function definitions encountered and generating the functions from
the  stack afterwards;  the exception  is the  -> form  which has  the same
limitations as Python's lambda, so is translated to a lambda. Thus,

f=t->function(a)return a+t;end;

is currently translated as

def xxx1(a):
  a+t

f=lambda t: xxx1

while it should be translated as

def f(t):
  def xxx1(a):
    a+t
  return xxx1

but  it is  very hard  to imagine  a general  scheme where this would occur
naturally.

Another  problem is to translate the operators * + - / which in GAP can act
on  lists as well as numbers.  A scheme which works well  for * and / is to
generate  functions GAPMul and GAPDiv  for them. But an  attempt to use the
same scheme for - by translating it to a function GAPMin provokes a 'parser
stack overflow' error in Python.

To use the parser, call to_py(f) where f is a CHEVIE filename.
"""
function lex(s::String) # massage GAP source so julia parser does not choke

s=readstring(s)

while true
  new=replace(s,r",\s*,"s,",missing,")
  if new==s break end
  s=new
end

while true
  new=replace(s,r"\[\s*,"s,"[missing,")
  if new==s break end
  s=new
end

while true
  new=replace(s,r",\s*\]"s,",missing]")
  if new==s break end
  s=new
end

while true #suppress comments makes some below changes harmless
  new=replace(s,r"^([^\"]*)#.*$"m,s"\1") 
  if new==s break end
  s=new
end

s=replace(s,r"([^<>:])=",s"\1==") 
s=replace(s,r":=","=")
s=replace(s,r"<>","!=")
s=replace(s,r"return\s*\n"," return")
s=replace(s,r";\s*fi\s*;"s," end;")
s=replace(s,r";\s*fi\s*"s," end;")
s=replace(s,r"\s*fi\s*;"s," end;")
s=replace(s,r";\s*od\s*;"s," end;")
s=replace(s,r";\s*od\s*"s," end;")
s=replace(s,r"\s*od\s*;"s," end;")
s=replace(s,r"\bthen\b"," ")
s=replace(s,r"\btype\b","type_")
s=replace(s,r"\bdo\b"," ")
s=replace(s,r"\bnot\b"," ! ")
s=replace(s,r"\band\b"," && ")
s=replace(s,r"\bor\b"," || ")
s=replace(s,r"\bmod\b"," % ")
s=replace(s,r"local\s*\w*\s*;"," ")
s=replace(s,r"local(\s*\w*\s*,)*\s*\w*\s*;"," ")
s=replace(s,r";","\n")
s=replace(s,r"elif","elseif")
s=replace(s,r"([^\w\"]\s*)((\s*\((\s*\d+\s*,)+\s*\d+\))+)",s"\1perm\"\2\"")
s=replace(s,r"\)\s*\[",")[")
s=replace(s,r"\n/","/\n")
s=replace(s,r"(?<!,)\n\s*\-","-\n")
s=replace(s,r"(?<!,)\n\s*\+","+\n")
s=replace(s,r"(?<!,)\n\s*\^","^\n")
s=replace(s,r"(?<!,)\n\s*\*","*\n")
s=replace(s,r"\n\s*\&\&","&&\n")
s=replace(s,r", *\-\n",",\n-")
s=replace(s,r"\[ *\-\n","[\n-")
s=replace(s,r"\^\-\n","^\n-")
s=replace(s,r"\)\s*\{","){")
s=replace(s,r"\)\s*\(",")(")
s=replace(s,r"\}\s*\{","}{")
s=replace(s,r"(\w)\s*\(",s"\1(")
s=replace(s,r"(\w+)\.(\d+)",s"\1[:\2]")
end

function myparse(s::String) # use julia parser
  s=lex(s)
  i=1
  l=[]
  while i<length(s)
    p=1
  try
    p=parse(s,i)
  catch
    write("test",s[i:end])
    rethrow()
  end
  # print(p)
    i=p[2]
    push!(l,p[1])
  end

  l
end

function_stack=[]
function_prefix=""
function_counter=0
indent=0

function cr(i=0)
  global indent
  indent+=i
  "\n"*" "^(4*indent)
end

function py(e) # try to output python code
  global function_stack, function_prefix, function_counter
  res=""
  if e isa Integer
    res= string(e)
  elseif e isa String
    res= repr(e)
  elseif e isa Char
    res= repr(e)
  elseif e isa Symbol
    if e == :missing res= "None"
    elseif e== :CHEVIE res="ChevieData"
    elseif e== :class res="class_"
    elseif e== :lambda res="lambda_"
    elseif e== :String res="str"
    elseif e== :ShallowCopy res="copy"
    elseif e== :ComplexReflectionGroup res="ReflectionGroup"
    elseif e== :IsFunc res="callable"
    else res= string(e)
    end
  elseif e.head== :.
    if e.args[2] isa Expr field=py(e.args[2])
    else                  field=repr(string(e.args[2].value))
    end
    if e.args[1]== :CHEVIE res= "ChevieData[$field]"
    else                   res= py(e.args[1])*"[$field]"
    end
  elseif e.head== :call 
    if e.args[1] == :(CHEVIE.AddData)
       res= "ChevieData["*py(e.args[3])*"]["*py(e.args[2])*"]="*py(e.args[4])
    elseif e.args[1] == :(CHEVIE.R)
       res= "ChevieData["*py(e.args[3])*"]["*py(e.args[2])*"]"
    elseif e.args[1] == :(CHEVIE.IndirectAddData)
       res= "for t in "*py(e.args[3])*":"*cr(1)*
         "ChevieData[t]["*py(e.args[2])*"]="*py(e.args[4])
       cr(-1)
    elseif e.args[1] in [:+, :(==), :%, :>, :<, :!=, :>=, :<=]
      if length(e.args)==2
        res= "$(e.args[1])$(py(e.args[2]))"
      else
        res= "$(py(e.args[2]))$(e.args[1])$(py(e.args[3]))"
      end
    elseif e.args[1]== :in
        res= py(e.args[2])*" in "*py(e.args[3])
    elseif e.args[1]== :-
      if length(e.args)==2
        res= "-$(py(e.args[2]))"
      elseif e.args[3] isa Expr && e.args[3].head==:vect
        res="["*py(e.args[2])*"-i for i in "*py(e.args[3])*"]"
      elseif e.args[2] isa Expr && e.args[2].head==:vect
        res="[i-"*py(e.args[3])*" for i in "*py(e.args[2])*"]"
      else
        res= py(e.args[2])*"-"*py(e.args[3])
      end
    elseif e.args[1]== :*
      res= "GAPMul("*py(e.args[2])*","*py(e.args[3])*")"
    elseif e.args[1]== :/
      res= "GAPDiv("*py(e.args[2])*","*py(e.args[3])*")"
    elseif e.args[1]== :^
      res= "$(py(e.args[2]))**$(py(e.args[3]))"
    elseif e.args[1]== :!
      res= "not "*py(e.args[2])
    elseif e.args[1]== :Product
       res= "prod($(py(e.args[2])))"
    elseif e.args[1]== :Unbind
       res= "del "*py(e.args[2])
    elseif e.args[1]== :Length
       res= "len($(py(e.args[2])))"
    elseif e.args[1]== :Print
       res= "print "*join(map(py,e.args[2:end]),",")*","
    elseif e.args[1]== :String
       if length(e.args)==2 res="str("*py(e.args[2])*")"
       else res= "just("*join(map(py,e.args[2:end]),",")*")"
       end
    elseif e.args[1]== :TransposedMat
       res= "Matrix("*py(e.args[2])*").transpose()"
    elseif e.args[1]== :Append
       res= py(e.args[2])*"+="*py(e.args[3])
    elseif e.args[1]== :Inherit
       res= py(e.args[2])*".update("*py(e.args[3])*")"
    elseif e.args[1]== :Add
       res= py(e.args[2])*".append("*py(e.args[3])*")"
    elseif e.args[1]== :Position
       res= py(e.args[2])*".index("*py(e.args[3])*")+1"
    elseif e.args[1]== :ApplyFunc
       res= py(e.args[2])*"(*"*py(e.args[3])*")"
    elseif e.args[1]== :List
       res= "map("*py(e.args[3])*","*py(e.args[2])*")"
    elseif e.args[1]== :Group
       res= "PermutationGroup(["*join(map(py,e.args[2:end]),",")*"])"
    elseif e.args[1] in [:E, :ER]
       res= "ER($(py(e.args[2])))"
    elseif e.args[1]== :rec
      inner=map(f->repr(string(f.args[1]))*":"*py(f.args[2]),e.args[2:end])
      cr(1)
      res= "{"*join(inner,","*cr())*"}"
      cr(-1)
    elseif e.args[1]==:IsBound && e.args[2] isa Expr
      b=e.args[2]
      if b.head==:. 
        if b.args[2] isa QuoteNode
          res=repr(string(b.args[2].value))*" in "*py(b.args[1])
        else
          res=py(b.args[2])*" in "*py(b.args[1])
        end
      elseif b.head==:ref 
        res=py(b.args[1])*"["*py(b.args[2])*"]==None"
      end
    else 
      res=py(e.args[1])*"("*join(map(py,e.args[2:end]),",")*")"
    end
  elseif e.head== :(=)
    if e.args[1] isa Symbol && e.args[2] isa Expr && e.args[2].head== :function
      res= "def "*py(e.args[1])*"("*
         join(map(py,e.args[2].args[1].args),",")*"):"*cr(1)*
         py(e.args[2].args[2])*cr(-1)*cr()
    elseif e.args[1] isa Expr && e.args[1].head== :curly
      res="for i,j in zip("*py(e.args[1].args[2])*","*py(e.args[2])*
        "):"*cr(1)*py(e.args[1].args[1])*"[i-1]=j"
      cr(-1)
    else
      res= "$(py(e.args[1]))=$(py(e.args[2]))"
    end
  elseif e.head== :block
    lines=filter(x->x.head!= :line,e.args)
    res=join(map(py,lines), cr())
  elseif e.head== :vect
    if length(e.args)>=1 && e.args[1] isa Expr && e.args[1].head==:call &&
       e.args[1].args[1]==:(..)
      res= "range("*py(e.args[1].args[2])*","*py(e.args[1].args[3])*"+1)"
    elseif length(e.args)>=2 && e.args[2] isa Expr && e.args[2].head==:call &&
       e.args[2].args[1]==:(..)
      res= "range("*py(e.args[1])*","*py(e.args[2].args[3])*"+1,"*
         py(e.args[2].args[2])*"-"*py(e.args[1])*")"
    else
      res= "[$(join(map(py,e.args),','))]"
    end
  elseif e.head== :function
    function_counter+=1
    fname=function_prefix*string(function_counter)
    push!(function_stack,(fname,e))
    res= fname
  elseif e.head== :->
    lines=filter(x->!(x isa Expr) || x.head!= :line,e.args[2].args)
    if length(lines)!= 1 error() end
    res="lambda "*py(e.args[1])*": "*py(lines[1])
  elseif e.head==:ref
    res= py(e.args[1])*"["*py(e.args[2])*"-1]"
  elseif e.head==:macrocall 
    if  e.args[1]==Symbol("@perm_str")
      res= "Permutation("*repr(e.args[2])*")"
    elseif e.args[1]==Symbol("@int128_str")
      res= e.args[2]
    end
  elseif e.head==:curly
    res= "["*py(e.args[1])*"[k-1] for k in "*py(e.args[2])*"]"
  elseif e.head==:return
    res= "return "*py(e.args[1])
  elseif e.head==:if
    if length(e.args)==2
      res="if "*py(e.args[1])*" :"*cr(1)*py(e.args[2])
      cr(-1)
    elseif length(e.args)==3
      res="if "*py(e.args[1])*" :"*cr(1)*py(e.args[2])*cr(-1)*
        "else:"*cr(1)*py(e.args[3])
      cr(-1)
    else error()
    end
  elseif e.head==:tuple
    if length(e.args)==0
      res="Permutation(\"()\")"
    elseif length(e.args)==1
      res=py(e.args[1])
    elseif length(e.args)==2
      res="Permutation(\"(%s,%s)\"%("*py(e.args[1])*","*py(e.args[2])*"))"
    else
      error("not implemented")
    end
  elseif e.head==:for
    res="for "*py(e.args[1].args[1])*" in "*py(e.args[1].args[2])*
      ":"*cr(1)*py(e.args[2])
    cr(-1)
  elseif e.head==:while
    res="while "*py(e.args[1])*":"*cr(1)*py(e.args[2])
    cr(-1)
  elseif e.head== :&&
    res= "$(py(e.args[1])) and $(py(e.args[2]))"
  elseif e.head== :||
    res= "$(py(e.args[1])) or $(py(e.args[2]))"
  else
    println("unimplemented case: head=$(e.head) ",e.args,"")
  end
  return res
end

function python(e) # handle the function stack
  global function_stack
  res=py(e)
  res=join(map(function_stack) do fns
        e=fns[2]
        res1="def $(fns[1])("*join(map(py,e.args[1].args),",")*"):"*cr(1)*
          py(e.args[2])*cr(-1)*cr()
        res1
       end)*res
  resize!(function_stack,0)
  res*cr()
end

function to_py(n)
  global function_prefix, function_counter
  l=myparse("$(homedir())/gap3-dev/pkg/chevie/tbl/"*n*".g")
  n=replace(n,r".*/","")
  function_prefix=n[max(1,end-4):end]
  function_counter=0
  open(n*".py","w")do f
    for e in l 
      write(f,"\n",python(e))
    end
  end
end

ChevieTbl=["cmp4_22", "cmplxg24", "cmplxg25", "cmplxg26", "cmplxg27",
"cmplxg29", "cmplxg31", "cmplxg32", "cmplxg33", "cmplxg34", "coxh3", "coxh4", 
"coxi", "exceptio", "weyla", "weylbc", "weyld", "weyle6", "weyle7", "weyle8", 
"weylf4", "weylg2"]

# cmplximp currently triggers parser error
