class File
def cut(s,c=/,/,max=80,n="  ")
# add \n+n in string s (contains no \n) after pattern c so no line exceeds max
  l=0
  while m=(c.match(s) or /.$/.match(s)) do
    w=m.end(0)
    if l+w>=max
      print "\n"+n
      l=n.length
    end
    l+=w
    print s[0...w]
    s=s[w..-1]
  end
  print "\n"
end
end

File.open("auto.g","wb"){ |out|
  out.print <<END_OF_STRING
CHEVIE.AUTO:=rec();
CHEVIE.AutoLoad:=function(flist,fields,file)local f,t;
  for t in fields do 
    if not IsBound(CHEVIE.AUTO.(t)) then CHEVIE.AUTO.(t):=rec();fi;
    for f in flist do CHEVIE.AUTO.(t).(f):=file;od;
  od;
end;
END_OF_STRING
for file in ARGV do
  p file
  auto=Hash.new
  File.open(file,"rb"){|f| f.read}.
    gsub(/\n\s*#[^\n]*/,""). # suppress commented lines
    gsub(/\s/,"").
  scan( /CHEVIE.(?:Indirect)?AddData
    \(((?:"\w*")|(?:\[[^\]]*\])),((?:"\w*")|(?:\[[^\]]*\])),/x)do |f,s|
    r=s.delete("[]").split(",")
    for f in f.delete("[]").split(",") do
      auto[f]=auto.fetch(f,[])+r
    end
  end
  a=Hash.new
  for f,r in auto do
    r=r.sort.join(",")
    a[r]=a.fetch(r,"")+f+","
  end
  for r,f in a do
    out.cut %{CHEVIE.AutoLoad([#{f.chop}],[#{r}],"#{file[0..-3]}");}
  end
end}
