#!/usr/bin/ruby
require 'tempfile'
require 'getoptlong'

def help
  print <<END
  Usage:
    maketest.rb [options] [pattern]
    test GAP3 manual files matching [pattern] --- all if no pattern
    options:
    -r         pass option -r to gap (that is, bypass .gaprc)
    -d xxx     use xxx as diff tool
END
end

$opt_r=false
$difftool="gvim -d -geometry 150x50"

begin
GetoptLong.new(
 [ "-r", GetoptLong::NO_ARGUMENT],
 [ "-d", GetoptLong::REQUIRED_ARGUMENT],
).each {|opt,arg|
  if opt=="--no gaprc" then $opt_r=true end
  if opt=="-d" then $difftool=arg end
}
rescue
  help();exit
end

SIZE_SCREEN="[72,70]" # gap manual examples should have this
files=[]
m=File.open("manual.tex","r"){|f| f.each{|l|
     files<<$~[1]+".tex" if /^[^%].*Include\{([^}]*)\}/=~l}}
files=files-["aboutgap.tex"]
files=files.grep(Regexp.new("^#{ARGV[0]}")) if ARGV.length>0 
p files
# test separately 0 (aboutgap) since it destroys lib

# outman=Tempfile.new("outman")
File.delete("outman") if File.exists?("outman"); outman=File.open("outman","w")
File.delete("err") if File.exists?("err"); err=File.open("err","w")
File.delete("script") if File.exists?("script"); spt=File.open("script","w")
spt.print "LogTo(\"out2\");SizeScreen("+SIZE_SCREEN+");;\n"
spt.print "InfoRead1:=Ignore;;InfoChevie:=Ignore;;InfoAlgebra:=Ignore;;\n"
#spt.print "CHEVIE.PrintSpets:=rec(GAP:=true);;\n"
files.each do |n| s=File.read(n)
  spt.print "###   "+n+"   ###\n"
  err.print "###   "+n+"   ###\n"
  outman.print "gap> ###   "+n+"   ###\n"
# s=s.gsub(/%.*$/,"")
  indent=4
  s=s.gsub(/\|'\\\|'\|/,"\\|").gsub(/\|\\#\|/,"#")
  s.scan(/[^\\]\|(([^|]|(\\\|))*?[^\\])\|/m){ |e| 
   sm=Regexp.last_match
   if /^\s*gap>/=~ e[0] 
     cnt=1
     e[0].split("\n",-1).each {|l| cnt+=1
       if /^\s*gap>/=~ l or /^\s*>\s/=~ l
	 l=l.split(/("(?:[^"]|\\")*")/)
	 (0...l.length).each{|i| 
	   if /(#|&)/=~l[i] and not /^"/=~l[i] then
	    l[i]=l[i].sub(/ *(#|&).*/,"")
	    l=l[0..i]
	   end
	 }
	 l=l.join("")
       end
       if /^\s*gap>/=~ l and not /quit;/=~ l
         spt.print l.sub(/\s*gap> /,""), "\n"
       elsif /^\s*>\s/=~ l
         spt.print l.sub(/\s*>\s/,""), "\n"
       end
       if /^\s*brk>/=~l
            outman.print "brk> "
       else
         l=l.sub(/&/,"#") if /^\s*&/=~l
         l=l.gsub(/\\\|/,"|")
         m=/^(\s*)gap>/.match(l)
         if m then
           if m[1].length!=4 then 
             err.print "!!! indent=#{m[1].length} on line ",
                  cnt+sm.pre_match.count("\n"),"\n"
           end
           indent=m[1].length
           l=l[indent..-1]
         else
           m=/^(\s*)/.match(l)
           if m[1].length<indent and not(m[1].length==0 and l.length==0) then
             err.print "*** indent=#{indent} but beg=#{m[1].inspect} on line ",
                  cnt+sm.pre_match.count("\n"),"\n"
           end
           l=l[indent..-1]
         end
         outman.print l,"\n"
       end
     }
#    spt.print "################\n"
   else
     err.print "|",e[0],"|\n"
   end
  }
  print n,"=>",s.length.to_s + "\n"
end
spt.close
outman.close
if $opt_r then
system "gap3 -r < "+spt.path
else
system "gap3 < "+spt.path
end
# out2=Tempfile.new("out2")
File.delete("outgap") if File.exists?("outgap"); outgap=File.open("outgap","w")
File.open("out2"){ |v| v.each {|l|
  l=l.sub(/\s*\n$/,"\n")
# outgap.print l unless /^#I/=~l
  outgap.print l unless /gap> InfoRead1/=~l
}}
File.delete("out2")
outgap.close
system $difftool+" "+outman.path+" "+outgap.path
