#!/usr/bin/ruby -w
#
# Script to convert the GAP manual to HTML
# usage convert.rb [-cs] <doc-directory> [<html-directory>]
#
# Caveats: 
#
#  1. This script assumes that the .toc file is up-to-date with the .tex
#  files and will almost certainly fail horribly if this is not true
#
#  2.  The  output  files  are  CxxxSxxx.htm,  (not  .html) plus index.htm,
#  theindex.htm  and biblio.htm. Not  all servers will  serve .htm files as
#  HTML without adjustments
#
#  3.   The  script  assumes  that  the   .tex  files  comply  with  GAP
#  conventions,  including  unwritten  ones.  It  tries  to  follow  the
#  behaviour of the on-line browser
#
#
# Global variables 
#
#  $dir  -- the full pathname of the input directory, including a trailing /
#  $odir -- the full pathname of the output directory, including a trailing /
#  $opt_c and $opt_s set by getopts()
#  $chapters -- the chapters data structure
#  $footer -- the trailer put on every page
#  $indexcount -- used within chapters to number the index anchors
#

# Mainly  diagnostic, prints  the chapters  data structure.  Also checks
# that each section has the correct back reference to its chapter

VERSION="11 Mar 2019"
$footer = "<P>\n<address>gap3-jm<br>#{VERSION}</address></body></html>"
# Printed at the bottom of every page

def pad(s)
  s.to_s.rjust(3,"0")
end

# The names of the section and chapter files are determined by this routine
def name2fn (name)
  sec=$sections_by_name[canonize(name)]
  return "badSectionlink.htm#"+name unless sec
  cnum=pad(sec.chapnum)
  snum=pad(sec.secnum)
  if $opt_c
    res="chap#{cnum}.htm"
    res<< "#SECT#{snum}" unless snum=="000"
    return res
  else
    return "c#{cnum}s#{snum}.htm"
  end
end

class String
  def h(c) # return self surrounded by html <c> element
    if c[0]==?/ then "<#{c}>"+self+"<#{c[1..-1]}>"  # complement use
    else "<#{c}>"+self+"</#{c}>"
    end
  end
  def style(c)
    "<span style=\"#{c}\">"+self+"</span>"
  end
  def clean # unescape most specials
    self.gsub(/\\(\^|_|"|:|'|\*|\{|\}| |<|>|#)/){m=$~;m[1]}.gsub(/\\cr/,"<BR>")
  end
end

MathMacs= { "in"=>"&isin;", "wedge"=>"&and;", "vee"=>"&or;", 
	 "cup"=>"&cup;", "bigcup"=>"<big>&cup;</big>",
	 "cap"=>"&cap;", "bigcap"=>"&<big>cap;</big>",
	 "otimes"=>"&otimes;", "bigotimes"=>"<big>&otimes;</big>", 
	 "oplus"=>"&oplus;", "bigoplus"=>"<big>&oplus;</big>", 
	 "backslash"=>" \\ ", "setminus"=>"\\","split"=>":",
	 "langle"=>"&lang;", "lt"=>" &lt;", "rangle"=>"&rang;", "gt"=>" &gt;",
	 "ge"=>" &ge;", "geq"=>" &ge;", 
	 "ne"=>" &ne;", "neq"=>" &ne;", 
	 "leq"=>" &le;", "le"=>" &le;", 
	 "forall"=>"&forall;", "exists"=>"&exist;",
	 "sqrt"=>"&radic;",
	 "cong"=>"&cong;",
	 "equiv"=>"&equiv;",
	 "lceil"=>"&lceil;", "rceil"=>"&rceil;",
	 "lfloor"=>"&lfloor;", "rfloor"=>"&rfloor;",
	 "emptyset"=>"&empty;",
	 "cdots"=>"...", 
	 "prime"=>"&prime;",
	 "underbrace"=>"",
	 "mid"=>"|", 
	 "rightarrow"=>"&rarr;","mapsto"=>"&rarr;",
	 "Rightarrow"=>"&rArr;",
	 "hookrightarrow"=>"&rarr;","to"=>"&rarr;",
	 "left"=>"","right"=>"", "times"=>"&times;",
	 "bullet"=>"&bull;","star"=>"&lowast;",
	 "cdot"=>".", "partial"=>"&part;", "notin"=>"&notin;",
	 "subset"=>"&sub;", "supset"=>"&sup;", 
	 "subseteq"=>"&sube;", "infty"=>"&infin;",
	 "colon"=>":","circ"=>"o","pm"=>"&plusmn;",
	 "dim"=>"dim".h("/I"), "cos"=>"cos".h("/I"),  "sin"=>"sin".h("/I"), 
	 "tan"=>"tan".h("/I"), "ker"=>"ker".h("/I"),
	 "inf"=>"inf".h("/I"), "sup"=>"sup".h("/I"),
	 "gcd"=>"gcd".h("/I"), "det"=>"det".h("/I"), "Stab"=>"Stab".h("/I"),
	 "GL"=>"GL".h("/I"), "SL"=>"SL".h("/I"), "GO"=>"GO".h("/I"),
	 "Aut"=>"Aut".h("/I"),
	 "log"=>"log".h("/I"), "bmod"=>" mod ".h("/I"),"mod"=>" mod ".h("/I"),
	 "ast"=>"*","tilde"=>"<sup>~</sup>","vdots"=>"&sdot;",
	 "Leftrightarrow"=>"&hArr;",
	 "supsetneq"=>"&sup;<sub>&ne;</sub>",
         "ltimes"=>"&#8905;",
         "rtimes"=>"&#8906;"
	 }

AlwaysMacs={
# mathbb
     "C"=>"&#x2102;",
     "F"=>"&#x1D53D;",
     "Q"=>"&#x211A;",
     "R"=>"&#x211D;",
     "Z"=>"&#x2124;",
# mathbf
     "bs"=>"s".h("strong"), "bb"=>"b".h("strong"), "bw"=>"w".h("strong"), 
     "bt"=>"t".h("strong"), 
     "bB"=>"B".h("strong"), "bG"=>"G".h("strong"), "bH"=>"H".h("strong"), 
     "bL"=>"L".h("strong"), "bP"=>"P".h("strong"), "bS"=>"S".h("strong"), 
     "bT"=>"T".h("strong"), "bU"=>"U".h("strong"), "bW"=>"W".h("strong"), 
# fraktur
     "Sym"=>"S".h("strong"),
# mathcal
     "ca"=>"a".style("font-family: cursive"),
     "cA"=>"A".style("font-family: cursive"),
     "cB"=>"B".style("font-family: cursive"),
     "cC"=>"C".style("font-family: cursive"),
     "cH"=>"H".style("font-family: cursive"),
     "gg"=>"g".style("font-family: cursive"),
     "gsl"=>"sl".style("font-family: cursive"),
     "ell"=>"l".style("font-family: cursive"),
# greek
     "alpha"=>"&alpha;", "beta"=>"&beta;",
     "Gamma"=>"&Gamma;","gamma"=>"&gamma;",
     "Delta"=>"&Delta;","delta"=>"&delta;",
     "epsilon"=>"&epsilon;", "varepsilon"=>"&epsilon;",
     "eta"=>"&eta;", "theta"=>"&theta;",
     "zeta"=>"&zeta;", "iota"=>"&iota;",
     "lambda"=>"&lambda;", "mu"=>"&mu;", "nu"=>"&nu;",
     "xi"=>"&xi;", "Pi"=>"&Pi;","pi"=>"&pi;","prod"=>"&prod;",
     "rho"=>"&rho;", "Phi"=>"&Phi;", "phi"=>"&phi;", "varphi"=>"&phi;",
     "chi"=>"&chi;", "sum"=>"&sum;","sigma"=>"&sigma;","Sigma"=>"&Sigma;",
     "tau"=>"&tau;", "Psi"=>"&Psi;", "psi"=>"&psi;",
     "Omega"=>"&Omega;","omega"=>"&omega;",
# package names
     "GAP"=>"GAP3".h("strong"), "CAS"=>"CAS".h("strong"),
     "ATLAS"=>"ATLAS".h("strong"), 
     "CHEVIE"=>"CHEVIE".h("strong"),
     "VKCURVE"=>"VKCURVE".h("strong"), "MAPLE"=>"MAPLE".h("strong"),
     "MeatAxe"=>"MeatAxe".h("strong"), "VE"=>"VE".h("strong"),
     "GUAVA"=>"GUAVA".h("strong"), "GRAPE"=>"GRAPE".h("strong"),
     "SISYPHOS"=>"Sisyphos".h("strong"), "LaTeX"=>"LaTeX".h("strong"),
     "Specht"=>"S"+"PECHT".h("small"),
# punctuation
     "par"=>"<P>", "smallskip"=>"<P>", "bigskip"=>"<P>", "medskip"=>"<P>",
     "hfill"=>"<P align=right>",
     "noindent"=>"",
     "item"=>"<P>&bull;",
     "quad"=>"&emsp;", "qquad"=>"&emsp;&emsp;",
     "ldots"=>"...", "dots"=>"...", "S"=>"&sect;", " "=>"&sp;"}

Accents={"127u"=>"&uuml;","127o"=>"&ouml;","127O"=>"&Ouml;",
    "127e"=>"&euml;","127a"=>"&auml;",
    "18e"=>"&egrave;","19e"=>"&eacute;","19a"=>"&aacute;"}

OneChar={"*"=>"*","{"=>"{","}"=>"}","|"=>"|","'"=>"'",":"=>":"," "=>" ",
  "<"=>"&lt;",">"=>"&gt;",";"=>" ","-"=>"","#"=>"#",","=>"","!"=>"","~"=>" "}

$trads=Hash.new;$key="zaaaaa"

# This  routine is called to process the  text of the section. The output file
# is assumed to be pre-selected. The input file is infile.
# 
# As  we process, we can be in :normal status, :math status (inside $..$)
# "verbatim" status inside a  multi-line example or "shortverb" status
# inside a short |..|
#
# We  separately track whether we are in bold  or tt, whether we are in a xxx:
# ....  paragraph and whether  we are reading  a cross-reference that is split
# across multiple lines
#
# Finally,  we track whether we  have already emitted a  <P> for this group of
# blank lines

# convert TeX text to html
# text assumed to be in state mode (:normal or :math)
def TeX2htm(text,state=:normal,fname="",sec="")
  sp=proc{|s| $trads[$key.succ!]=s; "\\#{$key} "}
  mathdisplay=proc{|s|
   sp[TeX2htm(s,:math).h("I").h("td").h("tr").h("table").h("center")]}
  text=text.gsub(/([^\\]|^)%.*/){m=$~;m[1]} # suppress tex comments
  unless state==:math
    text.gsub!(/\\$/,"\\ ")
    vb=/\|'\\\|'\|/      # a  | in verbatim
    vb2=/\|\s*\{\s*\\tt\s*\\\|\s*\}\s*\|/ # another way for | in verbatim
    vbb=/(?:#{vb2}|#{vb})/o
    text.gsub!(/([^\\]|^|\\\\)\|((?:[^|]|#{vbb})*)\|/om){ m=$~
       if m[2]=~/\n/
        m[1]+sp[m[2].gsub(vbb,"|").lines.map{|l| v=l.split("#");
	v[0].gsub!(/&/,"#");v[0].gsub!(/</,"&lt;");
	if v[1] then v[1]=TeX2htm(v[1]); v.join("#")
	else v[0]
	end}.join.h("pre")]
       else m[1]+sp[m[2].gsub(vbb,"|").gsub(/&/,"#").gsub(/</,"&lt;").clean.h("code")]
       end
    }
    text.gsub!(/\$\$(.*?)\$\$/m){m=$~; mathdisplay[m[1]]}
    text.gsub!(/\\\[(.*?)\\\]/m){m=$~; mathdisplay[m[1]]}
    text.gsub!(/([^\\]|^)\$(.*?)\$/m){m=$~
      m[1]+sp[TeX2htm(m[2],:math).h("I")]}
    text.gsub!(/([^\\]|^)\*([^*]*)([^\\])\*/){m=$~;
      m[1]+sp[(m[2]+m[3]).h("strong")]}
  end
  text.gsub!(/\\\\/){"\\cr"}
  text.gsub!(/([^\\])<([^>]*)>/m){m=$~;m[1]+sp[TeX2htm(m[2]).h("var")]}
  text.gsub!(/([^\\]|^)'((?:[^']|\\')*)([^\\])'/m){m=$~;
    m[1]+sp[(m[2]+m[3]).clean.h("code")]}
  text.gsub!(/([^\\]|^)\\([*{}|': <>;,#~-])/){m=$~;m[1]+sp[OneChar[m[2]]]}
  text.gsub!(/([^\\]|^)\\([*{}|': <>;,#~-])/){m=$~;m[1]+sp[OneChar[m[2]]]}
  unless state==:math
  text.gsub!(/^\s*$/,"\\par")
  text.gsub!(/^\s*\\vspace\{[^}]*\}\s*$/,"\\par")
  text.gsub!(/^([^:\n]*)([^\\]):(.*?)(?=\\par)/m){m=$~
    sp["<DL><DT>#{TeX2htm(m[1]+m[2])}:<DD>"]+"#{m[3]}"+sp["</DL>"]}
  text.gsub!(/\\cr/){sp["<BR>"]}
  # characters that can appear in a cross reference
  text.gsub!(/([^\\]|^)"([\w\s.$-]*)"/m){m=$~;
    m[1]+sp["<a href=\"#{name2fn m[2]}\">#{m[2]}</a>"]}
  text.gsub!(/\\"/){sp['"']}
  text.gsub!(/\\cite(?:\[([^\]\[]*)\])?\s*\{\s*((?:\w|,)+)\s*\}/m) { m=$~
     r="<A href=\"biblio.htm##{m[2].split(",")[0]}\">#{m[2].h("cite")}</a>"
     r<< ", "+m[1] if m[1]
     sp[r] }
  text.gsub!(/\\accent\s*([0-9]*)\s*([a-zA-Z])/m){m=$~
    a=Accents[m[1]+m[2]]
    if a then a
    else $err.print "** unhandled accent #{fname}:#{sec} #{m[0]}\n"
    end}
  text.gsub!(/\\psubsection\s*\{(.*?)\}/m){m=$~; m[1].h("h3")}
  text.gsub!(/\\psection\s*\{(.*?)\}/m){m=$~;m[1].h("h1").style("color:#ff0000")}
  text.gsub!(/\\index\s*\{(.*?)\}/m){m=$~
    $indexcount+=1 # emit an anchor and remember the index keys for later
    r=m[1].gsub(/\\([a-zA-Z]+) /){m=$~; $trads[m[1]] || "\\#{m[1]} "}
    $index[r]=$index[r].push(["#{fname}#I#{$indexcount}","#{sec.chapnum}.#{sec.secnum}"])
    "<A name = \"I#{$indexcount}\"></a>\n"
  }
  text.gsub!(/\\begin\{([^}]*)\}(.*?)\\end\{\1\}/m){m=$~
    case m[1]
    when "displaymath","equation" then mathdisplay[m[2]]
    when "center" then m[2].h("center")
    when "itemize","enumerate" then m[2]
    else "{beg-"+m[1]+"}".h("s")+m[2]+"{end-"+m[1]+"}".h("s")
    end}
  text.gsub!(/~/," ")
  end
#  nestedbraces=// # can define but not recognized inside another regexp
# text.gsub!(/\\(\s)/){m=$~;m[1]}
  text.gsub!(/\{\s*\\(?:sl|it)([^}]*)\}/m){m=$~;
    sp[TeX2htm(m[1],state).h("I")]}
  text.gsub!(/\{\s*\\cal([^}]*\{[^}]*\}[^}]*)\}/m){sp[$~[1].style("font-family: cursive")]}
  text.gsub!(/\{\s*\\cal([^}]*)\}/m){sp[$~[1].style("font-family: cursive")]}
  text.gsub!(/\{\s*\\tt([^}]*\{[^}]*\}[^}]*)\}/m){m=$~;sp[m[1].h("code")]}
  text.gsub!(/\{\s*\\tt([^}]*)\}/m){m=$~;sp[m[1].h("code")]}
  text.gsub!(/\{\s*\\bf([^}]*\{[^}]*\}[^}]*)\}/m){m=$~;sp[m[1].h("strong")]}
  text.gsub!(/\{\s*\\bf([^}]*)\}/m){m=$~;sp[m[1].h("strong")]}
  text.gsub!(/\{\s*\\em([^}]*\{[^}]*\}[^}]*)\}/m){m=$~;sp[m[1].h("em")]}
  text.gsub!(/\{\s*\\em([^}]*)\}/m){m=$~;sp[m[1].h("em")]}
  if state==:math
  text.gsub!(/</){sp["&lt;"]}
  text.gsub!(/>/){sp["&gt;"]}
  texarg=/([^{\\]|\{(?:(?:[^{}]|\{[^{}]*\})*)\}|\\[a-zA-Z]+)/
  text.gsub!(/\^\\prime/om){sp["'"]}
  text.gsub!(/\^\{\\prime\\prime\}/om){sp["''"]}
  text.gsub!(/\\not\s*\\equiv/om){sp["&#8802;"]}
  text.gsub!(/\\(?:text|mbox|hbox)\{([^{}]*)\}/){m=$~;
    sp[TeX2htm(m[1]).h("/i")]}
  text.gsub!(/\^\s*#{texarg}/om){m=$~;sp[TeX2htm(m[1],:math).h("sup")]}
  text.gsub!(/_\s*#{texarg}/om){m=$~;sp[TeX2htm(m[1],:math).h("sub")]}
  text.gsub!(/\\underline\s*#{texarg}/om){m=$~;sp[TeX2htm(m[1],:math).h("u")]}
  text.gsub!(/\\overline\s*#{texarg}/om){m=$~;
    sp['<span style="text-decoration: overline">'+
       TeX2htm(m[1],:math)+"<\/span>"]}
  text.gsub!(/\{\s*\\rm([^}]*)\}/m){m=$~;sp[m[1].h("I")]}
  text.gsub!(/\{\s*\\mathcal([^}]*)\}/m){m=$~;sp[m[1].h("u")]}
  text.gsub!(/\\frac#{texarg}#{texarg}/om){m=$~;
    n=TeX2htm(m[1],'math'); d=TeX2htm(m[2],'math')
    n="(#{n})" if n.length>1
    d="(#{d})" if n.length>1
    "#{n}/#{d}"}
  text.gsub!(/\\not\s*\\in/){sp["&notin;"]}
  text.gsub!(/\\not\s*=/){sp["&ne;"]}
  text.gsub!(/\\([a-zA-Z]+)/){m=$~;
   if MathMacs[m[1]] then sp[MathMacs[m[1]]] else m[0] end}
  text.gsub!(/\\pmod\s*#{texarg}/){m=$~;sp["(</i>mod<i> #{TeX2htm(m[1],:math)})"]}
  text.gsub!(/\\begin\{([^}]*)\}(.*?)\\end\{\1\}/m){m=$~
    case m[1]
    when "array","pmatrix" then 
      if m[1]=="array" then r=m[2].sub(/\{[cl]*\}/,"") else r=m[2] end
      r=r.split("\\cr").map{|l| l.split("&").map{|i| 
         i.h("I").h("td")}.join("").h("tr")}.join.h('table style="display:inline-table;"').h("td").h("/td")
      if m[1]=="pmatrix" then "("+r+")" else r end
    else "{beg-"+m[1]+"}".h("s")+m[2]+"{end-"+m[1]+"}".h("s")
    end
    }
  end
  text.gsub!(/\\(?:text|mbox|hbox)\{([^{}]*)\}/){m=$~;
    sp[TeX2htm(m[1]).h("/i")]}
  text.gsub!(/\\([a-zA-Z]+)/){m=$~;AlwaysMacs[m[1]] || m[0]}
  text.gsub!(/\{([^{}]*)\}/){m=$~; m[1]}
  text.gsub!(/\\([a-zA-Z]+) /){m=$~; $trads[m[1]] || "\\#{m[1]} "}
  text.gsub!(/\\([a-zA-Z]+) /){m=$~; $trads[m[1]] || "\\#{m[1]} "}
  text.gsub!(/\\([a-zA-Z]+)/){m=$~; $trads[m[1]] || "\\#{m[1]}"}
  text.scan(/\\([a-zA-Z]+)/){|m| $err.print "!!!!! #{m.inspect} state=#{state}\n"}
  text
end

# Basically the chapter file is read in one pass, using information previously 
# read from the .toc file to fill in next and previous pointers and the like
#

def navigation (sec)
  chap = sec.chapter
  cfname = name2fn chap.name
  res=""
  next_ref=proc{|n| "<a href =\"#{name2fn n.name}\">Next</a>"}
  prev_ref=proc{|n| "<a href =\"#{name2fn n.name}\">Previous</a> "}
  if sec.secnum == 0
    res<< prev_ref[$chapters[chap.number-1]]unless chap.number==1
    res<< "<a href = \"index.htm\">Up</a> "
    res<< next_ref[$chapters[chap.number+1]]if $chapters[chap.number+1]
  else
    res<< prev_ref[chap.sections[sec.secnum-1]]unless sec.secnum == 1
    res<< "<a href = \"#{cfname}\">Up</a> "
    res<< "<a href = \"index.htm\">Top</a> "
    res<< next_ref[chap.sections[sec.secnum+1]]if chap.sections[sec.secnum+1]
  end
  res<<  "<BR><a href = \"theindex.htm\">Index</a>\n"
  res
end

#
# Main program starts here
#
# Process option and sort out input and output directories   
#

def help
  print <<END
  Usage:
    convert.rb [options] <html-directory>

    -c, --per-chapter
    file-per-chapter  mode -- generates one HTML file chapxxx.htm for each
    chapter sections are level 2 headings and anchors chapxxx.htm#SECTxxx.
    This is intended for local browsing, especially under MS-DOS
 
    -s, --silent  silent running. Conversational messages are suppressed.

    <html-directory> defaults to the current directory
END
end

require 'getoptlong'

$opt_s=$opt_c=false
GetoptLong.new(
 [ "--help-chapter", "-h", GetoptLong::NO_ARGUMENT],
 [ "--per-chapter", "-c", GetoptLong::NO_ARGUMENT],
 [ "--silent",      "-s", GetoptLong::NO_ARGUMENT]).each {|opt,arg|
  if opt=="--silent" then $opt_s=true end
  if opt=="--per-chapter" then $opt_c=true end
  if opt=="--help " then help end
}

$err=File.open("convert.log","w")

help;exit if ARGV.length==0
$dir=File.expand_path(ARGV.shift)
raise "Can't read directory "+$dir unless FileTest.readable?($dir)
print  "Reading input from #{$dir}\n" unless $opt_s

# Scan the .tex and .toc files to get chapter names and numbers, section names
# and numbers and associated filenames Loads up chapters and sections_by_name

# used to standardize section names for use as hash indices.
def canonize (key)
  key.tr("A-Z","a-z").gsub(/\s/,"")
end

$chapters=[]
$sections_by_name=Hash.new

class Chapter
  attr_accessor :name, :number, :sections, :file
  def initialize(name, number, sections)
    @name=name
    @number=number
    @sections=sections
  end
  def inspect
    "Chapter #{number}[#{file}]=#{name}\n"
  end
end

class Section
  attr_accessor :name, :chapnum, :secnum, :chapter
  def initialize(name, chapnum, secnum, chapter)
    @name=name
    @chapnum=chapnum
    @secnum=secnum
    @chapter=chapter
  end
  def inspect
    "#{@name}:#{@chapnum}.#{@secnum}\n"
  end
end

File.open(man=File.join($dir,"manual.toc") ){|toc|
  where=proc{"#{File.basename(man)} line #{toc.lineno}: "}
  for line in toc do
    if line=~/\\numberline\s+\{(\d+)\}(.+)\}\s*\{\d+\}/
      chapno=$1.to_i
      raise "chapter number repeated" if $chapters[chapno] 
      $chapters[chapno] = Chapter.new($2, chapno, [])
      chap_as_sec = Section.new($2,chapno,0,$chapters[chapno])
      $chapters[chapno].sections[0] = chap_as_sec
      $sections_by_name[canonize(chap_as_sec.name)] = chap_as_sec
    elsif line=~/\\numberline\s+\{(\d+)\.(\d+)\}(.+)\}\s*\{\d+\}/
      chapno=$1.to_i
      secno=$2.to_i
      raise where[]+" section #{$2} in unknown chapter #{$1}" unless 
        $chapters[chapno] 
      raise "section number #{$2} repeated" if $chapters[chapno].sections[secno]
      sec = Section.new($3,chapno,secno,$chapters[chapno])
      $sections_by_name[canonize(sec.name)] = sec
      $chapters[chapno].sections[secno] = sec
    else 
      $err.print "Bad line: #{where[]} <#{line}>"
    end
  end
}
print  "Processed TOC file\n" unless $opt_s

File.open(File.join($dir,"manual.tex")){|tex|
  chapno = 0
  for line in tex do
    if line=~ /^[^%]*\\Include\{(.+)\}/  
      name=File.join($dir,$1+".tex")
      if not FileTest.exist?(name) or not FileTest.readable?(name) 
	$err.print "cannot read "+name+"\n"
      end
      chapno+=1
      begin
        $chapters[chapno].file = $1
      rescue
        raise "for #{name} chapno=#{chapno} missing\n"
      end
    end
  end
}

$odir=File.expand_path(ARGV.shift|| ".") 
raise "Can't write to directory "+$odir unless FileTest.writable?($odir) 
print  "Creating output in "+$odir+"\n" unless $opt_s

# OK go to work
$index=Hash.new{Array.new}

for chap in $chapters do
  next unless chap
  print  chap.inspect unless $opt_s
  $err.print  chap.inspect
  $indexcount = -1

  # loop, controlled by the list of sections that we expect.
  # Will fail, possibly messily if this does not match reality
   
  infile=File.open(File.join($dir,chap.file+".tex")) {|f| f.read}

  for sec in chap.sections do
    $err.print " Section #{sec.chapnum}.#{sec.secnum} #{sec.name}\n"
    # sort out what we are processing (chapter or section)
    if sec.secnum == 0
      num = chap.number
      name = chap.name
      re= "^\\\\Chapter\\{#{Regexp.escape name}\\}"
    else
      num = sec.chapnum.to_s + "." + sec.secnum.to_s
      name = sec.name
      re= "^\\\\Section\\{#{Regexp.escape name}\\}"
    end

    if not $opt_c or sec.secnum==0
      fname = name2fn sec.name
      if fname =~ /\#/ then raise "Filename contains #" end
      out=File.open(File.join($odir,fname),"w")
    # produce the header of the Web page
      out.print  "<html><head><title>GAP3 Manual: #{num} #{name}</title></head>\n"
      out.print  "<body bgcolor=\"ffffff\">\n<h1>#{num} #{name}</h1>\n<P>"
    else
      out.print "<A NAME=\"SECT#{pad(sec.secnum)}\">"+
        "#{num} #{TeX2htm(name)}".h("h2")+"</a>\n<P>"
    end

    # Look for the \Chapter or \Section line

    if infile !~/#{re}/
      raise "Missing chapter or section line matching "+re
    end

    
    infile=$'
    
    if  infile !~/^\%{16,}/ 
      out.print TeX2htm(infile,:normal,fname,sec)
    else
      infile=$'
      out.print TeX2htm($`,:normal,fname,sec)
    end

    # Here we have processed the whole section and start to attach footers
    # to it. If it is really a chapter then it gets a list of its sections

    if sec.secnum == 0
      out.print  "<P>\n<H3> Subsections</H3>\n<oL>\n"
      for subsec in chap.sections do
	next if subsec.secnum==0
	out.print  "<LI> <A HREF=\"#{name2fn subsec.name}\">"+
          TeX2htm(subsec.name)+"</a>\n"
      end
      out.print  "</ol>\n"
    end
    unless $opt_c
      out.print navigation(sec)
      out.print $footer
      out.close
    end
  end
  if $opt_c
    out.print navigation(chap.sections[0])
    out.print $footer
    out.close
  end
end

print "creating index.htm\n" unless $opt_s
File.open(File.join($odir,"index.htm"),"w"){|chaps|
chaps.print  <<END
<html><head><title>The GAP3 Manual -- Chapters</title></head>
<body bgcolor=\"ffffff\"><h1>The GAP3 Manual -- Chapters</h1><ol>
END
$chapters.compact.each { |chp|
  chaps.print  "<LI> <A HREF=\"#{name2fn chp.name}\">#{chp.name}</a>\n"
}
chaps.print  <<END
</ol>\n<P>\n<a href=\"biblio.htm\">References</a><P>
<a href=\"theindex.htm\">Index</a><P>
<a href=\"index.htm\">Up</a><P>
END
chaps.print  $footer
}

print "creating biblio.htm\n" unless $opt_s
File.open(File.join($odir,"biblio.htm"),"w"){|bib|
File.read(File.join($dir,"manual.bbl")).split("\n\n").each{|bl|
  v=proc{|s|
    TeX2htm(s.gsub(/\n/," ").gsub(/\{([A-Z])\}/,'\1')).
      gsub(/\[([a-zA-Z.]*)\]/,'\1').gsub(/\{([^}]*)\}/,'\1')
  }
  blk=/\\newblock(.*)/m
  case bl
  when /\\bibitem\[([^\]]*)\]\{([^}]*)\}(.*)#{blk}#{blk}#{blk}/m 
  then
    m=$~
    s="<A name=#$2>#{$2.h("b")}</A>"+v[m[3]]
    s+=(4..6).map{|i|v[m[i]]}.join("<br>\n   <sp>").h("blockquote")
    bib.print s.h("li")
  when /\\bibitem\[([^\]]*)\]\{([^}]*)\}(.*)#{blk}#{blk}/m then
    m=$~
    s="<A name=#$2>#{$2.h("b")}</A>"+v[m[3]]
    s+=(4..5).map{|i|v[m[i]]}.join("<br>\n   <sp>").h("blockquote")
    bib.print s.h("li")
  when /\\bibitem\[([^\]]*)\]\{([^}]*)\}(.*)#{blk}/m then
    m=$~
    s="<A name=#$2>#{$2.h("b")}</A>"+v[m[3]]
    s+=(4..4).map{|i|v[m[i]]}.join("<br>\n   <sp>").h("blockquote")
    bib.print s.h("li")
  end
  bib.print "\n"
}}

print "creating theindex.htm\n" unless $opt_s
File.open(File.join($odir,"theindex.htm"),"w"){|index|
index.print <<END
<html><head><title>The GAP3 Manual -- Index</title></head>
<body bgcolor=\"ffffff\"><h1>The GAP3 Manual -- Index</h1>
<P>
END
("A".."Z").each { |letter|
  index.print  "<a href=\"theindex.htm#L#{letter}\">#{letter}</A> "
}
index.print  "\n<ul>"
nextletter = ?A.ord
print "#{$index.keys.length} index entries\n"
for ent in $index.keys.sort { |a,b| a.downcase<=>b.downcase or a<=>b }do
  next if ent.length==0
#  print ent,"\n"
  letter = ent.upcase[0].ord || 32
  if nextletter <= ?Z.ord
    until letter <= nextletter 
      index.print  "<A name = \"L#{"%c"%nextletter}\"></a>"
      nextletter+=1
    end
  end
  index.print  "<LI>"+ent.sub(/!/," ")+" "
  for ref in $index[ent] do
    index.print  "<A HREF=\"#{ref[0]}\">#{ref[1]}</A> "
  end
  index.print  "\n"
end
index.print  "</ul>\n<P>\n"
index.print  "<a href=\"index.htm\">Up</a><P>\n"
index.print  $footer
}

print  "done\n" unless $opt_s
