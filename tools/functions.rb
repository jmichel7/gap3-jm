#!/usr/bin/ruby
ARGV.sort.each{|f|
  s=File.open(f){|h| h.read}
  list=[]
  s.scan(/\s+([\w]+)\s*:=\s*function/) { list<<$1}
  s.scan(/\s+([\w]+)\s*:=\s*\w+\s*->/) { list<<$1}
  s.scan(/\s+([\w]+)\s*:=\s*Dispatcher/) { list<<$1}
  s.scan(/\s+([\w]+)\s*:=\s*AttributeDispatcher/) { list<<$1}
  s.scan(/\s+([\w]+)\s*:=\s*OperationsRecord/) { list<<$1}
  print "#{f} defines:\n"
  list.sort.each{|f| print "  ",f,",\n"}
}
