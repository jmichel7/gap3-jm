CC          = gcc 
#COPTS       = -O2 -m32 -Wno-unused-result -g
COPTS       = -m32 -Wno-unused-result -g
LOPTS       = -m32 -g
CFLAGS      = $(COPTS)

PCQA_SRCS = ab_mem.c abelian.c change.c col_debug.c col_io.c \
  col_mem.c collect.c common.c complete.c const_debug.c const_io.c \
  const_mem.c construct.c ext_debug.c ext_io.c ext_mem.c extend.c flag.c \
  gbasis.c gbasis_basic.c gbasis_debug.c gbasis_io.c gbasis_mem.c \
  group.c group_io.c group_mem.c head.c head_mem.c hnf.c hnf_debug.c \
  hnf_mem.c homom.c init.c module.c parser.c polypre.c pcqa.c present.c \
  print.c quot_io.c quot_mem.c quotient.c random.c read.c reader.c rel_io.c \
  rel_mem.c relator.c ring.c ring_adv.c ring_basic.c ring_mem.c rule.c \
  rule_io.c rule_mem.c rule_plan.c save.c smith.c stack.c term.c util.c

PCQA_OBJS = pcqa.o ab_mem.o abelian.o change.o col_debug.o col_io.o \
  col_mem.o collect.o common.o complete.o const_debug.o const_io.o \
  const_mem.o construct.o ext_debug.o ext_io.o ext_mem.o extend.o flag.o \
  gbasis.o gbasis_basic.o gbasis_debug.o gbasis_io.o gbasis_mem.o \
  group.o group_io.o group_mem.o head.o head_mem.o hnf.o hnf_debug.o \
  hnf_mem.o homom.o init.o module.o parser.o polypre.o present.o \
  print.o quot_io.o quot_mem.o quotient.o random.o read.o reader.o rel_io.o \
  rel_mem.o relator.o ring.o ring_adv.o ring_basic.o ring_mem.o rule.o \
  rule_io.o rule_mem.o rule_plan.o save.o smith.o stack.o term.o util.o

pcqa : $(PCQA_OBJS)
	$(CC) -o ../bin/pcqa $(PCQA_OBJS) $(LOPTS) -lgmp

cleano:
	rm -f $(PCQA_OBJS)

pcqa.h : $(GMPLIB) ab.h col.h construct.h ext.h data.h def.h gbasis.h \
  group.h head.h io.h main.h module.h parser.h quot.h reader.h relator.h \
  ring.h rule.h smith.h term.h

ab_mem.o : ab_mem.c pcqa.h
abelian.o : abelian.c pcqa.h
change.o : change.c pcqa.h
col_debug.o : col_debug.c pcqa.h
col_io.o : col_io.c pcqa.h
col_mem.o : col_mem.c pcqa.h
collect.o : collect.c pcqa.h
common.o : common.c pcqa.h
complete.o : complete.c pcqa.h
const_debug.o : const_debug.c pcqa.h
const_io.o : const_io.c pcqa.h
const_mem.o : const_mem.c pcqa.h
construct.o : construct.c pcqa.h
ext_debug.o : ext_debug.c pcqa.h
ext_io.o : ext_io.c pcqa.h
ext_mem.o : ext_mem.c pcqa.h
extend.o : extend.c pcqa.h
flag.o : flag.c pcqa.h
gbasis.o : gbasis.c pcqa.h
gbasis_basic.o : gbasis_basic.c pcqa.h
gbasis_debug.o : gbasis_debug.c pcqa.h
gbasis_io.o : gbasis_io.c pcqa.h
gbasis_mem.o : gbasis_mem.c pcqa.h
group.o : group.c pcqa.h
group_io.o : group_io.c pcqa.h
group_mem.o : group_mem.c pcqa.h
head.o : head.c pcqa.h
head_mem.o : head_mem.c pcqa.h
hnf.o : hnf.c pcqa.h
hnf_debug.o : hnf_debug.c pcqa.h
hnf_mem.o : hnf_mem.c pcqa.h
homom.o : homom.c pcqa.h
init.o : init.c pcqa.h
module.o : module.c pcqa.h
parser.o : parser.c pcqa.h
polypre.o : polypre.c pcqa.h
pcqa.o : pcqa.c pcqa.h
present.o : present.c pcqa.h
print.o : print.c pcqa.h
quot_io.o : quot_io.c pcqa.h
quot_mem.o : quot_mem.c pcqa.h
quotient.o : quotient.c pcqa.h
random.o : random.c pcqa.h
read.o : read.c pcqa.h
reader.o : reader.c pcqa.h
rel_io.o : rel_io.c pcqa.h
rel_mem.o : rel_mem.c pcqa.h
relator.o : relator.c pcqa.h
ring.o : ring.c pcqa.h
ring_adv.o : ring_adv.c pcqa.h
ring_basic.o : ring_basic.c pcqa.h
ring_mem.o : ring_mem.c pcqa.h
rule.o : rule.c pcqa.h
rule_io.o : rule_io.c pcqa.h
rule_mem.o : rule_mem.c pcqa.h
rule_plan.o : rule_plan.c pcqa.h
save.o : save.c pcqa.h
smith.o : smith.c pcqa.h
stack.o : stack.c pcqa.h
term.o : term.c pcqa.h
util.o : util.c pcqa.h
