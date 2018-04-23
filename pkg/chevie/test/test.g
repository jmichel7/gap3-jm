#############################################################################
##
#A  test.g           CHEVIE library                            Jean Michel
#Y  Copyright (C) august 2017 -  University  Paris Diderot.
#
# Testing framework for CHEVIE. The framework is loaded by loading this file.
# The main functions defined here are:
#
# CHEVIE.Test(W) applies various tests to spets W
#
# CHEVIE.RegressionTest() applies CHEVIE.Test to all CHEVIE.TestObjs
#
CHEVIE.TestObjs:=[ # a selection of spets to be tested
  [CoxeterGroup],
  [CoxeterGroup,"A",1],
  [RootDatum,"psu",3],
  [CoxeterGroup,"A",2],
  [CoxeterGroup,"I",2,5],
  [CoxeterGroup,"B",2],
  [CoxeterGroup,"G",2],
  [CoxeterGroup,"A",3],
  [RootDatum,"2B2"],
  [RootDatum,"2G2"],
  [RootDatum,"2I",5],
  [RootDatum,"2I",8],
  [ComplexReflectionGroup,5],
  [RootDatum,"psu",4],
  [ComplexReflectionGroup,4],
  [ComplexReflectionGroup,9],
  [ComplexReflectionGroup,12],
  [ComplexReflectionGroup,13],
  [RootDatum,"3D4"],
  [RootDatum,"psu",5],
  [ComplexReflectionGroup,3,1,2],
  [ComplexReflectionGroup,3,3,3],
  [CoxeterGroup,"B",3],
  [CoxeterGroup,"C",3],
  [CoxeterGroup,"A",4],
  [ComplexReflectionGroup,7],
  [ComplexReflectionGroup,10],
  [ComplexReflectionGroup,22],
  [ComplexReflectionGroup,4,4,3],
  [CoxeterGroup,"H",3],
  [RootDatum,"pso-",8],
  [ComplexReflectionGroup,6],
  [RootDatum,"psu",6],
  [RootDatum,"2F4"],
  [ComplexReflectionGroup,8],
  [ComplexReflectionGroup,24],
  [CoxeterGroup,"D",4],
  [CoxeterGroup,"A",5],
  [ComplexReflectionGroup,3,3,4],
  [CoxeterGroup,"B",4],
  [CoxeterGroup,"C",4],
  [CoxeterGroup,"F",4],
  [ComplexReflectionGroup,15],
  [ComplexReflectionGroup,16],
  [ComplexReflectionGroup,20],
  [ComplexReflectionGroup,25],
  [RootDatum,"psu",7],
  [CoxeterGroup,"A",6],
  [CoxeterGroup,"D",5],
  [RootDatum,"pso-",10],
  [ComplexReflectionGroup,14],
  [CoxeterGroup,"E",6],
  [ComplexReflectionGroup,11],
  [ComplexReflectionGroup,29],
  [CoxeterGroup,"B",5],
  [CoxeterGroup,"C",5],
  [ComplexReflectionGroup,26],
  [CoxeterGroup,"A",7],
  [RootDatum,"2E6"],
  [RootDatum,"psu",8],
  [ComplexReflectionGroup,17],
  [CoxeterGroup,"H",4],
  [ComplexReflectionGroup,33],
  [RootDatum,"pso-",12],
  [CoxeterGroup,"D",6],
  [CoxeterGroup,"B",6],
  [CoxeterGroup,"C",6],
  [CoxeterGroup,"E",7],
  [CoxeterGroup,"D",7],
  [ComplexReflectionGroup,27],
  [ComplexReflectionGroup,18],
  [ComplexReflectionGroup,21],
  [ComplexReflectionGroup,32],
  [CoxeterGroup,"B",7],
  [CoxeterGroup,"C",7],
  [RootDatum,"pso-",14],
  [CoxeterGroup,"E",8],
  [ComplexReflectionGroup,34]];

CHEVIE.TestList:=[];

CHEVIE.AddTest:=function(arg)local res;
  res:=rec(name:=arg[1],test:=arg[2]);
  if Length(arg)>=3 then res.cond:=arg[3];else res.cond:=x->true;fi;
  if Length(arg)=4 then res.explanation:=arg[4];fi;
  Add(CHEVIE.TestList,res);
end;

CHEVIE.Title:=[];

# The next function is used for describing what is being tested
# Title[1] is the name of current test, Title[2] details the test if needed
CHEVIE.Testing:=function(arg)CHEVIE.Title[2]:=ApplyFunc(SPrint,arg);end;

# The next two functions are destined to be bound to ChevieErr, called
# if a problem is encountered when testing.
CHEVIE.Error:=function(arg)
  ApplyFunc(Error,Concatenation([Concatenation(CHEVIE.Title,[":"])],arg));
end;

CHEVIE.Warn:=function(arg)local s,t;
  s:=ApplyFunc(SPrint,arg);
  CHEVIE.errorDiagnosed:=true;
  t:=SPrint("\n****** WARNING!:",Concatenation(CHEVIE.Title),s);
  Cut(t,rec(after:=", =",before:="+"));
  if IsBound(CHEVIE.log) then
    t:=SPrint("[",Concatenation(CHEVIE.Title),"]",s);
    Cut(t,rec(after:=", =",before:="+",file:=CHEVIE.log));
  fi;
end;

ChevieErr:=CHEVIE.Warn;

# CHEVIE.Test([test or [tests],] W [,...])
# test W by all CHEVIE.TestList (default)
# or by given test or tests
CHEVIE.Test:=function(arg)local a,tests,i,W,ttest;
  if IsString(arg[1]) then tests:=Filtered(CHEVIE.TestList,x->x.name=arg[1]);
    W:=arg[2];
    arg:=arg{[3..Length(arg)]};
  elif IsList(arg[1]) then tests:=arg[1];W:=arg[2];
    arg:=arg{[3..Length(arg)]};
  else tests:=CHEVIE.TestList;W:=arg[1];
    arg:=arg{[2..Length(arg)]};
  fi;
  for i in [1..Length(tests)] do
    if Length(tests)>1 then InfoChevie(i,"/",Length(tests),":\c");fi;
    a:=tests[i];
    if a.cond(W) then
      CHEVIE.errorDiagnosed:=false;
      CHEVIE.Title[1]:=SPrint(a.name,"(",ReflectionName(W),")");
      if Length(arg)=0 then
        if IsBound(a.explanation) then InfoChevie2(a.explanation,"\n");fi;
        InfoChevie("Check[",CHEVIE.Title[1],"]\c");
        ttest:=Runtime();
      fi;
      ApplyFunc(a.test,Concatenation([W],arg));
      if Length(arg)=0 then
        if not CHEVIE.errorDiagnosed then InfoChevie("OK!");fi;
        InfoChevie(StringTime(Runtime()-ttest));
        InfoChevie("\n");
      fi;
      CHEVIE.Testing();
    else
      InfoChevie(a.name," not applicable\n");
    fi;
  od;
end;

# CHEVIE.RegressionTest(test or [tests])
# Apply all possible tests (or argument tests) to all TestObjs
CHEVIE.RegressionTest:=function(arg)local W,msg,i,tstart,objs,tests;
  objs:=CHEVIE.TestObjs;
  CHEVIE.log:="all.log";ChevieErr:=CHEVIE.Warn;
  if Length(arg)=0 then tests:=CHEVIE.TestList;else tests:=arg[1];fi;
  for i in [1..Length(objs)] do Print(i,"/",Length(objs),":");
    W:=objs[i];
    if IsList(W) then W:=ApplyFunc(W[1],W{[2..Length(W)]});fi;
    tstart:=Runtime();
    msg:=SPrint("======== diagnostics for ",ReflectionName(W)," =========\n");
    if IsBound(CHEVIE.log) then AppendTo(CHEVIE.log,msg);fi;
    Print(msg);
    CHEVIE.Test(tests,W);
    msg:=SPrint(Elapsed(tstart),"\n");
    if IsBound(CHEVIE.log) then AppendTo(CHEVIE.log,msg);fi;
    Print(msg);
  od;
  Unbind(CHEVIE.log);
end;

IsSpetsial:=W->UnipotentCharacters(W)<>false;

# now load various tests
ReadChv("test/check");
ReadChv("test/checkcharpara");
CHEVIE.AddTest("CharParams",CheckCharParams,W->not IsSpets(W));
ReadChv("test/checkbasic");
ReadChv("test/checkchar");
ReadChv("test/checkunideg");
ReadChv("test/checkunipclasses");
