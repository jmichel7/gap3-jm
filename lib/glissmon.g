# glissmon.g
#
#  Jean Michel apr 2016. 
# Hack to provide compatibility between monoid and gliss packages.
#
Transformation:=function(arg)
  if Length(arg)=1 then return ApplyFunc(Transformation_mon,arg);
  elif Length(arg)=2 then return ApplyFunc(Transformation_gliss,arg);
  else Error("should have 1 or 2 arguments");
  fi;
end;

IdentityTransformation:=function(x)
  if IsInt(x) then return IdentityTransformation_mon(x);
  else return IdentityTransformation_gliss(x);
  fi;
end;
