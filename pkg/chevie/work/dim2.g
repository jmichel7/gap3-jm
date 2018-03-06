# computes reflection group <x,y,z| x^ox=y^oy=z^oz=1, xyz=yzx=zxy>
ddim2:=function(ox,oy,oz)local b,x,y,z;
  b:=E(2*ox)*E(2*oy)*E(2*oz); x:=E(ox);y:=E(oy);z:=E(oz);
  return PermRootGroup([[x+y+b+b/z,1-x],[y-1,1],[z+b/x,b/x/y]],
                       [[0,1],[-1,0],[-1,y*(1+x/b)]]);
end;
