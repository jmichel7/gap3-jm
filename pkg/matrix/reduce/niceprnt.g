#############################################################################
##
#A  Matrix package                                               Anthony Pye 
##                                                      
#Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
##
##
PrintLayer := function(record,layerNmr)
    
    if layerNmr < 1 then
        Error(" Layers start from 1 ");
    fi;
    if layerNmr = 1 then
        if KernelFlag(record) = "unknown" then
            Print("#I  Layer Number = ",LayerNumberFlag(record),
                  "\n#I  Type = ", TypeFlag(record), "\n#I  Dimension = ",
                  DimensionFlag(record), "\n#I  Size = ",Size(record),"\n");
        else
            Print("#I  Layer Number = ", LayerNumberFlag(record),
                  "\n#I  Type = ", TypeFlag(record),"\n#I  Dimension = ",
                  DimensionQuotientFlag(record),"\n#I  Size = ",
                  SizeQuotientFlag(record),"\n");
        fi;
    elif KernelFlag(record) <> "unknown" then
        PrintLayer(KernelFlag(record),layerNmr-1);
    else
        Error(" Computation doesn't reach this layer ");
    fi;
    
end;

PrintLayers := function(record)
    
    if KernelFlag(record) = "unknown" then
        Print("#I  Layer Number = ",LayerNumberFlag(record),
              "\n#I  Type = ", TypeFlag(record), "\n#I  Dimension = ",
              DimensionFlag(record),"\n#I  Size = ", SizeFlag(record),"\n");
        return;
    else
        Print("#I  Layer Number = ",LayerNumberFlag(record),"\n#I  Type = ",
              TypeFlag(record), "\n#I  Dimension = ",
              DimensionQuotientFlag(record),"\n#I  Size = ", 
              SizeQuotientFlag(record),"\n");
    fi;
    PrintLayers(KernelFlag(record));
    
end;

PrintNumberOfLayers := function(record)
    
    if KernelFlag(record) = "unknown" then
        Print("#I  Number of layers is ",LayerNumberFlag(record),"\n");
    else
        PrintNumberOfLayers(KernelFlag(record));
    fi;
    
end;

PrintBasicInfo := function(record)
    
    Print("#I  Matrix group over field ",FieldFlag(record)," of dimension ",
          DimensionFlag(record)," has size ",SizeFlag(record),"\n");
    
end;

DisplayMatRecord := function(arg)
    local tmp, i;
    
    if Length(arg) = 1 then
        
        if   PrintLevelFlag(arg[1]) = 1 then
            PrintBasicInfo(arg[1]);
            PrintNumberOfLayers(arg[1]);
        elif PrintLevelFlag(arg[1]) = 2 then
            PrintBasicInfo(arg[1]);
            PrintLayers(arg[1]);
        elif PrintLevelFlag(arg[1]) = 3 then
            Print(arg[1],"\n");
        fi;
        
    else
        
        if   PrintLevelFlag(arg[1]) = 1 then
            PrintLayer(arg[1],arg[2]);
        elif PrintLevelFlag(arg[1]) = 2 then
            PrintLayer(arg[1],arg[2]);
        elif PrintLevelFlag(arg[1]) = 3 then
            
            tmp := arg[1];
            i := 1;
            while i <> arg[2] do
                if KernelFlag(tmp) = "unknown" then
                    Error(" Computation didn't reach such a layer ");
                fi;
                tmp := KernelFlag(tmp);
                i := i + 1;
            od;
            if KernelFlag(tmp) <> "unknown" then
                tmp := Copy(tmp);
                UndoKernelFlag(tmp);
            fi;
            Print(tmp,"\n");
        fi;
    fi;
    
end;

