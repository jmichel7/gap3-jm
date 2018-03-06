e:=(  2,  5)(  6,119)(  9, 11)( 12,117)( 13,116)( 15, 17)( 18,114)( 19,113)
( 20, 21)( 23,110)( 24,111)( 25,109)( 27,108)( 29,107)( 30,105)( 31,106)
( 33,104)( 34,103)( 35,102)( 36,101)( 38, 99)( 39, 98)( 40, 97)( 41, 96)
( 42, 93)( 43, 94)( 44, 92)( 45, 90)( 46, 88)( 47, 86)( 48, 87)( 49,126)
( 50, 82)( 51, 81)( 52,125)( 53, 76)( 54, 75)( 55,124)( 56, 69)( 57,123)
( 58,121)( 59,122)( 60,120)( 61,118)( 62,115)( 63,112)( 65, 68)( 72, 74)
( 78, 80)( 83, 84);

HF:=CoxeterSubCoset(CoxeterCoset(CoxeterGroup("E",7)),[1,2,3,4,5,7],e);

InitClassesCharTable:=function ( tbl ) local  x, order, initclasses;
   if not IsInt( tbl.centralizers[1] )  then return; fi;
   order := Maximum(tbl.centralizers);
   initclasses := List( tbl.centralizers, x->order/x);
   if not ForAll( initclasses, IsInt )  then
       Print( "#E InitClassesCharTable: not all centralizer orders divide", 
        " the group order\n" );
   fi;
   tbl.classes := initclasses;
   return initclasses;
end;

CharTable(HF);
