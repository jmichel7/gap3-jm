l := AllIrreducibleSolvableGroups();
failed := [];
nono := [ 94, 95, 96, 98, 99, 154, 156, 159, 186, 191, 318, 322,
          326, 328, 369, 370 ];

for i  in Difference( [ 1 .. Length(l) ], nono )  do
    grp := l[i];
    p := PermGroupRepresentation( grp, 100 );
    if p = false  then
        Add( failed, i );
    fi;
od;
for i  in nono  do
    grp := l[i];
    p := PermGroupRepresentation( grp, 1000 );
    if p = false  then
        Add( failed, i );
    fi;
od;

Print( "#I  'PermGroupRepresentation' for irreducible solvables: failed = ",
       Length(failed), "\n" );

#############################################################################

files := [ 
    "2j2d36.gap", "2j2d42.gap", "2m11d5.gap", "2m22d34.gap", "3m22_2.gap",
    "3m22_2b.gap", "3m22d12.gap", "a5.gap", "a5d4.gap", "a5xa5d25.gap",
    "a6c6d12.gap", "a6d16.gap", "a7c6d24.gap", "a7d20.gap", "a7d6.gap",
    "a8d20.gap", "a8d28.gap", "a8d64.gap", "a9d28.gap", "co2.gap",
    "co3d22.gap", "esd4.gap", "f29d18.gap", "f29d81.gap", "j1.gap",
    "j1d14.gap", "j1d27.gap", "l213.gap", "l217.gap", "l231.gap", "l281.gap",
    "l33d26.gap", "l34d63a.gap", "l34d63b.gap", "l35d124.gap",
    "l35d124b.gap", "l35d124c.gap", "m11d11.gap", "m11d24.gap", "m11d44.gap",
    "m11d44f2.gap", "m11d5.gap", "m11d55.gap", "m11wrm11.gap", "m12d120.gap",
    "m12d55.gap", "m22d154.gap", "m22d21.gap", "m22d30.gap", "m22d54.gap",
    "m23d11.gap", "m24d11.gap", "norm243.gap", "norm343.gap", "norm64.gap",
    "onand18.gap", "onand27.gap", "onand45.gap", "s22d21.gap", "sl25d6.gap",
    "sz8d65.gap", "u42d14.gap", "u42d30.gap", "u42d30ts.gap", "u42d58.gap",
    "u42d64.gap", "u42d81.gap"
];


for max  in [ 10^3, 10^4 ]  do
    failed := [];
    for file  in files  do
        ReadDataPkg( "matrix", "data", file );
        start := Runtime();
        p := PermGroupRepresentation( G, max );
        if p = false  then
            Print( "#I  \"", file, "\" failed with <max> = ", max,
                   ", time = ", Runtime()-start, "\n" );
            Add( failed, file );
        else
            Print( "#I  \"", file, "\" succeeded, degree = ",
                   Length(G.permDomain), ", time = ", Runtime()-start,
                   "\n" );
        fi;
    od;
    files := failed;
od;
Print( "#I  failed on: ", failed, "\n" );
