#!/bin/csh
#
# This is a modified copy of the file testNq from Werner Nickel

unset noclobber
foreach example (examples/*.fp)
    echo -n Testing $example 
    bin/Sq -p1 < $example | grep -v '#T' | fgrep -v runtime > $example:r.tst
    echo -n '  . . . . . . . . '
    diff $example:r.out $example:r.tst > /dev/null
    if( $status != 0 ) then
        echo Error 
        echo Please mail the file $example:r.tst
        echo "    " to alice@maths.uwa.edu.au
    else
        rm $example:r.tst
	echo  succeeded
    endif
end
