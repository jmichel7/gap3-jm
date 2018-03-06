#!/bin/csh
unset noclobber
foreach example (examples/*.[0-9]*)
    echo -n testing example `basename $example "\.[0-9]*"`:
    set class=$example:e
    bin/nq $example $class | grep -v "runtime\|size\|Machine\|Input" \
		> $example:r.tst
    diff $example:r.out $example:r.tst > /dev/null
    if( $status != 0 ) then
	echo " error."
	echo Please mail the file $example:r.tst 
	echo "    " to Werner.Nickel@math.rwth-aachen.de
    else
        echo " ok."
	rm $example:r.tst
    endif
end
