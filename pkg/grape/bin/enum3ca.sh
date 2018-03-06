#!/bin/sh
if test $1 
then gens=$1
else gens=abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ
fi
echo "0 1 1 0" >GRAPE_tcin
cat >GRAPE_tcfein
$GRAPE_BIN/tcfrontend3 <GRAPE_tcfein
cat GRAPE_tcfeout >>GRAPE_tcin
echo 'GROUP PRESENTATION'
cat GRAPE_tcfein 
echo "\input /galois/home/leonard/ca/lstyle" > GRAPE_coladj.tex
echo "\begin{document}" >> GRAPE_coladj.tex
echo "\begin{verbatim}" >> GRAPE_coladj.tex
cat GRAPE_tcfein >> GRAPE_coladj.tex
echo "\end{verbatim}" >> GRAPE_coladj.tex
$GRAPE_BIN/tcmainca3 | $GRAPE_BIN/coladjg3t $gens
cat GRAPE_tcout
echo "\begin{verbatim}" >> GRAPE_coladj.tex
echo " " >> GRAPE_coladj.tex
cat GRAPE_tcout >> GRAPE_coladj.tex
echo "\end{verbatim}" >> GRAPE_coladj.tex
echo "\end{document}" >> GRAPE_coladj.tex
if test $2
then mv GRAPE_coladj.g $2.g; mv GRAPE_coladj.tex $2.tex
fi
rm GRAPE_tcfein GRAPE_tcfeout GRAPE_tcin GRAPE_tcout
