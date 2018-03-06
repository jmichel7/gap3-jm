#!/bin/sh

#  This script converts the files from the cunningham project 
#  into GAP readable format

# # # # # # # # # extract primes from "n+" and "-" # # # # # # # # #

# remove exponent part
sed -e 's/^[0-9]*[LM]*//' [0-9]*[+-] > /tmp/primes.$$

# throw away lines containing a composite number
fgrep -v "C" /tmp/primes.$$ > /tmp/primesx.$$
mv /tmp/primesx.$$ /tmp/primes.$$

# throw away unspecified primes
sed -e 's/\.[ 	]*P[0-9]*//' /tmp/primes.$$ > /tmp/primesx.$$
mv /tmp/primesx.$$ /tmp/primes.$$

# and replace "." by ","
tr "." "," < /tmp/primes.$$ | tr -d " 	" > /tmp/primesx.$$
mv /tmp/primesx.$$ /tmp/primes.$$

# remove lines conataing only a "Pxyz"
egrep -v "P[0-9]" /tmp/primes.$$ > /tmp/primesx.$$
mv /tmp/primesx.$$ /tmp/primes1.$$
rm /tmp/primes.$$

# # # # # # # # # extract the primes from "primes" # # # # # # # # #

# throw away "There are ... main tables" lines
egrep -v "^    There (are|is) [0-9]+ additional [0-9]+-digit" primes > /tmp/primes.$$

# throw away "Appendix A" line
egrep -v "       Appendix A: Primes" /tmp/primes.$$ > /tmp/primesx.$$
mv /tmp/primesx.$$ /tmp/primes.$$

# now combine split lines
awk '
  $0 ~ /^                  / {
    printf( "%s", $0 );
    next;
  }

  {
    printf( "\n%s", substr($0,length(" 21  3,525L  P  ")) );
  }
  END {
    print
  }

' /tmp/primes.$$ | tr -d " " > /tmp/primesx.$$
mv /tmp/primesx.$$ /tmp/primes.$$
cat /tmp/primes1.$$ >> /tmp/primes.$$
rm /tmp/primes1.$$

# create GAP input
awk '
  BEGIN {
    printf( "%s", "Cunningham := [ 2" );
  }
  {
    if ( 0 < length($0) )
      printf( ",\n%s", $0 );
  }
  END {
    printf( "\n];\n" );
  }

' /tmp/primes.$$ > /tmp/gap.$$
rm /tmp/primes.$$

# now find new primes for GAP
(
  echo "Read(\"/tmp/gap.$$\");"
  echo "Cunningham := Difference( Cunningham, Primes );;"
  echo "Cunningham := Difference( Cunningham, Primes2 );;"
  echo "PrintTo(\"/tmp/gap.$$\",\"Cunningham := \",Cunningham,\";\n\");"
  echo "Print( \"CHECK WITH: Collected(List(Cunningham,IsPrimeInt)); \");"
) | gap -b
mv /tmp/gap.$$ primes.gap


