#!/bin/sh
execDir=`dirname $0`

base1=`basename $1`
base2=`basename $2`
diffFile1=`mktemp $base1.XXXXXX`
diffFile2=`mktemp $base2.XXXXXX`

head -1 $1 > $diffFile1
head -1 $2 > $diffFile2

diff -q --strip-trailing-cr $diffFile1 $diffFile2
status=$?
rm -f $diffFile1 $diffFile2

exit $status