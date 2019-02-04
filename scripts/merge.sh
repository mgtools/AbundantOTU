#!/bin/bash
rm -f $2
a=0
while read line
do
    file=$line
    cat ${file} >> $2 
    n=`grep ">" ${file} -c`
    b=${a}
    a=$((a+n))
    echo "${file} ${b} ${a}"
done < $1
