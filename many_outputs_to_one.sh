#!/bin/sh 
set -u
if [ $# -gt 0 ]; then 
 R=$1
 N=`grep 'r_ <-' "$R" | uniq | wc -l` 

 Name="${R%.R}"
 sed -n "/${Name} <-/,/vf_ <-/p" "$R" | sed -e "s/${Name} <-/${Name}_out <-/" -e "s/vf_ <- vector.*\$/r_ <-vector(len=${N})/"
 grep 'r_ <- ' "$R" | uniq | awk -v r_N="$N" '{print "    r_[" NR "]" " <- " $3}'
 echo "    return(r_)"
 echo "}"
fi
