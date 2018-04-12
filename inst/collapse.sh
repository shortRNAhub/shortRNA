#!/bin/bash
if [ $# -lt 1 ]
then
    echo "Usage: `basename $0` input [output] [minLength] [maxLength]"
    echo "By default, output is saved in basename.seqcounts, minLength=15 and maxLength=49."
    exit
fi
tmp=`basename $1 .gz`
bs=`dirname $1`/`basename $tmp .fastq`.seqcounts
o=${2:-$bs}
minLength=${3:-15}
maxLength=${4:-49}

if [[ $1 =~ \.gz$ ]]
then
    zcat $1 | awk '(NR%4==2)' | sort -S 2G | uniq -c | sed 's/^[ \t]*//' | awk -v minl="$minLength" -v maxl="$maxLength" -F" " '{if($1>1 && length($2)>=minl && length($2)<=maxl) print $2"\t"$1}' | sort > $o
else
    awk '(NR%4==2)' $1 | sort -S 2G | uniq -c | sed 's/^[ \t]*//' | awk -v minl="$minLength" -v maxl="$maxLength" -F" " '{if($1>1 && length($2)>=minl && length($2)<=maxl) print $2"\t"$1}' | sort > $o
fi
