#!/bin/bash
if [ $# -lt 1 ]
then
    echo "Usage: `basename $0` sample1.seqcounts [sampl2.seqcounts ...]"
    echo "Prints to stdout."
    exit
fi

printf "seq"

for f in $@; do
    printf "\t"`basename $f .seqcounts`
done
echo ""

`dirname $0`/full_outer_join.awk HEADER=0 $@ | awk '{keep=0; for(i=2;i<=NF;i++){ if($i != "NULL" && $i>3) keep++ }; if(keep>1) print $0}' | sed 's/NULL/0/g'
