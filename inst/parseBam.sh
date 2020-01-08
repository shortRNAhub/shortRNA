#!/bin/bash
if [ $# -lt 2 ]
then
    echo "Usage: `basename $0` align.bam elements.bed [output.srcs] [samtools] [bedtools]"
    exit
fi
plcheck ()
{
    for iter in $@; do
        if [ ! -r $iter ]; then
            echo "ERROR: $iter not found! Aborting..."
            exit 0
        else
            if [ ! -s $iter ]; then
                echo "ERROR: $iter is empty! Aborting..."
                exit 0
            fi
        fi
    done
}
plcheck $1
plcheck $2
tmpo=`dirname $1`/seqs.srcs
osrcs=${3:-$tmpo}
samtools=${4:-samtools}
bedtools=${5:-bedtools}

$samtools view $1 | paste - <($samtools view $1 | cut -f 6 | egrep -n -o '([0-9]+)M' | sed 's/M//g' | awk -F":" 'BEGIN { ln=0 } {if($1==ln){ s+=$2 }else{ if(ln!=0){ print s }; ln=$1; s=$2 } } END { print s }') |  awk -F "\t" '{if($3!="*"){ print $3"\t"$4"\t"($4+$NF-1)"\t"$10"\t"$6"\t"$NF }}' | $bedtools sort -i - | $bedtools intersect -wao -f 0.5 -a - -b $2 | awk -F "\t" '{
 if($7=="."){ 
    print $4"\t"$5"\t"$1":"$2"\tNA\tNA\tNA\tNA\tNA\t"$6;
 }else{
     print $4"\t"$5"\t"$1":"$2"\t"($2-$8)"\t"($9-$8)"\t"$10"\t"$11"\t"$12"\t"(1+$13);
 }
}' | sort | uniq > $osrcs
$samtools view -f 4 $1 | cut -f 10 | awk '{ print $1"\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"}'>> $osrcs
