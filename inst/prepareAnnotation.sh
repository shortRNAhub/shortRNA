#!/bin/bash

mature_miRNAs=ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz
rmsk=ftp://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/rmsk.txt.gz
gtRNAdbArchive=http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz
features=ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M14/gencode.vM14.chr_patch_hapl_scaff.annotation.gtf.gz
samtools=samtools

if [ $# -lt 3 ]
then
  echo "Usage: `basename $0` outputdir genome_fasta piRNA_precursors.gtf [gtRNAdbFile] [features.gtf.gz]"
  echo "Other sources files can be edited at the top of the script..."
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

mkdir -p $1
plcheck $1

fasta=`readlink -f $2`
piRNA=`readlink -f $3`
plcheck $fasta
plcheck $piRNA

if [ -z $5 ]; then
    echo "Fetching genomic features gtf and saving it in: $1/features.gtf.gz"
    wget -O $1/features.gtf.gz $features
    features=$1/features.gtf.gz
else
    features=`readlink -f $5`
fi
plcheck $features

if [ -z $4 ]; then
    cd $1
    echo "Fetching gtRNAdb archive:"
    wget http://gtrnadb.ucsc.edu/genomes/eukaryota/Mmusc10/mm10-tRNAs.tar.gz
    tar --one-top-level -xf mm10-tRNAs.tar.gz
    gtRNAdbBed=mm10-tRNAs/mm10-tRNAs.bed
else
    gtRNAdbBed=`readlink -f $4`
    cd $1
fi

plcheck $gtRNAdbBed

echo "
#####
# Preparing tRNA elements...
#####
"

# we fetch the genomic tRNA coordinates from gtRNAdb:
# we extract the spliced sequences, and then add the post-transcriptional modifications
bedtools getfasta -name -s -split -fi $fasta -bed $gtRNAdbBed | awk '
BEGIN{
split("ACGT", letters, "")
}
{
if(substr($0,1,1)==">"){
    split($0,a,"::")
    gsub("-","_",a[1])
}else{
    for(i=1; i<=4; i++){
        print a[1] letters[i] " " a[2]
        print letters[i] $0 "CCA"
    }
}
}
' > modified_tRNAs.fa

# we next create bed entries for each:
awk -F" " '{
if(substr($0,0,1)==">"){ 
    a=substr($1,2);
    b=a;
    gsub("_","-",b)
}else{
    print a"\t1\t"length($0)"\t"b"\t"substr(b,1,length(b)-1)"\ttRNA"
}
}' modified_tRNAs.fa > modified_tRNAs.artificialLocations.bed

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$4"\ttRNA"}' $gtRNAdbBed > tRNAs.genomicLocations.bed



echo "
#####
# Preparing miRNA elements...
#####
"

# we fetch reference miRNA sequences:
wget -qO- $mature_miRNAs | gunzip -c | grep -A 1 "Mus musculus" | grep -v "^--$" > mature.fa
# we create the 3' variations:
awk -F" " '
BEGIN{
    split("ACGT", letters, "")
}
{
if(substr($1,1,1)==">"){
    a=$1;
    gsub("-","-",a)
}else{
    gsub("U","T",$1)
    for(i=1; i<=4; i++){
        print a letters[i]
        print $1 letters[i]
    }
}
}
' mature.fa > modified_reference_miRNAs.fa

# we next create bed entries for each:
awk -F" " '
BEGIN{
    split("ACGT", letters, "")
}
{
if(substr($1,1,1)==">"){
    a=substr($1,2);
    b=a
    gsub("-","-",a)
}else{
    gsub("U","T",$1)
    for(i=1; i<=4; i++){
        print a letters[i] "\t" 1 "\t" length($1) "\t" b letters[i] "\t" b "\tmiRNA"
    }
}
}
' mature.fa > modified_reference_miRNAs.artificialLocations.bed



echo "
#####
# Preparing repeat masker elements...
#####
"

# we fetch the repeat masker elements:
# we'll rename tRNAs to be compatible with gtRNAdb

echo '
function revcomp(seq){
    ou=""
    for (i=length(seq);i>0;i--){
        ou=ou""a[substr(seq,i,1)]
    }
    return ou
}
BEGIN{
a["T"]="A";a["A"]="T";a["C"]="G";a["G"]="C";a["N"]="N";a["R"]="Y";a["Y"]="R"
}
{
if( $5=="tRNA" ){
    gname=substr($4,1,8)"-"revcomp(substr($4,10,3))
    print $1"\t"$2"\t"$3"\t"gname"\t"$6"\tputative_tRNA"
}else{
    if( $5!="Simple_repeat" && 
    $5!="Low_complexity" &&
    $5!="DNA" && $5!="DNA?" &&
    $5!="Satellite" && 
    $5!="Unknown" ){
        print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$5
    }
}
}' > clean_rmsk_table.awk

wget -qO- $rmsk | gunzip -c | cut -f 6-8,11-13 | awk -f clean_rmsk_table.awk > rmsk.bed

plcheck rmsk.bed



echo "
#####
# Preparing genes and piRNA precursors...
#####
"


# we prepare the following scripts for extracting the other elements

echo '
# For features with a single exon, the prepareElements.awk script will generate 
# duplicate entries (one for the precursor, one for the exon), so the following 
# script removes them (assuming they come one after the other, i.e. a sorted gtf)
BEGIN { lastp="" }
{ 
tmp=$1"\t"$2"\t"$3"\t"$4"\t"$5;
if(tmp != lastp && lastp!=""){
    print lastp"\tprecursor"
    lastp=""
}
if($6=="precursor"){
    lastp=tmp
}else{
    print $0
}
}
END { if(lastp!="") print lastp"\tprecursor" }
' > removeUnsplicedPrecursors.awk


echo '
{
    gid=".";
    gname=".";
    gtype=".";
    split($9,b,"%"); 
    for(i=1;i<=length(b);i++){ 
        if(b[i]=="gene_id") gid=b[i+1];
        if(b[i]=="gene_name") gname=b[i+1];
        if(b[i]=="gene_type") gtype=b[i+1];
        if(b[i]=="transcript_type" && gtype==".") gtype=b[i+1];
    }
    if(gname==".") gname=gid;
    sub("-miR-","-mir-",gname)
    if($3=="transcript" && gtype != "piRNA_precursor") gtype="precursor"
    print $1"\t"$4"\t"$5"\t"gid"\t"gname"\t"gtype
}' > prepareElements.awk

gunzip -c $features | grep -v "^#" | awk '{
if( $3=="transcript" || $3=="exon" ) print $0}' | sed 's/ //g' | sed 's/;//g' | tr '"' "%" | awk -F"\t" -f prepareElements.awk | awk -F"\t" -f removeUnsplicedPrecursors.awk | uniq > tmpelements.bed

# we assume to have the piRNA file:
# piRNA_precursors_zamore_mm10.gtf
awk '{if( $3=="transcript" || $3=="exon" ) print $0}' $piRNA | sed 's/ //g' | sed 's/;//g' | tr '"' "%" | awk -F"\t" -f prepareElements.awk | uniq >> tmpelements.bed

plcheck tmpelements.bed

cat tmpelements.bed *artificialLocations.bed | sortBed -i - | uniq > elements.bed

plcheck elements.bed

rm tmpelements.bed

rm *.awk

echo "
#####
# Now concatenating the artificial genome...
#####
"
cat $fasta modified_reference_miRNAs.fa modified_tRNAs.fa > customGenome.fa
$samtools faidx customGenome.fa
plcheck customGenome.fa.fai

echo "
#####
# Done!
#####

The annotation files were saved in:
$PWD

The following files will be needed to create the alignment indexes:
customGenome.fa
customGenome.fa.fai

The following file will be needed to run the read assignment pipeline:
elements.bed

The other files in the folder won't be needed anymore, but were left for conveniance.
"

