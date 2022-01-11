#pipeline in shell for ATACseq footprinting

organism=$1
fasta=$2
blacklist=$3
gtf=$4
motifs=$5
output=$6
macs=$7
mac="--nomodel --shift -100 --extsize 200 --broad"
