#!/bin/bash

#pipeline in shell for ATACseq footprinting



organism=$1
fasta=$2
blacklist=$3
gtf=$4
motifs=$5
output=$6
macs=$7
mac="--nomodel --shift -100 --extsize 200 --broad"

code_directory=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

data_path=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/
bam_path="/home/maurertm/labs/maurertm/atac_seq/data"

bam_directory=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/
WT_NT_files=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/ATAC_tet2_WT_NT*.bam
WT_LDL_files=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/ATAC_tet2_WT_LDL*.bam
Tet2_KO_NT_file=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/ATAC_tet2_KO_NT*.bam
Tet2_KO_LDL_file=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/ATAC_tet2_KO_LDL*.bam

FASTA_PATH=/oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/genomes/mm9/mm9_bgzip.fa.gz
BLACKLIST_PATH=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/blacklist/mm9-blacklist.bed.gz
GTF_PATH=/oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/genomes/mm9/GTF/gencode.vM1.annotation.gtf.gz
MOTIFS_DIR=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/TOBIAS_snakemake/data/individual_motifs
OUTPUT_DIR=/oak/stanford/groups/smontgom/kameronr/ATACseq/output


module load samtools/1.9

#Step 1. Combine bam files of replicate samples
#1a. sort each bam file
#1b. merge bam files with the same condition (for any conditions that have more than 1 replicate)
#1c. index all resulting files (should be 1 final sorted bam for each experimental condition)
#code for merging bams for single condition

cd $bam_path 
#bams=$bam_dir/.*
bam_file="$bam_path/BAMs" #give a path to a file to store the paths to the bams files in $bam_directory

find "${bam_path}/" -type f `#list all files in ${fastq_directory}` | \
        grep ".*\.bam$" `#only keep files with FASTQ in name (case insensitive)` | \
        grep -v ".sorted.bam" `#remove sorted BAMs` | \
	grep -v ".merged.bam" `#remove merged BAMs` | \
        sed -e 's/\_Rep.*//g' `#remove _Rep1/2_treat_rep1.bam extension` | \
        sort -u  `#sort and remove duplicate names` > ${bam_file}
            #| \
        #head -n -1 > ${bam_list} `#remove the last line and generate a list of unique FASTQs`
    bam_array_length=$(wc -l < ${bam_file}) #get the number of FASTQs 
    echo "$bam_array_length"

number_bams=$(wc -l < "${bam_file}") #get the number of files
array_length="$number_bams"


if ! [ -d "$bam_path/Logs" ]; then
    mkdir -p "$bam_path/Logs"
fi

sorted=$(find "$bam_path/" -type f | grep "sorted" | sort -u | wc -l)

if [ $sorted -le 1 ]; then
         sbatch -o "${bam_path}/Logs/%A_%a.log" `#put into log` \
        -a "1-${array_length}" `#initiate job array equal to the number of bam files` \
        -W `#indicates to the script not to move on until the sbatch operation is complete` \
            "${code_directory}/sort.sh" \
            $bam_path
        
        wait
    else
        echo "sorted files found"
fi


#bams=WT_NT_files
#bam=$OUTPUT_DIR
#bam=os.path.join(OUTPUTDIR, "mapping", "{condition}.bam"),
#bai = os.path.join(OUTPUTDIR, "mapping", "{condition}.bam.bai")
#output_dir=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/
#input_
#output_name= samtools sort bam -o merge_bams[-1] -@ threads -T temp_prefix

#step 2. Create coverage bigwig track

#get chromosomes available in fasta (fasta chroms) -- because it's needed for making the bigwig tracks.
#OUTPUT_DIR=/oak/stanford/groups/smontgom/kameronr/ATACseq/test
#this next line might not work? seems like a python command? https://www.biostars.org/p/173963/
mkdir $OUTPUT_DIR/flatfiles
faidx $FASTA_PATH -i chromsizes > $OUTPUT_DIR/flatfiles/chromsizes.txt
#samtools faidx input.fa
#cut -f1,2 input.fa.fai > sizes.genome
awk '{{ print $1\"\t\"0\"\t\"$2 }}' $OUTPUT_DIR/flatfiles/chromsizes.txt > $OUTPUT_DIR/flatfiles/chromsizes.bed



#Make bigwig tracks
bedtools genomecov -ibam {input.bam} -bg | sort -k1,1 -k2,2 -T $OUTPUTDIR/coverage > $OUTPUTDIR/coverage/{condition}_coverage.bg
bedGraphToBigWig $OUTPUTDIR/coverage/{condition}_coverage.bg $OUTPUT_DIR/flatfiles/chromsizes.txt $OUTPUTDIR/coverage/{condition}_coverage.bw





#step 3. Peak calling with MACS2 (get broad peaks (in future can also try narrow peaks)) (why broad peaks and not narrow peaks?). because ATACseq includes more broadpeaks by nature than ChIPseq?
#bams of condition and sample_id
gsize=2620345972 #is this correct?!
macs=$OUTPUTDIR/peak_calling/{condition}/{sample_id}_peaks.broadPeak
raw=$OUTPUTDIR/peak_calling/{condition}/{sample_id}_raw.bed
log=$OUTPUTDIR/logs/{condition}_{sample_id}_peak_calling.log

echo Running macs2 with .bam-file: {input}

macs2 callpeak -t {input} --name {sample_id} --outdir $OUTPUTDIR/peak_calling/{condition} --gsize $gsize --nomodel --shift -100 --extsize 200 --broad &> $OUTPUTDIR/logs/{condition}_{sample_id}_peak_calling.log
cp $OUTPUTDIR/peak_calling/{condition}/{sample_id}_peaks.broadPeak $OUTPUTDIR/peak_calling/{condition}/{sample_id}_raw.bed


#Step 4. Peak processing: reduce to genomic location columns and sort

#4a. reduce to genomic location columns, remove blacklisted regions, sort, then merge peaks per condition.
for each condition OUTPUTDIR/peak_calling/{condition}/sample_id_raw.bed
blacklist=$BLACKLIST
whitelist=$OUTPUT_DIR/flatfiles/chromsizes.bed

cat OUTPUTDIR/peak_calling/{condition}/sample_id_raw.bed | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -d 5 | bedtools substract -a - -b $blacklist -A | bedtools intersect -a - -b $whitelist -wa | awk '$1 !~ /[M]/' | awk '{{print $s\"\\t{wildcards.condition}\"}}' > $OUTPUTDIR/peak_calling/{condition}_union.bed
#excludes mitochondria chromosome (M)
#adds condition name to each peak

#4b. Union peaks across all conditions.
$OUTPUTDIR/peak_calling/{condition}_union.bed for each conditon 
temp?
$OUTPUTDIR/peak_calling/all_merged.tmp
cat $OUTPUTDIR/peak_calling/{condition}_union.bed | sort -k1,1 -k2,2n | bedtools merge -d 5 -c 4 -o distinct > $OUTPUTDIR/peak_calling/all_merged.tmp


#4c. get correct sorting of peak_names
$OUTPUTDIR/peak_calling/all_merged.tmp
peaks=$OUTPUTDIR/peak_calling/all_merged.bed

python script:
{
		out = open($OUTPUTDIR/peak_calling/all_merged.bed, "w")
		with open($OUTPUTDIR/peak_calling/all_merged.tmp) as f:
			for line in f:
				columns = line.rstrip().split("\t")

				#Sort according to condition names
				peak_ids = columns[3].split(",")
				columns[3] = ",".join(sorted(peak_ids, key= lambda x: CONDITION_IDS.index(x)))

				out.write("\t".join(columns) + "\n")
		out.close()
}

#4d. peak annotation. peaks per condition or across conditions, dependent on run info output

#need to make .config file for uropa to use.
python script:
{
		import json
		config = {"queries":[
					{"feature":"gene", "feature.anchor":"start", "distance":[10000,1000], "filter_attribute":"gene_biotype", "attribute_values":"protein_coding", "name":"protein_coding_promoter"},
					{"feature":"gene", "distance":1, "filter_attribute":"gene_biotype", "attribute_values":"protein_coding", "internals":0.1, "name":"protein_coding_internal"},
					{"feature":"gene", "feature.anchor":"start", "distance":[10000,1000], "name":"any_promoter"},
					{"feature":"gene", "distance":1, "internals":0.1, "name":"any_internal"},
					{"feature":"gene", "distance":[50000, 50000], "name":"distal_enhancer"},
					],
				"show_attributes":["gene_biotype", "gene_id", "gene_name"],
				"priority":"True"
				}

		config["gtf"] = $GTF_PATH
		config["bed"] = $OUTPUTDIR/peak_calling/all_merged.bed

		string_config = json.dumps(config, indent=4)

		config_file = open($OUTPUTDIR/peak_annotation/all_merged_annotated.config, "w")
		config_file.write(string_config)
		config_file.close()
}



$OUTPUTDIR/peak_annotation/all_merged_annotated.config
log=$OUTPUTDIR/logs/uropa.log
prefix=$OUTPUTDIR/peak_annotation/all_merged_annotated
finalhits=$OUTPUTDIR/peak_annotation/all_merged_annotated_finalhits.txt
finalhits_sub=$OUTPUTDIR/peak_annotation/all_merged_annotated_finalhits_sub.txt
peaks=$OUTPUTDIR/peak_annotation/all_merged_annotated.bed
header=$OUTPUTDIR/peak_annotation/all_merged_annotated_header.txt

uropa --input $OUTPUTDIR/peak_annotation/all_merged_annotated.config --prefix $OUTPUTDIR/peak_annotation/all_merged_annotated --threads $threads_count --log $OUTPUTDIR/logs/uropa.log; 
cut -f 1-4,7-13,16-19 $OUTPUTDIR/peak_annotation/all_merged_annotated_finalhits.txt > $OUTPUTDIR/peak_annotation/all_merged_annotated_finalhits_sub.txt;  #Get a subset of columns
head -n 1 $OUTPUTDIR/peak_annotation/all_merged_annotated_finalhits_sub.txt > $OUTPUTDIR/peak_annotation/all_merged_annotated_header.txt;  #header
tail -n +2 $OUTPUTDIR/peak_annotation/all_merged_annotated_finalhits_sub.txt > $OUTPUTDIR/peak_annotation/all_merged_annotated.bed #bedlines


#could later add expression information to each peaks if needed?




