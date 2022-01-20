#!/bin/bash

#pipeline in shell for ATACseq footprinting
TEMP=`getopt -o vdm: --long gsize:,extsize:,shifts:,broad:,nomodel \
    -n './submit_atacseq' -- "$@"`

   if [ $? != 0 ]; then
       echo "Unrecognized argument. Possible arguments: --gize, --extsize, shifts, --broad true, --broad false, --nomodel true, --nomodel false." >&2 ; exit 1 ; 
   fi
       eval set -- "$TEMP"

        gsize=2620345972
        extsize=200
        shifts=100
        broad=true
        nomodel=true
        
    while true; do
        case "$1" in
            --gsize ) gsize="$2"; shift 2 ;;
            --extsize ) extsize="$2"; shift 2 ;;
            --shifts ) shifts="$2"; shift 2 ;;
            --broad ) broad="$2"; shift 2 ;;
            --nomodel ) nomodel="$2"; shift 2 ;;
            -- ) shift; break ;;
            * ) break ;;
        esac
    done

    if [ $broad != true ] || [ $braod != false ]; then
        echo "broad can only be set to true or false"
        echo "example: --broad true"
        echo "example: --broad false"
        exit 1
    fi

    if [ $nomodel != true ] || [ $nomodel != false ]; then
        echo "broad can only be set to true or false"
        echo "example: --broad true"
        echo "example: --broad false"
        exit 1
    fi
    
bam_path=$1
genome_path=$2
genome_folder="$(dirname "${genome_path}")"
echo "$genome_folder"

code_directory=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

#organism=$1
#fasta=$2
#blacklist=$3
#gtf=$4
#motifs=$5
#output=$6
#macs=$7
#mac="--nomodel --shift -100 --extsize 200 --broad"

#FASTA_PATH=/oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/genomes/mm9/mm9_bgzip.fa.gz
#BLACKLIST_PATH=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/blacklist/mm9-blacklist.bed.gz
#GTF_PATH=/oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/genomes/mm9/GTF/gencode.vM1.annotation.gtf.gz
#MOTIFS_DIR=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/TOBIAS_snakemake/data/individual_motifs
#OUTPUT_DIR=/oak/stanford/groups/smontgom/kameronr/ATACseq/output

module load samtools/1.9

##################################################################################################################################
########################################---STEP 1: CREATE COVERAGE BED FILES---################################################### 
##################################################################################################################################
#get chromosomes available in fasta (fasta chroms) -- because it's needed for making the bigwig tracks.
#OUTPUT_DIR=/oak/stanford/groups/smontgom/kameronr/ATACseq/test
#this next line might not work? seems like a python command? https://www.biostars.org/p/173963/

if [ ! -f $genome_folder/chromsizes.bed ]; then
    samtools faidx mm9_bgzip.fa.gz
    cut -f1,2 mm9_bgzip.fa.gz.fai > sizes.genome

    awk '{{ print $1, 0, $2 }}' chromsizes.txt > chromsizes.bed
    echo "creation of chromosome-based coverage bed files complete"
else
    echo "chromosome-based coverage bed files have already been created"
fi


##################################################################################################################################
#####################################---STEP 2: SORT, MERGE, AND INDEX BAM FILES---############################################### 
##########################################---CREATE COVERAGE BIGWIG TRACK---######################################################
############################################---PEAK CALLING WITH MACS2---#########################################################
##################################################################################################################################
cd $bam_path 
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

bed_file=$(find "$bam_path/peak_calling" -type f | grep "raw.bed" | sort -u | wc -l)

if [ $bed_file -lt 1 ]; then
         sbatch -o "${bam_path}/Logs/%A_%a.log" `#put into log` \
        -a "1-${array_length}" `#initiate job array equal to the number of bam files` \
        -W `#indicates to the script not to move on until the sbatch operation is complete` \
            "${code_directory}/sort.sh" \
            $bam_path $gsize $extsize $shifts $broad $nomodel
        
        wait
    else
        echo "sorting, merging, and indexing of files already completed"
fi

##################################################################################################################################
#########################---STEP 4: PEAK PROCESSING: REDUCE GENOMIC LOCATION COLUMNS AND SORT---################################## 
##################################################################################################################################
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




