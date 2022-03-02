

#https://github.molgen.mpg.de/loosolab/TOBIAS_snakemake/blob/master/snakefiles/preprocessing.snake

#https://github.molgen.mpg.de/loosolab/TOBIAS_snakemake/blob/master/snakefiles/footprinting.snake


  organism: human                           #mouse/human
  fasta: data/genome_chr4.fa.gz             #.fasta-file containing organism genome
  blacklist: data/blacklist_chr4.bed        #.bed-file containing blacklisted regions
  gtf: data/genes_chr4.gtf                  #.gtf-file for annotation of peaks
  motifs: data/individual_motifs/*          #motifs (directory with files or individual files in MEME/JASPAR/PFM format)
  output: test_output                       #output directory 

  macs: "--nomodel --shift -100 --extsize 200 --broad"

#blacklist https://github.com/Boyle-Lab/Blacklist/tree/master/lists



#Step 1. Combine bam files of replicate samples
#[Merge sample bam files to condition bam file (using samtools merge)]

bam = os.path.join(OUTPUTDIR, "mapping", "{condition}.bam"),
bai = os.path.join(OUTPUTDIR, "mapping", "{condition}.bam.bai")
output_dir = os.path.dirname(output.bam)

output_dir=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/

#location where bam files are

#sort bam files first
#then merge the sorted bam files
#

merge_bams = []
for bam in input:
	prefix = os.path.basename(bam)
	temp_prefix = output_dir + prefix
	merge_bams.append(os.path.join(output_dir, prefix + ".sorted"))
	shell("samtools sort {0} -o {1} -@ {2} -T {3}".format(bam, merge_bams[-1], threads, temp_prefix))

if len(merge_bams) > 1:
	shell("samtools merge -@ {threads} {output.bam} " + " ".join(merge_bams))
	shell("samtools index {output.bam}")
else:
	shell("samtools sort -o {output.bam} {input}")
	shell("samtools index {output.bam}")

#Step 2. Create coverage bigwig track


#Step 3. Peak calling with MACS2 (get broad peaks?) (why broad peaks and not narrow peaks?). because ATACseq includes more broadpeaks by nature than ChIPseq?


#Step 4a. Peak processing: reduce to genomic location columns and sort

#Step 4b. merge peaks per condition

#Step 4c. remove blacklisted regions and add unique peak ids


#inputs
  WT_NT: [/oak/stanford/groups/sjaiswal/kameronr/ATACseq/ATAC_tet2_WT_NT*.bam]  #list of .bam-files
  WT_LDL: [/oak/stanford/groups/sjaiswal/kameronr/ATACseq/ATAC_tet2_WT_LDL*.bam]   #list of .bam-files
  Tet2_KO_NT: [/oak/stanford/groups/sjaiswal/kameronr/ATACseq/ATAC_tet2_KO_NT*.bam]  #list of .bam-files
  Tet2_KO_LDL: [/oak/stanford/groups/sjaiswal/kameronr/ATACseq/ATAC_tet2_KO_LDL*.bam]   #list of .bam-files

  organism: mouse                           #mouse/human
  fasta: /oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/genomes/mm9/mm9_bgzip.fa.gz    #.fasta-file containing organism genome
  blacklist: /oak/stanford/groups/sjaiswal/kameronr/ATACseq/blacklist/mm9-blacklist.bed.gz      #.bed-file containing blacklisted regions
  gtf: /oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/genomes/mm9/GTF/gencode.vM1.annotation.gtf.gz      #.gtf-file for annotation of peaks
  motifs: data/individual_motifs/*          #motifs (directory with files or individual files in MEME/JASPAR/PFM format)
  output: /oak/stanford/groups/sjaiswal/kameronr/ATACseq/TOBIAS_snakemake/tet2_pilot_output  

FASTA_PATH=/oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/genomes/mm9/mm9_bgzip.fa.gz
BLACKLIST_PATH=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/blacklist/mm9-blacklist.bed.gz
GTF_PATH=/oak/stanford/groups/sjaiswal/kameronr/sjaiswal_old/genomes/mm9/GTF/gencode.vM1.annotation.gtf.gz
MOTIFS_DIR=/oak/stanford/groups/sjaiswal/kameronr/ATACseq/TOBIAS_snakemake/data/individual_motifs
OUTPUT_DIR=/oak/stanford/groups/smontgom/kameronr/ATACseq/output


#get chromosomes available in fasta (fasta chroms) -- for what/why?
#OUTPUT_DIR=/oak/stanford/groups/smontgom/kameronr/ATACseq/test
module load samtools/1.9
#this next line might not work? seems like a python command? https://www.biostars.org/p/173963/
mkdir $OUTPUT_DIR/flatfiles
faidx $FASTA_PATH -i chromsizes > $OUTPUT_DIR/flatfiles/chromsizes.txt
#samtools faidx input.fa
#cut -f1,2 input.fa.fai > sizes.genome
awk '{{ print $1\"\t\"0\"\t\"$2 }}' $OUTPUT_DIR/flatfiles/chromsizes.txt > $OUTPUT_DIR/flatfiles/chromsizes.bed

#merge sample bam files to condition bam file
#if only one sample per condition, copy sample bam to merged bam file name for further processing
#Combine bam files of replicate samples
#join individual bamfiles from the same condition

$OUTPUTDIR/mapping/{condition}.bam.bai
output_dir=$(dirname $OUTPUTDIR/mapping/{condition}.bam)
threads_count=10
files=??? #the folder to all the bam files
for bam in files/*.bam:
	prefix=$(basename bam)
	temp_prefix=$ouput_dir/$prefix
	samtools sort bam -o $output_dir/$prefix.sorted -@ $threads_count -T $temp_prefix

#if more than 1 bam file
samtools merge -@ $threads_count $OUTPUTDIR/mapping/{condition}.bam *.bam
samtools index $OUTPUTDIR/mapping/{condition}.bam

#else:
samtools sort -o $OUTPUTDIR/mapping/{condition}.bam files/*.bam
samtools index $OUTPUTDIR/mapping/{condition}.bam



#create coverage bigwig track
bedtools genomecov -ibam {input.bam} -bg | sort -k1,1 -k2,2 -T $OUTPUTDIR/coverage > $OUTPUTDIR/coverage/{condition}_coverage.bg
bedGraphToBigWig $OUTPUTDIR/coverage/{condition}_coverage.bg $OUTPUT_DIR/flatfiles/chromsizes.txt $OUTPUTDIR/coverage/{condition}_coverage.bw


#peak calling
#GRCm37	2620345972
#GRCm38	2652783500
#bams of condition and sample_id
gsize=2620345972 #is this correct?!
macs=$OUTPUTDIR/peak_calling/{condition}/{sample_id}_peaks.broadPeak
raw=$OUTPUTDIR/peak_calling/{condition}/{sample_id}_raw.bed
log=$OUTPUTDIR/logs/{condition}_{sample_id}_peak_calling.log

echo Running macs2 with .bam-file: {input}

macs2 callpeak -t {input} --name {sample_id} --outdir $OUTPUTDIR/peak_calling/{condition} --gsize $gsize --nomodel --shift -100 --extsize 200 --broad &> $OUTPUTDIR/logs/{condition}_{sample_id}_peak_calling.log
cp $OUTPUTDIR/peak_calling/{condition}/{sample_id}_peaks.broadPeak $OUTPUTDIR/peak_calling/{condition}/{sample_id}_raw.bed



#peak processing
#1. reduce to genomic location columns, remove blacklisted regions, sort, then merge peaks per condition.
for each condition OUTPUTDIR/peak_calling/{condition}/sample_id_raw.bed
blacklist=$BLACKLIST
whitelist=$OUTPUT_DIR/flatfiles/chromsizes.bed

cat OUTPUTDIR/peak_calling/{condition}/sample_id_raw.bed | cut -f1-3 | sort -k1,1 -k2,2n | bedtools merge -d 5 | bedtools substract -a - -b $blacklist -A | bedtools intersect -a - -b $whitelist -wa | awk '$1 !~ /[M]/' | awk '{{print $s\"\\t{wildcards.condition}\"}}' > $OUTPUTDIR/peak_calling/{condition}_union.bed
#excludes mitochondria chromosome (M)
#adds condition name to each peak



#2. Union peaks across all conditions.
$OUTPUTDIR/peak_calling/{condition}_union.bed for each conditon 
temp?
$OUTPUTDIR/peak_calling/all_merged.tmp
cat $OUTPUTDIR/peak_calling/{condition}_union.bed | sort -k1,1 -k2,2n | bedtools merge -d 5 -c 4 -o distinct > $OUTPUTDIR/peak_calling/all_merged.tmp


#get correct sorting of peak_names
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



#3. peak annotation. peaks per condition or across conditions, dependent on run info output

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






















