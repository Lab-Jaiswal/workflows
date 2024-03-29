input_directory="/oak/stanford/groups/sjaiswal/Herra/30-909057548/00_fastq"
input_file_list=none
output_directory="/oak/stanford/groups/sjaiswal/Herra/30-909057548/output"
bam_extension=none
fastq_extension=fastq.gz
slurm_mode=true
slurm_runtime=3:00:00
slurm_account=sjaiswal
slurm_partition=batch
slurm_ncpus=1
slurm_memory=64G
slurm_jobname=CHIP_variant_call
assembly=GRCh38
run_mutect=true
mutect_bam_output=true
split_intervals=false
run_funcotator=true
run_varscan=true
run_annovar=true
varscan_min_coverage=10
varscan_min_var_freq=0.001
varscan_max_pvalue=0.1
run_haplotypecaller=false
realign=false

reference_genome="/oak/stanford/groups/sjaiswal/references/reference_genome/GRCh38.p12.genome.u2af1l5_mask.fa"
sequence_dictionary="/oak/stanford/groups/sjaiswal/references/reference_genome/GRCh38.p12.genome.u2af1l5_mask.fa.dict"
gnomad_reference="/oak/stanford/groups/sjaiswal/references/gnomad.genomes.v3.1.2.sites.maf05.vcf.bgz"
funcotator_sources="/oak/stanford/groups/sjaiswal/references/funcotator_dataSources.v1.6.20190124s"
interval_list="/oak/stanford/groups/sjaiswal/references/CHIP_exons_updated020223.interval_list"
transcript_list="/oak/stanford/groups/sjaiswal/references/chip_transcript_list.txt"
whitelist="/oak/stanford/groups/sjaiswal/references/variant_whitelist_cleaned.xlsx"
mpileup_interval_bed="/oak/stanford/groups/sjaiswal/references/CHIP_exons_updated020223.bed"
annovarroot="/labs/sjaiswal/tools/annovar/"
