# Parameters
input_directory="/oak/stanford/groups/sjaiswal/Herra/Twist_AGS_Batch_2a/FASTQ"
#input_file_list=
output_directory="/oak/stanford/scg/lab_sjaiswal/dnachun/projects/TWIST/new_test_output"
#bam_extension=
fastq_extension="fastq.gz"
#cloud_mode=false
slurm_mode=true
slurm_runtime="3:00:00"
slurm_account="sjaiswal"
slurm_ncpus="8"
slurm_memory="64G"
slurm_jobname="CHIP_variant_call"
assembly="GRCh38"
run_mutect=true
mutect_bam_output=true
#split_intervals=false
run_varscan=true
varscan_min_coverage="10"
varscan_min_var_freq="0.001"
varscan_max_pvalue="0.1"
run_haplotypecaller=true
#n_jobs=1
#realign=false

# References
reference_genome="/oak/stanford/groups/smontgom/maurertm/ADRC/Cloud_Testing_Folder/Params/GRCh38.p12.genome.u2af1l5_mask.fa"
gnomad_reference="/oak/stanford/groups/smontgom/maurertm/ADRC/Cloud_Testing_Folder/Params/gnomad.genomes.v3.1.2.sites.maf05.vcf.bgz"
funcotator_sources="/oak/stanford/groups/smontgom/maurertm/ADRC/Cloud_Testing_Folder/Params/funcotator_dataSources.v1.6.20190124s"
interval_list="/oak/stanford/groups/smontgom/maurertm/ADRC/Cloud_Testing_Folder/Params/CHIP_exons.interval_list"
#normal_bam=
#normal_pileups=
#sequence_dictionary=
#mpileup_interval_bed=
transcript_list="/oak/stanford/groups/smontgom/maurertm/ADRC/Cloud_Testing_Folder/Params/chip_transcript_list.txt"
annovarroot="/labs/sjaiswal/tools/annovar/"
germline_snps="/labs/sjaiswal/workflows/BWA_mutect_twist/twist_snps.bed"
