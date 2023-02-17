# Parameters
input_directory="$HOME/inputs"
output_directory="$HOME/outputs"
bam_extension="cram"
slurm_mode=false
assembly="GRCh38"
run_mutect=true
mutect_bam_output=true
run_varscan=false
run_haplotypecaller=false

# References
reference_genome="$HOME/references/cloud_references/GRCh38_full_analysis_set_plus_decoy_hla.fa"
gnomad_genomes="$HOME/references/gnomad.genomes.v3.1.2.sites.maf05.vcf.bgz"
funcotator_sources="$HOME/references/funcotator_dataSources.v1.6.20190124s"
interval_list="$HOME/references/CHIP_exons.interval_list"
transcript_list="$HOME/references/chip_transcript_list.txt"
