# Workflows
## BWA CHIP Workflow
### Set Up
Running this pipeline is done by:
1. Ensure that you have a directory containing only fastqs from your samples
  * If you do not have fastqs, but have bam files, please convert the files to fastqs using our [bam to fastq script](https://github.com/Lab-Jaiswal/workflows/tree/main/filetype_transformation/bam_to_fastq)
  * The fastq files must be aligned to the hg38 genome. If they are not, please convert them to bam files using our [bam to fastq script](https://github.com/Lab-Jaiswal/workflows/tree/main/filetype_transformation/bam_to_fastq) and then convert them back to fastqs using our [bam to fastq script](https://github.com/Lab-Jaiswal/workflows/tree/main/filetype_transformation/fastq_to_bam), which will align them to the hg38 genome.
2. Moving into the directory where you have housed the BWA_CHIP pipeline scripts:

    `cd /path/to/BWA/pipeline`
3. The script that submits the jobs on SCG for alignment and QC requires you to specify:
  * FASTQ filepath
  * output filepath
  * analysis_type (any combination of --mutect --varscan --haplotypecaller and --all)
 
    `./submit_BWA_CHIP FASTQ_DIRECTORY OUTPUT_DIRECTORY --analysis_type`
  * You may also specify optional arguments, as discussed below
  
### Mutect Optional Arguments

The BWA_CHIP workflow contains a variety of optional arguments, which allow great flexibility in calling CHIP.

#### Changing Default Values
Many of the genomes and sources in the pipeline have hardcoded default values thatn can be changed by specifying their corresponding optional argument
| Optional Argument    | Default Value                                                                              | 
| ---------------------|:------------------------------------------------------------------------------------------ |
| `--bwa_gref`          | "/oak/stanford/groups/sjaiswal/Herra/CHIP_Panel_AmpliSeq/GRCh38.p12.genome.u2af1l5_mask.fa"| 
| `--twist_snps`        | "/labs/sjaiswal/workflows/BWA_mutect_twist/twist_snps.bed"                                 | 
| `--funcotator_sources`|"/labs/sjaiswal/tools/funcotator/funcotator_dataSources.v1.6.20190124s"                     |   
| `--transcript_list`   | "/oak/stanford/groups/sjaiswal/Herra/CHIP_TWIST-PANEL_ATHEROMA/chip_transcript_list.txt"   | 
| `--assembly`          | "GRCh38"                                                                                   |

You can change any of these values by using the optional argument followed by your replacement value. For example:
* `--bwa_gref /path/to/new/bwa/gref` changes the bwa_gref variable from the default to the path you provided

#### Tumor Normal
In order to run mutect in [tumor-normal mode](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2), you need to 
add the flag --normal_sample followed *immediately* by the path to your normal sample

#### Filtering
| Optional Argument                               | Use Case                                                                                                                      | 
| ------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------| 
| `--twist`                                         | If you want the downstream R analysis to be run on a Twist panel (either the default one or one specified using --twist_snps) |
| `--filtered`                                      | If you want to only include frame shift, missense, nonsense, and splicesite mutations                                         | 
| `--intervals path/to/interval/file.interval_list` | If you would like to only call CHIP on specific intervals*                                                                    |   

* please note that the CHIP calls are currently only called on Chromosomes 1-22, X, and Y. If you need to analyze mutations outside of those chromosomes, please make an interval file and se --intervals 
#### Additional Arguments

| Optional Argument | Function                                                | 
|-------------------|:--------------------------------------------------------|
| `--log_name`        | To define the log name                                  |   
| `--no_funcotator`   | To indicate that funcotator annotation is not required  |   
| `--no_bam_out`      | To indicate that bam files do not need to be outputted  | 

### Varscan Optional Arguments
| Optional Argument        | Default Value | Example                                                 |   
| ------------------------ |:-------------:|:-------------------------------------------------------:| 
| `--min_coverage value`     |10            |`--min_coverage 100` to change the min_coverage to 100   |
| `--min_var_freq value`     |0.001         | `--min_var_freq 0.01` to change the min_var_freq to 0.01|  
| `--p_value value`          |0.1           | `--p_value .05` to change the p_value to 0.1            |  

# Create Panel of Normals using the UkBioBank
1. Use get_1000_youngest_with_crams_no_predispositions.R to get a list of the 1000 youngest with WES data and no CHIP predispositions (hx of malignant neoplasms, autoimmune disorders, or DMARD treatment
2. Upload the resulting file to DNANexus
3. dx run cloud_workstation -imax_session_length=10hr --ssh --brief -y --name "test"
4. Download the 1000 youngest file list using: dx download project-G5B07V8JPkg740v9GjfF9PzV:/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/top_1000_without_predis_with_crams
4. Get a list of all of the file paths using 
* dx ls --full project-G5B07V8JPkg740v9GjfF9PzV:/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/ | xargs -I {} bash -c 'array=( $(dx ls --full "project-G5B07V8JPkg740v9GjfF9PzV:{}") ); printf "%s\n" "${array[@]/#/{}/}"' > all_file_paths
* sed -e 's/\s/\\ /g' < all_file_paths > all_file_paths_fixed
5. Filter out all filepaths to only list the 1000 youngest using:
* grep -f top_1000_without_predis_with_crams < all_file_paths_fixed > filepaths_of_1000_youngest
* sed -e 's/\s/\\ /g' < filepaths_of_1000_youngest > filepaths_of_1000_youngest_fixed
* sed -e 's/^/project-G5B07V8JPkg740v9GjfF9PzV:/g' filepaths_of_1000_youngest_fixed > filepaths_of_1000_youngest_final
6. dx upload the filepaths_of_1000_youngest_final
* dx upload filepaths_of_1000_youngest_final --path  project-G5B07V8JPkg740v9GjfF9PzV:/Bulk/Exome\ sequences/Exome\ OQFE\ CRAM\ files/
7. dx download the 1000 youngest using the filepaths in step 5
* mkdir 1000_youngest
* cd 1000_youngest
* cat ../filepaths_of_1000_youngest_final | xargs -I % dx download %
8. upload the folder with the 1000 youngest

