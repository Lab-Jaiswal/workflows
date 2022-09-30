#!/bin/bash

#SBATCH --time=72:00:00
#SBATCH --account=sjaiswal
#SBATCH --cpus-per-task=8
#SBATCH --mem=124GB
#SBATCH --job-name=CHIP_variant_call

set -u

##################################################################################################################################
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
PARENT_DIRECTORY=$1
OUTPUT_DIRECTORY="$2"
MIN_COVERAGE="$3"
MIN_VAR_FREQ="$4"
P_VALUE="$5"
GET_MUTECT="$6"
GET_VARSCAN="$7"
GET_HAPLOTYPE="$8"
USE_BAM="$9"
INTERVALS_FILENAME="${10}"
NORMAL_SAMPLE_FILENAME="${11}"
CODE_DIRECTORY="${12}"
PARAMETER_FILE=${13}
BWA_GREF_FILENAME=${14}
TWIST_SNPS=${15}
ASSEMBLY=${16}
FUNCOTATOR_SOURCES_FILENAME=${17}    
TRANSCRIPT_LIST_FILENAME=${18}
FILTERED=${19}
RUN_FUNCOTATOR=${20}
BAM_OUT=${21}
NORMAL_PILEUPS_FILENAME="${22}"
MODE=${23}
DOCKER_IMAGE=${24}
CONTAINER_ENGINE=${25}
SPLIT_BY_CHR=${26}
SEQUENCE_DICT_FILENAME=${27}
CHR_INTERVALS_FILENAME=${28}
GNOMAD_GENOMES_FILENAME=${29}
RUN_MUTECT=${30}

if [[ $CONTAINER_ENGINE == "singularity" ]]; then
    gatk_command="singularity run instance://gatk_container gatk"
fi

if [[ $MODE == "slurm" ]]; then
    LINE_NUMBER=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM
else
    LINE_NUMBER=${TASK_ID}
fi

if [[ $LINE_NUMBER == 1 ]]; then
    echo "PARENT_DIRECTORY=$1 \
        OUTPUT_DIRECTORY=$2 \
        MIN_COVERAGE=$3 \
        MIN_VAR_FREQ=$4 \
        P_VALUE=$5 \
        GET_MUTECT=$6 \
        GET_VARSCAN=$7 \
        GET_HAPLOTYPE=$8 \
        USE_BAM=$9 \
        INTERVALS_FILENAME=${10} \
        NORMAL_SAMPLE_FILENAME=${11} \
        CODE_DIRECTORY=${12} \
        PARAMETER_FILE=${13} \
        BWA_GREF_FILENAME=${14} \
        TWIST_SNPS=${15} \
        ASSEMBLY=${16} \
        FUNCOTATOR_SOURCES_FILENAME=${17} \ 
        TRANSCRIPT_LIST_FILENAME=${18} \
        FILTERED=${19} \
        RUN_FUNCOTATOR=${20} \
        BAM_OUT=${21}
        NORMAL_PILEUPS_FILENAME=${22}
        MODE=${23}
        DOCKER_IMAGE=${24}
        CONTAINER_ENGINE=${25}
        SPLIT_BY_CHR=${26}
        SEQUENCE_DICT_FILENAME=${27}
        CHR_INTERVALS_FILENAME=${28}
        GNOMAD_GENOMES_FILENAME=${29}
        RUN_MUTECT=${30}" >> $PARAMETER_FILE
fi

if [[ $RUN_MUTECT == false ]]; then
    NORMAL_SAMPLE_FILENAME=false
    NORMAL_PILEUPS_FILENAME=false
fi

if [ $USE_BAM = false ]; then
    ARRAY_FILE="${PARENT_DIRECTORY}/fastq_files" #provide path to file containing list of fastq files
else 
    ARRAY_FILE="${PARENT_DIRECTORY}/bam_files" #provide path to file containing list of fastq files
fi

if [[ $USE_BAM == false ]]; then
    ARRAY_FILE="${PARENT_DIRECTORY}/normal_sample_path" #provide path to file containing list of fastq files
fi


ARRAY_PREFIX="$(sed "${LINE_NUMBER}q; d" "${ARRAY_FILE}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID

#CHR_INTERVALS="${CODE_DIRECTORY}/whole_genome_intervals.interval_list"
#CHR_INTERVALS="/oak/stanford/groups/smontgom/maurertm/ADRC/workflows/Interval_filtering/chr21_chr22.interval_list"

#Run script and save files in the same location of the fastq files

FILENAME=$(basename "${ARRAY_PREFIX}")
SAMPLE_NAME="${FILENAME}_${ASSEMBLY}"
if [ $USE_BAM = false ]; then
   R1="${ARRAY_PREFIX}_R1_001.fastq.gz"
   R2="${ARRAY_PREFIX}_R2_001.fastq.gz"
   READGROUP="@RG\tID:${FILENAME}\tLB:${FILENAME}\tPL:illumina\tSM:${FILENAME}"
   echo "$READGROUP"
   MUTECT_INPUT_FILENAME="${FILENAME}_${ASSEMBLY}"

   echo "FILENAME: $FILENAME"

   ${CODE_DIRECTORY}/fastq_to_bam.sh $SAMPLE_NAME $READGROUP $BWA_GREF_FILENAME $R1 $R2 $PARAMETER_FILE
else
    MUTECT_INPUT_FILENAME="${ARRAY_PREFIX}"
fi


OUTPUT_DIRECTORY=${OUTPUT_DIRECTORY}/${SAMPLE_NAME}
mkdir -p $OUTPUT_DIRECTORY

cd "${OUTPUT_DIRECTORY}" || exit

echo "MODE: $MODE !!!!!!!!!!!!!!!"
echo "RUN_FUNCOTATOR: $RUN_FUNCOTATOR"

if [[ $MODE == "slurm" ]]; then
    temp_dir="/tmp/sjaiswal"
    INPUT="${TMPDIR}/Inputs"
    OUTPUTS="${TMPDIR}/Outputs"
    PARAMS="${TMPDIR}/Params"

    mkdir -p $INPUT
    mkdir -p $OUTPUTS
    mkdir -p $PARAMS

    echo "temp_path is: " $temp_dir
    echo "copying FASTQs from the data path..." #copy data from data_path to the temp_path

    MUTECT_INPUT_DIRNAME=$(dirname "${MUTECT_INPUT_FILENAME}")
    MUTECT_INPUT_BASENAME=$(basename "${MUTECT_INPUT_FILENAME}")
    MUTECT_INPUT="${INPUT}/${MUTECT_INPUT_BASENAME}"

    if [[ $LINE_NUMBER == 1 ]]; then
        echo "MUTECT_INPUT_BASENAME: $MUTECT_INPUT_BASENAME
              MUTECT_INPUT: $MUTECT_INPUT" >> $PARAMETER_FILE
    fi

    rsync -vurPhlt "${MUTECT_INPUT_FILENAME}.bam" "${MUTECT_INPUT}.bam"
    rsync -vurPhlt "${MUTECT_INPUT_FILENAME}.bai" "${MUTECT_INPUT}.bai"
    MUTECT_INPUT_PATTERN=${MUTECT_INPUT_BASENAME%.*}
    rsync -vurPhlt $MUTECT_INPUT_DIRNAME/$MUTECT_INPUT_PATTERN* $INPUT/

    if [[ $RUN_MUTECT == true ]]; then
        BWA_GREF_DIRNAME=$(dirname "${BWA_GREF_FILENAME}")
        BWA_GREF_BASENAME=$(basename "${BWA_GREF_FILENAME}")
        BWA_GREF="${PARAMS}/${BWA_GREF_BASENAME}"

        rsync -vurPhlt "$BWA_GREF_FILENAME" $BWA_GREF
        BWA_GREF_PATTERN=${BWA_GREF_BASENAME%.*}
        rsync -vurPhlt $BWA_GREF_DIRNAME/$BWA_GREF_PATTERN* $PARAMS/
    
    else
        BWA_GREF=$BWA_GREF_FILENAME
    fi

    if [[ $SPLIT_BY_CHR == true ]]; then
        SEQUENCE_DICT_DIRNAME=$(dirname "${SEQUENCE_DICT_FILENAME}")
        SEQUENCE_DICT_BASENAME=$(basename "${SEQUENCE_DICT_FILENAME}")
        SEQUENCE_DICT="${PARAMS}/${SEQUENCE_DICT_BASENAME}"

        rsync -vurPhlt "$SEQUENCE_DICT_FILENAME" $SEQUENCE_DICT
        SEQUENCE_DICT_PATTERN=${SEQUENCE_DICT_BASENAME%.*}
        rsync -vurPhlt $SEQUENCE_DICT_DIRNAME/$SEQUENCE_DICT_PATTERN* $PARAMS/
    fi
    
    INTERVALS_FILE_DIRNAME=$(dirname "${INTERVALS_FILENAME}")
    INTERVALS_FILE_BASENAME=$(basename "${INTERVALS_FILENAME}")
    INTERVALS_FILE="${PARAMS}/${INTERVALS_FILE_BASENAME}"

    rsync -vurPhlt "$INTERVALS_FILENAME" $INTERVALS_FILE
    INTERVALS_FILE_PATTERN=${INTERVALS_FILE_BASENAME%.*}
    rsync -vurPhlt $INTERVALS_FILE_DIRNAME/$INTERVALS_FILE_PATTERN* $PARAMS/

    CHR_INTERVALS_DIRNAME=$(dirname "${CHR_INTERVALS_FILENAME}")
    CHR_INTERVALS_BASENAME=$(basename "${CHR_INTERVALS_FILENAME}")
    CHR_INTERVALS="${PARAMS}/${CHR_INTERVALS_BASENAME}"

    rsync -vurPhlt "$CHR_INTERVALS_FILENAME" $CHR_INTERVALS
    CHR_INTERVALS_PATTERN=${CHR_INTERVALS_BASENAME%.*}
    rsync -vurPhlt $CHR_INTERVALS_DIRNAME/$CHR_INTERVALS_PATTERN* $PARAMS/

    GNOMAD_GENOMES_DIRNAME=$(dirname "${GNOMAD_GENOMES_FILENAME}")
    GNOMAD_GENOMES_BASENAME=$(basename "${GNOMAD_GENOMES_FILENAME}")
    GNOMAD_GENOMES="${PARAMS}/${GNOMAD_GENOMES_BASENAME}"

    rsync -vurPhlt "$GNOMAD_GENOMES_FILENAME" $GNOMAD_GENOMES
    GNOMAD_GENOMES_PATTERN=${GNOMAD_GENOMES_BASENAME%.*}
    rsync -vurPhlt $GNOMAD_GENOMES_DIRNAME/$GNOMAD_GENOMES_PATTERN* $PARAMS/

    echo "GNOMAD_GENOMES: $GNOMAD_GENOMES"


else
    INPUT="${PARENT_DIRECTORY}"
    OUTPUTS="${OUTPUT_DIRECTORY}"
    PARAMS="${PARENT_DIRECTORY}/References"

    MUTECT_INPUT_BASENAME=$(basename "${MUTECT_INPUT_FILENAME}")
    MUTECT_INPUT=$INPUT/$MUTECT_INPUT_BASENAME
    mkdir -p $MUTECT_INPUT
    mv $MUTECT_INPUT_FILENAME $MUTECT_INPUT/

    BWA_GREF=$BWA_GREF_FILENAME
    SEQUENCE_DICT=$SEQUENCE_DICT_FILENAME
    INTERVALS_FILE=$INTERVALS_FILENAME
    CHR_INTERVALS=$CHR_INTERVALS_FILENAME
    GNOMAD_GENOMES=$GNOMAD_GENOMES_FILENAME
fi

##################################################################################################################################
#################################################---STEP 2: MUTECT.sh---########################################################## 
##################################################################################################################################       

#Note: will need changes to also work with Docker
if [[ $CONTAINER_ENGINE == "singularity" ]]; then
    singularity instance start -B $(readlink -f $TMPDIR) docker://broadinstitute/gatk:latest gatk_container
fi

if [ $GET_MUTECT = true ]; then
    echo "Mutect analysis requested"
    echo "$PARAMETER_FILE"

    if [[ $MODE != "slurm" ]]; then
         TRANSCRIPT_LIST=$TRANSCRIPT_LIST_FILENAME

         if [ $RUN_FUNCOTATOR = true ]; then
            FUNCOTATOR_SOURCES=$FUNCOTATOR_SOURCES_FILENAME
         else
            FUNCOTATOR_SOURCES=false
         fi
        
         if [ $NORMAL_SAMPLE_FILENAME != false ]; then
            NORMAL_SAMPLE=$NORMAL_SAMPLE_FILENAME
         else
            NORMAL_SAMPLE=false
        fi

        if [ $NORMAL_PILEUPS_FILENAME != false ]; then
            NORMAL_PILEUPS=$NORMAL_PILEUPS_FILENAME
        else
            NORMAL_PILEUPS=false
        fi
    fi

    if [[ $MODE == "slurm" ]]; then
        
        TRANSCRIPT_LIST_BASENAME=$(basename "${TRANSCRIPT_LIST_FILENAME}")
        TRANSCRIPT_LIST="${PARAMS}/${TRANSCRIPT_LIST_BASENAME}"

        rsync -vurPhlt "$TRANSCRIPT_LIST_FILENAME" $TRANSCRIPT_LIST


        NORMAL="${temp_dir}/Normal"
        mkdir -p $NORMAL
    
        if [[ $RUN_FUNCOTATOR = true ]] && [[ $RUN_MUTECT == true ]] ; then
            #FUNCTATOR_SOURCES_BASENAME=$(basename $FUNCOTATOR_SOURCES_FILENAME)
            #echo "^^^^^^^FUNCOTATOR_SOURCES_BASENAME: ${FUNCOTATOR_SOURCES_BASENAME}"
            echo "FUNCOTATOR_SOURCES_FILENAME: ${FUNCOTATOR_SOURCES_FILENAME} ^^^^^^^"
            FUNCOTATOR_SOURCES="$PARAMS/FUNCOTATOR_SOURCES"
            mkdir -p $FUNCOTATOR_SOURCES 
            rsync -vurPhlt $FUNCOTATOR_SOURCES_FILENAME/ $FUNCOTATOR_SOURCES
            echo "inside FUNC SOURCES"
            ls $FUNCOTATOR_SOURCES
        else
            FUNCOTATOR_SOURCES=false
        fi
    
        if [ $NORMAL_SAMPLE_FILENAME != false ]; then 
            NORMAL_SAMPLE_DIRNAME=$(dirname "${NORMAL_SAMPLE_FILENAME}")
            NORMAL_SAMPLE_BASENAME=$(basename "${NORMAL_SAMPLE_FILENAME}")
            NORMAL_SAMPLE="${NORMAL}/${NORMAL_SAMPLE_BASENAME}"

            rsync -vurPhlt "$NORMAL_SAMPLE_FILENAME" $NORMAL_SAMPLE
            NORMAL_SAMPLE_PATTERN=${NORMAL_SAMPLE_BASENAME%.*}
            echo "NORMAL_SAMPLE_PATTERN: $NORMAL_SAMPLE_PATTERN"
            rsync -vurPhlt $NORMAL_SAMPLE_DIRNAME/$NORMAL_SAMPLE_PATTERN* $NORMAL
        else
            NORMAL_SAMPLE=false
        fi
        
        if [ $NORMAL_PILEUPS_FILENAME != false ]; then
            NORMAL_PILEUPS_DIRNAME=$(dirname "${NORMAL_PILEUPS_FILENAME}")
            NORMAL_PILEUPS_BASENAME=$(basename "${NORMAL_PILEUPS_FILENAME}")
            NORMAL_PILEUPS="${NORMAL}/${NORMAL_PILEUPS_BASENAME}"

            rsync -vurPhlt "$NORMAL_PILEUPS_FILENAME" $NORMAL_PILEUPS
            NORMAL_PILEUPS_PATTERN=${NORMAL_PILEUPS_BASENAME%.*}
            rsync -vurPhlt /$NORMAL_PILEUPS_DIRNAME/$NORMAL_PILEUPS_PATTERN* $NORMAL/
        else
             NORMAL_PILEUPS=false
        fi
         
    fi
    
    
    if [ $SPLIT_BY_CHR = true ]; then
        num_intervals=$(grep -v "@" < $CHR_INTERVALS | wc -l)
        seq 1 ${num_intervals} | parallel -j8 --progress --ungroup "${CODE_DIRECTORY}/mutect_and_pileups.sh  $NORMAL_SAMPLE $INTERVALS_FILE $MUTECT_INPUT $SAMPLE_NAME $BWA_GREF $PARAMETER_FILE  $OUTPUTS $BAM_OUT $OUTPUT_DIRECTORY $MODE $SPLIT_BY_CHR $LINE_NUMBER $CHR_INTERVALS $GNOMAD_GENOMES $RUN_MUTECT ${gatk_command} {}"
        
        if ( [[ -d "${OUTPUTS}/pileups" ]] || [[ $MODE == "slurm" ]] ); then
              rsync -vurPhlt "${OUTPUTS}/pileups" "${OUTPUT_DIRECTORY}"
        fi
        # Add GatherPileupSummaries for pileups
        chr_pileups=$(find $OUTPUTS/pileups -type f | sort -V)
        #if [[ $(echo "${chr_pileups}" | wc -l) -lt ${num_intervals} ]] && [[ -d "${OUTPUT_DIRECTORY}/pileups" ]]; then
           
        if ( [[ -d "${OUTPUTS}/pileups" ]] || [[ $MODE == "slurm" ]] ); then
            rsync -vurPhlt "${OUTPUTS}/pileups" "${OUTPUT_DIRECTORY}"
        fi
    
        pileup_tables=$(find $OUTPUTS/pileups -type f | sed -e 's/^/-I /g' | tr '\n' ' ')
        echo "gatk GatherPileupSummaries --sequence-dictionary $SEQUENCE_DICT ${pileup_tables} -O ${OUTPUTS}/${SAMPLE_NAME}_pileups.table"
    
        "${gatk_command}" GatherPileupSummaries \
            --sequence-dictionary $SEQUENCE_DICT ${pileup_tables} -O ${OUTPUTS}/${SAMPLE_NAME}_pileups.table
    
    else
        ${CODE_DIRECTORY}/mutect_and_pileups.sh  $NORMAL_SAMPLE $INTERVALS_FILE $MUTECT_INPUT $SAMPLE_NAME $BWA_GREF  $PARAMETER_FILE  $OUTPUTS $BAM_OUT $OUTPUT_DIRECTORY $MODE $SPLIT_BY_CHR $LINE_NUMBER $CHR_INTERVALS $GNOMAD_GENOMES $RUN_MUTECT "${gatk_command}"    
    fi
    
    if [[ $RUN_MUTECT == false ]]; then
            echo "THIS IS THE RSYNC STEP"
            OVERALL_OUTPUT_FOLDER=$(dirname "${OUTPUT_DIRECTORY}")
            NORMAL_PILEUPS_FOLDER="${OVERALL_OUTPUT_FOLDER}/NORMAL_PILEUPS"
            mkdir -p $NORMAL_PILEUPS_FOLDER
            rsync -vurPhlt $OUTPUTS/ $NORMAL_PILEUPS_FOLDER
    else
            echo "FOR SOME REASON THIS RSYNC STEP WAS SKIPPED"
    fi
    
   
   echo "mutect_and_pileups.sh analysis complete"
   
   if [[ $RUN_MUTECT == true ]]; then 
  
       if [ $SPLIT_BY_CHR = true ]; then
           chr_vcfs=$(find $OUTPUTS/vcfs -type f | sort -V)
    
           #if [[ $(echo "${chr_vcfs}" | wc -l) -lt ${num_intervals} ]] && [[ -d "${OUTPUT_DIRECTORY}/vcfs" ]]; then
           if ( [[ -d "${OUTPUTS}/vcfs" ]] || [[ $MODE == "slurm" ]] ); then
              rsync -vurPhlt "${OUTPUTS}/vcfs" "${OUTPUT_DIRECTORY}"
           fi
    
           # Change to MergeVCFs
           echo "bcftools concat $(find $OUTPUTS/vcfs -type f | grep -E ".*.vcf$" | sort -V) > ${OUTPUTS}/${SAMPLE_NAME}_mutect2.vcf"
           vcf_files=$(find $OUTPUTS/vcfs -type f | grep -E ".*.vcf$" | sed -e 's/^/--INPUT /g' | tr '\n' ' ')
           "${gatk_command}" MergeVCFs \
              ${vcf_files} --OUTPUT ${OUTPUTS}/${SAMPLE_NAME}_mutect2.vcf
           #bcftools concat $(find $OUTPUTS/vcfs -type f | grep -E ".*.vcf$" | sort -V) > ${OUTPUTS}/${SAMPLE_NAME}_mutect2.vcf
    
           vcf_stats=$(find $OUTPUTS/vcfs -type f | grep -E ".*.vcf.stats$" | sed -e 's/^/--stats /g' | tr '\n' ' ')
          "${gatk_command}" MergeMutectStats \
              ${vcf_stats} -O "${OUTPUTS}/${SAMPLE_NAME}_mutect2.vcf.stats"
    
           echo "gatk MergeMutectStats ${vcf_stats} -O ${OUTPUTS}/${SAMPLE_NAME}_mutect2.vcf.stats"
           
           f1r2_files=$(find $OUTPUTS/f1r2 -type f | sed -e 's/^/-I /g' | tr '\n' ' ')
           "${gatk_command}" LearnReadOrientationModel \
              ${f1r2_files} -O "${OUTPUTS}/${SAMPLE_NAME}_mutect2_artifact_prior.tar.gz"
      else
        echo "running LearnReadOrientationModel"
        echo "LearnReadOrientationModel input: ${OUTPUTS}/f1r2/${SAMPLE_NAME}_f1r2.tar.gz"
        ${gatk_command} LearnReadOrientationModel \
            -I ${OUTPUTS}/f1r2/${SAMPLE_NAME}_f1r2.tar.gz -O "${OUTPUTS}/${SAMPLE_NAME}_mutect2_artifact_prior.tar.gz"
       fi
    
       rsync -vurPhlt "${OUTPUTS}/f1r2" "${OUTPUT_DIRECTORY}"
       echo "what is in f1r2:" 
       ls "$OUTPUT_DIRECTORY/f1r2"
      
    
      echo "$OUTPUTS folder contains:"
      tree -h $OUTPUTS
    
        ${CODE_DIRECTORY}/mutect_filter.sh $SAMPLE_NAME $BWA_GREF $PARAMETER_FILE $OUTPUTS $OUTPUT_DIRECTORY $NORMAL_SAMPLE $NORMAL_PILEUPS $MODE $LINE_NUMBER "${gatk_command}"
    
        #header for vcf header need to change
        #_mutect2_filter.vcf
    
        ${CODE_DIRECTORY}/funcotator.sh $SAMPLE_NAME $BWA_GREF $FUNCOTATOR_SOURCES $TRANSCRIPT_LIST $PARAMETER_FILE $FILTERED $OUTPUTS $RUN_FUNCOTATOR $OUTPUT_DIRECTORY $MODE $LINE_NUMBER "${gatk_command}"
        
       if [[ $MODE == "slurm" ]]; then 
            rsync -vurPhlt $OUTPUTS/ $OUTPUT_DIRECTORY
            echo "Results copied over to specified output directory"
       fi
    fi
else
    echo "No mutect analysis requested"
fi
            
##################################################################################################################################
############################################---STEP 3: HAPLOTYPECALLER.sh---######################################################
##################################################################################################################################    
if [[ $GET_HAPLOTYPE = true ]] && [[ $RUN_MUTECT = true ]]; then
    echo "Haplotypecaller analysis requested"
    ${CODE_DIRECTORY}/haplotypecaller.sh $SAMPLE_NAME $BWA_GREF $TWIST_SNPS $PARAMETER_FILE $MODE
    echo "Haplotypecaller analysis complete"

    if [[ $MODE == "slurm" ]]; then 
        rsync -vurPhlt $OUTPUTS/ $OUTPUT_DIRECTORY
        echo "Results copied over to specified output directory"
   fi

else
    echo "No HaplotypeCaller analysis requested"
fi

singularity instance stop gatk_container
singularity delete --force gatk_container

##################################################################################################################################
################################################---STEP 4: VARSCAN.sh---########################################################## 
##################################################################################################################################    
if [[ $GET_VARSCAN = true ]] && [[ $RUN_MUTECT = true ]]; then
    echo "Varscan analysis requested"
   ${CODE_DIRECTORY}/varscan.sh $SAMPLE_NAME $BWA_GREF $MIN_COVERAGE $MIN_VAR_FREQ $P_VALUE $PARAMETER_FILE $MODE
   echo "Varscan analysis complete"

    if [[ $MODE == "slurm" ]]; then 
        rsync -vurPhlt $OUTPUTS/ $OUTPUT_DIRECTORY
        echo "Results copied over to specified output directory"
   fi

else 
    echo "No Varscan analysis requested"
fi
