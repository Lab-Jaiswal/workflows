#!/bin/bash
echo "entering peak_calling.sh"

output_path=$1
output_temp_dir=$2
gsize=$3
extsize=$4
shifts=$5
broad=$6
nomodel=$7
whitelist=$8
parameter_file=$9
PREFIX=${10}

echo "peak_calling.sh used the following parameters:
$0 $1 $2 $3 $4 $5 $6 $7 $8 $9"

if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
	echo "#####################peak_calling.sh file#####################
        command given to peak_calling.sh $0 $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}
        The paramters given:            
            output_path=$1
            output_temp_dir=$2
            gsize=$3
            extsize=$4
            shifts=$5
            broad=$6
            nomodel=$7
            whitelist=$8
            parameter_file=$9
            PREFIX=${10}
            " >> $parameter_file
fi

###########################---STEP 6: PEAK CALLING WITH MACS2---########################################
if [ ! -f "$output_temp_dir/peak_calling/${PREFIX}/${PREFIX}_peaks.broadPeak" ]; then

      if [ ! -f "$output_temp_dir/peak_calling" ]; then
            mkdir "$output_temp_dir/peak_calling/"
      fi

      if [ ! -f "$output_temp_dir/peak_calling/${PREFIX}" ]; then
            mkdir "$output_temp_dir/peak_calling/${PREFIX}"
      fi
      
      module load macs2
      module load samtools/1.9

      echo "Running macs2 with .bam file: $output_temp_dir/${PREFIX}.merged.sorted.bam"

      if [ $broad == true ] && [ $nomodel == true ]; then
            macs2 callpeak -t "$output_temp_dir/${PREFIX}.merged.sorted.bam" --name ${PREFIX} --outdir "$output_temp_dir/peak_calling/${PREFIX}" --gsize $gsize --nomodel --shift -$shifts --extsize $extsize --broad 
      fi

      if [ $broad == true ] && [ $nomodel == false ]; then
            macs2 callpeak -t "$output_temp_dir/${PREFIX}.merged.sorted.bam" --name ${PREFIX} --outdir "$output_temp_dir/peak_calling/${PREFIX}" --gsize $gsize --shift -$shifts --extsize $extsize --broad 
      fi

      if [ $broad == false ] && [ $nomodel == true ]; then
            macs2 callpeak -t "$output_temp_dir/${PREFIX}.merged.sorted.bam" --name ${PREFIX} --outdir "$output_temp_dir/peak_calling/${PREFIX}" --gsize $gsize --nomodel --shift -$shifts --extsize $extsize
      fi
      
      echo "peak calling with macs2 complete"
else
      echo "peaking calling with macs2 already done"
fi
      
if [ ! -f "$output_temp_dir/peak_calling/${PREFIX}/${PREFIX}_raw.bed" ]; then
      
    cp "$output_temp_dir/peak_calling/${PREFIX}/${PREFIX}_peaks.broadPeak" "$output_temp_dir/peak_calling/${PREFIX}/${PREFIX}_raw.bed"
    rsync -vur "$output_temp_dir/" "$output_path"

fi


if [ $SLURM_ARRAY_TASK_ID -eq 1 ]; then
	echo "peak_calling.sh script complete
        " >> $parameter_file
fi
