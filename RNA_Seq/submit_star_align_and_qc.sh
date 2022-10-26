if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    echo "Format: ./submit_star_align_and_qc.sh [fastq_directory] [output_directory] --argument"
    echo "fastq_directory: path to raw .fastq or .fastq.gz files"
    echo "output_directory: path for BWA and mutect output"
    echo "argument: indicates which reference genome you would like your data aligned to (--human or --mouse)"
    exit 1
else
    TEMP=`getopt -o vdm: --long human,mouse,hsapiens,sapiens,mmusculus,musculus \
    -n 'submit_BWA_CHIP.sh' -- "$@"`

    if [ $? != 0 ] ; then echo "Unrecognized argument. Possible arguments: human or mouse. The argument represents the species of the reference genome you would like aligned to your FASTQ files." >&2 ; exit 1 ; fi
        eval set -- "$TEMP"

        get_human=false
        get_mouse=false
        
    while true; do
        case "$1" in
            --human | --hsapiens | --sapiens ) get_human=true; shift ;;
            --mouse | --mmusculus | --musculus ) get_mouse=true; shift ;;
            -- ) shift; break ;;
            * ) break ;;
        esac
    done

    #prevent the FastQ's from being aligned to both human and mouse genomes
    if [ $get_human = true ] && [ $get_mouse = true ]; then
        echo "Please pick one species to align the genome to (human OR mouse)"; exit 1
    fi
    
    #specify location of star_align_and_qc.sh
    code_directory=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

    fastq_directory=$1 #get directory path from second argument (first argument $0 is the path of this script)
    output_directory=$2
    fastq_file="${fastq_directory}/fastq_files" #give a path to a file to store the paths to the fastq files in $fastq_directory
    find "${fastq_directory}/" -type f | grep -E "fastq|fq" | grep -v "Undetermined" | sed -e 's/_R1.*$//g' | sed -e 's/_R2.*$//g' | sort -u > "${fastq_file}" #generate list of full paths to fastq files and save to the file in $fastq_list
    array_length=$(wc -l < "${fastq_file}") #get the number of files 

    mkdir -p $output_directory
    mkdir -p "${output_directory}/Logs"

    sbatch -o "${output_directory}/Logs/%A_%a.log" `#put standard output into log` \
        -a "1-${array_length}" `#initiate job array equal to the number of fastq files` \
        -W `#indicates to the script not to move on until the sbatch operation is complete` \
        "${code_directory}/star_align_and_qc.sh" `#specify STAR alignment script` \
        "${output_directory}" `#provide absolute path to fastq file` \
        "${code_directory}" \
        "${fastq_file}" \
        "${get_human}" \
        "${get_mouse}"

    wait
    
    if [ get_human = true ]; then
        export animal="human"; else
            export animal="mouse"
    fi

    module load R/4.0
    Rscript aggregate_counts.R $output_directory $animal \
        "$output_directory" > "$output_directory/Logs/OutFile.Rout" 2>&1
        
fi

