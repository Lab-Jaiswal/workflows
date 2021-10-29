#!/bin/bash
#Unrecognized argument. Possible arguments: mutect, varscan, haplotypecaller, all, min_coverage, min_var_freq, and p_value.
output_directory=$1
code_directory=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

chmod a+x CellRanger/submit_cellranger.sh
chmod a+x BWA_CHIP/submit_BWA_CHIP.sh
chmod a+x RNA_Seq/submit_star_align_and_qc.sh

mkdir -p "$output_directory/BWA_CHIP"
mkdir -p "$output_directory/RNA_Seq"
mkdir -p "$output_directory/CellRanger"
mkdir -p "$output_directory/BWA_CHIP/Mutect"
mkdir -p "$output_directory/BWA_CHIP/Varscan"
mkdir -p "$output_directory/BWA_CHIP/Haplotype"
mkdir -p "$output_directory/BWA_CHIP/All"
mkdir -p "$output_directory/BWA_CHIP/Varscan_Edited_Fields"
mkdir -p "$output_directory/RNA_Seq/Human"
mkdir -p "$output_directory/RNA_Seq/Mouse"
mkdir -p "$output_directory/CellRanger/Human"
mkdir -p "$output_directory/CellRanger/Mouse"
mkdir -p "$output_directory/CellRanger/Human_Nuclei"

rm test_case_summary.txt

./BWA_CHIP/submit_BWA_CHIP.sh /home/maurertm/labs/maurertm/mutectvarscanhaplotype $output_directory/BWA_CHIP --gkbhunjki >/dev/null 2>&1
Unrec_BWA=$?
./RNA_Seq/submit_star_align_and_qc.sh /labs/sjaiswal/maurertm/output/FASTQ $output_directory/RNA_Seq --htsryjth >/dev/null 2>&1
Unrec_RNASeq=$?
./CellRanger/submit_cellranger.sh /oak/stanford/groups/sjaiswal/jk/JG136/day5_10X_deepseq/raw_data $output_directory/CellRanger --tshegrth >/dev/null 2>&1
Unrec_CellRanger=$?

if [[ $Unrec_BWA -eq 1 ]]; then
    echo "Unrecognized BWA argument test case: pass" >> test_case_summary.txt; else 
    echo "Unrecognized BWA argument test case: fail" >> test_case_summary.txt
fi

if [[ $Unrec_RNASeq -eq 1 ]]; then
    echo "Unrecognized RNA_Seq argument test case: pass" >> test_case_summary.txt; else
    echo "Unrecognized RNA_Seq argument test case: fail" >> test_case_summary.txt
fi

if [[ $Unrec_CellRanger -eq 1 ]]; then
    echo "Unrecognized CellRanger argument test case: pass" >> test_case_summary.txt; else
    echo "Unrecognized CellRanger argument test case: fail" >> test_case_summary.txt
fi

./BWA_CHIP/submit_BWA_CHIP.sh /home/maurertm/labs/maurertm/mutectvarscanhaplotype $output_directory/BWA_CHIP >/dev/null 2>&1
NoArg_BWA=$?
./RNA_Seq/submit_star_align_and_qc.sh /labs/sjaiswal/maurertm/output/FASTQ $output_directory/RNA_Seq >/dev/null 2>&1
NoArg_RNASeq=$?
./CellRanger/submit_cellranger.sh /oak/stanford/groups/sjaiswal/jk/JG136/day5_10X_deepseq/raw_data $output_directory/CellRanger >/dev/null 2>&1
NoArg_CellRanger=$?

if [[ $NoArg_BWA -eq 1 ]]; then
    echo "Missing BWA argument test case: pass" >> test_case_summary.txt; else
    echo "Missing BWA argument test case: fail" >> test_case_summary.txt
fi

if [[ $NoArg_RNASeq -eq 1 ]]; then
    echo "Missing RNA_Seq argument test case: pass" >> test_case_summary.txt; else
    echo "Missing RNA_Seq argument test case: fail" >> test_case_summary.txt
fi

if [[ $NoArg_CellRanger -eq 1 ]] ; then
    echo "Missing CellRanger argument test case: pass" >> test_case_summary.txt; else
    echo "Missing CellRanger argument test case: fail" >> test_case_summary.txt
fi

./BWA_CHIP/submit_BWA_CHIP.sh /home/maurertm/labs/maurertm/mutectvarscanhaplotype $output_directory/BWA_CHIP --mutect --p_value=2 >/dev/null 2>&1
Mutect_Pval=$?
./BWA_CHIP/submit_BWA_CHIP.sh /home/maurertm/labs/maurertm/mutectvarscanhaplotype $output_directory/BWA_CHIP --haplotypecaller --p_value=2 >/dev/null 2>&1
Haplotype_Pval=$?
./BWA_CHIP/submit_BWA_CHIP.sh /home/maurertm/labs/maurertm/mutectvarscanhaplotype $output_directory/BWA_CHIP --mutect --min_coverage=2 >/dev/null 2>&1
Mutect_Mincov=$?
./BWA_CHIP/submit_BWA_CHIP.sh /home/maurertm/labs/maurertm/mutectvarscanhaplotype $output_directory/BWA_CHIP --haplotypecaller --min_coverage=2 >/dev/null 2>&1
Haplotype_Mincov=$?
./BWA_CHIP/submit_BWA_CHIP.sh /home/maurertm/labs/maurertm/mutectvarscanhaplotye $output_directory/BWA_CHIP --mutect --min_var_freq=2 >/dev/null 2>&1
Mutect_MVF=$?
./BWA_CHIP/submit_BWA_CHIP.sh /home/maurertm/labs/maurertm/mutectvarscanhaplotype $output_directory/BWA_CHIP --haplotypecaller --min_var_freq=2 >/dev/null 2>&1
Haplotype_MVF=$?

if ( [[ $Mutect_Pval -eq 1 ]] && \
    [[ $Haplotype_Pval -eq 1 ]] && \
    [[ $Mutect_Mincov -eq 1 ]] && \
    [[ $Haplotype_Mincov -eq 1 ]] && \
    [[ $Mutect_MVF -eq 1 ]] && \
    [[ $Haplotype_MVF -eq 1 ]] ); then
    echo "BWA_CHIP test case where user selects --mutect or --haplotypecaller and a varscan only argument (p_value, min_coverage, or min_var_freq): pass" >> test_case_summary.txt; else
    echo "BWA_CHIP test case where user selects --mutect or --haplotypecaller and a varscan only argument (p_value, min_coverage, or min_var_freq): fail" >> test_case_summary.txt
fi

./RNA_Seq/submit_star_align_and_qc.sh /labs/sjaiswal/maurertm/output/FASTQ $output_directory/RNA_Seq --human --mouse >/dev/null 2>&1
HumanMouse_RNASeq=$?
./CellRanger/submit_cellranger.sh /oak/stanford/groups/sjaiswal/jk/JG136/day5_10X_deepseq/raw_data $output_directory/CellRanger --human --mouse >/dev/null 2>&1
HumanMouse_CellRanger=$?
./CellRanger/submit_cellranger.sh /oak/stanford/groups/sjaiswal/jk/JG136/day5_10X_deepseq/raw_data $output_directory/CellRanger --human --human_nuclei >/dev/null 2>&1
HumanNuclei_CellRanger=$?
./CellRanger/submit_cellranger.sh /oak/stanford/groups/sjaiswal/jk/JG136/day5_10X_deepseq/raw_data $output_directory/CellRanger --mouse --human_nuclei >/dev/null 2>&1
MouseNuclei_CellRanger=$?
./CellRanger/submit_cellranger.sh /oak/stanford/groups/sjaiswal/jk/JG136/day5_10X_deepseq/raw_data $output_directory/CellRanger --mouse --human_nuclei --human >/dev/null 2>&1
All_CellRanger=$?

if [[ $HumanMouse_RNASeq -eq 1 ]]; then
    echo "RNA_Seq test case where user selects both Mouse AND Human genomes: pass" >> test_case_summary.txt; else
    echo "RNA_Seq test case where user selects both Mouse AND Human genomes: fail" >> test_case_summary.txt
fi

if ( [[ $HumanMouse_CellRanger -eq 1 ]] && \
    [[ $HumanNuclei_CellRanger -eq 1 ]] && \
    [[ $MouseNuclei_CellRanger -eq 1 ]] && \
    [[ $All_CellRanger -eq 1 ]] ); then
    echo "CellRanger test case where user selects multiple genome types: pass" >> test_case_summary.txt; else
    echo "CellRanger test case where user selects multiple genome types: fail" >> test_case_summary.txt
fi

cd "$code_directory/BWA_CHIP"
./submit_BWA_CHIP.sh /labs/sjaiswal/Herra/199011828_R6/FASTQ/ $output_directory/BWA_CHIP/Mutect --mutect
./submit_BWA_CHIP.sh /labs/sjaiswal/Herra/199011828_R6/FASTQ/ $output_directory/BWA_CHIP/Varscan --varscan
./submit_BWA_CHIP.sh /labs/sjaiswal/Herra/199011828_R6/FASTQ/ $output_directory/BWA_CHIP/Haplotype --haplotypecaller
./submit_BWA_CHIP.sh /labs/sjaiswal/Herra/199011828_R6/FASTQ/ $output_directory/BWA_CHIP/All --all
./submit_BWA_CHIP.sh /labs/sjaiswal/Herra/199011828_R6/FASTQ/ $output_directory/BWA_CHIP/Varscan_Edited_Fields --varscan --p_value=0.05, --Min_var_freq=0.002, --min_coverage=20

cd "$code_directory/RNA_Seq"
./submit_star_align_and_qc.sh /labs/sjaiswal/maurertm/output/FASTQ $output_directory/RNA_Seq/Human --human
./submit_star_align_and_qc.sh /labs/sjaiswal/maurertm/output/FASTQ $output_directory/RNA_Seq/Mouse --mouse

cd "$code_directory/CellRanger"
./submit_cellranger.sh /oak/stanford/groups/sjaiswal/jk/JG136/day5_10X_deepseq/raw_data $output_directory/CellRanger/Human --human
./submit_cellranger.sh /oak/stanford/groups/sjaiswal/jk/JG136/day5_10X_deepseq/raw_data $output_directory/CellRanger/Mouse --mouse
./submit_cellranger.sh /oak/stanford/groups/sjaiswal/jk/JG136/day5_10X_deepseq/raw_data $output_directory/CellRanger/Human_Nuclei --human_nuclei

