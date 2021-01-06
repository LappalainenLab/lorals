#!/bin/bash

# Align and process ONT reads with minimap2 in splice mode
# The alignment will happen to each copy of your genome
# Will align, sort and indec

set -eo pipefail

declare -a DEPDENCIES=(minimap2 samtools awk)
for DEP in ${DEPDENCIES[@]}; do $(command -v ${DEP} >/dev/null 2>&1) || (echo "Cannot find ${DEP}" >&2; exit 1); done

OUTDIR_DEFAULT="$(pwd -P)/aligned"
THREADS_DEFAULT=8
MIN_MAPQ_DEFAULT=10
INTERMEDIATE_DEFAULT="no"

function Usage() {
    echo -e "\
Usage: $(basename $0) -f \e[3minput.fq\e[0m -T \e[3mref.fa\e[0m [-r \e[3mrev.fq\e[0m] [-o \e[3m/path/to/outdir\
e[0m] [--csi] \n\
Where:  -f|--fastq is the path to the input or forward FASTQ file \n\
        -G|--reference is the path to the reference FASTA file \n\
        [-mapq|--min_mapq] is the minimum mapq score to keep an aligned read \n\
        [-intermediate|--keep_intermediate] is an option of whether to keep the intermediate files \n\
        [-r|--reverse] is an optional reverse FASTQ \n\
        [-t|--threads] is the number of threads to use \n\
        [-o|--outdir] is an optional path to an output directory \n\
            defaults to \e[1m${OUTDIR_DEFAULT}\e[0m \n\
        [--csi] will index the final BAM file with a CSI index \n\
            defaults to a BAI index
" >&2
    exit 1
}

function extension() {
    local fname="$1"
    local ext="$(echo ${fname} | rev | cut -f 1 -d '.' | rev)"
    $(grep -E 'gz|bz|zip' <(echo "${ext}") > /dev/null 2> /dev/null) && ext="$(echo ${fname} | rev | cut -f -2 -d '.' | rev)"
    echo ".${ext}"
}

[[ "$#" -lt 1 ]] && Usage
while [[ "$#" -ge 1 ]]; do
    case "$1" in
        -f|--fastq)
            FASTQ="$2"
            shift
            ;;
        -r|--reverse)
            REVERSE="$2"
            shift
            ;;
        -G|--reference)
            REFERENCE="$2"
            shift
            ;;
        -o|--outdir)
            OUTDIR="$2"
            shift
            ;;
        -mapq|--min_mapq)
            MIN_MAPQ="$2"
            shift
            ;;
        -intermediate|--keep_intermediate)
            INTERMEDIATE="$2"
            shift
            ;;
        --csi)
            INDEX='-c'
            ;;
        --threads)
            THREADS="$2"
            ;;
        *)
            Usage
            ;;
    esac
    shift
done

[[ -z "${FASTQ}" || -z "${REFERENCE}" ]] && Usage
[[ -f "${FASTQ}" ]] || (echo "Cannot find input FASTQ ${FASTQ}" >&2; exit 1)
#[[ -f "${REFERENCE}" ]] || (echo "Cannot find reference genome ${REFERENCE}" >&2; exit 1)

[[ -z "${INDEX}" ]] && INDEX='-b'
[[ -z "${REVERSE}" ]] && REVERSE='' || [[ -f "${REVERSE}" ]] || (echo "Cannot find reverse FASTQ ${REVERSE}" >&2; exit 1)
[[ -z "${OUTDIR}" ]] && OUTDIR="${OUTDIR_DEFAULT}"
[[ -z "${THREADS}" ]] && THREADS="${THREADS_DEFAULT}"
[[ -z "${MIN_MAPQ}" ]] && MIN_MAPQ="${MIN_MAPQ_DEFAULT}"
[[ -z "${INTERMEDIATE}" ]] && INTERMEDIATE="${INTERMEDIATE_DEFAULT}"

EXTENSION="$(extension ${FASTQ})"
NAME="$(basename ${FASTQ} $(extension ${FASTQ}))"

(set -x; mkdir ${NAME}_temp)
cd ${NAME}_temp
echo "Processing sample" ${NAME} >> ${NAME}_hap_aware.log

# Align to each haplotype using minimap2 and get sorted bam files
# Only keep primary alignments with MAPQ more than 10
(set -x; minimap2 -t ${THREADS} -ax splice -uf -k14 ${REFERENCE}_hap1.fa ${FASTQ} ${REVERSE} | \
     samtools view -q ${MIN_MAPQ} -F 2304 -Sb | \
     samtools sort -@ ${THREADS} - -o ${NAME}_reads_aln_sorted.hap1.bam)
(set -x; minimap2 -t ${THREADS} -ax splice -uf -k14 ${REFERENCE}_hap2.fa ${FASTQ} ${REVERSE} | \
     samtools view -q ${MIN_MAPQ} -F 2304 -Sb | \
     samtools sort -@ ${THREADS} - -o ${NAME}_reads_aln_sorted.hap2.bam)

(set -x; samtools index ${NAME}_reads_aln_sorted.hap1.bam)
(set -x; samtools index ${NAME}_reads_aln_sorted.hap2.bam)

# Select reads based on MAPQ score
(set -x; samtools view ${NAME}_reads_aln_sorted.hap1.bam | \
     cut -f1,5 | sort -k1,1 > ${NAME}_mapq_hap1.txt)
(set -x; samtools view ${NAME}_reads_aln_sorted.hap2.bam | \
     cut -f1,5 | sort -k1,1 > ${NAME}_mapq_hap2.txt)

(set -x; join -1 1 -2 1 ${NAME}_mapq_hap1.txt ${NAME}_mapq_hap2.txt -e 0 > ${NAME}_mapq_hap.txt)

(set -x; awk '$2 > $3' ${NAME}_mapq_hap.txt | cut -f1 -d' ' > ${NAME}_reads_from_hap1.txt)
echo "Reads from hap1" $(wc -l ${NAME}_reads_from_hap1.txt | cut -f1 -d' ') >> ${NAME}_hap_aware.log
(set -x; awk '$2 < $3' ${NAME}_mapq_hap.txt | cut -f1 -d' ' > ${NAME}_reads_from_hap2.txt)
echo "Reads from hap2" $(wc -l ${NAME}_reads_from_hap2.txt | cut -f1 -d' ') >> ${NAME}_hap_aware.log
(set -x; awk '$2 == $3' ${NAME}_mapq_hap.txt | cut -f1 -d' ' > ${NAME}_reads_equal.txt)
echo "Reads equal" $(wc -l ${NAME}_reads_equal.txt | cut -f1 -d' ') >> ${NAME}_hap_aware.log

(set -x; sort -R ${NAME}_reads_equal.txt > ${NAME}_temp.txt)
(set -x; head -n $(($(wc -l ${NAME}_reads_equal.txt | cut -f1 -d' ') / 2)) ${NAME}_reads_equal.txt >> ${NAME}_reads_from_hap1.txt)
(set -x; tail -n +$(($(wc -l ${NAME}_reads_equal.txt | cut -f1 -d' ') / 2 + 1)) ${NAME}_reads_equal.txt >> ${NAME}_reads_from_hap2.txt)

#Compile new bam file
(set -x; samtools view ${NAME}_reads_aln_sorted.hap1.bam | grep -wFf ${NAME}_reads_from_hap1.txt - > ${NAME}_reads_aln_sorted.hap1_select.sam)
(set -x; samtools view ${NAME}_reads_aln_sorted.hap2.bam | grep -wFf ${NAME}_reads_from_hap2.txt - > ${NAME}_reads_aln_sorted.hap2_select.sam)

(set -x; samtools view -H ${NAME}_reads_aln_sorted.hap1.bam > ${NAME}_header.txt)

cd ..

(set -x; cat ${NAME}_temp/${NAME}_header.txt ${NAME}_temp/${NAME}_reads_aln_sorted.hap1_select.sam ${NAME}_temp/${NAME}_reads_aln_sorted.hap2_select.sam | \
     samtools view -Sb | \
     samtools sort -@ ${THREADS} - -o ${NAME}_reads_aln_sorted.merged.bam)
(set -x; samtools index ${NAME}_reads_aln_sorted.merged.bam)


if ${INTERMEDIATE}=="no"; then
     echo "Deleting intermediate files"
     (set -x; rm -rf ${NAME}_temp)
fi
