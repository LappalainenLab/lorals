#!/bin/bash

# Using an aligned bam file it will correct phased haplotypes in vcf file and output a new vcf
set -eo pipefail

declare -a DEPENDENCIES=(parallel bcftools bgzip tabix hapcut2 gatk)
for DEP in ${DEPENDENCIES[@]}; do $(command -v ${DEP} > /dev/null 2> /dev/null) || (echo "Cannot find ${DEP}" >&2; exit 1); done

OUTDIR_DEFAULT="$(pwd -P)/vcf"
THREADS_DEFAULT=8

function Usage() {
    echo -e "\
Usage: $(basename $0) -b \e[3minput.bam\e[0m -G \e[3mreference.fasta\e[0m -V \e[3mVCF.vcf\e[0m \n\
Where:  -b|--bam is the path to the input bam file \n\
        -G|--reference is the path to the reference FASTA file \n\
        -V|--vcf is the existing vcf file \n\
        [-t|--threads] is the number of threads to use \n\
        [-o|--outdir] is an optional path to an output directory \n\
            defaults to \e[1m${OUTDIR_DEFAULT}\e[0m \n\
" >&2
    exit 1
}


function extension() {
    local fname="$1"
    local ext="$(echo ${fname} | rev | cut -f 1 -d '.' | rev)"
    $(grep -E 'gz|bz|zip' <(echo "${ext}") > /dev/null 2> /dev/null) && ext="$(echo ${fname} | rev | cut -f -2 -d '.' | r
ev)"
    echo ".${ext}"
}

[[ "$#" -lt 1 ]] && Usage

while [[ "$#" -ge 1 ]]; do
    case "$1" in
        -b|--bam)
            BAM="$2"
            shift
            ;;
        -G|--reference)
            REFERENCE="$2"
            shift
            ;;
        -V|--vcf)
            VCF="$2"
            shift
            ;;
        -o|--outdir)
            OUTDIR="$2"
            shift
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

[[ -z "${BAM}" || -z "${REFERENCE}" || -z "${VCF}" ]] && Usage
[[ -f "${BAM}" ]] || (echo "Cannot find input FASTQ ${FASTQ}" >&2; exit 1)
[[ -f "${REFERENCE}" ]] || (echo "Cannot find reference genome ${REFERENCE}" >&2; exit 1)
[[ -f "${VCF}" ]] || (echo "Cannot find VCF ${VCF}" >&2; exit 1)

[[ -z "${OUTDIR}" ]] && OUTDIR="${OUTDIR_DEFAULT}"
[[ -z "${THREADS}" ]] && THREADS="${THREADS_DEFAULT}"

EXTENSION="$(extension ${BAM})"
NAME="$(basename ${BAM} $(extension ${BAM}))"

# Use hapcut2 to extract reads containing het variants
(set -x; extractHAIRS --nf 1 --ont 1 --bam ${BAM} --VCF ${VCF} --out ${NAME}.extractHairs.txt --ref ${REFERENCE})

(set -x; HAPCUT2 --fragments ${NAME}.extractHairs.txt --VCF ${VCF} --output ${NAME}.hapcut.txt --nf 1)

(set -x; python HapCUT22VCF.py ${NAME}.hapcut.txt ${NAME}.vcf ${VCF})

(set -x; bgzip -c ${NAME}.vcf > ${NAME}.vcf.gz)
(set -x; tabix -p vcf ${NAME}.vcf.gz)

(set -x; bcftools consensus -f ${REFERENCE} -H 1 -o ${NAME}.hap1.fa ${NAME}.vcf.gz)
(set -x; bcftools consensus -f ${REFERENCE} -H 2 -o ${NAME}.hap2.fa ${NAME}.vcf.gz)

(set -x; gatk CreateSequenceDictionary -R ${NAME}.hap1.fa -O ${NAME}.hap1.dict)
(set -x; gatk CreateSequenceDictionary -R ${NAME}.hap2.fa -O ${NAME}.hap2.dict)
