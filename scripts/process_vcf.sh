#!/bin/bash

set -eo pipefail

declare -a DEPENDENCIES=(parallel bcftools bgzip tabix gatk)
for DEP in ${DEPENDENCIES[@]}; do $(command -v ${DEP} > /dev/null 2> /dev/null) || (echo "Cannot find ${DEP}" >&2; exit 1); done

OUTDIR_DEFAULT="$(pwd -P)/sample_vcfs"

function Usage() {
    echo -e "\
1. Filter VCF to only contain biallelic variants
2. Split a VCF containing records for multiple individuals into one VCF per individual and tabix the files \n\
3. For each sample create two fasta ref files for each haplotype
4. For each sample VCF only keep het variants
\n\
Usage: $(basename $0) -v|--vcf in.vcf -f|--fasta genome.fa [-o|--outdir /path/to/outdir] [-s|--samples list_samples.txt]\n\
Where:  -v|--vcf is an input VCF file to be split into sample VCF files \n\
        -f|--fasta is the genome reference fasta file \n\
        [-o|--outdir] is an optional path to an output directory to put sample VCF files \n\
            defaults to ${OUTDIR_DEFAULT} \n\
        [-s|--samples] is an optional file containing the list of samples to keep from the vcf \n\
            defaults to ${SAMPLES_DEFAULT} \n\
" >&2
    exit 1
}

export -f Usage

function SplitVCF() {
    local vcf="$1"
    local outdir="$2"
    local donor="$3"
    local outname="${outdir}/${donor}.vcf.gz"
    (set -x; bcftools view --force-samples -s "${donor}" "${vcf}" | bgzip -c > "${outname}")
    (set -x; tabix -p vcf "${outname}")
}

export -f SplitVCF

function generateFA() {
    local fasta="$1"
    local outdir="$2"
    local donor="$3"
    local hap="$4"
    local faname="${outdir}/${donor}_hap${hap}.fa"
    local dictname="${outdir}/${donor}_hap${hap}.dict"
    local vcfname="${outdir}/${donor}.vcf.gz"
    (set -x; bcftools consensus -f ${fasta} -H ${hap} -o ${faname} ${vcfname})
    (set -x; gatk CreateSequenceDictionary -R ${faname} -O ${dictname})
}

export -f generateFA

function keepHet() {
    local donor="$1"
    local outdir="$2"
    local samplename="${outdir}/${donor}.vcf.gz"
    local outname="${outdir}/${donor}_het.vcf.gz"
    (set -x; bcftools view -m2 -M2 -v snps -Ov -g het ${samplename} | bgzip -c > ${outname})
    (set -x; tabix -p vcf ${outname})
}

export -f keepHet

[[ "$#" -lt 1 ]] && Usage

while [[ "$#" -gt 1 ]]; do
    case "$1" in
        -v|--vcf)
            VCF="$2"
            shift
            ;;
        -f|--fasta)
            GENOME=$2
            shift
            ;;
        -o|--outdir)
            OUTDIR="$2"
            shift
            ;;
        -s|--samples)
            SAMPLES="$2"
            shift
            ;;
        *)
            Usage
            ;;
    esac
    shift
done

[[ -z "${VCF}" ]] && Usage
[[ -f "${VCF}" ]] || (echo "Cannot find VCF file ${VCF}" >&2; exit 1)
[[ -z "${GENOME}" ]] && Usage
[[ -f "${GENOME}" ]] || (echo "Cannot find genome fasta reference file ${GENOME}" >&2; exit 1)
[[ -z "${SAMPLES}" ]] && Usage
[[ -f "${SAMPLES}" ]] || (echo "Cannot find samples file ${SAMPLES}" >&2; exit 1)
[[ -z "${OUTDIR}" ]] && OUTDIR="${OUTDIR_DEFAULT}"

(set -x; bcftools norm -d all -Ou -N ${VCF} | \
        bcftools view -m2 -M2 -v snps -S ${SAMPLES} -Ov --force-samples | \
        bgzip -c > ${VCF/.vcf.gz/_biallelic.vcf.gz}
)
(set -x; tabix ${VCF/.vcf.gz/_biallelic.vcf.gz})

(set -x; mkdir -vp "${OUTDIR}")

#declare -a ALL_SAMPLES=($(set -x; bcftools query --list-samples "${VCF}"))
#[[ ${#ALL_SAMPLES[@]} -lt 1 ]] && (echo "Fewer than 1 sample found in ${VCF/.vcf.gz/_biallelic.vcf.gz}" >&2; exit 1)
declare -a NEW_SAMPLES=($(set -x; bcftools query --list-samples "${VCF/.vcf.gz/_biallelic.vcf.gz}"))
[[ ${#NEW_SAMPLES[@]} -lt 2 ]] && (echo "Fewer than 2 samples found in ${VCF}" >&2; exit 1)
echo "Splitting ${VCF} into ${#NEW_SAMPLES[@]} sample VCF files" >&2

parallel --verbose "SplitVCF ${VCF/.vcf.gz/_biallelic.vcf.gz} ${OUTDIR} {}" ::: ${NEW_SAMPLES[@]}

parallel --verbose "generateFA ${GENOME} ${OUTDIR} {} 1" ::: ${NEW_SAMPLES[@]}

parallel --verbose "generateFA ${GENOME} ${OUTDIR} {} 2" ::: ${NEW_SAMPLES[@]}

parallel --verbose "keepHet {} ${OUTDIR}" ::: ${NEW_SAMPLES[@]}
