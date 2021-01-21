=======
LORALS
=======

A Python package for allele-specific analyses in long-read data. Written by Dafni Glinos

Dependencies
============

LoRALS depends on Python 3, SAMTools, and Bedtools. All other dependencies are installed automatically by pip

Installation
============

To install, run the following command

.. code:: bash

    python3 setup.py install

Please note, to install system-wide, you must have sudo access.
If you are not installing system-wide, please install within your ``$PATH`` and ``$PYTHONPATH``
or use virtualenv

Background
============
The method behind LORALS can be found in this preprint:

Please cite this paper when using any part of this method.

.. image:: ../images/lorals_pipeline.svg
      :class: with-shadow
      :width: 200px

      Figure 1. Suggested pipeline of allelic read analysis using LORALS
      

The pipeline of analysis provided by LORALS is outlined on Figure 1. It is not necessary to run these steps in succession, as LORALS is modular and each step can accept different types of inputs. Please see specific functions for more information.

Usage
============

Calculate per variant read counts
------------------------------------

.. code:: bash

    calc_ase -b in.bam -f in.vcf [-o /path/to/output] [-w window]
                 [-t threhsold] [-c coverage] [-a allelic coverage]
                 [-q mapping quality]

Calculates the allelic coverage of each variant. It only requires a genome aligned bam file and a phased VCF file.

For additional options run with --help

.. code:: bash

    annotate_ase -i in_ase.tsv -b ref.bed [-f in.vcf] [-o /path/to/output]
                    [--blacklist blacklist.bed]
                    [--genotype genotype_warning.bed]
                    [--mapping multi_mapping.bed] [-c COVERAGE] [-n BINOMIAL]
                    [-m {mean,median,global}] [--other-threshold threshold]
                    [--indel-threshold threshold] [-v level]

Annotates the output of calc_ase based on five criteria and assigns it a gene. Make sure to provide a gene coordinates
file that does not contain introns if you want to avoid multiple genes assigned to a variant.

1. indel-threshold: Ratio of reads containing indels within the variant used for ASE to the total number of reads. If you don't want to use this flag you can set it to 0. The default is 0.2.
2. other-threshold: Ratio of REF and ALT containing reads to the total number of reads covering the variant site used for ASE.  If you don't want to use this flag you can set it to 0. The default is 0.8. 
3. blacklist (optional): The variant falls within the ENCODE blacklist region. The expected file is in BED format. For ease we provide one such file for hg38 which you can replace with any other file you like
4. mapping (optional): The variant falls  within a multi-mapping region. The expected file is in BED format. For ease we provide one such file for hg38 which you can replace with any other file you like
5. genotype (optional): The variant falls within a region that is potentially wrongly assumed to be heterozygous or where the imputed genotype is ambiguous. The expected file is in BED format.

The output file can be used as it is for allele specific expression, calculated per variant. If you want to carry allele specific expression
based on the exact reads assigned to a transcript please look into process_ase.

The most efficient way to carry out ASTS analysis is to first run calc_ase, followed by annotate_ase. You can then filter out the variants that get
tagged by different flags and use the filtered file as an input to calc_asts.

Calculate transcript counts assigned to each haplotype
--------------------------------------------------------

.. code:: bash

    calc_asts -m quant -b in.bam -i ase.tsv [-o /path/to/output]
              [-x transcripts.bam] [-w window] [-t threhsold] [-c coverage]
              [-a allelic coverage] [-q mapping quality] [-v level]

Calculates the number of reads containing the REF or ALT allele assigned to each transcript.
It requires the user to have aligned the reads to the relevant transcriptome and provide the alignments in BAM format.

.. code:: bash

    usage: process_asts -i asts.tsv [asts.tsv ...] -g genes.tsv [-o outdir]
                    [-r min reads per gene] [-t min reads per transcript]
                    [-v level]

Assigns a gene to each transcript. It then (1) adds up all the transcript counts per gene for the REF and the ALT allele and
performs a binomial test per gene, followed by FDR correction. This is the ASE final file. (2) It performs chi-square per gene
across the transripts, followed by fdr correction. This is the ASTS quant final file.

Note that chi-square test statistic is not reliable with low counts, we therefore set the default min. number of reads
for a transcript (-t) to 10.

It currently selects the top variant per gene based on the total number of reads. If you want to disable this function you should use X flag.

.. code:: bash

    calc_asts -m length -b in.bam -i ase.tsv [-o /path/to/output]
              [-w window] [-t threhsold] [-c coverage] [--raw-lengths]
              [-a allelic coverage] [-q mapping quality] [-v level]

In case the exact transcriptome is not readily available we provide this alternative ASTS analysis. Here the
distribution of the reads overlapping the REF allele are compared to the distribution of the reads overlapping the ALT
allele.

The user can either get a summary result where Kolmorogov-Smirnov test is performed or get the
lengths per variant to carry the test of their choice by. using the --raw-lengths option.

.. image:: ../images/pipeline_analysis.svg
      :class: with-shadow
      :width: 200px

      Figure 2. Statistical tests perfomed for different types of analysis using LORALS

Further investigation of specific genes/snps    
--------------------------------------------------------

.. code:: bash

    fetch_haplotype -b in.bam -t transcripts.bam -s snps.tsv [-o outdir]
                    [-w window size] [-m minimum matches] [-v level]

This script output the reads that overlap a specific SNP per haplotype and transcript. They can be useful for visualisation
using IGV or any other software.

Optional alignment steps
--------------------------------------------------------

.. code:: bash

    process_vcf.sh

We provide this script in order to obtain a per-individual VCF file, filtered to only
include heterozygous SNP variants. This script will perform these actions:

1. Filter VCF to only contain biallelic variants
2. Split a VCF containing records for multiple individuals into one VCF per individual and tabix the files
3. For each sample create two fasta ref files for each haplotype
4. For each sample VCF only keep het variants

.. code:: bash

    hap_aligner.sh

Aligns reads to each of the two genomes using minimap2, selects the best aligned read of the two based on the MAPQ score.
In case of ties it randomly selects an equal proportion from each of the two alignments.
It then converts the aligned minimap2 `sam` output to `bam` format and indexes the reads.

Alternatively, the user can align the reads themselves with their aligner of choice.

.. code:: bash

    make_new_vcf.sh

It uses an aligned bam file to correct the phased haplotypes in a vcf file.
This VCF file is then used to generate two haplotype specific genome references.
It requires bcftools, GATK, HAPCUT2 and `HapCUT2VCF.py <https://github.com/liangjiaoxue/PythonNGSTools/blob/master/HapCUT2VCF.py>`_ to be in your path.
