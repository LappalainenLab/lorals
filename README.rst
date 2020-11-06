LORALS
=======

A Python package for allele-specific analyses in long-read data. Written by Dafni Glinos

Dependencies
------------

LoRALS depends on Python 3, SAMTools, and Bedtools. All other dependencies are installed automatically by pip

Installation
------------

To install, run the following command

.. code:: bash

    python3 setup.py install

Please note, to install system-wide, you must have sudo access.
If you are not installing system-wide, please install within your ``$PATH`` and ``$PYTHONPATH``
or use virtualenv

Background
------------
The method behind LORALS can be found in this preprint:

It utilises the

Usage
------------

.. code:: bash

    calc_ase -b in.bam -f in.vcf [-o /path/to/output] [-w window]
                 [-t threhsold] [-c coverage] [-a allelic coverage]
                 [-q mapping quality]

Calculates the allelic coverage of each variant. It only requires a genome aligned bam file and a phased VCF file.

The most efficient way to carry out ASTS analysis is to first run this command and filter out the variants that get tagged by
different flags.

For additional options run with --help

.. code:: bash

    annotate_ase

Annotates the output of calc_asts based on five filters and assigns it a gene. Make sure to provide a gene coordinates
file that does not contain introns if you want to avoid multiple genes assigned to a variant.

1. Ratio of reads containing indels within the variant used for ASE to the total number of reads. If you don't want to use
this flag you can set it to 0.
2. Ratio of other alleles to the REF or the ALT that are found at the variant site used for ASE. If you don't want to use
this flag you can set it to 0.
Optional
3. The variant falls within the ENCODE blacklist region. The expected file is in BED format. For ease we provide one such
file for hg38 which you can replace with any other file you like
4. The variant falls  within a multi-mapping region. The expected file is in BED format. For ease we provide one such
file for hg38 which you can replace with any other file you like
5. The variant falls within a region that is potentially wrongly assumed to be heterozygous or where the imputed genotype
is ambiguous. The expected file is in BED format.

The output file can be used as it is for allele specific expression that is per variant. If you want to carry allele specific expression
based on the exact reads assigned to a transcript please look into process_ase

.. code:: bash

    calc_asts -m quant

Calculates the number of reads containing the REF or ALT allele assigned to each transcript.
It requires the user to have aligned the reads to the relevant transcriptome and provide the alignments in BAM format.

.. code:: bash

    process_ase

Assigns a gene to each transcript and adds up all the transcript counts per gene for the REF and the ALT allele and
performs a binomial test per gene, followed by FDR correction. It currently selects the top variant per gene based
on the total number of reads. If you want to disable this function you should use X flag

.. code:: bash

    process_asts

Assigns a gene to each transcript and performs chi-square per gene followed by fdr correction. It currently selects the
top variant per gene based on the total number of reads. If you want to disable this function you should use X flag.

Note that chi-square test statistic is not reliable with low counts, we therefore set the default min. number of reads
for a transcript to 10.

.. code:: bash

    calc_asts -m length

In case the exact transcriptome is not readily available we provide this alternative ASTS analysis. Here all the
distribution of the reads overlapping the REF allele are compared to the distribution of the reads overlapping the ALT
allele.

The user can either get a summary result where XX test is performed or get the lengths per variant to carry the test of
their choice.

.. code:: bash

    fetch_haplotype -b in.bam -t transcripts.bam -s snp.txt

This script output the reads that overlap a specific SNP per haplotype and transcript. They can be useful for visualisation
using IGV or any other software.


Optional alignment steps if high reference bias is observed.


.. code:: bash

    process_vcf

The pipeline requires the VCFs to only contain a single individual and for optimal performance to only
include heterozygous variants. We provide this script in order to obtain such a VCF.
This script will perform these actions:

1. Filter VCF to only contain biallelic variants
2. Split a VCF containing records for multiple individuals into one VCF per individual and tabix the files
3. For each sample create two fasta ref files for each haplotype
4. For each sample VCF only keep het variants

.. code:: bash

    make_new_vcf

It uses an aligned bam file to correct the phased haplotypes in a vcf file.
This VCF file is then used to generate two haplotype specific genome references.

.. code:: bash

    hap_aligner

Aligns reads to each of the two genomes using minimap2, selects the best aligned read of the two based on the MAPQ score.
In case of ties it randomly selects an equal proportion from each of the two alignments.
It then converts the aligned minimap2 `sam` output to `bam` format and indexes the reads.

Alternatively, the user can align the reads themselves with their aligner of choice.
