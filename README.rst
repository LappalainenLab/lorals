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

Optional alignment steps


.. code:: bash

    TBD

Generates two genome references based on a phased VCF file and a fasta genome reference file.

.. code:: bash

    TBD

Aligns reads to each of the two genomes using minimap2, selects the best aligned read of the two based on the MAPQ score.
In case of ties it randomly selects an equal proportion from each of the two alignments.
It then converts the aligned minimap2 `sam` output to `bam` format and indexes the reads.

Alternatively, the user can align the reads themselves with their aligner of choice.

.. code:: bash

    calc_ase -b in.bam -f in.vcf [-o /path/to/output] [-w window]
                 [-t threhsold] [-c coverage] [-a allelic coverage]
                 [-q mapping quality]

Calculates the allelic coverage of each variant. It only requires a genome aligned bam file and a phased VCF file.

In order to carry out ASTS analysis you should first run this command and filter out the variants that get tagged by
different flags.

For additional options run with --help

.. code:: bash

    annotate_ase

Annotates the output of calc_asts based on five different filters:
1. Ratio of reads containing indels within the variant used for ASE to the total number of reads
2. Ratio of other alleles to the REF or the ALT that are found at the variant used for ASE
Optional
3. The variant falls within ENCODE blacklist region. The expected file is in BED format. For ease we provide one such
file for hg38 which you can replace with whatever file you like
4. The variant falls  within a multi-mapping region. The expected file is in BED format. For ease we provide one such
file for hg38 which you can replace with whatever file you like
5. The variant falls within a region that is potentially wrongly assumed to be heterozygous. The expected file is in
BED format.

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
