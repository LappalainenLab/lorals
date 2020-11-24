:code:`calc_ase`
================

    Calculate per variant read counts

Arguments
---------

.. code::

    input/output options:
        -b|--bam        BAM file containing RNA-seq reads
        -f|--vcf        Genotype VCF for sample
        -o|--output     Name of output ASE file; defaults to lorals_out/ase.tsv

    ase options:
        -w|--window     Window around a variant to calculate number of matches;
                        defaults to 5
        -t|--threshold  Minimum number of matches in window around the variant;
                        defaults to 8

    utility options:
        -v|--verbosity  Verbosity level, choose from debug, info, warning, error, critical;
                        defaults to info

Input Formats
-------------

BAM file
    Standard `BAM file <http://samtools.github.io/hts-specs/SAMv1.pdf>`_

VCF file
    Bgzipped `VCF file <http://samtools.github.io/hts-specs/VCFv4.3.pdf>`_

Output
------

A tab-separated table with the following columns:

contig
    The chromosome or contig of the variant

position
    The position of the variant

variantID
    The ID given to the variant

refAllele
    The reference allele of the variant

altAllele
    The alternate alleles, joined by commas, of the variant

refCount
    The count of the reference allele for this variant

altCount
    The count of the alternate allele for this variant

totalCount
    The total expression (ref + alt) count for this variant

refIndelCount
    Ref indel count

altIndelCount
    alt indel count

otherBases
    The count of other alleles for this variant

rawDepth
    The total expression (ref + alt + other) count for this variant
