:code:`annotate_ase`
====================

    Annotate an ASE table with warning flags and gene information

Arguments
---------

.. code::

    input/output options:
        -i|--input      Input ASE table
        -b|--bed        Reference BED file
        -f|--vcf        Genotype VCF for sample
        -o|--output     Name of output ASE file; defaults to lorals_out/ase_annotated.tsv

    blacklist options:
        --blacklist blacklist.bed           Blacklist BED file; defaults to included file
        --genotype genotype_warning.bed     Genotype warning BED file
        --mapping multi_mapping.bed         BED file with multi-mapping regions;
                                            defaults to included file

    ase stats options:
        -c|--coverage       Minimum coverage for a site to be included; defaults to 20
        -n|--binomial-null  For binomial test, the null ref ratio to test against,
                            pass 'auto' to auto-calculate null ref ratio; defaults to 0.5
        -m|--method         Method for calculating biomial null, choose from mean, median,
                            global; defaults to mean
        --other-threshold   Threshold for issuing an other allele warning; defaults to 0.8
        --indel-threshold   Threshold for issuing an indel warning; defaults to 0.2

    utility options:
        -v|--verbosity  Verbosity level, choose from debug, info, warning, error, critical;
                        defaults to info

Input Formats
-------------

ASE table
    See the output of ``calc_ase``

BED files
    Standard `BED file format`_: requires at least a four-column BED file

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

GENOTYPE
    The genotype for this variant

GENE_ID
    A comma-separated list of gene identifiers for this variant

GENOTYPE_WARNING
    Was this variant flagged for being in a region prone to genotyping errors? ``1`` for yes, ``0`` for no

BLACKLIST
    Was this variant flagged for being in a blacklisted region? ``1`` for yes, ``0`` for no

MULTI_MAPPING
    Was this variant flagged for being in a region prone to multimapping? ``1`` for yes, ``0`` for no

OTHER_ALLELE_WARNING
    Was this variant flagged for having a high ratio of other alleles? ``1`` for yes, ``0`` for no

HIGH_INDEL_WARNING
    Was this variant flagged for having a high ratio of indels? ``1`` for yes, ``0`` for no

NULL_RATIO
    Null reference ratio used for the binomial test

BINOM_P
    `p`-value for the binomial test

BINOM_P_ADJUSTED
    FDR-adjusted `p`-value



.. _BED file format: http://genome.cse.ucsc.edu/FAQ/FAQformat.html#format1
