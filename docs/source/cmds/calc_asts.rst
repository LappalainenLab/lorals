:code:`calc_asts`
=================

    Calculate transcript counts assigned to each haplotype

Arguments
---------

.. code::

    input/output options:
        -b|--bam        BAM file containing RNA-seq reads
        -i|--input      Input ASE table
        -o|--output     Name of output ASTS file; defaults to lorals_out/asts.tsv

    asts options:
        -m|--mode           ASTS mode, choose from length, quant; defaults to length
        -x|--transcripts    BAM file aligned to transcriptome; used when 'mode' is
                            set to 'quant'
        -w|--window         Window around a variant to calculate number of matches;
                            defaults to 5
        -t--threshold       Minimum number of matches in window around the variant;
                            defaults to 8

    filter options:
        -c|--coverage           Minimum overall coverage; defaults to 10
        -a|--allelic-coverage   Minimum coverage per allele; defaults to 5
        -q|--mapq mapping       Minimum mapping quality; defaults to 10
        --filter                For annotated ASE tables, filter based on warnings provided.
                                Choose one or more from bl (blacklisted), gt (genotype),
                                mm (multimapping), other, or indel; pass '--filter' with no
                                extra arguments for all filters

    utility options:
        -v|--verbosity  Verbosity level, choose from debug, info, warning, error, critical;
                        defaults to info

Input Formats
-------------

ASE table
    See the output of ``calc_ase`` or ``annotate_ase``

BAM file
    Standard `BAM file <http://samtools.github.io/hts-specs/SAMv1.pdf>`_

Output
------

A tab-separated table with the following columns:

...
