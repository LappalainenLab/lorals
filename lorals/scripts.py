#!/usr/bin/env python

"""Command-line scripts included in LoRALS"""

from __future__ import division
from __future__ import print_function

__all__ = [] # type: List

import os
import sys
import time
import logging
import argparse
import subprocess

import pkg_resources

try:
    __LORALS_SETUP__
except NameError:
    import pandas
    import pybedtools
    from scipy.stats import binom_test
    if sys.version_info.major == 2:
        from . import ase
        from . import asts
        from . import maths
        from . import utils
        from . import annotate
        from . import features
        from . import fancy_logging
    else:
        from lorals import ase
        from lorals import asts
        from lorals import maths
        from lorals import utils
        from lorals import annotate
        from lorals import features
        from lorals import fancy_logging


VERSION = ''

def _common_parser(): # type: (None) -> argparse.ArgumentParser
    parser = argparse.ArgumentParser(add_help=False) # type: argparse.ArgumentParser
    parser.add_argument( # type: Verbosity level
        '-v',
        '--verbosity',
        dest='verbosity',
        type=str,
        required=False,
        default='info',
        choices=fancy_logging.LOG_LEVELS.keys(),
        metavar='level',
        help="Verbosity level, choose from %(choices)s; defaults to %(default)s"
    )
    parser.add_argument( # Version action
        '--version',
        action='version',
        version=VERSION
    )
    parser.add_argument( # Help action
        '-h',
        '--help',
        action=argparse._HelpAction,
        help=argparse.SUPPRESS
    )
    return parser


def _common_opts(parser, group=None, version=None): # type: (argparse.ArgumentParser, Optional[str], Optional[str]) -> None
    if group:
        group = parser.add_argument_group(title=group) # type: argparse._ArgumentGroup
    else:
        group = parser # type: argparse.ArgumentParser
    group.add_argument( # Verbosity
        '-v',
        '--verbosity',
        dest='verbosity',
        type=str,
        required=False,
        default='info',
        choices=fancy_logging.LOG_LEVELS.keys(),
        metavar='level',
        help="Verbosity level, choose from %(choices)s; defaults to %(default)s"
    )
    if version:
        group.add_argument( # Version action
            '--version',
            action='version',
            version=version
        )
    if not parser.add_help:
        parser.add_argument( # Help action
            '-h',
            '--help',
            action=argparse._HelpAction,
            help=argparse.SUPPRESS
        )
    return None


def _greeter(): # type: (None) -> None
    logging.info("Welcome to %s %s!", os.path.basename(sys.argv[0]), VERSION)
    logging.info("Author: Dafni Glinos (dglinos@nygenome.org)")


def calc_asts(*args): # type: (Optional[List[str, ...]]) -> None
    """Calculate ASTS"""
    # parser = _common_parser() # type: argparse.ArgumentParser
    parser = argparse.ArgumentParser(add_help=False)
    io_opts = parser.add_argument_group(title="input/output options") # type: argparse._ArgumentGroup
    io_opts.add_argument( # BAM file
        '-b',
        '--bam',
        dest='bam',
        type=str,
        required=True,
        metavar='in.bam',
        help="BAM file containing RNA-seq reads"
    )
    io_opts.add_argument( # VCF file
        '-f',
        '--vcf',
        dest='vcf',
        type=str,
        required=True,
        metavar='in.vcf',
        help="Genotype VCF for sample"
    )
    io_opts.add_argument( # Output prefix
        '-o',
        '--output',
        dest='out',
        type=str,
        required=False,
        default=os.path.join(os.getcwd(), 'lorals_out', 'asts'),
        metavar='/path/to/output',
        help="Directory and prefix of output files; defaults to %(default)s"
    )
    asts_opts = parser.add_argument_group(title="asts options") # type: argparse._ArgumentGroup
    asts_opts.add_argument( # ASTS mode
        '-m',
        '--mode',
        dest='mode',
        type=str,
        required=False,
        default='length',
        choices=('length', 'quant'),
        metavar='mode',
        # help="ASTS mode, choose from %(choices)s; defaults to %(default)s"
        help=argparse.SUPPRESS
    )
    asts_opts.add_argument( # Transcriptome-aligned BAM files
        '-x',
        '--transcripts',
        dest='flair',
        type=str,
        required=False,
        metavar='transcripts.bam',
        # help="BAM file aligned to transcriptome; used when 'mode' is set to 'quant'"
        help=argparse.SUPPRESS
    )
    asts_opts.add_argument( # Match window
        '-w',
        '--window',
        dest='window',
        type=int,
        required=False,
        default=5,
        metavar='window',
        help="Window around a variant to calculate number of matches; defaults to %(default)s"
    )
    asts_opts.add_argument( # Minimum match threshold
        '-t',
        '--threshold',
        dest='threshold',
        type=int,
        required=False,
        default=8,
        metavar="threhsold",
        help="Minimum number of matches in window around the variant; defaults to %(default)s"
    )
    filter_opts = parser.add_argument_group(title="filter options") # type: argparse._ArgumentGroup
    filter_opts.add_argument( # Minimum coverage
        '-c',
        '--coverage',
        dest='coverage',
        type=int,
        required=False,
        default=10,
        metavar='coverage',
        help="Minimum overall coverage; defaults to %(default)s"
    )
    filter_opts.add_argument( # Minimum coverage per allele
        '-a',
        '--allelic-coverage',
        dest='allelic',
        type=int,
        required=False,
        default=5,
        metavar='allelic coverage',
        help="Minimum coverage per allele; defaults to %(default)s"
    )
    filter_opts.add_argument( # Minimum mapping quality
        '-q',
        '--mapq',
        dest='mapq',
        type=int,
        required=False,
        default=10,
        metavar='mapping quality',
        help="Minimum mapping quality; defaults to %(default)s"
    )
    _common_opts(parser=parser, group='utility options', version=VERSION)
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    args = vars(parser.parse_args(*args)) # type: Dict[str, Any]
    fancy_logging.configure_logging(level=args['verbosity'])
    _ = utils.where('bedtools') # type: _
    _greeter()
    if args['mode'] == 'quant' and not args['flair']:
        parser.error("'-t|--transcripts' must be supplied when 'mode' is 'quant'")
    args['vcf'] = utils.fullpath(path=args['vcf']) # type: str
    args['bam'] = utils.fullpath(path=args['bam']) # type: str
    args['out'] = utils.fullpath(path=os.path.splitext(args['out'])[0]) # type: str
    if not os.path.isdir(os.path.dirname(args['out'])):
        os.makedirs(os.path.dirname(args['out']))
    #   Create a BED-like pilup from VCF
    bpileup = asts.vcf_bpileup(vcffile=args['vcf']) # type: str
    ifh = asts.bed_intersect(afile=bpileup, bfile=args['bam']) # type: pybedtools.bedtool.BedTool
    os.remove(bpileup)
    bpileup = asts.deduplicate_bpileup(bpileup=ifh.fn) # type: pandas.core.frame.DataFrame
    #   Calculate ASE
    ase_stats = {ase.allelic_stats(var=var, bamfile=args['bam'], window=args['window'], match_threshold=args['threshold']) for var in features.iter_var(df=bpileup)} # type: Set[ase.AllelicStat, ...]
    ase_stats = tuple(sorted(filter(None, ase_stats), key=lambda x: (x.chrom, x.position))) # type: Tuple[ase.AllelicStat, ...]
    with open('%s_ase.tsv' % args['out'], 'w') as ofile: # type: file
        logging.info("Saving raw ASE results to %s", ofile.name)
        ofile.write('\t'.join(ase.AllelicStat.HEADER))
        ofile.write('\n')
        ofile.flush()
        for stat in ase_stats: # type: ase.AllelicStat
            ofile.write(str(stat))
            ofile.write('\n')
            ofile.flush()
    ase_stats = ase.filter_stats( # type: Tuple[ase.AllelicStat, ...]
        stats=ase_stats,
        total_coverage=args['coverage'],
        allelic_coverage=args['allelic']
    )
    #   Calculate ASTS
    logging.info("Calculating ASTS")
    if args['mode'] == 'length':
        header = header = getattr(asts.LengthStat, '_fields') # type: Tuple[str, ...]
        # asts_stats = tuple((asts.asts_length(var=features.Bpileup.fromstat(stat=stat), bamfile=args['bam'], window=args['window'], match_threshold=args['threshold']) for stat in ase_stats)) # type: Tuple[asts.LengthStat, ...]
        asts_stats = tuple((asts.asts_length(var=stat, bamfile=args['bam'], window=args['window'], match_threshold=args['threshold']) for stat in ase_stats)) # type: Tuple[asts.LengthStat, ...]
    else:
        msg = "Unknown mode: %s" % args['mode']
        logging.critical(msg)
        raise ValueError(msg)
    with open("%(out)s_asts_%(mode)s.tsv" % {'out': args['out'], 'mode': args['mode']}, 'w') as ofile: # type: file
        logging.info("Saving ASTS results to %s", ofile.name)
        ofile.write('\t'.join(header))
        ofile.write('\n')
        ofile.flush()
        for stat in asts_stats: # type: Union[asts.LengthStat]
            ofile.write('\t'.join(map(str, stat)))
            ofile.write('\n')
            ofile.flush()
    logging.info("Done")


def annotate_ase(*args): # type: (Optional[List[str, ...]]) -> None
    """Annotate ASE"""
    parser = argparse.ArgumentParser( # type: argparse.ArgumentParser
        add_help=False,
        formatter_class=argparse.RawTextHelpFormatter
    )
    io_opts = parser.add_argument_group(title="input/output options") # type: argparse._ArgumentGroup
    io_opts.add_argument( # Input ASE file
        '-i',
        '--input',
        dest='input',
        type=str,
        required=True,
        metavar='in_ase.tsv',
        help="Input ASE table"
    )
    io_opts.add_argument( # BED file
        '-b',
        '--bed',
        dest='bed',
        type=str,
        required=True,
        metavar='ref.bed',
        help="Reference BED file"
    )
    io_opts.add_argument( # VCF file
        '-f',
        '--vcf',
        dest='vcf',
        type=str,
        required=False,
        metavar='in.vcf',
        help="Genotype VCF for sample"
    )
    io_opts.add_argument( # Output file
        '-o',
        '--output',
        dest='output',
        type=str,
        required=False,
        default='ase_annotated.tsv',
        metavar='ase_annotated.tsv',
        help="Name of output ASE file; defaults to %(default)s"
    )
    blacklist_opts = parser.add_argument_group(title='blacklist options') # type: argparse._ArgumentGroup
    blacklist_opts.add_argument( # Blacklist
        '--blacklist',
        dest='blacklist',
        type=str,
        required=False,
        default=pkg_resources.resource_filename('lorals', 'blacklists/hg38-blacklist.v2.bed'),
        metavar='blacklist.bed',
        help="Blacklist BED file; defaults to %(default)s"
    )
    blacklist_opts.add_argument( # Genotype warning
        '--genotype',
        dest='warning',
        type=str,
        required=False,
        default=pkg_resources.resource_filename('lorals', 'blacklists/GTEX_Q2AG_braincerebellarhemisphere_illumina_GT_warning.bed'),
        metavar='genotype_warning.bed',
        help="Genotype warning BED file; defaults to %(default)s"
    )
    blacklist_opts.add_argument( # Multimapping
        '--mapping',
        dest='mapping',
        type=str,
        required=False,
        default=pkg_resources.resource_filename('lorals', 'blacklists/wgEncodeCrgMapabilityAlign100mer.hg38.mappingover5.bed'),
        metavar='multi_mapping.bed',
        help="BED file with multi-mapping regions; defaults to %(default)s"
    )
    stats_opts = parser.add_argument_group(title="ase stats options") # type: argparse._ArgumentGroup
    stats_opts.add_argument( # Proportion cutoff
        '-p',
        '--proportion-cutoff',
        dest='cutoff',
        type=float,
        required=False,
        default=0.05,
        help="Maximum proportion of reads arising from non ref/alt read for variant to be used for genotype warning test; defaults to %(default)s"
    )
    stats_opts.add_argument( # Coverage
        '-c',
        '--coverage',
        dest='coverage',
        type=int,
        required=False,
        default=20,
        help="Minimum coverage for a site to be included; defaults to %(default)s"
    )
    binom_excl = stats_opts.add_mutually_exclusive_group()
    binom_excl.add_argument( # Optional binomial NULL
        '-n',
        '--binomial-null',
        dest='binomial',
        type=float,
        required=False,
        default=None,
        help="For binomail test, the null ref ratio to test against; defaults to auto-calculated null ref ratio"
    )
    binom_excl.add_argument( # Method for calculating binomial null
        '-m',
        '--method',
        type=str,
        required=False,
        default='mean',
        choices=('mean', 'median', 'global'),
        help="Method for calculating biomial null, choose from %(choices)s; defaults to %(default)s"
    )
    _common_opts(parser=parser, group='utility options', version=VERSION)
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    args = vars(parser.parse_args(*args)) # type: Dict[str, Any]
    fancy_logging.configure_logging(level=args['verbosity'])
    _ = utils.where("bedtools") # type: _
    _greeter()
    args['input'] = utils.fullpath(path=args['input'])
    if not os.path.exists(args['input']):
        msg = "Cannot find input file %s" % args['input'] # type: str
        logging.critical(msg)
        raise ValueError(msg)
    my_open = utils.find_open(filename=args['input']) # type: function
    ase_stats = list() # type: List[annotate.AnnotatedStat]
    with my_open(args['input'], 'rt') as ifile:
        logging.info("Reading in input ASE from %s", args['input'])
        read_start = time.time() # type: float
        for line in ifile: # type: str
            if line.startswith('#') or line.startswith(ase.AllelicStat.HEADER[0]):
                continue
            ase_stats.append(annotate.AnnotatedStat.fromstring(string=line))
    ase_stats = tuple(ase_stats) # type: Tuple[annotate.AnnotatedStats, ...]
    logging.debug("Reading input ASE took %s seconds", fancy_logging.fmttime(start=read_start))
    if args['vcf']:
        ase_stats = annotate.annotate_genotypes(stats=ase_stats, vcffile=args['vcf']) # type: Tuple[annotate.AnnotatedStat]
    ase_stats = annotate.annotate_genes(stats=ase_stats, bedfile=args['bed']) # type: Tuple[annotate.AnnotatedStat, ...]
    #   Add blacklist, genotype warning, and mappability annotations
    logging.info("Annotating blacklist status")
    blacklist = annotate.annotate_bed(stats=ase_stats, bedfile=args['blacklist']) # type: Tuple[str, ...]
    logging.info("Annotating genotype warning")
    gt_warning = annotate.annotate_bed(stats=ase_stats, bedfile=args['warning']) # type: Tuple[str, ...]
    logging.info("Annotating multi mapping")
    multi_mapping = annotate.annotate_bed(stats=ase_stats, bedfile=args['mapping']) # type: Tuple[str, ...]
    logging.info("Adding annotations")
    for i in range(len(ase_stats)): # type: int
        ase_stats[i].blacklisted = ase_stats[i].default in blacklist # type: bool
        ase_stats[i].warning = ase_stats[i].default in gt_warning # type: bool
        ase_stats[i].multi_mapping = ase_stats[i].default in multi_mapping # type: bool
    #   Get bias stats
    if args['binomial'] is None:
        logging.info("Calculating binomial null ratio")
        stats_bias = tuple(filter(
            lambda x: not x.blacklisted and not x.warning and not x.multi_mapping and not x.other_warning and not x.indel_warning and x.chrom not in ('chrX', 'chrY'),
            ase_stats
        ))
        bias_stats = annotate.bias_stats(
            stats=stats_bias,
            method=args['method'],
            coverage=args['coverage']
        )
        logging.info("Determining null ratios")
        for i, stat in enumerate(ase_stats): # type: int, annotate.AnnotatedStat
            ase_stats[i].null_ratio = bias_stats.get(stat.ref + stat.alts, float('nan'))
    else:
        for i in range(len(ase_stats)): # type: int
            ase_stats[i].null_ratio = args['binomial']
    #   Calculate binomial p-value
    logging.info("Calculating binomial p-value")
    for i, stat in enumerate(ase_stats): # type: i, annotate.AnnotatedStat
        ase_stats[i].pvalue = binom_test(stat.ref_count, stat.total_count, stat.null_ratio)
    #   Calculate q-values
    logging.info("Adjusting p-values")
    stats_cov = tuple(filter(
        lambda x: not x.multi_mapping and not x.blacklisted and not x.other_warning and not x.indel_warning and not x.warning and x.total_count >= args['coverage'] and x.chrom not in ('chrX', 'chrY'),
        ase_stats
    ))
    qvalues = maths.pvalue_adjust(tuple(x.pvalue for x in ase_stats)) # type: Tuple[float, ...]
    for i, stat in enumerate(ase_stats): # type: int, annotate.AnnotatedStat
        try:
            index = stats_cov.index(stat)
        except ValueError:
            pass
        else:
            ase_stats[i].qvalue = qvalues[index]
    #   Write the output file
    with open(args['output'], 'w') as ofile:
        logging.info("Writing output to %s", args['output'])
        ofile.write('\t'.join(annotate.AnnotatedStat.HEADER))
        ofile.write('\n')
        ofile.flush()
        for stat in ase_stats: # type: annotate.AnnotatedStat
            ofile.write(str(stat))
            ofile.write('\n')
            ofile.flush()



def fetch_haplotype(*args): # type: (Optional[List[str, ...]]) -> None
    """..."""
    pass
