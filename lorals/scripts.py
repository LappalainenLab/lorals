#!/usr/bin/env python3

"""Command-line scripts included in LoRALS"""

import os
import sys
import time
import logging
import argparse
import subprocess

from collections import defaultdict
from typing import Any, Callable, DefaultDict, Dict, List, Optional, Set, Tuple

import pkg_resources

try:
    __LORALS_SETUP__
except NameError:
    import pysam
    import pandas
    import pybedtools
    from scipy.stats import binom_test
    from lorals import ase
    from lorals import asts
    from lorals import maths
    from lorals import utils
    from lorals import annotate
    from lorals import features
    from lorals import fancy_logging


VERSION: str = ''

__all__: List = []

def _common_parser() -> argparse.ArgumentParser:
    parser: argparse.ArgumentParser = argparse.ArgumentParser(add_help=False)
    parser.add_argument( # Verbosity level
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


def _common_opts(
    parser: argparse.ArgumentParser,
    group: Optional[str]=None,
    version: Optional[str]=None
) -> None:
    if group:
        group: argparse._ArgumentGroup = parser.add_argument_group(title=group)
    else:
        group: argparse.ArgumentParser = parser
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


def _greeter() -> None:
    logging.info("Welcome to %s %s!", os.path.basename(sys.argv[0]), VERSION)
    logging.info("Author: Dafni Glinos (dglinos@nygenome.org)")


def _read_ase(filename: str) -> 'Tuple[annotate.AnnotatedStat, ...]':
    ase_stats: List[annotate.AnnotatedStat] = list()
    with open(filename, 'rt') as ifile:
        for line in ifile: # type: str
            if line.startswith(('#', ase.AllelicStat.HEADER[0], annotate.AnnotatedStat.HEADER[0])):
                continue
            stat: annotate.AnnotatedStat = annotate.AnnotatedStat.fromstring(string=line)
            line: List[str] = line.strip().split()
            if len(line) == len(annotate.AnnotatedStat.HEADER):
                stat.geno = line[12]
                stat.geno = line[13]
                stat.warning = line[14]
                stat.blacklisted = line[15]
                stat.multi_mapping = line[16]
                # stat.other_warning = line[17]
                # stat.indel_warning = line[18]
                stat.null_ratio = line[19]
                stat.pvalue = line[20]
                stat.qvalue = line[21]
            ase_stats.append(stat)
    return ase_stats


def calc_ase(*args: Optional[List[str]]) -> None:
    parser: argparse.ArgumentParser = argparse.ArgumentParser(add_help=False)
    io_opts: argparse._ArgumentGroup = parser.add_argument_group(title="input/output options")
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
        default=os.path.join(os.getcwd(), 'lorals_out', 'ase.tsv'),
        metavar='/path/to/output',
        help="Directory and prefix of output files; defaults to %(default)s"
    )
    ase_opts: argparse._ArgumentGroup = parser.add_argument_group(title="ase options")
    ase_opts.add_argument( # Match window
        '-w',
        '--window',
        dest='window',
        type=int,
        required=False,
        default=5,
        metavar='window',
        help="Window around a variant to calculate number of matches; defaults to %(default)s"
    )
    ase_opts.add_argument( # Minimum match threshold
        '-t',
        '--threshold',
        dest='threshold',
        type=int,
        required=False,
        default=8,
        metavar="threhsold",
        help="Minimum number of matches in window around the variant; defaults to %(default)s"
    )
    _common_opts(parser=parser, group='utility_options', version=VERSION)
    if not sys.argv[1:]:
        parser.print_help(file=sys.stderr)
        raise SystemExit(1)
    args: Dict[str, Any] = vars(parser.parse_args(*args))
    fancy_logging.configure_logging(level=args['verbosity'])
    _greeter()
    for key in ('bam', 'vcf', 'out'):
        args[key]: str = utils.fullpath(path=args[key])
    os.makedirs(os.path.dirname(args['out']), exist_ok=True)
    #   Create a BED-like pileup from the VCF
    bpileup: str = asts.vcf_bpileup(vcffile=args['vcf'])
    ifh: pybedtools.bedtool.BedTool = asts.bed_intersect(afile=bpileup, bfile=args['bam'])
    os.remove(bpileup)
    bpileup: pandas.core.frame.DataFrame = asts.deduplicate_bpileup(bpileup=ifh.fn)
    #   Calculate ASE
    ase_stats: Set[ase.AllelicStat] = set()
    for var in features.iter_var(df=bpileup): # type: features.Bpileup
        ase_stats.add(ase.allelic_stats(
            var=var,
            bamfile=args['bam'],
            window=args['window'],
            match_threshold=args['threshold']
        ))
    ase_stats: Tuple[ase.AllelicStat, ...] = tuple(sorted(filter(None, ase_stats), key=lambda x: (x.chrom, x.position)))
    with open(args['out'], 'w') as ofile: # type: file
        logging.info("Saving raw ASE results to %s", ofile.name)
        ofile.write('\t'.join(ase.AllelicStat.HEADER))
        ofile.write('\n')
        ofile.flush()
        for stat in ase_stats: # type: ase.AllelicStat
            ofile.write(str(stat))
            ofile.write('\n')
            ofile.flush()


def calc_asts(*args: Optional[List[str]]) -> None:
    """Calculate ASTS"""
    parser: argparse.ArgumentParser = argparse.ArgumentParser(add_help=False)
    io_opts: argparse._ArgumentGroup = parser.add_argument_group(title="input/output options")
    io_opts.add_argument( # BAM file
        '-b',
        '--bam',
        dest='bam',
        type=str,
        required=True,
        metavar='in.bam',
        help="BAM file containing RNA-seq reads"
    )
    # io_opts.add_argument( # VCF file
    #     '-f',
    #     '--vcf',
    #     dest='vcf',
    #     type=str,
    #     required=True,
    #     metavar='in.vcf',
    #     help="Genotype VCF for sample"
    # )
    io_opts.add_argument( # Input ASE matrix
        '-i',
        '--input',
        dest='input',
        type=str,
        required=True,
        metavar='ase.tsv',
        help="Input ASE file"
    )
    io_opts.add_argument( # Output file
        '-o',
        '--output',
        dest='out',
        type=str,
        required=False,
        default=os.path.join(os.getcwd(), 'lorals_out', 'asts.tsv'),
        #default='asts.tsv',
        metavar='/path/to/output',
        #metavar='asts.tsv',
        # help="Directory and prefix of output files; defaults to %(default)s"
        help="Name of output file; defaults to %(default)s"
    )
    asts_opts: argparse._ArgumentGroup = parser.add_argument_group(title="asts options")
    asts_opts.add_argument( # ASTS mode
        '-m',
        '--mode',
        dest='mode',
        type=str,
        required=False,
        default='length',
        choices=('length', 'quant'),
        metavar='mode',
        help="ASTS mode, choose from %(choices)s; defaults to %(default)s"
        # help=argparse.SUPPRESS
    )
    asts_opts.add_argument( # Transcriptome-aligned BAM files
        '-x',
        '--transcripts',
        dest='flair',
        type=str,
        required=False,
        metavar='transcripts.bam',
        help="BAM file aligned to transcriptome; used when 'mode' is set to 'quant'"
        # help=argparse.SUPPRESS
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
    filter_opts: argparse._ArgumentGroup = parser.add_argument_group(title="filter options")
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
        parser.print_help(file=sys.stderr)
        raise SystemExit(1)
    args: Dict[str, Any] = vars(parser.parse_args(*args))
    fancy_logging.configure_logging(level=args['verbosity'])
    _ = utils.where('bedtools')
    _greeter()
    if args['mode'] == 'quant' and not args['flair']:
        parser.error("'-x|--transcripts' must be supplied when 'mode' is 'quant'")
    for k in ('input', 'out', 'flair'): # type: str
        if args[k]:
            args[k]: str = utils.fullpath(path=args[k])
    # args['vcf'] = utils.fullpath(path=args['vcf']) # type: str
    args['bam'] = utils.fullpath(path=args['bam']) # type: str
    # args['out'] = utils.fullpath(path=os.path.splitext(args['out'])[0]) # type: str
    # os.makedirs(os.path.dirname(args['out']), exist_ok=True)
    # #   Create a BED-like pilup from VCF
    # bpileup: str = asts.vcf_bpileup(vcffile=args['vcf'])
    # ifh: pybedtools.bedtool.BedTool = asts.bed_intersect(afile=bpileup, bfile=args['bam'])
    # os.remove(bpileup)
    # bpileup: pandas.core.frame.DataFrame = asts.deduplicate_bpileup(bpileup=ifh.fn)
    # #   Calculate ASE
    # ase_stats: Set[ase.AllelicStat] = {ase.allelic_stats(var=var, bamfile=args['bam'], window=args['window'], match_threshold=args['threshold']) for var in features.iter_var(df=bpileup)}
    # ase_stats: Tuple[ase.AllelicStat, ...] = tuple(sorted(filter(None, ase_stats), key=lambda x: (x.chrom, x.position)))
    # with open('%s_ase.tsv' % args['out'], 'w') as ofile: # type: file
    #     logging.info("Saving raw ASE results to %s", ofile.name)
    #     ofile.write('\t'.join(ase.AllelicStat.HEADER))
    #     ofile.write('\n')
    #     ofile.flush()
    #     for stat in ase_stats: # type: ase.AllelicStat
    #         ofile.write(str(stat))
    #         ofile.write('\n')
    #         ofile.flush()
    # ase_stats: Tuple[ase.AllelicStat, ...] = ase.filter_stats(
    #     stats=ase_stats,
    #     total_coverage=args['coverage'],
    #     allelic_coverage=args['allelic']
    # )
    ase_stats: Tuple[annotate.AnnotatedStat, ...] = _read_ase(filename=args['input'])
    #   Calculate ASTS
    logging.info("Calculating ASTS")
    if args['mode'] == 'length':
        header: Tuple[str, ...] = getattr(asts.LengthStat, '_fields')
        asts_stats: Tuple[asts.LengthStat, ...] = tuple(asts.asts_length(var=stat, bamfile=args['bam'], window=args['window'], match_threshold=args['threshold']) for stat in ase_stats)
    elif args['mode'] == 'quant':
        header: Tuple[str, ...] = getattr(asts.QuantStat, '_fields')
        asts_stats: Tuple[asts.QuantStat, ...] = tuple()
        for stat in ase_stats: # type: ase.AllelicStat
            asts_stats += asts.asts_quant(
                var=stat,
                bamfile=args['bam'],
                trans_reads=asts.reads_dict(bamfile=args['flair']),
                window=args['window'],
                match_threshold=args['threshold']
            )
    else:
        msg = "Unknown mode: %s" % args['mode']
        logging.critical(msg)
        raise ValueError(msg)
    with open("%(out)s_asts_%(mode)s.tsv" % {'out': args['out'], 'mode': args['mode']}, 'w') as ofile: # type: file
        logging.info("Saving ASTS results to %s", ofile.name)
        ofile.write('\t'.join(header))
        ofile.write('\n')
        ofile.flush()
        for stat in asts_stats: # type: Union[asts.LengthStat, asts.QuantStat]
            ofile.write('\t'.join(map(str, stat)))
            ofile.write('\n')
            ofile.flush()
    logging.info("Done")


def annotate_ase(*args: Optional[List[str]]) -> None:
    """Annotate ASE"""
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        add_help=False,
        formatter_class=argparse.RawTextHelpFormatter
    )
    io_opts: argparse._ArgumentGroup = parser.add_argument_group(title="input/output options")
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
        default=os.path.join(os.getcwd(), 'lorals_out', 'ase_annotated.tsv'),
        #default='ase_annotated.tsv',
        metavar='/path/to/output',
        #metavar='ase_annotated.tsv',
        help="Name of output ASE file; defaults to %(default)s"
    )
    blacklist_opts: argparse._ArgumentGroup = parser.add_argument_group(title='blacklist options')
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
        # default=pkg_resources.resource_filename('lorals', 'blacklists/GTEX_Q2AG_braincerebellarhemisphere_illumina_GT_warning.bed'),
        metavar='genotype_warning.bed',
        # help="Genotype warning BED file; defaults to %(default)s"
        help="Genotype warning BED file"
    )
    blacklist_opts.add_argument( # Multimapping
        '--mapping',
        dest='mapping',
        type=str,
        required=False,
        default=pkg_resources.resource_filename('lorals', 'blacklists/wgEncodeCrgMapabilityAlign100mer.hg38.mappingover5.bed.gz'),
        metavar='multi_mapping.bed',
        help="BED file with multi-mapping regions; defaults to %(default)s"
    )
    stats_opts: argparse._ArgumentGroup = parser.add_argument_group(title="ase stats options")
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
    args: Dict[str, Any] = vars(parser.parse_args(*args))
    fancy_logging.configure_logging(level=args['verbosity'])
    _ = utils.where("bedtools")
    _greeter()
    args['input']: str = utils.fullpath(path=args['input'])
    if not os.path.exists(args['input']):
        msg = "Cannot find input file %s" % args['input'] # type: str
        logging.critical(msg)
        raise ValueError(msg)
    my_open: Callable = utils.find_open(filename=args['input'])
    ase_stats: List[annotate.AnnotatedStat] = list()
    with my_open(args['input'], 'rt') as ifile:
        logging.info("Reading in input ASE from %s", args['input'])
        read_start: float = time.time()
        for line in ifile: # type: str
            if line.startswith('#') or line.startswith(ase.AllelicStat.HEADER[0]):
                continue
            ase_stats.append(annotate.AnnotatedStat.fromstring(string=line))
    ase_stats: Tuple[annotate.AnnotatedStat, ...] = tuple(ase_stats)
    logging.debug("Reading input ASE took %s seconds", fancy_logging.fmttime(start=read_start))
    if args['vcf']:
        ase_stats: Tuple[annotate.AnnotatedStat, ...] = annotate.annotate_genotypes(stats=ase_stats, vcffile=args['vcf'])
    ase_stats: Tuple[annotate.AnnotatedStat, ...] = annotate.annotate_genes(stats=ase_stats, bedfile=args['bed'])
    #   Add blacklist, genotype warning, and mappability annotations
    logging.info("Annotating blacklist status")
    blacklist: Tuple[str, ...] = annotate.annotate_bed(stats=ase_stats, bedfile=args['blacklist'])
    if args['warning']:
        logging.info("Annotating genotype warning")
        gt_warning: Tuple[str, ...] = annotate.annotate_bed(stats=ase_stats, bedfile=args['warning'])
    else:
        logging.warning("No genotype warning BED file provided")
        gt_warning: Tuple = tuple()
    logging.info("Annotating multi mapping")
    multi_mapping: Tuple[str, ...] = annotate.annotate_bed(stats=ase_stats, bedfile=args['mapping'])
    logging.info("Adding annotations")
    for i in range(len(ase_stats)): # type: int
        ase_stats[i].blacklisted: bool = ase_stats[i].default in blacklist
        ase_stats[i].warning: bool = ase_stats[i].default in gt_warning
        ase_stats[i].multi_mapping: bool = ase_stats[i].default in multi_mapping
    #   Get bias stats
    if args['binomial'] is None:
        logging.info("Calculating binomial null ratio")
        stats_bias: Tuple[annotate.AnnotatedStat, ...] = tuple(filter(
            lambda x: not x.blacklisted and not x.warning and not x.multi_mapping and not x.other_warning and not x.indel_warning and x.chrom not in ('chrX', 'chrY'),
            ase_stats
        ))
        bias_stats: Dict[str, float] = annotate.bias_stats(
            stats=stats_bias,
            method=args['method'],
            coverage=args['coverage']
        )
        logging.info("Determining null ratios")
        for i, stat in enumerate(ase_stats): # type: int, annotate.AnnotatedStat
            ase_stats[i].null_ratio: float = bias_stats.get(stat.ref + stat.alts, float('nan'))
    else:
        for i in range(len(ase_stats)): # type: int
            ase_stats[i].null_ratio: float = args['binomial']
    #   Calculate binomial p-value
    logging.info("Calculating binomial p-value")
    for i, stat in enumerate(ase_stats): # type: i, annotate.AnnotatedStat
        ase_stats[i].pvalue: float = binom_test(stat.ref_count, stat.total_count, stat.null_ratio)
    #   Calculate q-values
    logging.info("Adjusting p-values")
    stats_cov: Tuple[annotate.AnnotatedStat, ...] = tuple(filter(
        lambda x: not x.multi_mapping and not x.blacklisted and not x.other_warning and not x.indel_warning and not x.warning and x.total_count >= args['coverage'],
        ase_stats
    ))
    qvalues: Tuple[float, ...] = maths.pvalue_adjust(tuple(x.pvalue for x in ase_stats))
    for i, stat in enumerate(ase_stats): # type: int, annotate.AnnotatedStat
        try:
            index: int = stats_cov.index(stat)
        except ValueError:
            pass
        else:
            ase_stats[i].qvalue: float = qvalues[index]
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


def fetch_haplotype(*args: Optional[List[str]]) -> None:
    """..."""
    # snp_header = ( # type: Tuple[str, ...]
    #     'CHR',
    #     'POS',
    #     'ALLELE1',
    #     'ALLELE2',
    # )
    parser = argparse.ArgumentParser(add_help=False) # type: argparse.ArgumentParser
    io_opts: argparse._ArgumentGroup = parser.add_argument_group(title="input/output options")
    io_opts.add_argument( # Input BAM file
        '-b',
        '--bam',
        dest='bam',
        type=str,
        required=True,
        metavar='in.bam',
        help="BAM file containing RNA-seq reads aligned to the genome"
    )
    io_opts.add_argument( # Transcript-aligned BAM file
        '-t',
        '--transcripts',
        dest='trans',
        type=str,
        required=True,
        metavar='transcripts.bam',
        help="BAM file with reads aligned to the transcriptome"
    )
    io_opts.add_argument( # SNP table
        '-s',
        '--snps',
        dest='snps',
        type=str,
        required=True,
        metavar='snps.tsv',
        help="Tab-delimited file of SNPs be used for plotting; should contain the following columns: chrom, position, reference allele, alternate allele(s). Any header should start with a '#'"
    )
    io_opts.add_argument( # Outdir
        '-o',
        '--outdir',
        dest='outdir',
        type=str,
        required=False,
        default=os.path.join(os.getcwd(), 'lorals_haplo_bams'),
        metavar='outdir',
        help="Output directory for the resulting BAM files; defaults to %(default)s"
    )
    filter_opts: argparse._ArgumentGroup = parser.add_argument_group(title='filtering options')
    filter_opts.add_argument( # Window size
        '-w',
        '--window',
        dest='window',
        type=int,
        required=False,
        default=5,
        metavar='window size',
        help="Window around the heterozygous variant to count number of matches and mismatches; defaults to %(default)s"
    )
    filter_opts.add_argument( # Minimum number of matches
        '-m',
        '--minimum-matches',
        dest='minmatch',
        type=int,
        required=False,
        default=8,
        metavar='minimum matches',
        help="Minimum number of matches within the window around the heterozygous variant; defaults to %(default)s"
    )
    _common_opts(parser=parser, group='utility options', version=VERSION)
    if not sys.argv[1:]:
        parser.print_help(file=sys.stderr)
        raise SystemExit(1)
    args: Dict[str, Any] = vars(parser.parse_args(*args))
    fancy_logging.configure_logging(level=args['verbosity'])
    _greeter()
    # read_hash: DefaultDict[str, List] = defaultdict(list)
    for k in ('bam', 'trans', 'snps', 'outdir'): # type: str
        args[k] = utils.fullpath(path=args[k]) # type: str
    os.makedirs(args['outdir'], exist_ok=True)
    my_open = utils.find_open(args['snps']) # type: Callable
    with my_open(args['snps'], 'rt') as sfile:
        #   TODO support multiple SNPs at once
        for line in sfile: # type: str
            if not line.startswith('#'):
                break
        line: List[str] = line.strip().split('\t')
        var: features.Bpileup = features.Bpileup(
            chrom=line[0],
            position=line[1],
            ref=line[2],
            alt=line[3]
        )
    with pysam.AlignmentFile(args['trans']) as flairfh:
        logging.info("Reading transcripts from %s", args['trans'])
        flair_reads: Dict[str, str] = {read.query_name: read.reference_name for read in flairfh}
    logging.info("getting qnames")
    qnames: Tuple[asts.Qname] = asts.qnames(
        var=var,
        bamfile=args['bam'],
        trans_reads=flair_reads,
        window=args['window'],
        min_matches=args['minmatch']
    )
    logging.info("splitting qnames")
    qnames: Dict[str, Dict[str, Tuple[str, ...]]] = asts.split_qnames(qnames=qnames)
    bamfh = pysam.AlignmentFile(args['bam'])
    outfh: Dict[str, Any] = dict()
    logging.info("prepping out bams")
    for vartype in qnames: # type: str
        # outfh[vartype] = dict()
        for transcript in qnames[vartype]: # type: str
            outname: str = os.path.join(
                args['outdir'],
                '%(ts)s_%(vt)s.bam' % {'ts': transcript, 'vt': vartype.lower()}
            )
            # outfh[vartype][transcript] = pysam.Samfile(
            #     outname,
            #     'wb',
            #     template=bamfh
            # )
            outfh[transcript] = pysam.Samfile(
                outname,
                'wb',
                template=bamfh
            )
    combined_qnames: Dict[str, Tuple[str, ...]] = dict()
    for d in qnames.items(): # type: Dict[str, Tuple[str, ...]]
        combined_qnames.update(d)
    for read in bamfh.fetch(until_eof=True): # type: Any
        if read.query_name in utils.unpack(combined_qnames.values()):
            transcript = utils.dictsearch(d=combined_qnames, query=read.query_name) # type: str
            outfh[transcript].write(read)
    bamfh.close()
