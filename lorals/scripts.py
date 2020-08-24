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
    if sys.version_info.major == 2:
        from . import ase
        from . import asts
        from . import utils
        from . import features
        from . import fancy_logging
        from .fancy_logging import fmttime
    else:
        from lorals import ase
        from lorals import asts
        from lorals import utils
        from lorals import features
        from lorals import fancy_logging
        from lorals.fancy_logging import fmttime


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
    io_opts.add_argument(
        '-o',
        '--output',
        dest='output',
        type=str,
        required=False,
        default=os.path.join(os.getcwd(), 'lorals_out', 'asts'),
        metavar='/path/to/output',
        help="Directory and prefix of output files; defaults to %(default)s"
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
    # filter_opts.add_argument
    _common_opts(parser=parser, group='utility options', version=VERSION)
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    args = vars(parser.parse_args(*args)) # type: Dict[str, Any]
    fancy_logging.configure_logging(level=args['verbosity'])
    _ = utils.where('bedtools') # type: _
    _greeter()
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
    ase_stats = {ase.allelic_stats(var=var, bamfile=args['bam']) for var in features.iter_var(df=bpileup)} # type: Set[ase.AllelicStat, ...]
    ase_stats = tuple(sorted(filter(None, ase_stats), key=lambda x: (x.contig, x.position))) # type: Tuple[ase.AllelicStat, ...]
    with open('%s_ase_raw.tsv' % args['out'], 'w') as ofile: # type: file
        logging.info("Saving raw ASE results to %s", ofile.name)
        raw_start = time.time() # type: float
        ofile.write('\t'.join(getattr(ase.AllelicStat, '_fields')))
        ofile.write('\n')
        ofile.flush()
        for stat in ase_stats: # type: ase.AllelicStat
            ofile.write('\t'.join(map(str, stat)))
            ofile.write('\n')
            ofile.flush()
    logging.debug("Writing raw ASE results took %s seconds", round(time.time() - raw_start, 3))
    ase_stats = ase.filter_stats( # type: Tuple[ase.AllelicStat, ...]
        stats=ase_stats,
        total_coverage=args['coverage'],
        allelic_coverage=args['allelic']
    )


def annotate_ase(*args): # type: (Optional[List[str, ...]]) -> None
    """Annotate ASE"""
    pass


def fetch_haplotype(*args): # type: (Optional[List[str, ...]]) -> None
    """..."""
    pass
