#!/usr/bin/env python

"""Command-line scripts included in laurels"""

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
    __LAURELS_SETUP__
except NameError:
    import pandas
    import pybedtools
    if sys.version_info.major == 2:
        from . import asts
        from . import fancy_logging
    else:
        from laurels import asts
        from laurels import fancy_logging


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
        version=''
    )
    parser.add_argument( # Help action
        '-h',
        '--help',
        action=argparse._HelpAction,
        help=argparse.SUPPRESS
    )
    return parser


def _greeter(): # type: (None) -> None
    logging.info("Welcome to %s %s!", os.path.basename(sys.argv[0]), VERSION)
    logging.info("Author: Dafni Glinos (dglinos@nygenome.org)")


def calc_asts(*args): # type: (Optional[List[str, ...]]) -> None
    """Calculate ASTS"""
    parser = _common_parser() # type: argparse.ArgumentParser
    parser.add_argument( # VCF file
        '-f',
        '--vcf',
        dest='vcf',
        type=str,
        required=True,
        metavar='in.vcf',
        help="Genotype VCF for sample"
    )
    parser.add_argument( # BAM file
        '-b',
        '--bam',
        dest='bam',
        type=str,
        required=True,
        metavar='in.bam',
        help="BAM file containing RNA-seq reads"
    )
    if not sys.argv[1:]:
        sys.exit(parser.print_help())
    args = vars(parser.parse_args(*args)) # type: Dict[str, Any]
    fancy_logging.configure_logging(level=args['verbosity'])
    _greeter()
    #   Create a BED-like pilup from VCF
    bpileup = asts.vcf_bpileup(vcffile=args['vcf']) # type: str
    #   Calculate ASE
    ifh = asts.bed_intersect(afile=bpileup, bfile=args['bam']) # type: pybedtools.bedtool.BedTool


def annotate_ase(*args): # type: (Optional[List[str, ...]]) -> None
    """Annotate ASE"""
    pass


def fetch_haplotype(*args): # type: (Optional[List[str, ...]]) -> None
    """..."""
    pass
