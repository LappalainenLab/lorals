#!/usr/bin/env python

"""ASTS-specific utilities"""

from __future__ import division
from __future__ import print_function

__all__ = [ # type: List[str, ...]
    "bam_to_bed",
    "bed_intersect",
    "deduplicate_bpileup",
    "vcf_bpileup",
]

import os
import sys
import time
import logging
import tempfile

from collections import namedtuple, defaultdict

if sys.version_info.major == 2:
    import cigar
    import utils
    from .fancy_logging import fmttime
else:
    from lorals import cigar
    from lorals import utils
    from lorals.fancy_logging import fmttime


import pysam
import pandas
import pybedtools

from scipy import stats

LengthStat = namedtuple( # type: type
    "LengthStat",
    (
        "D",
        "pvalue",
        "contig",
        "position",
        "refAllele",
        "altAllele",
    )
)

NullResult = stats.stats.KstestResult(statistic=utils.nan, pvalue=utils.nan) # type: scipy.stats.stats.KstestResult

def asts_length(var, bamfile, window=5, match_threshold=8, min_reads=20): # type: (features.Bpileup, str, int, int) -> LengthStat
    """Calculate read length based on SNP"""
    logging.info("Getting read lengths for %s", var.name)
    lengths_start = time.time() # type: float
    reads_completed = set() # type: Set[str, ...]
    lengths = defaultdict(list) # type: Mapping[str, List[int, ...]]
    bamfile = utils.fullpath(path=bamfile) # type: str
    bamfh = pysam.AlignmentFile(bamfile) # type: pysam.libcalcalignment.AlignmentFile
    for pile in bamfh.pileup(region=var.chrom, start=var.dummy, end=var.position): # type: pysam.libcalignedsegment.PileupColumn
        if pile.pos != var.dummy:
            continue
        for pile_read in pile.pileups: # type: pysam.libcalignedsegment.PileupRead
            if pile_read.alignment.query_name in reads_completed or not pile_read.query_position:
                continue
            logging.info("Processing read %s", pile_read.alignment.query_name)
            reads_completed.add(pile_read.alignment.query_name)
            pile_window = utils.window(position=pile_read.query_position, size=window) # type: slice
            count_m = cigar.Cigar(tuples=pile_read.alignment.cigartuples)[pile_window].count('M') # type: int
            if count_m < match_threshold:
                logging.warning("Skipping read %s", pile_read.alignment.query_name)
                continue
            if pile_read.alignment.query_sequence[pile_read.query_position] == var.ref:
                lengths['ref'].append(len(pile_read.alignment.query_sequence))
            elif pile_read.alignment.query_sequence[pile_read.query_position] in var.alt:
                lengths['alt'].append(len(pile_read.alignment.query_sequence))
    if len(lengths['ref']) >= min_reads and len(lengths['alt']) >= min_reads:
        logging.info("Running KS test")
        ks = stats.ks_2samp(lengths['ref'], lengths['alt']) # type: scipy.stats.stats.KstestResult
    else:
        logging.warning("Too few hits for KS test")
        ks = NullResult # type: scipy.stats.stats.KstestResult
    logging.info("Finished getting lengths for %s", str(var))
    logging.debug("Getting read lengths took %s seconds", fmttime(start=lengths_start))
    return LengthStat(
        D=ks.statistic,
        pvalue=ks.pvalue,
        contig=var.chrom,
        position=var.position,
        refAllele=var.ref,
        altAllele=','.join(var.alt)
    )


def asts_quant():
    """..."""
    pass


def bam_to_bed(bamfile, save=False): # type: (str, bool) -> pybedtools.bedtool.BedTool
    """Convert a BAM file to a BED file"""
    logging.info("Converting %s from BAM to BED", bamfile)
    bam_bed_start = time.time()
    bamfh = pybedtools.bedtool.BedTool(bamfile) # type: pybedtools.bedtool.BedTool
    bedfh = bamfh.bam_to_bed() # type: pybedtools.bedtool.BedTool
    if save:
        ofile = os.path.splitext(bamfile)[0] + '.bed' # type: str
        logging.info("Saving resulting BED file to %s", ofile)
        bedfh = bedfh.moveto(ofile) # type: pybedtools.bedtool.BedTool
    logging.debug("Converting from BAM to BED took %s seconds", fmttime(start=bam_bed_start))
    return bedfh


def bed_intersect(afile, bfile, ofile=None): # type: (str, str, Optional[str]) -> pybedtools.bedtool.BedTool
    """Intersect two bedfiles"""
    afile = utils.fullpath(path=afile) # type: str
    bfile = utils.fullpath(path=bfile) # type: str
    clean = set() # type: Set[str]
    logging.info("Connecting to %s", afile)
    if os.path.splitext(afile)[-1] == '.bam':
        afh = bam_to_bed(bamfile=afile) # type: pybedtools.bedtool.BedTool
        clean.add(afh.fn)
    else:
        afh = pybedtools.bedtool.BedTool(afile) # type: pybedtools.bedtool.BedTool
    logging.info("Connecting to %s", bfile)
    if os.path.splitext(bfile)[-1] == '.bam':
        bfh = bam_to_bed(bamfile=bfile) # type: pybedtools.bedtool.BedTool
        clean.add(bfh.fn)
    else:
        bfh = pybedtools.bedtool.BedTool(bfile) # type: pybedtools.bedtool.BedTool
    logging.info("Intersecting two files")
    intersect_start = time.time()
    ifh = afh.intersect(bfh) # type: pybedtools.bedtool.BedTool
    if isinstance(ofile, str):
        ofile = utils.fullpath(path=ofile) # type: str
        logging.info("Saving resulting BED file as %s", ofile)
        ifh = ifh.moveto(ofile) # type: pybedtools.bedtool.BedTool
    logging.info("Number of interesected regions: %s", len(ifh))
    logging.debug("Intersection took %s seconds", fmttime(start=intersect_start))
    if clean:
        logging.debug("Cleaning temporary files")
        for i in clean: # type: str
            os.remove(i)
    return ifh


def deduplicate_bpileup(bpileup): # type: (str) -> pandas.core.frame.DataFrame
    """Deduplicate a BED file"""
    bpileup = utils.fullpath(path=bpileup) # type: str
    dedup_start = time.time() # type: float
    bed = pandas.read_csv( # type: pandas.core.frame.DataFrame
        bpileup,
        sep="\t",
        names=['contig', 'dummy', 'position', 'variantID', 'refAllele', 'altAllele']
    )
    norig = len(bed) # type: int
    bed.drop_duplicates(keep='first', inplace=True)
    ndedup = len(bed) # type: int
    if norig == ndedup:
        logging.info("No duplicates in %s", bpileup)
    else:
        logging.info("Deduplicating %s from %s to %s", bpileup, norig, len(bed))
    logging.debug("Deduplication took %s seconds", fmttime(start=dedup_start))
    return bed


def vcf_bpileup(vcffile, pileup=None): # type: (str, Optional[str]) -> str
    """Create a BED-like pileup from a VCF file"""
    vcffile = utils.fullpath(path=vcffile) # type: str
    if pileup:
        pileup = utils.fullpath(path=pileup) # type: str
        pfile = utils.find_open(filename=pileup)(pileup, 'w+b')
    else:
        pfile = tempfile.NamedTemporaryFile(delete=False)
    my_open = utils.find_open(filename=vcffile) # type: function
    with my_open(vcffile, 'rt') as vfile:
        logging.info("Reading in VCF file %s", vcffile)
        pileup_start = time.time()
        for line in vfile: # type: str
            if line.startswith('#'):
                continue
            line = line.strip().split() # type: List[str, ...]
            out = ( # type: Tuple[Union[int, str], ...]
                line[0],
                int(line[1]) - 1,
                line[1],
                line[2],
                line[3],
                line[4]
            )
            pfile.write('\t'.join(map(str, out)))
            pfile.write('\n')
            pfile.flush()
    pfile.close()
    logging.info("Finished making pileup")
    logging.debug("Building pileup took %s seconds", fmttime(start=pileup_start))
    return pfile.name
