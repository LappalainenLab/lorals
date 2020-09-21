#!/usr/bin/env python

"""ASE-specific utilties"""

from __future__ import division
from __future__ import print_function

__all__ = [ # type: List[str, ...]
    'AllelicStat',
    'allelic_stats',
    'filter_stats',
]

import sys
import time
import logging

if sys.version_info.major == 2:
    from . import cigar
    from . import utils
    from .features import Bpileup
else:
    from lorals import cigar
    from lorals import utils
    from lorals.features import Bpileup


import pysam

class AllelicStat(Bpileup):

    HEADER = ( # type: Tuple[str, ...]
        'contig',
        'position',
        'variantID',
        'refAllele',
        'altAllele',
        'refCount',
        'altCount',
        'totalCount',
        'refIndelCount',
        'altIndelCount',
        'otherBases',
        'rawDepth',
    )

    @classmethod
    def fromstring(cls, string, sep='\t'): # type: (str) -> AllelicStat
        string = string.strip().split(sep) # type: List[str]
        if len(string) != len(AllelicStat.HEADER):
            raise ValueError("Incorrect number of columns")
        #   totalCount and rawDepth are calculated by AllelicStat
        return cls(
            chrom=string[0],
            position=string[1],
            name=string[2],
            ref=string[3],
            alt=string[4],
            ref_count=string[5],
            alt_count=string[6],
            ref_indel=string[8],
            alt_indel=string[9],
            other_count=string[10]
        )

    @classmethod
    def frombpileup(
            cls,
            var, # type: Bpileup
            ref_count, # type: int
            alt_count, # type: int
            other_count, # type: int
            ref_indel, # type: int
            alt_indel # type: int
    ): # type: (...) -> AllelicStat
        return cls(
            chrom=var.chrom,
            position=var.position,
            ref=var.ref,
            # alt = ','.join(var.alt),
            alt=var.alts,
            ref_count=ref_count,
            alt_count=alt_count,
            other_count=other_count,
            ref_indel=ref_indel,
            alt_indel=alt_indel,
            name=var.name
        )

    def __init__(
            self,
            chrom, # type: str
            position, # type: int
            ref, # type: str
            alt, # type: str
            ref_count, # type: int
            alt_count, # type: int
            other_count, # type: int
            ref_indel, # type: int
            alt_indel, # type: int
            name=None #type: Optional[str]
    ): # type: (...) -> None
        super(AllelicStat, self).__init__(
            chrom=chrom,
            position=position,
            ref=ref,
            alt=alt,
            name=name
        )
        self._rc = int(ref_count) # type: int
        self._ac = int(alt_count) # type: int
        self._oc = int(other_count) # type: int
        self._ric = int(ref_indel) # type: int
        self._aic = int(alt_indel) # type: int

    def __str__(self): # type: (None) -> str
        out = ( # type: Tuple[Union[int, str], ...]
            self.chrom,
            self.position,
            self.name,
            self.ref,
            self.alts,
            self.ref_count,
            self.alt_count,
            self.total_count,
            self.ref_indel,
            self.alt_indel,
            self.other_count,
            self.depth
        )
        return '\t'.join(map(str, out))

    def print_bpileup(self): # type: (None) -> str
        """Print this AllelicStat as if it was a Bpileup"""
        return super(AllelicStat, self).__str__()

    ref_count = property(fget=lambda self: self._rc, doc="Counts for reference allele")
    alt_count = property(fget=lambda self: self._ac, doc="Counts for alternate allele")
    other_count = property(fget=lambda self: self._oc, doc="Counts for other alleles")
    ref_indel = property(fget=lambda self: self._ric, doc="Reference indel counts")
    alt_indel = property(fget=lambda self: self._aic, doc="Alternate indel counts")
    total_count = property(fget=lambda self: self.ref_count + self.alt_count, doc="Total counts")
    depth = property(
        fget=lambda self: self.ref_count + self.alt_count + self.other_count,
        doc="Total depth"
    )


def allelic_stats(var, bamfile, window=5, match_threshold=8): # type: (features.Bpileup, str, int, int) -> Optional[AllelicStat]
    """..."""
    logging.info("Getting allelic reads for %s", repr(var))
    ase_start = time.time()
    reads_completed = set() # type: Set[str, ...]
    stats = dict.fromkeys(('keep_ref', 'keep_alt', 'indel_ref', 'indel_alt', 'other'), 0) # type: Dict[str, int]
    bamfile = utils.fullpath(path=bamfile) # type: str
    bamfh = pysam.AlignmentFile(bamfile) # type: pysam.libcalignmentfile.AlignmentFile
    for pile in bamfh.pileup(region=var.chrom, start=var.dummy, end=var.position): # type: pysam.libcalignedsegment.PileupColumn
        # if pile.pos != var.start:
        if pile.pos != var.dummy:
            continue
        for pile_read in pile.pileups: # type: pysam.libcalignedsegment.PileupRead
            if pile_read.alignment.query_name in reads_completed or not pile_read.query_position:
                continue
            logging.info("Processing read %s", pile_read.alignment.query_name)
            reads_completed.add(pile_read.alignment.query_name)
            pile_window = utils.window(position=pile_read.query_position, size=window) # type: slice
            count_m = cigar.Cigar(tuples=pile_read.alignment.cigartuples)[pile_window].count('M') # type: int
            if pile_read.alignment.query_sequence[pile_read.query_position] == var.ref:
                key = 'keep_ref' if count_m >= match_threshold else 'indel_ref' # type: str
            elif pile_read.alignment.query_sequence[pile_read.query_position] in var.alt:
                key = 'keep_alt' if count_m >= match_threshold else 'indel_alt' # type: str
            else:
                key = 'other' # type: str
            stats[key] += 1
        break
    logging.info("Finished getting allelic reads for %s", repr(var))
    logging.debug("Getting allelic reads took %s seconds", round(time.time() - ase_start, 3))
    if stats['keep_ref'] or stats['keep_alt']:
        return AllelicStat.frombpileup(
            var=var,
            ref_count=stats['keep_ref'],
            alt_count=stats['keep_alt'],
            other_count=stats['other'],
            ref_indel=stats['indel_ref'],
            alt_indel=stats['indel_alt']
        )
    else:
        logging.warning("No allelic reads for %s", repr(var))
        return None


def bias_stats(method, coverage): # type(...) -> ...
    pass


def filter_stats(stats, total_coverage=10, allelic_coverage=5, proportion_match=0.8, proportion_indel=0.8): # type: (Iterable[AllelicStat], int, int, float) -> Tuple[AllelicStat, ...]
    """Filter ase stats"""
    logging.info("Filtering ASE results")
    filter_start = time.time()
    stats = tuple(stats) # type: Tuple[AllelicStat, ...]
    stats = filter(lambda x: x.total_count >= total_coverage, stats) # type: Tuple[AllelicStat, ...]
    stats = filter(lambda x: x.ref_count >= allelic_coverage, stats) # type: Tuple[AllelicStat, ...]
    stats = filter(lambda x: x.alt_count >= allelic_coverage, stats) # type: Tuple[AllelicStat, ...]
    stats = filter(lambda x: x.total_count / x.depth >= proportion_match, stats) # type: Tuple[AllelicStat, ...]
    stats = filter(lambda x: x.total_count / (x.total_count + x.ref_indel + x.alt_indel) >= proportion_indel, stats) # type: Tuple[AllelicStat, ...]
    logging.debug("Filtering ASE results took %s seconds", round(time.time() - filter_start, 3))
    return stats
