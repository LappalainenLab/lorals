#!/usr/bin/env python

"""ASE-specific utilties"""

import time
import logging

from typing import Dict, Iterable, List, Optional, Set, Tuple

from . import cigar
from . import utils
from .features import Bpileup

import pysam

__all__: List[str] = [
    'AllelicStat',
    'allelic_stats',
    'filter_stats',
]

class AllelicStat(Bpileup):

    HEADER: Tuple[str, ...] = (
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
    def fromstring(cls, string: str, sep: str='\t') -> 'AllelicStat':
        string: str = string.strip().split(sep)
        # if len(string) != len(AllelicStat.HEADER):
        if len(string) < len(AllelicStat.HEADER):
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
            var: Bpileup,
            ref_count: int,
            alt_count: int,
            other_count: int,
            ref_indel: int,
            alt_indel: int
    ) -> 'AllelicStat':
        return cls(
            chrom=var.chrom,
            position=var.position,
            ref=var.ref,
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
            chrom: str,
            position: int,
            ref: str,
            alt: str,
            ref_count: int,
            alt_count: int,
            other_count: int,
            ref_indel: int,
            alt_indel: int,
            name: Optional[str]=None
    ) -> None:
        super(AllelicStat, self).__init__(
            chrom=chrom,
            position=position,
            ref=ref,
            alt=alt,
            name=name
        )
        self._rc: int = int(ref_count)
        self._ac: int = int(alt_count)
        self._oc: int = int(other_count)
        self._ric: int = int(ref_indel)
        self._aic: int = int(alt_indel)

    def __str__(self) -> str:
        out: Tuple = (
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

    def print_bpileup(self) -> str:
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


def allelic_stats(var: Bpileup, bamfile: str, window: int=5, match_threshold: int=8) -> AllelicStat:
    """..."""
    logging.info("Getting allelic reads for %s", repr(var))
    ase_start: float = time.time()
    reads_completed: Set[str] = set()
    stats: Dict[str, int] = dict.fromkeys(('keep_ref', 'keep_alt', 'indel_ref', 'indel_alt', 'other'), 0)
    bamfile: str = utils.fullpath(path=bamfile)
    bamfh = pysam.AlignmentFile(bamfile) # type: pysam.libcalignmentfile.AlignmentFile
    for pile in bamfh.pileup(region=var.chrom, start=var.dummy, end=var.position): # type: pysam.libcalignedsegment.PileupColumn
        if pile.pos != var.dummy:
            continue
        for pile_read in pile.pileups: # type: pysam.libcalignedsegment.PileupRead
            if pile_read.alignment.query_name in reads_completed or not pile_read.query_position:
                continue
            logging.info("Processing read %s", pile_read.alignment.query_name)
            reads_completed.add(pile_read.alignment.query_name)
            pile_window: slice = utils.window(position=pile_read.query_position, size=window)
            count_m: int = cigar.Cigar(tuples=pile_read.alignment.cigartuples)[pile_window].count('M')
            if pile_read.alignment.query_sequence[pile_read.query_position] == var.ref:
                key: str = 'keep_ref' if count_m >= match_threshold else 'indel_ref'
            elif pile_read.alignment.query_sequence[pile_read.query_position] in var.alt:
                key: str = 'keep_alt' if count_m >= match_threshold else 'indel_alt'
            else:
                key: str = 'other'
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


def filter_stats(
    stats: Iterable[AllelicStat],
    total_coverage: int=10,
    allelic_coverage: int=5,
    proportion_match: int=0.8,
    proportion_indel: int=0.8
) -> Tuple[AllelicStat, ...]:
    """Filter ase stats"""
    logging.info("Filtering ASE results")
    filter_start: float = time.time()
    stats: Tuple[AllelicStat, ...] = tuple(stats)
    stats: Iterable[AllelicStat] = filter(lambda x: x.total_count >= total_coverage, stats)
    stats: Iterable[AllelicStat] = filter(lambda x: x.ref_count >= allelic_coverage, stats)
    stats: Iterable[AllelicStat] = filter(lambda x: x.alt_count >= allelic_coverage, stats)
    stats: Iterable[AllelicStat] = filter(lambda x: x.total_count / x.depth >= proportion_match, stats)
    stats: Iterable[AllelicStat] = filter(lambda x: x.total_count / (x.total_count + x.ref_indel + x.alt_indel) >= proportion_indel, stats)
    logging.debug("Filtering ASE results took %s seconds", round(time.time() - filter_start, 3))
    return tuple(stats)
