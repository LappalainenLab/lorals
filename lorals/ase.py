#!/usr/bin/env python

"""ASE-specific utilties"""

from __future__ import division
from __future__ import print_function

__all__ = [ # type: List[str, ...]
    'allelic_stats',
]

import sys
import time
import logging

from collections import namedtuple

if sys.version_info.major == 2:
    from . import cigar
    from . import utils
else:
    from lorals import cigar
    from lorals import utils


import pysam
import pandas

AllelicStat = namedtuple(
    'AllelicStat',
    (
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
        'rawDepth'
    )
)

def allelic_stats(var, bamfile, window=5, match_threshold=8): # type: (features.Bpileup, str, int, int) -> pandas.core.frame.DataFrame
    """..."""
    logging.info("Getting allelic reads for %s", str(var))
    ase_start = time.time()
    reads_completed = set() # type: Set[str, ...]
    stats = dict.fromkeys(('keep_ref', 'keep_alt', 'indel_ref', 'indel_alt', 'other'), 0) # type: Dict[str, int]
    bamfile = utils.fullpath(path=bamfile) # type: str
    bamfh = pysam.AlignmentFile(bamfile) # type: pysam.libcalignmentfile.AlignmentFile
    for pile in bamfh.pileup(region=var.chrom, start=var.start, end=var.end): # type: pysam.libcalignedsegment.PileupColumn
        if pile.pos != var.start:
            continue
        for pile_read in pile.pileups: # type: pysam.libcalignedsegment.PileupRead
            if pile_read.alignment.query_name in reads_completed or not pile_read.query_position:
                continue
            logging.info("Processing read %s", pile_read.alignment.query_name)
            reads_completed.add(pile_read.alignment.query_name)
            pile_window = slice(pile_read.query_position - (window + 1), pile_read.query_position + window) # type: slice
            count_m = cigar.Cigar(tuples=pile_read.alignment.cigartuples)[pile_window].count('M') # type: int
            if pile_read.alignment.query_sequence[pile_read.query_position] == var.ref:
                key = 'keep_ref' if count_m >= match_threshold else 'indel_ref' # type: str
            elif pile_read.alignment.query_sequence[pile_read.query_position] in var.alt:
                key = 'keep_alt' if count_m >= match_threshold else 'indel_alt' # type: str
            else:
                key = 'other' # type: str
            stats[key] += 1
        break
    logging.info("Finished getting allelic reads for %s", str(var))
    logging.debug("Getting allelic reads took %s seconds", round(time.time() - ase_start, 3))
    if stats['keep_ref'] or stats['keep_alt']:
        return pandas.DataFrame([AllelicStat(
            contig=var.chrom,
            position=var.end,
            variantID=var.name,
            refAllele=var.ref,
            altAllele=','.join(var.alt),
            refCount=stats['keep_ref'],
            altCount=stats['keep_alt'],
            totalCount=stats['keep_ref'] + stats['keep_alt'],
            refIndelCount=stats['indel_ref'],
            altIndelCount=stats['indel_alt'],
            otherBases=stats['other'],
            rawDepth=stats['keep_ref'] + stats['keep_alt'] + stats['other']
        )])
    else:
        logging.warning("No allelic reads for %s", str(var))
        return pandas.DataFrame(columns=getattr(AllelicStat, '_fields'))
