#!/usr/bin/env python

"""ASTS-specific utilities"""

import os
import time
import logging
import tempfile

from collections import namedtuple, defaultdict, Counter
from typing import Callable, DefaultDict, Dict, Iterable, List, Optional, Tuple, Set

from . import cigar
from . import utils
from .features import Bpileup
from .fancy_logging import fmttime


import pysam
import pandas
import pybedtools

from scipy import stats

__all__: List[str] = [
    "bam_to_bed",
    "bed_intersect",
    "deduplicate_bpileup",
    "vcf_bpileup",
]

LengthStat = namedtuple(
    "LengthStat",
    (
        "contig",
        "position",
        "refAllele",
        "altAllele",
        "D",
        "pvalue"
    )
)

QuantStat = namedtuple(
    'QuantStat',
    (
        'contig',
        'position',
        'refAllele',
        'altAllele',
        'transcript',
        'refCount',
        'altCount',
    )
)

Qname = namedtuple('Qname', ('flairs', 'query', 'vartype'))

NullResult: stats.stats.KstestResult = stats.stats.KstestResult(statistic=utils.nan, pvalue=utils.nan)

class LengthSummary(Bpileup):

    HEADER: Tuple[str, ...] = (
        'contig',
        'position',
        'variantID',
        'refAllele',
        'altAllele',
        'type',
        'count',
    )

    def __init__(
        self,
        chrom: str,
        position: int,
        ref: str,
        alt: str,
        ref_lens: Iterable[int],
        alt_lens: Iterable[int]
    ):
        super().__init__(chrom=chrom, position=position, ref=ref, alt=alt)
        self._rlens: Tuple[int, ...] = tuple(int(x) for x in ref_lens)
        self._alens: Tuple[int, ...] = tuple(int(x) for x in alt_lens)

    def __str__(self) -> str:
        n: int = len(self._rlens) + len(self._alens)
        chroms: Tuple[str, ...] = (self.chrom,) * n
        poss: Tuple[str, ...] = (str(self.position),) * n
        ids: Tuple[str, ...] = (self.default,) * n
        refs: Tuple[str, ...] = (self.ref,) * n
        alts: Tuple[str, ...] = (self.alts,) * n
        types: Tuple[str, ...] = ('ref',) * len(self._rlens) + ('alt',) * len(self._alens)
        counts: Tuple[str, ...] = tuple(str(x) for x in self._rlens + self._alens)
        outs: Tuple[Tuple[str, ...], ...] = tuple(zip(chroms, poss, ids, refs, alts, types, counts))
        outs: Tuple[str, ...] = tuple('\t'.join(x) for x in outs)
        return '\n'.join(outs)


def asts_length(var: Bpileup, bamfile: str, window: int=5, match_threshold: int=8, min_reads: int=20) -> LengthStat:
    """Calculate read length based on SNP"""
    logging.info("Getting read lengths for %s", var.name)
    lengths_start: float = time.time()
    reads_completed: Set[str] = set()
    lengths: DefaultDict[str, List[int]] = defaultdict(list)
    bamfile: str = utils.fullpath(path=bamfile)
    bamfh = pysam.AlignmentFile(bamfile) # type: pysam.libcalcalignment.AlignmentFile
    for pile in bamfh.pileup(region=var.chrom, start=var.dummy, end=var.position, min_base_quality  = 0, max_depth = 100000): # type: pysam.libcalignedsegment.PileupColumn
        if pile.pos != var.dummy:
            continue
        for pile_read in pile.pileups: # type: pysam.libcalignedsegment.PileupRead
            if pile_read.alignment.query_name in reads_completed or not pile_read.query_position:
                continue
            logging.info("Processing read %s", pile_read.alignment.query_name)
            reads_completed.add(pile_read.alignment.query_name)
            pile_window: slice = utils.window(position=pile_read.query_position, size=window)
            count_m: int = cigar.Cigar(tuples=pile_read.alignment.cigartuples)[pile_window].count('M')
            if count_m < match_threshold:
                logging.warning("Skipping read %s", pile_read.alignment.query_name)
                continue
            if pile_read.alignment.query_sequence[pile_read.query_position] == var.ref:
                lengths['ref'].append(len(pile_read.alignment.query_sequence))
            elif pile_read.alignment.query_sequence[pile_read.query_position] in var.alt:
                lengths['alt'].append(len(pile_read.alignment.query_sequence))
    if len(lengths['ref']) >= min_reads and len(lengths['alt']) >= min_reads:
        logging.info("Running KS test")
        ks: stats.stats.KstestResult = stats.ks_2samp(lengths['ref'], lengths['alt'])
        raw: LengthSummary = LengthSummary(
            chrom=var.chrom,
            position=var.position,
            ref=var.ref,
            alt=var.alts,
            ref_lens=lengths['ref'],
            alt_lens=lengths['alt']
        )
    else:
        logging.warning("Too few hits for KS test")
        ks: stats.stats.KstestResult = NullResult
        raw = None
    logging.info("Finished getting lengths")
    logging.debug("Getting read lengths took %s seconds", fmttime(start=lengths_start))
    return LengthStat(contig=var.chrom,
        position=var.position,
        refAllele=var.ref,
        altAllele=var.alts,
        D=ks.statistic,
        pvalue=ks.pvalue
    ), raw


def asts_quant(
    var: Bpileup,
    bamfile: str,
    trans_reads: Dict[str, str],
    window: int=5,
    match_threshold: int=8,
    min_reads: int=10
):
    """..."""
    logging.info("Getting read quants for %s", var.name)
    flair_reads: DefaultDict[str, List[str]] = defaultdict(list)
    quants: List[QuantStat] = list()
    bamfile: str = utils.fullpath(path=bamfile)
    bamfh = pysam.AlignmentFile(bamfile)
    logging.debug("Building pileup on %s from %s to %s", var.chrom, var.dummy, var.position)
    x = var.name
    for pile in bamfh.pileup(region=var.chrom, start=var.dummy, end=var.position, min_base_quality  = 0, max_depth = 100000, truncate = True):
        for pile_read in pile.pileups: # type: pysam.libcalcalignedsegment.PileupRead
            qname: str = pile_read.alignment.query_name
            if pile.pos != var.dummy or not pile_read.query_position:
                continue
            logging.info("Processing read %s", qname)
            pile_window: slice = utils.window(position=pile_read.query_position, size=window)
            count_m: int = cigar.Cigar(tuples=pile_read.alignment.cigartuples)[pile_window].count('M')
            if pile_read.alignment.query_sequence[pile_read.query_position] == var.ref:
                key: str = 'ref' if count_m >= match_threshold else 'ref_indel'
            elif pile_read.alignment.query_sequence[pile_read.query_position] in var.alt:
                key: str = 'alt' if count_m >= match_threshold else 'alt_indel'
            else:
                key: str = 'unknown'
            logging.debug("Read %s classified as %s", qname, key)
            try:
                trans_id = trans_reads[qname]
            except KeyError:
                continue  
            if trans_reads[qname] == None or qname == None:
                continue
            else:
                flair_reads[key].extend([trans_id])             
    if len(flair_reads['ref']) >= min_reads or len(flair_reads['alt']) >= min_reads:
        logging.info("Assembling quant stats")
        counts: DefaultDict[str, Dict[str, int]] = defaultdict(dict)
        for key in ('ref', 'alt'):
            for tx, count in Counter(flair_reads[key]).items(): # type: str, int
                counts[tx][key] = count
        counts: Dict[str, Dict[str, int]] = dict(counts)
        for tx in counts: # type: str
            quants.append(QuantStat(
                contig=var.chrom,
                position=var.position,
                refAllele=var.ref,
                altAllele=var.alts,
                transcript=tx,
                refCount=counts[tx].get('ref', 0),
                altCount=counts[tx].get('alt', 0)
            ))
    else:
        logging.warning("Too few transcripts for %s", repr(var))
        logging.debug("Ref transcripts: %s", len(flair_reads['ref']))
        logging.debug("Alt transcripts: %s", len(flair_reads['alt']))

    return (tuple(quants))


def bam_to_bed(bamfile: str, save: bool=False) -> pybedtools.bedtool.BedTool:
    """Convert a BAM file to a BED file"""
    logging.info("Converting %s from BAM to BED", bamfile)
    bam_bed_start: float = time.time()
    bamfh: pybedtools.bedtool.BedTool = pybedtools.bedtool.BedTool(bamfile)
    bedfh: pybedtools.bedtool.BedTool = bamfh.bam_to_bed()
    if save:
        ofile: str = os.path.splitext(bamfile)[0] + '.bed'
        logging.info("Saving resulting BED file to %s", ofile)
        bedfh: pybedtools.bedtool.BedTool = bedfh.moveto(ofile)
    logging.debug("Converting from BAM to BED took %s seconds", fmttime(start=bam_bed_start))
    return bedfh


def bed_intersect(afile: str, bfile: str, ofile: Optional[str]=None, **kwargs: Optional) -> pybedtools.bedtool.BedTool:
    """Intersect two bedfiles"""
    afile: str = utils.fullpath(path=afile)
    bfile: str = utils.fullpath(path=bfile)
    clean: Set[str] = set()
    logging.info("Connecting to %s", afile)
    if os.path.splitext(afile)[-1] == '.bam':
        afh: pybedtools.bedtool.BedTool = bam_to_bed(bamfile=afile)
        clean.add(afh.fn)
    else:
        afh: pybedtools.bedtool.BedTool = pybedtools.bedtool.BedTool(afile)
    logging.info("Connecting to %s", bfile)
    if os.path.splitext(bfile)[-1] == '.bam':
        bfh: pybedtools.bedtool.BedTool = bam_to_bed(bamfile=bfile)
        clean.add(bfh.fn)
    else:
        bfh: pybedtools.bedtool.BedTool = pybedtools.bedtool.BedTool(bfile)
    logging.info("Intersecting two files")
    intersect_start: float = time.time()
    ifh: pybedtools.bedtool.BedTool = afh.intersect(bfh, **kwargs)
    if isinstance(ofile, str):
        ofile: str = utils.fullpath(path=ofile)
        logging.info("Saving resulting BED file as %s", ofile)
        ifh: pybedtools.bedtool.BedTool = ifh.moveto(ofile)
    logging.info("Number of interesected regions: %s", len(ifh))
    logging.debug("Intersection took %s seconds", fmttime(start=intersect_start))
    if clean:
        logging.debug("Cleaning temporary files")
        for i in clean: # type: str
            os.remove(i)
    return ifh


def deduplicate_bpileup(bpileup: str) -> pandas.core.frame.DataFrame:
    """Deduplicate a BED file"""
    bpileup: str = utils.fullpath(path=bpileup)
    dedup_start: float = time.time()
    bed: pandas.core.frame.DataFrame = pandas.read_csv(
        bpileup,
        sep="\t",
        names=['contig', 'dummy', 'position', 'variantID', 'refAllele', 'altAllele']
    )
    norig: int = len(bed)
    bed.drop_duplicates(keep='first', inplace=True)
    ndedup: int = len(bed)
    if norig == ndedup:
        logging.info("No duplicates in %s", bpileup)
    else:
        logging.info("Deduplicating %s from %s to %s", bpileup, norig, len(bed))
    logging.debug("Deduplication took %s seconds", fmttime(start=dedup_start))
    return bed


def qnames(
    var: Bpileup,
    bamfile: str,
    trans_reads: Dict[str, str],
    window: int=5,
    min_matches: int=8
) -> Tuple[Qname, ...]:
    flair_reads: List[Qname] = list()
    bamfile: str = utils.fullpath(path=bamfile)
    bamfh = pysam.AlignmentFile(bamfile) # type: pysam.libcalcalignment.AlignmentFile
    for pile in bamfh.pileup(region=var.chrom, start=var.dummy, end=var.position,min_base_quality  = 0, max_depth = 100000, truncate = True): # type: pysam.libcalignedsegment.PileupColumn
        if pile.pos != var.dummy:
            continue
        for pile_read in pile.pileups: # type: pysam.libcalignedsegment.PileupRead
            qname: str = pile_read.alignment.query_name
            if qname in reads_completed or not pile_read.query_position:
                continue
            pile_window: slice = utils.window(position=pile_read.query_position, size=window)
            count_m: int = cigar.Cigar(tuples=pile_read.alignment.cigartuples)[pile_window].count('M')
            flairs: Tuple[str, ...] = tuple(ref for query, ref in trans_reads.items() if query == qname)
            if not flairs:
                logging.warning("no flairs")
            if pile_read.alignment.query_sequence[pile_read.query_position] == var.ref:
                vartype: str = 'REF' if count_m >= min_matches else 'REF_INDEL'
            elif pile_read.alignment.query_sequence[pile_read.query_position] in var.alt:
                vartype: str = 'ALT' if count_m >= min_matches else 'ALT_INDEL'
            else:
                vartype: str = 'OTHER'
            flair_reads.append(Qname(flairs=''.join(flairs), query=qname, vartype=vartype))
    return tuple(filter(lambda x: x.vartype != 'OTHER', flair_reads))


def reads_dict(bamfile: str) -> Dict[str, str]:
    with pysam.libcalignmentfile.AlignmentFile(bamfile) as bamfh:
        reads: Dict[str, str] = {read.query_name: read.reference_name for read in bamfh}
    return reads


def split_qnames(qnames: Tuple[Qname, ...]) -> Dict[str, Dict[str, Tuple[str, ...]]]:
    qsplit: Dict = dict.fromkeys({q.vartype for q in qnames})
    for k in qsplit: # type: str
        qsplit[k]: DefaultDict[str, List[str]] = defaultdict(list)
    for q in qnames: # type: Qname
        qsplit[q.vartype][q.flairs].append(q.query)
    for k in qsplit: # type: str
        qsplit[k]: Dict[str, Tuple[str, ...]] = {f: tuple(q) for f, q, in qsplit[k].items()}
    return qsplit


def vcf_bpileup(vcffile: str, pileup: Optional[str]=None) -> str:
    """Create a BED-like pileup from a VCF file"""
    vcffile: str = utils.fullpath(path=vcffile)
    if pileup:
        pileup: str = utils.fullpath(path=pileup)
        pfile = utils.find_open(filename=pileup)(pileup, 'w')
    else:
        pfile = tempfile.NamedTemporaryFile(mode='w', delete=False)
    my_open: Callable = utils.find_open(filename=vcffile)
    with my_open(vcffile, 'rt') as vfile:
        logging.info("Reading in VCF file %s", vcffile)
        pileup_start: float = time.time()
        for line in vfile: # type: str
            if line.startswith('#'):
                continue
            line: List[str] = line.strip().split()
            out: Tuple = (
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
