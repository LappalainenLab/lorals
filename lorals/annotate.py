#!/usr/bin/env python3

import sys
import time
import logging

from numbers import Real
from collections import defaultdict
from typing import Callable, DefaultDict, Dict, Iterable, List, Optional, Set, Tuple

from . import ase
from . import asts
from . import maths
from . import utils
from . import features
from .fancy_logging import fmttime


class AnnotatedStat(ase.AllelicStat, features.GenoVar):

    HEADER: Tuple[str, ...] = (
        'CHR',
        'POS',
        'rsID',
        'REF_ALLELE',
        'ALT_ALLELE',
        'REF_COUNT',
        'ALT_COUNT',
        'TOTAL_COUNT',
        'REF_INDEL_COUNT',
        'ALT_INDEL_COUNT',
        'OTHER_ALLELE_COUNT',
        'RAW_TOTAL_COUNT',
        'GENOTYPE',
        'GENE_ID',
        'GENOTYPE_WARNING',
        'BLACKLIST',
        'MULTI_MAPPING',
        'OTHER_ALLELE_WARNING',
        'HIGH_INDEL_WARNING',
        'NULL_RATIO',
        'BINOM_P',
        'BINOM_P_ADJUSTED',
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
    ) ->  None:
        super(AnnotatedStat, self).__init__(
            chrom=chrom,
            position=position,
            ref=ref,
            alt=alt,
            ref_count=ref_count,
            alt_count=alt_count,
            other_count=other_count,
            ref_indel=ref_indel,
            alt_indel=alt_indel
        )
        self._gene: str = ""
        self._blacklist: bool = False
        self._warning: bool = False
        self._map: bool = False
        self._null: float = float('nan')
        self._pvalue: float = float('nan')
        self._qvalue: float = float('nan')

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
            self.depth,
            self.genos if self.genos else float("nan"),
            self.gene,
            int(self.warning),
            int(self.blacklisted),
            int(self.multi_mapping),
            int(self.other_warning),
            int(self.indel_warning),
            self.null_ratio,
            self.pvalue,
            self.qvalue,
        )
        return '\t'.join(map(str, out))

    @property
    def gene(self) -> str:
        """Gene name for this annotated stat"""
        return self._gene

    @gene.setter
    def gene(self, value: str) -> None:
        if value is None:
            value = "" # type: str
        self._gene = str(value)

    @property
    def blacklisted(self) -> bool:
        return self._blacklist

    @blacklisted.setter
    def blacklisted(self, value: bool) -> None:
        self._blacklist = bool(value) # type: bool

    @property
    def warning(self) -> bool:
        return self._warning

    @warning.setter
    def warning(self, value: bool):
        self._warning = bool(value) # type: bool

    @property
    def multi_mapping(self) -> bool:
        return self._map

    @multi_mapping.setter
    def multi_mapping(self, value: bool) -> None:
        self._map = bool(value)

    @property
    def null_ratio(self) -> float:
        return self._null

    @null_ratio.setter
    def null_ratio(self, value: float) -> None:
        self._null = float(value)

    @property
    def pvalue(self) -> float:
        return self._pvalue

    @pvalue.setter
    def pvalue(self, value: float) -> None:
        self._pvalue = float(value)

    @property
    def qvalue(self) -> float:
        return self._qvalue

    @qvalue.setter
    def qvalue(self, value: float) -> None:
        self._qvalue = float(value)

    ref_ratio = property(fget=lambda self: self.ref_count / self.total_count, doc="Reference Ratio")
    all_ratio = property(fget=lambda self: self.total_count / self.depth, doc="All ratio")
    indel_ratio = property(
        fget=lambda self: (self.ref_indel + self.alt_indel) / self.total_count,
        doc="Indel ratio"
    )
    other_warning = property(fget=lambda self: self.all_ratio < 0.8, doc="Other allele warning")
    indel_warning = property(fget=lambda self: self.indel_ratio >= 0.2, doc="High indel warning")


def annotate_bed(stats: Iterable[AnnotatedStat], bedfile: str) -> Tuple[str, ...]:
    mstats: Tuple[AnnotatedStat, ...] = tuple(stats[:])
    statfile: str = features.write_bed(bpileups=mstats, default=True)
    gfh = asts.bed_intersect(afile=statfile, bfile=bedfile, u=True)
    my_open: Callable = utils.find_open(filename=gfh.fn)
    with my_open(gfh.fn, 'rt') as gfile:
        regions: Tuple[str, ...] = tuple(line.strip().split('\t')[3] for line in gfile)
    return regions


def annotate_genes(stats: Iterable[AnnotatedStat], bedfile: str) -> Tuple[AnnotatedStat, ...]:
    logging.info("Annotating genes")
    gene_start: float = time.time()
    mstats: Tuple[AnnotatedStat, ...] = tuple(stats[:])
    statfile: str = features.write_bed(bpileups=mstats, default=True)
    gfh = asts.bed_intersect(afile=statfile, bfile=bedfile, loj=True) # type: pybedtools.bedtool.BedTool
    my_open: Callable = utils.find_open(filename=gfh.fn)
    dict_genes: DefaultDict[str, List[str]] = defaultdict(list)
    logging.info("Fetching gene names")
    with my_open(gfh.fn, 'rt') as gfile:
        for line in gfile: # type: str
            line: List[str] = line.strip().split('\t')
            dict_genes[line[3]].append(line[7])
    logging.info("Matching gene names")
    for i in range(len(mstats)): # type: int
        mstats[i].gene: str = ','.join(dict_genes.get(mstats[i].default, (None,)))
    logging.debug("Annotating genes took %s seconds", fmttime(start=gene_start))
    return mstats


def annotate_genotypes(stats: Iterable[AnnotatedStat], vcffile: str) -> Tuple[AnnotatedStat, ...]:
    mstats: Tuple[AnnotatedStat, ...] = tuple(stats[:])
    logging.info("Annotating genotypes")
    gt_start: float = time.time()
    statfile: str = features.write_bed(bpileups=mstats, default=True)
    gfh = asts.bed_intersect(afile=vcffile, bfile=statfile, u=True) # type: pybedtools.bedtool.BedTool
    genotypes: Set[features.GenoVar] = set()
    my_open: Callable = utils.find_open(filename=gfh.fn)
    logging.info("Fetching genotypes")
    with my_open(gfh.fn, 'rt') as gfile:
        for line in gfile: # type: str
            genotypes.add(features.GenoVar.fromvcf(vcf=line))
    genotypes: Tuple[features.GenoVar] = tuple(genotypes)
    logging.info("Matching genotypes")
    for index, stat in enumerate(mstats): # type: int, AnnotatedStat
        try:
            gt_index: int = genotypes.index(stat)
        except ValueError:
            pass
        else:
            mstats[index].geno = genotypes[gt_index].genos
    logging.debug("Annotating genotypes took %s seconds", fmttime(start=gt_start))
    return mstats


def bias_stats(stats: Iterable[AnnotatedStat], method: str, coverage: int) -> Dict[str, Real]:
    logging.info("Calculating bias stats")
    upper_percentile: float = maths.percentile(
        x=tuple(x.total_count for x in filter(lambda x: x.total_count >= coverage, stats)),
        probs=75
    )
    bias_stats: Dict[str, List[float, float, int]] = dict()
    ref_ratios: DefaultDict[str, List[float]] = defaultdict(list)
    ref_cuml, total_cuml = int(), int() # type: int, int
    all_ratios: List[float] = list()
    bias_out: Dict[str, Real] = dict()
    for stat in stats: # type: AnnotatedStat
        key: str = stat.ref + stat.alts
        if key not in bias_stats:
            bias_stats[key]: List[float, float, int] = [0, 0, 0]
        if stat.total_count > upper_percentile:
            ref_count: float = upper_percentile * (stat.ref_count / stat.total_count)
            alt_count: float = upper_percentile * (stat.alt_count / stat.total_count)
        else:
            ref_count, alt_count = stat.ref_count, stat.alt_count # type: int, int
        bias_stats[key][0] += ref_count
        bias_stats[key][1] += ref_count + alt_count
        bias_stats[key][2] += 1
        ref_ratios[key].append(stat.ref_ratio)
        ref_cuml += ref_count
        total_cuml += ref_count + alt_count
    bias_out['genomeWide']: float = ref_cuml / total_cuml
    if method == 'global':
        all_ratios: Tuple[float, ...] = tuple(x.ref_ratio for x in stats)
    logging.info("Determining bias")
    for key in bias_stats: # str
        if bias_stats[key][2] >= 100:
            if method == 'mean':
                bias_out[key]: Real = bias_stats[key][0] / bias_stats[key][1]
            elif method == 'median':
                bias_out[key]: Real = maths.median(ref_ratios[key])
            else:
                bias_out[key]: Real = maths.median(all_ratios)
        elif method == 'mean':
            bias_out[key]: Real = bias_out['genomeWide']
        else:
            bias_out[key]: Real = maths.median(all_ratios)
    return bias_out
