#!/usr/bin/env python3

from __future__ import division
from __future__ import print_function

import sys
import time
import logging

from collections import defaultdict

if sys.version_info.major == 2:
    from . import ase
    from . import asts
    from . import maths
    from . import utils
    from . import features
    from .fancy_logging import fmttime
else:
    from lorals import ase
    from lorals import asts
    from lorals import maths
    from lorals import utils
    from lorals import features
    from lorals.fancy_logging import fmttime


class AnnotatedStat(ase.AllelicStat, features.GenoVar):

    HEADER = ( # type: Tuple[str, ...]
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
        self._gene = "" # type: str
        self._blacklist = False # type: bool
        self._warning = False # type: bool
        self._map = False # type: bool
        self._null = float('nan') # type: float
        self._pvalue = float('nan') # type: float
        self._qvalue = float('nan') # type: float

    def __str__(self): # type: (None) -> str
        out = ( # type: Tuple[Union[str, int, float], ...]
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
    def gene(self): # type: (None) -> str
        """Gene name for this annotated stat"""
        return self._gene

    @gene.setter
    def gene(self, value): # type: (str) -> None
        if value is None:
            value = "" # type: str
        self._gene = str(value)

    @property
    def blacklisted(self): # type: (None) -> bool
        return self._blacklist

    @blacklisted.setter
    def blacklisted(self, value): # type: (bool) -> None
        self._blacklist = bool(value) # type: bool

    @property
    def warning(self): # type: (None) -> bool
        return self._warning

    @warning.setter
    def warning(self, value): # type: (bool) -> None
        self._warning = bool(value) # type: bool

    @property
    def multi_mapping(self): # type: (None) -> bool
        return self._map

    @multi_mapping.setter
    def multi_mapping(self, value): # type: (bool) -> None
        self._map = bool(value)

    @property
    def null_ratio(self): # type: (None) -> float
        return self._null

    @null_ratio.setter
    def null_ratio(self, value): # type: (float) -> None
        self._null = float(value)

    @property
    def pvalue(self): # type: (None) -> float
        return self._pvalue

    @pvalue.setter
    def pvalue(self, value): # type: (float) -> None
        self._pvalue = float(value)

    @property
    def qvalue(self): # type: (float) -> None
        return self._qvalue

    @qvalue.setter
    def qvalue(self, value): # type: (float) -> None
        self._qvalue = float(value)

    ref_ratio = property(fget=lambda self: self.ref_count / self.total_count, doc="Reference Ratio")
    all_ratio = property(fget=lambda self: self.total_count / self.depth, doc="All ratio")
    indel_ratio = property(
        fget=lambda self: (self.ref_indel + self.alt_indel) / self.total_count,
        doc="Indel ratio"
    )
    other_warning = property(fget=lambda self: self.all_ratio < 0.8, doc="Other allele warning")
    indel_warning = property(fget=lambda self: self.indel_ratio >= 0.2, doc="High indel warning")


def annotate_bed(stats, bedfile): # type: (Iterable[AnnotatedStat], str) -> Tuple[str, ...]
    mstats = tuple(stats[:]) # type: Tuple[AnnotatedStat, ...]
    statfile = features.write_bed(bpileups=mstats, default=True) # type: str
    gfh = asts.bed_intersect(afile=statfile, bfile=bedfile, u=True) # type: str
    my_open = utils.find_open(filename=gfh.fn) # type: function
    with my_open(gfh.fn, 'rt') as gfile:
        regions = tuple(line.strip().split('\t')[3] for line in gfile)
    return regions


def annotate_genes(stats, bedfile): # type: (Iterable[AnnotatedStat], str) -> Tuple[AnnotatedStat, ...]
    logging.info("Annotating genes")
    gene_start = time.time() # type: float
    mstats = tuple(stats[:]) # type: Tuple[AnnotatedStat, ...]
    statfile = features.write_bed(bpileups=mstats, default=True) # type: str
    gfh = asts.bed_intersect(afile=statfile, bfile=bedfile, loj=True) # type: pybedtools.bedtool.BedTool
    my_open = utils.find_open(filename=gfh.fn) # type: function
    dict_genes = defaultdict(list) # type: Mapping[str, List[str]]
    logging.info("Fetching gene names")
    with my_open(gfh.fn, 'rt') as gfile:
        for line in gfile: # type: str
            line = line.strip().split('\t') # type: List[str]
            dict_genes[line[3]].append(line[7])
    logging.info("Matching gene names")
    for i in range(len(mstats)): # type: int
        mstats[i].gene = ','.join(dict_genes.get(mstats[i].default, (None,)))
    logging.debug("Annotating genes took %s seconds", fmttime(start=gene_start))
    return mstats


def annotate_genotypes(stats, vcffile): # type: (Iterable[AnnotatedStat], str) -> Tuple[AnnotatedStat, ...]
    mstats = tuple(stats[:]) # type: Tuple[AnnotatedStat]
    logging.info("Annotating genotypes")
    gt_start = time.time() # type: float
    statfile = features.write_bed(bpileups=mstats, default=True) # type: str
    gfh = asts.bed_intersect(afile=vcffile, bfile=statfile, u=True) # type: pybedtools.bedtool.BedTool
    genotypes = set() # type: Set[features.GenoVar]
    my_open = utils.find_open(filename=gfh.fn) # type: function
    logging.info("Fetching genotypes")
    with my_open(gfh.fn, 'rt') as gfile:
        for line in gfile: # type: str
            genotypes.add(features.GenoVar.fromvcf(vcf=line))
    genotypes = tuple(genotypes) # type: Tuple[features.GenoVar, ...]
    logging.info("Matching genotypes")
    for index, stat in enumerate(mstats): # type: int, AnnotatedStat
        try:
            gt_index = genotypes.index(stat) # type: int
        except ValueError:
            pass
        else:
            mstats[index].geno = genotypes[gt_index].genos
    logging.debug("Annotating genotypes took %s seconds", fmttime(start=gt_start))
    return mstats


def bias_stats(stats, method, coverage): # type(...) -> ...
    logging.info("Calculating bias stats")
    upper_percentile = maths.percentile( # type: float
        x=tuple(x.total_count for x in filter(lambda x: x.total_count >= coverage, stats)),
        probs=75
    )
    bias_stats = dict() # type: Dict[str, List[float, float, int]]
    ref_ratios = defaultdict(list) # type: Mapping[str, List[float]]
    ref_cuml, total_cuml = int(), int() # type: int, int
    all_ratios = list() # type: List[float]
    bias_out = dict() # type: Dict[str, float]
    for stat in stats: # type: AnnotatedStat
        key = stat.ref + stat.alts # type: str
        if key not in bias_stats:
            bias_stats[key] = [0, 0, 0] # type: List[float, float, int]
        if stat.total_count > upper_percentile:
            ref_count = upper_percentile * (stat.ref_count / stat.total_count) # type: float
            alt_count = upper_percentile * (stat.alt_count / stat.total_count) # type: float
        else:
            ref_count, alt_count = stat.ref_count, stat.alt_count # type: int, int
        bias_stats[key][0] += ref_count
        bias_stats[key][1] += ref_count + alt_count
        bias_stats[key][2] += 1
        ref_ratios[key].append(stat.ref_ratio)
        ref_cuml += ref_count
        total_cuml += ref_count + alt_count
    bias_out['genomeWide'] = ref_cuml / total_cuml # type: float
    if method == 'global':
        all_ratios = tuple(x.ref_ratio for x in stats) # type: Tuple[float]
    logging.info("Determining bias")
    for key in bias_stats:
        if bias_stats[key][2] >= 100:
            if method == 'mean':
                bias_out[key] = bias_stats[key][0] / bias_stats[key][1]
            elif method == 'median':
                bias_out[key] = maths.median(ref_ratios[key])
            else:
                bias_out[key] = maths.median(all_ratios)
        elif method == 'mean':
            bias_out[key] = bias_out['genomeWide']
        else:
            bias_out[key] = maths.median(all_ratios)
    return bias_out
