#!/usr/bin/env python3

import sys
import time
import logging

from collections import defaultdict

if sys.version_info.major == 2:
    from . import ase
    from . import asts
    from . import utils
    from . import features
    from .fancy_logging import fmttime
else:
    from lorals import ase
    from lorals import asts
    from lorals import utils
    from lorals import features
    from lorals.fancy_logging import fmttime


class AnnotatedStat(ase.AllelicStat, features.GenoVar):

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

    @property
    def gene(self): # type: (None) -> str
        """Gene name for this annotated stat"""
        return self._gene

    @gene.setter
    def gene(self, value): # type: (str) -> None
        if value is None:
            value = "" # type: str
        self._gene = str(value)

    ref_ratio = property(fget=lambda self: self.ref_count / self.total_count, doc="Reference Ratio")
    all_ratio = property(fget=lambda self: self.total_count / self.depth, doc="All ratio")
    indel_ratio = property(
        fget=lambda self: (self.ref_indel + self.alt_indel) / self.total_count,
        doc="Indel ratio"
    )
    other_warning = property(fget=lambda self: self.all_ratio < 0.8, doc="Other allele warning")
    indel_warning = property(fget=lambda self: self.indel_ratio >= 0.2, doc="High indel warning")


def annotate_genes(stats, bedfile): # type: (Iterable[AnnotatedStat], str) -> Tuple[AnnotatedStat, ...]
    logging.info("Annotating genes")
    gene_start = time.time() # type: float
    mstats = tuple(stats[:]) # type: Tuple[AnnotatedStat]
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
