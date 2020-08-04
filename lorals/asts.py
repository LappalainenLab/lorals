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

if sys.version_info.major == 2:
    import utils
    from .fancy_logging import fmttime
else:
    from lorals import utils
    from lorals.fancy_logging import fmttime


import pandas
import pybedtools


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
