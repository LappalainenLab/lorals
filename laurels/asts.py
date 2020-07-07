#!/usr/bin/env python

"""ASTS-specific utilities"""

from __future__ import division
from __future__ import print_function

__all__ = [ # type: List[str, ...]
    "vcf_bpileup"
]

import os
import sys
import time
import logging
import tempfile

if sys.version_info.major == 2:
    from utils import find_open
else:
    from laurels.utils import find_open

import pybedtools


def bam_to_bed(bamfile, save=False): # type: (str, bool) -> pybedtools.bedtool.BedTool
    """Convert a BAM file to a BED file"""
    logging.info("Converting %s from BAM to BED")
    bam_bed_start = time.time()
    bamfh = pybedtools.bedtool.BedTool(bamfile) # type: pybedtools.bedtool.BedTool
    bedfh = bamfh.bam_to_bed() # type: pybedtools.bedtool.BedTool
    if save:
        ofile = os.path.splitext(bamfile)[0] + '.bed' # type: str
        logging.info("Saving resulting BED file to %s", ofile)
        bedfh = bedfh.moveto(ofile) # type: pybedtools.bedtool.BedTool
    logging.debug("Converting from BAM to BED took %s seconds", round(time.time() - bam_bed_start, 3))
    return bedfh


def bed_intersect(afile, bfile, ofile=None): # type: (str, str, Optional[str]) -> pybedtools.bedtool.BedTool
    """Intersect two bedfiles"""
    logging.info("Connecting to %s", afile)
    if os.path.splitext(afile)[-1] == '.bam':
        afh = bam_to_bed(bamfile=afile) # type: pybedtools.bedtool.BedTool
    else:
        afh = pybedtools.bedtool.BedTool(afile) # type: pybedtools.bedtool.BedTool
    logging.info("Connecting to %s", bfile)
    if os.path.splitext(bfile)[-1] == '.bam':
        bfh = bam_to_bed(bamfile=bfile) # type: pybedtools.bedtool.BedTool
    else:
        bfh = pybedtools.bedtool.BedTool(bfile) # type: pybedtools.bedtool.BedTool
    logging.info("Intersecting two files")
    intersect_start = time.time()
    ifh = afh.intersect(bfh) # type: pybedtools.bedtool.BedTool
    if isinstance(ofile, str):
        logging.info("Saving resulting BED file as %s", ofile)
        ifh = ifh.moveto(ofile) # type: pybedtools.bedtool.BedTool
    logging.info("Number of regions in afile: %s", len(afh))
    logging.info("Number of interesected regions: %s", len(ifh))
    logging.debug("Intersection took %s seconds", round(time.time() - intersect_start, 3))
    return ifh


def vcf_bpileup(vcffile, pileup=None): # type: (str, Optional[str]) -> str
    """Create a BED-like pileup from a VCF file"""
    if pileup:
        pfile = find_open(pileup)(pileup, 'w+b')
    else:
        pfile = tempfile.NamedTemporaryFile(delete=False)
    my_open = find_open(vcffile) # type: function
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
    logging.debug("Building pileup took %s seconds", round(time.time() - pileup_start, 3))
    return pfile.name
