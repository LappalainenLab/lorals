#!/usr/bin/env python

"""Some utility functions"""

__all__ = [ # type: List[str, ...]
    "dictsearch",
    "find_open",
    "unpack"
]

import os
import gzip

from collections import Counter

def cigar_parse(tuples):
    """
    arguments:
     <tuples> a CIGAR string tuple list in pysam format
    purpose:
     This function uses the pysam cigarstring tuples format and returns
     a list of tuples in the internal format, [(20, 'M'), (5, "I")], et
     cetera. The zeroth element of each tuple is the number of bases for the
     CIGAR string feature. The first element of each tuple is the CIGAR
     string feature type.
    There are several feature types in SAM/BAM files. See below:
     'M' - match
     'I' - insertion relative to reference
     'D' - deletion relative to reference
     'N' - skipped region from the reference
     'S' - soft clip, not aligned but still in sam file
     'H' - hard clip, not aligned and not in sam file
     'P' - padding (silent deletion from padded reference)
     '=' - sequence match
     'X' - sequence mismatch
     'B' - BAM_CBACK
    """
    # I used the map values from http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment
    psam_to_char = {0: 'M', 1: 'I', 2: 'D', 3: 'N', 4: 'S',
                    5: 'H', 6: 'P', 7: '=', 8: 'X', 9: 'B'}
    return [(value, psam_to_char[feature]) for feature, value in tuples]


def dictsearch(d, query): # type: (Dict[str, Iterable[Any]], Any) -> str
    """Search a dictionary 'd' by value 'query'"""
    for key, value in d.items(): # type: str, Iterable[Any]:
        if query in value:
            return key


def find_open(filename): # type: (str) -> function
    """Figure out which version of open (open vs gzip.open) to use"""
    extension = os.path.splitext(filename)[-1]
    if extension == '.gz':
        return gzip.open
    return open


def get_count_m(tuples, position, window): # type: (Iterable[], int) -> Counter
    """Count the number of matches in a window around the variant of interest"""
    mycigar = cigar_parse(tuples) # type:
    i = 0 # type: int
    s = '' # type: str
    while i < len(mycigar):
        mystring = mycigar[i][0] * mycigar[i][1]
        s += mystring
        i += 1
    wind = s[position - (window+1): position + window] # type: str
    return Counter(wind)


def unpack(collection): # type: (Iterable[Any]) -> Tuple[Any]
    """Unpack a series of nested lists, sets, or tuples"""
    result = [] # type: List
    for item in collection: # type: Any
        if hasattr(item, '__iter__') and not isinstance(item, str):
            result.extend(unpack(collection=item))
        else:
            result.append(item)
    return tuple(result)
