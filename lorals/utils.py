#!/usr/bin/env python

"""Some utility functions"""

import os
import gzip
import operator
import itertools

from collections import Counter
from typing import Any, Callable, Dict, Iterable, List, Tuple

__all__: List[str] = [
    "nan",
    "dictsearch",
    "find_open",
    "fullpath",
    "unpack",
    "where",
    "window",
]

nan: float = float('nan')

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


def dictsearch(d: Dict[str, Iterable], query: Any) -> str:
    """Search a dictionary 'd' by value 'query'"""
    for key, value in d.items(): # type: str, Iterable[Any]:
        if query in value:
            return key


def find_open(filename: str) -> Callable:
    """Figure out which version of open (open vs gzip.open) to use

    Standard open cannot handle gzipped files easily. To get around this, this
    function will return either the `gzip.open` function if the filename ends in
    ".gz" or the `open` function otherwise

    Arguments:
        filename (str): A file to find the correct open for

    Returns:
        function: An function capable of opening the file described by `filename`
    """
    extension = os.path.splitext(filename)[-1]
    if extension == '.gz':
        return gzip.open
    return open


def fullpath(path: str) -> str:
    """Find the full path to a file or directory"""
    return os.path.realpath(os.path.expanduser(path))


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


def unpack(collection: Iterable) -> Tuple:
    """Unpack a series of nested lists, sets, or tuples"""
    result: List = []
    for item in collection: # type: Any
        if hasattr(item, '__iter__') and not isinstance(item, str):
            result.extend(unpack(collection=item))
        else:
            result.append(item)
    return tuple(result)


def where(name: str, flags: int=os.X_OK) -> str:
    dirname, name = os.path.split(name) # type: str, str
    #   Figure out the paths we're working with
    paths: Tuple[str, str] = (dirname, os.getcwd())
    paths: Tuple[str, ...] = paths + tuple(os.environ.get('PATH', '').split(os.pathsep))
    paths: Tuple[str, ...] = tuple(filter(None, paths))
    #   Figure out all the extensions we need to test
    name, extensions = os.path.splitext(name) # type: str, str
    if extensions:
        extensions: Tuple[str] = (extensions,)
    else:
        extensions: List[str] = os.environ.get('PATHEXT', '').split(os.pathsep)
        extensions: Tuple[str, ...] = ('',) + tuple(filter(None, extensions))
    #   Get all the combinations of name + extension we're testing
    names: itertools.product = itertools.product([name], extensions)
    names: Tuple[str, ...] = tuple(map(lambda x: operator.add(*x), names))
    #   Find the first matching executable
    for path in paths: # type: str
        for exe in (os.path.join(path, e) for e in names): # type: str
            if os.access(exe, flags):
                return exe
    else:
        raise ValueError("Cannot find %(name)s%(ext)s" % {'name': name, 'ext': extensions[0]})


def window(position: int, size: int) -> slice:
    """Make a window slice"""
    return slice(position - (size + 1), position + size)
