#!/usr/bin/env python

__all__ = [ # type: List[str]
    'Bpileup',
    'iter_var',
]

import re
import sys

if sys.version_info.major == 2:
    from .ase import AllelicStat
else:
    from lorals.ase import AllelicStat


class Feature(object):

    """A genomic feature"""

    @staticmethod
    def numeric_chrom(feature): # type: (Feature) -> Union[int, str]
        """Get the numerical value of a Feature's chromosome"""
        try:
            return int(re.search(r'[a-zA-Z](\d+$)', feature.chrom).groups()[0])
        except AttributeError:
            return feature.chrom

    def __init__(self, chrom, start, end, strand='.'): # type: (str, int, int, str) -> None
        if strand not in {'-', '+', '.'}:
            raise ValueError("'strand' must be one of '-', '+', or '.'")
        self._chrom = chrom # type: str
        self._start = int(start) # type; int
        self._end = int(end) # type: int
        self._strand = strand # type: str
        if self._start > self._end:
            raise ValueError("'start' cannot be less than 'end'")

    def __repr__(self): # type: (None) -> str
        return self._chrom + ':' + '-'.join(map(str, (self.start, self.end)))

    def __hash__(self): # type: (None) -> int
        return hash(''.join(map(str, (self.chrom, self.start, self.end))))

    def __eq__(self, other): # type: (Any) -> bool
        if isinstance(other, Feature):
            return self.numeric_chrom(feature=self) == self.numeric_chrom(feature=other) and self.start == other.start and self.end == other.end
        return NotImplemented

    def __ne__(self, other): # type: (Any) -> bool
        return not self == other

    def __le__(self, other): # type: (Any) -> bool
        if isinstance(other, Feature):
            if self.chrom == other.chrom:
                return self.start <= other.start
            return self.numeric_chrom(feature=self) < self.numeric_chrom(feature=other)
        return NotImplemented

    def __lt__(self, other): # type: (Any) -> bool
        if isinstance(other, Feature):
            if self.chrom == other.chrom:
                return self.start < other.start
            return self.numeric_chrom(feature=self) < self.numeric_chrom(feature=other)
        return NotImplemented

    def __ge__(self, other): # type: (Any) -> bool
        return not self < other

    def __gt__(self, other): # type: (Any) -> bool
        return not self <= other

    def __contains__(self, item): # type: (Any) -> bool
        if isinstance(item, Feature):
            return self.numeric_chrom(feature=self) == self.numeric_chrom(feature=item) and self.start <= item.start and self.end >= item.end
        # elif isinstance(item, variants.Variant):
        #     return self.numeric_chrom(feature=self) == item.numeric_chrom(variant=item) and self.start <= item.position and self.end >= item.position
        return False

    chrom = property(fget=lambda self: self._chrom, doc="Chromsome of this feature")
    start = property(fget=lambda self: self._start, doc="Starting position of this feature")
    end = property(fget=lambda self: self._end, doc="Ending position of this feature")
    strand = property(fget=lambda self: self._strand, doc="Strand of this feature")


class Bpileup(Feature):

    """A BED-like pileup feature"""

    @classmethod
    def frominterval(cls, interval): # type: (pybedtools.cbedtools.Interval) -> Bpileup
        """Create a Bpileup from an Interval"""
        if interval.file_type != 'bed':
            raise TypeError("Must be a BED interval")
        elif len(interval.fields) is not 6:
            raise ValueError("wrong number of fields")
        return cls(
            chrom=interval.chrom,
            start=interval.start,
            end=interval.end,
            ref=interval.fields[4],
            alt=interval.fields[5],
            name=interval.fields[3]
        )

    @classmethod
    def fromstat(cls, stat): # type: (ase.AllelicStat) -> Bpileup
        """Create a Bpileup from an AllelicStat"""
        return cls(
            chrom=stat.contig,
            start=stat.position - 1,
            end=stat.position,
            ref=stat.refAllele,
            alt=stat.altAllele,
            name=stat.variantID
        )

    def __init__(self, chrom, start, end, ref, alt, name=None,): # type: (str, int, int, str, str, Optional[str]) -> None
        super(Bpileup, self).__init__(chrom=chrom, start=start, end=end)
        valid_chars = set('ACGT') # type: Set[str, ...]
        ref = str(ref) # type: str
        if not set(ref.upper()).issubset(valid_chars):
            raise ValueError("The reference allele must contain only %s" % ','.join(valid_chars))
        alt = str(alt) # type: str
        if not set(alt.upper()).issubset(valid_chars.union(set(',*'))):
            raise ValueError("The alternate variant must contain only %s" % ','.join(valid_chars))
        self._ref = ref # type: str
        self._alt = tuple(alt.split(',')) # type: Tuple[str, ...]
        if name:
            self._id = str(name) # type: str
        else:
            self._id = str(self) # type: str

    def __str__(self): # type: (None) -> str
        out = ( # type: Tuple[str, str, str]
            super(Bpileup, self).__repr__(),
            self.ref,
            ','.join(self.alt)
        )
        return '_'.join(out)

    def __repr__(self): # type: (None) -> str
        return self.name

    def __hash__(self): # type: (None) -> int
        return hash(''.join(map(
            str,
            (self.chrom, self.start, self.end, self.ref, ','.join(sorted(self.alt)))
        )))

    def __eq__(self, other): # type:(Any) -> bool
        if isinstance(other, Bpileup):
            return hash(self) == hash(other)
        return super(Bpileup, self).__eq__(other=other)

    ref = property(fget=lambda self: self._ref, doc="Reference allele")
    alt = property(fget=lambda self: self._alt, doc="Alternate variant")
    name = property(fget=lambda self: self._id, doc="ID")
    ID = name


def iter_var(df): # type: (pandas.core.frame.DataFrame) -> Iterator(Bpileup)
    """Yield Bpileups from a data frame"""
    if tuple(df.columns) != ('contig', 'dummy', 'position', 'variantID', 'refAllele', 'altAllele'):
        raise ValueError("Invalid data frame")
    for x in df.itertuples():
        yield Bpileup(
            chrom=x.contig,
            start=x.dummy,
            end=x.position,
            ref=x.refAllele,
            alt=x.altAllele,
            name=x.variantID
        )
