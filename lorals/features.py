#!/usr/bin/env python

__all__ = [ # type: List[str]
    'Bpileup',
    'iter_var',
    'write_bed',
]

import re
import sys
import tempfile

if sys.version_info.major == 2:
    from . import utils
else:
    from lorals import utils


class Bpileup(object):

    """A BED-like pileup feature"""

    @staticmethod
    def _valid_variant(var, mode='reference'): # type: (str, str) -> None
        if mode not in {'reference', 'alternate'}:
            mode = 'reference' # type: str
        valid_chars = set('ACGT') # type: Set[str]
        msg = 'The %s allele must contain only %s' % (mode, ','.join(valid_chars)) # type: str
        if mode == 'alternate':
            valid_chars.update(',*')
        if not set(var.upper()).issubset(valid_chars):
            raise ValueError(msg)

    @staticmethod
    def default_id(chrom, position, ref, alt): # type: (str, int, str, str) -> str
        Bpileup._valid_variant(var=ref, mode='reference')
        Bpileup._valid_variant(var=alt, mode='alternate')
        return "%(chr)s_%(pos)s_%(ref)s_%(alt)s" % {
            'chr': chrom,
            'pos': position,
            'ref': ref,
            'alt': alt
        }

    @classmethod
    def frominterval(cls, interval): # type: (pybedtools.cbedtools.Interval) -> Bpileup
        """Create a Bpileup from an Interval"""
        if interval.file_type != 'bed':
            raise TypeError("Must be a BED interval")
        elif len(interval.fields) is not 6:
            raise ValueError("wrong number of fields")
        return cls(
            chrom=interval.chrom,
            position=interval.end,
            ref=interval.fields[4],
            alt=interval.fields[5],
            name=interval.fields[3]
        )

    def __init__(self, chrom, position, ref, alt, name=None,): # type: (str, int, str, str, Optional[str]) -> None
        #   Validate data
        ref = str(ref) # type: str
        alt = str(alt) # type: str
        self._valid_variant(var=ref, mode='reference')
        self._valid_variant(var=alt, mode='alternate')
        #   Fill in the data
        self._chrom = str(chrom) # type: str
        self._pos = int(position) # type: int
        self._ref = ref # type: str
        self._alt = tuple(alt.split(',')) # type: Tuple[str, ...]
        if name:
            self._id = str(name) # type: str
        else:
            self._id = self._default_id() # type: str

    def __str__(self): # type: (None) -> str
        out = ( # type: Tuple[Union[int, str], ...]
            self.chrom,
            self.dummy,
            self.position,
            self.name
        )
        return '\t'.join(map(str, out))
        # return '%(chr)s:%(start)s-%(end)s_%(ref)s_%(alt)s' % {
        #     'chr': self.chrom,
        #     'start': self.dummy,
        #     'end': self.position,
        #     'ref': self.ref,
        #     'alt': self.alts
        # }

    def __repr__(self): # type: (None) -> str
        return self.name

    def __hash__(self): # type: (None) -> int
        return hash(self._default_id())

    def __eq__(self, other): # type:(Any) -> bool
        if isinstance(other, Bpileup):
            return hash(self) == hash(other)
        elif isinstance(other, str):
            return self.name == other or self._default_id() == other
        return NotImplemented

    def _default_id(self): # type: (None) -> str
        return Bpileup.default_id(
            chrom=self.chrom,
            position=self.position,
            ref=self.ref,
            alt=self.alts
        )

    def bed(self, default=False): # type(bool) -> str
        """Format this Bpileup to a bed4 format"""
        return "%(chr)s\t%(start)s\t%(end)s\t%(name)s" % {
            'chr': self.chrom,
            'start': self.dummy,
            'end': self.position,
            'name': self._default_id() if default else self.name
        }

    chrom = property(fget=lambda self: self._chrom, doc="Chromsome of this variant")
    position = property(fget=lambda self: self._pos, doc="Position of this variant")
    dummy = property(fget=lambda self: self.position - 1, doc="Dummy starting position")
    ref = property(fget=lambda self: self._ref, doc="Reference allele")
    alt = property(fget=lambda self: self._alt, doc="Alternate variant(s)")
    alts = property(fget=lambda self: ','.join(self.alt), doc="Alternate variant(s) as a string")
    name = property(fget=lambda self: self._id, doc="Variant name/ID")
    default = property(fget=lambda self: self._default_id(), doc="Standardized default ID")


class GenoVar(Bpileup):

    @staticmethod
    def _valid_genotype(genotype): # type: (str) -> None
        if genotype:
            genos = genotype.split(',') # type: List[str]
            for gt in genos: # type: str
                if not re.match(r'[01][/|][01]', gt):
                    raise ValueError("invalid genotypes")
        return None

    @classmethod
    def fromvcf(cls, vcf): # type: (str) -> Bpileup
        vcf = vcf.strip().split('\t') # type: List[str]
        gt_index = vcf[8].split(':').index('GT') # type: int
        gt = tuple(x.split(':')[gt_index] for x in vcf[9:]) # type: Tuple[str, ...]
        return cls(
            chrom=vcf[0],
            position=vcf[1],
            ref=vcf[3],
            alt=vcf[4],
            geno=','.join(gt),
            name=vcf[2]
        )

    def __init__(
        self,
        chrom, # type: str
        position, # type: int
        ref, # type: str
        alt, # type: str
        geno=None, # type: Optional[str]
        name=None # type: Optional[str]
    ): # type: (...) -> None
        super(GenoVar, self).__init__(chrom=chrom, position=position, ref=ref, alt=alt, name=name)
        geno = str(geno) if geno is not None else "" # type: str
        self._valid_genotype(genotype=geno)
        self._geno = tuple(geno.split(',')) # type Tuple[str, ...]

    @property
    def geno(self): # type: (None) -> Tuple[str, ...]
        return self._geno

    @geno.setter
    def geno(self, value): # type: (str) -> None
        self._valid_genotype(genotype=value)
        self._geno = tuple(value.split(',')) # type: Tuple[str, ...]

    genos = property(fget=lambda self: ','.join(self.geno), doc='')


def iter_var(df): # type: (pandas.core.frame.DataFrame) -> Iterator(Bpileup)
    """Yield Bpileups from a data frame"""
    if tuple(df.columns) != ('contig', 'dummy', 'position', 'variantID', 'refAllele', 'altAllele'):
        raise ValueError("Invalid data frame")
    for x in df.itertuples():
        yield Bpileup(
            chrom=x.contig,
            position=x.position,
            ref=x.refAllele,
            alt=x.altAllele,
            name=x.variantID
        )


def write_bed(bpileups, default=False, ofile=None): # type: (Iterable[Bpileup], bool, Optional[str]) -> str
    if ofile:
        ofh = utils.find_open(filename=ofile)(ofile, 'w+b')
    else:
        ofh = tempfile.NamedTemporaryFile(delete=False)
    for bp in bpileups: # type: Bpileup
        ofh.write(bp.bed(default=default))
        ofh.write('\n')
        ofh.flush()
    ofh.close()
    return ofh.name
