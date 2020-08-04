#!/usr/bin/env python

__all__ = [ # type: List[str, ...]
    "Cigar"
]

import operator

class Cigar(object):

    """blah"""

    oper_map = { # type: Dict[int, str]
        0: 'M',
        1: 'I',
        2: 'D',
        3: 'N',
        4: 'S',
        5: 'H',
        6: 'P',
        7: '=',
        8: 'X',
        9: 'B'
    }

    def __init__(self, tuples): # type: (Iterable[Tuple[int, int], ...]) -> None
        last = 0 # type: int
        self._starts = list() # type: List[int, ...]
        self._lengths = list() # type: List[int, ...]
        self._chars = list() # type: List[str, ...]
        for oper in tuples: # type: Tuple[int, int]
            if len(oper) != 2 or not all((isinstance(i, int) for i in oper)):
                raise ValueError("Each CIGAR tuple must have exactly two integer values")
            operation, length = oper # type: int, int
            if operation not in self.oper_map:
                raise ValueError("Invalid operation %s, must be one of %s" % (operation, ', '.join(map(str, self.oper_map))))
            self._starts.append(last)
            self._lengths.append(length)
            self._chars.append(self.oper_map.get(operation))
            last += length

    def __str__(self): # type: (None) -> str
        return self[:]

    def __len__(self): # type: (None) -> int
        return sum(self._lengths)

    def __getitem__(self, val): # type: (Union[int, slice]) -> str
        # return NotImplemented
        if not isinstance(val, (int, slice)):
            raise TypeError("Must pass a slice object")
        return ''.join(map(lambda x: operator.mul(*x), zip(self._lengths, self._chars)))[val]
