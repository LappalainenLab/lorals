#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

import time
import logging

from operator import itemgetter

def pvalue_adjust(pvalues): # type: (Iterable[float]) -> Iterable[float]
    """Adjust p-values using the Benjamin-Hochberg method
    Basically, it's like R's p.adjust(method = 'fdr')
    >>> x = tuple(map(lambda x: round(x, 1), tuple(x * 0.1 for x in range(0, 10))))
    >>> pvalue_adjust(x)
    (0.0, 0.5, 0.6666666666666667, 0.75, 0.8, 0.8333333333333334, 0.8571428571428571, 0.875, 0.888888888888889, 0.9)
    """
    num_pvalues = len(pvalues) # type: int
    logging.info("Adjusting %s p-values", num_pvalues)
    adjust_start = time.time() # type: float
    order = tuple(index for index, _ in sorted(enumerate(pvalues), key=itemgetter(1), reverse=True)) # type: Tuple[int]
    ranked_order = itemgetter(*(index for index, _ in sorted(enumerate(order), key=itemgetter(1)))) # type: operator.itemgetter
    adj_value = (num_pvalues / index * pvalues[o] for index, o in zip(range(num_pvalues, 0, -1), order)) # type: Iterable[float]
    adj_value = tuple(min(value, 1.0) for value in adj_value) # type: Tuple[float]
    adj_value = ranked_order(adj_value) # type: Tuple[float]
    logging.debug("Adjusting p-values took %s seconds", round(time.time() - adjust_start, 3))
    return adj_value
