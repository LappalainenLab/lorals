#!/usr/bin/env python

from __future__ import division
from __future__ import print_function

import math
import time
import logging

from operator import itemgetter

def median(x): # type: (Iterable[Union[int, float]]) -> Union[int, float]
    """Calculate the median
    >>> median([1, 2, 3])
    2
    >>> median([4, 3, 4, 1, 2])
    3
    """
    if isinstance(x, (int, float)):
        return x
    x_sorted = tuple(sorted(x)) # type: Tuple[Union[int, float]]
    midpoint = int(len(x_sorted) / 2) # type: int
    if len(x_sorted) % 2:
        return x_sorted[midpoint]
    elif len(x_sorted) == 2:
        lower, upper = x_sorted # type: Union[int, float], Union[int, float]
    else:
        lower, upper = x_sorted[midpoint:(midpoint + 2)] # type: Union[int, float], Union[int, float]
    return upper - ((upper - lower) / 2)


def percentile(x, probs): # type: (Iterable[Union[int, float]], Union[Tuple[Union[int, float]], int, float]) -> Union[int, float]
    """Calculate a percentile
    >>> percentile(x = (15, 20, 35, 40, 50), probs = 40)
    29.0
    >>> percentile(x = (15, 20, 35, 40, 50), probs = (5, 40, 95))
    (16.0, 29.0, 48.0)
    """
    if isinstance(probs, (int, float)):
        probs = (probs,) # type: Tuple[Union[int, float]]
    pcheck = (0 <= p <= 100 for p in probs) # type: Generator(bool)
    if not all(pcheck):
        raise ValueError("All percentile values must be between 0 and 100, inclusive")
    xs = tuple(sorted(x)) # type: Tuple[Union[int, float]]
    ranks = tuple(p / 100 * (len(xs) - 1) + 1 for p in probs) # type: Tuple[float]
    indexes = tuple(int(math.floor(r)) for r in ranks) # type: Tuple[int]
    percs = tuple(xs[i - 1] + ((r % 1) * (xs[i] - xs[i - 1])) for r, i in zip(ranks, indexes))
    if len(percs) == 1:
        return percs[0]
    return percs


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
